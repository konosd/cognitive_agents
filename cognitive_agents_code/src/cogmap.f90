! Cognitive Maps - simulate Brownian agents with repulsion
!
! Copyright (C) 2018 Maksym Romenskyy [m.romenskyy@imperial.ac.uk]
! Copyright (C) 2019 Louis Cochen [louis.cochen15@imperial.ac.uk],
!                    Klara Kaleb [klara.kaleb18@imperial.ac.uk],
!                    Alvin Ziqi Lu [ziqi.lu15@imperial.ac.uk],
!                    Marta Rudzite [marta.rudzite18@imperial.ac.uk],
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
! Please see the accompanying LICENSE file for further information or
! find a copy at http://www.gnu.org/licenses/gpl-3.0.html
!
! This code is written for compilation by gfortran 8.2.x, the use of a
! different compiler, or different version might cause unexpected behaviour.

program cognitive_maps

  use utilities
  use io_utils
  use interactions
  use ziggurat
  use biKDE

  ! ### CONSTANTS ###
  ! Simulation parameters
  ! integer, parameter    :: N_AGENT = 200
  ! integer, parameter    :: N_STEP = 2000
  integer, parameter    :: N_VTRAJ = 360
  integer, parameter    :: N_AGENT_STAT = 0 ! no of static / non-moving agents
  ! integer, parameter    :: N_VSTEP = 20
  real, parameter       :: DIAMETER = 2
  real, parameter       :: BOX_LENGTH = 80
  ! real, parameter       :: BOX_AREA = BOX_LENGTH ** 2.0
  real, parameter       :: TIMESTEP = 0.1
  character, parameter  :: VREPUL_TYPE = "S" ! Hard, Soft, sWarm, Ahc, None, stop Virtual
  character, parameter  :: REPUL_TYPE = "S" ! Hard, Soft, sWarm, Ahc
  real, parameter       :: FRICTION = 1.0 ! gamma
  real, parameter       :: TEMPERATURE = 1.0
  real, parameter       :: NOI = sqrt(2.0*FRICTION*TEMPERATURE/TIMESTEP)
  real, parameter       :: REPUL_STRENGTH = 20.0
  real, parameter       :: REPUL_EXP = 60.0
  ! real, parameter       :: PI = 3.14159265358979323846264338327950288419716939937510582097494459 ! from https://oeis.org/A000796
  ! real, parameter       :: TAU = 6.28318530717958647692528676655900576839433879875021164194988918 ! from https://oeis.org/A019692
  integer, parameter    :: INCR_COOR = 1
  ! real, parameter       :: empty = 0.00
  character, parameter  :: CogmapTYPE = 'L' ! <R>adiusGyration <L>ogFraction 

  ! TODO include analysis parameters

  ! ### VARIABLES ###
  integer           :: i, step, vtraj, vstep, vsteps, d ! iterators
  integer           :: n_agent, n_vstep, n_step ! , n_vtraj change constants to variables
  integer           :: start, incr, vrecord_agent
  real              :: r_gyration_norm, stat_radius, log_fraction
  character(len=90) :: path, dir_name, file_name, arg_n_vtraj
  character(len=90) :: arg_n_agent, arg_n_vstep, arg_n_step, arg_stat_r !, arg_repul
  character(len=90) :: box_length_char, n_step_char
  integer           :: coor_file_code = 67, vrecord_file_code = 86
  logical           :: stop_traj

  real, dimension(:,:), allocatable :: coordinates
  real, dimension(:,:), allocatable :: p_coordinates, c_coordinates, velocities
  real, dimension(:,:), allocatable :: v_coordinates, v_velocities
  real, dimension(:,:), allocatable :: rnumb, vrecord
  real, dimension(:,:), allocatable :: v_noise
  real, dimension(:,:), allocatable :: g_array, g_array_coor, vrecord2
  real, dimension(4) :: rand4
  real, dimension(2) :: f_langevin, f_repul, f_total, i_coordinates
  real, dimension(2) :: v_p_coordinates, v_c_coordinates, v_velocity, p_velocity
  real, dimension(2) :: gauss, mean_r_virt, r_gyr_comp
  real, dimension(:), allocatable :: r_gyration
  integer, dimension(:), allocatable :: path_size

  integer :: path_sum, y_index, x_index, index
  real, dimension(:), allocatable :: log_fractions, mean_log_fractions, v_coordinates2
  real, dimension(:,:), allocatable :: vr

  real, parameter :: incr_kde = 0.1
  character, parameter :: btype = 'S', kernel = 'U'

  ! ###SYNTHETIC###
  ! FIX THE SEEDS
  integer, allocatable :: seed(:)
  integer seed_size
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed(:) = 10
  call random_seed(put=seed)
  call zigset(10)
  deallocate(seed)

  ! COMMAND LINE ARGUMENTS
  call get_command_argument(1,arg_n_agent)
  call get_command_argument(2,arg_n_vstep)
  call get_command_argument(3,arg_n_step)
  call get_command_argument(4,file_name)
  ! call get_command_argument(5,arg_n_vtraj)
  ! call get_command_argument(5,arg_stat_r)     ! for circle SYNTHETIC
  ! call get_command_argument(5,arg_repul)      ! for adjustable repulsion

  read(arg_n_agent,*)N_AGENT   ! convert to int
  ! N_AGENT_STAT = N_AGENT-1      ! for circle + 1 moving agent sim only
  read(arg_n_vstep,*)N_VSTEP
  vsteps = N_VSTEP + 1
  read(arg_n_step,*)N_STEP
  ! read(arg_n_vtraj,*)n_vtraj
  ! read(arg_stat_r,*)stat_radius
  ! stat_radius = 10
  ! read(arg_repul,*)REPUL_EXP
  ! the 'file_name' argument stays a character

  ! directory name variables
  write(box_length_char,*)int(BOX_LENGTH)
  write(n_step_char,*)int(N_STEP)

  ! Size allocation
  allocate(coordinates(2,N_AGENT))
  allocate(p_coordinates(2,N_AGENT))
  allocate(c_coordinates(2,N_AGENT))
  allocate(velocities(2,N_AGENT))
  allocate(v_coordinates(2,N_VSTEP+1))
  allocate(rnumb(2,N_VSTEP+1))
  allocate(v_noise(2,N_VTRAJ))
  allocate(r_gyration(N_VTRAJ))
  ! allocate(vrecord(N_VTRAJ,N_VSTEP*2))

  allocate(vrecord2(N_VSTEP*N_VTRAJ,2))
  allocate(path_size(N_VTRAJ))

  allocate(vr(N_VTRAJ,2))
  allocate(v_coordinates2(2))
  allocate(log_fractions(N_VTRAJ*N_VSTEP))
  allocate(mean_log_fractions(N_VTRAJ))

  ! OUTPUT SETUP
  call init_directory(arg_n_agent, arg_n_vstep, &
                      box_length_char, REPUL_TYPE, n_step_char, dir_name)

  path = trim(dir_name)//'/'//trim(file_name)

  call write_info_file(path, N_STEP, TIMESTEP, FRICTION, TEMPERATURE, &
                       start, incr, INCR_COOR)

  call init_files(path, coor_file_code, vrecord_file_code)

  ! Initial conditions

  !call init_random_seed() ! TURN OFF IN SYNTHETIC

  coordinates = init_agents(N_AGENT, DIAMETER, BOX_LENGTH) ! TURN OFF IN SYNTHETIC

  ! ###SYNTHETIC###

  ! PLACE AGENTS
  !coordinates = place_agents(N_AGENT, 1.0) ! for set initial coor agents SYNTHETIC

  ! SYNTHETIC 2 HORIZONTAL
  ! coordinates(:,1) = [38,40]
  ! coordinates(:,2) = [42,40]
  ! SYNTHETIC 2 HORIZONTAL 2
  ! coordinates(:,1) = [02,40]
  ! coordinates(:,2) = [78,40]
  ! SYNTHETIC 2 VERTICAL
  ! coordinates(:,1) = [40,38]
  ! coordinates(:,2) = [40,42]
  ! SYNTHETIC 2 DIAGONAL
  ! coordinates(:,1) = [38,38]
  ! coordinates(:,2) = [42,42]
  ! SYNTHETIC 3 TRIANGLE
  ! coordinates(:,1) = [38,39]
  ! coordinates(:,2) = [42,39]
  ! coordinates(:,3) = [40,42]

  ! vrecord_agent = N_AGENT  ! virtual traj record particle number (<np)
  p_coordinates = coordinates ! past coordinates
  c_coordinates = coordinates ! corrected coordinates

  do i=1, N_AGENT
    velocities(1,i) = rnor() ! TURN OFF IN SYNTHETIC
    velocities(2,i) = rnor() ! TURN OFF IN SYNTHETIC
  end do

  ! velocities = 0

  do i=1, n_agent  ! print initial coordinates
    write(coor_file_code, *) c_coordinates(1,i),c_coordinates(2,i), &
        velocities(2,i), velocities(2,i)
  end do

  ! Main loop
  do step=1, N_STEP ! for every step in time

     ! used to print coord for non moving agents
      !do j=1, N_AGENT_STAT
       ! write(coor_file_code, *) coordinates(1,j),coordinates(2,j), &
        !    empty, empty
      !end do

    do i=N_AGENT_STAT+1, N_AGENT ! for every particle in that timestep
    ! do i=1, N_AGENT ! for every particle in that timestep

      ! vrecord2 = 0.0

      do vtraj=1, N_VTRAJ ! for every virt traj for that particle in that timestep

        ! Setup for the first virtual step
        v_coordinates(1,1) = coordinates(1,i)
        v_coordinates(2,1) = coordinates(2,i)
        v_c_coordinates(1) = c_coordinates(1,i)
        v_c_coordinates(2) = c_coordinates(2,i)
        v_p_coordinates(1) = p_coordinates(1,i)
        v_p_coordinates(2) = p_coordinates(2,i)
        v_velocity(1) = velocities(1,i)
        v_velocity(2) = velocities(2,i)

        ! Generate random number for the whole vtraj
        do vstep=1, N_VSTEP ! for every step in virtual trajectory
          rnumb(1,vstep) = rnor()
          rnumb(2,vstep) = rnor()
        end do

        ! Save the first vstep random numbers
        v_noise(1,vtraj) = rnumb(1,1)
        v_noise(2,vtraj) = rnumb(2,1)

        do vstep=1, N_VSTEP ! for every step in virtual trajectory

          if(VREPUL_TYPE.ne.'V')then
            call repulsions(c_coordinates, v_c_coordinates, i, N_AGENT, BOX_LENGTH, &
                            DIAMETER, REPUL_STRENGTH, REPUL_EXP, VREPUL_TYPE, f_repul)
          else
            stop_traj = virtual_repulsion(c_coordinates, i, v_c_coordinates, BOX_LENGTH, &
                                   DIAMETER)
            if(stop_traj)then
              vsteps = vstep + 1 ! record number of steps taken
              ! THIS IS A HACK
              v_coordinates(1,vstep:(N_VSTEP+1)) = i_coordinates(1)
              v_coordinates(2,vstep:(N_VSTEP+1)) = i_coordinates(2)
              exit
            end if
          end if        ! end if for equals 'V' check

          ! virtual langevin equation
          f_langevin(1) = -FRICTION*v_velocity(1) + NOI*rnumb(1,vstep)
          f_langevin(2) = -FRICTION*v_velocity(2) + NOI*rnumb(2,vstep)

          ! virtual total force
          f_total(1) = f_langevin(1) + f_repul(1)
          f_total(2) = f_langevin(2) + f_repul(2)

          !if((f_repul(1).ne.0.0).or.(f_repul(2).ne.0.0))then

           ! p_velocity = v_velocity

          !end if

          ! virtual integration
          ! call integration(TIMESTEP, v_coordinates(:,vstep), i_coordinates, &
          !                  v_p_coordinates, f_total, v_velocity)

          ! integrate
          i_coordinates(1) = 2*v_coordinates(1,vstep) - v_p_coordinates(2) + f_total(2)*TIMESTEP**2
          i_coordinates(2) = 2*v_coordinates(2,vstep) - v_p_coordinates(2) + f_total(2)*TIMESTEP**2

          ! measure new velocity
          v_velocity(1) = (i_coordinates(1)-v_coordinates(1,vstep)) / TIMESTEP
          v_velocity(2) = (i_coordinates(2)-v_coordinates(2,vstep)) / TIMESTEP

          ! update past coordinates
          v_p_coordinates(1) = v_coordinates(1,vstep)
          v_p_coordinates(2) = v_coordinates(2,vstep)

          ! if((f_repul(1).ne.0.0).or.(f_repul(2).ne.0.0))then
          !   print *, sqrt(sum(v_velocity**2))/sqrt(sum(p_velocity**2))
          ! end if

          v_coordinates(1,vstep+1) = i_coordinates(1)
          v_coordinates(2,vstep+1) = i_coordinates(2)

          v_c_coordinates(1) = modulo(v_coordinates(1,vstep), BOX_LENGTH)
          v_c_coordinates(2) = modulo(v_coordinates(2,vstep), BOX_LENGTH)

        end do

        path_size(vtraj) = vsteps

        if (CogmapTYPE=='R') then

          ! radius of gyration for trajectory i
          mean_r_virt(1) = sum(v_coordinates(1,:)) / vsteps
          mean_r_virt(2) = sum(v_coordinates(2,:)) / vsteps

          r_gyr_comp = 0

          do vstep=1, vsteps
            r_gyr_comp(1) = r_gyr_comp(1)+(v_coordinates(1,vstep)-mean_r_virt(1))**2
            r_gyr_comp(2) = r_gyr_comp(2)+(v_coordinates(2,vstep)-mean_r_virt(2))**2
          end do

          r_gyr_comp(1) = r_gyr_comp(1) / vsteps
          r_gyr_comp(2) = r_gyr_comp(2) / vsteps
          r_gyration(vtraj) = sqrt(sum(r_gyr_comp))

        else

          do vstep=1,N_VSTEP
            vrecord2((vtraj-1)*N_VSTEP+vstep,1)=MODULO(v_coordinates(1,vstep+1),BOX_LENGTH)
            vrecord2((vtraj-1)*N_VSTEP+vstep,2)=MODULO(v_coordinates(2,vstep+1),BOX_LENGTH)
          end do

         ! print *, 'vrecord2' ,vrecord2

        end if

        ! if (i.eq.vrecord_agent) then
        !   do vstep=1,N_VSTEP
        !     vrecord(vtraj,vstep)= modulo(v_coordinates(1,vstep+1),BOX_LENGTH)
        !     vrecord(vtraj,N_VSTEP+vstep)= modulo(v_coordinates(2,vstep+1),BOX_LENGTH)
        !   end do
        ! end if

      end do
      ! choosing which way to go

      gauss = 0

      if (CogmapTYPE=='R') then

        do vtraj=1, N_VTRAJ
        ! if (r_gyration(vtraj).ne.0.0)then
          r_gyration_norm = log(r_gyration(vtraj) / (sum(r_gyration)/size(r_gyration)))
          gauss(1) = gauss(1) + v_noise(1,vtraj)*r_gyration_norm
          gauss(2) = gauss(2) + v_noise(2,vtraj)*r_gyration_norm
        ! end if
        end do

      else

        do vstep=1, N_VSTEP

          vr(:,:) = 0

          do vtraj=1, N_VTRAJ

            vr(vtraj,:) = vrecord2((vtraj-1)*N_VSTEP+vstep,:)

          end do 

          ! at every level
          call get_biKDE(vr, incr_kde, btype, kernel, g_array)

          do vtraj=1, N_VTRAJ

            v_coordinates2 = vr(vtraj,:)

            x_index = MINLOC(abs(g_array(:,1)-v_coordinates2(1)),dim = 1)
            y_index = MINLOC(abs(g_array(:,2)-v_coordinates2(2)),dim = 1)

            index = x_index + y_index - 1 

            if(g_array(index,3).eq.0.0)then
                 log_fraction = 0.0
            else
                 log_fraction = log(g_array(index,3)/sum(g_array(:,3)))
            end if 

            log_fractions((vtraj-1)*N_VSTEP+vstep) = log_fraction

          end do

        end do

        ! take ensemble average log_fraction for every vtraj

        do vtraj=1, N_VTRAJ

            mean_log_fractions(vtraj) = sum(log_fractions((vtraj-1)*N_VSTEP+1:(vtraj-1)*N_VSTEP+N_VSTEP))/N_VSTEP

        end do 

        !call get_biKDE(vrecord2, incr_kde, btype, kernel, g_array)

        do vtraj=1, N_VTRAJ
          !call get_log_fraction(vrecord2((vtraj-1)*N_VSTEP+1:vtraj*N_VSTEP,:), g_array, path_size(vtraj), log_fraction)
          gauss(1) = gauss(1) + v_noise(1,vtraj)*mean_log_fractions(vtraj)
          gauss(2) = gauss(2) + v_noise(2,vtraj)*mean_log_fractions(vtraj)
        end do

      end if
      
      gauss(1) = gauss(1) / N_VTRAJ
      gauss(2) = gauss(2) / N_VTRAJ

      !print *, step, gauss(1), gauss(2) 


      ! real langevin
      ! f_langevin = -FRICTION*velocities(:,i) + 2*NOI*gauss
      f_langevin(1) = -FRICTION*velocities(1,i) + 2*NOI*gauss(1)
      f_langevin(2) = -FRICTION*velocities(2,i) + 2*NOI*gauss(2)

      ! call soft repulsion
      call repulsions(c_coordinates, c_coordinates(:,i), i, N_AGENT, BOX_LENGTH, &
                      DIAMETER, REPUL_STRENGTH, REPUL_EXP, REPUL_TYPE, f_repul)

      ! total force
      f_total(1) = f_langevin(1) + f_repul(1)
      f_total(2) = f_langevin(2) + f_repul(2)


      ! call integration(TIMESTEP, coordinates(:,i), i_coordinates, &
      !                  p_coordinates(:,i), f_total, velocities(:,i))

      ! integrate
      i_coordinates(1) = 2*coordinates(1,i) - p_coordinates(1,i) + f_total(1)*TIMESTEP**2
      i_coordinates(2) = 2*coordinates(2,i) - p_coordinates(2,i) + f_total(2)*TIMESTEP**2

      ! measure new velocity
      velocities(1,i) = (i_coordinates(1)-coordinates(1,i)) / TIMESTEP
      velocities(2,i) = (i_coordinates(2)-coordinates(2,i)) / TIMESTEP

      ! update past coordinates
      p_coordinates(1,i) = coordinates(1,i)   
      p_coordinates(2,i) = coordinates(2,i)

      coordinates(1,i) = i_coordinates(1)
      coordinates(2,i) = i_coordinates(2)

      c_coordinates(1,i) = modulo(coordinates(1,i), BOX_LENGTH)
      c_coordinates(2,i) = modulo(coordinates(2,i), BOX_LENGTH)

    end do

    do j=N_AGENT_STAT+1, N_AGENT
      write(coor_file_code, *) c_coordinates(1,j),c_coordinates(2,j), &
                                    velocities(1,j), velocities(2,j)
    end do

    !do vtraj=1,N_VTRAJ  ! print virt traj Marta way - Klara did small correction
     ! do vstep=1,N_VSTEP
      !  write(vrecord_file_code,*) vrecord(vtraj,vstep), vrecord(vtraj,(vstep+N_VSTEP-1))
      !end do
    !end do

     !do vtraj=1,N_VTRAJ ! print virt traj Alvin way
      !do vstep=1,N_VSTEP*2
         ! write(vrecord_file_code,*) vrecord(vtraj,vstep)
        !end do
     ! end do

  end do

  close(coor_file_code)
  close(vrecord_file_code)

end program cognitive_maps
