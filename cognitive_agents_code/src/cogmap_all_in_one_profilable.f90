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


  ! ### CONSTANTS ###
  ! Simulation parameters
  ! integer, parameter    :: N_AGENT = 200
  ! integer, parameter    :: N_STEP = 2000
  integer, parameter    :: N_VTRAJ = 20
  integer, parameter    :: N_AGENT_STAT = 0 ! no of static / non-moving agents
  ! integer, parameter    :: N_VSTEP = 20
  real, parameter       :: DIAMETER = 2
  real, parameter       :: BOX_LENGTH = 80
  ! real, parameter       :: BOX_AREA = BOX_LENGTH ** 2.0
  real, parameter       :: TIMESTEP = 0.1
  character, parameter  :: VREPUL_TYPE = "H" ! Hard, Soft, sWarm, Ahc, None, stop Virtual
  character, parameter  :: REPUL_TYPE = "S" ! Hard, Soft, sWarm, Ahc
  real, parameter       :: FRICTION = 1.0 ! gamma
  real, parameter       :: TEMPERATURE = 1.0
  real, parameter       :: NOI = sqrt(2.0*FRICTION*TEMPERATURE/TIMESTEP)
  real, parameter       :: REPUL_STRENGTH = 20.0
  real, parameter       :: REPUL_EXP = 60.0
  ! real, parameter       :: PI = 3.14159265358979323846264338327950288419716939937510582097494459 ! from https://oeis.org/A000796
  ! real, parameter       :: TAU = 6.28318530717958647692528676655900576839433879875021164194988918 ! from https://oeis.org/A019692
  integer, parameter    :: INCR_COOR = 1
  real, parameter       :: empty = 0.00
  character(len=2), parameter  :: CogmapTYPE = 'RG' ! 'LogFraction' 'RadiusGyration'

  ! TODO include analysis parameters

  ! ### VARIABLES ###
  integer           :: i, step, vtraj, vstep, vsteps, d ! iterators
  integer           :: n_agent, n_vstep, n_step !,N_AGENT_STAT ! , n_vtraj change constants to variables
  integer           :: start, incr, vrr
  real              :: r_gyration_norm, stat_radius, log_fraction
  character(len=90) :: path, dir_name, file_name, arg_n_vtraj
  character(len=90) :: arg_n_agent, arg_n_vstep, arg_n_step, arg_stat_r !, arg_repul
  character(len=90) :: box_length_char, n_step_char
  character         :: coor_file_code, vrecord_file_code
  logical           :: stop_traj

  real, dimension(:,:), allocatable :: coordinates
  real, dimension(:,:), allocatable :: p_coordinates, c_coordinates, velocities
  real, dimension(:,:), allocatable :: v_coordinates!, v_velocities
  real, dimension(:,:), allocatable :: rnumb, vrecord
  real, dimension(:,:), allocatable :: v_noise
  real, dimension(:,:), allocatable :: g_array, g_array_coor, vrecord2
  ! real, dimension(4) :: rand4
  real, dimension(2) :: f_langevin, f_repul, f_total, i_coordinates
  real, dimension(2) :: v_p_coordinates, v_c_coordinates, v_velocity
  real, dimension(2) :: gauss, mean_r_virt, r_gyr_comp
  real, dimension(:), allocatable :: r_gyration

  real, parameter :: incr_kde = 0.1
  character(len=90), parameter :: btype = 'S', kernel = 'U'

  ! ###SYNTHETIC###
  ! FIX THE SEEDS
  integer, allocatable :: seed(:)
  integer seed_size
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed(:) = 10
  call random_seed(put=seed)
  call zigset(10)


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
  read(arg_n_step,*)N_STEP
  ! read(arg_n_vtraj,*)n_vtraj
  ! read(arg_stat_r,*)stat_radius
  stat_radius = 10
  ! read(arg_repul,*)REPUL_EXP
  ! the 'file_name' argument stays a character

  ! directory name variables
  write(box_length_char,*)int(BOX_LENGTH)
  write(n_step_char,*)int(N_STEP)

  ! Size allocation
  allocate(coordinates(N_AGENT,2))
  allocate(p_coordinates(N_AGENT,2))
  allocate(c_coordinates(N_AGENT,2))
  allocate(velocities(N_AGENT,2))

  allocate(v_coordinates(N_VSTEP+1,2))
  allocate(rnumb(N_VSTEP+1,2))
  allocate(v_noise(N_VTRAJ,2))
  allocate(r_gyration(N_VTRAJ))
  allocate(vrecord(N_VTRAJ,N_VSTEP*2))

  allocate(vrecord2(N_VSTEP*N_VTRAJ,2))

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
  ! coordinates = place_agents(N_AGENT, stat_radius) ! for set initial coor agents SYNTHETIC

  ! SYNTHETIC 2 HORIZONTAL
  ! coordinates(1,:) = [38,40]
  ! coordinates(2,:) = [42,40]
  ! SYNTHETIC 2 HORIZONTAL 2
  ! coordinates(1,:) = [02,40]
  ! coordinates(2,:) = [78,40]
  ! SYNTHETIC 2 VERTICAL
  ! coordinates(1,:) = [40,38]
  ! coordinates(2,:) = [40,42]
  ! SYNTHETIC 2 DIAGONAL
  ! coordinates(1,:) = [38,38]
  ! coordinates(2,:) = [42,42]
  ! SYNTHETIC 3 TRIANGLE
  ! coordinates(1,:) = [38,39]
  ! coordinates(2,:) = [42,39]
  ! coordinates(2,:) = [40,42]

  vrr = N_AGENT  ! virtual traj record particle number (<np)
  p_coordinates = coordinates ! past coordinates
  c_coordinates = coordinates ! corrected coordinates

  do i=1, N_AGENT
    velocities(i,:) = [rnor(), rnor()] ! TURN OFF IN SYNTHETIC
    ! velocities(i,:) = 0
  end do

  do i=1, n_agent  ! print initial coordinates
    write(ichar(coor_file_code), *) c_coordinates(i,1),c_coordinates(i,2), &
        velocities(i,1), velocities(i,2)
  end do

  ! Main loop
  do step=1, N_STEP ! for every step in time

    ! ! used to print coord for non moving agents
     do j=1, N_AGENT_STAT
       write(ichar(coor_file_code), *) coordinates(j,1),coordinates(j,2), &
           empty, empty
     end do

    do i=N_AGENT_STAT+1, N_AGENT ! for every particle in that timestep

      vrecord2(:,:) = 0.0

      do vtraj=1, N_VTRAJ ! for every virt traj for that particle in that timestep

        ! Setup for the first virtual step
        v_coordinates(1,:) = coordinates(i,:)
        v_c_coordinates = c_coordinates(i,:)
        v_p_coordinates = p_coordinates(i,:)
        v_velocity = velocities(i,:)

        ! Generate random number for the whole vtraj
        do vstep=1, N_VSTEP ! for every step in virtual trajectory
          rnumb(vstep,:) = [rnor(), rnor()]
        end do

        ! Save the first vstep random numbers
        v_noise(vtraj,:) = rnumb(1,:)
        vsteps = N_VSTEP + 1

        do vstep=1, N_VSTEP ! for every step in virtual trajectory

          if(VREPUL_TYPE.ne.'V')then
            call repulsions(c_coordinates, v_c_coordinates, i, N_AGENT, BOX_LENGTH, &
                            DIAMETER, REPUL_STRENGTH, REPUL_EXP, VREPUL_TYPE, f_repul)
          end if

          if(VREPUL_TYPE.eq.'V')then
            stop_traj = virtual_repulsion(c_coordinates, i, v_c_coordinates, BOX_LENGTH, &
                                   DIAMETER)

            if(stop_traj)then
              vsteps = vstep
              ! THIS IS A HACK
              v_coordinates(vstep:(N_VSTEP+1),1) = i_coordinates(1)
              v_coordinates(vstep:(N_VSTEP+1),2) = i_coordinates(2)
              exit
            end if

          end if        ! end if for equals 'V' check

          ! virtual langevin equation
          f_langevin = -FRICTION*v_velocity + NOI*rnumb(vstep,:)

          ! virtual total force
          f_total = f_langevin + f_repul

          ! virtual integration
          call integration(TIMESTEP, v_coordinates(vstep,:), i_coordinates, &
                           v_p_coordinates, f_total, v_velocity)

          v_coordinates(vstep+1,:) = i_coordinates
          v_c_coordinates = modulo(v_coordinates(vstep,:), BOX_LENGTH)

        end do

        if (CogmapTYPE=='RG') then

          ! radius of gyration for trajectory i
          mean_r_virt = [ ((sum(v_coordinates(1:vsteps,d)) / vsteps) , d=1, 2) ]

          r_gyr_comp = 0

          do vstep=1, vsteps
            r_gyr_comp = r_gyr_comp+(v_coordinates(vstep,:)-mean_r_virt)**2
          end do

          r_gyr_comp = r_gyr_comp / vsteps
          r_gyration(vtraj) = sqrt(sum(r_gyr_comp))

        end if

        if (CogmapTYPE=='LF') then

          do vstep=1,N_VSTEP
            vrecord2((vtraj-1)*N_VSTEP+vstep,1)=MODULO(v_coordinates(vstep+1,1),BOX_LENGTH)
            vrecord2((vtraj-1)*N_VSTEP+vstep,2)=MODULO(v_coordinates(vstep+1,2),BOX_LENGTH)
          end do

         ! print *, 'vrecord2' ,vrecord2

        end if

        if (i.eq.vrr) then
          do vstep=1,N_VSTEP
            vrecord(vtraj,vstep)= modulo(v_coordinates(vstep+1,1),BOX_LENGTH)
            vrecord(vtraj,N_VSTEP+vstep)= modulo(v_coordinates(vstep+1,2),BOX_LENGTH)
          end do
        end if

      end do
      ! choosing which way to go

      gauss = 0

      if (CogmapTYPE=='RG') then

        do vtraj=1, N_VTRAJ
        if (r_gyration(vtraj).ne.0.0)then
          r_gyration_norm = log(r_gyration(vtraj) / (sum(r_gyration)/size(r_gyration)))
          gauss = gauss + v_noise(vtraj,:)*r_gyration_norm
        end if
      end do


      end if

      
      gauss = gauss / N_VTRAJ


      ! real langevin
      f_langevin = -FRICTION*velocities(i,:) + 2*NOI*gauss

      ! call soft repulsion
      call repulsions(c_coordinates, c_coordinates(i,:), i, N_AGENT, BOX_LENGTH, &
                      DIAMETER, REPUL_STRENGTH, REPUL_EXP, REPUL_TYPE, f_repul )

      ! total force
      f_total = f_langevin + f_repul


      call integration(TIMESTEP, coordinates(i,:), i_coordinates, &
                       p_coordinates(i,:), f_total, velocities(i,:))

      coordinates(i,:) = i_coordinates

      c_coordinates(i,:) = modulo(coordinates(i,:), BOX_LENGTH)

      if (modulo(step,INCR_COOR).eq.0)then
        write(ichar(coor_file_code), *) c_coordinates(i,1),c_coordinates(i,2), &
                                        velocities(i,1), velocities(i,2)
      end if

    end do

    !do vtraj=1,N_VTRAJ  ! print virt traj Marta way - Klara did small correction
     ! do vstep=1,N_VSTEP
      !  write(ichar(vrecord_file_code),*) vrecord(vtraj,vstep), vrecord(vtraj,(vstep+N_VSTEP-1))
      !end do
    !end do

     do vtraj=1,N_VTRAJ ! print virt traj Alvin way
      do vstep=1,N_VSTEP*2
         write(ichar(vrecord_file_code),*) vrecord(vtraj,vstep)
       end do
     end do

  end do

  contains

    subroutine init_files(path,coor_file_code, vrecord_file_code)

      ! INPUT
      character(len=90), intent(in) :: path
      ! OUTPUT
      character, intent(out) :: coor_file_code, vrecord_file_code ! just 2 for nwo
      character(len=90) :: coor_file, vrecord_file ! just 2 for nwo

      coor_file_code = 'c'
      vrecord_file_code = 'v'

      coor_file = trim(path)//"_coor.dat"
      vrecord_file = trim(path)//"_vrecord.dat"

      open(ichar(coor_file_code),file=coor_file,status="new", action="write")
      open(ichar(vrecord_file_code),file=vrecord_file,status="new", action="write")

    end subroutine init_files

    subroutine init_directory(arg_n_agent,arg_n_vstep,box_length_char,&
                              REPUL_TYPE,n_step_char,dir_name)
      ! INPUTS
      character(len=90), intent(in)  :: arg_n_agent, arg_n_vstep
      character(len=90), intent(in)  :: n_step_char, box_length_char
      character, intent(in)          :: REPUL_TYPE
      ! INTERNALS
      logical :: d_exists
      ! OUTPUT
      character(len=90), intent(out) :: dir_name

      dir_name=trim(adjustl("na"))//trim(adjustl(arg_n_agent))// &
      trim(adjustl("_vs"))//trim(adjustl(arg_n_vstep))// &
      trim(adjustl("_st"))//trim(adjustl(n_step_char))// &
      trim(adjustl("_lb"))//trim(adjustl(box_length_char))// &
      trim(adjustl("_rt"))//trim(adjustl(REPUL_TYPE))

      inquire(file=dir_name, exist=d_exists)

      if (.not.d_exists) then
        call system('mkdir '//dir_name)
      end if

    end subroutine init_directory

    subroutine write_info_file(path,nsteps,dt,FRICTION,temp,start,incr,incr_coor)

      ! INPUT
      character(len=90), intent(in) :: path
      integer, intent(in)           :: nsteps,start,incr,incr_coor
      real, intent(in)              :: dt, friction, temp
      ! INTERNALS
      character :: info_file
      character(len=90) :: info_file_name

      info_file_name = trim(path)//"_info.dat"

      open(ichar(info_file),file=info_file_name,status='new',action='write')
      write(ichar(info_file),*)'nsteps',nsteps
      write(ichar(info_file),*)'dt',dt
      write(ichar(info_file),*)'friction',friction
      write(ichar(info_file),*)'temp', temp
      write(ichar(info_file),*)'start',start
      write(ichar(info_file),*)'incr', incr
      write(ichar(info_file),*)'incr_coor',incr_coor
      close(ichar(info_file))

    end subroutine write_info_file

    pure function ext_coor(i, coordinates, BOX_LENGTH)

      ! INPUTS
      integer, intent(in) :: i
      real, dimension(:,:), intent(in) :: coordinates
      real, intent(in) :: BOX_LENGTH
      ! INTERNALS
      integer :: j, n, m
      real, dimension(:), allocatable :: coor_i
      ! OUTPUT
      real, dimension(:,:), allocatable :: ext_coor

      n = size(coordinates, 1)
      m = size(coordinates, 2)

      allocate(coor_i(m))

      coor_i = coordinates(i,:)

      allocate(ext_coor(4*n, m))

      ! transform 0 : copy coordinates four times
      do j=1, 4
        ext_coor((j-1)*n+1:j*n,:) = coordinates
      end do

      ! transform 1 : transform in x on copies 2 and 3
      ext_coor(n+1:n*3,1) = ext_coor(n+1:n*3,1) - sign(BOX_LENGTH, ext_coor(n+1:n*3,1)-coor_i(1))
      ! ext_coor = reshape([ext_coor, ext_x_coor], [size(ext_coor, 1)*2, size(ext_coor, 2)])

      ! transform 2 : transform in y on copies 3 and 4
      ext_coor(n*2+1:n*4,2) = ext_coor(n*2+1:n*4,2) - sign(BOX_LENGTH, ext_coor(n*2+1:n*4,2)-coor_i(2))

      ! ext_coor contains : [1] - original coordinates
      !                     [2] - transform x
      !                     [3] - transform x and y
      !                     [4] - transform y

      deallocate(coor_i)

    end function ext_coor


    subroutine init_random_seed()

      integer :: seed
      integer, dimension(8) :: time

      call date_and_time(VALUES=time)

      seed = time(4) + (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))

      call zigset(seed)

    end subroutine


    ! verlet integration
    subroutine integration(TIMESTEP, agent_coordinates, &
                           i_coordinates, agent_p_coordinates, &
                           f_total, velocity)

      ! INPUTS
      real, intent(in) :: TIMESTEP
      real, dimension(2), intent(in) :: f_total, agent_coordinates
      ! OUTPUTS
      real, dimension(2), intent(out) :: i_coordinates, agent_p_coordinates, velocity

      ! integrate
      i_coordinates = 2*agent_coordinates - agent_p_coordinates + f_total*TIMESTEP**2

      ! measure new velocity
      velocity = (i_coordinates-agent_coordinates) / TIMESTEP

      ! update past coordinates
      agent_p_coordinates = agent_coordinates

    end subroutine integration


    pure function distance(coor_i, coor_j, BOX_LENGTH)
      ! INPUTS
      real, dimension(2), intent(in) :: coor_i, coor_j
      real, intent(in)  :: BOX_LENGTH
      ! INTERNAL
      real              :: box, half_box_length
      ! OUPUT
      real, dimension(2) :: distance

      box = BOX_LENGTH
      half_box_length = box*0.5 ! to avoid calculating every time

      distance = coor_i - coor_j

      if (distance(1).gt.half_box_length) then
        distance(1) = distance(1) - box
      else if (distance(1).lt.-half_box_length) then
        distance(1) = distance(1) + box
      end if

      if (distance(2).gt.half_box_length) then
        distance(2) = distance(2) - box
      else if (distance(2).lt.-half_box_length) then
        distance(2) = distance(2) + box
      end if

      ! where (distance.lt.-half_box_length)
        ! distance = distance + box
      ! elsewhere (distance.gt.half_box_length)
        ! distance = distance - box
      ! end where

    end function distance


    function init_agents(N_AGENT, DIAMETER, BOX_LENGTH) result(agent_coor)

      ! INPUTS
      integer, intent(in)   :: N_AGENT
      real, intent(in)      :: DIAMETER, BOX_LENGTH
      ! INTERNALS
      real                  :: min_distance, allowed_length, radius, box, dia
      real                  :: theor_max, slow_regime ! safety checks
      real                  :: start_time, current_time ! timekeeping
      integer               :: i, j ! iterator
      real, dimension(2)    :: rand2, coor_vector, dist_vector
      ! OUTPUT
      real, dimension(N_AGENT,2) :: agent_coor

      box = BOX_LENGTH
      dia = DIAMETER
      radius = 0.5 * dia
      min_distance = 1.1 * dia
      allowed_length = box - dia

      ! Sanity checks
      theor_max = int(box*0.5)**2
      if (N_AGENT.ge.theor_max) then
        print *, "Theoretical maximum limit exceeded:", N_AGENT, "/", theor_max
        print *, "TERMINATING THE PROGRAM"
        stop
      end if

      slow_regime = theor_max * 2/3
      if (N_AGENT.ge.slow_regime) then
        print *, "Slow regime limit exceeded:", N_AGENT, "/",slow_regime
      end if

      ! Timer
      call cpu_time(start_time)

      ! Agent placement

      ! place first candidate
      call random_number(rand2)
      agent_coor(1,:) = allowed_length*rand2 + radius

      do i=2, N_AGENT ! for each agent

        candidate: do ! keep on creating candidates until exited

          ! create random x and y coordinates for the candidate agent
          call random_number(rand2)
          coor_vector = allowed_length*rand2 + radius

          ! calculate distance to each other agents
          do j=1, i
            if (i.eq.j) cycle ! not for the agent itself
            dist_vector = distance(coor_vector, agent_coor(j,:), box)
            if (norm2(dist_vector).le.MIN_DISTANCE) cycle candidate
          end do

          ! agent is not first and far from others, accept candidate
          agent_coor(i,:) = coor_vector

          ! if placing agents takes too long, stop program
          call cpu_time(current_time)
          if (current_time-start_time.ge.600) then
            print *, "Could not place agents in 600 seconds"
            print *, "TERMINATING THE PROGRAM"
            stop
          end if

          exit candidate

        end do candidate
      end do

    end function init_agents


    function place_agents(N_AGENT,radius) result(agent_coor)

      ! INPUTS
      integer, intent(in)   :: N_AGENT
      real, intent(in)      :: radius
      ! INTERNALS
      integer               :: i
      ! ! circle radius 30
      ! real,dimension(27)::x(1:27)=[70.00,68.19,62.98,55.00,45.21,34.79,&
      ! 25.00,17.02,11.81,10.00,11.81,17.02,25.00,34.79,45.21,55.00,62.98,&
      ! 68.19,40.00,42.00,38.00,40.00,42.00,38.00,40.00,42.00,38.00]
      ! real,dimension(27)::y(1:27)=[40.00,50.26,59.28,65.98,69.54,69.54,&
      ! 65.98,59.28,50.26,40.00,29.74,20.72,14.02,10.46,10.46,14.02,20.72,&
      ! 29.74,40.00,40.00,40.00,38.00,38.00,38.00,42.00,42.00,42.00]

       ! circle radius 10
       !real,dimension(28)::x(1:28)=[50.00,49.40,47.66,45.00,41.74,38.26,&
       !35.00,32.34,30.60,30.00,30.60,32.34,35.00,38.26,41.74,45.00,47.66,&
       !49.40,40.00,42.00,38.00,40.00,42.00,38.00,40.00,42.00,38.00, 40.00]
       !real,dimension(28)::y(1:28)=[40.00,43.42,46.43,48.66,49.85,49.85,&
       !48.66,46.43,43.42,40.00,36.58,33.57,31.34,30.15,30.15,31.34,33.57,&
       !36.58,40.00,40.00,40.00,38.00,38.00,38.00,42.00,42.00,42.00, 40.00]

      ! ! circle radius 5
      ! real,dimension(27)::x(1:27)=[45.00,44.70,43.83,42.50,40.87,39.13,&
      ! 37.50,36.17,35.30,35.00,35.30,36.17,37.50,39.13,40.87,42.50,43.83,&
      ! 44.70,40.00,42.00,38.00,40.00,42.00,38.00,40.00,42.00,38.00]
      ! real,dimension(27)::y(1:27)=[40.00,41.71,43.21,44.33,44.92,44.92,&
      ! 44.33,43.21,41.71,40.00,38.29,36.79,35.67,35.08,35.08,35.67,36.79,&
      ! 38.29,40.00,40.00,40.00,38.00,38.00,38.00,42.00,42.00,42.00]

     !  ! circle radius 5 + 1 dot
     ! real,dimension(19)::x(1:19)=(/45.00,44.70,43.83,42.50,40.87,39.13,&
     ! 37.50,36.17,35.30,35.00,35.30,36.17,37.50,39.13,40.87,42.50,43.83,&
     ! 44.70,40.00/)
     ! real,dimension(19)::y(1:19)=(/40.00,41.71,43.21,44.33,44.92,44.92,&
     ! 44.33,43.21,41.71,40.00,38.29,36.79,35.67,35.08,35.08,35.67,36.79,&
     ! 38.29,40.00/)

     ! 9 dot grid
     !real,dimension(9)::x(1:9)=(/38.00,40.00,42.00,38.00,40.00,42.00,&
     !38.00,40.00,42.00/)
     !real,dimension(9)::y(1:9)=(/42.00,42.00,42.00,40.00,40.00,40.00,&
     !38.00,38.00,38.00/)
     
     ! ! line + dot parallel
      !real,dimension(8)::x(1:8)=(/38.00,38.00,38.00,38.00,38.00,38.00,&
      !38.00,40.00/)
      !real,dimension(8)::y(1:8)=(/40.00,42.00,38.00,44.00,36.00,46.00,&
      !34.00,40.00/)

     ! ! line + dot in line
     ! real,dimension(8)::x(1:8)=(/40.00,42.00,38.00,44.00,36.00,46.00,&
     ! 34.00,32.00/)
     ! real,dimension(8)::y(1:8)=(/40.00,40.00,40.00,40.00,40.00,40.00,&
     ! 40.00,40.00/)

     ! ! 2 parallel lines + dot
     ! real,dimension(15)::x(1:15)=(/35.00,35.00,35.00,35.00,35.00,35.00,&
     ! 35.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,40.00/)
     ! real,dimension(15)::y(1:15)=(/40.00,42.00,38.00,44.00,36.00,46.00,&
     ! 34.00,40.00,42.00,38.00,44.00,36.00,46.00,34.00,40.00/)

     ! ! 2 parallel lines + dot     closer together
     ! real,dimension(15)::x(1:15)=(/37.00,37.00,37.00,37.00,37.00,37.00,&
     ! 37.00,43.00,43.00,43.00,43.00,43.00,43.00,43.00,40.00/)
     ! real,dimension(15)::y(1:15)=(/40.00,42.00,38.00,44.00,36.00,46.00,&
     ! 34.00,40.00,42.00,38.00,44.00,36.00,46.00,34.00,40.00/)

     ! ! 2 parallel lines + dot
     ! real,dimension(43)::x(1:43)=(/35.00,35.00,35.00,35.00,35.00,35.00,&
     ! 35.00,35.00,35.00,35.00,35.00,35.00,35.00,35.00,35.00,35.00,35.00,&
     ! 35.00,35.00,35.00,35.00,&
     ! 45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,&
     ! 45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,45.00,&
     ! 40.00/)
     ! real,dimension(43)::y(1:43)=(/40.00,42.00,44.00,46.00,48.00,50.00,&
     ! 52.00,54.00,56.00,58.00,60.00,38.00,36.00,34.00,32.00,30.00,28.00,&
     ! 26.00,24.00,22.00,20.00,&
     ! 40.00,42.00,44.00,46.00,48.00,50.00,52.00,54.00,56.00,58.00,60.00,&
     ! 38.00,36.00,34.00,32.00,30.00,28.00,26.00,24.00,22.00,20.00,&
     ! 40.00/)

     ! a square box + 1 particle starting in left lower corner (33.00) 
       real,dimension(41)::x(1:41)=(/30.00,30.00,30.00,30.00,30.00,30.00,&
       30.00,30.00,30.00,30.00,30.00,32.00,34.00,36.00,38.00,40.00,42.00,&
       44.00,46.00,48.00,50.00,32.00,34.00,36.00,38.00,40.00,42.00,44.00,&
       46.00,48.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,&
       50.00,33.00/)
       real,dimension(41)::y(1:41)=(/30.00,32.00,34.00,36.00,38.00,40.00,&
       42.00,44.00,46.00,48.00,50.00,30.00,30.00,30.00,30.00,30.00,30.00,&
       30.00,30.00,30.00,30.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,&
       50.00,50.00,50.00,32.00,34.00,36.00,38.00,40.00,42.00,44.00,46.00,&
       48.00,33.00/)

      real                :: angle
      ! real, parameter     :: radius = 5

      ! OUTPUT
      real, dimension(N_AGENT,2) :: agent_coor

    ! Agent placement

      ! ! for circle + 1 moving agent
      !
      ! angle = 6.28319 / (N_AGENT-1)
      ! do i=1, N_AGENT-1 ! for each agent
      !   agent_coor(i,:) = [(cos(angle*(i-1))*radius)+40, (sin(angle*(i-1))*radius)+40]
      ! end do
      ! agent_coor(N_AGENT,:) = [40.00, 40.00]

      ! do i=1, N_AGENT ! for each agent
      !   print*, agent_coor(i,:)
      ! end do
      !
      do i=1, N_AGENT ! for each agent
        agent_coor(i,:) = [x(i), y(i)]
      end do

      ! SYNTHETIC 2 HORIZONTAL
      ! coordinates(1,:) = [38,40]
      ! coordinates(2,:) = [42,40]
      ! SYNTHETIC 2 VERTICAL
      ! coordinates(1,:) = [40,38]
      ! coordinates(2,:) = [40,42]
      ! SYNTHETIC 3 TRIANGLE
      ! coordinates(1,:) = [38,39]
      ! coordinates(2,:) = [42,39]
      ! coordinates(2,:) = [40,42]

    end function place_agents


    subroutine repulsions(coordinates, coor_i, i, N_AGENT, BOX_LENGTH, DIAMETER, &
                          REPUL_STRENGTH, REPUL_EXP, REPUL_TYPE, repulsion)

      ! INPUTS
      real, dimension(N_AGENT,2), intent(in) :: coordinates
      real, dimension(2), intent(in) :: coor_i
      integer, intent(in)     :: i, N_AGENT
      real, intent(in)        :: BOX_LENGTH, DIAMETER
      real, intent(in)        :: REPUL_STRENGTH, REPUL_EXP
      character, intent(in)   :: REPUL_TYPE
      ! INTERNALS
      integer                 :: j
      real                    :: box, half_box_length, dia
      real                    :: dist, potential
      real, dimension(2)      :: dist_vector, unit_vector
      ! OUPUTS
      real, dimension(2), intent(out) :: repulsion

      box = BOX_LENGTH
      dia = DIAMETER
      half_box_length = box*0.5 ! to avoid calculating every time

      repulsion = 0

      do j=1, N_AGENT

        dist_vector = coor_i - coordinates(j,:)

        if (dist_vector(1).lt.-half_box_length) then
          dist_vector(1) = dist_vector(1) + box
        else if (dist_vector(1).gt.half_box_length) then
          dist_vector(1) = dist_vector(1) - box
        end if

        if (dist_vector(2).lt.-half_box_length) then
          dist_vector(2) = dist_vector(2) + box
        else if (dist_vector(2).gt.half_box_length) then
          dist_vector(2) = dist_vector(2) - box
        end if

        ! where (dist_vector.lt.-half_box_length)
        !   dist_vector = dist_vector + box
        ! elsewhere (dist_vector.gt.half_box_length)
        !   dist_vector = dist_vector - box
        ! end where

        ! dist_vector = distance(coor_i, coordinates(j,:), box)
        ! dist = norm2(dist_vector)
        dist = hypot(dist_vector(1), dist_vector(2))
        ! dist = sqrt(dist_vector(1)**2 + dist_vector(2)**2)

        if (dist.le.dia) then ! if contact between the agents

          if (j .eq. i) cycle

          ! decompose the dist into x and y components
          where (dist_vector.ne.0.0)
            unit_vector = dist_vector / dist
          elsewhere
            unit_vector = 0
          endwhere

          ! get the potential for type of repulsion
          select case (REPUL_TYPE) ! <S>oftcore, <H>ardcore, <N>ull, <A>djustable
            case ("S")
              potential = 1 - dist/dia
            case ("H")
              potential = 1
            case ("N") ! for testing purposes
              potential = 0
            case ("A")
              potential = (exp(REPUL_EXP*dist/dia)-exp(REPUL_EXP))/&
                          (1-exp(REPUL_EXP))
            ! Adjustable gives Soft for REPUL_EXP = 0, Hard for REPUL_EXP = +Inf
          end select

          repulsion = repulsion + REPUL_STRENGTH*unit_vector*potential

        end if

      end do

    end subroutine repulsions

    pure function virtual_repulsion(coordinates, i, coor_i, BOX_LENGTH, &
                                    DIAMETER) result(stop_traj)
      ! INPUTS
      real, dimension(:,:), intent(in)  :: coordinates
      integer, intent(in)               :: i
      real, dimension(:), intent(in)    :: coor_i
      real, intent(in)                  :: BOX_LENGTH, DIAMETER
      ! INTERNALS
      integer   :: j
      real      :: dist
      real, dimension(2)      :: dist_vector
      ! OUTPUTS
      logical   :: stop_traj

      stop_traj = .false.

      do j=1, size(coordinates, 1)
    
        dist_vector = distance(coor_i, coordinates(j,:), BOX_LENGTH)
        ! dist = norm2(dist_vector)
        dist = hypot(dist_vector(1), dist_vector(2))

          if (dist.le.DIAMETER) then
            if(i.eq.j) cycle
            stop_traj = .true.
            exit
          end if

      end do

    end function virtual_repulsion
    

end program cognitive_maps
