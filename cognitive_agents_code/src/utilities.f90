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

module utilities

  use ziggurat

  implicit none

  contains


    pure function ext_coor(i, coordinates, BOX_LENGTH, N_AGENT)

      ! INPUTS
      integer, intent(in) :: i, N_AGENT
      real, dimension(2,N_AGENT), intent(in) :: coordinates
      real, intent(in) :: BOX_LENGTH
      ! INTERNALS
      integer :: j, n
      real, dimension(2) :: coor_i
      ! OUTPUT
      real, dimension(2,4*N_AGENT) :: ext_coor

      coor_i = coordinates(:,i)
      n = N_AGENT

      ! transform 0 : copy coordinates four times
      do j=1, 4
        ext_coor(:,(j-1)*n+1:j*n) = coordinates
      end do

      ! transform 1 : transform in x on copies 2 and 3
      ext_coor(1,n+1:n*3) = ext_coor(1,n+1:n*3) - sign(BOX_LENGTH, ext_coor(1,n+1:n*3)-coor_i(1))
      ! ext_coor = reshape([ext_coor, ext_x_coor], [size(ext_coor, 1)*2, size(ext_coor, 2)])

      ! transform 2 : transform in y on copies 3 and 4
      ext_coor(2,n*2+1:n*4) = ext_coor(2,n*2+1:n*4) - sign(BOX_LENGTH, ext_coor(2,n*2+1:n*4)-coor_i(2))

      ! ext_coor contains : [1] - original coordinates
      !                     [2] - transform x
      !                     [3] - transform x and y
      !                     [4] - transform y

    end function ext_coor


    subroutine init_random_seed()
      INTEGER :: i, n, clock, zigseed
      INTEGER, DIMENSION(:), ALLOCATABLE :: ranseed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(ranseed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      ranseed=abs((clock+37*[ (i - 1, i = 1, n) ])**2)
      CALL RANDOM_SEED(PUT = ranseed)
      zigseed=sum(ranseed)
      call zigset(zigseed)

      DEALLOCATE(ranseed)
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
      i_coordinates(1) = 2*agent_coordinates(1) - agent_p_coordinates(1) + f_total(1)*TIMESTEP**2
      i_coordinates(2) = 2*agent_coordinates(2) - agent_p_coordinates(2) + f_total(2)*TIMESTEP**2

      ! measure new velocity
      velocity = (i_coordinates(1)-agent_coordinates(1)) / TIMESTEP
      velocity = (i_coordinates(2)-agent_coordinates(2)) / TIMESTEP

      ! update past coordinates
      agent_p_coordinates(1) = agent_coordinates(1)
      agent_p_coordinates(2) = agent_coordinates(2)

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
      real, dimension(2,N_AGENT) :: agent_coor

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
      agent_coor(:,1) = allowed_length*rand2 + radius

      do i=2, N_AGENT ! for each agent

        candidate: do ! keep on creating candidates until exited

          ! create random x and y coordinates for the candidate agent
          call random_number(rand2)
          coor_vector = allowed_length*rand2 + radius

          ! calculate distance to each other agents
          do j=1, i
            if (i.eq.j) cycle ! not for the agent itself
            dist_vector = distance(coor_vector, agent_coor(:,j), box)
            if (norm2(dist_vector).le.MIN_DISTANCE) cycle candidate
          end do

          ! agent is not first and far from others, accept candidate
          agent_coor(:,i) = coor_vector

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
     real,dimension(9)::x(1:9)=(/38.00,40.00,42.00,38.00,40.00,42.00,&
     38.00,40.00,42.00/)
     real,dimension(9)::y(1:9)=(/42.00,42.00,42.00,40.00,40.00,40.00,&
     38.00,38.00,38.00/)
     
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
       ! real,dimension(41)::x(1:41)=(/30.00,30.00,30.00,30.00,30.00,30.00,&
       ! 30.00,30.00,30.00,30.00,30.00,32.00,34.00,36.00,38.00,40.00,42.00,&
       ! 44.00,46.00,48.00,50.00,32.00,34.00,36.00,38.00,40.00,42.00,44.00,&
       ! 46.00,48.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,&
       ! 50.00,33.00/)
       ! real,dimension(41)::y(1:41)=(/30.00,32.00,34.00,36.00,38.00,40.00,&
       ! 42.00,44.00,46.00,48.00,50.00,30.00,30.00,30.00,30.00,30.00,30.00,&
       ! 30.00,30.00,30.00,30.00,50.00,50.00,50.00,50.00,50.00,50.00,50.00,&
       ! 50.00,50.00,50.00,32.00,34.00,36.00,38.00,40.00,42.00,44.00,46.00,&
       ! 48.00,33.00/)

      real                :: angle
      ! real, parameter     :: radius = 5

      ! OUTPUT
      real, dimension(2,N_AGENT) :: agent_coor

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
        agent_coor(1,i) = x(i)
        agent_coor(2,i) = y(i)
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

end module utilities
