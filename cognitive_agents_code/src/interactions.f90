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

module interactions

  use utilities, only: distance ! for the distance function

  implicit none

  contains

    subroutine repulsions(coordinates, coor_i, i, N_AGENT, BOX_LENGTH, DIAMETER, &
                          REPUL_STRENGTH, REPUL_EXP, REPUL_TYPE, repulsion)

      ! INPUTS
      real, dimension(2,N_AGENT), intent(in) :: coordinates
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

        dist_vector(1) = coor_i(1) - coordinates(1,j)
        dist_vector(2) = coor_i(2) - coordinates(2,j)

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

        ! dist_vector = distance(coor_i, coordinates(:,j), box)
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

          repulsion(1) = repulsion(1) + REPUL_STRENGTH*unit_vector(1)*potential
          repulsion(2) = repulsion(2) + REPUL_STRENGTH*unit_vector(2)*potential

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

      do j=1, size(coordinates, 2)
    
        dist_vector = distance(coor_i, coordinates(:,j), BOX_LENGTH)
        ! dist = norm2(dist_vector)
        dist = hypot(dist_vector(1), dist_vector(2))

          if (dist.le.DIAMETER) then
            if(i.eq.j) cycle
            stop_traj = .true.
            exit
          end if

      end do

    end function virtual_repulsion
    
end module interactions
