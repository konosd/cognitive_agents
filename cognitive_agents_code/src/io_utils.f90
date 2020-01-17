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

module io_utils
  
  implicit none

  contains

    subroutine init_files(path,coor_file_code, vrecord_file_code)

      ! INPUT
      character(len=90), intent(in) :: path
      integer, intent(in) :: coor_file_code, vrecord_file_code ! just 2 for nwo
      ! OUTPUT
      character(len=90) :: coor_file, vrecord_file ! just 2 for nwo

      coor_file = trim(path)//"_coor.dat"
      vrecord_file = trim(path)//"_vrecord.dat"

      open(coor_file_code, file=coor_file, status="new", action="write")
      open(vrecord_file_code, file=vrecord_file, status="new", action="write")

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
      integer :: info_file = 73
      character(len=90) :: info_file_name

      info_file_name = trim(path)//"_info.dat"

      open(info_file,file=info_file_name,status='new',action='write')
      write(info_file,*)'nsteps',nsteps
      write(info_file,*)'dt',dt
      write(info_file,*)'friction',friction
      write(info_file,*)'temp', temp
      write(info_file,*)'start',start
      write(info_file,*)'incr', incr
      write(info_file,*)'incr_coor',incr_coor
      close(info_file)

    end subroutine write_info_file

end module io_utils
