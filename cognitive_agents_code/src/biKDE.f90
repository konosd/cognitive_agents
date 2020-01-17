! Bivariate Kernel Density Estimation
!
! Copyright (C) 2019 Klara Kaleb [klara.kaleb18@imperial.ac.uk]
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

module biKDE

  use qsort_c_module 
  use stats 

  implicit none

  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510582097494459 ! from https://oeis.org/A000796

   contains

   subroutine get_biKDE(data, incr, btype, kernel, g_array)

    ! basi workflow

   	! INPUT
 		real, intent(in) :: incr
 		character(len=1) :: btype !  'S' only supported for now
    character(len=1) :: kernel ! 'Uniform', 'Normal' or custom 'Alvin'
    real, dimension(:,:) :: data 

 		! INTERNALS
 		real, dimension(2) :: bandwidth, coordinates
 		integer :: i, n

 		! OUTPUT
    real, dimension(:,:), allocatable, intent(out) :: g_array

    n = size(data(:,1))
   	
    ! set bandwidth - can be expanded
		select case (btype) 
            case ("S") 
              call set_bandwidth_s(data, bandwidth)
    end select

    call make_grid_array(data,incr,g_array)

    !print *, size(g_array(:,1))

		call density(data,n,bandwidth,g_array,PI,kernel)

   end subroutine get_biKDE

   subroutine make_grid_array(data, incr, g_array) 

   		! the grid size depends on the dimension of the data
      ! the resolution stays the same

   		! INPUT
   		real, intent(in) :: incr
      real, dimension(:,:), intent(in) :: data

   		! INTERNALS
      real, dimension(:), allocatable :: grid_x, grid_y
      real :: minx, miny, maxx, maxy, value
   		integer :: i, j, z, xsize, ysize
   		
   		! OUTPUT
   		real, dimension(:,:), allocatable, intent(out) :: g_array


   		minx = MINVAL(data(:,1))
   		miny = MINVAL(data(:,2))

   		maxx = MAXVAL(data(:,1))
   		maxy = MAXVAL(data(:,2))

      xsize = ceiling((maxx-minx)/incr) + 10  ! padding
      ysize = ceiling((maxy-miny)/incr) + 10  ! padding

      allocate(grid_x(xsize))
      allocate(grid_y(ysize))

      ! xgrid
      value = minx - 2*incr

      i = 1

      do while (value.le.(maxx+10*incr).and.i.le.xsize)
        grid_x(i) = value
        value = value + incr
        i = i + 1
      end do

      ! ygrid
      value = miny - 2*incr

      i = 1

      do while (value.le.(maxy+10*incr).and.i.le.ysize)
        grid_y(i) = value
        value = value + incr
        i = i + 1
      end do

      allocate(g_array(xsize*ysize,3))

      !make grid
      do i=1, xsize

        g_array(((i-1)*ysize+1):i*ysize,1) = grid_x(i)

        do j=1, ysize

          g_array(((i-1)*ysize+j),2) = grid_y(j)

        end do

      end do
        
      g_array(:,3) = 0.0

   end subroutine make_grid_array

   subroutine set_bandwidth_s(data, bandwidth)

    ! setting the bandwidth based on the Silverman's rule of thumb
   	! code translated from https://github.com/Neojume/pythonABC/blob/master/hselect.py

    ! INPUT
   	real, dimension(:,:), intent(in) :: data

   	! INTERNALS
   	real, dimension(2) :: iqr, A, sd
   	integer :: n 
    real, dimension(:,:), allocatable :: sorted_data
  
   	!OUTPUT
   	real, dimension(2), intent(out) :: bandwidth

   	n = size(data(:,1))

    allocate(sorted_data(n,2))
    sorted_data = data

    ! get IQR, which makes it less sensitive to outliers
    call QsortC(sorted_data(:,1))
    call QsortC(sorted_data(:,2))

 	  iqr(1) = sorted_data(floor(n*0.75),1) - sorted_data(floor(n*0.25),1)
   	iqr(2) = sorted_data(floor(n*0.75),2) - sorted_data(floor(n*0.25),2)

   	! get standard deviation

   	call get_sd(n,data,sd)

  	A(1) = MINVAL((/sd(1), (iqr(1)/1.349)/))
    A(2) = MINVAL((/sd(2), (iqr(2)/1.349)/))

   ! 0.9 better than 1.06 at handling multimodality
	  bandwidth = 0.9 * A * n ** (-0.2) 

   end subroutine set_bandwidth_s

   subroutine get_sd(n,data,sd)

   ! abstraction of the stats module functions

   ! INPUT 
   integer, intent(in) :: n
   real, dimension(:,:), intent(in) :: data

   ! INTERNALS
   real, dimension(2) :: ss, ss2
   real :: mean_x, var_x, std_x, mean_y, var_y, std_y
   integer :: i

   !OUTPUT
   real, dimension(2) ,intent(out) :: sd

    ss(:)   = 0.0 
    ss2(:)   = 0.0

    do i=1,n
    	call Sums(data(i,1), ss(1), ss2(1))
    	call Sums(data(i,2), ss(2), ss2(2))
    end do

    call Results(ss(1), ss2(1), n, mean_x, var_x, std_x)
    call Results(ss(2), ss2(2), n, mean_y, var_y, std_y)

    sd(1) = std_x
    sd(2) = std_y

   end subroutine get_sd

   subroutine density(data,n,bandwidth,g_array,PI,kernel)

    ! INPUT
    real, dimension(:,:), intent(in) :: data
    real, dimension(2), intent(in) :: bandwidth
    real, intent(in) :: PI
    integer, intent(in) :: n
    character(len=1), intent(in) :: kernel


    ! INTERNALS
    integer :: i, j, x, y, index, new_index, x_index, y_index, ysize
    real :: k_sum, ux, uy, guy, gux, rux, ruy, x_point
    real,dimension(2) :: xrange, yrange, data_point
    real,dimension(3) :: grid_point
    integer, dimension(4) :: indices


    ! OUTPUT
    real, dimension(:,:), intent(out) :: g_array

    if(kernel.eq.'A') then ! kernel used in the COGMAP project

      do i=1, size(data(:,1))

        data_point = data(i,:)

        ! first occurance of x
        x_index = MINLOC(abs(g_array(:,1)-data_point(1)),dim = 1)

        ! first occurance of y - index of y relative to x
        y_index = MINLOC(abs(g_array(:,2)-data_point(2)),dim = 1)

        index = x_index + y_index - 1

        indices(1) = index

        ! get xsize
        x_point = g_array(x_index,1)
        ysize = count((g_array(:,1)-x_point).eq.0.0)

       ! x axis
       if(g_array(index,1).lt.data_point(1)) then
          indices(2) = index + ysize
          x = 1
        else 
          indices(2) = index - ysize
          x = -1
       end if

       ! y axis
       if(g_array(index,2).lt.data_point(2)) then
         indices(3) = index + 1
         y = 1
        else 
          indices(3) = index - 1
          y = -1
       end if

       indices(4) = index + x*ysize + y 

       do j=1, 4
          if(g_array(indices(j),3).ne.1.0) then
             g_array(indices(j),3) =  1.0
          end if
       end do

      end do

    end if

    if(kernel.eq.'U') then

      do i=1, size(g_array(:,1))

          grid_point = g_array(i,:)

          xrange = (/ g_array(i,1)-bandwidth(1), g_array(i,1)+bandwidth(1) /)
          yrange = (/ g_array(i,2)-bandwidth(2), g_array(i,2)+bandwidth(2) /)

          do j=1, size(data(:,1))

            data_point = data(j,:)

            if((data_point(1).lt.xrange(2)).and.(data_point(1).gt.xrange(1)) &
              .and.(data_point(2).gt.yrange(1)).and.(data_point(2).lt.yrange(2))) then

              g_array(i,3) = g_array(i,3) + 1

            end if

          end do

      end do

      g_array(:,3) = g_array(:,3)/n 


    end if

    if(kernel.eq.'N') then
 
    do i=1, size(g_array(:,1)) 

    	k_sum = 0.0 ! sum of x.y

    	do j = 1, n 

        ux = (data(j,1)-g_array(i,1))/bandwidth(1)

        uy = (data(j,2)-g_array(i,2))/bandwidth(2)

        rux = -0.5 * (ux**2)
        ruy = -0.5 * (uy**2)

        if(rux.lt.-10.0) cycle ! this is to avoid underflow
        if(ruy.lt.-10.0) cycle

        gux = exp(rux)/sqrt(2*PI)
        guy = exp(ruy)/sqrt(2*PI)

        k_sum = k_sum + gux * guy

    	end do

    	g_array(i,3) = k_sum/n

    end do 

    end if

    end subroutine density

    subroutine get_log_fraction(v_coordinates, g_array, path_size,log_fraction)

    ! as used in the COGMAP project

    ! INPUT
    real, dimension(:,:), intent(in) :: v_coordinates,g_array
    integer, intent(in) :: path_size

    !INTERNALS
    real :: path_sum, x_point
    real, dimension(:), allocatable :: counter
    integer, dimension(4) :: indices
    integer :: i, j, x_index, y_index, index, new_index, ysize, x, y
  
    ! OUTPUT
    real,intent(out) :: log_fraction

    path_sum = 0.0

    allocate(counter(size(g_array(:,1))))
    counter(:) = 0.0

    ! for every point
    do i = 1, (path_size-1)

      x_index = MINLOC(abs(g_array(:,1)-v_coordinates(i,1)),dim = 1)
      y_index = MINLOC(abs(g_array(:,2)-v_coordinates(i,2)),dim = 1)

      index = x_index + y_index - 1 

      indices(1) = index

      x_point = g_array(x_index,1)
      ysize = count((g_array(:,1)-x_point).eq.0.0)

      ! x axis
       if(g_array(index,1).lt.v_coordinates(i,1)) then
          new_index = index + ysize
          x = 1
        else 
          new_index = index - ysize
          x = -1
       end if

      indices(2) = new_index

       ! y axis
       if(g_array(index,2).lt.v_coordinates(i,2)) then
         new_index = index +1
         y = 1
        else 
          new_index = index-1
          y = -1
       end if

       indices(3) = new_index

       indices(4) = index + x*ysize + y

      
       do j=1, 4
          if(counter(indices(j)).ne.1.0) then
             path_sum = path_sum + 1.0
             counter(indices(j)) = 1.0 ! counter to avoid duplicate counting
          end if
       end do

    end do


    log_fraction = log(path_sum/sum(g_array(:,3)))

    end subroutine get_log_fraction

    subroutine get_log_fraction2(v_coordinates, g_array, path_size, log_fraction)

    ! INPUT
    real, dimension(:,:), intent(in) :: v_coordinates,g_array
    integer, intent(in) :: path_size

    !INTERNALS
    real :: path_sum, x_point
    real, dimension(:), allocatable :: counter
    integer :: i, j, x_index, y_index, index, new_index, ysize, x, y
  
    ! OUTPUT
    real,intent(out) :: log_fraction

    path_sum = 0.0

    allocate(counter(size(g_array(:,1))))
    counter(:) = 0.0

    ! for every point
    do i = 1, (path_size-1)

      x_index = MINLOC(abs(g_array(:,1)-v_coordinates(i,1)),dim = 1)
      y_index = MINLOC(abs(g_array(:,2)-v_coordinates(i,2)),dim = 1)

      index = x_index + y_index - 1 

      path_sum = path_sum + g_array(index,3)

    end do

    if(path_sum.ne.0.0)then
      log_fraction = 0.0
    else
      log_fraction = log(path_sum/sum(g_array(:,3)))
    end if 


    end subroutine get_log_fraction2

end module biKDE