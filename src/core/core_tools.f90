! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher


! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.


!> \file core_tools.f90
!! \brief Core tools 

!> \page Core tools
!! mainly for writing data and managing units. 
!<


module core_tools

  public :: write_log
  public :: avg
  public :: mesh_allocater
  public :: find_unit
  public :: cut
  public :: pi
  public :: determinant, plane_coeff, plane_derivative, hyperplane_coeff, hyper_planeder
  
    !> generic procedure, can be overloaded with different vector/matrix types
  interface mean
    module procedure mean_single
    module procedure mean_array
    module procedure mean_matrix
  end interface mean


  contains
  
  
    function mean_array(input) result(output)
      use typy
      real(kind=rkind), dimension(:), intent(in) :: input
      real(kind=rkind) :: output
      
      integer(kind=ikind) :: i
      real(kind=rkind) :: suma
      
      if (ubound(input,1) < 0) then
        print *, "unallocated array in mean"
        print *, "exited from core_tools::mean_array"
        ERROR STOP
      end if
      
      
      output = sum(input)/ubound(input,1)
      
    end function mean_array
    
    function mean_single(input) result(output)
      use typy
      
      real(kind=rkind), intent(in) :: input
      real(kind=rkind) :: output
      
      output = input
      
    end function mean_single
    
    
    function mean_matrix(input) result(output)
      use typy
      real(kind=rkind), dimension(:,:), intent(in) :: input
      real(kind=rkind), dimension(:), allocatable :: output
      
      integer(kind=ikind) :: i
      
      if (.not. allocated(output)) allocate(output(ubound(input,2)))
      
      if (ubound(output,1) /= ubound(input,2)) then
        print *, "output array has different length than number of input array columns"
        print *, "exited from core_tools::mean_matrix"
        ERROR STOP
      end if
      
      
      
      if (ubound(input,1) < 0) then
        print *, "unallocated array in mean"
        print *, "exited from core_tools::mean_matrix"
        ERROR STOP
      end if
      
      do i=1, ubound(input,2)
        output(i) = sum(input(:,i)) / ubound(input,1)
      end do
      
    end function mean_matrix

 !> abreviation for adjustl(TRIM(ch))
  function cut(ch) result(out_ch)
    
    !> input character
    character(len=*), intent(in) :: ch
    !> output character without the leading free spaces and free spaces behind
    character(len=LEN(adjustl(TRIM(ch)))) :: out_ch    

    out_ch=adjustl(TRIM(ch))
    
  end function cut
  
    
!> writes data into out/DRUtES.log in specific format
  subroutine write_log(text, real1, int1, text2, real2, int2, text3, real3, int3, hidden)
    use typy
    use globals
    character(len=*), intent(in) :: text
    real(kind=rkind), intent(in), optional :: real1, real2, real3
    integer(kind=ikind), intent(in), optional :: int1, int2, int3
    character(len=*), intent(in), optional :: text2, text3
    !> character storing the time data
    character(len=10) :: timer
    !> character storing date info
    character(len=8) :: dater

    character(len=1024) :: ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8
    logical, optional :: hidden

    if (present(real1)) then
      write(unit=ch1, fmt=*) real1
    else
      write(unit=ch1, fmt=*) ""
    end if


    if (present(int1)) then
      write(unit=ch2, fmt=*) int1
    else
      write(unit=ch2, fmt=*) ""
    end if

    if (present(text2)) then
      write(unit=ch3, fmt=*) text2
    else
      write(unit=ch3, fmt=*) ""
    end if

    if (present(real2)) then
      write(unit=ch4, fmt=*) real2
    else
      write(unit=ch4, fmt=*) ""
    end if

    if (present(int2)) then
      write(unit=ch5, fmt=*) int2
    else
      write(unit=ch5, fmt=*) ""
    end if

    
    if (present(text3)) then
      write(unit=ch6, fmt=*) text3
    else
      write(unit=ch6, fmt=*) ""
    end if

    if (present(real3)) then
      write(unit=ch7, fmt=*) real3
    else
      write(unit=ch7, fmt=*) ""
    end if

    if (present(int3)) then
      write(unit=ch8, fmt=*) int3
    else
      write(unit=ch8, fmt=*) ""
    end if

    call date_and_time(dater, timer)
    write(unit=logfile, fmt = *) "---------------------------------------------------------------"
    write(unit=logfile, fmt = *)  timer(1:2), "hrs", " : ", timer(3:4), "min"," : ", timer(5:10), "s", "  -  ", dater(1:4), &
       "/", dater(5:6), "/", dater(7:8), "    ",  trim(text), trim(ch1), trim(ch2), trim(ch3), trim(ch4), trim(ch5), &
	    trim(ch6), trim(ch7), trim(ch8)
    if (.not. present(hidden)) then	    
      write(unit=terminal, fmt=*) trim(text), trim(ch1), trim(ch2), trim(ch3), trim(ch4), trim(ch5), &
		      trim(ch6), trim(ch7), trim(ch8)
    else 
      if (.not. hidden) then
	write(unit=terminal, fmt=*) trim(text), trim(ch1), trim(ch2), trim(ch3), trim(ch4), trim(ch5), &
	      trim(ch6), trim(ch7), trim(ch8)
      end if
    end if
      

 
!     if (present(time)) then
!       write(unit=logfile, fmt = *) "simulation time=", time
!     end if
    call flush(logfile)

  end subroutine write_log
  



   !>simple function, returns arithmetic average out of two reals
    function avg(a,b) result(c)
      use typy
      real(kind=rkind), intent(in) :: a
      real(kind=rkind), intent(in) :: b
      real(kind=rkind) :: c

      c = (a+b)/2.0_rkind

    end function avg

    !> allocates structures for nodes and elements
    subroutine mesh_allocater()
      use typy
      use global_objs
      use globals


      allocate(nodes%id(nodes%kolik))
      allocate(nodes%edge(nodes%kolik))
      allocate(nodes%data(nodes%kolik,drutes_config%dimen))
      allocate(nodes%boundary(nodes%kolik))
      allocate(nodes%el2integ(nodes%kolik))
      allocate(nodes%element(nodes%kolik))
      allocate(elements%length(elements%kolik,drutes_config%dimen + 1))
      allocate(elements%nvect_z(elements%kolik,drutes_config%dimen + 1))
      
      if (drutes_config%dimen > 1) allocate(elements%nvect_x(elements%kolik,drutes_config%dimen + 1))
      if (drutes_config%dimen > 2) allocate(elements%nvect_y(elements%kolik,drutes_config%dimen + 1))
      
      allocate(elements%id(elements%kolik))
      allocate(elements%data(elements%kolik, drutes_config%dimen+1))
      allocate(elements%areas(elements%kolik))
      allocate(elements%ders(elements%kolik, drutes_config%dimen+1, drutes_config%dimen))
      allocate(elements%border(elements%kolik))
      
      allocate(elements%material(elements%kolik))
      
      allocate(elements%neighbours(elements%kolik, ubound(elements%data,2)))
      
      allocate(nodes%boundary_order(nodes%kolik))

      

    end subroutine mesh_allocater


    !> with Fortran 2008 this function is depricated. Searches for the next available unit ID for attaching newly opened file. Use open(newunit instead.
    subroutine find_unit(iunit, start_id)
      integer, intent(out) :: iunit
      integer, intent(in), optional :: start_id
      logical :: op

      if (present(start_id)) then
        iunit = start_id*1 + 1
      else
        iunit = 100 ! the lower values are often reserved
      end if

      do
        inquire(unit=iunit,opened=op)
        if (.not.op) then
          exit
        end if 
        iunit = iunit + 1  
        if (iunit >= huge(iunit) - 1) then
          ERROR STOP "SYSTEM PANIC, unable to find any free unit"
          exit
        end if
      end do
  
  end subroutine find_unit


  function pi() result(value)
    use typy
    real(kind=rkind) :: value
    
    value = 4_ikind*atan(1.0_rkind)
  end function pi
  
  
     !> the matrix dimension must be exactly (2,2) or (3,3)
  function determinant(A) result(det)
    use typy


    real(kind=rkind) :: det
    real(kind=rkind), dimension(:,:), intent(in) :: A


    if (ubound(A,1) /= ubound(A,2)) then
      print *, "ERROR: input matrix is not a square matrix, called from fem_tools::determinant()"
      ERROR STOP
    end if


    select case(ubound(A,1))
      case(2)
            det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      case(3)
            det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2)&
                - A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
      case(4)
      det = A(1,1)*(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
                 -A(2,4)*A(3,3)*A(4,2) - A(2,3)*A(3,2)*A(4,4) - A(2,2)*A(3,4)*A(4,3) ) &
         -A(2,1)*(A(1,2)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,2) + A(1,4)*A(3,2)*A(4,3) &
                 -A(1,4)*A(3,3)*A(4,2) - A(1,3)*A(3,2)*A(4,4) - A(1,2)*A(3,4)*A(4,3) ) &
         +A(3,1)*(A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) &
                 -A(1,4)*A(2,3)*A(4,2) - A(1,3)*A(2,2)*A(4,4) - A(1,2)*A(2,4)*A(4,3) ) &  
         -A(4,1)*(A(1,2)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,2) + A(1,4)*A(2,2)*A(3,3) &
                 -A(1,4)*A(2,3)*A(3,2) - A(1,3)*A(2,2)*A(3,4) - A(1,2)*A(2,4)*A(3,3) )           
      
      case default
            print *, "ERROR: incorrect matrix dimension, called from fem_tools::determinant()"
            ERROR STOP
    end select


  end function determinant
  
  
 !> specify the x and y derivates of a plane defined by three points
  subroutine plane_derivative(a,b,c,xder,yder, element)
    use typy
    !> 1st point of the plane
    real(kind=rkind), dimension(:), intent(in) :: a
    !> 2nd point of the plane
    real(kind=rkind), dimension(:), intent(in) :: b
    !> 3rd point of the plane
    real(kind=rkind), dimension(:), intent(in) :: c
    !> resulting x derivate
    real(kind=rkind), intent(out) :: xder
    !> resulting y derivate
    real(kind=rkind), intent(out) :: yder
    
    !> element id
    integer(kind=ikind), intent(in), optional :: element
    !-------local variables--------------
    real(kind=rkind), dimension(3) :: u
    real(kind=rkind), dimension(3) :: v
    real(kind=rkind), dimension(3) :: n
    integer(kind=ikind) :: i
    real(kind=rkind) :: reps


    reps = 100*epsilon(reps)


    !check if the plane is not horizontal
    if (abs(a(3) - b(3)) < reps*(abs(a(3))-abs(b(3)))  .and.  &
        abs(a(3) - c(3)) < reps*(abs(a(3))-abs(c(3)))  .and.  &
        abs(b(3) - c(3)) < reps*(abs(b(3))-abs(c(3)))) then
      xder = 0.0_rkind
      yder = 0.0_rkind
      RETURN
    else
      CONTINUE
    end if

    !creates the plane vectors 
    do i=1,3
      u(i) = a(i) - b(i)
      v(i) = a(i) - c(i)
    end do

    ! the normal plane vector is defined as
    n(1) = u(2)*v(3) - v(2)*u(3)
    n(2) = u(3)*v(1) - v(3)*u(1)
    n(3) = u(1)*v(2) - v(1)*u(2)

    ! finally the derivate is as follows, the horizontality has been already checked
    ! the verticality check
    if (abs(n(3)) < 1e2*reps) then
      print *, "the mesh is wrong, base function can't be vertical"
      print *, abs(n(3))
      print *, element
      ERROR STOP
    end if 
    xder = -n(1)/n(3)
    yder = -n(2)/n(3)

    
   end subroutine plane_derivative 
   
   
   !> specify the x and y derivates of a plane defined by three points
  subroutine plane_derivativev2(a,b,c,xder,yder)
    use typy
    !> 1st point of the plane
    real(kind=rkind), dimension(:), intent(in) :: a
    !> 2nd point of the plane
    real(kind=rkind), dimension(:), intent(in) :: b
    !> 3rd point of the plane
    real(kind=rkind), dimension(:), intent(in) :: c
    !> resulting x derivate
    real(kind=rkind), intent(out) :: xder
    !> resulting y derivate
    real(kind=rkind), intent(out) :: yder
    
    
    real(kind=rkind), dimension(4) :: coeffs
    real(kind=rkind), dimension(3,3) :: pts
    
    pts(1,:) = a
    pts(2,:) = b
    pts(3,:) = c
    
    call plane_coeff(pts,coeffs)
    
    
    xder = -coeffs(1)/coeffs(3)
    
    yder = -coeffs(2)/coeffs(3)
    
  end subroutine plane_derivativev2
   
   
   subroutine plane_coeff(point, coeff)
      use typy
      
      real(kind=rkind), dimension(:,:), intent(in) :: point
      real(kind=rkind), dimension(:), intent(out) :: coeff      
      
      coeff(1) = (point(2,2)-point(1,2)) * (point(3,3)-point(1,3)) - (point(3,2) - point(1,2))*(point(2,3) - point(1,3)) 
      coeff(2) = (point(2,3)-point(1,3)) * (point(3,1)-point(1,1)) - (point(3,3) - point(1,3))*(point(2,1) - point(1,1)) 
      coeff(3) = (point(2,1)-point(1,1)) * (point(3,2)-point(1,2)) - (point(3,1) - point(1,1))*(point(2,2) - point(1,2))
      coeff(4) = -(coeff(1)*point(1,1) + coeff(2)*point(1,2) + coeff(3)*point(1,3))
      
      
    end subroutine plane_coeff
   
   
   subroutine hyper_planeder(A, B, C, D, dgdx, dgdy, dgdz)
    use typy
    
    real(kind=rkind), dimension(:), intent(in) :: A, B, C, D
    real(kind=rkind), intent(out) :: dgdx, dgdy, dgdz
    
    real(kind=rkind), dimension(5) :: cfs
    
    call hyperplane_coeff(A,B,C,D, cfs)
    
    
    if (abs(cfs(4)) < 1e3*epsilon(cfs(4))) then
      print *, "base function can't be vertical, mesh is wrong"
      print *, "exited from geom_tools::hyper_planeder"
      ERROR STOP
    end if
    
    dgdx = -cfs(1)/cfs(4)
    dgdy = -cfs(2)/cfs(4)
    dgdz = -cfs(3)/cfs(4)
  
  
  end subroutine hyper_planeder
  
  
    !> points have coordinates eg. \f \[A_x, A_y, A_z, A_f \] \f
  !!  vectors \f \vec{v}_1, \vec{v}_2, \vec{v}_3 \f
  !! are computed as 
  !! \f \vec{v}_1 = \vec{B} - \vec{A} \f
  !! \f \vec{v}_2 = \vec{C} - \vec{A} \f
  !! \f \vec{v}_3 = \vec{D} - \vec{A} \f
  !! plane coefficients computed from
  !! \f \mbox{det} \left( \begin{matrix} x - A_x & y - A_y & z - A_z \\ v_{1_x} &  v_{1_y} & v_{1_z} \\ v_{2_x} &  v_{2_y} & v_{2_z} \\ v_{3_x} &  v_{3_y} & v_{3_z} \end{matrix} \right) \f
  !<
  subroutine hyperplane_coeff(a,b,c,d,coeffs)
    use typy
    
    !> 1st point of the hyperplane
    real(kind=rkind), dimension(:), intent(in) :: a
    !> 2nd point of the hyperplane
    real(kind=rkind), dimension(:), intent(in) :: b
    !> 3rd point of the hyperplane
    real(kind=rkind), dimension(:), intent(in) :: c
    !> 4rd point of the hyperplane
    real(kind=rkind), dimension(:), intent(in) :: d
    
    !> hyperplane coefficients
    real(kind=rkind), dimension(:), intent(out) :: coeffs
    
    
    
    real(kind=rkind), dimension(3,3) :: YZF, XZF, XYF, XYZ
    real(kind=rkind), dimension(3,4) :: v
    real(kind=rkind) :: detYZF, detXZF, detXYF, detXYZ 
    integer(kind=ikind) :: i
    

    if (ubound(coeffs,1) /= 5) then
      print *, "hyperplane coefficients vector must have dimension 5"
      print *, "exited from geom_tools::hyperplane_coeff"
      ERROR STOP
    end if

    
    if (ubound(A,1) /= 4 .or. ubound(B,1) /= 4 .or. ubound(C,1) /= 4 .or. ubound(D,1) /= 4) then
      print *, "points for hyperplane must have 4 components"
      print *, "exited from geom_tools::hyperplane_coeff"
      ERROR STOP
    end if
    
    v(1,:) = B - A
    v(2,:) = C - A
    v(3,:) = D - A
    
    do i=1,3
      YZF(i,:) = [v(i,2), v(i,3), v(i,4)]
      XZF(i,:) = [v(i,1), v(i,3), v(i,4)] 
      XYF(i,:) = [v(i,1), v(i,2), v(i,4)]
      XYZ(i,:) = [v(i,1), v(i,2), v(i,3)]
    end do
    
    detYZF = determinant(YZF)
    detXZF = determinant(XZF)
    detXYF = determinant(XYF)
    detXYZ = determinant(XYZ)
    
    coeffs(1) = detYZF
    coeffs(2) = -detXZF
    coeffs(3) = detXYF
    coeffs(4) = -detXYZ
    
    coeffs(5) = -A(1)*detYZF + A(2)*detXZF - A(3)*detXYF + A(4)*detXYZ
    
  end subroutine hyperplane_coeff
  
  
 
  
  
  


end module core_tools
