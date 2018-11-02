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

!> \file dummy_procs.f90
!! \brief subroutines returning always zeroes for scalars, vectors and 2nd order tensors, and some logical functions returning always .TRUE.
!<

!>  Dummy procs.
!! DRUtES in general solves fully coupled CDRE equations. Let us consider simple heat equation
!! \f[ C\partial_t T = \kappa \Delta T. \f] 
!! this equation has just diffusion term, the other terms (such as convection) are zeroes. This is achieved by linking the convection pointer to this dummy_vector function.
!<

 


module dummy_procs

  public :: dummy_scalar, dummy_vector, dummy_tensor, dummy_logical, time_check_ok
  
  contains



    !> this function returns zero, in order to null zero terms from some particular PDE problems
    function dummy_scalar(pde_loc,i,quadpnt, r) result(null)
      use typy
      use pde_objs
      use global_objs
      class(pde_str), intent(in) :: pde_loc
      !>integer number
      integer(kind=ikind), intent(in) :: i
      type(integpnt_str), intent(in), optional :: quadpnt
      !> real number
      real(kind=rkind), dimension(:), intent(in), optional :: r
      real(kind=rkind) :: null

      null  = 0

    end function dummy_scalar

    !> this procedure returns vector of zeroes
    subroutine dummy_vector(pde_loc, i, quadpnt, r, vector_in, vector_out, scalar)
      use typy
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      !> integer num
      integer(kind=ikind), intent(in) :: i
      type(integpnt_str), intent(in), optional :: quadpnt
      !> real number
      real(kind=rkind), dimension(:), intent(in), optional :: r 
      !> input vector
      real(kind=rkind), dimension(:), intent(in), optional  :: vector_in
      !> output vector
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> scalar value to be supplied by zeroes
      real(kind=rkind), intent(out), optional :: scalar

      if (present(vector_out)) then
	vector_out = 0
      end if

      if (present(scalar)) then
	scalar = 0
      end if

    end subroutine dummy_vector
   
     !> this function returns zero matrix (tensor of second order) 
    subroutine dummy_tensor(pde_loc, i, quadpnt, r, tensor, scalar)
      use typy
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      !> integer number
      integer(kind=ikind), intent(in) :: i
      type(integpnt_str), intent(in), optional :: quadpnt
      !> real number
      real(kind=rkind), dimension(:), intent(in), optional :: r
      !> second order tensor to bye filled by zeroes
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor		
      !> scalar value to be supplied by zeroes
      real(kind=rkind), intent(out), optional :: scalar

      if (present(tensor)) then
	tensor = 0
      end if

      if (present(scalar)) then
	scalar = 0
      end if

    end subroutine dummy_tensor

    !> this function always returns .true.
    function dummy_logical() result(true)
      logical :: true

      true = .true.

    end function dummy_logical
    
    function time_check_ok(pde_loc) result(true)
      use typy
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      logical :: true
            
      true = .true.
      
    end function time_check_ok

end module dummy_procs
