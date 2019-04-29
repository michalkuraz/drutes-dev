
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

module kinfnc

  contains 
    subroutine kinconvect(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use geom_tools
      use debug_tools
      use kinglobs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
      real(kind=rkind), dimension(:), intent(in), optional :: vector_in
      !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
      !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
      !<
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
      real(kind=rkind), intent(out), optional :: scalar
      
      integer(kind=ikind) :: el
      
      el = quadpnt%element
      
      select case (drutes_config%dimen)
      
        case(1)
          vector_out(1) = 1.49_rkind * sign(1.0_rkind, watershed_el(el)%sx) * sqrt(abs( watershed_el(el)%sx))/manning(layer)
      
        case(2)
      
          vector_out(1) = 1.49_rkind * sign(1.0_rkind, watershed_el(el)%sx) * sqrt(abs( watershed_el(el)%sx))/manning(layer)
          
          vector_out(2) = 1.49_rkind * sign(1.0_rkind, watershed_el(el)%sy) * sqrt(abs( watershed_el(el)%sy))/manning(layer)
          
        case(3)
          print *, "kinematic wave has no sense for three-dimensions"
          print *, "exited from kinfnc::kinconvect"
          ERROR STOP
          
      end select
      
    end subroutine kinconvect
    
    
    
    subroutine kinbor(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code

      
      

      if (present(value)) then
        value = 0.0_rkind
      end if

      
      if (present(code)) then
        code = 1
      end if
  
    end subroutine kinbor
      
    subroutine kinematixinit(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs

      class(pde_str), intent(in out) :: pde_loc

      pde_loc%solution(:) = 0.0_rkind
     
      
    end subroutine kinematixinit
    
    
    function rainfall(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use kinglobs
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      integer(kind=ikind), save :: position = 2
      integer(kind=ikind) :: i
      
      
      if (position /= raindata(1)%series(1)%pos) then
        do i=position, raindata(1)%series(1)%pos
          if (raindata(1)%series(1)%data(i-1) < time .and. raindata(1)%series(1)%data(i) > time) then
            position = i
            EXIT
          else if (raindata(1)%series(1)%data(i) < time .and. i ==  raindata(1)%series(1)%pos) then
            position = raindata(1)%series(1)%pos
          end if
        end do
      end if
    
      val = -raindata(el2pt(quadpnt%element))%series(2)%data(position)
        
    
    end function rainfall
    


end module kinfnc
