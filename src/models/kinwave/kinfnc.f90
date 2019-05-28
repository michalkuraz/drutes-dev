
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

  public :: kinconvect, kinbor, kinematixinit, rainfall, kin_elast, getval_kinwave

  contains 
  
    !> specific function for kinematic wave equation, replaces pde_objs::getvalp1 surface runoff should be in [mm]
    function getval_kinwave(pde_loc, quadpnt) result(val)
      use typy
      use pde_objs
      use geom_tools
      use re_globals
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind), dimension(3) :: xyz
      integer(kind=ikind) :: D, layer
      
      val = getvalp1(pde_loc, quadpnt)
           
      if (quadpnt%preproc) then
        val = val*1e3
      end if
	
      
    end function getval_kinwave
  
  
  
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
      real(kind=rkind) :: hsurf, m
      
      el = quadpnt%element
      
      hsurf = max(0.0_rkind, pde_loc%getval(quadpnt))
      
      m = 5.0_rkind/3
  
      select case (drutes_config%dimen)
      
        case(1)
          vector_out(1) = -1.49_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                          sqrt(abs( watershed_el(el)%sx))/manning(layer)*m*hsurf**(m-1)
      
        case(2)
      
          vector_out(1) = -1.49_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                          sqrt(abs( watershed_el(el)%sx))/manning(layer)*m*hsurf**(m-1)
          
          vector_out(2) = -1.49_rkind * sign(1.0_rkind, watershed_el(el)%sy) * & 
                          sqrt(abs( watershed_el(el)%sy))/manning(layer)*m*hsurf**(m-1)
                          
          
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
      integer(kind=ikind), save :: position = 1
      integer(kind=ikind) :: i
      
      
       if (position < ubound(raindata(1)%series,1)) then
        do i=position, ubound(raindata(1)%series,1)-1
           if (raindata(1)%series(i,1) < time .and. raindata(1)%series(i+1,1) > time) then
            position = i
            EXIT
          else if (raindata(1)%series(i+1,1) < time .and. i == ubound(raindata(1)%series,1) ) then
            position = i
          end if
        end do
      end if
    
      val = raindata(el2pt(quadpnt%element))%series(position,2)
      
        
    
    end function rainfall
    
    
    function kin_elast(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use heat_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      end if


      E = 1.0_rkind
      

    end function kin_elast 
    


end module kinfnc
