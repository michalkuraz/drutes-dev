! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher, Copyright 2019  Michal Kuraz, Petr Mayer, Johanna Bloecher, Juliana Arbelaez

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

!> \file evap_fnc.f90
!! \brief This module contains subroutines that read input information from config files and additional input files
!<




module evap_fnc
  use pde_objs
  use typy
  use evap_globals
  use debug_tools
  use re_globals
  
  public :: difussion_hh, difussion_hT, convection
  public :: capacity_T, difussion_Th, difussion_TT, convection_T
  public :: 

  contains
    !!> Coefficents for modified Richards equation
    !!> Capacity water flow equation from RE equation 
    !! Difussion due to pressure gradient
    subroutine difussion_hh(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      real(kind=rkind) :: h
      real(kind=rkind), dimension(3,3) :: Klh, 
      integer(kind=ikind):: D
      
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      call mualem(pde_loc, layer, quadpnt,tensor(1:D,1:D) = Klh(1:D,1:D))
      tensor(1:D,1:D) = Klh(1:D,1:D) + 
      
      
      
      
    end subroutine difussion_hh
    !! Difussion due to temperature gradient
    subroutine difussion_hT(pde_loc, layer, quadpnt,  x, tensor, scalar)
        use typy
        use re_globals
        use pde_objs

        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> second order tensor of the unsaturated hydraulic conductivity
        real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
        real(kind=rkind) :: h
        !> relative hydraulic conductivity, (scalar value)
        real(kind=rkind), intent(out), optional :: scalar
    end subroutine difussion_hT
    !! Convection term for water flow 
    subroutine convection_h(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      real(kind=rkind) :: h
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
        
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if 
        !i can use pde_loc
        
      call dmualem_dh(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
    end subroutine convection_h
    
    !!> Coefficents for Heat equation
    !!> Capacity heat equation
    function capacity_T(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use re_globals
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      real(kind=rkind) :: h
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      call mualem(...)
      
      
      
    end subroutine capacity_T
    !! Difussion due to pressure gradient
    subroutine difussion_Th(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      real(kind=rkind) :: h
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
    end subroutine difussion_th
    !! Difussion due to temperature gradient
    subroutine difussion_TT(pde_loc, layer, quadpnt,  x, tensor, scalar)
        use typy
        use re_globals
        use pde_objs

        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> second order tensor of the unsaturated hydraulic conductivity
        real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
        real(kind=rkind) :: h
        !> relative hydraulic conductivity, (scalar value)
        real(kind=rkind), intent(out), optional :: scalar
    end subroutine difussion_hT
    !! Convection term for heat flow
    subroutine convection_T(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
        use typy
        use re_globals
        use pde_objs

        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> second order tensor of the unsaturated hydraulic conductivity
        real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
        real(kind=rkind) :: h
        !> relative hydraulic conductivity, (scalar value)
        real(kind=rkind), intent(out), optional :: scalar
        real(kind=rkind), dimension(3) :: Kvect
        
        integer(kind=ikind) :: D
        
        
        D = drutes_config%dimen
        
        !i need to use pde(order)
        
    call dmualem_dh(pde_loc, layer, quadpnt, x,  vector_out = Kvect(1:D))
    
     vector_out = kr()*C*pde_loc%getval(quadpnt)
    end subroutine convection_T
    
    
    
    !!> Thermal Properties of Liquid water
    function hydraulic_lT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_globals
      use re_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      real(kind=rkind) :: temperature, press_head
      
     press_head = pde(1)%getval(quadpnt)
     temperature = pde(2)%getval(quadpnt)
     
     call mualem(pde_loc, layer, quadpnt,  x, tensor, scalar)
     val = 
     
      
      
      
    end function hydraulic_lT
    
    !!> Isothermal Properties of water vapor
    function hydraulic_vh(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_globals
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      
      real(kind=rkind) :: temperature, press_head
      
      
      press_head = pde(1)%getval(quadpnt)
      temperature = pde(2)%getval(quadpnt)
      
      val = 
     
      
      
      
    end function hydraulic_vh
    
    !!> Thermal Properties of water vapor
    function hydraulic_vT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_globals
      
      
      
      
      
    end function hydraulic_vT

    !!> Water vapor time derivative
    function dtheta_vapordt(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      

      
      
      quadpnt_loc = quadpnt
      
      quadpnt_loc%column = 1
      
      
      
      
    end function dtheta_vapordt
    
    !!> Water vapor
    function theta_vapor(quadpnt, .....) result(val)
      
      presshead = pde(1)%getval(quadpnt)
      
      presshead_prev = pde(1)%getval(quadpnt)
      
      temp = pde(2)%getval(quadpnt)
    
    
    end function theta_vapor
    
    
    
    
    

end module evap_fnc
