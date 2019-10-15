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
  
  public :: difussion_hh, difussion_hT
  public :: capacity_T, difussion_Th, difussion_TT, convection_T
  public :: theta_vapor, dtheta_vapordt
  public :: hydraulic_lT
  public :: hydraulic_vh, hydraulic_vT

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
      
      
      real(kind=rkind), dimension(3,3) :: Klh, Kvh
      integer(kind=ikind):: D
      
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      call mualem(pde_loc, layer, quadpnt,tensor = Klh(1:D,1:D))
      
      Kvh(1:D,1:D) = hydraulic_vh(pde_loc, layer, quadpnt, x)* !!
      
      tensor(1:D,1:D) = Klh(1:D,1:D) + Kvh(1:D,1:D)
      
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
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      real(kind=rkind), dimension(3,3) :: KlT, KvT
      real(kind=rkind) :: Klt_scalar, Kvt_scalar
      integer(kind=ikind):: D, i
       
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::difussion_hT"
        ERROR STOP
      end if
      
      Klt_scalar = hydraulic_lT(pde_loc, layer, quadpnt) 
      Kvt_scalar = hydraulic_vT(pde_loc, layer, quadpnt)
      
      Klt = 0
      KvT = 0
      
      do i=1, D
        KlT(i,i) = Klt_scalar
        KvT(i,i) = Kvt_scalar 
      end do
      
      tensor(1:D,1:D) = Klh(1:D,1:D) + Kvh(1:D,1:D)
        
    end subroutine difussion_hT
    
    
    !! Convection term for water flow 
    subroutine convection_h(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use re_globals
      use pde_objs

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
      
      real(kind=rkind), dimension(3) :: Kvect
      
        
      if (.not. present(quadpnt) .or. present(vector_out)) then
        print *, "ERROR! output vector undefined, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if 
        
        
      call dmualem_dh(pde_loc, layer, quadpnt, x, vector_in, vector_out = Kvect(1:D))
    
      
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
      
      val = C_liq + C_vap + C_soil
      
    
      
    end subroutine capacity_T
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
      
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      real(kind=rkind) :: T, L, kappa
      real(kind=rkind), dimension(3,3) :: KlT, KvT
      integer(kind=ikind):: D
      
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if
      
      
      kappa = thermal_conduc(pde_loc, layer, quadpnt, x)
      L = latent_heat_wat(pde_loc, layer, quadpnt, x)
      
      KlT(1:D,1:D) = hydraulic_lT(pde_loc, layer, quadpnt, x)
      KvT(1:D,1:D) = hydraulic_vT(pde_loc, layer, quadpnt, x)*
      tensor(1:D,1:D) = kappa* + C_vap*T*KvT(1:D,1:D) + L*KvT(1:D,1:D)
      
    end subroutine difussion_TT
    !! Difussion due to pressure gradient
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
        
        !> relative hydraulic conductivity, (scalar value)
        real(kind=rkind), intent(out), optional :: scalar
          
        real(kind=rkind) :: T, L  
        real(kind=rkind), dimension(3,3) :: Klh, Kvh
        integer(kind=ikind):: D
        
        if (.not. present(quadpnt) .or. present(tensor)) then
          print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
          ERROR STOP
        end if
        
        D = drutes_config%dimen 
        T = pde(Heat_order)%getval(quadpnt)
        L = latent_heat_wat(pde_loc, layer, quadpnt, x)
        
        
        call mualem(pde_loc, layer, quadpnt,tensor = Klh(1:D,1:D))
        Kvh(1:D,1:D) = hydraulic_vh(pde_loc, layer, quadpnt, x)*
        tensor(1:D,1:D) = C_liq*T*Klh(1:D,1:D) + C_vap*T*Kvh(1:D,1:D) +  Kvh(1:D,1:D)*L
          
          
      end subroutine difussion_hT
      !! Convection term for heat flow
      subroutine convection_T(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
        use typy
        use re_globals
        use pde_objs

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
        
        real(kind=rkind), dimension(3) :: Kvect
        integer(kind=ikind) :: D
        real(kind=rkind) :: T
          
          

        
        if (.not. present(quadpnt) .or. present(vector_out)) then
          print *, "ERROR! output vector undefined, exited from evap_fnc::difussion_hh"
          ERROR STOP
        end if 
        
        D = drutes_config%dimen
        T = pde(Heat_order)%getval(quadpnt)
          
        call dmualem_dh(pde_loc, layer, quadpnt, x,  vector_out = Kvect(1:D))
      
        vector_out(1:D) = C_liq*T*Kvect(1:D)
        
        
      end subroutine convection_T
    
    
    
      !!> Thermal Properties of Liquid water
      function hydraulic_lT(pde_loc, layer, quadpnt, x) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use re_constitutive
        use evap_auxfnc
        
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> second order tensor of the unsaturated hydraulic conductivity
        real(kind=rkind), dimension(:,:) :: val
        
        real(kind=rkind) :: T, h
        real(kind=rkind), dimension(3,3) :: Klh
        integer(kind=ikind):: D
        
        if (.not. present(quadpnt) .or. present(tensor)) then
          print *, "ERROR! output tensor undefined, exited from evap_fnc::difussion_hh"
          ERROR STOP
        end if 
        
        
        D = drutes_config%dimen
        h = pde(RE_order)%getval(quadpnt)
        T = pde(Heat_order)%getval(quadpnt)
       
        call mualem(pde_loc, layer, quadpnt,  x, tensor = Klh(1:D,1:D) , scalar)
        val(1:D,1:D)  = Klh(1:D,1:D)*h*GwT*(1/gamma_0)*dsurf_tension_soilwat_dT(pde_loc, layer, quadpnt)  
        
      end function hydraulic_lT
      
      !!> Isothermal Properties of water vapor
      function hydraulic_vh(pde_loc, layer, quadpnt, x) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
        
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: val
        
        real(kind=rkind) :: rh_soil, rho_l,rho_sv,diff
        
        
       
        rh_soil = rh_soil(pde_loc, layer, quadpnt)
        rho_l = rho_l(pde_loc, layer, quadpnt, x) 
        rho_sv = rho_sv(pde_loc, layer, quadpnt, x) 
        diff = vapor_diff_soil(pde_loc, layer, quadpnt, x)
       
       
        val = (diff/rho_l)*rho_sv*((Molw*gravity)/(R_gas*T))*rh
    
      end function hydraulic_vh
      
      !!> Thermal Properties of water vapor
      function hydraulic_vT(pde_loc, layer, quadpnt, x) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
        
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: val
        
        real(kind=rkind) :: rh_soil, rho_l,drho_svdT,diff, enhancement_factor
        
        
        rh_soil= rh_soil(pde_loc, layer, quadpnt)
        rho_l = rho_l(pde_loc, layer, quadpnt) 
        diff = vapor_diff_soil(pde_loc, layer, quadpnt)
        drho_svdT = drho_sv_dT(pde_loc, layer, quadpnt)
        enhancement_factor = enhacement_factor(pde_loc, layer, quadpnt)
        
        val = (diff/rho_l)*enhancement_factor*drho_svdT*rh_soil
        
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
        
        if (.not. present(quadpnt)) then
         print *, "ERROR: you have not specified either integ point "
         print *, "exited from evap_auxfnc::enhacement_factor"
         ERROR stop
        end if


        quadpnt_loc = quadpnt
        
        quadpnt_loc%column = 1
        
        
        val = 
        
      end function dtheta_vapordt
      
      !!> Water vapor
      function theta_vapor(quadpnt, layer, quadpnt, x) result(val)
      
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: val
        
        real(kind=rkind) :: rh_soil, rho_l,rho_sv
      
      
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified either integ point "
          print *, "exited from evap_auxfnc::enhacement_factor"
          ERROR stop
        end if
        
        theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
        rh_soil = rh_soil(layer, quadpnt)
        rho_l = rho_l(pde_loc, layer, quadpnt, x) 
        rho_sv = rho_sv(pde_loc, layer, quadpnt, x) 
        
        
        val = (1 - theta_l)*rho_sv*rh_soil*(1.0_rkind/rho_l)
      
      end function theta_vapor
      
    
    
    
    

end module evap_fnc
