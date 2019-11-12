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
  public :: liquid_flux,vapor_flux, heatmod_flux

  contains
  
  
    !!> Coefficents for modified Richards equation
    !!> Capacity water flow equation from RE equation 
    !! Difussion due to pressure gradient
    subroutine difussion_hh(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
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
      
      
      real(kind=rkind), dimension(3,3) :: Klh, Kvh
      integer(kind=ikind):: D, i
      real(kind=rkind):: Kvh_scalar
      
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined or integ point, exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if
      
      if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      call mualem(pde_loc, layer, quadpnt, tensor = Klh(1:D,1:D))
      
      
      Kvh_scalar = hydraulic_vh( pde_loc,layer, quadpnt)
    
      
      Kvh = 0 
      do i=1, D
        Kvh(i,i) = Kvh_scalar 
      end do
      
      tensor(1:D,1:D) = Klh(1:D,1:D) + Kvh(1:D,1:D)
      
    end subroutine difussion_hh
    !! Difussion due to temperature gradient
    subroutine difussion_hT(pde_loc, layer, quadpnt,  x, tensor, scalar)
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
      real(kind=rkind) ::  KvT_scalar
      integer(kind=ikind):: D, i
       
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined or integ point, exited from evap_fnc::difussion_hT"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::difussion_hT"
        ERROR STOP
      end if
      
      KlT(1:D,1:D) = hydraulic_lT(pde_loc, layer, quadpnt) 
      KvT_scalar = hydraulic_vT(pde_loc, layer, quadpnt)
      

      KvT = 0
      
      do i=1, D
        
        KvT(i,i) = KvT_scalar 
      end do
      
      tensor(1:D,1:D) = KlT(1:D, 1:D) + KvT(1:D,1:D)
        
    end subroutine difussion_hT
    !! Convection term for water flow 
    subroutine convection_h(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use re_globals
      use pde_objs
      use re_constitutive

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
      integer(kind=ikind):: D
      
        
      if (.not. present(quadpnt) .or. present(vector_out)) then
        print *, "ERROR! output vector undefined or integ point, exited from evap_fnc::convection_h"
        ERROR STOP
      end if 
      
      if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::convection_h"
        ERROR STOP
      end if
        
      D = drutes_config%dimen
        
      call dmualem_dh(pde_loc, layer, quadpnt, x, vector_in, vector_out = Kvect(1:D))
      vector_out = Kvect(1:D)
      
    end subroutine convection_h
    
    !!> Coefficents for Heat equation
    !!> Capacity heat equation
    function capacity_T(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use re_globals
      use pde_objs
      use evap_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      real(kind=rkind) :: val 
  
      
      if (.not. present(quadpnt)) then
        print *, "ERROR! integ point undefined, exited from evap_fnc::capacity_T"
        ERROR STOP
      end if 
      
      if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::capacity_T"
        ERROR STOP
      end if
      
      val = C_liq + C_vap + C_soil

      
    end function capacity_T
    !! Difussion due to temperature gradient
    subroutine difussion_TT(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs
      use evap_auxfnc

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
      
      real(kind=rkind) :: T, L, kappa, KvT_scalar
      real(kind=rkind), dimension(3,3) :: KlT, KvT, kappa_tensor
      integer(kind=ikind):: D, i, j
      
      if (.not. present(quadpnt) .or. present(tensor)) then
        print *, "ERROR! output tensor undefined or integ point, exited from evap_fnc::difussion_TT"
        ERROR STOP
        
      end if
       if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::difussion_TT"
        ERROR STOP
      end if
      
      kappa = thermal_conduc(pde_loc, layer, quadpnt)
      L = latent_heat_wat(quadpnt)
      
      KlT(1:D,1:D) = hydraulic_lT(pde_loc, layer, quadpnt)
      
      KvT_scalar = hydraulic_vT(pde_loc, layer, quadpnt)
      
      KvT = 0 
      do i=1, D
        KvT(i,i) = KvT_scalar 
      end do
      
      kappa_tensor = 0
      
      do j=1, D
        kappa_tensor(i,i) = kappa
      end do
      
      tensor(1:D,1:D) = kappa_tensor + C_vap*T*(KvT(1:D,1:D)) + L*(KvT(1:D,1:D))
      
    end subroutine difussion_TT
    !! Difussion due to pressure gradient
    subroutine difussion_Th(pde_loc, layer, quadpnt,  x, tensor, scalar)
        use typy
        use re_globals
        use pde_objs
        use re_constitutive
        use evap_auxfnc

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
          
        real(kind=rkind) :: T, L , Kvh_scalar 
        real(kind=rkind), dimension(3,3) :: Klh, Kvh
        integer(kind=ikind):: D, i
        
        
        
        if (.not. present(quadpnt) .or. present(tensor)) then
          print *, "ERROR! output tensor undefined or integ point, exited from evap_fnc::difussion_Th"
          ERROR STOP
        end if
        
        if ( present(x) ) then
         print *, "This option is not implemented"
         print *, "exited from evap_fnc::difussion_Th"
        ERROR STOP
       end if
        
        D = drutes_config%dimen 
        T = pde(Heat_order)%getval(quadpnt)
        L = latent_heat_wat(quadpnt)
        
      
        call mualem(pde_loc, layer, quadpnt,tensor = Klh(1:D,1:D))
        
        Kvh_scalar = hydraulic_vh(pde_loc, layer, quadpnt)
        Kvh = 0 
        do i=1, D
          Kvh(i,i) = Kvh_scalar 
        end do
        
        tensor(1:D,1:D) = C_liq*T*Klh(1:D,1:D) + C_vap*T*Kvh(1:D,1:D) +  Kvh(1:D,1:D)*L
          
          
      end subroutine difussion_Th
      !! Convection term for heat flow
    subroutine convection_T(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
        use typy
        use re_globals
        use pde_objs
        use re_constitutive

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
          print *, "ERROR! output vector undefined or integ point, exited from evap_fnc::convection_T"
          ERROR STOP
         end if 
        
        if ( present(x) ) then
         print *, "This option is not implemented"
         print *, "exited from evap_fnc::convection_T"
         ERROR STOP
        end if
        
        D = drutes_config%dimen
        T = pde(Heat_order)%getval(quadpnt)
          
        call dmualem_dh(pde_loc, layer, quadpnt, x,  vector_out = Kvect(1:D))
      
        vector_out(1:D) = C_liq*T*Kvect(1:D)

    end subroutine convection_T
    

    
    subroutine liquid_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use re_constitutive
      
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      real(kind=rkind), intent(out), optional :: flux_length

      real(kind=rkind), dimension(3,3)  :: KlT
      integer :: D
      integer(kind=ikind)  :: i
      real(kind=rkind), dimension(3)  ::  q_liq
      real(kind=rkind), dimension(3)  :: vct
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable :: gradient
      type(integpnt_str) :: quadpnt_loc
      
      real(kind=rkind), dimension(:), allocatable, save :: gradT
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from evap_fnc::liquid_flux"
        ERROR stop
      end if
    
      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))
      
      if (present(quadpnt)) then
        call pde_loc%getgrad(quadpnt, gradT)
      else
        gradT = grad
      end if
      
      D = drutes_config%dimen
      
      if (present(x)) then
        call darcy_law(pde_loc, layer, x=x, flux = q_liq(1:D))
      end if
      
      if (present(quadpnt)) then
        call darcy_law(pde_loc, layer, quadpnt, flux = q_liq(1:D))
      end if
      
      KlT(1:D,1:D) = hydraulic_lT(pde_loc, layer, quadpnt) 
      
      vct(1:D) = q_liq(1:D) - matmul(-KlT(1:D,1:D), gradT(1:D))
      
      if (present(flux_length)) then
        select case(D)
          case(1)
                flux_length = vct(1)
          case(2)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
          case(3)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
        end select
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    
    end subroutine liquid_flux
    
    subroutine vapor_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evap_globals
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length

      real(kind=rkind), dimension(3,3)  :: Kvh, KvT
      integer  :: D
      integer(kind=ikind) :: i
      real(kind=rkind), dimension(3) :: gradh
      real(kind=rkind), dimension(3) :: vct
      real(kind=rkind) :: h, Kvh_scalar, KvT_scalar
      real(kind=rkind), dimension(:), allocatable :: gradient
      type(integpnt_str) :: quadpnt_loc
      
      
       real(kind=rkind), dimension(:), allocatable, save :: gradT
      
        if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from evap_fnc::liquid_flux"
        ERROR stop
      end if
    
      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))
      
      if (present(quadpnt)) then
        call pde_loc%getgrad(quadpnt, gradT)
      else
        gradT = grad
      end if
      
        if (present(quadpnt)) then
        quadpnt_loc=quadpnt
        quadpnt_loc%preproc=.true.
        h = pde_loc%getval(quadpnt_loc)
        call pde_loc%getgrad(quadpnt, gradient)
      else
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from re_constitutive::darcy_law"
          ERROR STOP
        end if
        h = x(1)
        allocate(gradient(ubound(grad,1)))
        gradient = grad
      end if
      D = drutes_config%dimen
      
      Kvh_scalar = hydraulic_vh(pde_loc, layer, quadpnt)
      KvT_scalar = hydraulic_vT(pde_loc, layer, quadpnt)
      
      Kvh = 0
      KvT = 0
      
       do i=1, D
        Kvh(i,i) = Kvh_scalar 
        KvT(i,i) =  KvT_scalar
      end do
      
      
      vct(1:D) = - matmul(-Kvh(1:D,1:D), gradient(1:D)) - matmul(-KvT(1:D,1:D), gradT(1:D))
      
      
       if (present(flux_length)) then
        select case(D)
          case(1)
                flux_length = vct(1)
          case(2)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
          case(3)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
        end select
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    end subroutine vapor_flux
    
    subroutine heatmod_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evap_globals
      use evap_auxfnc
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length

      real(kind=rkind), dimension(3,3)  :: K
      integer                           :: D
      integer(kind=ikind)               :: i
      real(kind=rkind), dimension(3)  :: vct, q_vap, q_liq
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind)::kappa,T, L
      real(kind=rkind), dimension(:), allocatable, save :: gradT
      
       if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from evap_fnc::liquid_flux"
        ERROR stop
      end if
    
      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))
      
      if (present(quadpnt)) then
        call pde_loc%getgrad(quadpnt, gradT)
      else
        gradT = grad
      end if
      
      D = drutes_config%dimen
      
      if (present(x)) then
        call vapor_flux(pde_loc, layer, x=x, flux = q_vap(1:D))
        call liquid_flux(pde_loc, layer, x=x,  flux=q_liq(1:D))
      end if
      
      if (present(quadpnt)) then
        call vapor_flux(pde_loc, layer, quadpnt, flux = q_vap(1:D))
        call liquid_flux(pde_loc, layer, quadpnt, flux=q_liq(1:D))
      end if
      
      
      kappa = thermal_conduc(pde_loc, layer, quadpnt)
      L = latent_heat_wat(quadpnt)
      T = pde(Heat_order)%getval(quadpnt)
      
      
      vct(1:D) =  gradT(1:D)*kappa  + C_liq*T*q_liq(1:D) + C_vap*T*q_vap(1:D) + L*q_vap(1:D)
      
      
       if (present(flux_length)) then
        select case(D)
          case(1)
                flux_length = vct(1)
          case(2)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
          case(3)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
        end select
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    end subroutine heatmod_flux
    
    
    
    
      !!> Thermal Properties of Liquid water
    function hydraulic_lT(pde_loc, layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use re_constitutive
        use evap_auxfnc
        
        
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> second order tensor of the unsaturated hydraulic conductivity
        real(kind=rkind), dimension(3, 3) :: val
        
        real(kind=rkind) :: T, h
        real(kind=rkind), dimension(3,3) :: Klh
        integer(kind=ikind):: D
        
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified  integ point "
          print *, "exited from evap_fnc::hydraulic_lT"
          ERROR stop
        end if
        
        
        D = drutes_config%dimen
        h = pde(RE_order)%getval(quadpnt)
        T = pde(Heat_order)%getval(quadpnt)
       
        call mualem(pde_loc, layer, quadpnt,  tensor = Klh(1:D,1:D))
        val(1:D,1:D)  = Klh(1:D,1:D)*h*GwT*(1/gamma_0)*dsurf_tension_soilwat_dT( quadpnt)  
        
    end function hydraulic_lT
      
      !!> Isothermal Properties of water vapor
    function hydraulic_vh(pde_loc, layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
        
       
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
     
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: val
        
        real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val,diff, T
        
        
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified  integ point "
          print *, "exited from evap_fnc::hydraulic_vh"
          ERROR stop
        end if
       
        rh_soil_val = rh_soil(layer, quadpnt)
        rho_l_val= rho_l( quadpnt) 
        rho_sv_val = rho_sv( quadpnt) 
        diff = vapor_diff_soil(pde_loc, layer, quadpnt)
        T = pde(Heat_order)%getval(quadpnt)
       
       
        val = (diff/rho_l_val)*rho_sv_val*((Molw*gravity)/(R_gas*T))*rh_soil_val
    
    end function hydraulic_vh
      
      !!> Thermal Properties of water vapor
    function hydraulic_vT(pde_loc, layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
        
        
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: val
        
        real(kind=rkind) :: rh_soil_val, rho_l_val,drho_svdT_val,diff, enhancement_factor_val
        
        
        rh_soil_val= rh_soil( layer, quadpnt)
        rho_l_val = rho_l(quadpnt) 
        diff = vapor_diff_soil(pde_loc, layer, quadpnt)
        drho_svdT_val = drho_sv_dT(quadpnt)
        enhancement_factor_val = enhancement_factor(pde_loc, layer, quadpnt)
        
        val = (diff/rho_l_val)*enhancement_factor_val*drho_svdT_val*rh_soil_val
        
    end function hydraulic_vT

        
    
      !!> Water vapor time derivative
    function dtheta_vapordt(pde_loc, layer, quadpnt_in, x)  result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_auxfnc
      use globals
      
      class(pde_str), intent(in) :: pde_loc
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt_in
       !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x    
      !> return value
      real(kind=rkind) :: val, theta_vapor_curr, theta_vapor_prev
      
      type(integpnt_str) :: quadpnt
  
      quadpnt = quadpnt_in
    
      quadpnt%column = 2
      theta_vapor_prev = theta_vapor(pde_loc,layer, quadpnt) 
    
      quadpnt%column = 1
      theta_vapor_curr= theta_vapor(pde_loc,layer, quadpnt) 
    
    !to be modified
!       val = sinkterm(pde(re_order, layer, quadpnt_in))
      val = 0
      val = val + (theta_vapor_curr - theta_vapor_prev)/ time_step !! check time
        
    end function dtheta_vapordt
    
    
    
    

      
      !!> Water vapor content
    function theta_vapor(pde_loc,layer, quadpnt) result(val)
      use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
      
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: val
        
        real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val, theta_l
      
      
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified either integ point "
          print *, "exited from evap_auxfnc::theta_vapor"
          ERROR stop
        end if
        
        theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
        rh_soil_val = rh_soil(layer, quadpnt)
        rho_l_val = rho_l(quadpnt) 
        rho_sv_val = rho_sv(quadpnt) 
        
        
        val = (1 - theta_l)*rho_sv_val*rh_soil_val*(1.0_rkind/rho_l_val)
      
    end function theta_vapor
      
      
      
      
    
  
end module evap_fnc
