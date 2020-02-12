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
!! \brief: this module contains the primary function to construct the PDE equations for liquid water, water vapor and heat flow
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
  public :: source_h, source_T

  contains
    !!> Coefficents for modified Richards equation
    !!> Capacity water flow equation from RE equation 
    !!> Difussion due to pressure gradient
    subroutine difussion_hh(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs
      use evap_globals
      use re_constitutive

      class(pde_str), intent(in) :: pde_loc
       !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      !> Klh: unsaturated non-thermal conductivity for water
      !> Klv: unsaturated thermal conductivity for water     
      real(kind=rkind), dimension(3,3) :: Klh, Kvh
      !> local variables
      integer(kind=ikind):: D, i
      !> Kvh_scalar: scalar value of unsaturated non-thermal conductivity for water
      real(kind=rkind):: Kvh_scalar
      
      if (present(x)) then
        print *, "ERROR! use quadpnt only"
        print *, "exited from evap_fnc::difussion_hh"
        ERROR STOP
      end if

      
      D = drutes_config%dimen !Dimension of the problem
      call mualem(pde_loc, layer, quadpnt, tensor = Klh(1:D,1:D))
      Kvh_scalar = hydraulic_vh( pde_loc,layer, quadpnt)
      !Filling the tensor of Kvh
      Kvh = 0 
      do i=1, D
        Kvh(i,i) = Kvh_scalar 
      end do
      tensor(1:D,1:D) = Klh(1:D,1:D) + Kvh(1:D,1:D)
      
    end subroutine difussion_hh
    !! Difussion due to temperature gradient: the same of normal RE
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
      !> Klt: total unsaturated non-thermal conductivity of liquid water
      !> Kvt: total unsaturated thermal conductivity for water vapor
      real(kind=rkind), dimension(3,3) :: KlT, KvT
      !> KvT_scalar: scalar value of unsaturated thermal conductivity for water
      real(kind=rkind) ::  KvT_scalar
      !> local variab
      integer(kind=ikind):: D, i
       
      if (.not. present(quadpnt) .or. .not. present(tensor)) then
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

    
    !! Source term of for water flow
    function source_h(pde_loc, layer, quadpnt_in, x)  result(val)
      use typy
      use global_objs
      use pde_objs
      use globals
      use evap_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt_in
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x    
      !> return value
      real(kind=rkind) :: val
      
      val = 0
      val = vgset(layer)%sinkterm
      val =  val + dtheta_vapordt(pde_loc, layer, quadpnt_in)
      
    end function source_h
    
    
    
    !!> Coefficents for Heat equation
    !!> Capacity heat equation
    function capacity_T(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use re_globals
      use pde_objs
      use evap_globals
      use evap_auxfnc

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> return value
      real(kind=rkind) :: val 
      !> rho_liq: liquid water density
      !> rho_vapor: water vapor density
      real(kind=rkind) :: rho_liq,rho_vapor
      
      rho_liq = rho_l( quadpnt)
      rho_vapor = rho_sv(quadpnt)*rh_soil(layer, quadpnt)
      
      if (.not. present(quadpnt)) then
        print *, "ERROR! integ point undefined, exited from evap_fnc::capacity_T"
        ERROR STOP
      end if 
      
      if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::capacity_T"
        ERROR STOP
      end if
      
      val = rho_liq*C_liq + rho_vapor*C_vap + rho_soil*C_soil
    end function capacity_T
    
    
    !! Difussion due to temperature gradient
    subroutine difussion_TT(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs
      use evap_auxfnc
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> temperature
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      !> T:temperature
      !> L: Specific latent heat
      !> rho_liq: liquid water density
      !> rho_vapor: water vapor density
      !> Kvh_scalar: scalar value of unsaturated non-thermal conductivity for water
      real(kind=rkind) :: T, L, kappa, KvT_scalar, rho_liq, rho_vapor
      !> Klt: total unsaturated non-thermal conductivity of liquid water
      !> Kvt: total unsaturated thermal conductivity for water vapor
      !> kappa_tensor: tensor of thermal conductivity
      real(kind=rkind), dimension(3,3) :: KlT, KvT, kappa_tensor
      !> local variables
      integer(kind=ikind):: D, i, j
      

      D = drutes_config%dimen
      
       if ( present(x) ) then
        print *, "This option is not implemented"
        print *, "exited from evap_fnc::difussion_TT"
        ERROR STOP
      end if
      
      T = pde(Heat_order)%getval(quadpnt)
      kappa = thermal_conduc(pde_loc, layer, quadpnt)
      L = latent_heat_wat(quadpnt)
      rho_liq = rho_l(quadpnt)
      rho_vapor = rho_sv(quadpnt)*rh_soil(layer, quadpnt)
      
      KlT(1:D,1:D) = hydraulic_lT(pde_loc, layer, quadpnt)
      
      KvT_scalar = hydraulic_vT(pde_loc, layer, quadpnt)
      
      KvT = 0 
      do i=1, D
        KvT(i,i) = KvT_scalar 
      end do
      
      kappa_tensor = 0
      
      do i=1, D
        kappa_tensor(i,i) = kappa
      end do
  
      
      tensor(1:D,1:D) = kappa_tensor(1:D, 1:D) + C_vap*rho_vapor*T*(KvT(1:D,1:D)) +  &
                      C_liq*rho_liq*T*(KlT(1:D,1:D)) + L*rho_liq*(KvT(1:D,1:D))
                    
      
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
        !> T:temperature
        !> L: Specific latent heat
        !> rho_liq: liquid water density
        !> rho_vapor: water vapor density
        !> Kvh_scalar: scalar value of unsaturated non-thermal conductivity for water
        real(kind=rkind) :: T, L , Kvh_scalar ,rho_liq, rho_vapor
        !> Klh: unsaturated non-thermal conductivity for water
        !> Klv: unsaturated thermal conductivity for water
        real(kind=rkind), dimension(3,3) :: Klh, Kvh
        !> local variables
        integer(kind=ikind):: D, i
        
        if (present(x)) then
          print *, "ERROR! use quadpnt only"
          print *, "exited from evap_fnc::difussion_Th"
          ERROR STOP
        end if
        
        
        D = drutes_config%dimen 
        T = pde(Heat_order)%getval(quadpnt)
        L = latent_heat_wat(quadpnt)
        rho_liq = rho_l( quadpnt)
        rho_vapor = rho_sv(quadpnt)*rh_soil(layer, quadpnt)
        
      
        call mualem(pde_loc, layer, quadpnt,tensor = Klh(1:D,1:D))
        
        Kvh_scalar = hydraulic_vh(pde_loc, layer, quadpnt)
        Kvh = 0 
        do i=1, D
          Kvh(i,i) = Kvh_scalar 
        end do
        
        tensor(1:D,1:D) = rho_liq*C_liq*T*Klh(1:D,1:D) + rho_vapor*C_vap*T*Kvh(1:D,1:D) +  Kvh(1:D,1:D)*L*rho_liq
          
          
    end subroutine difussion_Th
    !! Convection term for heat flow
    subroutine convection_T(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
        use typy
        use re_globals
        use pde_objs
        use re_constitutive
        use evap_auxfnc

        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in) :: layer
        type(integpnt_str), intent(in), optional :: quadpnt    
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
        real(kind=rkind), dimension(:), intent(in), optional :: vector_in
        !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. 
        !> it is the last column of the hydraulic conductivity second order tensor times  
        !>relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
        real(kind=rkind), dimension(:), intent(out), optional :: vector_out
        !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
        real(kind=rkind), intent(out), optional :: scalar
        
        !> Result from the derivative of mualem 
        real(kind=rkind), dimension(3) :: Kvect
        !> local variable
        integer(kind=ikind) :: D
        !> T:temperature
        !> rho_liq: liquid water density
        real(kind=rkind) :: T, rho_liq
          
        

        
        if ( present(x) ) then
         print *, "ERROR: use quadpnt only."
         print *, "exited from evap_fnc::convection_T"
         ERROR STOP
        end if
        
        D = drutes_config%dimen
        T = pde(Heat_order)%getval(quadpnt)
        rho_liq = rho_l( quadpnt)
          
        call dmualem_dh(pde_loc, layer, quadpnt, x,  vector_out = Kvect(1:D))
      
        vector_out(1:D) = C_liq*rho_liq*T*Kvect(1:D)

    end subroutine convection_T
    
    
    
    !> Source term heat flow
    function source_T(pde_loc, layer, quadpnt_in, x)  result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_auxfnc
      use globals
      use evap_globals
      use heat_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt_in
       !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x    
      !> return value
      real(kind=rkind) :: val
      !> Laten heat of vaporization, liquid water density
      real(kind =rkind):: latent_heat, rho_l_val
      
      val = 0
      val = heatpar(layer)%source
      latent_heat = latent_heat_wat(quadpnt_in) 
      rho_l_val = rho_l(quadpnt_in) 
      val =  val + dtheta_vapordt(pde_loc, layer, quadpnt_in)*latent_heat*rho_l_val 
      
      
    end function source_T


    !> Liquid water flux
    subroutine liquid_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use re_constitutive
      
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> Klt: total unsaturated non-thermal conductivity of liquid water
      real(kind=rkind), dimension(3,3)  :: KlT
      !> local variable
      integer(kind=ikind)  :: i, D
      !> Liquid water flux
      real(kind=rkind), dimension(3)  ::  q_liq
      !> resul of the modified flux of liquid water
      real(kind=rkind), dimension(3)  :: vct
      !> pressure head
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable :: gradient
      !> Gauss quadrature point structure local
      type(integpnt_str) :: quadpnt_loc
      !> Temperature gradient
      real(kind=rkind), dimension(:), allocatable, save :: gradT
      
      if (present(x)) then
        print *, "ERROR: use quadpnt only"
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
      
      if (present(quadpnt)) then
        call darcy_law(pde_loc, layer, quadpnt, flux = q_liq(1:D))
      end if
      
      KlT(1:D,1:D) = hydraulic_lT(pde_loc, layer, quadpnt) 
      
      vct(1:D) = q_liq(1:D) + matmul(-KlT(1:D,1:D), gradT(1:D))
      
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
    
    !> Water Vapour flux
    subroutine vapor_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evap_globals
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> KvT: unsaturated thermal conductivity for water
      !> Kvh: unsaturated non-thermal  conductivity for water
      real(kind=rkind), dimension(3,3)  :: Kvh, KvT
      !> local variables
      integer(kind=ikind) :: i, D
      !> pressure gradient
      real(kind=rkind), dimension(:), allocatable :: gradient
      !result of the vapor flux vector
      real(kind=rkind), dimension(3) :: vct
      !> h: pressure head
      !> Kvh_scalar: unsaturated non-thermal  conductivity for water
      !> KvT_scalar: unsaturated thermal conductivity for water
      real(kind=rkind) :: h, Kvh_scalar, KvT_scalar
      !> Gauss quadrature point structure local
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
      
      
      vct(1:D) =  matmul(-Kvh(1:D,1:D), gradient(1:D)) + matmul(-KvT(1:D,1:D), gradT(1:D))
      
      
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
    
    !> Modified heat flux
    subroutine heatmod_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evap_globals
      use evap_auxfnc
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> local variables
      integer(kind=ikind):: i,D
      !vct: result of te flux
      !q_vap water vapor liquid flux
      !q_liq liquid water flux
      real(kind=rkind), dimension(3)  :: vct, q_vap, q_liq
      !> Gauss quadrature point structure local
      type(integpnt_str) :: quadpnt_loc
      !> kappa:thermal hydraulic conductivity
      !> T:Temperature
      !> L:specif latent heat of vaporization
      real(kind=rkind)::kappa,T, L, rho_liq, rho_vap
      !> Temperature gradient
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
      
      print *, "ve fluxu", present(x), present(quadpnt)
      
      D = drutes_config%dimen
      
      if (present(x)) then
      ! this option is forbideen in your liquid flux!!!!!!!
        call vapor_flux(pde_loc, layer, x=x, flux = q_vap(1:D))
        call liquid_flux(pde_loc, layer, x=x,  flux=q_liq(1:D))
      end if
      
      if (present(quadpnt)) then
             print *, "ve fluxu 2"
        call vapor_flux(pde(1), layer, quadpnt, flux = q_vap(1:D))
               print *, "ve fluxu 3"
        call liquid_flux(pde(1), layer, quadpnt, flux=q_liq(1:D))
      end if
      
  print *, "ve fluxu 4"
      
      kappa = thermal_conduc(pde(1), layer, quadpnt)
      L = latent_heat_wat(quadpnt)
      T = pde(Heat_order)%getval(quadpnt)
      rho_liq = rho_l(quadpnt)
      rho_vap = rho_sv(quadpnt)*rh_soil(layer, quadpnt)
      
      
      print *, "ve fluxu 4.5"
      vct(1:D) =  gradT(1:D)*kappa  + C_liq*T*q_liq(1:D) *rho_liq + C_vap*T*q_vap(1:D)*rho_vap + L*q_vap(1:D)*rho_liq
      
      
       if (present(flux_length)) then
         flux_length = norm2(vct(1:D))
       end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
      print *, "ve fluxu 5"
      
    end subroutine heatmod_flux
    
    !> Thermal Properties of Liquid water
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
        !> unsaturated thermal conductivity for liquid water
        real(kind=rkind), dimension(3, 3) :: val
        !> T:temperature
        !> h: pressure head
        real(kind=rkind) :: T, h
        !unsaturated non'thermal conductivity
        real(kind=rkind), dimension(3,3) :: Klh
        !local variable: dimension
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
        val(1:D,1:D)  = Klh(1:D,1:D)*h*GwT*(1/gamma_0)*dsurf_tension_soilwat_dT(quadpnt)  
        
    end function hydraulic_lT
      
    !> Isothermal Properties of water vapor
    function hydraulic_vh(pde_loc, layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
        
       
        class(pde_str), intent(in) :: pde_loc
        !> MaterialID
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        !>unsaturated non-thermal conductuvity of water vapor
        real(kind=rkind) :: val
        !> T:temperature
        !> Rh: relatuive humidity of soil
        !> rho_l: liquid water density
        !> rho_sv: saturated vapor density
        !> D: diffusitivity
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
       
       
        val = (diff/rho_l_val)*rho_sv_val*((MolWat*gravity)/(R_gas*T))*rh_soil_val
    
    end function hydraulic_vh
      
    !> Thermal Properties of water vapor
    function hydraulic_vT(pde_loc, layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evap_globals
        use evap_auxfnc
        
        
        class(pde_str), intent(in) :: pde_loc
        !Material ID
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        !> unsaturated thermal hydraulic conductivity for water vapor
        real(kind=rkind) :: val
        !> Rh: relatuive humidity of soil
        !> rho_l: liquid water density
        !> rho_sv: saturated vapor density
        !> D: diffusitivity
        !> enhancement factor
        real(kind=rkind) :: rh_soil_val, rho_l_val,drho_svdT_val,diff, enhancement_factor_val
        
        
        rh_soil_val= rh_soil( layer, quadpnt)
        rho_l_val = rho_l(quadpnt) 
        diff = vapor_diff_soil(pde_loc, layer, quadpnt)
        drho_svdT_val = drho_sv_dT(quadpnt)
        enhancement_factor_val = enhancement_factor(pde_loc, layer, quadpnt)
        
        val = (diff/rho_l_val)*enhancement_factor_val*drho_svdT_val*rh_soil_val
        
    end function hydraulic_vT
    !> Water vapor time derivative
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
      !> val:return value
      !> vapor content current and previos
      real(kind=rkind) :: val, theta_vapor_curr, theta_vapor_prev
      !> Gauss quadrature point structure 
      type(integpnt_str) :: quadpnt
  
      quadpnt = quadpnt_in
    
      quadpnt%column = 2
      theta_vapor_prev = theta_vapor(pde_loc,layer, quadpnt) 
    
      quadpnt%column = 1
      theta_vapor_curr= theta_vapor(pde_loc,layer, quadpnt) 
      
      val =  (theta_vapor_curr - theta_vapor_prev)/ time_step 
        
    end function dtheta_vapordt
    
    !> Water vapor content
    function theta_vapor(pde_loc,layer, quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_globals
      use evap_auxfnc
      
      class(pde_str), intent(in) :: pde_loc
      !>material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt 
      !vapor volumetric content
      real(kind=rkind) :: val
      !> Rh: relatuive humidity of soil
      !> rho_l: liquid water density
      !> rho_sv: saturated vapor density
      !>theta_l: liquid water content
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
