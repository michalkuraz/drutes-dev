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

module evap_auxfnc
  use typy
  use global_objs
  use re_globals
  use evap_globals
  
  public :: rh_soil, rho_sv, drho_sv_dT,rho_l
  public :: latent_heat_wat, surf_tension_soilwat
  public :: dsurf_tension_soilwat_dT, thermal_conduc
  public :: vapor_diff_air,vapor_diff_soil, tortuosity
  public :: enhacement_factor
  

  contains
  
  !< Relative Humudity soil rh_soil [-]
  !Input: 
  function rh_soil(layer, quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind):: h,T
    
   
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    h = pde(RE_order)%getval(quadpnt)
    T = pde(Heat_order)%getval(quadpnt)

    
    val = exp ((h*MolW*grav)/R_gas*T)
    
  end fuction rh_soil
  
  !Saturated water vapor density rho_sv [kg/m^3]
  !Input: Temperature in 
  function rho_sv(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind):: T
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::rho_sv"
      ERROR stop
    end if
    
    
   
    T = pde(Heat_order)%getval(quadpnt)
   
    
    val = 1e-3 *exp(31.3716_rkind - (6014.79_rkind/T) - 7.92495e-3*T**3)/T

  end fuction rho_sv
  
  function drho_sv_dT(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind):: T
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::drho_sv_dT"
      ERROR stop
    end if
    
    
   
    T = pde(Heat_order)%getval(quadpnt)
   
    
    val = 1e-3 *exp(31.3716_rkind - (6014.79_rkind/T) - 7.92495e-3*T**3)/T

  end fuction drho_sv_dT
  
  !Liquid water density  rho_l [kg/m^3]
  !Input: Temperature in 
  function rho_l(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
     real(kind=rkind)::T
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::rho_l"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    val = 1000.0_rkind - 7.37e-3*(T - 4.0_rkind)**2 + 3.79e-5*(T -4.0_rkind)**3

  end fuction rho_l
  
  !Latent heat of evaporation of liquid water  [kg/m^3]
  !Input: Temperature in ºC
  function latent_heat_wat(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind)::T
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::latent_heat_wat"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    
    val = 2.501e-6 - 2369*T

  end fuction latent_heat_wat
  
  
  !Surface Tension Soil-Water  [g/s^2]
  !Input: Temperature in ºC
  function surf_tension_soilwat(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
     real(kind=rkind)::T
    
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::surf_tension_soilwat"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    val = 75.6_rkind - 0.1425_rkind*T - 2.38e-4*T**2

  end fuction surf_tension_soilwat
  
  
  function dsurf_tension_soilwat_dT(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
     real(kind=rkind)::T
    
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::dsurf_tension_soilwat_dT"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    val = - 0.1425_rkind - 4.76e-4*T

  end fuction dsurf_tension_soilwat_dT
  
  !thermal conductivity  []
  !Input: Liquid Water content
  function thermal_conduc(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind)::theta_l
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::thermal_conduc"
      ERROR stop
    end if
  
    theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    
    val = b1 + b2*theta_l + b3*theta_l**0.5

  end fuction thermal_conduc
  
  function vapor_diff_soil(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind)::theta_l, theta_sat, theta_air
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::vapor_diff_soil"
      ERROR stop
    end if
  
    theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    theta_air = 1 - theta_l
    theta_sat = vgset(layer)%ths
    
    val = tortuosity(theta_l)*theta_air*vapor_diff_air(pde_loc, layer, quadpnt)

  end fuction vapor_diff_soil
  
  
  function tortuosity(theta_l) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    
    real(kind=rkind),intent (in):: theta_l
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind):: theta_sat, theta_air
    
    
    
    theta_air = 1 - theta_l
    theta_sat = vgset(layer)%ths
    
    val = ((theta_air)**(7/3))/ (theta_sat**2)

  end fuction tortuosity
  
  function vapor_diff_air(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind)::T
    
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::vapor_diff_air"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    val =  2.12e-5 * (T/Tref)**2

  end fuction vapor_diff_air
  
  function enhacement_factor(pde_loc, layer, quadpnt, x) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    
    real(kind=rkind)::T
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::enhacement_factor"
      ERROR stop
    end if
    
    theta_sat = vgset(layer)%ths
    theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    
    val =  9.5_rkind + 3.0_rkind*(theta_l/theta_sat) -8.5_rkind *exp(-((1.0_rkind + &
    (2.6_rkind/f_c**0.5))*(theta_l/theta_sat))**4)
  end fuction enhacement_factor
  
  
end module evap_auxfnc
