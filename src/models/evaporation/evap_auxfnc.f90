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
  
  public :: rh_soil, rho_sv
  

  contains
  
  !< Relative Humudity soil rh_soil [-]
  !Input: 
  function rh_soil(pde_loc, layer, quadpnt, x) result(val)
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
    
     real(kind=rkind):: h,T
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      !quadpnt_loc=quadpnt
      !quadpnt_loc%preproc=.true.
      h = pde(1)%getval(quadpnt)
      T = pde(2)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from evap_auxfnc::rh_soil"
          ERROR STOP
      end if
        h = x(1)
        T = x(2)
    end if
    
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
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      !quadpnt_loc=quadpnt
      !quadpnt_loc%preproc=.true.
      T = pde(2)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from evap_auxfnc::rh_soil"
          ERROR STOP
      end if
        T = x(2)
    end if
    
    val = 1e-3 *exp(31.3716_rkind - (6014.79_rkind/T) - 7.92495e-3*T**3)/T

  end fuction rho_sv
  
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
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      !quadpnt_loc=quadpnt
      !quadpnt_loc%preproc=.true.
      T = pde(2)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from evap_auxfnc::rh_soil"
          ERROR STOP
      end if
        T = x(2)
    end if
    
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
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      !quadpnt_loc=quadpnt
      !quadpnt_loc%preproc=.true.
      T = pde(2)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from evap_auxfnc::rh_soil"
          ERROR STOP
      end if
        T = x(2)
    end if
    
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
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      !quadpnt_loc=quadpnt
      !quadpnt_loc%preproc=.true.
      T = pde(2)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from evap_auxfnc::rh_soil"
          ERROR STOP
      end if
        T = x(2)
    end if
    
    val = 75.6_rkind - 0.1425_rkind*T - 2.38e-4*T**2

  end fuction surf_tension_soilwat
  
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
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      quadpnt_loc=quadpnt
      quadpnt_loc%preproc=.true.
      theta_l =  pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    else
      if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from evap_auxfnc::rh_soil"
          ERROR STOP
      end if
        T = x(2) !!id what is here!!
    end if
    
    val = b1 + b2*theta_l + b3*theta_l**0.5

  end fuction thermal_conduc
  
  
end module evap_auxfnc
