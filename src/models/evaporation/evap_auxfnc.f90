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
  public :: enhancement_factor
  public :: evaporation
  public :: sensible_heat
  public :: get_calendar
  private :: set_february
  public :: get_datapos

  contains
  
  !< Relative Humudity soil rh_soil [-]
  !Input: pressure head: h [m]
  !Temperature: T [ºC]
  !Garvity : grvity  [m.s^-2]
  !Molecular weight of water: MolWat [kg mol^-1]
  !Universal Gas constant: R_gas [J mol^-1 K^-1]
  function rh_soil(layer, quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value: Relative Humudity soil rh_soil [-]
    real(kind=rkind):: val
    !> Pressure head: h, Temperature T in ºC and TemperatureT_abs in Kelvin
    real(kind=rkind):: h,T, T_abs
    
   
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not integ point "
      print *, "exited from evap_auxfnc::rh_soil"
      ERROR stop
    end if
    
    h = pde(RE_order)%getval(quadpnt)
    T = pde(Heat_order)%getval(quadpnt)
    T_abs = T + Tref ! T from ºC to Kelvin
    
    val = exp ((h*MolWat*gravity)/(R_gas*T_abs))
    
  end function rh_soil
  
  !Saturated water vapor density rho_sv [kg/m^3]
  !Input: Temperature in T [ºC]
  function rho_sv(quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
  
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value: Saturated water vapor density rho_sv [kg/m^3]
    real(kind=rkind):: val
    !>  Temperature T in ºC and TemperatureT_abs in Kelvin
    real(kind=rkind):: T, T_abs
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified integ point "
      print *, "exited from evap_auxfnc::rho_sv"
      ERROR stop
    end if
   
    T = pde(Heat_order)%getval(quadpnt)
    T_abs = T + Tref ! T from ºC to Kelvin
    
    val = 1e-3 *(exp(31.3716_rkind - (6014.79_rkind/T_abs) - 7.92495e-3*T_abs))/T_abs

  end function rho_sv
  
  !Derivative saturated water vapor density rho_sv [kg/m^3 K]
  !Input: Temperature in T [ºC]
  function drho_sv_dT( quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value : Derivative saturated water vapor density rho_sv [kg/m^3 K]
    real(kind=rkind):: val
    !>  Temperature T in ºC and TemperatureT_abs in Kelvin
    real(kind=rkind):: T, T_abs
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::drho_sv_dT"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    T_abs = T + Tref ! T from ºC to Kelvin
   
    
    val = exp(- 7.92495e-3*T_abs - (6014.79_rkind/T_abs))*((-3.33818e8*T_abs - 4.2122e10)*T_abs + 2.53357e14)*(1/T_abs**3)

  end function drho_sv_dT
  
  !Liquid water density  rho_l [kg/m^3]
  !Input: Temperature in ºC
  function rho_l( quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value:Liquid water density  rho_l [kg/m^3]
    real(kind=rkind):: val
    !>  Temperature T in ºC 
    real(kind=rkind)::T
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified integ point "
      print *, "exited from evap_auxfnc::rho_l"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    val = 1000.0_rkind - 7.37e-3*(T - 4.0_rkind)**2 + 3.79e-5*(T -4.0_rkind)**3

  end function rho_l
  
  !Specific Latent heat of evaporation of liquid water  [Jkg^-1]
  !Input: Temperature in ºC
  function latent_heat_wat(quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
  
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value: specific Latent heat of evaporation of liquid water  [Jkg^-1]
    real(kind=rkind):: val
    !>  Temperature T in ºC 
    real(kind=rkind)::T
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::latent_heat_wat"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt) !T is in ºC
    val = 2.501e-6 - 2369.2_rkind*T

  end function latent_heat_wat
  
  
  !Surface Tension Soil-Water  [g/s^2]
  !Input: Temperature in ºC
  function surf_tension_soilwat( quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value: Surface Tension Soil-Water  [g/s^2]
    real(kind=rkind):: val
    !>  Temperature T in ºC 
     real(kind=rkind)::T
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::surf_tension_soilwat"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    val = 75.6_rkind - 0.1425_rkind*T - 2.38e-4*T**2

  end function surf_tension_soilwat
  
  !Derivative Surface Tension Soil-Water  [g/s^2 ºC]
  !Input: Temperature in ºC  
  function dsurf_tension_soilwat_dT(quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    

    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value: Derivative Surface Tension Soil-Water  [g/s^2 ºC]
    real(kind=rkind):: val
    !>Temperature in ºC 
     real(kind=rkind)::T
    
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified either integ point "
      print *, "exited from evap_auxfnc::dsurf_tension_soilwat_dT"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    
    val = - 0.1425_rkind - 4.76e-4*T

  end function dsurf_tension_soilwat_dT
  
  !thermal conductivity  [Wm^-1 K^-1]]
  !Input: Liquid Water content [-]
  function thermal_conduc(pde_loc, layer, quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value: Thermal conductivity  [Wm^-1 K^-1]]
    real(kind=rkind):: val
    !>Liquid Water content [-]
    real(kind=rkind)::theta_l
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::thermal_conduc"
      ERROR stop
    end if
  
    theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    val = b1 + b2*theta_l + b3*theta_l**0.5

  end function thermal_conduc
  
  
  !Vapor Difussivity in soil  [m^2/s]
  !Input: Water content,  Saturated water content, Volumetric air content
  function vapor_diff_soil(pde_loc, layer, quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value: Vapor Difussivity in soil  [m^2/s]
    real(kind=rkind):: val
    !> Water content,Volumetric air content,
    real(kind=rkind)::theta_l, theta_air
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::vapor_diff_soil"
      ERROR stop
    end if
  
    theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    theta_air = 1 - theta_l
   
    
    val = tortuosity(theta_l, layer)*theta_air*vapor_diff_air(quadpnt)

  end function vapor_diff_soil
  
  !Tortuosity factor in gaseous phase [-]
  !Input: Volumetric liquid water  content [-]
  !Saturated water content [-]
  function tortuosity(theta_l, layer) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    !> Volumetric liquid water  content
    real(kind=rkind),intent (in):: theta_l
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value: Tortuosity factor in gaseous phase [-]
    real(kind=rkind):: val
    !> Volumetric air content,  Saturated water content 
    real(kind=rkind):: theta_sat, theta_air
    
    theta_air = 1 - theta_l
    theta_sat = vgset(layer)%ths
    
    val = ((theta_air)**(7/3))/ (theta_sat**2)

  end function tortuosity
  
  
  !Vapor Difussivity in air  [m^2/s]
  !Input: Temperature [ºC]
  function vapor_diff_air(quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> return value: Vapor Difussivity in air  [m^2/s]
    real(kind=rkind):: val
    !>Temperature T in ºC and TemperatureT_abs in Kelvin
    real(kind=rkind)::T, T_abs
    
     if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::vapor_diff_air"
      ERROR stop
    end if
  
    T = pde(Heat_order)%getval(quadpnt)
    T_abs = T + Tref ! T from ºC to Kelvin
    
    val =  2.12e-5 * (T_abs/Tref)**2

  end function vapor_diff_air
  
  
  !Enhacement Factor [-]
  !Input: Volumetric liquid water  content [-]
  !Saturated water content [-]
  function enhancement_factor(pde_loc, layer, quadpnt) result(val)
    use typy
    use global_objs
    use pde_objs
    use evap_globals
    
    class(pde_str), intent(in) :: pde_loc
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind):: val
    !> Volumetric liquid water  content, Saturated water content
    real(kind=rkind)::theta_l, theta_sat, tmp, const
    
    
    if (.not. present(quadpnt)) then
      print *, "ERROR: you have not specified  integ point "
      print *, "exited from evap_auxfnc::enhacement_factor"
      ERROR stop
    end if
    
    theta_sat = vgset(layer)%ths
    theta_l = pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
    const = 1 + (2.6_rkind/sqrt(f_c))
    tmp = exp(- (const * (theta_l/theta_sat))**4)
    
    val =  9.5_rkind + 3.0_rkind*(theta_l/theta_sat) -8.5_rkind *tmp
  end function enhancement_factor
  
    !> Evaporation rate [m/s]
  !> Input: Air Relative humiduty [-]
  function evaporation(layer, quadpnt, rh_air) result(val)
    use typy
    use evap_globals
      
    !>material ID  
    integer(kind=ikind), intent(in) :: layer
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt 
    !> Relative humidity of air 
    real(kind=rkind) :: rh_air
    !> Evaporation rate [m/s]
    real(kind=rkind) :: val
    !> Relative humidity soil
    !> liquid water density 
    !> saturated water vapor density
    real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val
      

    
    rh_soil_val = rh_soil(layer, quadpnt)
    rho_l_val = rho_l(quadpnt) 
    rho_sv_val = rho_sv(quadpnt) 
    
    val = (rh_soil_val*rho_sv_val  - rh_air* rho_sv_val )/(resistance*rho_l_val)
     
      
  
  end function evaporation
  
  
    !> Sensible heat[W/m^2]
  !> Input: Air temperature[K]
  function sensible_heat(quadpnt, temp_air) result(val)
    use typy
    use pde_objs
    use evap_globals
    use global_objs
    
  
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt 
    real(kind=rkind) :: temp_air
    real(kind=rkind) :: val
      
    real(kind=rkind) ::T
    
    T = pde(Heat_order)%getval(quadpnt)
    
    val = C_air*rho_air*((T - temp_air)/resistance)
  
  end function sensible_heat
  
  
  subroutine get_calendar(hour,day,month,year)
    use globals
    use typy
    use re_globals
    
    integer(kind=ikind), intent(out) :: hour, day, month, year
    integer(kind=ikind), dimension(12):: days_in_month   
    
    integer(kind=ikind) :: cum_hour, cum_day, days2end, oldmonth, newmonth, daystotal, day2count

    year = init_year
    
    days_in_month = set_february(year)

    cum_hour = int(time/3600)
    
    cum_day = int(time/86400)
    
    hour = cum_hour - cum_day*24

    days2end = days_in_month(month_in_year) - day_in_month
    
    
    if (days2end >= cum_day) then
      month = month_in_year
      day = day_in_month + cum_day
      year = init_year
      RETURN
    else
      oldmonth = month_in_year
      year = init_year
      cum_day = days_in_month(oldmonth) - days2end + cum_day
      daystotal = 0
      day2count = cum_day
      do        
        if (oldmonth < 12) then
          newmonth = oldmonth + 1
        else
          newmonth = 1
          year = year + 1
          days_in_month = set_february(year)
        end if
        
        daystotal = daystotal + days_in_month(oldmonth) + days_in_month(newmonth)

        if (daystotal > cum_day) then
          month = newmonth
          day = day2count - days_in_month(oldmonth)
          EXIT
        else
          oldmonth = newmonth
          day2count = cum_day - days_in_month(oldmonth)
        end if
      end do    
    end if
      
      
  end subroutine get_calendar  
  

  function set_february(year) result(days_in_month)
    use typy
    
    integer(kind=ikind), intent(in) :: year
    integer(kind=ikind), dimension(12) :: days_in_month
    
    if ((modulo(year,4_ikind) == 0 .and. modulo(year,100_ikind) /= 0) .or. &
           (modulo(year,4_ikind) == 0 .and. modulo(year,100_ikind) == 0 .and. modulo(year,400_ikind) == 0)) then
      days_in_month = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    else
      days_in_month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    end if
    
  end function set_february
  
    subroutine get_datapos(bcstr, datapos, dataprev, datainit)
      use typy
      use pde_objs
      use globals
      
      type(boundary_vals), intent(in) :: bcstr
      integer(kind=ikind), intent(out), optional :: dataprev
      integer(kind=ikind), intent(out) :: datapos
      integer(kind=ikind), intent(in out), optional :: datainit
      
      integer(kind=ikind) :: i, start, fin
      
      
      if (present(datainit)) then 
        start = datainit
      else
        start = 1
      end if
      
      fin = ubound(bcstr%series,1)
       do i = start, ubound(bcstr%series,1) - 1
         ! inside the table
         if (time > bcstr%series(i,1) .and. time < bcstr%series(i+1,1)) then
           datapos = i
           if (present(dataprev)) dataprev = i-1
           if (present(datainit)) datainit = datapos
           EXIT
         !above the upper row of the boundary, the table always starts with zero
         else if (i == fin-1 .and. time > bcstr%series(fin,1)) then
           print *, "Insufficient meteo data provided!!!"
           print *, "Actual simulation time is greater than final record in your meteo data file"
           print *, "Meteorological data doesn't overlap the simulation perion, exiting now..."
           print *, "exited from evap_bc::get_datapos"
           ERROR STOP
         else if (i == 1 .and.  bcstr%series(i,1) > time ) then
           print *, "Insufficient meteo data provided!!!"
           print *, "Your meteo data should provide at least one record prior to simulation start time"
           print *, "HINT: add negative value for the first time record"
           print *, "exited from evap_bc::get_datapos"
           ERROR STOP
          end if
       end do
       
     end subroutine get_datapos 
  
end module evap_auxfnc
