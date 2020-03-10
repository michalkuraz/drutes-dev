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


!> \file Re_evap_bc.f90
!! \brief This module contains subroutines that read input information from config files and additional input files
!<

module Re_evap_bc

  public :: evap_pm_bc
  public :: evap_datadt_bc
  public :: e_o
  public :: pressure_atm
  public :: wind_fcn
  public :: radiation_fcn
  public :: soilheat_fcn
  public :: num_day_fcn
  
  contains

  !> Defines dt of provide data for eveporation calculations
    subroutine evap_datadt_bc(evap_units, series)
      use typy
      use globals
      use core_tools
      use debug_tools


      real(kind=rkind), dimension(:,:), intent(in) :: series
      character(len=*), intent(out) :: evap_units
      real(kind=rkind) :: datascale
      
      integer :: i
      
      
      if (ubound(series,1)> 1) then
        datascale = series(2,1) - series(1,1)
      else
        datascale = series(1,1)
      end if
     

      
      select case(nint(datascale))
        case(0:3600)
          evap_units  = "hourly"
        case(3601:86400)
          evap_units  = "daily"
        case(86401:2678400)
          evap_units  = "monthly"
        case(2678401:31708800)
          evap_units  = "yearly"
        case default
          print *, "incorrect input data"
          print *, "exited from  Re_evap_bc::evap_datadt_bc"
          ERROR STOP
      end select

    end subroutine evap_datadt_bc
  
  
  
  !> Defines Neumann (flux) evaporation boundary condition using the Penman-Monteith Model
    subroutine evap_pm_bc(pde_loc, el_id, node_order, value, code, valarray) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use core_tools
      use geom_tools
      use debug_tools
      use evap_auxfnc
      use Re_evap_reader
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: valarray
      
      integer(kind=ikind) :: edge_id, i, D,num_day,hour,day, month,year,layer
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: xyz
      real(kind=rkind) :: tmax, tmin,tmean,tmean_prev,tmax_prev,tmin_prev,wind,solar,soil
      real(kind=rkind) ::  slope_vap,e_sat,e_act, Patm,gp, light,evap, rhmean,wind_elev
      real(kind=rkind) :: radiation, tmaxk,tmink,wind2,rain, theta
      logical, save :: run1st=.true.
      character(len=8), save :: evap_units 
      integer(kind=ikind), save :: datapos=1, dataprev=1
      integer, save :: fileid
      
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      
      
      if (run1st) then
      	call Re_evap_var()
        call evap_datadt_bc(evap_units, pde(re_order)%bc(edge_id)%series)
        run1st = .false.
        open(newunit=fileid, file="out/evap.rate", status="new", action="write")
      end if
      

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      quadpnt%column = 2
      layer = elements%material(el_id)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      Patm = pressure_atm(elevation) 
      gp = 0.665e-3*Patm !psychrometric constant [kPa/•C]
      
      
      if (present(value)) then

        if (pde(re_order)%bc(edge_id)%file) then
          call get_datapos(pde(re_order)%bc(edge_id), datapos, dataprev)
          
          
          call get_calendar(hour, day , month, year)
          print*, hour, day,month, year 
          tmax = pde(re_order)%bc(edge_id)%series(datapos,3)
          tmin = pde(re_order)%bc(edge_id)%series(datapos,2)
          tmax_prev = pde(re_order)%bc(edge_id)%series(dataprev,3)
          tmin_prev =  pde(re_order)%bc(edge_id)%series(dataprev,2)
          rhmean = pde(re_order)%bc(edge_id)%series(datapos,4)
          wind = pde(re_order)%bc(edge_id)%series(datapos,5)
          light = pde(re_order)%bc(edge_id)%series(datapos,6)
          solar = pde(re_order)%bc(edge_id)%series(datapos,7)
          
         
          
          rain = pde(re_order)%bc(edge_id)%series(datapos,8)
          theta =  pde(re_order)%mass(1)%val(pde(re_order), layer, quadpnt)
          
          tmean = (tmax+tmin)/2.0_rkind
          tmean_prev = (tmax_prev+tmin_prev)/2.0_rkind
          tmink = tmin + Tref
          tmaxk = tmax + Tref
          wind_elev = 10.0_rkind
          e_sat = ((e_o(tmax) + e_o(tmin))/2.0_rkind)
          e_act = ((e_o(tmax) + e_o(tmin))/2.0_rkind)*(rhmean/100.0_rkind)
          
         
          slope_vap = (4098.0_rkind*e_sat)/(tmean + Tref)**2.0_rkind

          
          !> num_day calculation
          num_day = num_day_fcn (day, month,evap_units)
          !> Net Radiation calculation             
          radiation = radiation_fcn(num_day,latitude,elevation,albedo,e_act,solar,tmink,tmaxk)
          !> Wind velocity calculation
          wind2 = wind*wind_fcn(wind_elev)
          !> Soil Flux calculation
          soil = soilheat_fcn(tmean,tmean_prev,radiation,hour,evap_units)
          !Evaporation rate [mm/day]
          evap = (0.408_rkind*(radiation - soil)*gp*(900.0_rkind/(tmean + Tref))*wind2*(e_sat - e_act))&
          / (slope_vap + gp*(1.0_rkind + 0.34_rkind*wind2))
          
          
          
          write(unit=fileid, fmt=*) time, evap, soil, radiation 
           
          !conversion to [m/s]
          evap = evap/(1e3*86400.0)   
          
          
           
          if ((rain - evap) >= 0) then
            value = rain - evap
          else
            value = rain - evap*theta**(2.0_rkind/3.0_rkind)
          end if


        else
          print *, "evaporation boundary must be time dependent, check record for the boundary", edge_id
          ERROR STOP
        end if
      end if
!		end if
      

      if (present(code)) then
        code = 2
      end if

    end subroutine evap_pm_bc
 
   
  !> Actual vapor pressure function
  !> val: Actual vapor pressure[kPa]
    function e_o(x) result (val) 
      use typy
        !x: Temperature in [°C]
      real (kind=rkind), intent(in) :: x
      real (kind=rkind) :: val
      val = 0.6108_rkind*exp(17.27_rkind*x/(x + 237.3_rkind))
    end function e_o
    
  !> Atmospheric  pressure function
  !> val: Atmospheric  pressure [kPa]
    function pressure_atm (x) result (val) 
      use typy
      !x: elevation [m]
      real (kind=rkind), intent(in) :: x
      real (kind=rkind) :: val
      val = 101.3_rkind*((293.0_rkind- 0.0065_rkind*x)/293.0_rkind)**5.26_rkind
    end function pressure_atm
    
  !> Wind velocity function
    function wind_fcn(x) result (val) 
      use typy
      !x: elevation of the wind measured [m]
      real (kind=rkind), intent(in) :: x
      real (kind=rkind) :: val
      val = (4.87_rkind/log(67.82_rkind*x - 5.42_rkind))
    end function wind_fcn
    
  !> Net radiation function
	!> val: Net Radiation [MJ/m2 day]
    function radiation_fcn(x,y,z,a,e,s,t_max,t_min) result (val) 
      use typy
      use core_tools
      real (kind=rkind), intent(in) :: y,z,a,e,s,t_max,t_min
      !x: # of the day in the year 
      !y: latitude [rad]
      !z: elevation [m]
      !a: albedo[-]
      !s: incoming solar radiation [MJ/m2 day]   
      !e: actual vapor pressure
      !t_min: air temperature min in [K]
      !t_max: air temperature max in [K]
      integer (kind=ikind), intent(in) :: x
      real (kind=rkind) :: val
      real(kind=rkind) :: omega, R_nl,R_so,R_ns,R_a,dr,delta
      
      !inverse distance Earth-sun [rad]
      dr = 1.0_rkind + 0.033_rkind*cos((2.0_rkind*pi()*x)/365.0_rkind)
      !solar declination [rad]
      delta = 0.409_rkind*sin(((2.0_rkind*pi()*x)/365.0_rkind) -1.39_rkind)
      !sunset hour angle [rad]
      omega = acos(-tan(y)*tan(delta))
      !Extraterrestrial radiation [MJ/m2 day]        
      R_a = ((24.0_rkind*60.0_rkind)/pi())*dr*0.0820_rkind*(omega*sin(y)*sin(delta)&
              + cos(y)*cos(delta)*sin(omega))
      !Clear-sky solar radiation [MJ/m2 day] 
      R_so = (0.75_rkind + z*2e-5)*R_a
      !net shortwave radiation MJ/m2 day]
      R_ns = (1.0_rkind-a)*s
      !Net longwave radiation [MJ/m2 day] 
      R_nl = 4.903e-9*((t_min**4.0_rkind + t_max**4.0_rkind)/2.0_rkind)*(0.34_rkind - &
               0.14_rkind*sqrt(e))*(1.35_rkind*(s/R_so) - 0.35_rkind)
      !Net radiation [MJ/m2 day]         
      val = R_ns - R_nl
    end function radiation_fcn
    
  !>Soil heat flux function
  !val: 
    function soilheat_fcn(t1,t2,r,h,u) result(val)
      use typy
      use core_tools
      !t1: temp mean current [°C]
      !t2: temp mean max [°C]
      !r: Net radiation [MJ/m2 day] 
      real (kind=rkind), intent(in) :: t1,t2,r
      !h: hour
      !u: evap_units (e.g hourly, daily, monthly,yearly)
      character(len=*), intent(in) :: u
      integer(kind =ikind), intent(in) :: h 
      real (kind=rkind) :: val
      
      select case(cut(u))
        case("hourly")
          if (h .LE.  20 .AND. h .GE. 6 ) then
            val = 0.1_rkind*r
          else 
            val = 0.5_rkind*r
          end if
        case("daily")
          val = 0.0_rkind
        case("monthly")
          val = 0.07_rkind*(t1- t2)
        case("yearly")
          val = 0.14_rkind*(t1- t2)
      end select
    end function soilheat_fcn
    
    !>Day of the year function
    !val: # of day in the year
    function num_day_fcn(x,y,z) result(val)
      use typy
      !z: evap_units (e.g hourly, daily, monthly,yearly)
      character(len=*), intent(in) :: z
      !x: day
      !y:month
      integer(kind =ikind), intent(in) :: x,y
      real (kind=rkind) :: val
      
      select case(z)
        case("hourly")
          val = nint(((275.0_rkind/9.0_rkind)*y-30.0_rkind + x)-2.0_rkind)
            if (y < 3) then
              val = val + 2
            end if
        case("daily")
          val = nint(((275.0_rkind/9.0_rkind)*y-30.0_rkind + x)-2.0_rkind)
        case("monthly")
          val = nint(30.4_rkind*y- 15.0_rkind)
            if (y < 3) then
              val = val + 2
            end if
        case("yearly")
          !>Assuming summer condition
          val = 183
      end select
    
    end function num_day_fcn
    
end module Re_evap_bc
