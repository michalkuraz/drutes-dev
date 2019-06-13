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

module Re_evap_bc

  public :: evap_pm_bc
  public :: evap_datadt_bc
  
  
 
  contains

  !> Defines dt of provide data for eveporation calculations
  subroutine evap_datadt_bc
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use re_globals
    

      integer(kind=ikind) :: evap_units
      real(kind=rkind) :: datadt
      logical, save :: run1st=.true.
      
      
      
        if (runf1st) then
        select case(time_units)
          case("s")
            datadt = (1.0_rkind/86400.0_rkind)*datadt
          case("min")
            datadt = (1.0_rkind/1440.0_rkind)*datadt
          case("hrs")
            datadt = (1.0_rkind/24.0_rkind)*datadt
          case("day")
            continue
          case("month")
            datadt = 30.0_rkind*datadt
          case("year")
            datadt = 365.0_rkind*datadt
          case default
            ERROR STOP
        end select
        
        select case(nint(time_units))
          case(1.0_rkind/24.0_rkind)
            evap_units  = "hourly"
          case(1)
            evap_units  = "daily"
          case(28.0_rkind:31.0_rkind)
            evap_units  = "monthly"
          case(365.0_rkind)
            evap_units  = "yearly"
          case default
            ERROR STOP
        end select
        run1st = .false.
      end if
      
      
      
  end subroutine evap_datadt_bc
  
  
  
  !> Defines Neumann (flux) evaporation boundary condition using the Penman-Monteith Model
  subroutine evap_pm_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use re_globals
      use evap_datadt_bc
      
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: value
      integer(kind=ikind), intent(out), optional :: code

      integer(kind=ikind) :: edge_id, i, j, D,num_day
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: xyz
      real(kind=rkind) :: tmax, tmin, tmean,rhmean, wind,solar,soil, slope_vap,e_soil,e_air, Patm,gp, light
      
      edge_id = nodes%edge(elements%data(el_id, node_order))

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      Patm = 
      gp = 
      
      if (present(value)) then
        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              tmax = pde_loc%bc(edge_id)%series(j,2)
              tmin = pde_loc%bc(edge_id)%series(j,3)
              rhmean = pde_loc%bc(edge_id)%series(j,4)
              wind = pde_loc%bc(edge_id)%series(j,5)
              solar = pde_loc%bc(edge_id)%series(j,6)
              light = pde_loc%bc(edge_id)%series(j,7)
              tmean = (tmax+tmin)/2.0_rkind
              e_soil = 0.6108_rkind*exp(17.27_rkinf*tmean/(tmean + 237.3_rkind))
              slope_vap = (4098.0_rkind*e_soil)/(tmean + 237.3_rkind)**2
              e_air = 0.6108_rkind*exp(17.27_rkind*tmean/(tmean + 237.3_rkind))*(rhmean/100.0_rkind)
              !num_day calculation
              dr = 1.0_rkind + 0.033_rkind*cos(2.0_rkind*3.14159265_rkind*J/365.0_rkind)
              delta = 0.409_rkind*sin((2.0_rkind*3.14159265_rkind*J/365.0_rkind) -1.39_rkind)
              omega = acos(-tan(phi)*tan(delta))
              R_a = (24*60/3.14159265)*dr*0.0820*(omega*sin(phi)*sin(delta) + cos(phi)*cos(delta)* sin(omega))
              R_so = (0.75 + z*2e-5)*R_a
              R_ns = (1-a)*solar
              tmink = tmin + 273.15_rkind
              tmaxk = tmax + 273.15_rkind
              R_nl = 4.903e-9*((tmink**4 + tmaxk**4)/2.0_rkind)*(0.34_rkind - 0.14_rkind*sqrt(e_air))*(1.35_rkind*(solar/R_so) - 0.35_rkind)
              radiation = R_ns - R_ln
              wind2 = wind*(4.87_rkind/log(67.82_rkind*z - 5.42_rkind))
              !Soil flux calculation
              value = 
              EXIT
            end if
          end do
        else
          print *, "evaporation boundary must be time dependent, check record for the boundary", edge_id
          ERROR STOP
        end if
      value = 
      end if
      

      if (present(code)) then
        code = 2
      end if

  end subroutine evap_pm_bc
  


end module Re_evap_bc
