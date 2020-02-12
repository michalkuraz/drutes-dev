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


!> \file evap_bc.f90
!! \brief This module contains subroutines that read input information from config files and additional input files
!<

module evap_bc
  use pde_objs
  use typy
  use evap_globals
  use debug_tools
  use re_globals

  public :: heat_robin
  public :: water_evap
  public :: evap4print
  

  contains

  !> implementation for Robin boundary condition
  !! solution is a scalar function \f[ p \f]
  !! \f ccoef = acoef  \pdv{p}{\vec{n}}  + bcoef p ||\vec{q} ||_2 \f]
  !!  valarray(1) = acoef
  !!  valarray(2) = bcoef 
  !!  valarray(3) = ccoef
  !< 
  subroutine heat_robin(pde_loc, el_id, node_order, val, code, valarray) 
    use typy
    use globals
    use global_objs
    use pde_objs
    use re_globals
    use core_tools
    use geom_tools
    use debug_tools
    use evap_fnc
    use evap_auxfnc
    use evap_globals
    use re_evap_bc
      
    class(pde_str), intent(in) :: pde_loc
    !>Node id and order
    integer(kind=ikind), intent(in)  :: el_id, node_order
    !>return value
    real(kind=rkind), intent(out), optional   :: val

    !> return type of boundary condition
    integer(kind=ikind), intent(out), optional :: code
    !> return value for Robin boundary
    real(kind=rkind), dimension(:), intent(out), optional :: valarray
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    
    
    
    !!!!-----local variables------------
    !Robin boundary coeficients
    real(kind = rkind) :: acoef, bcoef, ccoef
    type(integpnt_str) :: quadpnt
    !Layer ID
    integer(kind=ikind) :: layer
    !> coordinates
    real(kind=rkind), dimension(3) :: xyz
    !> T:temperature
    !> rho_l: liquid water density
    !> L: latent heat of vaporization
    !> kappa: thermal conductivity
    real(kind =rkind):: T, L , kappa, rho_liq, rho_vapor
    !> temperature maximum, mininum and mean
    !> solar: solar radiation 
    real(kind=rkind) :: tmax, tmin,tmean,solar
    !> evap: evaporation rate
    !> rh_air: air relative humidity
    !> e_act: Actual vapor pressure
    real(kind=rkind) ::  e_act, evap, rh_air
    !> temperature maximum, mininum in Kelvin 
    real(kind=rkind) :: tmaxk,tmink
    !> Hs: sensible heat 
    !> rad: solar radiation 
    !> soil heat flux
    real(kind = rkind):: Hs, rad, heat_soil_flux
    logical, save :: run1st=.true.
    !> time units for evaporation
    character(len=8), save :: evap_units 
    !> local variables
    integer(kind =ikind) :: D, i,  edge_id, datapos, dataprev
    !> Number of te day in ten year, hour,day and month
    integer(kind =ikind) :: num_day,hour, day, month,year
    !> liquid water and watwer vapor flux
    real(kind=rkind), dimension(3) :: q_vap, q_liq
      
      
      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      quadpnt%column = 2
      layer = elements%material(el_id)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      edge_id = nodes%edge(elements%data(el_id, node_order))
      
      
      if (run1st) then
        call evap_datadt_bc(evap_units, pde_loc%bc(edge_id)%series)
        run1st = .false.
      end if    
      if (present(valarray)) then
        if (pde_loc%bc(edge_id)%file) then
          do i = pde_loc%bc(edge_id)%series_pos, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time .and. i < ubound(pde_loc%bc(edge_id)%series,1)) then
              datapos = i + 1
              dataprev = i
              EXIT
            else if (pde_loc%bc(edge_id)%series(i,1) > time .and. i == ubound(pde_loc%bc(edge_id)%series,1)) then
              datapos = i
              dataprev = i-1 
              EXIT
            end if
          end do
      
          call get_calendar(hour, day , month, year)

          tmax = pde_loc%bc(edge_id)%series(datapos,3)
          tmin = pde_loc%bc(edge_id)%series(datapos,2)
          rh_air = pde_loc%bc(edge_id)%series(datapos,4)
          solar = pde_loc%bc(edge_id)%series(datapos,7)
          
          tmean = ((tmax+tmin)/2.0_rkind) + Tref
          tmink = tmin + Tref
          tmaxk = tmax + Tref
          e_act = ((e_o(tmax) + e_o(tmin))/2.0_rkind)*(rh_air/100.0_rkind)
      
          kappa = thermal_conduc(pde_loc, layer, quadpnt)
          L = latent_heat_wat(quadpnt)
          rho_liq = rho_l(quadpnt)
          rho_vapor = rho_sv(quadpnt)*rh_soil(layer, quadpnt)
          !temperature should be in K 
          Hs= sensible_heat(quadpnt, tmean)
          evap = evaporation(layer, quadpnt, rh_air)
          num_day = num_day_fcn (day, month,evap_units)
          rad = radiation_fcn(num_day,latitude,elevation,albedo,e_act,solar,tmink,tmaxk)
          
          call vapor_flux(pde_loc, layer , quadpnt=quadpnt, flux=q_vap(1:D))
          
          call liquid_flux(pde_loc, layer, quadpnt, flux=q_liq(1:D))
          
          heat_soil_flux = rad - Hs - L*evap*rho_liq
      
          ccoef =  - heat_soil_flux - L*norm2(q_vap(1:D))*rho_liq
          acoef = -kappa
          bcoef = C_liq*rho_liq*norm2(q_liq(1:D)) + C_vap*norm2(q_vap(1:D))*rho_vapor
          
          valarray(1) = acoef
          
          valarray(2) = bcoef
          
          valarray(3) = ccoef
              
  
        else
          print *, "evaporation boundary must be time dependent, check record for the boundary", edge_id
          ERROR STOP
        end if
      end if
      if (present(code)) then
        code = 5 ! should be 5? 1: Direchlet, 2: Neumann, 3: no boundary, 4: Dirichlet bc (node is not excluded from x vector),  5: Robin  
      end if

  end subroutine  heat_robin
  
  
  subroutine water_evap(pde_loc, el_id, node_order, value, code, valarray) 
    use typy
    use globals
    use global_objs
    use pde_objs
    use re_globals
    use core_tools
    use geom_tools
    use debug_tools
    use evap_fnc
    use evap_auxfnc
    use evap_globals
    
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in)  :: el_id, node_order
    real(kind=rkind), intent(out), optional   :: value
    integer(kind=ikind), intent(out), optional :: code
    real(kind=rkind), dimension(:), intent(out), optional :: valarray

    
    type(integpnt_str) :: quadpnt
    integer(kind=ikind) :: layer
    real(kind=rkind), dimension(3) :: xyz
    integer(kind=ikind) :: edge_id, i, datapos, D
    integer(kind=ikind), save :: datainit=1
    
  
    real(kind=rkind) ::  evap, rhmean, theta
    
  
  
    quadpnt%type_pnt = "ndpt"
    quadpnt%order = elements%data(el_id, node_order)
    quadpnt%column = 2
    layer = elements%material(el_id)
    D = drutes_config%dimen
    call getcoor(quadpnt, xyz(1:D))
    
    edge_id = nodes%edge(elements%data(el_id, node_order))
    
    if (present(value)) then
      if (pde_loc%bc(edge_id)%file) then
        call get_datapos(pde_loc%bc(edge_id), datapos, datainit=datainit)
        print *, datapos ; call wait()
        rhmean = pde_loc%bc(edge_id)%series(datapos,4)
        theta =  pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
        
        evap = evaporation(layer, quadpnt, rhmean)
        value = evap
      else
        print *, "evaporation boundary must be time dependent, check record for the boundary", edge_id
        ERROR STOP
      end if
    end if


    if (present(code)) then
      code = 2
    end if

  end subroutine water_evap
   
  
  
    !> For printing
    function evap4print(pde_loc,layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_globals
      use evap_auxfnc
      use evap_fnc
      use geom_tools
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      
      !-- local i/o variables----
      integer(kind=ikind) :: el_id, node_order
      !> Rhmean: relatuive humidity of air
      !> evap: evaporation arte
      real(kind=rkind) :: rhmean, evap
      integer(kind=ikind) :: edge_id, datapos
      integer(kind=ikind), save :: datainit=1
      
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point "
        print *, "exited from evap_bc::evap4print"
        ERROR stop
      end if
      
      

      
      select case(quadpnt%type_pnt)
       case ("obpt", "gqnd") 
          el_id = quadpnt%element
        case ("ndpt") 
          el_id = nodes%element(quadpnt%order)%data(1)
      end select
        
      if (isboundary(quadpnt)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))
        call get_datapos(pde_loc%bc(edge_id), datapos, datainit=datainit)
        rhmean = pde_loc%bc(edge_id)%series(datapos,4)
        evap = evaporation(layer, quadpnt, rhmean)
      else
        evap = dtheta_vapordt(pde_loc, layer, quadpnt, x)*elements%areas(el_id)
      end if 
      
    end function evap4print
    
  


end module evap_bc
