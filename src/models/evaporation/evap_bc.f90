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
  public :: water_neumann
  public :: evaporation
  public :: sensible_heat

  contains

  !> implementation for Robin boundary condition
  !! solution is a scalar function \f[ p \f]
  !! \f val = acoef  \pdv{p}{\vec{n}}  + bcoef p ||\vec{q} ||_2 \f]
  !<
  subroutine heat_robin(pde_loc, el_id, node_order, val, acoef, bcoef, code, valarray) 
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
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: val
      real(kind = rkind), intent(out), optional :: acoef, bcoef
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: valarray
    
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind), dimension(3) :: xyz
      real(kind =rkind):: T, L , kappa
      real(kind = rkind), dimension(3):: normal_vct
      
      real(kind=rkind) :: tmax, tmin,tmean,solar
      real(kind=rkind) ::  e_act, evap, rhmean
      real(kind=rkind) :: radiation, tmaxk,tmink
      real(kind = rkind):: Hs, rad, heat_soil_flux
      logical, save :: run1st=.true.
      character(len=8), save :: evap_units 
      integer(kind =ikind) :: D, i,  edge_id, datapos, dataprev
      integer(kind =ikind) :: num_day,hour,day, month
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
      if (present(acoef) .and. present(bcoef)) then
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
      
          day = pde_loc%bc(edge_id)%series(datapos,2)
          month = pde_loc%bc(edge_id)%series(datapos,3)
          tmax = pde_loc%bc(edge_id)%series(datapos,5)
          tmin = pde_loc%bc(edge_id)%series(datapos,6)
          rhmean = pde_loc%bc(edge_id)%series(datapos,7)
          solar = pde_loc%bc(edge_id)%series(datapos,10)
          
          tmean = (tmax+tmin)/2.0_rkind
          tmink = tmin + Tref
          tmaxk = tmax + Tref
          e_act = ((e_o(tmax) + e_o(tmin))/2.0_rkind)*(rhmean/100.0_rkind)
      
          kappa = thermal_conduc(pde_loc, layer, quadpnt)
          L = latent_heat_wat(quadpnt)
          Hs= sensible_heat(quadpnt, tmean)
          evap = evaporation(layer, quadpnt, rhmean)
          num_day = num_day_fcn (day, month,evap_units)
          rad = radiation_fcn(num_day,latitude,elevation,albedo,e_act,solar,tmink,tmaxk)
          
          call vapor_flux(pde_loc, layer , quadpnt=quadpnt, flux=q_vap(1:D))
          
          call liquid_flux(pde_loc, layer, quadpnt, flux=q_liq(1:D))
          
          heat_soil_flux = rad - Hs - L*evap
      
          normal_vct (1:D) = 0
          val =  - heat_soil_flux - L*norm2(q_vap(1:D))
          acoef = -kappa
          bcoef = C_liq*norm2(q_liq(1:D)) + C_vap*norm2(q_vap(1:D))
              
  
        else
          print *, "evaporation boundary must be time dependent, check record for the boundary", edge_id
          ERROR STOP
        end if
      end if
      if (present(code)) then
        code = 5 ! should be 5? 1: Direchlet, 2: Neumann, 3: no boundary, 4: Dirichlet bc (node is not excluded from x vector),  5: Robin  
      end if

  end subroutine  heat_robin
  
  
  subroutine water_neumann(pde_loc, el_id, node_order, value, code, valarray) 
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
      integer(kind=ikind) :: edge_id, i, datapos, dataprev, D
    
      real(kind=rkind) ::  evap, rhmean, rain, theta
      
    
      

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      quadpnt%column = 2
      layer = elements%material(el_id)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      
         if (present(value)) then
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
      
      
          rhmean = pde_loc%bc(edge_id)%series(datapos,7)
          rain = pde_loc%bc(edge_id)%series(datapos,11)
          theta =  pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
          
          evap = evaporation(layer, quadpnt, rhmean)
              
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


      if (present(code)) then
        code = 2
      end if

  end subroutine water_neumann
    
  function evaporation(layer, quadpnt, rh_air) result(val)
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
      
      
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: rh_air
      real(kind=rkind) :: val
        
      real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val
      
      
      
      rh_soil_val = rh_soil(layer, quadpnt)
      rho_l_val = rho_l(quadpnt) 
      rho_sv_val = rho_sv(quadpnt) 
      
      val = (rh_soil_val*rho_sv_val  - rh_air* rho_sv_val )/(resistance*rho_l_val)
  
  end function evaporation
  
   function sensible_heat(quadpnt, temp_air) result(val)
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
      
    
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt 
        real(kind=rkind) :: temp_air
      real(kind=rkind) :: val
        
      real(kind=rkind) ::T
      
      T = pde(Heat_order)%getval(quadpnt)
      
      val = C_air *((T -temp_air)/resistance)
  
  end function sensible_heat


end module evap_bc
