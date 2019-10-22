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

  subroutine heat _bc(pde_loc, el_id, node_order, value, code) 
       use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use core_tools
      use geom_tools
      use debug_tools
      use eva_fnc
      use evap_auxfnc
      use evap_globals
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: value
      integer(kind=ikind), intent(out), optional :: code

      
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: xyz
      real(kind =rkind):: T, L , kappa
      

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      quadpnt%column = 2
      layer = elements%material(el_id)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      
      T = pde(Heat_order)%getval(quadpnt)
      kappa = thermal_conduc(pde_loc, layer, quadpnt)
      L = latent_heat_wat(quadpnt)
      call vapor_flux(pde_loc, layer, quadpnt, x, grad,  flux = q_vap(1:D), flux_length)
      call liquid_flux(pde_loc, layer, quadpnt, x, grad,  flux= q_liq(1:D), flux_length)
      

      if (present(code)) then
        code = 2
      end if

  end subroutine heat_bc
  
  
  subroutine water _bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use core_tools
      use geom_tools
      use debug_tools
      use eva_fnc
      use evap_auxfnc
      use evap_globals
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: value
      integer(kind=ikind), intent(out), optional :: code

      
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: xyz
    
      

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      quadpnt%column = 2
      layer = elements%material(el_id)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      call vapor_flux(pde_loc, layer, quadpnt, x, grad,  flux = q_vap(1:D), flux_length)
      call liquid_flux(pde_loc, layer, quadpnt, x, grad,  flux= q_liq(1:D), flux_length)
      
      
      

      if (present(code)) then
        code = 2
      end if

  end subroutine water_bc
    
    
    

end module evap_bc
