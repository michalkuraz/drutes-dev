! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

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


module evappointers
  public :: REevap_proc, REevap_linker, REevap_linker2
  contains
  
    subroutine REevap_proc(processes) 
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 2
    
    end subroutine REevap_proc
    
    subroutine REevap_linker()
      use typy
      use global_objs
      use pde_objs
      use globals
      use heat_pointers
      use evapglob
      use evap_RE_constitutive
      use evapreader
      use evap_heat_constitutive
      use re_constitutive
      use evapbc4heat

      integer(kind=ikind) :: i


      call heat(pde(:))
      
      call evapread()

      pde(re_ord)%pde_fnc(re_ord)%dispersion => REdiffhh
      
      pde(re_ord)%pde_fnc(heat_ord)%dispersion => REdiffhT
      
      pde(re_ord)%pde_fnc(re_ord)%elasticity  =>  REcapacityhh
      
      pde(re_ord)%pde_fnc(heat_ord)%elasticity  =>  REcapacityhT
      
      pde(re_ord)%pde_fnc(re_ord)%convection => dKvhdh

      pde(re_ord)%flux  =>  totalflux
      
      pde(heat_ord)%pde_fnc(re_ord)%elasticity => heatcap_Th
      
      pde(heat_ord)%pde_fnc(heat_ord)%elasticity => heatcap_TT

      pde(heat_ord)%pde_fnc(heat_ord)%dispersion => heat_cond
          
      pde(heat_ord)%pde_fnc(re_ord)%dispersion => heatdiffTh
      
      pde(heat_ord)%pde_fnc(heat_ord)%convection => convection4heat
      
      pde(heat_ord)%pde_fnc(heat_ord)%zerord  => heatsrc_w_roots 
      
      pde(heat_ord)%flux  =>  heat_flux4evap
      
      deallocate(pde(re_ord)%mass)
      deallocate(pde(re_ord)%mass_name)
      
      allocate(pde(re_ord)%mass(3))
      allocate(pde(re_ord)%mass_name(3,2))
      
      pde(re_ord)%mass(1)%val => vangen
      pde(re_ord)%mass(2)%val => thetav
      pde(re_ord)%mass(3)%val => satdegree
      
      pde(re_ord)%mass_name(1,1) = "theta_l"
      pde(re_ord)%mass_name(1,2) = "theta_l [-]        "
      
      pde(re_ord)%mass_name(2,1) = "theta_v"
      pde(re_ord)%mass_name(2,2) = " theta_v [-]       "
      
      pde(re_ord)%mass_name(3,1) = "satdgree"
      pde(re_ord)%mass_name(3,2) = "sat degree [%]     "
      
      
      allocate(pde(re_ord)%fluxes(3))
      pde(re_ord)%fluxes(1)%name(1) = "liq_flux"
      pde(re_ord)%fluxes(1)%name(2) = "liquid flux [L/T]"
      pde(re_ord)%fluxes(1)%val => darcy4liq
      
      pde(re_ord)%fluxes(2)%name(1) = "vap_flux"
      pde(re_ord)%fluxes(2)%name(2) = "vapour flux [L/T]"
      pde(re_ord)%fluxes(2)%val => darcy4vap
      
      pde(re_ord)%fluxes(3)%name(1) = "total_flux"
      pde(re_ord)%fluxes(3)%name(2) = "total flux [L/T]"
      pde(re_ord)%fluxes(3)%val => totalflux
      

      
      do i=lbound(pde(heat_ord)%bc,1), ubound(pde(heat_ord)%bc,1)
        if (pde(heat_ord)%bc(i)%code == 3) then
          pde(heat_ord)%bc(i)%value_fnc => ebalance_flux
        end if
      end do
      
      do i=lbound(pde(re_ord)%bc,1), ubound(pde(re_ord)%bc,1)
        if (pde(re_ord)%bc(i)%code == 5) then
          pde(re_ord)%bc(i)%value_fnc => evaporation_bcflux
        end if
      end do   
          
    end subroutine REevap_linker
    
    
    
    subroutine REevap_linker2()
    
    end subroutine REevap_linker2
  

  


end module evappointers
