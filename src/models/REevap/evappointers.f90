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

      
      


      call heat(pde(:))
      
      call evapread()

      pde(re_ord)%pde_fnc(re_ord)%dispersion => REdiffhh
      
      pde(re_ord)%pde_fnc(heat_ord)%dispersion => REdiffhT
      
      pde(re_ord)%pde_fnc(re_ord)%elasticity  =>  REcapacityhh
      
      pde(re_ord)%pde_fnc(heat_ord)%elasticity  =>  REcapacityhT

      pde(re_ord)%flux  =>  totalflux
      
      pde(heat_ord)%pde_fnc(heat_ord)%elasticity => heatcap_TT

      pde(heat_ord)%pde_fnc(heat_ord)%dispersion => heat_cond
          
      pde(heat_ord)%pde_fnc(re_ord)%dispersion => heatdiffTh
      
      pde(heat_ord)%pde_fnc(heat_ord)%convection => convection4heat
      
      pde(heat_ord)%pde_fnc(heat_ord)%zerord  => heatsrc_w_roots 
      
      
!      pde(re_ord)%bc(102)%value_fnc => evap4bc
      
!      pde(heat_ord)%bc(102)%value_fnc => soil_heat_flux
      
      
          
    end subroutine REevap_linker
    
    
    
    subroutine REevap_linker2()
    
    end subroutine REevap_linker2
  

  


end module evappointers
