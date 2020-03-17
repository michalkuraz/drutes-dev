! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher
! Copyright 2019  Michal Kuraz, Petr Mayer, Johanna Bloecher, Juliana Arbelaez
! Copyright 2020  Michal Kuraz, Petr Mayer, Johanna Bloecher, Juliana Arbelaez
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

!> \file: vapour_pointers.f90
!! \brief: This module contains subroutines that linkes the functions to the PDE
!<



module vapour_pointers
  public :: vapour_processes
  public :: vapour_linker
  
  contains
    subroutine vapour_processes(processes)
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 2
      
    end subroutine vapour_processes
    
    
    subroutine vapour_linker()
      use typy
      use debug_tools
      use RE_pointers
      use pde_objs
      use evap_reader
      use evap_fnc
      use evap_globals
      use heat_pointers
      use heat_fnc
      use re_constitutive
      use heat_reader
      use evap_bc
      use RE_evap_reader
      
       
      call evap_var()
      
      call Re_evap_var()
      
      call RE_std(pde(1))
      
      
      pde(re_order)%pde_fnc(re_order)%dispersion => diffusion_hh
      
      pde(re_order)%flux => liquid_flux
      
      pde(re_order)%pde_fnc(heat_order)%dispersion => diffusion_hT
      
      pde(re_order)%pde_fnc(re_order)%zerord  =>  source_h
  
      
      call heatlinker(pde(2))
      
     
      deallocate(pde(re_order)%mass)
      
      allocate(pde(re_order)%mass(2))
      
      pde(re_order)%mass(1)%val => vangen
      
      pde(re_order)%mass(2)%val => evap4print
      
      deallocate(pde(re_order)%mass_name)
      
      allocate(pde(re_order)%mass_name(2,2))

      pde(re_order)%mass_name(1,1) = "theta"
      
      pde(re_order)%mass_name(1,2) = "theta [-]"
      
      pde(re_order)%mass_name(2,1) = "evap_rate_norm"
      
      pde(re_order)%mass_name(2,2) = "evaporation rate norm [L.T^{-1}]"
      
      call set_evapbc(pde(re_order))

      pde(heat_order)%pde_fnc(heat_order)%elasticity => capacity_T
      
      pde(heat_order)%pde_fnc(heat_order)%dispersion => diffusion_TT
      
      pde(heat_order)%pde_fnc(re_order)%dispersion => diffusion_Th
      
      pde(heat_order)%pde_fnc(heat_order)%convection => convection_T
      
      pde(heat_order)%pde_fnc(heat_order)%zerord  =>  source_T
       
      pde(heat_order)%initcond => heat_icond  
       
      pde(heat_order)%flux => heatmod_flux
      
      pde(heat_order)%print_mass = .false.
     
      call set_heatbc(pde(heat_order))
      
    
    end subroutine vapour_linker
    
    
    
    subroutine set_evapbc(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use evap_bc
      
      class(pde_str), intent(in out) :: pde_loc
      
      integer(kind=ikind) :: i
      
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
        select case(pde_loc%bc(i)%code)
          case(-1:4)
            continue
            !all done in re_pointers
          case(5)
            pde_loc%bc(i)%value_fnc => water_evap
          case default
            print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
            print *, "the incorrect boundary code specified is:", pde_loc%bc(i)%code
            ERROR stop
        end select
      end do
    end subroutine set_evapbc
    
    
    subroutine set_heatbc(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs 
      use heat_fnc
      use evap_bc 
      use re_constitutive
      
      class(pde_str), intent(in out) :: pde_loc
      
      integer(kind=ikind) :: i
      
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
        select case(pde_loc%bc(i)%code)
          case(0:2)
            CONTINUE
            !already linked from heatpointers
          case(3)
            pde_loc%bc(i)%value_fnc => heat_robin
          case default
            print *, "unrecognized bc option"
            print *, "exited from heat_pointers::heatlinker"
            ERROR STOP
        end select
      end do
      
      
    
    end subroutine set_heatbc
    

end module vapour_pointers

