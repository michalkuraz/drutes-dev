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

module evap_reader
  use typy
  
  
  public :: evap_var
  
  
  contains
    subroutine evap_var()
      use typy
      use globals
      use global_objs
      use core_tools
      use readtools
      use pde_objs
      use debug_tools
      use re_globals
      use evap_globals
       
      integer :: ierr
      pde(re_order)%problem_name(2) = "Richards' equation with vapour flow"
      open(newunit=file_vapor, file="drutes.conf/evaporation/vapor.conf", action="read", status="old", iostat = ierr)
      if (ierr /= 0) then
        print *, "missing vapor.conf file in drutes.conf/evaporation/vapor.conf"
        ERROR STOP
      end if 

        
      call fileread(b1, file_vapor, ranges=(/- huge(0.0_rkind), huge(0.0_rkind)/), & 
                      errmsg="Empirical regression parameters for the Thermal conductivity [Wm^-1 K^-1], see Chun&Horton, 1987")
      call fileread(b2, file_vapor, ranges=(/-huge(0.0_rkind), huge(0.0_rkind)/), & 
                      errmsg="Empirical regression parameters for the Thermal conductivity [Wm^-1 K^-1], see Chun&Horton, 1987")
      call fileread(b2, file_vapor, ranges=(/-huge(0.0_rkind), huge(0.0_rkind)/), & 
                      errmsg="Empirical regression parameters for the Thermal conductivity [Wm^-1 K^-1], see Chun&Horton, 1987")
      call fileread(resistance, file_vapor, ranges=(/0.0_rkind, huge(0.0_rkind)/), & 
                      errmsg="specify the aerodynamic resistance to water vapor flow and heat transfer")
      close(file_evap)	

    end subroutine evap_var
  
  
end module evap_reader
