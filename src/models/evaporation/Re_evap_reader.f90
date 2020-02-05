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

!> \file Re_evap_reader.f90
!! \brief This module contains subroutines that read input information from config files and additional input files
!<


module Re_evap_reader
  public :: Re_evap_var
  contains
    subroutine Re_evap_var()
      use typy
      use globals
      use global_objs
      use core_tools
      use readtools
      use pde_objs
      use debug_tools
      use re_globals
      
       
      integer :: ierr
      open(newunit=file_evap, file="drutes.conf/evaporation/evap.conf", action="read", status="old", iostat = ierr)
      if (ierr /= 0) then
        print *, "missing evap.conf file in drutes.conf/evaporation/evap.conf"
        ERROR STOP
      end if 

      call fileread(elevation, file_evap, ranges=(/0.0_rkind, huge(0.0_rkind)/), & 
                      errmsg="specify elevation above sea level in m. This can be only positive")
      call fileread(latitude, file_evap, ranges=(/- huge(0.0_rkind), huge(0.0_rkind)/), &
                      errmsg="specify latitude in radians. This can be negative or positive depend on hemisphere")
      call fileread(albedo, file_evap, ranges=(/0.0_rkind, 1.0_rkind/), &
                      errmsg="specify albedo or canopy reflection coefficient between 0 and 1")
                      
      call fileread(day_in_month, file_evap, ranges=(/1_ikind, 31_ikind/), &
                      errmsg="Days in a month can't be more than 31...:)")
      call fileread(month_in_year, file_evap, ranges=(/1_ikind, 12_ikind/), &
                      errmsg="Month in a year can't be more than 12...:)")          
      call fileread(init_year, file_evap, ranges=(/1_ikind, huge(0_ikind)/), &
                      errmsg="Year can be only positive:)")         
      close(file_evap)	

    end subroutine Re_evap_var
end module Re_evap_reader
