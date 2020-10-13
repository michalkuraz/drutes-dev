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


module evapreader
  public :: evapread
  
  contains
    
    subroutine evapread()
      use typy
      use globals
      use global_objs
      use readtools
      use evapglob
      
      integer(kind=ikind) :: layers, i
      real(kind=rkind), dimension(3) :: tmpdata
      integer :: fileid, ierr
      
      layers = maxval(elements%material(:))
      
      allocate(soil_heat_coef(layers))
      
      open(newunit=fileid, file="drutes.conf/evaporation/evap.conf", status="old", action="read", iostat=ierr)
      
      if (ierr /=0) then
        print *, "unable to open file drutes.conf/evaporation/evap.conf"
        print *, "exiting..."
        ERROR STOP
      end if
      
      do i=1, layers
        call fileread(tmpdata, fileid, checklen=.true.) 
        soil_heat_coef(i)%b1 = tmpdata(1)
        soil_heat_coef(i)%b2 = tmpdata(2)
        soil_heat_coef(i)%b3 = tmpdata(3)
      end do
        
      close(fileid)
    
    end subroutine evapread

end module evapreader
