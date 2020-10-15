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
      
      integer(kind=ikind) :: layers, i, counter, low, high
      integer :: evapconf, ierr, ebalancein, albedodat
      real(kind=rkind), dimension(:), allocatable :: tmpdata
      type(smartarray_real) :: datafiller
      character(len=12) :: ch
      logical :: success
      
      layers = maxval(elements%material(:))
      
      allocate(soil_heat_coef(layers))
      
      open(newunit=evapconf, file="drutes.conf/evaporation/evap.conf", status="old", action="read", iostat=ierr)
      
      if (ierr /=0) then
        print *, "unable to open file drutes.conf/evaporation/evap.conf"
        print *, "exiting..."
        ERROR STOP
      end if
      
      allocate(tmpdata(3))
      
      do i=1, layers
        call fileread(tmpdata, evapconf, checklen=.true.) 
        soil_heat_coef(i)%b1 = tmpdata(1)
        soil_heat_coef(i)%b2 = tmpdata(2)
        soil_heat_coef(i)%b3 = tmpdata(3)
      end do
  
      call fileread(albedo_conf%method, evapconf, ranges=(/1_ikind,2_ikind/))
      
      if (albedo_conf%method == 2) then
        deallocate(tmpdata)
        allocate(tmpdata(5))
        call fileread(tmpdata, evapconf, checklen=.true.,errmsg="Have you defined parameters for computing albedo from equation")
        albedo_conf%theta_min = tmpdata(1)
        albedo_conf%theta_max = tmpdata(2)
        albedo_conf%albd_min = tmpdata(3)
        albedo_conf%albd_max = tmpdata(4)
        albedo_conf%A = tmpdata(5)
        
        call fileread(zref, evapconf)
      else
        deallocate(tmpdata)
        allocate(tmpdata(1))
        call fileread(tmpdata, evapconf, checklen=.true.,errmsg="Have you commented out parameters for computing albedo from & 
                      equation. Note that you have defined an option for providing albedo from measured data")
        zref = tmpdata(1)
        
        open(newunit=albedodat, file="drutes.conf/evaporation/albedo.dat", status="old", action="read", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open file drutes.conf/evaporation/albedo.dat"
          print *, "exiting..."
          ERROR STOP
        end if
        
        deallocate(tmpdata)
        allocate(tmpdata(2))
        
        counter = 0
        
        call datafiller%clear()
        
        do 
          call comment(albedodat)
          read(unit=albedodat, fmt = *, iostat=ierr) tmpdata
          
          if (ierr /= 0) then
            EXIT
          else
            counter = counter + 1
          end if
          
          do i=1, 2
            call datafiller%fill(tmpdata(i))
          end do
        end do
        
        allocate(albedo_conf%albdat(counter,2))
        
        
        do i=1, counter
          low = (i-1)*2
          high = i*2
          albedo_conf%albdat(i,1) = datafiller%data(low+1)
          albedo_conf%albdat(i,2) = datafiller%data(low+2)
        end do
      end if
          
        
    
      
      counter = 0
      
      deallocate(tmpdata)
      allocate(tmpdata(5))
      call datafiller%clear()
      
      do 
        call comment(evapconf)
        read(unit=evapconf, fmt = *, iostat=ierr) tmpdata

        if (ierr /= 0) then
          EXIT
        else
          counter = counter + 1
        end if
   
        do i=1, 5
          call datafiller%fill(tmpdata(i))
        end do
      end do 
      
      if (counter == 0) then
        call file_error(evapconf, "incorrect inputs for energy balance equation")
      end if
      
      allocate(pars4ebalance(counter))
      
      if (datafiller%pos /= counter*5) then
        print *, datafiller%pos, counter*5, datafiller%data
        call file_error(evapconf, "You have defined incorrect number of columns for energy balance parameters...")
      end if
      
      do i=1, counter
        low = (i-1)*5
        high = i*5
        pars4ebalance(i)%day = datafiller%data(low+1)
        pars4ebalance(i)%latitude = datafiller%data(low+2)
        pars4ebalance(i)%sunsent_angle = datafiller%data(low+3)
        pars4ebalance(i)%solar_decl = datafiller%data(low+4)
        pars4ebalance(i)%noon = datafiller%data(low+5)
      end do
      
      close(evapconf)
      
      open(newunit=ebalancein, file="drutes.conf/evaporation/ebalance.in", status="old", action="read", iostat = ierr)
      
      if (ierr /= 0) then
        print *, "unable to open file drutes.conf/evaporation/ebalance.in"
        print *, "exiting..."
        ERROR STOP
      end if
      
      call datafiller%clear
      deallocate(tmpdata)
      allocate(tmpdata(7))
      
      counter = 0
      
      do
        call comment(ebalancein)
        read(unit=ebalancein, fmt = *, iostat=ierr) tmpdata
        
        if (ierr /= 0) then
          EXIT
        else
          counter = counter + 1
        end if
   
        do i=1, 7
          call datafiller%fill(tmpdata(i))
        end do
        

      end do
      
      
      allocate(meteo4evap(counter))
      
      if (datafiller%pos /= counter*7) then
        call file_error(ebalancein, "You have defined incorrect number of columns for energy balance time dependent parameters...")
      end if
      
      do i=1, counter
        low = (i-1)*7
        high = i*7
        meteo4evap%time = datafiller%data(low+1)
        meteo4evap(i)%inradiation = datafiller%data(low+2)
        meteo4evap(i)%extraterrad = datafiller%data(low+3)
        meteo4evap(i)%T_air = datafiller%data(low+4)
        meteo4evap(i)%wind_speed = datafiller%data(low+5)
        meteo4evap(i)%cloudiness = datafiller%data(low+6)
        meteo4evap(i)%relhum = datafiller%data(low+7)
      end do
      
    
    end subroutine evapread

end module evapreader
