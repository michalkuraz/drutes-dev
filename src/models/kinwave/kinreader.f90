
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



module kinreader

  private :: read_catchment, gen_catchment
  public :: kinematix

  contains
  
  
    subroutine read_catchment()
      use typy
      use kinglobs
      use readtools
      
      integer :: ierr, filenodes, fileel
      integer(kind=ikind) :: ndcounter, itmp, i
      
      open(newunit=filenodes, file="drutes.conf/kinwave/nodes.in", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open drutes.conf/kinwave/nodes.in, exiting..."
        error stop
      end if
      
      open(newunit=fileel, file="drutes.conf/kinwave/elements.in", status="old", action="read", iostat=ierr)
    
      if (ierr /= 0) then
        print *, "unable to open drutes.conf/kinwave/elements.in, exiting..."
        error stop
      end if
      
      ndcounter = 0
      do 
        call comment(filenodes)
        read(unit=filenodes, fmt=*, iostat=ierr) itmp
        if (ierr /= 0) then
          EXIT 
        end if
        if (itmp > ndcounter) then
          ndcounter = itmp
        end if
      end do
      
      allocate(watershed(ndcounter))
      
      close(filenodes)
      
      open(newunit=filenodes, file="drutes.conf/kinwave/nodes.in", status="old", action="read", iostat=ierr)
      
      print *, ndcounter
      
      do i=1, ndcounter
        call fileread(itmp, filenodes) 
        print *, itmp
      end do
      stop
    
    end subroutine read_catchment
    
    
    subroutine gen_catchment()
      use typy
      use kinglobs
      
!       allocate(
      
      
    
    
    end subroutine gen_catchment
    
    subroutine kinematix()
      use typy
      use pde_objs
      
      select case(drutes_config%mesh_type)
        case(1) 
          call gen_catchment()
          
        case(4)
          call read_catchment()
          
        case default
          print *, "for kinematic wave equation use Arc GIS data option only"
          print *, "  this is option number 4"
          print *, "exiting..."
          ERROR STOP
          
      end select
    
      
    
    end subroutine kinematix

end module kinreader
