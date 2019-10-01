
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
  public :: kininit

  contains
  
  
    subroutine read_catchment()
      use typy
      use kinglobs
      use readtools
      
      integer :: ierr, filenodes, fileel
      integer(kind=ikind) :: ndcounter, itmp, i
      real(kind=rkind), dimension(:), allocatable :: distance
      real(kind=rkind), dimension(:), allocatable :: tmpvals
      
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
      
      allocate(watershed_nd(ndcounter))
      
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
      use globals
      use geom_tools
      use debug_tools
      
      integer(kind=ikind) :: i
      real(kind=rkind), dimension(3) :: a,b,c
      
      allocate(watershed_nd(nodes%kolik))
      
      allocate(watershed_el(elements%kolik))
      
      
      select case(drutes_config%dimen)
        case(2)
          do i=1, nodes%kolik
            watershed_nd(i)%xyz(1:2) = nodes%data(i,:)
            watershed_nd(i)%xyz(3) = nodes%data(i,2)**2*0.001 + nodes%data(i,2)**2*0.02
          end do
          do i=1, elements%kolik
            a = watershed_nd(elements%data(i,1))%xyz
            b = watershed_nd(elements%data(i,2))%xyz
            c = watershed_nd(elements%data(i,3))%xyz
        
            call plane_derivative(a,b,c, watershed_el(i)%sx, watershed_el(i)%sy)
          end do
        case(1)
     
          do i=1, elements%kolik
            watershed_el(i)%sx = oneDslopes(elements%material(i))
          end do
          
      end select
          
    
      
      watershed_el(:)%cover = 1
      


    
    end subroutine gen_catchment
    
    

    
    subroutine kininit(pde_loc)
      use typy
      use pde_objs
      use kinglobs
      use readtools
      use core_tools
      use geom_tools
      use global_objs
      use debug_tools
      
      class(pde_str), intent(in out) :: pde_loc
      integer :: file_kinematix, ierr, filerain, frainpts
      integer(kind=ikind) :: n, i, counter, j, k
      character(len=512) :: msg
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      real(kind=rkind), dimension(:), allocatable :: pts, distance
      character(len=7), dimension(2) :: infnames
      
      if (drutes_config%dimen == 2) then
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
      end if
      
      pde_loc%problem_name(1) = "runoff"
      pde_loc%problem_name(2) = "Kinematic wave equation for real catchments"

      pde_loc%solution_name(1) = "runoff_depth" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "[L] " !popisek grafu

      pde_loc%flux_name(1) = "runoff_flux"  
      pde_loc%flux_name(2) = "runoff flux [L3.T-1.L-2]"
    
      allocate(pde_loc%mass_name(1,2))
      
      pde_loc%mass_name(1,1) = "runoff_height"
      pde_loc%mass_name(1,2) = "runoff elevation [m a.m.s.l]"
      
      open(newunit=file_kinematix, file="drutes.conf/kinwave/kinwave.conf", action="read", status="old", iostat = ierr)
      
      n = maxval(elements%material)
      
      write(msg, fmt=*) "You have specified incorrect number of subregions with different Manning cofficients,", new_line("a"),  &
      "according to mesh definitions, you should define", n, "different Manning values." 
      
      
      allocate(manning(n))
      
      allocate(oneDslopes(n))
      
      allocate(inf_model(n))
      
      call fileread(i, file_kinematix, ranges=(/n,n/), errmsg=cut(msg))
      
      do i=1, n
        call fileread(manning(i), file_kinematix, ranges=(/epsilon(manning(1)), huge(manning(1))/))
      end do
      
      call fileread(backwater, file_kinematix)
      
      if (drutes_config%dimen == 1) then
        do i=1, n
          call fileread(oneDslopes(i), file_kinematix, ranges=(/-huge(oneDslopes(1)), huge(oneDslopes(1))/))
        end do
        call gen_catchment()
      end if
      
      infnames(1) = "Ks"
      infnames(2) = "Schwarz"
      
      write(msg, *) "Set model for reducing rainfall by infiltration. ", &
            "The model must be set for each subregion, thus the number of subregions -> ", & 
            "number of lines with model type records.", new_line("a"), &
            "-------------------------------------", &
             new_line("a"), new_line("a"),&
            "    Ks - only saturated hydraulic conductivity is subtracted",   new_line("a"), new_line("a"),&
            "    Schwarz - Schwarz equation for unsteady infiltration is used,",  new_line("a"),& 
            "    see: D. Swartzendruber. ", & 
            "A Quasi-Solution of Richards Equation for the Downward Infiltration of Water into Soil. " , &  
            "Water Resour. Res., 23(5): 809–817, MAY 1987."
            
      do i=1, ubound(inf_model,1)
        call fileread(inf_model(i)%name, file_kinematix, cut(msg), options=infnames)
      end do
      
      do i=1, ubound(inf_model,1)
        if (allocated(tmp_array)) deallocate(tmp_array)
        select case(cut(inf_model(i)%name))
          case("Ks")
            write(msg,*) "For subregion:", i, "model Ks for reducing the rainfall intensity is defined", &
            "provide only a SINGLE value of Ks"
            allocate(tmp_array(1))
            call fileread(tmp_array, errmsg=msg, fileid=file_kinematix, checklen=.true.)
            inf_model(i)%Ks = tmp_array(1)
          case("Schwarz")
             write(msg,*) "For subregion:", i, &
             "model Schwarz (Schwarzendruber, 1987) for reducing the rainfall intensity is defined", &
             "provide THREE values for S, A, Ks"
             allocate(tmp_array(3))
             call fileread(tmp_array, errmsg=msg, fileid=file_kinematix, checklen=.true.)
             inf_model(i)%S = tmp_array(1)
             inf_model(i)%A = tmp_array(2)
             inf_model(i)%Ks = tmp_array(3)
        end select
      end do
            
      
      open(newunit=frainpts, file="drutes.conf/kinwave/rain.pts", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open file with rainfall data drutes.conf/kinwave/rain.pts"
        error stop
      end if
      
      
      counter = 0
      do
        call comment(frainpts)
        read(unit=frainpts, fmt=*, iostat=ierr) tmp
        
        if (ierr /= 0) then
          EXIT
        else
          counter = counter + 1
        end if
      end do
      
      close(frainpts)
      
      open(newunit=frainpts, file="drutes.conf/kinwave/rain.pts", status="old", action="read", iostat=ierr)
      
      allocate(raindata(counter))
      
      do i=1, ubound(raindata,1)
        allocate(raindata(i)%xy(drutes_config%dimen))
        call fileread(raindata(i)%xy, frainpts, checklen=.true.)
      end do
      
      allocate(el2pt(elements%kolik))
      
      allocate(pts(drutes_config%dimen))
      
      allocate(distance(ubound(raindata,1)))
      do i=1, elements%kolik
        distance = huge(distance)
        do j=1, ubound(distance,1)
          pts = 0
          do k=1, ubound(elements%data,2)
            pts = pts + nodes%data(elements%data(i,k),:)
          end do
          pts = pts/ubound(elements%data,2)
          distance(j) = dist(pts, raindata(j)%xy)
        end do
        el2pt(i) = minloc(distance,1)        
      end do
      
      
      open(newunit=filerain, file="drutes.conf/kinwave/rain.in", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then 
        print *, "unable to open drutes.conf/kinwave/rain.in, missing file"
        ERROR STOP
      end if
      
      if (allocated(tmp_array)) deallocate(tmp_array)
      allocate(tmp_array(ubound(raindata,1)+1))
      
      counter = 0
      do 
        call comment(filerain)
        read(unit=filerain, fmt=*, iostat=ierr) tmp_array
        if (ierr == 0) then
          counter = counter + 1
        else
          EXIT
        end if
      end do
      
      
      do i=1, ubound(raindata,1)
        allocate(raindata(i)%series(counter,2))
      end do
      
      close(filerain)
      
      open(newunit=filerain, file="drutes.conf/kinwave/rain.in", status="old", action="read", iostat=ierr)
      
      do i=1, ubound(raindata(1)%series,1)
        call fileread(tmp_array, filerain, checklen=.true.)
        do j=1, ubound(raindata,1)
          raindata(j)%series(i,1) = tmp_array(1)
          raindata(j)%series(i,2) = tmp_array(1+j)
        end do
      end do
      

    end subroutine kininit

end module kinreader
