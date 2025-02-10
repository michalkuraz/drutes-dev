
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

  private ::  gen_catchment, read_catchment_old, get_slopes
  public :: kininit

  contains
  
  
    subroutine get_slopes()
      use typy
      use kinglobs
      use globals
      use geom_tools
      use debug_tools
      use readtools
      use core_tools
      
      real(kind=rkind), dimension(3) :: a,b,c
      integer(kind=ikind) :: i, j, el
      integer(kind=ikind), dimension(3) :: nd
      integer(kind=ikind), dimension(:), allocatable :: ndpmt
      integer :: file_nodes, ierr
      real(kind=rkind) :: tmp, big
      


      allocate(watershed_el(elements%kolik))
      
      if (drutes_config%mesh_type /= 4) then
        open(newunit=file_nodes, file="drutes.conf/kinwave/nodes.in", action="read", status="old", iostat=ierr)
        
        if (ierr /= 0) then
          print *, "unable to open file drutes.conf/kinwave/nodes.in , check your inputs"
          ERROR STOP
        end if
        allocate(watershed_nd(nodes%kolik))
        
        do i=1, nodes%kolik
          call comment(file_nodes)
          read(unit=file_nodes, fmt=*, iostat=ierr) j, watershed_nd(i)%xyz
          if (abs(watershed_nd(i)%xyz(1)-nodes%data(i,1)) > 1e-3 .or. &
              abs(watershed_nd(i)%xyz(2)-nodes%data(i,2)) > 1e-3) then
            call file_error(file_nodes, "There is a mismatch between drutes.conf/kinwave/nodes.in and mesh file data!")
          end if
        end do
        
        do el=1, elements%kolik
          do j=1,3
            nd(j) = elements%data(el,j)
          end do
          
          a = watershed_nd(nd(1))%xyz
          b = watershed_nd(nd(2))%xyz
          c = watershed_nd(nd(3))%xyz
          
          call plane_derivative(a,b,c, watershed_el(el)%sx, watershed_el(el)%sy)
        end do
        

      else
      
        allocate(ndpmt(maxval(nodes4arcgis%id)))
        
        ndpmt = 0
        
        do i=1, ubound(nodes4arcgis%id,1)
          ndpmt(nodes4arcgis%id(i)) = i
        end do
      
      
        do i=1, ubound(elements4arcgis%data,1)
          if (elements4arcgis%elpermut(i) > 0) then
            el = elements4arcgis%elpermut(i)
            do j=1,3
              nd(j) = ndpmt(elements4arcgis%data(el,j))
            end do
              
            a = nodes4arcgis%data(nd(1),:)
            b = nodes4arcgis%data(nd(2),:)
            c = nodes4arcgis%data(nd(3),:)
            
            call plane_derivative(a,b,c, watershed_el(el)%sx, watershed_el(el)%sy)
          end if
        end do
      end if
      
      
      tmp = 0
      big = 0
      do i=1, elements%kolik
        tmp = tmp + sqrt(watershed_el(i)%sx**2 + watershed_el(i)%sy**2)
        if (sqrt(watershed_el(i)%sx**2 + watershed_el(i)%sy**2) > big) big = sqrt(watershed_el(i)%sx**2 + watershed_el(i)%sy**2)
      end do 
      
        
          
      
    end subroutine get_slopes
    
    
    subroutine read_catchment_old()
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
    
    end subroutine read_catchment_old
    
    
    subroutine gen_catchment()
      use typy
      use kinglobs
      use globals
      use geom_tools
      use debug_tools
      use core_tools
      
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
    
    

    
    subroutine kininit()
      use typy
      use pde_objs
      use kinglobs
      use readtools
      use core_tools
      use geom_tools
      use global_objs
      use debug_tools
      use core_tools
      
      integer :: file_kinematix, ierr, filerain, frainpts
      integer(kind=ikind) :: n, i, counter, j, k
      character(len=512) :: msg
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      real(kind=rkind), dimension(:), allocatable :: pts, distance
      character(len=7), dimension(2) :: infnames
      character(len=4) :: usefile
      
      if (drutes_config%dimen == 2) then
        if (drutes_config%mesh_type == 1) then
          call gen_catchment()
            
        else
          call get_slopes()
        end if
      end if
      
      pde(1)%problem_name(1) = "runoff"
      pde(1)%problem_name(2) = "Kinematic wave equation for real catchments"

      pde(1)%solution_name(1) = "runoff_depth" !nazev vystupnich souboru
      pde(1)%solution_name(2) = "[L] " !popisek grafu

      pde(1)%flux_name(1) = "runoff_flux"  
      pde(1)%flux_name(2) = "runoff flux [L3.T-1.L-2]"
    
      allocate(pde(1)%mass_name(1,2))
      
      pde(1)%mass_name(1,1) = "runoff_height"
      pde(1)%mass_name(1,2) = "runoff elevation [m a.m.s.l]"
      
      if (ubound(pde,1) == 3) then
        pde(2)%problem_name(1) = "debris_flow"
        pde(2)%problem_name(2) = "Transport of solutes for kinematic wave equation"

        pde(2)%solution_name(1) = "debris_concentration" !nazev vystupnich souboru
        pde(2)%solution_name(2) = "[M.L-3] " !popisek grafu

        pde(2)%flux_name(1) = "conc_flux"  
        pde(2)%flux_name(2) = "conc flux [M.L3.T-1.L-2]"
        
        allocate(pde(2)%mass_name(1,2))
        
        pde(2)%mass_name(1,1) = "solute_mass"
        pde(2)%mass_name(1,2) = "solute mass [M]"
      
        pde(3)%problem_name(1) = "concentration_soil"
        pde(3)%problem_name(2) = "Soil contamination"

        pde(3)%solution_name(1) = "soil_concentration" !nazev vystupnich souboru
        pde(3)%solution_name(2) = "[M.M-3] " !popisek grafu

        pde(3)%flux_name(1) = "na"  
        pde(3)%flux_name(2) = "na"
        
        allocate(pde(3)%mass_name(1,2))
        
        pde(3)%mass_name(1,1) = "solid_mass"
        pde(3)%mass_name(1,2) = "solid mass [M]"
      
        pde(3)%problem_name(1) = "concentration_soil"
        pde(3)%problem_name(2) = "Soil contamination"
        
      end if
        
        
      
      open(newunit=file_kinematix, file="drutes.conf/kinwave/kinwave.conf", action="read", status="old", iostat = ierr)
      
      n = maxval(elements%material)
      
      write(msg, fmt=*) "You have specified incorrect number of subregions with different Manning cofficients,", new_line("a"),  &
      "according to mesh definitions, you should define", n, "different Manning values." 
      
      
      call fileread(with_solutes, fileid=file_kinematix, errmsg="set y/n if you want to couple your model with solute transport")
      
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
      
      call fileread(icond_file, errmsg="set y/n for your initial condition input method  &
                  (n=use separate file with node wise defined values (drutes.conf/kinwave/initcond, y=define it for each material &
                  in drutes.conf/kinwave/kinwave.conf)", fileid=file_kinematix)
                

      if (with_solutes) then
        allocate(kinsols(n))
        deallocate(tmp_array)
        allocate(tmp_array(8))
        do i=1,n
          print *, i
          call fileread(tmp_array, fileid=file_kinematix, checklen=.true.)
          kinsols(i)%horb = tmp_array(1)
          if (icond_file) then
            if (tmp_array(2) /= -99) then
              call file_error(file_kinematix, errmsg="on line above you speficified initial condition to be supplied node-wise in &
                  drutes.conf/kinwave/initcond, if this is really what you want to do, then set of cs_init on material description &
                  value -99")
            end if
          end if
                  
          kinsols(i)%csinit = tmp_array(2)
          kinsols(i)%rhos = tmp_array(3)
          kinsols(i)%lambda_a = tmp_array(4)
          kinsols(i)%lambda_d = tmp_array(5)
          kinsols(i)%lost_me = tmp_array(6)
          kinsols(i)%ncs = tmp_array(7)
          kinsols(i)%diff =  tmp_array(8)
        end do
      end if
            
      
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
