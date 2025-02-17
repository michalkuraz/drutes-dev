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

!> \file postpro.f90
!! \brief Output files generator. Scilab friendly.
!<
module postpro
  use typy
  
  public :: make_print
  public :: write_obs , write_bcfluxes
  public :: get_RAM_use
  private :: print_scilab, print_pure, print_gmsh, write_bcfluxes_old
  


  contains
 
    subroutine make_print(behaviour, curtime, name)
      use typy
      use globals
      use global_objs
      use core_tools
      use geom_tools
      use pde_objs
      use debug_tools

      character(len=*), intent(in), optional                :: behaviour
      real(kind=rkind), intent(in), optional                :: curtime
      character(len=*), intent(in), optional                :: name
      logical                                               :: anime, op
      integer(kind=ikind)                                   :: mode, no_prints
      type :: filenames_str
      character(len=256), dimension(:), allocatable       :: names
      end type filenames_str
      type(filenames_str), dimension(:), allocatable        :: filenames
      integer(kind=ikind)                                   :: i, i_err, j, layer, proc, run, ii
      character(len=64)                                     :: forma
      integer, dimension(:,:), allocatable, save            :: ids
      integer, dimension(:,:), allocatable,  save           :: ids_obs
      integer, dimension(:,:), allocatable,  save           :: ids_anime
      character(len=5)                                      :: extension
      character(len=15)                                     :: prefix
      real(kind=rkind)                                      :: distance, avgval, val1, val2, val3, tmp, flux
      logical, save                                         :: first_run=.true., first_anime=.true., openedq
      type(integpnt_str)                                    :: quadpnt
      integer(kind=ikind), save                             :: anime_run, anime_dec
      integer                                               :: ierr
      integer(kind=ikind), dimension(:), allocatable, save  :: masscnt, fluxcnt, no_files
      character(len=256) :: name_of_file

      if (present(name)) then
        select case(name)
          case("obs")
            anime = .false.
            mode = 0
          case("avi")
          	if (drutes_config%dimen == 1) then
							print *, "animations are supported for 2 dimensional domains only"
							print *, "your setup is 1D"
							print *, "exited from postpro:: make_print"
							ERROR STOP
          	end if
            anime = .true.
            mode = -1
          case default
            print *, "this is a BUG, contact developer, called from postpro::make_print()"
            error STOP
        end select
      else
        anime = .false.
        mode = 0
      end if
      
      if (cut(observe_info%fmt) == "gmsh") then 
        mode = -1
      end if
      

      if (drutes_config%dimen < 2  .or. www) then
        extension = ".dat"
      else 
        select case(observe_info%fmt)
          case("pure")
            extension = ".dat"
          case("scil")
            extension = ".sci"
          case("gmsh")
            extension = ".msh"
        end select
      end if
	  
  
      call write_log("making output files for observation time")
      

      if (.not. anime) then
        postpro_run = postpro_run + 1

      else
        if (first_anime) then
          first_anime = .false.
          anime_run = 0
          ierr=system("mkdir out/anime")
        end if
        anime_run = anime_run + 1
      end if
      
      if (.not. anime) then
        do 
          if (postpro_run/(10**postpro_dec) < 1) then
            EXIT
          else
            postpro_dec = postpro_dec + 1
          end if
        end do
      end if

      i=1
      
      if (.not. allocated(masscnt)) then
        allocate(masscnt(ubound(pde,1)))
        allocate(fluxcnt(ubound(pde,1)))
        allocate(no_files(ubound(pde,1)))
        masscnt = 0
        fluxcnt = 0
    
        do proc=1, ubound(pde,1)
          if (allocated(pde(proc)%mass_name) ) then
            masscnt(proc) = ubound(pde(proc)%mass_name,1)
          end if

          if (allocated(pde(proc)%fluxes)) then
            fluxcnt(proc) = ubound(pde(proc)%fluxes,1)
          end if

      
          if (fluxcnt(proc) == 0) then
            no_files(proc) = masscnt(proc) + 3
          else
            no_files(proc) = masscnt(proc) + fluxcnt(proc) + 2
          end if
        end do
        

     
        allocate(ids(ubound(pde,1), maxval(no_files)))
        allocate(ids_obs(ubound(pde,1), maxval(no_files)))
        allocate(ids_anime(ubound(pde,1), maxval(no_files)))
        ids_obs = 0
        ids_anime = 0
      end if
      
      if (anime) then
        ids = ids_anime
      else
        ids = ids_obs
      end if
     

      allocate(filenames(ubound(pde,1)))
      
      do i=1, ubound(filenames,1)
        if (.not. allocated(pde(i)%fluxes)) then
          allocate(filenames(i)%names(3+ubound(pde(i)%mass_name,1)))
        else
          allocate(filenames(i)%names(2+ubound(pde(i)%mass_name,1)+ubound(pde(i)%fluxes,1)))
        end if
      end do
      

      if (anime) then
        prefix = "out/anime/"
        write(unit=forma, fmt="(a, I7, a)") "(a, a, a, a, a, I", 1," a)"
        run = anime_run
      else
        prefix = "out/"
        write(unit=forma, fmt="(a, I7, a)") "(a, a, a, a, a, I", postpro_dec," a)"
        run = postpro_run
      end if


      do proc=1, ubound(pde,1)

        write(unit=filenames(proc)%names(1), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", &
        trim(pde(proc)%solution_name(1)), "-",  run, trim(extension)


        write(unit=filenames(proc)%names(2), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", & 
                    trim(pde(proc)%solution_name(1)), &
                    "-el_avg-",  run,  trim(extension)

        if (ubound(pde(proc)%mass_name,1) > 0 ) then
          do i=1, ubound(pde(proc)%mass_name,1)
            write(unit=filenames(proc)%names(2+i), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_",  &
              cut(pde(proc)%mass_name(i,1)), "-", run,  trim(extension)
          end do
        end if

        if (.not. allocated(pde(proc)%fluxes)) then
        
          write(unit=filenames(proc)%names(3+ubound(pde(proc)%mass_name,1)), fmt=forma) trim(prefix),&
                       trim(pde(proc)%problem_name(1)), "_", &
                       trim(pde(proc)%flux_name(1)), "-", &
                      run,  trim(extension)
        else
          do i=1, ubound(pde(proc)%fluxes,1)
             write(unit=filenames(proc)%names(2+i+ubound(pde(proc)%mass_name,1)), fmt=forma) trim(prefix),&
                       trim(pde(proc)%problem_name(1)), "_", &
                       trim(pde(proc)%fluxes(i)%name(1)), "-", &
                      run,  trim(extension)
          end do
        end if
        


        if ( (.not. anime .and. mode == 0)  .or. &
          (anime_run == 1 .and. anime) .or. & 
         ( .not. anime .and. (mode == -1 .and. postpro_run == 0 ) ) ) then


          
          do i=1, ubound(filenames(proc)%names,1)
              open(newunit=ids(proc, i), file=trim(filenames(proc)%names(i)), action="write", status="replace", iostat=ierr)
          end do
          
          ids_obs = ids
          
          if (anime .and. anime_run == 1) then
            ids_anime = ids
          end if
        end if
        


        quadpnt%type_pnt = "ndpt"
        
        if (time > epsilon(time) ) then
          quadpnt%column = 3
        else
          quadpnt%column = 1
        end if
              
        if (present(curtime)) then
          quadpnt%globtime = .false.
          quadpnt%time4eval = curtime
        end if

   
        if (drutes_config%dimen == 1 .or. www) then

          call print_pure(ids(proc,:), proc, quadpnt)

        else
          select case(cut(observe_info%fmt))
            case("pure", "csv")
              call print_pure(ids(proc,1:no_files(proc)), proc, quadpnt)
            case("scil")
              call print_scilab(ids(proc,1:no_files(proc)), proc, quadpnt)
            case("gmsh")
              call print_gmsh(ids(proc,1:no_files(proc)), proc, quadpnt)
          end select
        end if
      end do


           
      if (.not. anime) then
        do proc=1, ubound(pde,1)
          do i=1,ubound(ids,2)
            if (i /= 3 .or. pde(proc)%print_mass) then
              if (mode == 0) then
                if (ids(proc,i) /= 0) close(ids(proc,i))
              else
                 if (ids(proc,i) /= 0) call flush(ids(proc,i))
              end if
            end if
          end do
        end do	
      else      
        do proc=1, ubound(pde,1)
          do i=1,ubound(ids,2)
           if (i /= 3 .or. pde(proc)%print_mass) then
              write(unit=ids(proc,i), fmt="(a,I6.6,a)" ) "xs2png(0, 'K-", anime_run , ".png');"
              write(unit=ids(proc,i), fmt=*) "clear"
              write(unit=ids(proc,i), fmt=*) "   "
              call flush(ids(proc,i))
            end if
          end do
        end do


      end if
      

      deallocate(filenames)
      
      first_run=.false.
      
      
    end subroutine make_print




    subroutine write_obs()
      use typy
      use globals
      use global_objs
      use geom_tools
      use pde_objs
      use debug_tools
      use core_tools
      
      integer(kind=ikind) :: i, layer, proc, D, j, massdim, printdim, fluxid
      real(kind=rkind) :: val
      real(kind=rkind), dimension(:), allocatable, save :: massval
      real(kind=rkind), dimension(3) :: advectval
      real(kind=rkind), dimension(:,:), allocatable :: extrafluxes
      type(integpnt_str) :: quadpnt
      character(len=1) :: csv
!      real(kind=rkind), dimension(:), allocatable, save :: volume
!      integer, dimension(:), allocatable, save :: volumeunits
      real(kind=rkind) :: tmp

      quadpnt%type_pnt = "obpt"
      quadpnt%column=3
      D = drutes_config%dimen
      
      if (cut(observe_info%fmt) == "csv" ) then
        csv = ","
      else
        csv = " "
      end if
      
      massdim=0
      if (.not. allocated(massval)) then
        do proc=1, ubound(pde,1)
          if (ubound(pde(proc)%mass,1) > massdim) massdim = ubound(pde(proc)%mass,1)
        end do
        allocate(massval(massdim))
      end if
      

      do proc=1, ubound(pde,1)
        do i =1, ubound(observation_array,1)
          layer = elements%material(observation_array(i)%element)

          quadpnt%order = i
          
          quadpnt%element = observation_array(i)%element
          
          quadpnt%preproc=.true.

          if (.not. allocated(pde(proc)%fluxes)) then
            call pde(proc)%flux(layer=layer, quadpnt=quadpnt,  vector_out=advectval(1:D))
            observation_array(i)%cumflux(proc) = observation_array(i)%cumflux(proc) + norm2(advectval(1:D))*dtprev
          else
            if (.not. allocated(extrafluxes)) allocate(extrafluxes(ubound(pde(proc)%fluxes,1),D+1))
            do fluxid=1, ubound(pde(proc)%fluxes,1)
              call pde(proc)%fluxes(fluxid)%val(pde(proc),layer=layer, quadpnt=quadpnt, vector_out=extrafluxes(fluxid,1:D))
              pde(proc)%fluxes(fluxid)%cumval = pde(proc)%fluxes(fluxid)%cumval + norm2(extrafluxes(fluxid,1:D))*dtprev
              extrafluxes(fluxid, D+1) = pde(proc)%fluxes(fluxid)%cumval
            end do
          end if
          
          
          
            
          
          if (pde(proc)%print_mass) then
            printdim  = ubound(pde(proc)%mass,1)
          else
            printdim = 0
          end if
          
          quadpnt%preproc = .true.
          val = pde(proc)%getval(quadpnt)
          
          do j=1, ubound(pde(proc)%mass,1)
            massval(j) = pde(proc)%mass(j)%val(pde(proc),layer, quadpnt)
          end do
          
          if (.not. allocated(pde(proc)%fluxes)) then
            write(unit=pde(proc)%obspt_unit(i), fmt=*) time, csv, val, csv,  massval(1:printdim), csv, advectval(1:D), &
                    csv, observation_array(i)%cumflux(proc)
          else
            fluxid = 1
            write(unit=pde(proc)%obspt_unit(i), fmt=*) time, csv, val, csv,  massval(1:printdim), csv, &
                (/ ( extrafluxes(fluxid,:), fluxid=1, ubound(extrafluxes,1)) /)
          end if

          call flush(pde(proc)%obspt_unit(i))
        end do
        if (allocated(extrafluxes)) deallocate(extrafluxes)
      end do
      
!      if (.not. allocated(volume)) then

!        allocate(volume(1))
!        allocate(volumeunits(ubound(pde,1)))
!        open(newunit=volumeunits(1), file="out/volume-re.dat", action="write", status="new")
          
!      end if  
      
!      volume(1) = 0
!      quadpnt%type_pnt = "ndpt"
      
!      do i=1, elements%kolik
!        tmp = 0
!        do j=1,3
!          quadpnt%order = elements%data(i,j)
!          tmp = tmp + pde(1)%mass(1)%val(pde(1),layer, quadpnt)
!        end do
!        tmp = tmp/3.0
!        volume(1) = volume(1) + tmp*elements%areas(i)
!      end do
      
      
!      write(volumeunits(1), fmt=*) time, volume(1)
      
!      call flush(volumeunits(1))
      
    end subroutine write_obs
    

   subroutine write_bcfluxes()
      use typy
      use globals
      use geom_tools
      use pde_objs
      use global_objs
      use core_tools
      use debug_tools
      
      
      integer(kind=ikind) :: bc, el, D, i, proc, surfel
      real(kind=rkind), dimension(:,:), allocatable, save :: flux
      real(kind=rkind), dimension(2,2) :: bcpt
      type(integpnt_str) :: quadpnt
      real(kind=rkind) :: length
      real(kind=rkind), dimension(:,:), allocatable, save :: integflux, cumflux
      integer(kind=ikind), dimension(:), allocatable, save :: surfnodes
      
      D = drutes_config%dimen
      
      quadpnt%type_pnt="xypt"

      quadpnt%column=2
      
      
      if (.not. allocated(flux)) then
        allocate(flux(ubound(pde,1), D))
        allocate(surfnodes(D))
    
        allocate(integflux(lbound(bcfluxes,1) : ubound(bcfluxes,1) , ubound(pde,1)))
        allocate(cumflux(lbound(bcfluxes,1) : ubound(bcfluxes,1) , ubound(pde,1)))
        cumflux = 0
      end if
      

      integflux = 0
      do bc = lbound(bcfluxes,1), ubound(bcfluxes,1)
        do surfel = 1, bcfluxes(bc)%pos
          el = bcfluxes(bc)%bcel(surfel)%element
          quadpnt%element = el
          surfnodes = bcfluxes(bc)%bcel(surfel)%surfnode(1:D)
          quadpnt%xy(1:D) = mean(nodes%data(surfnodes,:))
          
          do proc=1, ubound(pde,1)
            call pde(proc)%flux(elements%material(el), quadpnt, vector_out=flux(proc,:))
            integflux(bc, proc) = integflux(bc, proc) - &
                                  dot_product(flux(proc,:), bcfluxes(bc)%bcel(surfel)%n_out(1:D))* &
                                                                   bcfluxes(bc)%bcel(surfel)%area
          end do
        end do
      end do
            

      cumflux = cumflux + integflux*time_step          
      
      
      
      do bc = lbound(bcfluxes,1), ubound(bcfluxes,1)
        write(unit=bcfluxes(bc)%fileid, fmt=*) time, integflux(bc, :), cumflux(bc,:)
        call flush(bcfluxes(bc)%fileid)
      end do
    
    end subroutine write_bcfluxes
    
    
    subroutine write_bcfluxes_old()
      use typy
      use globals
      use geom_tools
      use pde_objs
      use global_objs
      use debug_tools
      
      integer(kind=ikind) :: bc, el, D, i, proc, surfel
      real(kind=rkind), dimension(:,:), allocatable, save :: flux, fluxtmp
      real(kind=rkind), dimension(2,2) :: bcpt
      type(integpnt_str) :: quadpnt
      real(kind=rkind) :: length
      real(kind=rkind), dimension(:,:), allocatable, save :: integflux, cumflux
      
      D = drutes_config%dimen
      
      quadpnt%type_pnt="ndpt"

      quadpnt%column=2
      
      
      if (.not. allocated(flux)) then
        allocate(flux(ubound(pde,1), D))
        allocate(fluxtmp(ubound(pde,1), D))
    
        allocate(integflux(lbound(bcfluxes,1) : ubound(bcfluxes,1) , ubound(pde,1)))
        allocate(cumflux(lbound(bcfluxes,1) : ubound(bcfluxes,1) , ubound(pde,1)))
        cumflux = 0
      end if
      

      integflux = 0
      do bc = lbound(bcfluxes,1), ubound(bcfluxes,1)
        do surfel = 1, bcfluxes(bc)%pos
          do proc=1, ubound(pde,1)
            flux = 0
            do i=1, D
              quadpnt%order=bcfluxes(bc)%bcel(surfel)%surfnode(i)
              call pde(proc)%flux(elements%material(el), quadpnt, vector_out=fluxtmp(proc,:))
              flux(proc,:) = flux(proc,:) + fluxtmp(proc,:)
            end do
            flux(proc,:) = flux(proc,:)/D
            integflux(bc, proc) = integflux(bc, proc) + &
                                  dot_product(flux(proc,:), bcfluxes(bc)%bcel(surfel)%n_out(1:D))* &
                                                            bcfluxes(bc)%bcel(surfel)%area
          end do
        end do
      end do
            

      cumflux = cumflux + integflux*time_step          
      
      
      
      do bc = lbound(bcfluxes,1), ubound(bcfluxes,1)
        write(unit=bcfluxes(bc)%fileid, fmt=*) time, integflux(bc, :), cumflux(bc,:)
        call flush(bcfluxes(bc)%fileid)
      end do
    
    end subroutine write_bcfluxes_old
 

    !> this subroutine parses /proc/[PID]/status file in order to get RAM consumption statistics, it works Linux only
    subroutine get_RAM_use()
      use typy
      use globals
      use core_tools
      integer :: PID, fileid, i, ierr
      integer(kind=ikind) :: bytes
      character(len=2) :: byte_unit
      character(len=7) :: ch
      character(len=256) :: filename, format, command

      PID = getpid()

      i = 0
      do
        i = i + 1
        if ((1.0*PID)/(10**i) < 1) then
        EXIT
        end if
      end do

      write(unit=format, fmt="(a, I7, a)") "(a, I", i, ", a)"

  !     write(unit=format, fmt = *) "(I.", i, ")"

      write(unit=filename, fmt=format) "/proc/", PID, "/status"
      
      write(unit=command, fmt=*) "cp ", cut(filename), " out/"
      
      ierr = system(cut(command))
      
      open(newunit=fileid, file=filename, action="read", status="old", iostat=ierr)

      if (ierr /= 0) then
        write(unit=terminal) "WARNING! this is not POSIX system, unable to get RAM consumption"
        RETURN
      end if

      do 
        read(unit=fileid, fmt=*, iostat=ierr) ch
        if (ch == "VmPeak:") then
          backspace fileid
          EXIT
        end if

        if (ierr /=0) then
          print *, "unable to fetch memory consumption from system files"
          RETURN
        end if
      end do

      read(unit=fileid, fmt=*, iostat=ierr) ch, bytes, byte_unit
      
      if (ierr /=0) then
        print *, "unable to fetch memory consumption from system files"
        RETURN
      end if

      call write_log(text="Peak RAM  consumption on image", int1=1_ikind*THIS_IMAGE(), text2="was:", int2=bytes, text3=byte_unit)

      do 
        read(unit=fileid, fmt=*) ch
        if (ch == "VmSwap:") then
          backspace fileid
          EXIT
        end if

        if (ierr /=0) then
          print *, "unable to fetch swap consumption from system files"
          RETURN
        end if
      end do

      read(unit=fileid, fmt=*, iostat=ierr) ch, bytes, byte_unit
      
      if (ierr /=0) then
        print *, "unable to fetch swap consumption from system files"
        RETURN
      end if

      call write_log(text="Peak SWAP consumption on image", int1=1_ikind*THIS_IMAGE(), text2="was:", int2=bytes, text3=byte_unit)


      close(fileid)

    end subroutine get_RAM_use
    
    subroutine print_scilab(ids, proc, quadpnt)
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use debug_tools
      use core_tools
      
      integer, dimension(:), intent(in) :: ids
      integer(kind=ikind), intent(in) :: proc
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind) :: i, j, layer, ii
      real(kind=rkind) :: tmp, totflux
      real(kind=rkind), dimension(3) :: flux
      real(kind=rkind), dimension(:,:), allocatable, save :: body
      real(kind=rkind), dimension(2) :: vct1, vct2
      real(kind=rkind), dimension(:), allocatable, save :: vct_tmp
      integer(kind=ikind) :: time_dec
      real(kind=rkind) :: curtime
      type(integpnt_str) :: qpntloc
      type(integpnt_str) :: quadpnt_loc
      
!      logical :: opened
!      character(len=256) :: filename
      
      
!      do i=1, ubound(ids,1)
!        INQUIRE(UNIT=ids(i), NAME=filename, OPENED=opened)
!        print *, i, ids(i), cut(filename), opened
!      end do
      
!      call wait()

      if (.not. allocated(body)) allocate(body(3, 7+ubound(pde(proc)%mass,1)))
      
      if (.not. allocated(vct_tmp)) allocate(vct_tmp(7+ubound(pde(proc)%mass,1)))
      
      quadpnt_loc = quadpnt
      if (.not. quadpnt%globtime) then
        curtime = quadpnt%time4eval
      else
        curtime = time
      end if
    
      do i=1, ubound(ids,1)
        write(unit=ids(i), fmt=*) "//", curtime
        write(unit=ids(i), fmt=*) "nt =", elements%kolik, ";"
        write(unit=ids(i), fmt=*) "x=zeros(nt,3);"
        write(unit=ids(i), fmt=*) "y=zeros(nt,3);"
        write(unit=ids(i), fmt=*) "z=zeros(nt,3);"
      end do
      



      quadpnt_loc%preproc=.true.
      do i=1, elements%kolik
        do j=1,ubound(elements%data,2)
          body(j,1:2) = nodes%data(elements%data(i,j),:)
          quadpnt_loc%order = elements%data(i,j)
          body(j,3) = pde(proc)%getval(quadpnt_loc)
        end do
      
        ! mass (constant over element)
        layer = elements%material(i)     
      
        qpntloc%element = i
        qpntloc%column = 3
        qpntloc%type_pnt = "gqnd"
        qpntloc%preproc=.true.
        do ii=1, ubound(pde(proc)%mass,1)
          tmp = 0
          do j=1, ubound(gauss_points%weight,1)
            qpntloc%order = j
            tmp = tmp + pde(proc)%mass(ii)%val(pde(proc),layer, qpntloc)*gauss_points%weight(j)
          end do
        
          body(:,4+ii) = tmp/gauss_points%area
        end do

          
        if (ubound(gauss_points%weight,1) > 1) then
          qpntloc%order = nint(ubound(gauss_points%weight,1)/2.0)
        else
          qpntloc%order = 1
        end if
        
        call pde(proc)%flux(layer, quadpnt_loc, vector_out=flux(1:drutes_config%dimen), scalar=totflux)

        body(:,5+ubound(pde(proc)%mass,1)) = totflux
        
        vct1 = body(3, 1:2) - body(1, 1:2) 
        vct2 = body(3, 1:2) - body(2, 1:2)
        
        if ( (vct1(1)*vct2(2)-vct1(2)*vct2(1)) > 0.0_rkind) then
          vct_tmp = body(3,:)
          body(3,:) = body(2,:)
          body(2,:) = vct_tmp
        else
          CONTINUE
        end if

        do j=1, ubound(ids,1)
        

          write(unit=ids(j), fmt=*) "x(", i, ",1) =", body(1,1), ";"
          write(unit=ids(j), fmt=*) "x(", i, ",2) =", body(2,1), ";"
          write(unit=ids(j), fmt=*) "x(", i, ",3) =", body(3,1), ";"
          
          write(unit=ids(j), fmt=*) "y(", i, ",1) =", body(1,2), ";"
          write(unit=ids(j), fmt=*) "y(", i, ",2) =", body(2,2), ";"
          write(unit=ids(j), fmt=*) "y(", i, ",3) =", body(3,2), ";"

          write(unit=ids(j), fmt=*) "z(", i, ",1) =", body(1,2+j), ";"
          write(unit=ids(j), fmt=*) "z(", i, ",2) =", body(2,2+j), ";"
          write(unit=ids(j), fmt=*) "z(", i, ",3) =", body(3,2+j), ";"


            
        end do

      end do
    
      time_dec = 0
      if (curtime>epsilon(curtime)) then
        do
          if (curtime*10.0_rkind**time_dec > 1) then
            EXIT 
          else
            time_dec = time_dec + 1
          end if
        end do
      end if
      


      do i=1, ubound(ids,1)
      

        write(unit=ids(i), fmt=*) "f=gcf();"
        write(unit=ids(i), fmt=*) "clf(f,'reset');"
        write(unit=ids(i), fmt=*) "f.color_map=jetcolormap(256);"
        write(unit=ids(i), fmt=*) "colorbar(min(z),max(z));"
!      write(unit=ids(i), fmt=*) "colorbar(0,0.5);"
        if (time_dec < 2) then
          write(unit=ids(i), fmt="(a,F10.2,a,a,a,a)")  "xtitle('$\mathbf{\LARGE t= ", curtime,"  ",&
          "} \quad \mbox{\Large ",   trim(time_units), "}$')"
        else
   !                                                             a	                 F10.2        a
         write(unit=ids(i), fmt="(a,F10.2,a,a,I16,a,a,a)")  "xtitle('$\mathbf{\LARGE t= ", curtime*10**time_dec,"  ",&
  !                    a            I16         a                         a                a
          "\times 10^{-", time_dec,"}} \quad \mbox{\Large ",   trim(time_units), "}$')"
        end if
        
        
        write(unit=ids(i), fmt=*) "plot3d1(x',y',z',alpha=0, theta=-90);"


      end do 
    
    end subroutine print_scilab
    
    
    subroutine print_pure(ids, proc, quadpnt)
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      use core_tools
      
      integer, dimension(:), intent(in) :: ids
      integer(kind=ikind), intent(in) :: proc
      type(integpnt_str),  intent(in out) :: quadpnt
      real(kind=rkind) :: curtime, tmp
      integer(kind=ikind) :: i, layer, j, velocity_id, ii, jj, D, iflux
      real(kind=rkind) ::  distance,  avgval
      type(integpnt_str) :: qpntloc
      real(kind=rkind), dimension(3) :: flux, flux_tmp
      real(kind=rkind), dimension(:,:), allocatable :: fluxes
      real(kind=rkind), dimension(:), allocatable :: massval
      real(kind=rkind), dimension(:), allocatable, save :: pts
      character(len=1) :: csv

      if (cut(observe_info%fmt) == "csv") then
        csv = ","
      else
        csv = " "
      end if
      

 
      if (.not. allocated(pts)) allocate(pts(ubound(nodes%data,2)))
      
      D = drutes_config%dimen
      
      do i=1, nodes%kolik
        quadpnt%type_pnt = "ndpt"
        quadpnt%order = i
        quadpnt%preproc=.true.
        write(unit=ids(1), fmt=*) i, csv,  nodes%data(i,:), csv,  pde(proc)%getval(quadpnt) 
      end do
      
      
      do i=1, elements%kolik
        avgval=0
        do ii=1, ubound(elements%data,2)
          quadpnt%type_pnt = "ndpt"
          quadpnt%preproc = .true.
          quadpnt%order = elements%data(i,ii)
          avgval = avgval + pde(proc)%getval(quadpnt) 
        end do
        avgval = avgval/ubound(elements%data,2)
        
   
        pts = 0
        do ii=1, ubound(elements%data,2)
          pts = pts + nodes%data(elements%data(i,ii),:)
        end do
        pts = pts/ubound(elements%data,2)

        write(unit=ids(2), fmt=*)  i, csv, pts, csv, avgval

        quadpnt%element = i
        quadpnt%column = 3
        quadpnt%type_pnt = "gqnd"
        quadpnt%preproc=.true. 
        
        layer = elements%material(i)
        
        if (allocated(pde(proc)%fluxes)) then
          if (.not. allocated(fluxes)) allocate(fluxes(ubound(pde(proc)%fluxes,1),D))
        end if

     
        ! 3 is for mass
        if (pde(proc)%print_mass) then
          do j=3, ubound(pde(proc)%mass,1) + 2
            tmp = 0
            do jj=1, ubound(gauss_points%weight,1)
              quadpnt%order = jj
              tmp = tmp + pde(proc)%mass(j-2)%val(pde(proc),layer, quadpnt)*gauss_points%weight(jj)
            end do
            tmp = tmp/gauss_points%area
            write(unit=ids(j), fmt=*)  i, csv, pts, csv, tmp
          end do
        end if
      
        
        
        if (.not. allocated(fluxes)) then 
          flux = 0
          do jj=1, ubound(gauss_points%weight,1)
            quadpnt%order = jj
            call pde(proc)%flux(layer, quadpnt, vector_out=flux_tmp(1:D))
            flux = flux + flux_tmp*gauss_points%weight(jj)
          end do
          flux = flux / gauss_points%area
        else
          do iflux=1, ubound(fluxes,1)
            fluxes(iflux,:) = 0 
              do jj=1, ubound(gauss_points%weight,1)
                quadpnt%order = jj
                call pde(proc)%fluxes(iflux)%val(pde(proc),layer, quadpnt, vector_out=flux_tmp(1:D))
                fluxes(iflux,:) = fluxes(iflux,:)  + flux_tmp(1:D)
              end do
          end do
          fluxes = fluxes / gauss_points%area
        end if
        
        if (.not. pde(proc)%print_mass) then
          velocity_id = 3
        else
          velocity_id = 3 + ubound(pde(proc)%mass,1)
        end if
        
        if (.not. allocated(fluxes)) then
          write(unit=ids(velocity_id), fmt=*) i, csv,  pts, csv, flux(1:drutes_config%dimen)
        else
          do iflux=1, ubound(fluxes,1)
            write(unit=ids(velocity_id-1+iflux), fmt=*) i, csv, pts, csv, fluxes(iflux,:)
          end do
        end if
      end do

      do i=1, ubound(ids,1)
        if (ids(i) /= 0) close(ids(i))
      end do

    end subroutine print_pure
    
    subroutine print_gmsh(ids, proc, quadpnt)
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      integer, dimension(:), intent(in) :: ids
      integer(kind=ikind), intent(in) :: proc
      type(integpnt_str),  intent(in out) :: quadpnt
      integer(kind=ikind) :: i, j, layer, tag1, tag2
      logical, dimension(:), allocatable, save :: printed 
      real(kind=rkind) :: curtime, totflux
      
      
      if (.not. allocated(printed)) then
        allocate(printed(ubound(pde,1)))
        printed = .false.
      end if
          
      if (.not. quadpnt%globtime) then
        curtime = quadpnt%time4eval
      else
        curtime = time
      end if
      
      select case(drutes_config%dimen)
        case(2)
          tag1 = 2
          tag2 = 14
        case(3)
          tag1 = 4
          tag2 = 36
      end select
      
      if (.not. printed(proc)) then
       do i=1, ubound(ids,1)
            write(unit=ids(i), fmt="(a)") "$MeshFormat"
            write(unit=ids(i), fmt="(a)") "2.2 0 8"
            write(unit=ids(i), fmt="(a)") "$EndMeshFormat"
            write(unit=ids(i), fmt="(a)") "$Nodes"
            
            write(unit=ids(i), fmt=*) nodes%kolik
            
            do j=1, nodes%kolik
              if (drutes_config%dimen == 2) then
                write(unit=ids(i), fmt=*) j,  nodes%data(j,:), "1.0"
              else
                write(unit=ids(i), fmt=*) j,  nodes%data(j,:)
              end if
            end do
            write(unit=ids(i), fmt="(a)") "$EndNodes"
            write(unit=ids(i), fmt="(a)") "$Elements"
            write(unit=ids(i), fmt=*) elements%kolik
            
            do j=1, elements%kolik
              write(unit=ids(i), fmt=*) j, tag1, "2 ", elements%material(j), tag2, elements%data(j,:)
            end do
            write(unit=ids(i), fmt="(a)") "$EndElements"
!            write(unit=ids(i), fmt="(a)") "$NodeData"
            

        end do
        printed(proc) = .true.
      end if
      

      
      
      
      do i=1, ubound(ids,1)
        if (i/=2) then
          write(unit=ids(i), fmt="(a)") "$NodeData"     
          write(unit=ids(i), fmt=*) "1"
          
          !write solution
          if (i==1) then
            write(unit=ids(i), fmt=*) ' " ', "time =  ", time, " ",trim(pde(proc)%solution_name(1)) ,' " '
          
          ! write mass
          else if (i >= 3 .and. i <= ubound(pde(proc)%mass_name,1)+2) then

            write(unit=ids(i), fmt=*) ' " ', "time =  ", time, " ", trim(pde(proc)%mass_name(i-2,1)), ' " '
          ! write fluxes   
          else if ( i > ubound(pde(proc)%mass_name,1)+2 ) then
            write(unit=ids(i), fmt=*) ' " ',  "time =  ", time, " ", trim(pde(proc)%flux_name(1)) , ' " '
          end if
          
          write(unit=ids(i), fmt=*) "1"
          write(unit=ids(i), fmt=*) curtime  
          write(unit=ids(i), fmt=*) "3"
          write(unit=ids(i), fmt=*) postpro_run
          write(unit=ids(i), fmt=*) "1"
            

          
          write(unit=ids(i), fmt=*) nodes%kolik
          
           !write solution
          if (i==1) then
            do j=1, nodes%kolik
              quadpnt%order = j
              quadpnt%preproc = .true.
              write(unit=ids(i), fmt=*) j, pde(proc)%getval(quadpnt)
            end do
            
          ! write mass  
          else if (i >= 3 .and. i <= ubound(pde(proc)%mass_name,1)+2) then
            do j=1, nodes%kolik
              quadpnt%order = j
              quadpnt%preproc = .true.
              layer = elements%material(nodes%element(j)%data(1))
              write(unit=ids(i), fmt=*) j, pde(proc)%mass(i-2)%val(pde(proc), layer, quadpnt)
            end do
          
          
          ! write fluxes   
          else if ( i > ubound(pde(proc)%mass_name,1)+2) then
            do j=1, nodes%kolik
              quadpnt%order = j
              quadpnt%preproc = .true.
              layer = elements%material(nodes%element(j)%data(1))
              call pde(proc)%flux(layer, quadpnt, scalar=totflux)
              write(unit=ids(i), fmt=*) j, totflux
            end do
          end if
          write(unit=ids(i), fmt="(a)") "$EndNodeData"
        end if
      end do
    
end subroutine print_gmsh

  
    
end module postpro
