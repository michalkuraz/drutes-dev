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
  public :: write_obs
  public :: get_RAM_use
  private :: print_scilab, print_pure, print_gmsh
  


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
      logical                                               :: anime
      integer(kind=ikind)                                   :: mode, no_prints
      character(len=256), dimension(:,:), allocatable       :: filenames
      integer(kind=ikind)                                   :: i, i_err, j, layer, proc, run
      character(len=64)                                     :: forma
      integer, dimension(:,:), pointer, save                :: ids
      integer, dimension(:,:), allocatable, target, save    :: ids_obs
      integer, dimension(:,:), allocatable, target, save    :: ids_anime
      character(len=5)                                      :: extension
      character(len=15)                                     :: prefix
      real(kind=rkind)                                      :: distance, avgval, val1, val2, val3, tmp, flux
      logical, save                                         :: first_run=.true.
      type(integpnt_str)                                    :: quadpnt
      integer(kind=ikind), save                             :: anime_run, anime_dec
      integer                                               :: ierr

      if (present(name)) then
        select case(name)
          case("obs")
            anime = .false.
            mode = 0
          case("avi")
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
        if (first_run) then
          anime_run = 0
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
  
      if (.not. allocated(ids_obs)) then
        allocate(ids(ubound(pde,1), 4))
        allocate(ids_obs(ubound(pde,1), 4))
        allocate(ids_anime(ubound(pde,1), 4))
      end if
      
      if (anime) then
        ids => ids_anime
      else
        ids => ids_obs
      end if
	

      allocate(filenames(ubound(pde,1),4))
 
      
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

        write(unit=filenames(proc,1), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", &
        trim(pde(proc)%solution_name(1)), "-",  run, trim(extension)


        write(unit=filenames(proc,2), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", & 
                    trim(pde(proc)%solution_name(1)), &
                    "-el_avg-",  run,  trim(extension)

        if (pde(proc)%print_mass) then
          write(unit=filenames(proc,3), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_",  &
            trim(pde(proc)%mass_name(1)), "-", run,  trim(extension)
        end if

        write(unit=filenames(proc,4), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", &
                     trim(pde(proc)%flux_name(1)), "-", &
                    run,  trim(extension)
      

        if ( (.not. anime .and. mode == 0)  .or. &
          (anime_run == 1 .and. anime) .or. & 
         ( .not. anime .and. (mode == -1 .and. postpro_run == 0 ) ) ) then

          do i=1, ubound(filenames,2)
            if (i /= 3 .or. pde(proc)%print_mass) then
              call find_unit(ids(proc, i), 6000)
              open(unit=ids(proc, i), file=trim(filenames(proc,i)), action="write", status="replace", iostat=ierr)
!               if (ierr /= 0) then
!                 call system("mkdir out/anime")
!                 open(unit=ids(proc, i), file=trim(filenames(proc,i)), action="write", status="replace", iostat=ierr)
!                 if (ierr /= 0) then
!                   print *, "unexpected system error, called from postpro::make_print()"
!                   error stop
!                 end if
!               end if
            end if
          end do  
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
          select case(observe_info%fmt)
            case("pure")
              call print_pure(ids(proc,:), proc, quadpnt)
            case("scil")
              call print_scilab(ids(proc,:), proc, quadpnt)
            case("gmsh")
              call print_gmsh(ids(proc,:), proc, quadpnt)
          end select
        end if
      end do


      if (.not. anime) then
        do proc=1, ubound(pde,1)
          do i=1,ubound(ids,2)
            if (i /= 3 .or. pde(proc)%print_mass) then
              if (mode == 0) then
                close(ids(proc,i))
              else
                call flush(ids(proc,i))
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
		
      integer(kind=ikind) :: i, layer, proc, D
      real(kind=rkind) :: val, massval
      real(kind=rkind), dimension(3) :: advectval
      type(integpnt_str) :: quadpnt

      
      quadpnt%type_pnt = "obpt"
      quadpnt%column=3
      D = drutes_config%dimen
      

      do proc=1, ubound(pde,1)
        do i =1, ubound(observation_array,1)
          layer = elements%material(observation_array(i)%element)

          quadpnt%order = i
          
          quadpnt%element = observation_array(i)%element
          
          quadpnt%preproc=.true.
          
          call pde(proc)%flux(layer=layer, quadpnt=quadpnt,  vector_out=advectval(1:D))

          val = pde(proc)%getval(quadpnt)

          massval = pde(proc)%mass(layer, quadpnt)
          
          observation_array(i)%cumflux(proc) = observation_array(i)%cumflux(proc) + &
              sqrt(dot_product(advectval(1:D), advectval(1:D)))*time_step
              
          write(unit=pde(proc)%obspt_unit(i), fmt="(50E24.12E3)") time, val, massval, advectval(1:D), &
                observation_array(i)%cumflux(proc)

          call flush(pde(proc)%obspt_unit(i))
          
        end do
      end do
    end subroutine write_obs
 

    !> this subroutine parses /proc/[PID]/status file in order to get RAM consumption statistics, it works Linux only
    subroutine get_RAM_use()
      use typy
      use globals
      use core_tools
      integer :: PID, fileid, i, ierr
      integer(kind=ikind) :: bytes
      character(len=2) :: byte_unit
      character(len=7) :: ch
      character(len=256) :: filename, format

      PID = getpid()

      call find_unit(fileid)


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
      
      open(unit=fileid, file=filename, action="read", status="old", iostat=ierr)

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

      read(unit=fileid, fmt=*) ch, bytes, byte_unit

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

      read(unit=fileid, fmt=*) ch, bytes, byte_unit

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
      
      integer, dimension(:), intent(in) :: ids
      integer(kind=ikind), intent(in) :: proc
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind) :: i, j, layer
      real(kind=rkind) :: tmp, totflux
      real(kind=rkind), dimension(3) :: flux
      real(kind=rkind), dimension(3,8) :: body
      real(kind=rkind), dimension(2) :: vct1, vct2
      real(kind=rkind), dimension(8) :: vct_tmp
      integer(kind=ikind) :: time_dec
      real(kind=rkind) :: curtime
      type(integpnt_str) :: qpntloc
      type(integpnt_str) :: quadpnt_loc

      
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
!         if (pde(proc)%print_mass) then
          qpntloc%element = i
          qpntloc%column = 2
          qpntloc%type_pnt = "gqnd"
          qpntloc%preproc=.true.
          tmp = 0
          do j=1, ubound(gauss_points%weight,1)
            qpntloc%order = j
            tmp = tmp + pde(proc)%mass(layer, qpntloc)*gauss_points%weight(j)
          end do
          
          body(:,5) = tmp/gauss_points%area
!         end if
          
        if (ubound(gauss_points%weight,1) > 1) then
          qpntloc%order = nint(ubound(gauss_points%weight,1)/2.0)
        else
          qpntloc%order = 1
        end if
        
        call pde(proc)%flux(layer, quadpnt_loc, vector_out=flux(1:drutes_config%dimen), scalar=totflux)

        body(:,6) = totflux
        
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
        
          ! 3 is for mass
          if (j /= 3 .or. pde(proc)%print_mass) then
            write(unit=ids(j), fmt=*) "x(", i, ",1) =", body(1,1), ";"
            write(unit=ids(j), fmt=*) "x(", i, ",2) =", body(2,1), ";"
            write(unit=ids(j), fmt=*) "x(", i, ",3) =", body(3,1), ";"
            
            write(unit=ids(j), fmt=*) "y(", i, ",1) =", body(1,2), ";"
            write(unit=ids(j), fmt=*) "y(", i, ",2) =", body(2,2), ";"
            write(unit=ids(j), fmt=*) "y(", i, ",3) =", body(3,2), ";"

            write(unit=ids(j), fmt=*) "z(", i, ",1) =", body(1,2+j), ";"
            write(unit=ids(j), fmt=*) "z(", i, ",2) =", body(2,2+j), ";"
            write(unit=ids(j), fmt=*) "z(", i, ",3) =", body(3,2+j), ";"
          end if  
            
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
      
        if (i /= 3 .or. pde(proc)%print_mass) then
          write(unit=ids(i), fmt=*) "f=gcf();"
          write(unit=ids(i), fmt=*) "clf(f,'reset');"
          write(unit=ids(i), fmt=*) "f.color_map=jetcolormap(256);"
          write(unit=ids(i), fmt=*) "colorbar(min(z),max(z));"
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
        end if

      end do 
    
    end subroutine print_scilab
    
    
    subroutine print_pure(ids, proc, quadpnt)
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      integer, dimension(:), intent(in) :: ids
      integer(kind=ikind), intent(in) :: proc
      type(integpnt_str),  intent(in out) :: quadpnt
      real(kind=rkind) :: curtime
      integer(kind=ikind) :: i, layer
      real(kind=rkind) ::  distance,  avgval
      type(integpnt_str) :: qpntloc
      real(kind=rkind), dimension(3) :: flux


    
      do i=1, nodes%kolik
        quadpnt%order = i
        quadpnt%preproc=.true.

        write(unit=ids(1), fmt=*) i,  nodes%data(i,:), pde(proc)%getval(quadpnt) 
        
        layer = elements%material(nodes%element(i)%data(1))
        ! 3 is for mass
        if (pde(proc)%print_mass) then
          write(unit=ids(3), fmt=*)  i, nodes%data(i,:), pde(proc)%mass(layer, quadpnt)
        end if
        
        call pde(proc)%flux(layer, quadpnt, vector_out=flux(1:drutes_config%dimen))

        write(unit=ids(4), fmt=*) i, nodes%data(i,:), flux
      end do

      close(ids(1))


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
      integer(kind=ikind) :: i
      logical, dimension(:), allocatable, save :: printed 
      real(kind=rkind) :: curtime
      
      
      if (.not. allocated(printed)) then
        allocate(printed(ubound(pde,1)))
        printed = .false.
      end if
          
      if (.not. quadpnt%globtime) then
        curtime = quadpnt%time4eval
      else
        curtime = time
      end if
      
      if (.not. printed(proc)) then
    
        write(unit=ids(1), fmt="(a)") "$MeshFormat"
        write(unit=ids(1), fmt="(a)") "2.2 0 8"
        write(unit=ids(1), fmt="(a)") "$EndMeshFormat"
        write(unit=ids(1), fmt="(a)") "$Nodes"
        
        write(unit=ids(1), fmt=*) nodes%kolik
        
        do i=1, nodes%kolik
          write(unit=ids(1), fmt=*) i,  nodes%data(i,:)
        end do
        write(unit=ids(1), fmt="(a)") "$EndNodes"
        write(unit=ids(1), fmt="(a)") "$Elements"
        write(unit=ids(1), fmt=*) elements%kolik
        do i=1, elements%kolik
          write(unit=ids(1), fmt=*) i,  elements%data(i,:)
        end do
        write(unit=ids(1), fmt="(a)") "$EndElements"
        write(unit=ids(1), fmt="(a)") "$NodeData"
        
        printed(proc) = .true.
      end if
      
      write(unit=ids(1), fmt="(a)") "$NodeData"     
      write(unit=ids(1), fmt=*) "1"
      write(unit=ids(1), fmt=*) ' " ', trim(pde(proc)%problem_name(2)), " ", trim(pde(proc)%solution_name(1)) , ' " '
      write(unit=ids(1), fmt=*) "1"
      write(unit=ids(1), fmt=*) curtime  !///tohle udává čas
      write(unit=ids(1), fmt=*) "3"
      write(unit=ids(1), fmt=*) postpro_run
      write(unit=ids(1), fmt=*) "1"
        

      
      write(unit=ids(1), fmt=*) nodes%kolik
      do i=1, nodes%kolik
        quadpnt%order = i
        write(unit=ids(1), fmt=*) i, pde(proc)%getval(quadpnt)
      end do
      write(unit=ids(1), fmt="(a)") "$EndNodeData"
    
    end subroutine print_gmsh

  
    
end module postpro
