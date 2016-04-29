module re_reader
  public :: res_read, red_read

  contains

    !> opens and reads water.conf/matrix.conf, input data for the Richards equation in single mode, 
    !! Richards equation with the dual porosity regime - matrix domain
    subroutine res_read(pde_loc)
      use typy
      use global_objs
      use pde_objs
      use globals
      use re_globals
      use core_tools
      use readtools
      
      class(pde_str), intent(in out) :: pde_loc
      integer :: ierr, i, j, filewww
      integer(kind=ikind) :: n
      character(len=1) :: yn
      character(len=4096) :: msg

      pde_loc%problem_name(1) = "RE_matrix"
      pde_loc%problem_name(2) = "Richards' equation"

      pde_loc%solution_name(1) = "press_head" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "h  [L]" !popisek grafu

      pde_loc%flux_name(1) = "flux"  
      pde_loc%flux_name(2) = "Darcian flow [L.T^{-1}]"

      pde_loc%mass_name(1) = "theta"
      pde_loc%mass_name(2) = "theta [-]"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !water.conf/matrix.conf
      call find_unit(file_waterm, 200)

      open(unit=file_waterm, file="drutes.conf/water.conf/matrix.conf", action="read", status="old", iostat = ierr)

      
      if (ierr /= 0) then
        print *, "missing drutes.conf/water.conf/matrix.conf file"
        ERROR STOP
      end if
      
      write(msg, *) "define method of evaluation of constitutive functions for the Richards equation", new_line("a"), &
	"   0 - direct evaluation (not recommended, extremely resources consuming due to complicated exponential functions)", &
	new_line("a"), &
	"   1 - function values are precalculated in program initialization and values between are linearly approximated"
      
      call fileread(drutes_config%fnc_method, file_waterm, ranges=(/0_ikind,1_ikind/),errmsg=msg)
      
      call fileread(maxpress, file_waterm, ranges=(/-huge(0.0_rkind), huge(0.0_rkind)/), &
	errmsg="set some positive nonzero limit for maximal suction pressure (think in absolute values) ")
	maxpress = abs(maxpress)
      
      call fileread(drutes_config%fnc_discr_length, file_waterm, ranges=(/tiny(0.0_rkind), maxpress/),  &
	errmsg="the discretization step for precalculating constitutive functions must be positive and smaller &
	then the bc")

      
      call fileread(n, file_waterm)
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/water.conf/matrix.conf  &
	the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
      backspace(file_waterm)
     
      call fileread(n, file_waterm, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
	errmsg=trim(msg))



 
      if (.not. allocated(vgmatrix)) then
	allocate (vgmatrix(n))
	do i=1, ubound(vgmatrix,1)
	  allocate(vgmatrix(i)%Ks_local(drutes_config%dimen))
	  allocate(vgmatrix(i)%Ks(drutes_config%dimen, drutes_config%dimen))
	  j = max(1,drutes_config%dimen-1)
	  allocate(vgmatrix(i)%anisoangle(j))
	end do
      end if


      do i = 1, ubound(vgmatrix,1)
        !nacteni pole pudnich charakteristik, podle schemat v deklaracich

        call comment(file_waterm)
        read(unit=file_waterm, fmt= *, iostat=ierr) vgmatrix(i)%alpha, vgmatrix(i)%n, vgmatrix(i)%m, &
                            vgmatrix(i)%Thr, vgmatrix(i)%Ths, vgmatrix(i)%Ss
        
        if (ierr /= 0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "HINT: check number of layers in matrix!"
          print *, "----------------------------------------"
          call file_error(file_waterm)
        end if
      end do




      do i = 1, ubound(vgmatrix,1)
          !nacteni pole pudnich charakteristik, podle schemat v deklaracich
          call comment(file_waterm)
          read(unit=file_waterm, fmt= *, iostat=ierr) vgmatrix(i)%anisoangle(:), vgmatrix(i)%Ks_local(:)
          if (ierr /= 0) then
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, "1st HINT: check number of records of anisothropy description in water.conf/matrix.conf!!"
            print *, "2nd HINT: for 3D problem you must specify exactly Kxx, Kyy and Kzz values."
            print *, "          Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem"
            print *, "3rd HINT: for 1D and 2D problem supply only 1 angle, for 3D problem supply 2 angles"
            print *, "----------------------------------------"
            call file_error(file_waterm)
          end if
        call set_tensor(vgmatrix(i)%Ks_local(:), vgmatrix(i)%anisoangle(:),  vgmatrix(i)%Ks)
      end do

      
      do i=1, ubound(vgmatrix,1)
	call fileread(vgmatrix(i)%sinkterm, file_waterm,  errmsg="Have you defined sink term for each layer?")
      end do
      
      if (.not. www) then
	do i=1, ubound(vgmatrix,1)
	  call comment(file_waterm)
	  read(unit=file_waterm, fmt= *, iostat=ierr) vgmatrix(i)%initcond, vgmatrix(i)%icondtype, &
							yn, vgmatrix(i)%rcza_set%val
	  select case(yn)
	    case("y")
	      vgmatrix(i)%rcza_set%use = .true.
	    case("n")
	      vgmatrix(i)%rcza_set%use = .false.
	    case default
	      write(msg, fmt=*) "type [y/n] value for using the retention curve zone approach at layer:", i
	      call file_error(file_waterm, msg)
	  end select
	  select case(vgmatrix(i)%icondtype)
	    case("H_tot", "hpres", "theta")
	      CONTINUE
	    case default
	      print *, "you have specified wrong initial condition type keyword"
	      print *, "the allowed options are:"
	      print *, "                        H_tot = total hydraulic head"
	      print *, "                        hpres = pressure head"
	      print *, "                        theta = water content"
	      call file_error(file_waterm)
	  end select
	  if (ierr /= 0) then
	    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	    print *, "HINT: check number of line records of initial conditions in water.conf/matrix.conf!"
	    print *, "----------------------------------------"
	    call file_error(file_waterm)
	  end if
	end do
      else
	do i=1, ubound(vgmatrix,1)
	  call comment(file_waterm)
	  read(unit=file_waterm, fmt= *, iostat=ierr) vgmatrix(i)%initcond, vgmatrix(i)%icondtype
          vgmatrix(i)%rcza_set%use = .false.
	  select case(vgmatrix(i)%icondtype)
	    case("H_tot", "hpres", "theta")
	      CONTINUE
	    case default
	      print *, "you have specified wrong initial condition type keyword"
	      print *, "the allowed options are:"
	      print *, "                        H_tot = total hydraulic head"
	      print *, "                        hpres = pressure head"
	      print *, "                        theta = water content"
	      call file_error(file_waterm)
	  end select
	  if (ierr /= 0) then
	    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	    print *, "HINT: check number of line records of initial conditions in water.conf/matrix.conf!"
	    print *, "----------------------------------------"
	    call file_error(file_waterm)
	  end if
	end do
      end if

   
	
	
      call fileread(n, file_waterm, ranges=(/1_ikind, huge(1_ikind)/), &
	errmsg="at least one boundary must be specified (and no negative values here)")
      

      call readbcvals(unitW=file_waterm, struct=pde_loc%bc, dimen=n, &
		      dirname="drutes.conf/water.conf/")
		      
       !!debuging
!        if (ubound(pde,1) == 3 .and. .not. allocated(pde(2)%bc)) then
! 	  allocate(pde(2)%bc(lbound(pde(1)%bc,1) : ubound(pde(1)%bc,1)))
! 	  allocate(pde(3)%bc(lbound(pde(1)%bc,1) : ubound(pde(1)%bc,1)))
! 	  pde(2)%bc = pde(1)%bc
! 	  pde(3)%bc = pde(1)%bc
! 	end if
	!!end debugging

		      
	close(file_waterm)	      

    end subroutine res_read


    subroutine red_read()
      use typy
      use global_objs
      use globals
      use core_tools

      integer :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !water.conf/fractures.conf
!       call find_unit(file_waterf, 200)
!       open(unit=file_waterf, file="drutes.conf/water.conf/fractures.conf", action="read", status="old", iostat = ierr)
!       if (ierr /= 0) then
!                 print *, "missing drutes.conf/water.conf/fractures.conf file"
!                 ERROR STOP
!       end if



!       call comment(file_waterf)
!       read(unit=file_waterf, fmt= *, iostat=ierr) retc_method_f ; if (ierr /= 0) call file_error(file_waterf)
! 
! 
!       call comment(file_waterf) 
!       read(unit=file_waterf, fmt = *,  iostat=ierr) n; if (ierr /= 0) call file_error(file_waterf)
!  
!       allocate (vgfractures(n))
!       do i=1, ubound(vgfractures,1)
!       allocate(vgfractures(i)%Ks_local(drutes_config%dimen))
!       allocate(vgfractures(i)%Ks(drutes_config%dimen, drutes_config%dimen))
!       j = max(1,drutes_config%dimen-1)
!       allocate(vgfractures(i)%anisoangle(j))
!       end do
! 
!       do i = 1, ubound(vgfractures,1)
!         !nacteni pole pudnich charakteristik, podle schemat v deklaracich
!         call comment(file_waterf)
!         read(unit=file_waterf, fmt= *, iostat=ierr) vgfractures(i)%alpha, vgfractures(i)%n, vgfractures(i)%m, vgfractures(i)%Thr, &
!                                                     vgfractures(i)%Ths, vgfractures(i)%Ss
!        
!         if (ierr /= 0) then
!           print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!           print *, "HINT: check number of layers in fractures!"
!           print *, "----------------------------------------"
!           call file_error(file_waterf)
!         end if
!         call set_tensor(vgfractures(i)%Ks_local(:), vgfractures(i)%anisoangle(:),  vgfractures(i)%Ks)
!       end do
!      
! 
!       do i = 1, ubound(vgfractures,1)
!         !nacteni pole pudnich charakteristik, podle schemat v deklaracich
!         call comment(file_waterf)
!         read(unit=file_waterf, fmt= *, iostat=ierr) vgfractures(i)%anisoangle, vgfractures(i)%Ks_local(:)
!         if (ierr /= 0) then
!           print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!           print *, "1st HINT: check number of records of anisothropy description in water.conf/fractures.conf!!"
!           print *, "2nd HINT: for 3D problem you must specify exactly Kxx, Kyy and Kzz values."
!           print *, "          Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem"
!           print *, "3rd HINT: for 1D and 2D problem supply only 1 angle, for 3D problem supply 2 angles"
!           print *, "  "
!           print *, "----------------------------------------"
!           call file_error(file_waterf)
!         end if
!       call set_tensor(vgfractures(i)%Ks_local(:), vgfractures(i)%anisoangle(:),  vgfractures(i)%Ks)
!       end do
! 
!       do i = 1, ubound(vgfractures,1)
!       call comment(file_waterf)
!       read(unit=file_waterf, fmt= *, iostat=ierr) vgfractures(i)%initcond, yn, vgfractures(i)%rcza_set%val
!       select case(yn)
!         case("y")
!           vgfractures(i)%rcza_set%use = .true.
!         case("n")
!           vgfractures(i)%rcza_set%use = .false.
!         case default
!           print *, "type [y/n] value for using the retention curve zone approach at layer:", i
!           call file_error(file_waterf)
!       end select
! 
!       if (ierr /= 0) then
!         print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!         print *, "HINT: check number of line records of initial conditions in water.conf/fractures.conf!"
!         print *, "----------------------------------------"
!         call file_error(file_waterf)
!       end if
!       end do
! 
!       allocate(dual_prop(ubound(vgfractures,1),2))
! 
!       do i = 1, ubound(dual_prop,1)
!       call comment(file_waterf)
!       read(unit=file_waterf, fmt= *, iostat=ierr) dual_prop(i,:)
!       if (ierr /= 0) then
!         print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!         print *, "HINT: check number of line records of dual permeability properties in matrix.conf/fractures.conf!"
!         print *, "----------------------------------------"
!         call file_error(file_waterf)
!       end if  
!       end do
! 
! 
!       call readbcvals(file_waterf, RICHARDS_dual%bc_sys2, dimen=(ubound(RICHARDS_dual%bc_sys1,1)-101),  &
!                         dirname="drutes.conf/water.conf/")


    end subroutine red_read

end module re_reader