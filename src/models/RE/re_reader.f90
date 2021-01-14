
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

!> \file re_reader.f90
!! \brief Config reader for the Richards equation.
!<


module re_reader
  public :: res_read
  private :: read_roots, REevapbc_read, read_shp_tab

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
      use debug_tools
      use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
      use, intrinsic :: iso_fortran_env, only: real32

      
      class(pde_str), intent(in out) :: pde_loc
      integer :: ierr,  filewww
      integer(kind=ikind) :: i, j
      integer(kind=ikind) :: n
      character(len=1) :: yn
      character(len=4096) :: msg
      real(kind=rkind), dimension(:), allocatable :: tmpdata

      real(kind=rkind) :: nan

      pde_loc%problem_name(1) = "RE_matrix"
      pde_loc%problem_name(2) = "Richards' equation"

      pde_loc%solution_name(1) = "press_head" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "h  [L]" !popisek grafu

      pde_loc%flux_name(1) = "flux"  
      pde_loc%flux_name(2) = "Darcian flow [L.T^{-1}]"
      
      deallocate(pde_loc%mass)
!      deallocate(pde_loc%mass_name)
      
      allocate(pde_loc%mass(2))
      allocate(pde_loc%mass_name(2,2))


      pde_loc%mass_name(1,1) = "theta"
      pde_loc%mass_name(1,2) = "theta [-]"
      
      pde_loc%mass_name(2,1) = "gwt"
      pde_loc%mass_name(1,2) = "gwt [y/n]"
      
      pde_loc%print_mass = .true.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !water.conf/matrix.conf

      open(newunit=file_waterm, file="drutes.conf/water.conf/matrix.conf", action="read", status="old", iostat = ierr)

      
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



 
      if (.not. allocated(vgset)) then
        allocate (vgset(n))
        do i=1, ubound(vgset,1)
          allocate(vgset(i)%Ks_local(drutes_config%dimen))
          allocate(vgset(i)%Ks(drutes_config%dimen, drutes_config%dimen))
          j = max(1,drutes_config%dimen-1)
          allocate(vgset(i)%anisoangle(j))
        end do
      end if

      write(msg, *) "HINT 1 : check number of layers in matrix", new_line("a"), &
         "   HINT 2 : have you specified all values in the following order: ", new_line("a"), &
         "         alpha   n   m   theta_r   theta_s   S_s "
      allocate(tmpdata(6))
      
      do i = 1, ubound(vgset,1)
      
        call comment(file_waterm)
        
        read(file_waterm, fmt = *, iostat=ierr) tmpdata
        
        if (ierr == 0) then 
          backspace file_waterm
          call fileread(tmpdata, errmsg=msg, fileid=file_waterm, checklen=.true., ranges=(/ 0.0_rkind, huge(0.0_rkind) /))
          vgset(i)%method = "vgfnc"
        else
          backspace file_waterm
          call fileread(vgset(i)%method, fileid=file_waterm, errmsg="set correct option for soil hydraulic functions", &
                      options=(/"vgfnc", "table"/))
          if (cut(vgset(i)%method) == "vgfnc") then
            backspace file_waterm
            tmpdata = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            read(unit=file_waterm, fmt = *) vgset(i)%method, tmpdata
            do j=1, ubound(tmpdata,1)
              if (isnan(tmpdata(j))) then
                write(msg, *) "not enough data for soil hydraulic properties for layer:", i
                call file_error(file_waterm, cut(msg))
              end if
            end do
          end if
        end if
        if (cut(vgset(i)%method) == "vgfnc") then
          vgset(i)%alpha=tmpdata(1)
          vgset(i)%n=tmpdata(2)
          vgset(i)%m=tmpdata(3)
          vgset(i)%thr=tmpdata(4)
          vgset(i)%ths=tmpdata(5)
          vgset(i)%Ss=tmpdata(6)
!          call check4tables()
        else
          call read_shp_tab(i)
        end if
      end do
                

      
     


      write(msg, *) "HINT: check number of records of anisothropy description in water.conf/matrix.conf!!", &
        new_line("a") ,  &
        "      for 3D problem you must specify exactly Kxx, Kyy and Kzz values.", new_line("a"), &
        "       Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem", &
        new_line("a") , &
        "       for 2D problem supply only 1 angle, for 3D problem supply 2 angles, and", new_line("a"), &
        "       for 1D problem the angle value defines the angle between the VERTICAL and the flow trajectory", new_line("a"), &
        "       (carefull some other softwares consider HORIZONTAL!!)"
        
      deallocate(tmpdata)
      
      select case(drutes_config%dimen)
        case(1,2)
          allocate(tmpdata(drutes_config%dimen+1))
        case(3)
          allocate(tmpdata(drutes_config%dimen+2))
      end select
      
      do i = 1, ubound(vgset,1)
        call fileread(tmpdata, file_waterm, errmsg=msg, checklen=.TRUE.)
        
        if (drutes_config%dimen > 1) then
          vgset(i)%anisoangle(:) = tmpdata(1:drutes_config%dimen-1)
        else
          vgset(i)%anisoangle(:) = tmpdata(1)
        end if
        
        select case(drutes_config%dimen)
          case(1)
            vgset(i)%Ks_local(:) = tmpdata(2)
          case(2)
            vgset(i)%Ks_local(:) = tmpdata(2:3)
          case(3)
            vgset(i)%Ks_local(:) = tmpdata(3:5)
        end select

        call set_tensor(vgset(i)%Ks_local(:), vgset(i)%anisoangle(:),  vgset(i)%Ks)
      end do

      allocate(roots(ubound(vgset,1)))
      
      call fileread(roots(1)%use, file_waterm)
      
      if (roots(1)%use) then
        call read_roots()
      end if
      
      if (.not. www) then
        do i=1, ubound(vgset,1)
          call comment(file_waterm)
          read(unit=file_waterm, fmt= *, iostat=ierr) vgset(i)%initcond, vgset(i)%icondtype, &
                    yn, vgset(i)%rcza_set%val
          select case(yn)
            case("y")
              vgset(i)%rcza_set%use = .true.
            case("n")
              vgset(i)%rcza_set%use = .false.
            case default
              write(msg, fmt=*) "type [y/n] value for using the retention curve zone approach at layer:", i
              call file_error(file_waterm, msg)
          end select
          select case(vgset(i)%icondtype)
            case("H_tot", "hpres", "theta","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        H_tot = total hydraulic head"
              print *, "                        hpres = pressure head"
              print *, "                        theta = water content"
              print *, "                        input = read from input file (drutes output file)"
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
        do i=1, ubound(vgset,1)
          call comment(file_waterm)
          read(unit=file_waterm, fmt= *, iostat=ierr) vgset(i)%initcond, vgset(i)%icondtype
                vgset(i)%rcza_set%use = .false.
          select case(vgset(i)%icondtype)
            case("H_tot", "hpres", "theta","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        H_tot = total hydraulic head"
              print *, "                        hpres = pressure head"
              print *, "                        theta = water content"
              print *, "                        input = read from input file (drutes output file)"
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

		      
      close(file_waterm)
      
!      !> Call reader if the user select atmospheric boundary
!       do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
!         if (pde_loc%bc(i)%code == 5) then 
!           call REevapbc_read()
!           exit
!         end if 
!       end do
       
       if (roots(1)%use) then
         call read_roots()
       end if
            

    end subroutine res_read
    
    
    subroutine read_roots()
      use typy
      use globals
      use readtools
      use re_globals
      
      integer :: fileid, ierr
      
      integer(kind=ikind) :: i
      real(kind=rkind), dimension(5) :: tmpdata
      character(len=4096) :: msg
      
      
      open(newunit=fileid, file="drutes.conf/water.conf/root4uptake.conf", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "file drutes.conf/water.conf/root4uptake.conf doesn't exist, exiting...."
        ERROR STOP
      end if
      
    
      do i=1, ubound(roots,1)
        call fileread(tmpdata, fileid, checklen=.true.)
        roots(i)%h_wet = tmpdata(1)
        roots(i)%h_start = tmpdata(2)
        roots(i)%h_end = tmpdata(3)
        roots(i)%h_dry = tmpdata(4)
        roots(i)%Smax = tmpdata(5)
        
        if (.not. roots(i)%h_wet >= roots(i)%h_start) then
          msg = "h_1 has to be greater or equal to h_2 -- see Feddes et al. (1978)"
          call file_error(fileid, msg)
        end if
        
        if (.not. roots(i)%h_start >= roots(i)%h_end) then
          msg = "h_2 has to be greater or equal to h_3 -- see Feddes et al. (1978)"
          call file_error(fileid, msg)
        end if
        
        if (.not. roots(i)%h_end >= roots(i)%h_dry) then
          msg = "h_3 has to be greater or equal to h_4 -- see Feddes et al. (1978)"
          call file_error(fileid, msg)
        end if
        
      end do
        
    
      close(fileid)

    end subroutine read_roots
    
    
    
    subroutine REevapbc_read()
      use typy
      use globals
      use readtools
      use re_globals
      
      integer :: fileid, ierr
      character(len=256), dimension(1) :: evapnames
      
      open(newunit=fileid, file="drutes.conf/water.conf/evap.conf", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "file drutes.conf/water.conf/evap.conf doesn't exist, exiting...."
        ERROR STOP
      end if
      
      evapnames(1) = "thorn"
      
      call fileread(evap_name, fileid, options=evapnames)
    
      
      close(fileid)

    
    end subroutine REevapbc_read
    
    
    subroutine read_shp_tab(layer)
      use typy
      use readtools
      use core_tools
      use re_globals
      use debug_tools
      
      integer(kind=ikind), intent(in) :: layer
      
      character(len=128)              :: filename
      integer(kind=ikind)             :: decimals=1
      character(len=64)               :: forma
      integer                         :: fileid, ierr
      real(kind=rkind), dimension(4)  :: tmpdata
      real(kind=rkind) :: tmpreal, step
      integer(kind=ikind) :: counter, i, tablesize, pos, j
      character(len=2048) :: msg
      real(kind=rkind), dimension(:,:), allocatable :: constable
    
      
      do 
        if (layer/(10**decimals) < 1) then
          EXIT
        else
          decimals = decimals + 1
        end if
      end do
      
      write(unit=forma, fmt="(a, I16, a)") "(a, I", decimals," a)"
      
      write(filename, fmt=forma) "drutes.conf/water.conf/SHP-", layer, ".in"
      
      open(newunit=fileid, file=cut(filename), status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open the file: ", cut(filename)
        print *, "check your configuration in: drutes.conf/water.conf where you assign soil hydraulic properties..."
        call file_error(file_waterm)
      end if
      
      call fileread(tmpdata(1:2), fileid, checklen=.true., ranges=(/0.0_rkind, 1.0_rkind /))
      
      if (tmpdata(1) >= tmpdata(2) ) then
        call file_error(fileid, "residual water content can't be greater or equal to saturated water content")
      else
        vgset(layer)%thr = tmpdata(1)
        vgset(layer)%ths = tmpdata(2)
      end if
      
      counter = 0
      do 
        call comment(fileid)
        read(fileid, fmt=*, iostat=ierr) tmpreal
        if (ierr == 0) then
          counter = counter + 1
        else
          EXIT
        end if
      end do
      
      close(fileid)
      open(newunit=fileid, file=cut(filename), status="old", action="read", iostat=ierr)
      
      call fileread(tmpdata(1:2), fileid, checklen=.true., ranges=(/0.0_rkind, 1.0_rkind /))
      
      !! constable(1) = h
      !! constable(2) = theta
      !! constable(3) = K
      !! constable(4) = C
      allocate(constable(counter,4))
      do i=1, counter
        call fileread(constable(i,:), fileid, checklen=.true.)
      end do
      
      if ( abs(constable(1,1)) > epsilon(tmpreal) ) then
        write(msg, fmt=*) "first value must be for h=0, in your data the first value is: ", constable(1,1)
        call file_error(fileid, cut(msg))
      end if
      
      if (abs(constable(1,2) - vgset(layer)%ths) > epsilon(tmpreal)) then
        write(msg, fmt = *) "Your saturated water content is:", vgset(layer)%ths, &
          "but in your data the water content for h=0 is: ", constable(1,2)
          
        call file_error(fileid, msg)
      end if 
        
      if (abs(constable(1,3) - 1.0_rkind) > epsilon(tmpreal)) then
        write(msg, fmt = *) "For h=0 the residual water content is 1.0, but in your data K_r(0)= ", constable(1,3)
        call file_error(fileid, msg)
      end if
      
      if (abs(constable(1,4)) > epsilon(tmpreal)) then
        write(msg, fmt = *) "For h=0 the rentention water capacity should be equal 0.0, but in your data C(0)=", constable(1,4)
        call file_error(fileid, msg)
      end if
       
      
      step = huge(step)
      
      tmpreal = 0
      
      do i=3, ubound(constable,1)
        tmpreal = tmpreal + abs(log10(-constable(i,1)) - log10(-constable(i-1,1)))
        if (abs(log10(-constable(i,1)) - log10(-constable(i-1,1))) < step) then
          step  = abs(log10(-constable(i,1)) - log10(-constable(i-1,1)))
        end if
      end do
      
      tablesize = int((log10(-constable(ubound(constable,1),1))-log10(-constable(2,1)))/step) + 1
      
      step = (log10(-constable(ubound(constable,1),1))-log10(-constable(2,1)))/tablesize
      
      allocate(vgset(layer)%logh(tablesize+1))
      allocate(vgset(layer)%Kr(0:tablesize+1))
      allocate(vgset(layer)%theta(0:tablesize+1))
      allocate(vgset(layer)%C(0:tablesize+1))
      allocate(vgset(layer)%dKdh(0:tablesize+1))
      
      vgset(layer)%logh(1) = log10(-constable(2,1))
      
      do i=2, tablesize+1
        vgset(layer)%logh(i) = vgset(layer)%logh(i-1) + step
      end do
      
      vgset(layer)%step4fnc = step
      
      vgset(layer)%theta(1) = constable(2,2)
      vgset(layer)%Kr(1) = constable(2,3)
      vgset(layer)%C(1) = constable(2,4)
      
      vgset(layer)%theta(0) = constable(1,2)
      vgset(layer)%Kr(0) = constable(1,3)
      vgset(layer)%C(0) = constable(1,4)

      pos = 1
      do i=2, ubound(vgset(layer)%logh,1)
        searchme: do j = pos, ubound(constable,1) - 1
          if (-10**(vgset(layer)%logh(i)) <= constable(j,1) .and. -10**(vgset(layer)%logh(i)) > constable(j+1,1)) then
            pos = j
            vgset(layer)%theta(i) = (constable(j+1,2)-constable(j,2))/(constable(j+1,1)- constable(j,1)) * &
                                     (-10**vgset(layer)%logh(i) - constable(j,1)) + constable(j,2)

            vgset(layer)%Kr(i) = (constable(j+1,3)-constable(j,3))/(constable(j+1,1)- constable(j,1)) * &
                                     (-10**vgset(layer)%logh(i) - constable(j,1)) + constable(j,3)
                                     
            vgset(layer)%C(i) = (constable(j+1,4)-constable(j,4))/(constable(j+1,1)- constable(j,1)) * &
                                     (-10**vgset(layer)%logh(i) - constable(j,1)) + constable(j,4)
            EXIT searchme
          end if
        end do searchme
      end do
      
      vgset(layer)%dKdh(0) = 0
      
      vgset(layer)%dKdh(1) = (vgset(layer)%Kr(0) - vgset(layer)%Kr(1))/10**vgset(layer)%logh(1)
        
      do i=2, ubound(vgset(layer)%dKdh,1)
        vgset(layer)%dKdh(i) = (vgset(layer)%Kr(i-1) - vgset(layer)%Kr(i))/(-10**vgset(layer)%logh(i-1) + 10**vgset(layer)%logh(i))
      end do
      
 
      
    
    end subroutine read_shp_tab
    
    
    subroutine check4tables()
      use typy
      use re_constitutive
      use global_objs
      use pde_objs
      
      integer(kind=ikind) :: i
      type(integpnt_str) :: quadpnt
      real(kind=rkind) :: hinit = -1e-3, logstep, Kr
      
      logstep = 0.05
      quadpnt%type_pnt="numb"
      do i=1, 200
        hinit = -10**(log10(-hinit) + logstep)
        quadpnt%this_is_the_value = hinit
        call  mualem(pde(1), 1_ikind, quadpnt, scalar=Kr)
        print *,  hinit, vangen(pde(1), 1_ikind, quadpnt), Kr, vangen_elast(pde(1), 1_ikind, quadpnt)
      end do
      
      stop
    
    end subroutine check4tables


   

end module re_reader
