
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

!> \file LTNE_reader.f90
!! \brief Config reader for the snow and freezing soil.
!<

module LTNE_read
  use typy
  use LTNE_globs
  
  public :: LTNE_reader
  
  contains

    !> opens and reads water.conf/matrix.conf, input data for the Richards equation in single mode, 
    !! Richards equation with the dual porosity regime - matrix domain
    subroutine LTNE_reader(pde_loc)
      use typy
      use global_objs
      use pde_objs
      use globals
      use ltne_globs
      use core_tools
      use readtools
      
      class(pde_str), intent(in out) :: pde_loc
      integer :: ierr, i, j, filewww,i_err
      integer(kind=ikind) :: n, tmp_int
      character(len=1) :: yn
      character(len=4096) :: msg
      real(kind=rkind), dimension(:), allocatable :: tmpdata

      select case (drutes_config%name)
        case ("LTNE")
          call find_unit(file_LTNE, 200)
          open(unit = file_LTNE, file="drutes.conf/LTNE/LTNE.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/LTNE/LTNE.conf"
            ERROR STOP
          end if
        case default
          print *, "procedure called because of unexpected problem name"
          print *, "exited from LTNE_read::LTNE_reader"
          error stop
      end select
      

      write(msg, *) "define method of evaluation of constitutive functions for the Richards equation", new_line("a"), &
        "   0 - direct evaluation (not recommended, extremely resources consuming due to complicated exponential functions)", &
        new_line("a"), &
        "   1 - function values are precalculated in program initialization and values between are linearly approximated"

      call fileread(drutes_config%fnc_method, file_LTNE, ranges=(/0_ikind,1_ikind/),errmsg=msg)
      
      call fileread(maxpress, file_LTNE, ranges=(/-huge(0.0_rkind), huge(0.0_rkind)/), &
        errmsg="set a positive nonzero limit for maximum suction pressure (think in absolute values) ")
        maxpress = abs(maxpress)
      
      call fileread(drutes_config%fnc_discr_length, file_LTNE, ranges=(/tiny(0.0_rkind), maxpress/),  &
        errmsg="the discretization step for precalculating constitutive functions must be positive and smaller &
        than the bc")
        
        
      allocate(LTNE_par(maxval(elements%material)))
      
      call fileread(n, file_LTNE)
      
      backspace(file_LTNE)
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/LTNE/LTNE.conf  &
        the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_LTNE, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
        errmsg = trim(msg))
	
      write(unit = msg, fmt = *) "Your porous material should be snow or soil?"
	
      do i=1, ubound(LTNE_par,1)
        call fileread(LTNE_par(i)%material, file_LTNE,&
        errmsg=trim(msg))
        select case(LTNE_par(i)%material)
            case("Snow", "Soil")
              CONTINUE
            case default
              print *, "you have specified the wrong material keyword"
              print *, "the allowed options are:"
              print *, "                        Snow"
              print *, "                        Soil"
              call file_error(file_LTNE)
          end select
      end do
      
     write(unit = msg, fmt = *) "Use freezing point depression? yes - 1 or no -0"
     call fileread(frz_pnt, file_LTNE, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(frz_pnt)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for freezing point depression"
              print *, "the allowed options are:"
              print *, "                        1 = yes"
              print *, "                        0 = no"
              call file_error(file_LTNE)
      end select
      
      if(frz_pnt > 0) then
       clap = .true.
      else
       clap = .false.
      end if
      
      write(unit = msg, fmt = *) "HINT 1: Is the diameter positive?", new_line("a"),&
        "   HINT 2 : Did you define mean diameters for each layer?"
        
     do i=1, ubound(LTNE_par,1)
        call fileread(LTNE_par(i)%diameter, file_LTNE, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg))
     end do
      
      write(unit = msg, fmt = *) "HINT 1: Are the heat capacities positive?", new_line("a"),&
        "   HINT 2 : Did you define heat capacities for each layer?"
      
      allocate(tmpdata(4))
  
      do i=1, ubound(LTNE_par,1)
        call fileread(r = tmpdata, fileid = file_LTNE, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg), checklen = .TRUE.)
        LTNE_par(i)%Cs = tmpdata(1)
        LTNE_par(i)%Cl = tmpdata(2)
        LTNE_par(i)%Ci = tmpdata(3)
        LTNE_par(i)%Ca = tmpdata(4)
      end do
      
    write(unit = msg, fmt = *) "HINT 1: Are the thermal conductivities positive?", new_line("a"),&
        "   HINT 2 : Did you define thermal conductivties for each layer?"
        
     do i=1, ubound(LTNE_par,1)
        call fileread(r = tmpdata, fileid = file_LTNE, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg), checklen = .TRUE.)
        LTNE_par(i)%Ls = tmpdata(1)
        LTNE_par(i)%Ll = tmpdata(2)
        LTNE_par(i)%Li = tmpdata(3)
        LTNE_par(i)%La = tmpdata(4)
      end do
      
      deallocate(tmpdata)
      
     write(unit = msg, fmt = *) "Does air change temperature? yes - 1 or no -0"
     call fileread(frz_pnt, file_LTNE, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(frz_pnt)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for heat flow handling of air"
              print *, "the allowed options are:"
              print *, "                        1 = yes"
              print *, "                        0 = no"
              call file_error(file_LTNE)
      end select
      
      if(frz_pnt > 0) then
       air = .true.
      else
       air = .false.
      end if
      
      call comment(file_LTNE)
      read(unit = file_LTNE, fmt= *, iostat=ierr) LTNE_par(1)%icondtype

      select case(LTNE_par(1_ikind)%icondtype)
            case("value","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        value = enter constant temp values"
              print *, "                        input = read from input file (drutes output file)"
              call file_error(file_LTNE)
      end select
      
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      select case(LTNE_par(1)%icondtype)
        case("value")
          do i=1, ubound(LTNE_par,1)
            call fileread(r = LTNE_par(i)%Tinit_l, fileid=file_LTNE, errmsg=trim(msg), ranges=(/-273.15_rkind, huge(0.0_rkind)/))
          end do
      end select
      
      write(unit=msg, fmt=*) "The number of boundaries should be greater than zero and smaller or equal the number of nodes"
      
      call fileread(n, file_LTNE, ranges=(/1_ikind, nodes%kolik/),&
        errmsg=trim(msg))
      
      call readbcvals(unitW=file_LTNE, struct=pde(2)%bc, dimen=n, &
          dirname="drutes.conf/LTNE/")
      
      do i=lbound(pde(2)%bc,1), ubound(pde(2)%bc,1)
        select case(pde(2)%bc(i)%code)
          case(3)
            call fileread(hc, file_LTNE)
        end select
      end do
 
       call comment(file_LTNE)
      read(unit = file_LTNE, fmt= *, iostat=ierr) LTNE_par(1)%icondtypeTs

      select case(LTNE_par(1_ikind)%icondtypeTs)
            case("value","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        value = enter constant temp values"
              print *, "                        input = read from input file (drutes output file)"
              call file_error(file_LTNE)
      end select
      
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      select case(LTNE_par(1)%icondtypeTs)
        case("value")
          do i=1, ubound(LTNE_par,1)
            call fileread(r = LTNE_par(i)%Tinit_s, fileid=file_LTNE, errmsg=trim(msg), ranges=(/-273.15_rkind, huge(0.0_rkind)/))
          end do
      end select
      
      write(unit=msg, fmt=*) "The number of boundaries should be greater than zero and smaller or equal the number of nodes"
      
      call fileread(n, file_LTNE, ranges=(/1_ikind, nodes%kolik/),&
        errmsg=trim(msg))
      
      call readbcvals(unitW=file_LTNE, struct=pde(3)%bc, dimen=n, &
          dirname="drutes.conf/LTNE/")
      
      do i=lbound(pde(3)%bc,1), ubound(pde(3)%bc,1)
        select case(pde(3)%bc(i)%code)
          case(3)
            call fileread(hc, file_LTNE)
        end select
      end do
 
 ! Water flow
 
       
     write(unit = msg, fmt = *) "Use qlt? yes - 1 or no -0"
     call fileread(tmp_int, file_ltne, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(tmp_int)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for qlt"
              print *, "the allowed options are:"
              print *, "                        1 = yes"
              print *, "                        0 = no"
              call file_error(file_ltne)
      end select
      
      if(tmp_int > 0) then
       qlt_log = .true.
      else
       qlt_log = .false.
      end if
      
      
      do i=1, ubound(LTNE_par,1)
        allocate(LTNE_par(i)%Ks_local(drutes_config%dimen))
        allocate(LTNE_par(i)%Ks(drutes_config%dimen, drutes_config%dimen))
        j = max(1,drutes_config%dimen-1)
        allocate(LTNE_par(i)%anisoangle(j))
      end do
      

      write(msg, *) "HINT 1 : check number of layers in matrix", new_line("a"), &
         "   HINT 2 : have you specified all values in the following order: ", new_line("a"), &
         "         alpha   n   m   theta_r   theta_s   "
      allocate(tmpdata(5))
      do i = 1, ubound(LTNE_par,1)
        call fileread(tmpdata, errmsg=msg, fileid=file_LTNE, checklen=.true.)
        LTNE_par(i)%alpha = tmpdata(1)
        LTNE_par(i)%n = tmpdata(2)
        LTNE_par(i)%m = tmpdata(3)
        LTNE_par(i)%thr = tmpdata(4)
        LTNE_par(i)%ths = tmpdata(5)
      end do
      deallocate(tmpdata)

     


      write(msg, *) "HINT: check number of records of anisothropy description in water.conf/matrix.conf!!", &
        new_line("a") ,  &
        "      for 3D problem you must specify exactly Kxx, Kyy and Kzz values.", new_line("a"), &
        "       Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem", &
        new_line("a") , &
        "       for 2D problem supply only 1 angle, for 3D problem supply 2 angles, and", new_line("a"), &
        "       for 1D problem the angle value defines the angle between the VERTICAL and the flow trajectory", new_line("a"), &
        "       (carefull some other softwares consider HORIZONTAL!!)"
        
      
      select case(drutes_config%dimen)
        case(1,2)
          allocate(tmpdata(drutes_config%dimen+1))
        case(3)
          allocate(tmpdata(drutes_config%dimen+2))
      end select
      
      do i = 1, ubound(LTNE_par,1)
        call fileread(tmpdata, file_LTNE, errmsg=msg, checklen=.TRUE.)
        
        if (drutes_config%dimen > 1) then
          LTNE_par(i)%anisoangle(:) = tmpdata(1:drutes_config%dimen-1)
        else
          LTNE_par(i)%anisoangle(:) = tmpdata(1)
        end if
        
        select case(drutes_config%dimen)
          case(1)
            LTNE_par(i)%Ks_local(:) = tmpdata(2)
          case(2)
            LTNE_par(i)%Ks_local(:) = tmpdata(2:3)
          case(3)
            LTNE_par(i)%Ks_local(:) = tmpdata(3:5)
        end select

        call set_tensor(LTNE_par(i)%Ks_local(:), LTNE_par(i)%anisoangle(:),  LTNE_par(i)%Ks)
      end do


        do i=1, ubound(LTNE_par,1)
        call comment(file_LTNE)
        read(unit = file_LTNE, fmt= *, iostat=ierr) LTNE_par(i)%initcond, LTNE_par(i)%icondtypeRE

          select case(LTNE_par(i)%icondtypeRE)
            case("H_tot", "hpres", "theta","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        H_tot = total hydraulic head"
              print *, "                        hpres = pressure head"
              print *, "                        theta = water content"
              print *, "                        input = read from input file (drutes output file)"
              call file_error(file_LTNE)
          end select
          if (ierr /= 0) then
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, "HINT: check number of line records of initial conditions in LTNE/LTNE.conf!"
            print *, "----------------------------------------"
            call file_error(file_LTNE)
          end if
        end do
   
	
	
      call fileread(n, file_LTNE, ranges=(/1_ikind, huge(1_ikind)/), &
      errmsg="at least one boundary must be specified (and no negative values here)")
      

      call readbcvals(unitW=file_LTNE, struct=pde(1)%bc, dimen=n, &
		      dirname="drutes.conf/LTNE/")

		      
      close(file_LTNE)	      

    end subroutine LTNE_reader


end module LTNE_read
