
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

!> \file freeze_reader.f90
!! \brief Config reader for the snow and freezing soil.
!<

module freeze_read
  use typy
  use freeze_globs
  use debug_tools
  
  public :: freeze_reader, read_frrate
  
  contains
    subroutine read_frrate
      use typy
      use global_objs
      use pde_objs
      use globals
      use core_tools
      use readtools
      use freeze_globs
      
      integer :: i_err
      character(len=4096) :: msg
      
      call find_unit(file_freeze, 200)
      open(unit = file_freeze, file="drutes.conf/freeze/freeze.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
        print *, "missing drutes.conf/freeze/freeze.conf"
        ERROR STOP
      end if
    
      
      
      write(unit = msg, fmt = *) "Use freezing rate? yes - 1 or no -0"
      call fileread(frz_pnt, file_freeze, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(frz_pnt)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for freezing rate"
              print *, "the allowed options are:"
              print *, "                        1 = yes (freezing rate)"
              print *, "                        0 = no (instantaneous freezing)"
              call file_error(file_freeze)
      end select
      
      if(frz_pnt > 0) then
       fr_rate = .true.
          wat = 1
          heat_proc = 2
          ice = 3
          heat_solid = 4

      else
       fr_rate = .false.
         wat = 1
         heat_proc = 2
         heat_solid = 3
      end if
      
    close(file_freeze)	  
    
    end subroutine read_frrate

    
    !> opens and reads water.conf/matrix.conf, input data for the Richards equation in single mode, 
    !! Richards equation with the dual porosity regime - matrix domain
    subroutine freeze_reader(pde_loc)
      use typy
      use global_objs
      use pde_objs
      use globals
      use core_tools
      use readtools
      use freeze_globs
      
      class(pde_str), intent(in out) :: pde_loc
      integer :: ierr, i, j, filewww,i_err
      integer(kind=ikind) :: n, tmp_int
      character(len=1) :: yn
      character(len=4096) :: msg
      real(kind=rkind), dimension(:), allocatable :: tmpdata

      select case (drutes_config%name)
        case ("freeze")
          call find_unit(file_freeze, 200)
          open(unit = file_freeze, file="drutes.conf/freeze/freeze.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/freeze/freeze.conf"
            ERROR STOP
          end if
        case ("LTNE")
          call find_unit(file_freeze, 200)
          open(unit = file_freeze, file="drutes.conf/freeze/freeze.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/freeze/freeze.conf"
            ERROR STOP
          end if
        case default
          print *, "procedure called because of unexpected problem name"
          print *, "exited from freeze_read::freeze_reader"
          error stop
      end select

      write(unit = msg, fmt = *) "Use freezing rate? yes - 1 or no -0"
      call fileread(frz_pnt, file_freeze, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(frz_pnt)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for freezing rate"
              print *, "the allowed options are:"
              print *, "                        1 = yes (freezing rate)"
              print *, "                        0 = no (instantaneous freezing)"
              call file_error(file_freeze)
      end select
      
      if(frz_pnt > 0) then
       fr_rate = .true.
        write(unit = msg, fmt = *) "beta should be positive"
        call fileread(beta, file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg))    
        call fileread(beta_melt, file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg))     

      else
       fr_rate = .false.
      end if
      
        
      allocate(freeze_par(maxval(elements%material)))
      
      call fileread(n, file_freeze)
      
      backspace(file_freeze)
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/freeze/freeze.conf  &
        the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_freeze, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
        errmsg = trim(msg))
	
      write(unit = msg, fmt = *) "Your porous material should be snow or soil?"
	
      do i=1, ubound(freeze_par,1)
        call fileread(freeze_par(i)%material, file_freeze,&
        errmsg=trim(msg))
        select case(freeze_par(i)%material)
            case("Snow", "Soil")
              CONTINUE
            case default
              print *, "you have specified the wrong material keyword"
              print *, "the allowed options are:"
              print *, "                        Snow"
              print *, "                        Soil"
              call file_error(file_freeze)
          end select
      end do
      
     write(unit = msg, fmt = *) "Use freezing point depression? yes - 1 or no -0"
     call fileread(frz_pnt, file_freeze, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(frz_pnt)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for freezing point depression"
              print *, "the allowed options are:"
              print *, "                        1 = yes"
              print *, "                        0 = no"
              call file_error(file_freeze)
      end select
      
      if(frz_pnt > 0) then
       clap = .true.
      else
       clap = .false.
      end if
      
     write(unit = msg, fmt = *) "Add factor should be real"
     call fileread(fac_add, file_freeze, ranges=(/-huge(0.0_rkind), huge(0.0_rkind)/), errmsg=trim(msg))
     write(unit = msg, fmt = *) "Scaling factor should be real"
     call fileread(fac_scale, file_freeze, ranges=(/-huge(0.0_rkind), huge(0.0_rkind)/), errmsg=trim(msg))

     write(unit = msg, fmt = *) "Impedance factor should be positive"
     call fileread(Omega, file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg))

     write(unit = msg, fmt = *) "Use qlt? yes - 1 or no -0"
     call fileread(tmp_int, file_freeze, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
      select case(tmp_int)
            case(1_ikind,0_ikind)
              CONTINUE
            case default
              print *, "you have specified wrong input for qlt"
              print *, "the allowed options are:"
              print *, "                        1 = yes"
              print *, "                        0 = no"
              call file_error(file_freeze)
      end select
      
      if(tmp_int > 0) then
       qlt_log = .true.
      else
       qlt_log = .false.
      end if


      write(msg, *) "HINT 1 : check number of layers in matrix", new_line("a"), &
         "   HINT 2 : have you specified all values in the following order: ", new_line("a"), &
         "         alpha   n   m   theta_r   theta_s   "
      allocate(tmpdata(5))
      do i = 1, ubound(freeze_par,1)
        call fileread(tmpdata, errmsg=msg, fileid=file_freeze, checklen=.true.)
        freeze_par(i)%alpha = tmpdata(1)
        freeze_par(i)%n = tmpdata(2)
        freeze_par(i)%m = tmpdata(3)
        freeze_par(i)%thr = tmpdata(4)
        freeze_par(i)%ths = tmpdata(5)
      end do
      deallocate(tmpdata)

     

 
      do i=1, ubound(freeze_par,1)
        allocate(freeze_par(i)%Ks_local(drutes_config%dimen))
        allocate(freeze_par(i)%Ks(drutes_config%dimen, drutes_config%dimen))
        j = max(1,drutes_config%dimen-1)
        allocate(freeze_par(i)%anisoangle(j))
      end do
      
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
      
      do i = 1, ubound(freeze_par,1)
        call fileread(tmpdata, file_freeze, errmsg=msg, checklen=.TRUE.)
        
        if (drutes_config%dimen > 1) then
          freeze_par(i)%anisoangle(:) = tmpdata(1:drutes_config%dimen-1)
        else
          freeze_par(i)%anisoangle(:) = tmpdata(1)
        end if
        
        select case(drutes_config%dimen)
          case(1)
            freeze_par(i)%Ks_local(:) = tmpdata(2)
          case(2)
            freeze_par(i)%Ks_local(:) = tmpdata(2:3)
          case(3)
            freeze_par(i)%Ks_local(:) = tmpdata(3:5)
        end select

        call set_tensor(freeze_par(i)%Ks_local(:), freeze_par(i)%anisoangle(:),  freeze_par(i)%Ks)
      end do
      deallocate(tmpdata)


      do i=1, ubound(freeze_par,1)
        call comment(file_freeze)
        read(unit = file_freeze, fmt= *, iostat=ierr) freeze_par(i)%initcond, freeze_par(i)%icondtypeRE

          select case(freeze_par(i)%icondtypeRE)
            case("H_tot", "hpres", "theta","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        H_tot = total hydraulic head"
              print *, "                        hpres = pressure head"
              print *, "                        theta = water content"
              print *, "                        input = read from input file (drutes output file)"
              call file_error(file_freeze)
          end select
          if (ierr /= 0) then
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, "HINT: check number of line records of initial conditions in freeze/freeze.conf!"
            print *, "----------------------------------------"
            call file_error(file_freeze)
          end if
        end do
   
	
	
      call fileread(n, file_freeze, ranges=(/1_ikind, huge(1_ikind)/), &
      errmsg="at least one boundary must be specified (and no negative values here)")
      

      call readbcvals(unitW=file_freeze, struct=pde(wat)%bc, dimen=n, &
		      dirname="drutes.conf/freeze/")

      if(fr_rate) then
        do i=1, ubound(freeze_par,1)
            call comment(file_freeze)
            read(unit = file_freeze, fmt= *, iostat=ierr) freeze_par(i)%iceini, freeze_par(i)%icondtypeIce

            select case(freeze_par(i)%icondtypeIce)
                case("theta")
                CONTINUE
                case default
                print *, "you have specified wrong initial condition type keyword"
                print *, "the allowed options are:"
                print *, "                        theta = ice content"
                call file_error(file_freeze)
            end select
            if (ierr /= 0) then
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, "HINT: check number of line records of initial conditions in freeze/freeze.conf!"
                print *, "----------------------------------------"
                call file_error(file_freeze)
            end if
        end do
        
       call fileread(n, file_freeze, ranges=(/1_ikind, huge(1_ikind)/), &
       errmsg="at least one boundary must be specified (and no negative values here)")
      

       call readbcvals(unitW=file_freeze, struct=pde(ice)%bc, dimen=n, &
		      dirname="drutes.conf/freeze/")
		      
      end if
	
	
		      
      close(file_freeze)	      

      select case (drutes_config%name)
        case ("freeze")
          call find_unit(file_freeze, 200)
          open(unit = file_freeze, file="drutes.conf/freeze/freeze_heat_LTE.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/freeze/freeze_heat_LTE.conf"
            ERROR STOP
          end if
        case ("LTNE")
          call find_unit(file_freeze, 200)
          open(unit = file_freeze, file="drutes.conf/freeze/freeze_heat_LTNE.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/freeze/freeze_heat_LTNE.conf.conf"
            ERROR STOP
          end if
        case default
          print *, "procedure called because of unexpected problem name"
          print *, "exited from freeze_read::freeze_reader"
          error stop
      end select
      select case (drutes_config%name)
        case ("freeze")
          CONTINUE
        case ("LTNE")
          write(unit = msg, fmt = *) "HINT 1: Is the diameter positive?", new_line("a"),&
        "   HINT 2 : Did you define mean diameters for each layer?"
        
         do i=1, ubound(freeze_par,1)
            call fileread(freeze_par(i)%diameter, file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg))
         end do
      end select

      allocate(tmpdata(4))
     write(unit = msg, fmt = *) "HINT 1: Are all capacities positive?", new_line("a"),&
        "   HINT 2 : Did you define heat capacities for each layer?"

      do i=1, ubound(freeze_par,1)
        call fileread(r = tmpdata, fileid = file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), errmsg=trim(msg), checklen = .TRUE.)
        freeze_par(i)%Cs = tmpdata(1)
        freeze_par(i)%Cl = tmpdata(2)
        freeze_par(i)%Ci = tmpdata(3)
        freeze_par(i)%Ca = tmpdata(4)
      end do
      
      deallocate(tmpdata)
      
      select case (drutes_config%name)
        case ("freeze")
        
          write(unit = msg, fmt = *) "HINT 1: Is the snow density positive?", new_line("a"),&
          "   HINT 2 : Did you define snow density for each layer?"
  
          do i=1, ubound(freeze_par,1)
            call fileread(freeze_par(i)%snow_density, fileid = file_freeze, ranges = (/0.0_rkind,huge(0.0_rkind)/),&
            errmsg = trim(msg))
          end do
      
         write(unit = msg, fmt = *) "HINT 1: Did you define all 8 values for each layer?", new_line("a"), &
          "   HINT 2 : Did you define values for each layer?"
          allocate(tmpdata(8))
          do i=1, ubound(freeze_par,1)
            call fileread(r = tmpdata, fileid = file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
            errmsg=trim(msg), checklen = .TRUE.)
            freeze_par(i)%C1 = tmpdata(1)
            freeze_par(i)%C2 = tmpdata(2)
            freeze_par(i)%C3 = tmpdata(3)
            freeze_par(i)%C4 = tmpdata(4)
            freeze_par(i)%C5 = tmpdata(5)
            freeze_par(i)%F1 = tmpdata(6)
            freeze_par(i)%F2 = tmpdata(7)
            freeze_par(i)%beta = tmpdata(8)
          end do
          deallocate(tmpdata)
        case ("LTNE")
            allocate(tmpdata(4))
            write(unit = msg, fmt = *) "HINT 1: Are the thermal conductivities positive?", new_line("a"),&
          "   HINT 2 : Did you define thermal conductivties for each layer?"
        
           do i=1, ubound(freeze_par,1)
             call fileread(r = tmpdata, fileid = file_freeze, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
             errmsg=trim(msg), checklen = .TRUE.)
             freeze_par(i)%Ls = tmpdata(1)
             freeze_par(i)%Ll = tmpdata(2)
             freeze_par(i)%Li = tmpdata(3)
             freeze_par(i)%La = tmpdata(4)
          end do
      
          deallocate(tmpdata)
      
          write(unit = msg, fmt = *) "Does air change temperature? yes - 1 or no -0"
         call fileread(frz_pnt, file_freeze, ranges=(/0_ikind,1_ikind/), errmsg = trim(msg))
         select case(frz_pnt)
           case(1_ikind,0_ikind)
            CONTINUE
           case default
             print *, "you have specified wrong input for heat flow handling of air"
             print *, "the allowed options are:"
              print *, "                        1 = yes"
              print *, "                        0 = no"
              call file_error(file_freeze)
         end select
        if(frz_pnt > 0) then
         air = .true.
        else
         air = .false.
        end if
        
     end select

      call comment(file_freeze)
      read(unit = file_freeze, fmt= *, iostat=ierr) freeze_par(1)%icondtype

      select case(freeze_par(1_ikind)%icondtype)
        case("value","input")
          CONTINUE
        case default
          print *, "you have specified wrong initial condition type keyword"
          print *, "the allowed options are:"
          print *, "                        value = enter constant temp values"
          print *, "                        input = read from input file (drutes output file)"
          call file_error(file_freeze)
      end select
      
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      select case(freeze_par(1)%icondtype)
        case("value")
          do i=1, ubound(freeze_par,1)
            call fileread(r = freeze_par(i)%Tinit, fileid=file_freeze, errmsg=trim(msg), &
            ranges=(/-273.15_rkind, huge(0.0_rkind)/))
          end do
      end select
      
      write(unit=msg, fmt=*) "The number of boundaries should be greater than zero and smaller or equal the number of nodes"
      
      call fileread(n, file_freeze, ranges=(/1_ikind, nodes%kolik/),&
        errmsg=trim(msg))
      
      call readbcvals(unitW=file_freeze, struct=pde(heat_proc)%bc, dimen=n, &
          dirname="drutes.conf/freeze/")
      
      do i=lbound(pde(heat_proc)%bc,1), ubound(pde(heat_proc)%bc,1)
        select case(pde(heat_proc)%bc(i)%code)
          case(3)
            call fileread(hc, file_freeze)
        end select
      end do
      
      
    select case (drutes_config%name)
      case ("freeze")
        CONTINUE
      case ("LTNE")
        call comment(file_freeze)
        read(unit = file_freeze, fmt= *, iostat=ierr) freeze_par(1)%icondtypeTs

      select case(freeze_par(1_ikind)%icondtypeTs)
            case("value","input")
              CONTINUE
            case default
              print *, "you have specified wrong initial condition type keyword"
              print *, "the allowed options are:"
              print *, "                        value = enter constant temp values"
              print *, "                        input = read from input file (drutes output file)"
              call file_error(file_freeze)
      end select
      
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      select case(freeze_par(1)%icondtypeTs)
        case("value")
          do i=1, ubound(freeze_par,1)
            call fileread(r = freeze_par(i)%Tinit_s, fileid=file_freeze, errmsg=trim(msg), &
            ranges=(/-273.15_rkind, huge(0.0_rkind)/))
          end do
      end select
      
      write(unit=msg, fmt=*) "The number of boundaries should be greater than zero and smaller or equal the number of nodes"
      
      call fileread(n, file_freeze, ranges=(/1_ikind, nodes%kolik/),&
        errmsg=trim(msg))
      
      call readbcvals(unitW=file_freeze, struct=pde(heat_solid)%bc, dimen=n, &
          dirname="drutes.conf/freeze/") 
      
      do i=lbound(pde(heat_solid)%bc,1), ubound(pde(heat_solid)%bc,1)
        select case(pde(heat_solid)%bc(i)%code)
          case(3)
            call fileread(hc, file_freeze)
        end select
      end do  
        
    end select
    end subroutine freeze_reader


end module freeze_read
