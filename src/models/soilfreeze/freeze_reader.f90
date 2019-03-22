
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
  
  public :: freeze_read
  
  contains

    !> opens and reads water.conf/matrix.conf, input data for the Richards equation in single mode, 
    !! Richards equation with the dual porosity regime - matrix domain
    subroutine freeze_read(pde_loc)
      use typy
      use global_objs
      use pde_objs
      use globals
      use freeze_globals
      use core_tools
      use readtools
      
      class(pde_str), intent(in out) :: pde_loc
      integer :: ierr, i, j, filewww
      integer(kind=ikind) :: n
      character(len=1) :: yn
      character(len=4096) :: msg
      real(kind=rkind), dimension(:), allocatable :: tmpdata

      
      pde_loc%problem_name(1) = "heat"
      pde_loc%problem_name(2) = "Heat conduction equation with convection (Sophocleous, 1979)"

      pde_loc%solution_name(1) = "temperature" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "T " !popisek grafu

      pde_loc%flux_name(1) = "heat_flux"  
      pde_loc%flux_name(2) = "heat flux [W.L-2]"
      
      allocate(pde_loc%mass_name(0,2))

      call find_unit(file_heat, 200)
      open(unit = file_heat, file="drutes.conf/heat/heat.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
        print *, "missing drutes.conf/heat/heat.conf file"
        ERROR STOP
      end if
     
      call fileread(with_richards, file_heat)

      allocate(heatpar(maxval(elements%material)))
      
      call fileread(n, file_heat)
      
      backspace(file_heat)
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/heat/heat.conf  &
        the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_heat, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
        errmsg=trim(msg))
	
      write(unit=msg, fmt=*) "HINT 1: Is the heat capacity (matrix/matrix) positive?", new_line("a"), &
        "   HINT 2 : Is the number of heat capacity values corresponding to the amount of layers?"
	
      do i=1, ubound(heatpar,1)
        call fileread(heatpar(i)%C, file_heat, ranges=(/0.0_rkind,huge(0.0_rkind)/),&
        errmsg=trim(msg))
      end do
      
      write(unit=msg, fmt=*) "HINT 1: Is the heat capacity (matrix/water) positive?", new_line("a"), &
        "   HINT 2 : Is the number of heat capacity values corresponding to the amount of layers?"
        
      do i=1, ubound(heatpar,1)
        call fileread(heatpar(i)%C_w, file_heat, ranges=(/0.0_rkind,huge(0.0_rkind)/),&
        errmsg=trim(msg))
      end do
      
      
      write(unit=msg, fmt=*) "HINT 1: Is the heat conductivity positive?", new_line("a"), &
	"   HINT 2 : Is the number of heat conductivity values corresponding to the amount of layers?"
			      

 
      write(unit=msg, fmt=*) "HINT 1: Are all values anisotropy defining anisotropical diffusivity positive? ", &
       new_line("a"), new_line("a"),  &
        "   HINT 2: Have you defined an EXACT NUMBER of values for anisotropy?", new_line("a"), &
        "     (e.g. for 2D define angle and the maximal and minimal value of diffusivity, in total 3 values", new_line("a"),&
         "      if you switch from 2D to 1D don't forget to erase the last value, otherwise this error is generated)", &
        new_line("a"), new_line("a"), &
        "   HINT 3: The number of lines with heat conductivity has to correspond to the number of materials & 
        defined by your mesh", new_line("a"), new_line("a")
      
      
      allocate(tmp_array(drutes_config%dimen + 1))
      do i=1, ubound(heatpar,1)
        allocate(heatpar(i)%lambda_loc(drutes_config%dimen))
        call fileread(r=tmp_array, fileid=file_heat, ranges=(/0.0_rkind, huge(tmp)/), errmsg=trim(msg), checklen=.TRUE.)
        heatpar(i)%anisoangle = tmp_array(1)
        heatpar(i)%lambda_loc = tmp_array(2:drutes_config%dimen + 1)
        allocate(heatpar(i)%lambda(drutes_config%dimen, drutes_config%dimen))
        call set_tensor(heatpar(i)%lambda_loc, (/heatpar(i)%anisoangle/), heatpar(i)%lambda)
      end do
      
      
      write(unit=msg, fmt=*) "Did you specify convection vector component for each coordinate (e.g. x,y,z)"
      do i=1, ubound(heatpar,1)
        allocate(heatpar(i)%convection(drutes_config%dimen))
        call fileread(r=heatpar(i)%convection, fileid=file_heat, errmsg=trim(msg))
      end do
        
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      do i=1, ubound(heatpar,1)
       call fileread(r=heatpar(i)%Tinit, fileid=file_heat, errmsg=trim(msg))
      end do
       
		    
      
      write(unit=msg, fmt=*) "Hint: The number of lines for the heat source has to be equal to the number of materials."
      do i=1, ubound(heatpar,1)
       call fileread(r=heatpar(i)%source, fileid=file_heat, errmsg=trim(msg))
      end do
      
      
      write(unit=msg, fmt=*) "The number of boundaries should be greater than zero and smaller or equal the number of nodes"
      call fileread(n, file_heat, ranges=(/1_ikind, nodes%kolik/),&
        errmsg=trim(msg))
      
      call readbcvals(unitW=file_heat, struct=pde_loc%bc, dimen=n, &
          dirname="drutes.conf/heat/")
      
      pde_loc%problem_name(1) = "RE_matrix"
      pde_loc%problem_name(2) = "Richards' equation"

      pde_loc%solution_name(1) = "press_head" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "h  [L]" !popisek grafu

      pde_loc%flux_name(1) = "flux"  
      pde_loc%flux_name(2) = "Darcian flow [L.T^{-1}]"
      
      allocate(pde_loc%mass_name(1,2))

      pde_loc%mass_name(1,1) = "theta"
      pde_loc%mass_name(1,2) = "theta [-]"
      
      pde_loc%print_mass = .true.

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
        call fileread(tmpdata, errmsg=msg, fileid=file_waterm, checklen=.true.)
        vgset(i)%alpha=tmpdata(1)
        vgset(i)%n=tmpdata(2)
        vgset(i)%m=tmpdata(3)
        vgset(i)%thr=tmpdata(4)
        vgset(i)%ths=tmpdata(5)
        vgset(i)%Ss=tmpdata(6)
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

      
      do i=1, ubound(vgset,1)
        call fileread(vgset(i)%sinkterm, file_waterm,  errmsg="Have you defined sink term for each layer?")
      end do
      
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

    end subroutine res_read


end module freeze_read
