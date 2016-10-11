!> @file
!! Načtení proměnných ze konfiguračního souboru modwater.conf/matrix.conf

!> @brief  Modul obsahuje subrutnyny pro čtení ze souboru
!!
!! Modul obsahuje subrutnyny pro čtení ze souboru modwater.conf/matrix.conf
!! @author Jakub Jerabek
!! @version TEST
!!	


module modRE_reader
  public :: modre_read
  public :: modheat_read
  public :: modsolute_read
  
  private :: mod_init_bc
  private :: mod_read
  
 contains

 !> @brief Subroutina načítá ze souboru modwater.conf/matrix.conf.
 !!
 !! @param[out] hyd_prop proměnna s datovým typem pro globální charakterystiky prostředí
 !! @param[out] saturation proměnná pro objevoou vlhkost a hydraulickou vodivost
 
  subroutine modre_read()
    use typy
    use modRE_globals
    use pde_objs
    use core_tools
    use readtools
    use debug_tools
    
    integer :: ierr
    
    pde(1)%problem_name(1) = "RE_matrix_Noborio"
    pde(1)%problem_name(2) = "Richards' equation (Noborio)"

    pde(1)%solution_name(1) = "volumetric_water_content" !nazev vystupnich souboru
    pde(1)%solution_name(2) = "\theta  [L^{3}/L^{3}]" !popisek grafu

    pde(1)%flux_name(1) = "flux"  
    pde(1)%flux_name(2) = "Noborio fluxes [L.T^{-1}]"

    pde(1)%mass_name(1) = "theta"
    pde(1)%mass_name(2) = "\theta [L^{3}/L^{3}]"
    
    call mod_read()
    
    call find_unit(file_reinitbc, 200)
    open(unit=file_reinitbc, file="drutes.conf/modwater.conf/init_bc_RE_conditions.conf", &
	status="old", action="read", iostat=ierr) 
    call mod_init_bc(file_reinitbc, 1)
    
  
  end subroutine modre_read
  
  
  subroutine modheat_read
    use typy
    use modRE_globals
    use core_tools
    use readtools
    use pde_objs
    
    integer ::  ierr
    
    pde(2)%problem_name(1) = "heat_transport"
    pde(2)%problem_name(2) = "Heat transport equation"

    pde(2)%solution_name(1) = "temperature" !nazev vystupnich souboru
    pde(2)%solution_name(2) = "T  [K]" !popisek grafu

    pde(2)%flux_name(1) = "flux"  
    pde(2)%flux_name(2) = "Noborio fluxes [L.T^{-1}]"

    pde(2)%mass_name(1) = "dummy"
    pde(2)%mass_name(2) = "[-]"
    
    call mod_read()
    
    call find_unit(file_heatinitbc, 200)
    open(unit=file_heatinitbc, file="drutes.conf/modwater.conf/init_bc_Heat_conditions.conf", &
	status="old", action="read", iostat=ierr) 
   call mod_init_bc(file_heatinitbc, 2)
  
  
  end subroutine modheat_read
  
  subroutine modsolute_read
    use typy
    use modRE_globals
    use core_tools
    use readtools
    use pde_objs
    
    integer ::  ierr
    
    pde(3)%problem_name(1) = "Solute_transport"
    pde(3)%problem_name(2) = "Solute transport"

    pde(3)%solution_name(1) = "concentration" !nazev vystupnich souboru
    pde(3)%solution_name(2) = "C  [L^{3}/L^{3}]" !popisek grafu

    pde(3)%flux_name(1) = "flux"  
    pde(3)%flux_name(2) = "Noborio fluxes [L.T^{-1}]"

    pde(3)%mass_name(1) = "dummy"
    pde(3)%mass_name(2) = "[-]"    
    
    call mod_read()
    
    call find_unit(file_soluteinitbc, 200)
    open(unit=file_soluteinitbc, file="drutes.conf/modwater.conf/init_bc_Solute_conditions.conf", &
	status="old", action="read", iostat=ierr) 

    call mod_init_bc(file_soluteinitbc, 3)
 
  end subroutine modsolute_read
  
  subroutine mod_init_bc(file_id, process)
    use modRE_globals
    use pde_objs
    use core_tools
    use readtools
    
    integer, intent(in) :: file_id, process
    integer :: i, ierr
    integer(kind=ikind) :: n
    character(len=1) :: yn
    
    
    do i=1, ubound(hyd_prop,1)
	  call comment(file_id)
	  read(unit=file_id, fmt= *, iostat=ierr) hyd_prop(i)%initcond(process), yn, &
							hyd_prop(i)%rcza_set(process)%val
	  select case(yn)
	    case("y")
	      hyd_prop(i)%rcza_set(process)%use = .true.
	    case("n")
	      hyd_prop(i)%rcza_set(process)%use = .false.
	    case default
	      print *, "type [y/n] value for using the retention curve zone approach at layer:", i
	      call file_error(file_id)
	  end select
	  if (ierr /= 0) then
	    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	    print *, "HINT: check number of line records of initial conditions in water.conf/matrix.conf!"
	    print *, "----------------------------------------"
	    call file_error(file_id)
	  end if
      end do
      
      call fileread(n, file_id)        

      call readbcvals(unitW=file_id, struct=pde(process)%bc, dimen=n, &
		      dirname="drutes.conf/modwater.conf/")


  end subroutine mod_init_bc
  
  subroutine mod_read()
    use typy
    use pde_objs
    use modRE_globals
    use RE_globals
    use readtools
    use core_tools
    use debug_tools 


    integer (kind=ikind) :: i, j, n
    integer :: ierr 
    real(kind=rkind), dimension(3) :: wrkreal

    
    if (.not. allocated(hyd_prop)) then
    call find_unit(file_modwaterm, 200) 
    
    open(unit=file_modwaterm, file="drutes.conf/modwater.conf/modmatrix.conf", &
	status="old", action="read", iostat=ierr) 

      if (ierr /= 0) then
        print *, "missing drutes.conf/water.conf/matrix.conf file"
        ERROR STOP
      end if
    call fileread(pde_common%timeint_method, file_modwaterm)
    
    call fileread(n, file_modwaterm)

    
    allocate (hyd_prop(1:n))
    
    
    do i =1, n
      call fileread(wrkreal(1:2), file_modwaterm)
      hyd_prop(i)%theta_s = wrkreal(1)
      hyd_prop(i)%theta_r = wrkreal(2)
    end do
    
    do i =1, n
      call fileread(wrkreal(1:2), file_modwaterm)
      hyd_prop(i)%hc_a 	= wrkreal(1)
      hyd_prop(i)%hc_b	= wrkreal(2)
    end do

    
    do i=1, ubound(hyd_prop,1)
      allocate(hyd_prop(i)%Ks_local(drutes_config%dimen))
      allocate(hyd_prop(i)%Ks(drutes_config%dimen, drutes_config%dimen))
      j = max(1,drutes_config%dimen-1)
      allocate(hyd_prop(i)%anisoangle(j))
    end do
    
    do i = 1, n
      !nacteni pole pudnich charakteristik, podle schemat v deklaracich
      call comment(file_modwaterm)
      read(unit=file_modwaterm, fmt= *, iostat=ierr) hyd_prop(i)%anisoangle(:), hyd_prop(i)%Ks_local(:)
      if (ierr /= 0) then
	print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	print *, "1st HINT: check number of records of anisothropy description in modwater.conf/modmatrix.conf!!"
	print *, "2nd HINT: for 3D problem you must specify exactly Kxx, Kyy and Kzz values."
	print *, "          Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem"
	print *, "3rd HINT: for 1D and 2D problem supply only 1 angle, for 3D problem supply 2 angles"
	print *, "----------------------------------------"
	call file_error(file_modwaterm)
      end if
      call set_tensor(hyd_prop(i)%Ks_local(:), hyd_prop(i)%anisoangle(:),  hyd_prop(i)%Ks)
    end do

    do i =1, n
      call fileread(wrkreal(1:2), file_modwaterm)
      hyd_prop(i)%specific_storage	= wrkreal(1)
      hyd_prop(i)%porosity 		= wrkreal(2)
    end do
    
    do i =1, n
      call fileread(hyd_prop(i)%heat_capacity_soil, file_modwaterm)
    end do
    
    do i=1, n
      allocate(hyd_prop(i)%apparent_T_cap_soil_loc(drutes_config%dimen))
      allocate(hyd_prop(i)%apparent_T_cap_soil(drutes_config%dimen, drutes_config%dimen))
    end do
    
    do i = 1, n
      call comment(file_modwaterm)
      read(unit=file_modwaterm, fmt= *, iostat=ierr) hyd_prop(i)%apparent_T_cap_soil_loc(:)
      if (ierr /= 0) then
	print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	print *, "1st HINT: check number of records of anisothropy description in modwater.conf/modmatrix.conf!!"
	print *, "2nd HINT: for 3D problem you must specify exactly Kxx, Kyy and Kzz values."
	print *, "          Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem"
	print *, "3rd HINT: for 1D and 2D problem supply only 1 angle, for 3D problem supply 2 angles"
	print *, "----------------------------------------"
	call file_error(file_modwaterm)
      end if
      call set_tensor(hyd_prop(i)%apparent_T_cap_soil_loc(:), hyd_prop(i)%anisoangle(:),  hyd_prop(i)%apparent_T_cap_soil)
    end do
    
    
    
    
    do i =1, n
      call fileread(hyd_prop(i)%specific_heat_soil_water, file_modwaterm)
    end do
    
    do i =1, n
      call fileread(hyd_prop(i)%kappa, file_modwaterm)
    end do
    
    do i =1, n
      call fileread(hyd_prop(i)%hdisp_select, file_modwaterm)
    end do
    
    do i=1, n
      allocate(hyd_prop(i)%hdisp_loc(drutes_config%dimen))
      allocate(hyd_prop(i)%hdisp(drutes_config%dimen, drutes_config%dimen))
    end do
    
    do i =1, n
      select case(hyd_prop(i)%hdisp_select)
      case (0)
        call comment(file_modwaterm)
	read(unit=file_modwaterm, fmt= *, iostat=ierr) hyd_prop(i)%hdisp_loc(:)
	if (ierr /= 0) then
	  print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	  print *, "1st HINT: check number of records of anisothropy description in modwater.conf/modmatrix.conf!!"
	  print *, "2nd HINT: for 3D problem you must specify exactly Kxx, Kyy and Kzz values."
	  print *, "          Obviously for 2D problem there must be exactly only Kxx and Kzz value, analogicaly for 1D problem"
	  print *, "3rd HINT: for 1D and 2D problem supply only 1 angle, for 3D problem supply 2 angles"
	  print *, "----------------------------------------"
	  call file_error(file_modwaterm)
	end if
	call set_tensor(hyd_prop(i)%hdisp_loc(:), hyd_prop(i)%anisoangle(:),  hyd_prop(i)%hdisp)
      case (1)
	continue
      case default
	print*, "Warning: incorect assigment in a type od hydrodynamic"
	print*, "dispersion coefficient. Section [2.3] in modwater.conf/modmatrix.conf"
	error stop
      end select
    end do
    
    do i =1, n
      call fileread(hyd_prop(i)%G_STVF_select, file_modwaterm)
    end do
    
    do i=1, n
      select case(hyd_prop(i)%G_STVF_select)
      case (0)
	print*, "Nedodelany STVF model"
	call wait("Presto pokracovat?")
      case (1)
	print*, "Nedodelany Gain faktor model"
	call wait("Presto pokracovat?")
      case (2)
	print*, "Nedodelany Gain faktor model jako linearni funkci"
	call wait("Presto pokracovat?")
      case (3)
	call fileread(hyd_prop(i)%gain_factor, file_modwaterm)
	call fileread(hyd_prop(i)%gain_factor_T, file_modwaterm)
      case default
	print*, "Warning: incorect assigment in a Gain factor model / STVF model."
	print*, "         Section [3] in modwater.conf/modmatrix.conf."
	error stop
      end select
    end do
    
    
    do i =1, n
      call fileread(wrkreal(1:3), file_modwaterm)
      hyd_prop(i)%viscosity_0		= wrkreal(1)
      hyd_prop(i)%molecular_weight	= wrkreal(2)
      hyd_prop(i)%molecular_weight_NaCl	= wrkreal(3)
    end do
    
    do i =1, n
      call fileread(wrkreal(1:3), file_modwaterm)
      hyd_prop(i)%gas_constant		= wrkreal(1)
      hyd_prop(i)%hydrated_iont_r	= wrkreal(2)
      hyd_prop(i)%water_molecul_r	= wrkreal(3)
    end do
    
    do i =1, n
      call fileread(wrkreal(1:3), file_modwaterm)
      hyd_prop(i)%half_of_solute_film	= wrkreal(1)
      hyd_prop(i)%clay_fraction		= wrkreal(2)
      hyd_prop(i)%density_soil_water	= wrkreal(3)
    end do

    close(unit=file_modwaterm)
    else 
      continue
    end if
    
  end subroutine mod_read
 
end module modRE_reader