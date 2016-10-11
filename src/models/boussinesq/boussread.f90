module boussread
  public :: boussreader
  
  contains
    subroutine boussreader(pde_loc)
      use typy
      use readtools
      use core_tools
      use boussglob
      use global_objs
      use pde_objs
      
      class(pde_str), intent(in out) :: pde_loc
      integer, dimension(5) :: fileid
      integer :: ierr
      integer(kind=ikind) :: i
      real(kind=rkind), dimension(2) :: inputs
      
      
      pde_loc%problem_name(1) = "bouss"
      pde_loc%problem_name(2) = "Boussinesq equation"

      pde_loc%solution_name(1) = "water_table" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "h  [L]" !popisek grafu

      pde_loc%flux_name(1) = "flux"  
      pde_loc%flux_name(2) = "Darcian specific flux [L^2.T^{-1}]"

      pde_loc%mass_name(1) = "not_available"
      pde_loc%mass_name(2) = "n_a"
      
      
      call find_unit(fileid(1), 200)
      open(unit=fileid(1), file="drutes.conf/boussinesq.conf/slope.in", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
	print *, "missing drutes.conf/boussinesq.conf/slope.in file"
	error stop
      end if
      
      call find_unit(fileid(2), 200)
      open(unit=fileid(2), file="drutes.conf/boussinesq.conf/cond.in", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
	print *, "missing drutes.conf/boussinesq.conf/cond.in file"
	error stop
      end if
      
      call find_unit(fileid(3), 200)
      open(unit=fileid(3), file="drutes.conf/boussinesq.conf/rain.in", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
	print *, "missing drutes.conf/boussinesq.conf/rain.in file"
	error stop
      end if
      
      
      
      call find_unit(fileid(4), 200)
      open(unit=fileid(4), file="drutes.conf/boussinesq.conf/porosity.in", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
	print *, "missing drutes.conf/boussinesq.conf/porosity.in file"
	error stop
      end if
      
      call find_unit(fileid(5), 200)
      open(unit=fileid(5), file="drutes.conf/boussinesq.conf/bouss.boundaries", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
	print *, "missing drutes.conf/boussinesq.conf/bouss.boundaries file"
	error stop
      end if
      
      
      call bouss_slopes(1)%clear()
      call bouss_slopes(2)%clear()
      do 
	call comment(fileid(1))
	read(unit=fileid(1), fmt=*, iostat=ierr) inputs
	if (ierr == 0) then
	  call bouss_slopes(1)%fill(inputs(1))
	  call bouss_slopes(2)%fill(inputs(2))
	else
	  close(fileid(1))
	  EXIT
	end if
      end do
      
       do 
	call comment(fileid(2))
	read(unit=fileid(2), fmt=*, iostat=ierr) inputs
	if (ierr == 0) then
	  call bouss_K(1)%fill(inputs(1))
	  call bouss_K(2)%fill(inputs(2))
	else
	  close(fileid(2))
	  EXIT
	end if
      end do
      
      do 
	call comment(fileid(3))
	read(unit=fileid(3), fmt=*, iostat=ierr) inputs
	if (ierr == 0) then
	  call bouss_rain(1)%fill(inputs(1))
	  call bouss_rain(2)%fill(inputs(2))
	else
	  close(fileid(3))
	  EXIT
	end if
      end do


      call fileread(bouss_por, fileid(4))
      
      close(fileid(4))
      
      call fileread(pde_common%timeint_method, fileid(5))
      
      call fileread(bouss_icond, fileid(5))
      
      call readbcvals(unitW=fileid(5), struct=pde_loc%bc, dimen=2_ikind, &
		      dirname="drutes.conf/boussinesq.conf/")
      
      
      
    
    end subroutine boussreader

end module boussread