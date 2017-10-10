module ADE_pointers
  public :: ADE
  public :: ADEkinsorb
  
  public :: ADE_processes
  
  contains
  
    subroutine ADE_processes(processes)
      use typy
      use readtools
      use ADE_globals
      use core_tools
      
      integer(kind=ikind) :: processes
      integer :: adeconf, ierr
      
      open(newunit=adeconf, file="drutes.conf/ADE/ADE.conf", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        call write_log("Unable to open drutes.conf/ADE/ADE.conf, exiting....")
        ERROR STOP
      end if

      
      call fileread(with_richards, adeconf)
      
      call fileread(no_solids, adeconf, ranges=(/1_ikind, huge(1_ikind)/))
      
      processes =  no_solids + 1
      
      if (with_richards) processes = processes + 1
      
    end subroutine ADE_processes
    
    
  
    subroutine ADE(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
      use ADE_globals
      use RE_pointers
      
      class(pde_str), intent(in out), dimension(:) :: pde_loc
      integer(kind=ikind) :: i, adepos
      
      if (with_richards) then
        adepos = 2
      else
        adepos = 1
      end if
      
      call ADE_read(pde_loc(adepos))
	    
      pde_loc(adepos)%pde_fnc(adepos)%dispersion => ADEdispersion
      
      pde_loc(adepos)%pde_fnc(adepos)%convection => ADE_convection

      pde_loc(adepos)%pde_fnc(adepos)%elasticity => ADE_tder_coef

      pde_loc(adepos)%mass => ADE_mass

      pde_loc(adepos)%pde_fnc(adepos)%reaction => ADE_reaction
            
      pde_loc(adepos)%pde_fnc(adepos)%zerord => ADE_zerorder
	  
      do i=lbound(pde_loc(adepos)%bc,1), ubound(pde_loc(adepos)%bc,1)
        select case(pde_loc(adepos)%bc(i)%code)
          case(1)
            pde_loc(adepos)%bc(i)%value_fnc => ADE_dirichlet
          case(2)
            pde_loc(adepos)%bc(i)%value_fnc => ADE_neumann
        end select
      end do    
	
      pde_loc(adepos)%flux => ADE_flux
      
      pde_loc(adepos)%initcond => ADE_icond
      
      if (with_richards) call REstdH(pde_loc(1))
      
      if (no_solids > 0) then 
        call ADEkinsorb(pde_loc(adepos+1:pde_common%processes))
      end if 
      
    
    end subroutine ADE
    
    subroutine ADEkinsorb(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
      
      class(pde_str), intent(in out), dimension(:) :: pde_loc  
      integer(kind=ikind) :: i
      
      call ADEcs_read(pde_loc)
      
!       pde_loc%pde_fnc(pde_loc%order)%elasticity => ADE_tder_cscs
!       
!       pde(pde_loc%order-1)%pde_fnc(pde_loc%order)%elasticity => ADE_tder_cscl
!       
!       pde_loc%pde_fnc(pde_loc%order-1)%reaction => ADE_cscl_react
!       
!       pde_loc%pde_fnc(pde_loc%order)%reaction => ADE_cscs_react
!       
!       allocate(pde_loc%bc(lbound(pde(pde_loc%order-1)%bc,1) : (ubound(pde(pde_loc%order-1)%bc,1) )  ))
!       
!       do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
!         pde_loc%bc(i)%code = 2
!         pde_loc%bc(i)%value_fnc => ADE_null_bc
!       end do 
!       
!       pde_loc%initcond => ADEcs_icond
      
    
    end subroutine ADEkinsorb
    
      
  

end module ADE_pointers
