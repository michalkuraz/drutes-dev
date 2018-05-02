module ADE_pointers
  use typy
  public :: ADE
  public :: ADEkinsorb
  
  public :: ADE_processes
  integer(kind=ikind), private :: adepos
  
  contains
  
    subroutine ADE_processes(processes)
      use typy
      use readtools
      use ADE_globals
      use core_tools
      use globals
      
      integer(kind=ikind), intent(out) :: processes
      integer :: adeconf, ierr
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      integer(kind=ikind) :: i
      character(len=4096) :: msg
      
      open(newunit=adeconf, file="drutes.conf/ADE/ADE.conf", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        call write_log("Unable to open drutes.conf/ADE/ADE.conf, exiting....")
        ERROR STOP
      end if

      call fileread(use_richards, adeconf)
      
      allocate(adepar(maxval(elements%material)))
      
      if (.not. use_richards) then
         allocate(tmp_array(2))
         do i=1, maxval(elements%material)
           call fileread(tmp_array, adeconf, errmsg="Convection has to be defined for each layer.", checklen=.true.)
           adepar(i)%convection = tmp_array(1)
           adepar(i)%water_cont = tmp_array(2)
         end do
      end if
      
      if (use_richards) then
        write(unit=msg, fmt=*) "HINT1: Have you commented out all lines with convection values? ", &
        "Since the convection is computed from the Richards equation."
      else
        write(unit=msg, fmt=*) "Is the number of lines with convection definition corresponding & 
          with the number of layers?"
      end if
      
      call fileread(use_sorption, adeconf, errmsg=msg)
      
      if (use_sorption) then
        call fileread(no_solids, adeconf, ranges=(/1_ikind, huge(1_ikind)/))
      else
        no_solids = 0
      end if
      
      processes =  no_solids + 1
      
      if (use_richards) processes = processes + 1
      
      
    end subroutine ADE_processes
    
    
  
    subroutine ADE()
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
      use ADE_globals
      use RE_pointers
      use re_constitutive
      
      integer(kind=ikind) :: i
      real(kind=rkind) :: r
      
      if (use_richards) then
        adepos = 2
      else
        adepos = 1
      end if
      
      
      call ADE_read(pde(adepos))
 
 
      pde(adepos)%pde_fnc(adepos)%dispersion => ADEdispersion
      
      pde(adepos)%pde_fnc(adepos)%convection => ADE_convection

      pde(adepos)%pde_fnc(adepos)%elasticity => ADE_tder_coef
      
      pde(adepos)%mass => ADE_mass

      pde(adepos)%pde_fnc(adepos)%reaction => ADE_reaction
            
            
      pde(adepos)%pde_fnc(adepos)%zerord => ADE_zerorder

	  
      do i=lbound(pde(adepos)%bc,1), ubound(pde(adepos)%bc,1)
        select case(pde(adepos)%bc(i)%code)
          case(0)
            pde(adepos)%bc(i)%value_fnc => re_null_bc
          case(1)
            pde(adepos)%bc(i)%value_fnc => ADE_dirichlet
          case(2)
            pde(adepos)%bc(i)%value_fnc => ADE_neumann
        end select
      end do    
	
      pde(adepos)%flux => ADE_flux
      
      pde(adepos)%initcond => ADE_icond

      
      if (use_richards) call REstdH(pde(1))


      if (use_sorption) then 
        call ADEkinsorb(adepos,no_solids+adepos)
      end if 
      
    
    end subroutine ADE
    
    subroutine ADEkinsorb(lb, tb)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
!       use debug_tools
      !>lb = low bound of pde, tp = top bound of pde to deal with
      integer(kind=ikind), intent(in) :: lb, tb
      integer(kind=ikind) :: i, j
      
      
      call ADEcs_read(lb, tb)
      
      do i=lb+1, tb
      
        pde(i)%pde_fnc(i)%elasticity => ADE_tder_cscs
      
        pde(lb)%pde_fnc(i)%elasticity => ADE_tder_cscl
      
        pde(i)%pde_fnc(lb)%reaction => ADE_cscl_react
      
        pde(i)%pde_fnc(i)%reaction => ADE_cscs_react
      
        allocate(pde(i)%bc(lbound(pde(lb)%bc,1) : (ubound(pde(lb)%bc,1) )  ))
        
        
        do j=lbound(pde(i)%bc,1), ubound(pde(i)%bc,1)
          pde(i)%bc(j)%code = 0
          pde(i)%bc(j)%value_fnc => ADE_null_bc
        end do 
      
        pde(i)%initcond => ADEcs_icond
      end do

    
    end subroutine ADEkinsorb
    
      
  

end module ADE_pointers
