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
      integer(kind=ikind) :: i
      real(kind=rkind) :: r
      
      if (use_richards) then
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
      
      if (use_richards) call REstdH(pde_loc(1))

      if (use_sorption) then 
        call ADEkinsorb(pde_loc(adepos:no_solids+adepos))
      end if 
      
    
    end subroutine ADE
    
    subroutine ADEkinsorb(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
      use debug_tools
      
      class(pde_str), intent(in out), dimension(:) :: pde_loc  
      integer(kind=ikind) :: i, j
      
      call ADEcs_read(pde_loc)
      
      do i=2, ubound(pde_loc,1)
      
        pde_loc(i)%pde_fnc(pde_loc(i)%order)%elasticity => ADE_tder_cscs
      
        pde_loc(1)%pde_fnc(pde_loc(i)%order)%elasticity => ADE_tder_cscl
      
        pde_loc(i)%pde_fnc(adepos)%reaction => ADE_cscl_react
      
        pde_loc(i)%pde_fnc(pde_loc(i)%order)%reaction => ADE_cscs_react
      
        allocate(pde_loc(i)%bc(lbound(pde_loc(1)%bc,1) : (ubound(pde_loc(1)%bc,1) )  ))
        
        
        do j=lbound(pde_loc(i)%bc,1), ubound(pde_loc(i)%bc,1)
          pde_loc(i)%bc(j)%code = 0
          pde_loc(i)%bc(j)%value_fnc => ADE_null_bc
        end do 
      
        pde_loc(i)%initcond => ADEcs_icond
      end do

    
    end subroutine ADEkinsorb
    
      
  

end module ADE_pointers
