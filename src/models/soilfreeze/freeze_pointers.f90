module freeze_pointers
  public :: freeze_processes, frz_pointers
  
  contains
  
    subroutine freeze_processes(processes)
      use typy
      use globals
      
      integer(kind=ikind), intent(out) :: processes
      
      
      select case (drutes_config%name)
        case ("freeze")
          processes = 2
          
        case ("LTNE")
          processes = 4
          
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::freeze_processes"
          error stop
      end select
        
    
    end subroutine freeze_processes
  
  
  
    subroutine frz_pointers()
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use freeze_globs
      use freeze_fnc
      
      call REstdH(pde(1))
      
      
    
    end subroutine frz_pointers
  
  



end module freeze_pointers
