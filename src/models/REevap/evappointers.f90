module evappointers
  public :: REevap_proc, REevap_linker
  contains
  
    subroutine REevap_proc(processes) 
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 2
    
    end subroutine REevap_proc
    
    subroutine REevap_linker()
      use typy
      use global_objs
      use pde_objs
      use globals
      use heat_pointers
      use re_pointers
      use evapglob
      use evap_heat_fnc
      use evap_RE_fnc
      
      
      call heat(pde(:))

      pde(re_ord)%pde_fnc(re_ord)%dispersion => evapdiffhh
      pde(re_ord)%pde_fnc(heat_ord)%dispersion => evapdiffhT
      pde(re_ord)%pde_fnc(re_ord)%zerord  =>  dtheta_vdt
          
    end subroutine REevap_linker
  

  


end module evappointers
