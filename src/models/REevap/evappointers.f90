module evappointers
  public :: REevap_proc, REevap_linker, REevap_linker2
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
      use evapglob
      use evap_RE_constitutive
      use evapreader
      use evap_heat_constitutive

      
      


      call heat(pde(:))
      
      call evapread()

      pde(re_ord)%pde_fnc(re_ord)%dispersion => REdiffhh
      
      pde(re_ord)%pde_fnc(heat_ord)%dispersion => REdiffhT
      
      pde(re_ord)%pde_fnc(re_ord)%elasticity  =>  REcapacityhh
      
      pde(re_ord)%pde_fnc(heat_ord)%elasticity  =>  REcapacityhT

      pde(re_ord)%flux  =>  totalflux
      
      pde(heat_ord)%pde_fnc(heat_ord)%elasticity => heatcap_TT

      pde(heat_ord)%pde_fnc(heat_ord)%dispersion => heat_cond
          
      pde(heat_ord)%pde_fnc(re_ord)%dispersion => heatdiffTh
      
      pde(heat_ord)%pde_fnc(heat_ord)%convection => convection4heat
      
      pde(heat_ord)%pde_fnc(heat_ord)%zerord  => heatsrc_w_roots 
      
      
!      pde(re_ord)%bc(102)%value_fnc => evap4bc
      
!      pde(heat_ord)%bc(102)%value_fnc => soil_heat_flux
      
      
          
    end subroutine REevap_linker
    
    
    
    subroutine REevap_linker2()
    
    end subroutine REevap_linker2
  

  


end module evappointers
