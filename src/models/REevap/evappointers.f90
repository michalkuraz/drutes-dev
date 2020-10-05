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

      
      


      call heat(pde(:))

      pde(re_ord)%pde_fnc(re_ord)%dispersion => REdiffhh
      
      pde(re_ord)%pde_fnc(heat_ord)%dispersion => REdiffhT
      
      pde(re_ord)%pde_fnc(re_ord)%elasticity  =>  REcapacityhh
      
      pde(re_ord)%pde_fnc(heat_ord)%elasticity  =>  REcapacityhT

!      pde(re_ord)%flux  =>  water_flux

!      pde(heat_ord)%pde_fnc(heat_ord)%dispersion => evapdiffTT
      
!      pde(heat_ord)%pde_fnc(re_ord)%dispersion => evapdiffTh
      
!      pde(heat_ord)%pde_fnc(heat_ord)%convection => evap_heatconvect
      
!      pde(heat_ord)%pde_fnc(heat_ord)%zerord  =>  Ldtheta_vdt
      
      
!      pde(re_ord)%bc(102)%value_fnc => evap4bc
      
!      pde(heat_ord)%bc(102)%value_fnc => soil_heat_flux
      
      
          
    end subroutine REevap_linker
    
    
    
    subroutine REevap_linker2()
    
    end subroutine REevap_linker2
  

  


end module evappointers
