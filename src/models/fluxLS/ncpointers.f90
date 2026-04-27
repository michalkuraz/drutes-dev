module ncpointers
  use typy
  
  contains
  
     subroutine nc_processes(processes)
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 1
      
    end subroutine nc_processes
    
    subroutine nclinker()
      use typy
      use globals
      use global_objs
      use pde_objs
      use ncglobvars
      use init_netcdf
      use lsconstitutive


      integer(kind=ikind) :: i
      
      call netcdf()
      

 
      pde(1)%pde_fnc(1)%dispersion => ADElsdisp
      
      pde(1)%pde_fnc(1)%convection => ADEls_convection

      pde(1)%pde_fnc(1)%elasticity => ADEls_tder_coef
      
      pde(1)%mass(1)%val => ADEls_mass



	  
      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
        select case(pde(1)%bc(i)%code)
          case(1)
            pde(1)%bc(i)%value_fnc => ADEls_dirichlet
        end select
      end do    
	
      pde(1)%flux => ncflux
      
      pde(1)%initcond => ADEls_icond

      
     
    
    end subroutine nclinker
    
    


end module ncpointers
