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
      nullify(pde(1)%stabilize_element)
      if (LSbank_noflow .or. (LSsupg .and. LSsupg_factor>0) .or. (LSshock .and. LSshock_factor>0)) &
        pde(1)%stabilize_element => adenc_element_corrections
      

 
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

    subroutine adenc_element_corrections(pde_loc,el_id,dt,quadpnt_in)
      use typy
      use pde_objs
      use global_objs
      use ncglobvars
      use ncsupg, only: adenc_supg_element
      use ncboundary, only: adenc_bank_element
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: el_id
      real(kind=rkind), intent(in) :: dt
      type(integpnt_str), intent(in), optional :: quadpnt_in
      if ((LSsupg .and. LSsupg_factor>0) .or. (LSshock .and. LSshock_factor>0)) &
        call adenc_supg_element(pde_loc,el_id,dt,quadpnt_in)
      call adenc_bank_element(pde_loc,el_id,dt,quadpnt_in)
    end subroutine adenc_element_corrections
    
    


end module ncpointers
