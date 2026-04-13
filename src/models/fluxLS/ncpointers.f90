module ncpointers
  use typy
  
  contains
  
     subroutine nc_processes(processes)
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 1
      
    end subroutine nc_processes
    
    subroutine nclinker(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ncglobvars
      use init_netcdf

      class(pde_str), intent(in out) :: pde_loc
      
      call netcdf()
    
    end subroutine nclinker
    
    


end module ncpointers
