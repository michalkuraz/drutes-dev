module bousspointers
  public :: boussi
  
  contains
  
    subroutine boussi(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use boussglob
      use boussfnc
      use boussread
	
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i
      
      call boussreader(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%dispersion => bouss_cond
      pde_loc%pde_fnc(pde_loc%order)%convection => bouss_adv
      pde_loc%pde_fnc(pde_loc%order)%elasticity => bouss_elast

      pde_loc%pde_fnc(pde_loc%order)%zerord => boussreact
	  
      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
	select case(pde_loc%bc(i)%code)
	  case(1,2)
	    pde_loc%bc(i)%value_fnc => bouss_bc
	end select
      end do
	      
	
      pde_loc%flux => darcy4bouss
      pde_loc%initcond => boussicond     
       
    end subroutine boussi


end module bousspointers