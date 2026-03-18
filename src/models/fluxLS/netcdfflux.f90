module netcdfflux

  contains

    subroutine ncflux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use globals 
      use global_objs
      use pde_objs
      use geom_tools
      use debug_tools
      use datetime
    
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
      
      real(kind=rkind), dimension(2) :: xy
      
      call getcoor(quadpnt, xy)
      
      
      
  
    end subroutine ncflux


end module netcdfflux
