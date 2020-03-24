module evap_heat_fnc

  public :: evapdiffTT !, diffTh, convectT, heat_flux, sourceT
  
  contains
    
    subroutine evapdiffTT(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs
      use evapglob
      use heat_fnc
      use evap_RE_fnc
      use heat_globals
      use evapextras
      use debug_tools
      use re_constitutive

      class(pde_str), intent(in) :: pde_loc
       !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      integer(kind=ikind) :: D, i
      real(kind=rkind) :: L
      real(kind=rkind) :: Kvt
      
      D = drutes_config%dimen
      
      L = latent_heat_wat(quadpnt)

      if (present(tensor)) then
        tensor =  heatpar(layer)%lambda
      end if
      
      Kvt = hydraulic_vT(layer, quadpnt)

      do i=1, D
        tensor(i,i) = tensor(i,i) + Kvt*L
      end do
    end subroutine evapdiffTT
    
    
    subroutine evapdiffTh(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs
      use evapglob
      use heat_fnc
      use evap_RE_fnc
      use heat_globals
      use evapextras
      use debug_tools
      use re_constitutive

      class(pde_str), intent(in) :: pde_loc
       !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar
      
      real(kind=rkind) :: L, Kv
      integer(kind=ikind) :: i


      L = latent_heat_wat(quadpnt)
      
      Kv = hydraulic_vh(layer, quadpnt)
      
      tensor = L*Kv
    
      do i=1, drutes_config%dimen
        tensor(i,i) = L*Kv
      end do
      
    end subroutine evapdiffTh
    
    
    subroutine evap_convect(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
      use evap_RE_fnc
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> input vector
      real(kind=rkind), dimension(:), intent(in), optional  :: vector_in
      !> output vector
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional               :: scalar
      
      real(kind=rkind), dimension(3) :: ql, qv
      
      real(kind=rkind) :: rhol, rhov

      call liquid_flux(pde(re_ord), layer, quadpnt, vector_out=ql(1:D))
      
      call vapor_flux(pde(re_ord), layer, quadpnt, vector_out=ql(1:D))

      rhol = 1000.0
      
      rhov = rho_sv(quadpnt)
    
      Cl = 

      

      
   end subroutine evap_convect
    
end module evap_heat_fnc
