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
    
    
    subroutine evap_heatconvect(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
      use evap_RE_fnc
      use evapextras
      use evapglob
      
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
      integer(kind=ikind) :: D
      
      D = drutes_config%dimen

      call liquid_flux(pde(re_ord), layer, quadpnt, flux=ql(1:D))
      
      call vapor_flux(pde(re_ord), layer, quadpnt, flux=qv(1:D))

      rhol = rho_l( quadpnt)
      
      rhov = rho_sv(quadpnt)*rh_soil(layer, quadpnt)
    
      if (present(vector_out)) then
        vector_out(1:D) = -(C_liq*rhol*ql(1:D) + C_vap*rhov*qv(1:D))
      end if
      
      if (present(scalar)) then
        scalar = norm2(-(C_liq*rhol*ql(1:D) + C_vap*rhov*qv(1:D)))
      end if
      
   end subroutine evap_heatconvect
   
   
   function Ldtheta_vdt(pde_loc, layer, quadpnt, x) result(val)
     use typy
     use global_objs
     use globals
     use pde_objs
     use evapglob
     use evapextras
     use evap_RE_fnc

     class(pde_str), intent(in) :: pde_loc
     !> value of the nonlinear function
     real(kind=rkind), dimension(:), intent(in), optional    :: x
     !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
     type(integpnt_str), intent(in), optional :: quadpnt
     !> material ID
     integer(kind=ikind), intent(in) :: layer
     !> return value
     real(kind=rkind)                :: val
     
     real(kind=rkind) :: L
     
     L = latent_heat_wat(quadpnt)
     
     val = L*dtheta_vdt(pde(re_ord), layer, quadpnt, x)
     
   end function Ldtheta_vdt  
     
    
end module evap_heat_fnc
