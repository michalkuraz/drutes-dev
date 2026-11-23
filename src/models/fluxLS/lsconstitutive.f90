module lsconstitutive

  use ncboundary, only: active_node_element
  use nchydroflow, only: hydro_value

contains

  subroutine ncflux(pde_loc, layer, quadpnt, x, grad, flux, flux_length)

    use typy
    use global_objs
    use pde_objs
    use geom_tools
    use ncglobvars, only: ncfluxdata

    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in)                          :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), dimension(:), optional                   :: x
    !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
    real(kind=rkind), dimension(:), intent(in), optional     :: grad
    real(kind=rkind), dimension(:), intent(out), optional    :: flux
    real(kind=rkind), intent(out), optional                  :: flux_length

    real(kind=rkind), dimension(2) :: xy, q
    real(kind=rkind) :: h, divq
    integer(kind=ikind) :: el

    el = hydraulic_element(quadpnt)
    q = 0.0_rkind

    if (ncfluxdata%activeel(el)) then

      call getcoor(quadpnt, xy)
      call hydro_value(int(el), xy, q, h, divq)

    end if

    if (present(flux)) flux = q
    if (present(flux_length)) flux_length = norm2(q)

  end subroutine ncflux

  function hydraulic_element(point) result(el)

    use typy, only: ikind
    use global_objs, only: integpnt_str
    type(integpnt_str), intent(in), optional :: point
    integer(kind=ikind) :: el

    if (.not. present(point)) error stop 'ADEnc hydraulic coefficients require an integration point'

    select case (point%type_pnt)

    case ('gqnd', 'obpt')

      el = point%element

    case ('ndpt')

      el = active_node_element(point%order)

    case default

      error stop 'Unsupported ADEnc hydraulic integration point type'

    end select

  end function hydraulic_element

  subroutine ADEls_convection(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)

    use typy
    use global_objs
    use pde_objs
    use ADE_globals
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

    if (present(vector_out)) then

      call pde(1)%flux(layer, quadpnt, vector_out=vector_out)

    end if

    if (present(scalar)) then

      call pde(1)%flux(layer, quadpnt, scalar=scalar)

    end if

  end subroutine ADEls_convection

  ! First-order loss of suspended inventory H*C. FEM supplies the factor C.
  function ADEls_deposition(pde_loc, layer, quadpnt, x) result(val)

    use typy
    use pde_objs
    use ncglobvars, only: adencpar, ncfluxdata
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), optional, intent(in) :: quadpnt
    real(kind=rkind), dimension(:), optional, intent(in) :: x
    real(kind=rkind) :: val
    integer(kind=ikind) :: el

    if (quadpnt%type_pnt == 'ndpt') then

      el = active_node_element(quadpnt%order)

    else

      el = quadpnt%element

    end if

    val = 0.0_rkind
    if (ncfluxdata%activeel(el)) then

      val = -adencpar(layer)%k_dep * pde_loc%pde_fnc(1)%elasticity(pde_loc, layer, quadpnt, x)

    end if

  end function ADEls_deposition

  function ADEls_tder_coef(pde_loc, layer, quadpnt, x) result(val)

    use typy
    use global_objs
    use pde_objs
    use ncglobvars, only: ncfluxdata
    use geom_tools

    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind)                :: val

    real(kind=rkind), dimension(2) :: xy, q
    real(kind=rkind) :: divq
    integer(kind=ikind) :: el

    el = hydraulic_element(quadpnt)
    val = 1.0_rkind

    if (ncfluxdata%activeel(el)) then

      call getcoor(quadpnt, xy)
      call hydro_value(int(el), xy, q, val, divq)

    end if

  end function ADEls_tder_coef

  function ADEls_mass(pde_loc, layer, quadpnt, x) result(val)

    use typy
    use global_objs
    use pde_objs

    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return value
    real(kind=rkind)                :: val

    val = pde(1)%getval(quadpnt)*ADEls_tder_coef(pde_loc, layer, quadpnt, x)

  end function ADEls_mass

  subroutine ADEls_icond(pde_loc)

    use typy
    use globals
    use global_objs
    use pde_objs
    use ncglobvars

    class(pde_str), intent(in out) :: pde_loc
    integer(kind=ikind) :: i, j, k, l, m, layer, D
    real(kind=rkind) :: value

    pde_loc%solution(:) = cinit_ls

  end subroutine ADEls_icond

  subroutine ADEls_dirichlet(pde_loc, el_id, node_order, value, code, array, bcpts)

    use typy
    use global_objs
    use pde_objs
    use ncglobvars, only: adenc_coefficient_time
    use ncconservative, only: adenc_boundary_value
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in)  :: el_id, node_order
    real(kind=rkind), intent(out), optional    :: value
    integer(kind=ikind), intent(out), optional :: code
    !> unused for this model (implementation for Robin boundary)
    real(kind=rkind), dimension(:), intent(out), optional :: array
    type(bcpts_str), intent(in), optional :: bcpts

    if (present(value)) value = adenc_boundary_value(pde_loc, el_id, node_order, adenc_coefficient_time())
    if (present(code)) code = 1

  end subroutine ADEls_dirichlet

  subroutine ADElsdisp(pde_loc, layer, quadpnt, x, tensor, scalar)

    use typy
    use global_objs
    use pde_objs
    use globals
    use ADE_globals
    use re_globals
    use debug_tools
    use ncglobvars
    use ncdispersion, only: ls_dispersion_tensor

    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    !> return tensor
    real(kind=rkind), dimension(:, :), intent(out), optional :: tensor
    !> relative scalar value of the nonlinear function
    real(kind=rkind), intent(out), optional                 :: scalar

    integer(kind=ikind) :: D
    real(kind=rkind), dimension(2) :: qvector

    if (present(quadpnt) .and. present(x)) then

      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from lsconstitutive::ADElsdisp"
      ERROR stop

    else if (.not. present(quadpnt) .and. .not. present(x)) then

      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from lsconstitutive::ADElsdisp"
      ERROR stop

    end if

    D = drutes_config%dimen

    if (D /= 2) then

      error stop "ADEnc dispersion requires two dimensions"

    end if

    if (.not. present(quadpnt)) then

      error stop "ADEnc dispersion requires an integration point"

    end if

    call pde(1)%flux(layer, quadpnt, vector_out=qvector)

    if (present(tensor)) then

      call ls_dispersion_tensor(qvector, adencpar(layer)%alpha_l, adencpar(layer)%alpha_t, tensor)

    end if

    if (present(scalar)) then

      ! Scalar interface reports the longitudinal coefficient, not the tensor.
      scalar = adencpar(layer)%alpha_l*norm2(qvector)

    end if

  end subroutine ADElsdisp

end module lsconstitutive
