! Optional conservative ADEnc local assembly using the existing post-capacity hook.

module ncconservative

  use typy
  use ncglobvars
  implicit none
  private
  public :: read_adenc_conservative, conservative_begin, conservative_end, conservative_element
  public :: adenc_boundary_history, adenc_boundary_value

contains

  subroutine read_adenc_conservative()

    use globals, only: drutes_config
    use ncconfig, only: open_transport_config
    use core_tools, only: write_log
    use ncbalance, only: balance_reset
    integer(kind=ikind) :: unit
    logical :: exists, banks, hydro

    LSconservative = .false.
    LSbalance = .false.
    LSstep_active = .false.
    LSclock_override = .false.
    LSstate_time = 0
    LSprevious_time = 0
    call balance_reset()

    if (allocated(LSdepth_old)) then

      deallocate (LSdepth_old, LSdepth_new)

    end if

    call open_transport_config(unit, exists, banks, LSconservative, LSbalance, hydro)

    if (exists) then

      close (unit)

    end if

    if (.not. (LSconservative .or. LSbalance)) then

      return

    end if

    if (drutes_config%run_from_backup) then

      error stop 'ADEnc conservative/audit restart from backup is not supported'

    end if

    if (drutes_config%dimen /= 2 .or. drutes_config%it_method /= 0 .or. .not. LSbank_noflow) then

      error stop 'ADEnc conservative/audit requires 2D standard Picard and no-flow active banks'

    end if

    call write_log('ADEnc optional conservative transport / inventory audit configured')

  end subroutine

  function adenc_boundary_value(pde_loc, el, node, t) result(value)

    use pde_objs, only: pde_str
    use globals, only: nodes, elements
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el, node
    real(kind=rkind), intent(in) :: t
    real(kind=rkind) :: value
    integer(kind=ikind) :: edge, j, selected

    edge = nodes%edge(elements%data(el, node))
    value = pde_loc%bc(edge)%value

    if (.not. pde_loc%bc(edge)%file) then

      return

    end if

    selected = 1

    do j = 1, size(pde_loc%bc(edge)%series, 1)

      if (pde_loc%bc(edge)%series(j, 1) > t) then

        exit

      end if

      selected = j

    end do

    value = pde_loc%bc(edge)%series(selected, 2) ! last value persists after final record

  end function

  function adenc_boundary_history(pde_loc, el, node, column) result(value)

    use pde_objs, only: pde_str
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el, node, column
    real(kind=rkind) :: value, t

    t = LSstate_time

    if (column == 4) then

      t = LSprevious_time

    end if

    if (LSstep_active .and. column /= 1 .and. column /= 4) then

      t = LStrial_time

    end if

    value = adenc_boundary_value(pde_loc, el, node, t)

  end function

  subroutine conservative_begin(pde_loc, t, dt)

    use pde_objs
    use globals, only: elements, gauss_points
    use global_objs, only: integpnt_str
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use nchydroflow, only: hydro_begin, hydro_clip_step
    class(pde_str), intent(in) :: pde_loc
    real(kind=rkind), intent(in) :: t
    real(kind=rkind), intent(in out) :: dt
    type(integpnt_str) :: point
    integer(kind=ikind) :: el, g, edge, j
    real(kind=rkind) :: event, next_day

    if (pde_common%timeint_method /= 1 .and. pde_common%timeint_method /= 2) then

      error stop 'ADEnc conservative/audit requires transient implicit Euler'

    end if

    if (size(pde) /= 1 .or. size(elements%data, 2) /= 3 .or. .not. pde_loc%diffusion) then

      error stop 'ADEnc conservative/audit requires a single 2D P1 equation'

    end if


    if (LSconservative) then

      do edge = lbound(pde_loc%bc, 1), ubound(pde_loc%bc, 1)

        if (pde_loc%bc(edge)%code /= 1) then

          error stop 'ADEnc conservative requires Dirichlet port records'

        end if

        if (.not. pde_loc%bc(edge)%file) then

          cycle

        end if

        if (.not. allocated(pde_loc%bc(edge)%series)) then

          error stop 'Missing ADEnc inlet series'

        end if

        if (size(pde_loc%bc(edge)%series, 2) /= 2 .or. size(pde_loc%bc(edge)%series, 1) < 1) then

          error stop 'ADEnc inlet series must have two columns'

        end if

        if (any(.not. ieee_is_finite(pde_loc%bc(edge)%series))) then

          error stop 'Nonfinite ADEnc inlet series'

        end if


        do j = 2, size(pde_loc%bc(edge)%series, 1)

          if (pde_loc%bc(edge)%series(j, 1) <= pde_loc%bc(edge)%series(j - 1, 1)) then

            error stop 'ADEnc inlet times must be strictly increasing'

          end if


        end do

      end do

    end if

    if (LSconservative) then

      next_day = (floor(t/86400.0_rkind) + 1.0_rkind)*86400.0_rkind
      dt = min(dt, next_day - t)

      do edge = lbound(pde_loc%bc, 1), ubound(pde_loc%bc, 1)

        if (.not. pde_loc%bc(edge)%file) then

          cycle

        end if


        do j = 1, size(pde_loc%bc(edge)%series, 1)

          event = pde_loc%bc(edge)%series(j, 1)

          if (event > t) then

            dt = min(dt, event - t)

          end if


        end do

      end do

    end if

    call hydro_clip_step(t, dt)

    if (.not. ieee_is_finite(dt) .or. dt <= 0 .or. t + dt <= t) then

      error stop 'Invalid ADEnc trial time step'

    end if

    LSstep_start = t
    LSstep_dt = dt
    ! Left limit at discontinuities: a step ending at pulse/day change uses
    ! the interval's forcing. Accepted H and BC history retain that SAME trace.
    LStrial_time = nearest(t + dt, -1.0_rkind)

    if (.not. LSconservative) then

      ! diagnostic-only preserves legacy timing
      LStrial_time = t

    end if


    if (.not. allocated(LSdepth_old)) then

      allocate (LSdepth_old(size(gauss_points%weight), elements%kolik), &
                LSdepth_new(size(gauss_points%weight), elements%kolik))
      LSdepth_old = 1
      LSdepth_new = 1
      point%type_pnt = 'gqnd'
      point%column = 1
      LSclock_override = .true.
      LSoverride_time = LSstate_time

      do el = 1, elements%kolik

        if (.not. ncfluxdata%activeel(el)) then

          cycle

        end if

        point%element = el

        do g = 1, size(gauss_points%weight)

          point%order = g
          LSdepth_old(g, el) = pde_loc%pde_fnc(1)%elasticity(pde_loc, elements%material(el), point)

          if (.not. ieee_is_finite(LSdepth_old(g, el)) .or. LSdepth_old(g, el) <= 0) then

            error stop 'Invalid ADEnc initial storage depth'

          end if


        end do

      end do

      LSclock_override = .false.

    end if

    call hydro_begin(LStrial_time, dt, t)
    LSstep_active = .true.

  end subroutine

  subroutine conservative_end(pde_loc, accepted)

    use pde_objs, only: pde_str
    use ncbalance, only: balance_finish
    use nchydroflow, only: hydro_end
    class(pde_str), intent(in) :: pde_loc
    logical, intent(in) :: accepted

    if (LSbalance) then

      call balance_finish(pde_loc, accepted, LSstep_start, LSstep_dt)

    end if


    if (accepted) then

      LSdepth_old = LSdepth_new
      LSprevious_time = LSstate_time
      LSstate_time = LStrial_time

    end if

    LSstep_active = .false.
    call hydro_end(accepted)

  end subroutine

  subroutine conservative_element(pde_loc, el, dt, oldcap, newcap, open_matrix, source)

    use pde_objs
    use globals
    use global_objs
    use ncboundary, only: adenc_open_element
    use nchydroflow, only: hydro_sources
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el
    real(kind=rkind), intent(in) :: dt
    real(kind=rkind), dimension(3, 3), intent(out) :: oldcap
    real(kind=rkind), dimension(3, 3), intent(out) :: newcap
    real(kind=rkind), dimension(3, 3), intent(out) :: open_matrix
    real(kind=rkind), intent(out) :: source
    type(integpnt_str) :: point
    real(kind=rkind), dimension(2) :: q
    real(kind=rkind) :: h
    real(kind=rkind) :: w
    real(kind=rkind) :: tmp
    real(kind=rkind), dimension(3) :: row_sum
    real(kind=rkind) :: water
    real(kind=rkind) :: lateral_load
    integer :: i, j, g

    if (.not. LSstep_active) then

      error stop 'ADEnc assembly without step preparation'

    end if

    oldcap = 0
    newcap = cap_mat
    open_matrix = 0
    ! time_integ already added the old-concentration capacity load to bside.
    ! Only the remaining load represents a physical source, not inventory.
    source = -sum(bside - matmul(newcap, elnode_prev))
    call hydro_sources(int(el), water, lateral_load)
    source = source + dt*elements%areas(el)*lateral_load
    point%type_pnt = 'gqnd'
    point%element = el
    point%column = 2

    do g = 1, size(gauss_points%weight)

      point%order = g
      h = pde_loc%pde_fnc(1)%elasticity(pde_loc, elements%material(el), point)

      if (.not. ieee_is_finite(h) .or. h <= 0) then

        error stop 'Invalid ADEnc storage depth'

      end if

      LSdepth_new(g, el) = h
      w = elements%areas(el)*gauss_points%weight(g)/gauss_points%area
      ! Explicit local P1 load: shared zerord assembly does not weight by Ni.
      ! Keep this ADEnc source local and audit its physical (positive) amount.
      bside = bside - dt*w*base_fnc(:, g)*lateral_load
      call pde_loc%pde_fnc(1)%convection(pde_loc, elements%material(el), point, vector_out=q)

      do j = 1, 3

        do i = 1, 3

          oldcap(i, j) = oldcap(i, j) - w*LSdepth_old(g, el)*base_fnc(i, g)*base_fnc(j, g)

          if (LSconservative) then

            ! Remove -Ni*q.grad(Nj), add +grad(Ni).q*Nj (negative DRUtES sign).
            tmp = base_fnc(i, g)*dot_product(elements%ders(el, j, 1:2), q) + &
                  dot_product(elements%ders(el, i, 1:2), q)*base_fnc(j, g)
            stiff_mat(i, j) = stiff_mat(i, j) + dt*w*tmp

          end if

        end do

      end do

    end do

    if (pde_common%timeint_method == 1) then

      row_sum = sum(oldcap, dim=2)
      oldcap = 0

      do i = 1, 3

        oldcap(i, i) = row_sum(i)

      end do

    end if

    if (LSconservative) then

      bside = bside + matmul(oldcap - newcap, elnode_prev)
      call adenc_open_element(el, dt, open_matrix)
      stiff_mat = stiff_mat + open_matrix

    end if

  end subroutine

end module
