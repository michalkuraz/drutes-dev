! Impermeable solute banks on the boundary of the fixed active ADEnc FE domain.
! With the existing advective volume form, K grad(C).n = (q.n) C.
! DRUtES uses negative mass/stiffness signs: the Robin edge addition is +dt*q.n*Ni*Nj.
! This removes absorbing banks; it does not repair interior flux jumps or H_t+div(q).

module ncboundary

  use typy
  implicit none
  private
  public :: read_adenc_banks, prepare_adenc_banks, adenc_bank_element, bank_edge_terms
  public :: active_node_element
  public :: adenc_open_element

contains

  subroutine read_adenc_banks()

    use ncglobvars, only: LSbank_noflow
    use globals, only: drutes_config
    use ncconfig, only: open_transport_config
    use core_tools, only: write_log
    logical :: exists, conservative, audit, hydro
    integer(kind=ikind) :: unit

    LSbank_noflow = .false.
    call open_transport_config(unit, exists, LSbank_noflow, conservative, audit, hydro)

    if (exists) then

      close (unit)

    end if

    if (.not. LSbank_noflow) then

      return

    end if

    if (drutes_config%dimen /= 2 .or. drutes_config%it_method /= 0) then

      error stop 'ADEnc no-flow banks require 2D and standard Picard (no Schwarz)'

    end if

    call write_log('ADEnc banks: zero TOTAL solute flux; active FE domain only')

  end subroutine read_adenc_banks

  subroutine prepare_adenc_banks(original_edge, close_exterior)

    use ncglobvars, only: LSbank_noflow, bank_edges, ncfluxdata, addedbc, open_edges
    use globals, only: nodes, elements
    use pde_objs, only: pde
    use core_tools, only: write_log
    integer(kind=ikind), dimension(:), intent(in) :: original_edge
    ! Manufactured closed-box tests may additionally close the external mesh.
    ! Production keeps external mesh boundaries separate from internal banks.
    logical, intent(in), optional :: close_exterior
    ! Hash undirected active-triangle edges. A second occurrence is internal.
    integer, dimension(:), allocatable :: head
    integer, dimension(:), allocatable :: next
    integer, dimension(:), allocatable :: lo
    integer, dimension(:), allocatable :: hi
    integer, dimension(:), allocatable :: owner
    integer, dimension(:), allocatable :: local_edge
    logical, dimension(:), allocatable :: participating
    logical, dimension(:, :), allocatable :: interface_edges
    integer :: e, k, a, b, bucket, item, item_count, buckets, nbanks
    integer, dimension(2, 3), parameter :: ends = reshape([1, 2, 2, 3, 3, 1], [2, 3])
    character(len=256) :: message

    if (.not. LSbank_noflow) then

      return

    end if

    if (size(elements%data, 2) /= 3) then

      error stop 'No-flow banks require P1 triangles'

    end if

    if (.not. any(ncfluxdata%activeel)) then

      error stop 'ADEnc active domain is empty'

    end if

    if (allocated(bank_edges)) then

      deallocate (bank_edges)

    end if

    allocate (bank_edges(3, elements%kolik))
    bank_edges = .false.

    if (allocated(pde(1)%assembly_mask)) then

      deallocate (pde(1)%assembly_mask)

    end if

    allocate (pde(1)%assembly_mask(elements%kolik))
    pde(1)%assembly_mask = ncfluxdata%activeel
    allocate (participating(nodes%kolik))
    participating = .false.
    buckets = 2*elements%kolik + 1
    allocate (head(buckets), next(3*elements%kolik), lo(3*elements%kolik), hi(3*elements%kolik), &
              owner(3*elements%kolik), local_edge(3*elements%kolik))
    head = 0
    item_count = 0

    do e = 1, elements%kolik

      if (.not. ncfluxdata%activeel(e)) then

        cycle

      end if

      participating(elements%data(e, :)) = .true.

      do k = 1, 3

        a = minval(elements%data(e, ends(:, k)))
        b = maxval(elements%data(e, ends(:, k)))
        bucket = int(modulo(31_8*int(a, 8) + int(b, 8), int(buckets, 8))) + 1
        item = head(bucket)

        do while (item /= 0)

          if (lo(item) == a .and. hi(item) == b) then

            exit

          end if

          item = next(item)

        end do

        if (item == 0) then

          item_count = item_count + 1
          item = item_count
          lo(item) = a
          hi(item) = b
          owner(item) = e
          local_edge(item) = k
          next(item) = head(bucket)
          head(bucket) = item
          bank_edges(k, e) = .true.

        else

          if (owner(item) == 0) then

            error stop 'Nonmanifold ADEnc active edge'

          end if

          bank_edges(local_edge(item), owner(item)) = .false.
          owner(item) = 0

        end if

      end do

    end do

    allocate (interface_edges(3, elements%kolik))
    interface_edges = .false.

    do e = 1, elements%kolik

      if (ncfluxdata%activeel(e)) then

        cycle

      end if


      do k = 1, 3

        a = minval(elements%data(e, ends(:, k)))
        b = maxval(elements%data(e, ends(:, k)))
        bucket = int(modulo(31_8*int(a, 8) + int(b, 8), int(buckets, 8))) + 1
        item = head(bucket)

        do while (item /= 0)

          if (lo(item) == a .and. hi(item) == b) then

            exit

          end if

          item = next(item)

        end do

        if (item == 0) then

          cycle

        end if

        if (owner(item) == 0) then

          error stop 'Nonmanifold active/inactive ADEnc interface'

        end if

        interface_edges(local_edge(item), owner(item)) = .true.

      end do

    end do

    if (allocated(open_edges)) then

      deallocate (open_edges)

    end if

    allocate (open_edges(3, elements%kolik))
    open_edges = bank_edges .and. .not. interface_edges

    if (present(close_exterior)) then

      if (close_exterior) then

        open_edges = .false.

      end if


    end if

    if (present(close_exterior)) then

      if (.not. close_exterior) then

        bank_edges = bank_edges .and. interface_edges

      end if


    else

      bank_edges = bank_edges .and. interface_edges

    end if

    ! Matching original physical IDs identify ports, not banks. Do not use the
    ! threshold-generated addedbc here, and do not infer edges from node pairs
    ! inside a triangle (all three nodes may lie on different bank edges).

    do e = 1, elements%kolik

      do k = 1, 3

        a = elements%data(e, ends(1, k))
        b = elements%data(e, ends(2, k))

        if (original_edge(a) > 100 .and. original_edge(a) == original_edge(b)) then

          bank_edges(k, e) = .false.
          open_edges(k, e) = .false.

        end if

      end do

    end do

    do a = 1, nodes%kolik

      if (participating(a)) then

        nodes%edge(a) = original_edge(a) ! bank nodes regain DOFs; preserve real ports

      else

        nodes%edge(a) = addedbc ! unused nodes remain eliminated, never assembled

      end if

    end do

    nbanks = count(bank_edges)
    write (message, *) 'ADEnc no-flow domain: elements=', count(ncfluxdata%activeel), &
      ' nodes=', count(participating), ' bank edges=', nbanks
    call write_log(trim(message))

  end subroutine prepare_adenc_banks

  function active_node_element(node_id) result(el)

    use globals, only: nodes
    use ncglobvars, only: LSbank_noflow, ncfluxdata
    integer(kind=ikind), intent(in) :: node_id
    integer(kind=ikind) :: el
    integer(kind=ikind) :: j

    el = nodes%element(node_id)%data(1)

    if (.not. LSbank_noflow) then

      return

    end if

    ! smartarray data includes spare, uninitialized capacity. Only pos entries
    ! are adjacent elements (especially important for inactive-only nodes).

    do j = 1, nodes%element(node_id)%pos

      if (ncfluxdata%activeel(nodes%element(node_id)%data(j))) then

        el = nodes%element(node_id)%data(j)
        return

      end if

    end do

  end function active_node_element

  pure subroutine bank_edge_terms(edge, length, normal_flux, dt, matrix)

    integer, dimension(2), intent(in) :: edge
    real(kind=rkind), intent(in) :: length
    real(kind=rkind), dimension(2), intent(in) :: normal_flux
    real(kind=rkind), intent(in) :: dt
    real(kind=rkind), dimension(3, 3), intent(out) :: matrix
    real(kind=rkind) :: s
    real(kind=rkind), dimension(3) :: basis
    integer :: g, i, j

    matrix = 0

    do g = 1, 2

      s = (1 + (2*g - 3)/sqrt(3.0_rkind))/2
      basis = 0
      basis(edge(1)) = 1 - s
      basis(edge(2)) = s

      do j = 1, 3

        do i = 1, 3

          matrix(i, j) = matrix(i, j) + dt*length/2*normal_flux(g)*basis(i)*basis(j)

        end do

      end do

    end do

  end subroutine bank_edge_terms

  subroutine adenc_bank_element(pde_loc, el_id, dt, quadpnt_in)

    use global_objs, only: integpnt_str
    use globals, only: stiff_mat, time, nodes, elements
    use pde_objs, only: pde_str
    use ncglobvars, only: LSbank_noflow, bank_edges, ncfluxdata, ora_di_ini
    use netcdfflux, only: ncflux_get_xy_cell
    use ncfluxarea, only: ncflux_active_width
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el_id
    real(kind=rkind), intent(in) :: dt
    type(integpnt_str), intent(in), optional :: quadpnt_in
    integer, dimension(2, 3), parameter :: ends = reshape([1, 2, 2, 3, 3, 1], [2, 3])
    integer :: k, g
    real(kind=rkind), dimension(2) :: a
    real(kind=rkind), dimension(2) :: b
    real(kind=rkind), dimension(2) :: center
    real(kind=rkind), dimension(2) :: point
    real(kind=rkind), dimension(2) :: normal
    real(kind=rkind) :: length
    real(kind=rkind) :: s
    real(kind=rkind) :: Q
    real(kind=rkind) :: W
    real(kind=rkind), dimension(2) :: qn
    real(kind=rkind), dimension(3, 3) :: matrix
    logical :: ok
    character(len=1024) :: message

    if (.not. LSbank_noflow) then

      return

    end if

    if (.not. ncfluxdata%activeel(el_id)) then

      return

    end if

    W = ncflux_active_width(el_id)

    if (W <= 0) then

      ! matches zero hydrological convection callback
      return

    end if

    center = sum(nodes%data(elements%data(el_id, :), 1:2), dim=1)/3

    do k = 1, 3

      if (.not. bank_edges(k, el_id)) then

        cycle

      end if

      a = nodes%data(elements%data(el_id, ends(1, k)), 1:2)
      b = nodes%data(elements%data(el_id, ends(2, k)), 1:2)
      length = norm2(b - a)

      if (length <= 0) then

        error stop 'Degenerate ADEnc bank edge'

      end if

      normal = [b(2) - a(2), a(1) - b(1)]/length

      if (dot_product(normal, center - (a + b)/2) > 0) then

        normal = -normal

      end if


      do g = 1, 2

        s = (1 + (2*g - 3)/sqrt(3.0_rkind))/2
        point = (1 - s)*a + s*b
        ! One-sided active-element trace at hydrological grid boundaries.
        point = point + 1.e-8_rkind*(center - point)
        call ncflux_get_xy_cell(point(1), point(2), ora_di_ini + int(time/86400.0_rkind)*24, Q, ok, message)
        qn(g) = 0

        if (ok .and. Q > 0) then

          qn(g) = Q/W*dot_product(ncfluxdata%fluxvct(el_id, :), normal)

        end if


      end do

      call bank_edge_terms(ends(:, k), length, qn, dt, matrix)
      stiff_mat = stiff_mat + matrix

    end do

  end subroutine adenc_bank_element

  ! Conservative form: J.n=0 is natural on banks; exterior outflow has
  ! zero diffusive flux and retains q.n*C. Unspecified inflow is an error.

  subroutine adenc_open_element(el_id, dt, matrix)

    use globals, only: nodes, elements
    use ncglobvars, only: open_edges, ncfluxdata, ora_di_ini, adenc_coefficient_time
    use netcdfflux, only: ncflux_get_xy_cell
    use ncfluxarea, only: ncflux_active_width
    use nchydroflow, only: LShydro, hydro_normal_flux
    integer(kind=ikind), intent(in) :: el_id
    real(kind=rkind), intent(in) :: dt
    real(kind=rkind), dimension(3, 3), intent(out) :: matrix
    integer, dimension(2, 3), parameter :: ends = reshape([1, 2, 2, 3, 3, 1], [2, 3])
    real(kind=rkind), dimension(2) :: a
    real(kind=rkind), dimension(2) :: b
    real(kind=rkind), dimension(2) :: center
    real(kind=rkind), dimension(2) :: normal
    real(kind=rkind), dimension(2) :: point
    real(kind=rkind) :: length
    real(kind=rkind) :: s
    real(kind=rkind) :: Q
    real(kind=rkind) :: W
    real(kind=rkind), dimension(2) :: qn
    real(kind=rkind), dimension(3, 3) :: part
    integer :: k, g
    logical :: ok
    character(len=1024) :: message

    matrix = 0

    if (.not. allocated(open_edges)) then

      error stop 'Conservative exterior topology not initialized'

    end if

    W = ncflux_active_width(el_id)

    if (W <= 0) then

      return

    end if

    center = sum(nodes%data(elements%data(el_id, :), 1:2), dim=1)/3

    do k = 1, 3

      if (.not. open_edges(k, el_id)) then

        cycle

      end if

      a = nodes%data(elements%data(el_id, ends(1, k)), 1:2)
      b = nodes%data(elements%data(el_id, ends(2, k)), 1:2)
      length = norm2(b - a)

      if (length <= 0) then

        error stop 'Degenerate conservative exterior edge'

      end if

      normal = [b(2) - a(2), a(1) - b(1)]/length

      if (dot_product(normal, center - (a + b)/2) > 0) then

        normal = -normal

      end if


      if (LShydro) then

        qn = hydro_normal_flux(int(el_id), k)

        if (any(qn < -1.e-12_rkind)) then

          error stop 'Reconstructed outlet became inflow'

        end if

        call bank_edge_terms(ends(:, k), length, qn, dt, part)
        matrix = matrix - part
        cycle

      end if

      do g = 1, 2

        s = (1 + (2*g - 3)/sqrt(3.0_rkind))/2
        point = (1 - s)*a + s*b
        point = point + 1.e-8_rkind*(center - point)
        call ncflux_get_xy_cell(point(1), point(2), &
                                ora_di_ini + int(adenc_coefficient_time()/86400.0_rkind)*24, Q, ok, message)
        qn(g) = 0

        if (ok .and. Q > 0) then

          qn(g) = Q/W*dot_product(ncfluxdata%fluxvct(el_id, :), normal)

        end if


        if (qn(g) < -1.e-12_rkind) then

          print *, 'Unspecified ADEnc exterior inflow at element/edge:', el_id, k
          error stop 'Label exterior inflow with an explicit Dirichlet inlet'

        end if

      end do

      call bank_edge_terms(ends(:, k), length, qn, dt, part)
      matrix = matrix - part ! negative equation sign: -dt integral Ni*Jn

    end do

  end subroutine

end module ncboundary
