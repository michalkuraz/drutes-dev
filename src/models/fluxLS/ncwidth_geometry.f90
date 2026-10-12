! Geometry for converting cell discharge to discharge per active transverse width.

module ncwidth_geometry

  use typy, only: rkind, ikind
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: river_contact_width
  public :: river_cell_width

contains

  ! Homogenized ADEnc width: full hydrological-cell span normal to flow [m].
  ! This is NOT a physical wetted-channel width or a projected face opening.
  ! Use it consistently for cell-averaged Q/W and H=Q/(W*v). A near-tangential
  ! neighbour contact must not collapse storage width over the entire cell.

  pure function river_cell_width(vertices, direction) result(width)

    real(kind=rkind), dimension(4, 2), intent(in) :: vertices
    real(kind=rkind), dimension(2), intent(in) :: direction
    real(kind=rkind) :: width
    real(kind=rkind), dimension(2) :: d
    real(kind=rkind), dimension(2) :: normal
    real(kind=rkind), dimension(4, 2) :: local
    real(kind=rkind), dimension(4) :: projections
    real(kind=rkind), dimension(2) :: centre
    real(kind=rkind) :: norm_d
    real(kind=rkind) :: tolerance
    real(kind=rkind) :: orientation
    integer(kind=ikind) :: i, next_i

    width = 0.0_rkind

    if (.not. all(ieee_is_finite(vertices))) then

      return

    end if

    if (.not. all(ieee_is_finite(direction))) then

      return

    end if

    norm_d = norm2(direction)

    if (.not. ieee_is_finite(norm_d) .or. norm_d <= tiny(1.0_rkind)) then

      return

    end if

    d = direction/norm_d
    normal = [-d(2), d(1)]
    centre = sum(vertices, dim=1)/4.0_rkind

    do i = 1, 4

      local(i, :) = vertices(i, :) - centre
      projections(i) = dot_product(local(i, :), normal)

    end do

    tolerance = max(1.0e-9_rkind, &
                    64.0_rkind*epsilon(1.0_rkind)*max(1.0_rkind, maxval(abs(vertices))))
    orientation = 0.0_rkind

    do i = 1, 4

      next_i = mod(i, 4) + 1
      orientation = orientation + cross2(local(i, :), local(next_i, :))

    end do

    if (abs(orientation) <= tolerance*maxval(abs(local))) then

      return

    end if

    width = maxval(projections) - minval(projections)

    if (width <= tolerance .or. .not. ieee_is_finite(width)) then

      width = 0.0_rkind

    end if


  end function river_cell_width

  ! Effective width, not the wetted width of a surveyed river cross-section.
  ! Retained as a contact diagnostic/legacy API, NOT ADEnc's storage width.
  ! Each shared edge contributes length * abs(dot(flow, outward_edge_normal)).
  ! Sum contacts on each side of the cell (branches are parallel openings),
  ! then limit the full transverse span by the narrower inlet/outlet.
  ! A missing upstream/downstream neighbour is an open channel end. A known
  ! neighbour touching only at a corner has zero opening, not an open end.
  ! An isolated cell retains its own span: no connecting width can be inferred.

  pure function river_contact_width(vertices, direction, neighbours) result(width)

    real(kind=rkind), dimension(4, 2), intent(in) :: vertices
    real(kind=rkind), dimension(2), intent(in) :: direction
    real(kind=rkind), dimension(:, :, :), intent(in) :: neighbours
    real(kind=rkind) :: width
    real(kind=rkind), dimension(2) :: d
    real(kind=rkind), dimension(2) :: normal
    real(kind=rkind), dimension(2) :: centre
    real(kind=rkind), dimension(2) :: other_centre
    real(kind=rkind), dimension(2) :: edge
    real(kind=rkind), dimension(2) :: outward
    real(kind=rkind), dimension(4) :: projections
    real(kind=rkind), dimension(4, 2) :: local
    real(kind=rkind) :: length, shared, overlap, total_shared, inlet, outlet, signed_width
    real(kind=rkind) :: orientation, tolerance, norm_d, along
    integer(kind=ikind) :: i, j, k, next_i, next_j
    logical :: touches, edge_touches, corner_inlet, corner_outlet, face_contact

    width = 0.0_rkind

    if (.not. all(ieee_is_finite(vertices))) then

      return

    end if

    if (.not. all(ieee_is_finite(direction))) then

      return

    end if

    if (size(neighbours, 1) /= 4 .or. size(neighbours, 2) /= 2) then

      return

    end if

    norm_d = norm2(direction)

    if (norm_d <= tiny(1.0_rkind)) then

      return

    end if

    d = direction/norm_d
    normal = [-d(2), d(1)]
    centre = sum(vertices, dim=1)/4.0_rkind
    ! Use local coordinates to avoid subtracting large UTM projections.

    do i = 1, 4

      local(i, :) = vertices(i, :) - centre
      projections(i) = dot_product(local(i, :), normal)

    end do

    tolerance = max(1.0e-9_rkind, &
                    64.0_rkind*epsilon(1.0_rkind)*max(1.0_rkind, maxval(abs(vertices))))
    orientation = 0.0_rkind

    do i = 1, 4

      next_i = mod(i, 4) + 1
      orientation = orientation + cross2(local(i, :), local(next_i, :))

    end do

    if (abs(orientation) <= tolerance*maxval(abs(local))) then

      return

    end if

    width = maxval(projections) - minval(projections)

    if (width <= tolerance) then

      width = 0.0_rkind
      return

    end if

    inlet = 0.0_rkind
    outlet = 0.0_rkind
    corner_inlet = .false.
    corner_outlet = .false.
    face_contact = .false.

    do k = 1, size(neighbours, 3)

      if (.not. all(ieee_is_finite(neighbours(:, :, k)))) then

        cycle

      end if

      touches = .false.
      total_shared = 0.0_rkind

      do i = 1, 4

        next_i = mod(i, 4) + 1
        edge = vertices(next_i, :) - vertices(i, :)
        length = norm2(edge)

        if (length <= tolerance) then

          cycle

        end if

        shared = 0.0_rkind

        do j = 1, 4

          next_j = mod(j, 4) + 1
          call shared_edge_overlap(vertices(i, :), vertices(next_i, :), &
                                   neighbours(j, :, k), neighbours(next_j, :, k), tolerance, overlap, edge_touches)
          shared = shared + overlap
          touches = touches .or. edge_touches

        end do

        shared = min(shared, length)

        if (shared <= tolerance) then

          cycle

        end if

        total_shared = total_shared + shared
        outward = sign(1.0_rkind, orientation)*[edge(2), -edge(1)]/length
        signed_width = shared*dot_product(d, outward)

        if (signed_width > tolerance) then

          outlet = outlet + signed_width

        end if

        if (signed_width < -tolerance) then

          inlet = inlet - signed_width

        end if


      end do

      if (total_shared > tolerance) then

        face_contact = .true.

      else if (touches) then

        other_centre = sum(neighbours(:, :, k), dim=1)/4.0_rkind
        along = dot_product(other_centre - centre, d)

        if (along < -tolerance) then

          corner_inlet = .true.

        end if

        if (along > tolerance) then

          corner_outlet = .true.

        end if


      end if

    end do

    if (inlet > tolerance) then

      width = min(width, inlet)

    else if (corner_inlet) then

      width = 0.0_rkind

    end if

    if (outlet > tolerance) then

      width = min(width, outlet)

    else if (corner_outlet) then

      width = 0.0_rkind

    end if

    ! Flow tangent to every real contact cannot pass through those contacts.

    if (face_contact .and. max(inlet, outlet) <= tolerance) then

      width = 0.0_rkind

    end if


  end function river_contact_width

  ! Overlap of collinear segments, including partial-face and point contacts.

  pure subroutine shared_edge_overlap(a, b, c, d, tolerance, length, touches)

    real(kind=rkind), dimension(2), intent(in) :: a
    real(kind=rkind), dimension(2), intent(in) :: b
    real(kind=rkind), dimension(2), intent(in) :: c
    real(kind=rkind), dimension(2), intent(in) :: d
    real(kind=rkind), intent(in) :: tolerance
    real(kind=rkind), intent(out) :: length
    logical, intent(out) :: touches
    real(kind=rkind), dimension(2) :: tangent
    real(kind=rkind) :: edge_length
    real(kind=rkind) :: lo
    real(kind=rkind) :: hi
    real(kind=rkind) :: tc
    real(kind=rkind) :: td

    length = 0.0_rkind
    touches = .false.
    tangent = b - a
    edge_length = norm2(tangent)

    if (edge_length <= tolerance) then

      return

    end if

    tangent = tangent/edge_length

    if (abs(cross2(tangent, c - a)) > tolerance) then

      return

    end if

    if (abs(cross2(tangent, d - a)) > tolerance) then

      return

    end if

    tc = dot_product(c - a, tangent)
    td = dot_product(d - a, tangent)
    lo = max(0.0_rkind, min(tc, td))
    hi = min(edge_length, max(tc, td))
    touches = hi >= lo - tolerance

    if (touches) then

      length = max(0.0_rkind, hi - lo)

    end if


  end subroutine shared_edge_overlap

  pure function cross2(a, b) result(value)

    real(kind=rkind), dimension(2), intent(in) :: a
    real(kind=rkind), dimension(2), intent(in) :: b
    real(kind=rkind) :: value

    value = a(1)*b(2) - a(2)*b(1)

  end function cross2

end module ncwidth_geometry
