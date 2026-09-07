! Geometry for converting cell discharge to discharge per active transverse width.
module ncwidth_geometry
  use typy, only: rkind
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: river_contact_width

contains

  ! Effective width, not the wetted width of a surveyed river cross-section.
  ! Each shared edge contributes length * abs(dot(flow, outward_edge_normal)).
  ! Sum contacts on each side of the cell (branches are parallel openings),
  ! then limit the full transverse span by the narrower inlet/outlet.
  ! A missing upstream/downstream neighbour is an open channel end. A known
  ! neighbour touching only at a corner has zero opening, not an open end.
  ! An isolated cell retains its own span: no connecting width can be inferred.
  pure function river_contact_width(vertices, direction, neighbours) result(width)
    real(kind=rkind), intent(in) :: vertices(4,2), direction(2), neighbours(:,:,:)
    real(kind=rkind) :: width
    real(kind=rkind) :: d(2), normal(2), centre(2), other_centre(2)
    real(kind=rkind) :: edge(2), outward(2), projections(4), local(4,2)
    real(kind=rkind) :: length, shared, overlap, total_shared, inlet, outlet, signed_width
    real(kind=rkind) :: orientation, tolerance, norm_d, along
    integer :: i, j, k, next_i, next_j
    logical :: touches, edge_touches, corner_inlet, corner_outlet, face_contact

    width = 0.0_rkind
    if (.not. all(ieee_is_finite(vertices))) return
    if (.not. all(ieee_is_finite(direction))) return
    if (size(neighbours,1) /= 4 .or. size(neighbours,2) /= 2) return
    norm_d = norm2(direction)
    if (norm_d <= tiny(1.0_rkind)) return
    d = direction/norm_d
    normal = [-d(2), d(1)]
    centre = sum(vertices, dim=1)/4.0_rkind
    ! Use local coordinates to avoid subtracting large UTM projections.
    do i = 1, 4
      local(i,:) = vertices(i,:) - centre
      projections(i) = dot_product(local(i,:), normal)
    end do
    tolerance = max(1.0e-9_rkind, &
      64.0_rkind*epsilon(1.0_rkind)*max(1.0_rkind, maxval(abs(vertices))))
    orientation = 0.0_rkind
    do i = 1, 4
      next_i = mod(i,4) + 1
      orientation = orientation + cross2(local(i,:), local(next_i,:))
    end do
    if (abs(orientation) <= tolerance*maxval(abs(local))) return
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
    do k = 1, size(neighbours,3)
      if (.not. all(ieee_is_finite(neighbours(:,:,k)))) cycle
      touches = .false.
      total_shared = 0.0_rkind
      do i = 1, 4
        next_i = mod(i,4) + 1
        edge = vertices(next_i,:) - vertices(i,:)
        length = norm2(edge)
        if (length <= tolerance) cycle
        shared = 0.0_rkind
        do j = 1, 4
          next_j = mod(j,4) + 1
          call shared_edge_overlap(vertices(i,:), vertices(next_i,:), &
            neighbours(j,:,k), neighbours(next_j,:,k), tolerance, overlap, edge_touches)
          shared = shared + overlap
          touches = touches .or. edge_touches
        end do
        shared = min(shared, length)
        if (shared <= tolerance) cycle
        total_shared = total_shared + shared
        outward = sign(1.0_rkind, orientation)*[edge(2), -edge(1)]/length
        signed_width = shared*dot_product(d, outward)
        if (signed_width > tolerance) outlet = outlet + signed_width
        if (signed_width < -tolerance) inlet = inlet - signed_width
      end do
      if (total_shared > tolerance) then
        face_contact = .true.
      else if (touches) then
        other_centre = sum(neighbours(:,:,k), dim=1)/4.0_rkind
        along = dot_product(other_centre - centre, d)
        if (along < -tolerance) corner_inlet = .true.
        if (along > tolerance) corner_outlet = .true.
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
    if (face_contact .and. max(inlet,outlet) <= tolerance) width = 0.0_rkind
  end function river_contact_width


  ! Overlap of collinear segments, including partial-face and point contacts.
  pure subroutine shared_edge_overlap(a, b, c, d, tolerance, length, touches)
    real(kind=rkind), intent(in) :: a(2), b(2), c(2), d(2), tolerance
    real(kind=rkind), intent(out) :: length
    logical, intent(out) :: touches
    real(kind=rkind) :: tangent(2), edge_length, lo, hi, tc, td

    length = 0.0_rkind
    touches = .false.
    tangent = b - a
    edge_length = norm2(tangent)
    if (edge_length <= tolerance) return
    tangent = tangent/edge_length
    if (abs(cross2(tangent,c-a)) > tolerance) return
    if (abs(cross2(tangent,d-a)) > tolerance) return
    tc = dot_product(c-a,tangent)
    td = dot_product(d-a,tangent)
    lo = max(0.0_rkind,min(tc,td))
    hi = min(edge_length,max(tc,td))
    touches = hi >= lo - tolerance
    if (touches) length = max(0.0_rkind,hi-lo)
  end subroutine shared_edge_overlap


  pure function cross2(a,b) result(value)
    real(kind=rkind), intent(in) :: a(2), b(2)
    real(kind=rkind) :: value
    value = a(1)*b(2) - a(2)*b(1)
  end function cross2

end module ncwidth_geometry
