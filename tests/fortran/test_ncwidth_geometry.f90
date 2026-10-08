program test_ncwidth_geometry
  use typy, only: rkind
  use ncwidth_geometry, only: river_contact_width
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  real(kind=rkind) :: cell(4,2), adjacent(4,2,3), d(2), rotated(4,2,3), transformed(4,2)
  real(kind=rkind) :: angle, rotation(2,2), offset(2), nan
  integer :: i, k, checks

  checks = 0
  cell = rectangle(0.0_rkind,2.0_rkind,0.0_rkind,2.0_rkind)
  adjacent(:,:,1) = rectangle(-2.0_rkind,0.0_rkind,0.0_rkind,2.0_rkind)
  adjacent(:,:,2) = rectangle(2.0_rkind,4.0_rkind,0.0_rkind,2.0_rkind)
  d = [1.0_rkind,0.0_rkind]
  call check("straight channel", river_contact_width(cell,d,adjacent(:,:,1:2)), 2.0_rkind)
  call check("direction normalization", river_contact_width(cell,12*d,adjacent(:,:,1:2)), 2.0_rkind)
  d = [1.0_rkind,1.0_rkind]
  call check("oblique full contacts", river_contact_width(cell,d,adjacent(:,:,1:2)), sqrt(2.0_rkind))
  call check("isolated cell span", river_contact_width(cell,d,adjacent(:,:,1:0)), 2*sqrt(2.0_rkind))

  adjacent(:,:,2) = rectangle(2.0_rkind,4.0_rkind,0.5_rkind,1.5_rkind)
  d = [1.0_rkind,0.0_rkind]
  call check("partial outlet", river_contact_width(cell,d,adjacent(:,:,1:2)), 1.0_rkind)
  call check("partial inlet under reversed flow", river_contact_width(cell,-d,adjacent(:,:,1:2)), 1.0_rkind)
  call check("open source end", river_contact_width(cell,d,adjacent(:,:,2:2)), 1.0_rkind)
  d = [1.0_rkind,1.0_rkind]
  call check("partial oblique contact", river_contact_width(cell,d,adjacent(:,:,1:2)), 1/sqrt(2.0_rkind))
  call check("clockwise vertices", &
    river_contact_width(cell([1,4,3,2],:),d,adjacent([1,4,3,2],:,1:2)), 1/sqrt(2.0_rkind))
  call check("neighbour order independent", &
    river_contact_width(cell,d,adjacent(:,:,[2,1])), 1/sqrt(2.0_rkind))

  ! Rotation/translation should not affect a physical width, including UTM-scale offsets.
  offset = [400000.0_rkind,5300000.0_rkind]
  do k = 1, 12
    angle = real(k,rkind)*0.37_rkind
    rotation(1,:) = [cos(angle),-sin(angle)]
    rotation(2,:) = [sin(angle), cos(angle)]
    do i = 1, 4
      transformed(i,:) = matmul(rotation,cell(i,:)) + offset
      rotated(i,:,1) = matmul(rotation,adjacent(i,:,1)) + offset
      rotated(i,:,2) = matmul(rotation,adjacent(i,:,2)) + offset
    end do
    call check("rotated UTM contact", &
      river_contact_width(transformed,matmul(rotation,d),rotated(:,:,1:2)), 1/sqrt(2.0_rkind))
  end do

  adjacent(:,:,2) = rectangle(2.0_rkind,4.0_rkind,0.0_rkind,1.0_rkind)
  adjacent(:,:,3) = rectangle(2.0_rkind,4.0_rkind,1.0_rkind,2.0_rkind)
  d = [1.0_rkind,0.0_rkind]
  call check("parallel branch openings are summed", river_contact_width(cell,d,adjacent), 2.0_rkind)
  adjacent(:,:,1) = rectangle(2.0_rkind,4.0_rkind,2.0_rkind,4.0_rkind)
  d = [1.0_rkind,1.0_rkind]
  call check("corner-only downstream", river_contact_width(cell,d,adjacent(:,:,1:1)), 0.0_rkind)
  call check("corner-only upstream", river_contact_width(cell,-d,adjacent(:,:,1:1)), 0.0_rkind)
  ! A real face on that side supplies an opening; a diagonal neighbour must not close it.
  adjacent(:,:,2) = rectangle(2.0_rkind,4.0_rkind,0.0_rkind,2.0_rkind)
  call check("face plus redundant corner", &
    river_contact_width(cell,d,adjacent(:,:,1:2)), sqrt(2.0_rkind))
  d = [0.0_rkind,1.0_rkind]
  call check("tangent to sole contact", river_contact_width(cell,d,adjacent(:,:,2:2)), 0.0_rkind)
  call check("zero direction", river_contact_width(cell,0*d,adjacent), 0.0_rkind)
  nan = ieee_value(0.0_rkind,ieee_quiet_nan)
  call check("nonfinite direction", river_contact_width(cell,[nan,1.0_rkind],adjacent), 0.0_rkind)
  cell(:,2) = 0.0_rkind
  call check("degenerate cell", river_contact_width(cell,d,adjacent), 0.0_rkind)
  print '(a,i0,a)', "PASS: ", checks, " river contact geometry checks"

contains
  function rectangle(x1,x2,y1,y2) result(vertices)
    real(kind=rkind), intent(in) :: x1,x2,y1,y2
    real(kind=rkind) :: vertices(4,2)
    vertices(:,1) = [x1,x2,x2,x1]
    vertices(:,2) = [y1,y1,y2,y2]
  end function rectangle

  subroutine check(label,actual,expected)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*), intent(in) :: label
    real(kind=rkind), intent(in) :: actual,expected
    if (.not. ieee_is_finite(actual)) error stop "Nonfinite width"
    if (abs(actual-expected) > 1.0e-7_rkind) then
      print *, "FAIL: ", label, "; got ", actual, "; expected ", expected
      error stop 1
    end if
    checks = checks + 1
  end subroutine check
end program test_ncwidth_geometry
