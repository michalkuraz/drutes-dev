! Direction-aligned mechanical dispersion for the depth-integrated ADEnc model.
module ncdispersion
  use typy
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: parse_ls_dispersivity, ls_dispersion_tensor
contains
  ! One legacy value means alpha_L = alpha_T. Two values mean alpha_L alpha_T.
  ! Parse one physical record only, never consuming the following Qmin record.
  subroutine parse_ls_dispersivity(line, alpha_l, alpha_t, ok)
    character(len=*), intent(in) :: line
    real(kind=rkind), intent(out) :: alpha_l, alpha_t
    logical, intent(out) :: ok
    integer :: i, first, last, count, ierr, hash
    real(kind=rkind) :: values(2)
    ok = .false.
    alpha_l = 0.0_rkind
    alpha_t = 0.0_rkind
    last = len_trim(line)
    hash = index(line, '#')
    if (hash > 0) last = hash-1
    i = 1
    count = 0
    do while (i <= last)
      if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
        i = i+1
        cycle
      end if
      first = i
      do while (i <= last)
        if (line(i:i) == ' ' .or. line(i:i) == achar(9)) exit
        i = i+1
      end do
      count = count+1
      if (count > 2) return
      ! Reject list-directed null/repeat/termination syntax and extra values.
      if (scan(line(first:i-1), ',/*') > 0) return
      read(line(first:i-1), *, iostat=ierr) values(count)
      if (ierr /= 0) return
      if (.not. ieee_is_finite(values(count))) return
      if (values(count) < 0.0_rkind) return
    end do
    if (count == 0) return
    alpha_l = values(1)
    alpha_t = values(1)
    if (count == 2) alpha_t = values(2)
    ok = .true.
  end subroutine parse_ls_dispersivity

  ! q is the depth-integrated flux [m2/s], not velocity [m/s].
  ! K = |q| [alpha_T I + (alpha_L-alpha_T) d d^T] [m3/s].
  ! Dividing by the existing depth/storage coefficient gives D [m2/s].
  pure subroutine ls_dispersion_tensor(q, alpha_l, alpha_t, tensor)
    real(kind=rkind), intent(in) :: q(:), alpha_l, alpha_t
    real(kind=rkind), intent(out) :: tensor(:,:)
    real(kind=rkind) :: magnitude, direction(size(q))
    integer :: i, j
    tensor = 0.0_rkind
    magnitude = norm2(q)
    if (magnitude <= 0.0_rkind) return
    direction = q/magnitude
    do i=1,size(q)
      do j=1,size(q)
        tensor(i,j) = magnitude*(alpha_l-alpha_t)*direction(i)*direction(j)
      end do
      tensor(i,i) = tensor(i,i)+magnitude*alpha_t
    end do
  end subroutine ls_dispersion_tensor
end module ncdispersion
