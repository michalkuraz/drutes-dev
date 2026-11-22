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

  subroutine parse_ls_dispersivity(line, alpha_l, alpha_t, ok, deposition, value_count)

    character(len=*), intent(in) :: line
    real(kind=rkind), intent(out) :: alpha_l, alpha_t
    logical, intent(out) :: ok
    real(kind=rkind), optional, intent(out) :: deposition
    integer(kind=ikind) :: i, first, last, count, ierr, hash
    real(kind=rkind), dimension(3) :: values

    integer(kind=ikind), optional, intent(out) :: value_count

    if (present(value_count)) value_count = 0
    ok = .false.
    if (present(deposition)) then

      deposition = 0.0_rkind

    end if
    alpha_l = 0.0_rkind
    alpha_t = 0.0_rkind
    last = len_trim(line)
    hash = index(line, '#')

    if (hash > 0) then

      last = hash - 1

    end if

    i = 1
    count = 0

    do while (i <= last)

      if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then

        i = i + 1
        cycle

      end if

      first = i

      do while (i <= last)

        if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then

          exit

        end if

        i = i + 1

      end do

      count = count + 1

      if (count > size(values) .or. (count == 3 .and. .not. present(deposition))) then

        return

      end if

      ! Reject list-directed null/repeat/termination syntax and extra values.

      if (scan(line(first:i - 1), ',/*') > 0) then

        return

      end if

      read (line(first:i - 1), *, iostat=ierr) values(count)

      if (ierr /= 0) then

        return

      end if

      if (.not. ieee_is_finite(values(count))) then

        return

      end if

      if (values(count) < 0.0_rkind) then

        return

      end if


    end do

    if (count == 0) then

      return

    end if

    alpha_l = values(1)
    alpha_t = values(1)

    if (count >= 2) then

      alpha_t = values(2)

    end if

    if (count == 3) then

      deposition = values(3)

    end if

    if (present(value_count)) value_count = count
    ok = .true.

  end subroutine parse_ls_dispersivity

  ! q is the depth-integrated flux [m2/s], not velocity [m/s].
  ! K = |q| [alpha_T I + (alpha_L-alpha_T) d d^T] [m3/s].
  ! Dividing by the existing depth/storage coefficient gives D [m2/s].

  pure subroutine ls_dispersion_tensor(q, alpha_l, alpha_t, tensor)

    real(kind=rkind), dimension(:), intent(in) :: q
    real(kind=rkind), intent(in) :: alpha_l
    real(kind=rkind), intent(in) :: alpha_t
    real(kind=rkind), dimension(:, :), intent(out) :: tensor
    real(kind=rkind) :: magnitude
    real(kind=rkind), dimension(size(q)) :: direction
    integer(kind=ikind) :: i, j

    tensor = 0.0_rkind
    magnitude = norm2(q)

    if (magnitude <= 0.0_rkind) then

      return

    end if

    direction = q/magnitude

    do i = 1, size(q)

      do j = 1, size(q)

        tensor(i, j) = magnitude*(alpha_l - alpha_t)*direction(i)*direction(j)

      end do

      tensor(i, i) = tensor(i, i) + magnitude*alpha_t

    end do

  end subroutine ls_dispersion_tensor

end module ncdispersion
