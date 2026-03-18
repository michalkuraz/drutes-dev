module datetime
  use typy

  ! Public API

  private :: timestamp_to_seconds
  public :: difftime, parse_iso_datetime

  type, public :: datetime_t
    integer(kind=ikind) :: year
    integer(kind=ikind) :: month
    integer(kind=ikind) :: day
    integer(kind=ikind) :: hour
    integer(kind=ikind) :: minute
    integer(kind=ikind) :: second
  end type datetime_t

contains



  !--------------------------------------------------------
  ! Parse ISO 8601 "YYYY-MM-DDTHH:MM:SSZ" into datetime_t
  !--------------------------------------------------------
  subroutine parse_iso_datetime(ts, dt)

    character(len=*), intent(in)  :: ts
    type(datetime_t), intent(out) :: dt

    ! Assumes ts has correct fixed format and length
    read(ts(1:4),  *) dt%year
    read(ts(6:7),  *) dt%month
    read(ts(9:10), *) dt%day
    read(ts(12:13),*) dt%hour
    read(ts(15:16),*) dt%minute
    read(ts(18:19),*) dt%second
  end subroutine parse_iso_datetime

  !--------------------------------------------------------
  ! Leap year test (Gregorian calendar)
  !--------------------------------------------------------
   function is_leap_year(y) result(isit)

    integer(kind=ikind), intent(in) :: y
    logical :: isit

    if (mod(y,4_ikind) == 0_ikind .and. &
                   (mod(y,100_ikind) /= 0_ikind .or. mod(y,400_ikind) == 0_ikind)) then
      isit = .true.
    else
      isit = .false.
    end if


  end function is_leap_year

  !--------------------------------------------------------
  ! Days before year y (from year 0), closed form
  !   days_before_year(1) = 0
  !   days_before_year(2) = days in year 1, etc.
  !--------------------------------------------------------
  function days_before_year(y) result(days)
    integer(kind=ikind), intent(in) :: y
    integer(kind=ikind) :: days

    integer(kind=ikind) :: y0
    y0 = y - 1_ikind
    days = 365_ikind*y0                      &
                     + y0/4_ikind - y0/100_ikind + y0/400_ikind
  end function days_before_year

  !--------------------------------------------------------
  ! Day-of-year (1..365/366) for a given date
  !--------------------------------------------------------
   function day_of_year(y, m, d) result(dy)

    integer(kind=ikind), intent(in) :: y, m, d
    integer(kind=ikind) :: dy

    integer(kind=ikind), dimension(12) :: mdays
    integer(kind=ikind) :: i

    mdays = (/ 31_ikind, 28_ikind, 31_ikind, 30_ikind, 31_ikind, 30_ikind, &
               31_ikind, 31_ikind, 30_ikind, 31_ikind, 30_ikind, 31_ikind /)
    if (is_leap_year(y)) mdays(2) = 29_ikind

    dy = d
    do i = 1_ikind, m-1_ikind
      dy = dy + mdays(i)
    end do

  end function day_of_year

    !--------------------------------------------------------
    ! Convert datetime_t to "seconds since year 0".
    ! Epoch is arbitrary; only differences matter.
    !--------------------------------------------------------
  function timestamp_to_seconds(dt) result(sec)
    use typy
    implicit none
    type(datetime_t), intent(in) :: dt
    integer(kind=ikind)          :: sec

    integer(kind=ikind) :: doy
    integer(kind=ikind) :: days_total

    ! Day-of-year and total day index
    doy        = day_of_year(dt%year, dt%month, dt%day)
    days_total = days_before_year(dt%year) + doy - 1_ikind   ! days since year 0, starting at 0

    ! Integer-only accumulation:
    ! sec = (((days*24 + hour)*60 + minute)*60 + second)
    sec = (((days_total * 24_ikind + dt%hour) * 60_ikind + dt%minute) * 60_ikind) + dt%second

  end function timestamp_to_seconds

  !--------------------------------------------------------
  ! Time difference between two datetime_t values:
  !
  !   dtime = t2 - t1  [in units given by optional 'unit']
  !
  ! If 'unit' is not present or equals "sec"  -> seconds (default)
  ! If 'unit' equals "min"                    -> minutes
  ! If 'unit' equals "hrs"                    -> hours
  !
  ! Any other 'unit' triggers error stop.
  !--------------------------------------------------------
  function difftime(t1, t2, unit) result(dtime)

    type(datetime_t), intent(in)  :: t1, t2
    real(kind=rkind)  :: dtime
    character(len=*), intent(in), optional :: unit

    real(kind=rkind) :: s1, s2, dsec, coeff
    character(len=8) :: u


    s1   = timestamp_to_seconds(t1)
    s2   = timestamp_to_seconds(t2)

    dsec = s2 - s1

    ! Default: seconds
    if (present(unit)) then
       u = adjustl(unit)
       select case (trim(u))
       case ("sec")
          coeff = 1.0_rkind
       case ("min")
          coeff = 1.0_rkind / 60.0_rkind
       case ("hrs")
          coeff = 1.0_rkind / 3600.0_rkind
       case default
          error stop "difftime: unsupported time unit '" // trim(unit) // "'"
       end select
    else
       coeff = 1.0_rkind
    end if

    dtime = dsec * coeff
  end function difftime


end module datetime
