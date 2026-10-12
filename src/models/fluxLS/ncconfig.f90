! Shared input records for the NetCDF transport model.

module ncconfig

  use typy, only: ikind
  implicit none
  private
  public :: open_transport_config, config_flag, config_end, reject_old_config

contains

  subroutine reject_old_config(path)

    character(len=*), intent(in) :: path
    logical :: exists

    inquire (file=path, exist=exists)

    if (exists) then

      print *, 'Obsolete configuration file: ', trim(path)
      error stop 'Use transport.conf and stabilization.conf; separate legacy files are no longer supported'

    end if

  end subroutine

  subroutine config_record(unit, line, status)

    integer(kind=ikind), intent(in) :: unit
    character(len=*), intent(out) :: line
    integer(kind=ikind), intent(out) :: status
    integer(kind=ikind) :: marker

    do

      read (unit, '(A)', iostat=status) line

      if (status /= 0) then

        return

      end if

      marker = index(line, '#')

      if (marker > 0) then

        line = line(:marker - 1)

      end if

      if (len_trim(line) > 0) then

        return

      end if


    end do

  end subroutine

  subroutine config_flag(unit, value)

    integer(kind=ikind), intent(in) :: unit
    logical, intent(out) :: value
    integer(kind=ikind) :: status
    character(len=1024) :: line

    call config_record(unit, line, status)

    if (status /= 0) then

      error stop 'Missing y/n flag in combined NetCDF configuration'

    end if


    select case (trim(adjustl(line)))

    case ('y')

      value = .true.

    case ('n')

      value = .false.

    case default

      error stop 'Combined NetCDF configuration flags must be a single y or n'

    end select

  end subroutine

  subroutine config_end(unit)

    integer(kind=ikind), intent(in) :: unit
    integer(kind=ikind) :: status
    character(len=1024) :: line

    call config_record(unit, line, status)

    if (status == 0) then

      error stop 'Unexpected trailing combined NetCDF configuration record'

    end if

    if (status > 0) then

      error stop 'Error reading combined NetCDF configuration'

    end if


  end subroutine

  ! The caller owns the open unit, positioned immediately after the four flags.

  subroutine open_transport_config(unit, exists, banks, conservative, audit, hydro)

    integer(kind=ikind), intent(out) :: unit
    logical, intent(out) :: exists, banks, conservative, audit, hydro
    integer(kind=ikind) :: status

    call reject_old_config('drutes.conf/netcdf/riverbank.conf')
    call reject_old_config('drutes.conf/netcdf/conservative.conf')
    call reject_old_config('drutes.conf/netcdf/hydroflow.conf')
    banks = .false.
    conservative = .false.
    audit = .false.
    hydro = .false.
    unit = 0
    inquire (file='drutes.conf/netcdf/transport.conf', exist=exists)

    if (.not. exists) then

      return

    end if

    open (newunit=unit, file='drutes.conf/netcdf/transport.conf', status='old', action='read', iostat=status)

    if (status /= 0) then

      error stop 'Cannot open drutes.conf/netcdf/transport.conf'

    end if

    call config_flag(unit, banks)
    call config_flag(unit, conservative)
    call config_flag(unit, audit)
    call config_flag(unit, hydro)

    if (.not. hydro) then

      call config_end(unit)

    end if


  end subroutine

end module ncconfig
