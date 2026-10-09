! Integrate the real ncfluxarea cache with small synthetic ncglobvars arrays.
! Does not call DRUtES main, initialize a live model, or touch configuration/output files.
program test_ncflux_widths
  use typy
  use ncglobvars
  use ncfluxarea, only: ncflux_prepare_widths, ncflux_active_width
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  integer(kind=ikind) :: row, col, cell, node_id, i, zero_width_count
  integer :: checks
  logical :: ok
  character(len=1024) :: errmsg

  checks = 0
  call check("unprepared cache",ncflux_active_width(1_ikind),0.0_rkind)
  ncfluxdata%initialized = .true.
  ncfluxdata%slice_loaded = .true.
  ncfluxdata%nlat = 3
  ncfluxdata%nlon = 3
  ncfluxdata%has_fill = .true.
  Qmin = 0.0_rkind
  ncelements%kolik = 9
  ncnodes%kolik = 36
  allocate(ncelements%data(9,4),ncnodes%data(36,2))
  allocate(ncfluxdata%qslice(3,3),ncfluxdata%activeel(7),ncfluxdata%fluxvct(7,2),el2ncgrid(7))
  do row = 1, 3
    do col = 1, 3
      cell = (row-1)*3 + col
      node_id = (cell-1)*4
      ncelements%data(cell,:) = node_id + [1,2,3,4]
      ncnodes%data(node_id+1:node_id+4,1) = real([col-1,col,col,col-1],rkind)
      ncnodes%data(node_id+1:node_id+4,2) = real([row-1,row-1,row,row],rkind)
    end do
  end do
  ncfluxdata%qslice = 0.0_rkind
  ncfluxdata%qslice(2,:) = [10.0_rkind,20.0_rkind,30.0_rkind]
  ncfluxdata%qslice(1,2) = ncfluxdata%fill_value
  ncfluxdata%qslice(3,2) = ieee_value(0.0_rkind,ieee_quiet_nan)
  el2ncgrid = [5,1,0,10,5,6,5]
  ncfluxdata%activeel = .true.
  ncfluxdata%activeel(5) = .false.
  do i = 1, 7
    ncfluxdata%fluxvct(i,:) = [1.0_rkind,1.0_rkind]
  end do
  ncfluxdata%fluxvct(6,:) = [1.0_rkind,0.0_rkind]
  ncfluxdata%fluxvct(7,:) = 0.0_rkind
  call prepare()
  if (zero_width_count /= 1) error stop "Incorrect count of blocked flowing elements"
  call check("full oblique cell span",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  call check("dry mapped cell",ncflux_active_width(2_ikind),0.0_rkind)
  call check("unmapped FE element",ncflux_active_width(3_ikind),0.0_rkind)
  call check("invalid mapping",ncflux_active_width(4_ikind),0.0_rkind)
  call check("inactive FE element",ncflux_active_width(5_ikind),0.0_rkind)
  call check("channel end",ncflux_active_width(6_ikind),1.0_rkind)
  call check("undefined flow direction",ncflux_active_width(7_ikind),0.0_rkind)
  call check("negative element index",ncflux_active_width(-1_ikind),0.0_rkind)
  call check("element index beyond cache",ncflux_active_width(8_ikind),0.0_rkind)

  ncnodes%data(:,2) = 3.0_rkind - ncnodes%data(:,2)
  ncfluxdata%fluxvct(:,2) = -ncfluxdata%fluxvct(:,2)
  call prepare()
  call check("descending latitude order",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  ncnodes%data(:,1) = 3.0_rkind - ncnodes%data(:,1)
  ncfluxdata%fluxvct(:,1) = -ncfluxdata%fluxvct(:,1)
  call prepare()
  call check("descending longitude order",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  ncnodes%data = 3.0_rkind - ncnodes%data
  ncfluxdata%fluxvct = -ncfluxdata%fluxvct
  call prepare()

  ! A later slice must not silently redefine the initial fixed channel geometry.
  ncfluxdata%qslice = 0.0_rkind
  ncfluxdata%qslice(2,2) = 20.0_rkind
  call check("width cache is fixed",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  call prepare()
  call check("explicit cache rebuild",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  Qmin = 21.0_rkind
  call prepare()
  call check("threshold excludes cell",ncflux_active_width(1_ikind),0.0_rkind)
  Qmin = 0.0_rkind

  ! A corner-only neighbour is not a finite flow opening, but it must NOT
  ! collapse the full-cell storage width. Connectivity is handled separately.
  ncfluxdata%qslice(3,3) = 10.0_rkind
  call prepare()
  call check("corner neighbour does not shrink storage width",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  ! Additional active neighbours must not redefine the cell storage volume.
  ncfluxdata%qslice(2,3) = 10.0_rkind
  call prepare()
  call check("face neighbour does not shrink storage width",ncflux_active_width(1_ikind),sqrt(2.0_rkind))
  ncfluxdata%fluxvct(1,:)=[1.0_rkind,1.e-8_rkind]
  call prepare()
  call check("near-tangential contacts retain cell-scale width",ncflux_active_width(1_ikind),1+1.e-8_rkind)

  ncfluxdata%nlon = 4
  call ncflux_prepare_widths(ok,errmsg)
  if (ok) error stop "Invalid array dimensions were accepted"
  call check("failed rebuild clears stale cache",ncflux_active_width(1_ikind),0.0_rkind)
  print '(a,i0,a)', "PASS: ", checks, " river width integration checks"
contains
  subroutine prepare()
    call ncflux_prepare_widths(ok,errmsg,zero_width_count)
    if (.not. ok) then
      print *, trim(errmsg)
      error stop 1
    end if
  end subroutine prepare
  subroutine check(label,actual,expected)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*), intent(in) :: label
    real(kind=rkind), intent(in) :: actual,expected
    if (.not. ieee_is_finite(actual)) error stop "Nonfinite width"
    if (abs(actual-expected) > 1.0e-10_rkind) then
      print *, "FAIL: ", label, "; got ", actual, "; expected ", expected
      error stop 1
    end if
    checks = checks + 1
  end subroutine check
end program test_ncflux_widths
