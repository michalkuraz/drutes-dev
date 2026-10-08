program test_adenc_boundaries
  use typy
  use globals, only: nodes, elements, time, end_time
  use ncglobvars, only: addedbc
  use pde_objs, only: pde
  use init_netcdf, only: set_adenc_boundary_id, read_adenc_boundaries
  use lsconstitutive, only: ADEls_dirichlet
  implicit none
  character(len=32) :: argument
  integer :: unit, mode
  real(kind=rkind) :: value

  call get_command_argument(1, argument)
  read(argument,*) mode
  allocate(nodes%edge(3), elements%data(1,3), pde(1))
  if (mode == 1) then
    nodes%edge = (/101_ikind, 0_ikind, 0_ikind/)
  else
    nodes%edge = (/101_ikind, 102_ikind, 0_ikind/)
  end if
  call set_adenc_boundary_id()
  if (addedbc /= 101_ikind+mode) error stop "incorrect inactive boundary ID"
  ! Exercise both a used inactive boundary and the all-active domain case.
  if (mode == 2) nodes%edge(3) = addedbc
  elements%data(1,:) = (/1_ikind, 2_ikind, 3_ikind/)
  end_time = 10.0_rkind
  time = 5.0_rkind
  open(newunit=unit, file="boundaries.conf", status="old", action="read")
  call read_adenc_boundaries(unit)
  close(unit)
  if (lbound(pde(1)%bc,1) /= 101 .or. ubound(pde(1)%bc,1) /= addedbc) &
    error stop "incorrect boundary array bounds"
  call ADEls_dirichlet(pde(1), 1_ikind, 1_ikind, value=value)
  if (abs(value-1.0_rkind) > 1.0e-12_rkind) error stop "incorrect first inlet value"
  if (mode == 2) then
    call ADEls_dirichlet(pde(1), 1_ikind, 2_ikind, value=value)
    if (abs(value-2.0_rkind) > 1.0e-12_rkind) error stop "incorrect second inlet value"
    call ADEls_dirichlet(pde(1), 1_ikind, 3_ikind, value=value)
    if (value /= 0.0_rkind) error stop "incorrect inactive boundary value"
  end if
  print *, "ADEnc boundary checks: OK"
end program test_adenc_boundaries
