! Exercise the real PDE dispersion callback without initializing/running a model.
module test_dispersion_flux
  use typy
  use global_objs
  use pde_objs
  implicit none
  real(kind=rkind) :: test_q(2) = [3.0_rkind, -4.0_rkind]
contains
  subroutine sample_flux(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:), vector_in(:)
    real(kind=rkind), intent(out), optional :: vector_out(:), scalar
    if (present(vector_out)) vector_out = test_q
    if (present(scalar)) scalar = norm2(test_q)
  end subroutine sample_flux
end module test_dispersion_flux

program test_adenc_dispersion_callback
  use typy
  use global_objs
  use pde_objs
  use globals
  use ncglobvars
  use lsconstitutive, only: ADElsdisp
  use test_dispersion_flux
  implicit none
  type(integpnt_str) :: point
  real(kind=rkind) :: tensor(2,2), scalar
  integer :: unit, ierr
  character(len=128) :: line
  logical :: ok
  ! Production comment skipping + line parser must not consume the Qmin row.
  block
    use ncdispersion, only: parse_ls_dispersivity
    use readtools, only: comment, fileread
    real(kind=rkind) :: threshold
    open(newunit=unit, status='scratch', action='readwrite')
    write(unit,'(A)') '# header'
    write(unit,'(A)') ''
    write(unit,'(A)') '200 0.2'
    write(unit,'(A)') '# Qmin'
    write(unit,'(A)') '300'
    rewind(unit)
    call comment(unit)
    read(unit,'(A)',iostat=ierr) line
    if (ierr/=0) error stop 'dispersion record read'
    call parse_ls_dispersivity(line,LSdisp,LSdisp_transverse,ok)
    if (.not. ok) error stop 'dispersion record parse'
    call fileread(threshold,unit)
    if (threshold/=300) error stop 'Qmin record shifted'
    close(unit)
  end block
  allocate(pde(1))
  pde(1)%flux => sample_flux
  drutes_config%dimen = 2
  point%type_pnt = 'gqnd'
  point%element = 1
  call ADElsdisp(pde(1),1_ikind,point,tensor=tensor,scalar=scalar)
  if (abs(scalar-1000)>1e-10_rkind) error stop 'scalar longitudinal coefficient'
  if (abs(tensor(1,1)-360.64_rkind)>1e-10_rkind) error stop 'Kxx'
  if (abs(tensor(2,2)-640.36_rkind)>1e-10_rkind) error stop 'Kyy'
  if (abs(tensor(1,2)+479.52_rkind)>1e-10_rkind) error stop 'Kxy'
  if (abs(tensor(2,1)-tensor(1,2))>1e-10_rkind) error stop 'Kyx'
  test_q = 0
  call ADElsdisp(pde(1),1_ikind,point,tensor=tensor,scalar=scalar)
  if (any(tensor/=0) .or. scalar/=0) error stop 'zero flux callback'
  print *, 'ADEnc dispersion callback checks passed'
end program test_adenc_dispersion_callback
