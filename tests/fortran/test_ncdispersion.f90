program test_ncdispersion
  use typy
  use ncdispersion
  implicit none
  real(kind=rkind) :: a, b, tensor(2,2), reverse(2,2), q(2), d(2), normal(2)
  logical :: ok
  integer :: i
  character(len=24) :: bad(9)
  call parse_ls_dispersivity('2000.0', a, b, ok)
  if (.not. ok .or. a /= 2000 .or. b /= 2000) error stop 'legacy parsing'
  call parse_ls_dispersivity('200.0'//achar(9)//'0.2 # comment', a, b, ok)
  if (.not. ok .or. a /= 200 .or. abs(b-0.2_rkind)>1e-12_rkind) error stop 'pair parsing'
  bad = [character(len=24) :: '', '-1 2', '1 -2', 'NaN 1', '1 Inf', '1 2 3', '1, 2', '2*1', '1 /']
  do i=1,size(bad)
    call parse_ls_dispersivity(trim(bad(i)), a, b, ok)
    if (ok) error stop 'invalid record accepted'
  end do
  call ls_dispersion_tensor([2.0_rkind,0.0_rkind],200.0_rkind,0.2_rkind,tensor)
  if (abs(tensor(1,1)-400)>1e-10_rkind .or. abs(tensor(2,2)-0.4_rkind)>1e-10_rkind) &
    error stop 'axis eigenvalues'
  if (tensor(1,2)/=0 .or. tensor(2,1)/=0) error stop 'axis off-diagonal'
  q = [3.0_rkind,-4.0_rkind]
  d = q/5
  normal = [-d(2),d(1)]
  call ls_dispersion_tensor(q,200.0_rkind,0.2_rkind,tensor)
  if (maxval(abs(matmul(tensor,d)-1000*d))>1e-10_rkind) error stop 'longitudinal eigenvector'
  if (maxval(abs(matmul(tensor,normal)-normal))>1e-10_rkind) error stop 'transverse eigenvector'
  if (abs(tensor(1,2)-tensor(2,1))>1e-10_rkind) error stop 'symmetry'
  call ls_dispersion_tensor(-q,200.0_rkind,0.2_rkind,reverse)
  if (maxval(abs(reverse-tensor))>1e-10_rkind) error stop 'reversed flow'
  call ls_dispersion_tensor(q,2000.0_rkind,2000.0_rkind,tensor)
  if (tensor(1,1)/=10000 .or. tensor(2,2)/=10000 .or. tensor(1,2)/=0) error stop 'isotropic regression'
  call ls_dispersion_tensor([0.0_rkind,0.0_rkind],200.0_rkind,0.2_rkind,tensor)
  if (any(tensor/=0)) error stop 'zero flow'
  call ls_dispersion_tensor(q,0.0_rkind,0.0_rkind,tensor)
  if (any(tensor/=0)) error stop 'zero dispersivity'
  print *, 'ADEnc dispersion checks passed'
end program test_ncdispersion
