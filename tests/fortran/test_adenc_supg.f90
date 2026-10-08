module supg_test_coefficients
  use typy
  use global_objs
  use pde_objs
  implicit none
  real(kind=rkind) :: qtest(2)=[2.0_rkind,1.0_rkind], depthtest=3.0_rkind
  real(kind=rkind) :: diffusiontest=.03_rkind
  real(kind=rkind), allocatable :: nodaltest(:)
contains
  function current_value(pde_loc,quadpnt) result(value)
    class(pde_str), intent(in) :: pde_loc
    type(integpnt_str), intent(in) :: quadpnt
    real(kind=rkind) :: value
    if (quadpnt%column/=2 .or. quadpnt%type_pnt/='ndpt') error stop 'shock iterate context'
    value=nodaltest(quadpnt%order)
  end function current_value
  subroutine convection(pde_loc,layer,quadpnt,x,vector_in,vector_out,scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:),vector_in(:)
    real(kind=rkind), intent(out), optional :: vector_out(:),scalar
    if (present(vector_out)) vector_out=qtest
    if (present(scalar)) scalar=norm2(qtest)
  end subroutine convection
  subroutine dispersion(pde_loc,layer,quadpnt,x,tensor,scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:)
    real(kind=rkind), intent(out), optional :: tensor(:,:),scalar
    if (present(tensor)) then
      tensor=0
      tensor(1,1)=diffusiontest
      tensor(2,2)=diffusiontest
    end if
    if (present(scalar)) scalar=diffusiontest
  end subroutine dispersion
  function storage(pde_loc,layer,quadpnt,x) result(value)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:)
    real(kind=rkind) :: value
    value=depthtest
  end function storage
end module supg_test_coefficients

program test_adenc_supg
  use typy
  use global_objs
  use globals
  use pde_objs
  use dummy_procs
  use ncsupg
  use ncglobvars, only: LSsupg,LSsupg_factor,LSshock,LSshock_factor
  use stiffmat
  use capmat
  use supg_test_coefficients
  implicit none
  real(kind=rkind) :: g(3,2), basis(3), tensor(2,2), mass(3,3), spatial(3,3), rhs(3)
  real(kind=rkind) :: tau,h,speed,d(2),test(3),rate,source,dt,residual,old(3),new(3)
  real(kind=rkind) :: before_cap(3,3),before_stiff(3,3),before_rhs(3),expected_cap(3,3)
  real(kind=rkind) :: expected_stiff(3,3),expected_rhs(3),integral_mass(3,3)
  integer :: i,j,l,method
  character(len=24) :: mode
  call get_command_argument(1,mode)
  if (trim(mode)=='strip-on' .or. trim(mode)=='strip-off' .or. trim(mode)=='strip-shock') then
    call strip_benchmark(trim(mode)/='strip-off')
    stop
  end if
  if (len_trim(mode)>0) then
    call read_adenc_supg()
    select case(trim(mode))
    case('absent','off')
      if (LSsupg) error stop 'SUPG should be off'
    case('on')
      if (.not. LSsupg .or. LSsupg_factor/=1) error stop 'SUPG should be on'
    case('shock')
      if (.not. LSshock .or. LSshock_factor/=1) error stop 'shock should be on'
    case('shock-absent','shock-off')
      if (LSshock) error stop 'shock should be off'
    case('invalid')
      error stop 'Invalid SUPG settings accepted'
    end select
    print *, 'SUPG configuration checks passed'
    stop
  end if
  g(:,1)=[-1.0_rkind,1.0_rkind,0.0_rkind]
  g(:,2)=[-1.0_rkind,0.0_rkind,1.0_rkind]
  basis=[0.2_rkind,0.3_rkind,0.5_rkind]
  tensor=0; tensor(1,1)=.03_rkind; tensor(2,2)=.03_rkind
  dt=.5_rkind; rate=-.2_rkind; source=.3_rkind
  speed=norm2(qtest/depthtest); d=qtest/norm2(qtest)
  h=2/sum(abs(matmul(g,d)))
  tau=supg_tau(qtest/depthtest,tensor/depthtest,g,dt,.true.,rate/depthtest)
  call close_scalar(tau,1/sqrt((2/dt)**2+(2*speed/h)**2+ &
       (4*.01_rkind/h**2)**2+(rate/depthtest)**2),'tau scales')
  if (tau>dt/2) error stop 'temporal tau bound'
  call close_scalar(supg_tau([0.0_rkind,0.0_rkind],tensor,g,dt,.true.,0.0_rkind),0.0_rkind,'zero velocity')
  call supg_terms(qtest,depthtest,tensor,g,basis,dt,.true.,rate,source,1.0_rkind,mass,spatial,rhs)
  test=tau*matmul(g,qtest/depthtest)
  old=[.2_rkind,.3_rkind,.4_rkind]; new=[.5_rkind,.7_rkind,.1_rkind]
  residual=depthtest*dot_product(basis,new-old)/dt+dot_product(matmul(g,qtest),new) &
           -rate*dot_product(basis,new)-source
  ! Independent strong-residual check covers signs, old-time RHS, source and reaction.
  if (maxval(abs(matmul(mass+spatial,new)-matmul(mass,old)-rhs+dt*test*residual))>1e-12_rkind) &
    error stop 'residual identity'
  if (maxval(abs(sum(mass,dim=1)))>1e-12_rkind) error stop 'test partition temporal'
  if (maxval(abs(sum(spatial,dim=1)))>1e-12_rkind) error stop 'test partition spatial'
  call supg_terms(qtest,depthtest,tensor,g,basis,dt,.false.,0.0_rkind,0.0_rkind,1.0_rkind,mass,spatial,rhs)
  if (any(mass/=0)) error stop 'steady mass correction'
  call supg_terms(qtest,depthtest,tensor,g,basis,dt,.true.,rate,source,0.0_rkind,mass,spatial,rhs)
  if (any(mass/=0) .or. any(spatial/=0) .or. any(rhs/=0)) error stop 'zero multiplier'
  call close_scalar(shock_diffusivity([1.0_rkind,0.0_rkind],g,[1.0_rkind,0.0_rkind], &
                    .2_rkind,1.0_rkind),.1_rkind,'residual viscosity')
  call close_scalar(shock_diffusivity([1.0_rkind,0.0_rkind],g,[1.0_rkind,0.0_rkind], &
                    10.0_rkind,1.0_rkind),.5_rkind,'advective cap')
  call close_scalar(shock_diffusivity(qtest,g,[0.0_rkind,0.0_rkind],1.0_rkind,1.0_rkind), &
                    0.0_rkind,'constant gradient')
  call close_scalar(shock_diffusivity(qtest,g,[1.0_rkind,0.0_rkind],0.0_rkind,1.0_rkind), &
                    0.0_rkind,'zero residual')
  call close_scalar(shock_diffusivity([0.0_rkind,0.0_rkind],g,[1.0_rkind,0.0_rkind], &
                    1.0_rkind,1.0_rkind),0.0_rkind,'zero flow shock')
  call close_scalar(shock_diffusivity([1.0_rkind,0.0_rkind],g,[10.0_rkind,0.0_rkind], &
                    2.0_rkind,1.0_rkind),.1_rkind,'concentration scale invariance')
  call shock_terms(qtest,depthtest,g,basis,new,old,dt,.true.,rate,source,1.0_rkind,spatial)
  residual=(depthtest*dot_product(basis,new-old)/dt+dot_product(qtest,matmul(transpose(g),new)) &
           -rate*dot_product(basis,new)-source)/depthtest
  tau=shock_diffusivity(qtest/depthtest,g,matmul(transpose(g),new),residual,1.0_rkind)
  if (maxval(abs(spatial+dt*depthtest*tau*matmul(g,transpose(g))))>1e-12_rkind) &
    error stop 'shock residual assembly identity'
  if (maxval(abs(sum(spatial,dim=1)))>1e-12_rkind) error stop 'shock partition'
  if (dot_product(new,matmul(spatial,new))>1e-12_rkind) error stop 'shock diffusion sign'
  call shock_terms(qtest,depthtest,g,basis,new,new,dt,.true.,rate,source,0.0_rkind,spatial)
  if (any(spatial/=0)) error stop 'shock disabled'

  allocate(pde(1))
  allocate(pde(1)%pde_fnc(1))
  if (associated(pde(1)%stabilize_element)) error stop 'default callback must be null'
  pde(1)%pde_fnc(1)%convection=>convection
  pde(1)%pde_fnc(1)%dispersion=>dispersion
  pde(1)%pde_fnc(1)%elasticity=>storage
  pde(1)%pde_fnc(1)%reaction=>dummy_scalar
  pde(1)%pde_fnc(1)%zerord=>dummy_scalar
  pde(1)%pde_fnc(1)%der_convect=>dummy_vector
  allocate(elements%ders(1,3,2),elements%material(1),elements%areas(1))
  elements%ders(1,:,:)=g; elements%material=1; elements%areas=.5_rkind
  allocate(stiff_mat(3,3),cap_mat(3,3),bside(3),elnode_prev(3))
  allocate(base_fnc(3,3),gauss_points%weight(3))
  base_fnc(:,1)=[2.0_rkind/3,1.0_rkind/6,1.0_rkind/6]
  base_fnc(:,2)=[1.0_rkind/6,2.0_rkind/3,1.0_rkind/6]
  base_fnc(:,3)=[1.0_rkind/6,1.0_rkind/6,2.0_rkind/3]
  gauss_points%weight=1.0_rkind/6; gauss_points%area=.5_rkind
  drutes_config%dimen=2
  elnode_prev=1
  LSsupg=.true.; LSsupg_factor=1
  do method=0,2
    pde_common%timeint_method=method
    call build_bvect(1_ikind,dt)
    call build_stiff_np(1_ikind,dt)
    select case(method)
    case(0)
      call steady_state_int(1_ikind)
    case(1)
      call impl_euler_np_diag(1_ikind)
    case(2)
      call impl_euler_np_nondiag(1_ikind)
    end select
    before_cap=cap_mat; before_stiff=stiff_mat; before_rhs=bside
    ! Null hook leaves every matrix untouched (other models and SUPG-off).
    call apply_element_stabilization(1_ikind,dt)
    if (any(cap_mat/=before_cap) .or. any(stiff_mat/=before_stiff) .or. any(bside/=before_rhs)) &
      error stop 'null hook changes assembly'
    expected_cap=before_cap; expected_stiff=before_stiff; expected_rhs=before_rhs
    integral_mass=0
    tau=supg_tau(qtest/depthtest,tensor/depthtest,g,dt,method/=0,0.0_rkind)
    test=tau*matmul(g,qtest/depthtest)
    do l=1,3
      do i=1,3
        do j=1,3
          if (method/=0) integral_mass(i,j)=integral_mass(i,j) &
            -gauss_points%weight(l)*test(i)*depthtest*base_fnc(j,l)
          expected_stiff(i,j)=expected_stiff(i,j)-gauss_points%weight(l)*dt*test(i)*dot_product(qtest,g(j,:))
        end do
      end do
    end do
    expected_cap=expected_cap+integral_mass
    expected_rhs=expected_rhs+matmul(integral_mass,elnode_prev)
    pde(1)%stabilize_element=>adenc_supg_element
    call apply_element_stabilization(1_ikind,dt)
    if (maxval(abs(cap_mat-expected_cap))>1e-12_rkind) error stop 'capacity assembly'
    if (maxval(abs(stiff_mat-expected_stiff))>1e-12_rkind) error stop 'stiffness assembly'
    if (maxval(abs(bside-expected_rhs))>1e-12_rkind) error stop 'old-time RHS assembly'
    if (maxval(abs(matmul(stiff_mat+cap_mat,elnode_prev)-bside))>1e-12_rkind) &
      error stop 'constant preservation'
    nullify(pde(1)%stabilize_element)
  end do
  ! Exercise the real hook with lagged nodal values and SUPG disabled.
  allocate(elements%data(1,3),nodaltest(3))
  elements%data(1,:)=[1,2,3]; nodaltest=new
  pde(1)%getval=>current_value
  LSsupg=.false.; LSshock=.true.; LSshock_factor=1
  elnode_prev=old; cap_mat=0; stiff_mat=0; bside=0
  expected_stiff=0
  do l=1,3
    call shock_terms(qtest,depthtest,g,base_fnc(:,l),new,old,dt,.true.,0.0_rkind,0.0_rkind, &
                     1.0_rkind,spatial)
    expected_stiff=expected_stiff+gauss_points%weight(l)*spatial
  end do
  pde(1)%stabilize_element=>adenc_supg_element
  call apply_element_stabilization(1_ikind,dt)
  if (maxval(abs(stiff_mat-expected_stiff))>1e-12_rkind) error stop 'shock real assembly'
  if (any(cap_mat/=0) .or. any(bside/=0)) error stop 'shock changes mass or RHS'
  print *, 'ADEnc SUPG checks passed'
contains
  ! Small synthetic steady convection-diffusion solve input, not DRUtES main.
  ! Python solves this assembled dense matrix and checks bounds and exact error.
  subroutine strip_benchmark(enabled)
    logical, intent(in) :: enabled
    integer, parameter :: nx=20, nn=2*(nx+1), ne=2*nx
    integer :: e,k,i,j,unit,conn(3),ix
    real(kind=rkind) :: coordinates(nn,2),local_xy(3,2),det,matrix(nn,nn),load(nn)
    qtest=[1.0_rkind,0.0_rkind]; depthtest=1; diffusiontest=.005_rkind
    allocate(pde(1)); allocate(pde(1)%pde_fnc(1))
    LSshock=trim(mode)=='strip-shock'; LSshock_factor=1
    if (LSshock) then
      allocate(nodaltest(nn))
      open(newunit=unit,file='iterate.dat',status='old',action='read')
      read(unit,*) nodaltest
      close(unit)
      pde(1)%getval=>current_value
    end if
    pde(1)%pde_fnc(1)%convection=>convection
    pde(1)%pde_fnc(1)%dispersion=>dispersion
    pde(1)%pde_fnc(1)%elasticity=>storage
    pde(1)%pde_fnc(1)%reaction=>dummy_scalar
    pde(1)%pde_fnc(1)%zerord=>dummy_scalar
    pde(1)%pde_fnc(1)%der_convect=>dummy_vector
    if (enabled) pde(1)%stabilize_element=>adenc_supg_element
    allocate(elements%ders(ne,3,2),elements%material(ne),elements%areas(ne),elements%data(ne,3))
    do ix=0,nx
      coordinates(ix+1,:)=[real(ix,rkind)/nx,0.0_rkind]
      coordinates(nx+2+ix,:)=[real(ix,rkind)/nx,1.0_rkind]
    end do
    do ix=1,nx
      elements%data(2*ix-1,:)=[ix,ix+1,ix+nx+2]
      elements%data(2*ix,:)=[ix,ix+nx+2,ix+nx+1]
    end do
    do e=1,ne
      local_xy=coordinates(elements%data(e,:),:)
      det=(local_xy(2,1)-local_xy(1,1))*(local_xy(3,2)-local_xy(1,2)) &
         -(local_xy(3,1)-local_xy(1,1))*(local_xy(2,2)-local_xy(1,2))
      elements%areas(e)=abs(det)/2
      elements%ders(e,1,:)=[local_xy(2,2)-local_xy(3,2),local_xy(3,1)-local_xy(2,1)]/det
      elements%ders(e,2,:)=[local_xy(3,2)-local_xy(1,2),local_xy(1,1)-local_xy(3,1)]/det
      elements%ders(e,3,:)=[local_xy(1,2)-local_xy(2,2),local_xy(2,1)-local_xy(1,1)]/det
    end do
    elements%material=1
    allocate(stiff_mat(3,3),cap_mat(3,3),bside(3),elnode_prev(3))
    allocate(base_fnc(3,3),gauss_points%weight(3))
    base_fnc(:,1)=[2.0_rkind/3,1.0_rkind/6,1.0_rkind/6]
    base_fnc(:,2)=[1.0_rkind/6,2.0_rkind/3,1.0_rkind/6]
    base_fnc(:,3)=[1.0_rkind/6,1.0_rkind/6,2.0_rkind/3]
    gauss_points%weight=1.0_rkind/6; gauss_points%area=.5_rkind
    drutes_config%dimen=2; pde_common%timeint_method=0
    LSsupg=enabled; LSsupg_factor=1; elnode_prev=0
    matrix=0; load=0
    do e=1,ne
      call build_bvect(int(e,ikind),1.0_rkind)
      call build_stiff_np(int(e,ikind),1.0_rkind)
      call steady_state_int(int(e,ikind))
      call apply_element_stabilization(int(e,ikind),1.0_rkind)
      conn=int(elements%data(e,:))
      do i=1,3
        load(conn(i))=load(conn(i))+bside(i)
        do j=1,3
          matrix(conn(i),conn(j))=matrix(conn(i),conn(j))+stiff_mat(i,j)+cap_mat(i,j)
        end do
      end do
    end do
    ! Dirichlet: C(0)=0, C(1)=1. Top and bottom natural zero diffusive flux.
    do k=1,4
      select case(k)
      case(1); i=1
      case(2); i=nx+2
      case(3); i=nx+1
      case(4); i=nn
      end select
      matrix(i,:)=0; matrix(i,i)=1
      load(i)=0
      if (k>=3) load(i)=1
    end do
    open(newunit=unit,file='benchmark-matrix.dat',status='new',action='write')
    do i=1,nn
      write(unit,*) matrix(i,:)
    end do
    close(unit)
    open(newunit=unit,file='benchmark-rhs.dat',status='new',action='write')
    write(unit,*) load
    close(unit)
    print *, 'SUPG strip assembled'
  end subroutine strip_benchmark

  subroutine close_scalar(actual,expected,label)
    real(kind=rkind), intent(in) :: actual,expected
    character(len=*), intent(in) :: label
    if (abs(actual-expected)>1e-12_rkind) then
      print *, label,actual,expected
      error stop 'scalar comparison'
    end if
  end subroutine close_scalar
end program test_adenc_supg
