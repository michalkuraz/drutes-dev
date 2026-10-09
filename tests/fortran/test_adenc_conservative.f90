! Small production-assembly tests; never calls main or deletes simulation outputs.
module conservative_test_coefficients
  use typy
  use pde_objs
  use global_objs
  use ncglobvars, only: adenc_coefficient_time
  implicit none
  logical :: variable=.true.
contains
  subroutine convection(pde_loc,layer,quadpnt,x,vector_in,vector_out,scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:),vector_in(:)
    real(kind=rkind), intent(out), optional :: vector_out(:),scalar
    real(kind=rkind) :: q(2)
    q=[.5_rkind,0.0_rkind]
    if (variable) q=q*(1+.2_rkind*quadpnt%element+.1_rkind*adenc_coefficient_time())
    if (present(vector_out)) vector_out=q
    if (present(scalar)) scalar=norm2(q)
  end subroutine
  subroutine dispersion(pde_loc,layer,quadpnt,x,tensor,scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:)
    real(kind=rkind), intent(out), optional :: tensor(:,:),scalar
    if (present(tensor)) then
      tensor=0; tensor(1,1)=.1_rkind; tensor(2,2)=.02_rkind
    end if
    if (present(scalar)) scalar=.1_rkind
  end subroutine
  function storage(pde_loc,layer,quadpnt,x) result(value)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:)
    real(kind=rkind) :: value
    value=1
    if (variable) value=1+.3_rkind*quadpnt%element+.2_rkind*adenc_coefficient_time()
    if (variable .and. adenc_coefficient_time()>=.5_rkind) value=value+.4_rkind
  end function
end module

program test_adenc_conservative
  use typy
  use globals
  use pde_objs
  use global_objs
  use ncglobvars
  use ncconservative
  use ncbalance
  use ncboundary
  use ncfluxarea
  use ncpointers, only: adenc_element_corrections
  use lsconstitutive, only: ADEls_dirichlet
  use dummy_procs
  use capmat
  use femmat
  use conservative_test_coefficients
  implicit none
  integer(kind=ikind) :: original(6)
  integer :: i,j,k,mode,method,ierr
  logical :: ok
  character(len=1024) :: message
  character(len=32) :: argument
  real(kind=rkind) :: A(4,4),old(4),new(4),weights(4),mass0,edge_matrix(3,3),probe_dt
  call get_command_argument(1,argument)
  drutes_config%dimen=2; drutes_config%it_method=0; drutes_config%run_from_backup=.false.
  LSbank_noflow=.true.
  if (trim(argument)=='schwarz') drutes_config%it_method=1
  if (trim(argument)=='3d') drutes_config%dimen=3
  if (trim(argument)=='backup') drutes_config%run_from_backup=.true.
  if (len_trim(argument)>0 .and. trim(argument)/='inflow') then
    call read_adenc_conservative()
    if (trim(argument)=='absent' .and. (LSconservative .or. LSbalance)) error stop 'legacy default changed'
    if (trim(argument)=='off' .and. (LSconservative .or. LSbalance)) error stop 'off flags changed'
    if (trim(argument)=='on' .and. .not.(LSconservative .and. LSbalance)) error stop 'not enabled'
    if (trim(argument)=='audit' .and. (LSconservative .or. .not.LSbalance)) error stop 'audit flags'
    print *, 'conservative configuration checks passed'
    stop
  end if
  allocate(pde(1),pde_common%xvect(4,4),pde_common%bvect(4))
  allocate(pde(1)%pde_fnc(1),pde(1)%permut(6),pde(1)%bc(101:103))
  nodes%kolik=6; elements%kolik=3
  allocate(nodes%data(6,3),nodes%edge(6),nodes%element(6),elements%data(3,3))
  allocate(elements%ders(3,3,2),elements%areas(3),elements%material(3))
  nodes%data=0
  nodes%data(:,1)=[0,1,1,0,2,2]+500000.0_rkind
  nodes%data(:,2)=[0,0,1,1,0,1]+5000000.0_rkind
  elements%data(1,:)=[1,2,3]; elements%data(2,:)=[1,3,4]; elements%data(3,:)=[2,3,5]
  call nodes%element(1)%fill(1_ikind); call nodes%element(1)%fill(2_ikind)
  call nodes%element(2)%fill(1_ikind); call nodes%element(2)%fill(3_ikind)
  call nodes%element(3)%fill(1_ikind); call nodes%element(3)%fill(2_ikind)
  call nodes%element(3)%fill(3_ikind); call nodes%element(4)%fill(2_ikind)
  call nodes%element(5)%fill(3_ikind); call nodes%element(6)%fill(3_ikind)
  elements%ders(1,:,1)=[-1,1,0]; elements%ders(1,:,2)=[0,-1,1]
  elements%ders(2,:,1)=[0,1,-1]; elements%ders(2,:,2)=[-1,0,1]
  elements%ders(3,:,:)=0; elements%areas=.5_rkind; elements%material=1
  allocate(ncfluxdata%activeel(3),ncfluxdata%fluxvct(3,2))
  ncfluxdata%activeel=[.true.,.true.,.false.]; addedbc=103
  original=0; nodes%edge=addedbc
  call prepare_adenc_banks(original,.true.)
  ! Cached hydrology: Q=2, W=4, hence q=.5. No NetCDF file needed.
  ncfluxdata%initialized=.true.; ncfluxdata%slice_loaded=.true.; ncfluxdata%current_time_index=1
  ncfluxdata%nlat=1; ncfluxdata%nlon=1; ncfluxdata%ntime=1; ncfluxdata%has_bounds=.true.
  allocate(ncfluxdata%time(1),ncfluxdata%lon(1),ncfluxdata%lat(1),ncfluxdata%qslice(1,1))
  allocate(ncfluxdata%lon_bnds(1,2),ncfluxdata%lat_bnds(1,2))
  ncfluxdata%time=0; ncfluxdata%lon=9; ncfluxdata%lat=45; ncfluxdata%qslice=2
  ncfluxdata%lon_bnds(1,:)=[-180.0_rkind,180.0_rkind]; ncfluxdata%lat_bnds(1,:)=[0.0_rkind,90.0_rkind]
  ncfluxdata%fluxvct(:,1)=1; ncfluxdata%fluxvct(:,2)=0
  allocate(ncnodes%data(4,2),ncelements%data(1,4),el2ncgrid(3))
  ncnodes%data(:,1)=[-1,3,3,-1]+500000.0_rkind
  ncnodes%data(:,2)=[-1,-1,3,3]+5000000.0_rkind
  ncelements%kolik=1; ncelements%data(1,:)=[1,2,3,4]; el2ncgrid=1; Qmin=0
  call ncflux_prepare_widths(ok,message)
  if (.not.ok) error stop 'width preparation'
  pde(1)%pde_fnc(1)%convection=>convection; pde(1)%pde_fnc(1)%dispersion=>dispersion
  pde(1)%pde_fnc(1)%elasticity=>storage; pde(1)%pde_fnc(1)%reaction=>dummy_scalar
  pde(1)%pde_fnc(1)%zerord=>dummy_scalar; pde(1)%pde_fnc(1)%der_convect=>dummy_vector
  pde(1)%getval=>getvalp1; pde(1)%stabilize_element=>adenc_element_corrections
  pde(1)%step_begin=>conservative_begin; pde(1)%step_end=>conservative_end
  pde(1)%boundary_history=>adenc_boundary_history
  do i=101,103
    pde(1)%bc(i)%code=1; pde(1)%bc(i)%file=.false.; pde(1)%bc(i)%value=0
    pde(1)%bc(i)%value_fnc=>ADEls_dirichlet
  end do
  pde(1)%permut=[1,2,3,4,0,0]
  allocate(stiff_mat(3,3),cap_mat(3,3),bside(3),elnode_prev(3))
  allocate(base_fnc(3,3),gauss_points%weight(3))
  base_fnc(:,1)=[2.0_rkind/3,1.0_rkind/6,1.0_rkind/6]
  base_fnc(:,2)=[1.0_rkind/6,2.0_rkind/3,1.0_rkind/6]
  base_fnc(:,3)=[1.0_rkind/6,1.0_rkind/6,2.0_rkind/3]
  gauss_points%weight=1.0_rkind/6; gauss_points%area=.5_rkind
  call spmatrix%init(4_ikind,4_ikind)
  ora_di_ini=0
  do method=1,2
    pde_common%timeint_method=method
    if (method==1) pde_common%time_integ=>impl_euler_np_diag
    if (method==2) pde_common%time_integ=>impl_euler_np_nondiag
    do mode=0,2
      call reset_case()
      LSsupg=mode>0; LSsupg_factor=2; LSshock=mode==2; LSshock_factor=1
      old=[.2_rkind,.8_rkind,.1_rkind,.4_rkind]
      weights=[2.9_rkind/6,1.3_rkind/6,2.9_rkind/6,1.6_rkind/6]
      mass0=dot_product(weights,old)
      do k=1,30
        time_step=.05_rkind
        call step(old,new,4)
        if (abs(balance_mass-mass0)>1.e-11_rkind) error stop 'variable-H closed mass drift'
        if (abs(balance_error)>1.e-11_rkind) error stop 'closed balance error'
        if (balance_free_residual>1.e-11_rkind) error stop 'free row residual'
        old=new
      end do
    end do
  end do
  ! Audit-only leaves legacy operator unchanged and exposes its mass drift.
  call reset_case(); LSconservative=.false.; nullify(pde(1)%boundary_history)
  LSsupg=.false.; LSshock=.false.; old=[.2_rkind,.8_rkind,.1_rkind,.4_rkind]
  mass0=dot_product(weights,old); time_step=.05_rkind
  call step(old,new,4)
  if (abs(balance_mass-mass0)<1.e-6_rkind) error stop 'audit masked legacy drift'
  if (abs(balance_error)<1.e-6_rkind) error stop 'audit forced conservation'
  ! Actual external outflow: left Dirichlet inlet, top/bottom tangential.
  call reset_case(); variable=.false.; LSsupg=.false.; LSshock=.false.
  pde(1)%boundary_history=>adenc_boundary_history
  elements%data(3,:)=[5,6,5] ! remove inactive neighbour on right boundary
  original=[101,0,0,101,0,0]
  if (trim(argument)=='inflow') original=0
  call prepare_adenc_banks(original)
  if (trim(argument)=='inflow') call adenc_open_element(2_ikind,.1_rkind,edge_matrix)
  call adenc_open_element(1_ikind,.1_rkind,edge_matrix)
  if (abs(sum(edge_matrix)+.05_rkind)>1.e-12_rkind) error stop 'exterior flux sign'
  if (trim(argument)=='inflow') error stop 'unspecified inflow did not fail'
  pde(1)%permut=[0,1,2,0,0,0]
  pde(1)%bc(101)%file=.true.; allocate(pde(1)%bc(101)%series(2,2))
  pde(1)%bc(101)%series(1,:)=[0.0_rkind,1.0_rkind]
  pde(1)%bc(101)%series(2,:)=[.2_rkind,0.0_rkind]
  old=0
  do k=1,20
    time_step=.1_rkind
    call step(old,new,2)
    if (abs(balance_error)>1.e-11_rkind) error stop 'pulse/outflow inventory budget'
    if (k==2 .and. pde(1)%boundary_history(2_ikind,1_ikind,1_ikind)/=1) &
      error stop 'pulse lost interval ending at discontinuity'
    if (k==3 .and. pde(1)%boundary_history(2_ikind,1_ikind,4_ikind)/=1) &
      error stop 'old Dirichlet trace overwritten'
    old=new
  end do
  if (maxval(new(1:2))<=0) error stop 'pulse did not propagate'
  ! Rejection changes neither accepted storage nor clock, and event clipping.
  call reset_case(); probe_dt=1; time=.1_rkind
  call pde(1)%step_begin(time,probe_dt)
  if (abs(probe_dt-.1_rkind)>1.e-14_rkind) error stop 'pulse clipping'
  LSdepth_new=99
  call pde(1)%step_end(.false.)
  if (LSstate_time/=0 .or. any(LSdepth_old/=1)) error stop 'rejected state committed'
  time=86399; probe_dt=10
  call pde(1)%step_begin(time,probe_dt)
  if (abs(probe_dt-1)>1.e-14_rkind) error stop 'daily forcing clipping'
  call pde(1)%step_end(.false.)
  call balance_reset()
  print *, 'ADEnc conservative checks passed: storage, modes, pulse, outflow, audit and rejection'
contains
  subroutine reset_case()
    call balance_reset()
    if (allocated(LSdepth_old)) deallocate(LSdepth_old,LSdepth_new)
    LSstate_time=0; LSprevious_time=0; LSstep_active=.false.; LSclock_override=.false.
    LSconservative=.true.; LSbalance=.true.; time=0
  end subroutine
  subroutine step(previous,next,n)
    real(kind=rkind), intent(in) :: previous(4)
    real(kind=rkind), intent(out) :: next(4)
    integer, intent(in) :: n
    pde_common%xvect(:,1)=previous; pde_common%xvect(:,2)=previous
    call pde(1)%step_begin(time,time_step)
    call assemble_mat(ierr)
    do j=1,n
      do i=1,n
        A(i,j)=spmatrix%get(int(i,ikind),int(j,ikind))
      end do
    end do
    next=0
    call dense_solve(A(:n,:n),pde_common%bvect(:n),next(:n))
    pde_common%xvect(:,3)=next
    call pde(1)%step_end(.true.)
    time=time+time_step
  end subroutine
  subroutine dense_solve(matrix,rhs,x)
    real(kind=rkind), intent(in) :: matrix(:,:),rhs(:)
    real(kind=rkind), intent(out) :: x(:)
    real(kind=rkind) :: a(size(x),size(x)),b(size(x)),f
    integer :: i,j,n
    n=size(x); a=matrix; b=rhs
    do j=1,n-1
      do i=j+1,n
        f=a(i,j)/a(j,j); a(i,j:)=a(i,j:)-f*a(j,j:); b(i)=b(i)-f*b(j)
      end do
    end do
    do i=n,1,-1
      x(i)=(b(i)-dot_product(a(i,i+1:),x(i+1:)))/a(i,i)
    end do
  end subroutine
end program
