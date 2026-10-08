module bank_test_coefficients
  use typy
  use pde_objs
  use global_objs
  implicit none
contains
  subroutine convection(pde_loc,layer,quadpnt,x,vector_in,vector_out,scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:),vector_in(:)
    real(kind=rkind), intent(out), optional :: vector_out(:),scalar
    if (quadpnt%element==3) error stop 'inactive element assembled'
    if (present(vector_out)) vector_out=[.5_rkind,0.0_rkind]
    if (present(scalar)) scalar=.5_rkind
  end subroutine
  subroutine dispersion(pde_loc,layer,quadpnt,x,tensor,scalar)
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind), intent(in), optional :: x(:)
    real(kind=rkind), intent(out), optional :: tensor(:,:),scalar
    if (present(tensor)) then
      tensor=0; tensor(1,1)=.1_rkind; tensor(2,2)=.2_rkind
      tensor(1,2)=.02_rkind; tensor(2,1)=.02_rkind
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
  end function
end module

program test_adenc_banks
  use typy
  use globals
  use pde_objs
  use ncglobvars
  use ncboundary
  use ncfluxarea
  use ncpointers, only: adenc_element_corrections
  use dummy_procs
  use capmat
  use femmat
  use bank_test_coefficients
  implicit none
  integer(kind=ikind) :: original(6)
  integer :: i,j,k,mode,method,ierr
  logical :: ok
  character(len=1024) :: message
  character(len=32) :: argument
  real(kind=rkind) :: edge_matrix(3,3),expected(3,3),A(4,4),old(4),new(4),weights(4),mass0
  call get_command_argument(1,argument)
  if (len_trim(argument)>0) then
    drutes_config%dimen=2; drutes_config%it_method=0
    if (trim(argument)=='schwarz') drutes_config%it_method=1
    if (trim(argument)=='3d') drutes_config%dimen=3
    call read_adenc_banks()
    if (trim(argument)=='on' .and. .not. LSbank_noflow) error stop 'banks not enabled'
    if (trim(argument)=='absent' .and. LSbank_noflow) error stop 'legacy default changed'
    print *, 'bank configuration checks passed'
    stop
  end if
  call bank_edge_terms([1,2],2.0_rkind,[3.0_rkind,3.0_rkind],.5_rkind,edge_matrix)
  expected=0; expected(1,1)=1; expected(2,2)=1; expected(1,2)=.5_rkind; expected(2,1)=.5_rkind
  if (maxval(abs(edge_matrix-expected))>1.e-12_rkind) error stop 'Robin edge integration/sign'
  call bank_edge_terms([1,2],2.0_rkind,[0.0_rkind,0.0_rkind],.5_rkind,edge_matrix)
  if (any(edge_matrix/=0)) error stop 'tangential flow changes natural bank'
  allocate(pde(1))
  allocate(pde(1)%pde_fnc(1))
  nodes%kolik=6; elements%kolik=3
  allocate(nodes%data(6,3),nodes%edge(6),nodes%element(6),elements%data(3,3))
  allocate(elements%ders(3,3,2),elements%areas(3),elements%material(3))
  nodes%data=0
  nodes%data(:,1)=[0,1,1,0,2,2]+500000.0_rkind
  nodes%data(:,2)=[0,0,1,1,0,1]+5000000.0_rkind
  elements%data(1,:)=[1,2,3]; elements%data(2,:)=[1,3,4]; elements%data(3,:)=[2,3,5]
  ! Use production smartarray fill, including spare capacity. Invalid sentinels
  ! in the unused tail must never be interpreted as adjacent element IDs.
  call nodes%element(1)%fill(1_ikind); call nodes%element(1)%fill(2_ikind)
  call nodes%element(2)%fill(3_ikind); call nodes%element(2)%fill(1_ikind)
  call nodes%element(3)%fill(1_ikind); call nodes%element(3)%fill(2_ikind)
  call nodes%element(3)%fill(3_ikind)
  call nodes%element(4)%fill(2_ikind)
  call nodes%element(5)%fill(3_ikind); call nodes%element(5)%fill(3_ikind)
  call nodes%element(5)%fill(3_ikind)
  call nodes%element(6)%fill(3_ikind)
  nodes%element(3)%data(nodes%element(3)%pos+1:)=-999999_ikind
  nodes%element(5)%data(nodes%element(5)%pos+1:)=-999999_ikind
  elements%ders(1,:,1)=[-1,1,0]; elements%ders(1,:,2)=[0,-1,1]
  elements%ders(2,:,1)=[0,1,-1]; elements%ders(2,:,2)=[-1,0,1]
  elements%ders(3,:,:)=0; elements%areas=.5_rkind; elements%material=1
  allocate(ncfluxdata%activeel(3),ncfluxdata%fluxvct(3,2))
  ncfluxdata%activeel=[.true.,.true.,.false.]; addedbc=103; LSbank_noflow=.true.
  original=0; nodes%edge=addedbc
  call prepare_adenc_banks(original)
  if (count(bank_edges)/=1 .or. .not. bank_edges(2,1)) error stop 'internal bank/exterior distinction'
  call prepare_adenc_banks(original,.true.)
  if (count(bank_edges)/=4) error stop 'incorrect active-domain perimeter'
  if (any(nodes%edge(1:4)/=0) .or. any(nodes%edge(5:6)/=addedbc)) error stop 'bank DOFs/unused nodes'
  if (any(pde(1)%assembly_mask .neqv. ncfluxdata%activeel)) error stop 'assembly mask'
  if (active_node_element(2_ikind)/=1) error stop 'inactive-first nodal trace'
  if (active_node_element(5_ikind)/=3) error stop 'inactive-only nodal trace reads spare capacity'
  if (active_node_element(6_ikind)/=3) error stop 'single inactive adjacency fallback'
  ! Opposite inlet/outlet edges must be preserved; diagonal must not be a bank.
  original=[101,102,102,101,0,0]
  call prepare_adenc_banks(original,.true.)
  if (count(bank_edges)/=2) error stop 'port edges treated as banks'
  if (any(nodes%edge(1:4)/=original(1:4))) error stop 'physical port IDs changed'
  original=0; call prepare_adenc_banks(original,.true.)
  ! Synthetic cached hydrology, no file/model main invoked. Width=4, Q=2 => qx=.5.
  ncfluxdata%initialized=.true.; ncfluxdata%slice_loaded=.true.; ncfluxdata%current_time_index=1
  ncfluxdata%nlat=1; ncfluxdata%nlon=1; ncfluxdata%ntime=1; ncfluxdata%has_bounds=.true.
  allocate(ncfluxdata%time(1),ncfluxdata%lon(1),ncfluxdata%lat(1),ncfluxdata%qslice(1,1))
  allocate(ncfluxdata%lon_bnds(1,2),ncfluxdata%lat_bnds(1,2))
  ncfluxdata%time=0; ncfluxdata%lon=9; ncfluxdata%lat=45; ncfluxdata%qslice=2
  ncfluxdata%lon_bnds(1,:)=[-180.0_rkind,180.0_rkind]; ncfluxdata%lat_bnds(1,:)=[0.0_rkind,90.0_rkind]
  ncfluxdata%fluxvct(:,1)=1; ncfluxdata%fluxvct(:,2)=0
  allocate(ncnodes%data(4,2),ncelements%data(1,4),el2ncgrid(3))
  ncnodes%data(:,1)=[-1,3,3,-1]+500000.0_rkind; ncnodes%data(:,2)=[-1,-1,3,3]+5000000.0_rkind
  ncelements%kolik=1; ncelements%data(1,:)=[1,2,3,4]; el2ncgrid=1; Qmin=0
  call ncflux_prepare_widths(ok,message)
  if (.not. ok) error stop 'synthetic width preparation'
  if (abs(ncflux_active_width(1_ikind)-4)>1.e-12_rkind) error stop 'synthetic width'
  pde(1)%pde_fnc(1)%convection=>convection; pde(1)%pde_fnc(1)%dispersion=>dispersion
  pde(1)%pde_fnc(1)%elasticity=>storage; pde(1)%pde_fnc(1)%reaction=>dummy_scalar
  pde(1)%pde_fnc(1)%zerord=>dummy_scalar; pde(1)%pde_fnc(1)%der_convect=>dummy_vector
  pde(1)%getval=>getvalp1; pde(1)%stabilize_element=>adenc_element_corrections
  allocate(pde(1)%permut(6)); pde(1)%permut=[1,2,3,4,0,0]
  allocate(pde_common%xvect(4,4),pde_common%bvect(4))
  allocate(stiff_mat(3,3),cap_mat(3,3),bside(3),elnode_prev(3))
  allocate(base_fnc(3,3),gauss_points%weight(3))
  base_fnc(:,1)=[2.0_rkind/3,1.0_rkind/6,1.0_rkind/6]
  base_fnc(:,2)=[1.0_rkind/6,2.0_rkind/3,1.0_rkind/6]
  base_fnc(:,3)=[1.0_rkind/6,1.0_rkind/6,2.0_rkind/3]
  gauss_points%weight=1.0_rkind/6; gauss_points%area=.5_rkind
  drutes_config%dimen=2; drutes_config%it_method=0
  call spmatrix%init(4_ikind,4_ikind)
  ora_di_ini=0; time=0; time_step=.05_rkind
  weights=[1.0_rkind/3,1.0_rkind/6,1.0_rkind/3,1.0_rkind/6]
  ! Exercise real assemble_mat, bank hook, hydrological trace and both capacity choices.
  do method=1,2
    pde_common%timeint_method=method
    if (method==1) pde_common%time_integ=>impl_euler_np_diag
    if (method==2) pde_common%time_integ=>impl_euler_np_nondiag
    do mode=0,2
      LSsupg=mode>0; LSsupg_factor=2; LSshock=mode==2; LSshock_factor=1
      old=[.2_rkind,.8_rkind,.1_rkind,.4_rkind]; mass0=dot_product(weights,old)
      do k=1,30
        pde_common%xvect(:,1)=old; pde_common%xvect(:,2)=old
        call assemble_mat(ierr)
        do j=1,4
          do i=1,4
            A(i,j)=spmatrix%get(int(i,ikind),int(j,ikind))
          end do
        end do
        if (maxval(abs(sum(A,dim=1)+weights))>1.e-12_rkind) error stop 'closed-domain operator balance'
        call dense_solve(A,pde_common%bvect,new)
        if (abs(dot_product(weights,new)-mass0)>1.e-11_rkind) error stop 'closed-domain mass drift'
        old=new
      end do
    end do
  end do
  print *, 'ADEnc no-flow banks: topology, ports, inactive mask, Robin signs, closed mass checks passed'
contains
  subroutine dense_solve(matrix,rhs,x)
    real(kind=rkind), intent(in) :: matrix(4,4),rhs(4)
    real(kind=rkind), intent(out) :: x(4)
    real(kind=rkind) :: a(4,4),b(4),f
    integer :: i,j
    a=matrix; b=rhs
    do j=1,3
      do i=j+1,4
        f=a(i,j)/a(j,j); a(i,j:)=a(i,j:)-f*a(j,j:); b(i)=b(i)-f*b(j)
      end do
    end do
    do i=4,1,-1
      x(i)=(b(i)-dot_product(a(i,i+1:),x(i+1:)))/a(i,i)
    end do
  end subroutine
end program
