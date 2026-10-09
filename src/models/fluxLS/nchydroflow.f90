! Optional compatible hydraulic reconstruction for conservative ADEnc.
! One oriented volume flux per FE edge, RT0 inside each triangle, P0 storage.
! This is a projection of supplied hydrology, NOT a hydrodynamic flow model.
module nchydroflow
  use typy
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: LShydro,read_adenc_hydro,hydro_filter,hydro_initialize,hydro_begin,hydro_end
  public :: hydro_value,hydro_normal_flux,hydro_project,hydro_rt0,hydro_velocity
  public :: hydro_project_outflow
  public :: hydro_sources,hydro_clip_step
  logical :: LShydro=.false.,ready=.false.,trial=.false.
  real(kind=rkind) :: tolerance=1.e-10_rkind,correction_limit=1.0_rkind
  integer :: max_iterations=20000,forcing_mode=1
  integer, allocatable :: ports(:,:),owners(:,:),edge_nodes(:,:),edge_of(:,:),orientation(:,:)
  integer, allocatable :: port_edge(:),port_kind(:),port_id(:)
  real(kind=rkind), allocatable :: port_flow(:),area(:),center(:,:),edge_length(:),normal(:,:)
  real(kind=rkind), allocatable :: depth_state(:),depth_trial(:),flux_state(:),flux_trial(:)
  real(kind=rkind), allocatable :: raw_depth(:,:),raw_q(:,:,:),cached_target(:),cached_flux(:)
  integer :: raw_hour=-huge(1),projected_hour=-huge(1)
  real(kind=rkind) :: raw_time=-huge(1.0_rkind),projected_time=-huge(1.0_rkind)
  logical :: projection_cached=.false.
  ! Explicit unresolved tributary/lateral input; never inferred from Q differences.
  type :: lateral_group
    integer, allocatable :: elements(:)
    real(kind=rkind), allocatable :: time(:),water(:),solute(:)
    real(kind=rkind) :: volume_area=0
  end type
  type(lateral_group), allocatable :: lateral(:)
  real(kind=rkind), allocatable :: water_state(:),water_trial(:),solute_state(:),solute_trial(:)
  integer, parameter :: ends(2,3)=reshape([1,2,2,3,3,1],[2,3]),opposite(3)=[3,1,2]
contains
  subroutine read_adenc_hydro()
    use ncglobvars, only: LSconservative,LSbank_noflow
    use readtools, only: fileread
    use globals, only: drutes_config
    integer :: unit,i,n,policy,status,a,b,kind,extra
    logical :: exists
    character(len=1024) :: line
    real(kind=rkind) :: flow
    LShydro=.false.; ready=.false.; trial=.false.
    if (allocated(lateral)) deallocate(lateral)
    if (allocated(ports)) deallocate(ports,port_kind,port_flow)
    inquire(file='drutes.conf/netcdf/hydroflow.conf',exist=exists)
    if (.not.exists) return
    open(newunit=unit,file='drutes.conf/netcdf/hydroflow.conf',status='old',action='read')
    call fileread(LShydro,unit)
    if (.not.LShydro) then
      close(unit)
      return
    end if
    if (.not.LSconservative .or. .not.LSbank_noflow .or. drutes_config%dimen/=2 .or. &
        drutes_config%it_method/=0) error stop 'Hydro reconstruction requires conservative 2D Picard/no-flow banks'
    call fileread(policy,unit)
    if (policy/=0 .and. policy/=1) error stop 'Hydroflow lateral source policy must be 0 or 1'
    call fileread(forcing_mode,unit)
    if (forcing_mode/=0 .and. forcing_mode/=1) error stop 'Hydroflow forcing mode must be 0 (daily) or 1 (linear)'
    call fileread(tolerance,unit)
    call fileread(max_iterations,unit)
    call fileread(correction_limit,unit)
    if (.not.ieee_is_finite(tolerance) .or. .not.ieee_is_finite(correction_limit)) &
      error stop 'Nonfinite hydroflow tolerance/correction limit'
    if (tolerance<=0 .or. tolerance>1.e-6_rkind .or. max_iterations<1 .or. correction_limit<=0) &
      error stop 'Invalid hydroflow tolerance/iterations/correction limit'
    call fileread(n,unit)
    if (n<2) error stop 'Hydroflow needs explicit inlet and outlet edges'
    allocate(ports(2,n),port_kind(n),port_flow(n))
    do i=1,n
      call record(unit,line,status)
      if (status/=0) error stop 'Missing hydroflow port record'
      if (scan(line,',/*')>0) error stop 'Invalid hydroflow port list syntax'
      read(line,*,iostat=status) a,b,kind,flow
      if (status/=0) error stop 'Hydroflow ports require node1 node2 kind discharge'
      read(line,*,iostat=status) a,b,kind,flow,extra
      if (status>=0) error stop 'Extra or invalid hydroflow port field'
      if (.not.ieee_is_finite(flow)) error stop 'Nonfinite hydroflow port discharge'
      if (a<1 .or. b<1 .or. a==b .or. (kind/=-1 .and. kind/=1)) error stop 'Invalid hydroflow port'
      if (flow<0 .or. (kind==1 .and. flow/=0)) &
        error stop 'Inlet uses positive m3/s or zero for Qrouted; outlet uses zero (free balanced discharge)'
      ports(:,i)=[min(a,b),max(a,b)]; port_kind(i)=kind; port_flow(i)=flow
    end do
    call record(unit,line,status)
    if (status==0) error stop 'Unexpected trailing hydroflow record'
    close(unit)
    if (policy==1) call read_lateral()
  end subroutine

  ! Each group has an independently supplied TOTAL discharge and concentration.
  ! Distribution over listed FE triangles is uniform per area, not Q per triangle.
  subroutine read_lateral()
    use globals, only: end_time
    integer :: unit,status,n,g,ne,nt,i,j,extra
    real(kind=rkind) :: t,q,c,unused
    character(len=1024) :: line
    open(newunit=unit,file='drutes.conf/netcdf/lateral.conf',status='old',action='read',iostat=status)
    if (status/=0) error stop 'Hydroflow policy 1 requires lateral.conf'
    if (.not.ieee_is_finite(end_time) .or. end_time<=0) error stop 'Invalid lateral-source simulation period'
    call record(unit,line,status)
    if (status/=0) error stop 'Missing lateral group count'
    read(line,*,iostat=status) n
    if (status/=0) error stop 'Invalid lateral group count'
    read(line,*,iostat=status) n,extra
    if (status>=0 .or. n<1 .or. scan(line,',/*')>0) error stop 'Invalid lateral group count'
    allocate(lateral(n))
    do g=1,n
      call record(unit,line,status)
      if (status/=0) error stop 'Missing lateral group sizes'
      read(line,*,iostat=status) ne,nt
      if (status/=0) error stop 'Invalid lateral group sizes'
      read(line,*,iostat=status) ne,nt,extra
      if (status>=0 .or. ne<1 .or. nt<2 .or. scan(line,',/*')>0) error stop 'Invalid lateral group sizes'
      allocate(lateral(g)%elements(ne),lateral(g)%time(nt),lateral(g)%water(nt),lateral(g)%solute(nt))
      do i=1,ne
        call record(unit,line,status)
        if (status/=0) error stop 'Missing lateral element'
        read(line,*,iostat=status) j
        if (status/=0) error stop 'Invalid lateral element'
        read(line,*,iostat=status) j,extra
        if (status>=0 .or. j<1 .or. scan(line,',/*')>0) error stop 'Invalid lateral element'
        if (any(lateral(g)%elements(:i-1)==j)) error stop 'Duplicate lateral element in one group'
        lateral(g)%elements(i)=j
      end do
      do i=1,nt
        call record(unit,line,status)
        if (status/=0) error stop 'Missing lateral time Q C record'
        read(line,*,iostat=status) t,q,c
        if (status/=0) error stop 'Invalid lateral time Q C record'
        read(line,*,iostat=status) t,q,c,unused
        if (status>=0 .or. scan(line,',/*')>0) error stop 'Extra or invalid lateral time Q C field'
        if (.not.all(ieee_is_finite([t,q,c]))) error stop 'Nonfinite lateral time Q C'
        if (t<0 .or. q<0 .or. c<0) error stop 'Lateral inputs require nonnegative time Q C (no withdrawals)'
        if (i==1) then
          if (t/=0) error stop 'Lateral series must start at simulation time zero'
        else
          if (t<=lateral(g)%time(i-1)) error stop 'Lateral times must strictly increase'
        end if
        lateral(g)%time(i)=t; lateral(g)%water(i)=q; lateral(g)%solute(i)=q*c
        if (.not.ieee_is_finite(q*c)) error stop 'Nonfinite lateral solute load'
      end do
      if (lateral(g)%time(nt)<end_time) error stop 'Lateral series does not cover configured simulation period'
    end do
    call record(unit,line,status)
    if (status==0) error stop 'Unexpected trailing lateral record'
    close(unit)
  end subroutine

  ! Knots are explicit events; no step crosses a change of interpolation segment.
  subroutine hydro_clip_step(t,dt)
    real(kind=rkind), intent(in) :: t
    real(kind=rkind), intent(in out) :: dt
    integer :: g,i
    if (.not.LShydro .or. .not.allocated(lateral)) return
    do g=1,size(lateral)
      if (t>=lateral(g)%time(size(lateral(g)%time))) error stop 'Lateral series does not cover trial interval'
      do i=2,size(lateral(g)%time)
        if (lateral(g)%time(i)>t) then
          dt=min(dt,lateral(g)%time(i)-t)
          exit
        end if
      end do
    end do
  end subroutine

  subroutine lateral_sample(t,water,solute)
    real(kind=rkind), intent(in) :: t
    real(kind=rkind), intent(out) :: water(:),solute(:)
    integer :: g,i,nt
    real(kind=rkind) :: fraction,q,load
    water=0; solute=0
    if (.not.allocated(lateral)) return
    do g=1,size(lateral)
      nt=size(lateral(g)%time)
      if (.not.ieee_is_finite(t)) error stop 'Invalid lateral forcing time'
      if (t<0 .or. t>lateral(g)%time(nt)) error stop 'Lateral series does not cover trial interval'
      i=1
      do while(i<nt-1)
        if (lateral(g)%time(i+1)>=t) exit
        i=i+1
      end do
      fraction=(t-lateral(g)%time(i))/(lateral(g)%time(i+1)-lateral(g)%time(i))
      q=(1-fraction)*lateral(g)%water(i)+fraction*lateral(g)%water(i+1)
      load=(1-fraction)*lateral(g)%solute(i)+fraction*lateral(g)%solute(i+1)
      water(lateral(g)%elements)=water(lateral(g)%elements)+q/lateral(g)%volume_area
      solute(lateral(g)%elements)=solute(lateral(g)%elements)+load/lateral(g)%volume_area
    end do
    if (.not.all(ieee_is_finite(water)) .or. .not.all(ieee_is_finite(solute))) &
      error stop 'Nonfinite distributed lateral source'
  end subroutine

  ! Depth-rate water source [m/s], solute source [concentration*m/s].
  ! Trial values are interval averages, shared by hydraulic balance and transport.
  subroutine hydro_sources(el,water,solute)
    integer, intent(in) :: el
    real(kind=rkind), intent(out) :: water,solute
    water=0; solute=0
    if (.not.LShydro) return
    if (.not.ready) error stop 'Lateral source before hydraulic initialization'
    if (trial) then
      water=water_trial(el); solute=solute_trial(el)
    else
      water=water_state(el); solute=solute_state(el)
    end if
  end subroutine

  subroutine record(unit,line,status)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    integer, intent(out) :: status
    integer :: comment
    do
      read(unit,'(A)',iostat=status) line
      if (status/=0) return
      comment=index(line,'#')
      if (comment>0) line=line(:comment-1)
      if (len_trim(line)>0) return
    end do
  end subroutine

  ! Preserve the existing empirical velocity law, including its width cancellation.
  function hydro_velocity(Q,W) result(v)
    use ncglobvars, only: vref,Qref
    real(kind=rkind), intent(in) :: Q,W
    real(kind=rkind) :: v,k,m
    m=3.0_rkind/5.0_rkind
    if (Q<=0 .or. W<=0 .or. Qref<=0 .or. vref<=0) error stop 'Invalid hydraulic velocity inputs'
    k=vref*W**(1-m)/Qref**(1-m)
    v=k*Q**(1-m)/W**(1-m)
  end function

  ! Opt-in P0 homogenization: ONE centroid hydro cell defines geometry, Q and H.
  ! Not exact cut-cell geometry. The approximation and removed elements are logged.
  subroutine hydro_filter()
    use globals, only: elements,nodes
    use ncglobvars, only: ncfluxdata,Qmin,ora_di_ini
    use ncfluxarea, only: ncflux_active_width
    use netcdfflux, only: ncflux_get_xy_cell
    use core_tools, only: write_log
    integer :: e,removed
    real(kind=rkind) :: xy(2),Q
    logical :: ok
    character(len=1024) :: message
    if (.not.LShydro) return
    raw_hour=-huge(1); projected_hour=-huge(1); projection_cached=.false.
    if (allocated(raw_depth)) deallocate(raw_depth,raw_q,cached_target,cached_flux)
    removed=0
    do e=1,elements%kolik
      if (.not.ncfluxdata%activeel(e)) cycle
      xy=sum(nodes%data(elements%data(e,:),1:2),dim=1)/3
      call ncflux_get_xy_cell(xy(1),xy(2),ora_di_ini,Q,ok,message)
      if (.not.ok .or. Q<Qmin .or. Q<=0 .or. ncflux_active_width(int(e,ikind))<=0) then
        ncfluxdata%activeel(e)=.false.; removed=removed+1
      end if
    end do
    write(message,*) 'ADEnc hydroflow P0 centroid domain; removed invalid/zero-width elements=',removed
    call write_log(trim(message))
  end subroutine

  subroutine hydro_initialize(original_edge)
    use globals, only: elements,nodes
    use ncglobvars, only: ncfluxdata,bank_edges,open_edges,addedbc
    integer(kind=ikind), intent(in) :: original_edge(:)
    integer, allocatable :: head(:),next(:),lo(:),hi(:),root(:)
    logical, allocatable :: inlet_node(:),has_in(:),has_out(:)
    integer :: e,k,a,b,bucket,item,nedge,buckets,p,n,r,e2
    real(kind=rkind) :: aa(2),bb(2),cc(2),side(2),qraw(elements%kolik,2),rate(elements%kolik)
    if (.not.LShydro) return
    raw_hour=-huge(1); projected_hour=-huge(1); projection_cached=.false.
    if (allocated(raw_depth)) deallocate(raw_depth,raw_q,cached_target,cached_flux)
    if (allocated(owners)) then
      deallocate(owners,edge_nodes,edge_of,orientation,port_edge,area,center,edge_length,normal, &
        depth_state,depth_trial,flux_state,flux_trial,port_id)
      deallocate(water_state,water_trial,solute_state,solute_trial)
    end if
    n=elements%kolik; buckets=2*n+1
    allocate(head(buckets),next(3*n),lo(3*n),hi(3*n),owners(3*n,2),edge_nodes(3*n,2))
    allocate(edge_of(3,n),orientation(3,n),area(n),center(n,2),edge_length(3*n),normal(3*n,2))
    head=0; owners=0; edge_of=0; orientation=0; area=0; nedge=0
    do e=1,n
      if (.not.ncfluxdata%activeel(e)) cycle
      aa=nodes%data(elements%data(e,1),1:2); bb=nodes%data(elements%data(e,2),1:2)
      cc=nodes%data(elements%data(e,3),1:2)
      side=bb-aa
      area(e)=abs(side(1)*(cc(2)-aa(2))-side(2)*(cc(1)-aa(1)))/2
      if (area(e)<=0) error stop 'Degenerate hydroflow triangle'
      center(e,:)=(aa+bb+cc)/3
      do k=1,3
        a=minval(elements%data(e,ends(:,k))); b=maxval(elements%data(e,ends(:,k)))
        bucket=int(modulo(31_8*int(a,8)+int(b,8),int(buckets,8)))+1
        item=head(bucket)
        do while(item/=0)
          if (lo(item)==a .and. hi(item)==b) exit
          item=next(item)
        end do
        if (item==0) then
          nedge=nedge+1; item=nedge
          lo(item)=a; hi(item)=b; next(item)=head(bucket); head(bucket)=item
          owners(item,1)=e; edge_nodes(item,:)=[a,b]; orientation(k,e)=1
          aa=nodes%data(a,1:2); bb=nodes%data(b,1:2)
          edge_length(item)=norm2(bb-aa)
          normal(item,:)=[bb(2)-aa(2),aa(1)-bb(1)]/edge_length(item)
          if (dot_product(normal(item,:),center(e,:)-(aa+bb)/2)>0) normal(item,:)=-normal(item,:)
        else
          if (owners(item,2)/=0) error stop 'Nonmanifold hydroflow edge'
          owners(item,2)=e; orientation(k,e)=-1
        end if
        edge_of(k,e)=item
      end do
    end do
    if (nedge==0) error stop 'Empty hydroflow domain'
    owners=owners(:nedge,:); edge_nodes=edge_nodes(:nedge,:)
    edge_length=edge_length(:nedge); normal=normal(:nedge,:)
    allocate(port_edge(size(port_kind)),port_id(size(port_kind)), &
      inlet_node(nodes%kolik),root(n),has_in(n),has_out(n))
    inlet_node=.false.; has_in=.false.; has_out=.false.; port_edge=0; port_id=0
    do e=1,n
      root(e)=e
    end do
    do item=1,nedge
      e2=owners(item,2)
      if (e2==0) cycle
      a=find_root(root,owners(item,1)); b=find_root(root,e2)
      root(b)=a
    end do
    do p=1,size(port_kind)
      do item=1,nedge
        if (all(edge_nodes(item,:)==ports(:,p))) exit
      end do
      if (item>nedge) error stop 'Hydroflow port is not an active FE edge'
      if (owners(item,2)/=0) error stop 'Hydroflow port must be on active-domain boundary'
      if (any(port_edge==item)) error stop 'Duplicate hydroflow port'
      port_edge(p)=item; r=find_root(root,owners(item,1))
      if (port_kind(p)==-1) then
        a=ports(1,p); b=ports(2,p)
        if (original_edge(a)<=100 .or. original_edge(a)/=original_edge(b) .or. &
            original_edge(a)==addedbc) error stop 'Hydroflow inlet must have a real concentration boundary ID'
        inlet_node([a,b])=.true.; has_in(r)=.true.
        port_id(p)=original_edge(a)
      else
        has_out(r)=.true.
      end if
    end do
    do p=1,size(port_kind)
      if (port_kind(p)/=-1) cycle
      do k=1,size(port_kind)
        if (port_kind(k)/=-1 .or. port_id(p)/=port_id(k)) cycle
        if ((port_flow(p)==0) .neqv. (port_flow(k)==0)) &
          error stop 'Do not mix fixed and Qrouted inlet edges under one boundary ID'
      end do
    end do
    do e=1,n
      if (area(e)<=0) cycle
      r=find_root(root,e)
      if (.not.has_in(r) .or. .not.has_out(r)) &
        error stop 'Every hydroflow component needs an explicit inlet and outlet'
    end do
    ! Every unlabelled exposed edge is a WATER bank, even on the exterior mesh.
    ! Only listed inlet nodes retain Dirichlet DOFs; outlet and bank nodes are free.
    do e=1,n
      if (area(e)<=0) cycle
      do k=1,3
        item=edge_of(k,e)
        bank_edges(k,e)=owners(item,2)==0; open_edges(k,e)=.false.
        do p=1,size(port_kind)
          if (port_edge(p)/=item) cycle
          bank_edges(k,e)=.false.; open_edges(k,e)=port_kind(p)==1
        end do
      end do
      do k=1,3
        a=elements%data(e,k)
        nodes%edge(a)=0
        if (inlet_node(a)) nodes%edge(a)=original_edge(a)
      end do
    end do
    do p=1,size(port_kind)
      if (port_kind(p)==1 .and. any(inlet_node(ports(:,p)))) &
        error stop 'Hydroflow inlet and outlet cannot share a node'
    end do
    allocate(depth_state(n),depth_trial(n),flux_state(nedge),flux_trial(nedge))
    allocate(raw_depth(n,2),raw_q(n,2,2),cached_target(n),cached_flux(nedge))
    allocate(water_state(n),water_trial(n),solute_state(n),solute_trial(n))
    if (allocated(lateral)) then
      do p=1,size(lateral)
        if (any(lateral(p)%elements>n)) error stop 'Lateral element outside FE mesh'
        if (any(area(lateral(p)%elements)<=0)) error stop 'Lateral source element is inactive'
        lateral(p)%volume_area=sum(area(lateral(p)%elements))
      end do
    end if
    call lateral_sample(0.0_rkind,water_state,solute_state)
    call raw_field(0.0_rkind,depth_state,qraw,rate)
    depth_trial=depth_state
    call project_field(qraw,area*(water_state-rate),flux_state)
    water_trial=water_state; solute_trial=solute_state
    flux_trial=flux_state; trial=.false.; ready=.true.
  end subroutine

  integer function find_root(root,i) result(r)
    integer, intent(in out) :: root(:)
    integer, intent(in) :: i
    integer :: a,b
    r=i
    do while(root(r)/=r)
      r=root(r)
    end do
    a=i
    do while(root(a)/=a)
      b=root(a); root(a)=r; a=b
    end do
  end function

  subroutine raw_field(t,depth,qraw,rate)
    use ncglobvars, only: ora_di_ini
    real(kind=rkind), intent(in) :: t
    real(kind=rkind), intent(out) :: depth(:),qraw(:,:)
    real(kind=rkind), intent(out), optional :: rate(:)
    integer :: hour,day
    real(kind=rkind) :: fraction
    if (.not.ieee_is_finite(t) .or. t<0) error stop 'Invalid hydroflow forcing time'
    day=floor(t/86400.0_rkind); hour=int(ora_di_ini)+day*24
    if (raw_hour/=hour) then
      call raw_sample(hour,raw_depth(:,1),raw_q(:,:,1))
      raw_depth(:,2)=raw_depth(:,1); raw_q(:,:,2)=raw_q(:,:,1)
      if (forcing_mode==1) call raw_sample(hour+24,raw_depth(:,2),raw_q(:,:,2))
      raw_hour=hour
    end if
    fraction=0
    if (forcing_mode==1) fraction=(t-day*86400.0_rkind)/86400.0_rkind
    depth=(1-fraction)*raw_depth(:,1)+fraction*raw_depth(:,2)
    qraw=(1-fraction)*raw_q(:,:,1)+fraction*raw_q(:,:,2)
    if (present(rate)) rate=(raw_depth(:,2)-raw_depth(:,1))/86400.0_rkind
    raw_time=t
  end subroutine

  subroutine raw_sample(hour,depth,qraw)
    use ncglobvars, only: ncfluxdata
    use netcdfflux, only: ncflux_get_xy_cell
    use ncfluxarea, only: ncflux_active_width
    integer, intent(in) :: hour
    real(kind=rkind), intent(out) :: depth(:),qraw(:,:)
    integer :: e
    real(kind=rkind) :: Q,W
    logical :: ok
    character(len=1024) :: message
    depth=1; qraw=0
    do e=1,size(depth)
      if (area(e)<=0) cycle
      call ncflux_get_xy_cell(center(e,1),center(e,2),int(hour,ikind),Q,ok,message)
      if (.not.ok .or. .not.ieee_is_finite(Q) .or. Q<=0) &
        error stop 'Hydroflow fixed domain became dry/missing: wetting/drying is not implemented'
      W=ncflux_active_width(int(e,ikind))
      if (W<=0) error stop 'Hydroflow width vanished'
      if (abs(norm2(ncfluxdata%fluxvct(e,:))-1)>1.e-10_rkind) &
        error stop 'Hydroflow requires unit channel directions'
      qraw(e,:)=Q/W*ncfluxdata%fluxvct(e,:)
      depth(e)=Q/(W*hydro_velocity(Q,W))
      if (.not.ieee_is_finite(depth(e)) .or. depth(e)<=0) error stop 'Invalid hydroflow storage'
    end do
  end subroutine

  subroutine hydro_begin(t,dt,start_time)
    real(kind=rkind), intent(in) :: t,dt
    real(kind=rkind), intent(in), optional :: start_time
    real(kind=rkind), allocatable :: qraw(:,:),target(:),water_old(:),solute_old(:)
    real(kind=rkind) :: start,span
    if (.not.LShydro) return
    if (.not.ready .or. dt<=0) error stop 'Hydroflow step before initialization/invalid dt'
    allocate(qraw(size(area),2),target(size(area)),water_old(size(area)),solute_old(size(area)))
    start=max(0.0_rkind,t-dt)
    ! The forcing endpoint is a left-limit timestamp; subtracting dt can put
    ! a knot-aligned start just BEFORE its knot. Use the actual solver start.
    if (present(start_time)) start=start_time
    ! Caller clips at knots. Exact interval mean of linear Q and Q*C.
    span=dt
    call hydro_clip_step(start,span)
    if (dt-span>32*epsilon(dt)*max(1.0_rkind,abs(t),dt)) &
      error stop 'Hydraulic trial crosses a lateral forcing knot'
    call lateral_sample(start,water_old,solute_old)
    call lateral_sample(t,water_trial,solute_trial)
    water_trial=(water_old+water_trial)/2; solute_trial=(solute_old+solute_trial)/2
    call raw_field(t,depth_trial,qraw)
    target=area*(water_trial-(depth_trial-depth_state)/dt)
    call project_field(qraw,target,flux_trial)
    trial=.true.
  end subroutine

  subroutine hydro_end(accepted)
    logical, intent(in) :: accepted
    if (.not.LShydro) return
    if (accepted) then
      depth_state=depth_trial; flux_state=flux_trial
      water_state=water_trial; solute_state=solute_trial
    end if
    trial=.false.
  end subroutine

  subroutine project_field(qraw,target,flux)
    use core_tools, only: write_log
    real(kind=rkind), intent(in) :: qraw(:,:),target(:)
    real(kind=rkind), intent(out) :: flux(:)
    real(kind=rkind) :: preferred(size(flux)),fixed_value(size(flux)),weight(size(flux)),q(2),error,change,total_length
    logical :: fixed(size(flux)),outflow(size(flux)),ok
    integer :: i,e,j,p,k,iterations,bounded,passes
    character(len=256) :: message
    if (projection_cached .and. projected_hour==raw_hour .and. &
        (forcing_mode==0 .or. projected_time==raw_time)) then
      if (all(target==cached_target)) then
        flux=cached_flux
        return
      end if
    end if
    outflow=.false.
    do i=1,size(flux)
      e=owners(i,1); j=owners(i,2); q=qraw(e,:)
      if (j>0) q=(q+qraw(j,:))/2
      preferred(i)=edge_length(i)*dot_product(q,normal(i,:))
      weight(i)=edge_length(i)**2/area(e)
      if (j>0) weight(i)=edge_length(i)**2/(area(e)+area(j))
      fixed(i)=j==0; fixed_value(i)=0
    end do
    do p=1,size(port_kind)
      i=port_edge(p)
      if (port_kind(p)==-1) then
        fixed_value(i)=-port_flow(p)
        if (port_flow(p)==0) then
          ! A real inlet ID represents one port. Split its Qrouted discharge
          ! across all listed edges, rather than injecting the full Q per edge.
          total_length=0
          do k=1,size(port_kind)
            if (port_kind(k)==-1 .and. port_id(k)==port_id(p)) &
              total_length=total_length+edge_length(port_edge(k))
          end do
          e=owners(i,1)
          fixed_value(i)=-norm2(qraw(e,:))*hydro_width(e)*edge_length(i)/total_length
        end if
      else
        fixed(i)=.false.
        outflow(i)=.true.
      end if
    end do
    call hydro_project_outflow(owners,preferred,weight,fixed,fixed_value,outflow,target,tolerance,max_iterations, &
      flux,ok,iterations,error,bounded,passes)
    if (.not.ok) error stop 'Hydroflow outflow-constrained solve failed: check water balance/ports/solver'
    change=norm2(flux-preferred)/max(norm2(preferred),tiny(1.0_rkind))
    write(message,*) 'ADEnc hydroflow: iterations=',iterations,' water residual=',error,' flux correction=',change
    call write_log(trim(message))
    write(message,*) 'ADEnc hydroflow: zero-bound outlets=',bounded,' projection passes=',passes
    call write_log(trim(message))
    if (change>correction_limit) error stop 'Hydroflow correction exceeds configured limit; check ports/domain/hydrology'
    do p=1,size(port_kind)
      if (port_kind(p)==1 .and. flux(port_edge(p))<-tolerance) &
        error stop 'Hydroflow outlet reversed: an inflow concentration is required'
    end do
    cached_target=target; cached_flux=flux; projected_hour=raw_hour; projected_time=raw_time
    projection_cached=.true.
  end subroutine

  function hydro_width(e) result(w)
    use ncfluxarea, only: ncflux_active_width
    integer, intent(in) :: e
    real(kind=rkind) :: w
    w=ncflux_active_width(int(e,ikind))
  end function

  ! Hydraulic outlets may discharge or close, but cannot supply unconfigured
  ! incoming water/solute. Solve the weighted projection with F_out >= 0.
  ! NOT post-solve clipping: binding an outlet re-solves ALL water balances.
  ! This monotone active set is specific to exterior, one-owner edge bounds:
  ! B W B^T is an M-matrix. Removing a negative outlet decreases the remaining
  ! potentials, so bound multipliers remain nonnegative (no release needed).
  ! Interior one-way bounds are NOT supported by this algorithm.
  subroutine hydro_project_outflow(owner,preferred,weight,fixed,fixed_value,outflow,target,rtol,maxiter, &
                                  flux,ok,iterations,error,bounded,passes)
    integer, intent(in) :: owner(:,:),maxiter
    real(kind=rkind), intent(in) :: preferred(:),weight(:),fixed_value(:),target(:),rtol
    logical, intent(in) :: fixed(:),outflow(:)
    real(kind=rkind), intent(out) :: flux(:),error
    logical, intent(out) :: ok
    integer, intent(out) :: iterations,bounded,passes
    logical :: active(size(preferred)),negative(size(preferred)),solved
    real(kind=rkind) :: values(size(preferred))
    integer :: niter,i
    ok=.false.; iterations=0; bounded=0; passes=0; error=huge(1.0_rkind)
    if (size(flux)/=size(preferred)) return
    flux=preferred
    if (size(outflow)/=size(flux) .or. size(fixed)/=size(flux) .or. size(fixed_value)/=size(flux)) return
    if (size(owner,1)/=size(flux) .or. size(owner,2)/=2) return
    do i=1,size(flux)
      if (.not.outflow(i)) cycle
      if (fixed(i) .or. owner(i,2)/=0) return
    end do
    active=fixed; values=fixed_value
    do passes=1,count(outflow)+1
      call hydro_project(owner,preferred,weight,active,values,target,rtol,maxiter,flux,solved,niter,error)
      iterations=iterations+niter
      if (.not.solved) return
      negative=outflow .and. .not.active .and. flux<0
      if (.not.any(negative)) then
        ok=.true.
        return
      end if
      where(negative)
        active=.true.
        values=0
      end where
      bounded=bounded+count(negative)
    end do
  end subroutine

  ! min sum((F-F0)**2 / weight), subject to B F=target and prescribed edges.
  ! Matrix-free Jacobi-PCG on B W B^T; no changes to DRUtES transport solver.
  ! Closed compatible components use one potential gauge; inconsistent ones fail.
  subroutine hydro_project(owner,preferred,weight,fixed,fixed_value,target,rtol,maxiter, &
                           flux,ok,iterations,error)
    integer, intent(in) :: owner(:,:),maxiter
    real(kind=rkind), intent(in) :: preferred(:),weight(:),fixed_value(:),target(:),rtol
    logical, intent(in) :: fixed(:)
    real(kind=rkind), intent(out) :: flux(:),error
    logical, intent(out) :: ok
    integer, intent(out) :: iterations
    real(kind=rkind) :: diagonal(size(target)),rhs(size(target)),lambda(size(target)),r(size(target))
    real(kind=rkind) :: z(size(target)),p(size(target)),ap(size(target)),sumr(size(target))
    real(kind=rkind) :: rz,rznew,alpha,denom,scale,tol,delta
    integer :: root(size(target)),i,a,b,j,gauges(size(target))
    logical :: anchored(size(target)),pinned(size(target))
    ok=.false.; error=huge(1.0_rkind); iterations=0; flux=preferred
    if (size(target)==0 .or. size(flux)==0 .or. size(owner,1)/=size(flux) .or. size(owner,2)/=2) return
    if (size(weight)/=size(flux) .or. size(fixed)/=size(flux) .or. size(fixed_value)/=size(flux)) return
    if (.not.ieee_is_finite(rtol)) return
    if (rtol<=0 .or. maxiter<1) return
    if (any(.not.ieee_is_finite(preferred)) .or. any(.not.ieee_is_finite(weight)) .or. &
        any(.not.ieee_is_finite(fixed_value)) .or. any(.not.ieee_is_finite(target))) return
    if (any(weight<=0)) return
    diagonal=0; anchored=.false.; pinned=.false.
    do i=1,size(target)
      root(i)=i
    end do
    do i=1,size(flux)
      a=owner(i,1); b=owner(i,2)
      if (a<1 .or. a>size(target) .or. b<0 .or. b>size(target) .or. a==b) return
      if (fixed(i)) then
        flux(i)=fixed_value(i)
        cycle
      end if
      diagonal(a)=diagonal(a)+weight(i)
      if (b>0) then
        diagonal(b)=diagonal(b)+weight(i)
        j=find_root(root,a); gauges(1)=find_root(root,b); root(gauges(1))=j
      end if
    end do
    do i=1,size(flux)
      if (.not.fixed(i) .and. owner(i,2)==0) anchored(find_root(root,owner(i,1)))=.true.
    end do
    call incidence(flux,ap)
    rhs=target-ap
    scale=max(1.0_rkind,maxval(abs(target)),maxval(abs(flux)))
    tol=rtol*scale
    sumr=0; gauges=0
    do i=1,size(target)
      a=find_root(root,i); sumr(a)=sumr(a)+rhs(i)
      gauges(a)=i
    end do
    do i=1,size(target)
      if (gauges(i)==0 .or. anchored(i)) cycle
      if (abs(sumr(i))>tol) return
      pinned(gauges(i))=.true.
    end do
    where(diagonal==0) pinned=.true.
    where(pinned) diagonal=1
    lambda=0; r=rhs; where(pinned) r=0
    z=r/diagonal; p=z; rz=dot_product(r,z)
    do while(norm2(r)>tol)
      if (iterations>=maxiter .or. rz<=0 .or. .not.ieee_is_finite(rz)) return
      call action(p,ap)
      denom=dot_product(p,ap)
      if (denom<=0 .or. .not.ieee_is_finite(denom)) return
      alpha=rz/denom; lambda=lambda+alpha*p; r=r-alpha*ap
      iterations=iterations+1
      ! Recompute occasionally instead of trusting accumulated PCG roundoff.
      if (mod(iterations,50)==0) then
        call action(lambda,ap); r=rhs-ap; where(pinned) r=0
      end if
      z=r/diagonal; rznew=dot_product(r,z)
      p=z+(rznew/rz)*p; rz=rznew
    end do
    do i=1,size(flux)
      if (fixed(i)) cycle
      delta=lambda(owner(i,1))
      if (owner(i,2)>0) delta=delta-lambda(owner(i,2))
      flux(i)=flux(i)+weight(i)*delta
    end do
    call incidence(flux,ap)
    error=maxval(abs(ap-target))/scale
    ok=all(ieee_is_finite(flux)) .and. error<=10*rtol
  contains
    subroutine incidence(f,v)
      real(kind=rkind), intent(in) :: f(:)
      real(kind=rkind), intent(out) :: v(:)
      integer :: k
      v=0
      do k=1,size(f)
        v(owner(k,1))=v(owner(k,1))+f(k)
        if (owner(k,2)>0) v(owner(k,2))=v(owner(k,2))-f(k)
      end do
    end subroutine
    subroutine action(x,v)
      real(kind=rkind), intent(in) :: x(:)
      real(kind=rkind), intent(out) :: v(:)
      integer :: k,a,b
      real(kind=rkind) :: d
      v=0
      do k=1,size(flux)
        if (fixed(k)) cycle
        a=owner(k,1); b=owner(k,2); d=0
        if (.not.pinned(a)) d=x(a)
        if (b>0) then
          if (.not.pinned(b)) d=d-x(b)
        end if
        if (.not.pinned(a)) v(a)=v(a)+weight(k)*d
        if (b>0) then
          if (.not.pinned(b)) v(b)=v(b)-weight(k)*d
        end if
      end do
      where(pinned) v=x
    end subroutine
  end subroutine

  pure subroutine hydro_rt0(vertices,volume,edge_flux,xy,q,divq)
    real(kind=rkind), intent(in) :: vertices(3,2),volume,edge_flux(3),xy(2)
    real(kind=rkind), intent(out) :: q(2),divq
    integer :: k
    q=0
    do k=1,3
      q=q+edge_flux(k)*(xy-vertices(opposite(k),:))/(2*volume)
    end do
    divq=sum(edge_flux)/volume
  end subroutine

  subroutine hydro_value(el,xy,q,h,divq)
    use globals, only: elements,nodes
    integer, intent(in) :: el
    real(kind=rkind), intent(in) :: xy(2)
    real(kind=rkind), intent(out) :: q(2),h,divq
    real(kind=rkind) :: f(3)
    if (.not.ready) error stop 'Hydroflow coefficients before initialization'
    q=0; h=1; divq=0
    if (area(el)<=0) return
    if (trial) then
      f=orientation(:,el)*flux_trial(edge_of(:,el)); h=depth_trial(el)
    else
      f=orientation(:,el)*flux_state(edge_of(:,el)); h=depth_state(el)
    end if
    call hydro_rt0(nodes%data(elements%data(el,:),1:2),area(el),f,xy,q,divq)
  end subroutine

  function hydro_normal_flux(el,k) result(qn)
    integer, intent(in) :: el,k
    real(kind=rkind) :: qn
    integer :: edge
    if (.not.ready) error stop 'Hydroflow boundary before initialization'
    edge=edge_of(k,el)
    if (trial) then
      qn=orientation(k,el)*flux_trial(edge)/edge_length(edge)
    else
      qn=orientation(k,el)*flux_state(edge)/edge_length(edge)
    end if
  end function
end module nchydroflow
