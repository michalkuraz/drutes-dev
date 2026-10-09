! Isolated hydraulic + real FEM assembly tests. Never calls main/time stepping.
program test_adenc_hydroflow
  use typy
  use globals
  use global_objs
  use pde_objs
  use ncglobvars
  use nchydroflow
  use ncfluxarea
  use ncboundary
  use ncconservative
  use ncbalance
  use ncsupg
  use ncdispersion
  use lsconstitutive
  use ncpointers, only: adenc_element_corrections
  use dummy_procs
  use femmat
  use capmat
  use netcdf
  implicit none
  integer :: owner(5,2),iterations,i,j,k,e,mode,method,ierr,dims(3),variable_id,unit
  integer(kind=ikind) :: original(6)
  real(kind=rkind) :: preferred(5),weights(5),fixed_value(5),target(3),f(5),error
  real(kind=rkind) :: vertices(3,2),q(2),h,divq,v,normal(2),xy(2),divk(2),fd(2),plus(2,2),minus(2,2)
  real(kind=rkind) :: current(3),previous(3),mass(3,3),stiffness(3,3),rhs(3),gradient(3,2),basis(3)
  real(kind=rkind) :: data(1,1,4),system_matrix(3,3),old(3),new(3),qsave(2),hsave,dt,total_in,total_out,expected_out
  real(kind=rkind) :: water,load,water_sum,load_sum,water_saved,load_saved,initial_lateral
  real(kind=rkind) :: offset(2)
  logical :: fixed(5),ok
  logical :: lateral_case,clean_case,load_case
  character(len=1024) :: message
  character(len=32) :: argument
  call get_command_argument(1,argument)
  lateral_case=argument=='lateral-uniform' .or. argument=='lateral-clean' .or. argument=='lateral-zero' .or. &
    argument=='lateral-short' .or. argument=='lateral-load' .or. argument=='lateral-overlap'
  clean_case=argument=='lateral-clean'
  load_case=argument=='lateral-load'
  drutes_config%dimen=2; drutes_config%it_method=0; drutes_config%run_from_backup=.false.
  end_time=172800
  LSbank_noflow=.true.; LSconservative=.true.
  if (len_trim(argument)>0 .and. trim(argument)/='production' .and. trim(argument)/='linear' .and. &
      .not.lateral_case) then
    if (argument=='legacy') LSconservative=.false.
    call read_adenc_hydro()
    if (argument=='absent' .and. LShydro) error stop 'hydroflow default changed'
    if (argument=='off' .and. LShydro) error stop 'hydroflow off changed'
    print *, 'Hydroflow config checks passed'
    stop
  end if
  ! Confluence: two prescribed inlets, one outlet. A shared edge has ONE flux.
  owner(:,1)=[1,1,2,2,3]; owner(:,2)=[0,2,0,3,0]
  preferred=[-.3_rkind,.8_rkind,.3_rkind,-.1_rkind,-.2_rkind]
  weights=[1.0_rkind,2.0_rkind,3.0_rkind,1.0_rkind,1.0_rkind]
  fixed=[.true.,.false.,.false.,.false.,.true.]
  fixed_value=[-2.0_rkind,0.0_rkind,0.0_rkind,0.0_rkind,-1.0_rkind]; target=0
  call hydro_project(owner,preferred,weights,fixed,fixed_value,target,1.e-12_rkind,100,f,ok,iterations,error)
  if (.not.ok .or. maxval(abs(f-[-2.0_rkind,2.0_rkind,3.0_rkind,-1.0_rkind,-1.0_rkind]))>1.e-11_rkind) &
    error stop 'Confluence water balance'
  target=[-.1_rkind,-.2_rkind,-.3_rkind]
  call hydro_project(owner,preferred,weights,fixed,fixed_value,target,1.e-12_rkind,100,f,ok,iterations,error)
  if (.not.ok .or. abs(f(3)-2.4_rkind)>1.e-11_rkind) error stop 'Transient water storage balance'
  ! Closed component compatibility and gauge; no outlet invented by the solver.
  owner(:,1)=[1,1,2,2,3]; owner(:,2)=[0,2,0,3,0]
  fixed=[.true.,.false.,.true.,.false.,.true.]; fixed_value=0; target=0
  call hydro_project(owner,preferred,weights,fixed,fixed_value,target,1.e-12_rkind,100,f,ok,iterations,error)
  if (.not.ok .or. maxval(abs(f))>1.e-11_rkind) error stop 'Closed component gauge'
  target(1)=1
  call hydro_project(owner,preferred,weights,fixed,fixed_value,target,1.e-12_rkind,100,f,ok,iterations,error)
  if (ok) error stop 'Incompatible closed water balance accepted'
  call check_outflow_bounds()
  ! RT0 edge traces and divergence, including a nonzero divergence.
  vertices(:,1)=[0,2,0]; vertices(:,2)=[0,0,1]
  f(1:3)=[-.5_rkind,1.2_rkind,.3_rkind]
  do k=1,3
    i=k; j=mod(k,3)+1
    xy=(vertices(i,:)+vertices(j,:))/2
    normal=[vertices(j,2)-vertices(i,2),vertices(i,1)-vertices(j,1)]
    call hydro_rt0(vertices,1.0_rkind,f(:3),xy,q,divq)
    if (abs(dot_product(q,normal)-f(k))>1.e-12_rkind) error stop 'RT0 normal trace'
    if (abs(divq-sum(f(:3)))>1.e-12_rkind) error stop 'RT0 divergence'
  end do
  ! Exact div(K) checked by numerical differentiation of the production tensor.
  xy=[.6_rkind,.2_rkind]; call hydro_rt0(vertices,1.0_rkind,f(:3),xy,q,divq)
  divk=rt0_dispersion_divergence(q,divq,2.0_rkind,.2_rkind); fd=0
  do i=1,2
    normal=0; normal(i)=1.e-6_rkind
    call ls_dispersion_tensor(q+.5_rkind*divq*normal,2.0_rkind,.2_rkind,plus)
    call ls_dispersion_tensor(q-.5_rkind*divq*normal,2.0_rkind,.2_rkind,minus)
    fd=fd+(plus(:,i)-minus(:,i))/(2.e-6_rkind)
  end do
  if (maxval(abs(divk-fd))>1.e-8_rkind) error stop 'RT0 dispersion divergence'
  ! Constant C residual: Hnew-Hold + dt*div(q)=0, including SUPG and shock.
  gradient(:,1)=[-1,1,0]; gradient(:,2)=[-1,0,1]; basis=1.0_rkind/3
  call supg_terms([1.0_rkind,.2_rkind],1.2_rkind,plus,gradient,basis,1.0_rkind,.true., &
    0.0_rkind,0.0_rkind,2.0_rkind,mass,stiffness,rhs,-.2_rkind,[0.0_rkind,0.0_rkind])
  current=1; previous=1
  if (maxval(abs(matmul(mass+stiffness,current)-matmul(mass,previous)/1.2_rkind))>1.e-12_rkind) &
    error stop 'SUPG constant-state residual'
  if (trim(argument)/='production' .and. trim(argument)/='linear' .and. .not.lateral_case) then
    print *, 'Hydroflow kernel checks passed'
    stop
  end if

  ! Four-triangle strip, two inlet edges sharing one ID, two explicit outlets.
  allocate(pde(1),pde_common%xvect(3,4),pde_common%bvect(3))
  allocate(pde(1)%pde_fnc(1),pde(1)%permut(6),pde(1)%bc(101:103))
  nodes%kolik=6; elements%kolik=4
  offset=[500000.0_rkind,5000000.0_rkind]
  ! Overlapping groups deliberately produce strongly nonuniform divergence.
  ! Use a well-conditioned local origin for this sub-metre synthetic fixture;
  ! all original production cases retain their realistic UTM coordinates.
  if (argument=='lateral-overlap') offset=1000
  allocate(nodes%data(6,2),nodes%edge(6),nodes%element(6),elements%data(4,3))
  allocate(elements%ders(4,3,2),elements%areas(4),elements%material(4))
  nodes%data=0; nodes%data(:,1)=[0,1,0,1,0,1]+offset(1)
  nodes%data(:,2)=[0.0_rkind,0.0_rkind,.5_rkind,.5_rkind,1.0_rkind,1.0_rkind]+offset(2)
  elements%data(1,:)=[1,2,4]; elements%data(2,:)=[1,4,3]
  elements%data(3,:)=[3,4,6]; elements%data(4,:)=[3,6,5]
  do e=1,4
    do i=1,3
      call nodes%element(elements%data(e,i))%fill(int(e,ikind))
    end do
  end do
  do e=1,4,2
    elements%ders(e,:,1)=[-1,1,0]; elements%ders(e,:,2)=[0,-2,2]
    elements%ders(e+1,:,1)=[0,1,-1]; elements%ders(e+1,:,2)=[-2,0,2]
  end do
  elements%areas=.25_rkind; elements%material=1
  allocate(ncfluxdata%activeel(4),ncfluxdata%fluxvct(4,2))
  ncfluxdata%activeel=.true.; ncfluxdata%fluxvct(:,1)=1; ncfluxdata%fluxvct(:,2)=0
  original=[101,0,101,0,101,0]; nodes%edge=original; addedbc=103
  ! Synthetic NetCDF with an actual second/third daily slice, read by real code.
  call nc_check(nf90_create('forcing.nc',nf90_clobber,netcdfID))
  call nc_check(nf90_def_dim(netcdfID,'lon',1,dims(1)))
  call nc_check(nf90_def_dim(netcdfID,'lat',1,dims(2)))
  call nc_check(nf90_def_dim(netcdfID,'time',4,dims(3)))
  call nc_check(nf90_def_var(netcdfID,'Qrouted',nf90_double,dims,variable_id))
  call nc_check(nf90_enddef(netcdfID))
  data(1,1,:)=[2.0_rkind,2.4_rkind,2.8_rkind,3.0_rkind]
  call nc_check(nf90_put_var(netcdfID,variable_id,data))
  ncfluxdata%initialized=.true.; ncfluxdata%slice_loaded=.true.; ncfluxdata%current_time_index=1
  ncfluxdata%q_varid=variable_id; ncfluxdata%nlat=1; ncfluxdata%nlon=1; ncfluxdata%ntime=4
  ncfluxdata%has_bounds=.true.
  allocate(ncfluxdata%time(4),ncfluxdata%lon(1),ncfluxdata%lat(1),ncfluxdata%qslice(1,1))
  allocate(ncfluxdata%lon_bnds(1,2),ncfluxdata%lat_bnds(1,2))
  ncfluxdata%time=[0,24,48,72]; ncfluxdata%lon=9; ncfluxdata%lat=45; ncfluxdata%qslice=2
  ncfluxdata%lon_bnds(1,:)=[-180.0_rkind,180.0_rkind]; ncfluxdata%lat_bnds(1,:)=[0.0_rkind,90.0_rkind]
  allocate(ncnodes%data(4,2),ncelements%data(1,4),el2ncgrid(4))
  ncnodes%data(:,1)=[-1,3,3,-1]+offset(1)
  ncnodes%data(:,2)=[-1,-1,3,3]+offset(2)
  ncelements%kolik=1; ncelements%data(1,:)=[1,2,3,4]; el2ncgrid=1; Qmin=0; ora_di_ini=0
  call ncflux_prepare_widths(ok,message)
  if (.not.ok) error stop 'Synthetic width preparation'
  call read_adenc_hydro()
  if (.not.LShydro) error stop 'Production hydro test config absent'
  call hydro_filter(); call prepare_adenc_banks(original); call hydro_initialize(original)
  pde(1)%pde_fnc(1)%convection=>ADEls_convection; pde(1)%flux=>ncflux
  pde(1)%pde_fnc(1)%dispersion=>ADElsdisp; LSdisp=.05_rkind; LSdisp_transverse=.005_rkind
  pde(1)%pde_fnc(1)%elasticity=>ADEls_tder_coef; pde(1)%pde_fnc(1)%reaction=>dummy_scalar
  pde(1)%pde_fnc(1)%zerord=>dummy_scalar; pde(1)%pde_fnc(1)%der_convect=>dummy_vector
  pde(1)%getval=>getvalp1; pde(1)%stabilize_element=>adenc_element_corrections
  pde(1)%step_begin=>conservative_begin; pde(1)%step_end=>conservative_end
  pde(1)%boundary_history=>adenc_boundary_history; pde(1)%permut=[0,1,0,2,0,3]
  do i=101,103
    pde(1)%bc(i)%code=1; pde(1)%bc(i)%file=.false.; pde(1)%bc(i)%value=1
    pde(1)%bc(i)%value_fnc=>ADEls_dirichlet
  end do
  allocate(stiff_mat(3,3),cap_mat(3,3),bside(3),elnode_prev(3))
  allocate(base_fnc(3,3),gauss_points%weight(3),gauss_points%point(3,2))
  base_fnc(:,1)=[2.0_rkind/3,1.0_rkind/6,1.0_rkind/6]
  base_fnc(:,2)=[1.0_rkind/6,2.0_rkind/3,1.0_rkind/6]
  base_fnc(:,3)=[1.0_rkind/6,1.0_rkind/6,2.0_rkind/3]
  gauss_points%point(:,1)=base_fnc(2,:); gauss_points%point(:,2)=base_fnc(3,:)
  gauss_points%weight=1.0_rkind/6; gauss_points%area=.5_rkind
  call spmatrix%init(3_ikind,3_ikind)
  total_in=hydro_normal_flux(2,3)*.5_rkind+hydro_normal_flux(4,3)*.5_rkind
  total_out=hydro_normal_flux(1,2)*.5_rkind+hydro_normal_flux(3,2)*.5_rkind
  expected_out=2
  if (argument=='linear' .or. lateral_case) expected_out=2-(2.4_rkind/(4*hydro_velocity(2.4_rkind,4.0_rkind))- &
    2.0_rkind/(4*hydro_velocity(2.0_rkind,4.0_rkind)))/86400
  initial_lateral=0
  do e=1,4
    call hydro_sources(e,water,load)
    initial_lateral=initial_lateral+elements%areas(e)*water
  end do
  if (lateral_case) then
    if (argument=='lateral-zero') then
      if (initial_lateral/=0) error stop 'Zero lateral source created water'
    else
      if (abs(initial_lateral-.2_rkind)>1.e-12_rkind) error stop 'Lateral Q repeated per element or group'
    end if
  end if
  expected_out=expected_out+initial_lateral
  if (abs(total_in+2)>1.e-10_rkind .or. abs(total_out-expected_out)>1.e-10_rkind) &
    error stop 'Qrouted inlet repeated per edge instead of shared per port'
  if (abs(hydro_normal_flux(1,1))+abs(hydro_normal_flux(4,2))>1.e-12_rkind) &
    error stop 'Water crossed a sealed bank'
  if (abs(hydro_normal_flux(1,3)+hydro_normal_flux(2,1))>1.e-12_rkind) &
    error stop 'Two traces on a shared edge disagree'
  do method=1,2
    pde_common%timeint_method=method
    if (method==1) pde_common%time_integ=>impl_euler_np_diag
    if (method==2) pde_common%time_integ=>impl_euler_np_nondiag
    do mode=0,2
      call hydro_initialize(original)
      if (allocated(LSdepth_old)) deallocate(LSdepth_old,LSdepth_new)
      LSstate_time=0; LSprevious_time=0; LSstep_active=.false.; time=0
      call balance_reset(); LSbalance=.true.
      LSsupg=mode>0; LSsupg_factor=2; LSshock=mode==2; LSshock_factor=1
      old=1; pde_common%xvect=1
      do k=1,3
        dt=600
        if (k==1) dt=86400
        time_step=dt; pde_common%xvect(:,1)=old; pde_common%xvect(:,2)=old
        call pde(1)%step_begin(time,time_step)
        if (lateral_case .and. k==1) then
          if (time_step/=3600) error stop 'Lateral knot not clipped'
          water_sum=0; load_sum=0
          do e=1,4
            call hydro_sources(e,water,load)
            water_sum=water_sum+elements%areas(e)*water
            load_sum=load_sum+elements%areas(e)*load
          end do
          if (argument/='lateral-zero') then
            if (abs(water_sum-.3_rkind)>1.e-12_rkind) error stop 'Lateral interval water integral'
            if (clean_case) then
              if (load_sum/=0) error stop 'Clean lateral water added solute'
            else if (load_case) then
              if (abs(load_sum-.5_rkind)>1.e-12_rkind) error stop 'Q*C load interpolation or integration'
            else
              if (abs(load_sum-water_sum)>1.e-12_rkind) error stop 'Matching lateral C1 lost solute'
            end if
          else
            if (water_sum/=0 .or. load_sum/=0) error stop 'Zero water added solute'
          end if
        end if
        call assemble_mat(ierr)
        ! assemble_mat's historical ierr argument is not assigned by that routine.
        ! Check the actual assembled equation and solved residual below instead.
        do j=1,3
          do i=1,3
            system_matrix(i,j)=spmatrix%get(int(i,ikind),int(j,ikind))
          end do
        end do
        if (.not.(clean_case .or. load_case) .and. maxval(abs(matmul(system_matrix,old)-pde_common%bvect))/ &
            max(1.0_rkind,maxval(abs(system_matrix)))>1.e-10_rkind) then
          print *, 'Constant residual, capacity method/mode/step:',method,mode,k
          print *, matmul(system_matrix,old)-pde_common%bvect
          error stop 'Production constant-state matrix residual'
        end if
        call dense_solve(system_matrix,pde_common%bvect,new)
        if (.not.(clean_case .or. load_case) .and. maxval(abs(new-1))>1.e-9_rkind) &
          error stop 'Uniform concentration not preserved'
        if (clean_case .and. k==1 .and. minval(new)>=.999_rkind) error stop 'Clean lateral water did not dilute'
        if (load_case .and. k==1 .and. maxval(new)<=1.001_rkind) error stop 'Lateral contaminant load missing'
        pde_common%xvect(:,3)=new; call pde(1)%step_end(.true.)
        if (abs(balance_error)>1.e-9_rkind) error stop 'Hydro transport inventory budget'
        time=time+time_step; old=new
      end do
      xy=offset+[.5_rkind,.25_rkind]
      call hydro_value(1,xy,qsave,hsave,divq)
      call hydro_sources(1,water_saved,load_saved)
      dt=600; call pde(1)%step_begin(172800.0_rkind,dt)
      call hydro_value(1,xy,q,h,divq)
      if (abs(h-hsave)<1.e-6_rkind) error stop 'New forcing storage not loaded'
      call pde(1)%step_end(.false.)
      call hydro_value(1,xy,q,h,divq)
      if (maxval(abs(q-qsave))>1.e-12_rkind .or. h/=hsave) error stop 'Rejected hydro state committed'
      call hydro_sources(1,water,load)
      if (water/=water_saved .or. load/=load_saved) error stop 'Rejected lateral source committed'
    end do
  end do
  if (argument=='lateral-short') then
    dt=600
    call pde(1)%step_begin(172801.0_rkind,dt)
    error stop 'Expired lateral forcing was accepted'
  end if
  call nc_check(nf90_close(netcdfID))
  print *, 'Hydroflow production checks passed: shared Qrouted inlet, storage, uniform C, modes and rejection'
contains
  subroutine check_outflow_bounds()
    integer :: own(8,2),nit,bound,passes,trial_case,subset,l,seed_size
    integer, allocatable :: seed(:)
    logical :: prescribed(8),outlets(8),reference_fixed(8),success
    real(kind=rkind) :: pref(8),w(8),values(8),t(4),answer(8),candidate(8),best(8)
    real(kind=rkind) :: err,objective,best_objective,randoms(20)
    ! A cascade: the second outlet reverses only AFTER the first is closed.
    own(:3,1)=1; own(:3,2)=0; pref(:3)=[-5.0_rkind,1.0_rkind,10.0_rkind]
    prescribed=.false.; outlets=.true.; values=0; w=1; t(1)=8
    call hydro_project_outflow(own(:3,:),pref(:3),w(:3),prescribed(:3),values(:3),outlets(:3),t(:1), &
      1.e-12_rkind,100,answer(:3),success,nit,err,bound,passes)
    if (.not.success .or. bound/=2 .or. passes/=3 .or. &
        maxval(abs(answer(:3)-[0.0_rkind,0.0_rkind,8.0_rkind]))>1.e-11_rkind) &
      error stop 'Outlet cascade was clipped instead of balanced'
    t(1)=0
    call hydro_project_outflow(own(:3,:),pref(:3),w(:3),prescribed(:3),values(:3),outlets(:3),t(:1), &
      1.e-12_rkind,100,answer(:3),success,nit,err,bound,passes)
    if (.not.success .or. maxval(abs(answer(:3)))>1.e-11_rkind) error stop 'Zero net outflow failed'
    t(1)=-1
    call hydro_project_outflow(own(:3,:),pref(:3),w(:3),prescribed(:3),values(:3),outlets(:3),t(:1), &
      1.e-12_rkind,100,answer(:3),success,nit,err,bound,passes)
    if (success) error stop 'Required incoming water silently accepted at outlet'
    ! Independent components cannot exchange water through the optimizer.
    own(:2,1)=[1,2]; own(:2,2)=0; t(:2)=[1.0_rkind,-1.0_rkind]
    call hydro_project_outflow(own(:2,:),pref(:2),w(:2),prescribed(:2),values(:2),outlets(:2),t(:2), &
      1.e-12_rkind,100,answer(:2),success,nit,err,bound,passes)
    if (success) error stop 'Incompatible separate outlet component accepted'
    own(1,:)=[1,2]
    call hydro_project_outflow(own(:2,:),pref(:2),w(:2),prescribed(:2),values(:2),outlets(:2),t(:2), &
      1.e-12_rkind,100,answer(:2),success,nit,err,bound,passes)
    if (success) error stop 'Interior one-way bound accepted'
    ! Compare 100 weighted graph problems to exhaustive active-subset search.
    ! This reference enumerates ALL possible outlet bounds, not the monotone rule.
    own(:,1)=[1,1,2,3,1,2,3,4]; own(:,2)=[0,2,3,4,4,0,0,0]
    prescribed=.false.; prescribed(1)=.true.; outlets=.false.; outlets(6:8)=.true.
    values=0; values(1)=-2
    call random_seed(size=seed_size); allocate(seed(seed_size)); seed=314159
    call random_seed(put=seed)
    do trial_case=1,100
      call random_number(randoms)
      pref=10*(randoms(:8)-.5_rkind); w=.1_rkind+3*randoms(9:16)
      t=.6_rkind*(randoms(17:20)-.5_rkind)
      call hydro_project_outflow(own,pref,w,prescribed,values,outlets,t,1.e-12_rkind,100, &
        answer,success,nit,err,bound,passes)
      if (.not.success .or. any(answer(6:8)<0)) error stop 'Feasible graph outflow projection failed'
      best_objective=huge(1.0_rkind); best=0
      do subset=0,7
        reference_fixed=prescribed
        do l=0,2
          reference_fixed(6+l)=btest(subset,l)
        end do
        call hydro_project(own,pref,w,reference_fixed,values,t,1.e-12_rkind,100,candidate,success,nit,err)
        if (.not.success) cycle
        if (any(candidate(6:8)<-1.e-10_rkind)) cycle
        objective=sum((candidate-pref)**2/w)
        if (objective<best_objective) then
          best_objective=objective; best=candidate
        end if
      end do
      if (maxval(abs(answer-best))>1.e-8_rkind) error stop 'Bounded projection not weighted minimizer'
    end do
    print *, 'Outflow bounds: cascade, zero flow, infeasibility, components and 100 reference graphs passed'
  end subroutine

  subroutine nc_check(status)
    integer, intent(in) :: status
    if (status/=nf90_noerr) error stop 'Synthetic NetCDF failure'
  end subroutine
  subroutine dense_solve(a,b,x)
    real(kind=rkind), intent(in) :: a(3,3),b(3)
    real(kind=rkind), intent(out) :: x(3)
    real(kind=rkind) :: aa(3,3),bb(3),factor
    integer :: i,j
    aa=a; bb=b
    do j=1,2
      do i=j+1,3
        factor=aa(i,j)/aa(j,j); aa(i,j:)=aa(i,j:)-factor*aa(j,j:); bb(i)=bb(i)-factor*bb(j)
      end do
    end do
    do i=3,1,-1
      x(i)=(bb(i)-dot_product(aa(i,i+1:),x(i+1:)))/aa(i,i)
    end do
  end subroutine
end program
