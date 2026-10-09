! Optional conservative ADEnc local assembly using the existing post-capacity hook.
module ncconservative
  use typy
  use ncglobvars
  implicit none
  private
  public :: read_adenc_conservative,conservative_begin,conservative_end,conservative_element
  public :: adenc_boundary_history,adenc_boundary_value
contains
  subroutine read_adenc_conservative()
    use globals, only: drutes_config
    use readtools, only: fileread
    use core_tools, only: write_log
    use ncbalance, only: balance_reset
    integer :: unit
    logical :: exists
    LSconservative=.false.; LSbalance=.false.; LSstep_active=.false.; LSclock_override=.false.
    LSstate_time=0; LSprevious_time=0
    call balance_reset()
    if (allocated(LSdepth_old)) deallocate(LSdepth_old,LSdepth_new)
    inquire(file='drutes.conf/netcdf/conservative.conf',exist=exists)
    if (.not. exists) return
    open(newunit=unit,file='drutes.conf/netcdf/conservative.conf',status='old',action='read')
    call fileread(LSconservative,unit)
    call fileread(LSbalance,unit)
    close(unit)
    if (.not. (LSconservative .or. LSbalance)) return
    if (drutes_config%run_from_backup) &
      error stop 'ADEnc conservative/audit restart from backup is not supported'
    if (drutes_config%dimen/=2 .or. drutes_config%it_method/=0 .or. .not. LSbank_noflow) &
      error stop 'ADEnc conservative/audit requires 2D standard Picard and no-flow active banks'
    call write_log('ADEnc optional conservative transport / inventory audit configured')
  end subroutine

  function adenc_boundary_value(pde_loc,el,node,t) result(value)
    use pde_objs, only: pde_str
    use globals, only: nodes,elements
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el,node
    real(kind=rkind), intent(in) :: t
    real(kind=rkind) :: value
    integer(kind=ikind) :: edge,j,selected
    edge=nodes%edge(elements%data(el,node))
    value=pde_loc%bc(edge)%value
    if (.not. pde_loc%bc(edge)%file) return
    selected=1
    do j=1,size(pde_loc%bc(edge)%series,1)
      if (pde_loc%bc(edge)%series(j,1)>t) exit
      selected=j
    end do
    value=pde_loc%bc(edge)%series(selected,2) ! last value persists after final record
  end function

  function adenc_boundary_history(pde_loc,el,node,column) result(value)
    use pde_objs, only: pde_str
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el,node,column
    real(kind=rkind) :: value,t
    t=LSstate_time
    if (column==4) t=LSprevious_time
    if (LSstep_active .and. column/=1 .and. column/=4) t=LStrial_time
    value=adenc_boundary_value(pde_loc,el,node,t)
  end function

  subroutine conservative_begin(pde_loc,t,dt)
    use pde_objs
    use globals, only: elements,gauss_points
    use global_objs, only: integpnt_str
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use nchydroflow, only: hydro_begin
    class(pde_str), intent(in) :: pde_loc
    real(kind=rkind), intent(in) :: t
    real(kind=rkind), intent(in out) :: dt
    type(integpnt_str) :: point
    integer(kind=ikind) :: el,g,edge,j
    real(kind=rkind) :: event,next_day
    if (pde_common%timeint_method/=1 .and. pde_common%timeint_method/=2) &
      error stop 'ADEnc conservative/audit requires transient implicit Euler'
    if (size(pde)/=1 .or. size(elements%data,2)/=3 .or. .not. pde_loc%diffusion) &
      error stop 'ADEnc conservative/audit requires a single 2D P1 equation'
    if (LSconservative) then
      do edge=lbound(pde_loc%bc,1),ubound(pde_loc%bc,1)
        if (pde_loc%bc(edge)%code/=1) error stop 'ADEnc conservative requires Dirichlet port records'
        if (.not. pde_loc%bc(edge)%file) cycle
        if (.not. allocated(pde_loc%bc(edge)%series)) error stop 'Missing ADEnc inlet series'
        if (size(pde_loc%bc(edge)%series,2)/=2 .or. size(pde_loc%bc(edge)%series,1)<1) &
          error stop 'ADEnc inlet series must have two columns'
        if (any(.not.ieee_is_finite(pde_loc%bc(edge)%series))) error stop 'Nonfinite ADEnc inlet series'
        do j=2,size(pde_loc%bc(edge)%series,1)
          if (pde_loc%bc(edge)%series(j,1)<=pde_loc%bc(edge)%series(j-1,1)) &
            error stop 'ADEnc inlet times must be strictly increasing'
        end do
      end do
    end if
    if (LSconservative) then
      next_day=(floor(t/86400.0_rkind)+1.0_rkind)*86400.0_rkind
      dt=min(dt,next_day-t)
      do edge=lbound(pde_loc%bc,1),ubound(pde_loc%bc,1)
        if (.not. pde_loc%bc(edge)%file) cycle
        do j=1,size(pde_loc%bc(edge)%series,1)
          event=pde_loc%bc(edge)%series(j,1)
          if (event>t) dt=min(dt,event-t)
        end do
      end do
    end if
    if (.not. ieee_is_finite(dt) .or. dt<=0 .or. t+dt<=t) error stop 'Invalid ADEnc trial time step'
    LSstep_start=t; LSstep_dt=dt
    ! Left limit at discontinuities: a step ending at pulse/day change uses
    ! the interval's forcing. Accepted H and BC history retain that SAME trace.
    LStrial_time=nearest(t+dt,-1.0_rkind)
    if (.not. LSconservative) LStrial_time=t ! diagnostic-only preserves legacy timing
    if (.not. allocated(LSdepth_old)) then
      allocate(LSdepth_old(size(gauss_points%weight),elements%kolik), &
        LSdepth_new(size(gauss_points%weight),elements%kolik))
      LSdepth_old=1; LSdepth_new=1
      point%type_pnt='gqnd'; point%column=1
      LSclock_override=.true.; LSoverride_time=LSstate_time
      do el=1,elements%kolik
        if (.not. ncfluxdata%activeel(el)) cycle
        point%element=el
        do g=1,size(gauss_points%weight)
          point%order=g
          LSdepth_old(g,el)=pde_loc%pde_fnc(1)%elasticity(pde_loc,elements%material(el),point)
          if (.not. ieee_is_finite(LSdepth_old(g,el)) .or. LSdepth_old(g,el)<=0) &
            error stop 'Invalid ADEnc initial storage depth'
        end do
      end do
      LSclock_override=.false.
    end if
    call hydro_begin(LStrial_time,dt)
    LSstep_active=.true.
  end subroutine

  subroutine conservative_end(pde_loc,accepted)
    use pde_objs, only: pde_str
    use ncbalance, only: balance_finish
    use nchydroflow, only: hydro_end
    class(pde_str), intent(in) :: pde_loc
    logical, intent(in) :: accepted
    if (LSbalance) call balance_finish(pde_loc,accepted,LSstep_start,LSstep_dt)
    if (accepted) then
      LSdepth_old=LSdepth_new
      LSprevious_time=LSstate_time; LSstate_time=LStrial_time
    end if
    LSstep_active=.false.
    call hydro_end(accepted)
  end subroutine

  subroutine conservative_element(pde_loc,el,dt,oldcap,newcap,open_matrix,source)
    use pde_objs
    use globals
    use global_objs
    use ncboundary, only: adenc_open_element
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el
    real(kind=rkind), intent(in) :: dt
    real(kind=rkind), intent(out) :: oldcap(3,3),newcap(3,3),open_matrix(3,3),source
    type(integpnt_str) :: point
    real(kind=rkind) :: q(2),h,w,tmp,row_sum(3)
    integer :: i,j,g
    if (.not. LSstep_active) error stop 'ADEnc assembly without step preparation'
    oldcap=0; newcap=cap_mat; open_matrix=0
    ! time_integ already added the old-concentration capacity load to bside.
    ! Only the remaining load represents a physical source, not inventory.
    source=-sum(bside-matmul(newcap,elnode_prev))
    point%type_pnt='gqnd'; point%element=el; point%column=2
    do g=1,size(gauss_points%weight)
      point%order=g
      h=pde_loc%pde_fnc(1)%elasticity(pde_loc,elements%material(el),point)
      if (.not. ieee_is_finite(h) .or. h<=0) error stop 'Invalid ADEnc storage depth'
      LSdepth_new(g,el)=h
      w=elements%areas(el)*gauss_points%weight(g)/gauss_points%area
      call pde_loc%pde_fnc(1)%convection(pde_loc,elements%material(el),point,vector_out=q)
      do j=1,3
        do i=1,3
          oldcap(i,j)=oldcap(i,j)-w*LSdepth_old(g,el)*base_fnc(i,g)*base_fnc(j,g)
          if (LSconservative) then
            ! Remove -Ni*q.grad(Nj), add +grad(Ni).q*Nj (negative DRUtES sign).
            tmp=base_fnc(i,g)*dot_product(elements%ders(el,j,1:2),q) + &
              dot_product(elements%ders(el,i,1:2),q)*base_fnc(j,g)
            stiff_mat(i,j)=stiff_mat(i,j)+dt*w*tmp
          end if
        end do
      end do
    end do
    if (pde_common%timeint_method==1) then
      row_sum=sum(oldcap,dim=2); oldcap=0
      do i=1,3
        oldcap(i,i)=row_sum(i)
      end do
    end if
    if (LSconservative) then
      bside=bside+matmul(oldcap-newcap,elnode_prev)
      call adenc_open_element(el,dt,open_matrix)
      stiff_mat=stiff_mat+open_matrix
    end if
  end subroutine
end module
