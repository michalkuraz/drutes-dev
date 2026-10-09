! Full production initialization/dry assembly only. Run in a NEW directory.
! Never calls main or solve_pde, never advances/accepts a transport time step.
program check_adenc_routing_init
  use typy
  use globals
  use pde_objs
  use ncglobvars
  use ncrouting
  use ncfluxarea, only: ncflux_active_width
  use drutes_init, only: parse_globals,init_measured,init_observe
  use manage_pointers, only: set_pointers
  use feminittools, only: feminit
  use femmat, only: assemble_mat
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: e,ierr,ncheck
  integer(kind=ikind) :: target,code
  integer(kind=ikind), allocatable :: upstream(:)
  real(kind=rkind) :: direction(2),width
  logical :: ok
  call parse_globals()
  call init_measured()
  call set_pointers()
  if (.not.routing_enabled) error stop 'mandatory routing not initialized'
  ncheck=0
  do e=1,elements%kolik
    if (.not.ncfluxdata%activeel(e)) cycle
    call routing_element(int(e,ikind),target,upstream,direction,code,ok)
    if (.not.ok) error stop 'active FE missing routing node'
    if (target>0) then
      if (norm2(direction-ncfluxdata%fluxvct(e,:))>1.e-12_rkind) error stop 'routing direction not used'
      ncheck=ncheck+1
    end if
    width=ncflux_active_width(int(e,ikind))
    if (.not.ieee_is_finite(width)) error stop 'nonfinite width'
    if (LSconservative .and. width<=0) error stop 'invalid active width after hydro filtering'
  end do
  if (ncheck==0) error stop 'no network directions checked'
  call init_observe()
  call feminit()
  if (any(.not.ncfluxdata%activeel(observation_array(:)%element))) error stop 'inactive observation point'
  time=0; time_step=init_dt
  if (associated(pde(1)%step_begin)) call pde(1)%step_begin(time,time_step)
  call assemble_mat(ierr)
  ! The legacy assemble_mat does not assign ierr: inspect actual matrices.
  if (.not.all(ieee_is_finite(pde_common%bvect))) error stop 'nonfinite assembled rhs'
  if (.not.all(ieee_is_finite(cap_mat)) .or. .not.all(ieee_is_finite(stiff_mat))) &
    error stop 'nonfinite local assembly'
  if (associated(pde(1)%step_end)) call pde(1)%step_end(.false.)
  print *, 'Routing production initialization/assembly passed; mapped FE=',ncheck
  print *, 'Observation FE indices=',observation_array(:)%element
end program
