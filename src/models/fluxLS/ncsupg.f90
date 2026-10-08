! Residual-based SUPG for the depth-integrated, P1, 2D ADEnc equation.
module ncsupg
  use typy
  implicit none
  private
  public :: read_adenc_supg, adenc_supg_element, supg_tau, supg_terms, shock_diffusivity, shock_terms
contains
  subroutine read_adenc_supg()
    use ncglobvars, only: LSsupg, LSsupg_factor, LSshock, LSshock_factor
    use readtools, only: fileread
    use core_tools, only: write_log
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer :: unit, ierr
    logical :: exists
    character(len=256) :: message
    LSsupg=.false.
    LSsupg_factor=1.0_rkind
    inquire(file='drutes.conf/netcdf/supg.conf',exist=exists)
    if (exists) then
      open(newunit=unit,file='drutes.conf/netcdf/supg.conf',status='old',action='read',iostat=ierr)
      if (ierr/=0) error stop 'Unable to open ADEnc supg.conf'
      call fileread(LSsupg,unit)
      call fileread(LSsupg_factor,unit)
      if (.not. ieee_is_finite(LSsupg_factor)) error stop 'SUPG factor must be finite'
      if (LSsupg_factor<0) error stop 'SUPG factor must be nonnegative'
      close(unit)
    end if
    write(message,*) 'ADEnc SUPG enabled=',LSsupg,' factor=',LSsupg_factor
    call write_log(trim(message))
    LSshock=.false.
    LSshock_factor=1.0_rkind
    inquire(file='drutes.conf/netcdf/shock.conf',exist=exists)
    if (exists) then
      open(newunit=unit,file='drutes.conf/netcdf/shock.conf',status='old',action='read',iostat=ierr)
      if (ierr/=0) error stop 'Unable to open ADEnc shock.conf'
      call fileread(LSshock,unit)
      call fileread(LSshock_factor,unit)
      if (.not. ieee_is_finite(LSshock_factor)) error stop 'Shock factor must be finite'
      if (LSshock_factor<0) error stop 'Shock factor must be nonnegative'
      close(unit)
    end if
    write(message,*) 'ADEnc shock capturing enabled=',LSshock,' factor=',LSshock_factor
    call write_log(trim(message))
  end subroutine read_adenc_supg

  ! Time/streamline/diffusion/reaction scales, in seconds. Gradients are P1.
  pure function supg_tau(velocity, diffusion, gradients, dt, transient, rate) result(tau)
    real(kind=rkind), intent(in) :: velocity(2), diffusion(2,2), gradients(3,2), dt, rate
    logical, intent(in) :: transient
    real(kind=rkind) :: tau, speed, h, d(2), parallel_diffusion, scales(4), maximum
    tau=0
    speed=norm2(velocity)
    if (speed<=0) return
    d=velocity/speed
    maximum=sum(abs(matmul(gradients,d)))
    if (maximum<=0) return
    h=2/maximum
    parallel_diffusion=max(0.0_rkind,dot_product(d,matmul(diffusion,d)))
    scales=[0.0_rkind,2*speed/h,4*parallel_diffusion/h**2,abs(rate)]
    if (transient) then
      if (dt<=0) return
      scales(1)=2/dt
    end if
    maximum=maxval(scales)
    if (maximum>0) tau=(1/maximum)/norm2(scales/maximum)
  end function supg_tau

  ! Return NEGATIVE-sign DRUtES corrections at one quadrature point.
  ! Strong residual: h C_t + q.grad(C) - reaction*C - source.
  ! P1 Hessian is zero; q and K are piecewise constant inside each hydro cell.
  pure subroutine supg_terms(q, depth, tensor, gradients, basis, dt, transient, &
                            reaction, source, factor, mass, stiffness, rhs)
    real(kind=rkind), intent(in) :: q(2), depth, tensor(2,2), gradients(3,2), basis(3)
    real(kind=rkind), intent(in) :: dt, reaction, source, factor
    logical, intent(in) :: transient
    real(kind=rkind), intent(out) :: mass(3,3), stiffness(3,3), rhs(3)
    real(kind=rkind) :: tau, test(3), convective_derivative(3)
    integer :: i,j
    mass=0; stiffness=0; rhs=0
    if (depth<=0 .or. factor<=0) return
    tau=factor*supg_tau(q/depth,tensor/depth,gradients,dt,transient,reaction/depth)
    test=tau*matmul(gradients,q/depth)
    convective_derivative=matmul(gradients,q)
    do i=1,3
      do j=1,3
        if (transient) mass(i,j)=-test(i)*depth*basis(j)
        stiffness(i,j)=-dt*test(i)*(convective_derivative(j)-reaction*basis(j))
      end do
      rhs(i)=-dt*test(i)*source
    end do
  end subroutine supg_terms

  ! Residual viscosity [length**2/time], capped by first-order advective diffusion.
  ! h is measured along grad(C), not the long edge of an anisotropic triangle.
  ! A zero gradient or zero flow gives zero viscosity; no concentration clipping.
  pure function shock_diffusivity(velocity, gradients, gradient, residual, factor) result(nu)
    real(kind=rkind), intent(in) :: velocity(2), gradients(3,2), gradient(2), residual, factor
    real(kind=rkind) :: nu, magnitude, speed, inverse_h, h
    nu=0
    magnitude=norm2(gradient); speed=norm2(velocity)
    if (magnitude<=tiny(1.0_rkind) .or. speed<=0 .or. factor<=0) return
    inverse_h=sum(abs(matmul(gradients,gradient/magnitude)))
    if (inverse_h<=0) return
    h=2/inverse_h
    nu=0.5_rkind*h*min(speed,factor*(abs(residual)/magnitude))
  end function shock_diffusivity

  ! Lag the viscosity at the current Picard iterate, diffuse the new unknown.
  ! Broken P1 residual shares the coefficient-jump limitation of SUPG.
  pure subroutine shock_terms(q,depth,gradients,basis,current,previous,dt,transient, &
                             reaction,source,factor,stiffness)
    real(kind=rkind), intent(in) :: q(2),depth,gradients(3,2),basis(3),current(3),previous(3)
    real(kind=rkind), intent(in) :: dt,reaction,source,factor
    logical, intent(in) :: transient
    real(kind=rkind), intent(out) :: stiffness(3,3)
    real(kind=rkind) :: gradient(2),residual,nu
    stiffness=0
    if (depth<=0 .or. factor<=0 .or. dt<=0) return
    gradient=matmul(transpose(gradients),current-current(1))
    residual=(dot_product(q,gradient)-reaction*dot_product(basis,current)-source)/depth
    if (transient) residual=residual+dot_product(basis,current-previous)/dt
    nu=shock_diffusivity(q/depth,gradients,gradient,residual,factor)
    stiffness=-dt*depth*nu*matmul(gradients,transpose(gradients))
  end subroutine shock_terms

  subroutine adenc_supg_element(pde_loc,el_id,dt,quadpnt_in)
    use global_objs
    use globals
    use pde_objs
    use ncglobvars, only: LSsupg, LSsupg_factor, LSshock, LSshock_factor
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el_id
    real(kind=rkind), intent(in) :: dt
    type(integpnt_str), intent(in), optional :: quadpnt_in
    type(integpnt_str) :: point, nodepoint
    integer :: l,j
    integer(kind=ikind) :: layer
    real(kind=rkind) :: q(2),tensor(2,2),depth,reaction,source,weight
    real(kind=rkind) :: temporal(3,3),spatial(3,3),forcing(3)
    real(kind=rkind) :: current(3),shock(3,3),supg_factor
    logical :: transient
    if (.not. ((LSsupg .and. LSsupg_factor>0) .or. (LSshock .and. LSshock_factor>0))) return
    if (size(pde)/=1 .or. drutes_config%dimen/=2 .or. size(stiff_mat,1)/=3) &
      error stop 'ADEnc SUPG requires a single 2D P1 equation'
    if (.not. pde_loc%diffusion) error stop 'ADEnc SUPG requires Galerkin diffusion assembly'
    if (pde_common%timeint_method<0 .or. pde_common%timeint_method>2) &
      error stop 'ADEnc SUPG supports steady state or implicit Euler only'
    transient=pde_common%timeint_method/=0
    if (dt<=0) error stop 'ADEnc SUPG requires positive assembly dt'
    if (present(quadpnt_in)) point=quadpnt_in
    point%element=el_id
    point%type_pnt='gqnd'
    point%column=2
    supg_factor=0
    if (LSsupg) supg_factor=LSsupg_factor
    if (LSshock .and. LSshock_factor>0) then
      nodepoint=point
      nodepoint%type_pnt='ndpt'
      do j=1,3
        nodepoint%order=elements%data(el_id,j)
        current(j)=pde_loc%getval(nodepoint)
      end do
    end if
    layer=elements%material(el_id)
    do l=1,size(gauss_points%weight)
      point%order=l
      call pde_loc%pde_fnc(1)%convection(pde_loc,layer,point,vector_out=q)
      if (norm2(q)<=0) cycle
      depth=pde_loc%pde_fnc(1)%elasticity(pde_loc,layer,point)
      call pde_loc%pde_fnc(1)%dispersion(pde_loc,layer,point,tensor=tensor)
      reaction=pde_loc%pde_fnc(1)%reaction(pde_loc,layer,point)
      source=pde_loc%pde_fnc(1)%zerord(pde_loc,layer,point)
      call supg_terms(q,depth,tensor,elements%ders(el_id,:,1:2),base_fnc(:,l), &
        dt,transient,reaction,source,supg_factor,temporal,spatial,forcing)
      if (LSshock .and. LSshock_factor>0) then
        call shock_terms(q,depth,elements%ders(el_id,:,1:2),base_fnc(:,l),current,elnode_prev, &
          dt,transient,reaction,source,LSshock_factor,shock)
        spatial=spatial+shock
      end if
      weight=gauss_points%weight(l)*elements%areas(el_id)/gauss_points%area
      cap_mat=cap_mat+weight*temporal
      stiff_mat=stiff_mat+weight*spatial
      bside=bside+weight*forcing
      if (transient) bside=bside+weight*matmul(temporal,elnode_prev)
    end do
  end subroutine adenc_supg_element
end module ncsupg
