! Residual-based SUPG for the depth-integrated, P1, 2D ADEnc equation.
module ncsupg
  use typy
  implicit none
  private
  public :: read_adenc_supg, adenc_supg_element, supg_tau, supg_terms, shock_diffusivity, shock_terms
  public :: rt0_dispersion_divergence
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
  ! Legacy strong residual: h C_t + q.grad(C) - reaction*C - source.
  ! Optional RT0 conservative mode adds div(q)*C - div(K).grad(C).
  ! P1 Hessian is zero. Legacy hydro-cell/interface jumps remain omitted.
  pure subroutine supg_terms(q, depth, tensor, gradients, basis, dt, transient, &
                            reaction, source, factor, mass, stiffness, rhs,divq,divtensor)
    real(kind=rkind), intent(in) :: q(2), depth, tensor(2,2), gradients(3,2), basis(3)
    real(kind=rkind), intent(in) :: dt, reaction, source, factor
    logical, intent(in) :: transient
    real(kind=rkind), intent(out) :: mass(3,3), stiffness(3,3), rhs(3)
    real(kind=rkind), intent(in), optional :: divq,divtensor(2)
    real(kind=rkind) :: tau, test(3), convective_derivative(3),rate,drift(2)
    integer :: i,j
    mass=0; stiffness=0; rhs=0
    if (depth<=0 .or. factor<=0) return
    rate=reaction; drift=q
    if (present(divq)) rate=rate-divq
    if (present(divtensor)) drift=drift-divtensor
    tau=factor*supg_tau(q/depth,tensor/depth,gradients,dt,transient,rate/depth)
    test=tau*matmul(gradients,q/depth)
    convective_derivative=matmul(gradients,drift)
    do i=1,3
      do j=1,3
        if (transient) mass(i,j)=-test(i)*depth*basis(j)
        stiffness(i,j)=-dt*test(i)*(convective_derivative(j)-rate*basis(j))
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
  ! Broken P1 residual: same optional RT0 divergence terms as SUPG.
  ! Does not add interelement diffusive jump residuals or a positivity limiter.
  pure subroutine shock_terms(q,depth,gradients,basis,current,previous,dt,transient, &
                             reaction,source,factor,stiffness,old_depth,divq,divtensor)
    real(kind=rkind), intent(in) :: q(2),depth,gradients(3,2),basis(3),current(3),previous(3)
    real(kind=rkind), intent(in) :: dt,reaction,source,factor
    logical, intent(in) :: transient
    real(kind=rkind), intent(out) :: stiffness(3,3)
    real(kind=rkind), intent(in), optional :: old_depth
    real(kind=rkind), intent(in), optional :: divq,divtensor(2)
    real(kind=rkind) :: gradient(2),residual,nu,rate,drift(2)
    stiffness=0
    if (depth<=0 .or. factor<=0 .or. dt<=0) return
    gradient=matmul(transpose(gradients),current-current(1))
    rate=reaction; drift=q
    if (present(divq)) rate=rate-divq
    if (present(divtensor)) drift=drift-divtensor
    residual=(dot_product(drift,gradient)-rate*dot_product(basis,current)-source)/depth
    if (transient) residual=residual+dot_product(basis,current-previous)/dt
    if (transient .and. present(old_depth)) &
      residual=residual+(1-old_depth/depth)*dot_product(basis,previous)/dt
    nu=shock_diffusivity(q/depth,gradients,gradient,residual,factor)
    stiffness=-dt*depth*nu*matmul(gradients,transpose(gradients))
  end subroutine shock_terms

  ! RT0 has grad(q)=b*I, div(q)=2*b. For the physical ADEnc tensor
  ! K=alpha_T*|q|*I+(alpha_L-alpha_T)*q*q^T/|q| this is exact away from q=0.
  pure function rt0_dispersion_divergence(q,divq,alpha_l,alpha_t) result(value)
    real(kind=rkind), intent(in) :: q(2),divq,alpha_l,alpha_t
    real(kind=rkind) :: value(2),speed
    value=0; speed=norm2(q)
    if (speed>tiny(1.0_rkind)) value=.5_rkind*divq*(2*alpha_l-alpha_t)*q/speed
  end function

  subroutine adenc_supg_element(pde_loc,el_id,dt,quadpnt_in)
    use global_objs
    use globals
    use pde_objs
    use ncglobvars, only: LSsupg, LSsupg_factor, LSshock, LSshock_factor, LSconservative,LSdepth_old
    use ncglobvars, only: LSdisp,LSdisp_transverse
    use nchydroflow, only: LShydro,hydro_value,hydro_sources
    use geom_tools, only: getcoor
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el_id
    real(kind=rkind), intent(in) :: dt
    type(integpnt_str), intent(in), optional :: quadpnt_in
    type(integpnt_str) :: point, nodepoint
    integer :: l,j
    integer(kind=ikind) :: layer
    real(kind=rkind) :: q(2),tensor(2,2),depth,reaction,source,weight,divq,divtensor(2),xy(2)
    real(kind=rkind) :: temporal(3,3),spatial(3,3),forcing(3)
    real(kind=rkind) :: current(3),shock(3,3),supg_factor
    real(kind=rkind) :: lateral_water,lateral_load
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
      call hydro_sources(int(el_id),lateral_water,lateral_load)
      source=source+lateral_load
      divq=0; divtensor=0
      if (LShydro) then
        call getcoor(point,xy)
        call hydro_value(int(el_id),xy,q,depth,divq)
        divtensor=rt0_dispersion_divergence(q,divq,LSdisp,LSdisp_transverse)
      end if
      call supg_terms(q,depth,tensor,elements%ders(el_id,:,1:2),base_fnc(:,l), &
        dt,transient,reaction,source,supg_factor,temporal,spatial,forcing,divq,divtensor)
      if (LSshock .and. LSshock_factor>0) then
        if (LSconservative) then
          call shock_terms(q,depth,elements%ders(el_id,:,1:2),base_fnc(:,l),current,elnode_prev, &
            dt,transient,reaction,source,LSshock_factor,shock,LSdepth_old(l,el_id),divq,divtensor)
        else
          call shock_terms(q,depth,elements%ders(el_id,:,1:2),base_fnc(:,l),current,elnode_prev, &
            dt,transient,reaction,source,LSshock_factor,shock)
        end if
        spatial=spatial+shock
      end if
      weight=gauss_points%weight(l)*elements%areas(el_id)/gauss_points%area
      cap_mat=cap_mat+weight*temporal
      stiff_mat=stiff_mat+weight*spatial
      bside=bside+weight*forcing
      if (transient) then
        if (LSconservative) then
          bside=bside+weight*matmul(temporal,elnode_prev)*LSdepth_old(l,el_id)/depth
        else
          bside=bside+weight*matmul(temporal,elnode_prev)
        end if
      end if
    end do
  end subroutine adenc_supg_element
end module ncsupg
