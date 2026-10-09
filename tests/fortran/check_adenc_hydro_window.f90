! Production hydraulic forcing-window test only. NEVER main or a transport solve.
! Run in a fresh directory with enabled hydroflow and the actual benchmark inputs.
program check_adenc_hydro_window
  use typy
  use globals
  use pde_objs
  use read_inputs, only: read_global,read_2dmesh_gmsh
  use init_netcdf, only: initialize_nc=>netcdf
  use ncglobvars
  use nchydroflow
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: ends(2,3)=reshape([1,2,2,3,3,1],[2,3])
  integer :: e,k,a,b,u,step
  real(kind=rkind), allocatable :: old_h(:),volume(:)
  real(kind=rkind) :: t,dt,xy(2),q(2),h,divq,edge_flux(3),balance,scale,max_error,min_out,total_out
  real(kind=rkind) :: vertices(3,2),residual_div
  real(kind=rkind) :: max_absolute,water_scale,water,solute,total_water,total_solute
  open(newunit=file_global,file='drutes.conf/global.conf',status='old',action='read')
  call read_global(); close(file_global)
  open(newunit=file_mesh,file='drutes.conf/mesh/mesh.msh',status='old',action='read')
  call read_2dmesh_gmsh(); close(file_mesh)
  allocate(pde(1)); pde_common%processes=1
  call initialize_nc()
  if (.not.LShydro) error stop 'Window preflight requires enabled hydroflow'
  allocate(old_h(elements%kolik),volume(elements%kolik)); old_h=0; volume=0
  do e=1,elements%kolik
    if (.not.ncfluxdata%activeel(e)) cycle
    vertices=nodes%data(elements%data(e,:),1:2)
    volume(e)=abs((vertices(2,1)-vertices(1,1))*(vertices(3,2)-vertices(1,2))- &
                  (vertices(2,2)-vertices(1,2))*(vertices(3,1)-vertices(1,1)))/2
    xy=sum(vertices,dim=1)/3
    call hydro_value(e,xy,q,old_h(e),divq)
  end do
  open(newunit=u,file='hydro-window.txt',status='new',action='write')
  write(u,'(A)') '# time[s] min_outlet[m3/s] total_outlet[m3/s] max_local_scaled_residual global_scaled_residual ' // &
    'lateral_water[m3/s] lateral_solute[concentration*m3/s]'
  t=0; step=0
  do while(t<end_time)
    dt=min(1800.0_rkind,end_time-t)
    call hydro_clip_step(t,dt)
    call hydro_begin(t+dt,dt,t)
    max_error=0; max_absolute=0; water_scale=1; min_out=huge(1.0_rkind); total_out=0
    total_water=0; total_solute=0
    do e=1,elements%kolik
      if (.not.ncfluxdata%activeel(e)) cycle
      vertices=nodes%data(elements%data(e,:),1:2); xy=sum(vertices,dim=1)/3
      call hydro_value(e,xy,q,h,divq)
      call hydro_sources(e,water,solute)
      total_water=total_water+volume(e)*water; total_solute=total_solute+volume(e)*solute
      if (any(.not.ieee_is_finite(q)) .or. .not.ieee_is_finite(h) .or. h<=0) &
        error stop 'Nonfinite/invalid hydraulic coefficients'
      do k=1,3
        a=ends(1,k); b=ends(2,k)
        edge_flux(k)=hydro_normal_flux(e,k)*norm2(vertices(b,:)-vertices(a,:))
        if (open_edges(k,e)) then
          min_out=min(min_out,edge_flux(k)); total_out=total_out+edge_flux(k)
          if (edge_flux(k)<0) error stop 'Outlet reversed in full-window test'
        end if
        if (bank_edges(k,e) .and. abs(edge_flux(k))>1.e-10_rkind) error stop 'Water crossed bank'
      end do
      balance=volume(e)*((h-old_h(e))/dt-water)+sum(edge_flux)
      scale=max(1.0_rkind,maxval(abs(edge_flux)))
      residual_div=volume(e)*divq-sum(edge_flux)
      max_error=max(max_error,abs(balance)/scale,abs(residual_div)/scale)
      max_absolute=max(max_absolute,abs(balance),abs(residual_div))
      water_scale=max(water_scale,scale,abs(volume(e)*(h-old_h(e))/dt),abs(volume(e)*water))
      old_h(e)=h
    end do
    ! PCG uses a global flux/target scale, NOT each tiny element's own flux.
    ! Match its configured Rhine rtol1e-10 acceptance (10*rtol), and also
    ! report the more stringent per-element normalization without hiding it.
    if (max_absolute/water_scale>1.e-9_rkind) then
      print *, 'Water residual absolute / global scale / local normalized:',max_absolute,water_scale,max_error
      error stop 'Local water/storage balance lost in window test'
    end if
    call hydro_end(.true.)
    t=t+dt; step=step+1
    write(u,'(7(ES24.16,1X))') t,min_out,total_out,max_error,max_absolute/water_scale,total_water,total_solute
    if (mod(step,48)==0) print *, 'Hydraulic window passed through day=',t/86400
  end do
  close(u)
  print *, 'Full hydraulic forcing window passed; end time=',t,' trials=',step
end program
