! Read/query test only; NEVER DRUtES main or a transport solve.
program test_ncrouting
  use typy
  use netcdf
  use ncglobvars
  use ncrouting
  use netcdfflux, only: ncflux_init
  implicit none
  character(len=1024) :: path,mode,message,forcing
  logical :: ok,found
  integer :: status,unit,i,j,heads,sinks,branches,links
  integer(kind=ikind) :: cell,target,code,other_code
  integer(kind=ikind), allocatable :: upstream(:),back(:)
  real(kind=rkind) :: direction(2),other(2)
  call get_command_argument(1,path)
  call get_command_argument(2,mode)
  if (trim(mode)=='real') then
    call get_command_argument(3,forcing)
    if (len_trim(forcing)==0) forcing='drutes.conf/netcdf/mRM_Fluxes_States.nc'
    status=nf90_open(trim(forcing),nf90_nowrite,netcdfID)
    if (status/=nf90_noerr) error stop 'cannot open real Qrouted'
    call ncflux_init(ok,message)
    if (.not.ok) error stop 'real Qrouted initialization failed'
  else
    allocate(ncfluxdata%lat(2),ncfluxdata%lon(2))
    ncfluxdata%lat=[51.0_rkind,50.875_rkind]; ncfluxdata%lon=[7.0_rkind,7.125_rkind]
    ncfluxdata%nlat=2; ncfluxdata%nlon=2
  end if
  if (trim(mode)=='required') then
    call read_adenc_routing()
    print *, 'Mandatory routing read passed'
    stop
  end if
  call routing_read(trim(path),ok,message)
  if (trim(mode)=='invalid') then
    if (ok) error stop 'invalid routing file accepted'
    call routing_cell(1_ikind,target,upstream,direction,code,found)
    if (found .or. size(upstream)/=0 .or. routing_enabled) error stop 'partial graph leaked'
    print *, trim(message)
    stop
  end if
  if (.not.ok) then
    print *, trim(message)
    error stop 'routing read failed'
  end if
  if (trim(mode)=='real') then
    open(newunit=unit,file='routing-analysis.csv',status='new',action='write')
    write(unit,'(A)') 'cell,downstream_cell,upstream_count,raw_fdir,unit_x,unit_y'
    heads=0; sinks=0; branches=0; links=0
    do cell=1,ncfluxdata%nlat*ncfluxdata%nlon
      call routing_cell(cell,target,upstream,direction,code,found)
      if (.not.found) cycle
      if (target==0) then
        sinks=sinks+1
        if (any(direction/=0)) error stop 'invented sink direction'
      else
        links=links+1
        if (abs(norm2(direction)-1)>1.e-12_rkind) error stop 'non-unit link'
        call routing_cell(target,code,back,other,other_code,found)
        if (.not.found .or. .not.any(back==cell)) error stop 'inconsistent downstream/upstream'
      end if
      if (size(upstream)==0) heads=heads+1
      if (size(upstream)>1) branches=branches+1
      call routing_cell(cell,target,upstream,direction,code,found)
      write(unit,'(I0,3(",",I0),2(",",ES24.16))') cell,target,size(upstream),code,direction
    end do
    close(unit)
    print *, 'Real routing checks passed: links, sinks, headwaters, confluences=',links,sinks,heads,branches
    stop
  end if
  call routing_cell(1_ikind,target,upstream,direction,code,ok)
  if (.not.ok .or. target/=2 .or. size(upstream)/=0 .or. direction(1)<.99_rkind) &
    error stop 'incorrect eastward headwater mapping'
  call routing_cell(2_ikind,target,upstream,direction,code,ok)
  if (.not.ok .or. target/=4 .or. direction(2)>-.99_rkind) error stop 'incorrect southward link'
  call routing_cell(4_ikind,target,upstream,direction,code,ok)
  if (.not.ok .or. target/=0 .or. size(upstream)/=2 .or. any(direction/=0)) error stop 'incorrect outlet'
  if (.not.all(upstream==[2_ikind,3_ikind])) error stop 'incorrect confluence'
  allocate(el2ncgrid(3)); el2ncgrid=[3_ikind,1_ikind,-1_ikind]
  call routing_element(1_ikind,target,upstream,direction,code,ok)
  if (.not.ok .or. target/=4 .or. direction(1)<.99_rkind) error stop 'incorrect FE mapping'
  other=routing_flow_direction(1_ikind,found)
  if (.not.found .or. norm2(other-direction)>1.e-12_rkind) error stop 'incorrect direction function'
  call routing_element(3_ikind,target,upstream,direction,code,ok)
  if (ok) error stop 'unmapped FE accepted'
  call routing_element(0_ikind,target,upstream,direction,code,ok)
  if (ok) error stop 'invalid FE accepted'
  allocate(ncfluxdata%activeel(3),ncfluxdata%fluxvct(3,2))
  ncfluxdata%activeel=.true.; ncfluxdata%fluxvct=0
  call routing_apply_directions(ok,message)
  if (.not.ok .or. ncfluxdata%fluxvct(1,1)<.99_rkind) error stop 'directions not applied'
  call routing_close()
  call routing_cell(1_ikind,target,upstream,direction,code,ok)
  if (ok) error stop 'closed cache accepted'
  ! Repeated reads must not keep stale topology or leak allocatable state.
  do j=1,3
    call routing_read(trim(path),ok,message)
    if (.not.ok) error stop 'repeated read failed'
  end do
  print *, 'Synthetic routing checks passed'
end program
