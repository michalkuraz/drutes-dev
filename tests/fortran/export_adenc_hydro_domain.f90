! Diagnostic only: export the centroid-filtered Rhine domain before choosing ports.
! Run in a fresh directory with hydroflow absent/off. Does not solve or run main.
program export_adenc_hydro_domain
  use typy
  use globals
  use pde_objs
  use read_inputs, only: read_global,read_2dmesh_gmsh
  use init_netcdf, only: initialize_nc=>netcdf
  use ncglobvars
  use ncfluxarea, only: ncflux_active_width
  use netcdfflux, only: ncflux_get_xy_cell
  use nchydroflow, only: LShydro,hydro_filter
  implicit none
  integer :: e,j,u,day,days,status
  integer(kind=ikind), allocatable :: original(:)
  real(kind=rkind) :: xy(2),discharge
  logical :: ok
  character(len=1024) :: message
  character(len=32) :: argument
  days=1
  call get_command_argument(1,argument)
  if (len_trim(argument)>0) then
    read(argument,*,iostat=status) days
    if (status/=0 .or. days<1 .or. days>100) error stop 'Diagnostic span must be 1..100 days'
  end if
  open(newunit=file_global,file='drutes.conf/global.conf',status='old',action='read')
  call read_global(); close(file_global)
  open(newunit=file_mesh,file='drutes.conf/mesh/mesh.msh',status='old',action='read')
  call read_2dmesh_gmsh(); close(file_mesh)
  original=nodes%edge
  allocate(pde(1)); pde_common%processes=1
  call initialize_nc()
  if (LShydro) error stop 'domain discovery requires hydroflow absent/off'
  LShydro=.true.
  call hydro_filter()
  open(newunit=u,file='nodes.txt',status='new',action='write')
  do j=1,nodes%kolik
    write(u,'(I0,1X,2(ES24.16,1X),I0)') j,nodes%data(j,1:2),original(j)
  end do
  close(u)
  open(newunit=u,file='active.txt',status='new',action='write')
  do e=1,elements%kolik
    if (.not.ncfluxdata%activeel(e)) cycle
    xy=sum(nodes%data(elements%data(e,:),1:2),dim=1)/3
    call ncflux_get_xy_cell(xy(1),xy(2),ora_di_ini,discharge,ok,message)
    if (.not.ok) error stop 'invalid active centroid after hydro_filter'
    write(u,'(4(I0,1X),6(ES24.16,1X))') e,elements%data(e,:),xy, &
      ncflux_active_width(int(e,ikind)),discharge,ncfluxdata%fluxvct(e,:)
  end do
  close(u)
  open(newunit=u,file='next.txt',status='new',action='write')
  do e=1,elements%kolik
    if (.not.ncfluxdata%activeel(e)) cycle
    xy=sum(nodes%data(elements%data(e,:),1:2),dim=1)/3
    call ncflux_get_xy_cell(xy(1),xy(2),ora_di_ini+24_ikind,discharge,ok,message)
    if (.not.ok) error stop 'missing next daily forcing slice'
    write(u,'(I0,1X,ES24.16)') e,discharge
  end do
  close(u)
  open(newunit=u,file='window.txt',status='new',action='write')
  do day=0,days
    do e=1,elements%kolik
      if (.not.ncfluxdata%activeel(e)) cycle
      xy=sum(nodes%data(elements%data(e,:),1:2),dim=1)/3
      call ncflux_get_xy_cell(xy(1),xy(2),ora_di_ini+int(24*day,ikind),discharge,ok,message)
      if (.not.ok) error stop 'missing diagnostic daily forcing slice'
      write(u,'(2(I0,1X),ES24.16)') day,e,discharge
    end do
  end do
  close(u)
  print *, 'Diagnostic domain exported, active triangles=',count(ncfluxdata%activeel)
end program
