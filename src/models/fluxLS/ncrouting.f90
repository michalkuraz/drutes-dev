! Static mRM L11 routing graph. Use explicit links, not Q ordering or DEM slopes.
! Restart dimension names/orientation can be misleading: map by lat/lon centres.
module ncrouting
  use typy
  use netcdf
  use ncglobvars, only: ncfluxdata, el2ncgrid
  use ncmesh, only: latlong2utm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: read_adenc_routing, routing_read, routing_cell, routing_element
  public :: routing_flow_direction
  public :: routing_apply_directions, routing_close, routing_enabled, routing_path
  character(len=*), parameter :: routing_path='drutes.conf/netcdf/mRM_restart_001.nc'
  logical, protected :: routing_enabled=.false.
  logical :: ready=.false.
  integer, allocatable :: cell_node(:), node_cell(:), downstream(:), codes(:), upstream_count(:)
  real(kind=rkind), allocatable :: centres(:,:)
contains
  ! Required for every ADEnc initialization. Single-domain fixed filename.
  subroutine read_adenc_routing()
    use core_tools, only: write_log
    logical :: exists, ok
    character(len=1024) :: message
    call routing_close()
    inquire(file=routing_path,exist=exists)
    if (.not.exists) error stop 'ADEnc requires drutes.conf/netcdf/mRM_restart_001.nc (mHM routing restart)'
    call routing_read(routing_path,ok,message)
    if (.not.ok) then
      print *, trim(message)
      error stop 'Invalid ADEnc routing network'
    end if
    write(message,'(A,I0,A,I0,A,I0)') 'ADEnc mRM routing: nodes=',size(downstream), &
      ', outlets=',count(downstream==0),', confluences=',count(upstream_count>1)
    call write_log(trim(message))
  end subroutine

  ! Read original restart or its routing-only extract. No dynamic restart state
  ! is restored. Errors invalidate the entire cache (no partial graph exposure).
  subroutine routing_read(path,ok,message)
    character(len=*), intent(in) :: path
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: nc, status, close_status, vid, nd, dims(2), nx, ny, fill
    integer :: i,j,k,id,n,ilat,ilon,c,target,source,head,tail,seen
    integer, allocatable :: ids(:,:), mask(:,:), fdir(:,:), from(:), to(:), degree(:), queue(:)
    real(kind=rkind), allocatable :: lat(:,:),lon(:,:)
    logical, allocatable :: linked(:)
    real(kind=rkind) :: x,y,dx,dy
    real(kind=rkind), parameter :: coordinate_tolerance=1.e-7_rkind
    call routing_close()
    ok=.false.; message='routing_read: invalid network'; nc=-1
    if (.not.allocated(ncfluxdata%lat) .or. .not.allocated(ncfluxdata%lon)) then
      message='routing_read: initialize Qrouted axes first'
      return
    end if
    if (size(ncfluxdata%lat)<1 .or. size(ncfluxdata%lon)<1) return
    if (any(.not.ieee_is_finite(ncfluxdata%lat)) .or. any(.not.ieee_is_finite(ncfluxdata%lon))) then
      message='routing_read: nonfinite Qrouted axes'; return
    end if
    status=nf90_open(path,nf90_nowrite,nc)
    if (status/=nf90_noerr) then
      message='routing_read: cannot open '//trim(path)//': '//trim(nf90_strerror(status))
      return
    end if
    status=nf90_inq_varid(nc,'L11_Id',vid)
    if (status/=nf90_noerr) goto 900
    status=nf90_inquire_variable(nc,vid,ndims=nd)
    if (status/=nf90_noerr) goto 900
    if (nd/=2) then
      message='routing_read: L11_Id must be a 2D grid'; goto 900
    end if
    status=nf90_inquire_variable(nc,vid,dimids=dims)
    if (status/=nf90_noerr) goto 900
    status=nf90_inquire_dimension(nc,dims(1),len=nx)
    if (status/=nf90_noerr) goto 900
    status=nf90_inquire_dimension(nc,dims(2),len=ny)
    if (status/=nf90_noerr) goto 900
    allocate(ids(nx,ny),mask(nx,ny),fdir(nx,ny),lat(nx,ny),lon(nx,ny))
    status=nf90_get_var(nc,vid,ids)
    if (status/=nf90_noerr) goto 900
    fill=-9999
    close_status=nf90_get_att(nc,vid,'_FillValue',fill)
    if (close_status/=nf90_noerr .and. close_status/=nf90_enotatt) then
      status=close_status; goto 900
    end if
    call int_grid('L11_domain_mask',mask)
    if (status/=nf90_noerr) goto 900
    call int_grid('L11_fDir',fdir)
    if (status/=nf90_noerr) goto 900
    call real_grid('L11_domain_lat',lat)
    if (status/=nf90_noerr) goto 900
    call real_grid('L11_domain_lon',lon)
    if (status/=nf90_noerr) goto 900
    n=count(mask==1)
    if (n<1 .or. any(mask/=0 .and. mask/=1)) then
      message='routing_read: empty or invalid L11_domain_mask'; goto 900
    end if
    allocate(cell_node(size(ncfluxdata%lat)*size(ncfluxdata%lon)),node_cell(n),downstream(n), &
      codes(n),centres(n,2),upstream_count(n),linked(n))
    cell_node=0; node_cell=0; downstream=0; codes=0; upstream_count=0; linked=.false.
    do j=1,ny
      do i=1,nx
        if (mask(i,j)/=1) cycle
        id=ids(i,j)
        if (id<1 .or. id>n) then
          message='routing_read: L11 IDs must cover 1..number of active nodes'; goto 900
        end if
        if (node_cell(id)/=0) then
          message='routing_read: duplicate L11 ID'; goto 900
        end if
        if (.not.ieee_is_finite(lat(i,j)) .or. .not.ieee_is_finite(lon(i,j))) then
          message='routing_read: nonfinite L11 coordinates'; goto 900
        end if
        ilat=minloc(abs(ncfluxdata%lat-lat(i,j)),dim=1)
        ilon=minloc(abs(ncfluxdata%lon-lon(i,j)),dim=1)
        if (abs(ncfluxdata%lat(ilat)-lat(i,j))>coordinate_tolerance .or. &
            abs(ncfluxdata%lon(ilon)-lon(i,j))>coordinate_tolerance) then
          message='routing_read: routing centres do not match Qrouted grid'; goto 900
        end if
        c=(ilat-1)*size(ncfluxdata%lon)+ilon
        if (cell_node(c)/=0) then
          message='routing_read: multiple routing nodes map to one Qrouted cell'; goto 900
        end if
        cell_node(c)=id; node_cell(id)=c; codes(id)=fdir(i,j)
        select case(codes(id))
        case(0,1,2,4,8,16,32,64,128)
        case default
          message='routing_read: invalid L11_fDir code'; goto 900
        end select
        call latlong2utm(lat(i,j),lon(i,j),x,y)
        centres(id,:)=[x,y]
      end do
    end do
    call int_vector('L11_fromN',from)
    if (status/=nf90_noerr) goto 900
    call int_vector('L11_toN',to)
    if (status/=nf90_noerr) goto 900
    if (size(from)/=size(to)) then
      message='routing_read: unequal from/to lengths'; goto 900
    end if
    do k=1,size(from)
      source=from(k); target=to(k)
      ! mRM has a trailing fill record; the sink is represented by fDir=0.
      if (source==fill .and. target==fill) cycle
      if (source<1 .or. source>n .or. target<1 .or. target>n) then
        message='routing_read: link references a missing node'; goto 900
      end if
      if (linked(source)) then
        message='routing_read: duplicate downstream link'; goto 900
      end if
      if (source==target) then
        message='routing_read: self-loop'; goto 900
      end if
      if (codes(source)==0) then
        message='routing_read: outlet has a downstream link'; goto 900
      end if
      c=node_cell(source); id=node_cell(target)
      ilat=(c-1)/size(ncfluxdata%lon); ilon=mod(c-1,size(ncfluxdata%lon))
      if (abs(ilat-(id-1)/size(ncfluxdata%lon))>1 .or. &
          abs(ilon-mod(id-1,size(ncfluxdata%lon)))>1) then
        message='routing_read: downstream cell is not a neighbouring hydrological cell'; goto 900
      end if
      downstream(source)=target; linked(source)=.true.
      upstream_count(target)=upstream_count(target)+1
      dx=centres(target,1)-centres(source,1); dy=centres(target,2)-centres(source,2)
      if (.not.ieee_is_finite(dx) .or. .not.ieee_is_finite(dy) .or. norm2([dx,dy])<=0) then
        message='routing_read: invalid projected link'; goto 900
      end if
    end do
    if (any((downstream==0) .neqv. (codes==0))) then
      message='routing_read: non-outlet node has no downstream link'; goto 900
    end if
    ! Kahn traversal validates the complete directed graph, including tributaries.
    allocate(degree(n),queue(n)); degree=upstream_count; head=1; tail=0; seen=0
    do id=1,n
      if (degree(id)/=0) cycle
      tail=tail+1; queue(tail)=id
    end do
    do while(head<=tail)
      source=queue(head); head=head+1; seen=seen+1; target=downstream(source)
      if (target==0) cycle
      degree(target)=degree(target)-1
      if (degree(target)/=0) cycle
      tail=tail+1; queue(tail)=target
    end do
    if (seen/=n) then
      message='routing_read: cycle in routing graph'; goto 900
    end if
    close_status=nf90_close(nc); nc=-1
    if (close_status/=nf90_noerr) then
      status=close_status; goto 900
    end if
    ready=.true.; routing_enabled=.true.; ok=.true.; message='routing_read: everything ok'
    return
900 continue
    if (status/=nf90_noerr) message='routing_read: NetCDF error: '//trim(nf90_strerror(status))
    if (nc/=-1) close_status=nf90_close(nc)
    call routing_close()
  contains
    subroutine check_grid(name,v)
      character(len=*), intent(in) :: name
      integer, intent(out) :: v
      integer :: d(2),rank
      status=nf90_inq_varid(nc,name,v)
      if (status/=nf90_noerr) return
      status=nf90_inquire_variable(nc,v,ndims=rank)
      if (status/=nf90_noerr) return
      if (rank/=2) then
        status=nf90_einval; return
      end if
      status=nf90_inquire_variable(nc,v,dimids=d)
      if (status/=nf90_noerr) return
      if (any(d/=dims)) status=nf90_einval
    end subroutine
    subroutine int_grid(name,a)
      character(len=*), intent(in) :: name
      integer, intent(out) :: a(:,:)
      integer :: v
      call check_grid(name,v)
      if (status==nf90_noerr) status=nf90_get_var(nc,v,a)
    end subroutine
    subroutine real_grid(name,a)
      character(len=*), intent(in) :: name
      real(kind=rkind), intent(out) :: a(:,:)
      integer :: v
      call check_grid(name,v)
      if (status==nf90_noerr) status=nf90_get_var(nc,v,a)
    end subroutine
    subroutine int_vector(name,a)
      character(len=*), intent(in) :: name
      integer, allocatable, intent(out) :: a(:)
      integer :: v,rank,d(1),length
      status=nf90_inq_varid(nc,name,v)
      if (status/=nf90_noerr) return
      status=nf90_inquire_variable(nc,v,ndims=rank)
      if (status/=nf90_noerr) return
      if (rank/=1) then
        status=nf90_einval; return
      end if
      status=nf90_inquire_variable(nc,v,dimids=d)
      if (status/=nf90_noerr) return
      status=nf90_inquire_dimension(nc,d(1),len=length)
      if (status/=nf90_noerr) return
      allocate(a(length)); status=nf90_get_var(nc,v,a)
    end subroutine
  end subroutine routing_read

  ! Returns DRUtES hydrological cell indices (NOT mRM IDs or FE indices).
  ! Outlet: downstream_cell=0, direction=0, ok=true. No hidden outlet direction.
  ! fdir is the raw mRM code, intentionally not decoded as an assumed GIS D8.
  subroutine routing_cell(cell,downstream_cell,upstream_cells,direction,fdir,ok)
    integer(kind=ikind), intent(in) :: cell
    integer(kind=ikind), intent(out) :: downstream_cell,fdir
    integer(kind=ikind), allocatable, intent(out) :: upstream_cells(:)
    real(kind=rkind), intent(out) :: direction(2)
    logical, intent(out) :: ok
    integer :: id,target,i,k
    downstream_cell=0; fdir=0; direction=0; ok=.false.
    allocate(upstream_cells(0))
    if (.not.ready) return
    if (cell<1 .or. cell>size(cell_node)) return
    id=cell_node(cell)
    if (id==0) return
    fdir=int(codes(id),ikind); target=downstream(id)
    if (target>0) then
      downstream_cell=int(node_cell(target),ikind)
      direction=centres(target,:)-centres(id,:)
      direction=direction/norm2(direction)
    end if
    deallocate(upstream_cells); allocate(upstream_cells(upstream_count(id))); k=0
    do i=1,size(downstream)
      if (downstream(i)/=id) cycle
      k=k+1; upstream_cells(k)=int(node_cell(i),ikind)
    end do
    ok=.true.
  end subroutine

  ! FE array index -> containing hydrological cell -> routing graph.
  subroutine routing_element(element,downstream_cell,upstream_cells,direction,fdir,ok)
    integer(kind=ikind), intent(in) :: element
    integer(kind=ikind), intent(out) :: downstream_cell,fdir
    integer(kind=ikind), allocatable, intent(out) :: upstream_cells(:)
    real(kind=rkind), intent(out) :: direction(2)
    logical, intent(out) :: ok
    downstream_cell=0; fdir=0; direction=0; ok=.false.; allocate(upstream_cells(0))
    if (.not.allocated(el2ncgrid)) return
    if (element<1 .or. element>size(el2ncgrid)) return
    call routing_cell(el2ncgrid(element),downstream_cell,upstream_cells,direction,fdir,ok)
  end subroutine

  ! Convenience function for an FE array index. Always inspect ok: an outlet
  ! legitimately has no graph exit vector, whereas an invalid lookup has ok=false.
  function routing_flow_direction(element,ok) result(direction)
    integer(kind=ikind), intent(in) :: element
    logical, intent(out) :: ok
    real(kind=rkind) :: direction(2)
    integer(kind=ikind) :: target,code
    integer(kind=ikind), allocatable :: upstream(:)
    call routing_element(element,target,upstream,direction,code,ok)
  end function

  ! Apply BEFORE width/hydraulic preparation. At actual mRM sinks retain the
  ! explicitly configured channel direction: a graph ending gives no exit vector.
  subroutine routing_apply_directions(ok,message)
    use core_tools, only: write_log
    integer :: e,nchanged,nsink
    integer(kind=ikind) :: cell,target,code
    integer(kind=ikind), allocatable :: upstream(:)
    real(kind=rkind) :: direction(2)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical :: found
    ok=.false.; message='routing_apply_directions: missing FE mapping/directions'
    if (.not.ready .or. .not.allocated(el2ncgrid)) return
    if (.not.allocated(ncfluxdata%activeel) .or. .not.allocated(ncfluxdata%fluxvct)) return
    if (size(ncfluxdata%activeel)/=size(el2ncgrid) .or. &
        size(ncfluxdata%fluxvct,1)/=size(el2ncgrid)) return
    if (size(ncfluxdata%fluxvct,2)/=2) return
    nchanged=0; nsink=0
    do e=1,size(el2ncgrid)
      if (.not.ncfluxdata%activeel(e)) cycle
      cell=el2ncgrid(e)
      ! Unmapped/masked fringe triangles retain their existing zero-width handling.
      call routing_cell(cell,target,upstream,direction,code,found)
      if (.not.found) cycle
      if (target==0) then
        nsink=nsink+1
        cycle
      end if
      ncfluxdata%fluxvct(e,:)=direction; nchanged=nchanged+1
    end do
    write(message,'(A,I0,A,I0)') 'ADEnc routing directions applied to FE=',nchanged, &
      '; mRM outlet FE retaining channel direction=',nsink
    call write_log(trim(message)); ok=.true.
  end subroutine

  subroutine routing_close()
    ready=.false.; routing_enabled=.false.
    if (allocated(cell_node)) deallocate(cell_node,node_cell,downstream,codes,centres,upstream_count)
  end subroutine
end module ncrouting
