! Impermeable solute banks on the boundary of the fixed active ADEnc FE domain.
! With the existing advective volume form, K grad(C).n = (q.n) C.
! DRUtES uses negative mass/stiffness signs: the Robin edge addition is +dt*q.n*Ni*Nj.
! This removes absorbing banks; it does not repair interior flux jumps or H_t+div(q).
module ncboundary
  use typy
  implicit none
  private
  public :: read_adenc_banks, prepare_adenc_banks, adenc_bank_element, bank_edge_terms
  public :: active_node_element
contains
  subroutine read_adenc_banks()
    use ncglobvars, only: LSbank_noflow
    use globals, only: drutes_config
    use readtools, only: fileread
    use core_tools, only: write_log
    logical :: exists
    integer :: unit, ierr
    LSbank_noflow=.false.
    inquire(file='drutes.conf/netcdf/riverbank.conf',exist=exists)
    if (exists) then
      open(newunit=unit,file='drutes.conf/netcdf/riverbank.conf',status='old',action='read',iostat=ierr)
      if (ierr/=0) error stop 'Cannot open ADEnc riverbank.conf'
      call fileread(LSbank_noflow,unit)
      close(unit)
    end if
    if (.not. LSbank_noflow) return
    if (drutes_config%dimen/=2 .or. drutes_config%it_method/=0) &
      error stop 'ADEnc no-flow banks require 2D and standard Picard (no Schwarz)'
    call write_log('ADEnc banks: zero TOTAL solute flux; active FE domain only')
  end subroutine read_adenc_banks

  subroutine prepare_adenc_banks(original_edge, close_exterior)
    use ncglobvars, only: LSbank_noflow, bank_edges, ncfluxdata, addedbc
    use globals, only: nodes, elements
    use pde_objs, only: pde
    use core_tools, only: write_log
    integer(kind=ikind), intent(in) :: original_edge(:)
    ! Manufactured closed-box tests may additionally close the external mesh.
    ! Production keeps external mesh boundaries separate from internal banks.
    logical, intent(in), optional :: close_exterior
    ! Hash undirected active-triangle edges. A second occurrence is internal.
    integer, allocatable :: head(:), next(:), lo(:), hi(:), owner(:), local_edge(:)
    logical, allocatable :: participating(:), interface_edges(:,:)
    integer :: e, k, a, b, bucket, item, item_count, buckets, nbanks
    integer, parameter :: ends(2,3)=reshape([1,2,2,3,3,1],[2,3])
    character(len=256) :: message
    if (.not. LSbank_noflow) return
    if (size(elements%data,2)/=3) error stop 'No-flow banks require P1 triangles'
    if (.not. any(ncfluxdata%activeel)) error stop 'ADEnc active domain is empty'
    if (allocated(bank_edges)) deallocate(bank_edges)
    allocate(bank_edges(3,elements%kolik)); bank_edges=.false.
    if (allocated(pde(1)%assembly_mask)) deallocate(pde(1)%assembly_mask)
    allocate(pde(1)%assembly_mask(elements%kolik))
    pde(1)%assembly_mask=ncfluxdata%activeel
    allocate(participating(nodes%kolik)); participating=.false.
    buckets=2*elements%kolik+1
    allocate(head(buckets),next(3*elements%kolik),lo(3*elements%kolik),hi(3*elements%kolik), &
             owner(3*elements%kolik),local_edge(3*elements%kolik))
    head=0; item_count=0
    do e=1,elements%kolik
      if (.not. ncfluxdata%activeel(e)) cycle
      participating(elements%data(e,:))=.true.
      do k=1,3
        a=minval(elements%data(e,ends(:,k))); b=maxval(elements%data(e,ends(:,k)))
        bucket=int(modulo(31_8*int(a,8)+int(b,8),int(buckets,8)))+1
        item=head(bucket)
        do while(item/=0)
          if (lo(item)==a .and. hi(item)==b) exit
          item=next(item)
        end do
        if (item==0) then
          item_count=item_count+1; item=item_count
          lo(item)=a; hi(item)=b; owner(item)=e; local_edge(item)=k
          next(item)=head(bucket); head(bucket)=item
          bank_edges(k,e)=.true.
        else
          if (owner(item)==0) error stop 'Nonmanifold ADEnc active edge'
          bank_edges(local_edge(item),owner(item))=.false.
          owner(item)=0
        end if
      end do
    end do
    allocate(interface_edges(3,elements%kolik)); interface_edges=.false.
    do e=1,elements%kolik
      if (ncfluxdata%activeel(e)) cycle
      do k=1,3
        a=minval(elements%data(e,ends(:,k))); b=maxval(elements%data(e,ends(:,k)))
        bucket=int(modulo(31_8*int(a,8)+int(b,8),int(buckets,8)))+1
        item=head(bucket)
        do while(item/=0)
          if (lo(item)==a .and. hi(item)==b) exit
          item=next(item)
        end do
        if (item==0) cycle
        if (owner(item)==0) error stop 'Nonmanifold active/inactive ADEnc interface'
        interface_edges(local_edge(item),owner(item))=.true.
      end do
    end do
    if (present(close_exterior)) then
      if (.not. close_exterior) bank_edges=bank_edges .and. interface_edges
    else
      bank_edges=bank_edges .and. interface_edges
    end if
    ! Matching original physical IDs identify ports, not banks. Do not use the
    ! threshold-generated addedbc here, and do not infer edges from node pairs
    ! inside a triangle (all three nodes may lie on different bank edges).
    do e=1,elements%kolik
      do k=1,3
        if (.not. bank_edges(k,e)) cycle
        a=elements%data(e,ends(1,k)); b=elements%data(e,ends(2,k))
        if (original_edge(a)>100 .and. original_edge(a)==original_edge(b)) bank_edges(k,e)=.false.
      end do
    end do
    do a=1,nodes%kolik
      if (participating(a)) then
        nodes%edge(a)=original_edge(a) ! bank nodes regain DOFs; preserve real ports
      else
        nodes%edge(a)=addedbc ! unused nodes remain eliminated, never assembled
      end if
    end do
    nbanks=count(bank_edges)
    write(message,*) 'ADEnc no-flow domain: elements=',count(ncfluxdata%activeel), &
      ' nodes=',count(participating),' bank edges=',nbanks
    call write_log(trim(message))
  end subroutine prepare_adenc_banks

  function active_node_element(node_id) result(el)
    use globals, only: nodes
    use ncglobvars, only: LSbank_noflow,ncfluxdata
    integer(kind=ikind), intent(in) :: node_id
    integer(kind=ikind) :: el
    integer :: j
    el=nodes%element(node_id)%data(1)
    if (.not. LSbank_noflow) return
    do j=1,size(nodes%element(node_id)%data)
      if (ncfluxdata%activeel(nodes%element(node_id)%data(j))) then
        el=nodes%element(node_id)%data(j)
        return
      end if
    end do
  end function active_node_element

  pure subroutine bank_edge_terms(edge, length, normal_flux, dt, matrix)
    integer, intent(in) :: edge(2)
    real(kind=rkind), intent(in) :: length,normal_flux(2),dt
    real(kind=rkind), intent(out) :: matrix(3,3)
    real(kind=rkind) :: s,basis(3)
    integer :: g,i,j
    matrix=0
    do g=1,2
      s=(1+(2*g-3)/sqrt(3.0_rkind))/2
      basis=0; basis(edge(1))=1-s; basis(edge(2))=s
      do j=1,3
        do i=1,3
          matrix(i,j)=matrix(i,j)+dt*length/2*normal_flux(g)*basis(i)*basis(j)
        end do
      end do
    end do
  end subroutine bank_edge_terms

  subroutine adenc_bank_element(pde_loc, el_id, dt, quadpnt_in)
    use global_objs, only: integpnt_str
    use globals, only: stiff_mat,time,nodes,elements
    use pde_objs, only: pde_str
    use ncglobvars, only: LSbank_noflow,bank_edges,ncfluxdata,ora_di_ini
    use netcdfflux, only: ncflux_get_xy_cell
    use ncfluxarea, only: ncflux_active_width
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: el_id
    real(kind=rkind), intent(in) :: dt
    type(integpnt_str), intent(in), optional :: quadpnt_in
    integer, parameter :: ends(2,3)=reshape([1,2,2,3,3,1],[2,3])
    integer :: k,g
    real(kind=rkind) :: a(2),b(2),center(2),point(2),normal(2),length,s,Q,W,qn(2),matrix(3,3)
    logical :: ok
    character(len=1024) :: message
    if (.not. LSbank_noflow) return
    if (.not. ncfluxdata%activeel(el_id)) return
    W=ncflux_active_width(el_id)
    if (W<=0) return ! matches zero hydrological convection callback
    center=sum(nodes%data(elements%data(el_id,:),1:2),dim=1)/3
    do k=1,3
      if (.not. bank_edges(k,el_id)) cycle
      a=nodes%data(elements%data(el_id,ends(1,k)),1:2)
      b=nodes%data(elements%data(el_id,ends(2,k)),1:2)
      length=norm2(b-a)
      if (length<=0) error stop 'Degenerate ADEnc bank edge'
      normal=[b(2)-a(2),a(1)-b(1)]/length
      if (dot_product(normal,center-(a+b)/2)>0) normal=-normal
      do g=1,2
        s=(1+(2*g-3)/sqrt(3.0_rkind))/2
        point=(1-s)*a+s*b
        ! One-sided active-element trace at hydrological grid boundaries.
        point=point+1.e-8_rkind*(center-point)
        call ncflux_get_xy_cell(point(1),point(2),ora_di_ini+int(time/86400.0_rkind)*24,Q,ok,message)
        qn(g)=0
        if (ok .and. Q>0) qn(g)=Q/W*dot_product(ncfluxdata%fluxvct(el_id,:),normal)
      end do
      call bank_edge_terms(ends(:,k),length,qn,dt,matrix)
      stiff_mat=stiff_mat+matrix
    end do
  end subroutine adenc_bank_element
end module ncboundary
