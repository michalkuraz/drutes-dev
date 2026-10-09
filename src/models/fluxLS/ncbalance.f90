! ADEnc inventory and discrete boundary reactions from UNSCALED local equations.
module ncbalance
  use typy
  implicit none
  private
  public :: balance_capture,balance_finish,balance_reset
  real(kind=rkind), allocatable :: matrices(:,:,:),rhs(:,:),boundary(:,:,:)
  real(kind=rkind), allocatable :: old_inventory(:),new_weights(:,:),source_amount(:)
  logical, allocatable :: seen(:)
  real(kind=rkind), public :: balance_error=0, balance_mass=0, balance_free_residual=0
  real(kind=rkind) :: initial_mass=0,cum_in=0,cum_out=0,cum_source=0
  integer :: output_unit=0
  logical :: started=.false.
contains
  subroutine balance_reset()
    if (output_unit/=0) close(output_unit)
    output_unit=0; started=.false.
    cum_in=0; cum_out=0; cum_source=0; balance_error=0; balance_mass=0; balance_free_residual=0
    if (allocated(seen)) deallocate(seen,matrices,rhs,boundary,old_inventory,new_weights,source_amount)
  end subroutine

  subroutine balance_capture(el,oldcap,newcap,open_matrix,source)
    use globals, only: elements,stiff_mat,cap_mat,bside,elnode_prev
    integer(kind=ikind), intent(in) :: el
    real(kind=rkind), intent(in) :: oldcap(3,3),newcap(3,3),open_matrix(3,3),source
    if (.not. allocated(seen)) then
      allocate(seen(elements%kolik),matrices(3,3,elements%kolik),rhs(3,elements%kolik), &
        boundary(3,3,elements%kolik),old_inventory(elements%kolik),new_weights(3,elements%kolik), &
        source_amount(elements%kolik))
      seen=.false.
    end if
    ! Replaced on every Picard assembly, never accumulated over iterations.
    matrices(:,:,el)=stiff_mat+cap_mat
    rhs(:,el)=bside
    boundary(:,:,el)=open_matrix
    old_inventory(el)=-sum(matmul(oldcap,elnode_prev))
    new_weights(:,el)=-sum(newcap,dim=1)
    source_amount(el)=source
    seen(el)=.true.
  end subroutine

  subroutine balance_finish(pde_loc,accepted,t,dt)
    use pde_objs, only: pde_str,pde_common
    use globals, only: elements,nodes
    use ncglobvars, only: ncfluxdata
    class(pde_str), intent(in) :: pde_loc
    logical, intent(in) :: accepted
    real(kind=rkind), intent(in) :: t,dt
    real(kind=rkind) :: c(3),residual(3),oldmass,net,source,step_in,step_out,negative,scale,flux
    real(kind=rkind), allocatable :: reaction(:)
    integer(kind=ikind) :: el,j,node,k
    integer :: ierr
    if (.not. accepted) then
      if (allocated(seen)) seen=.false.
      return
    end if
    if (.not. allocated(seen)) error stop 'No ADEnc balance assembly captured'
    if (any(ncfluxdata%activeel .and. .not. seen)) error stop 'Incomplete ADEnc balance assembly'
    allocate(reaction(nodes%kolik)); reaction=0
    balance_mass=0; oldmass=0; source=0; net=0; negative=0; balance_free_residual=0
    do el=1,elements%kolik
      if (.not. ncfluxdata%activeel(el)) cycle
      do j=1,3
        node=elements%data(el,j); k=pde_loc%permut(node)
        if (k>0) then
          c(j)=pde_common%xvect(k,3)
        else if (associated(pde_loc%boundary_history)) then
          c(j)=pde_loc%boundary_history(el,j,3_ikind)
        else
          call pde_loc%bc(nodes%edge(node))%value_fnc(pde_loc,el,j,c(j))
        end if
      end do
      balance_mass=balance_mass+dot_product(new_weights(:,el),c)
      negative=negative+dot_product(new_weights(:,el),max(-c,0.0_rkind))
      oldmass=oldmass+old_inventory(el); source=source+source_amount(el)
      net=net-sum(matmul(boundary(:,:,el),c)) ! outward amount over this step
      residual=matmul(matrices(:,:,el),c)-rhs(:,el)
      do j=1,3
        node=elements%data(el,j)
        reaction(node)=reaction(node)+residual(j)
      end do
    end do
    step_in=max(-net,0.0_rkind); step_out=max(net,0.0_rkind)
    do node=1,nodes%kolik
      if (pde_loc%permut(node)>0) then
        balance_free_residual=max(balance_free_residual,abs(reaction(node)))
      else
        ! Negative equation convention: reaction = dt * outward Dirichlet flux.
        flux=reaction(node)
        step_in=step_in+max(-flux,0.0_rkind); step_out=step_out+max(flux,0.0_rkind)
      end if
    end do
    if (.not. started) then
      initial_mass=oldmass
      open(newunit=output_unit,file='out/adenc_mass_balance.csv',status='replace',action='write',iostat=ierr)
      if (ierr/=0) error stop 'Cannot open ADEnc mass balance output'
      write(output_unit,'(a)') 'time_s,dt_s,inventory,initial_inventory,cumulative_in,cumulative_out,' // &
        'cumulative_source,error,relative_error,step_error,free_residual,negative_nodal_inventory'
      started=.true.
    end if
    cum_in=cum_in+step_in; cum_out=cum_out+step_out; cum_source=cum_source+source
    balance_error=balance_mass-initial_mass-cum_in+cum_out-cum_source
    scale=max(abs(initial_mass)+cum_in+abs(cum_source),tiny(1.0_rkind))
    write(output_unit,'(*(es24.16,:,","))') t+dt,dt,balance_mass,initial_mass,cum_in,cum_out,cum_source, &
      balance_error,balance_error/scale,balance_mass-oldmass-step_in+step_out-source, &
      balance_free_residual,negative
    flush(output_unit)
    seen=.false.
  end subroutine
end module
