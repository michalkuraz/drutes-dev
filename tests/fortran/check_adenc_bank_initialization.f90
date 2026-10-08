! Production Rhine initialization only. Never invokes main or solve_pde.
program check_adenc_bank_initialization
  use typy
  use globals
  use pde_objs
  use read_inputs, only: read_global,read_2dmesh_gmsh
  use init_netcdf, only: initialize_nc=>netcdf
  use ncglobvars
  use ncfluxarea, only: ncflux_active_width
  implicit none
  integer :: e,k,a,b,zero_widths
  integer, parameter :: ends(2,3)=reshape([1,2,2,3,3,1],[2,3])
  open(newunit=file_global,file='drutes.conf/global.conf',status='old',action='read')
  call read_global(); close(file_global)
  open(newunit=file_mesh,file='drutes.conf/mesh/mesh.msh',status='old',action='read')
  call read_2dmesh_gmsh(); close(file_mesh)
  allocate(pde(1)); pde_common%processes=1
  call initialize_nc()
  if (.not. LSbank_noflow) error stop 'bank preflight requires enabled riverbank.conf'
  if (.not. allocated(pde(1)%assembly_mask)) error stop 'missing assembly mask'
  if (any(pde(1)%assembly_mask .neqv. ncfluxdata%activeel)) error stop 'mask mismatch'
  if (size(pde(1)%bc(101)%series,2)/=2) error stop 'invalid inlet series columns'
  if (pde(1)%bc(addedbc)%code/=1 .or. pde(1)%bc(addedbc)%value/=0) error stop 'unused DOFs not eliminated'
  zero_widths=0
  do e=1,elements%kolik
    if (.not. ncfluxdata%activeel(e)) cycle
    if (ncflux_active_width(int(e,ikind))<=0) zero_widths=zero_widths+1
    do k=1,3
      if (.not. bank_edges(k,e)) cycle
      a=elements%data(e,ends(1,k)); b=elements%data(e,ends(2,k))
      if (nodes%edge(a)==addedbc .or. nodes%edge(b)==addedbc) then
        print *, 'PREFLIGHT invalid bank element/edge/nodes/IDs/addedbc:',e,k,a,b,nodes%edge(a),nodes%edge(b),addedbc
        error stop 'bank retains unused-node Dirichlet'
      end if
    end do
  end do
  print *, 'PREFLIGHT active elements / bank edges / unused nodes / zero widths:', &
    count(ncfluxdata%activeel),count(bank_edges),count(nodes%edge==addedbc),zero_widths
  print *, 'PREFLIGHT inlet nodes / reserved unused ID:',count(nodes%edge==101),addedbc
  print *, 'PREFLIGHT inlet series:',pde(1)%bc(101)%series
  print *, 'Rhine no-flow bank initialization passed (no time stepping)'
end program
