module decomp_vars
  use typy
  use globals
  use global_objs
  use sparsematrix

  type, private :: cluster
    !> contains ranges of defined material properties
    !! range(1) - low range
    !! range(2) - high range
    !<
    real(kind=rkind), dimension(2)    :: ranges
    integer(kind=ikind), dimension(:), allocatable :: elements
    !>number of elements at cluster
    integer(kind=ikind) :: el_count
    !> NDOF at cluster
    integer(kind=ikind) :: ndof_cl
    !> id of neibourhood elements, (only used for coarse meshes)
    integer(kind=ikind), dimension(:), allocatable :: neighbours
    !> the level is defined as follows
    !! the level is defined as an average exponent value from the nodal hydraulic conductivity
    !<
    integer(kind=ikind)   :: level
    !> if the cluster contains more then two orders of magnitude in conductivity, the cluster is treated saporately, treat_alone = .true., otherwise treat_alone = .false.
    logical :: treat_alone
    !> ID number of the subdomain, the subdomains are obtained from joining the clusters 
    integer(kind=ikind)   :: subdomain
    integer(kind=ikind)   :: prevsub
    !> if evaluated, then the subdomain for this cluster has been already set up
    logical :: evaluated
    !> number of iterations at current cluster
    integer :: iter_count
    !> number of iterations at current subdomain
    integer :: subdomiter
    logical :: active
    !> array of elements at the neighborhood cluster, the number of rows is the number of overlap levels, the lower column rank is zero, where zero is the list of the boundary nodes
    type(smartarray_int), dimension(:), allocatable :: overlap_el
  end type cluster

  type, private :: resvct_str
    real(kind=rkind), dimension(:), allocatable :: main
    real(kind=rkind), dimension(:), allocatable :: ext
  end type resvct_str
  

  type, public :: subdomain_str
    !> local matrix
    type(extsmtx) :: matrix
    !> extension of local matrix on nodes, that lie out of the subdomain but have a graph link to internal nodes in subdomain
    type(extsmtx) :: extmatrix    
    !> dtto b-vector
    real(kind=rkind), dimension(:), allocatable :: extbvect
    !> 
!     integer(kind=ikind), dimension(:), allocatable  
    !> solution in the subdomain
    real(kind=rkind), dimension(:,:), allocatable :: xvect
        !> solution in the subdomain
    real(kind=rkind), dimension(:,:), allocatable :: extxvect
    !> used for subcycling -- solution in the subdomain for previous time level in terms of the global time step
    real(kind=rkind), dimension(:), allocatable :: xvect_init
    !> subdomain b vector
    real(kind=rkind), dimension(:), allocatable :: bvect
    !> corrections of the solution on the subdomain (Schwarz(-Picard)) iteration increment)
    real(kind=rkind), dimension(:), allocatable :: ovect
    !> vector of residual, consist of two vectors, vector formed out from the main subdomain matrix and vector formed out from the extension of the matrix (extension for the exterior nodes that have graph connection with the interior subdomain nodes)
    type(resvct_str) :: resvct
    !> list of elements in the subdomain
    type(smartarray_int) :: elsub
    !> list of coarse elements in the subdomain
    integer(kind=ikind), dimension(:), allocatable :: coarse_el
    !> permut vector - local xvect to global xvect
    integer(kind=ikind), dimension(:), allocatable :: permut
    !> permut vector - extlocal (local nodes that lie out of subdomain but have a graph connection with the subdomain) 
    integer(kind=ikind), dimension(:), allocatable :: extpermut
    !> permut vector global xvect to local subdomain xvect, should be updated into sparse vector later
    integer(kind=ikind), dimension(:), allocatable :: invpermut
    !> permut vector global xvect to local xvect that lies out of subdomain, it is the x node with graph connections to the subdomain, should be updated into sparse vector later
    integer(kind=ikind), dimension(:), allocatable :: extinvpermut
    logical :: solved = .false.
    logical :: poor_res
    logical :: short_dt = .false.
    integer(kind=ikind) :: itcount
    !> NDOF at subdomain
    integer(kind=ikind) :: ndof
    !> numebr of exterior nodes with graph connection to the interior subdomain nodes apart the Dirichlet boundary
    integer(kind=ikind) :: extndof
    integer(kind=ikind) :: pcg_it
    real(kind=rkind) :: last_error
    logical :: critical
    !> an array with elements that should be formally added to the subdomain in order to be able to evaluate residual vector for a single subdomain only
    type(smartarray_int) :: elextra
    !> an array with nodes that should be formally added to the subdomain in order to be able to evaluate residual vector for a single subdomain only, constructed out of 'elextra' array
    type(smartarray_int) :: ndextra
    
    type(smartarray_int) :: disjoint	
    !> an array with elements from the overlap layers
    type(smartarray_int), dimension(:), allocatable :: olayers
    !> an array with nodes from the overlap layers -- global mesh indexes
    type(smartarray_int), dimension(:), allocatable :: ondlayers
    !> list of subdomain boundary clusters
    type(smartarray_int) :: bcluster
    !> time step could be different on different subdomains (subcycling)
    real(kind=rkind) :: time_step
    !> equals to the previous time level of the subdomain
    real(kind=rkind) :: timeprev
    logical :: time_increased
    !> equals to the time level, which is curretnly being processed
    real(kind=rkind) :: time
    real(kind=rkind) :: tmpval
    integer(kind=ikind) :: order
    logical :: finish
      contains 
	procedure :: returnval
  end type subdomain_str

  !> array map with a disjoint map of subdomains, 
  type, public :: disjoint_nd_str
    !> id number of subdomain at particular node with !!permutated!! index
    integer(kind=ikind), dimension(:), allocatable :: data
    !> if certain node is critical 
    class(smartarray_int), dimension(:), allocatable :: crit
  end type disjoint_nd_str

  
  type, private :: ddinfo_str
    integer(kind=ikind) :: number
    real(kind=rkind) :: t_change
    real(kind=rkind) :: cput_change
    type(disjoint_nd_str) :: nd_dsjnt
    integer(kind=ikind) :: ndofs_tot
    integer(kind=ikind) :: overlaps
    integer(kind=ikind), dimension(:), allocatable :: elincoarse
    integer(kind=ikind), dimension(:), allocatable :: coarseinsub
    !> list of nodes (geometrical) with a list of subdomains
    type(smartarray_int), dimension(:), allocatable :: nodesinsub
    !> list of nodes (geometrical) with a list of extended subdomains
    type(smartarray_int), dimension(:), allocatable :: nodesinextsub
  end type ddinfo_str
  


  type(cluster), dimension(:), allocatable :: ddcoarse_mesh

  type(node), public     :: coarse_nodes
  type(element), public  :: coarse_elements
  type(subdomain_str), dimension(:), allocatable :: subdomain
  type(ddinfo_str), public :: ddinfo
  
  type(extsmtx), public :: prolong_mtx
  type(extsmtx), public :: coarse_mtx
  !!!begin debug vars!!!!
  type(extsmtx), public :: coarse2
  type(extsmtx), public :: prolongT
  type(extsmtx), public :: tmpmat
  !!!end debug vars!!!!
  real(kind=rkind), dimension(:), allocatable, public :: bcoarse
  real(kind=rkind), dimension(:), allocatable, public :: xcoarse
  
  !> inner Ax=b criterion
  real(kind=rkind), public :: inner_criterion

  
  logical, public :: use_coarselev
  
  integer, public :: file_decomp
  integer, public :: file_ddinfo, file_dt
  
  private :: returnval
  
  contains
     function returnval(subdom, i, subtime) result(value)
	use typy
	use globals
	
	class(subdomain_str), intent(in) :: subdom
	!> local index
	integer(kind=ikind), intent(in) :: i
	real(kind=rkind), intent(in) :: subtime
	real(kind=rkind) :: value
	
	if (abs(subtime-subdom%time)  < epsilon(subtime)) then
	  value = subdom%xvect(i,2)
	else
	  value = (subdom%xvect(i,2) - subdom%xvect(i,1))/subdom%time_step*(subtime-subdom%time)+subdom%xvect(i,2)
	end if
	

     end function returnval



end module decomp_vars