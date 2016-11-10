module global_objs
  use typy
  use sparsematrix
  implicit none
  
  !> structure that specifies the point 
  type, public :: integpnt_str
    !> type_pnt value specifies the following
    !! type_pnt = gqnd -- Gauss quadrature node
    !! type_pnt = obpt -- observation point
    !! type_pnt = ndpt -- node from the mesh of geometrical discretization
    !<
    character(len=4) :: type_pnt
    !> this value must be supplied if type_pnt = obpt, or type_pnt = ndpt
    integer(kind=ikind) :: order
    !> these values (element, gqpt) has to be supplied if type_pnt = gqnd
    integer(kind=ikind) :: element
    !> this value specifies the vector of solution, that is used
    !! if column = 1 then previous time step is used
    !! if column = 2 then current iteration level is used
    !<
    integer(kind=ikind) :: column
    !> if true, then the values for quadpnt are returned from the local subdomain data, otherwise (default) the values are returned from the pde_common%xvect array
    logical :: ddlocal = .false.
    integer(kind=ikind) :: subdom
    !> if true then the values for quadpnt are returned from the exnteded subdomain data, 
    logical :: extended = .false.
    !> if .false. then the solution value will be evaluated for different time than the calculated time 
    logical :: globtime=.true.
    real(kind=rkind) :: time4eval
    logical :: debugstop=.false.
    !> use preprocessed values (some PDE problems e.g. Richards equation in total hydraulic head form distiguish between pressure head h and total hydraulic head H, where H is the solution, but e.g. the function for the water content (retention curve) requires h. In this case for evaluating the water content we need to use preprocessed values. The preprocessor could be created as a part of model setup by changing the default pointer pde%getval to your own routine. Because pde%getval is called both from your constitutive functions and from the FEM solver, you should be able to tell your own getval function, which value you want to get H or h? By default it is false.
    logical :: preproc=.false.
  end type integpnt_str
  
  type, public :: smartarray_int
    integer(kind=ikind), dimension(:), allocatable :: data
    integer(kind=ikind), dimension(:), allocatable :: info
    integer(kind=ikind) :: pos=0, infopos=0
    contains
      !> fills value into array, its dimension is automaticaly allocated
      procedure :: fill => ismartfill
      procedure :: clear => ismartclear
      !> disjoint fill -- it fills value into array, only if the same value does not exist in the array
      procedure :: nrfill => ismartfill_norepeat
      procedure :: exist => ismartexist
  end type smartarray_int
  
  type, public ::  smartarray_real
    real(kind=rkind), dimension(:), allocatable :: data
    integer(kind=ikind), dimension(:), allocatable :: info
    integer(kind=ikind) :: pos=0, infopos=0
    contains
      !> fills value into array, its dimension is automaticaly allocated
      procedure :: fill => rsmartfill
      procedure :: clear => rsmartclear
      !> disjoint fill -- it fills value into array, only if the same value does not exist in the array
      procedure :: nrfill => rsmartfill_norepeat
  end type smartarray_real  
  
  type, public, extends(smtx) :: extsmtx
    real(kind=rkind), dimension(:), allocatable :: weight
    logical :: weighted
    type(smartarray_int) :: rowsfilled
  end type extsmtx
  
  type, public :: dirglob_str
    character(len=4096) :: dir
    logical :: valid=.false.
  end type dirglob_str
  


  type, public :: version    
    !> number of version
    character(len=9) :: number
    !> beta release (b), stable release (s)
    character(len=7) :: reliability
  end type version
  

  !>structure with basic model configuration
  type, public :: configuration
    !> run dual yes or no
    character(len=1)    :: run_dual
    !> use damped newton for the modified pickard iteration y/n
    character(len=1)    :: damped_newton
    !> problem dimension 1 = 1D / 2 = 2D / 3 = 3D
    integer(kind=ikind) :: dimen
    !> mesh_type = 1 internal mesh generator (simple)
    !! mesh_type = 2 t3d mesh generator
    !< mesh type = 3 gmsh mesh generator
    integer(kind=ikind) :: mesh_type
    !> adapt time step to the observation time, or calculate the values of the observation time by linear approximation
    logical    :: adapt_observe
    !> parameter to decide if the code execution begins in some old backup
    logical    :: run_from_backup
    !> evaluate constitutive functions from table created at program init or directly?
    !! 0 - for direct evaluation
    !! 1 - to tabelarise the values and linearly approximate values between
    !<
    integer(kind=ikind) :: fnc_method
    !> length of interval in between the values in function table
    real(kind=rkind)    :: fnc_discr_length
    !>nonlinear iteration method
    !!0 - standard Picard method
    !! 1 - Schwarz-Picard method (dd-adaptivity)
    !!2 - no iterations
    !<
    integer(kind=ikind) :: it_method
    !> descriptor of the problem type (RE_std = standard Richards equation, RE_mod = modified Richards eq. (Noborio)
    character(len=256) :: name
    character(len=4096) :: fullname
  end type configuration


    
  !> structure to define observation points
  type, public :: observation
    real(kind=rkind), dimension(:), allocatable :: xyz
    integer(kind=ikind) :: element
    real(kind=rkind), dimension(:), allocatable :: cumflux
  end type observation
  
  type, public :: observe_time_str
    real(kind=rkind) :: value
    character(len=3) :: name
  end type observe_time_str
  
  type, public :: observe_info_str
    !> methods for observation time print 
    !! 1 - adjust time stepping to observation time values
    !! 2 - linearly interpolate solution between two consecutive solutions (recommended)
    !<
    integer(kind=ikind) :: method
    !> format of outputs for observation times
    !! pure - raw data are printed, just nodes id with FEM coefficients
    !! scil - scilab output 
    !! gmsh - gmsh output
    !<
    character(len=4) :: fmt
    logical :: anime
    integer(kind=ikind) :: nframes
    !> format of output files for observation
    character(len=4) :: output_fmt
  end type observe_info_str

  type, public, extends(observation) :: measured
    integer(kind=ikind) :: node
  end type measured
    



  !> mesh array type for node
  !! kolik = number of nodes
  !! id = id number of node, equal to position
  !! data = x,y coordinates
  !! bc = boundary condition for current node
  !! edge = boundary id number, if the node lies apart from any boundary, default value is 0
  !! results = final iteration at each time step will be copied into this vector
  !<
  type, public :: node
    integer(kind=ikind) ::  kolik
    integer(kind=ikind), dimension(:), allocatable :: id
    real(kind=rkind), dimension(:,:), allocatable  :: data
    integer(kind=ikind), dimension(:), allocatable :: edge
    real(kind=rkind), dimension(:,:), allocatable  :: results
    !> logical array.if true, then the particular node is the domain boundary node, the domain should be sufficiently Lipshitz type
    logical, dimension(:), allocatable             :: boundary
    !> array of elements el2integ(i)%data that are covered by basis function the basis function that originates from the node el2integ(i)
    type(smartarray_int), dimension(:), allocatable :: el2integ
    !> list of elements of geometrical discretization, where this node belongs
    type(smartarray_int), dimension(:), allocatable :: element
  end type node


  !> mesh array type for elements
  type, public :: element
    !> kolik = number of elements
    integer(kind=ikind) ::  kolik
    integer(kind=ikind), dimension(:), allocatable   :: id
    integer(kind=ikind), dimension(:,:), allocatable :: data
    real(kind=rkind), dimension(:), allocatable      :: areas
    !>
    !! the first value is the element identificator \n
    !! \n
    !! -------------------------------------------------- \n
    !! the second value is the basis function number \n
    !! 1D - (:,1,:) function [0,1] -> [1,0] 
    !!      (:,2,:) function [0,0] -> [1,1] \n
    !! 3D - (to be filled when implemented)
    !! 2D - (:,1,:) function [0,0,1] -> [1,0,0] -> [0,1,0] 
    !!      (:,2,:) function [0,0,0] -> [1,0,1] -> [0,1,0] 
    !!      (:,3,:) function [0,0,0] -> [1,0,0] -> [0,1,1] \n
    !! --------------------------------------------------- \n
    !! the third value defines derivatives with respect to particular axes x,y,z, and thus this value is equal to problem dimension \n
    !! (:,:,1) derivative with respect to x 
    !! (:,:,2) derivative with respect to y 
    !! (:,:,3) derivative with respect to z \n
    !! \n
    !<
    real(kind=rkind), dimension(:,:,:), allocatable  :: ders
    !> an array that stores data of neighbourhood elements (elements that shares an edge (two nodes))
    integer(kind=ikind), dimension(:,:), allocatable :: neighbours
    !> an array that contains coordinates of gravity centers of each element, the row number is the element number, and the column number is the coordinate -- x,y,z
    real(kind=rkind), dimension(:,:), allocatable    :: gc
    !> an array that contains lenght of the element edge, if the element edge is boundary edge, if the element 
    !! edge is apart the boundary, zero value is supplied
    !!
    !!         3
    !!       /|  
    !!      / |                   
    !!  l3 /  |  l2
    !!    /   | 
    !!   /____|
    !!  1      2
    !!     l1
    !!
    !! this figure explains the boundary length order in the lenght array, 
    !! thus the line between the nodes 1-2 has order length(1), the line between the nodes 3-1 (or 1-3) has order lenght(3)
    !< 
    real(kind=rkind), dimension(:,:), allocatable   :: length
    !> vertical component of the boundary normal vector
    real(kind=rkind), dimension(:,:), allocatable   :: nvect_z
    !> material = id number of material at current element, a constant material properties are required for each element
    integer(kind=ikind), dimension(:,:), allocatable   :: material
    !> domain id -- array, type smartarray_int, carries id number of subdomain, the subdomain split is of an overlap type
    type(smartarray_int), dimension(:), allocatable :: subdom
  end type element


  type, public :: integnodes
    real(kind=rkind), dimension(:,:), allocatable :: point
    real(kind=rkind), dimension(:), allocatable   :: weight
    real(kind=rkind)                              :: area
  end type integnodes

  
  private :: ismartfill, ismartclear, ismartfill_norepeat, rsmartfill, rsmartclear, rsmartfill_norepeat, ismartexist
  
  contains
    subroutine ismartfill(array,input, info)
      use typy
      class(smartarray_int), intent(in out) :: array
      integer(kind=ikind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: l
      integer(kind=ikind), dimension(:), allocatable :: itmp
      integer(kind=ikind), dimension(:), allocatable  :: logtmp
      
      ! check allocations
      if (.not. allocated(array%data)) then
        array%pos = 0
        allocate(array%data(1))
        if (present(info)) then
	  array%infopos = 0
	  allocate(array%info(1))
	end if
      end if
      
      array%pos = array%pos+1
      
      if (present(info)) then
	array%infopos = array%infopos + 1
        if (array%pos /= array%infopos) then
	    print *, "this is a bug in the code"
	    print *, "smartarray%data and smartarray%info has different amount of data written"
	    ERROR STOP
	end if
      end if
      
      if (ubound(array%data,1) < array%pos) then
        l = ubound(array%data,1)
        allocate(itmp(l))
        itmp = array%data
        deallocate(array%data)
        allocate(array%data(2*l))
        array%data(1:l) = itmp
        deallocate(itmp)
        if (present(info)) then
	  allocate (logtmp(l))
	  logtmp = array%info
	  deallocate(array%info)
	  allocate(array%info(2*l))
	  array%info(1:l) = logtmp
	  deallocate(logtmp)
	end if	  
      end if
      
      array%data(array%pos) = input
      
      if (present(info)) then
	array%info = info
      end if
      
      
    end subroutine ismartfill
   
    
    
    
    subroutine ismartfill_norepeat(array, input, info)
      use typy
      class(smartarray_int), intent(in out) :: array
      integer(kind=ikind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: i
      logical :: exist
      
      exist = .false.
      if (allocated(array%data)) then
	do i=1, array%pos
	  if (array%data(i) == input) then
	    exist = .true.
	    EXIT
	  end if
	end do
       end if
      
      if (.not. exist) then
	if (present(info)) then
	  call ismartfill(array, input, info)
	else
	  call ismartfill(array, input)
	end if
      end if
    
    
    end subroutine ismartfill_norepeat
    
    subroutine ismartclear(array, full)
      class(smartarray_int), intent(in out) :: array
      logical, intent(in), optional :: full
      
      if (present(full) .and. full .and. allocated(array%data)) then
        deallocate(array%data)
      end if

      array%pos = 0
      
    end subroutine ismartclear
    
    function ismartexist(array, value) result(exist)
      use typy
      class(smartarray_int), intent(in) :: array
      integer(kind=ikind), intent(in) :: value
      logical :: exist
      
      integer(kind=ikind) :: i
      
      exist = .false.
      
      do i=1, array%pos
	if (array%data(i) == value) then
	  exist = .true.
	  RETURN
	end if
      end do
      
    end function ismartexist
    

    subroutine rsmartfill(array,input, info)
      use typy
      class(smartarray_real), intent(in out) :: array
      real(kind=rkind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: l
      real(kind=rkind), dimension(:), allocatable :: rtmp
      integer(kind=ikind), dimension(:), allocatable  :: logtmp
      
      ! check allocations
      if (.not. allocated(array%data)) then
        array%pos = 0
        allocate(array%data(1))
        if (present(info)) then
	  array%infopos = 0
	  allocate(array%info(1))
	end if
      end if
      
      array%pos = array%pos+1
      
      if (present(info)) then
	array%infopos = array%infopos + 1
        if (array%pos /= array%infopos) then
	    print *, "this is a bug in the code"
	    print *, "smartarray%data and smartarray%info has different amount of data written"
	    ERROR STOP
	end if
      end if
      
      if (ubound(array%data,1) < array%pos) then
        l = ubound(array%data,1)
        allocate(rtmp(l))
        rtmp = array%data
        deallocate(array%data)
        allocate(array%data(2*l))
        array%data(1:l) = rtmp
        deallocate(rtmp)
        if (present(info)) then
	  allocate (logtmp(l))
	  logtmp = array%info
	  deallocate(array%info)
	  allocate(array%info(2*l))
	  array%info(1:l) = logtmp
	  deallocate(logtmp)
	end if	  
      end if
      
      array%data(array%pos) = input
      
      if (present(info)) then
	array%info = info
      end if
      
      
    end subroutine rsmartfill
   
    
    
    
    subroutine rsmartfill_norepeat(array, input, info)
      use typy
      class(smartarray_real), intent(in out) :: array
      real(kind=rkind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: i
      logical :: exist
      
      exist = .false.
      if (allocated(array%data)) then
	do i=1, array%pos
	  if (abs(array%data(i) - input) < epsilon(input)) then
	    exist = .true.
	    EXIT
	  end if
	end do
       end if
      
      if (.not. exist) then
	if (present(info)) then
	  call rsmartfill(array, input, info)
	else
	  call rsmartfill(array, input)
	end if
      end if
    
    
    end subroutine rsmartfill_norepeat
    
    subroutine rsmartclear(array, full)
      class(smartarray_real), intent(in out) :: array
      logical, intent(in), optional :: full
      
      if (present(full) .and. full .and. allocated(array%data)) then
        deallocate(array%data)
      end if

      array%pos = 0
      
    end subroutine rsmartclear



 
      


end module global_objs
