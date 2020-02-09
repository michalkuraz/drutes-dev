! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher


! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.



!> \file global_objs.f90
!! \brief main object definitions
!<

!> Main objects with global attribute are defined here.


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
    !> if true then the values for quadpnt are returned from the extended subdomain data, 
    logical :: extended = .false.
    !> if .false. then the solution value will be evaluated for different time than the calculated time 
    logical :: globtime=.true.
    real(kind=rkind) :: time4eval
    logical :: debugstop=.false.
    !> use preprocessed values (some PDE problems e.g. Richards equation in total hydraulic head form distiguish between pressure head h and total hydraulic head H, where H is the solution, but e.g. the function for the water content (retention curve) requires h. In this case for evaluating the water content we need to use preprocessed values. The preprocessor could be created as a part of model setup by changing the default pointer pde%getval to your own routine. Because pde%getval is called both from your constitutive functions and from the FEM solver, you should be able to tell your own getval function, which value you want to get H or h? By default it is false.
    logical :: preproc=.false.
    real(kind=rkind), dimension(2) :: xy
  end type integpnt_str
  
  !> smart array, you don't need to allocate -- usefull if you don't know how many data you will write in, 
  !! it is simply smart enough, it does all allocations automatically. By default the array data allocates on dimension 2, if the data 
  !! doesn't fit into this array, then the array is reallocated for a double of its original size
  !<
  type, public :: smartarray_int
    !> array with data
    integer(kind=ikind), dimension(:), allocatable :: data
    !> put there what ever you like, sometimes it can be usefull if you have this array
    integer(kind=ikind), dimension(:), allocatable :: info
    !> current dimension of your data, e.g. your smart array is declared as smart, then don't try ubound(smart%%data) NEVER!!
    !! for the array dimension see the value smart%pos
    !<
    integer(kind=ikind) :: pos=0, infopos=0
    contains
      !> fills value into array, its dimension is automaticaly allocated
      procedure :: fill => ismartfill
      !> clear the smart array, input parameter is logical FULL, this parameter is optional
      !! by default if FULL not passed or FULL is .false., then only the value smart%pos is set to zero
      !! if FULL .true., the also the array data is deallocated
      !! it turns out that FULL .false. is faster, but less memory efficient
      !<
      procedure :: clear => ismartclear
      !> disjoint fill -- it fills value into array, only if the same value does not already exist in the array
      procedure :: nrfill => ismartfill_norepeat
      !> return value is logical, if the value already exist in the array, then the reutrn value is .true. otherwise .false.
      procedure :: exist => ismartexist
  end type smartarray_int
  
  
  !> For comments and description see smartarray_int. 
  type, public ::  smartarray_real
    real(kind=rkind), dimension(:), allocatable :: data
    integer(kind=ikind), dimension(:), allocatable :: info
    integer(kind=ikind) :: pos=0, infopos=0
    contains
      !> fills value into array, its dimension is automaticaly allocated
      procedure :: fill => rsmartfill
      procedure :: clear => rsmartclear
      !> disjoint fill -- it fills value into array, only if the same value does not exist in the array. Since we deal with real numbers in this structure, then the real number is checked for absolute value of differece between two real numbers
      !! the value is checked as 
      !!  abs(array%data(i) - input) < max(abs(array%data(i)), abs(input))*(epsilon(input)
      !<
      procedure :: nrfill => rsmartfill_norepeat
  end type smartarray_real  
  
  !> polymorphic extension of sparse matrix
  type, public, extends(smtx) :: extsmtx
    real(kind=rkind), dimension(:), allocatable :: weight
    logical :: weighted
    type(smartarray_int) :: rowsfilled
  end type extsmtx
  
  type, public :: dirglob_str
    character(len=4096) :: dir
    logical :: valid=.false.
  end type dirglob_str
  

 !> DRUtES version structure
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
    integer(1) :: dimen
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
    !! 1 - to tabelarize the values and linearly approximate values between
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
    !> descriptor of the problem type (REstdH = standard Richards equation in total hydraulic head form 
    character(len=256) :: name
    !> full name of the partial differential equation problem -- e.g. dual permeability Richards equation in total hydraulic head form
    character(len=4096) :: fullname
    !> if .true. then DRUtES solves rotational symmetric flow
    logical :: rotsym=.false.
    !> if .true. then DRUtES checks for integral mass balance errors, requires more computational power
    logical :: check4mass=.false.
  end type configuration


    
  !> structure to define observation points
  type, public :: observation
    
    real(kind=rkind), dimension(:), allocatable :: xyz
    integer(kind=ikind) :: element
    real(kind=rkind), dimension(:), allocatable :: cumflux
    logical :: boundary = .false.
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
    !> format of output files for observation  real(kind=rkind), dimension(:,:), allocatable, private :: exp_data
    real(kind=rkind), dimension(:,:), allocatable, private :: model_data
    character(len=4) :: output_fmt
    logical :: isboundary = .true.
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
    !> logical array.if true, then the particular node is the domain boundary node, the domain should be obviously Lipshitz type
    logical, dimension(:), allocatable             :: boundary
    !> integer array, if the node is a boundary node, then this stores the order of boundary node on the boundary curve (for 2D)
    integer(kind=ikind), dimension(:), allocatable :: boundary_order
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
    !> vertical component of the inner boundary normal vector
    real(kind=rkind), dimension(:,:), allocatable   :: nvect_z
    !> material = id number of material at current element, a constant material properties are required for each element
    integer(kind=ikind), dimension(:), allocatable   :: material
    !> domain id -- array, type smartarray_int, carries id number of subdomain, the subdomain split is of an overlap type
    type(smartarray_int), dimension(:), allocatable :: subdom
    !> if element lies at domain boundary, then this array contains list of nodes at boundary
    type(smartarray_int), dimension(:), allocatable :: border
    !> list of elements at domain boundary
    type(smartarray_int) :: bcel
  end type element


  
  type, public :: integnodes
    real(kind=rkind), dimension(:,:), allocatable :: point
    real(kind=rkind), dimension(:), allocatable   :: weight
    real(kind=rkind)                              :: area
  end type integnodes

  
  private :: ismartfill, ismartclear, ismartfill_norepeat, rsmartfill, rsmartclear, rsmartfill_norepeat, ismartexist
  
  contains
    !> writes into smartarray integer vectors
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
            print *, "called from global_objs::rsmartfill"
            print *, "contact Michal -> michalkuraz@gmail.com"
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
   
    
    
    !> writes into smartarray integer vectors, only if the value is not present in the smartvector
    subroutine ismartfill_norepeat(array, input, info)
      use typy
      class(smartarray_int), intent(in out) :: array
      integer(kind=ikind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: i
!       logical :: exist
!       
!       exist = .false.
!       if (allocated(array%data)) then
!         do i=1, array%pos
!           if (array%data(i) == input) then
!             exist = .true.
!             EXIT
!           end if
!         end do
!        end if

      
      if (.not. ismartexist(array,input)) then
        if (present(info)) then
          call ismartfill(array, input, info)
        else
          call ismartfill(array, input)
        end if
      end if
    
    
    end subroutine ismartfill_norepeat
    
    !> clears (deallocates) smart vector
    subroutine ismartclear(array, full)
      class(smartarray_int), intent(in out) :: array
      !> if present and .true. the smart vector is completely deallocated
      logical, intent(in), optional :: full
      
      if (present(full)) then
         if (full .and. allocated(array%data)) then
           deallocate(array%data)
         end if
      end if

      array%pos = 0
      
    end subroutine ismartclear
    
    !> checks if the value is present in the vector
    function ismartexist(array, value) result(exist)
      use typy
      class(smartarray_int), intent(in) :: array
      integer(kind=ikind), intent(in) :: value
      logical :: exist
      
      integer(kind=ikind) :: i
      
      if (minval (abs(array%data(1:array%pos) - value)) == 0) then
        exist = .true.
      else
        exist = .false.
      end if
      
      
    end function ismartexist
    
    !> write real into smartarray vector
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
            print *, "called from global_objs::rsmartfill"
            print *, "contact Michal -> michalkuraz@gmail.com"
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
   
    
    
    !> writes real into smartarray vector, only if value sufficiently differ from the values already present. Analogical to smartfill_norepeat for integers
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
          if (abs(array%data(i) - input) < max(abs(array%data(i)), abs(input))*(epsilon(input))) then
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
    
    !> clears / deallocates real smartarray vector
    subroutine rsmartclear(array, full)
      class(smartarray_real), intent(in out) :: array
      logical, intent(in), optional :: full
      
      if (present(full) ) then
        if (full .and. allocated(array%data)) then
          deallocate(array%data)
        end if
      end if

      array%pos = 0
      
    end subroutine rsmartclear



 
      


end module global_objs
