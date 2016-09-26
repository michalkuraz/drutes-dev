module pde_objs
  use typy
  use sparsematrix 
  use global_objs
  use globals
  use decomp_vars
  implicit none
  
  
  type, public :: pde_fnc
    procedure(tensor_fnc), nopass, pointer           :: dispersion
    procedure(vector_fnc), nopass, pointer           :: convection
    !> derivative of the convection function, because \f[ \nabla .(a u) = (\nabla . a) u + a (\nabla . u) \f]
    procedure(vector_fnc), nopass, pointer           :: der_convect
    !> reaction
    procedure(scalar_fnc), nopass, pointer           :: reaction
    !> zero order reaction
    procedure(scalar_fnc), nopass, pointer           :: zerord
    procedure(scalar_fnc), nopass, pointer           :: elasticity
    logical                                          :: coupling
  end type pde_fnc
  

  type, public :: pde_common_str
    real(kind=rkind), dimension(:), allocatable :: bvect
    real(kind=rkind), dimension(:,:), allocatable :: xvect
    integer(kind=ikind), dimension(:), allocatable :: invpermut
    logical :: nonlinear
    procedure(timeint_fnc), nopass, pointer          :: time_integ
    integer :: timeint_method
    integer(kind=ikind)                              :: processes
    integer(kind=ikind)                              :: current_proc
    logical, dimension(:,:), allocatable             :: coupling
    !> id number of the current element (triangle), that is currently integrated
    integer(kind=ikind) :: current_el
    !> nonlinear solver method (standard Picard, Schwarz Picard, etc...)
    procedure(treat_pde_subrt), nopass, pointer       :: treat_pde
  end type pde_common_str
  
  !> data structure  to carry boundary value problem setup, allocation starts at 100
  type, public :: boundary_vals
    integer(kind=ikind) :: ID
    integer(kind=ikind) :: code
    logical             :: file
    real(kind=rkind)    :: value
    !>  value of the boundary condition, typically this is a function 
    procedure(bc_fnc), nopass, pointer    :: value_fnc
    !> in case of the Neumann boundary condition the value that is supplied into the system of equations contains only diffusion flux, and thus the convective flux should be subtracted 
    real(kind=rkind)    :: convect_value
    !> if file == .true., then series contains unsteady data
    !! series(:,1) = time
    !! series(:,2) = values
    !<
    real(kind=rkind), dimension(:,:), allocatable :: series
    !> current position in series data
    integer(kind=ikind) :: series_pos
    integer(kind=ikind) :: xy_count
    integer(kind=ikind) :: layer
    !> if icond4neumann == .true. then the initial condition is supplied into the Neumann boundary nodes, otherwise the initial condition in the Neumann nodes is evaluated for such value, that the Neumann condition is satisfied
    logical :: icond4neumann = .false.
  end type boundary_vals
  


  !> type definition for common quasilinear partial diferential equation in a format
  !! \f[ E(u)\frac{\partial u}{\partial t} = \nabla . D(u) \nabla u - \nabla . (c(u) u) + r(u) u \f]
  !<
  type, public :: PDE_str
    !> the first item is used for filename 
    !! the second item is printed inside the text files
    !>
    character(len=256), dimension(2)                 :: problem_name
    character(len=64), dimension(2)                  :: solution_name
    character(len=64), dimension(2)                  :: flux_name
    character(len=64), dimension(2)                  :: mass_name
    type(pde_fnc), dimension(:), allocatable         :: pde_fnc
    procedure(scalar_fnc), pass(pde_loc), pointer    :: mass
    procedure(vector_fnc), pass(pde_loc), pointer    :: flux
    procedure(time_check), pass(pde_loc), pointer    :: dt_check
    procedure(icond_fnc), pass(pde_loc), pointer     :: initcond
    !> procedure to be called after process change
    procedure(basic_subrt), nopass, pointer          :: process_change
    procedure(basic_subrt), nopass, pointer          :: read_parameters
    !> bc is allocated in read_inputs::readbcval
    type(boundary_vals), dimension(:), allocatable   :: bc
    integer(kind=ikind), dimension(:), allocatable   :: permut
    real(kind=rkind), dimension(:), allocatable      :: solution
    !> contains units of files opened for particular observation points
    integer, dimension(:), allocatable               :: obspt_unit
    !> procnodes(1) = lower id of the process base function
    !! procnodes(2) = upper id of the process base function
    integer, dimension(2) :: procbase_fnc
    integer(kind=ikind) :: order
    procedure(getval_str), pass(pde_loc), pointer :: getval
    contains 
      !> get vector of gradient of the solution
      procedure :: getgrad=>getgradp1
  end type PDE_str
  
  
  abstract interface
    function getval_str(pde_loc, quadpnt) result(value)
      use typy
      use global_objs
      import::pde_str
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: value
    end function getval_str
  end interface

  abstract interface
      subroutine matrix_solver(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
                  ll1,ll2,cond1,opcnt1,errcode1)
        use mtx
        use typy
        implicit none
        !> matice soustavy\n
        !! musi poskytovat getn, getm, mul (nasobeni vektorem)
        class(matrix), intent(in out) :: A
        !> vektor prave strany
        real(kind=rkind), dimension(:), intent(in) :: b
        !> aproximace reseni, postupne menena
        real(kind=rkind), dimension(:), intent(in out) :: x
        !> maximalni povoleny pocet iteraci, default = n ( Rozmer matice)
        integer(kind=ikind), intent(in), optional :: itmax1
        !> pozadovana relativni zmena rezidua, default = 1e-6
        real(kind=rkind), intent(in), optional :: reps1
        !> informacni podrobnost\n
        !> - 0 ... pracuj tise
        !! - 1 ... minimalni informace
        !! - 10 ... maximalni ukecanost
        integer, intent(in), optional :: ilev1
        !> skutecne provedeny pocet iteraci
        integer(kind=ikind), intent(out), optional :: itfin1
        !> skutecne dosazena relativni zmena residua
        real(kind=rkind), intent(out), optional :: repsfin1
        !> odhad nejvetsiho vlastniho cisla
        real(kind=rkind), intent(out), optional :: ll1
        !> odhad nejmensiho vlastniho cisla
        real(kind=rkind), intent(out), optional :: ll2
        !> odhad cisla podminenosti : cond1 = ll1/ll2
        real(kind=rkind), intent(out), optional :: cond1
        !> celkovy pocet provedenych operaci a cas behu
        type(tcount), intent(out), optional :: opcnt1
        !> kod pripadnr chyby
        !! - 0 ... OK
        !! - 1 ... matice neni ctvercova
        !! - 2 ... nesouhlasi b
        !! - 3 ... nesouhasi x
        !! - 4 ... ani jeden z vektoru nesouhlasi
        !! - 5 ... vycerpan povoleny pocet iteraci
        !! - 6 ... prestalo klesat residuum i energie
        integer, intent(out), optional :: errcode1
    end subroutine matrix_solver
  end interface 



 !> abstract interface for boundary value function
  abstract interface
    subroutine bc_fnc(pde_loc, element, node, value, code) 
      use typy
      import :: pde_str
      class(pde_str), intent(in) :: pde_loc
      !> element id
      integer(kind=ikind), intent(in)  :: element
      !> node id
      integer(kind=ikind), intent(in)  :: node
      !> return value
      real(kind=rkind), intent(out), optional    :: value
      !> return type of boundary condition
      integer(kind=ikind), intent(out), optional :: code
    end subroutine bc_fnc
  end interface 

  !> abstract interface for scalar function
  abstract interface
    function scalar_fnc(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      import :: pde_str
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
    end function scalar_fnc
  end interface 

  !> abstract interface for vector function
  abstract interface
    subroutine vector_fnc(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      import :: pde_str
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> input vector
      real(kind=rkind), dimension(:), intent(in), optional  :: vector_in
      !> output vector
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional               :: scalar
    end subroutine vector_fnc
  end interface 



  !> abstract interface for vector function
  abstract interface
    subroutine tensor_fnc(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      import :: pde_str
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
    end subroutine tensor_fnc 
  end interface


  abstract interface
    subroutine timeint_fnc(el_id, domain_id, quadpnt_in)
      use typy
      use global_objs
      integer(kind=ikind), intent(in) :: el_id
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str), intent(in out), optional :: quadpnt_in
    end subroutine timeint_fnc
  end interface

  abstract interface
    function logical_fnc() result(valid)
      logical :: valid
    end function logical_fnc
  end interface
  
  abstract interface
    function time_check(pde_loc) result(ok)
      use typy
      import :: pde_str
      class(pde_str), intent(in) :: pde_loc
      logical :: ok
    end function time_check
  end interface


  abstract interface
    subroutine icond_fnc(pde_loc)
      use typy
      import :: pde_str
      class(PDE_str), intent(in out) :: pde_loc
    end subroutine icond_fnc
  end interface


  abstract interface
    subroutine basic_subrt()
    end subroutine basic_subrt
  end interface


  abstract interface
    subroutine treat_pde_subrt(ierr, itcount, success)
      use typy
      integer, intent(out) :: ierr
      integer(kind=ikind), intent(out) :: itcount
      logical, intent(out) :: success
    end subroutine treat_pde_subrt
  end interface
  
  public :: getgradp1, getvalp1, do_nothing
  private :: getvalp1loc
  
  type(PDE_str), dimension(:), allocatable,  public :: PDE
  type(pde_common_str), public :: pde_common
  procedure(matrix_solver), pointer, public :: solve_matrix

  
  contains 
  

  
    subroutine getgradp1(pde_loc, quadpnt, grad)
      use typy
      use decomp_vars

      
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      type(integpnt_str), dimension(:), allocatable, save :: quadpntloc
      real(kind=rkind), dimension(:), allocatable, intent(out) :: grad
      real(kind=rkind), dimension(3) :: gradloc
      integer(kind=ikind), dimension(:), allocatable, save :: pts
      real(kind=rkind), dimension(3)    :: a,b,c
      real(kind=rkind) :: dx
      integer(kind=ikind) :: i, el, top, j

      
      if (.not. allocated(grad)) then
	allocate(grad(drutes_config%dimen))
      else if (ubound(grad,1) /= drutes_config%dimen ) then
	deallocate(grad)
	allocate(grad(drutes_config%dimen))
      end if
      
      if (.not. allocated(pts)) then
	allocate(pts(ubound(elements%data,2)))
	allocate(quadpntloc(ubound(elements%data,2)))
      end if
      
      select case(quadpnt%type_pnt)
	case("gqnd", "obpt")
	  top = 1
	case("ndpt")
	  top = nodes%element(quadpnt%order)%pos
	case default
	  print *, "RUNTIME ERROR: incorrect quadpnt type definition (value quadpnt%type_pnt)"
	  print *, "the value specified in code was:", quadpnt%type_pnt
	  print *, "exited from pde_objs::getgradp1"
	  ERROR STOP
	  
	  
      end select
      
      gradloc = 0
      do i=1, top  
	select case(quadpnt%type_pnt)
	  case("gqnd")
	    el = quadpnt%element
	  case("obpt")
	    el = observation_array(quadpnt%order)%element
	  case("ndpt")
	    el = nodes%element(quadpnt%order)%data(i)
	end select
	
	pts = elements%data(el,:)
	
	quadpntloc(:) = quadpnt
	quadpntloc(:)%type_pnt = "ndpt"
	select case(drutes_config%dimen)
	  case(1)
	    dx = nodes%data(pts(2),1) - nodes%data(pts(1),1)
	    quadpntloc(1)%order = pts(1)
	    quadpntloc(2)%order = pts(2)
	     gradloc(1) = gradloc(1) + (getvalp1(pde_loc,quadpntloc(2)) - getvalp1(pde_loc, quadpntloc(1)))/dx
	  case(2)
	    a(1:2) = nodes%data(pts(1),:)
	    b(1:2) = nodes%data(pts(2),:)
	    c(1:2) = nodes%data(pts(3),:)
	    
	    
	    quadpntloc(1)%order = pts(1)
	    quadpntloc(2)%order = pts(2)
	    quadpntloc(3)%order = pts(3)
	    
	    
	    
	    a(3) = getvalp1(pde_loc, quadpntloc(1))
	    b(3) = getvalp1(pde_loc, quadpntloc(2))
	    c(3) = getvalp1(pde_loc, quadpntloc(3))
	    call get2dderivative(a,b,c,grad(1), grad(2))

	    gradloc(1:2) = gradloc(1:2) + grad
	  case(3)
	end select
      end do
      
      grad = gradloc(1:drutes_config%dimen)/top
    
    end subroutine getgradp1
    
    function getvalp1(pde_loc, quadpnt) result(val)
      use typy
      use decomp_vars
      
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      real(kind=rkind) :: tmp
      
      type(integpnt_str) :: quadpntloc
      real(kind=rkind) :: valprev, timeprev, timeglob
      logical :: stopper=.false.

   
       val = getvalp1loc(pde_loc, quadpnt)


      if (quadpnt%globtime .or. quadpnt%column==1) then
	RETURN
      else
	quadpntloc = quadpnt
	quadpntloc%column = 4
	if (quadpnt%ddlocal) then
	  timeglob = subdomain(quadpnt%subdom)%time
	  timeprev = subdomain(quadpnt%subdom)%time - subdomain(quadpnt%subdom)%time_step
	else
	  timeglob = time
	  timeprev = time-time_step
	end if

	valprev = getvalp1loc(pde_loc, quadpntloc)
        tmp=val
	val = (val-valprev)/(timeglob-timeprev)*(quadpnt%time4eval-timeprev)+valprev

      end if
	
    
    end function getvalp1
    
    function getvalp1loc(pde_loc, quadpnt, stopme) result(val)
      use typy
      use decomp_vars

      
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      logical, intent(in), optional :: stopme
      real(kind=rkind) :: val
      
      integer(kind=ikind), dimension(:), allocatable, save :: pts, ppts
      real(kind=rkind), dimension(:), allocatable, save :: ndvals
      integer(kind=ikind) :: i, edge, el, j, order
      real(kind=rkind) :: xder, yder
      real(kind=rkind), dimension(3,3) :: a
      
            
      
     select case(quadpnt%type_pnt)
	case("gqnd", "obpt")
	    if (.not. allocated(pts) ) then
	      allocate(pts(ubound(elements%data,2)))
	      allocate(ppts(ubound(elements%data,2)))
	      allocate(ndvals(ubound(elements%data,2)))
	    end if
	    
	    
	    select case(quadpnt%type_pnt)
	      case("gqnd")
		el = quadpnt%element
	      case("obpt")
		el = observation_array(quadpnt%order)%element
	    end select
	      
	    pts = elements%data(el,:)
	    ppts = pde_loc%permut(pts)
	    
		
            if (quadpnt%ddlocal) then
	      if (.not. quadpnt%extended) then
		where (ppts /= 0)
		  ppts = subdomain(quadpnt%subdom)%invpermut(ppts)
		end where
	      else 
		where (ppts /= 0)
		  ppts = subdomain(quadpnt%subdom)%extinvpermut(ppts)
		end where
	      end if
	    end if


	    
	    do i=1, ubound(elements%data,2)
	
	      if (ppts(i) > 0) then

		if (.not. quadpnt%ddlocal) then
		  ndvals(i) = pde_common%xvect(ppts(i), quadpnt%column)
		else
		  if (.not. quadpnt%extended) then
		    ndvals(i) = subdomain(quadpnt%subdom)%xvect(ppts(i), quadpnt%column)
		  else
		    ndvals(i) = subdomain(quadpnt%subdom)%extxvect(ppts(i), quadpnt%column)
		  end if
		end if
	      else

		edge = nodes%edge(pts(i))

                
		call pde_loc%bc(edge)%value_fnc(pde_loc, el, i, ndvals(i))
		
	      end if
	    end do

       
	    select case(quadpnt%type_pnt)
	      case("gqnd")
		select case(drutes_config%dimen)
		  case(1)
		    val = (ndvals(2) - ndvals(1))*gauss_points%point(quadpnt%order,1) + ndvals(1)
		  case(2)
		    call get2dderivative((/0.0_rkind, 0.0_rkind, ndvals(1)/), (/1.0_rkind, 0.0_rkind, ndvals(2)/), (/0.0_rkind, &
				   1.0_rkind, ndvals(3)/), xder, yder)
		    val = ndvals(1) + xder*gauss_points%point(quadpnt%order,1) + yder*gauss_points%point(quadpnt%order,2)
		end select
	      case("obpt")
		select case(drutes_config%dimen)
		  case(1)
		     val = (ndvals(2) - ndvals(1))/(nodes%data(pts(2),1) - nodes%data(pts(1),1))*(observation_array(quadpnt%order)%xyz(1) - &
			   nodes%data(pts(1),1)) + ndvals(1)
		  case(2)
		    do i=1,3
		      do j=1,2
			a(i,j) = nodes%data(pts(i),j)
		      end do
		      a(i,3) = ndvals(i)
		    end do
		    	    
		    call get2dderivative(a(1,:), a(2,:), a(3,:), xder, yder)
		    
		    val = ndvals(1) + xder*(observation_array(quadpnt%order)%xyz(1) - a(1,1)) + &
			  yder * (observation_array(quadpnt%order)%xyz(2) - a(1,2))
		end select
	  end select
		    
	case("ndpt")
	    i = pde_loc%permut(quadpnt%order)

	    if (i > 0) then
	      if (.not. quadpnt%ddlocal) then
		val = pde_common%xvect(i, quadpnt%column)
	      else
	      
		if (.not. quadpnt%extended) then

		  i = subdomain(quadpnt%subdom)%invpermut(i)

		  val = subdomain(quadpnt%subdom)%xvect(i, quadpnt%column)

		else
		  i = subdomain(quadpnt%subdom)%extinvpermut(i)
		  val = subdomain(quadpnt%subdom)%extxvect(i, quadpnt%column)
		end if
	      end if
	    else
	      edge = nodes%edge(quadpnt%order)
	      el = nodes%element(quadpnt%order)%data(1)
	      do i=1, ubound(elements%data,2)
		if (elements%data(el,i) == quadpnt%order) then
		  order = i
		  EXIT
		end if
	      end do
	      call pde_loc%bc(edge)%value_fnc(pde_loc, el, order, val)
	    end if

	case default
	    print *, "RUNTIME ERROR: incorrect quadpnt type definition (value quadpnt%type_pnt)"
	    print *, "the value specified in code was:", quadpnt%type_pnt
	    print *, "exited from pde_objs::getvalp1"
	    ERROR STOP
      end select
	  

    end function getvalp1loc
    
      
    subroutine get2dderivative(a,b,c,xder,yder)
      use typy
      !> 1st point of the plane
      real(kind=rkind), dimension(:), intent(in) :: a
      !> 2nd point of the plane
      real(kind=rkind), dimension(:), intent(in) :: b
      !> 3rd point of the plane
      real(kind=rkind), dimension(:), intent(in) :: c
      !> resulting x derivate
      real(kind=rkind), intent(out) :: xder
      !> resulting y derivate
      real(kind=rkind), intent(out) :: yder
      !-------local variables--------------
      real(kind=rkind), dimension(3) :: u
      real(kind=rkind), dimension(3) :: v
      real(kind=rkind), dimension(3) :: n
      integer(kind=ikind) :: i
      real(kind=rkind) :: reps


      reps = epsilon(reps)


      !check if the plane is not horizontal
      if (abs(a(3) - b(3)) < reps*(abs(a(3))-abs(b(3)))  .and.  &
	  abs(a(3) - c(3)) < reps*(abs(a(3))-abs(c(3)))  .and.  &
	  abs(b(3) - c(3)) < reps*(abs(b(3))-abs(c(3)))) then
	xder = 0.0_rkind
	yder = 0.0_rkind
	RETURN
      else
	CONTINUE
      end if

      !creates the plane vectors 
      do i=1,3
	u(i) = a(i) - b(i)
	v(i) = a(i) - c(i)
      end do

      ! the normal plane vector is defined as
      n(1) = u(2)*v(3) - v(2)*u(3)
      n(2) = u(3)*v(1) - v(3)*u(1)
      n(3) = u(1)*v(2) - v(1)*u(2)

      ! finally the derivate is as follows, the horizontality has been already checked
      ! the verticality check
      if (abs(n(3)) < 1e2*reps) then
	print *, "the mesh is wrong, base function can't be vertical"
	ERROR STOP
      end if 
      xder = -n(1)/n(3)
      yder = -n(2)/n(3)

      
    end subroutine get2dderivative 
    
    subroutine do_nothing()
    
    end subroutine do_nothing
    
end module pde_objs