module RE_constitutive
  use typy
  use global_objs
  use re_globals

  public :: vangen, vangen_tab, inverse_vangen
  public :: vangen_elast, vangen_elast_tab
  public :: mualem, mualem_tab
  public :: dmualem_dh, dmualem_dh_tab
  public :: convection_rerot
  public :: darcy_law
  public :: domainswitch
  public :: init_zones
  public :: rcza_check
  public :: derKz
  public :: re_dirichlet_bc, re_neumann_bc, re_null_bc, re_dirichlet_height_bc, re_initcond
  public :: sinkterm
  private :: zone_error_val, secondder_retc,  rcza_check_old, setmatflux, intoverflow

 
  real(kind=rkind), dimension(:,:), pointer, public :: Ktab, dKdhtab, warecatab, watcontab
  real(kind=rkind), dimension(:,:,:), allocatable, target, public :: Ktab_all, dKdhtab_all, warecatab_all, watcontab_all
  type(soilpar), private :: van_gen_coeff
  integer(kind=ikind), private :: this_layer
  
  logical, private :: tabwarning = .false.
  character(len=4096), private :: tabmsg = "DRUtES requested values of the constitutional functions out of the range of &
					    your precalculated table, the function will be evaluated directly, it will & 
					    slow down your computation, is everything ok? If yes, increase the &
					    length of interval for precaculating the constitutive functions in your config files"   
					    
  logical, private :: intwarning = .false.


  contains
  
  
    subroutine intoverflow()
      use core_tools
      use postpro
      character(len=4096) :: msg
      
      write(msg, fmt=*) "You are probably experiencing oscilations", & 
	      "(probably convective dominance), I will print your data. Your computation (just for a case) is not", &
	      "interrupted, but the results are probably not reliable anymore.",  & 
	      "We will probably experience a crash soon.:) Check your model configuration to avoid", &
	      "convective dominant problems."
	      
      call write_log(msg)
      
      call make_print("separately")
	      
      
    end subroutine intoverflow

    !> \brief Van Genuchten relation \f[ \theta = f(pressure) \f]
    !!  \f[ \theta_e = \frac{1}{(1+(\alpha*h)^n)^m} \f]
    !! water content is considered as absolute value not the relative one \n
    !! see \f[ \theta_e = \frac{\theta - \theta_r}{\theta_s-\theta_r} \f]
    !<
    function vangen(pde_loc, layer, quadpnt, x) result(theta)
      use typy
      use re_globals
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting water content
      real(kind=rkind) :: theta

      real(kind=rkind) :: a,n,m, theta_e
      type(integpnt_str) :: quadpnt_loc
      

      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
	quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen"
	  ERROR STOP
	end if
	h = x(1)
      end if
      
      
      
      a = vgset(layer)%alpha
      n = vgset(layer)%n
      m = vgset(layer)%m
      

      if (h >=0.0_rkind) then
	  theta = vgset(layer)%Ths
	  RETURN
      else
	  theta_e = 1/(1+(a*(abs(h)))**n)**m
	  theta = theta_e*(vgset(layer)%Ths-vgset(layer)%Thr)+vgset(layer)%Thr
      end if

    end function vangen
    
    function inverse_vangen(pde_loc, layer, quadpnt, x) result(hpress)
      use typy
      use re_globals
      use pde_objs
      use core_tools
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> water content
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: theta
      !> resulting pressure head
      real(kind=rkind) :: hpress
      
      
      real(kind=rkind) :: a,n,m
      type(integpnt_str) :: quadpnt_loc
      
 

      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
      	quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	theta = pde_loc%getval(quadpnt_loc)
      else
	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen"
	  ERROR STOP
	end if
	theta = x(1)
      end if
      
      
      
      a = vgset(layer)%alpha
      n = vgset(layer)%n
      m = vgset(layer)%m
      
      if (abs(theta - vgset(layer)%Ths) < epsilon(theta)) then
	hpress = 0
      else
	if (theta >  vgset(layer)%Ths + 10*epsilon(theta)) then
	  call write_log("theta is greater then theta_s, exiting")
	  print *, "called from re_constitutive::inverse_vangen"
	  error stop
	else if (theta < 0) then
	  call write_log("theta is negative strange, exiting")
	  print *, "called from re_constitutive::inverse_vangen"
	  error stop 
	end if
	hpress = ((((vgset(layer)%Ths - vgset(layer)%Thr)/(theta-vgset(layer)%Thr))**(1.0_rkind/m)-1) &  
	**(1.0_rkind/n))/(-a)
      end if
      
    end function inverse_vangen

    function vangen_tab(pde_loc, layer, quadpnt, x) result(theta)
      use typy
      use re_globals
      use pde_objs
      use core_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting water content
      real(kind=rkind) :: theta

      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist
      type(integpnt_str) :: quadpnt_loc      

      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_tab"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
      	quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if

      if (h<0) then
	if ( h/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	  
	  pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (pos <= ubound(watcontab,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    theta = (watcontab(layer,pos+1)-watcontab(layer,pos))/drutes_config%fnc_discr_length*dist + watcontab(layer,pos)
	  else
	    if (present(quadpnt)) theta = vangen(pde_loc, layer, quadpnt)
	    if (present(x)) theta = vangen(pde_loc, layer, x=x)
	    if (.not. tabwarning) then
	      call write_log(trim(tabmsg))
	      tabwarning = .true.
	    end if
	  end if
	else
	
	  if (.not. intwarning) then
	    call intoverflow()
	    intwarning = .true.
	  end if
	
	  if (present(quadpnt)) theta = vangen(pde_loc, layer, quadpnt)
	  if (present(x)) theta = vangen(pde_loc, layer, x=x)
	  
	end if
	
      else
	theta = vgset(layer)%Ths	
      end if
      
      

    end function vangen_tab


    !> \brief so-called retention water capacity, it is a derivative to retention curve function
    !! \f E(h) = C(h) + \frac{\theta(h)}{\theta_s}S_s \f]
    !! where
    !! \f[ C(h) = \left\{ \begin{array}{l l}\frac{m n \alpha  (-h \alpha )^{-1+n}}{\left(1+(-h \alpha )^n\right)^{1+m}}(\theta_s - \theta_r) ,  & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ 0, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f]
    !! and 
    !! \f[ \theta(h) = \left\{ \begin{array}{l l} \frac{\theta_s -\theta_r}{(1+(-\alpha h)^n_{vg})^m_{vg}} + \theta_r,  & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ \theta_S, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f]
    !<
    function vangen_elast(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use re_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

      real(kind=rkind) :: C, a, m, n, tr, ts 
      type(integpnt_str) :: quadpnt_loc      
          
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen_elast"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_elast"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  ERROR STOP
	end if
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_elast"
	  ERROR STOP
	end if
	h = x(1)
      end if

      if (h < 0) then
	a = vgset(layer)%alpha
	n = vgset(layer)%n
	m = vgset(layer)%m
	tr = vgset(layer)%Thr
	ts = vgset(layer)%Ths
	C = a*m*n*(-tr + ts)*(-(a*h))**(-1 + n)*(1 + (-(a*h))**n)**(-1 - m)
      else
	E = vgset(layer)%Ss
	RETURN
      end if

      E = C + vangen(pde_loc, layer, x=(/h/))/vgset(layer)%Ths*vgset(layer)%Ss
      

    end function vangen_elast


    function vangen_elast_tab(pde_loc, layer, quadpnt, x) result(E)
      use typy
      use re_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist
      type(integpnt_str) :: quadpnt_loc     
         
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen_elast_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_elast_tab"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_elast_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if

      if (h<0) then
	if ( h/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	  pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (pos <= ubound(warecatab,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    E = (warecatab(layer,pos+1)-warecatab(layer,pos))/drutes_config%fnc_discr_length*dist + warecatab(layer,pos)
	  else
	    if (present(quadpnt)) E = vangen_elast(pde_loc, layer, quadpnt)
	    if (present(x)) E = vangen_elast(pde_loc, layer, x=x)
	    if (.not. tabwarning) then
	      call write_log(trim(tabmsg))
	      tabwarning = .true.
	    end if
	  end if
	else
	 
	  if (.not. intwarning) then
	    call intoverflow()
	    intwarning = .true.
	  end if
	  if (present(quadpnt)) E = vangen_elast(pde_loc, layer, quadpnt_loc)
	  if (present(x)) E = vangen_elast(pde_loc, layer, x=x)	  
	
	end if 
      else
	E = vgset(layer)%Ss	
      end if


    end function vangen_elast_tab



    !> \brief Mualem's fucntion for unsaturated hydraulic conductivity with van Genuchten's water content substitution
    !! \f[   K(h) = \left\{ \begin{array}{l l} K_s\frac{\left( 1- (-\alpha h)^{n_{vg}m_{vg}} \left( 1+ (-\alpha h)^{n_{vg}} \right)^{-m_{vg}} \right)^2}{\left(1+(-\alpha h)^{n_{vg}} \right)^{\frac{m_{vg}}{2}}},  &  \mbox{$\forall$  $h \in$ $(-\infty,0)$}\\ K_s, & \mbox{$\forall$   $h \in$ $\langle 0, +\infty)$}\\ \end{array} \right. \f]
    !<
    subroutine mualem(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      real(kind=rkind) :: h
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar

      real(kind=rkind) :: a,n,m, tmp
      type(integpnt_str) :: quadpnt_loc
      
  

      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::mualem"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::mualem"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::mualem"
	  ERROR STOP
	end if
	h = x(1)
      end if
      
      
      if (h >= 0) then
	tmp = 1
      else
	a = vgset(layer)%alpha
	n = vgset(layer)%n
	m = vgset(layer)%m

	tmp =  (1 - (-(a*h))**(m*n)/(1 + (-(a*h))**n)**m)**2/(1 + (-(a*h))**n)**(m/2.0_rkind)
      end if
	
      if (present(tensor)) then
	tensor = tmp* vgset(layer)%Ks
      end if

      if (present(scalar)) then
	scalar = tmp
      end if
    end subroutine mualem


    subroutine mualem_tab(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use re_globals
      use pde_objs
      use core_tools
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor		
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar

      real(kind=rkind) :: h
      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist, tmp
      type(integpnt_str) :: quadpnt_loc      
      
   
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::mualem_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::mualem_tab"
	ERROR stop
      end if
      
      

      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::mualem_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if
      
      


      
      if (h<0) then
	if (-h/drutes_config%fnc_discr_length < 0.1*huge(1) ) then
	  pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (pos <= ubound(Ktab,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    tmp = (Ktab(layer,pos+1)-Ktab(layer,pos))/drutes_config%fnc_discr_length*dist + Ktab(layer,pos)
	  else
	    if (present(quadpnt)) call mualem(pde_loc, layer, quadpnt, scalar = tmp)
	    if (present(x)) call mualem(pde_loc, layer, x=x, scalar = tmp)
	    if (.not. tabwarning) then
	      call write_log(trim(tabmsg))
	      tabwarning = .true.
	    end if
	  end if
	else
	  if (.not. intwarning) then
	    call intoverflow()
	    intwarning = .true.
	  end if
	  if (present(quadpnt)) call mualem(pde_loc, layer, quadpnt_loc, scalar = tmp)
	  if (present(x)) call mualem(pde_loc, layer, x=x, scalar = tmp)
	end if 
	 
      else
	tmp = 1
      end if

      if (present(tensor)) then
	tensor = tmp* vgset(layer)%Ks
      end if

      if (present(scalar)) then
	scalar = tmp
      end if
    end subroutine mualem_tab


    !> \brief derivative of Mualem's function with van Genuchten's substitution for water content function
    !! \f[ \begin{split} \frac{\textrm{d} K(h)}{\textrm{d} h}  = \left\{ \begin{array}{l l} K_s \frac{1}{2} \alpha (-\alpha h)^{(-1 + n)} (1 + (-\alpha h)^n)^{(-1 - m/ 2)} (1 - (-\alpha h)^{(m n)} \\ (1 + (-\alpha h)^n)^{-m})^2 m n +  2 (1 + (-\alpha h)^n)^{(-m/2)} (1 - (-\alpha h)^{(m n)}  \\ (1 + (-\alpha h)^n)^{-m}) (-\alpha (-\alpha h)^{(-1 + n + m n)}  \\ (1 + (-\alpha h)^n)^{(-1 - m)} m n + \alpha (-\alpha h)^{(-1 + m n)} \\ (1 + (-a h)^n)^{-m} m n),  \quad  \mbox{$\forall$  $h \in$ $(-\infty,0)$}\\ 0, \quad \mbox{$\forall$   $h \in$ $\langle 0, +\infty)$}\\ \end{array} \right. \end{split} \f]
    !<
    subroutine dmualem_dh(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
      real(kind=rkind), dimension(:), intent(in), optional :: vector_in
      !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
      !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
      !<
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
      real(kind=rkind), intent(out), optional :: scalar

      real(kind=rkind) :: a,n,m, tmp, h
      type(integpnt_str) :: quadpnt_loc     
      
         
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::dmualem_dh"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::dmualem_dh"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::dmualem_dh"
	  ERROR STOP
	end if
	h = x(1)
      end if
      

      if (h < 0) then
	a = vgset(layer)%alpha
	n = vgset(layer)%n
	m = vgset(layer)%m
	tmp = (    (a*(-(a*h))**(-1 + n)*(1 + (-(a*h))**n)**(-1 - m/2.0)* &
	    (1 - (-(a*h))**(m*n)/(1 + (-(a*h))**n)**m)**2*m*n)/2.0 + &
	 (2*(1 - (-(a*h))**(m*n)/(1 + (-(a*h))**n)**m)* &
	    (-(a*(-(a*h))**(-1 + n + m*n)*(1 + (-(a*h))**n)**(-1 - m)*m*n) + &
	      (a*(-(a*h))**(-1 + m*n)*m*n)/(1 + (-(a*h))**n)**m))/ &
	  (1 + (-(a*h))**n)**(m/2.0))
      else
	tmp = 0
      end if

      if (present(vector_out)) then
	! must be negative, because the commnon scheme of the CDE problem has negative convection, but RE has positive convection
	vector_out = -vgset(layer)%Ks(drutes_config%dimen,:) * tmp
      end if

      if (present(scalar)) then
	scalar = tmp
      end if

    end subroutine dmualem_dh
    
    
    function derKz(pde_loc,layer, quadpnt, x) result(dKdz)
      use typy
      use re_globals
      use pde_objs
      use geom_tools
      use debug_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: dKdz
      real(kind=rkind), dimension(3,3) :: Ktmp
      
      
      real(kind=rkind) :: dKdx
      real(kind=rkind), dimension(:,:), allocatable, save :: pts
      type(integpnt_str) :: quadpnt_loc
      integer(kind=ikind) :: layer_loc, vecino, i, D   


      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen_elast"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_elast"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  ERROR STOP
	end if
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_elast"
	  ERROR STOP
	end if
	h = x(1)
      end if
      
      
      if (.not. allocated(pts) ) then
	allocate(pts(drutes_config%dimen+1, ubound(elements%data,2)))
      end if
      
      quadpnt_loc%type_pnt = "ndpt"
      quadpnt_loc%element = quadpnt%element
      quadpnt_loc%column = 2
      D = drutes_config%dimen
     
      do i=1, ubound(elements%data,2)
	quadpnt_loc%order = elements%data(quadpnt%element,i)
	vecino = elements%neighbours(quadpnt%element,i)
	if (vecino > 0) then
	  layer_loc = elements%material(vecino,1)
	else
	  layer_loc = elements%material(quadpnt%element,1)
	end if 
	quadpnt_loc%order = elements%data(quadpnt%element, i)
	call pde_loc%pde_fnc(1)%dispersion(pde_loc, layer_loc, quadpnt_loc, tensor=Ktmp(1:D, 1:D))
	pts(D+1, i) = Ktmp(D,D)
	pts(1:D,i) = nodes%data(elements%data(quadpnt%element,i),1:D)
      end do
      
      
      
      select case(drutes_config%dimen)
	case(1)
	  if (pts(1,2)-pts(1,1) > 0) then
	    dKdz = (pts(2,2)-pts(2,1))/(pts(1,2)-pts(1,1))
	  else
	    dKdz = (pts(2,1)-pts(2,2))/(pts(1,1)-pts(1,2))
	  end if
	case(2)
	  call plane_derivative(pts(1,:), pts(2,:), pts(3,:), dKdx, dKdz)
      end select
      
      dKdz = -dKdz
	
      
    
    end function derKz


    !> \brief convection term for the Richards equation for axisymmetric problems
    !! the vertical convection is equal as for the Richards equatio in standard coordinates, the horizontal
    !! convection equals
    !! [\ \frac{1}{r}K(h)
    !<
    subroutine convection_rerot(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use geom_tools
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
      real(kind=rkind), dimension(:), intent(in), optional :: vector_in
      !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
      !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
      !<
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
      real(kind=rkind), intent(out), optional :: scalar
      
      real(kind=rkind), dimension(2) :: coord

      real(kind=rkind) :: a,n,m, tmp, h, K, flux
      integer :: proc
      type(integpnt_str) :: quadpnt_loc    
      

      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::dmualem_dh"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::dmualem_dh"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::dmualem_dh"
	  ERROR STOP
	end if
	h = x(1)
      end if
      

      select case (drutes_config%name)
	case("RE_rot")
	  if (drutes_config%fnc_method == 0) then
	    call dmualem_dh(pde_loc, layer, x=(/h/), vector_out = vector_out)
	  else
	    call dmualem_dh_tab(pde_loc, layer, x=(/h/), vector_out = vector_out)
	  end if
	case("RErotH")
	  vector_out = 0
      end select
      
      tmp = vector_out(2)
      
      call getcoor(quadpnt, coord)
      
      call pde_loc%flux(layer, quadpnt, vector_out=vector_out)
      
      vector_out(2) = tmp
      
      if (abs(coord(2)) < -12.5 ) then
	vector_out(1) =  1.0_rkind/coord(1)*vector_out(1)	
      else
	vector_out(1) = 0
      end if

    end subroutine convection_rerot
    
    function sinkterm(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      
      val = vgset(layer)%sinkterm
     
    
    
    end function sinkterm
    
    
    subroutine dmualem_dh_tab(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use pde_objs
      use re_globals
      use core_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
      real(kind=rkind), dimension(:), intent(in), optional :: vector_in
      !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  relative unsaturated 
      !!hydraulic conductivity derivative in respect to h (scalar value)     
      !<
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative unsaturated hydraulic conductivity derivative in respect to h
      real(kind=rkind), intent(out), optional :: scalar


      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist, tmp, h
      type(integpnt_str) :: quadpnt_loc     
      
    
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::dmualem_dh_tab"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::dmualem_dh_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if      
      
      if (h<0) then
	if (h/drutes_config%fnc_discr_length < 0.1*huge(1) ) then
	
	  pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (pos <= ubound(dKdhtab,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    tmp = (dKdhtab(layer,pos+1)-dKdhtab(layer,pos))/drutes_config%fnc_discr_length*dist + dKdhtab(layer,pos)
	  else
	    if (present(quadpnt)) call dmualem_dh(pde_loc, layer, quadpnt, scalar=tmp)
	    if (present(x)) call dmualem_dh(pde_loc, layer, x=x, scalar=tmp)
	    if (.not. tabwarning) then
	      call write_log(trim(tabmsg))
	      tabwarning = .true.
	    end if
	  end if
	else
          if (.not. intwarning) then
	    call intoverflow()
	    intwarning = .true.
	  end if
	  if (present(quadpnt))call dmualem_dh(pde_loc, layer, quadpnt_loc, scalar=tmp)
	  if (present(x))call dmualem_dh(pde_loc, layer, x=x, scalar=tmp)
	end if
	  
      else
	tmp = 0
      end if
      

      if (present(vector_out)) then
	! must be negative, because the commnon scheme of the CDE problem has negative convection, but RE has positive convection
	vector_out = -vgset(layer)%Ks(drutes_config%dimen,:) * tmp
      end if

      if (present(scalar)) then
	scalar = tmp
      end if

    end subroutine dmualem_dh_tab


    subroutine darcy_law(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length

      real(kind=rkind), dimension(3,3)  :: K
      integer                           :: D
      integer(kind=ikind)               :: i
      integer(kind=ikind), dimension(3) :: nablaz
      real(kind=rkind), dimension(3)  :: gradH
      real(kind=rkind), dimension(3)  :: vct
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable :: gradient
      type(integpnt_str) :: quadpnt_loc
      

      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
	print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
	ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::darcy_law"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
	call pde_loc%getgrad(quadpnt, gradient)
      else
        if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::darcy_law"
	  ERROR STOP
	end if
	h = x(1)
	allocate(gradient(ubound(grad,1)))
	gradient = grad
      end if
      
      D = drutes_config%dimen

      nablaz = 0
      nablaz(D) = 1
      
      gradH(1:D) = gradient(1:D) + nablaz(1:D)

      call pde_loc%pde_fnc(1)%dispersion(pde_loc, layer, x=(/h/), tensor=K(1:D, 1:D))
     
      
      vct(1:D) = matmul(-K(1:D,1:D), gradH(1:D))


      if (present(flux_length)) then
        select case(D)
          case(1)
                flux_length = vct(1)
          case(2)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
          case(3)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
        end select
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if

    end subroutine darcy_law




  !> creates a table of values of constitutive functions for the Richards equation to be linearly approximated 
    subroutine tabvalues(pde_loc, Kfnc, dKdhfnc, Cfnc, thetafnc)
      use typy
      use globals
      use pde_objs
      use printtools
      use core_tools

      class(pde_str), intent(in) :: pde_loc
      interface
	subroutine Kfnc(pde_loc, layer, quadpnt, x, tensor, scalar)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in)                           :: layer
	  type(integpnt_str), intent(in), optional                   :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional      :: x
	  real(kind=rkind), intent(out), dimension(:,:), optional   :: tensor
	  real(kind=rkind), intent(out), optional                   :: scalar 
	end subroutine Kfnc
      end interface

      interface
	subroutine dKdhfnc(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in) 		            :: layer
	  type(integpnt_str), intent(in), optional                  :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional      :: x
          real(kind=rkind), dimension(:), intent(in), optional      :: vector_in
          real(kind=rkind), dimension(:), intent(out), optional     :: vector_out
	  real(kind=rkind), intent(out), optional                   :: scalar 
	end subroutine dKdhfnc
      end interface

      interface
	function Cfnc(pde_loc, layer, quadpnt,  x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function Cfnc
      end interface


      interface
	function thetafnc(pde_loc, layer, quadpnt, x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc	  
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function thetafnc
      end interface


      integer(kind=ikind) :: i,j, n, k, domain_count	
      integer :: l
      integer(kind=ikind) :: maxcalls, counter
      character(len=1), dimension(2) :: domain
      real(kind=rkind) :: dx

      n = int(maxpress/drutes_config%fnc_discr_length)+1
      
      drutes_config%fnc_discr_length = 1.0_rkind*maxpress/n

      dx = drutes_config%fnc_discr_length

      domain(1) = "m"

      domain(2) = "f"


     domain_count = 1


      allocate(Ktab_all(domain_count, ubound(vgmatrix,1), n))
      allocate(warecatab_all(domain_count, ubound(vgmatrix,1),n))
      allocate(watcontab_all(domain_count, ubound(vgmatrix,1), n))
      allocate(dKdhtab_all(domain_count, ubound(vgmatrix,1), n))


      do k=1, domain_count
	call domainswitch(domain(k))
	maxcalls = ubound(vgset,1)*n
	counter = maxcalls
	call write_log(text="creating constitutive function table for domain: ", text2=domain(k))
	do i=1, ubound(vgset,1)
	  do j=1, n
	    if (this_image() == 1) then
	      counter = counter - 1
	      l = 100*(maxcalls - counter)/maxcalls
	      call progressbar(l)
	    end if
	    call Kfnc(pde_loc,i, x=(/-(j-1)*dx/), scalar=Ktab(i,j))
! 	    Ktab(i,j) = Kfnc(i, -(j-1)*dx)
! 	    dKdhtab(i,j) = dKdhfnc(i, -(j-1)*dx)
	    call dKdhfnc(pde_loc, i, x=(/-(j-1)*dx/), scalar=dKdhtab(i,j))
	    warecatab(i,j) = Cfnc(pde_loc, i, x=(/-(j-1)*dx/))
	    watcontab(i,j) = thetafnc(pde_loc, i, x=(/-(j-1)*dx/))
	  end do
	end do
      end do



    end subroutine tabvalues


    subroutine domainswitch(domain)
      use re_globals
 
      character(len=1), intent(in) :: domain


      if (domain == "m") then
	vgset => vgmatrix
	Ktab => Ktab_all(1,:,:)
	dKdhtab => dKdhtab_all(1,:,:)
	warecatab => warecatab_all(1,:,:)
	watcontab => watcontab_all(1,:,:)
      else if (domain == "f") then
	vgset => vgfractures
	Ktab => Ktab_all(2,:,:)
	dKdhtab => dKdhtab_all(2,:,:)
	warecatab => warecatab_all(2,:,:)
	watcontab => watcontab_all(2,:,:)
      else
	ERROR STOP "runtime error, invalid argument in RE_constitutive::domainswitch() procedure"
      end if
	


    end subroutine domainswitch


    subroutine init_zones(vg)
      use typy
      use globals
      use re_globals
      use geom_tools
      use printtools
      use debug_tools
      
      
      
      type(soilpar), intent(in out), dimension(:) :: vg
      integer(kind=ikind) :: i, j, k, l
      integer :: i_err
      logical :: increased, overzero
      real(kind=rkind) :: theta, htest, htest2, h, hprev, tmp1, tmp2, m, inflex_point
      integer, dimension(:,:), allocatable :: images_run
      real(kind=rkind), dimension(:,:,:), allocatable :: cotable[:]
      character(len=7) :: procentas
      integer(kind=ikind) :: n
      logical :: backward

      n = int(maxpress/drutes_config%fnc_discr_length)+1    



      do i=1, ubound(vg,1)
	if (vg(i)%rcza_set%use) then
         allocate(vg(i)%rcza_set%tab(n, 3))
	 vg(i)%rcza_set%tab = 0
	end if
      end do

      i_err = 0
      htest = 0.0
    

	
      allocate(images_run(NUM_IMAGES(),2))
      
      images_run(1,1) = 1
      images_run(1,2) = int(ubound(vgset,1)/NUM_IMAGES()) + modulo(ubound(vgset,1), NUM_IMAGES())

      allocate(cotable(images_run(1,2), n, 3)[*])

!       do i=2, NUM_IMAGES()
! 	images_run(i,1) = images_run(i-1,2) + 1
! 	if (images_run(i,1) + int(ubound(vgset,1)/NUM_IMAGES()) <= ubound(vgset,1) ) then
! 	  images_run(i,2) = images_run(i,1) + int(ubound(vgset,1)/NUM_IMAGES())
! 	else
! 	  images_run(i,2) = images_run(i,1)
! 	end if
!       end do


    

      layers: do j=images_run(THIS_IMAGE(),1), images_run(THIS_IMAGE(),2)
		if (vgset(j)%rcza_set%use) then
		  k = j - images_run(THIS_IMAGE(),1) + 1
		  !backward
		  write(unit=terminal, fmt=*) "preparing layer:", j
		  this_layer = j
		  call solve_bisect(-epsilon(h), -maxpress, secondder_retc, 0.001_rkind, i_err, inflex_point)

		  print *, "inflex point of the layer:", j, "found at:", inflex_point
		  

		  levels: do i=1, ubound(cotable,2)
			    if (this_image() == 1) then
			      call progressbar(int(50*i/ubound(cotable,2)))
			    end if
			    h = -drutes_config%fnc_discr_length*i
			    if (i < 2) then
			      htest = -drutes_config%fnc_discr_length
			      hprev = -drutes_config%fnc_discr_length
			    else
			      htest = cotable(k,i-1,2)[THIS_IMAGE()]
			      hprev = htest
			    end if
			    rough: do
				    if (zone_error_val(h, htest, vgset(j)%rcza_set%val) > 0) then
						  cotable(k,i,1)[THIS_IMAGE()] = h
						  cotable(k,i,2)[THIS_IMAGE()] = htest
						  exit rough
				    else
				      htest = htest - drutes_config%fnc_discr_length
				      if (htest < -2*maxpress) then
					cotable(k,i:ubound(cotable,2),2)[THIS_IMAGE()] = htest
					do l= i, ubound(cotable,2)
					  cotable(k,l,1)[THIS_IMAGE()] = -drutes_config%fnc_discr_length*l
					end do
					exit levels
				      end if
				    end if
			    end do rough
		  end do levels

    ! 		  !forward
		  levels2: do i=1, ubound(cotable,2)
			    if (this_image() == 1) then
			      call progressbar(int(50*i/ubound(cotable,2)+50))
			    end if
			    h = -drutes_config%fnc_discr_length*i
			    if (i < 2) then
			      htest = -drutes_config%fnc_discr_length
			      hprev = -drutes_config%fnc_discr_length
			      backward = .false.
			    else
			      if (cotable(k,i-1,3)[THIS_IMAGE()] > inflex_point-5*drutes_config%fnc_discr_length) then
				htest = h
				hprev = h
				backward = .false.
			      else
				htest = cotable(k,i-1,3)[THIS_IMAGE()]
				hprev = htest
				backward = .true.
			      end if
			    end if
			    rough2: do
				      if ((zone_error_val(h, htest, vgset(j)%rcza_set%val) > epsilon(0.0_rkind) .and. .not.(backward)) .or. &
					(zone_error_val(h, htest, vgset(j)%rcza_set%val) < 0 .and. backward))      then
				      cotable(k,i,3)[THIS_IMAGE()] = htest
				      exit rough2
				    else
					if (htest < 0) then
					  if (backward) then
					    l = -1
					  else 
					    l = 1
					  end if
					  htest = htest + l*drutes_config%fnc_discr_length
					else
					  cotable(k,i,3)[THIS_IMAGE()] = huge(1.0_rkind)
					  exit rough2
					end if
				      end if
				    end do rough2
		  end do levels2 
		  if (this_image() == 1) then
		    print *, " "
		    call flush(6)
		  end if
	  end if
	end do layers

	! sync data
	sync all
	
	  
	do i=1, ubound(images_run,1)
	  do j=images_run(i,1), images_run(i,2)
	    vg(j)%rcza_set%tab = cotable(j,:,:)[i]
	  end do
	end do



	 sync all 




	

	deallocate(images_run)
	deallocate(cotable)




      end subroutine init_zones



      
      
      function zone_error_val( h, htest, taylor_max) result(err_val)
	use typy
	use pde_objs

	!> starting point
	real(kind=rkind), intent(in) :: h
	!> point to be evaluated
	real(kind=rkind), intent(in) :: htest
	!> maximal error of taylor series
	real(kind=rkind), intent(in) :: taylor_max
	!> error value of first order taylor approximation starting at point h and evaluating at point htest
	real(kind=rkind) :: err_val

	err_val=abs(vangen(pde(1), this_layer, x=(/h/)) + vangen_elast(pde(1), this_layer, x=(/htest/))*(htest-h)- &
		vangen(pde(1), this_layer, x=(/htest/))) - taylor_max





      end function zone_error_val


      !second order derivative of the van Genuchten retention curve function
      pure &
      function secondder_retc(h) result(curvature)
	use typy
	use re_globals

	!> input pressure head, van genuchten's parameters are obtained from variable van_gen_coeff
	real(kind=rkind), intent(in) :: h
	!> resulting curvature
	real(kind=rkind) :: curvature

	!!!!--local variables------
	real(kind=rkind) :: a, n, m, tr, ts

	a = vgset(this_layer)%alpha
	n = vgset(this_layer)%n
	m = vgset(this_layer)%m
	tr = vgset(this_layer)%Thr
	ts = vgset(this_layer)%Ths


	if( h < 0.0_rkind) then
	  curvature = -(a**2*(-(a*h))**(-2 + n)*(1 + (-(a*h))**n)**(-1 - m)*m*(-1 + n)*n* &
	      (-tr + ts)) - a**2*(-(a*h))**(-2 + 2*n)*(1 + (-(a*h))**n)**(-2 - m)*(-1 - m)* &
	    m*n**2*(-tr + ts)
	else
	  curvature = 0.0_rkind
	end if

	

      end function secondder_retc


      function rcza_check_old() result(passed)
	use typy
	use globals
	use global_objs
        use re_globals
        use pde_objs

        integer(kind=ikind) :: process
	logical :: passed
	integer(kind=ikind) :: i, j, position
	real(kind=rkind) :: low, high, hprev, h

        process = pde_common%current_proc

        do i=1, ubound(vgset,1)
          if (vgset(i)%rcza_set%use) then
            do j=pde(process)%procbase_fnc(1), pde(process)%procbase_fnc(2)
              hprev = pde_common%xvect(j,1)
              h = pde_common%xvect(j,3)
              if ( h<0 ) then
                position = int(-hprev/drutes_config%fnc_discr_length)+1
                low = vgset(i)%rcza_set%tab(position,2)
                high = vgset(i)%rcza_set%tab(position,3)
!                 print *, low, h, high
                if (h < low .or. h > high) then
                  passed = .false.
                  RETURN
                end if
              end if
            end do
          end if
        end do


	passed = .true.


      end function rcza_check_old
      
      function rcza_check(pde_loc) result(passed)
	use typy
	use globals
	use global_objs
        use re_globals
        use pde_objs

        class(pde_str),intent(in) :: pde_loc
        integer(kind=ikind) :: el, pos, mat
	logical :: passed
	integer(kind=ikind) :: i, position, nd
	real(kind=rkind) :: low, high, hprev, h
	type(integpnt_str) :: quadpnt
	
        quadpnt%type_pnt = "ndpt"
               

	do i=1, ubound(pde_loc%solution,1)
	  el = nodes%element(i)%data(1)
	  mat = elements%material(el,1)
	    
	  quadpnt%order = i
	  quadpnt%column = 1
	  quadpnt%preproc=.true.
	  hprev = pde_loc%getval(quadpnt)
	  quadpnt%column = 2
	  h = pde_loc%getval(quadpnt)
	  if (vgset(mat)%rcza_set%use) then
	    if ( hprev<0 ) then
	      position = int(-hprev/drutes_config%fnc_discr_length)+1
	      low = vgset(mat)%rcza_set%tab(position,2)
	      high = vgset(mat)%rcza_set%tab(position,3)
	      if (h < low .or. h > high) then
		passed = .false.
		RETURN
	      end if
	      end if
	   end if
	 end do

 	passed = .true.


      end function rcza_check


      subroutine re_dirichlet_bc(pde_loc, el_id, node_order, value, code) 
	use typy
	use globals
	use global_objs
	use pde_objs

	class(pde_str), intent(in) :: pde_loc
	integer(kind=ikind), intent(in)  :: el_id, node_order
	real(kind=rkind), intent(out), optional    :: value
	integer(kind=ikind), intent(out), optional :: code

	integer(kind=ikind) :: edge_id, i, j

	
	edge_id = nodes%edge(elements%data(el_id, node_order))

	if (present(value)) then
	  if (pde_loc%bc(edge_id)%file) then
	    do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	      if (pde_loc%bc(edge_id)%series(i,1) > time) then
		if (i > 1) then
		  j = i-1
		else
		  j = i
		end if
		value = pde_loc%bc(edge_id)%series(j,2)
		EXIT
	      end if
	    end do
	  else
	    value = pde_loc%bc(edge_id)%value
	  end if
	end if

	if (present(code)) then
	  code = 1
	end if
      end subroutine re_dirichlet_bc


      subroutine re_null_bc(pde_loc, el_id, node_order, value, code) 
	use typy
	use globals
	use global_objs
	use pde_objs
	
	class(pde_str), intent(in) :: pde_loc
	integer(kind=ikind), intent(in)  :: el_id, node_order
	real(kind=rkind), intent(out), optional    :: value
	integer(kind=ikind), intent(out), optional :: code

	if (present(value)) then
	  value = 0.0_rkind
	end if

	if (present(code)) then
	  code = 2
	end if
	
      end subroutine re_null_bc



      subroutine re_dirichlet_height_bc(pde_loc, el_id, node_order, value, code) 
	use typy
	use globals
	use global_objs
	use pde_objs
	use debug_tools
	
	class(pde_str), intent(in) :: pde_loc
	integer(kind=ikind), intent(in)  :: el_id, node_order
	real(kind=rkind), intent(out), optional    :: value
	integer(kind=ikind), intent(out), optional :: code

	
	integer(kind=ikind) :: edge_id, i, j
	real(kind=rkind) :: tempval, node_height
	

	if (present(value)) then
	  edge_id = nodes%edge(elements%data(el_id, node_order))
	  node_height = nodes%data(elements%data(el_id, node_order), drutes_config%dimen)

	  if (pde_loc%bc(edge_id)%file) then
	    do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	      if (pde_loc%bc(edge_id)%series(i,1) > time) then
		if (i > 1) then
		  j = i-1
		else
		  j = i
		end if
		tempval = pde_loc%bc(edge_id)%series(j,2)
		EXIT
	      end if
	    end do
	  else
	    tempval =  pde_loc%bc(edge_id)%value
	  end if

	  
	  value = tempval - node_height
	end if

	
	if (present(code)) then
	  code = 1
	end if
	
      end subroutine re_dirichlet_height_bc





      subroutine re_neumann_bc(pde_loc, el_id, node_order, value, code) 
	use typy
	use globals
	use global_objs
	use pde_objs

	class(pde_str), intent(in) :: pde_loc
	integer(kind=ikind), intent(in)  :: el_id, node_order
	real(kind=rkind), intent(out), optional    :: value
	integer(kind=ikind), intent(out), optional :: code
	
	real(kind=rkind), dimension(3,3) :: K

	integer(kind=ikind) :: i, edge_id, j
	real(kind=rkind), dimension(3) :: gravflux, bcflux
	real(kind=rkind) :: bcval, gfluxval
	integer :: i1

      
	if (present(value)) then
	  edge_id = nodes%edge(elements%data(el_id, node_order))

	  i = pde_loc%permut(elements%data(el_id, node_order))

	  call pde_loc%pde_fnc(1)%dispersion(pde_loc, elements%material(el_id,1), x=(/pde_common%xvect(i,2)/), &
			      tensor=K(1:drutes_config%dimen, 1:drutes_config%dimen))

	  gravflux(1:drutes_config%dimen) = K(drutes_config%dimen, 1:drutes_config%dimen)*elements%nvect_z(el_id, node_order)
	  

	  if (pde_loc%bc(edge_id)%file) then
	    do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	      if (pde_loc%bc(edge_id)%series(i,1) > time) then
		if (i > 1) then
		  j = i-1
		else
		  j = i
		end if
		bcval = pde_loc%bc(edge_id)%series(j,2)
		EXIT
	      end if
	    end do
	  else
	    bcval = pde_loc%bc(edge_id)%value
	  end if
	  

	  


	  select case(drutes_config%dimen)
	    case(1)
	      value = bcval - gravflux(1)
	    case(2)
	      bcflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*bcval
	      bcflux(2) = elements%nvect_z(el_id, node_order)*bcval
	      bcflux = bcflux + gravflux
	      value = sqrt(bcflux(1)*bcflux(1) + bcflux(2)*bcflux(2))
  ! 	    print *, value, gravflux, el_id, elements%data(el_id,:) ; stop 
	  end select
	end if
	
       if (present(code)) then
	  code = 2
	end if

      end subroutine re_neumann_bc
      
      
    subroutine re_atmospheric(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      
      real(kind=rkind), dimension(3,3) :: K
      
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind) :: theta, bcval
      integer(kind=ikind) :: i, edge_id, j
  
  
      if (present(code)) then
	code = 2
      end if
      
      if (present(value)) then
	edge_id = nodes%edge(elements%data(el_id, node_order))

	i = pde_loc%permut(elements%data(el_id, node_order))

	call pde_loc%pde_fnc(1)%dispersion(pde_loc, elements%material(el_id,1), x=(/pde_common%xvect(i,2)/), &
			      tensor=K(1:drutes_config%dimen, 1:drutes_config%dimen))

        gravflux(1:drutes_config%dimen) = K(drutes_config%dimen, 1:drutes_config%dimen)*elements%nvect_z(el_id, node_order)
	

	if (pde_loc%bc(edge_id)%file) then
	  do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	    if (pde_loc%bc(edge_id)%series(i,1) > time) then
	      if (i > 1) then
		j = i-1
	      else
		j = i
	      end if
	      bcval = pde_loc%bc(edge_id)%series(j,2)
	      EXIT
	    end if
	  end do
	else
	  bcval = pde_loc%bc(edge_id)%value
	end if
	
	if (bcval < 0) then
	  quadpnt%type_pnt = "ndpt"
	  quadpnt%order = elements%data(el_id,node_order)
	  layer = elements%material(el_id,1)
	  theta =  pde_loc%mass(layer, quadpnt)
	  bcval = bcval*(theta*theta)**(1.0_rkind/3.0_rkind)
	end if
	  
	  
	 select case(drutes_config%dimen)
	    case(1)
	      value = bcval - gravflux(1)
	    case(2)
	      bcflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*bcval
	      bcflux(2) = elements%nvect_z(el_id, node_order)*bcval
	      bcflux = bcflux + gravflux
	      value = sqrt(bcflux(1)*bcflux(1) + bcflux(2)*bcflux(2))
  ! 	    print *, value, gravflux, el_id, elements%data(el_id,:) ; stop 
	  end select 
	
	

      end if
      
      
    end subroutine re_atmospheric


      subroutine re_initcond(pde_loc) 
        use typy
        use globals
        use global_objs
        use pde_objs
        use re_globals
        
        class(pde_str), intent(in out) :: pde_loc
        integer(kind=ikind) :: i, j, k,l, m, layer, D
        real(kind=rkind) :: value
        
        D = drutes_config%dimen
  
        
        do i=1, elements%kolik
          layer = elements%material(i,1)
          do j=1, ubound(elements%data,2)
            k = elements%data(i,j)
            l = nodes%edge(k)
            m = pde_loc%permut(k)
            if (m == 0) then
	      call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
              pde_loc%solution(k) =  value 
            else
              select case (vgset(layer)%icondtype)
		case("H_tot")
		  pde_loc%solution(k) = vgset(layer)%initcond - nodes%data(k,D)
		case("hpres")
		  pde_loc%solution(k) = vgset(layer)%initcond 
		case("theta")
		  value = inverse_vangen(pde_loc, layer, x=(/vgset(layer)%initcond/))
		  pde_loc%solution(k) = value
	      end select
            end if
          end do   
        end do
        
        



      end subroutine re_initcond
      
      subroutine setmatflux()
	use typy
	use globals
        use global_objs
        use pde_objs
        use debug_tools
        
        integer(kind=ikind) :: i, j, vecino, id, k, l, nd, ndvec
        
        id = maxval(nodes%edge) + 100
        do i=1, elements%kolik
	  do j=1, ubound(elements%neighbours,2)
	    vecino = elements%neighbours(i,j)
	    if (vecino /= 0) then
	      if (elements%material(i,1) /= elements%material(vecino,1)) then
		do k=1, ubound(elements%data,2)
		  nd = elements%data(i, k)
		  do l=1, ubound(elements%data,2)
		    ndvec = elements%data(vecino, l)
		    if (nd == ndvec) then
		      nodes%edge(nd) = -id
		    end if
		  end do
		end do
	      end if
	    end if
	  end do
	end do
	

      
      end subroutine setmatflux



end module RE_constitutive
