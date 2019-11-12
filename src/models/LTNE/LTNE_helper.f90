module LTNE_helper
  use pde_objs
  use typy
  use LTNE_globs
  use debug_tools
  use RE_constitutive
  use re_globals

  public :: iceswitch, rho_icewat, Q_reduction, hl, thetai, thetal
  public:: vangen_ltne, mualem_ltne, temp_l_initcond, temp_s_initcond, wat_init, getval_retotltne
      
      
  
  contains

      !> switch for freezing condition based on Clapeyron equation 
    function iceswitch(quadpnt) result(sw)
      use typy
      use global_objs
      use re_globals
      
      type(integpnt_str), intent(in) :: quadpnt
      logical :: sw
      
      real(kind=rkind) :: Tf
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.
      
      if(clap) then
        Tf = Tref*exp(pde(wat)%getval(quadpnt_loc)*grav/Lf)
        Tf = Tf - 273.15_rkind
      else
        Tf = 0
      end if
      


      if (pde(heat_proc)%getval(quadpnt_loc) > Tf) then
      !> melting
        sw = .FALSE.
      else
      !> freezing
        sw = .TRUE.
      end if
          
    end function iceswitch

    function rho_icewat(quadpnt) result(rho)

      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: rho
      
      integer(kind=ikind) :: layer, el
      real(kind=rkind) :: thl, thall, thice
      
      
      if (quadpnt%type_pnt == "ndpt" ) then
        el = nodes%element(quadpnt%order)%data(1)
      else
        el = quadpnt%element
      end if
      
      layer = elements%material(el)
      thl = vangen_ltne(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      thall = vangen_ltne(pde(wat), layer, quadpnt)
      thice = thall - thl
      rho = (thl * rho_wat + thice * rho_ice)/thall
       
    end function rho_icewat
    
    function Q_reduction(layer, quadpnt, x) result(val)

      use typy
      use global_objs
      use re_globals
      
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      integer(kind=ikind) :: el
      real(kind=rkind) :: thl, thall, thice, val

      if(present(quadpnt)) then
        thall = vangen_ltne(pde(wat), layer, quadpnt)
        thl = vangen_ltne(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      end if
      if(present(x)) then
        thall = vangen_ltne(pde(wat), layer,x = x)
        thl = vangen_ltne(pde(wat), layer, x = x)
      end if
      thice = thall - thl
      val = thice/(thall- LTNE_par(layer)%Thr)
       
    end function Q_reduction
    
    
    
    function hl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use LTNE_globs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val, T_melt
      real(kind=rkind) :: hw, temp
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.

      hw = pde(wat)%getval(quadpnt_loc)
      temp = pde(heat_proc)%getval(quadpnt)
      T_melt = Tref*exp(hw*grav/Lf)
      
      if(iceswitch(quadpnt)) then
        val = hw+Lf/grav*log((temp+273.15_rkind)/T_melt) !units
      else
        val = hw
      end if
    end function hl
    
    function thetai(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      
      real(kind=rkind) :: thl, thall
      
      thl = vangen_ltne(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      thall = vangen_ltne(pde(wat), layer, quadpnt)
      
      !val = (thall * rho_icewat(quadpnt) - thl * rho_wat)/rho_ice
      val = thall - thl
    end function thetai
    
    function thetal(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      
      val = vangen_ltne(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
    end function thetal
    
    

    
    
    function vangen_ltne(pde_loc, layer, quadpnt, x) result(theta)
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
        print *, "exited from LTNE_helper::vangen_ltne"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from LTNE_helper::vangen_ltne"
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
          print *, "exited from LTNE_helper::vangen_ltne"
          ERROR STOP
        end if
        h = x(1)
      end if
      
      
      
      a = LTNE_par(layer)%alpha
      n = LTNE_par(layer)%n
      m = LTNE_par(layer)%m
      

      if (h >=0.0_rkind) then
        theta = LTNE_par(layer)%Ths
        RETURN
      else
        theta_e = 1/(1+(a*(abs(h)))**n)**m
        theta = theta_e*(LTNE_par(layer)%Ths-LTNE_par(layer)%Thr)+LTNE_par(layer)%Thr
      end if

    end function vangen_ltne
    
    
    
        !> \brief so-called retention water capacity, it is a derivative to retention curve function
    !! \f E(h) = C(h) + \frac{\theta(h)}{\theta_s}S_s \f]
    !! where
    !! \f[ C(h) = \left\{ \begin{array}{l l}\frac{m n \alpha  (-h \alpha )^{-1+n}}{\left(1+(-h \alpha )^n\right)^{1+m}}(\theta_s - \theta_r) ,  & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ 0, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f]
    !! and 
    !! \f[ \theta(h) = \left\{ \begin{array}{l l} \frac{\theta_s -\theta_r}{(1+(-\alpha h)^n_{vg})^m_{vg}} + \theta_r,  & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ \theta_S, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f]
    !<
    function vangen_elast_ltne(pde_loc,layer, quadpnt, x) result(E)
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
        print *, "exited from LTNE_helper::vangen_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from LTNE_helper::vangen_elast"
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
        print *, "exited from LTNE_helper::vangen_elast"
        ERROR STOP
      end if
        h = x(1)
      end if

      if (h < 0) then
        a = LTNE_par(layer)%alpha
        n = LTNE_par(layer)%n
        m = LTNE_par(layer)%m
        tr = LTNE_par(layer)%Thr
        ts = LTNE_par(layer)%Ths
        C = a*m*n*(-tr + ts)*(-(a*h))**(-1 + n)*(1 + (-(a*h))**n)**(-1 - m)
      else
        E = 0
        RETURN
      end if

      E = C 
      

    end function vangen_elast_ltne
    
    
    
    
    !> \brief Mualem's fucntion for unsaturated hydraulic conductivity with van Genuchten's water content substitution
    !! \f[   K(h) = \left\{ \begin{array}{l l} K_s\frac{\left( 1- (-\alpha h)^{n_{vg}m_{vg}} \left( 1+ (-\alpha h)^{n_{vg}} \right)^{-m_{vg}} \right)^2}{\left(1+(-\alpha h)^{n_{vg}} \right)^{\frac{m_{vg}}{2}}},  &  \mbox{$\forall$  $h \in$ $(-\infty,0)$}\\ K_s,  \mbox{$\forall$   $h \in$ $\langle 0, +\infty)$}\\ \end{array} \right. \f]
    !<
    subroutine mualem_ltne(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use LTNE_globs
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
        a = LTNE_par(layer)%alpha
        n = LTNE_par(layer)%n
        m = LTNE_par(layer)%m

        tmp =  (1 - (-(a*h))**(m*n)/(1 + (-(a*h))**n)**m)**2/(1 + (-(a*h))**n)**(m/2.0_rkind)
      end if
	
      if (present(tensor)) then
        tensor = tmp* LTNE_par(layer)%Ks
      end if

      if (present(scalar)) then
        scalar = tmp
      end if
    end subroutine mualem_ltne
    
    subroutine wat_init(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use geom_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (LTNE_par(1_ikind)%icondtypeRE)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/LTNE/hini.in", correct_h = .true.)
      end select
      
      D = drutes_config%dimen
      do i=1, elements%kolik
        layer = elements%material(i)
        do j=1, ubound(elements%data,2)
          k = elements%data(i,j)
          l = nodes%edge(k)
          m = pde_loc%permut(k)
          if (m == 0) then
            call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
            pde_loc%solution(k) =  value 
          else
            select case (LTNE_par(layer)%icondtypeRE)
              case("H_tot")
                pde_loc%solution(k) = LTNE_par(layer)%initcond !+ nodes%data(k,1)
              case("hpres")
                pde_loc%solution(k) = LTNE_par(layer)%initcond + nodes%data(k,D)
              case("theta")
                value = inverse_vangen_ltne(pde_loc, layer, x=(/LTNE_par(layer)%initcond/))
                pde_loc%solution(k) = value + nodes%data(k,D)
            end select
          end if
        end do   
      end do
      

    end subroutine wat_init
    
    subroutine temp_l_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals
      use geom_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (LTNE_par(1_ikind)%icondtype)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/LTNE/Tini_l.in", correct_h = .false.)
        case("value")
          do i=1, elements%kolik
            layer = elements%material(i)
            do j=1, ubound(elements%data,2)
              k = elements%data(i,j)
              l = nodes%edge(k)
              m = pde_loc%permut(k)
              if (m == 0) then
                call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
                pde_loc%solution(k) = value 
              else
                pde_loc%solution(k) = LTNE_par(layer)%Tinit_l
              end if
            end do   
          end do
      end select
      
    end subroutine temp_l_initcond
    
    subroutine temp_s_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals
      use geom_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (LTNE_par(1_ikind)%icondtypeTs)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/LTNE/Tini_s.in", correct_h = .false.)
        case("value")
          do i=1, elements%kolik
            layer = elements%material(i)
            do j=1, ubound(elements%data,2)
              k = elements%data(i,j)
              l = nodes%edge(k)
              m = pde_loc%permut(k)
              if (m == 0) then
                call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
                pde_loc%solution(k) = value 
              else
                pde_loc%solution(k) = LTNE_par(layer)%Tinit_s
              end if
            end do   
          end do
      end select
    
    
      if(.not.air) then
        if(allocated(T_air))then
        else
            allocate(T_air(nodes%kolik))
        end if
        do i=1, elements%kolik
          do j=1, ubound(elements%data,2)
            k = elements%data(i,j)
            T_air(k) = pde_loc%solution(k) 
    end do   
        end do
      end if
    end subroutine temp_s_initcond
    
    subroutine Kliquid_temp(pde_loc, layer, quadpnt, x, T, tensor, scalar) 
      use typy
      use global_objs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind),intent(in), optional    :: T
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      real(kind=rkind), dimension(3,3) :: Klh, Klt
      integer(kind=ikind) :: D
      real(kind = rkind) :: h_l
      D = drutes_config%dimen

      
      if (present(tensor)) then
        call mualem_ltne(pde_loc, layer, x=(/hl(pde_loc, layer, quadpnt)/), tensor = Klt(1:D, 1:D))
        if(qlt_log) then
          Klt(1:D, 1:D) = 10**(-Omega*Q_reduction(layer, quadpnt))*Klt(1:D, 1:D)
        else
          Klt(1:D,1:D)= 0_rkind*Klt(1:D, 1:D)
        end if 
        if (present(quadpnt)) then
          h_l = hl(pde_loc, layer, quadpnt)
          tensor = Klt(1:D, 1:D)*gwt*h_l*surf_tens_deriv(pde_loc, layer, quadpnt)/surf_tens_ref
        else
          print *, "runtime error"
          print *, "exited from ltne_helper:: Kliquid_temp"
          ERROR STOP
        end if
      else
         print *, "ERROR! output tensor undefined, exited from ltne_helper:: Kliquid_temp"
      end if    

      
    end subroutine Kliquid_temp
    
    function surf_tens_deriv(pde_loc, layer, quadpnt, T) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
    
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), intent(in), optional    ::  T
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
    
      real(kind=rkind) :: temp
    
      temp = pde(heat_proc)%getval(quadpnt)
      if (present(T)) then
        temp = T
      end if
      
      val = -0.1425-4.76e-4*temp
      
    end function surf_tens_deriv
    
    !> specific function for Richards equation in H-form (total hydraulic head form), replaces pde_objs::getvalp1 in order to distinguish between H and h 
    function getval_retotltne(pde_loc, quadpnt) result(val)
      use typy
      use pde_objs
      use geom_tools
      use re_globals
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind), dimension(3) :: xyz
      integer(kind=ikind) :: D, layer
      

      if (quadpnt%preproc) then
        D = drutes_config%dimen
             
        call getcoor(quadpnt, xyz(1:D))
        if (drutes_config%dimen>1) then
          val = getvalp1(pde_loc, quadpnt) - xyz(D)
        else
          layer = get_layer(quadpnt)
          val = getvalp1(pde(wat), quadpnt) - xyz(D)*cos(4*atan(1.0_rkind)/180*LTNE_par(layer)%anisoangle(1))
        end if
      else
        val = getvalp1(pde_loc, quadpnt)
      end if
    end function getval_retotltne
    
        
    function inverse_vangen_ltne(pde_loc, layer, quadpnt, x) result(hpress)
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
        print *, "exited from LTNE_helper::inverse_vangen_ltne"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from LTNE_helper::inverse_vangen_ltne"
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
          print *, "exited from LTNE_helper::inverse_vangen_ltne"
          ERROR STOP
        end if
        theta = x(1)
      end if
      
      
      
      a = LTNE_par(layer)%alpha
      n = LTNE_par(layer)%n
      m = LTNE_par(layer)%m
      
      if (abs(theta - LTNE_par(layer)%Ths) < epsilon(theta)) then
        hpress = 0
      else
        if (theta >  LTNE_par(layer)%Ths + 10*epsilon(theta)) then
          call write_log("theta is greater then theta_s, exiting")
          print *, "called from LTNE_helper::inverse_vangen_ltne"
          error stop
        else if (theta < 0) then
          call write_log("theta is negative strange, exiting")
          print *, "called from LTNE_helper::inverse_vangen_ltne"
          error stop 
        end if
        hpress = ((((LTNE_par(layer)%Ths - LTNE_par(layer)%Thr)/(theta-LTNE_par(layer)%Thr))**(1.0_rkind/m)-1) &  
        **(1.0_rkind/n))/(-a)
      end if
      
    end function inverse_vangen_ltne
    
    
    subroutine LTNE_coolant_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
     

      integer(kind=ikind) :: i, edge_id, j
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, T
      integer :: i1
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      
      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))        
        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              quadpnt%type_pnt = "ndpt"
              quadpnt%column=1
              quadpnt%order = elements%data(el_id,node_order)
              T =  pde_loc%getval(quadpnt)
              bcval = -hc*(T-pde_loc%bc(edge_id)%series(j,2))
              value = bcval
              EXIT
            end if
          end do
        else
          quadpnt%type_pnt = "ndpt"
          quadpnt%column=1
          quadpnt%order = elements%data(el_id,node_order)
          T =  pde_loc%getval(quadpnt)
          bcval = -hc*(T-pde_loc%bc(edge_id)%value)
          value = bcval
        end if
      end if
      if (present(code)) then
        code = 2
      end if


    end subroutine LTNE_coolant_bc
    
    
end module LTNE_helper
