module freeze_helper
  use pde_objs
  use typy
  use freeze_globs
  use debug_tools
  use RE_constitutive

  public :: iceswitch,icefac, rho_icewat, Q_reduction, surf_tens_deriv, Kliquid_temp, hl, hcl, thetai, thetal, vf, mf
  public:: vangen_fr, mualem_fr, temp_initcond, temp_s_initcond, wat_initcond, getval_retotfr, ice_initcond, thetas
  public:: phase_ice, phase_wat, latent_heat_vf, latent_heat_mf
  public:: rho_wat, thetai_wat_eq
    
      
  
  contains
     function rho_wat(quadpnt) result(val)
      use typy
      use global_objs

      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value:Liquid water density  rho_l [kg/m^3]
      real(kind=rkind):: val
      !>  Temperature T in ºC 
      real(kind=rkind):: T
      
       if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified integ point "
        print *, "exited from evap_auxfnc::rho_l"
        ERROR stop
      end if
    
      T = pde(heat_proc)%getval(quadpnt)
      val = 1000.0_rkind - 7.37e-3*(T - 4.0_rkind)**2 + 3.79e-5*(T -4.0_rkind)**3

    end function rho_wat

      !> switch for freezing condition based on Clapeyron equation 
    function iceswitch(quadpnt) result(sw)
      use typy
      use global_objs
      
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

    function icefac(quadpnt) result(fac)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: Tf, fac, x
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
        fac = 0
      else
      !> freezing sigmoid function
        x = pde(heat_proc)%getval(quadpnt_loc)-Tf 
        fac = 1_rkind/(1_rkind+exp(x*fac_scale + fac_add))
      end if
          
    end function icefac
    
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
      
      if(fr_rate) then
        thl = vangen_fr(pde(wat), layer,  quadpnt)
        thice = thetai(pde(wat), layer, quadpnt)
        thall = thl+thice
      else
        thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
        thall = vangen_fr(pde(wat), layer, quadpnt)
        thice = thall - thl
      end if
      rho = (thl * rho_wat(quadpnt) + thice * rho_ice)/thall
       
    end function rho_icewat
    
    function Q_reduction(layer, quadpnt, x) result(val)

      use typy
      use global_objs
      use re_globals
      
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      integer(kind=ikind) :: el
      real(kind=rkind) :: thl, thice, val
      
      if(.not. present(quadpnt)) then
        print*, "Cant estimate Q_reduction without quadpnt"
        print *, "exited from freeze_helper::Q_reduction"
        ERROR STOP
      end if
      if(fr_rate) then
        thl = vangen_fr(pde(wat), layer, quadpnt)
      else
        thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      end if
      thice = thetai(pde(wat), layer, quadpnt)
      val = thice/(thice+thl*1.09)
      
       
    end function Q_reduction
    
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
        if(fr_rate) then 
          call mualem_fr(pde_loc, layer, quadpnt, tensor = Klt(1:D, 1:D))
        else
          call mualem_fr(pde_loc, layer, x=(/hl(pde_loc, layer, quadpnt)/), tensor = Klt(1:D, 1:D))
        end if
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
          print *, "exited from Kliquid_temp::freeze_helper"
          ERROR STOP
        end if
      else
         print *, "ERROR! output tensor undefined, exited from Kliquid_temp::freeze_helper"
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
    
    function hl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val, T_f, fac
      real(kind=rkind) :: hw, temp
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.

      hw = pde(wat)%getval(quadpnt_loc)
      temp = pde(heat_proc)%getval(quadpnt)
      T_f = Tref*exp(hw*grav/Lf)
      
      if(iceswitch(quadpnt)) then
        if(fr_rate) then
          val = hw
        else
          fac = icefac(quadpnt)
          val = hw+fac*Lf/grav*log((temp+273.15_rkind)/T_f) !units       
        end if        
      else
        val = hw
      end if
      
      if(isnan(val)) then
        print*, "hw is not a number! from freeze_helper::hl"
        print*, "hw", hw
        print*, "temp", temp
      end if
    end function hl
    
    function hcl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val, T_f
      real(kind=rkind) :: hw, temp
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.

      hw = pde(wat)%getval(quadpnt_loc)
      
      temp = pde(heat_proc)%getval(quadpnt)
      T_f = Tref*exp(hw*grav/Lf)
      
      if(iceswitch(quadpnt)) then
        val = hw+Lf/grav*log((temp+273.15_rkind)/T_f) !units       
      else
        val = hw
      end if
      
    end function hcl
    
    function thetai(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use debug_tools
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind) :: thl, thall, thi_tmp
      
      if(.not. present(quadpnt)) then
        print*, x
        print*, "Quadpnt needed"
        print *, "exited from freeze_helper::thetai"
        ERROR STOP
      end if
      
      if(fr_rate) then
        val =  pde(ice)%getval(quadpnt)
      else
        thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
        thall = vangen_fr(pde(wat), layer, quadpnt)
        val = thall - thl
        !val = rho_wat(quadpnt)/rho_ice*val
      end if

    end function thetai
    

    function thetai_wat_eq(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use debug_tools
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind) :: thi
      
      if(.not. present(quadpnt)) then
        print*, x
        print*, "Quadpnt needed"
        print *, "exited from freeze_helper::thetai"
        ERROR STOP
      end if
       thi = thetai(pde_loc, layer, quadpnt)
       val = thi*rho_ice/rho_wat(quadpnt)

    end function thetai_wat_eq
    
    function thetal(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      if(fr_rate) then
        val = vangen_fr(pde(wat), layer, quadpnt)
      else
        val = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      end if
    end function thetal
    
    
    function thetas(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
    
    
      select case (freeze_par(layer)%material)
        case ("Soil")
          val = freeze_par(layer)%ths
        case ("Snow")
          if(fr_rate) then
            if (present(x)) then
              print*, "EROOR thetas"
              ERROR STOP
            else
              val = 1-thetai(pde_loc, layer, quadpnt)
            end if
          else 
            val = freeze_par(layer)%ths
          end if
      end select
      
      if(val > 1.0) then
        val = 1.0
      end if
    end function thetas
    
    
    
        !> \brief Van Genuchten relation \f[ \theta = f(pressure) \f]
    !!  \f[ \theta_e = \frac{1}{(1+(\alpha*h)^n)^m} \f]
    !! water content is considered as absolute value not the relative one \n
    !! see \f[ \theta_e = \frac{\theta - \theta_r}{\theta_s-\theta_r} \f]
    !<
    function vangen_fr(pde_loc, layer, quadpnt, x) result(theta)
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

      real(kind=rkind) :: a,n,m, theta_e, ths
      type(integpnt_str) :: quadpnt_loc
      

!       if (present(quadpnt) .and. present(x)) then
!         print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
!         print *, "exited from freeze_helper::vangen_fr"
!         ERROR stop
!       else if (.not. present(quadpnt) .and. .not. present(x)) then
!         print *, "ERROR: you have not specified either integ point or x value"
!         print *, "exited from freeze_helper::vangen_fr"
!         ERROR stop
!       end if
      if (present(x)) then
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_helper::vangen_fr"
          ERROR STOP
        end if
        h = x(1)
        if (present(quadpnt)) then
          ths = thetas(pde_loc, layer, quadpnt)
        else 
          ths = thetas(pde_loc, layer, x = x)
        end if
      else
        if (present(quadpnt)) then
          quadpnt_loc=quadpnt
          quadpnt_loc%preproc=.true.
          h = pde_loc%getval(quadpnt_loc)
          ths = thetas(pde_loc, layer, quadpnt)
        else
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_helper::vangen_fr"
          ERROR STOP
        end if
      end if
      

      
      a = freeze_par(layer)%alpha
      n = freeze_par(layer)%n
      m = freeze_par(layer)%m
      
      if (h >=0.0_rkind) then
        theta = ths
        RETURN
      else
        theta_e = 1/(1+(a*abs(h))**n)**m
        theta = theta_e*(ths-freeze_par(layer)%Thr)+freeze_par(layer)%Thr
      end if

    end function vangen_fr
    
    
    
        !> \brief so-called retention water capacity, it is a derivative to retention curve function
    !! \f E(h) = C(h) + \frac{\theta(h)}{\theta_s}S_s \f]
    !! where
    !! \f[ C(h) = \left\{ \begin{array}{l l}\frac{m n \alpha  (-h \alpha )^{-1+n}}{\left(1+(-h \alpha )^n\right)^{1+m}}(\theta_s - \theta_r) ,  & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ 0, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f]
    !! and 
    !! \f[ \theta(h) = \left\{ \begin{array}{l l} \frac{\theta_s -\theta_r}{(1+(-\alpha h)^n_{vg})^m_{vg}} + \theta_r,  & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ \theta_S, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f]
    !<
    function vangen_elast_fr(pde_loc,layer, quadpnt, x) result(E)
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
      real(kind=rkind) :: trsh = 0   
          
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from freeze_helper::vangen_elast_fr"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from freeze_helper::vangen_elast_fr"
        ERROR stop
      end if
      
      if (present(x)) then
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_helper::vangen_elast_fr"
          ERROR STOP
        end if
        h = x(1)
        if (present(quadpnt)) then
          ts = thetas(pde_loc, layer, quadpnt)
        else 
          ts = thetas(pde_loc, layer, x = x)
        end if
      else
        if (present(quadpnt)) then
          quadpnt_loc=quadpnt
          quadpnt_loc%preproc=.true.
          h = pde_loc%getval(quadpnt_loc)
          ts = thetas(pde_loc, layer, quadpnt)
        else
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_helper::vangen_elast_fr"
          ERROR STOP
        end if
      end if

      if (h < 0) then
        a = freeze_par(layer)%alpha
        n = freeze_par(layer)%n
        m = freeze_par(layer)%m
        tr = freeze_par(layer)%Thr

        C = a*m*n*(-tr + ts)*(-(a*h))**(-1 + n)*(1 + (-(a*h))**n)**(-1 - m)
      else
        C = 0
      end if

      E = max(C, trsh)
      

    end function vangen_elast_fr
    
    
    
    
    !> \brief Mualem's fucntion for unsaturated hydraulic conductivity with van Genuchten's water content substitution
    !! \f[   K(h) = \left\{ \begin{array}{l l} K_s\frac{\left( 1- (-\alpha h)^{n_{vg}m_{vg}} \left( 1+ (-\alpha h)^{n_{vg}} \right)^{-m_{vg}} \right)^2}{\left(1+(-\alpha h)^{n_{vg}} \right)^{\frac{m_{vg}}{2}}},  &  \mbox{$\forall$  $h \in$ $(-\infty,0)$}\\ K_s,  \mbox{$\forall$   $h \in$ $\langle 0, +\infty)$}\\ \end{array} \right. \f]
    !<
    subroutine mualem_fr(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use freeze_globs
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
        a = freeze_par(layer)%alpha
        n = freeze_par(layer)%n
        m = freeze_par(layer)%m

        tmp =  (1 - (-(a*h))**(m*n)/(1 + (-(a*h))**n)**m)**2/(1 + (-(a*h))**n)**(m/2.0_rkind)
      end if
	
      if (present(tensor)) then
        tensor = tmp* freeze_par(layer)%Ks
      end if

      if (present(scalar)) then
        scalar = tmp
      end if
    end subroutine mualem_fr
    
    subroutine wat_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use geom_tools
      use debug_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (freeze_par(1_ikind)%icondtypeRE)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/freeze/hini.in", correct_h = .true.)
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
            select case (freeze_par(layer)%icondtypeRE)
              case("H_tot")
                pde_loc%solution(k) = freeze_par(layer)%initcond !+ nodes%data(k,1)
              case("hpres")
                pde_loc%solution(k) = freeze_par(layer)%initcond + nodes%data(k,D)
              case("theta")
                value = inverse_vangen_fr(pde_loc, layer, x=(/freeze_par(layer)%initcond/))
                pde_loc%solution(k) = value + nodes%data(k,D)
            end select
          end if
        end do   
      end do
      

    end subroutine wat_initcond
    
    subroutine ice_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use geom_tools
      use debug_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (freeze_par(1_ikind)%icondtypeIce)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/freeze/iceini.in", correct_h = .true.)
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
            select case (freeze_par(layer)%icondtypeIce)
              case("theta")
                pde_loc%solution(k) = freeze_par(layer)%iceini 
            end select
          end if
        end do   
      end do

    end subroutine ice_initcond
    
    subroutine temp_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals
      use geom_tools
      use debug_tools
      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (freeze_par(1)%icondtype)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/freeze/Tini.in", correct_h = .false.)
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
                pde_loc%solution(k) = freeze_par(layer)%Tinit
              end if
            end do   
          end do
      end select
    end subroutine temp_initcond
    
     subroutine temp_s_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals
      use geom_tools
      use debug_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
      D = drutes_config%dimen
      select case (freeze_par(1_ikind)%icondtypeTs)
        case("input")
          call map1d2dJ(pde_loc,"drutes.conf/freeze/Tini_s.in", correct_h = .false.)
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
                pde_loc%solution(k) = freeze_par(layer)%Tinit_s
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

    
    !> specific function for Richards equation in H-form (total hydraulic head form), replaces pde_objs::getvalp1 in order to distinguish between H and h 
    function getval_retotfr(pde_loc, quadpnt) result(val)
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
          val = getvalp1(pde_loc, quadpnt) - xyz(D)*cos(4*atan(1.0_rkind)/180*freeze_par(layer)%anisoangle(1))
        end if
        
        
      else
        val = getvalp1(pde_loc, quadpnt)
      end if
	
      
    end function getval_retotfr
    
        
    function inverse_vangen_fr(pde_loc, layer, quadpnt, x) result(hpress)
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
      
      
      real(kind=rkind) :: a,n,m, ths
      type(integpnt_str) :: quadpnt_loc
      
 

      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from freeze_helper::inverse_vangen_fr"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from freeze_helper::inverse_vangen_fr"
        ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
        quadpnt_loc%preproc=.true.
        theta = pde_loc%getval(quadpnt_loc)
        ths = thetas(pde_loc, layer, quadpnt)
      else
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_helper::inverse_vangen_fr"
          ERROR STOP
        end if
        theta = x(1)
        ths = freeze_par(layer)%ths
      end if
      
      
      
      a = freeze_par(layer)%alpha
      n = freeze_par(layer)%n
      m = freeze_par(layer)%m
      

      
      if (abs(theta - ths) < epsilon(theta)) then
        hpress = 0
      else
        if (theta >  ths + 10*epsilon(theta)) then
          call write_log("theta is greater then theta_s, exiting")
          print *, "called from freeze_helper::inverse_vangen_fr"
          error stop
        else if (theta < 0) then
          call write_log("theta is negative strange, exiting")
          print *, "called from freeze_helper::inverse_vangen_fr"
          error stop 
        end if
        hpress = ((((ths - freeze_par(layer)%Thr)/(theta-freeze_par(layer)%Thr))**(1.0_rkind/m)-1) &  
        **(1.0_rkind/n))/(-a)
      end if
      
    end function inverse_vangen_fr
    
    
    subroutine freeze_coolant_bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
     

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

    end subroutine freeze_coolant_bc
    
  !> Functions below are ONLY relevant for freezing rate  
    function vf(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use globals
      use global_objs
      use pde_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      real(kind=rkind) :: thcl, thl, T, T_f, hw, super_cooling
      real(kind=rkind) :: T_tomelt, T_dif, melt, ths
      real(kind=rkind) :: a,n,m, theta_e
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.
      ! Freezing rate
      hw = pde(wat)%getval(quadpnt_loc)
      T = pde(heat_proc)%getval(quadpnt)
      T_f = Tref*exp(hw*grav/Lf)-273.15
      super_cooling = T_f-T
      thl = vangen_fr(pde(wat), layer, quadpnt)
      ths = thetas(pde(wat), layer, quadpnt)
      a = freeze_par(layer)%alpha
      n = freeze_par(layer)%n
      m = freeze_par(layer)%m
      
      if (hcl(pde(wat), layer, quadpnt) >=0.0_rkind) then
        thcl = ths
        RETURN
      else
        theta_e = 1/(1+(a*(abs(hcl(pde(wat), layer, quadpnt))))**n)**m
        thcl = theta_e*(ths-freeze_par(layer)%Thr)+freeze_par(layer)%Thr
      end if

      thl = vangen_fr(pde(wat), layer, quadpnt)
      if(super_cooling < epsilon(super_cooling)) then
        val = 0
      else 
        val = beta*(thl-thcl)*(T_f-T)**(0.33333333_rkind)
      end if

    end function vf
        
    function phase_ice(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use debug_tools
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind) :: thl, thall, thi_tmp, T_tomelt,T
      select case (drutes_config%name)
        case ("LTNE")
          T_tomelt = pde(heat_solid)%getval(quadpnt)
        case("freeze")
          T = pde(heat_proc)%getval(quadpnt)
          T_tomelt = T
      end select
      val =  vf(pde_loc, layer, quadpnt)-&
       mf_react(pde_loc, layer, quadpnt)*T_tomelt+&
       mf(pde_loc, layer, quadpnt)

    end function phase_ice
    
    function mf(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use globals
      use global_objs
      use pde_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      real(kind=rkind) :: thcl, thl, T, T_f, hw, super_cooling
      real(kind=rkind) :: T_tomelt, T_dif, melt, ths
      real(kind=rkind) :: a,n,m, theta_e
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.
      ! Freezing rate
      hw = pde(wat)%getval(quadpnt_loc)
      T = pde(heat_proc)%getval(quadpnt)
      T_f = Tref*exp(hw*grav/Lf)-273.15
      ! melting
     ! select case (drutes_config%name)
     !   case ("LTNE")
     !     T_tomelt = pde(heat_solid)%getval(quadpnt)
     !   case("freeze")
      !    T_tomelt = T
      !end select
      
      if(T_tomelt >= T_f) then 
        !T_dif = T_tomelt-T_f
        T_dif = -T_f
        val =  freeze_par(layer)%Ci*T_dif/Lf/beta_melt
      else
        val = 0
      end if

    end function mf
    
    function mf_react(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use globals
      use global_objs
      use pde_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      real(kind=rkind) :: thcl, thl, T, T_f, hw, super_cooling
      real(kind=rkind) :: T_tomelt, T_dif, melt, ths
      real(kind=rkind) :: a,n,m, theta_e
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.

      ! melting
      select case (drutes_config%name)
        case ("LTNE")
          T_tomelt = pde(heat_solid)%getval(quadpnt)
        case("freeze")
          T = pde(heat_proc)%getval(quadpnt)
          T_tomelt = T
      end select
      
      if(T_tomelt >= T_f) then 
        T_dif = T_tomelt-T_f
        val =  freeze_par(layer)%Ci/Lf/beta_melt
      else
        val = 0
      end if
    end function mf_react
    
    function phase_wat(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use debug_tools
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind) :: thl, thall, thi_tmp, T, T_tomelt
      select case (drutes_config%name)
        case ("LTNE")
          T_tomelt = pde(heat_solid)%getval(quadpnt)
        case("freeze")
          T = pde(heat_proc)%getval(quadpnt)
          T_tomelt = T
      end select
      val =   mf_react(pde_loc, layer, quadpnt)*T_tomelt&
      -mf(pde_loc, layer, quadpnt)-vf(pde_loc, layer, quadpnt)
    end function phase_wat
    
    function latent_heat_vf(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use globals
      use global_objs
      use pde_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val

      val = Lf*rho_wat(quadpnt)*vf(pde_loc, layer, quadpnt)

    end function latent_heat_vf
    
    function latent_heat_mf(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use globals
      use global_objs
      use pde_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      val = mf_react(pde_loc, layer, quadpnt)-mf(pde_loc, layer, quadpnt)
      val = Lf*rho_wat(quadpnt)*val

    end function latent_heat_mf
    
    function latent_heat(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use globals
      use global_objs
      use pde_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      
      val = mf_react(pde_loc, layer, quadpnt)-mf(pde_loc, layer, quadpnt)
      val = Lf*rho_wat(quadpnt)*(vf(pde_loc, layer, quadpnt)-val)
    end function latent_heat
    
end module freeze_helper
