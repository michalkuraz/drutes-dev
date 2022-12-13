module freeze_helper
  use pde_objs
  use typy
  use freeze_globs
  use debug_tools
  use RE_constitutive

  public :: iceswitch,icefac, rho_icewat, Q_reduction, surf_tens_deriv, Kliquid_temp, hl, thetai, thetal
  public:: vangen_fr, mualem_fr, inverse_vangen_fr, temp_initcond, temp_s_initcond, wat_initcond, getval_retotfr
  public:: rho_wat, thetai_wat_eq, dhldT, hw_cl, theta_cl, thetai_eq, ice_initcond, thetas
  private:: linspace
      
  
  contains
    
     function T_fr(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> return value:Liquid water density  rho_l [kg/m^3]
      real(kind=rkind):: val

      
       if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified integ point "
        print *, "exited from freeze_helper::T_f"
        ERROR stop
      end if
    
      
      select case (freeze_par(layer)%material)
        case ("Soil")
          val = Tfk!*exp(hw_cl(pde(wat), layer, quadpnt)*grav/Lf)
        case ("Snow")
          val = Tfk
      end select

    end function T_fr
  
     function rho_wat(quadpnt) result(val)
      use typy
      use global_objs

      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value:Liquid water density  rho_l [kg/m^3]
      real(kind=rkind):: val
      !>  Temperature T in ºC 
      real(kind=rkind):: Temp_C
      
       if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified integ point "
        print *, "exited from freeze_helper::rho_l"
        ERROR stop
      end if
    
      Temp_C = pde(heat_proc)%getval(quadpnt)-273.15
      val = 1000.0_rkind - 7.37e-3*(Temp_C - 4.0_rkind)**2 + 3.79e-5*(Temp_C -4.0_rkind)**3

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
      
      Tf = Tfk
   
      

      if (pde(heat_proc)%getval(quadpnt_loc) > Tf) then
      !> melting
        sw = .FALSE.
      else
      !> freezing
        sw = .TRUE.
      end if
          
    end function iceswitch

    function icefac(quadpnt, Tf) result(fac)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind), intent(in) :: Tf
      real(kind=rkind) :: fac, pi, sin_T, sin_out
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.
      
 !     Tf = Tfk
      

      if (pde(heat_proc)%getval(quadpnt_loc) > Tf) then
      !> melting
        fac = 0
      else
        if(pde(heat_proc)%getval(quadpnt_loc) < Tr) then
          fac = 1 
        else
        !> freezing sign function
          pi = 4*atan(1.0_rkind)
          sin_T = (Tr - pde(heat_proc)%getval(quadpnt_loc))/(Tr-Tf)
          sin_out = sin(sin_T*pi+pi/2)
          fac = sin_out/2_rkind+0.5_rkind
        end if
      end if
    end function icefac

    function gaussianint(end, start, Tf) result(val)
      use typy
      use global_objs
      
      real(kind=rkind), intent(in) :: end, start, Tf
      real(kind=rkind), dimension(3) :: a, H, hh, aa, w
      real(kind=rkind) :: val
      integer(kind = ikind) :: i

      a = (/-(5.0_rkind/9.0_rkind)**0.5, 0.0_rkind, (3.0/5.0_rkind)**0.5/)
      H = (/5.0/9.0_rkind, 8_rkind/9.0_rkind, 5_rkind/9.0_rkind/)
      hh = (end-start)/2.0_rkind*H
      aa = (end - start)/2.0_rkind*a+(end + start)/2.0_rkind
      do  i = 1, 3
        w(i) = dhldT(T = aa(i), Tf = Tf)
      end do
      val = sum(hh*w)

    end function gaussianint
    
    function dhldT(T, Tf) result(val)
      use typy
      use global_objs
      
      real(kind=rkind), intent(in) :: T, Tf
      real(kind=rkind) :: fac, val    
      real(kind=rkind) :: pi, sin_T, sin_out
      
      !Tf = Tfk
      

      if (T > Tf) then
      !> melting
        fac = 0
      else
        if(T < Tr) then
            fac = 1
        else
        !> freezing sigmoid function
            pi = 4*atan(1.0_rkind)
            sin_T = (Tr - T)/(Tr-Tf)
            sin_out = sin(sin_T*pi+pi/2)
            fac = sin_out/2_rkind+0.5_rkind
        end if
      end if
      val = fac*Lf/grav/T
          
    end function dhldT
    
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
      
      thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      thall = vangen_fr(pde(wat), layer, quadpnt)
      thice = thall - thl
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
      select case (drutes_config%name)
        case ("freeze", "LTNE")
          thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
       case("ICENE")
          thl = vangen_fr(pde(wat), layer, quadpnt)
       end select 
      thice = thetai(pde(wat), layer, quadpnt)
      !val = thice/(thice+thl*1.09)
      val = thice/(thice+thl-freeze_par(layer)%Thr)      
       
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
        select case (drutes_config%name)
          case ("freeze", "LTNE")
            call mualem_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/), tensor = Klt(1:D, 1:D))
          case("ICENE")
            call mualem_fr(pde(wat), layer, quadpnt, tensor = Klt(1:D, 1:D))
        end select 
        if(qlt_log) then
          Klt(1:D, 1:D) = Klt(1:D, 1:D)!*10**(-Omega*Q_reduction(layer, quadpnt))
        else
          Klt(1:D,1:D)= 0_rkind*Klt(1:D, 1:D)
        end if 
        if (present(quadpnt)) then
          h_l = hl(pde(wat), layer, quadpnt)
          tensor = Klt(1:D, 1:D)*gwt*h_l*surf_tens_deriv(pde(heat_proc), layer, quadpnt)/surf_tens_ref
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
    
      real(kind=rkind) :: temp_C
    
      temp_C = pde(heat_proc)%getval(quadpnt)-273.15
      if (present(T)) then
        temp_C = T
      end if
      
      val = -0.1425-4.76e-4*temp_C
      
    end function surf_tens_deriv
    
    subroutine linspace(from, to, array)
        real(kind=rkind), intent(in) :: from, to
        real(kind=rkind), intent(out) :: array(:)
        real(kind=rkind) :: range
        integer :: n, i
        n = size(array)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if
        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end subroutine
    
    function hl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val, T_f, fac, dif, T1K, T2K, T_threshK
      real(kind=rkind) :: hw, temp, tempK, midtemp,  Tstart, diffx
      real(kind=rkind) :: integ, integ2, integ3, integ4
      real(kind=rkind) :: sin_T, sin_out, pi
      real(kind=rkind) :: T_threshK99, meanKs
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind), dimension(:), allocatable :: intpoints
      integer :: n, i
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.
      
      select case (drutes_config%name)
        case ("ICENE")
          val = pde(wat)%getval(quadpnt_loc)
        case ("freeze", "LTNE")
        hw = pde(wat)%getval(quadpnt_loc)
        if(hw > 0) then
          hw = 0
        end if
        temp = pde(heat_proc)%getval(quadpnt)
        T_f = T_fr(pde(wat), layer, quadpnt)
        if(temp < T_f) then
            if(temp  < Tr) then
              fac = 1 
            else
            !> freezing sign function
              pi = 4*atan(1.0_rkind)
              sin_T = (Tr -temp)/(Tr-T_f)
              sin_out = sin(sin_T*pi+pi/2)
              fac = sin_out/2_rkind+0.5_rkind
            end if
            tempK = temp
            T_threshK99 = Tr
            if(fac > 0.99_rkind) then
                dif = T_f-T_threshK99
                Tstart = T_threshK99
            else
                dif = T_f-TempK   
                Tstart = TempK
            end if
            meanKs = sum(freeze_par(layer)%Ks)/max(1,size(freeze_par(layer)%Ks))*8.64e+6
            if(meanKS < 10) then
                meanKS = 10
            end if
            if(meanKS > 25) then
                meanKs = 25
            end if
            diffx = (T_f-T_threshK99)/(meanKs*0.75)
            if((T_f-T_threshK99) < 1.0) then
                diffx = 0.1
            end if
            n = nint(dif/diffx)+1
            allocate(intpoints(n))
            call linspace(from=Tstart, to=T_f, array=intpoints)
            val = hw
            do i=1,n-1
              val = val + gaussianint(start = intpoints(i+1), end = intpoints(i), Tf = T_f)
            end do
            if(fac > 0.99_rkind) then
                val = val + Lf/grav*log(tempK/T_threshK99)
            end if
        else
            val = hw
        end if
        
        if(isnan(val)) then
            print*, "hw is not a number! from freeze_helper::hl"
            print*, "hw", hw
            print*, "temp", temp
            print*, pde(wat)%getval(quadpnt)
        end if
     end select
    end function hl
    
    
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
      select case (drutes_config%name)
        case ("freeze", "LTNE")
            thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
            thall = vangen_fr(pde(wat), layer, quadpnt)
            val = thall - thl
        case("ICENE")
          val = pde(ice)%getval(quadpnt)
      end select
      
      select case (freeze_par(layer)%material)
          case ("Soil")
            continue
          case ("Snow")
            continue
      end select
        !val = rho_wat(quadpnt)/rho_ice*val

    end function thetai
    
    
    function thetai_eq(pde_loc, layer, quadpnt, x) result(val)
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
      real(kind=rkind) :: T, theta_tot, theta_ice, Tf, minice
      real(kind=rkind) :: theta_l,theta_f, cp, th_air

      if(.not. present(quadpnt)) then
        print*, x
        print*, "Quadpnt needed"
        print *, "exited from freeze_helper::thetai"
        ERROR STOP
      end if
      

      
      select case (freeze_par(layer)%material)
        case ("Soil")
          Tf = T_fr(pde(wat), layer, quadpnt)
          T = pde(heat_proc)%getval(quadpnt)
          theta_ice = pde(ice)%getval(quadpnt)
          if(T .LE. Tf) then
            theta_ice = pde(ice)%getval(quadpnt)
            theta_tot = vangen_fr(pde_loc, layer, quadpnt) + theta_ice
            if(theta_tot > thetas(pde_loc, layer, quadpnt)) then
              theta_tot = thetas(pde_loc, layer, quadpnt)
            end if
            val = theta_tot - theta_cl(pde_loc, layer, quadpnt) 
          else 
            val = 0
          end if        
        case ("Snow")
          Tf = T_fr(pde(wat), layer, quadpnt)
          T = pde(heat_proc)%getval(quadpnt)
          theta_ice = pde(ice)%getval(quadpnt)
          theta_l = vangen_fr(pde(wat), layer, quadpnt)
          th_air = thetas(pde_loc, layer, quadpnt)-theta_l
          cp = freeze_par(layer)%Ca*th_air*rho_air 
          cp = cp+ rho_ice*freeze_par(layer)%Ci*theta_ice
          cp = cp+theta_l*freeze_par(layer)%Cl*rho_wat(quadpnt)
          if(T .LE. Tf) then
            !freezing
            theta_f = min(theta_l, cp*(Tf-T)/Lf/rho_wat(quadpnt))
            val = theta_ice + theta_f
          else 
          !melting
            minice = 0_rkind
            val = theta_ice - cp*(T-Tf)/Lf/rho_ice
            val = max(minice, val)      
          end if
       end select
    end function thetai_eq
    
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
       thi = thetai(pde(wat), layer, quadpnt)
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

      val = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
    end function thetal

    function hw_cl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      real(kind=rkind) :: theta_l
      real(kind=rkind) :: theta_ice
      real(kind=rkind) :: T
      real(kind=rkind) :: theta_tot, hw, hcl
      real(kind = rkind) :: ths
      real(kind = rkind) :: minice
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%preproc=.true.
      select case (drutes_config%name)
        case ("freeze", "LTNE")
          val = pde(wat)%getval(quadpnt_loc)
        case("ICENE")
		  minice = 0_rkind
		  theta_l = vangen_fr(pde(wat), layer, quadpnt)
		  
		 !if(T < T_fr(pde(wat), layer, quadpnt)) then
			theta_ice = max(pde(ice)%getval(quadpnt),minice)
			theta_tot = theta_l + theta_ice
			ths = thetas(pde(wat), layer, quadpnt)
			if(theta_tot > ths)then
			  theta_tot = ths
			end if
			hw = inverse_vangen_fr(pde(wat), layer, x = (/theta_tot/))
			val =  hw
		!  else 
		!	val = pde(wat)%getval(quadpnt_loc)
		!  end if
      end select 
    end function hw_cl
    
    function theta_cl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      real(kind=rkind) :: T, T_f
      real(kind=rkind) :: hcl
      
      
      select case (freeze_par(layer)%material)
        case ("Soil")
		  T = pde(heat_proc)%getval(quadpnt)
		  T_f = T_fr(pde(wat), layer, quadpnt)
		  if(T < T_f) then
			hcl = hw_cl(pde(wat), layer, quadpnt) + Lf/grav*log(T/T_f) 
			val = vangen_fr(pde(wat), layer, x = (/hcl/))
		  else 
			val = vangen_fr(pde(wat), layer, quadpnt)
		  end if
         case ("Snow")
        val = vangen_fr(pde(wat), layer, quadpnt)
      end select

    end function theta_cl
    
    
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
            select case (drutes_config%name)
              case ("freeze", "LTNE")
               val = freeze_par(layer)%ths
              case("ICENE")
               val = 1-pde(ice)%getval(quadpnt)
               if(val < freeze_par(layer)%Thr) then
                 print*, val
                 val = freeze_par(layer)%Thr
               end if
           end select
      end select
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
          ths = thetas(pde(wat), layer, quadpnt)
        else 
          ths = thetas(pde(wat), layer, x = x)
        end if
      else
        if (present(quadpnt)) then
          quadpnt_loc=quadpnt
          quadpnt_loc%preproc=.true.
          h = pde(wat)%getval(quadpnt_loc)
          ths = thetas(pde(wat), layer, quadpnt)
        else
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_helper::vangen_fr"
          ERROR STOP
        end if
      end if
      

      
      a = freeze_par(layer)%alpha
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
          ts = thetas(pde(wat), layer, quadpnt)
        else 
          ts = thetas(pde(wat), layer, x = x)
        end if
      else
        if (present(quadpnt)) then
          quadpnt_loc=quadpnt
          quadpnt_loc%preproc=.true.
          h = pde(wat)%getval(quadpnt_loc)
          ts = thetas(pde(wat), layer, quadpnt)
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
        h = pde(wat)%getval(quadpnt_loc)
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
                pde_loc%solution(k) = freeze_par(layer)%initcond + &
                nodes%data(k,D)*cos(4*atan(1.0_rkind)/180*freeze_par(layer)%anisoangle(1))
              case("theta")
                value = inverse_vangen_fr(pde_loc, layer, x=(/freeze_par(layer)%initcond/))
                pde_loc%solution(k) = value + nodes%data(k,D)*cos(4*atan(1.0_rkind)/180*freeze_par(layer)%anisoangle(1))
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
          call map1d2dJ(pde_loc,"drutes.conf/freeze/iceini.in", correct_h = .false.)
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
    
    
    subroutine freeze_coolant_bc(pde_loc, el_id, node_order, value, code, array, nvectin) 
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
      real(kind=rkind), dimension(:), intent(in), optional :: nvectin
     

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
    
    
    subroutine freeze_coolant_bc_bot(pde_loc, el_id, node_order, value, code, array, nvectin) 
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
      real(kind=rkind), dimension(:), intent(in), optional :: nvectin   

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
              bcval = -hcbot*(T-pde_loc%bc(edge_id)%series(j,2))
              value = bcval
              EXIT
            end if
          end do
        else
          quadpnt%type_pnt = "ndpt"
          quadpnt%column=1
          quadpnt%order = elements%data(el_id,node_order)
          T =  pde_loc%getval(quadpnt)
          bcval = -hcbot*(T-pde_loc%bc(edge_id)%value)
          value = bcval
        end if
      end if
      if (present(code)) then
        code = 2
      end if

    end subroutine freeze_coolant_bc_bot
    
end module freeze_helper
