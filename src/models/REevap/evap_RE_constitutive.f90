module evap_RE_constitutive
  public :: thetav, dthetav_dtemp, dthetav_dh
  private :: dens_satvap, dens_liquid, T2kelv, drelhumiddh, drelhumiddT, drho_svdT, invdrhol_dT, cond_vapour4h, vapour_diff
  
  contains
  
  
  
  
  
  
  
     !> Isothermal Properties of water vapor
    subroutine cond_vapour4h(layer, quadpnt, Kvh)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      
      !> MaterialID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpnt 
      !>unsaturated non-thermal conductuvity of water vapor
      real(kind=rkind), dimension(:,:), intent(out) :: Kvh 

      real(kind=rkind) :: val, tort, ths, theta, Da, D, T
      
      
      

      

      
!      rh_soil_val = rh_soil(layer, quadpnt)
!      rho_l_val= rho_l( quadpnt) 
!      rho_sv_val = rho_sv( quadpnt) 
!      diff = vapor_diff_soil(pde(heat_ord), layer, quadpnt)
      
!      T = pde(Heat_ord)%getval(quadpnt)
!      T = T + Tref
     
!      val = (diff/rho_l_val)*rho_sv_val*((MolWat*gravity)/(R_gas*T))*rh_soil_val
  
    end subroutine cond_vapour4h
    
    function thetav(pde_loc,layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use evapglob
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: theta, theta_s
      
      theta_s = vgset(layer)%ths
      
      theta = vangen(pde(re_ord), layer, quadpnt)
      
      val = (theta_s - theta)*relhumid(quadpnt) * dens_satvap(quadpnt) / dens_liquid(quadpnt)
    
    end function thetav
    
    !> derivative of the vapour water content with respect to temperature
    !! \dv{\theta_{v}}{T} = (\theta_{s} - \theta_{l})\left( \dv{H_{r}}{T} \frac{\rho_{sv}}{\rho_{l}} + H_{r}  \dv{\rho_{sv}}{T} \frac{1}{\rho_{l}} + H_{r} \rho_{sv} \dv{\frac{1}{\rho_{l}}}{T} \right) \f]
    !<
    function   dthetav_dtemp(pde_loc,layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use re_globals
      use evapglob
      use re_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta
      
      ths = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      
      val = (ths-theta)*(drelhumiddT(quadpnt)*dens_satvap(quadpnt)/dens_liquid(quadpnt) + &
            relhumid(quadpnt) * drho_svdT(quadpnt)/dens_liquid(quadpnt) + &
            relhumid(quadpnt)*dens_liquid(quadpnt)*invdrhol_dT(quadpnt))
    
    end function dthetav_dtemp
    
    
    function dthetav_dh(pde_loc,layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use re_globals
      use evapglob
      use re_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta, Cl
      
      ths = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      Cl = vangen_elast(pde(re_ord), layer, quadpnt)
      
      
      val = drelhumiddh(quadpnt)*ths*dens_satvap(quadpnt)/dens_liquid(quadpnt) - & 
              Cl*relhumid(quadpnt) * dens_satvap(quadpnt)/dens_liquid(quadpnt) - &
              drelhumiddh(quadpnt)*theta*dens_satvap(quadpnt)/dens_liquid(quadpnt) 
            
    
    end function dthetav_dh
    
    !> saturated vapour density
    !! \f[ \rho_{sv} = \num{e-3} \frac{\exp \left( 31.3716 - \frac{6014.79}{T} - \num{7.92495e-3}T \right)}{T} \f]
    !<
    function dens_satvap(quadpnt) result(val)
      use typy
      use evapglob
      use global_objs
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) :: T
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))
      
      val = 1e-3 *(exp(31.3716 - (6014.79/T) - 7.92495e-3*T))/T
      
    end function dens_satvap
    
    !> liquid density obtained from
    !! \f[ \rho_{l} =  1000 - \num{7.37e-3}(T - 3.98)^2 + \num{3.79e-5}(T - 3.98)^3 \f]
    !<
    function dens_liquid(quadpnt) result(val)
      use typy
      use evapglob
      use global_objs
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) :: T
      
      T = pde(heat_ord)%getval(quadpnt)
      
      val = 1000 - 7.37e-3*(T - 3.98)**2 + 3.79e-5*(T - 3.98)**3
    
    end function dens_liquid
    
    !> temperature in dg. C to kelvins
    function T2kelv(T) result(Tkelv)
      use typy
      real(kind=rkind), intent(in) :: T
      real(kind=rkind) :: Tkelv
      
      
      Tkelv = T + 273.15
    
    end function T2kelv
    
    !> relative humidity
    !! \f[ H_{r}= \left\{ \begin{array}{l l}\exp \left( \frac{hMg}{RT} \right) ,&  \mbox{if $h < 0$}\\ 1, & \mbox{ if $h \ge 0$}\\ \end{array} \right. \f]
    !<
    function relhumid(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) :: h, T
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))
      
      h = pde(re_ord)%getval(quadpnt)
      
      if (h >= 0) then
        val = 1
      else
        val = exp(h*MolWat*gravity/(R_gas*T))
      end if
      
      
    end function relhumid
    
    !> derivative of relative humidity with a respect to pressure head
    !! \f[ \dv{\theta_{v}}{h} = \dv{H_{r}}{h} \frac{\theta_{s} \rho_{s}}{\rho_{l}} -  C^{l}(h) \frac{H_{r}\rho_{sv}}{\rho_{l}} -  \dv{H_{r}}{h}  \frac{\theta_{l}\rho_{sv}}{\rho_{l}}. \f]
    !<
    function drelhumiddh(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) :: h, T
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))
      
      h = pde(re_ord)%getval(quadpnt)
      
      if (h < 0) then
        val = MolWat*gravity/(R_gas*T)*relhumid(quadpnt)
      else
        val = 0
      end if
    
    end function drelhumiddh
    
    
    !> derivative of relative humidity with a respect to temperature
    !! \f[ \dv{\theta_{v}}{h} = \dv{H_{r}}{h} \frac{\theta_{s} \rho_{s}}{\rho_{l}} -  C^{l}(h) \frac{H_{r}\rho_{sv}}{\rho_{l}} -  \dv{H_{r}}{h}  \frac{\theta_{l}\rho_{sv}}{\rho_{l}}. \f]
    !<
    function drelhumiddT(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) :: h, T
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))
      
      h = pde(re_ord)%getval(quadpnt)
      
      if (h < 0) then
        val = MolWat * gravity *h /(R_gas*T*T)*relhumid(quadpnt)
      else
        val = 0
      end if
    
    end function drelhumiddT
    
    !> Derivative of the saturated vapour density with a respect to temperature
    !! \f[ \dv{\rho_{sv}}{T} = -\dfrac{\left(317x^2+40000x-240591600\right)\exp\left(-\frac{317x}{40000}-\frac{601479}{100x}+\frac{78429}{2500}\right)}{40000000x^3} \f]
    !<
    function drho_svdT(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) ::  T   
    
      T = T2kelv(pde(heat_ord)%getval(quadpnt))    
    
      val = -((317.0*T*T+40000.0*T-240591600.0_rkind)*exp(-(317.0*T)/40000.0-601479.0/(100.0*T)+78429.0/2500.0))/(400000000.0*T*T*T)
      
    end function drho_svdT
    
    
    !>  derivative of the term \f[ \rho_{l}^{-1} \f] with a respect to temperature is given by
    !! \f[ 
    function invdrhol_dT(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) ::  T   
      
      T = pde(heat_ord)%getval(quadpnt)
      
      val = -(62500000000000.0*(50*T-199)*(56850.0*T-7596263.0))/(47375000.0*T**3-9778157500.0_rkind*T**2+75582816850.0_rkind*T &
              +1249851083567979.0_rkind)**2
              
    end function invdrhol_dT


    !> vapour soil diffusivity
    !! \f[ D = \tau(\theta_S-theta) D_a \f]
    !! \f[ \tau = \frac{(\theta_{s}-\theta_{l})^{\nicefrac{7}{3}}}{\theta_{s}^{2}} \f]
    !! \f[  D_{a} = \num{2.12e-5} \left(\frac{T}{273.15} \right)^{2}. \f]
    !<
    function vapour_diff(layer, quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use re_constitutive
      
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value: Vapor Difussivity in soil  [m^2/s]
      real(kind=rkind) :: val
      !> Water content,Volumetric air content,
      real(kind=rkind) :: theta_l, theta_air, Da, tort, T
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))    
  
      theta_air = vgset(layer)%Ths - vangen(pde(re_ord), layer, quadpnt)
      
      tort = ((theta_air)**(7.0/3.0))/ (vgset(layer)%Ths**2)
      
      Da = 2.12e-5*(T/273.15)**2
      
      val = tort*theta_air*Da

    end function vapour_diff
    
    


end module evap_RE_constitutive
