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


module evap_RE_constitutive
  public :: REdiffhh, REdiffhT, REcapacityhh, REcapacityhT, darcy4liq, darcy4vap, totalflux, cond_vapour4h
  private ::  drelhumiddh, drelhumiddT, drhosv_dT, invdrhol_dT, vapour_diff, cond_ht
              
  
  public :: dens_liquid, dens_satvap, relhumid, thetav, dthetav_dtemp, dthetav_dh, T2kelv
  
  contains
  
  
    subroutine REdiffhh(pde_loc, layer, quadpntin, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use re_constitutive
      use evapglob
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      
      type(integpnt_str) :: quadpnt
      
      integer(kind=ikind) :: D
      real(kind=rkind), dimension(3,3) :: Klh, Kvh
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_RE_consitutive::REdiffhh"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      if (present(scalar)) then
        print *, "runtime error evap_RE_consitutive::REdiffhh"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      if (present(tensor)) then
        call mualem(pde(re_ord), layer, quadpnt, tensor=Klh(1:D, 1:D))
        call cond_vapour4h(layer, quadpnt, Kvh(1:D, 1:D))
        tensor(1:D, 1:D) = Klh(1:D, 1:D) + Kvh(1:D, 1:D)
      end if
      
    end subroutine REdiffhh
    
    
    
    subroutine REdiffhT(pde_loc, layer, quadpntin, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
            
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: D
      real(kind=rkind), dimension(3,3) :: KlT, KvT
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_RE_consitutive::REdiffhT"
        ERROR STOP
      end if
      
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      if (present(scalar)) then
        print *, "runtime error evap_RE_consitutive::REdiffhT"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      if (present(tensor)) then
        call cond_hT(layer, quadpnt, KlT(1:D, 1:D))
        call cond_vT(layer, quadpnt, KvT(1:D, 1:D))
        tensor(1:D, 1:D) = Klt(1:D,1:D) + KvT(1:D, 1:D)
      end if
      
      
    end subroutine REdiffhT
    
    
    function REcapacityhh(pde_loc, layer, quadpntin, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use re_constitutive
      use evapglob
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      type(integpnt_str) :: quadpnt
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_RE_consitutive::REcapacityhh"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      val = vangen_elast(pde(re_ord),layer, quadpnt) + dthetav_dh(pde(re_ord),layer, quadpnt) 
    end function REcapacityhh
    
    
    
        
    function REcapacityhT(pde_loc, layer, quadpntin, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      type(integpnt_str) :: quadpnt
      
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_RE_consitutive::REcapacityhT"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      val = dthetav_dtemp(pde(re_ord),layer, quadpnt) 
      
      
    end function REcapacityhT
  
  
    !> hydraulic conductivity for vapour
    !! \f[  K_{vh} = \frac{D}{\rho_{w}} \rho_{sv} \frac{Mg}{RT}H_{r}, \f]
    !! where
    !! \f[ D = \tau (\theta_{s} - \theta_{l}) D_{a},  \f],
    !! and
    !! \f[ \tau = \frac{(\theta_{s}-\theta_{l})^{\nicefrac{7}{3}}}{\theta_{s}^{2}}. \f]
    !<
    subroutine cond_vapour4h(layer, quadpntin, Kvh)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use globals
      
      !> MaterialID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpntin 
      !>unsaturated non-thermal conductuvity of water vapor
      real(kind=rkind), dimension(:,:), intent(out) :: Kvh 

      real(kind=rkind) :: val, tort, ths, theta, Da, D, T
      integer(kind=ikind) :: i
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))

      val = vapour_diff(layer, quadpnt) / dens_liquid(quadpnt) * MolWat * gravity / (T*R_gas)*relhumid(quadpnt)*dens_satvap(quadpnt)
      
      Kvh = 0

      do i = 1, drutes_config%dimen
        Kvh(i,i) = val
      end do
  
    end subroutine cond_vapour4h
    
    
    !> cross factor term heat in liquid
    !! \f[ K_{lt} = K_{lh} \left( h G_{w} \frac{1}{\gamma_{0}} \dv{\gamma}{T} \right), \f]
    !! \f[  \dv{\gamma}{T} = -0.1425 - \num{4.76e-4} T. \f]
    !<
    subroutine cond_ht(layer, quadpntin, Klht)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use globals
      use re_constitutive
    
      !> MaterialID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpntin 
      !>unsaturated non-thermal conductuvity of water vapor
      real(kind=rkind), dimension(:,:), intent(out) :: Klht
      
      real(kind=rkind) :: dgamma, T, h
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      T = pde(heat_ord)%getval(quadpnt)
      
      h = pde(re_ord)%getval(quadpnt)
      
      dgamma = -0.1425 - 4.76e-4*T
      
      call mualem(pde(re_ord), layer, quadpnt, tensor=Klht)
      
      Klht = Klht * h * dgamma * 5 * 1.0_rkind/71.89 * dgamma 

    end subroutine cond_ht
    
    
    !> cross factor term heat in liquid for vapour
    !! \f[ K_{vt} = \frac{D}{\rho_{l}} \eta H_{r} \dv{\rho_{sv}}{T},  \f]
    !! \f[  \eta = 9.5 + 3\frac{\theta_{l}}{\theta_{s}} - 8.5 \exp \left( - \left( \left( 1 + \frac{2.6}{\sqrt{f_{c}}} \right) \frac{\theta_{l}}{\theta_{s}} \right)^{4} \right). \f]
    !! mass fraction of clay \f[ f_c = 0.04 \f]
    !<
    subroutine cond_vt(layer, quadpntin, Klvt)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use globals
      use re_constitutive
    
      !> MaterialID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpntin 
      !>unsaturated non-thermal conductuvity of water vapor
      real(kind=rkind), dimension(:,:), intent(out) :: Klvt
      
      real(kind=rkind) :: h, ths, theta, eta, val
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: i
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      h = pde(re_ord)%getval(quadpnt)
      
      ths = vgset(layer)%ths

      theta = vangen(pde(re_ord), layer, quadpnt)
      
      eta = 9.5 + 3*theta/ths - 8.5*exp(-((1+2.6/sqrt(fraction_clay)) * theta/ths)**4)
    
      val = vapour_diff(layer, quadpnt) / dens_liquid(quadpnt) * relhumid(quadpnt) * drhosv_dT(quadpnt) * eta      
      
      
      Klvt = 0
      
      do i=1, drutes_config%dimen
        Klvt(i,i) = val
      end do
      
      
    end subroutine cond_vt
  
    
    function thetav(pde_loc, layer, quadpntin, x) result(val)
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
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: theta, theta_s
      type (integpnt_str) :: quadpnt
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_RE_constitutive::thetav"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      theta_s = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      
      val = (theta_s - theta)*relhumid(quadpnt) * dens_satvap(quadpnt) / dens_liquid(quadpnt)
    
    end function thetav
    
    !> derivative of the vapour water content with respect to temperature
    !! \dv{\theta_{v}}{T} = (\theta_{s} - \theta_{l})\left( \dv{H_{r}}{T} \frac{\rho_{sv}}{\rho_{l}} + H_{r}  \dv{\rho_{sv}}{T} \frac{1}{\rho_{l}} + H_{r} \rho_{sv} \dv{\frac{1}{\rho_{l}}}{T} \right) \f]
    !<
    function  dthetav_dtemp(pde_loc,layer, quadpntin, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use re_globals
      use evapglob
      use re_constitutive
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      ths = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      
      val = (ths-theta)*(drelhumiddT(quadpnt)*dens_satvap(quadpnt)/dens_liquid(quadpnt) + &
            relhumid(quadpnt) * drhosv_dT(quadpnt)/dens_liquid(quadpnt) + &
            relhumid(quadpnt)*dens_liquid(quadpnt)*invdrhol_dT(quadpnt))
            
    
    end function dthetav_dtemp
    
    
    function dthetav_dh(pde_loc,layer, quadpntin, x) result(val)
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
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta, Cl
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
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
    function relhumid(quadpntin) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpntin
      real(kind=rkind) :: val
      
      real(kind=rkind) :: h, T
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
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
    function drelhumiddh(quadpntin) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpntin
      real(kind=rkind) :: val
      
      real(kind=rkind) :: h, T
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
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
    function drelhumiddT(quadpntin) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpntin
      real(kind=rkind) :: val
      
      real(kind=rkind) :: h, T
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
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
    function drhosv_dT(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) ::  T   
     
      T = pde(heat_ord)%getval(quadpnt)        
    
	  if (T > 0) then
		val = -((317.0*T*T+40000.0*T-240591600.0_rkind)*exp(-(317.0*T)/40000.0-601479.0/(100.0*T)+78429.0/2500.0))& 
				/(400000000.0*T*T*T)
	  else
		val = 0
	  end if
      
    end function drhosv_dT
    
    
    !>  derivative of the term \f[ \rho_{l}^{-1} \f] with a respect to temperature is given by
    !! \f[ -\dfrac{3c\left(x-r\right)^2-2b\left(x-r\right)}{\left(c\left(x-r\right)^3-b\left(x-r\right)^2+a\right)^2} \f]
    !<
    function invdrhol_dT(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) ::  x, a,b,c,r   
      
      
      x = pde(heat_ord)%getval(quadpnt)
      a=1000.0
      b=7.37e-3
      c=3.79e-5
      r=3.98
      
      val =-(3*c*(x-r)**2-2*b*(x-r))/(c*(x-r)**3-b*(x-r)**2+a)**2
              
    end function invdrhol_dT


    !> vapour soil diffusivity
    !! \f[ D = \tau(\theta_S-theta) D_a \f]
    !! \f[ \tau = \frac{(\theta_{s}-\theta_{l})^{\nicefrac{7}{3}}}{\theta_{s}^{2}} \f]
    !! \f[  D_{a} = \num{2.12e-5} \left(\frac{T}{273.15} \right)^{2}. \f]
    !<
    function vapour_diff(layer, quadpntin) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use re_constitutive
      
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value: Vapor Difussivity in soil  [m^2/s]
      real(kind=rkind) :: val
      !> Water content,Volumetric air content,
      real(kind=rkind) :: theta_l, theta_air, Da, tort, T
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      T = T2kelv(pde(heat_ord)%getval(quadpnt))    
  
      theta_air = vgset(layer)%Ths - vangen(pde(re_ord), layer, quadpnt)
      
      tort = ((theta_air)**(7.0/3.0))/ (vgset(layer)%Ths**2)
      
      Da = 2.12e-5*(T/273.15)**2
      
      val = tort*theta_air*Da

    end function vapour_diff
    
    !> derivative of vapour diffusivity with respect to \f[ h \f]
    !! \f[ \dv{K_{vh}}{h} = {\rho_{w}} \rho_{sv} \frac{Mg}{RT}H_{r} \dv{D}{h}. \f]
    !! \f[ \begin{split} D(h) &= \frac{\theta_s - \theta_l)^{7/3}}{\theta_s^2} \left(\theta_s - \theta_l) \right) . D_a = \\ &=  \frac{\left( \theta_s - \frac{\theta_s - \theta_r}{(1+(-\alpha h)^n)^m} +\theta_r \right)^{7/3}}{\theta_s^2} \left(\theta_s - \frac{\theta_s - \theta_r}{(1+(-\alpha h)^n)^m} +\theta_r\right) . D_a = \\ &= D_a \left(-\dfrac{\theta_s-\theta_r}{\left( 1 + \left(-\alpha h\right)^n \right)^m}+\theta_s+\theta_r \right)^\frac{10}{3} \end{split} \f]
    !! \f[ \dv{D}{h} = \dfrac{10 D_a mn\cdot\left(\theta_s - \theta_r \right)\left(-\alpha h \right)^n\left(\left(-\alpha h\right)^n+1\right)^{-m-1}\left(-\frac{\theta_s-\theta_r}{\left(\left(-\alpha h \right)^n+1\right)^m}+\theta_s + \theta_r \right)^\frac{7}{3}}{3h}. \f]
    !<
    subroutine dKvhdh(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use pde_objs
      use re_globals
      use core_tools
      use evapglob

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

      real(kind=rkind) :: a,n,m, ths, thr, h
      type(integpnt_str) :: quadpnt_loc 
      real(kind=rkind) ::  dDdh, absval, Da, T
      
      
!      quadpnt_loc = quadpnt
!      quadpnt_loc%preproc = .true.
      
!      h = pde(re_ord)%getval(quadpnt_loc)
      
!      a = vgset(layer)%alpha
!      n = vgset(layer)%n
!      m = vgset(layer)%m
      
!      T = T2kelv(pde(heat_ord)%getval(quadpnt_loc))  
      
!      Da = 2.12e-5*(T/273.15)**2  
      
!      dDdh = (10*Da*m*n*(ths-thr)*(-a*h)**n*((-a*h)**n+1)**(-m-1)*(-(ths-thr)/((-a*h)**n+1)**m+ths+thr)**(7.0_rkind/3))/(3*h)
	
!	  absval = dDdh / dens_liquid(quadpnt_loc) * MolWat * gravity / (T*R_gas)*relhumid(quadpnt_loc)*dens_satvap(quadpnt_loc)
      absval = 0
	  
	  if (present(scalar)) then
		scalar = absval
	  end if
	  
	  if (present(vector_out)) then
		vector_out = 0
		!should be uncommented and tested!
!		vector_out(drutes_config%dimen) = absval
	  end if
	  
	end subroutine dKvhdh
          
    
      !> Liquid water flux
      !! \f[ \vec{q} = - \mathbf{K}_ll (\nabla h + \nabla z) - \mathbf{K}_{lT} \nabla T \f]
      !<
    subroutine darcy4liq(pde_loc, layer, quadpntin, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use re_constitutive
      use evapglob
      use debug_tools
      
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      
      
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> Klt: total unsaturated non-thermal conductivity of liquid water
      real(kind=rkind), dimension(3,3)  :: KlT, Klh
      !> local variable
      integer(kind=ikind)  :: i, D
      !> Liquid water flux
      real(kind=rkind), dimension(3)  ::  q_liq
      !> resul of the modified flux of liquid water
      real(kind=rkind), dimension(3)  :: vct
      !> Temperature gradient
      real(kind=rkind), dimension(:), allocatable, save :: gradT, gradH
      type(integpnt_str) :: quadpnt
      
      if (.not. present(quadpntin)) then
        print *, "ERROR: use quadpnt only"
        print *, "exited from evap_RE_constitutive::darcy4liq"
        ERROR STOP
      end if
      
      quadpnt = quadpntin

      call pde(heat_ord)%getgrad(quadpnt, gradT)
      
      call pde(re_ord)%getgrad(quadpnt, gradH)
      
      quadpnt%preproc = .true.

      
!      if (quadpnt%type_pnt=="obpt") then
!        print *, gradH
!        call wait()
!      end if

      D = drutes_config%dimen
      
!      if (drutes_config%dimen == 1) then
!        gradH(D) = gradH(D) + cos(vgset(layer)%anisoangle(1))
!      else
!        gradH(D) = gradH(D) + 1
!      end if
      
      call mualem(pde(re_ord), layer, quadpnt, tensor=Klh(1:D, 1:D))
      
      call cond_ht(layer, quadpnt, KlT(1:D,1:D)) 
      
      vct(1:D) = matmul(-Klh(1:D, 1:D), gradH) + matmul(-KlT(1:D,1:D), gradT(1:D))
       
      
      if (present(flux_length)) then      
        flux_length = norm2(vct(1:D))
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    
    end subroutine darcy4liq
        
        
      !> Vapour flux
      !! \f[ \vec{q} = - \mathbf{K}_vh \nabla h  - \mathbf{K}_{vT} \nabla T \f]
      !<
    subroutine darcy4vap(pde_loc, layer, quadpntin, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evapglob
      use debug_tools
      use re_globals
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> KvT: unsaturated thermal conductivity for water
      !> Kvh: unsaturated non-thermal  conductivity for water
      real(kind=rkind), dimension(3,3)  :: Kvh, KvT
      !> local variables
      integer(kind=ikind) :: i, D
      !> pressure gradient
      real(kind=rkind), dimension(:), allocatable :: gradh,  gradT
      
      !result of the vapor flux vector
      real(kind=rkind), dimension(3) :: vct
      type(integpnt_str) :: quadpnt


      
      if (present(x)) then
        print *, "ERROR: use quadpnt only"
        print *, "exited from evap_fnc::liquid_flux"
        ERROR stop
      end if
      
      quadpnt = quadpntin
    
      call pde(heat_ord)%getgrad(quadpnt, gradT)

      call pde(re_ord)%getgrad(quadpnt, gradh)
      
      
      D = drutes_config%dimen
      quadpnt%preproc = .true.
      
      if (drutes_config%dimen == 1) then
        gradh(D) = gradh(D) - cos(vgset(layer)%anisoangle(1))
      else
        gradh(D) = gradh(D) - 1
      end if
      
      call cond_vapour4h(layer, quadpnt, Kvh(1:D, 1:D))
      call cond_vt(layer, quadpnt, KvT(1:D, 1:D))
      

      
      
      vct(1:D) =  matmul(-Kvh(1:D,1:D), gradh(1:D)) + matmul(-KvT(1:D,1:D), gradT(1:D))
      
      
       if (present(flux_length)) then
         flux_length = norm2(vct(1:D))
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    end subroutine darcy4vap
    
    
    subroutine totalflux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evapglob
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      
      integer(kind=ikind) :: D
      real(kind=rkind), dimension(3) :: vctvapour, vctliquid
      
      D = drutes_config%dimen
      
      call darcy4liq(pde(re_ord), layer, quadpnt, flux=vctliquid(1:D))
      
      call darcy4vap(pde(re_ord), layer, quadpnt, flux=vctvapour(1:D))
    
      
      if (present(flux)) flux(1:D) = vctliquid(1:D) + vctvapour(1:D)
      
      if (present(flux_length)) flux_length = norm2(vctliquid(1:D) + vctvapour(1:D)) 
      
      
    end subroutine totalflux
    
    
    
    
    


end module evap_RE_constitutive
