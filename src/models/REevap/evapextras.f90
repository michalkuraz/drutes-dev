module evapextras

  public :: rh_soil, rho_l


  contains
  
   !< Relative Humudity soil rh_soil [-]
    !Input: pressure head: h [m]
    !Temperature: T [ºC]
    !Garvity : grvity  [m.s^-2]
    !Molecular weight of water: MolWat [kg mol^-1]
    !Universal Gas constant: R_gas [J mol^-1 K^-1]
    function rh_soil(layer, quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value: Relative Humudity soil rh_soil [-]
      real(kind=rkind):: val
      !> Pressure head: h, Temperature T in ºC and TemperatureT_abs in Kelvin
      real(kind=rkind):: h,T, T_abs
      
     
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not integ point "
        print *, "exited from evap_auxfnc::rh_soil"
        ERROR stop
      end if
      
      h = pde(RE_ord)%getval(quadpnt)
      T = pde(Heat_ord)%getval(quadpnt)
      T_abs = T + Tref ! T from ºC to Kelvin
      
      val = exp ((h*MolWat*gravity)/(R_gas*T_abs))
      
    end function rh_soil
    
    
      
        !Liquid water density  rho_l [kg/m^3]
    !Input: Temperature in ºC
    function rho_l( quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value:Liquid water density  rho_l [kg/m^3]
      real(kind=rkind):: val
      !>  Temperature T in ºC 
      real(kind=rkind)::T
      
       if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified integ point "
        print *, "exited from evap_auxfnc::rho_l"
        ERROR stop
      end if
    
      T = pde(Heat_ord)%getval(quadpnt)
      
      val = 1000.0_rkind - 7.37e-3*(T - 4.0_rkind)**2 + 3.79e-5*(T -4.0_rkind)**3

    end function rho_l
    
      
!      Saturated water vapor density rho_sv [kg/m^3]
    !Input: Temperature in T [ºC]
    function rho_sv(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob 
      use re_globals
    
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value: Saturated water vapor density rho_sv [kg/m^3]
      real(kind=rkind):: val
      !>  Temperature T in ºC and TemperatureT_abs in Kelvin
      real(kind=rkind):: T, T_abs
      
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified integ point "
        print *, "exited from evap_auxfnc::rho_sv"
        ERROR stop
      end if
     
      T = pde(Heat_ord)%getval(quadpnt)
      T_abs = T + Tref ! T from ºC to Kelvin
      
      val = 1e-3 *(exp(31.3716_rkind - (6014.79_rkind/T_abs) - 7.92495e-3*T_abs))/T_abs

    end function rho_sv
    
  
    
    
    
      !Vapor Difussivity in soil  [m^2/s]
    !Input: Water content,  Saturated water content, Volumetric air content
    function vapor_diff_soil(pde_loc, layer, quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value: Vapor Difussivity in soil  [m^2/s]
      real(kind=rkind):: val
      !> Water content,Volumetric air content,
      real(kind=rkind)::theta_l, theta_air
      
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified  integ point "
        print *, "exited from evap_auxfnc::vapor_diff_soil"
        ERROR stop
      end if
    
      theta_l = pde(re_ord)%mass(1)%val(pde(re_ord), layer, quadpnt)
      theta_air = vgset(layer)%Ths - theta_l
     
      
      val = tortuosity(theta_l, layer)*theta_air*vapor_diff_air(quadpnt)

    end function vapor_diff_soil
    
    
      !Tortuosity factor in gaseous phase [-]
    !Input: Volumetric liquid water  content [-]
    !Saturated water content [-]
    function tortuosity(theta_l, layer) result(val)
      use typy
      use global_objs
      use pde_objs
      use evap_globals
      use re_globals
      
      !> Volumetric liquid water  content
      real(kind=rkind),intent (in):: theta_l
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value: Tortuosity factor in gaseous phase [-]
      real(kind=rkind):: val
      !> Volumetric air content,  Saturated water content 
      real(kind=rkind):: theta_sat, theta_air
      
      theta_air = vgset(layer)%Ths - theta_l
      theta_sat = vgset(layer)%Ths
      
      val = ((theta_air)**(7/3))/ (theta_sat**2)

    end function tortuosity
    
      !Vapor Difussivity in air  [m^2/s]
    !Input: Temperature [ºC]
    function vapor_diff_air(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value: Vapor Difussivity in air  [m^2/s]
      real(kind=rkind):: val
      !>Temperature T in ºC and TemperatureT_abs in Kelvin
      real(kind=rkind)::T, T_abs
      
       if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified  integ point "
        print *, "exited from evap_auxfnc::vapor_diff_air"
        ERROR stop
      end if
    
      T = pde(Heat_ord)%getval(quadpnt)
      T_abs = T + Tref ! T from ºC to Kelvin
      
      val =  2.12e-5 * (T_abs/Tref)**2

    end function vapor_diff_air
    
    
    !Derivative Surface Tension Soil-Water  [g/s^2 ºC]
    !Input: Temperature in ºC  
    function dsurf_tension_soilwat_dT(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      

      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value: Derivative Surface Tension Soil-Water  [g/s^2 ºC]
      real(kind=rkind):: val
      !>Temperature in ºC 
       real(kind=rkind)::T
      
      
       if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point "
        print *, "exited from evap_auxfnc::dsurf_tension_soilwat_dT"
        ERROR stop
      end if
    
      T = pde(Heat_ord)%getval(quadpnt) + Tref
      
      val = - 0.1425_rkind - 4.76e-4*T

    end function dsurf_tension_soilwat_dT
    
    
    
        !> Water vapor content
    function theta_vapor(pde_loc,layer, quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_constitutive
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      !>material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt 
      !vapor volumetric content
      real(kind=rkind) :: val
      !> Rh: relatuive humidity of soil
      !> rho_l: liquid water density
      !> rho_sv: saturated vapor density
      !>theta_l: liquid water content
      real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val, theta_l
      
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point "
        print *, "exited from evap_auxfnc::theta_vapor"
        ERROR stop
      end if
        
      theta_l = vangen(pde(re_ord), layer, quadpnt)
      rh_soil_val = rh_soil(layer, quadpnt)
      rho_l_val = rho_l(quadpnt) 
      rho_sv_val = rho_sv(quadpnt) 
        
      val = (vgset(layer)%Ths - theta_l)*rho_sv_val*rh_soil_val*(1.0_rkind/rho_l_val)

      
    end function theta_vapor

end module evapextras
