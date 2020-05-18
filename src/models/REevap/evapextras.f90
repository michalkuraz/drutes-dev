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
      use evapglob
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
    
    
      !Specific Latent heat of evaporation of liquid water  [Jkg^-1]
    !Input: Temperature in ºC
    function latent_heat_wat(quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
    
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value: specific Latent heat of evaporation of liquid water  [Jkg^-1]
      real(kind=rkind):: val
      !>  Temperature T in ºC 
      real(kind=rkind)::T
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified  integ point "
        print *, "exited from evap_auxfnc::latent_heat_wat"
        ERROR stop
      end if
    
      T = pde(Heat_ord)%getval(quadpnt) !T is in ºC
      val = 2.501e+6 - 2369.2_rkind*T
      val = val*rho_l(quadpnt)

    end function latent_heat_wat
    
    
      !Enhacement Factor [-]
    !Input: Volumetric liquid water  content [-]
    !Saturated water content [-]
    function enhancement_factor(pde_loc, layer, quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_constitutive
      use re_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind):: val
      !> Volumetric liquid water  content, Saturated water content
      real(kind=rkind)::theta_l, theta_sat, tmp, const
      
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified  integ point "
        print *, "exited from evap_auxfnc::enhacement_factor"
        ERROR stop
      end if
      
      theta_sat = vgset(layer)%ths
      theta_l = vangen(pde(re_ord), layer, quadpnt)
      const = 1 + (2.6_rkind/sqrt(fraction_clay))
      tmp = exp(- (const * (theta_l/theta_sat))**4)
      
      val =  9.5_rkind + 3.0_rkind*(theta_l/theta_sat) -8.5_rkind *tmp
    end function enhancement_factor
    
    
    !Derivative saturated water vapor density rho_sv [kg/m^3 K]
    !Input: Temperature in T [ºC]
    function drho_sv_dT( quadpnt) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> return value : Derivative saturated water vapor density rho_sv [kg/m^3 K]
      real(kind=rkind):: val
      !>  Temperature T in ºC and TemperatureT_abs in Kelvin
      real(kind=rkind):: T, T_abs
      
      if (.not. present(quadpnt)) then
        print *, "ERROR: you have not specified  integ point "
        print *, "exited from evap_auxfnc::drho_sv_dT"
        ERROR stop
      end if
    
      T = pde(Heat_ord)%getval(quadpnt)
      T_abs = T + Tref ! T from ºC to Kelvin
     
      
      val = exp(- 7.92495e-3*T_abs - (6014.79_rkind/T_abs))*((-3.33818e8*T_abs - 4.2122e10)*T_abs + 2.53357e14)*(1/T_abs**3)

    end function drho_sv_dT
    
    
     !> Isothermal Properties of water vapor
      function hydraulic_vh(layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evapglob
        use re_globals
        
        !> MaterialID
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        !>unsaturated non-thermal conductuvity of water vapor
        real(kind=rkind) :: val
        !> T:temperature
        !> Rh: relatuive humidity of soil
        !> rho_l: liquid water density
        !> rho_sv: saturated vapor density
        !> D: diffusitivity
        real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val,diff, T
        
        
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified  integ point "
          print *, "exited from evap_fnc::hydraulic_vh"
          ERROR stop
        end if
       
        
        rh_soil_val = rh_soil(layer, quadpnt)
        rho_l_val= rho_l( quadpnt) 
        rho_sv_val = rho_sv( quadpnt) 
        diff = vapor_diff_soil(pde(heat_ord), layer, quadpnt)
        
        T = pde(Heat_ord)%getval(quadpnt)
				T = T + Tref
       
        val = (diff/rho_l_val)*rho_sv_val*((MolWat*gravity)/(R_gas*T))*rh_soil_val
    
      end function hydraulic_vh
      
      
      
      
          !> Thermal Properties of Liquid water
      function hydraulic_lT(layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evapglob
        use re_constitutive

        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> unsaturated thermal conductivity for liquid water
        real(kind=rkind), dimension(3, 3) :: val
        !> T:temperature
        !> h: pressure head
        real(kind=rkind) :: T, h
        !unsaturated non'thermal conductivity
        real(kind=rkind), dimension(3,3) :: Klh
        !local variable: dimension
        integer(kind=ikind):: D
        
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified  integ point "
          print *, "exited from evap_fnc::hydraulic_lT"
          ERROR stop
        end if
        
        
        D = drutes_config%dimen
        h = pde(RE_ord)%getval(quadpnt)
        T = pde(Heat_ord)%getval(quadpnt)
       
        call mualem(pde(re_ord), layer, quadpnt,  tensor = Klh(1:D,1:D))
        val(1:D,1:D)  = Klh(1:D,1:D)*h*GwT*(1/gamma_0)*dsurf_tension_soilwat_dT(quadpnt)  
        
    end function hydraulic_lT
    
    
      
  
      
    !> Thermal Properties of water vapor
      function hydraulic_vT(layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evapglob
        
        !Material ID
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        !> unsaturated thermal hydraulic conductivity for water vapor
        real(kind=rkind) :: val
        !> Rh: relatuive humidity of soil
        !> rho_l: liquid water density
        !> rho_sv: saturated vapor density
        !> D: diffusitivity
        !> enhancement factor
        real(kind=rkind) :: rh_soil_val, rho_l_val,drho_svdT_val,diff, enhancement_factor_val
        
        
        rh_soil_val= rh_soil( layer, quadpnt)
        rho_l_val = rho_l(quadpnt) 
        diff = vapor_diff_soil(pde(re_ord), layer, quadpnt)
        drho_svdT_val = drho_sv_dT(quadpnt)
        enhancement_factor_val = enhancement_factor(pde(re_ord), layer, quadpnt)
        
        val = (diff/rho_l_val)*enhancement_factor_val*drho_svdT_val*rh_soil_val
        
      end function hydraulic_vT


end module evapextras
