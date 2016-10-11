!> @file
!! Parametrické funkce

!> @brief  Modul obsahuje funkce parametru vstupujících do výpočtu
!!
!! Module obsahuje funkce parametru vstupujících do výpočtu dle článku \cite noborio.
!! 
!! @author Jakub Jerabek
!! @version TEST

module modRE_parameter_functions


  public :: eVol 
  public :: viscosity
  public :: mol_conc
  public :: water_activity_polutant
  public :: surface_tension_polutant
  public :: derT_surface_tension_polutant
  public :: derC_surface_tension_polutant
  public :: osmotic_coef
  public :: relativ_humidity_soil
  public :: saturated_vapor_density
  public :: der_saturated_vapor_density
  public :: mechanics_enhancement_factor
  public :: molecular_D_vapor_water
  public :: rho_pure_water
  public :: latent_heat_v
  public :: spec_heat_solut
  public :: diff_coef_NaCl
 contains
  !> @brief Výpočet efektivní objemové vlhkosti.
  !!
  !! Do rovnice vstupuje objemová vlhkost @f$\theta@f$, vstupni hodnota do rovnic "z venku", reziduálni objemová vlhkost @f$\theta_r@f$ (@f$~\ref{par:the_r}@f$) a saturovaná objemová vlhkost @f$\theta_r@f$ (@f$\ref{par:the_r}@f$).
  !!
  !! @f[
  !! \theta_e=\frac{\theta - \theta_r}{\theta_s - \theta_r}
  !! \label{sec:evol}
  !! @f]
  !!
  !! @param[in] layer layer identificatr
  !! @param[in]  thet @f$\theta@f$, volumatic moisture
  !! @param[out]  thete @f$\theta_e@f$,   Nenasycená hydraulická vodivost plynné fáze
  
  
  
  function eVol(layer, theta) result(thete) 
    use typy
    use modRE_globals
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: theta
    real(kind=rkind) :: thete

    real(kind=rkind) :: the_s, the_r

    the_r = hyd_prop(layer)%theta_r
    the_s = hyd_prop(layer)%theta_s

    thete  = (theta - the_r)/(the_s - the_r)

  end function eVol
  
  
  
  !> @brief Výpočet kinematické viskosity solného roztoku.
  !!
  !! Vztah vyjedřije závislost kinematické viskosity solného roztoku na teplotě a koncentraci.
  !!
  !! @f[
  !! \nu = 10^{-1.469+200.93/(T+116.71)} \cdot 10^{-3} (1 + 9.359 \cdot 10^{-2} C + 9.320 \cdot 10^{-3} C^{2})
  !! \label{eq:vis}
  !! @f]
  !!
  !! @param[in] T @f$ T @f$, teplota ve @f$ ^{\circ}C @f$ (@f$\ref{par:tem}@f$)
  !! @param[in] C @f$ C @f$, koncentrace v @f$ mol/kg @f$ (@f$\ref{par:conc}@f$)
  !! @param[out] vis @f$ \nu @f$, kinematická viskosita
  
  
  function viscosity(T, C, layer) result(vis)
    
    use typy
    use modRE_globals
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: T, C
    real(kind=rkind) :: vis, wrkT, wrkC
    wrkT = T - 273.15_rkind
    wrkC = mol_conc(C, layer)
    vis = 10.0_rkind**(-1.469_rkind+200.93_rkind/(wrkT+116.71_rkind))    &
    * 0.001_rkind*(1.0_rkind + 9.359E-002_rkind*wrkC    & 
    + 9.320E-003_rkind*wrkC**2.0_rkind)
  end function viscosity
  
  !> @brief Přepočet koncentrace roztoku z jednotek @f$ kg/kg @f$ na @f$ mol/kg @f$ pomocí molární hmotnosti @f$ NaCl @f$ (@f$\ref{par:weiNaCl}@f$).
  !!
  !! @f[
  !!  C = \frac{C_{kg/kg}}{M_{salt}(1 - C_{kg/kg})}
  !!  \label{eq:mol_conc}
  !! @f]
  !!
  !! @param[in] conc @f$ C_{kg/kg} @f$, koncentrace @f$ kg/kg @f$ (@f$\ref{par:conc}@f$)
  !! @param[out] mol_conc @f$ C @f$, koncentrace @f$ mol/kg @f$
  
  

  function mol_conc(conc, layer) result(res)
    use typy
    use modRE_globals
    real(kind=rkind), intent(in) :: conc 
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind) :: res

    real(kind=rkind) :: Msalt
    
    Msalt = hyd_prop(layer)%molecular_weight_NaCl

    res = conc/(Msalt*(1.0_rkind-conc))

  end function mol_conc

  !> @brief Vztah pro určení aktivity rozpuštěného NaCl ve vodě.
  !!
  !! @f[
  !!  h_\pi = 1 - 3.1327 \cdot 10^{-2} C - 1.4553 \cdot 10^{-3} C^2
  !!  \label{eq:activityW}
  !! @f]
  !!
  !! @param[in] C @f$ C @f$, koncentrace @f$ mol/kg @f$ (@f$\ref{par:conc}@f$)
  !! @param[out] hpi @f$ h_{\pi} @f$, aktivity rozpuštěného NaCl ve vodě
  
  function water_activity_polutant(C, layer) result(hpi)
    use typy
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: C
    real(kind=rkind) :: hpi, wrkC
    wrkC = mol_conc(c, layer)
    hpi = 1.0_rkind - 3.1327e-2_rkind*wrkC - 1.4553e-3_rkind*wrkC**2.0_rkind
  end function water_activity_polutant
  
  function der_water_activity_polutant(C, layer) result(dhpi)
    use typy
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: C
    real(kind=rkind) :: hpi, wrkC
    real(kind=rkind) :: dhpi
    wrkC = mol_conc(c, layer)
    hpi = 1.0_rkind - 3.1327e-2_rkind*wrkC - 1.4553e-3_rkind*wrkC**2.0_rkind
    dhpi = (log(hpi)/wrkC)*hpi
  end function der_water_activity_polutant
  
  
  !> @brief Vztah pro určení povrchového napětí roztoku NaCl.
  !!
  !! @f[
  !!  \gamma = 7.5617 \cdot 10^{-2} - 1.3595 \cdot 10^{-4}T - 4.0815 \cdot 10^{-7}T^2 + 1.6342 \cdot 10^{-3}C
  !!  \label{eq:gamma}
  !! @f]
  !!
  !! @param[in] C @f$ C @f$, koncentrace @f$ mol/kg @f$  (@f$\ref{par:conc}@f$)
  !! @param[in] T @f$ T @f$, teplota @f$ ^\{circ}C @f$ (@f$\ref{par:tem}@f$)
  !! @param[out] gamma @f$ \gamma @f$, povrchové napětí roztoku NaCl
  function surface_tension_polutant(T, C, layer) result(gamma)
    use typy
    
    real(kind=rkind), intent(in) :: T, C
    real(kind=rkind) :: gamma, wrkT, wrkC
    integer(kind=ikind) :: layer
!     if ((T >= 50.0_rkind) .and. (T<100.0_rkind)) then
!       print*, "" //achar(27)//'[91m', "Špatně zadana teplota, musi byt v kelvinech!", "" //achar(27)//'[0m'
!     else if (T>= 100.0_rkind) then
!       print*, "" //achar(27)//'[91m', "Špatně zadana teplota, musi byt v kelvinech!" , "" //achar(27)//'[0m'
!       STOP
!     else
!       continue
!     end if
    wrkC = mol_conc(c, layer)
    wrkT = T - 273.15_rkind
    gamma = 7.5617e-2_rkind - 1.3595e-4_rkind*wrkT - 4.0815e-7_rkind*wrkT**2.0_rkind   &
	    + 1.6342e-3_rkind*wrkC
    
  end function surface_tension_polutant
  
  
  !> @brief Derivate vztah pro určení povrchového napětí roztoku NaCl podle @f$ T @f$.
  !!
  !! @f[
  !!  \frac{\mathrm{d} \gamma}{\mathrm{d} T} = - 1.3595 \cdot 10^{-4} - 2(4.0815 \cdot 10^{-7}T)
  !!  \label{eq:gammaT}  
  !! @f]
  !!
  !! @param[in] T @f$ T @f$, teplota @f$ ^\{circ}C @f$ (@f$\ref{par:tem}@f$)
  !! @param[out] gamma_derT @f$ \gamma @f$, povrchové napětí roztoku NaCl
  function derT_surface_tension_polutant(T) result(gamma_derT)
    use typy
    
    real(kind=rkind), intent(in) :: T
    real(kind=rkind) :: gamma_derT, wrkT
    wrkT = T - 273.15_rkind
    gamma_derT = - 1.3595e-4_rkind - 2.0_rkind*4.0815e-7_rkind*wrkT
    
   end function derT_surface_tension_polutant 
  
  !> @brief Derivate vztah pro určení povrchového napětí roztoku NaCl podle @f$ C @f$
  !!
  !! @f[
  !!  \frac{\mathrm{d} \gamma}{\mathrm{d} C} = 1.6342 \cdot 10^{-3}
  !!  \label{eq:gammaC} 
  !! @f]
  !! @param[in] C @f$C@f$ koncentrace @f$mol/kg@f$ (@f$\ref{par:conc}@f$)
  !! @param[out] gamma_derC @f$ \frac{\mathrm{d} \gamma}{\mathrm{d} C} @f$, povrchové napětí roztoku NaCl
  function derC_surface_tension_polutant(C) result(gamma_derC)
    use typy
    
    real(kind=rkind), intent(in) :: C
    real(kind=rkind) :: gamma_derC
    
    gamma_derC =  1.6342e-3_rkind
   end function derC_surface_tension_polutant 
   
  
  !> @brief Výpočet koeficientu vlivu osmotických sil
  !!
  !! Do rovnice vstupuje poloměr iontu polutantu @f$r_s@f$ (@f$\ref{par:rs}@f$), poloměr molekuly vody @f$r_w@f$ (@f$\ref{par:rw}@f$), poloviční tloušťka filmu polutantu na částěčkách půdy @f$b@f$ (@f$\ref{par:b}@f$).
  !!
  !! @f[
  !! \sigma = \frac{r_s - r_w}{b - r_w}
  !! \label{eq:osmotic}
  !! @f]
  !!
  !! @param[out] sigma @f$ \sigma @f$, koeficient efektivity osmotických sil
  
  function osmotic_coef(layer) result(sigma)	
    use typy
    use modRE_globals
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind) :: sigma
    real(kind=rkind) :: rs, rw, b
    
    rs = hyd_prop(layer)%hydrated_iont_r
    rw = hyd_prop(layer)%water_molecul_r
    b = hyd_prop(layer)%half_of_solute_film 
    sigma = (rs - rw)/(b - rw)
  end function osmotic_coef
  
  !> @brief Výpočet relativní vlhkosti půdního vzduchu
  !!
  !! @f[
  !! h_m = \exp{\left( M_w\frac{\Psi}{RT} \right)}
  !! \label{eq:relat_hum}
  !! @f]
  !!
  !! @param[in] pot @f$ \Psi @f$, potenciál (@f$\ref{eq:poten}@f$)
  !! @param[in] R @f$ R @f$, plynova konstanta (@f$\ref{par:gasCon}@f$)
  !! @param[in] Mw @f$ M_w @f$, molární hnotnost vody (@f$\ref{par:weiNaCl}@f$)
  !! @param[in]  T @f$ T @f$, teplota (@f$\ref{par:tem}@f$)
  !! @param[out] hm @f$ h_m @f$, relativní vlhkost půdních par 

  function relativ_humidity_soil(pot, R, Mw, T) result(hm)
    use typy
    real(kind=rkind), intent(in) :: pot, R, Mw, T
    real(kind=rkind) :: hm
    hm = exp(Mw*(pot/(R*T)))
  end function relativ_humidity_soil
  
  !> @brief Výpočet hustoty plně saturovanách půdních par 
  !!
  !! @f[
  !! \rho^{sat}_v = \frac{\exp{(31.3716 - 6014.79/T - 7.92495 \cdot 10^{-3} T})}{T} 10^{-3}
  !! \label{eq:vap_dens}
  !! @f]
  !!
  !! @param[in] T @f$ T @f$, teplota @f$ K @f$ (@f$\ref{par:tem}@f$)
  !! @param[out] rho_sat @f$ \rho^{sat}_v @f$, hustota plně saturovanách půdních par 
    

  function saturated_vapor_density(T) result(rho_sat)
    use typy
 
    real(kind=rkind), intent(in) :: T  ! musi byt v kelvinech!!!
    real(kind=rkind) :: rho_sat

    rho_sat = 10e-3_rkind*(exp(31.3716_rkind - 6014.79_rkind/T -    &
	      7.92495e-3_rkind*T))/T

  end function saturated_vapor_density
  
 
  !> @brief Výpočet hodnoty derivace hustoty plně saturovaného půdního vzduchu dle teploty
  !!
  !! @f[
  !! \frac{\mathrm{d} \rho^{sat}_v}{\mathrm{d} T } = \frac{\exp{(- 6014.79/T -0.00792495 T)} (2.53357 \cdot 10^{14} - 4.21224\cdot 10^{10}T -3.33818\cdot 10^{8}T^{2})}{T^3}
  !! \label{eq:der_vap_dens}
  !! @f]
  !!
  !! @param[in] T  @f$ T @f$, teplota (@f$\ref{par:tem}@f$)
  !! @param[out] der_rho_sat @f$ \frac{\mathrm{d} \rho_{sat}}{\mathrm{d} T } @f$ derivace ---
  
  
!     \begin{array}{c c}
!   
!   \frac{\mathrm{d} \rho^{sat}_v}{\mathrm{d} T } = \frac{\exp{(31.3716 - 6014.79/T -7.92495 \cdot 10^{-3}T)}(6014.79/T^2-7.92495 \cdot 10^{-3})T}{T^2}  \cdot 10^{-3} -\\
!   - \frac{\exp{(31.3716 - 6014.79/T - 7.92495 \cdot 10^{-3} T)} }{T^2} \cdot 10^{-3}
!   \end{array}  
!   \label{eq:der_vap_dens}
  function der_saturated_vapor_density(T) result(der_rho_sat)
    use typy

    real(kind=rkind), intent(in) :: T  ! musi byt v kelvinech!!!
    real(kind=rkind) :: der_rho_sat
!     derived by Wolfram Alpha
    der_rho_sat = (exp(-(6014.79_rkind/T) + T*(-0.00792495_rkind))   &
		   *(2.53357e14_rkind - 4.21224e10_rkind*T -3.33818e8_rkind*T*T))/(T*T*T)
		   
!     der_rho_sat = 10e-3_rkind*(exp(31.3716_rkind - 6014.79_rkind/T -  &
! 	7.92495e-3_rkind*T)*(6014.79_rkind/T**2.0_rkind-7.92495e-3_rkind)*T   &
! 	- exp(31.3716_rkind - 6014.79_rkind/T - 7.92495e-3_rkind*T))/T**2.0_rkind

  end function der_saturated_vapor_density
  
  !> @brief Výpočet faktoru mechanického zatížení.
  !!
  !! @f[
  !! \eta = 9.5 + 6 \theta - 8.5 \exp(-[(1 + 2.6/f_c^{0.5})\theta]^4)
  !! \label{eq:mech_enha}
  !! @f]
  !!
  !! @param[in] theta @f$ \theta @f$, objemová vlhkost 
  !! @param[in] clay_fraction @f$ f_c @f$, obsah jílu v půdě (@f$\ref{par:clay}@f$)
  !! @param[out] eta @f$ \eta @f$, parametr vlivu obsahu jílu v půdě (@f$\label{par:clay}@f$)

  function mechanics_enhancement_factor(theta, clay_fraction) result(eta)
    use typy

    real(kind=rkind), intent(in) :: theta, clay_fraction
    real(kind=rkind) :: eta

    eta = 9.5_rkind + 6.0_rkind*theta - 8.5_rkind*   &
	  exp(-((1.0_rkind + 2.6_rkind/(sqrt(clay_fraction)))*theta)**4.0_rkind)

  end function mechanics_enhancement_factor
  
  
  !> @brief Výpočet difůze půdního vzduchu
  !!
  !! @f[
  !! D_g = 229 \cdot 10^{-7} \left( \frac{T}{273.15}  \right)^{1.75}
  !! \label{eq:bin_gas}
  !! @f]
  !! 
  !! @param[in] T  @f$ T @f$, teplota (@f$\ref{par:tem}@f$)
  !! @param[out] Dg @f$ D_g @f$, koeficient difůze pudního vzduchu
  
  
  function molecular_D_vapor_water(T) result(Dg)
    use typy
    real(kind=rkind), intent(in) :: T
    real(kind=rkind) :: Dg
    
!     Dg = 229E-007_rkind*(T/273.15_rkind)**1.75_rkind
    Dg = 2.12e-5*(T/273.15_rkind)**2.0_rkind
  
  end function molecular_D_vapor_water
  
  
  
  
  
  
  
  
  
  
  function rho_pure_water(T) result(res)
    use typy
    real(kind=rkind), intent(in) :: T
    real(kind=rkind) :: res, wrkT
    wrkT = T - 273.15_rkind !prevod na stupne
    res = (1.0_rkind - (wrkT - 3.9863_rkind)**2.0_rkind*(wrkT + 288.9414_rkind)  &
           /(508929.2_rkind*(wrkT+68.12963_rkind)))*1e-3_rkind
  
  end function
  
  function latent_heat_v(T) result(res)
    use typy
    real(kind=rkind), intent(in) :: T
    real(kind=rkind) :: wrkT, res
    wrkT = T - 273.15_rkind
    res = 2.501e6_rkind - 2369.2_rkind*T
  end function
  
  function spec_heat_solut(C, layer) result(res)
    use typy
    real(kind=rkind), intent(in) :: C
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind) :: wrkC, res
    wrkC = mol_conc(C,layer)
    res = 4.18e3_rkind*(0.730_rkind + 0.270_rkind*exp(-0.268_rkind*C))
  end function spec_heat_solut
  
  function diff_coef_NaCl(T) result(res)
    use typy
    real(kind=rkind), intent(in) :: T
    real(kind=rkind) ::  wrkT, res
    wrkT = T - 273.15_rkind
    res = 7.26e-10_rkind + 2.63e-11_rkind*wrkT + 2.18e-13_rkind*wrkT*wrkT
  end function diff_coef_NaCl
  
end module modRE_parameter_functions