!musim nutne zmenit retencni hrivky predevsim okraje ifu
! dodelat kontrolu pro x(3)

!> @file
!! Konstitutivní vztahy

!> @brief  Modul obsahuje konstitutivní vztahy
!!
!! Modul obsahuje konstituční vztahy pro Richardsovu rovnici. Jedná se o vyjádření funkce nenasycené hydraulické vodivosti (NHV) pro kapalnou i plynnou fázy půdního prostředí. Dále pak jsou uvedeny vztahy pro NHV obou dvou fází za vlivu teploty média a koncentrace NaCl rozpuštěném v médiu. Rovnice a jsou převzaty z článku \cite noborio
!! 
!! @author Jakub Jerabek
!! @version TEST

module modRE_constitutive

  public :: con_1_HustonAndCass
  public :: inv_con_1_HustonAndCass
  public :: capillar_capacity
  public :: con_2_HydraulicConductivity
  public :: con_3_HydraulicConductivityInvTemperature
  public :: con_4_HydraulicConductivityInvConcentration
  public :: con_5_HydraulicConductivityVapor
  public :: con_6_HydraulicConductivityVaporInvTemperature
  public :: con_7_HydraulicConductivityVaporInvConcentration
  
  public :: diff_coef
  
  public :: T_dirichlet_bc
  public :: T_neumann_bc !jakub dodelat pochopit
  public :: C_dirichlet_bc
  public :: C_neumann_bc !jakub dodelat pochopit
  
  public :: modRE_initcond
  public :: modHeat_initcond
  public :: modSolute_initcond
  
  private :: initcond

 contains
 
  
    
  
   !> @brief Potenciál vyjádřen retenční křivkou dle Huston And Cass \cite huston.
  !!
  !! Do výpočtu vstupuje 
  !! efektivni objemová vlhkost dopočtena z objemové vlhkosti dle rovnice (@f$\ref{sec:evol}@f$).
  !!
  !! \f[   
  !! \Psi = \left\{
  !! \begin{array}{l l} 
  !! -( \frac{1 - \theta_e}{0.01223} )^{0.5} & \quad \theta_e \geq 0.56362\\
  !!  -(0.551 \theta_e^{-2.1} + 4.136488)  & \quad \theta_e < 0.56362
  !! \end{array} \right.
  !! \label{eq:poten}
  !! \f]
  !!
  !! @param[in] layer identifikátor vrstvy (resp. jejich parametrů)
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out] potential \f$ \Psi \f$,  potenciál
  
  function con_1_HustonAndCass(layer, theta) result(h)
    use typy
    use modRE_globals
    use modRE_parameter_functions
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: theta
    real(kind=rkind) :: h 
    
    real(kind=rkind) :: the_i, the_s, the_r
    real(kind=rkind) :: a, b
    
    the_s = hyd_prop(layer)%theta_s
    the_r = hyd_prop(layer)%theta_r

    a = hyd_prop(layer)%hc_a
    b = hyd_prop(layer)%hc_b
    
    the_i = (2.0_rkind*the_s*b)/(1.0_rkind + 2.0_rkind*b)
    if (theta>the_s) then
      print*, "Warning: theta > theta_s !!!"
    else if (theta <= the_s .and. theta >= the_i) then
      h = a*(1.0_rkind-theta/the_s)**0.5_rkind*(the_i/the_s)**(-b) / (1.0_rkind-the_i/the_s)**0.5_rkind
    else if (theta < the_i .and. theta >= the_r) then
      h = a*(theta/the_s)**(-b)
    else if (theta < the_r) then
      print*, "Warning: theta < theta_r !!!"
    end if
 
  end function con_1_HustonAndCass
  
  function inv_con_1_HustonAndCass(layer, h) result(theta)
    use typy
    use modRE_globals
    use modRE_parameter_functions
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h
    real(kind=rkind) :: theta 
    real(kind=rkind) :: h_i, the_i, h_s, the_s, h_r, the_r
    real(kind=rkind) :: a, b
    
    the_s = hyd_prop(layer)%theta_s
    the_r = hyd_prop(layer)%theta_r

    a = hyd_prop(layer)%hc_a
    b = hyd_prop(layer)%hc_b
    
    the_i = (2.0_rkind*the_s*b)/(1.0_rkind + 2.0_rkind*b)
    h_i = a*(2*b/(1+2*b))**(-b)
    h_r = con_1_HustonAndCass(layer, the_r)
!     print*, "theta", theta


    if (h >= 0.0_rkind) then 
      theta = the_s
    else if (h < 0.0_rkind .and. h >= h_i) then
      theta = the_s - (the_s*h*h*(1.0_rkind-the_i/the_s))/(a*a*(the_i/the_s)**(-2.0_rkind*b))
    else if (h < h_i .and. h >= h_r) then
      theta = the_s*(h/a)**(-1.0_rkind/b)
    else if (h < h_r) then  !jakub muze to nastat? aby bylo h mensi nez theta r v ty pude?
      print*, "Warning: h < h_r !!!"
      theta = the_s*(h/a)**(-1.0_rkind/b)
!       theta = the_r
    else 

    end if

  end function inv_con_1_HustonAndCass

  
  function capillar_capacity(layer, h) result(res)
    use typy
    use modRE_globals
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h
    real(kind=rkind) :: res
    real(kind=rkind) :: h_i, the_i, h_s, the_s, h_r, the_r
    real(kind=rkind) :: a, b
    
    the_s = hyd_prop(layer)%theta_s
    the_r = hyd_prop(layer)%theta_r

    a = hyd_prop(layer)%hc_a
    b = hyd_prop(layer)%hc_b
    
    the_i = (2.0_rkind*the_s*b)/(1.0_rkind + 2.0_rkind*b)
    h_i = a*(2*b/(1+2*b))**(-b)
    h_r = con_1_HustonAndCass(layer, the_r)
  
  
    if (h >= 0.0_rkind) then 
      res = 0.0_rkind
    else if (h < 0.0_rkind .and. h >= h_i) then
      res = -2.0_rkind*h*(the_i/the_s)**(2.0_rkind*b)*(1.0_rkind-the_i/the_s)*the_s/a**2.0_rkind
    else if (h < h_i .and. h >= h_r) then
      res = -1.0_rkind*the_s/(b*h*(h/a)**(1.0_rkind/b))
    else if (h < h_r) then  !jakub muze to nastat? aby bylo h mensi nez theta r v ty pude?
      print*, "Warning: h < h_r !!!"
      res = -1.0_rkind*the_s/(b*h*(h/a)**(1.0_rkind/b))
!       theta = the_r
    else 

    end if
!     print*, "h", h, "capillar_capacity", res
  end function capillar_capacity
  
  
  
  
  !> @brief Nenasycená hydraulická vodivost kapalné fáze.
  !!
  !! Do výpočtu vstupuje
  !! gain factor @f$ G_{K_{\Psi}} @f$ (@f$\ref{par:gain}@f$), 
  !! počáteční kinematická viskosita @f$ \nu_0 @f$ (@f$\ref{par:vis}@f$) a 
  !! kinematická viskosita závislá na teplotě a koncentraci @f$ \nu @f$ (@f$\ref{eq:vis}@f$), 
  !! dále pak nasycenná hydraulická vodivost @f$ K^0_S @f$ (@f$\ref{par:K0}@f$) a 
  !! efektivkí vlhkost @f$ \theta_e @f$ (@f$\ref{sec:evol}@f$).
  !!
  !! @f[   
  !! K_{l, \Psi} =  \left\{ \begin{array}{l l}
  !! (1-G_{K_{\Psi}}(1-\nu^0 / \nu)) K^0_S \theta_e^{11.39}  &for \quad 0.94558 < \theta_e \leq 0.56362   \\
  !! (1-G_{K_{\Psi}}(1-\nu^0 / \nu)) 0.71356 K^0_S \theta_e^{5.3817}  &for  \quad  0.59766 < \theta_e \leq 0.94558 \\
  !! (1-G_{K_{\Psi}}(1-\nu^0 / \nu)) 5.5444 K^0_S \theta_e^{8.9775}  &for \quad  0.0 < \theta_e \leq 0.59766 
  !! \end{array} \right.
  !! \label{eq:NHV}
  !! @f]
  !! @param[in] layer layer identifikátor
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out] hyd_conduc \f$ K_{l, \Psi} \f$,  nenasycená hydraulická vodivost

  function con_2_HydraulicConductivity(layer, h, t, c) result(hyd_conduc)
    use typy
    use pde_objs
!     use global_objs
    use globals
    use modRE_globals
    use modRE_parameter_functions    
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), allocatable, dimension(:,:) :: hyd_conduc
    
    real(kind=rkind) :: Gky
    real(kind=rkind) :: vis_0
    real(kind=rkind), dimension(3,3) :: Ks
    real(kind=rkind) :: theta, thete
    real(kind=rkind) :: vis
    integer :: d
    
    d = drutes_config%dimen
    allocate(hyd_conduc(1:d,1:d))
    theta = inv_con_1_HustonAndCass(layer, h)
    thete = evol(layer, theta)
    vis_0 = hyd_prop(layer)%viscosity_0
    Gky =  hyd_prop(layer)%gain_factor !jakub ten gain factor musim taky nejak posasit
    Ks(1:d,1:d) =  hyd_prop(layer)%Ks(1:d,1:d)
    vis = viscosity(t, c, layer)

    if (thete < 1.0_rkind .and. thete >= 0.94558_rkind) then
      hyd_conduc = (1.0_rkind - Gky*(1.0_rkind-vis_0/vis))*Ks(1:d,1:d)*thete**(11.39_rkind)
    else if (thete < 0.94558_rkind .and. thete >= 0.59766_rkind) then
      hyd_conduc = (1.0_rkind - Gky*(1.0_rkind-vis_0/vis))*0.71356_rkind*Ks(1:d,1:d)*thete**(5.3817_rkind)
    else if (thete <  0.59766_rkind .and. thete > 0.0_rkind) then
      hyd_conduc = (1.0_rkind - Gky*(1.0_rkind-vis_0/vis))*4.54444_rkind*Ks(1:d,1:d)*thete**(8.9775_rkind)
    else 
      hyd_conduc = KS(1:d,1:d)
    end if

  end function con_2_HydraulicConductivity

  function der_con_2_HydraulicConductivity(layer, h, t, c) result(dhyd_conduc)
    use typy
    use pde_objs
!     use global_objs
    use globals
    use modRE_globals
    use modRE_parameter_functions    
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), allocatable, dimension(:,:) ::  dhyd_conduc
    
    real(kind=rkind) :: Gky
    real(kind=rkind) :: vis_0
    real(kind=rkind), dimension(3,3) :: Ks
    real(kind=rkind) :: theta, thete
    real(kind=rkind) :: vis
    integer :: d
    
    d = drutes_config%dimen
    allocate(dhyd_conduc(1:d,1:d))
    theta = inv_con_1_HustonAndCass(layer, h)
    thete = evol(layer, theta)
    vis_0 = hyd_prop(layer)%viscosity_0
    Gky =  hyd_prop(layer)%gain_factor !jakub ten gain factor musim taky nejak posasit
    Ks(1:d,1:d) =  hyd_prop(layer)%Ks(1:d,1:d)
    vis = viscosity(t, c, layer)

    if (thete <= 1.0_rkind .and. thete >= 0.94558_rkind) then
      dhyd_conduc = 11.39_rkind*Ks(1:d,1:d)*thete**10.39_rkind*(1.0_rkind-Gky*(1.0_rkind-vis_0/vis))
    else if (thete < 0.94558_rkind .and. thete >= 0.59766_rkind) then
      dhyd_conduc = 3.840165852_rkind*Ks(1:d,1:d)*thete**4.3817_rkind*(1.0_rkind-Gky*(1.0_rkind-vis_0/vis))
    else if (thete <  0.59766_rkind .and. thete > 0.0_rkind) then
      dhyd_conduc = 40.79771009999999_rkind*Ks(1:d,1:d)*thete**7.977499999999999_rkind*(1.0_rkind-Gky*(1.0_rkind-vis_0/vis))
    else 
      dhyd_conduc = 0.0_rkind
    end if
  end function der_con_2_HydraulicConductivity
  
  
  
  
  
  !> @brief Nenasycená hydraulická vodivost kapalné fáze s vlivem teploty.
  !!
  !! Do výpočtu vstupuje 
  !! nenasycená hydraulická vodivost @f$K_{l, \Psi}@f$ (@f$\ref{eq:NHV}@f$), 
  !! potenciál z retenční křivky @f$\Psi@f$ (@f$\ref{eq:poten}@f$), 
  !! gain factor @f$G_{\Psi, T}@f$ (@f$\ref{par:gainT}@f$), 
  !! povrchové napětí vody s obsahem NaCl @f$\gamma_0@f$ (@f$\ref{eq:gamma}@f$), 
  !! derivace povrchového napětí vody s obsahem NaCl dle teploty @f$\frac{\partial \gamma}{\partial T}@f$ (@f$\ref{eq:gammaT}@f$), 
  !! koeficient vlivu osmotického efektu @f$\sigma@f$ (@f$\ref{eq:osmotic}@f$), 
  !! derivace osmotického potenciálu dle teploty @f$\frac{\partial \Psi_\pi}{\partial T}@f$ (@f$\ref{eq:der_osmo}@f$) (viz níže).
  !! 
  !! @f[
  !!	K_{l, T} = K_{l, \Psi} \left(\Psi G_{\Psi, T} \frac{1}{\gamma^0} \frac{\partial \gamma}{\partial T} + \sigma \frac{\partial \Psi_\pi}{\partial T}\right)
  !!	\label{eq:KlT}
  !! @f]
  !! 
  !! Rovnice pro derivace osmotického potenciálu dle teploty, kde do rovnice vstupuje 
  !! plynová konstanta @f$R@f$ (@f$\ref{par:gasCon}@f$), 
  !! molární hmotnost vody @f$M_w@f$ (@f$\ref{par:molW}@f$) a 
  !! "aktivita" vody @f$h_\pi@f$ (@f$\ref{eq:activityW}@f$).
  !! 
  !! @f[
  !! \frac{\partial \Psi_\pi}{\partial T} = \frac{R}{M_W} \ln{(h_\pi)}
  !! \label{eq:der_osmo}
  !! @f]
  !! 
  !! Člen @f$\frac{1}{\gamma^0} \frac{\partial \gamma}{\partial T}@f$ v rovnici (@f$\ref{eq:KlT}@f$) vyjadřije (podle \cite nassar) tzv.: koeficient teploty ovlivňující povrchová napětí. V tomto případě je @f$\gamma@f$ vyjádřeno rovnicí (@f$\ref{eq:gamma}@f$) a člen @f$\frac{\partial \gamma}{\partial T}@f$ je derivace této rovnice podle teploty (@f$\ref{eq:gammaT}@f$). 
  !! 
  !! @param[in] layer layer identifikátor
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out] hyd_conduc_inv_T \f$  K_{l,T} \f$,  nenasycená hydraulická vodivost s vlivem teploty

  function con_3_HydraulicConductivityInvTemperature(layer, h, t, c) result(hyd_conduc_inv_T)
    use typy
    use pde_objs
    use modRE_globals
    use modRE_parameter_functions    

    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), allocatable, dimension(:,:) :: hyd_conduc_inv_T

    real(kind=rkind), dimension(3,3) :: hyd_conduc
!     real(kind=rkind) :: potential
    real(kind=rkind) :: Gt
    real(kind=rkind) :: gamma
    real(kind=rkind) :: gamma_derT
    real(kind=rkind) :: sigma
    real(kind=rkind) :: osmotic_pot_der
    real(kind=rkind) :: R, Mw, hpi
    integer :: d
    
    d = drutes_config%dimen
    allocate(hyd_conduc_inv_T(1:d,1:d))
    hyd_conduc(1:d,1:d) = con_2_HydraulicConductivity(layer, h, t, c)

    Gt = hyd_prop(layer)%gain_factor_T
    R = hyd_prop(layer)%gas_constant 
    Mw = hyd_prop(layer)%molecular_weight

    sigma = osmotic_coef(layer)
    
    gamma = surface_tension_polutant(t, c, layer)
    gamma_derT = derT_surface_tension_polutant(t)
    hpi = water_activity_polutant(c, layer)
    osmotic_pot_der = (R/Mw)*log(hpi)     ! log je ve fortranu prirozeny logaritmus, dekadicy je log10

    hyd_conduc_inv_T = hyd_conduc(1:d,1:d)*(h*Gt*(1.0_rkind/gamma)*gamma_derT + sigma*osmotic_pot_der)

  end function con_3_HydraulicConductivityInvTemperature

  !> @brief Nenasycená hydraulická vodivost kapalné fáze s vlivem koncentrace.
  !!
  !! Do výpočtu vstupuje 
  !! nenasycená hydraulická vodivost @f$K_{l, \Psi}@f$ (@f$\ref{eq:NHV}@f$), 
  !! potenciál z retenční křivky @f$\Psi@f$ (@f$\ref{eq:poten}@f$), 
  !! povrchové napětí vody s obsahem NaCl @f$\gamma_0@f$ (@f$\ref{eq:gamma}@f$), 
  !! derivace povrchového napětí vody s obsahem NaCl dle koncentrace @f$\frac{\partial \gamma}{\partial C}@f$ (@f$\ref{eq:gammaC}@f$), 
  !! koeficient vlivu osmotických sil @f$\sigma@f$ (@f$\ref{eq:osmotic}@f$), 
  !! derivace osmotického potenciálu dle koncentrace @f$\frac{\partial \Psi_\pi}{\partial C}@f$ (@f$\ref{eq:der_osmoC}@f$) (viz níže).
  !!
  !! @f[
  !!	K_{l,C} = K_{l,\Psi} \left(\Psi G_{\Psi, T} \frac{1}{\gamma^0} \frac{\partial \gamma}{\partial C} + \sigma \frac{\partial \Psi_\pi}{\partial C}\right)
  !!	\label{eq:KlC}
  !! @f]
  !! 
  !! Rovnice pro derivaci osmotického potenciálu dle koncentrace, kde do rovnice vstupuje 
  !! plynová konstanta @f$R@f$ (@f$\ref{par:gasCon}@f$), 
  !! molární hmotnost vody @f$M_w@f$ (@f$\ref{par:molW}@f$), 
  !! teplota @f$T@f$ (@f$\ref{par:tem}@f$) a "aktivita" vody @f$h_\pi@f$ (@f$\ref{eq:activityW}@f$).
  !!  
  !! @f[
  !! \frac{\partial \Psi_\pi}{\partial C} = \frac{\ln{(h_\pi)}}{C} \frac{RT}{M_W}
  !! \label{eq:der_osmoC} 
  !! @f]
  !! 
  !! Rovnice (@f$\ref{eq:KlC}@f$) je podle \cite noborio v podstatě ronice (@f$\ref{eq:KlT}@f$), pouze je povrchové napětí @f$\gamma@f$ a osmotický potenciál @f$\Psi_\pi@f$ derivováno místo podle teploty podle koncentrace. Zde je PROBLEM, pokud zderivuju povrchové napětí @f$\gamma@f$ (rovnici (@f$\ref{eq:gamma}@f$)) podle C je výsledek s opačným znamínkem než deriva @f$\gamma@f$ podle @f$T@f$ a rovnice (@f$\ref{eq:KlC}@f$) vyjde záporně. Což je BLBOST!!
  !! 
  !! @param[in] layer layer identifikátor
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out] hyd_conduc_inv_C \f$ K_{l,C} \f$,  nenasycená hydraulická vodivost s vlivem koncentrace


  function con_4_HydraulicConductivityInvConcentration(layer, h, t, c) result(hyd_conduc_inv_C)
    use typy
    use pde_objs
    use modRE_globals
    use modRE_parameter_functions
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), allocatable, dimension(:,:) :: hyd_conduc_inv_C

!     real(kind=rkind) :: potential
    real(kind=rkind), dimension(3,3) :: hyd_conduc
    real(kind=rkind) :: gamma
    real(kind=rkind) :: gamma_derC
    real(kind=rkind) :: sigma
    real(kind=rkind) :: osmotic_pot_der
    real(kind=rkind) :: c_mol
    real(kind=rkind) :: R, Mw, hpi
    integer :: d
    d = drutes_config%dimen
    allocate(hyd_conduc_inv_C(1:d,1:d))
    hyd_conduc(1:d,1:d) = con_2_HydraulicConductivity(layer, h, t, c)
    
    R = hyd_prop(layer)%gas_constant 
    Mw = hyd_prop(layer)%molecular_weight

    sigma = osmotic_coef(layer)

    gamma = surface_tension_polutant(t, c, layer)
    gamma_derC =derC_surface_tension_polutant(c)
    hpi =  water_activity_polutant(c, layer)
    c_mol = mol_conc(c, layer)
    osmotic_pot_der = (log(hpi)/c_mol)*((R*t)/Mw)

    hyd_conduc_inv_C = hyd_conduc(1:d,1:d)*(h*(1.0_rkind/gamma)*gamma_derC + sigma * osmotic_pot_der)

  end function con_4_HydraulicConductivityInvConcentration
  
  
  !> @brief Nenasycená hydraulická vodivost plynné fáze
  !!
  !! Do rovnice vstupuje 
  !! koeficient difuze pudního vzduchu @f$D_g@f$ (@f$\ref{eq:bin_gas}@f$), 
  !! teplota @f$T@f$ (@f$\ref{par:tem}@f$), 
  !! poměr pórů naplněných vzduchem  @f$\phi_g@f$ (dopočteno z plně nasycene vlhkosti @f$\theta_s@f$ (@f$\ref{par:the_s}@f$) a aktualni objemové vlhkosti @f$\theta@f$; @f$\phi_g = \theta_s - \theta@f$), 
  !! relativní vlhkost vzduchu @f$h_m@f$ (@f$\ref{eq:relat_hum}@f$),  
  !! "aktivita" vody @f$h_\pi@f$ (@f$\ref{eq:activityW}@f$), 
  !! hustota půdního vzduchu @f$\rho^{sat}_v@f$ (@f$\ref{eq:vap_dens}@f$), 
  !! molární hmotnost vody @f$M_w@f$ (@f$\ref{par:molW}@f$), 
  !! plynová konstanta  @f$R@f$ (@f$\ref{par:gasCon}@f$).
  !!
  !! @f[
  !!	K_{v,\Psi} = D_g \left(\frac{T}{273.15}\right)^2 \phi^{5/3}_g h_m h_\pi \rho^{sat}_v \frac{M_w}{RT}
  !! @f]
  !! 
  !! @param[in] layer layer identifikátor
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out] hyd_conduc_vap \f$  K_{v,\Psi} \f$,  nenasycená hydraulická vodivost plynné fáze 

  function con_5_HydraulicConductivityVapor(layer, h, t, c) result(hyd_conduc_vap)
    use typy
    use pde_objs
    use modRE_globals
    use modRE_parameter_functions
    use debug_tools
    use readtools
    
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h
    real(kind=rkind), allocatable, dimension(:,:) :: hyd_conduc_vap

    real(kind=rkind) :: Dg
    real(kind=rkind) :: t
    real(kind=rkind) :: fi_air
    real(kind=rkind) :: hm
    real(kind=rkind) :: hpi
    real(kind=rkind) :: rho_sat
    real(kind=rkind) :: R
    real(kind=rkind) :: Mw
    real(kind=rkind), dimension(3) :: K_Psi_v_loc

    real(kind=rkind) :: c
    real(kind=rkind) :: theta
    integer :: d
    
    d = drutes_config%dimen
    allocate(hyd_conduc_vap(1:d,1:d))
    theta = inv_con_1_HustonAndCass(layer,h)
    
    Mw = hyd_prop(layer)%molecular_weight
    R = hyd_prop(layer)%gas_constant 
    Dg = molecular_D_vapor_water(t)
    hm = relativ_humidity_soil(theta, R, Mw, t)
    rho_sat = saturated_vapor_density(t)
    fi_air = hyd_prop(layer)%theta_s - theta
    hpi = water_activity_polutant(c, layer)

    K_Psi_v_loc(1:d) = Dg*(t/273.15_rkind)**2.0_rkind * fi_air**(5.0_rkind/3.0_rkind)  * &
				  hm*hpi*rho_sat*(Mw/(R*t))
!
    call set_tensor(K_Psi_v_loc(1:d), hyd_prop(layer)%anisoangle(:),hyd_conduc_vap)
				  
				      
!     print*, "h= ", h, "fi_air = ", fi_air, "theta", theta, "theta_s = ", the_s
  end function con_5_HydraulicConductivityVapor
  
  !> @brief Nenasycená hydraulická vodivost plynné fáze s vlivem teploty.
  !!
  !! Do rovnice vstupuje 
  !! koeficient difuze půdního vzduchu @f$D_g@f$ (@f$\ref{eq:bin_gas}@f$), 
  !! teplota @f$T@f$ (@f$\ref{par:tem}@f$), 
  !! poměr pórů naplněných vzduchem  @f$\phi_g@f$(dopočteno z plně nasycene vlhkosti @f$\theta_s@f$ (@f$\ref{par:the_s}@f$) a aktualni objemové vlhkosti @f$\theta@f$; @f$\phi_g = \theta_s - \theta@f$), 
  !! relativní vlhkosti vzduchu @f$h_m@f$ (@f$\ref{eq:relat_hum}@f$),  
  !! "aktivita" vody @f$h_\pi@f$ (@f$\ref{eq:activityW}@f$), 
  !! derivace hustoty půdního vzruchu podle teploty @f$\frac{\mathrm{d}\rho^{sat}_v}{\mathrm{d}T}@f$ (@f$\ref{eq:der_vap_dens}@f$), 
  !! faktor mechanického zatížení @f$\eta@f$ (@f$\ref{eq:mech_enha}@f$).
  !!
  !! @f[
  !!	K_{v,T} = D_g\left(\frac{T}{273.15} \right)^2 \phi^{5/3}_g h_m h_\pi \frac{\mathrm{d} \rho^{sat}_v}{\mathrm{d} T} \eta
  !! @f]
  !! 
  !! 
  !! @param[in] layer layer identifikátor
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out]  hyd_conduc_vap_inv_T \f$ K_{v,T} \f$,  nenasycená hydraulická vodivost plynné fáze s vlivem teploty


  function con_6_HydraulicConductivityVaporInvTemperature(layer, h, t, c) result(hyd_conduc_vap_inv_T)
    use typy
    use pde_objs
    use modRE_globals
    use modRE_parameter_functions
    use readtools
    
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), allocatable, dimension(:,:) :: hyd_conduc_vap_inv_T

    real(kind=rkind) :: Dg
    real(kind=rkind) :: fi_air
    real(kind=rkind) :: hm
    real(kind=rkind) :: hpi
    real(kind=rkind) :: der_rho_sat
    real(kind=rkind) :: ny
    real(kind=rkind), dimension(3) :: K_T_v_loc
!     real(kind=rkind) :: potential

    real(kind=rkind) :: R, Mw
    real(kind=rkind) :: clay_frac
    real(kind=rkind) :: theta
    integer :: d
    
    d = drutes_config%dimen
    allocate(hyd_conduc_vap_inv_T(1:d,1:d))
    
    theta = inv_con_1_HustonAndCass(layer, h)

    Mw = hyd_prop(layer)%molecular_weight
    R = hyd_prop(layer)%gas_constant 
    Dg = molecular_D_vapor_water(t)
    clay_frac = hyd_prop(layer)%clay_fraction




    der_rho_sat = der_saturated_vapor_density(t)
    hm = relativ_humidity_soil(theta, R, Mw, t)
    ny = mechanics_enhancement_factor(theta, clay_frac)
    hpi = water_activity_polutant(c, layer)
    fi_air = hyd_prop(layer)%theta_s - theta



    K_T_v_loc(1:d) = Dg*(t/273.15_rkind)**2.0_rkind * fi_air**(5.0_rkind/3.0_rkind)  * &
				hm*hpi*der_rho_sat*ny
    !
    call set_tensor(K_T_v_loc(1:d), hyd_prop(layer)%anisoangle(:),hyd_conduc_vap_inv_T)

  end function con_6_HydraulicConductivityVaporInvTemperature
  
  !> @brief Nenasycená hydraulická vodivost plynné fáze s vlivem koncentrace.
  !!
  !! Do rovnice vstupuje 
  !! koeficient difuze půdního vzduchu @f$D_g@f$ (@f$\ref{eq:bin_gas}@f$), 
  !! teplota @f$T@f$ (@f$\ref{par:tem}@f$), 
  !! poměr pórů naplněných vzduchem  @f$\phi_g@f$(dopočteno z plně nasycene vlhkosti @f$\theta_s@f$ (@f$\ref{par:the_s}@f$) a aktualni objemové vlhkosti @f$\theta@f$; @f$\phi_g = \theta_s - \theta@f$), 
  !! relativní vlhkosti vzduchu @f$h_m@f$ (@f$\ref{eq:relat_hum}@f$),  
  !! hustota půdního vzduchu @f$\rho^{sat}_v@f$ (@f$\ref{eq:vap_dens}@f$), 
  !! derivace "aktivity" vody podle koncentrace (viz níže).
  !!
  !! @f[
  !!	K_{v,C} = D_g\left(\frac{T}{273.15} \right)^2 \phi^{5/3}_g h_m \rho^{sat}_v  \frac{\mathrm{d} h_\pi}{\mathrm{d} C}
  !! @f]
  !! 
  !! Do rovnice derivace "aktivity" vody podle koncentrace vstupuje 
  !! "aktivita" vody @f$h_\pi@f$ (@f$\ref{eq:activityW}@f$), 
  !! koncentrace @f$C@f$ (@f$\ref{par:c}@f$).
  !! 
  !! @f[
  !! \frac{\mathrm{d} h_\pi}{\mathrm{d} C} = \frac{\ln{(h_\pi)}}{C} h_\pi
  !! @f]
  !! 
  !! @param[in] layer layer identifikátor
  !! @param[in] theta \f$ \theta \f$, objemová vlhkost
  !! @param[out] hyd_conduc_vap_inv_C \f$  K_{v,C} \f$,  nenasycená hydraulická vodivost plynné fáze s vlivem koncentrace


  function con_7_HydraulicConductivityVaporInvConcentration(layer, h, t, c) result(hyd_conduc_vap_inv_C)
    use typy
    use pde_objs
    use modRE_globals
    use modRE_parameter_functions
    use readtools
  
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), allocatable, dimension(:,:)  :: hyd_conduc_vap_inv_C

    real(kind=rkind) :: Dg
    real(kind=rkind) :: fi_air
    real(kind=rkind) :: hm
    real(kind=rkind) :: rho_sat
    real(kind=rkind) :: der_hpi
    real(kind=rkind), dimension(3) :: K_C_v_loc
    
    

    real(kind=rkind) :: R, Mw
    real(kind=rkind) :: hpi
    real(kind=rkind) :: theta
    integer :: d
    
    d = drutes_config%dimen
    allocate(hyd_conduc_vap_inv_C(1:d,1:d))
    
    theta = inv_con_1_HustonAndCass(layer, h)
    
    Mw = hyd_prop(layer)%molecular_weight
    R = hyd_prop(layer)%gas_constant 
    Dg = molecular_D_vapor_water(t)

    fi_air = hyd_prop(layer)%theta_s - theta
    hm = relativ_humidity_soil(theta, R, Mw, t)
    rho_sat = saturated_vapor_density(t)
    hpi = water_activity_polutant(c, layer)
    der_hpi = der_water_activity_polutant(c, layer)
!     der_hpi = 1.0_rkind - 3.1327e-2_rkind - 1.4553e-3_rkind*c_mol*2.0_rkind


    K_C_v_loc(1:d) = Dg*(t/273.15_rkind)**2.0_rkind * fi_air**(5.0_rkind/3.0_rkind)  * &
				  hm*rho_sat*der_hpi
				  
    call set_tensor(K_C_v_loc(1:d), hyd_prop(layer)%anisoangle(:),hyd_conduc_vap_inv_C)
  end function con_7_HydraulicConductivityVaporInvConcentration
  
!   
  function diff_coef(layer, h, t, c, gradh, gradt, gradc) result(res)
    use typy
    use globals
    use pde_objs
    use modRE_globals
    use modRE_parameter_functions
    use readtools
    
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), intent(in) :: h, t, c
    real(kind=rkind), dimension(:), intent(in) :: gradh, gradt, gradc
    real(kind=rkind), dimension(3,3) :: res
    real(kind=rkind), dimension(3) :: rec_loc

    integer :: i, d
    real(kind=rkind), dimension(3,3) :: K_Psi_l, K_T_l, K_C_l
    real(kind=rkind), dimension(3) :: vct
    real(kind=rkind) :: theta
    real(kind=rkind) :: flux_length
    real(kind=rkind) :: D0, kappa, rho_l
    
    d = drutes_config%dimen
    D0 = diff_coef_NaCl(t)
    theta = inv_con_1_HustonAndCass(layer, h)
    rho_l = hyd_prop(layer)%density_soil_water
    kappa = hyd_prop(layer)%kappa 
    
    K_Psi_l(1:d,1:d) = con_2_HydraulicConductivity(layer, h, t, c)
    K_T_l(1:d,1:d) = con_3_HydraulicConductivityInvTemperature(layer, h, t, c)
    K_C_l(1:d,1:d) = con_4_HydraulicConductivityInvConcentration(layer, h, t, c)
    
    vct = 0.0_rkind
    vct(d) = norm2(K_Psi_l(d,1:d))  !jakub vliv gravitace je ok?
    vct(1:d) = vct(1:d)  + matmul(K_Psi_l(1:d,1:d),gradh(1:d))
    vct(1:d) = vct(1:d)  + matmul(K_T_l(1:d,1:d), gradt(1:d))
    vct(1:d) = vct(1:d)  + matmul(K_C_l(1:d,1:d), gradc(1:d))
    
    flux_length = norm2(vct(1:d))
    
    rec_loc = 2.8_rkind*D0*theta**3.0_rkind + kappa*flux_length/(rho_l * theta)
    call set_tensor(rec_loc, hyd_prop(layer)%anisoangle, res)
    
  end function diff_coef
  
  
  subroutine T_dirichlet_bc(el_id, node_order, value, code) 
    use typy
    use globals
    use global_objs
    use pde_objs
    integer(kind=ikind), intent(in)  :: el_id, node_order
    real(kind=rkind), intent(out), optional    :: value
    integer(kind=ikind), intent(out), optional :: code
    integer(kind=ikind) :: edge_id, i, j
    
    if (present(value)) then
      edge_id = nodes%edge(elements%data(el_id, node_order))

      if (pde(2)%bc(edge_id)%file) then
	do i=1, ubound(pde(2)%bc(edge_id)%series,1)
	  if (pde(2)%bc(edge_id)%series(i,1) > time) then
	    if (i > 1) then
	      j = i-1
	    else
	      j = i
	    end if
	    value = pde(2)%bc(edge_id)%series(j,2)
	    EXIT
	  end if
	end do
      else
	value = pde(2)%bc(edge_id)%value
      end if
    end if

    if (present(code)) code = 1

  end subroutine T_dirichlet_bc
  
  
  subroutine T_neumann_bc(el_id, node_order, value, code)
    use typy 
    use globals
    use global_objs
    use pde_objs
    
    integer(kind=ikind), intent(in)  :: el_id, node_order
    real(kind=rkind), intent(out), optional    :: value
    integer(kind=ikind), intent(out), optional :: code
    
    real(kind=rkind), dimension(3,3) :: K
    
    integer(kind=ikind) :: i, edge_id, j, d
    real(kind=rkind) :: tmp
    real(kind=rkind), dimension(3) :: bcflux
    integer :: i1
    
    if (present(value)) then
      d = drutes_config%dimen
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      i = pde(2)%permut(elements%data(el_id, node_order))
      
      call pde(2)%pde_fnc(2)%dispersion(pde(2), elements%material(el_id,1),  &
					x=(/pde_common%xvect(i,2)/), tensor=K(1:d, 1:d))
      !jakub prej je tu neco blbe, grav flux tu asi nebude
  !     gravflux(1:d) = K(d, 1:d)*elements%nvect_z(el_id, node_order)
    
    
      if (pde(2)%bc(edge_id)%file) then
	do i=1, ubound(pde(2)%bc(edge_id)%series,1)
	  if (pde(2)%bc(edge_id)%series(i,1) > time) then
	    if (i > 1) then
	      j = i-1
	    else
	      j = i
	    end if
	      tmp = pde(2)%bc(edge_id)%series(j,2)
	      EXIT
	  end if
	end do
      else
	tmp = pde(2)%bc(edge_id)%value
      end if
    
      if (tmp >= 0) then
	i1 = 1
      else
	i1 = -1
      end if
      
      
      select case(d)
	case(1)
	  if (edge_id == 101) then
	    value = tmp! + gravflux(1)
	  else
	    value = tmp! - gravflux(1)
	  end if
	  RETURN
	case(2)
	  bcflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*tmp
	  bcflux(2) = elements%nvect_z(el_id, node_order)*tmp
  ! 	bcflux = bcflux - gravflux
	  value = i1*sqrt(bcflux(1)*bcflux(1) + bcflux(2)*bcflux(2))
      end select
    end if
    
    if (present(code)) code = 2

    
  end subroutine T_neumann_bc
  
  subroutine C_dirichlet_bc(el_id, node_order, value, code) 
    use typy
    use globals
    use global_objs
    use pde_objs
    integer(kind=ikind), intent(in)  :: el_id, node_order
    real(kind=rkind), intent(out), optional    :: value
    integer(kind=ikind), intent(out), optional :: code
    integer(kind=ikind) :: edge_id, i, j
    
    
    if (present(value)) then
      edge_id = nodes%edge(elements%data(el_id, node_order))

      if (pde(3)%bc(edge_id)%file) then
	do i=1, ubound(pde(3)%bc(edge_id)%series,1)
	  if (pde(3)%bc(edge_id)%series(i,1) > time) then
	    if (i > 1) then
	      j = i-1
	    else
	      j = i
	    end if
	    value = pde(3)%bc(edge_id)%series(j,2)
	    EXIT
	  end if
	end do
      else
	value = pde(3)%bc(edge_id)%value
      end if
    end if

    if (present(code)) code = 1
    
  end subroutine C_dirichlet_bc
  
  subroutine C_neumann_bc(el_id, node_order, value, code)
    use typy 
    use globals
    use global_objs
    use pde_objs
    
    integer(kind=ikind), intent(in)  :: el_id, node_order
    real(kind=rkind), intent(out), optional    :: value
    integer(kind=ikind), intent(out), optional :: code
    
    real(kind=rkind), dimension(3,3) :: K
    
    integer(kind=ikind) :: i, edge_id, j, d
    real(kind=rkind) :: tmp
    real(kind=rkind), dimension(3) ::  bcflux
    integer :: i1
    
    if (present(value)) then
      d = drutes_config%dimen
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      i = pde(3)%permut(elements%data(el_id, node_order))
      
      call pde(3)%pde_fnc(3)%dispersion(pde(3), elements%material(el_id,1),  &
					x=(/pde_common%xvect(i,2)/), tensor=K(1:d, 1:d))
      !jakub prej je tu neco blbe, grav flux tu asi nebude
      !gravflux(1:d) = K(d, 1:d)*elements%nvect_z(el_id, node_order)
    
    
      if (pde(3)%bc(edge_id)%file) then
	do i=1, ubound(pde(3)%bc(edge_id)%series,1)
	  if (pde(3)%bc(edge_id)%series(i,1) > time) then
	    if (i > 1) then
	      j = i-1
	    else
	      j = i
	    end if
	      tmp = pde(3)%bc(edge_id)%series(j,2)
	      EXIT
	  end if
	end do
      else
	tmp = pde(3)%bc(edge_id)%value
      end if
    
      if (tmp >= 0) then
      i1 = 1
      else
	i1 = -1
      end if
      
      
      select case(d)
	case(1)
	  if (edge_id == 101) then
	    value = tmp! + gravflux(1)
	  else
	    value = tmp! - gravflux(1)
	  end if
	  RETURN
	case(2)
	  bcflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*tmp
	  bcflux(2) = elements%nvect_z(el_id, node_order)*tmp
  ! 	bcflux = bcflux - gravflux
	  value = i1*sqrt(bcflux(1)*bcflux(1) + bcflux(2)*bcflux(2))
      end select
    end if
    
    if (present(code)) code = 2

    
  end subroutine C_neumann_bc
  
  subroutine modRE_initcond(pde_loc)
    use typy
    use globals
    use global_objs
    use pde_objs
!     use modRE_globals
    class(pde_str), intent(in out) :: pde_loc
    call initcond(pde_loc, 1_ikind)
  end subroutine modRE_initcond
  
  subroutine modHeat_initcond(pde_loc)
    use typy
    use globals
    use global_objs
    use pde_objs
!     use modRE_globals
    class(pde_str), intent(in out) :: pde_loc
    call initcond(pde_loc, 2_ikind)
  end subroutine modHeat_initcond
  
  subroutine modSolute_initcond(pde_loc)
    use typy
    use globals
    use global_objs
    use pde_objs
!     use modRE_globals
    class(pde_str), intent(in out) :: pde_loc
    call initcond(pde_loc, 3_ikind)
  end subroutine modSolute_initcond
  
  subroutine initcond(pde_loc, process) 
        use typy
        use globals
        use pde_objs
        use global_objs
        use modRE_globals
        
        class(pde_str), intent(in out) :: pde_loc
        integer(kind=ikind), intent(in) :: process
        integer(kind=ikind) :: i, j, k,l, m, layer
	
	  do i=1, elements%kolik
	    layer = elements%material(i,1)
	    do j=1, ubound(elements%data,2)
	      k = elements%data(i,j)
	      l = nodes%edge(k)
	      m = pde(process)%permut(k)
	      if (m == 0) then
		pde_loc%solution(k) =  pde(process)%bc(l)%value
	      else
		pde_loc%solution(k) = hyd_prop(layer)%initcond(process)
	      end if
	    end do   
	  end do
  end subroutine initcond
  
  
  subroutine STVF()
   use typy
   
  
  
  
  end subroutine
  
  

end module modRE_constitutive

!   subroutine dispersion_to_potential_vapor(layer, theta, tensor, scalar)
!     use typy 
!     use modRE_globals    
!     integer(kind=ikind), intent(in) :: layer
!     real(kind=rkind), intent(in) :: theta
!     real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
!     real(kind=rkind), intent(out), optional :: scalar        
!     real(kind=rkind) :: tmp
!     
!     tmp  = con_5_HydraulicConductivityVapor(layer, theta)
!     
!     if (present(scalar)) then
!       scalar = tmp
!     end if 
!   
!     if (present(tensor)) then
!       tensor = tmp*hyd_prop(layer)%Ks
!     end if 
!   
!   end subroutine dispersion_to_potential_vapor
  