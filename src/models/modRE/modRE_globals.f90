!> @file
!! Deklarace globalních proměnných a datových typů

!> @brief  Modul obsahuje deklarace a definice globalnich proměnných  
!!
!! Modul obsahuje deklarace a definice globalnich proměnných pro modifikovanou Richardsovu rovnici dle. \cite noborio.
!! @author Jakub Jerabek
!! @version TEST

module modRE_globals
  use typy
  use global_objs
!   integer (kind=ikind), public :: selector

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
!   !> Deklarace pointru na základní konstitutivní vztahy.
!   !!
!   !! Funkce s danou "sablonou"  scalar_fnc (z global_objs).
!   !!
!   type, public :: coupled_fnc
!     procedure(scalar_fnc), nopass, pointer :: potential !< Potenciál vyjádřen retenční křivkou dle \cite noborio.
!     procedure(scalar_fnc), nopass, pointer :: hydraulic_conduc !< Hydraulická vodivost kapalné fáze dle \cite noborio.
!     procedure(scalar_fnc), nopass, pointer :: hydraulic_conduc_invT !< Hydraulická vodivost kapalné fáze s vlivem teploty dle \cite noborio.
!     procedure(scalar_fnc), nopass, pointer :: hydraulic_conduc_invC !< Hydraulická vodivost kapalné fáze s vlivem koncentrace dle \cite noborio.
!     procedure(scalar_fnc), nopass, pointer :: hydraulic_conduc_vapor !< Hydraulická vodivost plynné fáze dle \cite noborio.
!     procedure(scalar_fnc), nopass, pointer :: hydraulic_conduc_vapor_invT !< Hydraulická vodivost plynné fáze s vlivem teploty dle \cite noborio.
!     procedure(scalar_fnc), nopass, pointer :: hydraulic_conduc_vapor_invC !< Hydraulická vodivost plynné fáze s vlivem koncentrace dle \cite noborio.
!   end type coupled_fnc


  type, public :: rcza
    logical                                         :: use
    real(kind=rkind)                                :: val
    real(kind=rkind), dimension(:,:), allocatable   :: tab
  end type rcza


!   type, public :: initial
!     real(kind=rkind) :: initcond
!     type(rcza)       :: rcza_set
!     real(kind=rkind) :: top, bottom
!   end type initial



  !> Datový typ globálních charakterystik prostředí
  type, public :: hyd_proper_vals
    real(kind=rkind) :: hc_a !< @f$ a @f$, parametr retenční čáry podle Huston a Cass @f$\label{par:hca}@f$
    real(kind=rkind) :: hc_b !< @f$ b @f$, parametr retenční čáry podle Huston a Cass  @f$\label{par:hcb}@f$
    real(kind=rkind) :: theta_r !< @f$ \theta_r @f$, reziduální objemová vlhkost @f$\label{par:the_r}@f$
    real(kind=rkind) :: theta_s !< @f$ \theta_s @f$, saturovaná objemová vlhkost @f$\label{par:the_s}@f$
!     real(kind=rkind) :: gain_factor !< @f$ G_{K, \Psi} @f$, Zjistit co to přesně  je \cite giakou @f$\label{par:sat} @f$
    !> @f$ G_{K_{\Psi}} @f$, gain faktor [-] @f$\label{par:gain} @f$
    !!
    !! Tento koeficient je vlastně rozšíření tzv. "surface-tension viscous-flow"  (STVF) modelu, který vyjadřuje ve vztahu pro nenasycenou hydraulickou vodivost vodivost (ale například i retenční křivku) vliv teploty na danou křivku pomocí povrchového napětí ("surface tension") a kinematické viskosity ("viscous flow"). Podle STVF modelu ovlivňuje  povrchové napětí retenční křivku a kinematická viskosita nenasycenou hydraulickou vodivost. Takzvaný gain faktor z tohoto modelu vycházi a zpřesňuje ho. Tento popis je převzat z \cite giakou. V tomto modelu je vyjádřen pouze konstantou. Pro hlinité písky je přibližně 1.0 \cite noborio.
    !!
    real(kind=rkind) :: gain_factor 
    real(kind=rkind) :: viscosity_0  !< @f$ \nu_0 @f$, počáteční kinematická viskosita [@f$N \: s \: m^{-2}@f$] @f$\label{par:vis} @f$
    real(kind=rkind) :: molecular_weight !< @f$ M_w @f$, molární hmotnost vody [@f$kg \:mol^{-1}@f$]@f$\label{par:molW}@f$
    real(kind=rkind) :: molecular_weight_NaCl !< @f$M_{salt}@f$, molární hmotnost NaCl [@f$kg \:mol^{-1}@f$] @f$\label{par:weiNaCl}@f$
    real(kind=rkind) :: gas_constant !< @f$ R @f$, plynová konstanta [@f$J \:mol^{-1} \:K^{-1}@f$]@f$\label{par:gasCon}@f$
    !> @f$ G_{\Psi, T} @f$, gain faktor pro vliv teploty [-] @f$\label{par:gainT}@f$
    !!
    !! Popis viz @f$\ref{par:gain}@f$. Tento gain factor se vztahuje k vztahu pro nenasycenou hydraulickou vodivost s vlivem teploty.
    !!
    real(kind=rkind) :: gain_factor_T 
    real(kind=rkind) :: hydrated_iont_r !< @f$ r_s @f$, poloměr hydratováného ionu [@f$m@f$] @f$\label{par:rs}@f$
    real(kind=rkind) :: water_molecul_r !< @f$ r_w @f$, poloměr molekuly vody [@f$m@f$] @f$\label{par:rw}@f$
    !> @f$ b @f$, síla filtru nečistot na částečkách [@f$-@f$] @f$\label{par:b}@f$
    !!
    !! !!! nastřelená hodnota !!!
    !!
    real(kind=rkind) :: half_of_solute_film 
    !> @f$ f_c @f$, podíl jílové frakce v půdě @f$\label{par:clay}@f$
    !!
    !! !!! nastřelená  hodnota !!! 
    !!
    real(kind=rkind) :: clay_fraction 
    real(kind=rkind) :: density_soil_water !< @f$ \rho_l @f$ hustota půdní vody
    real(kind=rkind) :: density_pure_water !< @f$ \rho_w @f$ hustota čisté vody
    real(kind=rkind) :: heat_capacity_soil !< @f$ C_{ps} @f$ tepalná kapacita půd
    real(kind=rkind) :: specific_heat_soil_water
    real(kind=rkind), dimension(:,:), allocatable :: apparent_T_cap_soil
    real(kind=rkind), dimension(:), allocatable :: apparent_T_cap_soil_loc
    real(kind=rkind), dimension(:,:), allocatable :: hdisp !< @f$ D @f$ hydrodynamická disperze
    real(kind=rkind), dimension(:), allocatable :: hdisp_loc
    real(kind=rkind) :: specific_storage
    real(kind=rkind) :: porosity
    real(kind=rkind) :: kappa
    integer :: G_STVF_select
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> diagonal values of the hydraulic conductivity tensor in local system of coordinates
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle of the anisothrophy axes with respect to global axes
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    
    
    real(kind=rkind), dimension(3) :: initcond
    type(rcza), dimension(3)       :: rcza_set
    real(kind=rkind), dimension(3) :: top, bottom
    
    integer :: hdisp_select
    
  end type hyd_proper_vals


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

!   type(coupled_fnc),  public :: main_func !< Deklarace pointru pro konstitunční vztahy.

!   type(satura_vals), dimension(:), pointer, public :: saturation !< Deklarace pointru pro hodnoty saturace a hydraulické vodivosti.
  type(hyd_proper_vals), dimension(:), allocatable, target, public :: hyd_prop !< Deklarace globálních charakterystik prostředí
  character(len=3), public  :: coor
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  integer, public :: file_modwaterm, file_reinitbc, file_heatinitbc, file_soluteinitbc
  
end module modRE_globals