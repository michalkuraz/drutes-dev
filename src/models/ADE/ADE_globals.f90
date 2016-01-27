module ade_globals

  !> ADE solute/material parameters array
  !! D is dispersivity(L)
  !! n is porosity
  !! kd is Freundlich isotherm linear constant
  !! expo is  Freundlich isotherm exponent
  !! rho is density
  !! diff is iont or molecular diffusion
  !! bd is the effective ionic or molecular diffusion coefficient (L^2/T) of the matrix near the interface
  !! coop is \f[ \frac{\beta}{\alpha^2} ratio, based on Gerke approach
  !! T_w - halflife  in water
  !! T_s - halflife on solid phase
  !<
  type, public :: soluteXsoil
    real(kind=rkind) :: Dz
    real(kind=rkind) :: Dx
    real(kind=rkind) :: n
    real(kind=rkind) :: kd
    real(kind=rkind) :: expo
    real(kind=rkind) :: rho
    real(kind=rkind) :: diff
    real(kind=rkind) :: bd
    real(kind=rkind) :: coop
    real(kind=rkind) :: T_w
    real(kind=rkind) :: T_s
    real(kind=rkind) :: icond
  end type soluteXsoil


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf/matrix.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(soluteXsoil), dimension(:), allocatable, public :: ade_par


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf/fractures.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(soluteXsoil), dimension(:), allocatable, public :: ade_par_f

  !> type of used sorption isotherm
  !! 0 - linear
  !! 1 - Friedrich exponential
  !! 2 - Langmuir
  !<
  integer(kind=ikind), public :: isotherm

end module ade_globals