module ade_globals
  use typy
  
  !> parameters of sorption model
  type, public :: sorption_str
    logical :: kinetic=.false.
    !> name="freu" - Freundlich isoterm, "langmu" - Langmuir isoterm
    character(len=6) :: name
    real(kind=rkind) :: adsorb
    real(kind=rkind) :: desorb
    !>the third parameter in sorption model -- either n exponent in Freundlich or csmax in Langmuir
    real(kind=rkind) :: third
    !> bulk density
    real(kind=rkind) :: bd
    !> ratio  of solid media, if single solid medium ratio=1, if more solid media the sum between sorption(:)%ratio has to be equal 1.0
    real(kind=rkind) :: ratio 
    real(kind=rkind) :: csinit
  end type sorption_str


  !> ADE solute/material parameters array
  !<
  type, public :: soluteXsoil
    real(kind=rkind) :: difmol
    real(kind=rkind), dimension(:), allocatable :: diff_loc
    real(kind=rkind) :: anisoangle
    real(kind=rkind), dimension(:,:), allocatable :: diff
    real(kind=rkind), dimension(:), allocatable :: orders, lambda 
    real(kind=rkind) :: convection
    real(kind=rkind) :: water_cont
    character(len=2) :: icondtype
    real(kind=rkind) :: cmax
    real(kind=rkind) :: cinit
  end type soluteXsoil


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(soluteXsoil), dimension(:), allocatable, public :: adepar
  
  !> structure of sorption parameters and solid medium parameters, the layers are defined in lines, if solid medium scattered into different media, then use solumns
  type(sorption_str), dimension(:,:), allocatable, public :: sorption

  !> type of used sorption isotherm
  !! 0 - linear
  !! 1 - Friedrich exponential
  !! 2 - Langmuir
  !<
  integer(kind=ikind), public :: isotherm
  
  logical, public :: use_richards
  
  integer(kind=ikind) :: no_solids
  
  logical, public :: kinsorb
  
  logical, public :: use_sorption
  
  integer, public :: file_contaminant
  
end module ade_globals
