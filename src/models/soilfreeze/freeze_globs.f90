module freeze_globs
  use typy
  
  !> gravity acceleration [kg.s^-2]
  real(kind=rkind), parameter :: grav = 9.8196
  
  !>latent heat for water [J.kg^-1]
  real(kind=rkind), parameter :: Lf = 334e3
  

end module freeze_globs
