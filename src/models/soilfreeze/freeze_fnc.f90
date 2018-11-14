module freeze_fnc


  contains
    
    function icedensity(quadpnt) result(rho)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in out) :: quadpnt
      real(kind=rkind) :: rho
      
    end function icedensity

end module freeze_fnc
