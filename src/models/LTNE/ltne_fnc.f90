module ltne_fnc


  contains

    pure & !! this is already defined in freeze_fnc.f90
      function  densityFluid(T) result(rho)
      use typy
       
      real(kind=rkind), intent (in) :: T
      real(kind=rkind) :: rho

               rho = (1.682208e-8*T*T*T - 6.05282462e-6*T*T + 2.36680033177935e-5*T + &
               0.999946997406686)*1e3
  end function densityFluid
    
    !> time derivative coefficient for solid temperature in solid equation
    !! \f[ \varepsilon_s \rho_s C_s \f
    !<
    function caps(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use heat_globals
      use pde_objs
      use core_tools
      use re_globals
      use ltne_globals

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

     
      real(kind=rkind) :: watcont
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from re_constitutive::vangen_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_elast"
        ERROR stop
      end if
      
      if (with_richards) then
        watcont = vgset(layer)%Ths
      else
        watcont = ltnepars(layer)%satwatcont
      end if
      
      E = (1-watcont)*ltnepars(layer)%densityS*ltnepars(layer)%heatCapS
      
      
    end function caps
        
    
   function TsinS(pde_loc, layer, quadpnt, x) result(val)
     use typy
     use global_objs
     use pde_objs
     use ADE_globals
     use re_globals
     use debug_tools
     use ltne_globals
     use heat_globals

     class(pde_str), intent(in) :: pde_loc
     !> value of the nonlinear function
     real(kind=rkind), dimension(:), intent(in), optional    :: x
     !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
     type(integpnt_str), intent(in), optional :: quadpnt
     !> material ID
     integer(kind=ikind), intent(in) :: layer
     !> return value
     real(kind=rkind)                :: val 
     
     real(kind=rkind) :: h, Pr, Re, T, watcont, thetas, S, As
     real(kind=rkind), dimension(3) :: q 
     integer(kind=ikind)            :: D

     if (.not. present(quadpnt)) then
       print *, "quadpnt must be present"
       print *, "called from TsinS::ltne_fnc"
       ERROR STOP
     end if
     
     Pr = heatCapL*viscL/CondL
     
     T = pde(Tl_order)%getval(quadpnt)
     
     D = drutes_config%dimen
     
     if (with_richards) then
       call pde(RE_order)%flux(layer,quadpnt, vector_out=q(1:D))
     else
       q(1:D) = ltnepars(layer)%q(1:D)
     end if
     
     if (with_richards) then
        watcont = pde(RE_order)%mass(1)%val(pde(RE_order), layer, quadpnt)
        thetas = vgset(layer)%Ths
     else
        watcont = ltnepars(layer)%actwatcont
        thetas =  ltnepars(layer)%satwatcont
     end if
   
     
     Re = densityFluid(T)*sqrt(dot_product(q(1:D), q(1:D)))/viscL
     
     h = CondL/ltnepars(layer)%grainDiameter*(betah + gammah*Pr**deltah*Re**epsilonh)
     
     S = watcont/thetas
     
     As = 6.0*(1.0-thetas)/ltnepars(layer)%grainDiameter
     
     val = -S*h*As
     
   
   end function TsinS
   
   function TlinS(pde_loc, layer, quadpnt, x) result(val)
     use typy
     use global_objs
     use pde_objs
     use ADE_globals
     use re_globals
     use debug_tools
     use ltne_globals

     class(pde_str), intent(in) :: pde_loc
     !> value of the nonlinear function
     real(kind=rkind), dimension(:), intent(in), optional    :: x
     !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
     type(integpnt_str), intent(in), optional :: quadpnt
     !> material ID
     integer(kind=ikind), intent(in) :: layer
     !> return value
     real(kind=rkind)                :: val 
     
     
     if (.not. present(quadpnt)) then
       print *, "quadpnt must be present"
       print *, "called from TsinS::ltne_fnc"
       ERROR STOP
     end if
     
     val = -TsinS(pde_loc, layer, quadpnt)
     
   end function TlinS
   
   
   function sourceinS(pde_loc, layer, quadpnt, x) result(val)
     use typy
     use global_objs
     use pde_objs
     use ADE_globals
     use re_globals
     use debug_tools
     use ltne_globals

    
     class(pde_str), intent(in) :: pde_loc
     !> value of the nonlinear function
     real(kind=rkind), dimension(:), intent(in), optional    :: x
     !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
     type(integpnt_str), intent(in), optional :: quadpnt
     !> material ID
     integer(kind=ikind), intent(in) :: layer
     !> return value
     real(kind=rkind)                :: val 
     
     val = ltnepars(layer)%Sh
     
   end function sourceinS

end module ltne_fnc
