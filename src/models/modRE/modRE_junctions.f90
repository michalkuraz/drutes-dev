module modRE_junctions 
!  e.g.: disp_R_to_heat means dispersion term in richards equation
!  with respect to heat change. The rest of the functions and subroutines 
!  headings are analoguous


!  water transport
  public :: mass_R
  public :: elas_R_to_potential
  public :: disp_R_to_potential
  public :: disp_R_to_heat
  public :: disp_R_to_solute
  public :: conv_R_to_gravity
!   heat transport
  public :: elas_H_to_heat
  public :: disp_H_to_potential
  public :: disp_H_to_heat
  public :: disp_H_to_solute
  public :: conv_H_to_heat
!   solute transport
  public :: elas_S_to_solute
  public :: disp_S_to_solute
  public :: conv_S_to_solute


  
 contains

!  
!  
!  
!  
!  WATER transport
!  
!  
!  
!  

function mass_R(pde_loc,layer, quadpnt, x) result(res)
    use typy
    use modRE_globals
    use modRE_parameter_functions
    use modRE_constitutive
    use pde_objs
    
    class(pde_str), intent(in) :: pde_loc 
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind) :: res, h
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::mass_R"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::mass_R"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
	print *, "ERROR: van Genuchten function is a function of a single variable h"
	print *, "       your input data has:", ubound(x,1), "variables"
	print *, "exited from modRE_junctions::R_mass"
	ERROR STOP
      end if
      h = x(1)
    end if
    
    res = inv_con_1_HustonAndCass(layer, h)

end function mass_R


function elas_R_to_potential(pde_loc,layer, quadpnt, x) result(res)
    use typy
    use modRE_globals
    use modRE_parameter_functions
    use modRE_constitutive
    use pde_objs
    
    class(pde_str), intent(in) :: pde_loc 
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind) :: h
    !> resulting system elasticity
    real(kind=rkind) :: res
    real(kind=rkind) :: Cw, poro, Ch
    real(kind=rkind) :: theta
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::elas_R_to_potential"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::elas_R_to_potential"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
	print *, "ERROR: Elasticity term in Richards equation"
	print *, "       is a function of a single variable h."
	print *, "       Your input data have:", ubound(x,1), "variable."
	print *, "Exited from modRE_junctions::elas_R_to_potential."
	ERROR STOP
      end if
      h = x(1)
    end if
    
    Cw = hyd_prop(layer)%specific_storage
    poro = hyd_prop(layer)%porosity
    
    theta = inv_con_1_HustonAndCass(layer, h)
    Ch = capillar_capacity(layer, h)
    
    
!     print*, Cw, poro, theta, Ch
    
    
    res = Cw * (theta/poro) + Ch
!     print*, "elas_R_to_potential", res
  end function elas_R_to_potential
 
 
 
  subroutine disp_R_to_potential(pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy
    use globals
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    use modRE_parameter_functions

    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point strucconv_R_to_gravityture (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor    
    real(kind=rkind), intent(out), optional :: scalar
    
    real(kind=rkind) :: h, t, c
    real(kind=rkind) :: rho_l, rho_w
    real(kind=rkind), dimension(3,3) :: K_Psi_l, K_Psi_v
    integer :: d
    d = drutes_config%dimen
    
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_R_to_potential"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_R_to_potential"
      ERROR stop
    end if
      
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Dispersion term with respect to potential"
	print *, "       in Richards equation is a function of three"
	print *, "       variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::disp_R_to_potential."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
    end if
!     print*, "h, t, c  ", h, t, c
    rho_l = hyd_prop(layer)%density_soil_water
    rho_w = rho_pure_water(t)    
    
    K_Psi_l(1:d,1:d) = con_2_HydraulicConductivity(layer, h, t, c)
    K_Psi_v(1:d,1:d) = con_5_HydraulicConductivityVapor(layer, h, t, c)
!     print*, rho_l, rho_w, K_Psi_l(1:d,1:d), K_Psi_v(1:d,1:d)
    if (present(tensor)) then
      tensor = K_Psi_l(1:d,1:d)/rho_l + K_Psi_v(1:d,1:d)/rho_w
!       print*, "disp_R_to_potential", tensor
    end if 
    
    if (present(scalar)) then
!       scalar = norm2(K_Psi_l(1:d,1:d)/rho_l + K_Psi_v(1:d,1:d)/rho_w)
    end if 
  
  
  end subroutine disp_R_to_potential
  
  subroutine disp_R_to_heat(pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy
    use globals
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    use modRE_parameter_functions
    use debug_tools
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor   
    real(kind=rkind), intent(out), optional :: scalar
    
    real(kind=rkind) :: h, t, c
    real(kind=rkind) :: rho_l, rho_w
    real(kind=rkind), dimension(3,3) :: K_T_l, K_T_v
    integer :: d
    d = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_R_to_heat"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_R_to_heat"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Dispersion term with respect to heat"
	print *, "       in Richards equation is a function of three"
	print *, "       variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::disp_R_to_heat."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
    end if
    rho_l = hyd_prop(layer)%density_soil_water
    rho_w = rho_pure_water(t)    
    
    K_T_l(1:d,1:d) = con_3_HydraulicConductivityInvTemperature(layer, h, t, c)
    K_T_v(1:d,1:d) = con_6_HydraulicConductivityVaporInvTemperature(layer, h, t, c)

    
    
    if (present(scalar)) then
!       scalar = norm2(K_T_l(1:d,1:d)/rho_l + K_T_v(1:d,1:d)/rho_w)
    end if 
  
    if (present(tensor)) then
      tensor = K_T_l(1:d,1:d)/rho_l + K_T_v(1:d,1:d)/rho_w
!       print*, "disp_R_to_heat", tensor
    end if 
  
  end subroutine disp_R_to_heat
  
  
  subroutine disp_R_to_solute(pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy
    use globals
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    use modRE_parameter_functions
!     use modRE_parameter_functions    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor   
    real(kind=rkind), intent(out), optional :: scalar
    
    real(kind=rkind) :: h, t, c
    real(kind=rkind) :: rho_l, rho_w
    real(kind=rkind), dimension(3,3) :: K_C_l, K_C_v
    integer :: d
    d = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_R_to_solute"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_R_to_solute"
      ERROR stop
    end if
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Dispersion term with respect to solute"
	print *, "       in Richards equation is a function of three"
	print *, "       variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::disp_R_to_solute."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
    end if
    rho_l = hyd_prop(layer)%density_soil_water
    rho_w = rho_pure_water(t)    

    K_C_l(1:d,1:d) = con_4_HydraulicConductivityInvConcentration(layer, h, t, c)
    K_C_v(1:d,1:d) = con_7_HydraulicConductivityVaporInvConcentration(layer, h, t, c)
    
    if (present(scalar)) then
!       scalar = norm2(K_C_l(1:d,1:d)/rho_l + K_C_v(1:d,1:d)/rho_w)
    end if 
  
    if (present(tensor)) then
      tensor = K_C_l(1:d,1:d)/rho_l +  K_C_v(1:d,1:d)/rho_w
!       print*, "disp_R_to_solute", tensor
    end if 
  
  end subroutine disp_R_to_solute
  
  subroutine conv_R_to_gravity(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
    use typy 
    use globals
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt    
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
    real(kind=rkind), dimension(:), intent(in), optional :: vector_in
    !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
    !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
    !<
    real(kind=rkind), dimension(:), intent(out), optional :: vector_out
    !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
    real(kind=rkind), intent(out), optional :: scalar
    
    real(kind=rkind) :: h, t, c
    real(kind=rkind) :: rho_l
!     real(kind=rkind) :: tmp
    real(kind=rkind), dimension(3,3) :: dK_Psi_l
    real(kind=rkind), dimension(3)  :: wrk
    integer :: d
    d = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from re_constitutive::conv_R_to_gravity"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from re_constitutive::conv_R_to_gravity"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Convection term with respect to gravity"
	print *, "       in Richards equation is a function of three"
	print *, "       variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::conv_R_to_gravity."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
    end if
    
    rho_l = hyd_prop(layer)%density_soil_water
    dK_Psi_l(1:d,1:d) = der_con_2_HydraulicConductivity(layer, h, t, c)/rho_l
    
    wrk = 0.0_rkind
!     wrk(d) = norm2(dK_Psi_l(d,1:d)) !jakub je totak ok?
    
    
    if (present(vector_out)) then
      ! must be negative, because the commnon scheme of the CDE problem has negative convection, but RE has positive convection
      vector_out(1:d) = -wrk(1:d)
!       print*, "conv_R_to_gravity", vector_out
    end if
   
    
    if (present(scalar)) then
!       scalar = -norm2(wrk)
    end if
  end subroutine conv_R_to_gravity
  
!  
!  
!  
!  
!  HEAT transport
!  
!  
!  
!  
  
  function elas_H_to_heat(pde_loc,layer, quadpnt, x) result(res)
    use typy
    use modRE_globals
    use pde_objs
    use modRE_constitutive
    class(pde_str), intent(in) :: pde_loc 
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind) :: res 
    
    res = hyd_prop(layer)%heat_capacity_soil
    
  end function elas_H_to_heat
  
  subroutine disp_H_to_potential(pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy 
    use pde_objs
!     use globals
    use modRE_globals    
    use modRE_constitutive
    use modRE_parameter_functions
    use debug_tools
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
    !> relative hydraulic conductivity, (scalar value)
    real(kind=rkind), intent(out), optional :: scalar   
    
    real(kind=rkind) :: h, t, c
    real(kind=rkind), dimension(3,3) :: wrk
    real(kind=rkind) :: Lw
    integer :: d
    d = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_H_to_potential"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_H_to_potential"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Dispersion term with respect to potential"
	print *, "       in heat transport equation is a function"
	print *, "        of three variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::disp_H_to_potential."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
    end if
    
    Lw =  latent_heat_v(t)

    wrk(1:d,1:d)  = con_5_HydraulicConductivityVapor(layer, h, t, c)

    if (present(scalar)) then
!       scalar = norm2(Lw*wrk(1:d,1:d))
    end if 
  
    if (present(tensor)) then
    !jakub asi je blbost delat z latentniho tepla
    !nebo neni? Je prostredi stejne anisotropni pro hydraulicky vlastnosti 
    !a tepelny?
      tensor = LW*wrk(1:d,1:d)
!       print*, "disp_H_to_potential", tensor
    end if 
!     call  printmtx(A=scalar, name="tensor")
!     print*, LW,tmp,hyd_prop(layer)%Ks
!     call  wait()
    
  end subroutine disp_H_to_potential
  
  subroutine disp_H_to_heat(pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy 
    use pde_objs
    use modRE_globals    
    use modRE_constitutive
    use debug_tools
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
    !> relative hydraulic conductivity, (scalar value)
    real(kind=rkind), intent(out), optional :: scalar
    
    real(kind=rkind), dimension(3,3) :: lambda
    integer :: d
    d = drutes_config%dimen
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_H_to_heat"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_H_to_heat"
      ERROR stop
    end if
    

    lambda(1:d, 1:d) = hyd_prop(layer)%apparent_T_cap_soil
    
    if (present(scalar)) then
!       scalar = norm2(lambda(1:d, 1:d))
    end if 
  
    if (present(tensor)) then
    !jakub asi je blbost delat z tepelny kapacity tenzor
    !nebo neni? Je prostredi stejne anisotropni pro hydraulicky vlastnosti 
    !a tepelny?
      tensor = lambda(1:d, 1:d)
      !print*, "disp_H_to_heat", tensor
    end if 
    
    
  end subroutine disp_H_to_heat
  
  
  
  subroutine disp_H_to_solute(pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy 
    use pde_objs
    use modRE_globals    
    use modRE_constitutive
    use modRE_parameter_functions
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
    !> relative hydraulic conductivity, (scalar value)
    real(kind=rkind), intent(out), optional :: scalar   
    
    real(kind=rkind) :: h, t, c
    real(kind=rkind), dimension(3,3) :: wrk
    real(kind=rkind) :: Lw
    integer :: d
    d = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_H_to_solute"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_H_to_solute"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Dispersion term with respect to solute"
	print *, "       in heat transport equation is a function"
	print *, "       of three variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::disp_H_to_solute."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
    end if
    
    wrk(1:d,1:d)  = con_7_HydraulicConductivityVaporInvConcentration(layer, h, t, c)

    Lw =  latent_heat_v(t)
    
    if (present(scalar)) then
!       scalar = norm2(Lw*wrk(1:d,1:d))
    end if 
  
    if (present(tensor)) then
    !jakub asi je blbost delat z latentniho tepla
    !nebo neni? Je prostredi stejne anisotropni pro hydraulicky vlastnosti 
    !a tepelny?
      tensor = LW*wrk(1:d,1:d)
!       print*, "disp_H_to_solute", tensor
    end if 
  
  end subroutine disp_H_to_solute
  
  subroutine conv_H_to_heat(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
    use typy 
    use globals
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    use modRE_parameter_functions
    use debug_tools
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt    
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
    real(kind=rkind), dimension(:), intent(in), optional :: vector_in
    !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
    !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
    !<
    real(kind=rkind), dimension(:), intent(out), optional :: vector_out
    !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
    real(kind=rkind), intent(out), optional :: scalar
    
    
    real(kind=rkind) :: h, t, c, Cpl
    real(kind=rkind), dimension(:), allocatable :: gradh, gradt, gradc
    real(kind=rkind), dimension(3,3) :: K_Psi_l, K_C_l, K_T_l
    real(kind=rkind), dimension(3) :: vct
    integer :: d
    
    d = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from re_constitutive::conv_H_to_heat"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from re_constitutive::conv_H_to_heat"
      ERROR stop
    end if
    
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
      call pde(1)%getgrad(quadpnt, gradh)
      call pde(2)%getgrad(quadpnt, gradt)
      call pde(3)%getgrad(quadpnt, gradc)
!       print*,  "gradh, gradt, gradc", gradh, gradt, gradc
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Convection term with respect to heat"
	print *, "       in heat transport equation is a function"
	print *, "       of three variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::conv_H_to_heat."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
      !jakub jak tu zadat gradienty ?
    end if
    
   ! taky to spočítá tensor pro danné K
    K_Psi_l(1:d,1:d) = con_2_HydraulicConductivity(layer, h, t, c)
    K_T_l(1:d,1:d) = con_3_HydraulicConductivityInvTemperature(layer, h, t, c)
    K_C_l(1:d,1:d) = con_4_HydraulicConductivityInvConcentration(layer, h, t, c)
    Cpl = spec_heat_solut(c,layer)
!     print*, h, t, c
!     print*, K_Psi_l(1:d,1:d), K_T_l(1:d,1:d), K_C_l(1:d,1:d) , Cpl
    vct = 0.0_rkind
!     vct(d) = norm2(K_Psi_l(d,1:d))  !jakub vliv gravitace je ok?
    vct(1:d) = vct(1:d)  + matmul(K_Psi_l(1:d,1:d),gradh(1:d))
    vct(1:d) = vct(1:d)  + matmul(K_T_l(1:d,1:d), gradt(1:d))
    vct(1:d) = vct(1:d)  + matmul(K_C_l(1:d,1:d), gradc(1:d))
    vct(1:d) = vct(1:d) * Cpl
    
    !jakub je to ok?
    if (present(scalar)) then
!       scalar = -norm2(vct(1:d))
    end if
    
    if (present(vector_out)) then
      ! must be negative, because the commnon scheme of the CDE problem has negative convection, but RE has positive convection
      vector_out = -vct(1:d)
      !print*, "conv_H_to_heat", vector_out
    end if

    
  end subroutine conv_H_to_heat
  

!  
!  
!  
!  SOLUTE transport
!  
!  
!  
! 

  
  function elas_S_to_solute(pde_loc,layer, quadpnt, x) result(res)
    use typy
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    
    class(pde_str), intent(in) :: pde_loc 
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind) :: res 
    real(kind=rkind) :: h
    real(kind=rkind) :: theta
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::elas_S_to_solute"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::elas_S_to_solute"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
    else
      h = x(1)
      if (ubound(x,1) /=1) then
	print *, "ERROR: Elasticity term in solute transport equation"
	print *, "       is a function of a single variable h."
	print *, "       Your input data have:", ubound(x,1), "variable."
	print *, "Exited from modRE_junctions::elas_S_to_solute."
	ERROR STOP
      end if
    end if
    
    theta = inv_con_1_HustonAndCass(layer, h)
    
    res = hyd_prop(layer)%density_soil_water*theta 
    
  end function elas_S_to_solute
  
  
  subroutine disp_S_to_solute (pde_loc, layer, quadpnt,  x, tensor, scalar)
    use typy 
    use pde_objs
    use modRE_globals   
    use modRE_constitutive
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt      
    !> second order tensor of the unsaturated hydraulic conductivity
    real(kind=rkind), dimension(:,:), intent(out), optional :: tensor    
    real(kind=rkind), intent(out), optional :: scalar      
    real(kind=rkind) :: h, t, c
    real(kind=rkind), dimension(:), allocatable :: gradh, gradt, gradc
    real(kind=rkind) :: theta,  rho_l
    real(kind=rkind), dimension(3,3) :: D
    integer :: i, dd
    dd = drutes_config%dimen
    
    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from modre_constitutive::disp_S_to_solute"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from modre_constitutive::disp_S_to_solute"
      ERROR stop
    end if
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
      call pde(1)%getgrad(quadpnt, gradh)
      call pde(2)%getgrad(quadpnt, gradt)
      call pde(3)%getgrad(quadpnt, gradc)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Dispersion term with respect to solute"
	print *, "       in solute transport equation is a function"
	print *, "       of single variables h."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::disp_S_to_solute."
	ERROR STOP
      end if
      h = x(1)
      t = x(1)
      c = x(1)
    end if
    
    rho_l = hyd_prop(layer)%density_soil_water
    theta = inv_con_1_HustonAndCass(layer, h)
    
    if (hyd_prop(layer)%hdisp_select == 0) then
      D(1:dd,1:dd) = hyd_prop(layer)%hdisp
    end if
    
    if (hyd_prop(layer)%hdisp_select == 1) then
      D(1:dd,1:dd) = diff_coef(layer, h, t, c, gradh, gradt, gradc)
    end if
    
    
    if (present(tensor)) then
      tensor = theta*rho_l*D(1:dd,1:dd)
            !print*, "disp_S_to_solute", tensor
    end if 
    
    if (present(scalar)) then
!       scalar = norm2(theta*rho_l*D(1:dd,1:dd))
    end if 
  
  end subroutine disp_S_to_solute
!   
  
 subroutine conv_S_to_solute(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
    use typy 
    use globals
    use pde_objs
    use modRE_globals
    use modRE_constitutive
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    type(integpnt_str), intent(in), optional :: quadpnt    
    !> pressure head
    real(kind=rkind), intent(in), dimension(:), optional :: x
    !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
    real(kind=rkind), dimension(:), intent(in), optional :: vector_in
    !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
    !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
    !<
    real(kind=rkind), dimension(:), intent(out), optional :: vector_out
    !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
    real(kind=rkind), intent(out), optional :: scalar
    
    integer :: d
    real(kind=rkind) :: h, t, c
    real(kind=rkind), dimension(:), allocatable :: gradh, gradt, gradc
    real(kind=rkind), dimension(3,3) :: K_Psi_l, K_C_l, K_T_l
    real(kind=rkind), dimension(3) :: vct


    if (present(quadpnt) .and. present(x)) then
      print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
      print *, "exited from re_constitutive::conv_S_to_solute"
      ERROR stop
    else if (.not. present(quadpnt) .and. .not. present(x)) then
      print *, "ERROR: you have not specified either integ point or x value"
      print *, "exited from re_constitutive::conv_S_to_solute"
      ERROR stop
    end if
    
    d = drutes_config%dimen
    
    if (present(quadpnt)) then
      h = pde(1)%getval(quadpnt)
      t = pde(2)%getval(quadpnt)
      c = pde(3)%getval(quadpnt)
      call pde(1)%getgrad(quadpnt, gradh)
      call pde(2)%getgrad(quadpnt, gradt)
      call pde(3)%getgrad(quadpnt, gradc)
    else
      if (ubound(x,1) /=3) then 
	print *, "ERROR: Convection term with respect to solute"
	print *, "       in solute transport equation is a function"
	print *, "       of three variables h, t, c."
	print *, "       Your input data have:", ubound(x,1), "variables"
	print *, "Exited from modRE_junctions::conv_S_to_solute."
	ERROR STOP
      end if
      h = x(1)
      t = x(2)
      c = x(3)
      !jakub jak tu zadat gradienty ?
    end if
    
   ! taky to spočítá tensor pro danné K
    K_Psi_l(1:d,1:d) = con_2_HydraulicConductivity(layer, h, t, c)
    K_T_l(1:d,1:d) = con_3_HydraulicConductivityInvTemperature(layer, h, t, c)
    K_C_l(1:d,1:d) = con_4_HydraulicConductivityInvConcentration(layer, h, t, c)
    
    vct = 0.0_rkind
!     vct(d) = norm2(K_Psi_l(d,1:d))  !jakub vliv gravitace je ok?
    vct(1:d) = vct(1:d)  + matmul(K_Psi_l(1:d,1:d),gradh(1:d))
    vct(1:d) = vct(1:d)  + matmul(K_T_l(1:d,1:d), gradt(1:d))
    vct(1:d) = vct(1:d)  + matmul(K_C_l(1:d,1:d), gradc(1:d))
    
    
    if (present(vector_out)) then
      ! must be negative, because the commnon scheme of the CDE problem has negative convection, but RE has positive convection
      vector_out = -vct(1:d)
    end if

    if (present(scalar)) then
!       scalar = -norm2(vct(1:d))
    end if
    
    
  end subroutine conv_S_to_solute
  

!   
! 
!   subroutine water_flux_liquid(layer, h, gradient,  flux, flux_length)
!       use typy
!       use globals
!       use modRE_globals
!       use global_objs
!       use modRE_constitutive
! 
!       integer(kind=ikind), intent(in)                          :: layer
!       real(kind=rkind), intent(in)                             :: h
!       !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
!       real(kind=rkind), dimension(:), intent(in), optional     :: gradient
!       real(kind=rkind), dimension(:), intent(out), optional    :: flux
!       real(kind=rkind), intent(out), optional                  :: flux_length
! 
!       real(kind=rkind), dimension(3,3)  :: K_Psi_l
!       integer                           :: D
!       integer(kind=ikind)               :: i
!       integer(kind=ikind), dimension(3) :: nablaz
!       real(kind=rkind), dimension(3)  :: gradtheta
!       real(kind=rkind), dimension(3)  :: vct
! 
!       if (.not. present(gradient)) then
!         print *, "I: interrupted from modre_constitutive::darcy_law"
!         print *, "ERROR! this function was called without gradient of h vector"
!         ERROR STOP
!       end if
! 
!       D = drutes_config%dimen
! 
!       nablaz = 0
!       nablaz(D) = 1
!       
!       gradtheta(1:D) = gradient(1:D) + nablaz(1:D)
!       
!       K_Psi_l(1:D,1:D) = hyd_prop(layer)%Ks(1:D,1:D)*con_2_HydraulicConductivity(layer, h)
! 
! !       call pde(1)%pde_fnc(1)%dispersion(layer, h, K(1:D, 1:D))
! 
!       vct(1:D) = matmul(-K_Psi_l(1:D,1:D), gradtheta(1:D))
! ! 
! ! 
!       if (present(flux_length)) then
!         select case(D)
!           case(1)
!                 flux_length = vct(1)
!           case(2)
!                 flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
!           case(3)
!                 flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
!         end select
!       end if
! ! 
! ! 
!       if (present(flux)) then
!         flux(1:D) = vct(1:D)
!       end if
!     
!   
!   end subroutine water_flux_liquid
  
  
end module modRE_junctions