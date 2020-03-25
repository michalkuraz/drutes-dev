  module REevapbc
  
    contains
  
  
    !> Evaporation rate [m/s]
    !> Input: Air Relative humiduty [-]
    function evaporation(layer, quadpnt) result(val)
      use typy
      use evapglob
      use evapextras
      use global_objs
        
      !>material ID  
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt 
      !> Evaporation rate [m/s]
      real(kind=rkind) :: val
      !> Relative humidity soil
      !> liquid water density 
      !> saturated water vapor density
      real(kind=rkind) :: rh_soil_val, rho_l_val, rho_sv_val_soil, rho_sv_val_air
      real(kind=rkind), dimension(3,3) :: Ks
      
      type(integpnt_str) :: localpnt
        

      
      rh_soil_val = rh_soil(layer, quadpnt)
      rho_l_val = rho_l(quadpnt) 
      rho_sv_val_soil = rho_sv(quadpnt) 
      
      !air temperature
      localpnt%type_pnt = "numb"
      
      localpnt%this_is_the_value = temp_air
      
      rho_sv_val_air = rho_sv(localpnt) 
            
      val = (rh_soil_val*rho_sv_val_soil  - (rel_air_hum/100.0_rkind)* rho_sv_val_air )/(resistance*rho_l_val)
      
      val = min(val, 0.0_rkind)
    
    end function evaporation
    
    
    subroutine evap4bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs

           
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(integpnt_str) :: quadpnt_loc
      
      integer(kind=ikind) :: layer
      
      if (present(value)) then
        quadpnt_loc%column = 2
        quadpnt_loc%type_pnt = "ndpt"
        quadpnt_loc%order = elements%data(el_id,node_order)
        layer = elements%material(el_id)
        value = evaporation(layer, quadpnt_loc)
      end if
      
      if (present(code)) then
        code = 2
      end if
      
    end subroutine evap4bc
    
  end module REevapbc
