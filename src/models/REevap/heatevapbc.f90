module heatevapbc

  contains
    !> [W/m^2] Rn
    function netradiation(quadpnt, layer) result(val)
      use typy
      use pde_objs
      use evapglob
      use global_objs
      use re_constitutive
      use re_globals
      
      !>material ID  
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure
      type(integpnt_str), intent(in), optional :: quadpnt 
      real(kind=rkind) :: val
      
      
      !---local I/O variables-----
      !emi_s: emissivity of the soil
      !emi_air: emissivity of the atmosphere
      real(kind = rkind):: emi_s, emi_a,e_a
      real(kind = rkind) :: T,theta 
      real(kind = rkind) :: R_nl, R_ns, temp_air_K

      
      theta = vangen(pde(re_ord), layer, quadpnt)
      T =  pde(Heat_ord)%getval(quadpnt)
      T =  T + Tref
      temp_air_K = temp_air + Tref
     
      
      !net shortwave radiation
      R_ns = (1.0_rkind-wave_albedo)*radiation
      !soil emissivity
      emi_s = min(0.90_rkind + 0.180_rkind*theta, 1.0_rkind)
      !atmospheric vapor pressure
      e_a = (0.611_rkind*exp((17.27_rkind*(temp_air_K -273.15_rkind))/(temp_air_K -35.85_rkind)))*(rel_air_hum/100)
      !atmosphere emissivity 
!      emi_a = 0.70_rkind + 5.95e-5*e_a*exp(1500/temp_air)
      emi_a = 1
      !Net longwave radiation
      R_nl = emi_s*emi_a*5.67e-8*temp_air_K**4 - emi_s*5.67e-8*T**4
      val = R_ns + R_nl
    
    
    end function netradiation
    

      !> Hs - Sensible heat [W/m^2]
    !> Input: Air temperature [K] 
    function sensibleheat(quadpnt) result(val)
      use typy
      use pde_objs
      use evapglob
      use global_objs
      
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in) :: quadpnt 
      real(kind=rkind) :: val
        
      real(kind=rkind) :: T
      
      T = pde(Heat_ord)%getval(quadpnt)
      
      val = C_air*rho_air*((T - temp_air)/resistance)
      
    end function sensibleheat

    !> heat flux for evaporation \f[ G = R_n - H_s - LE_v \f]
    subroutine soil_heat_flux(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      use REevapbc
      use evapextras
      use printtools
            
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      
      
      !------local I/O variables--------
      integer(kind=ikind) :: layer
      real(kind=rkind) :: Rn, Hs, Ev, L
      integer, save :: fileid
      logical, save :: is_opened = .false.
      type(integpnt_str) :: quadpnt_loc
      integer(kind=ikind) :: i
      
      if (.not. is_opened) then
        open(newunit = fileid, file = "out/surface_energy.out", action="write", status = "replace")
        call print_logo(fileid)
        write(unit=fileid, fmt=*) "# time            Rn              -Hs               -LEv                Ev       G"
        is_opened = .true.
      end if

            
      if (present(value)) then
        
        quadpnt_loc%column = 2
        quadpnt_loc%type_pnt = "ndpt"
        quadpnt_loc%order = elements%data(el_id,node_order)
        layer = elements%material(el_id)

        Rn = netradiation(quadpnt_loc, layer) 
        
        Hs = sensibleheat(quadpnt_loc)
        
        Ev = evaporation(layer, quadpnt_loc)
        
  
        L = latent_heat_wat(quadpnt_loc) 
          
        value = Rn - Hs - L*Ev
   
        write(unit=fileid, fmt=*) time, Rn, -Hs, L*Ev, Ev, value
      
      end if
      
      if (present(code)) then
        code = 2
      end if
    
    end subroutine soil_heat_flux

end module heatevapbc
