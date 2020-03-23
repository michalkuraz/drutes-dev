module evap_RE_fnc
  public :: evapdiffhh, evapdiffhT
  private :: hydraulic_vh
  
  
    contains
    
      subroutine evapdiffhh(pde_loc, layer, quadpnt,  x, tensor, scalar)
        use typy
        use re_globals
        use pde_objs
        use evapglob
        use re_constitutive

        class(pde_str), intent(in) :: pde_loc
         !> material ID
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), dimension(:), intent(in), optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt      
        !> second order tensor of the unsaturated hydraulic conductivity
        real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
        !> relative hydraulic conductivity, (scalar value)
        real(kind=rkind), intent(out), optional :: scalar
        
        real(kind=rkind), dimension(3,3) :: Ks
        real(kind=rkind) :: Kv
        
        integer(kind=ikind) :: D, i
        
        D = drutes_config%dimen
        
        call mualem(pde(re_ord), layer, quadpnt, tensor=Ks(1:D, 1:D))
        
        Kv = hydraulic_vh(pde(re_ord), layer, quadpnt)
        
        do i = 1,D
          Ks(i,i) = Kv + Ks(i,i)
        end do
        
        tensor(1:D, 1:D) = Ks(1:D, 1:D)
        
      end subroutine evapdiffhh
      
      
     subroutine evapdiffhT(pde_loc, layer, quadpnt,  x, tensor, scalar)
       use typy
       use re_globals
       use pde_objs
       use evapglob
       use re_constitutive
       use re_globals
       use evapextras

       class(pde_str), intent(in) :: pde_loc
        !> material ID
       integer(kind=ikind), intent(in) :: layer
       !> pressure head
       real(kind=rkind), dimension(:), intent(in), optional :: x
       !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
       type(integpnt_str), intent(in), optional :: quadpnt      
       !> second order tensor of the unsaturated hydraulic conductivity
       real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
       !> relative hydraulic conductivity, (scalar value)
       real(kind=rkind), intent(out), optional :: scalar
       
       real(kind=rkind), dimension(3,3) :: Ks
       real(kind=rkind) :: Kv, h, T
      
       integer(kind=ikind) :: D, i
      
       D = drutes_config%dimen
      
       call mualem(pde(re_ord), layer, quadpnt, tensor=Ks(1:D, 1:D))
      
        
       h = pde(re_ord)%getval(quadpnt)
       
       T = pde(heat_ord)%getval(quadpnt)
      
       call mualem(pde(re_ord), layer, quadpnt,  tensor = Ks(1:D,1:D))
       
       tensor(1:D,1:D)  = Ks(1:D,1:D)*h*GwT*(1/gamma_0)*dsurf_tension_soilwat_dT(quadpnt)  
      
      end subroutine evapdiffhT
      
      
      
      function dtheta_vdt(pde_loc, layer, quadpnt, x) result(val)
        use typy
        use global_objs
        use globals
        use pde_objs
        use evapglob
        use evapextras
        
        class(pde_str), intent(in) :: pde_loc
        !> value of the nonlinear function
        real(kind=rkind), dimension(:), intent(in), optional    :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt
        !> material ID
        integer(kind=ikind), intent(in) :: layer
        !> return value
        real(kind=rkind)                :: val
        
        type(integpnt_str) :: quadpnt_loc
        real(kind=rkind) :: thvap_prev, thvap_curr
        
        quadpnt_loc = quadpnt
        
        quadpnt_loc%column = 4
        
        thvap_prev = theta_vapor(pde(re_ord),layer, quadpnt_loc)
        
        quadpnt_loc%column = 3
        
        thvap_curr = theta_vapor(pde(re_ord),layer, quadpnt_loc)
        
        val = -(thvap_curr  - thvap_prev) / time_step

        
      end function dtheta_vdt
        
        
      
      
      
      !> Isothermal Properties of water vapor
      function hydraulic_vh(pde_loc, layer, quadpnt) result(val)
        use typy
        use global_objs
        use pde_objs
        use evapglob
        use re_globals
        use evapextras
        
       
        class(pde_str), intent(in) :: pde_loc
        !> MaterialID
        integer(kind=ikind), intent(in) :: layer
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt 
        !>unsaturated non-thermal conductuvity of water vapor
        real(kind=rkind) :: val
        !> T:temperature
        !> Rh: relatuive humidity of soil
        !> rho_l: liquid water density
        !> rho_sv: saturated vapor density
        !> D: diffusitivity
        real(kind=rkind) :: rh_soil_val, rho_l_val,rho_sv_val,diff, T
        
        
        if (.not. present(quadpnt)) then
          print *, "ERROR: you have not specified  integ point "
          print *, "exited from evap_fnc::hydraulic_vh"
          ERROR stop
        end if
       
        
        rh_soil_val = rh_soil(layer, quadpnt)
        rho_l_val= rho_l( quadpnt) 
        rho_sv_val = rho_sv( quadpnt) 
        diff = vapor_diff_soil(pde(heat_ord), layer, quadpnt)
        
        T = pde(Heat_ord)%getval(quadpnt)
				T = T + Tref
       
        val = (diff/rho_l_val)*rho_sv_val*((MolWat*gravity)/(R_gas*T))*rh_soil_val
    
      end function hydraulic_vh

end module evap_RE_fnc
