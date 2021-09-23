module freeze_fnc
  use pde_objs
  use typy
  use freeze_globs
  use debug_tools
  use freeze_helper
  
  public :: capacityhh,  diffhh, diffhT, Dirichlet_mass_bc
  public :: capacityTT, capacityTh, diffTT, convectTT, thermal_k, heat_flux_freeze
  public :: Dirichlet_Neumann_switch_bc, infiltration_bc
  public :: capacityTsTs, thermal_p, diffTsTs, heat_flux_s_LTNE
  public :: T_m, qsl_neg, qsl_pos
  public :: diffTh, dCpdhw, dCpdT
        
      
  
  contains
  
    !> Ice over time for pressure head terms
    function dice_dhw(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
    
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val, C_hl, C_hw

      C_hw = vangen_elast_fr(pde(wat), layer, quadpnt)
      C_hl = vangen_elast_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      val = C_hw-C_hl

    end function dice_dhw
    
    function dice_dT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
    
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      !> water capacity of liquid water
      real(kind=rkind) :: C_hl
      !> temperature in K
      real(kind=rkind) :: temp_K
      
      C_hl = vangen_elast_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      temp_K = pde(heat_proc)%getval(quadpnt)
      val = -C_hl*dhldT(temp_K)
    end function dice_dT
    !> Capacity term due to pressure head for flow model
    !> so pde(wat)
    function capacityhh(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
    
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      !> water caoacity of total water content
      real(kind = rkind) ::C_hw
      C_hw = vangen_elast_fr(pde(wat), layer, quadpnt)
      val = C_hw+dice_dhw(pde_loc, layer, quadpnt)*(rho_ice/rho_wat(quadpnt)-1)
    end function capacityhh
                 
    !> Capacity term due to temperature for flow model
    !> so pde(wat)
    function capacityhT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
    
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
    
      if (iceswitch(quadpnt)) then
        val = (rho_ice/rho_wat(quadpnt)-1)*dice_dT(pde_loc, layer, quadpnt)
      else
        val = 0
      end if
    end function capacityhT

    !> diffusion due to pressure head for flow model
    !> so pde(wat)
    subroutine diffhh(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use freeze_globs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      !> temperature in K
      real(kind=rkind) :: temp_K
      !> viscosity at given temperature
      real(kind=rkind):: u_temp
      
      temp_K = pde(heat_proc)%getval(quadpnt)
      u_temp=exp(ul_a+ul_b/(ul_c+temp_K))/1000_rkind
      
      if (present(tensor)) then
        if(present(quadpnt)) then 
          call mualem_fr(pde(wat), layer, x = (/hl(pde(wat), layer, quadpnt)/), tensor = tensor)
          if(isnan(tensor(1,1))) then
           print*, hl(pde(wat), layer, quadpnt)
           print*, "Q reduction", Q_reduction(layer, quadpnt)
           call wait()
          end if
          tensor = 10**(-Omega*Q_reduction(layer, quadpnt))*tensor*ul_20/u_temp
        end if
        if (present(x)) then
          print*, "quadpnt required"
          print *, "exited from diffhh::freeze_fnc"
          ERROR STOP
        end if
      else
        print *, "ERROR! output tensor undefined, exited from freeze_fnc::diffhh"
      end if
     
    end subroutine diffhh
    
    !> water flow: diffusion due to temperature for flow model
    !> so pde(wat)
    subroutine diffhT(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      !> Hydraulic conductivities
      real(kind=rkind), dimension(3,3) :: Klh, Klt
      integer(kind=ikind) :: D, i,j
      !> temperature in K
      real(kind=rkind) :: temp_K
      !> viscosity at given temperature
      real(kind=rkind):: u_temp
        
      temp_K = pde(heat_proc)%getval(quadpnt)
      
      u_temp=exp(ul_a+ul_b/(ul_c+temp_K))/1000_rkind
  
      D = drutes_config%dimen

      if (present(tensor)) then
        if (present(quadpnt)) then
          call Kliquid_temp(pde(wat), layer, quadpnt, tensor = Klt(1:D, 1:D))
          if(iceswitch(quadpnt)) then
              call mualem_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/), tensor = Klh(1:D, 1:D))
              Klh(1:D,1:D) = 10**(-Omega*Q_reduction(layer, quadpnt))*Klh(1:D, 1:D)*ul_20/u_temp
              tensor = (Klt(1:D, 1:D)*ul_20/u_temp + dhldT(temp_K)*Klh(1:D,1:D))
          else
            tensor = Klt(1:D, 1:D)*ul_20/u_temp
          end if
        end if
      else
         print *, "ERROR! output tensor undefined, exited from freeze_fnc::diffhT"
      end if   
      
    end subroutine diffhT
    
    
    !> heat: pde(heat_proc)
    
        
    function dCpdhw(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      !> capacities
      real(kind = rkind) :: C_hw, C_hl
      

      if(iceswitch(quadpnt)) then
        C_hw = vangen_elast_fr(pde(wat), layer, quadpnt)
        C_hl =  vangen_elast_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
        val = (freeze_par(layer)%Cl*rho_wat(quadpnt)-freeze_par(layer)%Ci*rho_ice)*C_hl
        val = val+(freeze_par(layer)%Ci*rho_ice-freeze_par(layer)%Ca*rho_air)*C_hw
      else
        C_hw = vangen_elast_fr(pde(wat), layer, quadpnt)
        val = (freeze_par(layer)%Cl*rho_wat(quadpnt)-freeze_par(layer)%Ca*rho_air)*C_hw
      end if

    end function dCpdhw
    
        
    function dCpdT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      !> total water capacity
      real(kind = rkind) :: C_hl
      !> temperature in K
      real(kind = rkind) :: temp_K
      
      if(.not.iceswitch(quadpnt)) then
        val = 0
      end if
      if(iceswitch(quadpnt)) then
        temp_K = pde(heat_proc)%getval(quadpnt)
        C_hl =  vangen_elast_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
        val = freeze_par(layer)%Cl*rho_wat(quadpnt)
        val = (val-freeze_par(layer)%Ci*rho_ice)*C_hl*dhldT(temp_K)
        
      end if
    end function dCpdT
    !> Capacity term due to pressure head for heat flow model

    function capacityTh(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      !> total differential of heat capacities over dhw
      real(kind=rkind) :: dcpdhw_val
      !> temperature in K
      real(kind = rkind):: temp_K
      
      val = dice_dhw(pde_loc, layer, quadpnt)*Lf*rho_ice
      temp_K = pde(heat_proc)%getval(quadpnt)
      dcpdhw_val = dCpdhw(pde_loc, layer, quadpnt)
      val = val+temp_K*dcpdhw_val
      
    end function capacityTh
    
    
    !> Capacity term due to temperature for heat flow model

    function capacityTT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      !> temperature in K
      real(kind=rkind) temp_K
      !> volumetric content of soil matrix
      real(kind=rkind) :: vol_soil
      !> volumetric content of air 
      real(kind=rkind):: th_air
      !> volumetric content of ice 
      real(kind=rkind)::thice
      !> heat capacity of the porous medium
      real(kind=rkind):: C_p
      !> density of the porous medium
      real(kind=rkind):: dens_p
      !> total volumetric content (LTNE)
      real(kind=rkind)::thtot
      !> saturated water content
      real(kind=rkind)::ths
      !> volumetric liquid water content
      real(kind=rkind)::theta
      !> water demsity
      real(kind=rkind) :: wat_dens
      
      thice = thetai(pde(wat), layer, quadpnt)
      ths = thetas(pde(wat), layer, quadpnt)
      wat_dens = rho_wat(quadpnt)
      vol_soil = 1_rkind - ths
      theta = vangen_fr(pde(wat), layer, x = (/hl(pde(wat), layer, quadpnt)/))
      th_air = ths-thetai(pde(wat), layer, quadpnt)-theta 
      
      if(th_air < 0) then
        th_air = 0
      end if
      temp_K = pde(heat_proc)%getval(quadpnt)
      select case (drutes_config%name)
        case ("LTNE")
            if(air) then
                thtot = theta+th_air+thice !+freeze_par(layer)%Ci*thice
                C_P = (freeze_par(layer)%Cl*theta&
                +freeze_par(layer)%Ca*th_air+freeze_par(layer)%Ci*thice)/(thtot)
                ! rho_ice*thice+
                dens_p = (wat_dens*theta&
                +rho_air*th_air)/thtot
            else
                thtot = theta + thice
                ! +freeze_par(layer)%Ci*thice
                C_P = (freeze_par(layer)%Cl*theta+freeze_par(layer)%Ci*thice)&
                /(thtot)
                !+rho_ice*thice
                dens_p = (wat_dens* theta+rho_ice*thice)/thtot
            end if
            val = C_P*dens_p
            val = val*thtot          
        case("freeze")
          val =  freeze_par(layer)%Cl*wat_dens*theta
          select case (freeze_par(layer)%material)
          case ("Soil")
            val = val + freeze_par(layer)%Cs*rho_soil*vol_soil + freeze_par(layer)%Ca*rho_air*th_air
            val = freeze_par(layer)%Ci*rho_ice*thice + val 
          case ("Snow")
            val = val + freeze_par(layer)%Ca*rho_air*th_air + freeze_par(layer)%Ci*rho_ice*thice
          end select
       end select

       val = val + temp_K*dCpdT(pde_loc, layer, quadpnt)
       val = val + dice_dT(pde_loc, layer, quadpnt)*Lf*rho_ice


    end function capacityTT
    
    !> dispersion for heat flow model

    subroutine diffTT(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      !> thermal conductivty
      real(kind=rkind), dimension(3) :: thermal_conduct
      integer(kind=ikind) :: D, i
      !> hydraulic conductivity due to thermal gradient
      real(kind=rkind), dimension(3,3)  ::Klt
      !> temperature in K
      real(kind=rkind)::temp_K
      !> saturated water content 
      real(kind = rkind),dimension(3):: ths
      
      temp_K = pde(heat_proc)%getval(quadpnt)

      D = drutes_config%dimen
      select case (drutes_config%name)
        case ("LTNE")
          ths = thetas(pde(wat), layer, quadpnt)
          if (present(tensor)) then
            do i= 1, D
              tensor(i,i) =  thermal_p(pde(heat_proc),layer, quadpnt)*ths(1)
            end do
          end if
        case("freeze")
          thermal_conduct = thermal_k(pde(heat_proc),layer, quadpnt)
          if (present(tensor)) then
            do i= 1, D
              tensor(i,i) =  thermal_conduct(i)
            end do
          end if
      end select 
      call pde(wat)%pde_fnc(heat_proc)%dispersion(pde(wat), layer, quadpnt, tensor = Klt(1:D, 1:D))
      tensor(1:D,1:D) = tensor(1:D,1:D) + Klt(1:D,1:D)*freeze_par(layer)%Cl*rho_wat(quadpnt)*temp_K
    end subroutine diffTT
    
    subroutine diffTh(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      real(kind=rkind), dimension(3) :: thermal_conduct, ths
      integer(kind=ikind) :: D
      !> tensor of liquid hydraulic conductivity
      real(kind=rkind), dimension(3,3)  ::Klh
      !> temperature in K
      real(kind=rkind)::temp_K

      temp_K = pde(heat_proc)%getval(quadpnt)

      D = drutes_config%dimen

      call pde(wat)%pde_fnc(wat)%dispersion(pde(wat), layer, quadpnt, tensor = Klh(1:D, 1:D))
      tensor(1:D,1:D) = Klh(1:D,1:D)*temp_K*freeze_par(layer)%Cl*rho_wat(quadpnt)
      
    end subroutine diffTh
    
    
    subroutine convectTT(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use freeze_globs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
      
      
      if (present(flux)) then
        call all_fluxes(pde(wat), layer, quadpnt,  flux = flux)
        flux = freeze_par(layer)%Cl *rho_wat(quadpnt)*flux
      end if
        
      if (present(flux_length)) then
        call all_fluxes(pde(wat), layer, quadpnt, flux_length = flux_length)
        flux_length = freeze_par(layer)%Cl *rho_wat(quadpnt)*flux_length
      end if
              
    end subroutine convectTT
    
    function thermal_k(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind), dimension(3) :: val 
      real(kind=rkind), dimension(3) :: flux
      real(kind = rkind) :: thl, thice, tk, F
      integer(kind = ikind) :: D, i
      D = drutes_config%dimen
      
      thice = thetai(pde(wat), layer, quadpnt)
      thl = vangen_fr(pde(wat), layer, quadpnt = quadpnt, x=(/hl(pde(wat), layer, quadpnt)/))

      select case (freeze_par(layer)%material)
        case ("Soil")
          call all_fluxes(pde(wat), layer, quadpnt, flux = flux)
          !> hansson changing campbell
          F = 1+ freeze_par(layer)%F1*thice**freeze_par(layer)%F2
          tk = freeze_par(layer)%C1 + freeze_par(layer)%C2*(thl+F*thice)-& 
          (freeze_par(layer)%C1-freeze_par(layer)%C4)*exp(-(freeze_par(layer)%C3*(thl+F*thice))**freeze_par(layer)%C5)
          do i = 1, D
            val(i) = tk + freeze_par(layer)%beta*freeze_par(layer)%Cl*rho_wat(quadpnt)*abs(flux(i))
          end do 
        case("Snow")
          tk = freeze_par(layer)%snow_density**2*2.5e-6-1.23e-4*freeze_par(layer)%snow_density+0.024
          do i = 1, D
            val(i) = tk
          end do 
      end select
      
    end function thermal_k
    
    subroutine all_fluxes(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length

      real(kind=rkind), dimension(3,3)  :: Klh, Klt
      integer                           :: D
      integer(kind=ikind), dimension(3) :: nablaz
      real(kind=rkind), dimension(3)  :: gradH
      real(kind=rkind), dimension(3)  :: vct
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable, save :: gradient, gradientT
      type(integpnt_str) :: quadpnt_loc
      

      if (.not. allocated(gradientT)) allocate(gradientT(drutes_config%dimen))
      if (.not. allocated(gradient)) allocate(gradient(drutes_config%dimen))

      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from freeze_fnc::all_fluxes"
        ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc = quadpnt
        quadpnt_loc%preproc=.true.
        call pde(wat)%getgrad(quadpnt, gradient)
        call pde(heat_proc)%getgrad(quadpnt, gradientT)
      else
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from freeze_fnc::all_fluxes"
          ERROR STOP
        end if
        h = x(1)
        allocate(gradient(ubound(grad,1)))
        gradient = grad
      end if
      
      D = drutes_config%dimen
      
      if(present(quadpnt)) then
        call pde(wat)%pde_fnc(wat)%dispersion(pde(wat), layer, quadpnt, tensor = Klh(1:D, 1:D))
        call pde(wat)%pde_fnc(heat_proc)%dispersion(pde(wat), layer, quadpnt, tensor = Klt(1:D, 1:D))
      end if
      
      vct(1:D) = matmul(-Klh(1:D,1:D), gradient(1:D))+matmul(-Klt(1:D,1:D), gradientT(1:D))

      if (present(flux_length)) then
        flux_length = norm2(vct(1:D))
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if

    end subroutine all_fluxes
    
    subroutine heat_flux_freeze(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use heat_globals
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
    

      real(kind=rkind), dimension(:), allocatable, save :: gradT, gradient
      real(kind=rkind), dimension(3,3) :: thermal_diff, K_Th
      integer(kind = ikind):: D
      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        print *, "exited from heat_fnc::heat_flux"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_flux"
        ERROR stop
      end if   

      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))

      if (present(quadpnt)) then
        call pde(heat_proc)%getgrad(quadpnt, gradT)
        call pde(wat)%getgrad(quadpnt, gradient)
      else
        gradT = grad
      end if
      D = drutes_config%dimen


      call diffTT(pde(heat_proc), layer, quadpnt, tensor = thermal_diff)
      call diffTh(pde(heat_proc), layer, quadpnt, tensor = K_th)

      if (present(flux)) then
        flux = -matmul(thermal_diff(1:D, 1:D), gradT)-matmul(K_th(1:D, 1:D), gradient)
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(flux(1:D))
      end if
    end subroutine heat_flux_freeze
  
    subroutine infiltration_bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
     
      real(kind=rkind), dimension(3,3)  :: Klh
      integer(kind=ikind) :: i, edge_id, j
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, h
      real(kind = rkind), save :: cumfilt
      integer :: i1, D
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      integer(kind = ikind) :: code_tmp = 2
      logical, save:: switch = .false.
      
      
      D = drutes_config%dimen

      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))
        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              bcval = pde_loc%bc(edge_id)%series(j,2)
              value = bcval
              EXIT
            end if
          end do
        else
          quadpnt%type_pnt = "ndpt"
          quadpnt%column = 2 ! otherwise column is random integer number
          quadpnt%order = elements%data(el_id,node_order)
          h = pde(wat)%getval(quadpnt)
          if(h .GE. -1e-1 .and. .not. switch) then
            switch = .true.
          end if
          if(switch) then
            layer = elements%material(el_id)
            call  diffhh(pde_loc, layer, quadpnt, tensor = Klh(1:D,1:D))
            if(Klh(1,1) >   pde_loc%bc(edge_id)%value) then
              bcval = pde_loc%bc(edge_id)%value
              switch = .false.
              code_tmp = 2
            else
              bcval = 0
              code_tmp = 4
            end if
          else
            bcval = pde_loc%bc(edge_id)%value
            code_tmp = 2
          end if
          value = bcval
        end if
      end if
      
      if (present(code)) then
        code = code_tmp
      end if

    end subroutine infiltration_bc
    
    subroutine Dirichlet_mass_bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
     

      integer(kind=ikind) :: i, edge_id, j
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, flux_length, infilt
      real(kind = rkind), save :: cumfilt
      integer :: i1
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      
      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))
        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              bcval = pde_loc%bc(edge_id)%series(j,2)
              value = bcval
              EXIT
            end if
          end do
        else
          quadpnt%type_pnt = "ndpt"
          quadpnt%column = 1 ! otherwise column is random integer number
          quadpnt%order = elements%data(el_id,node_order)
          if(time_step > 0_rkind) then
            layer = elements%material(el_id)
            call  all_fluxes(pde_loc, layer, quadpnt, flux_length = flux_length)
            !flux_length = 1e-3_rkind
          else
            flux_length = 0
            cumfilt = 0
          end if
          infilt = flux_length*time_step
          cumfilt = cumfilt + infilt
          bcval = pde_loc%bc(edge_id)%value-cumfilt
          if(bcval <0) then
            bcval = 0
          end if
          value = bcval
        end if
      end if
      if (present(code)) then
        if(bcval == 0) then
          code = 2
        else
          code = 4
        end if
      end if

    end subroutine Dirichlet_mass_bc
    
    subroutine Dirichlet_Neumann_switch_bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
     

      integer(kind=ikind) :: i, edge_id, j
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, flux_length, infilt
      integer :: i1
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer, code_tmp
      
      if (.not. allocated(pde_common%xvect) ) then
        if (present(value)) value = 0
        if (present(code)) code = 2
        RETURN
      end if
     
      edge_id = nodes%edge(elements%data(el_id, node_order))
      if (pde_loc%bc(edge_id)%file) then
        if (present(value)) then
          edge_id = nodes%edge(elements%data(el_id, node_order))
          i = pde_loc%permut(elements%data(el_id, node_order))
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              bcval = pde_loc%bc(edge_id)%series(j,2)
              value = bcval
              EXIT
            end if
          end do
       end if
        if (present(code)) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              code_tmp = pde_loc%bc(edge_id)%series(j,3)
              code = code_tmp
              EXIT
            end if
          end do
        end if 
      else
      end if
    end subroutine Dirichlet_Neumann_switch_bc
    

    ! Thermal properties LTNE
    function thermal_p(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val 
      real(kind=rkind), dimension(3) :: flux
      real(kind = rkind) :: thl, thice, tk, F, ths
      integer(kind = ikind) :: D, i
      D = drutes_config%dimen
      ths = thetas(pde(wat), layer, quadpnt)
      thice = thetai(pde(wat), layer, quadpnt)
      thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      if(air) then
        !thice*freeze_par(layer)%Li
        val = thl*freeze_par(layer)%Ll+(ths-thl)*freeze_par(layer)%La/ths
      else 
        val = (thl*freeze_par(layer)%Ll)/ vangen_fr(pde(wat), layer, quadpnt)
      end if
    end function thermal_p
    
        ! pde(heat_solid)
    function capacityTsTs(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: temp, vol_soil, th_air, ths
      ths = thetas(pde(wat), layer, quadpnt)

      val = (1-ths)*freeze_par(layer)%Cs

    end function capacityTsTs
    
    subroutine diffTsTs(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use freeze_globs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      real(kind=rkind), dimension(3) :: thermal_conduct, ths
      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen
      ths = thetas(pde(wat), layer, quadpnt)
      if (present(tensor)) then
        do i= 1, D
          select case (freeze_par(layer)%material)
            case ("Soil")
              tensor(i,i) =  freeze_par(layer)%Li * (1-ths(1))
            case("Snow")
              tensor(i,i) =  freeze_par(layer)%Ls * (1-ths(1))
          end select
        end do
      end if
      
      
    end subroutine diffTsTs
    
    subroutine heat_flux_s_LTNE(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use heat_globals
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
    

      real(kind=rkind), dimension(:), allocatable, save :: gradT
      real(kind=rkind), dimension(3,3) :: thermal_diff
      integer(kind = ikind):: D
      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        print *, "exited from heat_fnc::heat_flux"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_flux"
        ERROR stop
      end if   

      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))

      if (present(quadpnt)) then
        call pde_loc%getgrad(quadpnt, gradT)
      else
        gradT = grad
      end if
      
      D = drutes_config%dimen
      
      call diffTsTs(pde(heat_solid), layer, quadpnt, tensor = thermal_diff)
            
      if (present(flux)) then
        flux = -matmul(thermal_diff(1:D, 1:D), gradT) 
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(matmul(thermal_diff(1:D, 1:D), gradT))
      end if
    end subroutine heat_flux_s_LTNE
    
      
    function qsl_pos(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val, h, A, Re, Pr, thice, thl, Cp, up, densp, tp, flux_tmp, thair, thtot, ths
      real(kind=rkind), dimension(3) :: flux
      
      integer(kind=ikind) :: D

      D = drutes_config%dimen
      ths = thetas(pde(wat), layer, quadpnt)
      thice = thetai(pde(wat), layer, quadpnt)
      thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))

      thair = ths-thl
      if(air) then
        thtot = thice+thl+thair
        Cp = (thl*freeze_par(layer)%Cl+thice*freeze_par(layer)%Ci&
        +thair*freeze_par(layer)%Ca)/thtot
        up = (thl*ul+thair*ua)/(thl+thair) !+thice*ui
        densp = (thl*rho_wat(quadpnt)+thice*rho_ice+thair*rho_air)/thtot
      else
        thtot = thice+thl
        Cp = (thl*freeze_par(layer)%Cl+thice*freeze_par(layer)%Ci)/ thtot
        up = ul !+thice*ui
        densp = (thl*rho_wat(quadpnt)+thice*rho_ice)/ thtot
      end if 
      
      A = 6*(1-ths)/freeze_par(layer)%diameter
      tp = thermal_p(pde(heat_proc), layer, quadpnt)
      Pr = Cp*up/tp
      call all_fluxes(pde(wat), layer, quadpnt, flux_length = flux_tmp)
      Re = densp*abs(flux_tmp)*freeze_par(layer)%diameter/up
      h = tp/freeze_par(layer)%diameter*(2.4e-5+(285.6*Pr**2.7*Re**(0.33333333_rkind)))
      val = h * A

    end function qsl_pos
    
    function qsl_neg(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val, h, A, Re, Pr, thice, thl, Cp, up, densp
      real(kind=rkind), dimension(3) :: flux
      integer(kind=ikind) :: D

     val = - qsl_pos(pde_loc, layer, quadpnt)
    end function qsl_neg
    
    !> Mixture temperature
    function T_m(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val 
      real(kind = rkind):: thair, Thsol, thtot, ths
      real(kind=rkind), dimension(3) :: flux
      integer(kind=ikind) :: el_id
      ths = thetas(pde_loc, layer, quadpnt)

      if(air) then
        val = ths*pde(heat_proc)%getval(quadpnt)+(1-ths)*pde(heat_solid)%getval(quadpnt)     
      else
        thtot = vangen_fr(pde(wat), layer, quadpnt)
        thair = ths - thtot
        Thsol = (1-ths)
        if(quadpnt%element >75) then
          el_id = elements%kolik
        else
          el_id = quadpnt%element
        end if
        val = thtot*pde(heat_proc)%getval(quadpnt)+thsol*pde(heat_solid)%getval(quadpnt)+thair*T_air(el_id) 

      end if

    end function T_m

end module freeze_fnc
