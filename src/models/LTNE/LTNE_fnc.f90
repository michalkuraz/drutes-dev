module LTNE_fnc
  use pde_objs
  use typy
  use LTNE_globs
  use debug_tools
  use LTNE_helper
  
  public ::  all_fluxes_LTNE, Dirichlet_mass_bc, Dirichlet_Neumann_switch_bc
  public :: capacityhh, diffhh, capacityTlTl, capacityTlh, diffTlTl, convectTlTl, thermal_p
  public:: heat_flux_l_LTNE, heat_flux_s_LTNE, qsl_pos, qsl_neg, T_m

  
  procedure(scalar_fnc), pointer, public :: rwcap
      
      
  
  contains
    !> Capacity term due to pressure head for flow model
    !> so pde(1)
    function capacityhh(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use LTNE_globs
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
        val = rho_ice/rho_wat*rwcap(pde_loc, layer, quadpnt)
      else
        val = rwcap(pde_loc, layer, quadpnt)

      end if

    end function capacityhh


    !> diffusion due to pressure head for flow model
    !> so pde(1)
    subroutine diffhh(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use LTNE_globs
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
      
      if (present(tensor)) then
        if(present(quadpnt)) then 
          call mualem_ltne(pde_loc, layer, x = (/hl(pde(1), layer, quadpnt)/), tensor = tensor)
          tensor = 10**(-Omega*Q_reduction(layer, quadpnt))*tensor
        end if
        if (present(x)) then
          call mualem_ltne(pde_loc, layer, x = x, tensor = tensor)
          tensor = 10**(-Omega*Q_reduction(layer, x = x))*tensor
        end if
      else
        print *, "ERROR! output tensor undefined, exited from LTNE_fnc:diffhh"
      end if

    end subroutine diffhh
    
    subroutine diffhT(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use ltne_globs
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
      
      real(kind=rkind), dimension(3,3) :: Klh, Klt, E
      integer(kind=ikind) :: D, i,j
      real(kind=rkind) :: temp
      
      D = drutes_config%dimen

      temp = pde(2)%getval(quadpnt)+273.15_rkind
      if (present(tensor)) then
        if (present(quadpnt)) then
          call Kliquid_temp(pde_loc, layer, quadpnt, tensor = Klt(1:D, 1:D))
          call mualem_ltne(pde_loc, layer, x=(/hl(pde(1), layer, quadpnt)/), tensor = Klh(1:D, 1:D))
          Klh(1:D,1:D) = 10**(-Omega*Q_reduction(layer, quadpnt))*Klh(1:D, 1:D)
          if(iceswitch(quadpnt)) then
            tensor = (Klt(1:D, 1:D) + Lf/temp/grav*Klh(1:D,1:D))
          else
            tensor = Klt(1:D, 1:D)
          end if
        end if
      else
         print *, "ERROR! output tensor undefined, exited from Ltne_fnc::diffhT"
      end if   
      
      end subroutine diffhT
    
    !> heat: pde(2)
    !> Capacity term due to pressure head for heat flow model

    function capacityTlh(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use LTNE_globs
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
      
      if(.not.iceswitch(quadpnt)) then
        val = 0
      end if
      if(iceswitch(quadpnt)) then
        val = rwcap(pde_loc, layer, quadpnt)
      end if
      val = -val*Lf*rho_ice
      
      !val = 0
    end function capacityTlh
    
    !> Capacity term due to temperature for heat flow model

    function capacityTlTl(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use LTNE_globs
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
      
      real(kind=rkind) :: temp, vol_soil, th_air
      
      vol_soil = 1_rkind - LTNE_par(layer)%Ths
      th_air = LTNE_par(layer)%Ths-thetai(pde_loc, layer, quadpnt)-&
      vangen_ltne(pde_loc, layer, x = (/hl(pde(1), layer, quadpnt)/)) 
      if(th_air < 0) then
        if(abs(th_air) > epsilon(th_air)) then
          print*, th_air
          print*, epsilon(th_air)
          print *, "the volume of air is negative"
          print *, "exited from LTNE_fnc :: capacityTT"
          stop
        end if
      end if
      temp = pde(2)%getval(quadpnt)+ 273.15_rkind
      val =  ltne_par(layer)%Cl*rho_wat*vangen_ltne(pde_loc, layer, x = (/hl(pde(1), layer, quadpnt)/)) 
      if(air) then
        val = val + ltne_par(layer)%Ci*rho_ice*thetai(pde_loc, layer, quadpnt) + ltne_par(layer)%Ca*rho_air*th_air
      else
        val = val + ltne_par(layer)%Ci*rho_ice*thetai(pde_loc, layer, quadpnt) 
      end if
      val = val*LTNE_par(layer)%Ths
      if(iceswitch(quadpnt)) then
        val = val + Lf*rho_ice*Lf/temp/grav*rwcap(pde_loc, layer, x = (/hl(pde(1), layer, quadpnt)/))
      end if
      
    end function capacityTlTl
    
    !> dispersion for heat flow model

    subroutine diffTlTl(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use LTNE_globs
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
      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen
      
      if (present(tensor)) then
        do i= 1, D
          tensor(i,i) =  thermal_p(pde_loc,layer, quadpnt)*ltne_par(layer)%Ths
        end do
      end if
      
      
    end subroutine diffTlTl
    
    subroutine convectTlTl(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use LTNE_globs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
      
      
      if (present(flux)) then
          call all_fluxes_LTNE(pde_loc, layer, quadpnt,  flux = flux)
          flux = ltne_par(layer)%Cl*rho_wat*flux
          
        end if
        
        if (present(flux_length)) then
           call all_fluxes_LTNE(pde_loc, layer, quadpnt, flux_length = flux_length)
           flux_length = ltne_par(layer)%Cl *rho_wat*flux_length
        end if
              
    end subroutine convectTlTl
    
    ! PDE(3)
    function capacityTsTs(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use LTNE_globs
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
      
      real(kind=rkind) :: temp, vol_soil, th_air
      
      val = (1-ltne_par(layer)%Ths)*ltne_par(layer)%Cs
      
    end function capacityTsTs
    
    subroutine diffTsTs(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use LTNE_globs
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
      real(kind=rkind), dimension(3) :: thermal_conduct
      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen
      
      if (present(tensor)) then
        do i= 1, D
          tensor(i,i) =  ltne_par(layer)%Li * (1-ltne_par(layer)%Ths)
        end do
      end if
      
      
    end subroutine diffTsTs
    
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
      real(kind = rkind) :: thl, thice, tk, F
      integer(kind = ikind) :: D, i
      D = drutes_config%dimen
      
      thice = thetai(pde(1), layer, quadpnt)
      thl = vangen_ltne(pde(1), layer, x=(/hl(pde(1), layer, quadpnt)/))
      if(air) then
        val = thl*ltne_par(layer)%Ll+thice*ltne_par(layer)%Li+(ltne_par(layer)%ths-thl)*LTNE_par(layer)%La
      else 
        val = (thl*ltne_par(layer)%Ll+thice*ltne_par(layer)%Li)/ vangen_ltne(pde(1), layer, quadpnt)
      end if
      
    end function thermal_p
    
    subroutine all_fluxes_LTNE(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
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
      real(kind=rkind), dimension(:), allocatable :: gradient, gradientT
      type(integpnt_str) :: quadpnt_loc
            
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from LTNE_fnc::all_fluxes"
        ERROR stop
      end if

      if (present(quadpnt)) then
        quadpnt_loc = quadpnt
        quadpnt_loc%preproc=.true.
        h = hl(pde(1), layer, quadpnt)
        call pde(1)%getgrad(quadpnt, gradient)
        call pde(2)%getgrad(quadpnt, gradientT)
      else
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from LTNE_fnc::all_fluxes"
          ERROR STOP
        end if
        h = x(1)
        allocate(gradient(ubound(grad,1)))
        gradient = grad
      end if

      D = drutes_config%dimen
      if(iceswitch(quadpnt))then
        gradH(1:D) = gradient(1:D) + Lf/grav*gradientT(1:D)/(pde(2)%getval(quadpnt) + 273.15_rkind)
      else
        gradH(1:D) = gradient(1:D)
      end if
      
      if(present(quadpnt)) then
        call pde(1)%pde_fnc(1)%dispersion(pde_loc, layer, x=(/hl(pde(1), layer, quadpnt)/), tensor=Klh(1:D, 1:D))
        Klh(1:D, 1:D) = 10**(-Omega*Q_reduction(layer, quadpnt))*Klh(1:D, 1:D)
      end if
      
      vct(1:D) = matmul(-Klh(1:D,1:D), gradH(1:D))
      if (present(flux_length)) then
        select case(D)
          case(1)
                flux_length = vct(1)
          case(2)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
          case(3)
                flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
        end select

      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if

    end subroutine all_fluxes_LTNE
        
    subroutine heat_flux_l_LTNE(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
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
      !!!!! change here
      call diffTlTl(pde_loc, layer, quadpnt, tensor = thermal_diff)
            
      if (present(flux)) then
        flux = -matmul(thermal_diff(1:D, 1:D), gradT) 
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(matmul(thermal_diff(1:D, 1:D), gradT))
      end if
    end subroutine heat_flux_l_LTNE
    
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
      
      call diffTsTs(pde_loc, layer, quadpnt, tensor = thermal_diff)
            
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
      real(kind=rkind) :: val, h, A, Re, Pr, thice, thl, Cp, up, densp, tp, flux_tmp
      real(kind=rkind), dimension(3) :: flux
      
      integer(kind=ikind) :: D

      D = drutes_config%dimen

      thice = thetai(pde(1), layer, quadpnt)
      thl = vangen_ltne(pde(1), layer, x=(/hl(pde(1), layer, quadpnt)/))
      
      if(air) then
        Cp = thl*ltne_par(layer)%Cl+thice*ltne_par(layer)%Ci +(ltne_par(layer)%ths-thl)*LTNE_par(layer)%Ca
        up = thl*ul+(ltne_par(layer)%ths-thl)*ua !+thice*ui
        densp = thl*rho_wat+thice*rho_ice+(ltne_par(layer)%ths-thl)*rho_air
      else
       Cp = (thl*ltne_par(layer)%Cl+thice*ltne_par(layer)%Ci)/ vangen_ltne(pde(1), layer,  quadpnt)
       up = thl*ul !+thice*ui
       densp = (thl*rho_wat+thice*rho_ice)/ vangen_ltne(pde(1), layer,  quadpnt)
      end if 
      
      A = 6*(1-ltne_par(layer)%Ths)/LTNE_par(layer)%diameter
      tp = thermal_p(pde_loc, layer, quadpnt)
      Pr = Cp*up/tp
      call all_fluxes_LTNE(pde_loc, layer, quadpnt, flux_length = flux_tmp)
      Re = densp*abs(flux_tmp)*LTNE_par(layer)%diameter/up
      h = tp/LTNE_par(layer)%diameter*(2.4e-5+(285.6*Pr**2.7*Re**(0.33333333_rkind)))
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
    
    
        
    
    function T_m(pde_loc, layer, quadpnt, x) result(val)
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

      if(air) then
        val = ltne_par(layer)%ths*pde(2)%getval(quadpnt)+(1-Ltne_par(layer)%Ths)*pde(3)%getval(quadpnt)     
      else
        val = vangen_ltne(pde(1), layer, quadpnt)*pde(2)%getval(quadpnt)+(1-Ltne_par(layer)%Ths)*pde(3)%getval(quadpnt)     
      end if

    end function T_m
    
    subroutine Dirichlet_mass_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
     

      integer(kind=ikind) :: i, edge_id, j
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, flux_length, infilt
      integer :: i1
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      
      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))
        print*, edge_id
        call wait()
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
            call  all_fluxes_LTNE(pde_loc, layer, quadpnt, flux_length = flux_length)
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
    
    subroutine Dirichlet_Neumann_switch_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
     

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

end module LTNE_fnc
