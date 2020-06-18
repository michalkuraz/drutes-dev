module freeze_fnc
  use pde_objs
  use typy
  use freeze_globs
  use debug_tools
  use freeze_helper
  
  public :: capacityhh,  diffhh, diffhT, Dirichlet_mass_bc
  public :: capacityTT, capacityTh, diffTT, convectTT, thermal_k, heat_flux_freeze, Dirichlet_Neumann_switch_bc
  public :: capacityTsTs, thermal_p, diffTsTs, heat_flux_s_LTNE
  public:: T_m, qsl_neg, qsl_pos
        
      
  
  contains
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
    
      if (iceswitch(quadpnt)) then
        val = rho_ice/rho_wat(quadpnt)*vangen_elast_fr(pde_loc, layer, quadpnt)
      else
        val = vangen_elast_fr(pde_loc, layer, quadpnt)
      end if
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
      real(kind=rkind)                :: val, fac
    
      real(kind=rkind) :: temp
    
      temp = pde(heat_proc)%getval(quadpnt)+273.15_rkind
      if (iceswitch(quadpnt)) then
        fac = icefac(quadpnt)
        val = (1-rho_ice/rho_wat(quadpnt))*fac*Lf/temp/grav*&
        vangen_elast_fr(pde_loc, layer, x=(/hl(pde_loc, layer, quadpnt)/))
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
      
      if (present(tensor)) then
        if(present(quadpnt)) then 
          if(fr_rate) then
            call mualem_fr(pde_loc, layer, quadpnt, tensor = tensor)
          else
            call mualem_fr(pde_loc, layer, x = (/hl(pde(wat), layer, quadpnt)/), tensor = tensor)
          end if
          if(isnan(tensor(1,1))) then
           print*, hl(pde(wat), layer, quadpnt)
           print*, "Q reduction", Q_reduction(layer, quadpnt)
           call wait()
          end if
          tensor = 10**(-Omega*Q_reduction(layer, quadpnt))*tensor
        end if
        if (present(x)) then
          print*, "quadpnt required"
          print *, "exited from diffhh::freeze_fnc"
          ERROR STOP
          !call mualem_fr(pde_loc, layer, x = x, tensor = tensor)      
          !tensor = 10**(-Omega*Q_reduction(layer, x = x))*tensor
        end if
      else
        print *, "ERROR! output tensor undefined, exited from freeze_fnc::diffhh"
      end if
     
    end subroutine diffhh
    
    !> diffusion due to temperature for flow model
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
      
      real(kind=rkind), dimension(3,3) :: Klh, Klt, E
      integer(kind=ikind) :: D, i,j
      real(kind=rkind) :: temp, fac
      
      D = drutes_config%dimen

      temp = pde(heat_proc)%getval(quadpnt)+273.15_rkind
      if (present(tensor)) then
        if (present(quadpnt)) then
          call Kliquid_temp(pde_loc, layer, quadpnt, tensor = Klt(1:D, 1:D))
          if(fr_rate) then
            call mualem_fr(pde_loc, layer, quadpnt, tensor = Klh(1:D, 1:D))
          else
            call mualem_fr(pde_loc, layer, x=(/hl(pde(wat), layer, quadpnt)/), tensor = Klh(1:D, 1:D))
          end if
          Klh(1:D,1:D) = 10**(-Omega*Q_reduction(layer, quadpnt))*Klh(1:D, 1:D)
          if(iceswitch(quadpnt)) then
            if(.not. fr_rate) then
              fac = icefac(quadpnt)
              tensor = (Klt(1:D, 1:D) + fac*Lf/temp/grav*Klh(1:D,1:D))
            else
              tensor = Klt(1:D, 1:D)
            end if
          else
            tensor = Klt(1:D, 1:D)
          end if
        end if
      else
         print *, "ERROR! output tensor undefined, exited from freeze_fnc::diffhT"
      end if   

    end subroutine diffhT
    
    
    !> heat: pde(heat_proc)
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
      
      if(.not.iceswitch(quadpnt)) then
        val = 0
      end if
      if(iceswitch(quadpnt)) then
        val = vangen_elast_fr(pde_loc, layer, quadpnt)
      end if
      
      if(.not. fr_rate) then
        val = -val*Lf*rho_ice
      else
        val = 0
      end if
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
      real(kind=rkind)                :: val, fac
      
      real(kind=rkind) :: temp, vol_soil, th_air, thice, C_p, dens_p, thtot, ths, theta
      real(kind=rkind) :: wat_dens
      
      thice = thetai(pde_loc, layer, quadpnt)
      ths = thetas(pde_loc, layer, quadpnt)
      wat_dens = rho_wat(quadpnt)
      vol_soil = 1_rkind - ths
      if(fr_rate) then
        th_air = ths-&
        vangen_fr(pde_loc, layer, quadpnt)
        theta = vangen_fr(pde_loc, layer, quadpnt)
      else
        th_air = ths-thetai(pde_loc, layer, quadpnt)-&
        vangen_fr(pde_loc, layer, x = (/hl(pde(wat), layer, quadpnt)/)) 
        theta = vangen_fr(pde_loc, layer, x = (/hl(pde(wat), layer, quadpnt)/))
      end if

      if(th_air < 0) then
        th_air = 0
!         if(abs(th_air) > epsilon(th_air)) then
!           print*, th_air
!           print*, epsilon(th_air)
!           print *, "the volume of air is negative"
!           print *, "exited from freeze_fnc :: capacityTT"
!           stop
!         end if
      end if
      temp = pde(heat_proc)%getval(quadpnt)+ 273.15_rkind
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
            if(iceswitch(quadpnt)) then
              val = freeze_par(layer)%Ci*rho_ice*thice + val 
            end if
          case ("Snow")
            val = val + freeze_par(layer)%Ca*rho_air*th_air + freeze_par(layer)%Ci*rho_ice*thice
          end select
       end select

       if(.not. fr_rate) then
         if(iceswitch(quadpnt)) then
           fac = icefac(quadpnt)
           val = val + fac*Lf*rho_ice*Lf/temp/grav*&
           vangen_elast_fr(pde_loc, layer, x = (/hl(pde(wat), layer, quadpnt)/))
         end if
       end if


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
      real(kind=rkind), dimension(3) :: thermal_conduct, ths
      integer(kind=ikind) :: D, i
    
      ths = thetas(pde_loc, layer, quadpnt)
     
      D = drutes_config%dimen
      select case (drutes_config%name)
        case ("LTNE")
          if (present(tensor)) then
            do i= 1, D
              tensor(i,i) =  thermal_p(pde_loc,layer, quadpnt)*ths(1)
            end do
          end if
        case("freeze")
          thermal_conduct = thermal_k(pde_loc,layer, quadpnt)
          if (present(tensor)) then
            do i= 1, D
              tensor(i,i) =  thermal_conduct(i)
            end do
          end if
      end select 

    end subroutine diffTT
    
    
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
        call all_fluxes(pde_loc, layer, quadpnt,  flux = flux)
        flux = freeze_par(layer)%Cl *rho_wat(quadpnt)*flux
      end if
        
      if (present(flux_length)) then
        call all_fluxes(pde_loc, layer, quadpnt, flux_length = flux_length)
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
      if(fr_rate) then
        thl = vangen_fr(pde(wat), layer, quadpnt = quadpnt)
      else
        thl = vangen_fr(pde(wat), layer, quadpnt = quadpnt, x=(/hl(pde(wat), layer, quadpnt)/))
      end if
      select case (freeze_par(layer)%material)
        case ("Soil")
          call all_fluxes(pde_loc, layer, quadpnt, flux = flux)
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
      real(kind=rkind) :: h, fac
      real(kind=rkind), dimension(:), allocatable :: gradient, gradientT
      type(integpnt_str) :: quadpnt_loc
      

      
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
        h = hl(pde(wat), layer, quadpnt)
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

      if(iceswitch(quadpnt))then
       if(fr_rate) then
         gradH(1:D) = gradient(1:D)
       else
         fac = icefac(quadpnt)
         gradH(1:D) = gradient(1:D) + fac*Lf/grav*gradientT(1:D)/(pde(heat_proc)%getval(quadpnt) + 273.15_rkind)
       end if
      else
        gradH(1:D) = gradient(1:D)
      end if
      
      if(present(quadpnt)) then
        call pde(wat)%pde_fnc(wat)%dispersion(pde_loc, layer, quadpnt, tensor = Klh(1:D, 1:D))
        call pde(wat)%pde_fnc(heat_proc)%dispersion(pde_loc, layer, quadpnt, tensor = Klt(1:D, 1:D))
      end if
      
      vct(1:D) = matmul(-Klh(1:D,1:D), gradH(1:D))+matmul(-Klt(1:D,1:D), gradientT(1:D))

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
      call diffTT(pde_loc, layer, quadpnt, tensor = thermal_diff)
            
      if (present(flux)) then
        flux = -matmul(thermal_diff(1:D, 1:D), gradT) 
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(matmul(thermal_diff(1:D, 1:D), gradT))
      end if
    end subroutine heat_flux_freeze
    
    
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
      ths = thetas(pde_loc, layer, quadpnt)
      thice = thetai(pde(wat), layer, quadpnt)
      if(fr_rate) then
        thl = vangen_fr(pde(wat), layer, quadpnt)
      else
        thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      end if
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
      ths = thetas(pde_loc, layer, quadpnt)

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
      ths = thetas(pde_loc, layer, quadpnt)
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
      real(kind=rkind) :: val, h, A, Re, Pr, thice, thl, Cp, up, densp, tp, flux_tmp, thair, thtot, ths
      real(kind=rkind), dimension(3) :: flux
      
      integer(kind=ikind) :: D

      D = drutes_config%dimen
      ths = thetas(pde_loc, layer, quadpnt)
      thice = thetai(pde(wat), layer, quadpnt)
      if(fr_rate) then
        thl = vangen_fr(pde(wat), layer, quadpnt)
      else
        thl = vangen_fr(pde(wat), layer, x=(/hl(pde(wat), layer, quadpnt)/))
      end if
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
      tp = thermal_p(pde_loc, layer, quadpnt)
      Pr = Cp*up/tp
      call all_fluxes(pde_loc, layer, quadpnt, flux_length = flux_tmp)
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

    
    ! > Capacity iceswitch
    function cap_ice(pde_loc, layer, quadpnt, x) result(val)
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
      
      val = 1

    end function cap_ice
end module freeze_fnc
