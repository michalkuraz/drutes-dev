module evap_RE_fnc
  public :: evapdiffhh, evapdiffhT, dtheta_vdt
  
  
    contains
    
      subroutine evapdiffhh(pde_loc, layer, quadpnt,  x, tensor, scalar)
        use typy
        use re_globals
        use pde_objs
        use evapglob
        use re_constitutive
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
        real(kind=rkind) :: Kv
        
        integer(kind=ikind) :: D, i
        
        D = drutes_config%dimen
        
        call mualem(pde(re_ord), layer, quadpnt, tensor=Ks(1:D, 1:D))
        
        Kv = hydraulic_vh(layer, quadpnt)
        
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
       
       real(kind=rkind), dimension(3,3) :: Klt
       real(kind=rkind) :: KvT_scalar, h, T
      
       integer(kind=ikind) :: D, i
      
       D = drutes_config%dimen
      
       KlT(1:D,1:D) = hydraulic_lT(layer, quadpnt) 
       KvT_scalar = hydraulic_vT(layer, quadpnt)
      
       tensor(1:D, 1:D) = Klt(1:D, 1:D)
       
       do i=1, D
         tensor(i,i) = tensor(i,i) + KvT_scalar 
       end do
      
      
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
      
      quadpnt_loc%column = 1
      
      thvap_curr = theta_vapor(pde(re_ord),layer, quadpnt_loc)
      
      val = -(thvap_curr  - thvap_prev) / time_step

      
    end function dtheta_vdt
      
      
    !> total water flux
    subroutine water_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use re_constitutive
      use evapglob
      use evapextras
      
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> local variable
      integer(kind=ikind)  :: i, D
      !> Liquid water flux
      real(kind=rkind), dimension(3)  ::  q_liq, q_l, q_v
      !> resul of the modified flux of liquid water
      real(kind=rkind), dimension(3)  :: vct
      
      if (present(x)) then
        print *, "ERROR: use quadpnt only"
        print *, "exited from evap_RE_fnc::water_flux"
        ERROR stop
      end if
    
            
      D = drutes_config%dimen

      call liquid_flux(pde(re_ord), layer, quadpnt, flux=q_l(1:D))
      
      call vapor_flux(pde(re_ord), layer, quadpnt, flux=q_v(1:D))
      vct(1:D) = q_l(1:D) +q_v(1:D)
      
      if (present(flux_length)) then      
        flux_length = norm2(vct(1:D))
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    
    end subroutine water_flux

      !> Liquid water flux
    subroutine liquid_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use re_constitutive
      use evapglob
      use evapextras
      
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> Klt: total unsaturated non-thermal conductivity of liquid water
      real(kind=rkind), dimension(3,3)  :: KlT
      !> local variable
      integer(kind=ikind)  :: i, D
      !> Liquid water flux
      real(kind=rkind), dimension(3)  ::  q_liq
      !> resul of the modified flux of liquid water
      real(kind=rkind), dimension(3)  :: vct
      !> pressure head
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable :: gradient
      !> Gauss quadrature point structure local
      type(integpnt_str) :: quadpnt_loc
      !> Temperature gradient
      real(kind=rkind), dimension(:), allocatable, save :: gradT
      
      if (present(x)) then
        print *, "ERROR: use quadpnt only"
        print *, "exited from evap_fnc::liquid_flux"
        ERROR stop
      end if
    
      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))
      
      if (present(quadpnt)) then
        call pde(heat_ord)%getgrad(quadpnt, gradT)
      else
        gradT = grad
      end if
      
      D = drutes_config%dimen
      
      if (present(quadpnt)) then
        call darcy_law(pde(re_ord), layer, quadpnt, flux = q_liq(1:D))
      end if
      
      KlT(1:D,1:D) = hydraulic_lT(layer, quadpnt) 
      
      vct(1:D) = q_liq(1:D) + matmul(-KlT(1:D,1:D), gradT(1:D))
      
      if (present(flux_length)) then      
        flux_length = norm2(vct(1:D))
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    
    end subroutine liquid_flux
        
        
       !> Water Vapour flux
    subroutine vapor_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use evapglob
      use evapextras
       
      class(pde_str), intent(in) :: pde_loc
      !> Material ID
      integer(kind=ikind), intent(in) :: layer
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt   
      !> value of the nonlinear function
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional :: grad
      !> Vector of the flux
      real(kind=rkind), dimension(:), intent(out), optional :: flux
      !> lengh of the flux vector
      real(kind=rkind), intent(out), optional :: flux_length
      !> KvT: unsaturated thermal conductivity for water
      !> Kvh: unsaturated non-thermal  conductivity for water
      real(kind=rkind), dimension(3,3)  :: Kvh, KvT
      !> local variables
      integer(kind=ikind) :: i, D
      !> pressure gradient
      real(kind=rkind), dimension(:), allocatable :: gradient
      !result of the vapor flux vector
      real(kind=rkind), dimension(3) :: vct
      !> h: pressure head
      !> Kvh_scalar: unsaturated non-thermal  conductivity for water
      !> KvT_scalar: unsaturated thermal conductivity for water
      real(kind=rkind) :: h, Kvh_scalar, KvT_scalar
      !> Gauss quadrature point structure local
      type(integpnt_str) :: quadpnt_loc
      
      
     real(kind=rkind), dimension(:), allocatable, save :: gradT
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from evap_fnc::liquid_flux"
        ERROR stop
      end if
    
      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))
      
      if (present(quadpnt)) then
        call pde(heat_ord)%getgrad(quadpnt, gradT)
      else
        gradT = grad
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
        quadpnt_loc%preproc=.true.
        h = pde(re_ord)%getval(quadpnt_loc)
        call pde(re_ord)%getgrad(quadpnt, gradient)
      else
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from re_constitutive::darcy_law"
          ERROR STOP
        end if
        h = x(1)
        allocate(gradient(ubound(grad,1)))
        gradient = grad
      end if
      D = drutes_config%dimen
      
      Kvh_scalar = hydraulic_vh(layer, quadpnt)
      KvT_scalar = hydraulic_vT(layer, quadpnt)
      
      Kvh = 0
      KvT = 0
      
      do i=1, D
        Kvh(i,i) = Kvh_scalar 
        KvT(i,i) =  KvT_scalar
      end do
      
      
      vct(1:D) =  matmul(-Kvh(1:D,1:D), gradient(1:D)) + matmul(-KvT(1:D,1:D), gradT(1:D))
      
      
       if (present(flux_length)) then
         flux_length = norm2(vct(1:D))
      end if


      if (present(flux)) then
        flux(1:D) = vct(1:D)
      end if
      
    end subroutine vapor_flux
      
      
     

end module evap_RE_fnc
