module freeze_fnc
  use pde_objs
  use typy
  use freeze_globs
  private :: icerho, watrho, icewatrho, Kliquid_temp, hl, getwater_id, gettempice_id, gettempwat_id
  public :: iceswitch, capacityhh, capacityhT, diffhh, capacityTh, capacityTT, thetai


  
  procedure(tensor_fnc), pointer, public :: Kliquid
  procedure(scalar_fnc), pointer, public :: rwcap, theta
  
  contains

  
    pure function gettempice_id() result(id)
      use typy
      
      integer(kind=ikind) :: id
      
      select case(drutes_config%name)
        case("freeze")
          id = 2
        case("LTNE")
          !will be editted later, temporaly genreates error, so we won't forget to update it
          id = -1
      end select
      
    end function gettempice_id
  
  
    
    pure function gettempwat_id() result(id)
      use typy
      
      integer(kind=ikind) :: id
      
      select case(drutes_config%name)
        case("freeze")
          id = 2
        case("LTNE")
          !will be editted later, temporaly genreates error, so we won't forget to update it
          id = -1
      end select
      
    end function gettempwat_id
    
    
    pure function getwater_id() result(id)
      use typy
      
       integer(kind=ikind) :: id
       
       id = 1
       
     end function getwater_id
    
      
  

    function icerho(quadpnt) result(rho)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: rho
      real(kind=rkind) :: temp
      
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      
      quadpnt_loc%column = 1
      
      temp = pde(gettempice_id())%getval(quadpnt_loc)
      
      
      rho = exp(log(999.946997406686) + 5e-5*temp)

      
      
    end function icerho
   
   
    function watrho(quadpnt) result(rho)
      use typy
      use global_objs
    
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: rho
      
      real(kind=rkind) :: temp
      
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      
      quadpnt_loc%column = 1
      
      temp = pde(gettempwat_id())%getval(quadpnt_loc)
      
      if (temp>0) then
        rho = (1.682208e-8*temp*temp*temp - 6.05282462e-6*temp*temp + 2.36680033177935e-5*temp + &
               0.999946997406686)*1e3
      else
        rho = exp(log(999.946997406686) + 5e-5*temp)
      end if
        
    end function watrho
   
   
   
    function icewatrho(quadpnt) result(rho)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: rho
      
      integer(kind=ikind) :: layer, el
      real(kind=rkind) :: thl, thall
      
      
      if (quadpnt%type_pnt == "ndpt" ) then
        el = nodes%element(quadpnt%order)%data(1)
      else
        el = quadpnt%element
      end if
      
      layer = elements%material(el)
      
      thl = theta(pde(getwater_id()), layer, x=(/hl(quadpnt)/))
      thall = theta(pde(getwater_id()), layer, quadpnt)
      
      rho = (thl * watrho(quadpnt) + thall * icerho(quadpnt) - thl * icerho(quadpnt))/thall
        
     end function icewatrho
   
   
    function iceswitch(quadpnt) result(sw)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind) :: sw
      
      real(kind=rkind) :: Tf
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%column = 1
      
      Tf = Tref*exp(pde(getwater_id())%getval(quadpnt_loc)*grav/Lf) - Tref
      
      if (pde(gettempwat_id())%getval(quadpnt_loc) > Tf) then
        sw = 0
      else
        sw = 1
      end if
          
    end function iceswitch
   
    function hl(quadpnt) result(val)
      use typy
      use global_objs
      use freeze_globs
     
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind) :: hw, temp
      
      hw = pde(getwater_id())%getval(quadpnt)
      
      temp = pde(gettempwat_id())%getval(quadpnt)
      
      val = hw + iceswitch(quadpnt)*(Lf/grav*log((temp+Tref)/Tref) - hw)
     
    end function hl
    
    
    function thetai(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
     
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      
      real(kind=rkind) :: thl, thall
      
      
      thl = theta(pde(getwater_id()), layer, x=(/hl(quadpnt)/))
      thall = theta(pde(getwater_id()), layer, quadpnt)
      
      val = (thall * icewatrho(quadpnt) - thl * watrho(quadpnt))/icerho(quadpnt)
     
    end function thetai
    
   
   
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
    
    
      val = (1-iceswitch(quadpnt))*watrho(quadpnt)*rwcap(pde_loc, layer, x=(/hl(quadpnt)/))
    
      val = val + icerho(quadpnt)*(rwcap(pde_loc, layer, quadpnt) - &
            rwcap(pde_loc, layer, x=(/hl(quadpnt)/))*(1-iceswitch(quadpnt)))
          
      val = val * icewatrho(quadpnt) * grav 
      
   
    end function capacityhh
   
   
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
    
      real(kind=rkind) :: temp
    
      temp = pde(gettempwat_id())%getval(quadpnt)
    
    
      val = iceswitch(quadpnt) * rwcap(pde_loc, layer, x=(/hl(quadpnt)/)) * watrho(quadpnt)*watrho(quadpnt) * Lf/temp
    
      val = val - iceswitch(quadpnt) * rwcap(pde_loc, layer, x=(/hl(quadpnt)/)) *watrho(quadpnt)*icerho(quadpnt) * Lf/temp
    
      
    end function capacityhT
  
  
    subroutine diffhh(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
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
        call Kliquid(pde_loc, layer, quadpnt, tensor=tensor)
        tensor = icewatrho(quadpnt)*(1-iceswitch(quadpnt))*tensor
      else
        print *, "ERROR! output tensor undefined, exited from diffhh::freeze_fnc"
      end if
    
    
    end subroutine diffhh
  
  
    subroutine Kliquid_temp(pde_loc, layer, quadpnt, x, tensor, scalar) 
      use typy
      use global_objs
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
    
    
    end subroutine Kliquid_temp
  
      
      
      
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
      
      real(kind=rkind), dimension(3,3) :: Klt, E
      integer(kind=ikind) :: D, i,j
      real(kind=rkind) :: temp
      
      D = drutes_config%dimen
      
      E=0
      
      temp = pde(gettempwat_id())%getval(quadpnt)
      
      if (present(tensor)) then
        call Kliquid_temp(pde_loc, layer, quadpnt, tensor=Klt(1:D, 1:D))
        E(1,1) = iceswitch(quadpnt) * watrho(quadpnt) * watrho(quadpnt) * Lf/temp
        do i=2, D
          E(i,i) = E(1,1)
        end do
        tensor = E(1:D,1:D) - watrho(quadpnt) * Klt(1:D, 1:D)
      else
         print *, "ERROR! output tensor undefined, exited from diffhT::freeze_fnc"
      end if     
      
      
      
    end subroutine diffhT


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
      

      
      val = Lf*icerho(quadpnt)*rwcap(pde_loc, layer, quadpnt)
      val = val - Lf*icerho(quadpnt)*rwcap(pde_loc, layer, x=(/hl(quadpnt)/))*(1-iceswitch(quadpnt))
      
      val = val * icewatrho(quadpnt) * grav
      
      
    end function capacityTh
    
    
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
      
      real(kind=rkind) :: temp
      
      temp = pde(gettempwat_id())%getval(quadpnt)
      
      
      val = Ci*thetai(pde_loc, layer, quadpnt) + Cl*theta(pde_loc, layer, quadpnt) 
      
      val = val - Lf*icerho(quadpnt)*rwcap(pde_loc, layer, x=(/hl(quadpnt)/))*iceswitch(quadpnt)*Lf*watrho(quadpnt)/temp 
      
      
    end function capacityTT
    
    
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
      
      real(kind=rkind), dimension(3,3) :: Klt, E
      integer(kind=ikind) :: D, i,j
      
        
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
      
      call darcy_frozen(pde_loc, layer, quadpnt, flux)
      
      flux = flux * Cl
      
      
    end subroutine convectTT
      
              
    subroutine darcy_frozen(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
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

      real(kind=rkind), dimension(3,3)  :: K
      integer                           :: D
      integer(kind=ikind)               :: i
      integer(kind=ikind), dimension(3) :: nablaz
      real(kind=rkind), dimension(3)  :: gradH
      real(kind=rkind), dimension(3)  :: vct
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable :: gradient
      type(integpnt_str) :: quadpnt_loc
      

      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::darcy_law"
        ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
        quadpnt_loc%preproc=.true.
        h = pde_loc%getval(quadpnt_loc)
        call pde_loc%getgrad(quadpnt, gradient)
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

      nablaz = 0
      nablaz(D) = 1
      
      gradH(1:D) = gradient(1:D) + nablaz(1:D)

      call Kliquid(pde_loc, layer, x=(/h/), tensor=K(1:D, 1:D))
     
      
      vct(1:D) = matmul(-K(1:D,1:D), gradH(1:D))


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

    end subroutine darcy_frozen     

    
              


end module freeze_fnc
