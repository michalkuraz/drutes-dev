module freeze_fnc
  use pde_objs
  use typy
  use freeze_globs
  use debug_tools
  private :: icerho, watrho, icewatrho, Kliquid_temp, hl, getwater_id, gettempice_id, gettempwat_id 
  public :: iceswitch, capacityhh, capacityhT, diffhh, capacityTh, capacityTT, thetai, Kliquid, surf_tens_deriv


  
  !> double pointer see freeze_pointers for pointing of Kliquid_default, rwcap, theta    
  procedure(tensor_fnc), pointer, public :: Kliquid_default
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
   
   !> Function to calculate the density of water based on ITS-90 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909168/ for water 
    function watrho(quadpnt) result(rho)
      use typy
      use global_objs
    
      type(integpnt_str), intent(in) :: quadpnt
      ! > rho in kg/m^3
      real(kind=rkind) :: rho
      ! > temp in deg C
      real(kind=rkind) :: temp
      
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      
      quadpnt_loc%column = 1
      
      temp = pde(gettempwat_id())%getval(quadpnt_loc)
      
      
      if (temp > 0) then
        rho = (-3.821216e-10*temp*temp*temp*temp+ 6.943248e-8*temp*temp*temp -8.523829e-6*temp*temp + 6.32693e-5*temp + &
               0.99985308)*1e3
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
   
    !> switch for freezing condition based on Clapeyron equation 
    function iceswitch(quadpnt) result(sw)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind) :: sw
      
      real(kind=rkind) :: Tf
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%column = 1
      
      Tf = Tref*exp(pde(getwater_id())%getval(quadpnt_loc)*grav/Lf)
      Tf = Tf - 273.15_rkind

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

      val = hw + iceswitch(quadpnt)*(Lf/grav*log((temp+273.15_rkind)/Tref)) !units

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
  
    function surf_tens_deriv(pde_loc, layer, quadpnt, x) result(val)
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
    
      val = -0.1425-2.38e-4*temp
      
    end function surf_tens_deriv
  
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
        call Kliquid_default(pde_loc, layer, quadpnt, tensor=tensor)
        tensor = icewatrho(quadpnt)*(1-iceswitch(quadpnt))*tensor
      else
        print *, "ERROR! output tensor undefined, exited from diffhh::freeze_fnc"
      end if
    
    
    end subroutine diffhh
    
    subroutine Kliquid(pde_loc, layer, quadpnt, x, tensor, scalar) 
      use typy
      use global_objs
      use globals
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
      
      
      real(kind=rkind), dimension(3,3) :: K
      integer(kind=ikind) :: D
      real(kind=rkind) :: Ei, Ks, thi, thl
      
      
      D = drutes_config%dimen
      
      if (present(quadpnt)) then
        call Kliquid_default(pde_loc, layer, quadpnt, tensor = K(1:D, 1:D))
        thi =  thetai(pde_loc, layer, quadpnt)
        thl = theta(pde_loc, layer, quadpnt)
      else if (present(x)) then
        call Kliquid_default(pde_loc, layer, x=x, tensor=K(1:D, 1:D))
        thi =  thetai(pde_loc, layer, x=x)
        thl = theta(pde_loc, layer, x=x)
      else
        print *, "runtime error"
        print *, "exited from Kliquid::freeze_fnc"
        ERROR STOP
      end if
      !call Kliquid_default(pde_loc, layer, x=(/0.0_rkind/), scalar = Ks)
      
      !Ei = 5/4.0_rkind*(Ks-3)*(Ks-3) + 6 ! no idea what this is 
      
      !tensor = 0
      
      
      
   end subroutine Kliquid
    
  
  
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
      real(kind=rkind), dimension(3,3) :: Klh, Klt
      integer(kind=ikind) :: D

      D = drutes_config%dimen

      if (present(tensor)) then
        call Kliquid_default(pde_loc, layer, quadpnt, tensor = Klt(1:D, 1:D))

        if (present(quadpnt)) then
          tensor = Klt(1:D, 1:D)*gwt*pde(getwater_id())%getval(quadpnt)*surf_tens_deriv(pde_loc, layer, quadpnt)
        else if (present(x)) then
          tensor = Klt(1:D, 1:D)*gwt*pde(getwater_id())%getval(quadpnt)*surf_tens_deriv(pde_loc, layer, x = x)
        else
          print *, "runtime error"
          print *, "exited from Kliquid_temp::freeze_fnc"
          ERROR STOP
        end if
      else
         print *, "ERROR! output tensor undefined, exited from Kliquid_temp::freeze_fnc"
      end if    

      
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
      
      real(kind=rkind), dimension(3,3) :: Klh, Klt, E
      integer(kind=ikind) :: D, i,j
      real(kind=rkind) :: temp
      
      D = drutes_config%dimen
      
      E = 0
      
      temp = pde(gettempwat_id())%getval(quadpnt)
      print*, "bla"
      if (present(tensor)) then
        call Kliquid_temp(pde_loc, layer, quadpnt, tensor=Klt(1:D, 1:D))
        call Kliquid_default(pde_loc, layer, quadpnt, tensor = Klh(1:D, 1:D))
print*, Klt(1:D, 1:D), Klh(1:D, 1:D)
        E(1:D,1:D) = iceswitch(quadpnt) * watrho(quadpnt) * Lf/temp/grav*Klh(1:D,1:D)
        tensor = E(1:D,1:D) - watrho(quadpnt) * Klt(1:D, 1:D)
      else
         print *, "ERROR! output tensor undefined, exited from diffhT::freeze_fnc"
      end if     
      
      
      print*, tensor
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
      val = val - Lf*icerho(quadpnt)*rwcap(pde_loc, layer, x = (/hl(quadpnt)/))*(1-iceswitch(quadpnt))
      
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
      
      
      val = Ci*icerho(quadpnt)*thetai(pde_loc, layer, quadpnt) + Cl*watrho(quadpnt)*theta(pde_loc, layer, quadpnt) 
      
      val = val - Lf*icerho(quadpnt)*rwcap(pde_loc, layer, x = (/hl(quadpnt)/))*iceswitch(quadpnt)*Lf*watrho(quadpnt)/temp 
      
      
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
      
      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen

      
      if (present(tensor)) then
        do i= 1, D
          tensor(i,i) =  thermal_cond
        end do
      end if
      
      
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
          call pde(1)%flux(layer, quadpnt, vector_out = flux)
          flux = Cl *watrho(quadpnt)*flux
        end if
        
        if (present(flux_length)) then
           call pde(1)%flux(layer, quadpnt, scalar = flux_length)
           flux_length = Cl *watrho(quadpnt)*flux_length
        end if
              
    end subroutine convectTT
      
                 
              


end module freeze_fnc
