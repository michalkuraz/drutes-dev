module lsconstitutive

  contains 
  
  
    
  
    subroutine ncflux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use globals 
      use global_objs
      use pde_objs
      use geom_tools
      use debug_tools
      use datetime
      use ncglobvars
      use core_tools
      use netcdfflux
      use ncfluxarea
      use init_netcdf
    
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
      
      real(kind=rkind), dimension(2) :: xy, gradsl
      integer(kind=ikind) :: nowhrs, el, ncel
      real(kind=rkind) :: tmp, q, Wcell, Acell
      logical :: success
      character(len=1024) :: errmsg
      
      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt")
          if (.not. ncfluxdata%activeel(quadpnt%element)) then
            if (present(flux)) flux = 0
            if (present(flux_length)) flux_length = 0
            RETURN
          end if
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
            if (.not. ncfluxdata%activeel(el) )then
              if (present(flux)) flux = 0
              if (present(flux_length)) flux_length = 0
            RETURN
          end if
          CONTINUE
      end select
      
 
      
      call getcoor(quadpnt, xy)
      
      
      nowhrs = ora_di_ini + int(time/86400.0_rkind)*24
      

      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt")
          el = quadpnt%element


        case("numb") 
          print *, "unable to print convection value for quadpnt%type_pnt = numb "
          print *, "exited from lsconstitutive::ncflux"
          ERROR STOP
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
          
        case default
          print *, "incorrect quadpnt%type_pnt: ", quadpnt%type_pnt
          print *, "exited from lsconstitutive::ncflux"
          ERROR STOP
      end select
      
      
      ncel = el2ncgrid(el)      
      Acell = ncelements%areas(ncel)
      Wcell = sqrt(Acell)
      
      gradsl = ncelements%ders(ncel,:,1)
      

      
      if (norm2(gradsl) < 100*epsilon(tmp)) then  
        gradsl = 1.0_rkind/sqrt(2.0_rkind)
      else
        gradsl = gradsl/norm2(gradsl)
      end if
      
!      gradsl = [0.0_rkind, -1.0_rkind]

      call ncflux_get_xy(xy(1), xy(2),  nowhrs, q, success, errmsg)
      
      
      if (.not. success) then
        q=0.0_rkind
      end if
      
      if (present(flux)) then
        flux = -q*gradsl/Wcell
!flux = gradsl
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(-q*gradsl/Wcell)
      end if
      
  
    end subroutine ncflux
    
    function velocity(Q, w) result(v)
      use ncglobvars
      use typy
      
      real(kind=rkind), intent(in) :: Q, w
      real(kind=rkind) :: v
      
      real(kind=rkind) :: k, m=3.0_rkind/5.0_rkind
      
      k = vref*w**(1-m)/Qref**(1-m)
      
      v = k*Q**(1-m)/w**(1-m)
    
    end function velocity
    
    function heff(Q, w, v) result(h)
      use typy
      real(kind=rkind), intent(in) :: Q, w, v
      real(kind=rkind) :: h
      
      h = Q/(w*v)
      
    
    end function heff
      
    
    subroutine ADEls_convection(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> input vector
      real(kind=rkind), dimension(:), intent(in), optional  :: vector_in
      !> output vector
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional               :: scalar
    
     
      if (present(vector_out)) then
        call pde(1)%flux(layer, quadpnt, vector_out=vector_out)
      end if
        
      if (present(scalar)) then
        call pde(1)%flux(layer, quadpnt, scalar=scalar)
      end if
      
    end subroutine ADEls_convection
    
    
    function ADEls_tder_coef(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      use ncglobvars
      use netcdfflux
      use geom_tools
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind), dimension(2) :: xy
      integer(kind=ikind) :: nowhrs, ncell, el, ncel
          
      logical :: success
      character(len=1024) :: errmsg
      real(kind=rkind) :: q, v, Wcell, Acell
      
      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt")
          if (.not. ncfluxdata%activeel(quadpnt%element)) then
            val = 1
            RETURN
          end if
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
            if (.not. ncfluxdata%activeel(el) )then
              val = 1
              RETURN
            end if
      end select
      
      call getcoor(quadpnt, xy)
      
      nowhrs = ora_di_ini + int(time/86400.0_rkind)*24
      
      call ncflux_get_xy(xy(1), xy(2),  nowhrs, q, success, errmsg)
      
      if (q < 0) then
        val = 1
        RETURN
      end if
      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt")
          el = quadpnt%element

        case("numb") 
          print *, "unable to print convection value for quadpnt%type_pnt = numb "
          print *, "exited from lsconstitutive::ncflux"
          ERROR STOP
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
        case default
          print *, "incorrect quadpnt%type_pnt: ", quadpnt%type_pnt
          print *, "exited from lsconstitutive::ncflux"
          ERROR STOP
      end select
      
      
      ncel = el2ncgrid(el)
      
     
      Acell = ncelements%areas(ncel)
      Wcell = sqrt(Acell)
      
      v = velocity(q,Wcell)
      
      
      val = q/(Wcell*v)
      
      if (isnan(val)) then
        print *, Acell, Wcell, q, v
        stop
      end if
      
              
    end function ADEls_tder_coef
    
    
    function ADEls_mass(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      use debug_tools
      use netcdfflux
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      real(kind=rkind)                :: q, c
      integer(kind=ikind) :: el
      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt")
          el = quadpnt%element
        case("numb") 
          print *, "unable to print convection value for quadpnt%type_pnt = numb "
          print *, "exited from lsconstitutive::ncflux"
          ERROR STOP
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
        case default
          print *, "incorrect quadpnt%type_pnt",  quadpnt%type_pnt
          print *, "exited from lsconstitutive::ncflux"
          ERROR STOP
      end select
      
      c = pde(1)%getval(quadpnt)
      
      call pde(1)%flux(layer, quadpnt, scalar=q)
      
      val = c*q*elements%areas(el)
      

    end function ADEls_mass
    
    
    subroutine ADEls_icond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use ncglobvars

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value

      pde_loc%solution(:) = cinit_ls

    end subroutine ADEls_icond
    
    
    subroutine ADEls_dirichlet(pde_loc, el_id, node_order, value, code, array, bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts

      
      integer(kind=ikind) :: edge_id, i, j, proc
      real(kind=rkind) :: tempval
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      if (present(value)) then
        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              tempval = pde_loc%bc(edge_id)%series(j,2)
              EXIT
            end if
          end do
        else
          tempval =  pde_loc%bc(edge_id)%value
        end if
        value = tempval 
      end if

      
      if (present(code)) then
        code = 1
      end if
      

    end subroutine adeLS_dirichlet
    
    
    subroutine ADElsdisp(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use ADE_globals
      use re_globals
      use debug_tools
      use ncglobvars
      
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
      real(kind=rkind) :: q
      real(kind=rkind), dimension(2,2) :: identity
      
     
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from lsconstitutive::ADElsdisp"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from lsconstitutive::ADElsdisp"
        ERROR stop
      end if
     
      D = drutes_config%dimen
      
      identity = 0.0_rkind
      do i=1, D
        identity(i,i) = 1.0_rkind
      end do
      
      call pde(1)%flux(layer, quadpnt, scalar=q)
      
      
      if (present(tensor)) then
        tensor = identity*LSdisp*q
      end if
      
      if (present(scalar)) then
        scalar = LSdisp*q
      end if
 
    
    end subroutine ADElsdisp


end module lsconstitutive
