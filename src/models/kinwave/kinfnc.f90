
! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

module kinfnc

  public :: kinconvect, kinbor, kinematixinit, rainfall, kin_elast, getval_kinwave, kinematixinit4cs
  public :: kinflux, kinfluxcl, kincl_lostfactor

  contains 
  
    !> specific function for kinematic wave equation, replaces pde_objs::getvalp1 surface runoff should be in [mm]
    function getval_kinwave(pde_loc, quadpnt) result(val)
      use typy
      use pde_objs
      use geom_tools
      use re_globals
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      real(kind=rkind), dimension(3) :: xyz
      integer(kind=ikind) :: D, layer
      
      val = getvalp1(pde_loc, quadpnt)
           
      if (quadpnt%preproc) then
        val = val*1e3
      end if
      

	
      
    end function getval_kinwave
    
    
    subroutine kinflux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use kinglobs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
    

      real(kind=rkind) :: h, m
      integer(kind=ikind) :: el

      
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        print *, "exited from kinfnc::kinflux"
        ERROR stop
      else if (.not. present(x) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from kinfnc::kinflux"
        ERROR stop
      end if   

      
      if (present(quadpnt)) then
        h = pde_loc%getval(quadpnt)
        if (quadpnt%preproc) then
          h=h*1e-3
        end if
      else
        h = x(1)
      end if
      
      h=max(0.0_rkind, h)
      
      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt", "xypt")
          el = quadpnt%element
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
      end select
      
      m = 5.0_rkind/3

      
      
      if (present(flux)) then
        select case(drutes_config%dimen)
          case(1)
            flux(1) = -1.0_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
          case(2)
            flux(1) = -1.0_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
            flux(2) = -1.0_rkind * sign(1.0_rkind, watershed_el(el)%sy) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
        end select
      end if
          
      
      
      if (present(flux_length)) then
        select case(drutes_config%dimen)
          case(1)     
            flux_length = 1.0_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
          case(2)
            flux_length = norm2((/1.0_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m, &
                            1.0_rkind * sign(1.0_rkind, watershed_el(el)%sy) * & 
                            sqrt(abs( watershed_el(el)%sy))/manning(layer)*h**m/))
        end select          
       end if
       
       
    
    end subroutine kinflux
  
    subroutine kinfluxcl(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use kinglobs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
    

      real(kind=rkind) :: h, m, cl
      integer(kind=ikind) :: el

      
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        print *, "exited from kinfnc::kinflux"
        ERROR stop
      else if (.not. present(x) .and. .not. present(quadpnt)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from kinfnc::kinflux"
        ERROR stop
      end if   

      
      if (present(quadpnt)) then
        h = pde(1)%getval(quadpnt)
        cl = pde(2)%getval(quadpnt)
        if (quadpnt%preproc) then
          h=h*1e-3
        end if
      else
        h = x(1)
        cl = x(2)
      end if
      
      
      select case(quadpnt%type_pnt)
        case("gqnd", "obpt", "xypt")
          el = quadpnt%element
        case("ndpt")
          el = nodes%element(quadpnt%order)%data(1)
      end select
      
      m = 5.0_rkind/3

      
      
      if (present(flux)) then
        select case(drutes_config%dimen)
          case(1)
            flux(1) = -1.0_rkind * cl * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
          case(2)
            flux(1) = -1.0_rkind * cl * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
            flux(2) = -1.0_rkind * cl * sign(1.0_rkind, watershed_el(el)%sy) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
        end select
      end if
          
      
      
      if (present(flux_length)) then
        select case(drutes_config%dimen)
          case(1)     
            flux_length = 1.0_rkind * cl * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m
          case(2)
            flux_length = norm2((/cl*1.0_rkind * sign(1.0_rkind, watershed_el(el)%sx) * & 
                            sqrt(abs( watershed_el(el)%sx))/manning(layer)*h**m, &
                            cl*1.0_rkind * sign(1.0_rkind, watershed_el(el)%sy) * & 
                            sqrt(abs( watershed_el(el)%sy))/manning(layer)*h**m/))
        end select          
       end if
       
       
    
    end subroutine kinfluxcl
    
    
  
    subroutine kinconvect(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use geom_tools
      use debug_tools
      use kinglobs
      use core_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
      real(kind=rkind), dimension(:), intent(in), optional :: vector_in
      !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
      !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
      !<
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
      real(kind=rkind), intent(out), optional :: scalar
      
      integer(kind=ikind) :: el, D, i
      real(kind=rkind) :: hsurf, m
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind), dimension(:), allocatable, save :: ndvals, slopes, outflux
      real(kind=rkind), dimension(:,:), allocatable, save :: abc
      
      
      el = quadpnt%element
      
      hsurf = max(0.0_rkind, pde_loc%getval(quadpnt))
      
      hsurf = min(5e-3_rkind, hsurf)
      
      m = 5.0_rkind/3
      
      D = drutes_config%dimen
      
      if (.not. allocated(slopes)) allocate(slopes(drutes_config%dimen))
      if (.not. allocated(outflux)) allocate(outflux(drutes_config%dimen))
      
      if (backwater) then
        if (.not. allocated(ndvals)) allocate(ndvals(ubound(elements%data,2)))
        
        quadpnt_loc%type_pnt = "ndpt"
        quadpnt_loc%column = 3
        
        do i=1, ubound(ndvals,1)
          quadpnt_loc%order = elements%data(el,i)
          ndvals(i) = pde_loc%getval(quadpnt_loc)
        end do
        
        select case(D)
          case(1)
            slopes(1) = (ndvals(2) - ndvals(1))/(nodes%data(elements%data(el,2),1) - nodes%data(elements%data(el,1),1))
            slopes(1) = slopes(1) +  watershed_el(el)%sx
          case(2)
            if (.not. allocated(abc)) allocate(abc(ubound(elements%data,2),3))
            
            do i=1, ubound(elements%data,2)
              abc(i, 1:2) = nodes%data(elements%data(el,i),:)
              abc(i,3) = ndvals(i)
            end do
            
            call plane_derivative(abc(1,:), abc(2,:), abc(3,:), slopes(1), slopes(2))
            
            slopes(1) = slopes(1) + watershed_el(el)%sx
            slopes(2) = slopes(2) + watershed_el(el)%sy
            
          case(3)
            print *, "kinematic wave has no sense for three-dimensions"
            print *, "exited from kinfnc::kinconvect"
            ERROR STOP
          end select
        else
          select case(D)
            case(1)
              slopes(1) = watershed_el(el)%sx
            case(2)
              slopes(1) = watershed_el(el)%sx
              slopes(2) = watershed_el(el)%sy
            case(3)
            print *, "kinematic wave has no sense for three-dimensions"
            print *, "exited from kinfnc::kinconvect"
            ERROR STOP
          end select
        end if
        
        
        
        outflux = -1* sign(1.0_rkind, slopes(1:D)) * & 
                            sqrt(abs( slopes(1:D)))/manning(layer)*m*hsurf**(m-1)
                            
        if (present(vector_out)) vector_out(1:D) = outflux
                            
        if (present(scalar)) scalar = norm2(outflux)
  

        
    end subroutine kinconvect
    
    
     subroutine kinconvectcl(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use geom_tools
      use debug_tools
      use kinglobs
      use core_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> this argument is required by the global vector_fnc procedure pointer, unused in this procedure
      real(kind=rkind), dimension(:), intent(in), optional :: vector_in
      !> first order tensor of the unsaturated hydraulic conductivity derivative in respect to h. it is the last column of the hydraulic conductivity second order tensor times  
      !!relative unsaturated hydraulic conductivity derivative in respect to h (scalar value)
      !<
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative unsaturated hydraulic conductivity derivative in respect to h, scalar value
      real(kind=rkind), intent(out), optional :: scalar
      
      integer(kind=ikind) :: el, D, i
      real(kind=rkind) :: hsurf, m, cl
      type(integpnt_str) :: quadpnt_loc
      real(kind=rkind), dimension(:), allocatable, save :: ndvals, slopes, outflux
      real(kind=rkind), dimension(:,:), allocatable, save :: abc
      
      
      el = quadpnt%element
      
      hsurf = max(0.0_rkind, pde(1)%getval(quadpnt))
!      cl = max(0.0_rkind, pde(2)%getval(quadpnt))

      
      m = 5.0_rkind/3
      
      D = drutes_config%dimen
      
      if (.not. allocated(slopes)) allocate(slopes(drutes_config%dimen))
      if (.not. allocated(outflux)) allocate(outflux(drutes_config%dimen))
      
      if (backwater) then
        if (.not. allocated(ndvals)) allocate(ndvals(ubound(elements%data,2)))
        
        quadpnt_loc%type_pnt = "ndpt"
        quadpnt_loc%column = 3
        
        do i=1, ubound(ndvals,1)
          quadpnt_loc%order = elements%data(el,i)
          ndvals(i) = pde(1)%getval(quadpnt_loc)
        end do
        
        select case(D)
          case(1)
            slopes(1) = (ndvals(2) - ndvals(1))/(nodes%data(elements%data(el,2),1) - nodes%data(elements%data(el,1),1))
            slopes(1) = slopes(1) +  watershed_el(el)%sx
          case(2)
            if (.not. allocated(abc)) allocate(abc(ubound(elements%data,2),3))
            
            do i=1, ubound(elements%data,2)
              abc(i, 1:2) = nodes%data(elements%data(el,i),:)
              abc(i,3) = ndvals(i)
            end do
            
            call plane_derivative(abc(1,:), abc(2,:), abc(3,:), slopes(1), slopes(2))
            
            slopes(1) = slopes(1) + watershed_el(el)%sx
            slopes(2) = slopes(2) + watershed_el(el)%sy
            
          case(3)
            print *, "kinematic wave has no sense for three-dimensions"
            print *, "exited from kinfnc::kinconvect"
            ERROR STOP
          end select
        else
          select case(D)
            case(1)
              slopes(1) = watershed_el(el)%sx
            case(2)
              slopes(1) = watershed_el(el)%sx
              slopes(2) = watershed_el(el)%sy
            case(3)
            print *, "kinematic wave has no sense for three-dimensions"
            print *, "exited from kinfnc::kinconvect"
            ERROR STOP
          end select
        end if
        
        outflux = (1.0_rkind) * sign(1.0_rkind, slopes(1:D)) * & 
                            sqrt(abs( slopes(1:D)))/manning(layer)*hsurf**m
                            

        if (present(vector_out)) vector_out(1:D) = outflux
                            
        if (present(scalar)) scalar = norm2(outflux)
        
    end subroutine kinconvectcl
    
    
    
    subroutine kinbor(pde_loc, el_id, node_order, value, code, array, bcpts) 
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
      
      

      if (present(value)) then
        value = 0.0_rkind
      end if

      
      if (present(code)) then
        code = 1
      end if
  
    end subroutine kinbor
    
    
   subroutine kinborO(pde_loc, el_id, node_order, value, code, array, bcpts) 
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
      
      

      if (present(value)) then
        value = 0.0_rkind
      end if

      
      if (present(code)) then
        code = 0
      end if
  
    end subroutine kinborO
    
    
   subroutine kinborcs(pde_loc, el_id, node_order, value, code, array) 
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

      
      

      if (present(value)) then
        value = 0.0_rkind
      end if

      
      if (present(code)) then
        code = 0
      end if
  
    end subroutine kinborcs
      
    subroutine kinematixinit(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs

      class(pde_str), intent(in out) :: pde_loc

      pde_loc%solution(:) = 0.0_rkind
     
      
    end subroutine kinematixinit
    
    
   subroutine kinematixinit4cs(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use kinglobs
      use debug_tools
      use readtools

      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, el, layer
      integer :: file_id, ierr
      real(kind=rkind) :: tmp

      if (.not. icond_file) then
        do i=1, ubound(pde_loc%solution,1)
          el = nodes%element(i)%data(1)
          layer = elements%material(el)
          pde_loc%solution(i) = kinsols(layer)%csinit
        end do
      else
        open(newunit=file_id, file="drutes.conf/kinwave/initcond", action="read", status="old", iostat=ierr)
        if (ierr /= 0) then
          print *, "in file drutes.conf/kinwave/kinwave.conf you have defined to use node-wise data for initial c_s, and so"
          print *, "file drutes.conf/kinwave/initcond must be provided"
          print *, "ERROR: unable to open file drutes.conf/kinwave/initcond"
          ERROR STOP
        end if
        do i=1, nodes%kolik
          call fileread(pde_loc%solution(i), file_id, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
                errmsg="Are your data in drutes.conf/kinwave/initcond consistent with your mesh data?" )
        end do
        read(unit=file_id, fmt=*, iostat=ierr) tmp
        if (ierr == 0) then
          call file_error(file_id, errmsg="you have more values in drutes.conf/kinwave/initcond than the umber of nodes, & 
                Are your data in drutes.conf/kinwave/initcond consistent with your mesh data?")
        end if
      end if
      
          
   
      
    end subroutine kinematixinit4cs
    
    !> for model Ks
    !! \f[ q_{in} = K_s \f]
    !! for model Swartzendruber cumulative infiltration is given by
    !! \f[ I = \frac{S(1-\exp(-A \sqrt(t)))}{A} + K_s t \f]
    !! the infiltration speed is given by \f[ \dv{I}{t}, \f] and so
    !! \f[  q_{in} =  \dv{\frac{S(1-\exp(-A \sqrt(t)))}{A} + K_s t}{t} = \frac{S\exp(-A\sqrt{t})}{2\sqrt{t}}+K \f]
    !! surface runoff is then given by
    !! \f[ i = \mathrm{max}(0, q_{in})
    !<
    function rainfall(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use kinglobs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind) :: val
      real(kind=rkind) :: qin, S, A, K 
      integer(kind=ikind), save :: position = 1
      integer(kind=ikind) :: i
      
      
       if (position < ubound(raindata(1)%series,1)) then
        do i=position, ubound(raindata(1)%series,1)-1
           if (raindata(1)%series(i,1) < time .and. raindata(1)%series(i+1,1) > time) then
            position = i
            EXIT
          else if (raindata(1)%series(i+1,1) < time .and. i == ubound(raindata(1)%series,1) ) then
            position = i
          end if
        end do
      end if
    
      val = raindata(el2pt(quadpnt%element))%series(position,2)

      select case(inf_model(layer)%name)
        case("Ks")
          qin = inf_model(layer)%Ks
        case("Schwarz")
          S = inf_model(layer)%S
          A = inf_model(layer)%A
          K = inf_model(layer)%Ks
          
          if (time > epsilon(time)) then
            qin = S/(2*sqrt(time))*exp(-A*sqrt(time)) + K
          else
            qin = huge(qin)
          end if
      end select
      
      val = max(0.0_rkind, val - qin)
          
      

    
    end function rainfall
    
    
    
    function kincl_source(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use kinglobs
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      real(kind=rkind) :: h
      
      h = pde(1)%getval(quadpnt)
      
      val = h*kinsols(layer)%lambda_a
      
      
    end function kincl_source
    
    
    function kincl_lostfactor(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use kinglobs
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      real(kind=rkind) :: cs
      
      cs = max(.0_rkind,pde(3)%getval(quadpnt))
      
      val = kinsols(layer)%lost_me
      
      
    end function kincl_lostfactor
    
    function kincs_source(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use kinglobs
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
      
      real(kind=rkind) :: cs 
      
      
      cs = max(0.0_rkind, pde(3)%getval(quadpnt))
      
      
      val = -kinsols(layer)%horb*kinsols(layer)%rhos*kinsols(layer)%lambda_d * &
              (cs*kinsols(layer)%horb*kinsols(layer)%rhos)**(kinsols(layer)%ncs-1)
      

      
    end function kincs_source
    
    
    
    function kin_elast(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use kinglobs
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      end if


      E = 1.0_rkind
      

    end function kin_elast 
    
    
        
    function kin_clelast(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use heat_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      end if


      E = max(1e-7_rkind, pde(1)%getval(quadpnt))
!E = 1e-3
      

    end function kin_clelast
    
    function solmass(pde_loc,layer, quadpnt, x) result(M)
      use typy
      use heat_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> resulting system elasticity
      real(kind=rkind) :: M
      
      real(kind=rkind) :: h, cl

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      end if
      
      h = max(0.0_rkind, pde(1)%getval(quadpnt))
      
      cl = max(0.0_rkind, pde(2)%getval(quadpnt))
      
      if (quadpnt%preproc) h = h*1e-3

      M = h*cl
      

    end function solmass
    
    
    function csmass(pde_loc,layer, quadpnt, x) result(M)
      use typy
      use heat_globals
      use pde_objs
      use core_tools
      use kinglobs

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> resulting system elasticity
      real(kind=rkind) :: M
      
      real(kind=rkind) :: h, cs

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      end if
      
      h = kinsols(layer)%horb
      
      cs = max(0.0_rkind, pde(3)%getval(quadpnt))

      M = h*cs
      

    end function csmass
    
    
    function kin_cselast(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use kinglobs
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_elast"
        ERROR stop
      end if


      E = kinsols(layer)%horb*kinsols(layer)%rhos
      
      

    end function kin_cselast
    
        
!    function kin_csclelast(pde_loc,layer, quadpnt, x) result(E)
!      use typy
!      use heat_globals
!      use pde_objs
!      use core_tools

!      class(pde_str), intent(in) :: pde_loc 
!      integer(kind=ikind), intent(in) :: layer
!      !> pressure head
!      real(kind=rkind), intent(in), dimension(:),  optional :: x
!      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
!      type(integpnt_str), intent(in), optional :: quadpnt
!      real(kind=rkind) :: h
!      !> resulting system elasticity
!      real(kind=rkind) :: E

       
      
!      if (present(quadpnt) .and. present(x)) then
!        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
!        print *, "exited from heat_fnc::heat_elast"
!        ERROR stop
!      else if (.not. present(quadpnt) .and. .not. present(x)) then
!        print *, "ERROR: you have not specified either integ point or x value"
!        print *, "exited from heat_fnc::heat_elast"
!        ERROR stop
!      end if


!      E = pde(1)%getval(quadpnt)
      

!    end function kin_csclelast

    subroutine disp4kinwave(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use geom_tools
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
      

      integer(kind=ikind) :: D, i, nd1, nd2, el
      real(kind=rkind) :: height, convect
      
     
      D = drutes_config%dimen
      el = quadpnt%element
      height = 0
      do i=1, ubound(elements%data,2)
          
        if (i < ubound(elements%data,2)) then
          nd1 = elements%data(el,i)
          nd2 = elements%data(el,i+1)
        else
          nd1 = elements%data(el,i)
          nd2 = elements%data(el,1)
        end if
        
        height = max(height, dist(nodes%data(nd1,:), nodes%data(nd2,:)))
        
      end do
      
      call kinconvect(pde(1), layer, quadpnt, scalar=convect)
      
      if (present(tensor)) then
        tensor =  convect*height
      end if
      
    
    end subroutine disp4kinwave
    
    
    subroutine disp4kinwavecl(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use geom_tools
      use debug_tools
      use kinglobs
      
      
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

      integer(kind=ikind) :: D, i, nd1, nd2, el
      real(kind=rkind) :: height, convect
      
     
      D = drutes_config%dimen
      el = quadpnt%element
      height = 0
      do i=1, ubound(elements%data,2)
          
        if (i < ubound(elements%data,2)) then
          nd1 = elements%data(el,i)
          nd2 = elements%data(el,i+1)
        else
          nd1 = elements%data(el,i)
          nd2 = elements%data(el,1)
        end if
        
        height = max(height, dist(nodes%data(nd1,:), nodes%data(nd2,:)))
        
      end do
      
      call kinconvectcl(pde(2), layer, quadpnt, scalar=convect)
      
      if (present(tensor)) then
        tensor = 0
        do i=1, drutes_config%dimen
          tensor(i,i) = convect*kinsols(layer)%diff
        end do
      end if
    end subroutine disp4kinwavecl


end module kinfnc
