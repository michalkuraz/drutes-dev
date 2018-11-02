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

!> \file bousfnc.f90
!! \brief Boussinesq equation implementation
!<

!> Solves the problem of saturated flow over a sloping impermeable layer by an extended form of Dupuit-Forchheimer assumption, where the streamlines are parallel to the impermeable layer. The Boussinesq equation states as:
!!\f[
!!\frac{\theta_s}{K_s\textrm{cos} \alpha} \frac{\partial h_b}{\partial t} =  \frac{q(t)}{K} !!\left(1+ \textrm{tan} \alpha \frac{\partial h_b}{\partial x} \right) + \frac{\partial}{\partial !!x} \left( h_b \frac{\partial h_b}{\partial x} \right) - \textrm{tan} \alpha \frac{\partial h_b}!!{\partial x}, 
!!\f]
!! where \f$K_s\f$ is the saturated hydraulic conductivity [\f$L.T^{-1}\f$], \f$\theta_s\f$ is the maximal (saturated) water content (equal to the porosity) [-], \f$q(t)\f$ is the time dependent rainfall [\f$L.T^{-1}\f$], $\alpha$ is the angle formed by the impermeable layer and horizontal plane [\f$rad\f$], and \f$h_b\f$ is  solution of the Boussinesq equation [\f$L\f$] --  it is the vertical distance between the impermeable layer and the groundwater table.
!<

module boussfnc
  public :: bouss_elast
  public :: bouss_cond
  public :: bouss_adv
  public :: darcy4bouss
  public :: boussicond
  public :: bouss_bc
  public :: boussreact
  

    
  contains

    function bouss_elast(pde_loc,layer, quadpnt, x) result(E)
        use typy
        use boussglob
        use pde_objs
        use geom_tools

        class(pde_str), intent(in) :: pde_loc 
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), intent(in), dimension(:),  optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt
        real(kind=rkind) :: h
        !> resulting system elasticity
        real(kind=rkind) :: E
        real(kind=rkind), dimension(1) :: coord
        real(kind=rkind) :: cos_slope, slope
        integer(kind=ikind) :: i

        
        if (present(quadpnt) .and. present(x)) then
          print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
          print *, "exited from re_constitutive::vangen_elast"
          ERROR stop
        else if (.not. present(quadpnt) .and. .not. present(x)) then
          print *, "ERROR: you have not specified either integ point or x value"
          print *, "exited from re_constitutive::vangen_elast"
          ERROR stop
        end if
        
        if (present(quadpnt)) then
          h = pde_loc%getval(quadpnt)
        else
          if (ubound(x,1) /=1) then
            print *, "ERROR: van Genuchten function is a function of a single variable h"
            print *, "       your input data has:", ubound(x,1), "variables"
            ERROR STOP
          end if
          if (ubound(x,1) /=1) then
            print *, "ERROR: van Genuchten function is a function of a single variable h"
            print *, "       your input data has:", ubound(x,1), "variables"
            print *, "exited from re_constitutive::vangen_elast"
            ERROR STOP
          end if
          h = x(1)
        end if

        call getcoor(quadpnt, coord)
        
        do i=1, bouss_slopes(1)%pos - 1
          if (coord(1) >= bouss_slopes(1)%data(i) .and. coord(1) < bouss_slopes(1)%data(i+1)) then
            slope = bouss_slopes(2)%data(i)
            EXIT
          end if
          if (i == bouss_slopes(1)%pos - 1 .and. coord(1) >= bouss_slopes(1)%data(i+1)) then
            slope = bouss_slopes(2)%data(i+1)
          end if
        end do
        cos_slope = cos(atan(slope))
        E = bouss_por/cos_slope 
	

      end function bouss_elast
      
      
    subroutine bouss_cond(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use boussglob
      use pde_objs
      use geom_tools
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      real(kind=rkind) :: h
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar

      real(kind=rkind), dimension(1) :: coord
      integer(kind=ikind) :: i

      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from re_constitutive::mualem"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::mualem"
        ERROR stop
      end if
      
      if (present(quadpnt)) then
        h = pde_loc%getval(quadpnt)
      else
        if (ubound(x,1) /=1) then
          print *, "ERROR: van Genuchten function is a function of a single variable h"
          print *, "       your input data has:", ubound(x,1), "variables"
          print *, "exited from re_constitutive::mualem"
          ERROR STOP
        end if
        h = x(1)
      end if
      
      call getcoor(quadpnt, coord)
      

      
            do i=1, bouss_K(1)%pos - 1
        if (coord(1) >= bouss_K(1)%data(i) .and. coord(1) < bouss_K(1)%data(i+1)) then
          if (present(tensor)) then

            tensor =  bouss_K(2)%data(i) * h

          end if
          if (present(scalar)) then
            scalar =  bouss_K(2)%data(i) * h
          end if
          EXIT
        end if
        if (i == bouss_K(1)%pos - 1 .and. coord(1) >= bouss_K(1)%data(i+1)) then
          if (present(tensor)) then
            tensor =  bouss_K(2)%data(i+1) * h
          end if
          if (present(scalar)) then
            scalar =  bouss_K(2)%data(i+1) * h
          end if
        end if
        
      end do
	
	
      
      end subroutine bouss_cond
      
      
      subroutine bouss_adv(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
          use typy
          use globals
          use boussglob
          use pde_objs
          use geom_tools
          use global_objs
          use debug_tools

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
          real(kind=rkind), dimension(1) :: coord
          real(kind=rkind) :: rain, h, slope
          real(kind=rkind), dimension(1,1) :: K
          integer(kind=ikind) :: i


          
          if (present(quadpnt) .and. present(x)) then
            print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
            print *, "exited from re_constitutive::dmualem_dh"
            ERROR stop
          else if (.not. present(quadpnt) .and. .not. present(x)) then
            print *, "ERROR: you have not specified either integ point or x value"
            print *, "exited from re_constitutive::dmualem_dh"
            ERROR stop
          end if
          
          if (present(quadpnt)) then
            h = pde_loc%getval(quadpnt)
          else
            if (ubound(x,1) /=1) then
              print *, "ERROR: van Genuchten function is a function of a single variable h"
              print *, "       your input data has:", ubound(x,1), "variables"
              print *, "exited from re_constitutive::dmualem_dh"
              ERROR STOP
            end if
            h = x(1)
          end if
          
          
          
          
          do i=1, ubound(bouss_rain,1)-1
            if (i == 1 .and. time < bouss_rain(1)%data(i)) then
              rain = bouss_rain(2)%data(i) 
            end if
            if (bouss_rain(1)%data(i) <= time .and. bouss_rain(1)%data(i+1) > time) then
              rain = bouss_rain(2)%data(i) 
              EXIT
            end if
            if (i == bouss_rain(1)%pos - 1 .and. time >= bouss_rain(1)%data(i+1)) then
              rain = bouss_rain(2)%data(i+1)
            end if
          end do
          
          
          call getcoor(quadpnt, coord)
          
          do i=1, bouss_slopes(1)%pos - 1
            if (coord(1) >= bouss_slopes(1)%data(i) .and. coord(1) < bouss_slopes(1)%data(i+1)) then
              slope = bouss_slopes(2)%data(i)
              EXIT
            end if
            if (i == bouss_slopes(1)%pos - 1 .and. coord(1) >= bouss_slopes(1)%data(i+1)) then
              slope = bouss_slopes(2)%data(i+1)
            end if
          end do
          
          call pde_loc%pde_fnc(1)%dispersion(pde_loc, layer, quadpnt, tensor=K)

          if (present(vector_out)) then
            ! must be negative, because the commnon scheme of the CDE problem has negative convection, but RE has positive convection
            vector_out = -(rain*slope - slope*K(1,1))
          end if

        ! print *, vector_out, rain, slope ; stop

      end subroutine bouss_adv
      
      function boussreact(pde_loc, layer, quadpnt, x) result(react)
        use typy
        use re_globals
        use pde_objs
        use geom_tools
        use boussglob


        class(pde_str), intent(in) :: pde_loc 
        integer(kind=ikind), intent(in) :: layer
        !> pressure head
        real(kind=rkind), intent(in), dimension(:),  optional :: x
        !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
        type(integpnt_str), intent(in), optional :: quadpnt

        !> resulting reaction
        real(kind=rkind) :: react
        
        real(kind=rkind) :: rain
        integer(kind=ikind) :: i
        
        do i=1, ubound(bouss_rain,1)-1
          if (i == 1 .and. time < bouss_rain(1)%data(i)) then
            rain = bouss_rain(2)%data(i) 
          end if
          if (bouss_rain(1)%data(i) <= time .and. bouss_rain(1)%data(i+1) > time) then
            rain = bouss_rain(2)%data(i) 
            EXIT
          end if
          if (i == bouss_rain(1)%pos - 1 .and. time >= bouss_rain(1)%data(i+1)) then
            rain = bouss_rain(2)%data(i+1)
          end if
        end do	
        
        react = -rain
	 
      end function boussreact
      
      
      
      subroutine darcy4bouss(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
        use typy
        use pde_objs
        use global_objs
        use geom_tools
        use boussglob
        
        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in)                          :: layer
        type(integpnt_str), intent(in), optional :: quadpnt    
        real(kind=rkind), intent(in), dimension(:), optional                   :: x
        !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
        real(kind=rkind), dimension(:), intent(in), optional     :: grad
        real(kind=rkind), dimension(:), intent(out), optional    :: flux
        real(kind=rkind), intent(out), optional                  :: flux_length

        real(kind=rkind), dimension(1,1)  :: K
        integer                           :: D
        integer(kind=ikind)               :: i
        real(kind=rkind), dimension(:), allocatable, save  :: gradH
        real(kind=rkind), dimension(:), allocatable, save  :: vct
        real(kind=rkind) :: h, slope
        real(kind=rkind), dimension(1) :: coord
        

        D = drutes_config%dimen
          
        if (.not.(allocated(gradH))) then
          allocate(gradH(1:D))
          allocate(vct(1:D))
        end if

        if (present(quadpnt) .and. (present(grad) .or. present(x))) then
          print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
          ERROR stop
        else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
          print *, "ERROR: you have not specified either integ point or x value"
          print *, "exited from re_constitutive::darcy_law"
          ERROR stop
        end if
        
        if (present(quadpnt)) then
          h = pde_loc%getval(quadpnt)
          call pde_loc%getgrad(quadpnt, gradH)
        else
          if (ubound(x,1) /=1) then
            print *, "ERROR: van Genuchten function is a function of a single variable h"
            print *, "       your input data has:", ubound(x,1), "variables"
            print *, "exited from re_constitutive::darcy_law"
            ERROR STOP
          end if
          h = x(1)
          gradH(1:D) = grad
        end if
        
             call getcoor(quadpnt, coord)
        
        do i=1, bouss_slopes(1)%pos - 1
          if (coord(1) >= bouss_slopes(1)%data(i) .and. coord(1) < bouss_slopes(1)%data(i+1)) then
            slope = bouss_slopes(2)%data(i)
            EXIT
          end if
          if (i == bouss_slopes(1)%pos - 1 .and. coord(1) >= bouss_slopes(1)%data(i+1)) then
            slope = bouss_slopes(2)%data(i+1)
          end if
        end do
        
        
        gradH = gradH + slope
        
        call pde_loc%pde_fnc(1)%dispersion(pde_loc, layer, quadpnt, tensor=K(1:D, 1:D))
            
        
        vct(1:D) =  matmul(-K(1:D,1:D), gradH(1:D))


        if (present(flux_length)) then
          select case(D)
            case(1)
            flux_length = abs(vct(1))
            case(2)
            flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2))
            case(3)
            flux_length = sqrt(vct(1)*vct(1) + vct(2)*vct(2) + vct(3)*vct(3))
          end select
        end if


        if (present(flux)) then
          flux(1:D) = vct(1:D)
        end if

      end subroutine darcy4bouss
      
      subroutine boussicond(pde_loc)
        use typy
        use globals
        use boussglob
        use global_objs
              use pde_objs
              
              class(pde_str), intent(in out) :: pde_loc
              integer(kind=ikind) :: i, j, k, l, m, layer
              real(kind=rkind) :: value
        
             do i=1, elements%kolik
                layer = elements%material(i)
                do j=1, ubound(elements%data,2)
                  k = elements%data(i,j)
                  l = nodes%edge(k)
                  m = pde(1)%permut(k)
                  if (m == 0) then
              call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
                    pde_loc%solution(k) =  value
                  else
              pde_loc%solution(k) = bouss_icond
                  end if
                end do   
              end do
        
      end subroutine boussicond
      
      subroutine bouss_bc(pde_loc, el_id, node_order, value, code) 
        use typy
        use globals
        use global_objs
        use pde_objs
        use geom_tools

        class(pde_str), intent(in) :: pde_loc
        integer(kind=ikind), intent(in)  :: el_id, node_order
        real(kind=rkind), intent(out), optional   :: value
        integer(kind=ikind), intent(out), optional :: code
        
        integer(kind=ikind) :: edge_id, i, j
        
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
          value = pde_loc%bc(edge_id)%series(j,2)
          EXIT
              end if
            end do
          else
            value =  pde_loc%bc(edge_id)%value
          end if
        end if
        
        
        if (present(code)) then
          code = pde_loc%bc(edge_id)%code
        end if
	
      end subroutine bouss_bc


end module boussfnc
