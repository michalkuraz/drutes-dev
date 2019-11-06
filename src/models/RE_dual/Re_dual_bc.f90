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

!> \file Re_dual_bc.f90
!!  \brief This module contains different boundary definitions for the fracture and matrix domains in the dual permeability model
!<

module Re_dual_bc
  public :: dual_neumann_bc
  public :: dual_freedrainage
  public :: dual_atmospheric
  public :: inf_neumann_bc
  contains
  
   !> Defines Neumann (flux) boundary condition for the dual permeability model
   !< The Neumann boundary condition is assigned weighted according to area weights.
  subroutine dual_neumann_bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals
      !> takes the global pde_loc structure
      class(pde_str), intent(in) :: pde_loc 
      !> takes the element id and node order
      integer(kind=ikind), intent(in)  :: el_id, node_order
      !> optional value output of the assigned boundary 
      real(kind=rkind), intent(out), optional    :: value
      !> optional code output. Different codes for different boundaries. Code = 2 
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
     

      integer(kind=ikind) :: i, edge_id, j,layer
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, weight
      integer :: i1
      
      

      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))

        i = pde_loc%permut(elements%data(el_id, node_order))
        

        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              layer=elements%material(el_id)
	       
              select case (pde_loc%mfswitch)
                case("m")
                  bcval = pde_loc%bc(edge_id)%series(j,2)*exchange(layer)%weightm
                case("f")
                  bcval = pde_loc%bc(edge_id)%series(j,2)*exchange(layer)%weightf
              end select
              EXIT
            end if
          end do
        else
          layer=elements%material(el_id)
          select case (pde_loc%mfswitch)
            case("m")
              bcval = pde_loc%bc(edge_id)%value*exchange(layer)%weightm
            case("f")
              bcval = pde_loc%bc(edge_id)%value*exchange(layer)%weightf
          end select
        end if
    

        value = bcval

      end if
      
      if (present(code)) then
        code = 2
      end if


    end subroutine dual_neumann_bc
    
   !> Also defines Neumann type boundary. The weighting for matrix and fracture domains is user defined. 
  subroutine inf_neumann_bc(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals

      !> takes the global pde_loc structure
      class(pde_str), intent(in) :: pde_loc 
      !> takes the element id and node order
      integer(kind=ikind), intent(in)  :: el_id, node_order
      !> optional value output of the assigned boundary 
      real(kind=rkind), intent(out), optional    :: value
      !> optional code output. Different codes for different boundaries. Code = 2 
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array

      integer(kind=ikind) :: i, edge_id, j,layer
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, weight
      integer :: i1
      
      select case (pde_loc%mfswitch)
        case("m")
          weight=(1_rkind-infweight)
        case("f")
          weight=infweight
      end select

      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))

        i = pde_loc%permut(elements%data(el_id, node_order))
        

        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              layer=elements%material(el_id)
              bcval = pde_loc%bc(edge_id)%series(j,2)*weight
              EXIT
            end if
          end do
        else
          layer=elements%material(el_id)
          bcval = pde_loc%bc(edge_id)%value*weight
        end if

        value = bcval

      end if
      
      if (present(code)) then
        code = 2
      end if


    end subroutine inf_neumann_bc
    
  !> Assigns free drainage boundary (unit gradient) for dual permeability model
  subroutine dual_freedrainage(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals
      use dual_por
      
      !> takes the global pde_loc structure
      class(pde_str), intent(in) :: pde_loc 
      !> takes the element id and node order
      integer(kind=ikind), intent(in)  :: el_id, node_order
      !> optional value output of the assigned boundary 
      real(kind=rkind), intent(out), optional    :: value
      !> optional code output. Different codes for different boundaries. Code = 2 
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      
      
      real(kind=rkind), dimension(3,3) :: K
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer, D
      real(kind=rkind), dimension(3) :: gravflux
      
      
      if (present(value)) then

        quadpnt%type_pnt = "ndpt"
        quadpnt%column = 2
        quadpnt%order = elements%data(el_id, node_order)
        layer = elements%material(el_id)
        D = drutes_config%dimen
        select case (pde_loc%mfswitch)
          case("m")
            call pde(1)%pde_fnc(1)%dispersion(pde_loc, layer, quadpnt, tensor=K(1:D,1:D))          
          case("f")
            call pde(2)%pde_fnc(2)%dispersion(pde_loc, layer, quadpnt, tensor=K(1:D,1:D))         
        end select
	
	
      	select case(D)
          case(1)
	  
            value = K(1,1) * elements%nvect_z(el_id, node_order)
          
          case(2)	  
            gravflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*K(1,2)
            
            gravflux(2) = elements%nvect_z(el_id, node_order)*K(2,2)

            value = sqrt(gravflux(1)*gravflux(1) + gravflux(2)*gravflux(2))

        end select
      end if
      
      if (present(code)) then
        code = 2
      end if
      
    end subroutine dual_freedrainage
     
  !> Assigns a simple atmospheric boundary based on water content in the soil and potential evaporation and rain for dual permeability model
  !!  \f$rain - evap*\theta(h)^{\frac{2}{3}}\f$.
  subroutine dual_atmospheric(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals

      !> takes the global pde_loc structure
      class(pde_str), intent(in) :: pde_loc 
      !> takes the element id and node order
      integer(kind=ikind), intent(in)  :: el_id, node_order
      !> optional value output of the assigned boundary 
      real(kind=rkind), intent(out), optional    :: value
      !> optional code output. Different codes for different boundaries. Code = 2 
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      
      
      
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind) :: theta, rain, evap
      integer(kind=ikind) :: i, edge_id, j
      
      

      if (present(code)) then
	    code = 2
      end if
      
      if (present(value)) then
        edge_id = nodes%edge(elements%data(el_id, node_order))

        i = pde_loc%permut(elements%data(el_id, node_order))


        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
          j = i-1
              else
          j = i
              end if
              rain = pde_loc%bc(edge_id)%series(j,2)
              evap = pde_loc%bc(edge_id)%series(j,3)
              EXIT
            end if
          end do
        else
          print *, "atmospheric boundary must be time dependent, check record for the boundary", edge_id
          ERROR STOP
        end if

        quadpnt%type_pnt = "ndpt"
        quadpnt%column = 2
        quadpnt%order = elements%data(el_id,node_order)
        layer = elements%material(el_id)
        theta =  pde_loc%mass(1)%val(pde_loc, layer, quadpnt)
        select case (pde_loc%mfswitch)
          case("m")
            value = (rain - evap*theta**(2.0_rkind/3.0_rkind))*(1_rkind-infweight)  
          case("f")
            value = (rain - evap*theta**(2.0_rkind/3.0_rkind))*(infweight)
        end select


      end if
      
    end subroutine dual_atmospheric
     
end module Re_dual_bc
