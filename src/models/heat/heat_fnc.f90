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

!> \file heat_fnc.f90
!! \brief Heat transport equation with convection (Sophocleous, 1979)
!<

!> Heat equation is solved in the following form
!! \f[  C_T \frac{\partial T}{\partial t} = \overbrace{\kappa_T \Delta T}^{\mbox{ diffusion term}}  - \underbrace{\nabla \cdot (C_T^L \vec{q}_L)}_{\mbox{ convection term}} + \overbrace{S_H}^{\mbox{ reaction term}}. \f]
!! Convection term can be user defined or obtained from coupling with Richards equation.
!<


module heat_fnc
  public :: heat_conduct
  public :: heat_convect
  public :: heat_source
  public :: heat_dirichlet, heat_neumann
  public :: heat_flux
  public :: heat_icond
  public :: heat_elast
  
  contains
    subroutine heat_conduct(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use heat_globals
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

      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen

      
      if (present(tensor)) then
        tensor =  heatpar(layer)%lambda
      end if
      
    
    end subroutine heat_conduct
    
    
    subroutine heat_convect(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
      use re_total
      
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
 
      
      if (with_richards) then
        if (present(vector_out)) then
          call pde(1)%flux(layer, quadpnt, vector_out=vector_out)
          vector_out=heatpar(layer)%C_w*vector_out
        end if
        
        if (present(scalar)) then
           call pde(1)%flux(layer, quadpnt, scalar=scalar)
           scalar=heatpar(layer)%C_w*scalar
        end if
        
      else
        
        if (present(vector_out)) then
          vector_out = heatpar(layer)%C_w*heatpar(layer)%convection
        end if
      
        
        if (present(scalar)) then
          scalar = norm2(heatpar(layer)%C_w*heatpar(layer)%convection)
        end if
      end if
     
     
      
    end subroutine heat_convect
    
    function heat_elast(pde_loc,layer, quadpnt, x) result(E)
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


      E = heatpar(layer)%C
      

    end function heat_elast 

    function heat_source(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      
      val = heatpar(layer)%source
      
      
    end function heat_source
    
    subroutine heat_dirichlet(pde_loc, el_id, node_order, value, code, array, bcpts) 
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
      

    end subroutine heat_dirichlet
    
    
    subroutine heat_neumann(pde_loc, el_id, node_order, value, code, array,bcpts) 
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
        code = 2
      end if
      

    end subroutine heat_neumann
    
    
    
    
    subroutine heat_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
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
  
      
      if (present(flux)) then
        flux = -matmul(heatpar(layer)%lambda, gradT) 
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(matmul(heatpar(layer)%lambda, gradT))
      end if
    
    end subroutine heat_flux
    
    subroutine heat_icond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals
      use geom_tools

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
   
      do i=1, elements%kolik
        layer = elements%material(i)
        do j=1, ubound(elements%data,2)
          k = elements%data(i,j)
          l = nodes%edge(k)
          m = pde_loc%permut(k)
          if (m == 0) then
            call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
            pde_loc%solution(k) =  value 
          else
            if (.not. Tinitfile) then
              pde_loc%solution(k) = heatpar(layer)%Tinit
            else
              pde_loc%solution(k) = interpol_init1D(nodes%data(k,1), Tinitpoints)
            end if
              
          end if
        end do   
      end do

      
    
    end subroutine heat_icond
    
    

      
    subroutine heat_icondlin(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals
      
      class(pde_str), intent(in out) :: pde_loc
      real(kind=rkind) :: val1, val2
      integer(kind=ikind) :: i
      
      call pde_loc%bc(101)%value_fnc(pde_loc, 1_ikind, 1_ikind, val1)
      
      call pde_loc%bc(102)%value_fnc(pde_loc, elements%kolik, 2_ikind, val2)
      
      
      do i=1, ubound(pde_loc%solution,1)
        pde_loc%solution(i) = val1 + (val2-val1)/ubound(pde_loc%solution,1)*i
      end do
      
    
    end subroutine heat_icondlin
    
    

end module heat_fnc
