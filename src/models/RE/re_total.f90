
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

!> \file re_total.f90
!! \brief Total hydraulic head for the Richards equation.
!<

!> The Richards equation can be solved either in pressure head form $\f h\f$ in in total hydraulic head \f$ H \f$. Total hydraulic head requires relinking of pde%getval function. However, total hydraulic head is more convenient. Symmetric operator, can handle non-uniform material properties (different soil layers).
!<

module re_total
  public :: getval_retot
  public :: darcy4totH
  public :: retot_dirichlet_bc, retot_dirichlet_height_bc, retot_neumann_bc, retot_freedrainage
  public :: retot_atmospheric
  public :: retot_initcond
  public :: iconddebouss
  public :: retot_seepage
  public :: retot_well
  
  contains
  
    !> specific function for Richards equation in H-form (total hydraulic head form), replaces pde_objs::getvalp1 in order to distinguish between H and h 
    function getval_retot(pde_loc, quadpnt) result(val)
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
      

           
      if (quadpnt%preproc) then
      
        D = drutes_config%dimen
             
        call getcoor(quadpnt, xyz(1:D))

        if (drutes_config%dimen>1) then
          val = getvalp1(pde_loc, quadpnt) - xyz(D)
        else
          layer = get_layer(quadpnt)
          val = getvalp1(pde_loc, quadpnt) - xyz(D)*cos(4*atan(1.0_rkind)/180*vgset(layer)%anisoangle(1))
        end if
        
        
      else
        val = getvalp1(pde_loc, quadpnt)
      end if
	
      
    end function getval_retot
    
    
    subroutine darcy4totH(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use re_globals
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length

      real(kind=rkind), dimension(3,3)  :: K
      integer(kind=ikind)               :: D
      integer(kind=ikind)               :: i
      real(kind=rkind), dimension(:), allocatable, save  :: gradH
      real(kind=rkind), dimension(:), allocatable, save  :: vct
      real(kind=rkind) :: h
      type(integpnt_str) :: quadpnt_loc
      
      D = drutes_config%dimen

      if (.not.(allocated(gradH))) then
        allocate(gradH(1:D))
        allocate(vct(1:D))
      end if

      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
        print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
        print *, "exited from re_constitutive::darcy_law"
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
      
      
      
      call pde_loc%pde_fnc(pde_loc%order)%dispersion(pde_loc, layer, x=(/h/), tensor=K(1:D, 1:D))

      
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

    end subroutine darcy4totH
    
    subroutine retot_dirichlet_bc(pde_loc, el_id, node_order, value, code, array, bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts
      
      integer(kind=ikind) :: edge_id, i, j, D
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: xyz

      
      
      edge_id = nodes%edge(elements%data(el_id, node_order))

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      
      if (present(value)) then
        if (pde_loc%bc(edge_id)%file) then
          do i=1, ubound(pde_loc%bc(edge_id)%series,1)
            if (pde_loc%bc(edge_id)%series(i,1) > time) then
              if (i > 1) then
                j = i-1
              else
                j = i
              end if
              value = pde_loc%bc(edge_id)%series(j,2) + xyz(D)
              EXIT
            end if
          end do
        else
          value = pde_loc%bc(edge_id)%value + xyz(D)
        end if
      end if
      

      if (present(code)) then
        code = 1
      end if
      
      
    end subroutine retot_dirichlet_bc
    



    subroutine retot_dirichlet_height_bc(pde_loc, el_id, node_order, value, code, array,bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      use re_globals
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts

      
      integer(kind=ikind) :: edge_id, i, j
      real(kind=rkind) :: tempval, node_height
      
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      node_height = nodes%data(elements%data(el_id, node_order), drutes_config%dimen)

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
      
    end subroutine retot_dirichlet_height_bc
    
    subroutine retot_well(pde_loc, el_id, node_order, value, code, array,bcpts)
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional   :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts

      integer(kind=ikind) :: edge_id, i, j, D
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: xyz

      
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      

      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id, node_order)
      D = drutes_config%dimen
      call getcoor(quadpnt, xyz(1:D))
      
      
      if ( xyz(D) <= pde_loc%bc(edge_id)%value ) then
      
        if (present(value)) then
          value = pde_loc%bc(edge_id)%value 
        end if
        

        if (present(code)) then
          code = 1
        end if
      
      else
        if (present(value)) then
          value = 0
        end if
        

        if (present(code)) then
          code = 2
        end if
      end if
      
      
    end subroutine retot_well
    
    subroutine retot_seepage(pde_loc, el_id, node_order, value, code, array,bcpts)
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      use re_globals
      use geom_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts
      
      
      real(kind=rkind) :: solval, gradn
      real(kind=rkind), dimension(:), allocatable, save :: solgrad, nvect
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: i, el_vecino, el_vecino2, nd_tmp, nd, counter, status, nd_vecino, nd_vecino2, pos, test1, test2, nd3
      integer(kind=ikind) :: myel, edge_id
      real(kind=rkind), dimension(2,2) :: points
      real(kind=rkind), dimension(2) :: third
 

              
      
      if (.not. allocated(pde_common%xvect) ) then
        if (present(value)) value = 0
        if (present(code)) code = 2
        RETURN
      end if
     
      
      nd = elements%data(el_id, node_order)
      
      quadpnt%type_pnt = "ndpt"
      quadpnt%order = elements%data(el_id,node_order)
      quadpnt%column = 1
        
      quadpnt%preproc = .true.
        
      solval = pde_loc%getval(quadpnt)
        
      quadpnt%preproc = .false.
        
      call pde_loc%getgrad(quadpnt, solgrad)
      
 
      if (drutes_config%dimen > 1) then
      
        status = nodes%element(nd)%pos
        
        if (status == 3 .and. elements%border(el_id)%pos > 0) then
          status = 2
        end if
        
        select case(status)
          case(1)
            el_vecino = el_id
            el_vecino2 = el_id
          case(2)      
            do i=1, nodes%element(nd)%pos
              el_vecino = nodes%element(nd)%data(i)
              if (elements%border(el_vecino)%pos > 0 .and. el_vecino /= el_id) then
                EXIT
              end if
              if (i==nodes%element(nd)%pos) el_vecino=0
            end do
            el_vecino2 = el_id
          case default
            do i=1, nodes%element(nd)%pos
              el_vecino = nodes%element(nd)%data(i)
              if (elements%border(el_vecino)%pos > 0 .and. el_vecino /= el_id) then
                EXIT
              end if
              if (i==nodes%element(nd)%pos) el_vecino=0
            end do
            
            do i=1, nodes%element(nd)%pos
              el_vecino2 = nodes%element(nd)%data(i)
              if (elements%border(el_vecino2)%pos > 0 .and. el_vecino2 /= el_id .and. el_vecino2 /= el_vecino) then
                EXIT
              end if
              if (i==nodes%element(nd)%pos) el_vecino2=0
            end do
        end select
        
        
        if (el_vecino == 0 .or. el_vecino2 == 0) then
          print *, "bug in re_total::re_seepage "
          print *,  "bug keyword: all neighbours are zeroes (don't worry if you don't understand it)"
          print *, "contact Michal -> michalkuraz@gmail.com"
          ERROR STOP
        end if
        
        

        do i=1, elements%border(el_vecino)%pos
          if (elements%border(el_vecino)%data(i) == nd) then
            if (elements%border(el_vecino)%pos < 3 .or. nodes%element(nd)%pos == 1) then
              if (i < elements%border(el_vecino)%pos) then
                nd_vecino = elements%border(el_vecino)%data(i+1)
              else
                nd_vecino = elements%border(el_vecino)%data(1)
              end if
              EXIT
            else
              if (i < elements%border(el_vecino)%pos .and. i > 1) then
                test1 = i+1
                test2 = i-1
              else if (i == elements%border(el_vecino)%pos ) then
                test1 = 1
                test2 = i-1
              else if (i == 1) then
                test1 = 2
                test2 = elements%border(el_vecino)%pos
              end if
              
              if (nodes%element(elements%data(el_vecino,test1))%pos == 1) then
                nd_vecino = elements%data(el_vecino,test1)
              else if  (nodes%element(elements%data(el_vecino,test2))%pos == 1) then
                nd_vecino = elements%data(el_vecino,test2)
              else
                print *, "bug in re_total::re_seepage, bug keyword: unable to find neighbour vecino1 "
                print *, "(don't worry if you don't understand it)"
                print *, "contact Michal -> michalkuraz@gmail.com"
                ERROR STOP  
              end if   
              EXIT        
            end if
          end if
        end do
                
        do i=1, elements%border(el_vecino2)%pos
          if (elements%border(el_vecino2)%data(i) == nd) then
            if (elements%border(el_vecino2)%pos < 3 .or. nodes%element(nd)%pos == 1) then
              if (i < elements%border(el_vecino2)%pos) then
                nd_vecino2 = elements%border(el_vecino2)%data(i+1)
                if (nd_vecino2 == nd_vecino) then
                  pos = i-1
                  if (pos < 1) then
                    pos = elements%border(el_vecino)%pos
                  end if
                  nd_vecino2 = elements%border(el_vecino2)%data(pos)
                end if
              else
                nd_vecino2 = elements%border(el_vecino2)%data(1)
                if (nd_vecino2 == nd_vecino) then
                  pos = i-1
                  nd_vecino2 = elements%border(el_vecino2)%data(pos)
                end if
              end if
              EXIT
            else 
              if (i < elements%border(el_vecino2)%pos .and. i > 1) then
                test1 = i+1
                test2 = i-1
              else if (i == elements%border(el_vecino2)%pos ) then
                test1 = 1
                test2 = i-1
              else if (i == 1) then
                test1 = 2
                test2 = elements%border(el_vecino2)%pos
              end if
              
              if (nodes%element(elements%data(el_vecino2,test1))%pos == 1) then
                nd_vecino2 = elements%data(el_vecino2,test1)
              else if  (nodes%element(elements%data(el_vecino2,test2))%pos == 1) then
                nd_vecino2 = elements%data(el_vecino2,test2)
              else
                print *, "bug in re_total::re_seepage, bug keyword: unable to find neighbour vecino2 "
                print *, "(don't worry if you don't understand it)"
                print *, "contact Michal -> michalkuraz@gmail.com"
                ERROR STOP  
              end if   
              EXIT         
            end if
          end if
        end do
        
        points(1,:) = nodes%data(nd_vecino,:)
        points(2,:) = nodes%data(nd_vecino2,:)
        
        
        if (elements%border(el_id)%pos < 2) then
          myel = el_id
        else
          myel = elements%neighbours(el_id,1)
        end if
        
        do i=1, ubound(elements%data,2)
          nd3 =  elements%data(myel, i)
          if ( nd3 /= nd .and. nd3 /= nd_vecino .and. nd3 /= nd_vecino2) then
            third = nodes%data(nd3,:)
            EXIT
          end if
          
          if (i == ubound(elements%data,2)) then
            print *,  "bug in re_total::re_seepage, bug keyword: unable to find internal element node"
            print *, "(don't worry if you don't understand it)"
            print *, "contact Michal -> michalkuraz@gmail.com"
            ERROR STOP  
          end if
        end do  
      
      
        call getnormal(points, third, nvect)    
      

        gradn = solgrad(1)*nvect(1) + solgrad(2)*nvect(2)
      
      else
     
        if (node_order == 1) then
          gradn = -solgrad(1)
        else
          gradn = solgrad(1)
        end if
      
      end if
     
      edge_id = nodes%edge(elements%data(el_id, node_order))
      
      if (solval < pde_loc%bc(edge_id)%value .or. gradn > 0) then
        if (present(code)) code = 2
        if (present(value)) value = 0
      else 
        if (present(code)) code = 4
        if (present(value)) value = nodes%data(nd,drutes_config%dimen)
      end if 

      
    end subroutine retot_seepage
    




    subroutine retot_neumann_bc(pde_loc, el_id, node_order, value, code, array,bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use debug_tools

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts
     

      integer(kind=ikind) :: i, edge_id, j
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval
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
              bcval = pde_loc%bc(edge_id)%series(j,2)
              EXIT
            end if
          end do
        else
          bcval = pde_loc%bc(edge_id)%value
        end if
    
        value = bcval

      end if
      
      if (present(code)) then
        code = 2
      end if


    end subroutine retot_neumann_bc
    
    subroutine retot_atmospheric(pde_loc, el_id, node_order, value, code, array,bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts
      
      
      
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind) :: theta, rain, evap, h
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
        quadpnt%order = elements%data(el_id,node_order)
        quadpnt%column = 3
        layer = elements%material(el_id)
        theta =  pde_loc%mass(1)%val(pde_loc,layer, quadpnt)
        quadpnt%preproc = .true.
        quadpnt%column = 1
        h = pde_loc%getval(quadpnt)
        if ( h < h_crit_low .OR. theta < theta_crit ) evap = 0
        if ( h >= h_crit_high) rain = 0
        
        value = rain + evap

      end if
      
      
    end subroutine retot_atmospheric
    
    
    !> if h > 0 it becomes Dirichlet condition h=0, similar to seepage face
    !! solution for VUMOP runoff simulator
    !! the state when rain finishes and it starts to evaporate again is not implemented here
    !! use with a caution
    !<
    subroutine retot_atmosphericv2(pde_loc, el_id, node_order, value, code, array,bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts
      logical, save :: turned_on = .true.
      
      
      
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind) :: theta, rain, evap, h
      integer(kind=ikind) :: i, edge_id, j
      
      
      if (present(code)) code = 2



      if (.not. turned_on) then
        if (present(value)) value = 0
        if (present(code)) code = 4
        
        RETURN
      end if

!      print *, present(value)

      if (present(value) .and. turned_on) then
        quadpnt%preproc = .true.
        quadpnt%column = 1
        quadpnt%type_pnt = "ndpt"
        quadpnt%order = elements%data(el_id,node_order)
        h = pde_loc%getval(quadpnt)

        if ( h >= h_crit_high) then
          value = 0
          if (present(code)) code = 4
          turned_on = .false.
          RETURN
        end if

        if ( h < h_crit_low ) evap = 0
        
        
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
        quadpnt%order = elements%data(el_id,node_order)
        quadpnt%column = 3
        layer = elements%material(el_id)
!        theta =  pde_loc%mass(1)%val(pde_loc,layer, quadpnt)
        quadpnt%preproc = .false.
        quadpnt%column = 1
 
                
        value = rain + evap

      end if
      
      
    end subroutine retot_atmosphericv2
    
    subroutine retot_freedrainage(pde_loc, el_id, node_order, value, code, array,bcpts) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(:), intent(out), optional :: array
      type(bcpts_str), intent(in), optional :: bcpts
      
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
        call pde_loc%pde_fnc(pde_loc%order)%dispersion(pde_loc, layer, quadpnt, tensor=K(1:D,1:D))
        
        select case (drutes_config%dimen)
          case(1)
            gravflux(1) = K(1,1) * elements%nvect_z(el_id, node_order)*(-1)
          case(2)
            gravflux(1) = K(2,1) * elements%nvect_x(el_id, node_order)*(-1)
            gravflux(2) = K(2,2) * elements%nvect_z(el_id, node_order) * (-1)
          case(3)
            gravflux(1) = K(3,1) * elements%nvect_x(el_id, node_order)*(-1)
            gravflux(2) = K(3,2) * elements%nvect_y(el_id, node_order) * (-1)
            gravflux(3) = K(3,3) * elements%nvect_z(el_id, node_order)*(-1)

        end select
    
        value = norm2(gravflux(1:D))
      end if
      
      if (present(code)) then
        code = 2
      end if
      
    end subroutine retot_freedrainage
    
    


    subroutine retot_initcond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use geom_tools
      use core_tools
      use read_inputs

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value


      D = drutes_config%dimen
      if (cut(vgset(1_ikind)%icondtype) == "input") then
        call map1d2dJ(pde_loc,"drutes.conf/water.conf/hini.in", correct_h = .false.)
        RETURN  
      end if
      
      if (cut(vgset(1_ikind)%icondtype) == "ifile") then
        call read_icond(pde_loc, "drutes.conf/water.conf/RE_init.in")
        do i=1, nodes%kolik
          pde_loc%solution(i) = pde_loc%solution(i) + nodes%data(i,D)
        end do
        RETURN  
      end if
        
      
      D = drutes_config%dimen
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
            select case (vgset(layer)%icondtype)
              case("H_tot")
                pde_loc%solution(k) = vgset(layer)%initcond !+ nodes%data(k,1)
              case("hpres")
                pde_loc%solution(k) = vgset(layer)%initcond + nodes%data(k,D)
              case("theta")
                value = inverse_vangen(pde_loc, layer, x=(/vgset(layer)%initcond/))
                pde_loc%solution(k) = value + nodes%data(k,D)
            end select
          end if
        end do   
      end do
      
      
  

    end subroutine retot_initcond
    
    
    subroutine iconddebouss(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use core_tools
      use debug_tools
      
      class(pde_str), intent(in out) :: pde_loc
      
      real(kind=rkind), dimension(:,:), allocatable :: boussdata, slopesdata
      real(kind=rkind), dimension(:), allocatable :: ndheight
      
      integer(kind=ikind) :: i, counter, j
      
      integer, dimension(2) :: fileid
      real(kind=rkind) :: tmp, dx, slope
      integer :: ierr
      
           
      call find_unit(fileid(1), 200)
      open(unit=fileid(1), file="drutes.conf/boussinesq.conf/iconds/slopes.csv", action="read", status="old")
      
      call find_unit(fileid(2), 200)
      open(unit=fileid(2), file="drutes.conf/boussinesq.conf/iconds/bouss_water_table-1.dat", action="read", status="old")
      
      counter=0
      do 
        read(unit=fileid(1), fmt=*, iostat=ierr) tmp
        if (ierr == 0) then
          counter = counter+1
        else
          exit
        end if
      end do
      
      allocate(slopesdata(counter,3))
      
      close(fileid(1))
      open(unit=fileid(1), file="drutes.conf/boussinesq.conf/iconds/slopes.csv", action="read", status="old")
      
      
      do i=1, ubound(slopesdata,1)
        read(fileid(1), fmt=*) slopesdata(i,:)
      end do
      
      counter=0
      do 
        read(unit=fileid(2), fmt=*, iostat=ierr) tmp
        if (ierr == 0) then
          counter = counter+1
        else
          exit
        end if
      end do
      
      allocate(boussdata(counter,2))
      
      close(fileid(2))
      open(unit=fileid(2), file="drutes.conf/boussinesq.conf/iconds/bouss_water_table-1.dat", action="read", status="old")
      
      do i=1, ubound(boussdata,1)
        read(fileid(2), fmt=*) boussdata(i,:)
      end do
      boussdata(:,1) = 560 - boussdata(:,1)
      

      allocate(ndheight(nodes%kolik))
      do i=1, nodes%kolik
        do j=1, ubound(slopesdata,1) -1
          if (nodes%data(i,1) >=slopesdata(j,1) .and. nodes%data(i,1) < slopesdata(j+1,1))  then
            slope = slopesdata(j,2)
            dx = nodes%data(i,1) - slopesdata(j,1)
            ndheight(i) = dx*slope + slopesdata(j,3)
            EXIT
          end if
        end do
      end do
      

      
      do i=1, nodes%kolik
        do j=1, ubound(boussdata,1) - 1
          if (nodes%data(i,1) >=boussdata(j+1,1) .and. nodes%data(i,1) < boussdata(j,1) ) then!.or. &
          !j==ubound(boussdata,1) -1) then
            slope = (boussdata(j+1,2) - boussdata(j,2))/  (boussdata(j+1,1) - boussdata(j,1))
            pde_loc%solution(i) = boussdata(j+1,2) + (nodes%data(i,1)-boussdata(j+1,1))*slope + ndheight(i)
            EXIT
          end if
        end do
        if (nodes%data(i,1) > boussdata(1,1)) then
          print *, maxval(nodes%data(:,1)), maxval(boussdata(:,1))
          STOP
        end if
      end do
	  
    
    end subroutine iconddebouss

end module re_total
