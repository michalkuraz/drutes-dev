
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

!> \file heat_pointers.f90
!! \brief Pointer linker for heat equation. 
!! based on (Sopohocleous, 1979)
!! \f[ C\frac{\partial T}{\partial t} = \mathbf{\lambda} \Delta T - C_w \nabla \cdot \vec{q} T  \f]
!<

module heat_pointers
  public :: heat, heat_processes
  public :: heatlinker

  
  contains
  
    subroutine heat_processes(processes)
      use typy
      use readtools
      use heat_globals
      use core_tools
      use globals
      
      integer(kind=ikind) processes
      integer :: heatconf, ierr
      
      open(newunit=heatconf, file="drutes.conf/heat/heat.conf", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
        call write_log("Unable to open drutes.conf/heat/heat.conf, exiting....")
        ERROR STOP
      end if
      
      call fileread(with_richards, heatconf)
      
      close(heatconf)
      
      if (with_richards .or. cut(drutes_config%name) == "REevap" ) then
        processes = 2
      else
        processes = 1
      end if
      
    end subroutine heat_processes
  
    subroutine heat(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_fnc
      use heat_reader
      use re_pointers
      
      class(pde_str), dimension(:), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, pos
      
      if (ubound(pde,1) == 1) then
        call heatlinker(pde(1))
	    else
!        if (drutes_config%name == "REevap") then
!          call RE_std(pde(1))
!        else
          call REstdH(pde(1))
!        end if
	     call heatlinker(pde(2))
      end if
            
    
    end subroutine heat
    
    
    subroutine heatlinker(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_fnc
      use heat_reader
      use re_constitutive
!       use evap_bc
      
      
      class(pde_str),  intent(in out) :: pde_loc
      integer(kind=ikind) :: i
      
      call heat_read(pde_loc)
      
      pde_common%nonlinear = .false. 
      
      pde_loc%pde_fnc(pde_loc%order)%dispersion => heat_conduct
      
      pde_loc%pde_fnc(pde_loc%order)%convection => heat_convect

      pde_loc%pde_fnc(pde_loc%order)%elasticity => heat_elast
            
      pde_loc%pde_fnc(pde_loc%order)%zerord => heat_source
      
	  
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
        select case(pde_loc%bc(i)%code)
          case(1)
            pde_loc%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde_loc%bc(i)%value_fnc => heat_neumann
          case(0)
            pde_loc%bc(i)%value_fnc => re_null_bc
          case(3)
            CONTINUE
            !to be linked in vapour_pointers
          case default
            print *, "unrecognized bc option"
            print *, "exited from heat_pointers::heatlinker"
            ERROR STOP
        end select
      end do    
	
      pde_loc%flux => heat_flux
      
      pde_loc%initcond => heat_icond  
      
      pde_loc%print_mass = .false.
      
      deallocate(pde_loc%mass)
      
      allocate(pde_loc%mass(0))
      

    
    end subroutine heatlinker
    
    function heat_source(pde_loc, layer, quadpnt, x) result(val)
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
      !> return value
      real(kind=rkind)                :: val 
      
      integer(kind=ikind) :: n, i
      real(kind=rkind) :: theta
      
      if (pde_loc%order == 2) then
        theta = pde(1)%mass(1)%val(pde_loc, layer, quadpnt)
      else
        theta = adepar(layer)%water_cont
      end if
      
       val = 0.0_rkind
       do i=1, ubound(adepar(layer)%orders,1)
        if (abs(adepar(layer)%orders(i)) < 100*epsilon(1.0_rkind)) then
          val =  val + theta*adepar(layer)%lambda(i)
        end if
      end do
      
      
    end function heat_source
    
    
      
  

end module heat_pointers
