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


module evap_heat_constitutive

  public :: heatcap_TT, heatcap_hT, heat_cond, convection4heat, heatdiffTh,  heatsrc_w_roots, latentheat

  private :: water_cap, vapour_cap
  
  contains
  
  
  
    subroutine heat_cond(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use re_constitutive
      use evapglob
      
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
      
      real(kind=rkind) :: theta, val
      integer(kind=ikind) :: i
      
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::heat_cond"
        ERROR STOP
      end if
      
      theta = vangen(pde(re_ord), layer, quadpnt)
      
      val = soil_heat_coef(layer)%b1 + soil_heat_coef(layer)%b2*theta + soil_heat_coef(layer)%b3*sqrt(theta)
      
      if (present(tensor)) then
        tensor = 0
        do i=1, drutes_config%dimen
          tensor(i,i) = val
        end do
      end if
      
      if (present(scalar)) scalar = val
          
      
    end subroutine heat_cond
  
    function heatcap_TT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use heat_globals
      use re_constitutive
      use evap_RE_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta, theta_v
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::heatcap_TT"
        ERROR STOP
      end if
      
      ths = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      theta_v = thetav(pde(re_ord), layer, quadpnt) 
      
      val = heatpar(layer)%C*(1-ths) + water_cap(quadpnt)*theta + vapour_cap(quadpnt)*theta_v + &
            latentheat(quadpnt)*dthetav_dtemp(pde(re_ord), layer, quadpnt)
      
    end function heatcap_TT
    
    
    subroutine convection4heat(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use global_objs
      use pde_objs
      use evapglob
      use evap_RE_constitutive

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
      
      real(kind=rkind), dimension(:), allocatable, save :: ql, qv
      integer(kind=ikind) :: D
      
      D = drutes_config%dimen
      
      if (.not. allocated(ql)) then
        allocate(ql(D))
        allocate(qv(D))
      end if
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::convection4heat"
        ERROR STOP
      end if
      
      call darcy4liq(pde(re_ord), layer, quadpnt, flux=ql)
      call darcy4vap(pde(re_ord), layer, quadpnt, flux=qv)
      
!      ql =0
!      qv = 0
      
      if (present(vector_out)) then
        vector_out(1:D) = C_liq*dens_liquid(quadpnt)*ql + C_vap*relhumid(quadpnt)*dens_satvap(quadpnt)*qv
      end if
      
      if (present(scalar)) then
        scalar = norm2(C_liq*dens_liquid(quadpnt)*ql + C_vap*relhumid(quadpnt)*dens_satvap(quadpnt)*qv)
      end if
      
    end subroutine convection4heat
    
    
    subroutine heatdiffTh(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use evap_RE_constitutive
            
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
      
      integer(kind=ikind) :: D
      real(kind=rkind), dimension(3,3) :: Kvh
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::heatdiffTh"
        ERROR STOP
      end if
      
      D = drutes_config%dimen
      
      call cond_vapour4h(layer, quadpnt, Kvh(1:D, 1:D))
      
      if (present(tensor)) then
        tensor(1:D, 1:D) = latentheat(quadpnt)*Kvh(1:D, 1:D)
      end if
      
      if (present(scalar)) then
        scalar = latentheat(quadpnt)*Kvh(1,1)
      end if
      
    end subroutine heatdiffTh
    
    
    function heatsrc_w_roots(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
      use evapglob
      use re_constitutive
      use evap_RE_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      
      val = -C_liq*dens_liquid(quadpnt)*sinkterm(pde(re_ord), layer, quadpnt) +  heatpar(layer)%source
      
      
    end function heatsrc_w_roots
    
    
    function heatcap_hT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use evap_RE_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::heatcap_hT"
        ERROR STOP
      end if
    
      val = latentheat(quadpnt)*dthetav_dh(pde(re_ord), layer, quadpnt)
      
    end function heatcap_hT
    
    function water_cap(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use evap_RE_constitutive
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = pde(heat_ord)%getval(quadpnt) * dens_liquid(quadpnt)
    
    end function water_cap
    
    
    function vapour_cap(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use evap_RE_constitutive
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = C_vap * relhumid(quadpnt) * dens_satvap(quadpnt)
    
    end function vapour_cap
    
    
    function latentheat(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = 2.501e6 - 2369.2*pde(heat_ord)%getval(quadpnt)
    
    end function latentheat


end module evap_heat_constitutive
