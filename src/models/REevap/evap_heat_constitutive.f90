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

  public :: heatcap_TT, heatcap_Th, heat_cond, convection4heat, heatdiffTh,  heatsrc_w_roots, latentheat, heat_flux4evap

  private :: water_cap, vapour_cap
  
  contains
  
  
  
    subroutine heat_cond(pde_loc, layer, quadpntin, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use re_constitutive
      use evapglob
      use evap_RE_constitutive
      use debug_tools
      
            
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      
      real(kind=rkind) :: theta, val
      integer(kind=ikind) :: i, D
      
      real(kind=rkind), dimension(3,3) :: Klvt
      type(integpnt_str) :: quadpnt
      
      D = drutes_config%dimen
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_heat_constitutive::heat_cond"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.

      
      theta = vangen(pde(re_ord), layer, quadpnt)
      
      call cond_vt(layer, quadpnt, Klvt)
      
      val = soil_heat_coef(layer)%b1 + soil_heat_coef(layer)%b2*theta + soil_heat_coef(layer)%b3*sqrt(theta)
      
      
      if (present(tensor)) then
        tensor = 0
        do i=1, drutes_config%dimen
          tensor(i,i) = val
        end do
        tensor = tensor + latentheat(quadpnt)*dens_liquid(quadpnt)*Klvt(1:D, 1:D)
      end if
      
      if (present(scalar)) scalar = tensor(1,1)
          
      
    end subroutine heat_cond
  
    function heatcap_TT(pde_loc, layer, quadpntin, x) result(val)
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
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta, theta_v
      type(integpnt_str) :: quadpnt
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_heat_constitutive::heatcap_TT"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.

      ths = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      theta_v = thetav(pde(re_ord), layer, quadpnt) 
      
      val = heatpar(layer)%C*(1-ths) + water_cap(quadpnt)*theta + vapour_cap(quadpnt)*(ths-theta) + &
            dens_liquid(quadpnt)*latentheat(quadpnt)*dthetav_dtemp(pde(re_ord), layer, quadpnt)
      
    end function heatcap_TT
    
    
    function heatcap_Th(pde_loc, layer, quadpntin, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use evap_RE_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta, theta_v
      type(integpnt_str) :: quadpnt
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_heat_constitutive::heatcap_Th"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      val = dens_liquid(quadpnt)*latentheat(quadpnt)*dthetav_dh(pde(re_ord), layer, quadpnt)
      
    end function heatcap_Th
     
    subroutine heat_flux4evap(pde_loc, layer, quadpntin, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use evapglob
      use evap_RE_constitutive
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpntin    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
      real(kind=rkind), dimension(:), allocatable, save :: gradT, gradh
      
      real(kind=rkind), dimension(3,3) :: cond, Kvh, Klvt
      real(kind=rkind), dimension(3) :: qconvect, flux_loc
      integer(kind=ikind) :: D
      type(integpnt_str) :: quadpnt
      
      D = drutes_config%dimen
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_heat_constitutive::heat_flux4evap"
        ERROR STOP
      end if
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      call pde(heat_ord)%getgrad(quadpnt, gradT)
      
      call pde(re_ord)%getgrad(quadpnt, gradh)

      call heat_cond(pde(heat_ord), layer, quadpnt, tensor = cond(1:D, 1:D))
      
      call convection4heat(pde(heat_ord), layer, quadpnt, vector_out = qconvect(1:D))
      
      call cond_vapour4h(layer, quadpnt, Kvh(1:D, 1:D))
      
      call cond_vt(layer, quadpnt, Klvt(1:D, 1:D))
      
      flux_loc(1:D) = -matmul(cond(1:D, 1:D), gradT) +  qconvect(1:D)*pde(heat_ord)%getval(quadpnt)  &
                - latentheat(quadpnt)*dens_liquid(quadpnt)*matmul(Kvh(1:D, 1:D),gradh) &
                - latentheat(quadpnt)*dens_liquid(quadpnt)*matmul(Klvt(1:D, 1:D),gradT)
      
      if (present(flux)) then
        flux = flux_loc(1:D)
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(flux_loc(1:D))
      end if
    
    end subroutine heat_flux4evap
    
    subroutine convection4heat(pde_loc, layer, quadpntin, x, vector_in, vector_out, scalar)
      use typy
      use globals
      use global_objs
      use pde_objs
      use evapglob
      use evap_RE_constitutive

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpntin    
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
      type(integpnt_str) :: quadpnt
      
      D = drutes_config%dimen
      
      if (.not. allocated(ql)) then
        allocate(ql(D))
        allocate(qv(D))
      end if
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_heat_constitutive::convection4heat"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
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
    
    
    
    
    
    subroutine heatdiffTh(pde_loc, layer, quadpntin, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use evap_RE_constitutive
            
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      
      integer(kind=ikind) :: D
      real(kind=rkind), dimension(3,3) :: Kvh
      type(integpnt_str) :: quadpnt
      
      if (.not. present(quadpntin)) then
        print *, "runtime error evap_heat_constitutive::heatdiffTh"
        ERROR STOP
      end if
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      D = drutes_config%dimen
      
      call cond_vapour4h(layer, quadpnt, Kvh(1:D, 1:D))
      
      if (present(tensor)) then
        tensor(1:D, 1:D) = dens_liquid(quadpnt)*latentheat(quadpnt)*Kvh(1:D, 1:D)
      end if
      
      if (present(scalar)) then
        scalar = dens_liquid(quadpnt)*latentheat(quadpnt)*Kvh(1,1)
      end if
      
    end subroutine heatdiffTh
    
    
    function heatsrc_w_roots(pde_loc, layer, quadpntin, x) result(val)
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
      type(integpnt_str), intent(in), optional :: quadpntin
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
      val = -C_liq*dens_liquid(quadpnt)*sinkterm(pde(re_ord), layer, quadpnt) +  heatpar(layer)%source
      
    end function heatsrc_w_roots
    

    
    function water_cap(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use evap_RE_constitutive
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = C_liq * dens_liquid(quadpnt)
    
    end function water_cap
    
    
    function vapour_cap(quadpntin) result(val)
      use typy
      use global_objs
      use evapglob
      use evap_RE_constitutive
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpntin
      real(kind=rkind) :: val
      
      type(integpnt_str) :: quadpnt
      
      quadpnt = quadpntin
      quadpnt%preproc = .true.
      
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
