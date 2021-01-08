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


module evapbc4heat
  private :: RNterm, Hterm, Eterm, albedo_fnc, getmeteopos, resH
  public :: ebalance_flux, evaporation_bcflux
  
  contains
    function RNterm(quadpnt, layer) result(Rn)
      use typy
      use global_objs
      use evap_RE_constitutive
      use evapglob
      use pde_objs
      use re_constitutive
      
      
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind) :: Rn
      
      real(kind=rkind) :: Rns, Rnl, eps_s, Rld, Rlu, eps_a, ea, sigm = 5.67e-8, T_s, T_a, c, T_a_K, T_s_K 
      integer(kind=ikind) :: pos
      
      pos = getmeteopos()
      ! T_a for Boltzman law in kelvin
      T_a = meteo4evap(pos)%T_air
      T_a_K = T2kelv(T_a)
      c = meteo4evap(pos)%cloudiness
      
      Rns = (1-albedo_fnc(quadpnt, layer))*meteo4evap(pos)%inradiation

      eps_s = min(0.9 + 0.18*vangen(pde(re_ord),  layer, quadpnt), 1.0_rkind)
      
      ea = 0.611 * relhumid(quadpnt) * exp(17.27*T_a/(T_a_K - 35.85))
      
      eps_a = 1! 0.7 + 5.95e-5 * ea *exp (1500/T2kelv(T_a))
      
      Rld = ((1-0.84*c)*eps_a + 0.84*c) * sigm * T_a_K**4
      
      T_s = pde(heat_ord)%getval(quadpnt)

      T_s_K= T2kelv(T_s)
      
      Rlu = eps_s*sigm*T_s_K**4
      
      Rnl = eps_s*Rld - Rlu
      Rn = Rnl + Rns
    
    end function RNterm
    
    function getmeteopos() result(pos)
      use typy
      use globals
      use evapglob
      
      
      integer(kind=ikind) :: pos
      integer(kind=ikind) :: i
      
      if (time >= meteo4evap( ubound(meteo4evap,1))%time) then
        pos = ubound(meteo4evap,1)
        RETURN
      end if
      
      do i=meteo4evap(1)%datapos, ubound(meteo4evap,1) - 1
        if (time >= meteo4evap(i)%time .and. time <  meteo4evap(i+1)%time ) then
          meteo4evap(1)%datapos = i
          pos = i
          RETURN
        else if ( i == ubound(meteo4evap,1) - 1 .and. time >=   meteo4evap(i+1)%time) then
          meteo4evap(1)%datapos = i+1
          pos = i+1
        else if (i == 1 .and. time < meteo4evap(i)%time) then
          pos = 1
          RETURN
        end if
        
      end do
        

    
    end function getmeteopos
    
    function resH() result(val)
      use typy
      use evapglob
      
      real(kind=rkind) :: val
      real(kind=rkind) :: Ustar
      integer(kind=ikind) :: pos
      
      
      pos = getmeteopos()
      
      Ustar = meteo4evap(pos)%wind_speed * karman / (log((zref + 1e-3)/1e-3))
      
      val = 1/(karman * Ustar)*(log((zref+1e-3)/1e-3))
      
    end function resH
    
        
    function Hterm(quadpnt) result(H)
      use typy
      use global_objs
      use evapglob
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: H
      
      real(kind=rkind) :: Ts, Ta, rH, Ustar
      integer(kind=rkind) :: pos
      
      pos = getmeteopos()
      
      H = C_air  * (pde(heat_ord)%getval(quadpnt) - meteo4evap(pos)%T_air)/ resH()
      
    end function Hterm
    
    
    function Eterm(quadpnt, layer) result(E)
      use typy
      use global_objs
      use re_constitutive
      use re_globals
      use evap_RE_constitutive
      use pde_objs
      use evapglob
      
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind), intent(in) :: layer
      real(kind=rkind) :: E
      
      real(kind=rkind) :: rs, theta, ths
      integer(kind=ikind) :: pos
      type(integpnt_str) :: quad4atm
      
      theta = vangen(pde(re_ord), layer, quadpnt)
      ths = vgset(layer)%ths
      
      rs = max(-805 + 4140*(ths-theta), 0.0_rkind)
!      rs = 0
      
      pos = getmeteopos()
      quad4atm%type_pnt = "numb"
      quad4atm%this_is_the_value = meteo4evap(pos)%T_air
      
      E = min((dens_satvap(quad4atm)*meteo4evap(pos)%relhum - dens_satvap(quadpnt)*relhumid(quadpnt))/(resH() + rs), 0.0_rkind)
      
    end function Eterm
    
    subroutine ebalance_flux(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      use evap_heat_constitutive
      use printtools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      
      integer(kind=ikind) :: layer, nodeid, i
      type(integpnt_str) :: quadpnt_loc
      
      integer, save :: outfile
      integer :: ierr, D
      logical, save :: opened=.false., header_written = .false.
      logical, dimension(:), allocatable, save :: nodesmarker
      character(len=1), dimension(3) :: xyz = (/"x", "z", "y"/)
      integer(kind=ikind), save :: bc_nds = 0
      real(kind=rkind), dimension(:,:), allocatable, save :: ebalance_vals
      real(kind=rkind) R, H, LE
      
      D = drutes_config%dimen
      
      if (.not. allocated(nodesmarker)) then
        allocate(nodesmarker(nodes%kolik))
        nodesmarker = .false.
      end if

      
      if (present(code)) then
        code = 2
        if (.not. opened) then
          open(newunit=outfile, file="out/ebalance.out", status="replace", action="write", iostat=ierr)
          if (ierr /=0) then
            print *, "error opening file out/ebalance.out"
            ERROR STOP
          end if
          call print_logo(outfile)
          write(unit=outfile, fmt=*) "#----------------------------------------------------------------------------"
          write(unit=outfile, fmt=*) "#----------------------------------------------------------------------------"
          write(unit=outfile, fmt=*) "# node id                  | coordinates: ", xyz(1:D)
          write(unit=outfile, fmt=*) "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
          opened = .true.
        end if
        nodeid = elements%data(el_id, node_order)
        if (.not. nodesmarker(nodeid)) then
          write(unit=outfile, fmt=*) "#", nodeid, nodes%data(nodeid,:)
          nodesmarker(nodeid) = .true.
          bc_nds = bc_nds + 1
        end if
      end if
      
      
      if (present(value)) then
        if (.not. allocated(ebalance_vals)) allocate(ebalance_vals(bc_nds,3))
        if (time > epsilon(time) .and. itcount == 1 .and. bc_nds == 1) then
          if (.not. header_written) then
            write(unit=outfile, fmt=*) "#-----------------------------------------------------------------------------------------"
            write(unit=outfile, fmt=*) "#-----------------------------------------------------------------------------------------"
            write(unit=outfile, fmt=*) "#time                               R                           H                        LE"
            write(unit=outfile, fmt=*) "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
            header_written = .true.
          end if
          i = 1
          if (iter_succesfull) write(unit=outfile, fmt=*) time-dtprev, (/ ( ebalance_vals(i,:), i=1, ubound(ebalance_vals,1) ) /)
        end if
        quadpnt_loc%column = 2
        quadpnt_loc%type_pnt = "ndpt"
        quadpnt_loc%order = elements%data(el_id, node_order)
        layer = elements%material(el_id)
        R = Rnterm(quadpnt_loc, layer)
        H = Hterm(quadpnt_loc)
        LE = latentheat(quadpnt_loc)*Eterm(quadpnt_loc, layer)
        ebalance_vals(bc_nds,1) = R
        ebalance_vals(bc_nds,2) = H
        ebalance_vals(bc_nds,3) = LE
        value = R - H + LE
      end if
      
    
    end subroutine ebalance_flux
    
    
    
    
    subroutine evaporation_bcflux(pde_loc, el_id, node_order, value, code, array) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      use evap_heat_constitutive
      use evap_RE_constitutive
      use evapglob
      use printtools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      !> unused for this model (implementation for Robin boundary)
      real(kind=rkind), dimension(:), intent(out), optional :: array
      real(kind=rkind), dimension(3,3) :: K
      real(kind=rkind), dimension(3) :: gravflux, bcflux

      real(kind=rkind):: dens_wat
      integer(kind=ikind) :: layer, nodeid, D
      type(integpnt_str) :: quadpnt_loc
!      integer, save :: outfile
!      integer :: ierr
!      logical, save :: opened=.false.
!      logical, dimension(:), allocatable, save :: nodesmarker
!      character(len=1), dimension(3) :: xyz = (/"x", "z", "y"/)
!      integer(kind=ikind), save :: bc_nds = 0
!      real(kind=rkind), dimension(:), allocatable, save :: ebalance_vals
      
      
      layer = elements%material(el_id)
!      D = drutes_config%dimen
      
!      if (.not. allocated(nodesmarker)) then
!        allocate(nodesmarker(nodes%kolik))
!        nodesmarker = .false.
!      end if

      if (present(code)) then
        code = 2
!        if (.not. opened) then
!          open(newunit=outfile, file="out/ebalance.out", status="replace", action="write", iostat=ierr)
!          if (ierr /=0) then
!            print *, "error opening file out/ebalance.out"
!            ERROR STOP
!          end if
!          call print_logo(outfile)
!          write(unit=outfile, fmt=*) "#----------------------------------------------------------------------------"
!          write(unit=outfile, fmt=*) "#----------------------------------------------------------------------------"
!          write(unit=outfile, fmt=*) "# node id                  | coordinates: ", xyz(1:D)
!           write(unit=outfile, fmt=*) "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
!          opened = .true.
!        end if
!        nodeid = elements%data(el_id, node_order)
!        if (.not. nodesmarker(nodeid)) then
!          write(unit=outfile, fmt=*) "#", nodeid, nodes%data(nodeid,:)
!          nodesmarker(nodeid) = .true.
!          bc_nds = bc_nds + 1
!        end if
      end if
        
        
          
      
      if (present(value)) then
!        if (.not. allocated(ebalance_vals)) allocate(ebalance
        quadpnt_loc%column = 2
        quadpnt_loc%type_pnt = "ndpt"
        quadpnt_loc%order = elements%data(el_id, node_order)
        call pde_loc%pde_fnc(pde_loc%order)%dispersion(pde_loc, elements%material(el_id), quadpnt_loc, &
                  tensor=K(1:drutes_config%dimen, 1:drutes_config%dimen))
        gravflux(1:drutes_config%dimen) = K(drutes_config%dimen, 1:drutes_config%dimen)*elements%nvect_z(el_id, node_order)
        
        dens_wat = dens_liquid(quadpnt_loc)
        value = Eterm(quadpnt_loc, layer)/dens_wat
        
        select case(drutes_config%dimen)
         case(1)
           value = value + gravflux(1)
         case(2)
           bcflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*value
           bcflux(2) = elements%nvect_z(el_id, node_order)*value
           bcflux = bcflux + gravflux
           value = sqrt(bcflux(1)*bcflux(1) + bcflux(2)*bcflux(2))
       end select 
      end if
      

      
    end subroutine evaporation_bcflux
    
    function albedo_fnc(quadpnt, layer) result(val)
      use typy
      use global_objs
      use evapglob
      use pde_objs
      use globals
      use re_constitutive
      
      type(integpnt_str), intent(in) :: quadpnt
      integer(kind=ikind), intent(in) :: layer
      real(kind=rkind) :: val
      
      integer(kind=ikind) :: i
      real(kind=rkind) :: theta
      
      
      select case(albedo_conf%method)
        case(1)
          do i=albedo_conf%datapos, ubound(albedo_conf%albdat,1) - 1
            if (time >= albedo_conf%albdat(i,1) .and. time < albedo_conf%albdat(i+1,1)) then
              val = albedo_conf%albdat(i,2)
              albedo_conf%datapos = i
              EXIT
            else if ( i == ubound(albedo_conf%albdat,1) - 1 .and. time >= albedo_conf%albdat(i+1,1)) then
              val = albedo_conf%albdat(i+1,2)
              EXIT
            end if
          end do
        case(2)
          theta = vangen(pde(re_ord), layer, quadpnt)
          if (theta < albedo_conf%theta_min) val = albedo_conf%albd_max
          if (theta >= albedo_conf%theta_max) val = albedo_conf%albd_min
          if (theta >=  albedo_conf%theta_min .and. theta < albedo_conf%theta_max) val = albedo_conf%A - theta
        end select
        
      
    
    end function albedo_fnc


end module evapbc4heat
