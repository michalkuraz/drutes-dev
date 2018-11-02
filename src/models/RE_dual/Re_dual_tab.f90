
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

!> \file Re_dual_tab.f90
!>  \brief This module writes tabular versions of main functions to enhance computation. <br>
!! Tabular versions dispersion: mualem_tab; elastisticty: dual_ret; mass= vangen_tab 
!<
module dual_tab
  use typy
  use global_objs
  use dual_globals
  use dual_por
  use dual_coup


  public:: dual_mualem_tab
  public:: vangen_d_tab
  public:: dual_ret_cap_tab
  public:: dual_coupling_neg_tab,dual_coupling_tab
  public:: dual_tabvalues

  real(kind=rkind), dimension(:,:), allocatable, public :: Ktab_dm,watcontab_dm,warecatab_dm,couptab
  real(kind=rkind), dimension(:,:), allocatable, public :: Ktab_df,watcontab_df,warecatab_df
  
  contains 
  
!> Tabular versions of dispersion
  subroutine dual_mualem_tab(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use debug_tools
      use dual_por

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt      
      !> second order tensor of the unsaturated hydraulic conductivity
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor		
      !> relative hydraulic conductivity, (scalar value)
      real(kind=rkind), intent(out), optional :: scalar

      real(kind=rkind) :: h
      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist, tmp
      type(integpnt_str) :: quadpnt_loc 
      
      
      if (present(quadpnt) .and. present(x)) then
		print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
		print *, "exited from Re_dual::mualem_m_tab"
		ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
		print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from Re_dual::mualem_m_tab"
		ERROR stop
      end if
      
     
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	    quadpnt_loc%preproc=.true.
	    h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  		print *, "ERROR: van Genuchten function is a function of a single variable h"
	  		print *, "       your input data has:", ubound(x,1), "variables"
	  		print *, "exited from Re_dual::mualem_m_tab"
	  		ERROR STOP
		end if
		h = x(1)
      end if
      
      if (h<0) then
        if (-h/drutes_config%fnc_discr_length < 0.1*huge(1) ) then
          pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (abs(pos) <= ubound(Ktab_dm,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    select case (pde_loc%mfswitch)
              case("m")
	       tmp = (Ktab_dm(layer,pos+1)-Ktab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + Ktab_dm(layer,pos)
              case("f")
            tmp = (Ktab_df(layer,pos+1)-Ktab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + Ktab_df(layer,pos)
            end select
	  else
	    if (present(quadpnt)) call dual_mualem(pde_loc, layer, quadpnt, scalar = tmp)
	    if (present(x)) call dual_mualem(pde_loc, layer, x=x, scalar = tmp)
	  end if
 	else
 	  if (present(quadpnt)) call dual_mualem(pde_loc, layer, quadpnt, scalar = tmp)
 	  if (present(x)) call dual_mualem(pde_loc, layer, x=x, scalar = tmp)
  	end if 
      else
	if (present(quadpnt)) call dual_mualem(pde_loc, layer, quadpnt, scalar = tmp)
 	if (present(x)) call dual_mualem(pde_loc, layer, x=x, scalar = tmp)
      end if

      if (present(tensor)) then
        select case (pde_loc%mfswitch)
          case("m")
	    tensor = tmp*vgmatrix(layer)%Ks
          case("f")
             tensor = tmp*vgfracture(layer)%Ks
        end select
      end if

      if (present(scalar)) then
	scalar = tmp
      end if
  end subroutine dual_mualem_tab
!> Tabular versions of mass
  function vangen_d_tab(pde_loc, layer, quadpnt, x) result(theta)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use dual_por
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting water content
      real(kind=rkind) :: theta

      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist
      type(integpnt_str) :: quadpnt_loc

      
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from RE_dual::vangen_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from RE_dual::vangen_tab"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if

     if (h<0) then
        if ( h/drutes_config%fnc_discr_length < 0.1*huge(1)) then
          pos = int(-h/drutes_config%fnc_discr_length)+1
          if (abs(pos) <= ubound(watcontab_dm,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    select case (pde_loc%mfswitch)
              case("m")
	        theta = (watcontab_dm(layer,pos+1)-watcontab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + watcontab_dm(layer,pos)
	      case("f")
	        theta = (watcontab_df(layer,pos+1)-watcontab_df(layer,pos))/drutes_config%fnc_discr_length*dist&
	     + watcontab_df(layer,pos)
            end select

          else
            if (present(quadpnt)) theta = vangen_d(pde_loc, layer, quadpnt)
            if (present(x)) theta = vangen_d(pde_loc, layer, x=x)
          end if
        else
          if (present(quadpnt)) theta = vangen_d(pde_loc, layer, quadpnt)
          if (present(x)) theta = vangen_d(pde_loc, layer, x=x)
        end if
    else
      select case (pde_loc%mfswitch)
        case("m")
           theta = vgmatrix(layer)%Ths
        case("f")
           theta = vgfracture(layer)%Ths
        end select
    end if
      

    end function vangen_d_tab
!> Tabular versions of elastiocity
  function dual_ret_cap_tab(pde_loc, layer, quadpnt, x) result(E)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use debug_tools
      use dual_por

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist
      type(integpnt_str) :: quadpnt_loc      
      
   
      
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from re_constitutive::vangen_elast_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_elast_tab"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_elast_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if


     
      if (h<0) then !same as 1104
	if ( h/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	  pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (abs(pos) <= ubound(warecatab_dm,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    select case (pde_loc%mfswitch)
              case("m")
	        E = (warecatab_dm(layer,pos+1)-warecatab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
	         + warecatab_dm(layer,pos)
              case("f")
                E = (warecatab_df(layer,pos+1)-warecatab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
	      + warecatab_df(layer,pos)
            end select
	  else
	    if (present(quadpnt)) E = dual_ret_cap(pde_loc, layer, quadpnt)
	    if (present(x)) E = dual_ret_cap(pde_loc, layer, x=x)
	  end if
	else
	  if (present(quadpnt)) E = dual_ret_cap(pde_loc, layer, quadpnt)
	  if (present(x)) E = dual_ret_cap(pde_loc, layer, x=x)	  
	end if
	
      else
        select case (pde_loc%mfswitch)
          case("m")
            E = vgmatrix(layer)%Ss*exchange(layer)%weightm
          case("f")
            E = vgfracture(layer)%Ss*exchange(layer)%weightf
        end select
	
      end if


  end function dual_ret_cap_tab
!> Tabular versions of coupling term for case(1) and case(2) 
  function dual_coupling_tab(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use dual_por
	  
      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> resulting system elasticity
      real(kind=rkind)::beta,a,gam_par
      real(kind=rkind)::Ks,weightf,weightm,Ksm,Ksf
      real(kind=rkind):: Ka_f,Ka_m,Ka,ex_term,Ka_fm,Ka_mf,Ka_mm,Ka_ff,Ka_fw,Ka_mw
      real(kind=rkind):: hm,hf,one,hw
      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist
      type(integpnt_str) :: quadpnt_loc 

      beta=exchange(layer)%beta
      a=exchange(layer)%a
      gam_par=exchange(layer)%gam_par
      weightm=exchange(layer)%weightm
      weightf=exchange(layer)%weightf
      Ksm=vgmatrix(layer)%KS_local(1)
      Ksf=vgfracture(layer)%KS_local(1)     
      if (present(quadpnt) .and. present(x)) then
	print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
	print *, "exited from RE_dual: dual_coupling_tab"
	ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from RE_dual: dual_coupling_tab"
	ERROR stop
      end if
      
	if (present(quadpnt)) then
	  quadpnt_loc=quadpnt
	  quadpnt_loc%preproc=.true.
	  hm = pde(1)%getval(quadpnt_loc)
	  hf = pde(2)%getval(quadpnt_loc)
    else
	if (ubound(x,1) /=2) then
	  print *, "ERROR: the coupling term is a function of two variables hm and hf"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  ERROR STOP
	end if
    if (ubound(x,1) /=2) then
	  print *, "ERROR: the coupling term is a function of two variables hm and hf"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from RE_dual::dual_coupling_tab"
	  ERROR STOP
	end if
	   hm = x(1)
	   hf = x(2)
    end if

    hw=hm*weightm+hf*weightf
	select case(coup_model)
		case(1:2)
		 Ks=vgexchange(layer)%KS_local(1) 
		if (hm<0) then
		  if ( hm/drutes_config%fnc_discr_length < 0.1*huge(1)) then
			pos = int(-hm/drutes_config%fnc_discr_length)+1
			if (abs(pos) <= ubound(couptab,2)-1) then
			  dist = -hm - (pos - 1)*drutes_config%fnc_discr_length
			  Ka_m = (couptab(layer,pos+1)-couptab(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + couptab(layer,pos)
			else
			  Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
			end if
		  else
			Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
		  end if
		else
		  Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
		end if
	
		if (hf<0) then
		  if ( hf/drutes_config%fnc_discr_length < 0.1*huge(1)) then
			pos = int(-hf/drutes_config%fnc_discr_length)+1
			if (abs(pos) <= ubound(couptab,2)-1) then
			  dist = -hf - (pos - 1)*drutes_config%fnc_discr_length
			  Ka_f = (couptab(layer,pos+1)-couptab(layer,pos))/drutes_config%fnc_discr_length*dist&
			   + couptab(layer,pos)
			else
			  Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
			end if
		  else
			Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
		  end if
		else
		  Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
		end if
    case(4)
		if (hm<0) then
		  if ( hm/drutes_config%fnc_discr_length < 0.1*huge(1)) then
			pos = int(-hm/drutes_config%fnc_discr_length)+1
			if (abs(pos) <= ubound(couptab,2)-1) then
			  dist = -hm - (pos - 1)*drutes_config%fnc_discr_length
			  Ka_mm = (Ktab_dm(layer,pos+1)-Ktab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + Ktab_dm(layer,pos)
			  Ka_fm= (Ktab_df(layer,pos+1)-Ktab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + Ktab_df(layer,pos)
			else
			  call dual_mualem(pde_loc,layer,x=[hm],scalar=Ka_mm)
			  call dual_mualem(pde_loc,layer,x=[hm],scalar=Ka_fm)
			end if
		  else
			call  dual_mualem(pde_loc,layer,x=[hm],scalar=Ka_mm)
			call  dual_mualem(pde_loc,layer,x=[hm],scalar=Ka_fm)
		  end if
		else
			call  dual_mualem(pde_loc,layer,x=[hm],scalar=Ka_mm)
			call  dual_mualem(pde_loc,layer,x=[hm],scalar=Ka_fm)
		end if
    
		if (hf<0) then
		  if ( hf/drutes_config%fnc_discr_length < 0.1*huge(1)) then
			pos = int(-hf/drutes_config%fnc_discr_length)+1
			if (abs(pos) <= ubound(couptab,2)-1) then
			  dist = -hf - (pos - 1)*drutes_config%fnc_discr_length
			  Ka_ff = (Ktab_df(layer,pos+1)-Ktab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + Ktab_df(layer,pos)
			  Ka_mf = (Ktab_dm(layer,pos+1)-Ktab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + Ktab_dm(layer,pos)
			else
			  call dual_mualem(pde_loc,layer,x=[hf],scalar=Ka_mf)
			  call dual_mualem(pde_loc,layer,x=[hf],scalar=Ka_ff)
			end if
		  else
			  call dual_mualem(pde_loc,layer,x=[hf],scalar=Ka_mf)
			  call dual_mualem(pde_loc,layer,x=[hf],scalar=Ka_ff)
		  end if
		else
			call dual_mualem(pde_loc,layer,x=[hf],scalar=Ka_mf)
			call dual_mualem(pde_loc,layer,x=[hf],scalar=Ka_ff)
		end if
		Ka_mm=Ka_mm*Ksm
		Ka_fm=Ka_fm*Ksf
		Ka_mf=Ka_mf*Ksm
		Ka_ff=Ka_ff*Ksf
	case(5)

		if (hw<0) then
		  if ( hw/drutes_config%fnc_discr_length < 0.1*huge(1)) then
			pos = int(-hw/drutes_config%fnc_discr_length)+1
			if (abs(pos) <= ubound(couptab,2)-1) then
			  dist = -hw - (pos - 1)*drutes_config%fnc_discr_length
			  Ka_fw = (Ktab_df(layer,pos+1)-Ktab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + Ktab_df(layer,pos)
			  Ka_mw = (Ktab_dm(layer,pos+1)-Ktab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
			  + Ktab_dm(layer,pos)
			else
			  call dual_mualem(pde_loc,layer,x=[hw],scalar=Ka_mw)
			  call dual_mualem(pde_loc,layer,x=[hw],scalar=Ka_fw)
			end if
		  else
			  call dual_mualem(pde_loc,layer,x=[hw],scalar=Ka_mw)
			  call dual_mualem(pde_loc,layer,x=[hw],scalar=Ka_fw)
		  end if
		else
			call dual_mualem(pde_loc,layer,x=[hw],scalar=Ka_mw)
			call dual_mualem(pde_loc,layer,x=[hw],scalar=Ka_fw)
		end if
		Ka_fw=Ka_fw*Ksf
		Ka_mw=Ka_mw*Ksm
	case(3)
	
	case default
		print *, "ERROR! You have specified an unsupported coupling term model in dual.conf file"
		print *, "the incorrect value specified is:", coup_model
		ERROR stop
  end select
    
    
  select case(coup_model)
	 case(1)
	   Ka=0.5*(Ka_f+Ka_m)*Ks
	 case(2)
	   Ka=(Ka_f*Ka_m)**0.5*Ks
	 case(3)
           Ka=vgexchange(layer)%KS_local(1)
	 case(4)
	   Ka=minval((/Ka_ff,Ka_fm,Ka_mm,Ka_mf/))
	 case(5)
	   Ka=minval((/Ka_fw,Ka_mw/))
   end select
    
    if(abs(hm-hf)>max(hm,hf)*epsilon(hm)) then
      ex_term=beta/a**2*gam_par*Ka
    else
      ex_term=0.0_rkind
    end if
  end function dual_coupling_tab 
  
  function dual_coupling_neg_tab(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use dual_por
      
      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind):: ex_term
      
      if (present(quadpnt)) then
        ex_term=-dual_coupling_tab(pde_loc, layer, quadpnt)
      end if
      if (present(x)) then
        ex_term=-dual_coupling_tab(pde_loc, layer, x=x)
      end if
  end function dual_coupling_neg_tab 
  
! tabular function
 !> creates a table of values of constitutive functions for the Richards equation to be linearly approximated 
  subroutine dual_tabvalues(pde_loc, Kfnc, Cfnc, thetafnc,ex_K_fnc)
      use typy
      use globals
      use pde_objs
      use printtools
      use core_tools
      use dual_por
      class(pde_str), intent(in) :: pde_loc
      
      interface
	subroutine Kfnc(pde_loc, layer, quadpnt, x, tensor, scalar)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in)                           :: layer
	  type(integpnt_str), intent(in), optional                   :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional      :: x
	  real(kind=rkind), intent(out), dimension(:,:), optional   :: tensor
	  real(kind=rkind), intent(out), optional                   :: scalar 
	end subroutine Kfnc
      end interface

      interface
	function Cfnc(pde_loc, layer, quadpnt,  x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function Cfnc
      end interface


      interface
	function thetafnc(pde_loc, layer, quadpnt, x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc	  
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function thetafnc
      end interface

      interface
	function ex_K_fnc(pde_loc, layer, quadpnt, x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc	  
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function ex_K_fnc
      end interface

      integer(kind=ikind) :: i,j, n	
      integer :: l
      integer(kind=ikind) :: maxcalls, counter
      real(kind=rkind) :: dx
	
	
      n = int(maxpress/drutes_config%fnc_discr_length)+1

      drutes_config%fnc_discr_length = 1.0_rkind*maxpress/n
      dx = drutes_config%fnc_discr_length
     
      select case (pde_loc%mfswitch)
        case("m")
          if (.not. allocated(Ktab_dm)) then
            allocate(Ktab_dm(ubound(vgmatrix,1), n))
          end if
          if (.not. allocated(warecatab_dm)) then
            allocate(warecatab_dm(ubound(vgmatrix,1),n))
          end if
          if (.not. allocated(watcontab_dm)) then
            allocate(watcontab_dm(ubound(vgmatrix,1), n))
          end if
          if (.not. allocated(couptab)) then
            allocate(couptab(ubound(vgmatrix,1), n))
          end if
       case("f")
         if (.not. allocated(Ktab_df)) then
           allocate(Ktab_df(ubound(vgfracture,1), n))
         end if
         if (.not. allocated(warecatab_df)) then
           allocate(warecatab_df(ubound(vgfracture,1),n))
         end if
         if (.not. allocated(watcontab_df)) then
           allocate(watcontab_df(ubound(vgfracture,1), n))
         end if
      end select
      
	maxcalls = ubound(vgmatrix,1)*n
	counter = maxcalls
     select case (pde_loc%mfswitch)
       case("m")
	call write_log(text="creating constitutive function table for matrix")
	do i=1, ubound(vgmatrix,1)
	  do j=1, n
	    if (this_image() == 1) then
	      counter = counter - 1
	      l = 100*(maxcalls - counter)/maxcalls
	      call progressbar(l)
	    end if
	    call Kfnc(pde_loc,i, x=(/-(j-1)*dx/), scalar=Ktab_dm(i,j))
	    warecatab_dm(i,j) = Cfnc(pde_loc, i, x=(/-(j-1)*dx/))
	    watcontab_dm(i,j) = thetafnc(pde_loc, i, x=(/-(j-1)*dx/))
	    couptab(i,j) = ex_K_fnc(pde_loc, i, x=(/-(j-1)*dx/))
	  end do
	end do
      case("f")
      	call write_log(text="creating constitutive function table for fracture")
          do i=1, ubound(vgfracture,1)
	  do j=1, n
	    if (this_image() == 1) then
	      counter = counter - 1
	      l = 100*(maxcalls - counter)/maxcalls
	      call progressbar(l)
	    end if
	    call Kfnc(pde_loc,i, x=(/-(j-1)*dx/), scalar=Ktab_df(i,j))
	    warecatab_df(i,j) = Cfnc(pde_loc, i, x=(/-(j-1)*dx/))
	    watcontab_df(i,j) = thetafnc(pde_loc, i, x=(/-(j-1)*dx/))
	  end do
	end do
      end select
    end subroutine dual_tabvalues

end module dual_tab
