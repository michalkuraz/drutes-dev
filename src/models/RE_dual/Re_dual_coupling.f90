!> This module contains variants of coupling terms of the dual permeability model.
!! Boundary gets assigned mualem parameters in case(2) and (3).
!! Case(1) The exchange term is calculated based on arithmetic mean of \f$Ka=\frac{Ka(hm)+Ka(hf)}{2}\f$. <br>
!! Case(2) The exchange term is calculated based on geometric mean of \f$Ka=\sqrt{Ka(hm)*Ka(hf)}\f$. <br>
!! Case(3) The exchange boundary term is constant.<br>
!! Case (4) and (5) don't require further parameters. <br>
!! Case(4) The exchange term is the minimum of \f$Ka=min(Km(hm),Kf(hf),Km(hf),Kf(hf))\f$ <br>
!! Case(5) The exchange term is the minimum of \f$Ka=min(Km(hm*weight_m+hf*weight_f),Kf(hm*weight_m+hf*weight_f))\f$.
!<
module dual_coup

  use typy
  use global_objs
  use dual_globals
  use dual_por

  public:: dual_coupling,dual_coupling_neg
  public:: dual_coup_min,dual_coup_min_neg

  public:: dual_coupling_K
  contains 


!> Coupling defined for case(1), case (2) and case (3)
!! Case(1) The exchange term is calculated based on arithmetic mean of \f$Ka=\frac{Ka(hm)+Ka(hf)}{2}\f$. <br>
!! Case(2) The exchange term is calculated based on geometric mean of \f$Ka=\sqrt{Ka(hm)*Ka(hf)}\f$. <br>
!! Case(3) The exchange boundary term is constant.<br>
  function dual_coupling(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use global_objs
      use pde_objs
      use dual_globals
      use Re_dual_reader
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> vg and ex parameters, later from conf file
      real(kind=rkind)::beta,a,gam_par
      real(kind=rkind)::n,m,alpha,Ks
      real(kind=rkind)				  :: Ka_f,Ka_m,Ka,ex_term
      real(kind=rkind)				  :: hm,hf,one
      type(integpnt_str) :: quadpnt_loc     
        
       
            
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.

	hm = pde(1)%getval(quadpnt_loc)
	hf = pde(2)%getval(quadpnt_loc)
      else
	    if (ubound(x,1) /=2) then
	      print *, "ERROR: exchange term requires two variables h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      ERROR STOP
	    end if
      	if (ubound(x,1) /=2) then
	      print *, "ERROR: exchange term requires two variables h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      print *, "exited from RE_dual::dual_coupling"
	      ERROR STOP
	   end if
	   hm = x(1)
	   hf = x(2)
     end if
     
     alpha=vgexchange(layer)%alpha
     n=vgexchange(layer)%n
     m=vgexchange(layer)%m
     Ks=vgexchange(layer)%KS_local(1)
     one=1.0_rkind   
     if(hf<0) then
       Ka_f=(one-abs(alpha*hf)**(n*m)*(one+abs(alpha*hf)**n)**(-m))**2/(one+abs(alpha*hf)**n)**(m/2)
     else
       Ka_f=1.0_rkind
     end if
     if(hm<0) then
       Ka_m=(one-abs(alpha*hm)**(n*m)*(one+abs(alpha*hm)**n)**(-m))**2/(one+abs(alpha*hm)**n)**(m/2)
     else
       Ka_m=1.0_rkind
     end if  

     select case(coup_model)
	   case(1)
     	 Ka=0.5*(Ka_f+Ka_m)*Ks
       case(2)
         Ka=(Ka_f*Ka_m)**0.5_rkind*Ks
       case(3)
         Ka=vgexchange(layer)%KS_local(1)
    end select
    
     beta=exchange(layer)%beta
     a=exchange(layer)%a
     gam_par=exchange(layer)%gam_par
     if(abs(hm-hf)>max(hm,hf)*epsilon(hm)) then
       ex_term=beta/a**2*gam_par*Ka
     else
       ex_term=0.0_rkind
     end if

  end function dual_coupling
  !> neg coupling defined for case(1), case (2) and case (3) for coupled matrix, i.e. matrix part in fracture and fracture part in matrix. 
  function dual_coupling_neg(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use global_objs
      use pde_objs
      use dual_globals
      use Re_dual_reader
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> vg and ex parameters, later from conf file
      real(kind=rkind)   ::ex_term
          
            
      if (present(quadpnt)) then
        ex_term=-dual_coupling(pde_loc, layer, quadpnt)
      end if
      if (present(x)) then
        ex_term=-dual_coupling(pde_loc, layer, x=x)
      end if

     
  end function dual_coupling_neg

  !> The following coupling schemes use the minimum value of case(4) \f$Ka=min(Km(hm),Kf(hf),Km(hf),Kf(hf))\f$ and case(5) \f$Ka=min(Km(hm*weight_m+hf*weight_f),Kf(hm*weight_m+hf*weight_f))\f$. <br>
  !! This has the advantage that no extra parameters need to be defined describing the hydraulic properties at the boundary<br>
  !! This reduces the number of parameters to be estimated by 3: Kas, a, n <br>
  !! This also assures that the exchange term is always smaller than or equal the fracture and matrix conductivity.
  !<
    function dual_coup_min(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use global_objs
      use pde_objs
      use dual_globals
      use Re_dual_reader
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> vg and ex parameters, later from conf file
      real(kind=rkind)::beta,a,gam_par
      real(kind=rkind)::weightf,weightm,Ks,Ksm,Ksf
      real(kind=rkind)				  :: Ka_ff,Ka_mf,Ka_mm,Ka_fm,Ka,ex_term,Ka_m,Ka_f
      real(kind=rkind)				  :: hm,hf
      type(integpnt_str) :: quadpnt_loc     
        
      Ksm=vgmatrix(layer)%KS_local(1)
      Ksf=vgfracture(layer)%KS_local(1) 
            
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	quadpnt_loc%column=1
	hm = pde(1)%getval(quadpnt_loc)
	hf = pde(2)%getval(quadpnt_loc)
      else
	    if (ubound(x,1) /=2) then
	      print *, "ERROR: exchange term requires two variables h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      ERROR STOP
	    end if
      	if (ubound(x,1) /=2) then
	      print *, "ERROR: exchange term requires two variables h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      print *, "exited from RE_dual::dual_coupling"
	      ERROR STOP
	   end if
	   hm = x(1)
	   hf = x(2)
     end if
          
     select case(coup_model)
	   case(4)
		 call dual_mualem(pde_loc, layer, x=[hm],scalar=Ka_mm)
		 call dual_mualem(pde_loc, layer, x=[hf],scalar=Ka_mf)
		 call dual_mualem(pde_loc, layer, x=[hm],scalar=Ka_fm)
		 call dual_mualem(pde_loc, layer, x=[hf],scalar=Ka_ff)
		 Ka_ff=Ka_ff*Ksf
		 Ka_mf=Ka_mf*Ksm
		 Ka_mm=Ka_mm*Ksm
		 Ka_fm=Ka_fm*Ksf
		 Ka=minval((/Ka_ff,Ka_mf,Ka_mm,Ka_fm/))

       case(5)
         weightm=exchange(layer)%weightm
     	 weightf=exchange(layer)%weightf
     	 call dual_mualem(pde_loc, layer, x=[hm*weightm+hf*weightf],scalar=Ka_m)
     	 call dual_mualem(pde_loc, layer, x=[hm*weightm+hf*weightf],scalar=Ka_f)
     	 Ka_m=Ka_m*Ksm
     	 Ka_f=Ka_f*Ksf
   		 Ka=minval((/Ka_f,Ka_m/))
    end select
    
     beta=exchange(layer)%beta
     a=exchange(layer)%a
     gam_par=exchange(layer)%gam_par
     
     if(abs(hm-hf)>max(hm,hf)*epsilon(hm)) then
       ex_term=beta/a**2*gam_par*Ka
     else
       ex_term=0.0_rkind
     end if

  end function dual_coup_min
  
  !> neg coupling defined for case(4) and case (5) for coupled matrix, i.e. matrix part in fracture and fracture part in matrix. 
  function dual_coup_min_neg(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use global_objs
      use pde_objs
      use dual_globals
      use Re_dual_reader
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> vg and ex parameters

      real(kind=rkind)   ::ex_term
          
            
      if (present(quadpnt)) then
        ex_term=-dual_coup_min(pde_loc, layer, quadpnt)
      end if
      if (present(x)) then
        ex_term=-dual_coup_min(pde_loc, layer, x=x)
      end if

     
  end function dual_coup_min_neg
  

 !> coupling conductivity function for tab values 
  function dual_coupling_K(pde_loc,layer, quadpnt, x) result(Ka_c)
    use typy
    use global_objs
    use pde_objs
    use dual_globals
    use Re_dual_reader
    
    class(pde_str), intent(in) :: pde_loc
    !> value of the nonlinear function
    real(kind=rkind), dimension(:), intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    !> material ID
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind)::Ka_c
    real(kind=rkind)::n,m,alpha
    real(kind=rkind):: one,h
    type(integpnt_str) :: quadpnt_loc
          
     if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
     else
	 if (ubound(x,1) /=1) then
	   print *, "ERROR: van Genuchten function is a function of a single variable h"
	   print *, "       your input data has:", ubound(x,1), "variables"
	   print *, "exited from RE_dual::dual_coupling_K"
	   ERROR STOP
	 end if
	   h = x(1)
     end if  
    alpha=vgexchange(layer)%alpha
    n=vgexchange(layer)%n
    m=vgexchange(layer)%m
    one=1.0_rkind 
    if(h<0) then
      Ka_c=(one-(alpha*abs(h))**(n*m)*(one+(alpha*abs(h))**n)**(-m))**2/(one+(alpha*abs(h))**n)**(m/2)
    else
      Ka_c=1.0_rkind
    end if
  end function dual_coupling_K

end module dual_coup