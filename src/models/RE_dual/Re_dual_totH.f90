module dual_por
	use typy
    use global_objs
    use dual_globals


public:: dual_mualemm,dual_mualemf, dual_ret_capf, dual_ret_capm, dual_coupling
public:: vangen_d_f, vangen_d_m, dual_coupling_f
public:: dual_inicond_f,dual_inicond_m
public:: darcy_law_d
public:: dual_mualem_m_tab,dual_mualem_f_tab
public:: vangen_d_m_tab, vangen_d_f_tab
public:: dual_ret_capf_tab, dual_ret_capm_tab
public:: dual_coupling_f_tab,dual_coupling_tab,dual_coupling_K
public:: dual_tabvalues
public :: getval_retot_dual

real(kind=rkind), dimension(:,:), pointer, public :: Ktab_dm,watcontab_dm,warecatab_dm,couptab
real(kind=rkind), dimension(:,:), pointer, public :: Ktab_df,watcontab_df,warecatab_df
real(kind=rkind), dimension(:,:), allocatable, target, public ::Ktabd_tget_m,watcontabd_tget_m
real(kind=rkind), dimension(:,:), allocatable, target, public ::warecatabd_tget_m
real(kind=rkind), dimension(:,:), allocatable, target, public ::Ktabd_tget_f,watcontabd_tget_f
real(kind=rkind), dimension(:,:), allocatable, target, public ::warecatabd_tget_f
real(kind=rkind), dimension(:,:), allocatable, target, public ::couptab_all
contains 

  !> specific function for Richards equation in H-form (total hydraulic head form), replaces pde_objs::getvalp1 in order to distinguish between H and h 
 function getval_retot_dual(pde_loc, quadpnt) result(val)
    use typy
    use pde_objs
    use geom_tools
    use dual_globals

    class(pde_str), intent(in) :: pde_loc
    type(integpnt_str), intent(in) :: quadpnt
    real(kind=rkind) :: val

    real(kind=rkind), dimension(3) :: xyz
    integer(kind=ikind) :: D

    D = drutes_config%dimen
     
    call getcoor(quadpnt, xyz(1:D))


    if (get_pressh_vals) then

      val = getvalp1(pde_loc, quadpnt) 

      get_pressh_vals = .false.
      
    else
        val = getvalp1(pde_loc, quadpnt) - xyz(D)
    end if
        
      
 end function getval_retot_dual

 subroutine dual_inicond_m(pde_loc)
   use typy
   use globals
   use global_objs
   use pde_objs
   use dual_globals
     class(pde_str), intent(in out) :: pde_loc
     integer(kind=ikind) :: i, j, k,l, m, layer
     real(kind=rkind) :: value
     
      do i=1, elements%kolik
	    layer = elements%material(i,1)
	    do j=1, ubound(elements%data,2)
	      k = elements%data(i,j)
          pde_loc%solution(k)=vgmatrix(layer)%initcond
	    end do   
      end do
      
  end subroutine dual_inicond_m
  
   subroutine dual_inicond_f(pde_loc)
   use typy
   use globals
   use global_objs
   use pde_objs
   use dual_globals
     class(pde_str), intent(in out) :: pde_loc
     integer(kind=ikind) :: i, j, k,l, m, layer
     real(kind=rkind) :: value
     
      do i=1, elements%kolik
	    layer = elements%material(i,1)
	    do j=1, ubound(elements%data,2)
	      k = elements%data(i,j)
          pde_loc%solution(k)=vgfracture(layer)%initcond
	    end do   
      end do
      
  end subroutine dual_inicond_f
 
subroutine dual_mualemm(pde_loc, layer, quadpnt, x, tensor, scalar)
  use typy
  use global_objs
  use pde_objs
  use globals
  use debug_tools
  use dual_globals
  use Re_dual_reader
      
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
  !> vg parameters, later from conf file
  real(kind=rkind)::n,m,alpha, weight       
  real(kind=rkind) :: h,Kr,one
  
  get_pressh_vals = .true.
      
    if (present(quadpnt)) then
	  h = pde_loc%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
	    print *, "ERROR: van Genuchten function is a function of a single variable h"
	    print *, "       your input data has:", ubound(x,1), "variables"
	    print *, "exited from RE_dual::dual_mualemm"
	    ERROR STOP
	  end if
	h = x(1)
    end if
!     print *, "h in matrix"
!     print *, h
    alpha=vgmatrix(layer)%alpha
    n=vgmatrix(layer)%n
    m=vgmatrix(layer)%m
    weight=-exchange(layer)%weightm


    one=1.0_rkind    
    if(h < 0.0_rkind) then
      Kr=(one-(-alpha*h)**(n*m)*(one+(-alpha*h)**n)**(-m))**2/(one+(-alpha*h)**n)**(m/2)
    else
      Kr=1.0_rkind
    end if
        
    if (present(tensor)) then
		tensor=vgmatrix(layer)%KS*Kr*weight
    end if
      
    if (present(scalar)) then
      scalar=Kr*weight
    end if
!      print*, " weighted hydraulic conductivty of matrix"
!     print *, Kr
  end subroutine dual_mualemm 
     
subroutine dual_mualemf(pde_loc, layer, quadpnt, x, tensor, scalar)
  use typy
  use global_objs
  use pde_objs
  use globals
  use debug_tools
  use dual_globals
  use Re_dual_reader
       
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
  !> vg parameters, later from conf file
  real(kind=rkind)::n,m,alpha,weight
         
  real(kind=rkind) :: h,Kr,one
  
  get_pressh_vals = .true.
      
    if (present(quadpnt)) then
	  h = pde_loc%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
	    print *, "ERROR: van Genuchten function is a function of a single variable h"
	    print *, "       your input data has:", ubound(x,1), "variables"
	    print *, "exited from RE_dual::mualemf"
	    ERROR STOP
	  end if
	h = x(1)
    end if

    alpha=vgfracture(layer)%alpha
    n=vgfracture(layer)%n
    m=vgfracture(layer)%m
    weight=exchange(layer)%weightf

    one=1.0_rkind    
    
    if(h < 0.0_rkind) then
      Kr=(one-(-alpha*h)**(n*m)*(one+(-alpha*h)**n)**(-m))**2/(one+(-alpha*h)**n)**(m/2)
    else
      Kr=1.0_rkind
    end if
    
    if (present(tensor)) then
		tensor=vgfracture(layer)%KS*Kr*(-weight)
    end if
      
    if (present(scalar)) then
      scalar=Kr*(-weight)
    end if
!      print*, " weighted hydraulic conductivty of fracture"
!      print*, Kr
  end subroutine dual_mualemf
   
  function dual_ret_capm(pde_loc,layer,quadpnt,x) result(E)
      use typy
      use pde_objs
      use core_tools
      use dual_globals
      use Re_dual_reader
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> vg parameters, later from conf file
      real(kind=rkind)::n,m,alpha,thetaS,thetaR,weight
      real(kind=rkind)::E,h,C
      
      
      get_pressh_vals = .true.
      
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
	      print *, "exited from RE_dual::dual_ret_capm"
	      ERROR STOP
	    end if
	    h = x(1)
      end if
     thetaS=vgmatrix(layer)%ThS
     thetaR=vgmatrix(layer)%ThR
     alpha=vgmatrix(layer)%alpha
     n=vgmatrix(layer)%n
     m=vgmatrix(layer)%m
     weight=-exchange(layer)%weightm
     

     
     if(h<0.0_rkind) then
       C=(thetaS-thetaR)*alpha*m*n*abs(alpha*h)**(n-1)*(abs(alpha*h)**n+1)**(-m-1)
     else
       E=vgmatrix(layer)%SS*weight
       RETURN
     end if
     
      E=(C+vangen_d_m(pde_loc,layer,x=(/h/))/thetaS)*weight
!      print*,"weighted elasticity retention of matrix"
!       print*, E
  end function dual_ret_capm  
  
 function dual_ret_capf(pde_loc,layer,quadpnt,x) result(E)
      use typy
      use pde_objs
      use core_tools
      use dual_globals
      use Re_dual_reader
      class(pde_str), intent(in) :: pde_loc 
      !> pressure head
      integer(kind=ikind), intent(in) :: layer
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> vg parameters
      real(kind=rkind)::n,m,alpha,thetaS,thetaR,weight
      real(kind=rkind)::E,h,C
      
      
      get_pressh_vals = .true.
      
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
	      print *, "exited from RE_dual::dual_ret_capf"
	      ERROR STOP
	    end if
	    h = x(1)
      end if
     thetaS=vgfracture(layer)%ThS
     thetaR=vgfracture(layer)%ThR
     alpha=vgfracture(layer)%alpha
     n=vgfracture(layer)%n
     m=vgfracture(layer)%m
     weight=exchange(layer)%weightf
     
     if(h<0.0_rkind) then
       C=(thetaS-thetaR)*alpha*m*n*abs(alpha*h)**(n-1)*(abs(alpha*h)**n+1)**(-m-1)
     else
       E=vgfracture(layer)%SS*(-weight)
       RETURN
     end if     
      E=(C+vangen_d_f(pde_loc,layer,x=(/h/))/thetaS)*(-weight)
!      print*,"weighted elasticity retention of fracture "
!      print*, E
  end function dual_ret_capf
  
 function vangen_d_m(pde_loc,layer,quadpnt,x) result(theta)
   use typy
   use dual_globals
   use pde_objs
   class(pde_str), intent(in) :: pde_loc
   integer(kind=ikind), intent(in) :: layer
   !> pressure head
   real(kind=rkind), intent(in), dimension(:), optional :: x
   !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
   type(integpnt_str), intent(in), optional :: quadpnt
   real(kind=rkind) :: h
   !> resulting water content
   real(kind=rkind) :: theta

   real(kind=rkind) :: a,n,m,ths,thr, theta_e

      get_pressh_vals = .true.
     if (present(quadpnt)) then
	   h = pde_loc%getval(quadpnt)
     else
	 if (ubound(x,1) /=1) then
	   print *, "ERROR: van Genuchten function is a function of a single variable h"
	   print *, "       your input data has:", ubound(x,1), "variables"
	   print *, "exited from RE_dual::vangen_d_m"
	   ERROR STOP
	 end if
	   h = x(1)
     end if    
      
      a = vgmatrix(layer)%alpha
      n = vgmatrix(layer)%n
      m = vgmatrix(layer)%m
      Ths = vgmatrix(layer)%Ths
      Thr = vgmatrix(layer)%Thr
    
      if (h >=0.0_rkind) then
        theta = vgmatrix(layer)%Ths
        RETURN
      else
	    theta_e = 1/(1+(a*(abs(h)))**n)**m
	    theta = theta_e*(Ths-Thr)+Thr
      end if
 end function vangen_d_m
   
 function vangen_d_f(pde_loc,layer,quadpnt,x) result(theta)
   use typy
   use dual_globals
   use pde_objs
   class(pde_str), intent(in) :: pde_loc
   integer(kind=ikind), intent(in) :: layer
   !> pressure head
   real(kind=rkind), intent(in), dimension(:), optional :: x
   !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
   type(integpnt_str), intent(in), optional :: quadpnt
   real(kind=rkind) :: h
   !> resulting water content
   real(kind=rkind) :: theta,weight

   real(kind=rkind) :: a,n,m,ths,thr, theta_e

   
   get_pressh_vals = .true.
      
     if (present(quadpnt)) then
	   h = pde_loc%getval(quadpnt)
     else
	 if (ubound(x,1) /=1) then
	   print *, "ERROR: van Genuchten function is a function of a single variable h"
	   print *, "       your input data has:", ubound(x,1), "variables"
	   print *, "exited from RE_dual::vangen_d_f"
	   ERROR STOP
	 end if
	   h = x(1)
     end if    
      
      a = vgfracture(layer)%alpha
      n = vgfracture(layer)%n
      m = vgfracture(layer)%m
      Ths = vgfracture(layer)%Ths
      Thr = vgfracture(layer)%Thr
      if (h >=0.0_rkind) then
	    theta = vgfracture(layer)%Ths
 	    RETURN
      else
	    theta_e = 1/(1+(a*(abs(h)))**n)**m
	    theta = theta_e*(Ths-Thr)+Thr
      end if
 end function vangen_d_f      
 
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
      
        
       get_pressh_vals = .true.
            
      if (present(quadpnt)) then
	    hm = pde(1)%getval(quadpnt)
	    hf = pde(2)%getval(quadpnt)
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
     Ka_f=(one-(-alpha*hf)**(n*m)*(one+(-alpha*hf)**n)**(-m))**2/(one+(-alpha*hf)**n)**(m/2)
     Ka_m=(one-(-alpha*hm)**(n*m)*(one+(-alpha*hm)**n)**(-m))**2/(one+(-alpha*hm)**n)**(m/2)
     Ka=0.5*(Ka_f+Ka_m)*Ks
     beta=exchange(layer)%beta
     a=exchange(layer)%a
     gam_par=exchange(layer)%gam_par
     if(hf /= hm) then
       ex_term=beta/a**2*gam_par*Ka*(hf-hm)
     else
       ex_term=0.0_rkind
     end if
   !  print*, 'exchange term matrix'
   !  print*, ex_term
  end function dual_coupling
  
  function dual_coupling_f(pde_loc, layer, quadpnt, x) result(ex_term)
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
      real(kind=rkind)::n,m,alpha, Ks
      real(kind=rkind)				  :: Ka_f,Ka_m,Ka,ex_term
      real(kind=rkind)				  :: hm,hf,one
      
       get_pressh_vals = .true.
            
      if (present(quadpnt)) then
	    hm = pde(1)%getval(quadpnt)
	    hf = pde(2)%getval(quadpnt)
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
     Ks=vgexchange(layer)%Ks_local(1)
     beta=exchange(layer)%beta
     a=exchange(layer)%a
     gam_par=exchange(layer)%gam_par
    
     one=1.0_rkind   
     Ka_f=(one-(-alpha*hf)**(n*m)*(one+(-alpha*hf)**n)**(-m))**2/(one+(-alpha*hf)**n)**(m/2)
     Ka_m=(one-(-alpha*hm)**(n*m)*(one+(-alpha*hm)**n)**(-m))**2/(one+(-alpha*hm)**n)**(m/2)
     Ka=0.5*(Ka_f+Ka_m)*Ks
     if(hf /= hm) then
       ex_term=-beta/a**2*gam_par*Ka*(hf-hm)
     else
       ex_term=0
     end if
     
!print*, 'exchange term fracture'
   !  print*, ex_term
  end function dual_coupling_f
  
  subroutine darcy_law_d(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length

      real(kind=rkind), dimension(3,3)  :: K
      integer                           :: D
      integer(kind=ikind)               :: i
      integer(kind=ikind), dimension(3) :: nablaz
      real(kind=rkind), dimension(3)  :: gradH
      real(kind=rkind), dimension(3)  :: vct
      real(kind=rkind) :: h
      real(kind=rkind), dimension(:), allocatable :: gradient
      
      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
	print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
	ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::darcy_law"
	ERROR stop
      end if
      
      if (present(quadpnt)) then
	call pde_loc%getgrad(quadpnt, gradient)
	get_pressh_vals = .true.
	h = pde_loc%getval(quadpnt)
      else
        if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::darcy_law"
	  ERROR STOP
	end if
	h = x(1)
	allocate(gradient(ubound(grad,1)))
	gradient = grad
      end if
      
      D = drutes_config%dimen

      nablaz = 0
      nablaz(D) = 1
      
      gradH(1:D) = gradient(1:D) + nablaz(1:D)

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

    end subroutine darcy_law_d
  
  function dual_coupling_K(pde_loc,layer, quadpnt, x) result(Ka_c)
    use typy
    use global_objs
    use pde_objs
    use dual_globals
    use Re_dual_reader
    
    class(pde_str), intent(in) :: pde_loc
    integer(kind=ikind), intent(in) :: layer
    real(kind=rkind), dimension(:),intent(in), optional    :: x
    !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
    type(integpnt_str), intent(in), optional :: quadpnt
    real(kind=rkind)::Ka_c
    real(kind=rkind)::n,m,alpha
    real(kind=rkind):: one,h
          
     if (present(quadpnt)) then
           get_pressh_vals = .true.
	   h = pde_loc%getval(quadpnt)
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
    Ka_c=(one-(-alpha*h)**(n*m)*(one+(-alpha*h)**n)**(-m))**2/(one+(-alpha*h)**n)**(m/2)
  end function dual_coupling_K
  ! Tabular versions dispersion: mualem_tab; elastisticty: dual_ret; mass= vangen_tab 
  
 subroutine dual_mualem_m_tab(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use debug_tools

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
      
      get_pressh_vals = .true.
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
		h = pde_loc%getval(quadpnt)
      else
      	if (ubound(x,1) /=1) then
	  		print *, "ERROR: van Genuchten function is a function of a single variable h"
	  		print *, "       your input data has:", ubound(x,1), "variables"
	  		print *, "exited from re_constitutive::mualem_tab"
	  		ERROR STOP
		end if
		h = x(1)
      end if
      
      


      
 if (h<0) then
   if (-h/drutes_config%fnc_discr_length < 0.1*huge(1) ) then
      pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (pos <= ubound(Ktab_dm,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    tmp = (Ktab_dm(layer,pos+1)-Ktab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + Ktab_dm(layer,pos)
	  else
	  if (present(quadpnt)) call dual_mualemm(pde_loc, layer, quadpnt, scalar = tmp)
	  if (present(x)) call dual_mualemm(pde_loc, layer, x=x, scalar = tmp)
	  ! if (.not. tabwarning) then
! 	    call write_log(trim(tabmsg))
! 	    tabwarning = .true.
! 	  end if
	end if
 else
!  if (.not. intwarning) then
!    call intoverflow()
!    intwarning = .true.
!  end if
 
 if (present(quadpnt)) call dual_mualemm(pde_loc, layer, quadpnt, scalar = tmp)
 if (present(x)) call dual_mualemm(pde_loc, layer, x=x, scalar = tmp)
  end if 
 else
	tmp = 1
  end if

  if (present(tensor)) then
	tensor = tmp* vgmatrix(layer)%Ks
  end if

  if (present(scalar)) then
	scalar = tmp
  end if
end subroutine dual_mualem_m_tab

subroutine dual_mualem_f_tab(pde_loc, layer, quadpnt,  x, tensor, scalar)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use debug_tools

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
      
      get_pressh_vals = .true.
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
		h = pde_loc%getval(quadpnt)
      else
      	if (ubound(x,1) /=1) then
	  		print *, "ERROR: van Genuchten function is a function of a single variable h"
	  		print *, "       your input data has:", ubound(x,1), "variables"
	  		print *, "exited from re_constitutive::mualem_tab"
	  		ERROR STOP
		end if
		h = x(1)
      end if
      
      


      
 if (h<0) then
   if (-h/drutes_config%fnc_discr_length < 0.1*huge(1) ) then
      pos = int(-h/drutes_config%fnc_discr_length)+1
	  if (pos <= ubound(Ktab_df,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    tmp = (Ktab_df(layer,pos+1)-Ktab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + Ktab_df(layer,pos)
	  else
	  if (present(quadpnt)) call dual_mualemf(pde_loc, layer, quadpnt, scalar = tmp)
	  if (present(x)) call dual_mualemf(pde_loc, layer, x=x, scalar = tmp)
	  ! if (.not. tabwarning) then
! 	    call write_log(trim(tabmsg))
! 	    tabwarning = .true.
! 	  end if
	end if
 else
!  if (.not. intwarning) then
!    call intoverflow()
!    intwarning = .true.
!  end if
 
 if (present(quadpnt)) call dual_mualemf(pde_loc, layer, quadpnt, scalar = tmp)
 if (present(x)) call dual_mualemf(pde_loc, layer, x=x, scalar = tmp)
  end if 
 else
	tmp = 1
  end if

  if (present(tensor)) then
	tensor = tmp* vgfracture(layer)%Ks
  end if

  if (present(scalar)) then
	scalar = tmp
  end if
end subroutine dual_mualem_f_tab    

function vangen_d_m_tab(pde_loc, layer, quadpnt, x) result(theta)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      
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

      get_pressh_vals = .true.
      
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
	h = pde_loc%getval(quadpnt)
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
	  if (pos <= ubound(watcontab_dm,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    theta = (watcontab_dm(layer,pos+1)-watcontab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + watcontab_dm(layer,pos)
	  else
	    if (present(quadpnt)) theta = vangen_d_m(pde_loc, layer, quadpnt)
	    if (present(x)) theta = vangen_d_m(pde_loc, layer, x=x)
	   !  if (.not. tabwarning) then
! 	      call write_log(trim(tabmsg))
! 	      tabwarning = .true.
! 	    end if
	  end if
	else
	
! 	  if (.not. intwarning) then
! 	    call intoverflow()
! 	    intwarning = .true.
! 	  end if
	
	  if (present(quadpnt)) theta = vangen_d_m(pde_loc, layer, quadpnt)
	  if (present(x)) theta = vangen_d_m(pde_loc, layer, x=x)
	  
	end if
	
      else
	theta = vgmatrix(layer)%Ths	
      end if
      
      

    end function vangen_d_m_tab

function vangen_d_f_tab(pde_loc, layer, quadpnt, x) result(theta)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      
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

      get_pressh_vals = .true.
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
	h = pde_loc%getval(quadpnt)
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
	  if (pos <= ubound(watcontab_df,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    theta = (watcontab_df(layer,pos+1)-watcontab_df(layer,pos))/drutes_config%fnc_discr_length*dist&
	     + watcontab_df(layer,pos)
	  else
	    if (present(quadpnt)) theta = vangen_d_f(pde_loc, layer, quadpnt)
	    if (present(x)) theta = vangen_d_f(pde_loc, layer, x=x)
	   !  if (.not. tabwarning) then
! 	      call write_log(trim(tabmsg))
! 	      tabwarning = .true.
! 	    end if
	  end if
	else
	
! 	  if (.not. intwarning) then
! 	    call intoverflow()
! 	    intwarning = .true.
! 	  end if
	
	  if (present(quadpnt)) theta = vangen_d_f(pde_loc, layer, quadpnt)
	  if (present(x)) theta = vangen_d_f(pde_loc, layer, x=x)	  
	  end if
	
      else
	  theta = vgfracture(layer)%Ths	
    end if
         

end function vangen_d_f_tab
 
function dual_ret_capm_tab(pde_loc, layer, quadpnt, x) result(E)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use debug_tools

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
      
   
      get_pressh_vals = .true.
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
	h = pde_loc%getval(quadpnt)
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
	  if (pos <= ubound(warecatab_dm,2)-1) then
	    dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	    E = (warecatab_dm(layer,pos+1)-warecatab_dm(layer,pos))/drutes_config%fnc_discr_length*dist &
	    + warecatab_dm(layer,pos)
	  else
	    if (present(quadpnt)) E = dual_ret_capm(pde_loc, layer, quadpnt)
	    if (present(x)) E = dual_ret_capm(pde_loc, layer, x=x)
	  end if
	else
	  if (present(quadpnt)) E = dual_ret_capm(pde_loc, layer, quadpnt)
	  if (present(x)) E = dual_ret_capm(pde_loc, layer, x=x)	  
	end if
	
      else
	E = vgmatrix(layer)%Ss	
      end if


    end function dual_ret_capm_tab
 
function dual_ret_capf_tab(pde_loc, layer, quadpnt, x) result(E)
      use typy
      use dual_globals
      use pde_objs
      use core_tools
      use debug_tools

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
      
   
      get_pressh_vals = .true.
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
	h = pde_loc%getval(quadpnt)
      else
      	if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from re_constitutive::vangen_elast_tab"
	  ERROR STOP
	end if
	h = x(1)
      end if

     if (h<0) then
	  if ( h/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	    pos = int(-h/drutes_config%fnc_discr_length)+1
	    if (pos <= ubound(warecatab_df,2)-1) then
	      dist = -h - (pos - 1)*drutes_config%fnc_discr_length
	      E = (warecatab_df(layer,pos+1)-warecatab_df(layer,pos))/drutes_config%fnc_discr_length*dist &
	      + warecatab_df(layer,pos)
	    else
              if (present(quadpnt)) E = dual_ret_capf(pde_loc, layer, quadpnt)
              if (present(x)) E = dual_ret_capf(pde_loc, layer, x=x)
	    end if
	  else
	    if (present(quadpnt)) E = dual_ret_capf(pde_loc, layer, quadpnt)
	    if (present(x)) E = dual_ret_capf(pde_loc, layer, x=x)	  
	  end if 
    else
      E = vgfracture(layer)%Ss	
    end if
    
    end function dual_ret_capf_tab 

function dual_coupling_tab(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use dual_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> resulting system elasticity
      real(kind=rkind)::beta,a,gam_par
      real(kind=rkind)::Ks
      real(kind=rkind):: Ka_f,Ka_m,Ka,ex_term
      real(kind=rkind):: hm,hf,one
      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist

      beta=exchange(layer)%beta
      a=exchange(layer)%a
      gam_par=exchange(layer)%gam_par
      Ks=vgexchange(layer)%KS_local(1) 
      
      get_pressh_vals = .true.
      
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
	  hm = pde(1)%getval(quadpnt)
	  hf = pde(2)%getval(quadpnt)
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
    
    if (hm<0) then
	  if ( hm/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	    pos = int(-hm/drutes_config%fnc_discr_length)+1
	    if (pos <= ubound(couptab,2)-1) then
	      dist = -hm - (pos - 1)*drutes_config%fnc_discr_length
	      Ka_m = (couptab(layer,pos+1)-couptab(layer,pos))/drutes_config%fnc_discr_length*dist &
	      + couptab(layer,pos)
	    else
	      Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
	    end if
	  else
	    Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
      end if
    end if
    
    if (hf<0) then
	  if ( hf/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	    pos = int(-hf/drutes_config%fnc_discr_length)+1
	    if (pos <= ubound(couptab,2)-1) then
	      dist = -hf - (pos - 1)*drutes_config%fnc_discr_length
	      Ka_f = (couptab(layer,pos+1)-couptab(layer,pos))/drutes_config%fnc_discr_length*dist&
	       + couptab(layer,pos)
	    else
	      Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
	    end if
	  else
	    Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
      end if
    end if
    
    Ka=0.5*(Ka_f+Ka_m)*Ks
    if(hf /= hm) then
      ex_term=beta/a**2*gam_par*Ka*(hf-hm)
    else
      ex_term=0.0_rkind
    end if
  end function dual_coupling_tab 

function dual_coupling_f_tab(pde_loc, layer, quadpnt, x) result(ex_term)
      use typy
      use dual_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), dimension(:), intent(in), optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> resulting system elasticity
      real(kind=rkind)::beta,a,gam_par
      real(kind=rkind)::Ks
      real(kind=rkind):: Ka_f,Ka_m,Ka,ex_term
      real(kind=rkind):: hm,hf,one
      integer(kind=ikind) :: pos
      real(kind=rkind) :: res, dist

      
      get_pressh_vals = .true.
      
      beta=exchange(layer)%beta
      a=exchange(layer)%a
      gam_par=exchange(layer)%gam_par
      Ks=vgexchange(layer)%KS_local(1) 
      
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
	  hm = pde(1)%getval(quadpnt)
	  hf = pde(2)%getval(quadpnt)
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
    
    if (hm<0) then
	  if ( hm/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	    pos = int(-hm/drutes_config%fnc_discr_length)+1
	    if (pos <= ubound(couptab,2)-1) then
	      dist = -hm - (pos - 1)*drutes_config%fnc_discr_length
	      Ka_m = (couptab(layer,pos+1)-couptab(layer,pos))/drutes_config%fnc_discr_length*dist&
	       + couptab(layer,pos)
	    else
	      Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
	    end if
	  else	
	    Ka_m= dual_coupling_K(pde_loc,layer,x=[hm])
      end if
    end if
    if (hf<0) then
	  if ( hf/drutes_config%fnc_discr_length < 0.1*huge(1)) then
	    pos = int(-hf/drutes_config%fnc_discr_length)+1
	    if (pos <= ubound(couptab,2)-1) then
	      dist = -hf - (pos - 1)*drutes_config%fnc_discr_length
	      Ka_f = (couptab(layer,pos+1)-couptab(layer,pos))/drutes_config%fnc_discr_length*dist &
	      + couptab(layer,pos)
	    else
	      Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
	    end if
	  else
        Ka_f= dual_coupling_K(pde_loc,layer,x=[hf])
      end if
    end if
    
     Ka=0.5*(Ka_f+Ka_m)*Ks
    if(hf /= hm) then
      ex_term=-beta/a**2*gam_par*Ka*(hf-hm)
    else
      ex_term=0.0_rkind
    end if
  end function dual_coupling_f_tab 
! tabular function
 !> creates a table of values of constitutive functions for the Richards equation to be linearly approximated 
 subroutine dual_tabvalues(pde_loc, Kfnc, Cfnc, thetafnc_f, Kfnc_f, Cfnc_f, thetafnc,ex_K_fnc)
      use typy
      use globals
      use pde_objs
      use printtools
      use core_tools

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
	subroutine Kfnc_f(pde_loc, layer, quadpnt, x, tensor, scalar)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in)                           :: layer
	  type(integpnt_str), intent(in), optional                   :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional      :: x
	  real(kind=rkind), intent(out), dimension(:,:), optional   :: tensor
	  real(kind=rkind), intent(out), optional                   :: scalar 
	end subroutine Kfnc_f
      end interface

      interface
	function Cfnc_f(pde_loc, layer, quadpnt,  x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function Cfnc_f
      end interface


      interface
	function thetafnc_f(pde_loc, layer, quadpnt, x) result(val)	
	  use typy
	  use pde_objs
	  use global_objs
	  class(pde_str), intent(in) :: pde_loc	  
	  integer(kind=ikind), intent(in)                      :: layer
	  type(integpnt_str), intent(in), optional             :: quadpnt
	  real(kind=rkind), dimension(:), intent(in), optional :: x
	  real(kind=rkind)                                     :: val
	end function thetafnc_f
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
	
	!if (domainname == "matrix") then
		n = int(maxpress/drutes_config%fnc_discr_length)+1
    !else if (domainname == "fracture") then
     !   n = int(maxpress/drutes_config%fnc_discr_length)
    !end if
      drutes_config%fnc_discr_length = 1.0_rkind*maxpress/n
      dx = drutes_config%fnc_discr_length

      if (.not. allocated(Ktabd_tget_m)) then
        allocate(Ktabd_tget_m(ubound(vgmatrix,1), n))
      end if
      if (.not. allocated(warecatabd_tget_m)) then
        allocate(warecatabd_tget_m(ubound(vgmatrix,1),n))
      end if
      if (.not. allocated(watcontabd_tget_m)) then
        allocate(watcontabd_tget_m(ubound(vgmatrix,1), n))
      end if
      if (.not. allocated(couptab_all)) then
        allocate(couptab_all(ubound(vgmatrix,1), n))
      end if
      if (.not. allocated(Ktabd_tget_f)) then
        allocate(Ktabd_tget_f(ubound(vgmatrix,1), n))
      end if
      if (.not. allocated(warecatabd_tget_f)) then
        allocate(warecatabd_tget_f(ubound(vgmatrix,1),n))
      end if
      if (.not. allocated(watcontabd_tget_f)) then
        allocate(watcontabd_tget_f(ubound(vgmatrix,1), n))
      end if
      
	maxcalls = ubound(vgmatrix,1)*n
	counter = maxcalls
	couptab => couptab_all(:,:)
	!if (domainname == "matrix") then
	  Ktab_dm => Ktabd_tget_m(:,:)
	  warecatab_dm => warecatabd_tget_m(:,:)
	  watcontab_dm => watcontabd_tget_m(:,:)
    !else if (domainname == "fracture") then
	  Ktab_df => Ktabd_tget_f(:,:)
	  warecatab_df => warecatabd_tget_f(:,:)
	  watcontab_df => watcontabd_tget_f(:,:)
    !else
	!  ERROR STOP "runtime error, invalid argument in RE_pointers::domainname not correctly identified"
    !end if

	call write_log(text="creating constitutive function table for: matrix and fracture")
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
	    call Kfnc_f(pde_loc,i, x=(/-(j-1)*dx/), scalar=Ktab_df(i,j))
	    warecatab_df(i,j) = Cfnc_f(pde_loc, i, x=(/-(j-1)*dx/))
	    watcontab_df(i,j) = thetafnc_f(pde_loc, i, x=(/-(j-1)*dx/))
	    couptab(i,j) = ex_K_fnc(pde_loc, i, x=(/-(j-1)*dx/))
	  end do
	end do

    end subroutine dual_tabvalues

end module dual_por