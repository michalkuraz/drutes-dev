module dual_por

public:: dual_mualemm,dual_mualemf, dual_ret_capf, dual_ret_capm, dual_coupling
public:: vangen_d_f, vangen_d_m, dual_coupling_f
public:: dual_inicond_f,dual_inicond_m
public :: darcy_law_d


contains 

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
      
    if (present(quadpnt)) then
	  h = pde_loc%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
	    print *, "ERROR: van Genuchten function is a function of a single variable h"
	    print *, "       your input data has:", ubound(x,1), "variables"
	    print *, "exited from re_constitutive::mualem"
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
      
    if (present(quadpnt)) then
	  h = pde_loc%getval(quadpnt)
    else
      if (ubound(x,1) /=1) then
	    print *, "ERROR: van Genuchten function is a function of a single variable h"
	    print *, "       your input data has:", ubound(x,1), "variables"
	    print *, "exited from re_constitutive::mualem"
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
      
      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> vg parameters, later from conf file
      real(kind=rkind)::n,m,alpha,thetaS,thetaR,weight
      real(kind=rkind)::E,h,C
      
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
	      print *, "exited from re_constitutive::vangen_elast"
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
	      print *, "exited from re_constitutive::vangen_elast"
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

      
     if (present(quadpnt)) then
	   h = pde_loc%getval(quadpnt)
     else
	 if (ubound(x,1) /=1) then
	   print *, "ERROR: van Genuchten function is a function of a single variable h"
	   print *, "       your input data has:", ubound(x,1), "variables"
	   print *, "exited from re_constitutive::vangen"
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

      
     if (present(quadpnt)) then
	   h = pde_loc%getval(quadpnt)
     else
	 if (ubound(x,1) /=1) then
	   print *, "ERROR: van Genuchten function is a function of a single variable h"
	   print *, "       your input data has:", ubound(x,1), "variables"
	   print *, "exited from re_constitutive::vangen"
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
      
            
      if (present(quadpnt)) then
	    hm = pde(1)%getval(quadpnt)
	    hf = pde(2)%getval(quadpnt)
      else
	    if (ubound(x,1) /=1) then
	      print *, "ERROR: van Genuchten function is a function of a single variable h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      ERROR STOP
	    end if
      	if (ubound(x,1) /=1) then
	      print *, "ERROR: van Genuchten function is a function of a single variable h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      print *, "exited from re_constitutive::vangen_elast"
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
      
            
      if (present(quadpnt)) then
	    hm = pde(1)%getval(quadpnt)
	    hf = pde(2)%getval(quadpnt)
      else
	    if (ubound(x,1) /=1) then
	      print *, "ERROR: van Genuchten function is a function of a single variable h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      ERROR STOP
	    end if
      	if (ubound(x,1) /=1) then
	      print *, "ERROR: van Genuchten function is a function of a single variable h"
	      print *, "       your input data has:", ubound(x,1), "variables"
	      print *, "exited from re_constitutive::vangen_elast"
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
	h = pde_loc%getval(quadpnt)
	call pde_loc%getgrad(quadpnt, gradient)
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
      
end module dual_por