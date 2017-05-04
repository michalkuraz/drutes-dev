module dual_por
! This module contains governing functions of the dual permeability model
  use typy
  use global_objs
  use dual_globals


  public:: dual_mualemm,dual_mualemf, dual_ret_capf, dual_ret_capm
  public:: vangen_d_f, vangen_d_m
  public:: dual_inicond_f,dual_inicond_m
  public:: darcy_law_d 
  public :: getval_retot_dual

  contains 

    !> specific function for Richards equation in H-form (total hydraulic head form), replaces pde_objs::getvalp1 in order to distinguish between H and h 
    !> pointers point to this with pde_loc%getval
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
      

	  
      if (quadpnt%preproc) then
      
	D = drutes_config%dimen
      
	call getcoor(quadpnt, xyz(1:D))

	val = getvalp1(pde_loc, quadpnt) - xyz(D)
    ! 	
      else
	val = getvalp1(pde_loc, quadpnt)
      end if

      
	
    end function getval_retot_dual

   subroutine dual_inicond_m(pde_loc) 
      	use typy
       	use globals
       	use global_objs
       	use pde_objs
       	use dual_globals
       	use debug_tools
       	use geom_tools
       	!use RE_constitutive 
        
       class(pde_str), intent(in out) :: pde_loc
        integer(kind=ikind) :: i, j, k,l, m, layer, D,n
        real(kind=rkind) :: value
        
        D = drutes_config%dimen
  
       select case (vgmatrix(1_ikind)%icondtype)
            case("input")
		      call map1d2dJ(pde_loc,"drutes.conf/REdual/hinim.in")
        end select
        do i=1, elements%kolik
          layer = elements%material(i,1)

          do j=1, ubound(elements%data,2)
            k = elements%data(i,j)
            l = nodes%edge(k)
            m = pde_loc%permut(k)
            if (m == 0) then
	      	  call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
              pde_loc%solution(k) =  value 
            else
        	  select case (vgmatrix(layer)%icondtype)
				case("H_tot")
		  		  pde_loc%solution(k) = vgmatrix(layer)%initcond+disttozero ! nodes%data(k,D)
				case("hpres")
		  		  pde_loc%solution(k) = vgmatrix(layer)%initcond+nodes%data(k,D)
	      	  end select
            end if
          end do   
        end do
        

      end subroutine dual_inicond_m
      
    subroutine dual_inicond_f(pde_loc) 
      	use typy
       	use globals
       	use global_objs
       	use pde_objs
       	use dual_globals
        use RE_constitutive 
        use geom_tools
		class(pde_str), intent(in out) :: pde_loc
		
        integer(kind=ikind) :: i, j, k,l, m, layer, D,n
        real(kind=rkind) :: value
        
        D = drutes_config%dimen
          select case (vgfracture(1_ikind)%icondtype)
            case("input")
		     call map1d2dJ(pde_loc,"drutes.conf/REdual/hinif.in")
          end select

        do i=1, elements%kolik
          layer = elements%material(i,1)

          do j=1, ubound(elements%data,2)
            k = elements%data(i,j)
            l = nodes%edge(k)
            m = pde_loc%permut(k)
            if (m == 0) then
	      call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
              pde_loc%solution(k) =  value 
            else
        	  select case (vgfracture(layer)%icondtype)
				case("H_tot")
		  		  pde_loc%solution(k) = vgfracture(layer)%initcond+disttozero !- nodes%data(k,D)
				case("hpres")
		  		  pde_loc%solution(k) = vgfracture(layer)%initcond+nodes%data(k,D)
	      	  end select
            end if
          end do   
        end do

      end subroutine dual_inicond_f

    ! subroutine dual_inicond_m(pde_loc)
!       use typy
!       use globals
!       use global_objs
!       use pde_objs
!       use dual_globals
! 	class(pde_str), intent(in out) :: pde_loc
! 	integer(kind=ikind) :: i, j, k,l, m, layer
! 	real(kind=rkind) :: value
! 	
! 	  do i=1, elements%kolik
! 		layer = elements%material(i,1)
! 		do j=1, ubound(elements%data,2)
! 		  k = elements%data(i,j)
! 	      pde_loc%solution(k)=vgmatrix(layer)%initcond
! 		end do   
! 	  end do
! 	  
!       end subroutine dual_inicond_m
!     
!     subroutine dual_inicond_f(pde_loc)
!       use typy
!       use globals
!       use global_objs
!       use pde_objs
!       use dual_globals
!       class(pde_str), intent(in out) :: pde_loc
!       integer(kind=ikind) :: i, j, k,l, m, layer
!       real(kind=rkind) :: value
! 	
!       do i=1, elements%kolik
! 	layer = elements%material(i,1)
! 	do j=1, ubound(elements%data,2)
! 	  k = elements%data(i,j)
! 	  pde_loc%solution(k)=vgfracture(layer)%initcond
! 	end do   
!       end do
! 	  
!     end subroutine dual_inicond_f
! !  
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
      type(integpnt_str) :: quadpnt_loc
	

	  
	if (present(quadpnt)) then
	  quadpnt_loc=quadpnt
	  quadpnt_loc%preproc=.true.
	 h = pde_loc%getval(quadpnt_loc)
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
	weight=exchange(layer)%weightm


	one=1.0_rkind    
	if(h < 0.0_rkind) then
	  Kr=(one-(alpha*abs(h))**(n*m)*(one+(alpha*abs(h))**n)**(-m))**2/(one+(alpha*abs(h))**n)**(m/2)
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
	type(integpnt_str) :: quadpnt_loc
	
	
	    
	  if (present(quadpnt)) then
	    quadpnt_loc=quadpnt
	    quadpnt_loc%preproc=.true.
	   h = pde_loc%getval(quadpnt_loc)
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
	    Kr=(one-(alpha*abs(h))**(n*m)*(one+(alpha*abs(h))**n)**(-m))**2/(one+(alpha*abs(h))**n)**(m/2)
	  else
	    Kr=1.0_rkind
	  end if
	  
	  if (present(tensor)) then
		      tensor=vgfracture(layer)%KS*Kr*(weight)
	  end if
	    
	  if (present(scalar)) then
	    scalar=Kr*(weight)
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
      type(integpnt_str) :: quadpnt_loc
      
      
      
      if (present(quadpnt)) then
	quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
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
     weight=exchange(layer)%weightm
     

     
     if(h<0.0_rkind) then
       C=(thetaS-thetaR)*alpha*m*n*(alpha*abs(h))**(n-1)*((alpha*abs(h))**n+1)**(-m-1)
     else
       E=vgmatrix(layer)%Ss*weight
       RETURN
     end if
     
      E=C*weight
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
      type(integpnt_str) :: quadpnt_loc
      
      
      
      
      if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
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
       C=(thetaS-thetaR)*alpha*m*n*(alpha*abs(h))**(n-1)*((alpha*abs(h))**n+1)**(-m-1)
     else
       E=vgfracture(layer)%SS*(weight)
       RETURN
     end if     
      E=C*weight
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
   type(integpnt_str) :: quadpnt_loc
      
      
     if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
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
   type(integpnt_str) :: quadpnt_loc
   
      
     if (present(quadpnt)) then
        quadpnt_loc=quadpnt
	quadpnt_loc%preproc=.true.
	h = pde_loc%getval(quadpnt_loc)
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
        print *, "exited from Re_dual_totH::darcy_law"
	ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from Re_dual_totH::darcy_law"
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
	  print *, "exited from Re_dual_totH::darcy_law"
	  ERROR STOP
	end if
	h = x(1)
	gradH(1:D) = grad
      end if
      
      call pde_loc%pde_fnc(pde_loc%order)%dispersion(pde_loc, layer, x=(/h/), tensor=K(1:D, 1:D))
     
      
      vct(1:D) = matmul(-K(1:D,1:D), gradH(1:D))

	 !print*,vct
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