module re_analytical

  public :: tracy_bc, tracy

  contains

    !> Boundary condition for Tracy's analytical solution
    subroutine tracy_bc(bcval, x_coord, width, hinit)
      use typy
      use re_globals

      !> output boundary value
      real(kind=rkind), intent(out) :: bcval
      !> with - the width of the domain, coord - the x coordinate at the boundary, hihit - initial condition (should be very dry))
      real(kind=rkind), intent(in) :: width, x_coord, hinit

      real(kind=rkind) :: hbar, alpha
    
      alpha = vgset(1)%alpha
    
      hbar = 1 - exp(alpha*hinit)

      bcval = 1.0/alpha*log(exp(alpha*hinit) + hbar*sin(4*atan(1.0)*x_coord/width))


    end subroutine tracy_bc
  
   !> analytical solution to the transient 2D Richards' equation based on (Tracy, 2006)
    subroutine tracy(hinit, coord, t, width, length, h)
      use typy
      use re_constitutive

      !> initial state (must be constant)
      real(kind=rkind), intent(in) :: hinit
      !> point coordinates
      real(kind=rkind), dimension(:), intent(in) :: coord
      !> simulation time
      real(kind=rkind), intent(in)               :: t, width, length
      !> solution
      real(kind=rkind), intent(out)              :: h



!--------local variables--------------------
      real(kind=rkind) :: lambda
      real(kind=rkind) :: c
      real(kind=rkind) :: gamma
      real(kind=rkind) :: phi
      real(kind=rkind) :: hbar
      real(kind=rkind) :: beta
      real(kind=rkind) :: ho
      real(kind=rkind) :: hr
      real(kind=rkind) :: suma
      real(kind=rkind) :: hss
      real(kind=rkind) :: tmp
      real(kind=rkind) :: alpha
      real(kind=rkind) :: a
      real(kind=rkind) :: L, absval
      integer(kind=ikind) :: i
      
      if (abs(t) < epsilon(t)) then
        if (abs(coord(2)-length) < epsilon(length)) then
          call tracy_bc(h, coord(1), width, hinit)
        else
          h=hinit
        end if
        RETURN
      end if
  

      a = width

      L = length


      alpha = vgset(1)%alpha

      ho = 1-exp(alpha*hinit)

      beta = sqrt(alpha**2/4 + (4*atan(1.0)/a)**2)

      c = alpha*(vgset(1)%ths-vgset(1)%thr)/vgset(1)%Ks(1,1)

      hss = ho*sin(4*atan(1.0)*coord(1)/a)*exp(alpha/2*(L-coord(2)))*sinh(beta*coord(2))/sinh(beta*L)

      suma = 0

      i = 0

      do
        i = i+1
        tmp = suma
        lambda = i*4*atan(1.0)/L
        gamma = 1/c*(beta*beta + lambda*lambda)
        tmp = ((-1)**i)*lambda/gamma*sin(lambda*coord(2))*exp(-gamma*t)
        if (i==1) absval=abs(tmp)
	suma = suma + tmp
        if (abs(tmp) < absval*epsilon(tmp) ) then
          EXIT
        end if
      
      end do

      phi = 2*ho/(L*c)*sin(4*atan(1.0)*coord(1)/a)*exp(alpha/2*(L-coord(2)))*suma


      hbar = phi + hss


      h = 1/alpha*log(exp(alpha*hinit)+hbar)


  end subroutine tracy

end module re_analytical
