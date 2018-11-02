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

!> \file re_analytical.f90
!! \brief Tracy's analytical solution for the Richards equation.
!<

!>
!! Tracy's analytical solution for the Richards equation.
!! The computational domain $\Omega$ is a rectangle  $ a [L]  \times L [L] $. !!The initial condition is a so-called "very dry" medium defined by pressure head \f$h_r\f$ [L].
!! Boundary conditions for the pressure head are Dirichlet 
!!\f$h(x,t) = h_r \quad \forall (x,t) \in \Gamma_{1,2,4} \times [0,T)\f$ on the bottom edge \f$\Gamma_1\f$
!!and vertical edges \f$\Gamma_2\f$ and \f$\Gamma_4\f$, and on the 
!!top edge \f$\Gamma_3\f$ the value of \f$h\f$ is given by the relation
!!
!!\f[ h(x,t) = \frac{1}{\alpha}\log \left(e^{\alpha h_r} + \left( 1- e^{\alpha h_r} \right)\sin \frac{\pi x_1}{a} \right)  \quad \forall (x,t) \in \Gamma_{3} \times [0,T) \f]
!!
!!The analytical solution formulas state as
!!
!!
!!
!!\f[ 
!! h(x_1,x_3) =  \frac{1}{\alpha} \log \left( \textrm{e}^{\alpha h_r} + \frac{2h_0}{Lc} \textrm{sin} \left( \frac{\pi x_1}{a} \right) \textrm{e}^{\frac{\alpha}{2}(L-x_3)} \sum_{i=1}^{\infty}  (-1)^i \frac{\lambda_i}{\gamma} \textrm{sin}(\lambda_ix_3) \textrm{e}^{-\gamma t} +   \right.
!! \qquad \left.\vphantom{\int_t} + \bar{h_0} \textrm{sin} \left( \frac{\pi x_1}{a} \right) \textrm{e}^{\frac{\alpha}{2}(L-x_3)}\frac{\textrm{sinh}(\beta x_3)}{\textrm{sinh}(\beta L)}  \right)
!!\f]
!!
!!where
!!
!! \f[
!!\beta =  \sqrt{ \frac{\alpha^2}{4}  + \left(\frac{\pi}{a}\right)^2 }, \hspace{4mm} c = !!\frac{\alpha(\theta_s - \theta_r}{K_s} ,
!!\f]

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
