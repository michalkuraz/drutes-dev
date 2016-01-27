! Copyright 2008 Michal Kuraz, Petr Mayer


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

module band
  use typy
  implicit none
     !> pocet volani resice lin. rovnic
  ! pasova matice je ulozena po radcich, jen pas, deklarace je
  ! real(kind=kind_real) :: A(1:n,-p:q)
  ! n ... rad matice
  ! p ... pocet poddiagonal - sirka dolniho pasu
  ! q ... pocet naddiagonal - sirka horniho pasu
  ! k-ta digonala == A(:,k)
  
  public :: Axpy, Set_Laplace, axpyna, LDU, lusolve
  contains



  subroutine Axpy(A,x,y,n,p,q)
   use typy
    integer(kind=ikind), intent(in) :: n,p,q
    real(kind=rkind), dimension(1:,-p:), intent(in) :: A
    real(kind=rkind), dimension(1:), intent(in)     :: x
    real(kind=rkind), dimension(1:),intent(in out)  :: y
    ! y := Ax + y
    ! lokalni promenne
    integer(kind=ikind) :: i,l1,l2,k1,k2
    do i = -p,q 
       l1 = max(1_ikind,1-i)    ! dolni mez pro y a A
       l2 = min(n,n-i)    ! horni mez pro y a A
       k1 = max(1_ikind,1+i)    ! dolni mez pro x
       k2 = min(n,n+i)    ! horni mez pro x
       y(l1:l2) = y(l1:l2) + A(l1:l2,i)*x(k1:k2)
    end do
  end subroutine Axpy


  subroutine Set_Laplace(A,n,p,q)
   use typy
    integer(kind=ikind), intent(in)  :: p,q
    real(kind=rkind), dimension(1:,-p:), intent(in out) :: A
    integer(kind=ikind), intent(out) :: n
    integer(kind=ikind)              :: i
    n = p*q 
    A = 0  ! tohle ma tendenci generovat obri matici nul
    a(:,0) = 4
    a(:,(/ -p,-1_ikind,1_ikind,p /)) = -1
    do i=1,q-1
       a(i*p+1,-1) = 0
       a(i*p,1) = 0
    end do
  end subroutine Set_Laplace


  subroutine axpyna(a,n,p,q,x,y)   
   use typy
    ! y = a krat x plus y 
    integer(kind=ikind), intent(in) :: n,p,q 
    real(kind=rkind), dimension(1:,-p:), intent(in) :: a
    real(kind=rkind), dimension(1:), intent(in) :: x
    real(kind=rkind), dimension(1:), intent(inout) :: y
    
    integer(kind=ikind) :: i,l1,l2,k1,k2

    do i = -p,q
       l1 = max(1_ikind,1-i)    ! dolni mez pro y a A
       l2 = min(n,n-i)    ! horni mez pro y a A
       k1 = max(1_ikind,1+i)    ! dolni mez pro x
       k2 = min(n,n+i)    ! horni mez pro x
       y(l1:l2) = y(l1:l2) + a(l1:l2,i)*x(k1:k2)
    end do
    
  endsubroutine axpyna
  
  subroutine LDU(a,n,p,q,err)
   use typy
   use globals
    integer(kind=ikind), intent(in) :: n,p,q 
    real(kind=rkind), dimension(1:,-p:), intent(inout) :: a
    integer(kind=ikind), intent(out) :: err  
    
    integer(kind=ikind) :: i,j,k  
    real(kind=rkind), dimension(1:p+1,1:q+1) :: b
    err = 0
    solver_call = solver_call+1
    do i =1,n-1
       ! Sup s ni do separace
       b = 0
       do j = 0,min(p,n-i)
          b(j+1,:) = a(i+j,-j:q-j)
       enddo
       b(2:p+1,1) = b(2:p+1,1) / b(1,1)
       do j = 2,p+1
          do k = 2,q+1
             b(j,k) = b(j,k) - b(j,1)*b(1,k)
          enddo
       enddo
       b(1,2:q+1) = b(1,2:q+1) / b(1,1)
       ! Pujdes na svobodu
       do j = 0,min(p,n-i)
          a(i+j,-j:q-j) = b(j+1,:)
       enddo
    enddo
  end subroutine LDU
  
  
  subroutine lusolve (a,n,p,q,err,b,x)
   use typy
   use linAlg
   use globals
    integer(kind=ikind), intent(in) :: n,p,q
    real(kind=rkind), dimension(1:,-p:), intent(inout) :: a
    integer(kind=ikind), intent(out) :: err
    real(kind=rkind), dimension(1:), intent(in) :: b
    real(kind=rkind), dimension(1:), intent(out) :: x 
    integer(kind=ikind) :: i 
    err = 0
    x=b
    solver_call = solver_call+1
    do i=2,n 
       x(i)=x(i)-dot_product(x(max (i-p,1_ikind):i-1),a(i,max(1-i,-p):-1))
    enddo
    x=x/a(:,0)
    
    do i=n-1,1,-1
       x(i)=x(i)-dot_product(x(i+1:min (i+q,n)),a(i,1:min (n-i,q)))
    enddo     
  end subroutine lusolve
  


end module band

