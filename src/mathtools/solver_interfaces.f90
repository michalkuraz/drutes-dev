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

!> \file solver_interfaces.f90
!! \brief interfaces for different solvers, so that can be linked to pde_objs::solve_matrix
!<



module solver_interfaces
  use mtx
  public :: LDU_face
  public :: CG_face
  public :: CG_normal_face
  public :: jacobi_face
  public :: null_problem

  
 contains
 
    subroutine geteigenvals(A, l1, l2)
      use typy
      use mtx
      use matmod
      
      class(matrix), intent(in), target :: A
      real(kind=rkind), intent(out) :: l1, l2
      
      class(matrix), pointer :: asp
      type(symmul) :: asym
      
      asp => A
      call asym%setdata(asp)
      call estimeigvalues(asym,l1,l2)
      
      print *, "max. min. eigenvalues:", l1, l2, "conditioning:", l1/l2
    
    end subroutine geteigenvals
    
    
    subroutine LDU_face(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
                ll1,ll2,cond1,opcnt1,errcode1)
      use mtx
      use typy
      use solvers
      implicit none
      !> matice soustavy\n
      !! musi poskytovat getn, getm, mul (nasobeni vektorem)
      class(matrix), intent(in out) :: A
      !> vektor prave strany
      real(kind=rkind), dimension(:), intent(in) :: b
      !> aproximace reseni, postupne menena
      real(kind=rkind), dimension(:), intent(in out) :: x
      !> maximalni povoleny pocet iteraci, default = n ( Rozmer matice)
      integer(kind=ikind), intent(in), optional :: itmax1
      !> pozadovana relativni zmena rezidua, default = 1e-6
      real(kind=rkind), intent(in), optional :: reps1
      !> informacni podrobnost\n
      !> - 0 ... pracuj tise
      !! - 1 ... minimalni informace
      !! - 10 ... maximalni ukecanost
      integer, intent(in), optional :: ilev1
      !> skutecne provedeny pocet iteraci
      integer(kind=ikind), intent(out), optional :: itfin1
      !> skutecne dosazena relativni zmena residua
      real(kind=rkind), intent(out), optional :: repsfin1
      !> odhad nejvetsiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll1
      !> odhad nejmensiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll2
      !> odhad cisla podminenosti : cond1 = ll1/ll2
      real(kind=rkind), intent(out), optional :: cond1
      !> celkovy pocet provedenych operaci a cas behu
      type(tcount), intent(out), optional :: opcnt1
      !> kod pripadnr chyby
      !! - 0 ... OK
      !! - 1 ... matice neni ctvercova
      !! - 2 ... nesouhlasi b
      !! - 3 ... nesouhasi x
      !! - 4 ... ani jeden z vektoru nesouhlasi
      !! - 5 ... vycerpan povoleny pocet iteraci
      !! - 6 ... prestalo klesat residuum i energie
      integer, intent(out), optional :: errcode1
      
      integer(kind=ikind), dimension(:), allocatable, save :: p1, p2
      
  !! pivtype -- method of pivoting
  !! 0 ... no pivoting (not recommended)
  !! 1 ... full pivoting (use both permutation vector)
  !! 2 ... column pivoting (requires perm1 only)
  !! 3 ... row pivoting (requires perm2 only)
  !! 4 ... diagonal pivoting (symmetric matrix only)
  !! 5 ... diagonal pivoting with minimal degree (symmetric matrix only)
      
      
      
      
      if (.not. allocated(p1)) then
        allocate(p1(ubound(b,1)))
        allocate(p2(ubound(b,1)))
      end if
  
        call LDUd(A, pivtype=0, ilev=0, perm1=p1, perm2=p2)
        
        call LDUback(A, b, x, p1=p1, p2=p2)
        
        if (present(itfin1)) then 
          itfin1 = 1
        end if

    end subroutine LDU_face


    subroutine CG_face(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
                ll1,ll2,cond1,opcnt1,errcode1)
      use mtx
      use typy
      use sparsematrix
      use solvers
      implicit none
      !> matice soustavy\n
      !! musi poskytovat getn, getm, mul (nasobeni vektorem)
      class(matrix), intent(in out) :: A
      !> vektor prave strany
      real(kind=rkind), dimension(:), intent(in) :: b
      !> aproximace reseni, postupne menena
      real(kind=rkind), dimension(:), intent(in out) :: x
      !> maximalni povoleny pocet iteraci, default = n ( Rozmer matice)
      integer(kind=ikind), intent(in), optional :: itmax1
      !> pozadovana relativni zmena rezidua, default = 1e-6
      real(kind=rkind), intent(in), optional :: reps1
      !> informacni podrobnost\n
      !> - 0 ... pracuj tise
      !! - 1 ... minimalni informace
      !! - 10 ... maximalni ukecanost
      integer, intent(in), optional :: ilev1
      !> skutecne provedeny pocet iteraci
      integer(kind=ikind), intent(out), optional :: itfin1
      !> skutecne dosazena relativni zmena residua
      real(kind=rkind), intent(out), optional :: repsfin1
      !> odhad nejvetsiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll1
      !> odhad nejmensiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll2
      !> odhad cisla podminenosti : cond1 = ll1/ll2
      real(kind=rkind), intent(out), optional :: cond1
      !> celkovy pocet provedenych operaci a cas behu
      type(tcount), intent(out), optional :: opcnt1
      !> kod pripadnr chyby
      !! - 0 ... OK
      !! - 1 ... matice neni ctvercova
      !! - 2 ... nesouhlasi b
      !! - 3 ... nesouhasi x
      !! - 4 ... ani jeden z vektoru nesouhlasi
      !! - 5 ... vycerpan povoleny pocet iteraci
      !! - 6 ... prestalo klesat residuum i energie
      integer, intent(out), optional :: errcode1
      integer :: ilevel

      if (.not. present(ilev1) ) then
        ilevel = 0
      else
        ilevel = ilev1
      end if


      call CG(A=A, b=b,x=x,ilev1=ilevel,itmax1=itmax1,reps1=reps1)


    end subroutine CG_face



    subroutine CG_normal_face(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
                ll1,ll2,cond1,opcnt1,errcode1)
      use mtx
      use typy
      use sparsematrix
      use solvers
      implicit none
      !> matice soustavy\n
      !! musi poskytovat getn, getm, mul (nasobeni vektorem)
      class(matrix), intent(in out) :: A
      !> vektor prave strany
      real(kind=rkind), dimension(:), intent(in) :: b
      !> aproximace reseni, postupne menena
      real(kind=rkind), dimension(:), intent(in out) :: x
      !> maximalni povoleny pocet iteraci, default = n ( Rozmer matice)
      integer(kind=ikind), intent(in), optional :: itmax1
      !> pozadovana relativni zmena rezidua, default = 1e-6
      real(kind=rkind), intent(in), optional :: reps1
      !> informacni podrobnost\n
      !> - 0 ... pracuj tise
      !! - 1 ... minimalni informace
      !! - 10 ... maximalni ukecanost
      integer, intent(in), optional :: ilev1
      !> skutecne provedeny pocet iteraci
      integer(kind=ikind), intent(out), optional :: itfin1
      !> skutecne dosazena relativni zmena residua
      real(kind=rkind), intent(out), optional :: repsfin1
      !> odhad nejvetsiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll1
      !> odhad nejmensiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll2
      !> odhad cisla podminenosti : cond1 = ll1/ll2
      real(kind=rkind), intent(out), optional :: cond1
      !> celkovy pocet provedenych operaci a cas behu
      type(tcount), intent(out), optional :: opcnt1
      !> kod pripadnr chyby
      !! - 0 ... OK
      !! - 1 ... matice neni ctvercova
      !! - 2 ... nesouhlasi b
      !! - 3 ... nesouhasi x
      !! - 4 ... ani jeden z vektoru nesouhlasi
      !! - 5 ... vycerpan povoleny pocet iteraci
      !! - 6 ... prestalo klesat residuum i energie
      integer, intent(out), optional :: errcode1
      integer :: ilevel

      if (.not. present(ilev1) ) then
        ilevel = 0
      else
        ilevel = ilev1
      end if


      call CGnormal(A=A, b=b,x=x,ilev1=ilevel,itmax1=itmax1,reps1=reps1, itfin1=itfin1, repsfin1=repsfin1)


    end subroutine CG_normal_face

    
    subroutine jacobi_face(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
                ll1,ll2,cond1,opcnt1,errcode1)
      use mtx
      use typy
      use solvers
      implicit none
      !> matice soustavy\n
      !! musi poskytovat getn, getm, mul (nasobeni vektorem)
      class(matrix), intent(in out) :: A
      !> vektor prave strany
      real(kind=rkind), dimension(:), intent(in) :: b
      !> aproximace reseni, postupne menena
      real(kind=rkind), dimension(:), intent(in out) :: x
      !> maximalni povoleny pocet iteraci, default = n ( Rozmer matice)
      integer(kind=ikind), intent(in), optional :: itmax1
      !> pozadovana relativni zmena rezidua, default = 1e-6
      real(kind=rkind), intent(in), optional :: reps1
      !> informacni podrobnost\n
      !> - 0 ... pracuj tise
      !! - 1 ... minimalni informace
      !! - 10 ... maximalni ukecanost
      integer, intent(in), optional :: ilev1
      !> skutecne provedeny pocet iteraci
      integer(kind=ikind), intent(out), optional :: itfin1
      !> skutecne dosazena relativni zmena residua
      real(kind=rkind), intent(out), optional :: repsfin1
      !> odhad nejvetsiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll1
      !> odhad nejmensiho vlastniho cisla
      real(kind=rkind), intent(out), optional :: ll2
      !> odhad cisla podminenosti : cond1 = ll1/ll2
      real(kind=rkind), intent(out), optional :: cond1
      !> celkovy pocet provedenych operaci a cas behu
      type(tcount), intent(out), optional :: opcnt1
      !> kod pripadnr chyby
      !! - 0 ... OK
      !! - 1 ... matice neni ctvercova
      !! - 2 ... nesouhlasi b
      !! - 3 ... nesouhasi x
      !! - 4 ... ani jeden z vektoru nesouhlasi
      !! - 5 ... vycerpan povoleny pocet iteraci
      !! - 6 ... prestalo klesat residuum i energie
      integer, intent(out), optional :: errcode1
      
      call jacobi(a,b,x,itmax1, reps1)

    end subroutine jacobi_face
    
    
  subroutine null_problem(A,b)
    use typy
    use globals
    use sparsematrix
    class(extsmtx), intent(in out) :: A
    real(kind=rkind), dimension(:), intent(in out), optional :: b
    integer(kind=ikind) :: n_rows, m_rows

    
    n_rows = A%getn()
    m_rows = A%getm()
  
    call A%init(n_rows,m_rows)
    
    if (present(b)) then
      b = 0.0_rkind
    end if
    
    call A%rowsfilled%clear
    
  end subroutine null_problem

end module solver_interfaces
