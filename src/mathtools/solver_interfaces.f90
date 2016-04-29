module solver_interfaces
  public :: LDU_face
  public :: CG_face
  public :: CG_normal_face
  public :: jacobi_face
  public :: null_problem

    contains

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
        
        
         call LDUd(A, ilev=0)
         
         call LDUback(A, b, x)

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


        call CGnormal(A=A, b=b,x=x,ilev1=ilevel,itmax1=itmax1,reps1=reps1)


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