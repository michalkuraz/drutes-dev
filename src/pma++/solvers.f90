!> \file solvers.f90
!! \brief ruzne resice soustav linearnich rovnic
!!
!<

!> \page resice Resice soustav linearnich rovnic
!! --
!! Metody jsou ve dvou zakladnich tridach
!! - \subpage iterace "Iteracni metody"
!! - \subpage elim  "Prime metody"
!<

!> \page iterace "Iteracni metody"
!! - \subpage Jacobi "Jacobiova meoda"
!! - \subpage gs "Gauss-Seidelova metoda"
!! - \subpage sd "Metoda nejvetsiho spadu"
!! - \subpage cg "Metoda sdruzenych gradientu"
!! - \subpage cgn "Metoda sdruzenych gradientu pro normalni rovnice"
!! - \subpage pcg "Metoda predpodminenych sdruzenych gradientu"
!!
!<


!> \page elim "Prime metody"
!! \subpage ldu "LDU"
!<

!> resice soustav
module solvers
    implicit none
    public :: LDU
    public :: LDUd
    public :: Split
    public :: LDUback
    public :: jacobi
    public :: GS
    public :: SD
    public :: SDnormal
    public :: CG
    public :: CGnormal
    public :: PCG
    public :: PCGnormal
    private :: vycdet
    private :: esteig
    contains
    !> \page Jacobi Jacobiova metoda
    !!
    !! Jde o iteracni postup
    !!  \f$  x_{i+1} = x_i + M^{-1}*(b-A*x_i) \f$
    !! kde \f$ M \f$ je diagonala matice \f$ A \f$ . \n
    !! Oznacime-li residuum jako \f$ r_i = b-A*x_i \f$, pak \f$
    !! reps = || r_i || / || r_0 || \f$ .
    !<
    
    !> Jacobiova metoda
    subroutine jacobi(a,b,x,maxit1,reps1, it, repsfin, l1, l2, countop,&
        ilev1)
        use typy
        use mtx
        implicit none
        !> matice soustavy
        class(matrix), intent(in) :: a
        !> prava strana
        real(kind=rkind), dimension(:), intent(in) :: b
        !> pocitane reseni soustavy
        real(kind=rkind), dimension(:), intent(inout) :: x
        !> maximalni povoleny pocet iteraci
        integer(kind=ikind), intent(in),optional :: maxit1
        !> pozadovane relativni residuum
        real(kind=rkind), intent(in), optional :: reps1
        !> skutecny pocet iteraci
        integer(kind=ikind), intent(out), optional :: it
        !> skutecne relativni residuum
        real(kind=rkind), intent(out), optional :: repsfin
        !> odhad rychlosti konvergence pomoci mocninne metody
        real(kind=rkind), intent(out), optional :: l1
        !> odhad rychlosti konvergence pomoci odmocniny
        real(kind=rkind), intent(out), optional :: l2
        !> pocitadlo operaci
        type(tcount), intent(inout), optional :: countop
        !> podrobnost informaci
        integer, intent(in), optional :: ilev1
        
        integer(kind=ikind) :: maxit,it1, ilev
        real(kind=rkind) :: reps
        real(kind=rkind) :: r0,r1,r2, repsfin1
        real(kind=rkind) :: ll1,ll2
        real(kind=rkind), dimension(:), allocatable :: r,diag
        integer(kind=ikind) :: i,n
        type(tcount) :: countop1
        real(kind=rkind) :: t1, t2
        
        !print *,"jacobi zacina"
        call cpu_time(t1)
        !! vyresime nepovinne parametry
        n = a%getn()
        if (present(maxit1)) then
            maxit = maxit1
            else
            maxit = n
        end if
        if (present(reps1)) then
            reps = reps1
            else
            reps = 1.0e-6
        end if
        if (present(ilev1)) then
            ilev = ilev1
            else
            ilev = 0
        end if
        
        !print *, "nepovinne vyseny"
        allocate(r(1:a%getn()))
        allocate(diag(1:a%getn()))
        r=b-a%mul(x,countop1)
        r0 = norm2(r)
        countop1%ad = countop1%ad + a%getn()
        do i=1,a%getn()
            diag(i) = a%get(i,i)
        end do
        it1 = 0
        !print *,"pripraveno"
        do
            if ( it1 > maxit ) exit
            r1 = norm2(r)
            repsfin1 = r1/r0
            if (r1 /= 0) then
                ll1 = r1/r2
                ll2 = repsfin1**(1.0_rkind/it1)
            end if
            if ( ilev > 0 )  print"(i16,4ES18.10)", it1, r0,r1,ll1,ll2
            if ( repsfin1 < reps) exit
            do i=1, a%getn()
                x(i) = x(i) + r(i)/diag(i)
            end do
            r=b-a%mul(x,countop1)
            countop1%ad  = countop1%ad  + 2*a%getn()
            countop1%div = countop1%div +   a%getn()
            it1 = it1 + 1
            r2  = r1
        end do
        deallocate(r,diag)
        call cpu_time(t2)
        !print *,"jacobi konci"
        ! vyresime navratove hodnoty
        if (present(countop)) then
            countop%ad  = countop%ad  + countop1%ad
            countop%mul = countop%mul + countop1%mul
            countop%div = countop%div + countop1%div
            countop%time = t2-t1
        end if
        if (present(it)) it = it1
        if (present(repsfin)) repsfin = repsfin1
        if (present(l1)) l1 = ll1
        if (present(l2)) l2 = ll2
    end subroutine jacobi
    
    !> \page gs Gauss-Seidelova metoda
    !!
    !! Jde o iteracni postup
    !!  \f$  x_{i+1} = x_i + M^{-1}*(b-A*x_i) \f$
    !! kde \f$ M \f$ je dolni trojuhelnik matice \f$ A \f$ . \n
    !! Oznacime-li residuum jako \f$ r_i = b-A*x_i \f$, pak \f$
    !! reps = || r_i || / || r_0 || \f$ .
    !<
    
    !> Gauss-Seidelova metoda
    subroutine GS(a,b,x,maxit,reps, it, repsfin, l1, l2, countop, ilev)
        use typy
        use mtx
        implicit none
        !> matice soustavy
        class(matrix), intent(in) :: a
        !> prava strana
        real(kind=rkind), dimension(:), intent(in) :: b
        !> pocitane reseni soustavy
        real(kind=rkind), dimension(:), intent(inout) :: x
        !> maximalni povoleny pocet iteraci
        integer(kind=ikind), intent(in) :: maxit
        !> pozadovane relativni residuum
        real(kind=rkind), intent(in) :: reps
        !> skutecny pocet iteraci
        integer(kind=ikind), intent(out) :: it
        !> skutecne relativni residuum
        real(kind=rkind), intent(out) :: repsfin
        !> odhad rychlosti konvergence pomoci mocninne metody
        real(kind=rkind), intent(out) :: l1
        !> odhad rychlosti konvergence pomoci odmocniny
        real(kind=rkind), intent(out) :: l2
        !> pocitadlo operaci
        type(tcount), intent(inout), optional :: countop
        !> podrobnost informaci
        integer, intent(in), optional :: ilev
        
        real(kind=rkind) :: r0,r1,r2
        real(kind=rkind), dimension(:), allocatable :: r,diag
        integer(kind=ikind) :: i,j,k
        type(tcount) :: countop1
        real(kind=rkind) :: t1, t2
        
        call cpu_time(t1)
        allocate(r(1:a%getn()))
        allocate(diag(1:a%getn()))
        r=b-a%mul(x,countop1)
        r0 = norm2(r)
        countop1%ad = countop1%ad + a%getn()
        do i=1,a%getn()
            diag(i) = a%get(i,i)
        end do
        it = 0
        do
            !!!!!!!! predelet - tohle je Jacobiovka
            if ( it > maxit ) exit
            r1 = norm2(r)
            repsfin = r1/r0
            if (r1 /= 0) then
                l1 = r1/r2
                l2 = repsfin**(1.0_rkind/it)
            end if
            if ( ilev > 0 )  print"(i16,4ES18.10)", it, r0,r1,l1,l2
            if ( repsfin < reps) exit
            do i=1, a%getn()
                x(i) = x(i) + r(i)/diag(i)
            end do
            r=b-a%mul(x,countop1)
            countop1%ad  = countop1%ad  + 2*a%getn()
            countop1%div = countop1%div +   a%getn()
            it = it + 1
            r2 = r1
        end do
        deallocate(r,diag)
        call cpu_time(t2)
        if (present(countop)) then
            countop%ad  = countop%ad  + countop1%ad
            countop%mul = countop%mul + countop1%mul
            countop%div = countop%div + countop1%div
            countop%time = t2-t1
        end if
    end subroutine GS
    
    
    
    
    !> \page ldu "LDU rozklad"
    !<
    
    !> nedestruktivni LDU rozklad
    !!
    !! pro volbu pivota plati nasledujici pravidla\n
    !! 0 ... zadna\n
    !! 1 ... uplny vyber hlavniho prvku (nutne obe permutace )\n
    !! 2 ... sloupcovy vyber hlavniho prvku(nutno perm1)\n
    !! 3 ... radkovy vyber hlavniho prvku(nutno perm2)\n
    !! 4 ... diagonalni vyber hlavniho prvku(nutno aspon jedna = obe shodne)\n
    !! 5 ... diagonalni vyber hlavniho prvku pro minimal degree\n
    !!              (nutno aspon jedna = obe shodne) \n
    !! 6 ... rid se predepsanou permutaci\n
    !! Pro situace 1 az 6 musi byt pritomny permutacni vektory\n
    !! Pokud nejsou pozadavky splneny - zvoli se nejpodobnejsi volba\n
    !!
    !! pro podrobnost informaci plati\n
    !! 0 ... tichy chod\n
    !! 1 ... super podrobne\n
    !<
    subroutine LDU(A,LU,info,pivtype,ilev,perm1,perm2)
        use mtx
        use typy
        implicit none
        !> rozkladana matice
        class(matrix), intent(in) :: A
        !> rozlozena matice
        class(matrix), intent(inout) :: LU
        !> pocty operaci
        type(tcount), intent(inout), optional :: info
        !> typ pivotaze
        integer, intent(in), optional :: pivtype
        !> podrobnost informaci
        integer, intent(in), optional :: ilev
        !> permutace pro i
        integer(kind=ikind), dimension(:), allocatable, &
        intent(inout), optional :: perm1
        !> permutace pro j
        integer(kind=ikind), dimension(:), allocatable,&
        intent(inout), optional :: perm2
        
        integer :: pt1,il1
        integer(kind=ikind),dimension(:),allocatable :: p1,p2
        type(tcount) :: inf1
        
        if (present(pivtype)) then
            pt1 = pivtype
            else
            pt1 = 0
        end if
        if (present(ilev)) then
            il1 = ilev
            else
            il1 = 0
        end if
        allocate(p1(1:A%getn()))
        allocate(p2(1:A%getm()))
        select case(pt1)
            case (0)
            print *, "nula"
            case (1)
            case (2)
            case (3)
            case (4)
            case (5)
            case (6)
            if (present(perm1)) then
                else
            end if
            case default
            print *,"chybna volba"
        end select
        
        call LU%copy(A)
        print *,"zkopirovano"
        call LDUd(LU,inf1,pt1,il1,p1,p2)
        if ( present(perm1) ) perm1 = p1
        if ( present(perm2) ) perm2 = p2
        if ( present(info) )  info  = inf1
        if (il1 > 0) print *,"LDU konci"
    end subroutine LDU
    
    !> destruktivni LDU rozklad
    !!
    !! nahradi vstupni matici jejim LDU rozkladem\n
    !! A = L*D*U\n
    !! L ... dolni trojuhelnikova matice s jednickami na diagonale\n
    !! D ... diagonalni matice\n
    !! U ... horni trojuhelnikova matice s jednickami na diagonale\n
    !!
    !! pro volbu pivota plati nasledujici pravidla\n
    !! 0 ... zadna\n
    !! 1 ... uplny vyber hlavniho prvku (nutne obe permutace )\n
    !! 2 ... sloupcovy vyber hlavniho prvku(nutno perm1)\n
    !! 3 ... radkovy vyber hlavniho prvku(nutno perm2)\n
    !! 4 ... diagonalni vyber hlavniho prvku(nutno aspon jedna = obe shodne)\n
    !! 5 ... diagonalni vyber hlavniho prvku pro minimal degree\n
    !!              (nutno aspon jedna = obe shodne) \n
    !! 6 ... rid se predepsanou permutaci\n
    !! Pro situace 1 az 6 musi byt pritomny permutacni vektory\n
    !! Pokud nejsou pozadavky splneny - zvoli se nejpodobnejsi volba\n
    !!
    !! pro podrobnost informaci plati\n
    !! 0 ... tichy chod\n
    !! 1 ... super podrobne\n
    !<
    subroutine LDUd(A,info,pivtype,ilev,perm1,perm2)
        use mtx
        use typy
        use pmatools
        implicit none
        !> rozkladana matice, pote s rozkladem
        class(matrix), intent(inout) :: A
        !> pocty operaci
        type(tcount), intent(inout), optional :: info
        !> typ pivotaze
        integer, intent(in), optional :: pivtype
        !> podrobnost informaci
        integer, intent(in), optional :: ilev
        !> permutace pro i
        integer(kind=ikind), dimension(:), allocatable, &
        intent(inout), optional :: perm1
        !> permutace pro j
        integer(kind=ikind), dimension(:), allocatable,&
        intent(inout), optional :: perm2
        
        
        real(kind=rkind) :: t1,t2 !meze pro cas
        integer :: pt1
        integer :: il1
        type(tcount) :: inf1
        integer(kind=ikind) :: i,q,j,k
        integer(kind=ikind) :: piv1, piv2 ! souradnice pivota
        integer(kind=ikind), dimension(:), allocatable :: p1, pp1
        integer(kind=ikind), dimension(:), allocatable :: p2, pp2
        integer(kind=ikind),dimension(:), allocatable :: ri1
        integer(kind=ikind),dimension(:), allocatable :: ri2
        real(kind=rkind),dimension(:), allocatable :: r1
        real(kind=rkind),dimension(:), allocatable :: r2
        integer(kind=ikind) :: sr1
        integer(kind=ikind) :: sr2
        real(kind=rkind) :: wrk, wrk1
        logical, dimension(:), allocatable :: mp1
        logical, dimension(:), allocatable :: mp2
        !! data pro pivotaci
        integer(kind=ikind) :: pivnelem, wrkpiv, ii1
        integer(kind=ikind),dimension(:), allocatable :: pivjj
        real(kind=rkind), dimension(:), allocatable :: pivval
        integer(kind=ikind), dimension(:), allocatable :: degs, &
        dstart, dnext, dprev, dstartnext
        integer(kind=ikind) :: ip1,ip2,ip3,ip4, degsmin, lastmin
        
        call cpu_time(t1)
        lastmin = 0
        allocate(p1(1:A%getn()))
        allocate(p2(1:A%getm()))
        allocate(pp1(1:A%getn()))
        allocate(pp2(1:A%getm()))
        allocate(mp1(1:A%getn()))
        allocate(mp2(1:A%getm()))
        do i=1,A%getn()
            p1(i)=i
        end do
        do i=1,A%getm()
            p2(i)=i
        end do
        !! p1(1:i-1) jiz pouzite indexy pivotu
        !! p1(i:n) jeste nepouzite
        !! p2 totez pro sloupce
        !! pp1 kde je ktery index v p1, vlastne inverzni permutace
        pp1 = p1 ! na zacatku je to stejne
        pp2 = p2
        
        mp1 = .true.
        mp2 = .true.
        if ( present(pivtype)) then
            pt1 = pivtype
            else
            pt1 = 0
        end if
        if (present(ilev)) then
            il1 = ilev
            else
            il1 = 0
        end if
        
        select case(pt1)
            case(0)
            ! nic resit nemusim
            if (il1>0) print *,"bez vyberu hlavniho prvku"
            case (1)
            if (present(perm1)) then
                if (present(perm2)) then
                    ! tak tohle je ok
                    if (il1>0) print *,"s uplnym vyberem hlavniho prvku"
                    else
                    ! mam k dispozici je radkove permutace
                    ! menim na sloupcovy vyber
                    pt1 = 2
                    if (il1>0) print *,&
                    "menim uplnym vyber hlavniho prvku na sloupcovy"
                end if
                else
                if (present(perm2)) then
                    ! mam jen sloupcove permuatace
                    ! menim na radkovy vyber
                    pt1 = 3
                    if (il1>0) print *,&
                    "menim uplnym vyber hlavniho prvku na radkovy"
                    else
                    ! nemam zadne permutace
                    ! menim na zadny vyber
                    pt1 = 0
                    if (il1>0) print *,&
                    "menim uplnym vyber hlavniho prvku na bez vyberu"
                end if
            end if
            case (2)
            if (present(perm1)) then
                ! OK
                else
                ! menim na bez vyberu
                pt1 = 0
            end if
            case (3)
            if (present(perm2)) then
                ! OK
                else
                ! menim na bez vyberu
                pt1 = 0
            end if
            case (4)
            if (present(perm1) .or. present(perm2)) then
                ! OK
                else
                ! menim na bez vyberu
                pt1 = 0
            end if
            case (5)
            if (present(perm1) .or. present(perm2)) then
                ! OK
                else
                ! menim na bez vyberu
                pt1 = 0
            end if
            case (6)
            if (present(perm1)) p1=perm1
            if (present(perm2)) p2=perm2
            case default
            print *, "chybna volba pivota, menim na zadna pivotaz"
            pt1 = 0
        end select
        i = A%getn()
        q = A%getm()
        if (i < q) q = i
        !! spocitej stupne
        allocate(degs(1:A%getn()))
        allocate(dstart(1:A%getn()))
        allocate(dnext(1:A%getn()))
        allocate(dprev(1:A%getn()))
        allocate(dstartnext(1:A%getn()))
        degs   = -1
        dstart = -1
        dnext  = -1
        dprev  = -1
        dstartnext = -1
        degsmin = 0

        do i=1, A%getn()
            !print *,"i=",i
            call A%getrow(i,r1,ri1,sr1,mp1)
            !print *,"pred adddegree", sr1
            if (sr1 == 0) sr1 = sr1 + 1
            call adddegree(i,sr1)
            !print *, "sr1=", sr1
        end do
        !print *,"step3"
        do i=1,q-1
            if (il1 > 10) print *, "step=",i
            ! napred najdi pivota
            if (il1 > 10) print *,"hledam pivota"
            call findpivot
            if (il1 > 10) print *, "mam pivota", piv1, piv2
            if ( .not. mp1(piv1) ) stop "podivne mp1"
            if ( .not. mp2(piv2) ) stop "podivne mp2"
            mp1(piv1) = .false.
            mp2(piv2) = .false.
            !! tady zajistit, ze v piv je to usporadane
            ip1       = p1(i)
            ip2       = pp1(piv1)
            p1(i)     = piv1 ! tohle je nutne
            p1(ip2)   = ip1
            pp1(ip1)  = ip2
            pp1(piv1) = i
            
            ip1       = p2(i)
            ip2       = pp2(piv2)
            p2(i)     = piv2 ! aby to vubec mohlo bezet
            p2(ip2)   = ip1
            pp2(ip1)  = ip2
            pp2(piv2) = i
            if (il1>10) print *,"LDUd jdu pro radek "
            call a%getrow(piv1,r1,ri1,sr1,mp2)
            ! v r1 jsou cisla z piv. radku, v ri1 jsou sloupcove indexy
            if (il1 > 0 ) print *,"step=",i," pivi=",piv1, " pivj=", &
            piv2, "stupen=",sr1
            if (il1>10) print *,"LDUd jdu pro sloupec"
            call a%getcol(piv2,r2,ri2,sr2,mp1)
            ! v r2 jsou cisla z piv. sloupce v ri2 jsou radkove indexy
            wrk = a%get(piv1,piv2)
            wrk1 = wrk
            if (wrk==0) STOP "Pruser"
            if (il1 > 1) print *, "mam data",sr1,sr2,wrk
            do j=1,sr1
                r1(j) = r1(j)/wrk
            end do
            inf1%div = inf1%div+sr1
            if (il1 > 10) print *,"jdu eliminovat"
            do j=1,sr1 ! jdi pres sloupce
                do k=1,sr2 ! jdi pres radky
                    if (il1 > 10) print *,j,k,sr1,sr2  !!
                    wrk = a%get(ri2(k),ri1(j))         !! prezkoumat
                    wrk = wrk - r1(j)*r2(k)            !!
                    call a%set(wrk,ri2(k),ri1(j))      !!
                end do
            end do
            if (il1 > 10) print *, "krok hotov"
            wrk = wrk1
            do j=1,sr1
                call a%set(r1(j),piv1,ri1(j))
            end do
            do j=1,sr2
                call a%set(r2(j)/wrk,ri2(j),piv2)
            end do
            !! prepocitat stupne
            do j=1,sr1
                call A%getrow(ri1(j),r2,ri2,sr2,mp1)
                degs(ri1(j)) = sr2
            end do
            if (il1 > 10) print *," radek a sloupec ulozeny"
            inf1%ad  = inf1%ad  + sr1*sr2
            inf1%mul = inf1%mul + sr1*sr2
            inf1%div = inf1%div + sr2
        end do
        !! doplnit permutace
        do i=1, a%getn()
            if(mp1(i)) then
                p1(a%getn())=i
                exit
            end if
        end do
        do i=1, a%getm()
            if(mp2(i)) then
                p2(a%getm())=i
                exit
            end if
        end do
        
        call cpu_time(t2)
        if (present(info)) then
            info%ad   = info%ad  + inf1%ad
            info%mul  = info%mul + inf1%mul
            info%div  = info%div + inf1%div
            info%time = t2 - t1
        end if
        if ( present(perm1) ) perm1 = p1
        if ( present(perm2) ) perm2 = p2
        if (il1 > 10) print *,"LDUd konci"
        contains
        !> vlozi informaci o stupni do seznamu
        subroutine adddegree(row,deg)
            implicit none
            integer(kind=ikind) :: row, deg
            degs(row) = deg
            ! zjisti jestli existuje odpovidajici stupen
            if (dstart(deg) > 0) then
                ! pokud ano, zarad na zacatek seznamu
                dnext(row) = dstart(deg)
                dprev(row) = -1
                dprev(dstart(deg)) = row
                dstart(deg) = row
                else
                ! pokud ne, vytvor novy seznam a uprav ukazatel
                ! na nejnizsi stupen
                dstart(deg) = row
                dnext(row)  = -1
                dprev(row) = -1
                degsmin = deg
            end if
        end subroutine adddegree
        
        !> vyjme informaci o stupni ze zasobniku
        subroutine deldegree(row)
            implicit none
            integer(kind=ikind) :: row
            
            if (dprev(row) == -1 ) then
                ! pokud jsi na zacatku seznamu
                if (dnext(row) == -1) then
                    ! pokud jsi jediny, tak se odstran a ukld v seznamech
                    if ( degs(row) == degsmin) then
                        ! uvolnuje se nejnizi stupen
                        else
                        ! uvolnuje se uprostred
                    end if
                    else
                    ! jinak se jen odeber, nic pred, neco za
                    ! nutno upravit zacatk stupne
                    dstart(degs(row)) = dnext(row)
                    degs(row) = -1
                    dprev(dnext(row)) = -1
                    dnext(row) = -1
                    dprev(row) = -1
                end if
                else
                ! pokud ne, jen se odeber, neco je pred, mozna neco za
                dnext(dprev(row)) = dnext(row)
                if (dnext(row) > 0 ) then
                    ! neco ja za tebou
                    dprev(dnext(row)) = dprev(row)
                end if
                dnext(row) = -1
                dprev(row) = -1
                degs(row)  = -1
            end if
        end subroutine deldegree
        
        
        
        !> vybere pivota a zaznemena to v permutacich
        !!
        !! orientuje se podle promenne pt1 - ta se nastavuje podle pivtype
        !<
        subroutine findpivot
            implicit none
            integer(kind=ikind) :: wrk, i1, j1, mi, mj,i2,i3,md
            real(kind=rkind) :: w1, mx
            select case(pt1)
                case (0)
                piv1 = i
                piv2 = i
                case (1)
                if (il1 > 10) print *, "piv1 i=",i
                mx = 0
                mi = 0
                mj = 0
                !do i1 = i, a%getn()
                do i3 = i, a%getn()
                    i1 = p1(i3)
                    i2 = i1
                    if ( .not. mp1(i2) ) stop "probem v uv"
                    call a%getrow(i2,pivval,pivjj,&
                    pivnelem,mp2)
                    do j1=1,pivnelem
                        w1 = abs(pivval(j1))
                        if (w1 > mx) then
                            if (il1>10) print *, "menim"
                            mx = w1
                            mi = i2
                            mj = pivjj(j1)
                        end if
                    end do
                end do
                piv1 = mi
                piv2 = mj
                case (2)
                ! sloupcovy vyber
                piv2 = i
                call a%getcol(i,pivval,pivjj, &
                pivnelem,mp1)
                mx = 0
                mi = 0
                mj = i
                do i1=1,pivnelem
                    w1 = abs(pivval(i1))
                    if ( w1 > mx ) then
                        mx = w1
                        mi = pivjj(i1)
                    end if
                end do
                piv1 = mi
                case (3)
                piv1 = i
                call a%getrow(i,pivval,pivjj, &
                pivnelem,mp2)
                mx = 0
                mi = i
                mj = 0
                do i1=1,pivnelem
                    w1 = abs(pivval(i1))
                    if ( w1 > mx ) then
                        mx = w1
                        mj = pivjj(i1)
                    end if
                end do
                piv2 = mj
                case (4)
                continue
                case (5)
                ! symetricke minimal degree
                !! napred jem primitivne
                md = a%getn() + 1
                mi = 0
                mj = 0
                !                    do i1 = 1,a%getn()
                !                        if (.not. mp1(i1)) cycle
                !                        call a%getrow(i1,pivval,pivjj,&
                !                                      pivnelem,mp2)
                !                        if (pivnelem < md) then
                !                            md = pivnelem
                !                            mi = i1
                !                            mj = i1
                !                        end if
                !                    end do
                ! ted trochu lepe
                do j1=i, A%getn()
                    i1 = p1(j1)
                    if (.not. mp1(i1)) STOP "potiz v MD"
                    if (degs(i1) < md) then
                        md = degs(i1)
                        mi = i1
                        mj = i1
                        if (md <= lastmin) exit
                    end if
                end do
                piv1 = mi
                piv2 = mj
                lastmin = md
                case (6)
                piv1 = p1(i)
                piv2 = p2(i)
                case default
                stop "chybna volba pivotaze - tady bych nemel nikdy byt"
            end select
        end subroutine findpivot
    end subroutine LDUd
    
    !> rozdeli matici s LDU rozkladem na jednotlive komponenty
    !!
    !! pokud jsou dodany permutace, tak plati
    !<
    subroutine Split(LDU,L,D,U,p1,p2)
        use typy
        use mtx
        !> LDU dekompozice
        class(matrix), intent(in) :: LDU
        !> L cast rozkladu
        class(matrix), intent(out) :: L
        !> D cast rozkladu
        class(matrix), intent(out) :: D
        !> U cast rozkladu
        class(matrix), intent(out) :: U
        !> Leva permutace
        integer(Kind=ikind), dimension(:), intent(in) :: p1
        !> Prava permutace
        integer(Kind=ikind), dimension(:), intent(in) :: p2
        integer(kind=ikind) :: n,m,i,j,i1,j1
        
        !print *,"split zacina"
        call L%clone(LDU)
        call D%clone(LDU)
        call U%clone(LDU)
        
        n=LDU%getn()
        m=LDU%getm()
        do i=1,min(n,m)
            call L%set(1.0_rkind,p1(i),p2(i))
            call U%set(1.0_rkind,p1(i),p2(i))
            call D%set(LDU%get(p1(i),p2(i)),p2(i),p1(i))
            do j=i+1,m
                i1 = p1(i)
                j1 = p2(j)
                call U%set(LDU%get(i1,j1),i1,j1)
            end do
            do j=i+1,n
                i1 = p1(j)
                j1 = p2(i)
                call L%set(LDU%get(i1,j1),i1,j1)
            end do
        end do
    end subroutine Split
    
    !> Zpetny chod po LDU dekompozici.
    !! Predpoklada se, ze matice je ctvercova.
    !<
    subroutine LDUback(LDU,b,x,p1,p2,oc1,ilevel,errcode)
        use typy
        use mtx
        implicit none
        !> rozlozena matice
        class(matrix), intent(in) :: LDU
        !> prava strana
        real(kind=rkind), dimension(:), intent(in)  :: b
        real(kind=rkind), dimension(:), intent(out) :: x
        integer(kind=ikind), dimension(:), intent(in), optional :: p1
        integer(kind=ikind), dimension(:), intent(in), optional :: p2
        type(tcount), intent(out), optional :: oc1
        !> info level, nepovinny parametr, default = 0
        !! - 0 ... pracuj tise
        !! - 1 ... vse reportuj
        integer, intent(in), optional :: ilevel
        !> chybovy status, nepovinny. Pokud neni  pritomen,tak chyba 
        !! zpusobi konec programu
        !! - 0 ... vse OK
        !! - 1 ... chyba - blize neurcena
        !! - 2 ... lisi se delka b a pocet radku matice
        !! - 3 ... lisi se delka x a pocet sloupcu matice
        !! - 4 ... matice je singularni a soustava je konzistentni
        !! - 5 ... matice je singularni a soustava nekonzistentni
        integer, intent(out), optional :: errcode
        type(tcount) :: oc
        integer(kind=ikind), dimension(:), allocatable :: perm1, perm2
        integer(kind=ikind), dimension(:), allocatable :: iperm1, iperm2
        integer(kind=ikind) :: n,m,i,j, ip, ip2
        real(kind=rkind), dimension(1:LDU%getn()) :: diag
        real(kind=rkind), dimension(:), allocatable :: v
        integer(kind=ikind), dimension(:), allocatable :: ii
        integer(kind=ikind) :: nelem
        logical, dimension(1:LDU%getn()) :: mask
        integer :: ilv
        
        if ( present(ilevel)) then
            ilv = ilevel
            else
            ilv = 0
        end if
        
        n = LDU%getn()
        m = LDU%getm()
        ! vyresit optionaly
        if (present(p1)) then
            allocate(perm1(1:n))
            perm1 = p1
            else
            allocate(perm1(1:n))
            do i=1,n
                perm1(i) = i
            end do
        end if
        if (present(p2)) then
            allocate(perm2(1:m))
            perm2 = p2
            else
            allocate(perm2(1:m))
            do i=1,m
                perm2(i) = i
            end do
        end if
        ! konstrukce inverznich permutaci
        allocate(iperm1(1:n),iperm2(1:m))
        do i=1,n
            iperm1(perm1(i))=i
        end do
        do i=1,m
            iperm2(perm2(i))=i
        end do
        
        if (ubound(b,1) /= n) then
            if ( present(errcode) ) then
                else
                errcode = 2
                stop "LDUback: Pocet radku matice a prave strany se lisi"
            end if
        end if
        
        if (ubound(x,1) /= n) then
            if ( present(errcode) ) then
                else
                errcode = 3
                stop "LDUback: Pocet sloupcu matice a neznamych se lisi"
            end if
        end if
        
        ! 1. L y = b
        x = b
        if ( ilv==1) then
            print *, "LDUback: zacatek"
            print *, "ridici permutace"
            print *, perm1
            print *, perm2
            print *, "Prava strana"
            print *,x
        end if
        diag(perm2(1))=LDU%get(perm1(1),perm2(1))
        mask = .false.
        mask(perm2(1)) = .true.
        do i = 2,n
            ip = perm1(i)
            ip2 = perm2(i)
            diag(ip2)= LDU%get(ip,ip2)
            call LDU%getrow(ip,v,ii,nelem,mask)
            if (ilv==1) then
                print *, "i=",i, ip, nelem, mask
                print *, ii(1:nelem)
                print *, v(1:nelem)
            end if
            do j=1,nelem
                if ((ii(j)) == ip) then
                    ! jsem na diagonale
                    diag(ip2) = v(j)
                    else
                    if (ilv == 1) then
                        print *, j, x(ip), v(j), x(ii(j)), ii(j)
                    end if
                    x(ip) = x(ip) - v(j)*x(ii(j))
                     if (ilv == 1) then
                        print *, j, x(ip), v(j), x(ii(j)), ii(j)
                    end if
                end if
            end do
            if (ilv == 1) then
                print *,perm1
                print *, perm2
                print *,ip,mask
                print *,ii
                print *,x
            end if
            mask(ip2) = .true.
        end do
        if (ilv == 1) then
            print *, "po primem chodu -------------------------------------"
            print *,x
        end if
        ! 2. D z = y
        x = x/diag  ! tady zalezi na poradi, potreba prezkoumat
        if (ilv==1) then
            print *, "po deleni diagonalou --------------------------------"
            print *,x
            print *, diag
        end if
        ! 3. U x = z
        mask = .false.
        mask(perm2(n)) = .true.
        do i=n-1,1,-1
            ip = perm1(i)
            ip2 =perm2(i)
            call LDU%getrow(ip,v,ii,nelem,mask)
            do j=1,nelem
                if (ii(j) == ip2) then
                    ! jsem na diagonale
                    diag(ip) = v(j)
                    else
                    x(ip) = x(ip) - v(j)*x(ii(j))
                end if
            end do
            if (ilv==1) then
            print *,perm1
            print *,ip,mask
            print *,ii
            print *,x
            end if
            mask(ip2) = .true.
        end do
        
    end subroutine LDUback
    
    !> \page sd "Metoda nejvetsiho spadu"
    !! Metoda je zalozena na lokalni minimalizaci funkcionalu
    !! \f$ F(x)=x^{T}Ax-2x^{T}b  \f$ . Minimalizace je provadena ve smeru
    !! gradientu, coz je residuum. Vysledny postup je \n
    !! \f[ x_{k+1}=x_k+\frac{r_{k}^{T}r_k}{r_{k}^{T}Ar_k}r_k  \f]
    !! kde \f$ r_k=b-Ax_ k \f$ .
    !<
    
    !> metoda nejvetsiho spadu ( Steepest descent) podrobnosti viz \ref sd
    subroutine SD(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
        ll1,ll2,cond1,opcnt1, errcode1)
        use mtx
        use typy
        implicit none
        !> matice soustavy\n
        !! musi poskytovat getn, getm, mul (nasobeni vektorem)
        class(matrix), intent(in) :: A
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
        !> kod pripadne chyby
        !! - 0 ... OK
        !! - 1 ... matice neni ctvercova
        !! - 2 ... nesouhlasi b
        !! - 3 ... nesouhasi x
        !! - 4 ... ani jeden z vektoru nesouhlasi
        !! - 5 ... vypocet ukoncen pred dosazenim pozadovane presnosti na itmax
        !! - 6 ...  -||-, ale pro ukonceni poklesu energiie i residua
        integer, intent(out), optional :: errcode1
        
        integer(kind=ikind) :: itmax, cnt, n
        integer :: ilev,rcount, rcount1, errcode
        real(kind=rkind) :: reps, repsfin, r0, r1, r2,l1,l2,t1,t2
        real(kind=rkind) :: wrk, en, enmin, res2min
        real(kind=rkind), dimension(:), allocatable :: res
        type(tcount) :: opcnt
        
        call cpu_time(t1)
        errcode = 0
        n = A%getn()
        rcount = 0
        rcount1 = 0
        l1 = 1
        l2 = 1
        if (present(ilev1)) then
            ilev = ilev1
            else
            ilev = 0
        end if
        
        if (present(itmax1)) then
            itmax = itmax1
            else
            itmax = A%getn()
            if (ilev==10) print *,&
            "nastavuji default hodnotu pro itmax na",itmax
        end if
        if (A%getn() /= A%getm()) then
            if (present(errcode1)) then
                errcode1 = 1
                if (ilev>0) print *,"SD: Matice neni ctvercova"
                return
                else
                Stop "SD: Matice neni ctvercova"
            end if
        end if
        if (A%getn() /= ubound(b,1)) then
            if (present(errcode1)) then
                errcode1 = 2
                if (ilev>0) print *,"SD : nesouhlasi delka prave strany"
                else
                Stop "SD: Nesouhlasi delka prave strany"
            end if
        end if
        if (A%getm() /= ubound(x,1)) then
            if (present(errcode1)) then
                if (errcode1==0) then
                    errcode1 = 3
                    if (ilev>0) print *,"SD: nesouhlasi delka reseni"
                    else
                    errcode1 = 4
                    if (ilev>0) print *,"SD: nesouhlasi delky obou vektoru"
                end if
                else
                Stop "SD: Nesouhlasi delka vektoru neznamych"
            end if
        end if
        if (present(errcode1)) then
            if (errcode1>0) return
        end if
        if (present(reps1) ) then
            reps = reps1
            else
            reps = 1.0e-6
            if (ilev==10) print *,"nastavuji reps na defaultni hodnotu "&
            , reps
        end if
        ! ted uz je snad vyrizeno vsechno ze vstupu
        ! jdu pocitat
        cnt = 0
        do
            res = b - A%mul(x,opcnt)
            r1 = dot_product(res,A%mul(res,opcnt))
            r2 = dot_product(res,res)
            wrk = r2/r1
            ! tenhle usek je jen kvuli vypoctu energie
            en = -dot_product(x,b+res) ! je to vychozi energie
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            x = x + wrk * res
            opcnt%ad = opcnt%ad + 4*n-2 +2*n-1
            opcnt%mul = opcnt%mul + 3*n + n
            opcnt%div = opcnt%div + 1
            
            cnt = cnt + 1
            if (cnt==1) then
                enmin = en
                r0 = r2
                res2min = r2
                if (r0==0.0_rkind) then
                    repsfin=0
                    exit
                end if
            end if
            if (en < enmin) then
                rcount = 0
                enmin = en
                else
                rcount = rcount+1
            end if
            if (r2 < res2min) then
                res2min = r2
                rcount1 = 0
                else
                rcount1 = rcount1+1
            end if
            if (ilev>0)  print *, cnt,sqrt(r2), reps, en, enmin
            if (cnt == itmax) then
                errcode = 5
                exit
            end if
            repsfin = sqrt(r2/r0)
            if (repsfin < reps) exit
            if ((rcount1>5) .and. (rcount > 5)) then
                errcode = 6
                exit
            end if
        end do
        
        call cpu_time(t2)
        opcnt%time = t2-t1
        ! vyridit konec
        if (present(opcnt1)) call update_info(opcnt1,opcnt)
        if (present(errcode1)) errcode1 = errcode
        if (present(itfin1)) itfin1 = cnt
        if (present(repsfin1)) repsfin1 = repsfin
        if (present(ll1)) ll1 = l1  ! tohle zatim nepocita
        if (present(ll2)) ll2 = l2
        if (present(cond1)) cond1 = l1/l2
    end subroutine SD
    
    subroutine SDnormal
    end subroutine SDnormal
    
    !> \page cg Metoda sdruzenych gradientu
    !! metoda probiha podle nasledujiciho algoritmu
    !! -# \f$ r_0 = b-Ax_0\f$ , \f$ p_0 = r_0 \f$
    !! -# pro j=0, ... do konvergence
    !! -#   \f$ \alpha_j=\frac{(r_j , r_j)}{(p_j , Ap_j)} \f$
    !! -# \f$ x_{j+1} = x_j + \alpha_j p_j \f$
    !! -# \f$ r_{j+1} = r_j - \alpha_j Ap_j \f$
    !! -# \f$ \beta_j = \frac {(r_{j+1} , r_{j+1})} {(r_j , r_j)} \f$
    !! -# \f$ p_{j+1} = r_{j+1} + \beta_j *p_j\f$
    !! -# konec cyklu
    !!
    !! Odhad cisla podminenosti lze provest pomoci spocitanych
    !! udaju, viz. Saad, Iterative solution ... , str.181
    !! Vlastni cisla odhaduji vlastni cisla matice
    !! \f[
    !! \left(\begin{array}{cccccc}
    !! \frac{1}{\alpha_{0}} & \frac{\sqrt{\beta_{0}}}{\alpha_{0}} & 0 & 0 & \ldots & 0\\
    !! \frac{\sqrt{\beta_{0}}}{\alpha_{0}} & \frac{1}{\alpha_{1}}+\frac{\beta_{0}}{\alpha_{0}} & \frac{\sqrt{\beta_{1}}}{\alpha_{1}} & 0 & \ddots & \vdots\\
    !! 0 & \frac{\sqrt{\beta_{1}}}{\alpha_{1}} & \ddots & \ddots & \ddots & 0\\
    !! 0 & 0 & \ddots & \ddots & \ddots & 0\\
    !! \vdots & \ddots & \ddots & \ddots & \ddots & \frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}}\\
    !! 0 & \ldots & 0 & 0 & \frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} & \frac{1}{\alpha_{m-1}}+\frac{\beta_{m-2}}{\alpha_{m-2}}
    !! \end{array}\right)
    !! \f]
    !!
    !!
    !<
    
    !> metoda sdruzenych gradientu \ref cg
    subroutine CG(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
        ll1,ll2,cond1,opcnt1,errcode1)
        use mtx
        use typy
        implicit none
        !> matice soustavy\n
        !! musi poskytovat getn, getm, mul (nasobeni vektorem)
        class(matrix), intent(in) :: A
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
        !> kod pripadne chyby
        !! - 0 ... OK
        !! - 1 ... matice neni ctvercova
        !! - 2 ... nesouhlasi b
        !! - 3 ... nesouhasi x
        !! - 4 ... ani jeden z vektoru nesouhlasi
        !! - 5 ... vycerpan povoleny pocet iteraci
        !! - 6 ... prestalo klesat residuum i energie
        integer, intent(out), optional :: errcode1
        integer(kind=ikind) :: itmax, cnt, n
        integer :: ilev,rcount, rcount1, errcode
        real(kind=rkind) :: reps, repsfin, r0, r1, r2,r3,l1,l2,t1,t2
        real(kind=rkind) :: wrk, en, enmin, res2min
        real(kind=rkind), dimension(:), allocatable :: res,p,Ap,&
        alfa, beta
        type(tcount) :: opcnt
        
        call cpu_time(t1)
        errcode = 0
        n = A%getn()
        rcount = 0
        rcount1 = 0
        l1 = 1
        l2 = 1
        if (present(ilev1)) then
            ilev = ilev1
            else
            ilev = 0
        end if
        
        if (present(itmax1)) then
            itmax = itmax1
            else
            itmax = A%getn()
            if (ilev==10) print *,&
            "nastavuji default hodnotu pro itmax na",itmax
        end if
        allocate(alfa(0:itmax), beta(0:itmax))
        if (A%getn() /= A%getm()) then
            if (present(errcode1)) then
                errcode1 = 1
                if (ilev>0) print *,"CG: Matice neni ctvercova"
                return
                else
                Stop "CG: Matice neni ctvercova"
            end if
        end if
        if (A%getn() /= ubound(b,1)) then
            if (present(errcode1)) then
                errcode1 = 2
                if (ilev>0) print *,"CG : nesouhlasi delka prave strany"
                else
                Stop "CG: Nesouhlasi delka prave strany"
            end if
        end if
        if (A%getm() /= ubound(x,1)) then
            if (present(errcode1)) then
                if (errcode1==0) then
                    errcode1 = 3
                    if (ilev>0) print *,"CG: nesouhlasi delka reseni"
                    else
                    errcode1 = 4
                    if (ilev>0) print *,"CG: nesouhlasi delky obou vektoru"
                end if
                else
                Stop "CG: Nesouhlasi delka vektoru neznamych"
            end if
        end if
        if (present(errcode1)) then
            if (errcode1>0) return
        end if
        if (present(reps1) ) then
            reps = reps1
            else
            reps = 1.0e-6
            if (ilev==10) print *,"nastavuji reps na defaultni hodnotu "&
            , reps
        end if
        ! ted uz je snad vyrizeno vsechno ze vstupu
        ! jdu pocitat
        cnt = 0
        !! Step 1
        res = b - A%mul(x,opcnt)
        r2 = dot_product(res,res)
        r0 = r2
        p = res
        opcnt%ad  = n-1
        opcnt%mul = n
        !! Step 2
        do
            !! Step 3
            Ap =A%mul(p,opcnt)
            r1 = dot_product(p,Ap)
            wrk = r2/r1
            if (ilev>0) print *,"alfa=",wrk
            alfa(cnt) = wrk
            ! tenhle usek je jen kvuli vypoctu energie
            en = -dot_product(x,b+res) ! je to vychozi energie
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Step 4
            x = x + wrk*p
            !! Step 5
            res = res - wrk*Ap
            !! Step 6
            r3 = r2
            r2 = dot_product(res,res)
            wrk = r2/r3
            if (ilev>0) print *,"beta=",wrk
            beta(cnt) = wrk
            !! Step 7
            p = res + wrk*p
            opcnt%ad  = opcnt%ad  + 4*n-2 + 2*n-1
            opcnt%mul = opcnt%mul + 4*n   + n
            opcnt%div = opcnt%div + 1
            
            !! odhadnout vlastni cisla
            call esteig(l1,l2,cnt,alfa,beta)
            print *, "odhad vl. cisel:",l1,l2
            !! hotovo
            
            cnt = cnt + 1
            if (cnt==1) then
                enmin = en
                res2min = r2
                if (r0==0.0_rkind) then
                    repsfin=0
                    exit
                end if
            end if
            if (en < enmin) then
                rcount = 0
                enmin = en
                else
                rcount = rcount+1
            end if
            if (r2 < res2min) then
                res2min = r2
                rcount1 = 0
                else
                rcount1 = rcount1+1
            end if
            if (ilev>0)  print *, cnt,sqrt(r2), reps, en, enmin
            if (cnt == itmax) then
                errcode = 5
                exit
            end if
            repsfin = sqrt(r2/r0)
            if (repsfin < reps) exit
            if ((rcount1>5) .and. (rcount > 5)) then
                errcode = 6
                exit
            end if
        end do
        
        call cpu_time(t2)
        opcnt%time = t2-t1
        ! vyridit konec
        if (present(opcnt1)) opcnt1 = opcnt
        if (present(errcode1)) errcode1 = errcode
        if (present(itfin1)) itfin1 = cnt
        if (present(repsfin1)) repsfin1 = repsfin
        if (present(ll1)) ll1 = l1
        if (present(ll2)) ll2 = l2
        if (present(cond1)) cond1 = l1/l2
        
        
    end subroutine CG
    
    !> \page cgn Metoda sdruzenych gradientu pro normalni rovnice
    !! metoda probiha podle nasledujiciho algoritmu
    !! -# \f$ r_0 = A^T (b-Ax_0 )\f$ , \f$ p_0 = r_0 \f$
    !! -# pro j=0, ... do konvergence
    !! -#   \f$ \alpha_j=\frac{(r_j , r_j)}{(Ap_j , Ap_j)} \f$
    !! -# \f$ x_{j+1} = x_j + \alpha_j p_j \f$
    !! -# \f$ r_{j+1} = r_j - \alpha_j A^T Ap_j \f$
    !! -# \f$ \beta_j = \frac {(r_{j+1} , r_{j+1})} {(r_j , r_j)} \f$
    !! -# \f$ p_{j+1} = r_{j+1} + \beta_j *p_j\f$
    !! -# konec cyklu
    !!
    !! Odhad cisla podminenosti lze provest pomoci spocitanych
    !! udaju, viz. Saad, Iterative solution ... , str.181
    !! Vlastni cisla odhaduji vlastni cisla matice
    !! \f[
    !! \left(\begin{array}{cccccc}
    !! \frac{1}{\alpha_{0}} & \frac{\sqrt{\beta_{0}}}{\alpha_{0}} & 0 & 0 & \ldots & 0\\
    !! \frac{\sqrt{\beta_{0}}}{\alpha_{0}} & \frac{1}{\alpha_{1}}+\frac{\beta_{0}}{\alpha_{0}} & \frac{\sqrt{\beta_{1}}}{\alpha_{1}} & 0 & \ddots & \vdots\\
    !! 0 & \frac{\sqrt{\beta_{1}}}{\alpha_{1}} & \ddots & \ddots & \ddots & 0\\
    !! 0 & 0 & \ddots & \ddots & \ddots & 0\\
    !! \vdots & \ddots & \ddots & \ddots & \ddots & \frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}}\\
    !! 0 & \ldots & 0 & 0 & \frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} & \frac{1}{\alpha_{m-1}}+\frac{\beta_{m-2}}{\alpha_{m-2}}
    !! \end{array}\right)
    !! \f]
    !!
    !!
    !<
    
    !> metoda sdruzenych gradientu pro normalni rovnice \ref cgn
    subroutine CGnormal(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
        ll1,ll2,cond1,opcnt1,errcode1)
        use mtx
        use typy
        implicit none
        !> matice soustavy\n
        !! musi poskytovat getn, getm, mul (nasobeni vektorem)
        class(matrix), intent(in) :: A
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
        integer(kind=ikind) :: itmax, cnt, n
        integer :: ilev,rcount, rcount1, errcode
        real(kind=rkind) :: reps, repsfin, r0, r1, r2,r3,l1,l2,t1,t2
        real(kind=rkind) :: wrk, en, enmin, res2min
        real(kind=rkind), dimension(:), allocatable :: res,p,Ap,&
        alfa, beta,&
        res1, b1
        type(tcount) :: opcnt
        
        call cpu_time(t1)
        errcode = 0
        n = A%getn()
        rcount = 0
        rcount1 = 0
        l1 = 1
        l2 = 1
        if (present(ilev1)) then
            ilev = ilev1
            else
            ilev = 0
        end if
        
        if (present(itmax1)) then
            itmax = itmax1
            else
            itmax = A%getn()
            if (ilev==10) print *,&
            "nastavuji default hodnotu pro itmax na",itmax
        end if
        allocate(alfa(0:itmax), beta(0:itmax))
        if (A%getn() /= ubound(b,1)) then
            if (present(errcode1)) then
                errcode1 = 2
                if (ilev>0) print *,"CGnormal : nesouhlasi delka prave strany"
                else
                Stop "CGnormal: Nesouhlasi delka prave strany"
            end if
        end if
        if (A%getm() /= ubound(x,1)) then
            if (present(errcode1)) then
                if (errcode1==0) then
                    errcode1 = 3
                    if (ilev>0) print *,"CGnormal: nesouhlasi delka reseni"
                    else
                    errcode1 = 4
                    if (ilev>0) print *,"CGnormal: nesouhlasi delky obou vektoru"
                end if
                else
                Stop "CGnormal: Nesouhlasi delka vektoru neznamych"
            end if
        end if
        if (present(errcode1)) then
            if (errcode1>0) return
        end if
        if (present(reps1) ) then
            reps = reps1
            else
            reps = 1.0e-6
            if (ilev==10) print *,"nastavuji reps na defaultni hodnotu "&
            , reps
        end if
        ! ted uz je snad vyrizeno vsechno ze vstupu
        ! jdu pocitat
        cnt = 0
        !! Step 1
        b1 = A%mult(b)
        res = b - A%mul(x,opcnt)
        res = A%mult(res,opcnt)
        r2 = dot_product(res,res)
        r0 = r2
        p = res
        opcnt%ad  = n-1
        opcnt%mul = n
        !! Step 2
        do
            !! Step 3
            Ap =A%mul(p,opcnt)
            r1 = dot_product(Ap,Ap)
            wrk = r2/r1
            if (ilev>0) print *,"alfa=",wrk
            alfa(cnt) = wrk
            ! tenhle usek je jen kvuli vypoctu energie
            en = -dot_product(x,b1+res) ! je to vychozi energie
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Step 4
            x = x + wrk*p
            !! Step 5
            res = res - wrk*A%mult(Ap,opcnt)
            !! Step 6
            r3 = r2
            r2 = dot_product(res,res)
            wrk = r2/r3
            if (ilev>0) print *,"beta=",wrk
            beta(cnt) = wrk
            !! Step 7
            p = res + wrk*p
            opcnt%ad  = opcnt%ad  + 4*n-2 + 2*n-1
            opcnt%mul = opcnt%mul + 4*n   + n
            opcnt%div = opcnt%div + 1
            
            !! odhadnout vlastni cisla
!             call esteig(l1,l2,cnt,alfa,beta)
!             print *, "odhad vl. cisel:",l1,l2
            !! hotovo
            
            cnt = cnt + 1
            if (cnt==1) then
                enmin = en
                res2min = r2
                if (r0==0.0_rkind) then
                    repsfin=0
                    exit
                end if
            end if
            if (en < enmin) then
                rcount = 0
                enmin = en
                else
                rcount = rcount+1
            end if
            if (r2 < res2min) then
                res2min = r2
                rcount1 = 0
                else
                rcount1 = rcount1+1
            end if
            if (ilev>0)  print *, cnt,sqrt(r2), reps, en, enmin
            if (cnt == itmax) then
                errcode = 5
                exit
            end if
            repsfin = sqrt(r2/r0)
            if (repsfin < reps) exit
            if ((rcount1>5) .and. (rcount > 5)) then
                errcode = 6
                exit
            end if
        end do
        
        call cpu_time(t2)
        opcnt%time = t2-t1
        ! vyridit konec
        if (present(opcnt1)) opcnt1 = opcnt
        if (present(errcode1)) errcode1 = errcode
        if (present(itfin1)) itfin1 = cnt
        if (present(repsfin1)) repsfin1 = repsfin
        if (present(ll1)) ll1 = l1  ! tohle zatim nepocita
        if (present(ll2)) ll2 = l2
        if (present(cond1)) cond1 = l1/l2
        
        
        
    end subroutine CGnormal
    
    
    !> \page pcg Metoda predpodminenych sdruzenych gradientu
    !! metoda probiha podle nasledujiciho algoritmu
    !! -# \f$ r_0 = b-Ax_0\f$ , \f$ z_0 = M^{-1} r_0 \f$ , \f$ p_0 = z_0 \f$
    !! -# pro j=0, ... do konvergence
    !! -# \f$ \alpha_j=\frac{(r_j , z_j)}{(p_j , Ap_j)} \f$
    !! -# \f$ x_{j+1} = x_j + \alpha_j p_j \f$
    !! -# \f$ r_{j+1} = r_j - \alpha_j Ap_j \f$
    !! -# \f$ z_{j+1} = M^{-1} r_{j+1} \f$
    !! -# \f$ \beta_j = \frac {(r_{j+1} , z_{j+1})} {(r_j , z_j)} \f$
    !! -# \f$ p_{j+1} = r_{j+1} + \beta_j p_j\f$
    !! -# konec cyklu
    !!
    !! Odhad cisla podminenosti lze provest pomoci spocitanych
    !! udaju, viz. Saad, Iterative solution ... , str.181
    !! Vlastni cisla odhaduji vlastni cisla matice
    !! \f[
    !! \left(\begin{array}{cccccc}
    !! \frac{1}{\alpha_{0}} & \frac{\sqrt{\beta_{0}}}{\alpha_{0}} & 0 & 0 & \ldots & 0\\
    !! \frac{\sqrt{\beta_{0}}}{\alpha_{0}} & \frac{1}{\alpha_{1}}+\frac{\beta_{0}}{\alpha_{0}} & \frac{\sqrt{\beta_{1}}}{\alpha_{1}} & 0 & \ddots & \vdots\\
    !! 0 & \frac{\sqrt{\beta_{1}}}{\alpha_{1}} & \ddots & \ddots & \ddots & 0\\
    !! 0 & 0 & \ddots & \ddots & \ddots & 0\\
    !! \vdots & \ddots & \ddots & \ddots & \ddots & \frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}}\\
    !! 0 & \ldots & 0 & 0 & \frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} & \frac{1}{\alpha_{m-1}}+\frac{\beta_{m-2}}{\alpha_{m-2}}
    !! \end{array}\right)
    !! \f]
    !!
    !!
    !<
    
    !> metoda sdruzenych gradientu \ref pcg
    subroutine PCG(A,Minv,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
        ll1,ll2,cond1,opcnt1,errcode1)
        use mtx
        use typy
        implicit none
        !> matice soustavy\n
        !! musi poskytovat getn, getm, mul (nasobeni vektorem)
        class(matrix), intent(in) :: A
        !> matice predpodmineni - ve skutecnosti operator\n
        !! musi poskytovat mul (nasobeni vektorem)
        class(matrix), intent(in) :: Minv
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
        integer(kind=ikind) :: itmax, cnt, n
        integer :: ilev,rcount, rcount1, errcode
        real(kind=rkind) :: reps, repsfin, r0, r1, r2,r3,r4,l1,l2,t1,t2
        real(kind=rkind) :: wrk, en, enmin, res2min
        real(kind=rkind), dimension(:), allocatable :: res,p,Ap,&
        alfa, beta, z
        type(tcount) :: opcnt
        
        call cpu_time(t1)
        errcode = 0
        n = A%getn()
        rcount = 0
        rcount1 = 0
        l1 = 1
        l2 = 1
        if (present(ilev1)) then
            ilev = ilev1
            else
            ilev = 0
        end if
        
        if (present(itmax1)) then
            itmax = itmax1
            else
            itmax = A%getn()
            if (ilev==10) print *,&
            "nastavuji default hodnotu pro itmax na",itmax
        end if
        allocate(alfa(0:itmax), beta(0:itmax))
        if (A%getn() /= A%getm()) then
            if (present(errcode1)) then
                errcode1 = 1
                if (ilev>0) print *,"CG: Matice neni ctvercova"
                return
                else
                Stop "CG: Matice neni ctvercova"
            end if
        end if
        if (A%getn() /= ubound(b,1)) then
            if (present(errcode1)) then
                errcode1 = 2
                if (ilev>0) print *,"CG : nesouhlasi delka prave strany"
                else
                Stop "CG: Nesouhlasi delka prave strany"
            end if
        end if
        if (A%getm() /= ubound(x,1)) then
            if (present(errcode1)) then
                if (errcode1==0) then
                    errcode1 = 3
                    if (ilev>0) print *,"CG: nesouhlasi delka reseni"
                    else
                    errcode1 = 4
                    if (ilev>0) print *,"CG: nesouhlasi delky obou vektoru"
                end if
                else
                Stop "CG: Nesouhlasi delka vektoru neznamych"
            end if
        end if
        if (present(errcode1)) then
            if (errcode1>0) return
        end if
        if (present(reps1) ) then
            reps = reps1
            else
            reps = 1.0e-6
            if (ilev==10) print *,"nastavuji reps na defaultni hodnotu "&
            , reps
        end if
        ! ted uz je snad vyrizeno vsechno ze vstupu
        ! jdu pocitat
        cnt = 0
        !! Step 1
        res = b - A%mul(x,opcnt)
        z = Minv%mul(res)
        r4 = dot_product(res,res)
        r2 = dot_product(res,z)
        r0 = r2
        p = z
        opcnt%ad  = n-1
        opcnt%mul = n
        !! Step 2
        do
            !! Step 3
            Ap =A%mul(p,opcnt)
            r1 = dot_product(p,Ap)
            wrk = r2/r1
            if (ilev>0) print *,"alfa=",wrk
            alfa(cnt) = wrk
            ! tenhle usek je jen kvuli vypoctu energie
            en = -dot_product(x,b+res) ! je to vychozi energie
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Step 4
            x = x + wrk*p
            !! Step 5
            res = res - wrk*Ap
            z = Minv%mul(res)
            !! Step 6
            r3 = r2
            r4 = dot_product(res,res)
            r2 = dot_product(res,z)
            wrk = r2/r3
            if (ilev>0) print *,"beta=",wrk
            beta(cnt) = wrk
            !! Step 7
            p = res + wrk*p
            opcnt%ad  = opcnt%ad  + 4*n-2 + 2*n-1
            opcnt%mul = opcnt%mul + 4*n   + n
            opcnt%div = opcnt%div + 1
            
            !! odhadnout vlastni cisla
            !call esteig(l1,l2,cnt,alfa,beta)
            !print *, "odhad vl. cisel:",l1,l2
            !! hotovo
            
            cnt = cnt + 1
            if (cnt==1) then
                enmin = en
                res2min = r2
                if (r0==0.0_rkind) then
                    repsfin=0
                    exit
                end if
            end if
            if (en < enmin) then
                rcount = 0
                enmin = en
                else
                rcount = rcount+1
            end if
            if (r2 < res2min) then
                res2min = r2
                rcount1 = 0
                else
                rcount1 = rcount1+1
            end if
            if (ilev>0)  print *, cnt,sqrt(r2), reps, en, enmin
            if (cnt == itmax) then
                errcode = 5
                exit
            end if
            repsfin = sqrt(r2/r0)
            if (repsfin < reps) exit
            if ((rcount1>5) .and. (rcount > 5)) then
                errcode = 6
                exit
            end if
        end do
        
        call cpu_time(t2)
        opcnt%time = t2-t1
        ! vyridit konec
        if (present(opcnt1)) opcnt1 = opcnt
        if (present(errcode1)) errcode1 = errcode
        if (present(itfin1)) itfin1 = cnt
        if (present(repsfin1)) repsfin1 = repsfin
        if (present(ll1)) ll1 = l1
        if (present(ll2)) ll2 = l2
        if (present(cond1)) cond1 = l1/l2
        
        
    end subroutine PCG
    
    subroutine PCGnormal
    end subroutine PCGnormal
    
    function vycdet(x,a,b,n) result(y)
        use typy
        implicit none
        real(kind=rkind), intent(in) :: x
        real(kind=rkind), dimension(0:), intent(in) :: a
        real(kind=rkind), dimension(0:), intent(in) :: b
        integer(kind=ikind), intent(in) :: n
        real(kind=rkind) :: y, f1,f2,f3
        integer(kind=ikind) :: i
        
        f1 = 1/a(0)-x
        if (n==0) then
            y = f1
            return
        end if
        f2 = f1*(1/a(1)+b(0)/a(0)-x)-b(0)/a(0)/a(0)
        if (n==1) then
            y = f2
            return
        end if
        do i=2,n
            f3 = (1/a(i)+b(i-1)/a(i-1)-x)*f2-b(i-1)/a(i-1)/a(i-1)*f1
            f1 = f2
            f2 = f3
        end do
        y = f3
        
        
    end function vycdet
    
    subroutine esteig(l1,l2,n,alfa,beta)
        use typy
        implicit none
        real(kind=rkind), intent(in out) :: l1
        real(kind=rkind), intent(in out) :: l2
        integer(kind=ikind), intent(in) :: n
        real(kind=rkind), dimension(0:), intent(in) :: alfa
        real(kind=rkind), dimension(0:), intent(in) :: beta
        
        real(kind=rkind) :: x1,x2,x3,y1,y2,y3,d1,d2
        
        if (n==0) then
            l1 = 1/alfa(0)
            l2 = l1
            else
            !odhadneme nejmensi
            x1 = 0
            y1 = vycdet(x1,alfa,beta,n)
            x2 = l2
            y2 = vycdet(x2,alfa,beta,n)
            if (y1*y2>0) print *, "podivne"
            d1 =x2-x1
            d2 = 2*d1
            do while (d1 < d2)
                d2 = d1
                x3 = (x1+x2)/2
                y3 = vycdet(x3,alfa,beta,n)
                if (y3*y1 > 0) then
                    x1 = x3
                    y1 = y3
                    else
                    x2 = x3
                    y2 = y3
                end if
                d1 = x2 - x1
            end do
            l2 = x3
            !odhadneme nejvetsi
            x1 = l1
            y1 = vycdet(x1,alfa,beta,n)
            d1 = 0.00000000001
            x2 = x1+d1
            y2 = vycdet(x2,alfa,beta,n)
            do while (y1*y2>0)
                x1 = x2
                y1 = y2
                d1 = d1+d1
                x2 = x1+d1
                y2 = vycdet(x2,alfa,beta,n)
            end do
            d1 =x2-x1
            d2 = 2*d1
            do while (d1 < d2)
                d2 = d1
                x3 = (x1+x2)/2
                y3 = vycdet(x3,alfa,beta,n)
                if (y3*y1 > 0) then
                    x1 = x3
                    y1 = y3
                    else
                    x2 = x3
                    y2 = y3
                end if
                d1 = x2 - x1
            end do
            l1 = x3
        end if
        
    end subroutine esteig
    
    
end module solvers
