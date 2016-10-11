!> \file dd_tools.f90
!>\brief DD implementation
!!
!!  neco
!!
!!

!>
!! module for implementing some kind of Domain Decomposition methods
module dd_tools
    implicit none
    !definuje metody pro konstrukci oblasti a metod Domain Decomposition
    public :: pcg,  & !metoda sdruzenych gradientu
    Minv, & !provede aplikaci prepodminovace
    Precond1, Precond1_setup !predpodminovadla
    ! Precond1 ... realizuje Schwarzovskou Domain Decomposition
    public :: applySchwarzDD !realizuje jeden krok Schwarze
    public :: prepSchwarzDD !vytvori potrebna data pro Schwarzovu DD
    public :: DD_test !tester pro Domain Decomposition
    public :: NactiData
    public :: AddLevel
    public :: MTest
    public :: condestim



contains
    !>
    !> \brief Preconditioned conjugate gradients
    !!
    !!  resi soustavu A*x=b
    !!  pokud je matice oznacena jako normeq=.true.
    !!  resi ve skutecnosti soustavu AtA x = At*b
    !!  realizuje pcg
    !! \callgraph
    !! \callergraph
    !<
    subroutine pcg(A,bb,x,ptype,ilev,ierr, itmax, reps, info,normmul)
        use mytypes
        use typy
        implicit none
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! parametry
        !> system matrix - supposed be symmetric positive definite\n
        !> inout because of possibility of initalizing of preconditioner
        type(smtx), intent(in out)                  :: A
        !> right hand side
        real(kind=rkind), dimension(:), intent(in) :: bb
        !> solution, on input initial aproximation
        real(kind=rkind), dimension(:), intent(in out) :: x
        !> level of information, default = 0
        integer, intent(in), optional           :: ilev
        !> kind of preconditioner, default 0
        integer, intent(in), optional           :: ptype
        !> error message\n 0 .. OK\n
        !>                1 .. after itmax iterions is not founded sufficiently small relative error
        integer, intent(out), optional          :: ierr
        !> maximum allowed iterations, default = 500
        integer(kind=ikind), intent(in), optional           :: itmax
        !> wanted relative error, default=  1e-5
        real(kind=rkind), intent(in), optional  :: reps
        !> operations count
        type(info_type), intent(inout), optional  :: info
        !> zda pro normeq=true delat prenasobeni
        logical, optional, intent(in) :: normmul

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! lokalni promenne
        integer                                       :: ilev1, ptype1, itmax1
        integer                                       :: ierr1
        real(kind=rkind)                              :: reps1
        real(kind=rkind), dimension(:), allocatable :: ax,r,p,z,ap,b
        integer(kind=ikind)                           :: n,i
        real(kind=rkind)                              :: alfa, beta, rz, rzold
        real(kind=rkind)                              :: r0r0,rr, rp, enr, bnr
        type(info_type)                                :: info1
        logical                                       :: nm1
        real(kind=rkind)                              :: normp
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! vypocet - vyrizeni defaultu
        ierr1 = 0
        if ( .NOT. present(ptype)) then
            ptype1 = 0
        else
            ptype1 = ptype
        end if
        if (present(ilev)) then
            ilev1 = ilev
        else
            ilev1 = 0
        end if
        if (present(itmax)) then
            itmax1 = itmax
        else
            itmax1 = 500
        end if
        if (present(reps)) then
            reps1 = reps
        else
            reps1 = 1.0e-5
        end if
        if (present(normmul)) then
            nm1 = normmul
        else
            nm1 = .true.
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! vypocet opravdu pocitani
        ! matice muze byt principialne obdelnikova,
        ! pak se resi AtA - cili dulezity je pocet sloupcu
        n = A%col
        if (ilev1 > 0) then
            print *, "PCG zacina"
        end if
        allocate(ax(n),r(n),p(n),z(n),ap(n),b(n))
        if ( A%normeq .and. nm1) then
            call mulAtx(A,bb,b,info1)  ! prenasob transpozici
        else
            b = bb ! jen zkopiruj
        end if
        bnr = dot_product(b,b)
        !step1 r0=b-Ax0 z0=Minv r0 p0=z0
        if(ilev1>0) print *, "alokovano, n=", n
        call mulAx(A,x,ax,info1)
        if(ilev1>0) print *, "mam Ax"
        r = b - ax
        enr = dot_product(x,ax)-2*dot_product(b,x)
        if(ilev1>0) print *, "mam r"
        call Minv(A,r,z,ptype1,info1)
        if(ilev1>0) print *, "mam z"
        p = z
        if(ilev1>0) print *, "mam p"
        !step2 pro i=0,.....
        i = 0
        rz = dot_product(r,z)
        r0r0 = dot_product(r,r)
        if (ilev1 > 0) then
            print *, "jdu do cyklu, vstupni residuum a energie ",r0r0,enr
        end if
        do
            i = i+1
            !step3 alfaj=(rj,zj)/(Apj,pj)
            call mulAx(A,p,Ap,info1)
            rp = dot_product(r,p)
            alfa = dot_product(ap,p)
            if( alfa == 0) then
                print *, "pAp = 0 - podivne, koncim"
                exit
            end if
            alfa = rp/alfa
            !step4 xj+1 = xj+alfaj*pj
            x = x + alfa*p
            !step5 rj+1 = rj-alfaj*Apj, uziju primou definici
            call mulAx(A,x,ax,info1)
            r=b-ax
            enr = dot_product(x,ax)-2*dot_product(b,x)
            !r = r - alfa*ap
            !step6 zj+1 = Minv*rj+1
            call Minv(A,r,z,ptype1,info1)
            !step7 betaj+1 = (rj+1,zj+1)/(rj,zj)
            rzold = rz
            if (rz == 0) then
                print *,"soucin rz==0, koncim"
                exit
            end if
            rz = dot_product(r,z)
            rr = dot_product(r,r)
            beta = rz/rzold
            !step8 pj+1=zj+1+betaj*pj
            !print *," normy p a z pred opravou",norm2(p),norm2(z)
            p = z + beta*p
            ! zkusim normalizovat p
            normp = norm2(p)
            if (ilev1 > 0) print *, "norma p ",normp
            if (normp == 0 ) RETURN!stop "oprava nulove delky"
            p = p/normp
            rz = rz/normp
            normp = norm2(p)
            if (ilev1 > 0) print *, "norma p ",normp
            !step9 konec cyklu
            if ( ilev1 > 0) then
                print *, i, rz, rr, enr, rr/r0r0, bnr
            end if
            if ( i > itmax1) then
                ierr1 = 1
                exit
            end if
            if ( rr < reps1*reps1*r0r0) then
                exit
            end if
        end do
        deallocate(ax,r,p,z,ap,b)
        if (ilev1 > 0 ) then
            print *, "PCG konci"
        end if
        ! ted nasypeme pripadne navratove hodnoty
        if ( present(ierr) ) then
            ierr = ierr1
        end if
        if (present(info)) then
            info = info1
        end if
    end subroutine pcg

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! skutecna realizace predpodmineni
    ! ve skutecnosti jen rozskok
    !> \brief application of preconditioner
    !!
    !! if fact it is just call for appropriate method
    !<
    subroutine Minv(A,r,z,ptype,info)
        use mytypes
        implicit none
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !> Matrix of system\n
        !> inout because of possibility of initalizing of preconditioner
        type(smtx), intent(in out) :: A
        !> residual vector
        real(kind=rkind), dimension(:),intent(in) :: r
        !> preconditioned residual
        real(kind=rkind), dimension(:), intent(out) :: z
        !> info data
        type(info_type), intent(inout) :: info
        !> kind of preconditioner\n
        !>  0 ... do nothing, i.e. z=r\n
        !>  1 ... Schwarz type proconditioner
        integer, intent(in)          :: ptype

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! skutecny vypocet
        select case(ptype)
            case(0) ! zadne predpodmineni
                !print *, "case 0"
                z = r
            case(1) ! Schwarz DD
                !print *, "case 1"
                call Precond1(A,r,z)
            case default
                print *, "neznamy typ predpodmineni"
                stop "havarie v predpodmineni"
        end select
    end subroutine Minv

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Domain Decomposition Schwarzova typu s hrubou urovni
    !> \brief Schwarz type preconditioner
    !!
    !! do one step of multiplicative Schwarz type method
    !!   all necessary data are part of the matrix A
    !<
    subroutine Precond1(A,r,z)
        use mytypes
        implicit none
        !> system matrix
        type(smtx), intent(in out) :: A
        ! tohle je matice soustavy, musi byt spd
        !> residual vector
        real(kind=rkind), dimension(:),intent(in) :: r
        !> preconditioned residual vector
        real(kind=rkind), dimension(:), intent(out) :: z

        !napred pro jistotu zkonrolovat inicializaci
        !print *, "pc1"
        if ( .NOT.  associated(A%prcd) ) then
            !print *, "pc2"
            !musim inicializovat nic neinicializovano
            call Precond1_setup(A)
        elseif ( .NOT. associated(A%prcd%pc1) ) then
            !mozna uz nejake predpodmineni bylo, ale ne tohle
            !print *, "pc3"
            call Precond1_setup(A)
        elseif ( A%prcd%pc1%status /= 1) then
            ! jeste porad to neni uplne inicializovano
            !print *, "pc4"
            call Precond1_setup(A)
        end if
        !print *, "Precond1"
        !tak ted uz je to snad ok a muzu pocitat
        z = r
    end subroutine Precond1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !inicializace pro Domain Decomposition Schwarzova typu s hrubou urovni
    subroutine Precond1_setup(A,nd,coarse_level,overlap,smooth,symmetry_step,info)
        use mytypes
        use sparse_tools
        use pma_private_tools
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  parametry
        type(smtx), intent(in out)                    :: A
        integer(kind=ikind), intent(in), optional :: nd
        logical, intent(in), optional             :: coarse_level,&
            symmetry_step, &
            smooth
        integer(kind=ikind), intent(in), optional :: overlap
        type(info_type), intent(in out), optional :: info
        ! A            ... matice soustavy, musi byt spd
        ! nd           ... pocet podoblasti
        ! coarse_level ... zda delat hrubou uroven
        ! overlap      ... velikost prekryvu (pouze priblizna)
        ! smooth       ... zhladit prolongator?
        ! symmetry_step... delat i symetrizujici zpetny chod?
        ! info         ... udaje pro statistiku

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! lokalni promenne
        type(prcnd1), pointer :: wrk !pracovni odkaz na data o predpodminovaci
        integer(kind=ikind)   :: w1
        logical               :: w2
        type(info_type)       :: info1

        !print *, "aa"
        !napred pro jistotu zkontrolovat alokace
        if ( .NOT.  associated(A%prcd) ) then
            !musim alokovat zacatek, jen to, zbytek pozdeji
            !print *, "ab"
            allocate(A%prcd)
            A%prcd%pc1 => null()
            !Print *, "ac"
        end if
        if ( .NOT. associated(A%prcd%pc1) ) then
            !mozna uz nejake predpodmineni bylo, ale ne tohle
            !print *, "ad"
            allocate(A%prcd%pc1)
            wrk => A%prcd%pc1
            wrk%status = 0
            if (allocated(wrk%covering)) deallocate(wrk%covering)
            !print *, "ae"
        else
            !print *, "af"
            wrk => A%prcd%pc1
            !print *, "ag"
        end if
        if ( wrk%status /= 1) then
            ! jeste porad to neni uplne inicializovano
            ! a jdu na to
            !print *, "ah"
            if ( present(nd) ) then
                !print *, "ai"
                wrk%nd = nd
                !print *, "aj"
            else
                !print *, "ak"
                print *," zadej pocet podoblasti"
                read(unit=*,fmt=*) wrk%nd
                !print *, "al"
            end if
            print *, "pocet podoblasti=" ,wrk%nd
            if (present(overlap)) then
                wrk%overlap = overlap
            else
                print *, "zadej velikost prekryti"
                read(unit=*, fmt=*) wrk%overlap
            end if
            print *, "prekryv je ", wrk%overlap, " rad"
            if (present(coarse_level)) then
                wrk%coarse_level = coarse_level
            else
                call getyesno("provadet hrubou uroven?",wrk%coarse_level)
                !wrk%coarse_level=getyesno()
            end if
            if (wrk%coarse_level) then
                call getyesno("uzit zhlazeneho prolongatoru?",wrk%smooth)
            end if
            if ( present(symmetry_step)) then
                wrk%symmetry_step = symmetry_step
            else
                call getyesno("provadet zpetny chod pro zachovani symetrie?",&
                    wrk%symmetry_step)
            end if
            !ted opravdove inicializace
            !1. vytvorit disjunktni pokryti
            if (allocated(wrk%covering)) then
                deallocate(wrk%covering)
            end if
            allocate(wrk%covering(1:A%row))
            call disjoint_covering(A, wrk%nd, wrk%covering)
            !2. zkonstruovat opravdove oblasti
            !3. jde-li o zhlazeny prolongator - zkonstruovat
            !4. vybrat matice
            !5. rozlozit matice

            ! tohle az uplne na zaver
            wrk%status = 1
        end if
    end subroutine Precond1_setup




    !> implementuje jeden krok Schwarzovy metody
    subroutine applySchwarzDD(A,b,x)
        use typy
        use mytypes
        implicit none
        !> matice
        type(smtx), intent(inout) :: A
        !> prava strana
        real(kind=rkind), dimension(:), intent(in) :: b
        !> aproximace reseni
        real(kind=rkind), dimension(:), intent(out) :: x
        type(prcnd1), pointer :: wrk ! pouzijeme jako alias
        integer(kind=ikind) :: i, wn
        real(kind=rkind), dimension(:), allocatable :: wrkx,wrkb
        x = 0 ! vynulujeme
        if (.not. associated(A%prcd) ) then
            STOP "neni definovan predpodminovac"
        end if
        if ( .not. associated(A%prcd%pc1)) then
            STOP "neni definovan predpodminovac ve druhe urovni"
        end if
        wrk = A%prcd%pc1
        !prvni prubeh
        do i=1,wrk%nd
            if (wrk%domlist(i)%active) then
                wn = ubound(wrk%domlist(i)%domelem,1)
                allocate(wrkx(1:wn),wrkb(1:wn))
                wrkx = 0
                wrkb = b(wrk%domlist(i)%domelem) ! vezmi si svuj kus prave strany
                call pcg(A=wrk%domlist(i)%domA,bb=wrkb,x=wrkx)
                deallocate(wrkx,wrkb)
            end if
        end do
        if (wrk%coarse_level) then
            ! pro hrubou uroven
        end if
        if (wrk%symmetry_step) then
            ! symetrizace
            do i=wrk%nd,1,-1
                if (wrk%domlist(i)%active) then
                    wn = ubound(wrk%domlist(i)%domelem,1)
                    allocate(wrkx(1:wn),wrkb(1:wn))
                    wrkx = 0
                    wrkb = b(wrk%domlist(i)%domelem) ! vezmi si svuj kus prave strany
                    call pcg(A=wrk%domlist(i)%domA,bb=wrkb,x=wrkx)
                    deallocate(wrkx,wrkb)
                end if
            end do
        end if
    end subroutine applySchwarzDD

    !> pripravi Schwarzovsky DD predpodminovac
    !!
    !<
    subroutine prepSchwarzDD(A,covering,overlap)
        use typy
        use mytypes
        use pma_private_tools
        !> opracovavana matice
        type(smtx), intent(inout) :: A
        integer(kind=ikind), dimension(:), intent(in) :: covering
        integer(kind=ikind), intent(in) :: overlap
        integer(kind=ikind), dimension(:), allocatable :: wrklist, fullset

        integer(kind=ikind) :: nd,n,i,j,k,cnt


        print *, "jsem v priprave"
        !! spociyam domeny
        n = A%row
        allocate(wrklist(1:n))
        allocate(fullset(1:n))
        do i=1,n
            fullset(i) = i
        end do
        nd = 0
        do i=1,n
            if (covering(i)>nd) then
                nd = covering(i)
            end if
        end do
        print *, "pocet domen je ",nd
        !! ted to zacnu nastavovat
        allocate(A%prcd)
        allocate(A%prcd%pc1)
        A%prcd%pc1%nd = nd
        A%prcd%pc1%status = 1
        A%prcd%pc1%overlap = overlap
        A%prcd%pc1%coarse_level = .false. !tohle jen prozatim
        A%prcd%pc1%smooth = .false.       ! -||-
        A%prcd%pc1%symmetry_step = .true.
        allocate(A%prcd%pc1%covering(1:n))
        A%prcd%pc1%covering = covering
        allocate(A%prcd%pc1%domlist(1:nd))
        call getyesno("uzit normalnich rovnic?",A%prcd%pc1%normeq)
        !! ted spocitam jednotlive domeny i s prekryvem
        do j=1,nd
            wrklist = 0
            k = 0
            do i=1,n
                if (covering(i) == j) then
                    wrklist(i) = 1
                    k = k + 1
                end if
            end do
            print *," domena c.:",j," pocatecni velikost=",k
            do i=1,overlap
                call AddLevel(A,wrklist)
            end do
            k = 0
            do i=1,n
                k = k + wrklist(i)
            end do
            print *,"finalni velikost=",k
            !ted to uloz
            allocate(A%prcd%pc1%domlist(j)%domelem(1:k))
            if (A%prcd%pc1%normeq) then
                k = 0
                do i=1,n
                    if (wrklist(i) == 1) then
                        k = k + 1
                        A%prcd%pc1%domlist(j)%domelem(k) = i
                    end if
                end do
                call GetSubmatrix(A,fullset,A%prcd%pc1%domlist(j)%domelem,&
                    A%prcd%pc1%domlist(j)%domA)
                A%prcd%pc1%domlist(j)%domA%normeq = .true.
            else
                call GetSubmatrix(A,A%prcd%pc1%domlist(j)%domelem,&
                    A%prcd%pc1%domlist(j)%domelem,&
                    A%prcd%pc1%domlist(j)%domA)
                A%prcd%pc1%domlist(j)%domA%normeq = .false.
            end if
        end do
        !ted vyrobim zhlazeny prolongator
    end subroutine prepSchwarzDD

    !> prida jeden level sousedu
    subroutine AddLevel(A,wrklist)
        use typy
        use mytypes
        implicit none
        !> matice
        type(smtx), intent(in) :: A
        !> mapa charakteristicke funkce
        integer(kind=ikind), dimension(:), intent(inout) :: wrklist
        integer(kind=ikind) :: i,j
        do i=1,A%row
            if (wrklist(i) == 0) then
                do j=A%rows(i,1),A%rows(i,2)
                    if (wrklist(A%jj(j))==1) then
                        wrklist(A%jj(j)) = 2
                    end if
                end do
            end if
        end do
        do i=1,A%row
            if (wrklist(i) == 2) then
                wrklist(i) = 1
            end if
        end do
        if (.not. A%normeq) then
            return
        end if
        ! ted jeste transponovana matice
        do i=1,A%row
            if (wrklist(i) == 1) then
                do j=A%rows(i,1),A%rows(i,2)
                    if (wrklist(A%jj(j))==0) then
                        wrklist(A%jj(j)) = 2
                    end if
                end do
            end if
        end do
        do i=1,A%row
            if (wrklist(i) == 2) then
                wrklist(i) = 1
            end if
        end do
    end subroutine AddLevel

    !> nacte data pro deleni
    subroutine NactiData(A,deleni,prekryv)
        use typy
        use mytypes
        use pma_private_tools
        implicit none
        !> opracovavana matice
        type(smtx), intent(inout) :: A
        !> disjunktni deleni site
        integer(kind=ikind), dimension(:), pointer :: deleni
        !> velikost prekryvu
        integer(kind=ikind), intent(out) :: prekryv

        integer(kind=ikind) :: n, nx, ny, i
        character(len=200) :: nazev
        integer :: choice
        logical :: volba

        call read_mtx(A)
        call getyesno("vytvorit formalne matici normalnich rovnic?",volba)
        A%normeq = volba
        print *, " pripravim predpodmineni"
        n = A%row
        allocate(deleni(1:n))
        deleni = 0
        print *, " zadej soubor s rozdelenim"
        read(unit=*,fmt=*) nazev
        !! ted to  prectu
        open(unit=20,file=nazev,status="old",action="READ")
        read(unit=20,fmt=*, IOSTAT=choice) nazev
        read(unit=20,fmt=*, IOSTAT=choice) nazev
        do
            read(unit=20,fmt=*, IOSTAT=choice) nx,ny
            if ( choice < 0 ) then
                exit
            end if
            !print *, nx,ny
            if (nx > 0) then
                if (nx > n) then
                    print *," toto by se stat nemelo - prilis vysoke cislo uzlu"
                else
                    if (deleni(nx)>0) then
                        if (deleni(nx) /= ny) then
                            print *,"problem - ignoruji"
                        end if
                    else
                        deleni(nx)=ny
                    end if
                end if
            else
                print *," zaporny uzel - ignoruji"
            end if
        end do
        close(20)
        !!ted to zkotroluji
        nx = 0
        ny = 0
        do i=1,n
            if (deleni(i) > nx ) then
                nx = deleni(i)
            else if (deleni(i) == 0 ) then
                ny = ny + 1
            end if
        end do
        print *,"pocet domen:",nx," pocet chybnych cisel:",ny
        print *,"velikost prekryvu"
        read(unit=*,fmt=*, IOSTAT=choice) prekryv


    end subroutine NactiData


    !> tester pro domain decomposition
    subroutine DD_test()
        use typy
        use mytypes
        implicit none
        type(smtx) :: A
        integer(kind=ikind),dimension(:), pointer :: deleni
        integer(kind=ikind) :: prekryv

        !! zavolam
        call NactiData(A,deleni,prekryv)
        call prepSchwarzDD(A,deleni,prekryv)

    end subroutine DD_test

    !> odhadne cislo podminenosti
    !!
    !! je uzito mocninne metody
    !! pokud uz neni - provede se konstrukce normalni matice a spocitaji se singularni cisla
    !! pak se odmocni
    !<
    subroutine condestim(A,lmin,lmax,cnd)
        use typy
        use mytypes
        implicit none
        !> matice
        type(smtx), intent(inout) :: A
        !> odhad nejmensiho singularniho cisla
        real(kind=rkind), intent(out) :: lmin
        !> odhad nejvetsiho singularniho cisla
        real(kind=rkind), intent(out) :: lmax
        !> lmax / lmin
        real(kind=rkind), intent(out) :: cnd

        real(kind=rkind),dimension(1:A%col) :: v1,v2,v3
        real(kind=rkind) :: w1,w2,w3,lo1,lo2,d1,d2,d1o,d2o,di,dio,dmin
        integer(kind=ikind) :: it, rc
        logical :: stat_orig

        v1 = 0
        v2 = 0
        v1(1) = 1
        v2(1) = 1
        v3 = 0
        it = 0
        stat_orig = A%normeq
        A%normeq = .true.
        rc = 0
        do
            if (it > 1) then
                lo1 = lmax
                lo2 = lmin
            end if
            it = it + 1
            call mulAx(A,v1,v3)
            lmax = dot_product(v3,v3)/dot_product(v3,v1)
            v1 = v3/norm2(v3)
            call mulAx(A,v2,v3)
            v3 = lmax*v2-v3
            lmin = lmax  - dot_product(v3,v3)/dot_product(v3,v2)
            v2 = v3/norm2(v3)
            cnd = lmax/lmin
            print *, it, lmax,lmin,cnd, rc
            d1o = d1
            d2o = d2
            d1 = abs(lmax - lo1)
            d2 = abs(lmin - lo2)
            dio = di
            di = max(d1,d2)
            if (it > 2*(a%col+a%row)) then
                if (di < dmin) then
                    dmin = di
                    rc = 0
                else
                    rc = rc + 1
                    if (rc == 200) then
                        exit
                    end if
                end if
            else
                dmin = di
            end if
        end do
        A%normeq = stat_orig
        if (.not. stat_orig) then
            lmin = sqrt(lmin)
            lmax = sqrt(lmax)
            cnd  = sqrt(cnd)
        end if
    end subroutine condestim


  subroutine nonsymit(A,b,x,ptype,ilev,ierr, itmax, reps, info)
    use mytypes
    use typy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! parametry
    !> system matrix - symmetric positive definite\n
    !> inout because of possibility of initalizing of preconditioner
    type(smtx), intent(in out)                  :: A
    !> right hand side
    real(kind=rkind), dimension(:), pointer :: b
    !> solution, on input initial aproximation
    real(kind=rkind), dimension(:), intent(in out) :: x
    !> level of information, default = 0
    integer, intent(in), optional           :: ilev
    !> kind of preconditioner, default 0
    integer, intent(in), optional           :: ptype
    !> error message\n 0 .. OK\n
    !>                1 .. after itmax iterions is not founded sufficiently small relative error
    integer, intent(out), optional          :: ierr
    !> maximum allowed iterations, default = 500
    integer, intent(in), optional           :: itmax
    !> wanted relative error, default=  1e-5
    real(kind=rkind), intent(in), optional  :: reps
    !> operations count
    type(info_type), intent(inout), optional  :: info

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! lokalni promenne
    integer                                     :: ilev1, ptype1, itmax1
    integer                                     :: ierr1
    real(kind = rkind)                          :: reps1
    real(kind=rkind), dimension(:), allocatable :: ax,r,p,z,ap
    integer(kind=ikind)                         :: n,i
    real(kind=rkind)                            :: alfa, beta, rz, rzold
    real(kind=rkind)                            :: r0r0,rr, rp, rap, ra2
    type(info_type)                             :: info1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! vypocet - vyrizeni defaultu
    ierr1 = 0

    if ( .NOT. present(ptype)) then
        ptype1 = 0
    else
        ptype1 = ptype
    end if
    if (present(ilev)) then
        ilev1 = ilev
    else
        ilev1 = 0
    end if
    if (present(itmax)) then
        itmax1 = itmax
    else
        itmax1 = 50000
    end if
    if (present(reps)) then
        reps1 = reps
    else
        reps1 = 1.0e-7
    end if

    n = A%col
    allocate(ax(n),r(n),p(n),z(n),ap(n))

    do i=1,itmax

      call mulAx(A,x,ax,info1)
      r = b - ax
      call mulAx(A,r,ap, info1)
      rap = dot_product(r,ap)
      ra2 = dot_product(ap,ap)
      rr = dot_product(r,r)
      if (ilev == 1) then
	print *, "norma rezidua^2", rr
      end if
      if (i == 1) then
	r0r0 = rr
      end if
      !       if (rr <= reps1*reps1*r0r0 .or. abs(rr) <= epsilon(rr)) then
      if ( abs(rr) <= 1e-2 ) then
	EXIT
      end if
      if (i == itmax .and. ilev == 1) then
	print *, "solver to nedal"
	STOP
      end if

      alfa = rap/ra2
      x = x + alfa*r

    end do



    deallocate(ax,r,p,z,ap)
    if (present(ierr)) then
      ierr = ierr1
    end if
  end subroutine nonsymit

    !>
  !> \brief Preconditioned conjugate gradients
  !!
  !!
  !!
  !!  realizuje pcg s normalnimi rovnicemi
  !<
  subroutine pcg_normal(A,b,x,ptype,ilev,ierr, itmax, reps_rel, info, normul)
    use mytypes
    use typy
    use pma_private_tools
    implicit none
    !use globals
      ! parametry
      !> system matrix - symmetric positive definite\n
      !> inout because of possibility of initalizing of preconditioner
      type(smtx), intent(in out)                     :: A
      !> right hand side
      real(kind=rkind), dimension(:), intent(in)     :: b
      !> solution, on input initial aproximation
      real(kind=rkind), dimension(:), intent(in out) :: x
      !> level of information, default = 0
      integer, intent(in), optional                  :: ilev
      !> kind of preconditioner, default 0
      integer, intent(in), optional                  :: ptype
      !> error message\n 0 .. OK\n
      !>                1 .. after itmax iterions is not founded sufficiently small relative error
      integer, intent(out), optional                 :: ierr
      !> maximum allowed iterations, default = 500
      integer(kind=ikind), intent(in), optional      :: itmax
      !> maximal relative error
      real(kind=rkind), intent(in), optional         :: reps_rel

      !> operations count
      type(info_type), intent(inout), optional       :: info
      !> zda pro normeq=true delat prenasobeni
      logical, optional, intent(in) :: normul

            !> maximal absolute error
      real(kind=rkind)                               :: reps_abs

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! lokalni promenne
    integer                                     :: ilev1, ptype1, itmax1
    integer                                     :: ierr1
    real(kind=rkind), dimension(:), allocatable :: ax,r,p,z,ap, tmp
    integer(kind=ikind)                         :: n,i,m
    real(kind=rkind)                            :: alfa, beta, rz, rzold
    real(kind=rkind)                            :: r0r0,rr, rp
    type(info_type)                             :: info1
    real(kind=rkind)                            :: tmp2, reps1, reps2
    logical :: normstatus
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! vypocet - vyrizeni defaultu
    ierr1 = 0
    if ( .NOT. present(ptype)) then
        ptype1 = 0
    else
        ptype1 = ptype
    end if
    if (present(ilev)) then
        ilev1 = ilev
    else
        ilev1 = 0
    end if
    if (present(itmax)) then
        itmax1 = itmax
    else
        itmax1 = 50000
    end if

    if (present(reps_rel)) then
        reps1 = reps_rel
    else
        reps1 = 1.0e-7
    end if

!     if (present(reps_abs)) then
!         reps2 = reps_abs
!     else
        reps2 = 1.0e-15
!     end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! vypocet opravdu pocitani

    normstatus = A%normeq
    A%normeq   = .false.
    n = A%col
    m = A%row
    if (ilev1 > 0) then
       print *, "PCG zacina"
    end if
    allocate(ax(m),r(n),p(n),z(n),ap(m), tmp(m))
    !step1 r0=b-Ax0 z0=Minv r0 p0=z0
    !print *, "alokovano, n=", n
!     print *, ubound(ax,1)
    call mulAx(A,x,ax,info1)
    !print *, "mam Ax"
    tmp = b - ax
    call mulAtx(A, tmp, r, info1)
    !print *, "mam r"
    call Minv(A,r,z,ptype1,info1)
    !print *, "mam z"
    p = z
    !print *, "mam p"
    !step2 pro j=0,.....
    i = 0
    rz = dot_product(r,z)
    r0r0 = dot_product(r,r)
    if (ilev1 > 0) then
      print *, "jdu do cyklu"
    end if
    do
      i = i+1
      !step3 alfaj=(rj,zj)/(Apj,pj)
      call mulAx(A,p,Ap,info1)
      rp = dot_product(r,p)
      tmp2 = dot_product(ap,ap)
      if (sqrt(abs(tmp2)) < epsilon(tmp2)) then
! 	    it = 0
        exit
      end if
      alfa = rp/tmp2
      !step4 xj+1 = xj+alfaj*pj
      x = x + alfa*p
      !step5 rj+1 = rj-alfaj*Apj, uziju primou definici
      call mulAx(A,x,ax,info1)
      tmp=b-ax
      call mulAtx(A, tmp, r, info1)
      !step6 zj+1 = Minv*rj+1
      call Minv(A,r,z,ptype1,info1)
      !step7 betaj+1 = (rj+1,zj+1)/(rj,zj)
      rzold =rz
      rz = dot_product(r,z)
      rr = dot_product(r,r)
      beta = rz/rzold
      !step8 pj+1=zj+1+betaj*pj
      p = z + beta*p
      !step9 konec cyklu
      if ( ilev1 > 0) then
        print *, i, rz, dot_product(r,r),  maxval(abs(r)), reps1, reps2
      end if
      if ( i > itmax1) then
        ierr1 = 1
! 	call write_log("WARNING!!!, pcg iterations exceeded", time)
!         print *, abs(rr), reps2, maxval(r)
	print *, "pcg iterations exceeded!!"
!         if (present(it)) it = i
        exit
      end if
      !       if ( rr < reps1*reps1*r0r0 .or. abs(rr) <= epsilon(rr) ) then
      if ( maxval(abs(r)) <= reps2 .or. abs(rr/r0r0) <=reps1) then
!         print *, i, rz, dot_product(r,r),  maxval(abs(r)), reps1, reps2, abs(rr/r0r0)
        ierr1 = 0
! 	print *, "pcg iterations:", i
!         if (present(it)) it = i
        exit
      end if
    end do
    deallocate(ax,r,p,z,ap, tmp)
    if (ilev1> 0 ) then
      print *, "PCG normal konci"
    end if
    ! ted nasypeme pripadne navratove hodnoty
    if ( present(ierr) ) then
        ierr = ierr1
    end if
    if (present(info)) then
        info = info1
    end if
    A%normeq = normstatus
  end subroutine pcg_normal


    !> testy pro Michala
    subroutine MTest()
        use typy
        use mytypes
        use pma_private_tools
        implicit none
        type(smtx) :: A1
        type(smtx) :: A2
        real(kind=rkind), dimension(:), pointer :: b1=> null()
        real(kind=rkind), dimension(:), pointer :: b2=> null()
        real(kind=rkind), dimension(:), pointer :: x1=> null()
        real(kind=rkind), dimension(:), pointer :: x2=> null()
        logical :: volba
        real(kind=rkind) :: nf,n1,emin,emax, lmin, lmax, cnd, rs1,rs2
        integer(kind=ikind) :: i,j
        type(info_type) :: info


        print *, " testy zacinaji"
        call Laplace2D(A1,6_ikind,7_ikind)

        call GetSubmatrix(A1,(/ (i,i=1,30,2) /),(/ (i,i=1,30,3) /),A2)
        call print_matrix(A2)
        call GetSubmatrix(A1,(/ 5_ikind,4_ikind /),(/ 1_ikind,2_ikind /),A2)
        call print_matrix(A2)
        call Dump_Matrix(A2)
        allocate(x1(1:42),b1(1:42))
        do i= 1,5
            x1=0
            b1=0
            x1(i)=1
            call mulAsx(A1,x1,b1,info)
            print *,b1
            !pause
        end do
        deallocate(x1,b1)
        call read_mtx(A1)
        print *, " matice prectena"
        call matnorm(A1,nf,n1,emin,emax)
        print *, nf,n1,emin,emax
        call getyesno("cist pravou stranu?", volba)
        if (volba) then
            call read_vct(b1)
        else
            allocate(b1(1:A1%row))
            b1 = 10
        end if
        allocate(x1(1:A1%col))
        x1 = 0
        rs1 = resid(A1,b1,x1)
        print *," pocatecni residuum =",rs1
        call pcg(A1,b1,x1,ilev=1, reps=1.0e-15_rkind, itmax=A1%row)
        rs2 = resid(A1,b1,x1)
        print *," pocatecni a konecne residuum =", rs1,rs2
        call write_vct(x1,filename="vektor1.vct")
        call getyesno("pokracovat?",volba)
        A1%normeq = .true.
        x1 = 0
        do
            rs1 = resid(A1,b1,x1)
            print *," pocatecni residuum =",rs1
            call pcg(A1,b1,x1,ilev=1, reps=1.0e-15_rkind, itmax=A1%row)
            rs2 = resid(A1,b1,x1)
            print *," pocatecni a konecne residuum =", rs1,rs2
            call getyesno("pokracovat restartem?",volba)
            if (.not. volba) exit
        end do
        call write_vct(x1,filename="vektor2.vct")

        print *,"nasleduje odhad cisla podminenosti"
        print *, "provedeme to mocninou metodou"
        call condestim(A1,lmin,lmax,cnd)
        print *, lmin,lmax,cnd
        call getyesno("pokracovat? Nasleduje pcg_normal",volba)

        x1 = 0
        call pcg_normal(A1,b1,x1,ilev=1, reps_rel=1.0e-15_rkind, itmax=A1%row)
        call getyesno("pokracovat?",volba)


        ! ted testneme akce s Laplacem
        call Laplace2d(A1,7_ikind,7_ikind)
        deallocate(b1,x1)
        allocate(b1(1:a1%col),x1(1:A1%col))
        b1 = 10
        x1 = 0
        do
            print *," pocatecni residuum =", resid(A1,b1,x1)
            call pcg(A1,b1,x1,ilev=1, reps=1.0e-15_rkind, itmax=A1%row)
            print *," konecne residuum =", resid(A1,b1,x1)
            call getyesno("pokracovat restartem?",volba)
            if (.not. volba) exit
        end do
        print *,"nasleduje odhad cisla podminenosti"
        print *, "provedeme to mocninou metodou"
        call condestim(A1,lmin,lmax,cnd)
        print *, lmin,lmax,cnd
        call getyesno("pokracovat?",volba)


        ! ted testneme akce s Laplacem ale normalizovanym
        !call Laplace2d(A1,5_ikind,5_ikind)
        deallocate(b1,x1)
        allocate(b1(1:a1%col),x1(1:A1%col))
        b1 = 10
        x1 = 0
        a1%normeq = .true.
        do
            print *," pocatecni residuum =", resid(A1,b1,x1)
            call pcg(A1,b1,x1,ilev=1, reps=1.0e-15_rkind, itmax=A1%row)
            print *," konecne residuum =", resid(A1,b1,x1)
            call getyesno("pokracovat restartem?",volba)
            if (.not. volba) exit
        end do
        print *,"nasleduje odhad cisla podminenosti"
        print *, "provedeme to mocninou metodou"
        call condestim(A1,lmin,lmax,cnd)
        print *, lmin,lmax,cnd
        call getyesno("pokracovat?",volba)

        call GetSubmatrix(A1,(/ (i,i=1,30,2) /),(/ (i,i=1,30,2) /),A2)
        call print_matrix(A2)
        do j=1,1000
            call GetSubmatrix(A1,(/ (i,i=1,30,2) /),(/ (i,i=7,30,3) /),A2)
            call print_matrix(A2)
            print *,j
        end do
    end subroutine MTest


end module dd_tools
