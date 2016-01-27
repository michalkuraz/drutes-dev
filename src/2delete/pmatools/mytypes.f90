!>\file mytypes.f90
!!\brief dodefinovane typy
!!
!! typy urcene pro praci s ridkou matici a pro sledovani vypoctu\n
!! !!!! zmena definice typu smtx !!!!
!<

!>
!!\brief dodefinovane typy
!!
!!typy urcene pro praci s ridkou matici a pro sledovani vypoctu
!!
!<
module mytypes
    !kind hodnoty pro bezna cisla
    use typy
    implicit none
    !vektor je obycejne realne pole
    !tj. typicka deklarace:
    ! real(kind=rkind), dimension(:) :: neco
    ! pripadne podle potreby s atributy allocatable nebo pointer
    !

    !>  maximalni povoleny pocet kusu prazdneho mista
    integer, parameter, public :: efsmax=20

    !>
    !!
    !! Type for sparse matrices
    !!
    !! Elements are stored rowvise,
    !! preconditioner and possible LDU decomposition
    !! are stored as a special member. These special members are initially empty.
    !! Start and end of each row is stored in vector rows.
    !! Empty spaces are stored in emptysets.
    !! Indices of diagonal elements are in diag.
    !! when triplet is part of empty set, is set tu zero
    !!
    type, public :: smtx
        !> number of rows
        integer(kind=ikind) :: row = 0
        !> number of collumns
        integer(kind=ikind) :: col = 0
        !> number of nonzeroes elements
        integer(kind=ikind) :: nz = 0
        !> size of allocated arrays
        integer(kind=ikind) :: allocsize = 0     !delky alokovanych poli
        !> amount of free sets
        integer(kind=ikind) :: efs = 0
        !> nekompromisne delat normalni rovnice
        logical :: normeq = .false.
        !> values
        real(kind=rkind), dimension(:), allocatable :: vals
        !> indices of diagonal elements
        integer(kind=ikind), dimension(:), allocatable :: diag
        !> indices of diagonal elements
        real(kind=rkind), dimension(:), allocatable :: weight
        !> row indices
        integer(kind=ikind), dimension(:), allocatable :: ii
        !> column indices
        integer(kind=ikind), dimension(:), allocatable :: jj
        !> starts and ends of rows - elements of row(i) are in <rows(i,1),rows(i,2)>
        integer(kind=ikind), dimension(:,:), allocatable :: rows
        !> interval of empty spaces
        integer(kind=ikind), dimension(efsmax,2)     :: emptysets
        !> preconditioner - initially empty
        type(precond), pointer :: prcd => null()
        !> LDU decomposition - initially epmty
        type(smtx), pointer  :: decomp => null()
        !> permutations used for LDU decomposition (:,1) - left, (:,2) - right
        integer(kind=ikind), dimension(:,:), allocatable :: perm
    end type smtx

      !>type for some informations of computations
    type, public :: info_type
        !> stores informations about number of additive operations
        integer(kind=ikind) :: add = 0
        !> stores informations about number of mutiplicative operations
        integer(kind=ikind) :: mul = 0
        !> stores informations about number of divisions
        integer(kind=ikind) :: div = 0
        !> time of computation
        real(kind=rkind)    :: time = 0
    end type info_type

    !>typ pro predpodmineni
    type, public :: precond
        !> just pointer to structure containing necessary information
        type(prcnd1), pointer :: pc1
    end type precond

    !> necessary data for preconditioning by DD
    type, public :: prcnd1
        !> status of preconditioner \n 0 ... empty \n 1 ... done
        integer(kind=ikind) :: status = 0
        ! status - je-li inicializovano
        ! 0 ... neni
        ! 1 ... vse hotovo
        !> amount of domain
        integer(kind=ikind) :: nd
        !> size of overlap - counted as amount aditional levels
        integer(kind=ikind) :: overlap
        !> coarse level flag \n true ... use coarse level \n false ... do not use
        logical  :: coarse_level
        !> create smoothed prolongator
        logical  :: smooth
        !> provide additional sweep to get symmetric preconditioner
        logical  :: symmetry_step
        !> cost of preparing of preconditioner
        type(info_type) :: setup_info
        !> disjoint covering of domain - obsahuje hodnoty 1..nd
        integer(kind=ikind), dimension(:), allocatable :: covering
        !> seznam jednotlivych domen
        type(listaccess), dimension(:), pointer :: domlist
        !> prolongator
        type(smtx) :: P
        !> zda uzit normalnich rovnic
        logical :: normeq
    end type prcnd1

    type, public :: listaccess
      integer(kind=ikind), dimension(:), pointer :: domelem
      type(smtx) :: domA
      logical    :: active
    end type  listaccess

    !verejna jmena procedur
    public :: init_matrix          !initializace matice
    public :: Laplace2D     !vytvori matici diskretizace
    !public :: Laplace3D     !vytvori matici diskretizace
    public :: sort_matrix   !setridi prvky matice, hlavne interni pouziti
    public :: print_matrix_sparse  !vytiskne matici
    public :: print_matrix  !vytiskne matici
    public :: Dump_Matrix   !v zasade ladici tisk matice
    public :: read_mtx      !precte matici ze souboru
    public :: read_vct      !precte vector ze souboru
    public :: write_mtx     !precte matici ze souboru
    public :: write_vct     !precte vector ze souboru
    public :: mulAx         !znasobi matici vektorem, y=A*x (nebo At*A*x)
    public :: mulAsx        !znasobi matici vektorem, y=A*x
    public :: mulAtx        !znasobi matici vektorem, y=At*x
    public :: transp        !transponuje ridkou matici
    !public :: spy           !vytiskne strukturu nenul v matici
    public :: GetSubmatrix
    public :: Sparse_Test_Driver
    public :: resid         !spocte residuum
    public :: matnorm       !spocte Frobeniovu a radkovou normu matice
  !> tiskne matice
  interface print_matrix
    module procedure print_matrix_sparse
  end interface print_matrix


    public :: clear_matrix
    contains



    !>vytvori matici s rozmery (n,m) a alokuje misto pro nz nenulovych prvku
    !!
    !!diagonala se pocita mezi nenuly vzdycky
    !<
    subroutine init_matrix(A,n,m,nz)
        !implicit none
        !> vytvarena matice
        type(smtx), intent(out) :: A
        !> pocet radku vytvarene matice
        integer(kind=ikind), intent(in) :: n
        !> pocet sloupcu vytvarene matice - nepovinne, pokud nezadano = n
        integer(kind=ikind), intent(in), optional :: m
        !> potencialni pocet nenulovych prvku - nepovinne, pokud nezadano = n+m
        integer(kind=ikind), intent(in), optional :: nz

        !> zvolena hodnota poctu sloupcu
        integer(kind=ikind) :: lm
        !> zvoleny potencialni pocet nenul
        integer(kind=ikind) :: lnz
        !> delka diagonaly
        integer(kind=ikind) :: dsz

        integer(kind=ikind) :: i

        !pripadne doplneni defaultnich hodnot
        if (present(m)) then
            lm = m
        else
            lm = n
        end if

        if (present(nz)) then
            lnz = nz
            if ( lnz < n ) then
              lnz = n + 1
            end if
        else
            lnz = n + lm
        end if

        !implementace
        dsz   = min(n,lm)
        A%row = n
        A%col = lm
        A%nz  = n  !v kazdem radku prave jeden prvek - diagonala, nebo posledni sloupec
        !alokace vektoru
        allocate(A%ii(lnz))
        allocate(A%jj(lnz))
        allocate(A%diag(n))
        allocate(A%weight(n))
        allocate(A%vals(lnz))
        allocate(A%rows(lnz,2))

        !vycisteni predpominovadel a rozkladu
        A%prcd   => null(A%prcd)
        A%decomp => null(A%decomp)
        if (allocated(A%perm)) deallocate(A%perm)
        ! ted inicializace
        A%efs = 1
        A%allocsize = lnz
        A%ii = -1
        A%jj = -1
        A%vals = 0
        ! napred nastav nulovou diagonalu
        do i = 1,dsz
          A%diag(i)   = i
          A%vals(i)   = 0
          A%ii(i)     = i
          A%jj(i)     = i
          A%rows(i,1) = i
          A%rows(i,2) = i
        end do
        ! pokud je vic radku nez sloupcu, dosad 0 do posledniho sloupce
        do i = dsz+1, n
           A%vals(i)   = 0
           A%ii(i)     = i
           A%jj(i)     = lm
           A%rows(i,1) = i
           A%rows(i,2) = i
        end do
        ! ted nastav volnou diru
        A%emptysets(1,1) = A%rows(n,2) +1
        A%emptysets(1,2) = lnz
    end subroutine init_matrix


    subroutine sort_matrix(A)
    ! na vstupu je matice s nastavenymi hodnotami vals, ii a jj
    ! pokud na danem miste prvej neni, obsahuje jak v ii tak v jj 0
    ! procedura doplni ostatni vektory, pripadne prekontroluje diagonalu
    ! a doplni rozdeleni mezer
    ! rozmery dat(t.j. row, col, nz, alocsize) jsou take v poradku
        type(smtx), intent(in out) :: A
        !lokalni promenne
        integer(kind=ikind) :: i,j,li,lj,lli,llj,k,j1
        real(kind=rkind)    :: lval
        logical :: rep
        !kod
        !jestli mam spravnou delku rows
        if (allocated(A%rows)) then
            if (ubound(A%rows,1) < (A%row)) then
                deallocate(a%rows)
                allocate(A%rows(1:A%row,1:2))
            end if
        else
           allocate(A%rows(1:A%row,1:2))
        end if
        !print *, "sort matrix 2"
        !napred si spocitam delky jednotlivych radku
        do i=1, A%row
            !print *,"sort matrix 2.01",i
            A%rows(i,:) = 0
        end do
        !print *,"sort matrix 2.1"
        ! tady ocekovat delky
        do i=1, A%allocsize
            if (A%ii(i) > 0) then
              A%rows(A%ii(i),2)=A%rows(A%ii(i),2)+1
            end if
        end do
        !ted ty delky nastavim
        !ocekovat delky
        !print *, "sort matrix 3"
        A%rows(1,1)=1
        do i=2, A%row
            A%rows(i,1) = A%rows(i-1,2) + 1
            A%rows(i,2) = A%rows(i,1) + A%rows(i,2) - 1
        end do
        !print *,"sort matrix 4"
        !ted to opravdu setridim
        !1. premistim radky na misto a setridim
        do i=1, A%nz
          !print *, "i=",i
          do
            rep = .false.
            li = A%ii(i) !skutecny radkovy index prvku
            lj = A%jj(i) !skutecny sloupcvy index prvku
            if ( (i < A%rows(li,1)) .or. (i > A%rows(li,2))) then
              do j=A%rows(li,1), A%rows(li,2)
                !projdi radek a pripadne zarad
    naseljsem: if (A%ii(j) /= li) then
                  !a) tak to prohod
                  !print *, "prohazuji",i,j, li, lj
                  lval = A%vals(i)
                  lli  = A%ii(i)
                  llj  = A%jj(i)
                  A%vals(i) = A%vals(j)
                  A%ii(i)   = A%ii(j)
                  A%jj(i)   = A%jj(j)
                  A%vals(j) = lval
                  A%ii(j)   = lli
                  A%jj(j)   = llj
                  !b) a zatrid
                  do j1=j,A%rows(li,1)+1,-1
                      !print *, "zkousim zatridit"
                      if (A%jj(j1) < A%jj(j1-1)) then
                          !mimo porad1 - zamen
                          !print *, "zatriduji"
                          lval = A%vals(j1)
                          lli  = A%ii(j1)
                          llj  = A%jj(j1)
                          A%vals(j1)   = A%vals(j1-1)
                          A%ii(j1)     = A%ii(j1-1)
                          A%jj(j1)     = A%jj(j1-1)
                          A%vals(j1-1) = lval
                          A%ii(j1-1)   = lli
                          A%jj(j1-1)   = llj
                      else
                        ! uz v poradi, tak skonci
                        exit
                      end if
                  end do
                  rep = .true. ! neco jsem menil, tak znovu stejne i
                  exit !skoncime cyklus pres j - prave jsem vlozil
                end if  naseljsem
              end do
          end if
            if ( .not. rep) then
              exit
            end if
        end do
      end do
        !print *, "sort matrix 5"
        !ted to pro kontrolu dotridim
        do i=1,A%row
            do j=A%rows(i,1), A%rows(i,2)
                do j1=A%rows(i,1)+1, A%rows(i,2)
                    if (A%jj(j1) < A%jj(j1-1)) then
                        !mimo porad1 - zamen
                        !print *, "podivne!!! zatriduji"
                        lval = A%vals(j1)
                        lli  = A%ii(j1)
                        llj  = A%jj(j1)
                        A%vals(j1)   = A%vals(j1-1)
                        A%ii(j1)     = A%ii(j1-1)
                        A%jj(j1)     = A%jj(j1-1)
                        A%vals(j1-1) = lval
                        A%ii(j1-1)   = lli
                        A%jj(j1-1)   = llj
                    end if
                end do
            end do
        end do
        ! a nastavim diagonalu
        do i=1, A%nz
          if ( A%ii(i) == A%jj(i) ) then
            A%diag(A%ii(i)) = i
          end if
        end do
    end subroutine sort_matrix

    subroutine print_matrix_sparse(A)
    ! vytiskna matici na standartni vystup
        type(smtx), intent(in) :: A
        !lokalni promenne
        integer(kind=ikind) :: i
        real(kind=rkind) :: wrk
        print *, "tisknu matici"
        print *, "pocet radku = ", A%row
        print *, "pocet sloupcu=", A%col
        print *, "pocet nenul  =", A%nz
        print *, "alokovany prostor pro prvky=", A%allocsize

        !ted  prvky
        do i=1,A%nz
          wrk = A%vals(i) !nejaky problem v printu dlouhych cisel
         if (wrk /= 0) then
           print "(a,i8,a,i8,a,es25.15)","i=",A%ii(i), " j=", A%jj(i), " A(i,j)=",wrk
         end if
        end do
        !zacatky radku
        !print *, "zacatky radku"
        !do i=1,A%row
        !    print "(a,i8,a,i8,a,i8)","radek=",i," zacatek", A%rows(i,1)," konec",&
        !    A%rows(i,2)
        !end do
    end subroutine print_matrix_sparse

    !>precte matici ze souboru a pripadne ji prealokuje
    subroutine read_mtx(a,fil,filename,ilev,ierr)
        !use mytypes
        ! parametry subroutiny
        !> ctena matice - povinny parametr
        type(smtx), intent(in out)         :: a
        !> cislo kanalu - nepovinne. Je-li pritomno spolecne s filename,
        !! musi popisovat identicky objekt
        integer, intent(in), optional      :: fil
        !> jmeno souboru - nepovinne
        character(len=*), intent(in), optional :: filename !jmeno souboru
        !pokud se vyskytuji oba, musi popisovat shodny objekt
        integer, intent(in), optional  :: ilev !uroven podrobnosti informaci
        ! 0 ... nic (default)
        integer, intent(out), optional :: ierr !chybova hlaseni
        ! 0 .... vse je OK (default)

        !lokalni promenne
        character(len=100)  :: filname1
        integer             :: ilev1, fil1, ierr1
        integer(kind=ikind) :: i, mi, cnt,li,lj
        real(kind=rkind)    :: lvars
        integer(kind=ikind), dimension(:), allocatable ::ldiag

        !napred vyridit defaulty
        if ( .NOT. present(ilev) ) then
            ilev1 = 0
        else
            ilev1 = ilev
        end if
        ierr1 = 0
        ! casem je to treba upresnit
        if ( .NOT. present(fil)) then
            !nemame popisovac souboru
            if (.NOT. present(filename)) then
                !nemame ani jmeno, musime precist
                print *, "zadej jmeno souboru s matici"
                read(unit=*,fmt=*) filname1
                print *,"jmeno souboru=",filname1
            else
                filname1=filename
            end if
            ! ted uz mam jmeno, jen otevru
            fil1 = 101
            open(unit=fil1, file=filname1, status = "old", action="read")
        else
           fil1 = fil
        end if
        read(unit=fil1,fmt=*) a%row, a%nz
        a%col = a%row
        print *, "velikost matice=",a%col," pocet nenul=", a%nz
        ! ted uvolnime pamet
        if ( allocated(a%diag)) then
            deallocate(a%diag)
        end if
        if (allocated(a%ii)) then
            deallocate(a%ii)
        end if
        if (allocated(a%jj)) then
            deallocate(a%jj)
        end if
        if (allocated(a%rows)) then
            deallocate(a%rows)
        end if
        if (allocated(a%vals)) then
            deallocate(a%vals)
        end if
        !alokujeme nove prostory
        allocate(a%diag(1:a%col))
        allocate(a%ii(1:a%nz))
        allocate(a%jj(1:a%nz))
        allocate(a%rows(1:a%row,1:2))
        allocate(a%vals(1:a%nz))
        allocate(ldiag(0:a%col))
        a%allocsize=a%nz
        a%prcd => null(a%prcd)
        !ted to precteme
        mi = 1 ! standartni zacatek pole je 1
        cnt = 0
        do i=1,a%nz
            read(unit=fil1, fmt=*) li, lj, lvars
            !print *, li,lj,lvars
            if (li < mi ) then !pokud je jiny tak ho ted oznac
                mi = li
            end if
            if (lj < mi ) then
                mi = lj
            end if
            if (li == lj) then
                ldiag(li) = i
            end if
              cnt = cnt + 1
              a%ii(cnt) = li
              a%jj(cnt) = lj
              a%vals(cnt) = lvars
              !print *, li,lj,lvars
        end do
        if (mi == 0) then
            a%ii = a%ii + 1
            a%jj = a%jj + 1
            a%diag =  ldiag(0:a%col-1)
        else
            a%diag = ldiag(1:a%col)
        end if
        a%nz = cnt
        !nz je opravdovy pocet nenul
        !vcetne diagonaly
        print *, "jdu tridit"
        call sort_matrix(a)
        print *, "setrideno"
        if ( .NOT. present(fil)) then
            !zavru to jen kdyz jsem to otevrel
            close(fil1)
        end if
        deallocate(ldiag)
        if ( present(ierr) ) then
           ierr = ierr1
        end if
    end subroutine read_mtx


    !precte vektor ze souboru a pripadne ji prealokuje
    subroutine read_vct(v,fil,filename,ilev,ierr)
        !use mytypes
        ! parametry subroutiny
        real(kind=rkind), dimension(:), pointer :: v
        !cteny vektor
        integer, intent(in), optional      :: fil
        !popisovac otevreneho souboru s matici
        character(len=*), intent(in), optional :: filename !jmeno souboru
        !pokud se vyskytuji oba, musi popisovat shodny objekt
        integer, intent(in), optional  :: ilev !uroven podrobnosti informaci
        ! 0 ... nic (default)
        integer, intent(out), optional :: ierr !chybova hlaseni
        ! 0 .... vse je OK (default)


        !lokalni promenne
        character(len=100)  :: filname1
        integer             :: ilev1, fil1,ierr1
        integer(kind=ikind) :: i, mi, cnt,li,lj
        real(kind=rkind)    :: lvars

        !napred vyridit defaulty
        if ( .NOT. present(ilev) ) then
            ilev1 = 0
        else
            ilev1 = ilev
        end if
        ierr1 = 0
        !Print *, "jsem v a"
        ! casem je to treba upresnit
        if ( .NOT. present(fil)) then
            !nemame popisovac souboru
            if (.NOT. present(filename)) then
                !nemame ani jmeno, musime precist
                print *, "zadej jmeno souboru s vektorem"
                read(unit=*,fmt=*) filname1
                print *,"jmeno souboru=",filname1
            else
                filname1=filename
            end if
            ! ted uz mam jmeno, jen otevru
            fil1 = 101
            open(unit=fil1, file=filname1, status = "old", action="read")
        else
           fil1 = fil
        end if
        if (associated(v)) then
            deallocate(v)
        end if
        read(unit=fil1, fmt=*) cnt
        allocate(v(1:cnt))
        do i=1,cnt
            read(unit=fil1,fmt=*) v(i)
        end do
        if ( .NOT. present(fil)) then
            !zavru to jen kdyz jsem to otevrel
            close(fil1)
        end if
        if ( present(ierr)) ierr = ierr1
     end subroutine read_vct


    subroutine write_mtx(a,fil,filename,ilev,ierr)
        !use mytypes
        ! parametry subroutiny
        type(smtx), intent(in out)         :: a
            !ctena matice, ridka
        integer, intent(in), optional      :: fil
         !popisovac otevreneho souboru s matici
        character(len=*), intent(in), optional :: filename !jmeno souboru
        !pokud se vyskytuji oba, musi popisovat shodny objekt
        integer, intent(in), optional  :: ilev !uroven podrobnosti informaci
        ! 0 ... nic (default)
        integer, intent(out), optional :: ierr !chybova hlaseni
        ! 0 .... vse je OK (default)

        !lokalni promenne
        character(len=100)  :: filname1
        integer             :: ilev1, fil1, ierr1
        integer(kind=ikind) :: i,j, mi, cnt,li,lj
        real(kind=rkind)    :: lvars

       print *, "jdu zapsat matici"
        !napred vyridit defaulty
        if ( .NOT. present(ilev) ) then
            ilev1 = 0
        else
            ilev1 = ilev
        end if
        ierr1 = 0
        ! casem je to treba upresnit
        if ( .NOT. present(fil)) then
            !nemame popisovac souboru
            if (.NOT. present(filename)) then
                !nemame ani jmeno, musime precist
                print *, "zadej jmeno souboru s matici"
                read(unit=*,fmt=*) filname1
                print *,"jmeno souboru=",filname1
            else
                filname1=filename
            end if
            ! ted uz mam jmeno, jen otevru
            fil1 = 101
            open(unit=fil1, file=filname1, status = "unknown", action="write")
        else
           fil1 = fil
        end if
        write(unit=fil1,fmt=*) a%row, a%nz
        print *, "velikost matice=",a%row," pocet nenul=", a%nz
        !ted to vypiseme
        do i=1,a%row
            do j=a%rows(i,1),a%rows(i,2)
            write(unit=fil1, fmt=*) a%ii(j), a%jj(j), a%vals(j)
            end do
        end do
        if ( .NOT. present(fil)) then
            !zavru to jen kdyz jsem to otevrel
            close(fil1)
        end if
    end subroutine write_mtx


    !> zapise vektor ze souboru
    subroutine write_vct(v,fil,filename,ilev,ierr)
        !use mytypes
        ! parametry subroutiny
        real(kind=rkind), dimension(:), intent(in) :: v
        !cteny vektor
        integer, intent(in), optional      :: fil
        !popisovac otevreneho souboru s matici
        character(len=*), intent(in), optional :: filename !jmeno souboru
        !pokud se vyskytuji oba, musi popisovat shodny objekt
        integer, intent(in), optional  :: ilev !uroven podrobnosti informaci
        ! 0 ... nic (default)
        integer, intent(out), optional :: ierr !chybova hlaseni
        ! 0 .... vse je OK (default)

        !lokalni promenne
        character(len=100)  :: filname1
        integer             :: ilev1, fil1, ierr1
        integer(kind=ikind) :: i, mi, cnt,li,lj
        real(kind=rkind)    :: lvars

        !napred vyridit defaulty
        ierr1 = 0
        if ( .NOT. present(ilev) ) then
            ilev1 = 0
        else
            ilev1 = ilev
        end if
        ierr1 = 0
        ! casem je to treba upresnit
        if ( .NOT. present(fil)) then
            !nemame popisovac souboru
            if (.NOT. present(filename)) then
                !nemame ani jmeno, musime precist
                print *, "zadej jmeno souboru s vektorem"
                read(unit=*,fmt=*) filname1
                print *,"jmeno souboru=",filname1
            else
                filname1=filename
            end if
            ! ted uz mam jmeno, jen otevru
            fil1 = 101
            open(unit=fil1, file=filname1, status = "unknown", action="write")
        else
           fil1 = fil
        end if
        cnt = ubound(v,1)
        write(unit=fil1, fmt=*) cnt
        do i=1,cnt
            write(unit=fil1,fmt=*) v(i)
        end do
        if ( .NOT. present(fil)) then
            !zavru to jen kdyz jsem to otevrel
            close(fil1)
        end if
        if (present(ierr)) then
           ierr = ierr1
        end if

     end subroutine write_vct

     
         !> zapise vektor ze souboru
    subroutine write_vct_i(v,fil,filename,ilev,ierr)
        !use mytypes
        ! parametry subroutiny
        integer(kind=ikind), dimension(:), intent(in) :: v
        !cteny vektor
        integer, intent(in), optional      :: fil
        !popisovac otevreneho souboru s matici
        character(len=*), intent(in), optional :: filename !jmeno souboru
        !pokud se vyskytuji oba, musi popisovat shodny objekt
        integer, intent(in), optional  :: ilev !uroven podrobnosti informaci
        ! 0 ... nic (default)
        integer, intent(out), optional :: ierr !chybova hlaseni
        ! 0 .... vse je OK (default)

        !lokalni promenne
        character(len=100)  :: filname1
        integer             :: ilev1, fil1, ierr1
        integer(kind=ikind) :: i, mi, cnt,li,lj
        real(kind=rkind)    :: lvars

        !napred vyridit defaulty
        ierr1 = 0
        if ( .NOT. present(ilev) ) then
            ilev1 = 0
        else
            ilev1 = ilev
        end if
        ierr1 = 0
        ! casem je to treba upresnit
        if ( .NOT. present(fil)) then
            !nemame popisovac souboru
            if (.NOT. present(filename)) then
                !nemame ani jmeno, musime precist
                print *, "zadej jmeno souboru s vektorem"
                read(unit=*,fmt=*) filname1
                print *,"jmeno souboru=",filname1
            else
                filname1=filename
            end if
            ! ted uz mam jmeno, jen otevru
            fil1 = 101
            open(unit=fil1, file=filname1, status = "unknown", action="write")
        else
           fil1 = fil
        end if
        cnt = ubound(v,1)
        write(unit=fil1, fmt=*) cnt
        do i=1,cnt
            write(unit=fil1,fmt=*) v(i)
        end do
        if ( .NOT. present(fil)) then
            !zavru to jen kdyz jsem to otevrel
            close(fil1)
        end if
        if (present(ierr)) then
           ierr = ierr1
        end if

     end subroutine write_vct_i
  !> udela y=A*x
  !!
  !! pokud je A%normeq==.true. tak ve skutecnosti nasobi AtA
  !<
  subroutine mulAx(A,x,y,info)
    use typy
    implicit none
    !> matice
    type(smtx), intent(in) :: A
    !> vektor x
    real(kind=rkind), dimension(:), intent(in)  :: x
    !> vektor y
    real(kind=rkind), dimension(:), intent(out) :: y
    !> udaje o poctu operaci
    type(info_type), intent(inout), optional :: info
    integer(kind=ikind) :: i,j

    real(kind=rkind), dimension(1:A%row) :: wrk
    type(info_type) :: info1
     ! tady otestovat korektnost dat
     wrk = 0
     y = 0
     call mulAsx(A,x,wrk,info1)
     if (A%normeq) then
       call mulAtx(A,wrk,y,info1)
     else
       y = wrk
     end if
     if (present(info)) then
        info%add  = info%add  + info1%add
        info%mul  = info%mul  + info1%mul
        info%div  = info%div  + info1%div
        info%time = info%time + info1%time
     end if
  end subroutine mulAx

  subroutine mulAtx(A,x,y,info)
    implicit none
    type(smtx), intent(in) :: A
    real(kind=rkind), dimension(:), intent(in)  :: x
    real(kind=rkind), dimension(:), intent(out) :: y
    type(info_type), intent(inout) :: info
    integer(kind=ikind) :: i,j

     ! tady otestovat korektnost dat
    if (A%col == ubound(y,1)) then
        if (A%row == ubound(x,1)) then
            ! O.K
        else
            stop "mulAtx: pocet radku matice nesouhlasi s delkou x"
        end if
    else
        stop "mulAtx: pocet sloupcu neni totozny s delkou y"
    end if
     y = 0
     !do i=1,A%row
     !  y(A%jj(A%rows(i,1):A%rows(i,2))) = y(A%jj(A%rows(i,1):A%rows(i,2))) &
     !   + x(i)*A%vals(A%rows(i,1):A%rows(i,2))
     !end do
     do i=1,a%row
        if (a%rows(i,1)>a%rows(i,2)) then
!             print *,"prazdny radek pro i:",i
        else
            do j=a%rows(i,1), a%rows(i,2)
                y(a%jj(j)) = y(a%jj(j)) + a%vals(j)*x(i)
            end do
        end if
     end do

  end subroutine mulAtx

  subroutine mulAsx(A,x,y,info)
    use typy
    implicit none
    type(smtx), intent(in) :: A
    real(kind=rkind), dimension(:), intent(in)  :: x
    real(kind=rkind), dimension(:), intent(out) :: y
    type(info_type), intent(inout), optional :: info
    integer(kind=ikind) :: i,j

     ! tady otestovat korektnost dat
    if (A%row == ubound(y,1)) then
        if (A%col== ubound(x,1)) then
            ! O.K
        else
            stop "mulAx: pocet sloupcu matice nesouhlasi s delkou x"
        end if
    else
        stop "mulAtx: pocet radku matice neni totozny s delkou y"
    end if
     !!!$omp parallel do
     y = 0
     do i=1,A%row
        if (A%rows(i,1)>A%rows(i,2)) then
!             print *,"prazdny radek pro i:",i
        else
            do j = A%rows(i,1),A%rows(i,2)
                y(i) = y(i) + A%vals(j)*x(A%jj(j))
            end do
        end if
     end do

  end subroutine mulAsx


  subroutine mulasx_jinak(A, x, y)
    use typy

    type(smtx), intent(in) :: A
    real(kind=rkind), dimension(:), intent(in) :: x
    real(kind=rkind), dimension(:), intent(out) :: y

    integer(kind=ikind) :: i,j

    y = 0

    do i=1, ubound(A%vals,1)
      if (A%ii(i) > 0) then
        y(A%ii(i)) = x(A%jj(i))*A%vals(i) + y(A%ii(i))
      end if
    end do
  end subroutine mulasx_jinak

subroutine mulAx1(A,x,y,info)
     type(smtx), intent(in) :: A
     real(kind=rkind), dimension(:), intent(in)  :: x
     real(kind=rkind), dimension(:), intent(out) :: y
     type(info_type), intent(inout) :: info
     integer(kind=ikind) :: i

      y = 0

      do i=1,A%nz
!       print *, A%vals(i), i
         y(A%ii(i))=y(A%ii(i)) + A%vals(i)*x(A%jj(i))
      end do

   end subroutine mulAx1

   subroutine mulAtx1(A,x,y,info)
     type(smtx), intent(in) :: A
     real(kind=rkind), dimension(:), intent(in)  :: x
     real(kind=rkind), dimension(:), intent(out) :: y
     type(info_type), intent(inout) :: info
     integer(kind=ikind) :: i

      y = 0

      do i=1,A%nz
         y(A%jj(i))=y(A%jj(i)) + A%vals(i)*x(A%ii(i))
      end do
   end subroutine mulAtx1

  subroutine transp(A)
    use typy
    type(smtx), intent(in out) :: A
    integer(kind=ikind) :: wrk,i

    do i=1,ubound(A%ii,1)
      wrk     = A%ii(i)
      A%ii(i) = A%jj(i)
      A%jj(i) = wrk
    end do
    call sort_matrix(A)
  end subroutine transp


  subroutine Dump_Matrix(A)
      use typy
      type(smtx), intent(in) :: A
      integer(kind=ikind) :: i

      print *, "uplny vypis struktury matice"
      print *, "pocet radku    nc=",A%row
      print *, "pocet sloupcu col=",A%col
      print *, "pocet nenul    nz=",A%nz
      print *, "alokovane rozsahy allocsize=",A%allocsize
      print *, "pocet seznamu volnych bloku efs=",A%efs
      print *, "datove udaje matice"
      print *, "index","  vals  ","  ii  ","  jj  "
      do i=1,a%allocsize
         print *,i,A%vals(i),A%ii(i),A%jj(i)
      end do
      print *, " ted jeste radky"
      do i=1,A%row
        print *, i, A%rows(i,1), A%rows(i,2)
      end do

  end subroutine Dump_Matrix

  !> vytvori matici diskretizace ulohy ...
  !!
  !! vytvori a naalokuje matici podle parametru
  !<
  subroutine Laplace2D(A,nx,ny)
    use typy
    !> konstruovana matice - je pripadne realokovana
    type(smtx), intent(in out) :: A
    !> pocet uzlu ve smeru x
    integer(kind=ikind), intent(in) :: nx
    !> pocet uzlu ve smeru y
    integer(kind=ikind), intent(in) :: ny

    !! lokalni promenne
    integer(kind=ikind) :: n, nz, i, j, i1

    print *,"Laplace2D zacina"

    n = nx*ny
    nz = 5*n ! zatim jen odhad

    !! zkontroluj alokacni status a smaz
    call clear_matrix(A)

    !! realokuj
    call init_matrix(A,n,n,nz)
    !! a nasyp
    ! realne nasypu vals, ii a jj
    nz = 0
    do i=1,nx
      do j=1,ny
        i1 = (i-1)*ny + j
        nz = nz + 1
        A%vals(nz) = 4
        A%ii(nz)   = i1
        A%jj(nz)   = i1
        if (i > 1) then
          nz = nz + 1
          A%vals(nz) = -1
          A%ii(nz)   = i1
          A%jj(nz)   = i1-ny
        end if
        if (i < nx) then
          nz = nz + 1
          A%vals(nz) = -1
          A%ii(nz)   = i1
          A%jj(nz)   = i1+ny
        end if
        if (j > 1) then
          nz = nz + 1
          A%vals(nz) = -1
          A%ii(nz)   = i1
          A%jj(nz)   = i1-1
        end if
        if (j < ny) then
          nz = nz + 1
          A%vals(nz) = -1
          A%ii(nz)   = i1
          A%jj(nz)   = i1+1
        end if

      end do
    end do
    A%nz = nz
    ! nakonec zavolam dorovnani
    call sort_matrix(A)
  end subroutine Laplace2D

  !> smaze matici A a uvolni jeji vsechna data
  recursive subroutine clear_matrix(A)
    use typy
    !> mazana matice
    type(smtx), intent(inout) :: A

    A%row = 0
    A%col = 0
    A%nz  = 0
    A%allocsize = 0
    A%efs      = 0
    if (allocated(A%vals)) then
      deallocate(A%vals)
    end if
    if (allocated(A%ii)) then
      deallocate(A%ii)
    end if
    if (allocated(A%diag)) then
      deallocate(A%diag)
    end if
    if (allocated(A%jj)) then
      deallocate(A%jj)
    end if
      if (allocated(A%rows)) then
      deallocate(A%rows)
    end if
    ! smazni pripadnou dekompozici
    if (associated(A%decomp)) then
      call clear_matrix(A%decomp)
      if (allocated(A%perm)) then
        deallocate(A%perm)
      end if
    end if
    ! ted jeste predpodmineni
    if(associated(A%prcd)) then
      print *, "zatim jen nuluji - to je chybne"
      deallocate(A%prcd)
    end if


  end subroutine clear_matrix






  subroutine GetSubmatrix(A,si,sj,SA)
    use typy
    type(smtx), intent(in) :: A
    integer(kind=ikind), dimension(:), intent(in) :: si
    integer(kind=ikind), dimension(:), intent(in) :: sj
    type(smtx), intent(inout) :: SA
    logical, dimension(1:A%row) ::mi
    logical, dimension(1:A%col) ::mj
    integer(kind=ikind), dimension(1:A%row) :: imi
    integer(kind=ikind), dimension(1:A%col) :: imj
    integer(kind=ikind) :: i,j,nz,ni,nj

    mi = .false.
    mj = .false.
    ! nefunguje kdyz si nebo sj mimo meze
    mi(si) = .true.
    mj(sj) = .true.
    imi = 0
    imj = 0
    do i=1,ubound(si,1)
      imi(si(i)) = i
    end do
    do i=1,ubound(sj,1)
      imj(sj(i)) = i
    end do
    !spocitame pocet nenul submatice
    nz = 0
    do i=1,A%row
      if (mi(i)) then
        !radek tam potencialne patri
        do j=A%rows(i,1),A%rows(i,2)
          if (mj(A%jj(j))) then
            nz = nz+1
          end if
        end do
      end if
    end do
    print *, " v submatici je ",nz," prvku"
!     if (nz == 0) then
!       print *, sj; pause
!       print *, si; pause
!       call dump_matrix(A); stop
!     end if
        
    call clear_matrix(SA)
    ni = ubound(si,1)
    nj = ubound(sj,1)
    call init_matrix(SA,ni,nj,nz)
    ! ted naplnime
    !spocitame pocet nenul submatice
    nz = 0
    do i=1,A%row
      if (mi(i)) then
        !radek tam potencialne patri
        do j=A%rows(i,1),A%rows(i,2)
          if (mj(A%jj(j))) then
            nz = nz+1
            SA%vals(nz) = A%vals(j)
            SA%ii(nz) = imi(A%ii(j))
            SA%jj(nz) = imj(A%jj(j))
          end if
        end do
      end if
    end do
    SA%nz = nz

    ! usporadame
   call sort_matrix(SA)
   !call print_matrix(SA)
  end subroutine GetSubmatrix







  subroutine Sparse_Test_Driver()
      use typy
      !implicit none
      integer(kind=ikind) :: choice, n, nx, ny, prekryv, i
      integer :: ierr
      character(len=200) :: nazev
      integer(kind=ikind),dimension(:), pointer :: deleni
      type(smtx) :: A
      logical :: volba


      do
          print *, " Testovani ridkych matic"
          print *, " Vyber cinnost"
          print *, "0 .... konec"
          print *, "1 .... cteni matice"
          print *, "2 .... konstrukce prazdne matice"
          print *, "3 .... tisk matice"
          print *, "4 .... dump matice"
          print *, "5 .... manipulace s matici"
          print *, "6 .... konstrukce 2D Laplace"
          read(unit=*,fmt=*) choice
          select case(choice)
              case (0)
                  exit
              case (1)
                  print *, "jdeme cist"
                  call read_mtx(A)
              case (2)
                  print *, " zadej rad matice"
                  read(unit=*,fmt=*) choice
                  call init_matrix(A,choice)
              case (3)
                  call print_matrix(A)
              case (4)
                  call Dump_Matrix(A)
              case (5)
                   print *, "Budeme upravovat matici"
              case (6)
                   print *, "pocet uzlu ve smeru x"
                   read(unit=*, fmt=*) nx
                   print *, "pocet uzlu ve smeru y"
                   read(unit=*, fmt=*) ny
                   call Laplace2D(A,nx,ny)
              case default
                  print *, "chybna volba"
          end select
      end do
  end subroutine Sparse_Test_Driver

  !> spocita residuum
  function resid(A,b,x) result(rs)
    implicit none
    type(smtx), intent(in) :: A
    real(kind=rkind), dimension(:), intent(in) :: b
    real(kind=rkind), dimension(:), intent(in) :: x
    real(kind=rkind) :: rs
    type(info_type) :: if
    real(kind=rkind), dimension(lbound(b,1):ubound(b,1))  :: wrk

    call mulAsx(A,x,wrk,if)
    wrk = b-wrk
    rs = sqrt(dot_product(wrk,wrk))
  end function resid

  subroutine matnorm(A,nf,n1,emin, emax)
    use typy
    implicit none
    type(smtx) , intent(in) :: A
    real(kind=rkind), intent(out) :: nf
    real(kind=rkind), intent(out) :: n1
    real(kind=rkind), intent(out) :: emin
    real(kind=rkind), intent(out) :: emax
    integer(kind=ikind) :: i,j
    real(kind=rkind) :: w1, w2

    nf = 0
    emax = 0
    emin = abs(A%vals(A%rows(1,1)))
    n1 = emin
    do i=1,A%row
      w2 = 0
      do j=A%rows(i,1), A%rows(i,2)
        w1 = abs(A%vals(j))
        nf = nf + w1*w1
        w2 = w2 + w1
        if (w1 > emax ) emax = w1
        if (w1 < emin ) emin = w1
      end do
      if (w2 > n1 ) n1 = w2
    end do
    nf = sqrt(nf)
  end subroutine matnorm

end module mytypes
