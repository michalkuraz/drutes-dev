!> \file sparsematrix.f90
!! \brief implementace vlastnosti ridke matice
!!
!! xxx
!<

!> Ridka matice
!! --
!! \page ridka Ridka matice
!<
module sparsematrix
    use typy
    use mtx
    implicit none


    !> typ pro prvek ridke matice
    type, public :: elemt
        !>hodnota
        real(kind=rkind) :: vl = 0.0_rkind
        !> radek
        integer(kind=ikind) :: row = 0_ikind
        !> sloupec
        integer(kind=ikind) :: col = 0_ikind
        !> index dalsiho prvku v radku pokud zadny neni = -1
        integer(kind=ikind) :: nextinrow = -1_ikind
        !> index dalsiho prvku ve sloupci pokud zadny neni = -1
        integer(kind=ikind) :: nextincol = -1_ikind
    end type elemt

    !> ridka matice
    type, public, extends(matrix) :: smtx
        !> pocet nenul
        integer(kind=ikind), private :: nz = 0
        !> pocet alokovanych prvku
        integer(kind=ikind), private :: aloc = 0
        !> hodnoty prvku
        type(elemt), dimension(:), allocatable, private :: v
        !> zacatky radku - pokud prazdny = -1
        integer(kind=ikind), dimension(:), allocatable, private :: rowstart
        !> zacatky sloupcu - pokud prazdny = -1
        integer(kind=ikind), dimension(:), allocatable, private :: colstart
        !> index prvniho volneho prvku, pokud neni = -1
        integer(kind=ikind), private :: firstfree
    contains
        procedure :: init => initsp
        procedure :: get => getsp
        procedure, pass(a) :: set => setsp
        procedure :: print => sparseprint
        procedure :: dump => sparsedump
        procedure :: nonzero => sparsenz
        procedure :: getcol => getcolsparse
        procedure :: getrow => getrowsparse
        procedure, private :: findelem => felem
        procedure, private :: insert => insertsp
        procedure, private :: remove => removesp
    end type smtx

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! veci pro ridkou matici
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> najde prvek v ridke strukture
    function felem(a,i,j) result(res)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in) :: a
        !> souradnice
        integer(kind=ikind), intent(in) :: i
        !> souradnice
        integer(kind=ikind), intent(in) :: j
        integer(kind=ikind) :: res
        integer(kind=ikind) :: k

        res = -1
        if (a%nz == 0) return
        k = a%rowstart(i)
        do while (k > 0)
            if (a%v(k)%col == j) then
                res = k
                return
            end if
            k = a%v(k)%nextinrow
        end do
    end function felem


    !> ziska prvek z ridke matice
    function getsp(a,i,j) result (r)
        use typy
        !> matice
        class(smtx), intent(in) :: a
        !> souradnice
        integer(kind=ikind), intent(in) :: i,j
        integer(kind=ikind) :: index
        real(kind=rkind) :: r
        r = 0
        index = a%findelem(i,j)
        if (index > 0) r = a%v(index)%vl
    end function

    !> nastavi prvek matice, pokud neni tak vlozi
    subroutine setsp(r,a,i,j)
        implicit none
        !> hodnota
        real(kind=rkind), intent(in) :: r
        !> matice
        class(smtx), intent(in out)  :: a
        !> souradnice
        integer(kind=ikind), intent(in) :: i,j
        integer(kind=ikind) :: index

        index = a%findelem(i,j)
        if ( r == 0) then
            ! jde vlastne o mazani
            if ( index > 0 ) then
                call a%remove(i,j,index)
            end if
        else
            if ( index > 0 ) then
                !print *, "prvek uz existuje"
                a%v(index)%vl = r
            else
                !print *, "prvek neexistuje" 
                call a%insert(r,i,j)
            end if
        end if
        !call a%dump
        !pause
    end subroutine setsp

    !> prida neexistujici prvek do matice
    subroutine insertsp(a,r,i,j)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in out) :: a
        !> hodnota
        real(kind=rkind), intent(in) :: r
        !> souradnice
        integer(kind=ikind), intent(in) :: i,j

        type(elemt), dimension(:), allocatable :: temp
        integer(kind=ikind) :: k, w1

        if (a%aloc == 0 ) then
            ! prvotni alokace
            a%aloc=2*a%getn()
            allocate(a%v(1:a%aloc))
            allocate(a%colstart(1:a%getm()))
            allocate(a%rowstart(1:a%getn()))
            a%nz = 0
            do k=1,a%aloc
                a%v(k) = elemt(0.0,0,0,k+1,k+1)
            end do
            a%v(a%aloc) = elemt(0.0,0,0,-1,-1)
            a%colstart = -1
            a%rowstart = -1
            a%firstfree = 1
        end if
        ! neco uz je
        if ( a%firstfree == -1) then
            ! neni misto, musim ho udelat
            if (a%nz /= a%aloc) stop "insert: nekonzistentni data v matici"
            a%aloc = max(a%aloc+10,a%aloc + a%aloc/10)
            allocate(temp(1:a%aloc))
            temp(1:ubound(a%v,1)) = a%v
            call move_alloc(temp,a%v)
            ! dodat to do volnyho mista
            do k=a%nz+1,a%aloc
                a%v(k) = elemt(0.0,0,0,k+1,k+1)
            end do
            a%v(a%aloc) = elemt(0.0,0,0,-1,-1)
            a%firstfree = a%nz+1
        end if
        w1 = a%firstfree                 ! vyberu to z volnyho mista
        a%firstfree = a%v(w1)%nextincol  ! a oznacim novy volny misto
        a%v(w1) = elemt(r,i,j,a%rowstart(i),a%colstart(j))
        a%rowstart(i) = w1
        a%colstart(j) = w1
        a%nz = a%nz + 1
    end subroutine insertsp

    !> smazne existujici prvek z ridke matici
    subroutine removesp(a,i,j,index)
        use typy
        implicit none
        !> matice
        class(smtx), intent(inout) :: a
        !> souradnice mazaneho prvku
        integer(kind=ikind), intent(in) :: i
        !> souradnice mazaneho prvku
        integer(kind=ikind), intent(in) :: j
        !> pozice mazaneho prvku
        integer(kind=ikind), intent(in) :: index

        integer(kind=ikind) :: k

        !kontrola, ze to sedi
        if (a%v(index)%row /= i .or. a%v(index)%col /= j) &
            stop "remove: nekonzistentni data"
        !vlastni smazani
        ! napred vynechat mazany prvek ze seznamu
        ! napred v radku
        !call pockej( "jsem v mazani")
        k = a%rowstart(i)
        if ( k == index) then
            a%rowstart(i) = a%v(index)%nextinrow
        else
            do while (a%v(k)%nextinrow /= index)
                k = a%v(k)%nextinrow
            end do
            a%v(k)%nextinrow = a%v(index)%nextinrow
        end if
        ! potom sloupec
        k = a%colstart(j)
        if ( k == index) then
            a%colstart(j) = a%v(index)%nextincol
        else
            do while (a%v(k)%nextincol /= index)
                k = a%v(k)%nextincol
            end do
            a%v(k)%nextincol = a%v(index)%nextincol
        end if
        ! pridame do seznamu volneho mista
        a%v(index) = elemt(0.0,0,0,a%firstfree,a%firstfree)
        a%firstfree = index
        a%nz = a%nz - 1
    end subroutine removesp

    !> inicializace a vymazani matice
    subroutine initsp(a,n,m)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in out)  :: a
        !> rozmery
        integer(kind=ikind), intent(in) :: n,m
        call a%resize(n,m)
        a%nz = 0
        a%aloc = 0
        if (allocated(a%v)) deallocate(a%v)
        if (allocated(a%rowstart)) deallocate(a%rowstart)
        if (allocated(a%colstart)) deallocate(a%colstart)
    end subroutine initsp



    !> vytiskne ridkou matici
    subroutine sparseprint(a,ncol, width, caption)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in) :: a
        !> pocet sloupci tisku
        integer(kind=ikind), intent(in), optional :: ncol
        !> sirka cisla
        integer(kind=ikind), intent(in), optional :: width
        !> nadpis
        character(len=*), intent(in),    optional :: caption
        character(len=100) :: cpt,fmts
        integer(kind=ikind) :: nc
        integer(kind=ikind) :: wd
        integer(kind=ikind) :: i,j
        real(kind=rkind) :: wrk
        integer(kind=ikind), dimension(:), allocatable :: jj
        real(kind=rkind), dimension(:), allocatable :: v
        integer(kind=ikind) :: nelem

        if ( present(ncol) ) then
            nc = ncol
        else
            nc = 5
        end if

        if ( present(width) ) then
            wd = width
        else
            wd = 15
        end if

        if ( present(caption) ) then
            cpt = adjustl(caption)
        else
            cpt ="ridka matice"
        end if

        print *, trim(cpt)
        if (a%aloc == 0) then
            print *,"prazdna matice"
        else
            !print *, "printsmtx 1"
            do i=1,a%getn()
                call a%getrow(i,v,jj,nelem)
                !print *,"printsmtx 2"
                if ( nelem == 0) then
                    print *, "radek ",i," je prazdny"
                else
                    print *, "radek ",i," ma ",nelem," prvku"
                    do j = 1,nelem
                        print *,v(j),jj(j)
                    end do
                end if
            end do
        end if
        print *, " pocet radku=",a%getn(), " pocet sloupcu=", a%getm(), &
            " pocet nenulovych prvku=", a%nz

    end subroutine sparseprint


    !> ridky dump
    subroutine sparsedump(a)
        implicit none
        class(smtx), intent(in) :: a
        integer(kind=ikind) :: i

        print *,"ridky dump"
        print *, "pocet radku:", a%getn(), " pocet sloupcu:",a%getm(),&
            " alokovano:",a%aloc," pocet nenul:", a%nz
        if (a%nz > 0) then
            print *,"elementy matice:"
            do i=1,a%aloc
                print *, "index=",i
                print *, "value:",a%v(i)
            end do
            print *, "zacatky radku:"
            do i=1,a%getn()
                print *, i, a%rowstart(i)
            end do
            print *, "zacatky sloupcu"
            do i=1,a%getm()
                print *, i, a%colstart(i)
            end do
        end if
    end subroutine sparsedump

    !> vrati pocet nenul v matici
    function sparsenz(a) result(res)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in) :: a
        integer(kind=ikind) :: res
        res = a%nz
    end function sparsenz



    subroutine getcolsparse(a,j,v,ii,nelem,mask)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in) :: a
        !> identifikace sloupce
        integer(kind=ikind), intent(in) :: j
        !> hodnoty prvku
        real(kind=rkind), dimension(:), allocatable, intent(inout) :: v
        !> radkove indexy
        integer(kind=ikind), dimension(:), allocatable, intent(inout) :: ii
        !> pocet prvku vybraneho sloupce
        integer(kind=ikind), intent(out) :: nelem
        !> seznam povolenych indexu
        logical, dimension(:), optional, intent(in) :: mask
        integer(kind=ikind) :: i,ic
        real(kind=rkind) :: w
        integer(kind=ikind), dimension(:), allocatable :: iiw

        !! spocitame delku
        !print *, "jsem v getcol sparse"
        i = a%colstart(j)
        if (i== -1) then
            !! je to prazdny sloupec
            nelem = 0
            return
        else
            nelem = 0
            do
                ic = i
                i = a%v(i)%nextincol
                if ( present(mask) ) then
                    if ( mask(a%v(ic)%row) ) then
                        nelem=nelem+1
                        !print *,ic, a%v(ic)%row,mask(a%v(ic)%row),i
                    end if
                else
                    nelem = nelem + 1
                end if
                if (i == -1) exit
            end do
        end if
        !print *,"nelem=",nelem
        !! v nelem je pocet prvku sloupce
        !! prekontrolovat a pripadne pripravit prostor
        if (allocated(v) ) then
            if ( ubound(v,1) < nelem ) then
                deallocate(v)
                allocate(v(1:nelem))
            end if
        else
            allocate(v(1:nelem))
        end if
        !print *,"v hotovo"
        if (allocated(ii) ) then
            if ( ubound(ii,1) < nelem ) then
                deallocate(ii)
                allocate(ii(1:nelem))
            end if
        else
            !print *,"nekompromisni alokace ii"
            !print *,"nelem=",nelem
            allocate(ii(1:nelem))
            !print *, "hotovo"
            !!ii = iiw
        end if
        !print *, "realokovano"
        i = a%colstart(j)
        ic = 0
        !print *,"i=",i," ic=", ic
        do
            if ( present(mask)) then
                if (mask(a%v(i)%row)) then
                    !print *, "pridavam"
                    !print *,i, a%v(i)%row,mask(a%v(i)%row),&
                    !    a%v(i)%nextincol
                    ic = ic + 1
                    v(ic)  = a%v(i)%vl
                    ii(ic) = a%v(i)%row
                end if
            else
                ic = ic + 1
                v(ic)  = a%v(i)%vl
                ii(ic) = a%v(i)%row
            end if
            i = a%v(i)%nextincol
            !print *,"i=",i," ic=", ic
            if ( i==-1) exit
        end do
        if (ic /= nelem) then
            print *, "scol:konecne ic ",ic,"nelem=",nelem
            STOP
        end if
    end subroutine getcolsparse

    subroutine getrowsparse(a,i,v,jj,nelem,mask)
        use typy
        implicit none
        !> matice
        class(smtx), intent(in) :: a
        !> identifikace sloupce
        integer(kind=ikind), intent(in) :: i
        !> hodnoty prvku
        real(kind=rkind), dimension(:), allocatable, intent(inout) :: v
        !> radkove indexy
        integer(kind=ikind), dimension(:), allocatable, intent(inout) :: jj
        !> pocet prvku vybraneho sloupce
        integer(kind=ikind), intent(out) :: nelem
        !> seznam povolenych indexu
        logical, dimension(:), optional, intent(in) :: mask
        integer(kind=ikind) :: j,ic
        real(kind=rkind) :: w

        !! spocitame delku
        !!print *, "jsem v getrow sparse"
        j = a%rowstart(i)
        if (j == -1) then
            !! je to prazdny sloupec
            nelem = 0
            return
        else
            nelem = 0
            do
                ic = j
                j = a%v(j)%nextinrow
                if ( present(mask) ) then
                    if ( mask(a%v(ic)%col) ) nelem=nelem+1
                else
                    nelem = nelem + 1
                end if
                if (j == -1) exit
            end do
        end if
        !!print *,"nelem=",nelem
        !! v nelem je pocet prvku sloupce
        !! prekontrolovat a pripadne pripravit prostor
        if (allocated(v) ) then
            !!print *,"kontroluji v"
            if ( ubound(v,1) < nelem ) then
                deallocate(v)
                allocate(v(1:nelem))
            end if
        else
            !!print *, "nekompromisni alokace v"
            allocate(v(1:nelem))
        end if
        !!print *,"v hotovo"
        if (allocated(jj) ) then
            if ( ubound(jj,1) < nelem ) then
                deallocate(jj)
                allocate(jj(1:nelem))
            end if
        else
            allocate(jj(1:nelem))
        end if
        !!print *, "realokovano"
        j = a%rowstart(i)
        ic = 0
        !!print *,"j=",j," ic=", ic
        do
            if ( present(mask)) then
                if (mask(a%v(j)%col)) then
                    ic = ic + 1
                    v(ic)  = a%v(j)%vl
                    jj(ic) = a%v(j)%col
                end if
            else
                ic = ic + 1
                v(ic)  = a%v(j)%vl
                jj(ic) = a%v(j)%col
            end if
            j      = a%v(j)%nextinrow
            !!print *,"j=",j," ic=", ic
            if ( j==-1) exit
        end do
        !!print *," konecne ic=",ic
        !!print *,v,jj
    end subroutine getrowsparse


end module sparsematrix
