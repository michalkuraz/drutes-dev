!> \file mtx_int.f90
!! \brief Objektova implementace matic
!<


!> Module pro matice
!! ==
!! \page implem Implementace matic
!! Jde o implementaci abstraktniho typu. Zde je jednak nekolik malo
!! procedur, ktere je nutno pri konstrukci potomka implementovat. Dale
!! je zde rada nastroju, ktera je pomoci techto prostredku
!! implementovana, ale lze je v pripade vhodnosti reimplementovat v
!! potomkovi.
!! Zmeny v celkovem interface se nepredpokladaji.
!! Typ ma dve interni polozky - rozmery matice. Jsou zpristupneny
!! pomoci pristupovych funkci.
!! - \subpage plna "Plna matice"
!! - \subpage ridka "Ridka matice"
module mtx_int
    use typy
    implicit none
    private

    !> obecna matice
    type, abstract, public :: matrix_int
    !> pocet radku
    integer(kind=ikind), private :: n = 0
    !> pocet sloupcu
    integer(kind=ikind), private :: m = 0
    !> podrobnost informaci o vypoctech
    !! - 0 ... pracuj tise
    !! - 1 ... informuj jen o velmi dulezitych vecech
    !! - 10... maximalni ukecanost
    integer, private :: verbosity = 0
contains
    !> vrati pocet radku
    procedure, non_overridable :: getn => getnmatrix
    !> vrati pocet sloupcu
    procedure, non_overridable :: getm => getmmatrix
    !> nastavi podrobnost informaci
    procedure, non_overridable :: setverbosity
    !> inicializace
    procedure :: init => initmatrix
    !> vrati prvek matice
    procedure(geti), deferred :: get
    !> nastavi prvek matice
    procedure(seti), pass(a), deferred :: set
    !> pricte cislo k prvku matice
    procedure :: add => addmatrix
    !> vybere radek z matice, vraceny jsou nenulove prvky a jejich
    !! souadnice
    procedure :: getrow => getrowmatrix
    !> vybere sloupec z matice, vraceny jsou nenulove prvky a jejich
    !! souadnice
    procedure :: getcol => getcolmatrix
    !> vybere submatici z matice
    procedure :: submatrix => smatrix
    !> vytiskne matici
    procedure :: print => printmatrix
    !> dumpne celou strukturu matice
    procedure :: dump => dumpmatrix
    !> da pocet nenul v matici
    procedure :: nonzero => nzmatrix
    !> textova varianta spy
    procedure, non_overridable :: spy => spymatrix
    !> nasobi vektorem, vrati soucin matice krat vektor
    procedure :: mul => mulmatrix
    !> nasobi vektorem, vrati soucin transponovane matice krat vektor
    procedure :: mult => mulmatrixt
    !> znasobi dve matice
    procedure :: mulm => mulmatrixmatrix
    !> odecte argument od matice
    procedure :: subm => submatrixmatrix
    !> pricte argument k matici
    procedure :: addm => addmatrixmatrix
    !> udela kopii matice
    procedure :: copy => copymatrix
    !> udela prazdnou kopii matice
    procedure :: clone => clonematrix
    !> jen prenastavi velikost
    procedure :: resize => initmatrix
    !> spocita Frobeniovskou normu matice
    procedure :: normF => normFmatrix
    !> definovane prirazeni
    !procedure, private :: copy1
    generic :: assignment(=) => copy
    !> transpozice
    procedure :: transpose
    !> cteni ze souboru
    procedure :: Read
    !> zapis do souboru
    procedure :: Write
end type matrix_int


!> interface pro ziskani prvku
abstract interface
    !> funkce pro ziskani prvku
    function geti(a,i,j) result(r)
        use typy
        import :: matrix_int
        !> matice
        class(matrix_int), intent(in) :: a
        !> radkovy index
        integer(kind=ikind), intent(in) :: i
        !> sloupcovy index
        integer(kind=ikind), intent(in) :: j
        !> hodnota prvku
        integer(kind=ikind) :: r
    end function
end interface

!> interface pro nastaveni prvku
abstract interface
    !> rutina pro nastaveni prvku
    subroutine seti(r,a,i,j)
        use typy
        import :: matrix_int
        !> matice
        class(matrix_int), intent(in out) :: a
        !> souradnice
        integer(kind=ikind), intent(in) :: i,j
        !> hodnota
        integer(kind=ikind), intent(in) :: r
    end subroutine
end interface

public :: mtxtest_int


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! veci pro matrix - zcela univerzalni
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setverbosity(A,level)
    class(matrix_int), intent(in out) :: A
    integer, intent(in) :: level
    A%verbosity = max(0,min(level,10))
end subroutine setverbosity




!> pricte k prvku - kandidat reimplementace v potomkovi
subroutine addmatrix(a,r,i,j)
    use typy
    !> matice
    class(matrix_int), intent(in out) :: a
    !> souradnice
    integer(kind=ikind), intent(in) :: i,j
    !> hodnota
    integer(kind=ikind), intent(in) :: r
    integer(kind=ikind) :: v

    v=a%get(i,j)
    v = v + r
    call a%set(v,i,j)
end subroutine addmatrix

!> rutina pro inicializaci
subroutine initmatrix(a,n,m)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in out) :: a
    !> rozmery
    integer(kind=ikind), intent(in) :: n,m
    a%n = n
    a%m = m
    !print *,"initmatrix: konci"
end subroutine



!> vrati pocet radku
pure function getnmatrix(a) result(res)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    integer(kind=ikind) :: res
    res = a%n
end function getnmatrix


!> vrati pocet sloupcu
pure function getmmatrix(a) result(res)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    integer(kind=ikind) :: res
    res = a%m
end function getmmatrix

!> ziska radek z matice - univerzal - kandidat reimplementace v potomkovi
!!
!! vrati jen nenulove prvky z radku, v pripade potreby prealokuje prostor
!<
subroutine getrowmatrix(a,i,v,jj,nelem,mask)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    !> pozadovany radek
    integer(kind=ikind), intent(in) :: i
    !> nenulove hodnoty z radku
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: v
    !> sloupcove indexy nenulovych prvku radku
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: jj
    !> pocet nenulovych prvku v radku
    integer(kind=ikind), intent(out) :: nelem
    !> seznam povolenych indexu
    logical, dimension(:), optional, intent(in) :: mask

    integer(kind=ikind) :: j
    integer(kind=ikind) :: w

    ! dvoufazove
    ! spocitat kolik jich bude
    nelem = 0
    if (present(mask)) then
        do j = 1,a%m
            if ( a%get(i,j) /= 0 .and. mask(j)) then
                nelem = nelem + 1
            end if
        end do
    else
        do j = 1,a%m
            if ( a%get(i,j) /= 0 ) then
                nelem = nelem + 1
            end if
        end do
    end if
    ! overit misto
    if ( .not. allocated(v)) then
        allocate(v(1:nelem))
    else
        if (nelem > ubound(v,1)) then
            deallocate(v)
            allocate(v(1:nelem))
        end if
    end if
    if ( .not. allocated(jj)) then
        allocate(jj(1:nelem))
    else
        if (nelem > ubound(jj,1)) then
            deallocate(jj)
            allocate(jj(1:nelem))
        end if
    end if
    ! opravdu naplnit
    nelem = 0
    if (present(mask)) then
        do j = 1,a%m
            w = a%get(i,j)
            if ( w /= 0 .and. mask(j)) then
                nelem = nelem + 1
                v(nelem)  = w
                jj(nelem) = j
            end if
        end do
    else
        do j = 1,a%m
            w = a%get(i,j)
            if ( w /= 0) then
                nelem = nelem + 1
                v(nelem)  = w
                jj(nelem) = j
            end if
        end do
    end if
end subroutine getrowmatrix

!> ziska sloupec z matice  - kandidat reimplementace v potomkovi
subroutine getcolmatrix(a,j,v,ii,nelem,mask)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    !> identifikace sloupce
    integer(kind=ikind), intent(in) :: j
    !> hodnoty prvku
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: v
    !> radkove indexy
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: ii
    !> pocet prvku vybraneho sloupce
    integer(kind=ikind), intent(out) :: nelem
    !> seznam povolenych indexu
    logical, dimension(:), optional, intent(in) :: mask
    integer(kind=ikind) :: i
    integer(kind=ikind) :: w
    ! dvoufazove
    ! spocitat kolik jich bude
    nelem = 0
    if ( present(mask)) then
        do i = 1,a%n
            if ( a%get(i,j) /= 0 .and. mask(i)) then
                nelem = nelem + 1
            end if
        end do
    else
        do i = 1,a%n
            if ( a%get(i,j) /= 0) then
                nelem = nelem + 1
            end if
        end do
    end if
    ! overit misto
    if ( .not. allocated(v)) then
        allocate(v(1:nelem))
    else
        if (nelem > ubound(v,1)) then
            deallocate(v)
            allocate(v(1:nelem))
        end if
    end if
    if ( .not. allocated(ii)) then
        allocate(ii(1:nelem))
    else
        if (nelem > ubound(ii,1)) then
            deallocate(ii)
            allocate(ii(1:nelem))
        end if
    end if
    ! opravdu naplnit
    nelem = 0
    if (present(mask)) then
        do i = 1,a%n
            w = a%get(i,j)
            if ( w /= 0 .and. mask(i)) then
                nelem = nelem + 1
                v(nelem)  = w
                ii(nelem) = i
            end if
        end do
    else
        do i = 1,a%n
            w = a%get(i,j)
            if ( w /= 0 ) then
                nelem = nelem + 1
                v(nelem)  = w
                ii(nelem) = i
            end if
        end do
    end if
end subroutine getcolmatrix


!> napodoba matlabovske spy
subroutine spymatrix(a)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    integer(kind=ikind) :: i,j, nz
    character(len=a%getm()) :: radek
    character(len=8) :: fmts

    nz = 0
    write(fmts,fmt="(a2,i5,a1)") "(a",a%m,")"
    do i=1,a%getn()
        do j=1,a%getm()
            if (a%get(i,j) == 0) then
                radek(j:j) = '.'
            else
                radek(j:j) = 'X'
                nz = nz + 1
            end if
        end do
        print fmts,radek
    end do
    print *, " pocet radek=", A%getn()
    print *, " pocet sloupcu",A%getm()
    print *,"celkovy pocet nenul=",nz
end subroutine spymatrix

!> vytiskne matici
subroutine printmatrix(a,ncol, width, caption)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    !> pozadovany pocet sloupcu tisku, nepovinne
    integer(kind=ikind), intent(in), optional :: ncol
    !> pozadovane sirka sloupce, nepovinne
    integer(kind=ikind), intent(in), optional :: width
    !> nadpis, nepovinne
    character(len=*), intent(in),    optional :: caption
    character(len=100) :: cpt,fmts
    integer(kind=ikind) :: nc
    integer(kind=ikind) :: wd
    integer(kind=ikind) :: i,j, col1,col2
    real(kind=rkind), dimension(:), allocatable :: data

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
        cpt ="matice"
    end if

    allocate(data(1:nc))
    print *, trim(cpt)
    print *, " pocet radku=",a%n, " pocet sloupcu=", a%m
    !write(fmts,fmt="(i3,a2,i3,a2,i3,a2)") ncol,"(es",width,".",width-8,",a)"
    write(unit=fmts,fmt="(A1,I0,A2,I0,A1,I0,A1)")   "(",nc,"Es",wd+8,".",wd,")"
    !print *,fmts
    col1 = 1
    col2 = min(nc,a%getm())
    do
        print *,"tisknu sloupce ",col1," az ",col2
        do i=1,a%getn()
            do j=col1,col2
                data(j-col1+1) = a%get(i,j)
            end do
            print fmts, data(1:col2-col1+1)
        end do
        col1 = col2+1
        col2 = min(col2+nc,a%getm())
        if (col1>a%getm()) exit
    end do

end subroutine printmatrix

!> zakladni dump
subroutine dumpmatrix(a)
    implicit none
    !> dumpovana matice
    class(matrix_int), intent(in) :: a
    print *, "dump matice pro matrix - v zasade k nicemu"
    print *, "n=",a%n," m=",a%m
end subroutine dumpmatrix

!> vrati pocet nenul v matici - tady n*m
function nzmatrix(a) result(res)
    use typy
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    integer(kind=ikind) :: res
    res = a%n*a%m
end function nzmatrix


!> nasobi matici vektorem
function mulmatrix(a,x, count) result(y)
    !> matice
    class(matrix_int), intent(in) :: a
    !> vektor
    integer(kind=ikind), dimension(:), intent(in) :: x
    !> vysledek
    integer(kind=ikind), dimension(:), allocatable :: y
    !> pocet operaci
    type(tcount), intent(inout), optional :: count
    integer(kind=ikind) :: i,j,nelem
    integer(kind=ikind), dimension(:), allocatable :: v
    integer(kind=ikind), dimension(:), allocatable :: ii


    allocate(y(1:a%n))
    y = 0
    do i=1,a%n
        call a%getrow(i,v,ii,nelem)
        do j=1,nelem
            y(i) = y(i) + v(j)*x(ii(j))
        end do
        !do j=1,a%m
        !    y(i) = y(i) + a%get(i,j)*x(j)
        !end do
        if (present(count)) then
            count%ad  = count%ad  + nelem
            count%mul = count%mul + nelem
        end if
    end do
end function mulmatrix

!> nasobi transponovanou matici vektorem
function mulmatrixt(a,x, count) result(y)
    implicit none
    !> matice
    class(matrix_int), intent(in) :: a
    !> vektor
    integer(kind=ikind), dimension(:), intent(in) :: x
    !> vysledek
    integer(kind=ikind), dimension(:), allocatable :: y
    !> pocet operaci
    type(tcount), intent(inout), optional :: count
    integer(kind=ikind) :: i,j,nelem
    integer(kind=ikind), dimension(:), allocatable :: v
    integer(kind=ikind), dimension(:), allocatable :: ii

    allocate(y(1:a%m))
    y = 0
    do i=1,a%n
        call a%getrow(i,v,ii,nelem)
        do j = 1, nelem
            y(ii(j)) = y(ii(j)) + v(j)*x(i)
        end do
        !do j=1,a%m
        !    y(j) = y(j) + a%get(i,j)*x(i)
        !end do
        if (present(count)) then
            count%ad  = count%ad  + nelem
            count%mul = count%mul + nelem
        end if
    end do
end function mulmatrixt




!> zkopiruje matici  - kandidat reimplementace v potomkovi
subroutine copymatrix(A,B)
    !use pmatools
    implicit none
    !> cil
    class(matrix_int), intent(inout) :: A
    !> zdroj
    class(matrix_int), intent(in) :: B

    integer(kind=ikind) :: i,j
    integer(kind=ikind), dimension(:), allocatable :: v
    integer(kind=ikind), dimension(:), allocatable :: ii
    integer(kind=ikind) :: nelem

    call A%init(B%getn(),B%getm())
    do i=1,B%getn()
        call B%getrow(i,v,ii,nelem)
        do j=1,nelem
            call A%set(v(j),i,ii(j))
        end do
        !do j=1,B%getm()
        !    wrk = B%get(i,j)
        !    call A%set(wrk,i,j)
        !end do
    end do
end subroutine copymatrix

!> vytvori prazdnou kopii
subroutine clonematrix(A,B)
    !> cil
    class(matrix_int), intent(inout) :: A
    !> zdroj
    class(matrix_int), intent(in) :: B
    call A%init(B%getn(),B%getm())
end subroutine clonematrix

!> \brief  A = B*C  !!!!!prepsat
!!
!! \param A vysledna matice
!! \param B levy operand
!! \param C pravy operand
!<
subroutine mulmatrixmatrix(A,B,C,info)
    use typy
    implicit none
    class(matrix_int), intent(in out) :: A
    class(matrix_int), intent(in)     :: B
    class(matrix_int), intent(in)     :: C
    type(tcount),  intent(in out), optional :: info

    integer(kind=ikind) :: n,m,p,i,j,k
    integer(kind=ikind) :: wrk
    integer(kind=ikind), dimension(:), allocatable :: rd, sl
    integer(kind=likind), dimension(:), allocatable :: rdi, sli
    integer(kind=ikind) :: nrd, nsl
    type(tcount) :: inf1
    real :: t1,t2

    n = B%getn()
    m = C%getm()
    p = B%getm()
    i = C%getn()
    j = A%getn()
    k = A%getm()
    call cpu_time(t1)
    call A%init(j,k)
    if ((p/=i) .or. (j/=n) .or. (m/=k)) print *,"nesedi rozmery"
!    do i=1,n
!        do j=1,m
!            call A%set(0.0_rkind,i,j)
!            wrk = 0
!            do k=1,p
!                w1 = B%get(i,k)
!                w2 = C%get(k,j)
!                wrk = wrk + w1*w2
!            end do
!            call A%set(wrk,i,j)
!        end do
!    end do
    do k=1,p
        call B%GetCol(k,sl,sli,nsl)
        call C%GetRow(k,rd,rdi,nrd)
        !! pricist opravu
        do i=1,nrd
            do j=1,nsl
                wrk = sl(j)*rd(i)
                call A%add(wrk, sli(j), rdi(i))
            end do
        end do
        call cpu_time(t2)
        inf1%ad = inf1%ad + nrd*nsl
        inf1%mul = inf1%mul + nrd*nsl
        inf1%time = inf1%time + t2 - t1
    end do
    if (present(info)) call update_info(info,inf1)
end subroutine mulmatrixmatrix

subroutine addmatrixmatrix(A,B)
    use typy
    implicit none
    class(matrix_int), intent(in out) :: A
    class(matrix_int), intent(in)      :: B
    integer(kind=ikind)            :: i,j,n,m
    integer(kind=ikind)               :: w1,w2
    integer(kind=ikind), dimension(:), allocatable :: rb
    integer(kind=ikind), dimension(:), allocatable :: rib
    integer(kind=ikind) :: nrb
    logical, dimension(:), allocatable :: mask

    n = A%getn()
    m = A%getm()
    allocate(mask(1:m))
    mask = .true.
    do i=1,n
        call B%getrow(i,rb,rib,nrb,mask)
        do j=1,nrb
            call A%add(rb(j),i,rib(j))
        end do
    end do

end subroutine addmatrixmatrix

!> A = A - B
subroutine submatrixmatrix(A, B)
    use typy
    implicit none
    class(matrix_int), intent(in out) :: A
    class(matrix_int), intent(in) :: B
    integer(kind=ikind)            :: i,j,n,m
    integer(kind=ikind), dimension(:), allocatable :: rb
    integer(kind=ikind), dimension(:), allocatable :: rib
    integer(kind=ikind) :: nrb
    logical, dimension(:), allocatable :: mask

    n = A%getn()
    m = A%getm()
    do i=1,n
        call B%getrow(i,rb,rib,nrb)
        do j=1,nrb
            call A%add(-rb(j),i,rib(j))
        end do
    end do
end subroutine submatrixmatrix

!> \brief  spocte Frobeniovu normu matice
!!
!! \param A zkoumana matice
!! \return  hledana norma
!!
function normFmatrix(A) result(y)
    use typy
    implicit none
    class(matrix_int), intent(in) :: A
    real(kind=rkind) :: y
    integer(kind=ikind), dimension(:), allocatable :: rb
    integer(kind=ikind), dimension(:), allocatable :: rib
    integer(kind=ikind) :: nrb,i,j

    y = 0
    print *,"normF zacina"
    do i=1,A%getn()
        call A%getrow(i,rb,rib,nrb)
        do j=1,nrb
            y = y + rb(j)*rb(j)
        end do
    end do
    y = sqrt(y)
    print *,"normF konci"
end function normFmatrix

function smatrix(A,ii,jj) result(B)
    !! jen na chvili, kvuli pockej
    !use pmatools
    implicit none
    class(matrix_int), intent(in) :: A
    class(matrix_int), allocatable :: B
    integer(kind=ikind), dimension(:), intent(in) :: ii
    integer(kind=ikind), dimension(:), intent(in) :: jj
    integer(kind=ikind) :: n,m,i,j
    integer(kind=ikind) :: wrk
    integer(kind=ikind), dimension(:), allocatable :: v
    integer(kind=ikind), dimension(:), allocatable :: is, invi
    integer(kind=ikind) :: nelem
    logical, dimension(:), allocatable :: mask

    print *, "zacina konstrukce submatice"
    n = size(ii,1)
    m = size(jj,1)
    !call A%print
    !print *,"pred alokaci"
    allocate(B,source=A)
    !print *, "alokovano"
    call B%init(n,m)
    !call B%print(caption="matice B")
    !call A%print(caption="matice A")
    allocate(mask(1:A%getm()))
    allocate(invi(1:A%getm()))
    !print *,"maska vytvorena"
    mask = .false.
    mask(jj) = .true.
    invi = 0
    do i = 1,m
        invi(jj(i)) = i
    end do
    !print *,"maska vytvorena a naplnena"
    do i=1,n
        !print *,"i=",i
        call A%GetRow(ii(i),v,is,nelem)
        do j = 1,nelem
            if ( mask(is(j)) ) then
                ! index tam patri
                call B%set(v(j),i,invi(is(j)))
            end if
        end do
        !!!! dokoncit
        !do j=1,m
            !wrk = A%get(ii(i),jj(j))
            !if ( wrk /= 0.0 ) call B%set(wrk,i,j)
            !call A%print(caption="matice A")
            !call B%print(caption="matice B")
            !print *,i,j
        !end do
        !print *, "i hotovo"
    end do
    !call pockej("submatrix hotovo")
    !call B%print(caption="matice B")
end function smatrix


subroutine Read(A,name)
    implicit none
    class(matrix_int), intent(in out) :: A
    character(len=*), intent(in) :: name
    integer :: fil
    integer(kind=ikind) :: i,n,m,nz
    integer(kind=ikind) :: wrk

    open(newunit=fil,file=name)
    read(fil,*) n,m,nz
    print *, n, m, nz
    call A%init(n,m)
    do i=1,nz
        read(fil,*) n,m,wrk
        call A%set(wrk,n,m)
    end do
    close(fil)
    call A%print()
end subroutine Read

subroutine Write(A,name)
    implicit none
    class(matrix_int), intent(in) :: A
    character(len=*), intent(in) :: name
    integer :: fil
    integer(kind=ikind) :: n,m
    integer(kind=ikind) :: wrk

    open(newunit=fil,file=name)
    write(fil,*) A%getn(),A%getm(),A%nonzero()
    do n=1,A%getn()
        do m=1,A%getm()
        wrk = A%get(n,m)
        if (wrk /= 0.0) write(fil,*) n,m,wrk
        end do
    end do
    close(fil)
    call A%print()
end subroutine Write


subroutine mtxtest_int(A)
    implicit none
    class(matrix_int), intent(in) :: A
    class(matrix_int), allocatable :: B,C
    integer(kind=ikind), dimension(3) :: ii
    integer(kind=ikind), dimension(2) :: jj


    print *,"zacinaji testy matic"

    allocate(B,source=A)
    allocate(C,source=A)
    call B%setverbosity(10)
    call C%setverbosity(10)
    ii = (/ 1 ,3, 5 /)
    jj = (/2,4/)
    call B%init(5_ikind,7_ikind)
    call B%print()
    call B%set(2_ikind,3_ikind,5_ikind)
    call B%set(2_ikind,5_ikind,2_ikind)
    call B%set(-8_ikind,1_ikind,1_ikind)
    call B%set(9_ikind,2_ikind,2_ikind)
    call B%set(-2_ikind,2_ikind,7_ikind)
    call B%print
    C = B%submatrix(ii,jj)
    call C%print
end subroutine mtxtest_int

subroutine transpose(A,B)
    class(matrix_int), intent(in out) :: A
    class(matrix_int), intent(in)     :: B

    integer(kind=ikind) :: n,m,i, nelem, j
    integer(kind=ikind), dimension(:), allocatable :: ii
    integer(kind=ikind), dimension(:), allocatable    :: v

    n = B%getn()
    m = B%getm()
    call A%init(m,n)
    do i=1,n
        call B%getrow(i,v,ii,nelem)
        do j=1, nelem
            call A%set(v(j),ii(j),i)
        end do
    end do

end subroutine transpose

end module mtx_int
