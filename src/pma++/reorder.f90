module reorder
    use typy
    use mtx
    implicit none

    !> i-th level is in levlist from levstart(i) to levstart(i+1)-1
    !> nlev - amount of levels
    !> levstart(i) is defibed for i=1,nlev+1
    type :: vrstvy
        integer(kind=ikind), dimension(:), allocatable, private :: levlist
        integer(kind=ikind), dimension(:), allocatable, private :: levstart
        integer(kind=ikind), private :: nlev
        contains
        procedure, private :: init
        procedure, private :: addlevel
        procedure, private :: order
    end type vrstvy

    public :: CM
    public :: RCM
!    public :: testreord

    contains
    subroutine  CM(A,perm)
        class(matrix), intent(in) :: A
        integer(ikind),dimension(:), allocatable, intent(in out) :: perm
        type(vrstvy) :: levset
        logical, dimension(:), allocatable :: mapa
        integer(ikind) :: first, last
        ! find border node
        first = A%getn()/2
        allocate(mapa(1:A%getn()))
        mapa = .true.
        call levset%order(A,first,mapa)
!        print *, levset%nlev
        last = levset%levlist(levset%levstart(levset%nlev+1)-1)
        mapa = .true.
        call levset%order(A,last,mapa)
!        print *, levset%nlev
        first = levset%levlist(levset%levstart(levset%nlev+1)-1)

        ! order by levels
        mapa = .true.
!        print *, levset%nlev
        call levset%order(A,first,mapa)
        if ( allocated(perm)) deallocate(perm)
        allocate(perm(1:A%getn()))
        perm = levset%levlist

    end subroutine CM

    subroutine  RCM(A,perm)
        class(matrix), intent(in) :: A
        integer(ikind),dimension(:), allocatable, intent(in out) :: perm
        integer(ikind) :: i,w
        call CM(A,perm)
        do i=1,A%getn()/2  !! reversing order of Cuthill-McKae
            w = perm(i)
            perm(i) = perm(A%getn()+1-i)
            perm(A%getn()+1-i) = w
        end do
    end subroutine RCM

    subroutine init(ls,A,i)
        class(vrstvy), intent(inout) :: ls
        class(matrix),intent(in) :: A
        integer(ikind), intent(in) :: i
        ls%nlev = 1
        if ( allocated(ls%levlist)) deallocate(ls%levlist)
        if ( allocated(ls%levstart)) deallocate(ls%levstart)
        allocate(ls%levlist(1:A%getn()))
        allocate(ls%levstart(1:max(A%getn()+1,2)))
        ls%levstart(1) = 1
        ls%levlist(1) = i
        ls%levstart(2) = 2
    end subroutine init

    subroutine addlevel(ls,A,mapa)
        class(vrstvy), intent(in out) :: ls
        class(matrix), intent(in) :: A
        logical, dimension(:), intent(in out) :: mapa
        integer(ikind) :: i, j, lev1, lev2, lev3, nelem
        real(rkind), dimension(:),allocatable :: v
        integer(ikind), dimension(:),allocatable :: jj

        lev1 = ls%levstart(ls%nlev)
        lev2 = ls%levstart(ls%nlev+1)
        lev3 = lev2
        do i = lev1, lev2-1
            call A%getrow(ls%levlist(i),v,jj,nelem)
            do j = 1, nelem
                if ( mapa(jj(j)) ) then
                    ls%levlist(lev3) = jj(j)
                    mapa(jj(j)) = .false.
                    lev3 = lev3+1
                end if
            end do
        end do
        ls%nlev = ls%nlev+1
        ls%levstart(ls%nlev+1) = lev3

    end subroutine addlevel

    subroutine order(ls,A,start,mapa)
        class(vrstvy), intent(inout) :: ls
        class(matrix), intent(in) :: A
        integer(ikind), intent(in) :: start
        logical, dimension(:), intent(in out) :: mapa
        integer(ikind) :: nord
        call ls%init(A,start)
        mapa(start) = .false.
        do
            call ls%addlevel(A,mapa)
            nord = ls%levstart(ls%nlev+1) - 1
!            print *, "ocislovano = ",nord
            if ( nord == A%getn()) exit
        end do

    end subroutine order

!    subroutine testreord
!        use sparsematrix
!        use datasetup
!        use solvers

!        type(smtx) :: A, LDU1, LDU2, B1, B2, LDU3, LDU4
!        integer(ikind), dimension(:),allocatable :: perm
!        integer(ikind), dimension(:),allocatable :: perm1
!        real(rkind), dimension(:), allocatable :: x,xr,b, bb1, bb2
!        type(tcount) :: inf1,inf2
!        integer(ikind) :: i
!        print *, "test preusporadani"
!        call Laplace2D(A,15_ikind,15_ikind)
!        allocate(perm(1:A%getn()))
!        call A%spy
!        call CM(A,perm)
!        print *, perm
!        call RCM(A,perm1)
!        print *, perm1
!        call copyperm(source=A, dest=B1,permi=perm, permj=perm)
!        call B1%spy
!        call LDU(B1,LDU3,inf1,0)
!        call LDU3%spy
!        call copyperm(source=A, dest=B2,permi=perm1, permj=perm1)
!        call B2%spy
!        call LDU(B2,LDU4,inf2,0)
!        call LDU4%spy
!        ! ted reseni soustavy
!        allocate(x(1:a%getn()));
!        allocate(xr(1:a%getn()));
!        allocate(b(1:a%getn()));
!        allocate(bb1(1:a%getn()));
!        allocate(bb2(1:a%getn()));
!        x = 0
!        !> construct of test solution
!        do i=1,A%getn()
!            xr(i) = i
!        end do
!        !> compute test rhs
!        b = A%mul(xr)
!        !> permute rhs according to founded  permutation
!        bb1 = b(perm)
!        !> back substitution
!        call LDUback(LDU3,bb1,x)
!        !> back permutation of solution
!        x(perm) = x
!        print *, bb1
!        print *, ""
!        print *, x
!!        call LDU(A,LDU1,inf1,6,perm1=perm,ilev=100)
! !       call LDU(A,LDU2,inf2,6,perm1=perm1,ilev=100)
!        print *, inf1
!        print *, inf2
!    end subroutine testreord




end module reorder
