!> \file datasetup.f90
!! \brief konstruktory matic
!<

!> konstruktory matic
module datasetup
    public :: Laplace2D
    public :: Laplace3D
    public :: Hilbert

    contains
    !> diskretizace 2D Laplaceova operatoru
    !!
    !! sit je rovnomerna s krokem h=1/nx
    !<
    subroutine Laplace2D(a,nx,ny)
        use mtx
        use typy
        implicit none
        !> vytvarena matice
        class(matrix), intent(in out) :: a
        !> pocet uzlu ve smeru x
        integer(kind=ikind), intent(in) :: nx
        !> pocet uzlu ve smeru y
        integer(kind=ikind), intent(in) :: ny

        integer(kind=ikind) :: n,i,j,ii

        n = nx*ny
        call a%init(n,n)
        do i=1,nx
            do j=1,ny
                ii = (i-1)*ny+j
                !print *, i,j,ii
                call a%set(4.0_rkind,ii,ii)
                !print *,"a"
                if (j>1)  call a%set(-1.0_rkind,ii,ii-1)
                !print *,"b"
                if (j<ny) call a%set(-1.0_rkind,ii,ii+1)
                !print *,"c"
                if (i>1)  call a%set(-1.0_rkind,ii,ii-ny)
                !print *,"d"
                if (i<nx) call a%set(-1.0_rkind,ii,ii+ny)
                !print *,"e"
            end do
        end do
    end subroutine Laplace2D

    !> diskretizace 3D Laplaceova operatoru
    !!
    !! sit je rovnomerna s krokem h=1/(nx+1)
    !<
    subroutine Laplace3D(a,nx,ny,nz)
        use mtx
        use typy
        implicit none
        !> vytvarena matice
        class(matrix), intent(in out) :: a
        !> pocet uzlu ve smeru x
        integer(kind=ikind), intent(in) :: nx
        !> pocet uzlu ve smeru y
        integer(kind=ikind), intent(in) :: ny
        !> pocet uzlu ve smeru z
        integer(kind=ikind), intent(in) :: nz

        integer(kind=ikind) :: n,i,j,k,ii

        n = nx*ny*nz
        call a%init(n,n)
        !print *, "lap1"
        do i=1,nx
            do j=1,ny
                do k = 1,nz
                    ii = (i-1)*ny*nz+(j-1)*nz+k
                    !print *, i,j,k,ii
                    call a%set(6.0_rkind,ii,ii)
                    !print *,"a"
                    if (k>1)  call a%set(-1.0_rkind,ii,ii-1)
                    if (k<nz) call a%set(-1.0_rkind,ii,ii+1)
                    !print *, "a1"
                    if (j>1)  call a%set(-1.0_rkind,ii,ii-nz)
                    if (j<ny) call a%set(-1.0_rkind,ii,ii+nz)
                    !print *, "a2"
                    if (i>1)  call a%set(-1.0_rkind,ii,ii-ny*nz)
                    !print *, "a2a"
                    if (i<nx) call a%set(-1.0_rkind,ii,ii+ny*nz)
                    !print *, "a3"
                end do
            end do
        end do
    end subroutine Laplace3D




    !> \brief vytvori Hilbertovu matici
    !!
    !! \param a  vytvarena matice
    !! \param n  rozmer konstruovane matice
    !!
    !! v pripade potreby realokuje
    !!
    subroutine Hilbert(a,n)
        use typy
        use mtx
        implicit none
        class(matrix), intent(inout) :: a
        integer(kind=ikind),intent(in) :: n
        integer(kind=ikind) :: i,j
        !print *,"Hilbert:volam resize"
        call a%init(n,n)
        !print *,"Hilbert:po resize"
        do i=1,n
            do j=1,n
                call a%set(1.0_rkind/(i+j-1),i,j)
            end do
        end do

    end subroutine Hilbert
end module datasetup
