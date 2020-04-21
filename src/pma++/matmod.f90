module matmod
    use typy
    use mtx

    type, extends(matrix) :: symmul
        class(matrix), pointer :: v
        contains
        procedure :: get => getmm
        procedure, pass(A) :: set => setmm
        procedure :: mul => mulmm
        procedure setdata
    end type symmul

    contains
    subroutine setdata(A, M)
        class(symmul), intent(in out) :: A
        class(matrix), pointer :: M

        call A%resize(M%getm(), M%getm())
        A%V => M
    end subroutine setdata

    function getmm(A,i,j) result(v)
        class(symmul), intent(in) :: A
        integer(ikind), intent(in) :: i
        integer(ikind), intent(in) :: j
        real(rkind) :: v
        v = 0
        Stop "Nepripustna operace get"
    end function getmm

    subroutine setmm(r,A,i,j)
        class(symmul), intent(in out) :: A
        integer(ikind), intent(in) :: i
        integer(ikind), intent(in) :: j
        real(rkind),intent(in) :: r
        Stop "Nepripustna operace set"
    end subroutine setmm

    !> nasobi matici vektorem
    function mulmm(a,x, count) result(y)
        implicit none
        !> matice
        class(symmul), intent(in) :: a
        !> vektor
        real(kind=rkind), dimension(:), intent(in) :: x
        !> vysledek
        real(kind=rkind), dimension(:), allocatable :: y
        !> pocet operaci
        type(tcount), intent(inout), optional :: count
        !> pomocny vektor
        real(kind=rkind), dimension(:), allocatable :: wrk

        allocate(y(1:A%getn()));
        allocate(wrk(1:A%V%getn()))
        wrk = A%v%mul(x)
        y = A%v%mult(wrk)
    end function mulmm


end module matmod