!> \file pmatest.f90
!! \brief obsahuje modul pro obal pro maticovy tester
!!
!<

!> obal pro maticovy tester
module pmatest
    public :: mtxtester
contains
    !> maticovy tester
    subroutine mtxtester
        use mtx
        use fullmatrix
        use sparsematrix
        use mtxtester_module
        implicit none

        type(fullmtx) :: af
        type(smtx) :: as
        print *,"mtxtester zacina"
        call matrixtester(af)
    end subroutine mtxtester
end module pmatest
