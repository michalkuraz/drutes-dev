!> \file mtxtester.f90
!! \brief testovadlo matic


module mtxtester_module
    public :: matrixtester
contains

    !> \brief  testovadlo matic
    !!
    !! \param a matice, nemusi byt alokovana. jde jen o moznost polymorfie
    !!
    subroutine matrixtester(a)
        use pmatypy
        use pmatools
        use mtx
        use datasetup
        use solvers
        implicit none
        class(matrix), intent(in) :: a
        class(matrix), allocatable :: b
        class(matrix), allocatable :: c,c1
        class(matrix), allocatable :: L
        class(matrix), allocatable :: D
        class(matrix), allocatable :: U
        real(kind=rkind), dimension(:), allocatable :: x,x0,xr,rhs
        real(kind=rkind), dimension(:), allocatable :: xldu
        real(kind=rkind), dimension(:), allocatable :: xj
        real(kind=rkind), dimension(:), allocatable :: xgs,xcg
        real(kind=rkind), dimension(:), allocatable :: res
        integer(kind=ikind), dimension(:), allocatable :: p1,p2
        integer(kind=ikind) :: i,n
        type(tcount) :: opc1

        print *,"testy matic zacinaji"
        allocate(b,source=a)
        allocate(c,source=a)
        allocate(c1,source=a)
        allocate(l,source=a)
        allocate(d,source=a)
        allocate(u,source=a)
        !konstrukce
        call b%init(5_ikind,5_ikind)
        call b%print
        !manipulace s prvky
        call Hilbert(b,6_ikind)
        allocate(p1(1:6),p2(1:6))
        call b%print(ncol=6_ikind)
        call LDU(b,c)
        print *,"delam vektor"
        do i=1,6
            p1(i)=i
            p2(i)=i
        enddo
        print *,"volam split"
        call Split(c,L,D,U,p1,p2)
        call L%print(ncol=6_ikind,caption="L matice")
        call D%print(ncol=6_ikind,caption="D matice")
        call U%print(ncol=6_ikind,caption="U matice")

        call c1%clone(c)
        call c%mulm(L,D)
        call c1%mulm(c,U)
        call c1%print(ncol=6_ikind,caption="soucim LDU matice")
        call c1%subm(b)
        call c1%print(ncol=6_ikind,caption="rozdil LDU matice a originalu")
        print *, "norma rozdilu=",c1%normF()
        call pockej("jdu na Laplace")
        call Laplace2D(b,10_ikind,10_ikind)
        call B%Print()
        call B%Spy()
        n = B%getn()
        print *,"n=",n
        allocate(xr(1:n))
        xr = 2
        rhs = B%mul(xr)
        x = xr-xr
        x0 = x
        print *," vstupni data"
        print *,"                  i                x                " &
              // "            xr                            rhs"
        do i=1,n
            print *,i,x(i),xr(i),rhs(i)
        end do
        call pockej("jdu resit")
        ! napred eliminace
        ! potom iterace
        x = x0
        print *, "Jacobiova metoda"
        call jacobi(B,rhs,x,ilev1=1, reps1=1.0e-10_rkind, maxit1=1000*B%getn())
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        call pockej("Jacobiova metoda konci. prechod na nejvetsi spad")
        x=x0
        call SD(B,rhs,x, ilev1=1, reps1=1.0e-10_rkind, itmax1=1000*B%getn(),&
                 opcnt1=opc1)
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        print *,"celkovy pocet operaci"
        call print_info(opc1)

       call pockej("Druhe kolo nejvetsi spad")
        call SD(B,rhs,x, ilev1=1, reps1=1.0e-10_rkind, itmax1=1000*B%getn(),&
                 opcnt1=opc1)
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        print *,"celkovy pocet operaci"
        call print_info(opc1)

        call pockej("Prechod na sdruzene gradienty")
        x=x0
        call CG(B,rhs,x, ilev1=1, reps1=1.0e-10_rkind, itmax1=1000*B%getn(),&
                 opcnt1=opc1)
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        print *,"celkovy pocet operaci"
        call print_info(opc1)

       call pockej("Druhe kolo gradientu")
        call CG(B,rhs,x, ilev1=1, reps1=1.0e-10_rkind, itmax1=1000*B%getn(),&
                 opcnt1=opc1)
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        print *,"celkovy pocet operaci"
        call print_info(opc1)

        call pockej("CG konci")

        call pockej("Prechod na sdruzene gradienty-normal")
        x=x0
        call CGnormal(B,rhs,x, ilev1=1, reps1=1.0e-10_rkind,&
                 itmax1=1000*B%getn(),&
                 opcnt1=opc1)
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        print *,"celkovy pocet operaci"
        call print_info(opc1)

       call pockej("Druhe kolo gradientu-normal")
        call CGnormal(B,rhs,x, ilev1=1, reps1=1.0e-10_rkind,&
                 itmax1=1000*B%getn(),&
                 opcnt1=opc1)
        res = rhs - B%mul(x)
        do i=1,ubound(x,1)
            print *, i, x(i) ,xr(i), rhs(i), res(i)
        end do
        print *,"celkovy pocet operaci"
        call print_info(opc1)

        call pockej("CG-normal konci")


    end subroutine matrixtester
end module mtxtester_module
