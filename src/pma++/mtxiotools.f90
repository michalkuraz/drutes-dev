!> \file mtxiotools.f90
!! \brief V/V nastroje
!!
!<

!> V/V nastroje
!!
!! Matice: \n
!!
!! format je nezavisly na typu matice
!! tisknou se jen nenulove prvky
!! - prvni radek: pocet_radku pocet_nenul
!! - dalsi radky: radek sloupec hodnota
!! \n
!!
!! Vektor: \n
!! - prvni radek: pocet_prvku
!! - dalsi radky: hodnota
!<
module mtxiotools

    interface read_vector
        module procedure  read_vector_real
        module procedure  read_vector_int
    end interface

    interface write_vector
        module procedure  write_vector_real
        module procedure  write_vector_int
    end interface


    public :: read_matrix
    public :: write_matrix
    public :: read_vector
    public :: write_vector

    private
contains

    !> precte matici ze souboru
    subroutine read_matrix(A,name,fil,ilev,ierr,mode)
        use typy
        use mtx
        use pmatools
        implicit none
        !> matice
        class(matrix), intent(inout) :: A
        !> jmeno souboru, je-li zadano name i fil musi identifikovat totez
        character(len=*), intent(in), optional ::name
        !> identifikator kanalu, je-li zadano name i fil musi identifikovat totez
        integer, intent(in), optional :: fil
        !> stupen podrobnosti informaci
        integer, intent(in), optional :: ilev
        !> chybovy kod
        integer, intent(in), optional :: ierr
        !> operacni mod
        !> - 0 ... tise a v pripade problemu je ignoruj a zapis hlasku do ierr
        !> - 1 ... problem se soubory znamena interaktivni dotaz
        !> - 2 ... problem se soubory - uprednostni fil a pripadne zhavaruj
        !> - 3 ... problem se soubory - uprednostni name a pripadne zhavaruj
        !> - 4 ... problem se soubory - zhavaruj
        !>
        integer, intent(in), optional :: mode

        character(len=200) :: name1
        integer :: mode1
        integer :: fil1
        integer(kind=ikind) :: n, nz, i, j, k, maxi, maxj
        real(kind=rkind) :: val

        ! checknout parametry
        mode1 = 1
        call get_free_unit(fil1)
        call safe_open_read(name1,fil1,mode1)
        print *,"otevreno"
        read(fil1,*) n,nz
        print *,"n=",n," nz=",nz
        call A%init(n,n)
        print *,"realokovano"
        maxi = 0
        maxj = 0
        do k=1,nz
            if (i>maxi) maxi=i
            if (j>maxj) maxj=j
            read(fil1,*) i,j,val
            call A%set(val,i,j)
        end do

        call pockej("ctecka")
        call A%spy


    end subroutine read_matrix

    !> zapise matici ze souboru
    subroutine write_matrix(A,name,fil,ilev,ierr,mode)
        use typy
        use mtx
        implicit none
        !> matice
        class(matrix), intent(in) :: A
        !> jmeno souboru, je-li zadano name i fil musi identifikovat totez
        character(len=*), intent(in), optional ::name
        !> identifikator kanalu, je-li zadano name i fil musi identifikovat totez
        integer, intent(in), optional :: fil
        !> stupen podrobnosti informaci
        integer, intent(in), optional :: ilev
        !> chybovy kod
        integer, intent(in), optional :: ierr
        !> operacni mod
        !> - 0 ... tise a v pripade problemu je ignoruj a zapis hlasku do ierr
        !> - 1 ... problem se soubory znamena interaktivni dotaz
        !> - 2 ... problem se soubory - uprednostni fil a pripadne zhavaruj
        !> - 3 ... problem se soubory - uprednostni name a pripadne zhavaruj
        !> - 4 ... problem se soubory - zhavaruj
        !>
        integer, intent(in), optional :: mode
    end subroutine write_matrix


    !> precte vektor ze souboru
    subroutine read_vector_real(x,name,fil,ilev,ierr,mode,rlen)
        use typy
        use mtx
        implicit none
        !> vektor
        real(kind=rkind), dimension(:), intent(inout) :: x
        !> jmeno souboru, je-li zadano name i fil musi identifikovat totez
        character(len=*), intent(in), optional ::name
        !> identifikator kanalu, je-li zadano name i fil musi identifikovat totez
        integer, intent(in), optional :: fil
        !> stupen podrobnosti informaci
        integer, intent(in), optional :: ilev
        !> chybovy kod
        integer, intent(in), optional :: ierr
        !> operacni mod
        !> - 0 ... tise a v pripade problemu je ignoruj a zapis hlasku do ierr
        !> - 1 ... problem se soubory znamena interaktivni dotaz
        !> - 2 ... problem se soubory - uprednostni fil a pripadne zhavaruj
        !> - 3 ... problem se soubory - uprednostni name a pripadne zhavaruj
        !> - 4 ... problem se soubory - zhavaruj
        !>
        integer, intent(in), optional :: mode
        !> skutecna delka dat
        !> pokud delka vektoru a dat nesouhlasi a rlen neni zadano,
        !> dojde k ukonceni
        integer(kind=ikind), intent(out), optional :: rlen
    end subroutine read_vector_real

    !> zapise vektor do souboru
    subroutine write_vector_real(x,name,fil,ilev,ierr,mode)
        use typy
        use mtx
        implicit none
        !> vektor
        real(kind=rkind), dimension(:), intent(in) :: x
        !> jmeno souboru, je-li zadano name i fil musi identifikovat totez
        character(len=*), intent(in), optional ::name
        !> identifikator kanalu, je-li zadano name i fil musi identifikovat totez
        integer, intent(in), optional :: fil
        !> stupen podrobnosti informaci
        integer, intent(in), optional :: ilev
        !> chybovy kod
        integer, intent(in), optional :: ierr
        !> operacni mod
        !> - 0 ... tise a v pripade problemu je ignoruj a zapis hlasku do ierr
        !> - 1 ... problem se soubory znamena interaktivni dotaz
        !> - 2 ... problem se soubory - uprednostni fil a pripadne zhavaruj
        !> - 3 ... problem se soubory - uprednostni name a pripadne zhavaruj
        !> - 4 ... problem se soubory - zhavaruj
        !>
        integer, intent(in), optional :: mode
    end subroutine write_vector_real

    !> precte vektor ze souboru
    subroutine read_vector_int(x,name,fil,ilev,ierr,mode,rlen)
        use typy
        use mtx
        implicit none
        !> vektor
        integer(kind=ikind), dimension(:), intent(inout) :: x
        !> jmeno souboru, je-li zadano name i fil musi identifikovat totez
        character(len=*), intent(in), optional ::name
        !> identifikator kanalu, je-li zadano name i fil musi identifikovat totez
        integer, intent(in), optional :: fil
        !> stupen podrobnosti informaci
        integer, intent(in), optional :: ilev
        !> chybovy kod
        integer, intent(in), optional :: ierr
        !> operacni mod
        !> - 0 ... tise a v pripade problemu je ignoruj a zapis hlasku do ierr
        !> - 1 ... problem se soubory znamena interaktivni dotaz
        !> - 2 ... problem se soubory - uprednostni fil a pripadne zhavaruj
        !> - 3 ... problem se soubory - uprednostni name a pripadne zhavaruj
        !> - 4 ... problem se soubory - zhavaruj
        !>
        integer, intent(in), optional :: mode
        !> skutecna delka dat
        !> pokud delka vektoru a dat nesouhlasi a rlen neni zadano,
        !> dojde k ukonceni
        integer(kind=ikind), intent(out), optional :: rlen
    end subroutine read_vector_int

    !> zapise vektor do souboru
    subroutine write_vector_int(x,name,fil,ilev,ierr,mode)
        use typy
        use mtx
        use pmatools
        implicit none
        !> vektor
        integer(kind=ikind), dimension(:), intent(in) :: x
        !> jmeno souboru, je-li zadano name i fil musi identifikovat totez
        character(len=*), intent(in), optional ::name
        !> identifikator kanalu, je-li zadano name i fil musi identifikovat totez
        integer, intent(in), optional :: fil
        !> stupen podrobnosti informaci
        integer, intent(in), optional :: ilev
        !> chybovy kod
        integer, intent(in), optional :: ierr
        !> operacni mod
        !> - 0 ... tise a v pripade problemu je ignoruj a zapis hlasku do ierr
        !> - 1 ... problem se soubory znamena interaktivni dotaz
        !> - 2 ... problem se soubory - uprednostni fil a pripadne zhavaruj
        !> - 3 ... problem se soubory - uprednostni name a pripadne zhavaruj
        !> - 4 ... problem se soubory - zhavaruj
        !>
        integer, intent(in), optional :: mode

        !call safe_open_write(name1,fil1,mode1)
    end subroutine write_vector_int


end module mtxiotools
