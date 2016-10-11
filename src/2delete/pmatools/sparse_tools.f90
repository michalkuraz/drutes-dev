module sparse_tools
    use mytypes
   !definuje metody pro praci s ridkymi maticemi
    public :: lowest_degree,     &!nalezne vrchol s nejnizsim stupnem
              disjoint_covering, &!nalezne disjunktni pokryti
              add_level
    public :: LDU
    public :: Solve_LDU

contains
  subroutine lowest_degree(A, degree, ind)
    type(smtx), intent(in) :: A
    !matice reprezentujici graf, predpoklada se symetrie ve strukture
    !rozumny pozadavek -- S.P.D.
    integer(kind=ikind), intent(out) :: degree
    !nalezeny nejmensi stupen
    integer(kind=ikind), intent(out) :: ind
    !cislo prvniho vrcholu s nejnizsim stupnem
    !lokalni promenne
    integer(kind=ikind) :: i, wrk
    !telo procedury
    ind    = 1
    degree = A%rows(1,2)-A%rows(1,1)+1
    do i = 2, A%row
        if (A%rows(i,2)-A%rows(i,1)+1 < degree) then
            degree = A%rows(i,2)-A%rows(i,1)+1
            ind    = i
        end if
    end do
  end subroutine lowest_degree

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disjoint_covering(A,nd,covering)
    type(smtx), intent(in)                        :: A
    integer(kind=ikind), intent(in)               :: nd
    integer(kind=ikind), dimension(:), intent(in) :: covering

    !vypocet

  end subroutine disjoint_covering

  subroutine add_level(A,levels,nodes_avail, levstart)
    type(smtx), intent(in) :: A
    integer(kind=ikind), dimension(:), intent(inout) :: levels, &
                                                        nodes_avail, levstart
  end subroutine add_level


  !>LDU rozklad matice A
  !!
  !! tento rozklad je ulozen v polozce decomp matice\n
  !! tato polozka je v pripade potreby naalokovana
  !! matice L a U maji na diagonale jednicky
  !<
  subroutine LDU(A,perm,sz,info,ilev,ierr)
    use mytypes
    use typy
    !> rozkladana matice
    type(smtx), intent (in out) :: A
    !> permutacni vektor definujici poradi pri eliminaci
    integer(kind=ikind), intent(in), dimension(:) :: perm
    !> odhadovana velikost rozkladu, pripadne se upravi
    integer(kind=ikind), intent(in) :: sz
    !> informace o vypoctu
    type(info_type), intent(in out), optional :: info
    !> uroven zprav
    integer(kind=ikind), intent(in), optional :: ilev
    !> chybove hlaseni
    integer(kind=ikind), intent(out), optional :: ierr
  end subroutine LDU

  !> vyresi soustavu linearnich rovnic
  !!
  !! predpoklada se existence LDU rozkladu v matici\n
  !! pokud neni, je sestrojen a pouzije se Minimal Degree ordering
  !<
  subroutine Solve_LDU(A,b,x, info,ilev, ierr)
    use typy
    use mytypes
    !> system matrix
    type(smtx), intent(in out) :: A
    !> right hand side
    real(kind=rkind), dimension(:), intent(in) :: b
    !> solution
    real(kind=rkind), dimension(:), intent(out) :: x
    !> informace o vypoctu
    type(info_type), intent(in out), optional :: info
    !> uroven zprav
    integer(kind=ikind), intent(in), optional :: ilev
    !> chybove hlaseni
    integer(kind=ikind), intent(out), optional :: ierr

    ! lokalni promenne
    integer(kind=ikind)  :: i
    ! resim LDUx = b
    !postupne: 1. Ly = b
    !          2. Dz = y
    !          3. Ux = z
    !Step 1
    !x(perm(1)) = b(perm(1))
    do i=2, A%row
    end do
    !Step 2
    x = x/A%decomp%diag
    !Step 3
  end subroutine Solve_LDU


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sparse_tools
