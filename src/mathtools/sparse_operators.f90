! Sparse operators for the existing matrix/solver interfaces.
! Assembly remains separate from repeated matrix-vector products.
module sparse_operators

  use typy
  use mtx, only: matrix, initmatrix

  implicit none
  private

  public :: csr_matrix, diagonal_matrix

  ! Fixed-pattern snapshot: changing its source never silently changes this copy.
  ! Call from_matrix again after changing the source or its sparsity pattern.
  type, extends(matrix) :: csr_matrix

    private
    integer(kind=ikind), dimension(:), allocatable :: row_start, column
    real(kind=rkind), dimension(:), allocatable :: values

  contains

    procedure :: init => init_csr
    procedure :: resize => init_csr
    procedure :: clear => clear_csr
    procedure :: get => get_csr
    procedure, pass(a) :: set => set_csr
    procedure :: getrow => row_csr
    procedure :: nonzero => nonzero_csr
    procedure :: mul => multiply_csr
    procedure :: mult => transpose_csr
    procedure :: from_matrix

  end type csr_matrix

  ! Inverse Jacobi diagonal is represented by its values, not by linked entries.
  type, extends(matrix) :: diagonal_matrix

    private
    real(kind=rkind), dimension(:), allocatable :: values

  contains

    procedure :: init => init_diagonal
    procedure :: resize => init_diagonal
    procedure :: clear => clear_diagonal
    procedure :: get => get_diagonal
    procedure, pass(a) :: set => set_diagonal
    procedure :: getrow => row_diagonal
    procedure :: nonzero => nonzero_diagonal
    procedure :: mul => multiply_diagonal
    procedure :: mult => multiply_diagonal
    procedure :: set_values => set_diagonal_values

  end type diagonal_matrix

contains

  subroutine init_csr(a, n, m)

    class(csr_matrix), intent(inout) :: a
    integer(kind=ikind), intent(in) :: n, m

    if (n < 0 .or. m < 0) then

      error stop 'CSR: negative matrix dimension'

    end if

    call initmatrix(a, n, m)

    if (allocated(a%row_start)) then

      deallocate (a%row_start, a%column, a%values)

    end if

    allocate (a%row_start(n + 1), a%column(0), a%values(0))
    a%row_start = 1

  end subroutine init_csr

  subroutine from_matrix(a, source)

    class(csr_matrix), intent(inout) :: a
    class(matrix), intent(in) :: source
    integer(kind=ikind), dimension(:), allocatable :: offsets, indices
    real(kind=rkind), dimension(:), allocatable :: row_values
    integer(kind=ikind) :: i, n, count, position

    ! Two linear passes avoid allocating n*m for a dense matrix with zeros.
    n = source%getn()
    allocate (offsets(n + 1))
    offsets(1) = 1

    do i = 1, n

      call source%getrow(i, row_values, indices, count)
      offsets(i + 1) = offsets(i) + count

    end do

    call a%init(n, source%getm())
    a%row_start = offsets
    deallocate (a%column, a%values)
    allocate (a%column(offsets(n + 1) - 1), a%values(offsets(n + 1) - 1))

    do i = 1, n

      call source%getrow(i, row_values, indices, count)

      if (count /= offsets(i + 1) - offsets(i)) then

        error stop 'CSR: source pattern changed during conversion'

      end if

      if (count == 0) then

        cycle

      end if

      if (any(indices(:count) < 1) .or. any(indices(:count) > source%getm())) then

        error stop 'CSR: source column outside matrix'

      end if

      position = offsets(i)
      ! Preserve source row order to avoid unnecessary rounding changes.
      a%column(position:position + count - 1) = indices(:count)
      a%values(position:position + count - 1) = row_values(:count)

    end do

  end subroutine from_matrix

  subroutine clear_csr(a)

    class(csr_matrix), intent(inout) :: a

    if (allocated(a%values)) then

      a%values = 0

    end if

  end subroutine clear_csr

  function get_csr(a, i, j) result(value)

    class(csr_matrix), intent(in) :: a
    integer(kind=ikind), intent(in) :: i, j
    real(kind=rkind) :: value
    integer(kind=ikind) :: k

    call check_indices(a, i, j)
    value = 0

    do k = a%row_start(i), a%row_start(i + 1) - 1

      if (a%column(k) == j) then

        value = a%values(k)
        exit

      end if

    end do

  end function get_csr

  subroutine set_csr(r, a, i, j)

    real(kind=rkind), intent(in) :: r
    class(csr_matrix), intent(inout) :: a
    integer(kind=ikind), intent(in) :: i, j
    integer(kind=ikind) :: k

    call check_indices(a, i, j)

    do k = a%row_start(i), a%row_start(i + 1) - 1

      if (a%column(k) == j) then

        a%values(k) = r
        return

      end if

    end do

    if (r /= 0) then

      error stop 'CSR: new entry requires rebuilding the snapshot'

    end if

  end subroutine set_csr

  subroutine row_csr(a, i, v, jj, nelem, mask)

    class(csr_matrix), intent(in) :: a
    integer(kind=ikind), intent(in) :: i
    real(kind=rkind), dimension(:), allocatable, intent(inout) :: v
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: jj
    integer(kind=ikind), intent(out) :: nelem
    logical, dimension(:), intent(in), optional :: mask
    integer(kind=ikind) :: k, capacity

    if (i < 1 .or. i > a%getn()) then

      error stop 'CSR: row outside matrix'

    end if

    if (present(mask)) then

      if (size(mask, kind=ikind) /= a%getm()) then

        error stop 'CSR: column mask dimension mismatch'

      end if

    end if

    capacity = a%row_start(i + 1) - a%row_start(i)
    call reserve_row(v, jj, capacity)
    nelem = 0

    do k = a%row_start(i), a%row_start(i + 1) - 1

      if (a%values(k) == 0) then

        cycle

      end if

      if (present(mask)) then

        if (.not. mask(a%column(k))) then

          cycle

        end if

      end if

      nelem = nelem + 1
      v(nelem) = a%values(k)
      jj(nelem) = a%column(k)

    end do

  end subroutine row_csr

  function nonzero_csr(a) result(number)

    class(csr_matrix), intent(in) :: a
    integer(kind=ikind) :: number

    number = 0

    if (allocated(a%values)) then

      number = count(a%values /= 0, kind=ikind)

    end if

  end function nonzero_csr

  function multiply_csr(a, x, count) result(y)

    class(csr_matrix), intent(in) :: a
    real(kind=rkind), dimension(:), intent(in) :: x
    type(tcount), intent(inout), optional :: count
    real(kind=rkind), dimension(:), allocatable :: y
    integer(kind=ikind) :: i, k

    if (size(x, kind=ikind) /= a%getm()) then

      error stop 'CSR: product dimension mismatch'

    end if

    allocate (y(a%getn()))
    y = 0

    if (.not. allocated(a%row_start)) then

      return

    end if

    do i = 1, a%getn()

      do k = a%row_start(i), a%row_start(i + 1) - 1

        y(i) = y(i) + a%values(k)*x(a%column(k))

      end do

    end do

    if (present(count)) then

      count%ad = count%ad + size(a%values, kind=ikind)
      count%mul = count%mul + size(a%values, kind=ikind)

    end if

  end function multiply_csr

  function transpose_csr(a, x, count) result(y)

    class(csr_matrix), intent(in) :: a
    real(kind=rkind), dimension(:), intent(in) :: x
    type(tcount), intent(inout), optional :: count
    real(kind=rkind), dimension(:), allocatable :: y
    integer(kind=ikind) :: i, k

    if (size(x, kind=ikind) /= a%getn()) then

      error stop 'CSR: transpose product dimension mismatch'

    end if

    allocate (y(a%getm()))
    y = 0

    if (.not. allocated(a%row_start)) then

      return

    end if

    do i = 1, a%getn()

      do k = a%row_start(i), a%row_start(i + 1) - 1

        y(a%column(k)) = y(a%column(k)) + a%values(k)*x(i)

      end do

    end do

    if (present(count)) then

      count%ad = count%ad + size(a%values, kind=ikind)
      count%mul = count%mul + size(a%values, kind=ikind)

    end if

  end function transpose_csr

  subroutine init_diagonal(a, n, m)

    class(diagonal_matrix), intent(inout) :: a
    integer(kind=ikind), intent(in) :: n, m

    if (n < 0 .or. n /= m) then

      error stop 'Diagonal operator requires square nonnegative dimensions'

    end if

    call initmatrix(a, n, m)

    if (allocated(a%values)) then

      deallocate (a%values)

    end if

    allocate (a%values(n))
    a%values = 0

  end subroutine init_diagonal

  subroutine set_diagonal_values(a, values)

    class(diagonal_matrix), intent(inout) :: a
    real(kind=rkind), dimension(:), intent(in) :: values

    call a%init(size(values, kind=ikind), size(values, kind=ikind))
    a%values = values

  end subroutine set_diagonal_values

  subroutine clear_diagonal(a)

    class(diagonal_matrix), intent(inout) :: a

    if (allocated(a%values)) then

      a%values = 0

    end if

  end subroutine clear_diagonal

  function get_diagonal(a, i, j) result(value)

    class(diagonal_matrix), intent(in) :: a
    integer(kind=ikind), intent(in) :: i, j
    real(kind=rkind) :: value

    call check_indices(a, i, j)
    value = 0

    if (i == j) then

      value = a%values(i)

    end if

  end function get_diagonal

  subroutine set_diagonal(r, a, i, j)

    real(kind=rkind), intent(in) :: r
    class(diagonal_matrix), intent(inout) :: a
    integer(kind=ikind), intent(in) :: i, j

    call check_indices(a, i, j)

    if (i == j) then

      a%values(i) = r

    else if (r /= 0) then

      error stop 'Nonzero off-diagonal entry in diagonal operator'

    end if

  end subroutine set_diagonal

  subroutine row_diagonal(a, i, v, jj, nelem, mask)

    class(diagonal_matrix), intent(in) :: a
    integer(kind=ikind), intent(in) :: i
    real(kind=rkind), dimension(:), allocatable, intent(inout) :: v
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: jj
    integer(kind=ikind), intent(out) :: nelem
    logical, dimension(:), intent(in), optional :: mask

    call check_indices(a, i, i)
    call reserve_row(v, jj, 1_ikind)
    nelem = 0

    if (present(mask)) then

      if (size(mask, kind=ikind) /= a%getm()) then

        error stop 'Diagonal operator: column mask dimension mismatch'

      end if

      if (.not. mask(i)) then

        return

      end if

    end if

    if (a%values(i) /= 0) then

      nelem = 1
      v(1) = a%values(i)
      jj(1) = i

    end if

  end subroutine row_diagonal

  function nonzero_diagonal(a) result(number)

    class(diagonal_matrix), intent(in) :: a
    integer(kind=ikind) :: number

    number = 0

    if (allocated(a%values)) then

      number = count(a%values /= 0, kind=ikind)

    end if

  end function nonzero_diagonal

  function multiply_diagonal(a, x, count) result(y)

    class(diagonal_matrix), intent(in) :: a
    real(kind=rkind), dimension(:), intent(in) :: x
    type(tcount), intent(inout), optional :: count
    real(kind=rkind), dimension(:), allocatable :: y

    if (size(x, kind=ikind) /= a%getm()) then

      error stop 'Diagonal operator: product dimension mismatch'

    end if

    allocate (y(a%getn()))
    y = 0

    if (allocated(a%values)) then

      y = a%values*x

    end if

    if (present(count)) then

      count%mul = count%mul + size(x, kind=ikind)

    end if

  end function multiply_diagonal

  subroutine reserve_row(values, columns, capacity)

    real(kind=rkind), dimension(:), allocatable, intent(inout) :: values
    integer(kind=ikind), dimension(:), allocatable, intent(inout) :: columns
    integer(kind=ikind), intent(in) :: capacity

    if (allocated(values)) then

      if (size(values, kind=ikind) < capacity) then

        deallocate (values)

      end if

    end if

    if (.not. allocated(values)) then

      allocate (values(capacity))

    end if

    if (allocated(columns)) then

      if (size(columns, kind=ikind) < capacity) then

        deallocate (columns)

      end if

    end if

    if (.not. allocated(columns)) then

      allocate (columns(capacity))

    end if

  end subroutine reserve_row

  subroutine check_indices(a, i, j)

    class(matrix), intent(in) :: a
    integer(kind=ikind), intent(in) :: i, j

    if (i < 1 .or. i > a%getn() .or. j < 1 .or. j > a%getm()) then

      error stop 'Sparse operator: index outside matrix'

    end if

  end subroutine check_indices

end module sparse_operators
