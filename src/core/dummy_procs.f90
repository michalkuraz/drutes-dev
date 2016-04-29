module dummy_procs

  public :: dummy_scalar, dummy_vector, dummy_tensor, dummy_logical, time_check_ok
  
  contains

!   subroutine getder

    !> this function returns zero, in order to null zero terms from some particular PDE problems
    function dummy_scalar(pde_loc,i,quadpnt, r) result(null)
      use typy
      use pde_objs
      use global_objs
      class(pde_str), intent(in) :: pde_loc
      !>integer number
      integer(kind=ikind), intent(in) :: i
      type(integpnt_str), intent(in), optional :: quadpnt
      !> real number
      real(kind=rkind), dimension(:), intent(in), optional :: r
      real(kind=rkind) :: null

      null  = 0

    end function dummy_scalar

    !> this procedure returns vector of zeroes
    subroutine dummy_vector(pde_loc, i, quadpnt, r, vector_in, vector_out, scalar)
      use typy
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      !> integer num
      integer(kind=ikind), intent(in) :: i
      type(integpnt_str), intent(in), optional :: quadpnt
      !> real number
      real(kind=rkind), dimension(:), intent(in), optional :: r 
      !> input vector
      real(kind=rkind), dimension(:), intent(in), optional  :: vector_in
      !> output vector
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> scalar value to be supplied by zeroes
      real(kind=rkind), intent(out), optional :: scalar

      if (present(vector_out)) then
	vector_out = 0
      end if

      if (present(scalar)) then
	scalar = 0
      end if

    end subroutine dummy_vector
   
     !> this function returns zero matrix (tensor of second order) 
    subroutine dummy_tensor(pde_loc, i, quadpnt, r, tensor, scalar)
      use typy
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      !> integer number
      integer(kind=ikind), intent(in) :: i
      type(integpnt_str), intent(in), optional :: quadpnt
      !> real number
      real(kind=rkind), dimension(:), intent(in), optional :: r
      !> second order tensor to bye filled by zeroes
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor		
      !> scalar value to be supplied by zeroes
      real(kind=rkind), intent(out), optional :: scalar

      if (present(tensor)) then
	tensor = 0
      end if

      if (present(scalar)) then
	scalar = 0
      end if

    end subroutine dummy_tensor

    !> this function always returns .true.
    function dummy_logical() result(true)
      logical :: true

      true = .true.

    end function dummy_logical
    
    function time_check_ok(pde_loc) result(true)
      use typy
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      logical :: true
            
      true = .true.
      
    end function time_check_ok

end module dummy_procs