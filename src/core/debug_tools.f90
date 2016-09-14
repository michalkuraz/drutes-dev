module debug_tools
  private :: print_real_matrix, print_int_matrix, print_real_vector, print_int_vector, print_real_vector4
!   public :: sparse_gem_pig
!   public :: sparse_gem_pig_AtA
  public :: wait
  private :: print_quadpnt
    
  interface printmtx
    module procedure print_real_matrix
    module procedure print_int_matrix
    module procedure print_real_vector
    module procedure print_int_vector
    module procedure print_sparse_matrix
    module procedure print_smartmatrix_i
    module procedure print_smartarray_i
    module procedure print_logical_array
    module procedure print_real_vector4
    module procedure print_quadpnt
  end interface printmtx
 

  contains

    subroutine print_quadpnt(quadpnt, filunit, name)
      use global_objs
      use core_tools
      use globals
      
      type(integpnt_str), intent(in) :: quadpnt
      integer, intent(in), optional :: filunit
      character(len=*), intent(in), optional :: name
      
      integer :: filloc
      integer :: ierr
      logical :: op   
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      
      write(unit=filloc, fmt=*) "type:", quadpnt%type_pnt
      
      if (quadpnt%type_pnt == "ndpt" .or. quadpnt%type_pnt == "obpt") then
        write(unit=filloc, fmt=*) "order of point:", quadpnt%order
      end if
      
      if (quadpnt%type_pnt == "gqnd") then
        write(unit=filloc, fmt=*) "order of element:", quadpnt%element
      end if
      
      write(unit=filloc, fmt=*) "column:", quadpnt%column
      
      if (quadpnt%ddlocal) then
        write(unit=filloc, fmt=*) "using subdomain local data"
        write(unit=filloc, fmt=*) "subdomain id", quadpnt%subdom
        if (quadpnt%extended) then
          write(unit=filloc, fmt=*) "node from extended subdomain (see subcycling man)"
        else
	  write(unit=filloc, fmt=*) "node inside the subdomain"
	end if
      else
        write(unit=filloc, fmt=*) "using global data"
      end if
      
      if (.not. quadpnt%globtime) then
        write(unit=filloc, fmt=*) "not using the global time, solution will be interpolated"
        write(unit=filloc, fmt=*) "time 4 use:", quadpnt%time4eval
      else
        write(unit=filloc, fmt=*) "using the global time"
      end if
      
  
      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
      
    end subroutine print_quadpnt
  
    subroutine print_sparse_matrix(A, filunit, name)
      use sparsematrix
      use mtxiotools
      use globals
      use core_tools
      class(smtx), intent(in out) :: A
      integer, intent(in), optional :: filunit
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
	filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
	filloc = terminal
      end if
      

      call a%print()

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if

    end subroutine print_sparse_matrix

    subroutine print_real_matrix(a, filunit, name)
      use typy
      use globals
      use core_tools

      real(kind=rkind), dimension(:,:), intent(in out) :: a
      integer, intent(in), optional :: filunit     
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      
      do i=lbound(a,1), ubound(a,1)
        write(unit=filloc, fmt=*)  i, "|",  a(i,:)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
      
    end subroutine print_real_matrix

    subroutine print_int_matrix(a, filunit, name)
      use typy
      use globals
      use core_tools
      
      integer(kind=ikind), dimension(:,:), intent(in out) :: a
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      
      do i=lbound(a,1), ubound(a,1)
        write(unit=filloc, fmt=*)  i, "|",  a(i,:)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
      
    end subroutine print_int_matrix

  !>
    !! vytiskne vektor
    !<
    subroutine print_real_vector(V, filunit, name)
    ! vytiskne vektor, pocet sloupcu tisku je nc
      use typy
      use globals
      use core_tools
      
      !parametry
      real(kind=rkind), dimension(:), intent(in out) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
     

      do i=lbound(V,1),ubound(V,1)
	 write(unit=filloc, fmt=*) "radek:", i, "hodnota:", V(i)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  
    end subroutine print_real_vector

  !>
    !! vytiskne vektor
    !<
    subroutine print_real_vector4(V, filunit, name)
    ! vytiskne vektor, pocet sloupcu tisku je nc
      use typy
      use globals
      use core_tools
      
      !parametry
      real(4), dimension(:), intent(in out) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
     

      do i=lbound(V,1),ubound(V,1)
	 write(unit=filloc, fmt=*) "radek:", i, "hodnota:", V(i)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  
    end subroutine print_real_vector4

    subroutine print_int_vector(V, filunit, name)
    ! vytiskne vektor, pocet sloupcu tisku je nc
      use typy
      use globals
      use core_tools
      
      !parametry
      integer(kind=ikind), dimension(:), intent(in) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if

      do i=lbound(V,1),ubound(V,1)
	 write(unit=filloc, fmt=*) "radek:", i, "hodnota:", V(i)
      end do
  
      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  

    end subroutine print_int_vector
    
    subroutine print_smartmatrix_i(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      class(smartarray_int), dimension(:), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
            
      do i=lbound(array,1),ubound(array,1)
	if (.not. allocated(array(i)%data)) then
	  print *, "no values to print at row:", i
	else
	  write(unit=filloc, fmt=*) i, "|", array(i)%data(1:array(i)%pos)
	  write(unit=filloc, fmt=*) "-----------------------------------------------"
	end if
      end do
      
      call flush(filloc)
      
    end subroutine print_smartmatrix_i

    subroutine print_smartarray_i(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      class(smartarray_int), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      do i=1, array%pos
        write(unit=filloc, fmt=*)  i, "|", array%data(i)
      end do

      
      call flush(filloc)
      
    end subroutine print_smartarray_i
    
    subroutine print_logical_array(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      logical, dimension(:), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      do i=1, ubound(array,1)
        write(unit=filloc, fmt=*) i, "|", array(i)
      end do
      
    end subroutine print_logical_array




    subroutine wait(ch) 
      use globals       
      character(len=*), optional, intent(in) :: ch
      
      if (present(ch)) then
        print *, "-------------"
        print *, trim(ch)
        print *, "-------------"
      end if

      print *, "press [ENTER] to continue"
  
      call flush(terminal)
  
      read(*,*)


    end subroutine wait
    
!      !> this procedure bears a codename sparse_gem_pig because it creates full matrix out of sparse matrix and solves it on using full Gauss elimination, only for debugging
!   subroutine sparse_gem_pig(A,b,x,ptype,ilev,ierr, itmax, reps_rel, reps_abs, info, it)
!     use linalg
!     use typy
!     use sparsematrix
!     
!     type(smtx), intent(in out) :: a
!     real(kind=rkind), dimension(:), intent(in) :: b
!     real(kind=rkind), dimension(:), intent(in out) :: x
!   !> level of information, default = 0
!       integer, intent(in), optional                  :: ilev
!       !> kind of preconditioner, default 0
!       integer, intent(in), optional                  :: ptype
!       !> error message\n 0 .. OK\n
!       !>                1 .. after itmax iterions is not founded sufficiently small relative error
!       integer, intent(out), optional                 :: ierr
!       !> maximum allowed iterations, default = 500
!       integer(kind=ikind), intent(in), optional      :: itmax
!       !> maximal relative error
!       real(kind=rkind), intent(in), optional         :: reps_rel
!       !> maximal absolute error
!       real(kind=rkind), intent(in), optional         :: reps_abs
!       !> operations count
!       type(info_type), intent(inout), optional       :: info
!       !> iteration number
!       integer(kind=ikind), intent(out), optional     :: it
!     real(kind=rkind), dimension(:,:), allocatable :: matice
!     integer(kind=ikind) :: i
!     
!     allocate(matice(ubound(b,1), ubound(b,1)))
! 
!     matice = 0.0_rkind
! 
!     do i=1,ubound(a%vals,1)
!       if (a%ii(i) /=0 .or. a%jj(i) /= 0) then  public :: copy_mtx
!         matice(a%ii(i), a%jj(i)) =  matice(a%ii(i), a%jj(i)) + a%vals(i)
!       end if
!     end do
! 
!     call gem(matice, b, x)
! 
!   end subroutine sparse_gem_pig
! 
!   
!     !> this procedure bears a codename sparse_gem_pig because it creates full matrix out of sparse matrix and solves it on using full Gauss elimination, only for debugging
!     subroutine sparse_gem_pig_AtA(A,b,x,ptype,ilev,ierr, itmax, reps, info,normmul)
!         use mytypes
!         use typy
!         use linAlg
!         implicit none
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! parametry
!         !> system matrix - supposed be symmetric positive definite\n
!         !> inout because of possibility of initalizing of preconditioner
!         type(smtx), intent(in out)                  :: A
!         !> right hand side
!         real(kind=rkind), dimension(:), intent(in) :: b
!         !> solution, on input initial aproximation
!         real(kind=rkind), dimension(:), intent(in out) :: x
!         !> level of information, default = 0
!         integer, intent(in), optional           :: ilev
!         !> kind of preconditioner, default 0
!         integer, intent(in), optional           :: ptype
!         !> error message\n 0 .. OK\n
!         !>                1 .. after itmax iterions is not founded sufficiently small relative error
!         integer, intent(out), optional          :: ierr
!         !> maximum allowed iterations, default = 500
!         integer(kind=ikind), intent(in), optional           :: itmax
!         !> wanted relative error, default=  1e-5
!         real(kind=rkind), intent(in), optional  :: reps
!         !> operations count
!         type(info_type), intent(inout), optional  :: info
!         !> zda pro normeq=true delat prenasobeni
!         logical, optional, intent(in) :: normmul
!     real(kind=rkind), dimension(:,:), allocatable :: matice, matice2
!     real(kind=rkind), dimension(:), allocatable :: Atb
!     integer(kind=ikind) :: i
!     
!     allocate(matice(ubound(b,1), ubound(x,1)))
!     allocate(matice2(ubound(x,1), ubound(b,1)))
!     allocate(Atb(ubound(x,1)))
!     
!     matice = 0.0_rkind
!     matice2 = 0.0_rkind
! 
!     do i=1,ubound(a%vals,1)
!       if (a%ii(i) > 0 .or. a%jj(i) > 0) then
!         matice(a%ii(i), a%jj(i)) =  matice(a%ii(i), a%jj(i)) + a%vals(i)
!       end if
!     end do
! 
!     do i=1,ubound(a%vals,1)
!       if (a%ii(i) > 0 .or. a%jj(i) > 0 .and.  abs(a%vals(i)) > epsilon(a%vals(i))) then
!         matice2(a%jj(i), a%ii(i)) =  matice2(a%jj(i), a%ii(i)) + a%vals(i)
!       end if
!     end do
!     
! 
!     
!     matice = matmul(matice2,matice)
! 
!     
!     Atb = matmul(matice2,b)
!     
! !     call printmtx(matice)
!     
! !     call printmtx(Atb)
!     
!     call gem(matice, Atb, x)
!   end subroutine sparse_gem_pig_AtA
!   


end module debug_tools