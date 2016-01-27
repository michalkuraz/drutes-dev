module type_module
  implicit none

  type test_type
    procedure(fun_interface), nopass, pointer :: fun_ptr
    procedure(sub_interface), nopass, pointer :: sub_ptr    
  end type test_type

  abstract interface
    function fun_interface(n,x) result(f)
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision :: f(n)
    end function fun_interface
    subroutine  sub_interface(n,x,f)
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision :: f(n)
    end subroutine sub_interface
  end interface  

  contains

  subroutine  test_type_constructor(test, fun,sub)    
    interface
      function fun(n,x) result(f)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision :: f(n)
      end function fun
      subroutine  sub(n,x,f)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision :: f(n)
      end subroutine sub
    end interface  

    type(test_type) :: test
    test%fun_ptr => fun;
    test%sub_ptr => sub;
  end subroutine test_type_constructor
end module type_module
 
module funcs
  implicit none

  contains

  function fun1 (n,x ) result (f )
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision :: f(n)
    ! try this
    f = 2.0*x ;
    ! or this , both give WRONG results
    call  sub1(n,x,f)
  end function fun1

  subroutine sub1(n,x,f)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision :: f(n)
    f = 2.0*x ;
  end subroutine sub1
end module funcs

program main
  use type_module
  use funcs
  type(test_type) :: test
  integer :: n =2 ;
  double precision :: x(2), f(2)
  call test_type_constructor (test, fun1, sub1)
  x = (/-1.d0, 1.d0/);
  f = 0.d0;
  call  test%sub_ptr(n,x,f)
  print *, " f from call fun ", f;

  f = test%fun_ptr(n,x);
  print *, " f from call fun ", f;
  ! This seems to work
  print *, " size of returned value from fun_ptr ", size(test%fun_ptr(n,x));
end program