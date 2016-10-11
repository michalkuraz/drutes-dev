!> \file schwarz.f90
!! \brief implementace Schwarzovske domain decomposition
!!
!! podrobnejsi popis
!!
!<

!> implementace sch
!!
!! dalsi text
module schwarz
   public :: schwarz_driver
  ! public :: schwarz_setup
contains

   !> hlavni driver pro Schwarze
   !!
   !! zavola konstruktor matic, a zkonfiguruje 
   !! Schwarzovskou DD-metodu.
   !! - zjisti pocet domen
   !! - rozhodne o pouziti symetrizacniho kroku
   !! - rozhodne o konstrukci hrube ulohy. 
   !!   Pokud neni, tak zjisti velikost prekryvu
   !! - rozhodne o resicich pro jednotlive kroky   
   subroutine schwarz_driver
      use defs
      use mtx
      use uloha
      
      implicit none
      logical :: repeat
      type(Mtx_Type) :: A
      real(kind=rkind), dimension(:), pointer :: xright, b
       
      print *,"cau"
      call init(A,20_ikind)
      
      !ted skutecny driver
      repeat = .true.
      do while (repeat)
         print *, " Schwarzova domain decomposition"
         !1. vyber resenou matici
         call Set_Problem(A,b,xright)
         repeat = .false.
         call Print_Mtx(A)
      end do
   end subroutine schwarz_driver
end module schwarz
