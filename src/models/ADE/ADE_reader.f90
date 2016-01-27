module ADE_reader
  public :: ADEs_read, ADEd_read

  subroutine ADEs_read()
    use typy
    use globals
    use global_objs
    use core_tools

    integer :: ierr

    call find_unit(file_contaminantm, 200)
    open(unit = file_contaminantm, file="drutes.conf/contaminant.conf/matrix.conf", action="read", status="old", iostat=i_err)
    if (ierr /= 0) then
      print *, "missing drutes.conf/contaminant.conf/matrix.conf file"
      ERROR STOP
    end if

!       use typy
!       use globals
!       use core_tools
! 
!       integer :: ierr, n, i
! 
!       allocate(ade_par(ubound(vgmatrix,1)))
! 
!       call comment(file_contaminantm)
! 
!       read(unit=file_contaminantm, fmt=*, iostat=ierr) isotherm ; if (ierr /= 0) call file_error(file_waterf)
! 
!       call readADEvals(file_contaminantm, ade_par)
! 
!       select case(problem_type)
!       case(100,200)
!           call readbcvals(file_contaminantm, ADE_single%bc, dimen= (ubound(RICHARDS_single%bc,1)-101), &
!                         dirname="drutes.conf/contaminant.conf/")
!       case(110,210)
! !         call readbcvals(file_contaminantm, ADE_dual(1,1)%bc, dimen = (ubound(RICHARDS_dual(1,1)%bc,1)-101),  &
! !                       dirname="drutes.conf/contaminant.conf/")
!       end select

  end subroutine ADEs_read


  subroutine ADEd_read()
    use typy
    use globals
    use global_objs
    use core_tools

    integer :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !contaminant.conf/fractures.conf
    call find_unit(file_contaminantf, 200)
    open(unit=file_contaminantf, file="drutes.conf/contaminant.conf/fractures.conf", action="read", status="old", iostat = i_err)
    if (i_err /= 0) then
      print *, "missing drutes.conf/contaminant.conf/fractures.conf file"
      ERROR STOP
    end if


!       integer :: ierr, n, i
! 
!       allocate(ade_par_f(ubound(vgfractures,1)))
! 
!       call readADEvals(file_contaminantf, ade_par_f)
! 
!       call readbcvals(file_contaminantm, ADE_dual%bc_sys2, dimen=(ubound(RICHARDS_dual%bc_sys1,1)-101),  &
!                         dirname="drutes.conf/contaminant.conf/")
! 
!       call readbcvals(file_contaminantf, ADE_dual%bc_sys2, dimen=(ubound(RICHARDS_dual%bc_sys1,1)-101),  &
!                         dirname="drutes.conf/contaminant.conf/")

  end subroutine ADEd_read
  
    subroutine readADEvals(unitC, struct)
      use typy
      use globals
      use core_tools
      integer, intent(in) :: unitC
      type(soluteXsoil), dimension(:), intent(out) :: struct
      integer(kind=ikind) :: i
      integer :: i_err
      character(len=256) :: filename

      do i=1, ubound(struct,1)
        call comment(unitC)
        read(unit=unitC, fmt=*, iostat =i_err) struct(i)%Dz, struct(i)%Dx, struct(i)%diff, struct(i)%kd, struct(i)%expo, &
                struct(i)%rho, struct(i)%bd, struct(i)%coop, struct(i)%T_w, struct(i)%T_s, struct(i)%icond
        if (i_err /= 0) then
          inquire(unit=unitC, name=filename)
          print *, "!!!!!!!!!!!!!!!!!!!!"
          print *, "HINT: check number of line records and parameters in: ", trim(filename)
          print *, "!!!!!!!!!!!!!!!!!!!!"
          call file_error(unitC)
        end if
      end do

    end subroutine readADEvals

end module ADE_reader