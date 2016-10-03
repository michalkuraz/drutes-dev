! Copyright 2008 Michal Kuraz, Petr Mayer

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.


module fem_tools
  public :: in2global
  public :: icond4neumann


  contains

  !> procedure to put local stiffness matrix into the global stiffness matrix
  !<
  subroutine in2global(el_id, locmatrix, bvect, subprmt, debugme)
    use typy
    use sparsematrix
    use global_objs
    use pde_objs
    use globals
    use debug_tools
    
    !> number of element
    integer(kind=ikind), intent(in) :: el_id
    class(extsmtx), intent(in out)  :: locmatrix
    real(kind=rkind), dimension(:) :: bvect
    !> permutation array for submatrices, supplied only for domain decomposition problems
    integer(kind=ikind), dimension(:), intent(in), optional :: subprmt
    logical, intent(in), optional :: debugme
    
    integer(kind=ikind), dimension(:), allocatable, save :: bc
    integer(kind=ikind), dimension(:), allocatable, save :: n_row, m_col
    integer(kind=ikind) :: i,j,m, iproc, jproc, limits
    real(kind=rkind), dimension(:), allocatable, save :: bcval
    real(kind=rkind), dimension(:), allocatable, save :: surface
    integer(kind=ikind), dimension(:), allocatable, save :: fin
    real(kind=rkind), dimension(3,3) :: d
   

    if (.not. allocated(bc)) then	
      allocate(bc(ubound(stiff_mat,1)))
      allocate(n_row(ubound(stiff_mat,1)))
      allocate(m_col(ubound(stiff_mat,1)))
      allocate(bcval(ubound(stiff_mat,1)))
      allocate(surface(ubound(stiff_mat,1)))
      allocate(fin(ubound(pde,1)))
    end if
    
    surface = 0
    bc = 0
    n_row = 0
    m_col = 0
    bcval = 0
    limits = ubound(stiff_mat,1)/ubound(pde,1)
        

    do iproc=0, ubound(pde,1)-1
      do i=1, ubound(elements%data, 2) 	
	surface(i+iproc*limits) = elements%length(el_id,i)
      end do
    end do
    
    
    do iproc = 0,ubound(pde,1)-1
      do i=1, ubound(elements%data, 2)
        j = elements%data(el_id,i)
        n_row(i+iproc*limits) = pde(iproc+1)%permut(j)
        if (nodes%edge(j) > 100) then
          m = nodes%edge(j)
          call pde(iproc+1)%bc(m)%value_fnc(pde(iproc+1), el_id,i,bcval(i+iproc*limits),bc(i+iproc*limits))
        else
          bc(i) = 0
        end if
      end do
    end do
    
    

    

    if (present(subprmt)) then
      do i=1, ubound(m_col,1)
        if (n_row(i) > 0 ) then
          m_col(i) = subprmt(n_row(i))
        end if
      end do
    else
      m_col=n_row
    end if

 

    do i=1, ubound(stiff_mat,1)
      ! fill bside
      if (abs(bc(i)) /= 1) then
	select case(bc(i))
	  case(0,3)
              if (m_col(i) > 0) then
                bvect(m_col(i)) = bvect(m_col(i)) + bside(i)
              end if
	  case(2)
              if (m_col(i) > 0) then
                bvect(m_col(i)) = bvect(m_col(i)) - bcval(i)*surface(i)*time_step + bside(i)
              end if
	  case(4)
	    print *, "seepage face boundary not yet implemented"
	    print *, "the program was interupted from fem_tools::in2global procedure"
	    ERROR STOP
	end select
	! fill stiffness matrix
	do m=1, ubound(stiff_mat,1)
	  select case(bc(m))
	    case(-1,1)
                if (m_col(i) > 0) then
                  bvect(m_col(i)) = bvect(m_col(i)) - stiff_mat(i,m)*bcval(m)
                end if
	    case default

                if (n_row(i) > 0 .and. m_col(m) > 0) then
                  call locmatrix%add(stiff_mat(i,m), n_row(i), m_col(m))
                  if (drutes_config%it_method == 2 .or. drutes_config%it_method == 1) then
                    call locmatrix%rowsfilled%nrfill(n_row(i))
                  end if
                end if
                

              !**!
! 	      spmatrix%vals(g_row) = stiff_mat(i,m)
! 	      spmatrix%ii(g_row) = n(i)
! 	      spmatrix%jj(g_row) = n(m)
! 	      g_row = g_row + 1
	  end select
	end do
      else
	CONTINUE
      end if
    end do
	  

  
  end subroutine in2global




  subroutine icond4neumann(proc)
    use typy
    use globals
    use global_objs
    use pde_objs
    use core_tools
    use simplelinalg

    integer(kind=ikind), intent(in) :: proc
    integer(kind=ikind) :: i,j,k,l,m,n,D, o, kk, nn, code, z
   
    logical, dimension(:), allocatable, save :: evaluated
    real(kind=rkind) :: tmp, length, der
    real(kind=rkind), dimension(3) :: gradh, velocity
    real(kind=rkind), dimension(3,3) :: disp
    logical :: gofurther

    if (.not. allocated(evaluated)) then
      allocate(evaluated(nodes%kolik))
    end if

    evaluated = .false.

    D = drutes_config%dimen

    do i=1, elements%kolik
      do j=1, ubound(elements%data,2)
        k = elements%data(i,j)
        l = nodes%edge(k)
        if (.not.evaluated(k) .and. elements%length(i,j) > epsilon(tmp) .and. l >= lbound(pde(proc)%bc,1)) then
          if (pde(proc)%bc(l)%code == 2 .and. pde(proc)%bc(l)%icond4neumann) then
            do m=1, ubound(elements%data,2)
              o = nodes%edge(elements%data(i,m))
              if (o >= lbound(pde(proc)%bc,1)) then
                if (pde(proc)%bc(o)%code /= 2) then
                  gofurther = .true.
                else
                  gofurther = .false.
                end if
              else
                gofurther = .true.
              end if
              if (gofurther) then
		call pde(proc)%bc(l)%value_fnc(pde(proc), i, j, tmp, code)
                select case(drutes_config%dimen)
                    case(1)
                      velocity(1) = elements%nvect_z(i,j)*tmp
                    case(2)
                      velocity(1) = sqrt(1-elements%nvect_z(i,j)*elements%nvect_z(i,j))*tmp
                      velocity(2) = elements%nvect_z(i,j)*tmp
                end select

		kk = pde(proc)%permut(k)

                call pde(proc)%pde_fnc(proc)%dispersion(pde(proc), &
                elements%material(i,proc),x=(/pde_common%xvect(kk,1)/), & 
                                                          tensor=disp(1:D, 1:D))

		do z=1, D
		  if ( disp(z,z) < 10*epsilon(tmp)) then
		    call write_log("W: dispersion is nearly zero")
		    RETURN
		  end if
		end do

                call invert_matrix(disp(1:D, 1:D))

                gradh(1:D) = matmul(disp(1:D, 1:D), velocity(1:D))

                EXIT

              end if
            end do

            if (m==4) then
              print *, "ERROR: possibly strange mesh, contact developer, called from fem_tools::icond4neumann"
              print *, "element id is: ", i
              ERROR STOP
            end if
            n = elements%data(i,m)

	    nn = pde(proc)%permut(n)
            pde_common%xvect(kk,1) = pde_common%xvect(nn,1)
            do o=1,D
              pde_common%xvect(kk,1) = pde_common%xvect(kk,1) - (nodes%data(k,o)-nodes%data(n,o))*gradh(o)
            end do
            evaluated(k) = .true.
          end if
        end if
      end do
    end do

  end subroutine icond4neumann
  


end module fem_tools