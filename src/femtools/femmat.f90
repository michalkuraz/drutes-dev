! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher


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



!> \file femmat.f90
!! \brief Local FEM matrix assembler.
!<



module femmat
  public :: solve_picard
  private :: results_extractor
  private :: results_eraser
  public :: assemble_mat
  

  contains
  !>solves the problem using Picard method
    subroutine solve_picard(ierr, success)
      use typy
      use globals
      use global_objs
      use pde_objs
      use RE_constitutive
      use feminittools
      use fem_tools
      use linAlg
      use debug_tools
      use solvers
      use postpro
      use decomp_tools
      use decomp_vars
      use simplelinalg
      use debug_tools
      use mtxiotools
      use core_tools
      use solver_interfaces
      use global4solver

!       use decomposer


      integer, intent(out) :: ierr
      integer(kind=ikind) ::  fin, pcg_it
      logical, intent(out) :: success
      integer(kind=ikind) :: i, proc,j, m, l, k, top, bot
      real(kind=rkind) :: error, loc_error, reps_err, l1, l2
      logical :: dt_fine
      integer :: ierr_loc, ll
      real(kind=rkind), dimension(:), allocatable :: vcttmp
      real(kind=rkind) :: lambda_l, lambda_h, tmpxx=0, maxtime=0, solver_start, solver_end, picard_old
      
      integer, save :: pocitac=0
      character(len=100) :: soubor

      proc = ubound(pde,1)
      fin = maxval(pde(proc)%permut(:))
      itcount = 0
      allocate(vcttmp(ubound(pde_common%bvect,1)))
      
      picard_old = huge(picard_old)
      
      
      do

        itcount = itcount + 1

    
        call assemble_mat(ierr)

        if (record_solver_time) call cpu_time(solver_start)
        
        solver_error = .false.
        call solve_matrix(spmatrix, pde_common%bvect(1:fin), pde_common%xvect(1:fin,3),  itmax1=fin, &
            reps1=1e-14_rkind, itfin1=pcg_it, repsfin1=reps_err)
            
        if (solver_error) then
          ierr = 1
          success = .false.
          call write_log("Looks like a poor conditioning...:/")
          EXIT
        end if
        
        
        if (itersolver) then
          write(unit=file_itcg, fmt=*) time, pcg_it, reps_err
          call flush(file_itcg)
        end if
            
        if (record_solver_time) then
          call cpu_time(solver_end)
          write(unit=solver_time_file, fmt=*) time, solver_end - solver_start
          call flush(solver_time_file)
        end if
          
  
        if (pcg_it > 0.55*fin) then 
          ierr=-1
        else
          ierr=0
        end if
		  
    
        error = norm2(pde_common%xvect(1:fin,2)-pde_common%xvect(1:fin,3))/ubound(pde_common%xvect,1)
        
        
        write(unit=terminal, fmt=*) "Nonlinear solver convergence criterion:", error
        
        write(unit=file_picard, fmt=*) time, error
   
        if (itcount == 1 .or. error <= iter_criterion) then
          do proc=1, ubound(pde,1)
            top = pde(proc)%permut(1)
            bot = pde(proc)%permut(ubound(pde(proc)%permut,1))
            dt_fine = pde(proc)%dt_check()
          end do
        else
          dt_fine = .true.
        end if


        if (.not.(dt_fine)) then
          ierr = 2
          success = .false.
          RETURN
        end if
        
        if (itcount == 1) picard_old = error
        
        if (itcount >= 3) then
          if (error > picard_old) then
            ierr = 2
            success = .false.
            call results_eraser()
            RETURN
          end if
        end if
        
        if (isnan(error)) then
          ierr = 2
          success = .false.
          call results_eraser()
          RETURN
        end if
				



        if (error <= iter_criterion) then

          if (ierr /= -1) ierr = 0
          
          if (drutes_config%check4mass) then
            call make_print("separately")
            call do_masscheck()
            
          end if

          call results_extractor()

          success = .true.
          iter_succesfull = .true.
          RETURN
        else
          pde_common%xvect(:,2) = pde_common%xvect(:,3)
        end if

        if (itcount >= max_itcount) then
          ierr = 1
          success = .false.
          iter_succesfull = .false.
          call results_eraser()
          EXIT
        end if
      end do

    end subroutine solve_picard
    
    subroutine assemble_mat(ierr)
      use typy
      use globals
      use global_objs
      use pde_objs
      use capmat
      use stiffmat
      use feminittools
      use geom_tools
      use fem_tools
      use re_constitutive
      use linAlg
      use solver_interfaces
      use debug_tools

      integer, intent(out) :: ierr

      integer(kind=ikind) :: i,j,k,l, proc, ll, limits
            
      type(integpnt_str) :: quadpnt
      
      real(kind=rkind) :: value

      call null_problem(spmatrix, pde_common%bvect)
      
      limits = ubound(stiff_mat,1)/ubound(pde,1)

      do i=1, elements%kolik
        pde_common%current_el = i
        processes: do proc=1, ubound(pde,1)
		    do j=1+(proc-1)*limits, ubound(elements%data,2) + (proc-1)*limits
		      ll = j - (proc-1)*ubound(stiff_mat,1)/ubound(pde,1)
		      k = pde(proc)%permut(elements%data(i,ll))
		      if (k > 0) then
            elnode_prev(j) = pde_common%xvect(k,1)
		      else
            k = nodes%edge(elements%data(i,ll))
            call pde(proc)%bc(k)%value_fnc(pde(proc), i, ll, value)
            elnode_prev(j) = value
		      end if
		    end do
        end do processes

        call build_bvect(i, time_step)

        call build_stiff_np(i, time_step)

        call pde_common%time_integ(i)
 

        stiff_mat = stiff_mat + cap_mat 

        call in2global(i,spmatrix, pde_common%bvect)

      end do
      printedk=.false.

    end subroutine assemble_mat
    
    
    subroutine results_extractor()
      use typy
      use globals
      use global_objs
      use pde_objs

      integer(kind=ikind) :: i, j, proc, nd, edge
      real(kind=rkind) :: value

      
      do proc=1, ubound(pde,1)
        do i=1, elements%kolik
          do j=1, ubound(elements%data,2)
            nd = elements%data(i,j)
            if (pde(proc)%permut(nd) > 0) then
              pde(proc)%solution(nd) = pde_common%xvect(pde(proc)%permut(nd),3)
            else
              edge = nodes%edge(nd)
              call pde(proc)%bc(edge)%value_fnc(pde(proc), i, j, value)
              pde(proc)%solution(nd) = value
            end if
          end do
        end do
      end do

      pde_common%xvect(:,4) = pde_common%xvect(:,1)
      pde_common%xvect(:,1) = pde_common%xvect(:,3)
      pde_common%xvect(:,2) = pde_common%xvect(:,3)

    end subroutine results_extractor
    
    subroutine results_eraser()
      use pde_objs
    
      pde_common%xvect(:,2) = pde_common%xvect(:,1)
      pde_common%xvect(:,3) = 0
    
    end subroutine results_eraser
    


end module femmat
