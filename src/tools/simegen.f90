
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

!> \file simegen.f90
!! \brief Simple mesch generators
!<

!> 
!! Mesh generator for 1D - handles nonuniform discretization. For 2D creates uniform mesh, good for debugging, useless for real world computations.
!<


!>simple mesh generator
module simegen
  public :: simegen1D
  public :: simegen2D
  
  contains
  
    subroutine simegen1D()
      use typy
      use global_objs
      use globals
      use globals1d
      use core_tools
      use readtools
      use geom_tools
      use debug_tools


      integer(kind=ikind) :: pocet1
      integer(kind=ikind), dimension(:), allocatable :: pocty, sumy
      real(kind=rkind) :: dx1
      integer(kind=ikind) :: i,j, last, n
      logical, dimension(:), allocatable :: testmat
      logical :: err_mat = .false.
      real(kind=rkind), dimension(:), allocatable :: nd_tmp
      real(kind=rkind), dimension(2) :: A,B,pt


      allocate(pocty(ubound(deltax_1d,1)))
      allocate(sumy(0:ubound(deltax_1d,1)))
      

     
      pocet1 = 0
      do i = 1, ubound(deltax_1d,1)
        if (abs((int((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1)) - &
            ((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1)))) > 100*epsilon(dx1)) then
          n = abs(int((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1))) + 1
        else
          n = abs(int((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1)))
        end if    

        deltax_1d(i,1) = (deltax_1d(i,3) - deltax_1d(i,2))/(1.0_rkind*n)

        if (i == 1) then
          n = n + 1
        end if
        pocty(i) = n
        pocet1 = pocet1 + n
        sumy(i) = pocet1
      end do
      

      nodes%kolik=pocet1
      elements%kolik = nodes%kolik - 1
      call mesh_allocater()
      allocate(nd_tmp(nodes%kolik))
      
      sumy(0) = 0
      deltax_1d(0,:) = 0
      deltax_1d(0,3) = deltax_1d(1,2)

      
	
      do i=1, ubound(pocty,1)
        do j=1, pocty(i) 
          if (i == 1) then
            nd_tmp(j+sumy(i-1)) =  deltax_1d(i-1,3) + (j-1)*deltax_1d(i,1)
!            nodes%data(j+sumy(i-1),1) = deltax_1d(i-1,3) + (j-1)*deltax_1d(i,1)
          else
!            nodes%data(j+sumy(i-1),1) = deltax_1d(i-1,3) + j*deltax_1d(i,1)
             nd_tmp(j+sumy(i-1)) = deltax_1d(i-1,3) + j*deltax_1d(i,1)
          end if
        end do
      end do
      
      if (nd_tmp(1) < nd_tmp(2)) then
        nodes%data(:,1) = nd_tmp
      else
        do i=1,nodes%kolik
          nodes%data(i,1) = nd_tmp(nodes%kolik-i+1)
        end do
      end if
      
	
      elements%data(1,1) = 1 
      elements%data(1,2) = 2
      do i=2, elements%kolik
        elements%data(i,1) = elements%data(i-1,2)
        elements%data(i,2) = i+1
      end do
      
      A=0
      B=0
      pt=0
      do i=1, elements%kolik
        do j=1,ubound(materials_1D,1)
          pt(1) = avg(nodes%data(elements%data(i,1),1), nodes%data(elements%data(i,2),1))
          A(1) = materials_1D(j,2)
          B(1) =  materials_1D(j,3)
          if (inline(A,B,pt)) then
              elements%material(i) =  nint(materials_1D(j,1))
              EXIT
          end if
        end do
      end do

      
      if (minval(elements%material) == 0) then
        print *, "check file drutes.conf/mesh/drumesh1d.conf, seems like the material description does not cover the entire domain"
        ERROR STOP
      end if
      
      allocate(testmat(nint(maxval(materials_1D(:,1)))))
      
      testmat = .false.
      
      do i=1, ubound(elements%material,1)
        testmat(elements%material(i)) = .true.
      end do
      
      do i=1, ubound(testmat,1)
        if (.not. testmat(i)) then
          print *, "ERROR! material with ID", i, "undefined in drutes.conf/mesh.conf/drumesh1d.conf"
          err_mat = .true.
        end if
      end do
      
      if (err_mat) call file_error(file_mesh, "ERROR!! check your material description")
      

      nodes%edge = 0
      
      if (nodes%data(1,1) < nodes%data(nodes%kolik,1) ) then
        nodes%edge(1) = 101
        nodes%edge(ubound(nodes%edge,1)) = 102
      else
        nodes%edge(1) = 102
        nodes%edge(ubound(nodes%edge,1)) = 101
      end if


    end subroutine simegen1D

    subroutine  simegen2D(lx, ly, density, edges_xy, center_coor)
      use typy
      use globals
      use globals2D
      use geom_tools
      use core_tools
      use debug_tools

      !> mesh parameters
      real(kind=rkind), intent(in) :: lx, ly, density
      !> description of the mesh edges
      real(kind=rkind), dimension(:,:,:), intent(in) :: edges_xy
      real(kind=rkind), dimension(:), intent(in) :: center_coor
      real(kind=rkind) :: step1, step2, tmp1, tmp2, reps
      real(kind=rkind), dimension(2) :: a,b
      
      integer :: i,j, counter, m, count1, count2, k, l
      real(kind=rkind), dimension(2) :: c
      integer(kind=ikind), dimension(:,:), allocatable :: cpus
      
      

      count1 = int(lx/density)
      
      step1 = lx/count1
      
      count2 = int(ly/density)
      
      step2 = ly/count2
      
      nodes%kolik = (count1+1)*(count2+1)
      elements%kolik = 2*count1*count2
      

      call mesh_allocater()

      counter = 1_ikind
      
      do i=0, count2
        do j=0, count1
          nodes%data(counter,1) = j*step1
          nodes%data(counter,2) = i*step2
          counter = counter + 1
        end do
      end do
      

      j = 1_ikind
      m = 0_ikind
      do i=1, elements%kolik - 1, 2
        if (m >= count1) then
          j = j+1
          m = 0_ikind
        end if
        m = m + 1
        elements%data(i,1) = int(i/2) + j
        elements%data(i,2) = int(i/2) + j+1 + count1
        elements%data(i,3) = int(i/2) + j+1 + count1 + 1
        elements%data(i+1,1) = int(i/2) + j
        elements%data(i+1,2) = int(i/2) + j+1
        elements%data(i+1,3) = int(i/2) + j+1 + count1 + 1
      end do


      elements%material = 1
      
      do i=1, nodes%kolik
        nodes%data(i,:) = nodes%data(i,:) + center_coor
      end do
      

      reps = epsilon(reps)
      

      nodes%edge = 0_ikind
      
    

      do i=1, nodes%kolik
        do j=1, ubound(edges_xy,1)
          a = edges_xy(j,1,:)
          b = edges_xy(j,2,:)	  

          c = nodes%data(i,:)

          if (inline(a,b,c)) then
            nodes%edge(i) = j+100
          end if
        end do
      end do
	


  end subroutine simegen2D

end module simegen
