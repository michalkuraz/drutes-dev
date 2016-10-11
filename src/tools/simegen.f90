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


      !> vstup pozadovany krok deleni
      integer(kind=ikind) :: pocet1
      integer(kind=ikind), dimension(:), allocatable :: pocty, sumy
      real(kind=rkind) :: dx1
      integer(kind=ikind) :: i,j, last, n


      allocate(pocty(ubound(deltax_1d,1)))
      allocate(sumy(0:ubound(deltax_1d,1)))
      
      if (abs(length_1D-(deltax_1d(ubound(deltax_1d,1), 3)-deltax_1d(1,2))) > 10*epsilon(dx1)) then
       write(unit=terminal, fmt=*) " "
       write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "!!ERROR!!" //achar(27)//'[0m'
        write(unit=terminal, fmt=*)"incorrect definition in either:"
        write(unit=terminal, fmt=*)"                          - domain length definition"
        write(unit=terminal, fmt=*)"or"
        write(unit=terminal, fmt=*)"                          - the mesh description values"
        write(unit=terminal, fmt=*)" "
        write(unit=terminal, fmt=*)"in file drutes.conf/mesh/drumesh1d.conf"
        write(unit=terminal, fmt=*)" "
        write(unit=terminal, fmt=*)"the mesh description must cover (and cannot overlap) the entire domain length"
        write(unit=terminal, fmt=*)"-----------------------------------------------------------------------------"
        ERROR STOP
      end if


      if (abs(length_1D-(materials_1d(ubound(materials_1d,1), 2)-deltax_1d(1,2))) > 10*epsilon(dx1) .or. &
         abs(materials_1d(1,1) - deltax_1d(1,2)) >  10*epsilon(dx1) ) then
         
        write(unit=terminal, fmt=*) " "
        write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "!!ERROR!!" //achar(27)//'[0m' 
        write(unit=terminal, fmt=*)"incorrect definition in either:"
        write(unit=terminal, fmt=*)"                          - domain length definition"
        write(unit=terminal, fmt=*)"or"
        write(unit=terminal, fmt=*)"                          - the description of the material distribution"
        write(unit=terminal, fmt=*)" "
        write(unit=terminal, fmt=*)"in file drutes.conf/mesh/drumesh1d.conf"
        write(unit=terminal, fmt=*)" "
        write(unit=terminal, fmt=*)"the description of the material distribution must cover (and cannot overlap) &
	  the entire domain length"
        write(unit=terminal, fmt=*)"-----------------------------------------------------------------------------"
        ERROR STOP
      end if

     
      pocet1 = 0
      do i = 1, ubound(deltax_1d,1)
	if (abs((int((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1)) - &
		  ((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1)))) > 100*epsilon(dx1)) then
	  n = int((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1)) + 1
	else
	  n = int((deltax_1d(i,3) - deltax_1d(i,2))/deltax_1d(i,1))
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
      sumy(0) = 0
      deltax_1d(0,:) = 0
      deltax_1d(0,3) = deltax_1d(1,2)

      

	
      do i=1, ubound(pocty,1)
	do j=1, pocty(i) 
	  if (i == 1) then
	    nodes%data(j+sumy(i-1),1) = deltax_1d(i-1,3) + (j-1)*deltax_1d(i,1)
	  else
	    nodes%data(j+sumy(i-1),1) = deltax_1d(i-1,3) + j*deltax_1d(i,1)
	  end if
	end do
      end do
	
      elements%data(1,1) = 1
      elements%data(1,2) = 2
      do i=2, elements%kolik
	elements%data(i,1) = elements%data(i-1,2)
	elements%data(i,2) = i+1
      end do



      do i=1, elements%kolik
	do j=1,ubound(materials_1D,1)
	  if (avg(nodes%data(elements%data(i,1),1), nodes%data(elements%data(i,2),1)) > materials_1D(j,1) .and. &
	      avg(nodes%data(elements%data(i,1),1), nodes%data(elements%data(i,2),1)) < materials_1D(j,2)) then
	      
	      elements%material(i, :) = j
	      EXIT
	  end if
	end do
      end do
      if (minval(elements%material) == 0) then
	print *, "check file drutes.conf/mesh/drumesh1d.conf, seems like the material description does not cover the entire domain"
	ERROR STOP
      end if

      nodes%edge = 0

      nodes%edge(1) = 101

      nodes%edge(ubound(nodes%edge,1)) = 102
  
      


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

!       do i=1, elements%kolik
! 	if (nodes%data(elements%data(i,1),2) > 50) then
! 	  elements%material(i,:) = 1
! 	else
! 	  elements%material(i,:) = 2
! 	end if
!       end do
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
	



! 	nodes%data(2,:) = (/0,1/)
! 	nodes%data(3,:) = (/0,2/)
! 	nodes%data(4,:) = (/1,0/)
! 	nodes%data(5,:) = (/1,1/)
! 	nodes%data(6,:) = (/1,2/)
! 
! 	elements%data(1,:) = (/1,4,5/)
! 	elements%data(2,:) = (/1,5,2/)
! 	elements%data(3,:) = (/2,5,6/)
! 	elements%data(4,:) = (/2,6,3/)
! 
! 	nodes%edge = 0
! 	nodes%edge(1) = 101
! 	nodes%edge(4) = 101


! 	if (num_images() > 1) then
! 	  allocate(cpus(num_images(),2))
! 	  i = elements%kolik/num_images()
! 	  
! 	  do j=1, num_images()-1
! 	    cpus(j,1) = (j-1)*i+1
! 	    cpus(j,2) = j*i
! 	  end do
!     
! 	  
! 	  cpus(num_images(),1) = cpus(num_images()-1,2)+1
! 	  cpus(num_images(),2) = elements%kolik
!     
! 	  i=1
! 	  do j=1, elements%kolik
! 	    do
! 	      if (j >= cpus(i,1) .and. j <= cpus(i,2)) then
! 		elements%CPU(j) = i
! 		EXIT
! 	      else
! 		i = i+1
! 	      end if
! 	    end do
! 	  end do
!     
! 	  deallocate(cpus)
! 	else
! 	  elements%CPU = 1
! 	end if

    
  !     do i=1, nodes%kolik
  !       do j=lbound(waterbc_m,1), lbound(waterbc_m,1)
  ! 	do k=1, ubound(waterbc_m(j)%x,1)-1
  ! 	  a = (/waterbc_m(j)%x(k), waterbc_m(j)%y(k)/)
  ! 	  b = (/waterbc_m(j)%x(k+1), waterbc_m(j)%y(k+1)/)
  ! 	  	  print *, a,b
  ! 	  print *, nodes%data(i,:)
  ! 	  if (inline(a,b,nodes%data(i,:))) then
  ! 	      nodes%edge(i) = waterbc_m(j)%ID
  ! 	      EXIT
  ! 	  end if
  ! 	end do
  !       end do
  !     end do
		
	      
    
  !     do i=1, nodes%kolik
  !       if ( abs(nodes%data(i,2) - ly) < reps ) then 
  !         nodes%bc_c(i, :) = 1_ikind
  !       else
  !         nodes%bc_c(i,:) = 0_ikind
  !       end if
  !     end do
  !     
  !   



  !     do i=1, ubound(nodes%data,1)
  !       print *, i, nodes%data(i,:), "edge", nodes%edge(i)
  !     end do
  !   stop
  !     stop
  !     do i=1, nodes%kolik
  !       print *, i, nodes%data(i,:), nodes%bc(i,:)
  !     end do
  ! !     
  !     do i=1, elements%kolik
  !       print *, i, elements%data(i,:)
  !     end do
  !   stop 

  end subroutine simegen2D

end module simegen