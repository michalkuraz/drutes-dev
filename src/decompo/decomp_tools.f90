module decomp_tools
  public :: print_domains, print_elements_dd
  public :: get_residual
  public :: collect_matrices
  public :: addlevel
  public :: init_addlevel
  public :: build_coarsemtx
  public :: mulAx_dd
  private :: add_el
  public :: get_residual_old
  public :: getres_loc
   

  contains
  
    subroutine getres_loc(subdom, debugme)
      use typy 
      use globals
      use global_objs
      use pde_objs
      use decomp_vars
      use debug_tools
      
      type(subdomain_str), intent(in out) :: subdom
      logical, intent(in), optional :: debugme
      integer(kind=ikind) :: i, nd, domain_id, ndloc, ndpt
            
      real(kind=rkind), dimension(:), allocatable :: Arow
      logical, dimension(:), allocatable, save :: coltrue
      logical, dimension(:), allocatable, save :: extcoltrue
      real(kind=rkind), dimension(:), allocatable, save :: colvals, extcolvals
      integer(kind=ikind), dimension(:), allocatable, save :: jjvals, extjjvals
      integer(kind=ikind) :: numbers, row, rowperm, extrowperm, j
      
      
      
      subdom%resvct%main=subdom%bvect
      subdom%resvct%ext=subdom%extbvect
      

      ! build resvct main
      do i=1, subdom%matrix%rowsfilled%pos
        row = subdom%matrix%rowsfilled%data(i)
        

        
        if (subdom%invpermut(row) /= 0) then
          rowperm=subdom%invpermut(row)
          
          call subdom%matrix%getrow(i=row, v=colvals, jj=jjvals, nelem=numbers)
      
          
          subdom%resvct%main(rowperm) = subdom%resvct%main(rowperm) - dot_product(colvals(1:numbers), &
            subdom%xvect(jjvals(1:numbers), 2))
            
            
          call subdom%extmatrix%getrow(i=row, v=colvals, jj=jjvals, nelem=numbers)
          
          if (numbers > 0) then
	    subdom%resvct%main(rowperm) = subdom%resvct%main(rowperm) - dot_product(colvals(1:numbers), &
            subdom%extxvect(jjvals(1:numbers), 2))
          end if
          

        end if   
        
      end do
      
      ! build resvct ext  
      do i=1, subdom%extmatrix%rowsfilled%pos
        row = subdom%extmatrix%rowsfilled%data(i)
        
        if (subdom%extinvpermut(row) /= 0) then
          rowperm=subdom%extinvpermut(row)
          
          call subdom%extmatrix%getrow(i=row, v=colvals, jj=jjvals, nelem=numbers)      
                    
          subdom%resvct%ext(rowperm) = subdom%resvct%ext(rowperm) - dot_product(colvals(1:numbers), &
            subdom%extxvect(jjvals(1:numbers), 2))
            
          call subdom%matrix%getrow(i=row, v=colvals, jj=jjvals, nelem=numbers) 
          
          
          if (numbers > 0) then
	    subdom%resvct%ext(rowperm) = subdom%resvct%ext(rowperm) - dot_product(colvals(1:numbers), &
            subdom%xvect(jjvals(1:numbers), 2))
          end if         
        end if 
        
      end do    
      
     
            
    end subroutine getres_loc
  
    subroutine mulAx_dd(x,vct)
      use typy 
      use globals
      use global_objs
      use decomp_vars
      
      real(kind=rkind), dimension(:), intent(in) :: x
      real(kind=rkind), dimension(:), intent(out) :: vct
      integer(kind=ikind) :: i, j, fin, posglob, posloc, subfin
      real(kind=rkind), dimension(:), allocatable :: vectorin

      logical :: skip_solved=.false.
      
      vct = 0.0_rkind
      do i=1, ubound(subdomain,1)
        !loop for A.x
        if ( (.not. skip_solved) .or.  (.not. subdomain(i)%solved) ) then
          allocate(vectorin(subdomain(i)%ndof))
          vectorin = 0
          fin = subdomain(i)%disjoint%pos
          do j=1, fin
            posglob = subdomain(i)%disjoint%data(j)
            posloc = subdomain(i)%invpermut(posglob)
            vectorin(posloc) = x(posglob)
          end do
          vct = vct +  subdomain(i)%matrix%mul(vectorin)
          deallocate(vectorin)
        end if
      end do

    end subroutine mulAx_dd
    
    subroutine get_residual(resvct)
      use typy
      use global_objs
      use pde_objs
      use globals
      use decomp_vars
      use debug_tools

      real(kind=rkind), dimension(:), intent(out) :: resvct
      
      integer(kind=ikind) :: i, j, fin, posglob, posloc, subfin
      integer :: ierr
      logical :: skip_solved=.false.
     
      call mulAx_dd(pde_common%xvect(:,2), resvct)    
      

      !loop for b-(A.x) , where the vector (A.x) was calculated in the code above  
      do i=1, ubound(subdomain,1)
	if ( (.not. skip_solved) .or.  (.not. subdomain(i)%solved) ) then
	  fin = subdomain(i)%disjoint%pos
	  do j=1, fin
	    posglob = subdomain(i)%disjoint%data(j)
	    posloc = subdomain(i)%invpermut(posglob)
	    resvct(posglob) = subdomain(i)%bvect(posloc) - resvct(posglob) 
	  end do
	end if
      end do
      



    end subroutine get_residual
    

    subroutine get_residual_old(resvct)
      use typy
      use global_objs
      use pde_objs
      use globals
      use decomp_vars
      use debug_tools

      real(kind=rkind), dimension(:), intent(out) :: resvct
      
      integer(kind=ikind) :: i, j, fin, posglob, posloc, subfin
      integer :: ierr
      real(kind=rkind), dimension(:), allocatable :: vectorin
      type(extsmtx) :: mat
      logical :: skip_solved=.false.

      

      resvct = 0.0_rkind
      do i=1, ubound(subdomain,1)
        !loop for A.x
        if ( (.not. skip_solved) .or.  (.not. subdomain(i)%solved) ) then
          allocate(vectorin(subdomain(i)%ndof))
          vectorin = 0
          fin = subdomain(i)%disjoint%pos
          do j=1, fin
            posglob = subdomain(i)%disjoint%data(j)
            posloc = subdomain(i)%invpermut(posglob)
            vectorin(posloc) = pde_common%xvect(posglob,2)
          end do
          resvct = resvct +  subdomain(i)%matrix%mul(vectorin)
          deallocate(vectorin)
        end if
      end do
      
     
      

      !loop for b-(A.x) , where the vector (A.x) was calculated in the code above  
      do i=1, ubound(subdomain,1)
	if ( (.not. skip_solved) .or.  (.not. subdomain(i)%solved) ) then
	  fin = subdomain(i)%disjoint%pos
	  do j=1, fin
	    posglob = subdomain(i)%disjoint%data(j)
	    posloc = subdomain(i)%invpermut(posglob)
	    resvct(posglob) = subdomain(i)%bvect(posloc) - resvct(posglob) 
	  end do
	end if
      end do

!
    end subroutine get_residual_old
 
    
    subroutine collect_matrices(a)
      use typy
      use globals
      use global_objs
      use pde_objs
      use mtx
      use decomp_vars
      
      class(extsmtx), intent(out) :: a
      integer(kind=ikind) :: i, ii, jj, col
      real(kind=rkind) :: val
      
      i = ubound(pde,1)
      
      call a%init(maxval(pde(i)%permut(:)),maxval(pde(i)%permut(:)))
      
      do i=1, ubound(subdomain,1)
	do ii=1, subdomain(i)%matrix%getn()
	  do jj = 1, subdomain(i)%matrix%getm()
	    val = subdomain(i)%matrix%get(ii,jj)
	    if (abs(val) > epsilon(val) ) then
	      col = subdomain(i)%permut(jj)
	      call a%set(val, ii, col)
	    end if
	  end do
	end do
      end do
      
      if (.not. allocated(pde_common%bvect)) allocate(pde_common%bvect(maxval(pde(1)%permut(:))))
      
      do i=1, ubound(subdomain,1)
	pde_common%bvect(subdomain(i)%permut(1:subdomain(i)%ndof)) = subdomain(i)%bvect(1:subdomain(i)%ndof)
      end do
      
    
      	
    end subroutine collect_matrices
    
 
    subroutine addlevel()
      use typy
      use decomp_vars
      use debug_tools
      use globals
      use global_objs
      use pde_objs
      
      integer(kind=ikind) :: sub, coarse_el, j, i, vecino, el, elsub, k, level, l, nd, ii, vecsub, pnd, nd2
      integer(kind=ikind), dimension(:), allocatable, save :: weight
      !> an array that describes the element that was already inserted into subdomain
      logical, dimension(:), allocatable, save :: el_inserted
      logical, dimension(:), allocatable, save :: nd_inserted
      logical, dimension(:), allocatable, save :: nd_bellow
      logical :: insert
      integer(kind=ikind), dimension(:), allocatable :: indexes
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmprow, tmpvct
      
      if (.not. allocated(el_inserted)) then
	allocate(el_inserted(elements%kolik))
	allocate(nd_inserted(nodes%kolik))
	allocate(nd_bellow(nodes%kolik))
	allocate(weight(nodes%kolik))
      end if
      
      el_inserted = .false.
      nd_inserted = .false.
      nd_bellow = .false.
      
      weight = 0
      
      i = ubound(pde,1)
      call prolong_mtx%init(maxval(pde(i)%permut(:)), 1_ikind*ubound(subdomain,1))
            
      do sub=1, ubound(subdomain,1)
	do j=1, ubound(subdomain(sub)%coarse_el,1)
	  coarse_el = subdomain(sub)%coarse_el(j)
	  ddinfo%coarseinsub(coarse_el) = sub
	end do
      end do
      
     
     !searching for boundary clusters at subdomain boundary
      sub_loop: do sub = 1, ubound(subdomain,1)
	allocate(subdomain(sub)%olayers(0:ddinfo%overlaps))
	coarse_loop: do i=1,  ubound(subdomain(sub)%coarse_el,1)
	  coarse_el = subdomain(sub)%coarse_el(i)
	  coarse_vecino: do j=1, ubound(coarse_elements%neighbours,2)
	    vecino = coarse_elements%neighbours(coarse_el,j)
	    if (vecino > 0) then
	      if (ddinfo%coarseinsub(vecino) /= sub) then !nasel jsem krajni cluster
		call subdomain(sub)%bcluster%fill(coarse_el)
		EXIT coarse_vecino
	      end if
	    end if
	  end do coarse_vecino
	end do coarse_loop
      end do sub_loop
	
	      		
      sub_loop2: do sub=1, ubound(subdomain,1)
	j = subdomain(sub)%elsub%pos
	el_inserted(subdomain(sub)%elsub%data(1:j)) = .true.
	level_loop: do level=0, ddinfo%overlaps
	  do i=1, subdomain(sub)%bcluster%pos
	    coarse_el = subdomain(sub)%bcluster%data(i)
	    bound_loop: do k=1, ddcoarse_mesh(coarse_el)%overlap_el(level)%pos
	      el = ddcoarse_mesh(coarse_el)%overlap_el(level)%data(k)
	      elsub = ddinfo%coarseinsub(ddinfo%elincoarse(el)) 
	      if (level > 0) then
		if (elsub /= sub) then
		  insert = .true.
		else
		  insert = .false.
		end if
	      else
		insert = .false.
		do ii = 1, ubound(elements%neighbours,2)
		  vecino = elements%neighbours(el,ii)
		  if (vecino > 0) then
		    vecsub = ddinfo%coarseinsub(ddinfo%elincoarse(vecino))
		  else
		    vecsub = 0
		  end if
		  if (vecsub /= sub .and. vecsub /= 0) then
		    insert = .true.
		    EXIT
		  end if
	        end do
	      end if
	      if (level > 0 .and. insert .and. .not. el_inserted(el)) then
		call subdomain(sub)%elsub%fill(el)
		call subdomain(sub)%olayers(level)%fill(el)
		el_inserted(el) = .true.
	      else if (level == 0 .and. insert) then
		call subdomain(sub)%olayers(level)%fill(el)
	      end if
	    end do bound_loop
	  end do
	end do level_loop
	
	!-------
	! i f   s u b c y c l i n g 
	!-------
	!add an extra level including all neighbour nodes, here we start with adding an extra layer with elements, the neighbourhood elements are all elements that share at least a single node (different approach compared to the overlap construction, where the neighborhood elements are elements that share TWO nodes (edge))
	if (drutes_config%it_method == 2) then
	  do level=0, ddinfo%overlaps
	    do i=1, subdomain(sub)%olayers(level)%pos
	      el = subdomain(sub)%olayers(level)%data(i)
	      do j = 1, ubound(elements%data,2)
		nd = elements%data(el,j)
		do k=1, nodes%el2integ(nd)%pos
		  vecino = nodes%el2integ(nd)%data(k)
		  if (vecino > 0) then
		    if (.not. el_inserted(vecino)) then
		      call subdomain(sub)%elextra%fill(vecino)
		      el_inserted(vecino) = .true.
		    end if
		  end if
		end do
	      end do
	    end do
	  end do
	  j = subdomain(sub)%elextra%pos
	  el_inserted(subdomain(sub)%elextra%data(1:j)) = .false.
	end if
	
	
        j = subdomain(sub)%elsub%pos
	el_inserted(subdomain(sub)%elsub%data(1:j)) = .false.

      end do sub_loop2

	
            
      do sub=1, ubound(subdomain,1)
	allocate(subdomain(sub)%ondlayers(0:ddinfo%overlaps))
	do level=0, ddinfo%overlaps-1
	  do i=1, subdomain(sub)%olayers(level)%pos
	    el = subdomain(sub)%olayers(level)%data(i)
	    do j=1, ubound(elements%data,2)
	      nd = elements%data(el, j)
	      nd_bellow(nd) = .true.
	    end do
	  end do
	  do i=1, subdomain(sub)%olayers(level+1)%pos
	    el = subdomain(sub)%olayers(level+1)%data(i)
	    do j=1, ubound(elements%data,2)
	      nd = elements%data(el, j)
	      if (.not. nd_inserted(nd) .and. nd_bellow(nd) ) then
		call subdomain(sub)%ondlayers(level)%fill(nd)
		nd_inserted(nd) = .true.
              end if
            end do
          end do
        end do
        
        ! fill the last level
        level = ddinfo%overlaps
        do i=1, subdomain(sub)%olayers(level)%pos
	  el = subdomain(sub)%olayers(level)%data(i)
	  do j=1, ubound(elements%data,2)
	    nd = elements%data(el, j)
	    if (.not. nd_inserted(nd)) then
	      call subdomain(sub)%ondlayers(level)%fill(nd)
	      nd_inserted(nd) = .true.
	    end if
	  end do
	end do

	
       ! if subcycling, add nodes, that lie out of the domain, but still have a graph connection to the nodes, that lie inside the domain
	if (drutes_config%it_method == 2) then
	  do i = 1, subdomain(sub)%elextra%pos
	    el = subdomain(sub)%elextra%data(i)
	    do j=1, ubound(elements%data,2)
	      nd = elements%data(el, j)
	      call ddinfo%nodesinextsub(nd)%nrfill(sub)
	      if (.not. nd_inserted(nd)) then
		call subdomain(sub)%ndextra%fill(nd)
		nd_inserted(nd) = .true.
	      end if
	    end do
	  end do
	 
          j = subdomain(sub)%elextra%pos
          el_inserted(subdomain(sub)%elextra%data(1:j)) = .true.
	 
	  j = subdomain(sub)%ndextra%pos
	  do i = 1, j
	    nd = subdomain(sub)%ndextra%data(i)
	    do k = 1, nodes%el2integ(nd)%pos
	      el = nodes%el2integ(nd)%data(k)
	      if (.not. el_inserted(el)) then
		do l=1, ubound(elements%data,2)
		  nd2 = elements%data(el,l)
		  if (.not. nd_inserted(nd2)) then
		    call subdomain(sub)%ndextra%fill(nd2)
		    nd_inserted(nd2) = .true.
		  end if
		end do
	      end if
	    end do
	  end do
	  
	  j = subdomain(sub)%ndextra%pos
	  nd_inserted(subdomain(sub)%ndextra%data(1:j)) = .false.
	  
	  j = subdomain(sub)%elextra%pos
          el_inserted(subdomain(sub)%elextra%data(1:j)) = .false.
	end if
	  


	! set everywhere in nd_inserted .false.
	do level = 0, ddinfo%overlaps
	  j = subdomain(sub)%ondlayers(level)%pos
          if (allocated(subdomain(sub)%ondlayers(level)%data)) then
	    nd_inserted(subdomain(sub)%ondlayers(level)%data(1:j)) = .false.
	    nd_bellow(subdomain(sub)%ondlayers(level)%data(1:j)) = .false.
	  end if
	end do
	

	!create the prolongator matrix, set everywhere just ones 
	do i=1, subdomain(sub)%elsub%pos
	  el = subdomain(sub)%elsub%data(i)
	  do j=1, ubound(elements%data,2)
	    nd = elements%data(el, j)
	    pnd = pde(1)%permut(nd)
	    if (pnd > 0 .and. .not. nd_inserted(nd)) then
	      call prolong_mtx%set(1.0_rkind, pnd, sub)
	      nd_inserted(nd) = .true.
	    end if
	  end do
	end do
! 	if (.not. subdomain(sub)%critical) then
! 	  nd_inserted = .false.
! 	else
	  do i=1, subdomain(sub)%elsub%pos
	    el = subdomain(sub)%elsub%data(i)
            nd_inserted(elements%data(el,:)) = .false.
          end do
!         end if
      end do
      

      
      !adjust the prolongator matrix
      
      do sub=1, ubound(subdomain,1)
	do level = 1, ddinfo%overlaps
	  do j=1, subdomain(sub)%ondlayers(level)%pos
	    nd = subdomain(sub)%ondlayers(level)%data(j)
	    pnd = pde(1)%permut(nd)
	    if (pnd > 0) then
	      call prolong_mtx%set((ddinfo%overlaps-level)*1.0_rkind/ddinfo%overlaps, pnd, sub)
	    end if
	  end do
	end do
      end do
      
      do i=1, prolong_mtx%getn()
	call prolong_mtx%getrow(i=i, v=tmprow, jj=indexes, nelem=ii)
	if (ii > 1) then
	  tmp = 0.0_rkind
	  do j=1,ii
	    tmp = tmp + tmprow(j)
	  end do
	  tmprow = tmprow / tmp
	  do j=1, ii
	    call prolong_mtx%set(tmprow(j), i, indexes(j))
	  end do
	end if
      end do
      

      


	  
    
    end subroutine addlevel
    
    subroutine build_coarsemtx(resvct)
      use typy
      use globals
      use global_objs
      use pde_objs
      use decomp_vars
      use debug_tools
      
      
      real(kind=rkind), dimension(:), allocatable :: resvct
      real(kind=rkind), dimension(:,:), allocatable ::  Axcol
      real(kind=rkind), dimension(:), allocatable :: column
      integer(kind=ikind) :: fin, i, sub, j
      integer(kind=ikind), dimension(:), allocatable :: indexes
      real(kind=rkind), dimension(:), allocatable :: tmpdata
      real(kind=rkind), dimension(:), allocatable :: coarsecol
      real(kind=rkind) :: tmp
      
      fin = maxval(pde(1)%permut(:))
      
      allocate(Axcol(fin, ubound(subdomain,1)))
      allocate(column(fin))
                  
      do sub=1, ubound(subdomain,1)
	call prolong_mtx%getcol(sub, tmpdata, indexes, i)
	column = 0
	column(indexes(1:i)) = tmpdata(1:i)
	call mulAx_dd(column, Axcol(:,sub))
      end do
      
      i = ubound(subdomain,1)
      call coarse_mtx%init(i,i)
      
      allocate(coarsecol(ubound(subdomain,1)))
      
      do sub=1, ubound(subdomain,1)
	coarsecol = prolong_mtx%mulT(Axcol(:,sub))
	do i=1, ubound(coarsecol,1)
	  if (abs(coarsecol(i)) > epsilon(coarsecol(i))) then
	    call coarse_mtx%set(coarsecol(i), i, sub)
	  end if
	end do
      end do
      
      if (.not. allocated(bcoarse)) then
	allocate(bcoarse(ubound(subdomain,1)))
	allocate(xcoarse(ubound(subdomain,1)))
      else if (ubound(bcoarse,1) /= ubound(subdomain,1)) then
	deallocate(bcoarse)
	allocate(bcoarse(ubound(subdomain,1)))
	deallocate(xcoarse)
	allocate(xcoarse(ubound(subdomain,1)))
      end if
      
      bcoarse = prolong_mtx%mult(resvct)

    
      
    
    end subroutine build_coarsemtx
 
 
    subroutine init_addlevel()
      use typy
      use decomp_vars
      use global_objs
      use globals
      use debug_tools
      
      integer(kind=ikind) :: i, coarse_el, j, el, vecino, level
      logical :: border
      type(smartarray_int) :: tmp_arrays
      
      allocate(ddinfo%elincoarse(elements%kolik))      
	
      do coarse_el=1, ubound(ddcoarse_mesh,1)
	do i=1, ubound(ddcoarse_mesh(coarse_el)%elements,1)
	  el = ddcoarse_mesh(coarse_el)%elements(i)
	  ddinfo%elincoarse(el) = coarse_el
	end do
	allocate(ddcoarse_mesh(coarse_el)%overlap_el(0:ddinfo%overlaps))
      end do
	

      do coarse_el=1, ubound(ddcoarse_mesh,1)
	do i=1, ubound(ddcoarse_mesh(coarse_el)%elements,1)
	  el = ddcoarse_mesh(coarse_el)%elements(i)
	  border = .false.
	  vecino_loop: do j=1, ubound(elements%data,2)
			vecino = elements%neighbours(el, j)
			if (vecino > 0) then
			  if (ddinfo%elincoarse(el) /= ddinfo%elincoarse(vecino)) then
			    border = .true.
			    EXIT vecino_loop
			  end if
			end if
	  end do vecino_loop
	  if (border) then
	    call ddcoarse_mesh(coarse_el)%overlap_el(0)%fill(el)
	  end if
	end do
      end do
      
      
      do level = 1, ddinfo%overlaps
	do coarse_el=1, coarse_elements%kolik
	  border_loop: do i=1, ddcoarse_mesh(coarse_el)%overlap_el(level-1)%pos
			  el = ddcoarse_mesh(coarse_el)%overlap_el(level-1)%data(i)
			  do j=1, ubound(elements%data,2)	
			    vecino = elements%neighbours(el, j)
			    if (vecino > 0) then
			      if (coarse_el /= ddinfo%elincoarse(vecino)) then
				if (add_el(coarse_el, level, vecino)) then
				  call ddcoarse_mesh(coarse_el)%overlap_el(level)%nrfill(vecino)
				end if
			      end if
			    end if
			  end do 
          end do border_loop
        end do
      end do
      
 
 
      allocate(ddinfo%coarseinsub(coarse_elements%kolik))

    end subroutine init_addlevel
    
    function add_el(coarse_el, level, el) result(add)
      use typy
      use globals
      use global_objs
      use decomp_vars
      use debug_tools
      
      integer(kind=ikind), intent(in) :: coarse_el, level, el
      logical :: add
      
      integer(kind=ikind) :: i,j
      
      add = .true.
      
      do i=0, level-1
! 	call printmtx(ddcoarse_mesh(coarse_el)%overlap_el(i)) ; call wait()
	do j=1, ddcoarse_mesh(coarse_el)%overlap_el(i)%pos
! 	  print *, level-1, el ; stop
	  if (allocated(ddcoarse_mesh(coarse_el)%overlap_el(i)%data)) then
	    if (ddcoarse_mesh(coarse_el)%overlap_el(i)%data(j) == el) then
	      add = .false.
	      RETURN
	    end if
	  end if
	end do
      end do
    
    end function add_el



     subroutine print_domains(behavior)
      use typy 
      use decomp_vars
      use globals
      use global_objs
      use core_tools
      use geom_tools
      use pde_objs
      use debug_tools

      !> behavior specifies the behavior of this printing tool, it is a character
      !! behavior == "all_in_one" everythink is written into a single scilab exec file, scilab generates png files from all print levels
      !! behavior == "separately", each plot is written into different scilab exec file
      !<
      character(len=*), optional , intent(in)               :: behavior
      integer(kind=ikind)                                   :: mode
      character(len=256), dimension(-1:6)                   :: filenames
      integer(kind=ikind)                                   :: i, i_err, j, layer, dec, row, k, l, n, npoints, sub
      character(len=64)                                     :: forma
      integer, dimension(-1:6), save                        :: ids
      character(len=4)                                      :: extension
      real(kind=rkind)                                      :: distance, avgval, val1, val2, val3, tmp
      real(kind=rkind), dimension(:), allocatable, save     :: deriv
      real(kind=rkind), dimension(4,8) 			    :: body
      real(kind=rkind), dimension(2) 			    :: vct1, vct2
      real(kind=rkind), dimension(8) 			    :: vct_tmp
      type(smartarray_int)                                  :: nodes_domain
      type(smartarray_int)                                  :: nodes_domain_invprmt
      integer(kind=ikind), save                             :: decomp_print
      integer :: ierr


      if (present(behavior)) then
	select case(behavior)
	  case("all_in_one")
		mode = 0
	  case("separately")
		mode = 1
	  case("domains_only")
		mode = 2
	  case default
		print *, "incorrect argument of function decomp_tools::print_domains"
		ERROR STOP
	end select
      else
	mode = 0
      end if
      
      if (.not. allocated(deriv)) then
	allocate(deriv(drutes_config%dimen))
	decomp_print = 0
      end if

      decomp_print = decomp_print + 1

      dec = 1

      write(unit=forma, fmt="(a)") "(a,  I6.6 , a)"

      write(unit=filenames(0), fmt=forma) "out/dd/subiters-",  decomp_print, ".sci"

      write(unit=filenames(1), fmt=forma) "out/dd/subdomain-",  decomp_print, ".sci"
      
      write(unit=filenames(2), fmt=forma) "out/dd/petr/bside-", decomp_print, ".strana"

      write(unit=filenames(3), fmt=forma) "out/dd/domeny-", decomp_print, ".domeny"

      write(unit=filenames(4), fmt=forma) "out/dd/petr/matice-", decomp_print, ".matice"

      write(unit=filenames(5), fmt=forma) "out/dd/petr/prcd_matice-", decomp_print, ".matice"
      
      write(unit=filenames(6), fmt=forma) "out/dd/nodes-doms-", decomp_print, ".dat"
      
      write(unit=filenames(-1), fmt=forma) "out/dd/timestep-", decomp_print, ".sci" 

      npoints = ubound(coarse_elements%data,2)
      
      if (mode == 1 .or. decomp_print == 1) then
	do i=lbound(ids,1), ubound(ids,1)
	  call find_unit(ids(i), 6000)
	  open(unit=ids(i), file=trim(filenames(i)), action="write", status="replace", iostat=ierr)
	  if (ierr /= 0) then
	    i_err=system("mkdir out/dd/petr") 
	    if (i_err /= 0) then
	      i_err=system("mkdir out/dd/")
	      i_err=system("mkdir out/dd/petr")
	      open(unit=ids(i), file=trim(filenames(i)), action="write", status="replace", iostat=ierr)
	      if (i_err /=0 ) then
		print *, "unable to open and create directory out/dd (is it UN*X type os?), called  from decomp_tools::print_domains"
		ERROR STOP
	      end if
	    end if
	    open(unit=ids(i), file=trim(filenames(i)), action="write", status="replace", iostat=ierr)
	  end if
	end do  
      end if

     do l=-1,1	  
	write(unit=ids(l), fmt=*) "nt =", coarse_elements%kolik, ";"
	write(unit=ids(l), fmt=*) "x=zeros(nt,", npoints,");"
	write(unit=ids(l), fmt=*) "y=zeros(nt,", npoints ,");"
	write(unit=ids(l), fmt=*) "z=zeros(nt,", npoints, ");"


	do i=1, coarse_elements%kolik
          n = npoints + 1
	  do j=1,ubound(coarse_elements%data,2)
	    body(j,1:2) = coarse_nodes%data(coarse_elements%data(i,n-j),:)
	    body(j,4) = ddcoarse_mesh(i)%subdomain 
            body(j,3) = ddcoarse_mesh(i)%iter_count
            sub = ddinfo%coarseinsub(i)
            body(j,5) = subdomain(sub)%tmpval
	  end do
	  
	

	  
	  vct1 = body(3, 1:2) - body(1, 1:2) 
	  vct2 = body(3, 1:2) - body(2, 1:2)
	  
	  if ( (vct1(1)*vct2(2)-vct1(2)*vct2(1)) > 0.0_rkind) then
	    vct_tmp = body(3,:)
	    body(3,:) = body(2,:)
	    body(2,:) = vct_tmp
	  else
	    CONTINUE
	  end if


	  write(unit=ids(l), fmt=*) "x(", i, ",1) =", body(1,1), ";"
	  write(unit=ids(l), fmt=*) "x(", i, ",2) =", body(2,1), ";"
	  write(unit=ids(l), fmt=*) "x(", i, ",3) =", body(3,1), ";"
          if (npoints == 4) write(unit=ids(l), fmt=*) "x(", i, ",4) =", body(4,1), ";"


	  write(unit=ids(l), fmt=*) "y(", i, ",1) =", body(1,2), ";"
	  write(unit=ids(l), fmt=*) "y(", i, ",2) =", body(2,2), ";"
	  write(unit=ids(l), fmt=*) "y(", i, ",3) =", body(3,2), ";"
	  if (npoints == 4) write(unit=ids(l), fmt=*) "y(", i, ",4) =", body(4,2), ";"

	  if (l>=0) then
	      write(unit=ids(l), fmt=*) "z(", i, ",1) =", body(1,3+l), ";"
	      write(unit=ids(l), fmt=*) "z(", i, ",2) =", body(2,3+l), ";"
	      write(unit=ids(l), fmt=*) "z(", i, ",3) =", body(3,3+l), ";"
	    if (npoints == 4) write(unit=ids(l), fmt=*) "z(", i, ",4) =", body(4,3+l), ";"
	  else
	    write(unit=ids(l), fmt=*) "z(", i, ",1) =", body(1,5), ";"
	    write(unit=ids(l), fmt=*) "z(", i, ",2) =", body(2,5), ";"
	    write(unit=ids(l), fmt=*) "z(", i, ",3) =", body(3,5), ";"
	    if (npoints == 4) write(unit=ids(l), fmt=*) "z(", i, ",4) =", body(4,5), ";"
	  end if

	end do
	
    
	if (l > -1) then
	  write(unit=ids(l), fmt=*) "f=gcf();"
	  write(unit=ids(l), fmt=*) "clf(f,'reset');"
    ! !       write(unit=ids(1), fmt=*) "f.color_map=jetcolormap(max(z));"
	  write(unit=ids(l), fmt=*) "n=ceil(max(z));"
	  write(unit=ids(l), fmt=*)  "r=linspace(0,1,n)';"
	  write(unit=ids(l), fmt=*)  "g=linspace(1,0,n)';"
	  write(unit=ids(l), fmt=*)  " b=ones(r);"

	  write(unit=ids(l), fmt=*) "for i=1,n;"
	  write(unit=ids(l), fmt=*) "    if (i>1) then;"
	  write(unit=ids(l), fmt=*) "        if (b(i-1) > 0) then;"
	  write(unit=ids(l), fmt=*) "            b(i) = -1;"
	  write(unit=ids(l), fmt=*) "        else;"
	  write(unit=ids(l), fmt=*) "            b(i)=1;"
	  write(unit=ids(l), fmt=*) "        end;"
	  write(unit=ids(l), fmt=*) "    else;"
	  write(unit=ids(l), fmt=*) "        b(i) = 0;"
	  write(unit=ids(l), fmt=*) "    end;"
	  write(unit=ids(l), fmt=*) "end;"
	  write(unit=ids(l), fmt=*) "cmap=[r g b];"
	  write(unit=ids(l), fmt=*) "f=gcf();"
	  write(unit=ids(l), fmt=*) "f.color_map=cmap;"

	  write(unit=ids(l), fmt=*) "colorbar(floor(min(z)),ceil(max(z)),colminmax=[floor(min(z)) ceil(max(z))], fmt='%.0f');"
	  write(unit=ids(l), fmt=*) "plot3d1(x',y',z',alpha=0, theta=-90);"
	  ! syntax from scilab manual:  xs2png(0, 'foo.png')
	  write(unit=ids(l), fmt="(a,I6.6,a)" ) "xs2png(0, 'domains-", decomp_print , ".png');" 
	else
	  write(unit=ids(l), fmt=*) "f=gcf();"
	  write(unit=ids(l), fmt=*) "clf(f,'reset');"
	  write(unit=ids(l), fmt=*) "colorbar(min(z),max(z));"
	  write(unit=ids(l), fmt=*) "plot3d1(x',y',z',alpha=0, theta=-90);"
	  write(unit=ids(l), fmt="(a,I6.6,a)" ) "xs2png(0, 'domains_dt-", decomp_print , ".png');"
	end if

      end do

      select case(mode)
        case(0)
	  do l=lbound(ids,1), ubound(ids,1)
	    write(unit=ids(l), fmt=*) "clear;"
	    call flush(ids(l))
	  end do

        case(2)
	  do l=-1,1
	    close(ids(l))
	  end do

        case(1)

	  do i=1, ubound(pde_common%bvect,1)
	    write(unit=ids(2), fmt=*) pde_common%bvect(i)
	  end do

	  write(unit=ids(3), fmt=*) "#--node id (permutated in stiff. mtx)----node id (physical)---coordinates--------subdomain id----"
	  write(unit=ids(3), fmt=*) "#-----------------------------------------------------------------------------------------"

	  do i=1, coarse_elements%kolik
	    row = 0
	    do j=1,ubound(ddcoarse_mesh(i)%elements,1)
	      do k=1,3
		call nodes_domain%fill(pde(1)%permut(elements%data(ddcoarse_mesh(i)%elements(j),k)))
		call nodes_domain_invprmt%fill(elements%data(ddcoarse_mesh(i)%elements(j),k))
!                     nodes_domain(row) = pde(1)%permut(elements%data(ddcoarse_mesh(i)%elements(j),k), 2)
	      end do
	    end do  

	    do j=1, nodes_domain%pos
	      write(ids(3), fmt=*)  nodes_domain%data(j), nodes_domain_invprmt%data(j), &
	       nodes%data(nodes_domain_invprmt%data(j),:), ddcoarse_mesh(i)%subdomain
	    end do
	  end do

	  call nodes_domain%clear(full = .true.)
	  call nodes_domain_invprmt%clear(full = .true.)
	  
	  do i=1, nodes%kolik
	    write(unit = ids(6), fmt=*) i, pde(1)%solution(i)
	  end do
	      
	  do i=lbound(ids,1), ubound(ids,1)
	    close(ids(i))
	  end do
 
      end select


    end subroutine print_domains

     subroutine print_elements_dd(behavior)
      use typy 
      use decomp_vars
      use globals
      use global_objs
      use core_tools
      use geom_tools
      use debug_tools

      !> behavior specifies the behavior of this printing tool, it is a character
      !! behavior == "all_in_one" everythink is written into a single scilab exec file, scilab generates png files from all print levels
      !! behavior == "separately", each plot is written into different scilab exec file
      !<
      character(len=*), intent(in)                          :: behavior
      integer(kind=ikind)                                   :: mode
      character(len=256), dimension(0:1)                    :: filenames
      integer(kind=ikind)                                   :: i, i_err, j, layer, dec, row, k, l, n, npoints, el
      character(len=64)                                     :: forma
      integer, dimension(0:1), save                         :: ids
      character(len=4)                                      :: extension
      real(kind=rkind)                                      :: distance, avgval, val1, val2, val3, tmp
      real(kind=rkind), dimension(:), allocatable, save     :: deriv
      real(kind=rkind), dimension(4,8)                      :: body
      real(kind=rkind), dimension(2)                        :: vct1, vct2
      real(kind=rkind), dimension(8)                        :: vct_tmp
      integer(kind=ikind), dimension(:), allocatable        :: nodes_domain
      integer(kind=ikind), save                             :: decomp_print


      if (drutes_config%it_method /= 1) RETURN
      
      select case(behavior)
        case("all_in_one")
              mode = 0
              RETURN
        case("separately")
              mode = 1
        case("domains_only")
              mode = 2
        case default
              print *, "incorrect argument of function decomp_tools::print_domains"
              ERROR STOP
      end select
      
      if (.not. allocated(deriv)) then
        allocate(deriv(drutes_config%dimen))
        decomp_print = 0
      end if

      decomp_print = decomp_print + 1

      dec = 1

      allocate(nodes_domain(nodes%kolik))  

      write(unit=forma, fmt="(a)") "(a,  I6.6 , a)"

      write(unit=filenames(0), fmt=forma) "out/dd/subiters_el-",  decomp_print, ".sci"

      write(unit=filenames(1), fmt=forma) "out/dd/subdomain_el-",  decomp_print, ".sci"

      npoints = ubound(elements%data,2)
      
      if (mode == 1 .or. decomp_print == 1) then
        do i = lbound(ids,1), ubound(ids,1)
          call find_unit(ids(i), 6000)
          open(unit=ids(i), file=trim(filenames(i)), action="write", status="replace")
        end do
      end if

      do l=lbound(ids,1), ubound(ids,1)   
        write(unit=ids(l), fmt=*) "nt =", elements%kolik, ";"
        write(unit=ids(l), fmt=*) "x=zeros(nt,", npoints,");"
        write(unit=ids(l), fmt=*) "y=zeros(nt,", npoints ,");"
        write(unit=ids(l), fmt=*) "z=zeros(nt,", npoints, ");"
    

 
        do i=1, coarse_elements%kolik
          do k=1, ddcoarse_mesh(i)%el_count
            el = ddcoarse_mesh(i)%elements(k)
            n = npoints + 1
            do j=1,ubound(elements%data,2)
              body(j,1:2) = nodes%data(elements%data(el,n-j),:)
              body(j,4) = ddcoarse_mesh(i)%subdomain 
              body(j,3) = ddcoarse_mesh(i)%iter_count
            end do
            
            vct1 = body(3, 1:2) - body(1, 1:2) 
    	    vct2 = body(3, 1:2) - body(2, 1:2)
            
            if ( (vct1(1)*vct2(2)-vct1(2)*vct2(1)) > 0.0_rkind) then
	      vct_tmp = body(3,:)
    	      body(3,:) = body(2,:)
    	      body(2,:) = vct_tmp
    	    else
    	      CONTINUE
    	    end if


            write(unit=ids(l), fmt=*) "x(", el, ",1) =", body(1,1), ";"
            write(unit=ids(l), fmt=*) "x(", el, ",2) =", body(2,1), ";"
            write(unit=ids(l), fmt=*) "x(", el, ",3) =", body(3,1), ";"
            if (npoints == 4) write(unit=ids(l), fmt=*) "x(", el, ",4) =", body(4,1), ";"


            write(unit=ids(l), fmt=*) "y(", el, ",1) =", body(1,2), ";"
            write(unit=ids(l), fmt=*) "y(", el, ",2) =", body(2,2), ";"
            write(unit=ids(l), fmt=*) "y(", el, ",3) =", body(3,2), ";"
            if (npoints == 4) write(unit=ids(l), fmt=*) "y(", el, ",4) =", body(4,2), ";"

            write(unit=ids(l), fmt=*) "z(", el, ",1) =", body(1,3+l), ";"
            write(unit=ids(l), fmt=*) "z(", el, ",2) =", body(2,3+l), ";"
            write(unit=ids(l), fmt=*) "z(", el, ",3) =", body(3,3+l), ";"
            if (npoints == 4) write(unit=ids(l), fmt=*) "z(", el, ",4) =", body(4,3+l), ";"
        end do
      end do
      
    
      

      write(unit=ids(l), fmt=*) "f=gcf();"
      write(unit=ids(l), fmt=*) "clf(f,'reset');"
! !       write(unit=ids(1), fmt=*) "f.color_map=jetcolormap(max(z));"
      write(unit=ids(l), fmt=*) "n=max(z);"
      write(unit=ids(l), fmt=*)  "r=linspace(0,1,n)';"
      write(unit=ids(l), fmt=*)  "g=linspace(1,0,n)';"
      write(unit=ids(l), fmt=*)  " b=ones(r);"

      write(unit=ids(l), fmt=*) "for i=1,n;"
      write(unit=ids(l), fmt=*) "    if (i>1) then;"
      write(unit=ids(l), fmt=*) "        if (b(i-1) > 0) then;"
      write(unit=ids(l), fmt=*) "            b(i) = -1;"
      write(unit=ids(l), fmt=*) "        else;"
      write(unit=ids(l), fmt=*) "            b(i)=1;"
      write(unit=ids(l), fmt=*) "        end;"
      write(unit=ids(l), fmt=*) "    else;"
      write(unit=ids(l), fmt=*) "        b(i) = 0;"
      write(unit=ids(l), fmt=*) "    end;"
      write(unit=ids(l), fmt=*) "end;"
      write(unit=ids(l), fmt=*) "cmap=[r g b];"
      write(unit=ids(l), fmt=*) "f=gcf();"
      write(unit=ids(l), fmt=*) "f.color_map=cmap;"

      write(unit=ids(l), fmt=*) "colorbar(int(min(z)),int(max(z)),colminmax=[min(z) max(z)], fmt='%.0f');"
      write(unit=ids(l), fmt=*) "plot3d1(x',y',z',alpha=0, theta=-90);"
      ! syntax from scilab manual:  xs2png(0, 'foo.png')
      write(unit=ids(l), fmt="(a,I6.6,a)" ) "xs2png(0, 'domains-", decomp_print , ".png');" 
    end do

      select case(mode)
        case(0)
          do l=0,1
            write(unit=ids(l), fmt=*) "clear;"
            call flush(ids(l))
          end do

        case(2)
          do l=0,1
            close(ids(l))
          end do

        case(1)         
              do i=lbound(ids,1), ubound(ids,1)
                close(ids(i))
              end do
              

              deallocate(nodes_domain)
      end select



    end subroutine print_elements_dd
   

end module decomp_tools