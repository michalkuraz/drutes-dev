!> stiffness matrix procedure for non-linear problems linearized by Picard method
module stiffmat
  public :: build_stiff_np
  public :: build_bvect
  contains

     !> build local stifness matrix for nonlinear problems and Picard method
     subroutine build_stiff_np(el_id, dt, domain_id, quadpnt_in)
      use typy
      use globals
      use global_objs
      use pde_objs
      use feminittools
      use linAlg
      use debug_tools

      integer(kind=ikind), intent(in) :: el_id
      !> time step
      real(kind=rkind), intent(in) :: dt
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str), intent(in out), optional :: quadpnt_in

      
      integer(kind=ikind), dimension(:,:), allocatable, save :: layer
      !> sum of dispersive terms
      real(kind=rkind), dimension(1,1) :: dsum
      !> sum of convection terms
      real(kind=rkind) :: csum
      !> sum of reaction terms
      real(kind=rkind) :: rsum
      integer(kind=ikind) :: i,j,l, top, iproc, jproc, ii, jj, limits
      real(kind=rkind), dimension(1,3) :: u
      real(kind=rkind), dimension(3,1) :: v
      real(kind=rkind), dimension(1,3) :: w
      real(kind=rkind), dimension(3) :: conv
      real(kind=rkind), dimension(3,3) :: disp
      type(integpnt_str) :: quadpnt

      stiff_mat = 0

      if (.not. allocated(layer)) then
        allocate(layer(ubound(pde,1), ubound(pde,1)))
      end if

      do iproc=1,ubound(pde,1)
        do jproc=1,ubound(pde,1)
            layer(iproc, jproc) = elements%material(el_id, jproc)
        end do
      end do
      
      limits = ubound(stiff_mat,1)/ubound(pde,1)

      top = drutes_config%dimen
      
      if (present(quadpnt_in)) then
        quadpnt=quadpnt_in
      end if

      quadpnt%element = el_id
      quadpnt%column = 2
      quadpnt%type_pnt = "gqnd"
      
      if (present(domain_id)) then
	quadpnt%ddlocal = .true.
	quadpnt%subdom = domain_id
      end if
      

     
            
      do iproc=1,ubound(pde,1)
	do jproc=1, ubound(pde,1)
	  do i=1, ubound(stiff_mat,1)/ ubound(pde,1)
	    do j=1, ubound(stiff_mat,1)/ubound(pde,1)
	      dsum = 0
	      csum = 0
	      rsum = 0
	      
	      
	      v(1:top,1) = elements%ders(el_id,i,1:top)
	      u(1,1:top) = elements%ders(el_id,j,1:top)

	      do l=1, ubound(gauss_points%weight,1)
		quadpnt%order = l
	        call pde(iproc)%pde_fnc(jproc)%dispersion(pde(iproc), layer(iproc, jproc), &
		quadpnt, tensor=disp(1:top,1:top))
		w(:,1:top) =  matmul(u(:,1:top),disp(1:top,1:top))
		dsum = dsum - matmul(w(:,1:top) ,v(1:top,:))*gauss_points%weight(l)
	      end do

  
	      do l=1, ubound(gauss_points%weight,1)
		quadpnt%order = l
		call pde(iproc)%pde_fnc(jproc)%convection(pde(iproc), layer(iproc, jproc), quadpnt, &
		  vector_out=conv(1:top))
		  csum = csum - dot_product(u(1,1:top),conv(1:top))*base_fnc(i,l)*gauss_points%weight(l)
		  call pde(iproc)%pde_fnc(jproc)%der_convect(pde(iproc), layer(iproc, jproc), quadpnt, 	&
		  vector_out=conv(1:top))
		  w = base_fnc(i,l)*base_fnc(j,l)
		  csum = csum - dot_product(w(1,1:top), conv(1:top))*gauss_points%weight(l)
	      end do
	      

	      do l=1, ubound(gauss_points%weight,1)
	      	quadpnt%order = l
		rsum = rsum + pde(iproc)%pde_fnc(jproc)%reaction(pde(iproc),layer(iproc, jproc), &
		      quadpnt)*base_fnc(i,l)*base_fnc(j,l)*gauss_points%weight(l)
	      end do
	     
	      

	      ii = i + (iproc-1)*limits
	      jj = j + (jproc-1)*limits
	      
	      
	      
	      stiff_mat(ii,jj) = (dsum(1,1) + csum + rsum)

	    end do
	  end do 
	end do
      end do

      
     stiff_mat = stiff_mat/gauss_points%area*elements%areas(el_id)*dt
     
     
    end subroutine build_stiff_np
    
    
    subroutine build_bvect(el_id, dt, domain_id, quadpnt_in)    
      use typy
      use globals
      use global_objs
      use pde_objs
      use feminittools
      use linAlg
      use debug_tools
      
      integer(kind=ikind), intent(in) :: el_id
      real(kind=rkind), intent(in) :: dt
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str) , intent(in), optional :: quadpnt_in
      
      integer(kind=ikind) :: iproc, limits, ii, i, l
      real(kind=rkind) :: suma
      type(integpnt_str) :: quadpnt
      
      bside = 0
      
      if (present(quadpnt_in)) then
        quadpnt = quadpnt_in
      end if
      
      quadpnt%element = el_id
      quadpnt%column = 2
      quadpnt%type_pnt = "gqnd"
      
      limits = ubound(stiff_mat,1)/ubound(pde,1)

      if (present(domain_id)) then
	quadpnt%ddlocal = .true.
	quadpnt%subdom = domain_id
      end if
      
      do iproc = 1, ubound(pde,1)
	do i=1, ubound(stiff_mat,1)/ ubound(pde,1)
	  suma = 0
	  
	  do l=1, ubound(gauss_points%weight,1)
           quadpnt%order = l	  
	   suma = suma - pde(iproc)%pde_fnc(iproc)%zerord(pde(iproc), layer=elements%material(el_id, iproc), quadpnt=quadpnt) &
	   *gauss_points%weight(l)
	  end do
	  
	  ii = i + (iproc-1)*limits
	  
	  bside(ii) = suma*dt*elements%areas(el_id)
	  
	end do
      end do
    

    end subroutine build_bvect


end module stiffmat
