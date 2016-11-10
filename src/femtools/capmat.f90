module capmat


  !> implicit Euler for nonlinear problem and Picard iteration method, nondiagonal capacity matrix
  public :: impl_euler_np_nondiag
  !> implicit Euler for nonlinear problem and Picard iteration method, diagonal capacity matrix
  public :: impl_euler_np_diag
  !> steady state -> the capacity matrix is a zero matrix
  public :: steady_state_int
  !> routine that diagonalizes matrix - it lumps consistent capacity matrices
  private :: diagonalizer
  

 

  contains

    !> implicit euler method, nonlinear problem Picard iteration method - non-diagonal capacity matrix
    subroutine impl_euler_np_nondiag(el_id, domain_id, quadpnt_in)
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use feminittools
      use linAlg
      use debug_tools

      integer(kind=ikind), intent(in) :: el_id
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str), intent(in out), optional :: quadpnt_in
      
      integer(kind=ikind) :: i,j,l, m
      real(kind=rkind) :: tmp
      integer(kind=ikind) :: iproc, jproc, limits, ll, jj
      integer(kind=ikind), dimension(:,:), allocatable, save :: layer
      type(integpnt_str) :: quadpnt


      if (.not. allocated(layer)) then
        allocate(layer(ubound(pde,1), ubound(pde,1)))
      end if

      do iproc=1,ubound(pde,1)
        do jproc=1,ubound(pde,1)
            layer(iproc, jproc) = elements%material(el_id, jproc)
        end do
      end do


      cap_mat = 0
      
      if (present(quadpnt_in)) then
        quadpnt=quadpnt_in
      end if

      limits = ubound(cap_mat,1)/ubound(pde,1)
      quadpnt%column = 2
      quadpnt%element = el_id
      quadpnt%type_pnt = "gqnd"
      
      if (present(domain_id)) then
	quadpnt%ddlocal = .true.
	quadpnt%subdom = domain_id
      end if
      

      do iproc=1, ubound(pde,1)
        do jproc=1, ubound(pde,1)
          do l=1,  limits
            do j=1,  limits
              do i=1, ubound(gauss_points%weight,1)
		quadpnt%order = i
                tmp =  -pde(iproc)%pde_fnc(jproc)%elasticity(pde(iproc), layer(iproc, jproc), &
                quadpnt)*base_fnc(j,i)*base_fnc(l,i)*gauss_points%weight(i)
                ll = l + limits*(iproc-1)
                jj = j + limits*(jproc-1)
                cap_mat(ll,jj) = cap_mat(ll,jj) + tmp
              end do
            end do
          end do
        end do
      end do

      
      cap_mat = cap_mat/gauss_points%area*elements%areas(el_id)

      bside =  bside + matmul(cap_mat, elnode_prev)
      

    end subroutine impl_euler_np_nondiag



    !> implicit euler method, nonlinear problem Picard iteration method - diagonal capacity matrix
    subroutine impl_euler_np_diag(el_id, domain_id,  quadpnt_in)
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use feminittools
      use linAlg
      use debug_tools

      integer(kind=ikind), intent(in) :: el_id
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str), intent(in out), optional :: quadpnt_in
      
      
      integer(kind=ikind) :: i,j,l, m
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(2,2) :: locmat
      integer(kind=ikind) :: iproc, jproc, limits, ll, jj
      integer(kind=ikind), dimension(:,:), allocatable, save :: layer
      type(integpnt_str) :: quadpnt

      if (.not. allocated(layer)) then
        allocate(layer(ubound(pde,1), ubound(pde,1)))
      end if

      do iproc=1,ubound(pde,1)
        do jproc=1,ubound(pde,1)
            layer(iproc, jproc) = elements%material(el_id, jproc)
        end do
      end do

      cap_mat = 0
      locmat = 0
      limits = ubound(cap_mat,1)/ubound(pde,1)     
      quadpnt%column = 2
      quadpnt%element = el_id
      quadpnt%type_pnt = "gqnd"
      
      
      if (present(domain_id)) then
	quadpnt%ddlocal = .true.
	quadpnt%subdom = domain_id
      end if



      do iproc=1, ubound(pde,1)
        do jproc=1, ubound(pde,1)
          do l=1,  limits
            do j=1,  limits
              do i=1, ubound(gauss_points%weight,1)
		quadpnt%order = i
                tmp =  -pde(iproc)%pde_fnc(jproc)%elasticity(pde(iproc), layer(iproc, jproc), &
                quadpnt)*base_fnc(j,i)*base_fnc(l,i)*gauss_points%weight(i)
                ll = l + limits*(iproc-1)
                jj = j + limits*(jproc-1)
                cap_mat(ll,jj) = cap_mat(ll,jj) + tmp
              end do
            end do
          end do
        end do
      end do

      
     do iproc = 1, ubound(pde,1)
       do jproc = 1, ubound(pde,1) 
         ll = limits*(iproc-1)+1
         jj = limits*(jproc-1)+1
        call diagonalizer(cap_mat(ll:ll+limits-1, jj:jj+limits-1))
      end do
    end do


    
     cap_mat = cap_mat/gauss_points%area*elements%areas(el_id)
     

     bside = bside +   matmul(cap_mat, elnode_prev)
     
     
    end subroutine impl_euler_np_diag

    subroutine steady_state_int(el_id, domain_id, quadpnt_in)
      use typy
      use globals
      use global_objs
      use geom_tools

      integer(kind=ikind), intent(in) :: el_id
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str), intent(in out), optional :: quadpnt_in

      cap_mat = 0

    end subroutine steady_state_int

    !> routine that diagonalizes matrix - it lumps consistent capacity matrices
    subroutine diagonalizer(a)
      use typy
      
      real(kind=rkind), dimension(:,:), intent(in out) :: a

      integer(kind=ikind) :: i,j
      real(kind=rkind) :: tmp

      do i=1, ubound(a,1)
	tmp = 0.0_rkind
	do j=1, ubound(a,1)
	  tmp = tmp + a(i,j)
	end do
	
	a(i,i) = tmp

	if (i > 1) then
	  a(i,1:i-1) = 0
	end if

	if (i < ubound(a,1)) then
	  a(i, i+1:ubound(a,1)) = 0
	end if
      end do
    
    end subroutine diagonalizer



end module capmat