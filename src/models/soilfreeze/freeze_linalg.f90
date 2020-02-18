module freeze_linalg
public :: freeze_solver_face


contains

	subroutine freeze_solver_face(A,b,x,itmax1,reps1,ilev1,itfin1,repsfin1,&
				ll1,ll2,cond1,opcnt1,errcode1)
	  use mtx
	  use typy
	  use sparsematrix
	  use solvers
	  use global_objs
	  use pde_objs
	  use freeze_globs

	  implicit none
	  !> matice soustavy\n
	  !! musi poskytovat getn, getm, mul (nasobeni vektorem)
	  class(matrix), intent(in out) :: A
	  !> vektor prave strany
	  real(kind=rkind), dimension(:), intent(in) :: b
	  !> aproximace reseni, postupne menena
	  real(kind=rkind), dimension(:), intent(in out) :: x
	  !> maximalni povoleny pocet iteraci, default = n ( Rozmer matice)
	  integer(kind=ikind), intent(in), optional :: itmax1
	  !> pozadovana relativni zmena rezidua, default = 1e-6
	  real(kind=rkind), intent(in), optional :: reps1
	  !> informacni podrobnost\n
	  !> - 0 ... pracuj tise
	  !! - 1 ... minimalni informace
	  !! - 10 ... maximalni ukecanost
	  integer, intent(in), optional :: ilev1
	  !> skutecne provedeny pocet iteraci
	  integer(kind=ikind), intent(out), optional :: itfin1
	  !> skutecne dosazena relativni zmena residua
	  real(kind=rkind), intent(out), optional :: repsfin1
	  !> odhad nejvetsiho vlastniho cisla
	  real(kind=rkind), intent(out), optional :: ll1
	  !> odhad nejmensiho vlastniho cisla
	  real(kind=rkind), intent(out), optional :: ll2
	  !> odhad cisla podminenosti : cond1 = ll1/ll2
	  real(kind=rkind), intent(out), optional :: cond1
	  !> celkovy pocet provedenych operaci a cas behu
	  type(tcount), intent(out), optional :: opcnt1
	  !> kod pripadnr chyby
	  !! - 0 ... OK
	  !! - 1 ... matice neni ctvercova
	  !! - 2 ... nesouhlasi b
	  !! - 3 ... nesouhasi x
	  !! - 4 ... ani jeden z vektoru nesouhlasi
	  !! - 5 ... vycerpan povoleny pocet iteraci
	  !! - 6 ... prestalo klesat residuum i energie
	  integer, intent(out), optional :: errcode1
	  integer :: ilevel
	  
	  !local variables
	  integer(kind=ikind) :: size_water, size_heat, i, row_size, j
	  real(kind=rkind), dimension(:), allocatable :: values
	  integer(kind=ikind), dimension(:), allocatable :: col_indices
	  
	  type(extsmtx) :: mat4heat, mat4water
      integer(kind=ikind), dimension(:), allocatable, save :: pwat1, pwat2, pheat1, pheat2
	  
	  
	
	  if (.not. present(ilev1) ) then
		ilevel = 0
	  else
		ilevel = ilev1
	  end if
	  
	  ! first initialize small block matrix
	  
	  ! size of the small block matrix is
	  size_water = maxval(pde(wat)%permut(:))
	  size_heat = maxval(pde(heat_proc)%permut(:))
	  
	  call mat4heat%init(size_water, size_water)
	  
	  call mat4water%init(size_heat - size_water, size_heat - size_water)
	  
	  !set matrix for water
	  
	  do i=1, size_water
		call mat4water%getrow(i, values, col_indices, row_size)
		do j=1, row_size
		  call mat4water%set(values(j), i, col_indices(j))
		end do
	  end do

	
      if (.not. allocated(pwat1)) then
        allocate(pwat1(size_water))
        allocate(pwat2(size_water))
		allocate(pheat1(size_heat-size_water))
        allocate(pheat2(size_heat-size_water))
      end if
  
      call LDUd(mat4water, pivtype=0, ilev=0, perm1=pwat1, perm2=pwat2)
      
      call LDUback(mat4water, b(1:size_water), x(1:size_water), p1=pwat1, p2=pwat2)
       
	  

	end subroutine freeze_solver_face

end module freeze_linalg
