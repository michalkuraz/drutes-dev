!> restarted GMRES solver for linear algebra problem
module gmres_solver
!  use matrix_oper_internals
!  use matrix_oper
  implicit none
  public:: gmres
  public:: house, hdump, hsolve
  
  
  contains
  !c-------------------------------------------------------c
  !c     List of Routines:
  !c
  !c     gmres
  !c
  !c-------------------------------------------------------c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c     Non-symmetric solvers                             c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !> main GMRES subroutine
  subroutine gmres(nsize, x, d, nit,tol, prod, precond, nvec,  &
       tol2,iout,it, sum, not_converge)
    use typy

    !    implicit double precision (a-h,o-z)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c     preconditioned gmres.                               c
    !c+++++++++++++++++++++++++++++++++++++++++++++++++++++++c
    !c     Ref. Saad, Y. & Schultz M.H., 1986.   GMRES :         c
    !c     A Generalised                                         c
    !c     Minimal Residual Algorithm for Solving Non-symmetric  c
    !c     Linear Systems. SIAM J. Sci. Stat. Comput. V7 #3 July c
    !c     1986 pp 856-869.                                      c
    !c+++++++++++++++++++++++++++++++++++++++++++++++++++++++c
    !c     This version uses Fast Householder instead of Givensc
    !c     rotations to perform the QR factorisation in the      c
    !c     residual minimisation.                                c
    !c-------------------------------------------------------    c
    !c     solution of the m-vector linear system Ax = d or      c
    !c     x is feasible solution on input/      c
    !c     solution vector on output           c
    !c     d is the right hand vector for the    c
    !c     system Ax = d                       c
    !c     A is a non-singular non-symmetric matrixc
    !c     is implicitly defined in the segments prod(h,z,m) c
    !c     which compute the product of A with an m-vector       c
    !c     z to produce an m-vector h.                           c
    !c     this allows the user                                  c
    !c     to define their own (sparse) data structure.          c
    !c     the pre-conditioning matrix is defined by :           c
    !c     precond(z,r,m) must be                                c
    !c     written to solve the system  LRz = r.,                c
    !c     The preconditioned system is then                     c
    !c     (LR)**-1 A x = (LR)**-1 d                             c
    !c     where LR is an approximation to A                     c
    !c     m is the number of independent                        c
    !c     variables.                                            c
    !c     nit is the maximum number of                          c
    !c     iterations allowed.                                   c
    !c     tol is the stopping tolerance                         c
    !c     ibnd is an integer array. normally set ibnd(i)=1      c
    !c     for i=1,2,3,..,m. this is useful for finite           c
    !c     element problems where A is singular, but some        c
    !c     of the variable x(i) are known. usually in such       c
    !c     problems one would "overwrite the boundary            c
    !c     conditions" and overwrite A. This overwriting can     c
    !c     be completely avoided if one follows the following    c
    !c     rules :                                               c
    !c     (1) enter the routine with each of the known          c
    !c     values x(i) in their correct positions in the         c
    !c     solution vector x.                                    c
    !c     (2) the integer array ibnd should have the values    c
    !c     ibnd(i)=1 if x(i) is one of the unknown values  c
    !c     ibnd(i)=0 if x(i) is a known value.             c
    !c     the routine will then correctly find the solution     c
    !c     using the matrix A without overwriting.               c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c     Author : N.R.C Birkett                                c
    !c     Oxford University Computing Laboratory       c
    !c     11 Keble Road                                c
    !c     Oxford OX1 3QD                               c
    !c     email  : nrcb@uk.ac.oxford.na                         c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      include 'params.for'
    integer(kind=ikind) :: nsize, m, nit, nvec, iout
    integer(kind=ikind), intent(out) :: it
    real(kind=rkind), dimension(1:nsize), intent(inout) :: x
    real(kind=rkind), dimension(1:nsize), intent(in)    :: d
    real(kind=rkind), intent(out) :: sum
    real(kind=rkind), intent(in) :: tol, tol2
    !real, dimension(:), allocatable :: xold
    integer(kind=ikind), intent(inout) :: not_converge   ! not converge = 1, converge = 0

    integer(kind=ikind), parameter  :: mdim=61
    real(kind=rkind) :: normB, normX, normA, normR

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c     Version December 1992 : Economised storage            c
    !c     Continuous Householder updates => residual calculated c
    !c     from H matrix without need to calculate  x at each    c
    !c     sweep.                                                c
    !c     Max vectors allowed for is GMRES(mdim-1)              c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccc
    !ccc   Big vectors
    !ccc
    real(kind=rkind), dimension(:, :), allocatable ::  qd
    real(kind=rkind), dimension(:), allocatable ::  vk
    !ccc
    !ccc   Small vectors
    !ccc
    real(kind=rkind) :: amat(mdim,mdim),bvec(mdim),vv(mdim)
    real(kind=rkind) :: hmat(mdim,mdim),hvec(mdim)
    interface
       subroutine prod(b,x,nsize)
         use typy
         integer(kind=ikind), intent (in):: nsize
         real(kind=rkind), dimension(1:nsize), intent(in) :: x
         real(kind=rkind), dimension(1:nsize), intent(inout) :: b
       end subroutine prod

       subroutine precond(b,x,nsize)
         use typy
         integer(kind=ikind), intent (in):: nsize
         real(kind=rkind), dimension(1:nsize), intent(in) :: x
         real(kind=rkind), dimension(1:nsize), intent(inout) :: b
       end subroutine precond

    end interface

    !local variables:
    integer(kind=ikind) :: idump, ioption, ihist, iflag, i, j, kvec
    real(kind=rkind) :: dnorm, qtq, ss, sum0,rr, bi, hi
    integer :: file_history, file_iters

    m = nsize

    allocate( qd(1_ikind:nsize,1_ikind:mdim), vk(nsize), source = 0.0_rkind )
    !print*,'GMRES - allocate'


    !! TOL
    !allocate(xold(1:nsize) )
    !ccc
    !ccc nvec < mdim
    !ccc
    !cccc
    !cccc  ioption = 1 gives convergence if r'r < tol**2 *d'd
    !cccc  AND r'r < tol2**2
    !cccc  ioption = 2 gives convergence if r'r < tol**2 *r0'r0
    !cccc                              AND r'r < tol2**2
    !cccc
    !cccc  iout    = 1 writes residual at each sweep to standard output
    !cccc  iout    = 0 suppresses printing.
    !cccc
    !ccccc idump=1 to dump Upper Hessenberg matrix
    !ccccc in MATLAB format

    !j = 10
    !do i=1, 5 !nsize/j
    !   write(*,'(i5,20es14.6)') i*j,x((i-1)*j+1 : i*j)
    !enddo
    !write(*,*) '----------------'
    !do i=1, 5 !nsize/j
    !   write(*,'(i5,20es14.6)') i*j,d((i-1)*j+1 : i*j)
    !enddo
    !write(*,*) '----------------'

    !
    !print*,'^^^^^^^^^^^^^^^^^^^^^^'
    !iout = 1
    !print*,'^^^^^^^^^^^^^^^^^^^^^^'!,sparse(506)
    
    if(m.gt.nsize) stop 'Array error in gmres. Too many variables'
    idump=0
    ioption=2

    ihist=0
    if(ihist.eq.1) then
       open(newunit=file_history,file="out/history.cg",status='UNKNOWN')
    endif



    !cccc  h(1,1) contains the norm of the residual at the
    !cccc  start of nvec sweeps.
    !cccc  abs(h(i+1,i+1)) contains the norm of the residual
    !cccc  after i sweeps.
    !cccc


    !init
    not_converge = 0

    if(iout.eq.1) then
    !   write(*,400)nvec
    endif
400 format('GMRES3(',I4,')')
    if(nvec.gt.mdim-1) then
       print *,'Array bound error in gmres3'
       print *,'re-set mdim >',nvec
       print *,'Solving with ',mdim-1,' vectors'
       nvec=mdim-1
    endif
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c     This algorithm is GMRES(nvec). Increasing nvec(< mdim) c
    !c     will in general reduce the number of iterations to     c
    !c     convergence, but increase the execution time per step  c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !ccc
    !ccc   Initialisation
    !ccc   it

    iflag=0
    it=0

    !ccc
    !ccc   Loop back to  78  instead of 500
    !ccc   if non-recursive residual required
    !ccc
78  continue


    dnorm=0.0
    !ccc
    !ccc   Get new rhs (LU)**-1 d
    !ccc

    normB = norm2(d(:) )
    !normA = MatrixNormFrobenius()
    normA = 1.

    vk(1:nsize) = d(1:nsize)

    call precond(qd(1,1),vk,nsize)
    !write(*,'(a10, 200es14.6)') 'SOL   b:',   vk(1:nsize)
    !write(*,'(a10, 200es14.6)') 'SOL ilu b:', qd(1:nsize,1)

    !write(*,'(a9,20es12.4)') ' Pb:',dot_product(qd(:, 1), qd(:, 1))**0.5

    !ccc
    !ccc   (LU)**-1 Ax (arbitrary preconditioner defined by zsol)
    !ccc
    !write(*,'(a10, 200es14.6)') 'SOL x:', x(1:nsize)
    call prod(vk,x,nsize)
    !write(*,'(a10, 200es14.6)') 'SOL Ax:', vk(1:nsize)
    
    call precond(qd(1,2),vk,nsize)
    !write(*,'(a10, 200es14.6)') 'SOL iluAx:', qd(1:nsize,2)

    !ccc
    !ccc   The initial (preconditioned) residual is now in qd(*,1)
    !ccc

    qd(1:nsize,1) = qd(1:nsize,1) - qd(1:nsize,2)
    dnorm = dot_product(qd(1:m,1), qd(1:m,1))

    !write(*,'(a9,20es12.4)') 'P(Ax-b):',dot_product(qd(:, 1), qd(:, 1))**0.5, dnorm**0.5

    !ccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccc
    !ccc   MAIN GMRES(NVEC) LOOP       ccc
    !ccc   Loop to 78 for non-recursiveccc
    !ccc   residual.                   ccc
    !ccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccc
500 continue
    !c     write to history file
    if(ihist.eq.1) then
       call prod(vk,x,nsize)
       rr=0.0D0
       do  j=1,m
          rr=rr+(vk(j)-d(j))**2 !*ibnd(j)
       enddo
       write(file_history,*)it,rr
       ! 2000    format(i6,3(1x,1PD21.14))
    endif

    bvec(1:nvec+1) = 0.
    sum = dot_product(qd(1:m,1), qd(1:m,1))

    qtq=sqrt(sum)
    amat(1,1)=qtq
    hmat(1,1)=qtq
    bvec(1)=  qtq
    if(qtq > 0.) then
       ss=1.0/qtq
       qd(1:m,1) = ss * qd(1:m,1)
    endif




1000 format('Pass...',i6,' residual**2 = ',1pd10.3)
    if(it.eq.0) then
       if(ioption.eq.1) sum0=dnorm
       if(ioption.eq.2) sum0=sum
       if(iout.eq.1) then
          !write(*,1000)it,sum
          !print*,'#####################################'
          print*
          write(*,'(i5,3(es9.2,a2), 5(es9.2,a4))')  &
               it,sum,' < ',tol**2*sum0,' ; ',sum,' < ', tol2**2, &
               'sum=',sum,'su0=',sum0,'tol=',tol,'to2=',tol
          open(newunit=file_iters, file = "out/GMRES_iter", status='UNKNOWN', position='append')
          write(file_iters,'(x)')
          write(file_iters,'(i5,6es12.4)')it,sum,sum0,tol,tol2,tol**2*sum0,tol2**2
          close(file_iters)
       endif
    endif


    !write(*,'(a9,20es12.4)') 'p(Ax-b):',dot_product(qd(:, 1), qd(:, 1))**0.5, &
    !     sum**0.5, sum0**0.5, sum**0.5/ sum0**0.5

    !write(*,'(a6,i8, 3(es12.4,a3,es12.4, a3))') 'gm???1', it, &
    !     sum,'<?', tol**2*sum0, '|', sum, '<?', tol2**2, '|',sum0,'<?',1E-30
    
    if( (sum.le.tol**2*sum0 .and. sum.le.tol2**2) .or. sum .le. 1.E-30 ) then
    


       if(ihist.eq.1) then
          call prod(vk,x,nsize)
          rr=0.0
          do  j=1,m
             rr=rr+(vk(j)-d(j))**2 !*ibnd(j)
          enddo
          write(file_history,*)it,rr
          close(file_history)
       endif

       ! evaluation of the error using backward analysis
       !normX = VectorNorm(x(:) )
       normX = norm2(x(:) )

       call prod(vk,x,nsize)
       !do j=1, nsize
       !   write(*,'(a8, 2i5, 800es12.4)') '>)(I(J', j, nsize, vk(j) - d(j)
       !enddo
       
       normR = norm2( vk(1:nsize) - d(1:nsize) )

!       if(normB > 0. .or. normX > 0.) &
!            state%linSolver%backward  = normR/(normB+normA * normX)
       !write(*,'(a4,6es12.4)') 'gmr:', normR, normB, normA, normX,state%linSolver%backward

       !print*,'GMRES - deallocate'
       deallocate( qd, vk )

       return   ! END oF GMRES, sucessfull step
    endif

    kvec=nvec
    do i=1,nvec
       !ccc
       !ccc   w(i):=Aq(i)
       !ccc
       call prod(vk,qd(1,i),nsize)

       !print*
       !write(*,'(a10, 200es14.6)') 'SOL!Ax:', vk(1:nsize)

       call precond(qd(1,i+1),vk,nsize)

       !write(*,'(a10, 200es14.6)') 'SOL!iluAx:', qd(1:nsize, i+1)

       !ccc
       !ccc Gram-Schmidt Orthonormalise to get {q(1),q(2),..,q(nvec+1)}
       !ccc (It would be better to use modified  Gramm-Schmidt
       !ccc  if nvec large)
       do j=1,i

          amat(j,i+1) = dot_product(qd(1:m,j),qd(1:m,i+1) )

          qd(1:m,i+1) = qd(1:m,i+1) - amat(j,i+1) * qd(1:m,j)

       enddo

       !ccc
       !ccc Normalise q(i+1)
       !ccc

       qtq = dot_product(qd(1:m,i+1), qd(1:m,i+1))**0.5

       ss=1./qtq
       amat(i+1,i+1)=qtq
       qd(1:m,i+1) = ss * qd(1:m,i+1)
       !ccc
       !ccc Copy next column of H as Householder overwrites.
       !ccc

       hmat(1:i+1,i+1) = amat(1:i+1,i+1)
       call house(amat(1,2),vv,bvec,i-1,mdim)
       it=it+1

       !ccc
       !ccc     The 2-norm of the residual is available
       !ccc as abs(bvec(i+1) at each sweep without
       !ccc     needing to calculate r:=r-Av
       !ccc
       !ccc

       sum=bvec(i+1)**2

       !         if(iout.eq.1 .and. mod(it,100) .eq. 0) then
       if(iout.eq.1 ) then
          !if(mod(it,10) == 1)
          !write(*,'(i5,6es12.4)')it,sum,sum0,tol,tol2,tol**2*sum0,tol2**2
          write(*,'(i5,3(es9.2,a2), 5(es9.2,a4))')  &
               it,sum,' < ',tol**2*sum0,' ; ',sum,' < ', tol2**2, &
               'sum=',sum,'su0=',sum0,'tol=',tol,'to2=',tol
          open(newunit=file_iters, file = "out/GMRES_iter", status='UNKNOWN', position='append')
          write(file_iters,'(i5,6es12.4)')it,sum,sum0,tol,tol2,tol**2*sum0,tol2**2
          close(file_iters)
       endif

       !stop

       !write(*,'(a9,20es12.4)') 'Q(Ax-b):',dot_product(qd(:, i+1), qd(:,i+1))**0.5,&
       !     sum**0.5, sum0**0.5

       !write(*,'(a6,6es12.4)') 'gm???2',sum, tol**2*sum0, sum,tol2**2, sum0
       ! EARLIER ESCAPE FROM ONE RESTART CYCLE
       if( (sum.le.tol**2*sum0 .and. sum.le.tol2**2) .or. sum .le. 1.E-28  ) then
          !print*,'gm: ESCAPE ??'
          iflag=1
          kvec=i
          goto 5000
       endif
    enddo

    !ccc
    !ccc   At the end of nvec cycles update r,x by completing
    !ccc   the QR calculation.
    !ccc   Solve the least squares problem
    !ccc   min || rk - AQr0||2
    !ccc   which is the same as
    !ccc   min || (amat(1,1),0,0,0..0)T- amat r0 ||
    !ccc   where amat is nvec+1 X nvec upper Hessenberg
    !ccc   Use fast version of Householder for this
    !ccc
5000 continue
    if(idump.eq.1) then
       call hdump(hmat,kvec,mdim)
       idump=0
    endif
    call hsolve(amat(1,2),bvec,kvec+1,mdim)

    !ccc
    !ccc   x(k+1)=x(k)+v(k)
    !ccc   r(k+1)=r(k)-w(k)
    !ccc   where v(k) = Q r0
    !ccc   w(k) = W r0
    !ccc

    ! 3333 format(A2,20(1x,f10.3))

    !ccc
    !ccc   compute hvec = H'bvec, then W bvec = QH'bvec = Q hvec
    !ccc
    !      call zero(hvec,kvec+1)
    hvec(1:kvec+1) = 0.
    do  j=2,kvec+1
       bi=-bvec(j-1)
       hvec(1:j) = hvec(1:j) - bi*hmat(1:j,j)
    enddo
    !ccc
    !ccc   Restore previous residual to q(*,1)
    !ccc   by mutiplying by mod(r) so do i=1
    !ccc   as special case.
    !ccc
    i=1
    bi=bvec(i)
    hi=hmat(1,1)-hvec(1)
    !call aminvb(x,qd(1,i),-bi,m)
    x(1:m) = x(1:m) + bi*qd(1:m, i)

    qd(1:m,1) = hi * qd(1:m,1)

    !ccc
    !ccc   Do the rest
    !ccc
    do i=2,kvec
       bi=bvec(i)
       hi=hvec(i)

       x(1:m) = x(1:m) + bi * qd(1:m, i)

       qd(1:m,1) = qd(1:m,1) - hi * qd(1:m,i)

    enddo

    if(hvec(kvec+1) /= 0.) &      ! correction of the GMRES code???
         qd(1:m,1) = qd(1:m,1) - hvec(kvec+1) * qd(1:m,kvec+1)


    ! a priori expecting of non convergence, causes a limitation of time step
    !if(it .ge. nit*0.75 .and. not_converge .ne. 1)  not_converge = 1

    if(it.ge.nit) goto 9000
    !ccc
    !ccc   Loop back to 500
    !ccc   to use recursive residual
    !ccc   otherwise loop back to 78
    !ccc   (exact residual requires an
    !ccc   extra mat/vec multiply.
    if (iflag.eq.0) then
       goto 78
    else
       goto 500
    endif
9000 continue
!    if(state%space%adapt%adapt_method /= 'ALG2') &
!         print *,'Gmres  failed to converge in ',it,' iterations', ',   rez =',sum
       ! evaluation of the error using backward analysis
    normX = norm2(x(:) )

    call prod(vk,x,nsize)
    normR = norm2(vk(:) -d(:) )

!    state%linSolver%backward  = normR/(normB+normA * normX)
    !write(*,'(a4,6es12.4)') 'GMR:', normR, normB, normA, normX,state%linSolver%backward


    not_converge = 1
    if(ihist.eq.1) then
       close(file_history)
    endif

    !print*,'GMRES - deallocate 2'
    deallocate( qd, vk )

    !deallocate(xold )

  end subroutine gmres

    !> upper Hessenberg matrix
    subroutine house(A,V,B,L,MDIM)
      use typy
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL(kind=rkind) :: A(MDIM,*),B(*),V(*),hk,vb,d, S, SS
      integer(kind=ikind) :: mdim, k,l
!c======================================================================c
!c     Version for m X m-1 upper Hessenberg matrix                          c
!c     For solution of the least-squares problem min || Ax-b||  A is M X N  c
!c     M> = N, first use the routine QR to factor A. Then do the            c
!c     2 steps b:= Q'b, solve for x the equations Rx = b, and overwrites thec
!c     solution on b. This routine may be called any number of times        c
!c     Re-ordered routine written for GMRES algorithm :                     c
!c     (1) Apply Householder transformations to next column generated       c
!c     from Arnoldi ,                                                   c
!c     A(*,L+1):=P(L)P(L-1)..P(1)A(*,L+1)                                   c
!c     Calculate new P(L+1) to eliminate A(L+1,L)                           c
!c     and apply P(L+1) to vector b                                         c
!c     V is a work vector.                                                  c
!c======================================================================c
      DO  K=1,L
         HK=-V(K)*0.5/A(K,K)

         VB = A(K+1,K) * a(k+1,L+1)

         vb=vb+hk*A(k,L+1)
         D=2.0*VB/V(K)

         a(k+1,L+1) = a(k+1,L+1) - d * A(k+1,k)

         a(k,L+1)=a(k,L+1)-d*hk

         !write(*,'(a6,i5,5es12.4)') '??*',k, hk,vb,d, S, SS
         !print*, '??*',k, hk,vb,d, S, SS
      enddo
      K=L+1
      SS = dot_product(A(k:k+1,k), A(k:k+1,k))

      S=SQRT(SS)
      IF (A(K,K).LT. 0.0) S=-S
!c     Hkk = Vk/(-2rkk)
!C     STORE VTV IN V
      V(K)=2.0*SS+2.0*A(K,K)*S
!C     APPLY TRANSFORMATION,GET NEW A
      A(K,K)=-S
!c     and apply to get new right handside
      HK=-V(K)*0.5/A(K,K)

      VB= A(K+1,K)* b(k+1)

      vb=vb+hk*b(k)
      D=2.0*VB/V(K)
      b(k+1) = b(k+1) - d * A(k+1,k)

      b(k)=b(k)-d*hk
    end subroutine house


    subroutine hdump(hmat,m,mdim)
      use typy
      integer(kind=ikind) :: m, mdim, i,j
      real(kind=rkind) :: hmat(mdim,*)
      !ccc
      !ccc   Dump the upper hessenberg matrix in
      !ccc   matlab format
      !ccc
      open(37,file='hmat',status='UNKNOWN')
      do  i=1,m
         write(37,1000)(hmat(i,j),j=2,m+1)
1000     format(1000(1PE21.14,1x))
      enddo
      close(37)
    end subroutine hdump

    !> solve Rx = Q'b and overwrite the solution on b
    subroutine hsolve(A,B,M,MDIM)
      use typy
      real(kind=rkind) :: A(MDIM,*),B(*)
      integer(kind=ikind) :: m, mdim, n,k
      N=M-1
      !c=================================================================c
      !c     Now solve Rx = Q'b and overwrite the solution on b
      !c=================================================================c
      do  k=n,1,-1
         b(k)=b(k)/a(k,k)
         b(1:k-1) = b(1:k-1)  - b(k)* a(1:k-1,k)
      enddo
    end subroutine hsolve

  end module gmres_solver



