!>\file LinAlg.f90
!!\brief  zakladni utilitky linearni algebry
!!
!<

!>
!! obsahuje:
!!   -  I/O operace pro linearni algebru
!!   - reseni soustavy linearnich rovnic
!! 
!<
module LinAlg
   public :: gem 
   public :: print_matrix
   public :: print_vector
   public :: LUdecomp
   public :: LUmul
   
   private :: fmax
contains 

  !>
  !! vytiskne matici
  !<
  subroutine Print_Matrix(A,nc,width)
   ! vytiskne matici, pocet sloupcu tisku je nc
   use typy
   !vytiskne matici
   ! simuluje format tisku matic Matlabu

      !>matice k tisknuti
      real(kind=rkind), dimension(:,:), intent(in) :: A
      !> pocet tistenych sloupcu,nepovinny, default = 6
      integer(kind=ikind), intent(in), optional    :: nc
      !> sirka cisel
      integer, intent(in), optional                :: width
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !lokalni promenne
      integer(kind=ikind) :: nc1  !< interni pocet sloupcu pro tisk
      integer(kind=ikind) :: ncol !< pocet sloupcu matice
      integer(kind=ikind) :: nrow !< pocet radku matice
      integer(kind=ikind) :: i    !< 
      integer(kind=ikind) :: j    !<
      integer(kind=ikind) :: col1 !<
      integer(kind=ikind) :: col2 !<
      integer                  :: w1  !<
      character(len=200):: fm   !< pouzity formatovaci retezec 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !vlastni kod
      !vyridime defaulty
      if ( present(nc) ) then
        nc1 = nc
      else
        nc1 = 6
      end if
      
      if ( present(width) ) then
         w1 = width
      else
         w1 = 20
      end if
      !ted opravdu tisknem
      ncol = ubound(A,2)
      nrow = ubound(A,1)
      do i=1,ncol/nc1
        col1 = (i-1)*nc1+1
        col2 = i*nc1
        write(unit=fm,fmt="(A1,I0,A2,I0,A1,I0,A1)")   "(",nc1,"Es",w1+8,".",w1,")"
        print *, "tisknu sloupce : (",col1,",",col2,")"
        do j=1,nrow
          print fm, A(j,col1:col2)
        end do
        print *, "  " 
      end do
      if (modulo(ncol,nc1) /= 0 ) then
        col1 = (ncol/nc1)*nc1+1
        col2 = ncol      
        print *, "tisknu sloupce : (",col1,",",col2,")"
        do j=1,nrow
          print "(100E30.20)", A(j,col1:col2)
        end do
        print *, "  " 
      end if
   end subroutine Print_Matrix

  !>
  !! vytiskne vektor
  !<
   subroutine Print_Vector(V,nc)
   ! vytiskne vektor, pocet sloupcu tisku je nc
     use typy
     !parametry
      !> tisteny vektor
      real(kind=rkind), dimension(:), intent(in) :: V
      !> pocet tistenych sloupcu, default = 6
      integer(kind=ikind), intent(in), optional  :: nc 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !lokalni promenne
      integer(kind=ikind) :: nc1, ncol, i, col1, col2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !vlastni kod
      !vyridime defaulty
      if ( present(nc) ) then
        nc1 = nc
      else
        nc1 = 6
      end if
      !ted opravdu tisknem
      ncol = ubound(V,1)
      do i=1,ncol/nc1
        col1 = (i-1)*nc1+1
        col2 = i*nc1
        print *, "tisknu usek : (",col1,",",col2,")"
        print "(100E30.20)", V(col1:col2)
        print *, "  " 
      end do
      if (modulo(ncol,nc1) /= 0) then
        col1 = (ncol/nc1)*nc1+1
        col2 = ncol      
        print *, "tisknu usek : (",col1,",",col2,")"
        print "(100E30.20)", V(col1:col2)
        print *, "  "    
      end if   
  end subroutine Print_Vector


  pure &
  
  !>
  !! \brief Gaussova eliminace
  !!
  !! Gaussova eliminace
  !!
  !<
  subroutine gem(a,b,x)
    use typy
    real(kind=rkind), dimension(:,:), intent(in) :: a !> matice soustavy
    !realny 2rozmer pole je to parametr a vstupni parametr
    real(kind=rkind), dimension(:), intent(in) :: b !> prava strana
    real(kind=rkind), dimension(:), intent(out) :: x !> reseni
    real(kind=rkind), dimension(:), allocatable :: aradek
    !matice a pochazi z main, nemuzu ji tedy menit,
    !proto si definuju matici aa, do ktery zapisu matici a
    real(kind=rkind), dimension(:,:), allocatable :: aa
    real(kind=rkind) :: wrk


    integer :: asizeV, asizeH,i,j
    integer :: bsize, n,k, maxindex
    integer :: mm
    integer :: xsize

    asizeV=ubound(a,1)
    !priradi rozmer horizontalni 
    asizeH=ubound(a,2)
    !priradi rozmer vertikalni 
    bsize=ubound(b,1)
    xsize=ubound(x,1)
    !print *, asizeH , asizeV , bsize , xsize
    if (asizeH==xsize .AND. asizeV==bsize) then
      !    print *, "zatim ok"
      if (asizeH/=asizeV) then
        !   print *, "matice neni |_|"
      else 
        !    print *, "jdeme pocitat"
        allocate(aa(asizeV,asizeH+1))
        allocate(aradek(1:asizeH+1))
        aa(1:asizeH,1:asizeV)=a
        aa(:,asizeV+1)=b
        !    x=aa(:,asizeV+1) -pokus pro vypis vektoru vysledu
        n=asizeH
        do i=1,n-1
          !prvni hodnata UDAVA RADKY, DRUHA SLOUPEC
          !MAXINDEX je poloha max. cisla  z radky i az n a sloupec i
          maxindex = i
          wrk = abs(aa(i,i))
          do j=i+1,n
            if (abs(aa(j,i)) > wrk ) then
              wrk = abs(aa(j,i))
              maxindex  = j
            end if
          end do
          !jestlize se maxindex rovna i nic neni treba delat,
          !je-li ruzny, tzn. > nez i je treba vymenit i radek a maxindex radek
          !klasicka zamena prvku
          if (maxindex > i) then
            !a radek priradim zkoumany radek (pouze uschova)
            aradek=aa(i,:)
            !priradim radek s maximalni hodnoutou do reseneho radku...
            aa(i,:)=aa(maxindex,:)   
            aa(maxindex,:)=aradek
            if (aa(i,i) == 0) then
              !print *,"matice je singularni..."
              return
            endif
          endif
          !vydelim reseny radek hodnotu nejvetsiho prvku 
          !(ktery sme predtim umistili na diagonalu)
          wrk=aa(i,i)
          aa(i,i:)=aa(i,i:)/wrk
          do j= i+1,n
            !vezmu j-ty radek i-ty sloupec
            wrk=aa(j,i)
            !aa(j,:) - vybiram cely j-ty radek
            aa(j,:)=aa(j,:)-wrk*aa(i,:)
          end do
          !call printmatrix(aa)
        enddo
        aa(n,n+1)=aa(n,n+1)/aa(n,n)
        do i=n-1, 1, -1
          !dotproduct - funkce provedejici skalarni soucin
          aa(i,n+1)=(aa(i,n+1)- dot_product( aa(i,i+1:n),aa(i+1:n,n+1)))/aa(i,i)
        end do
        !aa(:,n+1) je vektor reseni soustavy
        x=aa(:,n+1)
        deallocate(aa) 
        deallocate(aradek)
      endif
    else
      ! print *, "co delas??"
    endif
  end subroutine gem
 
 
 
  !>
  !!  \brief  LU rozklad
  !!
  !! LU rozklad s uplnym vyberem hlavniho prvku\n
  !! ! realne behem rozkladu presouva prvky
  !! tedy L a U jsou opravdu trojuhelniky 
  !<
  subroutine LUdecomp(A,p,q)
    use typy
    !> ctvercova matice - vstup: rozkladana matice\n
    !! ..................vystup: LU faktory
    !< 
    Real(kind=rkind),dimension(:,:), intent(in out) :: A
    !> radkova permutace, vystupni parametr
    Integer(kind=ikind), dimension(:), intent(out)  :: p
    !> sloupcova permutace, vystupni parametr
    Integer(kind=ikind), dimension(:), intent(out)  :: q
    
    
    Integer(kind=ikind)  :: n,n1,n2,i,j,k,imax,jmax
    Real(kind=rkind) :: wrk

    ! zjistim rozmery matice
    n1 = ubound(A,1)
    n2 = ubound(A,2)
    !print *, "n1,n2=", n1,n2
    if ( n1 == n2) then ! je to ctvercova matice, jdu na to
      n = n1
      p = (/ (i,i=1,n) /) ! priprava zaznamu o permutaci
      q = p               ! napred zadna
      !print *,p
      do i=1,n-1
        call fmax(A,i,imax,jmax) ! najdi nejvetsi prvek
        n1     = p(imax)         ! zaznamenej presun
        p(imax)= p(i)
        p(i)   = n1
        n1     = q(jmax)
        q(jmax)= q(i)
        q(i)   = n1
        do j=1,n                 !a opravdu presun
          wrk       = A(i,j)
          A(i,j)    = A(imax,j)
          A(imax,j) = wrk 
        end do
        do j=1,n
          wrk       = A(j,i)
          A(j,i)    = A(j,jmax)
          A(j,jmax) = wrk 
        end do
        ! ted jiz opravdu rozkladej
        A(i+1:,i) = A(i+1:,i)/A(i,i)
        do j=i+1,n
          A(i+1:,j) = A(i+1:,j) - A(i,j)*A(i+1:,i)
        end do
      end do
    else
      print *, "neco divneho"
      stop  
    end if
  end subroutine LUdecomp


  subroutine fmax(A,i,imax,jmax)
    use typy
    Real(kind=rkind),dimension(:,:), intent(in) :: A
    Integer(kind=ikind), intent(in)             :: i
    Integer(kind=ikind), intent(out)            :: imax,jmax
  
    Real(kind=rkind)    :: maxv,wrk
    Integer(kind=ikind) :: j, k, n
 
    n = ubound(A,1)
    maxv =  abs (A(i,i))
    imax = i
    jmax = i
    do j=i,n
      do k=1,n
        wrk = abs(A(j,k))
        if ( wrk > maxv) then
          imax = j
          jmax = i
          maxv = wrk
        end if
      end do
    end do
  end subroutine fmax

  subroutine LUmul(LU, A)
    use typy
    Real(kind=rkind), dimension(:,:), intent(in) :: LU
    Real(kind=rkind), dimension(:,:), intent(out) :: A

    Real(kind=rkind)    :: wrk
    Integer(kind=ikind) :: i,j,k,n,kmax

    A = 0    
    n = ubound(A,1)
    do i=1,n
      do j=1,n
        wrk = 0
        if (i<=j) then
          kmax = i-1
          wrk = LU(i,j)
        else
          kmax = j
        end if
        do k=1,kmax
          wrk = wrk + LU(i,k)*LU(k,j)
        end do
        A(i,j) = wrk
      end do
    end do
  end subroutine LUmul


end module LinAlg
