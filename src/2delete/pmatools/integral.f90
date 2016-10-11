module  integral
  use typy

  ! data jen pro privatni pouziti  
  type, public :: formulka
    real(kind=rkind), dimension(:), pointer :: a,h
  end type formulka
  ! typ pro data o integracni formule
  ! a ... uzly
  ! h ... vahy

  type(formulka), dimension(:), pointer, private :: group => null()
  ! obsahuje Gaussovy-Legenderovy integracni formule
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! verejne sluzby
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  public :: init_integration_formulas
  ! inicializuje seznam formuli
  ! volani: call init_integration_formulas(n)
  ! vstup : n ... pozadovany pocet predpripravenych formuli
  public :: getform
  ! vrati Gaussovu-Legenderovu integracni formuli stupne n
  ! vstup : n .... stupen formule
  ! vystup: a .... uzly formule
  !         h .... vahy formule
  public :: integ1, integ, integs
  private :: setup, dohledej, integ2, horner
  integer, private :: prep_forms = 0
  
  contains

    subroutine init_integration_formulas(n,ierr)
      use typy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   parametry
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      integer(kind=ikind), intent(in)  :: n
      integer(kind=ikind), intent(out), optional :: ierr 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! n ... naximalni rad generovane formule
      ! ierr .. chybove hlaseni
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! lokalni promenne
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(kind=ikind)  :: i, ierr1
      real(kind=rkind), dimension(:), pointer :: wrk
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! vypocet
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(group(1:n))
      do i=1,n
        call setup(i,group(i)%a, group(i)%h, ierr1)
        if (ierr1 > 0 )then
          exit
        end if
      end do
      prep_forms = n
      if (present(ierr)) then
        ierr = ierr1
      end if
    end subroutine init_integration_formulas
        


    pure subroutine getform(n,a,h,ierr)
      use typy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   parametery 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rkind), dimension(:), intent(out) :: a,h
      integer(kind=ikind), intent(in) :: n
      integer(kind=ikind), intent(out), optional :: ierr
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! n ... pocet bodu kvadr. vzorce - vstupni parametr 
      ! a ... uzly kvadraturniho vzorce
      ! h ... vahy kvdraturniho vzorce
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! lokalni promenne
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(kind=ikind) :: ierr1
      real(kind=rkind), dimension(:), pointer :: a1,h1
      
      if ((n>=1) .and. (n<=prep_forms)) then
        a = group(n)%a
        h = group(n)%h
        ierr1 = 0
      else
        call setup(n,a1,h1,ierr1)
        a = a1
        h = h1
        deallocate(a1,h1)
      end if		
      if ( present(ierr) ) then
        ierr = ierr1
      end if
    end subroutine getform


    !--integrace s pevnou formuli, neadaptivni integrace
    function integ1(a,b,f,n) result(cislo)
      use typy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! parametry
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rkind), intent(in) :: a,b
      integer(kind=ikind), intent(in) :: n
      interface
        function f(x) result(y)
          use typy
          real(kind=rkind), intent(in) :: x
          real(kind=rkind) :: y
        end function
      end interface
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! a ... dolni mez intervalu integrace
      ! b ... horni mez intervalu integrace
      ! n ... stupen kvadratueniho vzorce
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rkind) :: cislo, wrk
      integer(kind=ikind) :: i,m
      real(kind=rkind), dimension(:), allocatable :: a_temp, h_temp
          cislo = 0

      m = (n+1)/2

      allocate(a_temp(m), h_temp(m))

      call getform(n,a_temp, h_temp)

      wrk = (b+a)/2
      a_temp = (b-a)/2*a_temp

      if (modulo(n,2_ikind)==0) then
          do i=1,m
              cislo = cislo + h_temp(i)*(f(wrk+a_temp(i))+f(wrk-a_temp(i)))
          end do
      else
          do i=1, m-1
              cislo = cislo + h_temp(i)*(f(wrk+a_temp(i))+f(wrk-a_temp(i)))
          end do
          cislo = cislo + h_temp(m)*f(wrk)
      end if

          cislo = cislo*(b-a)/2


      deallocate(a_temp, h_temp)


    end function integ1

  !--integrace s pevnou formuli, neadaptivni integrace
    pure function integ2(a,b,f,n) result(cislo)
      use typy
      real(kind=rkind), intent(in) :: a,b
      integer(kind=ikind), intent(in) :: n

      interface
          pure function f(x) result(y)
          use typy
          real(kind=rkind), intent(in) :: x
          real(kind=rkind) :: y
          end function
      end interface

      real(kind=rkind) :: cislo, wrk
      integer(kind=ikind) :: i,m
      real(kind=rkind), dimension(:), allocatable :: a_temp, h_temp
          cislo = 0

      m = (n+1)/2

      allocate(a_temp(m), h_temp(m))

      call getform(n,a_temp, h_temp)

      wrk = (b+a)/2
      a_temp = (b-a)/2*a_temp

      if (modulo(n,2_ikind)==0) then
          do i=1,m
              cislo = cislo + h_temp(i)*(f(wrk+a_temp(i))+f(wrk-a_temp(i)))
          end do
      else
          do i=1, m-1
              cislo = cislo + h_temp(i)*(f(wrk+a_temp(i))+f(wrk-a_temp(i)))
          end do
          cislo = cislo + h_temp(m)*f(wrk)
      end if

          cislo = cislo*(b-a)/2


      deallocate(a_temp, h_temp)


    end function integ2


    pure recursive subroutine integs(a,b,f,n,rel_chyba, i_err, counter, hloubka, ilvl,&
                                cislo)
      use typy
      !--parametry
      real(kind=rkind), intent(in) :: a,b
      real(kind=rkind), intent(in),     optional :: rel_chyba
      integer(kind=ikind), intent(in),  optional :: n
      integer(kind=ikind), intent(out), optional :: i_err
      integer(kind=ikind), intent(in ), optional :: counter
      integer(kind=ikind), intent(in ), optional :: hloubka
      integer(kind=ikind), intent(in ), optional :: ilvl
      interface
          pure function f(x) result(y)
          use typy
          real(kind=rkind), intent(in) :: x
          real(kind=rkind) :: y
          end function
      end interface
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! a         .... dolni mez intervalu integrace  *povinny argument*
      ! b         .... horni mez intervalu integrace  *povinny argument*
      ! f         .... integrovana funkce             *povinny argument*
      ! n         .... stupen uziteho Gaussova vzorce      -nepovinny argument
      !                default = 4
      ! rel_chyba .... pozadovana relativni chyba vysledku -nepovinny argument
      !                default = 1.0e-5              
      ! i_err     .... chybove hlaseni                     -nepovinny argument
      !                 0 ... O.K.
      !                -1 ... nepovedlo se
      !                POZOR pokud neni zadan, tak pri chybe dojde k ukonceni
      !                      programu
      ! counter   .... uroven ve ktere se nachazim         -nepovinny argument
      !                default = 0 
      ! hloubka   .... maximalni povoleny pocet urovni     -nepovinny argument
      !                default = 20 
      ! ilvl      .... podrobnost informaci                -nepovinny argument
      !                0 ... nic nepsat
      !                1 ... vypisovat kazdy interval
      !                default = 0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---- navratova hodnota
      real(kind=rkind), intent(out) :: cislo

      !---pomocny cisla
      real(kind=rkind) ::  wrk1, wrk2,rel_chyba_tmp, wrk3
      real(kind=rkind), dimension(:), allocatable :: a1,h1
      integer(kind=ikind) :: i,m, i_err1, i_err2, counter1, counter2, cnt,&
                            hloubka_tmp, n_tmp, ierr_tmp,ilvl_tmp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! vyridime nepovinne parametry
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (present(n))  then
          n_tmp = n 
      else
          n_tmp = 4
      end if

      if (present(rel_chyba))  then
          rel_chyba_tmp = rel_chyba 
      else
          rel_chyba_tmp = 1.0e-5_rkind
      end if

      if (present(ilvl))  then
          ilvl_tmp = ilvl 
      else
          ilvl_tmp = 0
      end if

      if (present(hloubka))  then
          hloubka_tmp = hloubka 
      else
          hloubka_tmp = 20
      end if

      if (present(counter))  then
          cnt = counter + 1
      else
          cnt = 0
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! vyrizeno
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if ( ilvl_tmp > 0 ) then
      !  print *, "interval", a,b, cnt
      end if
      allocate(a1(1:(n_tmp+1)/2),h1(1:(n+1)/2))
      call getform(n_tmp,a1,h1)
      wrk1 = integ2(a,b,f,n_tmp)
      wrk2 = integ2(a,(a+b)/2,f,n_tmp) + integ2((a+b)/2,b,f,n_tmp)
      if (abs(wrk1 - wrk2) < rel_chyba*(abs(wrk1) + abs(wrk2)) ) then
          cislo = wrk2
          ierr_tmp = 0
      else
          if (cnt ==  hloubka_tmp) then 
              cislo = wrk2
              ierr_tmp = -1   
          else
            call   integs(a,(a+b)/2,f,&
                        n         = n_tmp,&
                        rel_chyba = rel_chyba_tmp,&
                        i_err     = i_err1,&
                        counter   = cnt,&
                        hloubka   = hloubka_tmp,&
                        ilvl      = ilvl_tmp,&
                        cislo     = wrk2)
            call   integs((a+b)/2,b,f,&
                        n         = n_tmp,&
                        rel_chyba = rel_chyba_tmp,&
                        i_err     = i_err2,&
                        counter   = cnt,&
                        hloubka   = hloubka_tmp,&
                        ilvl      = ilvl_tmp,&
                        cislo     = wrk3)
            wrk2 = wrk2 + wrk3             
            if (i_err1 + i_err2 == 0) then
              cislo = wrk2
              ierr_tmp = 0
            else
              cislo = wrk2
              ierr_tmp = -1
            end if
          end if
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ted jeste vyresit chybove hlaseni
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (present(i_err)) then
        i_err = ierr_tmp
      else
        if ( ierr_tmp /= 0) then
      !   stop "chyba v integrovani"
        end if
      end if  
    end subroutine integs





    !--integrace s pohyblivou formuli
    function integ(a,b,f,n,rel_chyba, counter, hloubka, ilvl)&
      result(cislo)
      use typy
      !--parametry
      real(kind=rkind), intent(in) :: a,b
      real(kind=rkind), intent(in),     optional :: rel_chyba
      integer(kind=ikind), intent(in),  optional :: n
      integer(kind=ikind), intent(in ), optional :: counter
      integer(kind=ikind), intent(in ), optional :: hloubka
      integer(kind=ikind), intent(in ), optional :: ilvl
      interface
          pure function f(x) result(y)
          use typy
          real(kind=rkind), intent(in) :: x
          real(kind=rkind) :: y
          end function
      end interface
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! a         .... dolni mez intervalu integrace  *povinny argument*
        ! b         .... horni mez intervalu integrace  *povinny argument*
        ! f         .... integrovana funkce             *povinny argument*
        ! n         .... stupen uziteho Gaussova vzorce      -nepovinny argument
        !                default = 4
        ! rel_chyba .... pozadovana relativni chyba vysledku -nepovinny argument
        !                default = 1.0e-5              
        ! i_err     .... chybove hlaseni                     -nepovinny argument
        !                 0 ... O.K.
        !                -1 ... nepovedlo se
        !                POZOR pokud neni zadan, tak pri chybe dojde k ukonceni
        !                      programu
        ! counter   .... uroven ve ktere se nachazim         -nepovinny argument
        !                default = 0 
        ! hloubka   .... maximalni povoleny pocet urovni     -nepovinny argument
        !                default = 20 
        ! ilvl      .... podrobnost informaci                -nepovinny argument
        !                0 ... nic nepsat
        !                1 ... vypisovat kazdy interval
        !                default = 0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---- navratova hodnota
      real(kind=rkind) :: cislo
      integer(kind=ikind) :: i_err

      call integs(a=a,b=b,f=f,n=n,rel_chyba=rel_chyba,i_err=i_err,hloubka=hloubka,&
                  ilvl=ilvl, cislo=cislo)
    end function integ 








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! soukrome procedurky
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure subroutine setup(n,a,h,ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  slouzi k vytvorenu uzlu a koeficientu  pro 
    !  Gaussovu kvadraturu
      use typy
      use LinAlg


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! parrametry
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rkind), dimension(:), pointer :: a,h
      integer(kind=ikind), intent(in) :: n 
      integer(kind=ikind), intent(out), optional :: ierr 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    n ... pocet bodu kvadr. vzorce
      !    a ... uzly kvadraturniho
      !    h ... vahy kbadraturniho vzorce
      ! ierr ... chybova zprava - nepovinna

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! lokalni promenne
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rkind), allocatable, dimension(:) ::  p, b 
      ! p ... n-ty ortogonalni polynom,
      ! b ... prava strany soustavy, kterou budeme vyrabet
      real(kind=rkind) :: x,y, x0, y0, x1, y1, step
      integer(kind=ikind) :: i,k, ik,ik1, m, nn
      !indexy, m=n/2
      real(kind=rkind), allocatable, dimension(:,:) :: matice
      ! matice ... matice soustavy, pro polynom a vahy
      logical :: sud
      ! sud ... sudost/lichost n
      real(kind=rkind), allocatable, dimension(:,:) :: intervaly
      integer(kind=ikind) :: pocet, ierr1


      ierr1 = 0
      ! najdu n-ty ortogonalni polynom
      ! tvorim jen nenulove koeficienty
      if (modulo(n,2_ikind)==0) then ! je to sude
        sud=.true.
        m=n/2+1
        allocate(matice(n/2,n/2), b(n/2), p(n/2+1)) 
        !print *, "jsme v sude casti"
        do k=1,n/2
          do i=1,n/2
            matice(k,i)=2.0_rkind/(2*(i+k)-3)
          end do
          b(k)=-2.0_rkind/(n+2*k-1)
        end do
        call gem(matice,b,p(1:n/2))
        p(n/2+1)=1.0_rkind
        deallocate(matice,b)
      else 
        sud=.false.
        m=n/2+1
        if (n>1_ikind) then		
          allocate(matice(n/2,n/2), b(n/2), p(n/2+1))
          ! n/2 se rizne o pulku dolu
          do k=1,n/2
            do i=1,n/2
              matice(k,i)=2.0_rkind/(2*(i+k)-1)
            end do
            b(k)=-2.0_rkind/(n+2*k)
          end do
          call gem(matice,b,p(1:n/2))
          p(n/2+1)=1.0_rkind
          deallocate(matice,b)
            else !tim vyresime problem je-li p stupne jedna
          allocate(p(1))
          p(1)=1.0_rkind
        end if
      end if
      ! mam ortogonalni polynom

      !jdu hledat koereny
      if (sud) then
        allocate(a(m-1))
      else
        allocate(a(m))
        a(m)=0 ! v lichem pripade je 0 vzdy korenem
      end if
      !je naolokovano jdeme na koreny, najdu jejich ctverce
      select case(m)
        case(1)
          continue !je to nula a to uz je dano vyse
        case(2)
          a(1) = -p(1) !ve skutecnosti to je -p(1)/p(2) 
        case(3)
          a(1) = (-p(2)+sqrt(p(2)*p(2)-4*p(1)*p(3)))/(2*p(3))
          a(2) = p(1)/p(3)/a(1)
        case default
          allocate(intervaly(2,m-1))
          nn=10*m
          do 
            step = 1.0_rkind/nn
            pocet = 0
            x0 = step
            y0 = horner(p,x0)
            do i=1, nn
              x1 = x0
              y1 = y0
              x0 = x0+step
              y0 = horner(p,x0)
              if (y0*y1 <= 0) then
                pocet = pocet + 1
                intervaly(:,(m-pocet))= (/ x1, x0 /)
              end if
            end do	
            if (pocet == m-1) then 
              exit		
            else
              nn=5*nn
            end if
          end do	
          do i=1,m-1
            a(i)=dohledej(intervaly(1,i),intervaly(2,i),p)
          end do
      end select
      ! ted to odmocnime
      a = sqrt(a)
      ! a ted spocitam vahy
      if (modulo(n,2_ikind)==0) then !suda vetev
        allocate(matice(n/2,n/2))
        allocate(b(n/2))
        allocate(h(n/2))
        matice(1,:) = 1.0_rkind
        b = a*a
        do i=2, n/2
          matice(i,:) = matice(i-1,:)*b
        end do
        do i=1, n/2
          b(i) = 1.0_rkind/(2*i-1)
        end do
        call gem(matice,b,h)
      else  !licha vetev
        allocate(matice(n/2+1, n/2+1))
        allocate(b(n/2+1))
        allocate(h(n/2+1))
        matice(1,:) = 2.0_rkind
        matice(1,ubound(matice,2)) = 1.0_rkind
        b = a*a
        do i=2, ubound(matice,1)
          matice(i,:) = matice(i-1,:)*b
        end do
        do i=1, ubound(b,1)
          b(i) = 2.0_rkind/(2*i-1)
        end do
        call gem(matice,b,h)
      end if
      if (present(ierr) ) then
        ierr = ierr1
      end if
    end subroutine setup




    pure function dohledej(xl,xp,p) result (x)
      use typy

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111111111111111111111
      !parametry
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111111111111111111111
      real(kind=rkind), intent(in) :: xl,xp
      real(kind=rkind),dimension(:), intent(in) :: p
      real(kind=rkind) :: x
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! xl ... levy krajni bod intervalu s korenem
      ! xp ... pravy krajni bod intervalu s korenem
      ! p  ... polynom
      ! 
      ! x ... navratova hodnota - nalezeny koren 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rkind) :: x1,x2,x3,y1,y2,y3,y4,y5

      x1 = xl
      x2 = xp
      y1 = horner(p,x1)
      y2 = horner(p,x2)
      y4 = x2-x1
      do
        x  = (x2+x1)/2
        y3 = horner(p,x)
        if (y3*y1< 0) then
          x2 = x
          y2 = y3
        else
          x1 = x
          y1 = y3
        end if
        y5 = y4
        y4 = x2-x1
        if (y5==y4) then
          exit
        end if
      end do
    end function

    pure &
    !> simple horner's scheme (zero derivatives only)
    function horner(p,x) result(px)
      use typy
      integer(kind=ikind) :: i
      real(kind=rkind), intent(in) :: x
      real(kind=rkind) :: px
      real(kind=rkind), intent(in), dimension(:) :: p

      px=0

      do i=ubound(p,1), 1, -1
        px=px*x+p(i)
      end do

    end function horner


end module integral

