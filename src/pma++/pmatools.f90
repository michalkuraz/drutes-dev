!> \file pmatools.f90
!> \brief soubor s modulem utilitek

!> nektere nastroje
module pmatools
    character, private :: dirsep=""

    public :: getyesno
    public :: pockej
    public :: get_free_unit
    public :: safe_open_write
    public :: safe_open_read
    public :: read_line
    public :: FileToStr
    public :: FindSubstr

    interface ReadInt
        procedure ReadIntI, ReadIntL
    end interface

    contains
    !> vytiskne zpravu a ziska odpoved
    subroutine getyesno(message, choice)
        !> zprava
        character(len=*), intent(in) :: message
        !> odpoved
        logical, intent(out) :: choice
        character(len=1)     :: z
        do
            print *, message
            print *,"zvol y nebo n (y..yes, n .. no)"
            read(unit=*,fmt=*) z
            if (z == "y" .or. z == "Y") then
                choice = .true.
            exit
            else if (z == "n" .or. z == "N") then
                choice = .false.
            exit
            else
                print *, "chybna volba"
            end if
        end do
    end subroutine getyesno

    !> pozastavi vypocet a pripadne zastavi
    subroutine pockej(message)
        implicit none
        !> zprava
        character(len=*),optional, intent(in) :: message
        logical :: volba
        if (present(message)) then
        call getyesno(message  // "Preruseni, pokracovat?",volba)
        else
            call getyesno("Preruseni, pokracovat?",volba)
        end if
        if (.not. volba) stop "Program prerusen v pause"
    end subroutine pockej


    !> nalezne prvni volne cislo zarizeni
    !!
    !! v soucasnosti je lepsi uzit newunit
    !<
    subroutine get_free_unit(fil)
        !> hledane cislo zarizeni
        integer, intent(out) :: fil
        logical :: ex,op

        fil = 101
        do
            inquire(unit=fil,exist=ex,opened=op)
            if ( .not. ex) then
            fil = fil + 1 !jednotka vubec neexistuje, tak dalsi
            else
                if ( op ) then
                fil = fil + 1 !existuje, ale je otevrena
                else
                    exit  !existuje a neni otevrena
                end if
            end if
        end do
    end subroutine get_free_unit

    !> Otevre soubor k zapisu
    !!
    !! V pripade pokusu o prepis muze pozadat o schvaleni.
    !! Rozhoduje o tom hodnota parametru mode
    subroutine safe_open_write(name, fil, mode)
        !> jmeno oteviraneho souboru , muze se zmenit, pokud uz existuje
        character(len=*), intent(in out)  :: name
        character(len=200)            :: name1
        !> cislo jednotky
        integer, intent(in) :: fil
        !> mod chovani (neni-li zadan plati 0)
        !! - 0 neptej se a skonci
        !! - 1 pozadej o pomoc
        !! - 2 neptej se a prepis
        integer, intent(in), optional  :: mode
        integer :: mode1
        logical :: ex, choice

        if (present(mode)) then
        mode1 = mode
        else
            mode1 = 0
        end if
        name1 = name
        do
            inquire(file=name1, exist=ex)
            if ( ex ) then
                Md: if (mode1==0) then
                STOP "soubor uz existuje, neprepisuji"
                else if (mode1==1) then
                    call getyesno(&
                    " soubor s udanym jmenem existuje. Prepsat?",choice)
                    CHC: if (choice) then
                        !! tohle by byt melo
                        !open(unit=fil, file=name1,status="replace",&
                        !action="write")
                        !! zatim  nahradime
                        call unlink(name1)
                        open(unit=fil, file=name1,action="write")
                    exit
                    else
                        print *, "zadej jmeno souboru"
                        read(unit=*,fmt=*) name1
                    end if CHC
                else
                    !print *, " podivne"
                    open(unit=fil, file=name1,status="replace",&
                    action="write")
                    exit
                end if Md
            else
                open(unit=fil, file=name1,status="new", action="write")
                exit
            end if
        end do
        name = name1
    end subroutine safe_open_write

    !> Otevre soubor na cteni
    !!
    !! V pripade neexistence zkusit se zeptat
    !! Rozhoduje o tom hodnota parametru mode
    subroutine safe_open_read(name, fil, mode)
        !> jmeno oteviraneho souboru , muze se zmenit, pokud neexistuje
        character(len=*), intent(in out)  :: name
        character(len=200)            :: name1
        !> cislo jednotky
        integer, intent(in) :: fil
        !> mod chovani (neni-li zadan plati 0)
        !! - 0 neptej se a skonci
        !! - 1 pozadej o pomoc
        integer, intent(in), optional  :: mode
        integer :: mode1
        logical :: ex, choice

        if (present(mode)) then
        mode1 = mode
        else
            mode1 = 0
        end if
        name1 = name
        do
            inquire(file=name1, exist=ex)
            if ( .not. ex ) then
                Md: if (mode1==0) then
                STOP "soubor neexistuje - nelze pokracovat"
                else if (mode1==1) then
                    call getyesno(&
                    "soubor s udanym jmenem neexistuje. Zadat jine jmeno?",choice)
                    CHC: if (.not. choice) then
                    Stop "soubor neexistuje - nelze pokracovat"
                    else
                        print *, "zadej jmeno souboru"
                        read(unit=*,fmt=*) name1
                        print *,"jmeno:", name1
                    end if CHC
                else
                    open(unit=fil, file=name1,action="read")
                    print *,"otevreno"
                    exit
                end if Md
            else
                print *,"oteviram"
                open(unit=fil, file=name1, action="read")
                exit
            end if
        end do
        name = name1
    end subroutine safe_open_read

    !> precte jednu radku ze souboru pripojeneho k jednotce unit.
    !! prodpoklada se soubor otevreny s encoding='UTF-8'
    !! \return fdf
    !<
    function read_line(unit,ierr) result(txt)
        use typy
        implicit none

        !> cislo jednotky s pripojenym souborem
        integer, intent(in) :: unit
        !> nepovinna promenna pro chybovou zpravu
        !! - 0 ... OK
        !! - 1 ... konec souboru
        !! - 2 ... chyba
        integer, intent(out), optional :: ierr
        !> precteny obsah radku
        character(len=:,kind=chkind), allocatable :: txt

        character(len=2000, kind=chkind) :: wrk
        character(len=2048) :: msg
        integer :: ln,aln, chyba

        ln = 0
        txt = ""
        do
            read(unit=unit,advance='no',fmt="(a2000)",eor=1, err=2,end=3,&
            size=aln, iostat=chyba, iomsg=msg) wrk
            txt = txt // wrk
            ln = ln + aln
            !print *, txt
            cycle
            1 continue
            !print *, "konec radku"
            txt = txt // wrk(1:aln)
            ln = ln+aln
            !print *, "celkem ", ln, " znaku"
            !print *, "text: ",txt
            if (present(ierr) ) ierr = 0
            return
        end do
        3 continue
        !print *, "konec souboru"
        if (present(ierr)) ierr = 1
        !print *, msg
        return

        2 print *, "chybicka se vloudila"
        if (present(ierr)) ierr = 2
        print *, chyba, msg
    end function read_line

    function FileToStr(name) result(Ret)
        use typy
        implicit none
        character(Len=*), intent(in) :: name
        character(len=:), allocatable :: ret, wrk, name1
        integer :: ierr, fil

        name1= name
        call safe_open_read(name1,fil)
        ret = ""
        do
            wrk = read_line(fil,ierr)
            if ( ierr == 2 ) then
            print *,"chybicka se vloudila"
            else if (ierr == 1) then
                exit ! doslo se na konec souboru
            end if
            ret = ret // wrk
        end do
    end function FileToStr

    function FindSubstr(Line, Str, mode) result(Pos)
        use typy
        implicit none
        !> prohledavany radek
        character(kind=chkind, len=*), intent(in) :: Line
        !> hledany retezec
        character(kind=chkind, len=*), intent(in) :: Str
        !> kam unisti nalezeny inndex
        !! - 0 ... pred prvni znak nalezenoho retezce
        !! - 1 ... za posledni znak nalezeneho retezce
        integer, optional, intent(in) :: mode
        !> navratova hodnota
        !! - nezaporna hodnota nalezena pozice
        !! - -1 ... nalezeny retezec nebyl nalezen
        !! !!! pozor: nalezena pozice muze byt i o jedna nizsi, nebo
        !! vyssi nez jsou rozmery retezce


        integer :: Pos, md

        if (present(mode)) then
        md = mode
        else
            md = 0
        end if

        Pos = Index(Line,Str)
        if ( Pos == 0) then
            Pos = -1
            return
        end if
        if (md == 0 ) then
        Pos = Pos - 1
        else
            Pos = Pos + Len(Str)
        end if
    end function FindSubstr


    subroutine ReadDir(StrList,Prefix,Files)
        use typy
        implicit none
        type(StringList), intent(out) :: StrList
        character(len=*), intent(in), optional :: Prefix
        logical, optional, intent(in) :: Files
        integer :: fil, i, ierr, cnt, delka, status
        character(len=:), allocatable :: nm,pr1
        logical :: fl

        if (present(Prefix)) then
        pr1 = Prefix
        else
            pr1 = ""
        end if
        if (present(Files)) then
        fl = Files
        else
            fl = .true.
        end if

        if (fl) then
        call system("dir /B " // pr1 // "> lst.txt",status)
        else
            call system("dir /B /AD " // pr1 // "> lst.txt",status)
        end if
        print *,"status",status
        if (status==0) then
            open(NEWUNIT=fil, file="lst.txt")
            cnt = 0
            do
                nm = read_line(fil,ierr)
                if (ierr==1) exit
                cnt = cnt + 1
                print *,cnt, nm
            end do
            rewind(fil)
            StrList%pocet = cnt
            allocate(StrList%name(1:cnt))
            allocate(StrList%delka(1:cnt))
            do i=1,cnt
                nm = read_line(fil,ierr)
                delka = len(nm)
                StrList%name(i)(1:delka) = nm
                StrList%delka(i) = delka
            end do
        close(fil)
        else
            print *,"prvni pokus nevysel, zkusim linuxovou variantu"
            dirsep = "/"
            call system("ls -aF " // pr1 // "> lst.txt",status)
            if (status == 0) then
                open(NEWUNIT=fil, file="lst.txt")
                cnt = 0
                do
                    nm = read_line(fil,ierr)
                    if (ierr==1) exit
                    if (nm(len(nm):len(nm))/="/") then
                        if (fl) then
                            cnt = cnt + 1
                            !print *,cnt, nm
                        end if
                    else
                        if (.not. fl) then
                            cnt = cnt + 1
                            !print *,cnt, nm
                        end if
                    end if

                end do
                rewind(fil)
                StrList%pocet = cnt
                !print *, "pocet radku=",cnt
                allocate(StrList%name(1:cnt))
                allocate(StrList%delka(1:cnt))
                i = 0
                do
                    nm = read_line(fil,ierr)
                    delka = len(nm)

                    !print *,i,delka, nm
                    if(nm(delka:delka)/="/") then
                        if (fl) then
                            i = i+1
                            StrList%name(i)(1:delka) = nm
                            StrList%delka(i) = delka
                        end if
                    else
                        if (.not. fl) then
                            i = i+1
                            StrList%name(i)(1:delka-1) = nm
                            StrList%delka(i) = delka-1
                        end if
                    end if
                    if (i==cnt) exit
                end do
                close(fil)

            else
                stop "neuspech pri cteni adresare"
            end if
        end if
    end subroutine ReadDir




    subroutine ReadIntI(I,low1,high1,fil1)
        use typy
        implicit none
        integer, intent(out) :: I
        integer, intent(in), optional ::low1, high1
        integer ::low, high
        integer, intent(in), optional :: fil1
        integer :: fil

        fil = 5
        low = -huge(low)
        high = huge(high)
        if (present(fil1)) fil = fil1
        if (present(low1)) low = low1
        if (present(high1)) high = high1
        do
            read(unit=fil, fmt=*, err = 5) I
            if ((I<low).or. (I>high)) goto 5
            return
            5 print *, "ocekavam cele cislo v intervalu <",low,",",high,">"
        end do
    end subroutine ReadIntI

    subroutine ReadIntL(I,low1,high1,fil1)
        use typy
        implicit none
        integer(kind=likind), intent(out) :: I
        integer(kind=likind), intent(in), optional ::low1, high1
        integer(kind=likind) ::low, high
        integer, intent(in), optional :: fil1
        integer :: fil

        fil = 5
        low = -huge(low)
        high = huge(high)
        if (present(fil1)) fil = fil1
        if (present(low1)) low = low1
        if (present(high1)) high = high1
        do
            read(unit=fil, fmt=*, err = 5) I
            if ((I<low).or. (I>high)) goto 5
            return
            5 print *, "ocekavam cele cislo v intervalu <",low,",",high,">"
        end do
    end subroutine ReadIntL


    function dot_product(x,y) result(r)
        use typy
        implicit none
        real(kind=rkind), dimension(:), intent(in) :: x,y
        real(kind=rkind) :: r
        integer(kind=ikind) :: i
        r = 0
        do i = LBOUND(x,1), UBOUND(x,1)
            r = r + x(i)*y(i)
        end do
    end function dot_product



    subroutine tools_test
    end subroutine tools_test



end module pmatools
