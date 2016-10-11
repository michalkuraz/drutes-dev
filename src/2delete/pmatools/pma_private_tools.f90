
module pma_private_tools
  public :: get_free_unit, safe_open_write, getyesno

  contains
    !> nalezne prvni volne cislo zarizeni
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
                  call getyesno(" soubor s udanym jmenem existuje. Prepsat?", choice)
                  CHC: if (choice) then
                    open(unit=fil, file=name1,status="replace", action="write")
                    exit
                  else
                    print *, "zadej jmeno souboru"
                    read(unit=*,fmt=*) name1
                  end if CHC
                else
                  print *, " podivne"
                  open(unit=fil, file=name1,status="replace", action="write")
                  exit
                end if Md
            else
                open(unit=fil, file=name1,status="new", action="write")
                exit
            end if
        end do
        name = name1
    end subroutine safe_open_write


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

end module pma_private_tools