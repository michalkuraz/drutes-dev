! grider.f90
! 
! do grid 1D pridat geometry checks,

    if (abs(delta_x(ubound(delta_x,1),3) -length) > reps*abs(abs(delta_x(ubound(delta_x,1),3))-abs(length))) then
      print *, "check grid geometry info"
      STOP
    end if
    
    
    
    if (ubound(delta_x,1) > 1) then
      do i = 2, ubound(delta_x,1)
        if (abs(delta_x(i-1,3) - delta_x(i,2)) > reps*abs(abs(delta_x(i-1,3)) - abs(delta_x(i,2)))) then
          print *, "BaD c0nfiG FILE!!"
          print *, "check grid geometry info at mesh density line", i
	  print *, "the previous interval ends with", delta_x(i-1,3), "and the next one begins with", delta_x(i,2)
	  print *, "those values must be equal!"
          STOP
        end if
        if (abs(delta_x(i,2)) - abs(delta_x(i,3)) > reps) then
            print *, "BaD c0nfiG FILE!!"
            print *, "check grid geometry info at mesh density line", i
	    print *, "the lower interval value is higher then the upper interval value"
            STOP
        end if
      end do
    end if