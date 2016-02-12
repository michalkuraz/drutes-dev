module ADE_reader
  public :: ADE_read
  
  contains

    subroutine ADE_read()
      use typy
      use globals
      use global_objs
      use core_tools
      use ADE_globals
      use readtools

      integer :: i_err
      integer(kind=ikind) :: i, n
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      character(len=4096) :: msg

      call find_unit(file_contaminant, 200)
      open(unit = file_contaminant, file="drutes.conf/ADE/contaminant.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
	print *, "missing drutes.conf/ADE/contaminant.conf file"
	ERROR STOP
      end if

      allocate(adepar(maxval(elements%material)))
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/water.conf/matrix.conf  &
	the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_contaminant, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
	errmsg=trim(msg))
      
      write(unit=msg, fmt="(a)") "HINT 1: Is the molecular diffusion value positive?", new_line("a"), &
			      "HINT 2 : Is the number of molecular diffusion values corresponding to the amount of layers?"
			      
      do i=1, ubound(adepar,1)
	call fileread(adepar(i)%difmol, file_contaminant, ranges=(/0.0_rkind, 1.0_rkind*huge(tmp)/), &
	  errmsg=trim(msg))
      end do
 
      write(unit=msg, fmt="(a)") "HINT 1: Are all values anisotropy defining anisotropical diffusivity positive? ", new_line("a"), &
	"HINT 2 : Have you defined enough values for anisotropy &
	(e.g. for 2D define angle and the maximal and minimal value of diffusivity, in total 3 values)?", new_line("a"),&
	"HINT 3: The number of lines with diffusivity has to correspond to the number of materials & 
	defined by your mesh"
      
      
      allocate(tmp_array(drutes_config%dimen + 1))
      do i=1, ubound(adepar,1)
	allocate(adepar(i)%diff_loc(drutes_config%dimen))
	call fileread(r=tmp_array, fileid=file_contaminant, ranges=(/0.0_rkind, huge(tmp)/), errmsg=trim(msg))
	adepar(i)%anisoangle = tmp_array(1)
	adepar(i)%diff_loc = tmp_array(2:drutes_config%dimen + 1)
	allocate(adepar(i)%diff(drutes_config%dimen, drutes_config%dimen))
	call set_tensor(adepar(i)%diff_loc, (/adepar(i)%anisoangle/), adepar(i)%diff)
      end do
      
      
      

    end subroutine ADE_read				


  

end module ADE_reader