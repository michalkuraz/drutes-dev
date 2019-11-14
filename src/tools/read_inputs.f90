! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

!> \file read_inputs.f90
!! \brief Reader for input global and mesh files.
!<



module read_inputs
  public :: read_global
  public :: read_1dmesh_int, read_2dmesh_int, read_2dmesh_t3d, read_2dmesh_gmsh
  public :: read_scilab



  contains

    subroutine read_global()
      use globals
      use globals1D
      use globals2D
      use typy
      use core_tools
      use readtools
      use pde_objs
      use objfnc

      integer(kind=ikind) :: i, j, n, ierr
      real(kind=rkind), dimension(3) :: tmp
      character(len=4096) :: filename
      character(len=8192) :: msg
      integer :: local, global
      character(len=256), dimension(11) :: probnames
      character(len=2) :: dimensions
      
      if (.not. www) then
        local = file_global
        global = file_global
      else
        local = file_global
        global = file_wwwglob
      end if

      
       write(msg, *) "Incorrect option for problem type, the available options are:", new_line("a"),  new_line("a"), new_line("a"),&
        "   RE = Richards equation, primary solution is total hydraulic head H",&
        new_line("a"), new_line("a"),  &
        "   REstd = Richards equation, primary solution is pressure head h, use for homogeneous porous media only, & 
            and for EXPERIMENTAL purposes ONLY! , otherwise set RE", &
        new_line("a"), new_line("a"),  &
        "   boussi = Boussinesq equation for sloping land (1877)", &
        new_line("a"),  new_line("a"), &
        "   ADE = advection dispersion reaction equation (transport of solutes)",   &
        new_line("a"),  new_line("a"), &
        "   Re_dual = Richards equation dual porosity with total hydraulic head", &
        new_line("a"),  new_line("a"), &
        "   heat = Heat conduction equation (Sophoclea, 1979)", &
        new_line("a"),  new_line("a"), &
        "   LTNE = Local Thermal Non-Equilibrium heat transport model (unfinished yet)", &
        new_line("a"),  new_line("a"), &
        "   freeze = Richards equation with freezing/thawing processes (unfinished yet)", &
        new_line("a"),  new_line("a"), &
        "   kinwave = Kinematic wave equation for 2D catchments", &
                new_line("a"),  new_line("a"), &
        "   vapour = Richards equation coupled with heat equation, phase changes due evaporation, vapour flow ", &
        new_line("a"),  new_line("a"), new_line("a")
        
	

      probnames(1) = "REtest"
      probnames(2) = "RE"
      probnames(3) = "boussi"
      probnames(4) = "ADE"
      probnames(5) = "Re_dual" 
      probnames(6) = "heat"
      probnames(7) = "REstd" 
      probnames(8) = "LTNE"
      probnames(9) = "freeze"
      probnames(10) = "kinwave"
      probnames(11) = "vapour"
      
	
      call fileread(drutes_config%name, local, trim(msg), options=probnames)

      call fileread(dimensions, local, options=(/"1 ","2 ","2r","3 "/))
      
      select case(cut(dimensions))
        case("1","2")
          read(dimensions,'(I10)') drutes_config%dimen
        case("2r")
          drutes_config%dimen=2
          drutes_config%rotsym=.true.
        case("3")
          msg="3D no yet implemented, you can use only 1D, 2D and 2D for rotational symmetric flow (pseudo 3D)."
          call file_error(file_global, msg)
      end select
      
      if (drutes_config%name == "boussi" .and. drutes_config%dimen > 1) then
        write(msg, fmt=*) 'You have selected Boussinesq equation, Boussinesq equation originates from Dupuit &
		   approximation, and so it is assumed for 1D only!!! &
		   '//NEW_LINE('A')//'    But your domain was specified for: ', drutes_config%dimen, "D"
        call file_error(local, msg)
      end if
      

      write(unit=msg, fmt=*) "INCORRECT mesh type definition", new_line("a"), &
       "the available options are: 1 - internal mesh generator" , new_line("a"), &
       "   (very simple uniform meshes, for debuging only", new_line("a"), &
       "                           2 - t3d mesh generator", new_line("a"), &
       "                           3 - gmsh mesh generator"


      call fileread(drutes_config%mesh_type, local,msg,ranges=(/1_ikind,3_ikind/))
      
      call fileread(max_itcount, local, ranges=(/1_ikind, huge(1_ikind)/), &
      errmsg="maximal number of iterations must be positive, greater than 1, &
	      and smaller than the maximal number your computer can handle :)")
      
      call fileread(iter_criterion, local, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
      errmsg="iteration criterion must be positive, and smaller than the maximal number your computer can handle :)")
      
      call fileread(time_units, local)
      
      call fileread(init_dt, local, ranges=(/tiny(0.0_rkind), huge(0.0_rkind)/),&
      errmsg="initial time step must be positive, and smaller than the maximal number your computer can handle :)")
      
      call fileread(end_time, local, ranges=(/init_dt, huge(tmp(1))/), &
      errmsg="end time must be greater than the minimal time step, and smaller than the maximal number your computer can handle :)")
      
      
      call fileread(dtmin, local, ranges=(/0.0_rkind, init_dt/), &
      errmsg="minimal time step must be positive, and smaller than the initial time step")
      
      call fileread(dtmax, local, ranges=(/dtmin, huge(tmp(1))/), &
       errmsg="maximal time step must be greater than the minimal time step, &
		  and smaller than the maximal number your computer can handle :)")
      
      
      write(msg, fmt=*) "methods for observation time print could be only", new_line('a'), &
	     "	1 - adjust time stepping to observation time values",  new_line('a'), &
	     "	2 - linearly interpolate solution between two consecutive solutions (recommended)"
      
      call fileread(observe_info%method, local, ranges=(/1_ikind, 2_ikind/), errmsg=msg)
      
      write(msg, fmt=*) "set correct name for the observation time outputs format", new_line('a'), &
      	"	scil - scilab output files",  new_line('a'), &
      	"	pure - just raw data with nodes IDs and FEM coefficients", new_line('a'), &  
      	"       gmsh - gmsh output files"
      
      call fileread(observe_info%fmt, global, options=(/"scil", "pure", "gmsh"/))
      
      
      call fileread(observe_info%anime, global)
      
      if (observe_info%anime) then
        call fileread(observe_info%nframes, global, ranges=(/1_ikind, huge(1_ikind)/), &
        errmsg="number of frames for animation must be greater than zero &
        and smaller than the maximal number your computer can handle :)")
      else
	     observe_info%nframes = 0
      end if
      
      call fileread(n, local)

      allocate(observe_time(n + observe_info%nframes))
      
      if (.not. observe_info%anime) then
	
        write(msg, *) "HINT 1: check number of the observation time values", new_line('a'), &
	       "HINT 2: You have selected [n] (not) to create animation frames, have you commented out the required number of frames?" , & 
	         new_line('a')
      else
        write(msg, *) "HINT: check number of the observation time values"
      end if
      
      
      do i=1, n
      	call fileread(observe_time(i)%value, local, errmsg=msg)
      end do
      

      ! reads the observation points number
      call fileread(n, local,ranges=(/0_ikind, huge(1_ikind)/), &
	     errmsg="number of observation times cannot be negative :)")
      
      allocate(observation_array(n))   
      
       if (.not. observe_info%anime) then
	
	        write(msg, *) "HINT 1: check number of the observation point values", new_line('a'), &
             "   HINT 2: check number of the observation points coordinates" , new_line('a'), &
             "   HINT 3: You have selected [n] (not) to create animation frames, have you commented out the &
                required number of frames?" , &
		          new_line('a')
      else
        write(msg, *)  "HINT 1: check number of the observation point values", new_line('a'), &
         "HINT 2: check number of the observation points coordinates" , new_line('a')
      end if     
      
      do i=1,n
      	allocate(observation_array(i)%xyz(drutes_config%dimen))
      end do
      
      
      do i=1,n
	     call fileread(observation_array(i)%xyz(:), local, errmsg=msg, checklen=.true.)
      end do  
      !----
      
      ! reads coordinates with measured points
      call fileread(n, global, errmsg="Incorrect number of points with measurement data.")

      allocate(measured_pts(n))
      do i=1, n
	       allocate(measured_pts(i)%xyz(drutes_config%dimen))
      end do
      
      do i=1, n
	       call fileread(measured_pts(i)%xyz(:), global, errmsg="HINT: check coordinates of the points with measurement data")
      end do  
      
      write(msg, *) "set correct value for the terminal outputs", new_line("a"), "     0 - standard output", &
    	new_line("a"),  "     1 - everything goes to out/screen.log", new_line("a"), &
  	"    -1 - everything goes to /dev/null (use only on Linux based systems (I have no idea what will happen in MAC OS X))"
      call fileread(print_level, global, ranges=(/-1_ikind, 1_ikind/), errmsg=msg)
      

      call fileread(i=drutes_config%it_method, fileid=local, ranges=(/0_ikind,2_ikind/),&
	       errmsg="you have selected an improper iteration method")
      
      write(msg, *) "Define method of time integration", new_line("a"), &
      "   0 - steady state problem", &
           new_line("a"), &
      "   1 - unsteady problem with lumped (diagonal) capacity matrix (recommended)", new_line("a"), &
      "   2 - unsteady problem with consistent capacity matrix"
      
      
      call fileread(pde_common%timeint_method, global, ranges=(/0_ikind,2_ikind/), errmsg=msg)
      
      call fileread(objval%compute, global, errmsg="Set correct value [y/n] for evaluating objective function")
      
      call fileread(drutes_config%check4mass, global, & 
        errmsg="Set correct value [y/n] for evaluating integral mass balance accuracy")
      
      call fileread(drutes_config%run_from_backup, global, errmsg="specify [y/n] if you want to relaunch your computation")
      
      
      if (drutes_config%run_from_backup) then
      	call fileread(backup_file, global, errmsg="backup file not specified")
      end if
      
      call fileread(integ_method, global)
      
      if (integ_method/10 < 1 .or. integ_method/10 > 12 .or. modulo(integ_method,10_ikind)/=0  &
          .or. integ_method == 100 .or. integ_method == 110) then
        write(msg, *)  "The recognized integration method code names are : & 
          10, 20, 30, 40, 50, 60, 70, 80, 90, 120", new_line("a") ,&
          "10 - for single Gauss quadrature node, 120 - for 12 points Gauss quadrature", new_line("a"), &
         "Your unrecognized setup was:", integ_method
        call file_error(global, msg)
      end if
      
!       call fileread(debugmode, global,"specify [y/n] for debugging (development option, not recommended)" )
        

    end subroutine read_global
    
  




    subroutine read_1dmesh_int()
      use typy
      use globals
      use globals1d
      use core_tools
      use simegen
      use readtools
      use debug_tools

      integer :: ierr
      integer :: n, i, j
      character(len=1) :: ch
      character(len=4096) :: msg

      call comment(file_mesh)
      read(unit=file_mesh, fmt = *, iostat = ierr) length_1D ; if (ierr /= 0) call file_error(file_mesh)



      call comment(file_mesh)
      read(unit=file_mesh, fmt = *, iostat = ierr) n; if (ierr /= 0) call file_error(file_mesh)


      allocate(deltax_1D(0:n, 3))

      do i=1, n
        call comment(file_mesh)
        read(unit=file_mesh, fmt=*, iostat = ierr) deltax_1D(i,1), deltax_1D(i,2), deltax_1D(i,3)
        if (ierr /= 0) call file_error(file_mesh)
      end do
      
      
     if (abs(length_1D-(deltax_1d(ubound(deltax_1d,1), 3)-deltax_1d(1,2))) > 10*epsilon(deltax_1d(1,2))) then

        write(unit=msg, fmt=*)"incorrect definition in either:" , new_line("a"),  &
        "                          - domain length definition" ,  new_line("a") ,&
       "    or" , new_line("a"),  &
        "                          - the mesh density description values" ,  new_line("a"), &
        " " ,  new_line("a") ,&
        "in file drutes.conf/mesh/drumesh1d.conf" ,  new_line("a") ,&
        " " , new_line("a") , &
       "the mesh density description must cover (and cannot overlap) the entire domain length",  new_line("a"), &
        "-----------------------------------------------------------------------------"
        call file_error(file_mesh, msg)
      end if
      
      call comment(file_mesh)

      read(unit=file_mesh, fmt=*, iostat=ierr) n ; if (ierr /= 0) call file_error(file_mesh)

      allocate(materials_1D(n,3))

      do i=1,n
       call comment(file_mesh)
       read(unit=file_mesh, fmt=* , iostat=ierr)  materials_1D(i,:) ; if (ierr /= 0) call file_error(file_mesh)
      end do
      

      do i=1, ubound(materials_1D,1)-1
        if (abs(materials_1D(i,3)-materials_1D(i+1,2)) > epsilon(materials_1D(i+1,1))) then
          write(msg,*) "ERROR! material description does not cover the entire domain, check material id:", int(materials_1D(i+1,1))
          call file_error(file_mesh, msg)
        end if
      end do
      
      
      if (maxval(abs(materials_1D)) < length_1D) then
        write(msg, fmt=*) "ERROR! material description does not cover the entire domain", new_line("a"), new_line("a"), &
        "1. Check 1D domain length value",  new_line("a"), &
        "2. Check your number of materials, maybe it's wrong"
        call file_error(file_mesh, msg)
      end if
	
      call simegen1D()

    end subroutine read_1dmesh_int


    subroutine read_2dmesh_int()
      use typy
      use globals
      use globals2d
      use core_tools
      use simegen
      use readtools

      integer :: ierr
      real(kind=rkind) :: width, length, density
      integer(kind=ikind) :: edges, i, proc
      real(kind=rkind), dimension(:,:,:), allocatable :: xy
      real(kind=rkind), dimension(2) :: center_coor
      

      call comment(file_mesh)
      read(unit=file_mesh, fmt=*, iostat = ierr) width; if (ierr /= 0) call file_error(file_mesh)

      call comment(file_mesh)
      read(unit=file_mesh, fmt=*, iostat = ierr) length; if (ierr /= 0) call file_error(file_mesh)
      call comment(file_mesh)

      read(unit=file_mesh, fmt=*, iostat = ierr) density; if (ierr /= 0) call file_error(file_mesh)

      call comment(file_mesh)
      
      read(unit=file_mesh, fmt=*, iostat = ierr) center_coor; if (ierr /= 0) call file_error(file_mesh)    
      
      call comment(file_mesh)
            
      read(unit=file_mesh, fmt=*, iostat=ierr) edges; if(ierr /= 0) call file_error(file_mesh)

      allocate(xy(edges,2,2))
      
      

      do i=1, edges
        call comment(file_mesh)
        read(unit=file_mesh, fmt=*, iostat = ierr) xy(i,1,:); if(ierr /= 0) call file_error(file_mesh)
        read(unit=file_mesh, fmt=*, iostat = ierr) xy(i,2,:); if(ierr /= 0) call file_error(file_mesh)
      end do
	
      call simegen2d(width, length, density, xy, center_coor)


      deallocate(xy)
    end subroutine read_2dmesh_int

!> loads mesh data files from external mesh generator and allocate required mesh structures
    subroutine read_2dmesh_t3d()
      use typy
      use globals
      use globals2D
      use core_tools
      use readtools
      use debug_tools


      integer(kind=ikind) :: i, okro, l
      real :: ch

 
  
      call comment(file_mesh)
      read(unit=file_mesh, fmt=*) ch
      call comment(file_mesh)
      read(unit=file_mesh, fmt=*) nodes%kolik, ch, elements%kolik
      call comment(file_mesh) 

      call mesh_allocater()

      do i=1, nodes%kolik
        call comment(file_mesh)
        read(unit=file_mesh, fmt=*) nodes%id(i), nodes%data(i,:), ch, ch, ch, nodes%edge(i) 
        if (nodes%edge(i) > 100) then
          CONTINUE
        else
          nodes%edge(i) = 0
        end if
      end do
      
   
    
      do i=1, elements%kolik
        call comment(file_mesh)
        read(unit=file_mesh, fmt=*) elements%id(i), elements%data(i,:), ch, ch, elements%material(i)
        elements%material(i) = elements%material(i) ! - 10000
      end do
     

      close(file_mesh)

    end subroutine read_2dmesh_t3d
    
    
    subroutine read_2dmesh_gmsh()
      use typy
      use globals
      use globals2D
      use core_tools
      use readtools
      use pde_objs
      use debug_tools
      use geom_tools
      
      character(len=256) :: ch
      character(len=11) :: msh
      real(kind=rkind) :: mshfmt
      integer(kind=ikind) :: i, itmp, jtmp, itmp2, edge, nd1, nd2, i_err, k, l, j, n, id, el
      integer(kind=ikind), dimension(:), allocatable :: tmp_array
      logical :: use_filemat=.false., update
      integer :: file_mat
      real(kind=rkind), dimension(:,:), allocatable :: domain, smalldom
      real(kind=rkind), dimension(2) :: point
      integer(kind=ikind), dimension(:), allocatable :: excl_ids
      
      
      
      read(unit=file_mesh, fmt=*, iostat=i_err) msh
     
      if (i_err /= 0 .or. msh /= "$MeshFormat") then
        print *, "ERROR: file drutes.conf/mesh/mesh.msh has strange (unsupported) syntax"
        print *, "exiting from read_inputs::read_2dmesh_gmsh"
        ERROR STOP
      end if
      
      read(unit=file_mesh, fmt=*, iostat=i_err) mshfmt
      
      
      if (i_err /= 0) then 
        print *, "ERROR: file drutes.conf/mesh/mesh.msh has strange (unsupported) syntax"
        print *, "exiting from read_inputs::read_2dmesh_gmsh"
        ERROR STOP
      end if
      
      if (int(mshfmt) /= 2) then
        print *, "DRUtES supports gmsh version 2 only,"
        print *,  "launch gmsh with option -format msh2. "
        ERROR STOP
      end if

      
      do
        read(unit=file_mesh, fmt=*) ch
        if (trim(ch) == "$Nodes") EXIT
      end do
      
      read(unit=file_mesh, fmt=*) nodes%kolik

      do 
        read(unit=file_mesh, fmt=*) ch
        if (trim(ch) == "$Elements") EXIT
      end do

      read(unit=file_mesh, fmt=*) itmp
      
      elements%kolik = 0
      
      do i=1, itmp
        read(unit=file_mesh, fmt=*) ch, itmp2
        if (itmp2 == 2) then
          elements%kolik = elements%kolik + 1
        end if
      end do
      
      call mesh_allocater()

            
      do 
        itmp = ftell(file_mesh)
        backspace(unit=file_mesh, iostat=i_err)
        jtmp = ftell(file_mesh)
        if(itmp == jtmp) EXIT
      end do
  
  
      do
        read(unit=file_mesh, fmt=*) ch
        if (trim(ch) == "$Nodes") then
          read(unit=file_mesh, fmt=*) ch
          EXIT
        end if
      end do
      
      do i=1, nodes%kolik
        read(unit=file_mesh, fmt=*) itmp,  nodes%data(i,:)
      end do
      
      do 
        read(unit=file_mesh, fmt=*) ch
        if (trim(ch) == "$Elements") then
          read(unit=file_mesh, fmt=*) ch
          EXIT
        end if
      end do
      
      jtmp = 0
      
      nodes%edge = 0
      
      do 
        read(unit=file_mesh, fmt=*, iostat=i_err) ch, itmp
        if (ch == "$EndElements") then
          EXIT
        end if
        backspace(file_mesh, iostat = i_err)
        if (i_err /= 0) then
          print *, "Something wrong with your GMSH input file"
          print *, "called from read_inputs::read_2dmesh_gmsh"
          STOP
        end if

        
        select case(itmp)
          case(1)
            read(unit=file_mesh, fmt=*) k, l, i, edge, j,  nd1, nd2
            if (i > 2) then
              call write_log("number of tags for edge must be equal 1, in mesh file it equals", int1=i)
              call write_log("update your GMSH input file!")
              allocate(tmp_array(i-1))
              backspace(file_mesh)
              read(unit=file_mesh, fmt=*) k, l, i, edge,  tmp_array,  nd1, nd2
              edge = tmp_array(1)
              call write_log(text="tags with positions greater than 1 were ignored")
              deallocate(tmp_array)
            end if

            nodes%edge(nd1) = edge 
            nodes%edge(nd2) = edge 
            
          case(2)
            jtmp = jtmp + 1
            read(unit=file_mesh, fmt=*) k, l, i,  j, elements%material(jtmp), elements%data(jtmp,:)
            if (i /= 2) then
              call write_log("number of tags for element must be equal 2")
              call write_log("update your GMSH input file!")
              if (i > pde_common%processes + 1) then
          
              call write_log("the number of tags is greater than the number of solution components")
              
              if (allocated(tmp_array)) deallocate(tmp_array)
              
                allocate(tmp_array(i-1))
                
                backspace(file_mesh)
                read(unit=file_mesh, fmt=*) k, l, i, tmp_array, elements%material(jtmp),  elements%data(jtmp,:)
                elements%material(jtmp) = tmp_array(1)
                call write_log(text="tags with position greater than 2 were ignored")
                deallocate(tmp_array)
              else
                call write_log("the number of tags is lower than the number of solution components, I don't know what to do :(")
                      print *, "called from read_inputs::read_2dmesh_gmsh"
                error STOP
              end if
              end if
                  
            case default
              print *, "your GMSH input file contains unsupported element type"
              print *, "called from read_inputs::read_2dmesh_gmsh"
              STOP
         end select
	 
      end do
      

    
 
    end subroutine read_2dmesh_gmsh
    
    subroutine read_scilab(name, proc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use core_tools
      use debug_tools
      
      character(len=*), intent(in)  :: name
      integer(kind=ikind), intent(in) :: proc
      integer :: fileid, ierr
      character(len=1) :: ch
      integer(kind=ikind) :: i,j,k,l
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(3,3) :: elloc
      integer(kind=ikind), dimension(3) :: pmt
      


      call find_unit(fileid,20)
      
      open(unit=fileid, file=trim(name), action="read", status="old", iostat=ierr)
      
      if (ierr /=0) then
        print *, "unable to open input file with scilab data, called from read_inputs::read_scilab"
        print *, "the specified location was:", trim(name)
        error stop
      end if
      
      read(unit=fileid, fmt=*) ch, time
      
      read(unit=fileid, fmt=*) ch, ch, i
      
      if (i/=elements%kolik) then
        print *, "the mesh in the scilab file has different number of elements then the mesh "
        print *, i, "/=", elements%kolik
        print *, "that is loaded or gerenerated from drutes config files"
        print *, "called from read_inputs::read_scilab"
        error stop
      end if
      

      
      read(unit=fileid, fmt=*) ch
      read(unit=fileid, fmt=*) ch
      read(unit=fileid, fmt=*) ch
      
      do i=1,elements%kolik
        do j=1,3
          read(unit=fileid, fmt=*) ch, ch, ch, ch, elloc(j,1)
        end do
        
        do j=1,3
          read(unit=fileid, fmt=*) ch, ch, ch, ch, elloc(j,2)
        end do
        
        do j=1,3
          read(unit=fileid, fmt=*) ch, ch, ch, ch, elloc(j,3)
        end do
        
        pmt = 0
        do j=1,3
          l=elements%data(i,j)
          do k=1,3
            if (norm2(nodes%data(l,:)-elloc(k,1:2)) < 10*epsilon(tmp) ) then
              pmt(k) = j
            end if
          end do
        end do
        
        if (minval(pmt) == 0) then
          print *, "the scilab file has a different mesh (probably boundary conditions)"
          print *, "than the one defined from drutes config files"
          print *, "called from read_inputs::read_scilab"
          ERROR STOP
        end if
	  
        pde(proc)%solution(elements%data(i,:)) = elloc(pmt,3)  
	  

      end do
      
      do i=1,4
        read(unit=fileid, fmt=*) ch
      end do


      close(fileid)    
      
    
    end subroutine read_scilab




      

end module read_inputs

