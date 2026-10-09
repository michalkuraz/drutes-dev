module init_netcdf

  contains
  
    subroutine netcdf()
      use typy
      use netcdf
      use ncglobvars
      use datetime
      use globals
      use global_objs
      use nctools
      use debug_tools
      use ncdem
      use netcdfflux
      use core_tools
      use pde_objs
      use readtools
      use ncfluxarea
      use ncmesh
      use ncmap
      use geom_tools
      use ncdispersion, only: parse_ls_dispersivity
      use ncsupg, only: read_adenc_supg
      use ncboundary, only: read_adenc_banks,prepare_adenc_banks
      use ncconservative, only: read_adenc_conservative
      use nchydroflow, only: read_adenc_hydro,LShydro,hydro_filter,hydro_initialize
      use ncrouting, only: read_adenc_routing,routing_apply_directions
            
      integer :: ierr, filetmp, fileconf
      integer(kind=ikind) :: i, bccnt, j
      logical :: success
      real(kind=rkind) :: q, segment_distance, best_segment_distance
      character(len=1024) :: errmsg
      character(len=4096) :: dispersivity_line
      integer(kind=ikind), dimension(3) :: datearray
      integer(kind=ikind), allocatable :: original_edge(:)
      real(kind=rkind), dimension(2) :: xy, A, B, C
      real(kind=rkind), dimension(4,2) :: pts
      
      starttime%year = 2010
      starttime%month = 1
      starttime%day = 1
      starttime%hour = 0
      starttime%minute = 0
      starttime%second = 0
      
      
      
      
      open(newunit=fileconf, file="drutes.conf/netcdf/netcdf.conf", action="read", status="old")
      
!      read_int_array(r, fileid, ranges, errmsg, checklen, noexit)

      call fileread(datearray, fileconf)
      
      if (datearray(1) < 1900) then
        write(errmsg, *) "----------------------------------", new_line("a"), &
                        "W: the year of your simulation start is very old, check drutes.conf/netcdf/netcdf.conf", new_line("a"), &
                        "the year you defined is:", datearray(1), new_line("a"), "check your inputs! however simulation continues",&
                        new_line("a"), "----------------------------------"
        call write_log(cut(errmsg))
      end if
      
      starttime%year = datearray(1)
      
      if (datearray(2) < 1 .or. datearray(2) > 12) then
        write(errmsg, *) "incorrect month specified, month should be within 1-12, you specified", datearray(2)
        call file_error(fileconf, errmsg)
      end if
     
      starttime%month = datearray(2)
      
      starttime%day = datearray(3)
    
      
      call comment(fileconf)
      read(fileconf, '(A)', iostat=ierr) dispersivity_line
      if (ierr /= 0) call file_error(fileconf, "missing ADEnc dispersivity record")
      call parse_ls_dispersivity(dispersivity_line, LSdisp, LSdisp_transverse, success)
      if (.not. success) call file_error(fileconf, &
        "ADEnc dispersivity must be one or two finite nonnegative values: alpha_L [m] alpha_T [m]")
      write(errmsg, *) "ADEnc dispersivity [m]: alpha_L=", LSdisp, " alpha_T=", LSdisp_transverse
      call write_log(trim(errmsg))
      call read_adenc_supg()
      call read_adenc_banks()
      call read_adenc_conservative()
      call read_adenc_hydro()
      
      
      call fileread(Qmin, fileconf, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
                    errmsg="incorrect Qmin definition in drutes.conf/netcdf/netcdf.conf")
    
      call fileread(cinit_ls, fileconf, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
                    errmsg="incorrect c initial definition in drutes.conf/netcdf/netcdf.conf")

      call fileread(channel_count, fileconf, ranges=(/1_ikind, huge(1_ikind)/), &
                    errmsg="number of channels in drutes.conf/netcdf/netcdf.conf must be at least 1")

      
      
      call set_adenc_boundary_id()
      original_edge=nodes%edge
      
      
      ierr = nf90_open(path="drutes.conf/netcdf/mRM_Fluxes_States.nc", mode=nf90_nowrite, ncid=netcdfID)
      
      if (ierr /= nf90_noerr) ERROR STOP "Error opening NetCDF file: drutes.conf/netcdf/mRM_Fluxes_States.nc"
    
      call read_ncorigin(netcdfID, ncstart)
      

      
      call ncflux_init(success, errmsg)
   
      if (.not. success) then
        print *, cut(errmsg)
        print *, "unsupported structure of the file drutes.conf/netcdf/mRM_Fluxes_States.nc"
        print *, "is this correct output from mHM?"
        print *, "after succesfull mHM simulation you should copy "
        print *, "      mRM_Fluxes_States.nc -> drutes.conf/netcdf/mRM_Fluxes_States.nc"
        ERROR STOP
      end if
      
      call read_ncbounds(netcdfID, success, errmsg)
      
      if (.not. success) then
        print *, cut(errmsg)
        print *, "unsupported structure of the file drutes.conf/netcdf/mRM_Fluxes_States.nc"
        print *, "is this correct output from mHM?"
        print *, "after succesfull mHM simulation you should copy "
        print *, "      mRM_Fluxes_States.nc -> drutes.conf/netcdf/mRM_Fluxes_States.nc"
        ERROR STOP
      end if
    
      call read_adenc_routing()
      call getmeshalt()
      
      ora_di_ini = difftime(ncstart, starttime, "hrs")
      
    
      do i=1, nodes%kolik
      
        call ncflux_get_xy_cell(nodes%data(i,1), nodes%data(i,2), ora_di_ini, q, success, errmsg)
        if (.not. success) then
          nodes%edge(i) = addedbc
        else if (q < Qmin) then
          nodes%edge(i) = addedbc
        end if
        
!         if (q > 0.0) then
!           if (nint(nodealt(i)) == missing) then
!             write(errmsg, *) "W: your dem model doesn't conver the entire watershed, update drutes.conf/netcdf/dem.nc, node:", i, &
!               "will be deactivated"
!             call write_log(errmsg)
!             nodes%edge(i) = addedbc
!           end if
!         end if
        
        
      end do

      allocate(ncfluxdata%activeel(elements%kolik))
      
      do i=1, elements%kolik
        bccnt = 0
        do j=1, ubound(elements%data,2)
          if (nodes%edge(elements%data(i,j)) == addedbc) then
            bccnt = bccnt + 1
          end if
        end do
        if (bccnt == ubound(elements%data,2)) then
          ncfluxdata%activeel(i) = .false.
        else
          ncfluxdata%activeel(i) = .true.
        end if
      end do
      
!       ncfluxdata%activeel(:) = .true.
      
      allocate(ncfluxdata%cellarea(elements%kolik))
      allocate(ncfluxdata%fluxvct(elements%kolik,2))
      
      ncfluxdata%fluxvct = 0.0_rkind
      
      
      
      do i = 1, elements%kolik
        xy(1) = avgarr(nodes%data(elements%data(i,:),1))
        xy(2) = avgarr(nodes%data(elements%data(i,:),2))
        call ncflux_cell_area_xy(xy(1), xy(2), ncfluxdata%cellarea(i), success, errmsg)
      end do
      
      
      call getncmesh(netcdfID, ncnodes, ncelements, success, errmsg)

      if (.not. success) then
        print *, trim(errmsg)
        error stop
      end if
      
      call ncelslope()
      
      call mapel()
      
      call terrain_slopes()
      ! terrain_slopes historically pins entire triangles with missing DEM to
      ! addedbc. Restore participating DOFs only AFTER that legacy annotation.
      ! ADEnc flux direction below comes from channel polylines, not DEM slopes.
      call prepare_adenc_banks(original_edge)
      
      call readchannel()
      
      do i=1, elements%kolik
        if (ncfluxdata%activeel(i)) then 
          C(1) = mean_array(nodes%data(elements%data(i,:),1))
          C(2) = mean_array(nodes%data(elements%data(i,:),2))
          best_segment_distance = huge(1.0_rkind)
          channel: do j=1, channel_el%kolik
                    A = channel_nd%data(channel_el%data(j,1),:)
                    B = channel_nd%data(channel_el%data(j,2),:)
                    if (norm2(B-A) <= 10.0_rkind*epsilon(1.0_rkind)) cycle channel
                    segment_distance = point_segment_distance(A, B, C)
                    if (segment_distance < best_segment_distance) then
                      best_segment_distance = segment_distance
                      ncfluxdata%fluxvct(i,:) = unit_vector(A,B)
                    end if
          end do channel
        end if
      end do

      call routing_apply_directions(success,errmsg)
      if (.not.success) then
        print *, trim(errmsg)
        error stop 'Unable to apply mRM routing directions'
      end if
      call ncflux_prepare_widths(success, errmsg, bccnt)
      if (.not. success) then
        call write_log(trim(errmsg))
        error stop "Unable to compute active river widths"
      end if
      if (bccnt > 0) then
        write(errmsg, *) "W: active river FE elements with no positive width:", bccnt, &
          "; check corner-only/tangential contacts and flow directions."
        call write_log(trim(errmsg))
      end if
      if (LShydro) then
        call hydro_filter()
        call prepare_adenc_banks(original_edge)
        call hydro_initialize(original_edge)
      end if

      allocate(ncelements%areas(ncelements%kolik))
      
      do i=1, ncelements%kolik
        do j=1,4
          pts(j,:) = ncnodes%data(ncelements%data(i,j),1:2)
        end do

        ncelements%areas(i) = quad_area(pts(1,:), pts(2,:), pts(3,:), pts(4,:))
      end do
      
      
      
      call read_adenc_boundaries(fileconf)
      if (LSbank_noflow) then
        ! This ID now belongs ONLY to unused nodes; not the participating banks.
        pde(1)%bc(addedbc)%code=1
        pde(1)%bc(addedbc)%file=.false.
        pde(1)%bc(addedbc)%value=0
      end if
      close(fileconf)


      pde(1)%problem_name(1) = "ADE_in_watershed"
      pde(1)%problem_name(2) = "Advection-dispersion-reaction equation for watershed large scale"

      pde(1)%solution_name(1) = "solute_concentration" !nazev vystupnich souboru
      pde(1)%solution_name(2) = "c  [M/L^3]" !popisek grafu

      pde(1)%flux_name(1) = "conc_flux"  
      pde(1)%flux_name(2) = "concentration flux [M.L^{-2}.T^{-1}]"
      
      allocate(pde(1)%mass_name(1,2))

      pde(1)%mass_name(1,1) = "conc_in_river"
      pde(1)%mass_name(1,2) = "concetration [M/L^3]"
      if (LSconservative) then
        pde(1)%mass_name(1,1) = "solute_inventory_density"
        pde(1)%mass_name(1,2) = "H*c [M/L^2]"
        ! ncflux exports hydrological q, NOT C*q-K*grad(C); retain its filename.
        pde(1)%flux_name(2) = "depth-integrated water flux [L^2/T]"
      end if
      
      pde(1)%print_mass = .true.      

    end subroutine netcdf
    
    
    subroutine read_ncorigin(ncid, ncstart)
      use netcdf
      use datetime
      implicit none

      integer, intent(in) :: ncid
      type(datetime_t), intent(out) :: ncstart

      integer :: ierr, time_varid, p
      character(len=256) :: units
      character(len=256) :: refstr

      ! find variable "time"
      ierr = nf90_inq_varid(ncid, "time", time_varid)
      if (ierr /= nf90_noerr) then
        print *, "Cannot find variable 'time': ", trim(nf90_strerror(ierr))
        stop
      end if

      ! read attribute units, e.g. "hours since 1950-01-01 00:00:00"
      ierr = nf90_get_att(ncid, time_varid, "units", units)
      if (ierr /= nf90_noerr) then
        print *, "Cannot read attribute 'units': ", trim(nf90_strerror(ierr))
        stop
      end if

      units = trim(adjustl(units))

      ! locate the word "since"
      p = index(units, "since")
      if (p <= 0) then
        print *, "Unsupported time units format: ", trim(units)
        stop
      end if

      ! take everything after "since"
      refstr = adjustl(units(p+5:))
      refstr = trim(refstr)

      ! initialize defaults
      ncstart%year   = 0
      ncstart%month  = 0
      ncstart%day    = 0
      ncstart%hour   = 0
      ncstart%minute = 0
      ncstart%second = 0

      ! support formats:
      !   YYYY-MM-DD HH:MM:SS
      !   YYYY-MM-DDTHH:MM:SSZ
      !   YYYY-MM-DD
      if (len_trim(refstr) >= 19) then
        read(refstr(1:4),  *) ncstart%year
        read(refstr(6:7),  *) ncstart%month
        read(refstr(9:10), *) ncstart%day

        if (refstr(11:11) == 'T' .or. refstr(11:11) == ' ') then
          read(refstr(12:13), *) ncstart%hour
          read(refstr(15:16), *) ncstart%minute
          read(refstr(18:19), *) ncstart%second
        end if

      else if (len_trim(refstr) >= 10) then
        read(refstr(1:4),  *) ncstart%year
        read(refstr(6:7),  *) ncstart%month
        read(refstr(9:10), *) ncstart%day
      else
        print *, "Reference date in units is too short: ", trim(refstr)
        stop
      end if

    end subroutine read_ncorigin
    
    
    subroutine read_ncbounds(ncid, ok, errmsg)
      use netcdf
      use typy
      use ncglobvars
      

      integer, intent(in) :: ncid
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg

      integer :: ierr
      integer :: lat_bnds_varid, lon_bnds_varid
      real(kind=rkind), allocatable :: tmp_lat_bnds(:,:)
      real(kind=rkind), allocatable :: tmp_lon_bnds(:,:)

      ok = .false.
      errmsg = "read_ncbounds: unknown error"

      ncfluxdata%has_bounds = .false.

      ierr = nf90_inq_varid(ncid, "lat_bnds", lat_bnds_varid)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot find lat_bnds: " // trim(nf90_strerror(ierr))
        return
      end if

      ierr = nf90_inq_varid(ncid, "lon_bnds", lon_bnds_varid)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot find lon_bnds: " // trim(nf90_strerror(ierr))
        return
      end if

      if (allocated(ncfluxdata%lat_bnds)) deallocate(ncfluxdata%lat_bnds)
      if (allocated(ncfluxdata%lon_bnds)) deallocate(ncfluxdata%lon_bnds)

      ! ncdump shows:
      !   lat_bnds(lat,bnds)
      !   lon_bnds(lon,bnds)
      !
      ! For NetCDF Fortran, read into reversed shape:
      !   tmp_lat_bnds(bnds,lat)
      !   tmp_lon_bnds(bnds,lon)
      allocate(tmp_lat_bnds(2, ncfluxdata%nlat))
      allocate(tmp_lon_bnds(2, ncfluxdata%nlon))

      ierr = nf90_get_var(ncid, lat_bnds_varid, tmp_lat_bnds)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot read lat_bnds: " // trim(nf90_strerror(ierr))
        deallocate(tmp_lat_bnds, tmp_lon_bnds)
        return
      end if

      ierr = nf90_get_var(ncid, lon_bnds_varid, tmp_lon_bnds)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot read lon_bnds: " // trim(nf90_strerror(ierr))
        deallocate(tmp_lat_bnds, tmp_lon_bnds)
        return
      end if

      ! Store internally as bounds(index,1:2)
      allocate(ncfluxdata%lat_bnds(ncfluxdata%nlat, 2))
      allocate(ncfluxdata%lon_bnds(ncfluxdata%nlon, 2))

      ncfluxdata%lat_bnds = transpose(tmp_lat_bnds)
      ncfluxdata%lon_bnds = transpose(tmp_lon_bnds)

      deallocate(tmp_lat_bnds, tmp_lon_bnds)

      ncfluxdata%has_bounds = .true.
      ok = .true.
      errmsg = "everything ok"
      
    end subroutine read_ncbounds
    
    
    subroutine set_adenc_boundary_id()
      use typy
      use globals, only: nodes
      use ncglobvars, only: addedbc
      implicit none

      if (any(nodes%edge > 0_ikind .and. nodes%edge < 101_ikind) .or. maxval(nodes%edge) < 101_ikind) &
        error stop "ADEnc mesh boundary IDs must start at 101"
      addedbc = maxval(nodes%edge) + 1_ikind
    end subroutine set_adenc_boundary_id

    subroutine read_adenc_boundaries(fileconf)
      use typy
      use ncglobvars, only: addedbc
      use pde_objs, only: pde
      use readtools, only: readbcvals
      implicit none
      integer, intent(in) :: fileconf

      call readbcvals(unitW=fileconf, struct=pde(1)%bc, dimen=addedbc-100_ikind, &
        dirname="drutes.conf/netcdf/", highest_boundary_id=addedbc)
    end subroutine read_adenc_boundaries

    subroutine readchannel()
      use typy
      use ncglobvars
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use, intrinsic :: iso_fortran_env, only: iostat_end
      implicit none

      integer :: fileid, ierr, line_number
      integer(kind=ikind) :: channel_index, counter, i
      integer(kind=ikind) :: node_offset, element_offset
      integer(kind=ikind) :: total_nodes, total_elements
      integer(kind=ikind), dimension(:), allocatable :: node_counts
      real(kind=rkind), dimension(2) :: tmp
      real(kind=rkind), dimension(2) :: previous
      logical :: eof
      character(len=256) :: channel_file

      if (channel_count < 1_ikind) error stop "ADEnc requires at least one channel"
      allocate(node_counts(channel_count))
      total_nodes = 0_ikind
      total_elements = 0_ikind

      ! Count each polyline separately. This prevents the final point of one
      ! channel from being connected to the first point of the next channel.
      do channel_index = 1, channel_count
        call channel_filename(channel_index, channel_file)
        open(newunit=fileid, file=trim(channel_file), iostat=ierr, status="old", action="read")

        if (ierr /= 0) then
          print *, "unable to open file ", trim(channel_file)
          ERROR STOP
        end if

        counter = 0_ikind
        line_number = 0
        do
          call read_point(fileid, channel_file, line_number, tmp, eof)
          if (eof) exit
          if (counter > 0_ikind) then
            if (norm2(tmp-previous) <= 10.0_rkind*epsilon(1.0_rkind)) &
              call point_error(channel_file, line_number, "duplicate consecutive channel points")
          end if
          previous = tmp
          counter = counter + 1_ikind
        end do
        close(fileid)

        if (counter < 2_ikind) then
          print *, "file ", trim(channel_file), " doesn't contain enough values; at least two points are required"
          ERROR STOP
        end if

        node_counts(channel_index) = counter
        total_nodes = total_nodes + counter
        total_elements = total_elements + counter - 1_ikind
      end do

      if (allocated(channel_nd%data)) deallocate(channel_nd%data)
      if (allocated(channel_el%data)) deallocate(channel_el%data)
      allocate(channel_nd%data(total_nodes,2))
      allocate(channel_el%data(total_elements,2))
      channel_nd%kolik = total_nodes
      channel_el%kolik = total_elements

      node_offset = 0_ikind
      element_offset = 0_ikind

      do channel_index = 1, channel_count
        call channel_filename(channel_index, channel_file)
        open(newunit=fileid, file=trim(channel_file), iostat=ierr, status="old", action="read")

        if (ierr /= 0) then
          print *, "unable to reopen file ", trim(channel_file)
          ERROR STOP
        end if

        line_number = 0
        do i = 1, node_counts(channel_index)
          call read_point(fileid, channel_file, line_number, channel_nd%data(node_offset+i,:), eof)
          if (eof) call point_error(channel_file, line_number, "channel file changed during loading")
        end do
        close(fileid)

        do i = 1, node_counts(channel_index) - 1_ikind
          channel_el%data(element_offset+i,1) = node_offset + i
          channel_el%data(element_offset+i,2) = node_offset + i + 1_ikind
        end do

        node_offset = node_offset + node_counts(channel_index)
        element_offset = element_offset + node_counts(channel_index) - 1_ikind
      end do

      deallocate(node_counts)

    contains

      subroutine point_error(filename, line, message)
        character(len=*), intent(in) :: filename, message
        integer, intent(in) :: line
        write(*,'(a,a,a,i0,a,a)') "Invalid channel file ", trim(filename), " at line ", line, ": ", message
        error stop "Invalid channel geometry"
      end subroutine point_error

      subroutine read_point(unit, filename, line, point, eof)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename
        integer, intent(inout) :: line
        real(kind=rkind), intent(out) :: point(2)
        logical, intent(out) :: eof
        character(len=4096) :: record
        integer :: status, pos, last, first, count, comment_pos

        eof = .false.
        do
          read(unit,'(a)',iostat=status) record
          if (status == iostat_end) then
            eof = .true.
            return
          end if
          line = line + 1
          if (status /= 0) call point_error(filename, line, "unable to read line")
          if (len_trim(record) == len(record)) call point_error(filename, line, "line is too long")
          comment_pos = index(record, '#')
          if (comment_pos > 0) record(comment_pos:) = ' '
          last = len_trim(record)
          if (last == 0) cycle

          count = 0
          pos = 1
          do while (pos <= last)
            if (record(pos:pos) == ' ' .or. record(pos:pos) == achar(9)) then
              pos = pos + 1
              cycle
            end if
            first = pos
            do while (pos <= last)
              if (record(pos:pos) == ' ' .or. record(pos:pos) == achar(9)) exit
              pos = pos + 1
            end do
            count = count + 1
            if (count > 2) call point_error(filename, line, "expected exactly two coordinates")
            ! Disallow list-directed null/repeated values and slash termination.
            if (scan(record(first:pos-1), ',/*') > 0) &
              call point_error(filename, line, "expected a numeric coordinate")
            read(record(first:pos-1),*,iostat=status) point(count)
            if (status /= 0) call point_error(filename, line, "expected a numeric coordinate")
            if (.not. ieee_is_finite(point(count))) &
              call point_error(filename, line, "coordinates must be finite")
          end do
          if (count /= 2) call point_error(filename, line, "expected exactly two coordinates")
          return
        end do
      end subroutine read_point

      subroutine channel_filename(index, filename)
        integer(kind=ikind), intent(in) :: index
        character(len=*), intent(out) :: filename

        if (index == 1_ikind) then
          filename = "drutes.conf/netcdf/channel.dat"
        else
          write(filename, '("drutes.conf/netcdf/channel",I0,".dat")') index
        end if
      end subroutine channel_filename

    end subroutine readchannel 






end module init_netcdf
