module Re_dual_reader
  public :: Re_dual_readm,Re_dual_readf, Re_dual_var
  
  contains
    subroutine Re_dual_var()
     use typy
     use globals
     use global_objs
     use core_tools
     use dual_globals
     use readtools
     use pde_objs
     
      integer :: i_err, i, j, file_dual
      integer(kind=ikind) :: layers
	  real(kind=rkind)::inicond
	  	
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !input pars
      
      ! all materials in dual.conf
      call find_unit(file_dual, 300)
      open(unit=file_dual, file="drutes.conf/REdual/dual.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
        print *, "missing dual.conf file in drutes.conf/REdual"
        ERROR STOP
      end if 
      !number of materials
      call fileread(layers,file_dual,ranges=(/0_ikind,huge(0_ikind)/)&
      , errmsg="At least one layer")
      

      ! now that we know the number of layers/materials, let's allocate van genuchten matrices 
      ! vgmatrix, vgexchange and vgfracture all come from dual_globals, exchange contains specific exchange parameters
      
      if (.not. allocated(vgmatrix)) then
	    allocate (vgmatrix(layers))
	   	  do i=1, layers
	  	    allocate(vgmatrix(i)%Ks_local(drutes_config%dimen))
	  		allocate(vgmatrix(i)%Ks(drutes_config%dimen, drutes_config%dimen))	  
	  		if(drutes_config%dimen>1) then
	  		  j = max(1,drutes_config%dimen-1)
	  		  allocate(vgmatrix(i)%anisoangle(j))	
	  		else
	  		  allocate(vgmatrix(i)%anisoangle(1))
	  		end if		
		  end do
      end if
      
      if (.not. allocated(vgexchange)) then
	    allocate (vgexchange(layers))
	   	  do i=1, layers
	  	    allocate(vgexchange(i)%Ks_local(drutes_config%dimen))
	  		allocate(vgexchange(i)%Ks(drutes_config%dimen, drutes_config%dimen))	  
	  		if(drutes_config%dimen>1) then
	  		  j = max(1,drutes_config%dimen-1)
	  		  allocate(vgexchange(i)%anisoangle(j))	
	  		else
	  		  allocate(vgexchange(i)%anisoangle(1))
	  		end if			
		  end do
      end if
      
      if (.not. allocated(exchange)) then
	    allocate(exchange(layers))
      end if
      
      if (.not. allocated(vgfracture)) then
	    allocate (vgfracture(layers))
	   	  do i=1, layers
	  	    allocate(vgfracture(i)%Ks_local(drutes_config%dimen))
	  		allocate(vgfracture(i)%Ks(drutes_config%dimen, drutes_config%dimen))	  
	  		if(drutes_config%dimen>1) then
	  		  j = max(1,drutes_config%dimen-1)
	  		  allocate(vgfracture(i)%anisoangle(j))	
	  		else
	  		  allocate(vgfracture(i)%anisoangle(1))
	  		end if			
		  end do
      end if
      
     
      ! Matrix
      do i=1,layers
        call comment(file_dual)
        read(unit=file_dual, fmt= *, iostat=i_err) vgmatrix(i)%alpha, vgmatrix(i)%n, vgmatrix(i)%m, &
                            vgmatrix(i)%Ths, vgmatrix(i)%Thr, vgmatrix(i)%Ss
        if (i_err /= 0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "HINT: You don't seem to have defined a sufficient number of van Genuchten &
           parameters for all >> matrix << layers!"
          print *, "----------------------------------------"
          call file_error(file_dual)
        end if
      end do
      
    
      do i=1,layers
        call comment(file_dual)
        read(unit=file_dual, fmt= *, iostat=i_err) vgmatrix(i)%anisoangle, vgmatrix(i)%Ks_local(:)
        if(drutes_config%dimen>1) then
	  	  call set_tensor(vgmatrix(i)%Ks_local(:), vgmatrix(i)%anisoangle(:),  vgmatrix(i)%Ks)
	  	else
	  	  vgmatrix(i)%Ks(1,1)=vgmatrix(i)%Ks_local(1)
	  	end if	  
      end do


      ! Fracture
      do i=1,layers
        call comment(file_dual)
        read(unit=file_dual, fmt= *, iostat=i_err) vgfracture(i)%alpha, vgfracture(i)%n, vgfracture(i)%m, &
                            vgfracture(i)%Ths, vgfracture(i)%Thr, vgfracture(i)%Ss
                            
        if (i_err /= 0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "HINT: You don't seem to have defined a sufficient number of van Genuchten &
           parameters for all >> fracture << layers!"
          print *, "----------------------------------------"
          call file_error(file_dual)
        end if
      end do
      
      do i=1,layers
        call comment(file_dual)
        read(unit=file_dual, fmt= *, iostat=i_err) vgfracture(i)%anisoangle, vgfracture(i)%Ks_local(:)
        if(drutes_config%dimen>1) then
	  	  call set_tensor(vgfracture(i)%Ks_local(:), vgfracture(i)%anisoangle(:),  vgfracture(i)%Ks)
	  	else
	  	  vgfracture(i)%Ks(1,1)=vgfracture(i)%Ks_local(1)
	  	end if	  
      end do
      
      ! exchange
      do i=1,layers
        call comment(file_dual)
        read(unit=file_dual, fmt= *, iostat=i_err) vgexchange(i)%alpha, vgexchange(i)%n, vgexchange(i)%m, &
                            vgexchange(i)%Ths, vgexchange(i)%Thr, vgexchange(i)%Ss
        if (i_err /= 0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "HINT: You don't seem to have defined a sufficient number of van Genuchten &
           parameters for all >> exchange << layers!"
          print *, "----------------------------------------"
          call file_error(file_dual)
        end if
      end do
      
      do i=1,layers
      call comment(file_dual)
      read(unit=file_dual, fmt= *, iostat=i_err) vgexchange(i)%anisoangle, vgexchange(i)%Ks_local(:)
        if(drutes_config%dimen>1) then
	  	  call set_tensor(vgexchange(i)%Ks_local(:), vgexchange(i)%anisoangle(:),  vgexchange(i)%Ks)
	  	else
	  	  vgexchange(i)%Ks(1,1)=vgexchange(i)%Ks_local(1)
	  	end if	 
      end do
      
    
	  do i=1,layers
	    call comment(file_dual)
	    read(unit=file_dual, fmt= *, iostat=i_err) exchange(i)%beta, exchange(i)%a,&
	    exchange(i)%gam_par,exchange(i)%weightf,exchange(i)%weightm
        if (i_err /= 0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "HINT: You don't seem to have defined a sufficient number of Exchange &
           parameters for all >> exchange << layers!"
          print *, "----------------------------------------"
          call file_error(file_dual)
        end if
        if((exchange(i)%weightf+exchange(i)%weightm)/=1.0_rkind) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "In dual permeability fracture weight and matrix weight have to sum up to one!"
          print *, "----------------------------------------"
          call file_error(file_dual)
        end if
      end do
      
            
      ! initial conditions
      
      do i=1,layers
		call fileread(inicond, file_dual,ranges=(/-huge(0.0_rkind),0.0_rkind/),&
      	errmsg="Start unsaturated, please.")
        vgmatrix(i)%initcond=inicond
		vgfracture(i)%initcond=inicond
      end do

      close(file_dual)

    end subroutine Re_dual_var
    
    subroutine Re_dual_readm(pde_loc)
      use typy
      use globals
      use global_objs
      use core_tools
      use dual_globals
      use readtools
      use pde_objs
      class(pde_str), intent(in out) :: pde_loc
      integer :: i_err, file_bc
      integer(kind=ikind) :: dim_bc
      pde_loc%problem_name(1) = "Re_dual_totH_m"
      pde_loc%problem_name(2) = "Richards' equation"

      pde_loc%solution_name(1) = "total_head_m" !file output name
      pde_loc%solution_name(2) = "h  [L]" ! caption

      pde_loc%flux_name(1) = "flux_m"  
      pde_loc%flux_name(2) = "Darcian flow [L.T^{-1}]"

      pde_loc%mass_name(1) = "theta_m"
      pde_loc%mass_name(2) = "theta [-]"

    ! boundary values
    
      call find_unit(file_bc, 400)
      open(unit = file_bc,file = "drutes.conf/REdual/dual_bc.conf", action = "read", status = "old", iostat = i_err)
      if (i_err > 0) then
        print *, "something is wrong with your input"
      else if(i_err < 0) then
        print*, "end of file reached"
      end if 
      if (i_err /= 0) then
        print *, "missing boundary condition for dual permeability problem"
        ERROR STOP
      end if  
      call fileread(dim_bc, file_bc)
      call readbcvals(unitW=file_bc,struct=pde_loc%bc,dimen=dim_bc,dirname = "drutes.conf/REdual/dual_bc.conf")

      close(file_bc)
      
    end subroutine Re_dual_readm		

subroutine Re_dual_readf(pde_loc)
      use typy
      use globals
      use global_objs
      use core_tools
      use dual_globals
      use readtools
      use pde_objs
      class(pde_str), intent(in out) :: pde_loc
      integer :: i_err, file_bc
      integer(kind=ikind) :: dim_bc
      pde_loc%problem_name(1) = "Re_dual_totH_f"
      pde_loc%problem_name(2) = "Richards' equation"

      pde_loc%solution_name(1) = "total_head_f" !file output name
      pde_loc%solution_name(2) = "h  [L]" ! caption

      pde_loc%flux_name(1) = "flux_f"  
      pde_loc%flux_name(2) = "Darcian flow [L.T^{-1}]"

      pde_loc%mass_name(1) = "theta_f"
      pde_loc%mass_name(2) = "theta [-]"

    ! boundary values
    
      call find_unit(file_bc, 400)
      open(unit = file_bc,file = "drutes.conf/REdual/dual_bc.conf", action = "read", status = "old", iostat = i_err)
      if (i_err > 0) then
        print *, "something is wrong with your input"
      else if(i_err < 0) then
        print*, "end of file reached"
      end if 
      if (i_err /= 0) then
        print *, "missing boundary condition for dual permeability problem"
        ERROR STOP
      end if  
      call fileread(dim_bc, file_bc)
     
      call readbcvals(unitW=file_bc,struct=pde_loc%bc,dimen=dim_bc,dirname = "drutes.conf/dual/dual_bc.conf")

      close(file_bc)
      
      
    end subroutine Re_dual_readf	
        
end module Re_dual_reader