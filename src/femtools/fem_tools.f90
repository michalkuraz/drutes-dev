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


!> \file fem_tools.f90
!! \brief FEM tools
!<

!> Main FEM matrix assembling tools. Local matrix to global matrix procedures and mass checks.

module fem_tools
  public :: in2global
  public :: do_masscheck
  private :: get_bcflux, get_mass


  contains

  !> procedure to put local stiffness matrix into the global stiffness matrix
  !<
  subroutine in2global(el_id, locmatrix, bvect, subprmt, debugme)
    use typy
    use sparsematrix
    use global_objs
    use globals
    use pde_objs
    use globals
    use debug_tools
    use geom_tools
    
    !> number of element
    integer(kind=ikind), intent(in) :: el_id
    class(extsmtx), intent(in out)  :: locmatrix
    real(kind=rkind), dimension(:) :: bvect
    !> permutation array for submatrices, supplied only for domain decomposition problems
    integer(kind=ikind), dimension(:), intent(in), optional :: subprmt
    logical, intent(in), optional :: debugme
    
    integer(kind=ikind), dimension(:), allocatable, save :: bc
    integer(kind=ikind), dimension(:), allocatable, save :: n_row, m_col
    integer(kind=ikind) :: i,j,m, iproc, jproc, limits, fnc, locneu
    real(kind=rkind), dimension(:,:), allocatable, save :: bcval
    real(kind=rkind), dimension(:), allocatable, save :: surface
    integer(kind=ikind), dimension(:), allocatable, save :: fin
    integer(kind=ikind) :: space_dim, space_i, kk, kkk, extpt
    real(kind=rkind), dimension(3,3) :: d
    real(kind=rkind) :: tmp   
    real(kind=rkind), dimension(3) :: tmpdata
    integer(kind=ikind), dimension(:), allocatable, save :: neumann
    type(bcpts_str), save :: bcpts
   

    if (.not. allocated(bc)) then	
      allocate(bc(ubound(stiff_mat,1)))
      allocate(n_row(ubound(stiff_mat,1)))
      allocate(m_col(ubound(stiff_mat,1)))
      allocate(bcval(ubound(stiff_mat,1),3))
      allocate(surface(ubound(stiff_mat,1)*ubound(pde,1)))
      allocate(fin(ubound(pde,1)))
      allocate(neumann(drutes_config%dimen+1))
      allocate(bcpts%border(drutes_config%dimen))
    end if
    
    surface = 0
    bc = 0
    n_row = 0
    m_col = 0
    bcval = 0
    limits = ubound(stiff_mat,1)/ubound(pde,1)
    
        

    do iproc = 0,ubound(pde,1)-1
      do i=1, ubound(elements%data, 2)
        j = elements%data(el_id,i)
        n_row(i+iproc*limits) = pde(iproc+1)%permut(j)
        if (nodes%edge(j) > 100) then
          m = nodes%edge(j)
          call pde(iproc+1)%bc(m)%value_fnc(pde(iproc+1), el_id,i,code=bc(i+iproc*limits))
        else
          bc(i) = 0
        end if
      end do
    end do

    select case(drutes_config%dimen)
      case(1)
        surface = 1
         do iproc=0, ubound(pde,1)-1
          do i=1, limits
            if (bc(i+iproc*limits) == 2) then
              neumann = i
              EXIT
            end if
          end do
         end do  
         if (neumann(1) == 2) then
           extpt = 1
         else
           extpt = 2
         end if
      case(2)

        do iproc=0, ubound(pde,1)-1
          locneu = 1
          do i=1, limits
            if (bc(i+iproc*limits) == 2) then
              neumann(locneu) = i
              locneu = locneu + 1
            end if
          end do
          if (locneu == 3) then
            do i=1, limits
              if (bc(i+iproc*limits) == 2) then
                surface(i+iproc*limits) = dist(nodes%data(elements%data(el_id,neumann(1)),:), &
                                              nodes%data(elements%data(el_id,neumann(2)),:))/drutes_config%dimen
                 select case(neumann(1) + neumann(2))
                 	!1 + 2
                 	 case(3)
                      extpt = 3
                    !1 + 3
                 	 case(4)
                      extpt = 2
                   case(5)
                      extpt = 1
                 end select


              end if
            end do
          end if
        end do
      case(3)
        do iproc=0, ubound(pde,1)-1
          locneu = 1
          do i=1, limits
            if (bc(i+iproc*limits) == 2) then
              neumann(locneu) = i
              locneu = locneu + 1
            end if
          end do
          if (locneu == 4) then
            do i=1, limits
              if (bc(i+iproc*limits) == 2) then
                surface(i+iproc*limits) = triarea(nodes%data(elements%data(el_id,neumann(1)),:), &
                                                  nodes%data(elements%data(el_id,neumann(2)),:),  &
                                                  nodes%data(elements%data(el_id,neumann(3)),:))/drutes_config%dimen
                                              
                select case(neumann(1) + neumann(2) + neumann(3))
                  ! 1 + 2 + 3
                  case(6)
                  	extpt = 4
                  ! 1 + 3 + 4	
                  case(8)
                    extpt = 2
                  ! 1 + 2 + 4
                  case(7)
                    extpt = 3
                  !2 + 3 + 4
                  case(9)
                  	extpt = 1
                end select     
              end if
            end do
          end if
        end do
      

    end select


   bcpts%border = neumann
   bcpts%extpt = extpt
  



    do iproc = 0,ubound(pde,1)-1
      do i=1, limits
        select case(bc(i+iproc*limits))
          case(0)
            CONTINUE
          case(3)
            bcval(i+iproc*limits,1) = 0
          case(1,4)
            m = nodes%edge(elements%data(el_id, i))
            call pde(iproc+1)%bc(m)%value_fnc(pde(iproc+1), element=el_id,node=i,value=bcval(i+iproc*limits,1))
          case(2)
          	m = nodes%edge(elements%data(el_id, i))
            if (abs(surface(i + iproc*limits)) > 10*epsilon(1.0_rkind)) then
              call pde(iproc+1)%bc(m)%value_fnc(pde(iproc+1), element=el_id,node=i,value=bcval(i+iproc*limits,1), &
                                                bcpts=bcpts)
            else
              bcval(i+iproc*limits,1) = 0
            end if
          case(5)
            call pde(iproc+1)%bc(m)%value_fnc(pde(iproc+1), el_id,i,valarray=bcval(i+iproc*limits,:))
        end select
      end do
    end do

    

    if (present(subprmt)) then
      do i=1, ubound(m_col,1)
        if (n_row(i) > 0 ) then
          m_col(i) = subprmt(n_row(i))
        end if
      end do
    else
      m_col=n_row
    end if

 

    do i=1, ubound(stiff_mat,1)
      ! fill bside
      if (bc(i) /= 1) then
      select case(bc(i))
        case(0,3)
          if (m_col(i) > 0) then
            bvect(m_col(i)) = bvect(m_col(i)) + bside(i)
          end if
        case(2)
          if (m_col(i) > 0) then
            bvect(m_col(i)) = bvect(m_col(i)) - bcval(i,1)*surface(i)*time_step + bside(i)
          end if
        case(4)
          bvect(m_col(i)) = bcval(i,1)
        case(5)
          bvect(m_col(i)) = bcval(i,3)*surface(i)
          print *, bcval(i,3)
      end select
      
      ! fill stiffness matrix
      do m=1, ubound(stiff_mat,1)
        select case(bc(m))
          case(1)
            if (m_col(i) > 0) then
              bvect(m_col(i)) = bvect(m_col(i)) - stiff_mat(i,m)*bcval(m,1)
            end if
                    
          case(4)
            if (i == m) then
              call locmatrix%set(1.0_rkind, n_row(i), m_col(i))
              if (drutes_config%it_method == 2 .or. drutes_config%it_method == 1) then
                call locmatrix%rowsfilled%nrfill(n_row(i))
              end if
            else
              call locmatrix%add(stiff_mat(i,m), n_row(i), m_col(m))
            end if
          case(5)
            tmp = 0
            space_dim = ubound(stiff_mat,1)/ubound(pde,1)
            if (m > space_dim) then
              space_i = m - space_dim
            else
              space_i = m
            end if
            
            do fnc=1, drutes_config%dimen
              tmp = tmp + elements%ders(el_id, space_i, fnc) 
            end do
            
            if (m /= i) then
              tmp = tmp*bcval(m,1)*surface(i)
              tmp = 0
            else
              tmp = tmp*bcval(m,1)*surface(i) + bcval(m,2)*surface(i)
              tmp = 1
            end if
            
            call locmatrix%add(tmp, n_row(i), m_col(m))
            
        case default
                    
          if (bc(i) /= 4 ) then 
            if (n_row(i) > 0 .and. m_col(m) > 0) then
              call locmatrix%add(stiff_mat(i,m), n_row(i), m_col(m))
              if (drutes_config%it_method == 2 .or. drutes_config%it_method == 1) then
                call locmatrix%rowsfilled%nrfill(n_row(i))
              end if
            end if
          end if
                    
        end select
      end do
      else
        CONTINUE
      end if
    end do
	  

  
  end subroutine in2global
  
  
  !> Get difference between the volume change and boundary fluxes.
  subroutine do_masscheck()
    use typy
    use global_objs
    use globals
    use core_tools
    use pde_objs
    use printtools
    use debug_tools
    
    character(len=512) :: namef
    integer(kind=ikind) :: i, proc
    integer, dimension(:), allocatable, save :: fileids
    integer :: i_err
    real(kind=rkind) :: value
    real(kind=rkind), save :: totvolume
    real(kind=rkind), dimension(:), allocatable, save :: massdiff, totflux, totdiff
    
    
    if (.not. allocated(fileids)) then
    
      allocate(fileids(ubound(pde,1)))
      
      allocate(massdiff(ubound(pde,1)))
      
      allocate(totflux(ubound(pde,1)))
      
      allocate(totdiff(ubound(pde,1)))
      
      totdiff = 0
      
      do proc=1, ubound(pde,1)
        write(namef, *) "out/mass/", cut(pde(proc)%problem_name(1)), "-mass_check.dat"
      
        open(newunit=fileids(proc), file=cut(namef), action="write", status="replace", iostat=i_err)
        
        totvolume=0
        do i=1, elements%kolik
          totvolume = totvolume + elements%areas(i)
        end do
        
        if (i_err /= 0) then
          i_err=system("mkdir out/mass")
          
          if (i_err /= 0) then
            print *, "Unable to create directory out/mass, maybe disc full??"
            print *, "or contact Michal -> michalkuraz@gmail.com"
          end if
          
          open(newunit=fileids(proc), file=cut(namef), action="write", status="replace", iostat=i_err)
          
          if (i_err /= 0) then
            print *, "Unable to create new files in directory out/mass, maybe disc full??"
            print *, "or contact Michal -> michalkuraz@gmail.com"
          end if
          
        end if
        
        call print_logo(fileids(proc))
        
        write(unit=fileids(proc), fmt=*) " #      time                   time step              integral boundary flux [L3] & 
        volume change [L3]      difference [L3]            cum. difference [L3]        cum. difference [%]"
        
        write(unit=fileids(proc), fmt=*) " #--------------------------------------------------------------------------------&
        ---------------------------------------------------------------------------------------------- "
        
      end do
    end if
      
      
    do proc = 1, ubound(pde,1)
      massdiff(proc) =  get_mass(3_ikind, pde(proc)) -  get_mass(1_ikind, pde(proc))
      totflux(proc) = get_bcflux(pde(proc))
      totdiff(proc) = totdiff(proc) + massdiff(proc) - totflux(proc)
      write(unit=fileids(proc), fmt=*) time, time_step, totflux(proc), massdiff(proc), massdiff(proc) - totflux(proc), &
       totdiff(proc), totdiff/totvolume*100.0
    end do
    
    
      
    
   

  
  end subroutine do_masscheck
  
  !> Get boundary fluxes.
  function get_bcflux(pde_loc) result(val)
    use typy
    use pde_objs
    use global_objs
    use globals
    use geom_tools
    use integral
    use debug_tools
    
    
    class(pde_str), intent(in) :: pde_loc
    real(kind=rkind) :: val
    
    
    integer(kind=ikind) :: i, no_points
    real(kind=rkind), dimension(:), allocatable, save :: uzly, points
    real(kind=rkind), dimension(:), allocatable, save :: vahy, weights
    real(kind=rkind), dimension(2,2) :: bcpoints
    real(kind=rkind), dimension(2) :: thirdpt
    real(kind=rkind), dimension(:), allocatable, save :: nvect, flux
    integer(kind=ikind), dimension(3) :: ellocs
    integer(kind=ikind) :: el, j, k, pt
    type(integpnt_str) :: quadpnt_loc
    real(kind=rkind) :: locval
    
    
    if (.not. allocated(flux)) allocate(flux(drutes_config%dimen))
    
    quadpnt_loc%column=3
    
    if (drutes_config%dimen == 1) then
    
      quadpnt_loc%type_pnt = "ndpt"
      quadpnt_loc%order = 1
      
      call pde_loc%flux(elements%material(1), quadpnt_loc, vector_out=flux)
      
    
      
      locval = flux(1)
      quadpnt_loc%order = nodes%kolik
      
      call pde_loc%flux(elements%material(elements%kolik), quadpnt_loc, vector_out=flux)
      
      
      
      val = locval + flux(1)
      
    
    else 
    
    
      quadpnt_loc%type_pnt = "xypt"
      
      
      if (.not. allocated(uzly)) then
        no_points = int(integ_method/10_ikind)
        
        allocate(points(no_points))

        allocate(weights(no_points))
      
        
        if (modulo(no_points,2_ikind) == 0) then
          allocate(uzly(no_points/2))
          allocate(vahy(no_points/2))
        else
          allocate(uzly(no_points/2+1))
          allocate(vahy(no_points/2+1))
        end if

        call getform(no_points, uzly, vahy)


        if (modulo(no_points,2_ikind) == 0) then

          do i = 1, ubound(uzly,1) 
            points(i) = 0.5_rkind - uzly(i)/2.0_rkind
            points(i+ubound(uzly,1)) = 0.5_rkind + uzly(i)/2.0_rkind
            weights(i) = vahy(i)/2.0_rkind
            weights(i+ubound(uzly,1)) = vahy(i)/2.0_rkind
          end do

        else

          points(ubound(uzly,1)) = 0.5_rkind
          weights(ubound(uzly,1)) = vahy(ubound(vahy,1))/2.0_rkind

          do i = 1, ubound(uzly,1) - 1 
            points(i) = 0.5_rkind - uzly(i)/2.0_rkind
            points(i+ubound(uzly,1)) = 0.5_rkind + uzly(i)/2.0_rkind
            weights(i) = vahy(i)/2.0_rkind
            weights(i+ubound(uzly,1)) = vahy(i)/2.0_rkind
          end do

        end if
      end if
        
      val=0
      do i=1, elements%bcel%pos
        el = elements%bcel%data(i)
        quadpnt_loc%element=el
        do j=2,elements%border(el)%pos
          bcpoints(1,:) = nodes%data(elements%border(el)%data(j-1),:)
          bcpoints(2,:) = nodes%data(elements%border(el)%data(j),:)
          
          ellocs=1
          do k=1, ubound(elements%data,2)
            if (elements%border(el)%data(j-1) == elements%data(el,k) .or. elements%border(el)%data(j) == elements%data(el,k)) then
              ellocs(k) = 0
            end if
          end do
          
          thirdpt =  nodes%data(elements%data(el, maxloc(ellocs,1)),:)
          
          !returns outer normal vector
          call getnormal(bcpoints, thirdpt, nvect)
          
          nvect = -nvect
          
          locval = 0
          do pt=1, ubound(points,1)
            quadpnt_loc%xy(1) = bcpoints(1,1) + (bcpoints(2,1)-bcpoints(1,1))*points(pt)
            quadpnt_loc%xy(2) = bcpoints(1,2) + (bcpoints(2,2)-bcpoints(1,2))*points(pt)
            
            call pde_loc%flux(elements%material(el), quadpnt_loc, vector_out=flux)
            
            flux(1) = flux(1)*nvect(1)
            flux(2) = flux(2)*nvect(2)
            
!            print *, el
!            print *, flux
!            print *, "---"
!            call wait()
            
            locval = sqrt(flux(1)*flux(1) + flux(2)*flux(2))*weights(pt)
          end do
          
          locval = locval*dist(bcpoints(1,:), bcpoints(2,:))
          

        end do
        
        val = val + locval
        
      end do
    end if 
    
    
    
     val = val*time_step
    
    
  end function get_bcflux


  !> Get mass on domain.
  function get_mass(column, pde_loc) result(mass)
    use typy
    use pde_objs
    use global_objs
    use globals
    use geom_tools
    use debug_tools
    
    integer(kind=ikind), intent(in) :: column
    class(pde_str), intent(in) :: pde_loc
    real(kind=rkind) :: mass
    
    type(integpnt_str) :: quadpnt_loc
    real(kind=rkind) :: locmass
    integer(kind=ikind) :: i, pt
    
    quadpnt_loc%column=column
    quadpnt_loc%type_pnt="gqnd"
    
    mass = 0.0
    do i=1, elements%kolik
      quadpnt_loc%element = i
      locmass = 0.0
      do pt=1, ubound(gauss_points%weight,1)
       quadpnt_loc%order = pt
       locmass = locmass + pde_loc%mass(1)%val(pde_loc, elements%material(i), quadpnt_loc)*gauss_points%weight(pt)
      end do
      locmass = locmass/gauss_points%area*elements%areas(i)
      mass = mass + locmass
    end do
    
  end function get_mass
  


end module fem_tools
