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



!> \file stiffmat.f90
!! \briefLocal conductivity (stiffness) matrix assembler.
!<






!> stiffness matrix procedure for non-linear problems linearized by Picard method
module stiffmat
  public :: build_stiff_np
  public :: build_bvect
  contains
   
   
   !> build local stifness matrix for nonlinear problems and Picard method
     subroutine build_stiff_np(el_id, dt, domain_id, quadpnt_in)
      use typy
      use globals
      use global_objs
      use pde_objs
      use feminittools
      use linAlg
      use debug_tools
      use geom_tools

      integer(kind=ikind), intent(in) :: el_id
      !> time step
      real(kind=rkind), intent(in) :: dt
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str), intent(in out), optional :: quadpnt_in

      
      integer(kind=ikind), dimension(:,:), allocatable, save :: layer
      !> sum of dispersive terms
      real(kind=rkind), dimension(1,1) :: dsum
      !> sum of convection terms
      real(kind=rkind) :: csum
      !> sum of reaction terms
      real(kind=rkind) :: rsum
      integer(kind=ikind) :: i,j,l, top, iproc, jproc, ii, jj, limits
      real(kind=rkind), dimension(1,3) :: u
      real(kind=rkind), dimension(3,1) :: v
      real(kind=rkind), dimension(1,3) :: w
      real(kind=rkind), dimension(3) :: conv
      real(kind=rkind), dimension(3,3) :: disp
      type(integpnt_str) :: quadpnt
      real(kind=rkind) :: capacity, react

      stiff_mat = 0

      if (.not. allocated(layer)) then
        allocate(layer(ubound(pde,1), ubound(pde,1)))
      end if

      do iproc=1,ubound(pde,1)
        do jproc=1,ubound(pde,1)
            layer(iproc, jproc) = elements%material(el_id)
        end do
      end do
      
      limits = ubound(stiff_mat,1)/ubound(pde,1)

      top = drutes_config%dimen
      
      if (present(quadpnt_in)) then
        quadpnt=quadpnt_in
      end if

      quadpnt%element = el_id
      quadpnt%column = 2
      quadpnt%type_pnt = "gqnd"
      
      if (present(domain_id)) then
        quadpnt%ddlocal = .true.
        quadpnt%subdom = domain_id
      end if
      
      do iproc=1,ubound(pde,1)
        !use Galerkin FEM
        if (pde(iproc)%diffusion) then
          do jproc=1, ubound(pde,1)
            pde_block_column = jproc
            do i=1, ubound(stiff_mat,1)/ ubound(pde,1)
              do j=1, ubound(stiff_mat,1)/ubound(pde,1)
                dsum = 0
                csum = 0
                rsum = 0
                
                
                v(1:top,1) = elements%ders(el_id,i,1:top)
                u(1,1:top) = elements%ders(el_id,j,1:top)
                

                do l=1, ubound(gauss_points%weight,1)
                  quadpnt%order = l
                  call pde(iproc)%pde_fnc(jproc)%dispersion(pde(iproc), layer(iproc, jproc), &
                     quadpnt, tensor=disp(1:top,1:top))
                     

                  w(:,1:top) =  matmul(u(:,1:top),disp(1:top,1:top))
                  dsum = dsum - matmul(w(:,1:top) ,v(1:top,:))*gauss_points%weight(l)
                  
                  
                end do
                do l=1, ubound(gauss_points%weight,1)
                  quadpnt%order = l
                  call pde(iproc)%pde_fnc(jproc)%convection(pde(iproc), layer(iproc, jproc), quadpnt, &
                    vector_out=conv(1:top))
                  csum = csum - dot_product(u(1,1:top),conv(1:top))*base_fnc(i,l)*gauss_points%weight(l)
                  call pde(iproc)%pde_fnc(jproc)%der_convect(pde(iproc), layer(iproc, jproc), quadpnt, 	&
                     vector_out=conv(1:top))
                    w = base_fnc(i,l)*base_fnc(j,l)
                  csum = csum - dot_product(w(1,1:top), conv(1:top))*gauss_points%weight(l)

                end do


                do l=1, ubound(gauss_points%weight,1)
                  quadpnt%order = l
                  rsum = rsum + pde(iproc)%pde_fnc(jproc)%reaction(pde(iproc),layer(iproc, jproc), &
                    quadpnt)*base_fnc(i,l)*base_fnc(j,l)*gauss_points%weight(l)
                end do	      

                ii = i + (iproc-1)*limits
                jj = j + (jproc-1)*limits

                stiff_mat(ii,jj) = (dsum(1,1) + csum + rsum)*dt
 
              end do
            end do 
          end do
        
        !use least-square FEM (for convection only problems)
        else
           do i=1, ubound(stiff_mat,1)/ ubound(pde,1)
            do j=1, ubound(stiff_mat,1)/ubound(pde,1)
                csum = 0
                v(1:top,1) = elements%ders(el_id,i,1:top)
                u(1,1:top) = elements%ders(el_id,j,1:top)

                do l=1, ubound(gauss_points%weight,1)
                  quadpnt%order = l
                  call pde(iproc)%pde_fnc(iproc)%convection(pde(iproc), layer(iproc, iproc), quadpnt, &
                    vector_out=conv(1:top))
                  csum = csum + ((base_fnc(i,l) + dt*dot_product(conv(1:top), v(1:top,1)))* &
                          (base_fnc(j,l) + dt*dot_product(conv(1:top), u(1,1:top))))* &
                          gauss_points%weight(l)
                 end do
                 ii = i + (iproc-1)*limits
                 jj = j + (iproc-1)*limits
                
                stiff_mat(ii,jj) =  csum 
               end do
             end do
          
        end if
      end do

     
     stiff_mat = stiff_mat/gauss_points%area*elements%areas(el_id)
     

     

     
    end subroutine build_stiff_np
    
  
    
    
    subroutine build_bvect(el_id, dt, domain_id, quadpnt_in)    
      use typy
      use globals
      use global_objs
      use pde_objs
      use feminittools
      use linAlg
      use debug_tools
      
      integer(kind=ikind), intent(in) :: el_id
      real(kind=rkind), intent(in) :: dt
      !>subdomain number (inserted only if domain decomposition used and if only local data needed)
      integer(kind=ikind), intent(in), optional :: domain_id
      type(integpnt_str) , intent(in), optional :: quadpnt_in
      
      integer(kind=ikind) :: iproc, limits, ii, i, l, top, layer
      real(kind=rkind) :: suma, hp, source, cl, cs
      type(integpnt_str) :: quadpnt
      real(kind=rkind), dimension(3) :: conv
      real(kind=rkind), dimension(3,1) :: v
      real(kind=rkind) :: capacity, react
      
      bside = 0
      
      top = drutes_config%dimen
      
      if (present(quadpnt_in)) then
        quadpnt = quadpnt_in
      end if
      
      quadpnt%element = el_id
      quadpnt%column = 2
      quadpnt%type_pnt = "gqnd"
      
      limits = ubound(stiff_mat,1)/ubound(pde,1)

      if (present(domain_id)) then
        quadpnt%ddlocal = .true.
        quadpnt%subdom = domain_id
      end if
      
      do iproc = 1, ubound(pde,1)
        if (pde(iproc)%diffusion) then
          do i=1, ubound(stiff_mat,1)/ ubound(pde,1)
            suma = 0
            
            do l=1, ubound(gauss_points%weight,1)
              quadpnt%order = l	  
              suma = suma - pde(iproc)%pde_fnc(iproc)%zerord(pde(iproc), layer=elements%material(el_id), quadpnt=quadpnt)* &
              gauss_points%weight(l)
            end do
      
            ii = i + (iproc-1)*limits
            
            bside(ii) = suma*dt*elements%areas(el_id)/gauss_points%area

          end do
        else
        
          layer = elements%material(el_id)
          
          do i=1, ubound(stiff_mat,1)/ ubound(pde,1)
            v(1:top,1) = elements%ders(el_id,i,1:top)
            suma = 0

            do l=1,  ubound(gauss_points%weight,1)
              quadpnt%order = l
              quadpnt%column = 1
              hp = pde(iproc)%getval(quadpnt)
              quadpnt%column = 2
              source = pde(iproc)%pde_fnc(iproc)%zerord(pde(iproc), layer, quadpnt=quadpnt)
              call pde(iproc)%pde_fnc(iproc)%convection(pde(iproc), layer, quadpnt, &
                      vector_out=conv(1:top))
              suma = suma + (hp + dt*source) *(base_fnc(1,l)+dot_product(dt*conv(1:top), v(1:top,1)))* &
                    gauss_points%weight(l)  
            end do
            
            ii = i + (iproc-1)*limits
            bside(ii) = suma*elements%areas(el_id)/gauss_points%area
          end do
        end if
      end do
    
      

    end subroutine build_bvect


end module stiffmat
