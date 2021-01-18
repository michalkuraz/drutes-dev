
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

module kinpointer

  contains
  
    subroutine kinwaveprocs(number) 
      use typy
      use readtools
      
      integer(kind=ikind), intent(out) :: number
      integer :: fileid, ierr
      logical :: yes
      
      open(newunit=fileid, file="drutes.conf/kinwave/kinwave.conf",action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open drutes.conf/kinwave/kinwave.conf, exiting..."
        error stop
      end if
      
      call fileread(yes, fileid, errmsg="set y/n if you want to couple your model with solute transport")
      
      close(fileid)
      
      if (yes) then
        number = 3
      else
        number = 1
      end if
      
    end subroutine kinwaveprocs
    
    subroutine kinwavelinker()
      use typy
      use kinfnc
      use kinreader
      use pde_objs
      use debug_tools
      
      
      call kininit()
      
      pde(1)%pde_fnc(1)%convection => kinconvect
      
      allocate(pde(1)%bc(101:101))
      
      pde(1)%bc(101)%value_fnc => kinbor
      
!       pde_loc%bc(102)%value_fnc => kinbor

      if (drutes_config%dimen == 1) then
        nodes%edge(ubound(nodes%data)) = 101
        
        nodes%edge(1) = 0
      end if
    
      pde(1)%initcond => kinematixinit
      
      pde(1)%pde_fnc(1)%zerord => rainfall
      
      pde(1)%pde_fnc(1)%elasticity => kin_elast
      
!      pde(1)%diffusion = .false.
      
      pde(1)%getval => getval_kinwave
      
      pde(1)%flux => kinflux
      
      pde(1)%symmetric = .true.
      
      if (ubound(pde,1) == 3) then
      
        allocate(pde(2)%bc(101:101))
        
        allocate(pde(3)%bc(101:101))
      
        pde(2)%pde_fnc(2)%convection => kinconvectcl
        
        pde(2)%initcond =>  kinematixinit
        
!        pde(2)%pde_fnc(2)%reaction => kincl_source
        
!        pde(2)%pde_fnc(3)%reaction => kincs_source
        
        pde(2)%pde_fnc(2)%elasticity => kin_clelast
        
        pde(2)%pde_fnc(3)%elasticity => kin_cselast
        
        pde(2)%print_mass = .true.
        
        pde(2)%mass(1)%val => solmass
        
!        pde(2)%diffusion = .false.
        
        pde(2)%flux => kinfluxcl
        
        pde(2)%symmetric = .true.
        
        pde(2)%bc(101)%value_fnc => kinbor
        
        !--- soil
        pde(3)%initcond =>  kinematixinit4cs
        
        pde(3)%pde_fnc(2)%elasticity => kin_clelast
        
        pde(3)%pde_fnc(3)%elasticity => kin_elast
        
        pde(3)%pde_fnc(2)%reaction => kincl_source
        
        pde(3)%pde_fnc(3)%reaction => kincs_source
        
        pde(3)%bc(101)%value_fnc => kinborcs
        
      end if
        
        
        
!      pde(2)
      
    end subroutine kinwavelinker
    
    



end module kinpointer
