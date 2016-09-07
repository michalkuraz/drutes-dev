! Copyright 2008 Michal Kuraz, Petr Mayer


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

module globals
  use typy
  use global_objs
  use sparsematrix


  !> version id variable
  type(version), public :: version_id   
  !> file units ----------------------------------------------------------------------
  integer, public :: file_global,  file_mesh, file_itcg, file_wwwglob
  !> log file unit
  integer, public :: logfile
  !> debug file unit (put there whatever you like)
  integer, public :: debugfile
  !> terminal ID, the print_level sets its value
  integer, public :: terminal = 6
  logical, public :: terminal_assigned = .false.
  !> array of the gauss quadrature formula roots
  real(kind=rkind), dimension(:,:), allocatable, public :: gauss_integ

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!-----------------basic variables------------------------!!!!!!!!!!!!!!!!
  !> program configuration
  type(configuration), public :: drutes_config
  !> maximal number of iterations of Picard method for nonlinear problems
  integer(kind=ikind), public :: max_itcount
  !> Picard iteration threshold value
  real(kind=rkind), public :: iter_criterion
  !> an intial dt
  real(kind=rkind), public :: init_dt
  !>  simulation end time
  real(kind=rkind), public :: end_time
  !> minimal dt
  real(kind=rkind), public :: dtmin
  !> maximal dt
  real(kind=rkind), public :: dtmax
  !> time units
  character(len=5), public :: time_units
  !> maximal pressure to precalculate (tabelarise) the unsaturated hydraulic conductivity value
  real(kind=rkind), public :: maxpress
  !> number of outer boundaries
  integer(kind=ikind) :: outer_boundaries
  !> output print level
  !! 0 = standard output 
  !! 1 = screen goes to out/screen.log
  !! -1 = does goes to /dev/null, not suitable for m$ hell
  !<
  integer(kind=ikind), public :: print_level
  !> level of informations that are supplied to the user
  !! 0 the default level, only information about the time step solution is displayed
  !! 1 the detailed information about each iteration is displayed
  !! 2 extra detailed, even the information about the iterations of the conjugate gradient solver is displayed
  !<
  integer, public :: info_level
  !> observation time array
  type(observe_time_str), dimension(:), allocatable, public :: observe_time
  
  type(observe_info_str), public :: observe_info


  !> the number of observation points
  !! each row defines each point, the row is labeled in [observation] structure as [xyz]
  !! each row has an integer value, this value defines the element number to which this observation point belongs
  !<
  type(observation),  dimension(:), allocatable, public :: observation_array
  !> the number of points with observation data
  type(measured), dimension(:), allocatable, public :: measured_pts
  !> code defining the example type in format XYZ
  !! 10 - Richards equation in single regime and 1D
  !! 11 - Richards equation in dual regime and 1D
  !! 20 - Richards equation in single regime and 2D
  !! 21 - Richards equation in dual regime and 2D
  !! 100 - Richards equation and transport ADE equation in single regime and 1D
  !! 110 - Richards equation and transport ADE equation in dual regime and 1D
  !! 200 - Richards equation and transport ADE equation in single regime and 2D
  !! 210 - Richards equation and transport ADE equation in dual regime and 2D
  !<
!   integer, public :: problem_type
  !> the time period in hours between each backup
  real(kind=rkind), public :: backup_time
  !> stores the number of performed backups
  integer(kind=ikind), public :: backup_runs
  !> current directory name
  character(len=256), public :: dir_name
  !> CPU start time value
  real, public :: start_time

  !> number of linear algebra solver calls
  integer(kind=ikind), public :: solver_call
  type(node), public :: nodes
  type(element),  public :: elements
  !> local capacity matrix and local stiffness matrix
  real(kind=rkind), dimension(:,:), allocatable, public :: cap_mat, stiff_mat
  !> local matrix vector
  real(kind=rkind), dimension(:), allocatable, public :: bside 
  !> points for gauss quadrature on line \n
  !! component point is node's coordinate \n
  !! component weight is gaussian weight \n
  !! this array is derived for the unit length 1 and interval <0,1>, thus it differs from the init_gauss procedure, see init_gauss documentation
  !<
  type(integnodes),  public :: gauss_points
  !> code defining gauss quadrature integration method
      !!10 - 1 point formula, 
      !!20 - 2 point formula
      !!30 - 3 point formula
      !!40 - 4 point formula
      !!50 - 5 point formula
      !!60 - 6 point formula
      !!70 - 7 point formula
      !!90 - 9 point formula (in 2-dimensions derived from 3 point formula, ask Michal Kuraz :) )
      !!120 - 12 point formula
  !<
  integer(kind=ikind), public :: integ_method


  !> defines directory with global configuration files
  type(dirglob_str), public :: dirglob
  
  !> values of base function on gauss quadrature nodes
  real(kind=rkind), dimension(:,:), allocatable, public :: base_fnc
  !> the total time value
  real(kind=rkind), public :: time
  !> time step value
  real(kind=rkind), public :: time_step
  !> equals to the actual time level plus time step which is being processed
  real(kind=rkind), public :: time4solve
  !> previous time step
  real(kind=rkind), public :: dtprev

  !>minimal adjusted time step
  real(kind=rkind), public :: minimal_dt




  !> sparse matrix
  type(extsmtx), public, save :: spmatrix
  !> sparse matrix
  type(extsmtx), public, save :: spmatrix2

  !> current position in the global stiffness matrix
  integer(kind=ikind), public :: global_row_id

  
  logical, public :: debug_print=.false.
  
  logical, public :: debugmode

  logical, public :: coupled_problem


  !> the amount of postprocess runs
  integer(kind=ikind), public :: postpro_run
  !> the amount of decimal numbers in postproces runs
  integer(kind=ikind), public :: postpro_dec
 
  character(len=1024), public :: backup_file

  real(kind=rkind), dimension(:,:), allocatable :: integnode_now, integnode_prev
  real(kind=rkind), dimension(:), allocatable :: elnode_prev
  logical, public :: www = .false.
  real(kind=rkind), public :: cpu_max_time
  logical, public :: cpu_time_limit = .false.
  
  !cummulative iteration count
  integer(kind=ikind), public :: itcum

  



end module globals