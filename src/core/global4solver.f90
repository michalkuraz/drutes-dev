module global4solver

  character(len=64), public :: solver_name
  logical, public :: record_solver_time
  integer, public :: solver_time_file
    
  logical, public :: solver_error

end module global4solver
