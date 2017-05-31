program ssdcsolver
  ! Main program for the solution of the compressible Navier-Stokes
  ! equations discretized using a discontinous collocation method (DCM)

  ! Load modules
  use precision_vars
  use fileio
  use referencevariables
  use mpimod
  use navierstokes
  use timeinteg
  use variables, only: xg, ug, vg, wg, uhat, uold, du, dudt
  use controlvariables
  use physicsindependent
  use errorestimation
  use physics_driver
  use write_solution_file, only : create_solution_dir
  use restart_simulation, only : create_restart_dir
  
  ! Nothing is implicitly defined
  implicit none

  ! Error flag
  integer :: i_err

  ! Time parameters
  real(wp) :: time_start, time_end

  ! Index for command line argument
  integer :: i

  ! Start-up file name
  character(len=120) :: arg
  character(len=120) :: start_up_file_name

  continue

  ! Initialize mpi (mpimod)
  call mpiInit()

  ! Parse command line options
  do i = 1, command_argument_count(), 2
    call get_command_argument(i, arg)
     
    if (myprocid .eq. 0) then
      write(*,*)
      write(*,*) 'Parsing command line argument options'
      write(*,*) '==============================================================='
    end if
    select case (arg)
      case ('-sfn')
        call get_command_argument(i+1, arg)
        start_up_file_name = adjustl(arg)
        if (myprocid .eq. 0) then
          write(*,*) '  Start-up file name: ', start_up_file_name
          write(*,*) '==============================================================='
        end if
      case default
        if (myprocid .eq. 0) then
          write(*,*) '  Unrecognized command-line option: ', arg
          write(*,*) '  Exiting...'
          write(*,*)
          stop
        end if
    end select
  end do

  ! Default name of the input file
  if (command_argument_count() == 0) then
    start_up_file_name = 'SSDCstartup'
    if (myprocid .eq. 0) then
      write(*,*)
      write(*,*) 'Start-up file name: ', start_up_file_name
      write(*,*) '==============================================================='
    end if
  end if

  ! Wait 2 second before to continue
  call sleep(2)

  ! Start time
  time_start = MPI_WTIME()

  ! Read startup parameters
  call readSSDCstartup(start_up_file_name)

  ! If needed create directory for solution files
  if ((write_solution .eqv. .true.) .or. (write_aero_coeffs .eqv. .true.) .or. &
    & (write_dke_dt .eqv. .true.) .or. (write_errors_wall .eqv. .true.)) then
    if (myprocid .eq. 0) then
      call  create_solution_dir()
    endif
  endif

  ! If needed create directory for restarting file
  if (write_restart .eqv. .true.) then
    if (myprocid .eq. 0) then
      call  create_restart_dir()
    endif
  endif

  ! PDE independent setup: including
  ! read grid, define elements, parallel distribution, temporal scheme
  call physics_independent_setup()

  ! Initialize RK scheme 
  call Runge_Kutta_init()

  ! Initialize the physics dependent variables
  call Physics_init()

!  call cpu_time ( time_start )
  
  ! Call the physics solver
  call Physics_Solver()
  
!  call cpu_time ( time_end )
  
  ! Compute error
  call nse_calcerror(timeglobal) 

  ! End time
  time_end = MPI_WTIME()

  ! Output CPU time
  if (myprocid == 0) write(*,201) time_end-time_start
  201 format(' CPU time: ',F15.4)

  !  Master node writes useful message for locating the solution files and the
  !  restarting files
  if (myprocid .eq. 0) then

    if (write_solution .eqv. .true.) then
      write(*,*)
      write(*,*) 'Solution files written in ', write_solution_dir
      write(*,*)
    endif

    if (write_restart .eqv. .true.) then
      write(*,*)
      write(*,*) 'Restarting files written in ', write_restart_dir
      write(*,*)
    endif

  endif

  ! Finalize MPI and the hooks to PETSc
  call PetscFinalize(i_err)

end program ssdcsolver
