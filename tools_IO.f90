module tools_IO
  ! This module contains some useful subroutines and functions for IO
  ! operations

  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none

  private

  ! Public subroutine and functions
  public get_tag_proc
  public get_tag_time
  public get_tag_iter
  public get_write_restart_file_name
  public get_read_restart_file_name
  public get_file_unit
  public check_io_open_file
  public get_solution_proc_vtu_file_name
  public get_solution_pvtu_file_name
  public get_aerodynamic_coefficients_file_name
  public get_error_no_slip_wall_bc_file_name
  public get_dkinetic_energy_dt_file_name
  public get_enstrophy_file_name
  public get_time_space_errors_file_name
  public get_error_heat_entropy_flow_wall_bc_file_name
   

contains


  pure function get_tag_proc(proc_id)
    ! This function returns a string with the processor ID.

    use referencevariables, only : nprocs
    ! Nothing is defined implicitly
    implicit none

    integer, intent(in) :: proc_id
    character(60) :: get_tag_proc

    continue

    select case (nprocs)

      case(   1:   10)
        write(get_tag_proc,1001) proc_id
      case(  11:  100)
        write(get_tag_proc,1002) proc_id
      case( 101: 1000)
        write(get_tag_proc,1003) proc_id
      case(1001:10000)
        write(get_tag_proc,1004) proc_id
      case default
        write(get_tag_proc,1005) proc_id

    end select

    1001 format("",i1.1)
    1002 format("",i2.2)
    1003 format("",i3.3)
    1004 format("",i4.4)
    1005 format("",i5.5)

    return
  end function get_tag_proc

  !============================================================================
  
  !============================================================================

  pure function get_tag_time(global_time)
    ! This function returns a string with the global time.

    ! Nothing is defined implicitly
    implicit none

    real(wp), intent(in) :: global_time
    character(120) :: get_tag_time

    continue

    ! Assign string to get_tag_time
    write(get_tag_time,1001) global_time
    1001 format("",es20.13)

    ! Remove left white space from get_tag_time if any
    get_tag_time = adjustl(get_tag_time)

    return
  end function get_tag_time

  pure function get_tag_iter(iter_cnt)
    ! This function returns a string with the iteration (time-step) number.

    ! Nothing is defined implicitly
    implicit none

    integer, intent(in) :: iter_cnt
    character(120) :: get_tag_iter

    continue

    ! Assign string to get_tag_time
    write(get_tag_iter,1002) iter_cnt
    1002 format("",i8.8)

    ! Remove left white space from get_tag_time if any
    get_tag_iter = adjustl(get_tag_iter)

    return
  end function get_tag_iter

  !============================================================================
  
  !============================================================================

  pure function get_write_restart_file_name(tag_proc,tag_time)
    ! This function returns a string with the name of the file that is written 
    ! by each processor for the restarting procedure.

    ! Load modules
    use controlvariables, only : write_restart_dir, write_restart_common_name

    ! Nothing is defined implicitly
    implicit none

    character(Len=*), intent(in) :: tag_proc
    character(Len=*), intent(in) :: tag_time
    character(120) :: get_write_restart_file_name

    continue

    ! Construct the file name that has to be written by each processor
    get_write_restart_file_name =  trim(write_restart_dir) // '/' // trim(write_restart_common_name) &
      & // '_p' // trim(tag_proc) // '_t'// trim(tag_time) // ".dat"

    ! Remove left blank spaces from get_write_restart_file_name
    get_write_restart_file_name = adjustl(get_write_restart_file_name)

    return
  end function get_write_restart_file_name

  !============================================================================
  
  !============================================================================

  pure function get_read_restart_file_name(tag_proc)
    ! This function returns a string with the name of the file that is read by
    ! each processor for the restarting procedure.

    ! Load modules
    use controlvariables, only : read_restart_dir, read_restart_common_name, &
                                 & read_restart_time

    ! Nothing is defined implicitly
    implicit none

    character(Len=*), intent(in) :: tag_proc
    character(120) :: get_read_restart_file_name

    continue

    ! Construct the file name that has to be written by each processor
    get_read_restart_file_name = './' // trim(read_restart_dir) // '/' // trim(read_restart_common_name) &
      & // '_p' // trim(tag_proc) // '_t' // trim(read_restart_time) // ".dat"

    ! Remove left blank spaces from get_read_restart_file_name
    get_read_restart_file_name = adjustl(get_read_restart_file_name)

    return
  end function get_read_restart_file_name

  !============================================================================
  
  !============================================================================

  pure function get_solution_proc_vtu_file_name(tag_proc,tag_time)
    ! This function returns a string with the name of the file that contains 
    ! the solution written by each processor.

    ! Load modules
    use controlvariables, only : write_solution_dir, write_solution_common_name

    ! Nothing is defined implicitly
    implicit none

    character(Len=*), intent(in) :: tag_proc
    character(Len=*), intent(in) :: tag_time
    character(120) :: get_solution_proc_vtu_file_name

    continue

    ! Construct the file name that has to be written by each processor
    get_solution_proc_vtu_file_name = './' // trim(write_solution_dir) // '/' // trim(write_solution_common_name) &
      & // '_p' // trim(tag_proc) // '_t'// trim(tag_time) // ".vtu"

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_solution_proc_vtu_file_name = adjustl(get_solution_proc_vtu_file_name)

    return
  end function get_solution_proc_vtu_file_name

  !============================================================================
  
  !============================================================================

  pure function get_solution_pvtu_file_name(tag_proc,tag_time,tag_iter)
    ! This function returns a string with the name of the file that contains the 
    ! solution data for parallel visualization.

    ! Load modules
    use controlvariables, only : write_solution_dir, write_solution_common_name

    ! Nothing is defined implicitly
    implicit none

    character(Len=*), intent(in) :: tag_proc
    character(Len=*), intent(in) :: tag_time
    character(Len=*), intent(in) :: tag_iter
    character(120) :: get_solution_pvtu_file_name

    continue

    ! Construct the file name that has to be written by each processor
!    get_solution_pvtu_file_name = './' // trim(write_solution_dir) // '/' // trim(write_solution_common_name) &
!      & // '_p' // trim(tag_proc) // '_t'// trim(tag_time) // ".pvtu"

    get_solution_pvtu_file_name = './' // trim(write_solution_dir) // '/' // trim(write_solution_common_name) &
      & // '_' // trim(tag_iter) // ".pvtu"

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_solution_pvtu_file_name = adjustl(get_solution_pvtu_file_name)

    return
  end function get_solution_pvtu_file_name

  !============================================================================
  
  !============================================================================

  pure function get_aerodynamic_coefficients_file_name()
    ! This function returns a string with the name of the file that contains
    ! the aerodynamic coefficients. 

    ! Load modules
    use controlvariables, only : write_solution_dir
    
    ! Nothing is defined implicitly
    implicit none

    character(120) :: get_aerodynamic_coefficients_file_name
    character(120) :: file_name

    continue

    ! Define file name
    file_name = 'aerodynamic_coefficients.dat'

    get_aerodynamic_coefficients_file_name = './' // trim(write_solution_dir) // '/' // trim(file_name)

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_aerodynamic_coefficients_file_name = adjustl(get_aerodynamic_coefficients_file_name)

    return
  end function get_aerodynamic_coefficients_file_name

  !============================================================================
  
  !============================================================================

  pure function get_error_no_slip_wall_bc_file_name()
    ! This function returns a string with the name of file that contains the 
    ! no-slip wall boundary conditions error.

    ! Load modules
    use controlvariables, only : write_solution_dir
    
    ! Nothing is defined implicitly
    implicit none

    character(120) :: get_error_no_slip_wall_bc_file_name
    character(120) :: file_name

    continue

    ! Define file name
    file_name = 'error_no_slip_wall_bc.dat'

    get_error_no_slip_wall_bc_file_name = './' // trim(write_solution_dir) // '/' // trim(file_name)

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_error_no_slip_wall_bc_file_name = adjustl(get_error_no_slip_wall_bc_file_name)

    return
  end function get_error_no_slip_wall_bc_file_name

  !============================================================================
  
  !============================================================================

  pure function get_dkinetic_energy_dt_file_name()
    ! This function returns a string with the name of the file that contains
    ! the time derivative of the kinetic energy. 

    ! Load modules
    use controlvariables, only : write_solution_dir
    
    ! Nothing is defined implicitly
    implicit none

    character(120) :: get_dkinetic_energy_dt_file_name
    character(120) :: file_name

    continue

    ! Define file name
    file_name = 'dke_dt.dat'

    get_dkinetic_energy_dt_file_name = './' // trim(write_solution_dir) // '/' // trim(file_name)

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_dkinetic_energy_dt_file_name = adjustl(get_dkinetic_energy_dt_file_name)

    return
  end function get_dkinetic_energy_dt_file_name

  !============================================================================
  
  !============================================================================

  pure function get_enstrophy_file_name()
    ! This function returns a string with the name of the file that contains
    ! the enstrophy. 

    ! Load modules
    use controlvariables, only : write_solution_dir
    
    ! Nothing is defined implicitly
    implicit none

    character(120) :: get_enstrophy_file_name
    character(120) :: file_name

    continue

    ! Define file name
    file_name = 'enstrophy.dat'

    get_enstrophy_file_name = './' // trim(write_solution_dir) // '/' // trim(file_name)

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_enstrophy_file_name = adjustl(get_enstrophy_file_name)

    return
  end function get_enstrophy_file_name

  !============================================================================
  
  !============================================================================

  pure function get_error_heat_entropy_flow_wall_bc_file_name()
    ! This function returns a string with the name of file that contains the 
    ! no-slip wall boundary conditions error.

    ! Load modules
    use controlvariables, only : write_solution_dir
    
    ! Nothing is defined implicitly
    implicit none

    character(120) :: get_error_heat_entropy_flow_wall_bc_file_name
    character(120) :: file_name

    continue

    ! Define file name
    file_name = 'error_heat_entropy_flow_wall_bc.dat'

    get_error_heat_entropy_flow_wall_bc_file_name = './' // trim(write_solution_dir) // '/' // trim(file_name)

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_error_heat_entropy_flow_wall_bc_file_name = adjustl(get_error_heat_entropy_flow_wall_bc_file_name)

    return
  end function get_error_heat_entropy_flow_wall_bc_file_name


  !============================================================================
  
  !============================================================================

   function get_time_space_errors_file_name()
    ! This function returns a string with the name of the file that contains
    ! the space and time errors. 

    ! Load modules
    use controlvariables, only : casefile, write_solution_dir
    
    ! Nothing is defined implicitly
    implicit none

    character(300) :: get_time_space_errors_file_name
    character(300) :: file_name

    continue

    ! Define file name
    file_name = 'time_space_errors.' // trim(casefile) //'.data'

    get_time_space_errors_file_name = './' // trim(write_solution_dir) // '/' // trim(file_name)

    ! Remove left blank spaces from get_solution_proc_vtu_file_name
    get_time_space_errors_file_name = adjustl(get_time_space_errors_file_name)

    return
  end function get_time_space_errors_file_name


  !============================================================================
  
  !============================================================================


  function get_file_unit(max_unit)
    ! This function returns a unit number that is not in use

    ! Nothing is defined implicitly
    implicit none

    integer, intent(in) :: max_unit
    integer :: i_unit
    integer :: j
    integer :: iostat
    logical :: opened
    integer :: get_file_unit

    continue

    ! Check which unit is free in a descending order
    j = max_unit

    do i_unit = j,1,-1
      inquire(unit=i_unit,opened=opened,iostat=iostat)
      if (iostat .ne. 0) then 
        cycle
        if (.not. opened) then 
          exit
        endif
      endif
    enddo

    ! Assign free unit to get_file_unit
    get_file_unit = i_unit

    return
  end function get_file_unit

  !============================================================================
  
  !============================================================================

  subroutine check_io_open_file(io_status,proc_ID,message)
    ! This subroutine check whether io_status of the Fortran open instruction
    
    !Nothing is implicitly defined
    implicit none

    integer, intent(in) :: io_status
    integer, intent(in) :: proc_ID
    character(Len=*), intent(in) :: message

    continue

    ! Output at screen an error message if the io_status is "negative"
    if(io_status > 0) then
      write(*,*)
      write(*,*) trim(message) 
      write(*,*) 'Processor: ', proc_ID
      write(*,*) 'Exiting...'
      write(*,*)
      stop
    endif

    return
  end subroutine check_io_open_file

end module tools_IO
