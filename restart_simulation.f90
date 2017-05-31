! This module contains the necessary routines to write the restart and the
! solution files.

module restart_simulation

  ! Load modules
  use precision_vars
  use tools_IO

  ! Nothing is implicitly defined
  implicit none


contains

  !============================================================================

  !============================================================================
  ! create_restart_dir - Creates the restart directory specified by 
  !                      write_restart_dir.
  
  subroutine create_restart_dir()

    ! Load module
    use mpimod
    use referencevariables
    use controlvariables
      
    ! Nothing is implicitly defined
    implicit None
    
    logical :: lexist

    ! Check if the directory already exist
    inquire(file=write_restart_dir,exist=lexist)

    ! If it does not exist then create it
    if (.not.lexist) then
      call system('mkdir ' // write_restart_dir)
    endif
    
    return
  end subroutine create_restart_dir

  !============================================================================

  !============================================================================
  ! write_restart_file - Writes to a file (specified by "common_restart_name" 
  !                      and saved in the "restart" subdirectory) the values of 
  !                      the conserved variables in each solution point of each 
  !                      cell and the value of the time-averaged quantities if 
  !                      they have been calculated.
  !                      Each processor executes this routines. Therefore N 
  !                      files, with N = number of processors, are written. 
  !                      This implies that when the restarting procedure is 
  !                      used, the same number of processors adopted in the 
  !                      first part of the simulation must be used.

  subroutine write_restart_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use variables, only: ug, mean_vg, time_ave_prod_vel_comp

    ! Nothing is implicitly defined
    implicit none

    character(120) :: file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: elem_low, elem_high, i_elem, i_node
    integer :: i_unit, io_status
    integer :: max_unit
    character(120) :: message

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the restart solution.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Get restart file name
    file_name = get_write_restart_file_name(tag_proc,tag_time)

    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    
    if (write_restart_formatted) then
      ! Open file
      open(unit=i_unit,file=file_name,status='unknown',iostat=io_status)

      ! Check io_status
      call check_io_open_file(io_status,myprocid,message)

      ! Write Order of poly to restart file
      write(i_unit,*) npoly

      ! Write number of time steps executed
      write(i_unit,*) itimestep

      ! Writes errors needed by the time-step controller
      write(i_unit,*) err_space_lf
      write(i_unit,*) err_time_lf
      write(i_unit,*) err_time_tm1
      write(i_unit,*) err_time_tm2

      ! Write element ID and conserved variables to file
      do i_elem = elem_low,elem_high
        write(i_unit,*) i_elem
        do i_node = 1, nodesperelem
          write(i_unit,*) ug(:,i_node,i_elem)
        enddo
      enddo

      ! Write element ID and time-averaged quantities to file
      if (time_averaging) then
        do i_elem = elem_low,elem_high
          write(i_unit,*) i_elem
          do i_node = 1, nodesperelem
            write(i_unit,*) mean_vg(:,i_node,i_elem)
          enddo
        enddo

        do i_elem = elem_low,elem_high
          write(i_unit,*) i_elem
          do i_node = 1, nodesperelem
            write(i_unit,*) time_ave_prod_vel_comp(:,i_node,i_elem)
          enddo
        enddo

      endif

!      ! Format for writing the formatted output
!      1002 format(I8)
!      1003 format(5(es17.10,1x))
!      1004 format(6(es17.10,1x))

    else

      ! Open file
      open(unit=i_unit,file=file_name,form='UNFORMATTED',status='unknown', &
        & iostat=io_status)

      ! Check io_status
      call check_io_open_file(io_status,myprocid,message)

      ! Write Order of poly to restart file
      write(i_unit) npoly

      ! Write number of time steps executed
      write(i_unit) itimestep

      ! Writes errors needed by the time-step controller
      write(i_unit) err_space_lf
      write(i_unit) err_time_lf
      write(i_unit) err_time_tm1
      write(i_unit) err_time_tm2

      ! Write element ID and conserved variables to file
      do i_elem = elem_low,elem_high
        write(i_unit) i_elem
        do i_node = 1, nodesperelem
          write(i_unit) ug(:,i_node,i_elem)
        enddo
      enddo

      ! Write element ID and time-averaged quantitiesi to file
      if (time_averaging) then

        do i_elem = elem_low,elem_high
          write(i_unit) i_elem
          do i_node = 1, nodesperelem
            write(i_unit) mean_vg(:,i_node,i_elem)
          enddo
        enddo

        do i_elem = elem_low,elem_high
          write(i_unit) i_elem
          do i_node = 1, nodesperelem
            write(i_unit) time_ave_prod_vel_comp(:,i_node,i_elem)
          enddo
        enddo

      endif

    endif

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_restart_file

  !============================================================================

  !============================================================================
  ! read_restart_file - Reads a file (specified by "common_restart_name" 
  !                      and saved in the "restart" subdirectory) the values of 
  !                      the conserved variables in each solution point of each 
  !                      cell and the value of the time-averaged quantities if 
  !                      they have been calculated.
  !                      Each processor executes this routines. Therefore N 
  !                      files, with N = number of processors, are read. This 
  !                      implies that when the restarting procedure is used, the
  !                      same number of processors adopted in the first part of
  !                      the simulation must be used.

  subroutine read_restart_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use variables, only: ug, mean_vg, time_ave_prod_vel_comp
    use time_average
    use initcollocation
    use SSWENO_routines, only: Negative_Density_Removal

    ! Nothing is implicitly defined
    implicit none

    integer                                 :: nXA, nXB, k, nodesperelem_tmp
    real(wp), allocatable, dimension(:,:,:) :: ug_tmp
    real(wp), allocatable, dimension(:)     ::  XA,  XB, tmpA, tmpB, tmpfieldA

    logical :: lexist
    character(120) :: file_name
    character(120) :: tag_proc
    integer :: elem_low, elem_high, i_elem, i_node
    integer :: i_unit, io_status
    integer :: elem_ID
    integer :: max_unit
    integer :: npoly_tmp
    character(120) :: message

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for reading the restart solution.'

    ! Set max number of unit
    max_unit = 98

    ! Check if read_restart_dir exist
    if (myprocid == 0) then
      ! Check if the directory already exist
      INQUIRE(FILE=read_restart_dir,EXIST=lexist)

      ! If it does not exist stop execution of the code
      if (.not.lexist) then
        write(*,*)
        write(*,*) 'Error: read_restart_dir does not exist!'
        write(*,*) 'Check the SSDCstartup file!'
        write(*,*) 'Exiting...'
        write(*,*)
        stop
      endif
    endif

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get restart file name
    file_name = get_read_restart_file_name(tag_proc)

    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    if (read_restart_formatted) then
      
      ! Open file
      open(unit=i_unit,file=file_name, status='old',action='READ',&
        & iostat=io_status)
      
      ! Check io_status
      call check_io_open_file(io_status,myprocid,message)

      ! Read Order of poly in restart file
      read(i_unit,*) npoly_tmp
!     npoly_tmp = npoly

      nodesperelem_tmp = (npoly_tmp+1)**(ndim)
      allocate(ug_tmp(nequations,nodesperelem_tmp,elem_low:elem_high))

      ! Read number of time steps used in the previous run 
      read(i_unit,*) restart_time_steps

      ! Read errors needed by the time-step controller
      read(i_unit,*) err_space_lf
      read(i_unit,*) err_time_lf
      read(i_unit,*) err_time_tm1
      read(i_unit,*) err_time_tm2

      ! Read element ID and conserved variables from file
      do i_elem = elem_low,elem_high
        read(i_unit,*) elem_ID
        if (elem_ID .ne. i_elem) then
          write(*,*) 'Failure in reading the restarting solution. The', &
            & ' element ID written in the restarting file does not', &
            & ' correspond to the element ID defined in the loop.', &
            & ' Processor: ', myprocid
          write(*,*) 'Exiting...'
          write(*,*)
          stop
        endif

        do i_node = 1, nodesperelem_tmp
          read(i_unit,*) ug_tmp(:,i_node,i_elem)
        enddo
      enddo

      if ((read_restart_time_averaging) .and. (npoly == npoly_tmp)) then

        ! Read element ID and time-averaged quantities from file
        do i_elem = elem_low,elem_high
          read(i_unit,*) elem_ID
          if (elem_ID .ne. i_elem) then
            write(*,*) 'Failure in reading the restarting solution. The', &
              & ' element ID written in the restarting file does not', &
              & ' correspond to the element ID defined in the loop.',&
              & ' Processor: ', myprocid
            write(*,*) 'Exiting...'
            write(*,*)
            stop
          endif

          do i_node = 1, nodesperelem
            read(i_unit,*) mean_vg(:,i_node,i_elem)
          enddo
        enddo

        do i_elem = elem_low,elem_high
          read(i_unit,*) elem_ID
          if (elem_ID .ne. i_elem) then
            write(*,*) 'Failure in reading the restarting solution. The', &
              & ' element ID written in the restarting file does not', &
              & ' correspond to the element ID defined in the loop.', &
              & ' Processor: ', myprocid
            write(*,*) 'Exiting...'
            write(*,*)
            stop
          endif

          do i_node = 1, nodesperelem
            read(i_unit,*) time_ave_prod_vel_comp(:,i_node,i_elem)
          enddo
        enddo

        ! Compute Reynolds stress from the time-averaged quantities just read 
        call compute_reynolds_stress()

      endif

!      ! Format for reading the file
!      ! This format must be the same of that used to write the restarting file 
!      ! in the subroutine "write_restart_file_each_process"
!      1005 format(5(es17.10,1x))
!      1006 format(6(es17.10,1x))


    else

      ! Open file
      open(unit=i_unit,file=file_name,form='UNFORMATTED',status='old', &
        & action='READ',iostat=io_status)

      ! Check io_status
      call check_io_open_file(io_status,myprocid,message)

      ! Read Order of poly in restart file
      read(i_unit) npoly_tmp

      nodesperelem_tmp = (npoly_tmp+1)**(ndim)
      allocate(ug_tmp(nequations,nodesperelem_tmp,elem_low:elem_high))

      ! Read number of time steps used in the previous run 
      read(i_unit) restart_time_steps

      ! Read errors needed by the time-step controller
      read(i_unit) err_space_lf
      read(i_unit) err_time_lf
      read(i_unit) err_time_tm1
      read(i_unit) err_time_tm2

      ! Read element ID and conserved variables from file
      do i_elem = elem_low,elem_high
        read(i_unit) elem_ID
        if (elem_ID .ne. i_elem) then
          write(*,*) 'Failure in reading the restarting solution. The', &
            & ' element ID written in the restarting file does not', &
            & ' correspond to the element ID defined in the loop.', &
            & ' Processor: ', myprocid
          write(*,*) 'Exiting...'
          write(*,*)
          stop
        endif

        do i_node = 1, nodesperelem_tmp
          read(i_unit) ug_tmp(:,i_node,i_elem)
        enddo
      enddo

      if ((read_restart_time_averaging) .and. (npoly == npoly_tmp)) then

        if (npoly /= npoly_tmp) then
            write(*,*) 'Failure in reading the restarting solution. ', &
              & ' Interpolation not implemented for averaged quantities ', &
              & ' Processor: ', myprocid
          write(*,*) 'Exiting...'
          write(*,*)
          stop
        endif

        ! Read element ID and time-averaged quantities from file
        do i_elem = elem_low,elem_high
          read(i_unit) elem_ID
          if (elem_ID .ne. i_elem) then
            write(*,*) 'Failure in reading the restarting solution. The', &
              & ' element ID written in the restarting file does not', &
              & ' correspond to the element ID defined in the loop.', &
              & ' Processor: ', myprocid
            write(*,*) 'Exiting...'
            write(*,*)
            stop
          endif

          do i_node = 1, nodesperelem
            read(i_unit) mean_vg(:,i_node,i_elem)
          enddo
        enddo


        do i_elem = elem_low,elem_high
          read(i_unit) elem_ID
          if (elem_ID .ne. i_elem) then
            write(*,*) 'Failure in reading the restarting solution. The', &
              & ' element ID written in the restarting file does not', &
              & ' correspond to the element ID defined in the loop.', &
              & ' Processor: ', myprocid            
            write(*,*) 'Exiting...'
            write(*,*)
            stop
          endif

          do i_node = 1, nodesperelem
            read(i_unit) time_ave_prod_vel_comp(:,i_node,i_elem)
          enddo
        enddo

        ! Compute Reynolds stress from the time-averaged quantities just read 
        call compute_reynolds_stress()

      endif

    endif

    ! Interpolate to higher dimensional space if npoly /= npoly_tmp

    if    (npoly_tmp == npoly) then
      ug(:,:,:) = ug_tmp(:,:,:)
    else  
      nXA = npoly_tmp+1 ; allocate(XA(nXA), tmpA(nXA)) ;
      nXB = npoly    +1 ; allocate(XB(nXB), tmpB(nXB)) ;

      call Gauss_Lobatto_Legendre_points(nXA,XA,tmpA)
      call Gauss_Lobatto_Legendre_points(nXB,XB,tmpB)

      allocate(tmpfieldA( (nXA)**ndim ))

      do i_elem = elem_low,elem_high

        do k = 1,nequations
          tmpfieldA(:) = ug_tmp(k,:,i_elem)
          call ExtrpXA2XB(ndim,nXA,nXB,XA,XB,tmpfieldA(:),ug(k,:,i_elem))
        enddo

        do i_node = 1, nodesperelem
          call Negative_Density_Removal(ug(:,i_node,i_elem),nequations)
        enddo

      enddo

      if (read_restart_time_averaging) then
          mean_vg(:,:,:) = 0.0_wp
          time_ave_prod_vel_comp(:,:,:) = 0.0_wp
      endif

      deallocate(ug_tmp,tmpfieldA,tmpA,tmpB)


    endif


    ! Close file unit
    close(i_unit)

    ! Set timeglobal equal to the restarting value, i.e. the time associated to
    ! the solution used for restarting
    read(read_restart_time,*) timeglobal

    return
  end subroutine read_restart_file

  !============================================================================

end module restart_simulation

