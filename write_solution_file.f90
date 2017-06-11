! This module contains the necessary routine to write the solution to file
! using the new (XML) format for. Such files can be read in ParaView
! See the following links for more information regarding the software and
! the data formats
!
! - www.paraview.org 
! - http://paraview.org/Wiki/ParaView/Data_formats for more information 
!

module write_solution_file

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Subroutines in this module are generally private 
  private

  ! List of public subroutines
  public create_solution_dir
  public write_solution_vtu_file
  public write_solution_pvtu_file
  public write_aerodynamic_coefficients
  public write_error_no_slip_wall_bc
  public write_dkinetic_energy_dt
  public write_enstrophy
  public write_error_heat_entropy_flow_wall_bc
  public calculate_binary_sizes_vtu
  public Sods_Line_Plot

contains

  !============================================================================

  !============================================================================
  ! create_solution_dir - Creates solution directory specified by 
  !                       write_solution_dir.

  subroutine create_solution_dir()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
      
    ! Nothing is implicitly defined
    implicit none

    logical :: lexist

    continue
    
    ! Check if the directory already exist
    inquire(file=write_solution_dir,exist=lexist)

    ! If it does not exist then create it
    if (.not.lexist) then
      call system('mkdir ' // write_solution_dir)
    endif

    return
  end subroutine create_solution_dir

  !============================================================================

  !============================================================================
  ! calculate_binary_sizes_vtu - Computes the length (in bytes) of the binray
  !                              values for writing the raw binary .vtu file.
  subroutine calculate_binary_sizes_vtu()

    ! Load modules
    use binary_sizes
    use referencevariables
    use controlvariables, only : write_solution_formatted
    
    ! Nothing is implicilty defined
    implicit none

    real(wp) :: tmp
    integer :: iell, elem_high
    integer :: n_elems, n_subelems, n_vertices_sub_elem
    
    continue
    
    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high + 1 - iell

    ! Number of sub-elements
    n_subelems = n_elems*(nodesperedge-1)**ndim

    ! Number of vertices per sub-element
    n_vertices_sub_elem = 8

    ! Compute byte lengths
    if (write_solution_formatted .eqv. .false.) then
      n_bytes_scal   = nodesperelem*n_elems*sizeof(tmp)
      n_bytes_vec    = ndim*nodesperelem*n_elems*sizeof(tmp)
      n_bytes_xyz    = ndim*n_subelems*n_vertices_sub_elem*n_elems*sizeof(tmp)
      n_bytes_ien    = n_subelems*n_vertices_sub_elem*sizeof(n_elems)
      n_bytes_offset = n_elems*sizeof(n_elems)
      n_bytes_etype  = n_elems*sizeof(n_elems)
    end if

  end subroutine calculate_binary_sizes_vtu
  !============================================================================

  !============================================================================
  ! write_solution_vtu_file - This subroutine is the driver subroutine for 
  !                           writing the solution to the .vtu file.

  subroutine write_solution_vtu_file()
  
    ! Load modules
    use referencevariables
    use controlvariables
    use variables
!   use navierstokes, only : supersonicvortexFull

!   integer :: ielem, inode

!   real(wp), dimension(3)   :: ctmp
!   real(wp), dimension(5)   :: fvtmp
!   real(wp), dimension(5,3) :: phitmp

!   real(wp), dimension(1:nequations,nodesperelem,ihelems(1):ihelems(2)) :: exact
!   exact = 0.0_wp

!   do ielem = ihelems(1), ihelems(2)
!     ! Loop over nodes in each element
!     do inode = 1, nodesperelem
!       ! Use exact solution routine to initialize data
!       call supersonicvortexFull( &
!         exact(:,inode,ielem), &
!         phitmp, &
!         fvtmp, &
!         ctmp, &
!         xg(:,inode,ielem), &
!         0.0_wp, &
!         nequations, &
!         ndim, &
!         mut(inode,ielem)) ! (navierstokes)
!     enddo
!   enddo

    ! Write ASCII file
    if (write_solution_formatted) then

!     exact(:,:,:) = log10(abs(exact(:,:,:) - vg(:,:,:))+1.0e-16_wp)

      ! Write header file
      call write_header_ascii_vtu_file()

      ! Write XML string for appending point data
      call write_begin_point_data_ascii_vtu_file()
      
      ! Write density to file
      call write_scalar_ascii_vtu_file(vg(1,:,:),'Density')
!     call write_scalar_ascii_vtu_file(exact(1,:,:),'Density Error')

      ! Write temperature to file
      call write_scalar_ascii_vtu_file(vg(5,:,:),'Temperature')
!     call write_scalar_ascii_vtu_file(exact(5,:,:),'Temperature Error')

      ! Write x component of the momentum to file
      !call write_scalar_ascii_vtu_file(ug(2,:,:),'momentum x')

      ! Write y component of the momentum to file
      !call write_scalar_ascii_vtu_file(ug(3,:,:),'momentum y')

      ! Write z component of the momentum to file
      !call write_scalar_ascii_vtu_file(ug(4,:,:),'momentum z')
      
      ! Write momentum to file
      call write_vector_3d_ascii_vtu_file(ug(2:4,:,:),'Momentum') 

      ! Write total energy to file
      call write_scalar_ascii_vtu_file(ug(5,:,:),'Specific total energy')

      ! Write velocity to file
      call write_vector_3d_ascii_vtu_file(vg(2:4,:,:),'Velocity')
!     call write_vector_3d_ascii_vtu_file(exact(2:4,:,:),'Velocity Error')

      ! Write vorticity to file
      call write_vector_3d_ascii_vtu_file(omega(1:3,:,:),'Vorticity')

      
      ! Output time-averaged quantities if necessary
      if (time_averaging) then

        ! Write the mean of the density to file
        call write_scalar_ascii_vtu_file(mean_vg(1,:,:),'Mean density')

        ! Write the mean of the velocity to file
        call write_vector_3d_ascii_vtu_file(mean_vg(2:4,:,:),'Mean velocity')

        ! Write the mean of the temperature to file
        call write_scalar_ascii_vtu_file(mean_vg(5,:,:),'Mean temperature')

        ! Reynolds stress component <u'*u'>
        call write_scalar_ascii_vtu_file(reynolds_stress(1,:,:),'up up')

        ! Reynolds stress component <u'*v'>
        call write_scalar_ascii_vtu_file(reynolds_stress(2,:,:),'up vp')
      
        ! Reynolds stress component <u'*w'>
        call write_scalar_ascii_vtu_file(reynolds_stress(3,:,:),'up wp')
      
        ! Reynolds stress component <v'*v'>
        call write_scalar_ascii_vtu_file(reynolds_stress(4,:,:),'vp vp')
      
        ! Reynolds stress component <v'*w'>
        call write_scalar_ascii_vtu_file(reynolds_stress(5,:,:),'vp wp')
      
        ! Reynolds stress component <w'*w'>
        call write_scalar_ascii_vtu_file(reynolds_stress(6,:,:),'wp wp')

      endif

      ! Write specific entropy to file
      call write_scalar_ascii_vtu_file(specific_entropy(:,:),'Specific entropy')


      ! Write XML string for closing the point data
      call write_end_point_data_ascii_vtu_file()
      
      ! Writes nodes coordinates
      call write_nodes_discontinous_elements_ascii_vtu_file()
      
      ! Write connectivity list 
      call write_connectivity_discontinous_elements_ascii_vtu_file()
      
      ! Write footer file
      call write_footer_ascii_vtu_file()

      ! Write RAW BINARY file (data are appended at the end of the .vtu file) 
    else

      call write_xml_raw_binary_vtu_file()

    end if

    return
  end subroutine write_solution_vtu_file

  !============================================================================

  !============================================================================
  ! write_solution_pvtu_file - This subroutine is the driver subroutine for 
  !                            writing the .pvtu file.

  subroutine write_solution_pvtu_file()

    ! Load modules
    use variables

    ! Nothing is implicitly defined
    implicit none

    ! Write header file
    call write_parallel_vtu_file()

    return
  end subroutine write_solution_pvtu_file

  !============================================================================

  !============================================================================
  ! write_header_ascii_vtu_file - Writes the .vtu header file for unstructured 
  ! grids. ASCII file.

  subroutine write_header_ascii_vtu_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO
      
    ! Nothing is implicitly defined
    implicit none
    
    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_elems, n_subelems, n_nodes
    character(120) :: message

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the header of the &
      &  .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name = get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
    
    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='unknown',iostat=io_status)
    
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)
    
    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high+1-iell
    
    ! Number of sub-elements
    n_subelems = n_elems*(nodesperedge-1)**ndim

    ! Number of nodes or vertices necessary to write the solution has
    ! discontinuous field also within an element
    if(ndim == 2) then
      n_nodes = n_subelems*4 ! For quadrilateral elements (2D)
    else
      n_nodes = n_subelems*8 ! For hexahedral elements (3D)
    endif
  
    ! Write header file
    ! Use version="1.0" for XML. One could also use version 0.1
    write(i_unit,'(A)') '<?xml version="1.0"?>'
    
    write(i_unit,'(A)') '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">'
    
    write(i_unit,*) '<UnstructuredGrid>'
    
    write(i_unit,*) '<Piece NumberOfPoints="',n_nodes, &
      '" NumberOfCells="',n_subelems,'">'
    
    ! Close file unit
    close(i_unit)
    
    return
  end subroutine write_header_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_footer_ascii_vtu_file - Writes the .vtu footer file for unstructured 
  ! grids. ASCII file.

  subroutine write_footer_ascii_vtu_file()
     
    ! Load modules
    use referencevariables
    use controlvariables
    use tools_IO
      
    ! Nothing is implicitly defined
    implicit none
    
    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    character(120) :: message

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the footer of the .vtu &
      & file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name =  get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
    
    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
    
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write footer file
    write(i_unit,*)'</Piece>'
    write(i_unit,*)'</UnstructuredGrid>'
    write(i_unit,*)'</VTKFile>'
    
    ! Close file unit
    close(i_unit)
   
    return
  end subroutine write_footer_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_begin_point_data_ascii_vtu_file - Writes the XML instructions for 
  ! appending the point data in the .vtu file. ASCII file.

  subroutine write_begin_point_data_ascii_vtu_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    character(120) :: message

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for starting point data XML &
      & instruction in the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name =  get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
    
    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
    
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write XML instruction
    write(i_unit,*)'<PointData>'

    ! Close file unit
    close(i_unit)
   
    return
  end subroutine write_begin_point_data_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_end_point_data_ascii_vtu_file - Writes the XML instructions for 
  ! closing the point data in the .vtu file. ASCII file.

  subroutine write_end_point_data_ascii_vtu_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    character(120) :: message

    message = 'Failure in opening the file for closing the point data XML &
      & instruction in the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name =  get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
    
    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
    
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write XML instruction
    write(i_unit,*)'</PointData>'   

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_end_point_data_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_scalar_ascii_vtu_file - Writes a scalar field to the .vtu file for 
  ! unstructured grids. ASCII file.
  !
  ! Input parameters:
  ! var - 1D array of values
  ! var_name - name of the variable

  subroutine write_scalar_ascii_vtu_file(var,var_name)

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_elems, n_subelems_1d
    character(120) :: message
    character(Len=*),intent(in) :: var_name
    real(wp),intent(in) :: var(:,:)
    integer :: i_elem, i_row_face, i_subelem, i_slice
    integer :: local_node_1_idx, local_node_2_idx, local_node_3_idx, &
               & local_node_4_idx
    integer :: local_node_5_idx, local_node_6_idx, local_node_7_idx, &
               & local_node_8_idx

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high+1-iell
    
    ! Number of sub-elements in 1D
    n_subelems_1D = nodesperedge-1

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the scalar field'
    message = message // var_name // ' to the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name = get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
  
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write field header
    write(i_unit,*)'<DataArray type="Float64" Name="',var_name,'" format="ascii">'

    if(ndim == 2) then

      do i_elem = 1, n_elems
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the four nodes that form the quadrilateral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by four nodes
            local_node_1_idx = i_subelem &
              & + nodesperedge*(i_row_face-1) 

            local_node_2_idx = i_subelem + 1 &
              & + nodesperedge*(i_row_face-1)

            local_node_3_idx = i_subelem + nodesperedge + &
              & 1 + nodesperedge*(i_row_face-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1)

            ! Write scalar field of the four nodes
            write(i_unit,*) var(local_node_1_idx,i_elem)
            write(i_unit,*) var(local_node_2_idx,i_elem)
            write(i_unit,*) var(local_node_3_idx,i_elem)
            write(i_unit,*) var(local_node_4_idx,i_elem)

          enddo
        enddo
      enddo

    else

      do i_elem = 1, n_elems 
        do i_slice = 1, n_subelems_1D
          do i_row_face = 1, n_subelems_1D
            do i_subelem = 1, n_subelems_1d

              ! Construct indices of the eight nodes that form the hexahedral
              ! This algorithm should also work when we will have curved 
              ! elements because each subelement is made up by eight nodes
              local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1) 

              local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1)

              local_node_3_idx = i_subelem + nodesperedge + 1 + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

              local_node_4_idx = i_subelem + nodesperedge + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
          
              local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice

              local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice 

              local_node_7_idx = i_subelem + nodesperedge + 1 + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

              local_node_8_idx = i_subelem + nodesperedge + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

              ! Write scalar field of the eight nodes
              write(i_unit,*) var(local_node_1_idx,i_elem)
              write(i_unit,*) var(local_node_2_idx,i_elem)
              write(i_unit,*) var(local_node_3_idx,i_elem)
              write(i_unit,*) var(local_node_4_idx,i_elem)
          
              write(i_unit,*) var(local_node_5_idx,i_elem)
              write(i_unit,*) var(local_node_6_idx,i_elem)
              write(i_unit,*) var(local_node_7_idx,i_elem)
              write(i_unit,*) var(local_node_8_idx,i_elem)

            enddo
          enddo
        enddo
      enddo

    endif

     
    write(i_unit,*)'</DataArray>'
    
    ! Close file unit
    close(i_unit) 

    return
  end subroutine write_scalar_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_vector_2d_ascii_vtu_file - Writes 2D vector field to the .vtu file  
  ! for unstructured grids. ASCII file.
  !
  ! Input parameters:
  ! var - 2D array of values
  ! var_name - name of the variable

  subroutine write_vector_2d_ascii_vtu_file(var,var_name)
    ! This subroutine writes a 2D vector field to the .vtu file for unstructured
    ! grids

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_elems, n_subelems_1d
    character(120) :: message
    character(Len=*),intent(in) :: var_name
    real(wp),intent(in) :: var(:,:,:)
    integer :: i_elem, i_row_face, i_subelem, i_slice
    integer :: local_node_1_idx, local_node_2_idx, local_node_3_idx, &
               & local_node_4_idx
    integer :: local_node_5_idx, local_node_6_idx, local_node_7_idx, &
               & local_node_8_idx

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high+1-iell
    
    ! Number of sub-elements in 1D
    n_subelems_1D = nodesperedge-1

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the 2D vector field'
    message = message // var_name // ' to the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name = get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
    
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write 2D vector field
    write(i_unit,*)'<DataArray type="Float64" Name="',var_name,'" NumberOfComponents="2"  format="ascii">'
    
    if(ndim == 2) then
      
      do i_elem = 1, n_elems
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the four nodes that form the quadrilateral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by four nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) 
            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1)
            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1)
            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1)

            ! Write 2D vector field of the four nodes
            write(i_unit,*) var(1,local_node_1_idx,i_elem), &
              & var(2,local_node_1_idx,i_elem)

            write(i_unit,*) var(1,local_node_2_idx,i_elem), &
              & var(2,local_node_2_idx,i_elem) 

            write(i_unit,*) var(1,local_node_3_idx,i_elem), &
              & var(2,local_node_3_idx,i_elem)

            write(i_unit,*) var(1,local_node_4_idx,i_elem), &
              & var(2,local_node_4_idx,i_elem)

          enddo
        enddo
      enddo

    else

      do i_elem = 1, n_elems
        do i_slice = 1, n_subelems_1D
          do i_row_face = 1, n_subelems_1D
            do i_subelem = 1, n_subelems_1d

              ! Construct indices of the eight nodes that form the hexahedral
              ! This algorithm should also work when we will have curved 
              ! elements because each subelement is made up by eight nodes
              local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1) 
              
              local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1)
              
              local_node_3_idx = i_subelem + nodesperedge + 1 + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
              
              local_node_4_idx = i_subelem + nodesperedge + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
            
              local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice
              
              local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice 
              
              local_node_7_idx = i_subelem + nodesperedge + 1 &
                & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 
              
              local_node_8_idx = i_subelem + nodesperedge &
                & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

              ! Write 3D vector field of the eight nodes
              write(i_unit,*) var(1,local_node_1_idx,i_elem), &
                & var(2,local_node_1_idx,i_elem) 

              write(i_unit,*) var(1,local_node_2_idx,i_elem), &
                & var(2,local_node_2_idx,i_elem)

              write(i_unit,*) var(1,local_node_3_idx,i_elem), &
                & var(2,local_node_3_idx,i_elem)

              write(i_unit,*) var(1,local_node_4_idx,i_elem), &
                & var(2,local_node_4_idx,i_elem)

            
              write(i_unit,*) var(1,local_node_5_idx,i_elem), &
                & var(2,local_node_5_idx,i_elem)

              write(i_unit,*) var(1,local_node_6_idx,i_elem), &
                & var(2,local_node_6_idx,i_elem)

              write(i_unit,*) var(1,local_node_7_idx,i_elem), &
                & var(2,local_node_7_idx,i_elem)

              write(i_unit,*) var(1,local_node_8_idx,i_elem), &
                & var(2,local_node_8_idx,i_elem)

            enddo
          enddo
        enddo
      enddo
    
    endif
    
    write(i_unit,*)'</DataArray>'
    
    ! Close file unit
    close(i_unit) 
   
    return
  end subroutine write_vector_2d_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_vector_3d_ascii_vtu_file - Writes 3D vector field to the .vtu file 
  ! for unstructured grids. ASCII file.
  !
  ! Input parameters:
  ! var - 3D array of values
  ! var_name - name of the variable


  subroutine write_vector_3d_ascii_vtu_file(var,var_name)

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_elems, n_subelems_1d
    character(120) :: message
    character(Len=*),intent(in) :: var_name
    real(wp),intent(in) :: var(:,:,:)
    integer :: i_elem, i_row_face, i_subelem, i_slice
    integer :: local_node_1_idx, local_node_2_idx, local_node_3_idx, &
               & local_node_4_idx
    integer :: local_node_5_idx, local_node_6_idx, local_node_7_idx, &
               & local_node_8_idx

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high+1-iell
    
    ! Number of sub-elements in 1D
    n_subelems_1D = nodesperedge-1

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the 2D vector field'
    message = message // var_name // ' to the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name = get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
  
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write 3D vector field
    write(i_unit,*)'<DataArray type="Float64" Name="',var_name,'" NumberOfComponents="3"  format="ascii">'

    if(ndim == 2) then
      do i_elem = 1, n_elems
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the four nodes that form the quadrilateral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by four nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1)

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1)
            
            local_node_3_idx = i_subelem + nodesperedge + 1 & 
              & + nodesperedge*(i_row_face-1)
            
            local_node_4_idx = i_subelem + nodesperedge &
              & + nodesperedge*(i_row_face-1)

            ! Write 3D vector field of the four nodes
            write(i_unit,*) var(1,local_node_1_idx,i_elem), &
              & var(2,local_node_1_idx,i_elem), var(3,local_node_1_idx,i_elem)

            write(i_unit,*) var(1,local_node_2_idx,i_elem), &
              & var(2,local_node_2_idx,i_elem), var(3,local_node_2_idx,i_elem)

            write(i_unit,*) var(1,local_node_3_idx,i_elem), &
              & var(2,local_node_3_idx,i_elem), var(3,local_node_3_idx,i_elem)

            write(i_unit,*) var(1,local_node_4_idx,i_elem), &
              & var(2,local_node_4_idx,i_elem), var(3,local_node_4_idx,i_elem)

          enddo
        enddo
      enddo

    else

      do i_elem = 1, n_elems
        do i_slice = 1, n_subelems_1D
          do i_row_face = 1, n_subelems_1D
            do i_subelem = 1, n_subelems_1d

              ! Construct indices of the eight nodes that form the hexahedral
              ! This algorithm should also work when we will have curved 
              ! elements because each subelement is made up by eight nodes
              local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1) 
              local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1)
              local_node_3_idx = i_subelem + nodesperedge + 1 + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
              local_node_4_idx = i_subelem + nodesperedge + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
            
              local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice
              local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice 
              local_node_7_idx = i_subelem + nodesperedge + 1 &
                & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 
              local_node_8_idx = i_subelem + nodesperedge &
                & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice


              ! Write 3D vector field of the eight nodes
              write(i_unit,*) var(1,local_node_1_idx,i_elem), & 
                & var(2,local_node_1_idx,i_elem), &
                & var(3,local_node_1_idx,i_elem)

              write(i_unit,*) var(1,local_node_2_idx,i_elem), &
                & var(2,local_node_2_idx,i_elem), &
                & var(3,local_node_2_idx,i_elem)

              write(i_unit,*) var(1,local_node_3_idx,i_elem), &
                & var(2,local_node_3_idx,i_elem), &
                & var(3,local_node_3_idx,i_elem)

              write(i_unit,*) var(1,local_node_4_idx,i_elem), &
                & var(2,local_node_4_idx,i_elem), &
                & var(3,local_node_4_idx,i_elem)
            
              write(i_unit,*) var(1,local_node_5_idx,i_elem), &
                & var(2,local_node_5_idx,i_elem), &
                & var(3,local_node_5_idx,i_elem)

              write(i_unit,*) var(1,local_node_6_idx,i_elem), &
                & var(2,local_node_6_idx,i_elem), &
                & var(3,local_node_6_idx,i_elem)

              write(i_unit,*) var(1,local_node_7_idx,i_elem), &
                & var(2,local_node_7_idx,i_elem), &
                & var(3,local_node_7_idx,i_elem)

              write(i_unit,*) var(1,local_node_8_idx,i_elem), &
                & var(2,local_node_8_idx,i_elem), & 
                & var(3,local_node_8_idx,i_elem)

            enddo
          enddo
        enddo
      enddo

    endif

    write(i_unit,*)'</DataArray>'
    
    ! Close file unit
    close(i_unit) 
   
    return
  end subroutine write_vector_3d_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_nodes_discontinous_elements_ascii_vtu_file - Writes the nodes 
  ! coordinates of the nodes to the .vtu file. ASCII file.

  subroutine write_nodes_discontinous_elements_ascii_vtu_file()
    
    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO
    use variables, only : xg

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_subelems_1d
    character(120) :: message
    integer :: i_elem, i_subelem, i_row_face, i_slice
    integer :: local_node_1_idx, local_node_2_idx, local_node_3_idx, &
               & local_node_4_idx
    integer :: local_node_5_idx, local_node_6_idx, local_node_7_idx, &
               & local_node_8_idx


    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of sub-elements
    n_subelems_1D = nodesperedge-1

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the coordinates of &
      & the nodes to the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name = get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
    
    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
  
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write coordinates
    write(i_unit,*)'<Points>'

    write(i_unit,*)'<DataArray type="Float64"  NumberOfComponents="3"  format="ascii">'

    if(ndim == 2) then
      
      do i_elem = iell, elem_high
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the four nodes that form the quadrilateral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by four nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) 
            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) 
            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) 
            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1)

            ! Write coordinates of the four nodes
            write(i_unit,1003) xg(1,local_node_1_idx,i_elem), &
              & xg(2,local_node_1_idx,i_elem), xg(3,local_node_1_idx,i_elem)

            write(i_unit,1003) xg(1,local_node_2_idx,i_elem), &
              & xg(2,local_node_2_idx,i_elem), xg(3,local_node_2_idx,i_elem)

            write(i_unit,1003) xg(1,local_node_3_idx,i_elem), &
              & xg(2,local_node_3_idx,i_elem), xg(3,local_node_3_idx,i_elem)

            write(i_unit,1003) xg(1,local_node_4_idx,i_elem), &
              & xg(2,local_node_4_idx,i_elem), xg(3,local_node_4_idx,i_elem)

          enddo
        enddo
      enddo

    else

      do i_elem = iell, elem_high
        do i_slice = 1, n_subelems_1D
          do i_row_face = 1, n_subelems_1D
            do i_subelem = 1, n_subelems_1d

              ! Construct indices of the eight nodes that form the hexahedral
              ! This algorithm should also work when we will have curved 
              ! elements because each subelement is made up by eight nodes
              local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1) 

              local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * (i_slice-1)

              local_node_3_idx = i_subelem + nodesperedge + 1 + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

              local_node_4_idx = i_subelem + nodesperedge + &
                & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
            
              local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice

              local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
                & + nodesperedge**2 * i_slice 

              local_node_7_idx = i_subelem + nodesperedge + 1 &
                & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

              local_node_8_idx = i_subelem + nodesperedge &
                & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice


              ! Write coordinates of the eight nodes
              write(i_unit,1003) xg(1,local_node_1_idx,i_elem), &
                & xg(2,local_node_1_idx,i_elem), xg(3,local_node_1_idx,i_elem)

              write(i_unit,1003) xg(1,local_node_2_idx,i_elem), &
                & xg(2,local_node_2_idx,i_elem), xg(3,local_node_2_idx,i_elem)

              write(i_unit,1003) xg(1,local_node_3_idx,i_elem), &
                & xg(2,local_node_3_idx,i_elem), xg(3,local_node_3_idx,i_elem)

              write(i_unit,1003) xg(1,local_node_4_idx,i_elem), &
                & xg(2,local_node_4_idx,i_elem), xg(3,local_node_4_idx,i_elem)
            
              write(i_unit,1003) xg(1,local_node_5_idx,i_elem), &
                & xg(2,local_node_5_idx,i_elem), xg(3,local_node_5_idx,i_elem)

              write(i_unit,1003) xg(1,local_node_6_idx,i_elem), &
                & xg(2,local_node_6_idx,i_elem), xg(3,local_node_6_idx,i_elem)

              write(i_unit,1003) xg(1,local_node_7_idx,i_elem), &
                & xg(2,local_node_7_idx,i_elem), xg(3,local_node_7_idx,i_elem)

              write(i_unit,1003) xg(1,local_node_8_idx,i_elem), &
                & xg(2,local_node_8_idx,i_elem), xg(3,local_node_8_idx,i_elem)

            enddo
          enddo
        enddo
      enddo

    endif

    1003 format(3(f17.10,1x))

    write(i_unit,*)
    write(i_unit,*)'</DataArray>' 
    write(i_unit,*)'</Points>' 

   
    ! Close file unit
    close(i_unit)

    return
  end subroutine write_nodes_discontinous_elements_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_connectivity_discontinous_elements_ascii_vtu_file - Writes the cell 
  ! connectivity to the .vtu file. ASCII file.

  subroutine write_connectivity_discontinous_elements_ascii_vtu_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_elems, n_subelems_1d, n_subelems
    character(120) :: message
    integer :: i_elem, i_subelem, i_row_face, i_slice
    integer :: global_node_1_idx, global_node_2_idx, global_node_3_idx, &
               & global_node_4_idx
    integer :: global_node_5_idx, global_node_6_idx, global_node_7_idx, &
               & global_node_8_idx
    integer :: elem_type
    integer :: shift

    ! Define cell typ according to the paraview data format
    if(ndim==2) then
      elem_type = 9 ! Quadrilateral element
    else
      elem_type = 12 ! Hexahedral element
    endif 

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high+1-iell
    
    ! Number of sub-elements in 1D
    n_subelems_1D = nodesperedge-1

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the cell connectivity &
      & to the .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name =  get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='old',access='append', &
      & iostat=io_status)
  
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write connectivity to file
    write(i_unit,*)'<Cells>'

    write(i_unit,*)'<DataArray type="Int32" Name="connectivity" format="ascii">'

    ! Set to zero shift variable
    shift = 0
    if(ndim == 2) then
      do i_elem = 1, n_elems
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct global indices (global in the ParaView context) for 
            ! connectivity 
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by four nodes
            global_node_1_idx = 0 + shift
            global_node_2_idx = 1 + shift
            global_node_3_idx = 2 + shift
            global_node_4_idx = 3 + shift

            ! Write the ID of the four nodes
            write(i_unit,1004) global_node_1_idx, global_node_2_idx, &
              & global_node_3_idx, global_node_4_idx

            shift = shift + 4
  
          enddo
        enddo
      enddo

      write(i_unit,*)'</DataArray>' 
      write(i_unit,*)'<DataArray type="Int32"  Name="offsets"  format="ascii">'

      ! Number of sub-elements
      n_subelems = n_elems*(nodesperedge-1)**ndim

      shift = 4
      do i_subelem = 1, n_subelems 
        write(i_unit,1004) shift
        shift = shift + 4
      enddo

      write(i_unit,*)'</DataArray>' 
      write(i_unit,*)'<DataArray type="UInt8"  Name="types"  format="ascii">'
      write(i_unit,1004) (elem_type,i_elem = 1,n_subelems) 
      write(i_unit,*)'</DataArray>' 
      write(i_unit,*)'</Cells>'

    else

      do i_elem = 1, n_elems
        do i_slice = 1, n_subelems_1d
          do i_row_face = 1, n_subelems_1D
            do i_subelem = 1, n_subelems_1d

              ! Construct global indices (global in the ParaView context) for 
              ! connectivity 
              ! This algorithm should also work when we will have curved 
              ! elements because each subelement is made up by eight nodes
              global_node_1_idx = 0 + shift
              global_node_2_idx = 1 + shift
              global_node_3_idx = 2 + shift
              global_node_4_idx = 3 + shift
              global_node_5_idx = 4 + shift
              global_node_6_idx = 5 + shift
              global_node_7_idx = 6 + shift
              global_node_8_idx = 7 + shift

              ! Write the ID of the eight nodes
              write(i_unit,1004) global_node_1_idx, global_node_2_idx, &
                & global_node_3_idx, global_node_4_idx, global_node_5_idx, &
                & global_node_6_idx, global_node_7_idx, global_node_8_idx

              shift = shift + 8
  
            enddo
          enddo
        enddo
      enddo

      write(i_unit,*)'</DataArray>' 
      write(i_unit,*)'<DataArray type="Int32"  Name="offsets"  format="ascii">'

      ! Number of sub-elements
      n_subelems = n_elems*(nodesperedge-1)**ndim

      shift = 8
      do i_subelem = 1, n_subelems 
        write(i_unit,1004) shift
        shift = shift + 8
      enddo

      write(i_unit,*)'</DataArray>' 
      write(i_unit,*)'<DataArray type="UInt8"  Name="types"  format="ascii">'
      write(i_unit,1004) (elem_type,i_elem = 1,n_subelems) 
      write(i_unit,*)'</DataArray>' 
      write(i_unit,*)'</Cells>'

    endif

    1004 format(1x, 10i8)

    ! Close file unit
    close(i_unit) 

    return
   end subroutine write_connectivity_discontinous_elements_ascii_vtu_file

  !============================================================================

  !============================================================================
  ! write_parallel_vtu_file - Writes the .pvtu file for launching ParaView

  subroutine write_parallel_vtu_file()
    
    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO

     ! Nothing is implicitly defined
    implicit none

    character(120) :: pvtu_file_name
    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    character(120) :: tag_iter
    integer :: i_unit, io_status
    integer :: max_unit
    character(120) :: message
    integer :: proc_ID

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the .pvtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Get iteration tag for file name
    tag_iter = get_tag_iter(itimestep)

    ! Create .pvtu file name
    ! This is commented out because it would requires renaming the *.pvtu for
    ! creating movies
!    pvtu_file_name = get_solution_pvtu_file_name(tag_proc,tag_time)

    pvtu_file_name = get_solution_pvtu_file_name(tag_proc,tag_time,tag_iter)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
  
    ! Open file
    open(unit=i_unit,file=pvtu_file_name,status='unknown',iostat=io_status)

    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    ! Write XML string to the .pvtu file
    write(i_unit,'(A)') '<?xml version="1.0"?>'
  
    write(i_unit,'(A)') '<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian">'

    ! Overlapping between each file. 
    ! GhostLevel="0" means no overlapping
    write(i_unit,*)'<PUnstructuredGrid GhostLevel="0">'


    write(i_unit,*)'<PPointData>'

    ! Write ASCII .pvtu file
    if (write_solution_formatted) then
    
      ! Conserved variables
      write(i_unit,*)'<PDataArray type="Float64" Name="Density" format="ascii"/>'
!     write(i_unit,*)'<PDataArray type="Float64" Name="Density Error" format="ascii"/>'
      write(i_unit,*)'<PDataArray type="Float64" Name="temperature" format="ascii"/>'
!     write(i_unit,*)'<PDataArray type="Float64" Name="temperature Error" format="ascii"/>'
      write(i_unit,*)'<DataArray  type="Float64" Name="momentum" NumberOfComponents="3" format="ascii"/>'
      write(i_unit,*)'<PDataArray type="Float64" Name="specific total energy" format="ascii"/>'

      ! Velocity vector
      write(i_unit,*)'<DataArray  type="Float64" Name="velocity" NumberOfComponents="3" format="ascii"/>'
!     write(i_unit,*)'<DataArray  type="Float64" Name="velocity Error" NumberOfComponents="3" format="ascii"/>'

      ! Vorticity vector
      write(i_unit,*)'<DataArray  type="Float64" Name="vorticity" NumberOfComponents="3" format="ascii"/>'
    
      if (time_averaging) then
      
        ! Mean primitive variables
        write(i_unit,*)'<PDataArray type="Float64" Name="mean density" format="ascii"/>'
        write(i_unit,*)'<DataArray  type="Float64" Name="mean velocity" NumberOfComponents="3" format="ascii"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="mean temperature" format="ascii"/>'

        ! Reynolds stress
        write(i_unit,*)'<PDataArray type="Float64" Name="up up" format="ascii"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="up vp" format="ascii"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="up wp" format="ascii"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="vp vp" format="ascii"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="vp wp" format="ascii"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="wp wp" format="ascii"/>'

      endif

      ! Specific entropy
      write(i_unit,*)'<PDataArray type="Float64" Name="specific entropy" format="ascii"/>'

      ! Nodes
      write(i_unit,*)'</PPointData>'
      write(i_unit,*)'<PPoints>' 
      write(i_unit,*)'<PDataArray type="Float64"  NumberOfComponents="3"  format="ascii"/>'
      write(i_unit,*)'</PPoints>'

    ! Write BINARY .pvtu file  
    else
  
      ! Conserved variables
      write(i_unit,*)'<PDataArray type="Float64" Name="density" format="binary"/>'
      write(i_unit,*)'<PDataArray type="Float64" Name="temperature" format="binary"/>'
      write(i_unit,*)'<DataArray  type="Float64" Name="momentum" NumberOfComponents="3" format="binary"/>'
      write(i_unit,*)'<PDataArray type="Float64" Name="specific total energy" format="binary"/>'

      ! Velocity vector
      write(i_unit,*)'<DataArray  type="Float64" Name="velocity" NumberOfComponents="3" format="binary"/>'

      ! Vorticity vector
      write(i_unit,*)'<DataArray  type="Float64" Name="vorticity" NumberOfComponents="3" format="binary"/>'
    
      if (time_averaging) then
      
        ! Mean primitive variables
        write(i_unit,*)'<PDataArray type="Float64" Name="mean density" format="binary"/>'
        write(i_unit,*)'<DataArray  type="Float64" Name="mean velocity" NumberOfComponents="3" format="binary"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="mean temperature" format="binary"/>'

        ! Reynolds stress
        write(i_unit,*)'<PDataArray type="Float64" Name="up up" format="binary"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="up vp" format="binary"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="up wp" format="binary"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="vp vp" format="binary"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="vp wp" format="binary"/>'
        write(i_unit,*)'<PDataArray type="Float64" Name="wp wp" format="binary"/>'

      endif

      ! Specific entropy
      write(i_unit,*)'<PDataArray type="Float64" Name="specific entropy" format="binary"/>'

      ! Nodes
      write(i_unit,*)'</PPointData>'
      write(i_unit,*)'<PPoints>' 
      write(i_unit,*)'<PDataArray type="Float64"  NumberOfComponents="3"  format="binary"/>'
      write(i_unit,*)'</PPoints>'

    endif


    ! Write name of the .vtu that conatins the solution of each partition
    do proc_id = 0, nprocs-1

      ! Get processor tag for file name
      tag_proc = get_tag_proc(proc_id)

      ! Get time tag for file name
      tag_time = get_tag_time(timeglobal)

      ! Create the solution file name that has to be written by each processor
      vtu_file_name =  get_solution_proc_vtu_file_name(tag_proc,tag_time)
      
      ! ATTENTION: in order to open the files containing the solution, the
      ! .pvtu needs to have a path as ../write_solution_dir/vtu_file_name.
      ! Thus, we have to add a dot (.) in fornt of
      ! ./write_solution_dir/vtu_file_name
      write(i_unit,*)'  <Piece Source=".',trim(vtu_file_name),'"/>'

    enddo

    write(i_unit,*)' </PUnstructuredGrid>'
    write(i_unit,*)'</VTKFile>'
    
    ! Close file unit
    close(i_unit)  

  end subroutine write_parallel_vtu_file

  !============================================================================
  
  !============================================================================

  subroutine write_xml_raw_binary_vtu_file()

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables
    use tools_IO
    use binary_sizes
    use variables
      
    ! Nothing is implicitly defined
    implicit none
    
    character(120) :: vtu_file_name
    character(120) :: tag_proc
    character(120) :: tag_time
    integer :: i_unit, io_status
    integer :: max_unit
    integer :: iell, elem_high
    integer :: n_elems, n_subelems, n_nodes, n_vertices_sub_elem, n_subelems_1d
    character(120) :: message
    character(8) :: str_1, str_2, offset
    character(200) :: buffer
    character(1) :: lf
    integer :: ioff
    real(wp) :: tmp

    integer :: i_elem, i_row_face, i_subelem, i_slice
    integer :: local_node_1_idx, local_node_2_idx, local_node_3_idx, &
               & local_node_4_idx
    integer :: local_node_5_idx, local_node_6_idx, local_node_7_idx, &
               & local_node_8_idx
    
    integer :: global_node_1_idx, global_node_2_idx, global_node_3_idx, &
               & global_node_4_idx
    integer :: global_node_5_idx, global_node_6_idx, global_node_7_idx, &
               & global_node_8_idx

    integer :: elem_type
    integer :: shift

    continue 

    ! Define cell typ according to the paraview data format
    if(ndim == 2) then
      elem_type = 9 ! Quadrilateral element
    else
      elem_type = 12  ! Hexahedral element
    endif 


    ! Line feed character
    lf = char(10)

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the XML part of the &
      &  raw binary .vtu file.'

    ! Set max number of unit
    max_unit = 98

    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    ! Get time tag for file name
    tag_time = get_tag_time(timeglobal)

    ! Create the solution file name that has to be written by each processor
    vtu_file_name = get_solution_proc_vtu_file_name(tag_proc,tag_time)

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)
    
    ! Open file
    open(unit=i_unit,file=vtu_file_name,status='unknown',form='unformatted',access='stream',iostat=io_status)
    
    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)
    
    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements
    n_elems = elem_high + 1 - iell
    
    ! Number of sub-elements
    n_subelems = n_elems*(nodesperedge-1)**ndim

    ! Number of sub-elements in 1D
    n_subelems_1D = nodesperedge - 1

    ! Number of nodes or vertices necessary to write the solution has
    ! discontinuous field also within an element
    if(ndim == 2) then
      n_vertices_sub_elem = 4
      n_nodes = n_subelems*n_vertices_sub_elem ! For quadrilateral elements (2D)
    else
      n_vertices_sub_elem = 8
      n_nodes = n_subelems*n_vertices_sub_elem ! For hexahedral elements (3D)
    endif

    ! Compute byte lengths
    n_bytes_scal   = n_nodes*sizeof(tmp)
    n_bytes_vec    = ndim*n_nodes*sizeof(tmp)
    n_bytes_xyz    = ndim*n_nodes*sizeof(tmp)
    n_bytes_ien    = n_nodes*sizeof(n_elems)
    n_bytes_offset = n_subelems*sizeof(n_elems)
    n_bytes_etype  = n_subelems*sizeof(n_elems)

  
    ! Write header file
    ! Use version="1.0" for XML. One could also use version 0.1
    buffer = '<?xml version="1.0"?>'//lf
    write(i_unit) trim(buffer)

    buffer = '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">'//lf
    write(i_unit) trim(buffer)
   
    buffer = '  <UnstructuredGrid>'//lf
    write(i_unit) trim(buffer)

    write(str_1(1:8),'(i8)') n_nodes
    write(str_2(1:8),'(i8)') n_subelems

    buffer = '<Piece NumberOfPoints="'//str_1//'" NumberOfCells="'//str_2//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '<PointData> '//lf
    write(i_unit) trim(buffer)

    ! Density XML header    
    ! ==================
    ioff = 0
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Density" NumberOfComponents="1" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Temperature XML header
    ! ======================
    ioff = ioff + sizeof(n_elems) + n_bytes_scal
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Temperature" NumberOfComponents="1" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Momentum XML header
    ! ===================
    ioff = ioff + sizeof(n_elems) + n_bytes_scal
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Momentum" NumberOfComponents="3" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Specific total energy XML header
    ! ================================
    ioff = ioff + sizeof(n_elems) + n_bytes_vec
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Specific total energy" NumberOfComponents="1" &
      & format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Velocity XML header
    ! ===================
    ioff = ioff + sizeof(n_elems) + n_bytes_scal
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Vorticity XML header
    ! ====================
    ioff = ioff + sizeof(n_elems) + n_bytes_vec
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Vorticity" NumberOfComponents="3" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Specific entropy XML header
    ! ===========================
    ioff = ioff + sizeof(n_elems) + n_bytes_vec
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64" Name="Specific entropy" NumberOfComponents="1" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)


    buffer = '</PointData>'//lf
    write(i_unit) trim(buffer)

    buffer = '<Points>'//lf
    write(i_unit) trim(buffer)

    ! Nodes coordinates XML header
    ! ============================
    ioff = ioff + sizeof(n_elems) + n_bytes_scal
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Float64"  NumberOfComponents="3"  format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)
    
    buffer = '</Points>'//lf 
    write(i_unit) trim(buffer)

    ! Connectivity XML header
    ! =======================
    buffer = '<Cells>'//lf
    write(i_unit) trim(buffer)

    ioff = ioff + sizeof(n_elems) + n_bytes_xyz
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

    ! Offset XML header
    ! =================
    ioff = ioff + sizeof(n_elems) + n_bytes_ien
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="Int32"  Name="offsets"  format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)

    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)

  
    ! Elements type XML header
    ! ========================
    ioff = ioff + sizeof(n_elems) + n_bytes_offset
    write(offset(1:8),'(i8)') ioff
    buffer = '<DataArray type="UInt32"  Name="types"  format="appended" offset="'//offset//'">'//lf
    write(i_unit) trim(buffer)    
    
    buffer = '</DataArray>'//lf
    write(i_unit) trim(buffer)    

    buffer = '</Cells>'//lf
    write(i_unit) trim(buffer)  

    buffer = '</Piece>'//lf
    write(i_unit) trim(buffer)  
    
    buffer = '</UnstructuredGrid>'//lf
    write(i_unit) trim(buffer)  
    
    buffer = '<AppendedData encoding="raw">'//lf
    write(i_unit) trim(buffer)  
    
    buffer = '_'
    write(i_unit) trim(buffer)  
   
    ! Append binary data
    ! ==================
    ! Write density field
    write(i_unit) n_bytes_scal
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) ug(1,local_node_1_idx,i_elem)
            write(i_unit) ug(1,local_node_2_idx,i_elem)
            write(i_unit) ug(1,local_node_3_idx,i_elem)
            write(i_unit) ug(1,local_node_4_idx,i_elem)
        
            write(i_unit) ug(1,local_node_5_idx,i_elem)
            write(i_unit) ug(1,local_node_6_idx,i_elem)
            write(i_unit) ug(1,local_node_7_idx,i_elem)
            write(i_unit) ug(1,local_node_8_idx,i_elem)

          end do
        end do
      end do
    end do

    ! Write temperature field
    write(i_unit) n_bytes_scal
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) vg(5,local_node_1_idx,i_elem)
            write(i_unit) vg(5,local_node_2_idx,i_elem)
            write(i_unit) vg(5,local_node_3_idx,i_elem)
            write(i_unit) vg(5,local_node_4_idx,i_elem)
        
            write(i_unit) vg(5,local_node_5_idx,i_elem)
            write(i_unit) vg(5,local_node_6_idx,i_elem)
            write(i_unit) vg(5,local_node_7_idx,i_elem)
            write(i_unit) vg(5,local_node_8_idx,i_elem)

          end do
        end do
      end do
    end do

    ! Write momentum field
    write(i_unit) n_bytes_vec
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) ug(2,local_node_1_idx,i_elem), &
                        & ug(3,local_node_1_idx,i_elem), &
                        & ug(4,local_node_1_idx,i_elem)

            write(i_unit) ug(2,local_node_2_idx,i_elem), &
                        & ug(3,local_node_2_idx,i_elem), &
                        & ug(4,local_node_2_idx,i_elem) 

            write(i_unit) ug(2,local_node_3_idx,i_elem), &
                        & ug(3,local_node_3_idx,i_elem), &
                        & ug(4,local_node_3_idx,i_elem)

            write(i_unit) ug(2,local_node_4_idx,i_elem), &
                        & ug(3,local_node_4_idx,i_elem), &
                        & ug(4,local_node_4_idx,i_elem)
        
            write(i_unit) ug(2,local_node_5_idx,i_elem), &
                        & ug(3,local_node_5_idx,i_elem), &
                        & ug(4,local_node_5_idx,i_elem)

            write(i_unit) ug(2,local_node_6_idx,i_elem), &
                        & ug(3,local_node_6_idx,i_elem), &
                        & ug(4,local_node_6_idx,i_elem)

            write(i_unit) ug(2,local_node_7_idx,i_elem), &
                        & ug(3,local_node_7_idx,i_elem), &
                        & ug(4,local_node_7_idx,i_elem)

            write(i_unit) ug(2,local_node_8_idx,i_elem), &
                        & ug(3,local_node_8_idx,i_elem), &
                        & ug(4,local_node_8_idx,i_elem)
          end do
        end do
      end do
    end do

    ! Write total specific energy field
    write(i_unit) n_bytes_scal
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) ug(5,local_node_1_idx,i_elem)
            write(i_unit) ug(5,local_node_2_idx,i_elem)
            write(i_unit) ug(5,local_node_3_idx,i_elem)
            write(i_unit) ug(5,local_node_4_idx,i_elem)
        
            write(i_unit) ug(5,local_node_5_idx,i_elem)
            write(i_unit) ug(5,local_node_6_idx,i_elem)
            write(i_unit) ug(5,local_node_7_idx,i_elem)
            write(i_unit) ug(5,local_node_8_idx,i_elem)

          end do
        end do
      end do
    end do

    ! Write velocity field
    write(i_unit) n_bytes_vec
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) vg(2,local_node_1_idx,i_elem), &
                        & vg(3,local_node_1_idx,i_elem), &
                        & vg(4,local_node_1_idx,i_elem)

            write(i_unit) vg(2,local_node_2_idx,i_elem), &
                        & vg(3,local_node_2_idx,i_elem), &
                        & vg(4,local_node_2_idx,i_elem) 

            write(i_unit) vg(2,local_node_3_idx,i_elem), &
                        & vg(3,local_node_3_idx,i_elem), &
                        & vg(4,local_node_3_idx,i_elem)

            write(i_unit) vg(2,local_node_4_idx,i_elem), &
                        & vg(3,local_node_4_idx,i_elem), &
                        & vg(4,local_node_4_idx,i_elem)
        
            write(i_unit) vg(2,local_node_5_idx,i_elem), &
                        & vg(3,local_node_5_idx,i_elem), &
                        & vg(4,local_node_5_idx,i_elem)

            write(i_unit) vg(2,local_node_6_idx,i_elem), &
                        & vg(3,local_node_6_idx,i_elem), &
                        & vg(4,local_node_6_idx,i_elem)

            write(i_unit) vg(2,local_node_7_idx,i_elem), &
                        & vg(3,local_node_7_idx,i_elem), &
                        & vg(4,local_node_7_idx,i_elem)

            write(i_unit) vg(2,local_node_8_idx,i_elem), &
                        & vg(3,local_node_8_idx,i_elem), &
                        & vg(4,local_node_8_idx,i_elem)

          end do
        end do
      end do
    end do

    ! Write vorticity field
    write(i_unit) n_bytes_vec
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) omega(1,local_node_1_idx,i_elem), &
                        & omega(2,local_node_1_idx,i_elem), &
                        & omega(3,local_node_1_idx,i_elem)

            write(i_unit) omega(1,local_node_2_idx,i_elem), &
                        & omega(2,local_node_2_idx,i_elem), &
                        & omega(3,local_node_2_idx,i_elem) 

            write(i_unit) omega(1,local_node_3_idx,i_elem), &
                        & omega(2,local_node_3_idx,i_elem), &
                        & omega(3,local_node_3_idx,i_elem)

            write(i_unit) omega(1,local_node_4_idx,i_elem), &
                        & omega(2,local_node_4_idx,i_elem), &
                        & omega(3,local_node_4_idx,i_elem)
        
            write(i_unit) omega(1,local_node_5_idx,i_elem), &
                        & omega(2,local_node_5_idx,i_elem), &
                        & omega(3,local_node_5_idx,i_elem)

            write(i_unit) omega(1,local_node_6_idx,i_elem), &
                        & omega(2,local_node_6_idx,i_elem), &
                        & omega(3,local_node_6_idx,i_elem)

            write(i_unit) omega(1,local_node_7_idx,i_elem), &
                        & omega(2,local_node_7_idx,i_elem), &
                        & omega(3,local_node_7_idx,i_elem)

            write(i_unit) omega(1,local_node_8_idx,i_elem), &
                        & omega(2,local_node_8_idx,i_elem), &
                        & omega(3,local_node_8_idx,i_elem)

          end do
        end do
      end do
    end do

    ! Write specific entropy field
    write(i_unit) n_bytes_scal
    do i_elem = iell, elem_high 
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
        
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice 

            local_node_8_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write scalar field of the eight nodes
            write(i_unit) specific_entropy(local_node_1_idx,i_elem)
            write(i_unit) specific_entropy(local_node_2_idx,i_elem)
            write(i_unit) specific_entropy(local_node_3_idx,i_elem)
            write(i_unit) specific_entropy(local_node_4_idx,i_elem)
        
            write(i_unit) specific_entropy(local_node_5_idx,i_elem)
            write(i_unit) specific_entropy(local_node_6_idx,i_elem)
            write(i_unit) specific_entropy(local_node_7_idx,i_elem)
            write(i_unit) specific_entropy(local_node_8_idx,i_elem)

          end do
        end do
      end do
    end do

    ! Write coordinates
    write(i_unit) n_bytes_xyz
    do i_elem = iell, elem_high
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct indices of the eight nodes that form the hexahedral
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            local_node_1_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1) 

            local_node_2_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * (i_slice-1)

            local_node_3_idx = i_subelem + nodesperedge + 1 + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)

            local_node_4_idx = i_subelem + nodesperedge + &
              & nodesperedge*(i_row_face-1) + nodesperedge**2 * (i_slice-1)
          
            local_node_5_idx = i_subelem + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice

            local_node_6_idx = i_subelem + 1 + nodesperedge*(i_row_face-1) &
              & + nodesperedge**2 * i_slice 

            local_node_7_idx = i_subelem + nodesperedge + 1 &
              & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            local_node_8_idx = i_subelem + nodesperedge &
              & + nodesperedge*(i_row_face-1) + nodesperedge**2 * i_slice

            ! Write coordinates of the eight nodes
            write(i_unit) xg(1,local_node_1_idx,i_elem)
            write(i_unit) xg(2,local_node_1_idx,i_elem)
            write(i_unit) xg(3,local_node_1_idx,i_elem)

            write(i_unit) xg(1,local_node_2_idx,i_elem)
            write(i_unit) xg(2,local_node_2_idx,i_elem)
            write(i_unit) xg(3,local_node_2_idx,i_elem)

            write(i_unit) xg(1,local_node_3_idx,i_elem)
            write(i_unit) xg(2,local_node_3_idx,i_elem)
            write(i_unit) xg(3,local_node_3_idx,i_elem)

            write(i_unit) xg(1,local_node_4_idx,i_elem)
            write(i_unit) xg(2,local_node_4_idx,i_elem)
            write(i_unit) xg(3,local_node_4_idx,i_elem)
          
            write(i_unit) xg(1,local_node_5_idx,i_elem)
            write(i_unit) xg(2,local_node_5_idx,i_elem)
            write(i_unit) xg(3,local_node_5_idx,i_elem)

            write(i_unit) xg(1,local_node_6_idx,i_elem)
            write(i_unit) xg(2,local_node_6_idx,i_elem)
            write(i_unit) xg(3,local_node_6_idx,i_elem)

            write(i_unit) xg(1,local_node_7_idx,i_elem)
            write(i_unit) xg(2,local_node_7_idx,i_elem)
            write(i_unit) xg(3,local_node_7_idx,i_elem)

            write(i_unit) xg(1,local_node_8_idx,i_elem)
            write(i_unit) xg(2,local_node_8_idx,i_elem)
            write(i_unit) xg(3,local_node_8_idx,i_elem)

          enddo
        enddo
      enddo
    enddo

    ! Write connectivity
    shift = 0
    write(i_unit) n_bytes_ien
    do i_elem = 1, n_elems
      do i_slice = 1, n_subelems_1d
        do i_row_face = 1, n_subelems_1d
          do i_subelem = 1, n_subelems_1d

            ! Construct global indices (global in the ParaView context) for 
            ! connectivity 
            ! This algorithm should also work when we will have curved 
            ! elements because each subelement is made up by eight nodes
            global_node_1_idx = 0 + shift
            global_node_2_idx = 1 + shift
            global_node_3_idx = 2 + shift
            global_node_4_idx = 3 + shift
            global_node_5_idx = 4 + shift
            global_node_6_idx = 5 + shift
            global_node_7_idx = 6 + shift
            global_node_8_idx = 7 + shift

            ! Write the ID of the eight nodes
            write(i_unit) global_node_1_idx
            write(i_unit) global_node_2_idx
            write(i_unit) global_node_3_idx
            write(i_unit) global_node_4_idx
            write(i_unit) global_node_5_idx
            write(i_unit) global_node_6_idx
            write(i_unit) global_node_7_idx
            write(i_unit) global_node_8_idx

            shift = shift + 8

          enddo
        enddo
      enddo
    enddo

    ! Write shift
    write(i_unit) n_bytes_offset
    shift = 8
    do i_subelem = 1, n_subelems 
      write(i_unit) shift
      shift = shift + 8
    enddo
    
    ! Write element type
    write(i_unit) n_bytes_etype
    do i_elem = 1, n_subelems
      write(i_unit) elem_type
    end do

    buffer = lf//'</AppendedData>'//lf
    write(i_unit) trim(buffer)

    buffer = '</VTKFile>'//lf
    write(i_unit) trim(buffer)

    ! Close file
    close(i_unit)

    return
  end subroutine write_xml_raw_binary_vtu_file

  !============================================================================
  
  !============================================================================

  subroutine write_aerodynamic_coefficients(time_step_cnt,global_time)

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables, only : c_1_aero_coeff, c_2_aero_coeff, c_3_aero_coeff

    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time

    character(120) :: aero_coeff
    integer        :: i_unit, io_status
    integer        :: max_unit
    character(120) :: message
    
    continue

    aero_coeff = get_aerodynamic_coefficients_file_name()

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the aerodynamic coefficients.'

    ! Set max number of unit
    max_unit = 98

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=aero_coeff,status="unknown",position="append",action="write",iostat=io_status)

    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    write(i_unit,230) time_step_cnt, global_time, c_1_aero_coeff, c_2_aero_coeff, c_3_aero_coeff

    230 format(i8,1x,4(e17.10,1x))

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_aerodynamic_coefficients

  !============================================================================
 
  !============================================================================

  subroutine write_error_no_slip_wall_bc(time_step_cnt,global_time)

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables, only : l2_error_bc_no_slip_wall_u_1,   &
                               & l2_error_bc_no_slip_wall_u_2,   &
                               & l2_error_bc_no_slip_wall_u_3,   &
                               & linf_error_bc_no_slip_wall_u_1, &
                               & linf_error_bc_no_slip_wall_u_2, &
                               & linf_error_bc_no_slip_wall_u_3


    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time

    character(120) :: error_no_slip_wall_bc
    integer        :: i_unit, io_status
    integer        :: max_unit
    character(120) :: message
    
    continue

    error_no_slip_wall_bc = get_error_no_slip_wall_bc_file_name()

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the no-slip wall bc error.'

    ! Set max number of unit
    max_unit = 98

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=error_no_slip_wall_bc,status="unknown",position="append",action="write",iostat=io_status)

    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    write(i_unit,230) time_step_cnt, global_time,     &
                    & l2_error_bc_no_slip_wall_u_1,   &
                    & l2_error_bc_no_slip_wall_u_2,   &
                    & l2_error_bc_no_slip_wall_u_3,   &
                    & linf_error_bc_no_slip_wall_u_1, &
                    & linf_error_bc_no_slip_wall_u_2, &
                    & linf_error_bc_no_slip_wall_u_3

    230 format(i8,1x,7(e17.10,1x))

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_error_no_slip_wall_bc

  !============================================================================
  
  !============================================================================

  subroutine write_dkinetic_energy_dt(time_step_cnt,global_time)

    ! Load modules
    use mpimod
    use referencevariables
    use variables, only : dkinetic_energy_dt, kinetic_energy

    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time

    character(120) :: dke_dt
    integer        :: i_unit, io_status
    integer        :: max_unit
    character(120) :: message
    
    continue

    dke_dt = get_dkinetic_energy_dt_file_name()

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the dke/dt.'

    ! Set max number of unit
    max_unit = 98

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=dke_dt,status="unknown",position="append",action="write",iostat=io_status)

    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    write(i_unit,230) time_step_cnt, global_time, &
                    & dkinetic_energy_dt, kinetic_energy

    230 format(i8,1x,7(e17.10,1x))

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_dkinetic_energy_dt

  !============================================================================
  
  !============================================================================

  subroutine write_enstrophy(time_step_cnt,global_time)

    ! Load modules
    use mpimod
    use referencevariables
    use variables, only : enstrophy

    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time

    character(120) :: enstr
    integer        :: i_unit, io_status
    integer        :: max_unit
    character(120) :: message
    
    continue

    enstr = get_enstrophy_file_name()

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the enstrophy.'

    ! Set max number of unit
    max_unit = 98

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=enstr,status="unknown",position="append",action="write",iostat=io_status)

    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    write(i_unit,230) time_step_cnt, global_time, &
                    & enstrophy

    230 format(i8,1x,7(e17.10,1x))

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_enstrophy


  !============================================================================
 
  !============================================================================

  subroutine write_error_heat_entropy_flow_wall_bc(time_step_cnt,global_time)

    ! Load modules
    use mpimod
    use referencevariables
    use controlvariables, only : linf_error_heat_entropy_flow_wall_bc

    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time

    character(120) :: error_entropy_flow_wall_bc
    integer        :: i_unit, io_status
    integer        :: max_unit
    character(120) :: message
    
    continue

    error_entropy_flow_wall_bc = get_error_heat_entropy_flow_wall_bc_file_name()

    ! Define error message in case the open instruction returns an error
    message = 'Failure in opening the file for writing the heat entropy flow wall bc error.'

    ! Set max number of unit
    max_unit = 98

    ! Get free IO unit
    i_unit = get_file_unit(max_unit)

    ! Open file
    open(unit=i_unit,file=error_entropy_flow_wall_bc,status="unknown",position="append",action="write",iostat=io_status)

    ! Check io_status
    call check_io_open_file(io_status,myprocid,message)

    write(i_unit,230) time_step_cnt, global_time,     &
                    & linf_error_heat_entropy_flow_wall_bc

    230 format(i8,1x,2(e17.10,1x))

    ! Close file unit
    close(i_unit)

    return
  end subroutine write_error_heat_entropy_flow_wall_bc

  !============================================================================

  subroutine Sods_Line_Plot()


    use mpimod
    use referencevariables
    use variables,             only: xg, vg
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

    character(120)  :: output_file_name
    character(120)  :: tag_proc
    integer         :: i_elem, i_node

!   get_solution_proc_vtu_file_name = './' // trim(write_solution_dir) // '/' // trim(write_solution_common_name) &
!     & // '_p' // trim(tag_proc) // '_t'// trim(tag_time) // ".vtu"


    ! Get processor tag for file name
    tag_proc = get_tag_proc(myprocid)

    output_file_name = 'outrk6_p' // trim(tag_proc) // ".dat"
    write(*,*)'output_file_name = ',output_file_name
    open(unit=46,file=output_file_name)

      do i_elem = ihelems(1), ihelems(2) 
        do i_node = 1, nodesperelem

   !  X axis
          if(abs(dot_product(xg(2:3,i_node,i_elem),xg(2:3,i_node,i_elem))) <= 1.0e-12)then
            write(46,2) xg(1,i_node,i_elem),                        &
                        vg(1,i_node,i_elem),                        &
                        vg(2,i_node,i_elem),                        &
                        vg(1,i_node,i_elem)*vg(5,i_node,i_elem)
          endif

   !  X-Y 45
!         if((abs(xg(3,i_node,i_elem)-0.5_wp) <= 1.0e-10) .and.                          &          
!            (abs(xg(1,i_node,i_elem)-xg(2,i_node,i_elem)) <= 1.0e-10)  .and.            &
!            (sqrt(dot_product(xg(1:2,i_node,i_elem),xg(1:2,i_node,i_elem))) - 0.5_wp) <= 1.0e-10) then
!           write(46,2) sign(sqrt(dot_product(xg(1:2,i_node,i_elem),xg(1:2,i_node,i_elem))),xg(1,i_node,i_elem)) + 0.5_wp, &
!                       vg(1,i_node,i_elem),                                                     &
!                       sqrt(dot_product(vg(2:3,i_node,i_elem),vg(2:3,i_node,i_elem))),          &
!                       vg(1,i_node,i_elem)*vg(5,i_node,i_elem)
!         endif

        enddo
      enddo

    2 format(7(1x,e20.12))

      close(unit=46)

  end subroutine Sods_Line_Plot
  !============================================================================

end module  write_solution_file

