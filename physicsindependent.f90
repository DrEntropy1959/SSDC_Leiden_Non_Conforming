! This modules contains all the subroutine and functions which drive the
! preparation of the simulation independently of the physics. It can be
! considered as a pre-processing step. 

module physicsindependent

  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none

  ! Subroutines and functions in this module are usually private
  private

  ! Exceptions, i.e. public subroutines or functions
  public :: physics_independent_setup
  

contains

  !============================================================================
  
  !============================================================================
  ! physics_independent_setup - Sets up all the physics independent informations 
  ! needed to perform the simulation.

  subroutine physics_independent_setup()

    ! Load modules
    use polyinit
    use controlvariables
    use collocationvariables
    use referencevariables
    use initgrid
    use fileio
    use mpimod
    use write_solution_file
    use referencevariables,  only: ihelems, nfacesperelem, nelems
    use initcollocation,     only: element_properties

    ! Nothing is implicitly defined
    implicit none

    integer :: i_err

!-- DEBUG DAVID START
!   integer ielem, iface, i, n_pts_2d, inode, jnode
!   integer,  allocatable, dimension(:,:) :: kfacenodes_On
!   real(wp) :: temp(3,8)
!-- DEBUG DAVID END

    continue

    ! Initialize collocation approximation
!   i_err = rmapInit(npoly,ndim)

!   if (i_err > 0) then
!     write(*,*) 'An error occured in the function rmapInit!'
!     write(*,*) 'Exiting'
!     stop
!   end if

    ! Here we specify the qualities of each element type according to the CGNS
    ! standard. Such features are also used for all the other grid format
    ! supported in this code.
    call init_edge_2() 
    call init_quad_4() 
    call init_hex_8()  
    call init_elem_type()  

    ! Initialize the BC numbering of this code
    call init_code_bc_types()

    ! Read the grid and construct all the connectivity based on the input grid
    ! format
    if (grid_format == 'cgns') then

      ! Define boundary condition types
!      call init_cgns_bc()

      ! The master node compute the connectivity and calls metis to assign 
      ! the elements to each process.
      if (myprocid == 0) then
        ! Read only the necessary information from the datafile
!        call cgnsPrelimUnstructuredGrid(casefile)

        ! Check the number of processors and the number of elements
        call check_n_procs_n_elems(nprocs,nelems)

        ! Create connectivity from the original grid
        call E2EConnectivity_cgns() 
     
        ! Metis calculates the partitions
        call calculatepartitions() 
      
      end if

      ! Push element connectivity to all processes
      call distributeelements_cgns()
    
      ! Read vertices now from grid file using known connectivity information
!      call cgnsParallelUnstructuredGrid(casefile)
    
    
    else if (grid_format == 'aflr3') then

      ! The master node compute the connectivity and calls metis to assign 
      ! the elements to each process.
      if (myprocid == 0) then
         
        write(*,*) 'Master node reads the AFLR3 grid'
        write(*,*) '==============================================================='
        ! Read only the necessary information from the datafile
        call aflr3ReadUnstructuredGrid(casefile)
        ! Check the number of processors and the number of elements
        call check_n_procs_n_elems(nprocs,nelems) 

        write(*,*) 'Master node builds connectivity and orientation arrays'
!      write(*,*) '==============================================================='

        ! Create connectivity from the original grid
        call e2e_connectivity_aflr3()  
!-- DEBUG DAVID START
!do ielem = 1,size(ic2nh,2)
!write(*,*)'+++++++++++++++++++++++++++'
!  do iface = 1,6
!    temp = vx_Master(:,ic2nh(:,ielem))
!    write(*,*)'ielem = ',ielem,' iface = ',iface,&
!    NEW_LINE('A')//'Adjoining element face ID = ',ef2e(1,iface,ielem),&
!    NEW_LINE('A')//' Adjoining element ID = ', ef2e(2,iface,ielem)!,&
!    NEW_LINE('A')//' Adjoining element polynomial order = ',ef2e(4,iface,ielem),&
!    NEW_LINE('A')//' max r = ',maxval(abs(temp))
!  enddo
!write(*,*)'+++++++++++++++++++++++++++'
!enddo
!-- DEBUG DAVID END
        ! Establish face orientation between connected faces
        call face_orientation_aflr3()

!       write(*,*) 'Master node builds element orders if non-conforming'
!       write(*,*) '==============================================================='
        call set_element_orders_serial()    
  
        ! create e_edge2e connectivity
        call e_edge2e_connectivity()   
 !-- DAVID DEBUG START
!write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~'
!iface = 6
!ielem = 8
!  write(*,*)'original'//NEW_LINE('A'),&
!e_edge2e(1,1,:,iface,ielem),&
!e_edge2e(1,2,:,iface,ielem),&
!e_edge2e(1,3,:,iface,ielem),&
!e_edge2e(1,4,:,iface,ielem),&
!NEW_LINE('A'),e_edge2e(2,1,:,iface,ielem),&
!e_edge2e(2,2,:,iface,ielem),&
!e_edge2e(2,3,:,iface,ielem),&
!e_edge2e(2,4,:,iface,ielem)
!write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~'      
        ! Construct the vector of +1 and -1 for the LDG flip-flop
        call create_ldg_flip_flop_sign()

        write(*,*) 'Master node calls metis to subdivide the domain'
        write(*,*) '==============================================================='

        ! Metis calculates the partitions
        call calculatepartitions() 


        write(*,*) 'Master node distributes elements'
        write(*,*) '==============================================================='
!-- DAVID DEBUG START
!ielem = 523
!write(*,*)'original',&
!NEW_LINE('A')//'ox1 = ',vx_Master(1,ic2nh(1,ielem)),&
!NEW_LINE('A')//'oy1 = ',vx_Master(2,ic2nh(1,ielem)),&
!NEW_LINE('A')//'oz1 = ',vx_Master(3,ic2nh(1,ielem)),&
!NEW_LINE('A')//'ox2 = ',vx_Master(1,ic2nh(2,ielem)),&
!NEW_LINE('A')//'oy2 = ',vx_Master(2,ic2nh(2,ielem)),&
!NEW_LINE('A')//'oz2 = ',vx_Master(3,ic2nh(2,ielem)),&
!NEW_LINE('A')//'ox3 = ',vx_Master(1,ic2nh(3,ielem)),&
!NEW_LINE('A')//'oy3 = ',vx_Master(2,ic2nh(3,ielem)),&
!NEW_LINE('A')//'oz3 = ',vx_Master(3,ic2nh(3,ielem)),&
!NEW_LINE('A')//'ox4 = ',vx_Master(1,ic2nh(4,ielem)),&
!NEW_LINE('A')//'oy4 = ',vx_Master(2,ic2nh(4,ielem)),&
!NEW_LINE('A')//'oz4 = ',vx_Master(3,ic2nh(4,ielem)),&
!NEW_LINE('A')//'ox5 = ',vx_Master(1,ic2nh(5,ielem)),&
!NEW_LINE('A')//'oy5 = ',vx_Master(2,ic2nh(5,ielem)),&
!NEW_LINE('A')//'oz5 = ',vx_Master(3,ic2nh(5,ielem)),&
!NEW_LINE('A')//'ox6 = ',vx_Master(1,ic2nh(6,ielem)),&
!NEW_LINE('A')//'oy6 = ',vx_Master(2,ic2nh(6,ielem)),&
!NEW_LINE('A')//'oz6 = ',vx_Master(3,ic2nh(6,ielem)),&
!NEW_LINE('A')//'ox7 = ',vx_Master(1,ic2nh(7,ielem)),&
!NEW_LINE('A')//'oy7 = ',vx_Master(2,ic2nh(7,ielem)),&
!NEW_LINE('A')//'oz7 = ',vx_Master(3,ic2nh(7,ielem)),&
!NEW_LINE('A')//'ox8 = ',vx_Master(1,ic2nh(8,ielem)),&
!NEW_LINE('A')//'oy8 = ',vx_Master(2,ic2nh(8,ielem)),&
!NEW_LINE('A')//'oz8 = ',vx_Master(3,ic2nh(8,ielem))
!-- DAVID DEBUG END

      end if

      !-- pass number_of_possible_partners to all processes
      call mpi_bcast(number_of_possible_partners,1,mpi_integer,0,PETSC_COMM_WORLD,i_err)
      
      ! Push edge connectivity to all processes (has to be before distribute_elements_aflr3 
      !  because in that routine nelems is changed to the local number of elements)
      call distribute_e_edge2e()

      ! Push element connectivity to all processes
      call distribute_elements_aflr3()
!-- DEBUG DAVID START
!if((myprocid.EQ.0))then
!ielem = 253
!write(*,*)'after distribute',&
!NEW_LINE('A')//'adx1 = ',vx(1,e2v(1,ielem)),&
!NEW_LINE('A')//'ady1 = ',vx(2,e2v(1,ielem)),&
!NEW_LINE('A')//'adz1 = ',vx(3,e2v(1,ielem)),&
!NEW_LINE('A')//'adx2 = ',vx(1,e2v(2,ielem)),&
!NEW_LINE('A')//'ady2 = ',vx(2,e2v(2,ielem)),&
!NEW_LINE('A')//'adz2 = ',vx(3,e2v(2,ielem)),&
!NEW_LINE('A')//'adx3 = ',vx(1,e2v(3,ielem)),&
!NEW_LINE('A')//'ady3 = ',vx(2,e2v(3,ielem)),&
!NEW_LINE('A')//'adz3 = ',vx(3,e2v(3,ielem)),&
!NEW_LINE('A')//'adx4 = ',vx(1,e2v(4,ielem)),&
!NEW_LINE('A')//'ady4 = ',vx(2,e2v(4,ielem)),&
!NEW_LINE('A')//'adz4 = ',vx(3,e2v(4,ielem)),&
!NEW_LINE('A')//'adx5 = ',vx(1,e2v(5,ielem)),&
!NEW_LINE('A')//'ady5 = ',vx(2,e2v(5,ielem)),&
!NEW_LINE('A')//'adz5 = ',vx(3,e2v(5,ielem)),&
!NEW_LINE('A')//'adx6 = ',vx(1,e2v(6,ielem)),&
!NEW_LINE('A')//'ady6 = ',vx(2,e2v(6,ielem)),&
!NEW_LINE('A')//'adz6 = ',vx(3,e2v(6,ielem)),&
!NEW_LINE('A')//'adx7 = ',vx(1,e2v(7,ielem)),&
!NEW_LINE('A')//'ady7 = ',vx(2,e2v(7,ielem)),&
!NEW_LINE('A')//'adz7 = ',vx(3,e2v(7,ielem)),&
!NEW_LINE('A')//'adx8 = ',vx(1,e2v(8,ielem)),&
!NEW_LINE('A')//'ady8 = ',vx(2,e2v(8,ielem)),&
!NEW_LINE('A')//'adz8 = ',vx(3,e2v(8,ielem))
!endif
!-- DEBUG DAVID END

    end if

    call mpi_bcast(npoly_max,1,mpi_integer,0,PETSC_COMM_WORLD,i_err)
!-- DAVID DEBUG START
!write(*,*)'+++++++++++++++++++++++++++'
!if((myprocid.EQ.2))then
!iface = 6
!ielem = 6
!  write(*,*)'new'//NEW_LINE('A'),&
!e_edge2e(1,1,:,iface,ielem),&
!e_edge2e(1,2,:,iface,ielem),&
!e_edge2e(1,3,:,iface,ielem),&
!e_edge2e(1,4,:,iface,ielem),&
!NEW_LINE('A'),e_edge2e(2,1,:,iface,ielem),&
!e_edge2e(2,2,:,iface,ielem),&
!e_edge2e(2,3,:,iface,ielem),&
!e_edge2e(2,4,:,iface,ielem)
!endif
!write(*,*)'+++++++++++++++++++++++++++'
!-- DAVID DEBUG END

    ! Initialize collocation approximation
    i_err = rmapInit(npoly,ndim)

    if (i_err > 0) then
      write(*,*) 'An error occured in the function rmapInit!'
      write(*,*) 'Exiting'
      stop
    end if

    ! Assign element orders
    call set_element_orders()      
    
    if (myprocid == 0) then
      write(*,*) 'Each process constructs the coordinates of the collocated nodes'
      write(*,*) '==============================================================='
    end if

    ! Compute coordinates of the collocation points in each elements
    ! ==============================================================
    ! call perturbvertices(0.1_wp)
    ! Collocation points

    call calcnodes_LGL()
!-- uncomment to write solution to file in a way that can be read by Matlab
!   call parallel_write_grid_to_file('true.tf')
!-- end of write to file

    call calc_Gau_shell_pts_all_hexas()

    if (myprocid == 0) then
      write(*,*) 'Each process constructs the metrics'
    end if

    ! Calculate metrics
    call calcmetrics_LGL()

    if (myprocid == 0) then
      write(*,*) 'Each process finds the partner node of each collocated node'
      write(*,*) '==============================================================='
    end if

    call calc_Jacobian_Gau_shell_all_hexas()


    ! Setup collocated nodes connectivity
    call facenodesetup_LGL_Driver()
    call facenodesetup_Gau_Driver()

    if (myprocid == 0) then
      write(*,*) 'Each process finds the WENO partner node of each collocated node'
      write(*,*) '==============================================================='
    end if

    if(discretization == 'SSWENO') then
      ! Setup collocated nodes connectivity
      call facenodesetup_LGL_WENO()
    endif

    if (myprocid == 0) then
      write(*,*) 'Each process initializes PETSc ghost vectors'
      write(*,*) '==============================================================='
    end if

    ! Communicate grid values
    call PetscGridLocations_LGL()

    if(non_conforming .eqv. .true.) call PetscGridLocations_Gau()

    if (myprocid == 0) then
      write(*,*) 'Each process constructs the face-node connectivity'
      write(*,*) '==============================================================='
    end if

    ! Calculate connections
    call calculate_face_node_connectivity_LGL()

    call calculate_face_node_connectivity_Gau()

    if (myprocid == 0) then
      write(*,*) 'Each process constructs the normal vectors'
      write(*,*) '==============================================================='
    end if

    ! Calculate normals
    call calcfacenormals_LGL(.false.)

    call PetscNormals_LGL()

    call calcfacenormals_Gau()

    if((non_conforming .eqv. .true.).AND.(SAT_type.EQ."mod_metric")) call modify_metrics_nonconforming()
   
    call calcfacenormals_LGL(.true.)

    if (myprocid == 0) then
      write(*,*) 'Start actual computation'
      write(*,*) '==============================================================='
    end if

    ! Set binary sizes for writing the solution vtu files in raw binary vtu format
    call calculate_binary_sizes_vtu()

  end subroutine physics_independent_setup

  !============================================================================

  subroutine check_n_procs_n_elems(n_procs,n_elems)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_procs
    integer, intent(in) :: n_elems

    continue

    ! Check if the number of elements is at least equal to the number of
    ! processors passed to mpirun
    if (n_procs .gt. n_elems) then
      write(*,*)
      write(*,*) 'Error: Number of processors is larger than the number of &
      &  elements in the mesh!'
      write(*,*) 'Exiting...'
      write(*,*)
      stop
    endif

    return
  end subroutine

  !============================================================================
  
  !============================================================================

end module physicsindependent
