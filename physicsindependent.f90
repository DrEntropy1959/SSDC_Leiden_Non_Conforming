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
    use referencevariables
    use initgrid
    use fileio
    use mpimod
    use write_solution_file
    
    ! Nothing is implicitly defined
    implicit none

    integer :: i_err

    continue

    ! Initialize collocation approximation
    i_err = rmapInit(npoly,ndim)

    if (i_err > 0) then
      write(*,*) 'An error occured in the function rmapInit!'
      write(*,*) 'Exiting'
      stop
    end if

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

        write(*,*) 'Master node builds connectivity arrays'
!       write(*,*) '==============================================================='

        ! Create connectivity from the original grid
        call e2e_connectivity_aflr3()  

!       write(*,*) 'Master node builds element orders'
!       write(*,*) '==============================================================='
        call set_element_orders_serial()      

        ! Construct the vector of +1 and -1 for the LDG flip-flop
        call create_ldg_flip_flop_sign()

        write(*,*) 'Master node calls metis to subdivide the domain'
        write(*,*) '==============================================================='

        ! Metis calculates the partitions
        call calculatepartitions() 

        write(*,*) 'Master node distributes elements'
        write(*,*) '==============================================================='

      end if

      ! Push element connectivity to all processes
      call distribute_elements_aflr3()

      if (myprocid == 0) then
!       write(*,*) 'read element polynomial orders'
!       write(*,*) '==============================================================='
      endif

      ! Assign element orders
      call set_element_orders()      

    end if
    
    if (myprocid == 0) then
      write(*,*) 'Each process constructs the coordinates of the collocated nodes'
      write(*,*) '==============================================================='
    end if

    ! Compute coordinates of the collocation points in each elements
    ! ==============================================================
    ! call perturbvertices(0.1_wp)
    ! Collocation points
    call calcnodes_LGL()

    if (myprocid == 0) then
      write(*,*) 'Each process constructs the metrics'
    end if

    ! Calculate metrics
    call calcmetrics_LGL()

    if (myprocid == 0) then
      write(*,*) 'Each process finds the partner node of each collocated node'
      write(*,*) '==============================================================='
    end if

!   call calc_Gau_shell_pts_all_hexas()

    ! Setup collocated nodes connectivity
    call facenodesetup_LGL()
    call facenodesetup_Gau()

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

!   call Petsc_shell_Counter()

!   call PetscGridLocations_Gau()

    if (myprocid == 0) then
      write(*,*) 'Each process constructs the face-node connectivity'
      write(*,*) '==============================================================='
    end if

    ! Calculate connections
    call calculate_face_node_connectivity_LGL()
!   call calculate_face_node_connectivity_Gau()

    if (myprocid == 0) then
      write(*,*) 'Each process constructs the normal vectors'
      write(*,*) '==============================================================='
    end if

    ! Calculate normals
    call calcfacenormals_LGL()
!   call calcfacenormals_Gau()

    if (myprocid == 0) then
      write(*,*) 'Start actual computation'
      write(*,*) '==============================================================='
    end if

    ! Set binary sizes for writing the solution vtu files in raw binary vtu format
    call calculate_binary_sizes_vtu()

    return
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
