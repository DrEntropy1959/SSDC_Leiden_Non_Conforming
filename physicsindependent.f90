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
    use non_conforming,      only: h_refine, construct_h_refine_list
    use variables, only : parent_geo, nelems_to_refine
!-- DAVID DEBUG START
    use variables, only : ef2e, vx_master, ic2nh
!-- DAVID DEBUG END

    ! Nothing is implicitly defined
    implicit none

    integer :: i_err
!-- DAVID DEBUG START
    integer :: ielem, iface
!-- DAVID DEBUG END
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
        if(hrefine)then
          !-- if we are doing h-refinement we check after the refinement
        else
          call check_n_procs_n_elems(nprocs,nelems) 
        endif

        write(*,*) 'Master node builds connectivity and orientation arrays'
!      write(*,*) '==============================================================='

        ! Create connectivity from the original grid
        call e2e_connectivity_aflr3()  

        ! Establish face orientation between connected faces
        call face_orientation_aflr3()

!       write(*,*) 'Master node builds element orders if non-conforming'
!       write(*,*) '==============================================================='
        call set_element_orders_serial()    
 
        ! create e_edge2e connectivity
        call e_edge2e_connectivity()   


        ! Construct the vector of +1 and -1 for the LDG flip-flop
        call create_ldg_flip_flop_sign()

        ! h refine
        if(hrefine)then
          call construct_h_refine_list()
          call h_refine()

          write(*,*) 'h-refinement applied, number of hexahedron ',nelems
          write(*,*) '==============================================================='

          ! Check the number of processors and the number of elements
          call check_n_procs_n_elems(nprocs,nelems) 
        endif
        write(*,*) 'Master node calls metis to subdivide the domain'
        write(*,*) '==============================================================='

        ! Metis calculates the partitions
        call calculatepartitions() 

        write(*,*) 'Master node distributes elements'
        write(*,*) '==============================================================='
!-- DEBUG
!do ielem = 1,nelems
!write(*,*)"ielem = ",ielem,"ef2e(4,:,ielem) = ",ef2e(4,:,ielem)
!enddo
!-- DEBUG
      end if

      !-- pass number_of_possible_partners to all processes
      call mpi_bcast(number_of_possible_partners,1,mpi_integer,0,PETSC_COMM_WORLD,i_err)

      !-- pass maximum number of faces per element
      call mpi_bcast(nfacesperelem,1,mpi_integer,0,PETSC_COMM_WORLD,i_err)     

      if(hrefine)then
        !-- pass the number of elements being refined
        call mpi_bcast(nelems_to_refine,1,mpi_integer,0,PETSC_COMM_WORLD,i_err) 
       
        !-- brodcast the parent geometry information
        if(myprocid.NE.0)then
          allocate(parent_geo(3,8,nelems_to_refine))
        endif
        call mpi_bcast(parent_geo(:,:,:),3*8*nelems_to_refine,mpi_double,0,PETSC_COMM_WORLD,i_err)
      else
        ! Push edge connectivity to all processes (has to be before distribute_elements_aflr3 
        !  because in that routine nelems is changed to the local number of elements)
        call distribute_e_edge2e()
      endif

      ! Push element connectivity to all processes
      call distribute_elements_aflr3()
!-- DAVID DEBUG START
!do ielem = 1,nelems
!write(*,*)"============================================="
!  do iface = 1,6
!  write(*,*)"ielem = ",ielem, "iface = ",iface,"ef2e(1,iface,ielem) =&
!",ef2e(1,iface,ielem),"ef2e(2,iface,ielem) = ",ef2e(2,iface,ielem)
!enddo
!write(*,*)"============================================="
!enddo
!        call mpi_barrier(petsc_comm_world,i_err)
!        call PetscFinalize(i_err); stop
!-- DAVID DEBUG END
    end if

    call mpi_bcast(npoly_max,1,mpi_integer,0,PETSC_COMM_WORLD,i_err)


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

!   if(non_conforming .eqv. .true.) call PetscGridLocations_Gau()
    if(non_conforming .eqv. .true.) call Petsc_Gau_Mortar_Geometry_Data()

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
