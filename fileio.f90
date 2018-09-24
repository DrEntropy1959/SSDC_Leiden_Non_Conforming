module fileio
  use precision_vars
  implicit none

  private

  public readSSDCstartup
  public init_code_bc_types
  public aflr3ReadUnstructuredGrid
  
contains

  !============================================================================
  
  subroutine readSSDCstartup(file_name)

    ! Load module
    use referencevariables
    use controlvariables
    use nsereferencevariables
    use collocationvariables
    use SSWENOvariables

    ! Nothing is implicitly defined
    implicit none

    character(len=120), intent(in) :: file_name

    integer :: iunit, ierr


    namelist /PreProcParameters/ ndim, nequations, nprocs, npoly, npoly_max, npoly_DeltaF, &
      &casefile, grid_format, viscous, crossterms, RK_Method, IMEX_element, &
      &IMEX_penalty, physics, discretization, filter_solution, l01, l10, l00, &
      &periodic, periodic_distance, Entropy_Correction, variable_viscosity, &
      &heat_entropy_flow_wall_bc, flux_entropy_correction, alpha_ldg_flip_flop, &
      &Riemann_Diss, Riemann_Diss_BC, WENO_Bias, WENO_type, entropy_viscosity,  &
      &Grid_Topology, cylinder_x0, cylinder_x1, p_non_conforming, turbulent_viscosity, &
      &p_refine_strategy,entropy_flux,entropy_flux_BC,radius,origin,SAT_type,symmetric_metric,&
      hrefine,h_refine_strategy,check_negative_P_T

    namelist /GeneralParameters/ runcasename, solutionfile, outputfile, new, &
      & time_averaging, &
      & write_restart, write_restart_dir, write_restart_common_name, &
      & write_restart_formatted, write_restart_frequency, read_restart_dir, &
      & read_restart_common_name, read_restart_time, read_restart_formatted, &
      & read_restart_time_averaging, &
      & write_solution, write_solution_dir, write_solution_common_name, &
      & write_solution_formatted, write_solution_frequency, &
      & write_aero_coeffs, write_dke_dt, write_errors_wall, &
      & perturb_internal_vertices

    namelist /FlowParameters/ Re0, Re0inv, Pr0, Mach0, uniformFreeStreamAOA, &
      & membranelocation, referenceWaveSpeed, isothermalWallTemperature, &
      & InitialCondition

    namelist /AerodynamicParameters/ aero_coeffs_surface

    continue

    if (myprocid .eq. 0) then
      write(*,*) 'Reading the start-up file'
      write(*,*) '=========================================================================='
    end if


    ! Open the SSDCstartup file
    call get_unit(iunit) ! (precision_vars)
    open(unit = iunit, file = file_name, status = 'old', action = 'read', iostat = ierr)
    if(ierr /= 0)then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) "Fatal Error: the following start-up file does not exist: ", file_name
        write(*,*) 'Exiting...'
!      endif
      stop
    end if

    ! Read the namelist for PreProcParameters
    read(iunit,PreProcParameters,iostat=ierr)
    if((Grid_Topology == 'cylinder') .and. (npoly == 1)) then
      Grid_Topology = 'linear'
      if(myprocid == 0) write(*,*) 'Inconsistancy between Grid_Topology and npoly. Switching Topology to linear'
    endif
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read PreProcParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if

    if (myprocid .eq. 0) then
      write(*,*) 'Polynomial order: ', npoly
      write(*,*) '==============================================================='
      write(*,*) 'Over-collocation: ', npoly_DeltaF
      write(*,*) '==============================================================='
    end if

    ! Read the namelist for GeneralParameters
    read(iunit,GeneralParameters,iostat=ierr)
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read GeneralParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if 

    ! Read the namelist for GeneralParameters
    read(iunit,FlowParameters,iostat=ierr)
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read FlowParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if

    ! Read the namelist for AerodynamicParameters
    read(iunit,AerodynamicParameters,iostat=ierr)
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read AerodynamicParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if 

    close(iunit)

    
    ! Check consistency of the input parameters
    ! -----------------------------------------
    ! -----------------------------------------

    ! Check Reynolds number for viscous flow
    if (viscous .and. Re0 < 1e-14) then
      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file: You have not specified the', &
          & 'reference Reynolds number or you have set it to zero!' 
        write(*,*) 'Exiting...'
      endif
      stop
    endif

    ! Delay the code so that we can think 2 more seconds if we really want to
    ! overwrite the restart files
    if (write_restart) then
      if (myprocid==0) then
        write(*,*)
        write(*,*) 'WARNING: In 2 seconds the execution of the code will', & 
          & ' continue. If you did not change the restart_dir name or the', & 
          & ' common_restart_name, all the files saved in the default', & 
          & ' directory or in the directory specified in the start-up file', &
          & ' will be overwritten.', &
          & ' Press Ctrl+c to stop the code execution. GOOD LUCK!'
        write(*,*)
      endif

      call sleep(2)
    
    endif
    
    ! Print warning when the user decides to restart a simulation
    if (.not. new) then
      if (myprocid==0) then
        write(*,*)
        write(*,*) 'WARNING: You have chosen to restart the simulation.', &
          & ' Please, verify the details of the restarting files, i.e.,', &
          & ' read_restart_dir, read_restart_common_name, read_restart_time', &
          & ' and read_restart_formatted in the SSDCstartup file.', &
          & ' Press Ctrl+c to stop the code execution.'
        write(*,*)
      endif

      call sleep(2)

    endif

    ! Call run time modifiable startup parameters subroutine
    call read_run_time_modifiable_startup_parameters(file_name)

    return
  end subroutine readSSDCstartup

  !============================================================================
  
  !============================================================================
  
  subroutine read_run_time_modifiable_startup_parameters(file_name)
    
    ! Load modules
    use controlvariables

    ! Nothing is implicitly defined
    implicit none

    character(len=120), intent(in) :: file_name

    integer :: iunit, ierr

    namelist /RunParameters/ timestep, timemaximum , verbose, verbose_fac, &
                           & variabletimestep, CFL, Dt_by_CFL

    continue

    ! Open the SSDCstartup file
    call get_unit(iunit) ! (precision_vars)
    open(unit=iunit,file=file_name,status='old',action='read',iostat=ierr)
    if(ierr /= 0)then
      write(*,*)
      write(*,*) "Fatal Error (modifiable startup parameters): the following start-up file does not exist: ", file_name
      write(*,*) 'Exiting...'
      stop
    end if


    ! read the namelist for RunParameters
    read(iunit,RunParameters,iostat=ierr)
    if(ierr /= 0) then
      write(*,*) 'Error in startup file, cannot read RunParameters!', ierr
      write(*,*) 'Exiting...'
      stop
    end if

    close(iunit)

    return
  end subroutine read_run_time_modifiable_startup_parameters

  !============================================================================
  
  !============================================================================
  ! aflr3ReadUnstructuredGrid - Reads the AFLR3 (.ugrid) grid.

  subroutine aflr3ReadUnstructuredGrid(filein)

    use referencevariables
    use variables, only: vx_master, nqface, if2nq, ifacetag, ic2nh, &
                       & aflr3_bc_ids, periodic_face_data_x1, &
                       & periodic_face_data_x2, periodic_face_data_x3, &
                       & wall_face_data 


    !  nnodesg    = number of nodes
    !  ntface     = number of boundary faces with triangles
    !  nqface     = number of boundary faces with quads
    !  ntet       = number of tetrahedra
    !  npyr       = number of pyramids
    !  nprz       = number of prisms
    !  nhex       = number of hexahedra
    !  x,y,z      = node coordinates
    !  if2nt      = face-to-node pointers for each TRI boundary face (3 nodes
    !               are pointed to by each triangle)
    !  if2nq      = face-to-node pointers for each QUAD boundary face (4 nodes
    !               are pointed to by each quadralateral)
    !  ifacetag   = tag number that groups boundary faces (e.g., all faces
    !               with the number "1" correspond to the same boundary
    !               entity (like wing) with a given boundary condition)
    !  ic2nt      = cell-to-node pointers for tets (4 nodes are pointed to by
    !               each tetrahedral cell)
    !  ic2np      = cell-to-node pointers for pyramids (5 nodes are pointed to by
    !               each pyramidal cell)
    !  ic2nz      = cell-to-node pointers for prisms (6 nodes are pointed to by
    !               each prismatic cell)
    !  ic2nh      = cell-to-node pointers for hexes (8 nodes are pointed to by
    !               each hexahedral cell)
    !  For my use:
    !  ibcnum     = number of surface elements associated with each
    !               number (e.g., ibcnum(1) = number of elements labeled
    !               as BC type=1, ibcnum(n) = number of elements labeled
    !               as BC type=n)
    !  itag       = element numbers of all the surface triangles
    !               (e.g., itag(1,n) = nth triangle element number in group
    !               of BC type = 1, itag(m,n) = nth triangle element number in
    !               group of BC type = m)
    !

    implicit none                                         ! Nothing is implicitly defined

    character(120), intent(in) :: filein

    integer :: iunit
    integer :: ierr

    integer :: nnodesg, ntface, ntet, npyr, nprz, nhex

    integer , allocatable, dimension(:,:) :: if2nt
    integer , allocatable, dimension(:,:) :: ic2nt
    integer , allocatable, dimension(:,:) :: ic2np
    integer , allocatable, dimension(:,:) :: ic2nz

    integer :: i, j
    integer :: group_id
    integer :: i_face, i_group
    integer :: n_p_faces_x1_a, n_p_faces_x1_b, &
               n_p_faces_x2_a, n_p_faces_x2_b, &
               n_p_faces_x3_a, n_p_faces_x3_b
    integer :: n_p_faces_x1, n_p_faces_x2, n_p_faces_x3
    integer :: n_w_faces
    logical :: match

    continue

    ! Ensure that casefile exists
    call get_unit(iunit)
    open(unit=iunit,file=trim(filein)//'.b8.ugrid', &
!   open(unit=iunit,file=trim(filein)//'.ugrid',    &
      status='old',                                 &
      access='stream',                              &
     !convert='big_endian',                         &
      iostat=ierr)
    if(ierr /= 0)then
      write(*,*) "aflr3GridRead:", trim(filein)," doesn't exist."
      stop
    end if
    close(iunit)

    open(unit=iunit,file=trim(filein)//'.b8.ugrid', &
!   open(unit=iunit,file=trim(filein)//'.ugrid', &
      status='old',                                 &
      access='stream',                              &
      !convert='big_endian',                         &
      iostat=ierr)

    read(iunit) nnodesg, ntface, nqface, ntet, npyr, nprz, nhex

    ! Allocate vertices
    allocate(vx_master(1:3,1:nnodesg))

    ! Allocate faces
    if(ntface /= 0) allocate(if2nt(3,ntface))
    if(nqface /= 0) allocate(if2nq(4,nqface))

    if(ntface+nqface /= 0) allocate(ifacetag(ntface+nqface)) 

    ! Allocate elements
    if(ntet /= 0) allocate(ic2nt(4,ntet))
    if(npyr /= 0) allocate(ic2np(5,npyr))
    if(nprz /= 0) allocate(ic2nz(6,nprz))
    if(nhex /= 0) allocate(ic2nh(8,nhex))

    read(iunit) (vx_master(1,i),vx_master(2,i),vx_master(3,i),i=1,nnodesg), &
      ((if2nt(j,i),j=1,3),i=1,ntface),                                      &
      ((if2nq(j,i),j=1,4),i=1,nqface),                                      &
      (ifacetag(i),i=1,ntface+nqface),                                      &
      ((ic2nt(j,i),j=1,4),i=1,ntet),                                        &
      ((ic2np(j,i),j=1,5),i=1,npyr),                                        &
      ((ic2nz(j,i),j=1,6),i=1,nprz),                                        &
      ((ic2nh(j,i),j=1,8),i=1,nhex)

    ! Set number of nodes, vertices per element, and elements (serial)
    nvertices        = nnodesg
    nverticesperelem = 2**ndim
    nelems           = ntet + npyr + nprz + nhex

    ! Read BC information from the aflr3 (.b8.mapbc or .mapbc) grid file
    call aflr3_read_bc(filein)

    ! Re-set ifacetag to be the BC number. This allows to re-use the rest of the
    ! code without changing anything. In fact, ifacetag will be used then to set
    ! ef2e(2,i,j).
    ! =========================================================================
    do i_face = 1, size(ifacetag)

      ! Original face tag ID from AFLR3
      group_id = ifacetag(i_face)

      ! Initialize logical variable
      match = .false.

      do i_group = 1, size(aflr3_bc_ids(:,1))
        if (group_id == aflr3_bc_ids(i_group,1)) then
          ! Match found
          match = .true.

          ! Re-map ifacetage to the convention used in the code
          ifacetag(i_face) = aflr3_bc_ids(i_group,2)
        end if
      end do

      ! If the tag of the boundary face is not known throw an error
      if (match .eqv. .false.) then
        write(*,*) 'Unknown boundary face tag: ', group_id
        write(*,*) 'Please check the number of the boundary face tag.'
        write(*,*) 'Exting...'
        stop
      end if
    end do

    ! Count periodic and wall boundary faces
    ! =========================================================================
    ! Initialize number of periodic and wall boundary faces
    n_p_faces_x1_a = 0
    n_p_faces_x1_b = 0
    n_p_faces_x2_a = 0
    n_p_faces_x2_b = 0
    n_p_faces_x3_a = 0
    n_p_faces_x3_b = 0
    n_p_faces_x1   = 0
    n_p_faces_x2   = 0
    n_p_faces_x3   = 0
    n_w_faces = 0

    do i_face = 1, size(ifacetag)
      ! Count number of "periodic" faces in the x1 direction
      if (ifacetag(i_face) == 8) then 
        n_p_faces_x1_a = n_p_faces_x1_a + 1
        n_p_faces_x1   = n_p_faces_x1   + 1
      end if
      if (ifacetag(i_face) == 9) then  
        n_p_faces_x1_b = n_p_faces_x1_b + 1
        n_p_faces_x1   = n_p_faces_x1   + 1
      end if

      ! Count number of "periodic" faces in the x2 direction
      if (ifacetag(i_face) == 10) then 
        n_p_faces_x2_a = n_p_faces_x2_a + 1
        n_p_faces_x2   = n_p_faces_x2   + 1
      end if
      if (ifacetag(i_face) == 11) then  
        n_p_faces_x2_b = n_p_faces_x2_b + 1
        n_p_faces_x2   = n_p_faces_x2   + 1
      end if

      ! Count number of "periodic" faces in the x3 direction
      if (ifacetag(i_face) == 12) then 
        n_p_faces_x3_a = n_p_faces_x3_a + 1
        n_p_faces_x3   = n_p_faces_x3   + 1
      end if
      if (ifacetag(i_face) == 13) then  
        n_p_faces_x3_b = n_p_faces_x3_b + 1
        n_p_faces_x3   = n_p_faces_x3   + 1
      end if

      ! Count number of "wall" faces
      if (ifacetag(i_face) == 5 .or. ifacetag(i_face) == 6) then
        n_w_faces = n_w_faces + 1
      end if
    
    end do

    ! Check if the number of "periodic" faces in the x1 direction is correct
    if (n_p_faces_x1_a .ne. n_p_faces_x1_b) then
      write(*,*) 'Number of "periodic" faces in the x1_a plane does not match', &
        ' the number of "periodic" faces in the x1_b plane.'
      write(*,*) n_p_faces_x1_a, n_p_faces_x1_b
      write(*,*) 'Exiting...'
      stop
    end if
    
    ! Check if the number of "periodic" faces in the x2 direction is correct
    if (n_p_faces_x2_a .ne. n_p_faces_x2_b) then
      write(*,*) 'Number of "periodic" faces in the x2_a plane does not match', &
        ' the number of "periodic" faces in the x2_b plane.'
      write(*,*) n_p_faces_x2_a, n_p_faces_x2_b
      write(*,*) 'Exiting...'
      stop
    end if

    ! Check if the number of "periodic" faces in the x3 direction is correct
    if (n_p_faces_x3_a .ne. n_p_faces_x3_b) then
      write(*,*) 'Number of "periodic" faces in the x3_a plane does not match', &
        ' the number of "periodic" faces in the x3_b plane.'
      write(*,*) n_p_faces_x3_a, n_p_faces_x3_b
      write(*,*) 'Exiting...'
      stop
    end if

    ! Allocate memory for periodic element face data
    allocate(periodic_face_data_x1(4+nverticesperface,n_p_faces_x1)) ; periodic_face_data_x1 = 0 ;

    allocate(periodic_face_data_x2(4+nverticesperface,n_p_faces_x2)) ; periodic_face_data_x2 = 0 ;

    allocate(periodic_face_data_x3(4+nverticesperface,n_p_faces_x3)) ; periodic_face_data_x3 = 0 ;

    ! Allocate memory for wall element face data
    allocate(wall_face_data(2,n_w_faces)) ; wall_face_data = 0 ;

  end subroutine aflr3ReadUnstructuredGrid

  !============================================================================
  
  !============================================================================
  ! aflr3_read_bc - Reads the AFLR3 boundary conditions file and call the 
  ! subroutine that maps the AFLR3 boundary conditions numbering to the 
  ! numbering used in the code.

  subroutine aflr3_read_bc(file_in)
    
    ! Load modules
    use variables, only : aflr3_bc_ids 

    ! Nothing is implicitly defined
    implicit none

    character(120), intent(in) :: file_in
    integer :: i_unit
    integer :: i_err
    integer :: n_bc_groups
    integer :: i_group
    integer :: group_id, group_bc_id
    character(120) :: group_bc_name 
    character(120), allocatable, dimension(:) :: aflr3_bc_names

    continue

    ! Ensure that the  mapbc file exists
    ! =========================================================================
    call get_unit(i_unit)

    open(unit=i_unit,file=trim(file_in)//'.mapbc', status='old', &
      & action='READ', iostat=i_err)

    if(i_err /= 0)then
      write(*,*) "aflr3_read_bc:", trim(file_in)," doesn't exist."
      stop
    end if

    close(i_unit)

    ! Read mapbc file
    ! =========================================================================
    open(unit=i_unit,file=trim(file_in)//'.mapbc', status='old', &
      & action='READ', iostat=i_err)

    ! Read number of BC groups
    read(i_unit,*) n_bc_groups

    ! Allocate memory for AFLR3 BC info
    allocate(aflr3_bc_ids(n_bc_groups,2))
    allocate(aflr3_bc_names(n_bc_groups))

    do i_group = 1, n_bc_groups 
      ! Read group informations
      read(i_unit,*) group_id, group_bc_id, group_bc_name

      ! Assign information to the output array
      aflr3_bc_ids(i_group,1) = group_id
      aflr3_bc_ids(i_group,2) = group_bc_id
      aflr3_bc_names(i_group) = group_bc_name
    end do

    close(i_unit)

    ! Re-map the aflr3 group_bc_id to be the BC convention numbering used in
    ! this code
    ! =========================================================================
    call map_aflr3_bc_id(n_bc_groups,aflr3_bc_ids,aflr3_bc_names)

    return
  end subroutine aflr3_read_bc

  !============================================================================
  
  !============================================================================
  ! map_aflr3_bc_id - Maps the AFLR3 BC numbering to the convention used in the
  ! code.

  subroutine map_aflr3_bc_id(n_groups,bc_ids,bc_names)

    ! Load modules
    use variables, only : code_bc_names
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_groups
    integer, dimension(n_groups,2), intent(inout) :: bc_ids
    character(120), dimension(n_groups), intent(in) :: bc_names
    character(120) :: group_bc_name
    logical :: match
    integer :: i_group, i_code_bc

    continue

    do i_group = 1, n_groups
      ! Group BC name from AFLR3 grid
      group_bc_name = bc_names(i_group)

      ! Initialize logical variable to false
      match = .false.

      do i_code_bc = 1, size(code_bc_names)
        ! Check if the group BC name from AFLR3 matches one of the BC supported
        ! in the code
        if (group_bc_name == code_bc_names(i_code_bc)) then
          bc_ids(i_group,2) = i_code_bc
          match = .true.
        end if
      end do
      if (match .eqv. .false.) then
        write(*,*) 'Unknown BC type: ', group_bc_name
        write(*,*) 'Please check the name of the boundary condition or', &
          & ' implement this new type.'
        write(*,*) 'Exting...'
        stop
      end if
    end do

    return
  end subroutine map_aflr3_bc_id

  !============================================================================
  
  !============================================================================
  ! init_code_bc_types - Defines the IDs of the boundary conditions according to
  ! this code.

  subroutine init_code_bc_types()
    
    ! Load modules
    use variables, only : code_bc_names

    ! Nothing is implicitly defined
    implicit none

    integer, parameter :: n_code_bc_types = 16

    continue 

    ! Allocate memory for code BC names and IDs
    allocate(code_bc_names(n_code_bc_types))
    
    ! No BC
    code_bc_names(1) = 'null'

    ! Dirichlet BC. The actual BC subroutine will depends on the initial
    ! solution specific in the start-up file
    code_bc_names(2) = 'dirichlet'

    ! Uniform free-stream BC
    code_bc_names(3) = 'uniform_free_stream'

     ! Inviscid wall BC
    code_bc_names(4) = 'inviscid_wall'

    ! No slip adiabatic wall BC
    code_bc_names(5) = 'no_slip_adiabatic_wall'

    ! No slip wall entropy stable BC
    code_bc_names(6) = 'no_slip_wall_ss'

    ! Symmetry plane BC
    code_bc_names(7) = 'symmetry_plane'

    ! Periodic BC: 1a plane 
    code_bc_names(8) = 'periodic_x1_a'
    
    ! Periodic BC: 1b plane
    code_bc_names(9) = 'periodic_x1_b'

    ! Periodic BC: 2a plane 
    code_bc_names(10) = 'periodic_x2_a'
    
    ! Periodic BC: 2b plane
    code_bc_names(11) = 'periodic_x2_b'

    ! Periodic BC: 3a plane 
    code_bc_names(12) = 'periodic_x3_a'
    
    ! Periodic BC: 3b plane
    code_bc_names(13) = 'periodic_x3_b'

    ! Boundary Layer Inflow Profile
    code_bc_names(14) = 'BLayerProfile'

    ! Subsonic Pressure Boundary condition 
    code_bc_names(15) = 'SubsonicOutflow'

    ! Elevate the polynomial of the element on the boundaryBoundary condition 
    code_bc_names(16) = 'dirichlet_elevate'

    return
  end subroutine init_code_bc_types
          
  !============================================================================
  
end module fileio
