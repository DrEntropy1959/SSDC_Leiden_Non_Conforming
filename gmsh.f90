  !============================================================================
  ! aflr3ReadUnstructuredGrid - Reads the AFLR3 (.ugrid) grid.

  subroutine aflr3ReadUnstructuredGrid(filein)

    use referencevariables
    use variables, only: vx_master, e2v, nqface, if2nq, ifacetag, ic2nh, &
                       & aflr3_bc_ids, periodic_face_data, wall_face_data
    use controlvariables, only : grid_format

    implicit none

    character(120), intent(in) :: filein

    ! number of bases and zones
    integer :: nbases, nzones

    integer :: iunit
    ! error flag
    integer :: ierr

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

!   integer :: nqface
!   integer , allocatable, dimension(:,:) :: if2nq
!   integer , allocatable, dimension(:)   :: ifacetag
!   integer , allocatable, dimension(:,:) :: ic2nh

    integer :: nnodesg, ntface, ntet, npyr, nprz, nhex

    integer , allocatable, dimension(:,:) :: if2nt

    integer , allocatable, dimension(:,:) :: ic2nt
    integer , allocatable, dimension(:,:) :: ic2np
    integer , allocatable, dimension(:,:) :: ic2nz

    integer :: i, j
    integer :: nbocos
    integer :: group_id
    integer :: i_face, i_group
    integer :: n_p_faces
    integer :: n_w_faces

    character(120) :: header
    integer :: endianess


    real :: version_gmsh
    integer :: file_type_gmsh
    integer :: data_size_gmsh
    integer, allocatable, dimension(:) :: nodes_id_gmsh
    integer :: n_entities_gmsh
    integer, allocatable, dimension(:) :: elems_info_gmsh

    integer :: entity_id_gmsh, entity_type_gmsh, n_entity_tags_gmsh
    integer :: cnt_elems_gmsh, cnt_faces_gmsh


    ! If the meh has been generated with gmsh, then read it and create the same
    ! arrays used with the pure aflr3 format
    if (grid_format == 'gmsh') then
      ! Ensure that the casefile exists
      call get_unit(iunit)
      open(unit=iunit,file=trim(filein)//'.msh',status='old',iostat=ierr)
      if(ierr /= 0)then
        write(*,*) ".msh grid", trim(filein)," doesn't exist."
        stop
      end if
      close(iunit)

      ! If the casefile exists, open the mesh file and read it
      open(unit=iunit,file=trim(filein)//'.msh',status='old',form='formatted', &
        & iostat=ierr)

      ! Read starting section header: $MeshFormat
      read(iunit,*) header
      write(*,*) header
      
      ! Read version number, file_type and data_size
      read(iunit,*) version_gmsh, file_type_gmsh, data_size_gmsh
      write(*,*) version_gmsh, file_type_gmsh, data_size_gmsh
      
      ! Read ending section header: $EndMeshFormat
      read(iunit,*) header
      write(*,*) header

      ! Read starting section header: $Nodes
      read(iunit,*) header
      write(*,*) header

      ! Read number of nodes
      read(iunit,*) nnodesg
      write(*,*) nnodesg

      ! Allocate memory for node_ids. Gmsh assume that the node-number must be 
      ! a postive (non-zero) integer BUT it does not assume that node-numbers 
      ! do not necessarily have to form a dense nor an ordered sequence. 
      allocate(nodes_id_gmsh(nnodesg))
      nodes_id_gmsh = 0

      ! Allocate memory for vertex coordinates
      allocate(vx_master(1:3,1:nnodesg))
      vx_master = 0.0_wp

      ! Read node ID and coordinates
      read(iunit,*) (nodes_id_gmsh(i),vx_master(1,i),vx_master(2,i),vx_master(3,i), &
        & i = 1,nnodesg)
      
      ! Read ending section header: $EndNodes
      read(iunit,*) header
      write(*,*) header

      ! Read starting section header: $Elements
      read(iunit,*) header
      write(*,*) header

      ! Read number of elements (points, lines, quads, hexahedron, etc.)
      read(iunit,*) n_entities_gmsh
      write(*,*) n_entities_gmsh

      ! Allocate memory for element and faces
      ! This is an over estimation because the number of entities is the sum of
      ! number of elements + number of faces + number of lines + number of
      ! points
      ! We assume that only two tags per element are read in the .msh file 
      ! We assume that we have only hexahdreon elements with 8 nodes
      ! See gmhs tutorial for more information
      allocate(elems_info_gmsh(n_entities_gmsh,2,8))

      cnt_elems_gmsh = 0
      cnt_faces_gmsh = 0
      do i = 1, n_entities_gmsh
        
        read(iunit,*) entity_id_gmsh
        read(iunit,*) entity_type_gmsh
        read(iunit,*) n_entity_tags_gmsh
       
        ! Entity type = 5 means hexahedron with 8 nodes
        if (entity_type == 5) then
          cnt_elems_gmsh = cnt_elems_gmsh + 1
          elems_info_gmsh(cnt_elems_gmsh) = entity_id_gmsh
          do j = 1, n_entity_tags_gmsh
            read(iunit,*) elems_info_gmsh(cnt_elems_gmsh,1)
            read(iunit,*) elems_info_gmsh(cnt_elems_gmsh,2)
          end do
        ! Entity type = 3 means quadrilateral face
        else if (entity_type == 3) then
          cnt_faces = cnt_faces + 1
          faces_id_gmsh(cnt) = entity_id_gmsh
          do j = 1, n_entity_tags_gmsh
            read(iunit,*) faces_tags_gmsh(1)
            read(iunit,*) faces_tags_gmsh(2)
          end do
        else
          do j = 1, n_elem_tags_gmsh
            read(iunit,*) trash_int
            read(iunit,*) trash_int
          end do
        end if

        if (elem_type_gmsh == 3)
        read(iunit,*) 

      end do

      stop



    else

      ! ensure that casefile exists
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

      !write(*,*)iunit

      read(iunit) nnodesg, ntface, nqface, ntet, npyr, nprz, nhex
      !write(*,*) nnodesg, ntface, nqface, ntet, npyr, nprz, nhex


      ! allocate vertices
      allocate(vx_master(1:3,1:nnodesg))

      ! allocate faces
      if(ntface /= 0) allocate(if2nt(3,ntface))
      if(nqface /= 0) allocate(if2nq(4,nqface))

      if(ntface+nqface /= 0) allocate(ifacetag(ntface+nqface)) 

      ! allocate elements
      if(ntet /= 0) allocate(ic2nt(4,ntet))
      if(npyr /= 0) allocate(ic2np(5,npyr))
      if(nprz /= 0) allocate(ic2nz(6,nprz))
      if(nhex /= 0) allocate(ic2nh(8,nhex))

      read(iunit) (vx_master(1,i),vx_master(2,i),vx_master(3,i),i=1,nnodesg),                    &
        ((if2nt(j,i),j=1,3),i=1,ntface),                          &
        ((if2nq(j,i),j=1,4),i=1,nqface),                          &
        (ifacetag(i),i=1,ntface+nqface),                          &
        ((ic2nt(j,i),j=1,4),i=1,ntet),                            &
        ((ic2np(j,i),j=1,5),i=1,npyr),                            &
        ((ic2nz(j,i),j=1,6),i=1,nprz),                            &
        ((ic2nh(j,i),j=1,8),i=1,nhex)

      !do i = 1,nnodesg
      !  write(120,'(i6,1x,3(e15.8,1x))')i,vx(:,i)
      !enddo
      !do i = 1,nqface
      !  write(120,'(i6,1x,5(i5,1x))')i,if2nq(:,i),ifacetag(i)
      !enddo

      nvertices        = nnodesg
      nverticesperelem = 2**ndim
      nelems           = ntet + npyr + nprz + nhex

      ! Read BC information from the aflr3 (.b8.mapbc or .mapbc) grid file
      call aflr3_read_bc(filein)

    end if

    ! Re-set ifacetag to be the BC number. This allows to re-use the rest of the
    ! code without changing anything. In fact, ifacetag will be used then to set
    ! ef2e(2,i,j).
    n_p_faces = 0
    n_w_faces = 0

    do i_face = 1, size(ifacetag)

      group_id = ifacetag(i_face)

      do i_group = 1, size(aflr3_bc_ids(:,1))
        if (group_id == aflr3_bc_ids(i_group,1)) then
          ifacetag(i_face) = aflr3_bc_ids(i_group,2)
        end if
      end do

      ! Count number of "periodic" faces
      if (ifacetag(i_face) == 8 .or. ifacetag(i_face) == 9) then
        n_p_faces = n_p_faces + 1
      end if

      ! Count number of "wall" faces
      if (ifacetag(i_face) == 5 .or. ifacetag(i_face) == 6) then
        n_w_faces = n_w_faces + 1
      end if
    
    end do

    ! Check if the number of "periodic" faces is even
    if (MOD(n_p_faces,2) /= 0) then
      write(*,*) 'Number of "periodic" faces must be even.'
      write(*,*) 'Exiting...'
      stop
    end if

    write(*,*) 'Number of PERIODIC FACES', n_p_faces

    ! Allocate memory for periodic element face data
    allocate(periodic_face_data(4+nverticesperface,n_p_faces))
    periodic_face_data = 0

    ! Allocate memory for wall element face data
    allocate(wall_face_data(2,n_w_faces))
    wall_face_data = 0

    return
  end subroutine aflr3ReadUnstructuredGrid

  !============================================================================
 
