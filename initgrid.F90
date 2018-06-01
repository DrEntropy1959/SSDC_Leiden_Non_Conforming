module initgrid
  use precision_vars
  use iso_c_binding
  implicit none

  ! #include "finclude/petscsys.h"
  ! #include "finclude/petscvec.h"
  ! #include "finclude/petscdmda.h"
  ! #include "finclude/petscis.h"
  ! #include "finclude/petscmat.h"
  ! #include "finclude/petscksp.h"
  ! #include "finclude/petscpc.h"
  ! #include "finclude/petscsnes.h"

  private
  public init_edge_2
  public init_quad_4
  public init_hex_8
  public E2EConnectivity_cgns
  public e2e_connectivity_aflr3
  public face_orientation_aflr3
  public set_element_orders
  public set_element_orders_Serial
  public calculatepartitions

  public calcnodes_LGL
  public calcmetrics_LGL
  public facenodesetup_LGL_Driver
  public facenodesetup_LGL
  public facenodesetup_Gau_Driver
  public facenodesetup_Gau
  public facenodesetup_LGL_WENO
  public calculate_face_node_connectivity_LGL
  public calculate_face_node_connectivity_Gau
  public calcfacenormals_LGL
  public calcfacenormals_Gau
  public init_elem_type
  public create_ldg_flip_flop_sign
  public pert_int_vert
  public face_pairs
  public data_partner_element_serial
  public WENO_Adjoining_Data
  public Pencil_Coord
  public WENO_Intrp_Face_Nodes
  public Boundary_Vertex_2_Vertex_Connectivity
  public calc_Gau_shell_pts_all_hexas
  public calc_Jacobian_Gau_shell_all_hexas
  public modify_metrics_nonconforming
  public perturb_vertices_tg_vortex_1
  public transform_grid
  public write_matrix_to_file_matlab
  public e_edge2e_connectivity
  public map_face_orientation_k_On_2_k_Off
  public face_map

  integer, allocatable, dimension(:,:), target :: edge_2_faces
  integer, allocatable, dimension(:), target :: edge_2_facedirections
  integer, parameter :: edge_2_nfacesperelem = 2
  integer, parameter :: edge_2_nverticesperface = 1
  integer, allocatable, dimension(:,:), target :: quad_4_faces
  integer, allocatable, dimension(:), target :: quad_4_facedirections
  integer, parameter :: quad_4_nfacesperelem = 4
  integer, parameter :: quad_4_nverticesperface = 2
  integer, allocatable, dimension(:,:), target :: hex_8_faces
  integer, allocatable, dimension(:,:), target :: hex_8_faces_order
  integer, allocatable, dimension(:), target :: hex_8_facedirections
  integer, parameter :: hex_8_nfacesperelem = 6
  integer, parameter :: hex_8_nverticesperface = 4
  integer, pointer, dimension(:,:) :: eltypfaces
  integer, pointer, dimension(:,:) :: eltypfaces_Lexo
  integer, pointer, dimension(:) :: elfacedirections

  public elfacedirections

!-- DAVID DEBUG START
  public face_orientation_hex
!-- DAVID DEBUG END
  interface
    integer(c_int) function calcMetisPartitions(ne, nv, nps, xadj, adj, ncommon, &
        epart, npart, nnz) &
        bind(c,name='calcMetisPartitions') 
      ! this is the interface to the c function that calls metis
      use iso_c_binding, only: c_ptr, c_double, c_int, c_char
      implicit none
      integer(c_int), value :: ne
      integer(c_int), value :: nv
      integer(c_int), value :: nps
      type(c_ptr), value :: xadj
      type(c_ptr), value :: adj
      integer(c_int), value :: ncommon
      type(c_ptr), value :: epart
      type(c_ptr), value :: npart
      integer(c_int), value :: nnz
    end function calcMetisPartitions
  end interface

contains

  subroutine init_edge_2()    !   SERIAL Routine

    ! we require the following info about the edge_2 elements
    ! to fully define the connectivity. This is in compliance
    ! with the CGNS standard
    use referencevariables, only: nverticesperface, nfacesperelem
    implicit none

    ! number of vertices per face
    nverticesperface = edge_2_nverticesperface
    ! number of faces per element
    nfacesperelem = edge_2_nfacesperelem
    ! the local vertices that are used to construct each face
    allocate(edge_2_faces(nverticesperface,nfacesperelem))
    edge_2_faces(:,1) = (/ 1 /)
    edge_2_faces(:,2) = (/ 2 /)
    ! the outward signed computational direction of each face
    allocate(edge_2_facedirections(nfacesperelem))
    edge_2_facedirections(1) = -1
    edge_2_facedirections(2) = 1

  end subroutine init_edge_2   !  SERIAL Routine

  subroutine init_quad_4()     !  SERIAL Routine

    ! we require the following info about the quad_4 elements
    ! to fully define the connectivity. This is in compliance
    ! with the CGNS standard
    use referencevariables, only: nverticesperface, nfacesperelem
    implicit none

    ! number of vertices per face
    nverticesperface = quad_4_nverticesperface
    ! number of faces per element
    nfacesperelem = quad_4_nfacesperelem
    ! the local vertices that are used to construct each face
    allocate(quad_4_faces(nverticesperface,nfacesperelem))
    quad_4_faces(:,1) = (/ 1, 2 /)
    quad_4_faces(:,2) = (/ 2, 3 /)
    quad_4_faces(:,3) = (/ 3, 4 /)
    quad_4_faces(:,4) = (/ 4, 1 /)
    ! the outward signed computational direction of each face
    allocate(quad_4_facedirections(nfacesperelem))
    quad_4_facedirections(1) = -2
    quad_4_facedirections(2) = +1
    quad_4_facedirections(3) = +2
    quad_4_facedirections(4) = -1

  end subroutine init_quad_4   !  SERIAL Routine

  subroutine init_hex_8()   !  SERIAL Routine
!
!                        8-----------------7
!                       /.                /|
!                      / .               / |
!                     /  .              /  |
!       (6)--top------------->         /   |
!                   /    .       (4)  /    |
!                  /     .           /     |
!                 /      .          /      |
!                /       .         /       |
!               5-----------------6        |       
!               |   (5)  4........|...(3)..3
!               |       .         |       /
!               |      .          |      /
!               |     .           |     /  zeta ^     
!       (2)-front------->    (1)  |    /        |    / eta
!               |   .             |   /         |   /
!               |  .              |  /          |  /
!               | .               | /           | /
!               |.                |/            |/
!               1-----------------2             ----------> xi

    ! we require the following info about the hex_8 elements
    ! to fully define the connectivity. This is in compliance
    ! with the CGNS standard
    use referencevariables, only: nverticesperface, nfacesperelem
    implicit none

    ! number of vertices per face
    nverticesperface = hex_8_nverticesperface
    ! number of faces per element
    nfacesperelem = hex_8_nfacesperelem
    ! the local vertices that are used to construct each face
    allocate(hex_8_faces(nverticesperface,nfacesperelem))
    hex_8_faces(:,1) = (/ 1, 4, 3, 2 /)
    hex_8_faces(:,2) = (/ 1, 2, 6, 5 /)
    hex_8_faces(:,3) = (/ 2, 3, 7, 6 /)
    hex_8_faces(:,4) = (/ 3, 4, 8, 7 /)
    hex_8_faces(:,5) = (/ 1, 5, 8, 4 /)
    hex_8_faces(:,6) = (/ 5, 6, 7, 8 /)

    allocate(hex_8_faces_order(nverticesperface,nfacesperelem))
    hex_8_faces_order(:,1) = (/ 1, 2, 4, 3 /)
    hex_8_faces_order(:,2) = (/ 1, 2, 5, 6 /)
    hex_8_faces_order(:,3) = (/ 2, 3, 6, 7 /)
    hex_8_faces_order(:,4) = (/ 4, 3, 8, 7 /)
    hex_8_faces_order(:,5) = (/ 1, 4, 5, 8 /)
    hex_8_faces_order(:,6) = (/ 5, 6, 8, 7 /)

    ! the outward signed computational direction of each face
    allocate(hex_8_facedirections(nfacesperelem))
    hex_8_facedirections(1) = -3
    hex_8_facedirections(2) = -2
    hex_8_facedirections(3) = +1
    hex_8_facedirections(4) = +2
    hex_8_facedirections(5) = -1
    hex_8_facedirections(6) = +3

  end subroutine init_hex_8   !  SERIAL Routine

  function isort(iain,nl)
    integer, intent(in) :: nl
    integer, intent(in) :: iain(nl)

    integer :: isort(nl)

    integer :: iatmp(nl)

    integer :: i, imin(1), iadd

    iatmp = iain
    iadd = maxval(iain)+100
    do i = 1, nl
      imin = minloc(iatmp)
      isort(i) = iatmp(imin(1))
      iatmp(imin(1)) = iadd
    end do

  end function isort

  subroutine init_elem_type()   !  SERIAL Routine
    use referencevariables
    use variables, only: facenormalcoordinate

    implicit none
    integer :: i

    ! Limited to tensor product elements at this time
    if (ndim == 2) then
      eltypfaces       => quad_4_faces
      elfacedirections => quad_4_facedirections
      nverticesperface =  quad_4_nverticesperface
      nfacesperelem    =  quad_4_nfacesperelem
      allocate(facenormalcoordinate(nfacesperelem))
      do i = 1,nfacesperelem
        facenormalcoordinate(i) = quad_4_facedirections(i)
      enddo
    else if (ndim == 3) then
      eltypfaces       => hex_8_faces
      eltypfaces_Lexo  => hex_8_faces_order
      elfacedirections => hex_8_facedirections
      nverticesperface =  hex_8_nverticesperface
      nfacesperelem    =  hex_8_nfacesperelem
      allocate(facenormalcoordinate(nfacesperelem))
      do i = 1,nfacesperelem
        facenormalcoordinate(i) = hex_8_facedirections(i)
      enddo
    else
      write(*,*) 'error: unsupported dimension'
      stop
    end if

  end subroutine init_elem_type   !  SERIAL Routine

  subroutine E2EConnectivity_cgns()   !   SERIAL Routine

!   SERIAL ROUTINE
!     CGNS Element to element connectivity
!   SERIAL ROUTINE

    use referencevariables
    use variables, only: boundaryelems, ef2e,  &
                         iae2v, jae2v, nnze2v
!                        iae2e, jae2e, nnze2e, 
    implicit none
    integer :: i,j,k,izone,jzone,kzone
    integer :: ii,jj,kk,icount,ielem
    integer :: iface, jface

    integer, allocatable :: v2e(:,:), iav2e(:)
    integer, allocatable, dimension(:) :: ivtmp1, ivtmp2
    integer :: nnzv2e
    integer :: connectioncount, vertexcount

    !  routine called from myprocid == 0 ;  i.e. master works on entire grid
    !
    !  nelems     :    elements
    !  nelemzones :    zones = 2,  1 == interior; 2 == boundary 
    !
    !                   Dim,    Dim,         Dim
    !  ef2e       :    ( 2 ,nfaceperelem, nelements) 
    !             :  Two situation occur.  The face is either an 
    !                  (Interior face 
    !                      :  (1,j,k) = face ID of the adjoining element
    !                      :  (2,j,k) = Connected to Element 
    !                  (Boundary face 
    !                      :  (1,j,k) = Set to -11 
    !                      :  (2,j,k) = Boundary Face number
    !
    !         (Note:  redimensioned (3,:,:) to account for processor info)
    !
    ! iae2v,jae2v     :    Which vertices belong to each element


    ! calculate total elements, including boundary elements

    write(*,*) 'nelems     = ', nelems
    write(*,*) 'nelemzones = ', nelemzones
    write(*,*) 'nvertices  = ', nvertices

    ! allocate rows for vertex to element connectivity
    allocate(iav2e(nvertices+1))
    iav2e = 0

    ! loop through elements. Count the number of accesses for each vertex
    do izone = 1, nelemzones
      do k = boundaryelems(izone)%ielstart, boundaryelems(izone)%ielend
        do j = 1, boundaryelems(izone)%nbvpe
          iav2e(boundaryelems(izone)%belems(j,k)) = &
          iav2e(boundaryelems(izone)%belems(j,k)) + 1 
        end do
      end do
    end do

    ! allocate nonzeroes
    nnzv2e = sum(iav2e)
    allocate(v2e(2,nnzv2e))

    ! Count backwards such that row structure is correct upon populating v2e
    iav2e(1) = iav2e(1) + 1
    do i = 2, nvertices+1
      iav2e(i) = iav2e(i) + iav2e(i-1)
    end do
    do izone = 1, nelemzones
      do k = boundaryelems(izone)%ielstart, boundaryelems(izone)%ielend
        do j = 1, boundaryelems(izone)%nbvpe
          iav2e(boundaryelems(izone)%belems(j,k)) = &
          iav2e(boundaryelems(izone)%belems(j,k)) - 1
            v2e(:,iav2e(boundaryelems(izone)%belems(j,k))) = (/ k, izone /) 
        end do
      end do
    end do

    allocate(ivtmp1(nverticesperface),ivtmp2(nverticesperface))
    ! calculate element-to-element connectivity using shared nodes
    allocate(ef2e(2,2*ndim,1:nelems))
    ef2e = 0
    nnze2v = nverticesperelem*nelems
    allocate(iae2v(nelems+1))
    iae2v = 0
    allocate(jae2v(nnze2v))
    jae2v = 0

    connectioncount = 0
    k = 0
    vertexcount = 1
    do kzone = 1, nvolumesections
      izone = isectionvolume(kzone)

      do ielem = boundaryelems(izone)%ielstart, boundaryelems(izone)%ielend
        k = k+1
        iae2v(k) = vertexcount
        do j = 1, nverticesperelem
          jae2v(vertexcount) = boundaryelems(izone)%belems(j,k)
          vertexcount = vertexcount+1
        end do
        faceloop: do iface = 1, nfacesperelem
          do j = 1, nverticesperface
            ! populate vertices on face
            ivtmp1(j) = boundaryelems(izone)%belems(eltypfaces(j,iface),k)
          end do
          ivtmp1 = isort(ivtmp1,nverticesperface)
          ! search elements connected to each vertex until a match is found
          do j = 1, nverticesperface
            ! search elements connected to vertices and count number of shared vertices
            do i = iav2e(ivtmp1(j)), iav2e(ivtmp1(j)+1) - 1
              jzone = v2e(2,i)
              kk = v2e(1,i)
              ! don't check self!
              if (kk == k) cycle
              icount = 0
              ! check for every vertex
              do ii = 1,nverticesperface
                do jj = 1, boundaryelems(jzone)%nbvpe
                  if( ivtmp1(ii) == boundaryelems(jzone)%belems(jj,kk) ) then
                    icount = icount + 1
                    exit
                  end if
                end do
              end do
              ! check to see if nvertices per face match
              if (icount == nverticesperface) then
                ! element is connected to face
                ef2e(2,iface,k) = kk
                if (jzone /= izone) then
                  ! this is a boundary face
                  ef2e(1,iface,k) = -boundaryelems(jzone)%btype
                else
                  ! this is not a boundary face
                  connectioncount = connectioncount + 1
                  ! determine which face of kk I am connected to
                  do jface = 1, nfacesperelem
                    icount = 0
                    do jj = 1, nverticesperface
                      ivtmp2(jj) = boundaryelems(jzone)%belems(eltypfaces(jj,jface),kk)
                    end do
                    ivtmp2 = isort(ivtmp2,nverticesperface)
                    if (all(ivtmp1==ivtmp2)) then
                      ! this is the right face
                      ef2e(1,iface,k) = jface
                      exit
                    end if
                  end do
                end if
                cycle faceloop
              end if
            end do
          end do
          ! Search should always find a connection but didnt. Somethings wrong
          write(*,*) 'face search failed' ; write(*,*) izone, k, iface ; write(*,*) ivtmp1, ivtmp2
          stop
        end do faceloop

      end do
    end do

    iae2v(nelems+1) = vertexcount

    deallocate(ivtmp1,ivtmp2)

    deallocate(iav2e,v2e)


  end subroutine E2EConnectivity_cgns !   SERIAL Routine

  subroutine calculatepartitions()    !   SERIAL Routine

!   SERIAL Routine
    ! this subroutine contains the c-bindings for
    ! calling the metis library on the "grid", defined by
    ! the element-to-node connectivities read from the datafile.
!   SERIAL Routine

    use mpimod
    use referencevariables
    use variables, only: nnze2v, iae2v, iae2v_tmp, jae2v, jae2v_tmp, elempart
    use iso_c_binding, only: c_int, c_loc, c_ptr
    implicit none
    ! minimum number of vertices shared by connected elements
    integer :: ncommon
    ! number of partitions
    integer :: nparts
    ! error return
    integer(c_int) :: icerr
    ! c pointers for arrays being sent to metis call
    type(c_ptr) :: xadj, adj, eepart, nnpart
    ! these arrays are used on return from the metics call
    integer(c_int), allocatable, target :: epart(:), npart(:)
    ! these arrays contain the grid connectivity in CSR format
    integer(c_int), allocatable, target :: xadjtmp(:), jadjtmp(:)
    ! loop index
    integer :: i

    ! row pointer for element-to-vertex connectivity
    allocate(xadjtmp(0:nelems))
    xadjtmp(0:nelems) = iae2v(1:nelems+1)
    ! connected vertices (use c array indexing)
    allocate(jadjtmp(0:nnze2v-1))
    jadjtmp = jae2v-1

    ! allocate sizes for returned partition and vertex values
    ! on return from the c function these arrays contain
    ! the partition data
    allocate(epart(0:nelems-1))
    epart = 0
    allocate(npart(0:nvertices-1))
    npart = 0

    ! number of partitions is equivalent to number of processes
    nparts = nprocs
    ! minimum number of shared vertices for an element connection
    ncommon = nverticesperface

    if(nparts > 1) then

      ! assign c pointers
      xadj   = c_loc(xadjtmp(0))
       adj   = c_loc(jadjtmp(0))
      eepart = c_loc(epart(0))
      nnpart = c_loc(npart(0))

!      write(*,*) 'Before calling metis', xadjtmp(594)

      ! call c function that calls metis
      icerr = calcMetisPartitions(nelems, nvertices, nparts, xadj, adj, ncommon, &
        eepart, nnpart, nnze2v)

!        write(*,*) 'After if', xadjtmp(:)
    end if

    ! FIX FOR THE MEMORY PROBLEM
    ! ===========================
    iae2v = iae2v_tmp
    jae2v = jae2v_tmp
    deallocate(iae2v_tmp)
    deallocate(jae2v_tmp)

!    write(*,*) 'out of if', iae2v(595)
    
    ! assign the partition data to the global fortran array
    allocate(elempart(nelems))
    do i = 0, nelems-1
      elempart(i+1) = epart(i)
    end do

    ! Deallocate memory
    deallocate(epart,npart)
    deallocate(jadjtmp)
    deallocate(xadjtmp)

  end subroutine calculatepartitions

  subroutine TFI2D(xl2d,nk,pos)
    ! this subroutine applies TFI to a two-dimensional
    ! surface, assuming the edges and corners have been
    ! populated with data.

    implicit none
    integer, intent(in) :: nk
    real(wp), intent(inout) :: xl2d(3,nk,nk)
    real(wp), intent(in) :: pos(nk)

    integer :: i,j
    real(wp) :: r(2)

    do j = 1,nk
      do i = 1,nk
        r(1) = 0.5_wp*(pos(i)+1.0_wp)
        r(2) = 0.5_wp*(pos(j)+1.0_wp)
        xl2d(:,i,j) = (one-r(1))*xl2d(:,1,j) + r(1)*xl2d(:,nk,j) &
          &         + (one-r(2))*xl2d(:,i,1) + r(2)*xl2d(:,i,nk) &
          &         - (one-r(1))*(one-r(2))*xl2d(:, 1, 1) &
          &         - (one-r(1))*     r(2) *xl2d(:, 1,nk) &
          &         -      r(1) *(one-r(2))*xl2d(:,nk, 1) &
          &         -      r(1) *     r(2) *xl2d(:,nk,nk)
      end do
    end do

  end subroutine TFI2D

  subroutine TFI3D(xl3d,nk,pos)
    ! this subroutine applies TFI to a two-dimensional
    ! surface, assuming the edges and corners have been
    ! populated with data.

    implicit none
    integer, intent(in) :: nk
    real(wp), intent(inout) :: xl3d(3,nk,nk,nk)
    real(wp), intent(in) :: pos(nk)

    integer :: i,j,k
    real(wp) :: r(3)

    do k = 1,nk
      do j = 1,nk
        do i = 1,nk
          r(1) = 0.5_wp*(pos(i)+1.0_wp)
          r(2) = 0.5_wp*(pos(j)+1.0_wp)
          r(3) = 0.5_wp*(pos(k)+1.0_wp)
          xl3d(:,i,j,k) = &
            (1.0_wp-r(1))*(1.0_wp-r(2))*xl3d(:,1,1,k) &
            + (1.0_wp-r(1))*r(2)*xl3d(:,1,nk,k) &
            + r(1)*(1.0_wp-r(2))*xl3d(:,nk,1,k) &
            + r(1)*r(2)*xl3d(:,nk,nk,k) &
            + (1.0_wp-r(2))*(1.0_wp-r(3))*xl3d(:,i,1,1) &
            + (1.0_wp-r(2))*r(3)*xl3d(:,i,1,nk) &
            + r(2)*(1.0_wp-r(3))*xl3d(:,i,nk,1) &
            + r(2)*r(3)*xl3d(:,i,nk,nk) &
            + (1.0_wp-r(1))*(1.0_wp-r(3))*xl3d(:,1,j,1) &
            + (1.0_wp-r(1))*r(3)*xl3d(:,1,j,nk) &
            + r(1)*(1.0_wp-r(3))*xl3d(:,nk,j,1) &
            + r(1)*r(3)*xl3d(:,nk,j,nk) &
            - 2.0_wp*( (1.0_wp-r(1))*(1.0_wp-r(2))*(1.0_wp-r(3))*xl3d(:,1,1,1) &
            + r(1)*(1.0_wp-r(2))*(1.0_wp-r(3))*xl3d(:,nk,1,1) &
            + (1.0_wp-r(1))*r(2)*(1.0_wp-r(3))*xl3d(:,1,nk,1) &
            + (1.0_wp-r(1))*(1.0_wp-r(2))*r(3)*xl3d(:,1,1,nk) &
            + r(1)*r(2)*(1.0_wp-r(3))*xl3d(:,nk,nk,1) &
            + r(1)*(1.0_wp-r(2))*r(3)*xl3d(:,nk,1,nk) &
            + (1.0_wp-r(1))*r(2)*r(3)*xl3d(:,1,nk,nk) &
            + r(1)*r(2)*r(3)*xl3d(:,nk,nk,nk) )
        end do
      end do
    end do

  end subroutine TFI3D

  !============================================================================
  
  subroutine perturb_vertices_tg_vortex_1(p_scale)

    ! Load modules
    use referencevariables
    use variables, only: vx, e2v
    
    ! Nothing is implicitly defined
    implicit none
    
    real(wp), intent(in) :: p_scale
    integer :: ielem
    integer :: i, j
    real(wp) :: dx(3), dr
    integer,parameter :: seed = 86456
    real(wp) :: rand
    
    real(wp) :: diff_x, diff_y, diff_z
    real(wp), parameter :: toll = 1e-6
    real :: tmp
    
    continue
  
    call srand(seed)

    ! loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)

      do i = 1, nverticesperelem
        ! Compute the difference between each coordinate and pi.
        ! The domain of computation for the Taylor-Green vortex goes from 
        ! -pi to +pi in all the three directions
        diff_x = abs(abs(vx(1,e2v(i,ielem))) - pi)
        diff_y = abs(abs(vx(2,e2v(i,ielem))) - pi)
        diff_z = abs(abs(vx(3,e2v(i,ielem))) - pi )

        ! If the vertex lies on one of the boundary faces do not perturb it.
        ! This is necessary to periodic BC.
        if (diff_x .lt. toll .or. diff_y .lt. toll .or. diff_z .lt. toll) then
          cycle
        else
          dr = 1000.0_wp
          do j = 1, i-1
            dr = min( dr, magnitude( vx(:,e2v(j,ielem)) - vx(:,e2v(i,ielem)) ) )
          end do
          do j = i+1, nverticesperelem
            dr = min( dr, magnitude( vx(:,e2v(j,ielem)) - vx(:,e2v(i,ielem)) ) )
          end do
        
          tmp = rand(seed)
          dr = dr*p_scale*rand()
          dx = 0.0_wp
          dx(1:ndim) = dr
          vx(:,e2v(i,ielem)) = vx(:,e2v(i,ielem)) + dx
        end if
      end do
    end do

    return
  end subroutine perturb_vertices_tg_vortex_1

  !============================================================================

  subroutine pert_int_vert(p_scale)

    ! Load modules
    use referencevariables
    use variables, only: vx_master
    
    ! Nothing is implicitly defined
    implicit none
    
    real(wp), intent(in) :: p_scale
    integer :: j
    real(wp) :: dx(3)
    integer,parameter :: seed = 86456
    real(wp) :: rand
    
    integer :: n_tot_vertices, i_vertex

    real(wp), parameter :: toll = 1e-6
    real :: tmp
    real(wp) :: x, y, z
    
    continue
  
    call srand(seed)

    n_tot_vertices = size(vx_master(1,:))

    do i_vertex = 1, n_tot_vertices
      ! Compute the difference between each coordinate and pi.
      ! The domain of computation for the Taylor-Green vortex goes from 
      ! -pi to +pi in all the three directions
!      diff_x = abs(abs(vx_master(1,i_vertex)) - pi)
!      diff_y = abs(abs(vx_master(2,i_vertex)) - pi)
!      diff_z = abs(abs(vx_master(3,i_vertex)) - pi)

      x = abs(vx_master(1,i_vertex))
!      write(*,*) 'x', x
      y = abs(vx_master(2,i_vertex))
!      write(*,*) 'y', y
      z = abs(vx_master(3,i_vertex))
!      write(*,*) 'z', z

      ! If the vertex lies on one of the boundary faces do not perturb it.
      ! This is necessary to periodic BC.
!      if (diff_x .lt. toll .or. diff_y .lt. toll .or. diff_z .lt. toll) then
      if (abs(x-5.0_wp) .lt. toll) then
!        write(*,*) 'x', x
        cycle
      else if (abs(y-5.0_wp) .lt. toll) then
!        write(*,*) 'y', y
        cycle
      else if (abs(z-1.0_wp) .lt. toll) then
!        write(*,*) 'z', z
        cycle 
      else if (abs(z) .lt. toll) then
!        write(*,*) 'z', z
        cycle
      else
      !    dr = 1000.0_wp
      !    do j = 1, i-1
      !      dr = min( dr, magnitude( vx(:,e2v(j,ielem)) - vx(:,e2v(i,ielem)) ) )
      !    end do
      !    do j = i+1, nverticesperelem
      !      dr = min( dr, magnitude( vx(:,e2v(j,ielem)) - vx(:,e2v(i,ielem)) ) )
      !    end do
!       write(*,*) 'here'
        tmp = rand(seed)
!        dr = p_scale*rand()*vx_master(1,i_vertex)
        dx = 0.0_wp
        do j = 1, ndim
          dx(j) = p_scale*rand()*vx_master(j,i_vertex)
        end do
        vx_master(:,i_vertex) = vx_master(:,i_vertex) + dx
      end if
    end do

    return
  end subroutine pert_int_vert

  !============================================================================

  subroutine calcnodes_LGL()
    ! this subroutine populates the nodal grid locations for
    ! each element. Currently we assume that all elements are
    ! straight sided, but this can be remedied by incorporating
    ! CAD or analytical surface data. 
    use controlvariables, only: Grid_Topology, cylinder_x0, cylinder_x1, radius, origin, hrefine
    use referencevariables
    use variables, only: xg, vx, e2v, ef2e
    use initcollocation, only: element_properties, Gauss_Lobatto_Legendre_points
!-- DEBUG
    use collocationvariables, only : elem_props
!-- DEBUG

    implicit none
    ! indices
    integer :: ielem, inode, idir, iface
    integer :: i,j,k
    integer :: nE
    integer :: nodesperelem_max
    ! cartesian based grid coordinates
    real(wp), allocatable :: xl(:,:,:,:)

    real(wp), allocatable :: x_LGL_1d(:)

    ! high and low indices for each direction
    integer :: il(2,3)
    ! local grid distance
    real(wp)                :: dr
    real(wp), dimension(3)  :: dx
    real(wp), dimension(3)  :: x00,x01

    real(wp), dimension(:), allocatable  :: xi

    integer                                             :: i1d
    integer                                             :: nmin

    real(wp), parameter                                 :: tol = 1.0e-12_wp
    real(wp), allocatable                               :: x_LGL_1d_min(:), w_LGL_1d_min(:)

    real(wp), dimension(4,3)                            :: points_surf
    real(wp), dimension(4)                              :: r


    nodesperelem_max = (npoly_max+1)**ndim                     ! number of nodes in each element

    allocate(xg(3,1:nodesperelem_max,ihelems(1):ihelems(2)))   ! allocate global node matrix
    xg = 0.0_wp

    element_do : do ielem = ihelems(1), ihelems(2)             ! loop over volumetric elements

      call element_properties(ielem,       &                   !     ! nE is size of edge on element (varies with element)
                              n_pts_1d=nE, &
                              x_pts_1d=x_LGL_1d)
!write(*,*)"at top of calcnodes myprocid = ",myprocid," ielem = ",ielem, "nE&
!=",nE, "elem_props(1,ielem) = ",elem_props(1,ielem)
      if(allocated(xi)) deallocate(xi) ; allocate(xi(1:nE)) ;  xi = 0.0_wp

      if(allocated(xl)) deallocate(xl) ; allocate(xl(3,1:nE,1:nE,1:nE)) ;  xl = 0.0_wp

      ! low index is always 1
      il = 1
      ! set high index for each grid direction to nodesperedge
      do idir = 1,ndim
        il(2,idir) = nE
      end do

      ! reset local grid coordinates
      xl = 0.0_wp
      ! initialize corners to vertex values
        xl(:, 1, 1, 1) = vx(:,e2v(1,ielem))
        xl(:,nE, 1, 1) = vx(:,e2v(2,ielem))
      if (ndim > 1) then
        xl(:,nE,nE, 1) = vx(:,e2v(3,ielem))
        xl(:, 1,nE, 1) = vx(:,e2v(4,ielem))
      end if
      if (ndim > 2) then
        xl(:, 1, 1,nE) = vx(:,e2v(5,ielem))
        xl(:,nE, 1,nE) = vx(:,e2v(6,ielem))
        xl(:,nE,nE,nE) = vx(:,e2v(7,ielem))
        xl(:, 1,nE,nE) = vx(:,e2v(8,ielem))
      end if

      ! Initialize edges
      ! Account for Curved boundaries

      do iface = 1, nfacesperelem                 ! loop over 6 faces
        if (ef2e(1,iface,ielem) >= 0) then        ! Check if the face is a boundary
          cycle
        else
          
        endif
      enddo

      select case(Grid_Topology)

        case ('linear')

        ! Build the ``Bird cage'': 12 bounding edge connectors that define the Hexahedral Element
        do i = 1,nE                                 ! loop over nodes on edge
            dr = 0.5_wp*(x_LGL_1d(i)+1.0_wp)    ! distance in computational space
          if (ndim > 0) then
            dx = xl(:,nE, 1, 1)-xl(:, 1, 1, 1) ; xl(:, i, 1, 1) = xl(:, 1, 1, 1) + dr*dx ! xi_2 = 0, xi_3 = 0
          endif
          if (ndim > 1) then
            dx = xl(:,nE,nE, 1)-xl(:, 1,nE, 1) ; xl(:, i,nE, 1) = xl(:, 1,nE, 1) + dr*dx ! xi_2 = 1, xi_3 = 0
            dx = xl(:, 1,nE, 1)-xl(:, 1, 1, 1) ; xl(:, 1, i, 1) = xl(:, 1, 1, 1) + dr*dx ! xi_1 = 0, xi_3 = 0
            dx = xl(:,nE,nE, 1)-xl(:,nE, 1, 1) ; xl(:,nE, i, 1) = xl(:,nE, 1, 1) + dr*dx ! xi_1 = 1, xi_3 = 0
          end if
          if (ndim > 2) then
            dx = xl(:,nE, 1,nE)-xl(:, 1, 1,nE) ; xl(:, i, 1,nE) = xl(:, 1, 1,nE) + dr*dx ! xi_2 = 0, xi_3 = 1
            dx = xl(:,nE,nE,nE)-xl(:, 1,nE,nE) ; xl(:, i,nE,nE) = xl(:, 1,nE,nE) + dr*dx ! xi_2 = 1, xi_3 = 1
            dx = xl(:, 1,nE,nE)-xl(:, 1, 1,nE) ; xl(:, 1, i,nE) = xl(:, 1, 1,nE) + dr*dx ! xi_1 = 0, xi_3 = 1
            dx = xl(:,nE,nE,nE)-xl(:,nE, 1,nE) ; xl(:,nE, i,nE) = xl(:,nE, 1,nE) + dr*dx ! xi_1 = 1, xi_3 = 1
            dx = xl(:, 1, 1,nE)-xl(:, 1, 1, 1) ; xl(:, 1, 1, i) = xl(:, 1, 1, 1) + dr*dx ! xi_1 = 0, xi_2 = 0
            dx = xl(:, 1,nE,nE)-xl(:, 1,nE, 1) ; xl(:, 1,nE, i) = xl(:, 1,nE, 1) + dr*dx ! xi_1 = 0, xi_2 = 1
            dx = xl(:,nE,nE,nE)-xl(:,nE,nE, 1) ; xl(:,nE,nE, i) = xl(:,nE,nE, 1) + dr*dx ! xi_1 = 1, xi_2 = 1
            dx = xl(:,nE, 1,nE)-xl(:,nE, 1, 1) ; xl(:,nE, 1, i) = xl(:,nE, 1, 1) + dr*dx ! xi_1 = 1, xi_2 = 0
          end if
        end do
 
        case ('cylinder')

        x00 = cylinder_x0 ; x01 = cylinder_x1 ;

        xi(:) = 0.5_wp*(x_LGL_1d(:)+1.0_wp)    ! distance in computational space

        if (ndim > 0) then
          xl(:, :, 1, 1) = curved_connector_cylinder(nE,x00,x01,xl(:, 1, 1, 1),xl(:,nE, 1, 1),xi) ! xi_2 = 0, xi_3 = 0
        end if
        if (ndim > 1) then
          xl(:, :,nE, 1) = curved_connector_cylinder(nE,x00,x01,xl(:, 1,nE, 1),xl(:,nE,nE, 1),xi) ! xi_2 = 1, xi_3 = 0
          xl(:, 1, :, 1) = curved_connector_cylinder(nE,x00,x01,xl(:, 1, 1, 1),xl(:, 1,nE, 1),xi) ! xi_1 = 0, xi_3 = 0
          xl(:,nE, :, 1) = curved_connector_cylinder(nE,x00,x01,xl(:,nE, 1, 1),xl(:,nE,nE, 1),xi) ! xi_1 = 1, xi_3 = 0
        end if
        if (ndim > 2) then
          xl(:, :, 1,nE) = curved_connector_cylinder(nE,x00,x01,xl(:, 1, 1,nE),xl(:,nE, 1,nE),xi) ! xi_2 = 0, xi_3 = 1
          xl(:, :,nE,nE) = curved_connector_cylinder(nE,x00,x01,xl(:, 1,nE,nE),xl(:,nE,nE,nE),xi) ! xi_2 = 1, xi_3 = 1
          xl(:, 1, :,nE) = curved_connector_cylinder(nE,x00,x01,xl(:, 1, 1,nE),xl(:, 1,nE,nE),xi) ! xi_1 = 0, xi_3 = 1
          xl(:,nE, :,nE) = curved_connector_cylinder(nE,x00,x01,xl(:,nE, 1,nE),xl(:,nE,nE,nE),xi) ! xi_1 = 1, xi_3 = 1
          xl(:, 1, 1, :) = curved_connector_cylinder(nE,x00,x01,xl(:, 1, 1, 1),xl(:, 1, 1,nE),xi) ! xi_1 = 0, xi_2 = 0
          xl(:, 1,nE, :) = curved_connector_cylinder(nE,x00,x01,xl(:, 1,nE, 1),xl(:, 1,nE,nE),xi) ! xi_1 = 0, xi_2 = 1
          xl(:,nE,nE, :) = curved_connector_cylinder(nE,x00,x01,xl(:,nE,nE, 1),xl(:,nE,nE,nE),xi) ! xi_1 = 1, xi_2 = 1
          xl(:,nE, 1, :) = curved_connector_cylinder(nE,x00,x01,xl(:,nE, 1, 1),xl(:,nE, 1,nE),xi) ! xi_1 = 1, xi_2 = 0
        end if

        case ('parabola')

        xi(:) = 0.5_wp*(x_LGL_1d(:)+1.0_wp)    ! distance in computational space

        if (ndim > 0) then
          xl(:, :, 1, 1) = curved_connector_parabola(nE,xl(:, 1, 1, 1),xl(:,nE, 1, 1),xi) ! xi_2 = 0, xi_3 = 0
        end if
        if (ndim > 1) then
          xl(:, :,nE, 1) = curved_connector_parabola(nE,xl(:, 1,nE, 1),xl(:,nE,nE, 1),xi) ! xi_2 = 1, xi_3 = 0
          xl(:, 1, :, 1) = curved_connector_parabola(nE,xl(:, 1, 1, 1),xl(:, 1,nE, 1),xi) ! xi_1 = 0, xi_3 = 0
          xl(:,nE, :, 1) = curved_connector_parabola(nE,xl(:,nE, 1, 1),xl(:,nE,nE, 1),xi) ! xi_1 = 1, xi_3 = 0
        end if
        if (ndim > 2) then
          xl(:, :, 1,nE) = curved_connector_parabola(nE,xl(:, 1, 1,nE),xl(:,nE, 1,nE),xi) ! xi_2 = 0, xi_3 = 1
          xl(:, :,nE,nE) = curved_connector_parabola(nE,xl(:, 1,nE,nE),xl(:,nE,nE,nE),xi) ! xi_2 = 1, xi_3 = 1
          xl(:, 1, :,nE) = curved_connector_parabola(nE,xl(:, 1, 1,nE),xl(:, 1,nE,nE),xi) ! xi_1 = 0, xi_3 = 1
          xl(:,nE, :,nE) = curved_connector_parabola(nE,xl(:,nE, 1,nE),xl(:,nE,nE,nE),xi) ! xi_1 = 1, xi_3 = 1
          xl(:, 1, 1, :) = curved_connector_parabola(nE,xl(:, 1, 1, 1),xl(:, 1, 1,nE),xi) ! xi_1 = 0, xi_2 = 0
          xl(:, 1,nE, :) = curved_connector_parabola(nE,xl(:, 1,nE, 1),xl(:, 1,nE,nE),xi) ! xi_1 = 0, xi_2 = 1
          xl(:,nE,nE, :) = curved_connector_parabola(nE,xl(:,nE,nE, 1),xl(:,nE,nE,nE),xi) ! xi_1 = 1, xi_2 = 1
          xl(:,nE, 1, :) = curved_connector_parabola(nE,xl(:,nE, 1, 1),xl(:,nE, 1,nE),xi) ! xi_1 = 1, xi_2 = 0
        end if
      case('sphere')
        hrefine_if: if(hrefine)then
          if((ef2e(8,2,ielem).GE.1))then
          !-- refined element
          call curved_sphere_hrefine(nE,ielem,xl)
          else
            !-- populate the srufaces of the block
            call surface_nodes_sphere(nE,ielem,xl,.false.)
          endif
        else
          !-- populate the srufaces of the block
          call surface_nodes_sphere(nE,ielem,xl,.false.)
        endif hrefine_if
    end select

      ! build faces
      if(Grid_Topology.EQ.'sphere')then
      else
        if (ndim > 1) then
          ! xi_3 = 0
          call TFI2D(xl(:, :, :, 1),nE,x_LGL_1d)
        end if
        if (ndim > 2) then
          ! xi_3 = 1
          call TFI2D(xl(:, :, :,nE),nE,x_LGL_1d)
          ! xi_2 = 0
          call TFI2D(xl(:, :, 1, :),nE,x_LGL_1d)
          ! xi_2 = 1
          call TFI2D(xl(:, :,nE, :),nE,x_LGL_1d)
          ! xi_1 = 0
          call TFI2D(xl(:, 1, :, :),nE,x_LGL_1d)
          ! xi_1 = 1
          call TFI2D(xl(:,nE, :, :),nE,x_LGL_1d)
        end if
      endif

      ! build volumes
      if ((ndim > 2)) then
        if(hrefine.AND.(Grid_Topology.EQ.'sphere'))then
          if((ef2e(8,2,ielem).GE.1))then
            !-- do nothing the volume nodes have already been computed
          else
            call TFI3D(xl(:,:,:,:),nE,x_LGL_1d)
          endif
        else
          call TFI3D(xl(:,:,:,:),nE,x_LGL_1d)
        endif
      end if

      ! populate global coordinate matrix simply by packing
      ! in the typical manner
      inode = 0
      do k = il(1,3), il(2,3)
        do j = il(1,2), il(2,2)
          do i = il(1,1), il(2,1)
            inode = inode + 1
            xg(:,inode,ielem) = xl(:,i,j,k)
          end do
        end do
      end do
    end do element_do

    deallocate(xl)
    deallocate(xi)
    deallocate(x_LGL_1d)

  end subroutine calcnodes_LGL
  !================================================================================================
  !
  ! Purpose: This function populates the surface nodes of a hex for use on
  ! sphereical meshes 
  !
  ! inputs:
  !         nE = number of nodes in each direction
  !         ielem = element that is being populated
  !         xl = matrix that stores the nodal locations
  !
  ! Output:
  !         xl = populated with the nodal distribution
  !
  ! Notes:
  !================================================================================================
  subroutine surface_nodes_sphere(nE, ielem,xl,parent)

    use referencevariables, only: ndim, myprocid
    use variables, only: vx, e2v, ef2e, parent_geo
    use controlvariables, only: radius, origin
    use initcollocation, only: element_properties, Gauss_Lobatto_Legendre_points

!-- DEBUG
    use collocationvariables, only: elem_props
!-- DEBUG
    integer, intent(in) :: nE, ielem
    real(wp), intent(inout) :: xl(1:3,1:nE,1:nE,1:nE)
    logical, intent(in) :: parent

    integer :: iface, j, i1d, nmin, nE_temp   
    real(wp) :: xl_parent(1:3,1:nE,1:nE,1:nE)
    real(wp), dimension(4,3)                            :: points_surf
    real(wp), dimension(4)                              :: r
    real(wp), allocatable :: x_LGL_1d(:)
    real(wp)                :: dr
    real(wp), dimension(3)  :: dx
    real(wp), allocatable                               :: x_LGL_1d_min(:), w_LGL_1d_min(:)

    real(wp), parameter                                 :: tol = 1.0e-12_wp

!--debug vars
    integer :: i
    logical :: curv

    curv = .true.

    call element_properties(ielem,       &                   !     ! nE is size of edge on element (varies with element)
                              n_pts_1d=nE_temp, &
                              x_pts_1d=x_LGL_1d)
!write(*,*)"top of surf ielem ",ielem,"myprocid = ",myprocid," nE = ",nE, "nE_temp = ",nE_temp,&
!"elem_props(1,ielem) = ",elem_props(1,ielem)," elem_props(2,ielem) =",elem_props(2,ielem) 
    !-- build the nodal distribution on the parent element
    face_do : do iface = 1,2*ndim
      !-- determine the radius and ensure that all 4 points have the same radius 
      !-- these are stored in a counter-clockwise fashion relatlve to the coordinate axis
      ! for example
      !-- face 
      !xi_2    ^
      !        |
      !        |
      !        |
      !        ------> xi_1  
      !  
      ! r(4)__________r(3)
      !     |         |
      !     |         |
      !     |         |
      !     |         |
      ! r(1) ----------r(2)
  
      parent_if: if(parent)then
        if(iface.EQ.1)then
          points_surf(1,:) = parent_geo(1:3,1,ef2e(8,1,ielem)) 
          points_surf(2,:) = parent_geo(1:3,2,ef2e(8,1,ielem))
          points_surf(3,:) = parent_geo(1:3,3,ef2e(8,1,ielem))
          points_surf(4,:) = parent_geo(1:3,4,ef2e(8,1,ielem))
        elseif(iface.EQ.2)then
          points_surf(1,:) = parent_geo(1:3,1,ef2e(8,1,ielem))
          points_surf(2,:) = parent_geo(1:3,2,ef2e(8,1,ielem))
          points_surf(3,:) = parent_geo(1:3,6,ef2e(8,1,ielem))
          points_surf(4,:) = parent_geo(1:3,5,ef2e(8,1,ielem))
        elseif(iface.EQ.3)then
          points_surf(1,:) = parent_geo(1:3,2,ef2e(8,1,ielem))
          points_surf(2,:) = parent_geo(1:3,3,ef2e(8,1,ielem))
          points_surf(3,:) = parent_geo(1:3,7,ef2e(8,1,ielem))
          points_surf(4,:) = parent_geo(1:3,6,ef2e(8,1,ielem))
        elseif(iface.EQ.4)then
          points_surf(1,:) = parent_geo(1:3,4,ef2e(8,1,ielem))
          points_surf(2,:) = parent_geo(1:3,3,ef2e(8,1,ielem))
          points_surf(3,:) = parent_geo(1:3,7,ef2e(8,1,ielem))
          points_surf(4,:) = parent_geo(1:3,8,ef2e(8,1,ielem))
        elseif(iface.EQ.5)then
          points_surf(1,:) = parent_geo(1:3,1,ef2e(8,1,ielem))
          points_surf(2,:) = parent_geo(1:3,4,ef2e(8,1,ielem))
          points_surf(3,:) = parent_geo(1:3,8,ef2e(8,1,ielem))
          points_surf(4,:) = parent_geo(1:3,5,ef2e(8,1,ielem))
        elseif(iface.EQ.6)then
          points_surf(1,:) = parent_geo(1:3,5,ef2e(8,1,ielem))
          points_surf(2,:) = parent_geo(1:3,6,ef2e(8,1,ielem))
          points_surf(3,:) = parent_geo(1:3,7,ef2e(8,1,ielem))
          points_surf(4,:) = parent_geo(1:3,8,ef2e(8,1,ielem))
        endif
!--old code
!        if(iface.EQ.1)then
!          points_surf(1,:) = vx(:,ef2e(8,1,ielem))
!          points_surf(2,:) = vx(:,ef2e(8,2,ielem))
!          points_surf(3,:) = vx(:,ef2e(8,3,ielem))
!          points_surf(4,:) = vx(:,ef2e(8,4,ielem))
!        elseif(iface.EQ.2)then
!          points_surf(1,:) = vx(:,ef2e(8,1,ielem))
!          points_surf(2,:) = vx(:,ef2e(8,2,ielem))
!          points_surf(3,:) = vx(:,ef2e(8,6,ielem))
!          points_surf(4,:) = vx(:,ef2e(8,5,ielem))
!        elseif(iface.EQ.3)then
!          points_surf(1,:) = vx(:,ef2e(8,2,ielem))
!          points_surf(2,:) = vx(:,ef2e(8,3,ielem))
!          points_surf(3,:) = vx(:,ef2e(8,7,ielem))
!          points_surf(4,:) = vx(:,ef2e(8,6,ielem))
!        elseif(iface.EQ.4)then
!          points_surf(1,:) = vx(:,ef2e(8,4,ielem))
!          points_surf(2,:) = vx(:,ef2e(8,3,ielem))
!          points_surf(3,:) = vx(:,ef2e(8,7,ielem))
!          points_surf(4,:) = vx(:,ef2e(8,8,ielem))
!        elseif(iface.EQ.5)then
!          points_surf(1,:) = vx(:,ef2e(8,1,ielem))
!          points_surf(2,:) = vx(:,ef2e(8,4,ielem))
!          points_surf(3,:) = vx(:,ef2e(8,8,ielem))
!          points_surf(4,:) = vx(:,ef2e(8,5,ielem))
!        elseif(iface.EQ.6)then
!          points_surf(1,:) = vx(:,ef2e(8,5,ielem))
!          points_surf(2,:) = vx(:,ef2e(8,6,ielem))
!          points_surf(3,:) = vx(:,ef2e(8,7,ielem))
!          points_surf(4,:) = vx(:,ef2e(8,8,ielem))
!        endif
      else
        if(iface.EQ.1)then
          points_surf(1,:) = vx(:,e2v(1,ielem))
          points_surf(2,:) = vx(:,e2v(2,ielem))
          points_surf(3,:) = vx(:,e2v(3,ielem))
          points_surf(4,:) = vx(:,e2v(4,ielem))
        elseif(iface.EQ.2)then
          points_surf(1,:) = vx(:,e2v(1,ielem))
          points_surf(2,:) = vx(:,e2v(2,ielem))
          points_surf(3,:) = vx(:,e2v(6,ielem))
          points_surf(4,:) = vx(:,e2v(5,ielem))
        elseif(iface.EQ.3)then
          points_surf(1,:) = vx(:,e2v(2,ielem))
          points_surf(2,:) = vx(:,e2v(3,ielem))
          points_surf(3,:) = vx(:,e2v(7,ielem))
          points_surf(4,:) = vx(:,e2v(6,ielem))
        elseif(iface.EQ.4)then
          points_surf(1,:) = vx(:,e2v(4,ielem))
          points_surf(2,:) = vx(:,e2v(3,ielem))
          points_surf(3,:) = vx(:,e2v(7,ielem))
          points_surf(4,:) = vx(:,e2v(8,ielem))
        elseif(iface.EQ.5)then
          points_surf(1,:) = vx(:,e2v(1,ielem))
          points_surf(2,:) = vx(:,e2v(4,ielem))
          points_surf(3,:) = vx(:,e2v(8,ielem))
          points_surf(4,:) = vx(:,e2v(5,ielem))
        elseif(iface.EQ.6)then
          points_surf(1,:) = vx(:,e2v(5,ielem))
          points_surf(2,:) = vx(:,e2v(6,ielem))
          points_surf(3,:) = vx(:,e2v(7,ielem))
          points_surf(4,:) = vx(:,e2v(8,ielem))
        endif
      endif parent_if
      do j = 1,4
        r(j) = magnitude(points_surf(j,:)-origin(:)) 
      enddo

  !-- face 1
      face_if: if(iface.EQ.1)then
!       write(*,*)"myprocid",myprocid,"ielem = ",ielem,"nE = ",nE,"size(x_LGL_1d) = ",size(x_LGL_1d)
        do i1d = 1,nE                                 ! loop over nodes on edge
          dr = 0.5_wp*(x_LGL_1d(i1d)+1.0_wp)    ! distance in computational space
          dx = xl(:,nE, 1, 1)-xl(:, 1, 1, 1) ; xl(:, i1d, 1, 1) = xl(:, 1, 1, 1) + dr*dx ! xi_2 = 0, xi_3 = 0
          dx = xl(:,nE,nE, 1)-xl(:,nE, 1, 1) ; xl(:,nE, i1d, 1) = xl(:,nE, 1, 1) + dr*dx ! xi_1 = 1, xi_3 = 0
          dx = xl(:,nE,nE, 1)-xl(:, 1,nE, 1) ; xl(:, i1d,nE, 1) = xl(:, 1,nE, 1) + dr*dx ! xi_2 = 1, xi_3 = 0
          dx = xl(:, 1,nE, 1)-xl(:, 1, 1, 1) ; xl(:, 1, i1d, 1) = xl(:, 1, 1, 1) + dr*dx ! xi_1 = 0, xi_3 = 0
        enddo
  
if(curv)then 
        !-- check to see if the surface is on the sphere and if so move points onto sphere
        if ( (abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol)&
             .AND.(abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol) ) then
  
           call snap_surface_to_sphere(nE,x_LGL_1d,xl(1:3,1:nE,1:nE,1))
        else
          !-- check the four connectors to see if they lay on the sphere
          if( (abs(r(1)-r(2)).LE.tol) )then
            !-- connector at xi_2 = 0, xi_3 = 0
            if((abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,1)
            endif
          
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
 
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :, 1, 1) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(2,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin) 
          endif  
          if( (abs(r(2)-r(3)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_3 = 0
            if((abs(r(2)-radius).LE.tol).AND.(abs(r(3)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,2)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, nE, :, 1) = curved_connector_sphere(nE,nmin,points_surf(2,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(2),origin) 
          endif           
          if( (abs(r(3)-r(4)).LE.tol) )then
            !-- connector at xi_2 = 1, xi_3 = 0
            if((abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,3)
            endif
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
            xl(:, :,nE, 1) = curved_connector_sphere(nE,nmin,points_surf(4,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(3),origin)
          endif
          if( (abs(r(4)-r(1)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_3 = 0
            if((abs(r(4)-radius).LE.tol).AND.(abs(r(1)-radius).LE.tol))then
              !-- this connector is on the surface
               nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,4)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1, :, 1) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(4,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif 
          ! xi_3 = 0
          call TFI2D(xl(:, :, :, 1),nE,x_LGL_1d)
        endif
else
          ! xi_3 = 0
          call TFI2D(xl(:, :, :, 1),nE,x_LGL_1d)
endif
  !-- face 2   
      elseif(iface.EQ.2)then
        do i1d = 1,nE                                 ! loop over nodes on edge
          dr = 0.5_wp*(x_LGL_1d(i1d)+1.0_wp)    ! distance in computational space
          dx = xl(:,nE, 1, 1)-xl(:, 1, 1, 1) ; xl(:, i1d, 1, 1) = xl(:, 1, 1, 1) + dr*dx ! xi_2 = 0, xi_3 = 0
          dx = xl(:,nE, 1,nE)-xl(:,nE, 1, 1) ; xl(:,nE, 1, i1d) = xl(:,nE, 1, 1) + dr*dx ! xi_1 = 1, xi_2 = 0
          dx = xl(:,nE, 1,nE)-xl(:, 1, 1,nE) ; xl(:, i1d, 1,nE) = xl(:, 1, 1,nE) + dr*dx ! xi_2 = 0, xi_3 = 1
          dx = xl(:, 1, 1,nE)-xl(:, 1, 1, 1) ; xl(:, 1, 1, i1d) = xl(:, 1, 1, 1) + dr*dx ! xi_1 = 0, xi_2 = 0
        enddo
if(curv)then        
        !-- check to see if the surface is on the sphere and if so move points onto sphere
        if ( (abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol)&
             .AND.(abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol) ) then
          call snap_surface_to_sphere(nE,x_LGL_1d,xl(1:3,1:nE,1,1:nE))
        else
          !-- check to see if the connectors lay on the sphere
          if( (abs(r(1)-r(2)).LE.tol) )then
            !-- connector at xi_2 = 0, xi_3 = 0
            if((abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,1)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :, 1, 1) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(2,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif  
          if( (abs(r(2)-r(3)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_2 = 0
            if((abs(r(2)-radius).LE.tol).AND.(abs(r(3)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,2)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
            xl(:,nE, 1, :) = curved_connector_sphere(nE,nmin,points_surf(2,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(2),origin) 
          endif           
          if( (abs(r(3)-r(4)).LE.tol) )then
            !-- connector at xi_2 = 0, xi_3 = 1
            if((abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,3)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :, 1,nE) = curved_connector_sphere(nE,nmin,points_surf(4,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(3),origin)
          endif
          if( (abs(r(4)-r(1)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_2 = 0
            if((abs(r(4)-radius).LE.tol).AND.(abs(r(1)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,4)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1, 1, :) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(4,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif     
        
          ! xi_2 = 0
          call TFI2D(xl(:, :, 1, :),nE,x_LGL_1d)
        endif  
else
          ! xi_2 = 0
          call TFI2D(xl(:, :, 1, :),nE,x_LGL_1d)
endif  
  !-- face 3
      elseif(iface.EQ.3)then
        do i1d = 1,nE                                 ! loop over nodes on edge
          dr = 0.5_wp*(x_LGL_1d(i1d)+1.0_wp)    ! distance in computational space
          dx = xl(:,nE,nE, 1)-xl(:,nE, 1, 1) ; xl(:,nE, i1d, 1) = xl(:,nE, 1, 1) + dr*dx ! xi_1 = 1, xi_3 = 0
          dx = xl(:,nE,nE,nE)-xl(:,nE,nE, 1) ; xl(:,nE,nE, i1d) = xl(:,nE,nE, 1) + dr*dx ! xi_1 = 1, xi_2 = 1
          dx = xl(:,nE,nE,nE)-xl(:,nE, 1,nE) ; xl(:,nE, i1d,nE) = xl(:,nE, 1,nE) + dr*dx ! xi_1 = 1, xi_3 = 1
          dx = xl(:,nE, 1,nE)-xl(:,nE, 1, 1) ; xl(:,nE, 1, i1d) = xl(:,nE, 1, 1) + dr*dx ! xi_1 = 1, xi_2 = 0
        enddo
if(curv)then  
        !-- check to see if the surface is on the sphere and if so move points onto sphere
        if ( (abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol)&
             .AND.(abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol) ) then
          call snap_surface_to_sphere(nE,x_LGL_1d,xl(1:3,nE,1:nE,1:nE))
        else
          !-- check to see if the connectors lay on the sphere
          if( (abs(r(1)-r(2)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_3 = 0
            if((abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol))then 
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,1)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:,nE, :, 1) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(2,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif  
          if( (abs(r(2)-r(3)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_2 = 1
            if((abs(r(2)-radius).LE.tol).AND.(abs(r(3)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,2)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
        
            xl(:,nE,nE, :) = curved_connector_sphere(nE,nmin,points_surf(2,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(2),origin) 
          endif           
          if( (abs(r(3)-r(4)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_3 = 1
            if((abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,3)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:,nE, :,nE) = curved_connector_sphere(nE,nmin,points_surf(4,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(3),origin)
          endif
          if( (abs(r(4)-r(1)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_2 = 0
            if((abs(r(4)-radius).LE.tol).AND.(abs(r(1)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,4)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
        
            xl(:,nE, 1, :) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(4,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
           endif 
  
          ! xi_1 = 1
          call TFI2D(xl(:,nE, :, :),nE,x_LGL_1d)
        endif  
else
          ! xi_1 = 1
          call TFI2D(xl(:,nE, :, :),nE,x_LGL_1d)
endif
  !-- face 4
      elseif(iface.EQ.4)then
        do i1d = 1,nE                                 ! loop over nodes on edge
          dr = 0.5_wp*(x_LGL_1d(i1d)+1.0_wp)    ! distance in computational space
          dx = xl(:,nE,nE, 1)-xl(:, 1,nE, 1) ; xl(:, i1d,nE, 1) = xl(:, 1,nE, 1) + dr*dx ! xi_2 = 1, xi_3 = 0
          dx = xl(:,nE,nE,nE)-xl(:,nE,nE, 1) ; xl(:,nE,nE, i1d) = xl(:,nE,nE, 1) + dr*dx ! xi_1 = 1, xi_2 = 1
          dx = xl(:,nE,nE,nE)-xl(:, 1,nE,nE) ; xl(:, i1d,nE,nE) = xl(:, 1,nE,nE) + dr*dx ! xi_2 = 1, xi_3 = 1
          dx = xl(:, 1,nE,nE)-xl(:, 1,nE, 1) ; xl(:, 1,nE, i1d) = xl(:, 1,nE, 1) + dr*dx ! xi_1 = 0, xi_2 = 1
        enddo
if(curv)then  
        !-- check to see if the surface is on the sphere and if so move points onto sphere
        if ( (abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol)&
             .AND.(abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol) ) then
          call snap_surface_to_sphere(nE,x_LGL_1d,xl(1:3,1:nE,nE,1:nE))
        else
          !-- check to see if the connectors lay on the sphere
          if( (abs(r(1)-r(2)).LE.tol) )then
            !-- connector at xi_2 = 1, xi_3 = 0
            if((abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,1)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :,nE, 1) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(2,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif  
          if( (abs(r(2)-r(3)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_2 = 1
            if((abs(r(2)-radius).LE.tol).AND.(abs(r(3)-radius).LE.tol))then
              !-- this connector is on the surface            
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,2)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:,nE,nE, :) = curved_connector_sphere(nE,nmin,points_surf(2,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(2),origin) 
          endif           
          if( (abs(r(3)-r(4)).LE.tol) )then
            !-- connector at xi_2 = 1, xi_3 = 1
            if((abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,3)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :,nE,nE) = curved_connector_sphere(nE,nmin,points_surf(4,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(3),origin)
          endif
          if( (abs(r(4)-r(1)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_2 = 1
            if((abs(r(4)-radius).LE.tol).AND.(abs(r(1)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,4)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1,nE, :) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(4,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif
  
          ! xi_2 = 1
          call TFI2D(xl(:, :,nE, :),nE,x_LGL_1d)
  
        endif  
else
          ! xi_2 = 1
          call TFI2D(xl(:, :,nE, :),nE,x_LGL_1d)
endif  
  !-- face 5
      elseif(iface.EQ.5)then
        do i1d = 1,nE                                 ! loop over nodes on edge
          dr = 0.5_wp*(x_LGL_1d(i1d)+1.0_wp)    ! distance in computational space
          dx = xl(:, 1,nE, 1)-xl(:, 1, 1, 1) ; xl(:, 1, i1d, 1) = xl(:, 1, 1, 1) + dr*dx ! xi_1 = 0, xi_3 = 0
          dx = xl(:, 1,nE,nE)-xl(:, 1,nE, 1) ; xl(:, 1,nE, i1d) = xl(:, 1,nE, 1) + dr*dx ! xi_1 = 0, xi_2 = 1   
          dx = xl(:, 1,nE,nE)-xl(:, 1, 1,nE) ; xl(:, 1, i1d,nE) = xl(:, 1, 1,nE) + dr*dx ! xi_1 = 0, xi_3 = 1
          dx = xl(:, 1, 1,nE)-xl(:, 1, 1, 1) ; xl(:, 1, 1, i1d) = xl(:, 1, 1, 1) + dr*dx ! xi_1 = 0, xi_2 = 0
        enddo
 if(curv)then  
        !-- check to see if the surface is on the sphere and if so move points onto sphere
        if ( (abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol)&
             .AND.(abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol) ) then
          call snap_surface_to_sphere(nE,x_LGL_1d,xl(1:3,1,1:nE,1:nE))
        else
          !-- check to see if the connectors lay on the sphere
          if( (abs(r(1)-r(2)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_3 = 0
            if((abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,1)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1, :, 1) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(2,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif  
          if( (abs(r(2)-r(3)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_2 = 1
            if((abs(r(2)-radius).LE.tol).AND.(abs(r(3)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,2)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1,nE, :) = curved_connector_sphere(nE,nmin,points_surf(2,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(2),origin) 
          endif           
          if( (abs(r(3)-r(4)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_3 = 1
            if((abs(r(3)-radius).LE.toL).AND.(abs(r(4)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,3)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1, :,nE) = curved_connector_sphere(nE,nmin,points_surf(4,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(3),origin)
          endif
          if( (abs(r(4)-r(1)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_2 = 0
            if((abs(r(4)-radius).LE.tol).AND.(abs(r(1)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,4)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1, 1, :) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(4,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif
  
          ! xi_1 = 0
          call TFI2D(xl(:, 1, :, :),nE,x_LGL_1d)
        endif  
else
          ! xi_1 = 0
          call TFI2D(xl(:, 1, :, :),nE,x_LGL_1d)
endif  
  !-- face 6
      elseif(iface.EQ.6)then
        do i1d = 1,nE                                 ! loop over nodes on edge
          dr = 0.5_wp*(x_LGL_1d(i1d)+1.0_wp)    ! distance in computational space
          dx = xl(:,nE, 1,nE)-xl(:, 1, 1,nE) ; xl(:, i1d, 1,nE) = xl(:, 1, 1,nE) + dr*dx ! xi_2 = 0, xi_3 = 1
          dx = xl(:,nE,nE,nE)-xl(:,nE, 1,nE) ; xl(:,nE, i1d,nE) = xl(:,nE, 1,nE) + dr*dx ! xi_1 = 1, xi_3 = 1
          dx = xl(:,nE,nE,nE)-xl(:, 1,nE,nE) ; xl(:, i1d,nE,nE) = xl(:, 1,nE,nE) + dr*dx ! xi_2 = 1, xi_3 = 1
          dx = xl(:, 1,nE,nE)-xl(:, 1, 1,nE) ; xl(:, 1, i1d,nE) = xl(:, 1, 1,nE) + dr*dx ! xi_1 = 0, xi_3 = 1
        enddo
if(curv)then  
        !-- check to see if the surface is on the sphere and if so move points onto sphere
        if ( (abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol)&
             .AND.(abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol) ) then
          call snap_surface_to_sphere(nE,x_LGL_1d,xl(1:3,1:nE,1:nE,nE))
        else
          !-- check to see if the connectors lay on the sphere
          if( (abs(r(1)-r(2)).LE.tol) )then
            !-- connector at xi_2 = 0, xi_3 = 1
            if((abs(r(1)-radius).LE.tol).AND.(abs(r(2)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,1)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :, 1,nE) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(2,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif  
          if( (abs(r(2)-r(3)).LE.tol) )then
            !-- connector at xi_1 = 1, xi_3 = 1
            if((abs(r(2)-radius).LE.tol).AND.(abs(r(3)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,2)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:,nE, :,nE) = curved_connector_sphere(nE,nmin,points_surf(2,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(2),origin)
          endif           
          if( (abs(r(3)-r(4)).LE.tol) )then
            !-- connector at xi_2 = 1, xi_3 = 1
            if((abs(r(3)-radius).LE.tol).AND.(abs(r(4)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,3)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, :,nE,nE) = curved_connector_sphere(nE,nmin,points_surf(4,:),points_surf(3,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(3),origin)
          endif
          if( (abs(r(4)-r(1)).LE.tol) )then
            !-- connector at xi_1 = 0, xi_3 = 1
            if((abs(r(4)-radius).LE.tol).AND.(abs(r(1)-radius).LE.tol))then
              !-- this connector is on the surface
              nmin = nE
            else
              nmin = nmin_connector(ielem,iface,nE,4)
            endif
  
            if(allocated(x_LGL_1d_min)) deallocate(x_LGL_1d_min); allocate(x_LGL_1d_min(nmin)); x_LGL_1d_min = 0.0_wp;
            if(allocated(w_LGL_1d_min)) deallocate(w_LGL_1d_min); allocate(w_LGL_1d_min(nmin)); w_LGL_1d_min = 0.0_wp
  
            call Gauss_Lobatto_Legendre_points(nmin,x_LGL_1d_min,w_LGL_1d_min)  
      
            xl(:, 1, :,nE) = curved_connector_sphere(nE,nmin,points_surf(1,:),points_surf(4,:),&
                                                   x_LGL_1d,x_LGL_1d_min,r(1),origin)
          endif
  
          ! xi_3 = 1
          call TFI2D(xl(:, :, :,nE),nE,x_LGL_1d)
        endif
else
           ! xi_3 = 1
          call TFI2D(xl(:, :, :,nE),nE,x_LGL_1d) 
endif  
      endif face_if
    enddo face_do
  end subroutine surface_nodes_sphere
  !================================================================================================
  !
  ! Purpose: This function constructs the nodal distribution for a refined
  ! element based on the nodal distribution of the 
  ! inputs:
  !         nE = number of nodes in each direction
  !         xl = matrix that stores the nodal locations
  !
  ! Output:
  !         xl = populated with the nodal distribution
  !
  ! Notes:
  !================================================================================================
  subroutine curved_sphere_hrefine(nE,ielem,xl)

    use variables, only: ef2e, vx, parent_geo
    use referencevariables, only: ndim
    use initcollocation, only: element_properties, lagrange_basis_function_1d
    use non_conforming, only: Lagrange_interpolant_basis_1D

    integer, intent(in) :: nE, ielem
    real(wp), intent(inout) :: xl(1:3,1:nE,1:nE,1:nE)

    integer :: iface, i, j, k, ii, jj, kk, i1d,nE_temp   
    real(wp) :: xl_parent(1:3,1:nE,1:nE,1:nE)
    real(wp), dimension(nE) :: xi1d, eta1d, zeta1d
    real(wp) :: xiL, xiR, etaL, etaR, zetaL, zetaR, lxi, leta, lzeta
    real(wp), allocatable :: x_LGL_1d(:)

    call element_properties(ielem,       &                   !     ! nE is size of edge on element (varies with element)
                            n_pts_1d=nE_temp, &
                            x_pts_1d=x_LGL_1d)

    !-- populate parent element
    xl_parent = 0.0_wp
    ! set corner nodes
    xl_parent(:, 1, 1, 1) = parent_geo(1:3,1,ef2e(8,1,ielem))
    xl_parent(:,nE, 1, 1) = parent_geo(1:3,2,ef2e(8,1,ielem))
    xl_parent(:,nE,nE, 1) = parent_geo(1:3,3,ef2e(8,1,ielem))
    xl_parent(:, 1,nE, 1) = parent_geo(1:3,4,ef2e(8,1,ielem))
    xl_parent(:, 1, 1,nE) = parent_geo(1:3,5,ef2e(8,1,ielem))
    xl_parent(:,nE, 1,nE) = parent_geo(1:3,6,ef2e(8,1,ielem))
    xl_parent(:,nE,nE,nE) = parent_geo(1:3,7,ef2e(8,1,ielem))
    xl_parent(:, 1,nE,nE) = parent_geo(1:3,8,ef2e(8,1,ielem))

!-- old code
!    xl_parent(:, 1, 1, 1) = vx(:,ef2e(8,1,ielem))
!    xl_parent(:,nE, 1, 1) = vx(:,ef2e(8,2,ielem))
!    xl_parent(:,nE,nE, 1) = vx(:,ef2e(8,3,ielem))
!    xl_parent(:, 1,nE, 1) = vx(:,ef2e(8,4,ielem))
!    xl_parent(:, 1, 1,nE) = vx(:,ef2e(8,5,ielem))
!    xl_parent(:,nE, 1,nE) = vx(:,ef2e(8,6,ielem))
!    xl_parent(:,nE,nE,nE) = vx(:,ef2e(8,7,ielem))
!    xl_parent(:, 1,nE,nE) = vx(:,ef2e(8,8,ielem))

    call surface_nodes_sphere(nE,ielem,xl_parent,.true.) 

    ! build volumes
    call TFI3D(xl_parent(:,:,:,:),nE,x_LGL_1d)
    !-- populate the child block
    if(ef2e(8,2,ielem).EQ.1)then
      !-- domain [-1,0]X[-1,0]X[-1,0]
      xiL = -1.0_wp; xiR = 0.0_wp
      etaL = -1.0_wp; etaR = 0.0_wp
      zetaL = -1.0_wp; zetaR = 0.0_wp
    elseif(ef2e(8,2,ielem).EQ.2)then
      !-- domain [0,1]X[-1,0]X[-1,0]
      xiL = 0.0_wp; xiR = 1.0_wp
      etaL = -1.0_wp; etaR = 0.0_wp
      zetaL = -1.0_wp; zetaR = 0.0_wp
    elseif(ef2e(8,2,ielem).EQ.3)then
      !-- domain [0,1]X[0,1]X[-1,0]
      xiL = 0.0_wp; xiR = 1.0_wp
      etaL = 0.0_wp; etaR = 1.0_wp
      zetaL = -1.0_wp; zetaR = 0.0_wp
    elseif(ef2e(8,2,ielem).EQ.4)then
      !-- domain [-1,0]X[0,1]X[-1,0]
      xiL = -1.0_wp; xiR = 0.0_wp
      etaL = 0.0_wp; etaR = 1.0_wp
      zetaL = -1.0_wp; zetaR = 0.0_wp
    elseif(ef2e(8,2,ielem).EQ.5)then
      !-- domain [-1,0]X[-1,0]X[0,1]
      xiL = -1.0_wp; xiR = 0.0_wp
      etaL = -1.0_wp; etaR = 0.0_wp
      zetaL = 0.0_wp; zetaR = 1.0_wp
    elseif(ef2e(8,2,ielem).EQ.6)then
      !-- domain [0,1]X[-1,0]X[0,1]
      xiL = 0.0_wp; xiR = 1.0_wp
      etaL = -1.0_wp; etaR = 0.0_wp
      zetaL = 0.0_wp; zetaR = 1.0_wp
    elseif(ef2e(8,2,ielem).EQ.7)then
      !-- domain [0,1]X[0,1]X[0,1]
      xiL = 0.0_wp; xiR = 1.0_wp
      etaL = 0.0_wp; etaR = 1.0_wp
      zetaL = 0.0_wp; zetaR = 1.0_wp
    elseif(ef2e(8,2,ielem).EQ.8)then
      !-- domain [-1,0]X[0,1]X[0,1]
      xiL = -1.0_wp; xiR = 0.0_wp
      etaL = 0.0_wp; etaR = 1.0_wp
      zetaL = 0.0_wp; zetaR = 1.0_wp
    endif

    xi1d = (xiR-xiL)/2.0_wp*x_LGL_1d+(xiR+xiL)/2.0_wp
    eta1d = (etaR-etaL)/2.0_wp*x_LGL_1d+(etaR+etaL)/2.0_wp
    zeta1d = (zetaR-zetaL)/2.0_wp*x_LGL_1d+(zetaR+zetaL)/2.0_wp

    !-- loop over the three computational coordiantes and construct the nodes
    xl = 0.0_wp
    k_do: do k = 1,nE
      j_do: do j = 1,nE
        i_do: do i = 1,nE
          kk_do: do kk = 1,nE
            jj_do: do jj = 1 ,nE
              ii_do: do ii = 1, nE
                lxi = lagrange_basis_function_1d(xi1d(i), ii, x_LGL_1d, nE)
                leta = lagrange_basis_function_1d(eta1d(j), jj, x_LGL_1d, nE)
                lzeta = lagrange_basis_function_1d(zeta1d(k),kk, x_LGL_1d, nE)

                xl(1,i,j,k) = xl(1,i,j,k)+xl_parent(1,ii,jj,kk)*lxi*leta*lzeta
                xl(2,i,j,k) = xl(2,i,j,k)+xl_parent(2,ii,jj,kk)*lxi*leta*lzeta
                xl(3,i,j,k) = xl(3,i,j,k)+xl_parent(3,ii,jj,kk)*lxi*leta*lzeta
              enddo ii_do 
            enddo jj_do
          enddo kk_do
        enddo i_do
      enddo j_do
    enddo k_do
  end subroutine curved_sphere_hrefine
  !================================================================================================
  !
  ! Purpose: This function determines the minimum polynomial degree and hence number of nodes  of a 
  !          connector based on information in e_edge2e which is an array that contains information 
  !           on which elements touch a particular connector and what their polynomial degree is
  ! inputs:
  !         ielem, kelem, elem1 = element numbers
  !
  ! Output:
  !         nmin_connector = the number of nodes the interpolant through the connector should use
  !
  ! Notes:
  !================================================================================================
  function nmin_connector(ielem,iface,nE,connector)
    use variables, only: e_edge2e
    use controlvariables, only: SAT_type, hrefine
    use referencevariables, only: number_of_possible_partners
   
    implicit none
    integer, intent(in)    :: ielem, iface, nE, connector
    integer                :: nmin_connector, ipartner

   if(hrefine)then
     !-- do nothing, not using e_edge2e for now
     nmin_connector = nE
   else
     nmin_connector = nE
     do ipartner = 1,number_of_possible_partners
       if(e_edge2e(1,connector,ipartner,iface,ielem)>0)then
         nmin_connector = minval((/nmin_connector,e_edge2e(1,connector,ipartner,iface,ielem)/))
       endif
     enddo
   endif
   
   if((SAT_type.EQ."mod_SAT").OR.(hrefine))then
        nmin_connector = maxval((/floor((nmin_connector-1.0_wp)/2.0_wp),1/))+1
   endif
  end function nmin_connector
  !================================================================================================
  !
  ! Purpose: Constructs a patch of nodes on the sphere. It takes a straight patch with the nodal 
  !          locations of the four sides have been defined and then uses curved_connector_sphere
  !          to populate the interior nodes 
  !
  ! inputs:
  !         nE = number of points along the line segment to be computed
  !         xLGL = one dimensional nodal distribution in computational space
  !         xyz =  ending point of the circular arc
  !
  ! Output:
  !         xyz: a (3,nE,nE) matrix with the (x,y,z) locations of nodes of the patch on the sphere.
  !
  ! Notes:
  !================================================================================================
  subroutine snap_surface_to_sphere(nE,xLGL,xyz)
    use controlvariables, only:                     radius, origin
    implicit none
    integer, intent(in)                          :: nE
    real(wp), dimension(3,nE,nE), intent(inout)  :: xyz
    real(wp), dimension(nE),   intent(in) :: xLGL

    !-- local variables
    real(wp), dimension(3)                       :: x0_r, x1_r
    integer                                      :: i,j
    real(wp), parameter                          :: tol = 1.0e-12_wp

    do i = 1,nE 
      do j = 1,nE
        x0_r = xyz(1:3,1,j)
        x1_r = xyz(1:3,nE,j)
        xyz(1:3,i,j) = radius*(x0_r(1:3)+(x1_r(1:3)-x0_r(1:3))*(0.5_wp*xLGL(i)+0.5_wp))/&
                                       magnitude((x0_r(1:3)+(x1_r(1:3)-x0_r(1:3))*(0.5_wp*xLGL(i)+0.5_wp)))&
                                       +origin  
        if(abs(xyz(1,i,j)**2+xyz(2,i,j)**2+xyz(3,i,j)**2-radius**2).GE.tol)then
          write(*,*)'Error in initgrid: snap_surface_to_sphere point not on sphere with radius r = ',radius
        endif
      enddo
    enddo
  end subroutine snap_surface_to_sphere

  function curved_connector_sphere(nE,nmin,x0_vec,x1_vec,xLGL,xLGLmin,r,origin)
  !================================================================================================
  !
  ! Purpose: constructs the nodal distribution along a line segement on the sphere between two points
  !
  ! inputs:
  !         nE = number of points along the line segment to be computed
  !         nmin = number of points that will be used to construct the segment on the surface. These 
  !            points are then used to construct an interpolant and then this interpolant is evaluated 
  !            at the points xLGL to construct curved_connector_sphere
  !         x0_vec = starting point of the circular arc
  !         x1_vec =  ending point of the circular arc
  !         xLGL = Gauss Lobatto nodal distribution 
  !         r = radius of the sphere
  !         origin = origin of the sphere
  !
  ! Output:
  !       curved_connector_sphere: matrix containing the (x,y,z) positions of the interpolant 
  !         on the surface of the sphere constructed at the computational points xLGLmin, evaluated at the 
  !         computational points xLGL
  !
  ! Notes:
  !================================================================================================

    use initcollocation, only   : lagrange_basis_function_1d
    use referencevariables, only: ndim

    implicit none
    integer,                   intent(in) :: nE, nmin
    real(wp), dimension(ndim), intent(in) :: x0_vec, x1_vec
    real(wp), dimension(nE),   intent(in) :: xLGL, xLGLmin
    real(wp),                  intent(in) :: r
    real(wp), dimension(ndim), intent(in) :: origin

    !-- local variables
    real(wp), dimension(ndim)             :: x0_r, x1_r
    real(wp), parameter                   :: tol   = 1.0e-12_wp


    real(wp), dimension(ndim,nE)          :: curved_connector_sphere
    real(wp), dimension(ndim,nmin)        :: curved_connector_sphere_min

    integer                               :: i, j
    real(wp)                              :: xi, Lxi


    !-- position vectors of the first and last node relative to the origin of the sphere
    x0_r = x0_vec-origin
    x1_r = x1_vec-origin
    do i = 1,nmin
      curved_connector_sphere_min(1:3,i) = r*(x0_r(1:3)+(x1_r(1:3)-x0_r(1:3))*(0.5_wp*xLGLmin(i)+0.5_wp))/&
                                       magnitude((x0_r(1:3)+(x1_r(1:3)-x0_r(1:3))*(0.5_wp*xLGLmin(i)+0.5_wp)))&
                                       +origin 
    enddo

   !-- evaluate the interpolant through the points in curved_connector_sphere_temp at the 
   !   points xLGL in computational space
   curved_connector_sphere(1:3,1:nE) = 0.0_wp
   do i = 1,nE
     xi = xLGL(i)
     do j = 1,nmin
       Lxi = lagrange_basis_function_1d(xi, j, xLGLmin, nmin)
       curved_connector_sphere(1:3,i) = curved_connector_sphere(1:3,i)+Lxi*curved_connector_sphere_min(1:3,j)
     enddo
   enddo
  
  end function curved_connector_sphere
  
! =============================================================================

  subroutine facenodesetup_LGL_Driver()
    ! This subroutine calculates the partner node of each facial
    ! node on every volumetric element. Partner nodes of boundary
    ! faces are set to themselves.
    !  
    !   kfacenodes(nodesperface,nfacesperelem)  
    !      volumetric node index of face node  
    !      
    !   ifacenodes(nodesperface*nfacesperelem)  
    !      kfacenode flattened into a single vector
    !  
    use referencevariables, only: nfacesperelem
    use mpimod
    use variables, only: kfacenodes_LGL_p0, ifacenodes_LGL_p0 &
                       , kfacenodes_LGL_p1, ifacenodes_LGL_p1 &
                       , kfacenodes_LGL_p2, ifacenodes_LGL_p2 &
                       , kfacenodes, ifacenodes

    use collocationvariables, only: n_LGL_1d_p0, n_LGL_2d_p0 &
                                  , n_LGL_1d_p1, n_LGL_2d_p1 &
                                  , n_LGL_1d_p2, n_LGL_2d_p2 

    implicit none

    ! kfacenodes separates each face
    ! ifacenodes includes all faces

    allocate(kfacenodes_LGL_p0(n_LGL_2d_p0,nfacesperelem))
    allocate(ifacenodes_LGL_p0(n_LGL_2d_p0*nfacesperelem))

    allocate(kfacenodes_LGL_p1(n_LGL_2d_p1,nfacesperelem))
    allocate(ifacenodes_LGL_p1(n_LGL_2d_p1*nfacesperelem))

    allocate(kfacenodes_LGL_p2(n_LGL_2d_p2,nfacesperelem))
    allocate(ifacenodes_LGL_p2(n_LGL_2d_p2*nfacesperelem))

    call facenodesetup_LGL(n_LGL_1d_p0, n_LGL_2d_p0, kfacenodes_LGL_p0, ifacenodes_LGL_p0)
    call facenodesetup_LGL(n_LGL_1d_p1, n_LGL_2d_p1, kfacenodes_LGL_p1, ifacenodes_LGL_p1)
    call facenodesetup_LGL(n_LGL_1d_p2, n_LGL_2d_p2, kfacenodes_LGL_p2, ifacenodes_LGL_p2)

    allocate(kfacenodes(n_LGL_2d_p0,nfacesperelem))
    allocate(ifacenodes(n_LGL_2d_p0*nfacesperelem))

    kfacenodes(:,:) = kfacenodes_LGL_p0(:,:)
    ifacenodes(:)   = ifacenodes_LGL_p0(:)

    end subroutine facenodesetup_LGL_Driver

! =============================================================================

  subroutine facenodesetup_Gau_Driver()
    ! This subroutine calculates the partner node of each facial
    ! node on every volumetric element. Partner nodes of boundary
    ! faces are set to themselves.
    !  
    !   kfacenodes(nodesperface,nfacesperelem)  
    !      volumetric node index of face node  
    !      
    !   ifacenodes(nodesperface*nfacesperelem)  
    !      kfacenode flattened into a single vector
    !  
    use referencevariables, only: nfacesperelem
    use mpimod
    use variables, only: kfacenodes_Gau_p0, ifacenodes_Gau_p0 &
                       , kfacenodes_Gau_p1, ifacenodes_Gau_p1 &
                       , kfacenodes_Gau_p2, ifacenodes_Gau_p2 

    use collocationvariables, only: n_Gau_1d_p0, n_Gau_2d_p0 &
                                  , n_Gau_1d_p1, n_Gau_2d_p1 &
                                  , n_Gau_1d_p2, n_Gau_2d_p2 

    implicit none

    ! kfacenodes separates each face
    ! ifacenodes includes all faces

    allocate(kfacenodes_Gau_p0(n_Gau_2d_p0,nfacesperelem))
    allocate(ifacenodes_Gau_p0(n_Gau_2d_p0*nfacesperelem))

    allocate(kfacenodes_Gau_p1(n_Gau_2d_p1,nfacesperelem))
    allocate(ifacenodes_Gau_p1(n_Gau_2d_p1*nfacesperelem))

    allocate(kfacenodes_Gau_p2(n_Gau_2d_p2,nfacesperelem))
    allocate(ifacenodes_Gau_p2(n_Gau_2d_p2*nfacesperelem))

    call facenodesetup_Gau(n_Gau_1d_p0, n_Gau_2d_p0, kfacenodes_Gau_p0, ifacenodes_Gau_p0)
    call facenodesetup_Gau(n_Gau_1d_p1, n_Gau_2d_p1, kfacenodes_Gau_p1, ifacenodes_Gau_p1)
    call facenodesetup_Gau(n_Gau_1d_p2, n_Gau_2d_p2, kfacenodes_Gau_p2, ifacenodes_Gau_p2)

    end subroutine facenodesetup_Gau_Driver

! =============================================================================

  subroutine facenodesetup_LGL(n_LGL_1d, n_LGL_2d, kfacenodes, ifacenodes)

    ! This subroutine calculates the partner node of each facial
    ! node on every volumetric element. Partner nodes of boundary
    ! faces are set to themselves.
    !  
    !   kfacenodes(nodesperface,nfacesperelem)  
    !      volumetric node index of face node  
    !      
    !   ifacenodes(nodesperface*nfacesperelem)  
    !      kfacenode flattened into a single vector
    !  
    use referencevariables, only: nfacesperelem, ndim

    implicit none

    ! indices
    integer,                   intent(in   ) :: n_LGL_1d, n_LGL_2d
    integer,  dimension(:,:),  intent(inout) :: kfacenodes
    integer,  dimension(:  ),  intent(inout) :: ifacenodes

    integer :: i,j,k
    integer :: stride, stride1, stride2, ioffset

    real(wp), parameter :: nodetol = 1.0e-8_wp

    if (ndim == 2) then
      ! loop over every node on each face
      do i = 1, n_LGL_1d
        ! on face 1, the first n_LGL_2d nodes are just
        ! the first n_LGL_2d
        kfacenodes(i,1) = i
        ! on face 3 there is just an offset to where the
        ! counting starts
        j = (n_LGL_1d-1)*n_LGL_1d+i
        kfacenodes(i,3) = j
        ! onface 2, a stride and offset are required
        stride = n_LGL_1d
        j = n_LGL_1d + (i-1)*stride
        kfacenodes(i,2) = j
        ! on face 4, a stride is needed
        j = 1 + (i-1)*stride
        kfacenodes(i,4) = j 
      end do
    else if (ndim == 3) then
      k = 0
      do j = 1, n_LGL_1d
        do i = 1, n_LGL_1d
          k = k+1
          ! face 1 does not require an offset or a stride
          kfacenodes(k,1) = k
          ! on face 2, a stride is required
          ioffset = 1
          stride1 = 1
          stride2 = n_LGL_1d**2
          kfacenodes(k,2) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 3, offset and stride are needed
          ioffset = n_LGL_1d
          stride1 = n_LGL_1d
          stride2 = n_LGL_1d**2
          kfacenodes(k,3) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! face 4 requires an offset and a stride
          ioffset = 1+(n_LGL_1d-1)*n_LGL_1d
          stride1 = 1
          stride2 = n_LGL_1d**2
          kfacenodes(k,4) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 5 only a stride is required
          ioffset = 1
          stride1 = n_LGL_1d
          stride2 = n_LGL_1d**2
          kfacenodes(k,5) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 6 only an offset is required
          ioffset = (n_LGL_1d-1)*n_LGL_1d*n_LGL_1d
          kfacenodes(k,6) = ioffset+k
        end do
      end do
    else
      write(*,*) 'error: unsupported dimension', ndim
      stop
    end if
    ! the nodes on each face are packed into a 1D array
    k = 0
    ! loop over faces
    do j = 1, nfacesperelem
      ! loop over nodes on each face
      do i = 1, n_LGL_2d
        ! advance facial node index
        k = k+1
        ! map facial node index to volumetric node
        ifacenodes(k) = kfacenodes(i,j)
      end do
    end do

  end subroutine facenodesetup_LGL

! =============================================================================

  subroutine facenodesetup_LGL_WENO()
    ! This subroutine calculates the partner node of each facial
    ! node on every volumetric element. Partner nodes of boundary
    ! faces are set to themselves.
    !  
    !   kfacenodesWENO(nodesperface,nfacesperelem)  
    !      volumetric node index of face node  
    !      
    !   ifacenodesWENO(nodesperface*nfacesperelem)  
    !      kfacenode flattened into a single vector
    !  
    use referencevariables
    use mpimod
    use variables, only: kfacenodesWENO, ifacenodesWENO
    implicit none

    ! indices
    integer :: i,j,k
    integer :: stride, stride1, stride2, ioffset

    real(wp), parameter :: nodetol = 1.0e-8_wp

    ! local facial masks
    !
    ! kfacenodesWENO separates each face
    allocate(kfacenodesWENO(nodesperface,nfacesperelem))
    ! ifacenodesWENO includes all faces
    allocate(ifacenodesWENO(nodesperface*nfacesperelem))

    if (ndim == 2) then
      ! loop over every node on each face
      do i = 1, nodesperedge
        ! on face 1, the first nodesperface nodes are just
        ! the first nodesperface
        kfacenodesWENO(i,1) = i + nodesperedge
        ! on face 3 there is just an offset to where the
        ! counting starts
        j = (nodesperedge-2)*nodesperedge+i
        kfacenodesWENO(i,3) = j
        ! onface 2, a stride and offset are required
        stride = nodesperedge
        j = nodesperedge-1 + (i-1)*stride
        kfacenodesWENO(i,2) = j
        ! on face 4, a stride is needed
        j = 2 + (i-1)*stride
        kfacenodesWENO(i,4) = j 
      end do
    else if (ndim == 3) then
      k = 0
      do j = 1, nodesperedge
        do i = 1, nodesperedge
          k = k+1
          ! face 1 does not require an offset or a stride
          kfacenodesWENO(k,1) = k + nodesperface
          ! on face 2, a stride is required
          ioffset = nodesperedge + 1
          stride1 = 1
          stride2 = nodesperedge**2
          kfacenodesWENO(k,2) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 3, offset and stride are needed
          ioffset = nodesperedge - 1
          stride1 = nodesperedge
          stride2 = nodesperedge**2
          kfacenodesWENO(k,3) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! face 4 requires an offset and a stride
          ioffset = 1+(nodesperedge-2)*nodesperedge
          stride1 = 1
          stride2 = nodesperedge**2
          kfacenodesWENO(k,4) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 5 only a stride is required
          ioffset = 1 + 1
          stride1 = nodesperedge
          stride2 = nodesperedge**2
          kfacenodesWENO(k,5) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 6 only an offset is required
          ioffset = (nodesperedge-2)*nodesperedge*nodesperedge
          kfacenodesWENO(k,6) = ioffset+k
        end do
      end do
      !  Checked for accuracy:  MHC 09/02/2015
    else
      write(*,*) 'error: unsupported dimension', ndim
      stop
    end if
    ! the nodes on each face are packed into a 1D array
    k = 0
    ! loop over faces
    do j = 1, nfacesperelem
      ! loop over nodes on each face
      do i = 1, nodesperface
        ! advance facial node index
        k = k+1
        ! map facial node index to volumetric node
        ifacenodesWENO(k) = kfacenodesWENO(i,j)
      end do
    end do

  end subroutine facenodesetup_LGL_WENO

  !============================================================================
  ! calculate_face_node_connectivity - Sets the face-node connectivity for the
  ! collocation points.
  !============================================================================

  subroutine calculate_face_node_connectivity_LGL()
    
    ! Load modules
    use referencevariables
    use mpimod
    use variables, only: xg, xghst_LGL, ef2e, efn2efn,        &
      & jelems, periodic_elem_face_ids_x1,                    &
      & periodic_elem_face_ids_x2, periodic_elem_face_ids_x3
    use collocationvariables, only: elem_props
    use initcollocation,      only: element_properties

    ! Nothing is implicitly defined
    implicit none

    integer, allocatable, dimension(:,:) :: kfacenodes
    integer, allocatable, dimension(:)   :: ifacenodes

    integer ::  ielem, inode, jnode, iface, knode
    integer ::  i_low

    real(wp) :: x1(3), x2(3)
    real(wp), parameter :: nodetol = 1.0e-8_wp

    integer :: i_p_face, p_dir
    logical :: match_found, conforming_interface
    real(wp), dimension(2) :: x1_p, x2_p

    integer :: cnt_debug, ii
    integer :: n_LGL_1d, n_LGL_2d, nodesperface_max
    integer :: kelem, n_LGL_2d_Off

    continue

    cnt_debug = 0

    ! efn2efn contains the partner node information of every facenode in the domain

    nodesperface_max = (npoly_max+1)**(ndim-1)

    !-- old
    !allocate(efn2efn(4,nfacesperelem*nodesperface_max,ihelems(1):ihelems(2))) ; efn2efn = -1000 ;
    allocate(efn2efn(4,(2*ndim)*nodesperface_max,ihelems(1):ihelems(2))) ; efn2efn = -1000 ;

    ! Initialize position of the ghost point in the stack
    i_low = 0

    ! loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)
      
      call element_properties(ielem,&
                              n_pts_1d=n_LGL_1d,  &
                              n_pts_2d=n_LGL_2d,  &
                            kfacenodes=kfacenodes,&
                            ifacenodes=ifacenodes )

      ! Reset facial node index counter
      knode = 0
      
      ! Loop over faces
!-- old
!      do iface = 1, nfacesperelem
      do iface = 1, 2*ndim

        kelem = ef2e(2,iface,ielem)  
        call element_properties(kelem, n_pts_2d=n_LGL_2d_Off)

        knode = n_LGL_2d * (iface - 1)

        ! If on boundary, connect to self
        if (ef2e(1,iface,ielem) < 0) then
          
          ! Loop over nodes on the boundary face
          do inode = 1, n_LGL_2d
            
            ! Update facial node index counter
            knode = knode + 1
            
            ! The first index is the volumetric node index and the second index is the element index
            efn2efn(:,knode,ielem) = (/ ifacenodes(knode), ielem, 0 , 0 /)
          
          end do

        ! ==========================================================================================
        !                           Parallel       .and.               Conforming interface
        ! ==========================================================================================
!        else if ((ef2e(3,iface,ielem) /= myprocid) .and. (ef2e(4,iface,ielem) == elem_props(2,ielem))) then
        else if ((ef2e(3,iface,ielem) /= myprocid) .and. (ef2e(4,iface,ielem) == elem_props(2,ielem)) &
                 .and. (ef2e(9,iface,ielem) == 0)) then

          !  ==================================================
          !  Parallel periodic logic in x1, x2, x3 directions
          !  ==================================================

          match_found = .false.                                                      ! Initialize match_found
          
          if (size(periodic_elem_face_ids_x1(1,:)) /= 0) then                        ! Loop through the elements that owns a periodic face in the x1 direction

            do i_p_face = 1, size(periodic_elem_face_ids_x1(1,:))                    ! Check if the ielem owns a periodic face and if iface is a periodic face

              if (periodic_elem_face_ids_x1(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x1(2,i_p_face) == iface) then

                match_found = .true.                                                 ! There is a match: change logical value of match_found

                p_dir = periodic_elem_face_ids_x1(3,i_p_face)                        ! Get the direction of "periodicity"

                do inode = 1, n_LGL_2d                                               ! Loop over the nodes on the On-Element face
                  
                  knode = knode + 1                                                  ! Update the facial node index counter
                  
                  x1 = xg(:,ifacenodes(knode),ielem)                                 ! Save the coordinates of the facial node
                
                  x1_p(:) = Extract_Parallel_Invariant(p_dir,x1)                     !   Extact the two invariant coordinate locations between parallel faces

                  do jnode = 1, n_LGL_2d                                             ! Loop over the nodes on the Off-Element face
                    
                    x2 = xghst_LGL(:,i_low + jnode)                                  ! Coordinates of the jnode
                
                    x2_p(:) = Extract_Parallel_Invariant(p_dir,x2)                   !   Extact the two invariant coordinate locations between parallel faces

                    if (magnitude(x1_p-x2_p) <= nodetol) then                        ! Check distance between the two nodes
                      
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem)) ! Set the volumetric node index of the connected node; ef2e(2) gives the element of the neighbor

                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                   ! Set the element of the connected node

                      efn2efn(3,knode,ielem) = i_low + jnode                         ! Set the node index in the ghost array

                      exit                                                           ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do                                                             ! End do jnode
                
                end do                                                               ! End do inode 

                i_low = i_low + n_LGL_2d                                             ! Update the position in the ghost stack

              end if                                                                 ! End if match found

              if (match_found .eqv. .true.) exit                                     ! If a partner face has been found exit from the loop over the elements that own a periodic face

            end do                                                                   ! End do loop over the elements that own a periodic face

          end if                                                                     ! End if check periodic face in x1 direction

          ! Loop through the elements that owns a periodic face in the x2 direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x2(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x2(1,:))

              if (periodic_elem_face_ids_x2(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x2(2,i_p_face) == iface) then

                ! There is a match: change logical value of match_found
                match_found = .true.

                ! Get the direction of "periodicity"
                p_dir = periodic_elem_face_ids_x2(3,i_p_face)

                do inode = 1, n_LGL_2d                                               ! Loop over the nodes on the On-Element face
                  
                  knode = knode + 1                                                  ! Update the facial node index counter
                  
                  x1 = xg(:,ifacenodes(knode),ielem)                                 ! Save the coordinates of the facial node
                
                  x1_p(:) = Extract_Parallel_Invariant(p_dir,x1)                     !   Extact the two invariant coordinate locations between parallel faces

                  do jnode = 1, n_LGL_2d                                             ! Loop over the nodes on the Off-Element face
                    
                    x2 = xghst_LGL(:,i_low + jnode)                                  ! Coordinates of the jnode
                
                    x2_p(:) = Extract_Parallel_Invariant(p_dir,x2)                   !   Extact the two invariant coordinate locations between parallel faces

                    if (magnitude(x1_p-x2_p) <= nodetol) then                        ! Check distance between the two nodes
                      
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem)) ! Set the volumetric node index of the connected node; ef2e(2) gives the element of the neighbor

                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                   ! Set the element of the connected node

                      efn2efn(3,knode,ielem) = i_low + jnode                         ! Set the node index in the ghost array

                      exit                                                           ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do                                                             ! End do jnode
                
                end do                                                               ! End do inode 

                i_low = i_low + n_LGL_2d                                             ! Update the position in the ghost stack

              end if                                                                 ! End if match found

              if (match_found .eqv. .true.) exit                                     ! If a partner face has been found exit from the loop over the ! elements that own a periodic face

            end do                                                                   ! End do loop over the elements that own a periodic face

          end if                                                                     ! End if check periodic face in x2 direction

          ! Loop through the elements that owns a periodic face in the x3 direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x3(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x3(1,:))

              if (periodic_elem_face_ids_x3(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x3(2,i_p_face) == iface) then

                match_found = .true.                                                 ! There is a match: change logical value of match_found

                p_dir = periodic_elem_face_ids_x3(3,i_p_face)                        ! Get the direction of "periodicity"

                do inode = 1, n_LGL_2d                                               ! Loop over the nodes on the On-Element face
                  
                  knode = knode + 1                                                  ! Update the facial node index counter
                  
                  x1 = xg(:,ifacenodes(knode),ielem)                                 ! Save the coordinates of the facial node
                
                  x1_p(:) = Extract_Parallel_Invariant(p_dir,x1)                     !   Extact the two invariant coordinate locations between parallel faces

                  do jnode = 1, n_LGL_2d                                             ! Loop over the nodes on the Off-Element face
                    
                    x2 = xghst_LGL(:,i_low + jnode)                                  ! Coordinates of the jnode
                
                    x2_p(:) = Extract_Parallel_Invariant(p_dir,x2)                   !   Extact the two invariant coordinate locations between parallel faces

                    if (magnitude(x1_p-x2_p) <= nodetol) then                        ! Check distance between the two nodes
                      
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem)) ! Set the volumetric node index of the connected node; ef2e(2) gives the element of the neighbor

                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                   ! Set the element of the connected node

                      efn2efn(3,knode,ielem) = i_low + jnode                         ! Set the node index in the ghost array

                      exit                                                           ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do                                                             ! End do jnode
                
                end do                                                               ! End do inode 

                i_low = i_low + n_LGL_2d                                             ! Update the position in the ghost stack

              end if                                                                 ! End if match found

              if (match_found .eqv. .true.) exit                                     ! If a partner face has been found exit from the loop over the ! elements that own a periodic face

            end do                                                                   ! End do loop over the elements that own a periodic face

          end if                                                                     ! End if check periodic face in x3 direction
 
          !  ====================================================================
          !  Finished Parallel periodic logic:  Begin Parallel non-periodic logic
          !  ====================================================================

          if (match_found .eqv. .false.) then

            do inode = 1, n_LGL_2d                                                   ! Loop over the nodes on the face

              knode = knode + 1                                                      ! Update the facial node index counter

              x1 = xg(:,ifacenodes(knode),ielem)                                     ! Save the coordinates of the facial node
              
              do jnode = 1, n_LGL_2d                                                 ! Loop over the nodes on the Off-Element face

                x2 = xghst_LGL(:,i_low + jnode)                                      ! Coordinates of the jnode
                
                if (magnitude(x1-x2) <= nodetol) then                                ! Check the distance between the two nodes

                  cnt_debug = cnt_debug + 1                                          ! debugging the counter index
                  
                  efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))     ! Set the volumetric node index of the connected node, ef2e(2) gives the element of the neighbor
                  
                  efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                       ! Set the element of the connected node
                  
                  efn2efn(3,knode,ielem) = i_low + jnode                             ! Set the node index in the ghost array

                  exit                                                               ! Found a match, exit the loop
                
                end if             
              
              end do
             
              if (jnode > n_LGL_2d .and. myprocid==0) then                           ! Print information at screen if there is a problem and stop computation
                write(*,*) 'Connectivity error in face-node connectivity_LGL Parallel.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
                write(*,*) 'Node coordinates'
                write(*,*) x1(:)
                write(*,*) 'ghost node coordinates'

                do ii = 1,size(xghst_LGL,2)
!               do ii = i_low+1,i_low+n_LGL_2d
                  write(*,*)ii,xghst_LGL(:,ii)
                enddo
                write(*,*) 'Exiting...'
                stop
              end if

            end do

            i_low = i_low + n_LGL_2d                                                 ! Update the position in the ghost stack
          
          end if

        ! ==========================================================================================
        !                           Serial         .and.                  Conforming interface
        ! ==========================================================================================
!        else if ((ef2e(3,iface,ielem) == myprocid) .and. (ef2e(4,iface,ielem) == elem_props(2,ielem))) then
        else if ((ef2e(3,iface,ielem) == myprocid) .and. (ef2e(4,iface,ielem) == elem_props(2,ielem)) &
                 .and. (ef2e(9,iface,ielem) == 0) ) then

          match_found = .false.                                                      ! Initialize match_found

          if (size(periodic_elem_face_ids_x1(1,:)) /= 0) then

            do i_p_face = 1, size(periodic_elem_face_ids_x1(1,:))                    ! Check if the ielem owns a periodic face and if the iface is a periodic face

              if (periodic_elem_face_ids_x1(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x1(2,i_p_face) == iface) then

                match_found = .true.                                                 ! There is a match

                p_dir = periodic_elem_face_ids_x1(3,i_p_face)                        ! Get the direction of periodicity

                do inode = 1, n_LGL_2d                                               ! Loop over the nodes on the On-Element face
                  
                  knode = knode + 1                                                  ! Update the facial node index counter
                  
                  x1 = xg(:,ifacenodes(knode),ielem)                                 ! Save the coordinates of the facial node
                
                  x1_p(:) = Extract_Parallel_Invariant(p_dir,x1)                     !   Extact the two invariant coordinate locations between parallel faces

                  do jnode = 1, n_LGL_2d                                             ! Loop over the nodes on the Off-Element face
                    
                    x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)),ef2e(2,iface,ielem))
                
                    x2_p(:) = Extract_Parallel_Invariant(p_dir,x2)                   !   Extact the two invariant coordinate locations between parallel faces

                    if (magnitude(x1_p-x2_p) <= nodetol) then                        ! Check distance between the two nodes
                      
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem)) ! Set the volumetric node index of the connected node; ef2e(2) gives the element of the neighbor

                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                   ! Set the element of the connected node

                      efn2efn(4,knode,ielem) = jnode                                 ! Set the node index in the off element matching node

                      exit                                                           ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do                                                             ! End do jnode
                
                end do                                                               ! End do inode 

              end if                                                                 ! End if match found

              if (match_found .eqv. .true.) exit                                     ! If a partner face has been found exit from the loop over the elements that own a periodic face

            end do                                                                   ! End do loop over the elements that own a periodic face

          end if                                                                     ! End if check periodic face in x1 direction

          ! If the iface is not a periodic face  in the x1 direction, check if it is a periodic face in the x2 direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x2(1,:)) /= 0) then
           
            do i_p_face = 1, size(periodic_elem_face_ids_x2(1,:))                    ! Check if the ielem owns a periodic face and if the iface is a periodic face

              if (periodic_elem_face_ids_x2(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x2(2,i_p_face) == iface) then

                match_found = .true.                                                 ! There is a match

                p_dir = periodic_elem_face_ids_x2(3,i_p_face)                        ! Get the direction of periodicity

                do inode = 1, n_LGL_2d                                               ! Loop over the nodes on the On-Element face
                  
                  knode = knode + 1                                                  ! Update the facial node index counter
                  
                  x1 = xg(:,ifacenodes(knode),ielem)                                 ! Save the coordinates of the facial node
                
                  x1_p(:) = Extract_Parallel_Invariant(p_dir,x1)                     !   Extact the two invariant coordinate locations between parallel faces

                  do jnode = 1, n_LGL_2d                                             ! Loop over the nodes on the Off-Element face
                    
                    x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)),ef2e(2,iface,ielem))
                
                    x2_p(:) = Extract_Parallel_Invariant(p_dir,x2)                   !   Extact the two invariant coordinate locations between parallel faces

                    if (magnitude(x1_p-x2_p) <= nodetol) then                        ! Check distance between the two nodes
                      
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem)) ! Set the volumetric node index of the connected node; ef2e(2) gives the element of the neighbor

                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                   ! Set the element of the connected node

                      efn2efn(4,knode,ielem) = jnode                                 ! Set the node index in the off element matching node

                      exit                                                           ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do                                                             ! End do jnode
                
                end do                                                               ! End do inode 

              end if                                                                 ! End if match found

              if (match_found .eqv. .true.) exit                                     ! If a partner face has been found exit from the loop over the elements that own a periodic face

            end do                                                                   ! End do loop over the elements that own a periodic face

          end if                                                                     ! End if check periodic face in x2 direction

          ! If the iface is not a periodic face in the x2 direction, check if it is a periodic face in the x3 direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x3(1,:)) /= 0) then
           
            do i_p_face = 1, size(periodic_elem_face_ids_x3(1,:))                    ! Check if the ielem owns a periodic face and if the iface is a periodic face

              if (periodic_elem_face_ids_x3(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x3(2,i_p_face) == iface) then

                match_found = .true.                                                 ! There is a match

                p_dir = periodic_elem_face_ids_x3(3,i_p_face)                        ! Get the direction of periodicity

                do inode = 1, n_LGL_2d                                               ! Loop over the nodes on the On-Element face
                  
                  knode = knode + 1                                                  ! Update the facial node index counter
                  
                  x1 = xg(:,ifacenodes(knode),ielem)                                 ! Save the coordinates of the facial node
                
                  x1_p(:) = Extract_Parallel_Invariant(p_dir,x1)                     !   Extact the two invariant coordinate locations between parallel faces

                  do jnode = 1, n_LGL_2d                                             ! Loop over the nodes on the Off-Element face
                    
                    x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)),ef2e(2,iface,ielem))
                
                    x2_p(:) = Extract_Parallel_Invariant(p_dir,x2)                   !   Extact the two invariant coordinate locations between parallel faces

                    if (magnitude(x1_p-x2_p) <= nodetol) then                        ! Check distance between the two nodes
                      
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem)) ! Set the volumetric node index of the connected node; ef2e(2) gives the element of the neighbor

                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                   ! Set the element of the connected node

                      efn2efn(4,knode,ielem) = jnode                                 ! Set the node index in the off element matching node

                      exit                                                           ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do                                                             ! End do jnode
                
                end do                                                               ! End do inode 

              end if                                                                 ! End if match found

              if (match_found .eqv. .true.) exit                                     ! If a partner face has been found exit from the loop over the elements that own a periodic face

            end do                                                                   ! End do loop over the elements that own a periodic face

          end if                                                                     ! End if check periodic face in x3 direction

          !  ==================================================================
          !  Finished Serial periodic logic:  Begin Serial non-periodic logic
          !  ==================================================================

          if (match_found .eqv. .false.) then

            do inode = 1, n_LGL_2d                                                   ! Loop over the nodes on the face
 
              knode = knode + 1                                                      ! Update the facial node index counter

              x1 = xg(:,ifacenodes(knode),ielem)                                     ! Save coordinates of the facial ndoes

              do jnode = 1, n_LGL_2d                                                 ! Search for connected node on connected element face
                
                                                                                     ! ef2e(1) gives the face on the neighboring element and
                x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)),ef2e(2,iface,ielem)) ! Coordinates of the jnode
               
                if (magnitude(x1-x2) <= nodetol) then                                ! Check the distance between the two nodes

                  efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))     ! Set the volumetric node index of the connected node
                  
                  efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)                       ! Set the element of the connected node
                  
                  efn2efn(4,knode,ielem) = jnode                                     ! Set the index of the connected node
                 
                  exit
                
                end if
              
              end do                                                    ! End do jnode

              if (efn2efn(1,knode,ielem) < 0 .or. efn2efn(2,knode,ielem) < 0) then   ! Print information at screen if there is a problem and stop computation
                write(*,*) 'conforming_interface', conforming_interface
                write(*,*) 'Connectivity error in face-node connectivity_LGL Serial.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface
                write(*,*) 'ef2e'
                write(*,*) ef2e(:,iface,ielem)
                write(*,*)"ef2e(9,iface,ielem) = ",ef2e(9,iface,ielem)
                write(*,*) 'Node coordinates'
                write(*,*) x1
write(*,*)'knode = ',knode
write(*,*)'ifacenodes = ',ifacenodes
write(*,*)"orientation = ",ef2e(7,iface,ielem)

                write(*,*) 'Possible partner node coordinates'
                
                do jnode = 1, n_LGL_2d
                  x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)), ef2e(2,iface,ielem))
                  write(*,*) x2
                end do 

                write(*,*) 'Exiting...'

              end if

            end do ! End do inode
!           if(ef2e(7,iface,ielem) /= 0) write(*,*)ef2e(7,iface,ielem),efn2efn(4,knode+1-n_LGL_2d:knode,ielem)
          
          end if ! End if not a periodic face (match_found = .false.)
              
        else if (ef2e(4,iface,ielem) /= elem_props(2,ielem) .and. (ef2e(3,iface,ielem) /= myprocid)) then
!        else if (((ef2e(4,iface,ielem) /= elem_props(2,ielem)) .and. (ef2e(3,iface,ielem) /= myprocid)).or.&
!                 ((ef2e(9,iface,ielem) == 1) .and. (ef2e(3,iface,ielem) /= myprocid))) then!HACK

          kelem = ef2e(2,iface,ielem)
          call element_properties(kelem, n_pts_2d=n_LGL_2d_Off)
          i_low = i_low + n_LGL_2d_Off
          cycle

        end if ! End if type of face (boundary, off processor or on processor)
      
      end do ! End do loop over faces of the element

    end do ! End do loop elements owned by the processor

  end subroutine calculate_face_node_connectivity_LGL

pure function face_map(n_on,orientation,n_pts_1d)
  integer, intent(in) :: n_on, orientation, n_pts_1d
  
  !--local variables
  integer :: i_on,j_on, i_off, j_off
  integer :: a, b, c, d, e, f
  integer :: face_map
  
  !-- convert from n_on to (i,j) coordinates on on face
  i_on = mod(n_on-1,n_pts_1d) + 1 ; j_on = (n_on-i_on) / n_pts_1d + 1 ;
  
  if(orientation.EQ.0)then
    a = 1; b = 0; c = 0 
    d = 0; e = 1; f = 0
  elseif(orientation.EQ.1)then
    a = 0; b = 1; c = 0
    d = -1; e = 0; f = n_pts_1d+1
  elseif(orientation.EQ.2)then
    a = -1; b = 0; c = n_pts_1d+1
    d = 0; e = -1; f = n_pts_1d+1
  elseif(orientation.EQ.3)then
    a = 0; b = -1; c = n_pts_1d+1
    d = 1; e = 0; f = 0
  elseif(orientation.EQ.4)then
    a = 0; b = 1; c = 0 
    d = 1; e = 0; f = 0
  elseif(orientation.EQ.5)then
    a = -1; b = 0; c = n_pts_1d+1 
    d = 0; e = 1; f = 0
  elseif(orientation.EQ.6)then
    a = 0; b = -1; c = n_pts_1d+1 
    d = -1; e = 0; f = n_pts_1d+1
  elseif(orientation.EQ.7)then
    a = 1; b = 0; c = 0 
    d = 0; e = -1; f = n_pts_1d+1
  endif

  !-- local (i,j) on off face
  i_off = a*i_on + b*j_on + c
  j_off = d*i_on + e*j_on + f

  !-- convert to local numbering
  face_map = (j_off-1)*n_pts_1d+i_off
end function 
  !============================================================================
  
  pure function Extract_Parallel_Invariant(p_dir,x)

    ! Extract from x1 the two invaraint coordinates
    integer,                intent(in) :: p_dir
    real(wp), dimension(3), intent(in) :: x

    integer :: cnt_coord, i_coord

    real(wp), dimension(2)             :: Extract_Parallel_Invariant

    continue

        cnt_coord = 0

        do i_coord = 1,3

          if (i_coord /= p_dir) then
            cnt_coord = cnt_coord + 1
            Extract_Parallel_Invariant(cnt_coord) = x(i_coord)
          end if

        end do

   end function Extract_Parallel_Invariant

  !============================================================================
  
  subroutine calcfacenormals_LGL(toggle)
    ! this subroutine calculates the outward facing normals
    ! of each facial node
    use referencevariables
    use variables, only: kfacenodes, facenodenormal, r_x, ef2e, efn2efn, Jx_r, &
                       Jx_facenodenormal_LGL 
    use initcollocation,      only: element_properties
    use collocationvariables, only: elem_props

    implicit none

    logical, intent(in)   :: toggle
    ! indices
    integer :: ielem, kelem, inode, jnode, knode, lnode, iface, kface, idir
    integer :: i, face_shift, n_pts_1d_max
    integer :: n_LGL_2d, nodesperface_max

    real(wp) :: dx
    real(wp), dimension(3) :: w1, w2
   !real(wp), dimension(3) :: xg_target=(/1.5_wp,1.0_wp,0.0_wp/)
    logical                :: testing = .false., nonconforming_element

    ! number of nodes in each element

    nodesperface_max = (npoly_max+1)**(ndim-1)
    if(allocated(facenodenormal)) deallocate(facenodenormal)
    allocate(   facenodenormal    (3,nfacesperelem*nodesperface_max,ihelems(1):ihelems(2)))

    facenodenormal     = -1000000.0_wp
    if(toggle)then    
    else
      if(allocated(Jx_facenodenormal_LGL)) deallocate(Jx_facenodenormal_LGL)
      allocate(Jx_facenodenormal_LGL(3,nfacesperelem*nodesperface_max,ihelems(1):ihelems(2)))
      Jx_facenodenormal_LGL = -5000000.0_wp
    endif

    n_pts_1d_max = (npoly_max+1)**1
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)

       call element_properties(ielem,             &
                              n_pts_2d=n_LGL_2d,  &
                            kfacenodes=kfacenodes )

      knode = 0                                  ! reset facial node index counter
                                                 ! compute outward facing normals
      do iface = 1,nfacesperelem                 ! loop over faces
        face_shift = (iface-1)*n_pts_1d_max**2
        do inode = 1,n_LGL_2d                    ! loop over nodes on face

          knode = knode + 1                      ! update facial node index counter

          i = kfacenodes(inode,iface)            ! volumetric node index of facial node

          idir = abs(elfacedirections(iface))    ! unsigned direction of face in computational space

          dx = sign(1.0_wp,real(elfacedirections(iface),wp)) ! sign so normal is facing outward

                                                 ! outward facing normal divided by Jacobian 
             facenodenormal    (:,knode,ielem) =               r_x(idir,:,i,ielem) * dx

                                                 ! outward facing normal using metrics
          if(toggle)then
            !-- this checks if you are on the second computation of Jx_facenormal_LGL (toggle = .true.) and if you are on a boundary face
            !   in such a case you do not want to overwrite the metric terms
          else
            Jx_facenodenormal_LGL(1:3,knode,ielem) =     Jx_r(i,ielem)*r_x(idir,:,i,ielem) * dx
          endif

        end do

      end do

    end do


    ! testing facenodenormal calculations

    if(testing) then
      elem_loop:do ielem = ihelems(1), ihelems(2)

        call element_properties(ielem,             &
                               n_pts_2d=n_LGL_2d,  &
                             kfacenodes=kfacenodes )
        nonconforming_element = .false.
       !-- determine if this is an element with any face that is nonconforming
          do iface = 1,nfacesperelem
            if (ef2e(4,iface,ielem) /= elem_props(2,ielem)) then
              nonconforming_element = .true.
              exit
            endif
          enddo

        if(nonconforming_element)then
          !-- this is an element with a nonconforming face do nothing
        else
          face_loop:do iface = 1,nfacesperelem                      ! loop over faces

            !-- determine if the adjoining element is a fully conforming element 
            kelem = ef2e(2,iface,ielem)                   ! adjoining element

            nonconforming_element = .false.

            do kface = 1,nfacesperelem
              if (ef2e(4,kface,kelem) /= elem_props(2,kelem)) then
                nonconforming_element = .true.
                exit
              endif
            enddo

  !          if (ef2e(4,iface,ielem) == elem_props(2,ielem)) then
             if(nonconforming_element)then                    !-- nonconforming
               !-- adjoining elmenet has at least one nonconforming face do nothing
             else

               if (ef2e(3,iface,ielem) /= myprocid) then       !  Parallel conforming path

               else                                            !  Serial conforming path

!                kelem = ef2e(2,iface,ielem)                   ! adjoining element
                  kface = ef2e(1,iface,ielem)                   ! face on element

                  do i = 1,n_LGL_2d                             ! loop over nodes on face

                    if(ef2e(1,iface,ielem) > 0)then             !-- not a boundary

                      jnode =  n_LGL_2d*(iface-1) + i             ! Index in facial ordering
            
                      inode = kfacenodes(i,iface)                 ! Volumetric node index corresponding to facial node index

                      knode = efn2efn(1,jnode,ielem)              ! Volumetric node index of adjoining face point 

                      lnode = (kface-1)*n_LGL_2d + efn2efn(4,jnode,ielem)

                      w1(:) = facenodenormal(1:3,jnode,ielem)*Jx_r(inode,ielem)
                      w2(:) = facenodenormal(1:3,lnode,kelem)*Jx_r(knode,kelem)
 
                      if(magnitude(w1+w2) >= 1.0e-10_wp) then
                        write(*,*)'facenodenormals are incorrect'
                        write(*,*)'ielement ',ielem,' iface = ',iface,' face normal = ',w1
                        write(*,*)'kelement ',kelem,' kface = ',kface,' face normal = ',w2
                        write(*,*)'nonconforming_element = ',nonconforming_element
                      endif

                    endif!-- not a boundary
                  end do
               endif!--Parallel conforming path
             endif!--nonconforming
          end do face_loop
        endif
      end do elem_loop
    endif

  end subroutine calcfacenormals_LGL

  !============================================================================
  
  subroutine calcfacenormals_Gau()
    ! this subroutine calculates the outward facing normals
    ! of each facial node
    use referencevariables
    use initcollocation,      only: element_properties, ExtrpXA2XB_2D_neq, Gauss_Legendre_points
    use collocationvariables, only: Restrct_Gau_2_LGL_1d, elem_props
    use variables, only: Jx_facenodenormal_Gau, Jx_facenodenormal_LGL, xg_Gau_Shell, ef2e

    implicit none

    ! indices
    integer :: ielem, iface, knode, node_id
    integer :: n_pts_1d_max, n_pts_2d_max
    integer :: n_S_1d_On, n_S_1d_Off, n_S_1d_Mort, poly_val, istart_Mort, istart_On,& 
               iend_Mort, iend_On, n_s_2d_On

    logical                               :: testing = .false.
    real(wp), dimension(:),   allocatable :: x_S_1d_Mort, x_S_1d_On, w_S_1d_Mort
    real(wp), dimension(:,:), allocatable :: Intrp

    n_pts_1d_max = (npoly_max+1)**1
    n_pts_2d_max = (npoly_max+1)**2

    allocate(Jx_facenodenormal_Gau(3,nfacesperelem*n_pts_2d_max,ihelems(1):ihelems(2)))
    Jx_facenodenormal_Gau = 0.0_wp

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,              &
                              n_pts_1d=n_S_1d_On, &
                              x_pts_1d=x_S_1d_On, &
                              n_pts_2d = n_s_2d_On)

      knode = 0                                  !  reset facial node index counter
      node_id = 0

      do iface = 1,nfacesperelem                 ! loop over faces
        if(elem_props(2,ielem) == ef2e(4,iface,ielem)) then
          !-- conforming: do nothing
        elseif (ef2e(1,iface,ielem) < 0)then
          !-- boundary: do nothing
        elseif(elem_props(2,ielem) /= ef2e(4,iface,ielem))then
          !-- nonconforming interface

        n_S_1d_Off  = ef2e(4,iface,ielem)
        n_S_1d_Mort = max(n_S_1d_On,n_S_1d_Off)

        if(allocated(x_S_1d_Mort)) deallocate(x_S_1d_Mort) ; allocate(x_S_1d_Mort(n_S_1d_Mort)) ;
        if(allocated(w_S_1d_Mort)) deallocate(w_S_1d_Mort) ; allocate(w_S_1d_Mort(n_S_1d_Mort)) ;
        call Gauss_Legendre_points(n_S_1d_Mort,x_S_1d_Mort,w_S_1d_Mort)
      
       ! compute OUTWARD FACING normals on mortars
        call Shell_Metrics_Analytic(iface,n_pts_1d_max,n_S_1d_Mort,x_S_1d_Mort,     &
                              xg_Gau_shell(:,:,ielem),Jx_facenodenormal_Gau(:,:,ielem))

        ! Grab the correct restriction operator between mortar and on-element face
        if(allocated(Intrp)) deallocate(Intrp) ;
        if(n_S_1d_Mort == n_S_1d_On) then
          poly_val = n_S_1d_Mort - npoly
          allocate(Intrp(n_S_1d_On,n_S_1d_Mort)) ;  Intrp(:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_On,1:n_S_1d_Mort,poly_val,1) ;
        else
          poly_val = n_S_1d_On - npoly
          allocate(Intrp(n_S_1d_On,n_S_1d_Mort)) ;  Intrp(:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_On,1:n_S_1d_Mort,poly_val,2) ;
        endif

        ! Restrict all metric data from mortar back to element face
        istart_Mort = (iface-1)*n_pts_2d_max  + 1
        istart_On   = (iface-1)*n_s_2d_On  + 1
!        istart_On = istart_Mort
          iend_Mort = istart_Mort + n_S_1d_Mort**2 - 1
          iend_On   = istart_On + n_S_1d_On**2   - 1
!        iend_On = iend_Mort

        call ExtrpXA2XB_2D_neq(3, n_S_1d_Mort, n_S_1d_On,x_S_1d_Mort,x_S_1d_On, &
           Jx_facenodenormal_Gau(:,istart_Mort:iend_Mort,ielem),Jx_facenodenormal_LGL(:,istart_On:iend_On,ielem),&
                          Intrp, Intrp)
        endif
      end do

    end do

    ! testing facenodenormal calculations

    if(testing) then
!     do ielem = ihelems(1), ihelems(2)
!       knode = 0
!       ! loop over faces
!       do iface = 1,nfacesperelem
!         ! loop over nodes on face
!         kelem = ef2e(2,iface,ielem)
!         do inode = 1,n_Gau_2d_p1
!           knode = knode + 1
!           if(ef2e(1,iface,ielem) > 0)then
!             i = (ef2e(1,iface,ielem)-1)*n_Gau_2d_p1+efn2efn(4,knode,ielem)
!             wrk = facenodenormal(1:3,knode,ielem)*Jx_r(kfacenodes(inode,iface),ielem) &
!                 + facenodenormal(1:3, i ,kelem)*Jx_r(efn2efn(1,knode,ielem),kelem)
!             if(magnitude(wrk) >= 1.0e-10_wp) then
!               write(*,*)'facenodenormals are incorrect'
!               write(*,*)facenodenormal(1:2,knode,ielem)*Jx_r(kfacenodes(inode,iface),ielem) &
!                       , facenodenormal(1:2, i ,kelem)*Jx_r(efn2efn(1,knode,ielem),kelem)
!             endif
!           endif
!         end do
!       end do
!     end do
    endif

  end subroutine calcfacenormals_Gau

  !============================================================================

  pure function cross_product(a, b)
    
    ! Nothing is implicitly defined
    real(wp), dimension(3) :: cross_product
    real(wp), dimension(3), intent(in) :: a, b

    continue

    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product

  !============================================================================

  subroutine calcmetrics_LGL()
    ! This subroutine calculates the metric transformations
    ! between computational and physical space.
    use referencevariables
    use variables, only: xg, x_r, r_x, Jx_r, dx_min_elem
    use collocationvariables, only: nnzgrad, iagrad, jagrad, dagrad, pvol
    use initcollocation, only: element_properties
    use controlvariables, only: symmetric_metric
    use mpimod

    implicit none

    integer :: ielem, inode, jnode, idir, jdir
    integer :: i, j, l, ii, iE
    integer :: nodesperelem_max, n_LGL_3d
!   integer :: ierr

    real(wp), parameter :: tol = 1.0e-11_wp

    real(wp) :: test(3)
    real(wp) :: err,err_L2,err_Linf
    real(wp), dimension(:), allocatable :: err_max_proc

    integer :: s_tag, r_tag, m_size, m, &
               s_request_err_Linf, r_request_err_Linf, i_err

    integer :: s_status(mpi_status_size)
    integer :: r_status(mpi_status_size)

    continue 

    nodesperelem_max = (npoly_max+1)**ndim                    ! number of nodes in each element

    allocate(x_r(3,3,1:nodesperelem_max,ihelems(1):ihelems(2))) ;  x_r = 0.0_wp     ! dx/dr
    allocate(r_x(3,3,1:nodesperelem_max,ihelems(1):ihelems(2))) ;  r_x = 0.0_wp     ! dr/dx
    allocate(Jx_r(   1:nodesperelem_max,ihelems(1):ihelems(2))) ; Jx_r = 0.0_wp     ! J = |dx/dr|

    allocate(dx_min_elem(ihelems(1):ihelems(2))) ; dx_min_elem(:) = 0.0_wp ;

    err_L2 = 0.0_wp ; err_Linf = 0.0_wp                       ! Initialize metrics error norms

    elloop:do ielem = ihelems(1), ihelems(2)                  ! loop over volumetric elements

      call element_properties(ielem,         &
                          n_pts_3d=n_LGL_3d, &
                           nnzgrad=nnzgrad,  &
                            iagrad=iagrad,   &
                            jagrad=jagrad,   &
                            dagrad=dagrad,   &
                              pvol=pvol)

                                                              ! initialize dx/dr to identity and dr/dx to identity
      do idir = 1,3
        x_r(idir,idir,1:n_LGL_3d,ielem) = 1.0_wp
        r_x(idir,idir,1:n_LGL_3d,ielem) = 1.0_wp
      end do
                                                              ! calculate dx/dr

      x_r(:,1:ndim,1:n_LGL_3d,ielem) = 0.0_wp                 ! initialize to zero
      do inode = 1, n_LGL_3d                                  ! loop over every node in element

        do jdir = 1,ndim                                      ! loop over dimension of gradient

          do i = iagrad(inode), iagrad(inode+1)-1             ! loop over number of dependent nodes in gradient
                                                              ! column/node from gradient operator in CSR format in
            jnode = jagrad(jdir,i)                            ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
                                                              ! update gradient. MP: Well, actually this is the Jacobian of the transformation
            x_r(:,jdir,inode,ielem) = x_r(:,jdir,inode,ielem) + dagrad(jdir,i)*xg(:,jnode,ielem)

          end do

        end do

      end do

      ! calculate determinant
      do inode = 1,n_LGL_3d
        Jx_r(inode,ielem) = determinant3(x_r(:,:,inode,ielem))
      end do

      if (ndim < 3) then
                                                               ! inverse metrics (note that in 3D this is not sufficient to satisfy the GCL.
        do inode = 1,n_LGL_3d
          r_x(1,1,inode,ielem) =  x_r(2,2,inode,ielem)/Jx_r(inode,ielem)
          r_x(2,1,inode,ielem) = -x_r(2,1,inode,ielem)/Jx_r(inode,ielem)
          r_x(1,2,inode,ielem) = -x_r(1,2,inode,ielem)/Jx_r(inode,ielem)
          r_x(2,2,inode,ielem) =  x_r(1,1,inode,ielem)/Jx_r(inode,ielem)
        end do

      else

!        Metric data is stored as follows
!              --                                  --
!              | dxi_1/dx_1, dxi_1/dx_2, dxi_1/dx_3 |
!              |                                    |
!        r_x = | dxi_2/dx_1, dxi_2/dx_2, dxi_2/dx_3 |
!              |                                    |
!              | dxi_3/dx_1, dxi_3/dx_2, dxi_3/dx_3 |
!              --                                  --
!  
!       ifacenodenormal(idir,1:3) = r_x(idir,1:3)    i.e., the ith row of r_x
!

        ! Special formulas are used to ensure that GCLs are satisfied in each direction

        iE = ielem                                 !  New variable strictly cosmetic so formulae are shorter
        r_x(:,:,:,ielem) = 0.0_wp

        if(symmetric_metric) then        !   Symmetric form (.true) is taken from Sjogreen.Yee.Vinokur.LLNL_TR_637397.HOFD.Metrics.GCL.Moving.Meshes.pdf

          do inode = 1, n_LGL_3d

            jdir = 2
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(1,1,inode,iE) = r_x(1,1,inode,iE) + dagrad(jdir,i)*( + x_r(3,3,j,iE)*xg(2,j,iE) - x_r(2,3,j,iE)*xg(3,j,iE)) ! dxi_1/dx_1
              r_x(1,2,inode,iE) = r_x(1,2,inode,iE) + dagrad(jdir,i)*( + x_r(1,3,j,iE)*xg(3,j,iE) - x_r(3,3,j,iE)*xg(1,j,iE)) ! dxi_1/dx_2
              r_x(1,3,inode,iE) = r_x(1,3,inode,iE) + dagrad(jdir,i)*( + x_r(2,3,j,iE)*xg(1,j,iE) - x_r(1,3,j,iE)*xg(2,j,iE)) ! dxi_1/dx_3
            end do
            jdir = 3
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(1,1,inode,iE) = r_x(1,1,inode,iE) + dagrad(jdir,i)*( + x_r(2,2,j,iE)*xg(3,j,iE) - x_r(3,2,j,iE)*xg(2,j,iE)) ! dxi_1/dx_1
              r_x(1,2,inode,iE) = r_x(1,2,inode,iE) + dagrad(jdir,i)*( + x_r(3,2,j,iE)*xg(1,j,iE) - x_r(1,2,j,iE)*xg(3,j,iE)) ! dxi_1/dx_2
              r_x(1,3,inode,iE) = r_x(1,3,inode,iE) + dagrad(jdir,i)*( + x_r(1,2,j,iE)*xg(2,j,iE) - x_r(2,2,j,iE)*xg(1,j,iE)) ! dxi_1/dx_3
            end do
  
            jdir = 3
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(2,1,inode,iE) = r_x(2,1,inode,iE) + dagrad(jdir,i)*( + x_r(3,1,j,iE)*xg(2,j,iE) - x_r(2,1,j,iE)*xg(3,j,iE)) ! dxi_2/dx_1
              r_x(2,2,inode,iE) = r_x(2,2,inode,iE) + dagrad(jdir,i)*( + x_r(1,1,j,iE)*xg(3,j,iE) - x_r(3,1,j,iE)*xg(1,j,iE)) ! dxi_2/dx_2
              r_x(2,3,inode,iE) = r_x(2,3,inode,iE) + dagrad(jdir,i)*( + x_r(2,1,j,iE)*xg(1,j,iE) - x_r(1,1,j,iE)*xg(2,j,iE)) ! dxi_2/dx_3
            end do
            jdir = 1
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(2,1,inode,iE) = r_x(2,1,inode,iE) + dagrad(jdir,i)*( + x_r(2,3,j,iE)*xg(3,j,iE) - x_r(3,3,j,iE)*xg(2,j,iE)) ! dxi_2/dx_1
              r_x(2,2,inode,iE) = r_x(2,2,inode,iE) + dagrad(jdir,i)*( + x_r(3,3,j,iE)*xg(1,j,iE) - x_r(1,3,j,iE)*xg(3,j,iE)) ! dxi_2/dx_2
              r_x(2,3,inode,iE) = r_x(2,3,inode,iE) + dagrad(jdir,i)*( + x_r(1,3,j,iE)*xg(2,j,iE) - x_r(2,3,j,iE)*xg(1,j,iE)) ! dxi_2/dx_3
            end do
  
            jdir = 1
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(3,1,inode,iE) = r_x(3,1,inode,iE) + dagrad(jdir,i)*( + x_r(3,2,j,iE)*xg(2,j,iE) - x_r(2,2,j,iE)*xg(3,j,iE)) ! dxi_3/dx_1
              r_x(3,2,inode,iE) = r_x(3,2,inode,iE) + dagrad(jdir,i)*( + x_r(1,2,j,iE)*xg(3,j,iE) - x_r(3,2,j,iE)*xg(1,j,iE)) ! dxi_3/dx_2
              r_x(3,3,inode,iE) = r_x(3,3,inode,iE) + dagrad(jdir,i)*( + x_r(2,2,j,iE)*xg(1,j,iE) - x_r(1,2,j,iE)*xg(2,j,iE)) ! dxi_3/dx_3
            end do
            jdir = 2
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(3,1,inode,iE) = r_x(3,1,inode,iE) + dagrad(jdir,i)*( + x_r(2,1,j,iE)*xg(3,j,iE) - x_r(3,1,j,iE)*xg(2,j,iE)) ! dxi_3/dx_1
              r_x(3,2,inode,iE) = r_x(3,2,inode,iE) + dagrad(jdir,i)*( + x_r(3,1,j,iE)*xg(1,j,iE) - x_r(1,1,j,iE)*xg(3,j,iE)) ! dxi_3/dx_2
              r_x(3,3,inode,iE) = r_x(3,3,inode,iE) + dagrad(jdir,i)*( + x_r(1,1,j,iE)*xg(2,j,iE) - x_r(2,1,j,iE)*xg(1,j,iE)) ! dxi_3/dx_3
            end do

            r_x(:,:,inode,iE) = 0.5_wp * r_x(:,:,inode,iE)/Jx_r(inode,iE)

          end do

        else

          do inode = 1, n_LGL_3d

            jdir = 2
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(1,1,inode,iE) = r_x(1,1,inode,iE) + dagrad(jdir,i)*(                            - x_r(2,3,j,iE)*xg(3,j,iE)) ! dxi_1/dx_1
              r_x(1,2,inode,iE) = r_x(1,2,inode,iE) + dagrad(jdir,i)*(                            - x_r(3,3,j,iE)*xg(1,j,iE)) ! dxi_1/dx_2
              r_x(1,3,inode,iE) = r_x(1,3,inode,iE) + dagrad(jdir,i)*(                            - x_r(1,3,j,iE)*xg(2,j,iE)) ! dxi_1/dx_3
            end do
            jdir = 3
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(1,1,inode,iE) = r_x(1,1,inode,iE) + dagrad(jdir,i)*( + x_r(2,2,j,iE)*xg(3,j,iE)                           ) ! dxi_1/dx_1
              r_x(1,2,inode,iE) = r_x(1,2,inode,iE) + dagrad(jdir,i)*( + x_r(3,2,j,iE)*xg(1,j,iE)                           ) ! dxi_1/dx_2
              r_x(1,3,inode,iE) = r_x(1,3,inode,iE) + dagrad(jdir,i)*( + x_r(1,2,j,iE)*xg(2,j,iE)                           ) ! dxi_1/dx_3
            end do
  
            jdir = 3
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(2,1,inode,iE) = r_x(2,1,inode,iE) + dagrad(jdir,i)*(                            - x_r(2,1,j,iE)*xg(3,j,iE)) ! dxi_2/dx_1
              r_x(2,2,inode,iE) = r_x(2,2,inode,iE) + dagrad(jdir,i)*(                            - x_r(3,1,j,iE)*xg(1,j,iE)) ! dxi_2/dx_2
              r_x(2,3,inode,iE) = r_x(2,3,inode,iE) + dagrad(jdir,i)*(                            - x_r(1,1,j,iE)*xg(2,j,iE)) ! dxi_2/dx_3
            end do
            jdir = 1
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(2,1,inode,iE) = r_x(2,1,inode,iE) + dagrad(jdir,i)*( + x_r(2,3,j,iE)*xg(3,j,iE)                           ) ! dxi_2/dx_1
              r_x(2,2,inode,iE) = r_x(2,2,inode,iE) + dagrad(jdir,i)*( + x_r(3,3,j,iE)*xg(1,j,iE)                           ) ! dxi_2/dx_2
              r_x(2,3,inode,iE) = r_x(2,3,inode,iE) + dagrad(jdir,i)*( + x_r(1,3,j,iE)*xg(2,j,iE)                           ) ! dxi_2/dx_3
            end do
  
            jdir = 1
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(3,1,inode,iE) = r_x(3,1,inode,iE) + dagrad(jdir,i)*(                            - x_r(2,2,j,iE)*xg(3,j,iE)) ! dxi_3/dx_1
              r_x(3,2,inode,iE) = r_x(3,2,inode,iE) + dagrad(jdir,i)*(                            - x_r(3,2,j,iE)*xg(1,j,iE)) ! dxi_3/dx_2
              r_x(3,3,inode,iE) = r_x(3,3,inode,iE) + dagrad(jdir,i)*(                            - x_r(1,2,j,iE)*xg(2,j,iE)) ! dxi_3/dx_3
            end do
            jdir = 2
            do i = iagrad(inode), iagrad(inode+1)-1
              j = jagrad(jdir,i)
              r_x(3,1,inode,iE) = r_x(3,1,inode,iE) + dagrad(jdir,i)*( + x_r(2,1,j,iE)*xg(3,j,iE)                           ) ! dxi_3/dx_1
              r_x(3,2,inode,iE) = r_x(3,2,inode,iE) + dagrad(jdir,i)*( + x_r(3,1,j,iE)*xg(1,j,iE)                           ) ! dxi_3/dx_2
              r_x(3,3,inode,iE) = r_x(3,3,inode,iE) + dagrad(jdir,i)*( + x_r(1,1,j,iE)*xg(2,j,iE)                           ) ! dxi_3/dx_3
            end do

            r_x(:,:,inode,iE) = r_x(:,:,inode,iE)/Jx_r(inode,iE)

          end do

        endif   !  symmetric_metric portion of loop

      end if

!     Metric Test
! 
!     \frac{\partial}{\partial \xi^j} \left( J \frac{\partial \xi^j}{\partial x^i} \right) == 0 ; i,j = 1,3
!     ( Assumes Einstein Convention on ``j'' )
!
      do inode = 1,n_LGL_3d
        do i = 1, ndim
          test(i) = 0.0_wp
          do l = 1, ndim
            do ii = iagrad(inode), iagrad(inode+1)-1
              jnode = jagrad(l,ii)
              test(i) = test(i) + dagrad(l,ii) * r_x(l,i,jnode,ielem)*Jx_r(jnode,ielem)
            end do
          end do
        end do
        err = maxval(abs(test(:)))
        err_L2   = err_L2 + err*err
        err_Linf = max(err_Linf,err)
      end do

      !  Calculate a conservative estimate of the characteristic length of an element.  
      !  Needs further development

      dx_min_elem(ielem) = 0.0_wp
      do inode = 1,n_LGL_3d
        dx_min_elem(ielem) = dx_min_elem(ielem) + pvol(inode)*Jx_r(inode,ielem)
      end do

      dx_min_elem(ielem) = dx_min_elem(ielem)**(0.333333333333333333333333_wp)

      if(minval(Jx_r(:,ielem)) <= tol) then
!       write(*,*)' stopping !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!       write(*,*)'jacobian in element = ',ielem, 'on process = ',myprocid, 'is negative'
!       write(*,*)' stopping !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!       call PetscFinalize(ierr) ; stop
      endif
    end do elloop

!   Parallel reduction of errors using conventional MPI calls

    ! Reduce values on all processes to a single value
    if(myprocid == 0 )  then
      allocate(err_max_proc(0:nprocs-1))
      err_max_proc(:) = 0.0_wp ; err_max_proc(0) = err_Linf ;
    endif
    if(myprocid /= 0 ) then
      s_tag = 100 + myprocid
      m_size = 1
      call mpi_isend(err_Linf,m_size,mpi_double,0,s_tag,petsc_comm_world, &
        & s_request_err_Linf,i_err)
      
      call mpi_wait(s_request_err_Linf,s_status,i_err)
    else
      do m = 1, nprocs-1
        r_tag = 100 + m
        m_size = 1
        call mpi_irecv(err_max_proc(m),m_size,mpi_double,m,r_tag, &
          & petsc_comm_world,r_request_err_Linf,i_err)

        call mpi_wait(r_request_err_Linf,r_status,i_err)
      enddo
    endif

    ! Write at screen the L_inf of the metric error
    if(myprocid == 0 )  then
      write(*,*)'  Metric error: L_inf   ',maxval(err_max_proc(:))
      write(*,*) '==========================================================='

      deallocate(err_max_proc)
    endif

  end subroutine calcmetrics_LGL

  !============================================================================

  subroutine modify_metrics_nonconforming()

    use controlvariables,     only: verbose
    use referencevariables
    use collocationvariables, only: elem_props
    use initcollocation,      only: GCL_Triple_Qmat_Transpose, element_properties, GCL_Triple_HinvDQSEmat
    use eispack_module,       only: svd
    use unary_mod,            only: qsortd
    use variables,            only: r_x, Jx_r, ef2e
    use mpimod,               only: mpi_integer, mpi_double, mpi_status_size, mpi_sum, petsc_comm_world

    implicit none

    logical, parameter :: matu = .True.
    logical, parameter :: matv = .True.

    integer,  dimension(:),    allocatable :: ifacenodes
    real(wp), dimension(:,:),  allocatable :: qmat, dmat
    real(wp), dimension(:),    allocatable :: p_surf

    real(wp), dimension(:,:),  allocatable :: Amat, Hinv_3_mat, D_3_mat, Q_3_mat, S_3_mat, E_3_mat
    real(wp), dimension(:,:),  allocatable :: bvec
    real(wp), dimension(:,:),  allocatable :: a_t
    real(wp), dimension(:,:),  allocatable :: u, v
    real(wp), dimension(:),    allocatable :: w, work, pmat
    real(wp), dimension(:,:),  allocatable :: verror, verror2, verror3
    real(wp), dimension(:,:),  allocatable :: work3
    real(wp), dimension(:,:),  allocatable :: eye, wrk, diag, wI
    real(wp), dimension(:,:),  allocatable :: delta_a

    integer,  dimension(:),    allocatable :: perm
    integer :: ielem, inode, ierr
    integer :: i, iface
    integer :: n_pts_1d, n_pts_2d, n_pts_3d
    integer :: nm, m, n, icnt

    logical                                :: testing = .true., testing_metric_comp = .true.
    logical                                :: modify_metrics = .false.
    real(wp)                               :: t1, t2
    real(wp), parameter                    :: tol = 1.0e-12_wp

    real(wp) :: err,err_L2,err_Linf
    real(wp), dimension(:), allocatable :: err_max_proc

    integer :: s_tag, r_tag, m_size, s_request_err_Linf, r_request_err_Linf
    integer :: np_mods, ng_mods

    integer :: s_status(mpi_status_size)
    integer :: r_status(mpi_status_size)
!-- unccoment to write matrices to file
!   character(len=1024)                            :: numb

    np_mods = 0 ; ng_mods = 0 ;
    err_Linf = 0.0_wp

    ! loop over volumetric elements
    elloop:do ielem = ihelems(1), ihelems(2)

   
      modify_metrics = .false.
      icnt = 0
      do iface = 1,nfacesperelem
         if (ef2e(4,iface,ielem) /= elem_props(2,ielem)) then
            modify_metrics = .true.
            exit
         endif
      enddo

      if(modify_metrics .eqv. .false.) cycle  ! don't modify metrics if element is fully conforming

      np_mods = np_mods + 1                   !  Count the number of elements that are modified
      
!      call element_properties(ielem,&
!                              n_pts_1d=n_pts_1d,&
!                              n_pts_2d=n_pts_2d,&
!                              n_pts_3d=n_pts_3d,&
!                                  qmat=qmat,&
!                                  pmat=pmat,&
!                                p_surf=p_surf,&
!                            ifacenodes=ifacenodes)


      call element_properties(ielem,&
                              n_pts_1d=n_pts_1d,&
                              n_pts_2d=n_pts_2d,&
                              n_pts_3d=n_pts_3d,&
                                  qmat=qmat,&
                                  dmat=dmat,&
                                  pmat=pmat,&
                                p_surf=p_surf,&
                            ifacenodes=ifacenodes)

                     
      m = 1*n_pts_3d ; n = 3*n_pts_3d ; nm = max(m,n) ;

      if(allocated(Amat)) deallocate(Amat) ; allocate(Amat(nm,n)) ; Amat(:,:) = 0.0_wp
      if(allocated(   v)) deallocate(   v) ; allocate(   v(nm,n)) ;    v(:,:) = 0.0_wp
      if(allocated(   u)) deallocate(   u) ; allocate(   u(nm,n)) ;    u(:,:) = 0.0_wp
      if(allocated(   w)) deallocate(   w) ; allocate(   w(nm))   ;    w(:)   = 0.0_wp
      if(allocated(work)) deallocate(work) ; allocate(work(nm))   ; work(:)   = 0.0_wp
      if(allocated(  wI)) deallocate(  wI) ; allocate(  wI(m,m))  ;   wI(:,:) = 0.0_wp

      if(allocated(bvec)) deallocate(bvec) ; allocate( bvec(n_pts_3d,3)) ; bvec(:,:) = 0.0_wp
      if(allocated( a_t)) deallocate( a_t) ; allocate(a_t(3*n_pts_3d,3)) ;  a_t(:,:) = 0.0_wp
      if(allocated(delta_a)) deallocate(delta_a) ; allocate(delta_a(3*n_pts_3d,3)) ;  delta_a(:,:) = 0.0_wp

      call GCL_Triple_Qmat_Transpose(n_pts_1d, n_pts_3d, pmat, qmat, Amat)     

      call SVD(nm,n,m,transpose(Amat),w,matu,u,matv,v,ierr,work)

      !-- construct the pseudo inverse of w
      do i = 1,m
        if(w(i) >= tol) wI(i,i) = 1.0_wp / w(i)
      enddo

      if(testing) then
        if(allocated( eye)) deallocate( eye) ; allocate( eye(m,m))  ;  eye(:,:) = 0.0_wp
        if(allocated( wrk)) deallocate( wrk) ; allocate( wrk(m,m))  ;  wrk(:,:) = 0.0_wp
        if(allocated(diag)) deallocate(diag) ; allocate(diag(m,m))  ; diag(:,:) = 0.0_wp
        if(allocated(perm)) deallocate(perm) ; allocate(perm(m))    ; perm(:)   = 0
        do i = 1,m
          eye(i,i) = 1.0_wp
        enddo
          
        do i = 1,m
          diag(i,i) = w(i)
          !if(w(i) >= tol) wI(i,i) = 1.0_wp / w(i)
        enddo
        call qsortd(w(1:m),perm,m)
        wrk(:,:) = matmul(transpose(v(1:m,1:m)),v(1:m,1:m)) &
                 + matmul(v(1:m,1:m),transpose(v(1:m,1:m))) &
                 + matmul(transpose(u(1:n,1:m)),u(1:n,1:m)) 
        t1 = maxval(abs(3*eye(:,:) - wrk(:,:))) 
        t2 = maxval(abs(transpose(amat(1:m,1:n))-matmul(u(1:n,1:m),matmul(diag(1:m,1:m),transpose(v(1:m,1:m))) )))
        if(w(perm(2)) <= tol) write(*,*)'second singular mode',w(perm(1:2))
        if(t1+t2 > tol) then
          write(*,*)'second singular mode',w(perm(1:2))
          write(*,*)'error in u^T u', t1
          write(*,*)'error in A - u S v^T', t2
          write(*,*)'stopping'
          call PetscFinalize(ierr) ; Stop ;
          stop
         endif
      endif


      call Load_Mortar_Metric_Data(ielem,n_pts_1d,n_pts_2d,n_pts_3d, ifacenodes, p_surf, a_t, bvec)

      !-- test to see if bvec is orthogonal to the constant vector
      !   if it is not then it does not satisfy the assumptions behind the solution process for the metrics
      if( (abs(sum(bvec(:,1))).GE.1.0e-12_wp).OR.(abs(sum(bvec(:,2))).GE.1.0e-12_wp)&
          .OR.(abs(sum(bvec(:,3))).GE.1.0e-12_wp))then
        write(*,*)'======================================'
        write(*,*)'ielem = ',ielem
        write(*,*)'bvec is not orthogonal to the constant vector'
        write(*,*)'sum(bvec(:,1)) = ',sum(bvec(:,1))
        write(*,*)'sum(bvec(:,2)) = ',sum(bvec(:,2))
        write(*,*)'sum(bvec(:,3)) = ',sum(bvec(:,3))
        write(*,*)'======================================'
      endif

!     minimizing (a - a_t)(a - a_t) / 2 ; subject to M a = b
!     "t" subscript denotes "target metric values' :  i.e., a_t
!     a = a_t - M^* (M a_t - b)
      
      if(allocated(work3)) deallocate(work3) ; allocate(work3(1:m,1:3))   ; work3(:,:)   = 0.0_wp

      work3(1:m,1:3) = matmul(Amat(1:m,1:n),a_t(1:n,1:3)) - bvec(1:m,1:3)

      t1 = maxval(abs(work3)) ; 
      if(t1 >= tol .and. verbose ) write(*,*)'A a_t - bvec', t1

      delta_a = matmul( u(1:n,1:m),           &
                  matmul(wI(1:m,1:m),           &
                    matmul(transpose(v(1:m,1:m)), &
                      matmul(Amat(1:m,1:n),a_t(1:n,1:3)) - bvec(1:m,1:3))))

      err = maxval(abs(delta_a)) ; 
      err_L2   = err_L2 + err*err
      err_Linf = max(err_Linf,err)

      a_t(:,:) = a_t(:,:) - delta_a(:,:)

!     New metric coefficients.  Assumes that the Jacobian (Jx_r) remains unchanged

      do inode = 1,n_pts_3d

         r_x(1,1:3,inode,ielem) = a_t((1-1)*n_pts_3d+inode,1:3) / Jx_r(inode,ielem)
         r_x(2,1:3,inode,ielem) = a_t((2-1)*n_pts_3d+inode,1:3) / Jx_r(inode,ielem)
         r_x(3,1:3,inode,ielem) = a_t((3-1)*n_pts_3d+inode,1:3) / Jx_r(inode,ielem)

      enddo

      !-- testing to see if the discrete metric terms satisfy the correct GCL condition

      if(testing_metric_comp) then
        if(allocated(Hinv_3_mat  )) deallocate(Hinv_3_mat  ); allocate(Hinv_3_mat(m,m)); Hinv_3_mat(:,:) = 0.0_wp
        if(allocated(D_3_mat  )) deallocate(D_3_mat  ); allocate(D_3_mat(m,n)); D_3_mat(:,:) = 0.0_wp
        if(allocated(Q_3_mat  )) deallocate(Q_3_mat  ); allocate(Q_3_mat(m,n)); Q_3_mat(:,:) = 0.0_wp
        if(allocated(S_3_mat  )) deallocate(S_3_mat  ); allocate(S_3_mat(m,n)); S_3_mat(:,:) = 0.0_wp
        if(allocated(E_3_mat  )) deallocate(E_3_mat  ); allocate(E_3_mat(m,n)); E_3_mat(:,:) = 0.0_wp
        if(allocated(verror)) deallocate(verror); allocate(verror(m,3)); verror(:,:) = 0.0_wp
        if(allocated(verror2)) deallocate(verror2); allocate(verror2(m,3)); verror2(:,:) = 0.0_wp
        if(allocated(verror3)) deallocate(verror3); allocate(verror3(m,3)); verror3(:,:) = 0.0_wp


        call GCL_Triple_HinvDQSEmat(n_pts_1d, n_pts_3d, pmat,dmat,qmat,Hinv_3_mat,D_3_mat,Q_3_mat,S_3_mat,E_3_mat)

      endif
    
    end do elloop

    if(allocated(u))    deallocate(u) ;
    if(allocated(v))    deallocate(v) ;
    if(allocated(w))    deallocate(w) ;
    if(allocated(work)) deallocate(work) ;
    if(allocated(Amat)) deallocate(Amat) ;
    if(allocated(a_t))  deallocate(a_t) ;
    if(allocated(bvec)) deallocate(bvec) ;
    if(allocated(wI))   deallocate(wI) ;

!   Parallel reduction of errors using conventional MPI calls

    ! Reduce values on all processes to a single value
    if(myprocid == 0 )  then
      allocate(err_max_proc(0:nprocs-1))
      err_max_proc(:) = 0.0_wp ; err_max_proc(0) = err_Linf ;
    endif
    if(myprocid /= 0 ) then
      s_tag = 100 + myprocid
      m_size = 1
      call mpi_isend(err_Linf,m_size,mpi_double,0,s_tag,petsc_comm_world, &
        & s_request_err_Linf,ierr)
      
      call mpi_wait(s_request_err_Linf,s_status,ierr)
    else
      do m = 1, nprocs-1
        r_tag = 100 + m
        m_size = 1
        call mpi_irecv(err_max_proc(m),m_size,mpi_double,m,r_tag, &
          & petsc_comm_world,r_request_err_Linf,ierr)

        call mpi_wait(r_request_err_Linf,r_status,ierr)
      enddo
    endif

    ! Reduce values on all processes to a single value
    call MPI_reduce(np_mods, ng_mods,1, mpi_integer, mpi_sum, 0, petsc_comm_world, ierr)

    ! Write at screen the L_inf of the metric error
!   if(myprocid == 0 .and. verbose)  then
    if(myprocid == 0)  then
      write(*,*)'  Modified Elements and L_inf of modification', ng_mods, maxval(err_max_proc(:))
      write(*,*) '==========================================================='

      deallocate(err_max_proc)

    endif

  end subroutine modify_metrics_nonConforming

  !============================================================================

  ! e2e_connectivity_aflr3 - Constructs the element-to-element connectivity
  ! starting from the information read from the AFLR3 grid.

  subroutine face_orientation_aflr3()   !  SERIAL Routine

!  SERIAL Routine
!     face orientation between connected elements.  
!  SERIAL Routine

    ! Load modules
    use referencevariables, only : nfacesperelem, nelems
    use variables,          only : ef2e

    implicit none                               ! Nothing is implicitly defined

    integer :: iface, kface
    integer :: ielem, kelem

    do ielem = 1,nelems
       do iface = 1,nfacesperelem
         kelem = ef2e(2,iface,ielem)
         kface = ef2e(1,iface,ielem)
         if((ef2e( 1,iface,ielem) <= 0) .or. &  !  Test for boundary face which is by definition correctly oriented
            (ef2e(10,iface,ielem) <= 0) .or. &  !  Test for periodic face which is by definition correctly oriented
            (ielem == kelem))          then     !  Periodic connection with a single element thickness  (may be redundent)
           ef2e(7,iface,ielem) = 0 ;
           cycle
         else
           ef2e(7,iface,ielem) = face_orientation_hex(ielem,kelem,iface,kface)
         endif
       enddo
    enddo

  end subroutine face_orientation_aflr3

  !============================================================================================

  function face_orientation_hex(elem1,elem2,face1,face2)
     
    !  input : 
    !        elem1, elem2: two adjoining elements 
    !        face1, face2: local ordering for common faces between two adjoining elements 
    !  output:
    !        0 <= face_orientation_hex  <= 3 
    !        0: common ordering  
    !        1: Counter-clockwise rotation 090^o
    !        2: Counter-clockwise rotation 180^o
    !        3: Counter-clockwise rotation 270^o

    use referencevariables
    use variables, only: ic2nh, vx_master

    implicit none

    integer,                   intent(in   ) :: elem1 ,elem2
    integer,                   intent(in   ) :: face1 ,face2

    integer,  parameter                      :: nodetol = 1.0e-10_wp

    real(wp),  dimension(3,4)                :: xface1,xface2

    integer                                  :: j,j1,j2
    integer                                  :: k,kk,ierr
    integer                                  :: face_orientation_hex

!   integer,  dimension(4,8), parameter :: perm = reshape(    &
!                                               & (/1,2,3,4,  &
!                                               &   2,3,4,1,  &
!                                               &   3,4,1,2,  &
!                                               &   4,1,2,3,  &
!                                               &   4,3,2,1,  &
!                                               &   1,4,3,2,  &
!                                               &   2,1,4,3,  &
!                                               &   3,2,1,4/),&
!                                               &   (/4,8/) )
    integer,  dimension(4,8), parameter :: perm = reshape(    &
                                                & (/1,2,3,4,  &
                                                &   3,1,4,2,  &
                                                &   4,3,2,1,  &
                                                &   2,4,1,3,  &
                                                &   1,3,2,4,  &
                                                &   2,1,4,3,  &
                                                &   4,2,3,1,  &
                                                &   3,4,1,2/),&
                                                &   (/4,8/) )

       kk = hex_8_nverticesperface
                                                     !   face of element 1 
       do j = 1, kk                                  !   Sweep over vertices on face
         j1 = eltypfaces_Lexo(j,face1)               !   Which node in element
         j2 = ic2nh(j1,elem1)                        !   Which vertex is pointed to
         xface1(:,j) = vx_master(:,j2)               !   load face nodes with vertex coordinates
       enddo
                                                     !   face of element 2 
       do j = 1, kk                                  !   Sweep over vertices on face
         j1 = eltypfaces_Lexo(j,face2)               !   Which node in element
         j2 = ic2nh(j1,elem2)                        !   Which vertex is pointed to
         xface2(:,j) = vx_master(:,j2)               !   load face nodes with vertex coordinates
       enddo

       do k = 1,8
          if ( (magnitude(xface1(:,1) -xface2(:,perm(1,k))) <= nodetol)  .and. & ! Check distances between all nodes 
               (magnitude(xface1(:,2) -xface2(:,perm(2,k))) <= nodetol)  .and. &
               (magnitude(xface1(:,3) -xface2(:,perm(3,k))) <= nodetol)  .and. &
               (magnitude(xface1(:,4) -xface2(:,perm(4,k))) <= nodetol) ) then  
               face_orientation_hex = k-1
               exit
          endif
       end do
       
       if(k > 8) then
          write(*,*)'function face_orientation_hex: didnt find a face orientation that matched'
          write(*,*)'stopping'
          call PetscFinalize(ierr) ; stop
       endif
          
  end function face_orientation_hex

  !============================================================================================
  !
  ! Purpose: For a given face node number on the elemenet, k_On in local face ordering, this function
  !          returns the face node number of the adjoining element, k_Off in local face ordering.
  !
  ! Inputs k_On = the face node number in the local face ordering of the on element
  !        orientation = orientation number
  !                   0 = normal orientation
  !                   1 = counter-clock-wise rotation of theta = 3*pi/2
  !                   2 = counter-clock-wise rotation of theta = pi
  !                   3 = counter-clock-wise rotation of theta = pi/2
  !                   4 = counter-clock-wise rotation of theta = 3*pi/2 with a switch of axis
  !                   5 = switch of axis
  !                   6 = counter-clock-wise rotation of theta = pi/2 with a switch of axis
  !                   7 = counter-clock-wise rotation of theta = pi with a switch of axis
  !       n_1d_On = number of nodes in each direction on the On element
  !
  ! Outputs  map_face_orientation_k_On_2_k_Off = the face node number on the off 
  !          element associated with k_On
  ! 
  ! Notes: This function is constructed based on the idea that k_On map onto 
  !        (i_On,j_On) which can then be mapped to (i_Off,j_Off) which can then 
  !        be mapped to k_Off
  ! The map between k_On-->(i_On, j_On) and (i_Off, j_Off) --> k_Off are highlighted 
  ! are given below. Here the discussion focuses on the map between 
  ! (i_On, j_On)--> (i_Off,j_Off) was constructed 
  !  
  ! This maping can be thought of as mapping between integer sets and can be arrived 
  ! at by using a rotation about the point ( (n_1d_On-1)/2+1, (n_1d_On-1)/2 +1) 
  ! and a switch of coordinates
  !
  ! For the rotation we have
  !
  ! [i_Off; j_Off;1] =| a11     (1-a11)  0  |   |  1  0  (n1d-1)/2+1  | |  c  s  0  | |  1  0  -(n1d-1)/2-1  | [i_On; j_Off; 1]
  !                   | (1-a11)  a11     0  |   |  0  1  (n1d-1)/2+1  | |  -s c  0  | |  0  1  -(n1d-1)/2-1  | 
  !                   |  0       0       1  |   |  0  0       1       | |  0  0  1  | |  0  0       1        |
  !
  !  The first matrix will change the coordinate axis, i.e., i_Off = j_On, j_Off = i_On for a11 = 0, while for a11 = 1 it is an identity matrix 
  !  The next three matricies effect a rotation in the counter-clockwise direction of theta, where c = cos(theta), and s = sin(theta)
  !  about the point ( (n_1d_On-1)/2+1, (n_1d_On-1)/2 +1).
  !  The commented out code gives the resulting mappings for the various choices of orientation. 
  !  We can however remove the if statments in that version of the function by solving for a11, c, and s as degree 7 polynomials of orientation, 
  !  that is for example a11 = \sum_{i=0}^{7}a_{i}*orientation**i and we solve for the a_{i}'s so that they give the correct value for orientation.
  ! 
  ! Tests: Tested in MatLab script, see unit_tests/test_map_face_orientation_compact
  !        Tested in code with orientation = 0 
  !
  !=============================================================================

  function map_face_orientation_k_On_2_k_Off(k_On,orientation,n_1d_On)

     implicit none

     integer, intent(in)                         :: k_On, orientation, n_1d_On

     integer                                     :: map_face_orientation_k_On_2_k_Off

     real(wp)                                    :: i_On2, j_On2, i_Off2, j_Off2

     real(wp)                                    :: a11, a0, a1, a2, a3, a4, a5, a6 ,a7
     real(wp)                                    :: c, c0, c1, c2, c3, c4, c5, c6, c7
     real(wp)                                    :: s, s1, s2, s3, s4, s5, s6, s7

    !-- map between k_On and (i_On2, j_On2)
    i_On2 = k_On-floor((k_On-1.0_wp)/n_1d_On)*n_1d_On
    j_On2 = (k_On-i_On2)/n_1d_On+1.0_wp

    a0 = 1.0_wp
    a1 = 2341.0_wp/420.0_wp
    a2 = -931.0_wp/72.0_wp
    a3 = 791.0_wp/72.0_wp
    a4 = -161.0_wp/36.0_wp
    a5 = 337.0_wp/360.0_wp
    a6 = -7.0_wp/72.0_wp
    a7 = 1.0_wp/252.0_wp
    a11 = a0 + a1*orientation + a2*orientation**2 + a3*orientation**3 + a4*orientation**4&
          + a5*orientation**5 + a6*orientation**6 + a7*orientation**7

    c0 = 1.0_wp
    c1 = 335.0_wp/28.0_wp
    c2 = -2093.0_wp/72.0_wp
    c3 = 851.0_wp/36.0_wp
    c4 = -665.0_wp/72.0_wp
    c5 = 17.0_wp/9.0_wp
    c6 = -7.0_wp/36.0_wp
    c7 = 1.0_wp/126.0_wp
    c = c0 + c1*orientation + c2*orientation**2 + c3*orientation**3 + c4*orientation**4 + c5*orientation**5&
           + c6*orientation**6 + c7*orientation**7

    s1 = 49.0_wp/4.0_wp;
    s2 = -11837.0_wp/360.0_wp
    s3 = 5333.0_wp/180.0_wp
    s4 = -889.0_wp/72.0_wp
    s5 = 47.0_wp/18.0_wp
    s6 = -49.0_wp/180.0_wp
    s7 = 1.0_wp/90.0_wp
    s = s1*orientation + s2*orientation**2 + s3*orientation**3 + s4*orientation**4 + s5*orientation**5 &
        + s6*orientation**6 + s7*orientation**7

    i_Off2 = a11 * c * i_On2 + a11 * i_On2 * s - i_On2 * s - a11 * c * j_On2 + a11 * j_On2 * s + c * j_On2 &
            - n_1d_On * a11 * s + n_1d_On * s / 2.0_wp - a11 * s + s / 2.0_wp - n_1d_On * c / 2.0_wp &
            - c / 2.0_wp + n_1d_On / 2.0_wp + 0.5_wp
    j_Off2 = -a11 * c * i_On2 - a11 * i_On2 * s + c * i_On2 + a11 * c * j_On2 - a11 * j_On2 * s + j_On2 * s &
            + n_1d_On * a11 * s - n_1d_On * c / 2.0_wp + a11 * s - c / 2.0_wp - n_1d_On * s / 2.0_wp &
            - s / 2.0_wp + n_1d_On / 2.0_wp + 0.5_wp

    !-- this is k_Off
    map_face_orientation_k_On_2_k_Off = NINT((j_Off2-1)*n_1d_On+i_Off2)


!-- simple code with ifs
!    implicit none
!    
!    integer, intent(in)                         :: k_On, orientation, n_1d_On
!
!    integer                                     :: map_face_orientation_k_On_2_k_Off
!    integer                                     :: i_On, j_On, i_Off, j_Off
!
!    !-- map between k_On and (i_On, j_On)
!    i_On = k_On-floor((k_On-1.0_wp)/n_1d_On)*n_1d_On
!    j_On = (k_On-i_On)/n_1d_On+1
!
!    if(orientation.EQ.0)then
!      i_Off = i_On
!      j_Off = j_On
!    elseif(orientation.EQ.1)then
!      i_Off = -j_On+n_1d_On +1
!      j_Off = i_On    
!    elseif(orientation.EQ.2)then
!      i_Off = -i_On+n_1d_On+1
!      j_Off = -j_On+n_1d_On+1
!    elseif(orientation.EQ.3)then
!      i_Off = j_On
!      j_Off = -i_On+n_1d_On +1
!    elseif(orientation.EQ.4)then
!      i_Off = i_On;
!      j_Off = -j_On+n_1d_On +1;     
!    elseif(orientation.EQ.5)then
!      i_Off = j_On;
!      j_Off = i_On;
!    elseif(orientation.EQ.6)then
!      i_Off = -i_On+n_1d_On +1;
!      j_Off = j_On;
!    elseif(orientation.EQ.7)then
!      i_Off = -j_On+n_1d_On+1;
!      j_Off = -i_On+n_1d_On+1;
!    endif
!
!    !-- this is k_Off, i.e., the map from (i_Off, j_Off) and k_Off
!    map_face_orientation_k_On_2_k_Off = (j_Off-1)*n_1d_On+i_Off

  end function map_face_orientation_k_On_2_k_Off

  subroutine Load_Mortar_Metric_Data(ielem,n_LGL_1d,n_LGL_2d,n_LGL_3d,ifacenodes,p_surf,a_t,bvec)

    use referencevariables
    use variables, only: Jx_r, r_x, facenodenormal, Jx_facenodenormal_LGL, ef2e
    use collocationvariables, only: elem_props

    implicit none

    integer,                   intent(in   ) :: ielem,n_LGL_1d,n_LGL_2d,n_LGL_3d
    integer,   dimension(:),   intent(in   ) :: ifacenodes
    real(wp),  dimension(:),   intent(in   ) :: p_surf

    real(wp), dimension(:,:),  intent(inout) :: a_t
    real(wp), dimension(:,:),  intent(inout) :: bvec

    integer                                  :: i, iface, inode, jnode, knode

    real(wp), dimension(3)                   :: nx


    bvec(:,:) = 0.0_wp
    do iface = 1,nfacesperelem

       do i = 1,n_LGL_2d

            !knode = (iface-1)*(npoly_max+1)**2 + i
             knode = (iface-1)*n_LGL_2d + i
            ! Index in facial ordering
            jnode =  n_LGL_2d*(iface-1)+i
              
            ! Volumetric node index corresponding to facial node index
            inode = ifacenodes(jnode)
              
            ! Outward facing normal of facial node
            !if (ef2e(4,iface,ielem) == elem_props(2,ielem)) then !    Conforming interface-- ORIGINAL
            if ((ef2e(4,iface,ielem) == elem_props(2,ielem)).OR.(ef2e(1,iface,ielem) < 0)) then !    Conforming interface or a boundary
              nx(:) = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)
!             t1 = maxval(abs(nx(:) - Jx_facenodenormal_LGL(:,knode,ielem)))
!             if(t1 >= 1.0e-10_wp) write(*,*)'metric differences: iface',iface, t1
            else                                                 ! NonConforming interface
              nx(:) = Jx_facenodenormal_LGL(:,knode,ielem)
            endif

            bvec(inode,:) = bvec(inode,:) + p_surf(i)*nx(:)

       enddo

    enddo

!        Metric data is stored as follows
!              --                                  --
!              | dxi_1/dx_1, dxi_1/dx_2, dxi_1/dx_3 |
!              |                                    |
!        r_x = | dxi_2/dx_1, dxi_2/dx_2, dxi_2/dx_3 |
!              |                                    |
!              | dxi_3/dx_1, dxi_3/dx_2, dxi_3/dx_3 |
!              --                                  --
   
    do inode = 1,n_LGL_3d
      a_t((1-1)*n_LGL_3d + inode,1:3) =   Jx_r(inode,ielem)*r_x(1,1:3,inode,ielem)
      a_t((2-1)*n_LGL_3d + inode,1:3) =   Jx_r(inode,ielem)*r_x(2,1:3,inode,ielem)
      a_t((3-1)*n_LGL_3d + inode,1:3) =   Jx_r(inode,ielem)*r_x(3,1:3,inode,ielem)
    enddo

  end subroutine Load_Mortar_Metric_Data

  !============================================================================

  ! e2e_connectivity_aflr3 - Constructs the element-to-element connectivity
  ! starting from the information read from the AFLR3 grid.

  subroutine e2e_connectivity_aflr3()   !  SERIAL Routine

!  SERIAL Routine
!     Element to element connectivity using AFLR3 formated grids
!  SERIAL Routine

    ! Load modules
    use unary_mod, only: qsorti

    use referencevariables
    use variables, only: ef2e, iae2v, iae2v_tmp, jae2v, jae2v_tmp, nnze2v, &
                       & nnze2e, iae2e,jae2e, &
                       & nnze2e2e, iae2e2e,jae2e2e, &
                       & if2nq, ifacetag, ic2nh, nqface, vx_master, &
                       & periodic_face_data_x1, periodic_face_data_x2, &
                       & periodic_face_data_x3, wall_face_data

    ! Nothing is implicitly defined
    implicit none

    integer :: i,j,k, j1, j2, k1
    integer :: cnt,cntBC 
    integer :: ielem,jelem, kelem
    integer :: iface, jface, bigN, ave

    !integer, dimension(:),   allocatable :: iav2e_tmp
    integer, dimension(:,:), allocatable :: jav2e_tmp
    integer, dimension(:),   allocatable :: ivtmp1, ivtmp2, ivtmp3, ivtmp4
    integer, dimension(:),   allocatable :: jae2e2e_tmp

    integer :: nnzv2e
    !integer, allocatable :: v2e(:,:), iav2e(:), jav2e(:)
    integer, allocatable :: iav2e(:), jav2e(:)

    integer, dimension(:,:), allocatable :: test_conL, test_conR
    integer, dimension(:),   allocatable :: test_cBCL, test_cBCR
    integer, dimension(:),   allocatable :: stack, ordered

    logical                              :: testing = .true.

    ! Variables for periodic faces
    integer,  dimension(:),     allocatable :: ivtmp_iface, ivtmp_jface
    real(wp), dimension(:,:),   allocatable :: vx_iface,    vx_jface
    integer,  dimension(:,:),     allocatable :: same_coord 
    integer :: i_comp_face, i_vertex_iface, i_vertex_jface, i_coord, i_p_elem, &
      & j_p_elem, match, i_vertex_hexa, cnt_periodic_elem_x1, &
      & cnt_periodic_elem_x2, cnt_periodic_elem_x3

    integer :: cnt_wall_elem

    integer,  dimension(:),   allocatable :: iv_hexa_elem 

    integer :: jelem_iface, jface_iface, ielem_jface, iface_jface, p_face, &
      & i_check 

    !integer, allocatable, dimension(:) :: list_partner_faces

    real(wp), parameter  :: diff_toll = 1e-8
    integer,  parameter  :: qdim = 10           !  dimension of ef2e array

    continue

    !  routine called from myprocid == 0 ;  i.e. master works on entire grid
    !
    !  AFLR3 grid format includes the following variables
    !
    !  nhex       = number of hexahedra
    !  nnodesg    = number of nodes
    !  ic2nh      = cell-to-node pointers for hexes 
    !               (8 nodes are pointed to by each hexahedral cell)

    !  nqface     = number of boundary faces with quads
    !  if2nq      = face-to-node pointers for each QUAD boundary face (4 nodes
    !               are pointed to by each quadralateral)
    !  ifacetag   = tag number that groups boundary faces (e.g., all faces
    !               with the number "1" correspond to the same boundary
    !               entity (like wing) with a given boundary condition)

    !
    !  nelems     :    elements  =  nhex in this case
    !
    !                   Dim,    Dim,         Dim
    !  ef2e       :    (qdim, nfaceperelem, nelements) 
    !             :  Two situation occur.  The face is either an 
    !                  (Interior face 
    !                      :  (1,j,k) = Adjoining element face ID
    !                      :  (2,j,k) = Adjoining element ID
    !                      :  (3,j,k) = Adjoining element process ID
    !                      :  (4,j,k) = Adjoining element polynomial order
    !                      :  (5,j,k) = Number of Adjoining elements
    !                      :  (6,j,k) = HACK self polynomial order assigned to each face
    !                      :  (7,j,k) = face_orientation
    !                      :  (8,1,k) = which entries in parent_geo(1:3,1:8,ef2e(8,1,k)) to get the parent vertices
    !                      :  (8,2,k) = local element number for refined elements (1:8)  
    !                      :  (9,j,k) = if 0 then conforming face if 1
    !                      : (10,j,k) = carries the original BC tag for a periodic face
    !  !      do iface = 1, nfacesperelemnonconforming face
    !                  (Boundary face 
    !                      :  (1,j,k) = Set to -11 
    !                      :  (2,j,k) = -100000000
    !                      :  (3,j,k) = -100000000
    !                      :  (4,j,k) = -100000000
    !                      :  (5,j,k) = -100000000
    !                      :  (6,j,k) = -100000000
    !                      :  (7,j,k) = -100000000
    !
    ! iae2v,jae2v     :    Which vertices belong to each element

    if (ndim == 2) then
      write(*,*) '2D grid connectivity for the AFLR3 format is not yet', &
        & ' implemented'
      write(*,*) 'Exiting...'
      stop

    else if (ndim == 3) then
      ! Allocate memory for test
      allocate(test_conL(nelems,nfacesperelem)) ; test_conL(:,:) = 0
      allocate(test_conR(nelems,nfacesperelem)) ; test_conR(:,:) = 0
      allocate(test_cBCL(nqface)) ;               test_cBCL(:) = 0
      allocate(test_cBCR(nqface)) ;               test_cBCR(:) = 0

      ! Big number
      bigN = 1000

      ! Print at screen: number of vertices, number of hexhedral elements,
      ! number of boundary faces and number of vertices per element
      write(*,*) '  Number of vertices             = ', nvertices
      write(*,*) '  Number of hexahedrons          = ', nelems
      write(*,*) '  Number of boundary faces       = ', nqface
      write(*,*) '  Number of vertices per element = ', nverticesperelem
      write(*,*) '=========================================================================='

      ! iav2e,  jav2e_tmp: vessels to accumulate v2e pointer information
      ! First:  Count total accesses to each vertex from all elements

      allocate(iav2e(nvertices+1))        ;    iav2e(:)       = 0
      allocate(jav2e_tmp(nvertices,bigN)) ;    jav2e_tmp(:,:) = 0

      ! Count the number of accesses for each vertex over all elements
      ! iav2e    :  counts total number of touches to each vertex
      ! jav2e_tmp:  tracks elements that touch vertices
      do j = 1,nelems
        do i = 1,nverticesperelem
            iav2e(ic2nh(i,j)) = iav2e(ic2nh(i,j)) + 1
            jav2e_tmp(ic2nh(i,j),iav2e(ic2nh(i,j))) = j
        end do
      end do

      nnzv2e = sum(iav2e)
      allocate(jav2e(nnzv2e))

      ! Reassemble in standard CSR form
      cnt = 0
      do i = 1, nvertices
        do j = 1,iav2e(i)
          cnt = cnt + 1
          jav2e(cnt) = jav2e_tmp(i,j)
        enddo
      enddo

      do i = nvertices,1,-1
        iav2e(i+1) = iav2e(i)
      end do
      iav2e(1) = 1
      do i = 2,nvertices+1
        iav2e(i) = iav2e(i) + iav2e(i-1)
      enddo
      
      deallocate(jav2e_tmp)

      ! Calculate element-to-element connectivity using shared nodes
      allocate(ef2e(qdim,nfacesperelem,1:nelems))  ;   ef2e(:,:,:) = -1000000000

      allocate(ivtmp1(nverticesperface*bigN),ivtmp2(nverticesperface*bigN))
      allocate(ivtmp3(nverticesperface),     ivtmp4(nverticesperface))

      allocate(ivtmp_iface(nverticesperface))
      allocate(ivtmp_jface(nverticesperface))

      allocate(vx_iface(3,nverticesperface))
      allocate(vx_jface(3,nverticesperface))
      allocate(same_coord(nverticesperface,3))
      allocate(iv_hexa_elem(nverticesperelem))

                                                        ! Counter for the number of elements that have a periodic boundary face
      cnt_periodic_elem_x1 = 0                          ! x1 direction
      cnt_periodic_elem_x2 = 0                          ! x2 direction
      cnt_periodic_elem_x3 = 0                          ! x3 direction

      cnt_wall_elem = 0                                 ! Counter for the number of elements that have a wall boundary face

      cntBC = 0                                         ! Counter for the number of boundary conditions

      do ielem = 1,nelems

        faceloop: do iface = 1, nfacesperelem

          cnt = 0
          do j = 1, nverticesperface                    !   Sweep over vertices on face
            j1 = eltypfaces(j,iface)                    !   Which node in element
            j2 = ic2nh(j1,ielem)                        !   Which vertex is pointed to
            do k = iav2e(j2),iav2e(j2+1)-1              !   Sweep over all elements touching vertex &
              cnt = cnt + 1                             !   & and put them in a big stack
              ivtmp1(cnt) = jav2e(k)
            end do
          end do

          call qsorti(ivtmp1,ivtmp2,cnt)                ! Sort the stack from all facenodes
                                                       
          do j = nverticesperface,cnt                   ! search for duplicate entries = nverticesperface 

            ivtmp3(:) = ivtmp1(ivtmp2(j-nverticesperface+1:j)) ! test nverticesperface terms
            ave = sum(ivtmp3(:))/nverticesperface

            if(maxval(abs(ivtmp3(:)-ave)) == 0) then         !  Two adjacent terms match in stack
              if( ave /= ielem ) then
                ef2e(2,iface,ielem) = ave

                cnt = 0 ; ivtmp3(:) = 0                 ! Store vertices on local face
                do k = 1, nverticesperface              ! Sweep over vertices on face
                  cnt = cnt + 1
                  k1  = eltypfaces(k,iface)             !   Which node in element
                  ivtmp3(cnt) = ic2nh(k1,ielem)         !   Which vertex is pointed to
                enddo
                ivtmp3 = isort(ivtmp3,nverticesperface)

                do jface = 1,nfacesperelem              ! Sweep over all faces on adjoining element

                  cnt = 0 ; ivtmp4(:) = 0
                  do k = 1, nverticesperface            ! Sweep over vertices on face
                    cnt = cnt + 1
                    k1  = eltypfaces(k,jface)           !   Which node in element
                    ivtmp4(cnt) = ic2nh(k1,ave)         !   Which vertex is pointed to
                  enddo
                  ivtmp4 = isort(ivtmp4,nverticesperface)

                  if(maxval(abs(ivtmp3(:)-ivtmp4(:))) == 0)  then
                    ef2e(1,iface,ielem) = jface
                    test_conL(ielem,iface) = ielem*iface ; test_conR(ave ,jface) = ave*jface
                    cycle faceloop
                  endif

                enddo

              endif
            endif

          end do 

!===================
!         Begin periodic face and solid boundary face connectivity loop
!===================

          cnt = 0 ; ivtmp3(:) = 0                       ! Stack for vertices on face (iface) of element (ielem)
          do j = 1, nverticesperface
            cnt = cnt + 1                               ! index counter
            j1 = eltypfaces(j,iface)                    ! Which node in element
            ivtmp3(cnt) = ic2nh(j1,ielem)               ! Which vertex is pointed to
          end do
          ivtmp3 = isort(ivtmp3,nverticesperface)       ! Sort the stack

          do jface = 1,nqface                           ! Loop over all boundary faces

            ivtmp4(:) = 0 ; ivtmp4(:) = if2nq(:,jface)  ! Stack for vertices of BC face (jface) 
            ivtmp4 = isort(ivtmp4,nverticesperface)     ! Sort the stack

            if(maxval(abs(ivtmp3(:)-ivtmp4(:))) == 0)  then            ! Found a match between element face and BC face stack

              ef2e( 1,iface,ielem) = -ifacetag(jface)                  ! Record BC of face into ef2e(1) with a negative sign (denotes BC)
              ef2e(10,iface,ielem) = -ifacetag(jface)                  ! Record BC of face into ef2e(1) with a negative sign (denotes BC)

              if(ifacetag(jface) == 8 .or. ifacetag(jface) == 9) then  ! Build periodic faces data in the x1 direction

                cnt_periodic_elem_x1 = cnt_periodic_elem_x1 + 1        ! Counter for periodic connections in X1 direction

                periodic_face_data_x1(1,cnt_periodic_elem_x1) = ielem  ! Global ID of the element which owns a "periodic" boundary face
                periodic_face_data_x1(2,cnt_periodic_elem_x1) = iface  ! Get the local ID of the "periodic" boundary face of the ielem
                periodic_face_data_x1(3,cnt_periodic_elem_x1) = jface  ! Position in the if2nq stack of the iface

                do j = 1, nverticesperface
                  periodic_face_data_x1(4+j,cnt_periodic_elem_x1) = ivtmp4(j)  ! store vertices of face in X1 direction
                end do

              end if

              if(ifacetag(jface) == 10 .or. ifacetag(jface) == 11) then ! Build periodic faces data in the x2 direction

                cnt_periodic_elem_x2 = cnt_periodic_elem_x2 + 1         ! Counter for periodic connections in X1 direction

                periodic_face_data_x2(1,cnt_periodic_elem_x2) = ielem   ! Global ID of the element which owns a "periodic" boundary face
                periodic_face_data_x2(2,cnt_periodic_elem_x2) = iface   ! Get the local ID of the "periodic" boundary face of the ielem
                periodic_face_data_x2(3,cnt_periodic_elem_x2) = jface   ! Position in the if2nq stack of the iface

                do j = 1, nverticesperface
                  periodic_face_data_x2(4+j,cnt_periodic_elem_x2) = ivtmp4(j)  ! store vertices of face in X1 direction
                end do

              end if

              if(ifacetag(jface) == 12 .or. ifacetag(jface) == 13) then ! Build periodic faces data in the x3 direction

                cnt_periodic_elem_x3 = cnt_periodic_elem_x3 + 1         ! Counter for periodic connections in X1 direction

                periodic_face_data_x3(1,cnt_periodic_elem_x3) = ielem   ! Global ID of the element which owns a "periodic" boundary face
                periodic_face_data_x3(2,cnt_periodic_elem_x3) = iface   ! Get the local ID of the "periodic" boundary face of the ielem
                periodic_face_data_x3(3,cnt_periodic_elem_x3) = jface   ! Position in the if2nq stack of the iface

                do j = 1, nverticesperface
                  periodic_face_data_x3(4+j,cnt_periodic_elem_x3) = ivtmp4(j)  ! store vertices of face in X1 direction
                end do

              end if

              if(ifacetag(jface) == 5 .or. ifacetag(jface) == 6) then   ! Store data for computing the aerodynamic coefficients

                cnt_wall_elem = cnt_wall_elem + 1

                wall_face_data(1,cnt_wall_elem) = ielem                 ! Global ID of the element which owns a "wall" boundary face
                wall_face_data(2,cnt_wall_elem) = iface                 ! Get the local ID of the "wall" boundary face of the ielem
                
              end if

              cntBC = cntBC + 1 ; test_cBCL(cntBC) = jface

              cycle faceloop                                            ! Found a face match.  Cycle by one one on the face loop

            endif
          enddo

          write(*,*) 'Face search failed. element and face' ;           ! Search should always find a connection but didn't. Somethings wrong
          write(*,*) ielem,iface ; stop

        end do faceloop
      end do

      ! Test if the element are connected in a symmetric fashion
      ! ========================================================
      if(testing) then
        if(maxval(abs(test_conL(:,:)-test_conR(:,:))) > 0) then
          write(*,*)'Not symmetric connectivity' ; write(*,*)'Stopping' ; stop
        endif
        call qsorti(test_cBCL,test_cBCR,nqface)
        do i = 1,nqface
          if(i /= test_cBCL(test_cBCR(i))) then
            write(*,*)i,test_cBCL(i),test_cBCR(i)
            stop
          endif
        enddo
      endif

      ! ===================================================
      ! ===================================================
      ! Build ef2e for "periodic" faces in the x1 direction
      ! ===================================================
      ! ===================================================

      periodic_elem_x1_1 : do i_p_elem = 1, size(periodic_face_data_x1(1,:))   ! Loop over the element that owns a "periodic" face

        ielem = periodic_face_data_x1(1,i_p_elem)                              ! Get the global ID of the element that owns a periodic boundary face
        iface = periodic_face_data_x1(2,i_p_elem)                              ! Get the local (local for the element) ID of the periodic boundary face

        if (abs(ef2e(1,iface,ielem)) == 8) then                                ! Check if the face is a "periodic" face with tag equals to 8

          do i_vertex_iface = 1, nverticesperface                              ! Get the coordinate of the nodes of the iface
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
          end do

          partner_face_x1_1: do i_comp_face = 1, size(periodic_face_data_x1(3,:)) ! Partner face search, i.e. the corresponding "periodic" face with face tag 9
            
            jface = periodic_face_data_x1(3,i_comp_face)                        ! Get the face ID. This ID is given by the AFLR3 format
            
            if (ifacetag(jface) == 9) then                                      ! Partner face can only be a face with ifacetag = 9
                  
              ivtmp_jface(:) = if2nq(:,jface)                                   ! Get the ID of the nodes which form the boundary face

              do i_vertex_jface = 1, nverticesperface                           ! Get the coordinate of the nodes of the jface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

              same_coord(1:nverticesperface,1:3) = 0                                                    ! Set to zero the array that keeps track of the coordinates matches

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs by a shift.

              search_x1_1: do i_vertex_iface = 1, nverticesperface

                do i_vertex_jface = 1, nverticesperface

                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1
                    endif
                  end do
                  
                  if (sum(same_coord(i_vertex_iface,:)) == 2) then             ! Check if vertex on both iface and jface have two equal coordinates
                    if (i_vertex_iface .gt. 1) then                            ! Verify that the other possible partner nodes have the same two coordinates in common'
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x1_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x1_1
                          else
                            write(*,*) 'The periodic boundary conditions', &
                              & ' works only for two parallel boundary planes.', &
                              & ' Check x1 direction.'
                            write(*,*) 'Exiting...'
                            stop
                          end if
                        end do check_common_coord_x1_1
                      end do
                      cycle search_x1_1
                    end if
                    cycle search_x1_1
                  end if
                  
                end do

                if (sum(same_coord(i_vertex_iface,:)) .lt. 2) then              ! Exit from the companion_face loop because none of the nodes of the jface has two invariant coordinates
                  if (i_comp_face == size(periodic_face_data_x1(3,:))) then
                    write(*,*) 'No periodic partner face has been found for', &
                      & ' face', iface, 'of element', ielem, 'in the x1', &
                      & ' direction.'
                    write(*,*) 'Try to increase the diff_toll parameter and', &
                      & ' verify the input grid'
                    write(*,*) 'Exiting...'
                    stop
                  else
                    cycle partner_face_x1_1
                  end if
                endif

              end do search_x1_1
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                if (same_coord(1,i_coord) == 0) then                           ! Just pick up the first node for this check.
                  periodic_face_data_x1(4,i_p_elem) = i_coord
                end if
              end do     

              search_x1_2 : do j_p_elem = 1, size(periodic_face_data_x1(1,:))  ! Search the global ID of the element which owns the jface
                
                jelem = periodic_face_data_x1(1,j_p_elem)                      ! Get the global ID of the element

                iv_hexa_elem(:) = ic2nh(:,jelem)                               ! Get the ID of the nodes which form the jelem (nverticesperelem)

                ! We have nverticesperface so we should find nverticesperface
                ! matches to claim that we have found the right (not self) element
                match = 0
                do i_vertex_jface = 1, nverticesperface
                  do i_vertex_hexa = 1, nverticesperelem
                    if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
                      match = match + 1
                    endif
                  end do
                end do
                  
                ! Set in ef2e(2,iface,ielem) the ID of the adjoining element 
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

                  do p_face = 1, nfacesperelem                                 ! Set the ID of the adjoining face

                    cnt = 0 ; ivtmp4(:) = 0
                    do k = 1, nverticesperface            
                      cnt = cnt + 1
                      k1  = eltypfaces(k,p_face)           
                      ivtmp4(cnt) = ic2nh(k1,jelem)         
                    enddo
                    ivtmp4 = isort(ivtmp4,nverticesperface)

                    ivtmp_jface = isort(ivtmp_jface,nverticesperface)           ! Sort ivtmp_jface

                    if(maxval(abs(ivtmp_jface(:)-ivtmp4(:))) == 0)  then
                      ef2e(1,iface,ielem) = p_face
                      exit
                    endif

                  end do

                else
                  cycle search_x1_2
                endif
                exit
              enddo search_x1_2
              exit
            endif
          end do partner_face_x1_1

        endif

      end do periodic_elem_x1_1

      periodic_elem_x1_2 : do i_p_elem = 1, size(periodic_face_data_x1(1,:))

        ielem = periodic_face_data_x1(1,i_p_elem)                                 ! Get the global ID of the element that owns a boundary face
        iface = periodic_face_data_x1(2,i_p_elem)                                 ! Get the local (local for the element) ID of the boundary face

        if (abs(ef2e(1,iface,ielem)) == 9) then                                   ! Check if the face is a "periodic" face

          do i_vertex_iface = 1, nverticesperface                                 ! Get the coordinate of the nodes of the iface
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
          end do

          partner_face_x1_2: do i_comp_face = 1, size(periodic_face_data_x1(3,:)) ! Search the partner face, i.e. the corresponding "periodic" face
            
            jface = periodic_face_data_x1(3,i_comp_face)                          ! Get the face ID (the one given by the AFLR3 format)
            
            if (ifacetag(jface) == 8) then                                        ! Partner face can only be a face with ifacetag = 8
                  
              ivtmp_jface(:) = if2nq(:,jface)                                     ! Get the ID of the nodes which form the boundary face

              do i_vertex_jface = 1, nverticesperface                             ! Get the coordinate of the nodes of the jface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

              same_coord = 0                                                      ! Set to zero the array that keeps track of the coordinates matches

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs by a shift.

              search_x1_3: do i_vertex_iface = 1, nverticesperface
                do i_vertex_jface = 1, nverticesperface
                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1 
                    endif
                  end do

                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the same two coordinates in common
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x1_2: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x1_2
                          else
                            write(*,*) 'The periodic boundary conditions', &
                              & ' works only for two parallel boundary planes.', &
                              & ' Check x1 direction.'
                            stop
                          end if
                        end do check_common_coord_x1_2
                      end do

                      cycle search_x1_3
                    endif
                    cycle search_x1_3
                  endif
                  
                end do

                if (sum(same_coord(i_vertex_iface,:)) .lt. 2) then
                  ! Exit from the companion_face loop because none of the nodes
                  ! of the jface has two invariant coordinates
                  if (i_comp_face == size(periodic_face_data_x1(3,:))) then
                    write(*,*) 'No periodic partner face has been found for', &
                      & ' face', iface, 'of element', ielem, 'in the x1', &
                      & ' direction.'
                    write(*,*) 'Try to increase the diff_toll parameter and', &
                      & ' verify the input grid'
                    write(*,*) 'Exiting...'
                    stop
                  else
                    cycle partner_face_x1_2
                  end if
                endif

              end do search_x1_3
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x1(4,i_p_elem) = i_coord
                end if
              end do

              search_x1_4 : do j_p_elem = 1, size(periodic_face_data_x1(1,:)) ! Search the global ID of the element which owns the jface
                
                jelem = periodic_face_data_x1(1,j_p_elem)                     ! Get the global ID of the element
                iv_hexa_elem = ic2nh(:,jelem)                                 ! Get the ID of the nodes which form the jelem

                ! We have nverticesperface so we should find nverticesperface
                ! matches to claim that we have found the right (not self) element
                match = 0
                do i_vertex_jface = 1, nverticesperface
                  do i_vertex_hexa = 1, nverticesperelem
                    if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
                      match = match + 1
                    endif
                  end do
                end do
                  
                ! Set in ef2e(2,iface,ielem) the ID of the adjoining element 
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

                  ! Set the ID of the adjoining face
                  do p_face = 1, nfacesperelem              

                    cnt = 0 ; ivtmp4(:) = 0
                    do k = 1, nverticesperface            
                      cnt = cnt + 1
                      k1  = eltypfaces(k,p_face)           
                      ivtmp4(cnt) = ic2nh(k1,jelem)         
                    enddo
                    ivtmp4 = isort(ivtmp4,nverticesperface)

                    ! Sort ivtmp_jface
                    ivtmp_jface = isort(ivtmp_jface,nverticesperface)

                    if(maxval(abs(ivtmp_jface(:)-ivtmp4(:))) == 0)  then
                      ef2e(1,iface,ielem) = p_face
                      exit
                    endif

                  end do
                else
                  cycle search_x1_4
                endif
                  exit
              enddo search_x1_4
              exit
            endif
          end do partner_face_x1_2

        endif

      end do periodic_elem_x1_2

      ! ===================================================
      ! ===================================================
      ! Build ef2e for "periodic" faces in the x2 direction
      ! ===================================================
      ! ===================================================

      ! Loop over the element that owns a "periodic" face
      periodic_elem_x2_1 : do i_p_elem = 1, size(periodic_face_data_x2(1,:))

        ! Get the global ID of the element that owns a periodic boundary face
        ielem = periodic_face_data_x2(1,i_p_elem)

        ! Get the local (local for the element) ID of the periodic boundary face
        iface = periodic_face_data_x2(2,i_p_elem)

        ! Check if the face is a "periodic" face with tag equals to 8
        if (abs(ef2e(1,iface,ielem)) == 10) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

          ! Search the partner face, i.e. the corresponding "periodic" face with
          ! face tag equals to 9
          partner_face_x2_1: do i_comp_face = 1, size(periodic_face_data_x2(3,:))
            
            ! Get the face ID. This ID is given by the AFLR3 format
            jface = periodic_face_data_x2(3,i_comp_face)
            
            ! Partner face can only be a face with ifacetag = 9
            if (ifacetag(jface) == 11) then
                  
              ! Get the ID of the nodes which form the boundary face
              ivtmp_jface(:) = if2nq(:,jface)

              ! Get the coordinate of the nodes of the jface
              do i_vertex_jface = 1, nverticesperface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

              ! Set to zero the array that keeps track of the coordinates matches
              same_coord = 0

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs for 
              ! a shift.
              search_x2_1: do i_vertex_iface = 1, nverticesperface
                do i_vertex_jface = 1, nverticesperface
                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1
                    endif
                  end do
                  
                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common'
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x2_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x2_1
                          else
                            write(*,*) 'The periodic boundary conditions', &
                              & ' works only for two parallel boundary planes.', &
                              & ' Check x2 direction.'
                            
                            write(*,*) 'Exiting...'
                            
                            stop
                          end if
                        end do check_common_coord_x2_1
                      end do
                      
                      cycle search_x2_1
                    end if
                    cycle search_x2_1
                  end if
                  
                end do

                if (sum(same_coord(i_vertex_iface,:)) .lt. 2) then
                  ! Exit from the companion_face loop because none of the nodes
                  ! of the jface has two invariant coordinates
                  if (i_comp_face == size(periodic_face_data_x2(3,:))) then
                    write(*,*) 'No periodic partner face has been found for', &
                      & ' face', iface, 'of element', ielem, 'in the x2', &
                      & ' direction.'

                    write(*,*) 'Try to increase the diff_toll parameter and', &
                      & ' verify the input grid'
                    write(*,*) 'Exiting...'
                    stop
                  else
                    cycle partner_face_x2_1
                  end if
                endif

              end do search_x2_1
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x2(4,i_p_elem) = i_coord
                end if
              end do     

              ! Search the global ID of the element which owns the jface
              search_x2_2 : do j_p_elem = 1, size(periodic_face_data_x2(1,:))
                
                ! Get the global ID of the element
                jelem = periodic_face_data_x2(1,j_p_elem)

                ! Get the ID of the nodes which form the jelem
                iv_hexa_elem = ic2nh(:,jelem) 

                ! We have nverticesperface so we should find nverticesperface
                ! matches to claim that we have found the right (not self) 
                ! element
                match = 0
                do i_vertex_jface = 1, nverticesperface
                  do i_vertex_hexa = 1, nverticesperelem
                    if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
                      match = match + 1
                    endif
                  end do
                end do
                  
                ! Set in ef2e(2,iface,ielem) the ID of the adjoining element 
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

                  ! Set the ID of the adjoning face
                  do p_face = 1, nfacesperelem              

                    cnt = 0 ; ivtmp4(:) = 0
                    do k = 1, nverticesperface            
                      cnt = cnt + 1
                      k1  = eltypfaces(k,p_face)           
                      ivtmp4(cnt) = ic2nh(k1,jelem)         
                    enddo
                    ivtmp4 = isort(ivtmp4,nverticesperface)

                    ! Sort ivtmp_jface
                    ivtmp_jface = isort(ivtmp_jface,nverticesperface)

                    if(maxval(abs(ivtmp_jface(:)-ivtmp4(:))) == 0)  then
                      ef2e(1,iface,ielem) = p_face
                      exit
                    endif

                  end do

                else
                  cycle search_x2_2

                  endif

                  exit
                
                enddo search_x2_2
              
              exit
            
            endif

          end do partner_face_x2_1

        endif

      end do periodic_elem_x2_1

!      write(*,*) 'NOW THE SYMMETRIC PART'
      
      periodic_elem_x2_2 : do i_p_elem = 1, size(periodic_face_data_x2(1,:))

!        write(*,*) 'IN bc_elem loop'

        ! Get the global ID of the element that owns a boundary face
        ielem = periodic_face_data_x2(1,i_p_elem)

        ! Get the local (local for the element) ID of the boundary face
        iface = periodic_face_data_x2(2,i_p_elem)

        ! Check if the face is a "periodic" face
        if (abs(ef2e(1,iface,ielem)) == 11) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
          end do

          ! Search the partner face, i.e. the corresponding "periodic" face
          partner_face_x2_2: do i_comp_face = 1, size(periodic_face_data_x2(3,:))
            
            ! Get the face ID
            ! This ID is the one given by the AFLR3 format
            jface = periodic_face_data_x2(3,i_comp_face)
            
            ! Partner face can only be a face with ifacetag = 8
            if (ifacetag(jface) == 10) then
                  
              ! Get the ID of the nodes which form the boundary face
              ivtmp_jface(:) = if2nq(:,jface)

              ! Get the coordinate of the nodes of the jface
              do i_vertex_jface = 1, nverticesperface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

              ! Set to zero the array that keeps track of the coordinates
              ! matches
              same_coord = 0

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs for 
              ! a shift.
              search_x2_3: do i_vertex_iface = 1, nverticesperface
                do i_vertex_jface = 1, nverticesperface
                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1 
                    endif
                  end do

                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x2_2: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x2_2
                          else
                            write(*,*) 'The periodic boundary conditions', &
                              & ' works only for two parallel boundary planes.', &
                              & ' Check x2 direction.'
                            stop
                          end if
                        end do check_common_coord_x2_2
                      end do

                      cycle search_x2_3
                    endif
                    cycle search_x2_3
                  endif
                  
                end do

                if (sum(same_coord(i_vertex_iface,:)) .lt. 2) then
                  ! Exit from the companion_face loop because none of the nodes
                  ! of the jface has two invariant coordinates
                  if (i_comp_face == size(periodic_face_data_x2(3,:))) then
                    write(*,*) 'No periodic partner face has been found for', &
                      & ' face', iface, 'of element', ielem, 'in the x2', &
                      & ' direction.'
                    write(*,*) 'Try to increase the diff_toll parameter and', &
                      & ' verify the input grid'
                    write(*,*) 'Exiting...'
                    stop
                  else
                    cycle partner_face_x2_2
                  end if
                endif

              end do search_x2_3
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x2(4,i_p_elem) = i_coord
                end if
              end do

              ! Search the global ID of the element which owns the jface
              search_x2_4 : do j_p_elem = 1, size(periodic_face_data_x2(1,:))
                
                ! Get the global ID of the element
                jelem = periodic_face_data_x2(1,j_p_elem)

                ! Get the ID of the nodes which form the jelem
                iv_hexa_elem = ic2nh(:,jelem) 

                ! We have nverticesperface so we should find nverticesperface
                ! matches to claim that we have found the right (not self) 
                ! element
                match = 0
                do i_vertex_jface = 1, nverticesperface
                  do i_vertex_hexa = 1, nverticesperelem
                    if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
                      match = match + 1
                    endif
                  end do
                end do
                  
                ! Set in ef2e(2,iface,ielem) the ID of the adjoining element 
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

                  ! Set the ID of the adjoning face
                  do p_face = 1, nfacesperelem              

                    cnt = 0 ; ivtmp4(:) = 0
                    do k = 1, nverticesperface            
                      cnt = cnt + 1
                      k1  = eltypfaces(k,p_face)           
                      ivtmp4(cnt) = ic2nh(k1,jelem)         
                    enddo
                    ivtmp4 = isort(ivtmp4,nverticesperface)

                    ! Sort ivtmp_jface
                    ivtmp_jface = isort(ivtmp_jface,nverticesperface)

                    if(maxval(abs(ivtmp_jface(:)-ivtmp4(:))) == 0)  then
                      ef2e(1,iface,ielem) = p_face
                      exit
                    endif

                  end do

                else
                  cycle search_x2_4

                  endif

                  exit
                
                enddo search_x2_4
              

              exit
            
            endif

          end do partner_face_x2_2

        endif

      end do periodic_elem_x2_2


      ! ===================================================
      ! ===================================================
      ! Build ef2e for "periodic" faces in the x3 direction
      ! ===================================================
      ! ===================================================

      ! Loop over the element that owns a "periodic" face
      periodic_elem_x3_1 : do i_p_elem = 1, size(periodic_face_data_x3(1,:))

!        write(*,*) 'IN bc_elem loop'

        ! Get the global ID of the element that owns a periodic boundary face
        ielem = periodic_face_data_x3(1,i_p_elem)

        ! Get the local (local for the element) ID of the periodic boundary face
        iface = periodic_face_data_x3(2,i_p_elem)

        ! Check if the face is a "periodic" face with tag equals to 8
        if (abs(ef2e(1,iface,ielem)) == 12) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
          end do

          ! Search the partner face, i.e. the corresponding "periodic" face with
          ! face tag equals to 9
          partner_face_x3_1: do i_comp_face = 1, size(periodic_face_data_x3(3,:))
            
            ! Get the face ID. This ID is given by the AFLR3 format
            jface = periodic_face_data_x3(3,i_comp_face)
            
            ! Partner face can only be a face with ifacetag = 9
            if (ifacetag(jface) == 13) then
                  
              ! Get the ID of the nodes which form the boundary face
              ivtmp_jface(:) = if2nq(:,jface)

              ! Get the coordinate of the nodes of the jface
              do i_vertex_jface = 1, nverticesperface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

              ! Set to zero the array that keeps track of the coordinates
              ! matches
              same_coord = 0

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs for 
              ! a shift.
              search_x3_1: do i_vertex_iface = 1, nverticesperface
                do i_vertex_jface = 1, nverticesperface
                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1
                    endif
                  end do
                  
                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common'
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x3_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x3_1
                          else
                            write(*,*) 'The periodic boundary conditions', &
                              & ' works only for two parallel boundary planes.', &
                              & ' Check x3 direction.'
                            
                            write(*,*) 'Exiting...'
                            
                            stop
                          end if
                        end do check_common_coord_x3_1
                      end do
                      
                      cycle search_x3_1
                    end if
                    cycle search_x3_1
                  end if
                  
                end do

                if (sum(same_coord(i_vertex_iface,:)) .lt. 2) then
                  ! Exit from the companion_face loop because none of the nodes
                  ! of the jface has two invariant coordinates
                  if (i_comp_face == size(periodic_face_data_x3(3,:))) then
                    write(*,*) 'No periodic partner face has been found for', &
                      & ' face', iface, 'of element', ielem, 'in the x3', &
                      & ' direction.'

                    write(*,*) 'Try to increase the diff_toll parameter and', &
                      & ' verify the input grid'
                    write(*,*) 'Exiting...'
                    stop
                  else
                    cycle partner_face_x3_1
                  end if
                endif

              end do search_x3_1
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x3(4,i_p_elem) = i_coord
                end if
              end do     

              ! Search the global ID of the element which owns the jface
              search_x3_2 : do j_p_elem = 1, size(periodic_face_data_x3(1,:))
                
                ! Get the global ID of the element
                jelem = periodic_face_data_x3(1,j_p_elem)

                ! Get the ID of the nodes which form the jelem
                iv_hexa_elem = ic2nh(:,jelem) 

                ! We have nverticesperface so we should find nverticesperface
                ! matches to claim that we have found the right (not self) 
                ! element
                match = 0
                do i_vertex_jface = 1, nverticesperface
                  do i_vertex_hexa = 1, nverticesperelem
                    if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
                      match = match + 1
                    endif
                  end do
                end do
                  
                ! Set in ef2e(2,iface,ielem) the ID of the adjoining element 
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

                  ! Set the ID of the adjoning face
                  do p_face = 1, nfacesperelem              

                    cnt = 0 ; ivtmp4(:) = 0
                    do k = 1, nverticesperface            
                      cnt = cnt + 1
                      k1  = eltypfaces(k,p_face)           
                      ivtmp4(cnt) = ic2nh(k1,jelem)         
                    enddo
                    ivtmp4 = isort(ivtmp4,nverticesperface)

                    ! Sort ivtmp_jface
                    ivtmp_jface = isort(ivtmp_jface,nverticesperface)

                    if(maxval(abs(ivtmp_jface(:)-ivtmp4(:))) == 0)  then
                      ef2e(1,iface,ielem) = p_face
                      exit
                    endif

                  end do

                else
                  cycle search_x3_2

                  endif

                  exit
                
                enddo search_x3_2
              
              exit
            
            endif

          end do partner_face_x3_1

        endif

      end do periodic_elem_x3_1
      
      periodic_elem_x3_2 : do i_p_elem = 1, size(periodic_face_data_x3(1,:))

        ! Get the global ID of the element that owns a boundary face
        ielem = periodic_face_data_x3(1,i_p_elem)

        ! Get the local (local for the element) ID of the boundary face
        iface = periodic_face_data_x3(2,i_p_elem)

        ! Check if the face is a "periodic" face
        if (abs(ef2e(1,iface,ielem)) == 13) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
          end do

          ! Search the partner face, i.e. the corresponding "periodic" face
          partner_face_x3_2: do i_comp_face = 1, size(periodic_face_data_x3(3,:))
            
            ! Get the face ID
            ! This ID is the one given by the AFLR3 format
            jface = periodic_face_data_x3(3,i_comp_face)
            
            ! Partner face can only be a face with ifacetag = 8
            if (ifacetag(jface) == 12) then
                  
              ! Get the ID of the nodes which form the boundary face
              ivtmp_jface(:) = if2nq(:,jface)

              ! Get the coordinate of the nodes of the jface
              do i_vertex_jface = 1, nverticesperface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

              ! Set to zero the array that keeps track of the coordinates matches
              same_coord = 0

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs for 
              ! a shift.
              search_x3_3: do i_vertex_iface = 1, nverticesperface
                do i_vertex_jface = 1, nverticesperface
                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1 
                    endif
                  end do

                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x3_2: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x3_2
                          else
                            write(*,*) 'The periodic boundary conditions', &
                              & ' works only for two parallel boundary planes.', &
                              & ' Check x3 direction.'
                            stop
                          end if
                        end do check_common_coord_x3_2
                      end do

                      cycle search_x3_3
                    endif
                    cycle search_x3_3
                  endif
                  
                end do

                if (sum(same_coord(i_vertex_iface,:)) .lt. 2) then
                  ! Exit from the companion_face loop because none of the nodes
                  ! of the jface has two invariant coordinates
                  if (i_comp_face == size(periodic_face_data_x3(3,:))) then
                    write(*,*) 'No periodic partner face has been found for', &
                      & ' face', iface, 'of element', ielem, 'in the x3', &
                      & ' direction.'
                    write(*,*) 'Try to increase the diff_toll parameter and', &
                      & ' verify the input grid'
                    write(*,*) 'Exiting...'
                    stop
                  else
                    cycle partner_face_x3_2
                  end if
                endif

              end do search_x3_3
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x3(4,i_p_elem) = i_coord
                end if
              end do

              ! Search the global ID of the element which owns the jface
              search_x3_4 : do j_p_elem = 1, size(periodic_face_data_x3(1,:))
                
                ! Get the global ID of the element
                jelem = periodic_face_data_x3(1,j_p_elem)

                ! Get the ID of the nodes which form the jelem
                iv_hexa_elem = ic2nh(:,jelem) 

                ! We have nverticesperface so we should find nverticesperface
                ! matches to claim that we have found the right (not self) 
                ! element
                match = 0
                do i_vertex_jface = 1, nverticesperface
                  do i_vertex_hexa = 1, nverticesperelem
                    if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
                      match = match + 1
                    endif
                  end do
                end do
                  
                ! Set in ef2e(2,iface,ielem) the ID of the adjoining element 
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

                  ! Set the ID of the adjoning face
                  do p_face = 1, nfacesperelem              

                    cnt = 0 ; ivtmp4(:) = 0
                    do k = 1, nverticesperface            
                      cnt = cnt + 1
                      k1  = eltypfaces(k,p_face)           
                      ivtmp4(cnt) = ic2nh(k1,jelem)         
                    enddo
                    ivtmp4 = isort(ivtmp4,nverticesperface)

                    ! Sort ivtmp_jface
                    ivtmp_jface = isort(ivtmp_jface,nverticesperface)

                    if(maxval(abs(ivtmp_jface(:)-ivtmp4(:))) == 0)  then
                      ef2e(1,iface,ielem) = p_face
                      exit
                    endif

                  end do

                else
                  cycle search_x3_4

                  endif

                  exit
                
                enddo search_x3_4
              

              exit
            
            endif

          end do partner_face_x3_2

        endif

      end do periodic_elem_x3_2


! =======================================================
! Sanity check!  Check if the faces are connected in a symmetric fashion
! =======================================================

      if(testing) then

        ! x1 direction
        do i_p_elem = 1, size(periodic_face_data_x1(1,:))

          ! Get the global ID of the element that owns a boundary face
          ielem = periodic_face_data_x1(1,i_p_elem)

          ! Get the local (local for the element) ID of the boundary face
          iface = abs(periodic_face_data_x1(2,i_p_elem))

          jelem_iface = ef2e(2,iface,ielem)

          jface_iface = ef2e(1,iface,ielem)

          ielem_jface = ef2e(2,jface_iface,jelem_iface)

          iface_jface = ef2e(1,jface_iface,jelem_iface)

          if (iface /= iface_jface .or. ielem /= ielem_jface) then
            write(*,*) 'Periodic faces are not connected in a symmetric fashion'
            write(*,*) 'Check x1 direction'
            write(*,*) 'Exiting...'
            stop
          end if

        end do

        ! x2 direction
        do i_p_elem = 1, size(periodic_face_data_x2(1,:))

          ! Get the global ID of the element that owns a boundary face
          ielem = periodic_face_data_x2(1,i_p_elem)

          ! Get the local (local for the element) ID of the boundary face
          iface = abs(periodic_face_data_x2(2,i_p_elem))

          jelem_iface = ef2e(2,iface,ielem)

          jface_iface = ef2e(1,iface,ielem)

          ielem_jface = ef2e(2,jface_iface,jelem_iface)

          iface_jface = ef2e(1,jface_iface,jelem_iface)

          if (iface /= iface_jface .or. ielem /= ielem_jface) then
            write(*,*) 'Periodic faces are not connected in a symmetric fashion'
            write(*,*) 'Check x2 direction'
            write(*,*) 'Exiting...'
            stop
          end if

        end do

        ! x3 direction
        do i_p_elem = 1, size(periodic_face_data_x3(1,:))

          ! Get the global ID of the element that owns a boundary face
          ielem = periodic_face_data_x3(1,i_p_elem)

          ! Get the local (local for the element) ID of the boundary face
          iface = abs(periodic_face_data_x3(2,i_p_elem))

          jelem_iface = ef2e(2,iface,ielem)

          jface_iface = ef2e(1,iface,ielem)

          ielem_jface = ef2e(2,jface_iface,jelem_iface)

          iface_jface = ef2e(1,jface_iface,jelem_iface)

          if (iface /= iface_jface .or. ielem /= ielem_jface) then
            write(*,*) 'Periodic faces are not connected in a symmetric fashion'
            write(*,*) 'Check x3 direction'
            write(*,*) 'Exiting...'
            stop
          end if

        end do

      endif

      ! Construct iae2v and jae2v
      ! =========================
      nnze2v = nverticesperelem*nelems
      allocate(iae2v    (nelems+1))  ; iae2v     = 0
      allocate(iae2v_tmp(nelems+1))  ; iae2v_tmp = 0
      allocate(jae2v    (nnze2v  ))  ; jae2v     = 0
      allocate(jae2v_tmp(nnze2v  ))  ; jae2v_tmp = 0

      iae2v(1) = 1
      do j = 2,nelems+1
        iae2v(j) = iae2v(j-1) + nverticesperelem
      end do
      
      iae2v_tmp = iae2v

      cnt = 0
      do j = 1,nelems
        do i = 1,nverticesperelem
          cnt = cnt + 1
          jae2v(cnt) = ic2nh(i,j)
        enddo
      enddo

      jae2v_tmp = jae2v

      ! ====================================================================================================
      ! Construct iae2e and jae2e;  arrays that track connectivities between elements
      ! ====================================================================================================

      nnze2e = 0                                           !  count the size of all e2e connections (exclude BC's)
      do ielem = 1,nelems                                  !  element loop
        do iface = 1,nfacesperelem                         !  face loop
          if(ef2e(1,iface,ielem) > 0) nnze2e = nnze2e + 1  !  only add element connections to other elements
        enddo
      enddo

      allocate(iae2e(nelems+1)) ; iae2e = 0                ! allocate arrays based on initial counting
      allocate(jae2e(nnze2e))   ; jae2e = 0

      cnt = 1
      do ielem = 1, nelems
        iae2e(ielem) = cnt
        do iface = 1, nfacesperelem
          if ( ef2e(1,iface,ielem) > 0 ) then                  ! only add element connections to other elements
            jae2e(cnt) = ef2e(2,iface,ielem)
            cnt = cnt+1
          end if
        end do
      end do
      iae2e(nelems+1) = cnt

      ! ====================================================================================================
      ! Construct iae2e2e and jae2e2e;  arrays that track secondary-connectivities between elements (nearest-nearest neighbors)
      ! ====================================================================================================

      allocate(iae2e2e    (nelems+1))    ; iae2e2e     = 0
      allocate(jae2e2e_tmp(6*nnze2e))    ; jae2e2e_tmp = 0      !  tmp vector allocated to 6 times the size of nearest neighbor array
      allocate(  stack(100) ) ;   stack = 0
      allocate(ordered(100) ) ; ordered = 0

      nnze2e2e = 1
      do ielem = 1, nelems                                      !  Loop over all elements
         iae2e2e(ielem) = nnze2e2e
         cnt = 0
         do jelem = iae2e(ielem),iae2e(ielem+1) - 1             !  Loop over connected elements
            cnt = cnt + 1                                       !  update counter
            stack(cnt) = jae2e(jelem)                           !  add the connected element to the stack  
            do kelem = iae2e(jae2e(jelem)),iae2e(jae2e(jelem)+1) - 1       !  Loop over all connected elements to the connected element
              cnt = cnt + 1                                     !  update counter
              stack(cnt) = jae2e(kelem)                         !  add nearest-nearest neighbor to the stack
            enddo
         enddo
         call qsorti(stack,ordered,cnt)
         jae2e2e_tmp(nnze2e2e) = stack(ordered(1))
         nnze2e2e = nnze2e2e + 1
         do i = 2,cnt
           if(stack(ordered(i))-stack(ordered(i-1)) /= 0 ) then
              jae2e2e_tmp(nnze2e2e) = stack(ordered(i))
              nnze2e2e = nnze2e2e + 1
           endif
         enddo
      enddo
      iae2e2e(nelems+1) = nnze2e2e

      allocate(jae2e2e(nnze2e2e)) ; jae2e2e(1:nnze2e2e) = jae2e2e_tmp(1:nnze2e2e)

      deallocate(jae2e2e_tmp)
      deallocate(stack,ordered)

      ! Deallocate memory
      deallocate(iav2e,jav2e)
      deallocate(ivtmp1,ivtmp2,ivtmp3,ivtmp4)
      deallocate(test_conL,test_conR,test_cBCL,test_cBCR)

      deallocate(ivtmp_iface)
      deallocate(ivtmp_jface)

      deallocate(vx_iface)
      deallocate(vx_jface)
      deallocate(same_coord)
      deallocate(iv_hexa_elem)

    else
      write(*,*) 'Unknown number of spatial dimension of the problem:', ndim
      write(*,*) 'Exiting...'
      stop
    end if ! End if ndim == 3

  end subroutine e2e_connectivity_aflr3     !  SERIAL Routine

  !============================================================================
  
    subroutine set_element_orders_Serial()     !  Serial Routine

    ! Load modules
    use variables
    use collocationvariables
    use controlvariables,   only : p_non_conforming, p_refine_strategy
    use referencevariables, only : npoly, npoly_max, nfacesperelem
    use precision_vars, only     : magnitude
    use mpimod

    ! Nothing is implicitly defined
    implicit none
   
    integer,  parameter  :: seed = 86441
    real(wp), parameter  :: tol = 1.0e-9_wp

    integer,  dimension(:), allocatable :: irand_order, jrand_order

    integer :: ielem, nhex, iface
    integer :: i, j, nval, icnt, ierr, nx, nxb2

    real(wp)                                           :: r

    call srand(seed)

    npoly_max = -1000

    nhex = size(ic2nh,2)

    if(allocated(elem_props)) deallocate(elem_props) ; allocate(elem_props(2,1:nhex))
    elem_props(:,:) = -1000

    elem_props(1,:) = 1
    elem_props(2,:) = npoly+1

    !-- set ef2e(9,:,:) = 0 so that all faces are h-conforming
    ef2e(9,:,:) = 0

    !  adjust the element polynomials as per directives
    if(p_non_conforming .eqv. .true.) then

      select case(p_refine_strategy)

        case(1)   !  Quadrants set to different orders

          do ielem = 1,nhex
  
            ! NW quad  x <= +tol ; y >= -tol
            icnt = 0
            do j=1,8
              if((vx_master(1,ic2nh(j,ielem)) <= +tol) .and. (vx_master(2,ic2nh(j,ielem)) >= -tol)) icnt = icnt + 1
            enddo
            if(icnt == 8) elem_props(2,ielem) = npoly+2
            ! SE quad  x >= -tol ; y <= +tol
            icnt = 0
            do j=1,8
              if((vx_master(1,ic2nh(j,ielem)) >= -tol) .and. (vx_master(2,ic2nh(j,ielem)) <= +tol)) icnt = icnt + 1
            enddo
            if(icnt == 8) elem_props(2,ielem) = npoly+2
            ! NE quad  x, y >= -tol
            icnt = 0
            do j=1,8
              if((vx_master(1,ic2nh(j,ielem)) >= -tol) .and. (vx_master(2,ic2nh(j,ielem)) >= -tol)) icnt = icnt + 1
            enddo
            if(icnt == 8) elem_props(2,ielem) = npoly+3
   
          enddo

        case(2)   !  Upper half of elements elevated in order  :  default test

          do ielem = nhex/2+1,nhex
            elem_props(2,ielem) = npoly+2
          enddo

        case(3)   !  every fourth one elevated

          do ielem = 1,nhex,4
            elem_props(2,ielem) = npoly+2
          enddo

        case(4)  !  randomly elevate polynomials

          nval = nhex / 10
          if(nval == 0) nval = 1

          if(allocated(irand_order)) deallocate(irand_order) ; allocate(irand_order(1:nval))
          if(allocated(jrand_order)) deallocate(jrand_order) ; allocate(jrand_order(1:nval))

          do i = 1,nval
            r = rand()
            irand_order(i) = 1 + floor( r * nhex) 
          enddo

          jrand_order = isort(irand_order,nval)

          icnt = 1 ; irand_order(1) = jrand_order(1) ;

          do i = 2,nval
            if(jrand_order(i) /= jrand_order(i-1)) then
              icnt = icnt + 1
              irand_order(icnt) = jrand_order(i)
            endif
          enddo

          do ielem = 1,icnt
            elem_props(2,irand_order(ielem)) = npoly+2
          enddo
      case(5)
        elem_props(2,1) = npoly+1
        elem_props(2,2) = npoly+2
        elem_props(2,3) = npoly+2
      case(6)
        elem_props(2,1) = npoly+2
        elem_props(2,2) = npoly+1
        elem_props(2,3) = npoly+2
      case(7)
        elem_props(2,1) = npoly+1
        elem_props(2,2) = npoly+1
        elem_props(2,3) = npoly+2
      case(8)
        elem_props(2,1) = npoly+2
        elem_props(2,2) = npoly+1
        elem_props(2,3) = npoly+1
        elem_props(2,4) = npoly+1
        elem_props(2,5) = npoly+1
        elem_props(2,6) = npoly+1
        elem_props(2,7) = npoly+1
        elem_props(2,8) = npoly+1
      case(9)
        elem_props(2,1) = npoly+1
        elem_props(2,2) = npoly+2
        elem_props(2,3) = npoly+1
        elem_props(2,4) = npoly+1
        elem_props(2,5) = npoly+1
        elem_props(2,6) = npoly+1
        elem_props(2,7) = npoly+1
        elem_props(2,8) = npoly+1
      case(10)
        do ielem = 1,3**3
          elem_props(2,ielem) = npoly+1
        enddo
        elem_props(2,13) = npoly+2
        elem_props(2,14) = npoly+2
        elem_props(2,15) = npoly+2
        elem_props(2,11) = npoly+2
        elem_props(2,17) = npoly+2
        elem_props(2,5) = npoly+2
        elem_props(2,23) = npoly+2

        !elem_props(2,2) = npoly+2
        !elem_props(2,11) = npoly+2
      case(11)
        do ielem = 1,3**3
          elem_props(2,ielem) = npoly+2
        enddo
        elem_props(2,2) = npoly+1
        !elem_props(2,23) = npoly+2
      case(12)
        do ielem = 1,4**3
          elem_props(2,ielem) = npoly+1
        enddo
        elem_props(2,2) = npoly+1
      case(13)
        do ielem = 1,5**3
          elem_props(2,ielem) = npoly+1
        enddo
        elem_props(2,1) = npoly+2
        elem_props(2,3) = npoly+2
        elem_props(2,15) = npoly+2

      case(14)
        !-- sphere test case
        !-- set the surface elements to p+1
        do ielem = 1,1
          elem_props(2,ielem) = npoly+2
        enddo
        do ielem = 2,2
          elem_props(2,ielem) = npoly+1
        enddo
      case(15)
        !-- sphere test case
        !-- set the surface elements to p+1
        do ielem = 1,1
          elem_props(2,ielem) = npoly+2
        enddo
        do ielem = 2,3
          elem_props(2,ielem) = npoly+1
        enddo
      case(16)
        !-- sphere with origin at (0,0,0) and a radius = sqrt(3)
!-- HACK
        do ielem = 1,nhex
          elem_props(2,ielem) = npoly+1
        enddo
        do ielem = 2,156,2
          elem_props(2,ielem) = npoly+2
        enddo
        !-- THIS IS WHAT IT SHOULD BE
        !do ielem =1,nhex
        !  if(abs(sqrt(vx(1,e2v(1,ielem))**2+vx(2,e2v(1,ielem))**2+vx(3,e2v(1,ielem))**2)&
        !     -radius)<1.0e-12_wp)then
        !    elem_props(2,ielem) = npoly+2
        !   else 
        !     elem_props(2,ielem) = npoly+1
        !   endif
        !enddo
      case(17)
        !-- sphere with origin at (0,0,0) and a radius = sqrt(3)
        do ielem = 1,nhex
          elem_props(2,ielem) = npoly+1
          do iface = 1,nfacesperelem
             if(ef2e(1,iface,ielem) == -16) elem_props(2,ielem) = npoly + 2
          enddo
        enddo
      case(18)
        !-- set the surface elements to p+1
        do ielem = 1,1
          elem_props(2,ielem) = npoly+2
        enddo
        do ielem = 2,8
          elem_props(2,ielem) = npoly+1
        enddo
      case(19)
        !-- SnowFlake test case       
        elem_props(2,1) = npoly+2

        do ielem = 2,24
          elem_props(2,ielem) = npoly+1
        enddo

      case(20)

        !-- Plus sign on Cartesian mesh  (assumes equal dimensions in X and Y directions)

        !         |-------------|
        !         |      +      |
        !         |      +      |
        !         |      +      |
        !         |+ + + + + + +|
        !         |      +      |
        !         |      +      |
        !         |      +      |
        !         |-------------|

        nx = int(sqrt(nhex + 1.0e-10_wp))

        if(nx*nx - nhex /= 0) then
          write(*,*)'inappropriate mesh for elevating polynomial orders.  Not a square mesh '
          write(*,*)'stopping'
          call PetscFinalize(ierr) ; stop
        endif
        nxb2 = nx/2 + 1

        elem_props(2,:) = npoly+1

        ielem = 0 
        do i = 1,nx
          do j = 1,nx
            ielem = ielem + 1
            if((i == nxb2) .or. (j == nxb2)) elem_props(2,ielem) = npoly+2
          enddo
        enddo

      case(21)

        !  elevate two planar surfaces on the Taylor-Green Cube test cube
        !             ______________
        !            /+            /|
        !          / |+          / +|
        !         |-------------|  +|
        !         |+ |+        +|  +|
        !         |+ |+        +|  +|
        !         |+ |+        +|  +|
        !         |+ |+        +|  +|
        !         |+ |+        +|  +|
        !         |+ /+        +|  /
        !         |/           +|/
        !         |-------------|

        elem_props(2,:) = npoly+1

        nx = int((nhex + 1.0e-10_wp)**(1.0_wp/3.0_wp))

        if(nx*nx*nx - nhex /= 0) then
          write(*,*)'inappropriate mesh for elevating polynomial orders.  Not a cube mesh '
          write(*,*)'stopping'
          call PetscFinalize(ierr) ; stop
        endif

        icnt = 0
        do i = 1,nx
          do j = 1,nx
            icnt = icnt + 1
            elem_props(2,       icnt) = npoly+2
            elem_props(2,nhex+1-icnt) = npoly+2
          enddo
        enddo
      end select

    endif

    !  Serial ordering : ! ef2e(:,iface,ielem) 
    !                       :  (1,j,k) = Adjoining element face ID
    !                       :  (2,j,k) = Adjoining element ID
    !                       :  (3,j,k) = Adjoining element process ID
    !                       :  (4,j,k) = Adjoining element polynomial order
    !                       :  (5,j,k) = Number of Adjoining elements
    !                       :  (6,j,k) = HACK self polynomial order assigned to each face
    !                       :  (7,j,k) = face orientation of adjoining face (0,1,2,3)

    do ielem = 1,nhex

      ef2e(6,:,ielem) = elem_props(2,ielem)

      if(elem_props(2,ielem) >= npoly_max + 1) npoly_max = elem_props(2,ielem) - 1

      do iface = 1,nfacesperelem

        if(ef2e(2,iface,ielem) > 0) then
          ef2e(4,iface,ielem) = elem_props(2,ef2e(2,iface,ielem))
        else
          ef2e(4,iface,ielem) = elem_props(2,ielem)
        endif

      enddo

    enddo

    end subroutine set_element_orders_Serial    !  Serial Routine

  !============================================================================
  
    subroutine set_element_orders()     !  Parallel Routine

    ! Load modules
    use variables
    use collocationvariables
    use referencevariables, only : ihelems, nfacesperelem, myprocid
    use controlvariables, only : hrefine

    ! Nothing is implicitly defined
    implicit none
   
    integer :: ielem, iface, qface, ierr, elem_min, elem_max

    ! determine the lowest and highest element process is connected to (including self)

    elem_min = ihelems(1) ; elem_max = ihelems(2) ;

    do ielem = ihelems(1),ihelems(2)

       do iface = 1,nfacesperelem
         elem_min = min(ef2e(2,iface,ielem),elem_min)
         elem_max = max(ef2e(2,iface,ielem),elem_max)
       enddo

    enddo

    ! allocate size of elem_props for each process.  Includes self and all elements connected to faces
    if(allocated(elem_props)) deallocate(elem_props) ; allocate(elem_props(2,elem_min:elem_max))
    elem_props(1,:) = 1 ; elem_props(2,:) = -1000

    !  Parallel ordering : ! ef2e(:,iface,ielem) 
    !                         :  (1,j,k) = Adjoining element face ID
    !                         :  (2,j,k) = Adjoining element ID
    !                         :  (3,j,k) = Adjoining element process ID
    !                         :  (4,j,k) = Adjoining element polynomial order
    !                         :  (5,j,k) = Number of Adjoining elements
    !                         :  (6,j,k) = HACK self polynomial order assigned to each face

    ! ef2e has poly order stored in the parallel ordering 
    do ielem = ihelems(1),ihelems(2)

      !-HACK
      if(hrefine)then
        !-- not testing
      else
        qface = size(ef2e,2)
        if(sum(ef2e(6,:,ielem))/qface /= ef2e(6,1,ielem)) then
           write(*,*)'mpi bug in transfering ef2e:   Stopping'
           call PetscFinalize(ierr) ; stop ! Finalize MPI and the hooks to PETSc
        endif
      endif
      elem_props(2,ielem) = ef2e(6,1,ielem) 

      do iface = 1,nfacesperelem
         if( (elem_props(2,ef2e(2,iface,ielem)) /= -1000) .and.  &
             (elem_props(2,ef2e(2,iface,ielem)) /= ef2e(4,iface,ielem)) ) then
             write(*,*)'something wrong in set_element_orders: parallel.  Stopping'
             call PetscFinalize(ierr) ; stop ; ! Finalize MPI and the hooks to PETSc
         else
           elem_props(2,ef2e(2,iface,ielem)) = ef2e(4,iface,ielem)
         endif
      enddo
! write(*,*)"myprocid = ",myprocid,"from within set ielem = ",ielem,"elem_props(2,ielem) = ",elem_props(2,ielem),&
!           "ef2e(4,iface,ielem) = ",ef2e(4,nfacesperelem,ielem)
    enddo

    end subroutine set_element_orders    !  Parallel Routine

  !============================================================================

  subroutine create_ldg_flip_flop_sign()    !   SERIAL Routine

!   SERIAL ROUTINE
!     Creates flip flow operators on faces of elements.
!   SERIAL ROUTINE
    
    ! Load modules
    use collocationvariables
    use referencevariables
    use variables, only: ef2e 

    ! Nothing is implicitly defined
    implicit none
   
    integer, dimension(nfacesperelem) :: face_pairs

    logical :: hit_boundary, close_loop
    logical, parameter :: testing = .true.

    integer :: i, i_elem, i_face, face_id, elem_id, opposite_i_face, flip, &
             & partner_elem, partner_face, opposite_partner_face, init, &
             & face_id_bk, elem_id_bk, cnt

    ! Set the face pairs
    if (ndim == 2) then
      write(*,*)'ndim = 2 not supported', i_face
      write(*,*)'Check the subroutine create_ldg_flip_flop_sign'
      write(*,*)'Stopping...'
      stop                
    else if (ndim == 3) then
      face_pairs(1) = 6
      face_pairs(2) = 4
      face_pairs(3) = 5
      face_pairs(4) = 2
      face_pairs(5) = 3
      face_pairs(6) = 1
    end if

    init = 34

    ! Allocate array for the LDG flip-flop
    allocate(ldg_flip_flop_sign(nfacesperelem,nelems))
    ldg_flip_flop_sign(:,:) = init

    ! Initialize logical variables
    hit_boundary = .false.
    close_loop = .false.

    ! Loop through all the elements
    do i_elem = 1, nelems
       
      do i = 1, nfacesperelem
        
        hit_boundary = .false.
        close_loop = .false.
        
        i_face = i
        opposite_i_face = face_pairs(i_face)
        ! Check if the i_face has been already set
        if (ldg_flip_flop_sign(i_face,i_elem) /= init) then
          ! Check if the opposite face has been set
          if (ldg_flip_flop_sign(opposite_i_face,i_elem) /= init) then
            ldg_flip_flop_sign(i_face,i_elem) = -ldg_flip_flop_sign(opposite_i_face,i_elem)
          ! i_face and its opposite face have not been already set
          else
            ldg_flip_flop_sign(i_face,i_elem) = +1
          end if
        end if

        ! Check if the i_face is a boundary face 
        if (ef2e(1,i_face,i_elem) < 0) then
          ! We cannot do anything because face_id is a boundary face
          hit_boundary = .true.
        end if

        ! Prepare IDs for the do while loop 
        ! We want to move along a "connection line"
        face_id = i_face
        elem_id = i_elem
        
        do while ((hit_boundary .eqv. .false.) .and. (close_loop .eqv. .false.))

          flip = -1
          partner_elem = ef2e(2,face_id,elem_id)
          partner_face = ef2e(1,face_id,elem_id)
          opposite_partner_face = face_pairs(partner_face)

          ldg_flip_flop_sign(partner_face,partner_elem) = ldg_flip_flop_sign(face_id,elem_id)*flip
          ldg_flip_flop_sign(opposite_partner_face,partner_elem) = ldg_flip_flop_sign(face_id,elem_id)

          face_id_bk = face_id
          elem_id_bk = elem_id

          if (ef2e(1,opposite_partner_face,partner_elem) < 0) then
            hit_boundary = .true.
            ldg_flip_flop_sign(opposite_partner_face,partner_elem) = 0
          else
            face_id = opposite_partner_face
            elem_id = partner_elem
          end if

          if (elem_id == i_elem) then
            close_loop = .true.
          end if
        
        end do
      
      end do

    end do

    ! Check if everything sum up to zero
    cnt = 0
    
    if (testing) then
      ! Loop through all the elements
      do i_elem = 1, nelems
       
        do i_face = 1, nfacesperelem

          if (ef2e(1,i_face,i_elem) < 0) then
            cycle
          else
            partner_elem = ef2e(2,i_face,i_elem)
            partner_face = ef2e(1,i_face,i_elem)
            cnt = cnt + abs(ldg_flip_flop_sign(i_face,i_elem) + ldg_flip_flop_sign(partner_face,partner_elem))
          end if

        end do

      end do

      if (cnt /=0) then
        write(*,*) 'The sum of the flip-flop is not zero. Value:', cnt
        write(*,*) ' Check subroutine create_ldg_flip_flop_sign'
        write(*,*) 'Stopping...'
        stop
      end if

    end if

    ! Set to zero the LDG flip-flop for boundary faces
    ! ================================================
    ! Loop through all the elements
    do i_elem = 1, nelems
       
      do i_face = 1, nfacesperelem

        if (ef2e(1,i_face,i_elem) < 0) then
          ldg_flip_flop_sign(i_face,i_elem) = 0
        end if
      end do
    end do

  end subroutine create_ldg_flip_flop_sign

  !============================================================================

  pure function face_pairs(dir)
   
    ! Load modules

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: dir
    integer, dimension(2) :: face_pairs

    continue

    if (dir == 1) then
      face_pairs(1) = 5
      face_pairs(2) = 3
    else if (dir == 2) then
      face_pairs(1) = 2
      face_pairs(2) = 4
    else
      face_pairs(1) = 1
      face_pairs(2) = 6
    end if

  end function face_pairs

  !============================================================================

! SERIAL Routine

  pure subroutine data_partner_element_serial(node,face_id,elem_id,k_node,k_elem,k_face,i_node)

    ! Load modules
    use referencevariables, only : nodesperface
    use variables, only: efn2efn, ef2e, ifacenodes

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: node, face_id, elem_id
    integer, intent(out) :: k_node, k_elem, k_face, i_node

    integer :: j_node

    continue

    ! Index in facial ordering
    j_node = nodesperface*(face_id-1) + node

    ! Volumetric node index corresponding to facial node index
    i_node = ifacenodes(j_node)

    ! Volumetric index of partner node
    k_node = efn2efn(1,j_node,elem_id)

    ! Element index of partner node
    k_elem = efn2efn(2,j_node,elem_id)

    ! Face index of the partner face
    k_face = ef2e(1,face_id,elem_id)

  end subroutine data_partner_element_serial     !   SERIAL Routine

  !============================================================================

! pure function WENO_Adjoining_Data(k_node,k_face)     !   PARALLEL Routine
  function WENO_Adjoining_Data(k_node,k_face)
     !  Grab the data that lives at the first point off the surface of the
     !  adjoining element

     use referencevariables

     integer,   intent(in) :: k_node,k_face
     integer               :: WENO_Adjoining_Data

         if(k_face == 1) then
         WENO_Adjoining_Data = k_node + nodesperface
     elseif(k_face == 2) then
         WENO_Adjoining_Data = k_node + nodesperedge
     elseif(k_face == 3) then
         WENO_Adjoining_Data = k_node - 1
     elseif(k_face == 4) then
         WENO_Adjoining_Data = k_node - nodesperedge
     elseif(k_face == 5) then
         WENO_Adjoining_Data = k_node + 1
     elseif(k_face == 6) then
         WENO_Adjoining_Data = k_node - nodesperface
     else 
         write(*,*)'invalid face. Stopping'
         stop
     endif

  end function WENO_Adjoining_Data     !   PARALLEL Routine

  !============================================================================

  pure function Pencil_Coord(Ns,jdir,iface,i)     !   PARALLEL Routine

    integer, intent(in)  :: Ns,jdir,iface,i

    integer              :: Pencil_Coord

    Pencil_Coord = 0

    select case (jdir)
      case(1)
        Pencil_Coord =  1 + (mod((iface-0),Ns)-1)*Ns + int((iface-0)/Ns) * Ns*Ns + (i-1)
      case(2)
        Pencil_Coord =  1 + (mod((iface-1),Ns)-0)    + int((iface-1)/Ns) * Ns*Ns + (i-1) * Ns
      case(3)
        Pencil_Coord =            iface + (i-1) * Ns*Ns
    end select

  end function     !   PARALLEL Routine

  !============================================================================

  subroutine WENO_Intrp_Face_Nodes()   !   PARALLEL  Routine

    !  Preprocessing step to build physical locations needed for interpolation

    use referencevariables
    use interpolation
    use variables
    use SSWENOvariables
    use mpimod
    use petscvariables, only: xpetscWENO_partner, xlocpetscWENO_partner

    implicit none
    ! indices
    integer :: ielem, inode, jnode, jdir, ipen, face_id, iface, kface
    integer :: gnode, iloc, kelem, knode, nodespershell
    integer :: i, i_err
    integer :: extrnal_xi_cnt, extrnal_xi_sum
    integer, dimension(2)  :: faceLR

    real(wp), parameter                 :: sqrt5  = sqrt(5.0_wp)
    real(wp), parameter                 :: tol    = 1.0e-09_wp

    real(wp), dimension(3,nodesperedge) :: xgS
    real(wp), dimension(3)              :: xgL, xgR

    real(wp), dimension(3,8)            :: comp2phys_coeffs
    real(wp)                            :: tmp, rnorm, rnorm_max
    real(wp)                            :: XI_max, XI_max_Glob
    real(wp)                            :: del

!   allocate WENO face node arrays
    nodespershell = nfacesperelem*nodesperface
    allocate(xgWENO_self   (3,nodespershell,ihelems(1):ihelems(2))) ; xgWENO_self    = -10000.0_wp ;
    allocate(XIWENO_partner(3,nodespershell,ihelems(1):ihelems(2))) ; XIWENO_partner = -10000.0_wp ;

    ! Build the local interpolated physical value for each face element.

    ! loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)

      do jdir = 1,ndim            ! Directional loop

        faceLR = face_pairs(jdir)

        do ipen = 1,nodesperface

          !  Grab a pencil of data
          do i = 1,nodesperedge
            inode    = Pencil_Coord(nodesperedge,jdir,ipen,i)
            xgS(:,i) =   xg(:,inode,ielem)    !  Neqns
          enddo

          !  Extrapolate pencil coordinates to outside the element
          del = (1.0_wp - 1.0_wp/sqrt5) * WENO_Extrp
          xgL     = Poly3_Intrp(3,xgS,-1.0_wp,+1.0_wp,-1.0_wp-del)
          xgR     = Poly3_Intrp(3,xgS,-1.0_wp,+1.0_wp,+1.0_wp+del)

          face_id = faceLR(1) ; jnode = nodesperface*(face_id-1) + ipen ;
          xgWENO_self(:,jnode,ielem) = xgL

          face_id = faceLR(2) ; jnode = nodesperface*(face_id-1) + ipen ;
          xgWENO_self(:,jnode,ielem) = xgR

        enddo

      enddo

    enddo


    call UpdateComm1DGhostDataWENOGeom(xgWENO_self, xghstWENO_partner,             &
                                       xpetscWENO_partner, xlocpetscWENO_partner,  &
                                ndim, size(xgWENO_self,2), size(xghstWENO_partner,2))

     iloc = 0 
     do ielem = ihelems(1), ihelems(2)                                   ! element loop

       do iface = 1,nfacesperelem                                        ! Face Loop

         if (ef2e(3,iface,ielem) /= myprocid) then                       ! off process

           do ipen = 1, nodesperface                                     ! loop over facenodes
             iloc  = iloc + 1
             jnode = nodesperface*(iface-1) + ipen                       ! Shell node index
             inode = ifacenodes(jnode)                                   ! Volumetric node index corresponding to facial node index
             gnode = efn2efn(3,jnode,ielem)                              ! This is pointing to ghost stack not volumetric stack
             xgWENO_partner(:,jnode,ielem) = xghstWENO_partner(:,iloc) & ! off process partner data
                   - xghst_LGL(:,iloc) + xg(:,inode,ielem)               ! account for possibility of non-periodic domain.
           end do

         else                                                            ! On process

           do ipen = 1, nodesperface
             jnode = nodesperface*(iface-1)+ipen                         ! Index in facial ordering
             inode = ifacenodes(jnode)                                   ! Volumetric node index corresponding to facial node index
             if(ef2e(1,iface,ielem) < 0)then
               xgWENO_partner(:,jnode,ielem) = -1000000.0_wp             ! fill partner data with nonsense for BCs
             else
               kelem = ef2e(2,iface,ielem)                               ! Element index of partner node
               kface = ef2e(1,iface,ielem)                               ! Face of partner element
               knode = efn2efn(1,jnode,ielem)                            ! Volumetric index of partner node
               gnode = (kface-1)*nodesperface + ipen
               xgWENO_partner(:,jnode,ielem) = xgWENO_self(:,gnode,kelem) &  ! fill from on-process partner data
                     - xg(:,knode,kelem) + xg(:,inode,ielem)             ! account for possibility of periodic jumps
             end if
           end do

         end if

       end do

     end do

    ! DECODE physical position (xgWENO_partner) into computational coordinates

    extrnal_xi_cnt = 0 ; XI_max = 1.0_wp ; rnorm_max = 0.0_wp ;

    XIWENO_partner(:,:,:) = +0.5_wp

    do ielem = ihelems(1), ihelems(2)            ! loop over volumetric elements

      call WENO_Mapping_Coefs(xg(:,:,ielem),nodesperelem,comp2phys_coeffs)

      do iface = 1,nfacesperelem

        if (ef2e(1,iface,ielem) < 0) then
        else
          do ipen = 1,nodesperface
             jnode = nodesperface*(iface-1) + ipen

             XIWENO_partner(:,jnode,ielem) = xi_guess(xg(:,:,ielem),xgWENO_partner(:,jnode,ielem),nodesperelem)

             call WENO_xi_val(comp2phys_coeffs,XIWENO_partner(:,jnode,ielem),xgWENO_partner(:,jnode,ielem),rnorm)

             if(rnorm >= rnorm_max) rnorm_max = rnorm
             tmp = maxval(abs(XIWENO_partner(:,jnode,ielem)))
             if( tmp > 1.0_wp + tol) then
                if( tmp >= XI_max) XI_max = tmp
                extrnal_xi_cnt = extrnal_xi_cnt + 1
             endif
 
          enddo
        end if

      enddo

    enddo
    extrnal_xi_sum = 0 ;
    call mpi_allreduce(extrnal_xi_cnt,extrnal_xi_sum,1, MPI_INT,MPI_SUM,petsc_comm_world,i_err)

    XI_max_Glob = 0.0_wp ;
    call mpi_allreduce(XI_max,XI_max_Glob,1, MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,i_err)

    if((extrnal_xi_sum > 0) .and. (myprocid == 0)) then
     write(*,*)'number (and XI_max) of SSWENO extrapolants outside -1<=xi<=1 ',extrnal_xi_sum, XI_max_Glob
!    write(*,*)'max residual in nonlinear decode of XI',rnorm_max
    endif

  end subroutine WENO_Intrp_Face_Nodes       !   PARALLEL Routine


  !============================================================================

! pure function Geom_Neighbor_Data(k_node,k_face)

!    !  Grab data at first point orthogonal to surface of the adjoining element
!    !  Element A with 8 vertices is connected to X (and other three other pts) 
!    !  Y is the orthogonal connecting vertex.  
!    !   _______ _______
!    !  |       |       |
!    !  |       |       |
!    !  |   A   |   B   |
!    !  |_______X_______Y
!    !

!    use referencevariables

!    integer,   intent(in) :: k_node,k_face
!    integer               :: Geom_Neighbor_Data

!        if(k_face == 1) then
!        Geom_Neighbor_Data = k_node + nodesperface
!    elseif(k_face == 2) then
!        Geom_Neighbor_Data = k_node + nodesperedge
!    elseif(k_face == 3) then
!        Geom_Neighbor_Data = k_node - 1
!    elseif(k_face == 4) then
!        Geom_Neighbor_Data = k_node - nodesperedge
!    elseif(k_face == 5) then
!        Geom_Neighbor_Data = k_node + 1
!    elseif(k_face == 6) then
!        Geom_Neighbor_Data = k_node - nodesperface
!    endif

! end function Geom_Neighbor_Data

  subroutine Boundary_Vertex_2_Vertex_Connectivity()

    ! Load modules

    use referencevariables, only: nvertices
    use variables, only: if2nq, ifacetag, nqface
    use unary_mod, only: qsorti
!   use variables, only: iaBv2Bv, jaBv2Bv

    !  nqface     = number of boundary faces with quads
    !  if2nq      = face-to-node pointers for each QUAD boundary face (4 nodes)
    !  ifacetag   = tag number that groups boundary faces (e.g., all similar faces)

    ! Nothing is implicitly defined
    implicit none

    integer :: i,k,L,m
    integer :: k2,k3
    integer :: jV,jE
    integer :: icnt0,icnt1,icnt2,icnt3
    integer, dimension(2)  :: itmp

    integer, parameter   :: bigN = 20

    integer, dimension(5)                :: stack
    integer, dimension(bigN)             :: wrk_vec0, wrk_vec1, wrk_vec2, wrk_vec3

    integer                              :: n_bc_types

!   integer, dimension(2,4)              :: if2nq_off

    integer, dimension(:),   allocatable :: iaBv2Bv, ind, nqface_type
    integer, dimension(:,:), allocatable :: jaBv2Bv

    integer, dimension(4,2)              :: if2nqNeigh
    
    allocate(iaBv2Bv(nvertices+1))    ;    iaBv2Bv(:)   = 0
    allocate(jaBv2Bv(nvertices,bigN)) ;    jaBv2Bv(:,:) = 0
 

    ! Sort the BC's into local groups
    allocate(ind(nqface))    ;    ind(:) = 0

    call qsorti(ifacetag,ind,nqface)

    write(*,*)'sorted ifacetag'
    write(*,*)ifacetag(ind(:))

    n_bc_types = 1
    wrk_vec0(:) = 0   ; 

    icnt0 = 1
    do i = 2,nqface
      if(ifacetag(ind(i)) /= ifacetag(ind(i-1))) then
        wrk_vec0(n_bc_types) = icnt0
        n_bc_types = n_bc_types + 1
        icnt0 = 1
      endif
      icnt0 = icnt0 + 1
    enddo
    wrk_vec0(n_bc_types) = icnt0

    write(*,*)'wrk_vec0'   , wrk_vec0

    allocate(nqface_type(n_bc_types+1))
    nqface_type(1) = 1 
    do i = 2,n_bc_types+1
      nqface_type(i) = nqface_type(i-1) + wrk_vec0(i-1)
    enddo

    write(*,*)'n_bc_types',n_bc_types
    write(*,*)'ia_facetypes',nqface_type(:)

    ! Count number of accesses for each vertex over all element Bndry faces within their group
    ! iaBv2Bv :  counts total number of touches to each vertex
    ! jaBv2Bv :  tracks elements that touch vertices

    type_Loop:do k = 1,n_bc_types                                        !  loop over each face type

      iaBv2Bv(:) = 0 ; jaBv2Bv(:,:) = 0
      do jE = nqface_type(k),nqface_type(k+1)-1                          !  loop within a type
        stack(1:4) = if2nq(1:4,ind(jE)) ; stack(5) = if2nq(1,ind(jE)) ;  ! assemble stack of vertices in face
        do i = 1,4
            iaBv2Bv(stack(i+0)) = iaBv2Bv(stack(i+0)) + 1
            iaBv2Bv(stack(i+1)) = iaBv2Bv(stack(i+1)) + 1

            jaBv2Bv(stack(i+0),iaBv2Bv(stack(i+0))) = stack(i+1)
            jaBv2Bv(stack(i+1),iaBv2Bv(stack(i+1))) = stack(i+0)
        end do
      end do                  ! Individual vertices now have a nonunique stack of connected vertices

      write(*,*)'iaBv2Bv',iaBv2Bv(:)

      write(*,*)'connected to vertices'
      do jV = 1,nvertices
         if(iaBv2Bv(jV) == 0) cycle

         call remove_duplicates(iaBv2Bv(jV),jaBv2Bv(jV,:),icnt0,wrk_vec0)

         iaBv2Bv(jV) = icnt0
         jaBv2Bv(jV,1:icnt0) = wrk_vec0(1:icnt0)
!        write(*,*)jV,(jaBv2Bv(jV,L),L=1,iaBv2Bv(jV))
      enddo                !  Individual vertices now have a unique stack of connected verticew


      do jE = nqface_type(k),nqface_type(k+1)-1                !  loop within a type
        stack(1:4) = if2nq(1:4,ind(jE)) ; stack(5) = if2nq(1,ind(jE)) ;  ! assemble stack of vertices in face

        write(*,*)'loop around element face = ',jE
        write(*,*)'stack 1:4 = ',stack(1:4)
        do i = 1,4                                             !  loop around BC face on all four sides
          if2nqNeigh(:,:) = -10000
          icnt0 = 0                                            !  Vertices that touch face of BC quad
          do m = 0,1
            do L = 1,iaBv2Bv(stack(i+m)) 
              icnt0 = icnt0 + 1
              wrk_vec0(icnt0) = jaBv2Bv(stack(i+m),L)
            enddo
          enddo
!         write(*,*)'Before stack on side ',i
!         write(*,*)(wrk_vec0(L),L=1,icnt0)
          call remove_face_values(icnt0,wrk_vec0,if2nq(1:4,ind(jE)),icnt1,wrk_vec1) ! remove quad face vertices 
!         write(*,*)'looping around element = ',i
!         write(*,*)'after stack on side ',i
!         write(*,*)(wrk_vec1(L),L=1,icnt1)

          if(icnt1 >= 2) then
            icnt2 = 0
            wrk_vec2(:) = 100000
            do k2 = 1,icnt1
               do L = 1,iaBv2Bv(wrk_vec1(k2))
                 icnt2 = icnt2 + 1
                 wrk_vec2(icnt2) = jaBv2Bv(wrk_vec1(k2),L)
               enddo
            enddo
                 
!           wrk_vec3 = isort(wrk_vec2,icnt2)
            wrk_vec3 = isort(wrk_vec2,bigN)
            write(*,*)'stack on edge = ',i
            write(*,*)wrk_vec3(:)
            icnt3 = 0
            do k3 = 2,icnt2-1
              if(wrk_vec3(k3-0) == wrk_vec3(k3-1))then
                icnt3 = icnt3 + 1
                if2nqNeigh(i,icnt3) = wrk_vec3(k3)
              endif
            enddo

            ! Is if2nqNeigh stored in correct order?
            Test_Order:do m = 0,1
              do L = 1,iaBv2Bv(stack(i+m)) 
                if(if2nqNeigh(i,1+m) == jaBv2Bv(stack(i+m),L)) then
                  cycle Test_Order
                endif
              enddo
              itmp(:) = if2nqNeigh(i,:) ;
              if2nqNeigh(i,1) = itmp(2) ; 
              if2nqNeigh(i,2) = itmp(1) ;
            enddo Test_Order
          endif

        enddo
        write(*,*)'if2nqNeigh', jE
        write(*,*)1,if2nqNeigh(1,:)
        write(*,*)2,if2nqNeigh(2,:)
        write(*,*)3,if2nqNeigh(3,:)
        write(*,*)4,if2nqNeigh(4,:)
       

      enddo
      
    enddo type_Loop

  !============================================================================

  end subroutine Boundary_Vertex_2_Vertex_Connectivity

  !============================================================================

! subroutine High_Order_Boundary_Connectivity()

!   ! Load modules

!   use referencevariables
!   use variables, only: ef2e, if2nq, ifacetag, ic2nh, nqface

!   ! Nothing is implicitly defined
!   implicit none

!   integer :: i,j,k, j1, j2, k1
!   integer :: ielem,cnt,cntBC
!   integer :: iface, jface, bigN, ave

!   ! low and high volumetric element indices
!   integer :: neighbors(4,2)

!   ! low() and high(2) volumetric element index

!   do ielem = 1, nelems

!     ! loop over faces
!     do iface = 1, nfacesperelem
!        if (ef2e(1,iface,ielem) >= 0) then
!           cycle
!        else
!          do j = 1,nverticesperface                       ! Loop over all face nodes
!            i1 = eltypfaces(j,iface)                      ! Node number on iface
!            do kface = 1, nfacesperelem                   ! Find faces orthogonal to iface
!               if(kface /= iface) then                    ! Don't check myself
!                 do k = 1, nverticesperface               ! Loop over vertices on kface
!                    k1 = eltypfaces(k,kface)              ! Node number on kface
!                    if(i1 /= k1) then
!                      cycle
!                    else
!                      Lface = ef2e(1,kface,ielem)         ! Lface: Face of Adjoining elem
!                      Lelem = ef2e(2,kface,ielem)         ! Lelem:         Adjoining elem
!                      do m = 1, nverticesperface
!                        L1 = eltypfaces(m,Lface)
!                        if(magnitude(
!                      enddo
!                    endif
!                 enddo
!               endif
!            enddo
!        endif
!     enddo
!   enddo

  !============================================================================

!  end subroutine High_Order_Boundary_Connectivity()

  !============================================================================

  subroutine remove_face_values(nI,vecI,faceV,nO,vecO)

    ! Nothing is implicitly defined
    implicit none

    integer,               intent(in)    :: nI
    integer, dimension(:), intent(in)    :: vecI
    integer, dimension(4), intent(in)    :: faceV
    integer,               intent(inout) :: nO
    integer, dimension(:), intent(inout) :: vecO

    integer                              :: L

    nO = 0
    do L = 1,nI
      if( (abs(vecI(L)-faceV(1)) /= 0) .and. &
        & (abs(vecI(L)-faceV(2)) /= 0) .and. &
        & (abs(vecI(L)-faceV(3)) /= 0) .and. &
        & (abs(vecI(L)-faceV(4)) /= 0)) then
          nO = nO + 1
          vecO(nO) = vecI(L)
      endif
    enddo

  end subroutine remove_face_values

  !============================================================================

  subroutine remove_duplicates(nI,iaI,nO,iaO)

    ! Nothing is implicitly defined
    implicit none

    integer,               intent(in)    :: nI
    integer, dimension(:), intent(in)    :: iaI
    integer,               intent(inout) :: nO
    integer, dimension(:), intent(inout) :: iaO

    integer, dimension(nI)                :: iatmp

    integer :: i, icnt

    continue

    iatmp = isort(iaI,nI)
    icnt  = 1
    iaO(icnt) = iatmp(1)
    do i = 2,nI
      if(iaO(icnt) /= iatmp(i)) then
        icnt = icnt + 1
        iaO(icnt) = iatmp(i)
      endif
    enddo
    nO = icnt

  end subroutine remove_duplicates

!=======================================================================================

  function curved_connector_cylinder(nE,x00,x01,x1,x2,xLGL)

    use referencevariables, only: ndim

    implicit none
    integer,                   intent(in) :: nE
    real(wp), dimension(ndim), intent(in) :: x00,x01,x1,x2
    real(wp), dimension(nE),   intent(in) :: xLGL

    real(wp), parameter                   :: tol_o   = 1.0e-10_wp
    real(wp), parameter                   :: tol_r = 1.0e-06_wp


    real(wp), dimension(ndim)             :: dx
    real(wp)                              :: r1, r2
    real(wp), dimension(ndim)             :: origin_1, origin_2
    real(wp)                              :: theta1,theta2, theta
    real(wp)                              :: rmag, thetabar
    real(wp), dimension(ndim,nE)          :: curved_connector_cylinder

    integer                               :: i

    !  Eqn for line at centroid of cylinder:   x = t x00 + (1-t) x01
    !  Solve for t which minimizes distance to x_k, k=1,2 (i.e., the orthogonal projection)

    origin_1 = ortho_projection_to_line3D(x00,x01,x1)
    origin_2 = ortho_projection_to_line3D(x00,x01,x2)

    !  Distance from centroid to point (r1 and r2)
    r1 = magnitude(x1(:)-origin_1(:))
    r2 = magnitude(x2(:)-origin_2(:))

    !  Three cases can occur
    !  1)  r1 = r2 and not parallel to centroid of cylinder  (curve boundary along circumference of circle at radius $r$)
    !  2)  r1 = r2 but parallel with center of cylinder (linear connectivity)
    !  3)  r1 /= r2 (linear connectivity)

    if( (abs(r1-r2) <= tol_r) .and. (magnitude(origin_1-origin_2) <= tol_o))then

      ! General theta
      theta  = acos( dot_product(x1(:)-origin_1(:),x2(:)-origin_2(:))  &
              / magnitude(x1(:)-origin_1(:)) / magnitude(x2(:)-origin_2(:)))

      ! d_theta assumes centroid is along ``z'' axis
      theta1 = atan2(x1(2)-origin_1(2),x1(1)-origin_1(1))
      theta2 = atan2(x2(2)-origin_2(2),x2(1)-origin_2(1))

      if(abs(abs(theta) - abs(theta2-theta1)) >= tol_o) then
        write(*,*)'theta',theta
        write(*,*)'theta1',theta1
        write(*,*)'theta2',theta2
        theta2 = theta1 + theta
        write(*,*)'new theta2',theta2
      endif

      rmag  = abs(0.5_wp*(r1+r2))
      do i = 1,nE
        thetabar = theta1 + (theta2-theta1)*xLGL(i)
        ! d_theta assumes centroid is along ``z'' axis
        curved_connector_cylinder(:,i) = origin_1(:) + rmag*(/cos(thetabar),sin(thetabar),0.0_wp/)
      enddo

    else                            !     write(*,*)'straight'

      do i = 1,nE
        dx(:) = x2(:)-x1(:)
        curved_connector_cylinder(:,i) = x1(:) + xLGL(i)*dx(:)
      enddo

    endif


  end function curved_connector_cylinder

!=======================================================================================

  function curved_connector_parabola(nE,x0_vec,x1_vec,xLGL)

    use referencevariables, only: ndim

    implicit none
    integer,                   intent(in) :: nE
    real(wp), dimension(ndim), intent(in) :: x0_vec,x1_vec
    real(wp), dimension(nE),   intent(in) :: xLGL

    real(wp), parameter                   :: tol_r = 1.0e-06_wp

    real(wp), dimension(ndim)             :: xm_vec, dx_vec
    real(wp)                              :: delta, dydx, dxdy

    real(wp), parameter                   :: amp = 0.1_wp

    real(wp), dimension(ndim,nE)          :: curved_connector_parabola

    real(wp)                              :: aa,bb,cc,de
    real(wp)                              :: x0,xm,xp,x1
    real(wp)                              :: y0,ym,yp,y1
    real(wp)                              :: z0,zm,   z1

    integer                               :: i

    !  Eqn for parabola x = - a t^2  ; y = 2*a*t

    !  Four cases can occur
    !  1)  x0(1) /= x1(1) .and. x0(2) == x1(2) .and. x0(3) == x1(3)   :  Linear
    !  2)  x0(1) == x1(1) .and. x0(2) /= x1(2) .and. x0(3) == x1(3)   :  Linear
    !  3)  x0(1) == x1(1) .and. x0(2) == x1(2) .and. x0(3) /= x1(3)   :  Linear
    !  4)  x0(1) /= x1(1) .and. x0(2) /= x1(2) .and. x0(3) == x1(3)   :
    !  Parabola

    xm_vec(:) = 0.5_wp * (x0_vec(:) + x1_vec(:))

    x0 = x0_vec(1); y0 = x0_vec(2); z0 = x0_vec(3)
    x1 = x1_vec(1); y1 = x1_vec(2); z1 = x1_vec(3)
    xm = xm_vec(1); ym = xm_vec(2); zm = xm_vec(3)

    if( (abs(x1 - x0) >= tol_r) .and.  &
        (abs(y1 - y0) >= tol_r) .and.  &
        (abs(z1 - z0) <= tol_r) ) then

      delta   = magnitude(x0_vec(:) - x1_vec(:))

      dydx    = (y1 - y0) / (x1 - x0)
      dxdy    = 1.0_wp / dydx

      if(abs(dydx) >= abs(dxdy)) then   !  vertical   connectors
        xp     = xm +  sqrt((amp*delta)**2 / (1.0_wp + dxdy**2))
        yp     = ym - dxdy * (xp - xm)
        do i = 1,nE
          curved_connector_parabola(2,i) =  y0 + (y1 - y0) * xLGL(i)
          curved_connector_parabola(1,i) =  x0                                  &
          & + 1.0_wp * (-1.0_wp*x1 + 4.0_wp * xp - 3.0_wp*x0) * xLGL(i)**1  &
          & + 2.0_wp * (+1.0_wp*x1 - 2.0_wp * xp + 1.0_wp*x0) * xLGL(i)**2
          curved_connector_parabola(3,i) = zm
          de = ((y0 - y1)*(y0 - yp)*(y1 - yp))
          aa = (yp*(-x0 + x1) + y1*(x0 - xp) + y0*(-x1 + xp))/de
          bb = (yp**2*(x0 - x1) + y0**2*(x1 - xp) + y1**2*(-x0 + xp))/de
          cc = (yp*(y1*(y1 - yp)*x0 + y0*(-y0 + yp)*x1) + y0*(y0 - y1)*y1*xp)/de
!         write(*,*)'vert', maxval(abs(curved_connector_parabola(1,:)  &
!         - ( aa*curved_connector_parabola(2,:)**2 + bb*curved_connector_parabola(2,:) + cc ) ))
        enddo

      else                              !  horizontal connectors
        yp  = ym +  sqrt((amp*delta)**2 / (1.0_wp + dydx**2))
        xp  = xm - dydx * (yp - ym)
        do i = 1,nE
          curved_connector_parabola(1,i) =  x0 + (x1 - x0) * xLGL(i)
          curved_connector_parabola(2,i) =  y0                              &
          & + 1.0_wp * (-1.0_wp*y1 + 4.0_wp * yp - 3.0_wp*y0) * xLGL(i)**1  &
          & + 2.0_wp * (+1.0_wp*y1 - 2.0_wp * yp + 1.0_wp*y0) * xLGL(i)**2
          curved_connector_parabola(3,i) = zm
          de = ((x0 - x1)*(x0 - xp)*(x1 - xp))
          aa = (xp*(-y0 + y1) + x1*(y0 - yp) + x0*(-y1 + yp))/de
          bb = (xp**2*(y0 - y1) + x0**2*(y1 - yp) + x1**2*(-y0 + yp))/de
          cc = (xp*(x1*(x1 - xp)*y0 + x0*(-x0 + xp)*y1) + x0*(x0 - x1)*x1*yp)/de
!         write(*,*)'horz', maxval(abs(curved_connector_parabola(2,:)  &
!         - ( aa*curved_connector_parabola(1,:)**2 + bb*curved_connector_parabola(1,:) + cc ) ))
        enddo
      endif



    else                !     write(*,*)'linear'

      do i = 1,nE
        dx_vec(:) = x1_vec(:)-x0_vec(:)
        curved_connector_parabola(:,i) = x0_vec(:) + xLGL(i)*dx_vec(:)
      enddo

    endif

  end function curved_connector_parabola


  pure function ortho_projection_to_line3D(x00,x01,x)

    real(wp), dimension(3), intent(in) :: x00,x01,x

    real(wp)                           :: t
    real(wp), dimension(3)             :: ortho_projection_to_line3D

    t = dot_product(x01-x,x01-x00) / dot_product(x01-x00,x01-x00) ;

    ortho_projection_to_line3D = t*x00 + (1.0_wp - t)*x01 ;

  end function ortho_projection_to_line3D

  !============================================================================

  pure function LGL_pts_lexo_comp_hexa(x_1d,n_pts_1d)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_pts_1d
    real(wp), dimension(n_pts_1d), intent(in) :: x_1d
    real(wp), dimension(3,n_pts_1d**3) :: LGL_pts_lexo_comp_hexa
    real(wp) :: zeta, eta, xi

    integer :: i, j, k, node_id

    continue

    ! Initialize to zero the coordinates
    LGL_pts_lexo_comp_hexa(:,:) = 0.0_wp

    ! Set to zero the local node ID
    node_id = 0

    ! Set the coordinates of the Gauss points in computational space
    do i = 1, n_pts_1d
      zeta = x_1d(i)
      do j = 1, n_pts_1d
         eta = x_1d(j)
        do k = 1, n_pts_1d
            xi = x_1d(k)
          node_id = node_id + 1
          LGL_pts_lexo_comp_hexa(3,node_id) = zeta
          LGL_pts_lexo_comp_hexa(2,node_id) =  eta
          LGL_pts_lexo_comp_hexa(1,node_id) =   xi
        end do
      end do
    end do

  end function LGL_pts_lexo_comp_hexa

  !=================================================================================================

  subroutine calc_Gau_shell_pts_all_hexas()

    ! Load module
    use variables
    use referencevariables
    use collocationvariables
    use initcollocation, only : lagrange_basis_function_1d, Gauss_Legendre_points
    use initcollocation, only : element_properties

    ! Nothing is implicitly defined
    implicit none

    real(wp), allocatable, dimension(:,:) :: Gau_pts_comp_shell_one_face

    integer  :: ielem, i_Gau, i_LGL, j_LGL, k_LGL
    real(wp) :: l_xi, l_eta, l_zeta
    real(wp) :: xi_in, eta_in, zeta_in
    real(wp),  dimension(:), allocatable :: x_Gau_1d_Mort,w_Gau_1d_Mort
    real(wp),  dimension(:), allocatable :: x_LGL_1d_On

    integer :: l, ishift, iface
    integer :: n_Gau_shell_max, n_Gau_1d_max, n_Gau_2d_max, n_Gau_1d_Mort, n_Gau_2d_Mort
    integer :: n_LGL_1d_On, n_LGL_1d_Off

    continue

    ! Set the maximum size of buckets for shell data 
    n_Gau_1d_max = (npoly_max+1)**1
    n_Gau_2d_max = (npoly_max+1)**2
    n_Gau_shell_max = nfacesperelem * n_Gau_2d_max

    ! Allocate memory
    allocate(xg_Gau_Shell(3,n_Gau_shell_max,ihelems(1):ihelems(2))) ;  xg_Gau_Shell(:,:,:) = 0.0_wp
    allocate(Gau_pts_comp_shell_one_face(3,n_Gau_2d_max))           ;  Gau_pts_comp_shell_one_face(:,:) = 0.0_wp

    ! Loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,               &
                              n_pts_1d=n_LGL_1d_On, &
                              x_pts_1d=x_LGL_1d_On)

      do iface = 1,nfacesperelem

        n_LGL_1d_Off  = ef2e(4,iface,ielem)
        n_Gau_1d_Mort = max(n_LGL_1d_On, n_LGL_1d_Off)
        n_Gau_2d_Mort = n_Gau_1d_Mort**2

!-- MAKE SURE THIS IS NEEDED
!        if(elem_props(2,ielem) == ef2e(4,iface,ielem)) then      ! cycle for conforming interfaces
!
!            cycle
!        endif
!--

        if(allocated(x_Gau_1d_Mort)) deallocate(x_Gau_1d_Mort) ; allocate(x_Gau_1d_Mort(n_Gau_1d_Mort)) ;
        if(allocated(w_Gau_1d_Mort)) deallocate(w_Gau_1d_Mort) ; allocate(w_Gau_1d_Mort(n_Gau_1d_Mort)) ;

        call Gauss_Legendre_points(n_Gau_1d_Mort,x_Gau_1d_Mort,w_Gau_1d_Mort)

        Gau_pts_comp_shell_one_face = Shell_pts_lexo_comp_one_face(iface,n_Gau_2d_max, &
                                                                   n_Gau_1d_mort,x_Gau_1d_Mort)

        do i_Gau = 1, n_Gau_2d_Mort

          ishift = (iface-1) * n_Gau_2d_max + i_Gau
          l = 0
          xg_Gau_shell(:,ishift,ielem) = 0.0_wp
  
            xi_in = Gau_pts_comp_shell_one_face(1,i_Gau)
           eta_in = Gau_pts_comp_shell_one_face(2,i_Gau)
          zeta_in = Gau_pts_comp_shell_one_face(3,i_Gau)
  
          do k_LGL = 1, n_LGL_1d_On
  
            l_zeta = lagrange_basis_function_1d(zeta_in,k_LGL,x_LGL_1d_On,n_LGL_1d_On)
  
            do j_LGL = 1, n_LGL_1d_On
  
              l_eta  = lagrange_basis_function_1d(eta_in, j_LGL,x_LGL_1d_On,n_LGL_1d_On)
  
              do i_LGL = 1, n_LGL_1d_On
                l = l + 1
  
                l_xi   = lagrange_basis_function_1d(xi_in,  i_LGL,x_LGL_1d_On,n_LGL_1d_On)
  
                xg_Gau_shell(:,ishift,ielem) = xg_Gau_shell(:,ishift,ielem) + xg(:,l,ielem)*l_xi*l_eta*l_zeta
  
              end do
            end do
          end do
        end do

      end do

    end do

    deallocate(x_LGL_1d_On,x_Gau_1d_Mort,w_Gau_1d_Mort) ; 
    deallocate(Gau_pts_comp_shell_one_face) ;

  end subroutine calc_Gau_shell_pts_all_hexas

  !============================================================================

  subroutine calc_Jacobian_Gau_shell_all_hexas()

    ! Load module
    use variables
    use referencevariables
    use collocationvariables
    use initcollocation, only : lagrange_basis_function_1d, Gauss_Legendre_points
    use initcollocation, only : element_properties

    ! Nothing is implicitly defined
    implicit none

    real(wp), allocatable, dimension(:,:) :: Gau_pts_comp_shell_one_face

    integer  :: ielem, i_Gau, i_LGL, j_LGL, k_LGL
    real(wp) :: l_xi, l_eta, l_zeta
    real(wp) :: xi_in, eta_in, zeta_in
    real(wp),  dimension(:), allocatable :: x_Gau_1d_Mort,w_Gau_1d_Mort
    real(wp),  dimension(:), allocatable :: x_LGL_1d_On

    integer :: l, ishift, iface
    integer :: n_Gau_shell_max, n_Gau_1d_max, n_Gau_2d_max, n_Gau_1d_Mort, n_Gau_2d_Mort
    integer :: n_LGL_1d_On, n_LGL_1d_Off

    continue

    ! Set the maximum size of buckets for shell data 
    n_Gau_1d_max = (npoly_max+1)**1
    n_Gau_2d_max = (npoly_max+1)**2
    n_Gau_shell_max = nfacesperelem * n_Gau_2d_max

    ! Allocate memory
    allocate(Jx_r_Gau_Shell(n_Gau_shell_max,ihelems(1):ihelems(2))) ;  Jx_r_Gau_Shell(:,:) = 0.0_wp
    allocate(Gau_pts_comp_shell_one_face(3,n_Gau_2d_max))           ;  Gau_pts_comp_shell_one_face(:,:) = 0.0_wp

    ! Loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,               &
                              n_pts_1d=n_LGL_1d_On, &
                              x_pts_1d=x_LGL_1d_On)

      do iface = 1,nfacesperelem

        n_LGL_1d_Off  = ef2e(4,iface,ielem)
        n_Gau_1d_Mort = max(n_LGL_1d_On, n_LGL_1d_Off)
        n_Gau_2d_Mort = n_Gau_1d_Mort**2

        if(allocated(x_Gau_1d_Mort)) deallocate(x_Gau_1d_Mort) ; allocate(x_Gau_1d_Mort(n_Gau_1d_Mort)) ;
        if(allocated(w_Gau_1d_Mort)) deallocate(w_Gau_1d_Mort) ; allocate(w_Gau_1d_Mort(n_Gau_1d_Mort)) ;

        call Gauss_Legendre_points(n_Gau_1d_Mort,x_Gau_1d_Mort,w_Gau_1d_Mort)

        Gau_pts_comp_shell_one_face = Shell_pts_lexo_comp_one_face(iface,n_Gau_2d_max, &
                                                                   n_Gau_1d_mort,x_Gau_1d_Mort)

        do i_Gau = 1, n_Gau_2d_Mort

          ishift = (iface-1) * n_Gau_2d_max + i_Gau
          l = 0
          Jx_r_Gau_shell(ishift,ielem) = 0.0_wp
  
            xi_in = Gau_pts_comp_shell_one_face(1,i_Gau)
           eta_in = Gau_pts_comp_shell_one_face(2,i_Gau)
          zeta_in = Gau_pts_comp_shell_one_face(3,i_Gau)
  
          do k_LGL = 1, n_LGL_1d_On
  
            l_zeta = lagrange_basis_function_1d(zeta_in,k_LGL,x_LGL_1d_On,n_LGL_1d_On)
  
            do j_LGL = 1, n_LGL_1d_On
  
              l_eta  = lagrange_basis_function_1d(eta_in, j_LGL,x_LGL_1d_On,n_LGL_1d_On)
  
              do i_LGL = 1, n_LGL_1d_On
                l = l + 1
  
                l_xi   = lagrange_basis_function_1d(xi_in,  i_LGL,x_LGL_1d_On,n_LGL_1d_On)
  
                Jx_r_Gau_shell(ishift,ielem) = Jx_r_Gau_shell(ishift,ielem) + Jx_r(l,ielem)*l_xi*l_eta*l_zeta
  
              end do
            end do
          end do
        end do

      end do

    end do

    deallocate(x_LGL_1d_On,x_Gau_1d_Mort,w_Gau_1d_Mort) ; 
    deallocate(Gau_pts_comp_shell_one_face) ;

  end subroutine calc_Jacobian_Gau_shell_all_hexas

  !============================================================================

  pure function Shell_pts_lexo_comp_quad(x_pts_1d,n_pts_1d)

    ! Nothing is implicitly defined
    use referencevariables

    implicit none

    integer, intent(in) :: n_pts_1d
    real(wp), dimension(n_pts_1d), intent(in) :: x_pts_1d
    real(wp), dimension(3,nfacesperelem*n_pts_1d**2) :: Shell_pts_lexo_comp_quad
    real(wp) :: zeta, eta, xi

    integer :: i, j, k, node_id, iface
    integer :: ixL, ixH, iyL, iyH, izL, izH

    continue

    ! Initialize to zero the coordinates
    Shell_pts_lexo_comp_quad(:,:) = 0.0_wp

    ! Set to zero the local node ID
    node_id = 0

    do iface=1,nfacesperelem
      ! Set the coordinates of the Gauss points in computational space
          if(iface == 1) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = 1
        zeta = -1.0_wp
          do j = iyL, iyH
            eta = x_pts_1d(j)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_quad(3,node_id) = zeta
              Shell_pts_lexo_comp_quad(2,node_id) = eta
              Shell_pts_lexo_comp_quad(1,node_id) = xi
            end do
          end do
      elseif(iface == 2) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = 1
         izL = 1   ;  izH = n_pts_1d
         eta = -1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_quad(3,node_id) = zeta
              Shell_pts_lexo_comp_quad(2,node_id) = eta
              Shell_pts_lexo_comp_quad(1,node_id) = xi
            end do
        end do
      elseif(iface == 3) then
         ixL = 1   ;  ixH = 1
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = n_pts_1d
         xi  = +1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
          do j = iyL, iyH
            eta = x_pts_1d(j)
              node_id = node_id + 1
              Shell_pts_lexo_comp_quad(3,node_id) = zeta
              Shell_pts_lexo_comp_quad(2,node_id) = eta
              Shell_pts_lexo_comp_quad(1,node_id) = xi
          end do
        end do
      elseif(iface == 4) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = 1
         izL = 1   ;  izH = n_pts_1d
         eta = +1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_quad(3,node_id) = zeta
              Shell_pts_lexo_comp_quad(2,node_id) = eta
              Shell_pts_lexo_comp_quad(1,node_id) = xi
            end do
        end do
      elseif(iface == 5) then
         ixL = 1   ;  ixH = 1
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = n_pts_1d
          xi = -1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
          do j = iyL, iyH
            eta = x_pts_1d(j)
              node_id = node_id + 1
              Shell_pts_lexo_comp_quad(3,node_id) = zeta
              Shell_pts_lexo_comp_quad(2,node_id) = eta
              Shell_pts_lexo_comp_quad(1,node_id) = xi
          end do
        end do
      elseif(iface == 6) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = 1
        zeta = +1.0_wp
          do j = iyL, iyH
            eta = x_pts_1d(j)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_quad(3,node_id) = zeta
              Shell_pts_lexo_comp_quad(2,node_id) = eta
              Shell_pts_lexo_comp_quad(1,node_id) = xi
            end do
          end do
        endif
      enddo     !  face do loop 

    return
  end function Shell_pts_lexo_comp_quad

  !============================================================================

  pure function Shell_pts_lexo_comp_one_face(iface,n_max_2d,n_pts_1d,x_pts_1d)

    ! Nothing is implicitly defined
    use referencevariables

    implicit none

    integer,                       intent(in) :: iface,n_max_2d,n_pts_1d
    real(wp), dimension(n_pts_1d), intent(in) :: x_pts_1d

    real(wp), dimension(3,n_max_2d)           :: Shell_pts_lexo_comp_one_face

    real(wp) :: zeta, eta, xi

    integer :: i, j, k, node_id
    integer :: ixL, ixH, iyL, iyH, izL, izH

    continue

    ! Initialize to zero the coordinates
    Shell_pts_lexo_comp_one_face(:,:) = 0.0_wp

    ! Set to zero the local node ID
    node_id = 0

      ! Set the coordinates of the Gauss points in computational space
          if(iface == 1) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = 1
        zeta = -1.0_wp
          do j = iyL, iyH
            eta = x_pts_1d(j)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_one_face(3,node_id) = zeta
              Shell_pts_lexo_comp_one_face(2,node_id) = eta
              Shell_pts_lexo_comp_one_face(1,node_id) = xi
            end do
          end do
      elseif(iface == 2) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = 1
         izL = 1   ;  izH = n_pts_1d
         eta = -1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_one_face(3,node_id) = zeta
              Shell_pts_lexo_comp_one_face(2,node_id) = eta
              Shell_pts_lexo_comp_one_face(1,node_id) = xi
            end do
        end do
      elseif(iface == 3) then
         ixL = 1   ;  ixH = 1
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = n_pts_1d
         xi  = +1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
          do j = iyL, iyH
            eta = x_pts_1d(j)
              node_id = node_id + 1
              Shell_pts_lexo_comp_one_face(3,node_id) = zeta
              Shell_pts_lexo_comp_one_face(2,node_id) = eta
              Shell_pts_lexo_comp_one_face(1,node_id) = xi
          end do
        end do
      elseif(iface == 4) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = 1
         izL = 1   ;  izH = n_pts_1d
         eta = +1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_one_face(3,node_id) = zeta
              Shell_pts_lexo_comp_one_face(2,node_id) = eta
              Shell_pts_lexo_comp_one_face(1,node_id) = xi
            end do
        end do
      elseif(iface == 5) then
         ixL = 1   ;  ixH = 1
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = n_pts_1d
          xi = -1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
          do j = iyL, iyH
            eta = x_pts_1d(j)
              node_id = node_id + 1
              Shell_pts_lexo_comp_one_face(3,node_id) = zeta
              Shell_pts_lexo_comp_one_face(2,node_id) = eta
              Shell_pts_lexo_comp_one_face(1,node_id) = xi
          end do
        end do
      elseif(iface == 6) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = 1
        zeta = +1.0_wp
          do j = iyL, iyH
            eta = x_pts_1d(j)
            do k = ixL, ixH
              xi = x_pts_1d(k)
              node_id = node_id + 1
              Shell_pts_lexo_comp_one_face(3,node_id) = zeta
              Shell_pts_lexo_comp_one_face(2,node_id) = eta
              Shell_pts_lexo_comp_one_face(1,node_id) = xi
            end do
          end do
        endif

    return
  end function Shell_pts_lexo_comp_one_face

  !============================================================================

  subroutine facenodesetup_Gau(n_Gau_1d, n_Gau_2d, kfacenodes_Gau, ifacenodes_Gau)
    !  Establish pointers for grabbing the face nodes from an 
    !  element volumetic ordering.  Implicit are that a subset of the volume nodes 
    !  are on the face of the element
    !  
    !   kfacenodes(n_Gau_2d_p1,nfacesperelem)  
    !      volumetric node index of face node  
    !      
    !   ifacenodes(n_Gau_2d_p1*nfacesperelem)  
    !      kfacenode flattened into a single vector
    !  
    use referencevariables
    use mpimod

    implicit none

    integer,                              intent(in   ) :: n_Gau_1d, n_Gau_2d
    integer, dimension(:,:), allocatable, intent(inout) :: kfacenodes_Gau
    integer, dimension(:),   allocatable, intent(inout) :: ifacenodes_Gau

    ! indices
    integer :: i,j,k

    real(wp), parameter :: nodetol = 1.0e-8_wp

    ! local facial masks
    if(allocated(kfacenodes_Gau)) deallocate(kfacenodes_Gau) ; allocate(kfacenodes_Gau(n_Gau_2d,nfacesperelem))
    if(allocated(ifacenodes_Gau)) deallocate(ifacenodes_Gau) ; allocate(ifacenodes_Gau(n_Gau_2d*nfacesperelem))
     
    if (ndim == 2) then
    else if (ndim == 3) then
      do k = 1, n_Gau_2d
          ! face 1 does not require an offset or a stride
          kfacenodes_Gau(k,1) = (1-1)*n_Gau_2d + k
          kfacenodes_Gau(k,2) = (2-1)*n_Gau_2d + k
          kfacenodes_Gau(k,3) = (3-1)*n_Gau_2d + k
          kfacenodes_Gau(k,4) = (4-1)*n_Gau_2d + k
          kfacenodes_Gau(k,5) = (5-1)*n_Gau_2d + k
          kfacenodes_Gau(k,6) = (6-1)*n_Gau_2d + k
      end do
    else
      write(*,*) 'error: unsupported dimension', ndim
      stop
    end if
    ! the nodes on each face are packed into a 1D array
    k = 0
    ! loop over faces
    do j = 1, nfacesperelem
      ! loop over nodes on each face
      do i = 1, n_Gau_2d
        ! advance facial node index
        k = k+1
        ! map facial node index to volumetric node
        ifacenodes_Gau(k) = kfacenodes_Gau(i,j)
      end do
    end do

  end subroutine facenodesetup_Gau

  !============================================================================
  ! calculate_face_node_connectivity - Sets the face-node connectivity for the
  ! collocation points.
  !============================================================================

  subroutine calculate_face_node_connectivity_Gau()
    
    ! Load modules
    use referencevariables
    use collocationvariables
    use mpimod
    use variables,       only: ef2e, xg_Gau_shell, xgghst_Gau_Shell, efn2efn_Gau
    use initcollocation, only: element_properties

    ! Nothing is implicitly defined
    implicit none

    real(wp), parameter :: nodetol = 1.0e-8_wp

    integer ::  ielem, inode, jnode, iface, knode, kshell
    integer ::  i_low, ierr

    real(wp) :: x1(3), x2(3)
    real(wp) :: distance_min

    integer :: n_Gau_1d_max , n_Gau_2d_max , n_Gau_shell_max
    integer :: n_Gau_1d_Mort, n_Gau_2d_Mort
    integer :: n_LGL_1d_On  , n_LGL_1d_Off

    continue

    ! efn2efn_Gau contains the partner node information of every facenode in the domain

    ! efn2efn_Gau(1,k,m) = Volumetric ID of partner node  (not used)
    ! efn2efn_Gau(2,k,m) = element ID of partner node
    ! efn2efn_Gau(3,k,m) = Stack pointer for PETSC ghost array             (defined only for    parallel interface)
    ! efn2efn_Gau(4,k,m) = node value in shell coordinates of partner node (defined only if non-parallel interface)
    !  NOTE:  The definition of efn2efn_Gau(4) is different than efn2efn(4)

    ! Set the maximum size of buckets for shell data 
    n_Gau_1d_max = (npoly_max+1)**1
    n_Gau_2d_max = (npoly_max+1)**2
    n_Gau_shell_max = nfacesperelem * n_Gau_2d_max

    allocate(efn2efn_Gau(4,n_Gau_shell_max,ihelems(1):ihelems(2))) ; efn2efn_Gau = -1000 ;

    i_low = 0                                                    ! Initialize position of the ghost point in the stack

    elem_loop:do ielem = ihelems(1), ihelems(2)                  ! loop over elements

      call element_properties(ielem,         &
                     n_pts_1d=n_LGL_1d_On,   &
                     n_pts_2d=nodesperface,  &
                     n_pts_3d=nodesperelem)

      face_loop:do iface = 1, nfacesperelem                      ! loop over faces of element

        n_LGL_1d_Off  = ef2e(4,iface,ielem)
        n_Gau_1d_Mort = max(n_LGL_1d_On, n_LGL_1d_Off)
        n_Gau_2d_Mort = n_Gau_1d_Mort**2

        knode = (iface-1) * n_Gau_2d_max

        if(elem_props(2,ielem) == ef2e(4,iface,ielem)) then      ! cycle for conforming interfaces

            cycle

        else if((ef2e(3,iface,ielem) /= myprocid) .and. (elem_props(2,ielem) /= ef2e(4,iface,ielem))) then ! Off-process - NonConforming

          do inode = 1, n_Gau_2d_Mort                            ! Loop over the nodes on the mortar

            knode = (iface-1) * n_Gau_2d_max + inode             ! location in the on-element stack

            x1 = xg_Gau_shell(:,knode,ielem)                     ! Save the coordinates of the facial node
              
            distance_min = +10000000                             ! initial all distances to huge number

            do jnode = 1, n_Gau_2d_Mort                          ! Search for the connected node on face of the connected element

              x2 = xgghst_Gau_Shell(:,i_low + jnode)             ! ef2e(2) gives the element of the neighbor
              
              distance_min = min(distance_min,magnitude(x1-x2))  ! distance between two points

            enddo

            do jnode = 1,n_Gau_2d_Mort

              x2 = xgghst_Gau_Shell(:,i_low + jnode)             ! ef2e(2) gives the element of the neighbor
              
              if(magnitude(x1-x2) <= distance_min + nodetol) then
                
                efn2efn_Gau(1,knode,ielem) = -1000               ! Set the volumetric node index of the connected node
                
                efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem) ! Set the element of the connected node
                
                efn2efn_Gau(3,knode,ielem) = i_low + jnode       ! Set the node index in the ghost array
!               efn2efn_Gau(3,knode,ielem) = jnode

                exit
              
              end if
            
            end do
              
            if (jnode > n_Gau_2d_Mort .and. myprocid==1) then       ! Print information at screen if there is a problem and stop computation

              write(*,*) 'Connectivity error in face-node connectivity_Gau, Parallel.'
              write(*,*) 'Process ID, element ID, face ID, ef2e'
              write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
              write(*,*) 'Node coordinates'
              write(*,*) x1
              write(*,*) 'Possible partner node coordinates'
              
              do jnode = i_low+1, i_low+n_Gau_2d_Mort
                x2 = xgghst_Gau_Shell(:,jnode)
                write(*,*) x2
              end do 
              write(*,*) 'Exiting...'
              call PetscFinalize(ierr) ; stop
            end if

          end do

          i_low = i_low + n_Gau_2d_Mort                             ! Update the position in the ghost stack
          
        else if((ef2e(3,iface,ielem) == myprocid) .and. (elem_props(2,ielem) /= ef2e(4,iface,ielem))) then  ! On-process - Non-conforming

          do inode = 1, n_Gau_2d_Mort                               ! Loop over the nodes on the face

            knode = (iface-1) * n_Gau_2d_max + inode                ! location in the on-element stack

            x1 = xg_Gau_shell(:,knode,ielem)                        ! Save coordinates of the facial ndoes

            distance_min = +10000000                                ! initial all distances to huge number

            do jnode = 1, n_Gau_2d_Mort                             ! Search the for connected node on the face of the connected element
              
              kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_max + jnode ! ef2e(1) gives the face on the neighboring element 

              x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))       ! ef2e(2) gives the element
              
              distance_min = min(distance_min,magnitude(x1-x2))     ! distance between two points

            enddo

            do jnode = 1,n_Gau_2d_Mort

              kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_max + jnode ! ef2e(1) gives the face on the neighboring element 

              x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))       ! ef2e(2) gives the element

              if(magnitude(x1-x2) <= distance_min + nodetol) then

                efn2efn_Gau(1,knode,ielem) = -1000                  ! Set the volumetric node index of the connected node
                
                efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)    ! Set the element of the connected node
                
                efn2efn_Gau(4,knode,ielem) = kshell                 ! Set the index of the connected node
                
                exit
              
              end if
            
            end do ! End do jnode

            ! Print information at screen if there is a problem and stop computation

            if (efn2efn_Gau(2,knode,ielem) < 0) then
              write(*,*) 'Connectivity error in face-node connectivity of Gauss path.'
              write(*,*) 'Process ID, element ID, face ID, ef2e'
              write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
              write(*,*) 'Node coordinates'
              write(*,*) x1
              write(*,*) 'Possible partner node coordinates'
              
              do jnode = 1, n_Gau_2d_Mort
                kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_max + jnode
                x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))
                write(*,*) x2
              end do 
              write(*,*) 'Exiting...'
              call PetscFinalize(ierr) ; stop
            end if

          end do ! End do inode
          
        end if ! End if type of face (boundary, off processor or on processor)
      
      end do face_loop! End do loop over faces of the element

    end do elem_loop   !End do loop elements owned by the processor

  end subroutine calculate_face_node_connectivity_Gau

  !============================================================================

  subroutine Shell_Metrics_Analytic(iface, n_pts_1d_max, n_pts_1d, x_pts_1d,  &
                                    xg_Gau,Jx_facenodenormal_Gau)

    ! Nothing is implicitly defined
    use referencevariables
    use initcollocation,      only: D_lagrange_basis_function_1d, &
                                      lagrange_basis_function_1d

    implicit none

    integer,                       intent(in   ) :: iface, n_pts_1d_max, n_pts_1d
    real(wp), dimension(n_pts_1d), intent(in   ) :: x_pts_1d
    real(wp), dimension(:,:),      intent(in   ) :: xg_Gau
    real(wp), dimension(:,:),      intent(  out) :: Jx_facenodenormal_Gau
    real(wp) :: zeta, eta, xi
    real(wp) :: t1, t2, t3, t4

    integer  :: i, j, k, i1, j1, k1
    integer  :: ixL, ixH, iyL, iyH, izL, izH
    integer  :: n_Gau, node_id
    integer  :: ival, kval
    integer  :: ishift, jshift, kshift, face_shift

    real(wp), dimension(:), allocatable  :: x_Gau
    
    continue

    node_id = 0
    n_Gau = n_pts_1d
    face_shift = (iface-1)*n_pts_1d_max**2
    allocate(x_Gau(n_pts_1d)) ; x_Gau(:) = x_pts_1d

    ! Set to zero the local node ID
          if(iface == 1) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = 1
        zeta = -1.0_wp
          do j = iyL, iyH
            eta = x_pts_1d(j)
            do k = ixL, ixH
              xi = x_pts_1d(k)

              node_id = node_id + 1
              kval = mod(node_id-1,n_Gau)+1 ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval + face_shift
                 t1 =  t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 =  t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 =  t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t3 =  t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo   
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  (t1*t2 - t3*t4) * zeta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  (t1*t2 - t3*t4) * zeta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  (t1*t2 - t3*t4) * zeta
            end do
          end do
      elseif(iface == 2) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = 1
         izL = 1   ;  izH = n_pts_1d
         eta = -1.0_wp
        do i = izL, izH
          zeta = x_Gau(i)
            do k = ixL, ixH
              xi = x_Gau(k)

              node_id = node_id + 1
              kval = mod(node_id-1,n_Gau)+1 ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  (t1*t2 - t3*t4) * eta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  (t1*t2 - t3*t4) * eta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  (t1*t2 - t3*t4) * eta
            end do
        end do
      elseif(iface == 3) then
         ixL = 1   ;  ixH = 1
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = n_pts_1d
         xi  = +1.0_wp
        do i = izL, izH
          zeta = x_Gau(i)
            do j = iyL, iyH
              eta = x_Gau(j)

              node_id = node_id + 1
              ival = mod(node_id-1,n_Gau)+1 ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = node_id - ival + j1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + ival + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  (t1*t2 - t3*t4) * xi

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = node_id - ival + j1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + ival + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  (t1*t2 - t3*t4) * xi

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = node_id - ival + j1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + ival + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  (t1*t2 - t3*t4) * xi
          end do
        end do
      elseif(iface == 4) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = 1
         izL = 1   ;  izH = n_pts_1d
         eta = +1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
            do k = ixL, ixH
              xi = x_pts_1d(k)

              node_id = node_id + 1
              kval = mod(node_id-1,n_Gau)+1 ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  (t1*t2 - t3*t4) * eta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  (t1*t2 - t3*t4) * eta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  (t1*t2 - t3*t4) * eta
            end do
        end do
      elseif(iface == 5) then
         ixL = 1   ;  ixH = 1
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = n_pts_1d
          xi = -1.0_wp
        do i = izL, izH
          zeta = x_pts_1d(i)
          do j = iyL, iyH
            eta = x_pts_1d(j)

              node_id = node_id + 1
              ival = mod(node_id-1,n_Gau)+1 ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = node_id - ival + j1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + ival + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  (t1*t2 - t3*t4) * xi

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = node_id - ival + j1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + ival + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  (t1*t2 - t3*t4) * xi

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = node_id - ival + j1 + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + ival + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  (t1*t2 - t3*t4) * xi
          end do
        end do
      elseif(iface == 6) then
         ixL = 1   ;  ixH = n_pts_1d
         iyL = 1   ;  iyH = n_pts_1d
         izL = 1   ;  izH = 1
        zeta = +1.0_wp
          do j = iyL, iyH
            eta = x_pts_1d(j)
            do k = ixL, ixH
              xi = x_pts_1d(k)

              node_id = node_id + 1
              kval = mod(node_id-1,n_Gau)+1 ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval + face_shift
                 t1 =  t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 =  t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 =  t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t3 =  t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo   
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  (t1*t2 - t3*t4) * zeta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  (t1*t2 - t3*t4) * zeta

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  (t1*t2 - t3*t4) * zeta
            end do
          end do
      endif

    return
  end subroutine Shell_Metrics_Analytic
  
  subroutine transform_grid()
  !================================================================================================
  !
  ! Purpose: takes the grid and subjects it to a transformation
  !
  !================================================================================================
 
    !-- use statments
    use mpimod
    use referencevariables, only                 : myprocid, nprocs, ihelems, npoly_max, nfacesperelem
    use variables, only                          : xg, xg_Gau_shell, ef2e
    use precision_vars, only                     : pi
    use initcollocation, only                    : element_properties
  
    implicit none

    !-- local variables
    integer                                      :: m_size, iproc,i_err
    integer                                      :: s_tag, s_request_err, s_status
    integer                                      :: r_tag, r_request_err, r_status
    integer                                      :: inode, ielem
    integer                                      :: n_pts_3d
    integer                                      :: iface, ishift, i_Gau, n_Gau_2d_max
    integer                                      :: n_Gau_1d_Mort, n_Gau_2d_Mort
    integer                                      :: n_LGL_1d_On, n_LGL_1d_Off, order
    real(wp)                                     :: x_min, y_min, z_min
    real(wp)                                     :: x_max, y_max, z_max
    real(wp), allocatable,dimension(:,:,:)       :: xyz_minmax
    real(wp), dimension(3,2)                     :: xyz_loc_minmax
    real(wp)                                     :: x_old, y_old, z_old, x, y, z 
    real(wp)                                     :: xi, eta, zeta, factor 

    !-- set the order of the x, y, and z Taylor expansions of sin
    order = 1
    x_min = -1000
    !   Parallel determination of the max and min values

    ! setup matrix to store local max and min values
    if(myprocid == 0 )  then
      allocate(xyz_minmax(3,2,nprocs)); xyz_minmax = 0.0_wp
      !-- set the values on the root process
      xyz_minmax(1,1,1) = minval(xg(1,:,:)); xyz_minmax(2,1,1) = minval(xg(2,:,:)); xyz_minmax(3,1,1) = minval(xg(3,:,:))
      xyz_minmax(1,2,1) = maxval(xg(1,:,:)); xyz_minmax(2,2,1) = maxval(xg(2,:,:)); xyz_minmax(3,2,1) = maxval(xg(3,:,:))  
    endif
     
    !-- send the max min values from each process to the root process
    if(myprocid /= 0 ) then
      xyz_loc_minmax(1,1) = minval(xg(1,:,:)); xyz_loc_minmax(2,1) = minval(xg(2,:,:)); xyz_loc_minmax(3,1) = minval(xg(3,:,:))
      xyz_loc_minmax(1,2) = maxval(xg(1,:,:)); xyz_loc_minmax(2,2) = maxval(xg(2,:,:)); xyz_loc_minmax(3,2) = maxval(xg(3,:,:))

      s_tag = 100 + myprocid
      m_size = 3*2
      call mpi_isend(xyz_loc_minmax,m_size,mpi_double,0,s_tag,petsc_comm_world, &
        & s_request_err,i_err)

      call mpi_wait(s_request_err,s_status,i_err)
    
    else
    !-- root process
      do iproc = 1, nprocs-1
        r_tag = 100 + iproc
        m_size = 3*2
        call mpi_irecv(xyz_minmax(1:3,1:2,iproc+1),m_size,mpi_double,iproc,r_tag, &
          & petsc_comm_world,r_request_err,i_err)

        call mpi_wait(r_request_err,r_status,i_err)
      enddo

      !-- determine global min and max x, y, and z values
      x_min = minval(xyz_minmax(1,1,1:nprocs))
      y_min = minval(xyz_minmax(2,1,1:nprocs))
      z_min = minval(xyz_minmax(3,1,1:nprocs))

      x_max = maxval(xyz_minmax(1,2,1:nprocs))
      y_max = maxval(xyz_minmax(2,2,1:nprocs))
      z_max = maxval(xyz_minmax(3,2,1:nprocs))

    endif
    call mpi_bcast(x_min,1,mpi_double,0,PETSC_COMM_WORLD,i_err)
    call mpi_bcast(y_min,1,mpi_double,0,PETSC_COMM_WORLD,i_err)
    call mpi_bcast(z_min,1,mpi_double,0,PETSC_COMM_WORLD,i_err)
    call mpi_bcast(x_max,1,mpi_double,0,PETSC_COMM_WORLD,i_err)
    call mpi_bcast(y_max,1,mpi_double,0,PETSC_COMM_WORLD,i_err)
    call mpi_bcast(z_max,1,mpi_double,0,PETSC_COMM_WORLD,i_err)
    
    !-- apply the transformation on each local element
    do ielem = ihelems(1), ihelems(2)
      !-- obtain element properties
      call element_properties(ielem,n_pts_3d=n_pts_3d)

      do inode = 1,n_pts_3d
        !-- obtain old values (x,y,z)
        x_old = xg(1,inode,ielem)
        y_old = xg(2,inode,ielem)
        z_old = xg(3,inode,ielem)

        !-- compute new values of (x,y,z); here xi, eta, and zeta are variables used as parameters to vary between 0 and 1
        xi = (x_old-x_min)/(x_max-x_min)
        eta = (y_old-y_min)/(y_max-y_min)
        zeta = (z_old-z_min)/(z_max-z_min)

        !x = x_old+1.0_wp
        !y = y_old+1.0_wp
        !z = z_old+1.0_wp

        !x = x_old+abs(x_max-x_min)/10.0_wp*sin(pi*xi)*sin(pi*eta)*sin(pi*zeta)
        !y = y_old+abs(y_max-y_min)/10.0_wp*sin(pi*xi)*sin(pi*eta)*sin(pi*zeta)
        !z = z_old+abs(z_max-z_min)/10.0_wp*sin(pi*xi)*sin(pi*eta)*sin(pi*zeta)
        if(order.eq.1)then
          factor = pi*xi*pi*eta*pi*zeta
        elseif(order.eq.3)then
          factor = (pi*xi-(pi*xi)**3/6.0_wp)*(pi*eta-(pi*eta)**3/6.0_wp)*(pi*zeta-(pi*zeta)**3/6.0_wp)
        elseif(order.eq.5)then
          factor = (pi*xi-(pi*xi)**3/6.0_wp+(pi*xi)**5/120.0_wp)*(pi*eta-(pi*eta)**3/6.0_wp+&
                   (pi*eta)**5/120.0_wp)*(pi*zeta-(pi*zeta)**3/6.0_wp+(pi*zeta)**5/120.0_wp)
        endif

        x = x_old+abs(x_max-x_min)/10.0_wp*factor 
        y = y_old+abs(y_max-y_min)/10.0_wp*factor
        z = z_old+abs(z_max-z_min)/10.0_wp*factor
          
        if((abs(x-x_max)>(x_max-x_min)).or.(abs(y-y_max)>(y_max-y_min))&
           .or.(abs(z-z_max)>(z_max-z_min)))then
           write(*,*)'inode = ',inode
           write(*,*)'xi = ',xi,'eta = ',eta,'zeta = ',zeta
           write(*,*)'x_old = ',x_old,'y_old = ',y_old,'z_old=',z_old
           write(*,*)'x = ',x,'y = ',y,'z=',z
        endif
        !-- overwrite old values with new values
        xg(1,inode,ielem) = x; xg(2,inode,ielem) = y; xg(3,inode,ielem) = z
      enddo
    enddo
    
   !write(*,*)'YOU ARE NOT TRANSFORMING THE GUASS SHELL NODES'
    if(.true.)then
    !-- apply the transformation to the Gauss shell nodes on each element localy
    n_Gau_2d_max = (npoly_max+1)**2
    !-- loop over elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,n_pts_1d=n_LGL_1d_On)


      !-- loop over faces
      do iface = 1,nfacesperelem

        n_LGL_1d_Off  = ef2e(4,iface,ielem)
        n_Gau_1d_Mort = max(n_LGL_1d_On, n_LGL_1d_Off)
        n_Gau_2d_Mort = n_Gau_1d_Mort**2

        !-- loop over the face nodes
        do i_Gau = 1,n_Gau_2d_max  !n_Gau_2d_Mort
          ishift = (iface-1) * n_Gau_2d_max + i_Gau

          x_old = xg_Gau_shell(1,ishift,ielem)
          y_old = xg_Gau_shell(2,ishift,ielem)
          z_old = xg_Gau_shell(3,ishift,ielem)

          xi = (x_old-x_min)/(x_max-x_min)
          eta = (y_old-y_min)/(y_max-y_min)
          zeta = (z_old-z_min)/(z_max-z_min)

          !x = x_old+1.0_wp
          !y = y_old+1.0_wp
          !z = z_old+1.0_wp

          !x = x_old+abs(x_max-x_min)/10.0_wp*sin(pi*xi)*sin(pi*eta)*sin(pi*zeta)
          !y = y_old+abs(y_max-y_min)/10.0_wp*sin(pi*xi)*sin(pi*eta)*sin(pi*zeta)
          !z = z_old+abs(z_max-z_min)/10.0_wp*sin(pi*xi)*sin(pi*eta)*sin(pi*zeta)

          if(order.eq.1)then
            factor = pi*xi*pi*eta*pi*zeta
          elseif(order.eq.3)then
            factor = (pi*xi-(pi*xi)**3/6.0_wp)*(pi*eta-(pi*eta)**3/6.0_wp)*(pi*zeta-(pi*zeta)**3/6.0_wp)
          elseif(order.eq.5)then
            factor = (pi*xi-(pi*xi)**3/6.0_wp+(pi*xi)**5/120.0_wp)*(pi*eta-(pi*eta)**3/6.0_wp+&
                     (pi*eta)**5/120.0_wp)*(pi*zeta-(pi*zeta)**3/6.0_wp+(pi*zeta)**5/120.0_wp)
          endif

          x = x_old+abs(x_max-x_min)/10.0_wp*factor 
          y = y_old+abs(y_max-y_min)/10.0_wp*factor
          z = z_old+abs(z_max-z_min)/10.0_wp*factor
 
          !-- overwrite old values with new values
          xg_Gau_shell(1,ishift,ielem) = x; xg_Gau_shell(2,ishift,ielem) = y; xg_Gau_shell(3,ishift,ielem) = z
        enddo
      enddo
    enddo
    endif

    end subroutine transform_grid

  subroutine snap_to_sphere(ielem,points,origin,xyz)
  !================================================================================================
  !
  ! Purpose: takes the grid and subjects it to a transformation
  !
  ! Assumption: that the element being snaped comes from a soccer ball decomposition of the sphere
  !             and therefore, each surface on the sphere has sides of the same length
  !
  !================================================================================================
 
    !-- use statments
    use initcollocation, only                    : element_properties
    use precision_vars, only                     : pi
 
    implicit none

    !-- input variables
    integer, intent(in)                          :: ielem
    real(wp), intent(in), dimension(8,3)         :: points
    real(wp), intent(in), dimension(3)           :: origin
    real(wp), intent(inout),allocatable, dimension(:,:) :: xyz

    !-- local variables
    integer                                      :: i, j, count_min, count_max, inode
    integer                                      :: n_pts_1d, n_pts_2d, n_pts_3d
    integer, dimension(4)                        :: order
    integer                                      :: i_err
    real(wp),dimension(4,3)                      :: points_plane_r_min, points_plane_r_max
    real(wp),dimension(4,3)                      :: points_plane_r_min_cc, points_plane_r_max_cc
    real(wp), dimension(8)                       :: r_points, theta_points, phi_points
    real(wp)                                     :: r_max, r_min
    real(wp), allocatable, dimension(:)          :: r
    real(wp), allocatable, dimension(:)           :: x_pts_1d

    real(wp), parameter                          :: tol = 10**(-12)

    !-- determine the radius 
    do i = 1,8
      r_points(i) = sqrt( (points(i,1)-origin(1))**2 + (points(i,2)-origin(2))**2 + (points(i,3)-origin(3)) )
    enddo

    !-- determine the min and max radius
    r_min = minval(r_points); r_max = maxval(r_points)

    !-- determine theta and phi
    do i = 1,8
      theta_points(i) = acos( (points(i,3)-origin(3))/r_points(i) )

      !-- logic to bar against negative theta and theta greater than pi
      if( theta_points(i)<0.0_wp) then 
        theta_points(i) = abs(theta_points(i))
      endif
      if( theta_points(i)> pi ) then
        theta_points(i) = 2.0_wp*pi-theta_points(i)
      endif

      phi_points(i) = atan( (points(i,2)-origin(2))/(points(i,1)-origin(1)) )

      !-- logic to bar against negative phi
      if ( phi_points(i)< 0.0_wp) then
        phi_points(i) = 2.0_wp*pi+phi_points(i)
      endif

    enddo

    !-- determine which four points have the same radius
    count_min = 1
    count_max = 1
    do i = 1,8
      if (abs(r_points(1)-r_min).LE.tol) then
        points_plane_r_min(count_min,1) = r_min
        points_plane_r_min(count_min,2) = theta_points(i)
        points_plane_r_min(count_min,3) = phi_points(i) 
        count_min = count_min+1
      elseif(abs(r_points(1)-r_min).LE.tol) then
        points_plane_r_max(count_max,1) = r_max
        points_plane_r_max(count_max,2) = theta_points(i)
        points_plane_r_max(count_max,3) = phi_points(i) 
        count_max = count_max+1
      else
        write(*,*)' the 8 points furnished to initgrid: snap_to_sphere are incorrect',points
        call PetscFinalize(i_err); stop
      endif
    enddo
 
    !-- organize the points in each plane in a counter clockwise fashion where the first point has
    !   the minimum theta and phi values and the last point has the maximum theta and phi values

    !-- plane r_min
    do i = 1,4
      if((abs(points_plane_r_min(i,2)-minval(points_plane_r_min(1:4,2))) < tol)&
         .AND.(abs(points_plane_r_min(i,3)-minval(points_plane_r_min(1:4,3))) < tol))then
        !-- min theta, min phi point
        order(1) = i
      endif 
      if( (abs(points_plane_r_min(i,2)-maxval(points_plane_r_min(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_min(i,3)-minval(points_plane_r_min(1:4,3))) < tol)      )then
        !-- max theta, min phi point
        order(2) = i
      endif    
      if( (abs(points_plane_r_min(i,2)-maxval(points_plane_r_min(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_min(i,3)-maxval(points_plane_r_min(1:4,3))) < tol)      )then
        !-- max theta, max phi point
        order(3) = i
      endif 
      if( (abs(points_plane_r_min(i,2)-minval(points_plane_r_min(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_min(i,3)-maxval(points_plane_r_min(1:4,3))) < tol)      )then
        !-- min theta, max phi point
        order(4) = i
      endif 
    enddo

    !-- reorder
    do i = 1,4
      points_plane_r_min_cc(i,1:3) = points_plane_r_min(order(i),1:3)
    enddo

    !-- plane r_max
    do i = 1,4
      if( (abs(points_plane_r_max(i,2)-minval(points_plane_r_max(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_max(i,3)-minval(points_plane_r_max(1:4,3))) < tol)      )then
        !-- min theta, min phi point
        order(1) = i
      endif 
      if( (abs(points_plane_r_max(i,2)-maxval(points_plane_r_max(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_max(i,3)-minval(points_plane_r_max(1:4,3))) < tol)      )then
        !-- max theta, min phi point
        order(2) = i
      endif    
      if( (abs(points_plane_r_max(i,2)-maxval(points_plane_r_max(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_max(i,3)-maxval(points_plane_r_max(1:4,3))) < tol)      )then
        !-- max theta, max phi point
        order(3) = i
      endif 
      if( (abs(points_plane_r_max(i,2)-minval(points_plane_r_max(1:4,2))) < tol).AND.&
          &(abs(points_plane_r_max(i,3)-maxval(points_plane_r_max(1:4,3))) < tol)      )then
        !-- min theta, max phi point
        order(4) = i
      endif 
    enddo

    !-- reorder
    do i = 1,4
      points_plane_r_max_cc(i,1:3) = points_plane_r_max(order(i),1:3)
    enddo

    !-- obtain element properties
    call element_properties(ielem,n_pts_1d=n_pts_1d,n_pts_2d=n_pts_2d,&
                            n_pts_3d=n_pts_3d,x_pts_1d=x_pts_1d)

    !-- construct nodal distribution in the r direction
    if(allocated(r)) deallocate(r); allocate(r(n_pts_1d)); r = 0.0_wp
    do inode = 1, n_pts_1d
      r(i) = (r_max-r_min)/2.0_wp*x_pts_1d(i)+(r_max+r_min)/2.0_wp
    enddo
    
    !-- construct each plane of points
    allocate(xyz(n_pts_3d,3)); xyz = 0.0_wp
    
    do i = 1, n_pts_1d
      !-- set the new r value 
      do j = 1, 4
        points_plane_r_min_cc(j,1) = r(j)
      enddo
      call snap_to_sphere_patch(ielem,points_plane_r_min_cc,origin,n_pts_1d,xyz(n_pts_2d*(i-1)+1:i*n_pts_2d,1:3))
    enddo

    !-- deallocate statments
    deallocate(r)
  end subroutine snap_to_sphere

  subroutine snap_to_sphere_patch(ielem,points,origin,n_pts_1d,xyz)
  !================================================================================================
  !
  ! Purpose: takes the grid and subjects it to a transformation
  !
  !================================================================================================
 
    !-- use statments
    use initcollocation, only                    : element_properties
  
    implicit none

    !-- input variables
    integer, intent(in)                          :: ielem, n_pts_1d
    real(wp), intent(in), dimension(4,3)         :: points
    real(wp), intent(in), dimension(3)           :: origin
    real(wp), intent(inout),dimension(n_pts_1d**2,3) :: xyz

    !-- local variables
    integer                                      :: inode
    integer                                      :: i, j
    integer                                      :: i_err
    real(wp), allocatable, dimension(:)          :: x_pts_1d
    real(wp), dimension(4)                       :: r, theta_points, phi_points
    real(wp)                                     :: theta, phi, theta_min, theta_max, phi_min, phi_max
    real(wp)                                     :: xi, eta

    real(wp), parameter                          :: tol = 10**(-12)

    !-- determine the radius and ensure that all 4 points have the same radius 
    do i = 1,4
      r(i) = sqrt( (points(i,1)-origin(1))**2 + (points(i,2)-origin(2))**2 + (points(i,3)-origin(3)) )
    enddo

    if ( (abs(r(1)-r(2)).GE.tol).OR.(abs(r(1)-r(3)).GE.tol).OR.(abs(r(1)-r(4)).GE.tol) ) then
      write(*,*)'the four points that have been provided to initgrid: snap_to_sphere_patch do not have the same radius'
      call PetscFinalize(i_err); stop
    endif

    !-- determine min/max theta, phi
    do i = 1,4
      theta_points(i) = acos( (points(i,3)-origin(3))/r(1) )

      !-- logic to bar against negative theta and theta greater than pi
      if( theta_points(i).LE.0.0_wp) then 
        theta_points(i) = abs(theta_points(i))
      endif
      if( theta_points(i)> pi ) then
        theta_points(i) = 2.0_wp*pi-theta_points(i)
      endif

      phi_points(i) = atan( (points(i,2)-origin(2))/(points(i,1)-origin(1)) )

      !-- logic to bar against negative phi
      if ( phi_points(i)< 0.0_wp) then
        phi_points(i) = 2.0_wp*pi+phi_points(i)
      endif

    enddo

    theta_min = minval(theta_points); theta_max = maxval(theta_points)
    phi_min = minval(phi_points);     phi_max = maxval(phi_points)
                
    !-- obtain element properties
    call element_properties(ielem,x_pts_1d=x_pts_1d)

    inode = 1
    do i = 1,n_pts_1d
      do j = 1,n_pts_1d
        xi = x_pts_1d(i)
        eta = x_pts_1d(j)

        theta = (theta_max-theta_min)/2.0_wp*xi+(theta_max+theta_min)/2.0_wp
        phi = (phi_max-phi_min)/2.0_wp*eta+(phi_max+phi_min)/2.0_wp

        xyz(inode,1) = origin(1)+r(1)*cos(phi)*sin(theta)
        xyz(inode,2) = origin(2)+r(1)*sin(phi)*cos(theta)
        xyz(inode,3) = origin(3)+r(1)*cos(theta)
  
        inode = inode+1
      enddo
    enddo


    !-- deallocate statments
    deallocate(x_pts_1d)

  end subroutine snap_to_sphere_patch


subroutine write_matrix_to_file_matlab(A,n,m,file_name)
!==================================================================================================
!
! Purpose: writes a matrix A to file in Matlab format
!
! Comments:
!
! Additional documentation: THIS IS NOT SETUP FOR PARALLEL RUNS ONLY USE ONE PROCESS
!
! Unit tests: 
!
! Inputs: A (size(A) = [n,m], real): matrix to be writen to file
!         n (integer): number of rows
!         m (integer): number of rows
!         file_name (char(len=*): string with the name of the file
! 
!
! outputs: 
!
!
!==================================================================================================
  
  !-- variables
  use precision_vars, only                     : get_unit
  use referencevariables, only                 : myprocid

  implicit none
                     
  !-- input variables
  integer, intent(in)                            :: n, m
  real(wp),intent(in),dimension(n,m)             :: A
  character(len=*),optional,intent(in)           :: file_name 

  !-- local variables
  integer                                        :: iunit
  integer                                        :: i,j

  !-- file access variables
  integer                                        :: ios

  !-- writing to file
  print *,"writing to: ",file_name

  !-- open file on the master proc
  if (myprocid==0) then
    !-- obtain a free Fortarn unit
    call get_unit(iunit)
    
    !-- open the file
    open(UNIT=iunit,FILE=file_name,STATUS='NEW',FORM='FORMATTED',IOSTAT=ios)
    if(ios.NE.0) then 
      write(*,*)"File = ",file_name," not opened correctly in initgrd: write_matrix_to_file_matlab(), iostat = ",ios
      stop
    endif
    
    !-- write the root process nodal information to file
    write(iunit,*)'A = [...'
    do i = 1,n
        write(iunit,*)(A(i,j),j=1,m),";"   
    enddo
    write(iunit,*)"];"

    !-- close file
    close(UNIT=iunit)
  else
    write(*,*)'trying to write a matrix to file using initgrid:write_matrix_to_file_matlab, during a parallel run'
  endif
 
end subroutine write_matrix_to_file_matlab
  subroutine Amat2(x, N, pmat,pinv,qmat,dmat)                                
    ! Form the differentiation matrix dmat                                  
    ! This matrix is used to compute the derivative in the Gauss-Labotto points   
    !NOTE: HAD TO ADD THIS BECAUSE IN modify_metrics_nonconforming ONE OF THE VARIABLES IS Amat

    implicit none

    integer,                  intent(in)    :: N
    real(wp), dimension(N),   intent(in)    :: x
    real(wp), dimension(N),   intent(inout) :: pmat,pinv
    real(wp), dimension(N,N), intent(inout) :: qmat,dmat

    real(wp), allocatable, dimension(:)     :: wk, f, dex

    real(wp)                                :: err, errmax
    real(wp), parameter                     :: half = 0.50000000000000000000000000_wp
    logical                                 :: diagnostic = .false.

    integer                                 :: i,j,k


    ! Allocate memory
    allocate(wk(N))
    allocate(f(N))
    allocate(dex(N))

    do k = 1,N 
      wk(k) = 1.0_wp
      do i = 1,N 
        if(i/=k) wk(k) = wk(k) * (x(k) - x(i)) 
      enddo 
    enddo 

    do i = 1,N 
      do j = 1,N 
        if(j/=i) dmat(i,j)=wk(i)/ (wk(j) * (x(i) - x(j)) ) 
      enddo 
    enddo 
    do j = 1,N 
      dmat(j,j)=0.0_wp
      do k = 1,N 
        if(k/=j) dmat(j,j) = dmat(j,j) + 1.0_wp/(x(j) - x(k)) 
      enddo 
    enddo 

    pinv(1) = (-dmat(1,1)+dmat(N,N))
    pinv(N) = pinv(1)
    do i = 2,(N+1)/2
      pinv(i) = - (dmat(i,1)*pinv(1)) / dmat(1,i) 
      pinv(N+1-i) = pinv(i)
    enddo

    pmat(:) = 1.0_wp / pinv(:)

    do i = 1,N
      do j = 1,N
        qmat(i,j) = dmat(i,j)*pmat(i)
      enddo
    enddo

    errmax = (qmat(1,1) + half)**2 + (qmat(N,N) - half)**2
    do i = 1,N
      do j = 1,N
        if(((i /= 1) .or. (j /= 1)) .and.  ((i /= N) .or. (j /= N)))then 
          errmax = errmax + (qmat(i,j) + qmat(j,i))**2
        endif
      enddo
    enddo
    errmax = sqrt(errmax/N/N)

    !       differentiate exact data to test accuracy of discrete operator

    do i=1,N
      f(i)    =         x(i)**(N-1)
      dex(i)  = (N-1) * x(i)**(N-2)
    enddo 
    do i = 1,N 
      wk(i)=0.0_wp
      do j = 1,N 
        wk(i) = wk(i) + dmat(i,j) * f(j) 
      enddo 
    enddo 

    do i = 1, N
      err=abs(wk(i)-dex(i)) 
      if(err.gt.errmax)errmax=err 
    enddo 

    if(diagnostic) then
      do i = 1,N
         write(*,287)i,pmat(i)
      enddo
      do i = 1,N
         do j = i,N
           write(*,288)i,j,qmat(i,j)
         enddo
      enddo
      stop
    endif
 287  format('      d',i2,' = ',f22.17)
 288  format('      q',i2,1x,i2,' = ',f22.17)
    if(errmax > 1.0e-13_wp) then
      write(*,*)'x',x
      write(*,*)'dmat'
      do i = 1,N
        write(*,*)(dmat(i,j),j=1,N)
      enddo
      write(*,*)'roundoff alert',errmax 
    endif

    deallocate(f  )
    deallocate(dex)
    deallocate(wk )

    !   150 format( 25(f5.2,1x)) 
    !   160 format( 25(f6.3,1x))

    return 
  end subroutine Amat2

  subroutine e_edge2e_connectivity() ! Serial subroutine
  !================================================================================================
  !
  ! Purpose: This constructs the master e_edge2e 
  !  e_edge2e(1,iedge,ipartner,iface,ielem) = number of nodes in one dimension of the ith element (ipartner) touching iedge
  !  e_edge2e(2,iedge,ipartner,iface,ielem) = element number of the ith element (ipartner) touching iedge
  !  e_edge2e(3,iedge,ipartner,iface,ielem) = procid of the ith element (ipartner) touching iedge (not setup)
  !
  ! Notes:
  !-- the edge numbering is counter-clockwise starting from the origin i.e.
  !                 edge 3
  !               ----------
  !               |        |
  !        edge4  |        | edge 2
  !               |        |
  !               |        |
  !               ----------
  !               edge 1
  !
  !================================================================================================

    use mpimod
    use variables, only            : e_edge2e, vx_Master, ic2nh, iae2e2e, jae2e2e
    use referencevariables, only   : nelems, ndim, number_of_possible_partners
    use collocationvariables, only : elem_props
    use precision_vars, only       : magnitude

    implicit none
 
    integer                :: number_of_edges_per_face, number_of_faces
    integer                :: ielem, iface, connector

    real(wp), dimension(3) :: v1, v2
    real(wp)               :: tol = 1.0e-12_wp
    integer                :: kelem, partners_cnt, max_partners, kvertex, partner
    integer, allocatable, dimension(:,:,:,:,:) ::  e_edge2e_tmp
    logical                :: v1match, v2match
    integer                :: nvertex_per_element
    integer                :: ii

    integer                :: i_err

    !-- now we have all of the required information to build the global e_edge2e
    if(ndim.EQ.1)then
      write(*,*)'Logic in initgrid: e_edge2e_connectivity not setup for ndim = ',ndim
      call PetscFinalize(i_err); stop
    elseif(ndim.EQ.1)then
      number_of_edges_per_face = 1
    elseif(ndim.EQ.3)then
      number_of_edges_per_face = 4
    else
       write(*,*)'Error in initgrid: e_edge2e_connectivity ndim incorrect, ndim = ',ndim
      call PetscFinalize(i_err); stop   
    endif

    number_of_faces = ndim*2
    nvertex_per_element = 2**ndim

    
    max_partners = 0                                    !   approximate size of buckets for e_edge2e_tmp
    do ielem = 1,nelems
       max_partners = max(iae2e2e(ielem+1)-iae2e2e(ielem),max_partners)
    enddo 

    !-- allocate vector to store the matches. Worse case scenario each connector touches every element
    allocate(e_edge2e_tmp(3,number_of_edges_per_face,max_partners,number_of_faces,1:nelems)) 
    e_edge2e_tmp = -1000

    max_partners = 0

    Element_Loop: do ielem = 1,nelems
      Face_Loop: do iface = 1,number_of_faces
        Connector_Loop : do connector = 1,4

          !-- determine the verticies for the connector
          if(iface.EQ.1)then
            if(connector.EQ.1)then
              !-- connector at xi_2 = 0, xi_3 = 0
              v1 = vx_Master(:,ic2nh(1,ielem))
              v2 = vx_Master(:,ic2nh(2,ielem))
            elseif(connector.EQ.2)then
              !-- connector at xi_1 = 1, xi_3 = 0
              v1 = vx_Master(:,ic2nh(2,ielem))
              v2 = vx_Master(:,ic2nh(3,ielem))
            elseif(connector.EQ.3)then
              !-- connector at xi_2 = 1, xi_3 = 0
              v1 = vx_Master(:,ic2nh(3,ielem))
              v2 = vx_Master(:,ic2nh(4,ielem))
            elseif(connector.EQ.4)then
              !-- connector at xi_1 = 0, xi_3 = 0
              v1 = vx_Master(:,ic2nh(4,ielem))
              v2 = vx_Master(:,ic2nh(1,ielem))
            endif
          elseif(iface.EQ.2)then
            if(connector.EQ.1)then
              !-- connector at xi_2 = 0, xi_3 = 0
              v1 = vx_Master(:,ic2nh(1,ielem))
              v2 = vx_Master(:,ic2nh(2,ielem))
            elseif(connector.EQ.2)then
              !-- connector at xi_1 = 1, xi_2 = 0
              v1 = vx_Master(:,ic2nh(2,ielem))
              v2 = vx_Master(:,ic2nh(6,ielem))
            elseif(connector.EQ.3)then
              !-- connector at xi_2 = 0, xi_3 = 1
              v1 = vx_Master(:,ic2nh(6,ielem))
              v2 = vx_Master(:,ic2nh(5,ielem))
            elseif(connector.EQ.4)then
              !-- connector at xi_1 = 0, xi_2 = 0
              v1 = vx_Master(:,ic2nh(5,ielem))
              v2 = vx_Master(:,ic2nh(1,ielem))
            endif
          elseif(iface.EQ.3)then
            if(connector.EQ.1)then
              !-- connector at xi_1 = 1, xi_3 = 0
              v1 = vx_Master(:,ic2nh(2,ielem))
              v2 = vx_Master(:,ic2nh(3,ielem))
            elseif(connector.EQ.2)then
              !-- connector at xi_1 = 1, xi_2 = 1
              v1 = vx_Master(:,ic2nh(3,ielem))
              v2 = vx_Master(:,ic2nh(7,ielem))
            elseif(connector.EQ.3)then
              !-- connector at xi_1 = 1, xi_3 = 1
              v1 = vx_Master(:,ic2nh(7,ielem))
              v2 = vx_Master(:,ic2nh(6,ielem))
            elseif(connector.EQ.4)then
              !-- connector at xi_1 = 1, xi_2 = 0
              v1 = vx_Master(:,ic2nh(6,ielem))
              v2 = vx_Master(:,ic2nh(2,ielem))
            endif
          elseif(iface.EQ.4)then
            if(connector.EQ.1)then
              !-- connector at xi_2 = 1, xi_3 = 0
              v1 = vx_Master(:,ic2nh(4,ielem))
              v2 = vx_Master(:,ic2nh(3,ielem))
            elseif(connector.EQ.2)then
              !-- connector at xi_1 = 1, xi_2 = 1
              v1 = vx_Master(:,ic2nh(3,ielem))
              v2 = vx_Master(:,ic2nh(7,ielem))
            elseif(connector.EQ.3)then
              !-- connector at xi_2 = 1, xi_3 = 1
              v1 = vx_Master(:,ic2nh(7,ielem))
              v2 = vx_Master(:,ic2nh(8,ielem))
            elseif(connector.EQ.4)then
              !-- connector at xi_1 = 0, xi_2 = 1
              v1 = vx_Master(:,ic2nh(8,ielem))
              v2 = vx_Master(:,ic2nh(4,ielem))
            endif
          elseif(iface.EQ.5)then
            if(connector.EQ.1)then
              !-- connector at xi_1 = 0, xi_3 = 0
              v1 = vx_Master(:,ic2nh(1,ielem))
              v2 = vx_Master(:,ic2nh(4,ielem))
            elseif(connector.EQ.2)then
              !-- connector at xi_1 = 0, xi_2 = 1
              v1 = vx_Master(:,ic2nh(4,ielem))
              v2 = vx_Master(:,ic2nh(8,ielem))
            elseif(connector.EQ.3)then
              !-- connector at xi_1 = 0, xi_3 = 1
              v1 = vx_Master(:,ic2nh(8,ielem))
              v2 = vx_Master(:,ic2nh(5,ielem))
            elseif(connector.EQ.4)then
              !-- connector at xi_1 = 0, xi_2 = 0
              v1 = vx_Master(:,ic2nh(5,ielem))
              v2 = vx_Master(:,ic2nh(1,ielem))
            endif
          elseif(iface.EQ.6)then
            if(connector.EQ.1)then
              !-- connector at xi_2 = 0, xi_3 = 1
              v1 = vx_Master(:,ic2nh(5,ielem))
              v2 = vx_Master(:,ic2nh(6,ielem))
            elseif(connector.EQ.2)then
              !-- connector at xi_1 = 1, xi_3 = 1
              v1 = vx_Master(:,ic2nh(6,ielem))
              v2 = vx_Master(:,ic2nh(7,ielem))
            elseif(connector.EQ.3)then
              !-- connector at xi_2 = 1, xi_3 = 1
              v1 = vx_Master(:,ic2nh(7,ielem))
              v2 = vx_Master(:,ic2nh(8,ielem))
            elseif(connector.EQ.4)then
              !-- connector at xi_1 = 0, xi_3 = 1
              v1 = vx_Master(:,ic2nh(8,ielem))
              v2 = vx_Master(:,ic2nh(5,ielem))
            endif
          endif

          !-- counter for the number of elements touching a given connector
          partners_cnt = 0
!         !-- brute force check which elements match
!         do kelem = 1,nelems
          do ii = iae2e2e(ielem),iae2e2e(ielem+1)-1
            kelem = jae2e2e(ii)
            v1match = .false.
            v2match = .false.
            if(kelem.NE.ielem)then
              do kvertex = 1,8
                if(magnitude( v1-vx_Master(:,ic2nh(kvertex,kelem)) ) < tol)then
                  v1match = .true.
                endif
                if(magnitude( v2-vx_Master(:,ic2nh(kvertex,kelem)) )  < tol)then
                  v2match = .true.
                endif
              enddo
            endif
            if(v1match.AND.v2match)then
                 partners_cnt = partners_cnt+1
                e_edge2e_tmp(1,connector,partners_cnt,iface,ielem) = elem_props(2,kelem)
                e_edge2e_tmp(2,connector,partners_cnt,iface,ielem) = kelem
            endif
          enddo

          !-- keep track of the maximum number of neighbours (will be used to set number_of_possible_partners
          max_partners = max(partners_cnt,max_partners)

        enddo Connector_Loop
      enddo Face_Loop
    enddo Element_Loop

    number_of_possible_partners = max_partners

    allocate(e_edge2e(3,number_of_edges_per_face,number_of_possible_partners,number_of_faces,1:nelems))
    e_edge2e = -1000

    do ielem = 1,nelems
      do iface = 1,number_of_faces
        do connector = 1,4
          do partner = 1,max_partners
            e_edge2e(1:2,connector,partner,iface,ielem) = e_edge2e_tmp(1:2,connector,partner,iface,ielem)
          enddo
!          write(*,*)'=========================================================================&
!                                ================================'
!          write(*,*)'ielem = ', ielem, 'iface = ',iface,' connector = ',connector
!          write(*,*)'e_edge2e_tmp(1,connector,1:max_partners,iface,ielem) = ',e_edge2e_tmp(1,connector,1:max_partners,iface,ielem)
!          write(*,*)'e_edge2e_tmp(2,connector,1:max_partners,iface,ielem) = ',e_edge2e_tmp(2,connector,1:max_partners,iface,ielem)
!          write(*,*)'=========================================================================&
!                                ================================'
         enddo
      enddo
    enddo

  end subroutine e_edge2e_connectivity

end module initgrid

