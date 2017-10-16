module initgrid
  use precision_vars
  use iso_c_binding
  implicit none

  ! #include "finclude/petscsys.h"
  ! #include "finclude/petscvec.h"
  ! ! #include "finclude/petscdmda.h"
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
  public modify_metrics_nonconforming

  integer, allocatable, dimension(:,:), target :: edge_2_faces
  integer, allocatable, dimension(:), target :: edge_2_facedirections
  integer, parameter :: edge_2_nfacesperelem = 2
  integer, parameter :: edge_2_nverticesperface = 1
  integer, allocatable, dimension(:,:), target :: quad_4_faces
  integer, allocatable, dimension(:), target :: quad_4_facedirections
  integer, parameter :: quad_4_nfacesperelem = 4
  integer, parameter :: quad_4_nverticesperface = 2
  integer, allocatable, dimension(:,:), target :: hex_8_faces
  integer, allocatable, dimension(:), target :: hex_8_facedirections
  integer, parameter :: hex_8_nfacesperelem = 6
  integer, parameter :: hex_8_nverticesperface = 4
  integer, pointer, dimension(:,:) :: eltypfaces
  integer, pointer, dimension(:) :: elfacedirections

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

!    write(*,*) iae2v(595)

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

!        write(*,*) 'After if', xadjtmp(594)
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
    
    integer :: iell, ielh

    real(wp) :: diff_x, diff_y, diff_z
    real(wp), parameter :: toll = 1e-6
    real :: tmp
    
    continue
  
    iell = ihelems(1)
    ielh = ihelems(2)

    call srand(seed)

    do ielem = iell, ielh
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
    use controlvariables, only: Grid_Topology, cylinder_x0, cylinder_x1
    use referencevariables
    use variables, only: xg, vx, e2v, ef2e
    use collocationvariables, only: n_LGL_1d_p1, elem_props
    use initcollocation, only: JacobiP11
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
    real(wp), dimension(nodesperedge)  :: xi
    real(wp), dimension(3)  :: dx
    real(wp), dimension(3)  :: x00,x01

    ! number of nodes in each element
    nodesperelem_max = (npoly_max+1)**ndim

    ! allocate global node matrix
    allocate(xg(3,1:nodesperelem_max,ihelems(1):ihelems(2)))
    xg = 0.0_wp

    ! loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)

      ! nE is size of edge on element (varies with element)
      nE = elem_props(2,ielem)

      ! allocate local nodes

      if(allocated(xl)) deallocate(xl) ; allocate(xl(3,1:nE,1:nE,1:nE)) ;  xl = 0.0_wp

      if(allocated(x_LGL_1d)) deallocate(x_LGL_1d) ; allocate(x_LGL_1d(1:nE)) 

      call JacobiP11(nE-1,x_LGL_1D)

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

      end select

      ! build faces
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
      ! build volumes
      if (ndim > 2) then
        call TFI3D(xl(:,:,:,:),nE,x_LGL_1d)
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
    end do
    deallocate(xl)

  end subroutine calcnodes_LGL

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
                       , kfacenodes_Gau_p2, ifacenodes_Gau_p2 &
                       , kfacenodes, ifacenodes

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
      & periodic_elem_face_ids_x2, periodic_elem_face_ids_x3, &
      & kfacenodes_LGL_p0,ifacenodes_LGL_p0,                  &
      & kfacenodes_LGL_p1,ifacenodes_LGL_p1
    use collocationvariables, only: elem_props, n_LGL_1d_p0, n_LGL_1d_p1

    ! Nothing is implicitly defined
    implicit none

    integer, allocatable, dimension(:,:) :: kfacenodes
    integer, allocatable, dimension(:)   :: ifacenodes

    integer ::  ielem, inode, jnode, iface, knode
    integer ::  i_low

    integer :: iell, ielh
    real(wp) :: x1(3), x2(3)
    real(wp), parameter :: nodetol = 1.0e-8_wp

    integer :: i_p_face, p_dir, cnt_coord, i_coord
    logical :: match_found
    real(wp), dimension(2) :: x1_p, x2_p

    integer :: cnt_debug
    integer :: n_LGL_1d, n_LGL_2d, nodesperface_max


    continue

    cnt_debug = 0

    ! Low and High volumetric element index
    iell = ihelems(1) ; ielh = ihelems(2) ;

    ! efn2efn contains the partner node information of every facenode in the domain

    nodesperface_max = (npoly_max+1)**(ndim-1)

    allocate(efn2efn(4,nfacesperelem*nodesperface_max,iell:ielh))
    efn2efn = -1000

    ! Initialize position of the ghost point in the stack
    i_low = 0

    ! Loop over elements
    do ielem = iell, ielh
      
      n_LGL_1d = elem_props(2,ielem)**1 
      n_LGL_2d = elem_props(2,ielem)**2 

      if(allocated(kfacenodes)) deallocate(kfacenodes) ; allocate(kfacenodes(1:n_LGL_2d,1:nfacesperelem))
      if(allocated(ifacenodes)) deallocate(ifacenodes) ; allocate(ifacenodes(1:n_LGL_2d*nfacesperelem))

      if(n_LGL_1d == n_LGL_1d_p0) then
        kfacenodes(:,:) = kfacenodes_LGL_p0(:,:)
        ifacenodes(:)   = ifacenodes_LGL_p0(:)
      else
        kfacenodes(:,:) = kfacenodes_LGL_p1(:,:)
        ifacenodes(:)   = ifacenodes_LGL_p1(:)
      endif

      ! Reset facial node index counter
      knode = 0
      
      ! Loop over faces
      do iface = 1, nfacesperelem
        
        ! If on boundary, connect to self
        if (ef2e(1,iface,ielem) < 0) then
          
          ! Loop over nodes on the boundary face
          do inode = 1, n_LGL_2d
            
            ! Update facial node index counter
            knode = knode + 1
            
            ! The first index is the volumetric node index and the second
            ! index is the element index
            efn2efn(:,knode,ielem) = (/ ifacenodes(knode), ielem, 0 /)
          
          end do

        else if ((ef2e(3,iface,ielem) /= myprocid) .and. (ef2e(4,iface,ielem) == elem_props(2,ielem))) then ! A parallel interface

          ! Initialize match_found
          match_found = .false.
          
          ! Loop through the elements that owns a periodic face in the x1
          ! direction
          if (size(periodic_elem_face_ids_x1(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic
            ! face
            do i_p_face = 1, size(periodic_elem_face_ids_x1(1,:))

              if (periodic_elem_face_ids_x1(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x1(2,i_p_face) == iface) then

                ! There is a match: change logical value of match_found
                match_found = .true.

                ! Get the direction of "periodicity"
                p_dir = periodic_elem_face_ids_x1(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_LGL_2d
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg(:,ifacenodes(knode),ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0
                  
                  do i_coord = 1,3
                    
                    if (i_coord /= p_dir) then

                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  do jnode = 1, n_LGL_2d
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xghst_LGL(:,i_low + jnode)
                
                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0

                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,&
                        & ef2e(1,iface,ielem))

                      ! Set the element of the connected node
                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)

                      ! Set the node index in the ghost array
                      efn2efn(3,knode,ielem) = i_low + jnode

!                     efn2efn(4,knode,ielem) = jnode
                      
                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

                ! Update the position in the ghost stack
                i_low = i_low + n_LGL_2d

              end if ! End if match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit 
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if check periodic face in x1 direction


          ! Loop through the elements that owns a periodic face in the x2
          ! direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x2(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic
            ! face
            do i_p_face = 1, size(periodic_elem_face_ids_x2(1,:))

              if (periodic_elem_face_ids_x2(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x2(2,i_p_face) == iface) then

                ! There is a match: change logical value of match_found
                match_found = .true.

                ! Get the direction of "periodicity"
                p_dir = periodic_elem_face_ids_x2(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_LGL_2d
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg(:,ifacenodes(knode),ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0
                  
                  do i_coord = 1,3
                    
                    if (i_coord /= p_dir) then

                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  do jnode = 1, n_LGL_2d
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xghst_LGL(:,i_low + jnode)
                
                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0

                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,&
                        & ef2e(1,iface,ielem))

                      ! Set the element of the connected node
                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)

                      ! Set the node index in the ghost array
                      efn2efn(3,knode,ielem) = i_low + jnode

!                     efn2efn(4,knode,ielem) = jnode
                      
                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

                ! Update the position in the ghost stack
                i_low = i_low + n_LGL_2d

              end if ! End if match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit 
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if check periodic face in x2 direction


          ! Loop through the elements that owns a periodic face in the x3
          ! direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x3(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic
            ! face
            do i_p_face = 1, size(periodic_elem_face_ids_x3(1,:))

              if (periodic_elem_face_ids_x3(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x3(2,i_p_face) == iface) then

                ! There is a match: change logical value of match_found
                match_found = .true.

                ! Get the direction of "periodicity"
                p_dir = periodic_elem_face_ids_x3(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_LGL_2d
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg(:,ifacenodes(knode),ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0
                  
                  do i_coord = 1,3
                    
                    if (i_coord /= p_dir) then

                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  do jnode = 1, n_LGL_2d
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xghst_LGL(:,i_low + jnode)
                
                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0

                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,&
                        & ef2e(1,iface,ielem))

                      ! Set the element of the connected node
                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)

                      ! Set the node index in the ghost array
                      efn2efn(3,knode,ielem) = i_low + jnode

!                     efn2efn(4,knode,ielem) = jnode
                      
                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

                ! Update the position in the ghost stack
                i_low = i_low + n_LGL_2d

              end if ! End if match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit 
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if check periodic face in x3 direction

          if (match_found .eqv. .false.) then

            ! Loop over the nodes on the face
            do inode = 1, n_LGL_2d

              ! Update the facial node index counter
              knode = knode + 1

              ! Save the coordinates of the facial node
              x1 = xg(:,ifacenodes(knode),ielem)
              
              ! Search for the connected node on face of the connected element
              do jnode = 1, n_LGL_2d

                ! Coordinates of the jnode
                ! ef2e(2) gives the element of the neighbor
                x2 = xghst_LGL(:,i_low + jnode)
                
                ! Check the distance between the two nodes
                if (magnitude(x1-x2) <= nodetol) then
                  
                  ! Set the volumetric node index of the connected node
                  efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))
                  
                  ! Set the element of the connected node
                  efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)
                  
                  ! Set the node index in the ghost array
                  efn2efn(3,knode,ielem) = i_low + jnode

!                 efn2efn(4,knode,ielem) = jnode
                  
                  exit
                
                end if
              
              end do
              
              ! Print information at screen if there is a problem and stop
              ! computation
              if (jnode > n_LGL_2d .and. myprocid==1) then
                write(*,*) 'Connectivity error in face-node connectivity.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
                write(*,*) 'Node coordinates and ghost node coordinates'
                write(*,*) x1, xghst_LGL(:,i_low + 1:i_low + n_LGL_2d)
                write(*,*) 'Exiting...'
                stop
              end if

            end do

            ! Update the position in the ghost stack
            i_low = i_low + n_LGL_2d
          
          end if

        else if (ef2e(4,iface,ielem) == elem_props(2,ielem)) then ! Not a parallel interface

          ! Initialize match_found
          match_found = .false.

          if (size(periodic_elem_face_ids_x1(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if the iface is a
            ! periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x1(1,:))

              if (periodic_elem_face_ids_x1(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x1(2,i_p_face) == iface) then

                ! There is a match
                match_found = .true.

                ! Get the direction of periodicity
                p_dir = periodic_elem_face_ids_x1(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_LGL_2d
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg(:,ifacenodes(knode),ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0 

                  do i_coord = 1, 3
                    
                    if (i_coord /= p_dir) then
                      
                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  ! Search for the connected node on the face of the connected 
                  ! element
                  do jnode = 1,n_LGL_2d
                    ! Coordinates of the jnode
                    ! ef2e(1) gives the face on the neighboring element and
                    ! ef2e(2) gives the element
                    x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)), &
                      & ef2e(2,iface,ielem))

                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0
                    
                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check the distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))
                      
                      ! Set the element of the connected node
                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)
                      
                      ! Set the index of the connected node
                      efn2efn(4,knode,ielem) = jnode

                      cnt_debug = cnt_debug + 1

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

              end if ! Match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if periodic x1 direction


          ! If the iface is not a periodic face  in the x1 direction, check
          ! if it is a periodic face in the x2 direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x2(1,:)) /= 0) then
           
            ! Check if the ielem owns a periodic face and if the iface is a
            ! periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x2(1,:))

              if (periodic_elem_face_ids_x2(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x2(2,i_p_face) == iface) then

                ! There is a match
                match_found = .true.

                ! Get the direction of periodicity
                p_dir = periodic_elem_face_ids_x2(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_LGL_2d
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg(:,ifacenodes(knode),ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0 

                  do i_coord = 1, 3
                    
                    if (i_coord /= p_dir) then
                      
                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  ! Search for the connected node on the face of the connected 
                  ! element
                  do jnode = 1,n_LGL_2d
                    ! Coordinates of the jnode
                    ! ef2e(1) gives the face on the neighboring element and
                    ! ef2e(2) gives the element
                    x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)), &
                      & ef2e(2,iface,ielem))

                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0
                    
                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check the distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))
                      
                      ! Set the element of the connected node
                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)
                      
                      ! Set the index of the connected node
                      efn2efn(4,knode,ielem) = jnode

                      cnt_debug = cnt_debug + 1

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

              end if ! Match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit
              end if

            end do ! End do loop over the elements that own a periodic face
          
          end if ! End if periodic x2 direction


          ! If the iface is not a periodic face in the x2 direction, check
          ! if it is a periodic face in the x3 direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x3(1,:)) /= 0) then
           
            ! Check if the ielem owns a periodic face and if the iface is a
            ! periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x3(1,:))

              if (periodic_elem_face_ids_x3(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x3(2,i_p_face) == iface) then

                ! There is a match
                match_found = .true.

                ! Get the direction of periodicity
                p_dir = periodic_elem_face_ids_x3(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_LGL_2d
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg(:,ifacenodes(knode),ielem)
                
                  ! Extract from x1 the two invariant coordinates
                  cnt_coord = 0 

                  do i_coord = 1, 3
                    
                    if (i_coord /= p_dir) then
                      
                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  ! Search for the connected node on the face of the connected 
                  ! element
                  do jnode = 1,n_LGL_2d
                    ! Coordinates of the jnode
                    ! ef2e(1) gives the face on the neighboring element and
                    ! ef2e(2) gives the element
                    x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)), &
                      & ef2e(2,iface,ielem))

                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0
                    
                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check the distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))
                      
                      ! Set the element of the connected node
                      efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)
                      
                      ! Set the index of the connected node
                      efn2efn(4,knode,ielem) = jnode

                      cnt_debug = cnt_debug + 1

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

              end if ! Match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit
              end if

            end do ! End do loop over the elements that own a periodic face
          
          end if ! End if periodic x3 direction

          if (match_found .eqv. .false.) then

            ! Loop over the nodes on the face
            do inode = 1, n_LGL_2d

              ! Update the facial node index counter
              knode = knode + 1

              ! Save coordinates of the facial ndoes
              x1 = xg(:,ifacenodes(knode),ielem)
              ! Search the for connected node on the face of the connected 
              ! element
              
              do jnode = 1, n_LGL_2d
                
                ! Coordinates of the jnode
                ! ef2e(1) gives the face on the neighboring element and
                ! ef2e(2) gives the element
                x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)), &
                  & ef2e(2,iface,ielem))
               
                ! Check the distance between the two nodes
                if (magnitude(x1-x2) <= nodetol) then

                  ! Set the volumetric node index of the connected node
                  efn2efn(1,knode,ielem) = kfacenodes(jnode,ef2e(1,iface,ielem))
                  
                  ! Set the element of the connected node
                  efn2efn(2,knode,ielem) = ef2e(2,iface,ielem)
                  
                  ! Set the index of the connected node
                  efn2efn(4,knode,ielem) = jnode
                  
                  exit
                
                end if
              
              end do ! End do jnode

              ! Print information at screen if there is a problem and stop computation

              if (efn2efn(1,knode,ielem) < 0 .or. efn2efn(2,knode,ielem) < 0) then
                write(*,*) 'Connectivity error in face-node connectivity.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
                write(*,*) 'Node coordinates'
                write(*,*) x1
                write(*,*) 'Possible partner node coordinates'
                
                do jnode = 1, n_LGL_2d
                  x2 = xg(:,kfacenodes(jnode,ef2e(1,iface,ielem)), &
                    & ef2e(2,iface,ielem))
                  write(*,*) x2
                end do 

                write(*,*) 'Exiting...'
                stop
              end if

            end do ! End do inode
          
          end if ! End if not a periodic face (match_found = .false.)
              
        end if ! End if type of face (boundary, off processor or on processor)
      
      end do ! End do loop over faces of the element

    end do ! End do loop elements owned by the processor

    return
  end subroutine calculate_face_node_connectivity_LGL

  !============================================================================
  
  subroutine calcfacenormals_LGL()
    ! this subroutine calculates the outward facing normals
    ! of each facial node
    use referencevariables
    use variables, only: kfacenodes, facenodenormal, r_x, ef2e, efn2efn, Jx_r
    use collocationvariables, only: n_LGL_1d_p1, elem_props

    implicit none

    ! indices
    integer :: ielem, kelem, inode, iface, idir, knode
    integer :: i
    integer :: n_LGL_2d, nodesperface_max

    real(wp) :: dx
    real(wp), dimension(3) :: wrk
    !real(wp), dimension(3) :: xg_target=(/1.5_wp,1.0_wp,0.0_wp/)
    logical                :: testing = .false.

    ! number of nodes in each element

    nodesperface_max = (npoly_max+1)**(ndim-1)
    allocate(facenodenormal(3,nfacesperelem*nodesperface_max,ihelems(1):ihelems(2)))
    facenodenormal = 0.0_wp

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)

       n_LGL_2d = elem_props(2,ielem)**2

      ! reset facial node index counter
      knode = 0
      ! compute outward facing normals
      !
      ! loop over faces
      do iface = 1,nfacesperelem
        ! loop over nodes on face
        do inode = 1,n_LGL_2d
          ! update facial node index counter
          knode = knode + 1
          ! volumetric node index of facial node
          i = kfacenodes(inode,iface)
          ! unsigned direction of face in computational space
          idir = abs(elfacedirections(iface))
          ! sign so normal is facing outward
          dx = sign(1.0_wp,real(elfacedirections(iface),wp))
          ! outward facing normal using metrics
          facenodenormal(:,knode,ielem) = dx*r_x(idir,:,i,ielem)
        end do
      end do
    end do

    ! testing facenodenormal calculations

    if(testing) then
      do ielem = ihelems(1), ihelems(2)
        knode = 0
        ! loop over faces
        do iface = 1,nfacesperelem
          ! loop over nodes on face
          kelem = ef2e(2,iface,ielem)
          do inode = 1,n_LGL_2d
            knode = knode + 1
            if(ef2e(1,iface,ielem) > 0)then
              i = (ef2e(1,iface,ielem)-1)*n_LGL_2d+efn2efn(4,knode,ielem)
              wrk = facenodenormal(1:3,knode,ielem)*Jx_r(kfacenodes(inode,iface),ielem) &
                  + facenodenormal(1:3, i ,kelem)*Jx_r(efn2efn(1,knode,ielem),kelem)
              if(magnitude(wrk) >= 1.0e-10_wp) then
                write(*,*)'facenodenormals are incorrect'
                write(*,*)facenodenormal(1:2,knode,ielem)*Jx_r(kfacenodes(inode,iface),ielem) &
                        , facenodenormal(1:2, i ,kelem)*Jx_r(efn2efn(1,knode,ielem),kelem)
              endif
            endif
          end do
        end do
      end do
    endif

  end subroutine calcfacenormals_LGL

  !============================================================================
  
  subroutine calcfacenormals_Gau()
    ! this subroutine calculates the outward facing normals
    ! of each facial node
    use referencevariables
    use initcollocation,      only: JacobiP11, ExtrpXa2XB_2D_neq, Gauss_Legendre_points
    use collocationvariables, only: n_Gau_1d_p1, x_Gau_1d_p1, &
                                    n_LGL_1D_p0, n_LGL_1D_p1, n_Gau_2D_p1, &
                                    Rot_Gau_p0_2_LGL_p0_1d,   &
                                    Rot_Gau_p1_2_LGL_p1_1d,   &
                                    Int_Gau_p1_2_LGL_p0_1d,   &
                                    elem_props

    use variables, only: Jx_facenodenormal_Gau, xg_Gau_Shell
    use variables, only: Jx_facenodenormal_LGL, ef2e

    implicit none

    ! indices
    integer :: ielem, inode, iface, idir, knode, node_id
    integer :: n_LGL_1d, n_Gau_1d, n_pts_1d_max, n_pts_2d_max
    integer :: n_S_1d_On, n_S_1d_Off, n_S_1d_Mort

    real(wp) :: dx
    logical                               :: testing = .false.
    real(wp), dimension(:),   allocatable :: x_LGL_1d
    real(wp), dimension(:),   allocatable :: x_Gau_1d, w_Gau_1d
    real(wp), dimension(:),   allocatable :: x_S_1d_Mort, x_S_1d_On, w_S_1d_Mort
    real(wp), dimension(:,:), allocatable :: Intrp

!   real(wp), dimension(3) :: wrk

    n_pts_1d_max = (npoly_max+1)**1
    n_pts_2d_max = (npoly_max+1)**2

    allocate(Jx_facenodenormal_Gau(3,nfacesperelem*n_pts_2d_max,ihelems(1):ihelems(2)))
    Jx_facenodenormal_Gau = 0.0_wp

    allocate(Jx_facenodenormal_LGL(3,nfacesperelem*n_pts_2d_max,ihelems(1):ihelems(2)))
    Jx_facenodenormal_Gau = 0.0_wp

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)

      n_S_1d_On   = elem_props(2,ielem)
      if(allocated(x_S_1d_On )) deallocate(x_S_1d_On ) ; allocate(x_S_1d_On (n_S_1d_On )) ; x_S_1d_On (:) = 0.0_wp
      call JacobiP11(n_S_1d_On -1,x_S_1d_On)

      knode = 0                                  !  reset facial node index counter
      node_id = 0
      do iface = 1,nfacesperelem                 ! loop over faces

        n_S_1d_Off  = ef2e(4,iface,ielem)
        n_S_1d_Mort = max(n_S_1d_On,n_S_1d_Off)

        if(allocated(x_S_1d_Mort)) deallocate(x_S_1d_Mort) ; allocate(x_S_1d_Mort(n_S_1d_Mort)) ;
        if(allocated(w_S_1d_Mort)) deallocate(w_S_1d_Mort) ; allocate(w_S_1d_Mort(n_S_1d_Mort)) ;
        call Gauss_Legendre_points(n_S_1d_Mort,x_S_1d_Mort,w_S_1d_Mort)
      
       ! compute outward facing normals on Mortar
        call Shell_Metrics_Analytic(iface,n_pts_1d_max,n_S_1d_Mort,x_S_1d_Mort,     &
                              xg_Gau_shell(:,:,ielem),Jx_facenodenormal_Gau(:,:,ielem))

        ! loop over nodes on face
        do inode = 1,n_pts_2d_max
          ! update facial node index counter
          knode = knode + 1
          ! unsigned direction of face in computational space
          idir = abs(elfacedirections(iface))
          ! sign so normal is facing outward
          dx = sign(1.0_wp,real(elfacedirections(iface),wp))
          ! outward facing normal using metrics
          Jx_facenodenormal_Gau(:,knode,ielem) = dx*Jx_facenodenormal_Gau(:,knode,ielem)
        end do

        if(allocated(Intrp)) deallocate(Intrp) ;
!       if(n_S_1d_Mort == n_S_1d_On) then
!         allocate(Intrp(n_S_1d_On,n_S_1d_Mort)) ;  Intrp(:,:) = Rot_Gau_p1_2_LGL_p1_1d(:,:) ;
!       else
!         allocate(Intrp(n_S_1d_On,n_S_1d_Mort)) ;  Intrp(:,:) = Int_Gau_p1_2_LGL_p0_1d(:,:) ;
!       endif

!       call ExtrpXA2XB_2D_neq(3, n_S_1d_Mort, n_S_1d_On,x_S_1d_Mort,x_S_1d_On, &
!          Jx_facenodenormal_Gau(:,:,ielem),Jx_facenodenormal_LGL(:,:,ielem),Intrp)

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

    return
  end function cross_product


  subroutine calcmetrics_LGL()
    ! This subroutine calculates the metric transformations
    ! between computational and physical space.
    use referencevariables
    use variables, only: xg, x_r, r_x, Jx_r, dx_min_elem
    use collocationvariables, only: iagrad, jagrad, dagrad, pvol
    use mpimod

    implicit none
    ! indices
    integer :: ielem, inode, jnode, idir, jdir
    integer :: i,j,k,l,ii

    real(wp) :: test(3)
    real(wp) :: err,err_L2,err_Linf
    real(wp), dimension(:), allocatable :: err_max_proc

    integer :: s_tag, r_tag, m_size, m, &
               s_request_err_Linf, r_request_err_Linf, i_err

    integer :: s_status(mpi_status_size)
    integer :: r_status(mpi_status_size)

    continue 

    ! dx/dr
    allocate(x_r(3,3,1:nodesperelem,ihelems(1):ihelems(2)))
    x_r = 0.0_wp
    ! dr/dx
    allocate(r_x(3,3,1:nodesperelem,ihelems(1):ihelems(2)))
    r_x = 0.0_wp
    ! J = |dx/dr|
    allocate(Jx_r(1:nodesperelem,ihelems(1):ihelems(2)))
    Jx_r = 0.0_wp

    allocate(dx_min_elem(ihelems(1):ihelems(2)))
    dx_min_elem(:) = 0.0_wp

    ! Initialize metrics error norms
    err_L2 = 0.0_wp
    err_Linf = 0.0_wp

    ! loop over volumetric elements
    elloop:do ielem = ihelems(1), ihelems(2)
      ! initialize dx/dr to identity and dr/dx to identity
      do idir = 1,3
        x_r(idir,idir,:,ielem) = 1.0_wp
        r_x(idir,idir,:,ielem) = 1.0_wp
      end do
      ! calculate dx/dr
      !
      ! initialize to zero
      x_r(:,1:ndim,:,ielem) = 0.0_wp
      ! loop over every node in element
      do inode = 1, nodesperelem
        ! loop over dimension of gradient
        do jdir = 1,ndim
          ! loop over number of dependent nodes in gradient
          do i = iagrad(inode), iagrad(inode+1)-1
            ! column/node from gradient operator in CSR format in
            ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
            jnode = jagrad(jdir,i)
            ! update gradient. MP: Well, actually this is the Jacobian of the transformation
            x_r(:,jdir,inode,ielem) = x_r(:,jdir,inode,ielem) + dagrad(jdir,i)*xg(:,jnode,ielem)
          end do
        end do
      end do
      ! calculate determinant
      do inode = 1,nodesperelem
        Jx_r(inode,ielem) = determinant3(x_r(:,:,inode,ielem))
      end do
      if (ndim < 3) then
        ! inverse metrics (note that in 3D this is not sufficient to satisfy the GCL.
        do inode = 1,nodesperelem
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
        do inode = 1, nodesperelem
          ! dxi_1/dx_1
          j = 1 ; k = 1 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 3
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(2,2,jnode,ielem)*xg(3,jnode,ielem)
          end do
          jdir = 2
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(2,3,jnode,ielem)*xg(3,jnode,ielem)
          end do

          ! dxi_1/dx_2
          j = 1 ; k = 2 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 3
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(3,2,jnode,ielem)*xg(1,jnode,ielem)
          end do
          jdir = 2
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(3,3,jnode,ielem)*xg(1,jnode,ielem)
          end do

          ! dxi_1/dx_3
          j = 1 ; k = 3 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 3
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(1,2,jnode,ielem)*xg(2,jnode,ielem)
          end do
          jdir = 2
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(1,3,jnode,ielem)*xg(2,jnode,ielem)
          end do

          ! dxi_2/dx_1
          j = 2 ; k = 1 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 1
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(2,3,jnode,ielem)*xg(3,jnode,ielem)
          end do
          jdir = 3
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(2,1,jnode,ielem)*xg(3,jnode,ielem)
          end do

          ! dxi_2/dx_2
          j = 2 ; k = 2 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 1
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(3,3,jnode,ielem)*xg(1,jnode,ielem)
          end do
          jdir = 3
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(3,1,jnode,ielem)*xg(1,jnode,ielem)
          end do

          ! dxi_2/dx_3
          j = 2 ; k = 3 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 1
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(1,3,jnode,ielem)*xg(2,jnode,ielem)
          end do
          jdir = 3
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(1,1,jnode,ielem)*xg(2,jnode,ielem)
          end do

          ! dxi_3/dx_1
          j = 3 ; k = 1 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 2
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(2,1,jnode,ielem)*xg(3,jnode,ielem)
          end do
          jdir = 1
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(2,2,jnode,ielem)*xg(3,jnode,ielem)
          end do

          ! dxi_3/dx_2
          j = 3 ; k = 2 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 2
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(3,1,jnode,ielem)*xg(1,jnode,ielem)
          end do
          jdir = 1
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(3,2,jnode,ielem)*xg(1,jnode,ielem)
          end do

          ! dxi_3/dx_3
          j = 3 ; k = 3 ;
          r_x(j,k,inode,ielem) = 0.0_wp
          jdir = 2
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) + dagrad(jdir,i)*x_r(1,1,jnode,ielem)*xg(2,jnode,ielem)
          end do
          jdir = 1
          do i = iagrad(inode), iagrad(inode+1)-1
            jnode = jagrad(jdir,i)
            r_x(j,k,inode,ielem) = r_x(j,k,inode,ielem) - dagrad(jdir,i)*x_r(1,2,jnode,ielem)*xg(2,jnode,ielem)
          end do

          r_x(:,:,inode,ielem) = r_x(:,:,inode,ielem)/Jx_r(inode,ielem)

        end do
      end if

!     Metric Test
! 
!     \frac{\partial}{\partial \xi^j} \left( J \frac{\partial \xi^j}{\partial x^i} \right) == 0 ; i,j = 1,3
!     ( Assumes Einstein Convention on ``j'' )
!
      do inode = 1,nodesperelem
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

      !  Calculate a conservative estimate of the characteristic length of an
      !  element.  Needs further refinement

      dx_min_elem(ielem) = 0.0_wp
      do inode = 1,nodesperelem
        dx_min_elem(ielem) = dx_min_elem(ielem) + pvol(inode)*Jx_r(inode,ielem)
      end do

      dx_min_elem(ielem) = dx_min_elem(ielem)**(0.333333333333333333333333_wp)

    end do elloop

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

    return
  end subroutine calcmetrics_LGL

  !============================================================================

  subroutine modify_metrics_nonconforming()

    use referencevariables
!   use collocationvariables
    use initcollocation, only: GCL_Triple_Qmat_Transpose, element_properties
    use eispack_module,  only: svd
    use unary_mod,       only: qsortd

    implicit none

    logical, parameter :: matu = .True.
    logical, parameter :: matv = .True.

    integer,  dimension(:),    allocatable :: ifacenodes
    real(wp), dimension(:,:),  allocatable :: qmat
    real(wp), dimension(:),    allocatable :: p_surf

    real(wp), dimension(:,:),  allocatable :: Amat
    real(wp), dimension(:,:),  allocatable :: bvec
    real(wp), dimension(:,:),  allocatable :: a_t
    real(wp), dimension(:,:),  allocatable :: u, v
    real(wp), dimension(:),    allocatable :: w, work, pmat
    real(wp), dimension(:,:),  allocatable :: work3
    real(wp), dimension(:,:),  allocatable :: eye, wrk, diag, wI
    real(wp), dimension(:,:),  allocatable :: delta_a

    integer,  dimension(:),    allocatable :: perm
    integer :: ielem, ierr
    integer :: i
    integer :: n_pts_1d, n_pts_2d, n_pts_3d
    integer :: nm, m, n

    logical                                :: testing = .false.
    real(wp)                               :: t1, t2
    real(wp), parameter                    :: tol = 1.0e-12_wp

    ! loop over volumetric elements
    elloop:do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,&
                              n_pts_1d=n_pts_1d,&
                              n_pts_2d=n_pts_2d,&
                              n_pts_3d=n_pts_3d,&
                                  qmat=qmat,&
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
          if(w(i) >= tol) wI(i,i) = 1.0_wp / w(i)
        enddo
        call qsortd(w(1:m),perm,m)
        wrk(:,:) = matmul(transpose(v(1:m,1:m)),v(1:m,1:m)) &
                 + matmul(v(1:m,1:m),transpose(v(1:m,1:m))) &
                 + matmul(transpose(u(1:n,1:m)),u(1:n,1:m)) 
        t1 = maxval(abs(3*eye(:,:) - wrk(:,:))) 
        t2 = maxval(abs(transpose(amat)-matmul(u(1:n,1:m),matmul(diag(1:m,1:m),transpose(v(1:m,1:m))) )))
        if(w(perm(2)) <= tol) write(*,*)'second singular mode',w(perm(1:2))
        if(t1+t2 > tol) then
          write(*,*)'error in u^T u', t1
          write(*,*)'error in A - u S v^T', t2
          write(*,*)'stopping'
          stop
         endif
      endif

      call Load_Mortar_Metric_Data(ielem,n_pts_1d,n_pts_2d,n_pts_3d, ifacenodes, p_surf, a_t, bvec)

!     minimizing (a - a_t)(a - a_t) / 2 ; subject to M a = b
!     "t" subscript denotes "target metric values' :  i.e., a_t
!     a = a_t - M^* (M a_t - b)
      
      if(allocated(work3)) deallocate(work3) ; allocate(work3(1:m,1:3))   ; work3(:,:)   = 0.0_wp
      work3(1:m,1:3) = matmul(Amat(1:m,1:n),a_t(1:n,1:3)) - bvec(1:m,1:3)
      t1 = maxval(abs(work3))
      if(t1 >= tol) write(*,*)'A a_t - bvec', maxval(abs(work3))

      delta_a = matmul( u(1:n,1:m),           &
                  matmul(wI(1:m,1:m),           &
                    matmul(transpose(v(1:m,1:m)), &
                      matmul(Amat(1:m,1:n),a_t(1:n,1:3)) - bvec(1:m,1:3))))

!     write(*,*)'perturbation magnitude',maxval(abs(delta_a))
!     aa = a_t - delta_a
     
    end do elloop
    deallocate(u,v,w,work,Amat,a_t,bvec,wI)

  end subroutine modify_metrics_nonConforming

  !============================================================================

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

    integer                                  :: i, iface, inode, jnode

    real(wp), dimension(3)                   :: nx

    bvec(:,:) = 0.0_wp
    do iface = 1,nfacesperelem

       do i = 1,n_LGL_2d

            ! Index in facial ordering
            jnode =  n_LGL_2d*(iface-1)+i
              
            ! Volumetric node index corresponding to facial node index
            inode = ifacenodes(jnode)
              
            ! Outward facing normal of facial node
            if (ef2e(4,iface,ielem) == elem_props(2,ielem)) then !    Conforming interface
              nx(:) = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)
!             write(*,*)'metric diffeences',maxval(abs(nx(:) - Jx_facenodenormal_LGL(:,jnode,ielem)))
!             write(*,*)'metrics',nx(1),Jx_facenodenormal_LGL(1,jnode,ielem)
            else                                                 ! NonConforming interface
!             write(*,*)'ef2e(4)', iface, ef2e(1,iface,ielem),ef2e(4,iface,ielem), elem_props(2,ielem)
!             nx(:) = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)
              nx(:) = Jx_facenodenormal_LGL(:,jnode,ielem)
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
                       & if2nq, ifacetag, ic2nh, nqface, vx_master, &
                       & periodic_face_data_x1, periodic_face_data_x2, &
                       & periodic_face_data_x3, wall_face_data

    ! Nothing is implicitly defined
    implicit none

    integer :: i,j,k, j1, j2, k1
    integer :: ielem,cnt,cntBC 
    integer :: iface, jface, bigN, ave

    !integer, dimension(:),   allocatable :: iav2e_tmp
    integer, dimension(:,:), allocatable :: jav2e_tmp
    integer, dimension(:),   allocatable :: ivtmp1, ivtmp2, ivtmp3, ivtmp4

    integer :: nnzv2e
    !integer, allocatable :: v2e(:,:), iav2e(:), jav2e(:)
    integer, allocatable :: iav2e(:), jav2e(:)

    integer, dimension(:,:), allocatable :: test_conL, test_conR
    integer, dimension(:),   allocatable :: test_cBCL, test_cBCR

    logical                              :: testing = .true.

    ! Variables for periodic faces
    integer,  dimension(:),     allocatable :: ivtmp_iface, ivtmp_jface
    real(wp), dimension(:,:),   allocatable :: vx_iface,    vx_jface
    integer,  dimension(:,:),     allocatable :: same_coord 
    integer :: i_comp_face, i_vertex_iface, i_vertex_jface, i_coord, i_p_elem, &
      & j_p_elem, match, jelem, i_vertex_hexa, cnt_periodic_elem_x1, &
      & cnt_periodic_elem_x2, cnt_periodic_elem_x3

    integer :: cnt_wall_elem

    integer,  dimension(:),   allocatable :: iv_hexa_elem 

    integer :: jelem_iface, jface_iface, ielem_jface, iface_jface, p_face, &
      & i_check 

    !integer, allocatable, dimension(:) :: list_partner_faces

    real(wp), parameter  :: diff_toll = 1e-8
    integer,  parameter  :: qdim = 6             !  dimension of ef2e array

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
    !  ef2e       :    ( 6 ,nfaceperelem, nelements) 
    !             :  Two situation occur.  The face is either an 
    !                  (Interior face 
    !                      :  (1,j,k) = Adjoining element face ID
    !                      :  (2,j,k) = Adjoining element ID
    !                      :  (3,j,k) = Adjoining element process ID
    !                      :  (4,j,k) = Adjoining element polynomial order
    !                      :  (5,j,k) = Number of Adjoining elements
    !                      :  (6,j,k) = HACK self polynomial order assigned to each face
    !                  (Boundary face 
    !                      :  (1,j,k) = Set to -11 
    !                      :  (2,j,k) = -100000000
    !                      :  (3,j,k) = -100000000
    !                      :  (4,j,k) = -100000000
    !                      :  (5,j,k) = -100000000
    !                      :  (6,j,k) = -100000000
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
      allocate(ef2e(qdim,2*ndim,1:nelems))  ;   ef2e(:,:,:) = -1000000000

      allocate(ivtmp1(nverticesperface*bigN),ivtmp2(nverticesperface*bigN))
      allocate(ivtmp3(nverticesperface),     ivtmp4(nverticesperface))

      allocate(ivtmp_iface(nverticesperface))
      allocate(ivtmp_jface(nverticesperface))

      allocate(vx_iface(3,nverticesperface))
      allocate(vx_jface(3,nverticesperface))
      allocate(same_coord(nverticesperface,3))
      allocate(iv_hexa_elem(nverticesperelem))

      ! Counter for the number of elements that have a periodic boundary face
      cnt_periodic_elem_x1 = 0
      cnt_periodic_elem_x2 = 0
      cnt_periodic_elem_x3 = 0

      ! Counter for the number of elements that have a wall boundary face
      cnt_wall_elem = 0

      ! Counter for the number of boundary conditions
      cntBC = 0

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
          !         write(*,*)'elem',ielem,'face',iface,'touches element',ave,'on face',jface
                    test_conL(ielem,iface) = ielem*iface ; test_conR(ave ,jface) = ave*jface
                    cycle faceloop
                  endif

                enddo

              endif
            endif

          end do 

          cnt = 0 ; ivtmp3(:) = 0 
          do j = 1, nverticesperface
            cnt = cnt + 1           
            j1 = eltypfaces(j,iface)
            ivtmp3(cnt) = ic2nh(j1,ielem)    
          end do
          ivtmp3 = isort(ivtmp3,nverticesperface)
  !       write(*,*)'ielem',ielem,'face',iface,'ivtmp3',ivtmp3(:)

          do jface = 1,nqface
            ivtmp4(:) = 0 ; ivtmp4(:) = if2nq(:,jface)    
            ivtmp4 = isort(ivtmp4,nverticesperface)
            if(maxval(abs(ivtmp3(:)-ivtmp4(:))) == 0)  then
              ef2e(1,iface,ielem) = -ifacetag(jface)

              ! Build periodic faces data in the x1 direction
              if(ifacetag(jface) == 8 .or. ifacetag(jface) == 9) then

                cnt_periodic_elem_x1 = cnt_periodic_elem_x1 + 1

                ! Global ID of the element which owns a "periodic" boundary face
                periodic_face_data_x1(1,cnt_periodic_elem_x1) = ielem 

                ! Get the local ID of the "periodic" boundary face of the ielem
                periodic_face_data_x1(2,cnt_periodic_elem_x1) = iface
                
                ! Position in the if2nq stack of the iface
                periodic_face_data_x1(3,cnt_periodic_elem_x1) = jface

                do j = 1, nverticesperface
                  periodic_face_data_x1(4+j,cnt_periodic_elem_x1) = ivtmp4(j)
                end do
              end if


              ! Build periodic faces data in the x2 direction
              if(ifacetag(jface) == 10 .or. ifacetag(jface) == 11) then

                cnt_periodic_elem_x2 = cnt_periodic_elem_x2 + 1

                ! Global ID of the element which owns a "periodic" boundary face
                periodic_face_data_x2(1,cnt_periodic_elem_x2) = ielem 

                ! Get the local ID of the "periodic" boundary face of the ielem
                periodic_face_data_x2(2,cnt_periodic_elem_x2) = iface
                
                ! Position in the if2nq stack of the iface
                periodic_face_data_x2(3,cnt_periodic_elem_x2) = jface

                do j = 1, nverticesperface
                  periodic_face_data_x2(4+j,cnt_periodic_elem_x2) = ivtmp4(j)
                end do
              end if


              ! Build periodic faces data in the x3 direction
              if(ifacetag(jface) == 12 .or. ifacetag(jface) == 13) then

                cnt_periodic_elem_x3 = cnt_periodic_elem_x3 + 1

                ! Global ID of the element which owns a "periodic" boundary face
                periodic_face_data_x3(1,cnt_periodic_elem_x3) = ielem 

                ! Get the local ID of the "periodic" boundary face of the ielem
                periodic_face_data_x3(2,cnt_periodic_elem_x3) = iface
                
                ! Position in the if2nq stack of the iface
                periodic_face_data_x3(3,cnt_periodic_elem_x3) = jface

                do j = 1, nverticesperface
                  periodic_face_data_x3(4+j,cnt_periodic_elem_x3) = ivtmp4(j)
                end do
              end if


              ! Store data for computing the aerodynamic coefficients
              if(ifacetag(jface) == 5 .or. ifacetag(jface) == 6) then
                cnt_wall_elem = cnt_wall_elem + 1
                ! Global ID of the element which owns a "wall" boundary face
                wall_face_data(1,cnt_wall_elem) = ielem 

                ! Get the local ID of the "wall" boundary face of the ielem
                wall_face_data(2,cnt_wall_elem) = iface
                
              end if

              cntBC = cntBC + 1 ; test_cBCL(cntBC) = jface
              cycle faceloop
            endif
          enddo

          ! Search should always find a connection but didn't. Somethings wrong
          write(*,*) 'Face search failed. element and face' ; write(*,*) ielem,iface
          stop
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
      ! Build ef2e for "periodic" faces in the x1 direction
      ! ===================================================

      ! Loop over the element that owns a "periodic" face
      periodic_elem_x1_1 : do i_p_elem = 1, size(periodic_face_data_x1(1,:))

        ! Get the global ID of the element that owns a periodic boundary face
        ielem = periodic_face_data_x1(1,i_p_elem)

        ! Get the local (local for the element) ID of the periodic boundary face
        iface = periodic_face_data_x1(2,i_p_elem)

        ! Check if the face is a "periodic" face with tag equals to 8
        if (abs(ef2e(1,iface,ielem)) == 8) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
          end do

          ! Search the partner face, i.e. the corresponding "periodic" face with
          ! face tag equals to 9
          partner_face_x1_1: do i_comp_face = 1, size(periodic_face_data_x1(3,:))
            
            ! Get the face ID. This ID is given by the AFLR3 format
            jface = periodic_face_data_x1(3,i_comp_face)
            
            ! Partner face can only be a face with ifacetag = 9
            if (ifacetag(jface) == 9) then
                  
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
              search_x1_1: do i_vertex_iface = 1, nverticesperface
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
                        check_common_coord_x1_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            cycle check_common_coord_x1_1
                          else
                            !write(*,*) same_coord(1,:)
                            !write(*,*) same_coord(i_check,:)
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
                    cycle partner_face_x1_1
                  end if
                endif

              end do search_x1_1
             
              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x1(4,i_p_elem) = i_coord
                end if
              end do     

              ! Search the global ID of the element which owns the jface
              search_x1_2 : do j_p_elem = 1, size(periodic_face_data_x1(1,:))
                
                ! Get the global ID of the element
                jelem = periodic_face_data_x1(1,j_p_elem)

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
                !if (match == nverticesperface .and. jelem /= ielem) then
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

!                  write(*,*) 'ielem, jelem', ielem, jelem

!                  do i_vertex_jface = 1, nverticesperface
!                    do i_vertex_hexa = 1, nverticesperelem
!                      if (iv_hexa_elem(i_vertex_hexa) == ivtmp_jface(i_vertex_jface)) then
!                        ic2nh_mod(i_vertex_hexa,jelem) = periodic_face_data(3+i_vertex_jface,i_p_elem)
!                      endif
!                    end do
!                  end do
               
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
                  cycle search_x1_2

                  endif

                  exit
                
                enddo search_x1_2
              
              exit
            
            endif

          end do partner_face_x1_1

        endif

      end do periodic_elem_x1_1


!      write(*,*) 'NOW THE SYMMETRIC PART'
      
      periodic_elem_x1_2 : do i_p_elem = 1, size(periodic_face_data_x1(1,:))

!        write(*,*) 'IN bc_elem loop'

        ! Get the global ID of the element that owns a boundary face
        ielem = periodic_face_data_x1(1,i_p_elem)

        ! Get the local (local for the element) ID of the boundary face
        iface = periodic_face_data_x1(2,i_p_elem)

        ! Check if the face is a "periodic" face
        if (abs(ef2e(1,iface,ielem)) == 9) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
          end do

          ! Search the partner face, i.e. the corresponding "periodic" face
          partner_face_x1_2: do i_comp_face = 1, size(periodic_face_data_x1(3,:))
            
            ! Get the face ID
            ! This ID is the one given by the AFLR3 format
            jface = periodic_face_data_x1(3,i_comp_face)
            
            ! Partner face can only be a face with ifacetag = 8
            if (ifacetag(jface) == 8) then
                  
              ! Get the ID of the nodes which form the boundary face
              ivtmp_jface(:) = if2nq(:,jface)

              ! Get the coordinate of the nodes of the jface
              do i_vertex_jface = 1, nverticesperface
                vx_jface(1,i_vertex_jface) = vx_master(1,ivtmp_jface(i_vertex_jface))
                vx_jface(2,i_vertex_jface) = vx_master(2,ivtmp_jface(i_vertex_jface))
                vx_jface(3,i_vertex_jface) = vx_master(3,ivtmp_jface(i_vertex_jface))
              end do

!              write(*,*)
!              write(*,*) 'node coordinate jface'
!              do i_vertex_jface = 1, nverticesperface 
!                write(*,*) 'node', i_vertex_jface
!                write(*,*) 'x', vx_jface(1,i_vertex_jface) 
!                write(*,*) 'y', vx_jface(2,i_vertex_jface) 
!                write(*,*) 'z', vx_jface(3,i_vertex_jface) 
!              end do
!              write(*,*)
               

              ! Set to zero the array that keeps track of the coordinates
              ! matches
              same_coord = 0

              ! Check if the jface is a partner face.
              ! If the jface is the partner face of the iface then, two of its 
              ! node coordinate must be equal to the coordinates of the nodes 
              ! which form the iface. The remaining coordinate just differs for 
              ! a shift.
              search_x1_3: do i_vertex_iface = 1, nverticesperface
                do i_vertex_jface = 1, nverticesperface
                  same_coord(i_vertex_iface,:) = 0
                  do i_coord = 1, 3
!                    write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                    write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1 
                    endif
                  end do

                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
!                    write(*,*) 'Found node'
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common
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
                  !write(*,*) 'No the right face'
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

!                write(*,*) 'same_coord', same_coord(i_vertex_iface)

              end do search_x1_3
             
!              write(*,*) 'Match found' !, exiting from the companion_face loop'
!              write(*,*) 'ID vertices jface', ivtmp_jface

              ! Store the periodic direction (x or y or z) of each face
              do i_coord = 1, 3
                ! Just pick up the first node for this check.
                if (same_coord(1,i_coord) == 0) then
                  periodic_face_data_x1(4,i_p_elem) = i_coord
                end if
              end do

              ! Search the global ID of the element which owns the jface
              search_x1_4 : do j_p_elem = 1, size(periodic_face_data_x1(1,:))
                
                ! Get the global ID of the element
                jelem = periodic_face_data_x1(1,j_p_elem)

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


!      write(*,*) 'NOW THE SYMMETRIC PART'
      
      periodic_elem_x3_2 : do i_p_elem = 1, size(periodic_face_data_x3(1,:))

!        write(*,*) 'IN bc_elem loop'

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


      ! Check if the faces are connected in a symmetric fashion
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

      ! Cosntruct iae2v and jae2v
      ! =========================
      nnze2v = nverticesperelem*nelems
      allocate(iae2v(nelems+1))  ; iae2v = 0
      allocate(iae2v_tmp(nelems+1))  ; iae2v_tmp = 0
      allocate(jae2v(nnze2v))    ; jae2v = 0
      allocate(jae2v_tmp(nnze2v))    ; jae2v_tmp = 0

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

    return
  end subroutine e2e_connectivity_aflr3     !  SERIAL Routine

  !============================================================================
  
    subroutine set_element_orders_Serial()     !  Serial Routine

    ! Load modules
    use variables
    use collocationvariables
    use referencevariables, only : npoly, npoly_max, nfacesperelem
    use mpimod

    ! Nothing is implicitly defined
    implicit none
   
    integer :: ielem, j, nhex, iface

    npoly_max = -1000

    nhex = size(ic2nh,2)

    if(allocated(elem_props)) deallocate(elem_props) ; allocate(elem_props(2,1:nhex))
    elem_props(:,:) = -1000

    do ielem = 1,nhex

      elem_props(1,ielem) = 1
      elem_props(2,ielem) = npoly+1
      do j=1,8
        if(vx_master(1,ic2nh(j,ielem)) >= 100000.5_wp) elem_props(2,ielem) = npoly+2
      enddo
 
    enddo

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
    use referencevariables, only : ihelems

    ! Nothing is implicitly defined
    implicit none
   
    integer :: ielem, qdim

    if(allocated(elem_props)) deallocate(elem_props) ; allocate(elem_props(2,ihelems(1):ihelems(2)))

    do ielem = ihelems(1),ihelems(2)

      qdim = size(ef2e(:,1,1))
      elem_props(1,ielem) = 1
      if(sum(ef2e(6,:,ielem))/qdim /= ef2e(6,1,ielem)) then
         write(*,*)'mpi bug in transfering ef2e'
      endif
      elem_props(2,ielem) = ef2e(6,1,ielem)
 
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

    return
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

    return
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

    return
  end subroutine data_partner_element_serial     !   SERIAL Routine

  pure function WENO_Adjoining_Data(k_node,k_face)     !   PARALLEL Routine
! function WENO_Adjoining_Data(k_node,k_face)
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
     endif

  end function WENO_Adjoining_Data     !   PARALLEL Routine

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
    allocate(xgWENO_partner(3,nodespershell,ihelems(1):ihelems(2))) ; xgWENO_partner = -10000.0_wp ;
    allocate(XIWENO_partner(3,nodespershell,ihelems(1):ihelems(2))) ; XIWENO_partner = -10000.0_wp ;
    allocate(xghstWENO_partner(3,nghost)) ; xghstWENO_partner = 0.0_wp


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

    call UpdateComm1DGhostDataWENOGeom(xgWENO_self, xghstWENO_partner, &
                 xpetscWENO_partner, xlocpetscWENO_partner, 3, nodespershell, ihelems, nghost)

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
                   - xghst_LGL(:,iloc) + xg(:,inode,ielem)                   ! account for possibility of non-periodic domain.
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
!            write(*,*)'xgwenopartner',xgWENO_partner(:,jnode,ielem)
!            write(*,*)'XIwenopartner',XIWENO_partner(:,jnode,ielem)
             call WENO_xi_val(comp2phys_coeffs,XIWENO_partner(:,jnode,ielem),xgWENO_partner(:,jnode,ielem),rnorm)
             if(rnorm >= rnorm_max) rnorm_max = rnorm
             tmp = maxval(abs(XIWENO_partner(:,jnode,ielem)))
             if( tmp > 1.0_wp + tol) then
                if( tmp >= XI_max) XI_max = tmp
                extrnal_xi_cnt = extrnal_xi_cnt + 1
   !            write(*,*)ielem,iface,ipen,XIWENO_partner(:,jnode,ielem)
             endif
!            write(*,*)jnode,XIWENO_partner(1,jnode,ielem),XIWENO_partner(2,jnode,ielem),XIWENO_partner(3,jnode,ielem)
!            Test the (x0,y0,z0) -> (xi0,eta0,zeta0) decode step for a uniform grid.  Comment if non-uniform
 !           tol = 1.0e-10_wp
 !           do i = 1,3
 !             tmp = abs(XIWENO_partner(i,jnode,ielem))
 !             t1  = tmp - 1.0_wp      
 !             t2  = tmp - 1.0_wp/sqrt5
 !             if(min(t1,t2) >= tol) write(*,*)'error in XIWENO_partner',XIWENO_partner(i,jnode,ielem)
 !           enddo
!
          enddo
        end if

      enddo

    enddo
!   write(*,*)'extrnal_xi_cnt',myprocid,extrnal_xi_cnt
    extrnal_xi_sum = 0 ;
    call mpi_allreduce(extrnal_xi_cnt,extrnal_xi_sum,1, &
      & MPI_INT,MPI_SUM,petsc_comm_world,i_err)
!   if(myprocid == 0) write(*,*)'i_err',i_err

    call mpi_allreduce(XI_max,XI_max_Glob,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,i_err)
!   if(myprocid == 0) write(*,*)'i_err',i_err



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

    return
  end function LGL_pts_lexo_comp_hexa

  !=================================================================================================

  subroutine calc_Gau_shell_pts_all_hexas()

    ! Load module
    use variables
    use referencevariables
    use collocationvariables
    use initcollocation, only : lagrange_basis_function_1d, Gauss_Legendre_points, JacobiP11

    ! Nothing is implicitly defined
    implicit none

    real(wp), allocatable, dimension(:,:) :: Gau_pts_comp_shell_one_face

    integer  :: i_elem, i_Gau, i_LGL, j_LGL, k_LGL
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
    do i_elem = ihelems(1), ihelems(2)

      n_LGL_1d_On   = elem_props(2,i_elem)
      if(allocated(x_LGL_1d_On)) deallocate(x_LGL_1d_On) ; allocate(x_LGL_1d_On(n_LGL_1d_On))
      call JacobiP11(n_LGL_1d_On-1,x_LGL_1d_On)

      do iface = 1,nfacesperelem

        n_LGL_1d_Off  = ef2e(4,iface,i_elem)
        n_Gau_1d_Mort = max(n_LGL_1d_On, n_LGL_1d_Off)
        n_Gau_2d_Mort = n_Gau_1d_Mort**2

        if(allocated(x_Gau_1d_Mort)) deallocate(x_Gau_1d_Mort) ; allocate(x_Gau_1d_Mort(n_Gau_1d_Mort)) ;
        if(allocated(w_Gau_1d_Mort)) deallocate(w_Gau_1d_Mort) ; allocate(w_Gau_1d_Mort(n_Gau_1d_Mort)) ;

        call Gauss_Legendre_points(n_Gau_1d_Mort,x_Gau_1d_Mort,w_Gau_1d_Mort)

        Gau_pts_comp_shell_one_face = Shell_pts_lexo_comp_one_face(iface,n_Gau_2d_max, &
                                                                   n_Gau_1d_mort,x_Gau_1d_Mort)

        do i_Gau = 1, n_Gau_2d_Mort

          ishift = (iface-1) * n_Gau_2d_Mort + i_Gau
          l = 0
          xg_Gau_shell(:,ishift,i_elem) = 0.0_wp
  
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
  
                xg_Gau_shell(:,ishift,i_elem) = xg_Gau_shell(:,ishift,i_elem) + xg(:,l,i_elem)*l_xi*l_eta*l_zeta
  
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
    use variables, only: xg, ef2e,  &
      & xg_Gau_shell, xgghst_Gau_Shell, efn2efn_Gau, &
      & jelems, periodic_elem_face_ids_x1, &
      & periodic_elem_face_ids_x2, periodic_elem_face_ids_x3
    use initcollocation, only : lagrange_basis_function_1d, Gauss_Legendre_points, JacobiP11

    ! Nothing is implicitly defined
    implicit none

    integer ::  ielem, inode, jnode, iface, knode, kshell
    integer ::  i_low

    real(wp) :: x1(3), x2(3)
    real(wp), parameter :: nodetol = 1.0e-8_wp

    integer :: i_p_face, p_dir, cnt_coord, i_coord
    logical :: match_found
    real(wp), dimension(2) :: x1_p, x2_p

    integer :: n_Gau_1d_max , n_Gau_2d_max , n_Gau_shell_max
    integer :: n_Gau_1d_Mort, n_Gau_2d_Mort
    integer :: n_LGL_1d_On  , n_LGL_1d_Off
    integer :: cnt_debug

    continue

    cnt_debug = 0

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

    ! Initialize position of the ghost point in the stack
    i_low = 0

    ! Loop over elements
    do ielem = ihelems(1), ihelems(2)
      
      ! Reset facial node index counter

      n_LGL_1d_On   = elem_props(2,ielem)
      
      ! Loop over faces
      do iface = 1, nfacesperelem

        n_LGL_1d_Off  = ef2e(4,iface,ielem)
        n_Gau_1d_Mort = max(n_LGL_1d_On, n_LGL_1d_Off)
        n_Gau_2d_Mort = n_Gau_1d_Mort**2
        knode = (iface-1) * n_Gau_2d_max

        ! If on boundary, connect to self
        if (ef2e(1,iface,ielem) < 0) then
          
          ! Loop over nodes on the boundary face
          do inode = 1, n_Gau_2d_Mort
            
            ! Update facial node index counter
            knode = knode + 1
            
            ! The first index is the volumetric node index and the second
            ! index is the element index
            efn2efn_Gau(:,knode,ielem) = (/ -1000, ielem, 0 /)
          
          end do

        else if (ef2e(3,iface,ielem) /= myprocid) then ! A parallel interface

          ! Initialize match_found
          match_found = .false.
          
          ! Loop through the elements that owns a periodic face in the x1
          ! direction
          if (size(periodic_elem_face_ids_x1(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic
            ! face
            do i_p_face = 1, size(periodic_elem_face_ids_x1(1,:))

              if (periodic_elem_face_ids_x1(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x1(2,i_p_face) == iface) then

                ! There is a match: change logical value of match_found
                match_found = .true.

                ! Get the direction of "periodicity"
                p_dir = periodic_elem_face_ids_x1(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_Gau_2d_Mort
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node

                  x1 = xg_Gau_shell(:,knode,ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0
                  
                  do i_coord = 1,3
                    
                    if (i_coord /= p_dir) then

                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  do jnode = 1, n_Gau_2d_Mort
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xgghst_Gau_Shell(:,i_low + jnode)
                
                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0

                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn_Gau(1,knode,ielem) = -1000 

                      ! Set the element of the connected node
                      efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)

                      ! Set the node index in the ghost array
                      efn2efn_Gau(3,knode,ielem) = i_low + jnode

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

                ! Update the position in the ghost stack
                i_low = i_low + n_Gau_2d_Mort

              end if ! End if match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit 
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if check periodic face in x1 direction


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

                ! Loop over the nodes on the face
                do inode = 1, n_Gau_2d_Mort
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg_Gau_shell(:,knode,ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0
                  
                  do i_coord = 1,3
                    
                    if (i_coord /= p_dir) then

                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  do jnode = 1, n_Gau_2d_Mort
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xgghst_Gau_Shell(:,i_low + jnode)
                
                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0

                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn_Gau(1,knode,ielem) = -1000

                      ! Set the element of the connected node
                      efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)

                      ! Set the node index in the ghost array
                      efn2efn_Gau(3,knode,ielem) = i_low + jnode

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

                ! Update the position in the ghost stack
                i_low = i_low + n_Gau_2d_Mort

              end if ! End if match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit 
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if check periodic face in x2 direction


          ! Loop through the elements that owns a periodic face in the x3
          ! direction
          if (match_found .eqv. .false. .and. size(periodic_elem_face_ids_x3(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if iface is a periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x3(1,:))

              if (periodic_elem_face_ids_x3(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x3(2,i_p_face) == iface) then

                ! There is a match: change logical value of match_found
                match_found = .true.

                ! Get the direction of "periodicity"
                p_dir = periodic_elem_face_ids_x3(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_Gau_2d_Mort
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg_Gau_shell(:,knode,ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0
                  
                  do i_coord = 1,3
                    
                    if (i_coord /= p_dir) then

                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  do jnode = 1, n_Gau_2d_Mort
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xgghst_Gau_Shell(:,i_low + jnode)
                
                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0

                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn_Gau(1,knode,ielem) = -1000

                      ! Set the element of the connected node
                      efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)

                      ! Set the node index in the ghost array
                      efn2efn_Gau(3,knode,ielem) = i_low + jnode

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

                ! Update the position in the ghost stack
                i_low = i_low + n_Gau_2d_Mort

              end if ! End if match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit 
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if check periodic face in x3 direction

          if (match_found .eqv. .false.) then

            ! Loop over the nodes on the face
            do inode = 1, n_Gau_2d_Mort

              ! Update the facial node index counter
              knode = knode + 1

              ! Save the coordinates of the facial node
              x1 = xg_Gau_shell(:,knode,ielem)
              
              ! Search for the connected node on face of the connected element
              do jnode = 1, n_Gau_2d_Mort

                ! Coordinates of the jnode
                ! ef2e(2) gives the element of the neighbor
                x2 = xgghst_Gau_Shell(:,i_low + jnode)
                
                ! Check the distance between the two nodes
                if (magnitude(x1-x2) <= nodetol) then
                  
                  ! Set the volumetric node index of the connected node
                  efn2efn_Gau(1,knode,ielem) = -1000
                  
                  ! Set the element of the connected node
                  efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)
                  
                  ! Set the node index in the ghost array
                  efn2efn_Gau(3,knode,ielem) = i_low + jnode

                  exit
                
                end if
              
              end do
              
              ! Print information at screen if there is a problem and stop
              ! computation
              if (jnode > n_Gau_2d_Mort .and. myprocid==1) then
                write(*,*) 'Connectivity error in face-node connectivity.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
                write(*,*) 'Node coordinates and ghost node coordinates'
                write(*,*) x1, xgghst_Gau_Shell(:,i_low + 1:i_low + n_Gau_2d_Mort)
                write(*,*) 'Exiting...'
                stop
              end if

            end do

            ! Update the position in the ghost stack
            i_low = i_low + n_Gau_2d_Mort
          
          end if

        else ! Not a parallel interface

          ! Initialize match_found
          match_found = .false.

          if (size(periodic_elem_face_ids_x1(1,:)) /= 0) then

            ! Check if the ielem owns a periodic face and if the iface is a periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x1(1,:))

              if (periodic_elem_face_ids_x1(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x1(2,i_p_face) == iface) then

                ! There is a match
                match_found = .true.

                ! Get the direction of periodicity
                p_dir = periodic_elem_face_ids_x1(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_Gau_2d_Mort
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg_Gau_shell(:,knode,ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0 

                  do i_coord = 1, 3
                    
                    if (i_coord /= p_dir) then
                      
                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  ! Search for the connected node on the face of the connected 
                  ! element
                  do jnode = 1,n_Gau_2d_Mort
                    ! Coordinates of the jnode
                    ! ef2e(1) gives the face on the neighboring element and
                    ! ef2e(2) gives the element

                    kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_Mort + jnode
                    x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))

                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0
                    
                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check the distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn_Gau(1,knode,ielem) = -1000
                      
                      ! Set the element of the connected node
                      efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)
                      
                      ! Set the index of the connected node
                      efn2efn_Gau(4,knode,ielem) = kshell

                      cnt_debug = cnt_debug + 1


                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

              end if ! Match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit
              end if

            end do ! End do loop over the elements that own a periodic face

          end if ! End if periodic x1 direction

          ! If the iface is not a periodic face  in the x1 direction, check
          ! if it is a periodic face in the x2 direction
          if ((match_found .eqv. .false.) .and.              &
             (size(periodic_elem_face_ids_x2(1,:)) /= 0)) then

            ! Check if the ielem owns a periodic face and if the iface is a periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x2(1,:))

              if (periodic_elem_face_ids_x2(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x2(2,i_p_face) == iface) then

                ! There is a match
                match_found = .true.

                ! Get the direction of periodicity
                p_dir = periodic_elem_face_ids_x2(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_Gau_2d_Mort
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg_Gau_shell(:,knode,ielem)
                
                  ! Extract from x1 the two invaraint coordinates
                  cnt_coord = 0 

                  do i_coord = 1, 3
                    
                    if (i_coord /= p_dir) then
                      
                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  ! Search for the connected node on the face of the connected 
                  ! element
                  do jnode = 1,n_Gau_2d_Mort
                    ! Coordinates of the jnode
                    ! ef2e(1) gives the face on the neighboring element and
                    ! ef2e(2) gives the element
                    kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_Mort + jnode
                    x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))

                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0
                    
                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check the distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn_Gau(1,knode,ielem) = -1000
                      
                      ! Set the element of the connected node
                      efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)
                      
                      ! Set the index of the connected node
                      efn2efn_Gau(4,knode,ielem) = kshell

                      cnt_debug = cnt_debug + 1

                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

              end if ! Match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit
              end if

            end do ! End do loop over the elements that own a periodic face
          
          end if ! End if periodic x2 direction


          ! If the iface is not a periodic face in the x2 direction, check
          ! if it is a periodic face in the x3 direction
          if ((match_found .eqv. .false.) .and.              &
             (size(periodic_elem_face_ids_x3(1,:)) /= 0)) then
           
!           write(*,*)'periodic 3'
            ! Check if the ielem owns a periodic face and if the iface is a periodic face
            do i_p_face = 1, size(periodic_elem_face_ids_x3(1,:))

              if (periodic_elem_face_ids_x3(1,i_p_face) == jelems(ielem) .and. &
                & periodic_elem_face_ids_x3(2,i_p_face) == iface) then

                ! There is a match
                match_found = .true.

                ! Get the direction of periodicity
                p_dir = periodic_elem_face_ids_x3(3,i_p_face)

                ! Loop over the nodes on the face
                do inode = 1, n_Gau_2d_Mort
                  
                  ! Update the facial node index counter
                  knode = knode + 1
                  
                  ! Save the coordinates of the facial node
                  x1 = xg_Gau_shell(:,knode,ielem)
                
                  ! Extract from x1 the two invariant coordinates
                  cnt_coord = 0 

                  do i_coord = 1, 3
                    
                    if (i_coord /= p_dir) then
                      
                      cnt_coord = cnt_coord + 1
                      x1_p(cnt_coord) = x1(i_coord)
                    
                    end if
                  
                  end do

                  ! Search for the connected node on the face of the connected element
                  do jnode = 1,n_Gau_2d_Mort
                    ! Coordinates of the jnode
                    ! ef2e(1) gives the face on the neighboring element and
                    ! ef2e(2) gives the element
                    kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_Mort + jnode
                    x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))

                    ! Extract from x2 the two invaraint coordinates
                    cnt_coord = 0
                    
                    do i_coord = 1, 3
                      
                      if (i_coord /= p_dir) then
                        
                        cnt_coord = cnt_coord + 1
                        x2_p(cnt_coord) = x2(i_coord)
                      
                      end if
                    
                    end do

                    ! Check the distance between the two nodes
                    if (magnitude(x1_p-x2_p) <= nodetol) then
                      
                      ! Set the volumetric node index of the connected node
                      efn2efn_Gau(1,knode,ielem) = -1000
                      
                      ! Set the element of the connected node
                      efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)
                      
                      ! Set the index of the connected node
                      efn2efn_Gau(4,knode,ielem) = kshell

                      cnt_debug = cnt_debug + 1


                      exit ! partner jnode found; exit the jnode do loop
                    
                    end if
                  
                  end do ! End do jnode
                
                end do ! End do inode 

              end if ! Match found

              ! If a partner face has been found exit from the loop over the 
              ! elements that own a periodic face
              if (match_found .eqv. .true.) then
                exit
              end if

            end do ! End do loop over the elements that own a periodic face
          
          end if ! End if periodic x3 direction


          if (match_found .eqv. .false.) then

!           write(*,*)'non-periodic path'
            
            ! Loop over the nodes on the face
            do inode = 1, n_Gau_2d_Mort

              ! Update the facial node index counter
              knode = knode + 1

              ! Save coordinates of the facial ndoes
              x1 = xg_Gau_shell(:,knode,ielem)
              ! Search the for connected node on the face of the connected element
              
              do jnode = 1, n_Gau_2d_Mort
                
                ! Coordinates of the jnode
                ! ef2e(1) gives the face on the neighboring element and ef2e(2) gives the element
                kshell = (ef2e(1,iface,ielem)-1)*n_Gau_2d_Mort + jnode
                x2 = xg_Gau_shell(:,kshell,ef2e(2,iface,ielem))
                
                ! Check the distance between the two nodes
                if (magnitude(x1-x2) <= nodetol) then

                  ! Set the volumetric node index of the connected node
                  efn2efn_Gau(1,knode,ielem) = -1000
                  
                  ! Set the element of the connected node
                  efn2efn_Gau(2,knode,ielem) = ef2e(2,iface,ielem)
                  
                  ! Set the index of the connected node
                  efn2efn_Gau(4,knode,ielem) = kshell
                  
                  exit
                
                end if
              
              end do ! End do jnode

              ! Print information at screen if there is a problem and stop computation
!             if (efn2efn_Gau(2,knode,ielem) < 0) then
!               write(*,*) 'Connectivity error in face-node connectivity of Gauss path.'
!               write(*,*) 'Process ID, element ID, face ID, ef2e'
!               write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
!               write(*,*) 'Node coordinates'
!               write(*,*) x1
!               write(*,*) 'Possible partner node coordinates'
!               
!               do jnode = 1, n_Gau_2d_Mort
!                 x2 = xg(:,kfacenodes_Gau(jnode,ef2e(1,iface,ielem)),ef2e(2,iface,ielem))
!                 write(*,*) x2
!               end do 

!               write(*,*) 'Exiting...'
!               stop
!             end if

            end do ! End do inode
          
          end if ! End if not a periodic face (match_found = .false.)
              
        end if ! End if type of face (boundary, off processor or on processor)
      
      end do ! End do loop over faces of the element
    
    end do ! End do loop elements owned by the processor

    return
  end subroutine calculate_face_node_connectivity_Gau

  !============================================================================

  subroutine Shell_Metrics_Analytic(iface, n_pts_1d_max, n_pts_1d, x_pts_1d,  &
                                    xg_Gau,Jx_facenodenormal_Gau)

    ! Nothing is implicitly defined
    use referencevariables
    use initcollocation,      only: D_lagrange_basis_function_1d, &
                                      lagrange_basis_function_1d, &
                                      JacobiP11

!   use variables, only:  Jx_r, facenodenormal

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
    integer  :: ival, jval, kval
    integer  :: ishift, jshift, kshift, face_shift

    real(wp), dimension(:), allocatable  :: x_Gau
    
    continue

    node_id = 0
    n_Gau = n_pts_1d
    face_shift = (iface-1)*n_pts_1d_max
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
              kval = mod(node_id-1,n_Gau)+1 ; jval = 1 + (node_id-kval) / n_Gau ;

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
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  t1*t2 - t3*t4

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
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  t1*t2 - t3*t4

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
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  t1*t2 - t3*t4
            end do
          end do
          Jx_facenodenormal_Gau(:,:) =  Jx_facenodenormal_Gau(:,:) * zeta
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
              kval = mod(node_id-1,n_Gau)+1 ; ival = 1 + (node_id-kval) / n_Gau ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  t1*t2 - t3*t4
            end do
        end do
        Jx_facenodenormal_Gau(:,:) =  Jx_facenodenormal_Gau(:,:) * eta
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
              ival = mod(node_id-1,n_Gau)+1 ; jval = 1 + (node_id-ival) / n_Gau ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + ival + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = node_id - ival + i1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + ival + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = node_id - ival + i1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + ival + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = node_id - ival + i1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  t1*t2 - t3*t4
          end do
        end do
        Jx_facenodenormal_Gau(:,:) =  Jx_facenodenormal_Gau(:,:) * xi
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
              kval = mod(node_id-1,n_Gau)+1 ; ival = 1 + (node_id-kval) / n_Gau ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do i1 = 1, n_pts_1d
                 ishift = (i1-1)*n_Gau + kval + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t4 = t4 + D_lagrange_basis_function_1d(  xi,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  t1*t2 - t3*t4
            end do
        end do
        Jx_facenodenormal_Gau(:,:) =  Jx_facenodenormal_Gau(:,:) * eta
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
              ival = mod(node_id-1,n_Gau)+1 ; jval = 1 + (node_id-ival) / n_Gau ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + ival + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = node_id - ival + i1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
              enddo
              Jx_facenodenormal_Gau(1,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + ival + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = node_id - ival + i1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(3,ishift)
              enddo
              Jx_facenodenormal_Gau(2,node_id+face_shift) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + ival + face_shift
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
              enddo
              do i1 = 1, n_pts_1d
                 ishift = node_id - ival + i1 + face_shift
                 t2 = t2 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(2,ishift)
                 t3 = t3 + D_lagrange_basis_function_1d(zeta,i1,x_Gau,n_Gau) * xg_Gau(1,ishift)
              enddo
              Jx_facenodenormal_Gau(3,node_id+face_shift) =  t1*t2 - t3*t4
          end do
        end do
        Jx_facenodenormal_Gau(:,:) =  Jx_facenodenormal_Gau(:,:) * xi
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
              kval = mod(node_id-1,n_Gau)+1 ; jval = 1 + (node_id-kval) / n_Gau ;

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval
                 t1 =  t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
                 t4 =  t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1
                 t2 =  t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
                 t3 =  t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
              enddo   
              Jx_facenodenormal_Gau(1,node_id) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(3,jshift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1
                 t2 = t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(3,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
              enddo
              Jx_facenodenormal_Gau(2,node_id) =  t1*t2 - t3*t4

              t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
              do j1 = 1, n_pts_1d
                 jshift = (j1-1)*n_Gau + kval
                 t1 = t1 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(2,jshift)
                 t4 = t4 + D_lagrange_basis_function_1d( eta,j1,x_Gau,n_Gau) * xg_Gau(1,jshift)
              enddo
              do k1 = 1, n_pts_1d
                 kshift = node_id - kval + k1
                 t2 = t2 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(1,kshift)
                 t3 = t3 + D_lagrange_basis_function_1d(  xi,k1,x_Gau,n_Gau) * xg_Gau(2,kshift)
              enddo
              Jx_facenodenormal_Gau(3,node_id) =  t1*t2 - t3*t4
            end do
          end do
        Jx_facenodenormal_Gau(:,:) =  Jx_facenodenormal_Gau(:,:) * zeta
      endif

    return
  end subroutine Shell_Metrics_Analytic

  !============================================================================

end module initgrid

