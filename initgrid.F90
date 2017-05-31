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
  public E2EConnectivity
  public e2e_connectivity_aflr3
  public calculatepartitions
  public calcnodes
  public calcmetrics
  public facenodesetup
  public facenodesetupWENO
  public calculate_face_node_connectivity
  public calcfacenormals
  public init_elem_type
  public create_ldg_flip_flop_sign
  public pert_int_vert
  public face_pairs
  public data_partner_element_serial
  public WENO_Adjoining_Data
  public Pencil_Coord
  public WENO_Intrp_Face_Nodes
  public Boundary_Vertex_2_Vertex_Connectivity

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

  subroutine init_edge_2()
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

  end subroutine init_edge_2

  subroutine init_quad_4()
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

  end subroutine init_quad_4

  subroutine init_hex_8()
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

  end subroutine init_hex_8

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

  subroutine init_elem_type()
    use referencevariables
    use variables, only: boundaryelems, ef2e,  &
                         iae2v, jae2v, nnze2v, &
                         facenormalcoordinate
!                        iae2e, jae2e, nnze2e, 
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

  end subroutine init_elem_type

  subroutine E2EConnectivity()
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

!   nnze2e = connectioncount
!   !   local work (i.e. in initgrid)
!   allocate(iae2e(nelems+1))
!   iae2e = 0
!   allocate(jae2e(nnze2e))
!   jae2e = 0

!   ii = 1
!   do ielem = 1, nelems
!     iae2e(ielem) = ii
!     do j = 1, nfacesperelem
!       ! only add element connections to other elements
!       if ( ef2e(1,j,ielem) > 0 ) then
!         jae2e(ii) = ef2e(2,j,ielem)
!         ii = ii+1
!       end if
!     end do
!   end do
!   iae2e(nelems+1) = ii

    deallocate(ivtmp1,ivtmp2)

    deallocate(iav2e,v2e)


  end subroutine E2EConnectivity

  subroutine calculatepartitions()
    ! this subroutine contains the c-bindings for
    ! calling the metis library on the "grid", defined by
    ! the element-to-node connectivities read from the datafile.
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
    use collocationvariables, only: rcollocation
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
    use collocationvariables, only: rcollocation
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
  
  !============================================================================

  subroutine perturb_vertices_tg_vortex_1(p_scale)

    ! Load modules
    use referencevariables
    use variables, only: xg, vx, e2v
    use collocationvariables, only: rcollocation
    
    ! Nothing is implicitly defined
    implicit none
    
    real(wp), intent(in) :: p_scale
    integer :: ielem
    integer :: i, j, izone
    real(wp) :: dx(3), dr
    integer,parameter :: seed = 86456
    real(wp) :: rand
    
    integer :: low_elem, high_elem

    real(wp) :: diff_x, diff_y, diff_z
    real(wp), parameter :: toll = 1e-6
    real :: tmp
    
    continue
  
    low_elem = ihelems(1)
    high_elem = ihelems(2)

    call srand(seed)

    do ielem = low_elem, high_elem
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
  
  !============================================================================

  subroutine pert_int_vert(p_scale)

    ! Load modules
    use referencevariables
    use variables, only: vx_master
    use collocationvariables, only: rcollocation
    
    ! Nothing is implicitly defined
    implicit none
    
    real(wp), intent(in) :: p_scale
    integer :: ielem
    integer :: i, j, izone
    real(wp) :: dx(3), dr
    integer,parameter :: seed = 86456
    real(wp) :: rand
    
    integer :: n_tot_vertices, i_vertex

    real(wp) :: diff_x, diff_y, diff_z
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
  
  !============================================================================

  subroutine calcnodes()
    ! this subroutine populates the nodal grid locations for
    ! each element. Currently we assume that all elements are
    ! straight sided, but this can be remedied by incorporating
    ! CAD or analytical surface data. 
    use controlvariables, only: Grid_Topology, cylinder_x0, cylinder_x1
    use referencevariables
    use variables, only: xg, vx, e2v, ef2e
    use collocationvariables, only: rcollocation
    implicit none
    ! indices
    integer :: ielem, inode, idir, iface
    integer :: i,j,k
    integer :: nE
    ! cartesian based grid coordinates
    real(wp), allocatable :: xl(:,:,:,:)
    ! high and low indices for each direction
    integer :: il(2,3)
    ! local grid distance
    real(wp)                :: dr
    real(wp), dimension(nodesperedge)  :: xi
    real(wp), dimension(3)  :: dx, x0
    real(wp), dimension(3)  :: x00,x01

    integer,  dimension(6)  :: curved_faces
    
   
    ! nE is simply for convenience of presentation in the coding
    nE = nodesperedge
    ! number of nodes in each element
    nodesperelem = nE**ndim
    ! total number of nodes
    nnodes = nodesperelem*nelems

    ! allocate global node matrix
    allocate(xg(3,1:nodesperelem,ihelems(1):ihelems(2)))
    xg = 0.0_wp
    ! allocate local nodes
    allocate(xl(3,1:nE,1:nE,1:nE))
    xl = 0.0_wp

    ! low index is always 1
    il = 1
    ! set high index for each grid direction to nodesperedge
    do idir = 1,ndim
      il(2,idir) = nE
    end do

    ! loop over volumetric elements
    do ielem = ihelems(1), ihelems(2)
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
            dr = 0.5_wp*(rcollocation(i)+1.0_wp)    ! distance in computational space
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

        xi(:) = 0.5_wp*(rcollocation(:)+1.0_wp)    ! distance in computational space

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

      end select

      ! build faces
      if (ndim > 1) then
        ! xi_3 = 0
        call TFI2D(xl(:, :, :, 1),nE,rcollocation)
      end if
      if (ndim > 2) then
        ! xi_3 = 1
        call TFI2D(xl(:, :, :,nE),nE,rcollocation)
        ! xi_2 = 0
        call TFI2D(xl(:, :, 1, :),nE,rcollocation)
        ! xi_2 = 1
        call TFI2D(xl(:, :,nE, :),nE,rcollocation)
        ! xi_1 = 0
        call TFI2D(xl(:, 1, :, :),nE,rcollocation)
        ! xi_1 = 1
        call TFI2D(xl(:,nE, :, :),nE,rcollocation)
      end if
      ! build volumes
      if (ndim > 2) then
        call TFI3D(xl(:,:,:,:),nE,rcollocation)
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

  end subroutine calcnodes

! =============================================================================

! =============================================================================

  subroutine facenodesetup()
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
    use referencevariables
    use mpimod
    use variables, only: xg, kfacenodes, ifacenodes, ef2e, efn2efn, &
      & boundaryelems
    implicit none

    ! indices
    integer :: i,j,k
    integer :: stride, stride1, stride2, ioffset

    real(wp), parameter :: nodetol = 1.0e-8_wp

    ! local facial masks
    !
    ! kfacenodes separates each face
    allocate(kfacenodes(nodesperface,nfacesperelem))
    ! ifacenodes includes all faces
    allocate(ifacenodes(nodesperface*nfacesperelem))

    if (ndim == 2) then
      ! loop over every node on each face
      do i = 1, nodesperedge
        ! on face 1, the first nodesperface nodes are just
        ! the first nodesperface
        kfacenodes(i,1) = i
        ! on face 3 there is just an offset to where the
        ! counting starts
        j = (nodesperedge-1)*nodesperedge+i
        kfacenodes(i,3) = j
        ! onface 2, a stride and offset are required
        stride = nodesperedge
        j = nodesperedge + (i-1)*stride
        kfacenodes(i,2) = j
        ! on face 4, a stride is needed
        j = 1 + (i-1)*stride
        kfacenodes(i,4) = j 
      end do
    else if (ndim == 3) then
      k = 0
      do j = 1, nodesperedge
        do i = 1, nodesperedge
          k = k+1
          ! face 1 does not require an offset or a stride
          kfacenodes(k,1) = k
          ! on face 2, a stride is required
          ioffset = 1
          stride1 = 1
          stride2 = nodesperedge**2
          kfacenodes(k,2) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 3, offset and stride are needed
          ioffset = nodesperedge
          stride1 = nodesperedge
          stride2 = nodesperedge**2
          kfacenodes(k,3) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! face 4 requires an offset and a stride
          ioffset = 1+(nodesperedge-1)*nodesperedge
          stride1 = 1
          stride2 = nodesperedge**2
          kfacenodes(k,4) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 5 only a stride is required
          ioffset = 1
          stride1 = nodesperedge
          stride2 = nodesperedge**2
          kfacenodes(k,5) = ioffset+stride1*(i-1)+stride2*(j-1)
          ! on face 6 only an offset is required
          ioffset = (nodesperedge-1)*nodesperedge*nodesperedge
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
      do i = 1, nodesperface
        ! advance facial node index
        k = k+1
        ! map facial node index to volumetric node
        ifacenodes(k) = kfacenodes(i,j)
      end do
    end do

  end subroutine facenodesetup

! =============================================================================

  subroutine facenodesetupWENO()
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
    use variables, only: xg, kfacenodesWENO, ifacenodesWENO, ef2e, efn2efn, &
      & boundaryelems
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

  end subroutine facenodesetupWENO

  !============================================================================
  
  !============================================================================
  ! calculate_face_node_connectivity - Sets the face-node connectivity for the
  ! collocation points.
  subroutine calculate_face_node_connectivity()
    
    ! Load modules
    use referencevariables
    use mpimod
    use variables, only: xg, xghst, kfacenodes, ifacenodes, ef2e, efn2efn, &
      & boundaryelems, jelems, periodic_elem_face_ids_x1, &
      & periodic_elem_face_ids_x2, periodic_elem_face_ids_x3

    ! Nothing is implicitly defined
    implicit none

    integer ::  ielem, inode, jnode, iface, knode
    integer ::  i_low

    integer :: low_elem, high_elem
    real(wp) :: x1(3), x2(3)
    real(wp), parameter :: nodetol = 1.0e-8_wp

    integer :: i_p_face, p_dir, cnt_coord, i_coord
    logical :: match_found
    real(wp), dimension(2) :: x1_p, x2_p

    integer :: cnt_debug

    continue

    cnt_debug = 0

    ! Low volumetric element index
    low_elem = ihelems(1)

    ! High volumetric element index
    high_elem = ihelems(2)


    ! efn2efn contains the partner node information of every facenode in the 
    ! domain
    allocate(efn2efn(4,nfacesperelem*nodesperface,low_elem:high_elem))
    efn2efn = -1000

    ! Initialize position of the ghost point in the stack
    i_low = 0

    ! Loop over elements
    do ielem = low_elem, high_elem
      
      ! Reset facial node index counter
      knode = 0
      
      ! Loop over faces
      do iface = 1, nfacesperelem
        
        ! If on boundary, connect to self
        if (ef2e(1,iface,ielem) < 0) then
          
          ! Loop over nodes on the boundary face
          do inode = 1, nodesperface
            
            ! Update facial node index counter
            knode = knode + 1
            
            ! The first index is the volumetric node index and the second
            ! index is the element index
            efn2efn(:,knode,ielem) = (/ ifacenodes(knode), ielem, 0 /)
          
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
                do inode = 1, nodesperface
                  
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

                  do jnode = 1, nodesperface
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xghst(:,i_low + jnode)
                
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
                i_low = i_low + nodesperface

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
                do inode = 1, nodesperface
                  
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

                  do jnode = 1, nodesperface
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xghst(:,i_low + jnode)
                
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
                i_low = i_low + nodesperface

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
                do inode = 1, nodesperface
                  
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

                  do jnode = 1, nodesperface
                    
                    ! Coordinates of the jnode
                    ! ef2e(2) gives the element of the neighbor
                    x2 = xghst(:,i_low + jnode)
                
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
                i_low = i_low + nodesperface

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
            do inode = 1, nodesperface

              ! Update the facial node index counter
              knode = knode + 1

              ! Save the coordinates of the facial node
              x1 = xg(:,ifacenodes(knode),ielem)
              
              ! Search for the connected node on face of the connected element
              do jnode = 1, nodesperface

                ! Coordinates of the jnode
                ! ef2e(2) gives the element of the neighbor
                x2 = xghst(:,i_low + jnode)
                
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
              if (jnode > nodesperface .and. myprocid==1) then
                write(*,*) 'Connectivity error in face-node connectivity.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
                write(*,*) 'Node coordinates and ghost node coordinates'
                write(*,*) x1, xghst(:,i_low + 1:i_low + nodesperface)
                write(*,*) 'Exiting...'
                stop
              end if

            end do

            ! Update the position in the ghost stack
            i_low = i_low + nodesperface
          
          end if

        else ! Not a parallel interface

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
                do inode = 1, nodesperface
                  
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
                  do jnode = 1,nodesperface
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

!                      if (ielem == 1) then
!                        write(*,*) 'ielem', ielem
!                        write(*,*) 'x1_p', x1_p
!                        write(*,*) 'x2_p', x2_p
!                        write(*,*) 'x1', x1
!                        write(*,*) 'x2', x2
!                      end if

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
                do inode = 1, nodesperface
                  
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
                  do jnode = 1,nodesperface
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

!                      if (ielem == 1) then
!                        write(*,*) 'ielem', ielem
!                        write(*,*) 'x1_p', x1_p
!                        write(*,*) 'x2_p', x2_p
!                        write(*,*) 'x1', x1
!                        write(*,*) 'x2', x2
!                      end if

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
                do inode = 1, nodesperface
                  
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
                  do jnode = 1,nodesperface
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

!                      if (ielem == 1) then
!                        write(*,*) 'ielem', ielem
!                        write(*,*) 'x1_p', x1_p
!                        write(*,*) 'x2_p', x2_p
!                        write(*,*) 'x1', x1
!                        write(*,*) 'x2', x2
!                      end if

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
            do inode = 1, nodesperface

              ! Update the facial node index counter
              knode = knode + 1

              ! Save coordinates of the facial ndoes
              x1 = xg(:,ifacenodes(knode),ielem)
              ! Search the for connected node on the face of the connected 
              ! element
              
              do jnode = 1, nodesperface
                
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

              ! Print information at screen if there is a problem and stop
              ! computation
              if (efn2efn(1,knode,ielem) < 0 .or. efn2efn(2,knode,ielem) < 0) then
                write(*,*) 'Connectivity error in face-node connectivity.'
                write(*,*) 'Process ID, element ID, face ID, ef2e'
                write(*,*) myprocid, ielem, iface, ef2e(:,iface,ielem)
                write(*,*) 'Node coordinates'
                write(*,*) x1
!                write(*,*) 'Type of face', ef2e(1,iface,ielem)
!                write(*,*) 'Process ID', ef2e(3,iface,ielem)
                write(*,*) 'Possible partner node coordinates'
                
                do jnode = 1, nodesperface
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
  end subroutine calculate_face_node_connectivity


  !============================================================================
  
  !============================================================================
  
  subroutine calcfacenormals()
    ! this subroutine calculates the outward facing normals
    ! of each facial node
    use referencevariables
    use variables, only: kfacenodes, ifacenodes, &
      & boundaryelems, facenodenormal, r_x, xg, ef2e, efn2efn, Jx_r
    implicit none

    ! indices
    integer :: ielem, kelem, inode, iface, idir, knode
    integer :: i

    real(wp) :: dx
    real(wp), dimension(3) :: wrk
    !real(wp), dimension(3) :: xg_target=(/1.5_wp,1.0_wp,0.0_wp/)
    logical                :: testing = .false.

    allocate(facenodenormal(3,nfacesperelem*nodesperface,ihelems(1):ihelems(2)))
    facenodenormal = 0.0_wp

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! reset facial node index counter
      knode = 0
      ! compute outward facing normals
      !
      ! loop over faces
      do iface = 1,nfacesperelem
        ! loop over nodes on face
        do inode = 1,nodesperface
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
          do inode = 1,nodesperface
            knode = knode + 1
            if(ef2e(1,iface,ielem) > 0)then
              i = (ef2e(1,iface,ielem)-1)*nodesperface+efn2efn(4,knode,ielem)
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

  end subroutine calcfacenormals


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


  subroutine calcmetrics()
    ! This subroutine calculates the metric transformations
    ! between computational and physical space.
    use referencevariables
    use variables, only: xg, x_r, r_x, Jx_r, e2v, dx_min_elem
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
            ! update gradient. MP: Well, actually this is the Jacobian of the 
            ! transformation
            x_r(:,jdir,inode,ielem) = x_r(:,jdir,inode,ielem) &
              + dagrad(jdir,i)*xg(:,jnode,ielem)
          end do
        end do
      end do
      ! calculate determinant
      do inode = 1,nodesperelem
        Jx_r(inode,ielem) = determinant3(x_r(:,:,inode,ielem))
      end do
      if (ndim < 3) then
        ! inverse metrics (note that in 3D this is not sufficient to satisfy
        ! the GCL.
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
  end subroutine calcmetrics

  !============================================================================
  
  !============================================================================
  ! e2e_connectivity_aflr3 - Constructs the element-to-element connectivity
  ! starting from the information read from the AFLR3 grid.

  subroutine e2e_connectivity_aflr3()

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

    real(wp), parameter :: diff_toll = 1e-8

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
    !  ef2e       :    ( 2 ,nfaceperelem, nelements) 
    !             :  Two situation occur.  The face is either an 
    !                  (Interior face 
    !                      :  (1,j,k) = face ID of the adjoining element
    !                      :  (2,j,k) = Connected to Element 
    !                  (Boundary face 
    !                      :  (1,j,k) = Set to -11 
    !                      :  (2,j,k) = -100000000
    !
    !         (Note:  redimensioned (3,:,:) to account for processor info)
    !
    ! iae2v,jae2v     :    Which vertices belong to each element

    if (ndim == 2) then
      write(*,*) '2D grid connectivity for the AFLR3 format is not yet', &
        & ' implemented'
      write(*,*) 'Exting...'
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
      allocate(ef2e(3,2*ndim,1:nelems))  ;   ef2e(:,:,:) = -1000000000

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

     
!      do i = 1, nelems
!        iv_hexa_elem = ic2nh(:,i)
!
!        write(*,*)
!        write(*,*) 'Element', i
!        do j = 1, 8
!          write(*,*) 'Node:', j, iv_hexa_elem(j)
!          write(*,*) 'x', vx_master(1,iv_hexa_elem(j))
!          write(*,*) 'y', vx_master(2,iv_hexa_elem(j))
!          write(*,*) 'z', vx_master(3,iv_hexa_elem(j))
!        end do
!
!        write(*,*)
!      enddo


!      do i = 1, 8
!        write(*,*) vx_master(:,ic2nh(i,231))
!        write(*,*) vx_master(:,ic2nh(i,595))
!
!      end do
!
!      write(*,*) 'Partner element'
!      do i = 1, 8
!        write(*,*) vx_master(:,ic2nh(i,595))
!      end do
!      stop


      ! ===================================================
      ! ===================================================
      ! Build ef2e for "periodic" faces in the x1 direction
      ! ===================================================
      ! ===================================================

      ! Loop over the element that owns a "periodic" face
      periodic_elem_x1_1 : do i_p_elem = 1, size(periodic_face_data_x1(1,:))

!        write(*,*) 'IN bc_elem loop'

        ! Get the global ID of the element that owns a periodic boundary face
        ielem = periodic_face_data_x1(1,i_p_elem)

        ! Get the local (local for the element) ID of the periodic boundary face
        iface = periodic_face_data_x1(2,i_p_elem)

!        write(*,*) 'iface tag', abs(ef2e(1,iface,ielem))

        ! Check if the face is a "periodic" face with tag equals to 8
        if (abs(ef2e(1,iface,ielem)) == 8) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

!          write(*,*)
!          write(*,*) 'node coordinate iface'
!          do i_vertex_iface = 1, nverticesperface 
!            write(*,*) 'node', i_vertex_iface
!            write(*,*) 'x', vx_iface(1,i_vertex_iface) 
!            write(*,*) 'y', vx_iface(2,i_vertex_iface) 
!            write(*,*) 'z', vx_iface(3,i_vertex_iface) 
!          end do

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

!              write(*,*)
!              write(*,*) 'node coordinate jface'
!              do i_vertex_jface = 1, nverticesperface 
!                write(*,*) 'node', i_vertex_jface
!                write(*,*) 'x', vx_jface(1,i_vertex_jface) 
!                write(*,*) 'y', vx_jface(2,i_vertex_jface) 
!                write(*,*) 'z', vx_jface(3,i_vertex_jface) 
!              end do
!              write(*,*)

              !stop

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
!                    write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                    write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1
!                      write(*,*) 'i_vertex_iface', i_vertex_iface
!                      write(*,*) 'i_vertex_jface', i_vertex_jface
!                      write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                      write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
!                      write(*,*) 'HERE'
!                      write(*,*) 'i_coord', i_coord
                    endif
                  end do
                  
                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
!                    write(*,*) 'Found node'
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common'
                      !write(*,*) i_vertex_iface
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x1_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            !write(*,*) i_p_elem
                            !write(*,*) 'node 1',i_coord, same_coord(1,i_coord)
                            !write(*,*) 'other nodes', i_coord, same_coord(i_check,i_coord)
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
!                  write(*,*) 'Not the right face'
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

!                write(*,*) 'same_coord', same_coord(i_vertex_iface)

              end do search_x1_1
             
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
!                      write(*,*) 'ef2e 1', ef2e(1,iface,ielem)
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

!        write(*,*) 'iface tag', abs(ef2e(1,iface,ielem))

        ! Check if the face is a "periodic" face
        if (abs(ef2e(1,iface,ielem)) == 9) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x1(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

!          write(*,*)
!          write(*,*) 'node coordinate iface'
!          do i_vertex_iface = 1, nverticesperface 
!            write(*,*) 'node', i_vertex_iface
!            write(*,*) 'x', vx_iface(1,i_vertex_iface) 
!            write(*,*) 'y', vx_iface(2,i_vertex_iface) 
!            write(*,*) 'z', vx_iface(3,i_vertex_iface) 
!          end do

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
                !if (match == nverticesperface .and. jelem /= ielem) then
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

!                  write(*,*) 'ielem, jelem', ielem, jelem
               
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
!                      write(*,*) 'ef2e 1', ef2e(1,iface,ielem)
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

!        write(*,*) 'IN bc_elem loop'

        ! Get the global ID of the element that owns a periodic boundary face
        ielem = periodic_face_data_x2(1,i_p_elem)

        ! Get the local (local for the element) ID of the periodic boundary face
        iface = periodic_face_data_x2(2,i_p_elem)

!        write(*,*) 'iface tag', abs(ef2e(1,iface,ielem))

        ! Check if the face is a "periodic" face with tag equals to 8
        if (abs(ef2e(1,iface,ielem)) == 10) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

!          write(*,*)
!          write(*,*) 'node coordinate iface'
!          do i_vertex_iface = 1, nverticesperface 
!            write(*,*) 'node', i_vertex_iface
!            write(*,*) 'x', vx_iface(1,i_vertex_iface) 
!            write(*,*) 'y', vx_iface(2,i_vertex_iface) 
!            write(*,*) 'z', vx_iface(3,i_vertex_iface) 
!          end do

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

!              write(*,*)
!              write(*,*) 'node coordinate jface'
!              do i_vertex_jface = 1, nverticesperface 
!                write(*,*) 'node', i_vertex_jface
!                write(*,*) 'x', vx_jface(1,i_vertex_jface) 
!                write(*,*) 'y', vx_jface(2,i_vertex_jface) 
!                write(*,*) 'z', vx_jface(3,i_vertex_jface) 
!              end do
!              write(*,*)

              !stop
              ! Set to zero the array that keeps track of the coordinates
              ! matches
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
!                    write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                    write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1
!                      write(*,*) 'i_vertex_iface', i_vertex_iface
!                      write(*,*) 'i_vertex_jface', i_vertex_jface
!                      write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                      write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
!                      write(*,*) 'HERE'
!                      write(*,*) 'i_coord', i_coord
                    endif
                  end do
                  
                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
!                    write(*,*) 'Found node'
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common'
                      !write(*,*) i_vertex_iface
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x2_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            !write(*,*) i_p_elem
                            !write(*,*) 'node 1',i_coord, same_coord(1,i_coord)
                            !write(*,*) 'other nodes', i_coord, same_coord(i_check,i_coord)
                            cycle check_common_coord_x2_1
                          else
                            !write(*,*) same_coord(1,:)
                            !write(*,*) same_coord(i_check,:)
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
!                  write(*,*) 'Not the right face'
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

!                write(*,*) 'same_coord', same_coord(i_vertex_iface)

              end do search_x2_1
             
!              write(*,*) 'Match found' !, exiting from the companion_face loop'
!              write(*,*) 'ID vertices jface', ivtmp_jface

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
!                      write(*,*) 'ef2e 1', ef2e(1,iface,ielem)
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

!        write(*,*) 'iface tag', abs(ef2e(1,iface,ielem))

        ! Check if the face is a "periodic" face
        if (abs(ef2e(1,iface,ielem)) == 11) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x2(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

!          write(*,*)
!          write(*,*) 'node coordinate iface'
!          do i_vertex_iface = 1, nverticesperface 
!            write(*,*) 'node', i_vertex_iface
!            write(*,*) 'x', vx_iface(1,i_vertex_iface) 
!            write(*,*) 'y', vx_iface(2,i_vertex_iface) 
!            write(*,*) 'z', vx_iface(3,i_vertex_iface) 
!          end do

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
              search_x2_3: do i_vertex_iface = 1, nverticesperface
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
                  !write(*,*) 'No the right face'
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

!                write(*,*) 'same_coord', same_coord(i_vertex_iface)

              end do search_x2_3
             
!              write(*,*) 'Match found' !, exiting from the companion_face loop'
!              write(*,*) 'ID vertices jface', ivtmp_jface

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
                !if (match == nverticesperface .and. jelem /= ielem) then
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

!                  write(*,*) 'ielem, jelem', ielem, jelem
               
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
!                      write(*,*) 'ef2e 1', ef2e(1,iface,ielem)
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

!        write(*,*) 'iface tag', abs(ef2e(1,iface,ielem))

        ! Check if the face is a "periodic" face with tag equals to 8
        if (abs(ef2e(1,iface,ielem)) == 12) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

!          write(*,*)
!          write(*,*) 'node coordinate iface'
!          do i_vertex_iface = 1, nverticesperface 
!            write(*,*) 'node', i_vertex_iface
!            write(*,*) 'x', vx_iface(1,i_vertex_iface) 
!            write(*,*) 'y', vx_iface(2,i_vertex_iface) 
!            write(*,*) 'z', vx_iface(3,i_vertex_iface) 
!          end do

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

!              write(*,*)
!              write(*,*) 'node coordinate jface'
!              do i_vertex_jface = 1, nverticesperface 
!                write(*,*) 'node', i_vertex_jface
!                write(*,*) 'x', vx_jface(1,i_vertex_jface) 
!                write(*,*) 'y', vx_jface(2,i_vertex_jface) 
!                write(*,*) 'z', vx_jface(3,i_vertex_jface) 
!              end do
!              write(*,*)

              !stop
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
!                    write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                    write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
                    if (abs(vx_jface(i_coord,i_vertex_jface) - vx_iface(i_coord,i_vertex_iface)) .lt. diff_toll) then
                      same_coord(i_vertex_iface,i_coord) = 1
!                      write(*,*) 'i_vertex_iface', i_vertex_iface
!                      write(*,*) 'i_vertex_jface', i_vertex_jface
!                      write(*,*) 'vx_iface(i_coord,i_vertex_iface)', vx_iface(i_coord,i_vertex_iface)
!                      write(*,*) 'vx_jface(i_coord,i_vertex_jface)', vx_jface(i_coord,i_vertex_jface)
!                      write(*,*) 'HERE'
!                      write(*,*) 'i_coord', i_coord
                    endif
                  end do
                  
                  if (sum(same_coord(i_vertex_iface,:)) == 2) then
                    ! Found a vertex of the jface which has two coordinates
                    ! equal to two coordinates of the vertex on the iface
!                    write(*,*) 'Found node'
                    if (i_vertex_iface .gt. 1) then
                      ! Verify that the other possible partner nodes have the 
                      ! same two coordinates in common'
                      !write(*,*) i_vertex_iface
                      do i_check = 2, i_vertex_iface
                        check_common_coord_x3_1: do i_coord = 1, 3
                          if (same_coord(1,i_coord)-same_coord(i_check,i_coord) == 0) then
                            !write(*,*) i_p_elem
                            !write(*,*) 'node 1',i_coord, same_coord(1,i_coord)
                            !write(*,*) 'other nodes', i_coord, same_coord(i_check,i_coord)
                            cycle check_common_coord_x3_1
                          else
                            !write(*,*) same_coord(1,:)
                            !write(*,*) same_coord(i_check,:)
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
!                  write(*,*) 'Not the right face'
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

!                write(*,*) 'same_coord', same_coord(i_vertex_iface)

              end do search_x3_1
             
!              write(*,*) 'Match found' !, exiting from the companion_face loop'
!              write(*,*) 'ID vertices jface', ivtmp_jface

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
!                      write(*,*) 'ef2e 1', ef2e(1,iface,ielem)
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

!        write(*,*) 'iface tag', abs(ef2e(1,iface,ielem))

        ! Check if the face is a "periodic" face
        if (abs(ef2e(1,iface,ielem)) == 13) then

          ! Get the coordinate of the nodes of the iface
          do i_vertex_iface = 1, nverticesperface       
            vx_iface(1,i_vertex_iface) = vx_master(1,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(2,i_vertex_iface) = vx_master(2,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
            vx_iface(3,i_vertex_iface) = vx_master(3,periodic_face_data_x3(4+i_vertex_iface,i_p_elem))
!            write(*,*) 'ID vertices iface', periodic_face_data(3+i_vertex_iface,i_p_elem)
          end do

!          write(*,*)
!          write(*,*) 'node coordinate iface'
!          do i_vertex_iface = 1, nverticesperface 
!            write(*,*) 'node', i_vertex_iface
!            write(*,*) 'x', vx_iface(1,i_vertex_iface) 
!            write(*,*) 'y', vx_iface(2,i_vertex_iface) 
!            write(*,*) 'z', vx_iface(3,i_vertex_iface) 
!          end do

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
              search_x3_3: do i_vertex_iface = 1, nverticesperface
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
                  !write(*,*) 'No the right face'
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

!                write(*,*) 'same_coord', same_coord(i_vertex_iface)

              end do search_x3_3
             
!              write(*,*) 'Match found' !, exiting from the companion_face loop'
!              write(*,*) 'ID vertices jface', ivtmp_jface

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
                !if (match == nverticesperface .and. jelem /= ielem) then
                if (match == nverticesperface) then
                  ef2e(2,iface,ielem) = jelem

!                  write(*,*) 'ielem, jelem', ielem, jelem
               
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
!                      write(*,*) 'ef2e 1', ef2e(1,iface,ielem)
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

!      do i = 1, size(iae2v)
!       cnt_pack = cnt_pack + 1
!        write(93,*) iae2v(i)
!      end do

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

!      write(80,*) iae2v_tmp
!      stop

    else
      write(*,*) 'Unknown number of spatial dimension of the problem:', ndim
      write(*,*) 'Exting...'
      stop
    end if ! End if ndim == 3

    return
  end subroutine e2e_connectivity_aflr3

  !============================================================================
  
  !============================================================================

  subroutine create_ldg_flip_flop_sign()
    
    ! Load modules
    use collocationvariables
    use referencevariables
    use variables, only: ef2e 

    ! Nothing is implicitly defined
    implicit none
   
    integer :: low_elem, high_elem

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
  end subroutine data_partner_element_serial

  pure function WENO_Adjoining_Data(k_node,k_face)
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

  end function WENO_Adjoining_Data

  pure function Pencil_Coord(Ns,jdir,iface,i)

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

  end function

  !============================================================================

  subroutine WENO_Intrp_Face_Nodes()

    !  Preprocessing step to build physical locations needed for interpolation

    use referencevariables
    use interpolation
    use variables
    use SSWENOvariables
    use mpimod
    use petscvariables, only: xpetscWENO_partner, xlocpetscWENO_partner

    implicit none
    ! indices
    integer :: ielem, inode, jnode, jdir, ipen, face_id, shift, iface, kface
    integer :: gnode, ieq, iloc, kelem, knode, nodespershell
    integer :: i,j,k, petscsize, i_err
    integer :: extrnal_xi_cnt, extrnal_xi_sum
    integer, dimension(2)  :: faceLR

    real(wp), parameter                 :: sqrt5  = sqrt(5.0_wp)
    real(wp), parameter                 :: tol    = 1.0e-09_wp

    real(wp), dimension(3,nodesperedge) :: xgS
    real(wp), dimension(3)              :: xgL, xgR

    real(wp), dimension(3,8)            :: comp2phys_coeffs
    real(wp)                            :: tmp, t1,t2,t3, rnorm, rnorm_max
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
                   - xghst(:,iloc) + xg(:,inode,ielem)                   ! account for possibility of non-periodic domain.
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

  end subroutine WENO_Intrp_Face_Nodes

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

    use referencevariables, only: nvertices, nodesperface
    use variables, only: if2nq, ifacetag, ic2nh, nqface
    use unary_mod, only: qsorti
!   use variables, only: iaBv2Bv, jaBv2Bv

    !  nqface     = number of boundary faces with quads
    !  if2nq      = face-to-node pointers for each QUAD boundary face (4 nodes)
    !  ifacetag   = tag number that groups boundary faces (e.g., all similar faces)

    ! Nothing is implicitly defined
    implicit none

    integer :: i,j,k,L,m,n
    integer :: k1,k2,k3
    integer :: nI,nO
    integer :: jV,jE
    integer :: icnt0,icnt1,icnt2,icnt3
    integer :: iwrk
    integer, dimension(2)  :: itmp

    integer, parameter   :: bigN = 20

    integer, dimension(5)                :: stack
    integer, dimension(bigN)             :: shuffle, wrk_vec
    integer, dimension(bigN)             :: wrk_vec0, wrk_vec1, wrk_vec2, wrk_vec3

    integer                              :: n_bc_types, new_vert

    integer, dimension(2,4)              :: if2nq_off

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


    real(wp), dimension(ndim)             :: dx, vec1
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

  pure function ortho_projection_to_line3D(x00,x01,x)

    real(wp), dimension(3), intent(in) :: x00,x01,x

    real(wp), dimension(3)             :: xmin

    real(wp)                           :: t
    real(wp), dimension(3)             :: ortho_projection_to_line3D

    t = dot_product(x01-x,x01-x00) / dot_product(x01-x00,x01-x00) ;

    ortho_projection_to_line3D = t*x00 + (1.0_wp - t)*x01 ;

  end function ortho_projection_to_line3D

end module initgrid
