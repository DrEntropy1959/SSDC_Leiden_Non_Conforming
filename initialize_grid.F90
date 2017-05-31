module initialize_grid
  
  ! Load modules
  use precision_vars
  use iso_c_binding

  ! Nothing is implicitly defined
  implicit none

  ! Subroutines and functions in this module are usually private
  private
  
  ! Exceptions, i.e. public subroutines or functions
  public init_hex_aflr3
  public init_elem_type

  ! Global variables and parameters for this module 
  integer, allocatable, dimension(:,:), target :: hex_aflr3_faces
  integer, allocatable, dimension(:), target :: hex_aflr3_face_directions
  integer, parameter :: hex_aflr3_n_vertices_per_face = 4
  integer, parameter :: hex_aflr3_n_faces_per_elem = 6

  integer, pointer, dimension(:,:) :: elem_type_faces
  integer, pointer, dimension(:) :: elem_face_directions


  ! Interface to the C function that calls Metis
  interface
    
    integer(c_int) function calcMetisPartitions(ne, nv, nps, xadj, adj, ncommon, &
        epart, npart, nnz) &
        bind(c,name='calcMetisPartitions') 
      
      ! Load modules
      use iso_c_binding, only: c_ptr, c_double, c_int, c_char
      
      ! Nothing is implicitly defined
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

  !============================================================================
  
  !============================================================================
  ! init_hex_aflr3 - Sets the node face connectivity at the element level.
  
  subroutine init_hex_aflr3()
    
    ! Load modules
    use referencevariables, only: nverticesperface, nfacesperelem
    
    ! Nothing is implicitly defined
    implicit none

    continue

    
!
!                        7-----------------4
!                       /.                /|
!                      / .               / |
!                     /  .              /  |
!                    /   .             /   |
!                   /    .            /    |
!                  /     .           /     |
!                 /      .          /      |
!                /       .         /       |
!               5-----------------2        |
!               |        3........|........8
!               |       .         |       /
!               |      .          |      /
!               |     .           |     /  zeta ^      
!               |    .            |    /        |    / eta
!               |   .             |   /         |   /
!               |  .              |  /          |  /
!               | .               | /           | /
!               |.                |/            |/
!               1-----------------6             ----------> csi


    ! Number of vertices per face
    nverticesperface = 4
    
    ! Number of faces per element
    nfacesperelem = 8
    
    ! Local vertices used to construct each face
    allocate(hex_aflr3_faces(nverticesperface,nfacesperelem))

    hex_aflr3_faces(:,1) = (/ 1, 3, 8, 6 /)
    hex_aflr3_faces(:,2) = (/ 1, 6, 2, 5 /)
    hex_aflr3_faces(:,3) = (/ 6, 8, 4, 2 /)
    hex_aflr3_faces(:,4) = (/ 8, 3, 7, 4 /)
    hex_aflr3_faces(:,5) = (/ 1, 5, 7, 3 /)
    hex_aflr3_faces(:,6) = (/ 5, 2, 4, 7 /)

    ! Outward facing normal in each computational direction for each face
    allocate(hex_aflr3_face_directions(nfacesperelem))
    
    hex_aflr3_face_directions(1) = -3
    hex_aflr3_face_directions(2) = -2
    hex_aflr3_face_directions(3) = +1
    hex_aflr3_face_directions(4) = +2
    hex_aflr3_face_directions(5) = -1
    hex_aflr3_face_directions(6) = +3

    return
  end subroutine init_hex_aflr3

  !============================================================================
  
  !============================================================================
  ! init_elem_type - Sets up the properties of different element type.

  subroutine init_elem_type()
    
    ! Load modules
    use referencevariables
    use variables, only: facenormalcoordinate

    ! Nothing is implicitly defined
    implicit none
    integer :: i

    continue

    ! We are limited to tensor product elements at this time
    if (ndim == 3) then
      
      ! Set local node-face connectivity pointer 
      elem_type_faces => hex_aflr3_faces

      ! Set local outward faceing normal
      elem_face_directions => hex_aflr3_face_directions

      ! Number of vertices per face
      nverticesperface = hex_aflr3_n_vertices_per_face
      
      ! Number of face per element
      nfacesperelem = hex_aflr3_n_faces_per_elem
      
      ! Assign outward facing normal 
      allocate(facenormalcoordinate(nfacesperelem))
      do i = 1, nfacesperelem
        facenormalcoordinate(i) = hex_aflr3_face_directions(i)
      enddo
    else
      write(*,*) 'Error: unsupported dimension'
      Write(*,*) 'Exting...'
      stop
    end if

    return
  end subroutine init_elem_type

  !============================================================================

  !============================================================================

!  subroutine elem_to_elem_connectivity()
!    
!    ! Load modules
!    use referencevariables
!    use variables, only: boundaryelems, ef2e,  &
!                         iae2v, jae2v, nnze2v
!!                        iae2e, jae2e, nnze2e,
!
!    ! Nothing is implicitly defined
!    implicit none
!
!    integer :: i,j,k,izone,jzone,kzone
!    integer :: ii,jj,kk,icount,ielem
!    integer :: iface, jface
!
!    integer, allocatable :: v2e(:,:), iav2e(:)
!    integer, allocatable, dimension(:) :: ivtmp1, ivtmp2
!    integer :: nnzv2e
!    integer :: connectioncount, vertexcount
!
!    continue
!
!
!    !  routine called from myprocid == 0 ;  i.e. master works on entire grid
!    !
!    !  nelems     :    elements
!    !  nelemzones :    zones = 2,  1 == interior; 2 == boundary 
!    !
!    !                   Dim,    Dim,         Dim
!    !  ef2e       :    ( 2 ,nfaceperelem, nelements) 
!    !             :  Two situation occur.  The face is either an 
!    !                  (Interior face 
!    !                      :  (1,j,k) = face ID of the adjoining element
!    !                      :  (2,j,k) = Connected to Element 
!    !                  (Boundary face 
!    !                      :  (1,j,k) = Set to -11 
!    !                      :  (2,j,k) = Boundary Face number
!    !
!    !         (Note:  redimensioned (3,:,:) to account for processor info)
!    !
!    ! iae2v,jae2v     :    Which vertices belong to each element
!
!
!    ! calculate total elements, including boundary elements
!
!    write(*,*) 'nelems     = ', nelems
!    write(*,*) 'nelemzones = ', nelemzones
!    write(*,*) 'nvertices  = ', nvertices
!
!    ! allocate rows for vertex to element connectivity
!    allocate(iav2e(nvertices+1))
!    iav2e = 0
!
!    ! loop through elements. Count the number of accesses for each vertex
!    do izone = 1, nelemzones
!      do k = boundaryelems(izone)%ielstart, boundaryelems(izone)%ielend
!        do j = 1, boundaryelems(izone)%nbvpe
!          iav2e(boundaryelems(izone)%belems(j,k)) = &
!          iav2e(boundaryelems(izone)%belems(j,k)) + 1 
!        end do
!      end do
!    end do
!
!    ! allocate nonzeroes
!    nnzv2e = sum(iav2e)
!    allocate(v2e(2,nnzv2e))
!
!    ! Count backwards such that row structure is correct upon populating v2e
!    iav2e(1) = iav2e(1) + 1
!    do i = 2, nvertices+1
!      iav2e(i) = iav2e(i) + iav2e(i-1)
!    end do
!    do izone = 1, nelemzones
!      do k = boundaryelems(izone)%ielstart, boundaryelems(izone)%ielend
!        do j = 1, boundaryelems(izone)%nbvpe
!          iav2e(boundaryelems(izone)%belems(j,k)) = &
!          iav2e(boundaryelems(izone)%belems(j,k)) - 1
!            v2e(:,iav2e(boundaryelems(izone)%belems(j,k))) = (/ k, izone /) 
!        end do
!      end do
!    end do
!
!    allocate(ivtmp1(nverticesperface),ivtmp2(nverticesperface))
!    ! calculate element-to-element connectivity using shared nodes
!    allocate(ef2e(2,2*ndim,1:nelems))
!    ef2e = 0
!    nnze2v = nverticesperelem*nelems
!    allocate(iae2v(nelems+1))
!    iae2v = 0
!    allocate(jae2v(nnze2v))
!    jae2v = 0
!
!    connectioncount = 0
!    k = 0
!    vertexcount = 1
!    do kzone = 1, nvolumesections
!      izone = isectionvolume(kzone)
!      do ielem = boundaryelems(izone)%ielstart, boundaryelems(izone)%ielend
!        k = k+1
!        iae2v(k) = vertexcount
!        do j = 1, nverticesperelem
!          jae2v(vertexcount) = boundaryelems(izone)%belems(j,k)
!          vertexcount = vertexcount+1
!        end do
!        faceloop: do iface = 1, nfacesperelem
!          do j = 1, nverticesperface
!            ! populate vertices on face
!            ivtmp1(j) = boundaryelems(izone)%belems(eltypfaces(j,iface),k)
!          end do
!          ivtmp1 = isort(ivtmp1,nverticesperface)
!          ! search elements connected to each vertex until a match is found
!          do j = 1, nverticesperface
!            ! search elements connected to vertices and count number of shared vertices
!            do i = iav2e(ivtmp1(j)), iav2e(ivtmp1(j)+1) - 1
!              jzone = v2e(2,i)
!              kk = v2e(1,i)
!              ! don't check self!
!              if (kk == k) cycle
!              icount = 0
!              ! check for every vertex
!              do ii = 1,nverticesperface
!                do jj = 1, boundaryelems(jzone)%nbvpe
!                  if( ivtmp1(ii) == boundaryelems(jzone)%belems(jj,kk) ) then
!                    icount = icount + 1
!                    exit
!                  end if
!                end do
!              end do
!              ! check to see if nvertices per face match
!              if (icount == nverticesperface) then
!                ! element is connected to face
!                ef2e(2,iface,k) = kk
!                if (jzone /= izone) then
!                  ! this is a boundary face
!                  ef2e(1,iface,k) = -boundaryelems(jzone)%btype
!                else
!                  ! this is not a boundary face
!                  connectioncount = connectioncount + 1
!                  ! determine which face of kk I am connected to
!                  do jface = 1, nfacesperelem
!                    icount = 0
!                    do jj = 1, nverticesperface
!                      ivtmp2(jj) = boundaryelems(jzone)%belems(eltypfaces(jj,jface),kk)
!                    end do
!                    ivtmp2 = isort(ivtmp2,nverticesperface)
!                    if (all(ivtmp1==ivtmp2)) then
!                      ! this is the right face
!                      ef2e(1,iface,k) = jface
!                      exit
!                    end if
!                  end do
!                end if
!                cycle faceloop
!              end if
!            end do
!          end do
!          ! Search should always find a connection but didnt. Somethings wrong
!          write(*,*) 'face search failed' ; write(*,*) izone, k, iface ; write(*,*) ivtmp1, ivtmp2
!          stop
!        end do faceloop
!
!      end do
!    end do
!
!    iae2v(nelems+1) = vertexcount
!
!!   nnze2e = connectioncount
!!   !   local work (i.e. in initgrid)
!!   allocate(iae2e(nelems+1))
!!   iae2e = 0
!!   allocate(jae2e(nnze2e))
!!   jae2e = 0
!
!!   ii = 1
!!   do ielem = 1, nelems
!!     iae2e(ielem) = ii
!!     do j = 1, nfacesperelem
!!       ! only add element connections to other elements
!!       if ( ef2e(1,j,ielem) > 0 ) then
!!         jae2e(ii) = ef2e(2,j,ielem)
!!         ii = ii+1
!!       end if
!!     end do
!!   end do
!!   iae2e(nelems+1) = ii
!
!    deallocate(ivtmp1,ivtmp2)
!
!    deallocate(iav2e,v2e)
!
!  end subroutine elem_to_elem_connectivity


end module initialize_grid
