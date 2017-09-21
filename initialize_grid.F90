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

end module initialize_grid
