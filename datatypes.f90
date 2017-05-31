module datatypes
  use precision_vars
  implicit none

  type bcarray
    character(60) :: name
    ! boundary type
    integer :: btype = 0
    ! element type
    integer :: eltype = 0
    ! array of global element indices
    integer, allocatable, dimension(:,:) :: belems
    ! array of global node indices
    integer, allocatable, dimension(:) :: bnodes
    ! number of vertices per boundary element
    integer :: nbvpe
    ! number of boundary elements
    integer :: nbelems = 0
    ! starting and ending
    integer :: ielstart, ielend
    ! volumetric element dimension
    integer :: mvdim
  end type bcarray

end module datatypes
