module referencevariables
  
  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none

  integer :: nfacesperelem = 2
  integer :: ndim = 1
  integer :: nedges = 0
  integer :: nfaces = 0
  integer :: nelems = 0
  integer :: ndof = 0
  integer :: nnodes = 0
  integer :: nvertices = 0
  integer :: nequations = 1
  integer :: npoly = 1
  integer :: npoly_max = 1
  integer :: npoly_DeltaF = 0
  integer :: ndofperelem = 0
  integer :: nodesperelem = 0
  integer :: nodesperedge = 0
  integer :: nodesperface = 0
  integer :: nphysdim = 1
  integer :: nconstitutive = 0
  integer :: nverticesperelem = 0
  integer :: nvolumesections = 0
  integer, allocatable :: isectionvolume(:)
  integer :: nelemzones = 0
  integer :: nverticesperface = 0
  integer :: nghost = 0
  integer :: nghost_elem = 0
  integer :: nghost_Gau_shell = 0

  integer :: ihelems(2)
  integer :: nprocs
  integer :: myprocid

  logical :: periodic = .false.
  real(wp), dimension(3)  :: periodic_distance = 0.0_wp

end module referencevariables
