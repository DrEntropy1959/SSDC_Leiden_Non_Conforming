! This module contains the declaration of the working arrays used by the
! CSR sparsekit for constructing the element-wise Jacoibian matrix. 

module jacobian_matrix_implicit_ts_variables 

  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none

  integer, allocatable, dimension(:) :: iw_inviscid
  integer, allocatable, dimension(:) :: iw_viscous_1
  integer, allocatable, dimension(:) :: iw_viscous_2

end module jacobian_matrix_implicit_ts_variables
