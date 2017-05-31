module CSRlocalvariables
  
  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none


  ! CSR used for building the Jacobian matrix
  ! -----------------------------------------
  ! CSR of the first and second order differentiation matrices for all the
  ! elements owned by a processor

  ! Time terms
  integer,  allocatable, dimension(:) :: ia_0
  integer,  allocatable, dimension(:) :: ja_0
  real(wp), allocatable, dimension(:) ::  a_0
  integer,  allocatable, dimension(:) :: ka_0
  
  ! x1
  integer,  allocatable, dimension(:) :: ia_x1
  integer,  allocatable, dimension(:) :: ja_x1
  real(wp), allocatable, dimension(:) ::  a_x1
  integer,  allocatable, dimension(:) :: ka_x1
  
  ! x2
  integer,  allocatable, dimension(:) :: ia_x2
  integer,  allocatable, dimension(:) :: ja_x2
  real(wp), allocatable, dimension(:) ::  a_x2  
  integer,  allocatable, dimension(:) :: ka_x2

  ! x3
  integer,  allocatable, dimension(:) :: ia_x3
  integer,  allocatable, dimension(:) :: ja_x3
  real(wp), allocatable, dimension(:) ::  a_x3
  integer,  allocatable, dimension(:) :: ka_x3

  ! x11
  integer,  allocatable, dimension(:) :: ia_x11
  integer,  allocatable, dimension(:) :: ja_x11
  real(wp), allocatable, dimension(:) ::  a_x11
  integer,  allocatable, dimension(:) :: ka_x11

  ! x22
  integer,  allocatable, dimension(:) :: ia_x22
  integer,  allocatable, dimension(:) :: ja_x22
  real(wp), allocatable, dimension(:) ::  a_x22
  integer,  allocatable, dimension(:) :: ka_x22

  ! x33
  integer,  allocatable, dimension(:) :: ia_x33
  integer,  allocatable, dimension(:) :: ja_x33
  real(wp), allocatable, dimension(:) ::  a_x33
  integer,  allocatable, dimension(:) :: ka_x33

  ! x12
  integer,  allocatable, dimension(:) :: ia_x12
  integer,  allocatable, dimension(:) :: ja_x12
  real(wp), allocatable, dimension(:) ::  a_x12
  integer,  allocatable, dimension(:) :: ka_x12

  ! x13
  integer,  allocatable, dimension(:) :: ia_x13
  integer,  allocatable, dimension(:) :: ja_x13
  real(wp), allocatable, dimension(:) ::  a_x13
  integer,  allocatable, dimension(:) :: ka_x13

  ! x23
  integer,  allocatable, dimension(:) :: ia_x23
  integer,  allocatable, dimension(:) :: ja_x23
  real(wp), allocatable, dimension(:) ::  a_x23
  integer,  allocatable, dimension(:) :: ka_x23

  ! Inviscid penalty
  integer,  allocatable, dimension(:) :: ia_penI
  integer,  allocatable, dimension(:) :: ja_penI
  integer,  allocatable, dimension(:) :: ka_penI

  real(wp), allocatable, dimension(:) ::  a_penI

  integer,  allocatable, dimension(:) :: ia_penI_proc
  integer,  allocatable, dimension(:) :: ja_penI_proc
  integer,  allocatable, dimension(:) :: ka_penI_proc
  integer,  allocatable, dimension(:) :: la_penI_proc


  ! Viscous penalty
  integer,  allocatable, dimension(:) :: ia_penV
  integer,  allocatable, dimension(:) :: ja_penV
  real(wp), allocatable, dimension(:) ::  a_penV
  integer,  allocatable, dimension(:) :: ka_penV
  
  integer,  allocatable, dimension(:) :: ia_penV_proc
  integer,  allocatable, dimension(:) :: ja_penV_proc
  integer,  allocatable, dimension(:) :: ka_penV_proc
  integer,  allocatable, dimension(:) :: la_penV_proc

  real(wp), allocatable, dimension(:) ::  a_penV_proc

  ! "Sum" of all contributions 
  integer,  allocatable, dimension(:) :: iaS
  integer,  allocatable, dimension(:) :: jaS

  ! Number of rows in the processor
  integer :: nprows

  ! Number of global rows
  integer :: ngrows


  ! CSR of the first and second order differentiation matrices for one element
  ! --------------------------------------------------------------------------
  ! Time terms
  integer,  allocatable, dimension(:) :: ia_0_elem
  integer,  allocatable, dimension(:) :: ja_0_elem
  real(wp), allocatable, dimension(:) ::  a_0_elem
  
  ! x1
  integer, allocatable, dimension(:)  :: ia_x1_elem 
  integer, allocatable, dimension(:)  :: ja_x1_elem
  real(wp), allocatable, dimension(:) :: a_x1_elem
  
  ! x2
  integer, allocatable, dimension(:)  :: ia_x2_elem
  integer, allocatable, dimension(:)  :: ja_x2_elem
  real(wp), allocatable, dimension(:) :: a_x2_elem
  
  ! x3
  integer, allocatable, dimension(:)  :: ia_x3_elem
  integer, allocatable, dimension(:)  :: ja_x3_elem
  real(wp), allocatable, dimension(:) :: a_x3_elem

  ! x11
  integer, allocatable, dimension(:)  :: ia_x11_elem 
  integer, allocatable, dimension(:)  :: ja_x11_elem
  real(wp), allocatable, dimension(:) :: a_x11_elem

  ! x22
  integer, allocatable, dimension(:)  :: ia_x22_elem 
  integer, allocatable, dimension(:)  :: ja_x22_elem
  real(wp), allocatable, dimension(:) :: a_x22_elem

  ! x33
  integer, allocatable, dimension(:)  :: ia_x33_elem 
  integer, allocatable, dimension(:)  :: ja_x33_elem
  real(wp), allocatable, dimension(:) :: a_x33_elem

  ! x12
  integer, allocatable, dimension(:)  :: ia_x12_elem 
  integer, allocatable, dimension(:)  :: ja_x12_elem
  real(wp), allocatable, dimension(:) :: a_x12_elem

  ! x13
  integer, allocatable, dimension(:)  :: ia_x13_elem 
  integer, allocatable, dimension(:)  :: ja_x13_elem
  real(wp), allocatable, dimension(:) :: a_x13_elem

  ! x23
  integer, allocatable, dimension(:)  :: ia_x23_elem 
  integer, allocatable, dimension(:)  :: ja_x23_elem
  real(wp), allocatable, dimension(:) :: a_x23_elem

  real(wp), allocatable, dimension(:,:,:) :: inviscid_flux_jacobian_elem
  real(wp), allocatable, dimension(:,:,:) :: dwdu_elem
  real(wp), allocatable, dimension(:,:,:) :: hatc_elem
  real(wp), allocatable, dimension(:,:,:) :: dhatcdu_gradwj_elem

  integer,  allocatable, dimension(:)     :: ia_elem
  integer,  allocatable, dimension(:)     :: ja_elem
  real(wp), allocatable, dimension(:,:,:) :: dfdu_a_elem


!  Try something new  MHC

  integer,  allocatable, dimension(:)     :: ia_W1_elem
  integer,  allocatable, dimension(:)     :: ja_W1_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_W1_elem
  
  integer,  allocatable, dimension(:)     :: ia_W2_elem
  integer,  allocatable, dimension(:)     :: ja_W2_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_W2_elem
  
  integer,  allocatable, dimension(:)     :: ia_W3_elem
  integer,  allocatable, dimension(:)     :: ja_W3_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_W3_elem
  
  integer,  allocatable, dimension(:)     :: ia_W12_elem
  integer,  allocatable, dimension(:)     :: ja_W12_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_W12_elem

  integer,  allocatable, dimension(:)     :: ia_W123_elem
  integer,  allocatable, dimension(:)     :: ja_W123_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_W123_elem

  integer,  allocatable, dimension(:)     :: ia_dFvdU1_elem
  integer,  allocatable, dimension(:)     :: ja_dFvdU1_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_dFvdU1_elem

  integer,  allocatable, dimension(:)     :: ia_dFvdU2_elem
  integer,  allocatable, dimension(:)     :: ja_dFvdU2_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_dFvdU2_elem

  integer,  allocatable, dimension(:)     :: ia_dFvdUx_elem
  integer,  allocatable, dimension(:)     :: ja_dFvdUx_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_dFvdUx_elem

  integer,  allocatable, dimension(:)     :: ia_dFvdUy_elem
  integer,  allocatable, dimension(:)     :: ja_dFvdUy_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_dFvdUy_elem

  integer,  allocatable, dimension(:)     :: ia_dFvdUz_elem
  integer,  allocatable, dimension(:)     :: ja_dFvdUz_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_dFvdUz_elem

  integer,  allocatable, dimension(:)     :: ia_Containerxy
  integer,  allocatable, dimension(:)     :: ja_Containerxy
  real(wp), allocatable, dimension(:,:,:) ::  a_Containerxy

  integer,  allocatable, dimension(:)     :: ia_Containerxyz
  integer,  allocatable, dimension(:)     :: ja_Containerxyz
  real(wp), allocatable, dimension(:,:,:) ::  a_Containerxyz

  integer,  allocatable, dimension(:)     :: ia_diag_tmp
  integer,  allocatable, dimension(:)     :: ja_diag_tmp
  real(wp), allocatable, dimension(:,:,:) ::  a_diag_tmp

  integer,  allocatable, dimension(:,:)     :: ia_dFvdU_elem
  integer,  allocatable, dimension(:,:)     :: ja_dFvdU_elem
  real(wp), allocatable, dimension(:,:,:,:) ::  a_dFvdU_elem

!  Try something new  MHC

  ! "Operator" CSR for one element
  ! ------------------------------
  ! Diagonal CSR for matrix-matrix multiply
  integer, allocatable, dimension(:) :: ia_diag_elem
  integer, allocatable, dimension(:) :: ja_diag_elem

  ! Contribution 1 (first derivative)
  integer,  allocatable, dimension(:)     :: ia_1_matmul_elem
  integer,  allocatable, dimension(:)     :: ja_1_matmul_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_1_matmul_elem
  
  ! Contribution 2 (first derivative)
  integer,  allocatable, dimension(:)     :: ia_2_matmul_elem
  integer,  allocatable, dimension(:)     :: ja_2_matmul_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_2_matmul_elem

  ! Laplacian terms
  integer,  allocatable, dimension(:)     :: ia_lap_matmul_elem
  integer,  allocatable, dimension(:)     :: ja_lap_matmul_elem
  real(wp), allocatable, dimension(:,:,:) ::  a_lap_matmul_elem
  
  ! Cross-terms 
  integer,  allocatable, dimension(:) :: ia_ct_matmul_elem
  integer,  allocatable, dimension(:) :: ja_ct_matmul_elem
  real(wp), allocatable, dimension(:,:,:) :: a_ct_matmul_elem

  ! CSR of the first and second order differentiation matrices for all the
  ! elements owned by a processor
  ! ---------------------------------------------------------------------
  ! Time terms
  integer,  allocatable, dimension(:) :: ia_0_matmul
  integer,  allocatable, dimension(:) :: ja_0_matmul
  integer,  allocatable, dimension(:) :: ka_0_matmul

  ! x1
  integer,  allocatable, dimension(:) :: ia_x1_matmul
  integer,  allocatable, dimension(:) :: ja_x1_matmul
  integer,  allocatable, dimension(:) :: ka_x1_matmul
  
  ! x2
  integer,  allocatable, dimension(:) :: ia_x2_matmul
  integer,  allocatable, dimension(:) :: ja_x2_matmul
  integer,  allocatable, dimension(:) :: ka_x2_matmul

  ! x3
  integer,  allocatable, dimension(:) :: ia_x3_matmul
  integer,  allocatable, dimension(:) :: ja_x3_matmul
  integer,  allocatable, dimension(:) :: ka_x3_matmul

  ! x11
  integer,  allocatable, dimension(:) :: ia_x11_matmul
  integer,  allocatable, dimension(:) :: ja_x11_matmul
  integer,  allocatable, dimension(:) :: ka_x11_matmul
  
  ! x22
  integer,  allocatable, dimension(:) :: ia_x22_matmul
  integer,  allocatable, dimension(:) :: ja_x22_matmul
  integer,  allocatable, dimension(:) :: ka_x22_matmul
  
  ! x33
  integer,  allocatable, dimension(:) :: ia_x33_matmul
  integer,  allocatable, dimension(:) :: ja_x33_matmul
  integer,  allocatable, dimension(:) :: ka_x33_matmul

  ! x12
  integer,  allocatable, dimension(:) :: ia_x12_matmul
  integer,  allocatable, dimension(:) :: ja_x12_matmul
  integer,  allocatable, dimension(:) :: ka_x12_matmul
  
  ! x13
  integer,  allocatable, dimension(:) :: ia_x13_matmul
  integer,  allocatable, dimension(:) :: ja_x13_matmul
  integer,  allocatable, dimension(:) :: ka_x13_matmul

  ! x21
  integer,  allocatable, dimension(:) :: ia_x21_matmul
  integer,  allocatable, dimension(:) :: ja_x21_matmul
  integer,  allocatable, dimension(:) :: ka_x21_matmul
  
  ! x23
  integer,  allocatable, dimension(:) :: ia_x23_matmul
  integer,  allocatable, dimension(:) :: ja_x23_matmul
  integer,  allocatable, dimension(:) :: ka_x23_matmul

  ! x31
  integer,  allocatable, dimension(:) :: ia_x31_matmul
  integer,  allocatable, dimension(:) :: ja_x31_matmul
  integer,  allocatable, dimension(:) :: ka_x31_matmul

  ! x32
  integer,  allocatable, dimension(:) :: ia_x32_matmul
  integer,  allocatable, dimension(:) :: ja_x32_matmul
  integer,  allocatable, dimension(:) :: ka_x32_matmul


end module CSRlocalvariables

