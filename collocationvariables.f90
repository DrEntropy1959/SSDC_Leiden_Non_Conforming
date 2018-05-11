module collocationvariables
  use precision_vars
  implicit none

  ! Penalty Coefficients
  real(wp)                     :: l00  = +0.05_wp  !  IP      Penalty
  real(wp)                     :: l01  = -0.5_wp   !  Viscous Penalties: -0.5_wp
  real(wp)                     :: l10  = -0.5_wp   !  LDG     Penalties: -0.5_wp
  real(wp)                     :: sfix =  0.0001_wp   !  entropy fix

  real(wp)                     :: alpha_ldg_flip_flop = 0.0_wp
  
  real(wp), allocatable, dimension(:) :: face_sign
  integer, allocatable, dimension(:,:) :: ldg_flip_flop_sign

  real(wp), allocatable, dimension(:)     :: x_LGL_pts_1D
  real(wp), allocatable, dimension(:,:)   :: qmat
  real(wp), allocatable, dimension(:,:)   :: dmat
  real(wp), allocatable, dimension(:,:,:) :: gradmat
  real(wp), allocatable, dimension(:)     :: pmat
  real(wp), allocatable, dimension(:)     :: pinv
  real(wp), allocatable, dimension(:)     :: pvol


  real(wp), allocatable, dimension(:)     :: p_surf

  real(wp) , allocatable, dimension(:)    :: x_gll, w_gll

  ! Gradient matrix
  integer,  allocatable, dimension(:)     :: iagrad
  integer,  allocatable, dimension(:,:)   :: jagrad
  real(wp), allocatable, dimension(:,:)   :: dagrad
  real(wp), allocatable, dimension(:,:)   :: qagrad
  real(wp), allocatable, dimension(:,:)   :: fbar
  integer :: nnzgrad
 
  integer, allocatable, dimension (:,:,:) :: lgradindices

  ! Filter matrix
  integer,  allocatable, dimension(:)     :: ia_Filter
  integer,  allocatable, dimension(:,:)   :: ja_Filter
  real(wp), allocatable, dimension(:,:)   :: aa_Filter
  real(wp), allocatable, dimension(:,:,:) :: mat_Filter
  real(wp), allocatable, dimension(:,:)   :: Filter
  integer :: nnz_Filter

  ! Interpolation matrix
  integer,  allocatable, dimension(:)     :: ia_Intrpltn
  integer,  allocatable, dimension(:,:)   :: ja_Intrpltn
  real(wp), allocatable, dimension(:,:)   :: aa_Intrpltn
  real(wp), allocatable, dimension(:,:,:) :: mat_Intrpltn
  real(wp), allocatable, dimension(:,:)   :: Int_F2S
  integer :: nnz_Intrpltn

  ! Extrapolation matrix
  integer,  allocatable, dimension(:)     :: ia_Extrpltn
  integer,  allocatable, dimension(:,:)   :: ja_Extrpltn
  real(wp), allocatable, dimension(:,:)   :: aa_Extrpltn
  real(wp), allocatable, dimension(:,:,:) :: mat_Extrpltn
  real(wp), allocatable, dimension(:,:)   :: Ext_S2F
  integer :: nnz_Extrpltn

  ! Extrapolation matrix P^{-1}Q = P^{-1} Delta F = P^{-1} Delta Int_Coeff_S2F f
  real(wp), allocatable, dimension(:,:)   :: Ext_SSSCE_S2F

  integer :: N_Soln_Pts, N_Flux_Pts
  real(wp), allocatable, dimension(:)   :: X_Soln_Pts, X_Flux_Pts
  real(wp), allocatable, dimension(:)   :: pmat_Soln , pmat_Flux
  real(wp), allocatable, dimension(:)   :: pinv_Soln , pinv_Flux
  real(wp), allocatable, dimension(:,:) :: qmat_Soln , qmat_Flux
  real(wp), allocatable, dimension(:,:) :: dmat_Soln , dmat_Flux

  real(wp), allocatable, dimension(:,:) :: Dif_F2S

  integer, allocatable, dimension(:,:)    :: elem_props

  !  grid info LGL
  integer :: n_LGL_1d_p0, n_LGL_2d_p0, n_LGL_3d_p0
  integer :: n_LGL_1d_p1, n_LGL_2d_p1, n_LGL_3d_p1
  integer :: n_LGL_1d_p2, n_LGL_2d_p2, n_LGL_3d_p2
  
  real(wp), allocatable, dimension(:)     :: x_LGL_1d_p0, x_LGL_1d_p1, x_LGL_1d_p2
  real(wp), allocatable, dimension(:)     :: w_LGL_1d_p0, w_LGL_1d_p1, w_LGL_1d_p2

  !  grid info Gau
  integer :: n_Gau_1d_p0, n_Gau_2d_p0, n_Gau_3d_p0
  integer :: n_Gau_1d_p1, n_Gau_2d_p1, n_Gau_3d_p1
  integer :: n_Gau_1d_p2, n_Gau_2d_p2, n_Gau_3d_p2

  real(wp), allocatable, dimension(:)     :: x_Gau_1d_p0, x_Gau_1d_p1, x_Gau_1d_p2
  real(wp), allocatable, dimension(:)     :: w_Gau_1d_p0, w_Gau_1d_p1, w_Gau_1d_p2

  integer :: n_Gau_Shell

  real(wp), allocatable, dimension(:)     ::   pmat_p0,   pmat_p1,    pmat_p2
  real(wp), allocatable, dimension(:)     ::   pinv_p0,   pinv_p1,    pinv_p2
  real(wp), allocatable, dimension(:,:)   ::   qmat_p0,   qmat_p1,    qmat_p2
  real(wp), allocatable, dimension(:,:)   ::   dmat_p0,   dmat_p1,    dmat_p2
  real(wp), allocatable, dimension(:)     ::   pvol_p0,   pvol_p1,    pvol_p2
  real(wp), allocatable, dimension(:)     :: p_surf_p0, p_surf_p1,  p_surf_p2

  real(wp), allocatable, dimension(:,:,:) :: gradmat_p0, gradmat_p1, gradmat_p2
  
  ! Gradient matrix
  integer,  allocatable, dimension(:)     :: iagrad_p0, iagrad_p1, iagrad_p2
  integer,  allocatable, dimension(:,:)   :: jagrad_p0, jagrad_p1, jagrad_p2
  real(wp), allocatable, dimension(:,:)   :: dagrad_p0, dagrad_p1, dagrad_p2
  real(wp), allocatable, dimension(:,:)   :: qagrad_p0, qagrad_p1, qagrad_p2
  integer :: nnzgrad_p0, nnzgrad_p1, nnzgrad_p2

  real(wp), allocatable, dimension(:,:,:,:) :: Prolong_LGL_2_Gau_1d
  real(wp), allocatable, dimension(:,:,:,:) :: Restrct_Gau_2_LGL_1d

  real(wp), allocatable, dimension(:,:)   :: Rot_Gau_p0_2_LGL_p0_1d
  real(wp), allocatable, dimension(:,:)   :: Rot_LGL_p0_2_Gau_p0_1d

  real(wp), allocatable, dimension(:,:)   :: Ext_LGL_p0_2_Gau_p1_1d
  real(wp), allocatable, dimension(:,:)   :: Int_Gau_p1_2_LGL_p0_1d

  real(wp), allocatable, dimension(:,:)   :: Rot_Gau_p1_2_LGL_p1_1d
  real(wp), allocatable, dimension(:,:)   :: Rot_LGL_p1_2_Gau_p1_1d

  real(wp), allocatable, dimension(:,:)   :: Ext_LGL_p1_2_Gau_p2_1d
  real(wp), allocatable, dimension(:,:)   :: Int_Gau_p2_2_LGL_p1_1d

  real(wp), allocatable, dimension(:,:)   :: Rot_Gau_p2_2_LGL_p2_1d
  real(wp), allocatable, dimension(:,:)   :: Rot_LGL_p2_2_Gau_p2_1d

  real(wp), allocatable, dimension(:,:)   :: LGL_Coarse_2_LGL_Fine_1d
  real(wp), allocatable, dimension(:,:)   :: LGL_Fine_2_LGL_Coarse_1d


end module collocationvariables
