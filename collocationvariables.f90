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

  real(wp), allocatable, dimension(:)     :: rcollocation
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

end module collocationvariables
