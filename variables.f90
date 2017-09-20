module variables
  use precision_vars
  use datatypes
  implicit none

  real(wp), allocatable, dimension(:,:,:) :: xg
  real(wp), allocatable, dimension(:,:,:) :: ug
  real(wp), allocatable, dimension(:,:,:) :: vg
  real(wp), allocatable, dimension(:,:,:) :: wg

  real(wp), allocatable, dimension(:,:,:) :: xgWENO_self
  real(wp), allocatable, dimension(:,:,:) :: xgWENO_partner
  real(wp), allocatable, dimension(:,:,:) :: XIWENO_partner
  real(wp), allocatable, dimension(:,:,:) :: ugWENO_partner

  real(wp), allocatable, dimension(:,:) :: xghst
  real(wp), allocatable, dimension(:,:) :: xghstWENO_partner
  real(wp), allocatable, dimension(:,:) :: ughst
  real(wp), allocatable, dimension(:,:) :: ughstWENO
  real(wp), allocatable, dimension(:,:) :: ughstWENO_partner

  real(wp), allocatable, dimension(:,:) :: uelemghst
  real(wp), allocatable, dimension(:,:) :: velemghst
  real(wp), allocatable, dimension(:,:) :: welemghst
  real(wp), allocatable, dimension(:,:,:) :: r_x_ghst
  
  real(wp), allocatable, dimension(:,:,:,:) :: x_r
  real(wp), allocatable, dimension(:,:,:,:) :: r_x
  real(wp), allocatable, dimension(:,:)     :: Jx_r
  real(wp), allocatable, dimension(:)       :: dx_min_elem
  
  real(wp), allocatable, dimension(:,:,:) :: area_wall_faces

  real(wp), allocatable, dimension(:,:,:,:) :: divf
  real(wp), allocatable, dimension(:,:,:)   :: divf_S
  real(wp), allocatable, dimension(:,:,:,:) :: fg, fvg
  real(wp), allocatable, dimension(:,:,:,:) :: phig
  real(wp), allocatable, dimension(:,:,:,:) :: phig_err
  real(wp), allocatable, dimension(:,:,:)   :: phighst
  real(wp), allocatable, dimension(:,:)     :: mut
  real(wp), allocatable, dimension(:,:)     :: chig
  real(wp), allocatable, dimension(:,:,:)   :: omega
  real(wp), allocatable, dimension(:,:,:)   :: gsat

  real(wp), allocatable, dimension(:,:,:)   :: uold
  real(wp), allocatable, dimension(:,:,:)   :: uhat

  ! Variables used on the mortar interfaces 

  real(wp), allocatable, dimension(:,:)       ::  xgghst_Gau_Shell
  real(wp), allocatable, dimension(:,:)       ::  fvghst_Gau_Shell         ! (nq,nShell)

  real(wp), allocatable, dimension(:,:,:)     :: xg_Gau_Shell              ! (nq,nShell,nelem)

  real(wp), allocatable, dimension(:,:)       :: Jx_r_Gau_Shell            ! (nGau_Shell,nelem)

  real(wp), allocatable, dimension(:,:,:)     :: wg_Gau_Shell              ! (nq,nShell,nelem)
  real(wp), allocatable, dimension(:,:,:)     :: wg_Gau_Shell_tmp          ! (nq,nShell,nelem)
  real(wp), allocatable, dimension(:,:)       :: wgghst_Gau_Shell

  real(wp), allocatable, dimension(:,:,:,:)   :: phig_Gau_Shell            ! (nq,nd,nShell,nelem)
  real(wp), allocatable, dimension(:,:,:,:)   :: phig_Gau_Shell_tmp        ! (nq,nd,nShell,nelem)
  real(wp), allocatable, dimension(:,:,:)     ::  fvg_Gau_Shell            ! (nq,nShell,nelem)
  real(wp), allocatable, dimension(:,:,:)     ::  fvg_Gau_Shell_tmp        ! (nq,nShell,nelem)

  real(wp), allocatable, dimension(:,:,:)     :: phighst_Gau_Shell         ! (nq,nd,nShell)

  integer,  allocatable, dimension(:,:,:)     :: efn2efn_LGL               ! ( 4,nshell,nelem)
  integer,  allocatable, dimension(:,:,:)     :: efn2efn_Gau               ! ( 4,nshell,nelem)
  integer,  allocatable, dimension(:,:)       :: kfacenodes_Gau
  integer,  allocatable, dimension(:)         :: ifacenodes_Gau
  real(wp), allocatable, dimension(:,:,:)     :: facenodenormal_LGL_shell
  real(wp), allocatable, dimension(:,:,:)     :: facenodenormal_Gau_shell

  ! This variable are used with the low-storage-RK
  real(wp), allocatable, dimension(:,:,:)   :: du
  real(wp), allocatable, dimension(:,:,:)   :: dudt
  real(wp), allocatable, dimension(:,:)     :: dudt_S

  ! These variables are used with the IMEX-RK
  real(wp), allocatable, dimension(:,:,:)   :: uexp
  real(wp), allocatable, dimension(:,:,:)   :: non_lin_res
  real(wp), allocatable, dimension(:,:,:,:) :: Fexp
  real(wp), allocatable, dimension(:,:,:,:) :: Fimp

  ! Code BC names and IDs
  character(120), allocatable, dimension(:) :: code_bc_names

  ! AFLR3 variables
  integer :: nqface
  integer,  allocatable, dimension(:,:)   :: if2nq
  integer,  allocatable, dimension(:)     :: ifacetag
  integer,  allocatable, dimension(:,:)   :: ic2nh
  real(wp), allocatable, dimension(:,:)   :: vx_master
  integer,  allocatable, dimension(:,:)   :: aflr3_bc_ids
  integer,  allocatable, dimension(:,:)   :: periodic_face_data_x1, &
                                           & periodic_face_data_x2, &
                                           & periodic_face_data_x3

  integer                                 :: n_periodic_faces_x1, &
                                           & n_periodic_faces_x2, &
                                           & n_periodic_faces_x3

  integer,  allocatable, dimension(:,:)   :: periodic_elem_face_ids_x1, &
                                           & periodic_elem_face_ids_x2, &
                                           & periodic_elem_face_ids_x3

  integer,  allocatable, dimension(:,:)   :: wall_face_data
  integer                                 :: n_wall_faces
  integer,  allocatable, dimension(:,:)   :: wall_elem_face_ids


  real(wp), allocatable, dimension(:,:) :: vx
  integer, allocatable, dimension(:,:) :: e2v
  integer, allocatable, dimension(:) :: bctypes
  integer, allocatable, dimension(:,:,:) :: ef2e
  integer, allocatable, dimension(:) :: facenormalcoordinate
! integer, allocatable, dimension(:) :: iae2e
! integer, allocatable, dimension(:) :: jae2e
! integer :: nnze2e
  integer, allocatable, dimension(:) :: iae2v, iae2v_tmp
  integer, allocatable, dimension(:) :: jae2v, jae2v_tmp
  integer :: nnze2v
  integer, allocatable, dimension(:) :: elempart
  integer, allocatable, dimension(:) :: jelems
  integer, allocatable, dimension(:,:,:) :: efn2efn
  integer, allocatable, dimension(:,:) :: kfacenodes
  integer, allocatable, dimension(:)   :: ifacenodes
  integer, allocatable, dimension(:,:) :: kfacenodesWENO
  integer, allocatable, dimension(:)   :: ifacenodesWENO
  real(wp), allocatable, dimension(:,:,:) :: facenodenormal

  type(bcarray), allocatable, dimension(:) :: boundaryelems


  ! Time-averaged quantities
  real(wp), allocatable, dimension(:,:,:) :: mean_vg
  real(wp), allocatable, dimension(:,:,:) :: time_ave_prod_vel_comp
  real(wp), allocatable, dimension(:,:,:) :: reynolds_stress

  ! Entropy
  real(wp), allocatable, dimension(:,:) :: specific_entropy

  ! Gradient of the entropy variables in computational space
  real(wp), allocatable, dimension(:,:,:,:) :: grad_w_jacobian

  ! Kinetic energy
  real(wp), dimension(2) :: kinetic_energy

  ! Enstrophy
  real(wp) :: enstrophy

  ! Time derivative of the kinetic energy
  real(wp) :: dkinetic_energy_dt

end module variables
