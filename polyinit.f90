module polyinit
  use precision_vars
  implicit none

  private

  public rmapInit

contains

  function rmapInit(polyorder,ndimensions)
    ! This function initializes the polynomial
    ! approximation given a polynomial order
    ! and dimension of the grid.
    use referencevariables
    use initcollocation
    use collocationvariables
    use SSWENO_Routines
    use non_conforming, only: Rotate_xione_2_xitwo_and_back

    implicit none
    ! polynomial order
    integer, intent(in) :: polyorder
    integer, intent(in) :: ndimensions
    integer :: rmapInit
    integer :: n_LGL_1d, n_Gau_1d
    integer :: diff_poly
    integer :: i
    real(wp), dimension(:), allocatable  :: x_LGL_1d, w_LGL_1d, x_Gau_1d, w_Gau_1d

    rmapInit = 0

    ndim = ndimensions                    ! set dimension
    nfacesperelem = 2 * ndim              ! set number of faces per element

    npoly = polyorder                     ! set polynomial order

    nodesperedge = (npoly+1)              ! each element direction has npoly+1 nodes
    nodesperelem = nodesperedge**ndim     ! number of nodes per element
    nodesperface = nodesperedge**(ndim-1) ! number of nodes on each face

    n_LGL_1d_p0 = npoly + 1 
    n_LGL_2d_p0 = n_LGL_1d_p0**2
    n_LGL_3d_p0 = n_LGL_1d_p0**3

    n_LGL_1d_p1 = npoly + 2 
    n_LGL_2d_p1 = n_LGL_1d_p1**2
    n_LGL_3d_p1 = n_LGL_1d_p1**3

    n_LGL_1d_p2 = npoly + 3 
    n_LGL_2d_p2 = n_LGL_1d_p2**2
    n_LGL_3d_p2 = n_LGL_1d_p2**3

    n_Gau_1d_p0 = n_LGL_1d_p0
    n_Gau_2d_p0 = n_LGL_2d_p0
    n_Gau_3d_p0 = n_LGL_3d_p0

    n_Gau_1d_p1 = n_LGL_1d_p1
    n_Gau_2d_p1 = n_LGL_2d_p1
    n_Gau_3d_p1 = n_LGL_3d_p1

    n_Gau_1d_p2 = n_LGL_1d_p2
    n_Gau_2d_p2 = n_LGL_2d_p2
    n_Gau_3d_p2 = n_LGL_3d_p2

    n_Gau_shell = nfacesperelem*n_Gau_2d_p2

    ! initialize collocation matrices
    allocate(x_LGL_1d_p0(n_LGL_1d_p0))
    allocate(x_LGL_1d_p1(n_LGL_1d_p1))
    allocate(x_LGL_1d_p2(n_LGL_1d_p2))

    allocate(w_LGL_1d_p0(n_LGL_1d_p0))
    allocate(w_LGL_1d_p1(n_LGL_1d_p1))
    allocate(w_LGL_1d_p2(n_LGL_1d_p2))

    call JacobiP11(n_LGL_1D_p0-1,x_LGL_1d_p0)
    call JacobiP11(n_LGL_1D_p1-1,x_LGL_1d_p1)
    call JacobiP11(n_LGL_1D_p2-1,x_LGL_1d_p2)

    allocate(pmat_p0(n_LGL_1d_p0))
    allocate(pinv_p0(n_LGL_1d_p0))
    allocate(qmat_p0(n_LGL_1d_p0,n_LGL_1d_p0))
    allocate(dmat_p0(n_LGL_1d_p0,n_LGL_1d_p0))
    call Amat(x_LGL_1d_p0,n_LGL_1d_p0,pmat_p0,pinv_p0,qmat_p0,dmat_p0) 

    allocate(pmat_p1(n_LGL_1d_p1))
    allocate(pinv_p1(n_LGL_1d_p1))
    allocate(qmat_p1(n_LGL_1d_p1,n_LGL_1d_p1))
    allocate(dmat_p1(n_LGL_1d_p1,n_LGL_1d_p1))
    call Amat(x_LGL_1d_p1,n_LGL_1d_p1,pmat_p1,pinv_p1,qmat_p1,dmat_p1) 

    allocate(pmat_p2(n_LGL_1d_p2))
    allocate(pinv_p2(n_LGL_1d_p2))
    allocate(qmat_p2(n_LGL_1d_p2,n_LGL_1d_p2))
    allocate(dmat_p2(n_LGL_1d_p2,n_LGL_1d_p2))
    call Amat(x_LGL_1d_p2,n_LGL_1d_p2,pmat_p2,pinv_p2,qmat_p2,dmat_p2) 

    ! populate grad matrix
    call gradmatrix(n_LGL_1d_p0, n_LGL_2d_p0, n_LGL_3d_p0, nnzgrad_p0,   &
                    gradmat_p0, iagrad_p0, jagrad_p0, dagrad_p0, qagrad_p0, &
                    pmat_p0, qmat_p0, dmat_p0, pvol_p0, p_surf_p0)

    call gradmatrix(n_LGL_1d_p1, n_LGL_2d_p1, n_LGL_3d_p1, nnzgrad_p1,   &
                    gradmat_p1, iagrad_p1, jagrad_p1, dagrad_p1, qagrad_p1, &
                    pmat_p1, qmat_p1, dmat_p1, pvol_p1, p_surf_p1)

    call gradmatrix(n_LGL_1d_p2, n_LGL_2d_p2, n_LGL_3d_p2, nnzgrad_p2,   &
                    gradmat_p2, iagrad_p2, jagrad_p2, dagrad_p2, qagrad_p2, &
                    pmat_p2, qmat_p2, dmat_p2, pvol_p2, p_surf_p2)

    allocate(x_LGL_pts_1D(nodesperedge))
    call JacobiP11(polyorder,x_LGL_pts_1D)
    
    allocate(pmat(nodesperedge))
    allocate(pinv(nodesperedge))
    allocate(qmat(nodesperedge,nodesperedge))
    allocate(dmat(nodesperedge,nodesperedge))
    call Amat(x_LGL_pts_1D,nodesperedge,pmat,pinv,qmat,dmat) 

    ! populate grad matrix
    call gradmatrix(nodesperedge, nodesperface, nodesperelem, nnzgrad,   &
                    gradmat, iagrad, jagrad, dagrad, qagrad, &
                    pmat, qmat, dmat, pvol, p_surf)

    !  Interpolation matrices rotating LGL <-> Gau points of various orders

    allocate(x_Gau_1d_p0(n_Gau_1d_p0))
    allocate(w_Gau_1d_p0(n_Gau_1d_p0))

    allocate(x_Gau_1d_p1(n_Gau_1d_p1))
    allocate(w_Gau_1d_p1(n_Gau_1d_p1))

    allocate(x_Gau_1d_p2(n_Gau_1d_p2))
    allocate(w_Gau_1d_p2(n_Gau_1d_p2))

    call Gauss_Legendre_points(n_Gau_1d_p0,x_Gau_1d_p0,w_Gau_1d_p0)
    call Gauss_Legendre_points(n_Gau_1d_p1,x_Gau_1d_p1,w_Gau_1d_p1)
    call Gauss_Legendre_points(n_Gau_1d_p2,x_Gau_1d_p2,w_Gau_1d_p2)

    w_LGL_1d_p0(:) = pmat_p0(:) ; w_LGL_1d_p1(:) = pmat_p1(:) ; w_LGL_1d_p2(:) = pmat_p2(:)

    diff_poly = npoly_max - npoly

    allocate(Prolong_LGL_2_Gau_1d(npoly_max+2,npoly_max+2,1+diff_poly,2))
    allocate(Restrct_Gau_2_LGL_1d(npoly_max+2,npoly_max+2,1+diff_poly,2))

    do i = 0, diff_poly

      !  First:  n_LGL == n_Gau  i.e., mortar is same order as face

      n_LGL_1d=n_LGL_1d_p0 + i

      if(allocated(x_LGL_1d)) deallocate(x_LGL_1d) ; allocate(x_LGL_1d(n_LGL_1d))
      if(allocated(w_LGL_1d)) deallocate(w_LGL_1d) ; allocate(w_LGL_1d(n_LGL_1d))
      if(i==0)w_LGL_1d(:) = pmat_p0(:) ;
      if(i==1)w_LGL_1d(:) = pmat_p1(:) ;
      if(i==2)w_LGL_1d(:) = pmat_p2(:) ;
      call JacobiP11(n_LGL_1d-1,x_LGL_1d)

      n_Gau_1d=n_LGL_1d

      if(allocated(x_Gau_1d)) deallocate(x_Gau_1d) ; allocate(x_Gau_1d(n_Gau_1d))
      if(allocated(w_Gau_1d)) deallocate(w_Gau_1d) ; allocate(w_Gau_1d(n_Gau_1d))
      call Gauss_Legendre_points(n_Gau_1d,x_Gau_1d,w_Gau_1d)

      call Rotate_xione_2_xitwo_and_back(n_LGL_1d,n_Gau_1d, &
                                         x_LGL_1d,x_Gau_1d, &
                                         w_LGL_1d,w_Gau_1d, &
                                         Prolong_LGL_2_Gau_1d(1:n_Gau_1d,1:n_LGL_1d,1+i,1), &
                                         Restrct_Gau_2_LGL_1d(1:n_LGL_1d,1:n_Gau_1d,1+i,1) )

      n_Gau_1d=n_LGL_1d + 1

      if(allocated(x_Gau_1d)) deallocate(x_Gau_1d) ; allocate(x_Gau_1d(n_Gau_1d))
      if(allocated(w_Gau_1d)) deallocate(w_Gau_1d) ; allocate(w_Gau_1d(n_Gau_1d))
      call Gauss_Legendre_points(n_Gau_1d,x_Gau_1d,w_Gau_1d)

      call Rotate_xione_2_xitwo_and_back(n_LGL_1d,n_Gau_1d, &
                                         x_LGL_1d,x_Gau_1d, &
                                         w_LGL_1d,w_Gau_1d, &
                                         Prolong_LGL_2_Gau_1d(1:n_Gau_1d,1:n_LGL_1d,1+i,2), &
                                         Restrct_Gau_2_LGL_1d(1:n_LGL_1d,1:n_Gau_1d,1+i,2) )
    enddo

    allocate(Ext_LGL_p0_2_Gau_p1_1d(n_LGL_1d_p0,n_Gau_1d_p1)) ; Ext_LGL_p0_2_Gau_p1_1d = 0.0_wp ;
    allocate(Int_Gau_p1_2_LGL_p0_1d(n_Gau_1d_p1,n_LGL_1d_p0)) ; Int_Gau_p1_2_LGL_p0_1d = 0.0_wp ;

    allocate(Rot_LGL_p0_2_Gau_p0_1d(n_LGL_1d_p0,n_Gau_1d_p0)) ; Rot_LGL_p0_2_Gau_p0_1d = 0.0_wp ;
    allocate(Rot_Gau_p0_2_LGL_p0_1d(n_Gau_1d_p0,n_LGL_1d_p0)) ; Rot_Gau_p0_2_LGL_p0_1d = 0.0_wp ;

    allocate(Rot_LGL_p1_2_Gau_p1_1d(n_LGL_1d_p1,n_Gau_1d_p1)) ; Rot_LGL_p1_2_Gau_p1_1d = 0.0_wp ;
    allocate(Rot_Gau_p1_2_LGL_p1_1d(n_Gau_1d_p1,n_LGL_1d_p1)) ; Rot_Gau_p1_2_LGL_p1_1d = 0.0_wp ;

    call Rotate_xione_2_xitwo_and_back(n_LGL_1d_p0,n_Gau_1d_p0, &
                                       x_LGL_1d_p0,x_Gau_1d_p0, &
                                       w_LGL_1d_p0,w_Gau_1d_p0, &
                                       Rot_LGL_p0_2_Gau_p0_1d,  &
                                       Rot_Gau_p0_2_LGL_p0_1d )

    call Rotate_xione_2_xitwo_and_back(n_LGL_1d_p1,n_Gau_1d_p1, &
                                       x_LGL_1d_p1,x_Gau_1d_p1, &
                                       w_LGL_1d_p1,w_Gau_1d_p1, &
                                       Rot_LGL_p1_2_Gau_p1_1d,  &
                                       Rot_Gau_p1_2_LGL_p1_1d )

    call Rotate_xione_2_xitwo_and_back(n_LGL_1d_p0,n_Gau_1d_p1, &
                                       x_LGL_1d_p0,x_Gau_1d_p1, &
                                       w_LGL_1d_p0,w_Gau_1d_p1, &
                                       Ext_LGL_p0_2_Gau_p1_1d,  &
                                       Int_Gau_p1_2_LGL_p0_1d )

!   Allocate and build the 1D interpolation matrices and prolong (coarse_2_fine)
!   and restrict (fine_2_coarse) data at a non-conforming H-refinement interface 

    allocate(LGL_Coarse_2_LGL_Fine_1d(n_LGL_1d_p0,n_LGL_1d_p0)) ; LGL_Coarse_2_LGL_Fine_1d = 0.0_wp ;
    allocate(LGL_Fine_2_LGL_Coarse_1d(n_LGL_1d_p0,n_LGL_1d_p0)) ; LGL_Fine_2_LGL_Coarse_1d = 0.0_wp ;
                                       
    call Rotate_xione_2_xitwo_and_back(n_LGL_1d_p0,n_LGL_1d_p0,                 &
                                       x_LGL_1d_p0,0.5_wp*(1.0_wp+x_LGL_1d_p0), &
                                       w_LGL_1d_p0,0.5_wp*w_LGL_1d_p0,          &
                                       LGL_Coarse_2_LGL_Fine_1d,                &
                                       LGL_Fine_2_LGL_Coarse_1d )

!   Filter Matrix
    allocate(Filter(nodesperedge,nodesperedge))

    call Filter_GLL_2_GLL(nodesperedge,nodesperedge,Filter)
    call FilterMatrix()

!   Begin Conversion to staggered approach    

    N_Soln_Pts = npoly + 1
    N_Flux_Pts = npoly + 1 + npoly_DeltaF

    allocate(X_Soln_Pts(N_Soln_Pts))
    allocate(X_Flux_Pts(N_Flux_Pts))

    call JacobiP11(N_Soln_Pts-1,X_Soln_Pts)
    call JacobiP11(N_Flux_Pts-1,X_Flux_Pts)

    allocate(pmat_Soln(N_Soln_Pts))
    allocate(pinv_Soln(N_Soln_Pts))
    allocate(qmat_Soln(N_Soln_Pts,N_Soln_Pts))
    allocate(dmat_Soln(N_Soln_Pts,N_Soln_Pts))
    call Amat(X_Soln_Pts,N_Soln_Pts,pmat_Soln,pinv_Soln,qmat_Soln,dmat_Soln) 

    allocate(pmat_Flux(N_Flux_Pts))
    allocate(pinv_Flux(N_Flux_Pts))
    allocate(qmat_Flux(N_Flux_Pts,N_Flux_Pts))
    allocate(dmat_Flux(N_Flux_Pts,N_Flux_Pts))
    call Amat(X_Flux_Pts,N_Flux_Pts,pmat_Flux,pinv_Flux,qmat_Flux,dmat_Flux) 

    allocate(Int_F2S(N_Soln_Pts,N_Flux_Pts))
    allocate(Dif_F2S(N_Soln_Pts,N_Flux_Pts))

    allocate(Ext_S2F(N_Flux_Pts,N_Soln_Pts))

    call ComputeSolutiontoFluxExtrapolationMatrix(N_Soln_Pts,N_Flux_Pts, X_Soln_Pts,X_Flux_Pts, Ext_S2F)
    call ComputeFluxtoSolutionInterpolationMatrix(N_Soln_Pts,N_Flux_Pts, X_Soln_Pts,X_Flux_Pts, Int_F2S)

    Dif_F2S = matmul(Int_F2S,dmat_Flux)

!   Extrapolation matrix needed for SSSCE.  Takes Flux at Solution pts to Flux pts 
!
!   Extrapolation matrix P^{-1}Q f(U) = P^{-1} Delta F(U) = P^{-1} Delta Int_Coeff_S2F f(U)
!
    allocate(Ext_SSSCE_S2F(N_Flux_Pts+1,N_Flux_Pts))
    call Get_Ext_SSSCE_S2F(N_Flux_Pts,Ext_SSSCE_S2F)

    call WENO_Coeffs

  end function rmapInit

end module polyinit
