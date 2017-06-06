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

    implicit none
    ! polynomial order
    integer, intent(in) :: polyorder
    integer, intent(in) :: ndimensions
    integer :: rmapInit

    rmapInit = 0

    ! set polynomial order
    npoly = polyorder
    ! set dimension
    ndim = ndimensions
    ! set number of faces per element
    nfacesperelem = 2 * ndim
    ! each element direction has npoly+1 nodes
    nodesperedge = (npoly+1)
    ! number of nodes per element
    nodesperelem = nodesperedge**ndim
    ! number of nodes on each face
    nodesperface = nodesperedge**(ndim-1)

    n_LGL_1d_pL = npoly + 1 
    n_LGL_2d_pL = n_LGL_1d_pL**2
    n_LGL_3d_pL = n_LGL_1d_pL**3

    n_LGL_1d_pH = npoly + 2 
    n_LGL_2d_pH = n_LGL_1d_pH**2
    n_LGL_3d_pH = n_LGL_1d_pH**3

    ! initialize collocation matrices
    allocate(x_LGL_pts_1d_pL(n_LGL_1d_pL))
    allocate(x_LGL_pts_1d_pH(n_LGL_1d_pH))

    call JacobiP11(n_LGL_1D_pL-1,x_LGL_pts_1d_pL)
    call JacobiP11(n_LGL_1D_pH-1,x_LGL_pts_1d_pH)

    allocate(pmat_pL(n_LGL_1d_pL))
    allocate(pinv_pL(n_LGL_1d_pL))
    allocate(qmat_pL(n_LGL_1d_pL,n_LGL_1d_pL))
    allocate(dmat_pL(n_LGL_1d_pL,n_LGL_1d_pL))
    call Amat(x_LGL_pts_1d_pL,n_LGL_1d_pL,pmat_pL,pinv_pL,qmat_pL,dmat_pL) 

    allocate(pmat_pH(n_LGL_1d_pH))
    allocate(pinv_pH(n_LGL_1d_pH))
    allocate(qmat_pH(n_LGL_1d_pH,n_LGL_1d_pH))
    allocate(dmat_pH(n_LGL_1d_pH,n_LGL_1d_pH))
    call Amat(x_LGL_pts_1d_pH,n_LGL_1d_pH,pmat_pH,pinv_pH,qmat_pH,dmat_pH) 

    allocate(x_LGL_pts_1D(nodesperedge))
    call JacobiP11(polyorder,x_LGL_pts_1D)
    
    allocate(pmat(nodesperedge))
    allocate(pinv(nodesperedge))
    allocate(qmat(nodesperedge,nodesperedge))
    allocate(dmat(nodesperedge,nodesperedge))
    call Amat(x_LGL_pts_1D,nodesperedge,pmat,pinv,qmat,dmat) 

    ! populate grad matrix
    call gradmatrix() ! (initcollocation)

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

!    allocate(x_gll(N_Soln_Pts))
!    allocate(w_gll(N_Soln_Pts))
!    write(*,*) 'order', polyorder
!    call Gauss_Lobatto_Legendre_points(polyorder+1,x_gll,w_gll)

    call WENO_Coeffs

  end function rmapInit

end module polyinit
