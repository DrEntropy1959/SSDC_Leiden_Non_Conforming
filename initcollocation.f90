module initcollocation
  use precision_vars
  implicit none

  private

  public gradmatrix
  public jacobiP11
  public Amat

  public Gauss_Legendre_points
  public Gauss_Lobatto_Legendre_points
  public Extrapolate_GL_2_GLL
  public Interpolate_GLL_2_GL
  public Intrpltnmatrix
  public Extrpltnmatrix
  public ExtrpXA2XB_2D
  public ExtrpXA2XB_3D
  public Filter_GLL_2_GLL
  public FilterMatrix
  public Get_Ext_SSSCE_S2F
  public ComputeFluxtoSolutionInterpolationMatrix
  public ComputeSolutionToFluxExtrapolationMatrix
  public compute_gsat_f2s_matrix
  public element_properties

contains

  subroutine Amat(x, N, pmat,pinv,qmat,dmat)                                
    ! Form the differentiation matrix dmat                                  
    ! This matrix is used to compute the derivative in the Gauss-Labotto points   

    implicit none

    integer,                  intent(in)    :: N
    real(wp), dimension(N),   intent(in)    :: x
    real(wp), dimension(N),   intent(inout) :: pmat,pinv
    real(wp), dimension(N,N), intent(inout) :: qmat,dmat

    real(wp), allocatable, dimension(:)     :: wk, f, dex

    real(wp)                                :: err, errmax
    real(wp), parameter                     :: half = 0.50000000000000000000000000_wp
    logical                                 :: diagnostic = .false.

    integer                                 :: i,j,k


    ! Allocate memory
    allocate(wk(N))
    allocate(f(N))
    allocate(dex(N))

    do k = 1,N 
      wk(k) = 1.0_wp
      do i = 1,N 
        if(i/=k) wk(k) = wk(k) * (x(k) - x(i)) 
      enddo 
    enddo 

    do i = 1,N 
      do j = 1,N 
        if(j/=i) dmat(i,j)=wk(i)/ (wk(j) * (x(i) - x(j)) ) 
      enddo 
    enddo 
    do j = 1,N 
      dmat(j,j)=0.0_wp
      do k = 1,N 
        if(k/=j) dmat(j,j) = dmat(j,j) + 1.0_wp/(x(j) - x(k)) 
      enddo 
    enddo 

    pinv(1) = (-dmat(1,1)+dmat(N,N))
    pinv(N) = pinv(1)
    do i = 2,(N+1)/2
      pinv(i) = - (dmat(i,1)*pinv(1)) / dmat(1,i) 
      pinv(N+1-i) = pinv(i)
    enddo

    pmat(:) = 1.0_wp / pinv(:)

    do i = 1,N
      do j = 1,N
        qmat(i,j) = dmat(i,j)*pmat(i)
      enddo
    enddo

    errmax = (qmat(1,1) + half)**2 + (qmat(N,N) - half)**2
    do i = 1,N
      do j = 1,N
        if(((i /= 1) .or. (j /= 1)) .and.  ((i /= N) .or. (j /= N)))then 
          errmax = errmax + (qmat(i,j) + qmat(j,i))**2
        endif
      enddo
    enddo
    errmax = sqrt(errmax/N/N)

    !       differentiate exact data to test accuracy of discrete operator

    do i=1,N
      f(i)    =         x(i)**(N-1)
      dex(i)  = (N-1) * x(i)**(N-2)
    enddo 
    do i = 1,N 
      wk(i)=0.0_wp
      do j = 1,N 
        wk(i) = wk(i) + dmat(i,j) * f(j) 
      enddo 
    enddo 

    do i = 1, N
      err=abs(wk(i)-dex(i)) 
      if(err.gt.errmax)errmax=err 
    enddo 

    if(diagnostic) then
      do i = 1,N
         write(*,287)i,pmat(i)
      enddo
      do i = 1,N
         do j = i,N
           write(*,288)i,j,qmat(i,j)
         enddo
      enddo
      stop
    endif
 287  format('      d',i2,' = ',f22.17)
 288  format('      q',i2,1x,i2,' = ',f22.17)
    if(errmax > 1.0e-13_wp) then
      write(*,*)'x',x
      write(*,*)'dmat'
      do i = 1,N
        write(*,*)(dmat(i,j),j=1,N)
      enddo
      write(*,*)'roundoff alert',errmax 
    endif

    deallocate(f  )
    deallocate(dex)
    deallocate(wk )

    !   150 format( 25(f5.2,1x)) 
    !   160 format( 25(f6.3,1x))

    return 
  end subroutine Amat

  subroutine ComputeFluxtoSolutionInterpolationMatrix(ns,nf, solpnts,flxpnts, IntrpFlx2Sol)
    !   Build rotation matrices used to take flux data ==> solution data
    !   ns:   nbrSolPnts
    !   nf:   nbrFlxPnts

    implicit none 

    integer,                  intent( in) :: ns,nf
    real(wp), dimension(ns),  intent( in) :: solpnts
    real(wp), dimension(nf),  intent( in) :: flxpnts
        real(wp), dimension(ns,nf),  intent(inout) :: IntrpFlx2Sol

    integer      :: iFlx, iSol, iTerm

    do iSol = 1,ns
      do iFlx = 1,nf
            IntrpFlx2Sol(iSol,iFlx) = 1.0_wp
        do iTerm = 1,nf
          if(iTerm /= iFlx) then
                IntrpFlx2Sol(iSol,iFlx) = IntrpFlx2Sol(iSol,iFlx)*&
                (SOLPNTS(iSol)-FLXPNTS(iTerm))                &
                /(FLXPNTS(iFlx)-FLXPNTS(iTerm))
          endif
        enddo

      enddo

    enddo

    return
  end subroutine ComputeFluxtoSolutionInterpolationMatrix

  subroutine ComputeSolutionToFluxExtrapolationMatrix(ns,nf, solpnts,flxpnts, ExtrpSol2Flx)
    ! Build rotation matrices used to take solution data ==> flux data
    ! ns:   nbrSolPnts
    ! nf:   nbrFlxPnts

    implicit none 

    integer,                  intent( in) :: ns,nf
    real(wp), dimension(ns),  intent( in) :: solpnts
    real(wp), dimension(nf),  intent( in) :: flxpnts
        real(wp), dimension(nf,ns),  intent(inout) :: ExtrpSol2Flx

    integer      :: iFlx, iSol, iFactor

    do iFlx = 1,nf
      do iSol = 1,ns
            ExtrpSol2Flx(iFlx,iSol) = 1.0_wp
        do iFactor = 1,ns
          if(iFactor /= iSol) &
                ExtrpSol2Flx(iFlx,iSol) = ExtrpSol2Flx(iFlx,iSol)*&
            (FLXPNTS(iFlx)-SOLPNTS(iFactor))                &
            /(SOLPNTS(iSol)-SOLPNTS(iFactor)) 
        enddo
      enddo
    enddo

    return
  end subroutine ComputeSolutionToFluxExtrapolationMatrix

    subroutine element_properties(ielem, n_pts_1d, n_pts_2d, n_pts_3d, &
                      pinv, qmat, dmat, iagrad, jagrad, dagrad,        &
                      pmat, nnzgrad, pvol, p_surf)

    use collocationvariables, only : n_LGL_1d_pL, n_LGL_2d_pL, n_LGL_3d_pL, &
                                   & pinv_pL, pmat_pL, qmat_pL, dmat_pL, &
                                   & iagrad_pL, jagrad_pL, nnzgrad_pL, &
                                   & dagrad_pL, pvol_pL, p_surf_pL, &
                                   & n_LGL_1d_pH, n_LGL_2d_pH, n_LGL_3d_pH, &
                                   & pinv_pH, pmat_pH, qmat_pH, dmat_pH, &
                                   & iagrad_pH, jagrad_pH, nnzgrad_pH, &
                                   & dagrad_pH, pvol_pH, p_surf_pH, &
                                   & elem_props  !, &
                                  ! & qagrad_pL, qagrad_pH

    use referencevariables
    implicit none

    integer,                                   intent(in   ) :: ielem
    integer,                                   intent(  out) :: n_pts_1d, n_pts_2d, n_pts_3d
    real(wp), dimension(:),      allocatable,  intent(  out) :: pinv
    real(wp), dimension(:,:),    allocatable,  intent(  out) :: dmat, qmat

!   real(wp), dimension(:,:,:),  allocatable,  intent(inout) :: gradmat
    integer,  dimension(:),      allocatable,  intent(inout) :: iagrad
    integer,  dimension(:,:),    allocatable,  intent(inout) :: jagrad
!   real(wp), dimension(:,:),    allocatable,  intent(inout) :: dagrad, qagrad
    real(wp), dimension(:,:),    allocatable,  intent(inout) :: dagrad

    integer,                                  optional, intent(  out) :: nnzgrad
    real(wp), dimension(:),      allocatable, optional, intent(  out) :: pmat
    real(wp), dimension(:),      allocatable, optional, intent(inout) :: pvol
    real(wp), dimension(:),      allocatable, optional, intent(inout) :: p_surf
 
    if((elem_props(1,ielem) == 1) .and. (elem_props(2,ielem) == npoly+1)) then
      n_pts_1d = n_LGL_1d_pL
      n_pts_2d = n_LGL_2d_pL
      n_pts_3d = n_LGL_3d_pL
      nnzgrad  = nnzgrad_pL
    else
      n_pts_1d = n_LGL_1d_pH
      n_pts_2d = n_LGL_2d_pH
      n_pts_3d = n_LGL_3d_pH
      nnzgrad  = nnzgrad_pH
    endif

    if(allocated(pmat)) deallocate(pmat) ; allocate(pmat(n_pts_1d)) ;
    if(allocated(pinv)) deallocate(pinv) ; allocate(pinv(n_pts_1d)) ;
    if(allocated(qmat)) deallocate(qmat) ; allocate(qmat(n_pts_1d,n_pts_1d)) ;
    if(allocated(dmat)) deallocate(dmat) ; allocate(dmat(n_pts_1d,n_pts_1d)) ;


    if(allocated(iagrad)) deallocate(iagrad) ; allocate(iagrad(n_pts_3d+1)) ;
    if(allocated(jagrad)) deallocate(jagrad) ; allocate(jagrad(3,nnzgrad)) ;
    if(allocated(dagrad)) deallocate(dagrad) ; allocate(dagrad(3,nnzgrad)) ;

    if(allocated(pvol  )) deallocate(pvol  ) ; allocate(pvol(n_pts_3d)) ;
    if(allocated(p_surf)) deallocate(p_surf) ; allocate(p_surf(n_pts_2d)) ;

    if((elem_props(1,ielem) == 1) .and. (elem_props(2,ielem) == npoly+1)) then
      pmat(:)     = pmat_pL(:)
      pinv(:)     = pinv_pL(:)
      qmat(:,:)   = qmat_pL(:,:)
      dmat(:,:)   = dmat_pL(:,:)

      iagrad(:)   = iagrad_pL(:)
      jagrad(:,:) = jagrad_pL(:,:)
      dagrad(:,:) = dagrad_pL(:,:)
!     qagrad(:,:) = qagrad_pL(:,:)

      pvol(:)     = pvol_pL(:)
      p_surf(:)   = p_surf_pL(:)
    else
      pmat(:)     = pmat_pH(:)
      pinv(:)     = pinv_pH(:)
      qmat(:,:)   = qmat_pH(:,:)
      dmat(:,:)   = dmat_pH(:,:)

      iagrad(:)   = iagrad_pH(:)
      jagrad(:,:) = jagrad_pH(:,:)
      dagrad(:,:) = dagrad_pH(:,:)
!     qagrad(:,:) = qagrad_pH(:,:)

      pvol(:)     = pvol_pH(:)
      p_surf(:)   = p_surf_pH(:)
    endif

    end subroutine 

  ! ===================================================================================
  ! This subroutine calculates the multi-dimensional gradient operators 
  ! for an already set dimension and polynomial approximation.
  ! ===================================================================================

  subroutine gradmatrix(n_pts_1d, n_pts_2d, n_pts_3d, nnzgrad,   &
                        gradmat, iagrad, jagrad, dagrad, qagrad, &
                        pmat, qmat, dmat, pvol, p_surf)

!   use collocationvariables
    use referencevariables
    implicit none

    integer,                                   intent(in)    :: n_pts_1d, n_pts_2d, n_pts_3d
    integer,                                   intent(  out) :: nnzgrad
    real(wp), dimension(n_pts_1d),             intent(in)    :: pmat
    real(wp), dimension(n_pts_1d,n_pts_1d),    intent(in)    :: dmat, qmat

    real(wp), dimension(:,:,:),  allocatable,  intent(inout) :: gradmat
    integer,  dimension(:),      allocatable,  intent(inout) :: iagrad
    integer,  dimension(:,:),    allocatable,  intent(inout) :: jagrad
    real(wp), dimension(:,:),    allocatable,  intent(inout) :: dagrad, qagrad
    real(wp), dimension(:),      allocatable,  intent(inout) :: pvol
    real(wp), dimension(:),      allocatable,  intent(inout) :: p_surf

    integer :: i,j,k,jj,kk,inode, idir
    integer :: il(2,3), ix(3), ix_surf(2)
    integer, allocatable, dimension(:) :: stride
    integer :: icount

    ! high and low indices in each dimension
    il = 1
    do idir = 1,ndim
      il(2,idir) = n_pts_1d
    end do

    ! initialize gradient matrix for each direction
    allocate(gradmat(n_pts_3d,n_pts_3d,ndim))
    gradmat = 0.0_wp

    ! allocate csr storage
    ! ====================
    ! number of nonzeroes is the same in each direction
    nnzgrad = n_pts_3d*n_pts_1d
    
    ! row pointers
    allocate(iagrad(n_pts_3d+1))
    iagrad = 0
    
    ! column pointers
    allocate(jagrad(1:ndim,nnzgrad))
    jagrad = 0

    ! matrix values
    allocate(dagrad(1:ndim,nnzgrad))
    dagrad = 0.0_wp

    allocate(qagrad(1:ndim,nnzgrad))
    qagrad = 0.0_wp

    ! stride between columns in each direction
    allocate(stride(1:ndim))
    do idir = 1,ndim
      stride(idir) = n_pts_1d**(idir-1)
    end do

    ! set CSR pointer index counter
    icount = 1
    ! reset row index coutner
    inode = 0
    ! loop over third direction (will be 1 to 1 unless 3D)
    do k = il(1,3), il(2,3)
      ! loop over second dimension (1 to 1 unless 2D)
      do j = il(1,2), il(2,2)
        ! loop over first dimension
        do i = il(1,1), il(2,1)
          ! index vector
          ix = (/ i,j,k /)
          ! advance row index by 1
          inode = inode+1
          ! set CSR counter for row
          iagrad(inode) = icount
          ! each gradient requires n_pts_1d nonzero columns
          ! loop over n_pts_1d
          do kk = 1,n_pts_1d
            ! loop over directions
            do idir = 1,ndim
              ! column/node corresponding to coefficient in idir-direction
              jj = stride(idir)*(kk-ix(idir)) + inode
              jagrad(idir,icount) = jj
              ! set coefficient for dmat
              dagrad(idir,icount) = dmat(ix(idir),kk)
              ! set coefficient for qmat
              qagrad(idir,icount) = qmat(ix(idir),kk)
            end do
            ! advance CSR counter by 1
            icount = icount + 1
          end do
        end do
      end do
    end do
    ! set last pointer index for CSR
    iagrad(n_pts_3d+1) = icount

    ! gradmat is the same as dmat but is a full matrix
    ! loop over directions
    do idir = 1, ndim
      ! reset node counter
      inode = 0
      ! loop over third direction
      do k = il(1,3), il(2,3)
        ! loop over second direction
        do j = il(1,2), il(2,2)
          ! loop over first direction
          do i = il(1,1), il(2,1)
            ! index vector
            ix = (/ i,j,k /)
            ! advance node counter by 1
            inode = inode+1
            ! loop over coefficients in gradient operator
            do kk = 1,n_pts_1d
              ! column/node corresponding to coefficient
              jj = stride(idir)*(kk-ix(idir)) + inode
              ! set coefficient in larger gradient matrix
              gradmat(inode,jj,idir) = dmat(ix(idir),kk)
            end do
          end do
        end do
      end do
    end do

    ! pvol is used to integrate over the volume. It represents
    ! the local volume in computational space. It is a diagonal
    ! matrix, and thus set to a vector to save space.
    allocate(pvol(n_pts_3d))
    ! initialize to 1 because the length of every dimension not included is 1.
    pvol = 1.0_wp
    ! loop over directions
    do idir = 1, ndim
      ! reset node counter
      inode = 0
      ! loop over third direction
      do k = il(1,3), il(2,3)
        ! loop over second direction
        do j = il(1,2), il(2,2)
          ! loop over first direction
          do i = il(1,1), il(2,1)
            ! index vector
            ix = (/ i,j,k /)
            ! advance node counter by 1
            inode = inode+1
            ! update the product of the volume
            pvol(inode) = pvol(inode)*pmat(ix(idir))
          end do
        end do
      end do
    end do

    ! p_surf is used to integrate over the surface. It represents
    ! the local surface in computational space. It is a diagonal
    ! matrix, and thus set to a vector to save space.
    allocate(p_surf(n_pts_2d))
    ! initialize to 1 because the length of every dimension not included is 1.
    p_surf = 1.0_wp
    ! loop over directions
    do idir = 1, ndim-1
      ! reset node counter
      inode = 0
      ! loop over third direction
      !do k = il(1,3), il(2,3)
        ! loop over second direction
        do j = il(1,2), il(2,2)
          ! loop over first direction
          do i = il(1,1), il(2,1)
            ! index vector
            ix_surf = (/ i,j /)
            ! advance node counter by 1
            inode = inode+1
            ! update the product of the volume
            p_surf(inode) = p_surf(inode)*pmat(ix_surf(idir))
          end do
        end do
    end do

    deallocate(stride)

  end subroutine gradmatrix

  subroutine FilterMatrix()

    use collocationvariables
    use referencevariables
    implicit none

    integer :: i,j,k,jj,kk,inode, idir
    integer :: il(2,3), ix(3)
    integer, allocatable, dimension(:) :: stride
    integer :: icount

    ! high and low indices in each dimension
    il = 1
    do idir = 1,ndim
      il(2,idir) = nodesperedge
    end do

    ! initialize gradient matrix for each direction
    allocate(mat_Filter(nodesperelem,nodesperelem,ndim))
    mat_Filter = 0.0_wp

    ! allocate CSR storage
    ! ====================
    ! number of nonzeroes is the same in each direction
    nnz_Filter = nodesperelem*nodesperedge
    ! row pointers
    allocate(ia_Filter(nodesperelem+1))
    ia_Filter = 0
    ! column pointers
    allocate(ja_Filter(1:ndim,nnz_Filter))
    ja_Filter = 0
    ! matrix values
    allocate(aa_Filter(1:ndim,nnz_Filter))
    aa_Filter = 0.0_wp

    ! stride between columns in each direction
    allocate(stride(1:ndim))
    do idir = 1,ndim
      stride(idir) = nodesperedge**(idir-1)
    end do

    ! set CSR pointer index counter
    icount = 1
    ! reset row index counter
    inode = 0
    ! loop over third direction (will be 1 to 1 unless 3D)
    do k = il(1,3), il(2,3)
      ! loop over second dimension (1 to 1 unless 2D)
      do j = il(1,2), il(2,2)
        ! loop over first dimension
        do i = il(1,1), il(2,1)
          ! index vector
          ix = (/ i,j,k /)
          ! advance row index by 1
          inode = inode+1
          ! set CSR counter for row
          ia_Filter(inode) = icount
          ! each Filter requires nodesperedge nonzero columns
          ! loop over nodesperedge
          do kk = 1,nodesperedge
            ! loop over directions
            do idir = 1,ndim
              ! column/node corresponding to coefficient in idir-direction
              jj = stride(idir)*(kk-ix(idir)) + inode
              ja_Filter(idir,icount) = jj
              ! set coefficient for Fmat
              aa_Filter(idir,icount) = Filter(ix(idir),kk)
            end do
            ! advance CSR counter by 1
            icount = icount + 1
          end do
        end do
      end do
    end do
    ! set last pointer index for CSR
    ia_Filter(nodesperelem+1) = icount

    ! mat_Filter is the same as dmat but is a full matrix
    ! loop over directions
    do idir = 1, ndim
      ! reset node counter
      inode = 0
      ! loop over third direction
      do k = il(1,3), il(2,3)
        ! loop over second direction
        do j = il(1,2), il(2,2)
          ! loop over first direction
          do i = il(1,1), il(2,1)
            ! index vector
            ix = (/ i,j,k /)
            ! advance node counter by 1
            inode = inode+1
            ! loop over coefficients in Filter operator
            do kk = 1,nodesperedge
              ! column/node corresponding to coefficient
              jj = stride(idir)*(kk-ix(idir)) + inode
              ! set coefficient in larger Filter matrix
              mat_Filter(inode,jj,idir) = Filter(ix(idir),kk)
            end do
          end do
        end do
      end do
    end do

    deallocate(stride)

  end subroutine FilterMatrix


  subroutine Intrpltnmatrix()
    use collocationvariables
    use referencevariables
    implicit none

    integer :: i,j,k,jj,kk,inode, idir
    integer :: il(2,3), ix(3)
    integer, allocatable, dimension(:) :: stride
    integer :: icount

    ! high and low indices in each dimension
    il = 1
    do idir = 1,ndim
      il(2,idir) = nodesperedge
    end do

    ! initialize gradient matrix for each direction
    allocate(mat_Intrpltn(nodesperelem,nodesperelem,ndim))
    mat_Intrpltn = 0.0_wp

    ! allocate csr storage
    ! ====================
    ! number of nonzeroes is the same in each direction
    nnz_Intrpltn = (N_Soln_Pts**2) * N_Soln_Pts*N_Flux_Pts
    ! row pointers
    allocate(ia_Intrpltn(nodesperelem+1))
    ia_Intrpltn = 0
    ! column pointers
    allocate(ja_Intrpltn(1:ndim,nnz_Intrpltn))
    ja_Intrpltn = 0
    ! matrix values
    allocate(aa_Intrpltn(1:ndim,nnz_Intrpltn))
    aa_Intrpltn = 0.0_wp

    ! stride between columns in each direction
    allocate(stride(1:ndim))
    do idir = 1,ndim
      stride(idir) = N_Soln_Pts**(idir-1)
    end do

    ! set CSR pointer index counter
    icount = 1
    ! reset row index coutner
    inode = 0
    ! loop over third direction (will be 1 to 1 unless 3D)
    do k = il(1,3), il(2,3)
      ! loop over second dimension (1 to 1 unless 2D)
      do j = il(1,2), il(2,2)
        ! loop over first dimension
        do i = il(1,1), il(2,1)
          ! index vector
          ix = (/ i,j,k /)
          ! advance row index by 1
          inode = inode+1
          ! set CSR counter for row
          ia_Intrpltn(inode) = icount
          ! each Interpolation requires nodesperedge nonzero columns
          ! loop over nodesperedge
          do kk = 1,nodesperedge
            ! loop over directions
            do idir = 1,ndim
              ! column/node corresponding to coefficient in idir-direction
              jj = stride(idir)*(kk-ix(idir)) + inode
              ja_Intrpltn(idir,icount) = jj
              ! set coefficient for dmat
              aa_Intrpltn(idir,icount) = Int_F2S(ix(idir),kk)
            end do
            ! advance CSR counter by 1
            icount = icount + 1
          end do
        end do
      end do
    end do
    ! set last pointer index for CSR
    ia_Intrpltn(nodesperelem+1) = icount

    ! matgrad is the same as dmat but is a full matrix
    ! loop over directions
    do idir = 1, ndim
      ! reset node counter
      inode = 0
      ! loop over third direction
      do k = il(1,3), il(2,3)
        ! loop over second direction
        do j = il(1,2), il(2,2)
          ! loop over first direction
          do i = il(1,1), il(2,1)
            ! index vector
            ix = (/ i,j,k /)
            ! advance node counter by 1
            inode = inode+1
            ! loop over coefficients in gradient operator
            do kk = 1,nodesperedge
              ! column/node corresponding to coefficient
              jj = stride(idir)*(kk-ix(idir)) + inode
              ! set coefficient in larger gradient matrix
              mat_Intrpltn(inode,jj,idir) = Int_F2S(ix(idir),kk)
            end do
          end do
        end do
      end do
    end do

    deallocate(stride)

  end subroutine 

  subroutine Extrpltnmatrix()

    use collocationvariables
    use referencevariables
    implicit none

    integer :: i,j,k,jj,kk,inode, idir
    integer :: il(2,3), ix(3)
    integer, allocatable, dimension(:) :: stride
    integer :: icount

    ! high and low indices in each dimension
    il = 1
    do idir = 1,ndim
      il(2,idir) = nodesperedge
    end do

    ! initialize gradient matrix for each direction
    allocate(mat_Extrpltn(nodesperelem,nodesperelem,ndim))
    mat_Extrpltn = 0.0_wp

    ! allocate csr storage
    ! ====================
    ! number of nonzeroes is the same in each direction
    nnz_Extrpltn = nodesperelem * N_Flux_Pts
    ! row pointers
    allocate(ia_Extrpltn(nodesperelem+1))
    ia_Extrpltn = 0
    ! column pointers
    allocate(ja_Extrpltn(1:ndim,nnz_Extrpltn))
    ja_Extrpltn = 0
    ! matrix values
    allocate(aa_Extrpltn(1:ndim,nnz_Extrpltn))
    aa_Extrpltn = 0.0_wp

    ! stride between columns in each direction
    allocate(stride(1:ndim))

    do idir = 1,ndim
      stride(idir) = N_Soln_Pts**(idir-1)
    end do

    ! set CSR pointer index counter
    icount = 1
    ! reset row index coutner
    inode = 0
    ! loop over third direction (will be 1 to 1 unless 3D)
    do k = il(1,3), il(2,3)
      ! loop over second dimension (1 to 1 unless 2D)
      do j = il(1,2), il(2,2)
        ! loop over first dimension
        do i = il(1,1), il(2,1)
          ! index vector
          ix = (/ i,j,k /)
          ! advance row index by 1
          inode = inode+1
          ! set CSR counter for row
          ia_Extrpltn(inode) = icount
          ! each _Extrpltn requires N_Flux_Pts nonzero columns
          ! loop over N_Flux_Pts
          do kk = 1,N_Flux_Pts
            ! loop over directions
            do idir = 1,ndim
              ! column/node corresponding to coefficient in idir-direction
              jj = stride(idir)*(kk-ix(idir)) + inode
              ja_Extrpltn(idir,icount) = jj
              ! set coefficient for dmat
              aa_Extrpltn(idir,icount) = Ext_S2F(ix(idir),kk)
            end do
            ! advance CSR counter by 1
            icount = icount + 1
          end do
        end do
      end do
    end do
    ! set last pointer index for CSR
    ia_Extrpltn(nodesperelem+1) = icount

    ! mat_Extrpltn is the same as dmat but is a full matrix
    ! loop over directions
    do idir = 1, ndim
      ! reset node counter
      inode = 0
      ! loop over third direction
      do k = il(1,3), il(2,3)
        ! loop over second direction
        do j = il(1,2), il(2,2)
          ! loop over first direction
          do i = il(1,1), il(2,1)
            ! index vector
            ix = (/ i,j,k /)
            ! advance node counter by 1
            inode = inode+1
            ! loop over coefficients in _Extrpltn operator
            do kk = 1,N_Flux_Pts
              ! column/node corresponding to coefficient
              jj = stride(idir)*(kk-ix(idir)) + inode
              ! set coefficient in larger _Extrpltn matrix
              mat_Extrpltn(inode,jj,idir) = Ext_S2F(ix(idir),kk)
            end do
          end do
        end do
      end do
    end do

    deallocate(stride)

  end subroutine Extrpltnmatrix

!=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine ExtrpXA2XB_2D(NPtsA,NPtsB,XA,XB,fA,fB,Extrp)

  ! Extrapolate from Tensor product XA points to Tensor product XB points

  implicit none

  integer,                          intent(in   )  :: NPtsA, NPtsB
  real(wp), dimension(NPtsA),       intent(in   )  :: XA
  real(wp), dimension(NPtsB),       intent(in   )  :: XB
  real(wp), dimension(NPtsB,NptsA), intent(in   )  :: Extrp
  real(wp), dimension(:),           intent(in   )  :: fA
  real(wp), dimension(:),           intent(inout)  :: fB

  real(wp), allocatable, dimension(:,:)   :: F1

  integer                                 :: i,j,m,n
  integer                                 :: StrideY

    allocate(F1(NptsB,NptsA))

    ! Extrapolate in the xi direction;
    StrideY = NPtsA
    F1(:,:) = 0.0_wp
    do j = 1,NPtsA
      do i = 1,NPtsB
        do m = 1,NPtsA
          n = + (j-1)*StrideY + m
          F1(i,j) = F1(i,j) + Extrp(i,m)*FA(n)
        enddo
      enddo
    enddo

    ! Extrapolate in the eta direction;
    StrideY = NPtsB
    FB(:) = 0.0_wp
    do j = 1,NPtsB
      do i = 1,NPtsB
        do m = 1,NPtsA
          n = + (j-1)*StrideY + i
          FB(n) = FB(n) + Extrp(j,m)*F1(i,m)
        enddo
      enddo
    enddo

    deallocate(F1)

  end subroutine ExtrpXA2XB_2D

!=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine ExtrpXA2XB_3D(NPtsA,NPtsB,XA,XB,fA,fB,Extrp)

  ! Extrapolate from Tensor product XA points to Tensor product XB points

  implicit none

  integer,                          intent(in   )  :: NPtsA, NPtsB
  real(wp), dimension(NPtsA),       intent(in   )  :: XA
  real(wp), dimension(NPtsB),       intent(in   )  :: XB
  real(wp), dimension(NPtsB,NptsA), intent(in   )  :: Extrp
  real(wp), dimension(:),           intent(in   )  :: fA
  real(wp), dimension(:),           intent(inout)  :: fB


  real(wp), allocatable, dimension(:,:,:)  :: F1
  real(wp), allocatable, dimension(:,:,:)  :: F2

  integer                                 :: i,j,k,m,n
  integer                                 :: StrideY, StrideZ


    allocate(F1(NptsB,NptsA,NptsA))
    allocate(F2(NptsB,NptsB,NptsA))

    StrideY = NPtsA
    strideZ = NPtsA * NPtsA
    F1(:,:,:) = 0.0_wp
    do k = 1,NPtsA
      do j = 1,NPtsA
        do i = 1,NPtsB
          do m = 1,NPtsA
            n = (k-1)*strideZ + (j-1)*StrideY + m
            F1(i,j,k) = F1(i,j,k) + Extrp(i,m)*FA(n)
          enddo
        enddo
      enddo
    enddo

    F2(:,:,:) = 0.0_wp
    do k = 1,NPtsA
      do j = 1,NPtsB
        do i = 1,NPtsB
          do m = 1,NPtsA
            F2(i,j,k) = F2(i,j,k) + Extrp(j,m)*F1(i,m,k)
          enddo
        enddo
      enddo
    enddo

    StrideY = NPtsB
    strideZ = NPtsB * NPtsB
    FB(:) = 0.0_wp
    do k = 1,NPtsB
      do j = 1,NPtsB
        do i = 1,NPtsB
          do m = 1,NPtsA
            n = (k-1)*strideZ + (j-1)*StrideY + i
            FB(n) = FB(n) + Extrp(k,m)*F2(i,j,m)
          enddo
        enddo
      enddo
    enddo

    deallocate(F1)
    deallocate(F2)

  end subroutine ExtrpXA2XB_3D

!=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!   What follows is a HUGE amount of coefficient data
!   All numbers were generated using Mathematica in N[*,30]  format.
!=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine JacobiP11(p,x)
    !
    !     Subroutine specifies the Jacobi(1,1) Polynomials of degree "p"
    !     As reference  Jacobi(0,0) are Legendre
    !

    implicit none 

    integer,                    intent( in) :: p
    real(wp), dimension(p+1),   intent(out) :: x

    select case(p)

    case(01) !       P = 1 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) =  1.000000000000000000000000000000_wp
    case(02) !       P = 2 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) =  0.000000000000000000000000000000_wp
      x( 3) =  1.000000000000000000000000000000_wp
    case(03) !       P = 3 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.447213595499957939281834733746_wp
      x( 3) =  0.447213595499957939281834733746_wp
      x( 4) =  1.000000000000000000000000000000_wp
    case(04) !       P = 4 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.654653670707977143798292456247_wp
      x( 3) =  0.000000000000000000000000000000_wp
      x( 4) =  0.654653670707977143798292456247_wp
      x( 5) =  1.000000000000000000000000000000_wp
    case(05) !       P = 5 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.765055323929464692851002973959_wp
      x( 3) = -0.285231516480645096314150994041_wp
      x( 4) =  0.285231516480645096314150994041_wp
      x( 5) =  0.765055323929464692851002973959_wp
      x( 6) =  1.000000000000000000000000000000_wp
    case(06) !       P = 6 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.830223896278566929872032213967_wp
      x( 3) = -0.468848793470714213803771881909_wp
      x( 4) =  0.000000000000000000000000000000_wp
      x( 5) =  0.468848793470714213803771881909_wp
      x( 6) =  0.830223896278566929872032213967_wp
      x( 7) =  1.000000000000000000000000000000_wp
    case(07) !       P = 7 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.871740148509606615337445761221_wp
      x( 3) = -0.591700181433142302144510731398_wp
      x( 4) = -0.209299217902478868768657260345_wp
      x( 5) =  0.209299217902478868768657260345_wp
      x( 6) =  0.591700181433142302144510731398_wp
      x( 7) =  0.871740148509606615337445761221_wp
      x( 8) =  1.000000000000000000000000000000_wp
    case(08) !       P = 8 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.899757995411460157312345244418_wp
      x( 3) = -0.677186279510737753445885427091_wp
      x( 4) = -0.363117463826178158710752068709_wp
      x( 5) =  0.000000000000000000000000000000_wp
      x( 6) =  0.363117463826178158710752068709_wp
      x( 7) =  0.677186279510737753445885427091_wp
      x( 8) =  0.899757995411460157312345244418_wp
      x( 9) =  1.000000000000000000000000000000_wp
    case(09) !       P = 9 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.919533908166458813828932660822_wp
      x( 3) = -0.738773865105505075003106174860_wp
      x( 4) = -0.477924949810444495661175092731_wp
      x( 5) = -0.165278957666387024626219765958_wp
      x( 6) =  0.165278957666387024626219765958_wp
      x( 7) =  0.477924949810444495661175092731_wp
      x( 8) =  0.738773865105505075003106174860_wp
      x( 9) =  0.919533908166458813828932660822_wp
      x(10) =  1.000000000000000000000000000000_wp
    case(10) !       P = 10 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.934001430408059134332274136099_wp
      x( 3) = -0.784483473663144418622417816108_wp
      x( 4) = -0.565235326996205006470963969478_wp
      x( 5) = -0.295758135586939391431911515559_wp
      x( 6) =  0.000000000000000000000000000000_wp
      x( 7) =  0.295758135586939391431911515559_wp
      x( 8) =  0.565235326996205006470963969478_wp
      x( 9) =  0.784483473663144418622417816108_wp
      x(10) =  0.934001430408059134332274136099_wp
      x(11) =  1.000000000000000000000000000000_wp
    case(11) !       P = 11 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.944899272222882223407580138303_wp
      x( 3) = -0.819279321644006678348641581717_wp
      x( 4) = -0.632876153031860677662404854444_wp
      x( 5) = -0.399530940965348932264349791567_wp
      x( 6) = -0.136552932854927554864061855740_wp
      x( 7) =  0.136552932854927554864061855740_wp
      x( 8) =  0.399530940965348932264349791567_wp
      x( 9) =  0.632876153031860677662404854444_wp
      x(10) =  0.819279321644006678348641581717_wp
      x(11) =  0.944899272222882223407580138303_wp
      x(12) =  1.000000000000000000000000000000_wp

    case(12) ! P = 12
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.953309846642163911896905464755_wp
      x( 3) = -0.846347564651872316865925607099_wp
      x( 4) = -0.686188469081757426072759039566_wp
      x( 5) = -0.482909821091336201746937233637_wp
      x( 6) = -0.249286930106239992568673700374_wp
      x( 7) =  0.000000000000000000000000000000_wp
      x( 8) =  0.249286930106239992568673700374_wp
      x( 9) =  0.482909821091336201746937233637_wp
      x(10) =  0.686188469081757426072759039566_wp
      x(11) =  0.846347564651872316865925607099_wp
      x(12) =  0.953309846642163911896905464755_wp
      x(13) =  1.000000000000000000000000000000_wp
      
    case(13) ! P = 13
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.959935045267260901355100162015_wp
      x( 3) = -0.867801053830347251000220202908_wp
      x( 4) = -0.728868599091326140584672400521_wp
      x( 5) = -0.550639402928647055316622705859_wp
      x( 6) = -0.342724013342712845043903403642_wp
      x( 7) = -0.116331868883703867658776709736_wp
      x( 8) =  0.116331868883703867658776709736_wp
      x( 9) =  0.342724013342712845043903403642_wp
      x(10) =  0.550639402928647055316622705859_wp
      x(11) =  0.728868599091326140584672400521_wp
      x(12) =  0.867801053830347251000220202908_wp
      x(13) =  0.959935045267260901355100162015_wp
      x(14) =  1.000000000000000000000000000000_wp

    case(14) ! P = 14
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.96524592650383857279585139207_wp
      x( 3) = -0.885082044222976298825401631482_wp
      x( 4) = -0.763519689951815200704118475976_wp
      x( 5) = -0.606253205469845711123529938637_wp
      x( 6) = -0.420638054713672480921896938739_wp
      x( 7) = -0.215353955363794238225679446273_wp
      x( 8) =  0.000000000000000000000000000000_wp
      x( 9) =  0.215353955363794238225679446273_wp
      x(10) =  0.420638054713672480921896938739_wp
      x(11) =  0.606253205469845711123529938637_wp
      x(12) =  0.763519689951815200704118475976_wp
      x(13) =  0.885082044222976298825401631482_wp
      x(14) =  0.96524592650383857279585139207_wp
      x(15) =  1.000000000000000000000000000000_wp

    case(15) !  P = 15
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.969568046270217932952242738367_wp
      x( 3) = -0.89920053309347209299462826152_wp
      x( 4) = -0.792008291861815063931088270963_wp
      x( 5) = -0.652388702882493089467883219641_wp
      x( 6) = -0.486059421887137611781890785847_wp
      x( 7) = -0.299830468900763208098353454722_wp
      x( 8) = -0.101326273521949447843033005046_wp
      x( 9) =  0.101326273521949447843033005046_wp
      x(10) =  0.299830468900763208098353454722_wp
      x(11) =  0.486059421887137611781890785847_wp
      x(12) =  0.652388702882493089467883219641_wp
      x(13) =  0.792008291861815063931088270963_wp
      x(14) =  0.89920053309347209299462826152_wp
      x(15) =  0.969568046270217932952242738367_wp
      x(16) =  1.000000000000000000000000000000_wp

    case(16) ! P = 16
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.973132176631418314156979501874_wp
      x( 3) = -0.910879995915573595623802506398_wp
      x( 4) = -0.815696251221770307106750553238_wp
      x( 5) = -0.691028980627684705394919357372_wp
      x( 6) = -0.541385399330101539123733407504_wp
      x( 7) = -0.372174433565477041907234680735_wp
      x( 8) = -0.189511973518317388304263014753_wp
      x( 9) =  0.000000000000000000000000000000_wp
      x(10) =  0.189511973518317388304263014753_wp
      x(11) =  0.372174433565477041907234680735_wp
      x(12) =  0.541385399330101539123733407504_wp
      x(13) =  0.691028980627684705394919357372_wp
      x(14) =  0.815696251221770307106750553238_wp
      x(15) =  0.910879995915573595623802506398_wp
      x(16) =  0.973132176631418314156979501874_wp
      x(17) =  1.000000000000000000000000000000_wp

    case(17) ! P = 17
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.976105557412198542864518924342_wp
      x( 3) = -0.920649185347533873837854625431_wp
      x( 4) = -0.835593535218090213713646362328_wp
      x( 5) = -0.723679329283242681306210365302_wp
      x( 6) = -0.588504834318661761173535893194_wp
      x( 7) = -0.434415036912123975342287136741_wp
      x( 8) = -0.266362652878280984167665332026_wp
      x( 9) = -0.0897490934846521110226450100886_wp
      x(10) = +0.0897490934846521110226450100886_wp
      x(11) = +0.266362652878280984167665332026_wp
      x(12) = +0.434415036912123975342287136741_wp
      x(13) = +0.588504834318661761173535893194_wp
      x(14) = +0.723679329283242681306210365302_wp
      x(15) = +0.835593535218090213713646362328_wp
      x(16) = +0.920649185347533873837854625431_wp
      x(17) = +0.976105557412198542864518924342_wp
      x(18) = +1.000000000000000000000000000000_wp
    
    case(18) ! P = 18
      
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.97861176622208009515263406311_wp
      x( 3) = -0.928901528152586243717940258797_wp
      x( 4) = -0.852460577796646093085955970041_wp
      x( 5) = -0.751494202552613014163637489634_wp
      x( 6) = -0.628908137265220497766832306229_wp
      x( 7) = -0.488229285680713502777909637625_wp
      x( 8) = -0.333504847824498610298500103845_wp
      x( 9) = -0.169186023409281571375154153445_wp
      x(10) = 0.000000000000000000000000000000_wp
      x(11) = 0.169186023409281571375154153445_wp
      x(12) = 0.333504847824498610298500103845_wp
      x(13) = 0.488229285680713502777909637625_wp
      x(14) = 0.628908137265220497766832306229_wp
      x(15) = 0.751494202552613014163637489634_wp
      x(16) = 0.852460577796646093085955970041_wp
      x(17) = 0.928901528152586243717940258797_wp
      x(18) = 0.97861176622208009515263406311_wp
      x(19) = 1.000000000000000000000000000000_wp

    case default
      write(*,*) 'Jacobi11'
      write(*,*) 'LGL points only defined for p <= 18'
      write(*,*) 'Stopping'
      stop
    end select

    return

  end subroutine JacobiP11

  subroutine Gauss_Legendre_points(p,x,w)
    ! Subroutine provides coefficients for 
    ! Gauss_Legendre Polynomials  where "p"
    ! is the number of GL points
    !
    implicit none 

    integer,                  intent( in) :: p
    real(wp), dimension(p),   intent(out) :: x,w

    select case(p)

    !  Mathematica generator used to generate coefficients

    !    (*  Gauss Points  *)
    !    For[i=1,i<=10,i++,ans=N[Solve[LegendreP[i,x]==0,x],30]  ;
    !    Print[i];
    !    Print[Re[x/.ans]];
    !    Print[Re[N[2/((1-x^2)*(D[LegendreP[i,x],x])^2)/.ans,30]]] ] ;

    case(01) !     P =  1 ;
      x( 1) = +0.000000000000000000000000000000_wp

      w( 1) = +2.00000000000000000000000000000_wp
    case(02) !     P =  2 ;
      x( 1) = -0.577350269189625764509148780502_wp
      x( 2) = +0.577350269189625764509148780502_wp

      w( 1) = +1.000000000000000000000000000000_wp
      w( 2) = +1.000000000000000000000000000000_wp
    case(03) !     P =  3 ;
      x( 1) = -0.774596669241483377035853079956_wp
      x( 2) = +0.000000000000000000000000000000_wp
      x( 3) = +0.774596669241483377035853079956_wp

      w( 1) = +0.555555555555555555555555555556_wp
      w( 2) = +0.888888888888888888888888888889_wp
      w( 3) = +0.555555555555555555555555555556_wp
    case(04) !     P =  4 ;
      x( 1) = -0.861136311594052575223946488893_wp
      x( 2) = -0.339981043584856264802665759103_wp
      x( 3) = +0.339981043584856264802665759103_wp
      x( 4) = +0.861136311594052575223946488893_wp

      w( 1) = +0.34785484513745385737306394922_wp
      w( 2) = +0.65214515486254614262693605078_wp
      w( 3) = +0.65214515486254614262693605078_wp
      w( 4) = +0.34785484513745385737306394922_wp
    case(05) !     P =  5 ;
      x( 1) = -0.906179845938663992797626878299_wp
      x( 2) = -0.538469310105683091036314420700_wp
      x( 3) = +0.000000000000000000000000000000_wp
      x( 4) = +0.538469310105683091036314420700_wp
      x( 5) = +0.906179845938663992797626878299_wp

      w( 1) = +0.2369268850561890875142640407_wp
      w( 2) = +0.4786286704993664680412915148_wp
      w( 3) = +0.5688888888888888888888888889_wp
      w( 4) = +0.4786286704993664680412915148_wp
      w( 5) = +0.2369268850561890875142640407_wp
    case(06) !     P =  6 ;
      x( 1) = -0.932469514203152027812301554494_wp
      x( 2) = -0.661209386466264513661399595020_wp
      x( 3) = -0.238619186083196908630501721681_wp
      x( 4) = +0.238619186083196908630501721681_wp
      x( 5) = +0.661209386466264513661399595020_wp
      x( 6) = +0.932469514203152027812301554494_wp

      w( 1) = +0.171324492379170345040296142_wp
      w( 2) = +0.360761573048138607569833514_wp
      w( 3) = +0.467913934572691047389870344_wp
      w( 4) = +0.467913934572691047389870344_wp
      w( 5) = +0.360761573048138607569833514_wp
      w( 6) = +0.171324492379170345040296142_wp
    case(07) !     P =  7 ;
      x( 1) = -0.949107912342758524526189684048_wp
      x( 2) = -0.741531185599394439863864773281_wp
      x( 5) = -0.405845151377397166906606412077_wp
      x( 4) = +0.000000000000000000000000000000_wp
      x( 3) = +0.405845151377397166906606412077_wp
      x( 6) = +0.741531185599394439863864773281_wp
      x( 7) = +0.949107912342758524526189684048_wp

      w( 1) = +0.129484966168869693270611433_wp
      w( 2) = +0.279705391489276667901467771_wp
      w( 3) = +0.381830050505118944950369775_wp
      w( 4) = +0.417959183673469387755102040_wp
      w( 5) = +0.381830050505118944950369775_wp
      w( 6) = +0.279705391489276667901467771_wp
      w( 7) = +0.129484966168869693270611433_wp
    case(08) !     P =  8 ;
      x( 1) = -0.960289856497536231683560868569_wp
      x( 2) = -0.796666477413626739591553936476_wp
      x( 3) = -0.525532409916328985817739049189_wp
      x( 4) = -0.183434642495649804939476142360_wp
      x( 5) = +0.183434642495649804939476142360_wp
      x( 6) = +0.525532409916328985817739049189_wp
      x( 7) = +0.796666477413626739591553936476_wp
      x( 8) = +0.960289856497536231683560868569_wp

      w( 1) = +0.101228536290376259152531354_wp
      w( 2) = +0.22238103445337447054435599_wp
      w( 3) = + 0.313706645877887287337962202_wp
      w( 4) = +0.3626837833783619829651504493_wp
      w( 5) = +0.3626837833783619829651504493_wp
      w( 6) = +0.313706645877887287337962202_wp
      w( 7) = +0.22238103445337447054435599_wp
      w( 8) = +0.101228536290376259152531354_wp
    case(09) !     P =  9 ;
      x( 1) = -0.968160239507626089835576202904_wp
      x( 2) = -0.836031107326635794299429788070_wp
      x( 3) = -0.613371432700590397308702039341_wp
      x( 4) = -0.324253423403808929038538014643_wp
      x( 5) = +0.000000000000000000000000000000_wp
      x( 6) = +0.324253423403808929038538014643_wp
      x( 7) = +0.613371432700590397308702039341_wp
      x( 8) = +0.836031107326635794299429788070_wp
      x( 9) = +0.968160239507626089835576202904_wp

      w( 1) = +0.08127438836157441197189216_wp
      w( 2) = +0.18064816069485740405847203_wp
      w( 3) = +0.26061069640293546231874287_wp
      w( 4) = +0.312347077040002840068630407_wp
      w( 5) = +0.330239355001259763164525069287_wp
      w( 6) = +0.312347077040002840068630407_wp
      w( 7) = +0.26061069640293546231874287_wp
      w( 8) = +0.18064816069485740405847203_wp
      w( 9) = +0.08127438836157441197189216_wp
    case(10) !     P = 10 ;
      x( 1) = -0.973906528517171720077964012084_wp
      x( 2) = -0.865063366688984510732096688423_wp
      x( 3) = -0.679409568299024406234327365115_wp
      x( 4) = -0.433395394129247190799265943166_wp
      x( 5) = -0.148874338981631210884826001130_wp
      x( 6) = +0.148874338981631210884826001130_wp
      x( 7) = +0.433395394129247190799265943166_wp
      x( 8) = +0.679409568299024406234327365115_wp
      x( 9) = +0.865063366688984510732096688423_wp
      x(10) = +0.973906528517171720077964012084_wp

      w( 1) = +0.066671344308688137593568810_wp
      w( 2) = +0.149451349150580593145776340_wp
      w( 3) = +0.219086362515982043995534934_wp
      w( 4) = +0.2692667193099963550912269216_wp
      w( 5) = +0.29552422471475287017389299465_wp
      w( 6) = +0.29552422471475287017389299465_wp
      w( 7) = +0.2692667193099963550912269216_wp
      w( 8) = +0.219086362515982043995534934_wp
      w( 9) = +0.149451349150580593145776340_wp
      w(10) = +0.066671344308688137593568810_wp

    case(11) ! P = 11
      x( 1) = -0.978228658146056992803938001123_wp
      x( 2) = -0.887062599768095299075157769304_wp
      x( 3) = -0.730152005574049324093416252031_wp
      x( 4) = -0.519096129206811815925725669459_wp
      x( 5) = -0.269543155952344972331531985401_wp
      x( 6) = +0.0000000000000000000000000000000_wp
      x( 7) = +0.269543155952344972331531985401_wp
      x( 8) = +0.519096129206811815925725669459_wp
      x( 9) = +0.730152005574049324093416252031_wp
      x(10) = +0.887062599768095299075157769304_wp
      x(11) = +0.978228658146056992803938001123_wp

      w( 1) = +0.05566856711617366648275372_wp
      w( 2) = +0.1255803694649046246346943_wp
      w( 3) = +0.18629021092773425142609764_wp
      w( 4) = +0.233193764591990479918523705_wp
      w( 5) = +0.2628045445102466621806888699_wp
      w( 6) = +0.272925086777900630714483528336_wp
      w( 7) = +0.2628045445102466621806888699_wp
      w( 8) = +0.233193764591990479918523705_wp
      w( 9) = +0.18629021092773425142609764_wp
      w(10) = +0.1255803694649046246346943_wp
      w(11) = +0.05566856711617366648275372_wp
    
    case(12) ! P = 12
      x( 1) = -0.981560634246719250690549090149_wp
      x( 2) = -0.904117256370474856678465866119_wp
      x( 3) = -0.769902674194304687036893833213_wp
      x( 4) = -0.587317954286617447296702418941_wp
      x( 5) = -0.367831498998180193752691536644_wp
      x( 6) = -0.125233408511468915472441369464_wp
      x( 7) = +0.125233408511468915472441369464_wp
      x( 8) = +0.367831498998180193752691536644_wp
      x( 9) = +0.587317954286617447296702418941_wp
      x(10) = +0.769902674194304687036893833213_wp
      x(11) = +0.904117256370474856678465866119_wp
      x(12) = +0.981560634246719250690549090149_wp

      w( 1) = +0.047175336386511827194615961_wp
      w( 2) = +0.10693932599531843096025472_wp
      w( 3) = +0.16007832854334622633465253_wp
      w( 4) = +0.203167426723065921749064456_wp
      w( 5) = +0.2334925365383548087608498989_wp
      w( 6) = +0.24914704581340278500056243604_wp
      w( 7) = +0.24914704581340278500056243604_wp
      w( 8) = +0.2334925365383548087608498989_wp
      w( 9) = +0.203167426723065921749064456_wp
      w(10) = +0.16007832854334622633465253_wp
      w(11) = +0.10693932599531843096025472_wp
      w(12) = +0.047175336386511827194615961_wp

    case(13)
      x( 1) = -0.984183054718588149472829448807_wp
      x( 2) = -0.917598399222977965206547836501_wp
      x( 3) = -0.801578090733309912794206489583_wp
      x( 4) = -0.642349339440340220643984606996_wp
      x( 5) = -0.448492751036446852877912852128_wp
      x( 6) = -0.230458315955134794065528121098_wp
      x( 7) = +0.000000000000000000000000000000_wp
      x( 8) = +0.230458315955134794065528121098_wp
      x( 9) = +0.448492751036446852877912852128_wp
      x(10) = +0.642349339440340220643984606996_wp
      x(11) = +0.801578090733309912794206489583_wp
      x(12) = +0.917598399222977965206547836501_wp
      x(13) = +0.984183054718588149472829448807_wp

      w( 1) = +0.04048400476531587952002159_wp
      w( 2) = +0.09212149983772844791442178_wp
      w( 3) = +0.13887351021978723846360178_wp
      w( 4) = +0.17814598076194573828004669_wp
      w( 5) = +0.207816047536888502312523219_wp
      w( 6) = +0.226283180262897238412090186_wp
      w( 7) = +0.232551553230873910194589515269_wp
      w( 8) = +0.226283180262897238412090186_wp
      w( 9) = +0.207816047536888502312523219_wp
      w(10) = +0.17814598076194573828004669_wp
      w(11) = +0.13887351021978723846360178_wp
      w(12) = +0.09212149983772844791442178_wp
      w(13) = +0.04048400476531587952002159_wp

    case(14)
      x( 1) = -0.986283808696812338841597266704_wp
      x( 2) = -0.928434883663573517336391139378_wp
      x( 3) = -0.82720131506976499318979474265_wp
      x( 4) = -0.687292904811685470148019803019_wp
      x( 5) = -0.515248636358154091965290718551_wp
      x( 6) = -0.319112368927889760435671824168_wp
      x( 7) = -0.10805494870734366206624465022_wp
      x( 8) = +0.10805494870734366206624465022_wp
      x( 9) = +0.319112368927889760435671824168_wp
      x(10) = +0.515248636358154091965290718551_wp
      x(11) = +0.687292904811685470148019803019_wp
      x(12) = +0.82720131506976499318979474265_wp
      x(13) = +0.928434883663573517336391139378_wp
      x(14) = +0.986283808696812338841597266704_wp

      w( 1) = +0.03511946033175186303183288_wp
      w( 2) = +0.0801580871597602098056333_wp
      w( 3) = +0.1215185706879031846894148_wp
      w( 4) = +0.15720316715819353456960194_wp
      w( 5) = +0.18553839747793781374171659_wp
      w( 6) = +0.2051984637212956039659240657_wp
      w( 7) = +0.21526385346315779019587644332_wp
      w( 8) = +0.21526385346315779019587644332_wp
      w( 9) = +0.2051984637212956039659240657_wp
      w(10) = +0.18553839747793781374171659_wp
      w(11) = +0.15720316715819353456960194_wp
      w(12) = +0.1215185706879031846894148_wp
      w(13) = +0.0801580871597602098056333_wp
      w(14) = +0.03511946033175186303183288_wp

    case(15)
      x( 1) = -0.987992518020485428489565718587_wp
      x( 2) = -0.93727339240070590430775894771_wp
      x( 3) = -0.848206583410427216200648320774_wp
      x( 4) = -0.724417731360170047416186054614_wp
      x( 5) = -0.570972172608538847537226737254_wp
      x( 6) = -0.394151347077563369897207370981_wp
      x( 7) = -0.201194093997434522300628303395_wp
      x( 8) = +0.000000000000000000000000000000_wp
      x( 9) = +0.201194093997434522300628303395_wp
      x(10) = +0.394151347077563369897207370981_wp
      x(11) = +0.570972172608538847537226737254_wp
      x(12) = +0.724417731360170047416186054614_wp
      x(13) = +0.848206583410427216200648320774_wp
      x(14) = +0.93727339240070590430775894771_wp
      x(15) = +0.987992518020485428489565718587_wp

      w( 1) = +0.03075324199611726835462839_wp
      w( 2) = +0.0703660474881081247092674_wp
      w( 3) = +0.1071592204671719350118695_wp
      w( 4) = +0.1395706779261543144478048_wp
      w( 5) = +0.16626920581699393355320086_wp
      w( 6) = +0.186161000015562211026800562_wp
      w( 7) = +0.1984314853271115764561183264_wp
      w( 8) = +0.202578241925561272880620199968_wp
      w( 9) = +0.1984314853271115764561183264_wp
      w(10) = +0.186161000015562211026800562_wp
      w(11) = +0.16626920581699393355320086_wp
      w(12) = +0.1395706779261543144478048_wp
      w(13) = +0.1071592204671719350118695_wp
      w(14) = +0.0703660474881081247092674_wp
      w(15) = +0.03075324199611726835462839_wp

    case(16)
      x( 1) = -0.98940093499164993259615417345_wp
      x( 2) = -0.944575023073232576077988415535_wp
      x( 3) = -0.865631202387831743880467897712_wp
      x( 4) = -0.755404408355003033895101194847_wp
      x( 5) = -0.617876244402643748446671764049_wp
      x( 6) = -0.458016777657227386342419442984_wp
      x( 7) = -0.28160355077925891323046050146_wp
      x( 8) = -0.095012509837637440185319335425_wp
      x( 9) = +0.095012509837637440185319335425_wp
      x(10) = +0.28160355077925891323046050146_wp
      x(11) = +0.458016777657227386342419442984_wp
      x(12) = +0.617876244402643748446671764049_wp
      x(13) = +0.755404408355003033895101194847_wp
      x(14) = +0.865631202387831743880467897712_wp
      x(15) = +0.944575023073232576077988415535_wp
      x(16) = +0.98940093499164993259615417345_wp

      w( 1) = +0.0271524594117540948517806_wp
      w( 2) = +0.0622535239386478928628438_wp
      w( 3) = +0.0951585116824927848099251_wp
      w( 4) = +0.1246289712555338720524763_wp
      w( 5) = +0.1495959888165767320815017_wp
      w( 6) = +0.16915651939500253818931208_wp
      w( 7) = +0.182603415044923588866763668_wp
      w( 8) = +0.18945061045506849628539672321_wp
      w( 9) = +0.18945061045506849628539672321_wp
      w(10) = +0.182603415044923588866763668_wp
      w(11) = +0.16915651939500253818931208_wp
      w(12) = +0.1495959888165767320815017_wp
      w(13) = +0.1246289712555338720524763_wp
      w(14) = +0.0951585116824927848099251_wp
      w(15) = +0.0622535239386478928628438_wp
      w(16) = +0.0271524594117540948517806_wp
    
    case(17)
      
      x( 1) = -0.990575475314417335675434019941_wp
      x( 2) = -0.950675521768767761222716957896_wp
      x( 3) = -0.880239153726985902122955694488_wp
      x( 4) = -0.78151400389680140692523005552_wp
      x( 5) = -0.657671159216690765850302216643_wp
      x( 6) = -0.51269053708647696788624656863_wp
      x( 7) = -0.351231763453876315297185517095_wp
      x( 8) = -0.178484181495847855850677493654_wp
      x( 9) = +0.0000000000000000000000000000000_wp
      x(10) = +0.178484181495847855850677493654_wp
      x(11) = +0.351231763453876315297185517095_wp
      x(12) = +0.51269053708647696788624656863_wp
      x(13) = +0.657671159216690765850302216643_wp
      x(14) = +0.78151400389680140692523005552_wp
      x(15) = +0.880239153726985902122955694488_wp
      x(16) = +0.950675521768767761222716957896_wp
      x(17) = +0.990575475314417335675434019941_wp
      
      w( 1) = +0.02414830286854793196011_wp
      w( 2) = +0.05545952937398720112944_wp
      w( 3) = +0.085036148317179180883535_wp
      w( 4) = +0.111883847193403971094788_wp
      w( 5) = +0.13513636846852547328632_wp
      w( 6) = +0.15404576107681028808143159_wp
      w( 7) = +0.168004102156450044509970664_wp
      w( 8) = +0.1765627053669926463252709901_wp
      w( 9) = +0.179446470356206525458265644262_wp
      w(10) = +0.1765627053669926463252709901_wp
      w(11) = +0.168004102156450044509970664_wp
      w(12) = +0.15404576107681028808143159_wp
      w(13) = +0.13513636846852547328632_wp
      w(14) = +0.111883847193403971094788_wp
      w(15) = +0.085036148317179180883535_wp
      w(16) = +0.05545952937398720112944_wp
      w(17) = +0.02414830286854793196011_wp

    case(18)

      x( 1) = -0.991565168420930946730016004706_wp
      x( 2) = -0.95582394957139775518119589293_wp
      x( 3) = -0.892602466497555739206060591127_wp
      x( 4) = -0.803704958972523115682417455015_wp
      x( 5) = -0.691687043060353207874891081289_wp
      x( 6) = -0.559770831073947534607871548525_wp
      x( 7) = -0.411751161462842646035931793833_wp
      x( 8) = -0.251886225691505509588972854878_wp
      x( 9) = -0.0847750130417353012422618529358_wp
      x(10) = +0.0847750130417353012422618529358_wp
      x(11) = +0.251886225691505509588972854878_wp
      x(12) = +0.411751161462842646035931793833_wp
      x(13) = +0.559770831073947534607871548525_wp
      x(14) = +0.691687043060353207874891081289_wp
      x(15) = +0.803704958972523115682417455015_wp
      x(16) = +0.892602466497555739206060591127_wp
      x(17) = +0.95582394957139775518119589293_wp
      x(18) = +0.991565168420930946730016004706_wp
      
      w( 1) = 0.0216160135264833103133427_wp
      w( 2) = 0.049714548894969796453335_wp
      w( 3) = 0.07642573025488905652913_wp
      w( 4) = 0.100942044106287165562814_wp
      w( 5) = 0.1225552067114784601845191_wp
      w( 6) = 0.1406429146706506512047313_wp
      w( 7) = 0.154684675126265244925418_wp
      w( 8) = 0.1642764837458327229860537765_wp
      w( 9) = 0.16914238296314359184065647013_wp
      w(10) = 0.16914238296314359184065647013_wp
      w(11) = 0.1642764837458327229860537765_wp
      w(12) = 0.154684675126265244925418_wp
      w(13) = 0.1406429146706506512047313_wp
      w(14) = 0.1225552067114784601845191_wp
      w(15) = 0.100942044106287165562814_wp
      w(16) = 0.07642573025488905652913_wp
      w(17) = 0.049714548894969796453335_wp
      w(18) = 0.0216160135264833103133427_wp

    case default
      write(*,*)'only defined for p <= 18'
      write(*,*)'stopping'
      stop
    end select

    return
  end subroutine Gauss_Legendre_points

  subroutine Gauss_Lobatto_Legendre_points(p,x,w)
    ! Coefficients for Gauss_Lobatto_Legendre polynomials
    !
    ! Mathematica script that generates the coefficients
    !  (*  Gauss_Lobatto_Legendre Points  *)
    !  For[i=11,i<=11,i++,ans=N[Solve[(1-x^2)*D[LegendreP[i,x],x]==0,x],30] ;
    !  Print[i];
    !  Print[Re[x/.ans]];
    !  Print[Re[N[2/(i (i+1) (LegendreP[i,x])^2)/.ans,30]]] ] ;


    implicit none 

    integer,                    intent( in) :: p
    real(wp), dimension(p+1),   intent(out) :: x,w

    continue

    write(*,*) p

    select case(p)

    case(02) !       P = 1 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) =  1.000000000000000000000000000000_wp

      w( 1) =  1.00000000000000000000000000000_wp
      w( 2) =  1.00000000000000000000000000000_wp
    case(03) !       P = 2 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) =  0.000000000000000000000000000000_wp
      x( 3) =  1.000000000000000000000000000000_wp

      w( 1) =  0.33333333333333333333333333333_wp
      w( 2) =  1.33333333333333333333333333333_wp
      w( 3) =  0.33333333333333333333333333333_wp
    case(04) !       P = 3 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.447213595499957939281834733746_wp
      x( 3) =  0.447213595499957939281834733746_wp
      x( 4) =  1.000000000000000000000000000000_wp

      w( 1) =  0.16666666666666666666666666667_wp
      w( 2) =  0.83333333333333333333333333333_wp
      w( 3) =  0.83333333333333333333333333333_wp
      w( 4) =  0.16666666666666666666666666667_wp
    case(05) !       P = 4 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.654653670707977143798292456247_wp
      x( 3) =  0.000000000000000000000000000000_wp
      x( 4) =  0.654653670707977143798292456247_wp
      x( 5) =  1.000000000000000000000000000000_wp

      w( 1) =  0.10000000000000000000000000000_wp
      w( 2) =  0.5444444444444444444444444444_wp
      w( 3) =  0.711111111111111111111111111111_wp
      w( 4) =  0.5444444444444444444444444444_wp
      w( 5) =  0.10000000000000000000000000000_wp
    case(06) !       P = 5 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.765055323929464692851002973959_wp
      x( 3) = -0.285231516480645096314150994041_wp
      x( 4) =  0.285231516480645096314150994041_wp
      x( 5) =  0.765055323929464692851002973959_wp
      x( 6) =  1.000000000000000000000000000000_wp

      w( 1) =  0.06666666666666666666666666667_wp
      w( 2) =  0.3784749562978469803166128082_wp
      w( 3) =  0.55485837703548635301672052512_wp
      w( 4) =  0.55485837703548635301672052512_wp
      w( 5) =  0.3784749562978469803166128082_wp
      w( 6) =  0.06666666666666666666666666667_wp
    case(07) !       P = 6 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.830223896278566929872032213967_wp
      x( 3) = -0.468848793470714213803771881909_wp
      x( 4) =  0.000000000000000000000000000000_wp
      x( 5) =  0.468848793470714213803771881909_wp
      x( 6) =  0.830223896278566929872032213967_wp
      x( 7) =  1.000000000000000000000000000000_wp

      w( 1) =  0.0476190476190476190476190476_wp
      w( 2) =  0.2768260473615659480107004063_wp
      w( 3) =  0.4317453812098626234178710223_wp
      w( 4) =  0.487619047619047619047619047619_wp
      w( 5) =  0.4317453812098626234178710223_wp
      w( 6) =  0.2768260473615659480107004063_wp
      w( 7) =  0.0476190476190476190476190476_wp
    case(08) !       P = 7 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.871740148509606615337445761221_wp
      x( 3) = -0.591700181433142302144510731398_wp
      x( 4) = -0.209299217902478868768657260345_wp
      x( 5) =  0.209299217902478868768657260345_wp
      x( 6) =  0.591700181433142302144510731398_wp
      x( 7) =  0.871740148509606615337445761221_wp
      x( 8) =  1.000000000000000000000000000000_wp

      w( 1) =  0.0357142857142857142857142857_wp
      w( 2) =  0.21070422714350603938299207_wp
      w( 3) =  0.341122692483504364764240677_wp
      w( 4) =  0.4124587946587038815670529714_wp
      w( 5) =  0.4124587946587038815670529714_wp
      w( 6) =  0.341122692483504364764240677_wp
      w( 7) =  0.21070422714350603938299207_wp
      w( 8) =  0.0357142857142857142857142857_wp
    case(09) !       P = 8 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.899757995411460157312345244418_wp
      x( 3) = -0.677186279510737753445885427091_wp
      x( 4) = -0.363117463826178158710752068709_wp
      x( 5) =  0.000000000000000000000000000000_wp
      x( 6) =  0.363117463826178158710752068709_wp
      x( 7) =  0.677186279510737753445885427091_wp
      x( 8) =  0.899757995411460157312345244418_wp
      x( 9) =  1.000000000000000000000000000000_wp

      w( 1) =  0.0277777777777777777777777778_wp
      w( 2) =  0.16549536156080552504633972_wp
      w( 3) =  0.27453871250016173528070562_wp
      w( 4) =  0.346428510973046345115131532_wp
      w( 5) =  0.371519274376417233560090702948_wp
      w( 6) =  0.346428510973046345115131532_wp
      w( 7) =  0.27453871250016173528070562_wp
      w( 8) =  0.16549536156080552504633972_wp
      w( 9) =  0.0277777777777777777777777778_wp
    case(10) !       P = 9 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.919533908166458813828932660822_wp
      x( 3) = -0.738773865105505075003106174860_wp
      x( 4) = -0.477924949810444495661175092731_wp
      x( 5) = -0.165278957666387024626219765958_wp
      x( 6) =  0.165278957666387024626219765958_wp
      x( 7) =  0.477924949810444495661175092731_wp
      x( 8) =  0.738773865105505075003106174860_wp
      x( 9) =  0.919533908166458813828932660822_wp
      x(10) =  1.000000000000000000000000000000_wp

      w( 1) =  0.022222222222222222222222222_wp
      w( 2) =  0.13330599085107011112622717_wp
      w( 3) =  0.22488934206312645211945782_wp
      w( 4) =  0.292042683679683757875582257_wp
      w( 5) =  0.3275397611838974566565105279_wp
      w( 6) =  0.3275397611838974566565105279_wp
      w( 7) =  0.292042683679683757875582257_wp
      w( 8) =  0.22488934206312645211945782_wp
      w( 9) =  0.13330599085107011112622717_wp
      w(10) =  0.022222222222222222222222222_wp
    case(11) !       P = 10 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.934001430408059134332274136099_wp
      x( 3) = -0.784483473663144418622417816108_wp
      x( 4) = -0.565235326996205006470963969478_wp
      x( 5) = -0.295758135586939391431911515559_wp
      x( 6) =  0.000000000000000000000000000000_wp
      x( 7) =  0.295758135586939391431911515559_wp
      x( 8) =  0.565235326996205006470963969478_wp
      x( 9) =  0.784483473663144418622417816108_wp
      x(10) =  0.934001430408059134332274136099_wp
      x(11) =  1.000000000000000000000000000000_wp

      w( 1) =  0.018181818181818181818181818_wp
      w( 2) =  0.1096122732669948644614034_wp
      w( 3) =  0.1871698817803052041081415_wp
      w( 4) =  0.24804810426402831404008487_wp
      w( 5) =  0.286879124779008088679222403_wp
      w( 6) =  0.300217595455690693785931881170_wp
      w( 7) =  0.286879124779008088679222403_wp
      w( 8) =  0.24804810426402831404008487_wp
      w( 9) =  0.1871698817803052041081415_wp
      w(10) =  0.1096122732669948644614034_wp
      w(11) =  0.018181818181818181818181818_wp
    case(12) !       P = 11 ;
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.944899272222882223407580138303_wp
      x( 3) = -0.819279321644006678348641581717_wp
      x( 4) = -0.632876153031860677662404854444_wp
      x( 5) = -0.399530940965348932264349791567_wp
      x( 6) = -0.136552932854927554864061855740_wp
      x( 7) =  0.136552932854927554864061855740_wp
      x( 8) =  0.399530940965348932264349791567_wp
      x( 9) =  0.632876153031860677662404854444_wp
      x(10) =  0.819279321644006678348641581717_wp
      x(11) =  0.944899272222882223407580138303_wp
      x(12) =  1.000000000000000000000000000000_wp

      w( 1) = 0.015151515151515151515151515_wp
      w( 2) = 0.09168451741319613066834259_wp
      w( 3) = 0.15797470556437011516467106_wp
      w( 4) = 0.212508417761021145358302077_wp
      w( 5) = 0.2512756031992012802932444121_wp
      w( 6) = 0.27140524091069617700028833850_wp
      w( 7) = 0.27140524091069617700028833850_wp
      w( 8) = 0.2512756031992012802932444121_wp
      w( 9) = 0.212508417761021145358302077_wp
      w(10) = 0.15797470556437011516467106_wp
      w(11) = 0.09168451741319613066834259_wp
      w(12) = 0.015151515151515151515151515_wp

    case(13) ! P = 12
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.953309846642163911896905464755_wp
      x( 3) = -0.846347564651872316865925607099_wp
      x( 4) = -0.686188469081757426072759039566_wp
      x( 5) = -0.482909821091336201746937233637_wp
      x( 6) = -0.249286930106239992568673700374_wp
      x( 7) =  0.000000000000000000000000000000_wp
      x( 8) =  0.249286930106239992568673700374_wp
      x( 9) =  0.482909821091336201746937233637_wp
      x(10) =  0.686188469081757426072759039566_wp
      x(11) =  0.846347564651872316865925607099_wp
      x(12) =  0.953309846642163911896905464755_wp
      x(13) =  1.000000000000000000000000000000_wp
      
      w( 1) = 0.01282051282051282051282051_wp
      w( 2) = 0.077801686746818927793589_wp
      w( 3) = 0.1349819266896083491199148_wp
      w( 4) = 0.18364686520355009200749426_wp
      w( 5) = 0.220767793566110086085534008_wp
      w( 6) = 0.2440157903066763564585781484_wp
      w( 7) = 0.251930849333446736044138641541_wp
      w( 8) = 0.2440157903066763564585781484_wp
      w( 9) = 0.220767793566110086085534008_wp
      w(10) = 0.18364686520355009200749426_wp
      w(11) = 0.1349819266896083491199148_wp
      w(12) = 0.077801686746818927793589_wp
      w(13) = 0.01282051282051282051282051_wp

    case(14) ! P = 13
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.959935045267260901355100162015_wp
      x( 3) = -0.867801053830347251000220202908_wp
      x( 4) = -0.728868599091326140584672400521_wp
      x( 5) = -0.550639402928647055316622705859_wp
      x( 6) = -0.342724013342712845043903403642_wp
      x( 7) = -0.116331868883703867658776709736_wp
      x( 8) =  0.116331868883703867658776709736_wp
      x( 9) =  0.342724013342712845043903403642_wp
      x(10) =  0.550639402928647055316622705859_wp
      x(11) =  0.728868599091326140584672400521_wp
      x(12) =  0.867801053830347251000220202908_wp
      x(13) =  0.959935045267260901355100162015_wp
      x(14) =  1.000000000000000000000000000000_wp

      w( 1) = 0.01098901098901098901098901_wp
      w( 2) = 0.0668372844976812846340707_wp
      w( 3) = 0.1165866558987116515409967_wp
      w( 4) = 0.160021851762952142412821_wp
      w( 5) = 0.19482614937341611864033178_wp
      w( 6) = 0.219126253009770754871162524_wp
      w( 7) = 0.23161279446845705888962835729_wp
      w( 8) = 0.23161279446845705888962835729_wp
      w( 9) = 0.219126253009770754871162524_wp
      w(10) = 0.19482614937341611864033178_wp
      w(11) = 0.160021851762952142412821_wp
      w(12) = 0.1165866558987116515409967_wp
      w(13) = 0.0668372844976812846340707_wp
      w(14) = 0.01098901098901098901098901_wp

    case(15) ! P = 14
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.96524592650383857279585139207_wp
      x( 3) = -0.885082044222976298825401631482_wp
      x( 4) = -0.763519689951815200704118475976_wp
      x( 5) = -0.606253205469845711123529938637_wp
      x( 6) = -0.420638054713672480921896938739_wp
      x( 7) = -0.215353955363794238225679446273_wp
      x( 8) =  0.000000000000000000000000000000_wp
      x( 9) =  0.215353955363794238225679446273_wp
      x(10) =  0.420638054713672480921896938739_wp
      x(11) =  0.606253205469845711123529938637_wp
      x(12) =  0.763519689951815200704118475976_wp
      x(13) =  0.885082044222976298825401631482_wp
      x(14) =  0.96524592650383857279585139207_wp
      x(15) =  1.000000000000000000000000000000_wp

      w( 1) = 0.00952380952380952380952381_wp
      w( 2) = 0.0580298930286012490968806_wp
      w( 3) = 0.1016600703257180676036662_wp
      w( 4) = 0.1405116998024281094604468_wp
      w( 5) = 0.1727896472536009490520771_wp
      w( 6) = 0.196987235964613356092500347_wp
      w( 7) = 0.211973585926820920127430077_wp
      w( 8) = 0.217048116348815649514950214251_wp
      w( 9) = 0.211973585926820920127430077_wp
      w(10) = 0.196987235964613356092500347_wp
      w(11) = 0.1727896472536009490520771_wp
      w(12) = 0.1405116998024281094604468_wp
      w(13) = 0.1016600703257180676036662_wp
      w(14) = 0.0580298930286012490968806_wp
      w(15) = 0.00952380952380952380952381_wp

    case(16) !  P = 15
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.969568046270217932952242738367_wp
      x( 3) = -0.89920053309347209299462826152_wp
      x( 4) = -0.792008291861815063931088270963_wp
      x( 5) = -0.652388702882493089467883219641_wp
      x( 6) = -0.486059421887137611781890785847_wp
      x( 7) = -0.299830468900763208098353454722_wp
      x( 8) = -0.101326273521949447843033005046_wp
      x( 9) =  0.101326273521949447843033005046_wp
      x(10) =  0.299830468900763208098353454722_wp
      x(11) =  0.486059421887137611781890785847_wp
      x(12) =  0.652388702882493089467883219641_wp
      x(13) =  0.792008291861815063931088270963_wp
      x(14) =  0.89920053309347209299462826152_wp
      x(15) =  0.969568046270217932952242738367_wp
      x(16) =  1.000000000000000000000000000000_wp

      w( 1) = 0.0083333333333333333333333_wp
      w( 2) = 0.050850361005919905403245_wp
      w( 3) = 0.089393697325930800991052_wp
      w( 4) = 0.1242553821325140983495363_wp
      w( 5) = 0.1540269808071642808156449_wp
      w( 6) = 0.17749191339170412530107567_wp
      w( 7) = 0.1936900238252035843169135989_wp
      w( 8) = 0.20195830817822987148919912541_wp
      w( 9) = 0.20195830817822987148919912541_wp
      w(10) = 0.1936900238252035843169135989_wp
      w(11) = 0.17749191339170412530107567_wp
      w(12) = 0.1540269808071642808156449_wp
      w(13) = 0.1242553821325140983495363_wp
      w(14) = 0.089393697325930800991052_wp
      w(15) = 0.050850361005919905403245_wp
      w(16) = 0.0083333333333333333333333_wp

    case(17) ! P = 16
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.973132176631418314156979501874_wp
      x( 3) = -0.910879995915573595623802506398_wp
      x( 4) = -0.815696251221770307106750553238_wp
      x( 5) = -0.691028980627684705394919357372_wp
      x( 6) = -0.541385399330101539123733407504_wp
      x( 7) = -0.372174433565477041907234680735_wp
      x( 8) = -0.189511973518317388304263014753_wp
      x( 9) = +0.000000000000000000000000000000_wp
      x(10) = +0.189511973518317388304263014753_wp
      x(11) = +0.372174433565477041907234680735_wp
      x(12) = +0.541385399330101539123733407504_wp
      x(13) = +0.691028980627684705394919357372_wp
      x(14) = +0.815696251221770307106750553238_wp
      x(15) = +0.910879995915573595623802506398_wp
      x(16) = +0.973132176631418314156979501874_wp
      x(17) = +1.000000000000000000000000000000_wp

      w( 1) = 0.0073529411764705882352941_wp
      w( 2) = 0.044921940543254209647401_wp
      w( 3) = 0.079198270503687119190264_wp
      w( 4) = 0.110592909007028161375773_wp
      w( 5) = 0.1379877462019265590562016_wp
      w( 6) = 0.16039466199762153951632837_wp
      w( 7) = 0.177004253515657870436945745_wp
      w( 8) = 0.1872163396776192358920884829_wp
      w( 9) = 0.190661874753469433299407247028_wp
      w(10) = 0.1872163396776192358920884829_wp
      w(11) = 0.177004253515657870436945745_wp
      w(12) = 0.16039466199762153951632837_wp
      w(13) = 0.1379877462019265590562016_wp
      w(14) = 0.110592909007028161375773_wp
      w(15) = 0.079198270503687119190264_wp
      w(16) = 0.044921940543254209647401_wp
      w(17) = 0.0073529411764705882352941_wp
   
    case(18) ! P = 17
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.976105557412198542864518924342_wp
      x( 3) = -0.920649185347533873837854625431_wp
      x( 4) = -0.835593535218090213713646362328_wp
      x( 5) = -0.723679329283242681306210365302_wp
      x( 6) = -0.588504834318661761173535893194_wp
      x( 7) = -0.434415036912123975342287136741_wp
      x( 8) = -0.266362652878280984167665332026_wp
      x( 9) = -0.0897490934846521110226450100886_wp
      x(10) = +0.0897490934846521110226450100886_wp
      x(11) = +0.266362652878280984167665332026_wp
      x(12) = +0.434415036912123975342287136741_wp
      x(13) = +0.588504834318661761173535893194_wp
      x(14) = +0.723679329283242681306210365302_wp
      x(15) = +0.835593535218090213713646362328_wp
      x(16) = +0.920649185347533873837854625431_wp
      x(17) = +0.976105557412198542864518924342_wp
      x(18) = +1.000000000000000000000000000000_wp

      w( 1) = 0.0065359477124183006535948_wp
      w( 2) = 0.039970628810914066137599_wp
      w( 3) = 0.070637166885633664999223_wp
      w( 4) = 0.099016271717502802394424_wp
      w( 5) = 0.124210533132967100263396_wp
      w( 6) = 0.1454119615738022679830032_wp
      w( 7) = 0.16193951723760248926432671_wp
      w( 8) = 0.1732621094894562260106144038_wp
      w( 9) = 0.17901586343970308229381880694_wp
      w(10) = 0.17901586343970308229381880694_wp
      w(11) = 0.1732621094894562260106144038_wp
      w(12) = 0.16193951723760248926432671_wp
      w(13) = 0.1454119615738022679830032_wp
      w(14) = 0.124210533132967100263396_wp
      w(15) = 0.099016271717502802394424_wp
      w(16) = 0.070637166885633664999223_wp
      w(17) = 0.039970628810914066137599_wp
      w(18) = 0.0065359477124183006535948_wp
    
    case(19) ! P = 18
      
      x( 1) = -1.000000000000000000000000000000_wp
      x( 2) = -0.97861176622208009515263406311_wp
      x( 3) = -0.928901528152586243717940258797_wp
      x( 4) = -0.852460577796646093085955970041_wp
      x( 5) = -0.751494202552613014163637489634_wp
      x( 6) = -0.628908137265220497766832306229_wp
      x( 7) = -0.488229285680713502777909637625_wp
      x( 8) = -0.333504847824498610298500103845_wp
      x( 9) = -0.169186023409281571375154153445_wp
      x(10) = 0.000000000000000000000000000000_wp
      x(11) = 0.169186023409281571375154153445_wp
      x(12) = 0.333504847824498610298500103845_wp
      x(13) = 0.488229285680713502777909637625_wp
      x(14) = 0.628908137265220497766832306229_wp
      x(15) = 0.751494202552613014163637489634_wp
      x(16) = 0.852460577796646093085955970041_wp
      x(17) = 0.928901528152586243717940258797_wp
      x(18) = 0.97861176622208009515263406311_wp
      x(19) = 1.000000000000000000000000000000_wp

      w( 1) = 0.005847953216374269005848_wp
      w( 2) = 0.03579336518617647711543_wp
      w( 3) = 0.0633818917626297368517_wp
      w( 4) = 0.08913175709920708444801_wp
      w( 5) = 0.11231534147730504407091_wp
      w( 6) = 0.1322672804487507769260467_wp
      w( 7) = 0.14841394259593888500968064_wp
      w( 8) = 0.160290924044061241979910968_wp
      w( 9) = 0.1675565845271428672701372777_wp
      w(10) = 0.170001919284827234644672715617_wp
      w(11) = 0.1675565845271428672701372777_wp
      w(12) = 0.160290924044061241979910968_wp
      w(13) = 0.14841394259593888500968064_wp
      w(14) = 0.1322672804487507769260467_wp
      w(15) = 0.11231534147730504407091_wp
      w(16) = 0.08913175709920708444801_wp
      w(17) = 0.0633818917626297368517_wp
      w(18) = 0.03579336518617647711543_wp
      w(19) = 0.005847953216374269005848_wp

    case default
      write(*,*) 'LGL points only defined for p <= 18'
      write(*,*) 'Stopping'
      stop

    end select


    return
  end subroutine Gauss_Lobatto_Legendre_points

  subroutine Interpolate_GLL_2_GL(N_Soln_Pts,N_Flux_Pts,Int_F2S)
    !  Interpolate from Gauss_Lobatto_Legendre points
    !                to Gauss_        Legendre points
    !

    implicit none

    integer,                                    intent( in) :: N_Soln_Pts,N_Flux_Pts
    real(wp), dimension(N_Soln_Pts,N_Flux_Pts), intent(out) :: Int_F2S

    continue

    select case(N_Soln_Pts)
    case(2)

      Int_F2S(1,1) = +0.4553418012614795489212410569_wp
      Int_F2S(1,2) = +0.6666666666666666666666666667_wp
      Int_F2S(1,3) = -0.1220084679281462155879077236_wp
      Int_F2S(2,1) = -0.1220084679281462155879077236_wp
      Int_F2S(2,2) = +0.6666666666666666666666666667_wp
      Int_F2S(2,3) = +0.4553418012614795489212410569_wp

    case(3)

      Int_F2S(1,1) = +0.44364916731037084425896327_wp
      Int_F2S(1,2) = +0.6830127018922193233818615854_wp
      Int_F2S(1,3) = -0.1830127018922193233818615854_wp
      Int_F2S(1,4) = +0.05635083268962915574103673001_wp
      Int_F2S(2,1) = -0.1250000000000000000000000000_wp
      Int_F2S(2,2) = +0.6250000000000000000000000000_wp
      Int_F2S(2,3) = +0.6250000000000000000000000000_wp
      Int_F2S(2,4) = -0.1250000000000000000000000000_wp
      Int_F2S(3,1) = +0.05635083268962915574103673001_wp
      Int_F2S(3,2) = -0.1830127018922193233818615854_wp
      Int_F2S(3,3) = +0.6830127018922193233818615854_wp
      Int_F2S(3,4) = +0.44364916731037084425896327_wp

    case(4)

      Int_F2S(1,1) = +0.4389152966531085059730026146_wp
      Int_F2S(1,2) = +0.6887516501894637412985684673_wp
      Int_F2S(1,3) = -0.188740996194412276919596025_wp
      Int_F2S(1,4) = +0.09382253564559162985375479359_wp
      Int_F2S(1,5) = -0.03274848629375160020572985054_wp
      Int_F2S(2,1) = -0.1247624770988934354482036499_wp
      Int_F2S(2,2) = +0.6106019927354806547773027584_wp
      Int_F2S(2,3) = +0.6458838533372694197767388821_wp
      Int_F2S(2,4) = -0.1931761785705360259296260193_wp
      Int_F2S(2,5) = +0.06145280959667938682378802869_wp
      Int_F2S(3,1) = +0.06145280959667938682378802869_wp
      Int_F2S(3,2) = -0.1931761785705360259296260193_wp
      Int_F2S(3,3) = +0.6458838533372694197767388821_wp
      Int_F2S(3,4) = +0.6106019927354806547773027584_wp
      Int_F2S(3,5) = -0.1247624770988934354482036499_wp
      Int_F2S(4,1) = -0.03274848629375160020572985054_wp
      Int_F2S(4,2) = +0.09382253564559162985375479359_wp
      Int_F2S(4,3) = -0.188740996194412276919596025_wp
      Int_F2S(4,4) = +0.6887516501894637412985684673_wp
      Int_F2S(4,5) = +0.4389152966531085059730026146_wp

    case(5)

      Int_F2S(1,1) = +0.4365363738624496765746326981_wp
      Int_F2S(1,2) = +0.6914778893211668643913368159_wp
      Int_F2S(1,3) = -0.1902820486169034492210305453_wp
      Int_F2S(1,4) = +0.09917256452430220669718635795_wp
      Int_F2S(1,5) = -0.05839063727822323222977856403_wp
      Int_F2S(1,6) = +0.02148585818720793378765323735_wp
      Int_F2S(2,1) = -0.1244051944023349039114759978_wp
      Int_F2S(2,2) = +0.6037681555875741277558571256_wp
      Int_F2S(2,3) = +0.6541047177124134806412658927_wp
      Int_F2S(2,4) = -0.2010973282614030882103135116_wp
      Int_F2S(2,5) = +0.1049503907147767938791801325_wp
      Int_F2S(2,6) = -0.03732074135102641015451364136_wp
      Int_F2S(3,1) = +0.06250000000000000000000000000_wp
      Int_F2S(3,2) = -0.1946486423538422797658774615_wp
      Int_F2S(3,3) = +0.6321486423538422797658774615_wp
      Int_F2S(3,4) = +0.6321486423538422797658774615_wp
      Int_F2S(3,5) = -0.1946486423538422797658774615_wp
      Int_F2S(3,6) = +0.06250000000000000000000000000_wp
      Int_F2S(4,1) = -0.03732074135102641015451364136_wp
      Int_F2S(4,2) = +0.1049503907147767938791801325_wp
      Int_F2S(4,3) = -0.2010973282614030882103135116_wp
      Int_F2S(4,4) = +0.6541047177124134806412658927_wp
      Int_F2S(4,5) = +0.6037681555875741277558571256_wp
      Int_F2S(4,6) = -0.1244051944023349039114759978_wp
      Int_F2S(5,1) = +0.02148585818720793378765323735_wp
      Int_F2S(5,2) = -0.05839063727822323222977856403_wp
      Int_F2S(5,3) = +0.09917256452430220669718635795_wp
      Int_F2S(5,4) = -0.1902820486169034492210305453_wp
      Int_F2S(5,5) = +0.6914778893211668643913368159_wp
      Int_F2S(5,6) = +0.4365363738624496765746326981_wp

    case(6)

      Int_F2S(1,1) = +0.4351734281450880317824983441_wp
      Int_F2S(1,2) = +0.692995895766239736084917899_wp
      Int_F2S(1,3) = -0.1908635085315759876631423831_wp
      Int_F2S(1,4) = +0.1008503894174153852632139378_wp
      Int_F2S(1,5) = -0.06314645066888346283078122458_wp
      Int_F2S(1,6) = +0.04019745757287234107353289041_wp
      Int_F2S(1,7) = -0.01520721170115604371023946375_wp
      Int_F2S(2,1) = -0.124137237053192688765576656_wp
      Int_F2S(2,2) = +0.5999603880851736977319801094_wp
      Int_F2S(2,3) = +0.6583257378714987411110181478_wp
      Int_F2S(2,4) = -0.2035374890409166302752860395_wp
      Int_F2S(2,5) = +0.112061424419986019727185912_wp
      Int_F2S(2,6) = -0.0679896392766522253366170211_wp
      Int_F2S(2,7) = +0.02531681499410308580729554735_wp
      Int_F2S(3,1) = +0.06278424196534005266671621489_wp
      Int_F2S(3,2) = -0.194820244505407087643623716_wp
      Int_F2S(3,3) = +0.6251941541027388521887783395_wp
      Int_F2S(3,4) = +0.6410578198124032520958973673_wp
      Int_F2S(3,5) = -0.2034554337438759609226168988_wp
      Int_F2S(3,6) = +0.107833017018865265641115289_wp
      Int_F2S(3,7) = -0.0385935546500643740262665959_wp
      Int_F2S(4,1) = -0.0385935546500643740262665959_wp
      Int_F2S(4,2) = +0.107833017018865265641115289_wp
      Int_F2S(4,3) = -0.2034554337438759609226168988_wp
      Int_F2S(4,4) = +0.6410578198124032520958973673_wp
      Int_F2S(4,5) = +0.6251941541027388521887783395_wp
      Int_F2S(4,6) = -0.194820244505407087643623716_wp
      Int_F2S(4,7) = +0.06278424196534005266671621489_wp
      Int_F2S(5,1) = +0.02531681499410308580729554735_wp
      Int_F2S(5,2) = -0.0679896392766522253366170211_wp
      Int_F2S(5,3) = +0.112061424419986019727185912_wp
      Int_F2S(5,4) = -0.2035374890409166302752860395_wp
      Int_F2S(5,5) = +0.6583257378714987411110181478_wp
      Int_F2S(5,6) = +0.5999603880851736977319801094_wp
      Int_F2S(5,7) = -0.124137237053192688765576656_wp
      Int_F2S(6,1) = -0.01520721170115604371023946375_wp
      Int_F2S(6,2) = +0.04019745757287234107353289041_wp
      Int_F2S(6,3) = -0.06314645066888346283078122458_wp
      Int_F2S(6,4) = +0.1008503894174153852632139378_wp
      Int_F2S(6,5) = -0.1908635085315759876631423831_wp
      Int_F2S(6,6) = +0.692995895766239736084917899_wp
      Int_F2S(6,7) = +0.4351734281450880317824983441_wp

    case(7)

      Int_F2S(1,1) = +0.4343202773431659577187751906_wp
      Int_F2S(1,2) = +0.6939304865766944381798985655_wp
      Int_F2S(1,3) = -0.1911308324656548154046745156_wp
      Int_F2S(1,4) = +0.1015337067167920367717938512_wp
      Int_F2S(1,5) = -0.06484379890853231006885509719_wp
      Int_F2S(1,6) = +0.04433494178435378867428544948_wp
      Int_F2S(1,7) = -0.02948507959360312159185627151_wp
      Int_F2S(1,8) = +0.01134029854678402572063282758_wp
      Int_F2S(2,1) = -0.1239475950210044312855623761_wp
      Int_F2S(2,2) = +0.597614347259035618477522008_wp
      Int_F2S(2,3) = +0.660813468654549139195979822_wp
      Int_F2S(2,4) = -0.2045571191722765042771766864_wp
      Int_F2S(2,5) = +0.114501847692819048533023147_wp
      Int_F2S(2,6) = -0.07426343845739941112028028205_wp
      Int_F2S(2,7) = +0.0482341331750241442665210601_wp
      Int_F2S(2,8) = -0.01839564413074760379002669252_wp
      Int_F2S(3,1) = +0.06286527783698134326235830845_wp
      Int_F2S(3,2) = -0.1947323378660034080705554958_wp
      Int_F2S(3,3) = +0.6211128101503285960333783877_wp
      Int_F2S(3,4) = +0.6458262927229572593083844631_wp
      Int_F2S(3,5) = -0.206349172495029830176072195_wp
      Int_F2S(3,6) = +0.1157209965318159678767380533_wp
      Int_F2S(3,7) = -0.07101273159581146501042112764_wp
      Int_F2S(3,8) = +0.02656886471476153677618960593_wp
      Int_F2S(4,1) = -0.03906250000000000000000000000000_wp
      Int_F2S(4,2) = +0.1088400233988091816746151422_wp
      Int_F2S(4,3) = -0.2040293534695560887945905658_wp
      Int_F2S(4,4) = +0.6342518300707469071199754236_wp
      Int_F2S(4,5) = +0.6342518300707469071199754236_wp
      Int_F2S(4,6) = -0.2040293534695560887945905658_wp
      Int_F2S(4,7) = +0.1088400233988091816746151422_wp
      Int_F2S(4,8) = -0.03906250000000000000000000000000_wp
      Int_F2S(5,1) = +0.02656886471476153677618960593_wp
      Int_F2S(5,2) = -0.07101273159581146501042112764_wp
      Int_F2S(5,3) = +0.1157209965318159678767380533_wp
      Int_F2S(5,4) = -0.206349172495029830176072195_wp
      Int_F2S(5,5) = +0.6458262927229572593083844631_wp
      Int_F2S(5,6) = +0.6211128101503285960333783877_wp
      Int_F2S(5,7) = -0.1947323378660034080705554958_wp
      Int_F2S(5,8) = +0.06286527783698134326235830845_wp
      Int_F2S(6,1) = -0.01839564413074760379002669252_wp
      Int_F2S(6,2) = +0.0482341331750241442665210601_wp
      Int_F2S(6,3) = -0.07426343845739941112028028205_wp
      Int_F2S(6,4) = +0.114501847692819048533023147_wp
      Int_F2S(6,5) = -0.2045571191722765042771766864_wp
      Int_F2S(6,6) = +0.660813468654549139195979822_wp
      Int_F2S(6,7) = +0.597614347259035618477522008_wp
      Int_F2S(6,8) = -0.1239475950210044312855623761_wp
      Int_F2S(7,1) = +0.01134029854678402572063282758_wp
      Int_F2S(7,2) = -0.02948507959360312159185627151_wp
      Int_F2S(7,3) = +0.04433494178435378867428544948_wp
      Int_F2S(7,4) = -0.06484379890853231006885509719_wp
      Int_F2S(7,5) = +0.1015337067167920367717938512_wp
      Int_F2S(7,6) = -0.1911308324656548154046745156_wp
      Int_F2S(7,7) = +0.6939304865766944381798985655_wp
      Int_F2S(7,8) = +0.4343202773431659577187751906_wp

    case(8)

      Int_F2S(1,1) = +0.4337509518771576115969178631_wp
      Int_F2S(1,2) = +0.6945476305066941768399758643_wp
      Int_F2S(1,3) = -0.1912712021344270397406367577_wp
      Int_F2S(1,4) = +0.1018591949239170430002731713_wp
      Int_F2S(1,5) = -0.06559662289565308967195199342_wp
      Int_F2S(1,6) = +0.04596279483584462362322306265_wp
      Int_F2S(1,7) = -0.0330689164306377696949946885_wp
      Int_F2S(1,8) = +0.02260278446296199910060669063_wp
      Int_F2S(1,9) = -0.008786615145857555053413212347_wp
      Int_F2S(2,1) = -0.1238120368927635408664401186_wp
      Int_F2S(2,2) = +0.596064259301189625067849999_wp
      Int_F2S(2,3) = +0.6624131821149836839257083416_wp
      Int_F2S(2,4) = -0.2050650812627149573887624234_wp
      Int_F2S(2,5) = +0.1155679035853609635927666476_wp
      Int_F2S(2,6) = -0.07665717772194716717884773957_wp
      Int_F2S(2,7) = +0.05369956918745465138463441063_wp
      Int_F2S(2,8) = -0.03622275574301365905980214147_wp
      Int_F2S(2,9) = +0.01401213743145040052289302419_wp
      Int_F2S(3,1) = +0.06288225176545602955088532873_wp
      Int_F2S(3,2) = -0.1946011148368713594654555393_wp
      Int_F2S(3,3) = +0.6184921438459408613138572686_wp
      Int_F2S(3,4) = +0.6487340730633449342734310311_wp
      Int_F2S(3,5) = -0.2076237629864822564623020735_wp
      Int_F2S(3,6) = +0.1185665047806636646389569115_wp
      Int_F2S(3,7) = -0.07798725317277651670212834283_wp
      Int_F2S(3,8) = +0.05109465121326286490318156413_wp
      Int_F2S(3,9) = -0.01955749367253822205042614849_wp
      Int_F2S(4,1) = -0.03926405018612098137025949562_wp
      Int_F2S(4,2) = +0.1092499132738985359039482809_wp
      Int_F2S(4,3) = -0.2041411446146674369463796936_wp
      Int_F2S(4,4) = +0.6301407320459061771255630881_wp
      Int_F2S(4,5) = +0.6392144955921210112214494326_wp
      Int_F2S(4,6) = -0.2071631657066965952773878022_wp
      Int_F2S(4,7) = +0.1171189565094929775691855647_wp
      Int_F2S(4,8) = -0.07224778072100559641827932285_wp
      Int_F2S(4,9) = +0.02709204380707190819215994807_wp
      Int_F2S(5,1) = +0.02709204380707190819215994807_wp
      Int_F2S(5,2) = -0.07224778072100559641827932285_wp
      Int_F2S(5,3) = +0.1171189565094929775691855647_wp
      Int_F2S(5,4) = -0.2071631657066965952773878022_wp
      Int_F2S(5,5) = +0.6392144955921210112214494326_wp
      Int_F2S(5,6) = +0.6301407320459061771255630881_wp
      Int_F2S(5,7) = -0.2041411446146674369463796936_wp
      Int_F2S(5,8) = +0.1092499132738985359039482809_wp
      Int_F2S(5,9) = -0.03926405018612098137025949562_wp
      Int_F2S(6,1) = -0.01955749367253822205042614849_wp
      Int_F2S(6,2) = +0.05109465121326286490318156413_wp
      Int_F2S(6,3) = -0.07798725317277651670212834283_wp
      Int_F2S(6,4) = +0.1185665047806636646389569115_wp
      Int_F2S(6,5) = -0.2076237629864822564623020735_wp
      Int_F2S(6,6) = +0.6487340730633449342734310311_wp
      Int_F2S(6,7) = +0.6184921438459408613138572686_wp
      Int_F2S(6,8) = -0.1946011148368713594654555393_wp
      Int_F2S(6,9) = +0.06288225176545602955088532873_wp
      Int_F2S(7,1) = +0.01401213743145040052289302419_wp
      Int_F2S(7,2) = -0.03622275574301365905980214147_wp
      Int_F2S(7,3) = +0.05369956918745465138463441063_wp
      Int_F2S(7,4) = -0.07665717772194716717884773957_wp
      Int_F2S(7,5) = +0.1155679035853609635927666476_wp
      Int_F2S(7,6) = -0.2050650812627149573887624234_wp
      Int_F2S(7,7) = +0.6624131821149836839257083416_wp
      Int_F2S(7,8) = +0.596064259301189625067849999_wp
      Int_F2S(7,9) = -0.1238120368927635408664401186_wp
      Int_F2S(8,1) = -0.008786615145857555053413212347_wp
      Int_F2S(8,2) = +0.02260278446296199910060669063_wp
      Int_F2S(8,3) = -0.0330689164306377696949946885_wp
      Int_F2S(8,4) = +0.04596279483584462362322306265_wp
      Int_F2S(8,5) = -0.06559662289565308967195199342_wp
      Int_F2S(8,6) = +0.1018591949239170430002731713_wp
      Int_F2S(8,7) = -0.1912712021344270397406367577_wp
      Int_F2S(8,8) = +0.6945476305066941768399758643_wp
      Int_F2S(8,9) = +0.4337509518771576115969178631_wp

    case(9)

      Int_F2S(1,1) = +0.4333521587716858022453300286_wp
      Int_F2S(1,2) = +0.6949768560346215209878625803_wp
      Int_F2S(1,3) = -0.1913523645716244138581456982_wp
      Int_F2S(1,4) = +0.1020318140037896033797638652_wp
      Int_F2S(1,5) = -0.06597777939430220461145271576_wp
      Int_F2S(1,6) = +0.04673592039626910885491974758_wp
      Int_F2S(1,7) = -0.03458966059950067883084088777_wp
      Int_F2S(1,8) = +0.02571489140894872722846530336_wp
      Int_F2S(1,9) = -0.01790235717879486844964868094_wp
      Int_F2S(1,10) = +0.007010521128907403053746457642_wp
      Int_F2S(2,1) = -0.1237129041364125172952043243_wp
      Int_F2S(2,2) = +0.5949856802639829704280503694_wp
      Int_F2S(2,3) = +0.6635065413526627978184683299_wp
      Int_F2S(2,4) = -0.2053497125975989814115764656_wp
      Int_F2S(2,5) = +0.1161054653559125387356589479_wp
      Int_F2S(2,6) = -0.07777609873054110745472761113_wp
      Int_F2S(2,7) = +0.05596610033185344509561958376_wp
      Int_F2S(2,8) = -0.04097702098817907869019484586_wp
      Int_F2S(2,9) = +0.02830027388514305481181679986_wp
      Int_F2S(2,10) = -0.01104832473682312203791078397_wp
      Int_F2S(3,1) = +0.06287769045394333004646401915_wp
      Int_F2S(3,2) = -0.1944776290939593190706628057_wp
      Int_F2S(3,3) = +0.6167018550728947525177940489_wp
      Int_F2S(3,4) = +0.6506571572740953762611349526_wp
      Int_F2S(3,5) = -0.2082865093109048262173217336_wp
      Int_F2S(3,6) = +0.1198633155881047701430730178_wp
      Int_F2S(3,7) = -0.0807564516225138238941444646_wp
      Int_F2S(3,8) = +0.05719497218253428870285905079_wp
      Int_F2S(3,9) = -0.03884241952765412675806244985_wp
      Int_F2S(3,10) = +0.01506801898345957826886636452_wp
      Int_F2S(4,1) = -0.03935930221376297983881370487_wp
      Int_F2S(4,2) = +0.1094310811419253885279512941_wp
      Int_F2S(4,3) = -0.2041155004643057128209608407_wp
      Int_F2S(4,4) = +0.6274328626107246445527163505_wp
      Int_F2S(4,5) = +0.6423062402219947517647948098_wp
      Int_F2S(4,6) = -0.2085874097968236254873724543_wp
      Int_F2S(4,7) = +0.1201959176844107347940043709_wp
      Int_F2S(4,8) = -0.07959348581660181416519864295_wp
      Int_F2S(4,9) = +0.05237405573831092901643564479_wp
      Int_F2S(4,10) = -0.02008445910587231634355682737_wp
      Int_F2S(5,1) = +0.02734375000000000000000000000_wp
      Int_F2S(5,2) = -0.07283186257766806981930472846_wp
      Int_F2S(5,3) = +0.1177435823540570945963627472_wp
      Int_F2S(5,4) = -0.2074090936574362986857453131_wp
      Int_F2S(5,5) = +0.6351536238810472739086872944_wp
      Int_F2S(5,6) = +0.6351536238810472739086872944_wp
      Int_F2S(5,7) = -0.2074090936574362986857453131_wp
      Int_F2S(5,8) = +0.1177435823540570945963627472_wp
      Int_F2S(5,9) = -0.07283186257766806981930472846_wp
      Int_F2S(5,10) = +0.027343750000000000000000000000000_wp
      Int_F2S(6,1) = -0.02008445910587231634355682737_wp
      Int_F2S(6,2) = +0.05237405573831092901643564479_wp
      Int_F2S(6,3) = -0.07959348581660181416519864295_wp
      Int_F2S(6,4) = +0.1201959176844107347940043709_wp
      Int_F2S(6,5) = -0.2085874097968236254873724543_wp
      Int_F2S(6,6) = +0.6423062402219947517647948098_wp
      Int_F2S(6,7) = +0.6274328626107246445527163505_wp
      Int_F2S(6,8) = -0.2041155004643057128209608407_wp
      Int_F2S(6,9) = +0.1094310811419253885279512941_wp
      Int_F2S(6,10) = -0.03935930221376297983881370487_wp
      Int_F2S(7,1) = +0.01506801898345957826886636452_wp
      Int_F2S(7,2) = -0.03884241952765412675806244985_wp
      Int_F2S(7,3) = +0.05719497218253428870285905079_wp
      Int_F2S(7,4) = -0.0807564516225138238941444646_wp
      Int_F2S(7,5) = +0.1198633155881047701430730178_wp
      Int_F2S(7,6) = -0.2082865093109048262173217336_wp
      Int_F2S(7,7) = +0.6506571572740953762611349526_wp
      Int_F2S(7,8) = +0.6167018550728947525177940489_wp
      Int_F2S(7,9) = -0.1944776290939593190706628057_wp
      Int_F2S(7,10) = +0.06287769045394333004646401915_wp
      Int_F2S(8,1) = -0.01104832473682312203791078397_wp
      Int_F2S(8,2) = +0.02830027388514305481181679986_wp
      Int_F2S(8,3) = -0.04097702098817907869019484586_wp
      Int_F2S(8,4) = +0.05596610033185344509561958376_wp
      Int_F2S(8,5) = -0.07777609873054110745472761113_wp
      Int_F2S(8,6) = +0.1161054653559125387356589479_wp
      Int_F2S(8,7) = -0.2053497125975989814115764656_wp
      Int_F2S(8,8) = +0.6635065413526627978184683299_wp
      Int_F2S(8,9) = +0.5949856802639829704280503694_wp
      Int_F2S(8,10) = -0.1237129041364125172952043243_wp
      Int_F2S(9,1) = +0.007010521128907403053746457642_wp
      Int_F2S(9,2) = -0.01790235717879486844964868094_wp
      Int_F2S(9,3) = +0.02571489140894872722846530336_wp
      Int_F2S(9,4) = -0.03458966059950067883084088777_wp
      Int_F2S(9,5) = +0.04673592039626910885491974758_wp
      Int_F2S(9,6) = -0.06597777939430220461145271576_wp
      Int_F2S(9,7) = +0.1020318140037896033797638652_wp
      Int_F2S(9,8) = -0.1913523645716244138581456982_wp
      Int_F2S(9,9) = +0.6949768560346215209878625803_wp
      Int_F2S(9,10) = +0.4333521587716858022453300286_wp

    case(10)

      Int_F2S(1,1) = +0.4330619903678455256956031589_wp
      Int_F2S(1,2) = +0.695287598942437838401077166_wp
      Int_F2S(1,3) = -0.1914028590539965493874019236_wp
      Int_F2S(1,4) = +0.1021309085527247673737346675_wp
      Int_F2S(1,5) = -0.06618933391656590937161960221_wp
      Int_F2S(1,6) = +0.04714808609989238730270669018_wp
      Int_F2S(1,7) = -0.0353527917202562025323458273_wp
      Int_F2S(1,8) = +0.02711768311748524238594897763_wp
      Int_F2S(1,9) = -0.02061892653213852798212302305_wp
      Int_F2S(1,10) = +0.01454237858805149377737680863_wp
      Int_F2S(1,11) = -0.005724734445480065662957092683_wp
      Int_F2S(2,1) = -0.1236386497461782664234301893_wp
      Int_F2S(2,2) = +0.5942045592480217613641910321_wp
      Int_F2S(2,3) = +0.6642885472665671207844642283_wp
      Int_F2S(2,4) = -0.2055233542848130309956076597_wp
      Int_F2S(2,5) = +0.1164044994692509488090925847_wp
      Int_F2S(2,6) = -0.07836743250010144799583639749_wp
      Int_F2S(2,7) = +0.05708861383229933205217754244_wp
      Int_F2S(2,8) = -0.04308307397494072617744033502_wp
      Int_F2S(2,9) = +0.03245030620991884458156895327_wp
      Int_F2S(2,10) = -0.02276922534068969833259629152_wp
      Int_F2S(2,11) = +0.008945209820665162333416532264_wp
      Int_F2S(3,1) = +0.06286626714457784285860074426_wp
      Int_F2S(3,2) = -0.1943724234822238916607768267_wp
      Int_F2S(3,3) = +0.6154215269733453439403814119_wp
      Int_F2S(3,4) = +0.652002703500085656470772885_wp
      Int_F2S(3,5) = -0.2086710608907948062319904681_wp
      Int_F2S(3,6) = +0.1205413322331536787583869081_wp
      Int_F2S(3,7) = -0.08209557305608433674301913416_wp
      Int_F2S(3,8) = +0.05980976122659611270052386311_wp
      Int_F2S(3,9) = -0.04417313384970009427103380731_wp
      Int_F2S(3,10) = +0.03067143912905072378171832985_wp
      Int_F2S(3,11) = -0.01200083892800622960356390602_wp
      Int_F2S(4,1) = -0.03940705002597900068285893215_wp
      Int_F2S(4,2) = +0.1095138081036284164507524665_wp
      Int_F2S(4,3) = -0.2040503149156609544224848585_wp
      Int_F2S(4,4) = +0.6255418258514654735423135907_wp
      Int_F2S(4,5) = +0.6443900566575508980619143569_wp
      Int_F2S(4,6) = -0.2093481609035134862643443416_wp
      Int_F2S(4,7) = +0.1216370451704319797806665226_wp
      Int_F2S(4,8) = -0.08258447350066153374684293213_wp
      Int_F2S(4,9) = +0.05882328291144842827438571281_wp
      Int_F2S(4,10) = -0.04009317003577707783431604898_wp
      Int_F2S(4,11) = +0.01557715068706685684081446388_wp
      Int_F2S(5,1) = +0.02747673920879407320106116098_wp
      Int_F2S(5,2) = -0.07313577199799650719777553962_wp
      Int_F2S(5,3) = +0.1180505660639145872467343837_wp
      Int_F2S(5,4) = -0.2074618023807564574135738611_wp
      Int_F2S(5,5) = +0.6324348668489196986257613241_wp
      Int_F2S(5,6) = +0.638319990900739500110398289_wp
      Int_F2S(5,7) = -0.2089240882337212907095570151_wp
      Int_F2S(5,8) = +0.1209604142470982279433964763_wp
      Int_F2S(5,9) = -0.08039148237472842917410027053_wp
      Int_F2S(5,10) = +0.05302628202713949614514501513_wp
      Int_F2S(5,11) = -0.020355714309402898777489963_wp
      Int_F2S(6,1) = -0.020355714309402898777489963_wp
      Int_F2S(6,2) = +0.05302628202713949614514501513_wp
      Int_F2S(6,3) = -0.08039148237472842917410027053_wp
      Int_F2S(6,4) = +0.1209604142470982279433964763_wp
      Int_F2S(6,5) = -0.2089240882337212907095570151_wp
      Int_F2S(6,6) = +0.638319990900739500110398289_wp
      Int_F2S(6,7) = +0.6324348668489196986257613241_wp
      Int_F2S(6,8) = -0.2074618023807564574135738611_wp
      Int_F2S(6,9) = +0.1180505660639145872467343837_wp
      Int_F2S(6,10) = -0.07313577199799650719777553962_wp
      Int_F2S(6,11) = +0.02747673920879407320106116098_wp
      Int_F2S(7,1) = +0.01557715068706685684081446388_wp
      Int_F2S(7,2) = -0.04009317003577707783431604898_wp
      Int_F2S(7,3) = +0.05882328291144842827438571281_wp
      Int_F2S(7,4) = -0.08258447350066153374684293213_wp
      Int_F2S(7,5) = +0.1216370451704319797806665226_wp
      Int_F2S(7,6) = -0.2093481609035134862643443416_wp
      Int_F2S(7,7) = +0.6443900566575508980619143569_wp
      Int_F2S(7,8) = +0.6255418258514654735423135907_wp
      Int_F2S(7,9) = -0.2040503149156609544224848585_wp
      Int_F2S(7,10) = +0.1095138081036284164507524665_wp
      Int_F2S(7,11) = -0.03940705002597900068285893215_wp
      Int_F2S(8,1) = -0.01200083892800622960356390602_wp
      Int_F2S(8,2) = +0.03067143912905072378171832985_wp
      Int_F2S(8,3) = -0.04417313384970009427103380731_wp
      Int_F2S(8,4) = +0.05980976122659611270052386311_wp
      Int_F2S(8,5) = -0.08209557305608433674301913416_wp
      Int_F2S(8,6) = +0.1205413322331536787583869081_wp
      Int_F2S(8,7) = -0.2086710608907948062319904681_wp
      Int_F2S(8,8) = +0.652002703500085656470772885_wp
      Int_F2S(8,9) = +0.6154215269733453439403814119_wp
      Int_F2S(8,10) = -0.1943724234822238916607768267_wp
      Int_F2S(8,11) = +0.06286626714457784285860074426_wp
      Int_F2S(9,1) = +0.008945209820665162333416532264_wp
      Int_F2S(9,2) = -0.02276922534068969833259629152_wp
      Int_F2S(9,3) = +0.03245030620991884458156895327_wp
      Int_F2S(9,4) = -0.04308307397494072617744033502_wp
      Int_F2S(9,5) = +0.05708861383229933205217754244_wp
      Int_F2S(9,6) = -0.07836743250010144799583639749_wp
      Int_F2S(9,7) = +0.1164044994692509488090925847_wp
      Int_F2S(9,8) = -0.2055233542848130309956076597_wp
      Int_F2S(9,9) = +0.6642885472665671207844642283_wp
      Int_F2S(9,10) = +0.5942045592480217613641910321_wp
      Int_F2S(9,11) = -0.1236386497461782664234301893_wp
      Int_F2S(10,1) = -0.005724734445480065662957092683_wp
      Int_F2S(10,2) = +0.01454237858805149377737680863_wp
      Int_F2S(10,3) = -0.02061892653213852798212302305_wp
      Int_F2S(10,4) = +0.02711768311748524238594897763_wp
      Int_F2S(10,5) = -0.0353527917202562025323458273_wp
      Int_F2S(10,6) = +0.04714808609989238730270669018_wp
      Int_F2S(10,7) = -0.06618933391656590937161960221_wp
      Int_F2S(10,8) = +0.1021309085527247673737346675_wp
      Int_F2S(10,9) = -0.1914028590539965493874019236_wp
      Int_F2S(10,10) = +0.695287598942437838401077166_wp
      Int_F2S(10,11) = +0.4330619903678455256956031589_wp

    case default
      write(*,*)'Int_F2S only defined for p <= 10'
      write(*,*)'stopping'
      stop
    end select

    return
  end subroutine Interpolate_GLL_2_GL

  subroutine Extrapolate_GL_2_GLL(N_Soln_Pts,N_Flux_Pts,Ext_S2F)
    ! Interpolate from Gauss_Legendre points
    ! to Gauss_Lobatto_Legendre points
    !
    ! Prolongation operator used to take data from 
    ! the solutions points to the flux points

    implicit none

    integer,                    intent( in) :: N_Soln_Pts,N_Flux_Pts
    real(wp), dimension(N_Flux_Pts,N_Soln_Pts),   intent(out) :: Ext_S2F

    continue

    select case(N_Soln_Pts)

    case(2)

      Ext_S2F(1,1) = +1.366025403784438646763723171_wp
      Ext_S2F(1,2) = -0.3660254037844386467637231708_wp
      Ext_S2F(2,1) = +0.5000000000000000000000000000_wp
      Ext_S2F(2,2) = +0.5000000000000000000000000000_wp
      Ext_S2F(3,1) = -0.3660254037844386467637231708_wp
      Ext_S2F(3,2) = +1.366025403784438646763723171_wp

    case(3)

      Ext_S2F(1,1) = +1.478830557701236147529877567_wp
      Ext_S2F(1,2) = -0.6666666666666666666666666667_wp
      Ext_S2F(1,3) = +0.1878361089654305191367891_wp
      Ext_S2F(2,1) = +0.4553418012614795489212410569_wp
      Ext_S2F(2,2) = +0.6666666666666666666666666667_wp
      Ext_S2F(2,3) = -0.1220084679281462155879077236_wp
      Ext_S2F(3,1) = -0.1220084679281462155879077236_wp
      Ext_S2F(3,2) = +0.6666666666666666666666666667_wp
      Ext_S2F(3,3) = +0.4553418012614795489212410569_wp
      Ext_S2F(4,1) = +0.1878361089654305191367891_wp
      Ext_S2F(4,2) = -0.6666666666666666666666666667_wp
      Ext_S2F(4,3) = +1.478830557701236147529877567_wp

    case(4)

      Ext_S2F(1,1) = +1.526788125457266786984328278_wp
      Ext_S2F(1,2) = -0.8136324494869272605618980818_wp
      Ext_S2F(1,3) = +0.4007615203116504048002817771_wp
      Ext_S2F(1,4) = -0.1139171962819899312227119735_wp
      Ext_S2F(2,1) = +0.4400551811292643496078056205_wp
      Ext_S2F(2,2) = +0.7313898326544354607127893658_wp
      Ext_S2F(2,3) = -0.2313898326544354607127893658_wp
      Ext_S2F(2,4) = +0.05994481887073565039219437954_wp
      Ext_S2F(3,1) = -0.09232659844072882091060611425_wp
      Ext_S2F(3,2) = +0.5923265984407288209106061143_wp
      Ext_S2F(3,3) = +0.5923265984407288209106061143_wp
      Ext_S2F(3,4) = -0.09232659844072882091060611425_wp
      Ext_S2F(4,1) = +0.05994481887073565039219437954_wp
      Ext_S2F(4,2) = -0.2313898326544354607127893658_wp
      Ext_S2F(4,3) = +0.7313898326544354607127893658_wp
      Ext_S2F(4,4) = +0.4400551811292643496078056205_wp
      Ext_S2F(5,1) = -0.1139171962819899312227119735_wp
      Ext_S2F(5,2) = +0.4007615203116504048002817771_wp
      Ext_S2F(5,3) = -0.8136324494869272605618980818_wp
      Ext_S2F(5,4) = +1.526788125457266786984328278_wp

    case(5)

      Ext_S2F(1,1) = +1.551408049094313012813027985_wp
      Ext_S2F(1,2) = -0.8931583920000717373261767813_wp
      Ext_S2F(1,3) = +0.5333333333333333333333333333_wp
      Ext_S2F(1,4) = -0.2679416522233875093041099283_wp
      Ext_S2F(1,5) = +0.07635866179581290048392539171_wp
      Ext_S2F(2,1) = +0.4328680132620556028466721453_wp
      Ext_S2F(2,2) = +0.7635399510326310095392527997_wp
      Ext_S2F(2,3) = -0.292578010855927216089176134_wp
      Ext_S2F(2,4) = +0.1327228265446223581214880164_wp
      Ext_S2F(2,5) = -0.03655277998338175441823682734_wp
      Ext_S2F(3,1) = -0.08125124342860909342361422131_wp
      Ext_S2F(3,2) = +0.5642399653020522942338906198_wp
      Ext_S2F(3,3) = +0.6481335664114827716447316895_wp
      Ext_S2F(3,4) = -0.1734693948048192463512025327_wp
      Ext_S2F(3,5) = +0.04234710651989327389619444463_wp
      Ext_S2F(4,1) = +0.04234710651989327389619444463_wp
      Ext_S2F(4,2) = -0.1734693948048192463512025327_wp
      Ext_S2F(4,3) = +0.6481335664114827716447316895_wp
      Ext_S2F(4,4) = +0.5642399653020522942338906198_wp
      Ext_S2F(4,5) = -0.08125124342860909342361422131_wp
      Ext_S2F(5,1) = -0.03655277998338175441823682734_wp
      Ext_S2F(5,2) = +0.1327228265446223581214880164_wp
      Ext_S2F(5,3) = -0.292578010855927216089176134_wp
      Ext_S2F(5,4) = +0.7635399510326310095392527997_wp
      Ext_S2F(5,5) = +0.4328680132620556028466721453_wp
      Ext_S2F(6,1) = +0.07635866179581290048392539171_wp
      Ext_S2F(6,2) = -0.2679416522233875093041099283_wp
      Ext_S2F(6,3) = +0.5333333333333333333333333333_wp
      Ext_S2F(6,4) = -0.8931583920000717373261767813_wp
      Ext_S2F(6,5) = +1.551408049094313012813027985_wp

    case(6)

      Ext_S2F(1,1) = +1.565673200151071933093716587_wp
      Ext_S2F(1,2) = -0.9404628431763489290199939332_wp
      Ext_S2F(1,3) = +0.6169300554304887086169604506_wp
      Ext_S2F(1,4) = -0.3792277021146137546173491821_wp
      Ext_S2F(1,5) = +0.191800014038667954820175561_wp
      Ext_S2F(1,6) = -0.05471272432926591289350948294_wp
      Ext_S2F(2,1) = +0.4288872784717705086970442698_wp
      Ext_S2F(2,2) = +0.7818724265115144086727882193_wp
      Ext_S2F(2,3) = -0.3293010466673126297619681515_wp
      Ext_S2F(2,4) = +0.1822681490815033240243593419_wp
      Ext_S2F(2,5) = -0.08860455672505475974946608806_wp
      Ext_S2F(2,6) = +0.02487774932757914811724240857_wp
      Ext_S2F(3,1) = -0.07573814367451240643839685279_wp
      Ext_S2F(3,2) = +0.5500895646111271377218560971_wp
      Ext_S2F(3,3) = +0.6775684680130059693562304561_wp
      Ext_S2F(3,4) = -0.2204994810749714892041888841_wp
      Ext_S2F(3,5) = +0.09363726286655580552421633574_wp
      Ext_S2F(3,6) = -0.02505767074120501695971715203_wp
      Ext_S2F(4,1) = +0.03543368918328000190499362945_wp
      Ext_S2F(4,2) = -0.1505858006966852031663958385_wp
      Ext_S2F(4,3) = +0.615152111513405201261402209_wp
      Ext_S2F(4,4) = +0.615152111513405201261402209_wp
      Ext_S2F(4,5) = -0.1505858006966852031663958385_wp
      Ext_S2F(4,6) = +0.03543368918328000190499362945_wp
      Ext_S2F(5,1) = -0.02505767074120501695971715203_wp
      Ext_S2F(5,2) = +0.09363726286655580552421633574_wp
      Ext_S2F(5,3) = -0.2204994810749714892041888841_wp
      Ext_S2F(5,4) = +0.6775684680130059693562304561_wp
      Ext_S2F(5,5) = +0.5500895646111271377218560971_wp
      Ext_S2F(5,6) = -0.07573814367451240643839685279_wp
      Ext_S2F(6,1) = +0.02487774932757914811724240857_wp
      Ext_S2F(6,2) = -0.08860455672505475974946608806_wp
      Ext_S2F(6,3) = +0.1822681490815033240243593419_wp
      Ext_S2F(6,4) = -0.3293010466673126297619681515_wp
      Ext_S2F(6,5) = +0.7818724265115144086727882193_wp
      Ext_S2F(6,6) = +0.4288872784717705086970442698_wp
      Ext_S2F(7,1) = -0.05471272432926591289350948294_wp
      Ext_S2F(7,2) = +0.191800014038667954820175561_wp
      Ext_S2F(7,3) = -0.3792277021146137546173491821_wp
      Ext_S2F(7,4) = +0.6169300554304887086169604506_wp
      Ext_S2F(7,5) = -0.9404628431763489290199939332_wp
      Ext_S2F(7,6) = +1.565673200151071933093716587_wp

    case(7)

      Ext_S2F(1,1) = +1.574662499710550498744170263_wp
      Ext_S2F(1,2) = -0.9707266965061221906486769487_wp
      Ext_S2F(1,3) = +0.6721078619223617869349168557_wp
      Ext_S2F(1,4) = -0.4571428571428571428571428572_wp
      Ext_S2F(1,5) = +0.2840541467652299666802055797_wp
      Ext_S2F(1,6) = -0.1440701036120688469286845855_wp
      Ext_S2F(1,7) = +0.04111514886290592807521169289_wp
      Ext_S2F(2,1) = +0.426444057606557966271441516_wp
      Ext_S2F(2,2) = +0.7933203677297410295021739678_wp
      Ext_S2F(2,3) = -0.3528864105403762807549904859_wp
      Ext_S2F(2,4) = +0.2158983137048545341143881568_wp
      Ext_S2F(2,5) = -0.1286865254642843683714253231_wp
      Ext_S2F(2,6) = +0.06402978851334280215184631135_wp
      Ext_S2F(2,7) = -0.01811959154983568291343414282_wp
      Ext_S2F(3,1) = -0.07255034602202543015413657927_wp
      Ext_S2F(3,2) = +0.5418375676087438182598734699_wp
      Ext_S2F(3,3) = +0.6952323632956333509450346846_wp
      Ext_S2F(3,4) = -0.2499861307986277753659755555_wp
      Ext_S2F(3,5) = +0.1295303857640097375739371638_wp
      Ext_S2F(3,6) = -0.06089270689041354446770531462_wp
      Ext_S2F(3,7) = +0.01682886704267984320897213112_wp
      Ext_S2F(4,1) = +0.03187491392952970931288332024_wp
      Ext_S2F(4,2) = -0.1387186546654784268684681089_wp
      Ext_S2F(4,3) = +0.5978679304728864032998676313_wp
      Ext_S2F(4,4) = +0.6427099641774587864755941612_wp
      Ext_S2F(4,5) = -0.1910259060439313492214976048_wp
      Ext_S2F(4,6) = +0.07764844525055307727465976663_wp
      Ext_S2F(4,7) = -0.02035669312101820027303916573_wp
      Ext_S2F(5,1) = -0.02035669312101820027303916573_wp
      Ext_S2F(5,2) = +0.07764844525055307727465976663_wp
      Ext_S2F(5,3) = -0.1910259060439313492214976048_wp
      Ext_S2F(5,4) = +0.6427099641774587864755941612_wp
      Ext_S2F(5,5) = +0.5978679304728864032998676313_wp
      Ext_S2F(5,6) = -0.1387186546654784268684681089_wp
      Ext_S2F(5,7) = +0.03187491392952970931288332024_wp
      Ext_S2F(6,1) = +0.01682886704267984320897213112_wp
      Ext_S2F(6,2) = -0.06089270689041354446770531462_wp
      Ext_S2F(6,3) = +0.1295303857640097375739371638_wp
      Ext_S2F(6,4) = -0.2499861307986277753659755555_wp
      Ext_S2F(6,5) = +0.6952323632956333509450346846_wp
      Ext_S2F(6,6) = +0.5418375676087438182598734699_wp
      Ext_S2F(6,7) = -0.07255034602202543015413657927_wp
      Ext_S2F(7,1) = -0.01811959154983568291343414282_wp
      Ext_S2F(7,2) = +0.06402978851334280215184631135_wp
      Ext_S2F(7,3) = -0.1286865254642843683714253231_wp
      Ext_S2F(7,4) = +0.2158983137048545341143881568_wp
      Ext_S2F(7,5) = -0.3528864105403762807549904859_wp
      Ext_S2F(7,6) = +0.7933203677297410295021739678_wp
      Ext_S2F(7,7) = +0.426444057606557966271441516_wp
      Ext_S2F(8,1) = +0.04111514886290592807521169289_wp
      Ext_S2F(8,2) = -0.1440701036120688469286845855_wp
      Ext_S2F(8,3) = +0.2840541467652299666802055797_wp
      Ext_S2F(8,4) = -0.4571428571428571428571428572_wp
      Ext_S2F(8,5) = +0.6721078619223617869349168557_wp
      Ext_S2F(8,6) = -0.9707266965061221906486769487_wp
      Ext_S2F(8,7) = +1.574662499710550498744170263_wp

    case(8)

      Ext_S2F(1,1) = +1.580687063030955444233446524_wp
      Ext_S2F(1,2) = -0.9912041583117163203115854404_wp
      Ext_S2F(1,3) = +0.7101568903172424245822905786_wp
      Ext_S2F(1,4) = -0.5126556338013684782412776392_wp
      Ext_S2F(1,5) = +0.3537304181064418044090136305_wp
      Ext_S2F(1,6) = -0.2208713667044388962766533154_wp
      Ext_S2F(1,7) = +0.1121772102087164247782520863_wp
      Ext_S2F(1,8) = -0.03202042284583240317348642409_wp
      Ext_S2F(2,1) = +0.4248339008239181620443260752_wp
      Ext_S2F(2,2) = +0.8009492552175297325464600611_wp
      Ext_S2F(2,3) = -0.3688783929883293854593979579_wp
      Ext_S2F(2,4) = +0.2394216460585047339157113224_wp
      Ext_S2F(2,5) = -0.1583313163913489625334890832_wp
      Ext_S2F(2,6) = +0.09685305680621155434323960597_wp
      Ext_S2F(2,7) = -0.04867359312619565532994335557_wp
      Ext_S2F(2,8) = +0.01382544359970982047309333211_wp
      Ext_S2F(3,1) = -0.07052595114999431658169081719_wp
      Ext_S2F(3,2) = +0.5365659630759535327858573806_wp
      Ext_S2F(3,3) = +0.706731280921336700338138788_wp
      Ext_S2F(3,4) = -0.2696839436514560553695311882_wp
      Ext_S2F(3,5) = +0.1547218818991282519308350223_wp
      Ext_S2F(3,6) = -0.08911355120472101623547913679_wp
      Ext_S2F(3,7) = +0.04349756592378455137110072631_wp
      Ext_S2F(3,8) = -0.01219324581403164823923077509_wp
      Ext_S2F(4,1) = +0.02976388167620098400445055745_wp
      Ext_S2F(4,2) = -0.1316363505225929030947663617_wp
      Ext_S2F(4,3) = +0.5874579709267530432034089027_wp
      Ext_S2F(4,4) = +0.6597084752559567289011028323_wp
      Ext_S2F(4,5) = -0.2168837677479411214181349261_wp
      Ext_S2F(4,6) = +0.1073673192305480539007551673_wp
      Ext_S2F(4,7) = -0.04920813945769354446171600511_wp
      Ext_S2F(4,8) = +0.01343061063876875896489983316_wp
      Ext_S2F(5,1) = -0.01787323183289530362166581505_wp
      Ext_S2F(5,2) = +0.06917571098311708503476644465_wp
      Ext_S2F(5,3) = -0.1753151418600240654049516401_wp
      Ext_S2F(5,4) = +0.6240126627098022839918510105_wp
      Ext_S2F(5,5) = +0.6240126627098022839918510105_wp
      Ext_S2F(5,6) = -0.1753151418600240654049516401_wp
      Ext_S2F(5,7) = +0.06917571098311708503476644465_wp
      Ext_S2F(5,8) = -0.01787323183289530362166581505_wp
      Ext_S2F(6,1) = +0.01343061063876875896489983316_wp
      Ext_S2F(6,2) = -0.04920813945769354446171600511_wp
      Ext_S2F(6,3) = +0.1073673192305480539007551673_wp
      Ext_S2F(6,4) = -0.2168837677479411214181349261_wp
      Ext_S2F(6,5) = +0.6597084752559567289011028323_wp
      Ext_S2F(6,6) = +0.5874579709267530432034089027_wp
      Ext_S2F(6,7) = -0.1316363505225929030947663617_wp
      Ext_S2F(6,8) = +0.02976388167620098400445055745_wp
      Ext_S2F(7,1) = -0.01219324581403164823923077509_wp
      Ext_S2F(7,2) = +0.04349756592378455137110072631_wp
      Ext_S2F(7,3) = -0.08911355120472101623547913679_wp
      Ext_S2F(7,4) = +0.1547218818991282519308350223_wp
      Ext_S2F(7,5) = -0.2696839436514560553695311882_wp
      Ext_S2F(7,6) = +0.706731280921336700338138788_wp
      Ext_S2F(7,7) = +0.5365659630759535327858573806_wp
      Ext_S2F(7,8) = -0.07052595114999431658169081719_wp
      Ext_S2F(8,1) = +0.01382544359970982047309333211_wp
      Ext_S2F(8,2) = -0.04867359312619565532994335557_wp
      Ext_S2F(8,3) = +0.09685305680621155434323960597_wp
      Ext_S2F(8,4) = -0.1583313163913489625334890832_wp
      Ext_S2F(8,5) = +0.2394216460585047339157113224_wp
      Ext_S2F(8,6) = -0.3688783929883293854593979579_wp
      Ext_S2F(8,7) = +0.8009492552175297325464600611_wp
      Ext_S2F(8,8) = +0.4248339008239181620443260752_wp
      Ext_S2F(9,1) = -0.03202042284583240317348642409_wp
      Ext_S2F(9,2) = +0.1121772102087164247782520863_wp
      Ext_S2F(9,3) = -0.2208713667044388962766533154_wp
      Ext_S2F(9,4) = +0.3537304181064418044090136305_wp
      Ext_S2F(9,5) = -0.5126556338013684782412776392_wp
      Ext_S2F(9,6) = +0.7101568903172424245822905786_wp
      Ext_S2F(9,7) = -0.9912041583117163203115854404_wp
      Ext_S2F(9,8) = +1.580687063030955444233446524_wp

    case(9)

      Ext_S2F(1,1) = +1.584919424220149128068025478_wp
      Ext_S2F(1,2) = -1.00568288639079619501990709_wp
      Ext_S2F(1,3) = +0.7373969413834670307252522568_wp
      Ext_S2F(1,4) = -0.5532193350361340907810049376_wp
      Ext_S2F(1,5) = +0.4063492063492063492063492064_wp
      Ext_S2F(1,6) = -0.2822994943041909127772907332_wp
      Ext_S2F(1,7) = +0.1767099114311423569129627886_wp
      Ext_S2F(1,8) = -0.08981367941099661693425771149_wp
      Ext_S2F(1,9) = +0.02563991175815295059987074324_wp
      Ext_S2F(2,1) = +0.423715532505720062910591041_wp
      Ext_S2F(2,2) = +0.806288360284928822267319807_wp
      Ext_S2F(2,3) = -0.3802000947548687349690045422_wp
      Ext_S2F(2,4) = +0.2564061683483851333802961922_wp
      Ext_S2F(2,5) = -0.1804266047432212245673289268_wp
      Ext_S2F(2,6) = +0.1227167895317458342798143351_wp
      Ext_S2F(2,7) = -0.0759361971540055607499869455_wp
      Ext_S2F(2,8) = +0.03835080773093971874149604505_wp
      Ext_S2F(2,9) = -0.01091476174962405129319700581_wp
      Ext_S2F(3,1) = -0.0691542171337506724744805212_wp
      Ext_S2F(3,2) = +0.5329787316942739867125266852_wp
      Ext_S2F(3,3) = +0.7146585891936818333638018106_wp
      Ext_S2F(3,4) = -0.2834944482637475689718414771_wp
      Ext_S2F(3,5) = +0.172900877984813111263695014_wp
      Ext_S2F(3,6) = -0.110546780110469059887326981_wp
      Ext_S2F(3,7) = +0.06627980407827763444102961341_wp
      Ext_S2F(3,8) = -0.03291584832059875039695271434_wp
      Ext_S2F(3,9) = +0.009293290877519485949548570425_wp
      Ext_S2F(4,1) = +0.02839507284378788230308818698_wp
      Ext_S2F(4,2) = -0.127022692068743197886728552_wp
      Ext_S2F(4,3) = +0.5806281901680541369122801462_wp
      Ext_S2F(4,4) = +0.6710554026076937193718893093_wp
      Ext_S2F(4,5) = -0.2345364192925766462167564177_wp
      Ext_S2F(4,6) = +0.1285525906276219619237567395_wp
      Ext_S2F(4,7) = -0.07206479145855548387542401032_wp
      Ext_S2F(4,8) = +0.0346188199575036906372965471_wp
      Ext_S2F(4,9) = -0.009626173384786063169401949137_wp
      Ext_S2F(5,1) = -0.01637145867831334694615490729_wp
      Ext_S2F(5,2) = +0.0640357026803536216771136841_wp
      Ext_S2F(5,3) = -0.1657255047345989522487091287_wp
      Ext_S2F(5,4) = +0.6125133509676534148006841363_wp
      Ext_S2F(5,5) = +0.6403885815848292084879819807_wp
      Ext_S2F(5,6) = -0.1989122405227732243780003466_wp
      Ext_S2F(5,7) = +0.09537059572758986032030660447_wp
      Ext_S2F(5,8) = -0.04289588882555678506377648298_wp
      Ext_S2F(5,9) = +0.01159686180081620335055446004_wp
      Ext_S2F(6,1) = +0.01159686180081620335055446004_wp
      Ext_S2F(6,2) = -0.04289588882555678506377648298_wp
      Ext_S2F(6,3) = +0.09537059572758986032030660447_wp
      Ext_S2F(6,4) = -0.1989122405227732243780003466_wp
      Ext_S2F(6,5) = +0.6403885815848292084879819807_wp
      Ext_S2F(6,6) = +0.6125133509676534148006841363_wp
      Ext_S2F(6,7) = -0.1657255047345989522487091287_wp
      Ext_S2F(6,8) = +0.0640357026803536216771136841_wp
      Ext_S2F(6,9) = -0.01637145867831334694615490729_wp
      Ext_S2F(7,1) = -0.009626173384786063169401949137_wp
      Ext_S2F(7,2) = +0.0346188199575036906372965471_wp
      Ext_S2F(7,3) = -0.07206479145855548387542401032_wp
      Ext_S2F(7,4) = +0.1285525906276219619237567395_wp
      Ext_S2F(7,5) = -0.2345364192925766462167564177_wp
      Ext_S2F(7,6) = +0.6710554026076937193718893093_wp
      Ext_S2F(7,7) = +0.5806281901680541369122801462_wp
      Ext_S2F(7,8) = -0.127022692068743197886728552_wp
      Ext_S2F(7,9) = +0.02839507284378788230308818698_wp
      Ext_S2F(8,1) = +0.009293290877519485949548570425_wp
      Ext_S2F(8,2) = -0.03291584832059875039695271434_wp
      Ext_S2F(8,3) = +0.06627980407827763444102961341_wp
      Ext_S2F(8,4) = -0.110546780110469059887326981_wp
      Ext_S2F(8,5) = +0.172900877984813111263695014_wp
      Ext_S2F(8,6) = -0.2834944482637475689718414771_wp
      Ext_S2F(8,7) = +0.7146585891936818333638018106_wp
      Ext_S2F(8,8) = +0.5329787316942739867125266852_wp
      Ext_S2F(8,9) = -0.0691542171337506724744805212_wp
      Ext_S2F(9,1) = -0.01091476174962405129319700581_wp
      Ext_S2F(9,2) = +0.03835080773093971874149604505_wp
      Ext_S2F(9,3) = -0.0759361971540055607499869455_wp
      Ext_S2F(9,4) = +0.1227167895317458342798143351_wp
      Ext_S2F(9,5) = -0.1804266047432212245673289268_wp
      Ext_S2F(9,6) = +0.2564061683483851333802961922_wp
      Ext_S2F(9,7) = -0.3802000947548687349690045422_wp
      Ext_S2F(9,8) = +0.806288360284928822267319807_wp
      Ext_S2F(9,9) = +0.423715532505720062910591041_wp
      Ext_S2F(10,1) = +0.02563991175815295059987074324_wp
      Ext_S2F(10,2) = -0.08981367941099661693425771149_wp
      Ext_S2F(10,3) = +0.1767099114311423569129627886_wp
      Ext_S2F(10,4) = -0.2822994943041909127772907332_wp
      Ext_S2F(10,5) = +0.4063492063492063492063492064_wp
      Ext_S2F(10,6) = -0.5532193350361340907810049376_wp
      Ext_S2F(10,7) = +0.7373969413834670307252522568_wp
      Ext_S2F(10,8) = -1.00568288639079619501990709_wp
      Ext_S2F(10,9) = +1.584919424220149128068025478_wp

    case(10)

      Ext_S2F(1,1) = +1.588005378675122816844966727_wp
      Ext_S2F(1,2) = -1.01628796564473368992920028_wp
      Ext_S2F(1,3) = +0.7575227986514953950529923881_wp
      Ext_S2F(1,4) = -0.5836053892999149682347822469_wp
      Ext_S2F(1,5) = +0.4466023128802576369076345604_wp
      Ext_S2F(1,6) = -0.3308583679390711035409576795_wp
      Ext_S2F(1,7) = +0.2306924543937171515976488509_wp
      Ext_S2F(1,8) = -0.1446071081332395186303446385_wp
      Ext_S2F(1,9) = +0.07352805218733873737532497482_wp
      Ext_S2F(1,10) = -0.02099216577097245744328265545_wp
      Ext_S2F(2,1) = +0.422906646500601630347043083_wp
      Ext_S2F(2,2) = +0.8101708905783888422619855047_wp
      Ext_S2F(2,3) = -0.3884998090533988889790833141_wp
      Ext_S2F(2,4) = +0.2690248358902316298274456953_wp
      Ext_S2F(2,5) = -0.197180403931379974851729302_wp
      Ext_S2F(2,6) = +0.1429634694958448827223768335_wp
      Ext_S2F(2,7) = -0.09849039747561035080667709756_wp
      Ext_S2F(2,8) = +0.06130421194300179946328507552_wp
      Ext_S2F(2,9) = -0.03104480314892165820248224784_wp
      Ext_S2F(2,10) = +0.008845359201242088217835769435_wp
      Ext_S2F(3,1) = -0.06817916320872028949454869217_wp
      Ext_S2F(3,2) = +0.530420913183021595996662495_wp
      Ext_S2F(3,3) = +0.7203641017248800962345340663_wp
      Ext_S2F(3,4) = -0.2935512826577693195324813746_wp
      Ext_S2F(3,5) = +0.1863911099443084041345871251_wp
      Ext_S2F(3,6) = -0.1269308409904718762982907777_wp
      Ext_S2F(3,7) = +0.0846244719393531765242670906_wp
      Ext_S2F(3,8) = -0.0517056009439694137553909005_wp
      Ext_S2F(3,9) = +0.02591091043757951955135991139_wp
      Ext_S2F(3,10) = -0.007344619428211893360698943356_wp
      Ext_S2F(4,1) = +0.02745114698167567932200579559_wp
      Ext_S2F(4,2) = -0.1238297816101166321437315128_wp
      Ext_S2F(4,3) = +0.5758758007211883548452886259_wp
      Ext_S2F(4,4) = +0.6790521368343947694108685836_wp
      Ext_S2F(4,5) = -0.2471697515625376863451722647_wp
      Ext_S2F(4,6) = +0.1441120977223812084887518879_wp
      Ext_S2F(4,7) = -0.08964894253655170910095353426_wp
      Ext_S2F(4,8) = +0.05282645908124611948477050257_wp
      Ext_S2F(4,9) = -0.02595796307419274339445915506_wp
      Ext_S2F(4,10) = +0.007288797442512639432631071275_wp
      Ext_S2F(5,1) = -0.01538254787452211419504869678_wp
      Ext_S2F(5,2) = +0.06064160125376463087228028204_wp
      Ext_S2F(5,3) = -0.1593597433348712569245758194_wp
      Ext_S2F(5,4) = +0.6048289384800086407622696097_wp
      Ext_S2F(5,5) = +0.6514932860732885390492994445_wp
      Ext_S2F(5,6) = -0.2152200138196496137222418718_wp
      Ext_S2F(5,7) = +0.114169367062990596811981281_wp
      Ext_S2F(5,8) = -0.06269546622947128622158296326_wp
      Ext_S2F(5,9) = +0.02974064552429903841230431621_wp
      Ext_S2F(5,10) = -0.008216067135837174844685582238_wp
      Ext_S2F(6,1) = +0.01047049316709865825992295152_wp
      Ext_S2F(6,2) = -0.03901209887058674262121418876_wp
      Ext_S2F(6,3) = +0.0879660699823633548120794374_wp
      Ext_S2F(6,4) = -0.1877654519033348121949107574_wp
      Ext_S2F(6,5) = +0.6283409876244595417441225573_wp
      Ext_S2F(6,6) = +0.6283409876244595417441225573_wp
      Ext_S2F(6,7) = -0.1877654519033348121949107574_wp
      Ext_S2F(6,8) = +0.0879660699823633548120794374_wp
      Ext_S2F(6,9) = -0.03901209887058674262121418876_wp
      Ext_S2F(6,10) = +0.01047049316709865825992295152_wp
      Ext_S2F(7,1) = -0.008216067135837174844685582238_wp
      Ext_S2F(7,2) = +0.02974064552429903841230431621_wp
      Ext_S2F(7,3) = -0.06269546622947128622158296326_wp
      Ext_S2F(7,4) = +0.114169367062990596811981281_wp
      Ext_S2F(7,5) = -0.2152200138196496137222418718_wp
      Ext_S2F(7,6) = +0.6514932860732885390492994445_wp
      Ext_S2F(7,7) = +0.6048289384800086407622696097_wp
      Ext_S2F(7,8) = -0.1593597433348712569245758194_wp
      Ext_S2F(7,9) = +0.06064160125376463087228028204_wp
      Ext_S2F(7,10) = -0.01538254787452211419504869678_wp
      Ext_S2F(8,1) = +0.007288797442512639432631071275_wp
      Ext_S2F(8,2) = -0.02595796307419274339445915506_wp
      Ext_S2F(8,3) = +0.05282645908124611948477050257_wp
      Ext_S2F(8,4) = -0.08964894253655170910095353426_wp
      Ext_S2F(8,5) = +0.1441120977223812084887518879_wp
      Ext_S2F(8,6) = -0.2471697515625376863451722647_wp
      Ext_S2F(8,7) = +0.6790521368343947694108685836_wp
      Ext_S2F(8,8) = +0.5758758007211883548452886259_wp
      Ext_S2F(8,9) = -0.1238297816101166321437315128_wp
      Ext_S2F(8,10) = +0.02745114698167567932200579559_wp
      Ext_S2F(9,1) = -0.007344619428211893360698943356_wp
      Ext_S2F(9,2) = +0.02591091043757951955135991139_wp
      Ext_S2F(9,3) = -0.0517056009439694137553909005_wp
      Ext_S2F(9,4) = +0.0846244719393531765242670906_wp
      Ext_S2F(9,5) = -0.1269308409904718762982907777_wp
      Ext_S2F(9,6) = +0.1863911099443084041345871251_wp
      Ext_S2F(9,7) = -0.2935512826577693195324813746_wp
      Ext_S2F(9,8) = +0.7203641017248800962345340663_wp
      Ext_S2F(9,9) = +0.530420913183021595996662495_wp
      Ext_S2F(9,10) = -0.06817916320872028949454869217_wp
      Ext_S2F(10,1) = +0.008845359201242088217835769435_wp
      Ext_S2F(10,2) = -0.03104480314892165820248224784_wp
      Ext_S2F(10,3) = +0.06130421194300179946328507552_wp
      Ext_S2F(10,4) = -0.09849039747561035080667709756_wp
      Ext_S2F(10,5) = +0.1429634694958448827223768335_wp
      Ext_S2F(10,6) = -0.197180403931379974851729302_wp
      Ext_S2F(10,7) = +0.2690248358902316298274456953_wp
      Ext_S2F(10,8) = -0.3884998090533988889790833141_wp
      Ext_S2F(10,9) = +0.8101708905783888422619855047_wp
      Ext_S2F(10,10) = +0.422906646500601630347043083_wp
      Ext_S2F(11,1) = -0.02099216577097245744328265545_wp
      Ext_S2F(11,2) = +0.07352805218733873737532497482_wp
      Ext_S2F(11,3) = -0.1446071081332395186303446385_wp
      Ext_S2F(11,4) = +0.2306924543937171515976488509_wp
      Ext_S2F(11,5) = -0.3308583679390711035409576795_wp
      Ext_S2F(11,6) = +0.4466023128802576369076345604_wp
      Ext_S2F(11,7) = -0.5836053892999149682347822469_wp
      Ext_S2F(11,8) = +0.7575227986514953950529923881_wp
      Ext_S2F(11,9) = -1.01628796564473368992920028_wp
      Ext_S2F(11,10) = +1.588005378675122816844966727_wp

    case default
      write(*,*)'Ext_S2F only defined for p <= 10'
      write(*,*)'stopping'
      stop
    end select

    return
  end subroutine Extrapolate_GL_2_GLL


  subroutine Filter_GLL_2_GLL(N_Soln_Pts,N_Flux_Pts,Filter)
    ! Filter by Restricting from Gauss_Lobatto_Legendre N 
    ! down to Gauss_Lobatto_Legendre N-1
    !  Followed by Prolongation back to GLL N
    !

    implicit none

    integer,                       intent( in) :: N_Soln_Pts,N_Flux_Pts
    real(wp), dimension(N_Soln_Pts,N_Soln_Pts),   intent(out) :: Filter

    continue 

    select case(N_Soln_Pts)

    case(1)

      Filter(1,1) = 1.0_wp

    case(2)

      Filter(1,1) = 1.0_wp
      Filter(1,2) = 0.0_wp
      Filter(2,1) = 0.0_wp
      Filter(2,2) = 1.0_wp

    case(3)

      Filter(1,1) = +0.6666666666666666666666666667_wp
      Filter(1,2) = +0.6666666666666666666666666667_wp
      Filter(1,3) = -0.3333333333333333333333333333_wp
      Filter(2,1) = +0.1666666666666666666666666667_wp
      Filter(2,2) = +0.6666666666666666666666666667_wp
      Filter(2,3) = +0.1666666666666666666666666667_wp
      Filter(3,1) = -0.3333333333333333333333333333_wp
      Filter(3,2) = +0.6666666666666666666666666667_wp
      Filter(3,3) = +0.6666666666666666666666666667_wp

    case(4)

      Filter(1,1) = +0.7500000000000000000000000000_wp
      Filter(1,2) = +0.5590169943749474241022934172_wp
      Filter(1,3) = -0.5590169943749474241022934172_wp
      Filter(1,4) = +0.2500000000000000000000000000_wp
      Filter(2,1) = +0.1118033988749894848204586834_wp
      Filter(2,2) = +0.7500000000000000000000000000_wp
      Filter(2,3) = +0.2500000000000000000000000000_wp
      Filter(2,4) = -0.1118033988749894848204586834_wp
      Filter(3,1) = -0.1118033988749894848204586834_wp
      Filter(3,2) = +0.2500000000000000000000000000_wp
      Filter(3,3) = +0.7500000000000000000000000000_wp
      Filter(3,4) = +0.1118033988749894848204586834_wp
      Filter(4,1) = +0.2500000000000000000000000000_wp
      Filter(4,2) = -0.5590169943749474241022934172_wp
      Filter(4,3) = +0.5590169943749474241022934172_wp
      Filter(4,4) = +0.7500000000000000000000000000_wp

    case(5)

      Filter(1,1) = +0.800000000000000000000000000_wp
      Filter(1,2) = +0.4666666666666666666666666667_wp
      Filter(1,3) = -0.5333333333333333333333333333_wp
      Filter(1,4) = +0.4666666666666666666666666667_wp
      Filter(1,5) = -0.200000000000000000000000000_wp
      Filter(2,1) = +0.08571428571428571428571428571_wp
      Filter(2,2) = +0.800000000000000000000000000_wp
      Filter(2,3) = +0.2285714285714285714285714286_wp
      Filter(2,4) = -0.200000000000000000000000000_wp
      Filter(2,5) = +0.08571428571428571428571428571_wp
      Filter(3,1) = -0.07500000000000000000000000000_wp
      Filter(3,2) = +0.17500000000000000000000000000_wp
      Filter(3,3) = +0.800000000000000000000000000_wp
      Filter(3,4) = +0.17500000000000000000000000000_wp
      Filter(3,5) = -0.07500000000000000000000000000_wp
      Filter(4,1) = +0.08571428571428571428571428571_wp
      Filter(4,2) = -0.200000000000000000000000000_wp
      Filter(4,3) = +0.2285714285714285714285714286_wp
      Filter(4,4) = +0.800000000000000000000000000_wp
      Filter(4,5) = +0.08571428571428571428571428571_wp
      Filter(5,1) = -0.200000000000000000000000000_wp
      Filter(5,2) = +0.4666666666666666666666666667_wp
      Filter(5,3) = -0.5333333333333333333333333333_wp
      Filter(5,4) = +0.4666666666666666666666666667_wp
      Filter(5,5) = +0.800000000000000000000000000_wp

    case(6)

      Filter(1,1) = +0.8333333333333333333333333333_wp
      Filter(1,2) = +0.3971119470091982032443159701_wp
      Filter(1,3) = -0.4808232423993796985907724435_wp
      Filter(1,4) = +0.4808232423993796985907724435_wp
      Filter(1,5) = -0.3971119470091982032443159701_wp
      Filter(1,6) = +0.1666666666666666666666666667_wp
      Filter(2,1) = +0.06994948902188119805148663988_wp
      Filter(2,2) = +0.8333333333333333333333333333_wp
      Filter(2,3) = +0.2018000406940843945685840876_wp
      Filter(2,4) = -0.2018000406940843945685840876_wp
      Filter(2,5) = +0.1666666666666666666666666667_wp
      Filter(2,6) = -0.06994948902188119805148663988_wp
      Filter(3,1) = -0.05777128750923629100512185927_wp
      Filter(3,2) = +0.1376500107841259762897897786_wp
      Filter(3,3) = +0.8333333333333333333333333333_wp
      Filter(3,4) = +0.1666666666666666666666666667_wp
      Filter(3,5) = -0.1376500107841259762897897786_wp
      Filter(3,6) = +0.05777128750923629100512185927_wp
      Filter(4,1) = +0.05777128750923629100512185927_wp
      Filter(4,2) = -0.1376500107841259762897897786_wp
      Filter(4,3) = +0.1666666666666666666666666667_wp
      Filter(4,4) = +0.8333333333333333333333333333_wp
      Filter(4,5) = +0.1376500107841259762897897786_wp
      Filter(4,6) = -0.05777128750923629100512185927_wp
      Filter(5,1) = -0.06994948902188119805148663988_wp
      Filter(5,2) = +0.1666666666666666666666666667_wp
      Filter(5,3) = -0.2018000406940843945685840876_wp
      Filter(5,4) = +0.2018000406940843945685840876_wp
      Filter(5,5) = +0.8333333333333333333333333333_wp
      Filter(5,6) = +0.06994948902188119805148663988_wp
      Filter(6,1) = +0.1666666666666666666666666667_wp
      Filter(6,2) = -0.3971119470091982032443159701_wp
      Filter(6,3) = +0.4808232423993796985907724435_wp
      Filter(6,4) = -0.4808232423993796985907724435_wp
      Filter(6,5) = +0.3971119470091982032443159701_wp
      Filter(6,6) = +0.8333333333333333333333333333_wp

    case(7)

      Filter(1,1) = +0.8571428571428571428571428571_wp
      Filter(1,2) = +0.3444411917635988313750693971_wp
      Filter(1,3) = -0.4301554774778845456607836828_wp
      Filter(1,4) = +0.4571428571428571428571428571_wp
      Filter(1,5) = -0.4301554774778845456607836828_wp
      Filter(1,6) = +0.3444411917635988313750693971_wp
      Filter(1,7) = -0.1428571428571428571428571429_wp
      Filter(2,1) = +0.05925006576830365642710519047_wp
      Filter(2,2) = +0.8571428571428571428571428571_wp
      Filter(2,3) = +0.1784071823181250509991202571_wp
      Filter(2,4) = -0.1896002104585717005667366095_wp
      Filter(2,5) = +0.1784071823181250509991202571_wp
      Filter(2,6) = -0.1428571428571428571428571429_wp
      Filter(2,7) = +0.05925006576830365642710519047_wp
      Filter(3,1) = -0.04744369032556457732438972412_wp
      Filter(3,2) = +0.1143909286618041107482233084_wp
      Filter(3,3) = +0.8571428571428571428571428571_wp
      Filter(3,4) = +0.1518198090418066474380471172_wp
      Filter(3,5) = -0.1428571428571428571428571429_wp
      Filter(3,6) = +0.1143909286618041107482233084_wp
      Filter(3,7) = -0.04744369032556457732438972412_wp
      Filter(4,1) = +0.04464285714285714285714285714_wp
      Filter(4,2) = -0.1076378724261246348047091866_wp
      Filter(4,3) = +0.1344235867118389205189949009_wp
      Filter(4,4) = +0.8571428571428571428571428571_wp
      Filter(4,5) = +0.1344235867118389205189949009_wp
      Filter(4,6) = -0.1076378724261246348047091866_wp
      Filter(4,7) = +0.04464285714285714285714285714_wp
      Filter(5,1) = -0.04744369032556457732438972412_wp
      Filter(5,2) = +0.1143909286618041107482233084_wp
      Filter(5,3) = -0.1428571428571428571428571429_wp
      Filter(5,4) = +0.1518198090418066474380471172_wp
      Filter(5,5) = +0.8571428571428571428571428571_wp
      Filter(5,6) = +0.1143909286618041107482233084_wp
      Filter(5,7) = -0.04744369032556457732438972412_wp
      Filter(6,1) = +0.05925006576830365642710519047_wp
      Filter(6,2) = -0.1428571428571428571428571429_wp
      Filter(6,3) = +0.1784071823181250509991202571_wp
      Filter(6,4) = -0.1896002104585717005667366095_wp
      Filter(6,5) = +0.1784071823181250509991202571_wp
      Filter(6,6) = +0.8571428571428571428571428571_wp
      Filter(6,7) = +0.05925006576830365642710519047_wp
      Filter(7,1) = -0.1428571428571428571428571429_wp
      Filter(7,2) = +0.3444411917635988313750693971_wp
      Filter(7,3) = -0.4301554774778845456607836828_wp
      Filter(7,4) = +0.4571428571428571428571428571_wp
      Filter(7,5) = -0.4301554774778845456607836828_wp
      Filter(7,6) = +0.3444411917635988313750693971_wp
      Filter(7,7) = +0.8571428571428571428571428571_wp

    case(8)

      Filter(1,1) = +0.87500000000000000000000000000_wp
      Filter(1,2) = +0.3036166981166943199817710034_wp
      Filter(1,3) = -0.3863174574899937907080670461_wp
      Filter(1,4) = +0.4247949183584744204226588286_wp
      Filter(1,5) = -0.4247949183584744204226588286_wp
      Filter(1,6) = +0.3863174574899937907080670461_wp
      Filter(1,7) = -0.3036166981166943199817710034_wp
      Filter(1,8) = +0.12500000000000000000000000000_wp
      Filter(2,1) = +0.05146291392048065240393693195_wp
      Filter(2,2) = +0.87500000000000000000000000000_wp
      Filter(2,3) = +0.1590481764862919530173529542_wp
      Filter(2,4) = -0.1748894745387182039953789983_wp
      Filter(2,5) = +0.1748894745387182039953789983_wp
      Filter(2,6) = -0.1590481764862919530173529542_wp
      Filter(2,7) = +0.12500000000000000000000000000_wp
      Filter(2,8) = -0.05146291392048065240393693195_wp
      Filter(3,1) = -0.0404460106502039485171070478_wp
      Filter(3,2) = +0.09824067364486060467301644646_wp
      Filter(3,3) = +0.87500000000000000000000000000_wp
      Filter(3,4) = +0.1374500783366349857587300542_wp
      Filter(3,5) = -0.1374500783366349857587300542_wp
      Filter(3,6) = +0.12500000000000000000000000000_wp
      Filter(3,7) = -0.09824067364486060467301644646_wp
      Filter(3,8) = +0.0404460106502039485171070478_wp
      Filter(4,1) = +0.03678245507356665395837128874_wp
      Filter(4,2) = -0.08934214046449566511855378844_wp
      Filter(4,3) = +0.1136776361940815381619446243_wp
      Filter(4,4) = +0.87500000000000000000000000000_wp
      Filter(4,5) = +0.12500000000000000000000000000_wp
      Filter(4,6) = -0.1136776361940815381619446243_wp
      Filter(4,7) = +0.08934214046449566511855378844_wp
      Filter(4,8) = -0.03678245507356665395837128874_wp
      Filter(5,1) = -0.03678245507356665395837128874_wp
      Filter(5,2) = +0.08934214046449566511855378844_wp
      Filter(5,3) = -0.1136776361940815381619446243_wp
      Filter(5,4) = +0.12500000000000000000000000000_wp
      Filter(5,5) = +0.87500000000000000000000000000_wp
      Filter(5,6) = +0.1136776361940815381619446243_wp
      Filter(5,7) = -0.08934214046449566511855378844_wp
      Filter(5,8) = +0.03678245507356665395837128874_wp
      Filter(6,1) = +0.0404460106502039485171070478_wp
      Filter(6,2) = -0.09824067364486060467301644646_wp
      Filter(6,3) = +0.12500000000000000000000000000_wp
      Filter(6,4) = -0.1374500783366349857587300542_wp
      Filter(6,5) = +0.1374500783366349857587300542_wp
      Filter(6,6) = +0.87500000000000000000000000000_wp
      Filter(6,7) = +0.09824067364486060467301644646_wp
      Filter(6,8) = -0.0404460106502039485171070478_wp
      Filter(7,1) = -0.05146291392048065240393693195_wp
      Filter(7,2) = +0.12500000000000000000000000000_wp
      Filter(7,3) = -0.1590481764862919530173529542_wp
      Filter(7,4) = +0.1748894745387182039953789983_wp
      Filter(7,5) = -0.1748894745387182039953789983_wp
      Filter(7,6) = +0.1590481764862919530173529542_wp
      Filter(7,7) = +0.87500000000000000000000000000_wp
      Filter(7,8) = +0.05146291392048065240393693195_wp
      Filter(8,1) = +0.12500000000000000000000000000_wp
      Filter(8,2) = -0.3036166981166943199817710034_wp
      Filter(8,3) = +0.3863174574899937907080670461_wp
      Filter(8,4) = -0.4247949183584744204226588286_wp
      Filter(8,5) = +0.4247949183584744204226588286_wp
      Filter(8,6) = -0.3863174574899937907080670461_wp
      Filter(8,7) = +0.3036166981166943199817710034_wp
      Filter(8,8) = +0.8750000000000000000000000000_wp

    case(9)

      Filter(1,1) = +0.8888888888888888888888888889_wp
      Filter(1,2) = +0.2712074741356231093026211754_wp
      Filter(1,3) = -0.349309612744378161211946529_wp
      Filter(1,4) = +0.3923878528944693376236110678_wp
      Filter(1,5) = -0.4063492063492063492063492063_wp
      Filter(1,6) = +0.3923878528944693376236110678_wp
      Filter(1,7) = -0.349309612744378161211946529_wp
      Filter(1,8) = +0.2712074741356231093026211754_wp
      Filter(1,9) = -0.1111111111111111111111111111_wp
      Filter(2,1) = +0.04552116069696500200270664893_wp
      Filter(2,2) = +0.8888888888888888888888888889_wp
      Filter(2,3) = +0.1431088111325830709068651887_wp
      Filter(2,4) = -0.1607575545643158226034100274_wp
      Filter(2,5) = +0.1664773876917577216098986018_wp
      Filter(2,6) = -0.1607575545643158226034100274_wp
      Filter(2,7) = +0.1431088111325830709068651887_wp
      Filter(2,8) = -0.1111111111111111111111111111_wp
      Filter(2,9) = +0.04552116069696500200270664893_wp
      Filter(3,1) = -0.03534308407762068418285371681_wp
      Filter(3,2) = +0.08626777704769018178100951725_wp
      Filter(3,3) = +0.8888888888888888888888888889_wp
      Filter(3,4) = +0.1248137718829765788758878214_wp
      Filter(3,5) = -0.1292547074838699307258650215_wp
      Filter(3,6) = +0.1248137718829765788758878214_wp
      Filter(3,7) = -0.1111111111111111111111111111_wp
      Filter(3,8) = +0.08626777704769018178100951725_wp
      Filter(3,9) = -0.03534308407762068418285371681_wp
      Filter(4,1) = +0.03146294902168132320578809352_wp
      Filter(4,2) = -0.07679688239725259286392689588_wp
      Filter(4,3) = +0.0989127948470365326215231139_wp
      Filter(4,4) = +0.8888888888888888888888888889_wp
      Filter(4,5) = +0.1150644992792916962954535992_wp
      Filter(4,6) = -0.1111111111111111111111111111_wp
      Filter(4,7) = +0.0989127948470365326215231139_wp
      Filter(4,8) = -0.07679688239725259286392689588_wp
      Filter(4,9) = +0.03146294902168132320578809352_wp
      Filter(5,1) = -0.03038194444444444444444444444_wp
      Filter(5,2) = +0.07415829370895944394993547765_wp
      Filter(5,3) = -0.09551434723479090345639162901_wp
      Filter(5,4) = +0.1072935535258314595064561514_wp
      Filter(5,5) = +0.8888888888888888888888888889_wp
      Filter(5,6) = +0.1072935535258314595064561514_wp
      Filter(5,7) = -0.09551434723479090345639162901_wp
      Filter(5,8) = +0.07415829370895944394993547765_wp
      Filter(5,9) = -0.03038194444444444444444444444_wp
      Filter(6,1) = +0.03146294902168132320578809352_wp
      Filter(6,2) = -0.07679688239725259286392689588_wp
      Filter(6,3) = +0.0989127948470365326215231139_wp
      Filter(6,4) = -0.1111111111111111111111111111_wp
      Filter(6,5) = +0.1150644992792916962954535992_wp
      Filter(6,6) = +0.8888888888888888888888888889_wp
      Filter(6,7) = +0.0989127948470365326215231139_wp
      Filter(6,8) = -0.07679688239725259286392689588_wp
      Filter(6,9) = +0.03146294902168132320578809352_wp
      Filter(7,1) = -0.03534308407762068418285371681_wp
      Filter(7,2) = +0.08626777704769018178100951725_wp
      Filter(7,3) = -0.1111111111111111111111111111_wp
      Filter(7,4) = +0.1248137718829765788758878214_wp
      Filter(7,5) = -0.1292547074838699307258650215_wp
      Filter(7,6) = +0.1248137718829765788758878214_wp
      Filter(7,7) = +0.8888888888888888888888888889_wp
      Filter(7,8) = +0.08626777704769018178100951725_wp
      Filter(7,9) = -0.03534308407762068418285371681_wp
      Filter(8,1) = +0.04552116069696500200270664893_wp
      Filter(8,2) = -0.1111111111111111111111111111_wp
      Filter(8,3) = +0.1431088111325830709068651887_wp
      Filter(8,4) = -0.1607575545643158226034100274_wp
      Filter(8,5) = +0.1664773876917577216098986018_wp
      Filter(8,6) = -0.1607575545643158226034100274_wp
      Filter(8,7) = +0.1431088111325830709068651887_wp
      Filter(8,8) = +0.8888888888888888888888888889_wp
      Filter(8,9) = +0.04552116069696500200270664893_wp
      Filter(9,1) = -0.1111111111111111111111111111_wp
      Filter(9,2) = +0.2712074741356231093026211754_wp
      Filter(9,3) = -0.349309612744378161211946529_wp
      Filter(9,4) = +0.3923878528944693376236110678_wp
      Filter(9,5) = -0.4063492063492063492063492063_wp
      Filter(9,6) = +0.3923878528944693376236110678_wp
      Filter(9,7) = -0.349309612744378161211946529_wp
      Filter(9,8) = +0.2712074741356231093026211754_wp
      Filter(9,9) = +0.8888888888888888888888888889_wp

    case(10)

      Filter(1,1) = +0.900000000000000000000000000_wp
      Filter(1,2) = +0.2449238573168844040096643206_wp
      Filter(1,3) = -0.3181197949333032363630618483_wp
      Filter(1,4) = +0.3625178721881966622166893225_wp
      Filter(1,5) = -0.3839178200250072469868738773_wp
      Filter(1,6) = +0.3839178200250072469868738773_wp
      Filter(1,7) = -0.3625178721881966622166893225_wp
      Filter(1,8) = +0.3181197949333032363630618483_wp
      Filter(1,9) = -0.2449238573168844040096643206_wp
      Filter(1,10) = +0.100000000000000000000000000_wp
      Filter(2,1) = +0.04082901563591627464004621494_wp
      Filter(2,2) = +0.900000000000000000000000000_wp
      Filter(2,3) = +0.1298851808142631672025764087_wp
      Filter(2,4) = -0.1480124787187097911694958176_wp
      Filter(2,5) = +0.1567498667670791114056612122_wp
      Filter(2,6) = -0.1567498667670791114056612122_wp
      Filter(2,7) = +0.1480124787187097911694958176_wp
      Filter(2,8) = -0.1298851808142631672025764087_wp
      Filter(2,9) = +0.100000000000000000000000000_wp
      Filter(2,10) = -0.04082901563591627464004621494_wp
      Filter(3,1) = -0.03143469900103699202141903836_wp
      Filter(3,2) = +0.07699107732929192943730975095_wp
      Filter(3,3) = +0.900000000000000000000000000_wp
      Filter(3,4) = +0.1139564019473236157066726827_wp
      Filter(3,5) = -0.1206834111362039499809756011_wp
      Filter(3,6) = +0.1206834111362039499809756011_wp
      Filter(3,7) = -0.1139564019473236157066726827_wp
      Filter(3,8) = +0.100000000000000000000000000_wp
      Filter(3,9) = -0.07699107732929192943730975095_wp
      Filter(3,10) = +0.03143469900103699202141903836_wp
      Filter(4,1) = +0.02758484689220680377139658354_wp
      Filter(4,2) = -0.06756187104334961388467391578_wp
      Filter(4,3) = +0.08775285836615395500231262886_wp
      Filter(4,4) = +0.900000000000000000000000000_wp
      Filter(4,5) = +0.1059031428457963217241087054_wp
      Filter(4,6) = -0.1059031428457963217241087054_wp
      Filter(4,7) = +0.100000000000000000000000000_wp
      Filter(4,8) = -0.08775285836615395500231262886_wp
      Filter(4,9) = +0.06756187104334961388467391578_wp
      Filter(4,10) = -0.02758484689220680377139658354_wp
      Filter(5,1) = -0.02604724104587963697501472807_wp
      Filter(5,2) = +0.06379590749419519100600054136_wp
      Filter(5,3) = -0.08286142980093549029734717573_wp
      Filter(5,4) = +0.09442590400325344189102951936_wp
      Filter(5,5) = +0.900000000000000000000000000_wp
      Filter(5,6) = +0.100000000000000000000000000_wp
      Filter(5,7) = -0.09442590400325344189102951936_wp
      Filter(5,8) = +0.08286142980093549029734717573_wp
      Filter(5,9) = -0.06379590749419519100600054136_wp
      Filter(5,10) = +0.02604724104587963697501472807_wp
      Filter(6,1) = +0.02604724104587963697501472807_wp
      Filter(6,2) = -0.06379590749419519100600054136_wp
      Filter(6,3) = +0.08286142980093549029734717573_wp
      Filter(6,4) = -0.09442590400325344189102951936_wp
      Filter(6,5) = +0.100000000000000000000000000_wp
      Filter(6,6) = +0.900000000000000000000000000_wp
      Filter(6,7) = +0.09442590400325344189102951936_wp
      Filter(6,8) = -0.08286142980093549029734717573_wp
      Filter(6,9) = +0.06379590749419519100600054136_wp
      Filter(6,10) = -0.02604724104587963697501472807_wp
      Filter(7,1) = -0.02758484689220680377139658354_wp
      Filter(7,2) = +0.06756187104334961388467391578_wp
      Filter(7,3) = -0.08775285836615395500231262886_wp
      Filter(7,4) = +0.100000000000000000000000000_wp
      Filter(7,5) = -0.1059031428457963217241087054_wp
      Filter(7,6) = +0.1059031428457963217241087054_wp
      Filter(7,7) = +0.900000000000000000000000000_wp
      Filter(7,8) = +0.08775285836615395500231262886_wp
      Filter(7,9) = -0.06756187104334961388467391578_wp
      Filter(7,10) = +0.02758484689220680377139658354_wp
      Filter(8,1) = +0.03143469900103699202141903836_wp
      Filter(8,2) = -0.07699107732929192943730975095_wp
      Filter(8,3) = +0.100000000000000000000000000_wp
      Filter(8,4) = -0.1139564019473236157066726827_wp
      Filter(8,5) = +0.1206834111362039499809756011_wp
      Filter(8,6) = -0.1206834111362039499809756011_wp
      Filter(8,7) = +0.1139564019473236157066726827_wp
      Filter(8,8) = +0.900000000000000000000000000_wp
      Filter(8,9) = +0.07699107732929192943730975095_wp
      Filter(8,10) = -0.03143469900103699202141903836_wp
      Filter(9,1) = -0.04082901563591627464004621494_wp
      Filter(9,2) = +0.100000000000000000000000000_wp
      Filter(9,3) = -0.1298851808142631672025764087_wp
      Filter(9,4) = +0.1480124787187097911694958176_wp
      Filter(9,5) = -0.1567498667670791114056612122_wp
      Filter(9,6) = +0.1567498667670791114056612122_wp
      Filter(9,7) = -0.1480124787187097911694958176_wp
      Filter(9,8) = +0.1298851808142631672025764087_wp
      Filter(9,9) = +0.900000000000000000000000000_wp
      Filter(9,10) = +0.04082901563591627464004621494_wp
      Filter(10,1) = +0.100000000000000000000000000_wp
      Filter(10,2) = -0.2449238573168844040096643206_wp
      Filter(10,3) = +0.3181197949333032363630618483_wp
      Filter(10,4) = -0.3625178721881966622166893225_wp
      Filter(10,5) = +0.3839178200250072469868738773_wp
      Filter(10,6) = -0.3839178200250072469868738773_wp
      Filter(10,7) = +0.3625178721881966622166893225_wp
      Filter(10,8) = -0.3181197949333032363630618483_wp
      Filter(10,9) = +0.2449238573168844040096643206_wp
      Filter(10,10) = +0.900000000000000000000000000_wp

    case(11)

      Filter(1,1) = +0.9090909090909090909090909091_wp
      Filter(1,2) = +0.2232123665389235475487200701_wp
      Filter(1,3) = -0.2916799941563490570160773687_wp
      Filter(1,4) = +0.3357813846862732976367611448_wp
      Filter(1,5) = -0.3611088508639415832631989401_wp
      Filter(1,6) = +0.3694083694083694083694083694_wp
      Filter(1,7) = -0.3611088508639415832631989401_wp
      Filter(1,8) = +0.3357813846862732976367611448_wp
      Filter(1,9) = -0.2916799941563490570160773687_wp
      Filter(1,10) = +0.2232123665389235475487200701_wp
      Filter(1,11) = -0.09090909090909090909090909091_wp
      Filter(2,1) = +0.03702511172684604226488616819_wp
      Filter(2,2) = +0.9090909090909090909090909091_wp
      Filter(2,3) = +0.1187943280933708700626833367_wp
      Filter(2,4) = -0.1367557761218477350145593368_wp
      Filter(2,5) = +0.147071051036694643046239074_wp
      Filter(2,6) = -0.1504512476519458225366803025_wp
      Filter(2,7) = +0.147071051036694643046239074_wp
      Filter(2,8) = -0.1367557761218477350145593368_wp
      Filter(2,9) = +0.1187943280933708700626833367_wp
      Filter(2,10) = -0.09090909090909090909090909091_wp
      Filter(2,11) = +0.03702511172684604226488616819_wp
      Filter(3,1) = -0.02833400636139398780915135023_wp
      Filter(3,2) = +0.06956950674801232931766906977_wp
      Filter(3,3) = +0.9090909090909090909090909091_wp
      Filter(3,4) = +0.1046543507871240431440386691_wp
      Filter(3,5) = -0.1125482652528805301580182951_wp
      Filter(3,6) = +0.1151350099764581091927419946_wp
      Filter(3,7) = -0.1125482652528805301580182951_wp
      Filter(3,8) = +0.1046543507871240431440386691_wp
      Filter(3,9) = -0.09090909090909090909090909091_wp
      Filter(3,10) = +0.06956950674801232931766906977_wp
      Filter(3,11) = -0.02833400636139398780915135023_wp
      Filter(4,1) = +0.02461262948700683469200562117_wp
      Filter(4,2) = -0.06043227601994536170154754883_wp
      Filter(4,3) = +0.07896912787436791835093091116_wp
      Filter(4,4) = +0.9090909090909090909090909091_wp
      Filter(4,5) = +0.09776622185872302315549978217_wp
      Filter(4,6) = -0.1000132245821230108119593495_wp
      Filter(4,7) = +0.09776622185872302315549978217_wp
      Filter(4,8) = -0.09090909090909090909090909091_wp
      Filter(4,9) = +0.07896912787436791835093091116_wp
      Filter(4,10) = -0.06043227601994536170154754883_wp
      Filter(4,11) = +0.02461262948700683469200562117_wp
      Filter(5,1) = -0.02288634795337995090938723583_wp
      Filter(5,2) = +0.05619367476917907205809897631_wp
      Filter(5,3) = -0.07343038821032238810225466134_wp
      Filter(5,4) = +0.08453290566817554562537953262_wp
      Filter(5,5) = +0.9090909090909090909090909091_wp
      Filter(5,6) = +0.09299849327087726083814495828_wp
      Filter(5,7) = -0.09090909090909090909090909091_wp
      Filter(5,8) = +0.08453290566817554562537953262_wp
      Filter(5,9) = -0.07343038821032238810225466134_wp
      Filter(5,10) = +0.05619367476917907205809897631_wp
      Filter(5,11) = -0.02288634795337995090938723583_wp
      Filter(6,1) = +0.02237215909090909090909090909_wp
      Filter(6,2) = -0.05493116832793821677956782976_wp
      Filter(6,3) = +0.07178062356191402575005028995_wp
      Filter(6,4) = -0.08263370013763756934029668798_wp
      Filter(6,5) = +0.08886663126729812400617786416_wp
      Filter(6,6) = +0.9090909090909090909090909091_wp
      Filter(6,7) = +0.08886663126729812400617786416_wp
      Filter(6,8) = -0.08263370013763756934029668798_wp
      Filter(6,9) = +0.07178062356191402575005028995_wp
      Filter(6,10) = -0.05493116832793821677956782976_wp
      Filter(6,11) = +0.02237215909090909090909090909_wp
      Filter(7,1) = -0.02288634795337995090938723583_wp
      Filter(7,2) = +0.05619367476917907205809897631_wp
      Filter(7,3) = -0.07343038821032238810225466134_wp
      Filter(7,4) = +0.08453290566817554562537953262_wp
      Filter(7,5) = -0.09090909090909090909090909091_wp
      Filter(7,6) = +0.09299849327087726083814495828_wp
      Filter(7,7) = +0.9090909090909090909090909091_wp
      Filter(7,8) = +0.08453290566817554562537953262_wp
      Filter(7,9) = -0.07343038821032238810225466134_wp
      Filter(7,10) = +0.05619367476917907205809897631_wp
      Filter(7,11) = -0.02288634795337995090938723583_wp
      Filter(8,1) = +0.02461262948700683469200562117_wp
      Filter(8,2) = -0.06043227601994536170154754883_wp
      Filter(8,3) = +0.07896912787436791835093091116_wp
      Filter(8,4) = -0.09090909090909090909090909091_wp
      Filter(8,5) = +0.09776622185872302315549978217_wp
      Filter(8,6) = -0.1000132245821230108119593495_wp
      Filter(8,7) = +0.09776622185872302315549978217_wp
      Filter(8,8) = +0.9090909090909090909090909091_wp
      Filter(8,9) = +0.07896912787436791835093091116_wp
      Filter(8,10) = -0.06043227601994536170154754883_wp
      Filter(8,11) = +0.02461262948700683469200562117_wp
      Filter(9,1) = -0.02833400636139398780915135023_wp
      Filter(9,2) = +0.06956950674801232931766906977_wp
      Filter(9,3) = -0.09090909090909090909090909091_wp
      Filter(9,4) = +0.1046543507871240431440386691_wp
      Filter(9,5) = -0.1125482652528805301580182951_wp
      Filter(9,6) = +0.1151350099764581091927419946_wp
      Filter(9,7) = -0.1125482652528805301580182951_wp
      Filter(9,8) = +0.1046543507871240431440386691_wp
      Filter(9,9) = +0.9090909090909090909090909091_wp
      Filter(9,10) = +0.06956950674801232931766906977_wp
      Filter(9,11) = -0.02833400636139398780915135023_wp
      Filter(10,1) = +0.03702511172684604226488616819_wp
      Filter(10,2) = -0.09090909090909090909090909091_wp
      Filter(10,3) = +0.1187943280933708700626833367_wp
      Filter(10,4) = -0.1367557761218477350145593368_wp
      Filter(10,5) = +0.147071051036694643046239074_wp
      Filter(10,6) = -0.1504512476519458225366803025_wp
      Filter(10,7) = +0.147071051036694643046239074_wp
      Filter(10,8) = -0.1367557761218477350145593368_wp
      Filter(10,9) = +0.1187943280933708700626833367_wp
      Filter(10,10) = +0.9090909090909090909090909091_wp
      Filter(10,11) = +0.03702511172684604226488616819_wp
      Filter(11,1) = -0.09090909090909090909090909091_wp
      Filter(11,2) = +0.2232123665389235475487200701_wp
      Filter(11,3) = -0.2916799941563490570160773687_wp
      Filter(11,4) = +0.3357813846862732976367611448_wp
      Filter(11,5) = -0.3611088508639415832631989401_wp
      Filter(11,6) = +0.3694083694083694083694083694_wp
      Filter(11,7) = -0.3611088508639415832631989401_wp
      Filter(11,8) = +0.3357813846862732976367611448_wp
      Filter(11,9) = -0.2916799941563490570160773687_wp
      Filter(11,10) = +0.2232123665389235475487200701_wp
      Filter(11,11) = +0.9090909090909090909090909091_wp

    case(12)

      Filter(1,1) = 0.9166666666666666666666666667_wp
      Filter(1,2) = 0.2049928547073000458251347975_wp
      Filter(1,3) = -0.269081908317281034027143975_wp
      Filter(1,4) = 0.312089236395299859257177019_wp
      Filter(1,5) = -0.3393640888480010050670093064_wp
      Filter(1,6) = 0.3526954334134987744806119808_wp
      Filter(1,7) = -0.3526954334134987744806119808_wp
      Filter(1,8) = 0.3393640888480010050670093064_wp
      Filter(1,9) = -0.312089236395299859257177019_wp
      Filter(1,10) = 0.269081908317281034027143975_wp
      Filter(1,11) = -0.2049928547073000458251347975_wp
      Filter(1,12) = 0.08333333333333333333333333333_wp
      Filter(2,1) = 0.03387651952240042784054583799_wp
      Filter(2,2) = 0.9166666666666666666666666667_wp
      Filter(2,3) = 0.1093867022428215960798999739_wp
      Filter(2,4) = -0.1268699653137170137218837501_wp
      Filter(2,5) = 0.1379576901727312724170124716_wp
      Filter(2,6) = -0.1433771248259264568334520746_wp
      Filter(2,7) = 0.1433771248259264568334520746_wp
      Filter(2,8) = -0.1379576901727312724170124716_wp
      Filter(2,9) = 0.1268699653137170137218837501_wp
      Filter(2,10) = -0.1093867022428215960798999739_wp
      Filter(2,11) = 0.08333333333333333333333333333_wp
      Filter(2,12) = -0.03387651952240042784054583799_wp
      Filter(3,1) = -0.02580792030156141454547333401_wp
      Filter(3,2) = 0.06348527108010669883079408028_wp
      Filter(3,3) = 0.9166666666666666666666666667_wp
      Filter(3,4) = 0.09665248967838070482604474474_wp
      Filter(3,5) = -0.1050993762984146013315341516_wp
      Filter(3,6) = 0.1092280276351228451495659995_wp
      Filter(3,7) = -0.1092280276351228451495659995_wp
      Filter(3,8) = 0.1050993762984146013315341516_wp
      Filter(3,9) = -0.09665248967838070482604474474_wp
      Filter(3,10) = 0.08333333333333333333333333333_wp
      Filter(3,11) = -0.06348527108010669883079408028_wp
      Filter(3,12) = 0.02580792030156141454547333401_wp
      Filter(4,1) = 0.02225147052379737082031922677_wp
      Filter(4,2) = -0.05473670956930276899974491276_wp
      Filter(4,3) = 0.07184961781690950626555960986_wp
      Filter(4,4) = 0.9166666666666666666666666667_wp
      Filter(4,5) = 0.09061620023803975696361458486_wp
      Filter(4,6) = -0.09417590448574087557258874637_wp
      Filter(4,7) = 0.09417590448574087557258874637_wp
      Filter(4,8) = -0.09061620023803975696361458486_wp
      Filter(4,9) = 0.08333333333333333333333333333_wp
      Filter(4,10) = -0.07184961781690950626555960986_wp
      Filter(4,11) = 0.05473670956930276899974491276_wp
      Filter(4,12) = -0.02225147052379737082031922677_wp
      Filter(5,1) = -0.02046310930545988468472949172_wp
      Filter(5,2) = 0.05033749431256485232212128954_wp
      Filter(5,3) = -0.06607503002421908481928529138_wp
      Filter(5,4) = 0.07663579388897436321981419679_wp
      Filter(5,5) = 0.9166666666666666666666666667_wp
      Filter(5,6) = 0.08660694246572368695457447066_wp
      Filter(5,7) = -0.08660694246572368695457447066_wp
      Filter(5,8) = 0.08333333333333333333333333333_wp
      Filter(5,9) = -0.07663579388897436321981419679_wp
      Filter(5,10) = 0.06607503002421908481928529138_wp
      Filter(5,11) = -0.05033749431256485232212128954_wp
      Filter(5,12) = 0.02046310930545988468472949172_wp
      Filter(6,1) = 0.01968963526755619881060302716_wp
      Filter(6,2) = -0.04843481449970254471919531599_wp
      Filter(6,3) = 0.06357749558238312372590334669_wp
      Filter(6,4) = -0.07373907882664295668746386174_wp
      Filter(6,5) = 0.08018346158788411073276964394_wp
      Filter(6,6) = 0.9166666666666666666666666667_wp
      Filter(6,7) = 0.08333333333333333333333333333_wp
      Filter(6,8) = -0.08018346158788411073276964394_wp
      Filter(6,9) = 0.07373907882664295668746386174_wp
      Filter(6,10) = -0.06357749558238312372590334669_wp
      Filter(6,11) = 0.04843481449970254471919531599_wp
      Filter(6,12) = -0.01968963526755619881060302716_wp
      Filter(7,1) = -0.01968963526755619881060302716_wp
      Filter(7,2) = 0.04843481449970254471919531599_wp
      Filter(7,3) = -0.06357749558238312372590334669_wp
      Filter(7,4) = 0.07373907882664295668746386174_wp
      Filter(7,5) = -0.08018346158788411073276964394_wp
      Filter(7,6) = 0.08333333333333333333333333333_wp
      Filter(7,7) = 0.9166666666666666666666666667_wp
      Filter(7,8) = 0.08018346158788411073276964394_wp
      Filter(7,9) = -0.07373907882664295668746386174_wp
      Filter(7,10) = 0.06357749558238312372590334669_wp
      Filter(7,11) = -0.04843481449970254471919531599_wp
      Filter(7,12) = 0.01968963526755619881060302716_wp
      Filter(8,1) = 0.02046310930545988468472949172_wp
      Filter(8,2) = -0.05033749431256485232212128954_wp
      Filter(8,3) = 0.06607503002421908481928529138_wp
      Filter(8,4) = -0.07663579388897436321981419679_wp
      Filter(8,5) = 0.08333333333333333333333333333_wp
      Filter(8,6) = -0.08660694246572368695457447066_wp
      Filter(8,7) = 0.08660694246572368695457447066_wp
      Filter(8,8) = 0.9166666666666666666666666667_wp
      Filter(8,9) = 0.07663579388897436321981419679_wp
      Filter(8,10) = -0.06607503002421908481928529138_wp
      Filter(8,11) = 0.05033749431256485232212128954_wp
      Filter(8,12) = -0.02046310930545988468472949172_wp
      Filter(9,1) = -0.02225147052379737082031922677_wp
      Filter(9,2) = 0.05473670956930276899974491276_wp
      Filter(9,3) = -0.07184961781690950626555960986_wp
      Filter(9,4) = 0.08333333333333333333333333333_wp
      Filter(9,5) = -0.09061620023803975696361458486_wp
      Filter(9,6) = 0.09417590448574087557258874637_wp
      Filter(9,7) = -0.09417590448574087557258874637_wp
      Filter(9,8) = 0.09061620023803975696361458486_wp
      Filter(9,9) = 0.9166666666666666666666666667_wp
      Filter(9,10) = 0.07184961781690950626555960986_wp
      Filter(9,11) = -0.05473670956930276899974491276_wp
      Filter(9,12) = 0.02225147052379737082031922677_wp
      Filter(10,1) = 0.02580792030156141454547333401_wp
      Filter(10,2) = -0.06348527108010669883079408028_wp
      Filter(10,3) = 0.08333333333333333333333333333_wp
      Filter(10,4) = -0.09665248967838070482604474474_wp
      Filter(10,5) = 0.1050993762984146013315341516_wp
      Filter(10,6) = -0.1092280276351228451495659995_wp
      Filter(10,7) = 0.1092280276351228451495659995_wp
      Filter(10,8) = -0.1050993762984146013315341516_wp
      Filter(10,9) = 0.09665248967838070482604474474_wp
      Filter(10,10) = 0.9166666666666666666666666667_wp
      Filter(10,11) = 0.06348527108010669883079408028_wp
      Filter(10,12) = -0.02580792030156141454547333401_wp
      Filter(11,1) = -0.03387651952240042784054583799_wp
      Filter(11,2) = 0.08333333333333333333333333333_wp
      Filter(11,3) = -0.1093867022428215960798999739_wp
      Filter(11,4) = 0.1268699653137170137218837501_wp
      Filter(11,5) = -0.1379576901727312724170124716_wp
      Filter(11,6) = 0.1433771248259264568334520746_wp
      Filter(11,7) = -0.1433771248259264568334520746_wp
      Filter(11,8) = 0.1379576901727312724170124716_wp
      Filter(11,9) = -0.1268699653137170137218837501_wp
      Filter(11,10) = 0.1093867022428215960798999739_wp
      Filter(11,11) = 0.9166666666666666666666666667_wp
      Filter(11,12) = 0.03387651952240042784054583799_wp
      Filter(12,1) = 0.08333333333333333333333333333_wp
      Filter(12,2) = -0.2049928547073000458251347975_wp
      Filter(12,3) = 0.269081908317281034027143975_wp
      Filter(12,4) = -0.312089236395299859257177019_wp
      Filter(12,5) = 0.3393640888480010050670093064_wp
      Filter(12,6) = -0.3526954334134987744806119808_wp
      Filter(12,7) = 0.3526954334134987744806119808_wp
      Filter(12,8) = -0.3393640888480010050670093064_wp
      Filter(12,9) = 0.312089236395299859257177019_wp
      Filter(12,10) = -0.269081908317281034027143975_wp
      Filter(12,11) = 0.2049928547073000458251347975_wp
      Filter(12,12) = 0.9166666666666666666666666667_wp

    case(13)

      Filter(1,1) = 0.9230769230769230769230769231_wp
      Filter(1,2) = 0.1894953054991709300421621982_wp
      Filter(1,3) = -0.2495983789607200915631795105_wp
      Filter(1,4) = 0.2911358645588133904045865895_wp
      Filter(1,5) = -0.319206559769286565729778801_wp
      Filter(1,6) = 0.3355930160912697560936287713_wp
      Filter(1,7) = -0.340992340992340992340992341_wp
      Filter(1,8) = 0.3355930160912697560936287713_wp
      Filter(1,9) = -0.319206559769286565729778801_wp
      Filter(1,10) = 0.2911358645588133904045865895_wp
      Filter(1,11) = -0.2495983789607200915631795105_wp
      Filter(1,12) = 0.1894953054991709300421621982_wp
      Filter(1,13) = -0.07692307692307692307692307692_wp
      Filter(2,1) = 0.03122589104635891845684205556_wp
      Filter(2,2) = 0.9230769230769230769230769231_wp
      Filter(2,3) = 0.101321113228078246852204017_wp
      Filter(2,4) = -0.1181826982232131796997580261_wp
      Filter(2,5) = 0.1295776203363043782181859723_wp
      Filter(2,6) = -0.1362294824330045413994822512_wp
      Filter(2,7) = 0.1384212659371062012978626186_wp
      Filter(2,8) = -0.1362294824330045413994822512_wp
      Filter(2,9) = 0.1295776203363043782181859723_wp
      Filter(2,10) = -0.1181826982232131796997580261_wp
      Filter(2,11) = 0.101321113228078246852204017_wp
      Filter(2,12) = -0.07692307692307692307692307692_wp
      Filter(2,13) = 0.03122589104635891845684205556_wp
      Filter(3,1) = -0.02370672352902102540102067364_wp
      Filter(3,2) = 0.0584000666277108962961689308_wp
      Filter(3,3) = 0.9230769230769230769230769231_wp
      Filter(3,4) = 0.08972440685621789811407082899_wp
      Filter(3,5) = -0.09837544159430522710834530347_wp
      Filter(3,6) = 0.1034255410596984721081816964_wp
      Filter(3,7) = -0.1050895449944481818642648043_wp
      Filter(3,8) = 0.1034255410596984721081816964_wp
      Filter(3,9) = -0.09837544159430522710834530347_wp
      Filter(3,10) = 0.08972440685621789811407082899_wp
      Filter(3,11) = -0.07692307692307692307692307692_wp
      Filter(3,12) = 0.0584000666277108962961689308_wp
      Filter(3,13) = -0.02370672352902102540102067364_wp
      Filter(4,1) = 0.02032439312236731658326326776_wp
      Filter(4,2) = -0.05006790208950716154035143871_wp
      Filter(4,3) = 0.06594816249714277729146337237_wp
      Filter(4,4) = 0.9230769230769230769230769231_wp
      Filter(4,5) = 0.08433983490386245540131491655_wp
      Filter(4,6) = -0.08866941704607879405543714103_wp
      Filter(4,7) = 0.09009601107058065879334019996_wp
      Filter(4,8) = -0.08866941704607879405543714103_wp
      Filter(4,9) = 0.08433983490386245540131491655_wp
      Filter(4,10) = -0.07692307692307692307692307692_wp
      Filter(4,11) = 0.06594816249714277729146337237_wp
      Filter(4,12) = -0.05006790208950716154035143871_wp
      Filter(4,13) = 0.02032439312236731658326326776_wp
      Filter(5,1) = -0.01853708698088900324035423796_wp
      Filter(5,2) = 0.04566498248660745525309655423_wp
      Filter(5,3) = -0.06014874919408893340950508358_wp
      Filter(5,4) = 0.07015854097957957031643027646_wp
      Filter(5,5) = 0.9230769230769230769230769231_wp
      Filter(5,6) = 0.0808719200830157547592528398_wp
      Filter(5,7) = -0.08217306090229584120399454405_wp
      Filter(5,8) = 0.0808719200830157547592528398_wp
      Filter(5,9) = -0.07692307692307692307692307692_wp
      Filter(5,10) = 0.07015854097957957031643027646_wp
      Filter(5,11) = -0.06014874919408893340950508358_wp
      Filter(5,12) = 0.04566498248660745525309655423_wp
      Filter(5,13) = -0.01853708698088900324035423796_wp
      Filter(6,1) = 0.01763195143996782567154213739_wp
      Filter(6,2) = -0.04343523632062224888931493992_wp
      Filter(6,3) = 0.05721178446529134760141233313_wp
      Filter(6,4) = -0.06673281454234261909514229724_wp
      Filter(6,5) = 0.07316704929522622495799478581_wp
      Filter(6,6) = 0.9230769230769230769230769231_wp
      Filter(6,7) = 0.07816068517111278566086211551_wp
      Filter(6,8) = -0.07692307692307692307692307692_wp
      Filter(6,9) = 0.07316704929522622495799478581_wp
      Filter(6,10) = -0.06673281454234261909514229724_wp
      Filter(6,11) = 0.05721178446529134760141233313_wp
      Filter(6,12) = -0.04343523632062224888931493992_wp
      Filter(6,13) = 0.01763195143996782567154213739_wp
      Filter(7,1) = -0.01735276442307692307692307692_wp
      Filter(7,2) = 0.04274747614287937972630807401_wp
      Filter(7,3) = -0.05630588431633431753036569037_wp
      Filter(7,4) = 0.06567615694637294256197217009_wp
      Filter(7,5) = -0.0720085110417042936363075225_wp
      Filter(7,6) = 0.07570506515340167349377758415_wp
      Filter(7,7) = 0.9230769230769230769230769231_wp
      Filter(7,8) = 0.07570506515340167349377758415_wp
      Filter(7,9) = -0.0720085110417042936363075225_wp
      Filter(7,10) = 0.06567615694637294256197217009_wp
      Filter(7,11) = -0.05630588431633431753036569037_wp
      Filter(7,12) = 0.04274747614287937972630807401_wp
      Filter(7,13) = -0.01735276442307692307692307692_wp
      Filter(8,1) = 0.01763195143996782567154213739_wp
      Filter(8,2) = -0.04343523632062224888931493992_wp
      Filter(8,3) = 0.05721178446529134760141233313_wp
      Filter(8,4) = -0.06673281454234261909514229724_wp
      Filter(8,5) = 0.07316704929522622495799478581_wp
      Filter(8,6) = -0.07692307692307692307692307692_wp
      Filter(8,7) = 0.07816068517111278566086211551_wp
      Filter(8,8) = 0.9230769230769230769230769231_wp
      Filter(8,9) = 0.07316704929522622495799478581_wp
      Filter(8,10) = -0.06673281454234261909514229724_wp
      Filter(8,11) = 0.05721178446529134760141233313_wp
      Filter(8,12) = -0.04343523632062224888931493992_wp
      Filter(8,13) = 0.01763195143996782567154213739_wp
      Filter(9,1) = -0.01853708698088900324035423796_wp
      Filter(9,2) = 0.04566498248660745525309655423_wp
      Filter(9,3) = -0.06014874919408893340950508358_wp
      Filter(9,4) = 0.07015854097957957031643027646_wp
      Filter(9,5) = -0.07692307692307692307692307692_wp
      Filter(9,6) = 0.0808719200830157547592528398_wp
      Filter(9,7) = -0.08217306090229584120399454405_wp
      Filter(9,8) = 0.0808719200830157547592528398_wp
      Filter(9,9) = 0.9230769230769230769230769231_wp
      Filter(9,10) = 0.07015854097957957031643027646_wp
      Filter(9,11) = -0.06014874919408893340950508358_wp
      Filter(9,12) = 0.04566498248660745525309655423_wp
      Filter(9,13) = -0.01853708698088900324035423796_wp
      Filter(10,1) = 0.02032439312236731658326326776_wp
      Filter(10,2) = -0.05006790208950716154035143871_wp
      Filter(10,3) = 0.06594816249714277729146337237_wp
      Filter(10,4) = -0.07692307692307692307692307692_wp
      Filter(10,5) = 0.08433983490386245540131491655_wp
      Filter(10,6) = -0.08866941704607879405543714103_wp
      Filter(10,7) = 0.09009601107058065879334019996_wp
      Filter(10,8) = -0.08866941704607879405543714103_wp
      Filter(10,9) = 0.08433983490386245540131491655_wp
      Filter(10,10) = 0.9230769230769230769230769231_wp
      Filter(10,11) = 0.06594816249714277729146337237_wp
      Filter(10,12) = -0.05006790208950716154035143871_wp
      Filter(10,13) = 0.02032439312236731658326326776_wp
      Filter(11,1) = -0.02370672352902102540102067364_wp
      Filter(11,2) = 0.0584000666277108962961689308_wp
      Filter(11,3) = -0.07692307692307692307692307692_wp
      Filter(11,4) = 0.08972440685621789811407082899_wp
      Filter(11,5) = -0.09837544159430522710834530347_wp
      Filter(11,6) = 0.1034255410596984721081816964_wp
      Filter(11,7) = -0.1050895449944481818642648043_wp
      Filter(11,8) = 0.1034255410596984721081816964_wp
      Filter(11,9) = -0.09837544159430522710834530347_wp
      Filter(11,10) = 0.08972440685621789811407082899_wp
      Filter(11,11) = 0.9230769230769230769230769231_wp
      Filter(11,12) = 0.0584000666277108962961689308_wp
      Filter(11,13) = -0.02370672352902102540102067364_wp
      Filter(12,1) = 0.03122589104635891845684205556_wp
      Filter(12,2) = -0.07692307692307692307692307692_wp
      Filter(12,3) = 0.101321113228078246852204017_wp
      Filter(12,4) = -0.1181826982232131796997580261_wp
      Filter(12,5) = 0.1295776203363043782181859723_wp
      Filter(12,6) = -0.1362294824330045413994822512_wp
      Filter(12,7) = 0.1384212659371062012978626186_wp
      Filter(12,8) = -0.1362294824330045413994822512_wp
      Filter(12,9) = 0.1295776203363043782181859723_wp
      Filter(12,10) = -0.1181826982232131796997580261_wp
      Filter(12,11) = 0.101321113228078246852204017_wp
      Filter(12,12) = 0.9230769230769230769230769231_wp
      Filter(12,13) = 0.03122589104635891845684205556_wp
      Filter(13,1) = -0.07692307692307692307692307692_wp
      Filter(13,2) = 0.1894953054991709300421621982_wp
      Filter(13,3) = -0.2495983789607200915631795105_wp
      Filter(13,4) = 0.2911358645588133904045865895_wp
      Filter(13,5) = -0.319206559769286565729778801_wp
      Filter(13,6) = 0.3355930160912697560936287713_wp
      Filter(13,7) = -0.340992340992340992340992341_wp
      Filter(13,8) = 0.3355930160912697560936287713_wp
      Filter(13,9) = -0.319206559769286565729778801_wp
      Filter(13,10) = 0.2911358645588133904045865895_wp
      Filter(13,11) = -0.2495983789607200915631795105_wp
      Filter(13,12) = 0.1894953054991709300421621982_wp
      Filter(13,13) = 0.9230769230769230769230769231_wp

    case(14)
      Filter(1,1) = 0.9285714285714285714285714286_wp
      Filter(1,2) = 0.1761578734372195286126760827_wp
      Filter(1,3) = -0.2326575139773399279277843899_wp
      Filter(1,4) = 0.2725726687455749354649408189_wp
      Filter(1,5) = -0.3007573738470459669350127593_wp
      Filter(1,6) = 0.3189626763390876189015776231_wp
      Filter(1,7) = -0.327924551867495764972789611_wp
      Filter(1,8) = 0.327924551867495764972789611_wp
      Filter(1,9) = -0.3189626763390876189015776231_wp
      Filter(1,10) = 0.3007573738470459669350127593_wp
      Filter(1,11) = -0.2725726687455749354649408189_wp
      Filter(1,12) = 0.2326575139773399279277843899_wp
      Filter(1,13) = -0.1761578734372195286126760827_wp
      Filter(1,14) = 0.07142857142857142857142857143_wp
      Filter(2,1) = 0.02896288832724149791787501733_wp
      Filter(2,2) = 0.9285714285714285714285714286_wp
      Filter(2,3) = 0.09433807034147053944811814395_wp
      Filter(2,4) = -0.1105228847231078606855420829_wp
      Filter(2,5) = 0.1219512312525698099824178033_wp
      Filter(2,6) = -0.1293331252551389619894264374_wp
      Filter(2,7) = 0.1329669904569858913396781887_wp
      Filter(2,8) = -0.1329669904569858913396781887_wp
      Filter(2,9) = 0.1293331252551389619894264374_wp
      Filter(2,10) = -0.1219512312525698099824178033_wp
      Filter(2,11) = 0.1105228847231078606855420829_wp
      Filter(2,12) = -0.09433807034147053944811814395_wp
      Filter(2,13) = 0.07142857142857142857142857143_wp
      Filter(2,14) = -0.02896288832724149791787501733_wp
      Filter(3,1) = -0.02192940485396681624792303442_wp
      Filter(3,2) = 0.05408252254746087622153410442_wp
      Filter(3,3) = 0.9285714285714285714285714286_wp
      Filter(3,4) = 0.08368298967067060135720199886_wp
      Filter(3,5) = -0.09233602299470811111831425755_wp
      Filter(3,6) = 0.09792526327842488442097363487_wp
      Filter(3,7) = -0.1006766636324113647666546051_wp
      Filter(3,8) = 0.1006766636324113647666546051_wp
      Filter(3,9) = -0.09792526327842488442097363487_wp
      Filter(3,10) = 0.09233602299470811111831425755_wp
      Filter(3,11) = -0.08368298967067060135720199886_wp
      Filter(3,12) = 0.07142857142857142857142857143_wp
      Filter(3,13) = -0.05408252254746087622153410442_wp
      Filter(3,14) = 0.02192940485396681624792303442_wp
      Filter(4,1) = 0.01871809392998563208241099158_wp
      Filter(4,2) = -0.04616275470106154501911869584_wp
      Filter(4,3) = 0.06096867280202711401906171862_wp
      Filter(4,4) = 0.9285714285714285714285714286_wp
      Filter(4,5) = 0.07881446683326734846937730915_wp
      Filter(4,6) = -0.08358522670224506852313029316_wp
      Filter(4,7) = 0.0859337158932592346958172919_wp
      Filter(4,8) = -0.0859337158932592346958172919_wp
      Filter(4,9) = 0.08358522670224506852313029316_wp
      Filter(4,10) = -0.07881446683326734846937730915_wp
      Filter(4,11) = 0.07142857142857142857142857143_wp
      Filter(4,12) = -0.06096867280202711401906171862_wp
      Filter(4,13) = 0.04616275470106154501911869584_wp
      Filter(4,14) = -0.01871809392998563208241099158_wp
      Filter(5,1) = -0.01696397581567273269937830024_wp
      Filter(5,2) = 0.04183673066621062232323418931_wp
      Filter(5,3) = -0.05525515016624589437933157151_wp
      Filter(5,4) = 0.06473482624858631397007159593_wp
      Filter(5,5) = 0.9285714285714285714285714286_wp
      Filter(5,6) = 0.07575225178525944419853549642_wp
      Filter(5,7) = -0.07788065834143723534816993725_wp
      Filter(5,8) = 0.07788065834143723534816993725_wp
      Filter(5,9) = -0.07575225178525944419853549642_wp
      Filter(5,10) = 0.07142857142857142857142857143_wp
      Filter(5,11) = -0.06473482624858631397007159593_wp
      Filter(5,12) = 0.05525515016624589437933157151_wp
      Filter(5,13) = -0.04183673066621062232323418931_wp
      Filter(5,14) = 0.01696397581567273269937830024_wp
      Filter(6,1) = 0.01599572989192809709998235013_wp
      Filter(6,2) = -0.03944883266573506512128795448_wp
      Filter(6,3) = 0.05210137451272621484167206227_wp
      Filter(6,4) = -0.06103998299246691721631834786_wp
      Filter(6,5) = 0.06735167201088181634546784977_wp
      Filter(6,6) = 0.9285714285714285714285714286_wp
      Filter(6,7) = 0.07343549579245638786569354625_wp
      Filter(6,8) = -0.07343549579245638786569354625_wp
      Filter(6,9) = 0.07142857142857142857142857143_wp
      Filter(6,10) = -0.06735167201088181634546784977_wp
      Filter(6,11) = 0.06103998299246691721631834786_wp
      Filter(6,12) = -0.05210137451272621484167206227_wp
      Filter(6,13) = 0.03944883266573506512128795448_wp
      Filter(6,14) = -0.01599572989192809709998235013_wp
      Filter(7,1) = -0.015558581348273393484445433_wp
      Filter(7,2) = 0.03837073245616560350360446829_wp
      Filter(7,3) = -0.05067749200504896284038194492_wp
      Filter(7,4) = 0.05937181655991606076187705611_wp
      Filter(7,5) = -0.06551101293929272489447382239_wp
      Filter(7,6) = 0.06947649445638568615581082639_wp
      Filter(7,7) = 0.9285714285714285714285714286_wp
      Filter(7,8) = 0.07142857142857142857142857143_wp
      Filter(7,9) = -0.06947649445638568615581082639_wp
      Filter(7,10) = 0.06551101293929272489447382239_wp
      Filter(7,11) = -0.05937181655991606076187705611_wp
      Filter(7,12) = 0.05067749200504896284038194492_wp
      Filter(7,13) = -0.03837073245616560350360446829_wp
      Filter(7,14) = 0.015558581348273393484445433_wp
      Filter(8,1) = 0.015558581348273393484445433_wp
      Filter(8,2) = -0.03837073245616560350360446829_wp
      Filter(8,3) = 0.05067749200504896284038194492_wp
      Filter(8,4) = -0.05937181655991606076187705611_wp
      Filter(8,5) = 0.06551101293929272489447382239_wp
      Filter(8,6) = -0.06947649445638568615581082639_wp
      Filter(8,7) = 0.07142857142857142857142857143_wp
      Filter(8,8) = 0.9285714285714285714285714286_wp
      Filter(8,9) = 0.06947649445638568615581082639_wp
      Filter(8,10) = -0.06551101293929272489447382239_wp
      Filter(8,11) = 0.05937181655991606076187705611_wp
      Filter(8,12) = -0.05067749200504896284038194492_wp
      Filter(8,13) = 0.03837073245616560350360446829_wp
      Filter(8,14) = -0.015558581348273393484445433_wp
      Filter(9,1) = -0.01599572989192809709998235013_wp
      Filter(9,2) = 0.03944883266573506512128795448_wp
      Filter(9,3) = -0.05210137451272621484167206227_wp
      Filter(9,4) = 0.06103998299246691721631834786_wp
      Filter(9,5) = -0.06735167201088181634546784977_wp
      Filter(9,6) = 0.07142857142857142857142857143_wp
      Filter(9,7) = -0.07343549579245638786569354625_wp
      Filter(9,8) = 0.07343549579245638786569354625_wp
      Filter(9,9) = 0.9285714285714285714285714286_wp
      Filter(9,10) = 0.06735167201088181634546784977_wp
      Filter(9,11) = -0.06103998299246691721631834786_wp
      Filter(9,12) = 0.05210137451272621484167206227_wp
      Filter(9,13) = -0.03944883266573506512128795448_wp
      Filter(9,14) = 0.01599572989192809709998235013_wp
      Filter(10,1) = 0.01696397581567273269937830024_wp
      Filter(10,2) = -0.04183673066621062232323418931_wp
      Filter(10,3) = 0.05525515016624589437933157151_wp
      Filter(10,4) = -0.06473482624858631397007159593_wp
      Filter(10,5) = 0.07142857142857142857142857143_wp
      Filter(10,6) = -0.07575225178525944419853549642_wp
      Filter(10,7) = 0.07788065834143723534816993725_wp
      Filter(10,8) = -0.07788065834143723534816993725_wp
      Filter(10,9) = 0.07575225178525944419853549642_wp
      Filter(10,10) = 0.9285714285714285714285714286_wp
      Filter(10,11) = 0.06473482624858631397007159593_wp
      Filter(10,12) = -0.05525515016624589437933157151_wp
      Filter(10,13) = 0.04183673066621062232323418931_wp
      Filter(10,14) = -0.01696397581567273269937830024_wp
      Filter(11,1) = -0.01871809392998563208241099158_wp
      Filter(11,2) = 0.04616275470106154501911869584_wp
      Filter(11,3) = -0.06096867280202711401906171862_wp
      Filter(11,4) = 0.07142857142857142857142857143_wp
      Filter(11,5) = -0.07881446683326734846937730915_wp
      Filter(11,6) = 0.08358522670224506852313029316_wp
      Filter(11,7) = -0.0859337158932592346958172919_wp
      Filter(11,8) = 0.0859337158932592346958172919_wp
      Filter(11,9) = -0.08358522670224506852313029316_wp
      Filter(11,10) = 0.07881446683326734846937730915_wp
      Filter(11,11) = 0.9285714285714285714285714286_wp
      Filter(11,12) = 0.06096867280202711401906171862_wp
      Filter(11,13) = -0.04616275470106154501911869584_wp
      Filter(11,14) = 0.01871809392998563208241099158_wp
      Filter(12,1) = 0.02192940485396681624792303442_wp
      Filter(12,2) = -0.05408252254746087622153410442_wp
      Filter(12,3) = 0.07142857142857142857142857143_wp
      Filter(12,4) = -0.08368298967067060135720199886_wp
      Filter(12,5) = 0.09233602299470811111831425755_wp
      Filter(12,6) = -0.09792526327842488442097363487_wp
      Filter(12,7) = 0.1006766636324113647666546051_wp
      Filter(12,8) = -0.1006766636324113647666546051_wp
      Filter(12,9) = 0.09792526327842488442097363487_wp
      Filter(12,10) = -0.09233602299470811111831425755_wp
      Filter(12,11) = 0.08368298967067060135720199886_wp
      Filter(12,12) = 0.9285714285714285714285714286_wp
      Filter(12,13) = 0.05408252254746087622153410442_wp
      Filter(12,14) = -0.02192940485396681624792303442_wp
      Filter(13,1) = -0.02896288832724149791787501733_wp
      Filter(13,2) = 0.07142857142857142857142857143_wp
      Filter(13,3) = -0.09433807034147053944811814395_wp
      Filter(13,4) = 0.1105228847231078606855420829_wp
      Filter(13,5) = -0.1219512312525698099824178033_wp
      Filter(13,6) = 0.1293331252551389619894264374_wp
      Filter(13,7) = -0.1329669904569858913396781887_wp
      Filter(13,8) = 0.1329669904569858913396781887_wp
      Filter(13,9) = -0.1293331252551389619894264374_wp
      Filter(13,10) = 0.1219512312525698099824178033_wp
      Filter(13,11) = -0.1105228847231078606855420829_wp
      Filter(13,12) = 0.09433807034147053944811814395_wp
      Filter(13,13) = 0.9285714285714285714285714286_wp
      Filter(13,14) = 0.02896288832724149791787501733_wp
      Filter(14,1) = 0.07142857142857142857142857143_wp
      Filter(14,2) = -0.1761578734372195286126760827_wp
      Filter(14,3) = 0.2326575139773399279277843899_wp
      Filter(14,4) = -0.2725726687455749354649408189_wp
      Filter(14,5) = 0.3007573738470459669350127593_wp
      Filter(14,6) = -0.3189626763390876189015776231_wp
      Filter(14,7) = 0.327924551867495764972789611_wp
      Filter(14,8) = -0.327924551867495764972789611_wp
      Filter(14,9) = 0.3189626763390876189015776231_wp
      Filter(14,10) = -0.3007573738470459669350127593_wp
      Filter(14,11) = 0.2725726687455749354649408189_wp
      Filter(14,12) = -0.2326575139773399279277843899_wp
      Filter(14,13) = 0.1761578734372195286126760827_wp
      Filter(14,14) = 0.9285714285714285714285714286_wp

    case(15)

      Filter(1,1) = 0.9333333333333333333333333333_wp
      Filter(1,2) = 0.1645618933613750778294877144_wp
      Filter(1,3) = -0.2178103903674059459393184838_wp
      Filter(1,4) = 0.2560705499944623680495691075_wp
      Filter(1,5) = -0.2839633228400229675647994315_wp
      Filter(1,6) = 0.3031952783001183757767754366_wp
      Filter(1,7) = -0.3145171009116193712441774357_wp
      Filter(1,8) = 0.3182595182595182595182595183_wp
      Filter(1,9) = -0.3145171009116193712441774357_wp
      Filter(1,10) = 0.3031952783001183757767754366_wp
      Filter(1,11) = -0.2839633228400229675647994315_wp
      Filter(1,12) = 0.2560705499944623680495691075_wp
      Filter(1,13) = -0.2178103903674059459393184838_wp
      Filter(1,14) = 0.1645618933613750778294877144_wp
      Filter(1,15) = -0.06666666666666666666666666667_wp
      Filter(2,1) = 0.02700773765822274021972884958_wp
      Filter(2,2) = 0.9333333333333333333333333333_wp
      Filter(2,3) = 0.0882384882341697773419215975_wp
      Filter(2,4) = -0.1037382935437087529941208626_wp
      Filter(2,5) = 0.1150381039173082479541360069_wp
      Filter(2,6) = -0.1228292780331214659935029654_wp
      Filter(2,7) = 0.1274159302566867634601878215_wp
      Filter(2,8) = -0.1289320436457812866433675616_wp
      Filter(2,9) = 0.1274159302566867634601878215_wp
      Filter(2,10) = -0.1228292780331214659935029654_wp
      Filter(2,11) = 0.1150381039173082479541360069_wp
      Filter(2,12) = -0.1037382935437087529941208626_wp
      Filter(2,13) = 0.0882384882341697773419215975_wp
      Filter(2,14) = -0.06666666666666666666666666667_wp
      Filter(2,15) = 0.02700773765822274021972884958_wp
      Filter(3,1) = -0.02040510756602330401757284324_wp
      Filter(3,2) = 0.05036854702960972089055548688_wp
      Filter(3,3) = 0.9333333333333333333333333333_wp
      Filter(3,4) = 0.07837720675691629195481246765_wp
      Filter(3,5) = -0.08691453221034106124003652492_wp
      Filter(3,6) = 0.09280098400836430134527730964_wp
      Filter(3,7) = -0.0962663291318309916229999787_wp
      Filter(3,8) = 0.09741179555994341871326149872_wp
      Filter(3,9) = -0.0962663291318309916229999787_wp
      Filter(3,10) = 0.09280098400836430134527730964_wp
      Filter(3,11) = -0.08691453221034106124003652492_wp
      Filter(3,12) = 0.07837720675691629195481246765_wp
      Filter(3,13) = -0.06666666666666666666666666667_wp
      Filter(3,14) = 0.05036854702960972089055548688_wp
      Filter(3,15) = -0.02040510756602330401757284324_wp
      Filter(4,1) = 0.01735632795157645969304651008_wp
      Filter(4,2) = -0.04284285284268568358491483247_wp
      Filter(4,3) = 0.05670582849716381834465585576_wp
      Filter(4,4) = 0.9333333333333333333333333333_wp
      Filter(4,5) = 0.07392840836146231115021582192_wp
      Filter(4,6) = -0.07893535025319522283747809426_wp
      Filter(4,7) = 0.08188292924701699957871629215_wp
      Filter(4,8) = -0.08285724858934403135514977305_wp
      Filter(4,9) = 0.08188292924701699957871629215_wp
      Filter(4,10) = -0.07893535025319522283747809426_wp
      Filter(4,11) = 0.07392840836146231115021582192_wp
      Filter(4,12) = -0.06666666666666666666666666667_wp
      Filter(4,13) = 0.05670582849716381834465585576_wp
      Filter(4,14) = -0.04284285284268568358491483247_wp
      Filter(4,15) = 0.01735632795157645969304651008_wp
      Filter(5,1) = -0.01565147357762227884842159151_wp
      Filter(5,2) = 0.03863454188743585692794308052_wp
      Filter(5,3) = -0.05113580354650572425691030562_wp
      Filter(5,4) = 0.06011822170868298805220876297_wp
      Filter(5,5) = 0.9333333333333333333333333333_wp
      Filter(5,6) = 0.07118179330761204364604728819_wp
      Filter(5,7) = -0.07383984141942855811010665111_wp
      Filter(5,8) = 0.07471845661298467851181216646_wp
      Filter(5,9) = -0.07383984141942855811010665111_wp
      Filter(5,10) = 0.07118179330761204364604728819_wp
      Filter(5,11) = -0.06666666666666666666666666667_wp
      Filter(5,12) = 0.06011822170868298805220876297_wp
      Filter(5,13) = -0.05113580354650572425691030562_wp
      Filter(5,14) = 0.03863454188743585692794308052_wp
      Filter(5,15) = -0.01565147357762227884842159151_wp
      Filter(6,1) = 0.0146586862083158938530146222_wp
      Filter(6,2) = -0.03618391734946109555771016429_wp
      Filter(6,3) = 0.0478922124795989183387649042_wp
      Filter(6,4) = -0.05630486759339536655376844046_wp
      Filter(6,5) = 0.06243793866273897499818190048_wp
      Filter(6,6) = 0.9333333333333333333333333333_wp
      Filter(6,7) = 0.06915611234118979684572923645_wp
      Filter(6,8) = -0.06997899616464091051509078382_wp
      Filter(6,9) = 0.06915611234118979684572923645_wp
      Filter(6,10) = -0.06666666666666666666666666667_wp
      Filter(6,11) = 0.06243793866273897499818190048_wp
      Filter(6,12) = -0.05630486759339536655376844046_wp
      Filter(6,13) = 0.0478922124795989183387649042_wp
      Filter(6,14) = -0.03618391734946109555771016429_wp
      Filter(6,15) = 0.0146586862083158938530146222_wp
      Filter(7,1) = -0.01413101046513000904935767964_wp
      Filter(7,2) = 0.03488138755876799680372638905_wp
      Filter(7,3) = -0.04616821358543798908145437949_wp
      Filter(7,4) = 0.05427803432675017355595631391_wp
      Filter(7,5) = -0.06019033030148183815730123395_wp
      Filter(7,6) = 0.06426683475955467463112927555_wp
      Filter(7,7) = 0.9333333333333333333333333333_wp
      Filter(7,8) = 0.06745992874728731592793596248_wp
      Filter(7,9) = -0.06666666666666666666666666667_wp
      Filter(7,10) = 0.06426683475955467463112927555_wp
      Filter(7,11) = -0.06019033030148183815730123395_wp
      Filter(7,12) = 0.05427803432675017355595631391_wp
      Filter(7,13) = -0.04616821358543798908145437949_wp
      Filter(7,14) = 0.03488138755876799680372638905_wp
      Filter(7,15) = -0.01413101046513000904935767964_wp
      Filter(8,1) = 0.01396484375_wp
      Filter(8,2) = -0.03447121691993647870549327611_wp
      Filter(8,3) = 0.04562532102910993691795294412_wp
      Filter(8,4) = -0.05363977829473845502600837263_wp
      Filter(8,5) = 0.05948255151287590482680613092_wp
      Filter(8,6) = -0.06351112030798378086339680777_wp
      Filter(8,7) = 0.06588273256400620618347271479_wp
      Filter(8,8) = 0.9333333333333333333333333333_wp
      Filter(8,9) = 0.06588273256400620618347271479_wp
      Filter(8,10) = -0.06351112030798378086339680777_wp
      Filter(8,11) = 0.05948255151287590482680613092_wp
      Filter(8,12) = -0.05363977829473845502600837263_wp
      Filter(8,13) = 0.04562532102910993691795294412_wp
      Filter(8,14) = -0.03447121691993647870549327611_wp
      Filter(8,15) = 0.01396484375_wp
      Filter(9,1) = -0.01413101046513000904935767964_wp
      Filter(9,2) = 0.03488138755876799680372638905_wp
      Filter(9,3) = -0.04616821358543798908145437949_wp
      Filter(9,4) = 0.05427803432675017355595631391_wp
      Filter(9,5) = -0.06019033030148183815730123395_wp
      Filter(9,6) = 0.06426683475955467463112927555_wp
      Filter(9,7) = -0.06666666666666666666666666667_wp
      Filter(9,8) = 0.06745992874728731592793596248_wp
      Filter(9,9) = 0.9333333333333333333333333333_wp
      Filter(9,10) = 0.06426683475955467463112927555_wp
      Filter(9,11) = -0.06019033030148183815730123395_wp
      Filter(9,12) = 0.05427803432675017355595631391_wp
      Filter(9,13) = -0.04616821358543798908145437949_wp
      Filter(9,14) = 0.03488138755876799680372638905_wp
      Filter(9,15) = -0.01413101046513000904935767964_wp
      Filter(10,1) = 0.0146586862083158938530146222_wp
      Filter(10,2) = -0.03618391734946109555771016429_wp
      Filter(10,3) = 0.0478922124795989183387649042_wp
      Filter(10,4) = -0.05630486759339536655376844046_wp
      Filter(10,5) = 0.06243793866273897499818190048_wp
      Filter(10,6) = -0.06666666666666666666666666667_wp
      Filter(10,7) = 0.06915611234118979684572923645_wp
      Filter(10,8) = -0.06997899616464091051509078382_wp
      Filter(10,9) = 0.06915611234118979684572923645_wp
      Filter(10,10) = 0.9333333333333333333333333333_wp
      Filter(10,11) = 0.06243793866273897499818190048_wp
      Filter(10,12) = -0.05630486759339536655376844046_wp
      Filter(10,13) = 0.0478922124795989183387649042_wp
      Filter(10,14) = -0.03618391734946109555771016429_wp
      Filter(10,15) = 0.0146586862083158938530146222_wp
      Filter(11,1) = -0.01565147357762227884842159151_wp
      Filter(11,2) = 0.03863454188743585692794308052_wp
      Filter(11,3) = -0.05113580354650572425691030562_wp
      Filter(11,4) = 0.06011822170868298805220876297_wp
      Filter(11,5) = -0.06666666666666666666666666667_wp
      Filter(11,6) = 0.07118179330761204364604728819_wp
      Filter(11,7) = -0.07383984141942855811010665111_wp
      Filter(11,8) = 0.07471845661298467851181216646_wp
      Filter(11,9) = -0.07383984141942855811010665111_wp
      Filter(11,10) = 0.07118179330761204364604728819_wp
      Filter(11,11) = 0.9333333333333333333333333333_wp
      Filter(11,12) = 0.06011822170868298805220876297_wp
      Filter(11,13) = -0.05113580354650572425691030562_wp
      Filter(11,14) = 0.03863454188743585692794308052_wp
      Filter(11,15) = -0.01565147357762227884842159151_wp
      Filter(12,1) = 0.01735632795157645969304651008_wp
      Filter(12,2) = -0.04284285284268568358491483247_wp
      Filter(12,3) = 0.05670582849716381834465585576_wp
      Filter(12,4) = -0.06666666666666666666666666667_wp
      Filter(12,5) = 0.07392840836146231115021582192_wp
      Filter(12,6) = -0.07893535025319522283747809426_wp
      Filter(12,7) = 0.08188292924701699957871629215_wp
      Filter(12,8) = -0.08285724858934403135514977305_wp
      Filter(12,9) = 0.08188292924701699957871629215_wp
      Filter(12,10) = -0.07893535025319522283747809426_wp
      Filter(12,11) = 0.07392840836146231115021582192_wp
      Filter(12,12) = 0.9333333333333333333333333333_wp
      Filter(12,13) = 0.05670582849716381834465585576_wp
      Filter(12,14) = -0.04284285284268568358491483247_wp
      Filter(12,15) = 0.01735632795157645969304651008_wp
      Filter(13,1) = -0.02040510756602330401757284324_wp
      Filter(13,2) = 0.05036854702960972089055548688_wp
      Filter(13,3) = -0.06666666666666666666666666667_wp
      Filter(13,4) = 0.07837720675691629195481246765_wp
      Filter(13,5) = -0.08691453221034106124003652492_wp
      Filter(13,6) = 0.09280098400836430134527730964_wp
      Filter(13,7) = -0.0962663291318309916229999787_wp
      Filter(13,8) = 0.09741179555994341871326149872_wp
      Filter(13,9) = -0.0962663291318309916229999787_wp
      Filter(13,10) = 0.09280098400836430134527730964_wp
      Filter(13,11) = -0.08691453221034106124003652492_wp
      Filter(13,12) = 0.07837720675691629195481246765_wp
      Filter(13,13) = 0.9333333333333333333333333333_wp
      Filter(13,14) = 0.05036854702960972089055548688_wp
      Filter(13,15) = -0.02040510756602330401757284324_wp
      Filter(14,1) = 0.02700773765822274021972884958_wp
      Filter(14,2) = -0.06666666666666666666666666667_wp
      Filter(14,3) = 0.0882384882341697773419215975_wp
      Filter(14,4) = -0.1037382935437087529941208626_wp
      Filter(14,5) = 0.1150381039173082479541360069_wp
      Filter(14,6) = -0.1228292780331214659935029654_wp
      Filter(14,7) = 0.1274159302566867634601878215_wp
      Filter(14,8) = -0.1289320436457812866433675616_wp
      Filter(14,9) = 0.1274159302566867634601878215_wp
      Filter(14,10) = -0.1228292780331214659935029654_wp
      Filter(14,11) = 0.1150381039173082479541360069_wp
      Filter(14,12) = -0.1037382935437087529941208626_wp
      Filter(14,13) = 0.0882384882341697773419215975_wp
      Filter(14,14) = 0.9333333333333333333333333333_wp
      Filter(14,15) = 0.02700773765822274021972884958_wp
      Filter(15,1) = -0.06666666666666666666666666667_wp
      Filter(15,2) = 0.1645618933613750778294877144_wp
      Filter(15,3) = -0.2178103903674059459393184838_wp
      Filter(15,4) = 0.2560705499944623680495691075_wp
      Filter(15,5) = -0.2839633228400229675647994315_wp
      Filter(15,6) = 0.3031952783001183757767754366_wp
      Filter(15,7) = -0.3145171009116193712441774357_wp
      Filter(15,8) = 0.3182595182595182595182595183_wp
      Filter(15,9) = -0.3145171009116193712441774357_wp
      Filter(15,10) = 0.3031952783001183757767754366_wp
      Filter(15,11) = -0.2839633228400229675647994315_wp
      Filter(15,12) = 0.2560705499944623680495691075_wp
      Filter(15,13) = -0.2178103903674059459393184838_wp
      Filter(15,14) = 0.1645618933613750778294877144_wp
      Filter(15,15) = 0.9333333333333333333333333333_wp

    case(16)
      Filter(:,:) = 0.0_wp

    case(17)
      Filter(:,:) = 0.0_wp
    
    case(18)
      Filter(:,:) = 0.0_wp

    case default
      write(*,*)'Filter only defined for p <= 16'
      write(*,*)'stopping'
      stop
    end select


    return
  end subroutine Filter_GLL_2_GLL

  subroutine Get_Ext_SSSCE_S2F(N_Flux_Pts,coeff)

    implicit none
    integer,                  intent(in)  :: N_Flux_Pts
    real(wp), dimension(:,:), intent(out) :: coeff

    real(wp), dimension( 3, 2)  :: intrp_2
    real(wp), dimension( 4, 3)  :: intrp_3
    real(wp), dimension( 5, 4)  :: intrp_4
    real(wp), dimension( 6, 5)  :: intrp_5
    real(wp), dimension( 7, 6)  :: intrp_6
    !real(wp), dimension( 8, 7)  :: intrp_7
    !real(wp), dimension( 9, 8)  :: intrp_8
    !real(wp), dimension(10, 9)  :: intrp_9
    !real(wp), dimension(11,10)  :: intrp_10
    !real(wp) :: gaI, rootT, rootroe

    intrp_2(1,1) =  1.0_wp
    intrp_2(1,2) =  0.0_wp
    intrp_2(2,1) =  0.5_wp
    intrp_2(2,2) =  0.5_wp
    intrp_2(3,1) =  0.0_wp
    intrp_2(3,2) =  1.0_wp

    intrp_3(1,1) =  1.0_wp
    intrp_3(1,2) =  0.0_wp
    intrp_3(1,3) =  0.0_wp
    intrp_3(2,1) =  0.50_wp
    intrp_3(2,2) =  0.6666666666666666666666666667_wp
    intrp_3(2,3) = -0.1666666666666666666666666667_wp
    intrp_3(3,1) = -0.1666666666666666666666666667_wp
    intrp_3(3,2) =  0.6666666666666666666666666667_wp
    intrp_3(3,3) =  0.50_wp
    intrp_3(4,1) =  0.0_wp
    intrp_3(4,2) =  0.0_wp
    intrp_3(4,3) =  1.0_wp

    intrp_4(1,1) =  1.0_wp
    intrp_4(1,2) =  0.0_wp
    intrp_4(1,3) =  0.0_wp
    intrp_4(1,4) =  0.0_wp
    intrp_4(2,1) =  0.50_wp
    intrp_4(2,2) =  0.6741808286457895200852445143_wp
    intrp_4(2,3) = -0.2575141619791228534185778477_wp
    intrp_4(2,4) =  0.08333333333333333333333333333_wp
    intrp_4(3,1) = -0.1741808286457895200852445143_wp
    intrp_4(3,2) =  0.6741808286457895200852445143_wp
    intrp_4(3,3) =  0.6741808286457895200852445143_wp
    intrp_4(3,4) = -0.1741808286457895200852445143_wp
    intrp_4(4,1) =  0.08333333333333333333333333333_wp
    intrp_4(4,2) = -0.2575141619791228534185778477_wp
    intrp_4(4,3) =  0.6741808286457895200852445143_wp
    intrp_4(4,4) =  0.50_wp
    intrp_4(5,1) =  0.0_wp
    intrp_4(5,2) =  0.0_wp
    intrp_4(5,3) =  0.0_wp
    intrp_4(5,4) =  1.0_wp

    intrp_5(1,1) =  1.0_wp
    intrp_5(1,2) =  0.0_wp
    intrp_5(1,3) =  0.0_wp
    intrp_5(1,4) =  0.0_wp
    intrp_5(1,5) =  0.0_wp
    intrp_5(2,1) =  0.50_wp
    intrp_5(2,2) =  0.675650248872424000384302753_wp
    intrp_5(2,3) = -0.2666666666666666666666666667_wp
    intrp_5(2,4) =  0.1410164177942426662823639137_wp
    intrp_5(2,5) = -0.050_wp
    intrp_5(3,1) = -0.175650248872424000384302753_wp
    intrp_5(3,2) =  0.675650248872424000384302753_wp
    intrp_5(3,3) =  0.6837934774723223717367801587_wp
    intrp_5(3,4) = -0.2748098952665650380191440724_wp
    intrp_5(3,5) =  0.0910164177942426662823639137_wp
    intrp_5(4,1) =  0.0910164177942426662823639137_wp
    intrp_5(4,2) = -0.2748098952665650380191440724_wp
    intrp_5(4,3) =  0.6837934774723223717367801587_wp
    intrp_5(4,4) =  0.675650248872424000384302753_wp
    intrp_5(4,5) = -0.175650248872424000384302753_wp
    intrp_5(5,1) = -0.050_wp
    intrp_5(5,2) =  0.1410164177942426662823639137_wp
    intrp_5(5,3) = -0.2666666666666666666666666667_wp
    intrp_5(5,4) =  0.675650248872424000384302753_wp
    intrp_5(5,5) =  0.50_wp
    intrp_5(6,1) =  0.0_wp
    intrp_5(6,2) =  0.0_wp
    intrp_5(6,3) =  0.0_wp
    intrp_5(6,4) =  0.0_wp
    intrp_5(6,5) =  1.0_wp

    intrp_6(1,1) =  1.0_wp
    intrp_6(1,2) =  0.0_wp
    intrp_6(1,3) =  0.0_wp
    intrp_6(1,4) =  0.0_wp
    intrp_6(1,5) =  0.0_wp
    intrp_6(1,6) =  0.0_wp
    intrp_6(2,1) =  0.50_wp
    intrp_6(2,2) =  0.6760943957546446186823019514_wp
    intrp_6(2,3) = -0.2690791513536898670183019099_wp
    intrp_6(2,4) =  0.1496456432117444549514198094_wp
    intrp_6(2,5) = -0.08999422094603253994875318422_wp
    intrp_6(2,6) =  0.03333333333333333333333333333_wp
    intrp_6(3,1) = -0.1760943957546446186823019514_wp
    intrp_6(3,2) =  0.6760943957546446186823019514_wp
    intrp_6(3,3) =  0.6859746879547401163848011176_wp
    intrp_6(3,4) = -0.286670943709068812772138847_wp
    intrp_6(3,5) =  0.1573571433670279030027575803_wp
    intrp_6(3,6) = -0.05666088761269920661541985088_wp
    intrp_6(4,1) =  0.0929847555990452483359999585_wp
    intrp_6(4,2) = -0.2789594435537853647208010761_wp
    intrp_6(4,3) =  0.6859746879547401163848011176_wp
    intrp_6(4,4) =  0.6859746879547401163848011176_wp
    intrp_6(4,5) = -0.2789594435537853647208010761_wp
    intrp_6(4,6) =  0.0929847555990452483359999585_wp
    intrp_6(5,1) = -0.05666088761269920661541985088_wp
    intrp_6(5,2) =  0.1573571433670279030027575803_wp
    intrp_6(5,3) = -0.286670943709068812772138847_wp
    intrp_6(5,4) =  0.6859746879547401163848011176_wp
    intrp_6(5,5) =  0.6760943957546446186823019514_wp
    intrp_6(5,6) = -0.1760943957546446186823019514_wp
    intrp_6(6,1) =  0.03333333333333333333333333333_wp
    intrp_6(6,2) = -0.08999422094603253994875318422_wp
    intrp_6(6,3) =  0.1496456432117444549514198094_wp
    intrp_6(6,4) = -0.2690791513536898670183019099_wp
    intrp_6(6,5) =  0.6760943957546446186823019514_wp
    intrp_6(6,6) =  0.50_wp
    intrp_6(7,1) =  0.0_wp
    intrp_6(7,2) =  0.0_wp
    intrp_6(7,3) =  0.0_wp
    intrp_6(7,4) =  0.0_wp
    intrp_6(7,5) =  0.0_wp
    intrp_6(7,6) =  1.0_wp

    select case (N_Flux_Pts)

      case(2)
        coeff(:,:) = intrp_2(:,:)
      case(3)
        coeff(:,:) = intrp_3(:,:)
      case(4)
        coeff(:,:) = intrp_4(:,:)
      case(5)
        coeff(:,:) = intrp_5(:,:)
      case(6)
        coeff(:,:) = intrp_6(:,:)
!         case(7)
!           coeff(:,:) = intrp_7(:,:)
!         case(8)
!           coeff(:,:) = intrp_8(:,:)
!         case(9)
!           coeff(:,:) = intrp_9(:,:)
!         case(10)
!           coeff(:,:) = intrp_10(:,:)
      case default

    end select

  end subroutine Get_Ext_SSSCE_S2F


  subroutine compute_gsat_f2s_matrix(Int_F2S,ns, nf, gsat_f2s_1d)
    ! Build rotation matrices used to take solution data ==> flux data
    ! ns:   nbrSolPnts
    ! nf:   nbrFlxPnts

    implicit none 

    real(wp), dimension(ns,nf),  intent(in) :: Int_F2S
    integer,  intent(in)                       :: ns,nf
    real(wp), dimension(ns,nf),  intent(inout) :: gsat_f2s_1d

    real(wp), dimension(nf,nf) :: mat

    continue

    mat = 0.0_wp
    mat(1,1) = 1.0_wp
    mat(nf,nf) = 1.0_wp

!    gsat_f2s_1d = matmul(Int_F2S,mat)
    gsat_f2s_1d = Int_F2S


    return
  end subroutine compute_gsat_f2s_matrix

end module initcollocation
