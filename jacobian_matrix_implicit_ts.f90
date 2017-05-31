! This module contains the necessary routines to compute the Jacobian matrix 
! of the residual of the Navier-Stokes equations discretized with a
! discontinuous collocation method and an implicit time-stepping scheme.
! The Jacobian matrix is build for a single element.

module jacobian_matrix_implicit_ts

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Subroutines and functions in this module are usually private
  private

  ! Exceptions, i.e. public subroutines or functions
  public compute_residual_jacobian_element
  public communicate_data_jacobian

contains
  
  !============================================================================

  !============================================================================
  ! communicate_data_jacobian - Communicates data across parallel interfaces.

  subroutine communicate_data_jacobian()

    ! Load modules
    use controlvariables, only: IMEX_penalty
    use nsereferencevariables, only : viscous
    use variables, only : ug, ughst, uelemghst, velemghst, welemghst, r_x, &
      & r_x_ghst 
    use referencevariables, only : nequations, nodesperelem, ihelems, nghost, &
      & nghost_elem
    use navierstokes, only : primitivevariables, entropyvariables
    use petscvariables, only:  upetsc, ulocpetsc, uelempetsc, uelemlocpetsc, &
      & r_x_petsc,r_x_loc_petsc
    use mpimod, only : UpdateComm1DGhostData, UpdateComm1DElementGhostData, &
      & UpdateComm2DGeomGhostData 

    ! Nothing is implicitly defined
    implicit none

    integer :: i

    ! Exchange conservative variables of at the ineterfaces
    call UpdateComm1DGhostData(ug,ughst,upetsc,ulocpetsc,nequations, &
      & nodesperelem,ihelems,nghost)
    
    if (IMEX_penalty == 'implicit')  then
      if (viscous) then
        ! Exchange conservative variables of the adjoining elements
        call UpdateComm1DElementGhostData(ug,uelemghst,uelempetsc, &
          & uelemlocpetsc,nequations,nodesperelem,ihelems,nghost_elem)

        ! Compute primitive and entropy variables of the adjoining elements
        do i = 1, nghost_elem
          call primitivevariables(uelemghst(:,i),velemghst(:,i),nequations)    
          call entropyvariables(velemghst(:,i),welemghst(:,i),nequations)    
        enddo

        ! Exchange the geometrical data r_x (Jacobian of the transformation)
        call UpdateComm2DGeomGhostData(r_x,r_x_ghst,r_x_petsc,r_x_loc_petsc,3, &
          & nodesperelem,ihelems,nghost)
      endif
    endif

    return
  end subroutine
   
  !============================================================================

  !============================================================================
  ! compute_residual_jacobian_element - Drives the computation of the residual 
  ! Jacobian matrix of one element. ( dResidual/d(Conserved Variable) ) 
  !
  ! Input parameters:
  ! dt  - time-step.
  ! a_kk - diagonal coefficient of the IMEX's Butcher tableau.
  ! elem_id - global element ID number.
  ! elems_counter - local element ID on the process (it is just a local counter  
  !                 of the element on the process).
  !
  ! Output parameter:
  ! dfdu_a_elem - residual Jacobi matrix of one element.
 
  subroutine compute_residual_jacobian_element(dt,a_kk,elem_id,elems_counter)

    ! Load modules
    use CSRlocalvariables
    use nsereferencevariables
    use controlvariables, only : imex_element, imex_penalty

    use variables
    use referencevariables

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt, a_kk
    integer, intent(in) :: elem_id, elems_counter
    real(wp) :: dt_times_a_kk
    integer :: i_loc

    ! Compute dt*a_kk
    dt_times_a_kk = dt*a_kk

    ! Initialize to zero the element-wise residual Jacobian matrix 
    dfdu_a_elem = 0.0_wp

    if(IMEX_element == 'implicit') then

      ! Add time contribution
      call compute_diag_ts_residual_jacobian_element(elems_counter)

      ! Add contribution of the on element inviscid flux
      call compute_inviscid_residual_jacobian_element(dt_times_a_kk,elem_id, &
        & elems_counter)

      ! Add contribution of the on element viscous flux
      if (viscous) then
        call compute_viscous_residual_jacobian_element(dt_times_a_kk,elem_id, &
          & elems_counter)
      endif

      if(IMEX_penalty == 'implicit')  then
        ! Add contribution of the penalty inviscid flux
        call compute_inviscid_residual_sat_jacobian_element(dt_times_a_kk, &
          & elem_id,elems_counter)
  
        ! Add contribution of the viscous penalty flux
        i_loc = 0
        if (viscous) then
          call compute_viscous_residual_sat_jacobian_element(dt_times_a_kk, &
            & elem_id,elems_counter,i_loc)
        endif

      endif

    endif
    
    !write(*,*) 'elem_id', elem_id
    !do i=1,nodesperelem
    !  write(*,*) 'proc ID, node ID, x, y', myprocid, i, xg(1,i,elem_id), xg(2,i,elem_id)
    !enddo

    return
  end subroutine compute_residual_jacobian_element

  !============================================================================

  !============================================================================
  ! compute_diag_ts_residual_jacobian_element - Computes the diagonal 
  ! contribution of the time-stepping algorithm to the residual Jacobian matrix 
  ! of the element. 
  !
  ! Input parameters:
  ! elems_counter - local element ID on the processor.
  !
  ! Output parameter:
  ! dfdu_a_elem - time contribution to the residual Jacobian matrix of the 
  !               element (added to current value).

  subroutine compute_diag_ts_residual_jacobian_element(elems_counter)

    ! Load modules
    use CSRlocalvariables
    use referencevariables

    ! Nothing is implicitly defined
    implicit none
    integer, intent(in) :: elems_counter

    real(wp),dimension(nequations,nequations) :: I_n_eq
    integer :: i_node, i, pos


    ! Assemble at the element level the contribution of the identity matrix 
    ! (arising from the time discretization) to the residual Jacobian matrix.
    ! -------------------------------------------------------------------------
    ! Identity matrix used for diagonal terms arising from the time 
    ! discretization
    I_n_eq = 0.0_wp     
    do i = 1, nequations
      I_n_eq(i,i) = one 
    enddo
    
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
      do i = ia_0_matmul(pos),ia_0_matmul(pos+1)-1
        dfdu_a_elem(:,:,ka_0_matmul(i)) = dfdu_a_elem(:,:,ka_0_matmul(i)) &
          & + I_n_eq(:,:)
      enddo   
    enddo

    return
  end subroutine compute_diag_ts_residual_jacobian_element

  !============================================================================
  
  !============================================================================
  ! compute_inviscid_residual_jacobian_element - Drives the computation of the 
  ! contribution of the inviscid flux on the element to the residual Jacobian 
  ! matrix of the element. 
  !
  ! Input parameters:
  ! dt_times_a_kk - time-step multiplied by the diagonal coefficient of the RK
  !                 IMEX scheme.
  ! elem_id - global element ID number.
  ! elems_counter - local element ID on the process.
  !
  ! Output parameter:
  ! dfdu_a_elem - contribution of the inviscid flux to the residual Jacobian 
  !               matrix of the element (added to current value).

  subroutine compute_inviscid_residual_jacobian_element(dt_times_a_kk,elem_id, &
      & elems_counter)

    ! Load modules
    use CSRlocalvariables
    use variables
    use referencevariables
    use unary_mod
    use jacobian_matrix_implicit_ts_variables, only : iw_inviscid
    use precision_vars
    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer, intent(in) :: elem_id, elems_counter
    integer :: i_err
    integer :: i_node, cnt, i, pos
    integer :: dir_div

    ! INVISCID CONTRIBUTION
    ! ---------------------
    ! Compute the element contribution to the residual Jacobian matrix arising 
    ! from the component x1 of the divergence
    ! -------------------------------------------------------------------------
    
    ! Direction of the derivative
    dir_div = 1 
      
    ! Compute the node values of the Euler flux Jacobian matrix
    call compute_inviscid_flux_jacobian_element(elem_id,dir_div)

    ! Compute ([a_x1_elem] tensor product [I_5]) times the 
    ! inviscid_flux_jacobian_elem matrix
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
     a_1_matmul_elem = 0.0_wp
    iw_inviscid = 0
    call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
      & a_x1_elem,ja_x1_elem,ia_x1_elem,inviscid_flux_jacobian_elem, &
      & ja_diag_elem,ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem, &
      & ia_1_matmul_elem,size(a_x1_elem),iw_inviscid,i_err)

    ! Check error status
    call check_csr_matmul_inviscid_error(i_err,dir_div)

    cnt = 0
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
     do i = ia_x1_matmul(pos), ia_x1_matmul(pos+1) - 1
        cnt = cnt + 1
        dfdu_a_elem(:,:,ka_x1_matmul(i)) = dfdu_a_elem(:,:,ka_x1_matmul(i)) &
          & + a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
      enddo
    enddo


    ! Compute the element contribution to the residual Jacobian matrix arising 
    ! from the component x2 of the divergence
    ! -------------------------------------------------------------------------
    ! Direction of the derivative
    dir_div = 2
     
    ! Compute the node values of the Euler flux Jacobian matrix
    call compute_inviscid_flux_jacobian_element(elem_id,dir_div)

    ! Compute ([a_x2_elem] tensor product [I_5]) times the 
    ! inviscid_flux_jacobi_elem matrix
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
    a_1_matmul_elem = 0.0_wp
    iw_inviscid = 0
    call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
      & a_x2_elem,ja_x2_elem,ia_x2_elem,inviscid_flux_jacobian_elem, &
      & ja_diag_elem,ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem, &
      & ia_1_matmul_elem,size(a_x2_elem),iw_inviscid,i_err)

    ! Check error status
    call check_csr_matmul_inviscid_error(i_err,dir_div)

    cnt = 0
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
      do i = ia_x2_matmul(pos), ia_x2_matmul(pos+1)-1
        cnt = cnt + 1
        dfdu_a_elem(:,:,ka_x2_matmul(i)) = dfdu_a_elem(:,:,ka_x2_matmul(i)) &
          & + a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
      enddo
    enddo
 
    ! Compute the element contribution to the residual Jacobian matrix arising 
    ! from the component x3 of the divergence
    ! -------------------------------------------------------------------------
    if (ndim .eq. 3) then

      ! Direction of the derivative
      dir_div = 3
     
      ! Compute the node values of the Euler flux Jacobian matrix
      call compute_inviscid_flux_jacobian_element(elem_id,dir_div)

      ! Compute ([a_x3_elem] tensor product [I_5]) times the 
      ! inviscid_flux_jacobian_elem matrix
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
       a_1_matmul_elem = 0.0_wp
      iw_inviscid = 0
      call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
        & a_x3_elem,ja_x3_elem,ia_x3_elem,inviscid_flux_jacobian_elem, &
        & ja_diag_elem,ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem, &
        & ia_1_matmul_elem,size(a_x3_elem),iw_inviscid,i_err)

      ! Check error status
      call check_csr_matmul_inviscid_error(i_err,dir_div)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x3_matmul(pos), ia_x3_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x3_matmul(i)) = dfdu_a_elem(:,:,ka_x3_matmul(i)) &
            & + a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo

    endif
    
    return
  end subroutine compute_inviscid_residual_jacobian_element

  !============================================================================

  !============================================================================
  ! compute_viscous_residual_jacobian_element - Drives the computation of the 
  ! contribution of the viscous flux on the element to the residual Jacobian 
  ! matrix of the element. 
  !
  ! Input parameters:
  ! dt_times_a_kk - time-step multiplied by the diagonal coefficient of the RK
  !                 IMEX scheme.
  ! elem_id - global element ID number.
  ! elems_counter - local element ID on the process.
  !
  ! Output parameter:
  ! dfdu_a_elem - contribution of the viscous flux to the residual Jacobian 
  !               matrix of the element (added to current value).

  subroutine compute_viscous_residual_jacobian_element(dt_times_a_kk,elem_id, &
      & elems_counter)

    ! Load modules
    use CSRlocalvariables
    use variables
    use referencevariables
    use unary_mod
    use nsereferencevariables

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer, intent(in) :: elem_id, elems_counter
    
    ! Add first viscous contribution, i.e. [D_xi] d/dU([\hat{C}]) [D_xj] {W}, 
    ! i,j= 1,2,3
    call compute_hatc_contribution_element(dt_times_a_kk,elem_id,elems_counter)

    ! Add second viscous contribution, i.e. [D_xi] [\hat{C}] [D_xj] [dWdU],
    ! i,j = 1,2,3
    call compute_fix_hatc_contribution_element(dt_times_a_kk,elem_id, &
      & elems_counter)
    
!   call compute_Jacobian_Matrix_Multiply(dt_times_a_kk,elem_id, &
!     & elems_counter)
    
    return
  end subroutine compute_viscous_residual_jacobian_element

  !============================================================================

  !============================================================================
  ! compute_hatc_contribution_element - Drives the computation of the first 
  ! contribution of the viscous flux on the element to the residual Jacobian 
  ! matrix of the element, i.e. [D_xi] d/dU([\hat{C}]) [D_xj] {W}, i,j= 1,2,3.
  !
  ! Input parameters:
  ! dt_times_a_kk - time-step multiplied by the diagonal coefficient of the RK
  !                 IMEX scheme.
  ! elem_id - global element ID number.
  ! elems_counter - local element ID on the process.
  !
  ! Output parameter:
  ! dfdu_a_elem - contribution of the viscous flux to the residual Jacobian 
  !               matrix of the element (added to current value).

  subroutine compute_hatc_contribution_element(dt_times_a_kk,elem_id, &
      & elems_counter)

    ! Load modules
    use CSRlocalvariables
    use variables
    use referencevariables
    use unary_mod
    use nsereferencevariables
    use jacobian_matrix_implicit_ts_variables, only : iw_viscous_1, iw_viscous_2

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer, intent(in) :: elem_id, elems_counter
    integer :: i_err
    integer :: i_node, cnt, i, pos 
    integer :: dir_div, dir_grad
    
    ! [D_xi] d/dU([\hat{C}]) [D_xj] {W}, i,j= 1,2,3
    ! ===========================================
    ! ===========================================

    ! [D_x1] d/dU([\hat{C}]) [D_x1] {W}
    ! --------------------------------
    ! Construct d/dU([\hat{C}]) [D_x1] {W} for the "[D_x1][D_x1]" term of one 
    ! element
    dir_div  = 1 
    dir_grad = 1

    dhatcdu_gradwj_elem(:,:,:) = 0.0_wp

    call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

    ! Compute ([a_x1_elem] tensor product [I_5]) times 
    ! d/dU([\hat{C}]) [D_x1] {W}, i.e. [D_x1] d/dU([\hat{C}]) [D_x1] {W}
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
     a_1_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
      & a_x1_elem,ja_x1_elem,ia_x1_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
      & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
      & size(a_x1_elem),iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    cnt = 0
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
      do i = ia_x1_matmul(pos), ia_x1_matmul(pos+1)-1
        cnt = cnt + 1
        dfdu_a_elem(:,:,ka_x1_matmul(i)) = dfdu_a_elem(:,:,ka_x1_matmul(i)) &
          & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
      enddo
    enddo

    ! [D_x2] d/dU([\hat{C}]) [D_x2] {W}
    ! --------------------------------
    ! Construct d/dU([\hat{C}]) [D_x2] {W} for the "[D_x2][D_x2]" term of one 
    ! element
    dir_div = 2
    dir_grad = 2

    call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

    ! Compute (a_x2_elem] tensor product [I_5]) times 
    ! d/dU([\hat{C}]) [D_x2] {W}, i.e. [D_x2] d/dU([\hat{C}]) [D_x2] {W}
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
    a_1_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
      & a_x2_elem,ja_x2_elem,ia_x2_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
      & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
      & size(a_x2_elem),iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    cnt = 0
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
      do i = ia_x2_matmul(pos), ia_x2_matmul(pos+1)-1
        cnt = cnt + 1
        dfdu_a_elem(:,:,ka_x2_matmul(i)) = dfdu_a_elem(:,:,ka_x2_matmul(i)) &
          & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
      enddo
    enddo

    ! CROSS TERMS
    ! -----------
    if (crossterms) then
      ! [D_x1] d/dU([\hat{C}]) [D_x2] {W}
      ! --------------------------------
      ! Construct d/dU([\hat{C}]) [D_x2] {W} for the "[D_x1][D_x2]" term of one 
      ! element
      dir_div = 1 
      dir_grad = 2

      call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

      ! Compute ([a_x1_elem] tensor product [I_5]) times 
      ! d/dU([\hat{C}]) [D_x2] {W}, i.e. [D_x1] d/dU([\hat{C}]) [D_x2] {W}
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
        & a_x1_elem,ja_x1_elem,ia_x1_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
        & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & size(a_x1_elem),iw_viscous_1,i_err)
      
      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x1_matmul(pos), ia_x1_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x1_matmul(i)) = dfdu_a_elem(:,:,ka_x1_matmul(i)) &
            & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo

      ! [D_x2] d/dU([\hat{C}]) [D_x1] {W}
      ! --------------------------------
      ! Construct d/dU([\hat{C}]) [D_x1] {W} for the "[D_x2][D_x1]" term of one 
      ! element
      dir_div = 2
      dir_grad = 1

      call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

      ! Compute [a_x2_elem] tensor product [I_5]) times 
      ! d/dU([\hat{C}]) [D_x1] {W}, i.e. [D_x2] d/dU([\hat{C}]) [D_x1] {W}
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
        & a_x2_elem,ja_x2_elem,ia_x2_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
        & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & size(a_x2_elem),iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x2_matmul(pos), ia_x2_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x2_matmul(i)) = dfdu_a_elem(:,:,ka_x2_matmul(i)) &
            & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo

    endif

    if (ndim .eq. 3) then
      ! [D_x3] d/dU([\hat{C}]) [D_x3] {W}
      ! --------------------------------
      ! Construct d/dU([\hat{C}]) [D_x3] {W} for the "[D_x3][D_x3]" term of one 
      ! element
      dir_div = 3
      dir_grad = 3

      call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

      ! Compute ([a_x3_elem] tensor product [I_5]) times 
      ! d/dU([\hat{C}]) [D_x3] {W}, i.e. [D_x3] d/dU([\hat{C}]) [D_x3] {W}
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
        & a_x3_elem,ja_x3_elem,ia_x3_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
        & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & size(a_x3_elem),iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x3_matmul(pos), ia_x3_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x3_matmul(i)) = dfdu_a_elem(:,:,ka_x3_matmul(i)) &
            & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo


      ! 3D CROSS TERMS
      ! --------------
      if (crossterms) then
        ! [D_x1] d/dU([\hat{C}]) [D_x3] {W}
        ! --------------------------------
        ! Construct d/dU([\hat{C}]) [D_x3] {W} for the "[D_x1][D_x3]" term of 
        ! one element
        dir_div = 1 
        dir_grad = 3

        call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x1_elem] tensor product [I_5]) times 
        ! d/dU([\hat{C}]) [D_x3] {W}, i.e. [D_x1] d/dU([\hat{C}]) [D_x3] {W}
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
         a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
          & a_x1_elem,ja_x1_elem,ia_x1_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x1_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x1_matmul(pos), ia_x1_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x1_matmul(i)) = dfdu_a_elem(:,:,ka_x1_matmul(i)) &
              & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

        ! [D_x2] d/dU([\hat{C}]) [D_x3] {W}
        ! --------------------------------
        ! Construct d/dU([\hat{C}]) [D_x3] {W} for the "[D_x2][D_x3]" term of 
        ! one element
        dir_div = 2
        dir_grad = 3

        call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x2_elem] tensor product [I_5]) times 
        ! d/dU([\hat{C}]) [D_x3] {W}, i.e. [D_x2] d/dU([\hat{C}]) [D_x3] {W}
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
         a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
          & a_x2_elem,ja_x2_elem,ia_x2_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x2_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x2_matmul(pos), ia_x2_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x2_matmul(i)) = dfdu_a_elem(:,:,ka_x2_matmul(i)) &
              & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

        ! [D_x3] d/dU([\hat{C}]) [D_x1] {W}
        ! --------------------------------
        ! Construct d/dU([\hat{C}]) [D_x1] {W} for the "[D_x3][D_x1]" term of 
        ! one element
        dir_div = 3 
        dir_grad = 1

        call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x3_elem] tensor product [I_5]) times 
        ! d/dU([\hat{C}]) [D_x1] {W}, i.e. [D_x3] d/dU([\hat{C}]) [D_x1] {W}
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x3_matmul(pos), ia_x3_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x3_matmul(i)) = dfdu_a_elem(:,:,ka_x3_matmul(i)) &
              & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

        ! [D_x3] d/dU([\hat{C}]) [D_x2] {W}
        ! --------------------------------
        ! Construct d/dU([\hat{C}]) [D_x2] {W} for the "[D_x3][D_x2]" term of 
        ! one element
        dir_div = 3
        dir_grad = 2

        call compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x3_elem] tensor product [I_5]) times 
        ! d/dU([\hat{C}]) [D_x2] {W}, i.e. [D_x3] d/dU([\hat{C}]) [D_x2] {W}
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,dhatcdu_gradwj_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x3_matmul(pos), ia_x3_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x3_matmul(i)) = dfdu_a_elem(:,:,ka_x3_matmul(i)) &
              & - a_1_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

      endif
    endif
    
    return
  end subroutine compute_hatc_contribution_element

  !============================================================================

  !============================================================================
  ! compute_fix_hatc_contribution_element - Drives the computation of the second 
  ! contribution of the viscous flux on the element to the residual Jacobian 
  ! matrix of the element, i.e. [D_xi] [\hat{C}] [D_xj] [dWdU], i,j = 1,2,3.
  !
  ! Input parameters:
  ! dt_times_a_kk - time-step multiplied by the diagonal coefficient of the RK
  !                 IMEX scheme.
  ! elem_id - global element ID number.
  ! elems_counter - local element ID on the process.
  !
  ! Output parameter:
  ! dfdu_a_elem - contribution of the viscous flux to the residual Jacobian 
  !               matrix of one element (added to current value).

  subroutine compute_fix_hatc_contribution_element(dt_times_a_kk,elem_id, &
      & elems_counter)

    ! Load modules
    use CSRlocalvariables
    use variables
    use referencevariables
    use unary_mod
    use nsereferencevariables
    use jacobian_matrix_implicit_ts_variables, only : iw_viscous_1, iw_viscous_2

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer, intent(in) :: elem_id, elems_counter
    integer :: i_err
    integer :: i_node, cnt, i, pos 
    integer :: dir_div, dir_grad

    ! [D_xi] [\hat{C}] [D_xj] [dW/DU], i,j= 1,2,3
    ! ===========================================
    ! ===========================================
    ! Compute the Jacobian matrix of the transformation between entropy and 
    ! conservative variables of one element
    call compute_dwdu_element(elem_id)

    ! LAPLACIAN TERMS
    ! ---------------
    ! [D_x1] [\hat{C}] [D_x1] [dW/dU]
    ! ------------------------------
    ! Construct [\hat{C}] for the "[D_x1][D_x1]" term of one element
    dir_div = 1 
    dir_grad = 1

    call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

    ! Compute ([a_x1_elem] tensor product [I_5]) times [\hat{C}], i.e. 
    ! [D_x1] [\hat{C}]
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
    a_1_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
      & a_x1_elem,ja_x1_elem,ia_x1_elem,hatC_elem,ja_diag_elem,ia_diag_elem, &
      & a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem,size(a_x1_elem), &
      & iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    ! Compute ([a_x1_elem] tensor product [I_5]) times dwdu_elem, i.e. 
    ! [D_x1] [dW/dU]
    ia_2_matmul_elem = 0
    ja_2_matmul_elem = 0
    a_2_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
      & a_x1_elem,ja_x1_elem,ia_x1_elem,dwdu_elem,ja_diag_elem,ia_diag_elem, &
      & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,size(a_x1_elem), &
      & iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    ! Compute [D_x1] [\hat{C}] times [D_x1] [dW/dU]
    ia_lap_matmul_elem = 0
    ja_lap_matmul_elem = 0
    a_lap_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
      & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
      & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,a_lap_matmul_elem, &
      & ja_lap_matmul_elem,ia_lap_matmul_elem, &
      & size(a_x11_elem),iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    cnt = 0
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
      do i = ia_x11_matmul(pos), ia_x11_matmul(pos+1)-1
        cnt = cnt + 1
        dfdu_a_elem(:,:,ka_x11_matmul(i)) = dfdu_a_elem(:,:,ka_x11_matmul(i)) &
          & - a_lap_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
      enddo
    enddo

    ! [D_x2] [\hat{C}] [D_x2] [dW/dU]
    ! ------------------------------
    ! Construct [\hat{C}] for the "[D_x2][D_x2]" term of one element
    dir_div = 2 
    dir_grad = 2

    call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

    ! Compute ([a_x2_elem] tensor product [I_5]) times [\hat{C}], i.e. 
    ! [D_x2] [\hat{C}]
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
    a_1_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
      & a_x2_elem,ja_x2_elem,ia_x2_elem,hatC_elem,ja_diag_elem,ia_diag_elem, &
      & a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem,size(a_x2_elem), &
      & iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    ! Compute ([a_x2_elem] tensor product [I_5]) times dwdu_elem, i.e. 
    ! [D_x2] [dW/dU]
    ia_2_matmul_elem = 0
    ja_2_matmul_elem = 0
    a_2_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
      & a_x2_elem,ja_x2_elem,ia_x2_elem,dwdu_elem,ja_diag_elem,ia_diag_elem, &
      & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,size(a_x2_elem), &
      & iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    ! Compute [D_x2] [\hat{C}] times [D_x2] [dW/dU]
    ia_lap_matmul_elem = 0
    ja_lap_matmul_elem = 0
    a_lap_matmul_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
      & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
      & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,a_lap_matmul_elem, &
      & ja_lap_matmul_elem,ia_lap_matmul_elem, &
      & size(a_x22_elem),iw_viscous_1,i_err)

    ! Check error status
    call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    cnt = 0
    do i_node = 1, nodesperelem
      pos = (elems_counter - 1)*nodesperelem + i_node
      do i = ia_x22_matmul(pos), ia_x22_matmul(pos+1)-1
        cnt = cnt + 1
        dfdu_a_elem(:,:,ka_x22_matmul(i)) = dfdu_a_elem(:,:,ka_x22_matmul(i)) &
          & - a_lap_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
      enddo
    enddo

    ! CROSS TERMS
    ! -----------
    if (crossterms) then

      ! [D_x1] [\hat{C}] [D_x2] [dW/dU]
      ! ------------------------------
      ! Construct [\hat{C}] for the "[D_x1][D_x2]" term of one element
      dir_div = 1 
      dir_grad = 2

      call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

      ! Compute ([a_x1_elem] tensor product [I_5]) times [\hat{C}], i.e. 
      ! [D_x1] [\hat{C}]
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
        & a_x1_elem,ja_x1_elem,ia_x1_elem,hatC_elem,ja_diag_elem,ia_diag_elem, &
        & a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem,size(a_x1_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      ! Compute ([a_x2_elem] tensor product [I_5]) times dwdu_elem, i.e. 
      ! [D_x2] [dW/dU]
      ia_2_matmul_elem = 0
      ja_2_matmul_elem = 0
      a_2_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
        & a_x2_elem,ja_x2_elem,ia_x2_elem,dwdu_elem,ja_diag_elem,ia_diag_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,size(a_x2_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      ! Compute [D_x1] [\hat{C}] times [D_x2] [dW/dU]
      ia_ct_matmul_elem = 0
      ja_ct_matmul_elem = 0
      a_ct_matmul_elem = 0.0_wp
      iw_viscous_2 = 0
      call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
        & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,a_ct_matmul_elem, &
        & ja_ct_matmul_elem,ia_ct_matmul_elem,size(a_x12_elem),iw_viscous_2, &
        & i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x12_matmul(pos), ia_x12_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x12_matmul(i)) = dfdu_a_elem(:,:,ka_x12_matmul(i)) &
            & - a_ct_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo

      ! [D_x2] [\hat{C}] [D_x1] [dW/dU]
      ! ------------------------------
      ! Construct [\hat{C}] for the "[D_x2][D_x1]" term of one element
      dir_div = 2 
      dir_grad = 1

      call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

      ! Compute ([a_x2_elem] tensor product [I_5]) times [\hat{C}], i.e. 
      ! [D_x2] [\hat{C}]
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
        & a_x2_elem,ja_x2_elem,ia_x2_elem,hatC_elem,ja_diag_elem,ia_diag_elem, &
        & a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem,size(a_x2_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)
   
      ! Compute ([a_x1_elem] tensor product [I_5]) times dwdu_elem, i.e. 
      ! [D_x1] [dW/dU]
      ia_2_matmul_elem = 0
      ja_2_matmul_elem = 0
      a_2_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
        & a_x1_elem,ja_x1_elem,ia_x1_elem,dwdu_elem,ja_diag_elem,ia_diag_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,size(a_x1_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      ! Compute [D_x2] [\hat{C}] times [D_x1] [dW/dU]
      ia_ct_matmul_elem = 0
      ja_ct_matmul_elem = 0
      a_ct_matmul_elem = 0.0_wp
      iw_viscous_2 = 0
      call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
        & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,a_ct_matmul_elem, &
        & ja_ct_matmul_elem,ia_ct_matmul_elem,size(a_x12_elem),iw_viscous_2, &
        & i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x12_matmul(pos), ia_x12_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x21_matmul(i)) = dfdu_a_elem(:,:,ka_x21_matmul(i)) &
            & - a_ct_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo

    endif

    ! 3D LAPLACIAN TERM
    ! -----------------
    if (ndim .eq. 3) then

      ! [D_x3] [\hat{C}] [D_x3] [dW/dU]
      ! ------------------------------
      ! Construct [\hat{C}] for the "[D_x3][D_x3]" term of one element
      dir_div = 3 
      dir_grad = 3

      call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

      ! Compute ([a_x3_elem] tensor product [I_5]) times [\hat{C}], i.e. 
      ! [D_x3] [\hat{C}]
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
        & a_x3_elem,ja_x3_elem,ia_x3_elem,hatC_elem,ja_diag_elem,ia_diag_elem, &
        & a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem,size(a_x3_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)
 
      ! Compute ([a_x3_elem] tensor product [I_5]) times dwdu_elem, i.e. 
      ! [D_x3] [dW/dU]
      ia_2_matmul_elem = 0
      ja_2_matmul_elem = 0
      a_2_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
        & a_x3_elem,ja_x3_elem,ia_x3_elem,dwdu_elem,ja_diag_elem,ia_diag_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,size(a_x3_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      ! Compute [D_x3] [\hat{C}] times [D_x3] [dW/dU]
      ia_lap_matmul_elem = 0
      ja_lap_matmul_elem = 0
      a_lap_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1, &
        & size(ia_2_matmul_elem)-1,nequations,a_1_matmul_elem, &
        & ja_1_matmul_elem,ia_1_matmul_elem,a_2_matmul_elem, &
        & ja_2_matmul_elem,ia_2_matmul_elem,a_lap_matmul_elem, &
        & ja_lap_matmul_elem,ia_lap_matmul_elem,size(a_x33_elem), &
        & iw_viscous_1,i_err)

      ! Check error status
      call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

      cnt = 0
      do i_node = 1, nodesperelem
        pos = (elems_counter - 1)*nodesperelem + i_node
        do i = ia_x33_matmul(pos), ia_x33_matmul(pos+1)-1
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_x33_matmul(i)) = dfdu_a_elem(:,:,ka_x33_matmul(i)) &
            & - a_lap_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
        enddo
      enddo

      ! 3D CROSS TERMS
      ! --------------
      if (crossterms) then

        ! [D_x1] [\hat{C}] [D_x3] [dW/dU]
        ! ------------------------------
        ! Construct [\hat{C}] for the "[D_x1][D_x3]" term of one element
        dir_div = 1 
        dir_grad = 3

        call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x1_elem] tensor product [I_5]) times [\hat{C}], i.e. 
        ! [D_x1] [\hat{C}]
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
          & a_x1_elem,ja_x1_elem,ia_x1_elem,hatC_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x1_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)
    
        ! Compute ([a_x3_elem] tensor product [I_5]) times dwdu_elem, i.e. 
        ! [D_x3] [dW/dU]
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,dwdu_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        ! Compute [D_x1] [\hat{C}] times [D_x3] [dW/dU]
        ia_ct_matmul_elem = 0
        ja_ct_matmul_elem = 0
        a_ct_matmul_elem = 0.0_wp
        iw_viscous_2 = 0
        call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
          & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
          & size(a_x13_elem),iw_viscous_2,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x13_matmul(pos), ia_x13_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x13_matmul(i)) = dfdu_a_elem(:,:,ka_x13_matmul(i)) &
              & - a_ct_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

        ! [D_x2] [\hat{C}] [D_x3] [dW/dU]
        ! ------------------------------
        ! Construct [\hat{C}] for the "[D_x2][D_x3]" term of one element
        dir_div = 2 
        dir_grad = 3

        call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x2_elem] tensor product [I_5]) times [\hat{C}], i.e. 
        ! [D_x2] [\hat{C}]
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
          & a_x2_elem,ja_x2_elem,ia_x2_elem,hatC_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x2_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)
    
        ! Compute ([a_x3_elem] tensor product [I_5]) times dwdu_elem, i.e. 
        ! [D_x3] [dW/dU]
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,dwdu_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        ! Compute [D_x2] [\hat{C}] times [D_x3] [dW/dU]
        ia_ct_matmul_elem = 0
        ja_ct_matmul_elem = 0
        a_ct_matmul_elem = 0.0_wp
        iw_viscous_2 = 0
        call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
          & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
          & size(a_x23_elem),iw_viscous_2,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x23_matmul(pos), ia_x23_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x23_matmul(i)) = dfdu_a_elem(:,:,ka_x23_matmul(i)) &
              & - a_ct_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

        ! [D_x3] [\hat{C}] [D_x1] [dW/dU]
        ! ------------------------------
        ! Construct [\hat{C}] for the "[D_x3][D_x1]" term of one element
        dir_div = 3 
        dir_grad = 1

        call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x3_elem] tensor product [I_5]) times [\hat{C}], i.e. 
        ! [D_x3] [\hat{C}]
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,hatC_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)
    
        ! Compute ([a_x1_elem] tensor product [I_5]) times dwdu_elem, i.e. 
        ! [D_x1] [dW/dU]
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
          & a_x1_elem,ja_x1_elem,ia_x1_elem,dwdu_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x1_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        ! Compute [D_x3] [\hat{C}] times [D_x1] [dW/dU]
        ia_ct_matmul_elem = 0
        ja_ct_matmul_elem = 0
        a_ct_matmul_elem = 0.0_wp
        iw_viscous_2 = 0
        call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
          & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
          & size(a_x13),iw_viscous_2,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x13_matmul(pos), ia_x13_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x31_matmul(i)) = dfdu_a_elem(:,:,ka_x31_matmul(i)) &
              & - a_ct_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

        ! [D_x3] [\hat{C}] [D_x2] [dW/dU]
        ! ------------------------------
        ! Construct [\hat{C}] for the "[D_x3][D_x2]" term of one element
        dir_div = 3 
        dir_grad = 2

        call compute_matrix_hatc_element(elem_id,dir_div,dir_grad)

        ! Compute ([a_x3_elem] tensor product [I_5]) times [\hat{C}], i.e. 
        ! [D_x3] [\hat{C}]
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,hatC_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)
    
        ! Compute ([a_x2_elem] tensor product [I_5]) times dwdu_elem, i.e. 
        ! [D_x2] [dW/dU]
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
          & a_x2_elem,ja_x2_elem,ia_x2_elem,dwdu_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x2_elem),iw_viscous_1,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        ! Compute [D_x3] [\hat{C}] times [D_x2] [dW/dU]
        ia_ct_matmul_elem = 0
        ja_ct_matmul_elem = 0
        a_ct_matmul_elem = 0.0_wp
        iw_viscous_2 = 0
        call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
          & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
          & size(a_x23_elem),iw_viscous_2,i_err)

        ! Check error status
        call check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

        cnt = 0
        do i_node = 1, nodesperelem
          pos = (elems_counter - 1)*nodesperelem + i_node
          do i = ia_x23_matmul(pos), ia_x23_matmul(pos+1)-1
            cnt = cnt + 1
            dfdu_a_elem(:,:,ka_x32_matmul(i)) = dfdu_a_elem(:,:,ka_x32_matmul(i)) &
              & - a_ct_matmul_elem(:,:,cnt)/Jx_r(i_node,elem_id)*dt_times_a_kk
          enddo
        enddo

      endif
    endif

    return
  end subroutine compute_fix_hatc_contribution_element

  !============================================================================
 
  !============================================================================

  subroutine compute_Jacobian_Matrix_Multiply(dt_times_a_kk,elem_id, &
      & elems_counter)

    ! Load modules
    use CSRlocalvariables
    use variables
    use referencevariables
    use unary_mod
    use nsereferencevariables
    use jacobian_matrix_implicit_ts_variables, only : iw_viscous_1, iw_viscous_2

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer , intent(in) :: elem_id, elems_counter
    integer :: ncol, nrow
    integer :: i_err
    integer :: i_node, i, j
    integer :: dir_div, dir_grad
    integer :: npe, nq, n1, n12, n123, nFin
    integer, parameter :: bigN = 50000

    integer, dimension(bigN) :: iwork

    npe  = nodesperelem
    nq   = nequations
    n1   = nodesperelem * nodesperedge
    n12  = n1 + n1 - nodesperelem
    n123 = n1 + n1 + n1 - nodesperelem - nodesperelem
    nFin = nodesperelem * nodesperedge * ((ndim-1)*nodesperedge - 1)

    ! Transformation Jacobian matrix between entropy and conservative variables: 1 element

    call compute_dwdu_element(elem_id)

    nrow  = npe ; 
    ncol  = npe ;

!===========================================================================================
! F^{I}_dir = [D_xk] { [dF_k/dU]
! F^{v}_dir = [D_xk] { [d\hat{kd1}/dU].dWdx1 + [d\hat{kd2}/dU].dWdx2 + [d\hat{kd3}/dU].dWdx3 } 
!           + [D_xk] { [ \hat{kd1}]   [D_x1] + [ \hat{kd2}]   [D_x2] + [ \hat{kd3}]   [D_x3] } [dW/dU]
!===========================================================================================

do dir_div = 1,ndim

    a_dFvdU1_elem(:,:,:) = 0.0_wp

    ! Compute the node values of the Euler flux Jacobian matrix
    call compute_inviscid_flux_jacobian_element(elem_id,dir_div)

    ! Store in block diagonal container used to accumulate first component of viscous flux jacobian
    a_diag_tmp(:,:,:) = - inviscid_flux_jacobian_elem(:,:,:)

    dir_grad = 1 ! [\hat{C11}] [D_x1] 

    call compute_hatc_and_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

    a_dFvdU1_elem(:,:,:) = a_dFvdU1_elem(:,:,:) + dhatcdu_gradwj_elem(:,:,:) !d{\hat{C11}}/dU*dw

    ia_W1_elem = 0 ; ja_W1_elem = 0 ; a_W1_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Bl_mu_b_sc(nrow, ncol, nequations,           &
      & hatC_elem, ja_diag_elem ,ia_diag_elem,          & ! a
      & a_x1_elem, ja_x1_elem   ,ia_x1_elem,            & ! b
      & a_W1_elem, ja_W1_elem   ,ia_W1_elem,            & ! c
      & n1 ,iw_viscous_1,i_err)

    dir_grad = 2 ! [\hat{C12}] [D_x2] 

    call compute_hatc_and_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

    a_dFvdU1_elem(:,:,:) = a_dFvdU1_elem(:,:,:) + dhatcdu_gradwj_elem(:,:,:) !d{\hat{C11}}/dU*dw

    ia_W2_elem = 0 ; ja_W2_elem = 0 ; a_W2_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Bl_mu_b_sc(nrow, ncol, nequations,           &
      & hatC_elem, ja_diag_elem, ia_diag_elem,          & ! a
      & a_x2_elem, ja_x2_elem  , ia_x2_elem,            & ! b
      & a_W2_elem, ja_W2_elem  , ia_W2_elem,            & ! c
      & n1 ,iw_viscous_1     ,i_err)

    ! [\hat{C11}] [D_x1] + [\hat{C12}] [D_x2] 
    ia_W12_elem = 0 ; ja_W12_elem = 0 ; a_W12_elem = 0.0_wp
    call a_Bl_pl_b_Bl1(nrow, ncol, nequations,          &
      & a_W1_elem  ,ja_W1_elem  ,ia_W1_elem ,           & ! a
      & a_W2_elem  ,ja_W2_elem  ,ia_W2_elem ,           & ! b
      & a_W12_elem ,ja_W12_elem ,ia_W12_elem,           & ! c
      & n12,i_err)

    if(ndim .eq. 3) then

      dir_grad = 3 ! [\hat{C13}] [D_x3] 

      call compute_hatc_and_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)

      a_dFvdU1_elem(:,:,:) = a_dFvdU1_elem(:,:,:) + dhatcdu_gradwj_elem(:,:,:) !d{\hat{C11}}/dU*dw

      ia_W3_elem = 0 ; ja_W3_elem = 0 ; a_W3_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Bl_mu_b_sc(nrow, ncol, nequations,         &
        & hatC_elem, ja_diag_elem , ia_diag_elem,       & ! a
        & a_x3_elem, ja_x3_elem   , ia_x3_elem,         & ! b
        & a_W3_elem, ja_W3_elem   , ia_W3_elem,         & ! c
        & n1,iw_viscous_1,i_err)

      ! [\hat{C11}] [D_x1] + [\hat{C12}] [D_x2] + [\hat{C13}] [D_x3] 
      ia_W123_elem = 0 ; ja_W123_elem = 0 ; a_W123_elem = 0.0_wp
      call a_Bl_pl_b_Bl1(nrow, ncol, nequations,        &
        & a_W12_elem ,ja_W12_elem ,ia_W12_elem ,        & ! a
        & a_W3_elem  ,ja_W3_elem  ,ia_W3_elem  ,        & ! b
        & a_W123_elem,ja_W123_elem,ia_W123_elem,        & ! c
        & n123,i_err)
    else
      ia_w123_elem(:) = ia_w12_elem(:)
      do i = 1, npe
         do j = ia_w123_elem(i),ia_w123_elem(i+1)-1
           ja_w123_elem(j)     = ja_w12_elem(j)
            a_w123_elem(:,:,j) =  a_w12_elem(:,:,j)
        enddo
      enddo
    endif

    ! { [\hat{C11}] [D_x1] + [\hat{C12}] [D_x2] + [\hat{C13}] [D_x3] } [dW/dU]
    ia_dFvdU2_elem = 0 ; ja_dFvdU2_elem = 0 ; a_dFvdU2_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_Bl_mu_b_Bl(nrow, ncol, nequations,                              &
      & a_W123_elem       ,ja_W123_elem       ,ia_W123_elem,               & ! a
      & dwdu_elem         ,ja_diag_elem       ,ia_diag_elem,               & ! b
      & a_dFvdU2_elem     ,ja_dFvdU2_elem     ,ia_dFvdU2_elem,             & ! c
      & N123,iw_viscous_1,i_err)

    
    a_diag_tmp(:,:,:) = a_diag_tmp(:,:,:) + a_dFvdU1_elem(:,:,:)

    ! { F^{v1} + F^{v2} }  :  In place
    iw_viscous_1 = 0
    call a_Bl_pl_b_Bl_diag(nrow, nequations, 1,                            &
      & a_dFvdU2_elem     ,ja_dFvdU2_elem     ,ia_dFvdU2_elem,             & ! a
      & a_diag_tmp        ,                                                & ! diag
      & a_dFvdU2_elem     ,ja_dFvdU2_elem     ,ia_dFvdU2_elem,             & ! b
      & iw_viscous_1)

     ia_dFvdU_elem(:    ,dir_div) = ia_dFvdU2_elem(:)
     ja_dFvdU_elem(:    ,dir_div) = ja_dFvdU2_elem(:)
      a_dFvdU_elem(:,:,:,dir_div) =  a_dFvdU2_elem(:,:,:)

enddo

!                   __
!   The divergence: \/ . ( F^{I} - F^{v} )
!

    !  d / dx1  (F^{v1})
    ia_dFvdUx_elem = 0 ; ja_dFvdUx_elem = 0 ; a_dFvdUx_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_sc_mu_b_Bl(nrow, ncol, nequations,                             &
      & a_x1_elem             ,ja_x1_elem         ,ia_x1_elem        ,    & ! a
      & a_dFvdU_elem(:,:,:,1) ,ja_dFvdU_elem(:,1) ,ia_dFvdU_elem(:,1),    & ! b
      & a_dFvdUx_elem         ,ja_dFvdUx_elem     ,ia_dFvdUx_elem    ,    & ! c
      & nFin,iw_viscous_1,i_err)

    call csort_block(nequations, nrow,                                    &
      & a_dFvdUx_elem, ja_dFvdUx_elem, ia_dFvdUx_elem,                    &
      & iwork, .true.)

    !  d / dx2  (F^{v2})
    ia_dFvdUy_elem = 0 ; ja_dFvdUy_elem = 0 ; a_dFvdUy_elem = 0.0_wp
    iw_viscous_1 = 0
    call a_sc_mu_b_Bl(nrow, ncol, nequations,                             &
      & a_x2_elem             ,ja_x2_elem         ,ia_x2_elem        ,    & ! a
      & a_dFvdU_elem(:,:,:,2) ,ja_dFvdU_elem(:,2) ,ia_dFvdU_elem(:,2),    & ! b
      & a_dFvdUy_elem         ,ja_dFvdUy_elem     ,ia_dFvdUy_elem    ,    & ! c
      & nFin,iw_viscous_1,i_err)

    call csort_block(nequations, nrow,                                    &
      & a_dFvdUy_elem, ja_dFvdUy_elem, ia_dFvdUy_elem,                    &
      & iwork, .true.)

    ia_containerxy = 0 ; ja_containerxy = 0 ; a_containerxy = 0.0_wp
    call a_Bl_pl_b_Bl1(nrow, ncol, nequations,                            &
      & a_dFvdUx_elem         ,ja_dFvdUx_elem     ,ia_dFvdUx_elem    ,    & ! a
      & a_dFvdUy_elem         ,ja_dFvdUy_elem     ,ia_dFvdUy_elem    ,    & ! b
      & a_containerxy         ,ja_containerxy     ,ia_containerxy    ,    & ! c
      & bigN,i_err)

    if(ndim .eq. 3) then
      !  d / dx3  (F^{v3})
      ia_dFvdUz_elem = 0 ; ja_dFvdUz_elem = 0 ; a_dFvdUz_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_sc_mu_b_Bl(nrow, ncol, nequations,                             &
        & a_x3_elem             ,ja_x3_elem         ,ia_x3_elem        ,    & ! a
        & a_dFvdU_elem(:,:,:,3) ,ja_dFvdU_elem(:,3) ,ia_dFvdU_elem(:,3),    & ! b
        & a_dFvdUz_elem         ,ja_dFvdUz_elem     ,ia_dFvdUz_elem    ,    & ! c
        & nFin,iw_viscous_1,i_err)

      call csort_block(nequations, nrow,                                    &
        & a_dFvdUz_elem, ja_dFvdUz_elem, ia_dFvdUz_elem,                    &
        & iwork, .true.)

      ia_containerxyz = 0 ; ja_containerxyz = 0 ; a_containerxyz = 0.0_wp
      call a_Bl_pl_b_Bl1(nrow, ncol, nequations,                            &
        & a_containerxy         ,ja_containerxy     ,ia_containerxy    ,    & ! a
        & a_dFvdUz_elem         ,ja_dFvdUz_elem     ,ia_dFvdUz_elem    ,    & ! b
        & a_containerxyz        ,ja_containerxyz    ,ia_containerxyz   ,    & ! c
        & bigN,i_err)
    endif

    do i_node = 1,nrow
       do j = ia_containerxyz(i_node),ia_containerxyz(i_node+1)-1
         dfdu_a_elem(:,:,j) = - a_containerxyz(:,:,j) / Jx_r(i_node,elem_id)*dt_times_a_kk
       enddo
    enddo

    call a_Bl_pl_b_bl_scal(nrow, nequations,                                &
        & dfdu_a_elem        ,ja_containerxyz    ,ia_containerxyz      ,    & ! a
        & 1.0_wp, iw_viscous_1)

    return
  end subroutine compute_Jacobian_Matrix_Multiply

  !============================================================================

  !============================================================================
  ! compute_inviscid_flux_jacobian_element - Drives the computation of the 
  ! inviscid flux Jacobian matrix of the element using the computational space 
  ! coordinates ( dFlux/d(Conserved Variable) ).
  !
  ! Input parameters:
  ! elem_id - global element ID number.
  ! dir_div - direction of the "div" differentiation used for the contribution 
  !           of the inviscid flux to the residual Jacobian matrix.
  !
  ! Output parameter:
  ! inviscid_flux_jacobi_elem - inviscid flux Jacobian matrix of the element.

  subroutine compute_inviscid_flux_jacobian_element(elem_id,dir_div)
    
    ! Load modules
    use CSRlocalvariables
    use variables, only : ug, r_x, Jx_r
    use referencevariables
    use navierstokes

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: elem_id, dir_div
    integer :: i_node
    real(wp), dimension(3) :: n_div
    real(wp), dimension(nequations) :: v_i_node

    ! Set to zero all the elements of the block Jacobian matrix
    inviscid_flux_jacobian_elem = 0.0_wp

    ! Loop over the nodes
    do i_node = 1, nodesperelem
      ! Transform conserved variables to primitive variables
      call primitivevariables(ug(:,i_node,elem_id),v_i_node(:),nequations)

      ! Contravariant vector at the i_node
      n_div = r_x(dir_div,:,i_node,elem_id)*Jx_r(i_node,elem_id)  

      ! Euler flux Jacobi matrix of the i_node
      inviscid_flux_jacobian_elem(:,:,i_node) = &
        & inviscid_flux_jacobian_node(v_i_node,n_div,nequations) 
    enddo

    return    
  end subroutine compute_inviscid_flux_jacobian_element

  !============================================================================
 
  !============================================================================
  ! inviscid_flux_jacobian_node - Computes the inviscid flux Jacobian matrix of
  ! one node using the computational space coordinates. 
  ! ( dFlux/d(Conserved Variable) ) 
  !
  ! Input parameters:
  ! v_in  - primitive variables at the node.
  ! n_in - contravariant vector.
  ! n_eq - number of equations.
  !
  ! Output parameter:
  ! inviscid_flux_jacobi_node - inviscid flux Jacobian matrix of the node 
 
  pure function inviscid_flux_jacobian_node(v_in,n_in,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, &
                                     & Pr0, Mach0
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq
    real(wp), intent(in) :: n_in(3)
    real(wp), intent(in) :: v_in(n_eq)
    real(wp) :: inviscid_flux_jacobian_node(n_eq,n_eq)
    real(wp) :: Ke, Un, M2inv, H

    ! Kinetic energy
    Ke = half*dot_product(v_in(2:4),v_in(2:4))

    ! Normal velocity
    Un = dot_product(v_in(2:4),n_in)

    ! 1/M^2
    M2inv = 1.0_wp / Mach0 / Mach0

    ! Enthalpy
    H = (v_in(5) + gm1M2*Ke)

    ! Continuity equation
    inviscid_flux_jacobian_node(1,1) = 0.0_wp
    inviscid_flux_jacobian_node(1,2) = n_in(1)
    inviscid_flux_jacobian_node(1,3) = n_in(2)
    inviscid_flux_jacobian_node(1,4) = n_in(3)
    inviscid_flux_jacobian_node(1,5) = 0.0_wp

    ! Momentum equation (mixed up to show patterns) 
    inviscid_flux_jacobian_node(2,1) = -v_in(2)*Un  + gm1*Ke*n_in(1)
    inviscid_flux_jacobian_node(3,1) = -v_in(3)*Un  + gm1*Ke*n_in(2)
    inviscid_flux_jacobian_node(4,1) = -v_in(4)*Un  + gm1*Ke*n_in(3)

    inviscid_flux_jacobian_node(2,2) = +v_in(2)*n_in(1) - gm1*v_in(2)*n_in(1) + Un
    inviscid_flux_jacobian_node(2,3) = +v_in(2)*n_in(2) - gm1*v_in(3)*n_in(1)
    inviscid_flux_jacobian_node(2,4) = +v_in(2)*n_in(3) - gm1*v_in(4)*n_in(1)

    inviscid_flux_jacobian_node(3,2) = +v_in(3)*n_in(1) - gm1*v_in(2)*n_in(2)
    inviscid_flux_jacobian_node(3,3) = +v_in(3)*n_in(2) - gm1*v_in(3)*n_in(2) + Un
    inviscid_flux_jacobian_node(3,4) = +v_in(3)*n_in(3) - gm1*v_in(4)*n_in(2)

    inviscid_flux_jacobian_node(4,2) = +v_in(4)*n_in(1) - gm1*v_in(2)*n_in(3)
    inviscid_flux_jacobian_node(4,3) = +v_in(4)*n_in(2) - gm1*v_in(3)*n_in(3)
    inviscid_flux_jacobian_node(4,4) = +v_in(4)*n_in(3) - gm1*v_in(4)*n_in(3) + Un

    inviscid_flux_jacobian_node(2,5) =  n_in(1)*M2inv
    inviscid_flux_jacobian_node(3,5) =  n_in(2)*M2inv
    inviscid_flux_jacobian_node(4,5) =  n_in(3)*M2inv

    ! Energy equation
    inviscid_flux_jacobian_node(5,1) = + Un*(-v_in(5) + (gm1-1.0_wp)*gm1M2*Ke)
    inviscid_flux_jacobian_node(5,2) = + n_in(1)*H - gm1*gm1M2*v_in(2)*Un
    inviscid_flux_jacobian_node(5,3) = + n_in(2)*H - gm1*gm1M2*v_in(3)*Un
    inviscid_flux_jacobian_node(5,4) = + n_in(3)*H - gm1*gm1M2*v_in(4)*Un
    inviscid_flux_jacobian_node(5,5) = gamma0 * Un

    return
  end function inviscid_flux_jacobian_node

  !============================================================================

  !============================================================================
  ! compute_dwdu_element - Drives the computation of the Jacobian matrix of the 
  ! transformation between entropy and conservative variables of the element.
  ! ( d(Entropy Variables)/d(Conserved Variable) 
  !
  ! Input parameters:
  ! elem_id - global element ID number. 
  !
  ! Output parameter:
  ! dwdu_elem - Jacobi matrix of the transformation between entropy and
  !             conservative variables of the element.

  subroutine compute_dwdu_element(elem_id)
    
    ! Load modules
    use CSRlocalvariables
    use variables, only : ug
    use referencevariables
    use navierstokes, only : primitivevariables, dwdu, dudv, dvdw

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: elem_id
    integer :: i_node
    real(wp), dimension(nequations) :: v_i_node

!    real(wp), dimension(nequations,nequations) :: err

    ! Set to zero all the elements of the block Jacobi matrix
    dwdu_elem = 0.0_wp

    ! Loop over the nodes
    do i_node = 1, nodesperelem
      ! Transform conserved variables to primitive variables
      call primitivevariables(ug(:,i_node,elem_id),v_i_node(:),nequations)

      ! Jacobi matrix of the i_node
      dwdu_elem(:,:,i_node) = dwdu(v_i_node,nequations)

!      err = matmul(dwdu(v_i_node,nequations),matmul(dudv(v_i_node,nequations),dvdw(v_i_node,nequations))) - &
!        eye(nequations)
!      write(*,*) maxval(abs(err))
    enddo

    return    
  end subroutine compute_dwdu_element

  !=============================================================================

  !============================================================================
  ! compute_matrix_hatc_element - Drives the computation of the matrix [\hat{C}] 
  ! of the element using the computational space coordinates. 
  !
  ! Input parameters:
  ! elem_id - global element ID number. 
  ! dir_div - direction of the "div" differentiation used for the contribution 
  !           of the viscous flux to the residual Jacobi matrix 
  !           (see subroutine compute_residual_jacobi_element()).
  ! dir_grad - direction of the "grad" differentiation used for the contribution
  !            of the viscosu flux to the residual Jacobi matrix
  !            (see subroutine compute_residual_jacobi_element()).
  !
  ! Output parameter:
  ! hatC_elem - matrix [\hat{C}] of the element.

  subroutine compute_matrix_hatc_element(elem_id,dir_div,dir_grad)
    
    ! Load modules
    use CSRlocalvariables
    use variables, only : ug, r_x, Jx_r
    use referencevariables
    use navierstokes

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: elem_id, dir_div, dir_grad
    integer :: i_node
    real(wp), dimension(3) :: n_div, n_grad
    real(wp), dimension(nequations) :: v_i_node

!    real(wp), dimension(nequations,nequations,nodesperelem) :: test_hatc_1 

    ! Set to zero all the elements of the matrix \hat{C}
    hatC_elem = 0.0_wp

    ! Loop over the nodes
    do i_node = 1, nodesperelem
      ! Transform conserved variables to primitive variables
      call primitivevariables(ug(:,i_node,elem_id),v_i_node(:),nequations)

      ! Contravariant vector at the i_node
      n_div = r_x(dir_div,:,i_node,elem_id)*Jx_r(i_node,elem_id)  
      n_grad = r_x(dir_grad,:,i_node,elem_id)

      ! Matrix \hat{C} at the i_node
      hatc_elem(:,:,i_node) = &
        & matrix_hatc_node(v_i_node,n_div,n_grad,nequations)
    enddo

    return    
  end subroutine compute_matrix_hatc_element

  !============================================================================
 
  !============================================================================
  ! compute_dhatcdu_gradwj_element - Computes the matrix dhat{[C]}/dU*grad(W)_j 
  ! of the element using the computational space coordinates. The subscript "j"
  ! refers to the j-th component of the gradient.
  !
  ! Input parameters:
  ! elem_id - global element ID number 
  ! dir_div - direction of the "div" differentiation used for the contribution 
  !           of the viscous flux to the residual Jacobi matrix 
  !           (see subroutine compute_residual_jacobi_element()).
  ! dir_grad - direction of the "grad" differentiation used for the contribution
  !            of the viscosu flux to the residual Jacobi matrix
  !            (see subroutine compute_residual_jacobi_element()).
  !
  ! Output parameter:
  ! dhatcdu_gradwj_elem - matrix dhat{[C]}/dU * grad(W)_j of the element.

  subroutine compute_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)
    
    ! Load modules
    use CSRlocalvariables
    use variables, only : ug, r_x, Jx_r, grad_w_jacobian
    use navierstokes, only : primitivevariables
    use referencevariables

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: elem_id, dir_div, dir_grad
    integer :: i_node
    real(wp), dimension(3) :: n_div, n_grad
    real(wp), dimension(nequations) :: v_i_node

    ! Set to zero all the elements of the matrix dhat{[C]}/dU * grad(W)_j
    dhatcdu_gradwj_elem = 0.0_wp

    ! Loop over the nodes
    do i_node = 1, nodesperelem
      ! Transform conserved variables to primitive variables
      call primitivevariables(ug(:,i_node,elem_id),v_i_node(:),nequations)

      ! Contravariant vector at the i_node
      n_div = r_x(dir_div,:,i_node,elem_id)*Jx_r(i_node,elem_id)  
      n_grad = r_x(dir_grad,:,i_node,elem_id)

      ! Matrix dhat{[C]}/dU * grad(W)_j at the i_node
      dhatcdu_gradwj_elem(:,:,i_node) = &
        & dhatcdu_gradwj_node(v_i_node,&
        & grad_w_jacobian(:,dir_grad,i_node,elem_id),n_div,n_grad,nequations)

    enddo

    return    
  end subroutine compute_dhatcdu_gradwj_element

  !============================================================================

  !============================================================================
  
  subroutine compute_hatc_and_dhatcdu_gradwj_element(elem_id,dir_div,dir_grad)
    
    ! Load modules
    use CSRlocalvariables
    use variables, only : ug, r_x, Jx_r, grad_w_jacobian
    use navierstokes, only : primitivevariables
    use referencevariables

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: elem_id, dir_div, dir_grad
    integer :: i_node
    real(wp), dimension(3) :: n_div, n_grad
    real(wp), dimension(nequations) :: v_i_node

    ! Set to zero all the elements of the matrix dhat{[C]}/dU * grad(W)_j
    dhatcdu_gradwj_elem = 0.0_wp

    ! Loop over the nodes
    do i_node = 1, nodesperelem
      ! Transform conserved variables to primitive variables
      call primitivevariables(ug(:,i_node,elem_id),v_i_node(:),nequations)

      ! Contravariant vector at the i_node
      n_div = r_x(dir_div,:,i_node,elem_id)*Jx_r(i_node,elem_id)  
      n_grad = r_x(dir_grad,:,i_node,elem_id)
 
      call hatc_dhatcdu_gradwj_node(v_i_node,                        &
                      & grad_w_jacobian(:,dir_grad,i_node,elem_id),  &
                      & n_div,n_grad,nequations,                     &
                      & hatc_elem(:,:,i_node),                       &
                      & dhatcdu_gradwj_elem(:,:,i_node))

    enddo

    return    
  end subroutine compute_hatc_and_dhatcdu_gradwj_element

  !============================================================================

  !============================================================================
  ! dhatcdu_gradwj_node - Computes the matrix dhat{[C]}/dU * grad(W)_j 
  ! of one node using the computational space coordinates. The subscript "j"
  ! refers to the j-th component of the gradient.
  !
  ! Input parameters:
  ! v_in  - primitive variables at the node.
  ! gradw_j - component of the gradient in the direction used for n_j 
  !           (see below).
  ! n_i - divergence contravariant vector.
  ! n_j - gradient contravariant vector.
  ! n_eq - number of equations. 
  !
  ! Output parameter:
  ! dhatcdu_gradwj_node - matrix d\hatcdu_gradw_node of the node.

  function dhatcdu_gradwj_node(vin,dw,ni,nj,nq)
      
    ! Load modules
    use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, &
      & Pr0, Mach0, k0, mu0
    use navierstokes, only: dWdU, dVdU
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: nq
    real(wp), intent(in) :: ni(3), nj(3)
    real(wp), intent(in) :: vin(nq),dw(nq)
    real(wp) :: dhatcdu_gradwj_node(nq,nq)
    real(wp), dimension(nq,nq) :: mat
    real(wp) :: con1, con2, con3, con4
    real(wp) :: u, v, w, T

    ! Set all elements to zero
    mat = 0.0_wp

    ! Dereference primitive variables 
    u = vin(2)
    v = vin(3) 
    w = vin(4) 
    T = vin(5)

    con1 = Re0inv*mu0*T/gm1M2/3.0_wp
    con2 = gm1M2
    con3 = gm1M2
    con4 = Re0inv*k0*T*T/Pr0 

    ! Momentum equation
    mat(2,2) = con1*( 4*ni(1)*nj(1) + 3*ni(2)*nj(2) + 3*ni(3)*nj(3))
    mat(2,3) = con1*( 3*ni(2)*nj(1) - 2*ni(1)*nj(2)                )
    mat(2,4) = con1*( 3*ni(3)*nj(1)                 - 2*ni(1)*nj(3))

    mat(3,2) = con1*(-2*ni(2)*nj(1) + 3*ni(1)*nj(2)                )
    mat(3,3) = con1*( 3*ni(1)*nj(1) + 4*ni(2)*nj(2) + 3*ni(3)*nj(3))
    mat(3,4) = con1*(                 3*ni(3)*nj(2) - 2*ni(2)*nj(3))

    mat(4,2) = con1*(-2*ni(3)*nj(1)                 + 3*ni(1)*nj(3))
    mat(4,3) = con1*(               - 2*ni(3)*nj(2) + 3*ni(2)*nj(3))
    mat(4,4) = con1*( 3*ni(1)*nj(1) + 3*ni(2)*nj(2) + 4*ni(3)*nj(3))

    mat(2,5) = con2*(mat(2,2)*u + mat(2,3)*v + mat(2,4)*w) 
    mat(3,5) = con2*(mat(3,2)*u + mat(3,3)*v + mat(3,4)*w)
    mat(4,5) = con2*(mat(4,2)*u + mat(4,3)*v + mat(4,4)*w)

    ! Energy equation
    mat(5,2) = con2*(u*mat(2,2) + v*mat(3,2) + w*mat(4,2)) 
    mat(5,3) = con2*(u*mat(2,3) + v*mat(3,3) + w*mat(4,3))
    mat(5,4) = con2*(u*mat(2,4) + v*mat(3,4) + w*mat(4,4))

    mat(5,5) = con3*(mat(5,2)*u + mat(5,3)*v + mat(5,4)*w) &
             + con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3))

    dhatcdu_gradwj_node = 0.0_wp

    dhatcdu_gradwj_node(2,2) = con2 * mat(2,2) * dW(5) 
    dhatcdu_gradwj_node(3,2) = con2 * mat(3,2) * dW(5) 
    dhatcdu_gradwj_node(4,2) = con2 * mat(4,2) * dW(5) 

    dhatcdu_gradwj_node(2,3) = con2 * mat(2,3) * dW(5) 
    dhatcdu_gradwj_node(3,3) = con2 * mat(3,3) * dW(5) 
    dhatcdu_gradwj_node(4,3) = con2 * mat(4,3) * dW(5) 

    dhatcdu_gradwj_node(2,4) = con2 * mat(2,4) * dW(5) 
    dhatcdu_gradwj_node(3,4) = con2 * mat(3,4) * dW(5) 
    dhatcdu_gradwj_node(4,4) = con2 * mat(4,4) * dW(5) 

    dhatcdu_gradwj_node(5,2) = con2 * ( mat(2,2)*dW(2) + mat(2,3)*dW(3) + &
      & mat(2,4)*dW(4) + (mat(2,5)+mat(5,2))*dW(5) ) 
    dhatcdu_gradwj_node(5,3) = con2 * ( mat(3,2)*dW(2) + mat(3,3)*dW(3) + &
      & mat(3,4)*dW(4) + (mat(3,5)+mat(5,3))*dW(5) ) 
    dhatcdu_gradwj_node(5,4) = con2 * ( mat(4,2)*dW(2) + mat(4,3)*dW(3) + &
      & mat(4,4)*dW(4) + (mat(4,5)+mat(5,4))*dW(5) ) 

    dhatcdu_gradwj_node(2,5) = ( mat(2,2)*dW(2) + mat(2,3)*dW(3) + &
      & mat(2,4)*dW(4) + mat(2,5)*dW(5) ) / T
    dhatcdu_gradwj_node(3,5) = ( mat(3,2)*dW(2) + mat(3,3)*dW(3) + &
      & mat(3,4)*dW(4) + mat(3,5)*dW(5) ) / T
    dhatcdu_gradwj_node(4,5) = ( mat(4,2)*dW(2) + mat(4,3)*dW(3) + &
      & mat(4,4)*dW(4) + mat(4,5)*dW(5) ) / T
    dhatcdu_gradwj_node(5,5) = ( mat(5,2)*dW(2) + mat(5,3)*dW(3) + &
      & mat(5,4)*dW(4) + mat(5,5)*dW(5) ) / T + &
      & (con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3)))*dW(5)/T

    dhatcdu_gradwj_node = matmul(dhatcdu_gradwj_node,dVdU(vin,nq))

    return
  end function dhatcdu_gradwj_node

  !============================================================================
  
  !============================================================================
  ! Reformulate function dhatcdu_gradwj_node as a subroutine call
  !============================================================================

  subroutine hatc_dhatcdu_gradwj_node(vin,dw,ni,nj,nq,hatc,dhatc_du)
      
    ! Load modules
    use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, &
                                   & Pr0, Mach0, k0, mu0
    use navierstokes, only: dWdU, dVdU
    
    ! Nothing is implicitly defined
    implicit none

    integer,                    intent(in)    :: nq
    real(wp), dimension(3),     intent(in)    :: ni, nj
    real(wp), dimension(nq),    intent(in)    :: vin,dw
    real(wp), dimension(nq,nq), intent(inout) :: hatc, dhatc_du

    real(wp), dimension(nq,nq)                ::       dhatc_dv

    real(wp) :: con1, con2, con3, con4
    real(wp) :: u, v, w, T

    ! Set all elements to zero
    hatc = 0.0_wp

    ! Dereference primitive variables 
    u = vin(2)
    v = vin(3) 
    w = vin(4) 
    T = vin(5)

    con1 = Re0inv*mu0*T/gm1M2/3.0_wp
    con2 = gm1M2
    con3 = gm1M2
    con4 = Re0inv*k0*T*T/Pr0 

    ! Momentum equation
    hatc(2,2) = con1*( 4*ni(1)*nj(1) + 3*ni(2)*nj(2) + 3*ni(3)*nj(3))
    hatc(2,3) = con1*( 3*ni(2)*nj(1) - 2*ni(1)*nj(2)                )
    hatc(2,4) = con1*( 3*ni(3)*nj(1)                 - 2*ni(1)*nj(3))

    hatc(3,2) = con1*(-2*ni(2)*nj(1) + 3*ni(1)*nj(2)                )
    hatc(3,3) = con1*( 3*ni(1)*nj(1) + 4*ni(2)*nj(2) + 3*ni(3)*nj(3))
    hatc(3,4) = con1*(                 3*ni(3)*nj(2) - 2*ni(2)*nj(3))

    hatc(4,2) = con1*(-2*ni(3)*nj(1)                 + 3*ni(1)*nj(3))
    hatc(4,3) = con1*(               - 2*ni(3)*nj(2) + 3*ni(2)*nj(3))
    hatc(4,4) = con1*( 3*ni(1)*nj(1) + 3*ni(2)*nj(2) + 4*ni(3)*nj(3))

    hatc(2,5) = con2*(hatc(2,2)*u + hatc(2,3)*v + hatc(2,4)*w) 
    hatc(3,5) = con2*(hatc(3,2)*u + hatc(3,3)*v + hatc(3,4)*w)
    hatc(4,5) = con2*(hatc(4,2)*u + hatc(4,3)*v + hatc(4,4)*w)

    ! Energy equation
    hatc(5,2) = con2*(u*hatc(2,2) + v*hatc(3,2) + w*hatc(4,2)) 
    hatc(5,3) = con2*(u*hatc(2,3) + v*hatc(3,3) + w*hatc(4,3))
    hatc(5,4) = con2*(u*hatc(2,4) + v*hatc(3,4) + w*hatc(4,4))

    hatc(5,5) = con3*(hatc(5,2)*u + hatc(5,3)*v + hatc(5,4)*w) &
             + con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3))

    dhatc_dv = 0.0_wp

    dhatc_dv(2,2) = con2 * hatc(2,2) * dW(5) 
    dhatc_dv(3,2) = con2 * hatc(3,2) * dW(5) 
    dhatc_dv(4,2) = con2 * hatc(4,2) * dW(5) 

    dhatc_dv(2,3) = con2 * hatc(2,3) * dW(5) 
    dhatc_dv(3,3) = con2 * hatc(3,3) * dW(5) 
    dhatc_dv(4,3) = con2 * hatc(4,3) * dW(5) 

    dhatc_dv(2,4) = con2 * hatc(2,4) * dW(5) 
    dhatc_dv(3,4) = con2 * hatc(3,4) * dW(5) 
    dhatc_dv(4,4) = con2 * hatc(4,4) * dW(5) 

    dhatc_dv(5,2) = con2 * ( hatc(2,2)*dW(2) + hatc(2,3)*dW(3) + &
      & hatc(2,4)*dW(4) + (hatc(2,5)+hatc(5,2))*dW(5) ) 
    dhatc_dv(5,3) = con2 * ( hatc(3,2)*dW(2) + hatc(3,3)*dW(3) + &
      & hatc(3,4)*dW(4) + (hatc(3,5)+hatc(5,3))*dW(5) ) 
    dhatc_dv(5,4) = con2 * ( hatc(4,2)*dW(2) + hatc(4,3)*dW(3) + &
      & hatc(4,4)*dW(4) + (hatc(4,5)+hatc(5,4))*dW(5) ) 

    dhatc_dv(2,5) = ( hatc(2,2)*dW(2) + hatc(2,3)*dW(3) + &
      & hatc(2,4)*dW(4) + hatc(2,5)*dW(5) ) / T
    dhatc_dv(3,5) = ( hatc(3,2)*dW(2) + hatc(3,3)*dW(3) + &
      & hatc(3,4)*dW(4) + hatc(3,5)*dW(5) ) / T
    dhatc_dv(4,5) = ( hatc(4,2)*dW(2) + hatc(4,3)*dW(3) + &
      & hatc(4,4)*dW(4) + hatc(4,5)*dW(5) ) / T
    dhatc_dv(5,5) = ( hatc(5,2)*dW(2) + hatc(5,3)*dW(3) + &
      & hatc(5,4)*dW(4) + hatc(5,5)*dW(5) ) / T + &
      & (con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3)))*dW(5)/T

    dhatc_du = matmul(dhatc_dv,dVdU(vin,nq))

    return
  end subroutine hatc_dhatcdu_gradwj_node

  !============================================================================

  !============================================================================
  ! check_csr_matmul_inviscid_error - Checks error in the CSR matrix 
  ! multiplications for the inviscid term contributions. 
  !
  ! Input parameters:
  ! i_err  - integer, if largere than 0 than there was an error.
  ! dir_div - direction of the "div" differentiation used for the contribution 
  !           of the viscous flux to the residual Jacobian matrix. 

  subroutine check_csr_matmul_inviscid_error(i_err,dir_div)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: i_err
    integer, intent(in) :: dir_div

    ! Check error status
    if(i_err > 0) then
      write(*,*) 'Failure in jacobian_matrix_implicit_ts.f90, &
        & compute_residual_jacobian_element(); i_err > 0 in csr matmul viscous.'
      write(*,*) 'dir_div: ', dir_div 
      write(*,*) 'Stopping'
      stop
    endif
    
    return
  end subroutine check_csr_matmul_inviscid_error

  !============================================================================

  !============================================================================
  ! check_csr_matmul_viscous_error - Check error in the CSR matrix 
  ! multiplications for the viscous term contributions. 
  !
  ! Input parameters:
  ! i_err  - integer, if larger than 0 than there was an error.
  ! dir_div - direction of the "div" differentiation used for the contribution 
  !           of the viscous flux to the residual Jacobian matrix. 
  ! dir_grad - direction of the "grad" differentiation used for the contribution
  !            of the viscosu flux to the residual Jacobian matrix.

  subroutine check_csr_matmul_viscous_error(i_err,dir_div,dir_grad)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: i_err
    integer, intent(in) :: dir_div, dir_grad

    ! Check error status
    if(i_err > 0) then
      write(*,*) 'Failure in jacobian_matrix_implicit_ts.f90, &
        & compute_residual_jacobian_element(); i_err > 0 in csr matmul viscous.'
      write(*,*) 'dir_div: ', dir_div 
      write(*,*) 'dir_grad: ', dir_grad
      write(*,*) 'Stopping'
      stop
    endif

    return
  end subroutine check_csr_matmul_viscous_error

  !============================================================================

  !============================================================================
  ! eye_matrix - Construct the identity matrix of dimension "n".
  !
  ! Input parameter:
  ! n - dimension of the square matrix.

  pure function eye_matrix(n)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) ::n
    real, dimension(n,n) :: eye_matrix
    integer :: i
    
    ! Initialize to zero all the element of the matrix
    eye_matrix = 0.0_wp

    ! Set to 1 the diagonal element of the matrix
    do i = 1, n
      eye_matrix(i,i) = 1.0_wp
    enddo
    
    return
  end function eye_matrix

  !============================================================================
  
  !============================================================================
  ! compute_inviscid_residual_sat_jacobian_element - Drives the computation of
  ! the contribution of the residual Jacobian matrix of the inviscid penalty 
  ! terms to the element using complex variables.
  !
  ! Input parameters:
  ! dt_times_a_kk - time-step multiplied by the diagonal coefficient of the RK
  !                 IMEX scheme.
  ! elem_id - global element ID number.
  !
  ! Output parameters:
  ! dfdu_a_elem - contribution of the inviscid penalty to the residual Jacobian 
  !               matrix of the element (added to current value).

  subroutine compute_inviscid_residual_sat_jacobian_element(dt_times_a_kk, &
      & elem_id,elems_counter)

    ! Load modules
    use variables, only : ifacenodes, efn2efn, Jx_r, facenodenormal, ef2e, &
      & kfacenodes, ug, ughst, vg, xg
    use referencevariables, only : nequations, nfacesperelem, nodesperface 
    use CSRlocalvariables
    use referencevariables
    use collocationvariables, only : pinv

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer,  intent(in) :: elem_id, elems_counter
    
    integer :: i
    integer :: i_node, j_node, k_node, g_node, k_elem, i_face

    real(wp), dimension(3) :: n_v
    real(wp), dimension(nequations,nequations,2) :: elem_id_i_node_inviscid_sat
    real(wp), dimension(nequations) :: u_l, u_r
    integer :: cnt

    ! Set counter to access the correct element in ka_penI_proc
    cnt = (elems_counter - 1)*nodesperface*nfacesperelem*2

    ! Loop over each face
    do i_face = 1, nfacesperelem

      ! Check to see if face is a boundary face
      if (ef2e(1,i_face,elem_id) < 0) then

        ! Loop over each node on the face
        do i = 1, nodesperface
          
          ! Volumetric node index corresponding to facial node index
          i_node = kfacenodes(i,i_face)
          
          ! Volumetric index of partner node
          k_node = (i_face - 1)*nodesperface + i

          ! Outward facing normal of facial node
          n_v = Jx_r(i_node,elem_id)*facenodenormal(:,k_node,elem_id)

          ! Compute numerical Jacobian of the boundary node
          elem_id_i_node_inviscid_sat = &
            & jacobian_inviscid_bc_sat_node(elem_id,i_face,i_node,n_v, &
            & nequations)
          
          ! Update counter and add to dfdu_a_elem the contrubution of the 
          ! i_node itself
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_penI_proc(cnt)) = dfdu_a_elem(:,:,ka_penI_proc(cnt)) &
            & - elem_id_i_node_inviscid_sat(:,:,1)/Jx_r(i_node,elem_id) &
            & *dt_times_a_kk
          
          ! Update the counter 
          cnt = cnt + 1

        enddo

      ! Check if the face is off-process
      else if (ef2e(3,i_face,elem_id) /= myprocid) then

        ! Loop over each node on the face
        do i = 1, nodesperface

          ! Index in facial ordering
          j_node = nodesperface*(i_face - 1) + i
            
          ! Volumetric node index corresponding to facial node index
          i_node = ifacenodes(j_node)
          
          ! Volumetric index of partner node
          g_node = efn2efn(3,j_node,elem_id)

          ! Element index of partner node
          k_elem = efn2efn(2,j_node,elem_id)
          
          ! Outward facing normal of facial node
          n_v = Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          ! Left state
          u_l = ug(:,i_node,elem_id)
          
          ! Right state
          u_r = ughst(:,g_node)

          ! Compute numerical Jacobian of the i_node and the g_node on the 
          ! i_face for the elem_id
          elem_id_i_node_inviscid_sat = &
            & jacobian_inviscid_sat_node(u_l,u_r,n_v,nequations)
            
          ! Update the counter and add to dfdu_a_elem the contrubution of the 
          ! i_node itself
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_penI_proc(cnt)) = dfdu_a_elem(:,:,ka_penI_proc(cnt)) &
            & + elem_id_i_node_inviscid_sat(:,:,1)/Jx_r(i_node,elem_id) &
            & * dt_times_a_kk
           
          ! Update the counter and add to dfdu_a_elem the contrubution of the 
          ! g_node
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_penI_proc(cnt)) = dfdu_a_elem(:,:,ka_penI_proc(cnt)) &
            & + elem_id_i_node_inviscid_sat(:,:,2)/Jx_r(i_node,elem_id) &
            & *dt_times_a_kk

        enddo

      else

        ! Loop over each node on the face
        do i = 1, nodesperface

          ! Index in facial ordering
          j_node = nodesperface*(i_face - 1) + i
          
          ! Volumetric node index corresponding to facial node index
          i_node = ifacenodes(j_node)
        
          ! Volumetric index of partner node
          k_node = efn2efn(1,j_node,elem_id)

          ! Element index of partner node
          k_elem = efn2efn(2,j_node,elem_id)
        
          ! Outward facing normal of facial node
          n_v = Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          ! Left state
          u_l = ug(:,i_node,elem_id)
          
          ! Right state
          u_r = ug(:,k_node,k_elem)

          ! Compute numerical Jacobian of the i_node and the k_node on the i_face
          ! for the elem_id element
          elem_id_i_node_inviscid_sat = &
            & jacobian_inviscid_sat_node(u_l,u_r,n_v,nequations)

          ! Update the counter and add to dfdu_a_elem the contrubution of the 
          ! i_node itself
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_penI_proc(cnt)) = dfdu_a_elem(:,:,ka_penI_proc(cnt)) &
            & + elem_id_i_node_inviscid_sat(:,:,1)/Jx_r(i_node,elem_id) &
            & * dt_times_a_kk

          ! Update the counter and add to dfdu_a_elem the contrubution of the 
          ! k_node
          cnt = cnt + 1
          dfdu_a_elem(:,:,ka_penI_proc(cnt)) = dfdu_a_elem(:,:,ka_penI_proc(cnt)) &
            & + elem_id_i_node_inviscid_sat(:,:,2)/Jx_r(i_node,elem_id) &
            & *dt_times_a_kk

          !if (elem_id == 5) then
          !  if (k_elem == 6) then
          !    if (i_node == 2) then                
          !      stop
          !    endif
          !  endif
          !endif


        enddo

      endif

    enddo

    return
  end subroutine compute_inviscid_residual_sat_jacobian_element
  
  !============================================================================ 
  
  !============================================================================ 
  ! compute_viscous_residual_sat_jacobian_element - Drives the computation of
  ! the contribution of the residual Jacobian matrix of the viscous penalty 
  ! terms.  The matrix is analytic
  !
  ! Input parameters:
  ! dt_times_a_kk - time-step multiplied by the diagonal coefficient of the RK
  !                 IMEX scheme.
  ! elem_id - global element ID number.
  ! elems_counter - number of elements already done.
  ! i_loc - position in the ghost stack. 
  !
  ! Output parameters:
  ! dfdu_a_elem - contribution of the viscous penalty to the residual Jacobian 
  !               matrix of the element (added to current value).

  subroutine compute_viscous_residual_sat_jacobian_element(dt_times_a_kk,&
      & elem_id,elems_counter,i_loc)

    ! Load modules
    use variables, only : ifacenodes, efn2efn, efn2efn, r_x, Jx_r, & 
      & facenodenormal, ef2e, kfacenodes, ug, ughst, uelemghst, vg, wg, &
      & grad_w_jacobian, velemghst, welemghst, r_x_ghst 
    use referencevariables, only : nequations, nfacesperelem, nodesperface 
    use CSRlocalvariables
    use referencevariables
    use nsereferencevariables
    use navierstokes, only : dWdU, primitivevariables, entropyvariables
    use collocationvariables, only: iagrad,jagrad,dagrad,dmat,pinv,l01

    ! Nothing is implicitly defined
    implicit none

    real(wp), intent(in) :: dt_times_a_kk
    integer, intent(in) :: elem_id, elems_counter
    integer, intent(inout) :: i_loc
    
    integer :: i_face
    integer :: i,j, jdir
    integer :: j_node, i_node, k_node, g_node
    integer :: k_elem
    integer :: jnode
    real(wp), dimension(3)     :: n_div, n_grad
    real(wp)                   :: dt_akk_pinv_l01

    integer :: cnt1, cnt2
    
    real(wp), dimension(nequations,nequations) :: elem_id_i_node_viscous_bc_sat

    real(wp), parameter :: strength = 2.0_wp

    real(wp), dimension(nequations) :: v_ghost, w_ghost

    integer :: l
    real(wp), allocatable :: phi_tmp(:,:)

    real(wp), dimension(nequations,3) :: grad_w_ghost_node

    integer :: pos

    ! Allocate array for temporary gradient
    allocate(phi_tmp(nequations,ndim))

    ! Compute dt*a_kk*pinv(1)*l01
    dt_akk_pinv_l01 = dt_times_a_kk * pinv(1) * l01

    ! Counters
    cnt1 = (elems_counter-1) * nodesperface*nfacesperelem * 2
    cnt2 = (elems_counter-1) * nodesperface*nfacesperelem * 2 * ndim * nodesperedge

    ! Loop over each face
    do i_face = 1, nfacesperelem

      ! Check to see if face is a boundary face
      if (ef2e(1,i_face,elem_id) < 0) then

        ! Loop over each node on the face
        do i = 1, nodesperface
          
          ! Volumetric node index corresponding to facial node index
          i_node = kfacenodes(i,i_face)
          
          ! Index in facial ordering
          j_node = nodesperface*(i_face - 1) + i
         
          ! Facial index corresponding to face and node on face indices
          k_node = i_node
          
          ! Adjoining element: self
          k_elem = elem_id
          
          ! Vector needed for reusing on element subroutines
          n_div  = Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)
          
          ! Update cnt1 counter
          cnt1 = cnt1 + 1

          ! ON element Contributions
          ! ------------------------
          do jdir = 1,ndim
            n_grad  = r_x(jdir,:,i_node,elem_id)

            call hatc_dhatcdu_gradwj_node(vg(:,i_node,elem_id), &
              & grad_w_jacobian(:,jdir,i_node,elem_id), &
              & n_div,n_grad,nequations, &
              & hatc_elem(:,:,i_node), &
              & dhatcdu_gradwj_elem(:,:,i_node))

            ! Matrix dhat{[C]}/dU * grad(W)_j at the i_node
            dfdu_a_elem(:,:,ka_penI_proc(cnt1)) = dfdu_a_elem(:,:,ka_penI_proc(cnt1)) &
              & - dhatcdu_gradwj_elem(:,:,i_node) / Jx_r(i_node,elem_id) * &
              & dt_akk_pinv_l01*strength

            ! Matrix \hat{[C]} * grad_j * dWdU
            do j = iagrad(i_node), iagrad(i_node+1)-1

              jnode = jagrad(jdir,j)
              cnt2 = cnt2 + 1

              dfdu_a_elem(:,:,ka_penV_proc(cnt2)) = dfdu_a_elem(:,:,ka_penV_proc(cnt2)) &
                & - matmul(hatc_elem(:,:,i_node),dWdU(vg(:,jnode,elem_id),nequations)) * &
                & dagrad(jdir,j) / Jx_r(i_node,elem_id) * dt_akk_pinv_l01*strength
            
            end do

          end do

          ! Update cnt1 counter
          cnt1 = cnt1 + 1

          ! OFF element contributions
          ! -------------------------
          ! Compute numerical Jacobian of the i_node boundary node
          elem_id_i_node_viscous_bc_sat = &
            & jacobian_viscous_bc_sat_node(elem_id,i_face,i_node,n_div, &
            & nequations)
         
          ! Add to dfdu_a_elem the viscous contribution of the i_node itself
          ! It arises from the normal viscous flux (analytical)
          ! Sign to be checked !!!!!!!!!!
          dfdu_a_elem(:,:,ka_penI_proc(cnt1)) = dfdu_a_elem(:,:,ka_penI_proc(cnt1)) &
            & + elem_id_i_node_viscous_bc_sat(:,:)/Jx_r(i_node,elem_id) &
            & * dt_akk_pinv_l01     
         
          do jdir = 1,ndim
        
            ! Matrix \hat{[C]} * grad_j * dWdU
            do j = iagrad(k_node), iagrad(k_node+1)-1

              jnode = jagrad(jdir,j)
              cnt2 = cnt2 + 1

              dfdu_a_elem(:,:,ka_penV_proc(cnt2)) = dfdu_a_elem(:,:,ka_penV_proc(cnt2))             &
              & + matmul(hatc_elem(:,:,k_node),dWdU(vg(:,jnode,k_elem),nequations)) * dagrad(jdir,j) &
              & / Jx_r(i_node,elem_id) * dt_akk_pinv_l01

            enddo

          enddo

        enddo

        ! Check if the face is off-process
      else if (ef2e(3,i_face,elem_id) /= myprocid) then

        ! Loop over each node on the face
        do i = 1, nodesperface
          
          ! Index in facial ordering
          j_node = nodesperface*(i_face - 1) + i
            
          ! Volumetric node index corresponding to facial node index
          i_node = ifacenodes(j_node)
          
          ! Index of partner node in the ghost stack
          g_node = efn2efn(3,j_node,elem_id)

          ! Element index of partner node
          k_elem = efn2efn(2,j_node,elem_id)
          
          ! Vector needed for reusing on element subroutines
          n_div  = Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          ! Update cnt1 counter
          cnt1 = cnt1 + 1

          ! ON element Contributions
          ! ------------------------
          do jdir = 1,ndim

            n_grad = r_x(jdir,:,i_node,elem_id)

            call hatc_dhatcdu_gradwj_node(vg(:,i_node,elem_id),      &
                      & grad_w_jacobian(:,jdir,i_node,elem_id),  &
                      & n_div,n_grad,nequations,                     &
                      & hatc_elem(:,:,i_node),                       &
                      & dhatcdu_gradwj_elem(:,:,i_node))

            ! Matrix dhat{[C]}/dU * grad(W)_j at the i_node
            dfdu_a_elem(:,:,ka_penI_proc(cnt1)) = dfdu_a_elem(:,:,ka_penI_proc(cnt1))                &
             & - dhatcdu_gradwj_elem(:,:,i_node) / Jx_r(i_node,elem_id) * dt_akk_pinv_l01

            ! Matrix \hat{[C]} * grad_j * dWdU
            do j = iagrad(i_node), iagrad(i_node+1)-1

              jnode = jagrad(jdir,j)
              cnt2 = cnt2 + 1
              
              dfdu_a_elem(:,:,ka_penV_proc(cnt2)) = dfdu_a_elem(:,:,ka_penV_proc(cnt2))              &
                & - matmul(hatc_elem(:,:,i_node),dWdU(vg(:,jnode,elem_id),nequations)) * dagrad(jdir,j) &
                & / Jx_r(i_node,elem_id) * dt_akk_pinv_l01

            end do

          end do

          ! Update cnt1 counter
          cnt1 = cnt1 + 1

          ! Off element contributions
          ! -------------------------
          call primitivevariables(ughst(:,g_node),v_ghost,nequations)
          call entropyvariables(v_ghost,w_ghost,nequations)

          ! Compute point-wise gradient

          ! Initialize ghost computational gradient to zero
          grad_w_ghost_node = 0.0_wp
          
          ! Get the volumetric index of the adjoining node
          k_node = efn2efn(1,j_node,elem_id)

          ! Loop over number of dependent elements in gradient
          do l = iagrad(k_node), iagrad(k_node+1)-1
            ! Loop over dimensions
            do jdir = 1,ndim
              ! Column/node from gradient operator in CSR format in
              ! the jdir-direction corresponding to the coefficient 
              ! dagrad(jdir,i)
              jnode = jagrad(jdir,l) + i_loc*nodesperelem
              ! Update gradient using coefficient and entropy variables at 
              ! appropriate node
              grad_w_ghost_node(:,jdir) = grad_w_ghost_node(:,jdir) + &
                & dagrad(jdir,l)*welemghst(:,jnode) 
            end do   
          end do

          !if (elem_id == 4) then
          !  if (i_node == 2) then
          !    write(70,*) 'i_node', i_node
          !    write(70,*) 'k_node', k_node
          !    write(70,*) 'grad_w_ghost_node(:,jdir)', grad_w_ghost_node(:,1)
          !    write(70,*) 'grad_w_ghost_node(:,jdir)', grad_w_ghost_node(:,2)
          !  endif
          !endif


          !if (myprocid == 1) then
          !  if (i_node == 4) then
          !    write(70,*) 'elem ID', elem_id
          !    write(70,*) 'i_node', i_node
          !    write(70,*) 'k_node', k_node
          !    do jdir = 1,ndim
          !      write(70,*) 'jdir', jdir
          !      write(70,*) 'grad_w', grad_w_ghost_node(:,jdir)
          !    enddo
          !    stop
          !  endif
          !endif

          ! Matrix dhat{[C]}/dU * grad(W)_j at the i_node
          do jdir = 1,ndim

            pos = i + i_loc*nodesperface
            n_grad  = r_x_ghst(jdir,:,pos)

            !if (myprocid == 0) then
            !  if (elem_id == 2) then
            !    if (i_node == 2) then
            !      write(70,*) 'elem ID', elem_id
            !      write(70,*) 'i_node', i_node
            !      write(70,*) 'k_node', k_node
            !      write(70,*) 'jdir', jdir
            !      write(70,*) 'r_x', r_x_ghst(jdir,:,pos)
            !    endif
            !  endif
            !endif

            call hatc_dhatcdu_gradwj_node(v_ghost,       &
                      & grad_w_ghost_node(:,jdir),   &
                      & n_div,n_grad,nequations,                     &
                      & hatc_elem(:,:,k_node),                       &
                      & dhatcdu_gradwj_elem(:,:,k_node))

            !if (myprocid == 1) then
            !  if (i_node == 4) then
            !    write(70,*) 'elem ID', elem_id
            !    write(70,*) 'i_node', i_node
            !    write(70,*) 'k_node', k_node
            !    write(70,*) 'jdir', jdir
            !    write(70,*) 'hatc_elem',hatc_elem(:,:,k_node)
            !    write(70,*) 'hatc_elem',dhatcdu_gradwj_elem(:,:,k_node)
            !endif
          !endif

            dfdu_a_elem(:,:,ka_penI_proc(cnt1)) = dfdu_a_elem(:,:,ka_penI_proc(cnt1)) &
              & + dhatcdu_gradwj_elem(:,:,k_node) / Jx_r(i_node,elem_id) * dt_akk_pinv_l01
            
            !if (myprocid == 0) then
            !  if (i_node == 1) then
            !    write(70,*) 'elem ID', elem_id
            !    write(70,*) 'i_node', i_node
            !    write(70,*) 'k_node', k_node
            !    write(70,*) 'jdir', jdir
            !    write(70,*) 'dfdu_a_elem(:,:,ka_penI_proc(cnt1))', dfdu_a_elem(:,:,ka_penI_proc(cnt1))
            !  endif
            !endif

            ! Matrix \hat{[C]} * grad_j * dWdU
            do j = iagrad(k_node), iagrad(k_node+1)-1

              jnode = jagrad(jdir,j) + i_loc*nodesperelem
              cnt2 = cnt2 + 1

              dfdu_a_elem(:,:,ka_penV_proc(cnt2)) = dfdu_a_elem(:,:,ka_penV_proc(cnt2)) &
                & + matmul(hatc_elem(:,:,k_node),dWdU(velemghst(:,jnode),nequations)) * dagrad(jdir,j) &
                & / Jx_r(i_node,elem_id) * dt_akk_pinv_l01
            end do

            !if (myprocid == 0) then
            !  if (elem_id == 2) then
            !    if (i_node == 2) then                
            !      if(jdir == ndim) then
            !        stop
            !      endif
            !    endif
            !  endif
            !endif


          end do

        enddo

        ! Update the stack coefficient position
        i_loc = i_loc + 1

      else

        ! Loop over each node on the face
        do i = 1, nodesperface

          ! Index in facial ordering
          j_node = nodesperface*(i_face - 1) + i
          
          ! Volumetric node index corresponding to facial node index
          i_node = ifacenodes(j_node)
        
          ! Volumetric index of partner node
          k_node = efn2efn(1,j_node,elem_id)

          ! Element index of partner node
          k_elem = efn2efn(2,j_node,elem_id)
        
          ! Vector needed for reusing on element subroutines
          n_div  = Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          ! Update cnt1 counter
          cnt1 = cnt1 + 1

          ! ON element contributions
          ! ------------------------
          do jdir = 1,ndim

            ! Vector for the computation of the gradient
            n_grad = r_x(jdir,:,i_node,elem_id)

            call hatc_dhatcdu_gradwj_node(vg(:,i_node,elem_id),      &
                      & grad_w_jacobian(:,jdir,i_node,elem_id),      &
                      & n_div,n_grad,nequations,                     &
                      & hatc_elem(:,:,i_node),                       &
                      & dhatcdu_gradwj_elem(:,:,i_node))

            ! Matrix dhat{[C]}/dU * grad(W)_j at the i_node
            dfdu_a_elem(:,:,ka_penI_proc(cnt1)) = dfdu_a_elem(:,:,ka_penI_proc(cnt1))                &
             & - dhatcdu_gradwj_elem(:,:,i_node) / Jx_r(i_node,elem_id) * dt_akk_pinv_l01 

            ! Matrix \hat{[C]} * grad_j * dWdU
            do j = iagrad(i_node), iagrad(i_node+1)-1

              jnode = jagrad(jdir,j)
              cnt2 = cnt2 + 1

              dfdu_a_elem(:,:,ka_penV_proc(cnt2)) = dfdu_a_elem(:,:,ka_penV_proc(cnt2))              &
             & - matmul(hatc_elem(:,:,i_node),dWdU(vg(:,jnode,elem_id),nequations)) * dagrad(jdir,j) &
             & / Jx_r(i_node,elem_id) * dt_akk_pinv_l01
             
            end do

          end do

          ! Update cnt1 counter
          cnt1 = cnt1 + 1

          ! Off element contributions
          ! -------------------------
          do jdir = 1,ndim

            n_grad = r_x(jdir,:,k_node,k_elem)

            !if (elem_id == 2) then
            !  if (i_node == 4) then
            !    if (k_elem == 3) then
            !      write(70,*) 'i_node', i_node
            !      write(70,*) 'k_node', k_node
            !      write(70,*) 'jdir', jdir
            !      write(70,*) 'r_x', r_x(jdir,:,k_node,k_elem)
            !    endif
            !  endif
            !endif

            !write(70,*) 'elem ID', elem_id
            !write(70,*) 'i_node', i_node
            !write(70,*) 'k_node', k_node
            !write(70,*) 'jdir', jdir
            !write(70,*) 'r_x', r_x(jdir,:,k_node,k_elem)



            !if (elem_id == 5) then
            !  if (i_node == 2) then
            !    write(70,*) 'proc', myprocid
            !    write(70,*) 'i_node', i_node
            !    write(70,*) 'k_node', k_node
            !    write(70,*) 'grad_w_jacobian(:,jdir,k_node,k_elem)', grad_w_jacobian(:,jdir,k_node,k_elem)
            !  endif
            !endif


          
            call hatc_dhatcdu_gradwj_node(vg(:,k_node,k_elem),       &
                      & grad_w_jacobian(:,jdir,k_node,k_elem),       &
                      & n_div,n_grad,nequations,                     &
                      & hatc_elem(:,:,k_node),                       &
                      & dhatcdu_gradwj_elem(:,:,k_node))

            !write(70,*) 'elem ID', elem_id
            !write(70,*) 'i_node', i_node
            !write(70,*) 'k_node', k_node
            !write(70,*) 'jdir', jdir
            !write(70,*) 'hatc_elem',hatc_elem(:,:,k_node)
            !write(70,*) 'hatc_elem',dhatcdu_gradwj_elem(:,:,k_node)

            ! Matrix dhat{[C]}/dU * grad(W)_j at the i_node
            dfdu_a_elem(:,:,ka_penI_proc(cnt1)) = dfdu_a_elem(:,:,ka_penI_proc(cnt1))               &
              & + dhatcdu_gradwj_elem(:,:,k_node) / Jx_r(i_node,elem_id) * dt_akk_pinv_l01

            !write(70,*) 'elem ID', elem_id
            !write(70,*) 'i_node', i_node
            !write(70,*) 'k_node', k_node
            !write(70,*) 'jdir', jdir
            !write(70,*) 'dfdu_a_elem(:,:,ka_penI_proc(cnt1))', dfdu_a_elem(:,:,ka_penI_proc(cnt1))

            ! Matrix \hat{[C]} * grad_j * dWdU
            do j = iagrad(k_node), iagrad(k_node+1)-1

              jnode = jagrad(jdir,j)
              cnt2 = cnt2 + 1

              dfdu_a_elem(:,:,ka_penV_proc(cnt2)) = dfdu_a_elem(:,:,ka_penV_proc(cnt2))             &
             & + matmul(hatc_elem(:,:,k_node),dWdU(vg(:,jnode,k_elem),nequations)) * dagrad(jdir,j) &
             & / Jx_r(i_node,elem_id) * dt_akk_pinv_l01

            end do

            !if (elem_id == 2) then
            !  if (i_node == 4) then
            !    if (k_elem == 3) then
            !      if(jdir == ndim) then
            !        stop
            !      endif
            !    endif
            !  endif
            !endif

          end do

        end do

      endif

    enddo

    return
  end subroutine compute_viscous_residual_sat_jacobian_element
  
  !============================================================================
  
  !============================================================================
  ! jacobian_viscous_bc_sat_node - Computes the contributions to the viscous 
  ! residual Jacobian matrix of the BC penalty term for one node on the element 
  ! face. 
  !
  ! Input parameters:
  ! elem_id - global element ID.
  ! i_face - face ID.
  ! i_node - left ID node.
  ! k_node - right ID node.
  ! n_eq - number of equations.
  !
  ! Output parameters:
  ! jacobian_inviscid_penalty - Jacobian matrices of the inviscid penalty.

  function jacobian_viscous_bc_sat_node(elem_id,i_face,i_node,n_v,n_eq)
 
    ! Load modules
    use nsereferencevariables
    use referencevariables
    use variables, only : ug, phig, xg, ef2e, vg
    use controlvariables, only : timelocal
    use collocationvariables, only : pinv
    use navierstokes
    use CSRlocalvariables
    
    ! Nothing is implicitly defined
    implicit none

    integer,                intent(in) :: elem_id, i_face, i_node 
    real(wp), dimension(3), intent(in) :: n_v
    integer,                intent(in) :: n_eq

    real(wp), dimension(n_eq,n_eq)   :: jacobian_viscous_bc_sat_node

    ! Here, according to the BC we call the function tha implements the
    ! analytical contribution of the normal viscous flux.
    select case(ef2e(1,i_face,elem_id))
      case default
        hatc_elem = 0.0_wp
        jacobian_viscous_bc_sat_node(:,:) = 0.0_wp
    end select

    return
  end function jacobian_viscous_bc_sat_node

  !============================================================================
  
  !============================================================================
  ! jacobian_inviscid_bc_sat_node - Computes the contributions to the inviscid 
  ! residual Jacobian matrix of the BC penalty term for one node on the element 
  ! face. 
  !
  ! Input parameters:
  ! elem_id - global element ID.
  ! i_face - face ID.
  ! i_node - left ID node.
  ! k_node - right ID node.
  ! n_eq - number of equations.
  !
  ! Output parameters:
  ! jacobian_inviscid_penalty - Jacobian matrices of the inviscid penalty.

  function jacobian_inviscid_bc_sat_node(elem_id,i_face,i_node,n_v,n_eq)
 
    ! Load modules
    use nsereferencevariables
    use referencevariables
    use variables, only : ug, phig, xg, ef2e, vg
    use controlvariables, only : timelocal
    use collocationvariables, only : pinv
    use navierstokes
    
    ! Nothing is implicitly defined
    implicit none

    integer,                intent(in) :: elem_id, i_face, i_node 
    real(wp), dimension(3), intent(in) :: n_v
    integer,                intent(in) :: n_eq

    real(wp), dimension(n_eq,n_eq,2)   :: jacobian_inviscid_bc_sat_node

    complex(wp), dimension(n_eq)       :: u_i_node!, u_star
    complex(wp), dimension(n_eq)       :: v_i_node, v_star
    complex(wp), dimension(n_eq)       :: w_i_node, w_star
    
    complex(wp), dimension(n_eq,n_eq)  :: s_inv, s_mat
    complex(wp), dimension(n_eq)       :: v_avg, ev, f_star

    real(wp), dimension(n_eq,n_eq)     :: dfndu
 
    real(wp), dimension(n_eq)          :: zero_v, ev_abs
    real(wp)                           :: ev_max, s_fix

    integer             :: comp

    real(wp), parameter :: eps = 1e-15_wp
    real(wp), parameter :: eps_inv = 1.0_wp/eps

    ! Parameter for the entropy fix
    s_fix = 0.0001_wp

    ! Initialize Jacobian contributions
    jacobian_inviscid_bc_sat_node = 0.0_wp

    ! Vector of zeros for the initialization of the imaginary parts
    zero_v = 0.0_wp

    ! Initialize u_i_node
    u_i_node = cmplx(real(ug(:,i_node,elem_id),wp),zero_v,wp) 

    ! Compute analitically dfn/du
    dfndu = inviscid_flux_jacobian_node(vg(:,i_node,elem_id),n_v,n_eq)

    do comp = 1, n_eq
      ! Perturb the imaginary part of one component of u_i_node
      u_i_node(comp) = cmplx(real(ug(comp,i_node,elem_id),wp),eps,wp)

      ! Compute the perturbed primitive variables at the i_node
      v_i_node = u_to_v_complex(u_i_node,n_eq)
      
      ! Compute the perturbed entropy variables at the i_node
      w_i_node = v_to_w_complex(v_i_node,n_eq)

      ! Set v_star
      v_star = v_i_node

      select case(ef2e(1,i_face,elem_id))
        case(-03)
          call inviscid_wall_complex(v_star,n_v,n_eq)
        case(-04)
          !boundary_condition_complex => no_slip_wall_isothermal_complex
        case(-06)
          !boundary_condition_complex => subsonic_outflow_complex
        case(-08)
          !boundary_condition_complex => uniform_free_stream_complex
        case(-09)
          !boundary_condition_complex => uniform_free_stream_complex
        case(-11)
          !boundary_condition_complex => uniform_free_stream_complex
        case(-12)  !  BC Dirichlet handles all the exact solution data
          if (InitialCondition == 'ExactSolutionIsentropicVortex') then
            call isentropic_vortex_full_complex(v_star,n_v, &
              & xg(:,i_node,elem_id),timelocal,n_eq)
          elseif(InitialCondition  == 'ExactSolutionViscousShock') then
            call viscous_shock_full_complex(v_star,n_v,xg(:,i_node,elem_id), &
              & timelocal,n_eq)
          endif
        case(-13)
          !boundary_condition_complex => symmetry_plane_complex
        case default
          write(*,*) 'Something wrong with the computation of the penalty term &
            & arising from the inviscid BC.'
          write(*,*) 'Stopping...'
          stop
      end select

      ! Compute w_star 
      w_star = v_to_w_complex(v_star,n_eq)

      ! Compute Roe-averaged variables
      call roe_average_complex(v_i_node,v_star,v_avg,n_eq)
      
      ! Compute the characteristic decomposition using the Roe-averaged state
      call v_to_characteristic_decomposition_complex(v_avg,n_eq,s_inv,s_mat, &
        & ev,n_v)
     
      ! Get the maximum eigenvalue modules
      ev_max = maxval(abs(ev(:)))

      ! Apply entropy fix
      ev_abs(:) = sqrt(conjg(ev(:))*ev(:) + s_fix*ev_max)

      ! Compute the entropy consitent flux
      f_star = v_to_entropy_consistent_flux_complex(v_i_node,v_star,n_v,n_eq)

      ! Apply the entropy consistent flux
      f_star = f_star + half * matmul(s_mat,ev_abs*matmul(transpose(s_mat), &
        & w_i_node - w_star))

      ! Set one column of the Jacobian matrix
      jacobian_inviscid_bc_sat_node(:,comp,1) = aimag((-f_star)*pinv(1))*eps_inv
      
      ! Restore the perturbed component of v_i_node
      u_i_node(comp) = cmplx(real(ug(comp,i_node,elem_id),wp),0.0_wp,wp)

      ! Compute the perturbed primitive variables at the i_node
      v_i_node = u_to_v_complex(u_i_node,n_eq)

      ! Restore the perturbed entropy variables at the i_node
      w_i_node = v_to_w_complex(v_i_node,n_eq)

    enddo
    
    ! Add analytical contribution of the "standard" normal flux
    jacobian_inviscid_bc_sat_node(:,:,1) = jacobian_inviscid_bc_sat_node(:,:,1) + &
      & dfndu*pinv(1)

    return
  end function jacobian_inviscid_bc_sat_node

  !============================================================================
  
  !============================================================================
  ! viscous_shock_full_complex - Sets the ghost values using the analytical
  ! expression for a viscous shock. Complex variables are used.

  subroutine viscous_shock_full_complex(v_in_out,n_v,x_in,t_in,n_eq_in)
    
    ! Load modules
    use nsereferencevariables
    use navierstokes, only : rhalf

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq_in
    complex(wp), intent(inout) :: v_in_out(n_eq_in)
    real(wp), intent(in) :: n_v(3)
    real(wp), intent(in) :: x_in(3), t_in
    
    real(wp) :: Uinf
    real(wp) :: x0,y0,r0
    real(wp) :: f, xp
    real(wp) :: alph, vf, mdot, Mach, wave, M2, uhat, Rloc, gloc
    real(wp) :: cc, ss
    real(wp) :: theta
    real(wp) :: mu, kappa
    
    integer :: ns

    ns = n_eq_in - 4
    theta = uniformFreeStreamAOA*pi/180._wp
    r0 = membranelocation
    y0 = zero
    x0 = 0.0_wp
    Uinf = 1.0_wp
    M2 = Mach0*Mach0 !Uinf*Uinf
    mach = sqrt(M2)
    mdot = abs(Uinf)

    cc = cos(theta)
    ss = sin(theta)

    xp = x_in(1)*cc + x_in(2)*ss
    
    ! Initialize v_in_out (input data is trashed because the analytic BC is 
    ! imposed)
    v_in_out = cmplx(0.0_wp,0.0_wp,wp)

    ! set mixture variable
    Rloc = 1.0_wp
    gloc = gamma0
    
    wave = referencewavespeed

    mu = 1.0_wp
    kappa = 1.0_wp

    vf = (gm1+two/M2)/(gamma0+one)
    
    alph = (four/three*mu/Re0)/mdot*(two*gamma0)/(gamma0+one)
    
    call rhalf(xp-r0-wave*t_in,f,alph,vf)
    
    uhat = f*Uinf
    
    ! Density
    v_in_out(1) = abs(mdot)/f
    
    ! Velocity
    v_in_out(2) = (uhat + wave)*cc
    v_in_out(3) = (uhat + wave)*ss

    ! Temperature
    v_in_out(5) = (Uinf*Uinf + gm1M2*half*(Uinf*Uinf-uhat*uhat))/Rloc
    
    return
  end subroutine viscous_shock_full_complex

  !============================================================================
 
  !============================================================================
  ! isentropic_vortex_full_complex - Sets the ghost values using the analytical
  ! expression for a isontropic vortex. Complex variables are used.

  pure subroutine isentropic_vortex_full_complex(v_in_out,n_v_in,x_in,t_in, &
      & n_eq_in)
    
    ! Load modules
    use nsereferencevariables
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq_in
    complex(wp), intent(inout) :: v_in_out(n_eq_in)
    real(wp), intent(in) :: n_v_in(3)
    real(wp), intent(in) :: x_in(3), t_in
    
    real(wp) :: epsvortex
    real(wp) :: u_inf
    real(wp) :: x0,y0
    real(wp) :: f
    real(wp) :: alpha, rin2
   
    ! Center of the vortex at t = t_0
    y0 = 0.0_wp
    x0 = 0.0_wp
    
    ! Module of the far field velocity
    u_inf = referencewavespeed

    ! Vortex width
    epsvortex = 5.0_wp
    
    ! Angle of attach of the free-stream that carries the vortex 
    alpha = uniformFreeStreamAOA*pi/180._wp

    ! Vortex radius squared
    rin2 = ((x_in(1)-x0)-u_inf*t_in*cos(alpha))**2 + ((x_in(2)-y0)- &
      & u_inf*t_in*sin(alpha))**2

    f = 1.0_wp - rin2
    
    ! Initialize v_in_out (input data is trashed because the analytic BC is 
    ! imposed)
    v_in_out = cmplx(0.0_wp,0.0_wp,wp)
    
    ! Density
    v_in_out(1) = (one-gm1M2*epsvortex*epsvortex/(eight*pi*pi)*exp(f))** &
      & (one/gm1)

    ! Velocity vector (2D)
    v_in_out(2) = u_inf*cos(alpha) + epsvortex/(two*pi)*((x_in(2)-y0)- &
      & u_inf*t_in*sin(alpha))*exp(f*half)
    v_in_out(3) = u_inf*sin(alpha) - epsvortex/(two*pi)*((x_in(1)-x0)- &
      & u_inf*t_in*cos(alpha))*exp(f*half)
    
    ! Temperature
    v_in_out(5) = (one-gm1M2*epsvortex*epsvortex/ &
      & (eight*pi*pi)*exp(f))

    return
  end subroutine isentropic_vortex_full_complex

  !============================================================================

 !============================================================================
  ! inviscid_wall_complex - Computes the ghost cell values for the inviscid wall
  ! boundary condition. Complex variables are used.
  !
  ! Input parameters:
  ! v_in_out - primitive variables at the i_node (intent(inout)).
  ! n_v - normal vector.
  ! n_eq - number of equations.
  !
  ! Output parameters:
  ! v_in_out - primitive variables at the ghost node.

  pure subroutine inviscid_wall_complex(v_in_out,n_v_in,n_eq_in)

    ! Load modules
    use nsereferencevariables

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq_in
    complex(wp), intent(inout) :: v_in_out(n_eq_in)
    real(wp), intent(in) :: n_v_in(3)

    real(wp), parameter  :: Mirror = -1.0_wp
    complex(wp), dimension(3)   :: v_normal, v_tangent, n_v_magnitude

    ! Magnitude of the normal vector
    n_v_magnitude = n_v_in/Magnitude(n_v_in)

    ! Normal velocity
    v_normal = (v_in_out(2)*n_v_magnitude(1) + v_in_out(3)*n_v_magnitude(2) + &
      & v_in_out(4)*n_v_magnitude(3))*n_v_magnitude
    
    ! Tangential velocity
    v_tangent = v_in_out(2:4) - v_normal

    ! Density
    v_in_out(1)   = v_in_out(1)
    
    ! Velocity vector
    v_in_out(2:4) = v_tangent + Mirror*v_normal
    
    ! Temperature
    v_in_out(5)   = v_in_out(5)

    return
  end subroutine inviscid_wall_complex

  !============================================================================
 
  !============================================================================
  ! jacobian_inviscid_sat_node - Computes the contributions to the residual 
  ! Jacobian matrix of the inviscid penalty term for one node on the element 
  ! face. Two contributions are computed: One for the i_node and one for the 
  ! k_node which lives on the same position but belongs to the adjacent element. 
  !
  ! Input parameters:
  ! elem_id - global element ID.
  ! i_node - left ID node.
  ! k_node - right ID node.
  ! n_eq - number of equations.
  !
  ! Output parameters:
  ! jacobian_inviscid_penalty - Jacobian matrices of the inviscid penalty.

  function jacobian_inviscid_sat_node(u_l,u_r,n_v,n_eq)
 
    ! Load modules
    use nsereferencevariables
    use referencevariables
    use variables, only : ug, vg
    use collocationvariables, only: pinv
    use navierstokes, only : primitivevariables
    
    ! Nothing is implicitly defined
    implicit none
    
    !integer, intent(in) :: i_elem, k_elem
    !integer, intent(in) :: i_node, k_node
    real(wp), dimension(3), intent(in) :: n_v
    integer, intent(in) :: n_eq
    real(wp), dimension(n_eq), intent(in) :: u_l, u_r


    real(wp), dimension(n_eq,n_eq,2) :: jacobian_inviscid_sat_node

    real(wp), parameter :: eps = 1e-15_wp
    real(wp), parameter :: eps_inv = 1.0_wp/eps
    
    complex(wp), dimension(n_eq) :: u_i_node, u_k_node
    complex(wp), dimension(n_eq) :: v_i_node, v_k_node
    complex(wp), dimension(n_eq) :: w_i_node, w_k_node
    
    integer :: comp
    real(wp), dimension(n_eq) :: zero_v

    complex(wp), dimension(n_eq) :: v_avg
    
    complex(wp), dimension(n_eq,n_eq) :: s_inv, s_mat
    complex(wp), dimension(n_eq) :: ev

    real(wp) :: ev_max
    real(wp), dimension(n_eq) :: ev_abs
    real(wp) :: s_fix

    complex(wp), dimension(n_eq) :: f_star

    real(wp), dimension(n_eq,n_eq) :: dfndu
    real(wp), dimension(n_eq) :: v_l

    ! Parameter for the entropy fix
    s_fix = 0.0001_wp

    ! Initialize Jacobian contributions
    jacobian_inviscid_sat_node = 0.0_wp

    ! Vector of zero for the initialization
    zero_v(:) = 0.0_wp
    
    u_i_node = cmplx(real(u_l,wp), zero_v,wp) 

    u_k_node = cmplx(real(u_r,wp), zero_v,wp)
    v_k_node = u_to_v_complex(u_k_node,n_eq)
    w_k_node = v_to_w_complex(v_k_node,n_eq)

    call primitivevariables(u_l,v_l,n_eq)
    
    ! Compute analitically dfn/du
    dfndu = inviscid_flux_jacobian_node(v_l,n_v,n_eq)

    do comp = 1, n_eq
      ! Perturb the imaginary part of one component of u_i_node
      u_i_node(comp) = cmplx(real(u_l(comp),wp),eps,wp)

      ! Compute the perturbed primitive variables at the i_node
      v_i_node(:) = u_to_v_complex(u_i_node,n_eq)

      ! Compute the perturbed entropy variables at the i_node
      w_i_node(:) = v_to_w_complex(v_i_node,n_eq)
      
      ! Compute Roe-averaged variables
      call roe_average_complex(v_i_node,v_k_node,v_avg,n_eq)
      
      ! Compute the characteristic decomposition using the Roe-averaged state
      call v_to_characteristic_decomposition_complex(v_avg,n_eq,s_inv,s_mat, &
        & ev,n_v)
     
      ! Get the maximum eigenvalue modules
      ev_max = maxval(cdabs(ev(:)))

      ! Apply entropy fix
      ev_abs(:) = sqrt(conjg(ev(:))*ev(:) + s_fix*ev_max)

      ! Compute the entropy consitent flux
      f_star = v_to_entropy_consistent_flux_complex(v_i_node,v_k_node,n_v,n_eq)

      ! Apply the entropy consistent flux
      f_star = f_star + half * matmul(s_mat,ev_abs*matmul(transpose(s_mat), &
        & w_i_node - w_k_node))

      ! Set one column of the Jacobian matrix
      jacobian_inviscid_sat_node(:,comp,1) = aimag(f_star)*eps_inv*pinv(1)

      ! Restore the perturbed component of v_i_node
      u_i_node(comp) = cmplx(real(u_l(comp),wp),0.0_wp,wp)

    enddo

    jacobian_inviscid_sat_node(:,:,1) = jacobian_inviscid_sat_node(:,:,1) - dfndu*pinv(1)
    
    u_k_node = cmplx(real(u_r,wp), zero_v,wp) 

    u_i_node = cmplx(real(u_l,wp), zero_v,wp)
    v_i_node = u_to_v_complex(u_i_node,n_eq)
    w_i_node = v_to_w_complex(v_i_node,n_eq)

    do comp = 1, n_eq
      ! Perturb the imaginary part of one component of u_k_node
      u_k_node(comp) = cmplx(real(u_r(comp),wp),eps,wp)

      ! Compute the perturbed primitive variables at the k_node
      v_k_node(:) = u_to_v_complex(u_k_node,n_eq)
      
      ! Compute the perturbed entropy variables at the k_node
      w_k_node(:) = v_to_w_complex(v_k_node,n_eq)
      
      ! Compute Roe-averaged variables
      call roe_average_complex(v_i_node,v_k_node,v_avg,n_eq)
      
      ! Compute the characteristic decomposition using the Roe-averaged state
      call v_to_characteristic_decomposition_complex(v_avg,n_eq,s_inv,s_mat, &
        & ev,n_v)
     
      ! Get the maximum eigenvalue modules
      ev_max = maxval(cdabs(ev(:)))

      ! Apply entropy fix
      ev_abs(:) = sqrt(conjg(ev(:))*ev(:) + s_fix*ev_max)

      ! Compute the entropy consitent flux
      f_star = v_to_entropy_consistent_flux_complex(v_i_node,v_k_node,n_v,n_eq)

      ! Apply the entropy consistent flux
      f_star = f_star + half * matmul(s_mat,ev_abs*matmul(transpose(s_mat), &
        & w_i_node - w_k_node))

      ! Set one column of the Jacobian matrix
      jacobian_inviscid_sat_node(:,comp,2) = aimag(f_star)*eps_inv*pinv(1)

      ! Restore the perturbed component of v_i_node
      u_k_node(comp) = cmplx(real(u_r(comp),wp),0.0_wp,wp)

    enddo

    return
  end function jacobian_inviscid_sat_node

  !============================================================================

  !============================================================================
  ! roe_average_complex - Computes the Roe averaged given the primitive
  ! variables on the left and right state. Complex data type are used.
  !
  ! Input parameters:
  ! v_l - left state (primitive variables).
  ! v_r - right state (primitive variables).
  ! n_eq - number of equations.
  !
  ! Output parameters:
  ! v_avg - Roe averaged variable (primitive variables).

  pure subroutine roe_average_complex(v_l,v_r,v_avg,n_eq)
    
    ! Load modules
    use nsereferencevariables, only : gm1M2

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq
    complex(wp), dimension(n_eq), intent(in) :: v_l, v_r
    complex(wp), dimension(n_eq), intent(out) :: v_avg
    complex(wp) :: roe, roe_p1_inv
    complex(wp) :: ht_l, ht_r, ht_avg

    ! Initialize the Roe average to zero
    v_avg = cmplx(real(0.0_wp,wp),0.0_wp,wp)
    
    ! Density ratio
    roe = sqrt(v_r(1)/v_l(1))
    
    ! Save for efficiency
    roe_p1_inv = 1.0_wp/(roe + 1.0_wp)

    ! Roe averaged density
    v_avg(1) = roe*v_l(1)

    ! Roe averaged velocity (2:4) 
    v_avg(2:4) = (roe*v_r(2:4) + v_l(2:4))*roe_p1_inv

    ! Roe averaged total enthalpy
    ht_l = v_l(5) + gm1M2*0.5_wp*(v_l(2)*v_l(2) + v_l(3)*v_l(3) + &
      & v_l(4)*v_l(4))
    ht_r = v_r(5) + gm1M2*0.5_wp*(v_r(2)*v_r(2) + v_r(3)*v_r(3) + &
      & v_r(4)*v_r(4))
    ht_avg = (roe*ht_r + ht_l)*roe_p1_inv
    
    ! Roe averaged temperature
    v_avg(5) = ht_avg - gm1M2*0.5_wp*(v_avg(2)*v_avg(2) + v_avg(3)*v_avg(3) + &
      & v_avg(4)*v_avg(4))

    return
  end subroutine roe_average_complex

  !============================================================================

  !============================================================================
  ! v_to_characteristic_decomposition_complex - Calculates the left and the 
  ! right eigenvetor matrices and the eigenvalues of the flux Jacobian in the 
  ! normal direction given a set of primitive variables. The resulting right 
  ! eigenvectors have the property that [R][R]^T = du/dw. Complex data types are
  ! used.
  !                  
  ! Input parameters:
  ! v_in - primitive variables.
  ! n_eq - number of equations.
  ! Jx - contravariant normal vector.
  !
  ! Output parameters:
  ! e_l - matrix of the left eigenvectors of the flux Jacobian.
  ! e_r - matrix of the right eigenvectors of the flux Jacobian.
  ! ev - eiegenvalues of the flux Jacobian of the flux Jacobian.

  pure subroutine v_to_characteristic_decomposition_complex(v_in,n_eq,e_l,e_r, &
      & ev,n_v)
    
    ! Load modules
    use nsereferenceVariables, only: gM2, gm1M2, gamma0, gm1, gm1og

    ! Nothing is implicitly defined
    implicit none

    integer,intent(in) :: n_eq
    complex(wp), dimension(n_eq), intent(in) :: v_in
    complex(wp), dimension(n_eq,n_eq), intent(out) :: e_l, e_r
    complex(wp), dimension(n_eq), intent(out) :: ev
    real(wp), intent(in) :: n_v(3)

    complex(wp) :: rho_inv, T_inv
    complex(wp) :: cav, cav_inv, cav2, cav2_inv
    complex(wp) :: u_n
    complex(wp) :: mat_tmp(n_eq,n_eq), tmp, tmp2 
    complex(wp) :: g1, g2, g1_inv, g2_inv
    
    real(wp) :: lgm1, lgm1og, lg_inv
    real(wp) :: l_tmp, l_tmp_inv
    real(wp) :: sqrt2, sqrt2_inv

    ! Initialize to zero the left and the right eigenvector matrices and the
    ! eigenvalues vector
    e_l = cmplx(real(0.0_wp,wp),0.0_wp,wp)
    e_r = cmplx(real(0.0_wp,wp),0.0_wp,wp)
    ev  = cmplx(real(0.0_wp,wp),0.0_wp,wp)

    ! Magnitude of the normal vector  
    l_tmp = sqrt(dot_product(n_v,n_v))

    ! Inverse magnitude of normal vector
    l_tmp_inv = 1.0_wp/l_tmp

    ! Inverse of density
    rho_inv = 1.0_wp/v_in(1)

    ! Inverse of temperature
    T_inv = 1.0_wp/v_in(5)

    ! Convenience parameters
    sqrt2 = sqrt(2.0_wp)
    sqrt2_inv = 1.0_wp/sqrt2
    lg_inv = 1.0_wp/gamma0
    lgm1 = gm1
    lgm1og = gm1og

    ! Speed of sound
    cav2 = gamma0*v_in(5)/gM2
    cav = sqrt(cav2)

    ! Inverse speed of sound
    cav2_inv = 1.0_wp/cav2
    cav_inv = 1.0_wp/cav

    ! Normal velocity vector
    u_n = v_in(2)*n_v(1) + v_in(3)*n_v(2) + v_in(4)*n_v(3)

    ! Other convenience parameters
    g1 = sqrt(0.5_wp*v_in(5)*v_in(1)*cav2_inv/gm1M2)
    g1_inv = 1.0_wp/g1
    g2 = sqrt(1.0_wp/gm1M2)
    g2_inv = 1.0_wp/g2

    ! calculate eigenvalues
    ev = u_n
    ev(1) = u_n - cav*l_tmp
    ev(2) = u_n + cav*l_tmp

    ! Calculate the right eigenvectors of the primtive flux jacobian 
    ! dv/du*df/du*du/dv
    ! First eigenvector
    e_r(1,1) = g1
    e_r(2:4,1) = -g1*l_tmp_inv*n_v(1:3)*cav*rho_inv
    e_r(5,1) = g1*lgm1*v_in(5)*rho_inv

    ! Second eigenvector
    e_r(1,2)   = g1
    e_r(2:4,2) = g1 * l_tmp_inv * n_v(1:3) * cav * rho_inv
    e_r(5,2)   = g1 * lgm1      * v_in(5)        * rho_inv

    tmp2 = sqrt(lgm1*v_in(5)*v_in(1)*cav2_inv)
    tmp = sqrt(v_in(5)*rho_inv)

    ! Third eigenvector
    e_r(1,3) = g2*n_v(1)*l_tmp_inv*tmp2
    e_r(3,3) = -g2*n_v(3)*l_tmp_inv*tmp
    e_r(4,3) = g2*n_v(2)*l_tmp_inv*tmp
    e_r(5,3) = -g2*n_v(1)*v_in(5)*rho_inv*l_tmp_inv*tmp2

    ! Fourth eigenvector
    e_r(1,4) = +g2*n_v(2)*l_tmp_inv*tmp2
    e_r(2,4) = +g2*n_v(3)*l_tmp_inv*tmp
    e_r(4,4) = -g2*n_v(1)*l_tmp_inv*tmp
    e_r(5,4) = -g2*n_v(2)*v_in(5)*rho_inv*l_tmp_inv*tmp2

    ! Fifth eigenvector
    e_r(1,5) = g2*n_v(3)*l_tmp_inv*tmp2
    e_r(2,5) = -g2*n_v(2)*l_tmp_inv*tmp
    e_r(3,5) = g2*n_v(1)*l_tmp_inv*tmp
    e_r(5,5) = -g2*n_v(3)*v_in(5)*rho_inv*l_tmp_inv*tmp2

    ! Transform to eigenvectors of df/du by similarity
    ! dv/du*df/du*du/dv = S \Lambda S^{-1} so
    ! df/du = R \Lambda R^{-1} = du/dv*S \Lambda S^{-1}*dv/du
    mat_tmp = v_to_dudv_complex(v_in,n_eq)
    e_r = matmul(mat_tmp,e_r)

    ! Calculate the left eigenvectors of the primtive flux jacobian 
    ! dv/du*df/du*du/dv
    tmp = sqrt(v_in(1)*T_inv)
    tmp2 = 1.0_wp/sqrt(v_in(5)*rho_inv*lgm1*cav2_inv)

    ! First column
    e_l(1:2,1) = g1_inv*0.5_wp*lg_inv
    e_l(3:5,1) = g2_inv*lgm1og*rho_inv*n_v(1:3)*l_tmp_inv*tmp2

    ! Second column
    e_l(1,2) = -0.5_wp*g1_inv*v_in(1)*cav_inv*n_v(1)*l_tmp_inv
    e_l(2,2) =  0.5_wp*g1_inv*v_in(1)*cav_inv*n_v(1)*l_tmp_inv
    e_l(4,2) =  g2_inv*n_v(3)*l_tmp_inv*tmp
    e_l(5,2) = -g2_inv*n_v(2)*l_tmp_inv*tmp

    ! Third column
    e_l(1,3) = -0.5_wp*g1_inv*v_in(1)*cav_inv*n_v(2)*l_tmp_inv
    e_l(2,3) =  0.5_wp*g1_inv*v_in(1)*cav_inv*n_v(2)*l_tmp_inv
    e_l(3,3) = -g2_inv*n_v(3)*l_tmp_inv*tmp
    e_l(5,3) =  g2_inv*n_v(1)*l_tmp_inv*tmp

    ! Fourth column
    e_l(1,4) = -0.5_wp*g1_inv*v_in(1)*cav_inv*n_v(3)*l_tmp_inv
    e_l(2,4) =  0.5_wp*g1_inv*v_in(1)*cav_inv*n_v(3)*l_tmp_inv
    e_l(3,4) =  g2_inv*n_v(2)*l_tmp_inv*tmp
    e_l(4,4) = -g2_inv*n_v(1)*l_tmp_inv*tmp

    ! Fifth column
    e_l(1:2,5) = g1_inv*0.5_wp*lg_inv*v_in(1)*T_inv
    e_l(3:5,5) = -g2_inv*T_inv*lg_inv*n_v(1:3)*l_tmp_inv*tmp2

    ! Transform to eigenvectors of df/du by similarity
    ! dv/du*df/du*du/dv = S \Lambda S^{-1} so
    ! df/du = R \Lambda R^{-1} = du/dv*S \Lambda S^{-1}*dv/du
    mat_tmp = v_to_dvdu_complex(v_in,n_eq)
    e_l = matmul(e_l,mat_tmp)

    return
  end subroutine v_to_characteristic_decomposition_complex

  !============================================================================

  !============================================================================
  ! v_to_dudv_complex - Calculates the Jacobian of the conservative variables 
  ! with respect to the primitive variables using the primitive variables. 
  ! Complex data type are used.
  !
  ! Input parameters:
  ! v_in - input primitive variables.
  ! n_eq - number of equations.
  ! 
  ! Output parameters:
  ! v_to_dudv_complex - Jacobian of the conservative variables with respect to  
  !                     the primitive variables.

  pure function v_to_dudv_complex(v_in,n_eq)
    
    ! Load modules
    use nsereferencevariables, only : gm1M2, gm1og, gm1

    ! Nothing implicitly defined
    implicit none

    integer, intent(in) :: n_eq
    complex(wp), intent(in) :: v_in(n_eq)

    complex(wp) :: v_to_dudv_complex(n_eq,n_eq)
    complex(wp) :: ht

    ! Initialize to zero the Jacobian
    v_to_dudv_complex = cmplx(real(0.0_wp,wp),0.0_wp,wp)

    ! Total enthalpy
    ht = v_in(5) + gm1M2*0.5_wp*(v_in(2)*v_in(2) + v_in(3)*v_in(3) + &
      & v_in(4)*v_in(4))

    ! Continuity
    v_to_dudv_complex(1,1) = cmplx(real(1.0_wp,wp),0.0_wp,wp)
    
    ! Momentum
    v_to_dudv_complex(2:4,1) = v_in(2:4)
    v_to_dudv_complex(2,2) = v_in(1)
    v_to_dudv_complex(3,3) = v_in(1)
    v_to_dudv_complex(4,4) = v_in(1)
    
    ! Energy
    v_to_dudv_complex(5,1) = ht - gm1og*v_in(5)
    v_to_dudv_complex(5,2:4) = gm1M2*v_in(1)*v_in(2:4)
    v_to_dudv_complex(5,5) = gm1og*v_in(1)/gm1

    return
  end function v_to_dudv_complex


  !============================================================================
  ! v_to_dvdu_complex - Calculates the Jacobian of the primitive variables with
  ! respect to the conserved variables. Complex variables are used.
  !
  ! Input parameters:
  ! v_in - input primitive variables.
  ! n_eq - number of equations.
  ! 
  ! Output parameters:
  ! v_to_dvdu_complex - Jacobian of the primitive variables with respect to the 
  !                     conserved variables using complex variables.

  pure function v_to_dvdu_complex(v_in,n_eq)
    
    ! Load modules
    use nsereferencevariables, only : gm1M2, gM2, gm1, gamma0

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq
    complex(wp), intent(in) :: v_in(n_eq)
    complex(wp) :: v_to_dvdu_complex(n_eq,n_eq)

    complex(wp) :: rho_inv
    complex(wp) :: ht

    ! Initialize to zero the Jacobian
    v_to_dvdu_complex = cmplx(real(0.0_wp,wp),0.0_wp,wp)

    ! Total enthalpy
    ht = v_in(5) + gm1M2*0.5_wp*(v_in(2)*v_in(2) + v_in(3)*v_in(3) + &
      & v_in(4)*v_in(4))
    
    ! dv/du
    rho_inv = 1.0_wp/v_in(1)

    ! Continuity
    v_to_dvdu_complex(1,1)   = cmplx(real(1.0_wp,wp),0.0_wp,wp)
    
    ! Momentum
    v_to_dvdu_complex(2:4,1) = -rho_inv*v_in(2:4)
    v_to_dvdu_complex(2,2)   = +rho_inv
    v_to_dvdu_complex(3,3)   = +rho_inv
    v_to_dvdu_complex(4,4)   = +rho_inv
    
    ! Energy
    v_to_dvdu_complex(5,1)   = +rho_inv*(-v_in(5) + 0.5_wp*gM2*gm1* &
      & (v_in(2)*v_in(2) + v_in(3)*v_in(3) + v_in(4)*v_in(4)))
    v_to_dvdu_complex(5,2:4) = -rho_inv*v_in(2:4)*gm1*gM2
    v_to_dvdu_complex(5,5)   = +rho_inv*gamma0

    return
  end function v_to_dvdu_complex

  !============================================================================

  !============================================================================
  ! v_to_normal_flux_complex - Calculates the convective (or inviscid) flux in 
  ! the normal direction. Complex variables are used.
  !
  ! Input parameters:
  ! v_in - input primitive variables.
  ! n_v - normal direction.
  ! n_eq - number of equations.
  ! 
  ! Output parameters:
  ! v_to_normal_flux_complex - convective (or inviscid) flux in the normal 
  !                            direction.

  function v_to_normal_flux_complex(v_in,n_v,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gm1M2, gm1og, gM2

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq
    real(wp), dimension(3), intent(in) :: n_v
    complex(wp), intent(in) :: v_in(n_eq)

    complex(wp) :: v_to_normal_flux_complex(n_eq)
    complex(wp) :: u_n
    complex(wp) :: p
    
    ! Initialize to zero the normal flux
    v_to_normal_flux_complex(:) = cmplx(real(0.0_wp,wp),0.0_wp,wp)
    
    ! Normal mass flux
    u_n = v_in(1)*(v_in(2)*n_v(1) + v_in(3)*n_v(2) + v_in(4)*n_v(3))

    ! Calculate pressure (R = 1)
    p = v_in(1)*v_in(5)

    ! Mass flux (\rho u \cdot n)
    v_to_normal_flux_complex(1) = u_n

    ! Momentum flux (\rho u \cdot n) u + p n /(gamma_0 M_0^2)
    v_to_normal_flux_complex(2:4) = u_n*v_in(2:4) + p*n_v/gM2
    
    ! Energy flux (\rho u \cdot n) H
    !v_to_normal_flux_complex(5) = u_n*(v_in(5) &
    !  + gm1M2*0.5_wp*(v_in(2)*v_in(2) + v_in(3)*v_in(3) + v_in(4)*v_in(4)))
    
    v_to_normal_flux_complex(5) = u_n*v_in(5) + u_n*cmplx(real(gm1M2,wp),0.0_wp,wp)* &
      & cmplx(real(0.5_wp,wp),0.0_wp,wp)*v_in(2)*v_in(2) + &
      & u_n*cmplx(real(gm1M2,wp),0.0_wp,wp)*cmplx(real(0.5_wp,wp),0.0_wp,wp)*v_in(3)*v_in(3) + &
      & u_n*cmplx(real(gm1M2,wp),0.0_wp,wp)*cmplx(real(0.5_wp,wp),0.0_wp,wp)*v_in(4)*v_in(4)


    return
  end function v_to_normal_flux_complex

  !============================================================================

  !============================================================================
  ! v_to_entropy_consistent_flux_complex - Calculates the normal entropy 
  ! consistent flux based on the left and the right states in primitive 
  ! variables. It is consistent with the non-dimensionalization employed herein
  ! and follows directly from the work of Ismail and Roe, 
  ! DOI: 10.1016/j.jcp.2009.04.021. Complex variables are used.
  !
  ! Input parameters:
  ! v_l - left state (primitive variables).
  ! v_r - right state (primitive variables).
  ! Jx - normal contravariant vector.
  ! n_eq - number of equations.
  !
  ! Output parameters:
  ! v_to_entropy_consistent_flux_complex - normal entropy consistent flux.
  
  pure function v_to_entropy_consistent_flux_complex(v_l,v_r,n_v,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gM2, gM2I, gm1og, gp1og, gm1M2, gamma0
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_eq
    
    complex(wp), intent(in), dimension(n_eq) :: v_l, v_r
    real(wp), intent(in), dimension(3) :: n_v
  
    complex(wp), dimension(n_eq) :: v_to_entropy_consistent_flux_complex
    
    complex(wp), dimension(n_eq) :: v_hat
    complex(wp) :: root_T_l, root_T_r, root_T_l_inv, root_T_r_inv
    complex(wp) :: T_inv_avg, T_inv_avg_inv
    real(wp), parameter :: sdiv = 1e-015_wp, sdiv2 = 1e-030_wp
    real(wp), parameter :: seps = 1.0e-04_wp
    complex(wp) :: xi, gs, us, ut, f_tmp
    complex(wp) :: s1, s2, al
    complex(wp) :: mdot, p, T

    ! Sqrt[temperature] and Sqrt[temperature]^{-1} are used a lot
    root_T_l = sqrt(v_l(5)) 
    root_T_r = sqrt(v_r(5))
    root_T_l_inv = one/root_T_l
    root_T_r_inv = one/root_T_r
    T_inv_avg = root_T_l_inv + root_T_r_inv
    T_inv_avg_inv = one/T_inv_avg
  
    ! Velocity
    v_hat(2:4) = (v_l(2:4)*root_T_l_inv + v_r(2:4)*root_T_r_inv)*T_inv_avg_inv

    ! Pressure
    f_tmp = v_l(1)*root_T_l + v_r(1)*root_T_r
    p = f_tmp*T_inv_avg_inv

    ! Logarithmic averages used in density and temperature
    ! Calculate s1
    xi = (root_T_l*v_l(1))/(root_T_r*v_r(1))
    gs = (xi - one)/(xi + one)
    us = gs*gs
    al = exp(-us/seps)
    s1 = half/(one + us*(third + us*(fifth + us*(seventh + us*ninth))))
    ut = log(xi)
    s1 = f_tmp*(al*s1 + (one - al)*gs*ut/(ut*ut + sdiv2))

    ! Calculate s2
    xi = root_T_l_inv/root_T_r_inv
    gs = (xi - one)/(xi + one)
    us = gs*gs
    al = exp(-us/seps)
    s2 = half/(one + us*(third + us*(fifth + us*(seventh + us*ninth))))
    ut = log(xi)
    s2 = T_inv_avg*(al*s2 + (one - al)*gs*ut/(ut*ut + sdiv2))

    ! Density
    v_hat(1) = half*T_inv_avg*s1

    ! Temperature
    T = half*(gm1og*p + gp1og*s1/s2)/v_hat(1)

    ! Total enthalpy
    v_hat(5) = gm1M2*half*(v_hat(2)*v_hat(2) + v_hat(3)*v_hat(3) + &
      & v_hat(4)*v_hat(4)) + T

    ! Normal mass flow rate
    mdot = v_hat(1)*(v_hat(2)*n_v(1) + v_hat(3)*n_v(2) + &
      & v_hat(4)*n_v(3))

    v_to_entropy_consistent_flux_complex(1) = mdot
    v_to_entropy_consistent_flux_complex(2) = mdot*v_hat(2) + n_v(1)*P*gM2I
    v_to_entropy_consistent_flux_complex(3) = mdot*v_hat(3) + n_v(2)*P*gM2I
    v_to_entropy_consistent_flux_complex(4) = mdot*v_hat(4) + n_v(3)*P*gM2I
    v_to_entropy_consistent_flux_complex(5) = mdot*v_hat(5)

    return
  end function v_to_entropy_consistent_flux_complex

  !============================================================================
  
  !============================================================================
  ! u_to_v_complex - Calculates the primitive variables from the conservative 
  ! variables in one node.
  !
  ! Input parameters:
  ! u_in - conservative variables.
  ! n_eq - number of equations.

  function u_to_v_complex(u_in,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gm1M2, gm1og
    
    ! Nothing is implicitly defined
    implicit none
    
    integer, intent(in) :: n_eq
    complex(wp), intent(in) :: u_in(n_eq)
    complex(wp) :: u_to_v_complex(n_eq)
    real(wp) :: rs

    ! Specific gas constant
    rs = 1.0_wp

    u_to_v_complex = cmplx(real(0.0_wp,wp),0.0_wp,wp) 

    ! Density
    u_to_v_complex(1) = u_in(1)
    
    ! Velocity
    u_to_v_complex(2) = u_in(2)/u_in(1)
    u_to_v_complex(3) = u_in(3)/u_in(1)
    u_to_v_complex(4) = u_in(4)/u_in(1)
    
    ! Temperature
    u_to_v_complex(5) = (u_in(5)/u_in(1) - &
      & cmplx(real(gm1M2,wp),0.0_wp,wp)*cmplx(real(0.5_wp,wp),0.0_wp,wp)*(u_to_v_complex(2)*u_to_v_complex(2) + &
      & u_to_v_complex(3)*u_to_v_complex(3) + &
      & u_to_v_complex(4)*u_to_v_complex(4)))/(cmplx(real(1.0_wp,wp),0.0_wp,wp)-cmplx(real(gm1og,wp),0.0_wp,wp)* &
      & cmplx(real(rs,wp),0.0_wp,wp))

    return
  end function u_to_v_complex

  !============================================================================
  
  !============================================================================
  ! v_to_w_complex - Calculates the entropy variables from the primitive 
  ! variables in one node. Such variables correspond to the entropy-entropy flux
  ! pair (S,F^i) = (-rho*s,-rho*u_i*s). Complex variables are used.
  !
  ! Input parameters:
  ! v_in - primitive variables.
  ! n_eq - number of equations.

  pure function v_to_w_complex(v_in,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gm1M2, gm1og
    
    ! Nothing is implicitly defined
    implicit none
    
    integer, intent(in) :: n_eq
    complex(wp), dimension(n_eq), intent(in) :: v_in
    complex(wp), dimension(n_eq) :: v_to_w_complex

    ! w_1 = h/T - s - (gamma_0 - 1) M_0^2 u_k u_k/(2T)
    v_to_w_complex(1) = 1.0_wp-0.5_wp*gm1M2*(v_in(2)*v_in(2) + &
      & v_in(3)*v_in(3) + v_in(4)*v_in(4))/v_in(5) - &
      & v_to_specific_entropy_complex(v_in,n_eq)

    ! w_{k+1} = (gamma_0 - 1) M_0^2 u_k/T, k = 1,2,3
    v_to_w_complex(2:4) = gm1M2*v_in(2:4)/v_in(5)
    
    ! w_5 = -1/T
    v_to_w_complex(5) = -1.0_wp/v_in(5)

    return
  end function v_to_w_complex

  !============================================================================

  !============================================================================
  ! v_to_specific_entropy_complex - Calculates the specific thermodynamic  
  ! entropy from the primitive variables. Complex variables are used.
  !
  ! Input parameters:
  ! v_in - primitive variables.
  ! n_eq - number of equations.
  ! 
  ! Ouput parameters:
  ! specific_entropy_complex - specific thermodynamic entropy.
  
  pure function v_to_specific_entropy_complex(v_in,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gm1og
    
    ! Nothing is implicitly defined
    implicit none
    
    integer, intent(in) :: n_eq
    complex(wp), intent(in) :: v_in(n_eq)

    real(wp) :: rs
    complex(wp) :: v_to_specific_entropy_complex

    ! Nondimensional gas constant (R)
    rs = 1.0_wp

    ! Specific thermodynamic entropy
    v_to_specific_entropy_complex = (1.0_wp-gm1og)*rs*log(v_in(5)) &
      - gm1og*rs*log(v_in(1))

    return
  end function v_to_specific_entropy_complex

  !============================================================================

  !============================================================================

  subroutine test_complexification(inode,ielem,nx,tmp_cnt,i_face,bc,knode,kelem)
 
    use variables
    use referencevariables
    use nsereferencevariables
    use controlvariables, only: verbose, timelocal
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol
    use navierstokes
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer, intent(in) :: inode, ielem, tmp_cnt, i_face
    real(wp), dimension(3), intent(in) :: nx

    ! indices
    integer :: knode
    integer :: kelem
    integer :: i

    ! reconstructed flux
    real(wp), allocatable :: fstar(:), fstarV(:)
    ! local normal flux
    real(wp), allocatable :: fn(:), fnV(:)
    ! boundary conservative and primitive states
    real(wp), allocatable :: ustar(:), vstar(:), wstar(:), phistar(:,:)
    ! right and left eigenvector matrices
    real(wp), allocatable, dimension(:,:) :: smat, sinv
    ! eigenvalues
    real(wp), allocatable, dimension(:) :: ev, evabs
    ! average state
    real(wp), allocatable, dimension(:) :: vav
    ! Lax-Freidrich max Eigenvalue
    real(wp) :: evmax
    ! penalty parameter for interfaces
    real(wp) :: Sfix , strength, eps

    real(wp), allocatable :: f_unpert(:), f_pert(:)
    integer :: comp

    real(wp), dimension(nequations,nequations) :: jacob
    real(wp) :: u_bk

    logical, intent(in):: bc

    ! allocate local arrays
    allocate(ustar(nequations))
    allocate(vstar(nequations))
    allocate(wstar(nequations))
    allocate(fstar(nequations))
    allocate(fstarV(nequations))
    allocate(fn(nequations))
    allocate(fnV(nequations))
    allocate(phistar(nequations,3))

    allocate(ev(nequations))
    allocate(evabs(nequations))
    allocate(vav(nequations))
    allocate(smat(nequations,nequations))
    allocate(sinv(nequations,nequations))

    allocate(f_unpert(nequations))
    allocate(f_pert(nequations))
    
    Sfix     = 0.0001_wp
    strength = 2.0_wp

    eps = 1.e-8_wp

    !write(*,*) 'nx', nx
    !stop

    if (bc .eqv. .true.) then
      do comp = 1, nequations
        
        ! compute the boundary state
        call primitivevariables( ug(:,inode,ielem), &
          vg(:,inode,ielem), &
          nequations ) ! (navierstokes)

        call entropyvariables(vg(:,inode,ielem),wg(:,inode,ielem),nequations)

        vstar(:) = vg(:,inode,ielem)

        call isentropicVortexFull( &
          Vx  = vstar, &
          phi = phig(:,:,inode,ielem), &
          fv = fnV, &
          Jx = nx(:), &
          xin = xg(:,inode,ielem), &
          tin = timelocal, &
          neqin = nequations, &
          nd    = ndim, &
          mut = mut(inode,ielem) )

        ! ==  Eigen values/vectors
        ! calculate the conserved variables
        call conservativevariables( vstar, ustar, nequations) ! (navierstokes)
        call entropyvariables(vstar, wstar, nequations) ! (navierstokes)
        ! compute the roe average of the solution and boundary states
        call roeavg( vg(:,inode,ielem),&
          vstar, &
          Vav, &
          nequations ) ! (navierstokes)

        ! compute the characteristic decomposition in the normal direction
        call CharacteristicDecomp( vav, &
          nequations, &
          sinv, &
          smat, &
          ev, &
          nx ) ! (navierstokes)
        evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)
        ! ==  Eigen values/vectors

        ! ==  Fluxes
        fn = normalflux( &
          vg(:,inode,ielem), &
          nx, &
          nequations ) ! (navierstokes)

        fstar = EntropyConsistentFlux(vg(:,inode,ielem), &
          vstar, &
          nx, &
          nequations ) ! (navierstokes)

        fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-wstar) )

        f_unpert = 0.0_wp
        !f_unpert = pinv(1) * (fn-fstar)
        f_unpert = pinv(1) * fn
        !f_unpert = vav
        !f_unpert = vg(:,inode,ielem)
        !f_unpert = vstar
        !f_unpert = fn
        !f_unpert = pinv(1)*fn
        !f_unpert = -pinv(1)*fstar
        !f_unpert =  wg(:,inode,ielem)
        !f_unpert =  wstar
        !f_unpert = sinv(:,5)

        u_bk = ug(comp,inode,ielem)
        ug(comp,inode,ielem) = ug(comp,inode,ielem) + eps

        vg(:,inode,ielem) = 0.0_wp

        ! compute the boundary state
        call primitivevariables( ug(:,inode,ielem), &
          vg(:,inode,ielem), &
          nequations ) ! (navierstokes)

        call entropyvariables(vg(:,inode,ielem),wg(:,inode,ielem),nequations)

        vstar = 0.0_wp
        vstar(:) = vg(:,inode,ielem)

        call isentropicVortexFull( &
          Vx  = vstar, &
          phi = phig(:,:,inode,ielem), &
          fv = fnV, &
          Jx = nx(:), &
          xin = xg(:,inode,ielem), &
          tin = timelocal, &
          neqin = nequations, &
          nd    = ndim, &
          mut   = mut(inode,ielem) )

        ! ==  Eigen values/vectors
        ! calculate the conserved variables
        ustar = 0.0_wp
        call conservativevariables( vstar, ustar, nequations) ! (navierstokes)
        
        wstar = 0.0_wp
        call entropyvariables(vstar, wstar, nequations) ! (navierstokes)
        ! compute the roe average of the solution and boundary states
        
        vav = 0.0_wp
        call roeavg( vg(:,inode,ielem),&
          vstar, &
          Vav, &
          nequations ) ! (navierstokes)
        ! compute the characteristic decomposition in the normal direction
        
        
        sinv = 0.0_wp
        smat = 0.0_wp
        ev = 0.0_wp
        call CharacteristicDecomp( vav, &
          nequations, &
          sinv, &
          smat, &
          ev, &
          nx ) ! (navierstokes)
        evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)
        ! ==  Eigen values/vectors

        fn = 0.0_wp
        ! ==  Fluxes
        fn = normalflux( &
          vg(:,inode,ielem), &
          nx, &
          nequations ) ! (navierstokes)

        fstar = 0.0_wp
        fstar = EntropyConsistentFlux(vg(:,inode,ielem), &
          vstar, &
          nx, &
          nequations ) ! (navierstokes)

        fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-wstar) )

        f_pert = 0.0_wp
        !f_pert = pinv(1) * (fn-fstar)
        f_pert = pinv(1) * fn
        !f_pert = vav
        !f_pert = vg(:,inode,ielem)
        !f_pert = vstar
        !f_pert = fn
        !f_pert = pinv(1)*fn
        !f_pert = -pinv(1)*fstar
        !f_pert = wg(:,inode,ielem)
        !f_pert = wstar
        !f_pert = sinv(:,5)
        
        jacob(:,comp) = (f_pert-f_unpert)/eps
        
        ug(comp,inode,ielem) = u_bk

        call primitivevariables( ug(:,inode,ielem), &
          vg(:,inode,ielem), &
          nequations )


        call entropyvariables(vg(:,inode,ielem),wg(:,inode,ielem),nequations)
        

      enddo

    else

      do comp = 1, nequations


        call roeavg( vg(:,inode,ielem),&
          vg(:,knode,kelem), &
          Vav, &
          nequations ) ! (navierstokes)
        ! compute the characteristic decomposition using the averaged state
        call CharacteristicDecomp( vav, &
          nequations, &
          sinv, &
          smat, &
          ev, &
          nx ) ! (navierstokes)
        evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)

 ! ==  Eigen values/vectors


  ! ==  Fluxes
        fn = normalflux( &
          vg(:,inode,ielem), &
          nx, &
          nequations ) ! (navierstokes)

        fstar = EntropyConsistentFlux(vg(:,inode,ielem), &
          vg(:,knode,kelem), &
          nx, &
          nequations ) ! (navierstokes)

        fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-wg(:,knode,kelem)) )

        f_unpert = 0.0_wp
        !f_unpert = pinv(1) * (fn-fstar)
        f_unpert = pinv(1) * fstar
        !f_unpert = vav
        !f_unpert = vg(:,inode,ielem)
        !f_unpert = vstar
        !f_unpert = fn
        !f_unpert = pinv(1)*fn
        !f_unpert = -pinv(1)*fstar
        !f_unpert =  wg(:,inode,ielem)
        !f_unpert =  wstar
        !f_unpert = sinv(:,5)

        u_bk = ug(comp,inode,ielem)
        ug(comp,inode,ielem) = ug(comp,inode,ielem) + eps

        vg(:,inode,ielem) = 0.0_wp

        ! compute the boundary state
        call primitivevariables( ug(:,inode,ielem), &
          vg(:,inode,ielem), &
          nequations ) ! (navierstokes)

        call entropyvariables(vg(:,inode,ielem),wg(:,inode,ielem),nequations)

        ! compute the roe average of the solution and boundary states
        vav = 0.0_wp
        call roeavg( vg(:,inode,ielem),&
          vg(:,knode,kelem), &
          Vav, &
          nequations ) ! (navierstokes)
        ! compute the characteristic decomposition in the normal direction
        
        !write(*,*) 'vav', vav 
        !write(*,*) 'nx', nx 
        sinv = 0.0_wp
        smat = 0.0_wp
        ev = 0.0_wp
        call CharacteristicDecomp( vav, &
          nequations, &
          sinv, &
          smat, &
          ev, &
          nx ) ! (navierstokes)
        evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)
        ! ==  Eigen values/vectors
        
        !write(*,*) 'smat'
        !do l = 1,nequations
        !  write(*,*) smat(l,:)
        !enddo
        !write(*,*) 'eig', ev
        !write(*,*) 'sinv'
        !do l = 1,nequations
        !  write(*,*) sinv(l,:)
        !enddo
        !write(*,*) 'abs eig', evabs

        fn = 0.0_wp
        ! ==  Fluxes
        
        fn = normalflux( &
          vg(:,inode,ielem), &
          nx, &
          nequations ) ! (navierstokes)

        fstar = EntropyConsistentFlux(vg(:,inode,ielem), &
          vg(:,knode,kelem), &
          nx, &
          nequations ) ! (navierstokes)
 
        fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-wg(:,knode,kelem)) )

        !write(*,*) 'i_face', i_face
        !write(*,*) 'inode', inode
        !        write(*,*) 'wg inode',wg(:,inode,ielem)
        !write(*,*) 'wg knode',wg(:,knode,kelem)

        !write(*,*) 'eig', ev
 
 
        f_pert = 0.0_wp
        !f_pert = pinv(1) * (fn-fstar)
        f_pert = pinv(1) * fstar
        !f_pert = vav
        !f_pert = vg(:,inode,ielem)
        !f_pert = vstar
        !f_pert = fn
        !f_pert = pinv(1)*fn
        !f_pert = -pinv(1)*fstar
        !f_pert = wg(:,inode,ielem)
        !f_pert = wstar
        !f_pert = sinv(:,5)
        
        jacob(:,comp) = (f_pert-f_unpert)/eps
        
        ug(comp,inode,ielem) = u_bk

        call primitivevariables( ug(:,inode,ielem), &
          vg(:,inode,ielem), &
          nequations )


        call entropyvariables(vg(:,inode,ielem),wg(:,inode,ielem),nequations)
      

      enddo
 
    endif

    write(*,*) 'Fréchet derivative'
    do i = 1, 5
      write(*,*) jacob(i,:)
    enddo

    !mat = inviscid_flux_jacobian_node(vg(:,inode,ielem),nx,nequations)
    !norm = maxval(abs(jacob-mat*pinv(1)))
    !write(*,*) 'Norm', norm

    if (tmp_cnt == nodesperface .and. i_face == 4) then
      stop
    endif

  end subroutine test_complexification

  !============================================================================
 

end module jacobian_matrix_implicit_ts


