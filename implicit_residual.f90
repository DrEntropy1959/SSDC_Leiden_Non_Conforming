! This module contains all the subroutines to compute the residual of the
! implicit time integration scheme and the stage-value predictor.

module implicit_residual

  ! Load modules
  use precision_vars
  use time_integ_coeff 
  
  ! Nothing is implicitly defined
  implicit none

contains

  !============================================================================

  !============================================================================
  ! compute_implicit_residual_imex - Computes the residual of the implicit
  ! portion of the RK-IMEX scheme for all the elements.
  !
  ! Input parameters:
  ! irkstep - stage number of the RK-IMEX scheme.
  ! dt  - time-step.
  !
  ! Output parameter:
  ! non_lin_res - residual of the implicit portion of the RK-IMEX scheme.

  subroutine compute_implicit_residual_imex(irkstep,dt) 

    ! Load modules
    use variables, only: ug, uexp, Fimp, non_lin_res
    use controlvariables
    use navierstokes

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: irkstep
    real(wp), intent(in) :: dt

    ! Calculate implicit residual, i.e. Fimp
    call nse_calcrhsimplicit(irkstep,timelocal)

    ! Create nonlinear RK_IMEX stage residual, i.e. 
    ! rnonlin = u^k - uexp - dt * A_ii f(u^k) * Fimp(:,:,:,istep)
    ! Note: the function evaluation at each stage are stored in Fexp and 
    ! Fimp
    non_lin_res(:,:,:) = ug(:,:,:) - uexp(:,:,:) &
      - dt * arkimp(irkstep,irkstep) * Fimp(:,:,:,irkstep) 
    !non_lin_res(:,:,:) = -Fimp(:,:,:,irkstep) 
    
    return
  end subroutine compute_implicit_residual_imex

  !============================================================================

  !============================================================================
  ! rk_imex_svp - Computes the stage-value prediction of the solution of the
  ! RK-IMEX scheme for all the elements.
  !
  ! Input parameters:
  ! irkstep - stage number of the RK-IMEX scheme.
  ! dt  - time-step.
  !
  ! Output parameter:
  ! ug - stage-value prediction. It is overwritten in ug. 

  subroutine rk_imex_svp(irkstep,dt)
    
    ! Load modules
    use referencevariables
    use variables
    use navierstokes

    ! Nothing is implicitly defined
    integer, intent(in) :: irkstep
    real(wp), intent(in) :: dt
    integer :: i

    ! Set vector of conservative variables equal to the old vector of the
    ! conservative variables
    ug = uold

    ! Compute the stage-value prediction
    if (irkstep<3) return
    do i = 1,irkstep-1
      ug(:,:,:) = ug(:,:,:) &
        & + dt*svp   (irkstep,i)*Fimp(:,:,:,i) &
        & + dt*arkexp(irkstep,i)*Fexp(:,:,:,i)
    end do
    
    return
  end subroutine rk_imex_svp

  !============================================================================

end module implicit_residual
