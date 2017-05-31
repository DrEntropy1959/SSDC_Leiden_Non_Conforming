! This module contains the necessary routines to compute the time-averaged
! quantities for for Euler and Navier-Stokes equations.  

module time_average

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Subroutines in this module are generally private 
  private

  ! Exceptions, i.e. public subroutines or functions
  public compute_time_averaged_quantities
  public compute_reynolds_stress

contains

  !============================================================================

  !============================================================================
  ! compute_time_average_quantities - Drives the computation of the 
  !                                   time-averaged quantities. 
  !
  ! Input parameters:
  ! time_step_cnt  - time-step counter.

  subroutine compute_time_averaged_quantities(time_step_cnt)

    ! Load modules
    use variables

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: time_step_cnt

    ! Compute sum of time-averaged of the primitive variables
    call compute_time_averaged_primitive_variables(time_step_cnt)

    ! Compute sum of the time-averaged of the product of the velocity components
    call compute_time_averaged_product_velocity_components(time_step_cnt)

    ! Compute sum of Reynolds stress
    call compute_reynolds_stress()

    return
  end subroutine compute_time_averaged_quantities

  !============================================================================

  !============================================================================
  ! compute_time_averaged_primitive_variables - Computes the time-average of the
  !                                             primitive variables.
  !
  ! Input parameters:
  ! time_step_cnt  - time-step counter.

  subroutine compute_time_averaged_primitive_variables(time_step_cnt)

    ! Load modules
    use referencevariables
    use controlvariables, only : restart_time_steps
    use variables, only : vg, mean_vg
      
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: time_step_cnt

    ! Compute <rho>, <u>, <v>, <w>, <T>   
    mean_vg = (mean_vg * (restart_time_steps + time_step_cnt - 1) + vg) / &
      & (restart_time_steps + time_step_cnt)

    return
  end subroutine compute_time_averaged_primitive_variables

  !============================================================================
  
  !============================================================================
  ! compute_time_averaged_product_velocity_components - Computes the 
  !                                                     time-average of the
  !                                                     product of the velocity
  !                                                     components.
  !
  ! Input parameters:
  ! time_step_cnt  - time-step counter.

  subroutine compute_time_averaged_product_velocity_components(time_step_cnt)

    ! Load modules
    use referencevariables
    use controlvariables, only : restart_time_steps
    use variables, only : vg, time_ave_prod_vel_comp
      
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: time_step_cnt
    
    integer :: elem_low, elem_high
    integer :: i_elem, i_node


    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)


    ! Compute time-averaged of the product of the velocity components 
    forall(i_node = 1:nodesperelem, i_elem = elem_low:elem_high)

      ! <u*u> 
      time_ave_prod_vel_comp(1,i_node,i_elem) = (time_ave_prod_vel_comp(1,i_node,i_elem) * &
        (restart_time_steps + time_step_cnt -1) + vg(2,i_node,i_elem)*vg(2,i_node,i_elem)) / &
        (restart_time_steps + time_step_cnt)
      
      ! <u*v>
      time_ave_prod_vel_comp(2,i_node,i_elem) = (time_ave_prod_vel_comp(2,i_node,i_elem) * &
        (restart_time_steps + time_step_cnt -1) + vg(2,i_node,i_elem)*vg(3,i_node,i_elem)) / &
        (restart_time_steps + time_step_cnt)
   
      ! <u*w>
      time_ave_prod_vel_comp(3,i_node,i_elem) = (time_ave_prod_vel_comp(3,i_node,i_elem) * & 
        (restart_time_steps + time_step_cnt -1) + vg(2,i_node,i_elem)*vg(4,i_node,i_elem)) / &
        (restart_time_steps + time_step_cnt)

      ! <v*v> 
      time_ave_prod_vel_comp(4,i_node,i_elem) = (time_ave_prod_vel_comp(4,i_node,i_elem) * & 
        (restart_time_steps + time_step_cnt -1) + vg(3,i_node,i_elem)*vg(3,i_node,i_elem)) / &
        (restart_time_steps + time_step_cnt)

      
      ! <v*w> 
      time_ave_prod_vel_comp(5,i_node,i_elem) = (time_ave_prod_vel_comp(5,i_node,i_elem) * & 
        (restart_time_steps + time_step_cnt -1) + vg(3,i_node,i_elem)*vg(4,i_node,i_elem)) / &
        (restart_time_steps + time_step_cnt)

      ! <w*w> 
      time_ave_prod_vel_comp(6,i_node,i_elem) = (time_ave_prod_vel_comp(6,i_node,i_elem) * & 
        (restart_time_steps + time_step_cnt -1) + vg(4,i_node,i_elem)*vg(4,i_node,i_elem)) / &
        (restart_time_steps + time_step_cnt)

    end forall

    return
  end subroutine compute_time_averaged_product_velocity_components

  !============================================================================
  
  !============================================================================
  ! compute_reynolds_stress - Computes the Reynolds stresses.

  subroutine compute_reynolds_stress()

    ! Load modules
    use referencevariables
    use variables, only : vg, mean_vg, time_ave_prod_vel_comp, &
                          reynolds_stress
      
    ! Nothing is implicitly defined
    implicit none

    integer :: elem_low, elem_high
    integer :: i_elem, i_node


    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)


    ! Compute Reynolds stress
    forall(i_node = 1:nodesperelem, i_elem = elem_low:elem_high)

      ! <u'*u'> 
      reynolds_stress(1,i_node,i_elem) = time_ave_prod_vel_comp(1,i_node,i_elem) - &
        mean_vg(2,i_node,i_elem)*mean_vg(2,i_node,i_elem)
      
      ! <u'*v'> 
      reynolds_stress(2,i_node,i_elem) = time_ave_prod_vel_comp(2,i_node,i_elem) - &
        mean_vg(2,i_node,i_elem)*mean_vg(3,i_node,i_elem)
      
      ! <u'*w'> 
      reynolds_stress(3,i_node,i_elem) = time_ave_prod_vel_comp(3,i_node,i_elem) - &
        mean_vg(2,i_node,i_elem)*mean_vg(4,i_node,i_elem)

      ! <v'*v'> 
      reynolds_stress(4,i_node,i_elem) = time_ave_prod_vel_comp(4,i_node,i_elem) - &
        mean_vg(3,i_node,i_elem)*mean_vg(3,i_node,i_elem)
      
      ! <v'*w'> 
      reynolds_stress(5,i_node,i_elem) = time_ave_prod_vel_comp(5,i_node,i_elem) - &
        mean_vg(3,i_node,i_elem)*mean_vg(4,i_node,i_elem)
      
      ! <w'*w'> 
      reynolds_stress(6,i_node,i_elem) = time_ave_prod_vel_comp(6,i_node,i_elem) - &
        mean_vg(4,i_node,i_elem)*mean_vg(4,i_node,i_elem)

    end forall

    return
  end subroutine compute_reynolds_stress

  !============================================================================

end module time_average

