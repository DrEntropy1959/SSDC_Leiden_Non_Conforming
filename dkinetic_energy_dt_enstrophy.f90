! This module contains the subroutines and functions to compute the integral of
! the time derivative of the kinetic energy (dissipation) and the enstrophy
! =============================================================================

module dkinetic_energy_dt_enstrophy
  
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  private

  public compute_enstrophy
  public compute_dkinetic_energy_dt

contains

  subroutine compute_enstrophy()
    
    ! Load modules
    use variables, only: Jx_r, vg, omega, enstrophy
    use collocationvariables, only: pvol
    use controlvariables
    use referencevariables
    use navierstokes
    use mpimod

    ! Nothing is implicitly defined
    implicit none

    integer :: i_node, i_elem

    integer :: low_elem, high_elem

    real(wp) :: local_enstrophy, enstrophy_sum
    
    integer :: i_err

    continue

    ! Low volumetric element index
    low_elem = ihelems(1)
    
    ! High volumetric element index
    high_elem = ihelems(2)

    ! Initialize local_enstrophy to zero
    local_enstrophy = 0.0_wp

!    ! Compute the vorticity field for all the elements own by a process
!    call compute_vorticity_field_elements()
    
    ! Loop over the elements own by a process
    do i_elem = low_elem, high_elem
      
      ! Loop over each node in the element
      do i_node = 1, nodesperelem
        
        ! Calculate the integral contribution to the enstrophy
        ! pvol*J is the volumetric integration weight.
        local_enstrophy = local_enstrophy &
          & + pvol(i_node)*Jx_r(i_node,i_elem)*vg(1,i_node,i_elem)*dot_product(omega(:,i_node,i_elem),omega(:,i_node,i_elem))/2.0_wp
!          & + pvol(i_node)*Jx_r(i_node,i_elem)*dot_product(omega(:,i_node,i_elem),omega(:,i_node,i_elem))/2.0_wp

      end do
    end do

    ! MPI all reduce
    call mpi_allreduce(local_enstrophy,enstrophy_sum,1, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)

    ! Set the global value of the enstrophy
    ! Volume is (2*Pi)
    enstrophy = enstrophy_sum/((2*pi)**3)

    return
  end subroutine compute_enstrophy

  !============================================================================

  !============================================================================

  subroutine compute_dkinetic_energy_dt()
    
    ! Load modules
    use variables, only: Jx_r, ug, vg, dudt, dkinetic_energy_dt, kinetic_energy
    use collocationvariables, only: pvol
    use controlvariables
    use referencevariables
    use navierstokes
    use mpimod

    ! Nothing is implicitly defined
    implicit none

    integer :: i_node, i_elem

    integer :: low_elem, high_elem

!    real(wp), dimension(nodesperelem) :: ke_elem
    
    real(wp) :: drho_dt_times_ke_elem, vel_dot_dmomentum_dt

    real(wp) :: local_integrand_dke_dt, local_integrand_ke

    real(wp), dimension(2) :: local_ke, ke_sum
    
    integer :: i_err

    continue

    ! Low volumetric element index
    low_elem = ihelems(1)
    
    ! High volumetric element index
    high_elem = ihelems(2)

    ! Initialize local_ke to zero
    local_ke(:) = 0.0_wp

    ! Compute the vorticity field for all the elements own by a process
    call compute_vorticity_field_elements()
    
    ! Loop over the elements own by a process
    do i_elem = low_elem, high_elem

!      ! Compute kinetic energy in the element with ID i_elem
!      ke_elem = kinetic_energy_element(i_elem) 

      ! Loop over each node in the element
      do i_node = 1, nodesperelem

        ! First piece of integrand of dke_dt: drho/dt*ke
        ! Explicit RK scheme: drho/dt = dudt(1,i_node,i_elem)
        ! drho/dt = dudt(1,i_node,i_elem) is the residual of the continuity
        ! equation
        ! ===================================================================
        drho_dt_times_ke_elem = - 0.5_wp*dudt(1,i_node,i_elem)*dot_product(vg(2:4,i_node,i_elem),vg(2:4,i_node,i_elem))
!       drho_dt_times_ke_elem = - 0.5_wp*dudt(1,i_node,i_elem)*dot_product(ug(2:4,i_node,i_elem),ug(2:4,i_node,i_elem))&
!                             &  /(ug(1,i_node,i_elem)**2)


        ! Second piece of the integrand of dke_dt: vec{v}\dot d(rho\vec{v})/dt
        ! Explicit RK scheme: d(rho\vec{v})/dt = dudt(2:4,i_node,i_elem)
        ! d(rho\vec{v})/dt = dudt(2:4,i_node,i_elem) is the residual of the
        ! momentum equation
        ! ====================================================================
        vel_dot_dmomentum_dt = dot_product(vg(2:4,i_node,i_elem),dudt(2:4,i_node,i_elem))
!       vel_dot_dmomentum_dt = dot_product(ug(2:4,i_node,i_elem),dudt(2:4,i_node,i_elem))/(ug(1,i_node,i_elem))

        ! Sum first and second contribution
        ! =================================
        local_integrand_dke_dt = drho_dt_times_ke_elem + vel_dot_dmomentum_dt
!       local_integrand = vel_dot_dmomentum_dt
        
        ! Calculate the integral contribution to the dke/dt
        ! pvol*J is the volumetric integration weight.
        local_ke(2) = local_ke(2) &
          & + pvol(i_node)*Jx_r(i_node,i_elem)*local_integrand_dke_dt

        ! Calculate the node value of ke
!       local_integrand_ke = 0.5_wp*dot_product(vg(2:4,i_node,i_elem),vg(2:4,i_node,i_elem))   !   No density
        local_integrand_ke = 0.5_wp*dot_product(ug(2:4,i_node,i_elem),vg(2:4,i_node,i_elem))   !   With density


        ! Calculate the integral contribution to ke
        local_ke(1) = local_ke(1) + pvol(i_node)*Jx_r(i_node,i_elem)*local_integrand_ke

      end do
    end do

    ! MPI all reduce
    call mpi_allreduce(local_ke,ke_sum,2, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)

    ! Set the global value of dke/dt
    ! Volume is (2*Pi)
    dkinetic_energy_dt = -1.0_wp*ke_sum(2)/((2*pi)**3)

    ! Set global value of ke
    kinetic_energy = ke_sum(1)/((2*pi)**3)

    return
  end subroutine compute_dkinetic_energy_dt

  !============================================================================

  !============================================================================


end module dkinetic_energy_dt_enstrophy
