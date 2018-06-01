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
    use initcollocation,       only: element_properties
    use variables, only: Jx_r, vg, omega, enstrophy
    use collocationvariables, only: pvol
    use controlvariables
    use referencevariables
    use navierstokes
    use mpimod

    ! Nothing is implicitly defined
    implicit none

    integer  :: inode, ielem
    integer  :: n_pts_3d

    real(wp) :: local_enstrophy, enstrophy_sum
    
    integer  :: i_err

    continue

    local_enstrophy = 0.0_wp                                    ! Initialize local_enstrophy to zero  

!   call compute_vorticity_field_elements()                     ! Vorticity field for all the elements own by a process
    
    do ielem = ihelems(1), ihelems(2)                           ! Loop over the elements own by a process
      
      call element_properties(ielem, n_pts_3d=n_pts_3d, pvol=pvol)

      do inode = 1, n_pts_3d                                ! Loop over each node in the element
        
        ! Calculate the integral contribution to the enstrophy
        ! pvol*J is the volumetric integration weight.
        local_enstrophy = local_enstrophy &
          & + pvol(inode)*Jx_r(inode,ielem)*vg(1,inode,ielem)*dot_product(omega(:,inode,ielem),omega(:,inode,ielem))/2.0_wp
!         & + pvol(inode)*Jx_r(inode,ielem)*dot_product(omega(:,inode,ielem),omega(:,inode,ielem))/2.0_wp

      end do
    end do

    call mpi_allreduce(local_enstrophy,enstrophy_sum,1,  &       ! MPI all reduce
                       MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)

    enstrophy = enstrophy_sum/((2*pi)**3)                        ! Set the global value of the enstrophy;  ! Volume is (2*Pi)

  end subroutine compute_enstrophy

  !============================================================================

  subroutine compute_dkinetic_energy_dt()
    
    ! Load modules
    use initcollocation,       only: element_properties
    use variables, only: Jx_r, ug, vg, dudt, dkinetic_energy_dt, kinetic_energy
    use collocationvariables, only: pvol
    use controlvariables
    use referencevariables
    use navierstokes
    use mpimod

    ! Nothing is implicitly defined
    implicit none

    integer :: inode, ielem
    integer :: n_pts_3d

    real(wp) :: drho_dt_times_ke_elem, vel_dot_dmomentum_dt

    real(wp) :: local_integrand_dke_dt, local_integrand_ke

    real(wp), dimension(2) :: local_ke, ke_sum
    
    integer :: i_err

    continue

    local_ke(:) = 0.0_wp                                        ! Initialize local_ke to zero

    call compute_vorticity_field_elements()                     ! Compute the vorticity field for all the elements own by a process
    
    do ielem = ihelems(1), ihelems(2)                           ! Loop over the elements own by a process

      call element_properties(ielem, n_pts_3d=n_pts_3d, pvol=pvol)

!      ke_elem = kinetic_energy_element(ielem)                  ! Compute kinetic energy in the element with ID ielem

      do inode = 1, n_pts_3d                                    ! Loop over each node in the element

        ! First piece of integrand of dke_dt: drho/dt*ke
        ! Explicit RK scheme: drho/dt = dudt(1,inode,ielem)
        ! drho/dt = dudt(1,inode,ielem) is the residual of continuity equation
        ! ===================================================================
        drho_dt_times_ke_elem = - 0.5_wp*dudt(1,inode,ielem)*dot_product(vg(2:4,inode,ielem),vg(2:4,inode,ielem))
!       drho_dt_times_ke_elem = - 0.5_wp*dudt(1,inode,ielem)*dot_product(ug(2:4,inode,ielem),ug(2:4,inode,ielem))&
!                             &  /(ug(1,inode,ielem)**2)


        ! Second piece of the integrand of dke_dt: vec{v}\dot d(rho\vec{v})/dt
        ! Explicit RK scheme: d(rho\vec{v})/dt = dudt(2:4,inode,ielem)
        ! d(rho\vec{v})/dt = dudt(2:4,inode,ielem) is the residual of momentum equation
        ! ====================================================================
        vel_dot_dmomentum_dt = dot_product(vg(2:4,inode,ielem),dudt(2:4,inode,ielem))
!       vel_dot_dmomentum_dt = dot_product(ug(2:4,inode,ielem),dudt(2:4,inode,ielem))/(ug(1,inode,ielem))

        ! Sum first and second contribution
        ! =================================
        local_integrand_dke_dt = drho_dt_times_ke_elem + vel_dot_dmomentum_dt
!       local_integrand = vel_dot_dmomentum_dt
        
        ! Calculate the integral contribution to the dke/dt
        ! pvol*J is the volumetric integration weight.
        local_ke(2) = local_ke(2) + pvol(inode)*Jx_r(inode,ielem)*local_integrand_dke_dt

        ! Calculate the node value of ke
!       local_integrand_ke = 0.5_wp*dot_product(vg(2:4,inode,ielem),vg(2:4,inode,ielem))   !   No density
        local_integrand_ke = 0.5_wp*dot_product(ug(2:4,inode,ielem),vg(2:4,inode,ielem))   !   With density


        ! Calculate the integral contribution to ke
        local_ke(1) = local_ke(1) + pvol(inode)*Jx_r(inode,ielem)*local_integrand_ke

      end do
    end do

    call mpi_allreduce(local_ke,ke_sum,2,MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)  ! MPI all reduce

    dkinetic_energy_dt = -1.0_wp*ke_sum(2)/((2*pi)**3)               ! Set the global value of dke/dt; Volume is (2*Pi)

    kinetic_energy = ke_sum(1)/((2*pi)**3)                           ! Set global value of ke

  end subroutine compute_dkinetic_energy_dt


end module dkinetic_energy_dt_enstrophy
