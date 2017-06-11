! This module contains the necessary routines to compute the aerodynamic
! coefficients such as lift and drag coefficients. 
! NOTE: These subroutines are going to build the coefficients assuming that all
! the solid wall faces belong to one body.


module aerodynamic_coefficients

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Subroutines in this module are generally private 
  private

  ! Exceptions, i.e. public subroutines or functions
  public compute_aerodynamic_coefficients

  contains

  !============================================================================

  !============================================================================
  ! compute_aerodynamic_coefficients - Drives the computation of the aerodynamic
  ! coefficients. 
  !
  ! Forces:
  !
  ! h_force = -p*nx + mu/Re*[ nx*(2*dUxdx - 2* dUxdx) + ny*(dUydx + dUxdy)     + nz*(dUzdx + dUxdz)  ] 
  ! v_force = -p*ny + mu/Re*[ nx*(dUydx + dUxdy)      + ny*(2*dUydy - 2*dUydy) + nz*(dUzdy + dUydz)  ]
  ! s_force = -p*nz + mu/Re*[ nx*(dUzdx + dUxdz)      + ny*(dUzdy + dUydz)     + nz*(2*dUzdz - dUzdz)]


  subroutine compute_aerodynamic_coefficients()

    ! Load modules
    use variables
    use referencevariables
    use navierstokes, only : dVdW, normalviscousflux
    use nsereferencevariables, only : mu0, gM2, aero_coeffs_surface
    use collocationvariables, only : pmat_Soln
    use mpimod
    use controlvariables

    ! Nothing is implicitly defined
    implicit none

    integer :: iell, elem_high
    integer :: i_elem, i_node, j_node
    integer :: n_wall_elems
    integer :: i, j, elem_id, face_id

    real(wp), dimension(3) :: normal

    real(wp) :: p

    real(wp), dimension(nodesperedge,nodesperedge) :: h_force_face_nodes, &
                                                    & v_force_face_nodes, &
                                                    & s_force_face_nodes
    integer :: cnt

    real(wp) :: h_force, v_force, s_force
    
    real(wp), dimension(nequations) :: n_visc_flux
    
    real(wp), dimension(3) :: area

    real(wp) :: h_force_sum, v_force_sum, s_force_sum

    integer :: i_err


    continue

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements whcih own a wall face 
    n_wall_elems = size(wall_elem_face_ids(1,:))

    ! Initialize to zero the forces
    h_force = 0.0_wp
    v_force = 0.0_wp
    s_force = 0.0_wp

    h_force_sum = 0.0_wp
    v_force_sum = 0.0_wp
    s_force_sum = 0.0_wp

    ! Initialize to zero the area
    area = 0.0_wp
    area_sum = 0.0_wp

    cnt = 0
    h_force_face_nodes = 0.0_wp
    v_force_face_nodes = 0.0_wp
    s_force_face_nodes = 0.0_wp

    do i_elem = 1, n_wall_elems
      elem_id = wall_elem_face_ids(1,i_elem)
      face_id = wall_elem_face_ids(2,i_elem)

      cnt = 0
      do j = 1, nodesperedge
       
        do i = 1, nodesperedge

          cnt = cnt + 1

          ! Volumetric node index corresponding to face and node on face indices
          i_node = kfacenodes(cnt,face_id)
         
          ! Facial index corresponding to face and node on face indices
          j_node = (face_id-1)*nodesperface + cnt
         
          ! Outward facing normal of facial node
          normal = -1.0_wp*Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          ! Compute the pressure at the node
          ! Note: the factor 1/gM2 comes from the non-dimensionalization
          p = vg(1,i_node,elem_id)*vg(5,i_node,elem_id)/gM2

          ! Contribution of the pressure to the horizontal force
          h_force_face_nodes(i,j) = -p*normal(1)

          ! Contribution of the pressure to the vertical force
          v_force_face_nodes(i,j) = -p*normal(2)

          ! Contribution of the pressure to the side force
          s_force_face_nodes(i,j) = -p*normal(3)

! 1st approach: not used
! ----------------------
!          ! Compute the gradient of the primitive variables
!          ! grad_prim = dprim/dentr * grad_entr
!          grad_prim = MatMul(dVdW(vg(:,i_node,elem_id),nequations),phig(:,:,i_node,elem_id))
!
!          ! Divergence of the velocity
!          divergence = grad_prim(2,1) + grad_prim(3,2) + grad_prim(4,3)
!
!
!          ! Add contribution of the viscous stress to the horizontal force
!          h_force_face_nodes(i,j) = h_force_face_nodes(i,j) - mu0*Re0inv*( normal(1)*(grad_prim(2,1) + grad_prim(2,1)  &
!                                                                         & - 2.0_wp/3.0_wp*divergence) + &
!                                                                         & normal(2)*(grad_prim(2,2) + grad_prim(3,1)) + &
!                                                                         & normal(3)*(grad_prim(2,3) + grad_prim(4,1)))
!
!          ! Add contribution of the viscous stress to the vertical force
!          v_force_face_nodes(i,j) = v_force_face_nodes(i,j) - mu0*Re0inv*( normal(1)*(grad_prim(3,1) + grad_prim(2,2)) + &
!                                                                         & normal(2)*(grad_prim(3,2) + grad_prim(3,2)  &
!                                                                         & - 2.0_wp/3.0_wp*divergence) + &
!                                                                         & normal(3)*(grad_prim(3,3) + grad_prim(4,2)))
!
!          ! Add contribution of the viscous stress to the side force
!          s_force_face_nodes(i,j) = s_force_face_nodes(i,j) - mu0*Re0inv*( normal(1)*(grad_prim(4,1) + grad_prim(2,3)) + &
!                                                                         & normal(2)*(grad_prim(4,2) + grad_prim(3,3)) + &
!                                                                         & normal(3)*(grad_prim(4,3) + grad_prim(4,3)  - &
!                                                                         & 2.0_wp/3.0_wp*divergence))  

! 2nd approach: used

          n_visc_flux = normalviscousflux(vg(:,i_node,elem_id),phig(:,:,i_node,elem_id), normal(:), &
                                &   nequations,mut(i_node,elem_id))
          
          ! Add contribution of the viscous stress to the horizontal force
          h_force_face_nodes(i,j) = h_force_face_nodes(i,j) - n_visc_flux(2)
          

          ! Add contribution of the viscous stress to the vertical force
          v_force_face_nodes(i,j) = v_force_face_nodes(i,j) - n_visc_flux(3)

          ! Add contribution of the viscous stress to the side force
          s_force_face_nodes(i,j) = s_force_face_nodes(i,j) - n_visc_flux(4) 

        end do
      end do

      cnt = 0

      ! Integrate the pointwise forces over the face uisng a tensor product
      ! approach
      do j = 1, nodesperedge
        do i = 1, nodesperedge

          cnt = cnt + 1

          ! Volumetric node index corresponding to face and node on face indices
          i_node = kfacenodes(cnt,face_id)
         
          ! Facial index corresponding to face and node on face indices
          j_node = (face_id-1)*nodesperface + cnt
         
          ! Outward facing normal of facial node
          normal = -1.0_wp*Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          
          h_force = h_force + pmat_Soln(i)*pmat_Soln(j)*h_force_face_nodes(i,j)

          v_force = v_force + pmat_Soln(i)*pmat_Soln(j)*v_force_face_nodes(i,j)

          s_force = s_force + pmat_Soln(i)*pmat_Soln(j)*s_force_face_nodes(i,j)

!          ! Compute wet area
!          area(1) = area(1) + pmat_Soln(i)*pmat_Soln(j)*abs(normal(1))
!          area(2) = area(2) + pmat_Soln(i)*pmat_Soln(j)*abs(normal(2))
!          area(3) = area(3) + pmat_Soln(i)*pmat_Soln(j)*abs(normal(3))

        end do
      end do

    end do



    ! Combines values of the 3 forces from all the processes and distribute the
    ! result back to all processes
    call mpi_allreduce(h_force,h_force_sum,1, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)

    call mpi_allreduce(v_force,v_force_sum,1, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)
    
    call mpi_allreduce(s_force,s_force_sum,1, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)
        
!    write(*,*) 'Before MPI'
!    write(*,*) 'Area', area(1), area(2)
!    write(*,*)
!
!    call mpi_allreduce(area,area_sum,3, &
!      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,i_err)


    ! Aerodynamic coefficients
    c_1_aero_coeff = h_force_sum/(1.0_wp/2.0_wp*aero_coeffs_surface) 
    c_2_aero_coeff = v_force_sum/(1.0_wp/2.0_wp*aero_coeffs_surface) 
    c_3_aero_coeff = s_force_sum/(1.0_wp/2.0_wp*aero_coeffs_surface) 

!    write(*,*) 'After MPI'
!    write(*,*) 'Area', area_sum
!    write(*,*)
!
!    stop

    return
  end subroutine compute_aerodynamic_coefficients

  !============================================================================

  !============================================================================

end module aerodynamic_coefficients

