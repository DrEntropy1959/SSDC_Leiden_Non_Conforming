! This module contains the necessary routines to compute the error for the
! no-slip boundary conditions at the solid wall.


module error_heat_entropy_flow_wall_bc

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Subroutines in this module are generally private 
  private

  ! Exceptions, i.e. public subroutines or functions
  public compute_heat_entropy_flow_wall_bc_error

  contains

  !============================================================================

  !============================================================================
  ! compute_bc_no_slip_wall_error - Drives the computation of the no-slip wall
  ! boundary condition error 
  !
  subroutine compute_heat_entropy_flow_wall_bc_error()

    ! Load modules
    use variables
    use referencevariables
    use navierstokes, only : dVdW
    use nsereferencevariables, only : mu0
    use mpimod
    use controlvariables

    ! Nothing is implicitly defined
    implicit none

    integer :: iell, elem_high
    integer :: i_elem, i_node, j_node
    integer :: n_wall_elems
    integer :: i, j, elem_id, face_id

    !real(wp), dimension(2,nodesperedge,nodesperedge) :: l2_error_u_1, l2_error_u_2, l2_error_u_3
    !real(wp), dimension(2) :: l2_error_sum_u_1, l2_error_sum_u_2, l2_error_sum_u_3

    real(wp) :: linf_error_entropy_flow
    real(wp) :: linfmax_error_entropy_flow

    real(wp) :: diff

    integer :: cnt

    integer :: i_err

    real(wp), dimension(3) :: nx, unit_normal

    continue

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements whcih own a wall face 
    n_wall_elems = size(wall_elem_face_ids(1,:))

    ! Initialize to zero the errors
    linf_error_entropy_flow = 0.0_wp

 
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
          nx = Jx_r(i_node,elem_id)*facenodenormal(:,j_node,elem_id)

          ! Unit normal direction
          unit_normal = nx/Magnitude(nx)

          ! Compute error for the entropy flow
          diff = vg(5,i_node,elem_id)*dot_product(unit_normal(:),phig(5,:,i_node,elem_id)) - heat_entropy_flow_wall_bc
          linf_error_entropy_flow = max(linf_error_entropy_flow,abs(diff))

        end do
      end do
    end do

    ! Finish max norm
    ! ---------------
    call mpi_allreduce(linf_error_entropy_flow,linfmax_error_entropy_flow,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,i_err)

    linf_error_heat_entropy_flow_wall_bc = linfmax_error_entropy_flow


    return
  end subroutine compute_heat_entropy_flow_wall_bc_error

  !============================================================================

  !============================================================================

end module error_heat_entropy_flow_wall_bc
