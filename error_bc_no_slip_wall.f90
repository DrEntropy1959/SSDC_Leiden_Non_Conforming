! This module contains the necessary routines to compute the error for the
! no-slip boundary conditions at the solid wall.


module error_bc_no_slip_wall

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Subroutines in this module are generally private 
  private

  ! Exceptions, i.e. public subroutines or functions
  public compute_bc_no_slip_wall_error

  contains

  !============================================================================

  !============================================================================
  ! compute_bc_no_slip_wall_error - Drives the computation of the no-slip wall
  ! boundary condition error 
  !
  subroutine compute_bc_no_slip_wall_error()

    ! Load modules
    use variables
    use referencevariables
    use navierstokes, only : dVdW
    use nsereferencevariables, only : mu0
    use mpimod
    use controlvariables

    ! Nothing is implicitly defined
    implicit none

    integer :: elem_low, elem_high
    integer :: i_elem, i_node
    integer :: n_wall_elems
    integer :: i, j, elem_id, face_id

    !real(wp), dimension(2,nodesperedge,nodesperedge) :: l2_error_u_1, l2_error_u_2, l2_error_u_3
    !real(wp), dimension(2) :: l2_error_sum_u_1, l2_error_sum_u_2, l2_error_sum_u_3

    real(wp) :: linf_error_u_1, linf_error_u_2, linf_error_u_3
    real(wp) :: linfmax_error_u_1, linfmax_error_u_2, linfmax_error_u_3

    real(wp) :: diff

    integer :: cnt

    integer :: i_err

    continue

    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Number of elements whcih own a wall face 
    n_wall_elems = size(wall_elem_face_ids(1,:))

    ! Initialize to zero the errors
    linf_error_u_1 = 0.0_wp
    linf_error_u_2 = 0.0_wp
    linf_error_u_3 = 0.0_wp

    cnt = 0


    do i_elem = 1, n_wall_elems
      elem_id = wall_elem_face_ids(1,i_elem)
      face_id = wall_elem_face_ids(2,i_elem)

      cnt = 0
      do j = 1, nodesperedge
       
        do i = 1, nodesperedge

          cnt = cnt + 1
          
          ! Volumetric node index corresponding to face and node on face indices
          i_node = kfacenodes(cnt,face_id)

          ! Compute error for the three velocity components
          diff = vg(2,i_node,elem_id) - 0.0_wp
          linf_error_u_1 = max(linf_error_u_1,abs(diff))

        
          diff = vg(3,i_node,elem_id) - 0.0_wp
          linf_error_u_2 = max(linf_error_u_2,abs(diff))
        
        
          diff = vg(4,i_node,elem_id) - 0.0_wp
          linf_error_u_3 = max(linf_error_u_3,abs(diff))

        end do
      end do
    end do

    ! Finish max norm
    ! ---------------
    call mpi_allreduce(linf_error_u_1,linfmax_error_u_1,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,i_err)

    call mpi_allreduce(linf_error_u_2,linfmax_error_u_2,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,i_err)

    call mpi_allreduce(linf_error_u_3,linfmax_error_u_3,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,i_err)

    linf_error_bc_no_slip_wall_u_1 = linfmax_error_u_1
    linf_error_bc_no_slip_wall_u_2 = linfmax_error_u_2
    linf_error_bc_no_slip_wall_u_3 = linfmax_error_u_3


    return
  end subroutine compute_bc_no_slip_wall_error

  !============================================================================

  !============================================================================

end module error_bc_no_slip_wall
