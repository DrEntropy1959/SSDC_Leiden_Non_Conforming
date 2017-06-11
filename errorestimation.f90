module errorestimation
  ! this module contains the necessary routines to simulate the
  ! Navier-Stokes equations using the SSDC method. 
  use precision_vars
  implicit none

  private

  public calcembeddedspatialerror
  public calcembeddedtemporalerror

contains

  subroutine calcembeddedspatialerror()
    ! this subroutine calculates the embedded error approximation
    ! using the solution, ug, and the embedded solution uhat.
    use variables, only: gsat, Jx_r
    use collocationvariables, only: pvol
    use controlvariables
    use referencevariables
    use mpimod
    implicit none

    ! indices
    integer :: inode, ielem

    ! low and high volumetric element indices
    integer :: iell, ielh

    ! different error estimates
    real(wp) :: l2(2), linf
    real(wp) :: l2sum(2), linfmax
    ! local error
    real(wp), allocatable :: ex(:)
    integer :: ierr

    allocate(ex(nequations))

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! initialize errors to zero
    l2 = 0.0_wp
    linf = 0.0_wp
    ! loop over elements
    do ielem = iell, ielh
      ! loop over each index in the element
      do inode = 1, nodesperelem
        ! compute the local embedded error
        ex = gsat(:,inode,ielem)/Jx_r(inode,ielem)
        ! calculate linf contribution
        linf = max(linf,Jx_r(inode,ielem)*maxval(abs(ex)))
        ! calculate the integral contribution of l2 error. pvol*J is the volumetric
        ! integration weight.
        l2(1) = l2(1) + pvol(inode)*Jx_r(inode,ielem)*dot_product(ex,ex)
        l2(2) = l2(2) + pvol(inode)*Jx_r(inode,ielem)
      end do
    end do
    ! finish l2 error
    l2(1) = sqrt(l2(1))

    call mpi_allreduce(l2,l2sum,2, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)

    l2sum(1) = l2sum(1)/l2sum(2)

    call mpi_allreduce(linf,linfmax,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,ierr)

    !   if(verbose .and. myprocid == 0) write(*,101) l2sum(1), linfmax
    !   101 format('embedded l2 error: ',ES12.5,1X,'embedded linf error: ',ES12.5,1X)
    err_space_lf = linfmax
    err_space_l2 = l2sum(1)

    deallocate(ex)

  end subroutine calcembeddedspatialerror

  subroutine calcembeddedtemporalerror()
    ! this subroutine calculates the embedded error approximation
    ! using the solution, ug, and the embedded solution uhat.
    use variables, only: ug,uhat,Jx_r
    use collocationvariables, only: pvol
    use controlvariables
    use referencevariables
    use mpimod
    implicit none

    ! indices
    integer :: inode, ielem

    ! low and high volumetric element indices
    integer :: iell, ielh

    ! different error estimates
    real(wp) :: l2, linf
    real(wp) :: l2sum, linfmax
    ! local error
    real(wp), allocatable :: ex(:)
    integer :: ierr

    allocate(ex(nequations))

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! initialize errors to zero
    l2 = 0.0_wp
    linf = 0.0_wp
    ! loop over elements
    do ielem = iell, ielh
      ! loop over each index in the element
      do inode = 1, nodesperelem
        ! compute the local embedded error
        ex = ug(:,inode,ielem) - uhat(:,inode,ielem)
        ! calculate linf contribution
        linf = max(linf,maxval(abs(ex)))
        ! calculate the integral contribution of l2 error. pvol*J is the volumetric
        ! integration weight.
        l2 = l2 + pvol(inode)*Jx_r(inode,ielem)*dot_product(ex,ex)
      end do
    end do
    ! finish l2 error
    l2 = sqrt(l2)

    call mpi_allreduce(l2,l2sum,1, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)

    call mpi_allreduce(linf,linfmax,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,ierr)

    !   if(verbose .and. myprocid == 0) write(*,101) l2sum, linfmax
    !   101 format('l2 error: ',ES12.5,1X,'linf error: ',ES12.5,1X,'S: ',ES12.5)
    err_time_lf = linfmax
    err_time_l2 = l2sum

    deallocate(ex)

  end subroutine calcembeddedtemporalerror

end module errorestimation
