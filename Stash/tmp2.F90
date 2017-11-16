    allocate(xghst_LGL(ndim,nghost)) ; xghst_LGL = -100.0_wp ;

    call PetscComm1DDataSetup (xg,xghst_LGL,xpetsc,xlocpetsc,size(xg,1), size(xg,2),  &
                               nelems, nodesperproc, size(xghst_LGL,2))

    call UpdateComm1DGhostData(xg,xghst_LGL,xpetsc,xlocpetsc,size(xg,1), size(xg,2),  &
                                       nodesperproc, size(xghst_LGL,2))

  end subroutine PetscGridLocations_LGL

!============================================================================

  subroutine PetscComm2DDataSetup(Zin, Zghstin, Zpetscin, Zlocin, nq, nd, nk, ne, ngh)

    ! this routine allocates the ghost data for Navier Stokes computations

    use referencevariables,  only: ihelems, nfacesperelem, myprocid, nelems
    use variables,           only: ef2e, kfacenodes
    use initcollocation,     only: element_properties
    implicit none

    ! Arguments
    ! =========
    integer,  intent(in) :: nq, nd, nk, ne, ngh
    real(wp), intent(in) :: Zin(nq,nd,nk,ne)
    real(wp), intent(in) :: Zghstin(nq,nd,ngh)

    Vec Zpetscin
    Vec Zlocin

    PetscErrorCode ierrpetsc
    PetscScalar xinit

    integer :: nodesperface
    integer :: ntotu,ntotv,ntotG
    integer :: ielem, iloc, idir, iface
    integer :: i, kelem, kface, ieq
    integer, allocatable :: iyu(:)

    xinit = 0.0_wp

    ntotu = nq * nd * nk * nelems                        ! length of on process data including padding (for 1D solution vector)
    ntotv = nq * nd * np                                 ! length of data being used (excludes padding)
    ntotG = nq * nd * ngh                                ! length of ghost data on process

    ! allocate memory for ghost locations
    allocate(iyu(ntotG))

    iloc = 0                                             ! ghost data counter
    elem_loop:do ielem = ihelems(1), ihelems(2)          ! loop over elements

      face_loop:do iface = 1,nfacesperelem               ! loop over faces on element

        if (ef2e(3,iface,ielem) == myprocid) cycle       ! cycle if neighbor is ON process

        kelem = ef2e(2,iface,ielem)                      ! element of neighbor
        kface = ef2e(1,iface,ielem)                      ! face of neighbor

        call element_properties(kelem,              &    ! Get off-element properties
                       n_pts_2d=nodesperface,       & 
                     kfacenodes=kfacenodes)

        do i = 1, nodesperface                           ! loop over nodes on neighbor face

          do idir = 1,nd                                 ! loop over gradient directions

            do ieq = 1,nq                                ! loop over equations

              iloc = iloc+1                              ! advance position in ghost array
                                                         ! set position of ghost in global vector containing  entropy variable gradient data
              iyphi(iloc) = nq * nd * nk * (kelem-1)      &
                          + (kfacenodes(i,kface)-1)*nd*nq &
                          + nq*(idir-1)                   &
                          + ieq

            end do
          end do

        end do

      end do face_loop

    end do elem_loop

    ! use C indexing
    iyu = iyu-1

    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntotu, petsc_decide, ntotG, iyu, Zpetscin, ierrpetsc)
    ! initialize to zero
    call VecSet(Zpetscin, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(Zpetscin,ierrpetsc)
    call VecAssemblyEnd  (Zpetscin,ierrpetsc)

    deallocate(iyu)

    ! create container for local vector that contains on process plus ghost data for solution gradients
    call VecCreateSeq(petsc_comm_self, ntotu + ntotG, Zlocin, ierrpetsc)

  end subroutine PetscComm2DDataSetup

  !===========================================================================================
  !===========================================================================================
  !  End of Communication setup routines                                       ===============
  !===========================================================================================
  !===========================================================================================

  subroutine UpdateComm2DGhostData(Zin, Zghstin, Zpetscin, Zlocin, nq, nd, nk, ngh)

    ! this subroutine communicates solution data across processes and updates the array ughst.

    use referencevariables,   only: ihelems, nfacesperelem, myprocid, nelems
    use variables,            only: ef2e
    use initcollocation,      only: element_properties

    implicit none

    ! Arguments
    ! =========
    integer,  intent(in) :: nq, nd, nk, ngh
    real(wp), intent(in) :: Zin(nq,nd,nk,ihelems(1):ihelems(2))
    real(wp), intent(inout) :: Zghstin(nq,nd,ngh)

    Vec Zpetscin
    Vec Zlocin

    PetscErrorCode ierrpetsc

    integer :: ntotu, ntotv, ntotw, ntot
    integer :: ielem, iloc, inode, iface, nodesperface, nodesperelem, idir
    integer :: kelem
    integer :: i, ieq
    integer,  allocatable :: iyu(:)
    real(wp), allocatable ::  yu(:)

    real(wp), pointer :: xx_Z(:)

                                   ! length of arrays for filling global vectors with data
    ntotu = nq * nd * nk * nelem   ! length of on process data
    ntotv = nq * nd * np           ! length of on-process data (no padding, just data)
    ntotw = nq * nd * nk           ! length of incoming block of data WITH PADDING of zeros

    do ielem = ihelems(1), ihelems(2)                       ! loop over elements

      call element_properties(ielem, n_pts_3d=nodesperelem) ! size is element dependent

      ntot = nq * nd * nodesperelem                         ! bucket size
      if(allocated(iyu)) deallocate(iyu) ; allocate(iyu(ntot)) ; iyu = 0
      if(allocated( yu)) deallocate( yu) ; allocate( yu(ntot)) ;  yu = 0.0_wp

      do inode = 1, nodesperelem                            ! loop over nodes

        do idir = 1,nd                                      ! loop over directions
          do ieq = 1, nq                                    ! loop over equations
                                                            ! update gradient data
            yu(nq*nd*(inode-1)+nq*(idir-1)+ieq) = Zin(ieq,idir,inode,ielem)
                                                            ! update global location of gradient data in 1D vector
           iyu(nq*nd*(inode-1)+nq*(idir-1)+ieq) = ntotu*(ielem-1) &! shift over previous elements
                                                + nq*nd*(inode-1) &! shift over previous nodes
                                                + nq*(idir-1)     &! shift over previous directions
                                                + ieqn            &! shift over previous eqns
                                                - 1                ! shift from fortran to C indexing
          end do
        end do

      end do
                                                            ! set values in petsc vector
      call VecSetValues(Zpetscin,ntot,iyu,yu,insert_values,ierrpetsc)

    end do

    call VecAssemblyBegin(Zpetscin,ierrpetsc)               ! assemble Petsc vector
    call VecAssemblyEnd  (Zpetscin,ierrpetsc)
                                                            ! update ghost values
    call VecGhostUpdateBegin(Zpetscin, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd  (Zpetscin, insert_values, scatter_forward, ierrpetsc)

    call VecGhostGetLocalForm(Zpetscin, Zlocin, ierrpetsc)  ! get local data including ghost points

    call VecGetArrayF90(Zlocin, xx_Z, ierrpetsc)            ! use fortran pointer for convenience

    iloc = 0
    do ielem = ihelems(1), ihelems(2)                       ! loop over elements

      do iface = 1,nfacesperelem                            ! loop over faces

        if (ef2e(3,iface,ielem) == myprocid) cycle          ! cycle if neighbor is on process

        kelem = ef2e(2,iface,ielem)                         ! element of neighbor

        call element_properties(kelem, n_pts_2d=nodesperface)

        do i = 1, nodesperface

          iloc = iloc+1                                     ! update ghost node index

          do idir = 1, nd                                   ! loop over directions
            do ieq = 1,nq                                   ! loop over equations
              Zghstin(ieq,idir,iloc) = xx_Z(ntotu+nd*nq*(iloc-1)+nq*(idir-1)+ieq) ! fill ghost values
            end do
          end do

        end do

      end do

    end do

    call VecRestoreArrayF90(Zlocin,xx_Z,ierrpetsc)          ! release pointer
    call VecGhostRestoreLocalForm(Zpetscin,Zlocin,ierrpetsc)
    if(associated(xx_Z)) deallocate(xx_Z)

  end subroutine UpdateComm2DGhostData

