module mpimod

  use precision_vars
  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscis.h"
#include "finclude/petscvec.h90"

  integer :: mpi_real_wp, mpi_complex_wp
  integer :: mpi_default_wp

contains

  subroutine mpiInit()
    use referencevariables, only: nprocs, myprocid
    implicit none
    integer :: ierr

    ! initialize petsc library
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    if (wp==sp) then
      mpi_real_wp = mpi_real
      mpi_complex_wp = mpi_complex
      mpi_default_wp = mpi_real
    else
      mpi_real_wp = mpi_double_precision
      mpi_complex_wp = mpi_double_complex
      mpi_default_wp = mpi_double_precision
    end if

    ! get the total number of MPI processes in the specified 
    ! communicator. If the communicator is MPI_COMM_WORLD, then it 
    ! represents the number of MPI tasks available in the application. By 
    ! default PETSC_COMM_WORLD and MPI_COMM_WORLD are identical unless you 
    ! PETSc is run on ONLY a subset of MPI_COMM_WORLD.
    call mpi_comm_size(petsc_comm_world, nprocs, ierr)

    ! get the rank of the calling MPI process within the specified 
    ! communicator. Initially, each process will be assigned a unique 
    ! integer rank between 0 and number of tasks - 1 within the communicator
    ! PETSC_COMM_WORLD. 		
    call mpi_comm_rank(petsc_comm_world, myprocid, ierr)

  end subroutine mpiInit

  function checkNegTemp(tempflag)
    logical, intent(inout) :: tempflag

    integer :: checkNegTemp

    integer :: iconv, jconv, ierr

    iconv = 0
    if (tempflag) iconv = -1
    call mpi_allreduce(iconv,jconv,1,MPI_INTEGER,MPI_MIN,petsc_comm_world,ierr)
    checkNegTemp = jconv

  end function checkNegTemp

  subroutine distributeelements_cgns()
    ! this subroutine actually tells each process what elements
    ! it owns, which will be read in from the data file, and it
    ! specifies the connectivity of each element.
    use referencevariables
    use variables, only: e2v, elempart, jelems, iae2v, jae2v, ef2e, boundaryelems
    implicit none
    integer :: i, iproc, ierr
    integer :: j, ii, jj
    integer :: netmp, npetmp

    integer, allocatable :: melemsonproc(:), ef2etmp1(:,:,:), ef2etmp2(:,:,:)
    integer, allocatable :: i2jelems(:), j2ielems(:)
    integer, allocatable :: icnt(:)
    integer, allocatable :: e2vtmp(:,:)
    integer :: rtag, stag
    integer :: irsreq(2)
    integer, allocatable :: irsstatus(:,:)
    integer :: msgsize

    allocate(irsstatus(MPI_STATUS_SIZE,2))

    allocate(melemsonproc(0:2*nprocs-1))
    allocate(icnt(0:nprocs-1))

    ! process 0 (i.e. the master) orders elements for scattering
    if (myprocid == 0) then
      melemsonproc = 0
      ! count number of elements on each process and put it in
      ! a container as the upper index
      do i = 1, nelems
        melemsonproc(2*elempart(i)+1) = melemsonproc(2*elempart(i)+1)+1
      end do
      ! now make the upper index additive, such that all odd indices contain
      ! the upper element bound of a given process
      do i = 3,2*nprocs-1,2
        melemsonproc(i) = melemsonproc(i-2)+melemsonproc(i)
      end do
      ! the first process starts with element 1
      melemsonproc(0) = 1
      ! now determine the first element on the other processes
      do i = 2,2*(nprocs-1),2
        melemsonproc(i) = melemsonproc(i-1)+1
      end do

      allocate(i2jelems(1:nelems))
      allocate(j2ielems(1:nelems))
      ! we now determine the mapping from the old
      ! element indices to the new element indices
      icnt = 0
      do i = 1, nelems
        iproc = elempart(i)
        i2jelems(melemsonproc(2*iproc)+icnt(iproc)) = i
        j2ielems(i) = melemsonproc(2*iproc)+icnt(iproc)
        icnt(iproc) = icnt(iproc) + 1
      end do
    else
      allocate(i2jelems(0))
    end if

    ! create a barrier synchronization in the group. Each task, when 
    ! reaching the MPI_Barrier call, blocks until all tasks in the group 
    ! reach the same MPI_Barrier call. 
    call mpi_barrier(petsc_comm_world,ierr)

    ! the lower and upper ordered-element bounds for each process is communicated
    call mpi_scatter(melemsonproc, 2, mpi_integer, &
      ihelems, 2, mpi_integer, 0, PETSC_COMM_WORLD, ierr)

    ! the original-to-ordered-element mapping is communicated to each process
    allocate(jelems(ihelems(1):ihelems(2)))
    netmp = ihelems(2)-ihelems(1)+1
    call mpi_scatterv(i2jelems, icnt, melemsonproc(0:2*(nprocs-1):2)-1, mpi_integer, &
      jelems, netmp, mpi_integer, 0, PETSC_COMM_WORLD, ierr)

    
    ! the ordered-element to vertex connectivity is allocated
    allocate(e2v(2**ndim,ihelems(1):ihelems(2)))

    ! post non-blocking receives for each process
    ! to receive the element to vertex connectivity
    rtag = 100*nprocs+myprocid
    msgsize = netmp*(2**ndim)
    call mpi_irecv(e2v, msgsize, mpi_integer, 0, &
      rtag, PETSC_COMM_WORLD, irsreq(1), ierr)

    ! ordered-element face-to-face connectivity array is allocated
    allocate(ef2etmp2(3,2*ndim,ihelems(1):ihelems(2)))
    ef2etmp2 = 0

    ! post non-blocking receives for each process to
    ! receive the ordered-element face-to-face connectivity
    rtag = 200*nprocs+myprocid
    msgsize = netmp*2*ndim*3
    call mpi_irecv(ef2etmp2, msgsize, mpi_integer, 0, &
      rtag, PETSC_COMM_WORLD, irsreq(2), ierr)

    ! on process 0 (i.e. the master) use blocking sends
    if (myprocid == 0) then
      do iproc = 0, nprocs-1
        ! create temporary array for element-to-vertex connectivity on process iproc
        allocate( e2vtmp(2**ndim,melemsonproc(2*iproc):melemsonproc(2*iproc+1)) )
        ! create temporary array for face-to-face connectivity on process iproc
        allocate( ef2etmp1(3,2*ndim,melemsonproc(2*iproc):melemsonproc(2*iproc+1)) ) 
        ! number of elements on process iproc
        npetmp = melemsonproc(2*iproc+1)-melemsonproc(2*iproc)+1
        ! loop over elements on process iproc
        do i = melemsonproc(2*iproc), melemsonproc(2*iproc+1)
          ! the original element index
          ii = i2jelems(i)
          ! fill the e2v connectivity
          do j = 1,2**ndim
            jj = iae2v(ii)+j-1
            e2vtmp(j,i) = jae2v(jj)
          end do
          ! loop over faces
          do j = 1,2*ndim
            ! face of neighbor
            ef2etmp1(1,j,i) = ef2e(1,j,ii)
            ! original element index of neighbor
            jj = ef2e(2,j,ii)
            if (jj <= nelems) then
              ! neighbor element (with new element index)
              ef2etmp1(2,j,i) = j2ielems(jj)
              ! process of neighbor
              ef2etmp1(3,j,i) = elempart(jj)
            else
              ! self is neighbor (new element index)
              ef2etmp1(2,j,i) = i
              ! process of neighbor
              ef2etmp1(3,j,i) = iproc
            end if
          end do
        end do
        ! use blocking send to transfer e2v data
        stag = 100*nprocs+iproc
        msgsize = npetmp*(2**ndim)
        call mpi_send(e2vtmp(:,:), msgsize, mpi_integer, iproc, &
          stag, petsc_comm_world, ierr)
        ! use blocking send to transfer ef2e data
        stag = 200*nprocs+iproc
        msgsize = npetmp*2*ndim*3
        call mpi_send(ef2etmp1, msgsize, mpi_integer, iproc, &
          stag, petsc_comm_world, ierr)
        deallocate(e2vtmp)
        deallocate(ef2etmp1)
      end do
      deallocate(ef2e)
      deallocate(boundaryelems)
      deallocate(isectionvolume)
      deallocate(j2ielems)
    end if

    ! wait for receives to finish
    call MPI_Wait(irsreq(1), irsstatus(:,1), ierr)
    call MPI_Wait(irsreq(2), irsstatus(:,2), ierr)

    ! wait for other processes
    call mpi_barrier(PETSC_COMM_WORLD,ierr)

    allocate(ef2e(3,2*ndim,ihelems(1):ihelems(2)))
    ! assign global ef2e
    ef2e = ef2etmp2
    deallocate(ef2etmp2)

    deallocate(i2jelems)
    deallocate(icnt)
    deallocate(melemsonproc)
    deallocate(irsstatus)

  end subroutine distributeelements_cgns

  !============================================================================
  
  !============================================================================
  ! distribute_elements_aflr3 - Distributes elements and connectivity
  ! informations to each process.
  
  subroutine distribute_elements_aflr3()
    
    ! Load modules
    use referencevariables
    use variables, only: e2v, elempart, jelems, iae2v, jae2v, ef2e, vx_master, &
                       & vx, periodic_face_data_x1, n_periodic_faces_x1, &
                       & periodic_face_data_x2, n_periodic_faces_x2, &
                       & periodic_face_data_x3, n_periodic_faces_x3, &
                       & periodic_elem_face_ids_x1, periodic_elem_face_ids_x2, &
                       & periodic_elem_face_ids_x3, wall_face_data, &
                       & n_wall_faces, wall_elem_face_ids

    use collocationvariables, only : ldg_flip_flop_sign

    ! Nothing is implicitly defined
    implicit none

    integer :: i, i_proc, m
    integer :: i_elem, j, ii, jj, kk
    integer :: netmp, npetmp

    integer, allocatable :: melemsonproc(:), ef2etmp1(:,:,:), ef2etmp2(:,:,:)

    integer, allocatable :: ldg_flip_flop_sign_tmp1(:,:), ldg_flip_flop_sign_tmp2(:,:)

    integer, allocatable :: i2jelems(:), j2ielems(:)
    integer, allocatable :: icnt(:)
    integer, allocatable :: e2vtmp(:,:)
    integer :: s_tag, r_tag
    integer, allocatable :: i_s_r_status(:,:)
    integer :: m_size
    integer :: i_err

    integer :: max_n_vert_per_proc, n_vert_per_proc, j_vert_pos, j_vert_elem, &
             & k_vert_proc

    integer, allocatable, dimension(:) :: vert_list_proc

    integer :: i_p_face, p_elem, p_face, p_dir
    integer, allocatable, dimension(:) :: cnt_p_faces_x1, cnt_p_faces_x2, &
                                        & cnt_p_faces_x3
    integer, allocatable, dimension(:,:) :: p_data_x1, p_data_x2, p_data_x3

    integer :: i_w_face, w_elem, w_face
    integer, allocatable, dimension(:) :: cnt_w_faces
    integer, allocatable, dimension(:,:) :: w_data
    integer :: cnt_pack, cnt_unpack

    integer, allocatable, dimension(:) :: tmp_periodic_elem_face_ids_x1, &
                                          tmp_periodic_elem_face_ids_x2, &
                                          tmp_periodic_elem_face_ids_x3


    integer, allocatable, dimension(:) :: tmp_wall_elem_face_ids

    integer :: s_status(mpi_status_size)
    integer :: r_status(mpi_status_size)

    integer, allocatable, dimension(:) :: s_request_periodic_x1, &
                                          s_request_periodic_x2, &
                                          s_request_periodic_x3, &
                                          s_request_wall

    integer :: r_request_periodic_x1, r_request_periodic_x2, &
               r_request_periodic_x3, r_request_wall

    logical :: test_request_periodic_x1, test_request_periodic_x2, &
               test_request_periodic_x3, test_request_wall

    integer, allocatable, dimension(:) :: s_request_n_vertices, &
                                          s_request_vertices, & 
                                          s_request_e2v, &
                                          s_request_ef2e, &
                                          s_request_ldg_flip_flop_sign


    integer :: r_request_n_vertices, r_request_vertices, r_request_e2v, &
               r_request_ef2e, r_request_ldg_flip_flop_sign

    logical :: test_request_n_vertices, test_request_vertices, &
               test_request_e2v, test_request_ef2e, test_request_ldg_flip_flop_sign


    real(wp), allocatable, dimension(:) :: packed_vertices
    
    integer, allocatable, dimension(:) :: packed_e2v, packed_ef2e, packed_ldg_flip_flop_sign

    real(wp), allocatable, dimension(:) :: tmp_vx
    
    integer, allocatable, dimension(:) :: tmp_e2v
    
    integer, allocatable, dimension(:) :: tmp_ef2etmp2, &
                                        & tmp_ldg_flip_flop_sign_tmp2

    integer            :: n_elems

    integer, parameter :: qdim = 6   !  dimension of first array in ef2e see initgrid.F90 
                                     !  for definitions of each of the qdim 

    continue

    ! Number of vertices per element 
    nverticesperelem = 2**ndim

    ! Allocate memory for MPI
    allocate(i_s_r_status(mpi_status_size,4))

    ! Allocate memory for lower and upper index of the element ID own by each 
    ! process
    allocate(melemsonproc(0:2*nprocs-1))
    melemsonproc = 0

    ! Allocate memory for counter array
    allocate(icnt(0:nprocs-1))

    ! Allocate memory for counter array
    allocate(cnt_p_faces_x1(0:nprocs-1))
    allocate(cnt_p_faces_x2(0:nprocs-1))
    allocate(cnt_p_faces_x3(0:nprocs-1))

    ! Allocate memory for counter array
    allocate(cnt_w_faces(0:nprocs-1))


    ! Process 0 (i.e. the master process) orders the elements for scattering
    ! ==========================================================================
    if (myprocid == 0) then
      ! Count the number of elements on each process and put it in a container 
      ! as the upper index
      do i_elem = 1, nelems
        melemsonproc(2*elempart(i_elem) + 1) = melemsonproc(2*elempart(i_elem) + 1) + 1
      end do

      ! Make the upper index additive, such that all odd indices contain the
      ! upper element bound of a given process
      do i_proc = 3, 2*nprocs-1, 2
        melemsonproc(i_proc) = melemsonproc(i_proc-2) + melemsonproc(i_proc)
      end do

      ! The first process starts with element ID 1
      melemsonproc(0) = 1

      ! Determine the first element ID on the other processes
      do i_proc = 2, 2*(nprocs-1), 2
        melemsonproc(i_proc) = melemsonproc(i_proc-1) + 1
      end do

      ! Determine the mapping from the old element indices to the new element indices
      ! ========================================================================

      ! Allocate data for indices map
      allocate(i2jelems(1:nelems))
      i2jelems = 0
      
      allocate(j2ielems(1:nelems))
      j2ielems = 0

      ! Allocate memory for periodic face data, 1D arrays for mpi_isend
      allocate(p_data_x1(3*size(periodic_face_data_x1(1,:)),0:nprocs-1))
      p_data_x1 = 0
      
      allocate(p_data_x2(3*size(periodic_face_data_x2(1,:)),0:nprocs-1))
      p_data_x2 = 0
      
      allocate(p_data_x3(3*size(periodic_face_data_x3(1,:)),0:nprocs-1))
      p_data_x3 = 0

      ! Allocate memory for wall faces
      allocate(w_data(2*size(wall_face_data(1,:)),0:nprocs-1))
      w_data = 0
      
      ! Initialize counters to zero
      icnt = 0
      cnt_p_faces_x1 = 0
      cnt_p_faces_x2 = 0
      cnt_p_faces_x3 = 0
      cnt_w_faces = 0

      ! Loop over the elements
      do i_elem = 1, nelems
        ! Process ID
        i_proc = elempart(i_elem)

        ! Original numbering of the element for each processor
        i2jelems(melemsonproc(2*i_proc)+icnt(i_proc)) = i_elem
        
        ! New numbering of the elements for each processor
        j2ielems(i_elem) = melemsonproc(2*i_proc) + icnt(i_proc)
        
        ! Update processor counter
        icnt(i_proc) = icnt(i_proc) + 1

        ! Original ID of the element which owns a "periodic" face, x1 direction
        do i_p_face = 1, size(periodic_face_data_x1(1,:))
          
          ! Get the global ID of the element that owns a periodic boundary face
          p_elem = periodic_face_data_x1(1,i_p_face)

          ! Get the local (local for the element) ID of the periodic boundary face
          p_face = periodic_face_data_x1(2,i_p_face)

          ! Get the periodic direction of the face
          p_dir = periodic_face_data_x1(4,i_p_face) 


          if (p_elem == i_elem) then
            ! The i_elem owns a periodic boundary face
            cnt_p_faces_x1(i_proc) = cnt_p_faces_x1(i_proc) + 1
            p_data_x1(1+(cnt_p_faces_x1(i_proc)-1)*3,i_proc) = p_elem
            p_data_x1(2+(cnt_p_faces_x1(i_proc)-1)*3,i_proc) = p_face
            p_data_x1(3+(cnt_p_faces_x1(i_proc)-1)*3,i_proc) = p_dir
          end if
        end do

        ! Original ID of the element which owns a "periodic" face, x2 direction
        do i_p_face = 1, size(periodic_face_data_x2(1,:))
          
          ! Get the global ID of the element that owns a periodic boundary face
          p_elem = periodic_face_data_x2(1,i_p_face)

          ! Get the local (local for the element) ID of the periodic boundary 
          ! face
          p_face = periodic_face_data_x2(2,i_p_face)

          ! Get the periodic direction of the face
          p_dir = periodic_face_data_x2(4,i_p_face) 

          if (p_elem == i_elem) then
            ! The i_elem owns a periodic boundary face
            cnt_p_faces_x2(i_proc) = cnt_p_faces_x2(i_proc) + 1
            p_data_x2(1+(cnt_p_faces_x2(i_proc)-1)*3,i_proc) = p_elem
            p_data_x2(2+(cnt_p_faces_x2(i_proc)-1)*3,i_proc) = p_face
            p_data_x2(3+(cnt_p_faces_x2(i_proc)-1)*3,i_proc) = p_dir
          end if
        end do

        ! Original ID of the element which owns a "periodic" face, x3 direction
        do i_p_face = 1, size(periodic_face_data_x3(1,:))
          
          ! Get the global ID of the element that owns a periodic boundary face
          p_elem = periodic_face_data_x3(1,i_p_face)

          ! Get the local (local for the element) ID of the periodic boundary 
          ! face
          p_face = periodic_face_data_x3(2,i_p_face)

          ! Get the periodic direction of the face
          p_dir = periodic_face_data_x3(4,i_p_face) 

          if (p_elem == i_elem) then
            ! The i_elem owns a periodic boundary face
            cnt_p_faces_x3(i_proc) = cnt_p_faces_x3(i_proc) + 1
            p_data_x3(1+(cnt_p_faces_x3(i_proc)-1)*3,i_proc) = p_elem
            p_data_x3(2+(cnt_p_faces_x3(i_proc)-1)*3,i_proc) = p_face
            p_data_x3(3+(cnt_p_faces_x3(i_proc)-1)*3,i_proc) = p_dir
          end if
        end do

        ! Original ID of the element which owns a "wall" face
        do i_w_face = 1, size(wall_face_data(1,:))
          
          ! Get the global ID of the element that owns a wall boundary face
          w_elem = wall_face_data(1,i_w_face)

          ! Get the local (local for the element) ID of the wall boundary face
          w_face = wall_face_data(2,i_w_face)

          if (w_elem == i_elem) then
            ! The i_elem owns a wall boundary face
            cnt_w_faces(i_proc) = cnt_w_faces(i_proc) + 1
            w_data(1+(cnt_w_faces(i_proc)-1)*2,i_proc) = j2ielems(i_elem)
            w_data(2+(cnt_w_faces(i_proc)-1)*2,i_proc) = w_face
          end if
        end do

      end do

      ! Deallocate memory
      deallocate(periodic_face_data_x1,periodic_face_data_x2,periodic_face_data_x3)
      deallocate(wall_face_data)
    
    else
      allocate(i2jelems(0))
    end if

    ! Create a barrier synchronization in the group. Each task, when reaching 
    ! the MPI_Barrier call, blocks until all tasks in the group reach the same 
    ! MPI_Barrier call. This forces all the other processes to wait for the
    ! master process. 
    call mpi_barrier(petsc_comm_world,i_err)

    ! The lower and upper ordered-element bounds for each process is 
    ! communicated
    call mpi_scatter(melemsonproc,2,mpi_integer,ihelems,2,mpi_integer,0, &
      & petsc_comm_world,i_err)

    ! Set number of elements on each process 
    nelems = ihelems(2) - ihelems(1) + 1

    ! The original-to-ordered-element mapping is communicated to each process
    ! Each process will receive the global ID of the elements (from the original 
    ! ugrid). Such IDs are stored in jelems.
    allocate(jelems(ihelems(1):ihelems(2)))
    netmp = ihelems(2) - ihelems(1) + 1
    call mpi_scatterv(i2jelems,icnt,melemsonproc(0:2*(nprocs-1):2)-1, &
      & mpi_integer,jelems,netmp,mpi_integer,0,petsc_comm_world,i_err)

    ! Create a barrier synchronization in the group. 
    call mpi_barrier(petsc_comm_world,i_err)


    ! ==========================================================================
    ! Send the number of periodic and wall faces owned by each process
    ! ==========================================================================
    if (myprocid == 0) then

      ! Allocate memory for send requests and send status
      allocate(s_request_periodic_x1(nprocs))
      allocate(s_request_periodic_x2(nprocs))
      allocate(s_request_periodic_x3(nprocs))
      allocate(s_request_wall(nprocs))

      ! Special treatment for the master node (i.e., myprocid=0). There is no
      ! need to send data. We already know their values. Hence we just set them.
      ! ========================================================================
      ! Number of periodic faces in the x1 direction
      n_periodic_faces_x1 = cnt_p_faces_x1(myprocid)
      
      ! Number of periodic faces in the x2 direction
      n_periodic_faces_x2 = cnt_p_faces_x2(myprocid)
      
      ! Number of periodic faces in the x3 direction
      n_periodic_faces_x3 = cnt_p_faces_x3(myprocid)

      ! Number of wall faces
      n_wall_faces = cnt_w_faces(myprocid)

      ! General treatment for all the processes with ID different from zero
      ! ========================================================================
      ! Loop over the number of processes
      do i_proc = 1, nprocs-1

        ! Transfer the number of periodic faces in the x1 direction
        ! ======================================================================
        s_tag = 301*nprocs + i_proc
        m_size = 1
        call mpi_isend(cnt_p_faces_x1(i_proc),m_size,mpi_integer,i_proc, &
          & s_tag,petsc_comm_world,s_request_periodic_x1(i_proc),i_err)
        
        call mpi_wait(s_request_periodic_x1(i_proc),s_status,i_err)


        ! Transfer the number of periodic faces in the x2 direction
        ! ======================================================================
        s_tag = 302*nprocs + i_proc
        m_size = 1
        call mpi_isend(cnt_p_faces_x2(i_proc),m_size,mpi_integer,i_proc, &
          & s_tag,petsc_comm_world,s_request_periodic_x2(i_proc),i_err)

        call mpi_wait(s_request_periodic_x2(i_proc),s_status,i_err)


        ! Transfer the number of periodic faces in the x3 direction
        ! ======================================================================
        s_tag = 303*nprocs + i_proc
        m_size = 1
        call mpi_isend(cnt_p_faces_x3(i_proc),m_size,mpi_integer,i_proc, &
          & s_tag,petsc_comm_world,s_request_periodic_x3(i_proc),i_err)

        call mpi_wait(s_request_periodic_x3(i_proc),s_status,i_err)

        
        ! Transfer the number of wall faces
        ! ======================================================================
        s_tag = 1000*nprocs + i_proc
        m_size = 1
        call mpi_send(cnt_w_faces(i_proc),m_size,mpi_integer,i_proc, &
          & s_tag,petsc_comm_world,s_request_wall(i_proc),i_err)

        call mpi_wait(s_request_wall(i_proc),s_status,i_err)


      end do ! End do loop over the number of processes

    end if ! End if myprocid == 0

    ! Create a barrier synchronization in the group. 
    call mpi_barrier(petsc_comm_world,i_err)


    ! ==========================================================================
    ! The number of periodic and/or wall faces are received by the processes
    ! ==========================================================================
    if (myprocid > 0) then
      ! Receive number of periodic faces in x1 direction
      ! This number is used to allocate the memory for the periodic face data
      ! ========================================================================
      r_tag = 301*nprocs + myprocid
      m_size = 1
      call mpi_irecv(n_periodic_faces_x1,m_size,mpi_integer,0,r_tag, &
        & petsc_comm_world,r_request_periodic_x1,i_err)
      
      call mpi_wait(r_request_periodic_x1,r_status,i_err)
      call mpi_test(r_request_periodic_x1,test_request_periodic_x1, &
        & r_status,i_err) 

      ! Receive number of periodic faces in x2 direction
      ! This number is used to allocate the memory for the periodic face data
      ! ========================================================================
      r_tag = 302*nprocs + myprocid
      m_size = 1
      call mpi_irecv(n_periodic_faces_x2,m_size,mpi_integer,0,r_tag, &
        & petsc_comm_world,r_request_periodic_x2,i_err)

      call mpi_wait(r_request_periodic_x2,r_status,i_err)
      call mpi_test(r_request_periodic_x2,test_request_periodic_x2, &
        & r_status,i_err) 

      ! Receive number of periodic faces in x3 direction
      ! This number is used to allocate the memory for the periodic face data
      ! ========================================================================
      r_tag = 303*nprocs + myprocid
      m_size = 1
      call mpi_irecv(n_periodic_faces_x3,m_size,mpi_integer,0,r_tag, &
        & petsc_comm_world,r_request_periodic_x3,i_err)

      call mpi_wait(r_request_periodic_x3,r_status,i_err)
      call mpi_test(r_request_periodic_x3,test_request_periodic_x3, &
        & r_status,i_err) 
      
      ! Receive number of wall faces
      ! This number is used to allocate the memory for the wall face data
      ! ========================================================================
      r_tag = 1000*nprocs + myprocid
      m_size = 1
      call mpi_irecv(n_wall_faces,m_size,mpi_integer,0,r_tag, &
        & petsc_comm_world,r_request_wall,i_err)

      call mpi_wait(r_request_wall,r_status,i_err)
      call mpi_test(r_request_wall,test_request_wall,r_status,i_err) 

    end if ! End if myprocid > 0


    ! Create a barrier synchronization in the group. 
    call mpi_barrier(petsc_comm_world,i_err)

    ! ==========================================================================
    ! At this point if a process owns periodic and/or baoundary face we send the
    ! data with the following mpi calls.
    ! ==========================================================================
    if (myprocid == 0) then

      ! Special treatment for the master node (i.e., myprocid=0). There is no
      ! need to send data. We already know their values. Hence we just set them.
      ! ========================================================================
      ! p_data_x1 has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array periodic_elem_face_ids_x1
      ! ========================================================================     
      ! Allocate memory for the periodic face data in the x1 direction
      allocate(periodic_elem_face_ids_x1(3,n_periodic_faces_x1))

      ! Unpack p_data_x1
      if (n_periodic_faces_x1 /= 0) then
        periodic_elem_face_ids_x1(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_periodic_faces_x1
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 3
            periodic_elem_face_ids_x1(jj,ii) = p_data_x1(jj+(cnt_unpack-1)*3,myprocid)
          end do
        end do
      end if

      ! p_data_x2 has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array periodic_elem_face_ids_x2
      ! ========================================================================     
      ! Allocate memory for the periodic face data in the x2 direction
      allocate(periodic_elem_face_ids_x2(3,n_periodic_faces_x2))
     
      ! Unpack 1D array p_data_x2
      if (n_periodic_faces_x2 /= 0) then
        periodic_elem_face_ids_x2(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_periodic_faces_x2
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 3
            periodic_elem_face_ids_x2(jj,ii) = p_data_x2(jj+(cnt_unpack-1)*3,myprocid)
          end do
        end do
      end if

      ! p_data_x3 has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array periodic_elem_face_ids_x3
      ! ========================================================================     
      ! Allocate memory for the periodic face data in the x3 direction
      allocate(periodic_elem_face_ids_x3(3,n_periodic_faces_x3))
      
      ! Unpack 1D array p_data_x3
      if (n_periodic_faces_x3 /= 0) then
        periodic_elem_face_ids_x3(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_periodic_faces_x3
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 3
            periodic_elem_face_ids_x3(jj,ii) = p_data_x3(jj+(cnt_unpack-1)*3,myprocid)
          end do
        end do
      end if

      ! w_data has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array wall_elem_face_ids
      ! ========================================================================     
      ! Allocate memory for the wall face data
      allocate(wall_elem_face_ids(2,n_wall_faces))
      
      ! Unpack 1D array w_data
      if (n_wall_faces /= 0) then
        wall_elem_face_ids(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_wall_faces
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 2
           wall_elem_face_ids(jj,ii) = w_data(jj+(cnt_unpack-1)*2,myprocid)
          end do
        end do
      end if

      ! General treatment for all the processes with ID different from zero
      ! ========================================================================
      ! Loop over the number of processes
      do i_proc = 1, nprocs-1

        ! Transfer p_data_x1 to the process
        ! ======================================================================
        if (cnt_p_faces_x1(i_proc) /= 0) then
          s_tag = 401*nprocs + i_proc
          m_size = cnt_p_faces_x1(i_proc)*3 !size(p_data_x1(:,i_proc))
          call mpi_isend(p_data_x1(1:m_size,i_proc),m_size,mpi_integer,i_proc, &
            & s_tag,petsc_comm_world,s_request_periodic_x1(i_proc),i_err)
        
          call mpi_wait(s_request_periodic_x1(i_proc),s_status,i_err)
        end if

        ! Transfer p_data_x2 to the process
        ! ======================================================================
        if (cnt_p_faces_x2(i_proc) /= 0) then
          s_tag = 402*nprocs + i_proc
          m_size = cnt_p_faces_x2(i_proc)*3
          call mpi_isend(p_data_x2(:,i_proc),m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_periodic_x2(i_proc),i_err)
        
          call mpi_wait(s_request_periodic_x2(i_proc),s_status,i_err)
        end if

        ! Transfer p_data_x3 to the process
        ! ======================================================================
        if (cnt_p_faces_x3(i_proc) /= 0) then
          s_tag = 403*nprocs + i_proc
          m_size = cnt_p_faces_x3(i_proc)*3
          call mpi_send(p_data_x3(:,i_proc),m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_periodic_x3(i_proc),i_err)
        
          call mpi_wait(s_request_periodic_x3(i_proc),s_status,i_err)
        end if

        ! Transfer w_data to the process
        ! ======================================================================
        if (cnt_w_faces(i_proc) /= 0) then
          s_tag = 1100*nprocs + i_proc
          m_size = 2*cnt_w_faces(i_proc)
          call mpi_send(w_data(:,i_proc),m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_wall(i_proc),i_err)
        
          call mpi_wait(s_request_wall(i_proc),s_status,i_err)
        end if
     
      end do ! End loop over the processes
  
      ! Deallocate memory
      deallocate(p_data_x1, p_data_x2, p_data_x3)
      deallocate(cnt_p_faces_x1, cnt_p_faces_x2, cnt_p_faces_x3)
      deallocate(w_data)
      deallocate(cnt_w_faces)
      deallocate(s_request_periodic_x1,s_request_periodic_x2, s_request_periodic_x3)

    end if ! End if myprocid == 0

    
    ! Each process (except the master node) allocates and receive periodic 
    ! and/or wall faces data
    ! ==========================================================================
    if (myprocid > 0) then
      ! p_data_x1 has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array periodic_elem_face_ids_x1
      ! ========================================================================     
      ! Allocate memory for the periodic face data in the x1 direction
      allocate(periodic_elem_face_ids_x1(3,n_periodic_faces_x1))
      
      if (n_periodic_faces_x1 /= 0) then
        ! Allocate memory for the 1D temporary array
        allocate(tmp_periodic_elem_face_ids_x1(3*n_periodic_faces_x1))
       
        ! Receive 1D temporary array
        r_tag = 401*nprocs + myprocid
        m_size = 3*n_periodic_faces_x1
        call mpi_irecv(tmp_periodic_elem_face_ids_x1,m_size,mpi_integer,0,r_tag, &
          & petsc_comm_world,r_request_periodic_x1,i_err)

        call mpi_wait(r_request_periodic_x1,r_status,i_err)
        test_request_periodic_x1 = .false.
        call mpi_test(r_request_periodic_x1,test_request_periodic_x1, &
          & r_status,i_err)

        ! Unpack 1D temporary array into a 2D permanent array
        periodic_elem_face_ids_x1(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_periodic_faces_x1
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 3
            periodic_elem_face_ids_x1(jj,ii) = tmp_periodic_elem_face_ids_x1(jj+(cnt_unpack-1)*3)
          end do
        end do

        deallocate(tmp_periodic_elem_face_ids_x1)
      end if

      ! p_data_x2 has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array periodic_elem_face_ids_x2
      ! ========================================================================     
      ! Allocate memory for the periodic face data in the x2 direction
      allocate(periodic_elem_face_ids_x2(3,n_periodic_faces_x2))
      
      if (n_periodic_faces_x2 /= 0) then 
        ! Allocate memory for the 1D temporary array
        allocate(tmp_periodic_elem_face_ids_x2(3*n_periodic_faces_x2))

        ! Receive 1D temporary array
        r_tag = 402*nprocs + myprocid
        m_size = 3*n_periodic_faces_x2
        call mpi_irecv(tmp_periodic_elem_face_ids_x2,m_size,mpi_integer,0,r_tag, &
          & petsc_comm_world,r_request_periodic_x2,i_err)

        call mpi_wait(r_request_periodic_x2,r_status,i_err)
        test_request_periodic_x2 = .false.
        call mpi_test(r_request_periodic_x2,test_request_periodic_x2, &
          & r_status,i_err)

        ! Unpack 1D temporary array into a 2D permanent array
        periodic_elem_face_ids_x2(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_periodic_faces_x2
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 3
            periodic_elem_face_ids_x2(jj,ii) = tmp_periodic_elem_face_ids_x2(jj+(cnt_unpack-1)*3)
          end do
        end do

        deallocate(tmp_periodic_elem_face_ids_x2)

      end if

      ! p_data_x3 has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array periodic_elem_face_ids_x3
      ! ========================================================================     
      ! Allocate memory for the periodic face data in the x3 direction
      allocate(periodic_elem_face_ids_x3(3,n_periodic_faces_x3))
      
      if (n_periodic_faces_x3 /= 0) then 
        ! Allocate memory for 1D temporary array
        allocate(tmp_periodic_elem_face_ids_x3(3*n_periodic_faces_x3))

        ! Receive 1D temporary array
        r_tag = 403*nprocs + myprocid
        m_size = 3*n_periodic_faces_x3
        call mpi_irecv(tmp_periodic_elem_face_ids_x3,m_size,mpi_integer,0,r_tag, &
          & petsc_comm_world,r_request_periodic_x3,i_err)

        call mpi_wait(r_request_periodic_x3,r_status,i_err)
        test_request_periodic_x3 = .false.
        call mpi_test(r_request_periodic_x3,test_request_periodic_x3, &
          & r_status,i_err)

        ! Unpack 1D temporary array into a 2D permanent array
        periodic_elem_face_ids_x3(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_periodic_faces_x3
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 3
            periodic_elem_face_ids_x3(jj,ii) = tmp_periodic_elem_face_ids_x3(jj+(cnt_unpack-1)*3)
          end do
        end do

        deallocate(tmp_periodic_elem_face_ids_x3)
      end if

      ! w_data has been "packed". Now we need to receive and unpack it; and 
      ! set the 2D array wall_elem_face_ids
      ! ========================================================================     
      ! Allocate memory for the wall face data
      allocate(wall_elem_face_ids(2,n_wall_faces))

      if (n_wall_faces /= 0) then
        ! Allocate memory for 1D temporary array
        allocate(tmp_wall_elem_face_ids(2*n_wall_faces))

        ! Receive 1D temporary array
        r_tag = 1100*nprocs + myprocid
        m_size = 2*n_wall_faces
        call mpi_irecv(tmp_wall_elem_face_ids,m_size,mpi_integer,0,r_tag, &
          & petsc_comm_world,r_request_wall,i_err)

        call mpi_wait(r_request_wall,r_status,i_err)
        test_request_wall = .false.
        call mpi_test(r_request_wall,test_request_wall, &
          & r_status,i_err)

        ! Unpack 1D temporary array into a 2D permanent array
        wall_elem_face_ids(:,:) = 0
        cnt_unpack = 0
        do ii = 1, n_wall_faces
          cnt_unpack = cnt_unpack + 1
          do jj = 1, 2
            wall_elem_face_ids(jj,ii) = tmp_wall_elem_face_ids(jj+(cnt_unpack-1)*2)
          end do
        end do

        deallocate(tmp_wall_elem_face_ids)

      end if

    end if ! End if myprocid > 0 

    ! Create a barrier synchronization in the group. 
    call mpi_barrier(petsc_comm_world,i_err)

    ! =========================================================================
    ! Master node sends to each process the following information:
    !
    ! 1) e2v, connectivity using the original global node numbering (from the 
    ! serial read)
    !
    ! 2) ef2e, connectivity in the new global numbering after metis
    ! has partitioned the grid
    !
    ! 3) coordinates of the vertex (or nodes)
    ! =========================================================================

    if (myprocid == 0) then
      ! Allocate memory for the requests
      allocate(s_request_n_vertices(nprocs))
      allocate(s_request_vertices(nprocs))
      allocate(s_request_e2v(nprocs))
      allocate(s_request_ef2e(nprocs))
      allocate(s_request_ldg_flip_flop_sign(nprocs))
      
      do i_proc = 0, nprocs-1
        ! Allocate memory for temporary array for element-to-vertex connectivity on the process
        allocate(e2vtmp(2**ndim,melemsonproc(2*i_proc):melemsonproc(2*i_proc+1)))
        
        ! Allocate memory for temporary array for face-to-face connectivity on the process
        allocate(ef2etmp1(qdim,2*ndim,melemsonproc(2*i_proc):melemsonproc(2*i_proc+1))) 

        ! Allocate memory for temporary array for the LDG flip-flop sign on the process
        allocate(ldg_flip_flop_sign_tmp1(2*ndim,melemsonproc(2*i_proc):melemsonproc(2*i_proc+1))) 
        
        ! Number of elements on the process
        npetmp = melemsonproc(2*i_proc+1)-melemsonproc(2*i_proc)+1
        
        ! Number of vertices per process. This is a conservative estimate 
        ! because most of the nodes are shared among the element excpet the 
        ! boundary nodes.
        max_n_vert_per_proc = nverticesperelem*(melemsonproc(2*i_proc+1) - &
          & melemsonproc(2*i_proc) + 1)

        ! Construct the local e2v array (local for each process) 
        allocate(vert_list_proc(max_n_vert_per_proc))
        vert_list_proc = max_n_vert_per_proc + 100
    
        ! Initialize number of vertices per process
        n_vert_per_proc = 0

        ! Vertex position in the array
        j_vert_pos = 1
        
        cnt_pack = 0 

        ! Loop over all elements on process
        do i = melemsonproc(2*i_proc), melemsonproc(2*i_proc+1)
         
          cnt_pack = cnt_pack + 1
          
          ! Original element index
          ii = i2jelems(i)

          ! Fill the e2v connectivity
          do j = 1, 2**ndim
            jj = iae2v(ii) + j - 1
            e2vtmp(j,i) = jae2v(jj)
          end do
          
          ! Loop over all faces
          do j = 1, 2*ndim
            ! Face of neighbor
            ef2etmp1(1,j,i) = ef2e(1,j,ii)
            ef2etmp1(4,j,i) = ef2e(4,j,ii)
            ef2etmp1(5,j,i) = ef2e(5,j,ii)
            ef2etmp1(6,j,i) = ef2e(6,j,ii)
            
            ! Original element index of neighbor
            jj = ef2e(2,j,ii)
            
            ! Internal faces have ef2e(2,j,ii)>0
            if (jj > 0) then
              ! Neighbor element (with new element index)
              ef2etmp1(2,j,i) = j2ielems(jj)
              
              ! Process of neighbor
              ef2etmp1(3,j,i) = elempart(jj)
            
            ! Boundary face
            else
              ! Self is neighbor (new element index)
              ef2etmp1(2,j,i) = i
              
              ! Process of neighbor
              ef2etmp1(3,j,i) = i_proc
            end if
          end do

          ! Loop over all faces
          do j = 1, 2*ndim
            ldg_flip_flop_sign_tmp1(j,i) = ldg_flip_flop_sign(j,ii)
          end do


          ! Loop over all the vertices of the element
          v_loop_m: do j_vert_elem = 1, nverticesperelem
            
            m = e2vtmp(j_vert_elem,i)
            
            ! Loop over the list of known vertices
            do k_vert_proc = 1, n_vert_per_proc
              
              if (m == vert_list_proc(k_vert_proc)) then
              
                cycle v_loop_m
              
              end if
            
            end do
            
            ! Add vertex to the list
            vert_list_proc(j_vert_pos) = m
            
            ! Update position and counter
            j_vert_pos = j_vert_pos + 1
            n_vert_per_proc = n_vert_per_proc + 1
          
          end do v_loop_m

        end do ! En do loop over the element

        ! Special treatment for the master node. There is no need to use mpi
        ! since the data are already available in memory.
        ! =====================================================================
        if (i_proc == 0) then
          ! Set number of vertices
          nvertices = n_vert_per_proc

          ! Set vertex coordinates
          allocate(vx(3,nvertices))
          do ii = 1, size(vert_list_proc(1:n_vert_per_proc))
            do jj = 1, 3
              vx(jj,ii) = vx_master(jj,vert_list_proc(ii))
            end do
          end do

          ! Element-to-vertex connectivity (e2v) 
          allocate(e2v(2**ndim,ihelems(1):ihelems(2)))
          e2v = e2vtmp

          ! Element-face-to-element connectivity (ef2e)
          allocate(ef2etmp2(qdim,2*ndim,ihelems(1):ihelems(2)))
          ef2etmp2 = ef2etmp1

          ! Sign of the LDG flip-flop
          allocate(ldg_flip_flop_sign_tmp2(2*ndim,ihelems(1):ihelems(2)))
          ldg_flip_flop_sign_tmp2 = ldg_flip_flop_sign_tmp1
          
          ! Deallocate memory
          deallocate(e2vtmp)
          deallocate(ef2etmp1)
          deallocate(vert_list_proc)
          deallocate(ldg_flip_flop_sign_tmp1)

        else 
          ! General procedure for all the other processes. 
          ! ===================================================================
          ! Transfer number of vertices
          ! ===================================================================
          s_tag = 25*nprocs + i_proc
          m_size = 1
          call mpi_isend(n_vert_per_proc,m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_n_vertices(i_proc),i_err)
        
          call mpi_wait(s_request_n_vertices(i_proc),s_status,i_err)

          ! Transfer vertex coordinates
          ! ===================================================================
          ! Prepare 1D vector
          allocate(packed_vertices(3*n_vert_per_proc))
          cnt_pack = 0
          do ii = 1, n_vert_per_proc
            do jj = 1, 3
              cnt_pack = cnt_pack + 1
              packed_vertices(cnt_pack) = vx_master(jj,vert_list_proc(ii))
            end do
          end do

          ! Transfer packed array
          s_tag = 5001*nprocs + i_proc
          m_size = 3*n_vert_per_proc
          call mpi_isend(packed_vertices,m_size,mpi_double,i_proc,s_tag, &
            & petsc_comm_world,s_request_vertices(i_proc),i_err)
        
          call mpi_wait(s_request_vertices(i_proc),s_status,i_err)

          ! Transfer e2v data
          ! ===================================================================
          ! Prepare 1D vector
          allocate(packed_e2v(npetmp*(2**ndim)))
          cnt_pack = 0
          do ii = melemsonproc(2*i_proc), melemsonproc(2*i_proc+1)
            do jj = 1, 2**ndim
              cnt_pack = cnt_pack + 1
              packed_e2v(cnt_pack) = e2vtmp(jj,ii)
            end do
          end do

          ! Send 1D array
          s_tag = 100*nprocs + i_proc
          m_size = npetmp*(2**ndim)
          call mpi_isend(packed_e2v,m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_e2v(i_proc),i_err)

          call mpi_wait(s_request_e2v(i_proc),s_status,i_err)

          ! Transfer ef2e data
          ! ===================================================================
          ! Prepare 1D vector
          allocate(packed_ef2e(qdim*npetmp*(2*ndim)))
          cnt_pack = 0
          do ii = melemsonproc(2*i_proc), melemsonproc(2*i_proc+1)
            do jj = 1, 2*ndim
              do kk = 1, qdim
                cnt_pack = cnt_pack + 1
                packed_ef2e(cnt_pack) = ef2etmp1(kk,jj,ii)
              end do
            end do
          end do

          ! Send 1D array
          s_tag = 200*nprocs + i_proc
          m_size = npetmp*2*ndim*qdim
          call mpi_send(packed_ef2e,m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_ef2e(i_proc),i_err)

          call mpi_wait(s_request_ef2e(i_proc),s_status,i_err)

          ! Transfer ldg_flip_flop_sign
          ! ===================================================================
          ! Prepare 1D vector
          allocate(packed_ldg_flip_flop_sign(npetmp*(2*ndim)))
          cnt_pack = 0
          do ii = melemsonproc(2*i_proc), melemsonproc(2*i_proc+1)
            do jj = 1, 2*ndim
              cnt_pack = cnt_pack + 1
              packed_ldg_flip_flop_sign(cnt_pack) = ldg_flip_flop_sign_tmp1(jj,ii)
            end do
          end do

          
          ! Send 1D array
          s_tag = 327*nprocs + i_proc
          m_size = npetmp*2*ndim
          call mpi_send(packed_ldg_flip_flop_sign,m_size,mpi_integer,i_proc,s_tag, &
            & petsc_comm_world,s_request_ldg_flip_flop_sign(i_proc),i_err)

          call mpi_wait(s_request_ldg_flip_flop_sign(i_proc),s_status,i_err)

          ! Deallocate memory
          deallocate(e2vtmp)
          deallocate(ef2etmp1)
          deallocate(vert_list_proc)
          deallocate(packed_vertices)
          deallocate(packed_e2v)
          deallocate(packed_ef2e)
          deallocate(ldg_flip_flop_sign_tmp1)
          deallocate(packed_ldg_flip_flop_sign)
        
        end if ! End if not master node

      end do ! End do processors
 
      ! Deallocate memory
      deallocate(ef2e)
      deallocate(j2ielems)
      deallocate(vx_master)
      deallocate(s_request_n_vertices)
      deallocate(s_request_vertices)
      deallocate(s_request_e2v)
      deallocate(s_request_ef2e)
      deallocate(ldg_flip_flop_sign)
      deallocate(s_request_ldg_flip_flop_sign)
    
    end if ! End if myprocid == 0


    ! All the processes receive the number of vertices that they own
    ! These valus are used to allocate memory for the vertex coordinates
    ! =========================================================================
    if (myprocid > 0) then
      ! Receive number of vertices. This number is used to allocate the memory for
      ! the vertex coordinates array.
      r_tag = 25*nprocs + myprocid
      m_size = 1
      call mpi_irecv(nvertices,m_size,mpi_integer,0,r_tag,petsc_comm_world, &
        & r_request_n_vertices,i_err)

      call mpi_wait(r_request_n_vertices,r_status,i_err)
      test_request_n_vertices = .false.
      call mpi_test(r_request_n_vertices,test_request_n_vertices, &
        & r_status,i_err)

      do while (test_request_n_vertices .neqv. .true.)
        call mpi_test(r_request_n_vertices,test_request_n_vertices, &
          & r_status,i_err)
      end do

      ! All the processes receive the vertex coordinates
      ! =========================================================================
      ! Allocate memory for vertex coordinates
      allocate(tmp_vx(3*nvertices))
      
      ! Receive vertex coordinates
      r_tag = 5001*nprocs + myprocid
      m_size = 3*nvertices
      call mpi_irecv(tmp_vx,m_size,mpi_double,0,r_tag, petsc_comm_world, &
        & r_request_vertices,i_err)

      call mpi_wait(r_request_vertices,r_status,i_err)
      test_request_vertices = .false.
      call mpi_test(r_request_vertices,test_request_vertices, &
        & r_status,i_err)


      ! Allocate memory for 2D array
      allocate(vx(3,nvertices))
      vx(:,:) = 0.0_wp
      cnt_unpack = 0
      do ii = 1, nvertices
        do jj = 1, 3
          cnt_unpack = cnt_unpack + 1
          vx(jj,ii) = tmp_vx(cnt_unpack)
        end do
      end do

      ! Deallocate memory for tmp_vx
      deallocate(tmp_vx)


      ! All the processes receive e2v (element-to-vertex)
      ! =======================================================================
      ! Allocate memory for 1D temporary array
      n_elems = ihelems(2) - ihelems(1) + 1
      allocate(tmp_e2v(2**ndim*n_elems))
      tmp_e2v = 0

      ! Receive element to vertex connectivity
      r_tag = 100*nprocs + myprocid
      m_size = 2**ndim*n_elems
      call mpi_irecv(tmp_e2v,m_size,mpi_integer,0,r_tag,petsc_comm_world, &
        & r_request_e2v,i_err)

      call mpi_wait(r_request_e2v,r_status,i_err)
      test_request_e2v = .false.
      call mpi_test(r_request_e2v,test_request_e2v, &
        & r_status,i_err)

      ! Unpack 1D array into the e2v 2D array
      allocate(e2v(2**ndim,ihelems(1):ihelems(2)))
      cnt_unpack = 0
      do ii = ihelems(1),ihelems(2)
        do jj = 1, 2**ndim
          cnt_unpack = cnt_unpack + 1
          e2v(jj,ii) = tmp_e2v(cnt_unpack)
        end do
      end do
      
      ! Deallocate memory for tmp_e2v
      deallocate(tmp_e2v)


      ! All the processes receive ef2etmp2
      ! =======================================================================
      allocate(tmp_ef2etmp2(qdim*2*ndim*(ihelems(2)-ihelems(1)+1)))
      tmp_ef2etmp2 = 0

      ! Receive the ordered-element face-to-face connectivity
      r_tag = 200*nprocs + myprocid
      m_size = netmp*2*ndim*qdim
      call mpi_irecv(tmp_ef2etmp2,m_size,mpi_integer,0,r_tag,petsc_comm_world, &
        & r_request_ef2e,i_err)

      call mpi_wait(r_request_ef2e,r_status,i_err)
      test_request_ef2e = .false.
      call mpi_test(r_request_ef2e,test_request_ef2e, &
        & r_status,i_err)

      ! Unpack 1D array into the e2v 2D array
      allocate(ef2etmp2(qdim,2*ndim,ihelems(1):ihelems(2)))
      cnt_unpack = 0
      do ii = ihelems(1),ihelems(2)
        do jj = 1, 2*ndim
          do kk = 1, qdim
            cnt_unpack = cnt_unpack + 1
            ef2etmp2(kk,jj,ii) = tmp_ef2etmp2(cnt_unpack)
          end do
        end do
      end do

      ! Deallocate memory for tmp_ef2etmp2
      deallocate(tmp_ef2etmp2)


      ! All the processes receive ldg_flip_flop_sign_tmp2 
      ! =======================================================================
      ! Allocate memory for 1D temporary array
      n_elems = ihelems(2) - ihelems(1) + 1
      allocate(tmp_ldg_flip_flop_sign_tmp2(2*ndim*n_elems))
      tmp_ldg_flip_flop_sign_tmp2 = 0

      ! Receive element to vertex connectivity
      r_tag = 327*nprocs + myprocid
      m_size = 2*ndim*n_elems
      call mpi_irecv(tmp_ldg_flip_flop_sign_tmp2,m_size,mpi_integer,0,r_tag,petsc_comm_world, &
        & r_request_ldg_flip_flop_sign,i_err)

      call mpi_wait(r_request_ldg_flip_flop_sign,r_status,i_err)
      test_request_ldg_flip_flop_sign = .false.
      call mpi_test(r_request_ldg_flip_flop_sign,test_request_ldg_flip_flop_sign, &
        & r_status,i_err)

      ! Unpack 1D array into the ldg_flip_flop_sign 2D array
      allocate(ldg_flip_flop_sign_tmp2(2*ndim,ihelems(1):ihelems(2)))
      cnt_unpack = 0
      do ii = ihelems(1),ihelems(2)
        do jj = 1, 2*ndim
          cnt_unpack = cnt_unpack + 1
          ldg_flip_flop_sign_tmp2(jj,ii) = tmp_ldg_flip_flop_sign_tmp2(cnt_unpack)
        end do
      end do

!      write(889,*) tmp_ldg_flip_flop_sign_tmp2
      
      ! Deallocate memory for tmp_e2v
      deallocate(tmp_ldg_flip_flop_sign_tmp2)

    end if ! End if myprocid > 0

    ! Allocate memory for ef2e connectivity
    allocate(ef2e(qdim,2*ndim,ihelems(1):ihelems(2)))
    
    ! Assign ef2e to each processor
    ef2e = ef2etmp2

    ! Allocate memory for ldg_flip_flop_sign
    allocate(ldg_flip_flop_sign(2*ndim,ihelems(1):ihelems(2)))
    
    ! Assign ldg_flip_flop_sign to each processor
    ldg_flip_flop_sign = ldg_flip_flop_sign_tmp2

    
    ! Deallocate memory
    deallocate(ef2etmp2)
    deallocate(i2jelems)
    deallocate(icnt)
    deallocate(melemsonproc)
    deallocate(i_s_r_status)
    deallocate(ldg_flip_flop_sign_tmp2)

    ! =========================================================================
    ! Construct the local e2v array (local for each process)
    ! =========================================================================
    if(allocated(vert_list_proc)) deallocate(vert_list_proc)
    allocate(vert_list_proc(nvertices))
    vert_list_proc = nvertices + 100
    
    ! Number of vertices per process
    n_vert_per_proc = 0

    ! Vertex position
    j_vert_pos = 1

    ! Loop over all the elements 
    do i_elem = ihelems(1), ihelems(2)
     
     ! Loop over all the vertices of the element
      v_loop: do j_vert_elem = 1, nverticesperelem
        
        i = e2v(j_vert_elem,i_elem)
        
        ! Loop over the list of known vertices
        do k_vert_proc = 1, n_vert_per_proc
          
          if (i == vert_list_proc(k_vert_proc)) then
          
            ! Point to local e2v entry
            e2v(j_vert_elem,i_elem) = k_vert_proc
            cycle v_loop
          
          end if
        
        end do
        
        ! Add vertex to the list
        vert_list_proc(j_vert_pos) = i
        
        ! Assign new local e2v entry
        e2v(j_vert_elem,i_elem) = j_vert_pos
        j_vert_pos = j_vert_pos + 1
        n_vert_per_proc = n_vert_per_proc + 1
      
      end do v_loop
    end do

    ! Deallocate memory
    deallocate(vert_list_proc)

    ! Wait for other processes
    call mpi_barrier(petsc_comm_world,i_err)

    return
  end subroutine distribute_elements_aflr3

  !============================================================================

  subroutine distributeelements()
    ! Load modules
    use referencevariables
    use variables, only: e2v, elempart, jelems, iae2v, jae2v, ef2e, vx_master, &
                       & vx, periodic_elem_face_ids_x1, periodic_elem_face_ids_x2, &
                       & periodic_elem_face_ids_x3

    ! Nothing is implicitly defined
    implicit none

    integer :: i, i_proc, m
    integer :: i_elem, j, ii, jj
    integer :: netmp, npetmp

    integer, allocatable :: melemsonproc(:), ef2etmp1(:,:,:), ef2etmp2(:,:,:)
    integer, allocatable :: i2jelems(:), j2ielems(:)
    integer, allocatable :: icnt(:)
    integer, allocatable :: e2vtmp(:,:)
    integer :: s_tag, r_tag
    integer :: irsreq(4)
    integer, allocatable :: i_s_r_status(:,:)
    integer :: m_size
    integer :: i_err

    integer :: max_n_vert_per_proc, n_vert_per_proc, j_vert_pos, j_vert_elem, &
             & k_vert_proc

    integer, allocatable :: vert_list_proc(:)

    continue

    allocate(periodic_elem_face_ids_x1(3,0))
    allocate(periodic_elem_face_ids_x2(3,0))
    allocate(periodic_elem_face_ids_x3(3,0))


    ! Number of vertices per element 
    nverticesperelem = 2**ndim

    ! Allocate memory for MPI status array
    allocate(i_s_r_status(mpi_status_size,4))

    ! Allocate memory for lower and upper index of the element ID own by each 
    ! process
    allocate(melemsonproc(0:2*nprocs-1))
    melemsonproc = 0

    ! Allocate memory for counter array
    allocate(icnt(0:nprocs-1))

    ! Process 0 (i.e. the master process) orders the elements for scattering
    if (myprocid == 0) then
      ! Count the number of elements on each process and put it in a container 
      ! as the upper index
      do i_elem = 1, nelems
        melemsonproc(2*elempart(i_elem) + 1) = melemsonproc(2*elempart(i_elem) + 1) + 1
      end do

      ! Make the upper index additive, such that all odd indices contain the
      ! upper element bound of a given process
      do i_proc = 3, 2*nprocs-1, 2
        melemsonproc(i_proc) = melemsonproc(i_proc-2) + melemsonproc(i_proc)
      end do

      ! The first process starts with element ID 1
      melemsonproc(0) = 1

      ! Determine the first element ID on the other processes
      do i_proc = 2, 2*(nprocs-1), 2
        melemsonproc(i_proc) = melemsonproc(i_proc-1) + 1
      end do

      ! Determine the mapping from the old element indices to the new element 
      ! indices
      allocate(i2jelems(1:nelems))
      i2jelems = 0
      
      allocate(j2ielems(1:nelems))
      j2ielems = 0
      
      icnt = 0
      do i_elem = 1, nelems
        i_proc = elempart(i_elem)
        i2jelems(melemsonproc(2*i_proc)+icnt(i_proc)) = i_elem
        j2ielems(i_elem) = melemsonproc(2*i_proc) + icnt(i_proc)
        icnt(i_proc) = icnt(i_proc) + 1
      end do
    else
      allocate(i2jelems(0))
    end if

    ! Create a barrier synchronization in the group. Each task, when reaching 
    ! the MPI_Barrier call, blocks until all tasks in the group reach the same 
    ! MPI_Barrier call. 
    call mpi_barrier(petsc_comm_world,i_err)

    ! The lower and upper ordered-element bounds for each process is 
    ! communicated
    call mpi_scatter(melemsonproc,2,mpi_integer,ihelems,2,mpi_integer,0, &
      & petsc_comm_world,i_err)

    ! Set number of elements 
    nelems = ihelems(2) - ihelems(1) + 1

    ! The original-to-ordered-element mapping is communicated to each process
    ! Each process will receive the global ID of the elements (from the original 
    ! ugrid). Such IDs are stored in jelems.
    allocate(jelems(ihelems(1):ihelems(2)))
    netmp = ihelems(2) - ihelems(1) + 1
    call mpi_scatterv(i2jelems,icnt,melemsonproc(0:2*(nprocs-1):2)-1, &
      & mpi_integer,jelems,netmp,mpi_integer,0,petsc_comm_world,i_err)


    ! Process 0 (i.e. the master process) use blocking sends for sending to each
    ! process the following information:
    !
    ! 1) e2v, connectivity using the original global node numbering (from the 
    ! serial read)
    !
    ! 2) ef2e, connectivity in the new global numbering after metis
    ! has partitioned the grid
    !
    ! 3) coordinates of the vertex (or nodes)
    if (myprocid == 0) then
      
      do i_proc = 0, nprocs-1
        ! Allocate memory for temporary array for element-to-vertex connectivity
        ! on the process
        allocate(e2vtmp(2**ndim,melemsonproc(2*i_proc):melemsonproc(2*i_proc+1)))
        
        ! Allocate memory for temporary array for face-to-face connectivity on 
        ! the process
        allocate(ef2etmp1(3,2*ndim,melemsonproc(2*i_proc):melemsonproc(2*i_proc+1))) 
        
        ! Number of elements on the process
        npetmp = melemsonproc(2*i_proc+1)-melemsonproc(2*i_proc)+1
        
        ! Number of vertices per process. This is a conservative estimate 
        ! because most of the nodes are shared among the element excpet the 
        ! boundary nodes.
        max_n_vert_per_proc = nverticesperelem*(melemsonproc(2*i_proc+1) - &
          & melemsonproc(2*i_proc) + 1)

        ! Construct the local e2v array (local for each process) 
        allocate(vert_list_proc(max_n_vert_per_proc))
        vert_list_proc = max_n_vert_per_proc + 100
    
        ! Initialize number of vertices per process
        n_vert_per_proc = 0

        ! Vertex position in the array
        j_vert_pos = 1
        
        ! Loop over all elements on process
        do i = melemsonproc(2*i_proc), melemsonproc(2*i_proc+1)
          ! Original element index
          ii = i2jelems(i)
          
          ! Fill the e2v connectivity
          do j = 1, 2**ndim
            jj = iae2v(ii) + j - 1
            e2vtmp(j,i) = jae2v(jj)
          end do
          
          ! Loop over all faces
          do j = 1, 2*ndim
            ! Face of neighbor
            ef2etmp1(1,j,i) = ef2e(1,j,ii)
            
            ! Original element index of neighbor
            jj = ef2e(2,j,ii)
            
            ! Internal faces have ef2e(2,j,ii)>0
            if (jj > 0) then
              ! Neighbor element (with new element index)
              ef2etmp1(2,j,i) = j2ielems(jj)
              
              ! Process of neighbor
              ef2etmp1(3,j,i) = elempart(jj)
            
            ! Boundary face
            else
              ! Self is neighbor (new element index)
              ef2etmp1(2,j,i) = i
              
              ! Process of neighbor
              ef2etmp1(3,j,i) = i_proc
            end if
          end do

          ! Loop over all the vertices of the element
          v_loop_m: do j_vert_elem = 1, nverticesperelem
            
            m = e2vtmp(j_vert_elem,i)
            
            ! Loop over the list of known vertices
            do k_vert_proc = 1, n_vert_per_proc
              
              if (m == vert_list_proc(k_vert_proc)) then
              
                cycle v_loop_m
              
              end if
            
            end do
            
            ! Add to vertex to the list
            vert_list_proc(j_vert_pos) = m
            
            ! Update position and counter
            j_vert_pos = j_vert_pos + 1
            n_vert_per_proc = n_vert_per_proc + 1
          
          end do v_loop_m

        end do ! En do loop over the element

        ! Transfer number of vertices
        s_tag = 25*nprocs + i_proc
        m_size = 1
        call mpi_send(n_vert_per_proc,m_size,mpi_integer,i_proc,s_tag, &
          & petsc_comm_world,i_err)

        ! Transfer vertex coordinates
        s_tag = 50*nprocs + i_proc
        m_size = 3*n_vert_per_proc
        call mpi_send(vx_master(1:3,vert_list_proc(1:n_vert_per_proc)),m_size, &
          & mpi_double,i_proc,s_tag,petsc_comm_world,i_err)

        ! Transfer e2v data
        s_tag = 100*nprocs + i_proc
        m_size = npetmp*(2**ndim)
        call mpi_send(e2vtmp(:,:),m_size,mpi_integer,i_proc,s_tag, &
          & petsc_comm_world,i_err)

        ! Transfer ef2e data
        s_tag = 200*nprocs + i_proc
        m_size = npetmp*2*ndim*3
        call mpi_send(ef2etmp1,m_size,mpi_integer,i_proc,s_tag, &
          & petsc_comm_world,i_err)

        ! Deallocate memory
        deallocate(e2vtmp)
        deallocate(ef2etmp1)
        deallocate(vert_list_proc)
      
      end do ! End do processors
      
      ! Deallocate memory
      deallocate(ef2e)
      deallocate(j2ielems)
      deallocate(vx_master)
    
    end if ! End if master process

    ! Receive number of vertices. This number is used to allocate the memory for
    ! the vertex coordinates array.
    r_tag = 25*nprocs + myprocid
    m_size = 1
    call mpi_irecv(nvertices,m_size,mpi_integer,0,r_tag,petsc_comm_world, &
      & irsreq(1),i_err)

    ! Waits for the MPI receive to complete
    call mpi_wait(irsreq(1),i_s_r_status(:,1),i_err)

    ! Allocate memory for vertex coordinates
    allocate(vx(3,nvertices))
    vx(:,:) = 0.0_wp

    ! Receive vertex coordinate
    r_tag = 50*nprocs + myprocid
    m_size = 3*nvertices
    call mpi_irecv(vx,m_size,mpi_double,0,r_tag, petsc_comm_world,irsreq(2), &
      & i_err)

    ! Allocate memory for ordered-element to vertex connectivity 
    allocate(e2v(2**ndim,ihelems(1):ihelems(2)))

    ! Receive element to vertex connectivity
    r_tag = 100*nprocs + myprocid
    m_size = netmp*(2**ndim)
    call mpi_irecv(e2v,m_size,mpi_integer,0,r_tag,petsc_comm_world,irsreq(3), &
      & i_err)

    ! Allocate memory for ordered-element face-to-face connectivity
    allocate(ef2etmp2(3,2*ndim,ihelems(1):ihelems(2)))
    ef2etmp2 = 0

    ! Receive the ordered-element face-to-face connectivity
    r_tag = 200*nprocs + myprocid
    m_size = netmp*2*ndim*3
    call mpi_irecv(ef2etmp2,m_size,mpi_integer,0,r_tag,petsc_comm_world, &
      & irsreq(4),i_err)

    ! Wait for receives to finish
    call MPI_Wait(irsreq(2),i_s_r_status(:,2),i_err)
    call MPI_Wait(irsreq(3),i_s_r_status(:,3),i_err)
    call MPI_Wait(irsreq(4),i_s_r_status(:,4),i_err)

    ! Wait for other processes to arrive here
    call mpi_barrier(petsc_comm_world,i_err)

    ! Allocate memory for ef2e connectivity
    allocate(ef2e(3,2*ndim,ihelems(1):ihelems(2)))
    
    ! Assign global ef2e
    ef2e = ef2etmp2
    
    ! Deallocate memory
    deallocate(ef2etmp2)
    deallocate(i2jelems)
    deallocate(icnt)
    deallocate(melemsonproc)
    deallocate(i_s_r_status)

    ! Construct the local e2v array (local for each process)
    if(allocated(vert_list_proc)) deallocate(vert_list_proc)
    allocate(vert_list_proc(nvertices))
    vert_list_proc = nvertices + 100
    
    ! Number of vertices per process
    n_vert_per_proc = 0

    ! Vertex position
    j_vert_pos = 1

    ! Loop over all the elements 
    do i_elem = ihelems(1), ihelems(2)
     
      ! Loop over all the vertices of the element
      v_loop: do j_vert_elem = 1, nverticesperelem
        
        i = e2v(j_vert_elem,i_elem)
        
        ! Loop over the list of known vertices
        do k_vert_proc = 1,n_vert_per_proc
          
          if (i == vert_list_proc(k_vert_proc)) then
          
            ! Point to local e2v entry
            e2v(j_vert_elem,i_elem) = k_vert_proc
            cycle v_loop
          
          end if
        
        end do
        
        ! Add vertex to the list
        vert_list_proc(j_vert_pos) = i
        
        ! Assign new local e2v entry
        e2v(j_vert_elem,i_elem) = j_vert_pos
        j_vert_pos = j_vert_pos + 1
        n_vert_per_proc = n_vert_per_proc + 1
      
      end do v_loop
    end do

    ! Deallocate memory
    deallocate(vert_list_proc)

    ! Wait for other processes
    call mpi_barrier(petsc_comm_world,i_err)

    return
  end subroutine distributeelements

  !============================================================================
  
  subroutine PetscGridLocations_LGL()
    ! Initialize the global and ghost arrays for the grid
    ! It is run in parallel by all processes
    use referencevariables
    use variables, only: xg, xghst_LGL, ef2e, kfacenodes
    use petscvariables, only: xpetsc
    implicit none

    PetscErrorCode ierrpetsc
    integer :: ntot
    real(wp), allocatable :: yy(:)
    integer, allocatable :: iy(:)
    integer :: ielem, inode, iloc, idir, iface
    integer :: i, ii, jj
    Vec xlocal
    real(wp), pointer, dimension(:) :: xx_v
    PetscScalar xinit

    ! Set number of ghost points to zero. These are LGL points.
    nghost = 0

    ! count necessary ghost points for each process
    ! loop over elements
  
    do ielem = ihelems(1), ihelems(2)

      ! loop over faces
      do iface = 1, nfacesperelem
        ! if face neighbor is off process, then add ghost nodes
        if (ef2e(3,iface,ielem) /= myprocid) nghost = nghost + nodesperface
      end do
    end do
    ! nghost_elem is the total storage needed for the entire adjoining element
    ! it is needed for the jacobian of the viscous penaltys, but for no other operations 
    nghost_elem = nghost * nodesperedge

    ! allocate ghost indices for physical grid locations
    allocate(iy(nghost*3))
    iloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! if neighbor does not need off process data then cycle
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element index of neighbor
        ii = ef2e(2,iface,ielem)
        ! connected face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on face of neighbor
        do i = 1, nodesperface
          ! loop over physical directions
          do idir = 1,3
            ! advance bookkeeping index
            iloc = iloc+1
            ! global index of physical coordinate in 1D parallel array
            iy(iloc) = 3*nodesperelem*(ii-1) + (kfacenodes(i,jj)-1)*3 + idir
          end do
        end do
      end do
    end do

    ! actual ghost cells
    allocate(xghst_LGL(3,nghost))
    xghst_LGL = -100._wp

    ! total number of local values in physical coordinate array
    ntot = 3*nelems*nodesperelem

    ! create global vector with ghosts for grid locations (this is a PETSc feature)
    !
    ! use C array ordering
    iy = iy-1
    ! call to create parallel vector xpetsc with ghost points
    call VecCreateGhost(petsc_comm_world, ntot, petsc_decide, nghost*3, iy, xpetsc, ierrpetsc)
    
    ! initialize to zero
    xinit = 0.0_wp
    call VecSet(xpetsc, xinit, ierrpetsc)
    deallocate(iy)

    ! assemble vector
    call VecAssemblyBegin(xpetsc,ierrpetsc)
    call VecAssemblyEnd  (xpetsc,ierrpetsc)

    ! temporary arrays for filling vector with data
    allocate(iy(3*nodesperelem))
    allocate(yy(3*nodesperelem))

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      yy = 0.0_wp
      ! loop over nodes
      do inode = 1, nodesperelem
        ! loop over physical direction
        do idir = 1,3
          ! update temporary value for physical coordinate
          yy(3*(inode-1)+idir) = xg(idir,inode,ielem)
          ! location in global vector for physical coordinate
          iy(3*(inode-1)+idir) = 3*nodesperelem*(ielem-1)+3*(inode-1)+idir
        end do
      end do
      ! use C indexing
      iy = iy-1
      ! set values in petsc vector
      call VecSetValues(xpetsc,3*nodesperelem,iy,yy,insert_values,ierrpetsc)
    end do

    ! parallel assembly
    call VecAssemblyBegin(xpetsc,ierrpetsc)
    call VecAssemblyEnd  (xpetsc,ierrpetsc)

    ! update local ghost values
    call VecGhostUpdateBegin(xpetsc, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd  (xpetsc, insert_values, scatter_forward, ierrpetsc)

    ! create a sequential PETSc vector for the local physical coordinates
    call VecCreateSeq(petsc_comm_self, ntot+nghost*3, xlocal, ierrpetsc)
    ! fill local vector
    call VecGhostGetLocalForm(xpetsc, xlocal, ierrpetsc)
    ! use fortran pointer to access data conveniently
    call VecGetArrayF90(xlocal, xx_v, ierrpetsc)
    iloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! cycle if a neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over nodes on neighboring face
        do i = 1, nodesperface
          ! update location in ghost array
          iloc = iloc+1
          ! loop over direction
          do idir = 1,3
            ! fill ghost data
            xghst_LGL(idir,iloc) = xx_v(ntot+3*(iloc-1)+idir)
          end do
        end do
      end do
    end do

    ! release pointer
    call VecRestoreArrayF90(xlocal,xx_v,ierrpetsc) 
    call VecGhostRestoreLocalForm(xpetsc,xlocal,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)
    call VecDestroy(xlocal, ierrpetsc)

    deallocate(iy,yy)

  end subroutine PetscGridLocations_LGL

  !============================================================================
  
  subroutine PetscGridLocations_Gau()
    ! Initialize the global and ghost arrays for the grid
    ! It is run in parallel by all processes
    use referencevariables
    use variables, only: xg_Gau_shell, xgghst_Gau_shell
    use petscvariables, only: xpetsc_shell, xlocpetsc_shell
    use collocationvariables, only: n_Gau_2d_p1
    implicit none

    allocate(xgghst_Gau_shell(ndim,nghost_Gau_shell))

!   call PetscComm1D_shell_DataSetup_Gau(xg_Gau_shell,xgghst_Gau_shell,xpetsc_shell,xlocpetsc_shell, &
!                           & ndim, n_Gau_2d_p1, nelems, nghost_Gau_shell)

!   call UpdateComm1D_shell_GhostData_Gau(xg_Gau_shell,xgghst_Gau_shell,xpetsc_shell,xlocpetsc_shell, &
!                           & ndim, n_Gau_2d_p1, ihelems, nghost_Gau_shell)


  end subroutine PetscGridLocations_Gau
  

  !============================================================================

  subroutine Petsc_shell_Counter()
    ! Initialize the global and ghost arrays for the grid
    ! It is run in parallel by all processes
    use referencevariables
    use variables, only: ef2e
    use collocationvariables, only: n_Gau_2d_p1, elem_props
    implicit none

    integer :: ielem, iface

    ! Set number of ghost points to zero. These are LGL points.
    nghost_Gau_shell = 0

    ! count necessary ghost points for each process
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1, nfacesperelem
        if(ef2e(4,iface,ielem) /= elem_props(2,ielem)) then
           write(*,*)'found dissimilar element orders'
           ! if face neighbor is off process, then add ghost nodes
           if (ef2e(3,iface,ielem) /= myprocid) nghost_Gau_shell = nghost_Gau_shell + n_Gau_2d_p1
        endif
      end do
    end do

  end subroutine Petsc_shell_Counter

  !============================================================================

  subroutine PetscComm1DDataSetup(vin, vghstin, vpetscin, vlocin, nq, nk, ne, ngh)
    ! this routine allocates the ghost data for Navier Stokes computations
    use referencevariables
    use variables, only: ef2e, kfacenodes
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nk, ne, ngh
    real(wp), intent(in) :: vin(nq,nk,ne)
    real(wp), intent(in) :: vghstin(nq,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntot
    integer, allocatable :: iyu(:)
    !real(wp), allocatable :: yu(:)
    integer :: ielem, iloc, jloc, iface
    integer :: i, ii, jj, ieq
    PetscScalar xinit

    xinit = 0.0_wp

    ! allocate memory for ghost locations
    allocate(iyu(ngh*nq))
    iloc = 0
    jloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor face
        do i = 1, nodesperface
          ! loop over equations
          do ieq = 1, nq
            ! advance position in ghost array
            jloc = jloc+1
            ! set position of ghost in global vector containing solution data
            iyu(jloc) = nq*nodesperelem*(ii-1) + (kfacenodes(i,jj)-1)*nq + ieq
          end do
        end do
      end do
    end do

    ! total length of on process data for 1D vector of solution
    ntot = nq*nelems*nodesperelem

    ! create global vector for solution with ghosts
    ! 
    ! use C indexing
    iyu = iyu-1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntot, petsc_decide, &
      ngh*nq, iyu, vpetscin, ierrpetsc)
    ! initialize to zero
    call VecSet(vpetscin, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd  (vpetscin,ierrpetsc)

    deallocate(iyu)

    ! create container for local vector that contains on process plus ghost data for solution
    call VecCreateSeq(petsc_comm_self, ntot+ngh*nq, vlocin, ierrpetsc)

  end subroutine PetscComm1DDataSetup

!============================================================================

  subroutine PetscComm1DElementDataSetup(vin, vghstin, vpetscin, vlocin, nq, n_pts_3d, ne, ngh)
    ! this routine allocates the ghost data for Navier Stokes computations
    ! ngh = number of ghost points in adjoining elements.  ngh = nghost*nodesperedge
    ! n_pts_3d is either n_LGL_pts_hex or n_Gau_pts_hex
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, n_pts_3d, ne, ngh
    real(wp), intent(in) :: vin(nq,n_pts_3d,ne)
    real(wp), intent(in) :: vghstin(nq,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntot
    integer, allocatable :: iyu(:)
    integer :: ielem, jloc, iface
    integer :: i, ii, jj, ieq
    PetscScalar xinit

    xinit = 0.0_wp

    ! allocate memory for ghost locations
    allocate(iyu(ngh*nq))
    jloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor elements
        do i = 1, n_pts_3d
          ! loop over equations
          do ieq = 1, nq
            ! advance position in ghost array
            jloc = jloc+1
            ! set position of ghost in global vector containing solution data
            iyu(jloc) = nq * n_pts_3d * (ii-1) + (i-1)*nq + ieq
          end do
        end do
      end do
    end do

    ! total length of on process data for 1D vector of solution
    ntot = nq*nelems * n_pts_3d

    ! create global vector for solution with ghosts
    ! 
    ! use C indexing
    iyu = iyu-1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntot, petsc_decide, &
      ngh*nq, iyu, vpetscin, ierrpetsc)
    ! initialize to zero
    call VecSet(vpetscin, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd(vpetscin,ierrpetsc)

    deallocate(iyu)

    ! create container for local vector that contains on process plus ghost data for solution
    call VecCreateSeq(petsc_comm_self, ntot+ngh*nq, vlocin, ierrpetsc)

  end subroutine PetscComm1DElementDataSetup

  subroutine PetscComm2DDataSetup(vin, vghstin, vpetscin, vlocin, nq, nd, nk, ne, ngh)
    ! this routine allocates the ghost data for Navier Stokes computations
    use referencevariables
    use variables, only: ef2e, kfacenodes
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nd, nk, ne, ngh
    real(wp), intent(in) :: vin(nq,nd,nk,ne)
    real(wp), intent(in) :: vghstin(nq,nd,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntotphi
    !real(wp), allocatable :: yphi(:)
    integer, allocatable :: iyphi(:)
    integer :: ielem, iloc, idir, iface
    integer :: i, ii, jj, ieq

    PetscScalar xinit

    xinit = 0.0_wp

    ! allocate memory for ghost locations
    allocate(iyphi(ngh*nq*nd))
    iloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor face
        do i = 1, nodesperface
          ! loop over gradient directions
          do idir = 1,nd
            ! loop over equations
            do ieq = 1,nq
              ! advance position in ghost array
              iloc = iloc+1
              ! set position of ghost in global vector containing
              ! entropy variable gradient data
              iyphi(iloc) = nd*nq*nodesperelem*(ii-1) &
                + (kfacenodes(i,jj)-1)*nd*nq + nq*(idir-1) + ieq
            end do
          end do
        end do
      end do
    end do

    ! total length of on process data for 1D vector of entropy variable gradients
    ntotphi = nd*nq*nelems*nodesperelem

    ! create global vector for entropy variable gradients with ghosts
    ! 
    ! use C indexing
    iyphi = iyphi-1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntotphi, petsc_decide, ngh*nq*nd, iyphi, vpetscin, ierrpetsc)
    ! initialize to zero
    call VecSet(vpetscin, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd(vpetscin,ierrpetsc)

    deallocate(iyphi)

    ! create container for local vector that contains on process plus ghost data for gradients
    call VecCreateSeq(petsc_comm_self, ntotphi+ngh*nq*nd, vlocin, ierrpetsc)

    return
  end subroutine PetscComm2DDataSetup

  !============================================================================

  subroutine PetscComm2DGeomDataSetup(vin,vghstin,vpetscin,vlocin,nd,nk,ne,ngh)
    ! this routine allocates the ghost data for Navier Stokes computations
    use referencevariables
    use variables, only: ef2e, kfacenodes
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nd, nk, ne, ngh
    real(wp), intent(in) :: vin(nd,nd,nk,ne)
    real(wp), intent(in) :: vghstin(nd,nd,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntot_r_x
    integer, allocatable :: iy_r_x(:)
    integer :: ielem, iloc, jloc, idir, iface
    integer :: i, ii, jj, icomp

    PetscScalar xinit

    xinit = 0.0_wp

    ! allocate memory for ghost locations
    allocate(iy_r_x(ngh*nd*nd))
    iloc = 0
    jloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor face
        do i = 1, nodesperface
          ! loop over gradient directions
          do idir = 1,nd
            ! loop over equations
            do icomp = 1,nd
              ! advance position in ghost array
              iloc = iloc + 1
              ! set position of ghost in global vector containing
              ! the contravariant vectors
              iy_r_x(iloc) = nd*nd*nodesperelem*(ii-1) &
                + (kfacenodes(i,jj)-1)*nd*nd + nd*(idir-1) + icomp
            end do
          end do
        end do
      end do
    end do

    ! total length of on process data for 1D vector of contravariant vectors
    ntot_r_x = nd*nd*nelems*nodesperelem

    ! create global vector for contravariant vectors with ghosts
    ! 
    ! use C indexing
    iy_r_x = iy_r_x - 1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntot_r_x, petsc_decide, ngh*nd*nd,iy_r_x,vpetscin,ierrpetsc)
    ! initialize to zero
    call VecSet(vpetscin,xinit,ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd(vpetscin,ierrpetsc)

    deallocate(iy_r_x)

    ! create container for local vector that contains on process plus ghost data 
    ! for contravariant vectors
    call VecCreateSeq(petsc_comm_self,ntot_r_x+ngh*nd*nd,vlocin,ierrpetsc)

    return
  end subroutine PetscComm2DGeomDataSetup

  !============================================================================
  !============================================================================
  !  End of Communication setup routines
  !============================================================================
  !============================================================================

  subroutine UpdateComm1DGhostData(vin, vghstin, vpetscin, vlocin, nq, nk, ih, ngh)
    ! this subroutine communicates solution data across processes and updates
    ! the array ughst.
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer,  intent(in) :: nq, nk, ngh
    integer,  intent(in) :: ih(2)
    real(wp), intent(in) :: vin(nq,nk,ih(1):ih(2))
    real(wp), intent(inout) :: vghstin(nq,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntotu
    real(wp), allocatable ::  yu(:)
    integer,  allocatable :: iyu(:)
    integer :: ielem, iloc, inode, iface
    integer :: i, ieq

    real(wp), pointer :: xx_v(:)

    ! length of arrays for filling global vectors with data
    ntotu = nequations*nodesperelem
    allocate(iyu(ntotu))
    allocate( yu(ntotu))

    ! loop over elements
    do ielem = ih(1), ih(2)
      yu = 0.0_wp
      ! loop over nodes
      do inode = 1, nodesperelem
        ! loop over variables
        do ieq = 1, nq
          ! update temporary solution values
          yu(nq*(inode-1)+ieq) = vin(ieq,inode,ielem)
          ! update global location of solution values
          iyu(nq*(inode-1)+ieq) = ntotu*(ielem-1)+nq*(inode-1)+ieq-1
        end do
      end do
      ! set values in petsc vector
      call VecSetValues(vpetscin,ntotu,iyu,yu,insert_values,ierrpetsc)
    end do

    ! assemble vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd  (vpetscin,ierrpetsc)
    ! update ghost values
    call VecGhostUpdateBegin(vpetscin, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd  (vpetscin, insert_values, scatter_forward, ierrpetsc)

    ! total length of on process data
    ntotu = nelems*nodesperelem*nq
    ! get local data including ghost points
    call VecGhostGetLocalForm(vpetscin, vlocin, ierrpetsc)
    ! use fortran pointer for convenience
    call VecGetArrayF90(vlocin, xx_v, ierrpetsc)
    iloc = 0
    ! loop over elements
    do ielem = ih(1), ih(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over nodes
        do i = 1, nodesperface
          ! update ghost node index
          iloc = iloc+1
          ! loop over equations
          do ieq = 1,nq
            ! fill ghost data
            vghstin(ieq,iloc) = xx_v(ntotu+nq*(iloc-1)+ieq)
          end do
        end do
      end do
    end do

    ! release pointer
    call VecRestoreArrayF90(vlocin,xx_v,ierrpetsc)
    call VecGhostRestoreLocalForm(vpetscin,vlocin,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)

  end subroutine UpdateComm1DGhostData

  !============================================================================
  
  !============================================================================
  
  subroutine UpdateComm1DElementGhostData(uin, ughstin, upetscin, ulocin, nq, nk, ih, ngh)
    ! this subroutine communicates solution data across processes and updates
    ! the array ughst.
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nk, ngh
    integer, intent(in) :: ih(2)
    real(wp), intent(in) :: uin(nq,nk,ih(1):ih(2))
    real(wp), intent(inout) :: ughstin(nq,ngh)
    Vec upetscin
    Vec ulocin

    PetscErrorCode ierrpetsc
    integer :: ntotu
    real(wp), allocatable :: yu(:)
    integer, allocatable :: iyu(:)
    integer :: ielem, inode, iloc, iface
    integer :: i, ieq

    real(wp), pointer :: xx_v(:)

    ! length of arrays for filling global vectors with data
    ntotu = nequations*nodesperelem
    allocate(iyu(ntotu))
    allocate( yu(ntotu))

    ! loop over elements
    do ielem = ih(1), ih(2)
      yu = 0.0_wp
      ! loop over nodes
      do inode = 1, nodesperelem
        ! loop over variables
        do ieq = 1, nq
          ! update temporary solution values
           yu(nq*(inode-1)+ieq) = uin(ieq,inode,ielem)
          ! update global location of solution values
          iyu(nq*(inode-1)+ieq) = ntotu*(ielem-1)+nq*(inode-1)+ieq-1
        end do
      end do
      ! set values in petsc vector
      call VecSetValues(upetscin,ntotu,iyu,yu,insert_values,ierrpetsc)
    end do

    ! assemble vector
    call VecAssemblyBegin(upetscin,ierrpetsc)
    call VecAssemblyEnd  (upetscin,ierrpetsc)
    ! update ghost values
    call VecGhostUpdateBegin(upetscin, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd  (upetscin, insert_values, scatter_forward, ierrpetsc)

    ! total length of on process data
    ntotu = nelems*nodesperelem*nq
    ! get local data including ghost points
    call VecGhostGetLocalForm(upetscin, ulocin, ierrpetsc)
    ! use fortran pointer for convenience
    call VecGetArrayF90(ulocin, xx_v, ierrpetsc)
    iloc = 0
    ! loop over elements
    do ielem = ih(1), ih(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over nodes
        do i = 1, nodesperelem
          ! update ghost node index
          iloc = iloc+1
          ! loop over equations
          do ieq = 1,nq
            ! fill ghost data
            ughstin(ieq,iloc) = xx_v(ntotu+nq*(iloc-1)+ieq)
          end do
        end do
      end do
    end do

    ! release pointer
    call VecRestoreArrayF90(ulocin,xx_v,ierrpetsc)
    call VecGhostRestoreLocalForm(upetscin,ulocin,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)

  end subroutine UpdateComm1DElementGhostData

  !============================================================================
  
  subroutine UpdateComm2DGhostData(vin, vghstin, vpetscin, vlocin, nq, nd, nk, ih, ngh)
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nd, nk, ngh
    integer, intent(in) :: ih(2)
    real(wp), intent(in) :: vin(nq,nd,nk,ih(1):ih(2))
    real(wp), intent(inout) :: vghstin(nq,nd,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntotphi
    real(wp), allocatable :: yphi(:)
    integer, allocatable :: iyphi(:)
    integer :: ielem, inode, iloc, idir, iface
    integer :: i, ieq

    real(wp), pointer :: xx_v(:)

    ! length of temporary arrays for filling in global data
    ntotphi = nd*nq*nodesperelem
    allocate(iyphi(ntotphi))
    allocate(yphi(ntotphi))

    ! loop over elements
    do ielem = ih(1), ih(2)
      yphi = 0.0_wp
      ! loop over nodes
      do inode = 1, nodesperelem
        ! loop over gradient direction
        do idir = 1,3
          ! loop over variables
          do ieq = 1, nq
            ! update gradient data
            yphi(nq*nd*(inode-1)+nq*(idir-1)+ieq) = vin(ieq,idir,inode,ielem)
            ! update global location of gradient data in 1D vector
            iyphi(nq*nd*(inode-1)+nq*(idir-1)+ieq) = ntotphi*(ielem-1) + nq*nd*(inode-1) + nq*(idir-1)+ieq-1
          end do
        end do
      end do
      ! set values in petsc vector
      call VecSetValues(vpetscin,ntotphi,iyphi,yphi,insert_values,ierrpetsc)
    end do

    ! assemble parallel vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd(vpetscin,ierrpetsc)
    ! update ghost data
    call VecGhostUpdateBegin(vpetscin, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd(vpetscin, insert_values, scatter_forward, ierrpetsc)

    ! total length of on process data
    ntotphi = nd*nq*nelems*nodesperelem
    ! get local vector with ghost data
    call VecGhostGetLocalForm(vpetscin, vlocin, ierrpetsc)
    ! use fortran 90 pointers for convenience
    call VecGetArrayF90(vlocin, xx_v, ierrpetsc)
    iloc = 0
    ! loop over elements
    do ielem = ih(1), ih(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over indices
        do i = 1, nodesperface
          ! update location of ghost node
          iloc = iloc+1
          ! loop over directions
          do idir = 1, nd
            ! loop over variables
            do ieq = 1,nq
              ! fill ghost values
              vghstin(ieq,idir,iloc) = xx_v(ntotphi+nd*nq*(iloc-1)+nq*(idir-1)+ieq)
            end do
          end do
        end do
      end do
    end do

    ! release pointer
    call VecRestoreArrayF90(vlocin,xx_v,ierrpetsc) 
    call VecGhostRestoreLocalForm(vpetscin,vlocin,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)

  end subroutine UpdateComm2DGhostData

  !============================================================================

  subroutine UpdateComm2DGeomGhostData(vin,vghstin,vpetscin,vlocin,nd,nk,ih,ngh)
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nd, nk, ngh
    integer, intent(in) :: ih(2)
    real(wp), intent(in) :: vin(nd,nd,nk,ih(1):ih(2))
    real(wp), intent(inout) :: vghstin(nd,nd,ngh)
    Vec vpetscin
    Vec vlocin

    PetscErrorCode ierrpetsc
    integer :: ntot_r_x
    real(wp), allocatable :: y_r_x(:)
    integer, allocatable :: iy_r_x(:)
    integer :: ielem, inode, iloc, idir, iface
    integer :: i, icomp
    !Vec xlocal
    !PetscScalar xinit
    real(wp), pointer :: xx_v(:)

    ! length of temporary arrays for filling in global data
    ntot_r_x = nd*nd*nodesperelem
    allocate(iy_r_x(ntot_r_x))
    allocate(y_r_x(ntot_r_x))

    ! loop over elements
    do ielem = ih(1), ih(2)
      y_r_x = 0.0_wp
      ! loop over nodes
      do inode = 1, nodesperelem
        ! loop over the directions
        do idir = 1,3
          ! loop over variables
          do icomp = 1, nd
            ! update geometrical data data
            y_r_x(nd*nd*(inode-1)+nd*(idir-1)+icomp) = vin(icomp,idir,inode,ielem)
            ! update global location of geometrical data in 1D vector
            iy_r_x(nd*nd*(inode-1)+nd*(idir-1)+icomp) = ntot_r_x*(ielem-1) &
              + nd*nd*(inode-1) + nd*(idir-1)+icomp-1
          end do
        end do
      end do
      ! set values in petsc vector
      call VecSetValues(vpetscin,ntot_r_x,iy_r_x,y_r_x,insert_values,ierrpetsc)
    end do

    ! assemble parallel vector
    call VecAssemblyBegin(vpetscin,ierrpetsc)
    call VecAssemblyEnd(vpetscin,ierrpetsc)
    ! update ghost data
    call VecGhostUpdateBegin(vpetscin, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd(vpetscin, insert_values, scatter_forward, ierrpetsc)

    ! total length of on process data
    ntot_r_x = nd*nd*nelems*nodesperelem
    ! get local vector with ghost data
    call VecGhostGetLocalForm(vpetscin, vlocin, ierrpetsc)
    ! use fortran 90 pointers for convenience
    call VecGetArrayF90(vlocin, xx_v, ierrpetsc)
    iloc = 0
    ! loop over elements
    do ielem = ih(1), ih(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over indices
        do i = 1, nodesperface
          ! update location of ghost node
          iloc = iloc+1
          ! loop over directions
          do idir = 1, nd
            ! loop over variables
            do icomp = 1,nd
              ! fill ghost values
              vghstin(icomp,idir,iloc) = xx_v(ntot_r_x+nd*nd*(iloc-1)+nd*(idir-1)+icomp)
            end do
          end do
        end do
      end do
    end do

    ! release pointer
    call VecRestoreArrayF90(vlocin,xx_v,ierrpetsc) 
    call VecGhostRestoreLocalForm(vpetscin,vlocin,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)

  end subroutine UpdateComm2DGeomGhostData

  !============================================================================

  subroutine PetscComm1DDataSetupWENO(vinWENO, vghstinWENO, vpetscinWENO, vlocinWENO, nq, nk, ne, ngh)
    ! this routine allocates the ghost data for Navier Stokes computations
    use referencevariables
    use variables, only: ef2e, kfacenodesWENO
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nk, ne, ngh
    real(wp), intent(in) :: vinWENO(nq,nk,ne)
    real(wp), intent(in) :: vghstinWENO(nq,ngh)
    Vec vpetscinWENO
    Vec vlocinWENO

    PetscErrorCode ierrpetsc
    integer :: ntot
    integer, allocatable :: iyu(:)
    !real(wp), allocatable :: yu(:)
    integer :: ielem, iloc, jloc, iface
    integer :: i, ii, jj, ieq
    PetscScalar xinit

    xinit = 0.0_wp

    ! allocate memory for ghost locations
    allocate(iyu(ngh*nq))
    iloc = 0
    jloc = 0
    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor face
        do i = 1, nodesperface
          ! loop over equations
          do ieq = 1, nq
            ! advance position in ghost array
            jloc = jloc+1
            ! set position of ghost in global vector containing solution data
            iyu(jloc) = nq*nodesperelem*(ii-1) + (kfacenodesWENO(i,jj)-1)*nq + ieq
          end do
        end do
      end do
    end do

    ! total length of on process data for 1D vector of solution
    ntot = nq*nelems*nodesperelem

    ! create global vector for solution with ghosts
    ! 
    ! use C indexing
    iyu = iyu-1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntot, petsc_decide, &
      ngh*nq, iyu, vpetscinWENO, ierrpetsc)
    ! initialize to zero
    call VecSet(vpetscinWENO, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(vpetscinWENO,ierrpetsc)
    call VecAssemblyEnd(vpetscinWENO,ierrpetsc)

    deallocate(iyu)

    ! create container for local vector that contains on process plus ghost data for solution
    call VecCreateSeq(petsc_comm_self, ntot+ngh*nq, vlocinWENO, ierrpetsc)

  end subroutine PetscComm1DDataSetupWENO

  !============================================================================
  
  subroutine UpdateComm1DGhostDataWENO(vinWENO, vghstinWENO, vpetscinWENO, vlocinWENO, nq, nk, ih, ngh)
    ! this subroutine communicates solution data across processes and updates
    ! the array ughst.
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nk, ngh
    integer, intent(in) :: ih(2)
    real(wp), intent(in) :: vinWENO(nq,nk,ih(1):ih(2))
    real(wp), intent(inout) :: vghstinWENO(nq,ngh)
    Vec vpetscinWENO
    Vec vlocinWENO

    PetscErrorCode ierrpetsc
    integer :: ntotu
    real(wp), allocatable :: yu(:)
    integer, allocatable :: iyu(:)
    integer :: ielem, iloc, inode, iface
    integer :: i, ieq
    !Vec xlocal
    !PetscScalar xinit
    real(wp), pointer :: xx_v(:)

    ! length of arrays for filling global vectors with data
    ntotu = nequations*nodesperelem
    allocate(iyu(ntotu))
    allocate( yu(ntotu))

    ! loop over elements
    do ielem = ih(1), ih(2)
      yu = 0.0_wp
      ! loop over nodes
      do inode = 1, nodesperelem
        ! loop over variables
        do ieq = 1, nq
          ! update temporary solution values
           yu(nq*(inode-1)+ieq) = vinWENO(ieq,inode,ielem)
          ! update global location of solution values
          iyu(nq*(inode-1)+ieq) = ntotu*(ielem-1)+nq*(inode-1)+ieq-1
        end do
      end do
      ! set values in petsc vector
       call VecSetValues(vpetscinWENO,ntotu,iyu,yu,insert_values,ierrpetsc)
    end do

     ! assemble vector
     call VecAssemblyBegin(vpetscinWENO,ierrpetsc)
     call VecAssemblyEnd  (vpetscinWENO,ierrpetsc)
     ! update ghost values
     call VecGhostUpdateBegin(vpetscinWENO, insert_values, scatter_forward, ierrpetsc)
     call VecGhostUpdateEnd  (vpetscinWENO, insert_values, scatter_forward, ierrpetsc)
 
     ! total length of on process data
     ntotu = nelems*nodesperelem*nq
     ! get local data including ghost points
     call VecGhostGetLocalForm(vpetscinWENO, vlocinWENO, ierrpetsc)
     ! use fortran pointer for convenience
     call VecGetArrayF90(vlocinWENO, xx_v, ierrpetsc)
     iloc = 0
     ! loop over elements
     do ielem = ih(1), ih(2)
       ! loop over faces
       do iface = 1,nfacesperelem
         ! do nothing if neighbor is on process
         if (ef2e(3,iface,ielem) == myprocid) cycle
         ! loop over nodes
         do i = 1, nodesperface
           ! update ghost node index
           iloc = iloc+1
           ! loop over equations
           do ieq = 1,nq
             ! fill ghost data
             vghstinWENO(ieq,iloc) = xx_v(ntotu+nq*(iloc-1)+ieq)
           end do
         end do
       end do
     end do
 
     ! release pointer
     call VecRestoreArrayF90(vlocinWENO,xx_v,ierrpetsc)
     call VecGhostRestoreLocalForm(vpetscinWENO,vlocinWENO,ierrpetsc)
     if(associated(xx_v)) deallocate(xx_v)

  end subroutine UpdateComm1DGhostDataWENO

  !============================================================================
  !  The next three routines allocate the PETSC buckets for exchanging the
  !  partner node geometrical data
  !============================================================================

  subroutine PetscComm1DDataSetupWENOGeom(xinWENO_self, xghstinWENO_partner, &
                              xpetscWENO_partner, xlocpetscWENO_partner,     &
                              nq, nshell, ne, ngh)
    !
    ! Allocates ghost containers for SSWENO Navier Stokes computations
    !
    ! the geometric location of partner data: (x), is needed from the adjoining elements
    ! the (x) locations in the following cartoon are transferred from the element where they are 
    ! generated, to the adjoining element.  It is currently assumed that the adjoining element contains
    ! the physical location (x).  This may not be the case for severly skewed elements.
    !
    !      ------------------------------
    !      |\ (x)                       |
    !      (x)\     X            x      |
    !      |    \                       |
    !      |      \                     |
    !      |        \(x)                |
    !      |   x  (x) \    x       x    |
    !      |            \               |
    !      |          x   \--------------
    !      |              |
    !      |              |
    !      |              |
    !      |              |
    !      |   x          |
    !      |          x   |
    !      ----------------
    !
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nshell, ne, ngh
    real(wp), intent(in) :: xinWENO_self(nq,nshell,ne)
    real(wp), intent(in) :: xghstinWENO_partner(nq,ngh)
    Vec xpetscWENO_partner
    Vec xlocpetscWENO_partner

    PetscErrorCode ierrpetsc
    integer :: ntot
    integer,  allocatable :: iyu(:)
   !real(wp), allocatable ::  yu(:)
    integer :: ielem, iloc, jloc, iface
    integer :: i, ii, jj, ieq
    PetscScalar xinit

    xinit = 0.0_wp

    ! Data lives on element in a container with a prescribed storage structure
    ! iyu points to location in local array that contains data to be transferred.

    ! allocate memory for ghost pointers
    allocate(iyu(ngh*nq))
    iloc = 0 ; jloc = 0 ;

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor face
        do i = 1, nodesperface
          ! loop over equations
          do ieq = 1, nq
            ! advance position in ghost array
            jloc = jloc+1
            ! set position of ghost in global vector containing solution data
            iyu(jloc) = nq*nodesperface*nfacesperelem*(ii-1) &
                       +nq*nodesperface*(jj-1)               &
                       +nq*(i-1) + ieq    
          end do
        end do
      end do
    end do


    ! total length of on process data for 1D vector of solution
    ntot = nq*nelems*nodesperface*nfacesperelem

    ! create global vector for solution with ghosts
    ! 
    ! use C indexing
    iyu = iyu-1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntot, petsc_decide, &
      ngh*nq, iyu, xpetscWENO_partner, ierrpetsc)
    ! initialize to zero
    call VecSet          (xpetscWENO_partner, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(xpetscWENO_partner,ierrpetsc)
    call VecAssemblyEnd  (xpetscWENO_partner,ierrpetsc)

    deallocate(iyu)

    ! create container for local vector that contains on process plus ghost data for solution
    call VecCreateSeq(petsc_comm_self, ntot+ngh*nq, xlocpetscWENO_partner, ierrpetsc)

  end subroutine PetscComm1DDataSetupWENOGeom

  !============================================================================
  
  subroutine UpdateComm1DGhostDataWENOGeom(xinWENO_self, xghstinWENO_partner,  &
                                xpetscWENO_partner, xlocpetscWENO_partner, &
                                nq, nshell, ih, ngh)
    ! this subroutine communicates solution data across processes and updates
    ! the array ughst.
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer,  intent(in)    :: nq, nshell, ngh
    integer,  intent(in)    :: ih(2)
    real(wp), intent(in)    :: xinWENO_self(nq,nshell,ih(1):ih(2))
    real(wp), intent(inout) :: xghstinWENO_partner(nq,ngh)
    Vec xpetscWENO_partner
    Vec xlocpetscWENO_partner

    PetscErrorCode ierrpetsc
    integer :: ntotu
    real(wp), allocatable ::  yu(:)
    integer,  allocatable :: iyu(:)
    integer :: ielem, iloc, iface
    integer :: i, ieq, ipen, jnode
    real(wp), pointer :: xx_v(:)

    ! length of arrays for filling global vectors with data
    ntotu = nq * nshell
    allocate(iyu(ntotu))
    allocate( yu(ntotu))

    do ielem = ih(1), ih(2)                       ! loop over elements

      iyu = 0 ; yu = 0.0_wp ;

      do iface = 1, nfacesperelem                 ! loop over 6 faces

        do ipen = 1, nodesperface                 ! loop over nodes on face

          jnode = (iface-1)*nodesperface + ipen
          do ieq = 1, nq                          ! loop over variables
            ! update temporary solution values
             yu(nq*(jnode-1)+ieq) = xinWENO_self(ieq,jnode,ielem)
            ! update global location of solution values
            iyu(nq*(jnode-1)+ieq) = ntotu*(ielem-1)+nq*(jnode-1)+ieq-1
          end do

        end do
      end do
      ! set values in petsc vector
       call VecSetValues(xpetscWENO_partner,ntotu,iyu,yu,insert_values,ierrpetsc)
    end do

    ! assemble vector
    call VecAssemblyBegin(xpetscWENO_partner,ierrpetsc)
    call VecAssemblyEnd  (xpetscWENO_partner,ierrpetsc)
    ! update ghost values
    call VecGhostUpdateBegin(xpetscWENO_partner, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd  (xpetscWENO_partner, insert_values, scatter_forward, ierrpetsc)

    ! get local data including ghost points
    call VecGhostGetLocalForm(xpetscWENO_partner, xlocpetscWENO_partner, ierrpetsc)
    ! use fortran pointer for convenience
    call VecGetArrayF90(xlocpetscWENO_partner, xx_v, ierrpetsc)

    ! total length of on process data WITHOUT ghost data
    ntotu = nq * nshell * nelems
    iloc = 0
    ! loop over elements
    do ielem = ih(1), ih(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over nodes
        do i = 1, nodesperface
          ! update ghost node index
          iloc = iloc+1
          ! loop over equations
          do ieq = 1,nq
            ! fill ghost data
            xghstinWENO_partner(ieq,iloc) = xx_v(ntotu+nq*(iloc-1)+ieq)
          end do
        end do
      end do
    end do
 
    ! release pointer
    call VecRestoreArrayF90(xlocpetscWENO_partner,xx_v,ierrpetsc)
    call VecGhostRestoreLocalForm(xpetscWENO_partner,xlocpetscWENO_partner,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)


  end subroutine UpdateComm1DGhostDataWENOGeom

  !============================================================================

  subroutine PetscCommShellDataSetup(shellin, shellghstin, shellpetscin, shelllocin, nq, nk, ne, ngh)
    !
    ! Allocate the shell ghost data for Navier Stokes computations
    ! Data is stored in SHELL coordinates (dimension => nk = nodesperface*nfacesperelem )
    !
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer, intent(in) :: nq, nk, ne, ngh
    real(wp), intent(in) :: shellin(nq,nk,ne)
    real(wp), intent(in) :: shellghstin(nq,ngh)
    Vec shellpetscin
    Vec shelllocin

    PetscErrorCode ierrpetsc
    integer :: ntot
    integer, allocatable :: iyu(:)
    !real(wp), allocatable :: yu(:)
    integer :: ielem, iloc, jloc, iface
    integer :: i, ii, jj, ieq, nshell
    PetscScalar xinit

    xinit  = 0.0_wp
    nshell = nk

    ! allocate memory for ghost locations
    allocate(iyu(ngh*nq))
    iloc = 0 ; jloc = 0

    ! loop over elements
    do ielem = ihelems(1), ihelems(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! element of neighbor
        ii = ef2e(2,iface,ielem)
        ! face of neighbor
        jj = ef2e(1,iface,ielem)
        ! loop over nodes on neighbor face
        do i = 1, nodesperface
          ! loop over equations
          do ieq = 1, nq
            ! advance position in ghost array
            jloc = jloc+1
            ! set position of ghost in global vector containing solution data
            iyu(jloc) = nq*nshell      *(ii-1)  &
                       +nq*nodesperface*(jj-1)  &
                       +nq*(i-1) + ieq
          end do
        end do
      end do
    end do

    ! total length of on process data for 1D vector of solution
    ntot = nq*nelems*nshell

    ! create global vector for solution with ghosts
    ! 
    ! use C indexing
    iyu = iyu-1
    ! call to petsc to create global vector with ghosts
    call VecCreateGhost(petsc_comm_world, ntot, petsc_decide, &
      ngh*nq, iyu, shellpetscin, ierrpetsc)
    ! initialize to zero
    call VecSet(shellpetscin, xinit, ierrpetsc)
    ! assemble parallel vector
    call VecAssemblyBegin(shellpetscin,ierrpetsc)
    call VecAssemblyEnd(shellpetscin,ierrpetsc)

    deallocate(iyu)

    ! create container for local vector that contains on process plus ghost data for solution
    call VecCreateSeq(petsc_comm_self, ntot+ngh*nq, shelllocin, ierrpetsc)

  end subroutine PetscCommShellDataSetup

  !============================================================================
  
  subroutine UpdateCommShellGhostData(shellin, shellghstin,  &
                                shellpetsc, shelllocpetsc, &
                                nq, nshell, ih, ngh)
    ! this subroutine communicates solution data across processes and updates
    ! the array ughst.
    use referencevariables
    use variables, only: ef2e
    implicit none

    ! Arguments
    ! =========
    integer,  intent(in)    :: nq, nshell, ngh
    integer,  intent(in)    :: ih(2)
    real(wp), intent(in)    :: shellin(nq,nshell,ih(1):ih(2))
    real(wp), intent(inout) :: shellghstin(nq,ngh)
    Vec shellpetsc
    Vec shelllocpetsc

    PetscErrorCode ierrpetsc
    integer :: ntotu
    real(wp), allocatable ::  yu(:)
    integer,  allocatable :: iyu(:)
    integer :: ielem, iloc, iface
    integer :: i, ieq, ipen, jnode
    real(wp), pointer :: xx_v(:)

    ! length of arrays for filling global vectors with data
    ntotu = nq * nshell
    allocate(iyu(ntotu))
    allocate( yu(ntotu))

    do ielem = ih(1), ih(2)                       ! loop over elements

      iyu = 0 ; yu = 0.0_wp ;

      do iface = 1, nfacesperelem                 ! loop over 6 faces

        do ipen = 1, nodesperface                 ! loop over nodes on face

          jnode = (iface-1)*nodesperface + ipen
          do ieq = 1, nq                          ! loop over variables
            ! update temporary solution values
             yu(nq*(jnode-1)+ieq) = shellin(ieq,jnode,ielem)
            ! update global location of solution values
            iyu(nq*(jnode-1)+ieq) = ntotu*(ielem-1)+nq*(jnode-1)+ieq-1
          end do

        end do
      end do
      ! set values in petsc vector
       call VecSetValues(shellpetsc,ntotu,iyu,yu,insert_values,ierrpetsc)
    end do

    ! assemble vector
    call VecAssemblyBegin(shellpetsc,ierrpetsc)
    call VecAssemblyEnd  (shellpetsc,ierrpetsc)
    ! update ghost values
    call VecGhostUpdateBegin(shellpetsc, insert_values, scatter_forward, ierrpetsc)
    call VecGhostUpdateEnd  (shellpetsc, insert_values, scatter_forward, ierrpetsc)

    ! get local data including ghost points
    call VecGhostGetLocalForm(shellpetsc, shelllocpetsc, ierrpetsc)
    ! use fortran pointer for convenience
    call VecGetArrayF90(shelllocpetsc, xx_v, ierrpetsc)

    ! total length of on process data WITHOUT ghost data
    ntotu = nq * nshell * nelems
    iloc = 0
    ! loop over elements
    do ielem = ih(1), ih(2)
      ! loop over faces
      do iface = 1,nfacesperelem
        ! do nothing if neighbor is on process
        if (ef2e(3,iface,ielem) == myprocid) cycle
        ! loop over nodes
        do i = 1, nodesperface
          ! update ghost node index
          iloc = iloc+1
          ! loop over equations
          do ieq = 1,nq
            ! fill ghost data
            shellghstin(ieq,iloc) = xx_v(ntotu+nq*(iloc-1)+ieq)
          end do
        end do
      end do
    end do
 
    ! release pointer
    call VecRestoreArrayF90(shelllocpetsc,xx_v,ierrpetsc)
    call VecGhostRestoreLocalForm(shellpetsc,shelllocpetsc,ierrpetsc)
    if(associated(xx_v)) deallocate(xx_v)


  end subroutine UpdateCommShellGhostData
  
end module mpimod
