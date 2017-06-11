! This module is the code API for PETSc nonlinear solver SNES.
! Ths subroutines implemented here perform the following tasks:
!
! - initialization of the PETSc solver
! - assembling of matrices and vectors to solve the nonlinear system of 
!   algebraic equations
! - monitoring PETSc performances
! - finalization of the PETSc solver

module petsc_snes_solver

  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none

  ! Subroutines and functions in this module are usually private
  private

  ! Exceptions, i.e. public subroutines and functions
  public star_petsc_perf_mon
  public finish_petsc_perf_mon
  public setup_petsc_solver
  public solve_implicit_petsc
  public finalize_petsc_solver


! C PETSc header files to be included 
#include "finclude/petscsys.h"
#include "finclude/petscpc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"
#include "finclude/petscksp.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

  ! PETSc variables
  ! ===============
  ! Petsc event
  PetscLogEvent :: petsc_user_event

  ! SNES context
  SNES :: petsc_snes

  ! KSP context
  KSP :: petsc_ksp

  ! Preconditioner
  PC :: petsc_pc
  
  ! Solution and nonlinear residual vectors
  Vec :: petsc_x_vec, petsc_f_vec
  
  !  Residual Jacobian and preconditioning matrices
  Mat :: petsc_pc_mat, petsc_a_mat
  
  ! Matrix structure
  MatStructure :: petsc_str

  ! This array is used extensively within this module. To avoid over-head in the
  ! allocation we preallocate memory inside the "setup_petsc_solver" routine.
  integer, allocatable, dimension(:) :: index_list

contains

  !============================================================================

  !============================================================================
  !
  ! start_petsc_perf_mon - Register PETSc event for logging operation. 
  !                    
  subroutine star_petsc_perf_mon()

    ! Nothing is implicitly defined
    implicit none
    
    ! PETSc variables
    PetscErrorCode :: i_err

    ! Register event for logging operation 
    call PetscLogEventRegister('rhs',0,petsc_user_event,i_err)
    call check_petsc_error(i_err)

    ! Print message at screen 
    write(*,*) 'Event registered'

    ! Login the beginning of a user event 
    call PetscLogEventBegin(petsc_user_event,0,0,0,0,i_err)
    call check_petsc_error(i_err)

    ! Print message at screen
    write(*,*) 'Event begins'

    return
  end subroutine star_petsc_perf_mon

  !============================================================================

  !============================================================================
  !
  ! finish_petsc_perf_mon - Counts user's flops and ends PETSc event. 
  !             
  subroutine finish_petsc_perf_mon()

    ! Nothing is implicitly defined
    implicit none
    
    ! PETSc variables
    PetscErrorCode :: i_err
    
    ! Code variables
    real(wp) :: user_flops

    ! Add floating point operation to the global counter
    call PetscLogFlops(user_flops,i_err); call check_petsc_error(i_err);

    ! Login the end of a user event
    call PetscLogEventEnd(petsc_user_event,0,0,0,0,i_err)
    call check_petsc_error(i_err)

    ! Print message at screen
    write(*,*) 'Event ends'

    return
  end subroutine finish_petsc_perf_mon

  !============================================================================

  !============================================================================
  !
  ! setup_petsc_solver - Creates and sets up the nonlinear solver in PETSc. 
  !
  subroutine setup_petsc_solver()

    ! Load modules
    use referenceVariables
    use controlvariables
    use CSRlocalvariables

    ! Nothing is implicitly defined
    implicit none

    ! PETSc variables
    PetscErrorCode :: i_err

    ! Code variables
    integer :: i_brow
    integer :: block_size
    integer :: iell, ielh
    integer, allocatable, dimension(:) :: d_nnz, o_nnz
    real(wp) :: d_tol_ksp
    integer :: max_iters_ksp

    ! Allocate memory for local PETSc vectors
    call VecCreate(PETSC_COMM_WORLD,petsc_x_vec,i_err);
    call check_petsc_error(i_err)
    call VecSetSizes(petsc_x_vec,nprows,ngrows,i_err)
    call check_petsc_error(i_err)
    call VecSetFromOptions(petsc_x_vec,i_err)
    call check_petsc_error(i_err)
    call VecSet(petsc_x_vec,0.0_wp,i_err)
    call check_petsc_error(i_err)
    call VecDuplicate(petsc_x_vec,petsc_f_vec,i_err)
    call check_petsc_error(i_err)
    call VecAssemblyBegin(petsc_x_vec,i_err)           
    call check_petsc_error(i_err)
    call VecAssemblyEnd(petsc_x_vec,i_err)
    call check_petsc_error(i_err)
    call VecAssemblyBegin(petsc_f_vec,i_err)
    call check_petsc_error(i_err)
    call VecAssemblyEnd(petsc_f_vec,i_err)
    call check_petsc_error(i_err)

    ! Set linear and nonlinear solver options
    if (IMEX_Element == 'explicit' .and. IMEX_Penalty == 'explicit') then
      call PetscOptionsSetValue('-pc_type','none',i_err)
      call check_petsc_error(i_err)
    else 
      call PetscOptionsSetValue('-pc_type','asm',i_err)
      call PetscOptionsSetValue('-sub_pc_type','ilu',i_err)
      call PetscOptionsSetValue('-sub_pc_factor_levels','1',i_err)
    end if

    ! Set SNES tolerances and maximum number of function evaluation
    call PetscOptionsSetValue('-snes_max_funcs','200',i_err)
    call check_petsc_error(i_err)
    call PetscOptionsSetValue('-snes_atol','1e-8',i_err)
    call check_petsc_error(i_err)
    call PetscOptionsSetValue('-snes_rtol','1e-25',i_err)
    call check_petsc_error(i_err)
    call PetscOptionsSetValue('-snes_stol','1e-25',i_err)
    call check_petsc_error(i_err)

    ! Set matrix free operator    
    call PetscOptionsSetValue('-snes_mf_operator','',i_err)
    call check_petsc_error(i_err)
    call PetscOptionsSetValue('-mat_mffd_type','wp',i_err)
    call check_petsc_error(i_err)
    call PetscOptionsSetValue('-mat_mffd_compute_normu','',i_err)
    call check_petsc_error(i_err)

    ! Create nonlinear solver context
    call SNESCreate(PETSC_COMM_WORLD,petsc_snes,i_err); 
    call check_petsc_error(i_err)

    ! Set nonlinear function    
    call SNESSetFunction(petsc_snes,petsc_f_vec,form_function_snes, &
      & PETSC_NULL_OBJECT,i_err)
    call check_petsc_error(i_err)
    
    ! Set up matrix-free Jacobian
    call MatCreateSNESMF(petsc_snes,petsc_a_mat,i_err)
    call check_petsc_error(i_err)

    ! Local portion of the residual Jacobi matrix
    call MatCreate(PETSC_COMM_WORLD,petsc_pc_mat,i_err)
    call check_petsc_error(i_err)

    !call MatSetType(petsc_pc_mat,MATMPIBAIJ,i_err)
    call MatSetType(petsc_pc_mat,MATBAIJ,i_err)
    call check_petsc_error(i_err)

    call MatSetSizes(petsc_pc_mat,nprows,nprows,ngrows,ngrows,i_err)
    call check_petsc_error(i_err)

    block_size = 5
    call MatSetBlockSize(petsc_pc_mat, block_size) 
    call check_petsc_error(i_err)
    
    call MatSetFromOptions(petsc_pc_mat,i_err)
    call check_petsc_error(i_err)

    call MatSetUp(petsc_pc_mat,i_err)
    call check_petsc_error(i_err)

    ! Option 1: let PETSc compute it with a finite differencing approach
    ! without taking advantage of the sparsity. This is very very slow and 
    ! expensive.
    ! --------------------------------------------------------------------
    ! call SNESSetJacobian(petsc_snes,petsc_jacobi_mat,petsc_jacobi_mat, &
    !                    & SNESComputeJacobianDefault,PETSC_NULL_OBJECT,i_err)
    ! call check_petsc_error(i_err)


    ! Option 2: provide the function for computing the Jacobian
    ! ---------------------------------------------------------
    ! Pre-allocate storage for Jacobi matrix.
    ! Here the pre-allocation is perfect because we pre-allocate exactly the 
    ! number of nonzero block terms per row, corresponding to both diagonal and 
    ! off-diagonal submatrices. By doing that it is possible to improve 
    ! drastically the performance of the code.
    !
    ! d_nnz = array containing the number of block nonzeros in the various 
    !         block rows of the in diagonal portion of the local submatrix 
    ! o_nnz = array containing the number of nonzeros in the various block rows 
    !         of the off-diagonal portion of the local submatrix
    allocate(o_nnz(size(iaS)-1),d_nnz(size(iaS)-1))

    ! Initialize the number of nonzero blocks
    d_nnz = 0
    o_nnz = 0

    ! Compute the number of nonzero blocks of the local submatrix
    do i_brow = 1, size(iaS)-1
      d_nnz(i_brow) = iaS(i_brow+1) - iaS(i_brow)
    enddo

    ! Preallocate memory for the local submatrix
    call MatMPIBAIJSetPreallocation(petsc_pc_mat,block_size,0,d_nnz,0,o_nnz, &
      & i_err)
    call check_petsc_error(i_err)

    ! Use column oriented layout
    call MatSetOption(petsc_pc_mat,MAT_ROW_ORIENTED,PETSC_FALSE,i_err)
    call check_petsc_error(i_err)

    ! Deallocate memory for d_nnz and o_nnz
    deallocate(d_nnz)
    deallocate(o_nnz)

    ! Residual Jacobi matrix is equal to preconditioning matrix
    !petsc_jacobi_mat = petsc_pc_mat

    ! Set function for the calculation of the residual Jacobi matrix
    call SNESSetJacobian(petsc_snes,petsc_a_mat,petsc_pc_mat, &
      & form_jacobian_matrix,PETSC_NULL_OBJECT,i_err)
    call check_petsc_error(i_err)

    ! Set SNES options
    ! Recompute the Jacobian matrix every N iterations
    !call SNESSetLagJacobian(petsc_snes,-2)
    !call SNESSetLagPreconditioner(petsc_snes,1)

    call SNESKSPSetUseEW(petsc_snes,PETSC_TRUE,i_err)
    call check_petsc_error(i_err)

    ! Get linear solver 
    call SNESGetKSP(petsc_snes,petsc_ksp,i_err)
    call check_petsc_error(i_err)
    
    ! Set KSP to GMRES
    call KSPSetType(petsc_ksp,KSPGMRES,i_err) ! another option is KSPLGMRES
    call check_petsc_error(i_err)

    ! Set some KSP parameter values
    ! The divergence tolerance (amount residual can increase before 
    ! KSPDefaultConverged() concludes that the method is diverging) 
    d_tol_ksp = 1.0e01_wp
    
    ! Maximum number of iterations to use 
    max_iters_ksp = 200
    call KSPSetTolerances(petsc_ksp,PETSC_DEFAULT_DOUBLE_PRECISION, &
      & PETSC_DEFAULT_DOUBLE_PRECISION,d_tol_ksp,max_iters_ksp,i_err)
    call check_petsc_error(i_err)
    
    ! Get preconditioner 
    call KSPGetPC(petsc_ksp,petsc_pc,i_err)
    call check_petsc_error(i_err)

    ! Set preconditioner options
    call PCSetFromOptions(petsc_pc,i_err)
    call check_petsc_error(i_err)
    
    ! Set the overlap between a pair of subdomains for the additive Schwarz 
    ! preconditioner
    call PCASMSetOverlap(petsc_pc,2,i_err)
    call check_petsc_error(i_err)

    ! Set various SNES and KSP parameters from user options
    call SNESSetFromOptions(petsc_snes,i_err)
    call check_petsc_error(i_err)

    return
  end subroutine setup_petsc_solver

  !============================================================================

  !============================================================================
  !
  ! form_jacobian_matrix - Forms the action of the Jacobian matrix on a vector
  !                        and set the preconditioning matrix equal to the 
  !                        analytical residual Jacobian, i.e. PC = dRHS/dCons.
  !
  ! Input parameters:
  ! snes  - the SNES context
  ! x_vec - input vector
  ! dummy - optional user-defined context (not used here)
  !
  ! Output parameter:
  ! jacobi_mat - residual Jacobi matrix
  ! pc_mat - preconditioning matrix
  ! flag - flag indicating matrix structure 

  subroutine form_jacobian_matrix(snes,x_vec,a_mat,pc_mat,flag,dummy,i_err)
    
    ! Load modules
    use controlvariables
    use referencevariables 
    use time_integ_coeff
    use CSRlocalvariables
    use jacobian_matrix_implicit_ts, only : compute_residual_jacobian_element

    ! Nothing is implicitly defined
    implicit none

    ! PETSc variables
    SNES, intent(in) :: snes
    Vec, intent(in) :: x_vec
    Mat, intent(out) :: a_mat, pc_mat
    MatStructure, intent(out) :: flag
    PetscObject, intent(in):: dummy
    PetscErrorCode, intent(inout) :: i_err

    double  precision info(MAT_INFO_SIZE)
    double  precision mal, nz_a
    
    ! Code variables
    integer :: iell, ielh
    integer :: shift
    integer :: i_node, i_elem
    integer :: nodes_counter, blocks_counter
    integer :: nnz_blocks_per_block_row, i
    real(wp) :: a_kk
    logical :: assemble_jacobi_mat 

    ! Get the diagonal coefficient of the IMEX scheme
    a_kk = arkimp(current_stage_imex,current_stage_imex)

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    ielh = ihelems(2)

    ! Compute shift of the diagonal portion of each process using C indexing
    shift = (iell-1)*nodesperelem-1

    ! Initialize node_counter
    nodes_counter = 1

    ! ON ELEMENT contribution to the residual Jacobi submatrix
    ! ========================================================
    do i_elem = iell, ielh

      ! Compute element contribution
      call compute_residual_jacobian_element(timestep,a_kk,i_elem)

      ! Initialize element block counter
      blocks_counter = 1

      do i_node = 1, nodesperelem

        ! Number of nonzero blocks per block row
        nnz_blocks_per_block_row = iaS(nodes_counter + 1) - iaS(nodes_counter)

        call MatSetValuesBlocked(pc_mat,1,shift+nodes_counter, &
             & nnz_blocks_per_block_row, &
             & jaS(iaS(nodes_counter):iaS(nodes_counter+1)-1)+shift, &
             & dfdu_a_elem(:,:,blocks_counter:nnz_blocks_per_block_row), &
             & INSERT_VALUES,i_err)
        call check_petsc_error(i_err)
        
        ! Update counters
        blocks_counter = blocks_counter + nnz_blocks_per_block_row
        nodes_counter = nodes_counter + 1
    
      enddo
  
      ! Check Jacobian matrix of the element
      !call check_jacobian_matrix_element(i_elem)
      !stop

    enddo

    ! Parallel assembly of the preconditioning matrix
    call MatAssemblyBegin(pc_mat,MAT_FINAL_ASSEMBLY,i_err)
    call check_petsc_error(i_err)
    call MatAssemblyEnd(pc_mat,MAT_FINAL_ASSEMBLY,i_err)
    call check_petsc_error(i_err)
    call MatSetOption(pc_mat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,i_err)
    call check_petsc_error(i_err)
    call MatSetOption(pc_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,i_err)
    call check_petsc_error(i_err)

    ! Parallel assembly of the residual Jacobian matrix if necessary
    !if(a_mat .ne. pc_mat ) then
    !  call MatAssemblyBegin(a_mat,MAT_FINAL_ASSEMBLY,i_err)
    !  call check_petsc_error(i_err)
    !  call MatAssemblyEnd(a_mat,MAT_FINAL_ASSEMBLY,i_err)
    !  call check_petsc_error(i_err)
    !endif

    call MatMFFDComputeJacobian(snes,x_vec,a_mat,pc_mat,PETSC_NULL_OBJECT, &
      & PETSC_NULL_OBJECT,i_err)
    call check_petsc_error(i_err)

    ! Print at screen the pc_mat matrix assembled by PETSc 
    !call MatView(pc_mat,PETSC_VIEWER_STDOUT_SELF, i_err)
    !call check_petsc_error(i_err)

    ! Print at screen some efficiency information
    !call MatGetInfo(pc_mat,MAT_LOCAL,info,i_err)
    !call check_petsc_error(i_err)
    !mal = info(MAT_INFO_MALLOCS)
    !nz_a = info(MAT_INFO_NZ_ALLOCATED)
    !write(*,*) 'Additional calloc calls', mal
    !write(*,*) 'Number of nonzero', nz_a

    ! SAT contribution to the residual Jacobi matrix
    ! ==============================================
    ! Here we will call MatSetValuesBlocked() with the the flag ADD_VALUES

    return
  end subroutine form_jacobian_matrix

  !============================================================================

  !============================================================================
  !
  ! form_initial_guess - Forms initial approximation
  !
  ! Input parameters:
  ! x_vec - input vector
  !
  ! Output parameter:
  ! x_vec - function vector

  subroutine form_initial_guess(snes,x_vec)
    use variables
    use referencevariables
    use mpimod
    use CSRlocalvariables
    use nsereferencevariables

    ! PETSc variables
    SNES :: snes
    Vec :: x_vec
    PetscErrorCode :: i_err

    ! Code variables
    integer :: m, i_elem, i_node, i_eq
    integer ::  iell, ielh
    integer :: len_local
    real(wp) :: error
    ! Temporary variable for calculating the error
    !real(wp), allocatable, dimension(:,:) :: tmp_x


    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    ielh = ihelems(2)

    ! Number of local solution values and indices
    len_local = (ielh-iell+1)*nodesperelem*nequations

    !allocate(tmp_x(len_local,1))

    ! Store sequence of integer numbers in index_list
    index_list = (/(m, m=1, len_local)/)

    ! Add contribution previous elements/nodes and use C indexing
    index_list = index_list + (iell-1)*nodesperelem*nequations - 1

    !tmp_x = reshape(ug,(/len_local,1/))

    ! Set values in the portion of the vector owned by the process
    call VecSetValues(x_vec,len_local,index_list,reshape(ug,(/len_local,1/)),INSERT_VALUES,i_err) 
    call check_petsc_error(i_err)

    ! Parallel assembly of initial guess
    call VecAssemblyBegin(x_vec,i_err)
    call check_petsc_error(i_err)
    call VecAssemblyEnd(x_vec,i_err)
    call check_petsc_error(i_err)

    ! Deallocate memory
    !deallocate(tmp_x)

    ! Check if reshape is working correctly
    !error = 0.0_wp
    !do i_elem = iell, ielh
    !  do i_node = 1,nodesperelem
    !    do i_eq = 1,nequations
    !      error = error + ug(i_eq,i_node,i_elem) - tmp_x((i_node-1)*5+(i_elem-iell)*nodesperelem*5+i_eq,1) 
    !    enddo
    !  enddo
    !enddo
    !
    !write(*,*) 'Error in reshape', error

    return
  end subroutine form_initial_guess

  !============================================================================

  !============================================================================
  !
  ! form_function_snes - Evaluates nonlinear function, F(x).
  !
  ! Input parameters:
  ! snes  - the SNES context
  ! x_vec - input vector
  ! dummy - optional user-defined context, as set by SNESSetFunction()
  !            (not used here)
  !
  ! Output parameter:
  ! f_vec - function vector

  subroutine form_function_snes(snes,x_vec,f_vec,dummy,i_err)
    use variables, only : ug, non_lin_res
    use controlvariables
    use referencevariables 
    use time_integ_coeff
    use implicit_residual, only: compute_implicit_residual_imex

    ! PETSc variables
    SNES, intent(in) :: snes
    Vec, intent(in):: x_vec
    Vec, intent(out) :: f_vec
    PetscObject, intent(in) :: dummy
    PetscErrorCode, intent(inout) :: i_err
    PetscScalar, pointer, dimension(:) :: x_vec_pointer

    ! Code variables  
    integer :: m
    integer :: iell, ielh
    integer :: len_local

    ! Low volumetric element index
    iell = ihelems(1)

    ! Low volumetric element index
    ielh = ihelems(2)

    ! Number of local solution values and indices
    len_local = (ielh-iell+1)*nodesperelem*nequations

    ! Store sequence of integer numbers in index_list
    index_list =(/(m, m=1, len_local)/)

    ! Add contribution of previous elements/nodes and use C indexing
    index_list = index_list + (iell-1)*nodesperelem*nequations - 1

    ! Get local vector, i.e. portion of the vector owned by the processor
    call VecGetArrayF90(x_vec,x_vec_pointer,i_err)

    ! Set solution value based on input vector x_vec
    ! Actually the portion of the global vector owned by the process
    ug = reshape(x_vec_pointer,(/nequations,nodesperelem,ielh-iell+1/))

    ! Release x_vec_pointer pointer
    call VecRestoreArrayF90(x_vec,x_vec_pointer,i_err)
    call check_petsc_error(i_err)

    ! Deallocate x_vec_pointer pointer
    if(associated(x_vec_pointer)) deallocate(x_vec_pointer)

    ! Calculate total nonlinear residual for the implicit IMEX RK step
    call compute_implicit_residual_imex(current_stage_imex,timestep)

    ! Set values in the portion of the vector owned by the process
    !call VecSetValues(f_vec,len_local,index_list,reshape(non_lin_res,(/len_local,1/)),INSERT_VALUES,& 
    !  & i_err)
    call VecSetValues(f_vec,len_local,index_list,non_lin_res,INSERT_VALUES,& 
      & i_err)

    ! Parallel assembly of the nonlinear residual
    call VecAssemblyBegin(f_vec,i_err)
    call check_petsc_error(i_err)
    call VecAssemblyEnd(f_vec,i_err)
    call check_petsc_error(i_err)

    return
  end subroutine form_function_snes

  !============================================================================

  !============================================================================
  !
  ! solve_implicit_petsc - Solves nonlinear system of equations with PETSc SNES. 
  !
  ! Output parameters:
  ! converged - whether SNES has converged

  subroutine solve_implicit_petsc(converged)
    use variables, only : ug
    use referenceVariables
    use controlvariables
    use time_integ_coeff
    use implicit_residual, only: rk_imex_svp

    ! PETSc variables
    PetscErrorCode :: i_err
    PetscInt :: newton_iter, linear_iter
    SNESConvergedReason :: conv
    
    ! Code variables
    logical, intent(out) :: converged
    integer :: m
    integer :: iell, ielh
    real(wp), pointer, dimension(:) :: x_vec_pointer
    integer :: len_local

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    ielh = ihelems(2)

    ! Number of local solution values and indices
    len_local = (ielh-iell+1)*nodesperelem*nequations

    ! Store sequence of integer numbers in index_list
    index_list = (/(m, m=1, len_local)/)

    ! Add contribution previous elements/nodes and use C indexing
    index_list = index_list + (iell-1)*nodesperelem*nequations - 1

    ! Call RK-IMEX SVP routine
    call rk_imex_svp(current_stage_imex,timestep)

    ! Set values in the portion of the vector owned by the process
    !call VecSetValues(petsc_x_vec,len_local,index_list, &
    ! & reshape(ug,(/len_local,1/)),INSERT_VALUES,i_err) 
    call VecSetValues(petsc_x_vec,len_local,index_list, &
      & ug(:,:,:),INSERT_VALUES,i_err) 
    call check_petsc_error(i_err)

    ! Parallel assembly of initial guess
    call VecAssemblyBegin(petsc_x_vec,i_err)
    call check_petsc_error(i_err)
    call VecAssemblyEnd(petsc_x_vec,i_err)
    call check_petsc_error(i_err)

    ! Solve nonlinear system
    call SNESSolve(petsc_snes,PETSC_NULL_OBJECT,petsc_x_vec,i_err)
    call check_petsc_error(i_err)

    ! Check if the solver has converged
    call SNESGetConvergedReason(petsc_snes,conv,i_err)
    call check_petsc_error(i_err)

    ! If SNES did not converge write a message at screen
    if(conv<0) then
      write(*,*) 'SNES did not converge', conv, myprocid
      ! This logical variable is passed back to the timeinteg.f90 so that 
      ! a smaller time step is taken      
      converged = .false.
      return
    end if

    ! Set logical variable to .true. (the above if condition was not 
    ! satisfied)
    converged = .true.

    ! Get number of nonlinear (Newton) iterations
    call SNESGetIterationNumber(petsc_snes,newton_iter,i_err)
    call check_petsc_error(i_err)

    ! Get number of linear iterations
    call SNESGetLinearSolveIterations(petsc_snes,linear_iter,i_err)
    call check_petsc_error(i_err)

    ! Write at screen the number of Newton and linear iterations
    if(myprocid .eq.0) then
      write(*,100) "RK stage = ",current_stage_imex,"nonlinear iterations = ",newton_iter, &
        & 'linear iterations = ',linear_iter
    endif
    100 format(A,i1.1,1X,A,i5.1,1X,A,i5.1,1X,i2.1)

    ! Get local vector, i.e. portion of the vector owned by the processor
    call VecGetArrayF90(petsc_x_vec,x_vec_pointer,i_err)
    call check_petsc_error(i_err)

    ! Assign to the solution array ug(:,:,:) the new computed solution
    ug = reshape(x_vec_pointer,(/nequations,nodesperelem,ielh-iell+1/))

    ! Release x_vec_local pointer
    call VecRestoreArrayF90(petsc_x_vec,x_vec_pointer,i_err)
    call check_petsc_error(i_err)

    ! Deallocate x_vec_pointer pointer
    if(associated(x_vec_pointer)) deallocate(x_vec_pointer)

    return
  end subroutine solve_implicit_petsc 

  !============================================================================

  !============================================================================
  ! check_petsc_error - Checks the content of i_err for possible PETSc errors
  !
  ! Input parameters:
  ! i_err - output of a call to a PETSc routine

  subroutine check_petsc_error(i_err)

    ! Nothing is implicitly defined
    implicit none

    ! PETSc variables
    PetscErrorCode, intent(in) :: i_err 

    ! Print error message at screen
    if (i_err /= 0) then
      write(*,*) 'PetscError: ',i_err
      stop
    end if

    return
  end subroutine check_petsc_error

  !============================================================================

  !============================================================================
  !
  ! finalize_petsc_solver - Finalises PESTc solver, i.e. destroy vectors and 
  !                         matrices and deallocate memory

  subroutine finalize_petsc_solver()

    ! PETSc variables
    PetscErrorCode :: i_err 

    ! Deallocate PETSc vectors and matrices 
    call VecDestroy(petsc_x_vec,i_err)
    call check_petsc_error(i_err)

    call VecDestroy(petsc_f_vec,i_err)
    call check_petsc_error(i_err)

    call MatDestroy(petsc_a_mat,i_err)
    call check_petsc_error(i_err)

    call MatDestroy(petsc_pc_mat,i_err)
    call check_petsc_error(i_err)

    call SNESDestroy(petsc_snes,i_err)
    call check_petsc_error(i_err)

    ! Deallocate other vectors and matrices
    deallocate(index_list)

    ! call PetscLogPrintSummary(PETSC_COMM_WORLD,'petscsummary.log',i_err)

    return
  end subroutine finalize_petsc_solver

  !============================================================================

  !============================================================================
  !
  ! check_jacobian_matrix_element - Verifies the correctnes of the element-wise
  !                                 Jacobian matrix.
  !
  subroutine check_jacobian_matrix_element(elem_id)
    
    ! Load modules
    use controlvariables
    use referencevariables 
    use time_integ_coeff
    use CSRlocalvariables
    use variables
    use unary_mod
    use navierstokes

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: elem_id
    integer :: length
    real(wp), allocatable, dimension(:) :: w_elem, u_elem, div_fv_elem, res
    real(wp), allocatable, dimension(:,:) :: time_der
    real(wp), allocatable, dimension(:,:) :: tmp_sol
    real(wp), allocatable, dimension(:) :: sol_node, v_i_node

    ! Code variables
    integer :: iell, ielh
    integer :: i_node, j_dir, i
    real(wp) :: a_kk

    logical :: testing_inviscid
    logical :: testing_2nd_visc_c

    testing_inviscid = .false.
    testing_2nd_visc_c = .false.

    ! Get the diagonal coefficient of the IMEX scheme
    a_kk = arkimp(current_stage_imex,current_stage_imex)

    ! Low volumetric element index
    iell = ihelems(1)

    ! High volumetric element index
    ielh = ihelems(2)

    ! Length of the solution vector for one element
    length = nequations*nodesperelem


    if (testing_2nd_visc_c) then
      ! Allocate memory 
      allocate(w_elem(length))
      w_elem = 0.0_wp

      allocate(u_elem(length))
      u_elem = 0.0_wp

      allocate(div_fv_elem(length))
      div_fv_elem = 0.0_wp

      allocate(res(length))
      res = 0.0_wp

      allocate(time_der(1:nequations,nodesperelem))
      time_der = 0.0_wp

      allocate(tmp_sol(1:nequations,nodesperelem))
      tmp_sol = 0.0_wp

      allocate(sol_node(nequations))
      sol_node = 0.0_wp

      allocate(v_i_node(nequations))
      v_i_node = 0.0_wp


      ! Compute dfdu_a_elem times u_elem
      do i_node =1, nodesperelem
        call conserved_to_primitive(ug(:,i_node,elem_id),v_i_node(:),nequations)

        sol_node = wg(:,i_node,elem_id)

        sol_node = matmul(dudw(v_i_node,nequations),sol_node)

        tmp_sol(:,i_node) = sol_node
      enddo

      ! Reshape tmp_sol for matrix product
      w_elem = reshape(tmp_sol,(/length/))
        
      call a_bl_mu_x(nequations,nodesperelem,w_elem,div_fv_elem,dfdu_a_elem, &
        & ja_elem,ia_elem)

      divf = 0.0_wp

!     HACK
!     call Flux_Divergence(timestep,elem_id)
!     HACK

      ! Loop over all nodes in the element
      do i_node = 1, nodesperelem
        ! Add the contribution from the flux divergence in each direction
        do j_dir = 2,2
          time_der(:,i_node) = time_der(:,i_node) + divf(:,j_dir,i_node,elem_id)/Jx_r(i_node, elem_id)
        end do
      end do

      !res = a_kk*timestep*reshape(time_der,(/length/)) 
      res = reshape(time_der,(/length/)) 
      
      write(*,*) 'max val res', maxval(abs(res-div_fv_elem))
      do i = 1,length
        write(*,*) i, res(i),div_fv_elem(i)
      enddo

      stop

    endif

    return
  end subroutine check_jacobian_matrix_element

  !============================================================================


end module petsc_snes_solver



