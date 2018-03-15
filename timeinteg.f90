! This module contains all the subroutines and functions that drives the
! integration in time of the nonlinear system of ODEs arising from the spatial
! discretization of the PDEs.

module timeinteg
  
  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

contains

  !============================================================================

  !============================================================================
  ! Runge_Kutta_init - Set the Runge-Kutta scheme using the options passed in 
  !                    the SSDCstartup file and initialize its coefficients. 
  
  subroutine Runge_Kutta_init()
    
    ! Load modules
    use controlvariables
    use time_integ_coeff 
    
    ! Nothing is implicitly defined
    implicit none

    ! Set the RK method based on the input parameters given in the SSDCstartup
    ! file
    select case (RK_Method)
    case ('Williamson_Low_Storage_45')
      call ls_rk_initialization()

    case ('Kraaij_LS_RK_35')
      call ls_TVDrk_initialization()

    case ('heun_method')

    case ('IMEX_RK_34')

    case ('IMEX_RK_46')
      call rk_imex_46_initialization()

    case ('Steady_State')
    case ('Frequency_Domain')
    case ('Kraai')
    case default
      write(*,*) 'runge_kutta_init'
      write(*,*)'Not a valid Temporal Scheme'
      write(*,*)'Check RK_Method in setup file'
      write(*,*)'Stopping'
      stop
    end select

    return
  end subroutine Runge_Kutta_init

  !============================================================================

  !============================================================================
  ! Physics_Solver - Calls the subroutine that actually drives the time
  !                  integration of the ODEs system. 

  subroutine Physics_Solver()

    ! Load modules
    use controlvariables

    ! Nothing is implicitly defined
    implicit none

    ! Call the appropriate subroutine to drive the time integration
    select case (RK_Method)
    case ('Williamson_Low_Storage_45')
      call LowStorage_RK() 

    case ('Kraaij_LS_RK_35')
      call LowStorage_RK() 

    case ('heun_method')
      call heun_method()

    case ('IMEX_RK_34')
      call rk_imex()
    
    case ('IMEX_RK_46')
      call rk_imex()
    
    case ('Steady_State')
    case ('Frequency_Domain')
    case ('Kraai')
    case default
      write(*,*) 'Physics_Solver'
      write(*,*)'Not a valid Temporal Scheme'
      write(*,*)'Check RK_Method in setup file'
      write(*,*)'Stopping'
      stop
    end select

    return
  end subroutine Physics_Solver

  !============================================================================

  !============================================================================
  ! LowStorage_RK - Advances the solution in time using a pre-determined 
  !                 low-storage explicit RK schemes.

  subroutine LowStorage_RK()
    
    ! Load modules
    use controlvariables
    use referencevariables
    use physics_driver
    use errorestimation
    use navierstokes
    use mpimod
    use variables, only:  ug, uhat, uold, du
    use time_integ_coeff 
    use dkinetic_energy_dt_enstrophy
    use tools_IO

    implicit none                            ! Nothing is implicitly defined

    integer :: irkstep
    integer :: i_err
!   character(len=300) :: time_space_error_file_name
!   integer :: i_unit, io_status
!   logical :: negTemp = .false.
!   character(120) :: message

    continue

    ! Set time step counter to zero
    itimestep = 0

    ! Compute primitive variables, entropy variables and CNG gradients
    call nse_reconcilestates()
    ! Compute interesting quantities and write output if needed
!   call post_process_ts_0(itimestep,timeglobal)

    ! Time-step loop
    timeloop:do
      ! Check whether the global time has hit maximum
      if (timeglobal >= timemaximum) exit

      ! Update time-step index counter
      itimestep = itimestep + 1

      ! Limit time-step such that we don't exceed max
      timestep = min(timestep,timemaximum-timeglobal)

      ! Reset du (note that this is for LSRK)
      du = 0.0_wp

      ! Set uhat and uold to current solution
      uhat = ug
      uold = ug 

      ! Advance solution one timestep using a predetermined RK scheme
      rkloop:do irkstep = 1,rksteps 

        ! check for negative temperatures 
!       if (checkNegTemp(negTemp) < 0) then ! (mpimod)
!         itimestep = itimestep - 1 ; timestep = 0.5_wp*timestep ; ug = uold ;
!         cycle timeloop
!       end if
        ! update stage time
        timelocal = timeglobal + crk(irkstep)*timestep

!       call Calc_Entropy_Viscosity()   !  Decodes dudt_S into mut

        ! calculate time derivative 
        call physics_timederivative()

        if (irkstep == 1) then
          if (write_dke_dt) then 
            call compute_dkinetic_energy_dt()
          end if
        end if

        ! Integrate in time
        call LSRK(irkstep)

        ! Update primitive, entropy variables and CNG gradients
        call nse_reconcilestates()

      end do rkloop

      ! Update global time
      timeglobal = timeglobal + timestep

      ! Calculate embedded temporal error estimate
      call calcembeddedtemporalerror() 

      ! Calculate embedded spatial error estimate
      call calcembeddedspatialerror()  

      ! Evaluate next dt
      call timestepcontroller(itimestep)

      ! Compute time-averaged quantities and write output if needed
      if(.not. non_conforming) call post_process_ts(itimestep,timeglobal)

      ! Call Solution filter to clip highest mode
      if (filter_solution) call Solution_Filter()

      ! Output to file this timestep
      if (myprocid == 0) then

        ! If verbose is set to true then output more at screen
        if (verbose) then
        
          ! Write meaning of variables
          if(mod(itimestep-1,20) == 0 ) write(*,100)

          ! Write values
          write(*,101) itimestep, timeglobal, timestep, err_time_lf,  &
            & err_space_lf
        endif
      endif

    end do timeloop

    ! Format for writing at screen
    100 format('     step   |&
               &          timeG          |&
               &            dt           |&
               &        Linf_time        |&
               &        Linf_space ')
!    101 format(i8,3x, 4(ES12.5,1x))
    101 format(i8,3x, 4(ES25.15,1x))
    
    ! Close unit used for max error norm 
!   close(unit=i_unit)

    ! Wait for other processes
    call mpi_barrier(PETSC_COMM_WORLD,i_err)

    return
  end subroutine LowStorage_RK

  !============================================================================

  !============================================================================
  ! heun_method - Advances the solution in time using Heun's method (2nd order
  ! RK method) with the embedded Euler method

  subroutine heun_method()
    
    ! Load modules
    use controlvariables
    use referencevariables
    use physics_driver
    use errorestimation
    use navierstokes
    use mpimod
    use variables, only:  ug, uhat, uold, du, dudt
    use time_integ_coeff 
    use dkinetic_energy_dt_enstrophy
    use tools_IO

    ! Nothing is implicitly defined
    implicit none

!   character(len=300) :: time_space_error_file_name
    integer :: i_err
!   integer :: i_unit
!   character(120) :: message

    continue

    ! Set time step counter to zero
    itimestep = 0

    ! Compute primitive variables, entropy variables and CNG gradients
    call nse_reconcilestates()

    ! Compute interesting quantities and write output if needed
    call post_process_ts_0(itimestep,timeglobal)

    ! Time-step loop
    do
      ! Check whether the global time has hit maximum
      if (timeglobal >= timemaximum) exit

      ! Update time-step index counter
      itimestep = itimestep + 1

      ! Limit time-step such that we don't exceed max
      timestep = min(timestep,timemaximum-timeglobal)

      ! Time at the first stage (c_1 = 0)
      timelocal = timeglobal

      ! Compute dke/dt
      call compute_dkinetic_energy_dt()
      
      ! Set uhat and uold to current solution
      uhat = ug
      uold = ug 

      ! Calculate time derivative at t = t_n
      call physics_timederivative()

      ! First stage (intermediate value)
      ug(:,:,:) = uold(:,:,:) + timestep*dudt(:,:,:)
      uhat = ug

      ! Store dudt for the final update
      du(:,:,:) = dudt(:,:,:)

      ! Update primitive, entropy variables and CNG gradients
      call nse_reconcilestates()

      ! Time at the second stage (c_2 = 1)
      timelocal = timeglobal + timestep

      ! Calculate time derivative at t_n+1 using the intermediate value for the
      ! solution
      call physics_timederivative()

      ! Update solution
      ug(:,:,:) =  uold(:,:,:) + timestep/2.0_wp*(du(:,:,:) + dudt(:,:,:))

      ! Update primitive, entropy variables and CNG gradients
      call nse_reconcilestates()

      ! Update global time
      timeglobal = timeglobal + timestep

      ! Calculate embedded temporal error estimate
      call calcembeddedtemporalerror() 

      ! Calculate embedded spatial error estimate
      call calcembeddedspatialerror()  

      ! Evaluate next dt
      call timestepcontroller(itimestep)

      ! Compute time-averaged quantities and write output if needed
      call post_process_ts(itimestep,timeglobal)

      ! Call Solution filter to clip highest mode
      if (filter_solution) then
        call Solution_Filter()
      endif

      ! Output to file this timestep
      if (myprocid == 0) then
!        write(i_unit,101) itimestep, timeglobal, timestep, err_time_lf, &
!          & err_space_lf

        ! If verbose is set to true then output more at screen
        if (verbose) then
        
          ! Write meaning of variables
          if(mod(itimestep-1,20) == 0 ) then
            write(*,100)
          !  write(*,*) 'here', myprocid 
          endif

          ! Write values
          write(*,101) itimestep, timeglobal, timestep, err_time_lf,  &
            & err_space_lf
        endif
      endif

    end do

    ! Format for writing at screen
    100 format('    step   |   timeG    |     dt     |  Linf_time | Linf_space ')
!    101 format(i8,3x, 4(ES12.5,1x))
    101 format(i8,3x, 4(ES25.15,1x))
    
    ! Close unit used for max error norm 
!   close(unit=i_unit)

    ! Wait for other processes
    call mpi_barrier(PETSC_COMM_WORLD,i_err)

    return
  end subroutine heun_method

  !============================================================================

  !============================================================================
  ! rk_imex_middle - Computes the explicit portion of the explicit RK IMEX
  !                  scheme, i.e. computes uexp.
  !
  ! Input parameters:
  ! istep - current stage.
  ! dt - current time-step.

  subroutine rk_imex_middle(istep,dt)  

    ! Load modules
    use variables
    use time_integ_coeff 

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: istep
    real(wp), intent(in) :: dt
    integer :: i

    ! set uexp equal to uold, i.e. solution at the current time level
    uexp = uold

    ! Calculate residual contribution from previous stages, i.e. "explicit" 
    ! part of u^k: 
    ! uexp = u^k = u^n + sum_j^{k-1} (a_{kj} F(u^j) + \hat{a}_{kj} G(u^j)
    do i = 1,istep-1
      uexp(:,:,:) = uexp(:,:,:) + dt * (arkexp(istep,i)*Fexp(:,:,:,i) &
        &  + arkimp(istep,i)*Fimp(:,:,:,i))
    end do

    return
  end subroutine rk_imex_middle

  !============================================================================

  !============================================================================
  ! rk_imex_final - Compute the solution at the next time level using the RK
  !                 IMEX weights b and \hat{b}.
  !
  ! Input parameters:
  ! n_stage - number of stage.
  ! dt - current time step.
  
  subroutine rk_imex_final(n_stages,dt)

    ! Load modules
    use variables
    use time_integ_coeff 

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_stages
    real(wp), intent(in) :: dt
    integer :: i

    ! Initialize ug and uhat
    ug = uold
    uhat = uold

    ! Update ug and uhat
    do i = 1,n_stages
      ug(:,:,:) = ug(:,:,:) + dt * brk (i) *(Fimp(:,:,:,i) + Fexp(:,:,:,i))
      uhat(:,:,:) = uhat(:,:,:) + dt * brkh(i) *(Fimp(:,:,:,i) + Fexp(:,:,:,i))
    end do

    return
  end subroutine rk_imex_final
  
  !============================================================================
  ! LSRK - Integrates the system of ODEs in time using an explicit fourth-order 
  !        accurate Runge-Kutta scheme.
  !============================================================================

  subroutine LSRK(nstep)

    ! Load modules
    use SSWENO_routines,    only: Negative_Density_Removal
    use time_integ_coeff,   only: alsrk, brk, brkh
    use variables,          only: ug, uhat, du, dudt
    use referencevariables, only: ihelems, nodesperelem
    use controlvariables,   only: timestep
    use initcollocation,    only: element_properties
    
    ! Nothing is implicitly defined
    implicit none

    integer,  intent(in) :: nstep
    
    integer              :: inode, ielem

    
    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

       do inode = 1, nodesperelem
          dU(:,inode,ielem) = alsrk(nstep)*dU(:,inode,ielem) + timestep   * dUdt(:,inode,ielem)
          ug(:,inode,ielem) =              ug(:,inode,ielem) + brk (nstep)* dU  (:,inode,ielem)
        uhat(:,inode,ielem) =            uhat(:,inode,ielem) + brkh(nstep)* dU  (:,inode,ielem)
       enddo

!     call Negative_Density_Removal(nodesperelem,ielem,ug(:,:,ielem))

    end do

  end subroutine LSRK

  !============================================================================
  ! timestepcontroller - Adjusts the current timestep based on the error 
  !                      committed in the two previous timesteps. There is no 
  !                      check for excessive error and no rejection of steps.
  ! 
  ! Input parameters:
  ! time_step_cnt - time-step counter. 

  subroutine timestepcontroller(time_step_cnt)

    ! Load modules
    use controlvariables
    use navierstokes

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp) :: predicted_time_step, tol
    real(wp) :: al, be, ga
    real(wp), save :: err_time_tm0
    real(wp), save :: dt0, dt1, dt2

    integer,  parameter :: p = 3

    if(Dt_by_CFL .eqv. .true.) then
      call compute_explicit_timestep(timestep)
      return
    endif

    ! If fixed time-step is selected then exit from this subroutine
    if(variabletimestep .eqv. .false.) then
      return
    endif

    ! tol = err_space_lf / 100000
    tol = max(err_space_lf / 100000.0_wp, 1.0e-12_wp)
    err_time_tm0 = err_time_lf

    !     if(er1 <= 1.0e-12) er1 = er0
    !     if(er2 <= 1.0e-12) er2 = er0
    if(new .and. time_step_cnt == 1) then
      err_time_tm1 = err_time_tm0
      err_time_tm2 = err_time_tm0
    endif
    al  = 0.49_wp/p
    be  = 0.34_wp/p
    ga  = 0.10_wp/p

    ! Estimate the new time-step
    predicted_time_step = (0.9_wp*timestep) * ((tol/err_time_tm0)**al) * &
      & ((err_time_tm1/tol)**be) * ((tol/err_time_tm2)**ga)
    
    ! Select the appropriate time-step
    if(predicted_time_step < 0.9_wp*timestep) then
      timestep = predicted_time_step
    else
      timestep = min(predicted_time_step,1.1*timestep)
    endif

    ! Shuffle time data to get ready for next time step
    err_time_tm2 = err_time_tm1
    err_time_tm1 = err_time_tm0
    dt2 = dt1
    dt1 = dt0

    return
  end subroutine timestepcontroller

  !============================================================================

  !============================================================================
  ! rk_imex - Advances with an RK IMEX scheme the nonlinear system of ODEs 
  !           arising from the spatial discretization of the PDEs.
  
  subroutine rk_imex()
    
    ! Load modules
    use controlvariables
    use referencevariables
    use errorestimation
    use initialize_CSR
    use navierstokes
    use time_integ_coeff 
    use mpimod
    use petsc_snes_solver,  only : setup_petsc_solver, solve_implicit_petsc

    use variables, only: ug, uold, uexp, Fexp, Fimp

    ! Nothing is implicitly defined
    implicit none

    logical  :: converged = .false.
    integer  :: irkstep
    integer :: i_err

    ! Set up storage for implicit method
    call csr_term_footprints()
!   write(*,*) 'csr_term_footprints'
    call csr_get_pointers()
!   write(*,*) 'csr_combine_pointers'
    call csr_initialize_jacobian()
!   write(*,*) 'csr_get_pointers'
    call csr_combine_pointers()
!   write(*,*) 'csr_initialize_jacobian'
    call csr_combine_pointers_element()
!   write(*,*) 'csr_combine_pointers_element'
    call csr_on_element_operator() 
!   write(*,*) 'csr_on_element_operator'
!   call csr_on_element_Matrix_Multiply()
!   write(*,*) 'csr_on_element_Matrix_Multiply'

    ! Test the CSR format
    !call csr_testing()

    ! Set-up PETSc solver
    call setup_petsc_solver()

    ! Initialize to zero the number of time step execuyted
    itimestep = 0

    ! Both the implicit and explicit portions of the IMEX scheme are exercised 
    if (IMEX_element == 'implicit' .or. IMEX_penalty == 'implicit') then

      ! Start the time loop
      rk_imex_time_loop: do

        ! Check to see if the global time has hit the maximum value defined in the
        ! startip file
        if (timeglobal >= timemaximum) exit
        
        ! Update the time-step counter
        itimestep = itimestep + 1

        ! Limit the time-step such that we don't exceed the maximum value defined 
        ! in the startip file
        timestep = min(timestep,timemaximum-timeglobal)
        
        ! The first stage is explicit and different from the other stages
        irkstep = 1
        timelocal = timeglobal

        ! Update uold
        uold = ug

        ! Initialize to zero both explicit and implicit residuals
        Fexp = 0.0_wp
        Fimp = 0.0_wp
        
        ! Calculate the implicit residual contribution for the first stage
        call nse_calcrhsimplicit(irkstep,timelocal) 

        ! Calculate the explicit residual contribution for the first stage
        call nse_calcrhsexplicit(irkstep,timelocal) 

        ! Set the flag to calculate Jacobian matrix for preconditioning in the NK 
        ! solver
        calcjacobian = .false.

        ! Start implicit stages
        rk_imex_stage_loop: do irkstep = 2, narksteps

          ! The variables current_stage_imex is defined in time_integ_coeff.f90 
          ! and used by the nonlinear solver to compute the nonlinear residual	
          current_stage_imex = irkstep

          ! Update flow time of current time step
          timelocal = timeglobal + crk(irkstep)*timestep

          ! Calculate residual contribution from previous stages
          call rk_imex_middle(irkstep,timestep) 

          ! Solve the nonlinear system of ODEs
          call solve_implicit_petsc(converged) 

          ! Fault tolerance, recompute time-step
          if (.not.converged) then
            itimestep = itimestep - 1

            ! Half time-step
            timestep = 0.5_wp*timestep

            ! Reset the solution to what it was at the beginning of the time step
            ug = uold
            if(myprocid==0) write(*,*) 'Time step diverged: recomputing with  &
              &half time-step, time-step = ',timestep

            ! exit the imex loop because a new time step is used
            cycle rk_imex_time_loop

          end if

          ! Calculate the implicit residual contribution for the current stage
          call nse_calcrhsimplicit(irkstep,timelocal)

          ! Calculate the explicit residual contribution for the current stage
          call nse_calcrhsexplicit(irkstep,timelocal)

          ! Flush open unit
          call flush()

        end do rk_imex_stage_loop

        ! Compute the fully integrated solution for the current time step
        call rk_imex_final(narksteps,timestep)

        ! Update the time-step and stage time
        timeglobal = timeglobal + timestep

        ! Calculate the embedded errors
        call calcembeddedtemporalerror() 
        call calcembeddedspatialerror()  
        
        ! Evaluate next dt
        call timestepcontroller(itimestep)

        ! Compute time-averaged quantities and write output if needed
        call post_process_ts(itimestep,timeglobal)
        
        ! Output that this timestep is complete
        if (verbose .and. myprocid == 0)then
          if(mod(itimestep-1,20) == 0 ) then  
            write(*,100)
          endif
          write(*,101) itimestep, timeglobal, timestep, err_time_lf, err_space_lf
        endif

      end do rk_imex_time_loop

    ! Only the explicit portion of the IMEX scheme is exercised
    else if (IMEX_element == 'explicit' .and. IMEX_penalty == 'explicit') then

      ! Start the time loop
      exp_rk_imex_time_loop: do

        ! Check to see if the global time has hit the maximum value defined in the
        ! startip file
        if (timeglobal >= timemaximum) exit
        
        ! Update the time-step counter
        itimestep = itimestep + 1

        ! Limit the time-step such that we don't exceed the maximum value defined 
        ! in the startip file
        timestep = min(timestep,timemaximum-timeglobal)
        
        ! The first stage is explicit and different from the other stages
        irkstep = 1
        timelocal = timeglobal

        ! Update uold
        uold = ug

        ! Initialize to zero both explicit and implicit residuals
        Fexp = 0.0_wp
        Fimp = 0.0_wp

        ! Calculate the explicit residual contribution for the first stage
        call nse_calcrhsexplicit(irkstep,timelocal) 

        ! Start implicit stages
        exp_rk_imex_stage_loop: do irkstep = 2, narksteps

          ! The variables current_stage_imex is defined in time_integ_coeff.f90 
          ! and used by the nonlinear solver to compute the nonlinear residual	
          current_stage_imex = irkstep

          ! Update flow time of current time step
          timelocal = timeglobal + crk(irkstep)*timestep

          ! Calculate residual contribution from previous stages
          call rk_imex_middle(irkstep,timestep)
          
          ! Set ug to the new stage value
          ug = uexp

          ! Calculate the explicit residual contribution for the current stage
          call nse_calcrhsexplicit(irkstep,timelocal)

          ! Flush open unit
          call flush()

        end do exp_rk_imex_stage_loop

        ! Compute the fully integrated solution for the current time step
        call rk_imex_final(narksteps,timestep)

        ! Update the time-step and stage time
        timeglobal = timeglobal + timestep

        ! Calculate the embedded errors
        call calcembeddedtemporalerror() 
        call calcembeddedspatialerror()  
        
        ! Evaluate next dt
        call timestepcontroller(itimestep)

        ! Compute time-averaged quantities and write output if needed
        call post_process_ts(itimestep,timeglobal)
        
        ! Output that this timestep is complete
        if (verbose .and. myprocid == 0)then
          if(mod(itimestep-1,20) == 0 ) then  
            write(*,100)
          endif
          write(*,101) itimestep, timeglobal, timestep, err_time_lf, err_space_lf
        endif

      end do exp_rk_imex_time_loop

    endif

    100 format('     step      |&
               &         timeG         |&
               &           dt          |&
               &       Linf_time       |&
               &      Linf_space ')
    101 format(i8,3x, 4(ES12.5,1x))

    ! Wait for other processes
    call mpi_barrier(PETSC_COMM_WORLD,i_err)

    return
  end subroutine rk_imex

  !============================================================================

  !============================================================================
  ! post_process_ts - Computes time-averaged quantities and write output if 
  !                   needed.
  !
  ! Input parameters:
  ! time_step_cnt - time-step counter.

  subroutine post_process_ts(time_step_cnt,global_time)
    
    ! Load modules
    use controlvariables
    use referencevariables
    use restart_simulation
    use write_solution_file
    use navierstokes
    use time_average
    use aerodynamic_coefficients
    use error_bc_no_slip_wall
    use dkinetic_energy_dt_enstrophy
    use error_heat_entropy_flow_wall_bc
    use nsereferencevariables, only : viscous, InitialCondition

      
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time
    logical :: done_vorticity

    ! Set done_vorticity
    done_vorticity = .false.

    ! Compute time-averaged quantities if needed
    if (time_averaging) then
      ! Compute time average quantities
      call compute_time_averaged_quantities(time_step_cnt)
    endif

    ! Write solution file if needed
    if (write_solution .eqv. .true.) then
      if (MOD(time_step_cnt,write_solution_frequency)==0 .or. &
        & global_time >= timemaximum) then

        ! Compute gradient of the entropy variables without LFG penalty
        ! This is used to compute the vorticity below
        if (viscous .eqv. .false.) then
          call compute_gradient_entropy_variables()
        end if

        ! Compute vorticity field
        call compute_vorticity_field_elements()

        ! Set done_vorticity
        done_vorticity = .true.

        ! Compute specific entropy
        call compute_specific_entropy_elements()

!-- DAVID DEBUG START
        ! Write solution file 
!        call write_solution_vtu_file()
        write(*,*)' timeinteg: post_process_ts you have turned off write_solution_vtu_file'
!-- DAVID DEBUG END
      
        ! Master node writes the .pvtu file
        if (myprocid .eq. 0) then
          call write_solution_pvtu_file()
        endif

      endif
      if (abs(global_time - timemaximum) <= 1.0e-18_wp) then
        if(InitialCondition  == 'SodsProblem') then
          call Sods_Line_Plot()
        endif
!       if(InitialCondition  == 'ExactSolutionViscousShock') then
!         call Sods_Line_Plot()
!       endif
      endif
    endif

    ! Write restart file for each processor if needed
    if (write_restart .eqv. .true.) then
      if (MOD(time_step_cnt,write_restart_frequency)==0 .or. &
        & global_time >= timemaximum) then
      
      ! Write restart file
      call write_restart_file()        
      
      endif
    endif

    ! Write time derivative of the kinetic energy and the enstrophy
    if (write_dke_dt .eqv. .true.) then
      ! Compute time derivative of the kinetic energy
!     This is a bug.  Calling here the ug and dudt are out of sync 
!                     (i.e., ug has been updated while dudt is still at the final stage)
!     call compute_dkinetic_energy_dt()

      if (done_vorticity) then
        ! Compute enstrophy
        call compute_enstrophy()
      else

        ! Compute gradient of the entropy variables without LFG penalty
        ! This is used to compute the vorticity below
        if (viscous .eqv. .false.) then
          call compute_gradient_entropy_variables()
        end if

        ! Compute vorticity field
        call compute_vorticity_field_elements()
        
        ! Compute enstrophy
        call compute_enstrophy()
      end if

      if (myprocid .eq. 0) then
        call write_dkinetic_energy_dt(time_step_cnt,global_time)
        call write_enstrophy(time_step_cnt,global_time)
      end if
    end if

    ! Write velocity and heat entropy flow error at the wall if required
    if (write_errors_wall .eqv. .true.) then
      ! Compute no-slip wall BC error
      call compute_bc_no_slip_wall_error()

      ! Compute heat entropy flow error
      call compute_heat_entropy_flow_wall_bc_error()

      if (myprocid .eq. 0) then
        call write_error_no_slip_wall_bc(time_step_cnt,global_time)
        call write_error_heat_entropy_flow_wall_bc(time_step_cnt,global_time)
      end if
    end if

    ! Write aerodynamic coefficients and error at the wall for the no-slip BC
    if (write_aero_coeffs .eqv. .true.) then
      ! Compute aerodynamic coefficients
      call compute_aerodynamic_coefficients()

      if (myprocid .eq. 0) then
        call write_aerodynamic_coefficients(time_step_cnt,global_time)
      end if
    end if

    
    return
  end subroutine post_process_ts

  !============================================================================
  
  !============================================================================
  ! post_process_ts_0 - Computes quantities and write output if at t = 0, i.e.,
  ! at the initial state

  subroutine post_process_ts_0(time_step_cnt,global_time)
    
    ! Load modules
    use controlvariables
    use referencevariables
    use restart_simulation
    use write_solution_file
    use navierstokes
    use time_average
    use aerodynamic_coefficients
    use error_bc_no_slip_wall
    use dkinetic_energy_dt_enstrophy
    use error_heat_entropy_flow_wall_bc
      
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in)  :: time_step_cnt
    real(wp), intent(in) :: global_time
    logical :: done_vorticity

    ! Set done_vorticity
    done_vorticity = .false.


    ! Write solution file if needed
    if (write_solution .eqv. .true.) then

      ! Compute vorticity field
      call compute_vorticity_field_elements()

      done_vorticity = .true.

      ! Compute specific entropy
      call compute_specific_entropy_elements()

      ! Write solution file
      call write_solution_vtu_file()
      
      ! Master node writes the .pvtu file
      if (myprocid .eq. 0) then
        call write_solution_pvtu_file()
      end if

    end if

    ! Write time derivative of the kinetic energy
    if (write_dke_dt .eqv. .true.) then
      ! Compute time derivative of the kinetic energy
      call compute_dkinetic_energy_dt()

      if (done_vorticity) then
        ! Compute enstrophy
        call compute_enstrophy()
      else
        ! Compute vorticity field
        call compute_vorticity_field_elements()

        ! Compute enstrophy
        call compute_enstrophy()
      end if

      if (myprocid .eq. 0) then
        call write_dkinetic_energy_dt(time_step_cnt,global_time)
        call write_enstrophy(time_step_cnt,global_time)
      end if
    end if

    ! Write velocity and heat entropy flow error at the wall if required
    if (write_errors_wall .eqv. .true.) then
      ! Compute no-slip wall BC error
      call compute_bc_no_slip_wall_error()

      ! Compute heat entropy flow error
      call compute_heat_entropy_flow_wall_bc_error()

      if (myprocid .eq. 0) then
        call write_error_no_slip_wall_bc(time_step_cnt,global_time)
        call write_error_heat_entropy_flow_wall_bc(time_step_cnt,global_time)
      end if
    end if

    ! Write aerodynamic coefficients and error at the wall for the no-slip BC
    if (write_aero_coeffs .eqv. .true.) then
      ! Compute aerodynamic coefficients
      call compute_aerodynamic_coefficients()

      if (myprocid .eq. 0) then
        call write_aerodynamic_coefficients(time_step_cnt,global_time)
      end if
    end if
    
    return
  end subroutine post_process_ts_0

  !============================================================================
  
  !============================================================================

end module timeinteg
