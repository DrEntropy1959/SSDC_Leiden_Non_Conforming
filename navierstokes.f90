module navierstokes
  ! This module contains the necessary routines to simulate the
  ! Navier-Stokes equations using the SSDC method. 
  use precision_vars
  implicit none

  ! This procedure pointer is called for the initial condition and
  ! for imposing boundary data and has the same interface as
  ! isentropicVortexFull. It may point to any subroutine with
  ! the same interface.
  procedure(isentropicVortexFull), pointer :: InitialSubroutine => null()
  procedure(isentropicVortexFull), pointer :: BoundaryCondition => null()

  private

  public dUdV, dVdU, dVdW, dWdV, dWdU, dUdW
  public conserved_to_primitive
  public primitive_to_conserved
  public primitive_to_entropy
  public entropy_to_primitive
  public EntropyConsistentFlux
  public CharacteristicDecomp
  public normalflux
  public matrix_hatc_node
  public normalviscousflux
  public roeavg
  public viscousflux3D
  public viscousflux
  public isentropicVortexFull
  public supersonicVortexFull
! public viscous_dissipation_2pt


  public nse_calcinitialcondition
  public nse_initializesemidiscretization
  public nse_calc_dudt_LSRK
  public nse_calcembeddederror
  public nse_calcerror
  public nse_calcrhsexplicit
  public nse_calcrhsimplicit
  public nse_calcembeddedspatialerror
  public nse_communicationsetup
  public Navier_Stokes_init
  public UniformFreeStream

  public Flux_Divergence
  public Solution_Filter
  public compute_vorticity_field_elements
  public compute_specific_entropy_elements
  public nse_reconcilestates
  public rhalf
  public compute_explicit_timestep
  public kinetic_energy_element
  public compute_gradient_entropy_variables
  public Calc_Entropy_Viscosity


contains

  subroutine Navier_Stokes_init()
    use controlvariables
    use initialize_CSR
    use collocationvariables
    use nsereferencevariables
    use initgrid

    ! calculate initial condition
    timeglobal = 0.0_wp

    ! If the flow is inviscid then set the penalty interfaces parameters l01, 
    ! l10 and l00 to zero
    if (.not. viscous) then
      l01 = 0.0_wp
      l10 = 0.0_wp
      l00 = 0.0_wp
    end if

    call nse_calcinitialcondition() ! (navierstokes)
    ! allocate memory for semidiscretization
    call nse_initializesemidiscretization() !(navierstokes)
    ! set up parallel communication
    call nse_communicationsetup()

    ! set up extrapolation points from adjoining elements to bridge the interfaces
    if(discretization == 'SSWENO') call WENO_Intrp_Face_Nodes()

    return
  end subroutine Navier_Stokes_init

  subroutine nse_calcinitialcondition()
    ! This routine sets the nondimensional parameters for the problem;
    ! sets the initial condition routine and calculate the initial data;
    ! and allocates the memory for arrays required to describe the solution.
    ! The grid and connectivities must already be populated.

    ! Load modules
    use variables, only: vg, wg, ug, xg, omega, &
                         specific_entropy, mean_vg, time_ave_prod_vel_comp, &
                         reynolds_stress, mut, Log_Ave_Counter
    use initcollocation, only: element_properties
    use referencevariables
    use nsereferencevariables
    use controlvariables
    use restart_simulation
    use mpimod
    
    ! Nothing is implicitly defined
    implicit none

    !
    ! Local Variables
    ! ===============
    ! Indices
    integer :: inode,ielem, nodesperelem_max
    ! Output unit
    integer, allocatable, dimension(:) :: iunit
    ! Normals are needed in the call to InitialSubroutine 
    ! but since we don't care about fv in this routine
    ! we can set them to arbitrary values
    real(wp) :: ctmp(3)
    ! We don't care about viscous flux, but need a vector for it.
    real(wp), allocatable :: fvtmp(:), phitmp(:,:)

!   real(wp), dimension(:,:), pointer :: my_pointer => null()

    integer :: i_err

    nodesperelem_max = (npoly_max+1)**(ndim)

    ! Arbitrary value for normals
    ctmp = 0.0_wp

    ! We are limited to the calorically perfect Navier-Stokes equations
    nequations = 5

    ! Set the flow parameters
    call set_Flow_parameters()

    ! Set initial and boundary condition procedure
    call set_InitialCondition(InitialCondition)

    ! Allocate memory for flow state data
    !
    ! Conserved variables (denisty, momentum, total energy)
    allocate(ug(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    ug = 0.0_wp
    
    ! Primitive variables (density, velocity, temperature)
    allocate(vg(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    vg = 0.0_wp
    
    ! Entropy variables (see primitive_to_entropy subroutine)
    allocate(wg(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    wg = 0.0_wp
    
    ! Vorticity field (\nabla \times velocity)
    allocate(omega(1:3,1:nodesperelem_max,ihelems(1):ihelems(2)))
    omega = 0.0_wp

    ! Entropy field
    allocate(specific_entropy(1:nodesperelem_max,ihelems(1):ihelems(2)))
    specific_entropy = 0.0_wp

    allocate(mut(1:nodesperelem_max,ihelems(1):ihelems(2)))
    mut = 0.0_wp

    ! Allocate memory if time averaging is required
    if (time_averaging) then

      ! Time-average of primitive variables
      allocate(mean_vg(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      mean_vg = 0.0_wp

      ! Time-average of the product of the velocity components
      allocate(time_ave_prod_vel_comp(6,1:nodesperelem_max,ihelems(1):ihelems(2)))
      time_ave_prod_vel_comp = 0.0_wp

      ! Reynolds stresses
      allocate(reynolds_stress(6,1:nodesperelem_max,ihelems(1):ihelems(2)))
      reynolds_stress = 0.0_wp
    
    endif

    ! We have a different output unit for each equation.
    ! if this were a real code we would replace this with something more useful.
    allocate(iunit(nequations))

    ! For our call to InitialSubroutine we need the viscous fluxes
    allocate(fvtmp(nequations))
    allocate(phitmp(nequations,ndim))

    Log_Ave_Counter = 0

    if (new .eqv. .true.) then

      ! Loop over elements
      do ielem = ihelems(1),ihelems(2)

        call element_properties(ielem, n_pts_3d=nodesperelem)

        ! Loop over nodes in each element
        do inode = 1, nodesperelem
          ! Use exact solution routine to initialize data
          call InitialSubroutine( &
            Vx = vg(:,inode,ielem), &
            phi = phitmp, &
            fv = fvtmp, &
            Jx = ctmp, &
            xin = xg(:,inode,ielem), &
            tin = 0.0_wp, &
            neqin = nequations, &
            nd = ndim, &
            mut = mut(inode,ielem)) ! (navierstokes)
          ! Calculate conservative variables from primitive variables
          call primitive_to_conserved( &
            vin = vg(:,inode,ielem), &
            uout = ug(:,inode,ielem), &
            nq = nequations ) ! (navierstokes)
          ! Calculate entropy variables from primitive variables
          call primitive_to_entropy( &
            vin = vg(:,inode,ielem), &
            wout = wg(:,inode,ielem), &
            nq = nequations ) ! (navierstokes)
        end do
      end do

    else

      ! Read solution from restart file
      call read_restart_file()
 
      ! Loop over elements
      do ielem = ihelems(1),ihelems(2)

        call element_properties(ielem, n_pts_3d=nodesperelem)

        ! Loop over nodes in each element
        do inode = 1, nodesperelem
          ! Calculate primitive variables from conservative variables
          call conserved_to_primitive( &
            uin = ug(:,inode,ielem), &
            vout = vg(:,inode,ielem), &
            nq = nequations)
          ! Calculate entropy variables from primitive variables
          call primitive_to_entropy( &
            vin = vg(:,inode,ielem), &
            wout = wg(:,inode,ielem), &
            nq = nequations )
        enddo
      enddo
    endif

    ! Deallocate temporary viscous flux array
    deallocate(fvtmp,phitmp)
    deallocate(iunit)

    ! Wait for other processes
    call mpi_barrier(PETSC_COMM_WORLD,i_err)

    return
  end subroutine nse_calcinitialcondition

  !============================================================================
  
  !============================================================================
  ! nse_communicationsetup - Allocates the ghost data for both explicit and
  ! implicit Navier-Stokes computations. The implicit part sets the ghost data
  ! for the viscous penalty calculation which requires the knowledge of the
  ! internal solution in the adjoining element.
  
  subroutine nse_communicationsetup()
    
    ! Load modules
    use referencevariables
    use initgrid
    use variables
!   use variables, only: ug, ughst, ughstWENO, uelemghst, phig, phighst, ef2e, & 
!     & kfacenodes, r_x, r_x_ghst, xgWENO_partner,xghstWENO_partner
    use controlvariables, only: IMEX_penalty, discretization
    use petscvariables
!   use petscvariables, only: upetsc,     ulocpetsc,     &
!                             upetscWENO, ulocpetscWENO, &
!                             xpetscWENO_partner, xlocpetscWENO_partner, &
!                             phipetsc,   philocpetsc,   &
!                             uelempetsc, uelemlocpetsc, &
!                             r_x_petsc,  r_x_loc_petsc
    use mpimod, only: PetscComm0DDataSetup, PetscComm1DDataSetup, PetscComm1DDataSetupWENO,    &
      & PetscComm1D_LGL_Shell_DataSetup, &
      & PetscComm1DElementDataSetup, PetscComm2DDataSetup, PetscComm2DGeomDataSetup, &
      & PetscComm1DDataSetupWENOGeom, PetscCommShellDataSetup

    use nsereferencevariables, only : viscous
    
    ! Nothing is implicitly defined
    implicit none

    integer :: nshell

    ! Setup the PETSc parallel communication for exchanging the conservative variables at the parallel face element 
    call PetscComm1DDataSetup(ug,ughst,upetsc,ulocpetsc, size(ug,1),size(ug,2), nelems, size(ughst,2))

    ! Setup the PETSc parallel communication for exchanging the conservative variables at the parallel face element 
    call PetscComm0DDataSetup(mut,mutghst,mutpetsc,mutlocpetsc, size(mut,1), nelems, size(mutghst))

    ! PETSc parallel communication setup for exchanging the conservative variables at the parallel face element 
    if(discretization == 'SSWENO')then

      nshell = nfacesperelem*nodesperface

      call PetscComm1DDataSetupWENO(ug,ughstWENO,upetscWENO,ulocpetscWENO,nequations, &
      & nodesperelem,nelems,nghost)

!     call PetscComm1D_LGL_Shell_DataSetup(ugWENO_partner,ughstWENO_partner,upetscWENO_Shell, &
!     & ulocpetscWENO_Shell,nequations,nshell,nelems,nghost)

      call PetscCommShellDataSetup(ugWENO_partner,ughstWENO_partner,upetscWENO_Shell, &
      & ulocpetscWENO_Shell,nequations,nshell,nelems,nghost)

      allocate(xgWENO_partner(3,nshell,ihelems(1):ihelems(2))) ; xgWENO_partner = -10000.0_wp ;
      allocate(xghstWENO_partner(3,nghost)) ; xghstWENO_partner = 0.0_wp

      call PetscComm1DDataSetupWENOGeom(xgWENO_partner,xghstWENO_partner,xpetscWENO_partner, &
      & xlocpetscWENO_partner,3,nshell,nelems,nghost)

    endif

    ! Setup the PETSc parallel communication for exchanging the gradient of the 
    ! entropy variables at the parallel face element 
    call PetscComm2DDataSetup(phig,phighst,phipetsc,philocpetsc,                    &
                     size(phig,1),size(phig,2),size(phig,3), nelems, size(phighst,3))

    ! Setup the PETSc parallel communication for exchanging the conservative 
    ! variables in the adjoining parallel element and the geometrical data r_x
    ! at the parallel interfaces
    if (IMEX_penalty == 'implicit')  then
      if (viscous) then
        ! Conservative variables in the ghost element
        call PetscComm1DElementDataSetup(ug,uelemghst,uelempetsc, &
          & uelemlocpetsc,nequations,nodesperelem,nelems,nghost_elem)

        ! Geometrical data at the parallel interfaces
        call PetscComm2DGeomDataSetup(r_x,r_x_ghst,r_x_petsc,r_x_loc_petsc, &
          & 3,nodesperelem,nelems,nghost)

      endif
    endif

  end subroutine nse_communicationsetup

  !============================================================================

  pure subroutine primitive_to_conserved(vin,uout,nq)
    ! this routine calculates the conserved variables
    ! (density, momentum, total energy)
    ! from the primitive variables (density, velocity, temperature).
    ! Velocity is nondimensionalized by reference velocity, density
    ! is nondimensionalized by reference density, and temperature
    ! is nondimensionalized by reference temperature. 
    ! the specific heat is nondimensionalized by reference specific
    ! heat and the specific gas constant is nondimensionalized by
    ! the reference specific gas constant.

    use nsereferencevariables, only: gm1M2, gamI
    implicit none

    integer,  intent(in ) :: nq                                                        ! number of equations
    real(wp), intent(in ) :: vin(nq)                                                   ! primitive variables
    real(wp), intent(out) :: uout(nq)                                                  ! conserved variables

    uout(1)   = vin(1)                                                                 ! density

    uout(2:4) = vin(1)*vin(2:4)                                                        ! momentum

    uout(5)   = vin(1)*( gamI * vin(5) + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) ) ! energy

  end subroutine primitive_to_conserved

  !============================================================================

  pure subroutine conserved_to_primitive(uin,vout,nq)
    ! calculate the primitive variables from the conserved variables as subroutine call

    use nsereferencevariables, only: gm1M2, gamma0
    implicit none

    integer,  intent(in ) :: nq                                                            ! number of equations
    real(wp), intent(in ) :: uin(nq)                                                       ! conserved variables
    real(wp), intent(out) :: vout(nq)                                                      ! primitive variables

    vout(1)   = uin(1)                                                                     ! density

    vout(2:4) = uin(2:4)/uin(1)                                                            ! velocity

    vout(5)   = (uin(5)/uin(1) - gm1M2*0.5_wp*dot_product(vout(2:4),vout(2:4)) ) * gamma0  ! temperature

  end subroutine conserved_to_primitive

  !============================================================================

  pure function conserved_to_primitive_F(uin,nq)
    ! calculate the primitive variables from the conserved variables as function call

    use nsereferencevariables, only: gm1M2, gamma0
    implicit none

    integer,  intent(in) :: nq                                                             ! number of equations
    real(wp), intent(in) :: uin(nq)                                                        ! primitive variables

    real(wp) :: conserved_to_primitive_F(nq)                                               ! output primitive variables

    conserved_to_primitive_F(1) = uin(1)                                                   ! density

    conserved_to_primitive_F(2:4) = uin(2:4)/uin(1)                                        ! velocity

    conserved_to_primitive_F(5) = ( uin(5)/uin(1) &
      - gm1M2*0.5_wp*dot_product(conserved_to_primitive_F(2:4),conserved_to_primitive_F(2:4)) ) * gamma0 ! temperature

  end function conserved_to_primitive_F

  !============================================================================

  pure function conserved_to_entropy(uin,nq)

    use nsereferencevariables, only: gm1M2, gamma0
    implicit none

    integer, intent(in) :: nq                                                             ! number of equations
    real(wp), intent(in) :: uin(nq)                                                       ! primitive variables

    real(wp)             :: vin(nq)
    real(wp)             :: conserved_to_entropy(nq)
    real(wp)             :: Tinv

    vin(1) = uin(1)                                                                       ! density

    vin(2:4) = uin(2:4)/uin(1)                                                            ! velocity

    vin(5) = ( uin(5)/uin(1) - gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) ) * gamma0     ! temperature

    Tinv = 1.0_wp/vin(5)

    conserved_to_entropy(1) = 1.0_wp-0.5_wp*gm1M2*dot_product(vin(2:4),vin(2:4)) * Tinv & ! w_1 = h/T - s - (gamma_0 - 1) M_0^2 u_k u_k/(2T)
      -specificentropy(vin,nq)

    conserved_to_entropy(2:4) = gm1M2*vin(2:4) * Tinv                                     ! w_{k+1} = (gamma_0 - 1) M_0^2 u_k/T, k = 1,2,3

    conserved_to_entropy(5) = - Tinv                                                      ! w_5 = -1/T

  end function conserved_to_entropy

  !============================================================================

  pure function specificentropy(vin,nq)
    ! calculate the specific thermodynamic entropy using primitive variable vector

    use nsereferencevariables, only: gm1og, gamI
    implicit none

    integer,  intent(in) :: nq            ! number of equations
    real(wp), intent(in) :: vin(nq)       ! primitive variables

    real(wp) :: specificentropy           ! output thermodynamic specific entropy

    specificentropy = gamI * log(vin(5)) - gm1og*log(vin(1))

  end function specificentropy

  !============================================================================

  pure subroutine primitive_to_entropy(vin,wout,nq)
    ! this routine calculates the entropy variables corresponding
    ! to the entropy--entropy flux pair (S,F^i) = (-rho*s,-rho*u_i*s)
    ! using the primitive variables.
    use nsereferencevariables, only: gm1M2
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! entropy variables
    real(wp), intent(out) :: wout(nq)

    ! w_1 = h/T - s - (gamma_0 - 1) M_0^2 u_k u_k/(2T)
    wout(1) = 1.0_wp-0.5_wp*gm1M2*dot_product(vin(2:4),vin(2:4))/vin(5) &
      -specificentropy(vin,nq)
    ! w_{k+1} = (gamma_0 - 1) M_0^2 u_k/T, k = 1,2,3
    wout(2:4) = gm1M2*vin(2:4)/vin(5)
    ! w_5 = -1/T
    wout(5) = -1.0_wp/vin(5)

    return
  end subroutine primitive_to_entropy

  !============================================================================

  pure subroutine entropy_to_primitive(win,vout,nq)
    use nsereferencevariables, only:gm1M2, gm1og
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: win(nq)
    ! entropy variables
    real(wp), intent(out) :: vout(nq)

    real(wp)              :: w5inv

    continue

    w5inv = 1.0_wp / win(5)

    ! Temperature
    vout(5) = - w5inv

    ! Velocity components
    vout(2:4) = -1.0_wp/gm1M2*win(2:4)*w5inv

    ! Density
    vout(1) = exp(- 1.0_wp*(win(2)**2 + win(3)**2 + win(4)**2 &
                & - 2.0_wp*gm1M2*(win(1)-1.0_wp)*win(5) &
                & + 2.0_wp*gm1M2*(gm1og-1.0_wp)*win(5)*log(-w5inv)) &
                & /(2.0_wp*gm1M2*gm1og*win(5)))
    
    return
  end subroutine entropy_to_primitive

!===================================================================================================

  pure function normal_Entropy_flux(vin,nx,nq)

    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! normal vector
    real(wp), intent(in) :: nx(3)
    ! primitive variables
    real(wp), intent(in) :: vin(nq)

    ! output normal flux
    real(wp) :: normal_Entropy_flux

    ! calculate normal mass flux
    normal_Entropy_flux = vin(1)*dot_product(vin(2:4),nx) * specificentropy(vin(:),nq)

  end function normal_Entropy_flux

!===================================================================================================

  pure function normalflux(vin,nx,nq)
    ! this function calculates the convective flux in the normal
    ! direction. Note that because we nondimensionalize the pressure
    ! by the reference pressure that there is a nondimensional parameter in the flux.

    use nsereferencevariables, only: gm1M2, gM2
    implicit none

    integer,  intent(in) :: nq                                                  ! number of equations

    real(wp), intent(in) :: nx(3)                                               ! normal vector

    real(wp), intent(in) :: vin(nq)                                             ! primitive variables

    real(wp) :: normalflux(nq)                                                  ! output normal flux

    real(wp) :: un, p                                                           ! local variables:  Normal mass flux and Pressure

    un = vin(1)*dot_product(vin(2:4),nx(:))                                     ! calculate normal mass flux

    p = vin(1)*vin(5)                                                           ! calculate pressure (R = 1)

    normalflux(1) = un                                                          ! mass flux (\rho u \cdot n)

    normalflux(2:4) = un*vin(2:4) + p*nx(:)/gM2                                 ! momentum flux (\rho u \cdot n) u + p n /(gamma_0 M_0^2)

    normalflux(5) = un*( vin(5) + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) ) ! energy flux (\rho u \cdot n) H

  end function normalflux

!===================================================================================================

  pure function normalviscousflux(vin,phi,nx,nq,mut)
    ! this function calculates the viscous flux in 
    ! the normal direction based on the primitive
    ! variables and the gradients of the entropy variables,
    ! penalized with an LDC/LDG methodology. The use
    ! of entropy variables ensures stability.
    use nsereferencevariables, only: Re0inv
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: nq
    ! contravariant vector
    real(wp), intent(in) :: nx(3)
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! entropy variable gradients
    real(wp), intent(in) :: phi(nq,3)

    real(wp), intent(in) :: mut

    ! output normal viscous flux
    real(wp) :: normalviscousflux(nq)

    ! physical space viscous flux
    real(wp) :: fvl(nq,3)
    ! direction index
    integer :: idir

    continue

    fvl(:,:) = ViscousFlux(vin,phi,nq,mut)

    normalviscousflux = 0.0_wp

    do idir = 1,3
      normalviscousflux = normalviscousflux + Re0inv*fvl(:,idir)*nx(idir)
    end do

    return
  end function normalviscousflux

!===================================================================================================

  pure function viscousflux3D(vin,phi,r_x,nq,nd,mut)
    ! this function calculates the viscous flux in 
    ! the three computational space directions based on the primitive
    ! variables and the gradients of the entropy variables,
    ! penalized with an LDC/LDG methodology. The use
    ! of entropy variables ensures stability.
    use nsereferencevariables, only: Re0inv
    
    ! Nothing is implicitly defined
    implicit none
    
    integer,  intent(in) :: nq, nd
    ! contravariant vector
    real(wp), intent(in) :: r_x(3,3)
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! entropy variable gradients
    real(wp), intent(in) :: phi(nq,3)

    real(wp), intent(in) :: mut

    ! output viscous flux in the computational space directions
    real(wp) :: viscousflux3D(nq,nd)

    ! physical space viscous flux
    real(wp) :: fvl(nq,3)
    ! direction index
    integer :: idir, jdir

    continue

    fvl(:,:) = ViscousFlux(vin,phi,nq,mut)

    viscousflux3D = 0.0_wp

    do jdir = 1,nd
      do idir = 1,3
        viscousflux3D(:,jdir) = viscousflux3D(:,jdir) + Re0inv*fvl(:,idir)*r_x(jdir,idir)
      end do
    end do

    return
  end function viscousflux3D

  !============================================================================

  pure function ViscousFlux(vin,phi,nq,mut)

    use nsereferencevariables, only: gm1M2, Pr0, gm1M2I

    use controlvariables, only : variable_viscosity

    implicit none

    integer,  intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! entropy variable gradients
    real(wp), intent(in) :: phi(nq,3)

    real(wp), intent(in) :: mut

    real(wp), dimension(nq,3) :: ViscousFlux

    ! thermal conductivity (normalized by kappa0), Kinetic energy and constants
    real(wp)             :: kappa_Pr0, KE2, t1,t2,t3
    real(wp)             :: u_phi51,v_phi52,w_phi53
    real(wp), parameter  :: third = 1.0_wp/3.0_wp

    ! dynamic viscosity (normalized by mu0)
    real(wp) :: mu

    continue

    ! Set dynamic viscosity
    if (variable_viscosity .eqv. .true.) then
      mu = sutherland_law(vin(5)) + mut
    else
      mu = 1.0_wp + mut
    end if

    ! Initialize viscous flux to zero
    ViscousFlux(:,:) = 0.0_wp

    u_phi51 = vin(2)*phi(5,1)
    v_phi52 = vin(3)*phi(5,2)
    w_phi53 = vin(4)*phi(5,3)

    kappa_Pr0 = 1.0_wp/Pr0
    KE2       = + vin(2)*vin(2) + vin(3)*vin(3) + vin(4)*vin(4)

    t1 = mu*vin(5)
    t2 = t1*third
    t3 = kappa_Pr0*vin(5)*vin(5)

    ! momentum
    ViscousFlux(2,1) = t2*(4.0_wp*(u_phi51+phi(2,1)*gm1M2I) - 2.0_wp*( v_phi52 + w_phi53 + (phi(3,2)+phi(4,3))*gm1M2I ))

    ViscousFlux(3,2) = t2*(4.0_wp*(v_phi52+phi(3,2)*gm1M2I) - 2.0_wp*( u_phi51 + w_phi53 + (phi(2,1)+phi(4,3))*gm1M2I ))

    ViscousFlux(4,3) = t2*(4.0_wp*(w_phi53+phi(4,3)*gm1M2I) - 2.0_wp*( u_phi51 + v_phi52 + (phi(2,1)+phi(3,2))*gm1M2I ))

    ViscousFlux(2,2) = t1*(phi(5,1)*vin(3)+phi(5,2)*vin(2)+(phi(3,1)+phi(2,2))*gm1M2I)
    ViscousFlux(2,3) = t1*(phi(5,1)*vin(4)+phi(5,3)*vin(2)+(phi(4,1)+phi(2,3))*gm1M2I)
    ViscousFlux(3,3) = t1*(phi(5,2)*vin(4)+phi(5,3)*vin(3)+(phi(4,2)+phi(3,3))*gm1M2I)

    ViscousFlux(3,1) = ViscousFlux(2,2)
    ViscousFlux(4,1) = ViscousFlux(2,3)
    ViscousFlux(4,2) = ViscousFlux(3,3)

    ! energy
    ViscousFlux(5,1) = t2*vin(2)*(gm1M2*( v_phi52 + w_phi53 ) &
      +  4.0_wp*phi(2,1)-2.0_wp*(phi(3,2)+phi(4,3)) ) &
      + t1*( vin(3)*(phi(3,1)+phi(2,2))+vin(4)*(phi(4,1)+phi(2,3)) ) &
      + ( t1*gm1M2*(KE2+third*vin(2)*vin(2)) + t3 )*phi(5,1)

    ViscousFlux(5,2) = t2*vin(3)*(gm1M2*( u_phi51 + w_phi53 ) &
      +  4.0_wp*phi(3,2)-2.0_wp*(phi(2,1)+phi(4,3)) ) &
      + t1*( vin(2)*(phi(3,1)+phi(2,2))+vin(4)*(phi(4,2)+phi(3,3)) ) &
      + ( t1*gm1M2*(KE2+third*vin(3)*vin(3)) + t3 )*phi(5,2)

    ViscousFlux(5,3) = t2*vin(4)*(gm1M2*( u_phi51 + v_phi52 ) &
      + 4.0_wp*phi(4,3)-2.0_wp*(phi(3,2)+phi(2,1)) ) &
      + t1*( vin(2)*(phi(4,1)+phi(2,3))+vin(3)*(phi(4,2)+phi(3,3)) ) &
      + ( t1*gm1M2*(KE2+third*vin(4)*vin(4)) + t3 )*phi(5,3)


  end function ViscousFlux

!===================================================================================================

  pure function HoneinMoinFlux(vl,vr,Jx,neqin)
    ! this function calculates the normal Kinetic Energy Preserving
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the nondimensionalization employed
    ! herein. It follows loosely the work of A.E. Honein, P. Moin 
    ! Journal of Computational Physics 201 (2004) 531–545

    use nsereferenceVariables, only: gm1og, gm1M2,gM2I,gamI

    implicit none
    ! Arguments
    ! =========
    ! number of equations
    integer, intent(in) :: neqin
    ! left and right states
    real(wp), intent(in), dimension(neqin) :: vl, vr
    ! metrics scaled by jacobian
    real(wp), intent(in), dimension(3)       :: Jx
  
    ! Function
    ! ========
    real(wp), dimension(neqin) :: HoneinMoinFlux

    ! Local Variables
    ! ===============

    ! Average Velocity
    real(wp), dimension(3)     :: vav

    real(wp) :: unl, unr, pl, pr, P, En, mdot, Tav

    ! normal velocity
    unl = dot_product(Jx,vl(2:4))
    unr = dot_product(Jx,vr(2:4))

    ! mass flux
    mdot = 0.5_wp*(vl(1)*unl + vr(1)*unr)

    ! Velocity average
    vav(:) = 0.5_wp * (vl(2:4) + vr(2:4))

    ! pressure average
    pl = vl(1)*vl(5) ; pr = vr(1)*vr(5)
    P  = 0.5_wp * (pl + pr)

    ! Temperature average
    Tav = 0.5_wp * (vl(5)+vr(5))

    ! Internal energy
    En = Tav * gamI  + 0.5_wp * gm1M2 * dot_product(vav,vav)

    HoneinMoinFlux(1) = mdot 
    HoneinMoinFlux(2) = mdot*vav(1) + Jx(1) * P * gM2I
    HoneinMoinFlux(3) = mdot*vav(2) + Jx(2) * P * gM2I
    HoneinMoinFlux(4) = mdot*vav(3) + Jx(3) * P * gM2I
    HoneinMoinFlux(5) = mdot*En + gm1og*0.5_wp*(pl*unr+pr*unl)

  end function HoneinMoinFlux

!===================================================================================================

! pure function EntropyConsistentFlux(vl,vr,Jx,neq)
       function EntropyConsistentFlux(vl,vr,Jx,neq)

    ! this function calculates the normal entropy consistent
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the nondimensionalization employed
    ! herein and follows directly from the work of Ismail and Roe,
    ! DOI: 10.1016/j.jcp.2009.04.021

    use nsereferencevariables, only: gM2I, gm1og, gp1og, gm1M2

    implicit none

    ! Arguments
    ! =========
    integer,  intent(in)                   :: neq      ! number of equations

    real(wp), intent(in), dimension(neq)   :: vl, vr   ! left and right states

    real(wp), intent(in), dimension(3)     :: Jx       ! metrics scaled by jacobian
  
    ! Function
    ! ========
    real(wp), dimension(neq) :: EntropyConsistentFlux

    ! Local Variables
    ! ===============

    real(wp), dimension(neq) :: vhat           ! temporary variables

    real(wp) :: root_Tl  , root_Tr             !   sqrt[T_i], i=L,R
    real(wp) :: root_Tl_I, root_Tr_I           ! 1/sqrt[T_i], i=L,R

    real(wp) :: tinvav, tinvavinv              ! average of inverse temperature and its inverse

    real(wp) :: s1, s2                         ! Logarithmic averages of density and temperature

    real(wp) :: mdot, P, T                     ! normal mass flux (mdot), Pressure,  Temperature

    continue

    root_Tl   = sqrt(vl(5))    ; root_Tr   = sqrt(vr(5))            ! Sqrt[T_i] 
    root_Tl_I = 1.0_wp/root_Tl ; root_Tr_I = 1.0_wp/root_Tr         ! Sqrt[T_i]^{-1}
    tinvav    = root_Tl_I   + root_Tr_I
    tinvavinv = 1.0_wp/tinvav
  
    vhat(2:4) = (vl(2:4)*root_Tl_I + vr(2:4)*root_Tr_I)*tinvavinv   ! velocity

    P = (vl(1)*root_Tl + vr(1)*root_Tr) * tinvavinv                 ! pressure

    s1 = Logarithmic_Average(root_Tl*vl(1),root_Tr*vr(1))           ! Logarithmic average of rho/Sqrt(T)

    s2 = Logarithmic_Average(root_Tl_I    ,root_Tr_I    )           ! Logarithmic average of   1/Sqrt(T)

    vhat(1) = 0.5_wp *tinvav * s1                                   ! density

    T = 0.5_wp * (gm1og * P + gp1og * s1/s2) / vhat(1)              ! temperature

    vhat(5) = T + 0.5_wp * gm1M2 * dot_product(vhat(2:4),vhat(2:4)) ! total enthalpy

    mdot = vhat(1)*dot_product(vhat(2:4),Jx)                        ! normal mass flow rate

    EntropyConsistentFlux(1) = mdot
    EntropyConsistentFlux(2) = mdot*vhat(2) + Jx(1) * P * gM2I
    EntropyConsistentFlux(3) = mdot*vhat(3) + Jx(2) * P * gM2I
    EntropyConsistentFlux(4) = mdot*vhat(4) + Jx(3) * P * gM2I
    EntropyConsistentFlux(5) = mdot*vhat(5)

  end function EntropyConsistentFlux

!===================================================================================================

       function Logarithmic_Average(a,b)
! pure function Logarithmic_Average(a,b)

!   use variables, only: Log_Ave_Counter

    implicit none

    real(wp),  intent(in) :: a,b

!  Logarithmic expansion valid to us = 0.00000894, or ratio =  1.006
!   real(wp), dimension(0:1), parameter :: c = (/1.0_wp,3.0_wp/5.0_wp/)
!   real(wp), dimension(0:1), parameter :: d = (/2.0_wp,8.0_wp/15.0_wp/)
!   real(wp),                 parameter :: eps = 8.9e-06_wp

!  Logarithmic expansion valid to us = 0.00276, or ratio =  1.111
!   real(wp), dimension(0:2), parameter :: c = (/1.0_wp,10.0_wp/9.0_wp,5.0_wp/21.0_wp /)
!   real(wp), dimension(0:2), parameter :: d = (/2.0_wp,14.0_wp/9.0_wp,128.0_wp/945.0_wp /)
!   real(wp),                 parameter :: eps = 2.7e-03_wp

!  Logarithmic expansion valid to us = 0.0226, or ratio =  1.353
    real(wp), dimension(0:3), parameter :: c = (/1.0_wp,21.0_wp/13.0_wp,105.0_wp/143.0_wp,35.0_wp/429.0_wp /)
    real(wp), dimension(0:3), parameter :: d = (/2.0_wp,100.0_wp/39.0_wp,566.0_wp/715.0_wp,512.0_wp/15015.0_wp /)
    real(wp),                 parameter :: eps = 2.0e-02_wp

!  Logarithmic expansion valid to us = 0.0677, or ratio =  1.703
!   real(wp), dimension(0:4), parameter :: c = (/1.0_wp,36.0_wp/17.0_wp,126.0_wp/85.0_wp, &
!                                                84.0_wp/221.0_wp,63.0_wp/2431.0_wp/)
!   real(wp), dimension(0:4), parameter :: d = (/2.0_wp,182.0_wp/51.0_wp,166.0_wp/85.0_wp, &
!                                              2578.0_wp/7735.0_wp,32768.0_wp/3828825.0_wp/)
!   real(wp),                 parameter :: eps = 6.0e-02_wp

!  Logarithmic expansion valid to us = 0.135, or ratio =  2.165
!   real(wp), dimension(0:5), parameter :: c = (/1.0_wp, 55.0_wp/21.0_wp, 330.0_wp/133.0_wp, &
!                                       330.0_wp/323.0_wp, 55.0_wp/323.0_wp, 33.0_wp/4199.0_wp/)
!   real(wp), dimension(0:5), parameter :: d = (/2.0_wp, 32.0_wp/7.0_wp, 3092.0_wp/855.0_wp, &
!                       7808.0_wp/6783.0_wp, 17926.0_wp/142443.0_wp, 131072.0_wp/61108047.0_wp/)
!   real(wp),                 parameter :: eps = 1.1e-01_wp

!  Logarithmic expansion valid to us = 0.221, or ratio =  2.776
!   real(wp), dimension(0:6), parameter :: c = (/1.0_wp,78.0_wp/25.0_wp,429.0_wp/115.0_wp,  &
!                                                1716.0_wp/805.0_wp, 1287.0_wp/2185.0_wp,   &
!                                                2574.0_wp/37145.0_wp, 429.0_wp/185725.0_wp/)
!   real(wp), dimension(0:6), parameter :: d = (/2.0_wp,418.0_wp/75.0_wp,3324.0_wp/575.0_wp,  &
!                                               55116.0_wp/20125.0_wp,399118.0_wp/688275.0_wp,&
!                                               1898954.0_wp/42902475.0_wp,                   &
!                                               2097152.0_wp/3904125225.0_wp/)
!   real(wp),                 parameter :: eps = 2.0e-01_wp

!  Logarithmic expansion valid to us = 0.315, or ratio =  3.560
!   real(wp), dimension(0:7), parameter :: c = (/1.0_wp,105.0_wp/29.0_wp,455.0_wp/87.0_wp,  &
!                                                1001.0_wp/261.0_wp,1001.0_wp/667.0_wp,    &
!                                                1001.0_wp/3335.0_wp,1001.0_wp/38019.0_wp, &
!                                                143.0_wp/215441.0_wp/)
!   real(wp), dimension(0:7), parameter :: d = (/2.0_wp,572.0_wp/87.0_wp,3674.0_wp/435.0_wp,  &
!                                                3256.0_wp/609.0_wp, 31054.0_wp/18009.0_wp,   &
!                                        86644.0_wp/330165.0_wp, 3622802.0_wp/244652265.0_wp, &
!                                                8388608.0_wp/62386327575.0_wp /)
!   real(wp),                 parameter :: eps = 3.0e-01_wp



    real(wp)              :: xi, gs, us, ave

    real(wp)              :: Logarithmic_Average

    xi  = a/b
    gs  = (xi-1.0_wp)/(xi+1.0_wp)
    us  = gs*gs
    ave = a + b 

!   if(                        (us <= 1.0e-8_wp)) Log_Ave_Counter( 1) = Log_Ave_Counter( 1) + 1
!   if((1.0e-8_wp <= us) .and. (us <= 1.0e-7_wp)) Log_Ave_Counter( 2) = Log_Ave_Counter( 2) + 1
!   if((1.0e-7_wp <= us) .and. (us <= 1.0e-6_wp)) Log_Ave_Counter( 3) = Log_Ave_Counter( 3) + 1
!   if((1.0e-6_wp <= us) .and. (us <= 1.0e-5_wp)) Log_Ave_Counter( 4) = Log_Ave_Counter( 4) + 1
!   if((1.0e-5_wp <= us) .and. (us <= 1.0e-4_wp)) Log_Ave_Counter( 5) = Log_Ave_Counter( 5) + 1
!   if((1.0e-4_wp <= us) .and. (us <= 1.0e-3_wp)) Log_Ave_Counter( 6) = Log_Ave_Counter( 6) + 1
!   if((1.0e-3_wp <= us) .and. (us <= 1.0e-2_wp)) Log_Ave_Counter( 7) = Log_Ave_Counter( 7) + 1
!   if((1.0e-2_wp <= us) .and. (us <= 1.0e-1_wp)) Log_Ave_Counter( 8) = Log_Ave_Counter( 8) + 1
!   if((1.0e-1_wp <= us) .and. (us <= 1.0e-0_wp)) Log_Ave_Counter( 9) = Log_Ave_Counter( 9) + 1
!   if((      eps <= us)                        ) Log_Ave_Counter(10) = Log_Ave_Counter(10) + 1

    if(us <= eps) then
      Logarithmic_Average =                                               &
!     ave * (c(0)-us*c(1)) / &
!           (d(0)-us*d(1))
!     ave * (c(0)-us*(c(1)-us*c(2))) / &
!           (d(0)-us*(d(1)-us*d(2)))
      ave * (c(0)-us*(c(1)-us*(c(2)-us*c(3)))) / &
            (d(0)-us*(d(1)-us*(d(2)-us*d(3))))
!     ave * (c(0)-us*(c(1)-us*(c(2)-us*(c(3)-us*c(4))))) / &
!           (d(0)-us*(d(1)-us*(d(2)-us*(d(3)-us*d(4)))))
!     ave * (c(0)-us*(c(1)-us*(c(2)-us*(c(3)-us*(c(4)-us*c(5)))))) / &
!           (d(0)-us*(d(1)-us*(d(2)-us*(d(3)-us*(d(4)-us*d(5))))))
!     ave * (c(0)-us*(c(1)-us*(c(2)-us*(c(3)-us*(c(4)-us*(c(5)-us*c(6))))))) / &
!           (d(0)-us*(d(1)-us*(d(2)-us*(d(3)-us*(d(4)-us*(d(5)-us*d(6)))))))
!     ave * (c(0)-us*(c(1)-us*(c(2)-us*(c(3)-us*(c(4)-us*(c(5)-us*(c(6)-us*c(7)))))))) / &
!           (d(0)-us*(d(1)-us*(d(2)-us*(d(3)-us*(d(4)-us*(d(5)-us*(d(6)-us*d(7))))))))
    else
      Logarithmic_Average = ave * gs / log(xi) 
    endif

    end function Logarithmic_Average

!===================================================================================================

  pure function Exp_Series(x)

    implicit none

    real(wp),  intent(in)   :: x

    real(wp), dimension(0:4), parameter :: a = (/139230.0_wp  ,16380.0_wp  ,      &
                                                    238.875_wp,    0.875_wp,      &
                                                  7.10227272727272727272727e-4_wp/)
    real(wp), dimension(0:4), parameter :: b = (/-69615.0_wp   ,-2388.75_wp,      &
                                                    -17.0625_wp,-   0.03125_wp,   &
                                                   -7.891414141414141e-6_wp/)

    real(wp)  :: Exp_Series

    real(wp)  :: t1, t2, x2

      x2 = x*x

      t1 = (a(0) + x2*(a(1) + x2*(a(2) + x2*(a(3) + x2*a(4)))))   ;
      t2 = (b(0) + x2*(b(1) + x2*(b(2) + x2*(b(3) + x2*b(4)))))*x ;

      Exp_Series = (t1 - t2) / (t1 + t2)

  end function Exp_Series

!===================================================================================================

  subroutine EntropyConsistentFlux_Vectors(vl,vr,neq,fx,fy,fz)

    !  Checked for accuracy :  10/25/2017

    ! this function calculates the normal entropy consistent
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the nondimensionalization employed
    ! herein and follows directly from the work of Ismail and Roe,
    ! DOI: 10.1016/j.jcp.2009.04.021

    use nsereferencevariables, only: gM2I, gm1og, gp1og, gm1M2

    implicit none

    ! Arguments
    ! =========
    integer,  intent(in)                   :: neq      ! number of equations

    real(wp), intent(in), dimension(neq)   :: vl, vr   ! left and right states

    real(wp), intent(out), dimension(neq)  :: fx,fy,fz ! left and right states

    ! Local Variables
    ! ===============

    real(wp), dimension(neq) :: vhat           ! temporary variables

    real(wp) :: root_Tl  , root_Tr             !   sqrt[T_i], i=L,R
    real(wp) :: root_Tl_I, root_Tr_I           ! 1/sqrt[T_i], i=L,R

    real(wp) :: tinvav, tinvavinv              ! average of inverse temperature and its inverse

    real(wp) :: s1, s2                         ! Logarithmic averages of density and temperature

    real(wp) :: P, T                           ! normal mass flux (mdot), Pressure,  Temperature

    continue

    root_Tl   = sqrt(vl(5))    ; root_Tr   = sqrt(vr(5))            ! Sqrt[T_i] 
    root_Tl_I = 1.0_wp/root_Tl ; root_Tr_I = 1.0_wp/root_Tr         ! Sqrt[T_i]^{-1}
    tinvav    = root_Tl_I   + root_Tr_I
    tinvavinv = 1.0_wp/tinvav

    vhat(2:4) = (vl(2:4)*root_Tl_I + vr(2:4)*root_Tr_I)*tinvavinv   ! velocity

    P = (vl(1)*root_Tl + vr(1)*root_Tr) * tinvavinv                 ! pressure

    s1 = Logarithmic_Average(root_Tl*vl(1),root_Tr*vr(1))           ! Logarithmic average of rho/Sqrt(T)

    s2 = Logarithmic_Average(root_Tl_I    ,root_Tr_I    )           ! Logarithmic average of   1/Sqrt(T)

    vhat(1) = 0.5_wp *tinvav * s1                                   ! density

    T = 0.5_wp * (gm1og * P + gp1og * s1/s2) / vhat(1)              ! temperature

    vhat(5) = T + 0.5_wp * gm1M2 * dot_product(vhat(2:4),vhat(2:4)) ! total enthalpy

    fx(1) = vhat(1)*vhat(2)
    fx(2) = vhat(1)*vhat(2)*vhat(2) + P * gM2I
    fx(3) = vhat(1)*vhat(2)*vhat(3)
    fx(4) = vhat(1)*vhat(2)*vhat(4)
    fx(5) = vhat(1)*vhat(2)*vhat(5)

    fy(1) = vhat(1)*vhat(3)
    fy(2) = vhat(1)*vhat(3)*vhat(2)
    fy(3) = vhat(1)*vhat(3)*vhat(3) + P * gM2I
    fy(4) = vhat(1)*vhat(3)*vhat(4)
    fy(5) = vhat(1)*vhat(3)*vhat(5)

    fz(1) = vhat(1)*vhat(4)
    fz(2) = vhat(1)*vhat(4)*vhat(2)
    fz(3) = vhat(1)*vhat(4)*vhat(3)
    fz(4) = vhat(1)*vhat(4)*vhat(4) + P * gM2I
    fz(5) = vhat(1)*vhat(4)*vhat(5)

  end subroutine EntropyConsistentFlux_vectors

!===================================================================================================

  pure function Entropy_KE_Consistent_Flux_AVX_512(vl,vr,Jx)

    ! this function calculates the normal entropy consistent
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the strange nondimensionalization employed herein 
    ! The logarithmic means follow Ismail and Roe, DOI: 10.1016/j.jcp.2009.04.021
    ! while the KE / S flux follows Chandrashekar

    ! Logarithmic mean for arbitrary variable ``a'' => DOI: 10.1016/j.jcp.2009.04.021
    ! \tilde{a} = (aL - aR) / (ln(aL) - ln(aR))
    ! Define xi = aL/aR ; f(xi) = ((xi-1)/ (xi+1)) ; u(xi) = f * f
    ! \tilde{a} = (aL + aR) / ln(xi) * f(xi) 
    ! \tilde{a} uses the asymptotic expansion of Ln() for values of () -> 1
    
    ! Both the ``rho'' or ``T'' logarithmic means are formed

    use nsereferencevariables, only: gM2I, gm1M2, gamI

    implicit none
    ! Arguments
    ! =========
    real(wp), intent(in), dimension(8,5) :: vl, vr                  ! left and right states
    real(wp), intent(in), dimension(8,3) :: Jx                      ! metrics scaled by jacobian

!   Logarithmic expansion valid to us = 0.135, or ratio =  2.165
    real(wp), dimension(0:5), parameter :: c = (/1.0_wp, 55.0_wp/21.0_wp, 330.0_wp/133.0_wp, &
                                        330.0_wp/323.0_wp, 55.0_wp/323.0_wp, 33.0_wp/4199.0_wp/)

    real(wp), dimension(0:5), parameter :: d = (/2.0_wp, 32.0_wp/7.0_wp, 3092.0_wp/855.0_wp, &
                        7808.0_wp/6783.0_wp, 17926.0_wp/142443.0_wp, 131072.0_wp/61108047.0_wp/)

    real(wp),                 parameter :: eps = 1.1e-01_wp

    ! Function
    ! ========
    real(wp), dimension(8,5)             :: Entropy_KE_Consistent_Flux_AVX_512

    ! Local temporary Variables
    ! ===============
    real(wp), dimension(8,5) :: vave

    real(wp), dimension(8)   :: mdot, P, rhotil, Btil, bL, bR, Keave  ! normal mass flux (mdot), Pressure, Logave density \& temperature
    real(wp), dimension(8)   :: xi, gs, us, ave, a, b

    integer                  :: i

    vave(:,1:5) = 0.5_wp * (vl(:,1:5) + vr(:,1:5))                  ! Average Density & velocity

    bL(:) = 1.0_wp/vl(:,5) ; bR(:) = 1.0_wp/vr(:,5)                 ! inverse Temperature (defined as B)

    P(:)  = 2.0_wp * gM2I * vave(:,1) / (bL(:) + bR(:))             ! Pressure

    Keave(:) = 0.25_wp * ( + vL(:,2)*vL(:,2)+vL(:,3)*vL(:,3)+vL(:,4)*vL(:,4)  &    ! Average KE 
                           + vR(:,2)*vR(:,2)+vR(:,3)*vR(:,3)+vR(:,4)*vR(:,4) )

    xi(:)  = vL(:,1)/vR(:,1)
    gs(:)  = (xi(:)-1.0_wp)/(xi(:)+1.0_wp)
    us(:)  = gs(:)*gs(:)
    ave(:) = a(:) + b(:)

    rhotil(:) = ave(:) * (c(0)-us(:)*(c(1)-us(:)*(c(2)-us(:)*(c(3)-us(:)*(c(4)-us(:)*c(5)))))) / &
                         (d(0)-us(:)*(d(1)-us(:)*(d(2)-us(:)*(d(3)-us(:)*(d(4)-us(:)*d(5))))))

    if(maxval(us(:)) >= eps) then
      do i = 1,8
        if(us(i) >= eps) then 
          rhotil(i) = ave(i) * gs(i) / log(xi(i))
        endif
      enddo
    endif
 
    xi(:)  = bL(:)/bR(:)
    gs(:)  = (xi(:)-1.0_wp)/(xi(:)+1.0_wp)
    us(:)  = gs(:)*gs(:)
    ave(:) = a(:) + b(:)

    Btil(:) =   ave(:) * (c(0)-us(:)*(c(1)-us(:)*(c(2)-us(:)*(c(3)-us(:)*(c(4)-us(:)*c(5)))))) / &
                         (d(0)-us(:)*(d(1)-us(:)*(d(2)-us(:)*(d(3)-us(:)*(d(4)-us(:)*d(5))))))

    if(maxval(us(:)) >= eps) then
      do i = 1,8
        if(us(i) >= eps) then 
          Btil(i) = ave(i) * gs(i) / log(xi(i))
        endif
      enddo
    endif

    mdot(:) = rhotil(:)*(vave(:,2)*Jx(:,1)+vave(:,3)*Jx(:,2)+vave(:,4)*Jx(:,3))     ! normal mass flow rate

    Entropy_KE_Consistent_Flux_AVX_512(:,1) = mdot(:)
    Entropy_KE_Consistent_Flux_AVX_512(:,2) = mdot(:)*vave(:,2) + Jx(:,1) * P(:)
    Entropy_KE_Consistent_Flux_AVX_512(:,3) = mdot(:)*vave(:,3) + Jx(:,2) * P(:)
    Entropy_KE_Consistent_Flux_AVX_512(:,4) = mdot(:)*vave(:,4) + Jx(:,3) * P(:)
    Entropy_KE_Consistent_Flux_AVX_512(:,5) = mdot(:)*(gamI/Btil(:) - gm1M2*Keave(:)) &
                                            + gm1M2 * ( Entropy_KE_Consistent_Flux_AVX_512(:,2) * vave(:,2) + &
                                                        Entropy_KE_Consistent_Flux_AVX_512(:,3) * vave(:,3) + &
                                                        Entropy_KE_Consistent_Flux_AVX_512(:,4) * vave(:,4) )

  end function Entropy_KE_Consistent_Flux_AVX_512

!===================================================================================================

! pure function Entropy_KE_Consistent_Flux(vl,vr,Jx,nq)
       function Entropy_KE_Consistent_Flux(vl,vr,Jx,nq)

    ! this function calculates the normal entropy consistent
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the strange nondimensionalization employed herein 
    ! The logarithmic means follow Ismail and Roe, DOI: 10.1016/j.jcp.2009.04.021
    ! while the KE / S flux follows Chandrashekar

!   Logarithmic mean for arbitrary variable ``a'' => DOI: 10.1016/j.jcp.2009.04.021
!   \tilde{a} = (aL - aR) / (ln(aL) - ln(aR))
!   Define xi = aL/aR ; f(xi) = ((xi-1)/ (xi+1)) ; u(xi) = f * f
!   \tilde{a} = (aL + aR) / ln(xi) * f(xi) 
!   \tilde{a} uses the asymptotic expansion of Ln() for values of () -> 1
!   ln(xi)    = 2 f(1 + u / 3 + u u / 5 + u u u / 7 + u u u u / 9 + ... ) 
!   Define  F =    (1 + u / 3 + u u / 5 + u u u / 7 + u u u u / 9 )
!   \tilde{a} = (aL + aR)/2 * 1/(al F + (1-al) ln(xi)/2/f) 
!
!   Both the ``rho'' or ``T'' logarithmic means are formed

    use nsereferencevariables, only: gM2I, gm1M2, gamI!, gm1og

    implicit none
    ! Arguments
    ! =========
    ! number of equations
    integer,  intent(in)                :: nq
    ! left and right states
    real(wp), intent(in), dimension(nq) :: vl, vr
    ! metrics scaled by jacobian
    real(wp), intent(in), dimension(3)  :: Jx

!   real(wp),             dimension(nq) :: wl, wr
!   real(wp)                            :: err
  
    ! Function
    ! ========
    real(wp), dimension(nq) :: Entropy_KE_Consistent_Flux
  
    ! Local temporary Variables
    ! ===============
    real(wp), dimension(nq) :: vave

    real(wp) :: mdot, P, rhotil, Btil, bL, bR, Keave                ! normal mass flux (mdot), Pressure, Logave density and temperature

    vave(1:5) = 0.5_wp * (vl(1:5) + vr(1:5))                        ! Average Density & velocity

    bL = 1.0_wp/vl(5) ; bR = 1.0_wp/vr(5)                           ! inverse Temperature (defined as B)

    P  = 2.0_wp * gM2I * vave(1) / (bL + bR)                        ! Pressure

    Keave = 0.25_wp * ( + vL(2)*vL(2)+vL(3)*vL(3)+vL(4)*vL(4)  &    ! Average KE 
                        + vR(2)*vR(2)+vR(3)*vR(3)+vR(4)*vR(4) )

    rhotil = Logarithmic_Average(vl(1),vr(1))                       ! logarithmic average of density

      Btil = Logarithmic_Average(   bL,   bR)                       ! logarithmic average of  B = 1/T

      mdot = rhotil*(vave(2)*Jx(1)+vave(3)*Jx(2)+vave(4)*Jx(3))     ! normal mass flow rate

    Entropy_KE_Consistent_Flux(1) = mdot
    Entropy_KE_Consistent_Flux(2) = mdot*vave(2) + Jx(1) * P
    Entropy_KE_Consistent_Flux(3) = mdot*vave(3) + Jx(2) * P
    Entropy_KE_Consistent_Flux(4) = mdot*vave(4) + Jx(3) * P
    Entropy_KE_Consistent_Flux(5) = mdot*(gamI/Btil - gm1M2*Keave) &
                                  + gm1M2 * ( Entropy_KE_Consistent_Flux(2) * vave(2) + &
                                              Entropy_KE_Consistent_Flux(3) * vave(3) + &
                                              Entropy_KE_Consistent_Flux(4) * vave(4) )

!   !  Check if it satisfies the entropy shuffle   DelW^T f = Delpsi
!   call primitive_to_entropy(vL,wL,nq)
!   call primitive_to_entropy(vR,wR,nq)

!   err = dot_product(wL-wR,Entropy_KE_Consistent_Flux) - gm1og* (           &
!                   ((vl(2)-vR(2))*vave(1) + (vl(1)-vR(1))*vave(2))*Jx(1) +  &
!                   ((vl(3)-vR(3))*vave(1) + (vl(1)-vR(1))*vave(3))*Jx(2) +  &
!                   ((vl(4)-vR(4))*vave(1) + (vl(1)-vR(1))*vave(4))*Jx(3) )
!   if(abs(err) >= 1.0e-14_wp) then
!     write(*,*)'error in Chandreshekar',err
!     write(*,*)'alr, alB',alr,alb
!     write(*,*)'rho in Chandreshekar',vl(1),vr(1)
!     write(*,*)'  u in Chandreshekar',vl(2),vr(2)
!     write(*,*)'  v in Chandreshekar',vl(3),vr(3)
!     write(*,*)'  w in Chandreshekar',vl(4),vr(4)
!     write(*,*)'  T in Chandreshekar',vl(5),vr(5)
!   endif

    return
  end function Entropy_KE_Consistent_Flux

!===================================================================================================

! pure function diabolical_flux(vL,vR,wL,wR,nx,nq)
  function diabolical_flux(vL,vR,wL,wR,nx,nq)

    use nsereferenceVariables, only: gm1og, gm1M2,gM2I,gamI

      implicit none
      integer,                 intent(in) :: nq
      real(wp),                intent(in) :: nx(3)
      real(wp), dimension(nq), intent(in) :: vL,vR
      real(wp), dimension(nq), intent(in) :: wL,wR

      real(wp), dimension(nq)             :: diabolical_flux

      ! matrices for minimization 
      real(wp), dimension(nq,3)  :: Flux_Vec
      real(wp), dimension(2,nq)  :: delta
      real(wp), dimension(2,3)   :: A
      real(wp), dimension(3)     :: x
      real(wp), dimension(2)     :: b

      ! Average Velocity
      real(wp), dimension(3)     :: vA

      real(wp) :: unL, unR
      real(wp) :: mdotL, mdotR, mdotA
      real(wp) :: pL, pR, pA, tA, En, keL, keR

      keL = dot_product(vL(2:4),vL(2:4)) * 0.5_wp
      keR = dot_product(vR(2:4),vR(2:4)) * 0.5_wp

      unL = dot_product(nx,vL(2:4))         ;
      unR = dot_product(nx,vR(2:4))         ;

      mdotL = vL(1)*unL ; mdotR = vR(1)*unR ;

      pL  = vL(1)*vL(5) ; pR  = vR(1)*vR(5) ;

      mdotA = 0.5_wp * (mdotL   + mdotR  )  ;
      vA(:) = 0.5_wp * (vL(2:4) + vR(2:4))  ;
      pA    = 0.5_wp * (pL + pR)            ;
      tA    = 0.5_wp * (vL(5)+vR(5))        ;

      En = tA * gamI  + 0.5_wp * gm1M2 * dot_product(vA,vA)

      Flux_Vec(1,1) = mdotL
      Flux_Vec(2,1) = mdotL * vL(2) + nx(1) * pL * gM2I
      Flux_Vec(3,1) = mdotL * vL(3) + nx(2) * pL * gM2I
      Flux_Vec(4,1) = mdotL * vL(4) + nx(3) * pL * gM2I
      Flux_Vec(5,1) = mdotL * ( vL(5) + gm1M2 * keL)

      Flux_Vec(1,2) = mdotR
      Flux_Vec(2,2) = mdotR * vR(2) + nx(1) * pR * gM2I
      Flux_Vec(3,2) = mdotR * vR(3) + nx(2) * pR * gM2I
      Flux_Vec(4,2) = mdotR * vR(4) + nx(3) * pR * gM2I
      Flux_Vec(5,2) = mdotR * ( vR(5) + gm1M2 * keR)

      Flux_Vec(1,3) = mdotA
      Flux_Vec(2,3) = mdotA * vA(1) + nx(1) * pA * gM2I
      Flux_Vec(3,3) = mdotA * vA(2) + nx(2) * pA * gM2I
      Flux_Vec(4,3) = mdotA * vA(3) + nx(3) * pA * gM2I
      Flux_Vec(5,3) = mdotA * En + gm1og*0.5_wp*(pL*unR+pR*unL)

      delta(1,1)    = - (keL   - keR  )
      delta(1,2)    =   (vL(2) - vR(2))
      delta(1,3)    =   (vL(3) - vR(3))
      delta(1,4)    =   (vL(4) - vR(4))
      delta(1,5)    =   0.0_wp

      delta(2,:)    =  wL(:) - wR(:)

      b(1)          =  (   unL -   unR) * pA * gM2I
      b(2)          =  ( mdotL - mdotR)

      A = matmul(delta,Flux_Vec)

!     Minimize somehow

      diabolical_flux(:) = matmul(Flux_Vec,x)

  end function diabolical_flux

!===================================================================================================

  pure function dUdV(Vin,neq)
    ! Checked MHC 08_09_13
    ! this function calculates the jacobian of the
    ! conserved variables with respect to the primitive
    ! variables
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neq
    ! primitive variables
    real(wp), intent(in) :: Vin(neq)

    ! output jacobian
    real(wp) :: dUdV(neq,neq)

    ! local convenience
    real(wp) :: ht

    dUdV = 0.0_wp

    ! total enthalpy
    ht = vin(5) + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4))

    ! continuity
    dUdV(1,1) = 1.0_wp
    ! momentum
    dUdV(2:4,1) = vin(2:4)
    dUdV(2,2)   = vin(1)
    dUdV(3,3)   = vin(1)
    dUdV(4,4)   = vin(1)
    ! energy
    dUdV(5,1)   = ht - gm1og*vin(5)
    dUdV(5,2:4) = gm1M2*vin(1)*vin(2:4)
    dUdV(5,5)   = gm1og*vin(1)/gm1

    return
  end function dUdV

  !============================================================================

  pure function dVdU(Vin,neq)
    ! Checked MHC 08_09_13
    ! this function calculates the jacobian of the
    ! primitive variables with respect to the conserved
    ! variables
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neq
    ! primitive variables
    real(wp), intent(in) :: Vin(neq)

    ! output jacobian
    real(wp) :: dVdU(neq,neq)

    real(wp) :: rhoinv
    real(wp) :: ht

    dVdU = 0.0_wp

    ! total enthalpy
    ht = vin(5) + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4))
    ! dv/du
    rhoinv = 1.0_wp/vin(1)

    dVdU(1,1)   = +1.0_wp
    dVdU(2:4,1) = -rhoinv*vin(2:4)
    dVdU(2,2)   = +rhoinv
    dVdU(3,3)   = +rhoinv
    dVdU(4,4)   = +rhoinv
    dVdU(5,1)   = +rhoinv*(-vin(5)+0.5_wp*gM2*gm1*dot_product(vin(2:4),vin(2:4)))
    dVdU(5,2:4) = -rhoinv*vin(2:4)*gm1*gM2
    dVdU(5,5)   = +rhoinv*gamma0

    return
  end function dVdU

  !============================================================================

  pure function dWdU(Vin,neqin)
    ! Checked 11-1-2013
    ! this function calculates the jacobian of the
    ! entropy variables with respect to the conserved
    ! variables using the chain rule, dw/du = dw/dv*dv/du.
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neqin
    ! primitive variables
    real(wp), intent(in) :: Vin(neqin)

    ! output jacobian
    real(wp) :: dWdU(neqin,neqin)

!   ! Local work jacobian
    real(wp) :: dWdV(neqin,neqin)

    real(wp) :: ht

    dWdV = 0.0_wp
    dWdU = 0.0_wp

    ! first calculate dW/dV
    ht = vin(5) + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4))

    dWdV(1,1) = +gm1og/vin(1)
    dWdV(1,2) = -gm1M2*vin(2)/vin(5)
    dWdV(2,2) = +gm1M2/vin(5)
    dWdV(1,3) = -gm1M2*vin(3)/vin(5)
    dWdV(3,3) = +gm1M2/vin(5)
    dWdV(1,4) = -gm1M2*vin(4)/vin(5)
    dWdV(4,4) = +gm1M2/vin(5)
    dWdV(1,5) = (gm1M2*dot_product(vin(2:4),vin(2:4))-ht)/(vin(5)*vin(5)) &
      & + 1.0_wp/vin(5) - 1.0_wp/(vin(5)*gamma0)
    dWdV(2:4,5) = -gm1M2*vin(2:4)/(vin(5)*vin(5))
    dWdV(5,5) = 1.0_wp/(vin(5)*vin(5))

    dWdU = matmul(dWdV,dVdU(vin,neqin))

    return
  end function dWdU

  !============================================================================

  function dUdW(Vin,neqin)
    ! this function calculates the jacobian of the
    ! conserved variables with respect to the entropy
    ! variables using the chain rule, du/dw = du/dv*dv/dw.

    implicit none
    ! number of equations
    integer, intent(in) :: neqin
    ! primitive variables
    real(wp), intent(in) :: Vin(neqin)

    ! output jacobian
    real(wp) :: dUdW(neqin,neqin)

    dUdW = 0.0_wp

    dUdW = matmul(dUdV(Vin,neqin),dVdW(Vin,neqin))

    return
  end function dUdW

  !============================================================================

  pure function dWdV(Vin,neqin)
    ! Checked MHC 08_09_13
    ! this function calculates the jacobian of the
    ! primitive variables with respect to the Entropy variables.
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neqin
    ! primitive variables
    real(wp), intent(in) :: Vin(neqin)

    ! output jacobian
    real(wp) :: dWdV(neqin,neqin)

    real(wp) :: TInv, T2Inv

    TInv  = 1.0_wp / vin(5)
    T2Inv = TInv*TInv

    dWdV  = 0.0_wp

    dWdV(1,1) = +gm1og/vin(1)
    dWdV(1,2) = -gm1M2*vin(2) * TInv
    dWdV(1,3) = -gm1M2*vin(3) * TInv
    dWdV(1,4) = -gm1M2*vin(4) * TInv
    dWdV(1,5) = +(-vin(5)/gamma0 + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4))) * T2Inv

    dWdV(2,2) = +gm1M2 * TInv
    dWdV(3,3) = +gm1M2 * TInv
    dWdV(4,4) = +gm1M2 * TInv

    dWdV(2,5) = -gm1M2 * vin(2) * T2Inv
    dWdV(3,5) = -gm1M2 * vin(3) * T2Inv
    dWdV(4,5) = -gm1M2 * vin(4) * T2Inv
    dWdV(5,5) =                   T2Inv

    return
  end function dWdV

  !============================================================================

  pure function dVdW(Vin,neqin)
    ! Checked MHC 08_09_13
    ! this function calculates the jacobian of the
    ! primitive variables with respect to the Entropy variables.
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neqin
    ! primitive variables
    real(wp), intent(in) :: Vin(neqin)

    ! output jacobian
    real(wp) :: dVdW(neqin,neqin)

    real(wp) :: gogm1, gm1I

    gogm1 = 1.0_wp/gm1og
    gm1I  = 1.0_wp/gm1

    dVdW = 0.0_wp

    ! first calculate dv/dW
    dVdW(1,1) = gogm1*vin(1)
    dVdW(1,2) = gogm1*vin(1)*vin(2)
    dVdW(1,3) = gogm1*vin(1)*vin(3)
    dVdW(1,4) = gogm1*vin(1)*vin(4)
    dVdW(1,5) = vin(1)*(gm1I*vin(5) + gM2*0.5_wp*dot_product(vin(2:4),vin(2:4)))

    dVdW(2,2) = vin(5)/gm1M2
    dVdW(2,5) = vin(5)*vin(2)

    dVdW(3,3) = vin(5)/gm1M2
    dVdW(3,5) = vin(5)*vin(3)

    dVdW(4,4) = vin(5)/gm1M2
    dVdW(4,5) = vin(5)*vin(4)

    dVdW(5,5) = vin(5)*vin(5)

    return
  end function dVdW

  !============================================================================

  pure subroutine roeavg(VLin,VRin,Vav,nq)
    ! this subroutine calculates the roe average
    use nsereferencevariables
    implicit none
    ! Arguments
    ! =========
    !
    ! number of equations
    integer, intent(in) :: nq
    ! left and right states
    real(wp), intent(in), dimension(nq) :: VLin, VRin
    ! average state
    real(wp), dimension(nq), intent(out) :: Vav

    ! local variables
    ! =================
    !
    real(wp) :: roe, roep1i
    real(wp) :: htl, htr, htav

    Vav = zero
    ! density ratio
    roe = sqrt(VRin(1)/VLin(1))
    ! save for efficiency
    roep1i = one/(roe+one)
    ! roe averaged density
    Vav(1) = roe*VLin(1)
    ! roe averaged velocity (2:4) 
    Vav(2:4) = (roe*VRin(2:4)+VLin(2:4))*roep1i
    ! roe averaged total enthalpy
    htl = vlin(5) + gm1M2*0.5_wp*dot_product(vlin(2:4),vlin(2:4))
    htr = vrin(5) + gm1M2*0.5_wp*dot_product(vrin(2:4),vrin(2:4))
    htav = (roe*htr+htl)*roep1i
    ! corresponding roe averaged temperature
    vav(5) = htav - gm1M2*0.5_wp*dot_product(vav(2:4),vav(2:4))

    return
  end subroutine roeavg

  pure subroutine CharacteristicDecomp(Vav,nq,Le,Re,ev,Jx)
    ! this subroutine calculates the left and right eigenvector
    ! matrices and the eigenvalues of the flux jacobian 
    ! in the normal direction. The resulting right eigenvectors
    ! have the special magic property that R R^T = du/dw.
    use nsereferenceVariables, only: gM2, gm1M2, gamma0, gm1, gm1og
    implicit none
    ! number of equations
    integer,intent(in) :: nq
    ! input primitive variables
    real(wp), dimension(nq), intent(in) :: Vav
    ! eigenvalues of df/du
    real(wp), dimension(nq), intent(out) :: ev
    ! left and right eigenvectors of df/du, respectively.
    real(wp), dimension(nq,nq), intent(out) :: Le, Re
    ! normal vector
    real(wp), intent(in) :: Jx(3)

    real(wp) :: rhoinv, lgm1, lgm1og, Tinv, lginv
    real(wp) :: cav, cavi, cav2, cav2i
    real(wp) :: un
    real(wp) :: mattmp(nq,nq), tmp, tmp2, ltmp, ltmpi
    real(wp) :: sqrt2, sqrt2i
    real(wp) :: g1, g2, g1i, g2i

    ! magnitude of normal  
    ltmp = sqrt(dot_product(jx,jx))
    ! inverse magnitude of normal
    ltmpi = one/ltmp
    ! inverse of density
    rhoinv = one/vav(1)
    ! inverse of temperature
    Tinv = one/vav(5)

    ! convenience
    sqrt2 = sqrt(2.0_wp)
    sqrt2i = 1.0_wp/sqrt2
    lginv = one/gamma0
    lgm1 = gm1
    lgm1og = gm1og

    ! speed of sound
    cav2 = gamma0*vav(5)/gM2
    cav = sqrt(cav2)
    ! inverse speed of sound
    cav2i = one/cav2
    cavi = one/cav

    ! normal velocity
    un = dot_product(vav(2:4),Jx)

    ! more convenience
    g1 = sqrt(0.5_wp*vav(5)*vav(1)*cav2i/gm1M2)
    g1i = 1.0_wp/g1
    g2 = sqrt(1.0_wp/gm1M2)
    g2i = 1.0_wp/g2

    ! calculate eigenvalues
    ev(1)   = un - cav*ltmp
    ev(2)   = un + cav*ltmp
    ev(3:5) = un

    ! start by calculating the right eigenvectors 
    ! of the primtive flux jacobian dv/du*df/du*du/dv
    Re(:,:) = 0.0_wp;

    ! First Vector
    Re(1,1) = g1
    Re(2:4,1) = -g1*ltmpi*Jx(1:3)*cav*rhoinv
    Re(5,1) = g1*lgm1*vav(5)*rhoinv

    ! Second Vector
    Re(1,2) = g1
    Re(2:4,2) = g1*ltmpi*Jx(1:3)*cav*rhoinv
    Re(5,2) = g1*lgm1*vav(5)*rhoinv

    tmp2 = sqrt(lgm1*vav(5)*vav(1)*cav2i)
    tmp = sqrt(vav(5)*rhoinv)
    ! Third Vector
    Re(1,3) = g2*Jx(1)*ltmpi*tmp2
    Re(3,3) = -g2*Jx(3)*ltmpi*tmp
    Re(4,3) = g2*Jx(2)*ltmpi*tmp
    Re(5,3) = -g2*Jx(1)*vav(5)*rhoinv*ltmpi*tmp2

    ! Fourth Vector
    Re(1,4) = g2*Jx(2)*ltmpi*tmp2
    Re(2,4) = g2*Jx(3)*ltmpi*tmp
    Re(4,4) = -g2*Jx(1)*ltmpi*tmp
    Re(5,4) = -g2*Jx(2)*vav(5)*rhoinv*ltmpi*tmp2

    ! Fifth Vector
    Re(1,5) = g2*Jx(3)*ltmpi*tmp2
    Re(2,5) = -g2*Jx(2)*ltmpi*tmp
    Re(3,5) = g2*Jx(1)*ltmpi*tmp
    Re(5,5) = -g2*Jx(3)*vav(5)*rhoinv*ltmpi*tmp2

    ! transform to eigenvectors of df/du by similarity
    ! dv/du*df/du*du/dv = S \Lambda S^{-1} so
    ! df/du = R \Lambda R^{-1} = du/dv*S \Lambda S^{-1}*dv/du
    mattmp = dUdV(Vav,nq)
    Re = matmul(mattmp,Re)

    ! now calculate the left eigenvectors 
    ! of the primtive flux jacobian dv/du*df/du*du/dv
    Le(:,:) = 0.0_wp

    tmp = sqrt(vav(1)*Tinv)
    tmp2 = 1.0_wp/sqrt(vav(5)*rhoinv*lgm1*cav2i)

    ! First Column
    Le(1:2,1) = g1i*0.5_wp*lginv
    Le(3:5,1) = g2i*lgm1og*rhoinv*Jx(1:3)*ltmpi*tmp2

    ! Second Column
    Le(1,2) = -0.5_wp*g1i*vav(1)*cavi*Jx(1)*ltmpi
    Le(2,2) =  0.5_wp*g1i*vav(1)*cavi*Jx(1)*ltmpi
    Le(4,2) =  g2i*Jx(3)*ltmpi*tmp
    Le(5,2) = -g2i*Jx(2)*ltmpi*tmp

    ! Third Column
    Le(1,3) = -0.5_wp*g1i*vav(1)*cavi*Jx(2)*ltmpi
    Le(2,3) =  0.5_wp*g1i*vav(1)*cavi*Jx(2)*ltmpi
    Le(3,3) = -g2i*Jx(3)*ltmpi*tmp
    Le(5,3) =  g2i*Jx(1)*ltmpi*tmp

    ! Fourth Column
    Le(1,4) = -0.5_wp*g1i*vav(1)*cavi*Jx(3)*ltmpi
    Le(2,4) =  0.5_wp*g1i*vav(1)*cavi*Jx(3)*ltmpi
    Le(3,4) =  g2i*Jx(2)*ltmpi*tmp
    Le(4,4) = -g2i*Jx(1)*ltmpi*tmp

    ! Fifth Column
    Le(1:2,5) = g1i*0.5_wp*lginv*vav(1)*Tinv
    Le(3:5,5) = -g2i*Tinv*lginv*Jx(1:3)*ltmpi*tmp2

    ! transform to eigenvectors of df/du by similarity
    ! dv/du*df/du*du/dv = S \Lambda S^{-1} so
    ! df/du = R \Lambda R^{-1} = du/dv*S \Lambda S^{-1}*dv/du
    mattmp = dVdU(Vav,nq)
    Le = matmul(Le,mattmp)

    return
  end subroutine CharacteristicDecomp
  
  !============================================================================

  subroutine nse_initializesemidiscretization()
    ! This routine allocates memory for the other variables
    ! in the discretization of space and time.
    use variables
    use referencevariables
    use controlvariables
    use nsereferencevariables, only : viscous
    implicit none

    integer :: nshell, nodesperelem_max

    nodesperelem_max = (npoly_max+1)**ndim

    allocate(  ughst(1:nequations,1:nghost)) ;   ughst = 0.0_wp           ! ghost cells for solution

    allocate(mutghst(             1:nghost)) ; mutghst = 0.0_wp           ! turbulent viscosity
    
    if(discretization == 'SSWENO') then
      allocate(ughstWENO(1:nequations,1:nghost))             ; ughstWENO          = 0.0_wp
      allocate(ughstWENO_partner(1:nequations,1:nghost))     ; ughstWENO_partner  = 0.0_wp

      nshell  = nfacesperelem*nodesperface
      allocate(ugWENO_partner(1:nequations,nshell,ihelems(1):ihelems(2)))  ; ugWENO_partner = 0.0_wp
    endif
    
    if (IMEX_penalty == 'implicit')  then
      if (viscous) then
        allocate(uelemghst(1:nequations,nghost_elem))
        uelemghst = 0.0_wp
      
        allocate(velemghst(1:nequations,nghost_elem))
        velemghst = 0.0_wp
      
        allocate(welemghst(1:nequations,nghost_elem))
        welemghst = 0.0_wp

        ! ghost points for exchanging r_x at the parallel interfaces
        allocate(r_x_ghst(3,3,1:nghost))
        r_x_ghst = 0.0_wp
      endif
    endif

    ! stores solution at T^n only used if step is rejected 
    allocate(uold(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    uold = 0.0_wp
    ! used in RK error estimation
    allocate(uhat(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    uhat = 0.0_wp


    ! penalty terms
    allocate(gsat(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    gsat = 0.0_wp
    ! flux divergence vectors -- one for each direction
    allocate(divf(1:nequations,ndim,1:nodesperelem_max,ihelems(1):ihelems(2)))
    divf = 0.0_wp

    ! entropy flux divergence used for error estimation
    allocate(divf_S(ndim,1:nodesperelem_max,ihelems(1):ihelems(2)))
    divf_S = 0.0_wp

    ! convective flux in each computational direction
    allocate(fg(1:nequations,ndim,1:nodesperelem_max,ihelems(1):ihelems(2)))
    fg = 0.0_wp
    ! viscous flux in each computational direction
    allocate(fvg(1:nequations,ndim,1:nodesperelem_max,ihelems(1):ihelems(2)))
    fvg = 0.0_wp
    ! variable gradients used for LDC/LDG approximation
    allocate(phig(1:nequations,3,1:nodesperelem_max,ihelems(1):ihelems(2)))
    phig = 0.0_wp
    allocate(phig_err(1:nequations,3,1:nodesperelem_max,ihelems(1):ihelems(2)))
    phig_err = 0.0_wp

    ! Gradient of the entropy variables in computational space for the
    ! calculation of the residual Jacobian matrix
    allocate(grad_w_jacobian(1:nequations,3,1:nodesperelem_max,ihelems(1):ihelems(2)))
    grad_w_jacobian = 0.0_wp
    
    ! ghost points for LDC/LDG approximation
    allocate(phighst(1:nequations,3,1:nghost))
    phighst = 0.0_wp
    ! shock sensor
    allocate(chig(1:nodesperelem_max,ihelems(1):ihelems(2)))
    chig = 0.0_wp
    ! artificial viscosity

    select case(RK_Method)

    case('Williamson_Low_Storage_45')

      ! Stores update
      allocate(du(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      du   = 0.0_wp
      ! local time derivative of conserved variables
      allocate(dudt(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      dudt = 0.0_wp

      ! local time derivative of entropy equation
      allocate(dudt_S(1:nodesperelem_max,ihelems(1):ihelems(2)))
      dudt_S = 0.0_wp

    case('Kraaij_LS_RK_35')

      ! Stores update
      allocate(du(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      du   = 0.0_wp
      ! local time derivative of conserved variables
      allocate(dudt(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      dudt = 0.0_wp

    case ('heun_method')
      ! Stores update
      allocate(du(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      du   = 0.0_wp
      ! local time derivative of conserved variables
      allocate(dudt(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      dudt = 0.0_wp

    case('IMEX_RK_46')
      allocate(Fimp(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2),6))
      allocate(Fexp(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2),6))
      allocate(uexp(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      allocate(non_lin_res(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    case('IMEX_RK_34')
      allocate(Fimp(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2),4))
      allocate(Fexp(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2),4))
      allocate(uexp(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
      allocate(non_lin_res(1:nequations,1:nodesperelem_max,ihelems(1):ihelems(2)))
    case default
      write(*,*)'Not a valid Temporal Scheme'
      write(*,*)'Check RK_Method in setup file'
      write(*,*)'Stopping'
    end select

    return
  end subroutine nse_initializesemidiscretization

  !============================================================================
  

  pure function dilatation(Axi,Jmat)
    real(wp), intent(in), dimension(3,3) :: Axi, Jmat

    real(wp) :: dilatation

    real(wp) :: Ax(3,3)
    integer :: dir,i

    do dir = 1,3
      do i = 1,3
        Ax(i,dir) = dot_product(Axi(i,:),Jmat(:,dir))
      end do
    end do

    dilatation = Ax(1,1)+Ax(2,2)+Ax(3,3)

    return
  end function dilatation

  pure function vorticity(Axi,Jmat)
    real(wp), intent(in), dimension(3,3) :: Axi, Jmat

    real(wp) :: vorticity(3)

    real(wp) :: Ax(3,3)

    integer :: dir,i

    do dir=1,3
      do i=1,3
        Ax(i,dir) = dot_product(Axi(i,:),Jmat(:,dir))
      end do
    end do

    vorticity(1) = Ax(3,2)-Ax(2,3)
    vorticity(2) = Ax(1,3)-Ax(3,1)
    vorticity(3) = Ax(2,1)-Ax(1,2)

    return
  end function vorticity

  !============================================================================

  subroutine nse_artificialviscosity()
    ! this subroutine is called to update the primitive
    ! and entropy variables. For viscous calculations, the
    ! entropy variable gradients are also updated according
    ! to the LDG approach. 
    use variables
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: iagrad,jagrad,dagrad
    use initcollocation, only: element_properties
    implicit none

    ! loop indices
    integer :: inode,ielem, jdir, idir
    integer :: jnode
    integer :: i

    ! LDC/LDG coefficient
    ! temporary arrays for phi and delta phi
    real(wp), allocatable :: phitmp(:,:)
    real(wp) :: theta, omega_loc(3), omegamag, sos, lmag
    real(wp) :: eta, delta
    real(wp), allocatable :: lambda(:), muvec(:)

    ! phitmp is calculated in computational space
    allocate(phitmp(nequations,3))
    ! dphi is calculated at faces
    allocate(lambda(nequations))
    allocate(muvec(nequations))
     mut = 0.0_wp
    chig = 0.0_wp

    ! loop over all elements
    do ielem = ihelems(1),ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      delta = 0.0_wp
      ! loop over every node in element
      do inode = 1, nodesperelem
        ! reinitialize computational gradient to zero
        phitmp = 0.0_wp
        ! loop over number of dependent elements in gradient
        do i = iagrad(inode), iagrad(inode+1)-1
          ! loop over dimensions
          do jdir = 1,ndim
            ! column/node from gradient operator in CSR format in
            ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
            jnode = jagrad(jdir,i)
            ! update gradient using coefficient and entropy variables at appropriate node
            phitmp(:,jdir) = phitmp(:,jdir) + dagrad(jdir,i) * wg(:,jnode,ielem)
          end do
        end do
        theta = dilatation(phitmp(2:4,:),r_x(:,:,inode,ielem))
        omega_loc = vorticity(phitmp(2:4,:),r_x(:,:,inode,ielem))
        omegamag = magnitude(omega_loc)
        sos = sqrt(gamma0*vg(5,inode,ielem)/gM2)
        lmag = 0.0_wp
        do idir = 1,ndim
          lmag = lmag + magnitude(x_r(:,idir,inode,ielem))*2.0_wp
        end do
        eta = -theta-max(5.0_wp*omegamag,0.05_wp*sos/lmag)
        chig(inode,ielem) = 0.5_wp*(eta+abs(eta))/(abs(eta)+1.0e-12_wp)
        delta =  max(delta, chig(inode,ielem))
      end do
    end do
    deallocate(phitmp,lambda,muvec)

    return
  end subroutine nse_artificialviscosity

  !============================================================================

  subroutine nse_reconcilestates()
    ! this subroutine is called to update the primitive
    ! and entropy variables. For viscous calculations, the
    ! entropy variable gradients are also updated according
    ! to the LDG approach. 
    use variables
    use SSWENOvariables,      only: WENO_type
    use controlvariables,     only: discretization
    use initcollocation,      only: element_properties
    use referencevariables
    use nsereferencevariables

    use mpimod, only: UpdateComm0DGhostData,UpdateComm1DGhostData, UpdateComm2DGhostData,     &
                      UpdateComm1DGhostDataWENO, UpdateCommShellGhostData
    use petscvariables, only: upetsc, phipetsc, ulocpetsc, philocpetsc, &
                              mutpetsc, mutlocpetsc, &
                              upetscWENO,  ulocpetscWENO,  &
                              upetscWENO_Shell, ulocpetscWENO_Shell
    implicit none

    ! loop indices
    integer :: ielem
    integer :: inode
    integer :: nshell

    call UpdateComm1DGhostData(ug, ughst, upetsc, ulocpetsc, size(ug,1), size(ug,2), size(ughst,2))

    ! Update ghost value in the interior plane for WENO p = 3
    if((discretization == 'SSWENO')) then

      call UpdateComm1DGhostDataWENO(ug, ughstWENO, upetscWENO, ulocpetscWENO, &
                                     size(ug,1), size(ug,2), size(ughstWENO,2))

      if((WENO_type == 'Neighbr_WENO')) call Planar_Interpolation_WENO()

      nshell = nodesperface*nfacesperelem

      call UpdateCommShellGhostData(ugWENO_partner, ughstWENO_partner, upetscWENO_Shell, &
                        ulocpetscWENO_Shell, size(ugWENO_partner,1),  &
                        size(ugWENO_partner,2), size(ughstWENO_partner,2))

    endif

    ! loop over all elements
    do ielem = ihelems(1),ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      do inode = 1,nodesperelem   ! loop over nodes and compute primitive and entropy variables
 
        call conserved_to_primitive(ug(:,inode,ielem), vg(:,inode,ielem), nequations ) ! (navierstokes)

        call primitive_to_entropy  (vg(:,inode,ielem), wg(:,inode,ielem), nequations ) ! (navierstokes)

      end do

    end do

    if (viscous) then

      if(turbulent_viscosity) call UpdateComm0DGhostData(mut, mutghst, mutpetsc, mutlocpetsc, size(mut,1), size(mutghst))

      call viscous_gradients()

      call UpdateComm2DGhostData(phig, phighst, phipetsc, philocpetsc,  &
                              size(phig,1), size(phig,2), size(phig,3), size(phighst,3))

    end if

    return
  end subroutine nse_reconcilestates

  !============================================================================

  subroutine nse_calc_dudt_LSRK(tin)

    ! Calculates the time derivative of the conserved variables at all nodes.

    use variables
    use referencevariables
    use initcollocation,       only: element_properties
    use nsereferencevariables, only: entropy_viscosity
    use collocationvariables,  only: iagrad,jagrad,dagrad,nnzgrad, &
                                     pinv, pmat, qmat, dmat, pvol, p_surf
    implicit none

    ! local time of evaluation (for RK schemes, this is the stage time)
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode,ielem, jdir
    integer :: n_pts_1d, n_pts_2d, n_pts_3d

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,           &
                           n_pts_1d=n_pts_1d,  &
                           n_pts_2d=n_pts_2d,  &
                           n_pts_3d=n_pts_3d,  &
                               pinv=pinv,      &
                               qmat=qmat,      &
                               dmat=dmat,      &
                            nnzgrad=nnzgrad,   &
                             iagrad=iagrad,    &
                             jagrad=jagrad,    &
                             dagrad=dagrad,    &
                               pmat=pmat,      &
                               pvol=pvol,      &
                             p_surf=p_surf)
          
      !  Calculate the elementwise Divergence  \/ * (F - Fv)

      call Flux_Divergence(tin, n_pts_1d, n_pts_2d, n_pts_3d,  &             !  result returned in divf(:,:,:,:)
                           pinv, qmat, dmat, iagrad, jagrad, dagrad, ielem)

      !  Form the elementwise SAT_Penalties

      call SAT_Penalty(tin, ielem, n_pts_1d, n_pts_2d, pinv )

      ! compute time derivative
        
      do inode = 1, n_pts_3d                                                 ! loop over all nodes in the element

          dudt(:,inode,ielem) =  ( - divf(:,1,inode,ielem) &                 ! Thus this is the dudt of u and NOT J u*
                                 & - divf(:,2,inode,ielem) &
                                 & - divf(:,3,inode,ielem) &
                                 & + gsat(:  ,inode,ielem) ) / Jx_r(inode,ielem) 


      end do

      ! Entropy equation
      if( entropy_viscosity .eqv. .true.) then

        call Entropy_Flux_Divergence(tin,ielem)   !  result in divF_S

        do inode = 1, nodesperelem
          ! reset the time derivative to zero
          dudt_S(inode,ielem) = 0.0_wp
          ! add the contribution from the flux divergence in each direction
  
          do jdir = 1,ndim
            dudt_S(inode,ielem) = dudt_S(inode,ielem) - divf_S(jdir,inode,ielem)
          end do

          dudt_S(inode,ielem) = (dudt_S(inode,ielem))/Jx_r(inode,ielem) + dot_product(wg(:,inode,ielem),dudt(:,inode,ielem))

        end do

      endif

    end do

    if( entropy_viscosity .eqv. .true.) write(*,*)'max entropy error',maxval(abs(dudt_S(:,:)))
         
    return
  end subroutine nse_calc_dudt_LSRK

  subroutine nse_calcembeddedspatialerror()
    ! this subroutine calculates the embedded error approximation
    ! using the solution, ug, and the embedded solution uhat.
    use variables, only: gsat, Jx_r
    use collocationvariables, only: pvol
    use controlvariables, only: verbose
    use initcollocation, only: element_properties
    use referencevariables
    use mpimod
    implicit none

    ! indices
    integer :: inode,ielem

    ! different error estimates
    real(wp) :: l2(2), linf
    real(wp) :: l2sum(2), linfmax
    ! local error
    real(wp), allocatable :: ex(:)
    integer :: ierr

    ! update primitives and entropy variables
    call nse_reconcilestates()

    allocate(ex(nequations))

    ! initialize errors to zero
      l2 = 0.0_wp
    linf = 0.0_wp

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem, pvol=pvol)
          
      ! loop over each index in the element
      do inode = 1, nodesperelem
        ! compute the local embedded error
        ex = gsat(:,inode,ielem)/Jx_r(inode,ielem)
        ! calculate linf contribution
        linf = max(linf,Jx_r(inode,ielem)*maxval(abs(ex)))
        ! calculate the integral contribution of l2 error. pvol*J is the volumetric
        ! integration weight.
        l2(1) = l2(1) + pvol(inode)*Jx_r(inode,ielem)*dot_product(ex,ex)/nequations
        l2(2) = l2(2) + pvol(inode)*Jx_r(inode,ielem)
      end do
    end do

    call mpi_allreduce(l2,l2sum,2, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)

    l2sum(1) = sqrt(l2sum(1)/l2sum(2))

    call mpi_allreduce(linf,linfmax,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,ierr)

    if(verbose .and. myprocid == 0) then
      write(*,101) l2sum(1), linfmax
    endif
    101 format('embedded l2 error: ',ES12.5,1X,'embedded linf error: ',ES12.5,1X)

    deallocate(ex)

    return
  end subroutine nse_calcembeddedspatialerror

  subroutine nse_calcembeddederror()
    ! this subroutine calculates the embedded error approximation
    ! using the solution, ug, and the embedded solution uhat.
    use variables, only: ug,uhat,vg,Jx_r
    use collocationvariables, only: pvol
    use controlvariables, only: verbose
    use initcollocation,  only: element_properties
    use referencevariables
    use mpimod
    implicit none

    ! indices
    integer :: inode,ielem

    ! different error estimates
    real(wp) :: l2(2) 
    real(wp) :: linf, sglob
    real(wp) :: l2sum(2), linfmax, sglobsum
    ! local error
    real(wp), allocatable :: ex(:)
    integer :: ierr

    ! update primitives and entropy variables
    call nse_reconcilestates()

    allocate(ex(nequations))

    ! initialize errors to zero
    l2 = 0.0_wp

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem, pvol=pvol)

      do inode = 1, nodesperelem
        ! compute the local embedded error
        ex = ug(:,inode,ielem) - uhat(:,inode,ielem)
        ! calculate linf contribution
        linf = max(linf,maxval(abs(ex)))
        ! calculate the integral contribution of l2 error. pvol*J is the volumetric
        ! integration weight.
        l2(1) = l2(1) + pvol(inode)*Jx_r(inode,ielem)*dot_product(ex,ex)/nequations
        l2(2) = l2(2) + pvol(inode)*Jx_r(inode,ielem)
        ! calculate the integral contribution to the global entropy 
        sglob = sglob - pvol(inode)*Jx_r(inode,ielem)*vg(1,inode,ielem) * &
          specificentropy(vg(:,inode,ielem),nequations)
      end do
    end do

    call mpi_allreduce(l2,l2sum,2, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)

    l2sum(1) = sqrt(l2sum(1)/l2sum(2))

    call mpi_allreduce(linf,linfmax,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,ierr)

    call mpi_allreduce(sglob,sglobsum,1, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)

    if(verbose .and. myprocid == 0) write(*,101) l2sum, linfmax, sglobsum
    101 format('l2 error: ',ES12.5,1X,'linf error: ',ES12.5,1X,'S: ',ES12.5)

    deallocate(ex)

    return
  end subroutine nse_calcembeddederror

  subroutine nse_calcerror(tin)
    ! This subroutine calculates the L1, L2 and Linfinity errors
    ! when the exact solution is known.
    use variables, only: ug, vg, xg, Jx_r, phig, mut!, Log_Ave_Counter
    use collocationvariables, only: pvol
    use initcollocation,  only: element_properties
    use referencevariables
    use mpimod
    implicit none

    ! Time of evaluation
    real(wp), intent(in) :: tin
    ! Indices
    integer :: inode, ielem, i
!   integer, dimension(10) :: Log_Ave_Counter_sum

    ! Errors
    real(wp) :: l1(2),    l2(2),    linf
    real(wp) :: l1sum(2), l2sum(2), linfmax
    ! Local error and temporary array
    real(wp), allocatable :: uex(:), vex(:), ftmp(:)
    ! Output unit
!   integer, allocatable, dimension(:) :: iunit
    ! Temporary for arbitrary normal vector
    real(wp) :: ctmp(3)
    integer :: ierr

    ! Global error norms file name
!   character(len=23) :: file_name_err

    ! Assign global error norms file name
!    file_name_err = 'global_error_norms.data'

    ! Arbitrary normal vector
    ctmp = 0.0_wp

    ! Update primitive variables
    call nse_reconcilestates()

    allocate(uex(nequations))
    allocate(vex(nequations))
    allocate(ftmp(nequations))

    ! Initialize error to zero
    l1 = 0.0_wp; l1sum = 0.0_wp ; l2 = 0.0_wp; l2sum = 0.0_wp ; linf = 0.0_wp ;

    ! Loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem, pvol=pvol)

      ! Loop over every node in the element
      do inode = 1, nodesperelem
        ! Compute the exact solution
        call InitialSubroutine( &
          vex,phig(:,:,inode,ielem),ftmp,ctmp,xg(:,inode,ielem),tin,nequations,ndim,mut(inode,ielem))
        call primitive_to_conserved(vex,uex,nequations)
        uex = uex - ug(:,inode,ielem)
        ! calculate the local error
        vex = vex - vg(:,inode,ielem)
        ! calculate the local contribution to the linf error
        linf = max(linf,maxval(abs(uex)))
        ! compute the addition of the error integral to the l2
        ! note that pvol*J is the volumetric integration weight
        l2(1) = l2(1) + pvol(inode)*Jx_r(inode,ielem)*dot_product(uex,uex)/nequations
        l2(2) = l2(2) + pvol(inode)*Jx_r(inode,ielem)
        ! compute the addition of the error integral to the l2
        ! note that pvol*J is the volumetric integration weight
        l1(1) = l1(1) + pvol(inode)*Jx_r(inode,ielem)*sum(abs(uex))/nequations
        l1(2) = l1(2) + pvol(inode)*Jx_r(inode,ielem)
      end do
    end do

    ! Sum all the L1 terms across all the processes
    call mpi_allreduce(l1,l1sum,2, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)


    ! Complete L1 norm
    l1sum(1) = l1sum(1)/l1sum(2)

    ! Sum all the L2 terms across all the processes
    call mpi_allreduce(l2,l2sum,2, &
      & MPI_DEFAULT_WP,MPI_SUM,petsc_comm_world,ierr)

    ! Complete L2 norm
    l2sum(1) = sqrt(l2sum(1)/l2sum(2))

    ! Sum all the LInf terms across all the processes
    call mpi_allreduce(linf,linfmax,1, &
      & MPI_DEFAULT_WP,MPI_MAX,petsc_comm_world,ierr)

    ! Reduce counter information for Log_Averages 
!   call mpi_allreduce(Log_Ave_Counter,Log_Ave_Counter_sum,8, &
!     & MPI_INT,MPI_SUM,petsc_comm_world,ierr)

    if (myprocid == 0) then

      write(*,*) 'P',nodesperedge-1
      write(*,*) 'L1   error: ', l1sum(1)
      write(*,*) 'L2   error: ', l2sum(1)
      write(*,*) 'Linf error: ', linfmax
      write(40,230)i,l1sum(1),l2sum(1),linfmax
      230  format('P',I2,1x,'  l1 error:  ', e17.10,1x,  &
                            '  l2 error:  ', e17.10,1x,  &
                            'linf error:  ', e17.10,1x)
!     write(*,'(10(f9.4,1x))') 1.0_wp*Log_Ave_Counter(:)/sum(Log_Ave_Counter(:))
    end if

    close(unit=40)
    deallocate(uex)
    deallocate(vex)
    deallocate(ftmp)

    return
  end subroutine nse_calcerror

  subroutine isentropicVortexFull(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut
    
    real(wp) :: epsvortex
    real(wp) :: Uinf
    real(wp) :: x0,y0
    real(wp) :: f
    real(wp) :: alpha, rin2
   

    y0 = +0.00_wp
    x0 = -0.25_wp
    Uinf = referencewavespeed

    epsvortex = 5.0_wp
    
    alpha = uniformFreeStreamAOA*pi/180._wp

    rin2 = ((xin(1)-x0)-Uinf*tin*cos(alpha))**2 + ((xin(2)-y0)-Uinf*tin*sin(alpha))**2

    f = 1.0_wp-rin2
    
    Vx = 0.0_wp
    ! density
!   tt = (one-gm1M2*epsvortex*epsvortex/(eight*pi*pi)*exp(f))
!   Vx(1) = (tt)**(one/gm1)
    Vx(1) = (one-gm1M2*epsvortex*epsvortex/(eight*pi*pi)*exp(f))**(one/gm1)
!     ! velocity
    Vx(2) = Uinf*cos(alpha) + epsvortex/(two*pi)*((xin(2)-y0)-Uinf*tin*sin(alpha))*exp(f*half)
    Vx(3) = Uinf*sin(alpha) - epsvortex/(two*pi)*((xin(1)-x0)-Uinf*tin*cos(alpha))*exp(f*half)
    ! Temperature
    Vx(5) = (one-gm1M2*epsvortex*epsvortex/ &
      & (eight*pi*pi)*exp(f))

    fv = 0.0_wp

    if(abs(Vx(4)) .gt. 1e-10) then
      write(*,*) 'rin2', rin2
      write(*,*) 'angle', uniformFreeStreamAOA
      write(*,*) 'pi', pi
      write(*,*) 'f', f
      write(*,*) 'Vx(1)', Vx(1)
      write(*,*) 'Vx(2)', Vx(2)
      write(*,*) 'Vx(3)', Vx(3)
      write(*,*) 'Vx(4)', Vx(4)
      write(*,*) 'Vx(5)', Vx(5)
    endif

    return
  end subroutine isentropicVortexFull

  subroutine supersonicvortexFull(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)

    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    real(wp)             :: r, ratio, alpha

    r = sqrt((xin(1)**2 + xin(2)**2))

    ratio = 1.0_wp / r

    alpha = atan2(xin(2),xin(1))

    Vx = 0.0_wp
    ! density
    Vx(1) = (1.0_wp + 0.5_wp*gm1M2*(1.0_wp - ratio**2) )**(1.0_wp/gm1)
!     ! velocity
    Vx(2) = - ratio*sin(alpha)
    Vx(3) = + ratio*cos(alpha)
    ! Temperature
    Vx(5) = Vx(1)**(gm1)

    fv = 0.0_wp

    return
  end subroutine supersonicvortexFull

  subroutine Potential_Flow_Around_cylinder(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
!
!   Assumes velocity profile satisfies potential flow around cylinder
!   T = constant
!   compressible Bernoulli's law to solve for density
!
    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    real(wp)             :: r, ratio, alpha
    real(wp)             :: v_rad, v_cir, v_mag
    
    r = sqrt((xin(1)**2 + xin(2)**2))

    ratio = 1.0_wp / r 

    alpha = atan2(xin(2),xin(1))
    v_rad = + (1.0_wp - ratio**2)*cos(alpha)
    v_cir = - (1.0_wp + ratio**2)*sin(alpha)


    v_mag = v_rad*v_rad + v_cir*v_cir

    Vx = 0.0_wp
    ! density
    Vx(1) = (1.0_wp + 0.5_wp*gM2) / (1.0_wp + 0.5_wp*gM2*v_mag)
!     ! velocity
    Vx(2) = v_rad * cos(alpha) - v_cir * sin(alpha)
    Vx(3) = v_rad * sin(alpha) + v_cir * cos(alpha)
    ! Temperature
    Vx(5) = 1.0_wp

    fv = 0.0_wp

    return
  end subroutine Potential_Flow_Around_cylinder

  subroutine ShockVortexInteraction(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut
    
    real(wp) :: epsvortex,rv
    real(wp) :: Uinf
    real(wp) :: x0,y0
    real(wp) :: f
    real(wp) :: alpha, rin2, gp1, M2
    real(wp) :: rho,u,v,p
   
    ! vortex radius
    rv  = 0.05_wp

    gp1 = gamma0 + 1.0_wp
    M2  = Mach0*Mach0

    y0   =  0.0_wp
    x0   = -0.25_wp
    Uinf = 1.0_wp

    if(xin(1).le.0.0_wp)then
      p   = 1.0_wp
      rho = 1.0_wp
      u   = Uinf
      v   = 0.0_wp
    else
      p   = (7.0_wp*M2-one)/(6.0_wp)
      rho = (6.0_wp*M2)/(M2+5.0_wp)
      u   = Uinf * (M2+5.0_wp) / (6.0_wp*M2)
      v   = 0.0_wp
    endif

    epsvortex = 2.7206_wp
    
    alpha = 0.0_wp

    rin2 = ((xin(1)-x0)-Uinf*tin*cos(alpha))**2 + ((xin(2)-y0)-Uinf*tin*sin(alpha))**2

    f = 1.0_wp-rin2/(rv*rv)
    
    Vx = 0.0_wp

    ! density
    if(xin(1).le.0.0_wp)then
      Vx(1) = (rho-gm1M2*epsvortex*epsvortex/(eight*pi*pi)*exp(f))**(one/gm1)
    else
      Vx(1) = rho
    endif

!     ! velocity
    if(xin(1).le.0.0_wp)then
      Vx(2) = u + epsvortex/(two*pi)*((xin(2)-y0)-Uinf*tin*sin(alpha))*exp(f*half)/rv
      Vx(3) = v - epsvortex/(two*pi)*((xin(1)-x0)-Uinf*tin*cos(alpha))*exp(f*half)/rv
    else
      Vx(2) = u
      Vx(3) = v
    endif

    ! Temperature
    if(xin(1).le.0.0_wp)then
      Vx(5) = p/rho - gm1M2*epsvortex*epsvortex/(eight*pi*pi)*exp(f)
    else
      Vx(5) = p/rho
    endif

    fv = 0.0_wp

    if(abs(Vx(4)) .gt. 1e-10) then
      write(*,*) 'rin2', rin2
      write(*,*) 'angle', uniformFreeStreamAOA
      write(*,*) 'pi', pi
      write(*,*) 'f', f
      write(*,*) 'Vx(1)', Vx(1)
      write(*,*) 'Vx(2)', Vx(2)
      write(*,*) 'Vx(3)', Vx(3)
      write(*,*) 'Vx(4)', Vx(4)
      write(*,*) 'Vx(5)', Vx(5)
    endif

    return
  end subroutine ShockVortexInteraction

  subroutine exact_Riemann(p1,ro1,p2,ro2,g,xb,x0,y1,y2,y3,t)

      implicit none
      integer   :: it

      real(wp), intent(in)    :: p1,ro1,p2,ro2,g,xb,x0,t
      real(wp), intent(out)   :: y1,y2,y3
      real(wp)                :: eps,z,p,p0,c1,c2,a1,a2,u,x

      !  adjust for initial position of diaphram
      x  = xb - x0
      
      eps=1.0_wp

      it=0
      c1=sqrt(g*p1/ro1)
      c2=sqrt(g*p2/ro2)
      p0=(p1*ro2*c2+p2*ro1*c1)/(ro1*c1+ro2*c2)
      p=0.0_wp
      do while(eps.gt.1.e-10_wp)
        a1=(g-1.0_wp)/(2*g)*ro1*c1* &
     &            (1.0_wp-p/p1)/(1.0_wp-(p/p1)**((g-1)/(2*g)))
        a2=sqrt(ro2*((g+1)/2*p+(g-1)/2*p2))
        p=(a2*p1+a1*p2)/(a1+a2)
        eps=abs(p-p0)
        p0=p
        it=it+1        

      enddo

      z= (p/p2 - 1.0_wp)
      u=(p1-p2)/(a1+a2)

!     Left state                                                        
       if(x.le.-c1*t)then

        y1=ro1
        y2=p1
        y3=0.0_wp
       endif

!     Expansion wave                                                    
       if(x.gt.-c1*t  .and.  x.le.((g+1)/2*u-c1)*t)then
        y1=ro1*c1**(2.0_wp/(1-g))* &
     &        ((g-1)/(g+1)*(2*c1/(g-1)-x/t))**(2.0_wp/(g-1))
        y2=p1/ro1**g*y1**g
        y3=  2.0_wp/(g+1.0_wp) * (c1 + x/t)
       endif

!     Contact discontinuity                                             
       if(x.gt.((g+1)/2*u-c1)*t  .and.  x.le.u*t)then
        y1=g*p/(c1-(g-1)/2*u)**2
        y2=p1/ro1**g*y1**g
        y3=c2 * z / (g * sqrt(1.0_wp + 0.5_wp * (g+1.0_wp)/g * z) )
       endif

!     Shock wave                                                        
       if(x.gt.u*t  .and.  x.le.a2/ro2*t)then
        y1=ro2*a2/(a2-ro2*u)
        y2=p1/ro1**g*(g*p/(c1-(g-1)/2*u)**2)**g
        y3=c2 * z / (g * sqrt(1.0_wp + 0.5_wp * (g+1.0_wp)/g * z) )
       endif

!     Right state                                                       
       if(x.gt.a2/ro2*t)then
        y1=ro2
        y2=p2
        y3=0.0_wp
       endif

  end subroutine exact_Riemann

!-----

  subroutine SodsProblemICBC(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)

    ! Load modules
    use nsereferencevariables
  
    implicit none

    integer,  intent(in)    :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out)   :: fv(neqin)
    real(wp), intent(in )   :: phi(neqin,3)
    real(wp), intent(in)    :: Jx(3)
    real(wp), intent(in)    :: xin(3)
    real(wp), intent(in)    :: tin, mut

    ! Local Variables
    ! ===============

!   real(wp) ::   U1,   U2
!   real(wp) ::   T1,   T2
    real(wp) :: rho1, rho2
    real(wp) ::   P1,   P2

    real(wp) :: cc, ss, xp
    real(wp) :: theta
    real(wp) :: den,pres,vel

!   Use a one-dimensional exact solution rotated into the orientation of the uniformfreestream flow
    theta = uniformFreeStreamAOA*pi/180._wp

    cc = cos(theta)
    ss = sin(theta)

    xp = xin(1)*cc+xin(2)*ss

!   Initial Left state
    rho1 = 1.0_wp
    P1   = 1.0_wp
!   U1   = 0.0_wp
!   T1   = P1 / rho1

!   Initial Right state
    rho2 = 0.125_wp
    P2   = 0.1_wp
!   U2   = 0.0_wp
!   T2   = P2 / rho2

    call exact_Riemann(P1,rho1,P2,rho2,gamma0,xp,membranelocation,den,pres,vel,tin)

    Vx = zero

    Vx(1) = den
    Vx(2) = vel*cc
    Vx(3) = vel*ss
    Vx(4) = 0.0_wp
    Vx(5) = pres/den

    fv(:) = zero

    return
  end subroutine SodsProblemICBC

  subroutine viscousShockFull(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    
    ! Load modules
    use nsereferencevariables
    
    use controlvariables, only : variable_viscosity
  
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut
    
    real(wp) :: gs(3,neqin)
    real(wp) :: Uinf
    real(wp) :: x0,y0,r0
    real(wp) :: f, xp
    real(wp) :: alph, vf, mdot, Mach, wave, M2, uhat, Rloc, gloc
    real(wp) :: tautmp(3,3)
    real(wp) :: u_xp
    integer :: meq, dir
    real(wp) :: cc, ss
    real(wp) :: theta
    real(wp) :: mu, kappa
    
    integer :: ns

    continue

    ns = neqin-4
    theta = uniformFreeStreamAOA*pi/180._wp
    r0 = membranelocation
    y0 = zero
    x0 = 0.0_wp
    Uinf = 1.0_wp
    M2 = Mach0*Mach0 !Uinf*Uinf
    mach = sqrt(M2)
    mdot = abs(Uinf)

    cc = cos(theta)
    ss = sin(theta)

    xp = xin(1)*cc+xin(2)*ss
    
    Vx = zero
    
    ! set mixture variable
    Rloc = 1.0_wp
    gloc = gamma0
    
    wave = referencewavespeed

    ! Set dynamic viscosity
    if (variable_viscosity .eqv. .true.) then
      mu = sutherland_law(Vx(5))
    else
      mu = 1.0_wp
    end if

    ! Set heat conductivity
    kappa = 1.0_wp

    vf = (gm1+two/M2)/(gamma0+one)
    
    alph = (four/three*mu/Re0)/mdot*(two*gamma0)/(gamma0+one)
    
    call rhalf(xp-r0-wave*tin,f,alph,vf)
    
    uhat = f*Uinf
    
    ! density
    Vx(1) = abs(mdot)/f
    ! velocity
    Vx(2) = (uhat+wave)*cc
    Vx(3) = (uhat+wave)*ss
    ! Temperature
    Vx(5) = (Uinf*Uinf + gm1M2*half*(Uinf*Uinf-uhat*uhat))/Rloc
    
    gs = zero
    
    tautmp = zero
    u_xp = (f-one)*(f-vf)/(alph*f)*Uinf
    tautmp(1,1) = mu/three*(four*u_xp*cc*cc-two*ss*ss*u_xp)
    tautmp(2,1) = mu*u_xp*(cc*ss+ss*cc)
    tautmp(1,2) = tautmp(2,1)
    tautmp(2,2) = mu/three*(four*u_xp*ss*ss-two*cc*cc*u_xp)
    tautmp(3,3) = -mu*two/three*u_xp
    
    do dir = 1,3
      do meq = 1,3
        gs(dir,meq+1) = tautmp(meq,dir)/Re0
      end do
      gs(dir,5) = dot_product(tautmp(:,dir),Vx(2:4))*gm1M2/Re0
    end do
    gs(1,5) = gs(1,5) - kappa/(Re0*Pr0)*gm1M2*uhat*u_xp*cc
    gs(2,5) = gs(2,5) - kappa/(Re0*Pr0)*gm1M2*uhat*u_xp*ss
!     gs(1,5) = tautmp(1,1)*gm1M2/Re0*wave
!     gs(1,5) = u_x*wave*kND(mix)/(Re0*Pr0)*gm1M2
!     gs(1,2) = tautmp(1,1)/Re0

    do meq = 1,neqin
      fv(meq) = dot_product(gs(:,meq),Jx)
    end do
    
    return
  end subroutine viscousShockFull

  !============================================================================
  
  !============================================================================

  subroutine kelvin_helmoholtz(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)

    ! Load module 
    use nsereferencevariables

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut
    
    real(wp) :: eps
    real(wp) :: rho_1, rho_2
    real(wp) :: u_1, u_2
    real(wp) :: v_1, v_2
    real(wp) :: w_1, w_2
    real(wp) :: p_1, p_2


    ! Set up the primitive variables of the two states
    ! ================================================
    ! Perturbation strength
    eps = 0.1_wp

    ! Density
    rho_1 = 2.0_wp
    rho_2 = 1.0_wp

    ! Velocity components
    u_1 = -0.5_wp + eps*sin(2.0_wp*pi*xin(1))
    u_2 = +0.5_wp + eps*sin(2.0_wp*pi*xin(1))
    
    v_1 = +0.0_wp + eps*sin(2.0_wp*pi*xin(2))
    v_2 = +0.0_wp + eps*sin(2.0_wp*pi*xin(2))
    
    w_1 = +0.0_wp
    w_2 = +0.0_wp
    
    p_1 = 2.5_wp
    p_2 = 2.5_wp
    
    ! Set primitive variables (density,velocity,temperature) used by the code 
    Vx = 0.0_wp
    
    if (xin(2) .gt. 0.25_wp .and. xin(2) .lt. 0.75_wp) then
      ! Density
      Vx(1) = rho_1  
     
      ! Velocity
      Vx(2) = u_1
      Vx(3) = v_1
      Vx(4) = w_1

      ! Temperature (T = p/(R*T), with R = one)
      Vx(5) = p_1 / (Vx(1)*one)
    else
      ! Density
      Vx(1) = rho_2  
     
      ! Velocity
      Vx(2) = u_2
      Vx(3) = v_2
      Vx(4) = w_2

      ! Temperature (T = p/(R*T), with R = one)
      Vx(5) = p_2 / (Vx(1)*one)

    end if

    fv = 0.0_wp

    return
  end subroutine kelvin_helmoholtz

  !============================================================================
  
  !============================================================================


  subroutine rhalf(xin,v,alpha,vf)
    !   solve equation by interval halving
    !
    !       f = EXP(2*(1-VF)*X/ALPH)-(V-1)**2/ABS(VF-V)**(2*VF)
    !       f = EXP(-2*(1-VF)*X/ALPH)-ABS(VF-V)**(2*VF)/(V-1)**2
    
    real(kind=dp),intent(in) :: xin,alpha,vf
    real(kind=dp),intent(out) :: v
    real(kind=dp) :: eps,vH,vL,vm,fh,fm,fL,tmp,tmp1,tmp2
    integer :: i
  
    eps  = 1.0d-12
  
    vh   = one
    vl   = vf
    vm   = (vh+vl)*half
    fh   = exp(two*(one-vf)*xin/alpha)-(vh-one)**2/abs(vf-vh)**(two*vf)
    if(xin.lt.zero) then
      fm   = exp(two*(one-vf)*xin/alpha)-(vm-one)**2/abs(vf-vm)**(two*vf)
    else
      fm   = -exp(-two*(one-vf)*xin/alpha)+abs(vf-vm)**(two*vf)/(vm-one)**2
    endif
    fl   = exp(-two*(one-vf)*xin/alpha)-abs(vf-vl)**(two*vf)/(vl-one)**2
  
    if(fm.lt.zero) then
  
      fl   = fm 
      vl   = vm
  
      do i = 1,200
          vm = (vh+vl)*half
          fm = exp(2.0_wp*(1.0_wp-vf)*xin/alpha)-(vm-1.0_wp)**2/abs(vf-vm)**(2.0_wp*vf)
          tmp = (vm - vf) * (one - vm)
          if(tmp.lt.zero)tmp = -one*tmp
          if(tmp.le.zero)tmp = 1.0d-14
          tmp1 = (one-vm)
          if(tmp1.lt.zero)tmp1=1.0d-14
          tmp2 = (vm-vf)
          if(tmp2.lt.zero)tmp2=1.0d-14
          fm=xin-alpha*(log(tmp)+(vf+one)*log(tmp1/tmp2)/(one-vf))*half
          if(fm.lt.zero)then
            fl = fm
            vl = vm
          else
            fh = fm
            vh = vm
          endif
          if(abs(fm).lt.eps) exit
      enddo
  
    else
  
      fh   = fm 
      vh   = vm
  
      do i = 1,100
          vm = (vh+vl)*half
          fm =-exp(-two*(one-vf)*xin/alpha)+abs(vf-vm)**(two*vf)/(vm-one)**2
          if(fm.lt.zero)then
            fl = fm
            vl = vm
          else
            fh = fm
            vh = vm
          endif
          if(abs(fm).lt.eps) exit
      enddo
  
    endif
    v = vm
  
    return
  end subroutine rhalf

!=====================================================================================================

  subroutine Flux_Divergence(tin, N_S, N_S_2d, N_S_3d, pinv, qmat, dmat, iagrad, jagrad, dagrad, ielem)

    ! This subroutine calculates elementwise the Divergence of the Conservative Flux

    use variables
    use referencevariables
    use nsereferencevariables

    implicit none

    integer ,                     intent(in)  :: N_S, N_S_2d, N_S_3d, ielem
    real(wp),                     intent(in)  :: tin
    real(wp), dimension(N_S),     intent(in)  :: pinv
    real(wp), dimension(N_S,N_S), intent(in)  :: qmat, dmat
    integer,  dimension(:),       intent(in)  :: iagrad
    integer,  dimension(:,:),     intent(in)  :: jagrad
    real(wp), dimension(:,:),     intent(in)  :: dagrad

    ! indices
    integer :: inode, jdir
    integer :: jnode
    integer :: i

    call Flux_Div_Pencil(ielem, N_S, N_S_2d, pinv, qmat, dmat)    !  result in divF

    if (viscous) then
      ! loop over all nodes in element
      do inode = 1, N_S_3d

        ! calculate viscous flux
        fvg (:,:,inode,ielem) = 0.0_wp

        fvg(:,1:ndim,inode,ielem) = Jx_r(inode,ielem) * &
            viscousflux3D( vg(:,inode,ielem), &
            phig(:,:,inode,ielem), &
            r_x(:,:,inode,ielem), &
            nequations, &
            ndim, &
            mut(inode,ielem)) ! (navierstokes)
      end do
      !
      ! calculate divergence of the flux
      ! 

      ! loop over all nodes in the element
      do inode = 1, N_S_3d
        ! loop over all nonzero columns in CSR corresponding to this row
        do i = iagrad(inode), iagrad(inode+1)-1
          ! loop over each direction
          do jdir = 1,ndim
            ! column/node from gradient operator in CSR format in
            ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
            jnode = jagrad(jdir,i)

            divf(:,jdir,inode,ielem) = divf(:,jdir,inode,ielem) - dagrad(jdir,i) * fvg(:,jdir,jnode,ielem)

          end do
        end do
      end do
    endif

    return
  end subroutine Flux_Divergence

  !============================================================================

  subroutine Entropy_Flux_Divergence(tin,ielem)
    ! This subroutine calculates elementwise 
    ! the Divergence of the entropy Flux
    use variables
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: iagrad,jagrad,dagrad
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem
    real(wp), intent(in) :: tin

    ! indices
    integer :: jdir
    integer :: inode, jnode
    integer :: i

    real(wp), dimension(:,:,:), allocatable :: fvg_err
    real(wp), dimension(:,:,:), allocatable :: divfV_err

   
    call Entropy_Inviscid_Flux_Div(ielem)   !   output in divf_S

    if (viscous) then

      allocate(fvg_err(nequations,ndim,nodesperelem))
      allocate(divfV_err(nequations,ndim,nodesperelem))
      ! loop over all nodes in element
      do inode = 1,nodesperelem

        ! calculate viscous flux
        fvg_err(:,:,inode) = 0.0_wp

        fvg_err(:,1:ndim,inode) = Jx_r(inode,ielem) * &
            viscousflux3D( vg(:,inode,ielem), &
            phig_err(:,:,inode,ielem), &
            r_x(:,:,inode,ielem), &
            nequations, &
            ndim, &
            mut(inode,ielem))

      end do

      divfV_err(:,:,:) = 0.0_wp
      ! loop over all nodes in the element
      do inode = 1,nodesperelem

        ! loop over all nonzero columns in CSR corresponding to this row
        do i = iagrad(inode), iagrad(inode+1)-1
          ! loop over each direction
          do jdir = 1,ndim
            ! column/node from gradient operator in CSR format in
            ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
            jnode = jagrad(jdir,i)

            divfV_err(:,jdir,inode) = divfV_err(:,jdir,inode) - dagrad(jdir,i) * fvg_err(:,jdir,jnode)
          end do
        end do

        do jdir = 1,ndim
          divf_S(jdir,inode,ielem) = divf_S(jdir,inode,ielem) - dot_product(wg(:,inode,ielem),divfV_err(:,jdir,inode))
        enddo

      end do

      deallocate(fvg_err,divfV_err)

    endif

    return
  end subroutine Entropy_Flux_Divergence


  !===========================================================================================
  ! sat_penalty - Calculates both inviscid and viscous penalty according to the SAT procedure.
  !===========================================================================================
  
  subroutine SAT_Penalty(tin, ielem, n_S_1d, n_S_2d, pinv )
    
    ! Load modules
    use variables
    use referencevariables
    use nsereferencevariables
    use controlvariables, only: heat_entropy_flow_wall_bc, Riemann_Diss_BC,  &
                                entropy_flux_BC, SAT_type
    use collocationvariables, only: l01, l00, Sfix, elem_props,              &
                                 Restrct_Gau_2_LGL_1d, Prolong_LGL_2_Gau_1d, &
                                 LGL_Fine_2_LGL_Coarse_1d, LGL_Coarse_2_LGL_Fine_1d

    use initcollocation,  only: ExtrpXA2XB_2D_neq, &
                                JacobiP11, Gauss_Legendre_points, element_properties

    use initgrid

    ! Nothing is implicitly defined
    implicit none

    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in)                     :: ielem, n_S_1d, n_S_2d
    real(wp), intent(in)                     :: tin
    real(wp), dimension(n_S_1d),  intent(in) :: pinv

    ! indices
    integer :: inode, jnode, knode, lnode, gnode
    integer :: kelem
    integer :: iface, kface
    integer :: i,j,k

    real(wp), allocatable :: fstar(:), fstarV(:)               ! reconstructed fluxs
    real(wp), allocatable :: fRoeI(:), fLLF(:)
    real(wp), allocatable :: fn(:), fV_Off(:)                  ! Inviscid and viscous fluxes

    real(wp), allocatable, dimension(:,:) :: smat, sinv        ! right and left eigenvector matrices
    real(wp), allocatable, dimension(:)   :: ev, evabs         ! eigenvalues
    real(wp), allocatable, dimension(:)   :: Vav               ! average state

    real(wp) :: nx(3), nx_On(3), nx_Off(3)                     ! normal vector
    real(wp) :: evmax                                          ! Lax-Freidrich max Eigenvalue
    real(wp) :: Jx_r_Ave                                       ! Average of Jacobian between and On- and Off- elements at interface

    integer                               :: n_S_1d_On , n_S_1d_Off, n_S_1d_Mort, n_S_1d_max
    integer                               :: n_S_2d_On , n_S_2d_Off, n_S_2d_Mort, n_S_2d_max
    integer                               :: poly_val
    integer                               :: nghst_volume, nghst_shell, nghst_Mortar
 
    integer,  allocatable, dimension(:,:) :: Eye, Jay

    real(wp), allocatable, dimension(:)   :: x_S_1d_On, x_S_1d_Off
    real(wp), allocatable, dimension(:)   :: x_S_1d_Mort, w_S_1d_Mort

    real(wp), allocatable, dimension(:)   :: mut_2d_Off
    real(wp), allocatable, dimension(:)   :: Jx_r_2d_Mort

    real(wp), allocatable, dimension(:,:) :: Intrp_On_2_Off_x1, Intrp_On_2_Off_x2
    real(wp), allocatable, dimension(:,:) :: Intrp_Off_2_On_x1, Intrp_Off_2_On_x2

    real(wp), allocatable, dimension(:,:) :: Extrp_Off, Extrp_On
    real(wp), allocatable, dimension(:,:) ::            Intrp_On
    real(wp), allocatable, dimension(:,:) ::            Intrp_Off
    real(wp), allocatable, dimension(:,:) ::            IOn2Off, IOff2On
    real(wp), allocatable, dimension(:,:) :: wg_Mort_On,wg_Mort_Off
    real(wp), allocatable, dimension(:,:) :: vg_2d_On,  vg_2d_Off
    real(wp), allocatable, dimension(:,:) ::            nx_2d_Off
    real(wp), allocatable, dimension(:,:) :: wg_2d_On,  wg_2d_Off
    real(wp), allocatable, dimension(:,:) :: nx_Off_ghst

    integer,  allocatable, dimension(:,:) :: kfacenodes_On, kfacenodes_Off
    integer,  allocatable, dimension(:)   :: ifacenodes_On, ifacenodes_Off
    integer,  allocatable, dimension(:)   :: cnt_Mort_Off

    real(wp), allocatable, dimension(:,:,:) :: phig_2d_On, phig_2d_Off

    logical,  parameter :: testing                = .false.
    logical             :: nonconforming_element  = .false.

    real(wp), parameter :: mirror          = -1.0_wp
    real(wp), parameter :: no_slip         = -0.0_wp


    real(wp), parameter :: Cevmax          =  1.0_wp
    real(wp), parameter :: deltaU          =  0.1_wp
    real(wp), parameter :: stab_strength   =  1.1_wp
    real(wp), parameter :: bc_pen_strength =  1.0_wp
    real(wp), parameter :: LocalLaxF_factor=  2.0_wp

    real(wp), dimension(nequations,nequations) :: matrix_ip

    real(wp), dimension(nequations)       ::   ug_On,   ug_Off
    real(wp), dimension(nequations)       ::   vg_On,   vg_Off
    real(wp), dimension(nequations)       ::   wg_On,   wg_Off

    real(wp), dimension(nequations,3)     :: phig_On, phig_Off
    real(wp), dimension(nequations)       :: SAT_Pen

    real(wp), dimension(nequations)       :: vg_Off_NoSlip
    real(wp), dimension(nequations)       :: wg_Off_NoSlip
    real(wp), dimension(nequations)       :: prim_ref
    
    real(wp), dimension(nequations,3)     :: phig_On_Normal, phig_On_Tangent

    real(wp), dimension(3)                :: unit_normal, normal_vel, tangent_vel, ref_vel_vector

    integer                               :: i_err
 
    ! allocate local arrays
    allocate(fRoeI(nequations))
    allocate(fLLF(nequations))
    allocate(fstar(nequations))
    allocate(fstarV(nequations))
    allocate(fn(nequations))
    allocate(fV_Off(nequations))

    allocate(ev(nequations))
    allocate(evabs(nequations))
    allocate(vav(nequations))
    allocate(smat(nequations,nequations))
    allocate(sinv(nequations,nequations))

    ! initialize penalty
    gsat(:,:,ielem) = 0.0_wp

    n_S_1d_max  = (npoly_max+1)**1
    n_S_2d_max  = (npoly_max+1)**2

    ! establish ``on-element'' face information
    call element_properties(ielem,&
                   n_pts_1d=n_S_1d_On,&
                   n_pts_2d=n_S_2d_On,&
                   x_pts_1d=x_S_1d_On,&
                 kfacenodes=kfacenodes_On,&
                 ifacenodes=ifacenodes_On)

    nghst_volume = nelem_ghst(1,ielem)                  ! point at beginning of ``ielem'' data in ghost stack 
    nghst_shell  = nelem_ghst(2,ielem)                  ! point at beginning of ``ielem'' data in ghost stack 
    nghst_Mortar = nelem_ghst(3,ielem)                  ! point at beginning of ``ielem'' data in ghost stack 
       
    nonconforming_element = .false.
    do iface = 1,nfacesperelem
       if (ef2e(4,iface,ielem) /= elem_props(2,ielem)) then
          nonconforming_element = .true.
          exit
       endif
    enddo

    ! Loop over each face
    faceloop:do iface = 1,nfacesperelem

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       Boundary Contribution to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      if (ef2e(1,iface,ielem) < 0) then


        if (abs(ef2e(1,iface,ielem)) /= 6) then
 
          call set_boundary_conditions(ef2e(1,iface,ielem),InitialCondition)   ! Specify the Boundary Condition procedure on the face

          do i = 1, n_S_2d_On                                                  ! Loop over On_element face nodes

            inode = kfacenodes_On(i,iface)                                     ! Volumetric node index corresponding to face and node on face indices
          
            jnode = (iface-1)* n_S_2d_On + i                                   ! Facial index corresponding to face and node on face indices
          
            nx_On = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)            ! On_Element outward facing face node normal
            if(SAT_type == "mod_metric")then
              nx_Off = Jx_facenodenormal_LGL(:,jnode,ielem)
            elseif(SAT_type == "mod_SAT")then
              nx_Off = nx_On
            else
              write(*,*)'In navierstokes: SAT_Penalty you have chosen an incorrect value of SAT_type = ',&
                SAT_type,' ending computation'
                call PetscFinalize(i_err); stop
            endif
            
              ug_On(:)   =  ug(:,inode,ielem)                                  ! Boundary state on-element
            phig_On(:,:) = phig(:,:,inode,ielem)                               ! dw/dx_j  on-element

            call conserved_to_primitive(ug_On, vg_On, nequations)              ! Rotate into primitive
            call primitive_to_entropy  (vg_On, wg_On, nequations)              ! Rotate into entropy

              vg_Off(:)   =   vg_On(:)                                         !  Preload Off registers in case BC does not overwrite the data
            phig_Off(:,:) = phig_On(:,:)
            
            call BoundaryCondition(vg_Off,phig_Off,fV_Off,nx_Off,xg(:,inode,ielem),tin, &   ! Boundary condition populates vg_Off and fV_Off
              & nequations,ndim,mut(inode,ielem))

            call primitive_to_conserved(vg_Off, ug_Off, nequations)            ! Rotate into conserved
            call primitive_to_entropy  (vg_Off, wg_Off, nequations)            ! Rotate into entropy

                                                                               ! ==  Eigen values/vectors
            call roeavg( vg_On, vg_Off, Vav, nequations )         

            call CharacteristicDecomp( vav, nequations, sinv, smat, ev, nx_On ) 

            evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)
                                                                               ! ==  Eigen values/vectors

                                                                               ! ==  Fluxes
            fn(:) = normalflux( vg_On, nx_On, nequations )                     ! (Euler Flux)

            select case(Riemann_Diss_BC)
              case('LocalLaxF')
                fLLF  = half * ( normalflux( vg_On , nx_On , nequations )  &
                    &        +   normalflux( vg_Off, nx_Off, nequations )  &
                             +   LocalLaxF_factor*evmax*(ug_On - ug_Off) )
                fstar = fLLF
              case('Roe')
                select case(entropy_flux_BC)
                  case('Ismail_Roe'   ) 
!                    fstar = EntropyConsistentFlux     (vg_On, vg_Off, nx_On , nequations ) ! (Entropy Flux)
                     fstar = EntropyConsistentFlux     (vg_On, vg_Off, nx_Off, nequations ) ! (Entropy Flux)
                  case('Chandrashekar') 
!                    fstar = Entropy_KE_Consistent_Flux(vg_On, vg_Off, nx_On , nequations ) ! (Entropy Flux)
                     fstar = Entropy_KE_Consistent_Flux(vg_On, vg_Off, nx_Off, nequations ) ! (Entropy Flux)
                end select
                if(SAT_type == "mod_metric")then
                  fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg_On - wg_Off) )
                elseif(SAT_type == "mod_SAT")then
                  fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg_On - wg_Off) )
                else
                  write(*,*)'navierstokes: SAT_Penalty you have chosen an incorrect value of SAT_type = ',SAT_type
                  call PetscFinalize(i_err); stop
                endif
            end select

            fstarV(:) = normalviscousflux( vg_On, phig_On, nx_On, nequations,mut(inode,ielem)) - fV_Off(:)

            matrix_ip = 0.5_wp * pinv(1) * ( matrix_hatc_node(vg_On ,nx_On,nx_On,nequations)    &  ! IP penalty contribution, i.e. M (u-v), where M is S.P.D. and defined as
                                           + matrix_hatc_node(vg_Off,nx_On,nx_On,nequations)) / &  ! M = pinv(1) (c_ii_side_1 + c_ii_side_2)/2, in the normal direction
                                             Jx_r(inode,ielem)


            gsat(:,inode,ielem) = gsat(:,inode,ielem) + &
              & pinv(1)* ( (fn - fstar) + l01*fstarV*bc_pen_strength ) - &
              & pinv(1)* l00 * matmul(matrix_ip,wg_On - wg_Off)

          end do
        
        else         ! Entropy Stable Solid Wall BC

          do i = 1, n_S_2d_On                                                       ! Loop over each node on the face

            inode = kfacenodes_On(i,iface)                                          ! Volumetric node index corresponding to facenode
          
            jnode = (iface-1)* n_S_2d_On + i                                        ! Facial index corresponding to face and node on face indices

!           nx_On = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)                 ! On_Element outward facing face node normal
!           if(SAT_type == "mod_metric")then
!             nx_Off = Jx_facenodenormal_LGL(:,jnode,ielem)
!           elseif(SAT_type == "mod_SAT")then
!             nx_Off = nx_On
!           else
!             write(*,*)'In navierstokes: SAT_Penalty you have chosen an incorrect value of SAT_type = ',&
!               SAT_type,' ending computation'
!               call PetscFinalize(i_err); stop
!           endif
          
            nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)                    ! Outward facing normal of facial node

            unit_normal = nx/Magnitude(nx)                                          ! Unit normal direction

              ug_On(:)   =   ug(:,  inode,ielem)                                    ! Compute the boundary state starting with the conserved variables
            phig_On(:,:) = phig(:,:,inode,ielem)                                    ! Grad of entropy variables on-element

            call conserved_to_primitive(ug_On, vg_On, nequations)                   ! Rotate into primitives
            call primitive_to_entropy  (vg_On, wg_On, nequations)                   ! Rotate into entropies

            normal_vel = dot_product(vg_On(2:4),unit_normal)*unit_normal            ! Normal velocity

            tangent_vel = vg_On(2:4) - normal_vel                                   ! Tangent velocity of interior velocity field

            vg_Off(1) = vg_On(1)                                                    ! primitive variables with flipped velocity vector
            vg_Off(2) = tangent_vel(1) - normal_vel(1)
            vg_Off(3) = tangent_vel(2) - normal_vel(2)
            vg_Off(4) = tangent_vel(3) - normal_vel(3)
            vg_Off(5) = vg_On(5)                                                    ! Thermodynamic variables are identical (Density and temperature)

            call primitive_to_conserved(vg_Off,ug_Off,nequations)                   ! Off element conserved variables for imposing nonpenetration condition
            call primitive_to_entropy  (vg_Off,wg_Off,nequations)                   ! Off element entropy   variables for imposing nonpenetration condition

                                                                                                    ! ==  Eigen values/vectors
            call roeavg(vg_On,vg_Off,vav,nequations)                                ! Compute the roe average state of the primitive variables

            call CharacteristicDecomp(vav,nequations,sinv,smat,ev,nx)               ! Compute characteristic decomposition

            evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)  ! ==  Eigen values/vectors

            fn = normalflux( vg_On, nx, nequations )                                ! (Euler Flux)

!           select case(Riemann_Diss_BC)
!             case('LocalLaxF')
                fLLF  = half * ( normalflux( vg_On , nx, nequations )  &
                    &        +   normalflux( vg_Off, nx, nequations )  &
                             +   LocalLaxF_factor*evmax*(ug_On - ug_Off) )
                fstar = fLLF
!             case('Roe')
!               fstar = EntropyConsistentFlux(vg(:,inode,ielem), vg_Off, nx, nequations ) ! (Entropy Flux)
!!              fstar = Entropy_KE_Consistent_Flux(vg(:,inode,ielem), vg_Off, nx, nequations ) ! (Entropy Flux)
!               fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg_On - wg_Off) )
!           end select

            vg_Off_NoSlip(1) = vg_On(1)                                             ! no-slip primitive variables 4 no-slip isothermal wall BC
            vg_Off_NoSlip(2) = 0.0_wp
            vg_Off_NoSlip(3) = 0.0_wp
            vg_Off_NoSlip(4) = 0.0_wp
            vg_Off_NoSlip(5) = vg_On(5)

            call primitive_to_entropy(vg_Off_NoSlip,wg_Off_NoSlip,nequations)       ! Rotate into entropy variables for no-slip isothermal wall BC

            phig_Off(1,:) = phig_On(1,:)                                            ! Gradients of the primitive variables
            phig_Off(2,:) = phig_On(2,:)
            phig_Off(3,:) = phig_On(3,:)
            phig_Off(4,:) = phig_On(4,:)

            phig_On_Normal(5,:) = dot_product(phig_On(5,:),unit_normal(:))*unit_normal(:)            ! Normal component of grad_entropy_variables

            phig_On_Tangent(5,:) = phig_On(5,:) - phig_On_Normal(5,:)                                !  Tangent component of grad_entropy_variables
  
            phig_Off(5,:) = phig_On_Tangent(5,:) + heat_entropy_flow_wall_bc*unit_normal(:)/vg_On(5) !  add thermal condition to phig_off

            fstarV(:) = normalviscousflux(vg_On,phig_On ,nx,nequations,mut(inode,ielem))             !  normal viscous flux on-element

            fV_Off(:) = normalviscousflux(vg_On,phig_Off,nx,nequations,mut(inode,ielem))             ! normal viscous flux arising from the ghost point
            
            ref_vel_vector(:) = U0                                                                   ! c_ii reference state matrix
            prim_ref(1) = rho0
            prim_ref(2) = Magnitude(ref_vel_vector)
            prim_ref(3) = Magnitude(ref_vel_vector)
            prim_ref(4) = Magnitude(ref_vel_vector)
            prim_ref(5) = T0

            matrix_ip = pinv(1) * matrix_hatc_node(prim_ref,nx,nx,nequations) / Jx_r(inode,ielem)    ! IP penalty contribution, i.e. M (u-v), where M is S.P.D. and defined as                                           

            matrix_ip(5,:) = 0.0_wp ; matrix_ip(:,5) = 0.0_wp ;                                      ! Zero out last row and column

            gsat(:,inode,ielem) = gsat(:,inode,ielem) &                                              ! Compute the penalty term
              & + pinv(1)* (fn - fstar) &
              & - pinv(1)* (fstarV-fV_Off) &
              & - pinv(1)* l00 * matmul(matrix_ip,wg_On - wg_Off_NoSlip)

          end do
        end if

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       Conforming interface:  polynomial orders match and h conforming 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!      else if (elem_props(2,ielem) == ef2e(4,iface,ielem)) then
      else if ((elem_props(2,ielem) == ef2e(4,iface,ielem)) .and. (ef2e(9,iface,ielem) == 0)) then

        if (ef2e(3,iface,ielem) /= myprocid) then 
        
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !       Off Processor Contributions to gsat:  Conforming Interface straddles parallel partition
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          kelem = ef2e(2,iface,ielem)                                                  ! adjoining element
          kface = ef2e(1,iface,ielem)                                                  ! face on element

          do i = 1,  n_S_2d_On

            jnode =  n_S_2d_On*(iface-1) + i                                           ! Index in facial ordering
            
            inode = ifacenodes_On(jnode)                                               ! Volumetric node index corresponding to facial node index
            
            gnode = efn2efn(3,jnode,ielem)                                             ! Index in Petsc ghost stack (not volumetric stack)
  
            ug_On(:)  = ug   (:,inode,ielem)                                           ! On -element solution vector
            ug_Off(:) = ughst(:,gnode)                                                 ! Off-element solution vector from PETSc ghost registers

            call conserved_to_primitive(ug_On (:), vg_On (:), nequations)              ! Rotate to conserved variables
            call conserved_to_primitive(ug_Off(:), vg_Off(:), nequations)
  
            phig_On (:,:) = phig   (:,:,inode,ielem)                                   ! On -element dW/dx_j
            phig_Off(:,:) = phighst(:,:,gnode)                                         ! Off-element dW/dx_j from PETSC ghost registers

            Jx_r_Ave = 0.5_wp * (abs(Jx_r(inode,ielem)) + abs(Jx_r_ghst_LGL(gnode)) )  ! Average of Jacobian at the face

            nx_On = + Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)                  ! Outward facing normal of facial node
            if(nonconforming_element) then
              nx_Off = - nxghst_LGL_Shell(:,nghst_shell + i)                           ! Outward facing normal in Petsc ghost registers
            else
              nx_Off = + Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)               ! Outward facing normal of facial node
            endif

            SAT_Pen(:) =  SAT_Inv_Vis_Flux( nequations,iface,ielem,    &               ! Build SAT for inviscid and viscous contributions
                                          & vg_On,vg_Off,              &
                                          & phig_On,phig_Off,          &
                                          & nx_On,nx_Off,Jx_r_Ave,     &
                                          & pinv(1), mut(inode,ielem))

            gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * SAT_Pen(:)           !  Update global SAT

          end do

          nghst_volume = nghst_volume + n_S_2d_On                                      !  Track position in Ghost stack (n_S_2d_On=n_S_2d_Off)
          nghst_shell  = nghst_shell  + n_S_2d_On                                      !  Track position in Ghost shell stack (On and Off are the same)

        else if (ef2e(3,iface,ielem) == myprocid) then 

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !       On-Processor Contributions to gsat:  Conforming Interface is on-process
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          kelem = ef2e(2,iface,ielem)                                                  ! adjoining element
          kface = ef2e(1,iface,ielem)                                                  ! face on element

          do i = 1, n_S_2d_On
        
            jnode =  n_S_2d_On*(iface-1) + i                                           ! Index in facial ordering
              
            inode = ifacenodes_On(jnode)                                               ! Volumetric node index corresponding to facial node index
              
            knode = efn2efn(1,jnode,ielem)                                             ! Volumetric index of partner node
              
            vg_On(:)  = vg(:,inode,ielem)                                              ! On -element primitive variables
            vg_Off(:) = vg(:,knode,kelem)                                              ! Off-element primitive variables

            phig_On (:,:) = phig(:,:,inode,ielem)                                      ! On -element dW/dx_j
            phig_Off(:,:) = phig(:,:,knode,kelem)                                      ! Off-element dW/dx_j from PETSC ghost registers

            Jx_r_Ave = 0.5_wp * (abs(Jx_r(inode,ielem)) + abs(Jx_r(knode,kelem)))      ! Average of Jacobian at the face
 
            nx_On  = + Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)                 ! Outward facing normal of facial node (On -Element)
            if(nonconforming_element) then
              nx_Off = + Jx_facenodenormal_LGL(:,jnode,ielem)                          ! Initially conforming interfaces on nonconforming element were the same. Just use self.
            else
              nx_Off = + Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)               ! Outward facing normal of facial node (Off-Element)
            endif

            SAT_Pen(:) =  SAT_Inv_Vis_Flux( nequations,iface,ielem,    &               ! Build SAT for inviscid and viscous contributions
                                          & vg_On,vg_Off,              &
                                          & phig_On,phig_Off,          &
                                          & nx_On,nx_Off,Jx_r_Ave,     &
                                          & pinv(1), mut(inode,ielem))

            gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * SAT_Pen(:)           ! Update global SAT

          end do
        endif

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       NON-CONFORMING p refinement Contributions to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!      else if (elem_props(2,ielem) /= ef2e(4,iface,ielem)) then
      else if ((elem_props(2,ielem) /= ef2e(4,iface,ielem)) .and. (ef2e(9,iface,ielem) == 0)) then

        kface       = ef2e(1,iface,ielem)
        kelem       = ef2e(2,iface,ielem)

        call element_properties(kelem,&
                       n_pts_1d=n_S_1d_Off,&
                       n_pts_2d=n_S_2d_Off,&
                       x_pts_1d=x_S_1d_Off,&
                     kfacenodes=kfacenodes_Off,&
                     ifacenodes=ifacenodes_Off)

        n_S_1d_Mort = max(n_S_1d_On,n_S_1d_Off)
        n_S_2d_Mort = (n_S_1d_Mort)**2
        if(allocated(x_S_1d_Mort)) deallocate(x_S_1d_Mort) ; allocate(x_S_1d_Mort(n_S_1d_Mort)) ;
        if(allocated(w_S_1d_Mort)) deallocate(w_S_1d_Mort) ; allocate(w_S_1d_Mort(n_S_1d_Mort)) ;
        call Gauss_Legendre_points(n_S_1d_Mort,x_S_1d_Mort,w_S_1d_Mort)

        allocate(  nx_2d_Off (ndim,n_S_2D_Off ))
        allocate( mut_2d_Off (     n_S_2D_Off ))
        allocate(Jx_r_2d_Mort(     n_S_2D_Mort))

        allocate(vg_2d_On    (nequations,n_S_2d_On  ))
        allocate(wg_2d_On    (nequations,n_S_2d_On  ))

        allocate(vg_2d_Off   (nequations,n_S_2d_Off ))
        allocate(wg_2d_Off   (nequations,n_S_2d_Off ))

        allocate( wg_Mort_On (nequations,n_S_2d_Mort))
        allocate( wg_Mort_Off(nequations,n_S_2d_Mort))

        allocate(cnt_Mort_Off(           n_S_2d_Mort))

        allocate(phig_2d_On  (nequations,ndim,n_S_2d_On ))
        allocate(phig_2d_Off (nequations,ndim,n_S_2d_Off))

        if(n_S_1d_Mort == n_S_1d_On) then
          poly_val = n_S_1d_Mort - npoly
           allocate(Intrp_On (n_S_1d_On  ,n_S_1d_Mort)) ; 
                    Intrp_On (:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_On  ,1:n_S_1d_Mort,poly_val,1) ;
           allocate(Extrp_On (n_S_1d_Mort,n_S_1d_On  )) ; 
                    Extrp_On (:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_On  ,poly_val,1) ;
          poly_val = n_S_1d_Off  - npoly
           allocate(Extrp_Off(n_S_1d_Mort,n_S_1d_Off )) ; 
                    Extrp_Off(:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_Off ,poly_val,2) ;
           if(allocated(Intrp_Off)) deallocate(Intrp_Off)   
             allocate(Intrp_Off (n_S_1d_Off  ,n_S_1d_Mort)) ; 
             Intrp_Off (:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_Off  ,1:n_S_1d_Mort,poly_val,2) ;
      
                  
        else
          poly_val = n_S_1d_On - npoly
           allocate(Intrp_On (n_S_1d_On  ,n_S_1d_Mort)) ; 
                    Intrp_On (:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_On  ,1:n_S_1d_Mort,poly_val,2) ;
           allocate(Extrp_On (n_S_1d_Mort,n_S_1d_On  )) ; 
                    Extrp_On (:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_On  ,poly_val,2) ;
          poly_val = n_S_1d_Mort - npoly
           allocate(Extrp_Off(n_S_1d_Mort,n_S_1d_Off )) ; 
                    Extrp_Off(:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_Off ,poly_val,1) ;
           if(allocated(Intrp_Off)) deallocate(Intrp_Off)
             allocate(Intrp_Off (n_S_1d_Off  ,n_S_1d_Mort)) ; 
             Intrp_Off (:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_Off  ,1:n_S_1d_Mort,poly_val,1) ;

        endif

        !QUESTION: why do I aonly have to do this for these and Intrp_Off
        if(allocated(IOn2Off)) deallocate(IOn2Off)
          allocate(IOn2Off(1:n_S_1d_Off,1:n_S_1d_On)) 
          IOn2Off = matmul(Intrp_Off,Extrp_On)
        if(allocated(IOff2On)) deallocate(IOff2On)
          allocate(IOff2On(1:n_S_1d_On,1:n_S_1d_Off))
          IOff2On = matmul(Intrp_On,Extrp_Off)

!=========
!       face data:  On_element, Off_element and On_Mortar
!=========
                                                                             ! ============================================
                                                                             ! On_Element face data same for serial or parallel
                                                                             ! ============================================
        On_Elem_0:do i = 1, n_S_2d_On                                             ! On_Element Loop over 2D LGL points
        
         jnode =  n_S_2d_On*(iface-1) + i                                         ! Index in facial ordering
         inode = ifacenodes_On(jnode)                                             ! Volumetric node index corresponding to facial node index

           vg_2d_On(:,  i) =   vg(:,  inode,ielem)                                ! On-element face data
         phig_2d_On(:,:,i) = phig(:,:,inode,ielem)                                ! Viscous derivatives

         call primitive_to_entropy(vg_2d_On(:,i),wg_2d_On(:,i),nequations)        ! Rotate into entropy variables and store as face plane data

        enddo  On_Elem_0                                                          ! End off-element loop

                                                                             ! ============================================
        if (ef2e(3,iface,ielem) /= myprocid) then                            ! ====== Parallel NON-CONFORMING data ========
                                                                             ! ============================================
          Off_Elem_0:do k = 1, n_S_2d_Off                                         ! Off-element loop over data

           
                 ug_Off(:    ) =           ughst(:,  nghst_volume + k)            ! conserved variable    in Petsc ghost registers
            phig_2d_Off(:,:,k) =         phighst(:,:,nghst_volume + k)            ! Viscous derivatives   in Petsc ghost registers
             mut_2d_Off(    k) =         mutghst(    nghst_volume + k)            ! Turbulent viscosity   in Petsc ghost registers

              nx_2d_Off(:,  k) =  nxghst_LGL_Shell(:,nghst_shell  + k)            ! Outward facing normal in Petsc ghost registers

            call conserved_to_primitive(ug_Off   (:  ),vg_2d_Off(:,k),nequations) ! Rotate into primitive variables and store as face plane data
            call primitive_to_entropy  (vg_2d_Off(:,k),wg_2d_Off(:,k),nequations) ! Rotate into entropy   variables and store as face plane data

          enddo Off_Elem_0                                                        ! End off-element loop

          nghst_volume = nghst_volume + n_S_2d_Off                                !  Keep track of position in Ghost volume stack
          nghst_shell  = nghst_shell  + n_S_2d_Off                                !  Keep track of position in Ghost shell  stack

          On_Mortar_0:do j = 1, n_S_2d_Mort

            jnode =  n_S_2d_max*(iface-1) + j                                     ! Index in facial ordering

            lnode = efn2efn_Gau(3,jnode,ielem)

            cnt_Mort_Off(j) = efn2efn_Gau(3,jnode,ielem) - nghst_Mortar

            Jx_r_2d_Mort(j) = 0.5_wp * (abs(Jx_r_Gau_shell(jnode,ielem)) + abs(Jx_r_ghst_Gau_shell(lnode)))
     
          enddo On_Mortar_0

          nghst_Mortar = nghst_Mortar + n_S_2d_Mort
                                                                             ! ============================================
        else                                                                 ! ====== Serial NON-CONFORMING data ==========
                                                                             ! ============================================
          Off_Elem_1:do k = 1, n_S_2d_Off                                         ! Off-element loop over data
 
            lnode =  n_S_2d_Off*(kface-1) + k                                     ! Index in facial ordering
            knode = ifacenodes_Off(lnode)                                         ! Volumetric node index corresponding to facial node index

              vg_2d_Off(:,  k) =   vg(:,  knode,kelem)                            ! volumetric node data from off element face
            phig_2d_Off(:,:,k) = phig(:,:,knode,kelem)                            ! Viscous derivatives
             mut_2d_Off(    k) =  mut(    knode,kelem)                            ! Outward facing normal of facial node

              nx_2d_Off(:,  k) = Jx_facenodenormal_LGL(:,lnode,kelem)             ! Outward facing normal of facial node

            call primitive_to_entropy(vg_2d_Off(:,k),wg_2d_Off(:,k),nequations)   ! Rotate into entropy variables and store as face plane data

          enddo Off_Elem_1                                                        ! End off-element loop

          On_Mortar_1:do j = 1, n_S_2d_Mort

            jnode =  n_S_2d_max*(iface-1) + j                                     ! Index in facial ordering

            cnt_Mort_Off(j) = efn2efn_Gau(4,jnode,ielem) - n_S_2d_max*(kface-1)   ! Correct for face orientation and shift back to 1:n_S_2d_Mort

            lnode = efn2efn_Gau(4,jnode,ielem)

            Jx_r_2d_Mort(j) = 0.5_wp * (abs(Jx_r_Gau_shell(jnode,ielem)) + abs(Jx_r_Gau_shell(lnode,kelem)))
     
          enddo On_Mortar_1

        endif

        ! Extrapolate On_element and Off_element entropy variables wg ==>  Mortar
        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_On ,n_S_1d_Mort,x_S_1d_On ,x_S_1d_Mort, &
                               wg_2d_On ,wg_Mort_On ,Extrp_On ,Extrp_On )
        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                               wg_2d_Off,wg_Mort_Off,Extrp_Off,Extrp_Off)

!=========
!       Inviscid interface SATs (skew-symmetric portion + Upwind Entropy Stable dissipation)
!=========

        if(SAT_type == "mod_metric")then                                         !-- modified metric approach

          call Inviscid_SAT_Non_Conforming_Interface_Mod_Metric(ielem, iface, kface, ifacenodes_On ,n_S_2d_max, &
                                                     n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                     n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                     n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,            &
                                                     pinv,                                           &
                                                     vg_2d_On,  vg_2d_Off, wg_Mort_On, wg_Mort_Off,  &
                                                     cnt_Mort_Off, Intrp_On, Extrp_Off)

        elseif(SAT_type == "mod_SAT")then                                        !-- modified SAT approach

          if(allocated(nx_Off_ghst)) deallocate(nx_Off_ghst)
          allocate(nx_Off_ghst(3,n_S_2d_Off))

          nx_Off_ghst = -nx_2d_Off

          if(.true.)then
            call Inviscid_SAT_Non_Conforming_Interface_Mod_SAT(ielem, iface, ifacenodes_On ,n_S_2d_max, &
                                                        n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                        n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                        pinv, nx_Off_ghst,                              &
                                                        vg_2d_On,  vg_2d_Off, wg_2d_On, wg_2d_Off,      &
                                                        IOn2Off, IOn2Off,                               &
                                                        IOff2On, IOff2On )

          else
            !-- old version
            call Inviscid_SAT_Non_Conforming_Interface_Mod_SAT_OLD(ielem, iface, ifacenodes_On ,n_S_2d_max, &
                                                        n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                        n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                        n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,            &
                                                        pinv,                                           &
                                                        vg_2d_On,  vg_2d_Off, wg_Mort_On, wg_Mort_Off,  &
                                                        cnt_Mort_Off, Intrp_On, Extrp_Off, nx_Off_ghst)
          endif
        else
          write(*,*)'In navierstokes: SAT_Penalty you have chosen an incorrect value of SAT_type = ',&
            SAT_type,' ending computation'
          call PetscFinalize(i_err); stop
        endif

! ========
!       Viscous interface SATs 
! ========

        if(viscous .eqv. .true.) then
          call Viscous_SAT_Non_Conforming_Interface(ielem, kelem, iface, kface,                    &
                                                    ifacenodes_On, ifacenodes_Off, n_S_2d_max,     &
                                                    n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,           &
                                                    n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,           &
                                                    n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,           &
                                                    pinv,                                          &
                                                      vg_2d_On,   vg_2d_Off,                       &
                                                    phig_2d_On, phig_2d_Off,                       &
                                                    wg_Mort_On, wg_Mort_Off,                       &
                                                    nx_2d_Off, mut_2d_Off, Jx_r_2d_Mort,           &
                                                    cnt_Mort_Off, Intrp_On, Extrp_Off)
        endif

        deallocate(vg_2d_On,  vg_2d_Off  )
        deallocate(wg_2d_On,  wg_2d_Off  )
        deallocate(wg_Mort_On,wg_Mort_Off,cnt_Mort_Off)
        deallocate(nx_2d_Off,mut_2d_Off  )
        deallocate(Extrp_Off,Extrp_On)
        deallocate(Intrp_On)
        deallocate(x_S_1d_Off,x_S_1d_Mort)
        deallocate(phig_2d_On, phig_2d_Off)
        deallocate(Jx_r_2d_Mort)

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       NON-CONFORMING h refinement: note that ghost faces are skipped over since ef2e(9,iface,ielem) = -1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      elseif(ef2e(9,iface,ielem) == 1)then

        kface       = ef2e(1,iface,ielem)
        kelem       = ef2e(2,iface,ielem)

        call element_properties(kelem,&
                       n_pts_1d=n_S_1d_Off,&
                       n_pts_2d=n_S_2d_Off,&
                       x_pts_1d=x_S_1d_Off,&
                     kfacenodes=kfacenodes_Off,&
                     ifacenodes=ifacenodes_Off)

        allocate(  nx_2d_Off (ndim,n_S_2D_Off ))
        allocate( mut_2d_Off (     n_S_2D_Off ))

        allocate(vg_2d_On    (nequations,n_S_2d_On  ))
        allocate(wg_2d_On    (nequations,n_S_2d_On  ))

        allocate(vg_2d_Off   (nequations,n_S_2d_Off ))
        allocate(wg_2d_Off   (nequations,n_S_2d_Off ))

        allocate(phig_2d_On  (nequations,ndim,n_S_2d_On ))
        allocate(phig_2d_Off (nequations,ndim,n_S_2d_Off))

!         Sub-face ordering

!                face 1      face 2      face 3      face 4      face 5       face 6
!               ---------   ---------   ---------   ---------   ---------   ---------
!               | 19| 13|   | 20| 14|   | 21| 15|   | 22| 16|   | 23| 17|   | 24| 18|  
!               ---------   ---------   ---------   ---------   ---------   ---------
!               |  1|  7|   |  2|  8|   |  3|  9|   |  4| 10|   |  5| 11|   |  6| 12|
!               ---------   ---------   ---------   ---------   ---------   ---------
!               ---------   ---------   ---------   ---------   ---------   ---------
!               | +-| --|   | 20| 14|   | 21| 15|   | 22| 16|   | 23| 17|   | 24| 18|  
!               ---------   ---------   ---------   ---------   ---------   ---------
!               | ++| -+|   |  2|  8|   |  3|  9|   |  4| 10|   |  5| 11|   |  6| 12|
!               ---------   ---------   ---------   ---------   ---------   ---------
!
        allocate(Eye(n_S_1d_On,n_S_1d_On)); Eye(:,:) = 0 ;                                  !  Initialize matrix to zero
        allocate(Jay(n_S_1d_On,n_S_1d_On)); Jay(:,:) = 0 ;                                  !  Initialize matrix to zero

        do i = 1,n_S_1d_On ; Eye(i,i) = 1 ; Jay(i,n_S_1d_On+1-i) = 1 ; enddo                !  Set diagonal and per-diagonal to 1

        allocate(Intrp_On_2_Off_x1(n_S_1d_On,n_S_1d_On)) ; Intrp_On_2_Off_x1(:,:) = 0.0_wp  ! Intialize Interpolation matrices
        allocate(Intrp_Off_2_On_x1(n_S_1d_On,n_S_1d_On)) ; Intrp_Off_2_On_x1(:,:) = 0.0_wp  ! Intialize Interpolation matrices
        allocate(Intrp_On_2_Off_x2(n_S_1d_On,n_S_1d_On)) ; Intrp_On_2_Off_x2(:,:) = 0.0_wp  ! Intialize Interpolation matrices
        allocate(Intrp_Off_2_On_x2(n_S_1d_On,n_S_1d_On)) ; Intrp_Off_2_On_x2(:,:) = 0.0_wp  ! Intialize Interpolation matrices

        if(ef2e(10,iface,ielem) == 0) then                                                  !  The Coarse Side
              
           select case(iface)
             case( 1: 6)
                Intrp_On_2_Off_x1(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
             case( 7:12)
                Intrp_On_2_Off_x1(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
             case(13:18)
                Intrp_On_2_Off_x1(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
             case(19:24)
                Intrp_On_2_Off_x1(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(eye,LGL_Coarse_2_LGL_Fine_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(Jay,LGL_Coarse_2_LGL_Fine_1d(:,:))
           end select

        elseif(ef2e(10,iface,ielem) == 1) then                                              !  The Fine side

           select case(iface)
             case( 1: 6)
                Intrp_On_2_Off_x1(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
             case( 7:12)
                Intrp_On_2_Off_x1(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
             case(13:18)
                Intrp_On_2_Off_x1(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
             case(19:24)
                Intrp_On_2_Off_x1(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_On_2_Off_x2(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x1(:,:) = matmul(eye,LGL_Fine_2_LGL_Coarse_1d(:,:))
                Intrp_Off_2_On_x2(:,:) = matmul(Jay,LGL_Fine_2_LGL_Coarse_1d(:,:))
           end select

        endif

!=========
!       face data:  On_element, Off_element
!=========
                                                                             ! ============================================
                                                                             ! On_Element face data same for serial or parallel
                                                                             ! ============================================
        On_Elem_1:do i = 1, n_S_2d_On                                             ! On_Element Loop over 2D LGL points
        
         jnode =  n_S_2d_On*(iface-1) + i                                         ! Index in facial ordering
         inode = ifacenodes_On(jnode)                                             ! Volumetric node index corresponding to facial node index

           vg_2d_On(:,  i) =   vg(:,  inode,ielem)                                ! On-element face data
         phig_2d_On(:,:,i) = phig(:,:,inode,ielem)                                ! Viscous derivatives

         call primitive_to_entropy(vg_2d_On(:,i),wg_2d_On(:,i),nequations)        ! Rotate into entropy variables and store as face plane data

        enddo  On_Elem_1                                                          ! End off-element loop

                                                                             ! ============================================
        if (ef2e(3,iface,ielem) /= myprocid) then                            ! ====== Parallel NON-CONFORMING data ========
                                                                             ! ============================================
          Off_Elem_2:do k = 1, n_S_2d_Off                                         ! Off-element loop over data

           
                 ug_Off(:    ) =           ughst(:,  nghst_volume + k)            ! conserved variable    in Petsc ghost registers
            phig_2d_Off(:,:,k) =         phighst(:,:,nghst_volume + k)            ! Viscous derivatives   in Petsc ghost registers
             mut_2d_Off(    k) =         mutghst(    nghst_volume + k)            ! Turbulent viscosity   in Petsc ghost registers

              nx_2d_Off(:,  k) =  nxghst_LGL_Shell(:,nghst_shell  + k)            ! Outward facing normal in Petsc ghost registers

            call conserved_to_primitive(ug_Off   (:  ),vg_2d_Off(:,k),nequations) ! Rotate into primitive variables and store as face plane data
            call primitive_to_entropy  (vg_2d_Off(:,k),wg_2d_Off(:,k),nequations) ! Rotate into entropy   variables and store as face plane data

          enddo Off_Elem_2                                                        ! End off-element loop

          nghst_volume = nghst_volume + n_S_2d_Off                                !  Keep track of position in Ghost volume stack
          nghst_shell  = nghst_shell  + n_S_2d_Off                                !  Keep track of position in Ghost shell  stack

                                                                             ! ============================================
        else                                                                 ! ====== Serial NON-CONFORMING data ==========
                                                                             ! ============================================
          Off_Elem_3:do k = 1, n_S_2d_Off                                         ! Off-element loop over data
 
            lnode =  n_S_2d_Off*(kface-1) + k                                     ! Index in facial ordering
            knode = ifacenodes_Off(lnode)                                         ! Volumetric node index corresponding to facial node index

              vg_2d_Off(:,  k) =   vg(:,  knode,kelem)                            ! volumetric node data from off element face
            phig_2d_Off(:,:,k) = phig(:,:,knode,kelem)                            ! Viscous derivatives
             mut_2d_Off(    k) =  mut(    knode,kelem)                            ! Outward facing normal of facial node

              nx_2d_Off(:,  k) = Jx_facenodenormal_LGL(:,lnode,kelem)             ! Outward facing normal of facial node

            call primitive_to_entropy(vg_2d_Off(:,k),wg_2d_Off(:,k),nequations)   ! Rotate into entropy variables and store as face plane data

          enddo Off_Elem_3                                                        ! End off-element loop

        endif

!=========
!       Inviscid interface SATs (skew-symmetric portion + Upwind Entropy Stable dissipation)
!=========

        if(SAT_type == "mod_SAT")then                                         !-- modified metric approach

          if(allocated(nx_Off_ghst)) deallocate(nx_Off_ghst)
          allocate(nx_Off_ghst(3,n_S_2d_Off))

          nx_Off_ghst = -nx_2d_Off

          call Inviscid_SAT_Non_Conforming_Interface_Mod_SAT(ielem, iface, ifacenodes_On ,n_S_2d_max, &
                                                      n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                      n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                      pinv, nx_Off_ghst,                              &
                                                      vg_2d_On,  vg_2d_Off, wg_2d_On, wg_2d_Off,      &
                                                      Intrp_On_2_Off_x1, Intrp_On_2_Off_x2,           &
                                                      Intrp_Off_2_On_x1, Intrp_Off_2_On_x2)
        else
          write(*,*)'In navierstokes: SAT_Penalty you have chosen an incorrect value of SAT_type = ',&
            SAT_type,' ending computation'
          call PetscFinalize(i_err); stop
        endif

! ========
!       Viscous interface SATs 
! ========

        if(viscous .eqv. .true.) then
          call Viscous_SAT_Non_Conforming_Interface(ielem, kelem, iface, kface,                    &
                                                    ifacenodes_On, ifacenodes_Off, n_S_2d_max,     &
                                                    n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,           &
                                                    n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,           &
                                                    n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,           &
                                                    pinv,                                          &
                                                      vg_2d_On,   vg_2d_Off,                       &
                                                    phig_2d_On, phig_2d_Off,                       &
                                                    wg_Mort_On, wg_Mort_Off,                       &
                                                    nx_2d_Off, mut_2d_Off, Jx_r_2d_Mort,           &
                                                    cnt_Mort_Off, Intrp_On, Extrp_Off)
        endif

        deallocate(Eye,Jay)
        deallocate(Intrp_On_2_Off_x1,Intrp_On_2_Off_x2)
        deallocate(Intrp_Off_2_On_x1,Intrp_Off_2_On_x2)

      end if                                                        !  Main ``if'' conditional (interface / BC type) in SAT_Pen routine

    end do faceloop                                                 !  Main ``face'' conditional in SAT_Pen routine
    

    deallocate(x_S_1d_On)                                           ! Deallocate memory

    deallocate(fRoeI,fLLF)
    deallocate(fstar,fstarV)
    deallocate(fn,fV_Off)
    deallocate(sinv,smat)
    deallocate(ev,evabs,vav)

    return
  end subroutine SAT_Penalty

  !============================================================================
  
  subroutine Inviscid_SAT_Non_Conforming_Interface_Mod_Metric(ielem, iface, kface, ifacenodes_On ,n_S_2d_max, &
                                                     n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                     n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                     n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,            &
                                                     pinv,                                           &
                                                     vg_2d_On,  vg_2d_Off, wg_Mort_On, wg_Mort_Off,  &
                                                     cnt_Mort_Off, Intrp_On, Extrp_Off)

    use referencevariables,   only: nequations, ndim
    use variables,            only: facenodenormal, Jx_r, Jx_facenodenormal_Gau, gsat
    use initcollocation,      only: ExtrpXA2XB_2D_neq, ExtrpXA2XB_2D_neq_k

    implicit none

    integer,                    intent(in) :: ielem, iface, kface
    integer,                    intent(in) :: n_S_1d_On, n_S_1d_Off, n_S_1d_Mort
    integer,                    intent(in) :: n_S_2d_On, n_S_2d_Off, n_S_2d_Mort, n_S_2d_max
    real(wp),  dimension(:),    intent(in) :: x_S_1d_On, x_S_1d_Off, x_S_1d_Mort
    real(wp),  dimension(:),    intent(in) :: pinv
    integer,   dimension(:),    intent(in) :: ifacenodes_On
    integer,   dimension(:),    intent(in) :: cnt_Mort_Off
    real(wp),  dimension(:,:),  intent(in) :: vg_2d_On,   vg_2d_Off
    real(wp),  dimension(:,:),  intent(in) :: wg_Mort_On, wg_Mort_Off
    real(wp),  dimension(:,:),  intent(in) :: Intrp_On, Extrp_Off
    
    real(wp),  dimension(ndim)             :: nx
    real(wp),  dimension(nequations)       :: fn, fstar
    real(wp),  dimension(nequations)       :: vg_On, vg_Off

    real(wp), allocatable, dimension(:,:) :: FxA, FyA, FzA
    real(wp), allocatable, dimension(:,:) :: FxB, FyB, FzB
    real(wp), allocatable, dimension(:,:) :: FC_Mort_On
    real(wp), allocatable, dimension(:,:) :: Up_diss_Mort, Up_diss_On

    integer                               :: i, j, k, l
    integer                               :: ival, jval
    integer                               :: inode, jnode

    continue

      allocate(FxA(nequations,n_S_2D_Off ), FyA(nequations,n_S_2D_Off ), FzA(nequations,n_S_2D_Off ))
      allocate(FxB(nequations,n_S_2D_Mort), FyB(nequations,n_S_2D_Mort), FzB(nequations,n_S_2D_Mort))
      allocate( FC_Mort_On (nequations,n_S_2D_Mort))
      allocate(Up_diss_Mort(nequations,n_S_2d_Mort))
      allocate(Up_diss_On  (nequations,n_S_2d_On  ))

!=========
!         Skew-symmetric matrix portion of SATs (Entropy Stable through Mortar)
!=========

      On_Element_1:do i = 1, n_S_2d_On                                       ! On_Element Loop over 2D LGL points
        
        Off_Element_1:do k = 1, n_S_2d_Off                                   ! Off_Element Loop over 2D LGL points

          call EntropyConsistentFlux_Vectors(vg_2d_On(:,i), vg_2d_Off(:,k), nequations, FxA(:,k), FyA(:,k), FzA(:,k)) ! (Entropy Flux vectors)

        enddo Off_Element_1                                                  ! End Off_Element

        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                               FxA,FxB,Extrp_Off,Extrp_Off)                  ! Extrapolate f^S_x
        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                               FyA,FyB,Extrp_Off,Extrp_Off)                  ! Extrapolate f^S_y
        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                               FzA,FzB,Extrp_Off,Extrp_Off)                  ! Extrapolate f^S_z

        On_Mortar_1:do j = 1, n_S_2d_Mort                                    ! Mortar loop over 2D Gauss points
  
          jnode =  n_S_2d_max*(iface-1) + j                                  ! Index in Mortar facial ordering
  
          nx(:) = Jx_facenodenormal_Gau(:,jnode,ielem)                       ! Outward facing normal on mortar

              l = cnt_Mort_Off(j)                                            ! Correct for face orientation and shift back to 1:n_S_2d_Mort

          FC_Mort_On(:,j) = FxB(:,l)*nx(1) + FyB(:,l)*nx(2) + FzB(:,l)*nx(3) ! Outward facing component of f^S on mortar

        enddo On_Mortar_1                                                    ! End Mortar loop

        ival = mod(i-1,n_S_1d_On) + 1 ; jval = (i-ival) / n_S_1d_On + 1 ;    ! Decode On-element point coordinates (i,j) from planar coordinates

        call ExtrpXA2XB_2D_neq_k(nequations,n_S_1d_Mort,n_S_1d_On,ival,jval,x_S_1d_Mort,x_S_1d_On, &
                                 FC_Mort_On,fstar,Intrp_On,Intrp_On)         ! Restrict planar data to the (ival,jval) point

        jnode =  n_S_2d_On*(iface-1) + i                                     ! Index in facial ordering
              
        inode = ifacenodes_On(jnode)                                         ! Volumetric node index corresponding to facial node index
              
        nx(:) = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)              ! Outward facing normal of facial node

        fn(:) = normalflux(vg_2d_On(:,i), nx(:), nequations)                 ! One point flux based on vg_On and nx

        gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1)*(fn - fstar)     ! SAT penalty:  subtract the on-element contribution and replace with penalty 

      enddo On_Element_1

!=========
!         Inviscid interface dissipation (Entropy Stable Upwinding of SATs)
!=========

      On_Mortar_2:do j = 1, n_S_2d_Mort                                      ! Mortar loop over data
  
        jnode =  n_S_2d_max*(iface-1) + j                                    ! Index in facial ordering (bucket is padded so n_S_2d_max is needed)
  
        nx(:) = Jx_facenodenormal_Gau(:,jnode,ielem)                         ! Outward facing normal of facial node

            l = cnt_Mort_Off(j)                                              ! Correct for face orientation and shift back to 1:n_S_2d_Mort

       call entropy_to_primitive(wg_Mort_On (:,j),vg_On (:),nequations)      ! Entropy -> primitive variables:  On_element
       call entropy_to_primitive(wg_Mort_Off(:,l),vg_Off(:),nequations)      ! Entropy -> primitive variables: Off_element

       Up_diss_Mort(:,j) = SAT_Vis_Diss(nequations,vg_On(:),vg_Off(:),nx(:)) ! Viscous dissipation on Mortar based on L-R states

      enddo On_Mortar_2                                                      ! End mortar loop

                                                                                 ! Restrict data plane from Mortar to on-face plane
      call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Mort ,n_S_1d_On ,x_S_1d_Mort ,x_S_1d_On, &
                             Up_diss_Mort,Up_diss_On, Intrp_On , Intrp_On)

      On_Elem_2: do i = 1, n_S_2d_On                                         ! On-element loop: Begin
        
        jnode =  n_S_2d_On*(iface-1) + i                                     ! Index in facial ordering
              
        inode = ifacenodes_On(jnode)                                         ! Volumetric node index corresponding to facial node index
              
        gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * Up_diss_On(:,i)! On-element Viscous penalty contribution

      end do On_Elem_2                                                       ! On-element loop: end

      deallocate(FxA,FyA,FzA)
      deallocate(FxB,FyB,FzB)
      deallocate(FC_Mort_On)
      deallocate(Up_diss_On,Up_diss_Mort)

  end subroutine Inviscid_SAT_Non_Conforming_Interface_Mod_Metric
  !============================================================================
  !
  ! Purpose: Constructs the inviscid SAT for the modified SAT approach where the 
  !          symmetric portion of the sat is for the xil computational direction and 
  !          xm Cartesian direction
  !\begin{equation*}
  !  P^{-1}E_{xil,ielem}\circ\[Jdxildx]_{ielem}\circ F^{SC}(q_ielem,q_ielem)1_{ielem}
  !  -\frac{1}{2}P^{-1}\left([Jdxildxm]_{ielem} Evtok + Evtok [Jdxildx_m]_{kelem}\right)\circ 
  !  F^{SC}(q_{ielem},q_{kelem}) 1_{kelem}
  !\end{equation*}
  !
  ! The upwinding contribution is given as (TO BE COMPLETED)
  !
  ! Inputs
  !       ielem: current element you are on
  !       iface: current face you are on
  !       kface: face of adjoining element
  !       ifacenodes_On: ifacenodes(nodesperface*nfacesperelem) kfacenode flattened into a single vector volumetric node index of face node
  !       n_S_2d_max: maximum number of nodes on a face
  !       n_S_1d_On: number of nodes in each direction for ielem
  !       n_S_2d_On: number of nodes on each face for ielem
  !       x_S_1d_On: one-dimensional computational coordinates for ielem
  !       n_S_1d_Off: number of nodes in each direction for kelem
  !       n_S_2d_Off: number of nodes on each face for kelem
  !       x_S_1d_Off: one-dimensional computational coordinates for kelem
  !       n_S_1d_Mort: number of nodes in each direction on the Morter
  !       n_S_2d_Mort: number of nodes on the mortar face
  !       x_S_1d_Mort: one-dimensional computational coordiantes on the Mortar
  !       pinv: vector with the norm weights in one-dimension
  !       vg_2d_On: two dimensional array (dim,n_S_2d_On) holding the primative variables over the iface of element ielem 
  !       vg_2d_Off: two dimensional array (dim,n_S_2d_Off) holding the primative variables over the iface of element kelem 
  !       wg_Mort_On: two dimensional array (dim,n_S_2d_Mort) holding the projection of the entropy variables form ielem to the mortar
  !       wg_Mort_Off:two dimensional array (dim,n_S_2d_Mort) holding the projection of the entropy variables form kelem to the mortar
  !       cnt_Mort_Off: vector size n_S_2d_Mort gives the correspondance between the jth mortar nod for the ielem to the node on the kelem
  !       Intrp_On: from the mortar to ielem
  !       Extrp_Off: from kelem to the mortar
  !
  ! Outputs
  ! 
  ! Notes: 
  !
  !=============================================================================
  subroutine Inviscid_SAT_Non_Conforming_Interface_Mod_SAT(ielem, iface, ifacenodes_On ,n_S_2d_max,  &
                                                     n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                     n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                     pinv, nx_Off_ghst,                              &
                                                     vg_2d_On,  vg_2d_Off, wg_2d_On, wg_2d_Off,      &
                                                     IOn2Off_x1, IOn2Off_x2,                         &
                                                     IOff2On_x1, IOff2On_x2 )

    use referencevariables,   only: nequations, ndim
    use variables,            only: facenodenormal, Jx_r, gsat, ef2e
    use initcollocation,      only: ExtrpXA2XB_2D_neq, element_properties
    use initgrid,             only: map_face_orientation_k_On_2_k_Off

    implicit none

    integer,                    intent(in) :: ielem, iface
    integer,                    intent(in) :: n_S_1d_On, n_S_1d_Off
    integer,                    intent(in) :: n_S_2d_On, n_S_2d_Off, n_S_2d_max
    integer,   dimension(:),    intent(in) :: ifacenodes_On

    real(wp),  dimension(:),    intent(in) :: x_S_1d_On, x_S_1d_Off
    real(wp),  dimension(:),    intent(in) :: pinv
    real(wp),  dimension(:,:),  intent(in) :: nx_Off_ghst
    real(wp),  dimension(:,:),  intent(in) :: vg_2d_On,   vg_2d_Off
    real(wp),  dimension(:,:),  intent(in) :: wg_2d_On, wg_2d_Off
    real(wp),  dimension(:,:),  intent(in) :: IOn2Off_x1,IOn2Off_x2
    real(wp),  dimension(:,:),  intent(in) :: IOff2On_x1,IOff2On_x2
    
    real(wp),  dimension(ndim)             :: nx, nx_On, nx_Off, nx_Ave
    real(wp),  dimension(nequations)       :: fn, fstar
    real(wp),  dimension(nequations)       :: vg_On, vg_Off

    real(wp), allocatable, dimension(:,:) :: FA, FB
    real(wp), allocatable, dimension(:,:) :: Up_diss_On, Up_diss_Off

    real(wp), allocatable, dimension(:,:) :: wg_On2Off, wg_Off2On

    integer                               :: i, j, k, l, orientation
    integer                               :: inode, jnode, kelem

    continue
     
      allocate(wg_On2Off(nequations,1:n_S_2d_Off))
      allocate(wg_Off2On(nequations,1:n_S_2d_On))

      !-- interpolate/extrapolate wg from On to Off and vice versa
      call ExtrpXA2XB_2D_neq(nequations, n_S_1d_On , n_S_1d_Off, x_S_1d_On , x_S_1d_Off, &
                             wg_2d_On ,wg_On2Off, IOn2Off_x1, IOn2Off_x2) 
      call ExtrpXA2XB_2D_neq(nequations, n_S_1d_Off, n_S_1d_On , x_S_1d_Off, x_S_1d_On , &
                             wg_2d_Off,wg_Off2On, IOff2On_x1, IOff2On_x2)

      allocate(FA (nequations,n_S_2D_Off ))
      allocate(FB (nequations,n_S_2D_On  ))
      allocate(Up_diss_On  (nequations,n_S_2d_On  ))
      allocate(Up_diss_Off (nequations,n_S_2d_Off ))

      orientation = ef2e(7,iface,ielem)                                        !-- orientation of the off element relative to the on element

!=========
!         Skew-symmetric matrix portion of SATs (Entropy Conservative without Mortar)
!=========
      kelem = ef2e(2,iface,ielem)


      On_Element_1:do i = 1, n_S_2d_On                                         ! On element Loop over 2D LGL points
 
        jnode =  n_S_2d_On*(iface-1) + i                                       ! Index in facial ordering
              
        inode = ifacenodes_On(jnode)                                           ! Volumetric node index corresponding to facial node index

        nx_On = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)                ! On element outward facing normals

        fn(:) = normalflux(vg_2d_On(:,i), nx_On(:), nequations)                ! One point flux based on vg_On and nx

        Off_Element_1:do k = 1, n_S_2d_Off                                     ! Off_Element Loop over 2D LGL points
         
          nx_Off(:) = nx_Off_ghst(:,k)                                         ! On element outward facing normals

          nx_Ave(:) = 0.5_wp*(nx_On(:)+nx_Off(:))                              ! Average of normal(i) and normal(j)

!         FA(:,k) = EntropyConsistentFlux     (vg_2d_On(:,i), vg_2d_Off(:,k), nx_Ave(:),nequations) ! Entropy conservative flux
          FA(:,k) = Entropy_KE_Consistent_Flux(vg_2d_On(:,i), vg_2d_Off(:,k), nx_Ave(:),nequations) ! Entropy conservative flux

        enddo Off_Element_1                                                    ! End Off_Element

        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_On,x_S_1d_Off,x_S_1d_On, &
                               FA,FB, IOff2On_x1, IOff2On_x2)  ! Extrapolate f^S_x

        l =  map_face_orientation_k_On_2_k_Off(i,orientation,n_S_1d_On)        ! Correct for face orientation and shift back to 1:n_S_2d_On

        fstar = FB(:,l)

        gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1)*(fn - fstar)       ! SAT penalty:  subtract the on-element contribution and replace with penalty 

      enddo On_Element_1

!=========
!         Inviscid interface dissipation (Entropy Stable Upwinding of SATs) without mortar
!=========

      !=========
      !   On element contribution to the dissipation: R_ON^T*Lambda_On*Ron*(wg_ON-wg_Off2ON)
      !=========
      On_Element_2:do j = 1, n_S_2d_On                                         ! On element loop over data

        jnode =  n_S_2d_On*(iface-1) + j                                       ! Index in facial ordering (bucket is padded so n_S_2d_max is needed)

        inode = ifacenodes_On(jnode)                                           ! Volumetric node index corresponding to facial node index

        nx = + Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)                 ! Outward facing normal of facial node (On-Element)


        l = map_face_orientation_k_On_2_k_Off(j,orientation,n_S_1d_On)         ! Correct for face orientation and shift back to 1:n_S_2d_On


        call entropy_to_primitive(wg_2d_On (:,j),vg_On (:),nequations)         ! Entropy -> primitive variables:  On_element
        call entropy_to_primitive(wg_Off2On(:,l),vg_Off(:),nequations)         ! Entropy -> primitive variables: Off_element

        gsat(:,inode,ielem) = gsat(:,inode,ielem) + &
         0.50_wp * pinv(1) * SAT_Vis_Diss(nequations,vg_On(:),vg_Off(:),nx(:)) ! On-element Viscous penalty contribution

      enddo On_Element_2                                                       ! End On element loop

      !=========
      !   Off element contribution to the dissipation: R_Off^T*Lambda_Off*R_Off*(wg_On2Off-wg_Off)
      !   construct 2D plane data in the ordering of the On element
      !=========
      Off_Element_2: do j = 1, n_S_2d_Off

        l = map_face_orientation_k_On_2_k_Off(j,orientation,n_S_1d_Off)        ! Correct for face orientation

        nx(:) = nx_Off_ghst(:,l)                                               ! Outward facing normal of facial node


        call entropy_to_primitive(wg_On2Off (:,l),vg_On (:),nequations)        ! Entropy -> primitive variables:  On_element
        call entropy_to_primitive(wg_2d_Off(:,l),vg_Off(:),nequations)         ! Entropy -> primitive variables: Off_element

       Up_diss_Off(:,l) = SAT_Vis_Diss(nequations,vg_On(:),vg_Off(:),nx(:))    ! Viscous dissipation on Mortar based on L-R states

      enddo Off_Element_2
 
                                                                                      ! Restrict data plane from Mortar to on-face plane
      call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off ,n_S_1d_On ,x_S_1d_Off ,x_S_1d_On, &
                             Up_diss_Off,Up_diss_On, IOff2On_x1, IOff2On_x2)

      !-- add dissipation to sat

      On_Elem_2: do i = 1, n_S_2d_On                                           ! On-element loop: Begin
        
        jnode =  n_S_2d_On*(iface-1) + i                                       ! Index in facial ordering
              
        inode = ifacenodes_On(jnode)                                           ! Volumetric node index corresponding to facial node index
             
        gsat(:,inode,ielem) = gsat(:,inode,ielem) + 0.50_wp*pinv(1) * Up_diss_On(:,i)! On-element Viscous penalty contribution

      end do On_Elem_2                                                         ! On-element loop: end

      deallocate(FA,FB,Up_diss_On,Up_diss_Off)

  end subroutine Inviscid_SAT_Non_Conforming_Interface_Mod_SAT

  !============================================================================
  !
  ! Purpose: Constructs the inviscid SAT for the modified SAT approach where the 
  !          symmetric portion of the sat is for the xil computational direction and 
  !          xm Cartesian direction
  !\begin{equation*}
  !  P^{-1}E_{xil,ielem}\circ\[Jdxildx]_{ielem}\circ F^{SC}(q_ielem,q_ielem)1_{ielem}
  !  -\frac{1}{2}P^{-1}\left([Jdxildxm]_{ielem} Evtok + Evtok [Jdxildx_m]_{kelem}\right)\circ 
  !  F^{SC}(q_{ielem},q_{kelem}) 1_{kelem}
  !\end{equation*}
  !
  ! The upwinding contribution is given as (TO BE COMPLETED)
  !
  ! Inputs
  !       ielem: current element you are on
  !       iface: current face you are on
  !       kface: face of adjoining element
  !       ifacenodes_On: ifacenodes(nodesperface*nfacesperelem) kfacenode flattened into a single vector volumetric node index of face node
  !       n_S_2d_max: maximum number of nodes on a face
  !       n_S_1d_On: number of nodes in each direction for ielem
  !       n_S_2d_On: number of nodes on each face for ielem
  !       x_S_1d_On: one-dimensional computational coordinates for ielem
  !       n_S_1d_Off: number of nodes in each direction for kelem
  !       n_S_2d_Off: number of nodes on each face for kelem
  !       x_S_1d_Off: one-dimensional computational coordinates for kelem
  !       n_S_1d_Mort: number of nodes in each direction on the Morter
  !       n_S_2d_Mort: number of nodes on the mortar face
  !       x_S_1d_Mort: one-dimensional computational coordiantes on the Mortar
  !       pinv: vector with the norm weights in one-dimension
  !       vg_2d_On: two dimensional array (dim,n_S_2d_On) holding the primative variables over the iface of element ielem 
  !       vg_2d_Off: two dimensional array (dim,n_S_2d_Off) holding the primative variables over the iface of element kelem 
  !       wg_Mort_On: two dimensional array (dim,n_S_2d_Mort) holding the projection of the entropy variables form ielem to the mortar
  !       wg_Mort_Off:two dimensional array (dim,n_S_2d_Mort) holding the projection of the entropy variables form kelem to the mortar
  !       cnt_Mort_Off: vector size n_S_2d_Mort gives the correspondance between the jth mortar nod for the ielem to the node on the kelem
  !       Intrp_On: from the mortar to ielem
  !       Extrp_Off: from kelem to the mortar
  !
  ! Outputs
  ! 
  ! Notes: 
  !
  !=============================================================================
  subroutine Inviscid_SAT_Non_Conforming_Interface_Mod_SAT_OLD(ielem, iface, ifacenodes_On ,n_S_2d_max, &
                                                     n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                     n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                     n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,            &
                                                     pinv,                                           &
                                                     vg_2d_On,  vg_2d_Off, wg_Mort_On, wg_Mort_Off,  &
                                                     cnt_Mort_Off, Intrp_On, Extrp_Off,nx_Off_ghst)

    use referencevariables,   only: nequations, ndim
    use variables,            only: facenodenormal, Jx_r, Jx_facenodenormal_Gau, gsat, ef2e
    use initcollocation,      only: ExtrpXA2XB_2D_neq, ExtrpXA2XB_2D_neq_k, element_properties

    implicit none

    integer,                    intent(in) :: ielem, iface
    integer,                    intent(in) :: n_S_1d_On, n_S_1d_Off, n_S_1d_Mort
    integer,                    intent(in) :: n_S_2d_On, n_S_2d_Off, n_S_2d_Mort, n_S_2d_max
    real(wp),  dimension(:),    intent(in) :: x_S_1d_On, x_S_1d_Off, x_S_1d_Mort
    real(wp),  dimension(:),    intent(in) :: pinv
    integer,   dimension(:),    intent(in) :: ifacenodes_On
    integer,   dimension(:),    intent(in) :: cnt_Mort_Off
    real(wp),  dimension(:,:),  intent(in) :: vg_2d_On,   vg_2d_Off
    real(wp),  dimension(:,:),  intent(in) :: wg_Mort_On, wg_Mort_Off
    real(wp),  dimension(:,:),  intent(in) :: Intrp_On, Extrp_Off, nx_Off_ghst
    
    real(wp),  dimension(ndim)             :: nx, nx_On, nx_Off, nx_Ave
    real(wp),  dimension(nequations)       :: fn, fstar
    real(wp),  dimension(nequations)       :: vg_On, vg_Off

    real(wp), allocatable, dimension(:,:) :: FA , FB
    real(wp), allocatable, dimension(:,:) :: FC_Mort_On
    real(wp), allocatable, dimension(:,:) :: Up_diss_Mort, Up_diss_On

    integer                               :: i, j, k, l
    integer                               :: ival, jval
    integer                               :: inode, jnode, kelem

    continue

      allocate(FA (nequations,n_S_2D_Off ))
      allocate(FB (nequations,n_S_2D_Mort))

      allocate( FC_Mort_On (nequations,n_S_2D_Mort))
      allocate(Up_diss_Mort(nequations,n_S_2d_Mort))
      allocate(Up_diss_On  (nequations,n_S_2d_On  ))

!=========
!         Skew-symmetric matrix portion of SATs (Entropy Stable through Mortar)
!=========
      kelem = ef2e(2,iface,ielem)


      On_Element_1:do i = 1, n_S_2d_On                                       ! On_Element Loop over 2D LGL points
 
        jnode =  n_S_2d_On*(iface-1) + i                                     ! Index in facial ordering
              
        inode = ifacenodes_On(jnode)                                         ! Volumetric node index corresponding to facial node index

        nx_On = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)              ! On element outward facing normals

        fn(:) = normalflux(vg_2d_On(:,i), nx_On(:), nequations)              ! One point flux based on vg_On and nx

        Off_Element_1:do k = 1, n_S_2d_Off                                   ! Off_Element Loop over 2D LGL points
         
          nx_Off(:) = nx_Off_ghst(:,k)                                       ! On element outward facing normals

          nx_Ave(:) = 0.5_wp*(nx_On(:)+nx_Off(:))                            ! Average of normal(i) and normal(j)

!         FA(:,k) = EntropyConsistentFlux     (vg_2d_On(:,i), vg_2d_Off(:,k), nx_Ave(:),nequations) ! Entropy conservative flux
          FA(:,k) = Entropy_KE_Consistent_Flux(vg_2d_On(:,i), vg_2d_Off(:,k), nx_Ave(:),nequations) ! Entropy conservative flux

        enddo Off_Element_1                                                  ! End Off_Element

        call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                               FA,FB,Extrp_Off, Extrp_Off)                   ! Extrapolate f^S_x

        On_Mortar_1:do j = 1, n_S_2d_Mort                                    ! Mortar loop over 2D Gauss points
  
          l = cnt_Mort_Off(j)                                                 ! Correct for face orientation and shift back to 1:n_S_2d_Mort

          FC_Mort_On(:,j) = FB(:,l)                                          ! Outward facing component of f^S on mortar

        enddo On_Mortar_1                                                    ! End Mortar loop

        ival = mod(i-1,n_S_1d_On) + 1 ; jval = (i-ival) / n_S_1d_On + 1 ;    ! Decode On-element point coordinates (i,j) from planar coordinates

        call ExtrpXA2XB_2D_neq_k(nequations,n_S_1d_Mort,n_S_1d_On,ival,jval,x_S_1d_Mort,x_S_1d_On, &
                                 FC_Mort_On,fstar,Intrp_On,Intrp_On)         ! Restrict planar data to the (ival,jval) point

        gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1)*(fn - fstar)     ! SAT penalty:  subtract the on-element contribution and replace with penalty 

      enddo On_Element_1

!=========
!         Inviscid interface dissipation (Entropy Stable Upwinding of SATs) using the mortar
!=========
if(.true.)then
      On_Mortar_2:do j = 1, n_S_2d_Mort                                      ! Mortar loop over data
  
        jnode =  n_S_2d_max*(iface-1) + j                                    ! Index in facial ordering (bucket is padded so n_S_2d_max is needed)
  
        nx(:) = Jx_facenodenormal_Gau(:,jnode,ielem)                         ! Outward facing normal of facial node

            l = cnt_Mort_Off(j)                                              ! Correct for face orientation and shift back to 1:n_S_2d_Mort

       call entropy_to_primitive(wg_Mort_On (:,j),vg_On (:),nequations)      ! Entropy -> primitive variables:  On_element
       call entropy_to_primitive(wg_Mort_Off(:,l),vg_Off(:),nequations)      ! Entropy -> primitive variables: Off_element

       Up_diss_Mort(:,j) = SAT_Vis_Diss(nequations,vg_On(:),vg_Off(:),nx(:)) ! Viscous dissipation on Mortar based on L-R states

      enddo On_Mortar_2                                                      ! End mortar loop

                                                                                 ! Restrict data plane from Mortar to on-face plane
      call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Mort ,n_S_1d_On ,x_S_1d_Mort ,x_S_1d_On, &
                             Up_diss_Mort,Up_diss_On, Intrp_On, Intrp_On )

      On_Elem_2: do i = 1, n_S_2d_On                                         ! On-element loop: Begin
        
        jnode =  n_S_2d_On*(iface-1) + i                                     ! Index in facial ordering
              
        inode = ifacenodes_On(jnode)                                         ! Volumetric node index corresponding to facial node index
             
        gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * Up_diss_On(:,i)! On-element Viscous penalty contribution

      end do On_Elem_2                                                       ! On-element loop: end

      deallocate(Up_diss_On,Up_diss_Mort)
endif

      deallocate(FC_Mort_On,FA,FB)

  end subroutine Inviscid_SAT_Non_Conforming_Interface_Mod_SAT_OLD

  !============================================================================
  !       NONCONFORMING VISCOUS interface SATs 
  !============================================================================
  
  subroutine Viscous_SAT_Non_Conforming_Interface(ielem, kelem, iface, kface,                     &
                                                  ifacenodes_On, ifacenodes_Off, n_S_2d_max,      &
                                                  n_S_1d_On  ,n_S_2d_On  ,x_S_1d_On  ,            &
                                                  n_S_1d_Off ,n_S_2d_Off ,x_S_1d_Off ,            &
                                                  n_S_1d_Mort,n_S_2d_Mort,x_S_1d_Mort,            &
                                                  pinv,                                           &
                                                    vg_2d_On,   vg_2d_Off,                        &
                                                  phig_2d_On, phig_2d_Off,                        &
                                                  wg_Mort_On, wg_Mort_Off,                        &
                                                  nx_2d_Off, mut_2d_Off, Jx_r_2d_Mort,            &
                                                  cnt_Mort_Off, Intrp_On, Extrp_Off)

    use referencevariables,   only: nequations, ndim
    use variables,            only: facenodenormal, Jx_r, Jx_facenodenormal_Gau, gsat, mut, ef2e
    use collocationvariables, only: l00, l01, ldg_flip_flop_sign, alpha_ldg_flip_flop
    use initcollocation,      only: ExtrpXA2XB_2D_neq
    use initgrid,             only: map_face_orientation_k_On_2_k_Off

    implicit none

    integer,                    intent(in) :: ielem, kelem, iface, kface
    integer,                    intent(in) :: n_S_1d_On, n_S_1d_Off, n_S_1d_Mort
    integer,                    intent(in) :: n_S_2d_On, n_S_2d_Off, n_S_2d_Mort, n_S_2d_max
    real(wp),  dimension(:),    intent(in) :: x_S_1d_On, x_S_1d_Off, x_S_1d_Mort
    integer,   dimension(:),    intent(in) :: ifacenodes_On, ifacenodes_Off
    integer,   dimension(:),    intent(in) :: cnt_Mort_Off
    real(wp),  dimension(:),    intent(in) :: pinv
    real(wp),  dimension(:,:),  intent(in) :: vg_2d_On,   vg_2d_Off
    real(wp),  dimension(:,:),  intent(in) ::             nx_2d_Off
    real(wp),  dimension(:,:),  intent(in) :: wg_Mort_On, wg_Mort_Off
    real(wp),  dimension(:,:,:),intent(in) :: phig_2d_On, phig_2d_Off
    real(wp),  dimension(:),    intent(in) ::              mut_2d_Off
    real(wp),  dimension(:),    intent(in) ::             Jx_r_2d_Mort
    real(wp),  dimension(:,:),  intent(in) :: Intrp_On,   Extrp_Off
    
    real(wp),  dimension(ndim)             :: nx, nx_Off
    real(wp),  dimension(nequations)       :: fV_Del, SAT_Pen
    real(wp),  dimension(nequations)       ::   vg_On,   vg_Off
    real(wp),  dimension(nequations,ndim)  :: phig_On, phig_Off

    real(wp), dimension(nequations,nequations) :: matrix_ip

    real(wp), allocatable, dimension(:,:) :: fV_2d_Mort, fV_Mort_On, fV_Mort_Off, fV_2d_On, fV_2d_Off
    real(wp), allocatable, dimension(:,:) :: IP_2d_On, IP_2d_Mort
    real(wp), allocatable, dimension(:,:) :: IOn2Off, IOff2On

    integer                               :: i, j, k, l, orientation
    integer                               :: inode, jnode

    real(wp)                                   ::  mut_Off
    real(wp)                                   :: l01_ldg_flip_flop

    continue

      allocate(fV_Mort_On (nequations,n_S_2d_Mort))
      allocate(fV_Mort_Off(nequations,n_S_2d_Mort))
      allocate(fV_2d_Mort (nequations,n_S_2d_Mort))

      allocate(fV_2d_On   (nequations,n_S_2d_On  ))
      allocate(fV_2d_Off  (nequations,n_S_2d_Off ))

      allocate(IP_2d_Mort (nequations,n_S_2d_Mort))
      allocate(IP_2d_On   (nequations,n_S_2d_On  ))

      allocate(IOff2On(1:n_S_1d_On,1:n_S_1d_Off)) ; IOff2On = matmul(Intrp_On,Extrp_Off) ;

      orientation = ef2e(7,iface,ielem)

      !  =======
      !  LDG viscous dissipation: Connects Off with On interfaces through mortar
      !  =======

      Off_Elem_3:do k = 1, n_S_2d_Off

          nx_Off(:  ) =   nx_2d_Off(:,  k)         ! Outward facing normal of facial node 

          vg_Off(:  ) =   vg_2d_Off(:,  k)         ! primitive variables of facial node

        phig_Off(:,:) = phig_2d_Off(:,:,k)         ! viscous gradients of facial node

         mut_Off      =  mut_2d_Off(    k)         ! Turbulent viscosity of facial node

        fV_2d_Off(:,k) = normalviscousflux(vg_Off(:), phig_Off(:,:), nx_Off(:), nequations, mut_Off)

      enddo Off_Elem_3

      call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_On,x_S_1d_Off,x_S_1d_On, & ! Extrapolate f^S_x
                            fV_2d_Off,fV_2d_On, IOff2On, IOff2On )                   ! fV_2d_On still has wrong ordering and orientation

      !  =======
      !  IP dissipation: formed on the common mortar between element interfaces
      !  =======

      On_Mortar_4:do j = 1, n_S_2d_Mort

        jnode =  n_S_2d_max*(iface-1) + j                        ! Index in facial ordering

        nx(:) = Jx_facenodenormal_Gau(:,jnode,ielem)             ! Outward facing normal of facial node

            l = cnt_Mort_Off(j)                                  ! Correct for face orientation and shift back to 1:n_S_2d_Mort

        call entropy_to_primitive(wg_Mort_On (:,j),vg_On (:),nequations)
        call entropy_to_primitive(wg_Mort_Off(:,l),vg_Off(:),nequations)

        matrix_ip = 0.5_wp * pinv(1) * (matrix_hatc_node(vg_On (:),nx,nx,nequations) &
                                     +  matrix_hatc_node(vg_Off(:),nx,nx,nequations)) / Jx_r_2d_Mort(j)

        IP_2d_Mort(:,j) = - l00*matmul(matrix_ip,wg_Mort_On(:,j)-wg_Mort_Off(:,l))
          
      enddo On_Mortar_4

      call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Mort,n_S_1d_On,x_S_1d_Mort,x_S_1d_On, &
                                IP_2d_Mort(:,:),IP_2d_On(:,:),Intrp_On,Intrp_On)

      On_Elem_4:do i = 1, n_S_2d_On                              !  Onface sweep 
    
        jnode =  n_S_2d_On*(iface-1) + i                         ! Index in facial ordering
          
        inode = ifacenodes_On(jnode)                             ! Volumetric node index corresponding to facial node index

          vg_On (:  ) =   vg_2d_On(:,  i)                        ! On-element primitive face data
        phig_On (:,:) = phig_2d_On(:,:,i)                        ! On-element viscous derivatives (LDG)

        nx(:) = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)  ! Outward facing normal of facial node

        l01_ldg_flip_flop = l01*(1.0_wp - ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)

        l =  map_face_orientation_k_On_2_k_Off(i,orientation,n_S_1d_On)        ! Correct for face orientation and shift back to 1:n_S_2d_On

        fV_Del(:) = normalviscousflux(vg_On (:), phig_On (:,:), nx, nequations, mut(inode,ielem)) &
                  - (-1 * fV_2d_On(:,l))                         ! -1 accounts for opposite direction of outward faceing normal

        SAT_Pen(:) = + l01_ldg_flip_flop*fV_Del(:) + IP_2d_On(:,i)

        gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * SAT_Pen(:)

      enddo On_Elem_4

      deallocate(fV_2d_Mort, fV_Mort_Off, fV_Mort_On, fV_2d_On, fV_2d_Off)
      deallocate(IP_2d_On, IP_2d_Mort)
      deallocate(IOff2On)

  end subroutine Viscous_SAT_Non_Conforming_Interface

  !============================================================================
  
  subroutine nse_calcrhsexplicit(irk,tin)

    ! This subroutine calculates the time derivative of the conserved variables at every node.

    use variables
    use referencevariables
    use controlvariables
    use initcollocation, only: element_properties
    use mpimod
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer,  intent(in) :: irk
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode,ielem, jdir

!   integer :: i_err

    ! update the primitive and entropy variables and
    ! the LDG/LDC gradients
    call nse_reconcilestates() ! (navierstokes)

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      !                                        __
      !  Calculate the elementwise Divergence  \/ * (F - Fv)
       
!     HACK
!     if(IMEX_element == 'explicit') call Flux_Divergence(tin,ielem)    !  result in divF

      !  Form the elementwise SAT_Penalties

!     if(IMEX_penalty == 'explicit') call SAT_Penalty(tin,ielem)        !  result in gsat
!     HACK
      !  
      ! compute time derivative
      ! 
      ! reset the time derivative to zero
      Fexp(:,:,ielem,irk) = 0.0_wp
      ! loop over all nodes in the element
      do inode = 1, nodesperelem
        ! add the contribution from the flux divergence in each direction
        if(IMEX_element == 'explicit')then
          do jdir = 1,ndim
            Fexp(:,inode,ielem,irk) = Fexp(:,inode,ielem,irk) - divf(:,jdir,inode,ielem)/Jx_r(inode,ielem)
          end do
        endif
        ! add the contribution from the boundary and interface penalties
        if(IMEX_penalty == 'explicit')then
          Fexp(:,inode,ielem,irk) = Fexp(:,inode,ielem,irk) + gsat(:,inode,ielem)/Jx_r(inode,ielem)
        endif
      end do

    end do

    return
  end subroutine nse_calcrhsexplicit

  !============================================================================

  subroutine nse_calcrhsimplicit(irk,tin)
    ! This subroutine calculates the time derivative of the
    ! conserved variables at every node.
    use variables
    use referencevariables
    use controlvariables
    use initcollocation, only: element_properties

    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer,  intent(in) :: irk
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode,ielem, jdir

    ! update the primitive and entropy variables and
    ! the LDG/LDC gradients
    call nse_reconcilestates() ! (navierstokes)

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      !                                        __
      !  Calculate the elementwise Divergence  \/ * (F - Fv)
       
!     HACK
!     if(IMEX_element == 'implicit') call Flux_Divergence(tin,ielem)    !  result in divF

      !  Form the elementwise SAT_Penalties

!     if(IMEX_penalty == 'implicit') call SAT_Penalty(tin,ielem)        !  result in gsat
!     HACK
      !  
      ! compute the "IMPLICIT" flux Fimp
      ! 
      ! reset the time derivative to zero
      Fimp(:,:,ielem,irk) = 0.0_wp
      ! loop over all nodes in the element
      do inode = 1, nodesperelem
        ! add the contribution from the flux divergence in each direction
        if(IMEX_element == 'implicit')then
          do jdir = 1,ndim
            Fimp(:,inode,ielem,irk) = Fimp(:,inode,ielem,irk) - divf(:,jdir,inode,ielem)/Jx_r(inode,ielem)
          end do
        endif
        ! add the contribution from the boundary and interface penalties
        if(IMEX_penalty == 'implicit')then
          Fimp(:,inode,ielem,irk) = Fimp(:,inode,ielem,irk) + gsat(:,inode,ielem)/Jx_r(inode,ielem)
        endif
      end do

    end do

    return
  end subroutine nse_calcrhsimplicit

  !============================================================================

  subroutine set_Flow_parameters()

    use nsereferencevariables

    ! specific heat ratio
    gamma0 = 1.4_wp
    ! molecular weight
    MW0    = 28.97_wp
    ! reference pressure
    p0 = patm
    ! reference temperature
    t0 = tref
    ! calculate reference density
    rho0 = p0*MW0/(Ru*T0)
    ! calculate reference speed of sound
    csound0 = sqrt(gamma0*Ru/MW0*T0)
    ! calculate reference velocity
    U0 = Mach0*csound0

    gamI  = 1.0_wp / gamma0
    gm1   = gamma0 - one
    gm1M2 = gm1*Mach0*Mach0
    gm1M2I= 1.0_wp / gm1M2

    gm1og = gm1/gamma0
    gp1og = (gamma0 + one)/gamma0
    gM2   = u0*u0/(Ru/MW0*T0)
    gM2I  = one / gM2

    ! Where viscous routines are called whether or
    ! not the problem is viscous, we need the viscous
    ! fluxes to be zero. This is accomplished by always
    ! multiplying by Re0inv instead of dividing by Re0.
    Re0inv = 0.0_wp
    if (viscous) Re0inv = 1.0_wp/Re0

    return
  end subroutine set_Flow_parameters

  subroutine set_InitialCondition(bfuncin)
! this subroutine sets the pointer for the exactSolution procedure
! it will default to uniformFreeStream
    character(180), intent(in) :: bfuncin

    integer :: ispline

    select case(trim(bfuncin))
      case('UniformFreeStream')
        InitialSubroutine => UniformFreeStream
      case('ExactSolutionViscousShock')
        InitialSubroutine => viscousShockFull
      case('ExactSolutionIsentropicVortex')
        InitialSubroutine => isentropicVortexFull
      case('ExactSolutionSupersonicVortex')
        InitialSubroutine => supersonicvortexFull
      case('ShockVortex')
        InitialSubroutine => ShockVortexInteraction
      case('SodsProblem')
        InitialSubroutine => SodsProblemICBC
      case('QuiescentFluid')
        InitialSubroutine => UniformFreeStream
      case('taylorGreenVortex')
        InitialSubroutine => taylorGreen
      case('Constant_rhoS')
        InitialSubroutine => Constant_rhoS
      case('Potential_cylinder')
        InitialSubroutine => Potential_Flow_Around_cylinder
      case('PreserveFreeStream')
        InitialSubroutine => UniformFreeStream
      case('VortexKicker')
        InitialSubroutine => isentropicVortexFull
      case('GFIT')
        InitialSubroutine => GFIT
      case('BLayerProfile')
        ispline = 1
        InitialSubroutine => BLayerProfile
!     case('BLayerProfile-01':'BLayerProfile-99')
!       read(bfuncin(15:16),*) ispline
!       InitialSubroutine => BLayerProfileS
      case('Kelvin-Helmoholtz')
        InitialSubroutine => kelvin_helmoholtz
      case default
        InitialSubroutine => UniformFreeStream
        write(*,*)'no match for InitialCondition specificiation'
        write(*,*)'using default UniformFreeStream'
     end select

     return
  end subroutine set_InitialCondition

  pure subroutine UniformFreeStream(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    real(wp) :: alpha

    alpha = uniformFreeStreamAOA*pi/180._wp

    Vx = zero
    ! density
    Vx(1) = one
    ! velocity
    Vx(2) = cos(alpha)
    Vx(3) = sin(alpha)
    Vx(4) = zero
    ! Temperature
    Vx(5) = one
    
    fv = zero

    return
  end subroutine UniformFreeStream

  pure subroutine GFIT(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    real(wp) :: alpha

    alpha = uniformFreeStreamAOA*pi/180._wp

    Vx = zero
    ! density
    Vx(1) = one
    ! velocity
    Vx(2) = cos(alpha)
    Vx(3) = sin(alpha)
    Vx(4) = zero
    if(xin(2) <= zero) Vx(2:4) = zero
    ! Temperature
    Vx(5) = one
    
    fv = zero

    return
  end subroutine GFIT

  pure subroutine QuiescentFluid(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    Vx = zero
    ! density
    Vx(1) = one
    ! velocity
    Vx(2) = zero
    Vx(3) = zero
    Vx(4) = zero
    ! Temperature
    Vx(5) = one
    
    fv = zero

    return
  end subroutine QuiescentFluid

  pure subroutine taylorGreen(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut
    real(wp) :: p
    
    real(wp) :: fcorr
    
    fcorr = two*pi/6.2831853071795862_wp
    

    p = one + gM2*( (cos(two*xin(1))+cos(two*xin(2))) * (cos(two*xin(3))+two) )/16.0_wp

    Vx    = zero

    Vx(5) = one

    Vx(2) = +sin(xin(1))*cos(xin(2))*cos(xin(3))
    Vx(3) = -cos(xin(1))*sin(xin(2))*cos(xin(3))
    Vx(4) = zero

    ! Density
    Vx(1) = p/Vx(5)

    ! Viscous flux
    fv = zero

    return
  end subroutine taylorGreen 

  pure subroutine InviscidWall(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer,  intent(in)    :: neqin, nd

    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out)   :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in)    :: Jx(3)
    real(wp), intent(in)    :: xin(3)
    real(wp), intent(in)    :: tin, mut

    real(wp), parameter  :: Mirror = -1.0_wp
    real(wp), dimension(3)   :: VNormal, Vtangent, JxNormal

    JxNormal = Jx(:) / Magnitude(Jx)

    VNormal(:)   = dot_product(Vx(2:4),JxNormal(:)) * JxNormal(:)
    VTangent(:)  = Vx(2:4) - VNormal(:)

    Vx(1)   = Vx(1)
    Vx(2:4) = VTangent(:) + 1.0_wp*Mirror*VNormal(:)
    Vx(5)   = Vx(5)

    fv = 0.0_wp

    return
  end subroutine InviscidWall

  pure subroutine SymmetryPlane(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    real(wp), parameter  :: Mirror = -1.0_wp
    real(wp), dimension(3)   :: VNormal, Vtangent, JxNormal

    JxNormal = Jx(:) / Magnitude(Jx)

    VNormal(:)   = dot_product(Vx(2:4),JxNormal(:)) * JxNormal(:)
    VTangent(:)  = Vx(2:4) - VNormal(:)

    Vx(1)   = Vx(1)
    Vx(2:4) = VTangent(:) + Mirror*VNormal(:)
    Vx(5)   = Vx(5)

    fv = Mirror*normalviscousflux(Vx, phi, Jx, neqin,mut)
    
    return
  end subroutine SymmetryPlane

  !============================================================================
  
  pure subroutine NoSlipWallAdiabatic(prim_int,grad_entr_int,f_v,normal,pos,time, &
      & n_eq,n_d,mut)
    
    ! Load modules
    use nsereferencevariables
    
    ! Nothing is implicitly defined
    implicit none
    
    integer,  intent(in)    :: n_eq, n_d
    real(wp), intent(out)   :: f_v(n_eq)
    real(wp), intent(inout) :: prim_int(n_eq)
    real(wp), intent(in)    :: grad_entr_int(n_eq,3)
    real(wp), intent(in)    :: normal(3)
    real(wp), intent(in)    :: pos(3)
    real(wp), intent(in)    :: time, mut

    real(wp), parameter  :: mirror = -1.0_wp
    real(wp), parameter  :: no_slip = -0.0_wp

    real(wp), dimension(3)    :: unit_normal, normal_vel, tangent_vel
    real(wp), dimension(n_eq,3) :: grad_prim_int, grad_prim_int_normal, &
      & grad_prim_int_tangent

    integer                         :: i

    ! Unit normal direction
    unit_normal(:) = normal(:) / Magnitude(normal(:))

    normal_vel(:)   = dot_product(prim_int(2:4),unit_normal) * unit_normal(:)
    tangent_vel(:)  = prim_int(2:4) - normal_vel(:)

    ! Set primitive variables in the ghost node
    ! ======================================== 
    ! Density
    prim_int(1) =  prim_int(1)
    
    ! Velocity vector
    prim_int(2:4) = no_slip*tangent_vel(:) + mirror*normal_vel(:)

    ! Temperature
    prim_int(5) = prim_int(5)
    
   
    ! Grad(V) = dVdW Grad(W)  = dvdw phi
    grad_prim_int = MatMul(dVdW(prim_int,n_eq),grad_entr_int)


    ! Set gradient of the primitive variables in the ghost node
    ! =========================================================

    ! Normal component of grad_prim_int
    do i = 1, n_eq
      grad_prim_int_normal(i,:)  = dot_product(grad_prim_int(i,:),unit_normal(:))*unit_normal(:)
    enddo

    grad_prim_int_tangent = grad_prim_int - grad_prim_int_normal

    ! Adiabatic condition
    grad_prim_int_normal(5:5,1:3) = 0.0_wp

    ! Compute normal viscous flux arising from the ghost point
    f_v   =  normalviscousflux(prim_int,MatMul(dWdV(prim_int,n_eq), &
      & grad_prim_int_tangent+grad_prim_int_normal),normal,n_eq,mut)


    return
  end subroutine NoSlipWallAdiabatic

  !============================================================================

  subroutine SubsonicOutflow(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
   use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut

    real(wp)                 :: rho, p

    rho = Vx(1)
    p   = vx(1) * vx(5) 

    Vx(1)   =  (p0 / p)**(1.0_wp/gamma0) * rho
    Vx(2:4) = Vx(2:4) 
    Vx(5)   =   p0 / rho

    fv = 0.0_wp

    return
 end subroutine SubsonicOutflow

! subroutine UniformFreeStream_AcousticPulse(Vx,phi,fv,Jx,xin,tin,neqin,nd)
!   use nsereferencevariables
!   integer, intent(in) :: neqin, nd
!   real(wp), intent(inout) :: Vx(neqin)

!   alpha = uniformFreeStreamAOA*pi/180._wp
!   omega = 1.0_wp
!   eps   = 0.01_wp

!   p   = rho0*T0 * (1.0_wp + eps * sin(2*pi*omega*tin)
!   Vx(1)   =  (p0 / p)**(1.0_wp/gamma0) * rho0

!   Vx = zero
!   ! density
!   Vx(1) = one
!   ! velocity
!   Vx(2) = cos(alpha)
!   Vx(3) = sin(alpha)
!   Vx(4) = zero
!   ! Temperature
!   Vx(5) = one
!   
!   fv = zero

!   return

!   return
! end subroutine UniformFreeStream_AcousticPulse

  !============================================================================

  pure subroutine Constant_rhoS(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin,mut

    real(wp)             :: c1,rs,mag

    real(wp) :: alpha

    alpha = uniformFreeStreamAOA*pi/180._wp
    
    c1 = 1.0_wp
    rs = 1.0_wp
    mag = sqrt(dot_product(xin(:),xin(:)))
    vx(1)   = 1.0_wp * ( 1.0_wp + 0.1_wp*sin(mag - 0.1_wp*tin))
    Vx(2) = cos(alpha)
    Vx(3) = sin(alpha)
    vx(4)   = 0.0_wp 
    vx(5)   = exp( (-c1 - gm1og*vx(1)*log(vx(1)))/((-1 + gm1og)*vx(1)*rs) )

    fv = zero

    return
  end subroutine Constant_rhoS

  !============================================================================

  subroutine BLayerProfile(Vx,phi,fv,Jx,xin,tin,neqin,nd,mut)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin, mut
    
    Vx = Vx
    fv = zero

  end subroutine BLayerProfile

  !============================================================================
  ! compute_explicit_timestep- given a CFL, computes the max timestsp 
  !============================================================================

  subroutine compute_explicit_timestep(dt_global)

    ! Load modules
    use variables
    use controlvariables
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: pmat
    use initcollocation, only: element_properties
    use mpimod 

    implicit none

    real(wp), intent(inout) :: dt_global
    real(wp), dimension(:), allocatable :: dt_min_proc

    integer :: iell, ielh
    integer :: s_tag, r_tag, m_size, &
               s_request_dt_min, r_request_dt_min, i_err

    integer :: s_status(mpi_status_size)
    integer :: r_status(mpi_status_size)

    integer :: ielem,  inode, n_pts_1d
    integer :: i, j, k, m

    real(wp)               :: Lngth, a0, dt0, dt_min, dt_global_max
    real(wp)               :: tI, tV
    real(wp)               :: eigtot, dtN, con1
    real(wp), dimension(3) :: sq, ucon, eig


    continue

    ! Compute gradient of the velocity components
    ! ===========================================
    dt_global_max = dt_global*1.1_wp
    dt_min  = 100.0_wp

    ! Low and high  volumetric element index
    iell  = ihelems(1) ;  ielh = ihelems(2)

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_1d=n_pts_1d, pmat=pmat)

      inode = 0
      do k = 1,n_pts_1d                                                          !  z-dir

        do j = 1,n_pts_1d                                                        !  y-dir

          do i = 1,n_pts_1d                                                      !  x-dir

            inode = inode + 1                                                    !  Accumulate volume index

            Lngth = (pmat(i)*pmat(j)*pmat(k))**(third)                           !  approximate scale of computational volume

            a0  = sqrt(abs(gamma0*vg(5,inode,ielem)/gM2))                        ! Speed of sound (in strange nondimensionalization)

            sq(1) = magnitude(r_x(1,:, inode,ielem))                             !  Sqrt( d(xi_1)/dx_j . d(xi_1)/dx_j )
            sq(2) = magnitude(r_x(2,:, inode,ielem))                             !  Sqrt( d(xi_2)/dx_j . d(xi_2)/dx_j )
            sq(3) = magnitude(r_x(3,:, inode,ielem))                             !  Sqrt( d(xi_3)/dx_j . d(xi_3)/dx_j )

            ucon(1) = dot_product(r_x(1,:, inode,ielem),vg(2:4, inode,ielem))    !  d(xi_1)/dx_j . U
            ucon(2) = dot_product(r_x(2,:, inode,ielem),vg(2:4, inode,ielem))    !  d(xi_2)/dx_j . U
            ucon(3) = dot_product(r_x(3,:, inode,ielem),vg(2:4, inode,ielem))    !  d(xi_3)/dx_j . U
!    ========                                                                    !  Simplified version of timestep constraint 
            tI  =  sum(abs(ucon(:))) + a0 * ( sum(sq) )                          !  Inviscid scale

            tV  =  Re0Inv * magnitude(sq)                                        !  viscous  scale

            dt0 = CFL / (tI / Lngth + tV / Lngth / Lngth)                        !  Dt_max
!    ========                                                                    !  Directional version of timestep constraint 
            eig(:) = abs(ucon(:)) + a0 * sq(:)                                   !  Magnitude of contravariant velocity + scales speed of sound

            con1 = 4.0_wp * Re0inv / gm1M2 / 3.0_wp                              !  viscous  scaling

            eig(1) = (eig(1) + con1 * sq(1) * sq(1) / pmat(i) ) / pmat(i)        !  Inviscid + viscous scaling in xi_1 direction
            eig(2) = (eig(2) + con1 * sq(2) * sq(2) / pmat(j) ) / pmat(j)        !  Inviscid + viscous scaling in xi_2 direction
            eig(3) = (eig(3) + con1 * sq(3) * sq(3) / pmat(k) ) / pmat(k)        !  Inviscid + viscous scaling in xi_3 direction

            eigtot = magnitude(eig(:))                                           !  Magnitude of eigenvalue vector

            dtN = CFL / eigtot                                                   !  Dt_max

            dt_min  = min(dt_min,dt0)                                            !  minimum of simplified  version
            dt_min  = min(dt_min,dtN)                                            !  minimum of directional version

          enddo

        enddo

      enddo

    enddo

    if(myprocid == 0 )  then                                                     ! Reduce values on all processes to a single value
      allocate(dt_min_proc(0:nprocs-1))
      dt_min_proc(:) = 10000.0_wp ; dt_min_proc(0) = dt_min ;
    end if

    if(myprocid /= 0 ) then
      s_tag = 100 + myprocid
      m_size = 1
      call mpi_isend(dt_min,m_size,mpi_double,0,s_tag,petsc_comm_world, &
        & s_request_dt_min,i_err)
      
      call mpi_wait(s_request_dt_min,s_status,i_err)
    else
      do m = 1, nprocs-1
        r_tag = 100 + m
        m_size = 1
        call mpi_irecv(dt_min_proc(m),m_size,mpi_double,m,r_tag, &
          & petsc_comm_world,r_request_dt_min,i_err)

        call mpi_wait(r_request_dt_min,r_status,i_err)
      end do
    end if

    if(myprocid == 0) then
      dt_global = minval(dt_min_proc(:))
      deallocate(dt_min_proc)
      if(dt_global >= dt_global_max)dt_global = dt_global_max
    end if

    call mpi_barrier(petsc_comm_world,i_err)                                     ! Create a barrier synchronization in the group. 

    m_size = 1                                                                   ! Broadcast dt_global
    call mpi_bcast(dt_global,m_size,mpi_default_wp,0,petsc_comm_world,i_err)

  end subroutine compute_explicit_timestep

!============================================================================
    ! compute_vorticity_field_elements - Computes the vorticity field for all the elements
!============================================================================

  subroutine compute_vorticity_field_elements()
    ! This subroutine computes the vorticity field given the velocity field.
    ! This means that the primitive variable must be already computed.

    ! Load modules
    use variables
    use referencevariables
    use initcollocation, only: element_properties

    ! Nothing is implicitly defined
    implicit none

    real(wp), dimension(5,3) :: GradV

    integer :: ielem,  inode

    continue
    
    ! Compute gradient of the velocity components
    ! ===========================================

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      ! loop over all elements
      do  inode = 1, nodesperelem

        !  __           __
        !  \/ V  = dVdW \/ W  = dVdW phi
        !  
        GradV(:,:)    = MatMul(dVdW(vg(:, inode,ielem),nequations),phig(:,:, inode,ielem))

        ! Compute the vorticity
        omega(:, inode,ielem) = 0.0_wp 

        ! Note:  GradV(1,:) is gradient of density
        omega(1, inode,ielem) = GradV(4,2) - GradV(3,3)
        omega(2, inode,ielem) = GradV(2,3) - GradV(4,1)
        omega(3, inode,ielem) = GradV(3,1) - GradV(2,2)
      
      end do
    end do

    return
  end subroutine compute_vorticity_field_elements

  !============================================================================
  
  !============================================================================
  ! kinetic_energy_element - Computes the kinetic energy for one element

  pure function kinetic_energy_element(ielem)
    ! This subroutine computes the kinetic energy in one element given the 
    ! velocity field. This means that the primitive variable must be already computed.

    ! Load modules
    use variables
    use referencevariables

    ! Nothing is implicitly defined
    implicit none

    real(wp), dimension(nodesperelem) :: kinetic_energy_element

    integer, intent(in) :: ielem
    
    integer :: inode

    continue 

    ! Compute kinetic energy
    ! ======================
    do inode = 1, nodesperelem

      kinetic_energy_element(inode) = 0.0_wp 

      kinetic_energy_element(inode) = 0.5_wp*dot_product(vg(2:4,inode,ielem),vg(2:4,inode,ielem))
      
    end do

    return
  end function kinetic_energy_element

  !============================================================================
  ! compute_primitive_variables_elements - Computes the primitive variables for
  !                                        all elements using the conservative 
  !                                        variables. 

  subroutine compute_primitive_variables_elements()

    ! Load modules
    use variables
    use referencevariables
    use initcollocation, only: element_properties

    ! Nothing is defined implicitly
    implicit none

    integer :: ielem, inode

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      ! Loop over nodes in each element
      do inode = 1, nodesperelem

        ! Calculate primitive variables from conservative variables
        call conserved_to_primitive(ug(:,inode,ielem),vg(:,inode,ielem),nequations)
      
      enddo
    enddo

    return
  end subroutine compute_primitive_variables_elements

  !============================================================================

  !============================================================================
  ! compute_entropy_variables_elements - Computes the entropy variables using 
  !                                      for all elements using the primitive 
  !                                      variables. 

  subroutine compute_entropy_variables_elements()
 
    ! Load modules
    use variables
    use referencevariables
    use initcollocation, only: element_properties

    ! Nothing is defined implicitly
    implicit none

    integer :: ielem, inode

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      ! Loop over nodes in each element
      do inode = 1, nodesperelem
            
        ! Calculate entropy variables from primitive variables
        call primitive_to_entropy( &
        vin = vg(:,inode,ielem), &
        wout = wg(:,inode,ielem), &
        nq = nequations )

      enddo
    enddo

    return
  end subroutine compute_entropy_variables_elements

  !============================================================================
  ! compute_specific_entropy_elements - Computes the specific entropy using the 
  ! primitive variables 

  subroutine compute_specific_entropy_elements()
 
    ! Load modules
    use variables
    use referencevariables
    use initcollocation, only: element_properties

    ! Nothing is defined implicitly
    implicit none

    integer :: ielem, inode

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      ! Loop over nodes in each element
      do inode = 1, nodesperelem
            
        ! Calculate specific entropy
        specific_entropy(inode,ielem) = specificentropy(vg(:,inode,ielem),nequations)

      enddo
    enddo

    return
  end subroutine compute_specific_entropy_elements


  subroutine Solution_Filter()
    ! This subroutine calculates elementwise 
    ! the Divergence of the Conservative Flux
    use variables
    use referencevariables
    use controlvariables
    use collocationvariables, only: ia_Filter,ja_Filter,aa_Filter
    use initcollocation, only: element_properties

    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)

    ! indices
    integer :: inode, jdir, ielem, nodesperelem_max
    integer :: jnode
    integer :: i

    real(wp),  parameter :: al = 0.01_wp
    real(wp),  allocatable, dimension(:,:) :: tmp88,tmp89
    real(wp)             :: one_al

    one_al = 1.0_wp - al

    nodesperelem_max = (npoly_max+1)**ndim
    
    allocate(tmp88(1:nequations,1:nodesperelem_max))
    allocate(tmp89(1:nequations,1:nodesperelem_max))

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem, n_pts_3d=nodesperelem)

      ! loop over all nodes in the element
      tmp89(:,1:nodesperelem) = ug(:,1:nodesperelem,ielem)
      ! loop over each direction
      do jdir = 1,ndim
        tmp88(:,:) = 0.0_wp
        do inode = 1,nodesperelem
          ! column/node from gradient operator in CSR format in
          do i = ia_Filter(inode), ia_Filter(inode+1)-1
            jnode = ja_Filter(jdir,i)
            tmp88(:,inode) = tmp88(:,inode) + aa_Filter(jdir,i) * tmp89(:,jnode)
          end do
        end do
        tmp89(:,1:nodesperelem) = tmp88(:,1:nodesperelem)
      end do
      ug(:,1:nodesperelem,ielem) = one_al * ug(:,1:nodesperelem,ielem)  + al*tmp88(:,1:nodesperelem)
    end do
    deallocate(tmp88) ; deallocate(tmp89)

    return
  end subroutine Solution_Filter
 
!============================================================================

  subroutine Flux_Div_Pencil(ielem, N_S, N_S_2d, pinv, qmat, dmat)

    ! This subroutine calculates elementwise the Divergence of the Conservative Flux

    use variables
    use referencevariables
    use controlvariables, only: discretization, Entropy_Correction
    use SSWENOvariables
    use initgrid

    implicit none

    integer ,                     intent(in)   :: ielem, N_S, N_S_2d

    real(wp), dimension(N_S),     intent(in)   :: pinv
    real(wp), dimension(N_S,N_S), intent(in)   :: qmat, dmat

    ! indices
    integer :: inode, jdir, ipen
    integer :: i, k

    real(wp), dimension(nequations,N_S)        :: ugS, vgS, fgS, d_fgS
    real(wp), dimension(         3,N_S)        :: JnS

    integer,  dimension(2)                     :: faceLR
    real(wp), dimension(nequations)            :: uLL,uL,uR,uRR
    real(wp), parameter                        :: tol1 = 1.0e-10_wp

    integer                                    :: jnode,gnode
    integer                                    :: k_node,k_elem,k_face
    integer                                    :: inb


    select case(discretization)

    case('SpecColl')

      do jdir = 1,ndim            ! Directional loop

        do ipen = 1, N_S_2d

          !  Grab a pencil of data
          fgS(:,:) = zero
          do i = 1,N_S
            inode    = Pencil_Coord(N_S,jdir,ipen,i)
            ugS(:,i) =   ug(     :,inode,ielem)    !  Neqns
            JnS(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)

            call conserved_to_primitive(ugS(:,i),vgS(:,i),nequations)
            fgS(:,i) = normalflux( vgS(:,i), JnS(:,i), nequations )

          enddo

          ! Differentiate Pencil of Fluxes on to solution points
          do k = 1,nequations
            d_fgS(k,1:N_S) = matmul(dmat,fgS(k,1:N_S))
          enddo

          do i = 1,N_S
            inode   = Pencil_Coord(N_S,jdir,ipen,i)
            divf(:,jdir,inode,ielem) = d_fgS(:,i)
          enddo

        end do

      end do

    case('SSDC')

      do jdir = 1,ndim            ! Directional loop

        do ipen = 1, N_S_2d

          !  Grab a pencil of data
          do i = 1,N_S
            inode    = Pencil_Coord(N_S,jdir,ipen,i)  
            ugS(:,i) =   ug(     :,inode,ielem)    !  Neqns
            JnS(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
          enddo

          !  Extrapolate data from Solution points to Flux points

          if( .not. Entropy_Correction) then
            call SS_Euler_Dspec(N_S, nequations, pinv, qmat, dmat, ugS, JnS, d_fgS)
          else
            call SS_Stabilized_Euler_Dspec(N_S,nequations, pinv, qmat, ugS, JnS, d_fgS)
          endif

          do i = 1,N_S
            inode   = Pencil_Coord(N_S,jdir,ipen,i)
            divf(:,jdir,inode,ielem) = d_fgS(:,i)
          enddo

        end do

      end do

    case('SSWENO')

      do jdir = 1,ndim            ! Directional loop

        faceLR = face_pairs(jdir)

        do ipen = 1,N_S_2d

           inb  =  0     !  Assumes the pencil is an interior one.

           !  Left interface partner nodes

           if     (ef2e(1,faceLR(1),ielem) < 0) then            ! Partner is a BC (i.e. no partner) 
             inb   = -1
             inode = Pencil_Coord(N_S,jdir,ipen,1)  
             uL(:) = ug(:, inode,ielem)
            uLL(:) = ug(:, inode,ielem)

          else if (ef2e(3,faceLR(1),ielem) /= myprocid) then    ! Partner is off-process
            jnode  = nodesperface*(faceLR(1)-1)+ipen
            gnode  = efn2efn(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack

             uL(:) = ughst(:,gnode)
!           uLL(:) = ughstWENO(:,gnode)
            uLL(:) = ughstWENO_partner(:,gnode)

          else                                              ! Partner is  on-process
             call data_partner_element_serial(ipen,faceLR(1),ielem,k_node,k_elem,k_face, inode)
            gnode  = nodesperface*(k_face-1) + ipen

             uL(:) = ug(:,k_node,k_elem)
!           uLL(:) = ug(:,WENO_Adjoining_Data(k_node,k_face),k_elem)
            uLL(:) = ugWENO_partner(:,gnode,k_elem)

!           t1 = abs(maxval(uLLT(:)-uLL(:)))
!           if(t1 >= tol1) write(*,*)t1
          endif

!          !  Right interface partner nodes

          if     (ef2e(1,faceLR(2),ielem) < 0) then            ! Partner is a BC (no partner) 
             inb   = +1
             inode = Pencil_Coord(N_S,jdir,ipen,N_S)  
             uR(:) = ug(:, inode,ielem)
            uRR(:) = ug(:, inode,ielem)

          else if (ef2e(3,faceLR(2),ielem) /= myprocid) then    ! Partner is off-process
            jnode  = nodesperface*(faceLR(2)-1)+ipen
            gnode  = efn2efn(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack

             uR(:) = ughst(:,gnode)
!           uRR(:) = ughstWENO(:,gnode)
            uRR(:) = ughstWENO_partner(:,gnode)

          else                                              ! Partner is  on-process
             call data_partner_element_serial(ipen,faceLR(2),ielem,k_node,k_elem,k_face, inode)
            gnode  = nodesperface*(k_face-1) + ipen

             uR(:) = ug(:,k_node,k_elem)
!           uRR(:) = ug(:,WENO_Adjoining_Data(k_node,k_face),k_elem)
            uRR(:) = ugWENO_partner(:,gnode,k_elem)
!           t1 = abs(maxval(uRRT(:)-uRR(:)))
!           if(t1 >= tol1) write(*,*)t1
          endif

          !  Grab a pencil of data
          do i = 1,N_S
            inode    = Pencil_Coord(N_S,jdir,ipen,i)  
            ugS(:,i) =   ug(     :,inode,ielem)    !  Neqns
            JnS(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
          enddo
             
          ! First  conditional catches elements that have BC on both faces
          ! Second conditional catches elements that are periodic onto themselves (1 wide)
          if ( (ef2e(1,faceLR(1),ielem) <      0) .and. (ef2e(1,faceLR(2),ielem) <      0)    .or. &
             & (ef2e(2,faceLR(1),ielem) == ielem) .and. (ef2e(2,faceLR(2),ielem) == ielem)  ) then
            call SS_Stabilized_Euler_Dspec(N_S,nequations, pinv, qmat, ugS, JnS, d_fgS)
          else
            if(WENO_type == 'Element_WENO') then
              uLL = uL ; uRR = uR ;     !  Ensures consistency of smoothness parameters tau_i
            endif
            call SSWENO4(inb, ugS, uLL, uL, uR, uRR, WENO_Extrp, 1.0_wp, WENO_Extrp, JnS, d_fgS)
          endif

          do i = 1,N_S
            inode   = Pencil_Coord(N_S,jdir,ipen,i)
            divf(:,jdir,inode,ielem) = d_fgS(:,i)
          enddo

        enddo

      enddo

    end select

    return
  end subroutine Flux_Div_Pencil

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 
  subroutine Entropy_Inviscid_Flux_Div(ielem)
    ! This subroutine calculates elementwise 
    ! the Divergence of the Conservative Flux
    use variables
    use referencevariables
    use collocationvariables
    use initgrid

    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem

    ! indices
    integer :: inode, jnode, jdir
    integer :: i

    real(wp), dimension(nequations) :: ugS, vgS
    real(wp), dimension(         3) :: JnS

        do inode = 1, nodesperelem
          divf_S(:,inode,ielem) = 0.0_wp
          ! loop over number of dependent elements in gradient
          do i = iagrad(inode), iagrad(inode+1)-1
            ! loop over dimensions
            do jdir = 1,ndim
              ! column/node from gradient operator in CSR format in
              ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
              jnode = jagrad(jdir,i)
              ! update gradient using coefficient and entropy variables at appropriate node
              JnS(:) = Jx_r(jnode,ielem)*r_x(jdir,:,jnode,ielem)
              ugS(:) =  ug(:,jnode,ielem)
              call conserved_to_primitive(ugS(:),vgS(:),nequations)
              divf_S(jdir,inode,ielem) = divf_S(jdir,inode,ielem) + dagrad(jdir,i) * normal_Entropy_flux(vgS(:),JnS(:),nequations)
            end do
          end do
        end do

  end subroutine Entropy_Inviscid_Flux_Div

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine SS_Euler_Dspec(N_S, N_q, pinv, qmat, dmat, ugS, JnS, dfn)

    use variables
    use referencevariables
    use controlvariables, only: entropy_flux

    implicit none 

    integer,                      intent(in ) :: N_S, N_q
    real(wp), dimension(N_S),     intent(in ) :: pinv
    real(wp), dimension(N_S,N_S), intent(in ) :: qmat, dmat
    real(wp), dimension(N_q,N_S), intent(in ) :: ugS
    real(wp), dimension(  3,N_S), intent(in ) :: JnS

    real(wp), dimension(N_q,N_S), intent(out) :: dfn

    real(wp), dimension(N_q,N_S,N_S)          :: SSFlux

    real(wp), dimension(3)                    :: nx

    integer                                   :: i,j,k,l

    real(wp), dimension(N_q,N_S)              :: vgS
    real(wp), dimension(N_q,0:N_S)            :: fnS

    logical                                   :: flux = .false.

    !  Rotate pencil from conserved variables to primitive variables
    do i=1,N_S
      call conserved_to_primitive(ugS(:,i),vgS(:,i),N_q)
    enddo

    SSFlux(:,:,:) = 0.0_wp ;

    !  EntropyConsistentFlux is symmetric w.r.t. the Left and Right states
    !  Only the upper half of the matrix are calculated
    !  Interior Diagonals are not calculated because d[[i,i]] = 0 for i /= 1, i/= N_F

    if(.not. flux) then

      select case(entropy_flux)
        case('Ismail_Roe')
          do i=1,N_S-1
              do j=i+1,N_S
                        nx(:) = 0.5_wp *(JnS(:,i) + JnS(:,j))
                SSFlux(:,i,j) = 2.0_wp * EntropyConsistentFlux     (vgS(:,i),vgS(:,j),nx,N_q)
                SSFlux(:,j,i) = SSFlux(:,i,j)
            enddo
          enddo
        case('Chandrashekar') 
          do i=1,N_S-1
              do j=i+1,N_S
                        nx(:) = 0.5_wp *(JnS(:,i) + JnS(:,j))
                SSFlux(:,i,j) = 2.0_wp * Entropy_KE_Consistent_Flux(vgS(:,i),vgS(:,j),nx,N_q)
                SSFlux(:,j,i) = SSFlux(:,i,j)
            enddo
          enddo
      end select

      SSFlux(:,  1,  1) = 2.0_wp * normalflux( vgS(:,  1), JnS(:,  1), N_q )
      SSFlux(:,N_S,N_S) = 2.0_wp * normalflux( vgS(:,N_S), JnS(:,N_S), N_q )

      do i=1,N_S
        dfn(:,i) = zero
        do j=1,N_S
          dfn(:,i) = dfn(:,i) + dmat(i,j)*SSFlux(:,i,j)
        enddo
      enddo

    else

      select case(entropy_flux)
        case('Ismail_Roe')
          do i=1,N_S-1
            do j=i+1,N_S
                      nx(:) = 0.5_wp *(JnS(:,i) + JnS(:,j))
              SSFlux(:,i,j) = 2.0_wp * qmat(i,j)*EntropyConsistentFlux     (vgS(:,i),vgS(:,j),nx,N_q) 
            enddo
  
            fnS(:,i) = 0.0_wp
            do k=i+1,N_S
              do l=1,i
                fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
              enddo
            enddo
          enddo
        case('Chandrashekar') 
          do i=1,N_S-1
            do j=i+1,N_S
                      nx(:) = 0.5_wp *(JnS(:,i) + JnS(:,j))
              SSFlux(:,i,j) = 2.0_wp * qmat(i,j)*Entropy_KE_Consistent_Flux(vgS(:,i),vgS(:,j),nx,N_q)
            enddo
  
            fnS(:,i) = 0.0_wp
            do k=i+1,N_S
              do l=1,i
                fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
              enddo
            enddo
          enddo
      end select

      fnS(:,  0) = normalflux  (vgS(:,  1),JnS(:,  1),N_q)
      fnS(:,N_S) = normalflux  (vgS(:,N_S),JnS(:,N_S),N_q)

      do k = 1,N_S
        dfn(:,k) = pinv(k) * (fnS(:,k) - fnS(:,k-1))
      enddo

    endif

  end subroutine SS_Euler_Dspec

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine SS_Stabilized_Euler_Dspec(N_S, N_q, pinv, qmat, ugS,JnS, dfn)

    use variables
    use referencevariables
    use controlvariables, only: flux_entropy_correction, entropy_flux

    implicit none 

    integer,                      intent(in ) :: N_S, N_q
    real(wp), dimension(N_S),     intent(in ) :: pinv
    real(wp), dimension(N_S,N_S), intent(in ) :: qmat
    real(wp), dimension(N_q,N_S), intent(in ) :: ugS
    real(wp), dimension(  3,N_S), intent(in ) :: JnS

    real(wp), dimension(N_q,N_S), intent(out) :: dfn

    real(wp), dimension(N_q,N_S,N_S)          :: SSFlux, DSFlux

    real(wp), dimension(N_q,1:N_S)            :: vgS, wgS

    real(wp), dimension(N_q,0:N_S)            :: fnS, fnD

    real(wp), dimension(3)                    :: nx

    real(wp), parameter                       :: cc2  = 1.0e-24_wp

    integer                                   :: i,j,k,l

    real(wp)                                  :: bbS,dsS,deltaS

    !  Rotate pencil from conserved variables to primitive variables
    do i=1,N_S
      call conserved_to_primitive(ugS(:,i),vgS(:,i),N_q)
      call primitive_to_entropy  (vgS(:,i),wgS(:,i),N_q)
    enddo

    !  Form Entropy Fluxes 0:N_S
       fnS(:,:)   = 0.0_wp ;    fnD(:,:)   = 0.0_wp ;
    SSFlux(:,:,:) = 0.0_wp ; DSFlux(:,:,:) = 0.0_wp ;

    select case(entropy_flux)

      case('Ismail_Roe')

        do i=1,N_S-1
          do j=i+1,N_S
                    nx(:) = 0.5_wp * (JnS(:,i) + JnS(:,j))
            SSFlux(:,i,j) = 2.0_wp * qmat(i,j)*EntropyConsistentFlux     (vgS(:,i),vgS(:,j),nx,N_q) 
            
            select case (flux_entropy_correction)
              case ('normal')
                DSFlux(:,i,j) = 2.0_wp * qmat(i,j)*normalflux    (0.5_wp*(vgS(:,i)+vgS(:,j)),nx,N_q)
              case ('Honein-Moin')
                DSFlux(:,i,j) = 2.0_wp * qmat(i,j)*HoneinMoinFlux(vgS(:,i),vgS(:,j),nx,N_q)
              case default
                write(*,*) 'The flux selected to be use for the entropy correction is unknown.'
                write(*,*) 'Check the subroutine SS_Stabilized_Euler_Dspec()'
                write(*,*) 'Exting....'
                stop
            end select
          enddo

          fnS(:,i) = 0.0_wp ; fnD(:,i) = 0.0_wp ;
          do k=i+1,N_S
            do l=1,i
              fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
              fnD(:,i) = fnD(:,i) + DSFlux(:,l,k)
            enddo
          enddo
        enddo

      case('Chandrashekar')

        do i=1,N_S-1
          do j=i+1,N_S
                    nx(:) = 0.5_wp * (JnS(:,i) + JnS(:,j))
            SSFlux(:,i,j) = 2.0_wp * qmat(i,j)*Entropy_KE_Consistent_Flux(vgS(:,i),vgS(:,j),nx,N_q)
            
            select case (flux_entropy_correction)
              case ('normal')
                DSFlux(:,i,j) = 2.0_wp * qmat(i,j)*normalflux    (0.5_wp*(vgS(:,i)+vgS(:,j)),nx,N_q)
              case ('Honein-Moin')
                DSFlux(:,i,j) = 2.0_wp * qmat(i,j)*HoneinMoinFlux(vgS(:,i),vgS(:,j),nx,N_q)
              case default
                write(*,*) 'The flux selected to be use for the entropy correction is unknown.'
                write(*,*) 'Check the subroutine SS_Stabilized_Euler_Dspec()'
                write(*,*) 'Exting....'
                stop
            end select
          enddo

          fnS(:,i) = 0.0_wp ; fnD(:,i) = 0.0_wp ;
          do k=i+1,N_S
            do l=1,i
              fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
              fnD(:,i) = fnD(:,i) + DSFlux(:,l,k)
            enddo
          enddo
        enddo
      case default 

    end select

    fnS(:,  0) = normalflux(vgS(:,  1),JnS(:,  1),N_q) ; fnD(:,  0) = fnS(:,  0)
    fnS(:,N_S) = normalflux(vgS(:,N_S),JnS(:,N_S),N_q) ; fnD(:,N_S) = fnS(:,N_S)

    !  Entropy Correction 
    do i = 1,N_S-1
!         bb(:) = (wgS(:,i+1) - wgS(:,i+0))*(fnS(:,i) - fnD(:,i))
!         ds(:) = sqrt(bb(:)*bb(:) + cc2)
!      delta(:) = (ds(:) - bb(:))/(two*ds(:))
!      fnD(:,i) = fnD(:,i) + delta(:)*(fnS(:,i) - fnD(:,i))  !  Entropy Correction
         bbS = dot_product(wgS(:,i+1)-wgS(:,i+0),fnS(:,i)-fnD(:,i))
         dsS = sqrt(bbS*bbS + cc2)
      deltaS = (dsS - bbS)/(two*dsS)
      fnD(:,i) = fnD(:,i) + deltaS * (fnS(:,i) - fnD(:,i))
    enddo

    do i = 1,N_S
      dfn(:,i) = pinv(i) * (fnD(:,i) - fnD(:,i-1))
    enddo

  end subroutine SS_Stabilized_Euler_Dspec

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !============================================================================
  ! matrix_hatc_node - Computes the matrix \hat{[C]} of one node using the 
  ! computational space coordinates.
  !
  ! Input parameters:
  ! v_in  - primitive variables at the node.
  ! n_i - divergence contravariant vector.
  ! n_j - gradient contravariant vector.
  ! n_eq - number of equations. 
  !
  ! Output parameter:
  ! matrix_hatc_node - matrix \hat{C} of the node.

  pure function matrix_hatc_node(v_in,n_i,n_j,n_eq)
    
    ! Load modules
    use nsereferencevariables, only: gm1M2, Re0inv, Pr0, k0

    use controlvariables, only : variable_viscosity
      
    ! Nothing is implicitly defined
    implicit none
    
    integer,  intent(in) :: n_eq
    real(wp), intent(in) :: n_i(3), n_j(3)
    real(wp), intent(in) :: v_in(n_eq)
    real(wp), dimension(n_eq,n_eq) :: mat
    real(wp) :: con1, con2, con3, con4
    real(wp) :: u, v, w, T
    real(wp) :: mu

    real(wp) :: matrix_hatc_node(n_eq,n_eq)

    ! Initialize all elements of mat to zero
    mat = 0.0_wp

    ! Dereference variables
    u = v_in(2)
    v = v_in(3) 
    w = v_in(4) 
    T = v_in(5)

    ! Set dynamic viscosity
    if (variable_viscosity .eqv. .true.) then
      mu = sutherland_law(T)
    else
      mu = 1.0_wp
    end if

    ! Define some useful constants
    con1 = Re0inv * mu * T / gm1M2 / 3.0_wp
    con2 = gm1M2
    con3 = gm1M2
    con4 = Re0inv * k0 * T * T / Pr0 

    ! Momentum equation
    mat(2,2) = con1*( 4.0_wp*n_i(1)*n_j(1) + 3.0_wp*n_i(2)*n_j(2) + 3.0_wp*n_i(3)*n_j(3))
    mat(2,3) = con1*( 3.0_wp*n_i(2)*n_j(1) - 2.0_wp*n_i(1)*n_j(2)                )
    mat(2,4) = con1*( 3.0_wp*n_i(3)*n_j(1)                        - 2.0_wp*n_i(1)*n_j(3))

    mat(3,2) = con1*(-2.0_wp*n_i(2)*n_j(1) + 3.0_wp*n_i(1)*n_j(2)                )
    mat(3,3) = con1*( 3.0_wp*n_i(1)*n_j(1) + 4.0_wp*n_i(2)*n_j(2) + 3.0_wp*n_i(3)*n_j(3))
    mat(3,4) = con1*(                        3.0_wp*n_i(3)*n_j(2) - 2.0_wp*n_i(2)*n_j(3))

    mat(4,2) = con1*(-2.0_wp*n_i(3)*n_j(1)                        + 3.0_wp*n_i(1)*n_j(3))
    mat(4,3) = con1*(                      - 2.0_wp*n_i(3)*n_j(2) + 3.0_wp*n_i(2)*n_j(3))
    mat(4,4) = con1*( 3.0_wp*n_i(1)*n_j(1) + 3.0_wp*n_i(2)*n_j(2) + 4.0_wp*n_i(3)*n_j(3))

    mat(2,5) = con2*(mat(2,2)*u + mat(2,3)*v + mat(2,4)*w) 
    mat(3,5) = con2*(mat(3,2)*u + mat(3,3)*v + mat(3,4)*w)
    mat(4,5) = con2*(mat(4,2)*u + mat(4,3)*v + mat(4,4)*w)

    ! Energy equation
    mat(5,2) = con2*(u*mat(2,2) + v*mat(3,2) + w*mat(4,2)) 
    mat(5,3) = con2*(u*mat(2,3) + v*mat(3,3) + w*mat(4,3))
    mat(5,4) = con2*(u*mat(2,4) + v*mat(3,4) + w*mat(4,4))

    mat(5,5) = con3*(mat(5,2)*u + mat(5,3)*v + mat(5,4)*w) &
             + con4*(n_i(1)*n_j(1) + n_i(2)*n_j(2) + n_i(3)*n_j(3)) 

    ! Jacobian matrix wrt to the conserved variables 
    matrix_hatc_node = mat(:,:)

    return
  end function matrix_hatc_node

  !============================================================================
  ! set_boundary_conditions - Set the pointer to the boundary condition subroutine. 

  subroutine set_boundary_conditions(bc_type,initial_condition)
    
    ! Nothing is implicitly defined
    implicit none

    integer,                     intent(in)    :: bc_type
    character(180),              intent(in)    :: initial_condition
    
    continue

    ! Call the the BC function specified by bc_type
    select case(bc_type)
      case(-2)  !  BCDirichlet handles ALL the exact solution data
        if( initial_condition  == 'ExactSolutionIsentropicVortex') then
          BoundaryCondition => isentropicVortexFull

        else if(initial_condition  == 'ExactSolutionSupersonicVortex') then
          BoundaryCondition => supersonicvortexFull

        else if(initial_condition  == 'ShockVortex') then
          BoundaryCondition => ShockVortexInteraction

        else if(initial_condition  == 'ExactSolutionViscousShock') then
          BoundaryCondition => viscousShockFull
        
        else if(initial_condition  == 'Constant_rhoS') then
          BoundaryCondition => Constant_rhoS

        else if(initial_condition  == 'PreserveFreeStream') then
          BoundaryCondition => UniformFreeStream

        else if(initial_condition  == 'SodsProblem') then
          BoundaryCondition => SodsProblemICBC

        else if(initial_condition == 'Potential_cylinder') then
          BoundaryCondition => Potential_Flow_Around_cylinder

        endif
      case(-3)
        BoundaryCondition => UniformFreeStream
      case(-4)
        BoundaryCondition => InviscidWall
      !case(-5)
      !  BoundaryCondition => no_penetration_BC
      case(-5)
        BoundaryCondition => NoSlipWallAdiabatic 
      case(-7)
        BoundaryCondition => SymmetryPlane
!     case(-8)
!       BoundaryCondition => Periodic_1
!     case(-9)
!       BoundaryCondition => Periodic_2
      case(-10)
        BoundaryCondition => BLayerProfile
      case(-11)
        BoundaryCondition => SubsonicOutflow
      case(-16)

        if( initial_condition  == 'ExactSolutionIsentropicVortex') then
          BoundaryCondition => isentropicVortexFull

        else if(initial_condition  == 'ExactSolutionViscousShock') then
          BoundaryCondition => viscousShockFull

        else if(initial_condition  == 'ExactSolutionSupersonicVortex') then
          BoundaryCondition => supersonicvortexFull

        else if(initial_condition  == 'PreserveFreeStream') then
          BoundaryCondition => UniformFreeStream

        endif

      case default
        BoundaryCondition => UniformFreeStream
    end select
    
    return
  end subroutine set_boundary_conditions
  
  !============================================================================

  pure function sutherland_law(T_in)

    ! Load modules
    use nsereferencevariables, only: Tref

    ! Nothing is implicitly defined
    implicit none

    real(wp) :: sutherland_law
    real(wp), intent(in) :: T_in 
    real(wp) :: s_c

    continue

    ! Sutherland constant divided by the free stream reference temperature in
    ! Rankine
!   s_c = 198.6_wp/(Tref*9.0_wp/5.0_wp)
    s_c = 198.6_wp/(Tref*1.8_wp)

    ! Compute new non-dimensional dynamic viscosity
!   sutherland_law = T_in**(3.0_wp/2.0_wp)*(1.0_wp + s_c)/(T_in + s_c)
    sutherland_law = sqrt(T_in*T_in*T_in)*(1.0_wp + s_c)/(T_in + s_c)

    return
  end function

  !============================================================================

  subroutine compute_gradient_entropy_variables()
    
    ! Load modules
    use variables
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: nnzgrad, iagrad, jagrad, dagrad
    use initcollocation,      only: element_properties
    
    ! Nothing is implicitly defined
    implicit none

    ! loop indices
    integer :: ielem, jdir, idir
    integer :: inode, jnode
    integer :: i

    ! Temporary arrays for phi
    real(wp), allocatable :: phitmp(:,:)

    ! phitmp is calculated in computational space
    allocate(phitmp(nequations,ndim))

    ! loop over all elements
    do ielem = ihelems(1), ihelems(2)

      call element_properties(ielem,        &
                      n_pts_3d=nodesperelem,&
                       nnzgrad=nnzgrad,     &
                        iagrad=iagrad,      &
                        jagrad=jagrad,      &
                        dagrad=dagrad) 

    ! compute computational gradients of the entropy variables
    ! initialize phi
    phig(:,:,:,ielem) = 0.0_wp
    ! loop over every node in element
    do inode = 1, nodesperelem
      ! reinitialize computational gradient to zero
      phitmp(:,:) = 0.0_wp
      ! loop over number of dependent elements in gradient
      do i = iagrad(inode), iagrad(inode+1)-1
        ! loop over dimensions
        do jdir = 1,ndim
          ! column/node from gradient operator in CSR format in
          ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
          jnode = jagrad(jdir,i)
          ! update gradient using coefficient and entropy variables at appropriate node
          phitmp(:,jdir) = phitmp(:,jdir) + dagrad(jdir,i)*wg(:,jnode,ielem) 
        end do
      end do
      ! transform to physical space using dxi_jdir/dx_idir
      do jdir = 1,ndim
        do idir = 1,ndim
          phig(:,idir,inode,ielem) = phig(:,idir,inode,ielem) + phitmp(:,jdir)*r_x(jdir,idir,inode,ielem)
        end do
      end do

    end do
      
    end do
    deallocate(phitmp)

    return
  end subroutine compute_gradient_entropy_variables

!============================================================================

      SUBROUTINE SSWENO4(inb, uint, uLL, uL, uR, uRR, dxL, dx, dxR, nxint, dfbar)

      use SSWENOvariables
      use nsereferencevariables
      use precision_vars
      use collocationvariables
      use controlvariables, only: WENO_Bias, entropy_flux
      use interpolation

      implicit none
      integer, parameter   :: nq = 5, ixd = 4

      real(wp), parameter :: sqrt5  = sqrt(5.0_wp)
      real(wp), parameter :: theta  = 3.0_wp

      real(wp), parameter :: c1     = 15.0_wp/4.0_wp
      real(wp), parameter :: c2     = 15.0_wp/4.0_wp*sqrt5
      real(wp), parameter :: x66F   = 1.0_wp * (3.0_wp/178.0_wp*(19.0_wp+sqrt5))! +0.357911_wp  (Original value)
      real(wp), parameter :: x66E   = 2.0_wp * (3.0_wp/178.0_wp*(19.0_wp+sqrt5))! +0.715823_wp  (Creates an Element-wise WENO)

      real(wp), dimension(4), parameter :: tau_coeff   = (/-c1,c2,-c2,c1/)
      real(wp), dimension(5), parameter :: tau_coeffT0 = (/0.0_wp,-c1,c2,-c2,c1/)
      real(wp), dimension(5), parameter :: tau_coeffT1 = (/-c1,c2,-c2,c1,0.0_wp/)

      real(wp), dimension(5)   :: tau_cL4,tau_cR4

      integer,                       intent(in)    :: inb
      real(wp), dimension(nq, 4),    intent(in)    :: uint
      real(wp), dimension(nq),       intent(in)    :: uLL, uL, uR, uRR 
      real(wp),                      intent(in)    :: dxL,dx,dxR
      real(wp), dimension( 3,1:4),   intent(in)    :: nxint
      real(wp), dimension(nq,1:4),   intent(inout) :: dfbar

      real(wp), dimension(3)           :: nLL, nL, nR, nRR
      real(wp), dimension(3)           :: nxL,nxR, nxave
      real(wp), dimension(nq,1:4)      :: vint, wint
      real(wp), dimension(nq, 6)       :: uin, qin, vin
      real(wp), dimension(3,  6)       :: nin
      real(wp), dimension(6)           :: uhat, metr, cav2
      real(wp), dimension(nq,ixd,ixd)  :: SSFlux!, DSFlux
      real(wp), dimension(nq,0:ixd)    :: fbarW, fbarC, fnS!, fnK

      ! Interpolation coefs and target values
      real(wp), dimension(5,6)         :: Tar
      real(wp), dimension(ixd+1,2)     :: IM2n, IM1n, IP1n, IP2n

      ! Beta's - smoothness indicators
      real(wp), dimension(nq,6)        :: beta
      real(wp), dimension(nq)          :: betam, tau, betad

      ! split weights
      real(wp), dimension(nq,6)        :: wpos, wneg, tbratp1, weights
      real(wp), dimension(nq,6)        :: fcp, fcm, fh, fhC
      real(wp), dimension(6)           :: dh, dhp, dhm
      real(wp), dimension(nq)          :: wpsum, wmsum
      real(wp)                         :: lambda, eps, Lsigp, Lsigm

      real(wp), dimension(nq,nq)       :: Smat,Sinv
      real(wp), dimension(nq)          :: ev, vL, vR, vav

      real(wp), dimension(nq)          :: bb ,ds ,delta 
      real(wp), parameter              :: cc2 = 1.0e-24_wp
      real(wp), parameter              :: cc1 = 1.0e-12_wp

      real(wp)                         :: x66

      integer                          :: i,j,k,l,n
      integer                          :: iL, iH

      if    (WENO_type == 'Neighbr_WENO') then
        x66 = x66F     ! +0.357911_wp  (Full WENO involving neighbors)
      elseif(WENO_type == 'Element_WENO') then
        x66 = x66E     ! +0.715823_wp  (Creates an Element-wise WENO)
      else
        write(*,*)'Invalid WENO_type'
        write(*,*)'running Element_WENO'
        x66 = x66E
      endif

      nLL(:) = Poly3_Intrp(3,nxint,-1.0_wp,+1.0_wp,-2.0_wp+1.0_wp/sqrt5)
       nL(:) = nxint(:,  1)
       nR(:) = nxint(:,ixd)
      nRR(:) = Poly3_Intrp(3,nxint,-1.0_wp,+1.0_wp,+2.0_wp-1.0_wp/sqrt5)

      IP1n(4,1) =        + ((5.0_wp + sqrt(5.0_wp))*dx/(24.0_wp*dxL))
      IP1n(4,2) = 1.0_wp - ((5.0_wp + sqrt(5.0_wp))*dx/(24.0_wp*dxL))
      IP2n(4,1) = IP1n(4,1)
      IP2n(4,2) = IP1n(4,2)

      IM1n(2,1) = 1.0_wp - ((5.0_wp + sqrt(5.0_wp))*dx/(24.0_wp*dxR))
      IM1n(2,2) =        + ((5.0_wp + sqrt(5.0_wp))*dx/(24.0_wp*dxR))
      IM2n(2,1) = IM1n(2,1)
      IM2n(2,2) = IM1n(2,2)

      Tar = zero

      if(inb.eq.-1)then
        Tar(2,1)= 3.0_wp/89.0_wp *(19.0_wp + sqrt5);
        Tar(2,2)= 1.0_wp/178.0_wp *(95.0_wp - 17.0_wp*sqrt5);
        Tar(2,3)= -2.0_wp/(31.0_wp+11.0_wp*sqrt5);
        Tar(2,4:6) = zero
      else
        Tar(2,1)= dxL*(-3.0_wp*(-5+sqrt(5.0_wp))+(-25.0_wp+6.0_wp*sqrt(5.0_wp))*x66)/  &
                  (6.0_wp*(-5.0_wp+sqrt(5.0_wp))*dxL - 5.0_wp*dx);
        Tar(2,2)= -Tar(2,1);
        Tar(2,3)= 3.0_wp/89.0_wp *(19.0_wp + sqrt5) - x66;
        Tar(2,4)= x66
        Tar(2,5)= 1.0_wp/178.0_wp *(95.0_wp - 17.0_wp*sqrt5);
        Tar(2,6)= 1.0_wp/178.0_wp *(-31.0_wp+11.0_wp*sqrt5);
      endif

        Tar(3,1)=  1.0_wp/6.0_wp *(8.0_wp - 3.0_wp*sqrt5);
        Tar(3,2)= -5.0_wp/3.0_wp+sqrt5;
        Tar(3,3)=  1.0_wp/6.0_wp *(8.0_wp - 3.0_wp*sqrt5);
        Tar(3,4:6) = zero

      if(inb.eq.+1)then
        Tar(4,1)= -2.0_wp/(31.0_wp+11.0_wp*sqrt5);
        Tar(4,2)=  1.0_wp/178.0_wp *(95.0_wp - 17.0_wp*sqrt5);
        Tar(4,3)=  3.0_wp/89.0_wp *(19.0_wp + sqrt5);
        Tar(4,4:6) = zero
      else
        Tar(4,1)= 1.0_wp/178.0_wp *(-31.0_wp+11.0_wp*sqrt5);
        Tar(4,2)= 1.0_wp/178.0_wp *(95.0_wp - 17.0_wp*sqrt5);
        Tar(4,3)= x66
        Tar(4,4)= 3.0_wp/89.0_wp *(19.0_wp + sqrt5) - x66;
        Tar(4,5)= dxR*(3.0_wp*(-5+sqrt(5.0_wp))+(25.0_wp-6.0_wp*sqrt(5.0_wp))*x66)/   &
                  (6.0_wp*(-5.0_wp+sqrt(5.0_wp))*dxR - 5.0_wp*dx);
        Tar(4,6)= -Tar(4,5);
      endif

      ! calculate coefficients for the 2nd smoothness indicator tau

      tau_cL4(1) = 5.0_wp*(-1.0_wp+sqrt(5.0_wp))*dx**3/(2.0_wp*(dx+dxL)*    &
      ((-5.0_wp+3.0_wp*sqrt(5.0_wp))*dx**2+(-19.0_wp+9.0_wp*sqrt(5.0_wp))*dx*dxL+(-11.0_wp &
      +5.0_wp*sqrt(5.0_wp))*dxL**2))
      tau_cL4(2) = -1.0_wp*(5.0_wp+sqrt(5.0_wp))/4.0_wp
      tau_cL4(3) = 5.0_wp*(1.0_wp+sqrt(5.0_wp))*dxL/(4.0_wp*(dx+dxL))
      tau_cL4(4) = 5.0_wp*(-1.0_wp+sqrt(5.0_wp))*dxL/(-4.0_wp*dx+2.0_wp*(-3.0_wp+sqrt(5.0_wp))*dxL)
      tau_cL4(5) = 5.0_wp*dxL/(10.0_wp*dx-(-5.0_wp+sqrt(5.0_wp))*dxL)
      tau_cfL    = 0.5_wp*(tau_coeffT0+tau_cL4)

      tau_cR4(1) = 5.0_wp*dxR/(10.0_wp*dx-(-5.0_wp+sqrt(5.0_wp))*dxR)
      tau_cR4(2) = 5.0_wp*(-1.0_wp+sqrt(5.0_wp))*dxR/(-4.0_wp*dx+2.0_wp*(-3.0_wp+sqrt(5.0_wp))*dxR)
      tau_cR4(3) = 5.0_wp*(1.0_wp+sqrt(5.0_wp))*dxR/(4.0_wp*(dx+dxR))
      tau_cR4(4) = -1.0_wp*(5.0_wp+sqrt(5.0_wp))/4.0_wp
      tau_cR4(5) = 5.0_wp*(-1.0_wp+sqrt(5.0_wp))*dx**3/(2.0_wp*(dx+dxR)*    &
      ((-5.0_wp+3.0_wp*sqrt(5.0_wp))*dx**2+(-19.0_wp+9.0_wp*sqrt(5.0_wp))*dx*dxR+(-11.0_wp &
      +5.0_wp*sqrt(5.0_wp))*dxR**2))
      tau_cfR    = 0.5_wp*(tau_coeffT1+tau_cR4)

      ! calculate scaling parameter

      eps = Ceps*dx0**3

      ! Interior flux points in element
      do i = 1,ixd-1

        uin(:,:) = 0.0_wp; vin(:,:) = 0.0_wp; nin(:,:) = 0.0_wp
        if(i.eq.1.and.inb.ne.-1)then
          uin(1:nq,1)   = uLL(1:nq)
          uin(1:nq,2)   =  uL(1:nq)
          uin(1:nq,3:6) = uint(1:nq,1:4)
          call conserved_to_primitive(uin(:,3),vL(:),nq)   !  primitives
          call conserved_to_primitive(uin(:,4),vR(:),nq)   !  primitives
          call roeavg(vL, vR, Vav, nq)   

          nin(:,1)   = nLL(:)
          nin(:,2)   =  nL(:)
          nin(:,3:6) = nxint(:,1:4)
          nxL(:)     = nxint(:,3)
          nxR(:)     = nxint(:,4)
          nxave(:)   = 0.5_wp*(nxL(:) + nxR(:))
          iL = 1 ;  iH = 6 ;
        elseif(i.eq.3.and.inb.ne.+1)then
          uin(1:nq,1:4) = uint(1:nq,1:4)
          uin(1:nq,5)   = uR(1:nq)
          uin(1:nq,6)   = uRR(1:nq)
          call conserved_to_primitive(uin(:,3),vL(:),nq)   !  primitives
          call conserved_to_primitive(uin(:,4),vR(:),nq)   !  primitives
          call roeavg(vL, vR, Vav, nq)   

          nin(:,1:4) = nxint(:,1:4)
          nin(:,5)   = nR(:)
          nin(:,6)   = nRR(:)
          nxL(:)     = nxint(:,3)
          nxR(:)     = nxint(:,4)
          nxave(:)   = 0.5_wp*(nxL(:) + nxR(:))
          iL = 1 ;  iH = 6 ;
        else
          uin(1:nq,1:4) = uint(1:nq,1:4)
          uin(1:nq,5:6) = zero
          call conserved_to_primitive(uin(:,i  ),vL(:),nq)   !  primitives
          call conserved_to_primitive(uin(:,i+1),vR(:),nq)   !  primitives
          call roeavg(vL, vR, Vav, nq)   

          nin(:,1:4) = nxint(:,1:4)
          nin(:,5:6) = zero
          nxL(:)     = nxint(:,i  )
          nxR(:)     = nxint(:,i+1)
          nxave(:)   = 0.5_wp*(nxL(:) + nxR(:))
          iL = 1 ;  iH = 4 ;
        endif
        
        uhat(:) = 0.0_wp ; metr(:) = 0.0_wp ; cav2(:) = 0.0_wp ;
        do j = iL,iH
          call conserved_to_primitive(uin(:,j),vin(:,j),nq)    ! primitives

          uhat(j)  =  abs(dot_product(vin(2:4,j),nin(:,j)))! normal velocity  |u.n|
          metr(j)  = sqrt(dot_product(nin( : ,j),nin(:,j)))! computational coordinate scaling
          cav2(j)  = sqrt(gamma0*vin(5,j)/gM2) * metr(j)   ! Speed of sound * metric scaling
        enddo

        lambda     = maxval(uhat(iL:iH) + cav2(iL:iH))             !  Max eigenvalue in pencil

        do j = 1,6
          fcp(:,j) = normalflux(vin(:,j),nin(:,j),nq) + lambda * uin(:,j)
          fcm(:,j) = normalflux(vin(:,j),nin(:,j),nq) - lambda * uin(:,j)
        enddo

        call CharacteristicDecomp(Vav,nq, Sinv, Smat, ev, nxave )

        qin(:,:) = 0.0_wp
        do j = 1,6
          qin(:,j) = matmul(Sinv,uin(:,j))
        enddo

        ! Smoothness indicators
        if(i.eq.1.and.inb.ne.-1)then
         beta(:,1) = (qin(:,2) - qin(:,1)  )**2
         beta(:,2) = (qin(:,3) - qin(:,1)  )**2
         beta(:,3) = (qin(:,4) - qin(:,2)  )**2
         beta(:,4) = (qin(:,4) - qin(:,3)  )**2
         beta(:,5) = (qin(:,5) - qin(:,4)  )**2
         beta(:,6) = (qin(:,6) - qin(:,5)  )**2
         tau(:)    = ( qin(:,1)*tau_cfL(1)+0.5_wp*(qin(:,2)+qin(:,3))*tau_cfL(2)         &
                   +   qin(:,4)*tau_cfL(3)+ qin(:,5)*tau_cfL(4)+ qin(:,6)*tau_cfL(5))**2
         n         = 6
        elseif(i.eq.3.and.inb.ne.+1)then
         beta(:,1) = (qin(:,2) - qin(:,1)  )**2
         beta(:,2) = (qin(:,3) - qin(:,2)  )**2
         beta(:,3) = (qin(:,4) - qin(:,3)  )**2
         beta(:,4) = (qin(:,5) - qin(:,3)  )**2
         beta(:,5) = (qin(:,6) - qin(:,4)  )**2
         beta(:,6) = (qin(:,6) - qin(:,5)  )**2
         tau(:)    = (qin(:,1)*tau_cfR(1)+qin(:,2)*tau_cfR(2)+qin(:,3)*tau_cfR(3)  &
                   + 0.5_wp*(qin(:,4)+qin(:,5))*tau_cfR(4) + qin(:,6)*tau_cfR(5))**2
         n         = 6
        else
         beta(:,1) = (qin(:,2) - qin(:,1)  )**2
         beta(:,2) = (qin(:,3) - qin(:,2)  )**2
         beta(:,3) = (qin(:,4) - qin(:,3)  )**2
         beta(:,4) = 0.0_wp

         do k=1,nq
           tau(k)    = (dot_product(tau_coeff(1:4),qin(k,1:4)))**2 ! Second smoothness indicator
         enddo
         n         = 3
        endif

!    positive and negative propagating fluxes

        ! split target weights
        dh  = Tar(i+1,:)
        dhp = 0.5_wp*(dh+theta*abs(dh))
        dhm = dhp-dh
        Lsigp = sum(dhp) ; dhp = dhp/Lsigp
        Lsigm = sum(dhm) ; dhm = dhm/Lsigm

!    Positive propagating speed

        ! candidate fluxes for fbarC(i) -- positive propagating speed
        if(i.eq.1.and.inb.ne.-1)then
!         fh(1:nq,1) = IP1(5-i,1)*fcp(1:nq,1) + IP1(5-i,2)*fcp(1:nq,2)
!         fh(1:nq,2) = IP2(5-i,1)*fcp(1:nq,1) + IP2(5-i,2)*fcp(1:nq,3)
          fh(1:nq,1) = IP1n(5-i,1)*fcp(1:nq,1) + IP1n(5-i,2)*fcp(1:nq,2)
          fh(1:nq,2) = IP2n(5-i,1)*fcp(1:nq,1) + IP2n(5-i,2)*fcp(1:nq,3)
          fh(1:nq,3) = IM2(i+1,1)*fcp(1:nq,2) + IM2(i+1,2)*fcp(1:nq,4)
          fh(1:nq,4) = IM1(i+1,1)*fcp(1:nq,3) + IM1(i+1,2)*fcp(1:nq,4)
          fh(1:nq,5) = IC (i+1,1)*fcp(1:nq,4) + IC (i+1,2)*fcp(1:nq,5)
          fh(1:nq,6) = IP1(i+1,1)*fcp(1:nq,5) + IP1(i+1,2)*fcp(1:nq,6)
        elseif(i.eq.3.and.inb.ne.+1)then
          fh(1:nq,1) = IM1(i+1,1)*fcp(1:nq,1) + IM1(i+1,2)*fcp(1:nq,2)
          fh(1:nq,2) = IC (i+1,1)*fcp(1:nq,2) + IC (i+1,2)*fcp(1:nq,3)
          fh(1:nq,3) = IP1(i+1,1)*fcp(1:nq,3) + IP1(i+1,2)*fcp(1:nq,4)
          fh(1:nq,4) = IP2(i+1,1)*fcp(1:nq,3) + IP1(i+1,2)*fcp(1:nq,5)
          fh(1:nq,5) = IM2n(5-i,1)*fcp(1:nq,4) + IM2n(5-i,2)*fcp(1:nq,6)
          fh(1:nq,6) = IM1n(5-i,1)*fcp(1:nq,5) + IM1n(5-i,2)*fcp(1:nq,6)
!         fh(1:nq,5) = IM2(5-i,1)*fcp(1:nq,4) + IM2(5-i,2)*fcp(1:nq,6)
!         fh(1:nq,6) = IM1(5-i,1)*fcp(1:nq,5) + IM1(5-i,2)*fcp(1:nq,6)
        else
          fh(1:nq,1) = IM1(i+1,1)*fcp(1:nq,1) + IM1(i+1,2)*fcp(1:nq,2)
          fh(1:nq,2) = IC (i+1,1)*fcp(1:nq,2) + IC (i+1,2)*fcp(1:nq,3)
          fh(1:nq,3) = IP1(i+1,1)*fcp(1:nq,3) + IP1(i+1,2)*fcp(1:nq,4)
          fh(1:nq,4) = 0.0_wp
          fh(1:nq,5:6) = 0.0_wp
        endif

        do j = 1,n-1
          tbratp1(1:nq,j) = (1.0_wp+tau(1:nq)/(beta(1:nq,j)+eps))       
        enddo

        betam(1:nq) = 0.0_wp
        do k = 1,n
          betam(1:nq) = betam(1:nq) + beta(1:nq,k)*beta(1:nq,k)*beta(1:nq,k)*beta(1:nq,k)
        enddo
        betam = sqrt(sqrt(betam/n))

        betad(1:nq)     = WENO_Bias*beta(1:nq,n) + (1.0_wp-WENO_Bias)*betam(1:nq)

        tbratp1(1:nq,n) = (1.0_wp+tau(1:nq)/(betad(1:nq)+eps))

        ! calculate split weights
        wpsum(1:nq) = 0.0_wp ; wmsum(1:nq) = 0.0_wp
        wpos(:,:)   = 0.0_wp ;   wneg(:,:) = 0.0_wp
        do j = 1,n                                       ! positive & negative weights
          wneg(1:nq,j) = dhm(j)*tbratp1(1:nq,j)   ; 
          wpos(1:nq,j) = dhp(j)*tbratp1(1:nq,j)   ; 

          wmsum(1:nq)  = wmsum(1:nq) + wneg(1:nq,j) ;
          wpsum(1:nq)  = wpsum(1:nq) + wpos(1:nq,j) ; 
        end do
        do k = 1,nq
          wneg(k,:) = wneg(k,:)/wmsum(k) ;
          wpos(k,:) = wpos(k,:)/wpsum(k) ; 
        enddo
        weights(:,:) = 0.0_wp ; 
        weights(:,:) = wpos(:,:)*Lsigp-wneg(:,:)*Lsigm
        
        fhC(:,:) = 0.0_wp
        do j = 1,6
          fhC(:,j) = matmul(Sinv,fh(:,j))
        enddo

        fbarC(:,:) = 0.0_wp
        do k = 1,nq
          fbarC(k,i) = 0.5_wp*dot_product(weights(k,:),fhC(k,:))
        enddo

!       split target weights
        dh  = Tar(i+1,:)
        dhp = 0.5_wp*(dh+theta*abs(dh))
        dhm = dhp-dh
        Lsigp = sum(dhp) ; dhp = dhp/Lsigp
        Lsigm = sum(dhm) ; dhm = dhm/Lsigm
        
!       Negative propagating speed

        ! candidate fluxes for fbarC(i) -- negative propagating speed
        if(i.eq.1.and.inb.ne.-1)then
          fh(1:nq,1) = IP1(5-i,1)*fcm(1:nq,1) + IP1(5-i,2)*fcm(1:nq,2)
          fh(1:nq,2) = IP2(5-i,1)*fcm(1:nq,1) + IP2(5-i,2)*fcm(1:nq,3)
          fh(1:nq,3) = IM2(i+1,1)*fcm(1:nq,2) + IM2(i+1,2)*fcm(1:nq,4)
          fh(1:nq,4) = IM1(i+1,1)*fcm(1:nq,3) + IM1(i+1,2)*fcm(1:nq,4)
          fh(1:nq,5) = IC (i+1,1)*fcm(1:nq,4) + IC (i+1,2)*fcm(1:nq,5)
          fh(1:nq,6) = IP1(i+1,1)*fcm(1:nq,5) + IP1(i+1,2)*fcm(1:nq,6)
        elseif(i.eq.3.and.inb.ne.+1)then
          fh(1:nq,1) = IM1(i+1,1)*fcm(1:nq,1) + IM1(i+1,2)*fcm(1:nq,2)
          fh(1:nq,2) = IC (i+1,1)*fcm(1:nq,2) + IC (i+1,2)*fcm(1:nq,3)
          fh(1:nq,3) = IP1(i+1,1)*fcm(1:nq,3) + IP1(i+1,2)*fcm(1:nq,4)
          fh(1:nq,4) = IP2(i+1,1)*fcm(1:nq,3) + IP1(i+1,2)*fcm(1:nq,5)
          fh(1:nq,5) = IM2(5-i,1)*fcm(1:nq,4) + IM2(5-i,2)*fcm(1:nq,6)
          fh(1:nq,6) = IM1(5-i,1)*fcm(1:nq,5) + IM1(5-i,2)*fcm(1:nq,6)
        else
          fh(1:nq,1) = IM1(i+1,1)*fcm(1:nq,1) + IM1(i+1,2)*fcm(1:nq,2)
          fh(1:nq,2) = IC (i+1,1)*fcm(1:nq,2) + IC (i+1,2)*fcm(1:nq,3)
          fh(1:nq,3) = IP1(i+1,1)*fcm(1:nq,3) + IP1(i+1,2)*fcm(1:nq,4)
          fh(1:nq,4) = 0.0_wp
          fh(1:nq,5:6) = 0.0_wp
        endif

        ! biasing the stencil in the upwind direction

        betad(1:nq)     = WENO_Bias*beta(1:nq,1) + (1.0_wp-WENO_Bias)*betam(1:nq)

        tbratp1(1:nq,1) = (1.0_wp + tau(1:nq)/(betad(1:nq)+eps))
        do j = 2,n
          tbratp1(1:nq,j) = (1.0_wp + tau(1:nq)/(beta(1:nq,j)+eps))
        enddo

        ! calculate split weights
        wpsum(1:nq) = 0.0_wp ; wmsum(1:nq) = 0.0_wp
        wpos(:,:)   = 0.0_wp ;   wneg(:,:) = 0.0_wp
        do j = 1,n
          wneg(1:nq,j) = dhm(j)*tbratp1(1:nq,j)   ; 
          wpos(1:nq,j) = dhp(j)*tbratp1(1:nq,j)   ; 

          wmsum(1:nq)  = wmsum(1:nq)+wneg(1:nq,j) ;
          wpsum(1:nq)  = wpsum(1:nq)+wpos(1:nq,j) ; 
        end do
        do k = 1,nq
          wneg(k,:) = wneg(k,:)/wmsum(k) ;
          wpos(k,:) = wpos(k,:)/wpsum(k) ; 
        enddo
        weights(:,:) = 0.0_wp
        weights(:,:) = wpos(:,:)*Lsigp-wneg(:,:)*Lsigm

        fhC(:,:) = 0.0_wp
        do j = 1,6
          fhC(:,j) = matmul(Sinv,fh(:,j))
        enddo

        do k = 1,nq
          fbarC(k,i) = fbarC(k,i) + 0.5_wp*dot_product(weights(k,:),fhC(k,:))
        enddo

        fbarW(:,i) = matmul(Smat,fbarC(:,i))

      end do

      call conserved_to_primitive(uint(:,1),vL(:),nq)   !  primitives
      call conserved_to_primitive(uint(:,4),vR(:),nq)   !  primitives
      fbarW(:,  0) = normalflux(vL(:),nxint(:,1),nq)
      fbarW(:,ixd) = normalflux(vR(:),nxint(:,4),nq)

      do i=1,ixd
        call conserved_to_primitive    (uint(:,i),vint(:,i),nq)
        call primitive_to_entropy      (vint(:,i),wint(:,i),nq)
!       call KineticEnergyVariables(vint(:,i),Kint(:,i),nq)
      enddo

      !  Form Entropy Fluxes 0:N_F
      fnS(:,:) = 0.0_wp ; !fnK(:,:) = 0.0_wp ;
      do i=1,ixd-1
  
        select case(entropy_flux)
    
          case('Ismail_Roe')

            do j=i+1,ixd
                   nxave(:) = 0.5_wp * (nxint(:,i) + nxint(:,j))
              SSFlux(:,i,j) = two  * qmat(i,j)*EntropyConsistentFlux     (vint(:,i),vint(:,j),nxave,nq)
            enddo
!         DSFlux(:,i,j) = two  * qmat(i,j)*HoneinMoinFlux(vint(:,i),vint(:,j),nxave,nq)
  
          case('Chandrashekar')

            do j=i+1,ixd
                   nxave(:) = 0.5_wp * (nxint(:,i) + nxint(:,j))
              SSFlux(:,i,j) = two  * qmat(i,j)*Entropy_KE_Consistent_Flux(vint(:,i),vint(:,j),nxave,nq)
            enddo

        end select

        fnS(:,i) = 0.0_wp ; !fnK(:,i) = 0.0_wp ;
        do k=i+1,ixd
          do l=1,i
            fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
!           fnK(:,i) = fnK(:,i) + DSFlux(:,l,k)
          enddo
        enddo
  
      enddo

      fnS(:,  0) = normalflux(vint(:,  1),nxint(:,  1),nq) 
      fnS(:,ixd) = normalflux(vint(:,ixd),nxint(:,ixd),nq)
!     fnK(:,  0) = fnS(:,  0) ;
!     fnK(:,ixd) = fnS(:,ixd) ;

      !  Entropy Correction 
      do i = 1,ixd-1
  
!     Component-wise correction
              bb(:) = (wint(:,i+1) - wint(:,i+0))*(fnS(:,i) - fbarW(:,i))
              ds(:) = sqrt(bb(:)*bb(:) + cc2)
           delta(:) = (ds(:) - bb(:))/(two*ds(:))
         fbarW(:,i) = fbarW(:,i) + delta(:)*(fnS(:,i) - fbarW(:,i))  !  Entropy
  
!     Vector contraction correction
!               bbS = dot_product(wint(:,i+1)-wint(:,i+0),fnS(:,i)-fbarW(:,i))
!               dsS = sqrt(bbS*bbS + cc2)
!            deltaS = (dsS - bbS)/(two*dsS)
!        fbarW(:,i) = fbarW(:,i) + deltaS * (fnS(:,i) - fbarW(:,i))
  
!     Vector Entropy / KE correction
!         Wmat(1,:) = wint(:,i+1) - wint(:,i+0) 
!         Wmat(2,:) = Kint(:,i+1) - Kint(:,i+0) 
! 
!             Hmat  = matmul(Wmat,Transpose(Wmat))
!                a  = Hmat(1,1)+cc1 ; b = Hmat(1,2) ; c = Hmat(2,2)+cc1 ; r  = sqrt( (0.5_wp*(a-c))**2 + b*b + cc1)
! 
!          Lam(:,:) = 0.0_wp
!          Lam(1,1) = 1.0_wp / (0.5_wp*(a+c) + r + cc1)
!          Lam(2,2) = 1.0_wp / (0.5_wp*(a+c) - r + cc1)
! 
!               x1  = sqrt(2 + (c - a) / r) ; x2  = sqrt(2 - (c - a) / r) ;
!         Tmat(1,1) = + b / (r*x1)  ; Tmat(1,2) = 0.5_wp*x1  ;
!         Tmat(2,1) = - b / (r*x2)  ; Tmat(2,2) = 0.5_wp*x2  ;
!             Hinv  = matmul(Transpose(Tmat),matmul(Lam,Tmat))
  
!          bSK(1)   = dot_product(wint(:,i+1)-wint(:,i+0),fnS(:,i)-fbarW(:,i))
!          bSK(2)   = dot_product(Kint(:,i+1)-Kint(:,i+0),fnK(:,i)-fbarW(:,i))
!          bSK(:)   = 0.5_wp*(bSK(:) - sqrt(bSK(:)*bSK(:) + cc2))
  
!        fbarW(:,i) = fbarW(:,i) + matmul(Transpose(Wmat),matmul(Hinv,bSK))
  
      enddo
  
      do i = 1,ixd
          dfbar(:,i) = pinv(i) * (fbarW(:,i) - fbarW(:,i-1))
      enddo
  
    end subroutine SSWENO4

!============================================================================

    subroutine Planar_Interpolation_WENO()

      use precision_vars
      use referencevariables
      use variables
      use interpolation
      use initcollocation,      only: element_properties

      implicit none

      integer    :: ielem, iface, ipen, jnode
      integer    :: elem_face_nodes

      elem_face_nodes = nodesperface*nfacesperelem

      ! loop over all elements
      do ielem = ihelems(1), ihelems(2)

        call element_properties(ielem, n_pts_3d=nodesperelem)

        do iface = 1,nfacesperelem

          if (ef2e(1,iface,ielem) < 0) then
          else

            do ipen = 1,nodesperface

              ! Index in facial ordering
              jnode = nodesperface*(iface-1) + ipen
  
              call Extrp_XiEtaZeta_neq(ndim,nequations,nodesperedge,XIWENO_partner(:,jnode,ielem), &
                                  ug(:,:,ielem), ugWENO_partner(:,jnode,ielem))
      
            end do

          end if

        end do

      end do

    end subroutine Planar_Interpolation_WENO

!============================================================================

    subroutine Calc_Entropy_Viscosity

      use precision_vars
      use nsereferencevariables
      use referencevariables
      use variables
      use interpolation
      use initcollocation,      only: element_properties

      real(wp), parameter :: c1 = 0.5_wp   !  Guermond's magic constant number one
      real(wp), parameter :: c2 = 1.0_wp   !  Guermond's magic constant number one

      real(wp) :: ev, evmax

      real(wp) :: t1,t2
      integer  :: ielem, inode

      ! loop over all elements
      do ielem = ihelems(1), ihelems(2)

        call element_properties(ielem, n_pts_3d=nodesperelem)

        ! calcualte the max - max eigenvalue on element (times density)
        evmax = 0.0_wp
        do inode = 1,nodesperelem
     
          ev = vg(1,inode,ielem)*(sqrt(abs(gamma0*vg(5,inode,ielem)/gM2))   &
                                + sqrt(dot_product(vg(2:4,inode,ielem),vg(2:4,inode,ielem))))
          if( ev >= evmax ) evmax = ev
  
        enddo

        t1 = c1 * evmax * dx_min_elem(ielem) 
        t2 = c2 * dx_min_elem(ielem)*dx_min_elem(ielem) * maxval(abs(dudt_S(:,ielem)))
  
        ! currently the mut is assumed constant over each element
        mut(:,ielem) = min(t1,t2)

      enddo

    end subroutine Calc_Entropy_Viscosity
      
  !============================================================================

    function SAT_Inv_Vis_Flux(neq,iface,ielem,vg_On,vg_Off,phig_On,phig_Off,nx_On,nx_Off,Jx_r,pinv,mut)

      use precision_vars
  
      use controlvariables,     only: Riemann_Diss, entropy_flux_BC
      use collocationvariables, only: l01, l00, Sfix, ldg_flip_flop_sign, alpha_ldg_flip_flop


      integer,                       intent(in) :: neq, iface, ielem
      real(wp),  dimension(neq),     intent(in) ::   vg_On,   vg_Off
      real(wp),  dimension(neq,3),   intent(in) :: phig_On, phig_Off
      real(wp),  dimension(3),       intent(in) :: nx_On, nx_Off
      real(wp),                      intent(in) :: Jx_r, pinv, mut

      real(wp), parameter :: Cevmax          =  1.0_wp
      real(wp), parameter :: deltaU          =  0.1_wp
      real(wp), parameter :: LocalLaxF_factor=  2.0_wp

      real(wp),  dimension(neq,neq)             :: smat, sinv, matrix_ip

      real(wp),  dimension(neq)                 :: fLLF
      real(wp),  dimension(neq)                 :: wg_On, wg_Off
      real(wp),  dimension(neq)                 :: ug_On, ug_Off
      real(wp),  dimension(neq)                 :: vav, ev, evabs
      real(wp),  dimension(neq)                 :: fn, fstar, fV_Del
      real(wp)                                  :: evmax
      real(wp)                                  :: l01_ldg_flip_flop

      real(wp),  parameter                      :: tol = 2.0e-12_wp
      real(wp)                                  :: UavAbs, switch

      real(wp),  dimension(neq)                 :: SAT_Inv_Vis_Flux

         call primitive_to_entropy (vg_On (:),wg_On (:),neq)
         call primitive_to_entropy (vg_Off(:),wg_Off(:),neq)

         call roeavg( vg_On (:), vg_Off(:), Vav, neq )   
!-- DAVID DEBUG START
         !-- I have arbitrarily set nx to nx_On do not know if that is the correct choice
         call CharacteristicDecomp( vav, neq, sinv, smat, ev, nx_On )      
!-- DAVID DEBUG END
         evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)

         select case(Riemann_Diss)

           case('Roe')

             select case(entropy_flux_BC)
               case('Ismail_Roe')
                 fstar = EntropyConsistentFlux     (vg_On(:), vg_Off(:), nx_Off, neq )
               case('Chandrashekar') 
                 fstar = Entropy_KE_Consistent_Flux(vg_On(:), vg_Off(:), nx_Off, neq )
             end select

             fstar = fstar + 0.5_wp * matmul(smat,evabs*matmul(transpose(smat), wg_On(:)-wg_Off(:)) )

           case('LocalLaxF')
             call primitive_to_conserved(vg_On (:),ug_On (:),neq)
             call primitive_to_conserved(vg_Off(:),ug_Off(:),neq)
             fLLF  = half * ( normalflux( vg_On (:), nx_On, neq )    &
                          +   normalflux( vg_Off(:), nx_Off, neq )    &
                          +   LocalLaxF_factor*evmax*(ug_On(:)-ug_Off(:)) )

             fstar = fLLF

           case('RoeLF')

             call primitive_to_conserved(vg_On (:),ug_On (:),neq)
             call primitive_to_conserved(vg_Off(:),ug_Off(:),neq)

             UavAbs = 0.5_wp*(ug_On(1)+ug_Off(1)) ; switch = 0.0_wp
             if(UavAbs <= 1.0e-7_wp)then
               switch = 1.0_wp
             else
               switch = 1.0_wp/deltaU * abs((ug_On(1)-ug_Off(1))/UavAbs)
             endif
             if(switch.gt.1.0_wp)switch = 1.0_wp

             select case(entropy_flux_BC)
               case('Ismail_Roe')
                 fstar = EntropyConsistentFlux     (vg_On(:), vg_Off(:), nx_Off, neq )
               case('Chandrashekar') 
                 fstar = Entropy_KE_Consistent_Flux(vg_On(:), vg_Off(:), nx_Off, neq )
             end select

             fstar = fstar + 0.5_wp * matmul(smat,evabs*matmul(transpose(smat), wg_On(:)-wg_Off(:)) )

             fLLF  = half * ( normalflux( vg_On (:), nx_Off, neq )    &
                          +   normalflux( vg_Off(:), nx_Off, neq )    &
                          +   LocalLaxF_factor*evmax*(ug_On(:)- ug_Off(:))  )

             fstar = switch*fLLF + (1.0_wp-switch) * fstar

         end select

                                                                                           
         l01_ldg_flip_flop = l01*(1.0_wp - ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop) ! Flip-flop sign convention

         fV_Del = normalviscousflux(vg_On (:), phig_On (:,:), nx_On , neq, mut)  &              ! Add LDG viscous flux dissipation
              & - normalviscousflux(vg_Off(:), phig_Off(:,:), nx_Off, neq, mut)

         matrix_ip = 0.5_wp * pinv * ( matrix_hatc_node(vg_On (:),nx_On,nx_Off,neq) &           ! Average of c_ii_L and cii_R matrix  appropriately 
                                     + matrix_hatc_node(vg_Off(:),nx_On,nx_Off,neq)) / Jx_r     ! scaled by inverse norm and size of elements

         fn = normalflux( vg_On (:), nx_On, neq )                                               ! (Euler Flux)

         SAT_Inv_Vis_Flux = + (fn - fstar)                 &                                    ! Inviscid penalty
                            + l01_ldg_flip_flop * fV_Del   &                                    ! Visous penalty
                            - l00*matmul(matrix_ip,wg_On (:)-wg_Off(:))                         ! Interior penalty approach

     end function

  !============================================================================

    function SAT_Vis_Diss(neq,vg_On,vg_Off,nx)

      use precision_vars

      use collocationvariables, only: Sfix

      integer,                       intent(in) :: neq
      real(wp),  dimension(neq),     intent(in) ::   vg_On,   vg_Off
      real(wp),  dimension(3),       intent(in) :: nx

      real(wp),  dimension(neq,neq)             :: smat,sinv


      real(wp),  dimension(neq)                 :: vav, ev, evabs
      real(wp)                                  :: evmax

!     real(wp),  dimension(neq)                 ::   ug_On, ug_Off
      real(wp),  dimension(neq)                 ::   wg_On, wg_Off
 

      real(wp),  dimension(neq)                 :: SAT_Vis_Diss

!        call primitive_to_conserved(vg_On (:),ug_On (:),neq)
!        call primitive_to_conserved(vg_Off(:),ug_Off(:),neq)
         call primitive_to_entropy(vg_On (:),wg_On (:),neq)
         call primitive_to_entropy(vg_Off(:),wg_Off(:),neq)

         call roeavg( vg_On (:), vg_Off(:), Vav, neq )   

         call CharacteristicDecomp( vav, neq, sinv, smat, ev, nx )      

         evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)

!        SAT_Vis_Diss = - 0.5_wp * matmul(smat,evabs*matmul(          sinv , ug_On (:)-ug_Off(:)) )
         SAT_Vis_Diss = - 0.5_wp * matmul(smat,evabs*matmul(transpose(smat), wg_On (:)-wg_Off(:)) )

          return
     end function

  !============================================================================

      subroutine viscous_gradients()

        use referencevariables
        use nsereferencevariables
        use variables,            only: ef2e, efn2efn, efn2efn_Gau,       &
                                      & phig, phig_err, grad_w_jacobian,  &
                                      & vg, wg, ughst,                    &
                                      & r_x, facenodenormal, nelem_ghst
        use collocationvariables, only: nnzgrad,iagrad,jagrad,dagrad,pinv,l10, &
                                      & ldg_flip_flop_sign, alpha_ldg_flip_flop,&
                                      & elem_props, &
                                      & Restrct_Gau_2_LGL_1d, Prolong_LGL_2_Gau_1d
        use initcollocation,      only: ExtrpXA2XB_2D_neq, &
                                        & Gauss_Legendre_points, element_properties

        implicit none

        ! temporary arrays for phi and delta phi
        real(wp), dimension(:),   allocatable :: dphi, ug_Off, vg_On, vg_Off, wg_On, wg_Off
        real(wp), dimension(:,:), allocatable :: phig_tmp1, phig_tmp2, phig_err1
        real(wp), dimension(:),   allocatable :: x_S_1d_Mort, w_S_1d_Mort
        real(wp), dimension(:),   allocatable :: x_S_1d_On , x_S_1d_Off
        real(wp), dimension(:,:), allocatable :: wg_2d_Mort_On, wg_2d_Mort_Off, wg_2d_On, wg_2d_Off

        real(wp), allocatable, dimension(:,:) :: Extrp_Off, Extrp_On
        real(wp), allocatable, dimension(:,:) ::            Intrp_On

        integer,  allocatable, dimension(:)   :: ifacenodes_On, ifacenodes_Off

        ! normal direction at face
        real(wp) :: nx(3)

        ! loop indices
        integer :: jdir, idir
        integer :: i,j,k,l
        integer :: n_S_1d_Mort, n_S_2d_Mort
        integer :: n_S_1d_On, n_S_1d_Off
        integer :: n_S_2d_On, n_S_2d_Off
        integer :: n_S_1d_max, n_S_2d_max
        integer :: n_S_3d_On
        integer :: inode, jnode, knode, gnode, lnode
        integer :: ielem, iface
        integer :: kelem, kface
        integer :: poly_val
        integer :: nghst_volume, nghst_Mortar

        real(wp) :: l10_ldg_flip_flop

        allocate(phig_tmp1(nequations,ndim),phig_tmp2(nequations,ndim),phig_err1(nequations,ndim))           ! phig_tmp1 is calculated in computational space
        allocate(dphi(nequations),ug_Off(nequations)) 
        allocate(vg_On(nequations),vg_Off(nequations),wg_On(nequations),wg_Off(nequations))                  ! tmp variables: vg_On, vg_Off, wg_On and wg_Off

        ! loop over all elements
        element_Loop:do ielem = ihelems(1), ihelems(2)

          call element_properties(ielem,              &
                                n_pts_1d=n_S_1d_On,   &
                                x_pts_1d=x_S_1d_On,   &
                                n_pts_2d=n_S_2d_On,   &
                                n_pts_3d=n_S_3d_On,   &
                                    pinv=pinv,        &
                                 nnzgrad=nnzgrad,     &
                                  iagrad=iagrad,      &
                                  jagrad=jagrad,      &
                                  dagrad=dagrad,      &
                              ifacenodes=ifacenodes_On)

                                                                      ! compute computational gradients of the entropy variables
                 phig    (:,:,:,ielem) = 0.0_wp                       ! initialize phi
                 phig_err(:,:,:,ielem) = 0.0_wp                       ! initialize phi_err used for error estimate
          grad_w_jacobian(:,:,:,ielem) = 0.0_wp

          do inode = 1, n_S_3d_On                                     ! loop over every node in element

            phig_tmp1(:,:) = 0.0_wp ; phig_tmp2(:,:) = 0.0_wp ;

            do i = iagrad(inode), iagrad(inode+1)-1                   ! loop over number of dependent elements in gradient

              do jdir = 1,ndim                                        ! loop over dimensions
                                                                      ! column/node from gradient operator in CSR format in
                jnode = jagrad(jdir,i)                                ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
                phig_tmp1(:,jdir) = phig_tmp1(:,jdir) + dagrad(jdir,i)*wg(:,jnode,ielem) ! update gradient using coefficient and entropy variables at appropriate node
                phig_tmp2(:,jdir) = phig_tmp2(:,jdir) + dagrad(jdir,i)*vg(:,jnode,ielem) 
              end do
            end do

            grad_w_jacobian(:,:,inode,ielem) = phig_tmp1(:,:)         ! Store gradient of the entropy variables in computational space

            phig_err1(:,:) = 0.0_wp                                   ! transform to physical space using dxi_jdir/dx_idir
            do jdir = 1,ndim
              do idir = 1,ndim
                phig(:,idir,inode,ielem) = phig(:,idir,inode,ielem) + phig_tmp1(:,jdir)*r_x(jdir,idir,inode,ielem)
                phig_err1(:,idir)        = phig_err1(:,idir)        + phig_tmp2(:,jdir)*r_x(jdir,idir,inode,ielem)
              end do
            end do

            phig_err(:,:,inode,ielem) = matmul(dWdV(vg(:,inode,ielem),nequations),phig_err1(:,:))  &
                                      - phig(:,:,inode,ielem)
  
          end do

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! LDC/LDG penalty on phig
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            
          nghst_volume = nelem_ghst(1,ielem)                  ! point at beginning of ``ielem'' data in ghost stack 
          nghst_Mortar = nelem_ghst(3,ielem)                  ! point at beginning of ``ielem'' data in ghost stack 

          face_loop:do iface = 1, nfacesperelem               ! loop over faces
  
            if (ef2e(1,iface,ielem) < 0) then                 ! face if cycle for BC or different face-types

              cycle
                                               !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                                               !       CONFORMING INTERFACES:  polynomial orders match and h conforming
                                               !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            else if ((elem_props(2,ielem) == ef2e(4,iface,ielem)) .and. (ef2e(9,iface,ielem) == 0)) then
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
              if (ef2e(3,iface,ielem) /= myprocid) then         !       Off-Processor Contributions to gsat:  Conforming Interface
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                do i = 1,n_S_2d_On

                  jnode = n_S_2d_On*(iface-1) + i

                  inode = ifacenodes_On(jnode)                       ! corresponding volumetric node for face node

                  wg_On(:) = wg(:,inode,ielem)

                  gnode = efn2efn(3,jnode,ielem)                     ! volumetric node of partner node

                  kelem = efn2efn(2,jnode,ielem)                     ! volumetric element of partner node

                  nx = facenodenormal(:,jnode,ielem)                 ! outward facing normal of facial node
  
                  call conserved_to_primitive(ughst(:,gnode),vg_Off(:),nequations)
                  call primitive_to_entropy(vg_Off(:),wg_Off(:),nequations)
    
                                                                     ! LDC/LDG penalty value
                  l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
                  dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg_On(:) - wg_Off(:))
                                                                     ! add LDC/LDG penalty to each physical gradient using the normal
                  do jdir = 1,ndim
                    phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
                  end do
  
                end do

                nghst_volume = nghst_volume + n_S_2d_On               !  Keep track of position in Ghost stack (n_S_2d_On=n_S_2d_Off)
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
              else                                              !       ON-Processor Contributions to gsat:  Conforming Interface
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                do i = 1,n_S_2d_On

                  jnode = n_S_2d_On*(iface-1) + i                      ! shell coordinate counter

                  inode = ifacenodes_On(jnode)                         ! corresponding volumetric node for face node

                  knode = efn2efn(1,jnode,ielem)                       ! volumetric node of partner node

                  kelem = efn2efn(2,jnode,ielem)                       ! volumetric element of partner node

                  nx = facenodenormal(:,jnode,ielem)                   ! outward facing normal of facial node
    
                  wg_On (:) = wg(:,inode,ielem)
                  wg_Off(:) = wg(:,knode,kelem)
                                                                       ! LDC/LDG penalty value
                  l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
                  dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg_On(:) - wg_Off(:))
                                                                       ! add LDC/LDG penalty to each physical gradient using the normal
                  do jdir = 1,ndim
                    phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
                  end do
                end do

              endif                                   
                                               !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                                               !       NON-CONFORMING INTERFACES:  polynomial orders do NOT match 
                                               !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            else if ((elem_props(2,ielem) /= ef2e(4,iface,ielem)) .and. (ef2e(9,iface,ielem) == 0)) then 

              n_S_1d_max  = (npoly_max+1)**1
              n_S_2d_max  = (npoly_max+1)**2
              kface       = ef2e(1,iface,ielem)
              kelem       = ef2e(2,iface,ielem)

              call element_properties(kelem,&
                             n_pts_1d=n_S_1d_Off,&
                             n_pts_2d=n_S_2d_Off,&
                             x_pts_1d=x_S_1d_Off,&
                           ifacenodes=ifacenodes_Off)

              n_S_1d_Mort = max(n_S_1d_On,n_S_1d_Off)
              n_S_2d_Mort = (n_S_1d_Mort)**2
              if(allocated(x_S_1d_Mort)) deallocate(x_S_1d_Mort) ; allocate(x_S_1d_Mort(n_S_1d_Mort)) ;
              if(allocated(w_S_1d_Mort)) deallocate(w_S_1d_Mort) ; allocate(w_S_1d_Mort(n_S_1d_Mort)) ;
              call Gauss_Legendre_points(n_S_1d_Mort,x_S_1d_Mort,w_S_1d_Mort)

              allocate(wg_2d_Mort_On (nequations,n_S_2d_Mort))
              allocate(wg_2d_Mort_Off(nequations,n_S_2d_Mort))

              allocate(wg_2d_On      (nequations,n_S_2d_On  ))
              allocate(wg_2d_Off     (nequations,n_S_2d_Off ))

              if(n_S_1d_Mort == n_S_1d_On) then
                poly_val = n_S_1d_Mort - npoly
                 allocate(Intrp_On (n_S_1d_On  ,n_S_1d_Mort)) ; 
                          Intrp_On (:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_On  ,1:n_S_1d_Mort,poly_val,1) ;
                 allocate(Extrp_On (n_S_1d_Mort,n_S_1d_On  )) ; 
                          Extrp_On (:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_On  ,poly_val,1) ;
                poly_val = n_S_1d_Off  - npoly
                 allocate(Extrp_Off(n_S_1d_Mort,n_S_1d_Off )) ; 
                          Extrp_Off(:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_Off ,poly_val,2) ;
              else
                poly_val = n_S_1d_On - npoly
                 allocate(Intrp_On (n_S_1d_On  ,n_S_1d_Mort)) ; 
                          Intrp_On (:,:) = Restrct_Gau_2_LGL_1d(1:n_S_1d_On  ,1:n_S_1d_Mort,poly_val,2) ;
                 allocate(Extrp_On (n_S_1d_Mort,n_S_1d_On  )) ; 
                          Extrp_On (:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_On  ,poly_val,2) ;
                poly_val = n_S_1d_Mort - npoly
                 allocate(Extrp_Off(n_S_1d_Mort,n_S_1d_Off )) ; 
                          Extrp_Off(:,:) = Prolong_LGL_2_Gau_1d(1:n_S_1d_Mort,1:n_S_1d_Off ,poly_val,1) ;
              endif

                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
              if (ef2e(3,iface,ielem) /= myprocid) then         !       Parallel Contributions to gsat:  Non-Conforming Interface
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                Off_Elem_0:do k = 1, n_S_2d_Off
    
                  ug_Off(:) = ughst(:,  nghst_volume + k)            ! conserved variable    in Petsc ghost registers

                  call conserved_to_primitive(ug_Off(:),vg_Off(:),nequations)
                  call primitive_to_entropy  (vg_Off(:),wg_Off(:),nequations)

                  wg_2d_Off(:,k) = wg_Off(:)
       
                enddo Off_Elem_0

                nghst_volume = nghst_volume + n_S_2d_Off                !  Keep track of position in Ghost stack (n_S_2d_On=n_S_2d_Off)

                call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                                       wg_2d_Off(:,:),wg_2d_Mort_Off(:,:),Extrp_Off,Extrp_Off)

                On_Mortar_0:do j = 1, n_S_2d_Mort
          
                  jnode =  n_S_2d_max*(iface-1) + j                     ! Index in facial ordering
          
                  l = efn2efn_Gau(3,jnode,ielem) - nghst_Mortar
        
                  wg_2d_Mort_On(:,j) = wg_2d_Mort_Off(:,l)

                enddo On_Mortar_0

                nghst_Mortar = nghst_Mortar + n_S_2d_Mort
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
              else                                              !       Serial Contributions to gsat:  Non-Conforming Interface
                                                                !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                Off_Elem_1:do k = 1, n_S_2d_Off
    
                  lnode =  n_S_2d_Off*(kface-1) + k                     ! Index in facial ordering
       
                  knode = ifacenodes_Off(lnode)                         ! Volumetric node index corresponding to facial node index
       
                  vg_Off(:)      = vg(:,knode,kelem)
                  call primitive_to_entropy(vg_Off,wg_Off,nequations)
       
                  wg_2d_Off(:,k) = wg_Off(:)
       
                enddo Off_Elem_1


                call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Off,n_S_1d_Mort,x_S_1d_Off,x_S_1d_Mort, &
                                       wg_2d_Off(:,:),wg_2d_Mort_Off(:,:),Extrp_Off,Extrp_Off)

                On_Mortar_1:do j = 1, n_S_2d_Mort
          
                  jnode =  n_S_2d_max*(iface-1) + j                     ! Index in facial ordering
          
                  l = efn2efn_Gau(4,jnode,ielem) - n_S_2d_max*(kface-1) ! Index in off-element local facial ordering
    
                  wg_2d_Mort_On(:,j) = wg_2d_Mort_Off(:,l)
        
                enddo On_Mortar_1

              endif

              call ExtrpXA2XB_2D_neq(nequations,n_S_1d_Mort,n_S_1d_On,x_S_1d_Mort,x_S_1d_On, &
                                     wg_2d_Mort_On(:,:),wg_2d_On(:,:),Intrp_On,Intrp_On)

              On_Elem:do i = 1,n_S_2d_On                                  ! On-Element face loop
    
                jnode = n_S_2d_On*(iface-1) + i                           ! Index in facial ordering

                inode = ifacenodes_On(jnode)                              ! corresponding volumetric node for face node
    
                nx = facenodenormal(:,jnode,ielem)                        ! outward facing normal of facial node

                vg_On(:) =   vg(:,inode,ielem)
                call primitive_to_entropy(vg_On,wg_On,nequations)
                                                                          ! LDC/LDG penalty value
                l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
                dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg_On(:) - wg_2d_On(:,i))
    
                                                                          ! add LDC/LDG penalty to each physical gradient using the normal
                do jdir = 1,ndim
                  phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
                end do

              end do On_Elem

              deallocate(x_S_1d_Mort, w_S_1d_Mort)
              deallocate(wg_2d_Mort_On, wg_2d_Mort_Off)
              deallocate(wg_2d_On, wg_2d_Off)
              deallocate(Intrp_On,Extrp_On, Extrp_Off)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       NON-CONFORMING h refinement: note that ghost faces are skipped over since ef2e(9,iface,ielem) = -1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            elseif(ef2e(9,iface,ielem) == 1)then

            end if          !  face cycle for BC / Off / On-process face-types

          end do face_Loop

        end do element_Loop

        deallocate(dphi,vg_On,vg_Off,wg_On,wg_Off)
        deallocate(phig_tmp1,phig_tmp2,phig_err1)

      end subroutine viscous_gradients

!=========================================================================================

      subroutine viscous_gradients_New()

        use variables
        use referencevariables
        use nsereferencevariables
        use collocationvariables, only: nnzgrad,iagrad,jagrad,dagrad,pinv,l10, &
                                      & ldg_flip_flop_sign, alpha_ldg_flip_flop
        use initcollocation,      only: element_properties

        implicit none

        ! temporary arrays for phi and delta phi
        real(wp), dimension(:),    allocatable :: dphi, vg_Off, wg_On, wg_Off
        real(wp), dimension(:,:),  allocatable :: phig_tmp1, phig_tmp2, phig_err1
        real(wp), dimension(:,:,:),allocatable :: phig_tmp3, phig_test

        ! normal direction at face
        real(wp) :: nx(3)

        ! loop indices
        integer :: ielem, jdir, idir, i
        integer :: inode, jnode, knode, gnode
        integer :: kelem
        integer :: iface

        real(wp) :: l10_ldg_flip_flop

        integer :: jdir_face
        real(wp) :: dx,t1

        ! phig_tmp is calculated in computational space
        allocate(phig_tmp1(nequations,ndim),phig_tmp2(nequations,ndim),phig_err1(nequations,ndim))
        ! tmp variables: dphi, vg_Off, wg_On and wg_Off
        allocate(dphi(nequations),vg_Off(nequations),wg_On(nequations),wg_Off(nequations))
        allocate(phig_tmp3(nequations,ndim,nodesperelem))
        allocate(phig_test(nequations,ndim,nodesperelem))

        ! loop over all elements
        do ielem = ihelems(1), ihelems(2)

          call element_properties(ielem,              &
                                n_pts_2d=nodesperface,&
                                n_pts_3d=nodesperelem,&
                                 nnzgrad=nnzgrad,     &
                                  iagrad=iagrad,      &
                                  jagrad=jagrad,      &
                                  dagrad=dagrad,      &
                              ifacenodes=ifacenodes)

          ! compute computational gradients of the entropy variables
          ! initialize phi
          phig    (:,:,:,ielem) = 0.0_wp
          phig_err(:,:,:,ielem) = 0.0_wp
          grad_w_jacobian(:,:,:,ielem) = 0.0_wp
          ! loop over every node in element
          do inode = 1, nodesperelem

            phig_tmp1(:,:) = 0.0_wp ; phig_tmp2(:,:) = 0.0_wp ;

            ! loop over number of dependent elements in gradient
            do i = iagrad(inode), iagrad(inode+1)-1
              ! loop over dimensions
              do jdir = 1,ndim
                ! column/node from gradient operator in CSR format in
                ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
                jnode = jagrad(jdir,i)
                ! update gradient using coefficient and entropy variables at appropriate node
                phig_tmp1(:,jdir) = phig_tmp1(:,jdir) + dagrad(jdir,i)*wg(:,jnode,ielem) 
                phig_tmp2(:,jdir) = phig_tmp2(:,jdir) + dagrad(jdir,i)*vg(:,jnode,ielem) 
              end do
            end do
  
            phig_tmp3(:,:,inode) = phig_tmp1(:,:) 

            ! Store gradient of the entropy variables in computational space
            grad_w_jacobian(:,:,inode,ielem) = phig_tmp1(:,:)

            ! transform to physical space using dxi_jdir/dx_idir
            do jdir = 1,ndim
              do idir = 1,ndim
                phig(:,idir,inode,ielem) = phig(:,idir,inode,ielem) + phig_tmp1(:,jdir)*r_x(jdir,idir,inode,ielem)
                phig_err1(:,idir)        = phig_err1(:,idir)        + phig_tmp2(:,jdir)*r_x(jdir,idir,inode,ielem)
              end do
            end do

            phig_err(:,:,inode,ielem) = matmul(dWdV(vg(:,inode,ielem),nequations),phig_err1(:,:))  &
                                      - phig(:,:,inode,ielem)
  
          end do
           
          ! LDC/LDG penalty on phig
          ! 
          ! loop over only faces
          do iface = 1, nfacesperelem

            ! sign so normal is facing outward
            dx = sign(1.0_wp,real(facenormalcoordinate(iface),wp))
            jdir_face = abs(facenormalcoordinate(iface))

            if (ef2e(1,iface,ielem) < 0) then
              cycle
            else if (ef2e(3,iface,ielem) /= myprocid) then
              do i = 1,nodesperface
                jnode = nodesperface*(iface-1)+i
                ! corresponding volumetric node for face node
                inode = ifacenodes(jnode)
                ! volumetric node of partner node
                gnode = efn2efn(3,jnode,ielem)
                ! volumetric element of partner node
                kelem = efn2efn(2,jnode,ielem)
                ! outward facing normal of facial node
                nx = facenodenormal(:,jnode,ielem)
  
                wg_On(:) = wg(:,inode,ielem)
  
                call conserved_to_primitive(ughst(:,gnode),vg_Off(:),nequations)
                call primitive_to_entropy(vg_Off(:),wg_Off(:),nequations)
  
                ! LDC/LDG penalty value
                l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
                dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg_On(:) - wg_Off(:))
                ! add LDC/LDG penalty to each physical gradient using the normal
                do jdir = 1,ndim
                  phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
                end do

                phig_tmp3(:,jdir_face,inode) = phig_tmp3(:,jdir_face,inode) + dphi(:)*dx

              end do
            else
              do i = 1,nodesperface
                jnode = nodesperface*(iface-1)+i
                ! corresponding volumetric node for face node
                inode = ifacenodes(jnode)
                ! volumetric node of partner node
                knode = efn2efn(1,jnode,ielem)
                ! volumetric element of partner node
                kelem = efn2efn(2,jnode,ielem)
                ! outward facing normal of facial node
                nx = facenodenormal(:,jnode,ielem)
  
                wg_On (:) = wg(:,inode,ielem)
                wg_Off(:) = wg(:,knode,kelem)
  
                ! LDC/LDG penalty value
                l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
                dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg_On(:) - wg_Off(:))
                ! add LDC/LDG penalty to each physical gradient using the normal
                do jdir = 1,ndim
                  phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
                end do

                phig_tmp3(:,jdir_face,inode) = phig_tmp3(:,jdir_face,inode) + dphi(:)*dx

              end do
            end if
          end do

          do inode = 1, nodesperelem
            phig_test(:,:,inode) = 0.0_wp
            do jdir = 1,ndim
              do idir = 1,ndim
                phig_test(:,idir,inode) = phig_test(:,idir,inode) + phig_tmp3(:,jdir,inode)*r_x(jdir,idir,inode,ielem)
              end do
            end do
          end do
          t1 = maxval(abs(phig_test(:,:,:) - phig(:,:,:,ielem)))
          if(t1 >= 1.0e-15) write(*,*)'error in reconcilestates',t1
!         This line runs the new path through the code
          phig(:,:,:,ielem) = phig_test(:,:,:)

        end do 

        deallocate(phig_tmp1,phig_tmp2,phig_err1)
        deallocate(dphi,vg_Off,wg_On,wg_Off)
        deallocate(phig_tmp3,phig_test)

      end subroutine viscous_gradients_New

 end module navierstokes
