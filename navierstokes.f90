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
    use variables, only: vg, wg, ug, xg, r_x, boundaryelems, omega, phig, &
                         specific_entropy, mean_vg, time_ave_prod_vel_comp, &
                         reynolds_stress, mut
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
    integer :: inode,ielem
    ! High and low element indices for volume elements
    integer :: iell, ielh
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

    ! Arbitrary value for normals
    ctmp = 0.0_wp

    ! Low volumetric element index
    iell = ihelems(1)
    ! High volumetric element index
    ielh = ihelems(2)

    ! We are limited to the calorically perfect Navier-Stokes equations
    nequations = 5
    ! Total degrees of freedom
    ndof = nequations*nnodes

    ! Set the flow parameters
    call set_Flow_parameters()

    ! Set initial and boundary condition procedure
    call set_InitialCondition(InitialCondition)

    ! Allocate memory for flow state data
    !
    ! Conserved variables (denisty, momentum, total energy)
    allocate(ug(1:nequations,nodesperelem,iell:ielh))
    ug = 0.0_wp
    
    ! Primitive variables (density, velocity, temperature)
    allocate(vg(1:nequations,nodesperelem,iell:ielh))
    vg = 0.0_wp
    
    ! Entropy variables (see primitive_to_entropy subroutine)
    allocate(wg(1:nequations,nodesperelem,iell:ielh))
    wg = 0.0_wp
    
    ! Vorticity field (\nabla \times velocity)
    allocate(omega(1:3,nodesperelem,iell:ielh))
    omega = 0.0_wp

    ! Entropy field
    allocate(specific_entropy(nodesperelem,iell:ielh))
    specific_entropy = 0.0_wp

    allocate(mut(nodesperelem,iell:ielh))
    mut = 0.0_wp

    ! Allocate memory if time averaging is required
    if (time_averaging) then

      ! Time-average of primitive variables
      allocate(mean_vg(1:nequations,nodesperelem,iell:ielh))
      mean_vg = 0.0_wp

      ! Time-average of the product of the velocity components
      allocate(time_ave_prod_vel_comp(6,nodesperelem,iell:ielh))
      time_ave_prod_vel_comp = 0.0_wp

      ! Reynolds stresses
      allocate(reynolds_stress(6,nodesperelem,iell:ielh))
      reynolds_stress = 0.0_wp
    
    endif

    ! We have a different output unit for each equation.
    ! if this were a real code we would replace this with
    ! something more useful.
    allocate(iunit(nequations))

    ! For our call to InitialSubroutine  we need
    ! the viscous fluxes
    allocate(fvtmp(nequations))
    allocate(phitmp(nequations,ndim))

    if (new .eqv. .true.) then
      ! Loop over elements
      do ielem = iell, ielh
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



      !1000 FORMAT(8(e25.15,1x))
      !stop
    else

      ! Read solution from restart file
      call read_restart_file()
 
      ! Loop over elements
      do ielem = iell, ielh
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
    use mpimod, only: PetscComm1DDataSetup, PetscComm1DDataSetupWENO,    &
      & PetscComm1DElementDataSetup, PetscComm2DDataSetup, PetscComm2DGeomDataSetup, &
      & PetscComm1DDataSetupWENOGeom, PetscCommShellDataSetup

    use nsereferencevariables, only : viscous
    
    ! Nothing is implicitly defined
    implicit none

    integer :: nshell

    ! Setup the PETSc parallel communication for exchanging the conservative
    ! variables at the parallel face element 
    call PetscComm1DDataSetup(ug,ughst,upetsc,ulocpetsc,nequations, &
      & nodesperelem,nelems,nghost)

    ! Setup the PETSc parallel communication for exchanging the conservative
    ! variables at the parallel face element 
    if(discretization == 'SSWENO')then

      nshell = nfacesperelem*nodesperface

      call PetscComm1DDataSetupWENO(ug,ughstWENO,upetscWENO,ulocpetscWENO,nequations, &
      & nodesperelem,nelems,nghost)

      call PetscComm1DDataSetupWENOGeom(xgWENO_partner,xghstWENO_partner,xpetscWENO_partner, &
      & xlocpetscWENO_partner,3,nshell,nelems,nghost)

      call PetscCommShellDataSetup(ugWENO_partner,ughstWENO_partner,upetscWENO_Shell, &
      & ulocpetscWENO_Shell,nequations,nshell,nelems,nghost)

    endif

    ! Setup the PETSc parallel communication for exchanging the gradient of the 
    ! entropy variables at the parallel face element 
    call PetscComm2DDataSetup(phig,phighst,phipetsc,philocpetsc,nequations,3, &
      & nodesperelem,nelems,nghost)

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

    return
  end subroutine nse_communicationsetup

  !============================================================================

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
    use nsereferencevariables, only: gm1M2, gm1og
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! conserved variables
    real(wp), intent(out) :: uout(nq)

    ! density
    uout(1) = vin(1)
    ! momentum
    uout(2:4) = vin(1)*vin(2:4)
    ! energy
    uout(5) = vin(1)*( (1.0_wp-gm1og)*vin(5) &
      + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) )

    return
  end subroutine primitive_to_conserved

  pure subroutine conserved_to_primitive(uin,vout,nq)
    ! this routine calculates the primitive variables
    ! from the conserved variables
    use nsereferencevariables, only: gm1M2, gm1og
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! conserved variables
    real(wp), intent(in) :: uin(nq)
    ! primitive variables
    real(wp), intent(out) :: vout(nq)

    ! density
    vout(1) = uin(1)
    ! velocity
    vout(2:4) = uin(2:4)/uin(1)
    ! temperature
    vout(5) = ( uin(5)/uin(1) - gm1M2*0.5_wp*dot_product(vout(2:4),vout(2:4)) )/(1.0_wp-gm1og)

    return
  end subroutine conserved_to_primitive

  pure function conserved_to_primitive_F(uin,nq)

    use nsereferencevariables, only: gm1og, gm1M2

    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: uin(nq)

    ! output primitive variables
    real(wp) :: conserved_to_primitive_F(nq)

    ! density
    conserved_to_primitive_F(1) = uin(1)
    ! velocity
    conserved_to_primitive_F(2:4) = uin(2:4)/uin(1)
    ! temperature
    conserved_to_primitive_F(5) = ( uin(5)/uin(1) &
      - gm1M2*0.5_wp*dot_product(conserved_to_primitive_F(2:4),conserved_to_primitive_F(2:4)) )/(1.0_wp-gm1og)

    return
  end function conserved_to_primitive_F

  pure function conserved_to_entropy(uin,nq)

    use nsereferencevariables, only: gm1og, gm1M2, gamma0

    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: uin(nq)

    real(wp)             :: vin(nq)

    real(wp)             :: conserved_to_entropy(nq)
    real(wp)             :: Tinv

    ! density
    vin(1) = uin(1)
    ! velocity
    vin(2:4) = uin(2:4)/uin(1)
    ! temperature
    vin(5) = ( uin(5)/uin(1) - gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) ) * gamma0

    Tinv = 1.0_wp/vin(5)

    ! w_1 = h/T - s - (gamma_0 - 1) M_0^2 u_k u_k/(2T)
    conserved_to_entropy(1) = 1.0_wp-0.5_wp*gm1M2*dot_product(vin(2:4),vin(2:4)) * Tinv &
      -specificentropy(vin,nq)
    ! w_{k+1} = (gamma_0 - 1) M_0^2 u_k/T, k = 1,2,3
    conserved_to_entropy(2:4) = gm1M2*vin(2:4) * Tinv
    ! w_5 = -1/T
    conserved_to_entropy(5) = - Tinv

    return

  end function conserved_to_entropy

  pure function specificentropy(vin,nq)
    ! this function calculates the specific thermodynamic
    ! entropy using the primitive variable vector
    use nsereferencevariables, only: gm1og
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: vin(nq)

    ! specific gas constant
    real(wp) :: rs

    ! output thermodynamic specific entropy
    real(wp) :: specificentropy

    rs = 1.0_wp

    specificentropy = (1.0_wp-gm1og)*rs*log(vin(5)) &
      - gm1og*rs*log(vin(1))

    return
  end function specificentropy

  pure subroutine primitive_to_entropy(vin,wout,nq)
    ! this routine calculates the entropy variables corresponding
    ! to the entropy--entropy flux pair (S,F^i) = (-rho*s,-rho*u_i*s)
    ! using the primitive variables.
    use nsereferencevariables, only:gm1M2, gm1og
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

  pure subroutine KineticEnergyVariables(vin,keout,nq)
    ! this routine calculates the kinetic energy variables corresponding to
    ! the entropy function keout = 1/2 r (u_i u_i)
    ! using the primitive variables.
    implicit none
    ! number of equations
    integer, intent(in)   :: nq
    ! primitive variables
    real(wp), intent(in)  :: vin(nq)
    ! entropy variables
    real(wp), intent(out) :: keout(nq)

    ! ke_1 = -1/2 u_i u_i 
    keout(1)   = - 0.5_wp*dot_product(vin(2:4),vin(2:4))
    ! ke_{2:4} = u_i
    keout(2:4) = vin(2:4)
    ! no fifth variable because KE not dependent on rho E .
    keout(5)   = 0.0_wp

    return
  end subroutine KineticEnergyVariables

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
    ! by the reference pressure that there is a nondimensional parameter
    ! in the flux.
    use nsereferencevariables, only: gm1M2, gm1og, gM2
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! normal vector
    real(wp), intent(in) :: nx(3)
    ! primitive variables
    real(wp), intent(in) :: vin(nq)

    ! output normal flux
    real(wp) :: normalflux(nq)

    ! local variables
    !
    ! normal mass flux
    real(wp) :: un
    ! pressure
    real(wp) :: p

    ! calculate normal mass flux
    un = vin(1)*dot_product(vin(2:4),nx)
    ! calculate pressure (R = 1)
    p = vin(1)*vin(5)
    ! mass flux (\rho u \cdot n)
    normalflux(1) = un
    ! momentum flux (\rho u \cdot n) u + p n /(gamma_0 M_0^2)
    normalflux(2:4) = un*vin(2:4) + p*nx/gM2
    ! energy flux (\rho u \cdot n) H
    normalflux(5) = un*( vin(5) &
      + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) )

    return
  end function normalflux

!===================================================================================================

  pure function normalviscousflux(vin,phi,nx,nq,mut)
    ! this function calculates the viscous flux in 
    ! the normal direction based on the primitive
    ! variables and the gradients of the entropy variables,
    ! penalized with an LDC/LDG methodology. The use
    ! of entropy variables ensures stability.
    use nsereferencevariables, only: gm1M2, gm1og, gM2, Re0inv, Pr0
    
    use controlvariables, only : variable_viscosity

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

    ! thermal conductivity (normalized by kappa0)
    real(wp) :: kappa
    ! dynamic viscosity (normalized by mu0)
    real(wp) :: mu

    continue

    fvl(:,:) = ViscousFlux(vin,phi,nq,mut)

    normalviscousflux = 0.0_wp

    do idir = 1,3
      normalviscousflux = normalviscousflux + Re0inv*fvl(:,idir)*nx(idir)
    end do

    return
  end function normalviscousflux

!===================================================================================================

  pure function viscousflux3D(vin,phi,Jx,nq,nd,mut)
    ! this function calculates the viscous flux in 
    ! the three computational space directions based on the primitive
    ! variables and the gradients of the entropy variables,
    ! penalized with an LDC/LDG methodology. The use
    ! of entropy variables ensures stability.
    use nsereferencevariables, only: Re0inv
    
    use controlvariables, only : variable_viscosity

    ! Nothing is implicitly defined
    implicit none
    
    integer,  intent(in) :: nq, nd
    ! contravariant vector
    real(wp), intent(in) :: Jx(3,3)
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
        viscousflux3D(:,jdir) = viscousflux3D(:,jdir) + Re0inv*fvl(:,idir)*Jx(jdir,idir)
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

  pure function EntropyConsistentFlux(vl,vr,Jx,neqin)
    ! this function calculates the normal entropy consistent
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the nondimensionalization employed
    ! herein and follows directly from the work of Ismail and Roe,
    ! DOI: 10.1016/j.jcp.2009.04.021
    use nsereferencevariables, only: gM2, gM2I, gm1og, gp1og, gm1M2, gamma0
    implicit none
    ! Arguments
    ! =========
    ! number of equations
    integer,  intent(in)                   :: neqin
    ! left and right states
    real(wp), intent(in), dimension(neqin) :: vl, vr
    ! metrics scaled by jacobian
    real(wp), intent(in), dimension(3)     :: Jx
  
    ! Function
    ! ========
    real(wp), dimension(neqin) :: EntropyConsistentFlux
  
    ! Local Variables
    ! ===============
    ! temporary variables (dimension is legacy from different code)
    real(wp), dimension(neqin) :: vhat

    ! inverse and square roots of temperature states
    real(wp) :: root_Tl_I, root_Tr_I, root_Tl, root_Tr
    ! average of inverse temperature and its inverse
    real(wp) :: tinvav, tinvavinv

    ! prevent division by zero
    real(wp), parameter :: sdiv = 1e-015_wp, sdiv2 = 1e-030_wp
    ! log degenerates in multiple cases (see reference for use)
    real(wp), parameter :: seps = 1.0e-04_wp
    real(wp) :: xi, gs, us, ut, ftmp
    real(wp) :: s1, s2, al

    ! normal mass flux (mdot), Pressure,  Temperature
    real(wp) :: mdot, P, T

    ! Sqrt[temperature] and Sqrt[temperature]^{-1} are used a lot
    root_Tl   = sqrt(vl(5)) ; root_Tr   = sqrt(vr(5))
    root_Tl_I = one/root_Tl ; root_Tr_I = one/root_Tr
    tinvav    = root_Tl_I   + root_Tr_I
    tinvavinv = one/tinvav
  
    ! velocity
    vhat(2:4) = (vl(2:4)*root_Tl_I + vr(2:4)*root_Tr_I)*tinvavinv

    ! pressure
    ftmp = vl(1)*root_Tl + vr(1)*root_Tr
    P    = ftmp*tinvavinv

    ! logarithmic averages used in density and temperature
    ! calculate s1
    xi = (root_Tl*vl(1))/(root_Tr*vr(1))
    gs = (xi-one)/(xi+one)
    us = gs*gs
    al = exp(-us/seps)
    s1 = half / (one + us*(third + us*(fifth + us*(seventh + us*ninth) ) ) )
    ut = log(xi)
    s1 = ftmp * (al * s1 + (one-al) * gs * ut / (ut*ut + sdiv2) )

    ! calculate s2
    xi = root_Tl_I/root_Tr_I
    gs = (xi-one)/(xi+one)
    us = gs*gs
    al = exp(-us/seps)
    s2 = half / (one + us*(third + us*(fifth + us*(seventh + us*ninth) ) ) )
    ut = log(xi)
    s2 = tinvav * (al * s2 + (one-al) * gs * ut / (ut*ut + sdiv2) )

    ! density
    vhat(1) = half *tinvav*s1

    ! temperature
    T = half * (gm1og * P + gp1og*s1/s2) / vhat(1)

    ! total enthalpy
    vhat(5) = gm1M2*half*dot_product(vhat(2:4),vhat(2:4)) + T

    ! normal mass flow rate
    mdot = vhat(1)*dot_product(vhat(2:4),Jx)

    EntropyConsistentFlux(1) = mdot
    EntropyConsistentFlux(2) = mdot*vhat(2) + Jx(1) * P * gM2I
    EntropyConsistentFlux(3) = mdot*vhat(3) + Jx(2) * P * gM2I
    EntropyConsistentFlux(4) = mdot*vhat(4) + Jx(3) * P * gM2I
    EntropyConsistentFlux(5) = mdot*vhat(5)

!   !  Check if it satisfies the entropy shuffle   DelW^T f = Delpsi
!   call primitive_to_entropy(vL,wL,neqin)
!   call primitive_to_entropy(vR,wR,neqin)

!   err = dot_product(wL-wR,EntropyConsistentFlux) - gm1og* (           &
!                   ((vl(2)-vR(2))*vhat(1) + (vl(1)-vR(1))*vhat(2))*Jx(1) +  &
!                   ((vl(3)-vR(3))*vhat(1) + (vl(1)-vR(1))*vhat(3))*Jx(2) +  &
!                   ((vl(4)-vR(4))*vhat(1) + (vl(1)-vR(1))*vhat(4))*Jx(3) )
!   if(abs(err) >= 1.0e-05_wp) write(*,*)'error in Ismail_Roe',err

    return
  end function EntropyConsistentFlux

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

    use nsereferencevariables, only: gM2, gM2I, gm1og, gp1og, gm1M2, gamma0, gamI
    implicit none
    ! Arguments
    ! =========
    ! number of equations
    integer,  intent(in)                :: nq
    ! left and right states
    real(wp), intent(in), dimension(nq) :: vl, vr
    ! metrics scaled by jacobian
    real(wp), intent(in), dimension(3)  :: Jx
  
    ! Function
    ! ========
    real(wp), dimension(nq) :: Entropy_KE_Consistent_Flux
  
    ! Local temporary Variables
    ! ===============
    real(wp), dimension(nq) :: vave
    real(wp), dimension(nq) :: wL,wR

    ! prevent division by zero and 
    real(wp), parameter :: sdiv2 = 1e-030_wp
    ! log degenerates in multiple cases (toggle between asymptotic expression and actual logarithm)
    real(wp), parameter :: seps = 1.0e-02_wp
    real(wp) :: xi, fs, us, ut, al, F
    real(wp) :: alr,alB
    real(wp) :: err


    ! normal mass flux (mdot), Pressure, Logave density and temperature
    real(wp) :: mdot, P, rhotil, Btil, bL, bR, Bave, Keave

    ! Average Density & velicity
    vave(1:5) = half * (vl(1:5) + vr(1:5))

    ! average inverse Temperature (defined as B) and Pressure
    bL = 1.0_wp/vl(5) ; bR = 1.0_wp/vr(5) ; Bave   = half * (bL + bR)
    P  = gM2I * vave(1) / Bave

    ! Average KE 
    ! Keave = half*half*(dot_product(vL(2:4),vL(2:4)) + dot_product(vR(2:4),vR(2:4)))
      Keave = half*half*(vL(2)*vL(2)+vL(3)*vL(3)+vL(4)*vL(4) + vR(2)*vR(2)+vR(3)*vR(3)+vR(4)*vR(4))

    ! logarithmic average of density
    xi = vl(1)/vr(1)
    fs = (xi-one)/(xi+one)
    us = fs * fs
    alr= exp(-us/seps)
    F  = half / (one + us*(third + us*(fifth + us*(seventh + us*ninth) ) ) )
!   F  = (one + us*(third + us*(fifth + us*(seventh + us*ninth) ) ) )
    ut = log(xi)

    rhotil = 2.0_wp * vave(1) * (alr * F + (one-alr) * fs * ut / (ut*ut + sdiv2) )
!   rhotil = vave(1) / (al * F + half*(one-al) * fs * ut / (fs*fs + sdiv2) )

    ! logarithmic average of  B = 1 / temperature
    xi = bL/bR
    fs = (xi-one)/(xi+one)
    us = fs * fs
    alB= exp(-us/seps)
    F  = half / (one + us*(third + us*(fifth + us*(seventh + us*ninth) ) ) )
!   F  = (one + us*(third + us*(fifth + us*(seventh + us*ninth) ) ) )
    ut = log(xi)

    Btil   = 2.0_wp * Bave * (alB * F + (one-alB) * fs * ut / (ut*ut + sdiv2) )
!   Btil   = Bave / (al * F + half*(one-al) * fs * ut / (fs*fs + sdiv2) )

    ! normal mass flow rate
    ! mdot = rhotil*dot_product(vave(2:4),Jx)
      mdot = rhotil*(vave(2)*Jx(1)+vave(3)*Jx(2)+vave(4)*Jx(3))

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
!   if(abs(err) >= 2.0e-05_wp) then
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

    use nsereferenceVariables, only: gm1og, gm1M2,gM2I,gamI,gM2

      implicit none
      integer,                 intent(in) :: nq
      real(wp),                intent(in) :: nx(3)
      real(wp), dimension(nq), intent(in) :: vL,vR
      real(wp), dimension(nq), intent(in) :: wL,wR
      real(wp), parameter                 ::    ep = 1.0e-10_wp
      real(wp), parameter                 :: third = 1.0_wp/3.0_wp

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
      real(wp) :: pL, pR, pA, tA, En, keL, keR, den
      real(wp) :: norm_x

      integer  :: i,j

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

  pure function dUdV(Vin,neqin)
    ! Checked MHC 08_09_13
    ! this function calculates the jacobian of the
    ! conserved variables with respect to the primitive
    ! variables
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neqin
    ! primitive variables
    real(wp), intent(in) :: Vin(neqin)

    ! output jacobian
    real(wp) :: dUdV(neqin,neqin)

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

  pure function dVdU(Vin,neqin)
    ! Checked MHC 08_09_13
    ! this function calculates the jacobian of the
    ! primitive variables with respect to the conserved
    ! variables
    use nsereferencevariables
    implicit none
    ! number of equations
    integer, intent(in) :: neqin
    ! primitive variables
    real(wp), intent(in) :: Vin(neqin)

    ! output jacobian
    real(wp) :: dVdU(neqin,neqin)

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
 !! subroutine CharacteristicDecomp(Vav,nq,Le,Re,ev,Jx)
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

    integer :: iell, ielh
    integer :: nshell

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! ghost cells for solution
    allocate(ughst(1:nequations,nghost)) ; ughst = 0.0_wp
    
    if(discretization == 'SSWENO') then
      allocate(ughstWENO(1:nequations,nghost))             ; ughstWENO          = 0.0_wp
      allocate(ughstWENO_partner(1:nequations,nghost))     ; ughstWENO_partner  = 0.0_wp

      nshell  = nfacesperelem*nodesperface
      allocate(ugWENO_partner(1:nequations,nshell,iell:ielh))  ; ugWENO_partner = 0.0_wp
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
        allocate(r_x_ghst(3,3,nghost))
        r_x_ghst = 0.0_wp
      endif
    endif

    ! stores solution at T^n only used if step is rejected 
    allocate(uold(1:nequations,nodesperelem,iell:ielh))
    uold = 0.0_wp
    ! used in RK error estimation
    allocate(uhat(1:nequations,nodesperelem,iell:ielh))
    uhat = 0.0_wp


    ! penalty terms
    allocate(gsat(1:nequations,nodesperelem,iell:ielh))
    gsat = 0.0_wp
    ! flux divergence vectors -- one for each direction
    allocate(divf(1:nequations,ndim,nodesperelem,iell:ielh))
    divf = 0.0_wp

    ! entropy flux divergence used for error estimation
    allocate(divf_S(ndim,nodesperelem,iell:ielh))
    divf_S = 0.0_wp

    ! convective flux in each computational direction
    allocate(fg(1:nequations,ndim,nodesperelem,iell:ielh))
    fg = 0.0_wp
    ! viscous flux in each computational direction
    allocate(fvg(1:nequations,ndim,nodesperelem,iell:ielh))
    fvg = 0.0_wp
    ! variable gradients used for LDC/LDG approximation
    allocate(phig(1:nequations,3,nodesperelem,iell:ielh))
    phig = 0.0_wp
    allocate(phig_err(1:nequations,3,nodesperelem,iell:ielh))
    phig_err = 0.0_wp

    ! Gradient of the entropy variables in computational space for the
    ! calculation of the residual Jacobian matrix
    allocate(grad_w_jacobian(1:nequations,3,nodesperelem,iell:ielh))
    grad_w_jacobian = 0.0_wp
    
    ! ghost points for LDC/LDG approximation
    allocate(phighst(1:nequations,3,nghost))
    phighst = 0.0_wp
    ! shock sensor
    allocate(chig(nodesperelem,iell:ielh))
    chig = 0.0_wp
    ! artificial viscosity

    select case(RK_Method)

    case('Williamson_Low_Storage_45')

      ! Stores update
      allocate(du(1:nequations,nodesperelem,iell:ielh))
      du   = 0.0_wp
      ! local time derivative of conserved variables
      allocate(dudt(1:nequations,nodesperelem,iell:ielh))
      dudt = 0.0_wp

      ! local time derivative of entropy equation
      allocate(dudt_S(nodesperelem,iell:ielh))
      dudt_S = 0.0_wp

    case('Kraaij_LS_RK_35')

      ! Stores update
      allocate(du(1:nequations,nodesperelem,iell:ielh))
      du   = 0.0_wp
      ! local time derivative of conserved variables
      allocate(dudt(1:nequations,nodesperelem,iell:ielh))
      dudt = 0.0_wp

    case ('heun_method')
      ! Stores update
      allocate(du(1:nequations,nodesperelem,iell:ielh))
      du   = 0.0_wp
      ! local time derivative of conserved variables
      allocate(dudt(1:nequations,nodesperelem,iell:ielh))
      dudt = 0.0_wp

    case('IMEX_RK_46')
      allocate(Fimp(1:nequations,nodesperelem,iell:ielh,6))
      allocate(Fexp(1:nequations,nodesperelem,iell:ielh,6))
      allocate(uexp(1:nequations,nodesperelem,iell:ielh))
      allocate(non_lin_res(1:nequations,nodesperelem,iell:ielh))
    case('IMEX_RK_34')
      allocate(Fimp(1:nequations,nodesperelem,iell:ielh,4))
      allocate(Fexp(1:nequations,nodesperelem,iell:ielh,4))
      allocate(uexp(1:nequations,nodesperelem,iell:ielh))
      allocate(non_lin_res(1:nequations,nodesperelem,iell:ielh))
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
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,l10
    implicit none

    ! loop indices
    integer :: inode,ielem, jdir, idir
    integer :: jnode
    integer :: i

    ! low and high volumetric element indices
    integer :: iell, ielh
    ! LDC/LDG coefficient
    ! temporary arrays for phi and delta phi
    real(wp), allocatable :: phitmp(:,:)
    real(wp) :: theta, omega_loc(3), omegamag, sos, lmag
    real(wp) :: eta, delta
    real(wp), allocatable :: lambda(:), muvec(:)

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! phitmp is calculated in computational space
    allocate(phitmp(nequations,3))
    ! dphi is calculated at faces
    mut = 0.0_wp
    chig = 0.0_wp
    allocate(lambda(nequations))
    allocate(muvec(nequations))
    ! LDC/LDG coefficient according to CNG-revisited
!   l10 = -0.5_wp
    ! loop over elements
    do ielem = iell, ielh
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
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol,l10, &
                                  & ldg_flip_flop_sign, alpha_ldg_flip_flop, &
                                  & face_sign

    use mpimod, only: UpdateComm1DGhostData, UpdateComm2DGhostData,     &
                      UpdateComm1DGhostDataWENO, UpdateCommShellGhostData
    use petscvariables, only: upetsc, phipetsc, ulocpetsc, philocpetsc, &
                              upetscWENO,  ulocpetscWENO,  &
                              upetscWENO_Shell, ulocpetscWENO_Shell
    implicit none

    ! loop indices
    integer :: ielem, jdir, idir
    integer :: inode, jnode, knode, gnode
    integer :: kelem
    integer :: iface
    integer :: i
    integer :: nshell
!  HACK TESTING new Viscous Path
!   integer :: jdir_face
!   real(wp) :: dx,t1
!  HACK TESTING new Viscous Path
    real(wp) :: l10_ldg_flip_flop


    ! low and high volumetric element indices
    integer :: iell, ielh
    ! normal direction at face
    real(wp) :: nx(3)
    ! LDC/LDG coefficient
    ! temporary arrays for phi and delta phi
    real(wp), allocatable :: dphi(:), wtmp(:)
    real(wp), allocatable :: phitmp(:,:), phig_err1(:,:), phig_err2(:,:)
!  HACK Testing
!   real(wp), allocatable :: phig_tmp3(:,:,:), phig_test(:,:,:)
!  HACK Testing

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    call UpdateComm1DGhostData(ug, ughst, upetsc, ulocpetsc, &
      nequations, nodesperelem, ihelems, nghost)

    ! Update ghost value in the interior plane for WENO p = 3
    if((discretization == 'SSWENO')) then

      call UpdateComm1DGhostDataWENO(ug, ughstWENO, upetscWENO, ulocpetscWENO, &
        nequations, nodesperelem, ihelems, nghost)

      if((WENO_type == 'Neighbr_WENO')) call Planar_Interpolation_WENO()

      nshell = nodesperface*nfacesperelem

      call UpdateCommShellGhostData(ugWENO_partner, ughstWENO_partner, upetscWENO_Shell, &
                        ulocpetscWENO_Shell, nequations, nshell, ihelems, nghost)

    endif

    ! loop over elements
    do ielem = iell, ielh
      ! loop over every node in element
      do inode = 1,nodesperelem
        ! compute primitive variables
        call conserved_to_primitive( ug(:,inode,ielem), &
          vg(:,inode,ielem), &
          nequations ) ! (navierstokes)
        ! compute entropy variables
        call primitive_to_entropy( vg(:,inode,ielem), &
          wg(:,inode,ielem), &
          nequations ) ! (navierstokes)
      end do
    end do

    if (viscous) then
      ! phitmp is calculated in computational space
      allocate(phig_err1(nequations,ndim),phig_err2(nequations,ndim))
!  HACK Testing
!     allocate(phig_tmp3(nequations,ndim,nodesperelem))
!     allocate(phig_test(nequations,ndim,nodesperelem))
!  HACK Testing
      allocate(phitmp(nequations,ndim))
      ! dphi is calculated at faces
      allocate(dphi(nequations))
      allocate(wtmp(nequations))
      ! LDC/LDG coefficient according to CNG-revisited
      ! loop over elements
      do ielem = iell, ielh
        ! compute computational gradients of the entropy variables
        !
        ! initialize phi
        phig(:,:,:,ielem) = 0.0_wp
        phig_err(:,:,:,ielem) = 0.0_wp
!  HACK TESTING new viscous path
!       phig_test(:,:,:) = 0.0_wp
!       phig_tmp3(:,:,:) = 0.0_wp
!  HACK Testing New viscous path
        grad_w_jacobian(:,:,:,ielem) = 0.0_wp
        ! loop over every node in element
        do inode = 1, nodesperelem
          ! reinitialize computational gradient to zero
          phitmp(:,:) = 0.0_wp
          phig_err1(:,:) = 0.0_wp
          phig_err2(:,:) = 0.0_wp
          ! loop over number of dependent elements in gradient
          do i = iagrad(inode), iagrad(inode+1)-1
            ! loop over dimensions
            do jdir = 1,ndim
              ! column/node from gradient operator in CSR format in
              ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
              jnode = jagrad(jdir,i)
              ! update gradient using coefficient and entropy variables at appropriate node
              phitmp(:,jdir)    = phitmp(:,jdir)    + dagrad(jdir,i)*wg(:,jnode,ielem) 
              phig_err1(:,jdir) = phig_err1(:,jdir) + dagrad(jdir,i)*vg(:,jnode,ielem) 
            end do
          end do

!  HACK Testing
!         phig_tmp3(:,:,inode) = phitmp(:,:) 
!  HACK Testing

          ! Store gradient of the entropy variables in computational space
          do i = 1, ndim
            grad_w_jacobian(:,i,inode,ielem) = phitmp(:,i)
          enddo
          ! transform to physical space using dxi_jdir/dx_idir
          do jdir = 1,ndim
            do idir = 1,ndim
              phig(:,idir,inode,ielem) = phig(:,idir,inode,ielem) + phitmp(:,jdir)*r_x(jdir,idir,inode,ielem)
              phig_err2(:,idir) = phig_err2(:,idir) + phig_err1(:,jdir)*r_x(jdir,idir,inode,ielem)
            end do
          end do

          phig_err(:,:,inode,ielem) = matmul(dWdV(vg(:,inode,ielem),nequations),phig_err2(:,:))

        end do
         
        ! LDC/LDG penalty on phig
        ! 
        ! loop over only faces
        do iface = 1, nfacesperelem
!  HACK Testing
!         ! sign so normal is facing outward
!         dx = sign(1.0_wp,real(facenormalcoordinate(iface),wp))
!         jdir_face = abs(facenormalcoordinate(iface))
!  HACK Testing
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

              call conserved_to_primitive(ughst(:,gnode),dphi,nequations)
              call primitive_to_entropy(dphi,wtmp,nequations)

              ! LDC/LDG penalty value
              l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
              dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg(:,inode,ielem) - wtmp)
              ! add LDC/LDG penalty to each physical gradient using the normal
              do jdir = 1,ndim
                phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
              end do
!  HACK Testing
!             phig_tmp3(:,jdir_face,inode) = phig_tmp3(:,jdir_face,inode) + dphi(:)*dx
!  HACK Testing
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

              ! LDC/LDG penalty value
              l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
              dphi(:) = l10_ldg_flip_flop*pinv(1)*(wg(:,inode,ielem) - wg(:,knode,kelem))
              ! add LDC/LDG penalty to each physical gradient using the normal
              do jdir = 1,ndim
                phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
              end do
!  HACK Testing
!             phig_tmp3(:,jdir_face,inode) = phig_tmp3(:,jdir_face,inode) + dphi(:)*dx
!  HACK Testing
            end do
          end if
        end do
!  HACK Testing
!       phig_test(:,:,:) = 0.0_wp
!       do inode = 1, nodesperelem
!         do jdir = 1,ndim
!           do idir = 1,ndim
!             phig_test(:,idir,inode) = phig_test(:,idir,inode) + phig_tmp3(:,jdir,inode)*r_x(jdir,idir,inode,ielem)
!           end do
!         end do
!       end do
!       t1 = maxval(abs(phig_test(:,:,:) - phig(:,:,:,ielem)))
!       if(t1 >= 1.0e-15) write(*,*)'error in reconcilestates',t1
!       This line runs the new path through the code
!       phig(:,:,:,ielem) = phig_test(:,:,:)
!  HACK Testing
      end do
      deallocate(phitmp,dphi,wtmp)
      deallocate(phig_err1,phig_err2)
!  HACK Testing
!     deallocate(phig_tmp3,phig_test)
!  HACK Testing
    end if

!   write(*,*)'error in phig',maxval(abs(phig-phig_err))
!   stop

    call UpdateComm2DGhostData(phig, phighst, phipetsc, philocpetsc, &
      nequations, 3, nodesperelem, ihelems, nghost)

    return
  end subroutine nse_reconcilestates

  
  !============================================================================

  subroutine nse_calc_dudt_LSRK(tin)
    ! This subroutine calculates the time derivative of the
    ! conserved variables at every node.
    use variables
    use referencevariables
    use controlvariables, only: verbose
    use nsereferencevariables, only: entropy_viscosity
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode,ielem, jdir

    ! low and high volumetric element indices
    integer :: iell, ielh
    real(wp) :: t1

    ! update the primitive and entropy variables and
    ! the LDG/LDC gradients
    !call nse_reconcilestates() IMPORTANT NOTE: This is now called in timeinteg.f90 

    ! low : high volumetric element index
    iell = ihelems(1) ;  ielh = ihelems(2)

    ! loop over all elements
    do ielem = iell, ielh
          
      !  Calculate the elementwise Divergence  \/ * (F - Fv)

      call Flux_Divergence(tin,ielem)           !  result in divF

      !  Form the elementwise SAT_Penalties

      call SAT_Penalty(tin,ielem)        !  result in gsat
      !  
      ! compute time derivative
      ! 
      ! loop over all nodes in the element
      do inode = 1, nodesperelem
        ! reset the time derivative to zero
        dudt(:,inode,ielem) = 0.0_wp
        ! add the contribution from the flux divergence in each direction

        do jdir = 1,ndim
          dudt(:,inode,ielem) = dudt(:,inode,ielem) - divf(:,jdir,inode,ielem)
        end do
        ! add the contribution from the boundary and interface penalties
        dudt(:,inode,ielem) = (dudt(:,inode,ielem) + gsat(:,inode,ielem))/Jx_r(inode,ielem) ! Thus this is the dudt of u and NOT u*J

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
    use variables, only: gsat, Jx_r, boundaryelems
    use collocationvariables, only: pvol
    use controlvariables, only: verbose
    use referencevariables
    use mpimod
    implicit none

    ! indices
    integer :: inode,ielem

    ! low and high volumetric element indices
    integer :: iell, ielh

    ! different error estimates
    real(wp) :: l2(2), linf
    real(wp) :: l2sum(2), linfmax
    ! local error
    real(wp), allocatable :: ex(:)
    integer :: ierr

    ! update primitives and entropy variables
    call nse_reconcilestates()

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
    use variables, only: ug,uhat,vg,Jx_r,boundaryelems
    use collocationvariables, only: pvol
    use controlvariables, only: verbose
    use referencevariables
    use mpimod
    implicit none

    ! indices
    integer :: inode,ielem

    ! low and high volumetric element indices
    integer :: iell, ielh

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

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! initialize errors to zero
    l2 = 0.0_wp
    linf = 0.0_wp
    sglob = 0.0_wp
    ! loop over elements
    do ielem = iell, ielh
      ! loop over each index in the element
      do inode = 1, nodesperelem
        ! compute the local embedded error
        ex = ug(:,inode,ielem) &
          - uhat(:,inode,ielem)
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
    use variables, only: uhat, ug, vg, xg, Jx_r, boundaryelems, phig, mut
    use collocationvariables, only: pvol
    use referencevariables
    use mpimod
    implicit none

    ! Time of evaluation
    real(wp), intent(in) :: tin
    ! Indices
    integer :: inode,ielem, i

    ! Low and high volumetric element indices
    integer :: iell, ielh

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
    character(len=23) :: file_name_err

    ! Assign global error norms file name
!    file_name_err = 'global_error_norms.data'

    ! Arbitrary normal vector
    ctmp = 0.0_wp

    ! Update primitive variables
    call nse_reconcilestates()

    allocate(uex(nequations))
    allocate(vex(nequations))
    allocate(ftmp(nequations))

    ! Low volumetric element index
    iell = ihelems(1)
    ! High volumetric element index
    ielh = ihelems(2)

    ! Initialize error to zero
    l2 = 0.0_wp; l2sum = 0.0_wp
    linf = 0.0_wp
    ! Loop over all elements
    do ielem = iell, ielh
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

    if (myprocid == 0) then
!      open(unit=40,file=file_name_err,Access='append',Status='unknown')
!      i = nodesperedge-1
      write(*,*) 'P',nodesperedge-1
      write(*,*) 'L1   error: ', l1sum(1)
      write(*,*) 'L2   error: ', l2sum(1)
      write(*,*) 'Linf error: ', linfmax
      write(40,230)i,l1sum(1),l2sum(1),linfmax
      230  format('P',I2,1x,'  l1 error:  ', e17.10,1x,  &
                            '  l2 error:  ', e17.10,1x,  &
                            'linf error:  ', e17.10,1x)
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
      integer   :: it,j,i

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
    use controlvariables, only : variable_viscosity
  
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

    integer :: i

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

    real(wp) :: x0,y0
    real(wp) :: f
    real(wp) :: alpha, rin2
   

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
          if(tmp.eq.zero)tmp = 1.0d-14
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

  subroutine Flux_Divergence(tin,ielem)
    ! This subroutine calculates elementwise 
    ! the Divergence of the Conservative Flux
    use variables
    use referencevariables
    use nsereferencevariables
    use controlvariables, only: verbose, discretization
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode, jdir
    integer :: jnode
    integer :: i

    real(wp), dimension(nequations) :: t1
    ! normal vector
!   real(wp) :: nx(3)
     
    call Flux_Div_Pencil(ielem)    !  result in divF

    !divf = 0.0_wp

    if (viscous) then
      ! loop over all nodes in element
      do inode = 1,nodesperelem

        ! calculate viscous flux
        fvg (:,:,inode,ielem) = 0.0_wp

        fvg(:,1:ndim,inode,ielem) = Jx_r(inode,ielem) * &
            viscousflux3D( vg(:,inode,ielem), &
            phig(:,:,inode,ielem), &
            r_x(:,:,inode,ielem), &
            nequations, &
            ndim, &
            mut(inode,ielem))
      end do

      !
      ! calculate divergence of the flux
      ! 

      ! loop over all nodes in the element
      do inode = 1,nodesperelem

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
    use controlvariables, only: verbose, discretization
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem
    real(wp), intent(in) :: tin

    ! indices
    integer :: idir,  jdir
    integer :: inode, jnode
    integer :: i
    real(wp), dimension(3)  :: n_i,n_j
    real(wp), dimension(:,:),   allocatable :: fvg_S
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


  !============================================================================
  ! sat_penalty - Calculates both inviscid and viscous penalty according to the
  ! SAT procedure.
  
  subroutine SAT_Penalty(tin,ielem)
    
    ! Load modules
    use variables
    use referencevariables
    use nsereferencevariables
    use controlvariables, only: verbose, heat_entropy_flow_wall_bc, Riemann_Diss,&
                             & Riemann_Diss_BC
    use collocationvariables, only: iagrad, jagrad, dagrad, gradmat, pinv, &
                                  & pvol, l01, l00, ldg_flip_flop_sign, &
                                  & alpha_ldg_flip_flop, Sfix
    use initgrid

    ! Nothing is implicitly defined
    implicit none

    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode, jnode, knode, gnode
    integer :: kelem
    integer :: iface
    integer :: i,j,k

    ! reconstructed flux
    real(wp), allocatable :: fstar(:), fstarV(:)
    real(wp), allocatable :: fRoeI(:), fLLF(:)
    ! local normal flux
    real(wp), allocatable :: fn(:), fnV(:)
    ! boundary conservative and primitive states
    real(wp), allocatable :: ustar(:), vstar(:), wstar(:), phistar(:,:)
    ! right and left eigenvector matrices
    real(wp), allocatable, dimension(:,:) :: smat, sinv
    ! eigenvalues
    real(wp), allocatable, dimension(:)   :: ev, evabs
    ! average state
    real(wp), allocatable, dimension(:)   :: Vav
    ! normal vector
    real(wp) :: nx(3)
    ! Lax-Freidrich max Eigenvalue
    real(wp) :: evmax
    ! penalty parameter for interfaces
    real(wp) :: t1 

    logical,  parameter :: testing         = .false.

    real(wp), parameter :: mirror          = -1.0_wp
    real(wp), parameter :: no_slip         = -0.0_wp

!  Original
!   real(wp), parameter :: Cevmax          =  4.0_wp
!  Nail's recent modification
    real(wp), parameter :: Cevmax          =  1.0_wp
    real(wp), parameter :: deltaU          =  0.1_wp
    real(wp), parameter :: stab_strength   =  1.1_wp
    real(wp), parameter :: bc_pen_strength =  1.0_wp
    real(wp), parameter :: LocalLaxF_factor=  2.0_wp

    real(wp), dimension(nequations,nequations) :: hatc_side_1, hatc_side_2, &
                                                & matrix_ip

    real(wp), dimension(nequations)    :: tmpr
    real(wp), dimension(nequations)    :: w_side_1, w_side_2
!   Original Version
!   real(wp), dimension(nequations)    :: UavAbs, switch
!   Nails Version
    real(wp) :: UavAbs, switch


    real(wp), dimension(nequations)    ::   vg_On,   vg_Off
    real(wp), dimension(nequations)    ::   ug_On,   ug_Off
    real(wp), dimension(nequations,3)  :: phig_On, phig_Off
    real(wp), dimension(nequations)    :: SAT_Pen
    real(wp), dimension(nequations)    :: f_viscous_normal_ghost

    real(wp), dimension(nequations)    :: prim_ghost_adiabatic
    real(wp), dimension(nequations)    :: entr_ghost_adiabatic
    real(wp), dimension(nequations)    :: prim_ref
    
    real(wp), dimension(nequations,3)  :: grad_prim_int, grad_prim_ghost, grad_entr_ghost
    real(wp), dimension(nequations,3)  :: grad_prim_int_normal,grad_prim_int_tangent
    real(wp), dimension(nequations,3)  :: grad_entr_int_normal,grad_entr_int_tangent

    real(wp), dimension(3)             :: unit_normal, normal_vel, tangent_vel, ref_vel_vector
    
    real(wp) :: l01_ldg_flip_flop


    ! allocate local arrays
    allocate(ustar(nequations))
    allocate(vstar(nequations))
    allocate(wstar(nequations))
    allocate(fRoeI(nequations))
    allocate(fLLF(nequations))
    allocate(fstar(nequations))
    allocate(fstarV(nequations))
    allocate(fn(nequations))
    allocate(fnV(nequations))
    allocate(phistar(nequations,3))

    allocate(ev(nequations))
    allocate(evabs(nequations))
    allocate(vav(nequations))
    allocate(smat(nequations,nequations))
    allocate(sinv(nequations,nequations))

    ! initialize penalty
    gsat(:,:,ielem) = 0.0_wp
       
       
    ! Loop over each face
    do iface = 1,nfacesperelem

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       Boundary Contribution to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      if (ef2e(1,iface,ielem) < 0) then

        if (abs(ef2e(1,iface,ielem)) /= 6) then
 
          ! Specify the Boundary Condition procedure on the face
          call set_boundary_conditions(ef2e(1,iface,ielem),InitialCondition)

          ! Loop over each node on the face
          do i = 1,nodesperface

            ! Volumetric node index corresponding to face and node on face indices
            inode = kfacenodes(i,iface)
          
            ! Facial index corresponding to face and node on face indices
            jnode = (iface-1)*nodesperface + i
          
            ! Outward facing normal of facial node
            nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)

            ! Compute the boundary state
            call conserved_to_primitive( ug(:,inode,ielem), vg(:,inode,ielem), nequations ) ! (navierstokes)

            vstar(:) = vg(:,inode,ielem)
            phistar = phig(:,:,inode,ielem)
            
            call BoundaryCondition(vstar,phistar,fnV,nx,xg(:,inode,ielem),tin, &
              & nequations,ndim,mut(inode,ielem))

            call primitive_to_conserved( vstar, ustar, nequations) 
            call primitive_to_entropy(vstar, wstar, nequations) 

            ! ==  Eigen values/vectors
            call roeavg( vg(:,inode,ielem), vstar, Vav, nequations )         

            call CharacteristicDecomp( vav, nequations, sinv, smat, ev, nx ) 

            evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)
            ! ==  Eigen values/vectors

            ! ==  Fluxes
            fn = normalflux( vg(:,inode,ielem), nx, nequations )                     ! (Euler Flux)

            select case(Riemann_Diss_BC)
              case('LocalLaxF')
                fLLF  = half * ( normalflux( vg(:,inode,ielem), nx, nequations )  &
                    &        +   normalflux( vstar            , nx, nequations )  &
                             +   LocalLaxF_factor*evmax*(ug(:,inode,ielem)-ustar) )
                fstar = fLLF
              case('Roe')
                fstar = EntropyConsistentFlux(vg(:,inode,ielem), vstar, nx, nequations ) ! (Entropy Flux)
 !              fstar = Entropy_KE_Consistent_Flux(vg(:,inode,ielem), vstar, nx, nequations ) ! (Entropy Flux)
                fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-wstar(:)) )
            end select

            fstarV = normalviscousflux( vg(:,inode,ielem), phig(:,:,inode,ielem), nx, nequations,mut(inode,ielem)) &
              & - fnV(:)

            ! ==  Fluxes

            ! Compute the IP penalty contribution, i.e. M (u-v), where M, in the
            ! case without flip-flop is a positive matrix defined as:
            ! M = pinv(1) (c_ii_side_1 + c_ii_side_2)/2, in the normal direction
            ! ------------------------------------------------------------------

            ! c_ii_side_1 matrix ! c_ii_side_2 matrix
            hatc_side_1 = matrix_hatc_node(vg(:,inode,ielem),nx,nx,nequations)
            hatc_side_2 = matrix_hatc_node(vstar,nx,nx,nequations)

            ! IP penalty matrix
            matrix_ip = 0.5_wp*(hatc_side_1 + hatc_side_2)*pinv(1)

            call primitive_to_entropy(vg(:,inode,ielem),w_side_1,nequations)
              
            call primitive_to_entropy(vstar,w_side_2,nequations)

            ! Add the LDG and IP terms to the penalty
            !l01_ldg_flip_flop = l01*(1.0_wp - ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)

            gsat(:,inode,ielem) = gsat(:,inode,ielem) + &
              & pinv(1)* ( (fn - fstar) + l01*fstarV*bc_pen_strength ) - &
              & pinv(1)* l00*matmul(matrix_ip,w_side_1 - wstar)

          end do
        
        else         ! Entropy Stable Solid Wall BC

          ! Loop over each node on the face
          do i = 1,nodesperface

            ! Volumetric node index corresponding to face and node on face indices
            inode = kfacenodes(i,iface)
          
            ! Facial index corresponding to face and node on face indices
            jnode = (iface-1)*nodesperface + i
          
            ! Outward facing normal of facial node
            nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)

            ! Unit normal direction
            unit_normal = nx/Magnitude(nx)

            ! Normal velocity
            normal_vel = dot_product(vg(2:4,inode,ielem),unit_normal)*unit_normal

            ! Tangent velocity
            tangent_vel = vg(2:4,inode,ielem) - normal_vel

            vstar(1) = vg(1,inode,ielem)
            vstar(2) = tangent_vel(1) - normal_vel(1)
            vstar(3) = tangent_vel(2) - normal_vel(2)
            vstar(4) = tangent_vel(3) - normal_vel(3)
!  HACK HUGE  (need an entropy proof for this change
!           vstar(2) = -vg(2,inode,ielem)
!           vstar(3) = -vg(3,inode,ielem)
!           vstar(4) = -vg(4,inode,ielem)
!  HACK HUGE  (need an entropy proof for this change
            vstar(5) = vg(5,inode,ielem) 

            ! Compute the entropy variables in the ghost node for imposing the
            ! nonpenetration condition
            call primitive_to_conserved( vstar, ustar, nequations) 
            call primitive_to_entropy(vstar,wstar,nequations) 

            ! Compute the roe average state of the primitive variables
            call roeavg(vg(:,inode,ielem),vstar,vav,nequations)

            ! Compute characteristic decomposition
            call CharacteristicDecomp(vav,nequations,sinv,smat,ev,nx)

            evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)
            ! ==  Eigen values/vectors

            ! ==  Fluxes
            fn = normalflux( vg(:,inode,ielem), nx, nequations )                     ! (Euler Flux)

!           select case(Riemann_Diss_BC)
!             case('LocalLaxF')
                fLLF  = half * ( normalflux( vg(:,inode,ielem), nx, nequations )  &
                    &        +   normalflux( vstar            , nx, nequations )  &
                             +   LocalLaxF_factor*evmax*(ug(:,inode,ielem)-ustar) )
                fstar = fLLF
!             case('Roe')
!               fstar = EntropyConsistentFlux(vg(:,inode,ielem), vstar, nx, nequations ) ! (Entropy Flux)
!!              fstar = Entropy_KE_Consistent_Flux(vg(:,inode,ielem), vstar, nx, nequations ) ! (Entropy Flux)
!               fstar = fstar + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-wstar(:)) )
!           end select

            ! Compute the boundary data in primitive variables for imposing the no-slip isothermal wall BC
            prim_ghost_adiabatic(1) = vg(1,inode,ielem)
            prim_ghost_adiabatic(2) = 0.0_wp
            prim_ghost_adiabatic(3) = 0.0_wp
            prim_ghost_adiabatic(4) = 0.0_wp
            prim_ghost_adiabatic(5) = vg(5,inode,ielem)

            ! Compute the entropy variables in the ghost node for imposing the no-slip isothermal wall BC
            call primitive_to_entropy(prim_ghost_adiabatic,entr_ghost_adiabatic,nequations)

            fstarV = normalviscousflux( vg(:,inode,ielem), phig(:,:,inode,ielem), nx, nequations,mut(inode,ielem)) 

            ! Grad(V) = dVdW Grad(W)  = dvdw phi
            grad_prim_int = MatMul(dVdW(vg(:,inode,ielem),nequations),phig(:,:,inode,ielem))

            ! Set gradient of the primitive variables in the ghost node
            ! =========================================================
            grad_entr_ghost(1,:) = phig(1,:,inode,ielem) !grad_prim_int(1,:) 
            grad_entr_ghost(2,:) = phig(2,:,inode,ielem) !grad_prim_int(2,:) 
            grad_entr_ghost(3,:) = phig(3,:,inode,ielem) !grad_prim_int(3,:) 
            grad_entr_ghost(4,:) = phig(4,:,inode,ielem) !grad_prim_int(4,:) 

            ! Normal component of grad_prim_int(V)
            do j = 1, nequations
              grad_entr_int_normal(j,:) = dot_product(phig(j,:,inode,ielem),unit_normal(:))*unit_normal(:) !dot_product(grad_prim_int(j,:),unit_normal(:))*unit_normal(:)
            end do

            grad_entr_int_tangent(:,:) = phig(:,:,inode,ielem) - grad_entr_int_normal!grad_prim_int(:,:) - grad_prim_int_normal(:,:)
  
            grad_entr_ghost(5,:) = grad_entr_int_tangent(5,:) + heat_entropy_flow_wall_bc*unit_normal(:)/vg(5,inode,ielem)

            !grad_prim_ghost(5,:) = 0.0_wp !grad_prim_int(5,:) - grad_prim_int_normal(5,:)

            ! Compute normal viscous flux arising from the ghost point
            !f_viscous_normal_ghost = normalviscousflux(vg(:,inode,ielem),MatMul(dWdV(vg(:,inode,ielem),nequations,mut(inode,ielem)), &
            !  & grad_prim_ghost),nx,nequations)
            f_viscous_normal_ghost = normalviscousflux(vg(:,inode,ielem),grad_entr_ghost,nx,nequations,mut(inode,ielem))      
            
            ! c_ii_side_1 matrix
            hatc_side_1 = matrix_hatc_node(vg(:,inode,ielem),nx,nx,nequations)

            ! c_ii_side_2 matrix
            ref_vel_vector(:) = U0
            prim_ref(1) = rho0
            prim_ref(2) = Magnitude(ref_vel_vector)
            prim_ref(3) = Magnitude(ref_vel_vector)
            prim_ref(4) = Magnitude(ref_vel_vector)
            prim_ref(5) = T0
            hatc_side_2 = matrix_hatc_node(prim_ref,nx,nx,nequations)

            hatc_side_2(5,:) = 0.0_wp
            hatc_side_2(:,5) = 0.0_wp

            ! IP penalty matrix
            matrix_ip = hatc_side_2!*pinv(1)

            call primitive_to_entropy(vg(:,inode,ielem),w_side_1,nequations)
              
            ! Compute the penalty term
            gsat(:,inode,ielem) = gsat(:,inode,ielem) &
              & + pinv(1)*(fn - fstar) &
              & - pinv(1)*(fstarV-f_viscous_normal_ghost) &
              & - pinv(1)*1.0_wp*matmul(matrix_ip,w_side_1-entr_ghost_adiabatic)

          end do

        end if

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       Off Processor Contributions to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      else if (ef2e(3,iface,ielem) /= myprocid) then
        
        ! This is a parallel interface
        do i = 1, nodesperface
          jnode = nodesperface*(iface-1)+i
          
          ! Volumetric node index corresponding to facial node index
          inode = ifacenodes(jnode)
          
          ! Index in ghost
          gnode = efn2efn(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack
          
          ! Element index of partner node
          kelem = efn2efn(2,jnode,ielem)
          
          ! Outward facing normal of facial node
          nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)

          ug_On(:)  = ug(:,inode,ielem)
          ug_Off(:) = ughst(:,gnode)
          call conserved_to_primitive(ug_On (:), vg_On (:), nequations)
          call conserved_to_primitive(ug_Off(:), vg_Off(:), nequations)

          phig_On (:,:) = phig(:,:,inode,ielem)
          phig_Off(:,:) = phighst(:,:,gnode)

          SAT_Pen(:) =  SAT_Inv_Vis_Flux( nequations,iface,ielem,    &
                                        & vg_On,vg_Off,              &
                                        & phig_On,phig_Off,          &
                                        & nx,Jx_r(inode,ielem),      &
                                        & pinv(1), mut(inode,ielem))

          gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * SAT_Pen(:)

        end do
      
      else

        do i = 1,nodesperface
        
          ! Index in facial ordering
          jnode = nodesperface*(iface-1)+i
            
          ! Volumetric node index corresponding to facial node index
          inode = ifacenodes(jnode)
            
          ! Volumetric index of partner node
          knode = efn2efn(1,jnode,ielem)
            
          ! Element index of partner node
          kelem = efn2efn(2,jnode,ielem)
          
          ! Outward facing normal of facial node
          nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)
          
          vg_On(:)  = vg(:,inode,ielem)
          vg_Off(:) = vg(:,knode,kelem)

          phig_On (:,:) = phig(:,:,inode,ielem)
          phig_Off(:,:) = phig(:,:,knode,kelem)

          SAT_Pen(:) =  SAT_Inv_Vis_Flux( nequations,iface,ielem,    &
                                        & vg_On,vg_Off,              &
                                        & phig_On,phig_Off,          &
                                        & nx,Jx_r(inode,ielem),      &
                                        & pinv(1), mut(inode,ielem))

          gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv(1) * SAT_Pen(:)

        end do
      end if
    end do
    
    
    ! Deallocate memory
    deallocate(ustar)
    deallocate(vstar)
    deallocate(wstar)
    deallocate(fRoeI)
    deallocate(fLLF)
    deallocate(fstar)
    deallocate(fstarV)
    deallocate(phistar)
    deallocate(fn)
    deallocate(fnV)
    deallocate(sinv)
    deallocate(smat)
    deallocate(ev)
    deallocate(evabs)
    deallocate(vav)

    return
  end subroutine SAT_Penalty

  !============================================================================
  
  subroutine nse_calcrhsexplicit(irk,tin)
    ! This subroutine calculates the time derivative of the
    ! conserved variables at every node.
    use variables
    use referencevariables
    use controlvariables
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol
    use mpimod
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer,  intent(in) :: irk
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode,ielem, jdir

    ! low and high volumetric element indices
    integer :: iell, ielh

!   integer :: i_err

    ! update the primitive and entropy variables and
    ! the LDG/LDC gradients
    call nse_reconcilestates() ! (navierstokes)

    ! low : high volumetric element index
    iell = ihelems(1) ;  ielh = ihelems(2)

    ! loop over all elements
    do ielem = iell, ielh
      !                                        __
      !  Calculate the elementwise Divergence  \/ * (F - Fv)
       
      if(IMEX_element == 'explicit') call Flux_Divergence(tin,ielem)    !  result in divF

      !  Form the elementwise SAT_Penalties

      if(IMEX_penalty == 'explicit') call SAT_Penalty(tin,ielem)        !  result in gsat
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


  subroutine nse_calcrhsimplicit(irk,tin)
    ! This subroutine calculates the time derivative of the
    ! conserved variables at every node.
    use variables
    use referencevariables
    use controlvariables
    use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer,  intent(in) :: irk
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode,ielem, jdir

    ! low and high volumetric element indices
    integer :: iell, ielh

    ! update the primitive and entropy variables and
    ! the LDG/LDC gradients
    call nse_reconcilestates() ! (navierstokes)

    ! low : high volumetric element index
    iell = ihelems(1) ;  ielh = ihelems(2)

    ! loop over all elements
    do ielem = iell, ielh
      !                                        __
      !  Calculate the elementwise Divergence  \/ * (F - Fv)
       
      if(IMEX_element == 'implicit') call Flux_Divergence(tin,ielem)    !  result in divF

      !  Form the elementwise SAT_Penalties

      if(IMEX_penalty == 'implicit') call SAT_Penalty(tin,ielem)        !  result in gsat
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
    
    real(wp)                 :: rho, p

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
    use collocationvariables, only: pmat, pvol
    use mpimod 

    implicit none

    real(wp), intent(inout) :: dt_global
    real(wp), dimension(:), allocatable :: dt_min_proc

    integer :: elem_low, elem_high
    integer :: s_tag, r_tag, m_size, &
               s_request_dt_min, r_request_dt_min, i_err

    integer :: s_status(mpi_status_size)
    integer :: r_status(mpi_status_size)

    integer :: i_elem, i_node
    integer :: i, j, k, m

    real(wp)               :: Lngth, a0, dt0, dt_min, dt_global_max
    real(wp)               :: tI, tV
    real(wp), dimension(3) :: xvec, sq, ucon

    continue

    ! Low and high  volumetric element index
    elem_low  = ihelems(1) ;  elem_high = ihelems(2)

    ! Compute gradient of the velocity components
    ! ===========================================
    dt_global_max = dt_global*1.1_wp
    dt_min = 100.0_wp

    do i_elem = elem_low, elem_high
      do i_node = 1, nodesperelem
                
          Lngth = pvol(i_node)**(third)

          a0  = sqrt(abs(gamma0*vg(5,i_node,i_elem)/gM2))

          sq(1) = magnitude(r_x(1,:,i_node,i_elem))
          sq(2) = magnitude(r_x(2,:,i_node,i_elem))
          sq(3) = magnitude(r_x(3,:,i_node,i_elem))

          ucon(1) = dot_product(r_x(1,:,i_node,i_elem),vg(2:4,i_node,i_elem))
          ucon(2) = dot_product(r_x(2,:,i_node,i_elem),vg(2:4,i_node,i_elem))
          ucon(3) = dot_product(r_x(3,:,i_node,i_elem),vg(2:4,i_node,i_elem))

          tI  =  sum(abs(ucon(:))) + a0 * ( sum(sq) )
          tV  =  Re0Inv * magnitude(sq)

          dt0 = CFL / (tI / Lngth + tV / Lngth / Lngth)

          dt_min = min(dt_min,dt0)

      enddo
    enddo

    ! Reduce values on all processes to a single value
    if(myprocid == 0 )  then
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

    ! Create a barrier synchronization in the group. 
    call mpi_barrier(petsc_comm_world,i_err)

    ! Broadcast dt_global
    m_size = 1
    call mpi_bcast(dt_global,m_size,mpi_default_wp,0,petsc_comm_world,i_err)

    return
  end subroutine compute_explicit_timestep

  !============================================================================
  
  !============================================================================
  ! compute_vorticity_field_elements - Computes the vorticity field for all the
  ! elements

  subroutine compute_vorticity_field_elements()
    ! This subroutine computes the vorticity field given the velocity field.
    ! This means that the primitive variable must be already computed.

    ! Load modules
    use variables
    use referencevariables
    use collocationvariables, only: iagrad,jagrad,dagrad

    ! Nothing is implicitly defined
    implicit none

    real(wp), dimension(5,3) :: GradV

    integer :: elem_low, elem_high
    
    integer :: i_elem, i_node

    continue
    
    ! Low volumetric element index
    elem_low  = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Compute gradient of the velocity components
    ! ===========================================
    do i_elem = elem_low, elem_high
      do i_node = 1, nodesperelem

        !  __           __
        !  \/ V  = dVdW \/ W  = dVdW phi
        !  
        GradV(:,:)    = MatMul(dVdW(vg(:,i_node,i_elem),nequations),phig(:,:,i_node,i_elem))

        ! Compute the vorticity
        omega(:,i_node,i_elem) = 0.0_wp 

        ! Note:  GradV(1,:) is gradient of density
        omega(1,i_node,i_elem) = GradV(4,2) - GradV(3,3)
        omega(2,i_node,i_elem) = GradV(2,3) - GradV(4,1)
        omega(3,i_node,i_elem) = GradV(3,1) - GradV(2,2)
      
      end do
    end do

    return
  end subroutine compute_vorticity_field_elements

  !============================================================================
  
  !============================================================================
  ! kinetic_energy_element - Computes the kinetic energy for one element

  pure function kinetic_energy_element(elem_id)
    ! This subroutine computes the kinetic energy in one element given the 
    ! velocity field. This means that the primitive variable must be already computed.

    ! Load modules
    use variables
    use referencevariables
    use collocationvariables, only: iagrad,jagrad,dagrad

    ! Nothing is implicitly defined
    implicit none

    real(wp), dimension(nodesperelem) :: kinetic_energy_element

    integer, intent(in) :: elem_id
    
    integer :: i_node

    continue 

    ! Compute kinetic energy
    ! ======================
    do i_node = 1, nodesperelem

      kinetic_energy_element(i_node) = 0.0_wp 

      kinetic_energy_element(i_node) = 0.5_wp*dot_product(vg(2:4,i_node,elem_id),vg(2:4,i_node,elem_id))
      
    end do

    return
  end function kinetic_energy_element

  !============================================================================

  !============================================================================
  ! compute_primitive_variables_elements - Computes the primitive variables for
  !                                        all elements using the conservative 
  !                                        variables. 

  subroutine compute_primitive_variables_elements()

    ! Load modules
    use variables
    use referencevariables

    ! Nothing is defined implicitly
    implicit none

    integer :: elem_low, elem_high
    integer :: i_elem, i_node

    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Loop over elements
    do i_elem = elem_low, elem_high
      ! Loop over nodes in each element
      do i_node = 1, nodesperelem

        ! Calculate primitive variables from conservative variables
        call conserved_to_primitive( &
          uin = ug(:,i_node,i_elem), &
          vout = vg(:,i_node,i_elem), &
          nq = nequations)
      
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

    ! Nothing is defined implicitly
    implicit none

    integer :: elem_low, elem_high
    integer :: i_elem, i_node

    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Loop over elements
    do i_elem = elem_low, elem_high
      ! Loop over nodes in each element
      do i_node = 1, nodesperelem
            
        ! Calculate entropy variables from primitive variables
        call primitive_to_entropy( &
        vin = vg(:,i_node,i_elem), &
        wout = wg(:,i_node,i_elem), &
        nq = nequations )

      enddo
    enddo

    return
  end subroutine compute_entropy_variables_elements

  !============================================================================

  !============================================================================
  ! compute_specific_entropy_elements - Computes the specific entropy using the 
  ! primitive variables 

  subroutine compute_specific_entropy_elements()
 
    ! Load modules
    use variables
    use referencevariables

    ! Nothing is defined implicitly
    implicit none

    integer :: elem_low, elem_high
    integer :: i_elem, i_node


    ! Low volumetric element index
    elem_low = ihelems(1)

    ! High volumetric element index
    elem_high = ihelems(2)

    ! Loop over elements
    do i_elem = elem_low, elem_high
      ! Loop over nodes in each element
      do i_node = 1, nodesperelem
            
        ! Calculate specific entropy
        specific_entropy(i_node,i_elem) = specificentropy(vg(:,i_node,i_elem),&
                                                      nequations)

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
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)

    ! indices
    ! high and low element indices for volume elements
    integer :: iell, ielh
    integer :: inode, jdir, ielem
    integer :: jnode
    integer :: i

    real(wp),  parameter :: al = 0.01_wp
    real(wp),  allocatable, dimension(:,:) :: tmp88,tmp89
    real(wp)             :: one_al

    one_al = 1.0_wp - al

    allocate(tmp88(1:nequations,1:nodesperelem))
    allocate(tmp89(1:nequations,1:nodesperelem))

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! loop over elements
    do ielem = iell, ielh
      ! loop over all nodes in the element
      tmp89(:,:) = ug(:,:,ielem)
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
        tmp89(:,:) = tmp88(:,:)
      end do
      ug(:,:,ielem) = one_al * ug(:,:,ielem)  + al*tmp88(:,:)
    end do
    deallocate(tmp88) ; deallocate(tmp89)

    return
  end subroutine Solution_Filter
 
  subroutine Flux_Div_Pencil(ielem)
    ! This subroutine calculates elementwise 
    ! the Divergence of the Conservative Flux
    use variables
    use referencevariables
    use controlvariables, only: discretization, Entropy_Correction
    use collocationvariables
    use SSWENOvariables
    use initgrid

    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem

    ! indices
    integer :: inode, jdir, iface, ipen
    integer :: i, k, l

!   real(wp), dimension(N_Soln_Pts) :: tmpS
!   real(wp), dimension(N_Flux_Pts) :: tmpF
!   real(wp), dimension(nequations) :: tmpE
!   real(wp), dimension(         3) :: tmpD

    real(wp), dimension(nequations,N_Soln_Pts) :: ugS, wgS, d_fgS
!   real(wp), dimension(         3,N_Soln_Pts) :: nxS
    real(wp), dimension(         3,N_Soln_Pts) :: JnS

    real(wp), dimension(nequations,N_Flux_Pts) :: ugF, wgF, vgF, fgF
!   real(wp), dimension(         3,N_Flux_Pts) :: nxF
    real(wp), dimension(         3,N_Flux_Pts) :: JnF
!   real(wp), dimension(           N_Flux_Pts) :: JxF

    integer,  dimension(2)                     :: faceLR
    real(wp), dimension(nequations)            :: uLL,uL,uR,uRR
    real(wp), dimension(nequations)            :: uLLT,uRRT
    real(wp)                                   :: dxL,dx,dxR
    real(wp)                                   :: t1, scl
    real(wp), parameter                        :: tol1 = 1.0e-10_wp

    integer                                    :: jnode,gnode
    integer                                    :: i_node,k_node,k_elem,k_face
    integer                                    :: inb


    select case(discretization)

    case('SpecColl')

      do jdir = 1,ndim            ! Directional loop

        do ipen = 1,nodesperface

          !  Grab a pencil of data
          do i = 1,N_Soln_Pts
            inode    = Pencil_Coord(N_Soln_Pts,jdir,ipen,i)
            ugS(:,i) =   ug(     :,inode,ielem)    !  Neqns
            JnS(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
          enddo

          !  Extrapolate data from Solution points to Flux points
          ugF(:,:) = zero  
          do k = 1,nequations
            ugF(k,1:N_Flux_Pts) = matmul(Ext_S2F,ugS(k,1:N_Soln_Pts))
          enddo
          JnF(:,:) = zero  
          do l = 1,ndim
            JnF(l,1:N_Flux_Pts) = matmul(Ext_S2F,JnS(l,1:N_Soln_Pts))
          enddo

          ! Build Inviscid fluxes on Pencil
          fgF(:,:) = zero
          do i = 1,N_Flux_Pts

            call conserved_to_primitive(ugF(:,i),vgF(:,i),nequations)

            fgF(:,i) = normalflux( vgF(:,i), JnF(:,i), nequations )

          end do

          ! Differentiate Pencil of Fluxes on to solution points
          do k = 1,nequations
            d_fgS(k,1:N_Soln_Pts) = matmul(Dif_F2S,fgF(k,1:N_Flux_Pts))
          enddo

          do i = 1,N_Soln_Pts
            inode   = Pencil_Coord(N_Soln_Pts,jdir,ipen,i)
            divf(:,jdir,inode,ielem) = d_fgS(:,i)
          enddo

        end do

      end do

    case('SSDC')

      do jdir = 1,ndim            ! Directional loop

        do ipen = 1,nodesperface

          !  Grab a pencil of data
          do i = 1,N_Soln_Pts
            inode    = Pencil_Coord(N_Soln_Pts,jdir,ipen,i)  
            ugS(:,i) =   ug(     :,inode,ielem)    !  Neqns
            wgS(:,i) =   wg(     :,inode,ielem)    !  Neqns
            JnS(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
          enddo

          !  Extrapolate data from Solution points to Flux points
          if(N_Soln_Pts == N_Flux_Pts) then
            ugF(:,:) = ugS(:,:)
            wgF(:,:) = wgS(:,:)
            JnF(:,:) = JnS(:,:)
          else
            ugF(:,:) = zero  
            wgF(:,:) = zero  
            do k = 1,nequations
              ugF(k,1:N_Flux_Pts) = matmul(Ext_S2F,ugS(k,1:N_Soln_Pts))
              wgF(k,1:N_Flux_Pts) = matmul(Ext_S2F,wgS(k,1:N_Soln_Pts))
            enddo
            JnF(:,:) = zero  
            do l = 1,ndim
              JnF(l,1:N_Flux_Pts) = matmul(Ext_S2F,JnS(l,1:N_Soln_Pts))
            enddo
          endif

          if( .not. Entropy_Correction) then
            call SS_Euler_Dspec(N_Soln_Pts,N_Flux_Pts,nequations, wgF, JnF, d_fgS)
          else
            call SS_Stabilized_Euler_Dspec(N_Soln_Pts,N_Flux_Pts,nequations, ugF, JnF, d_fgS)
          endif

          do i = 1,N_Soln_Pts
            inode   = Pencil_Coord(N_Soln_Pts,jdir,ipen,i)
            divf(:,jdir,inode,ielem) = d_fgS(:,i)
          enddo

        end do

      end do

    case('SSWENO')

      do jdir = 1,ndim            ! Directional loop

        faceLR = face_pairs(jdir)

        do ipen = 1,nodesperface

           inb  =  0     !  Assumes the pencil is an interior one.

           !  Left interface partner nodes

           if     (ef2e(1,faceLR(1),ielem) < 0) then            ! Partner is a BC (i.e. no partner) 
             inb   = -1
            i_node = Pencil_Coord(N_Soln_Pts,jdir,ipen,1)  
             uL(:) = ug(:,i_node,ielem)
            uLL(:) = ug(:,i_node,ielem)

          else if (ef2e(3,faceLR(1),ielem) /= myprocid) then    ! Partner is off-process
            jnode  = nodesperface*(faceLR(1)-1)+ipen
            gnode  = efn2efn(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack

             uL(:) = ughst(:,gnode)
!           uLL(:) = ughstWENO(:,gnode)
            uLL(:) = ughstWENO_partner(:,gnode)

          else                                              ! Partner is  on-process
             call data_partner_element_serial(ipen,faceLR(1),ielem,k_node,k_elem,k_face,i_node)
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
            i_node = Pencil_Coord(N_Soln_Pts,jdir,ipen,N_Soln_Pts)  
             uR(:) = ug(:,i_node,ielem)
            uRR(:) = ug(:,i_node,ielem)

          else if (ef2e(3,faceLR(2),ielem) /= myprocid) then    ! Partner is off-process
            jnode  = nodesperface*(faceLR(2)-1)+ipen
            gnode  = efn2efn(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack

             uR(:) = ughst(:,gnode)
!           uRR(:) = ughstWENO(:,gnode)
            uRR(:) = ughstWENO_partner(:,gnode)

          else                                              ! Partner is  on-process
             call data_partner_element_serial(ipen,faceLR(2),ielem,k_node,k_elem,k_face,i_node)
            gnode  = nodesperface*(k_face-1) + ipen

             uR(:) = ug(:,k_node,k_elem)
!           uRR(:) = ug(:,WENO_Adjoining_Data(k_node,k_face),k_elem)
            uRR(:) = ugWENO_partner(:,gnode,k_elem)
!           t1 = abs(maxval(uRRT(:)-uRR(:)))
!           if(t1 >= tol1) write(*,*)t1
          endif

          !  Grab a pencil of data
          do i = 1,N_Soln_Pts
            inode    = Pencil_Coord(N_Soln_Pts,jdir,ipen,i)  
            ugS(:,i) =   ug(     :,inode,ielem)    !  Neqns
            wgS(:,i) =   wg(     :,inode,ielem)    !  Neqns
            JnS(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
          enddo
             
          ! First  conditional catches elements that have BC on both faces
          ! Second conditional catches elements that are periodic onto themselves (1 wide)
          if ( (ef2e(1,faceLR(1),ielem) <      0) .and. (ef2e(1,faceLR(2),ielem) <      0)    .or. &
             & (ef2e(2,faceLR(1),ielem) == ielem) .and. (ef2e(2,faceLR(2),ielem) == ielem)  ) then
            call SS_Stabilized_Euler_Dspec(N_Soln_Pts,N_Flux_Pts,nequations, ugS, JnS, d_fgS)
          else
            if(WENO_type == 'Element_WENO') then
              uLL = uL ; uRR = uR ;     !  Ensures consistency of smoothness parameters tau_i
            endif
            call SSWENO4(inb, ugS, uLL, uL, uR, uRR, WENO_Extrp, 1.0_wp, WENO_Extrp, JnS, d_fgS)
          endif

          do i = 1,N_Soln_Pts
            inode   = Pencil_Coord(N_Soln_Pts,jdir,ipen,i)
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
    integer :: i, k, l

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

  subroutine SS_Euler_Dspec(N_S,N_F,N_q, wgF,JnF, dfn)

    use variables
    use referencevariables
    use controlvariables, only: discretization
    use collocationvariables

    implicit none 

    integer,                      intent(in ) :: N_S, N_F, N_q
    real(wp), dimension(N_q,N_F), intent(in ) :: wgF
    real(wp), dimension(  3,N_F), intent(in ) :: JnF

    real(wp), dimension(N_q,N_S), intent(out) :: dfn

    real(wp), dimension(N_q,N_F    )          :: vgF
    real(wp), dimension(N_q,N_F,N_F)          :: SSFlux


    real(wp), dimension(3)                    :: nx
    real(wp), dimension(N_q,N_F)              :: dfnT

    integer                                   :: i,j,k,l
    logical                                   :: flux = .true.
    real(wp), dimension(N_q,0:N_F)            :: fnS


    !  Rotate pencil from conserved variables to primitive variables
    do i=1,N_F
      !call conserved_to_primitive(ugF(:,i),vgF(:,i),N_q)
      call entropy_to_primitive(wgF(:,i),vgF(:,i),N_q)
    enddo

    !  EntropyConsistentFlux is symmetric w.r.t. the Left and Right states
    !  Only the upper half of the matrix are calculated
    !  Interior Diagonals are not calculated because d[[i,i]] = 0 for i /= 1, i/= N_F
    if(.not. flux) then

      SSFlux = 0.0_wp

      do i=1,N_F-1
          do j=i+1,N_F
                    nx(:) = half*(JnF(:,i) + JnF(:,j))
            SSFlux(:,i,j) = EntropyConsistentFlux(vgF(:,i),vgF(:,j),nx,N_q) 
!           SSFlux(:,i,j) = Entropy_KE_Consistent_Flux(vgF(:,i),vgF(:,j),nx,N_q) 
            SSFlux(:,j,i) = SSFlux(:,i,j)
        enddo
      enddo

!       SSFlux(:,  1,  1) = EntropyConsistentFlux(vgF(:,  1),vgF(:,  1),JxF(  1)*nxF(:,  1),N_q)
!       SSFlux(:,N_F,N_F) = EntropyConsistentFlux(vgF(:,N_F),vgF(:,N_F),JxF(N_F)*nxF(:,N_F),N_q)
!       To machine precision, the EntropyConsistenTFlux and the conventional flux are the same
!       when the left and right states are identical

      SSFlux(:,  1,  1) = normalflux( vgF(:,  1), JnF(:,  1), N_q )
      SSFlux(:,N_F,N_F) = normalflux( vgF(:,N_F), JnF(:,N_F), N_q )

      do i=1,N_F
        dfnT(:,i) = zero
        do j=1,N_F
          dfnT(:,i) = dfnT(:,i) + dmat_Flux(i,j)*SSFlux(:,i,j)
        enddo
      enddo
  
      do i=1,N_q
        dfn(i,:) = two * matmul(Int_F2S,dfnT(i,:))
      enddo

    else

      fnS(:,:) = 0.0_wp ;
      do i=1,N_F-1
  
        do j=i+1,N_F
                  nx(:) = half*(JnF(:,i) + JnF(:,j))
          SSFlux(:,i,j) = 2.0_wp* qmat_Flux(i,j)*EntropyConsistentFlux(vgF(:,i),vgF(:,j),nx,N_q) 
!         SSFlux(:,i,j) = 2.0_wp* qmat_Flux(i,j)*diabolical_flux(vgF(:,i), vgF(:,j), wgF(:,i), wgF(:,j),nx,N_q)
!         SSFlux(:,i,j) = 2.0_wp* qmat_Flux(i,j)*Entropy_KE_Consistent_Flux(vgF(:,i),vgF(:,j),nx,N_q) 
        enddo
  
        fnS(:,i) = 0.0_wp
        do k=i+1,N_F
          do l=1,i
            fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
          enddo
        enddo

      enddo

      fnS(:,  0) = normalflux  (vgF(:,  1),JnF(:,  1),N_q)
      fnS(:,N_F) = normalflux  (vgF(:,N_F),JnF(:,N_F),N_q)

      !dfnT(:.:) = 0.0_wp
      do k = 1,N_F
        dfnT(:,k) = pinv_Flux(k) * (fnS(:,k) - fnS(:,k-1))
      enddo

      do i=1,N_q
        dfn(i,:) = matmul(Int_F2S,dfnT(i,:))
      enddo


    endif

    return
  end subroutine SS_Euler_Dspec

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine SS_Stabilized_Euler_Dspec(N_S,N_F,N_q, ugF,JnF, dfn)

    use variables
    use referencevariables
    use controlvariables, only: discretization, flux_entropy_correction
    use collocationvariables

    implicit none 

    integer,                      intent(in ) :: N_S, N_F, N_q
    real(wp), dimension(N_q,N_F), intent(in ) :: ugF
    real(wp), dimension(  3,N_F), intent(in ) :: JnF

    real(wp), dimension(N_q,N_S), intent(out) :: dfn

    real(wp), dimension(N_q,N_F,N_F)          :: SSFlux, DSFlux

    real(wp), dimension(N_q,1:N_F)            :: vgF, wgF

    real(wp), dimension(N_q,0:N_F)            :: fnS, fnD

    real(wp), dimension(3)                    :: nx

    real(wp), dimension(N_q)                  :: bb, ds, delta

    real(wp), parameter                       :: cc2  = 1.0e-24_wp

    integer                                   :: i,j,k,l

    real(wp)                                  :: bbS,dsS,deltaS



    !  Rotate pencil from conserved variables to primitive variables
    do i=1,N_F
      call conserved_to_primitive(ugF(:,i),vgF(:,i),N_q)
      call primitive_to_entropy  (vgF(:,i),wgF(:,i),N_q)
    enddo

    !  Form Entropy Fluxes 0:N_F
    fnS(:,:) = 0.0_wp ; fnD(:,:) = 0.0_wp ;
    do i=1,N_F-1

      do j=i+1,N_F
                nx(:) = half * (JnF(:,i) + JnF(:,j))
        SSFlux(:,i,j) = two  * qmat(i,j)*EntropyConsistentFlux(vgF(:,i),vgF(:,j),nx,N_q) 
!       SSFlux(:,i,j) = two  * qmat(i,j)*Entropy_KE_Consistent_Flux(vgF(:,i),vgF(:,j),nx,N_q) 
        
        select case (flux_entropy_correction)
          case ('normal')
            DSFlux(:,i,j) = two  * qmat(i,j)*normalflux    (half*(vgF(:,i)+vgF(:,j)),nx,N_q)
          case ('Honein-Moin')
            DSFlux(:,i,j) = two  * qmat(i,j)*HoneinMoinFlux(vgF(:,i),vgF(:,j),nx,N_q)
          case default
            write(*,*) 'The flux selected to be use for the entropy correction is unknown.'
            write(*,*) 'Check the subroutine SS_Stabilized_Euler_Dspec()'
            write(*,*) 'Exting....'
            stop
        end select
      enddo

      fnS(:,i) = 0.0_wp ; fnD(:,i) = 0.0_wp ;
      do k=i+1,N_F
        do l=1,i
          fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
          fnD(:,i) = fnD(:,i) + DSFlux(:,l,k)
        enddo
      enddo

    enddo

    fnS(:,  0) = normalflux(vgF(:,  1),JnF(:,  1),N_q) ; fnD(:,  0) = fnS(:,  0)
    fnS(:,N_F) = normalflux(vgF(:,N_F),JnF(:,N_F),N_q) ; fnD(:,N_F) = fnS(:,N_F)

    !  Entropy Correction 
    do i = 1,N_F-1
!         bb(:) = (wgF(:,i+1) - wgF(:,i+0))*(fnS(:,i) - fnD(:,i))
!         ds(:) = sqrt(bb(:)*bb(:) + cc2)
!      delta(:) = (ds(:) - bb(:))/(two*ds(:))
!      fnD(:,i) = fnD(:,i) + delta(:)*(fnS(:,i) - fnD(:,i))  !  Entropy Correction
         bbS = dot_product(wgF(:,i+1)-wgF(:,i+0),fnS(:,i)-fnD(:,i))
         dsS = sqrt(bbS*bbS + cc2)
      deltaS = (dsS - bbS)/(two*dsS)
      fnD(:,i) = fnD(:,i) + deltaS * (fnS(:,i) - fnD(:,i))
    enddo

    do i = 1,N_F
      dfn(:,i) = pinv(i) * (fnD(:,i) - fnD(:,i-1))
    enddo

    return
  end subroutine SS_Stabilized_Euler_Dspec

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

! pure function Pencil_Coord(Ns,jdir,iface,i)
!   
!   integer, intent(in)  :: Ns,jdir,iface,i

!   integer              :: Pencil_Coord

!   Pencil_Coord = 0

!   select case (jdir)
!     case(1)
!       Pencil_Coord =  1 + (mod((iface-0),Ns)-1)*Ns + int((iface-0)/Ns) * Ns*Ns + (i-1)
!     case(2)
!       Pencil_Coord =  1 + (mod((iface-1),Ns)-0)    + int((iface-1)/Ns) * Ns*Ns + (i-1) * Ns
!     case(3)
!       Pencil_Coord =            iface                                          + (i-1) * Ns*Ns
!   end select

! end function

  !============================================================================

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
  
  !============================================================================

  subroutine compute_gradient_entropy_variables()
    
    ! Load modules
    use variables
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: iagrad, jagrad, dagrad, gradmat, pinv, pvol
    
    ! Nothing is implicitly defined
    implicit none

    ! loop indices
    integer :: ielem, jdir, idir
    integer :: inode, jnode, knode, gnode
    integer :: kelem
    integer :: iface
    integer :: i

    ! low and high volumetric element indices
    integer :: iell, ielh
    
    ! normal direction at face
    real(wp) :: nx(3)
    
    ! Temporary arrays for phi
    real(wp), allocatable :: phitmp(:,:)

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    ! phitmp is calculated in computational space
    allocate(phitmp(nequations,ndim))
    
    do ielem = iell, ielh
    ! compute computational gradients of the entropy variables
    !
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
          phig(:,idir,inode,ielem) = phig(:,idir,inode,ielem) &
            + phitmp(:,jdir)*r_x(jdir,idir,inode,ielem)
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
      use controlvariables, only: WENO_Bias
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
      real(wp), dimension(nq,1:4)      :: vint, wint, Kint
      real(wp), dimension(nq, 6)       :: uin, qin, vin
      real(wp), dimension(3,  6)       :: nin
      real(wp), dimension(6)           :: uhat, metr, cav2
      real(wp), dimension(nq,ixd,ixd)  :: SSFlux, DSFlux
      real(wp), dimension(nq,0:ixd)    :: fbarW, fbarC, fnS, fnK

      ! Interpolation coefs and target values
      real(wp), dimension(5,6)         :: Tar, TarP, TarM
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

      real(wp)                         :: a00,w0,w1
      real(wp), dimension(nq)          :: bb ,ds ,delta 
      real(wp)                         :: bbS,dsS,deltaS
      real(wp), parameter              :: cc2 = 1.0e-24_wp
      real(wp), parameter              :: cc1 = 1.0e-12_wp

      real(wp)                         :: x66
      real(wp)                         :: a, b, c, r, x1, x2
      real(wp), dimension(2,2)         :: Hmat, Hinv, Lam, Tmat
!     real(wp), dimension(2,2)         :: eye
      real(wp), dimension(2)           :: bSK
      real(wp), dimension(2,nq)        :: Wmat

      integer                          :: i,j,k,l,n

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
        endif
        
        uhat(:) = 0.0_wp ; metr(:) = 0.0_wp ; cav2(:) = 0.0_wp ;
        do j = 1,6
          call conserved_to_primitive(uin(:,j),vin(:,j),nq)    ! primitives

          uhat(j)  =  abs(dot_product(vin(2:4,j),nin(:,j)))! normal velocity  |u.n|
          metr(j)  = sqrt(dot_product(nin( : ,j),nin(:,j)))! computational coordinate scaling
          cav2(j)  = sqrt(gamma0*vin(5,j)/gM2) * metr(j)   ! Speed of sound * metric scaling
        enddo

        lambda     = maxval(uhat(:) + cav2(:))             !  Max eigenvalue in pencil

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
        dhp = half*(dh+theta*abs(dh))
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
        dhp = half*(dh+theta*abs(dh))
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
      fnS(:,:) = 0.0_wp ; fnK(:,:) = 0.0_wp ;
      do i=1,ixd-1
  
        do j=i+1,ixd
               nxave(:) = half * (nxint(:,i) + nxint(:,j))
          SSFlux(:,i,j) = two  * qmat(i,j)*EntropyConsistentFlux(vint(:,i),vint(:,j),nxave,nq)
!         SSFlux(:,i,j) = two  * qmat(i,j)*Entropy_KE_Consistent_Flux(vint(:,i),vint(:,j),nxave,nq)
!         DSFlux(:,i,j) = two  * qmat(i,j)*HoneinMoinFlux(vint(:,i),vint(:,j),nxave,nq)
        enddo
  
        fnS(:,i) = 0.0_wp ; fnK(:,i) = 0.0_wp ;
        do k=i+1,ixd
          do l=1,i
            fnS(:,i) = fnS(:,i) + SSFlux(:,l,k)
            fnK(:,i) = fnK(:,i) + DSFlux(:,l,k)
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

      implicit none

      integer    :: iell, ielh, ielem, iface, ipen, jnode, kval
      integer    :: elem_face_nodes


      ! low and high volumetric element index
      iell = ihelems(1) ; ielh = ihelems(2) ;

      elem_face_nodes = nodesperface*nfacesperelem

      do ielem = iell, ielh

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

      real(wp), parameter :: c1 = 0.5_wp   !  Guermond's magic constant number one
      real(wp), parameter :: c2 = 1.0_wp   !  Guermond's magic constant number one

      real(wp) :: ev, evmax

      real(wp) :: t1,t2
      integer  :: ielem, iell, ielh, inode

      ! low and high volumetric element index
      iell = ihelems(1) ; ielh = ihelems(2) ;

      do ielem = iell, ielh

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

    function SAT_Inv_Vis_Flux(neq,iface,ielem,vg_On,vg_Off,phig_On,phig_Off,nx,Jx_r,pinv,mut)

      use precision_vars
  
      use controlvariables,     only: Riemann_Diss
      use collocationvariables, only: l01, l00, Sfix, ldg_flip_flop_sign, alpha_ldg_flip_flop


      integer,                       intent(in) :: neq, iface, ielem
      real(wp),  dimension(neq),     intent(in) ::   vg_On,   vg_Off
      real(wp),  dimension(neq,3),   intent(in) :: phig_On, phig_Off
      real(wp),  dimension(3),       intent(in) :: nx
      real(wp),                      intent(in) :: Jx_r, pinv, mut

      real(wp), parameter :: Cevmax          =  1.0_wp
      real(wp), parameter :: deltaU          =  0.1_wp
      real(wp), parameter :: LocalLaxF_factor=  2.0_wp

      real(wp),  dimension(neq,neq)             :: smat,sinv
      real(wp),  dimension(neq,neq)             :: hatc_On, hatc_Off, matrix_ip
      real(wp),  dimension(neq,neq)             :: mattmp


      real(wp),  dimension(neq)                 :: fLLF
      real(wp),  dimension(neq)                 :: wg_On, wg_Off
      real(wp),  dimension(neq)                 :: ug_On, ug_Off
      real(wp),  dimension(neq)                 :: vav, ev, evabs
      real(wp),  dimension(neq)                 :: fn, fstar, fstarV
      real(wp)                                  :: evmax
      real(wp)                                  :: l01_ldg_flip_flop

      real(wp),  parameter                      :: tol = 2.0e-12_wp
      real(wp),  dimension(neq)                 ::   tmpr,tmp2
      real(wp)                                  :: UavAbs, switch
      real(wp)                                  :: t1
      logical                                   :: testing = .false.
 
      integer                                   :: k


      real(wp),  dimension(neq)                 :: SAT_Inv_Vis_Flux

         call primitive_to_entropy (vg_On (:),wg_On (:),neq)
         call primitive_to_entropy (vg_Off(:),wg_Off(:),neq)

         call roeavg( vg_On (:), vg_Off(:), Vav, neq )   
         call CharacteristicDecomp( vav, neq, sinv, smat, ev, nx )      
         evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)

         select case(Riemann_Diss)

           case('Roe')

             fstar = EntropyConsistentFlux(vg_On(:), vg_Off(:), nx, neq ) &
                 & + half * matmul(smat,evabs*matmul(transpose(smat), wg_On(:)-wg_Off(:)) )

           case('LocalLaxF')
             call primitive_to_conserved(vg_On (:),ug_On (:),neq)
             call primitive_to_conserved(vg_Off(:),ug_Off(:),neq)
             fLLF  = half * ( normalflux( vg_On (:), nx, neq )    &
                          +   normalflux( vg_Off(:), nx, neq )    &
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

             fstar = EntropyConsistentFlux(vg_On(:), vg_Off(:), nx, neq ) &
                 & + half * matmul(smat,evabs*matmul(transpose(smat), wg_On(:)-wg_Off(:)) )

             fLLF  = half * ( normalflux( vg_On (:), nx, neq )    &
                          +   normalflux( vg_Off(:), nx, neq )    &
                          +   LocalLaxF_factor*evmax*(ug_On(:)- ug_Off(:))  )

             fstar = switch*fLLF + (1.0_wp-switch) * fstar

         end select

         ! Add the LDG
         l01_ldg_flip_flop = l01*(1.0_wp - ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
         fstarV = normalviscousflux(vg_On (:), phig_On (:,:), nx, neq, mut)  &
              & - normalviscousflux(vg_Off(:), phig_Off(:,:), nx, neq, mut)

         ! Compute the IP penalty contribution, 
         ! c_ii_L matrix    ! cii_R matrix
         hatc_On  = matrix_hatc_node(vg_On (:),nx,nx,neq)
         hatc_Off = matrix_hatc_node(vg_Off(:),nx,nx,neq)

         matrix_ip = 0.5_wp*(hatc_On + hatc_Off) * pinv / Jx_r

         fn = normalflux( vg_On (:), nx, neq )                                  ! (Euler Flux)
         SAT_Inv_Vis_Flux = + (fn - fstar) + l01_ldg_flip_flop*fstarV     &
                            - l00*matmul(matrix_ip,wg_On (:)-wg_Off(:))

          return
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

      real(wp),  dimension(neq)                 ::   ug_On,   ug_Off
 

      real(wp),  dimension(neq)                 :: SAT_Vis_Diss

         call primitive_to_conserved(vg_On (:),ug_On (:),neq)
         call primitive_to_conserved(vg_Off(:),ug_Off(:),neq)

         call roeavg( vg_On (:), vg_Off(:), Vav, neq )   

         call CharacteristicDecomp( vav, neq, sinv, smat, ev, nx )      

         evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)

         SAT_Vis_Diss = half * matmul(smat,evabs*matmul(          sinv , ug_On (:)-ug_Off(:)) )
!                     = half * matmul(smat,evabs*matmul(transpose(smat), wg_On (:)-wg_Off(:)) )

          return
     end function

 end module navierstokes
