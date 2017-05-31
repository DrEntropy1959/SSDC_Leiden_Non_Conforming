module navierstokes_Physics

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

  public nse_calcinitialcondition
  public nse_initializesemidiscretization
  public nse_calc_dudt_LGL
  public nse_calc_dudt_LGL_Dense
  public nse_calc_dudt_Gau
  public nse_calc_dudt_Generalized
  public nse_calcembeddederror
  public nse_calcerror
  public nse_calcrhsexplicit
  public nse_calcrhsimplicit
  public nse_calcembeddedspatialerror
  public nse_communicationsetup
  public Navier_Stokes_init
  public dUdV, dVdU, dVdW, dWdV, dWdU, dUdW
  public Flux_Divergence_LGL
  public Flux_Divergence_Gau
  public compute_vorticity_field_elements
  public compute_specific_entropy_elements
  public nse_reconcilestates_LGL
  public nse_reconcilestates_Gau
  public roeavg
  public EntropyConsistentFlux
  public CharacteristicDecomp
  public normalflux
  public isentropicVortexFull
  public viscousShockFull

  public primitive_to_conserved
  public conserved_to_primitive

  public primitive_to_entropy
  public entropy_to_primitive

  public conserved_to_entropy
! public entropy_to_conserved

  public rhalf
  public matrix_hatc_node
  public normalviscousflux
  public compute_explicit_timestep
  public kinetic_energy_element
  public supersonicvortexFull
  public uniformfreestream

  public compute_primitive_conservative_variables_LGL_pts_hexa
  public nse_calcerror_gauss

contains

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

  !============================================================================

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

  end subroutine conserved_to_primitive

  !============================================================================

  pure subroutine conserved_to_entropy(uin,wout,nq)
    ! this routine calculates the primitive variables
    ! from the conserved variables
    use nsereferencevariables, only: gm1M2, gm1og
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! conserved variables
    real(wp), intent(in) :: uin(nq)
    ! entropy variables
    real(wp), intent(out) :: wout(nq)
    ! primitive variables
    real(wp)              :: vtmp(nq)

    ! density
    vtmp(1) = uin(1)
    ! velocity
    vtmp(2:4) = uin(2:4)/uin(1)
    ! temperature
    vtmp(5) = ( uin(5)/uin(1) - gm1M2*0.5_wp*dot_product(vtmp(2:4),vtmp(2:4)) )/(1.0_wp-gm1og)

    ! w_1 = h/T - s - (gamma_0 - 1) M_0^2 u_k u_k/(2T)
    wout(1) = 1.0_wp-0.5_wp*gm1M2*dot_product(vtmp(2:4),vtmp(2:4))/vtmp(5) - specificentropy(vtmp,nq)
    ! w_{k+1} = (gamma_0 - 1) M_0^2 u_k/T, k = 1,2,3
    wout(2:4) = gm1M2*vtmp(2:4)/vtmp(5)
    ! w_5 = -1/T
    wout(5) = -1.0_wp/vtmp(5)


  end subroutine conserved_to_entropy

  !============================================================================

  pure function specificentropy(vin,nq)
    ! this function calculates the specific thermodynamic
    ! entropy using the primitive variable vector
    use nsereferencevariables, only: gm1og
    implicit none
    ! number of equations
    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: vin(nq)

    ! output thermodynamic specific entropy
    real(wp) :: specificentropy

    specificentropy = (1.0_wp-gm1og)*log(vin(5)) - gm1og*log(vin(1))

    return
  end function specificentropy

  !============================================================================

  pure subroutine primitive_to_entropy(vin,wout,nq)
    ! this routine calculates the entropy variables corresponding
    ! to the entropy--entropy flux pair (S,F^i) = (-rho*s,-rho*u_i*s)
    ! using the primitive variables.
    use nsereferencevariables, only:gm1M2
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

!   real(wp) :: w_test(nq)

    continue 

    ! Temperature
    vout(5) = -1.0_wp/win(5)

    ! Velocity components
    vout(2:4) = -1.0_wp/gm1M2*win(2:4)/win(5)

    ! Density
    vout(1) = exp(-1.0_wp*(win(2)**2 + win(3)**2 + win(4)**2 &
                & - 2.0_wp*gm1M2*(win(1)-1.0_wp)*win(5) &
                & + 2.0_wp*gm1M2*(gm1og-1.0_wp)*win(5)*log(-1/win(5))) &
                & /(2.0_wp*gm1M2*gm1og*win(5)))

!    call primitive_to_entropy(vout,w_test,nq)

!    if (abs(maxval(w_test(:) - win(:))) .gt. 1e-12) then
!      write(*,*) w_test(:) - win(:)
!    end if
    
    return
  end subroutine entropy_to_primitive

  !============================================================================

  pure function normalflux(vin,nx,nq)
    ! this function calculates the convective flux in the normal
    ! direction. Note that because we nondimensionalize the pressure
    ! by the reference pressure that there is a nondimensional parameter
    ! in the flux.
    use nsereferencevariables, only: gm1M2, gM2
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

  !============================================================================

  pure function normalviscousflux(vin,phi,nx,nq)
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

    ! output normal viscous flux
    real(wp) :: normalviscousflux(nq)

    ! physical space viscous flux
    real(wp) :: fvl(nq,3)
    ! direction index
    integer :: idir

    continue

    fvl(:,:) = ViscousFlux(vin,phi,nq)

    normalviscousflux = 0.0_wp

    do idir = 1,3
      normalviscousflux = normalviscousflux + Re0inv*fvl(:,idir)*nx(idir)
    end do


    return
  end function normalviscousflux

  !============================================================================

  pure function viscousflux3D(vin,phi,Jx,nq,nd)
    ! this function calculates the viscous flux in 
    ! the three computational space directions based on the primitive
    ! variables and the gradients of the entropy variables,
    ! penalized with an LDC/LDG methodology. The use
    ! of entropy variables ensures stability.
    use nsereferencevariables, only: Re0inv
    
    ! Nothing is implicitly defined
    implicit none
    
    integer, intent(in) :: nq, nd
    ! contravariant vector
    real(wp), intent(in) :: Jx(3,3)
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! entropy variable gradients
    real(wp), intent(in) :: phi(nq,3)

    ! output viscous flux in the computational space directions
    real(wp) :: viscousflux3D(nq,nd)

    ! physical space viscous flux
    real(wp) :: fvl(nq,3)
    ! direction index
    integer :: idir, jdir

    continue

    fvl(:,:) = ViscousFlux(vin,phi,nq)

    viscousflux3D = 0.0_wp

    do jdir = 1,nd
      do idir = 1,3
        viscousflux3D(:,jdir) = viscousflux3D(:,jdir) + Re0inv*fvl(:,idir)*Jx(jdir,idir)
      end do
    end do

    return
  end function viscousflux3D

  !============================================================================

  pure function ViscousFlux(vin,phi,nq)

    use nsereferencevariables, only: gm1M2, Pr0, gm1M2I

    use controlvariables, only : variable_viscosity

    implicit none

    integer, intent(in) :: nq
    ! primitive variables
    real(wp), intent(in) :: vin(nq)
    ! entropy variable gradients
    real(wp), intent(in) :: phi(nq,3)

    real(wp), dimension(nq,3) :: ViscousFlux

    ! thermal conductivity (normalized by kappa0)
    real(wp)             :: kappa_Pr0, KE2, t1,t2,t3
    real(wp)             :: u_phi51,v_phi52,w_phi53
    real(wp), parameter  :: third = 1.0_wp/3.0_wp

    ! dynamic viscosity (normalized by mu0)
    real(wp) :: mu

    continue

    ! Set dynamic viscosity
    if (variable_viscosity .eqv. .true.) then
      mu = sutherland_law(vin(5))
    else
      mu = 1.0_wp
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

  end function viscousflux

  !============================================================================

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
    P  = 0.5_wp * (pl + pr) * gM2I

    ! Temperature average
    Tav = 0.5_wp * (vl(5)+vr(5))

    ! Internal energy
    En = Tav * gamI  + 0.5_wp * gm1M2 * dot_product(vav,vav)

    HoneinMoinFlux(1) = mdot 
    HoneinMoinFlux(2) = mdot*vav(1) + Jx(1) * P
    HoneinMoinFlux(3) = mdot*vav(2) + Jx(2) * P
    HoneinMoinFlux(4) = mdot*vav(3) + Jx(3) * P
    HoneinMoinFlux(5) = mdot*En + gm1og*0.5_wp*(pl*unr+pr*unl)

  end function HoneinMoinFlux

  !============================================================================

  pure function EntropyConsistentFlux(vl,vr,Jx,neqin)
    ! this function calculates the normal entropy consistent
    ! flux based on left and right states of primitive variables.
    ! it is consistent with the nondimensionalization employed
    ! herein and follows directly from the work of Ismail and Roe,
    ! DOI: 10.1016/j.jcp.2009.04.021
    use nsereferencevariables, only: gM2I, gm1og, gp1og, gm1M2
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

    return
  end function EntropyConsistentFlux

  !============================================================================

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

  !============================================================================

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
    integer :: k
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
!   mattmp = 0.0_wp
!   do k = 1,nq
!     mattmp(k,k) = 1.0_wp
!   enddo
!   write(*,*)'dudv.dvdu - eye',maxval(abs(matmul(dUdV(Vav,nq),dVdU(Vav,nq))-mattmp))

    return
  end subroutine CharacteristicDecomp
  
  !============================================================================
  
  subroutine nse_reconcilestates_LGL()
    ! this subroutine is called to update the primitive
    ! and entropy variables. For viscous calculations, the
    ! entropy variable gradients are also updated according
    ! to the LDG approach. 
    
    ! Load modules
    use variables
    use initgrid, only: load_shell_stacks_LGL
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: iagrad_LGL,jagrad_LGL,dagrad_LGL,pinv_LGL,l10, &
                                  & ldg_flip_flop_sign, alpha_ldg_flip_flop, &
                                  & n_LGL_shell, n_LGL_pts_2d, n_LGL_pts_3d

    use mpimod, only: UpdateComm1DGhostData, UpdateComm2DGhostData, &
                    & UpdateComm1D_shell_GhostData, UpdateComm2D_shell_GhostData
    use petscvariables, only: upetsc, phipetsc, ulocpetsc, philocpetsc, &
                            & wglocpetsc_shell, wgpetsc_shell, &
                            & philocpetsc_shell, phipetsc_shell
    
    ! Nothing is implicitly defined
    implicit none

    ! loop indices
    integer :: ielem, jdir, idir
    integer :: inode, jnode, knode, gnode1
    integer :: kelem
    integer :: iface
    integer :: i
    real(wp) :: l10_ldg_flip_flop


    ! low and high volumetric element indices
    integer :: iell, ielh
    ! normal direction at face
    real(wp) :: nx(3)
    ! LDC/LDG coefficient
    ! temporary arrays for phi and delta phi
    real(wp), allocatable :: phitmp(:,:), dphi(:)
    real(wp), allocatable :: wg_On(:), wg_Off(:)

    call UpdateComm1DGhostData(ug,ughst, upetsc, ulocpetsc, &
                            & nequations, n_LGL_pts_3d, ihelems, nghost)

    call UpdateComm1D_shell_GhostData(wg_LGL_shell,wgghst_LGL_shell,wgpetsc_shell, &
                           wglocpetsc_shell,nequations, n_LGL_shell, ihelems, nghost_shell)

    ! low volumetric element index
    iell = ihelems(1)
    ! high volumetric element index
    ielh = ihelems(2)

    if (viscous) then
      ! phitmp is calculated in computational space
      allocate(phitmp(nequations,ndim))
      ! dphi is calculated at faces
      allocate(dphi(nequations))
      allocate(wg_On(nequations),wg_Off(nequations))
      ! LDC/LDG coefficient according to CNG-revisited
      ! loop over elements
      do ielem = iell, ielh
        ! compute computational gradients of the entropy variables
        !
        ! initialize phi
        phig(:,:,:,ielem) = 0.0_wp
        grad_w_jacobian(:,:,:,ielem) = 0.0_wp
        ! loop over every node in element
        do inode = 1, n_LGL_pts_3d
          ! reinitialize computational gradient to zero
          phitmp(:,:) = 0.0_wp
          ! loop over number of dependent elements in gradient
          do i = iagrad_LGL(inode), iagrad_LGL(inode+1)-1
            ! loop over dimensions
            do jdir = 1,ndim
              ! column/node from gradient operator in CSR format in
              ! the jdir-direction corresponding to the coefficient dagrad_LGL(jdir,i)
              jnode = jagrad_LGL(jdir,i)
              ! update gradient using coefficient and entropy variables at appropriate node
              phitmp(:,jdir) = phitmp(:,jdir) + dagrad_LGL(jdir,i)*wg(:,jnode,ielem) 
            end do
          end do
          ! Store gradient of the entropy variables in computational space
          do i = 1, ndim
            grad_w_jacobian(:,i,inode,ielem) = phitmp(:,i)
          enddo
          ! transform to physical space using dxi_jdir/dx_idir
          do jdir = 1,ndim
            do idir = 1,ndim
              phig(:,idir,inode,ielem) = phig(:,idir,inode,ielem) &
                + phitmp(:,jdir)*r_x(jdir,idir,inode,ielem)
            end do
          end do
        end do
         
        ! LDC/LDG penalty on phig
        ! 
        ! loop over only faces
        do iface = 1, nfacesperelem
          if (ef2e(1,iface,ielem) < 0) then
            cycle
          else if (ef2e(3,iface,ielem) /= myprocid) then
            do i = 1,n_LGL_pts_2d
              jnode = n_LGL_pts_2d*(iface-1)+i
              ! corresponding volumetric node for face node
              inode = ifacenodes(jnode)
              ! volumetric node of partner node
              gnode1 = efn2efn_LGL(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack
!             gnode2 = efn2efn(3,jnode,ielem)        ! This is pointing to ghost stack not volumetric stack
              ! outward facing normal of facial node
              nx = facenodenormal_LGL_shell(:,jnode,ielem)

              wg_On(:)  = wg_LGL_shell(:,jnode,ielem)
              wg_Off(:) = wgghst_LGL_shell(:,gnode1)

              ! LDC/LDG penalty value
              l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
              dphi(:) = l10_ldg_flip_flop*pinv_LGL(1)*(wg_On - wg_Off)
              ! add LDC/LDG penalty to each physical gradient using the normal
              do jdir = 1,ndim
                phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
              end do
            end do
          else
            do i = 1,n_LGL_pts_2d
              jnode = n_LGL_pts_2d*(iface-1)+i
              ! corresponding volumetric node for face node
              inode = ifacenodes(jnode)

              ! shell ordering of partner node
              knode = efn2efn_LGL(4,jnode,ielem)

              ! volumetric element of partner node
              kelem = efn2efn_LGL(2,jnode,ielem)

!             nx = facenodenormal(:,jnode,ielem)
              nx = facenodenormal_LGL_shell(:,jnode,ielem)

              wg_On(:)  = wg_LGL_shell(:,jnode,ielem)
              wg_Off(:) = wg_LGL_shell(:,knode,kelem)

              ! LDC/LDG penalty value
              l10_ldg_flip_flop = l10*(1.0_wp + ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)
              dphi(:) = l10_ldg_flip_flop*pinv_LGL(1)*(wg_On - wg_Off)
              ! add LDC/LDG penalty to each physical gradient using the normal
              do jdir = 1,ndim
                phig(:,jdir,inode,ielem) = phig(:,jdir,inode,ielem) + dphi(:)*nx(jdir)
              end do
            end do
          end if
        end do

        do jdir = 1,ndim
          call load_shell_stacks_LGL(phig(:,jdir,:,ielem),phig_LGL_shell_tmp(:,jdir,:,ielem))
        enddo

      end do
      deallocate(phitmp,dphi,wg_On,wg_Off)
    end if

    phig_LGL_shell = phig_LGL_shell_tmp

    call UpdateComm2DGhostData(phig, phighst, phipetsc, philocpetsc, &
                            & nequations, 3, n_LGL_pts_3d, ihelems, nghost)

    call UpdateComm2D_shell_GhostData(phig_LGL_shell, phighst_LGL_shell, phipetsc_shell, philocpetsc_shell, &
                            & nequations, 3, n_LGL_shell, ihelems, nghost_shell)

    return
  end subroutine nse_reconcilestates_LGL

  !============================================================================

  subroutine isentropicVortexFull(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin
    
    real(wp) :: epsvortex
    real(wp) :: Uinf
    real(wp) :: x0,y0
    real(wp) :: f
    real(wp) :: alpha, rin2
   

    y0 = 0.0_wp
    x0 = 0.0_wp
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

  subroutine supersonicvortexFull(Vx,phi,fv,Jx,xin,tin,neqin,nd)

    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin

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

  subroutine Potential_Flow_Around_cylinder(Vx,phi,fv,Jx,xin,tin,neqin,nd)
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
    real(wp), intent(in) :: tin

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

  subroutine ShockVortexInteraction(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    implicit none
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in ) :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3)
    real(wp), intent(in) :: tin
    
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

  subroutine SodsProblemICBC(Vx,phi,fv,Jx,xin,tin,neqin,nd)

    ! Load modules
    use nsereferencevariables
  
    implicit none

    integer,  intent(in)    :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out)   :: fv(neqin)
    real(wp), intent(in )   :: phi(neqin,3)
    real(wp), intent(in)    :: Jx(3)
    real(wp), intent(in)    :: xin(3)
    real(wp), intent(in)    :: tin

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

  subroutine viscousShockFull(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    
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
    real(wp), intent(in) :: xin(3), tin
    
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

!   Check inviscid and viscous terms  MHC, 19/03/2017
!    fnV(:)/Jx_r(inode,ielem) * pinv_LGL(1) = viscousflux1D * pinvFD(1)

    
    return
  end subroutine viscousShockFull

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

  !============================================================================

  subroutine Flux_Divergence_LGL(tin,ielem)
    ! This subroutine calculates elementwise 
    ! the Divergence of the Conservative Flux
    use variables
    use referencevariables
    use nsereferencevariables
    use collocationvariables, only: iagrad_LGL,jagrad_LGL,dagrad_LGL, n_LGL_pts_3d
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem
    real(wp), intent(in) :: tin

    ! indices
    integer :: inode, jdir
    integer :: jnode
    integer :: i

    call Flux_Div_Pencil(ielem)    !  result in divF

    if (viscous) then
      ! loop over all nodes in element
      do inode = 1,n_LGL_pts_3d

        ! calculate viscous flux
        fvg (:,:,inode,ielem) = 0.0_wp

        fvg(:,1:ndim,inode,ielem) = Jx_r(inode,ielem) * &
            viscousflux3D( vg(:,inode,ielem), &
            phig(:,:,inode,ielem), &
            r_x(:,:,inode,ielem), &
            nequations, &
            ndim ) ! (navierstokes)
      end do

      !
      ! calculate divergence of the flux
      ! 

      ! loop over all nodes in the element
      do inode = 1,n_LGL_pts_3d

        ! loop over all nonzero columns in CSR corresponding to this row
        do i = iagrad_LGL(inode), iagrad_LGL(inode+1)-1
          ! loop over each direction
          do jdir = 1,ndim
            ! column/node from gradient operator in CSR format in
            ! the jdir-direction corresponding to the coefficient dagrad_LGL(jdir,i)
            jnode = jagrad_LGL(jdir,i)
            divf(:,jdir,inode,ielem) = divf(:,jdir,inode,ielem)                 & 
             - dagrad_LGL(jdir,i) * fvg(:,jdir,jnode,ielem)
          end do
        end do
      end do
    endif

    return
  end subroutine Flux_Divergence_LGL

  !============================================================================
  ! sat_penalty - Calculates both inviscid and viscous penalty according to the SAT procedure.
  
  subroutine SAT_Penalty_LGL(tin,ielem)
    
    ! Load modules
    use variables
    use referencevariables
    use nsereferencevariables
    use controlvariables, only: heat_entropy_flow_wall_bc
    use collocationvariables, only: pinv_LGL, &
                                  & l01, l00, Sfix, ldg_flip_flop_sign, &
                                  & alpha_ldg_flip_flop,  &
                                  & n_LGL_pts_2d
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
    integer :: i, j

    ! reconstructed flux
    real(wp), allocatable :: fstar(:), fstarV(:)
    ! local normal flux
    real(wp), allocatable :: fn(:), fnV(:)
    ! boundary conservative and primitive states
    real(wp), allocatable :: phistar(:,:)
    ! right and left eigenvector matrices
    real(wp), allocatable, dimension(:,:) :: smat, sinv
    ! eigenvalues
    real(wp), allocatable, dimension(:) :: ev, evabs
    ! average state
    real(wp), allocatable, dimension(:) :: vav
    ! normal vector
    real(wp) :: nx(3)
    ! Lax-Freidrich max Eigenvalue
    real(wp) :: evmax
    ! penalty parameter for interfaces
    real(wp) :: t1

    real(wp), dimension(nequations)            :: tmpr
    logical                                    :: testing = .false.

    real(wp), dimension(nequations,nequations) :: hatc_side_1, hatc_side_2, &
                                                & matrix_ip

    real(wp), parameter :: mirror          = -1.0_wp
    real(wp), parameter :: no_slip         = -0.0_wp
    real(wp), parameter :: tol             =  2.0e-12_wp
    real(wp), parameter :: stab_strength   = 1.1
    real(wp), parameter :: Cevmax          = 1.0
    real(wp), parameter :: bc_pen_strength = 2.0_wp

    real(wp), dimension(3)                :: unit_normal

    real(wp), dimension(nequations)       :: prim_ghost_euler
    real(wp), dimension(nequations)       :: entr_ghost_euler

    real(wp), dimension(3)                :: normal_vel, tangent_vel

    real(wp), dimension(nequations)       :: f_viscous_normal_ghost

    real(wp), dimension(nequations)       :: prim_ghost_adiabatic
    real(wp), dimension(nequations)       :: entr_ghost_adiabatic
    real(wp), dimension(nequations)       :: prim_ref

    real(wp), dimension(nequations)       ::   ug_On,  ug_Off
    real(wp), dimension(nequations)       ::   vg_On,  vg_Off
    real(wp), dimension(nequations)       ::   wg_On,  wg_Off
    real(wp), dimension(nequations,3)     :: phig_On,phig_Off
    real(wp), dimension(nequations)       :: tmpvec
    real(wp), dimension(nequations)       ::   SAT_Pen

    integer                               :: kknode, kkelem
    
    real(wp)    :: grad_prim_int(nequations,3)
    real(wp)    :: grad_entr_ghost(nequations,3)
    real(wp)    :: grad_entr_int_normal(nequations,3),grad_entr_int_tangent(nequations,3)
    
    real(wp), dimension(3) :: ref_vel_vector

    real(wp) :: l01_ldg_flip_flop
    real(wp) :: tmp


    ! allocate local arrays
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
    gsat (:,:,ielem) = 0.0_wp
    gsatI(:,:,ielem) = 0.0_wp
    gsatV(:,:,ielem) = 0.0_wp
       
       
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
          do i = 1,n_LGL_pts_2d

            ! Volumetric node index corresponding to face and node on face indices
            inode = kfacenodes(i,iface)
          
            ! Facial index corresponding to face and node on face indices
            jnode = (iface-1)*n_LGL_pts_2d + i
          
            ! Outward facing normal of facial node
            nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)

            wg_On(:)  = wg_LGL_shell(:,jnode,ielem)
            call entropy_to_primitive(wg_On (:),vg_On(:) ,nequations)
            call primitive_to_conserved(vg_On (:),ug_On(:) ,nequations)

            vg_Off(:) = vg_On(:)
            phistar   = phig(:,:,inode,ielem)
            
            call BoundaryCondition(vg_Off(:),phistar,fnV,nx,xg(:,inode,ielem),tin, nequations,ndim)
            call primitive_to_entropy  (vg_Off(:),wg_Off(:), nequations) 
            call primitive_to_conserved(vg_Off(:),ug_Off(:), nequations) 

            call roeavg( vg_On(:), vg_Off(:), Vav, nequations )         

            call CharacteristicDecomp( vav, nequations, sinv, smat, ev, nx ) 

            evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)
!           evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)

            fn = normalflux( vg_On(:), nx, nequations )                     ! (Euler Flux)

            fstar = EntropyConsistentFlux(vg_On(:), vg_Off(:), nx, nequations ) & ! (Entropy Flux)
                & + half * matmul(smat,evabs*matmul(          sinv , ug_On(:)-ug_Off(:)) )
!               & + half * matmul(smat,evabs*matmul(transpose(smat), wg_On(:)-wg_Off(:)) )

            fstarV = normalviscousflux( vg_On(:), phig(:,:,inode,ielem), nx, nequations) &
                 & - fnV(:)

            ! Compute the IP penalty contribution, i.e. M (u-v), where M, in the
            ! c_ii_side_1 matrix and c_ii_side_2 matrix
            hatc_side_1 = matrix_hatc_node(vg_On (:),nx,nx,nequations)
            hatc_side_2 = matrix_hatc_node(vg_Off(:),nx,nx,nequations)

            ! IP penalty matrix
            matrix_ip = 0.5_wp*(hatc_side_1 + hatc_side_2)*pinv_LGL(1) / Jx_r(inode,ielem)

            gsat(:,inode,ielem) = gsat(:,inode,ielem) &
              & + pinv_LGL(1)* ( (fn - fstar) + l01*fstarV*bc_pen_strength ) !&
!             & - pinv_LGL(1)* l00*matmul(matrix_ip,wg_On(:) - wg_Off(:))

          end do
        
        else  !  Entropy stable solid wall BC's 

          ! Loop over each node on the face
          do i = 1,n_LGL_pts_2d

            ! Volumetric node index corresponding to face and node on face indices
            inode = kfacenodes(i,iface)
          
            ! Facial index corresponding to face and node on face indices
            jnode = (iface-1)*n_LGL_pts_2d + i
          
            ! Outward facing normal of facial node
            nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)

            ! Unit normal direction
            unit_normal = nx/Magnitude(nx)

            ! Normal velocity
            normal_vel = dot_product(vg(2:4,inode,ielem),unit_normal)*unit_normal

            ! Tangent velocity
            tangent_vel = vg(2:4,inode,ielem) - normal_vel

            prim_ghost_euler(1) = vg(1,inode,ielem)
            prim_ghost_euler(2) = tangent_vel(1) - normal_vel(1)
            prim_ghost_euler(3) = tangent_vel(2) - normal_vel(2)
            prim_ghost_euler(4) = tangent_vel(3) - normal_vel(3)
            prim_ghost_euler(5) = vg(5,inode,ielem) 

            ! Compute the entropy variables in the ghost node for imposing the
            ! nonpenetration condition
            call primitive_to_entropy(prim_ghost_euler,entr_ghost_euler,nequations) 

            ! Compute the roe average state of the primitive variables
            call roeavg(vg(:,inode,ielem),prim_ghost_euler,vav,nequations)

            ! Compute characteristic decomposition
            call CharacteristicDecomp(vav,nequations,sinv,smat,ev,nx)

            evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)
!           evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)
            ! ==  Eigen values/vectors

            ! ==  Fluxes
            fn = normalflux( vg(:,inode,ielem), nx, nequations )                     ! (Euler Flux)

            fstar = EntropyConsistentFlux(vg(:,inode,ielem), prim_ghost_euler, nx, nequations ) & ! (Entropy Flux)
                  &      + half * matmul(smat,evabs*matmul(transpose(smat), wg(:,inode,ielem)-entr_ghost_euler) )
    
            ! Entropy consistent flux + upwinding

            ! Compute the boundary data in primitive variables for imposing the
            ! no-slip isothermal wall BC
            prim_ghost_adiabatic(1) = vg(1,inode,ielem)
            prim_ghost_adiabatic(2) = 0.0_wp
            prim_ghost_adiabatic(3) = 0.0_wp
            prim_ghost_adiabatic(4) = 0.0_wp
            prim_ghost_adiabatic(5) = vg(5,inode,ielem)

            ! Compute the entropy variables in the ghost node for imposing the
            ! no-slip isothermal wall BC
            call primitive_to_entropy(prim_ghost_adiabatic,entr_ghost_adiabatic,nequations)

            fstarV = normalviscousflux( vg(:,inode,ielem), phig(:,:,inode,ielem), nx, nequations) 

            ! Grad(V) = dVdW Grad(W)  = dvdw phi
            grad_prim_int = MatMul(dVdW(vg(:,inode,ielem),nequations),phig(:,:,inode,ielem))

            ! Set gradient of the primitive variables in the ghost node
            ! =========================================================
            grad_entr_ghost(1,:) = phig(1,:,inode,ielem) !grad_prim_int(1,:) ! grad(rho_ghost) = grad(rho_in)
            grad_entr_ghost(2,:) = phig(2,:,inode,ielem) !grad_prim_int(2,:) ! grad(u_ghost) = grad(u_in)
            grad_entr_ghost(3,:) = phig(3,:,inode,ielem) !grad_prim_int(3,:) ! grad(v_ghost) = grad(v_in)
            grad_entr_ghost(4,:) = phig(4,:,inode,ielem) !grad_prim_int(4,:) ! grad(w_ghost) = grad(w_in)

            ! Normal component of grad_prim_int(V)
            do j = 1, nequations
              grad_entr_int_normal(j,:) = dot_product(phig(j,:,inode,ielem),unit_normal(:))*unit_normal(:) 
            end do

            grad_entr_int_tangent(:,:) = phig(:,:,inode,ielem) - grad_entr_int_normal
  
            grad_entr_ghost(5,:) = grad_entr_int_tangent(5,:) + heat_entropy_flow_wall_bc*unit_normal(:)/vg(5,inode,ielem)

            ! Compute normal viscous flux arising from the ghost point
            f_viscous_normal_ghost = normalviscousflux(vg(:,inode,ielem),grad_entr_ghost,nx,nequations)      
            
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
            matrix_ip = hatc_side_2 * pinv_LGL(1) / Jx_r(inode,ielem)

            call primitive_to_entropy(vg(:,inode,ielem),wg_On(:),nequations)
              
            ! Compute the penalty term
            gsatI(:,inode,ielem) = gsatI(:,inode,ielem) + pinv_LGL(1)*(fn - fstar)
            gsatV(:,inode,ielem) = gsatV(:,inode,ielem)       &
              & - pinv_LGL(1)*(fstarV-f_viscous_normal_ghost) &
              & - pinv_LGL(1)*(1.0_wp*matmul(matrix_ip,wg_On(:)-entr_ghost_adiabatic))
            gsat(:,inode,ielem) = gsat(:,inode,ielem) &
              & + pinv_LGL(1)*(fn - fstar) &
              & - pinv_LGL(1)*(fstarV-f_viscous_normal_ghost) &
              & - pinv_LGL(1)*(1.0_wp*matmul(matrix_ip,wg_On(:)-entr_ghost_adiabatic))

          end do

        end if


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       Off Processor Contributions to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      else if (ef2e(3,iface,ielem) /= myprocid) then
        
        ! This is a parallel interface
        do i = 1, n_LGL_pts_2d
          jnode = n_LGL_pts_2d*(iface-1)+i
          
          ! Volumetric node index corresponding to facial node index
          inode = ifacenodes(jnode)
          
          ! Index in ghost
          gnode = efn2efn(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack
          
          ! Element index of partner node
          kelem = efn2efn(2,jnode,ielem)
          
          ! Outward facing normal of facial node
          nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)

          kknode = efn2efn_LGL(3,jnode,ielem)  ! This is pointing to ghost stack not volumetric stack
          wg_On (:) = wg_LGL_shell(:,jnode,ielem)
          wg_Off(:) = wgghst_LGL_shell(:,kknode)

          call entropy_to_primitive(wg_On (:) ,vg_On (:) ,nequations)
          call entropy_to_primitive(wg_Off(:) ,vg_Off(:) ,nequations)

          phig_On (:,:) = phig(:,:,inode,ielem)
          phig_Off(:,:) = phighst(:,:,gnode)

          SAT_Pen =  SAT_Inv_Vis_Flux( nequations,iface,ielem,    &
                                     & vg_On,vg_Off,wg_On,wg_Off, &
                                     & phig_On,phig_Off,          &
                                     & nx,Jx_r(inode,ielem),pinv_LGL(1))

          gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv_LGL(1) * SAT_Pen

        end do
      
      else
        do i = 1,n_LGL_pts_2d
        
          ! Index in facial ordering
          jnode = n_LGL_pts_2d*(iface-1)+i
            
          ! Volumetric node index corresponding to facial node index
          inode = ifacenodes(jnode)
            
          ! Volumetric index of partner node
          knode = efn2efn(1,jnode,ielem)
            
          ! Element index of partner node
          kelem = efn2efn(2,jnode,ielem)
          
          ! Outward facing normal of facial node
          nx = Jx_r(inode,ielem)*facenodenormal(:,jnode,ielem)
          
          kknode = efn2efn_LGL(4,jnode,ielem)
          kkelem = efn2efn_LGL(2,jnode,ielem)
          wg_On(:)  = wg_LGL_shell(:,jnode,ielem)
          wg_Off(:) = wg_LGL_shell(:,kknode,kkelem)

          call entropy_to_primitive(wg_On (:) ,vg_On (:) ,nequations)
          call entropy_to_primitive(wg_Off(:) ,vg_Off(:) ,nequations)

          phig_On (:,:) = phig(:,:,inode,ielem)
          phig_Off(:,:) = phig(:,:,knode,kelem)

          SAT_Pen =  SAT_Inv_Vis_Flux( nequations,iface,ielem,    &
                                     & vg_On,vg_Off,wg_On,wg_Off, &
                                     & phig_On,phig_Off,          &
                                     & nx,Jx_r(inode,ielem),pinv_LGL(1))

          gsat(:,inode,ielem) = gsat(:,inode,ielem) + pinv_LGL(1) * SAT_Pen

        end do
      end if
    end do
    
    
    ! Deallocate memory
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
  end subroutine SAT_Penalty_LGL

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
      case default
        InitialSubroutine => UniformFreeStream
        write(*,*)'no match for InitialCondition specificiation'
        write(*,*)'using default UniformFreeStream'
     end select

     return
  end subroutine set_InitialCondition

  pure subroutine UniformFreeStream(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin

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

  pure subroutine GFIT(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin

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

  pure subroutine QuiescentFluid(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin

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

  pure subroutine taylorGreen(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin
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

  pure subroutine InviscidWall(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer,  intent(in)    :: neqin, nd

    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out)   :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in)    :: Jx(3)
    real(wp), intent(in)    :: xin(3)
    real(wp), intent(in)    :: tin

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

  pure subroutine SymmetryPlane(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin

    real(wp), parameter  :: Mirror = -1.0_wp
    real(wp), dimension(3)   :: VNormal, Vtangent, JxNormal

    JxNormal = Jx(:) / Magnitude(Jx)

    VNormal(:)   = dot_product(Vx(2:4),JxNormal(:)) * JxNormal(:)
    VTangent(:)  = Vx(2:4) - VNormal(:)

    Vx(1)   = Vx(1)
    Vx(2:4) = VTangent(:) + Mirror*VNormal(:)
    Vx(5)   = Vx(5)

    fv = Mirror*normalviscousflux(Vx, phi, Jx, neqin)
    
    return
  end subroutine SymmetryPlane

  !============================================================================
  
  pure subroutine NoSlipWallAdiabatic(prim_int,grad_entr_int,f_v,normal,pos,time, &
      & n_eq,n_d)
    
    ! Load modules
    use nsereferencevariables
    
    ! Nothing is implicitly defined
    implicit none
    
    integer,  intent(in)    :: n_eq, n_d
    real(wp), intent(out)   :: f_v(n_eq)
    real(wp), intent(inout) :: prim_int(n_eq)
    real(wp), intent(in)    :: grad_entr_int(n_eq,3)
    real(wp), intent(in)    :: normal(3)
    real(wp), intent(in)    :: pos(3), time

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
      & grad_prim_int_tangent+grad_prim_int_normal),normal,n_eq)


    return
  end subroutine NoSlipWallAdiabatic

  !============================================================================

  subroutine SubsonicOutflow(Vx,phi,fv,Jx,xin,tin,neqin,nd)
   use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin

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

  pure subroutine Constant_rhoS(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin

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

  subroutine BLayerProfile(Vx,phi,fv,Jx,xin,tin,neqin,nd)
    use nsereferencevariables
    integer, intent(in) :: neqin, nd
    real(wp), intent(inout) :: Vx(neqin)
    real(wp), intent(out) :: fv(neqin)
    real(wp), intent(in)    :: phi(neqin,3)
    real(wp), intent(in) :: Jx(3)
    real(wp), intent(in) :: xin(3), tin
    
    Vx = Vx
    fv = zero

  end subroutine BLayerProfile

  !============================================================================
  ! compute_explicit_timestep- given a CFL, computes the max timestsp 

  subroutine compute_explicit_timestep(dt_global)

    ! Load modules
    use variables
    use controlvariables
    use collocationvariables, only: pvol_Gau_pts_hex, n_Gau_pts_3d
    use referencevariables
    use nsereferencevariables
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
    integer :: m

    real(wp)               :: Lngth, a0, dt0, dt_min, dt_global_max
    real(wp)               :: tI, tV
    real(wp), dimension(3) :: sq, ucon

    continue

    ! Low and high  volumetric element index
    elem_low  = ihelems(1) ;  elem_high = ihelems(2)

    ! Compute gradient of the velocity components
    ! ===========================================
    dt_global_max = dt_global*1.1_wp
    dt_min = 100.0_wp

    do i_elem = elem_low, elem_high
      do i_node = 1, n_Gau_pts_3d
                
          call conserved_to_primitive(ug_Gau_pts_hex(:,i_node,i_elem),vg_Gau_pts_hex(:,i_node),nequations)

          Lngth = pvol_Gau_pts_hex(i_node)**(third)
          a0  = sqrt(abs(gamma0*vg_Gau_pts_hex(5,i_node)/gM2))

          sq(1) = magnitude(r_x_Gau_pts_hex(1,:,i_node,i_elem))
          sq(2) = magnitude(r_x_Gau_pts_hex(2,:,i_node,i_elem))
          sq(3) = magnitude(r_x_Gau_pts_hex(3,:,i_node,i_elem))

          ucon(1) = dot_product(r_x_Gau_pts_hex(1,:,i_node,i_elem),vg_Gau_pts_hex(2:4,i_node))
          ucon(2) = dot_product(r_x_Gau_pts_hex(2,:,i_node,i_elem),vg_Gau_pts_hex(2:4,i_node))
          ucon(3) = dot_product(r_x_Gau_pts_hex(3,:,i_node,i_elem),vg_Gau_pts_hex(2:4,i_node))

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

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine Flux_Div_Pencil(ielem)
    ! This subroutine calculates elementwise 
    ! the Divergence of the Conservative Flux
    use variables
    use referencevariables
    use controlvariables, only: discretization, Entropy_Correction
    use collocationvariables
    implicit none
    ! local time of evaluation (for RK schemes, this is the stage time)
    integer , intent(in) :: ielem

    ! indices
    integer :: inode, jdir, iface
    integer :: i, k

    real(wp), dimension(nequations,n_LGL_pts_1d) :: wgF, vgF, fgF, d_fgF
    real(wp), dimension(         3,n_LGL_pts_1d) :: JnF, xgF

    select case(discretization)

    case('SpecColl')

      do jdir = 1,ndim            ! Directional loop

        do iface = 1,n_LGL_pts_2d

          fgF(:,:) = zero
          !  Grab a pencil of data and build Inviscid fluxes on Pencil
          do i = 1,n_LGL_pts_1d
            inode    = Pencil_Coord(n_LGL_pts_1d,jdir,iface,i)  
            wgF(:,i) =   wg(     :,inode,ielem)  ; call entropy_to_primitive(wgF(:,i),vgF(:,i),nequations)
            JnF(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
            fgF(:,i) = normalflux( vgF(:,i), JnF(:,i), nequations )
          enddo

          ! Differentiate Pencil of Fluxes on to solution points
          do k = 1,nequations
            d_fgF(k,1:n_LGL_pts_1d) = matmul(dmat_LGL,fgF(k,1:n_LGL_pts_1d))
          enddo

          do i = 1,n_LGL_pts_1d
            inode   = Pencil_Coord(n_LGL_pts_1d,jdir,iface,i)
            divf(:,jdir,inode,ielem) = d_fgF(:,i)
          enddo

        end do

      end do

    case('SSDC')

      do jdir = 1,ndim            ! Directional loop

        do iface = 1,n_LGL_pts_2d

          !  Grab a pencil of data
          do i = 1,n_LGL_pts_1d
            inode    = Pencil_Coord(n_LGL_pts_1d,jdir,iface,i)  
            wgF(:,i) =   wg(     :,inode,ielem)    !  Neqns
            JnF(:,i) = Jx_r(inode,ielem)*r_x(jdir,:,inode,ielem)
            xgF(:,i) = xg(:,inode,ielem)
          enddo

          if( .not. Entropy_Correction) then
!           call SS_Euler_Dspec_Test(n_LGL_pts_1d,nequations, wgF, JnF, d_fgF,xgF)
            call SS_Euler_Dspec(n_LGL_pts_1d,nequations, wgF, JnF, d_fgF)
          else
            call SS_Stabilized_Euler_Dspec(n_Gau_pts_1d,n_LGL_pts_1d,nequations, wgF, JnF, d_fgF)
          endif

          do i = 1,n_LGL_pts_1d
            inode   = Pencil_Coord(n_LGL_pts_1d,jdir,iface,i)
            divf(:,jdir,inode,ielem) = d_fgF(:,i)
          enddo

        end do

      end do

    end select

    return
  end subroutine Flux_Div_Pencil

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
    
    integer, intent(in) :: n_eq
    real(wp), intent(in) :: n_i(3), n_j(3)
    real(wp), intent(in) :: v_in(n_eq)
    real(wp) :: matrix_hatc_node(n_eq,n_eq)
    real(wp), dimension(n_eq,n_eq) :: mat
    real(wp) :: con1, con2, con3, con4
    real(wp) :: u, v, w, T
    real(wp) :: mu

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

    ! Jacobi matrix wrt to the conserved variables 
    matrix_hatc_node = mat(:,:)

    return
  end function matrix_hatc_node

  !============================================================================
  
  !============================================================================
  ! set_boundary_conditions - Set the pointer to the boundary condition 
  ! subroutine. 

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

    function SAT_Inv_Vis_Flux(neq,iface,ielem,vg_On,vg_Off,wg_On,wg_Off,phig_On,phig_Off,nx,Jx_r,pinv)


    use collocationvariables, only: l01, l00, Sfix, ldg_flip_flop_sign, alpha_ldg_flip_flop

      integer,                       intent(in) :: neq, iface, ielem
      real(wp),  dimension(neq),     intent(in) ::   vg_On,   vg_Off
      real(wp),  dimension(neq),     intent(in) ::   wg_On,   wg_Off
      real(wp),  dimension(neq,3),   intent(in) :: phig_On, phig_Off
      real(wp),  dimension(3),       intent(in) :: nx
      real(wp),                      intent(in) :: Jx_r, pinv

      real(wp),  dimension(neq,neq)             :: smat,sinv
      real(wp),  dimension(neq,neq)             :: hatc_On, hatc_Of, matrix_ip
      real(wp),  dimension(neq,neq)             :: mattmp


      real(wp),  dimension(neq)                 :: vav, ev, evabs
      real(wp),  dimension(neq)                 :: fn, fstar, fstarV
      real(wp)                                  :: evmax
      real(wp)                                  :: l01_ldg_flip_flop

      real(wp),  parameter                      :: tol = 2.0e-12_wp
      real(wp),  dimension(neq)                 ::   ug_On,   ug_Off
      real(wp),  dimension(neq)                 ::   tmpr,tmp2
      real(wp)                                  :: t1
      logical                                   :: testing = .false.
 
      integer                                   :: k


      real(wp),  dimension(neq)                 :: SAT_Inv_Vis_Flux

         call primitive_to_conserved(vg_On (:),ug_On (:),neq)
         call primitive_to_conserved(vg_Off(:),ug_Off(:),neq)

         call roeavg( vg_On (:), vg_Off(:), Vav, neq )   

         call CharacteristicDecomp( vav, neq, sinv, smat, ev, nx )      

         evmax = maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax)
!        evmax = Cevmax*maxval( abs(ev(:)) )  ; evabs(:) = sqrt(ev(:)*ev(:) + Sfix*evmax*evmax)

         fn = normalflux( vg_On (:), nx, neq )                                  ! (Euler Flux)

         ! Entropy consistent flux + upwinding
!  HACK
!        fstar = EntropyConsistentFlux(vg_On (:), vg_Off(:), nx, neq )  ! (Entropy Flux)
         fstar = EntropyConsistentFlux(vg_On (:), vg_Off(:), nx, neq ) & ! (Entropy Flux)
               & + half * matmul(smat,evabs*matmul(          sinv , ug_On (:)-ug_Off(:)) )
!              & + half * matmul(smat,evabs*matmul(transpose(smat), wg_On (:)-wg_Off(:)) )
!        fstar = 0.5_wp *( normalflux( vg_On(:),nx,neq) + normalflux( vg_Off(:),nx,neq) )  &
!              & + half * matmul(smat,evabs*matmul(          sinv , ug_On (:)-ug_Off(:)) )
!        fstar = 0.5_wp *( normalflux( vg_On(:),nx,neq) + normalflux( vg_Off(:),nx,neq) )  &
!            & + 0.5_wp * evmax * (ug_On (:)-ug_Off(:)) 
!  HACK

         fstarV = normalviscousflux(vg_On (:), phig_On (:,:), nx, neq)  &
              & - normalviscousflux(vg_Off(:), phig_Off(:,:), nx, neq)

         ! Compute the IP penalty contribution, i.e. M (u-v), where M, in the
         ! c_ii_L matrix    ! cii_R matrix
         hatc_On = matrix_hatc_node(vg_On (:),nx,nx,neq)
         hatc_Of = matrix_hatc_node(vg_Off(:),nx,nx,neq)

         ! IP penalty matrix
         matrix_ip = 0.5_wp*(hatc_On + hatc_Of) * pinv / Jx_r

         ! Add the LDG and IP terms to the penalty
         l01_ldg_flip_flop = l01*(1.0_wp - ldg_flip_flop_sign(iface,ielem)*alpha_ldg_flip_flop)

         SAT_Inv_Vis_Flux =   +  fn - fstar                           &
                            & + l01_ldg_flip_flop*fstarV              &
                            & - l00*matmul(matrix_ip,wg_On (:)-wg_Off(:))

!  TESTING
         if( testing ) then     !  This conditional test the condition dF = A dU  for Roe flux

!           mattmp = 0.0_wp
!           do k = 1,neq
!             mattmp(k,k) = 1.0_wp
!           enddo
!           write(*,*)'dudv.dvdu - eye',maxval(abs(matmul(dUdW(Vav,neq),dWdU(Vav,neq))-mattmp))

!           call primitive_to_conserved(vg_On (:),tmpr(:),neq)
!           call conserved_to_primitive(tmpr(:),tmp2(:),neq)
!           t1 = sqrt(dot_product(tmp2-vg_On,tmp2-vg_On)/neq)
!           if(t1 >= tol) write(*,*)'v_2_u u_2_v ', t1
!           call primitive_to_entropy(vg_On (:),tmpr(:),neq)
!           call entropy_to_primitive(tmpr(:),tmp2(:),neq)
!           t1 = sqrt(dot_product(tmp2-vg_On,tmp2-vg_On)/neq)
!           if(t1 >= tol) write(*,*)'v_2_w w_2_v ', t1
            
            tmpr(:) = (+ normalflux( vg_On (:), nx, neq )                    &
                       - normalflux( vg_Off(:), nx, neq )  )                 &
                       - matmul(smat,ev*matmul(sinv,(ug_On (:) - ug_Off(:))))
!                      - matmul(smat,ev*matmul(Transpose(smat),(wg_On (:) - wg_Off(:))))
            t1 = sqrt(dot_product(tmpr,tmpr)/neq)
            if(t1 >= tol) write(*,*)'Roe_unit error', t1
            t1 = maxval(abs(matmul(smat,transpose(smat)) - dUdW(vav,neq)))
            if(t1 >= tol) write(*,*)'smat.smat^T - dqdw', t1

            t1 = maxval(abs(matmul(          sinv ,ug_On (:)-ug_Off(:)) -  &
                            matmul(Transpose(smat),wg_On (:)-wg_Off(:)))) 
            if(t1 >= tol) write(*,*)'Sinv dq - Smat^T dw', t1
            t1 = maxval(abs(matmul(dUdW(vav,neq),wg_On (:)-wg_Off(:)) - (ug_On (:) - ug_Off(:))))
            if(t1 >= tol) write(*,*)'dqdw dw - dq', t1
            t1 = maxval(abs(matmul(dWdU(vav,neq),ug_On (:)-ug_Off(:)) - (wg_On (:) - wg_Off(:))))
            if(t1 >= tol) write(*,*)'dWdq dq - dw', t1
          endif

          return
     end function

  !============================================================================

    function SAT_Vis_Diss(neq,vg_On,vg_Off,nx)


    use collocationvariables, only: Sfix

      integer,                       intent(in) :: neq, iface, ielem
      real(wp),  dimension(neq),     intent(in) ::   vg_On,   vg_Off
      real(wp),  dimension(neq,3),   intent(in) :: phig_On, phig_Off
      real(wp),  dimension(3),       intent(in) :: nx
      real(wp),                      intent(in) :: Jx_r

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

  !============================================================================

end module navierstokes_Physics

