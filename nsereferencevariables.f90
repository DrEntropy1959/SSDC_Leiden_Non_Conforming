module nsereferencevariables
  
  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined
  implicit none

  ! Universal gas constant
  real(wp), parameter :: Ru = 8314.472_wp ! (J/kmol-K)
  
  ! Reference temperature for enthalpy calculation
  real(wp), parameter :: Tref = 298.15_wp ! (K)
  
  ! Atmospheric Pressure
  real(wp), parameter :: Patm = 101325.0_wp ! (Pa)
  
  ! Max number of iterations in newton solve
  integer, parameter :: nmaxnewton = 10
  
  ! Newton iteration tolerance for temperature solve
  real(wp), parameter :: temperaturetol = 1.0e-10_wp

  ! Reference values
  ! ================
  real(wp), parameter :: mu0 = 1.0_wp
  real(wp), parameter :: lambda0 = -2.0_wp/3.0_wp*mu0
  real(wp), parameter :: k0  = 1.0_wp
  real(wp) :: MW0
  real(wp) :: cp0
  real(wp) :: U0
  real(wp) :: rho0
  real(wp) :: p0
  real(wp) :: T0
  real(wp) :: gamma0
  real(wp) :: csound0
  real(wp) :: Mach0
  real(wp) :: L0
  real(wp) :: aero_coeffs_surface

  logical :: negTemp = .false.

  integer :: sutherlandpower = 1
  real(wp) :: Re0 = 0.0_wp
  real(wp) :: Re0inv
  real(wp) :: Pr0
  real(wp) :: Sc0
  real(wp) :: Le0

  real(wp) :: vcross, vnoslip
  real(wp) :: sigmanoslip = 1.0_wp

  real(wp) :: membranelocation = 0.0_wp
  real(wp) :: uniformFreestreamAOA = 0.0_wp
  real(wp) :: referenceWaveSpeed = 0.0_wp
  real(wp) :: isothermalWallTemperature = 1.0_wp
  character(180) :: InitialCondition = 'UniformFreeStream'
  logical  ::  entropy_viscosity = .false.

  real(wp) :: rkerrortol

  logical :: viscous    = .false.
  logical :: crossterms = .false.

  ! functions of reference variables (often used in nondimensionalization)
  real(wp) :: gm1, gm1og, gp1og, gm1M2, Msqrtgam, Msqrtgami, gM2, gM2I, gamI, cpT0i, gm1M2I

end module nsereferencevariables
