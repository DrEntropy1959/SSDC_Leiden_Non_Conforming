module SSWENOvariables

        use precision_vars

        implicit none

        ! Penalty Coefficients
        real(wp), parameter    :: bias = 0.75_wp  !  WENO Stencil biasing: 0: maximum biasing,  1: no biasing
        real(wp), parameter    :: Ceps = 1.00_wp  !  scaling of the eps value in stencil biasing
! Original
!       real(wp)               :: dx0  = 0.10_wp  !  Average physical grid spacing in x direction Vol / N_Elem

!       real(wp)               :: WENO_Extrp = 1.0_wp  !  Scale the Extrapolation into adjoining element  (1.0 is full distance)
! NAIL's Version    Note that  dx0 causes EXTREME STIFFNESS as it is decreased
        real(wp)               :: dx0  = 0.01_wp  !  Average physical grid spacing in x direction Vol / N_Elem
                                                  !  dx0 causes EXTREME STIFFNESS as it is decreased
  
        real(wp)               :: WENO_Extrp = 0.1_wp  !  Nails modification that destroys accuracy of smooth flows. (Bugs)

        character(120)         :: WENO_type = 'Neighbr_WENO'   !  'Element_WENO'

        real(wp), allocatable, dimension(:)   :: tau_cfL, tau_cfR
        real(wp), allocatable, dimension(:,:) :: IM2,IM1,IC,IP1,IP2

end module SSWENOvariables
