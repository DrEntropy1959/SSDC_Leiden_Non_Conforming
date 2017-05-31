module SSWENO_Routines

        use precision_vars
        implicit none

        private
        public  WENO_Coeffs
        public Negative_Density_Removal

  contains

        subroutine WENO_Coeffs

        use SSWENOvariables

        implicit none

        real(wp), parameter :: sqrt5  = sqrt(5.0_wp)
        real(wp), parameter :: c1     = 15.0_wp/4.0_wp
        real(wp), parameter :: c2     = 15.0_wp/4.0_wp*sqrt5

        real(wp), dimension(5), parameter :: tau_coeffT0 = (/0.0_wp,-c1,c2,-c2,c1/)
        real(wp), dimension(5), parameter :: tau_coeffT1 = (/-c1,c2,-c2,c1,0.0_wp/)
        real(wp), dimension(5)            :: tau_cL4,tau_cR4

!       WENO smoothness indicators 
        allocate(tau_cfL(5),tau_cfR(5))

!       WENO Interpolation coefficients

        allocate(IM2(5,2),IM1(5,2),IC(5,2),IP1(5,2),IP2(5,2))

        IM2(1,1)= 1.0_wp
        IM2(2,1)= 0.69849716760417542931628443_wp
        IM2(2,2)= 0.30150283239582457068371557_wp
        IM2(3,1)=-0.80901699437494742410229342_wp
        IM2(3,2)= 1.8090169943749474241022934_wp
        IM2(4,1)=-2.3165311563540702775208713_wp
        IM2(4,2)= 3.3165311563540702775208713_wp
        IM2(5,1)=-2.6180339887498948482045868_wp
        IM2(5,2)= 3.6180339887498948482045868_wp

        IM1(1,1)= 1.0_wp
        IM1(2,1)= 0.69849716760417542931628443_wp
        IM1(2,2)= 0.30150283239582457068371557_wp
        IM1(3,1)=-0.80901699437494742410229342_wp
        IM1(3,2)= 1.8090169943749474241022934_wp
        IM1(4,1)=-2.3165311563540702775208713_wp
        IM1(4,2)= 3.3165311563540702775208713_wp
        IM1(5,1)=-2.6180339887498948482045868_wp
        IM1(5,2)= 3.6180339887498948482045868_wp

        IC(1,1)= 1.6180339887498948482045868_wp
        IC(1,2)=-0.61803398874989484820458683_wp
        IC(2,1)= 1.4316949906249123735038224_wp
        IC(2,2)=-0.43169499062491237350382236_wp
        IC(3,1)= 0.5_wp
        IC(3,2)= 0.5_wp
        IC(4,1)=-0.43169499062491237350382236_wp
        IC(4,2)= 1.4316949906249123735038224_wp
        IC(5,1)=-0.61803398874989484820458683_wp
        IC(5,2)= 1.6180339887498948482045868_wp

        IP1(1,1)= 3.6180339887498948482045868_wp
        IP1(1,2)=-2.6180339887498948482045868_wp
        IP1(2,1)= 3.3165311563540702775208713_wp
        IP1(2,2)=-2.3165311563540702775208713_wp
        IP1(3,1)= 1.8090169943749474241022934_wp
        IP1(3,2)=-0.80901699437494742410229342_wp
        IP1(4,1)= 0.30150283239582457068371557_wp
        IP1(4,2)= 0.69849716760417542931628443_wp
        IP1(5,2)= 1.0_wp

        IP2(1,1)= 3.6180339887498948482045868_wp
        IP2(1,2)=-2.6180339887498948482045868_wp
        IP2(2,1)= 3.3165311563540702775208713_wp
        IP2(2,2)=-2.3165311563540702775208713_wp
        IP2(3,1)= 1.8090169943749474241022934_wp
        IP2(3,2)=-0.80901699437494742410229342_wp
        IP2(4,1)= 0.30150283239582457068371557_wp
        IP2(4,2)= 0.69849716760417542931628443_wp
        IP2(5,2)= 1.0_wp

        end subroutine WENO_Coeffs

    pure subroutine Negative_Density_Removal(uin,nq)
      ! this routine calculates the primitive variables
      ! from the conserved variables
      use nsereferencevariables, only: gm1M2, gm1og
      implicit none
      ! number of equations
      integer, intent(in) :: nq
      ! conserved variables
      real(wp), intent(inout) :: uin(nq)
  
      real(wp)                :: vin(nq)
  
      ! density
      vin(1) = abs(uin(1))
      ! velocity
      vin(2:4) = uin(2:4)/vin(1)
      ! temperature
      vin(5) = abs(( uin(5)/vin(1) - gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)))/(1.0_wp-gm1og) )

      ! density
      uin(1) = vin(1)
      ! momentum
      uin(2:4) = vin(1)*vin(2:4)
      ! energy
      uin(5) = vin(1)*( (1.0_wp-gm1og)*vin(5) + gm1M2*0.5_wp*dot_product(vin(2:4),vin(2:4)) )

      return
    end subroutine Negative_Density_Removal

end module SSWENO_Routines
