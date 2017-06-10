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

    subroutine Negative_Density_Removal(nnodes,ielem,uin)

      use unary_mod,             only: qsortd
      use referencevariables,    only: nequations
      use nsereferencevariables, only: gm1M2, gm1og
      use collocationvariables,  only: pvol
      use variables,             only: Jx_r

      implicit none

      integer,                                intent(in   ) :: nnodes, ielem
      real(wp), dimension(nequations,nnodes), intent(inout) :: uin
  
      real(wp), parameter                                   :: tol = 1.0e-10_wp

      integer,  dimension(nnodes)                           :: ind
      integer                                               :: i

      real(wp), dimension(nequations,nnodes)                :: vin
      real(wp), dimension(nnodes)                           :: d_rho
      real(wp)                                              :: xnumer, xdenom, xlam, volume
      real(wp)                                              :: rho_min, T_min
  

      vin(1,:) =  uin(1,:)
      vin(2,:) =  uin(2,:)/vin(1,:)
      vin(3,:) =  uin(3,:)/vin(1,:)
      vin(4,:) =  uin(4,:)/vin(1,:)
      vin(5,:) = (uin(5,:)/vin(1,:)  &
               - gm1M2*0.5_wp*(vin(2,:)*vin(2,:)+vin(3,:)*vin(3,:)+vin(4,:)*vin(4,:)))/(1.0_wp-gm1og)

      rho_min = minval(vin(1,:))
        T_min = minval(vin(5,:))
      if((rho_min >= tol) .and. (T_min >= tol))  return

      volume = 0.0_wp
      do i = 1,nnodes
        volume = volume + pvol(i)*Jx_r(i,ielem)
      enddo

      if(rho_min < tol) then

        call qsortd(uin(1,:),ind,nnodes)

        xnumer = 0.0_wp ; xdenom = volume ;
        do i = 1,nnodes-1
          d_rho(ind(i)) = - (uin(1,ind(i)) - tol)
          xnumer = xnumer + (uin(1,ind(i)) - tol) * pvol(ind(i)) * Jx_r(ind(i),ielem)
          xdenom = xdenom -                         pvol(ind(i)) * Jx_r(ind(i),ielem)
          xlam   = xnumer / xdenom
          if( uin(1,ind(i+1)) + xlam >= tol) then
            uin(1,ind(i+1:nnodes)) = xlam
            exit
          endif
        enddo
      
        uin(1,:) = uin(1,:) + d_rho(:)
      endif

      if(T_min < tol) then
        write(*,*)'not finished.  Stopping'
        stop
      endif


      return
    end subroutine Negative_Density_Removal

end module SSWENO_Routines
