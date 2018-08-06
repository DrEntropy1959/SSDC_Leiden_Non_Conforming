module SSWENO_Routines

        use precision_vars
        implicit none

        private
        public  WENO_Coeffs
        public  Negative_Density_Removal
        public  Positivity_Preserving_Limiter_Shu

  contains

        subroutine WENO_Coeffs

        use SSWENOvariables

        implicit none

        real(wp), parameter :: sqrt5  = sqrt(5.0_wp)
        real(wp), parameter :: c1     = 15.0_wp/4.0_wp
        real(wp), parameter :: c2     = 15.0_wp/4.0_wp*sqrt5

        real(wp), dimension(5), parameter :: tau_coeffT0 = (/0.0_wp,-c1,c2,-c2,c1/)
        real(wp), dimension(5), parameter :: tau_coeffT1 = (/-c1,c2,-c2,c1,0.0_wp/)

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

!=====================================================================================================

    subroutine Positivity_Preserving_Limiter_Shu(nnodes,ielem,uin,Jac)

      use unary_mod,             only: qsortd
      use referencevariables,    only: nequations
      use nsereferencevariables, only: gm1M2, gamma0, gamI
      use collocationvariables,  only: pvol

      implicit none

      integer,                                intent(in   ) :: nnodes, ielem
      real(wp), dimension(           nnodes), intent(in   ) :: Jac
      real(wp), dimension(nequations,nnodes), intent(inout) :: uin
  
      real(wp), parameter                                   :: tol_r = 1.0e-4_wp
      real(wp), parameter                                   :: tol_p = 1.0e-4_wp

      integer                                               :: i

      real(wp), dimension(nequations,nnodes)                :: vin
      real(wp), dimension(nnodes)                           :: p
      real(wp)                                              :: volume
      real(wp)                                              :: rho_min, P_min, p_ave
      real(wp)                                              :: theta_r,theta_p
      real(wp), dimension(nequations)                       :: u_ave,v_ave
  
      continue

      do i = 1,nnodes
        vin(  1,i) =  uin(  1,i)                              !  conserved variables => primitive variables
        vin(2:4,i) =  uin(2:4,i)/vin(1,i)
        vin(  5,i) = (uin(  5,i)/vin(1,i) - gm1M2*0.5_wp*dot_product(vin(2:4,i),vin(2:4,i)) ) * gamma0
              p(:) =  vin(  1,i)*vin(5,i)                     !  pressure
      enddo

      rho_min = minval(vin(1,:))                              !  minimum value of density
        P_min = minval(  p(  :))                              !  minimum value of pressure

      if((rho_min >= tol_r) .and. (P_min >= tol_p)) return    !  density and pressure positive => return

      u_ave(:) = 0.0_wp ; v_ave(:) = 0.0_wp                   !  Initialize averages of conserved and primitive variableszero
      volume   = 0.0_wp                                       !  Initialize volume
      do i = 1,nnodes                                         !  Total volume and conserved variable in element
        volume   = volume   + pvol(i)*Jac(i)
        u_ave(:) = u_ave(:) + pvol(i)*Jac(i)*uin(:,i)
      enddo
      u_ave(:) = u_ave(:) / volume                            !  normalize by total volume

      v_ave(  1) =  u_ave(  1)
      v_ave(2:4) =  u_ave(2:4)/v_ave(1)
      v_ave(  5) = (u_ave(  5)/v_ave(1) - gm1M2*0.5_wp*dot_product(v_ave(2:4),v_ave(2:4)) ) * gamma0

      p_ave      =  v_ave(1)*v_ave(5)
                                                              !  ==============================================
      if(rho_min <= tol_r) then                               !  begin density correction  ====================
                                                              !  ==============================================
        theta_r = min(1.0_wp , minval( (u_ave(1) - tol_r) / (u_ave(1) - uin(1,:)) ) )

        do i = 1,nnodes
          vin(  1,i) =  u_ave(1) + theta_r*( uin(1,i) - u_ave(1) )!  conserved variables => primitive variables

          vin(2:4,i) =  uin(2:4,i)/vin(1,i)
          vin(  5,i) = (uin(  5,i)/vin(1,i) - gm1M2*0.5_wp*dot_product(vin(2:4,i),vin(2:4,i)) ) * gamma0

            p(    i) =  vin(1,i)*vin(5,i)                       !  new pressure
        enddo

          P_min    = minval(p(:))                               !  minimum value of pressure based on new corrected density
  
        if(P_min >= tol_p) then                                 !  pressure >= tol  => decode conserved variables
        
          do i = 1,nnodes
            uin(  1,i) = vin(1,i)                               ! density

            uin(2:4,i) = vin(1,i)*vin(2:4,i)                    ! momentum

            uin(  5,i) = vin(1,i)*( gamI * vin(5,i) + gm1M2*0.5_wp*dot_product(vin(2:4,i),vin(2:4,i)) ) ! energy
          enddo

          return                                                ! Everything positive => return

        endif
                                                                !  ==============================================
      endif                                                     !  end density correction    ====================
                                                                !  ==============================================

                                                                !  ==============================================
      theta_p = 1.0_wp                                          !  Begin pressure correction ====================
      do i = 1,nnodes                                           !  ==============================================
        if(p(i) <= tol_p) theta_p = min(theta_p, (p_ave - tol_p) / (p_ave - p(i)))
      enddo

      do i = 1,nnodes

        uin(  1,i) = theta_p*(vin(  1,i) - v_ave(  1)) + v_ave(  1)

        uin(2:5,i) = theta_p*(uin(2:5,i) - u_ave(2:5)) + u_ave(2:5)
                                                                !  ==============================================
      enddo                                                     !  End pressure correction ====================
                                                                !  ==============================================

    end subroutine Positivity_Preserving_Limiter_Shu
    
!=====================================================================================================

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
      real(wp), dimension(nnodes)                           :: d_rho, p
      real(wp)                                              :: xnumer, xdenom, xlam, volume
      real(wp)                                              :: rho_min, P_min
  

      vin(1,:) =  uin(1,:)                                  !  conserved variables => primitive variables
      vin(2,:) =  uin(2,:)/vin(1,:)
      vin(3,:) =  uin(3,:)/vin(1,:)
      vin(4,:) =  uin(4,:)/vin(1,:)
      vin(5,:) = (uin(5,:)/vin(1,:)  &
               - gm1M2*0.5_wp*(vin(2,:)*vin(2,:)+vin(3,:)*vin(3,:)+vin(4,:)*vin(4,:)))/(1.0_wp-gm1og)
      p(:)     = vin(1,:)*vin(5,:)      

      rho_min = minval(vin(1,:))                            !  minimum value of density
        P_min = minval(  p(  :))                            !  minimum value of pressure
      if((rho_min >= tol) .and. (P_min >= tol))  return     !  Nothing negative then return

      volume = 0.0_wp                                       !  Total volume of element
      do i = 1,nnodes
        volume = volume + pvol(i)*Jx_r(i,ielem)
      enddo

      if(rho_min < tol) then                                !  First work on the negative density

        call qsortd(uin(1,:),ind,nnodes)                    !  Sort the density into ascending order
        d_rho(:) = 0.0_wp                                   !  Set perturbation to zero

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

      if(P_min < tol) then
        write(*,*)'not finished.  Stopping'
        stop
      endif

    end subroutine Negative_Density_Removal

end module SSWENO_Routines
