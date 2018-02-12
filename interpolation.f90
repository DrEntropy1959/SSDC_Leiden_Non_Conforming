module interpolation

  use precision_vars
  use tools

  implicit none

  private

  public Poly3_Intrp
  public Extrp_XiEtaZeta_neq
  public WENO_Mapping_Coefs
  public WENO_xi_val
  public xi_guess

contains 


  function useSpline2(ain,rvec,rin,npts,ndim)
    ! this function calculates the interpolated value from a cubic spline
    ! with given coefficients.
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rvec(1:npts) 
    real(wp), intent(in) :: rin
    real(wp), intent(in) :: ain(4,ndim,npts-1) ! coefficients for splines

    real(wp) :: useSpline2(ndim)


    useSpline2 = useSpline(ain,rvec,rin,npts,ndim)

  end function useSpline2

  function useSpline2_s(ain,rvec,rin,npts,ndim)
    ! this function calculates the interpolated value from a cubic spline
    ! with given coefficients.
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rvec(1:npts) 
    real(wp), intent(in) :: rin
    real(wp), intent(in) :: ain(4,ndim,npts-1) ! coefficients for splines

    real(wp) :: useSpline2_s(ndim)


    useSpline2_s = useSpline_s(ain,rvec,rin,npts,ndim)

    !     if (rin > rvec(npts)) then
    !       useSpline2_s = 0.0_wp
    !     end if

  end function useSpline2_s

  function useSpline2_s2(ain,rvec,rin,npts,ndim)
    ! this function calculates the interpolated value from a cubic spline
    ! with given coefficients.
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rvec(1:npts) 
    real(wp), intent(in) :: rin
    real(wp), intent(in) :: ain(4,ndim,npts-1) ! coefficients for splines

    real(wp) :: useSpline2_s2(ndim)

    useSpline2_s2 = useSpline_s2(ain,rvec,rin,npts,ndim)

    !     if (rin > rvec(npts)) then
    !       useSpline2_s2 = 0.0_wp
    !     end if

  end function useSpline2_s2

  function useSpline(ain,rvec,rin,npts,ndim)
    ! this function calculates the interpolated value from a cubic spline
    ! with given coefficients.
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rvec(1:npts) 
    real(wp), intent(in) :: rin
    real(wp), intent(in) :: ain(4,ndim,npts-1) ! coefficients for splines

    real(wp) :: useSpline(ndim)

    real(wp) :: rx

    real(wp) :: dpp(npts)
    integer :: imin(1)
    integer :: dir

    useSpline = 0.0_wp

    rx = rin
    rx = min(rin,rvec(npts))
    rx = max(rx,rvec(1))
    dpp = rx-rvec
    imin = min(minloc(dpp,dpp>=zero),npts-1)
    if(imin(1)<1)then
      write(*,*) "error in spline fit", imin(1)
      stop
    end if
    do dir=1,ndim
      useSpline(dir) = ain(1,dir,imin(1))+ain(2,dir,imin(1))*rx &
        & + ain(3,dir,imin(1))*rx*rx + ain(4,dir,imin(1))*rx*rx*rx
    end do

  end function useSpline

  function useSpline_s(ain,rvec,rin,npts,ndim)
    ! this function calculates the interpolated value from a cubic spline
    ! with given coefficients.
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rvec(1:npts) 
    real(wp), intent(in) :: rin
    real(wp), intent(in) :: ain(4,ndim,npts-1) ! coefficients for splines

    real(wp) :: useSpline_s(ndim)

    real(wp) :: rx

    real(wp) :: dpp(npts)
    integer :: imin(1)
    integer :: dir

    useSpline_s = 0.0_wp

    rx = rin
    rx = min(rin,rvec(npts))
    rx = max(rx,rvec(1))
    dpp = rin-rvec
    imin = min(minloc(dpp,dpp>=zero),npts-1)
    if (imin(1) < 1) then
      write(*,*) "error in spline fit", imin(1)
      stop
    end if
    do dir = 1,ndim
      useSpline_s(dir) = ain(2,dir,imin(1))+two*ain(3,dir,imin(1))*rx &
        & + three*ain(4,dir,imin(1))*rx*rx
    end do

  end function useSpline_s

  function useSpline_s2(ain,rvec,rin,npts,ndim)
    ! this function calculates the interpolated value from a cubic spline
    ! with given coefficients.
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rvec(1:npts) 
    real(wp), intent(in) :: rin
    real(wp), intent(in) :: ain(4,ndim,npts-1) ! coefficients for splines

    real(wp) :: useSpline_s2(ndim)

    real(wp) :: rx

    real(wp) :: dpp(npts)
    integer :: imin(1)
    integer :: dir

    useSpline_s2 = 0.0_wp

    rx = rin
    rx = min(rin,rvec(npts))
    rx = max(rx,rvec(1))
    dpp = rin-rvec
    imin = min(minloc(dpp,dpp>=zero),npts-1)
    if (imin(1) < 1) then
      write(*,*) "error in spline fit", imin(1)
      stop
    end if
    do dir = 1,ndim
      useSpline_s2(dir) = 2.0_wp*ain(3,dir,imin(1)) &
        & + 6.0_wp*ain(4,dir,imin(1))*rx
    end do

  end function useSpline_s2

  subroutine splinefit(aout,rin,xin,npts,ndim)
    ! this subroutine fits cubic splines through
    ! a vector function as a function of parameter r
    integer, intent(in) :: npts, ndim
    real(wp), intent(in) :: rin(npts) ! parameter function
    real(wp), intent(in) :: xin(ndim,npts) ! vector values
    real(wp), intent(out) :: aout(4,ndim,npts-1) ! coefficients for splines

    integer :: dir, i
    real(wp) :: b(npts-2), DD(npts-2), LL(npts-3), UU(npts-3)
    real(wp) :: dx(npts-1)
    real(wp) :: ypp(npts)

    do i=1,npts-1
      dx(i) = rin(i+1)-rin(i)
    end do
    ! set tridiagonal matrix (only a function of grid spacing)
    i=1 
    !     DD(i) = three*half*dx(i)+two*dx(i+1)
    !     UU(i) = dx(i+1)
    DD(i) = (dx(i)+dx(i+1))*(dx(i)+two*dx(i+1))/dx(i+1)
    UU(i) = (dx(i+1)**2-dx(i)**2)/dx(i+1)
    do i=2,npts-3
      LL(i-1) = dx(i)
      DD(i) = two*(dx(i)+dx(i+1))
      UU(i) = dx(i+1)
    end do
    i=npts-2
    !     LL(i-1) = dx(i)
    !     DD(i) = two*dx(i) + three*half*dx(i+1)
    LL(i-1) = (dx(i)**2-dx(i+1)**2)/dx(i)
    DD(i) = (dx(i)+dx(i+1))*(two*dx(i)+dx(i+1))/dx(i)
    ! loop over dimensions for rhs and solve
    do dir=1,ndim
      ! set rhs
      do i=1,npts-2
        b(i) = six*((xin(dir,i+2)-xin(dir,i+1))/dx(i+1) &
          & - (xin(dir,i+1)-xin(dir,i))/dx(i))
      end do
      ! solve equation for y''
      ypp = 0.0_wp
      call TDMA(ypp(2:npts-1),b,LL,DD,UU,npts-2)
      ypp(1)=ypp(2)*(dx(1)+dx(2))/dx(2)-ypp(3)*dx(1)/dx(2)
      ypp(npts)=ypp(npts-1)*(dx(npts-2)+dx(npts-1))/dx(npts-2)-ypp(npts-2)*dx(npts-1)/dx(npts-2)
      ! calculate coefficients
      do i=1,npts-1
        ! order 0 term
        aout(1,dir,i) = (dx(i)*rin(i)*ypp(i+1) - dx(i)*rin(i+1)*ypp(i) &
          & - rin(i)**3*ypp(i+1)/dx(i) + rin(i+1)**3*ypp(i)/dx(i))/six &
          & - rin(i)*xin(dir,i+1)/dx(i) + rin(i+1)*xin(dir,i)/dx(i)
        ! order 1 term
        aout(2,dir,i) = dx(i)/six*(ypp(i)-ypp(i+1)) &
          & + half/dx(i)*(rin(i)*rin(i)*ypp(i+1)-rin(i+1)*rin(i+1)*ypp(i)) &
          & + (xin(dir,i+1)-xin(dir,i))/dx(i)
        ! order 2 term
        aout(3,dir,i) = half/dx(i)*(rin(i+1)*ypp(i)-rin(i)*ypp(i+1))
        ! order 3 term
        aout(4,dir,i) = (ypp(i+1)-ypp(i))/(six*dx(i))
      end do
    end do

  end subroutine splinefit

!============================================================================

      function Poly3_Intrp(nq,u,xL,xR,x)

        implicit none
        integer,                    intent(in)  :: nq
        real(wp),                   intent(in)  :: xL,xR,x
        real(wp), dimension(nq,4),  intent(in)  :: u

        real(wp), dimension(nq)                 :: Poly3_Intrp
        real(wp)                                :: dx,dx2,dx3
        real(wp)                                :: sqrt5

        dx  = xR - xL
        dx2 = dx*dx
        dx3 = dx2*dx

        sqrt5 = sqrt(5.0_wp)

!       Poly3_Intrp(:) =  (2.0_wp*dx3*u(:,1) + dx2*(-12.0_wp*u(:,1) + 5.0_wp*(1.0_wp + sqrt5)*u(:,2) + 5.0_wp*u(:,3)  &
!                      & - 5.0_wp*sqrt5*u(:,3) + 2.0_wp*u(:,4))*(x - xL) - 5.0_wp*dx*(-4.0_wp*u(:,1) + u(:,2) +       &
!                      &   3.0_wp*sqrt5*u(:,2) + u(:,3) - 3.0_wp*sqrt5*u(:,3) + 2.0_wp*u(:,4))*(x - xL)*(x - xL)     &
!                      & -10.0_wp*(u(:,1) - sqrt5*u(:,2) + sqrt5*u(:,3) - u(:,4))*(x-xL)*(x-xL)*(x-xL))/(2.0_wp*dx3)

        Poly3_Intrp(:) = - (  (-1 +       x)*(-1 + 5*x**2)*u(:,1))/8.0_wp &
                         + (5*(-1 + sqrt5*x)*(-1 +   x**2)*u(:,2))/8.0_wp &
                         - (5*(+1 + sqrt5*x)*(-1 +   x**2)*u(:,3))/8.0_wp &
                         + (  (+1 +       x)*(-1 + 5*x**2)*u(:,4))/8.0_wp

      end function

  subroutine Extrp_XiEtaZeta_neq(Ndim,neq,NPts,xi,FA,F)
  !  
  ! Extrapolate data in the computational space 
  ! Assuses existing tensor product grid and an arbitrary location in (xi,eta,zeta) stace
  ! Assume tensor product distributions.  i.e. the same in each direction
  !  
  
  ! Load modules

  ! Nothing is implicitly defined
  implicit none

  integer,                    intent(in   )  :: Ndim, neq, NPts
  real(wp), dimension(Ndim),  intent(in   )  :: xi
  real(wp), dimension(:,:),   intent(in   )  :: FA
  real(wp), dimension(:),     intent(inout)  :: F

  real(wp), dimension(neq,NPts,NPts)         :: F1
  real(wp), dimension(neq,NPts)              :: F2

  real(wp), dimension(neq,NPts)              :: u

  integer                                    :: j,k,m,n
  integer                                    :: StrideY, StrideZ

  select case(ndim)

  case(2)

    ! Extrapolate in the xi direction;
!   StrideY = NPts
!   F1(:,:,:,:) = 0.0_wp
!   do j = 1,NPts
!     do i = 1,NPts
!       do m = 1,NPts
!         n = + (j-1)*StrideY + m
!         u(:,m) = FA(:,n)
!       enddo
!       F1(:,i,j,1) = F1(:,i,j,1) + Poly3_Intrp(neq,u,-1.0_wp,1.0_wp,xi(1))
!     enddo
!   enddo

    ! Extrapolate in the eta direction;
!   StrideY = NPts
!   F(:) = 0.0_wp
!   do j = 1,NPts
!     do i = 1,NPts
!       do m = 1,NPts
!         u(:,m) = F1(:,i,m,1)
!       enddo
!       F(:)        = F(:)        + Poly3_Intrp(neq,u,-1.0_wp,1.0_wp,xi(2))
!     enddo
!   enddo

  case(3)

    !  Poly3_Intrp(3,  u  ,-1.0_wp,+1.0_wp,1/sqrt5)
    !  Poly3_Intrp(nq,data, xL    ,   xR  ,  xx   )

    StrideY = NPts
    strideZ = NPts * NPts

    F1(:,:,:) = 0.0_wp
    do k = 1,NPts
      do j = 1,NPts
        do m = 1,NPts
          n = (k-1)*strideZ + (j-1)*StrideY + m
          u(:,m) = FA(:,n)
        enddo
        F1(:,j,k) = F1(:,j,k) + Poly3_Intrp(neq,u,-1.0_wp,1.0_wp,xi(1))
      enddo
    enddo

    F2(:,:) = 0.0_wp
    do k = 1,NPts
      do m = 1,NPts
        u(:,m) = F1(:,m,k)
      enddo
      F2(:,k) = F2(:,k) + Poly3_Intrp(neq,u,-1.0_wp,1.0_wp,xi(2))
    enddo

    F(:) = Poly3_Intrp(neq,F2,-1.0_wp,1.0_wp,xi(3))

  case default

    write(*,*)'NDim must be either 2 or three'
    write(*,*)'stopping'
    stop
  end select

  return

  end subroutine Extrp_XiEtaZeta_neq

  subroutine WENO_Mapping_Coefs(xg,nodesperelem,comp2phys_coeffs)

    integer, parameter  :: wp = 8 
    integer,                             intent(in)     :: nodesperelem

    real(wp), dimension(3,nodesperelem), intent(in)     :: xg
    real(wp), dimension(3,8),            intent(inout)  :: comp2phys_coeffs

    real(wp), dimension(8,8)                 :: Jac, vec_test, eye
    real(wp), dimension(3,8)                 :: xtrain
    real(wp)                                 :: t1
    integer                                  :: i

!   real(wp),  dimension(8,8)                :: JacT= reshape(                    &
!                                                   & (/+1,-1,-1,-1,+1,+1,+1,-1,  &
!                                                   &   +0,+1,+0,+0,-1,-1,+0,+1,  &
!                                                   &   +0,+0,+0,+0,+1,+0,+0,-1,  &
!                                                   &   +0,+0,+1,+0,-1,+0,-1,+1,  &
!                                                   &   +0,+0,+0,+1,+0,-1,-1,+1,  &
!                                                   &   +0,+0,+0,+0,+0,+1,+0,-1,  &
!                                                   &   +1,+0,+0,+0,+0,+0,+0,+1,  &
!                                                   &   +0,+0,+0,+0,+0,+0,+1,-1/),&
!                                                   &      (/8,8/) ) 

!    integer,  dimension(8),   parameter     :: corners  = reshape(                    &
!                                                        & (/1,4,16,13,49,52,64,61/),  &
!                                                        &   (/8/) )

!    real(wp),  dimension(8,8), parameter     :: vec_testT= reshape(                    &
!                                                        & (/+1,+1,+1,+1,+1,+1,+1,+1,  &
!                                                        &   -1,+1,+1,-1,-1,+1,+1,-1,  &
!                                                        &   -1,-1,+1,+1,-1,-1,+1,+1,  &
!                                                        &   -1,-1,-1,-1,+1,+1,+1,+1,  &
!                                                        &   -1,-1,+1,-1,-1,-1,+1,-1,  &
!                                                        &   -1,-1,-1,-1,-1,+1,+1,-1,  &
!                                                        &   -1,-1,-1,-1,-1,-1,+1,+1,  &
!                                                        &   -1,-1,-1,-1,-1,-1,+1,-1/),&
!                                                        &      (/8,8/) )

!    Jac      = Transpose(JacT) / 2.0_wp
!    vec_test = Transpose(vec_testT)

    real(wp),  dimension(8,8), parameter :: vec_testT= reshape(                    &
                                                     & (/+1,+1,+1,+1,+1,+1,+1,+1,  &
                                                     &   -1,+1,-1,+1,-1,+1,-1,+1,  &
                                                     &   -1,-1,+1,+1,-1,-1,+1,+1,  &
                                                     &   -1,-1,-1,-1,+1,+1,+1,+1,  &
                                                     &   +1,-1,-1,+1,+1,-1,-1,+1,  &
                                                     &   +1,-1,+1,-1,-1,+1,-1,+1,  &
                                                     &   +1,+1,-1,-1,-1,-1,+1,+1,  &
                                                     &   -1,+1,+1,-1,+1,-1,-1,+1/),&
                                                     &      (/8,8/) ) 

     integer,  dimension(8),   parameter     :: corners  = reshape(                    &
                                                         & (/1,4,13,16,49,52,61,64/),  &
                                                         &   (/8/) )

!    Jac      = Transpose(JacT) 
!    vec_test = Transpose(Jac ) / 8 

     eye(:,:) = 0.0_wp
     do i = 1,8
       eye(i,i) = 1.0_wp
     enddo

     vec_test = Transpose(vec_testT)
     Jac      = Transpose(vec_test) / 8.0_wp


     do i = 1,8
       xtrain(:,i) = xg(:,corners(i))
     enddo
     

     comp2phys_coeffs = matmul(xtrain,Jac) ;
     t1 = maxval(abs(xtrain - matmul(comp2phys_coeffs,vec_test)))
     if(t1 >= 1.0e-12) then
       write(*,*)'eye error ',maxval(matmul(Jac,Vec_test)-eye(:,:))
       write(*,*)'mapping error',t1
       write(*,*)'somethings wrong in the computational 2 physical mapping routine'
       write(*,*)'WENO_Mapping_Coefs'
       write(*,*)'stopping'
       stop
     endif

  end subroutine WENO_Mapping_Coefs

  subroutine WENO_xi_val(comp2phys_coeffs,xi,x,rnorm)

    !  xi enters with the initial estimate of the computational variables (xi,et,zt)
    !  x  is the target location in physical space  (x0,y0,z0)
    !  Routine implements modified Newtons method with line search

    integer, parameter  :: wp = 8 

    real(wp), dimension(3,8), intent(inout)  :: comp2phys_coeffs
    real(wp), dimension(3),   intent(inout)  :: xi, x
    real(wp),                 intent(  out)  :: rnorm

    real(wp), dimension(3)                   :: res,xiNew, dxi, xi0
    real(wp), dimension(3,3)                 :: WENO_Jac, WENO_JacI

    real(wp)                                 :: al, rnormt
    integer                                  :: i,j,k, ierr


     xi0(:) = xi(:)

     call WENOIntrpFunc(comp2phys_coeffs,xi,x,res)
     rnorm  = sqrt(dot_product(res(:),res(:)))
     ierr   = 1

     do k = 1,100
        xi(:) = xi0(:) + 0.1*(k-1)*rand()
       do i = 1,100
  
         call WENOIntrpFunc(comp2phys_coeffs,xi,x,res)
         rnorm = sqrt(dot_product(res(:),res(:)))
  
         call WENO_Intrp_Jac(comp2phys_coeffs,xi,WENO_Jac)
         call Inv_3x3(WENO_Jac,WENO_JacI)
  
         dxi = MatMul(WENO_JacI,res)
  
          al = 1.0_wp
         do j = 1,20    !   under-relax the value of the parameter alpha
  
           xiNew(:) = xi(:) - al*dxi
  
             call WENOIntrpFunc(comp2phys_coeffs,xiNew,x,res)
  
             rnormt = sqrt(dot_product(res(:),res(:)))
  
           if(rnormt >= rnorm) then
               al = 0.5_wp * al
           else
             exit
           endif
  
         enddo
  
         xi(:) = xi(:) - al*dxi
  
         rnorm = rnormt  
         
!        write(*,*)k,i,j,rnorm
         if(rnorm <= 1.0e-13_wp) then
           ierr = 0
           return
         endif
        
       enddo
     enddo
     if(ierr == 1) write(*,*)'did not converge to computational point xi,et,ze for ', x

  end subroutine WENO_xi_val

  subroutine WENOIntrpFunc(comp2phys_coeffs,xi,x,residual)

    !  xi is the current iterate of the computational variables (xi,et,zt)
    !  x  is the target location in physical space  (x0,y0,z0)

    integer, parameter  :: wp = 8 

    real(wp), dimension(3,8), intent(inout)  :: comp2phys_coeffs
    real(wp), dimension(3),   intent(inout)  :: xi,x

    real(wp), dimension(3),   intent(  out)  :: residual

    real(wp), dimension(8)  :: xivec 

    xivec(1) = 1.0_wp
    xivec(2) = xi(1)
    xivec(3) =       xi(2)
    xivec(4) =             xi(3)
    xivec(5) = xi(1)*xi(2)
    xivec(6) = xi(1)*      xi(3)
    xivec(7) =       xi(2)*xi(3)
    xivec(8) = xi(1)*xi(2)*xi(3)
    
    residual(:) = Matmul(comp2phys_coeffs,xivec) - x(:)

  end subroutine WENOIntrpFunc

  subroutine WENO_Intrp_Jac(comp2phys_coeffs,xi,Jac)

    integer, parameter  :: wp = 8 

    real(wp), dimension(3,8), intent(in)     :: comp2phys_coeffs
    real(wp), dimension(3),   intent(in)     :: xi
    real(wp), dimension(3,3), intent(out)    :: Jac

    real(wp), dimension(3,8)                 :: xi_derivativesT


        xi_derivativesT=    reshape(  (/0.0_wp,0.0_wp,0.0_wp,   &
                                      & 1.0_wp,0.0_wp,0.0_wp,   &
                                      & 0.0_wp,1.0_wp,0.0_wp,   &
                                      & 0.0_wp,0.0_wp,1.0_wp,   &
                                      & xi(2) ,xi(1) ,0.0_wp,   &
                                      & xi(3) ,0.0_wp,xi(1) ,   &
                                      & 0.0_wp,xi(3) ,xi(2) ,   &
                        xi(2)*xi(3),xi(1)*xi(3),xi(1)*xi(2)  /),&
                                      &      (/3,8/) )

       Jac = Matmul(comp2phys_coeffs,transpose(xi_derivativesT))

  end subroutine WENO_Intrp_Jac

  subroutine Inv_3x3(J,JI)

    integer, parameter  :: wp = 8 

    real(wp), dimension(3,3), intent(in)     :: J
    real(wp), dimension(3,3), intent(out)    :: JI
    real(wp), dimension(3,3)                 :: JIT

    real(wp)                                 :: det

       det = - J(1,3)*J(2,2)*J(3,1)  + J(1,2)*J(2,3)*J(3,1)  &
             + J(1,3)*J(2,1)*J(3,2)  - J(1,1)*J(2,3)*J(3,2)  &
             - J(1,2)*J(2,1)*J(3,3)  + J(1,1)*J(2,2)*J(3,3)

        JIT = reshape(                                                                                 &
    (/ -J(2,3)*J(3,2) + J(2,2)*J(3,3),+J(1,3)*J(3,2) - J(1,2)*J(3,3),-J(1,3)*J(2,2) + J(1,2)*J(2,3),   &
       +J(2,3)*J(3,1) - J(2,1)*J(3,3),-J(1,3)*J(3,1) + J(1,1)*J(3,3),+J(1,3)*J(2,1) - J(1,1)*J(2,3),   &
       -J(2,2)*J(3,1) + J(2,1)*J(3,2),+J(1,2)*J(3,1) - J(1,1)*J(3,2),-J(1,2)*J(2,1) + J(1,1)*J(2,2) /),&
        (/3,3/))

       if(abs(det) <= 1.0e-15) then
         write(*,*)'mapping jacobian is zero:  too close to a corner'
        write(*,*)'stopping'
         stop
       endif

       JI = transpose(JIT) / det

  end subroutine Inv_3x3

  !============================================================================
  !  Find the 8 closest neighbors to a point x0.
  !  Assumes within element, not adjoining.
  !============================================================================
  
  function xi_guess(xg,x0,nodesperelem)

    use unary_mod,      only: qsortd
    use precision_vars, only: magnitude

    integer, parameter  :: wp = 8 

    real(wp), dimension(3,nodesperelem), intent(in   ) :: xg
    real(wp), dimension(3),              intent(in   ) :: x0
    integer,                             intent(in   ) :: nodesperelem

    integer,  dimension(nodesperelem)                  :: ind
!   integer,  dimension(8)                             :: neighbors

    real(wp), dimension(nodesperelem)                  :: distance

    real(wp), parameter                                :: sqrt5I = 1.0_wp/sqrt(5.0_wp)
    real(wp), dimension(4)                 :: xi

    integer                                            :: i ,j ,k
    integer                                            :: i0,j0,k0
    integer, parameter                                 :: ns = 4

    real(wp), dimension(3)                             :: xi_guess

    xi = reshape((/-1.0_wp,-sqrt5I,+sqrt5I,+1.0_wp/),(/4/))

    do i = 1,nodesperelem
      distance(i) = magnitude(x0(:)-xg(:,i))
    enddo
      
    call qsortd(distance,ind,nodesperelem)
    
    i0 = 0 ; j0 = 0 ; k0 = 0 ;
    do k = 1,ns
      do j = 1,ns
        do i = 1,ns
          if(ind(1) == (k-1)*ns*ns + (j-1)*ns + i) then
            i0 = i ; j0 = j ; k0 = k
          endif
        enddo
      enddo
    enddo
    
    xi_guess(1) = xi(i0)
    xi_guess(2) = xi(j0)
    xi_guess(3) = xi(k0)

  end function xi_guess

end module interpolation
