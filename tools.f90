module tools
  use precision_vars
  implicit none

  private

  public TDMA

contains

  function distance(x1,x2)
    real(wp), intent(in) :: x1(:), x2(:)

    real(wp) :: distance

    distance = sqrt(dot_product(x2-x1,x2-x1))

  end function distance

  function LS_interpolate(order,npoints,xi,yi,xin)
    ! This subroutine performs a least squares interpolation
    ! to the point xin, given the order, and a set of 
    ! neigboring points.
    integer, intent(in) :: order
    integer, intent(in) :: npoints
    real(wp), intent(in) :: xi(npoints)
    real(wp), intent(in) :: yi(npoints)
    real(wp), intent(in) :: xin

    ! Output:
    real(wp) :: LS_interpolate

    integer :: i
    real(wp) :: ftmp
    real(wp) :: ai(0:order)
    ! First calculate the polynomial coefficients
    call poly_fit(order,npoints,xi,yi,ai)

    ftmp=zero

    ! Now calculate the interpolation value
    do i=order,1,-1
      ftmp=ftmp+ai(order-i)*(xin**(i))  
    end do
    LS_interpolate=ftmp+ai(order)

  end function LS_interpolate

  pure function polynomial(order,ai,xin)
    integer, intent(in) :: order
    real(wp), intent(in) :: ai(0:order)
    real(wp), intent(in) :: xin

    real(wp) :: polynomial

    integer :: i
    real(wp) :: tmp
    tmp=zero

    do i=0,order-1
      tmp=tmp+ai(i)*(xin**(order-i))
    end do
    polynomial=tmp+ai(order)

  end function polynomial

  subroutine poly_fit(order,npoints,xi,yi,ai)
    ! This subroutine performs a polynomial curve fit
    ! of a given order to a set of data. The curve fit
    ! are calculate by solving the system Csq*ai=D
    integer, intent(in) :: order
    integer, intent(in) :: npoints
    real(wp), intent(in) :: xi(npoints)
    real(wp), intent(in) :: yi(npoints)
    real(wp), intent(out) :: ai(order+1)

    ! Local Variables
    real(wp) :: C(npoints,order+1)
    real(wp) :: Ct(order+1,npoints)
    real(wp) :: Csq(order+1,order+1)
    real(wp) :: D(order+1)
    integer :: i,j

    if(npoints<order+1) then
      write(*,*) "poly_fit error: not enough data points for fit"
      stop
    end if

    C=one
    ! Poplulate C
    do i=1,npoints
      do j=0,order-1
        C(i,j+1)=xi(i)**(order-j)
      end do
    end do

    Ct=transpose(C)
    Csq=matmul(Ct,C)
    D=matmul(Ct,yi)

    ! Call gauss elimination
    call gauss_elim(Csq,ai,D,order+1)

  end subroutine poly_fit

  pure subroutine gauss_elim(A,xx,b,n)
    ! This subroutine solves a linear system A*xx=b using
    ! gauss elimination with partial pivoting.
    integer, intent(in) :: n
    real(wp), intent(in) :: A(n,n)
    real(wp), intent(out) :: xx(n)
    real(wp), intent(in) :: b(n)

    ! Local variables
    real(wp) :: AB(n,n), Bb(n)
    integer :: i, k, p
    integer :: ipivot
    real(wp) :: vtemp(n), btemp
    logical :: swap
    real(wp) :: mult

    AB=A
    Bb=b

    vtemp=zero

    ! Forward elimination
    do k=1,n-1
      swap=.false.
      ! partial pivoting
      ipivot=k
      do p=k+1,n
        if(abs(AB(p,k))>abs(AB(ipivot,k)))then
          swap=.true.
          ipivot=p
        end if
      end do
      if(swap)then
        vtemp(:)=AB(k,:)
        AB(k,:)=AB(ipivot,:)
        AB(ipivot,:)=vtemp(:)
        btemp=BB(k)
        BB(k)=BB(ipivot)
        BB(ipivot)=btemp
      end if
      ! end pivoting
      ! eliminate column k
      do i=k+1,n
        mult=AB(i,k)/AB(k,k)
        AB(i,k+1:n)=AB(i,k+1:n)-mult*AB(k,k+1:n)
        Bb(i)=Bb(i)-mult*Bb(k)
      end do
    end do
    ! End Forward elimination

    ! Begin backword substitution
    xx(n)=bb(n)/AB(n,n)
    do i=n-1,1,-1
      xx(i)=(Bb(i)-sum(AB(i,i+1:n)*xx(i+1:n)))/AB(i,i)
    end do

  end subroutine gauss_elim

  function quadsolve(a,b,c)
    ! This function solves a quadratic equation
    ! ax^2 + bx +c = 0, for x. Both possible values
    ! are returned.
    real(wp), intent(in) :: a, b, c
    real(wp) :: quadsolve(2)
    real(wp) :: radicand

    radicand=b*b-four*a*c

    if(radicand<=zero)then
      quadsolve=zero
      return
    end if

    quadsolve(1)=(-b+sqrt(radicand))*half/a
    quadsolve(2)=(-b-sqrt(radicand))*half/a

  end function quadsolve

  pure function dist2d(x1,x2,y1,y2)
    ! This function uses the distance formula in 2 dimensions
    ! to calculate the distance between two points.
    real(wp), intent(in) :: x1, x2
    real(wp), intent(in) :: y1, y2
    real(wp) :: dist2d

    dist2d=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

  end function dist2d

  subroutine addindex(indexptr,ij)
    type(iptr) :: indexptr
    integer :: ij(2)

    integer :: nel(2)
    integer, allocatable :: itmp(:,:)

    if(associated(indexptr%p))then
      nel=size(indexptr%p,DIM=1)
      allocate(itmp(nel(1),2))
      itmp(:,:)=indexptr%p(:,:)
      deallocate(indexptr%p)
      allocate(indexptr%p(1:nel(1)+1,2))
      indexptr%p(1:nel(1),:)=itmp(1:nel(1),:)
      indexptr%p(nel(1)+1,:)=ij(:)
      deallocate(itmp)
    else
      allocate(indexptr%p(1,2))
      indexptr%p(1,:)=ij(:)
    end if

  end subroutine addindex

  !***********************************************************************

  subroutine TDMA(xx,b,L,D,U,n)
    ! This subroutine uses the Thomas Algorithm to solve a tridiagonal
    ! system of equations
    !
    ! Arguments
    ! =========
    !
    ! number of equations to solve
    integer, intent(in) :: n
    real(wp) :: xx(1:n)        ! Solution
    real(wp) :: b(1:n)        ! b vector in Ax=b
    real(wp) :: L(2:n)        ! Lower Vectorization of A
    real(wp) :: D(1:n)        ! Diagonal Vectorization of A
    real(wp) :: U(1:n-1)      ! Upper Vectorization of A
    !
    ! Local Variables
    ! ===============
    !
    ! Loop index 
    integer :: i
    ! Multiplication Factor
    real(wp) :: m

    ! Elementary row operation
    do i=2,n
      m=L(i)/D(i-1)
      D(i)=D(i)-m*U(i-1)
      b(i)=b(i)-m*b(i-1)
    end do

    ! Backward Substitution
    xx(n)=b(n)/D(n)
    do i=n-1,1,-1
      xx(i)=(b(i)-U(i)*xx(i+1))/D(i)
    end do

  end subroutine TDMA

  subroutine tdma_per(xx,dt,at,bt,ct,n)
    ! This subroutine uses the Thomas Algorithm to solve a tridiagonal
    ! system of equations for a periodic system
    !
    ! Arguments
    ! =========
    !
    ! number of equations to solve
    integer, intent(in) :: n
    real(wp) :: xx(1:n)       ! Solution
    real(wp) :: dt(1:n)        ! b vector in Ax=b
    real(wp) :: at(1:n)        ! Lower Vectorization of A
    real(wp) :: bt(1:n)        ! Diagonal Vectorization of A
    real(wp) :: ct(1:n)        ! Upper Vectorization of A
    !
    ! Local Variables
    ! ===============
    !
    ! Loop index 
    integer :: i
    ! Multiplication Factor
    real(wp) :: m
    real(wp) :: x1(1:n)
    real(wp) :: x2(1:n)
    real(wp) :: xstar
    real(wp) :: vt(1:n)

    vt=zero
    xstar=zero
    x1=zero; x2=zero;

    vt(1)=at(1)
    vt(n)=ct(n)
    bt(1)=bt(1)-at(1)
    bt(n)=bt(n)-ct(n)

    ! Forward Elimination
    do i=2,n
      m=at(i)/bt(i-1)
      bt(i)=bt(i)-m*ct(i-1)
      dt(i)=dt(i)-m*dt(i-1)
      vt(i)=vt(i)-m*vt(i-1)
    end do

    ! Backward Substitutions
    ! The two matrix systems are M* X1 = dt  and M* X2 = vt
    x1(n)=dt(n)/bt(n)
    x2(n)=vt(n)/bt(n)
    do i=n-1,1,-1
      x1(i)=(dt(i)-ct(i)*x1(i+1))/bt(i)
      x2(i)=(vt(i)-ct(i)*x2(i+1))/bt(i)
    end do

    !Calculate supplementary unknown
    xstar=(x1(1)+x1(n))/(one+x2(1)+x2(n))

    forall(i=1:n)
      xx(i)=x1(i)-x2(i)*xstar
    end forall

  end subroutine tdma_per

  !***********************************************************************

  subroutine PENT(xx,b,LL,L,D,U,UU,n)
    ! This subroutine solves a linear system Ax=b, where A is a  
    ! Pentadiagonal matrix. Diagonal Dominance must be satisfied
    ! in the matrix in order for the solver to give an accurate
    ! solution. There is  no checking algorithm, an incorrect
    ! answer can be given for a system where diagonal dominance
    ! is not satisfied.

    integer, intent(in) :: n
    real(wp) :: xx(1:n)           ! Solution
    real(wp) :: b(1:n)           ! b vector in Ax=b
    real(wp) :: LL(3:n)          ! Second Lower Vectorization of A
    real(wp) :: L(2:n)           ! Lower Vectorization of A
    real(wp) :: D(1:n)           ! Diagonal Vectorization of D
    real(wp) :: U(1:n-1)           ! Upper Vectorization of A
    real(wp) :: UU(1:n-2)        ! Second Upper Vectorization of A

    integer :: i
    real(wp) :: m                ! multiplication factor

    ! Peform Gauss Elimination for pentediagonal Matrix
    ! Start by setting the LL vector to zero.
    do i=3,n-1
      m=LL(i)/L(i-1)
      L(i)=L(i)-m*D(i-1)
      D(i)=D(i)-m*U(i-1)
      U(i)=U(i)-m*UU(i-1)
      b(i)=b(i)-m*b(i-1)
    end do
    i=n
    m=LL(i)/L(i-1)
    L(i)=L(i)-m*D(i-1)
    D(i)=D(i)-m*U(i-1)
    b(i)=b(i)-m*b(i-1)
    do i=2,n-1
      m=L(i)/D(i-1)
      D(i)=D(i)-m*U(i-1)
      U(i)=U(i)-m*UU(i-1)
      b(i)=b(i)-m*b(i-1)
    end do
    i=n
    m=L(i)/D(i-1)
    D(i)=D(i)-m*U(i-1)
    b(i)=b(i)-m*b(i-1)

    ! Backward Substitution
    xx(n)=b(n)/D(n)
    xx(n-1)=(b(n-1)-U(n-1)*xx(n))/D(n-1)
    do i=n-2,1,-1
      xx(i)=(b(i)-U(i)*xx(i+1)-UU(i)*xx(i+2))/D(i)
    end do

  end subroutine PENT

end module tools
