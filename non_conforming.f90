module non_conforming

  use precision_vars
  use initcollocation

  !-- Nothing is implicitely defined
  implicit none

  !-- private subroutines functions etc
  private

  !-- public subroutines functions etc
  public Lagrange_interpolant_basis_1D
  public Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo
  public Vandermonde_1D_monomial
  public Rotate_xione_2_xitwo_and_back
contains
!==================================================================================================
!
! Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo()
!
! Purpose: Constructs the Vandermonde matrix from the Lagrange basis functions on the nodal set 
!          xione evaluated at the nodes of xitwo
! 
! Comments: The computation of the interpolation basis functions is based on 
!           "Barycentric Lagrange Interpolation", SIAM REview, Vol 46, No 3 pp. 501-517 
!           (see docs/barycentric.pdf). The notation in this subroutine is consistent with that 
!           in the paper and specifically with the section 3 "An improved Lagrange Formula".
!
! Additional documentation: docs/barycentric.pdf
!
! Unit tests: unit_tests/Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo_program 
!             compile using sh compile_Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo
!
! Inputs: nxione (int): number of nodes in the nodal set xione
!         nxitwo (int): number of nodes in the nodal set xitwo
!         xione (real(wp)) size (nxioneX1): nodal set that is used to construct the Lagrange basis functions
!         xitwo (real(wp)) size (nxitwoX1): nodal set at which the Lagrange basis functions will be evaluated at
!
! Outputs:
!         Vandermonde (real(wp)) size (nxitwoXnxione): the above mentioned Vandermonde matrix
!
!===================================================================================================
   subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo(nxione,nxitwo,xione,xitwo,Vandermonde)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nxione, nxitwo
     real(wp), intent(in)                        :: xione(nxione), xitwo(nxitwo)
     real(wp), intent(inout)                     :: Vandermonde(nxitwo,nxione)
     
     !-- local variables
     integer                                     :: i, j, k
     real(wp)                                    :: wj(nxione),l
 
     !-- construct the bary centric weights, wj,
     do j = 0,nxione-1 
       wj(j+1) = 1.0_wp
       do k = 0,nxione-1
         if (k.NE.j) then
           wj(j+1) = wj(j+1)*(xione(j+1)-xione(k+1))
         end if
       end do
       wj(j+1) = 1.0_wp/wj(j+1)
     end do
     !-- construct the Vandermonde matrix V(j,k) = l_k(xitwo(j)), where l_k is the Lagrange 
     !-- interpolation basis construted from xione
     do j = 1,nxione
       do k = 1,nxitwo
         !-- construct l
         l = 1.0_wp
         do i = 1, nxione
           l = l*(xitwo(j)-xione(i))
         end do
         !-- check to see if the nodal location on the second set of nodes matches the node 
         !-- location associated with the kth Lagrange basis function 
         if (abs(xitwo(j)-xione(k)) <=2*epsilon(1.0_wp)) then
           Vandermonde(j,k) = 1.0_wp
         else
           Vandermonde(j,k)  = l*wj(k)/(xitwo(j)-xione(k))
         end if
       end do
     end do

   end subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo

!==================================================================================================
!
! Vandermonde_1D_monomial()
!
! Purpose: Constructs the Vandermonde matrix from the monomials on the node set xi 
! 
! Comments: 
!
! Additional documentation:
!
! Unit tests: unit_tests/Vandermonde_1D_monomial.f90 
!             compile using sh unit_tests/comiplie_Vandermonde_1D_monomial.sh
!
! Inputs: n (int): number of nodes in the nodal set xi
!         p (int): the highest degree monomial that will be evaluated
!         xi (real(wp)) size (nX1): nodal on which the monomials will be evaluated
!
! Outputs:
!         Vandermonde (real(wp)) size (nX(p+1)): the above mentioned Vandermonde matrix
!
!===================================================================================================
   subroutine Vandermonde_1D_monomial(n,p,xi,Vandermonde)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: n, p
     real(wp), intent(in)                        :: xi(n)
     real(wp), intent(inout)                     :: Vandermonde(n,(p+1))
     
     !-- local variables
     integer                                     :: i, j
 
     !-- construct the Vandermonde matrix
     do i = 1,n
       do j = 1,p+1
         Vandermonde(i,j) = xi(i)**(j-1)
       end do
     end do
   end subroutine Vandermonde_1D_monomial

!==================================================================================================
!
! Rotate_xione_2_xitwo_and_back()
!
! Purpose: This constructs the interpolants from xione to xitwo and back.
!
! Comments: This is done by constructing the Vandermonde matrix of the Lagrange basis functions, 
!           constructed from xione, evaluated on xitwo, Ixione2xitwo. The interpolant back 
!           is given by Pxione^-1*Ixione2xitwo^T*Pxitwo, where Pxione and Pxitwo are the 
!           diagonal norm matrices on xione and xitwo, respectively. 
!
! Additional documentation
!
! Inputs: nxione (int): number of nodes in the one-dimensional nodal distribution xione
!         nxitwo (int): number of nodes in the one-dimensional nodal distribution xitwo
!         xione (real(wp)) size (nxioneX1): first nodal distribution
!         xitwo (real(wp)) size (nxitwoX1): second nodal distribution
!         Bxione (real(wp) size (nxioneX1): quadrature weights on the first nodal distribution
!         Bxitwo (real(wp)) size (nxitwoX1): quadrature weights on the second nodal distribution
!
! Outputs:
!         xione2xitwo (real(wp)) size(nxioneXnxitwo): interpolation matrix from xione to xitwo
!         xitwo2xione (real(wp)) size(nxitwoXnxione): interpolation matrix from xitwo to xione
!
!===================================================================================================
   subroutine Rotate_xione_2_xitwo_and_back(nxione,nxitwo,xione,xitwo,Bxione,Bxitwo,xione2xitwo,xitwo2xione)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nxione,nxitwo
     real(wp), intent(in)                        :: xione(nxione),xitwo(nxitwo),Bxione(nxione),Bxitwo(nxitwo)
     real(wp), intent(inout)                     :: xione2xitwo(nxitwo,nxione), xitwo2xione(nxione,nxitwo)

     !-- local variables
     integer                                     :: i,j
     !-- construct the Vandermonde matrix from the Lagrange basis funcitons on xione evaluated at
     !-- the nodal locations of xitwo

     call Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo(nxione,nxitwo,xione,xitwo,xione2xitwo)
    
     !-- construct the interpolant from xitwo to xione
     do i = 1, nxione
       do j = 1,nxitwo
         xitwo2xione(i,j) = Bxione(i)*xione2xitwo(j,i)*Bxitwo(j)
       end do 
     end do 
   end subroutine Rotate_xione_2_xitwo_and_back

!==================================================================================================
!
! Rotate_LGL_2_Gau_and_back_I()
!
! Purpose:
!
! Additional documentation
!
! Inputs:
!
! Outputs:
!
!===================================================================================================
   subroutine Rotate_LGL_2_Gau_and_back_I()
   end subroutine Rotate_LGL_2_Gau_and_back_I

!==================================================================================================
!
! Lagrange_interpolant_basis_1D()
!
! Purpose: Evaluates the jth Lagrange basis function, constructed from the nodes xivec at the 
!          point xi.
!
! Additional documentation
!
! Inputs: 
!       n (int): number of nodes in the one-dimensional nodal distribution xivec
!       xivec (real(wp) vector of size [nX1]): one-dimensional nodal distribution
!       xi_eval (real(wp)): location at which the jth Lagrange basis function is evaluated
!       j (int): which Lagrange basis (numbered 0 through n-1)
!
! Outputs:
! 
!        lj_eval (real(wp)): evaluation of the jth Lagrange basis at xi_eval
!
!===================================================================================================
   subroutine Lagrange_interpolant_basis_1D(n,xivec,xi_eval,j,lj_eval)

     !-- arguments
     integer, intent(in)                         :: n,j
     real(wp), intent(in)                        :: xivec(1:n)
     real(wp), intent(in)                        :: xi_eval
     real(wp), intent(inout)                     :: lj_eval
    
     !-- local variables
     integer                                     :: k
     real(wp)                                    :: wj,l,xij

     !-- construct wj
     wj = 1.0_wp
     xij = xivec(j+1)
     do k = 0,n-1
       if (k.NE.j) then
         wj = wj*(xij-xivec(k+1))
       end if
     end do
     wj = 1.0_wp/wj

     !-- construct l
     l = 0.0_wp
     do k = 0,n-1
       l = l*(xi_eval)
     end do
     lj_eval = 1.0_wp
   end subroutine Lagrange_interpolant_basis_1D

!==================================================================================================
!
! Derivative_Lagrange_interpolant_basis_1D()
!
! Purpose:
!
! Additional documentation
!
! Inputs:
!
! Outputs:
!
!===================================================================================================
   !subroutine Derivative_Lagrange_interpolant_basis_1D()
   !end subroutine Derivative_Lagrange_interpolant_basis_1D
end module non_conforming
