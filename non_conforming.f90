module non_conforming

  use precision_vars

  !-- Nothing is implicitely defined
  implicit none

  !-- private subroutines functions etc
  private

  !-- public subroutines functions etc
  public Lagrange_interpolant_basis_1D

  

contains
!==================================================================================================
!
! Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo()
!
! Purpose: Constructs the Vandermonde matrix from the Lagrange basis functions on the nodal set 
!          xione evaluated at the nodes of xitwo
!
! Additional documentation
!
! Unit tests
!
! Inputs: nxione (int): number of nodes in the nodal set xione
!         nxitwo (int): number of nodes in the nodal set xitwo
!         xione (real(wp)): nodal set that is used to construct the Lagrange basis functions
!         xitwo (real(wp)): nodal set at which the Lagrange basis functions will be evaluated at
!
! Outputs:
!         Vandermonde (real(wp) (nxitwoXnxione): the above mentioned Vandermonde matrix
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
     integer                                     :: j, k
     real(wp)                                    :: wj(nxione),l

     !-- construct the wjs
     do j = 0,nxione-1 
       wj(j+1) = 1.0_wp
       do k = 0,nxione-1
         if (k.NE.j) then
           wj(j+1) = wj(j+1)*(xione(j+1)-xione(k+1))
         end if
       end do
       wj(j+1) = 1.0_wp/wj(j+1)
     end do

     !-- construct the Vandermonde matrix
     do j = 1,nxione
       do k = 1,nxitwo
         Vandermonde(j,k)  = l*wj(k)/(xitwo(j)-xione(k))
       end do
     end do

   end subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo

!==================================================================================================
!
! Rotate_LGL_H_2_Gau_H_2_LGL_H()
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
   subroutine Rotate_LGL_H_2_Gau_H_2_LGL_H()

   end subroutine Rotate_LGL_H_2_Gau_H_2_LGL_H

!==================================================================================================
!
! Rotate_LGL_L_2_Gau_H_2_LGL_L()
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
   subroutine Rotate_LGL_L_2_Gau_H_2_LGL_L()
   end subroutine Rotate_LGL_L_2_Gau_H_2_LGL_L

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
     real(wp)                                    :: wj,l,lj,xij,xi

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
   subroutine Derivative_Lagrange_interpolant_basis_1D()
   end subroutine Derivative_Lagrange_interpolant_basis_1D
end module non_conforming
