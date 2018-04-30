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
  public Rotate_GL_2_G_and_back_I
  public h_refine_boundary
contains
!==================================================================================================
!
!
! Purpose: This subroutine h refines the mesh at the boundaries
! 
! Comments: hardcoded for 3D
!
! Outputs:
!         
!
!===================================================================================================
   subroutine h_refine_boundary()

     use referencevariables, only : nelems, nfaces, ndim
     use variables, only : ef2e, e_edge2e
     use variables, only : vx_master, e2v, vx, iae2e, jae2e, iae2v, jae2v, ic2nh, xg

     !-- local variables
     integer :: ielem, iface, bc_element_count, max_partners
     integer, allocatable, dimension(:,:,:) :: ef2e_temp                                            !-- (7 ,nfaceperelem, nelements)
     integer, allocatable, dimension(:,:,:,:,:) :: e_edge2e_temp                                    !-- (3,number_of_edges_per_face,max_partners,nfaceperelem,nelems)

     integer :: i_err


     !-- first we loop over all elements to determine how many boundary elements we have
     bc_element_count = 0
     e_loop_bc : do ielem = 1,nelems

       f_loop_bc: do iface = 1,nfaces

         !-- check if the element has a face with a boundary condition
         f_if_bc: if (ef2e(1,iface,ielem) < 0) then
           bc_element_count = bc_element_count+1
         endif f_if_bc

       end do f_loop_bc

     enddo e_loop_bc

     !-- allocate the necessary temporary arrays
     if(ndim.EQ.1)then
       call PetscFinalize(i_err); stop
     elseif(ndim.EQ.2)then
       nelems = nelems-bc_element_count+4*bc_element_count 
     elseif(ndim.EQ.3)then
       nelems = nelems-bc_element_count+8*bc_element_count                                           
     endif
     allocate(ef2e_temp(7,nfaces,nelems))
     max_partners = size(e_edge2e(1,1,:,1,1))
     allocate(e_edge2e_temp(3,2**ndim,max_partners,nfaces,nelems))


   end subroutine h_refine_boundary
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
!           the barycentric weights, wj, are computed as wj = 1/(Pi_{k=0^{n-1,K!=j}(x_{j}-x_{k})}
!           the Lagrange basis functions are constructed as l_{j} =l*w_{j}/(x-x_{j}) where 
!           l = Pi_{k=0}^{n-1}(x-x_{j})
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
           l = l*(xitwo(k)-xione(i))
         end do
         !-- check to see if the nodal location on the second set of nodes matches the node 
         !-- location associated with the kth Lagrange basis function 
         if (abs(xitwo(k)-xione(j)) <=2.0_wp*epsilon(1.0_wp)) then
           Vandermonde(k,j) = 1.0_wp
         else
           Vandermonde(k,j)  = l*wj(j)/(xitwo(k)-xione(j))
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
         xitwo2xione(i,j) = 1.0_wp/Bxione(i)*(xione2xitwo(j,i)*Bxitwo(j))
       end do 
     end do 
   end subroutine Rotate_xione_2_xitwo_and_back

!==================================================================================================
!
! Rotate_GL_2_G_and_back_I()
!
! Purpose: This gives back the interpolant from GL to G and back
!
! Comments: The interpolants were solved for in the Maple script maple/Construct_1D_interpolants_Identity.mw
!
! Additional documentation
!
! Inputs: nGL (int): number of nodes in the one-dimensional Gauss Lobatto (GL) nodal distribution
!         nG (int): number of nodes in the one-dimensional Gauss (G) nodal distribution
!         xGL (real(wp)) size (nxioneX1): GL nodal distribution
!         xG (real(wp)) size (nxitwoX1): G nodal distribution
!         BxiGL (real(wp) size (nxioneX1): quadrature weights on the GL nodal distribution
!         BxiG (real(wp)) size (nxitwoX1): quadrature weights on the G nodal distribution
!
! Outputs:
!         IGL2G (real(wp)) size(nGXnGL): interpolation matrix from the GL nodes to the G nodes
!         IG2GL (real(wp)) size(nGLXnG): interpolation matrix from the G nodes to GL nodes
!
!===================================================================================================
   subroutine Rotate_GL_2_G_and_back_I(nGL,nG,xiGL,xiG,BxiGL,BxiG,IGL2G,IG2GL)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nGL,nG
     real(wp), intent(in)                        :: xiGL(nGL),xiG(nG),BxiGL(nGL),BxiG(nG)
     real(wp), intent(inout)                     :: IGL2G(nG,nGL), IG2GL(nGL,nG)

     !-- local variables
     

    if (nGL==2.AND.nG==2) then
      IGL2G(1,1) = 1.000000000000000000000000000000_wp
      IGL2G(1,2) = 0.000000000000000000000000000000_wp
      IGL2G(2,1) = 0.000000000000000000000000000000_wp
      IGL2G(2,2) = 1.000000000000000000000000000000_wp
      
      IG2GL(1,1) = 1.000000000000000000000000000000_wp
      IG2GL(1,2) = 0.000000000000000000000000000000_wp
      IG2GL(2,1) = 0.000000000000000000000000000000_wp
      IG2GL(2,2) = 1.000000000000000000000000000000_wp
    else if (nGL==3.AND.nG==3) then
      IGL2G(1,1) = 0.343146490609516399718000303683_wp
      IGL2G(1,2) = 1.088303688022450577599852472590_wp
      IGL2G(1,3) = -0.431450178631966977317852776273_wp
      IGL2G(2,1) = 0.430189805014031610999907795369_wp
      IGL2G(2,2) = 0.139620389971936778000184409261_wp
      IGL2G(2,3) = 0.430189805014031610999907795369_wp
      IGL2G(3,1) = -0.431450178631966977317852776273_wp
      IGL2G(3,2) = 1.088303688022450577599852472590_wp
      IGL2G(3,3) = 0.343146490609516399718000303683_wp
      
      IG2GL(1,1) = 0.571910817682527332863333839473_wp
      IG2GL(1,2) = 1.147172813370750962666420787650_wp
      IG2GL(1,3) = -0.719083631053278295529754627123_wp
      IG2GL(2,1) = 0.453459870009354407333271863580_wp
      IG2GL(2,2) = 0.093080259981291185333456272841_wp
      IG2GL(2,3) = 0.453459870009354407333271863580_wp
      IG2GL(3,1) = -0.719083631053278295529754627123_wp
      IG2GL(3,2) = 1.147172813370750962666420787650_wp
      IG2GL(3,3) = 0.571910817682527332863333839473_wp      
    else if (nGL==4.AND.nG==4) then 
      IGL2G(1,1) = 0.670133597061877309375642466273_wp
      IGL2G(1,2) = 0.382690211445752349043404473625_wp
      IGL2G(1,3) = -0.059634895378013858655525747055_wp
      IGL2G(1,4) = 0.006811086870384200236478807157_wp
      IGL2G(2,1) = -0.124994022375115610291702237717_wp
      IGL2G(2,2) = 1.094392949182378686886803453950_wp
      IGL2G(2,3) = 0.011123163321311394153889248044_wp
      IGL2G(2,4) = 0.019477909871425529251009535715_wp
      IGL2G(3,1) = 0.019477909871425529251009535715_wp
      IGL2G(3,2) = 0.011123163321311394153889248052_wp
      IGL2G(3,3) = 1.094392949182378686886803453960_wp
      IGL2G(3,4) = -0.124994022375115610291702237717_wp
      IGL2G(4,1) = 0.006811086870384200236478807157_wp
      IGL2G(4,2) = -0.059634895378013858655525747055_wp
      IGL2G(4,3) = 0.382690211445752349043404473625_wp
      IGL2G(4,4) = 0.670133597061877309375642466273_wp
      
      IG2GL(1,1) = 1.398655311764185408399936481650_wp
      IG2GL(1,2) = -0.489085476472273963755402275747_wp
      IG2GL(1,3) = 0.076214547296997107993785971074_wp
      IG2GL(1,4) = 0.014215617411091447361679823021_wp
      IG2GL(2,1) = 0.159744773085697986032413787639_wp
      IG2GL(2,2) = 0.856443671190025127617787284033_wp
      IG2GL(2,3) = 0.008704700480085614686859237835_wp
      IG2GL(2,4) = -0.024893144755808728337060309506_wp
      IG2GL(3,1) = -0.024893144755808728337060309506_wp
      IG2GL(3,2) = 0.008704700480085614686859237829_wp
      IG2GL(3,3) = 0.856443671190025127617787284040_wp
      IG2GL(3,4) = 0.159744773085697986032413787639_wp
      IG2GL(4,1) = 0.014215617411091447361679823021_wp
      IG2GL(4,2) = 0.076214547296997107993785971074_wp
      IG2GL(4,3) = -0.489085476472273963755402275747_wp
      IG2GL(4,4) = 1.398655311764185408399936481650_wp
    else if (nGL==5.AND.nG==5) then
      IGL2G(1,1) = 0.470502966472490507460316236133_wp
      IGL2G(1,2) = 0.803126567593681665261854299464_wp
      IGL2G(1,3) = -0.491470918189824019737687496142_wp
      IGL2G(1,4) = 0.369914191164023953735669915576_wp
      IGL2G(1,5) = -0.152072807040372106720152955027_wp
      IGL2G(2,1) = 0.071767876791755725608492888402_wp
      IGL2G(2,2) = 0.529448395544771508275461087673_wp
      IGL2G(2,3) = 0.689001782387354883935218360338_wp
      IGL2G(2,4) = -0.492612611092600584063108759494_wp
      IGL2G(2,5) = 0.202394556368718466243936423086_wp
      IGL2G(3,1) = -0.187499999999999999999999999998_wp
      IGL2G(3,2) = 0.437500000000000000000000000001_wp
      IGL2G(3,3) = 0.500000000000000000000000000006_wp
      IGL2G(3,4) = 0.437500000000000000000000000001_wp
      IGL2G(3,5) = -0.187499999999999999999999999998_wp
      IGL2G(4,1) = 0.202394556368718466243936423086_wp
      IGL2G(4,2) = -0.492612611092600584063108759494_wp
      IGL2G(4,3) = 0.689001782387354883935218360338_wp
      IGL2G(4,4) = 0.529448395544771508275461087673_wp
      IGL2G(4,5) = 0.071767876791755725608492888402_wp
      IGL2G(5,1) = -0.152072807040372106720152955027_wp
      IGL2G(5,2) = 0.369914191164023953735669915576_wp
      IGL2G(5,3) = -0.491470918189824019737687496142_wp
      IGL2G(5,4) = 0.803126567593681665261854299464_wp
      IGL2G(5,5) = 0.470502966472490507460316236133_wp
      
      IG2GL(1,1) = 1.114748022560237464834855958820_wp
      IG2GL(1,2) = 0.343501634534003810651995244455_wp
      IG2GL(1,3) = -1.066666666666666666666666666660_wp
      IG2GL(1,4) = 0.968718374310688038674062097476_wp
      IG2GL(1,5) = -0.360301364738262647494246634073_wp
      IG2GL(2,1) = 0.349498057896440617267939081403_wp
      IG2GL(2,2) = 0.465445435697663304036304550688_wp
      IG2GL(2,3) = 0.457142857142857142857142857140_wp
      IG2GL(2,4) = -0.433062586135970603458639415191_wp
      IG2GL(2,5) = 0.160976235399009539297252925968_wp
      IG2GL(3,1) = -0.163747509950278330491814509669_wp
      IG2GL(3,2) = 0.463747509950278330491814509666_wp
      IG2GL(3,3) = 0.400000000000000000000000000005_wp
      IG2GL(3,4) = 0.463747509950278330491814509666_wp
      IG2GL(3,5) = -0.163747509950278330491814509669_wp
      IG2GL(4,1) = 0.160976235399009539297252925968_wp
      IG2GL(4,2) = -0.433062586135970603458639415191_wp
      IG2GL(4,3) = 0.457142857142857142857142857140_wp
      IG2GL(4,4) = 0.465445435697663304036304550688_wp
      IG2GL(4,5) = 0.349498057896440617267939081403_wp
      IG2GL(5,1) = -0.360301364738262647494246634073_wp
      IG2GL(5,2) = 0.968718374310688038674062097476_wp
      IG2GL(5,3) = -1.066666666666666666666666666660_wp
      IG2GL(5,4) = 0.343501634534003810651995244455_wp
      IG2GL(5,5) = 1.114748022560237464834855958820_wp     
    else if (nGL==6.AND.nG==6) then

    else if (nGL==7.AND.nG==7) then
    else if (nGL==8.AND.nG==8) then
    else if (nGL==9.AND.nG==9) then
    else if (nGL==10.AND.nG==10) then
    else if (nGL==11.AND.nG==11) then
    else if (nGL==12.AND.nG==12) then
    else if (nGL==13.AND.nG==13) then
    else if (nGL==15.AND.nG==15) then
    else if (nGL==16.AND.nG==16) then
    else if (nGL==17.AND.nG==17) then
    end if 
   end subroutine Rotate_GL_2_G_and_back_I

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
