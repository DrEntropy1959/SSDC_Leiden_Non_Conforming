program Rotate_xione_2_xitwo_and_back_program
!==================================================================================================
!
! Purpose: this tests the subroutine Rotate_xione_2_xitwo_and_back against 
!          Vandermonde matrices constructed using the Maple script Construct_1D_interpolants
!==================================================================================================
  !-- load modules
  use non_conforming
  use precision_vars
  use initcollocation

  !-- Nothing is implicitely definted
  implicit none

  !-- local variables
  integer                                        :: nxione,nxitwo
  real(wp),allocatable                           :: xiGL(:), wxiGL(:), xiG(:), wxiG(:),&
                                                    VandermondeGL(:,:), VandermondeG(:,:),&
                                                    IGL2G(:,:), IG2GL(:,:)

  !-- loop through all possible combinations of GL to G mortar and check the interpolation accuracy
  do nxione = 2,18
    do nxitwo = nxione,18
      allocate(xiGL(nxione))
      allocate(xiG(nxitwo))
      allocate(wxiGL(nxione))
      allocate(wxiG(nxitwo))
      allocate(VandermondeGL(nxione,nxione))
      allocate(VandermondeG(nxitwo,nxione))
      allocate(IGL2G(nxitwo,nxione))
      allocate(IG2GL(nxione,nxitwo))
      call Gauss_Lobatto_Legendre_points(nxione,xiGL,wxiGL)
      call Gauss_Legendre_points(nxitwo,xiG,wxiG)
      call Vandermonde_1D_monomial(nxione,nxione-1,xiGL,VandermondeGL)
      call Vandermonde_1D_monomial(nxitwo,nxione-1,xiG,VandermondeG)

      call  Rotate_xione_2_xitwo_and_back(nxione,nxitwo,xiGL,xiG,wxiGL,wxiG,IGL2G,IG2GL)
      if (maxval(abs(matmul(IGL2G,VandermondeGL)-VandermondeG))>10.0_wp**(-14)) then
        write(*,*)"Error above tolerance 10^-14, nxione = ",nxione," nxitwo = ",nxitwo 
        write(*,*)"Max error IGL2G = ",maxval(abs(matmul(IGL2G,VandermondeGL)-VandermondeG))
      end if
      if (maxval(abs(matmul(IG2GL,VandermondeG(1:nxitwo,1:nxione-1))-VandermondeGL(1:nxione,1:nxione-1)))>10.0_wp**(-14)) then
        write(*,*)"Error above tolerance 10^-14, nxione = ",nxione," nxitwo = ",nxitwo 
        write(*,*)"Max error IG2GL = ",maxval(abs(matmul(IG2GL,VandermondeG(1:nxitwo,1:nxione-1))-VandermondeGL(1:nxione,1:nxione-1)))
      endif

      !-- deallocate statments
      deallocate(xiGL)
      deallocate(xiG)
      deallocate(wxiGL)
      deallocate(wxiG)
      deallocate(VandermondeGL)
      deallocate(VandermondeG)
      deallocate(IGL2G)
      deallocate(IG2GL)
    end do
  end do
end program Rotate_xione_2_xitwo_and_back_program
