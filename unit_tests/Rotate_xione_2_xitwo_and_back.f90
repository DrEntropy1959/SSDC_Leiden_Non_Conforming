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
  integer                                        :: i,j
  real(wp),allocatable                           :: xiGL(:), wxiGL(:), xiG(:), wxiG(:),&
                                                    VandermondeGL(:,:), VandermondeG(:,:),&
                                                    IGL2G(:,:), IG2GL(:,:)

  !-- loop through all possible combinations of GL to G mortar and check the interpolation accuracy
  do i = 2,18
    do j = 2,18
      allocate(xiGL(i),xiG(j),wxiGL(i),wxiG(j),VandermondeGL(i,i-1),VandermondeG(j,i-1))
      allocate(IGL2G(j,i),IG2GL(i,j))
      call Gauss_Lobatto_Legendre_points(i,xiGL,wxiGL)
      call Gauss_Legendre_points(j,xiG,wxiG)
      call Vandermonde_1D_monomial(i,i-1,xiGL,VandermondeGL)
      call Vandermonde_1D_monomial(j,i-1,xiG,VandermondeG)

      write(*,*)"Size of IGL2G = ",size(IGL2G(:,1)),"x",size(IGL2G(1,:))
      write(*,*)"Size of IGL2G = ",size(IG2GL(:,1)),"x",size(IG2GL(1,:))
      write(*,*)"Size of VandermondeGL = ",size(VandermondeGL(:,1)),"x",size(VandermondeGL(1,:))
      write(*,*)"Size of VandermondeG = ",size(VandermondeG(:,1)),"x",size(VandermondeG(1,:))

      call  Rotate_xione_2_xitwo_and_back(i,j,xiGL,xiG,wxiGL,wxiG,IGL2G,IG2GL)
      write(*,*)"Max error IGL2G = ",maxval(abs(matmul(IGL2G,VandermondeGL)-VandermondeG))
      write(*,*)"Max error IG2GL = ",maxval(abs(matmul(IG2GL,VandermondeG)-VandermondeGL))
      deallocate(xiGL,xiG,wxiGL,wxiG,VandermondeGL,VandermondeG)
      deallocate(IGL2G,IG2GL)
    end do
  end do
end program Rotate_xione_2_xitwo_and_back_program
