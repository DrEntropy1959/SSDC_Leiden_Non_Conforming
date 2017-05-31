module Jacobian_NavierStokes
  ! This module contains the necessary routines to compute the Jacobi matrix 
  ! of the residual computed with an single diagonal implicit RK scheme for 
  ! the Navier-Stokes equations discretized with a DC method. 
  use precision_vars
  implicit none

  private

  public Element_Jacobian, SAT_Jacobian

contains

  subroutine Element_Jacobian(dt,Akk)

    use CSRlocalvariables
    use variables
    use controlvariables
    use collocationvariables
    use nsereferencevariables
    use referencevariables
    use navierstokes
    use Jacobi_Matrix_Variables

    implicit none

    ! primitive variables
    real(wp),                       intent(in) :: dt,Akk

    integer                                    :: iell , ielh, jdir, idir, jnode
    integer                                    :: inode, ielem, j_node, k_node, j_node_help
    integer                                    :: i, j, k
    integer                                    :: m1, m2, n, n_help
    integer                                    :: nnz, nnzt, cnt, N_tot_proc
    real(wp), dimension(3)                     :: n_div, n_grad
    real(wp), dimension(nequations)            :: vin, vk
    real(wp), dimension(nequations,nequations) :: a_loc
    real(wp), dimension(nequations,nequations) :: t1,t2
    real(wp), dimension(nequations,nequations) :: eye
    real(wp) :: Jx
    logical                                    :: testing = .false.
    real(wp) :: tmp 


    ! Low and high element ID owned by a process
    iell = ihelems(1) ;  ielh = ihelems(2)

    ! Set to zero all the elements in the CSR Jacobi matrix
    dfdu_a = 0.0_wp

    ! Nodes counter
    cnt = 1
    do ielem = iell, ielh
      do inode = 1,nodesperelem

        !  Testing
        if(testing)then
          eye(:,:) = zero
          do i = 1,nequations
            eye(i,i) = 1.0_wp
          enddo
          call primitivevariables(ug(:,inode,ielem),vin(:),nequations)
          t1 = 0.0_wp
          t1 = t1 + abs(eye(:,:) - matmul(dUdV(vin,nequations),dVdU(vin,nequations)))
          t1 = t1 + abs(eye(:,:) - matmul(dWdV(vin,nequations),dVdW(vin,nequations)))
          t1 = t1 + abs(eye(:,:) - matmul(dWdU(vin,nequations),dUdW(vin,nequations)))
          t2 = matmul(dVdW(vin,nequations),dWdU(vin,nequations))
          t1 = t1 + abs(eye(:,:) - matmul(dUdV(vin,nequations),t2))
          if(maxval(t1) >= 1.0e-12)write(*,*)'||dudv.dvdu-eye||',maxval(t1)
        endif
          
        !  Testing

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Time terms
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        a_loc = 0.0_wp     
        do i = 1,nequations
          a_loc(i,i) = one != I: Identity matrix used for diagonal terms arising from the time discretization
        enddo

        do n = ia_0(cnt),ia_0(cnt+1)-1
          do m2 = 1,nequations
            do m1 = 1,nequations
              nnzt  = (iaS(cnt)-1)*nequations*nequations         &
                + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations  &
                + (ka_0(n)-iaS(cnt))*nequations + m1
              dfdu_a(nnzt) = dfdu_a(nnzt) + a_loc(m2,m1) 
            enddo
          enddo
        enddo

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Inviscid Terms
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        ! xi-direction
        ! -------------

        !go to 1111


        ! Assign the value to dfdu_a
        do n = ia_x1(cnt),ia_x1(cnt+1)-1
        
          j_node = ja_x1(n) - (ielem-iell)*nodesperelem

          ! Transform conserved variables to primitive variables
          call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)

          ! Contravariant vector
          n_div = r_x(1,:,j_node,ielem)*Jx_r(j_node,ielem)  

          ! Jacobi matrix of the node
          a_loc = dt*Akk*Eul_Jac_Normal(vin,n_div,nequations)/Jx_r(inode,ielem)

          do m1 = 1,nequations
            do m2 = 1,nequations
              nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                + (ka_x1(n)-iaS(cnt))*nequations + m1
              dfdu_a(nnzt) = dfdu_a(nnzt)  + a_loc(m2,m1) * a_x1(n)

              !write(71,*) dfdu_a(nnzt)
            enddo
          enddo

        enddo


        
        !1111 continue

        !go to 1112

        ! eta-direction
        ! -------------

        do n = ia_x2(cnt),ia_x2(cnt+1)-1

          j_node = ja_x2(n) - (ielem-iell)*nodesperelem

          ! Transform conserved variables to primitive variables
          call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)

          ! Contravariant vector
          n_div = r_x(2,:,j_node,ielem)*Jx_r(j_node,ielem) 
          
          ! Jacobi matrix of the node
          a_loc = dt*Akk*Eul_Jac_Normal(vin,n_div,nequations)/Jx_r(inode,ielem)

          ! Assign the value to dfdu_a
          do m2 = 1,nequations
            do m1 = 1,nequations
              nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                + (ka_x2(n)-iaS(cnt))*nequations + m1
              dfdu_a(nnzt) = dfdu_a(nnzt) + a_loc(m2,m1) * a_x2(n)

              !write(71,*) dfdu_a(nnzt)
            enddo
          enddo

        enddo


        !1112 continue

        ! zeta-direction
        ! --------------

        if(ndim == 3) then

          do n = ia_x3(cnt),ia_x3(cnt+1)-1

            j_node = ja_x3(n) - (ielem-iell)*nodesperelem

            ! Transform conserved variables to primitive variables
            call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)
  
            ! Contravariant vector
            n_div = r_x(3,:,j_node,ielem)*Jx_r(j_node,ielem) 
  
            ! Jacobi matrix of the node
            a_loc = dt*Akk*Eul_Jac_Normal(vin,n_div,nequations)/Jx_r(inode,ielem)
  
            ! Assign the value to dfdu_a
            do m2 = 1,nequations
              do m1 = 1,nequations
                nnzt  = (iaS(cnt)-1)*nequations*nequations       &
                  + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations&
                  + (ka_x3(n)-iaS(cnt))*nequations + m1
                dfdu_a(nnzt) = dfdu_a(nnzt) + a_loc(m2,m1) * a_x3(n)

                !write(71,*) dfdu_a(nnzt)
              enddo
            enddo

          enddo

        endif ! End if 3D


        !1112 continue

        
        ! Update nodes counter
        cnt = cnt + 1

      enddo ! End loop over the nodes

    enddo ! End loop over the elements


    write(71,*) 'Old'
    do i = 1 , size(dfdu_a)
			if (dfdu_a(i) .gt. 0.0_wp) then
      	write(71,*) dfdu_a(i)
			endif
    enddo


    dfdu_a = 0.0_wp

    stop


    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !  Viscous Terms
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if(viscous) then

      go to 10001 

      ! Nodes counter
      cnt = 1

       !go to 10001
      ! Loop over elements
      do ielem = iell, ielh
        grad_u = 0.0_wp
        grad_v = 0.0_wp
        grad_w = 0.0_wp

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  First contribution: dFVisc/dPrim * dPrim/dCons
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        ! Loop over every node in element
        do inode = 1,nodesperelem
          ! compute primitive variables
          call primitivevariables(ug(:,inode,ielem),vin(:),nequations) ! (navierstokes)

          ! Compute gradient of the velocity components
          !
          ! Initialize gradients
          grad_u_tmp = 0.0_wp
          grad_v_tmp = 0.0_wp
          grad_w_tmp = 0.0_wp

          ! Loop over number of dependent nodes in gradient
          do i = iagrad(inode), iagrad(inode+1)-1
            ! Loop over dimensions
            do jdir = 1,ndim
              ! Column/node from gradient operator in CSR format in
              ! the jdir-direction corresponding to the coefficient dagrad(jdir,i)
              jnode = jagrad(jdir,i)
              ! Update gradient using coefficient and primitive variables at appropriate node
              grad_u_tmp(jdir) = grad_u_tmp(jdir) + dagrad(jdir,i)*vg(2,jnode,ielem)
              grad_v_tmp(jdir) = grad_v_tmp(jdir) + dagrad(jdir,i)*vg(3,jnode,ielem)

              if(ndim == 3) then
                grad_w_tmp(jdir) = grad_w_tmp(jdir) + dagrad(jdir,i)*vg(4,jnode,ielem)
              endif
            end do
          end do

          ! Transform gradients to physical space using dxi_jdir/dx_idir (inverse of the Jacobi matrix)
          do jdir = 1,ndim
            do idir = 1,ndim
              grad_u(idir) = grad_u(idir) + grad_u_tmp(jdir)*r_x(jdir,idir,inode,ielem)
              grad_v(idir) = grad_v(idir) + grad_v_tmp(jdir)*r_x(jdir,idir,inode,ielem)

              if(ndim == 3) then
                grad_w(idir) = grad_w(idir) + grad_w_tmp(jdir)*r_x(jdir,idir,inode,ielem)
              endif
            end do
          end do

          
          ! xi-direction
          ! -------------

          !  go to 10004
          ! Assign the value to dfdu_a
          do n = ia_x1(cnt),ia_x1(cnt+1)-1
        
            j_node = ja_x1(n) - (ielem-iell)*nodesperelem

            ! Transform conserved variables to primitive variables
            call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)

            ! Contravariant vector
            n_div = r_x(1,:,j_node,ielem)*Jx_r(j_node,ielem)  

            ! Jacobi matrix of the node
            a_loc = dt*Akk*Vis_Jac_Contrib_1(vin,n_div,nequations)/Jx_r(inode,ielem)

            do m1 = 1,nequations
              do m2 = 1,nequations
                nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                  + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                  + (ka_x1(n)-iaS(cnt))*nequations + m1
                dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1) * a_x1(n)
              enddo
            enddo

          enddo

          
          ! eta-direction
          ! -------------
        
          do n = ia_x2(cnt),ia_x2(cnt+1)-1

            j_node = ja_x2(n) - (ielem-iell)*nodesperelem

            ! Transform conserved variables to primitive variables
            call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)

            ! Contravariant vector
            n_div = r_x(2,:,j_node,ielem)*Jx_r(j_node,ielem) 
          
            ! Jacobi matrix of the node
            a_loc = dt*Akk*Vis_Jac_Contrib_1(vin,n_div,nequations)/Jx_r(inode,ielem)

            ! Assign the value to dfdu_a
            do m2 = 1,nequations
              do m1 = 1,nequations
                nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                  + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                  + (ka_x2(n)-iaS(cnt))*nequations + m1
                dfdu_a(nnzt) = dfdu_a(nnzt) - a_loc(m2,m1) * a_x2(n)
              enddo
            enddo

          enddo


          ! zeta-direction
          ! --------------

          if(ndim == 3) then

            do n = ia_x3(cnt),ia_x3(cnt+1)-1

              j_node = ja_x3(n) - (ielem-iell)*nodesperelem

              ! Transform conserved variables to primitive variables
              call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)
  
              ! Contravariant vector
              n_div = r_x(3,:,j_node,ielem)*Jx_r(j_node,ielem) 
  
              ! Jacobi matrix of the node
              a_loc = dt*Akk*Vis_Jac_Contrib_1(vin,n_div,nequations)/Jx_r(inode,ielem)
  
              ! Assign the value to dfdu_a
              do m2 = 1,nequations
                do m1 = 1,nequations
                  nnzt  = (iaS(cnt)-1)*nequations*nequations       &
                    + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations&
                    + (ka_x3(n)-iaS(cnt))*nequations + m1
                  dfdu_a(nnzt) = dfdu_a(nnzt) - a_loc(m2,m1) * a_x3(n)
                enddo
              enddo

            enddo

          endif ! End if 3D

          ! Update nodes counter
          cnt = cnt + 1 

        end do ! End loop over the nodes

      enddo ! End loop over the elements
    
      10001 continue


      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !  Second contribution: dFVisc/d(gradient(Prim)) * d(gradient(Prim))/dPrim * dPrim/dCons
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Nodes counter
      cnt = 1

      ! Loop over elements
      do ielem = iell, ielh

        ! Loop over every node in element
        do inode = 1,nodesperelem


          !write(*,*) 'inode', inode

         !go to 10002

          ! Laplacian derivatives 
          ! ---------------------

          ! xi-direction
          ! -------------

          ! Assign the value to dfdu_a
          do n = ia_x1(cnt),ia_x1(cnt+1)-1
        
            j_node = ja_x1(n) - (ielem-iell)*nodesperelem 

            !write(*,*) 'j_node', j_node
          
            a_loc = 0.0_wp

            ! We assume that a_x1 and a_x1_transpose have the same dimension for
            ! one node
            n_help = ia_x1(cnt)


            do k = ia_x1_transpose(cnt),ia_x1_transpose(cnt+1)-1
              
              k_node = ja_x1_transpose(k) - (ielem-iell)*nodesperelem

              ! Pseudo-loop over the j_node
              ! ===========================
              j_node_help = ja_x1(n_help) - (ielem-iell)*nodesperelem 
              call primitivevariables(ug(:,j_node_help,ielem),vin(:),nequations) 
              ! First contravariant vectors multiplied by the determinat of the
              ! Jacobi matrix evaluated at the j_node. The determinant of the 
              ! Jacobi matrix has to be added according to the Flux_Divergence in
              ! navierstokes.f90 (tip: look at how the divergence is computed,
              ! i.e. the divergence at the inode is just the summation of the
              ! viscous flux at the jnode multipled by the coefficient which
              ! appears in dagrad).
              n_div  = r_x(1,:,j_node,ielem)*abs(Jx_r(j_node,ielem))


              ! Second contravariant vector 
              n_grad = r_x(1,:,k_node,ielem)

              ! Accumulation of terms
              call primitivevariables(ug(:,k_node,ielem),vk(:),nequations)
              
              !a_loc = a_loc + matmul(Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations),dWdU(vk,nequations)) * a_x1_transpose(k)
              !a_loc = a_loc + Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations) * a_x1_transpose(k)
              a_loc = a_loc + a_x1(n_help)*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)*a_x1_transpose(k)
              !a_loc = Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)

              !a_loc = dWdU(vk,nequations)
              
              !write(*,*) 'dwdu', a_loc

              !write(*,*) 'a_x1', a_x1(n_help)
              !write(*,*) 'a_x1_transpose', a_x1_transpose(n_help)
              !write(*,*) 'a_loc'
              !write(*,*) a_loc(1,1:2), a_loc(1,5)
              !write(*,*) a_loc(2,1:2), a_loc(2,5)
              !write(*,*) a_loc(5,1:2), a_loc(5,5)

              ! Check the values of the differentiation matrix against the 1D mma
              ! code.
              !write(*,*) 'a_x1_transpose(k)', a_x1_transpose(k)

              ! Check the values of DWDU against the 1D mma code
              !write(*,*)'DWDU',dWdU(vk,nequations)

              !write(*,*) 'n_div',n_div
              !write(*,*) 'n_grad',n_grad

              n_help = n_help + 1

            enddo




            ! Element (inode,j_node) of the Jacobi matrix
            !a_loc = dt * akk * a_x1(n) * a_loc / Jx_r(inode,ielem)

            a_loc =  a_loc / abs(Jx_r(inode,ielem))

            write(*,*) 'a_loc'
            write(*,*) a_loc(1,1:2), a_loc(1,5)
            write(*,*) a_loc(2,1:2), a_loc(2,5)
            write(*,*) a_loc(5,1:2), a_loc(5,5)

            !write(*,*) 'Jx_r(inode,ielem)', Jx_r(inode,ielem)
            

            ! Check C11 values against the 1D mma code. THIS ARE CORRECT
            !write(*,*) 'C11_22', 4.0_wp/3.0_wp/gm1M2*Re0inv*vin(5)
            !write(*,*) 'C11_23', 4.0_wp/3.0_wp*mu0*Re0inv*vin(5)*vin(2)
            !write(*,*) 'C11_32', 4.0_wp/3.0_wp*mu0*Re0inv*vin(5)*vin(2)
            !write(*,*) 'C11_35', Re0inv*(k0*vin(5)**2/Pr0+4.0_wp/3.0_wp*gm1M2*mu0*vin(5)*vin(2)**2)

            ! Check the values of the differentiation matrix against the 1D mma
            ! code.
            !write(*,*) 'a_x1(n)', a_x1(n)

            !if(abs(xg(2,inode,ielem)) <= 1.0e-10) then
            !write(*,*) 'inode, j_node', inode, j_node
            !do m1=1,5  
            !  write(*,*) 'a_loc', a_loc(m1,:)
            !enddo
            !write(*,*) 'a_locu', a_loc(2,1:2),a_loc(2,5)
            !write(*,*) 'a_loct', a_loc(5,1:2),a_loc(5,5)
            !endif


            do m1 = 1,nequations
              do m2 = 1,nequations
                nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                  + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                  + (ka_x11(n)-iaS(cnt))*nequations + m1
                dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
              enddo
            enddo

          enddo

          if (inode .eq. 3) then
            stop
          endif



          goto 10002 
          ! eta-direction
          ! --------------

          ! Assign the value to dfdu_a
          do n = ia_x2(cnt),ia_x2(cnt+1)-1
        
            j_node = ja_x2(n) - (ielem-iell)*nodesperelem

            a_loc =0.0_wp

            do k = ia_x2_transpose(j_node),ia_x2_transpose(j_node+1)-1
              
              k_node = ja_x2_transpose(k) - (ielem-iell)*nodesperelem 


              write(*,*) 'inode, j_node, k,k_node',inode, j_node,k_node
              
              ! Transform conserved variables to primitive variables
              call primitivevariables(ug(:,j_node,ielem),vin(:),nequations)

              ! Contravariant vectors
              n_div  = r_x(2,:,j_node,ielem)*Jx_r(j_node,ielem)
              n_grad = r_x(2,:,k_node,ielem)

              ! Accumulation of terms
              a_loc = a_loc + matmul(Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations),dWdU(vin,nequations)) * a_x2_transpose(k)
            enddo

            ! Element (inode,j_node) of the Jacobi matrix
            a_loc = dt * akk * a_x2(n) * a_loc / Jx_r(inode,ielem)

            do m1 = 1,nequations
              do m2 = 1,nequations
                nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                  + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                  + (ka_x22(n)-iaS(cnt))*nequations + m1
                dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
              enddo
            enddo

          enddo

          
          ! zeta-direction
          ! --------------

          if(ndim == 3) then

            ! Assign the value to dfdu_a
            do n = ia_x3(cnt),ia_x3(cnt+1)-1
        
              j_node = ja_x3(n) - (ielem-iell)*nodesperelem

              a_loc =0.0_wp

              do k = ia_x3_transpose(cnt),ia_x3_transpose(cnt+1)-1
              
                k_node = ja_x3_transpose(k) - (ielem-iell)*nodesperelem 

                ! Transform conserved variables to primitive variables
                call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                ! Contravariant vectors
                n_div  = r_x(3,:,j_node,ielem)*Jx_r(j_node,ielem)
                n_grad = r_x(3,:,k_node,ielem)

                ! Jacobi matrix of the node
                a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x3(n)*a_x3_transpose(k)
              enddo

              do m1 = 1,nequations
                do m2 = 1,nequations
                  nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                    + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                    + (ka_x33(n)-iaS(cnt))*nequations + m1
                  dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                enddo
              enddo

            enddo

          endif ! End if 3D

          10002 continue


          ! Second order cross derivatives 
          ! ------------------------------

          if(crossterms) then

            !goto 10001
            !  Dx12

            ! Assign the value to dfdu_a
            do n = ia_x1(cnt),ia_x1(cnt+1)-1
        
              j_node = ja_x1(n) - (ielem-iell)*nodesperelem

              a_loc =0.0_wp
  
              do k = ia_x2_transpose(cnt),ia_x2_transpose(cnt+1)-1

                k_node = ja_x2_transpose(k) - (ielem-iell)*nodesperelem 

                ! Transform conserved variables to primitive variables
                call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                ! Contravariant vectors
                n_div  = r_x(1,:,j_node,ielem)*Jx_r(j_node,ielem)
                n_grad = r_x(2,:,k_node,ielem)

                ! Jacobi matrix of the node
                a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x1(n)*a_x2_transpose(k)
              enddo

              do m1 = 1,nequations
                do m2 = 1,nequations
                  nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                    + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                    + (ka_x12(n)-iaS(cnt))*nequations + m1
                  dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                enddo
              enddo

            enddo

            !  Dx21

            ! Assign the value to dfdu_a
            do n = ia_x2(cnt),ia_x2(cnt+1)-1
        
              j_node = ja_x2(n) - (ielem-iell)*nodesperelem

              do k = ia_x1_transpose(cnt),ia_x1_transpose(cnt+1)-1
              
                k_node = ja_x1_transpose(k) - (ielem-iell)*nodesperelem 

                ! Transform conserved variables to primitive variables
                call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                ! Contravariant vectors
                n_div  = r_x(2,:,j_node,ielem)*Jx_r(j_node,ielem)
                n_grad = r_x(1,:,k_node,ielem)

                ! Jacobi matrix of the node
                a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x2(n)*a_x1_transpose(k)
              enddo

              do m1 = 1,nequations
                do m2 = 1,nequations
                  nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                    + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                    + (ka_x12(n)-iaS(cnt))*nequations + m1
                  dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                enddo
              enddo

            enddo

            !10001 continue



            if(ndim == 3) then

              ! Dx13

              ! Assign the value to dfdu_a
              do n = ia_x1(cnt),ia_x1(cnt+1)-1
        
                j_node = ja_x1(n) - (ielem-iell)*nodesperelem

                do k = ia_x3_transpose(cnt),ia_x3_transpose(cnt+1)-1
              
                  k_node = ja_x3_transpose(k) - (ielem-iell)*nodesperelem 

                  ! Transform conserved variables to primitive variables
                  call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                  ! Contravariant vectors
                  n_div  = r_x(1,:,j_node,ielem)*Jx_r(j_node,ielem)
                  n_grad = r_x(3,:,k_node,ielem)

                  ! Jacobi matrix of the node
                  a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x1(n)*a_x3_transpose(k)
                enddo


                do m1 = 1,nequations
                  do m2 = 1,nequations
                    nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                      + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                      + (ka_x13(n)-iaS(cnt))*nequations + m1
                    dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                  enddo
                enddo

              enddo

             ! Dx23

             ! Assign the value to dfdu_a
              do n = ia_x2(cnt),ia_x2(cnt+1)-1
        
                j_node = ja_x2(n) - (ielem-iell)*nodesperelem

                do k = ia_x3_transpose(cnt),ia_x3_transpose(cnt+1)-1
              
                  k_node = ja_x3_transpose(k) - (ielem-iell)*nodesperelem 

                  ! Transform conserved variables to primitive variables
                  call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                  ! Contravariant vectors
                  n_div  = r_x(2,:,j_node,ielem)*Jx_r(j_node,ielem)
                  n_grad = r_x(3,:,k_node,ielem)

                  ! Jacobi matrix of the node
                  a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x2(n)*a_x3_transpose(k)

                enddo

                do m1 = 1,nequations
                  do m2 = 1,nequations
                    nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                      + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                      + (ka_x23(n)-iaS(cnt))*nequations + m1
                    dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                  enddo
                enddo

              enddo

             !  Dx31

             ! Assign the value to dfdu_a
              do n = ia_x3(cnt),ia_x3(cnt+1)-1
        
                j_node = ja_x3(n) - (ielem-iell)*nodesperelem

                do k = ia_x1_transpose(cnt),ia_x1_transpose(cnt+1)-1
              
                  k_node = ja_x1_transpose(k) - (ielem-iell)*nodesperelem 

                  ! Transform conserved variables to primitive variables
                  call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                  ! Contravariant vectors
                  n_div  = r_x(3,:,j_node,ielem)*Jx_r(j_node,ielem)
                  n_grad = r_x(1,:,k_node,ielem)

                  ! Jacobi matrix of the node
                  a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x3(n)*a_x1_transpose(k)
                enddo

                do m1 = 1,nequations
                  do m2 = 1,nequations
                    nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                      + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                      + (ka_x13(n)-iaS(cnt))*nequations + m1
                    dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                  enddo
                enddo

              enddo


              ! Dx32

              ! Assign the value to dfdu_a
              do n = ia_x3(cnt),ia_x3(cnt+1)-1
        
                j_node = ja_x3(n) - (ielem-iell)*nodesperelem

                do k = ia_x2_transpose(cnt),ia_x2_transpose(cnt+1)-1
  
                  k_node = ja_x2_transpose(k) - (ielem-iell)*nodesperelem 

                  ! Transform conserved variables to primitive variables
                  call primitivevariables(ug(:,k_node,ielem),vin(:),nequations)

                  ! Contravariant vectors
                  n_div  = r_x(3,:,j_node,ielem)*Jx_r(j_node,ielem)
                  n_grad = r_x(2,:,k_node,ielem)

                  ! Jacobi matrix of the node
                  a_loc = dt*Akk*Vis_Jac_Contrib_2_Wvar(vin,n_div,n_grad,nequations)/Jx_r(inode,ielem)*a_x2(n)*a_x3_transpose(k)
                enddo

                do m1 = 1,nequations
                  do m2 = 1,nequations
                    nnzt  = (iaS(cnt)-1)*nequations  * nequations      &
                      + (m2-1)*(iaS(cnt+1)-iaS(cnt)) * nequations      &
                      + (ka_x23(n)-iaS(cnt))*nequations + m1
                    dfdu_a(nnzt) = dfdu_a(nnzt)  - a_loc(m2,m1)
                  enddo
                enddo

              enddo

            endif ! End if 3D

          endif ! Second order cross derivatives 

          ! Update nodes counter
          cnt = cnt + 1 

        enddo ! End loop over the nodes

      enddo ! End loop over the elements

    endif ! End viscous contribution


  end subroutine Element_Jacobian



    subroutine     SAT_Jacobian
    end subroutine SAT_Jacobian

    function Eul_Jac_Normal(vin,nin,nq)
      ! 3D Euler flux Jacobi matrix. This routine also works in 2D
      ! because the number of equations is still 5. However, in 2D, the momentum in the z direction is not 
      ! calculated.

      use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, Pr0, Mach0
      implicit none
      ! Problem sizes
      integer, intent(in) :: nq

      ! Contravariant normal vector multiplied by the Jacobi matrix determinant
      real(wp), intent(in) :: nin(3)

      ! Primitive variables
      real(wp), intent(in) :: vin(nq)

      ! Coefficients of the Flux Jacobi matrix (in computational space coordinates)
      real(wp) :: Eul_Jac_Normal(nq,nq)

      ! Local work variables
      real(wp) :: Ke, Un, M2inv, H


      ! Normal vector multiplied by the determinant Jx of the transformation's Jacobi matrix
      ! Note: Flux_comput_space =  det|J|*J^(-1)*Flux_phys_space

      ! Kinetic energy
      Ke = half*dot_product(vin(2:4),vin(2:4))

      ! Normal velocity
      Un = dot_product(vin(2:4),nin)

      ! 1/M^2
      M2inv = 1.0_wp / Mach0 / Mach0

      ! Enthalpy
      H = (vin(5) + gm1M2*Ke)

      ! Continuity equation
      Eul_Jac_Normal(1,1) = 0.0_wp
      Eul_Jac_Normal(1,2) = nin(1)
      Eul_Jac_Normal(1,3) = nin(2)
      Eul_Jac_Normal(1,4) = nin(3)
      Eul_Jac_Normal(1,5) = 0.0_wp

      ! Momentum equation (mixed up to show patterns) 
      Eul_Jac_Normal(2,1) = -vin(2)*Un  + gm1*Ke*nin(1)
      Eul_Jac_Normal(3,1) = -vin(3)*Un  + gm1*Ke*nin(2)
      Eul_Jac_Normal(4,1) = -vin(4)*Un  + gm1*Ke*nin(3)

      Eul_Jac_Normal(2,2) = +vin(2)*nin(1) - gm1*vin(2)*nin(1) + Un
      Eul_Jac_Normal(2,3) = +vin(2)*nin(2) - gm1*vin(3)*nin(1)
      Eul_Jac_Normal(2,4) = +vin(2)*nin(3) - gm1*vin(4)*nin(1)

      Eul_Jac_Normal(3,2) = +vin(3)*nin(1) - gm1*vin(2)*nin(2)
      Eul_Jac_Normal(3,3) = +vin(3)*nin(2) - gm1*vin(3)*nin(2) + Un
      Eul_Jac_Normal(3,4) = +vin(3)*nin(3) - gm1*vin(4)*nin(2)

      Eul_Jac_Normal(4,2) = +vin(4)*nin(1) - gm1*vin(2)*nin(3)
      Eul_Jac_Normal(4,3) = +vin(4)*nin(2) - gm1*vin(3)*nin(3)
      Eul_Jac_Normal(4,4) = +vin(4)*nin(3) - gm1*vin(4)*nin(3) + Un

      Eul_Jac_Normal(2,5) =  nin(1)*M2inv
      Eul_Jac_Normal(3,5) =  nin(2)*M2inv
      Eul_Jac_Normal(4,5) =  nin(3)*M2inv

      ! Energy equation
      Eul_Jac_Normal(5,1) = + Un*(-vin(5) + (gm1-1.0_wp)*gm1M2*Ke)
      Eul_Jac_Normal(5,2) = + nin(1)*H - gm1*gm1M2*vin(2)*Un
      Eul_Jac_Normal(5,3) = + nin(2)*H - gm1*gm1M2*vin(3)*Un
      Eul_Jac_Normal(5,4) = + nin(3)*H - gm1*gm1M2*vin(4)*Un
      Eul_Jac_Normal(5,5) = gamma0 * Un

    end function Eul_Jac_Normal


    function Vis_Jac_Contrib_1(vin,nin,nq)
      ! First contribution of the viscous flux to the 3D Jacobi matrix. This routine also works in 2D
      ! because the number of equations is still 5. However, in 2D, the momentum in the z direction is not 
      ! calculated.
      use nsereferencevariables, only: gm1M2, Re0inv
      use navierstokes, only: dVdU
      use Jacobi_Matrix_Variables
      implicit none
      ! Problem sizes
      integer, intent(in) :: nq

      ! Contravariant normal vector
      real(wp), intent(in) :: nin(3)

      ! Primitive variables
      real(wp), intent(in) :: vin(nq)

      ! Coefficients of the Flux Jacobi matrix (function of the computational space coordinates)
      real(wp) :: Vis_Jac_Contrib_1(nq,nq)

      ! Matrices of the coefficients
      real(wp), dimension(nq,nq) :: mat


      ! Set all elements to zero
      mat = 0.0_wp

      ! Only the Jacobi matrix of the energy equation is nonzero	
      mat(5,1) = gm1M2*Re0inv/(3.0_wp*vin(1))*(                                                 &
        - 2.0_wp*grad_u(1)*(2.0_wp*nin(1)*vin(2) -        nin(2)*vin(3) -        nin(3)*vin(4)) &
        - 3.0_wp*grad_u(2)*(       nin(1)*vin(3) +        nin(2)*vin(2)                       ) &
        - 3.0_wp*grad_u(3)*(       nin(1)*vin(4)                        +        nin(3)*vin(2)) &

        - 3.0_wp*grad_v(1)*(       nin(1)*vin(3) +        nin(2)*vin(2)                       ) &
        - 2.0_wp*grad_v(2)*(      -nin(1)*vin(2) + 2.0_wp*nin(2)*vin(3) -        nin(3)*vin(4)) &
        - 3.0_wp*grad_v(3)*(                              nin(2)*vin(4) +        nin(3)*vin(3)) &

        - 3.0_wp*grad_w(1)*(       nin(1)*vin(4)                        +        nin(3)*vin(2)) &
        - 3.0_wp*grad_w(2)*(                     +        nin(2)*vin(4) +        nin(3)*vin(3)) &
        - 2.0_wp*grad_w(3)*(      -nin(1)*vin(2) -        nin(2)*vin(3) + 2.0_wp*nin(3)*vin(4)) )

      mat(5,2) = gm1M2*Re0inv/(3.0_wp*vin(1))*( &
        + 4.0_wp*grad_u(1)*nin(1) &
        + 3.0_wp*grad_u(2)*nin(2) &
        + 3.0_wp*grad_u(3)*nin(3) &
        + 3.0_wp*grad_v(1)*nin(2) &
        - 2.0_wp*grad_v(2)*nin(1) &
        + 3.0_wp*grad_w(1)*nin(3) &
        - 2.0_wp*grad_w(3)*nin(1))

      mat(5,3) = gm1M2*Re0inv/(3.0_wp*vin(1))*( &
        - 2.0_wp*grad_u(1)*nin(2) &
        + 3.0_wp*grad_u(2)*nin(1) &
        + 3.0_wp*grad_v(1)*nin(1) &
        + 4.0_wp*grad_v(2)*nin(2) &
        + 3.0_wp*grad_v(3)*nin(3) &
        + 3.0_wp*grad_w(2)*nin(3) &
        - 2.0_wp*grad_w(3)*nin(2))

      mat(5,4) = gm1M2*Re0inv/(3.0_wp*vin(1))*( &
        - 2.0_wp*grad_u(1)*nin(3) &
        + 3.0_wp*grad_u(3)*nin(1) &
        - 2.0_wp*grad_v(2)*nin(3) &
        + 3.0_wp*grad_v(3)*nin(2) &
        + 3.0_wp*grad_w(1)*nin(1) &
        + 3.0_wp*grad_w(2)*nin(2) &
        + 4.0_wp*grad_w(3)*nin(3))


      ! Jacobi matrix wrt to the conserved variables 
      Vis_Jac_Contrib_1 = matmul(mat,dVdU(vin,nq))


    end function Vis_Jac_Contrib_1


!    function Vis_Jac_Contrib_1_Wvar(vin,nin,nq)
!      ! First contribution of the viscous flux to the 3D Jacobi matrix. This routine also works in 2D
!      ! because the number of equations is still 5. However, in 2D, the momentum in the z direction is not 
!      ! calculated.
!      use nsereferencevariables, only: gm1M2, Re0inv
!      use navierstokes, only: dVdU
!      use Jacobi_Matrix_Variables
!      implicit none
!      ! Problem sizes
!      integer, intent(in) :: nq
!
!      ! Contravariant normal vector
!      real(wp), intent(in) :: nin(3)
!
!      ! Primitive variables
!      real(wp), intent(in) :: vin(nq)
!
!      ! Coefficients of the Flux Jacobi matrix (function of the computational space coordinates)
!      real(wp) :: Vis_Jac_Contrib_1_Wvar(nq,nq)
!
!      ! Matrices of the coefficients
!      real(wp), dimension(nq,nq) :: mat
!
!
!      ! Set all elements to zero
!      mat = 0.0_wp
!
!      mat(2,1) = (mu0*(3.0_wp*nin(2)*(((phi(2,2) + phi(3,1) + gm1M2*(phi(5,2)vin(2) + phi(5,1)vin(3)))*(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0 & 
!                 - 2.0_wp*vin(5)))/gm1M2 - 2.0_wp*phi(5,2)*vin(2)*vin(5) - 2.0_wp*phi(5,1)*vin(3)*vin(5)) + 3.0_wp*nin(3)*(((phi(2,3) + phi(4,1) + &
!              gm1M2*(phi(5,3)vin(2) + phi(5,1)vin(4)))*
!            (gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0 - 2.0_wp*vin(5)))/
!          gm1M2 - 2.0_wp*phi(5,3)*vin(2)*vin(5) - 
!         2.0_wp*phi(5,1)*vin(4)*vin(5)) + 
!      nin(1)*(-8*phi(5,1)*vin(2)*vin(5) + 4*phi(5,2)*vin(3)*vin(5) + 
!         4*phi(5,3)*vin(4)*vin(5) + 
!         (2.0_wp*(-2.0_wp*phi(2,1) + phi(3,2) + phi(4,3) - 
!              gm1M2*(2.0_wp*phi(5,1)vin(2) - phi(5,2)vin(3) - phi(5,3)vin(4)))*
!            (-(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0) + 2.0_wp*vin(5)))/
!          gm1M2)))/(6.0_wp*vin(1)) 
!
!
!      mat(2,2) = (mu0*(-(M**2*gamma0**2*vin(2)*
!         (3*phi(5,3)*nin(3)*vin(2) + 3*phi(5,2)*nin(2)*vin(2) + 
!           4.0_wp*phi(5,1)*nin(1)*vin(2) + 3*phi(5,1)*nin(2)*vin(3) - 
!           2.0_wp*phi(5,2)*nin(1)*vin(3) + 3*phi(5,1)*nin(3)*vin(4) - 
!           2.0_wp*phi(5,3)*nin(1)*vin(4))) + 
!      gamma0*vin(2)*(-3.0_wp*phi(2,3)*nin(3) - 3.0_wp*phi(4,1)*nin(3) - 
!         3.0_wp*(phi(2,2) + phi(3,1))*nin(2) + 
!         2.0_wp*(-2.0_wp*phi(2,1) + phi(3,2) + phi(4,3))*nin(1) + 
!         M**2*(3.0_wp*phi(5,3)*nin(3)*vin(2) + 
!            3.0_wp*phi(5,2)*nin(2)*vin(2) + 4.0_wp*phi(5,1)*nin(1)*vin(2) + 
!            3.0_wp*phi(5,1)*nin(2)*vin(3) - 2.0_wp*phi(5,2)*nin(1)*vin(3) + 
!            3.0_wp*phi(5,1)*nin(3)*vin(4) - 2.0_wp*phi(5,3)*nin(1)*vin(4))) + 
!      (3.0_wp*phi(5,3)*nin(3) + 3.0_wp*phi(5,2)*nin(2) + 
!         4.0_wp*phi(5,1)*nin(1))*vin(5)))/(3.0_wp*vin(1))
!
!
!      mat(2,3) = -(mu0*((3.0_wp*phi(2,3)*nin(3) + 3.0_wp*phi(4,1)*nin(3) + 
!          3.0_wp*(phi(2,2) + phi(3,1))*nin(2) + 4.0_wp*phi(2,1)*nin(1) - 
!          2.0_wp*phi(3,2)*nin(1) - 2.0_wp*phi(4,3)*nin(1) + 
!          gm1M2*phi(5,3)*(3.0_wp*nin(3)*vin(2) - 2.0_wp*nin(1)*vin(4)))*gamma0*vin(3) + 
!       phi(5,1)*(gm1M2*vin(3)*(4.0_wp*nin(1)*vin(2) + 3.0_wp*nin(3)*vin(4))*gamma0 - 
!          3.0_wp*nin(2)*(-(gm1M2*vin(3)**2*gamma0) + vin(5))) + 
!       phi(5,2)*(3.0_wp*gm1M2*nin(2)*vin(2)*vin(3)*gamma0 + 
!          2.0_wp*nin(1)*(-(gm1M2*vin(3)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!
!
!      mat(2,4) = -(mu0*((3.0_wp*(phi(2,3) + phi(4,2))*nin(3) + 3.0_wp*phi(2,2)*nin(2) + 
!          3.0_wp*phi(3,1)*nin(2) + 4*phi(2,1)*nin(1) - 
!          2.0_wp*phi(3,2)*nin(1) - 2.0_wp*phi(4,3)*nin(1) + 
!          gm1M2*phi(5,2)*(3.0_wp*nin(2)*vin(2) - 2.0_wp*nin(1)*vin(3)))*gamma0*vin(4) + 
!       phi(5,1)*(gm1M2*(4*nin(1)*vin(2) + 3.0_wp*nin(2)*vin(3))*vin(4)*gamma0 - 
!          3.0_wp*nin(3)*(-(gm1M2*vin(4)**2*gamma0) + vin(5))) + 
!       phi(5,3)*(3.0_wp*gm1M2*nin(3)*vin(2)*vin(4)*gamma0 + 
!          2.0_wp*nin(1)*(-(gm1M2*vin(4)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!      
!      
!      mat(2,5) = (mu0*(3.0_wp*phi(2,3)*nin(3) + 3.0_wp*phi(4,1)*nin(3) + 
!          3.0_wp*(phi(2,2) + phi(3,1))*nin(2) + 
!          2.0_wp*(2.0_wp*phi(2,1) - phi(3,2) - phi(4,3))*nin(1) + 
!          gm1M2*(3.0_wp*phi(5,3)*nin(3)*vin(2) + 3.0_wp*phi(5,2)*nin(2)*vin(2) + 
!         4.0_wp*phi(5,1)*nin(1)*vin(2) + 3.0_wp*phi(5,1)*nin(2)*vin(3) - 
!         2.0_wp*phi(5,2)*nin(1)*vin(3) + 3.0_wp*phi(5,1)*nin(3)*vin(4) - 
!         2.0_wp*phi(5,3)*nin(1)*vin(4)))*gamma0)/(3.0_wp.*gm1M2*vin(1))
!
!
!      mat(3,1) = (mu0*(3.0_wp*nin(1)*(((phi(2,2) + phi(3,1) + 
!              gm1M2*(phi(5,2)*vin(2) + phi(5,1)*vin(3)))*
!            (gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0 - 2.0_wp*vin(5)))/
!          gm1M2 - 2.0_wp*phi(5,2)*vin(2)*vin(5) - 
!         2.0_wp*phi(5,1)*vin(3)*vin(5)) + 
!      3.0_wp*nin(3)*(((phi(3,3) + phi(4,2) + 
!              gm1M2*(phi(5,3)*vin(3) + phi(5,2)*vin(4)))*
!            (gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0 - 2.0_wp*vin(5)))/
!          gm1M2 - 2.0_wp*phi(5,3)*vin(3)*vin(5) - 
!         2.0_wp*phi(5,2)*vin(4)*vin(5)) + 
!      nin(2)*(4.0_wp*phi(5,1)*vin(2)*vin(5) - 8*phi(5,2)*vin(3)*vin(5) + 
!         4.0_wp*phi(5,3)*vin(4)*vin(5) + 
!         (2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3) + 
!              gm1M2*(phi(5,1)*vin(2) - 2.0_wp*phi(5,2)*vin(3) + phi(5,3)*vin(4)))*
!            (-(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0) + 2.0_wp*vin(5)))/
!          gm1M2)))/(6.0_wp*vin(1))
!
!
!       mat(3,2) = -(mu0*((3.0_wp*phi(3,3)*nin(3) + 3.0_wp*phi(4,2)*nin(3) - 
!          2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*nin(2) + 
!          3.0_wp*(phi(2,2) + phi(3,1))*nin(1) + 
!          gm1M2*phi(5,3)*(3.0_wp*nin(3)*vin(3) - 2.0_wp*nin(2)*vin(4)))*gamma0*vin(2) + 
!       phi(5,1)*(3.0_wp*gm1M2*nin(1)*vin(2)*vin(3)*gamma0 + 
!          2.0_wp*nin(2)*(-(gm1M2*vin(2)**2*gamma0) + vin(5))) + 
!       phi(5,2)*(gm1M2*vin(2)*(4.0_wp*nin(2)*vin(3) + 3.0_wp*nin(3)*vin(4))*gamma0 - 
!          3.0_wp*nin(1)*(-(gm1M2*vin(2)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!
!      mat(3,3) = (mu0*(-(M**2*gamma0**2*vin(3)*
!         (-2.0_wp*phi(5,1)*nin(2)*vin(2) + 3.0_wp*phi(5,2)*nin(1)*vin(2) + 
!           3.0_wp*phi(5,3)*nin(3)*vin(3) + 4.0_wp*phi(5,2)*nin(2)*vin(3) + 
!           3.0_wp*phi(5,1)*nin(1)*vin(3) + 3.0_wp*phi(5,2)*nin(3)*vin(4) - 
!           2.0_wp*phi(5,3)*nin(2)*vin(4))) + 
!      gamma0*vin(3)*(-3.0_wp*phi(3,3)*nin(3) - 3.0_wp*phi(4,2)*nin(3) + 
!         2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*nin(2) - 
!         3.0_wp*(phi(2,2) + phi(3,1))*nin(1) + 
!         M**2*(-2.0_wp*phi(5,1)*nin(2) + 3.0_wp*phi(5,2)*nin(1))*vin(2) + 
!         M**2*(3.0_wp*phi(5,3)*nin(3)*vin(3) + 
!            4.0_wp*phi(5,2)*nin(2)*vin(3) + 3.0_wp*phi(5,1)*nin(1)*vin(3) + 
!            3.0_wp*phi(5,2)*nin(3)*vin(4) - 2.0_wp*phi(5,3)*nin(2)*vin(4))) + 
!      (3.0_wp*phi(5,3)*nin(3) + 4.0_wp*phi(5,2)*nin(2) + 
!         3.0_wp*phi(5,1)*nin(1))*vin(5)))/(3.0_wp*vin(1))
!
!      mat(3,4) = -(mu0*((3.0_wp*(phi(3,3) + phi(4,2))*nin(3) - 2.0_wp*phi(2,1)*nin(2) + 
!          4.0_wp*phi(3,2)*nin(2) - 2.0_wp*phi(4,3)*nin(2) + 
!          3.0_wp*(phi(2,2) + phi(3,1))*nin(1) - 
!          gm1M2*phi(5,1)*(2.0_wp*nin(2)*vin(2) - 3.0_wp*nin(1)*v))*gamma0*vin(4) + 
!       phi(5,2)*(gm1M2*(3.0_wp*nin(1)*vin(2) + 4.0_wp*nin(2)*vin(3))*vin(4)*gamma0 - 
!          3.0_wp*nin(3)*(-(gm1M2*vin(4)**2*gamma0) + vin(5))) + 
!       phi(5,3)*(3.0_wp*gm1M2*nin(3)*vin(3)*vin(4)*gamma0 + 
!          2.0_wp*nin(2)*(-(gm1M2*vin(4)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!
!      mat(3,5) = (mu0*(3.0_wp*phi(3,3)*nin(3) + 3.0_wp*phi(4,2)*nin(3) - 
!      2.0_wp*phi(2,1)*nin(2) + 4.0_wp*phi(3,2)*nin(2) - 
!      2.0_wp*phi(4,3)*nin(2) + 3.0_wp*(phi(2,2) + phi(3,1))*nin(1) - 
!      gm1M2*phi(5,1)*(2.0_wp*nin(2)*vin(2) - 3.0_wp*nin(1)*vin(3)) + 
!      gm1M2*phi(5,2)*(3.0_wp*nin(1)*vin(2) + 4.0_wp*nin(2)*vin(3) + 3.0_wp*nin(3)*vin(4)) + 
!      gm1M2*phi(5,3)*(3.0_wp*nin(3)*vin(3) - 2.0_wp*nin(2)*vin(4)))*gamma0)/
!  (3.0_wp*M**2*(-1 + gamma0)*vin(1))
!
!      mat(4,1) = (mu0*(3.0_wp*nin(1)*(((phi(2,3) + phi(4,1) + 
!              gm1M2*(phi(5,3)*vin(2) + phi(5,1)*vin(4)))*
!            (gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0 - 2.0_wp*vin(5)))/
!          gm1M2 - 2.0_wp*phi(5,3)*vin(2)*vin(5) - 
!         2.0_wp*phi(5,1)*vin(4)*vin(5)) + 
!      3.0_wp*nin(2)*(((phi(3,3) + phi(4,2) + 
!              gm1M2*(phi(5,3)*vin(3) + phi(5,2)*vin(4)))*
!            (gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0 - 2.0_wp*vin(5)))/
!          gm1M2 - 2.0_wp*phi(5,3)*vin(3)*vin(5) - 
!         2.0_wp*phi(5,2)*vin(4)*vin(5)) + 
!      nin(3)*(4.0_wp*phi(5,1)*vin(2)*vin(5) + 4.0_wp*phi(5,2)*vin(3)*vin(5) - 
!         8*phi(5,3)*vin(4)*vin(5) + 
!         (2.0_wp*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3) + 
!              gm1M2*(phi(5,1)*vin(2) + phi(5,2)*vin(3) - 2.0_wp*phi(5,3)*vin(4)))*
!            (-(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0) + 2.0_wp*vin(5)))/
!          gm1M2)))/(6.0_wp*vin(1))
!
!      mat(4,2) = -(mu0*(-((2.0_wp*phi(2,1)*nin(3) + 2.0_wp*phi(3,2)*nin(3) - 
!            4.0_wp*phi(4,3)*nin(3) - 
!            3.0_wp*(phi(3,3) + phi(4,2))*nin(2) - 
!            3.0_wp*(phi(2,3) + phi(4,1))*nin(1) + 
!            gm1M2*phi(5,2)*(2.0_wp*nin(3)*vin(3) - 3.0_wp*nin(2)*vin(4)))*gamma0*vin(2)) \
!+ phi(5,1)*(3.0_wp*gm1M2*nin(1)*vin(2)*vin(4)*gamma0 + 
!          2.0_wp*nin(3)*(-(gm1M2*vin(2)**2*gamma0) + vin(5))) + 
!       phi(5,3)*(gm1M2*vin(2)*(3.0_wp*nin(2)*vin(3) + 4.0_wp*nin(3)*vin(4))*gamma0 - 
!          3.0_wp*nin(1)*(-(gm1M2*vin(2)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!
!      mat(4,3) = -(mu0*(-((2.0_wp*phi(2,1)*nin(3) + 2.0_wp*phi(3,2)*nin(3) - 
!            4.0_wp*phi(4,3)*nin(3) - 
!            3.0_wp*(phi(3,3) + phi(4,2))*nin(2) - 
!            3.0_wp*(phi(2,3) + phi(4,1))*nin(1) + 
!            gm1M2*phi(5,1)*(2.0_wp*nin(3)*vin(2) - 3.0_wp*nin(1)*vin(4)))*gamma0*vin(3)) \
!+ phi(5,2)*(3.0_wp*gm1M2*nin(2)*vin(3)*vin(4)*gamma0 + 
!          2.0_wp*nin(3)*(-(gm1M2*vin(3)**2*gamma0) + vin(5))) + 
!       phi(5,3)*(gm1M2*vin(3)*(3.0_wp*nin(1)*vin(2) + 4.0_wp*nin(3)*vin(4))*gamma0 - 
!          3.0_wp*nin(2)*(-(gm1M2*vin(3)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!
!      mat(4,4) = (mu0*(M**2*gamma0**2*vin(4)*(2.0_wp*nin(3)*
!          (phi(5,1)*vin(2) + phi(5,2)*vin(3)) - 
!         3.0_wp*(phi(5,2)*nin(2) + phi(5,1)*nin(1))*vin(4)) + 
!      gamma0*vin(4)*(-3.0_wp*(phi(3,3) + phi(4,2))*nin(2) - 
!         3.0_wp*(phi(2,3) + phi(4,1))*nin(1) + 
!         2.0_wp*nin(3)*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3) - 
!            M**2*(phi(5,1)*vin(2) + phi(5,2)*vin(3))) + 
!         3.0_wp*M**2*(phi(5,2)*nin(2) + phi(5,1)*nin(1))*vin(4)) + 
!      3.0_wp*(phi(5,2)*nin(2) + phi(5,1)*nin(1))*vin(5) + 
!      phi(5,3)*(-3.0_wp*gm1M2*(nin(1)*vin(2) + nin(2)*vin(3))*vin(4)*gamma0 + 
!         4.0_wp*nin(3)*(-(gm1M2*vin(4)**2*gamma0) + vin(5)))))/(3.0_wp*vin(1))
!
!      mat(4,5) = -(mu0*(2.0_wp*phi(2,1)*nin(3) + 2.0_wp*phi(3,2)*nin(3) - 
!       4.0_wp*phi(4,3)*nin(3) - 3.0_wp*(phi(3,3) + phi(4,2))*nin(2) - 
!       3.0_wp*(phi(2,3) + phi(4,1))*nin(1) - 
!       gm1M2*phi(5,3)*(3.0_wp*nin(1)*vin(2) + 3.0_wp*nin(2)*vin(3) + 4.0_wp*nin(3)*vin(4)) + 
!       gm1M2*phi(5,2)*(2.0_wp*nin(3)*vin(3) - 3.0_wp*nin(2)*vin(4)) + 
!       gm1M2*phi(5,1)*(2.0_wp*nin(3)*vin(2) - 3.0_wp*nin(1)*vin(4)))*gamma0)/
!  (3.0_wp*M**2*(-1 + gamma0)*vin(1))
!
!      mat(5,1) = (nin(1)*(2.0_wp*mu0*(-4.0_wp*phi(2,1) + 2.0_wp*phi(3,2) + 2.0_wp*phi(4,3) - 
!          gm1M2*(8*phi(5,1)*vin(2) + phi(5,2)*vin(3) + phi(5,3)*vin(4)))*vin(2)*
!        vin(5) - 2.0_wp*mu0*(3.0_wp*phi(2,2) + 3.0_wp*phi(3,1) + 
!          gm1M2*(phi(5,2)*vin(2) + 6.0_wp*phi(5,1)*vin(3)))*vin(3)*vin(5) - 
!       2.0_wp*mu0*(3.0_wp*phi(2,3) + 3.0_wp*phi(4,1) + 
!          gm1M2*(phi(5,3)*vin(2) + 6.0_wp*phi(5,1)*vin(4)))*vin(4)*vin(5) - 
!       ((-(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0) + 2.0_wp*vin(5))*
!          (2.0_wp*(2.0_wp*phi(2,1) - phi(3,2) - phi(4,3))*mu0*Pr0*vin(2) + 
!            mu0*Pr0*(gm1M2*(phi(5,2)*vin(2)*vin(3) + phi(5,3)*vin(2)*vin(4) + 
!                  phi(5,1)*(4.0_wp*vin(2)**2 + 3.0_wp*(v**2 + w**2))) + 
!               3.0_wp*(phi(2,2) + phi(3,1))*vin(3) + 
!               3.0_wp*(phi(2,3) + phi(4,1))*vin(4)) + 
!            6.0_wp*phi(5,1)*k0*vin(5)))/Pr0) + 
!    nin(2)*(-2.0_wp*mu0*(3.0_wp*phi(2,2) + 3.0_wp*phi(3,1) + 
!          gm1M2*(6.0_wp*phi(5,2)*vin(2) + phi(5,1)*vin(3)))*vin(2)*vin(5) + 
!       2.0_wp*mu0*(2.0_wp*phi(2,1) - 4.0_wp*phi(3,2) + 2.0_wp*phi(4,3) - 
!          gm1M2*(phi(5,1)*vin(2) + 8*phi(5,2)*vin(3) + phi(5,3)*vin(4)))*vin(3)*
!        vin(5) - 2.0_wp*mu0*(3.0_wp*phi(3,3) + 3.0_wp*phi(4,2) + 
!          gm1M2*(phi(5,3)*vin(3) + 6.0_wp*phi(5,2)*vin(4)))*vin(4)*vin(5) - 
!       ((-(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0) + 2.0_wp*vin(5))*
!          (mu0*Pr0*(gm1M2*(phi(5,1)*vin(2)*vin(3) + phi(5,3)*vin(3)*vin(4) + 
!                  phi(5,2)*(3.0_wp*vin(2)**2 + 4.0_wp*vin(3)**2 + 3.0_wp*vin(4)**2)) + 
!               3.0_wp*phi(2,2)*vin(2) + 3.0_wp*phi(3,1)*vin(2) - 
!               2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*vin(3) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(4)) + 
!            6.0_wp*phi(5,2)*k0*vin(5)))/Pr0) + 
!    nin(3)*(-2.0_wp*mu0*(3.0_wp*phi(2,3) + 3.0_wp*phi(4,1) + 
!          gm1M2*(6.0_wp*phi(5,3)*vin(2) + phi(5,1)*vin(4)))*vin(2)*vin(5) - 
!       2.0_wp*mu0*(3.0_wp*phi(3,3) + 3.0_wp*phi(4,2) + 
!          gm1M2*(6.0_wp*phi(5,3)*vin(3) + phi(5,2)*vin(4)))*vin(3)*vin(5) + 
!       2.0_wp*mu0*(2.0_wp*phi(2,1) + 2.0_wp*phi(3,2) - 4.0_wp*phi(4,3) - 
!          gm1M2*(phi(5,1)*vin(2) + phi(5,2)*vin(3) + 8*phi(5,3)*vin(4)))*vin(4)*
!        vin(5) - ((-(gm1M2*dot_product(vin(2:4),vin(2:4))*gamma0) + 2.0_wp*vin(5))*
!          (mu0*Pr0*(gm1M2*(3.0_wp*phi(5,3)*(u**2 + v**2) + 
!                  (phi(5,1)*vin(2) + phi(5,2)*vin(3))*vin(4) + 4.0_wp*phi(5,3)*vin(4)**2) \
!+ 3.0_wp*phi(2,3)*vin(2) + 3.0_wp*phi(4,1)*vin(2) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(3) - 
!               2.0_wp*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3))*vin(4)) + 
!            6.0_wp*phi(5,3)*k0*vin(5)))/Pr0))/(6.0_wp*vin(1))
!
!      mat(5,2) = (3.0_wp*nin(2)*(2.0_wp*gm1M2*phi(5,2)*vin(T)*vin(2)*mu0 + 
!       (gm1M2*phi(5,1)*vin(T)*vin(3)*mu0)/3. - 
!       (gm1M2*vin(2)*gamma0*(6.0_wp*phi(5,2)*vin(T)*k0 + 
!            Pr0*(3.0_wp*phi(2,2)*vin(2) + 3.0_wp*phi(3,1)*vin(2) - 
!               2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*vin(3) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(4) + 
!               M**2*(phi(5,1)*vin(2)*vin(3) + phi(5,3)*vin(3)*vin(4) + 
!                  phi(5,2)*(3.0_wp*vin(2)**2 + 4.0_wp*vin(3)**2 + 3.0_wp*vin(4)**2))*(-1 + gamma0)\
!)*mu0))/(3.0_wp*Pr0) + (phi(2,2) + phi(3,1))*mu0*vin(5)) + 
!    3.0_wp*nin(3)*(2.0_wp*gm1M2*phi(5,3)*vin(T)*vin(2)*mu0 + 
!       (gm1M2*phi(5,1)*vin(T)*vin(4)*mu0)/3. - 
!       (gm1M2*vin(2)*gamma0*(6.0_wp*phi(5,3)*vin(T)*k0 + 
!            Pr0*(3.0_wp*phi(2,3)*vin(2) + 3.0_wp*phi(4,1)*vin(2) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(3) - 
!               2.0_wp*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3))*vin(4) + 
!               M**2*(3.0_wp*phi(5,3)*(u**2 + v**2) + 
!                  (phi(5,1)*vin(2) + phi(5,2)*vin(3))*vin(4) + 4.0_wp*phi(5,3)*vin(4)**2\
!)*(-1 + gamma0))*mu0))/(3.0_wp*Pr0) + (phi(2,3) + phi(4,1))*mu0*vin(5)\
!) + nin(1)*(8*gm1M2*phi(5,1)*vin(T)*vin(2)*mu0 + 
!       gm1M2*vin(T)*(phi(5,2)*vin(3) + phi(5,3)*vin(4))*mu0 - 
!       (gm1M2*vin(2)*gamma0*(6.0_wp*phi(5,1)*vin(T)*k0 + 
!            2.0_wp*(2.0_wp*phi(2,1) - phi(3,2) - phi(4,3))*Pr0*vin(2)*mu0 + 
!            Pr0*(3.0_wp*(phi(2,2) + phi(3,1))*vin(3) + 
!               3.0_wp*(phi(2,3) + phi(4,1))*vin(4) + 
!               M**2*(phi(5,2)*vin(2)*vin(3) + phi(5,3)*vin(2)*vin(4) + 
!                  phi(5,1)*(4.0_wp*vin(2)**2 + 3.0_wp*(v**2 + w**2)))*(-1 + gamma0))*
!             mu0))/Pr0 + 2.0_wp*(2.0_wp*phi(2,1) - phi(3,2) - phi(4,3))*mu0*
!        vin(5)))/(3.0_wp*vin(1))
!
!      mat(5,3) = (3.0_wp*nin(1)*((gm1M2*phi(5,2)*vin(T)*vin(2)*mu0)/3. + 
!       2.0_wp*gm1M2*phi(5,1)*vin(T)*vin(3)*mu0 - 
!       (gm1M2*vin(3)*gamma0*(6.0_wp*phi(5,1)*vin(T)*k0 + 
!            2.0_wp*(2.0_wp*phi(2,1) - phi(3,2) - phi(4,3))*Pr0*vin(2)*mu0 + 
!            Pr0*(3.0_wp*(phi(2,2) + phi(3,1))*vin(3) + 
!               3.0_wp*(phi(2,3) + phi(4,1))*vin(4) + 
!               M**2*(phi(5,2)*vin(2)*vin(3) + phi(5,3)*vin(2)*vin(4) + 
!                  phi(5,1)*(4.0_wp*vin(2)**2 + 3.0_wp*(v**2 + w**2)))*(-1 + gamma0)\
!)*mu0))/(3.0_wp*Pr0) + (phi(2,2) + phi(3,1))*mu0*vin(5)) + 
!    3.0_wp*nin(3)*(2.0_wp*gm1M2*phi(5,3)*vin(T)*vin(3)*mu0 + 
!       (gm1M2*phi(5,2)*vin(T)*vin(4)*mu0)/3. - 
!       (gm1M2*vin(3)*gamma0*(6.0_wp*phi(5,3)*vin(T)*k0 + 
!            Pr0*(3.0_wp*phi(2,3)*vin(2) + 3.0_wp*phi(4,1)*vin(2) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(3) - 
!               2.0_wp*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3))*vin(4) + 
!               M**2*(3.0_wp*phi(5,3)*(u**2 + v**2) + 
!                  (phi(5,1)*vin(2) + phi(5,2)*vin(3))*vin(4) + 4.0_wp*phi(5,3)*vin(4)**2\
!)*(-1 + gamma0))*mu0))/(3.0_wp*Pr0) + (phi(3,3) + phi(4,2))*mu0*vin(5)\
!) + nin(2)*(8*gm1M2*phi(5,2)*vin(T)*vin(3)*mu0 + 
!       gm1M2*vin(T)*(phi(5,1)*vin(2) + phi(5,3)*vin(4))*mu0 - 
!       (gm1M2*vin(3)*gamma0*(6.0_wp*phi(5,2)*vin(T)*k0 + 
!            Pr0*(3.0_wp*phi(2,2)*vin(2) + 3.0_wp*phi(3,1)*vin(2) - 
!               2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*vin(3) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(4) + 
!               M**2*(phi(5,1)*vin(2)*vin(3) + phi(5,3)*vin(3)*vin(4) + 
!                  phi(5,2)*(3.0_wp*vin(2)**2 + 4.0_wp*vin(3)**2 + 3.0_wp*vin(4)**2))*(-1 + gamma0))*
!             mu0))/Pr0 - 2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*mu0*
!        vin(5)))/(3.0_wp*vin(1))
!
!      mat(5,4) = (3.0_wp*nin(1)*((gm1M2*phi(5,3)*vin(T)*vin(2)*mu0)/3. + 
!       2.0_wp*gm1M2*phi(5,1)*vin(T)*vin(4)*mu0 - 
!       (gm1M2*vin(4)*gamma0*(6.0_wp*phi(5,1)*vin(T)*k0 + 
!            2.0_wp*(2.0_wp*phi(2,1) - phi(3,2) - phi(4,3))*Pr0*vin(2)*mu0 + 
!            Pr0*(3.0_wp*(phi(2,2) + phi(3,1))*vin(3) + 
!               3.0_wp*(phi(2,3) + phi(4,1))*vin(4) + 
!               M**2*(phi(5,2)*vin(2)*vin(3) + phi(5,3)*vin(2)*vin(4) + 
!                  phi(5,1)*(4.0_wp*vin(2)**2 + 3.0_wp*(v**2 + w**2)))*(-1 + gamma0)\
!)*mu0))/(3.0_wp*Pr0) + (phi(2,3) + phi(4,1))*mu0*vin(5)) + 
!    3.0_wp*nin(2)*((gm1M2*phi(5,3)*vin(T)*vin(3)*mu0)/3. + 
!       2.0_wp*gm1M2*phi(5,2)*vin(T)*vin(4)*mu0 - 
!       (gm1M2*vin(4)*gamma0*(6.0_wp*phi(5,2)*vin(T)*k0 + 
!            Pr0*(3.0_wp*phi(2,2)*vin(2) + 3.0_wp*phi(3,1)*vin(2) - 
!               2.0_wp*(phi(2,1) - 2.0_wp*phi(3,2) + phi(4,3))*vin(3) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(4) + 
!               M**2*(phi(5,1)*vin(2)*vin(3) + phi(5,3)*vin(3)*vin(4) + 
!                  phi(5,2)*(3.0_wp*vin(2)**2 + 4.0_wp*vin(3)**2 + 3.0_wp*vin(4)**2))*(-1 + gamma0)\
!)*mu0))/(3.0_wp*Pr0) + (phi(3,3) + phi(4,2))*mu0*vin(5)) + 
!    nin(3)*(gm1M2*vin(T)*(phi(5,1)*vin(2) + phi(5,2)*vin(3))*mu0 + 
!       8*gm1M2*phi(5,3)*vin(T)*vin(4)*mu0 - 
!       (gm1M2*vin(4)*gamma0*(6.0_wp*phi(5,3)*vin(T)*k0 + 
!            Pr0*(3.0_wp*phi(2,3)*vin(2) + 3.0_wp*phi(4,1)*vin(2) + 
!               3.0_wp*(phi(3,3) + phi(4,2))*vin(3) - 
!               2.0_wp*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3))*vin(4) + 
!               M**2*(3.0_wp*phi(5,3)*(u**2 + v**2) + 
!                  (phi(5,1)*vin(2) + phi(5,2)*vin(3))*vin(4) + 4.0_wp*phi(5,3)*vin(4)**2)*
!                (-1 + gamma0))*mu0))/Pr0 - 
!       2.0_wp*(phi(2,1) + phi(3,2) - 2.0_wp*phi(4,3))*mu0*vin(5)))/(3.0_wp*vin(1))
!
!      mat(5,5) = (gamma0*(mu0*Pr0*(gm1M2*phi(5,3)*
!          (3.0_wp*nin(3)*(u**2 + v**2) + (nin(1)*vin(2) + nin(2)*vin(3))*vin(4) + 
!            4.0_wp*nin(3)*vin(4)**2) + 
!         gm1M2*phi(5,2)*(3.0_wp*nin(2)*vin(2)**2 + nin(1)*vin(2)*v + 4.0_wp*nin(2)*vin(3)**2 + 
!            nin(3)*vin(3)*vin(4) + 3.0_wp*nin(2)*vin(4)**2) + 
!         gm1M2*phi(5,1)*(4.0_wp*nin(1)*vin(2)**2 + nin(2)*vin(2)*vin(3) + 3.0_wp*nin(1)*vin(3)**2 + 
!            nin(3)*vin(2)*vin(4) + 3.0_wp*nin(1)*vin(4)**2) + 
!         3.0_wp*phi(2,2)*(nin(2)*vin(2) + nin(1)*vin(3)) + 
!         3.0_wp*phi(3,1)*(nin(2)*vin(2) + nin(1)*vin(3)) - 
!         2.0_wp*phi(4,3)*(nin(1)*vin(2) + nin(2)*vin(3) - 
!            2.0_wp*nin(3)*vin(4)) - 
!         2.0_wp*phi(3,2)*(nin(1)*vin(2) - 2.0_wp*nin(2)*vin(3) + 
!            nin(3)*vin(4)) - 
!         2.0_wp*phi(2,1)*(-2.0_wp*nin(1)*vin(2) + nin(2)*vin(3) + 
!            nin(3)*vin(4)) + 
!         3.0_wp*phi(3,3)*(nin(3)*vin(3) + nin(2)*vin(4)) + 
!         3.0_wp*phi(4,2)*(nin(3)*vin(3) + nin(2)*vin(4)) + 
!         3.0_wp*phi(2,3)*(nin(3)*vin(2) + nin(1)*vin(4)) + 
!         3.0_wp*phi(4,1)*(nin(3)*vin(2) + nin(1)*vin(4))) + 
!      6.0_wp*(phi(5,3)*nin(3) + phi(5,2)*nin(2) + phi(5,1)*nin(1))*
!       k0*vin(5)))/(3.0_wp*Pr0*vin(1))
!      
!      
!      ! Jacobi matrix wrt to the conserved variables 
!      Vis_Jac_Contrib_1_Wvar = matmul(mat,dVdU(vin,nq))
!
!
!    end function Vis_Jac_Contrib_1_Wvar



    function Vis_Jac_Contrib_2(vin,ni,nj,nq)
      ! Second contribution of the viscous flux to the 3D Jacobi matrix. This routine also works in 2D
      ! because the number of equations is still 5. However, in 2D, the momentum in the z direction is not 
      ! calculated.
      use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, Pr0, Mach0, k0
      use navierstokes, only: dVdU
      use Jacobi_Matrix_Variables
      implicit none
      ! Problem sizes
      integer, intent(in) :: nq

      ! Contravariant normal vector
      real(wp), intent(in) :: ni(3), nj(3)

      ! Primitive variables
      real(wp), intent(in) :: vin(nq)

      ! Coefficients of the Flux Jacobi matrix (function of the computational space coordinates)
      real(wp) :: Vis_Jac_Contrib_2(nq,nq)

      ! Matrices of the coefficients
      real(wp), dimension(nq,nq) :: mat

      ! Work constants
      real(wp) :: con1, con2

      ! Set all elements to zero
      mat = 0.0_wp

      con1 = Re0inv/3.0_wp
      con2 = Re0inv*gm1M2

      ! Momentum equation
      mat(2,2) = con1*( 4*ni(1)*nj(1) + 3*ni(2)*nj(2) + 3*ni(3)*nj(3))
      mat(2,3) = con1*( 3*ni(2)*nj(1) - 2*ni(1)*nj(2)                )
      mat(2,4) = con1*( 3*ni(3)*nj(1)                 - 2*ni(1)*nj(3))

      mat(3,2) = con1*(-2*ni(2)*nj(1) + 3*ni(1)*nj(2)                )
      mat(3,3) = con1*( 3*ni(1)*nj(1) + 4*ni(2)*nj(2) + 3*ni(3)*nj(3))
      mat(3,4) = con1*(                 3*ni(3)*nj(2) - 2*ni(2)*nj(3))

      mat(4,2) = con1*(-2*ni(3)*nj(1)                 + 3*ni(1)*nj(3))
      mat(4,3) = con1*(               - 2*ni(3)*nj(2) + 3*ni(2)*nj(3))
      mat(4,4) = con1*( 3*ni(1)*nj(1) + 3*ni(2)*nj(2) + 4*ni(3)*nj(3))

      ! Energy equation
      mat(5,2) = con2*(mat(2,2)*vin(2) + mat(3,2)*vin(3) + mat(4,2)*vin(4)) 
      mat(5,3) = con2*(mat(2,3)*vin(2) + mat(3,3)*vin(3) + mat(4,3)*vin(4))
      mat(5,4) = con2*(mat(2,4)*vin(2) + mat(3,4)*vin(3) + mat(4,4)*vin(4))
      mat(5,5) = Re0inv*gamma0*k0/Pr0*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3)) 


      ! Jacobi matrix wrt to the conserved variables 
      Vis_Jac_Contrib_2 = matmul(mat,dVdU(vin,nq))

    end function Vis_Jac_Contrib_2

    function Vis_Jac_Contrib_2_Wvar(vin,ni,nj,nq)
      ! Second contribution of the viscous flux to the 3D Jacobi matrix. This routine also works in 2D
      ! because the number of equations is still 5. However, in 2D, the momentum in the z direction is not 
      ! calculated.
      use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, Pr0, Mach0, k0, mu0
      use navierstokes, only: dWdU
      use Jacobi_Matrix_Variables
      implicit none
      ! Problem sizes
      integer, intent(in) :: nq

      ! Contravariant normal vector
      real(wp), intent(in) :: ni(3), nj(3)

      ! Primitive variables
      real(wp), intent(in) :: vin(nq)

      ! Coefficients of the Flux Jacobi matrix (function of the computational space coordinates)
      real(wp) :: Vis_Jac_Contrib_2_Wvar(nq,nq)

      ! Matrices of the coefficients
      real(wp), dimension(nq,nq) :: mat

      ! Work constants
      real(wp) :: con1, con2, con3, con4
      real(wp) :: u, v, w, T

      ! Set all elements to zero
      mat = 0.0_wp

      u   = vin(2) ; v   = vin(3) ; w   = vin(4) ; T   = vin(5)

      con1 = Re0inv * mu0 * T / gm1M2 / 3.0_wp
      con2 = gm1M2
      con3 = gm1M2
      con4 = Re0inv * k0 * T * T / Pr0 

      ! Momentum equation
      mat(2,2) = con1*( 4*ni(1)*nj(1) + 3*ni(2)*nj(2) + 3*ni(3)*nj(3))
      mat(2,3) = con1*( 3*ni(2)*nj(1) - 2*ni(1)*nj(2)                )
      mat(2,4) = con1*( 3*ni(3)*nj(1)                 - 2*ni(1)*nj(3))

      mat(3,2) = con1*(-2*ni(2)*nj(1) + 3*ni(1)*nj(2)                )
      mat(3,3) = con1*( 3*ni(1)*nj(1) + 4*ni(2)*nj(2) + 3*ni(3)*nj(3))
      mat(3,4) = con1*(                 3*ni(3)*nj(2) - 2*ni(2)*nj(3))

      mat(4,2) = con1*(-2*ni(3)*nj(1)                 + 3*ni(1)*nj(3))
      mat(4,3) = con1*(               - 2*ni(3)*nj(2) + 3*ni(2)*nj(3))
      mat(4,4) = con1*( 3*ni(1)*nj(1) + 3*ni(2)*nj(2) + 4*ni(3)*nj(3))

      mat(2,5) = con2*(mat(2,2)*u + mat(2,3)*v + mat(2,4)*w) 
      mat(3,5) = con2*(mat(3,2)*u + mat(3,3)*v + mat(3,4)*w)
      mat(4,5) = con2*(mat(4,2)*u + mat(4,3)*v + mat(4,4)*w)

      ! Energy equation
      mat(5,2) = con2*(u*mat(2,2) + v*mat(3,2) + w*mat(4,2)) 
      mat(5,3) = con2*(u*mat(2,3) + v*mat(3,3) + w*mat(4,3))
      mat(5,4) = con2*(u*mat(2,4) + v*mat(3,4) + w*mat(4,4))

      mat(5,5) = con3*(mat(5,2)*u + mat(5,3)*v + mat(5,4)*w) &
               + con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3)) 


      ! Jacobi matrix wrt to the conserved variables 
!     Vis_Jac_Contrib_2_Wvar = matmul(mat,dWdU(vin,nq))
      Vis_Jac_Contrib_2_Wvar = mat(:,:)

    end function Vis_Jac_Contrib_2_Wvar

    function Vis_Jac_Contrib_1_Wvar(vin,dw,ni,nj,nq)
      ! Second contribution of the viscous flux to the 3D Jacobi matrix. This routine also works in 2D
      ! because the number of equations is still 5. However, in 2D, the momentum in the z direction is not 
      ! calculated.
      use nsereferencevariables, only: gamma0, gm1, gm1M2, gm1og, gM2, Re0inv, Pr0, Mach0, k0, mu0
      use navierstokes, only: dWdU, dVdU
      use Jacobi_Matrix_Variables
      implicit none
      ! Problem sizes
      integer, intent(in) :: nq

      ! Contravariant normal vector
      real(wp), intent(in) :: ni(3), nj(3)

      ! Primitive variables
      real(wp), intent(in) :: vin(nq),dw(nq)

      ! Coefficients of the Flux Jacobi matrix (function of the computational space coordinates)
      real(wp) :: Vis_Jac_Contrib_1_Wvar(nq,nq)

      ! Tmp Vector to hold C_{i,k} dw_k
      real(wp), dimension(nq)    :: tmp

      ! Matrices of the coefficients
      real(wp), dimension(nq,nq) :: mat

      ! Work constants
      real(wp) :: con1, con2, con3, con4
      real(wp) :: u, v, w, T

      ! Set all elements to zero
      mat = 0.0_wp

      u   = vin(2) ; v   = vin(3) ; w   = vin(4) ; T   = vin(5)

      con1 = Re0inv * mu0 * T / gm1M2 / 3.0_wp
      con2 = gm1M2
      con3 = gm1M2
      con4 = Re0inv * k0 * T * T / Pr0 

      ! Momentum equation
      mat(2,2) = con1*( 4*ni(1)*nj(1) + 3*ni(2)*nj(2) + 3*ni(3)*nj(3))
      mat(2,3) = con1*( 3*ni(2)*nj(1) - 2*ni(1)*nj(2)                )
      mat(2,4) = con1*( 3*ni(3)*nj(1)                 - 2*ni(1)*nj(3))

      mat(3,2) = con1*(-2*ni(2)*nj(1) + 3*ni(1)*nj(2)                )
      mat(3,3) = con1*( 3*ni(1)*nj(1) + 4*ni(2)*nj(2) + 3*ni(3)*nj(3))
      mat(3,4) = con1*(                 3*ni(3)*nj(2) - 2*ni(2)*nj(3))

      mat(4,2) = con1*(-2*ni(3)*nj(1)                 + 3*ni(1)*nj(3))
      mat(4,3) = con1*(               - 2*ni(3)*nj(2) + 3*ni(2)*nj(3))
      mat(4,4) = con1*( 3*ni(1)*nj(1) + 3*ni(2)*nj(2) + 4*ni(3)*nj(3))

      mat(2,5) = con2*(mat(2,2)*u + mat(2,3)*v + mat(2,4)*w) 
      mat(3,5) = con2*(mat(3,2)*u + mat(3,3)*v + mat(3,4)*w)
      mat(4,5) = con2*(mat(4,2)*u + mat(4,3)*v + mat(4,4)*w)

      ! Energy equation
      mat(5,2) = con2*(u*mat(2,2) + v*mat(3,2) + w*mat(4,2)) 
      mat(5,3) = con2*(u*mat(2,3) + v*mat(3,3) + w*mat(4,3))
      mat(5,4) = con2*(u*mat(2,4) + v*mat(3,4) + w*mat(4,4))

      mat(5,5) = con3*(mat(5,2)*u + mat(5,3)*v + mat(5,4)*w) &
               + con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3)) 


      Vis_Jac_Contrib_1_Wvar = 0.0_wp

      Vis_Jac_Contrib_1_Wvar(2,2) = con2 * mat(2,2) * dW(5) 
      Vis_Jac_Contrib_1_Wvar(3,2) = con2 * mat(3,2) * dW(5) 
      Vis_Jac_Contrib_1_Wvar(4,2) = con2 * mat(4,2) * dW(5) 

      Vis_Jac_Contrib_1_Wvar(2,3) = con2 * mat(2,3) * dW(5) 
      Vis_Jac_Contrib_1_Wvar(3,3) = con2 * mat(3,3) * dW(5) 
      Vis_Jac_Contrib_1_Wvar(4,3) = con2 * mat(4,3) * dW(5) 

      Vis_Jac_Contrib_1_Wvar(2,4) = con2 * mat(2,4) * dW(5) 
      Vis_Jac_Contrib_1_Wvar(3,4) = con2 * mat(3,4) * dW(5) 
      Vis_Jac_Contrib_1_Wvar(4,4) = con2 * mat(4,4) * dW(5) 

      Vis_Jac_Contrib_1_Wvar(5,2) = con2 * ( mat(2,2)*dW(2) + mat(2,3)*dW(3) + mat(2,4)*dW(4) + (mat(2,5)+mat(5,2))*dW(5) ) 
      Vis_Jac_Contrib_1_Wvar(5,3) = con2 * ( mat(3,2)*dW(2) + mat(3,3)*dW(3) + mat(3,4)*dW(4) + (mat(3,5)+mat(5,3))*dW(5) ) 
      Vis_Jac_Contrib_1_Wvar(5,4) = con2 * ( mat(4,2)*dW(2) + mat(4,3)*dW(3) + mat(4,4)*dW(4) + (mat(4,5)+mat(5,4))*dW(5) ) 

      Vis_Jac_Contrib_1_Wvar(2,5) = ( mat(2,2)*dW(2) + mat(2,3)*dW(3) + mat(2,4)*dW(4) + mat(2,5)*dW(5) ) / T
      Vis_Jac_Contrib_1_Wvar(3,5) = ( mat(3,2)*dW(2) + mat(3,3)*dW(3) + mat(3,4)*dW(4) + mat(3,5)*dW(5) ) / T
      Vis_Jac_Contrib_1_Wvar(4,5) = ( mat(4,2)*dW(2) + mat(4,3)*dW(3) + mat(4,4)*dW(4) + mat(4,5)*dW(5) ) / T
      Vis_Jac_Contrib_1_Wvar(5,5) = ( mat(5,2)*dW(2) + mat(5,3)*dW(3) + mat(5,4)*dW(4) + mat(5,5)*dW(5) ) / T + &
        & 2.0_wp*(con4*(ni(1)*nj(1) + ni(2)*nj(2) + ni(3)*nj(3)))*dW(5)/T
     
      Vis_Jac_Contrib_1_Wvar = matmul(Vis_Jac_Contrib_1_Wvar,dVdU(vin,nq))

    end function Vis_Jac_Contrib_1_Wvar


  end module Jacobian_NavierStokes
