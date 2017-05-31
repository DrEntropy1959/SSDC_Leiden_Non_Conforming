module matrix_mod

  ! module contains:

  use math_mod
  use matvec_module

  implicit none

contains

  !=============================================================================80

  subroutine CSR_Term_Pointers(CSR,                                        &
      bc_type,viscous_flow,cross_terms,             & 
      Dxi,Deta,Dzeta,                               &
      D2xi,D2eta,D2zeta,                            &
      DCxi,DCeta,DCzeta)

    !  Subroutine forms pointers to ia, ja from each class of terms 

    !  Classes include  1)  \Xi      Convective  operator: d  / dx
    !                   2)  \Eta     Convective  operator: d  / dy
    !                   3)  \Zeta    Convective  operator: d  / dz
    !                   4)  \Xi      Dissipation operator: d^2/dx^2
    !                   5)  \Eta     Dissipation operator: d^2/dy^2
    !                   6)  \Zeta    Dissipation operator: d^2/dz^2
    !                   7)  \XiEta   Dissipation operator: d^2/dxdy
    !                   8)  \XiZeta  Dissipation operator: d^2/dxdz
    !                   9)  \EtaZeta Dissipation operator: d^2/dydz

    type(diff_op),                  intent(in)    :: Dxi, Deta, Dzeta
    type(diff_op),                  intent(in)    :: D2xi,D2eta,D2zeta
    type(diff_op),                  intent(in)    :: DCxi,DCeta,DCzeta

    type(CSR_Matrix),               intent(inout) :: CSR

    integer,  dimension(:,:),       intent(in)    :: bc_type
    integer,                        intent(in)    :: viscous_flow,cross_terms

    integer,  dimension(9)                        :: statusX

    continue

    statusX(:) = 0

    !  0)   Time  Convective  operator: d  / dt     :  Local pointers

    call CSR_Get_Pointers_time(CSR%ia0,CSR%ja0,CSR%ka0,CSR%a0, bc_type, CSR%nt0, CSR%nz0)

    !
    !  1)  \Xi    Convective  operator: d  / dXi    :  Local pointers
    !  2)  \Eta   Convective  operator: d  / dEta   :  Local pointers
    !  3)  \Zeta  Convective  operator: d  / dZeta  :  Local pointers
    !
    call CSR_Get_Pointers_Xi  (CSR%ia1,CSR%ja1,CSR%ka1,CSR%a1, Dxi   , bc_type, CSR%nt1, CSR%nz1)
    call CSR_Get_Pointers_Eta (CSR%ia2,CSR%ja2,CSR%ka2,CSR%a2, Deta  , bc_type, CSR%nt2, CSR%nz2)
    if(nz(1) > 1) call CSR_Get_Pointers_Zeta(CSR%ia3,CSR%ja3,CSR%ka3,CSR%a3, Dzeta , bc_type, CSR%nt3, CSR%nz3)

    !
    !  4)  \Xi    Dissipation  operator: d^2  / dXi^2    :  Local pointers
    !  5)  \Eta   Dissipation  operator: d^2  / dEta^2   :  Local pointers
    !  6)  \Zeta  Dissipation  operator: d^2  / dZeta^2  :  Local pointers
    !
    if(viscous_flow == 1) then
      call CSR_Get_Pointers_Xi  (CSR%ia4,CSR%ja4,CSR%ka4,CSR%a4, D2xi   , bc_type, CSR%nt4, CSR%nz4)
      call CSR_Get_Pointers_Eta (CSR%ia5,CSR%ja5,CSR%ka5,CSR%a5, D2eta  , bc_type, CSR%nt5, CSR%nz5)
      if(nz(1) > 1) call CSR_Get_Pointers_Zeta(CSR%ia6,CSR%ja6,CSR%ka6,CSR%a6, D2zeta , bc_type, CSR%nt6, CSR%nz6)

      !
      !    7)  \Xi    D1 Derivative operator: d / dXi     for Cross-terms   :  Local pointers
      !    8)  \Eta   D1 Derivative operator: d / dEta    for Cross-terms  :  Local pointers
      !    9)  \Zeta  D1 Derivative operator: d / dZeta   for Cross-terms :  Local pointers
      !
      if(cross_terms == 1) then
        call CSR_Get_Pointers_Xi  (CSR%ia7,CSR%ja7,CSR%ka7,CSR%a7, DCxi   , bc_type, CSR%nt7, CSR%nz7)
        call CSR_Get_Pointers_Eta (CSR%ia8,CSR%ja8,CSR%ka8,CSR%a8, DCeta  , bc_type, CSR%nt8, CSR%nz8)
        if(nz(1) > 1) call CSR_Get_Pointers_Zeta(CSR%ia9,CSR%ja9,CSR%ka9,CSR%a9, DCzeta , bc_type, CSR%nt9, CSR%nz9)
      endif
    endif

    statusX(1) = CSR%N_tot - CSR%nt1 ; statusX(2) = CSR%N_tot - CSR%nt2 ; if(nz(1) > 1 ) statusX(3) = CSR%N_tot - CSR%nt3 ;
    if(viscous_flow == 1) then
      statusX(4) = CSR%N_tot - CSR%nt4 ; statusX(5) = CSR%N_tot - CSR%nt5 ; if(nz(1) > 1 ) statusX(6) = CSR%N_tot - CSR%nt6 ;
      if(cross_terms == 1) then
        statusX(7) = CSR%N_tot - CSR%nt7 ; statusX(8) = CSR%N_tot - CSR%nt8 ; if(nz(1) > 1 ) statusX(9) = CSR%N_tot - CSR%nt9 ;
      endif
    endif
    if(sum(statusX) /= 0) then ; write(*,*)'somethings up with CSR_TERM_POINTERS; Stopping';stop;endif

      call CSR_Global_Pointers(CSR, bc_type,viscous_flow,cross_terms) 

      return
    end subroutine CSR_Term_Pointers

    !=============================================================================80

    subroutine CSR_Get_Pointers_time(ia,ja,ka,a, bc_type, n_tot, nnz)

      integer,  dimension(:),         pointer       :: ia
      integer,  dimension(:),         pointer       :: ja
      integer,  dimension(:),         pointer       :: ka
      real(dp), dimension(:),         pointer       ::  a

      integer,  dimension(:,:),       intent(in)    :: bc_type

      integer,                        intent(inout) :: n_tot,nnz

      integer,  dimension(4)                        :: statusX
      integer                                       :: i,j,k,L,n
      integer                                       :: cnt
      integer                                       :: iell, ielh

      continue

      iell = ihelems(1) ;  ielh = ihelems(2)

      n_tot = 0
      nnz = 0
      do L = iell, ielh
        n_tot = n_tot + nodesperelem
        nnz =   nnz + nodesperelem
      enddo

      allocate(ia(n_tot+1),stat=statusX(1))
      allocate(ja(  nnz  ),stat=statusX(2))
      allocate(ka(  nnz  ),stat=statusX(3))
      allocate( a(  nnz  ),stat=statusX(4))
      if(sum(statusX) > 0) then
        write(*,*)'trouble allocating CSR in CSR_Get_Pointers'
        write(*,*)'\Xi Convective Operator'
        stop
      endif

      nnz   = 1
      cnt   = 1
      ia(cnt) = nnz

      do n=iell,ielh
        do i=1,nodesperelem
          ia(cnt+1) = ia(cnt) + 1
          ja(cnt)   = cnt
          a(cnt)   = one
          cnt    = cnt + 1 
        enddo
      enddo

      return
    end subroutine CSR_Get_Pointers_time

    !=============================================================================80

    subroutine CSR_Get_Pointers_Xi(ia,ja,ka,a, Dxi, bc_type, n_tot, nnz)

      type(diff_op),                  intent(in)    :: Dxi

      integer,  dimension(:),         pointer       :: ia
      integer,  dimension(:),         pointer       :: ja
      integer,  dimension(:),         pointer       :: ka
      real(dp), dimension(:),         pointer       ::  a

      integer,  dimension(:,:),       intent(in)    :: bc_type

      integer,                        intent(inout) :: n_tot,nnz

      integer,  dimension(4)                        :: statusX
      integer                                       :: c2,c3,c4,c5,c6,c7,c8
      integer                                       :: i,j,k,L,n
      integer                                       :: cnt

      continue

      n_tot = 0
      nnz = 0
      do L = 1, nblock
        n_tot = n_tot + nodesperelem
        nnz =   nnz + (nx(L)-2*c4) * ny(L) * nz(L) * c6    &
          + (      2*c4) * ny(L) * nz(L) * c5
      enddo

      write(*,*)'n_tot,nnz,xi   direction',n_tot,nnz

      allocate(ia(n_tot+1),stat=statusX(1))
      allocate(ja(  nnz  ),stat=statusX(2))
      allocate(ka(  nnz  ),stat=statusX(3))
      allocate( a(  nnz  ),stat=statusX(4))
      if(sum(statusX) > 0) then
        write(*,*)'trouble allocating CSR in CSR_Get_Pointers'
        write(*,*)'\Xi Convective Operator'
        stop
      endif

      nnz   = 1
      cnt   = 1
      ia(cnt) = nnz

      N_tot = 0
      do n=iell,ielh
        do i = 1,e
        end do
        N_tot = N_tot + nx(n)*ny(n)*nz(n)
      end do
      nnz = nnz - 1    !  subtract off the last outer

      return
    end subroutine CSR_Get_Pointers_Xi

    !=============================================================================80

    subroutine CSR_Get_Pointers_Eta(ia,ja,ka,a, Deta, bc_type, n_tot, nnz)

      type(diff_op),                  intent(in)    :: Deta

      integer,  dimension(:),         pointer       :: ia
      integer,  dimension(:),         pointer       :: ja
      integer,  dimension(:),         pointer       :: ka
      real(dp), dimension(:),         pointer       ::  a

      integer,  dimension(:,:),       intent(in)    :: bc_type

      integer,                        intent(inout) :: n_tot,nnz

      integer,  dimension(4)                        :: statusX
      integer                                       :: c2,c3,c4,c5,c6,c7,c8
      integer                                       :: i,j,k,L,n
      integer                                       :: cnt

      continue

      c4=  Deta%rows
      c5=  Deta%cols
      c6=2*Deta%width+1
      c7=  Deta%width+1

      n_tot = 0
      nnz = 0
      do L = 1, nblock
        n_tot = n_tot + nx(L) * ny(L) * nz(L)
        nnz =   nnz + nx(L) * (ny(L)-2*c4) * nz(L) * c6    &
          + nx(L) * (      2*c4) * nz(L) * c5
      enddo
      write(*,*)'n_tot,nnz,eta  direction',n_tot,nnz

      allocate(ia(n_tot+1),stat=statusX(1))
      allocate(ja(  nnz  ),stat=statusX(2))
      allocate(ka(  nnz  ),stat=statusX(3))
      allocate( a(  nnz  ),stat=statusX(4))
      if(sum(statusX) > 0) then
        write(*,*)'trouble allocating CSR in CSR_Get_Pointers'
        write(*,*)'\Eta Convective Operator'
        stop
      endif

      nnz   = 1
      cnt   = 1
      ia(cnt) = nnz

      n_tot = 0
      do n=1,nblock
        c2=Deta%size(n)-Deta%cols
        c3=Deta%size(n)-Deta%rows
        do k=1,nz(n)
          ! left corner
          do j=1,c4               !  Sweep over rows
            c8 = (k-1)*nx(n)*ny(n)
            do i=1,nx(n)
              do L=1,c5               !  Sweep over columns
                ja(nnz) =  N_tot + c8 + (L-1)*nx(n) + i
                a(nnz) = Deta%block1(j,L,n)
                nnz  = nnz + 1
              end do
              cnt  = cnt + 1
              ia(cnt) = nnz
            end do
          end do

          ! interior
          do j=c4+1,c3            !  Sweep over rows
            c8 = (k-1)*nx(n)*ny(n)
            do i=1,nx(n)
              do L=1,c6               !  Sweep over columns
                c8 = (k-1)*nx(n)*ny(n) + (j-1)*nx(n) + (L-c7)*nx(n) + i
                ja(nnz) =   N_tot + c8
                a(nnz) = Deta%stencil(L,n)
                nnz  = nnz + 1
              end do
              cnt  = cnt + 1
              ia(cnt) = nnz
            end do
          end do

          ! right corner
          do j=1,c4              !  Sweep over rows
            c8 = (k-1)*nx(n)*ny(n)
            do i=1,nx(n)
              do L=1,c5               !  Sweep over columns
                ja(nnz) =   N_tot + c8 + (L + c2-1)*nx(n) + i
                a(nnz) = Deta%block2(j,L,n)
                nnz  = nnz + 1
              end do
              cnt  = cnt + 1
              ia(cnt) = nnz
            end do
          end do

        end do
        N_tot = N_tot + nx(n)*ny(n)*nz(n)
      end do
      nnz = nnz - 1    !  subtract off the last outer


      return
    end subroutine CSR_Get_Pointers_Eta

    !=============================================================================80

    subroutine CSR_Get_Pointers_Zeta(ia,ja,ka,a, Dzeta, bc_type, n_tot, nnz)

      type(diff_op),                  intent(in)    :: Dzeta

      integer,  dimension(:),         pointer       :: ia
      integer,  dimension(:),         pointer       :: ja
      integer,  dimension(:),         pointer       :: ka
      real(dp), dimension(:),         pointer       ::  a

      integer,  dimension(:,:),       intent(in)    :: bc_type

      integer,                        intent(inout) :: n_tot,nnz

      integer,  dimension(4)                        :: statusX
      integer                                       :: c2,c3,c4,c5,c6,c7,c8
      integer                                       :: i,j,k,L,n
      integer                                       :: cnt

      continue

      c4=  Dzeta%rows
      c5=  Dzeta%cols
      c6=2*Dzeta%width+1
      c7=  Dzeta%width+1

      n_tot = 0
      nnz = 0
      do L = 1, nblock
        n_tot = n_tot + nx(L) * ny(L) * nz(L)
        nnz =   nnz + nx(L) * ny(L) * (nz(L)-2*c4) * c6    &
          + nx(L) * ny(L) * (      2*c4) * c5
      enddo
      write(*,*)'n_tot,nnz,zeta direction',n_tot,nnz

      allocate(ia(n_tot+1),stat=statusX(1))
      allocate(ja(  nnz  ),stat=statusX(2))
      allocate(ka(  nnz  ),stat=statusX(3))
      allocate( a(  nnz  ),stat=statusX(4))
      if(sum(statusX) > 0) then
        write(*,*)'trouble allocating CSR in CSR_Get_Pointers'
        write(*,*)'\Zeta Convective Operator'
        stop
      endif

      nnz   = 1
      cnt   = 1
      ia(cnt) = nnz

      n_tot = 0
      do n=1,nblock
        c2=Dzeta%size(n)-Dzeta%cols
        c3=Dzeta%size(n)-Dzeta%rows
        ! left corner
        do k=1,c4               !  Sweep over rows
          do j=1,ny(n)
            do i=1,nx(n)
              do L=1,c5               !  Sweep over columns
                ja(nnz) = N_tot + (L-1)*nx(n)*ny(n) + (j-1)*nx(n) + i
                a(nnz) = Dzeta%block1(k,L,n)
                nnz  = nnz + 1
              end do
              cnt  = cnt + 1
              ia(cnt) = nnz
            end do
          end do
        end do

        ! interior
        do k=c4+1,c3            !  Sweep over rows
          do j=1,ny(n)
            do i=1,nx(n)
              do L=1,c6               !  Sweep over columns
                c8 = (k-1)*nx(n)*ny(n) + (j-1)*nx(n) + (L-c7)*nx(n)*ny(n) + i
                ja(nnz) = N_tot + c8
                a(nnz) = Dzeta%stencil(L,n)
                nnz  = nnz + 1
              end do
              cnt  = cnt + 1
              ia(cnt) = nnz
            end do
          end do
        end do

        ! right corner
        do k=1,c4              !  Sweep over rows
          do j=1,ny(n)
            do i=1,nx(n)
              c8 =  (c2-1)*nx(n)*ny(n) + (j-1)*nx(n) + i
              do L=1,c5               !  Sweep over columns
                ja(nnz) = N_tot + c8 + nx(n)*ny(n)*L
                a(nnz) = Dzeta%block2(k,L,n)
                nnz  = nnz + 1
              end do
              cnt  = cnt + 1
              ia(cnt) = nnz
            end do
          end do
        end do
        N_tot = N_tot + nx(n)*ny(n)*nz(n)
      end do
      nnz = nnz - 1    !  subtract off the last outer

      return 
    end subroutine CSR_Get_Pointers_Zeta

    !=============================================================================80

    subroutine CSR_Global_Pointers(CSR, bc_type,viscous_flow,cross_terms) 

      !  Subroutine forms pointers to ia, ja by combining classes of terms into
      !     a single data structure.

      use unary

      type(CSR_Matrix),               intent(inout) :: CSR

      integer,  dimension(:,:),       intent(in)    :: bc_type
      integer,                        intent(in)    :: viscous_flow,cross_terms

      integer,  dimension(10)                       :: statusX
      integer                                       :: i,j,k

      integer                                       :: ierr, cnt

      continue

      CSR%ntS = CSR%N_tot
      statusX = 0
      allocate(CSR%ndegr(CSR%nt1),stat=statusX(1))
      allocate(CSR%iwrk (CSR%nt2),stat=statusX(2))
      if( associated(CSR%ia ) ) deallocate(CSR%ia ) ; allocate(CSR%ia (CSR%n_tot+1),stat=statusX(3)) ;
      if( associated(CSR%iaS) ) deallocate(CSR%iaS) ; allocate(CSR%iaS(CSR%n_tot+1),stat=statusX(4)) ;

      if( associated(CSR%atmp) ) deallocate(CSR%atmp) ; allocate(CSR%atmp(1),stat=statusX(5)) ;

      if(sum(statusX) > 0) then
        write(*,*)'failure allocating ndegr in CSR_Global_Pointers.  Stopping'
        stop
      endif

      !=-=-=-=-=-=-=-=-=  Begin  =-=-=-=-=-=-=-=-=

      !     Calculate the number of nonzeros in the matrix sum A0 + A1  = A ;  Store in nnz and allocate ja array

      call aplbdg(CSR%nt0,CSR%nt1,CSR%ja0,CSR%ia0,CSR%ja1,CSR%ia1,CSR%ndegr,CSR%nnz,CSR%iwrk)

      if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
      if(sum(statusX) > 0) then
        write(*,*)'failure allocating ndegr in CSR_Global_Pointers.  Stopping'
        stop
      endif

      !     Now perform sum A0 + A1  = A ;  Store in nnz and allocate ja array


      call aplb1(CSR%nt0,CSR%nt1,0,CSR%a0,CSR%ja0,CSR%ia0,CSR%a1,CSR%ja1,CSR%ia1,CSR%atmp,CSR%ja,CSR%ia,CSR%nnz,ierr) 
      if(ierr /= 0) then ;  write(*,*)'ierr ',ierr,' : stopping after aplb1 terms 0-1' ; stop ; endif

        write(*,*)'nnz after combining 0-1',CSR%nnz
        if( associated(CSR%jaS) ) deallocate(CSR%jaS) ; allocate(CSR%jaS(CSR%nnz    ),stat=statusX(2)) ;

        CSR%iaS(1) = CSR%ia(1)
        do i = 1,CSR%n_tot
          CSR%iaS(i+1) = CSR%ia(i+1)
          do j = CSR%ia(i),CSR%ia(i+1)-1
            CSR%jaS(j) = CSR%ja(j)
          enddo
        enddo

        !     Calculate the number of nonzeros in the matrix sum A1 + A2  = A ;  Store in nnz and allocate ja array

        call aplbdg(CSR%ntS,CSR%nt2,CSR%jaS,CSR%iaS,CSR%ja2,CSR%ia2,CSR%ndegr,CSR%nnz,CSR%iwrk)

        if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
        if(sum(statusX) > 0) then
          write(*,*)'failure allocating ndegr in CSR_Global_Pointers.  Stopping'
          stop
        endif

        !     Now perform sum A1 + A2  = A ;  Store in nnz and allocate ja array


        call aplb1(CSR%ntS,CSR%nt2,0,CSR%aS,CSR%jaS,CSR%ia1,CSR%a2,CSR%ja2,CSR%ia2,CSR%atmp,CSR%ja,CSR%ia,CSR%nnz,ierr) 
        if(ierr /= 0) then ;  write(*,*)'ierr ',ierr,' : stopping after aplb1 terms 1-2' ; stop ; endif

          write(*,*)'nnz after combining 0-2',CSR%nnz
          if( associated(CSR%jaS) ) deallocate(CSR%jaS) ; allocate(CSR%jaS(CSR%nnz    ),stat=statusX(2)) ;

          CSR%iaS(1) = CSR%ia(1)
          do i = 1,CSR%n_tot
            CSR%iaS(i+1) = CSR%ia(i+1)
            do j = CSR%ia(i),CSR%ia(i+1)-1
              CSR%jaS(j) = CSR%ja(j)
            enddo
          enddo

          if(nz(1) > 1) then  !  Convective terms in z direction

            call aplbdg(CSR%ntS,CSR%nt3,CSR%jaS,CSR%iaS,CSR%ja3,CSR%ia3,CSR%ndegr,CSR%nnz,CSR%iwrk)

            if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
            if(sum(statusX) > 0) then
              write(*,*)'failure allocating ndegr in CSR_Global_Pointers.  Stopping'
              stop
            endif

            call aplb1(CSR%ntS,CSR%nt3,0,CSR%aS,CSR%jaS,CSR%iaS,CSR%a3,CSR%ja3,CSR%ia3,CSR%atmp,CSR%ja,CSR%ia,CSR%nnz,ierr) 
            if(ierr /= 0) then ;  write(*,*)'ierr ',ierr,' : stopping after aplb1 terms 0-3' ; stop ; endif

              write(*,*)'nnz after combining 0-3',CSR%nnz
              if( associated(CSR%jaS) ) deallocate(CSR%jaS) ; allocate(CSR%jaS(CSR%nnz    ),stat=statusX(2)) ;

              CSR%iaS(1) = CSR%ia(1)
              do i = 1,CSR%n_tot
                CSR%iaS(i+1) = CSR%ia(i+1)
                do j = CSR%ia(i),CSR%ia(i+1)-1
                  CSR%jaS(j) = CSR%ja(j)
                enddo
              enddo

            endif

            !=-=-=-=-=-=-=-=-=

            if(viscous_flow == 1) then

              call aplbdg(CSR%ntS,CSR%nt4,CSR%jaS,CSR%iaS,CSR%ja4,CSR%ia4,CSR%ndegr,CSR%nnz,CSR%iwrk)

              if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
              if(sum(statusX) > 0) then
                write(*,*)'failure allocating ndegr in CSR_Global_Pointers.  Stopping'
                stop
              endif


              call aplb1(CSR%ntS,CSR%nt4,0,CSR%aS,CSR%jaS,CSR%iaS,CSR%a4,CSR%ja4,CSR%ia4,CSR%atmp,CSR%ja,CSR%ia,CSR%nnz,ierr) 
              if(ierr /= 0) then ;  write(*,*)'ierr ',ierr,' : stopping after aplb1 terms 0-4' ; stop ; endif

                write(*,*)'nnz after combining 0-4',CSR%nnz
                if( associated(CSR%jaS) ) deallocate(CSR%jaS) ; allocate(CSR%jaS(CSR%nnz    ),stat=statusX(2)) ;

                CSR%iaS(1) = CSR%ia(1)
                do i = 1,CSR%n_tot
                  CSR%iaS(i+1) = CSR%ia(i+1)
                  do j = CSR%ia(i),CSR%ia(i+1)-1
                    CSR%jaS(j) = CSR%ja(j)
                  enddo
                enddo

                call aplbdg(CSR%ntS,CSR%nt5,CSR%jaS,CSR%iaS,CSR%ja5,CSR%ia5,CSR%ndegr,CSR%nnz,CSR%iwrk)

                if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
                if(sum(statusX)>0) then;write(*,*)'Allocation failure (ja) CSR_Global_Pointers. Stopping';stop;endif

                  call aplb1(CSR%ntS,CSR%nt5,0,CSR%aS,CSR%jaS,CSR%iaS,CSR%a5,CSR%ja5,CSR%ia5,CSR%atmp,CSR%ja,CSR%ia,CSR%nnz,ierr) 

                  if( associated(CSR%jaS) ) deallocate(CSR%jaS) ; allocate(CSR%jaS(CSR%nnz    ),stat=statusX(2)) ;
                  write(*,*)'nnz after combining 0-5',CSR%nnz

                  CSR%iaS(1) = CSR%ia(1)
                  do i = 1,CSR%n_tot
                    CSR%iaS(i+1) = CSR%ia(i+1)
                    do j = CSR%ia(i),CSR%ia(i+1)-1
                      CSR%jaS(j) = CSR%ja(j)
                    enddo
                  enddo

                  if(nz(1) > 1) then

                    call aplbdg(CSR%ntS,CSR%nt6,CSR%jaS,CSR%iaS,CSR%ja6,CSR%ia6,CSR%ndegr,CSR%nnz,CSR%iwrk)

                    if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
                    if(sum(statusX)>0) then;write(*,*)'Allocation failure (ja) CSR_Global_Pointers. Stopping';stop;endif

                      call aplb1(CSR%ntS,CSR%nt6,0,CSR%aS,CSR%jaS,CSR%iaS,CSR%a6,CSR%ja6,CSR%ia6,CSR%atmp,CSR%ja,CSR%ia,CSR%nnz,ierr) 

                      if( associated(CSR%jaS) ) deallocate(CSR%jaS) ; allocate(CSR%jaS(CSR%nnz    ),stat=statusX(2)) ;
                      write(*,*)'nnz after combining 0-6',CSR%nnz

                      CSR%iaS(1) = CSR%ia(1)
                      do i = 1,CSR%n_tot
                        CSR%iaS(i+1) = CSR%ia(i+1)
                        do j = CSR%ia(i),CSR%ia(i+1)-1
                          CSR%jaS(j) = CSR%ja(j)
                        enddo
                      enddo
                    endif

                  endif

                  !     Finished with accumulation of terms.

                  write(*,*)'n_tot, nnz',CSR%N_tot,CSR%nnz

                  if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz    ),stat=statusX(1)) ;
                  if(sum(statusX) > 0) then
                    write(*,*)'failure allocating ndegr in CSR_Global_Pointers.  Stopping'
                    stop
                  endif

                  CSR%ia(:) = CSR%iaS(:)
                  CSR%ja(:) = CSR%jaS(:)

                  CSR%ka0(1) = 1 ; 
                  CSR%ka1(1) = 1 ; CSR%ka2(1) = 1 ; if(nz(1) > 1) CSR%ka3(1) = 1 ;
                  if(viscous_flow == 1 ) then
                    CSR%ka4(1) = 1 ; CSR%ka5(1) = 1 ; if(nz(1) > 1) CSR%ka6(1) = 1 ;
                    if(cross_terms == 1 ) then
                      CSR%ka7(1) = 1 ; CSR%ka8(1) = 1 ; if(nz(1) > 1) CSR%ka9(1) = 1 ;
                    endif
                  endif
                  cnt = 1
                  do i = 1,CSR%n_tot
                    do j = CSR%ia(i),CSR%ia(i+1)-1

                      do k = CSR%ia0(i),CSR%ia0(i+1)-1
                        if(CSR%ja(j) == CSR%ja0(k)) CSR%ka0(k) = cnt
                      enddo

                      do k = CSR%ia1(i),CSR%ia1(i+1)-1
                        if(CSR%ja(j) == CSR%ja1(k)) CSR%ka1(k) = cnt
                      enddo
                      do k = CSR%ia2(i),CSR%ia2(i+1)-1
                        if(CSR%ja(j) == CSR%ja2(k)) CSR%ka2(k) = cnt
                      enddo
                      if(nz(1) > 1)then
                        do k = CSR%ia3(i),CSR%ia3(i+1)-1
                          if(CSR%ja(j) == CSR%ja3(k)) CSR%ka3(k) = cnt
                        enddo
                      endif
                      if(viscous_flow == 1 ) then
                        do k = CSR%ia4(i),CSR%ia4(i+1)-1
                          if(CSR%ja(j) == CSR%ja4(k)) CSR%ka4(k) = cnt
                        enddo
                        do k = CSR%ia5(i),CSR%ia5(i+1)-1
                          if(CSR%ja(j) == CSR%ja5(k)) CSR%ka5(k) = cnt
                        enddo
                        if(nz(1) > 1)then
                          do k = CSR%ia6(i),CSR%ia6(i+1)-1
                            if(CSR%ja(j) == CSR%ja6(k)) CSR%ka6(k) = cnt
                          enddo
                        endif
                        if(cross_terms == 1 ) then
                          do k = CSR%ia7(i),CSR%ia7(i+1)-1
                            if(CSR%ja(j) == CSR%ja7(k)) CSR%ka7(k) = cnt
                          enddo
                          do k = CSR%ia8(i),CSR%ia8(i+1)-1
                            if(CSR%ja(j) == CSR%ja8(k)) CSR%ka8(k) = cnt
                          enddo
                          if(nz(1) > 1)then
                            do k = CSR%ia9(i),CSR%ia9(i+1)-1
                              if(CSR%ja(j) == CSR%ja9(k)) CSR%ka9(k) = cnt
                            enddo
                          endif
                        endif
                      endif
                      cnt = cnt + 1
                    enddo
                  enddo


                  return
                end subroutine CSR_Global_Pointers

                !=============================================================================80

                subroutine Assemble_Jacobian_matrix(CSR, u,    du,    ddu,                &
                    uMean, duMdx,duMdy, duMdz,                     &
                    DivU, Gradrho, GradP, GradDivU, DivGradU,                   &
                    wrk,wrk1,wrk2,wrk3,                              &
                    metrics,                                         &
                    bc_type,temp_buffer,temp_buffer_all,ph,viscous_flow,       &
                    cross_terms,Dxi,Deta,Dzeta,D2xi,D2eta,D2zeta) 

                  type(CSR_Matrix),               intent(inout) :: CSR

                  type(physical_parameters),      intent(in)    :: ph
                  type(diff_op),                  intent(in)    :: Dxi, Deta, Dzeta
                  type(diff_op),                  intent(in)    :: D2xi,D2eta,D2zeta
                  type(temporary_buffer)                        :: temp_buffer
                  type(boundary_buffer),          intent(inout) :: temp_buffer_all
                  type(metric_fields),            intent(in)    :: metrics


                  complex(dp), dimension(:,:,:,:,:), intent(inout) :: u,     du    , ddu

                  real(dp),    dimension(:,:,:,:,:), intent(inout) :: uMean, duMdx, duMdy,  duMdz, wrk
                  real(dp),    dimension(:,:,:,:,:), intent(inout) :: Gradrho, GradP, GradDivU
                  real(dp),    dimension(:,:,:,:,:), intent(inout) :: DivGradU
                  real(dp),    dimension(:,:,:,:),   intent(inout) :: DivU, wrk1,wrk2,wrk3

                  integer,     dimension(:,:),       intent(in)    :: bc_type

                  integer,                           intent(in)    :: viscous_flow, cross_terms

                  complex(dp), dimension(nvar,nvar)                :: ca_loc
                  real(dp),    dimension(nvar,nvar)                ::  a_loc

                  real(dp)                                      :: Gam, Re, Ma, Pr, Gm1
                  real(dp)                                      :: mu, kap, con1, con2, norm, tmp
                  real(dp)                                      :: r1, r2, r4

                  integer                                       :: i,j,k,L,m
                  integer                                       :: cnt, m1,m2,n
                  integer                                       :: nnz,nnzt

                  integer,  dimension(3)                        :: statusX

                  logical                                       :: check_CSR = .false.
                  continue

                  Gam = ph%gamma ; Re = ph%Re  ; Ma = ph%Mach ; 
                  Pr  = ph%Pr    ; mu = ph%mu  ; kap= ph%kappa ;
                  Gm1 = Gam - 1
                  con1 = one  / Re ; con2 = Gam / Pr

                  !     allocate the matrix 'a'

                  if( associated(CSR%ia ) ) deallocate(CSR%ia ) ; allocate(CSR%ia (CSR%n_tot * nvar +    1),stat=statusX(1)) ;
                  if( associated(CSR%ja ) ) deallocate(CSR%ja ) ; allocate(CSR%ja (CSR%nnz   * nvar * nvar),stat=statusX(2)) ;
                  if( associated(CSR%a  ) ) deallocate(CSR%a  ) ; allocate(CSR%a  (CSR%nnz   * nvar * nvar),stat=statusX(3)) ;

                  CSR%ia(:) = 0
                  CSR%ja(:) = 0
                  CSR% a(:) = c_zero

                  !     Expand accounting for nvar /= 1  ;  Low dimensional data is still in iaS and jaS.

                  CSR%ia(1) = 1
                  do i = 1,CSR%n_tot
                    do j = 1,nvar
                      n = (i-1)*nvar + j
                      CSR%ia(n+1) = CSR%ia(n) + (CSR%iaS(i+1)-CSR%iaS(i))*nvar
                    enddo
                  enddo

                  nnz = 1
                  do i = 1,CSR%n_tot
                    do m2 = 1,nvar
                      do n = CSR%iaS(i),CSR%iaS(i+1) - 1
                        do m1 = 1,nvar

                          nnzt = (CSR%iaS(i)-1)*nvar*nvar + (m2-1)*(CSR%iaS(i+1)-CSR%iaS(i))* nvar + (n-CSR%iaS(i))*nvar + m1
                          if(nnzt-nnz /= 0)then ; write(*,*)'shuffle error in expanding ia, ja: stopping';stop;endif; 

                            CSR%ja(nnz) = (CSR%jaS(n)-1)*nvar +  m1
                            nnz  =  nnz + 1
                            cnt  =  cnt + 1
                          enddo
                        enddo
                      enddo
                    enddo
                    CSR%n_tot = CSR%n_tot * nvar
                    CSR%nnz   = CSR%nnz   * nvar

                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !                             Time Terms
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                    cnt  = 1
                    do L = 1,nblock
                      do k = 1,nz(L)
                        do j = 1,ny(L)
                          do i = 1,nx(L)

                            ca_loc(:,:) = c_zero

                            ca_loc(1,1) = ph%freeK * c_eye
                            ca_loc(2,2) = ph%freeK * c_eye
                            ca_loc(3,3) = ph%freeK * c_eye
                            ca_loc(4,4) = ph%freeK * c_eye
                            ca_loc(5,5) = ph%freeK * c_eye

                            do m2 = 1,nvar
                              do n = CSR%ia0(cnt),CSR%ia0(cnt+1)-1
                                do m1 = 1,nvar
                                  nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                    + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                    + (CSR%ka0(n)-CSR%iaS(cnt))*nvar + m1
                                  CSR%a(nnzt) = CSR%a(nnzt) + ca_loc(m2,m1)
                                enddo
                              enddo
                            enddo

                            cnt = cnt + 1

                          enddo
                        enddo
                      enddo
                    enddo

                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !                             Inviscid Terms
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                    cnt  = 1
                    do L = 1,nblock
                      do k = 1,nz(L)
                        do j = 1,ny(L)
                          do i = 1,nx(L)

                            a_loc(:,:) = zero

                            a_loc(1,1) =     uMean(i,j,k,2,L)  *  metrics%dxidx(i,j,k,L)
                            a_loc(2,2) =     uMean(i,j,k,2,L)  *  metrics%dxidx(i,j,k,L)
                            a_loc(3,3) =     uMean(i,j,k,2,L)  *  metrics%dxidx(i,j,k,L)
                            a_loc(4,4) =     uMean(i,j,k,2,L)  *  metrics%dxidx(i,j,k,L)
                            a_loc(5,5) =     uMean(i,j,k,2,L)  *  metrics%dxidx(i,j,k,L)

                            a_loc(1,2) =     uMean(i,j,k,1,L)  *  metrics%dxidx(i,j,k,L)
                            a_loc(2,5) = one/uMean(i,j,k,1,L)  *  metrics%dxidx(i,j,k,L)
                            a_loc(5,2) = gam*uMean(i,j,k,5,L)  *  metrics%dxidx(i,j,k,L)

                            do m2 = 1,nvar
                              do n = CSR%ia1(cnt),CSR%ia1(cnt+1)-1
                                do m1 = 1,nvar
                                  nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                    + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                    + (CSR%ka1(n)-CSR%iaS(cnt))*nvar + m1
                                  CSR%a(nnzt) = CSR%a(nnzt)  + a_loc(m2,m1) * CSR%a1(n)
                                enddo
                              enddo
                            enddo

                            cnt = cnt + 1

                          enddo
                        enddo
                      enddo
                    enddo

                    cnt  = 1
                    do L = 1,nblock
                      do k = 1,nz(L)
                        do j = 1,ny(L)
                          do i = 1,nx(L)

                            a_loc(:,:) = zero

                            a_loc(1,1) =     uMean(i,j,k,3,L)  *  metrics%detady(i,j,k,L)
                            a_loc(2,2) =     uMean(i,j,k,3,L)  *  metrics%detady(i,j,k,L)
                            a_loc(3,3) =     uMean(i,j,k,3,L)  *  metrics%detady(i,j,k,L)
                            a_loc(4,4) =     uMean(i,j,k,3,L)  *  metrics%detady(i,j,k,L)
                            a_loc(5,5) =     uMean(i,j,k,3,L)  *  metrics%detady(i,j,k,L)

                            a_loc(1,3) =     uMean(i,j,k,1,L)  *  metrics%detady(i,j,k,L)
                            a_loc(3,5) = one/uMean(i,j,k,1,L)  *  metrics%detady(i,j,k,L)
                            a_loc(5,3) = gam*uMean(i,j,k,5,L)  *  metrics%detady(i,j,k,L)

                            do m2 = 1,nvar
                              do n = CSR%ia2(cnt),CSR%ia2(cnt+1)-1
                                do m1 = 1,nvar
                                  nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                    + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                    + (CSR%ka2(n)-CSR%iaS(cnt))*nvar + m1
                                  CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a2(n)
                                enddo
                              enddo
                            enddo

                            cnt = cnt + 1


                          enddo
                        enddo
                      enddo
                    enddo

                    if(nz(1) > 1) then

                      cnt  = 1
                      do L = 1,nblock
                        do k = 1,nz(L)
                          do j = 1,ny(L)
                            do i = 1,nx(L)

                              a_loc(:,:) = zero

                              a_loc(1,1) =     uMean(i,j,k,3,L)  *  metrics%dzetadz(i,j,k,L)
                              a_loc(2,2) =     uMean(i,j,k,3,L)  *  metrics%dzetadz(i,j,k,L)
                              a_loc(3,3) =     uMean(i,j,k,3,L)  *  metrics%dzetadz(i,j,k,L)
                              a_loc(4,4) =     uMean(i,j,k,3,L)  *  metrics%dzetadz(i,j,k,L)
                              a_loc(5,5) =     uMean(i,j,k,3,L)  *  metrics%dzetadz(i,j,k,L)

                              a_loc(1,4) =     uMean(i,j,k,1,L)  *  metrics%dzetadz(i,j,k,L)
                              a_loc(4,5) = one/uMean(i,j,k,1,L)  *  metrics%dzetadz(i,j,k,L)
                              a_loc(5,4) = gam*uMean(i,j,k,5,L)  *  metrics%dzetadz(i,j,k,L)

                              do m2 = 1,nvar
                                do n = CSR%ia3(cnt),CSR%ia3(cnt+1)-1
                                  do m1 = 1,nvar
                                    nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                      + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                      + (CSR%ka3(n)-CSR%iaS(cnt))*nvar + m1
                                    CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a3(n)
                                  enddo
                                enddo
                              enddo

                              cnt = cnt + 1

                            enddo
                          enddo
                        enddo
                      enddo
                    endif

                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !                               Viscous terms
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                    if ((viscous_flow/=0)) then

                      !------------------------------------------------------------------------------
                      !       calculate mean flow gradients dUdx, dUdy , dUdz
                      !------------------------------------------------------------------------------

                      call diff_xi  (Dxi,  uMean, temp_buffer_all, wrk, bc_type)
                      duMdx(:,:,:,:,:) = zero
                      do L = 1,nblock
                        do m = 1,nvar
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)
                                duMdx(i,j,k,m,L) = duMdx(i,j,k,m,L) + wrk(i,j,k,m,L)*metrics%dxidx  (i,j,k,L)
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                      call diff_eta (Deta, uMean, temp_buffer_all, wrk, bc_type)
                      duMdy(:,:,:,:,:) = zero
                      do L = 1,nblock
                        do m = 1,nvar
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)
                                duMdy(i,j,k,m,L) = duMdy(i,j,k,m,L) + wrk(i,j,k,m,L)*metrics%detady (i,j,k,L)
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                      if(nz(1) > 1) then
                        call diff_zeta(Dzeta,uMean, temp_buffer_all, wrk, bc_type)
                        duMdz(:,:,:,:,:) = zero
                        do L = 1,nblock
                          do m = 1,nvar
                            do k = 1,nz(L)
                              do j = 1,ny(L)
                                do i = 1,nx(L)
                                  duMdz(i,j,k,m,L) = duMdz(i,j,k,m,L) + wrk(i,j,k,m,L)*metrics%dzetadz(i,j,k,L)
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo
                      endif

                      !        Div U , Grad P, and Grad rho

                      do L = 1,nblock
                        do k = 1,nz(L)
                          do j = 1,ny(L)
                            do i = 1,nx(L)

                              DivU   (i,j,k,L)   = + duMdx(i,j,k,2,L) + duMdy(i,j,k,3,L) + duMdz(i,j,k,4,L)

                              Gradrho(i,j,k,1,L) = + duMdx(i,j,k,1,L)
                              Gradrho(i,j,k,2,L) = + duMdy(i,j,k,1,L)
                              Gradrho(i,j,k,3,L) = + duMdz(i,j,k,1,L)
                              GradP  (i,j,k,1,L) = + duMdx(i,j,k,5,L)
                              GradP  (i,j,k,2,L) = + duMdy(i,j,k,5,L)
                              GradP  (i,j,k,3,L) = + duMdz(i,j,k,5,L)
                            enddo
                          enddo
                        enddo
                      enddo

                      call diff_xi  (Dxi,   DivU, temp_buffer, wrk1, bc_type)
                      call diff_eta (Deta,  DivU, temp_buffer, wrk2, bc_type)
                      call diff_zeta(Dzeta, DivU, temp_buffer, wrk3, bc_type)

                      do L = 1,nblock
                        do m = 1,nvar
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)
                                GradDivU(i,j,k,1,L) = wrk1(i,j,k,L)*metrics%dxidx  (i,j,k,L)
                                GradDivU(i,j,k,2,L) = wrk2(i,j,k,L)*metrics%detady (i,j,k,L)
                                GradDivU(i,j,k,3,L) = wrk3(i,j,k,L)*metrics%dzetadz(i,j,k,L)
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo

                      DivGradU(:,:,:,:,:) = zero

                      call diff_xi  (D2xi,  uMean, temp_buffer_all, wrk, bc_type)
                      do L = 1,nblock
                        do m = 1,nvar
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)
                                DivGradU(i,j,k,m,L) = DivGradU(i,j,k,m,L) + wrk(i,j,k,m,L)*metrics%dxidx  (i,j,k,L)**2
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                      call diff_eta (D2eta, uMean, temp_buffer_all, wrk, bc_type)
                      do L = 1,nblock
                        do m = 1,nvar
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)
                                DivGradU(i,j,k,m,L) = DivGradU(i,j,k,m,L) + wrk(i,j,k,m,L)*metrics%detady (i,j,k,L)**2
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                      call diff_zeta(D2zeta,uMean, temp_buffer_all, wrk, bc_type)
                      do L = 1,nblock
                        do m = 1,nvar
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)
                                DivGradU(i,j,k,m,L) = DivGradU(i,j,k,m,L) + wrk(i,j,k,m,L)*metrics%dzetadz(i,j,k,L)**2
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo

                      !       lower order terms

                      cnt  = 1
                      do L = 1,nblock
                        do k = 1,nz(L)
                          do j = 1,ny(L)
                            do i = 1,nx(L)

                              r2 = one / (uMean(i,j,k,1,L) * uMean(i,j,k,1,L))

                              a_loc(1,1) =      DivU(i,j,k,  L)
                              a_loc(1,2) =     duMdx(i,j,k,1,L)
                              a_loc(1,3) =     duMdy(i,j,k,1,L)
                              a_loc(1,4) =     duMdz(i,j,k,1,L)
                              a_loc(1,5) = zero

                              a_loc(2,1) = -r2*duMdx(i,j,k,5,L)
                              a_loc(2,2) =     duMdx(i,j,k,2,L)
                              a_loc(2,3) =     duMdy(i,j,k,2,L)
                              a_loc(2,4) =     duMdz(i,j,k,2,L)
                              a_loc(2,5) = zero

                              a_loc(3,1) = -r2*duMdy(i,j,k,5,L)
                              a_loc(3,2) =     duMdx(i,j,k,3,L)
                              a_loc(3,3) =     duMdy(i,j,k,3,L)
                              a_loc(3,4) =     duMdz(i,j,k,3,L)
                              a_loc(3,5) = zero

                              a_loc(4,1) = -r2*duMdz(i,j,k,5,L)
                              a_loc(4,2) =     duMdx(i,j,k,4,L)
                              a_loc(4,3) =     duMdy(i,j,k,4,L)
                              a_loc(4,4) =     duMdz(i,j,k,4,L)
                              a_loc(4,5) = zero

                              a_loc(5,1) = zero
                              a_loc(5,2) =     duMdx(i,j,k,5,L)
                              a_loc(5,3) =     duMdy(i,j,k,5,L)
                              a_loc(5,4) =     duMdz(i,j,k,5,L)
                              a_loc(5,5) = gam* DivU(i,j,k,  L)

                              do m2 = 1,nvar
                                do n = CSR%ia0(cnt),CSR%ia0(cnt+1)-1
                                  do m1 = 1,nvar
                                    nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                      + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                      + (CSR%ka0(n)-CSR%iaS(cnt))*nvar + m1
                                    CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1)
                                  enddo
                                enddo
                              enddo

                              cnt = cnt + 1

                            enddo
                          enddo
                        enddo
                      enddo

                      !------------------------------------------------------------------------------
                      !       Lower order convective terms arrising from viscosity (and nonuniform mean flow) 
                      !------------------------------------------------------------------------------

                      cnt  = 1
                      do L = 1,nblock
                        do k = 1,nz(L)
                          do j = 1,ny(L)
                            do i = 1,nx(L)

                              tmp = con1 * metrics%dxidx(i,j,k,L)

                              a_loc(:,:) = zero

                              a_loc(5,1) =  tmp * (                                                     &
                                + con2 * two * kap * (+ one * uMean(i,j,k,1,L)*duMdx(i,j,k,5,L)       &
                                - two * uMean(i,j,k,1,L)*duMdx(i,j,k,1,L))/     &
                                uMean(i,j,k,1,L)**3)

                              a_loc(5,2) =  tmp * (                                     &
                                - four_by_three * gm1 * mu * (+ two * duMdx(i,j,k,2,L)                &
                                - one * duMdy(i,j,k,3,L)                &
                                - one * duMdz(i,j,k,4,L))  )

                              a_loc(5,3) =  tmp * (                                     &
                                -           two * gm1 * mu * (+ one * duMdy(i,j,k,2,L)                &
                                + one * duMdx(i,j,k,3,L))  )

                              a_loc(5,4) =  tmp * (                                     &
                                -           two * gm1 * mu * (+ one * duMdz(i,j,k,2,L)                &
                                + one * duMdx(i,j,k,4,L))  )

                              a_loc(5,5) =  tmp * (                                     &
                                + con2 * two * kap * (  duMdx(i,j,k,1,L))/ uMean(i,j,k,1,L)**2 )

                              do m2 = 1,nvar
                                do n = CSR%ia1(cnt),CSR%ia1(cnt+1)-1
                                  do m1 = 1,nvar
                                    nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                      + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                      + (CSR%ka1(n)-CSR%iaS(cnt))*nvar + m1
                                    CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a1(n)
                                  enddo
                                enddo
                              enddo

                              cnt = cnt + 1

                            enddo
                          enddo
                        enddo
                      enddo

                      cnt  = 1
                      do L = 1,nblock
                        do k = 1,nz(L)
                          do j = 1,ny(L)
                            do i = 1,nx(L)

                              tmp = con1 * metrics%detady(i,j,k,L)

                              a_loc(:,:) = zero

                              a_loc(5,1) = +  tmp * ( + con2 * two * kap * (+ one * uMean(i,j,k,1,L)*duMdy(i,j,k,5,L)   &
                                - two * uMean(i,j,k,1,L)*duMdy(i,j,k,1,L))/ &
                                uMean(i,j,k,1,L)**3      )

                              a_loc(5,2) = +  tmp * ( - two * gm1 * mu * (+ one * duMdy(i,j,k,2,L)             &
                                + one * duMdx(i,j,k,3,L)))

                              a_loc(5,3) = +  tmp * ( - four_by_three * gm1 * mu * (- one * duMdx(i,j,k,2,L)   &
                                + two * duMdy(i,j,k,3,L)   &
                                - one * duMdz(i,j,k,4,L)))

                              a_loc(5,4) = +  tmp * ( - two * gm1 * mu * (+ one * duMdz(i,j,k,3,L)      &
                                + one * duMdy(i,j,k,4,L)))

                              a_loc(5,5) = +  tmp * ( + con2 * two * kap * ( duMdy(i,j,k,1,L))/ uMean(i,j,k,1,L)**2 )

                              do m2 = 1,nvar
                                do n = CSR%ia2(cnt),CSR%ia2(cnt+1)-1
                                  do m1 = 1,nvar
                                    nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                      + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                      + (CSR%ka2(n)-CSR%iaS(cnt))*nvar + m1
                                    !                     CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a2(n)
                                  enddo
                                enddo
                              enddo

                              cnt = cnt + 1

                            enddo
                          enddo
                        enddo
                      enddo

                      if(nz(1) > 1) then

                        cnt  = 1
                        do L = 1,nblock
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)

                                a_loc(:,:) = zero

                                tmp = con1 * metrics%dzetadz(i,j,k,L)

                                a_loc(5,1) =                    +  tmp  * (                               &
                                  +   con2 * two * kap * (+ one * uMean(i,j,k,1,L)*duMdy(i,j,k,5,L)     &
                                  - two * uMean(i,j,k,1,L)*duMdy(i,j,k,1,L))/   &
                                  uMean(i,j,k,1,L)**3   )

                                a_loc(5,2) =                    +  tmp  * (                               &
                                  -          two * gm1 * mu * (+ one * duMdz(i,j,k,2,L)                 &
                                  + one * duMdx(i,j,k,4,L)))

                                a_loc(5,3) =                    +  tmp  * (                               &
                                  -          two * gm1 * mu * (+ one * duMdz(i,j,k,3,L)                 &
                                  + one * duMdy(i,j,k,4,L)))

                                a_loc(5,4) =                    +  tmp  * (                               &
                                  -four_by_three * gm1 * mu * (- one * duMdx(i,j,k,2,L)                 &
                                  - one * duMdy(i,j,k,3,L)                 &
                                  + two * duMdz(i,j,k,4,L)))
                                a_loc(5,5) =                    +  tmp  * (                               &

                                  +   con2 * two * kap * (                         duMdz(i,j,k,1,L))/   &
                                    uMean(i,j,k,1,L)**2   )

                                  do m2 = 1,nvar
                                    do n = CSR%ia3(cnt),CSR%ia3(cnt+1)-1
                                      do m1 = 1,nvar
                                        nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                          + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                          + (CSR%ka3(n)-CSR%iaS(cnt))*nvar + m1
                                        !                       CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a3(n)
                                      enddo
                                    enddo
                                  enddo

                                  cnt = cnt + 1

                                enddo
                              enddo
                            enddo
                          enddo
                        endif

                        !----------------------------------__---__-------------------------------------
                        !       Second order viscous terms \/ * \/ \phi
                        !------------------------------------------------------------------------------


                        !       D^2 \phi / D xi^2
                        cnt  = 1
                        do L = 1,nblock
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)

                                tmp = con1 * metrics%dxidx(i,j,k,L)*metrics%dxidx(i,j,k,L)

                                a_loc(:,:) = zero

                                a_loc(2,2) = - (four_by_three * mu / uMean(i,j,k,1,L)    ) * tmp
                                a_loc(3,3) = - (                mu / uMean(i,j,k,1,L)    ) * tmp
                                a_loc(4,4) = - (                mu / uMean(i,j,k,1,L)    ) * tmp
                                a_loc(5,5) = - (        con2  *kap / uMean(i,j,k,1,L)    ) * tmp

                                a_loc(5,1) = + (        con2  *kap * uMean(i,j,k,5,L)    &
                                  / uMean(i,j,k,1,L)**2 ) * tmp

                                do m2 = 1,nvar
                                  do n = CSR%ia4(cnt),CSR%ia4(cnt+1)-1
                                    do m1 = 1,nvar
                                      nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                        + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                        + (CSR%ka4(n)-CSR%iaS(cnt))*nvar + m1
                                      CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a4(n)
                                    enddo
                                  enddo
                                enddo

                                cnt = cnt + 1

                              enddo
                            enddo
                          enddo
                        enddo

                        !       D^2 \phi / D xi eta
                        if(cross_terms == 1) then
                          cnt  = 1
                          do L = 1,nblock
                            do k = 1,nz(L)
                              do j = 1,ny(L)
                                do i = 1,nx(L)

                                  tmp = con1 * metrics%dxidx(i,j,k,L)*metrics%detady(i,j,k,L) * third

                                  a_loc(:,:) = zero

                                  a_loc(2,3) = - ( mu / uMean(i,j,k,1,L) ) * tmp
                                  a_loc(3,2) = - ( mu / uMean(i,j,k,1,L) ) * tmp

                                  do m2 = 1,nvar
                                    do n = CSR%ia7(cnt),CSR%ia7(cnt+1)-1
                                      do m1 = 1,nvar
                                        nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                          + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                          + (CSR%ka7(n)-CSR%iaS(cnt))*nvar + m1
                                        CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a7(n)
                                      enddo
                                    enddo
                                  enddo

                                  cnt = cnt + 1

                                enddo
                              enddo
                            enddo
                          enddo
                        endif

                        !       D^2 \phi / D eta^2
                        cnt  = 1
                        do L = 1,nblock
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)

                                tmp = con1 * metrics%detady(i,j,k,L)*metrics%detady(i,j,k,L)

                                a_loc(:,:) = zero

                                a_loc(2,2) = - (                mu / uMean(i,j,k,1,L)    ) * tmp
                                a_loc(3,3) = - (four_by_three * mu / uMean(i,j,k,1,L)    ) * tmp
                                a_loc(4,4) = - (                mu / uMean(i,j,k,1,L)    ) * tmp
                                a_loc(5,5) = - (        con2  *kap / uMean(i,j,k,1,L)    ) * tmp

                                a_loc(5,1) = + (        con2  *kap * uMean(i,j,k,5,L)    &
                                  / uMean(i,j,k,1,L)**2 ) * tmp

                                do m2 = 1,nvar
                                  do n = CSR%ia5(cnt),CSR%ia5(cnt+1)-1
                                    do m1 = 1,nvar
                                      nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                        + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                        + (CSR%ka5(n)-CSR%iaS(cnt))*nvar + m1
                                      CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a5(n)
                                    enddo
                                  enddo
                                enddo

                                cnt = cnt + 1

                              enddo
                            enddo
                          enddo
                        enddo

                        if(nz(1) > 1) then     !     3D path

                          if(cross_terms == 1) then
                            !           D^2 \phi / D xi zeta
                            cnt  = 1
                            do L = 1,nblock
                              do k = 1,nz(L)
                                do j = 1,ny(L)
                                  do i = 1,nx(L)

                                    tmp = con1 * metrics%dxidx(i,j,k,L)*metrics%dzetadz(i,j,k,L) * third

                                    a_loc(:,:) = zero

                                    a_loc(2,4) = - ( mu / uMean(i,j,k,1,L) ) * tmp
                                    a_loc(4,2) = - ( mu / uMean(i,j,k,1,L) ) * tmp

                                    do m2 = 1,nvar
                                      do n = CSR%ia8(cnt),CSR%ia8(cnt+1)-1
                                        do m1 = 1,nvar
                                          nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                            + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                            + (CSR%ka8(n)-CSR%iaS(cnt))*nvar + m1
                                          CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a8(n)
                                        enddo
                                      enddo
                                    enddo

                                    cnt = cnt + 1

                                  enddo
                                enddo
                              enddo
                            enddo


                            !           D^2 \phi / D eta zeta
                            cnt  = 1
                            do L = 1,nblock
                              do k = 1,nz(L)
                                do j = 1,ny(L)
                                  do i = 1,nx(L)

                                    tmp = con1 * metrics%detady(i,j,k,L)*metrics%dzetadz(i,j,k,L) * third

                                    a_loc(:,:) = zero

                                    a_loc(3,4) = - ( mu / uMean(i,j,k,1,L) ) * tmp
                                    a_loc(4,3) = - ( mu / uMean(i,j,k,1,L) ) * tmp

                                    do m2 = 1,nvar
                                      do n = CSR%ia9(cnt),CSR%ia9(cnt+1)-1
                                        do m1 = 1,nvar
                                          nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                            + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                            + (CSR%ka9(n)-CSR%iaS(cnt))*nvar + m1
                                          CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a9(n)
                                        enddo
                                      enddo
                                    enddo

                                    cnt = cnt + 1
                                  enddo
                                enddo
                              enddo
                            enddo
                          endif

                          !         D^2 \phi / D zeta^2
                          cnt  = 1
                          do L = 1,nblock
                            do k = 1,nz(L)
                              do j = 1,ny(L)
                                do i = 1,nx(L)

                                  tmp = con1 * metrics%dzetadz(i,j,k,L)*metrics%dzetadz(i,j,k,L)

                                  a_loc(:,:) = zero

                                  a_loc(2,2) = - (                mu / uMean(i,j,k,1,L)    ) * tmp
                                  a_loc(3,3) = - (                mu / uMean(i,j,k,1,L)    ) * tmp
                                  a_loc(4,4) = - (four_by_three * mu / uMean(i,j,k,1,L)    ) * tmp
                                  a_loc(5,5) = - (        con2  *kap / uMean(i,j,k,1,L)    ) * tmp

                                  a_loc(5,1) = + (        con2  *kap * uMean(i,j,k,5,L)    &
                                    / uMean(i,j,k,1,L)**2 ) * tmp

                                  do m2 = 1,nvar
                                    do n = CSR%ia6(cnt),CSR%ia6(cnt+1)-1
                                      do m1 = 1,nvar
                                        nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                          + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                          + (CSR%ka6(n)-CSR%iaS(cnt))*nvar + m1
                                        CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a6(n)
                                      enddo
                                    enddo
                                  enddo

                                  cnt = cnt + 1

                                enddo
                              enddo
                            enddo
                          enddo

                        endif


                        !     lower order viscous terms (C^* \phi)

                        cnt  = 1
                        do L = 1,nblock
                          do k = 1,nz(L)
                            do j = 1,ny(L)
                              do i = 1,nx(L)

                                r1 = one / uMean(i,j,k,1,L)
                                r2 = r1*r1
                                r4 = r2*r2

                                tmp = con1 

                                a_loc(:,:) = zero

                                a_loc(2,1) =   con1 * mu * r2 * ( DivGradU(i,j,k,2,L) + third*GradDivU(i,j,k,1,L))
                                a_loc(3,1) =   con1 * mu * r2 * ( DivGradU(i,j,k,3,L) + third*GradDivU(i,j,k,2,L))
                                a_loc(4,1) =   con1 * mu * r2 * ( DivGradU(i,j,k,4,L) + third*GradDivU(i,j,k,3,L))

                                a_loc(5,1) =  + con2 * kap* r4 * ( + six*uMean(i,j,k,5,L)    * (                                                   &
                                  + Gradrho(i,j,k,1,L)*Gradrho(i,j,k,1,L)     &
                                  + Gradrho(i,j,k,2,L)*Gradrho(i,j,k,2,L)     &
                                  + Gradrho(i,j,k,3,L)*Gradrho(i,j,k,3,L) )   &
                                  + one*uMean(i,j,k,1,L)**2 * DivGradU(i,j,k,5,L)                 &
                                  - two*uMean(i,j,k,1,L) *    ( two * (                           &
                                  + GradP(i,j,k,1,L)*Gradrho(i,j,k,1,L)       &
                                  + GradP(i,j,k,2,L)*Gradrho(i,j,k,2,L)       &
                                  + GradP(i,j,k,3,L)*Gradrho(i,j,k,3,L) ))    &
                                  + one*uMean(i,j,k,5,L)**2 * DivGradU(i,j,k,1,L))
                                a_loc(5,5) =  + con2 * kap* r4 * ( - two*uMean(i,j,k,1,L) * (                                                      &
                                  + Gradrho(i,j,k,1,L)*Gradrho(i,j,k,1,L)     &
                                  + Gradrho(i,j,k,2,L)*Gradrho(i,j,k,2,L)     &
                                  + Gradrho(i,j,k,3,L)*Gradrho(i,j,k,3,L) )   &
                                  + one*uMean(i,j,k,1,L)**2 * DivGradU(i,j,k,1,L))

                                do m2 = 1,nvar
                                  do n = CSR%ia0(cnt),CSR%ia0(cnt+1)-1
                                    do m1 = 1,nvar
                                      nnzt  = (CSR%iaS(cnt)-1)*nvar*nvar                      &
                                        + (m2-1)*(CSR%iaS(cnt+1)-CSR%iaS(cnt)) * nvar     &
                                        + (CSR%ka0(n)-CSR%iaS(cnt))*nvar + m1
                                      CSR%a(nnzt) = CSR%a(nnzt) + a_loc(m2,m1) * CSR%a0(n)
                                    enddo
                                  enddo
                                enddo

                                cnt = cnt + 1

                              enddo
                            enddo
                          enddo
                        enddo

                      endif    !  end viscous terms (conditional)

                      if(check_CSR) then
                        nnz = 1
                        cnt = 0
                        do i = 1,CSR%n_tot
                          do j = CSR%ia(i),CSR%ia(i+1)-1
                            if(abs(CSR%a(nnz)) <= 1.0e-8)cnt = cnt + 1
                            nnz = nnz + 1 
                          enddo
                        enddo
                        if(cnt /= 0) then 
                          write(*,*)'Error in pointers for Amat (a) matris  :  ', nnz,cnt   
                          write(*,*)'Stopping in Assemble_Jacobian_Matrix'
                          !          write(180,*)'Here comes ja and a'
                          !          do i = 1,CSR%n_tot
                          !               write(180,'(i5,7x,32(i4  ,1x))')i,(CSR%ja(j),j=CSR%ia(i),CSR%ia(i+1)-1)
                          !               write(180,'(i5,7x,32(f4.1,1x))')i,(CSR% a(j),j=CSR%ia(i),CSR%ia(i+1)-1)
                          !          enddo
                          stop
                        endif
                      endif
                      180 continue

                      return
                    end subroutine Assemble_Jacobian_matrix

                    !=============================================================================80

                    subroutine testing(CSR, x,y,z, uu, u, du, ddu, res, cwrk, v1,v2,v3,r1,w1,w2, tmp_b_all,  &
                        bc_type,viscous_flow, order,                  & 
                        Dxi,Deta,Dzeta,D2xi,D2eta,D2zeta)


                      implicit none

                      type(diff_op),                  intent(in)    :: Dxi, Deta, Dzeta
                      type(diff_op),                  intent(in)    :: D2xi,D2eta,D2zeta

                      type(CSR_Matrix),               intent(inout) :: CSR
                      type(boundary_buffer)                         :: tmp_b_all

                      complex(dp), dimension(:,:,:,:,:), intent(inout) :: uu, u, du, ddu, res, cwrk
                      complex(dp), dimension(:),         intent(inout) :: v1,v2,v3,r1
                      complex(dp), dimension(:),         intent(inout) :: w1,w2

                      real(dp),    dimension(:,:,:,:),   intent(inout) :: x,y,z

                      integer,  dimension(:,:),       intent(in)    :: bc_type
                      integer,                        intent(in)    :: viscous_flow, order

                      integer                                       :: n_tot,cnt
                      integer                                       :: i,j,k,L,n

                      integer                                       :: p

                      real(dp), parameter                           :: eps = 1.0e-10
                      real(dp)                                      :: norm
                      character(len=1)                              :: direction = 'x'
                      character(len=10)                             :: output_file = 'plotter.dat'
                      character(len=10)                             :: title = 'error'

                      logical                                       :: check_design_order = .true.

                      continue

                      p = 1
                      select case(order)
                      case(21) ; p = 1
                      case(42) ; p = 2  ; case(43) ; p = 3
                      case(63) ; p = 3  ; case(65) ; p = 5
                      case(84) ; p = 4
                      end select

                      if(check_design_order)then
                        do cnt = 1,3
                          if(cnt == 1)direction='x'
                          if(cnt == 2)direction='y'
                          if(cnt == 3)direction='z'
                          n_tot =  0
                          do n=1,nblock
                            do k=1,nz(n)
                              do j=1,ny(n)
                                do i=1,nx(n)
                                  do L=1,nvar
                                    select case(direction)
                                    case('x')
                                      uu(i,j,k,L,n) =                x(i,j,k,n) ** (p-0) * c_one
                                      du(i,j,k,L,n) =  (p-0)       * x(i,j,k,n) ** (p-1) * c_one
                                      ddu(i,j,k,L,n) =  (p-0)*(p-1) * x(i,j,k,n) ** (p-2) * c_one
                                    case('y')
                                      uu(i,j,k,L,n) =                y(i,j,k,n) ** (p-0) * c_one
                                      du(i,j,k,L,n) =  (p-0)       * y(i,j,k,n) ** (p-1) * c_one
                                      ddu(i,j,k,L,n) =  (p-0)*(p-1) * y(i,j,k,n) ** (p-2) * c_one
                                    case('z')
                                      uu(i,j,k,L,n) =                z(i,j,k,n) ** (p-0) * c_one
                                      du(i,j,k,L,n) =  (p-0)       * z(i,j,k,n) ** (p-1) * c_one
                                      ddu(i,j,k,L,n) =  (p-0)*(p-1) * z(i,j,k,n) ** (p-2) * c_one
                                    end select 
                                    if(p==1) ddu(i,j,k,L,n) =  0.0_dp * c_one
                                  enddo
                                enddo
                              enddo
                            enddo
                            n_tot = n_tot + nx(n)*ny(n)*nz(n)
                          enddo

                          n_tot =  0
                          do n=1,nblock
                            do k=1,nz(n)
                              do j=1,ny(n)
                                do i=1,nx(n)
                                  n_tot = n_tot + 1
                                  select case(direction)
                                  case('x')
                                    v1(n_tot) =                x(i,j,k,n) ** (p  ) * c_one
                                    v2(n_tot) =  (p-0)       * x(i,j,k,n) ** (p-1) * c_one
                                    v3(n_tot) =  (p-0)*(p-1) * x(i,j,k,n) ** (p-2) * c_one
                                  case('y')
                                    v1(n_tot) =                y(i,j,k,n) ** (p  ) * c_one
                                    v2(n_tot) =  (p-0)       * y(i,j,k,n) ** (p-1) * c_one
                                    v3(n_tot) =  (p-0)*(p-1) * y(i,j,k,n) ** (p-2) * c_one
                                  case('z')
                                    v1(n_tot) =                z(i,j,k,n) ** (p  ) * c_one
                                    v2(n_tot) =  (p-0)       * z(i,j,k,n) ** (p-1) * c_one
                                    v3(n_tot) =  (p-0)*(p-1) * z(i,j,k,n) ** (p-2) * c_one
                                  end select 
                                  if(p==1) v3(n_tot) =  0.0_dp * c_one
                                enddo
                              enddo
                            enddo
                          enddo
                          select case(direction)
                          case('x')
                            write(*,*)' '
                            call diff_xi  (Dxi  ,uu, tmp_b_all, cwrk, bc_type)      
                            call L2_norm(cwrk- du,norm) ; write(*,*)'D1 SBP-Original error  ', direction,' ',norm
                            call a_mat_vec (n_tot,v1,r1, CSR%a1,CSR%ja1,CSR%ia1)
                            call L2_norm(r1-v2,norm) ; write(*,*)'D1 CSR error norm      ', direction,' ',norm

                            if(viscous_flow == 1)then
                              call diff_xi  (D2xi ,uu, tmp_b_all, cwrk, bc_type)      
                              call L2_norm(cwrk-ddu,norm) ; write(*,*)'D2 SBP-Original error  ', direction,' ',norm
                              call a_mat_vec (n_tot,v1,r1, CSR%a4,CSR%ja4,CSR%ia4)
                              call L2_norm(r1-v3,norm) ; write(*,*)'D2 CSR error norm      ', direction,' ',norm
                            endif

                          case('y')
                            write(*,*)' '
                            call diff_eta (Deta  ,uu, tmp_b_all, cwrk, bc_type)      
                            call L2_norm(cwrk- du,norm) ; write(*,*)'D1 SBP-Original error  ', direction,' ',norm
                            call a_mat_vec (n_tot,v1,r1, CSR%a2,CSR%ja2,CSR%ia2)
                            call L2_norm(r1-v2,norm) ; write(*,*)'D1 CSR error norm      ', direction,' ',norm

                            if(viscous_flow == 1)then
                              call diff_eta (D2eta ,uu, tmp_b_all, cwrk, bc_type)      
                              call L2_norm(cwrk-ddu,norm) ; write(*,*)'D2 SBP-Original error  ', direction,' ',norm
                              call a_mat_vec (n_tot,v1,r1, CSR%a5,CSR%ja5,CSR%ia5)
                              call L2_norm(r1-v3,norm) ; write(*,*)'D2 CSR error norm      ', direction,' ',norm
                            endif

                          case('z')
                            if(nz(1) > 1) then
                              write(*,*)' '
                              call diff_zeta(Dzeta ,uu, tmp_b_all, cwrk, bc_type)      
                              call L2_norm(cwrk- du,norm) ; write(*,*)'D1 SBP-Original error  ', direction,' ',norm
                              call a_mat_vec (n_tot,v1,r1, CSR%a3,CSR%ja3,CSR%ia3)
                              call L2_norm(r1-v2,norm) ; write(*,*)'D1 CSR error norm      ', direction,' ',norm

                              if(viscous_flow == 1)then
                                call diff_zeta(D2zeta,uu, tmp_b_all, cwrk, bc_type)      
                                call L2_norm(cwrk-ddu,norm) ; write(*,*)'D2 SBP-Original error  ', direction,' ',norm
                                call a_mat_vec (n_tot,v1,r1, CSR%a6,CSR%ja6,CSR%ia6)
                                call L2_norm(r1-v3,norm) ; write(*,*)'D2 CSR error norm      ', direction,' ',norm
                              endif

                            endif
                          end select
                        enddo
                      endif

                      !     Testing residual

                      cnt = 1
                      do n=1,nblock
                        do k=1,nz(n)
                          do j=1,ny(n)
                            do i=1,nx(n)
                              do L=1,nvar    !   note that variables at each gridpoint are number first.
                                w1(cnt) = u(i,j,k,L,n)
                                w2(cnt) = c_zero
                                cnt  = cnt + 1
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo

                      call a_mat_vec (cnt, w1, w2, CSR%a,CSR%ja,CSR%ia)

                      cnt = 1
                      do n=1,nblock
                        do k=1,nz(n)
                          do j=1,ny(n)
                            do i=1,nx(n)
                              do L=1,nvar
                                if(abs(res(i,j,k,L,n)) >= 1.0e-4) &
                                  write(180,'(4(i3,1x),6(f12.5,1x))')i,j,k,L,res(i,j,k,L,n),w2(cnt)
                                !               write(180,'(4(i3,1x),6(f12.5,1x))')i,j,k,L,res(i,j,k,L,n),w2(cnt)
                                w1(cnt) = res(i,j,k,L,n)-w2(cnt)
                                !               if(abs(res(i,j,k,L,n)-w2(cnt)) >= 1.0e-5) &
                                !               write(*,'(4(i3,1x),6(f12.5,1x))')i,j,k,L,res(i,j,k,L,n),w2(cnt)
                                cnt  = cnt + 1
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo

                      call L2_norm(w1 ,norm) ; write(*,*)'L2_norm inconsistency',norm
                      call L2_norm(res,norm) ; write(*,*)'L2_norm of residual  ',norm

                      do l=1,nblock
                        write (22,200) trim(output_file)
                        write (22,*)'variables = x,y,rho,u,v,w,p'
                        write(22,201) title,nx(l),ny(l)

                        do j=1,ny(l)
                          do i=1,nx(l)

                            write (22,'(7(1x,e16.9))') x(i,j,1,l),y(i,j,1,l),Log10(abs(res(i,j,1,1,l)+eps)),&
                              Log10(abs(res(i,j,1,2,l)+eps)),Log10(abs(res(i,j,1,3,l)+eps)),               &
                              Log10(abs(res(i,j,1,4,l)+eps)),Log10(abs(res(i,j,1,5,l)+eps))
                            !               write (22,'(7(1x,e16.9))') x(i,j,1,l),y(i,j,1,l),abs(u(i,j,1,1,l)),&
                            !                    abs(u(i,j,1,2,l)),abs(u(i,j,1,3,l)),abs(u(i,j,1,4,l)),abs(u(i,j,1,5,l))
                          end do
                        end do
                      end do

                      200 format('title = "',a30,'"',(f20.16))
                      201 format('zone t = "',a30,'", i = ',i5,', j = ',i5,', f = point')



                      return
                    end subroutine testing

                  end module matrix_mod

