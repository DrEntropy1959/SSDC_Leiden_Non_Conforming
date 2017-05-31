! This modules contains all the subroutines to construct the CSR formats needed
! for building the analytical Jacobian matrix.

module initialize_CSR

  ! Load modules
  use precision_vars
  
  ! Nothing is implicitly defined 
  implicit none

  ! Subroutines and functions in this module are usually private
  private

  ! Exceptions, i.e. public subroutines or functions
  public csr_term_footprints, csr_get_pointers, csr_combine_pointers, &
    & csr_testing, csr_initialize_jacobian, csr_combine_pointers_element, &
    & csr_on_element_operator, csr_on_element_Matrix_Multiply


contains

  !============================================================================

  !============================================================================
  ! csr_term_footprints - Constructs the CSR format of the differentiation 
  ! matrices for one element. This is the standard reference way to cosntruct 
  ! the CSR for one element. This "bricks" will then be used to construct the 
  ! CSR arrays which contains the CSR of all the elements owned by a processor. 

  subroutine csr_term_footprints()

    ! Load modules
    use CSRlocalvariables
    use collocationvariables
    use referencevariables
    use nsereferencevariables
    use unary_mod

    ! Nothing is implicitly defined
    implicit none

    integer, allocatable, dimension(:) :: n_degr, iw, iw1
    integer :: job, nrow, ncol, ncolb, i_err, nnz
    logical :: values
    integer :: cnt, i
    character(120) :: message

    ! Subroutine name for possible error message
    message = 'csr_term_footprints'

    ! Allocate CSR storage for tensor product cells
    ! ---------------------------------------------
    ! Number of nnz terms in the first derivatives
    nnz = nodesperelem*nodesperedge 

    ! x1 direction
    allocate(ia_x1_elem(nodesperelem+1))
    ia_x1_elem  = 0

    allocate(ja_x1_elem(nnz))
    ja_x1_elem = 0

    allocate(a_x1_elem(nnz)) 
    a_x1_elem = 0.0_wp
    
    ! x2 direction
    allocate(ia_x2_elem(nodesperelem+1))
    ia_x2_elem  = 0
    
    allocate(ja_x2_elem(nnz))
    ja_x2_elem = 0

    allocate(a_x2_elem(nnz))
    a_x2_elem = 0.0_wp

    ! Set values of the CSR format in x1 and x2 directions
    ia_x1_elem(:) = iagrad(:)
    ia_x2_elem(:) = iagrad(:)

    ja_x1_elem(:) = jagrad(1,:)
    ja_x2_elem(:) = jagrad(2,:)

    a_x1_elem(:) = dagrad(1,:)
    a_x2_elem(:) = dagrad(2,:)

    ! If 3D set the values of the CSR format in the x3 direction
    if(ndim == 3) then
      allocate(ia_x3_elem(nodesperelem+1))
      ia_x3_elem  = 0 

      allocate(ja_x3_elem(nnz))
      ja_x3_elem = 0 
      
      allocate(a_x3_elem(nnz))
      a_x3_elem = 0.0_wp 

      ia_x3_elem(:) = iagrad(:)
      ja_x3_elem(:) = jagrad(3,:)
      a_x3_elem(:) = dagrad(3,:)
    endif

    ! If viscous set the values of the CSR format
    if(viscous) then
      allocate(ia_x11_elem(nodesperelem+1))
      ia_x11_elem  = 0 
      
      allocate(ia_x22_elem(nodesperelem+1))
      ia_x22_elem  = 0 
      
      ! The footprint of DiDj is  [1 <= i,j <= 3]
      ! Viscous second derivative CSR
      job = 1
      values = .true.
      nrow = nodesperelem
      ncol = nodesperelem
      ncolb = nodesperelem

      allocate(n_degr(nrow))
      n_degr = 0

      allocate(iw(ncolb))
      iw = 0

      ! amubdg() gets number of nonzeros in each row
      ! CSR stride and column pointers for second-derivative terms
      ! DxDx
      call amubdg(nrow,ncol,ncolb,ja_x1_elem,ia_x1_elem,ja_x1_elem,ia_x1_elem, &
        & n_degr,nnz,iw) 
      
      allocate(ja_x11_elem(nnz)) 
      ja_x11_elem = 0
      
      allocate(a_x11_elem(nnz))
      a_x11_elem = 0.0_wp

      iw = 0
      call amub(nrow,ncol,job,a_x1_elem,ja_x1_elem,ia_x1_elem,a_x1_elem, &
        & ja_x1_elem,ia_x1_elem,a_x11_elem,ja_x11_elem,ia_x11_elem,nnz,iw,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)
      
      allocate(iw1(2*(ia_x11_elem(nrow+1)-ia_x11_elem(1))))
      iw1 = 0
      call csort(nrow,a_x11_elem,ja_x11_elem,ia_x11_elem,iw1,values)

      ! DyDy
      n_degr = 0
      iw = 0
      call amubdg(nrow,ncol,ncolb,ja_x2_elem,ia_x2_elem,ja_x2_elem,ia_x2_elem, &
        & n_degr,nnz,iw)

      allocate(ja_x22_elem(nnz))
      ja_x22_elem = 0
      
      allocate(a_x22_elem(nnz))
      a_x22_elem = 0.0_wp
      
      iw = 0
      call amub(nrow,ncol,job,a_x2_elem,ja_x2_elem,ia_x2_elem,a_x2_elem, &
        & ja_x2_elem,ia_x2_elem,a_x22_elem,ja_x22_elem,ia_x22_elem,nnz,iw,i_err)
      
      ! Check for error
      call check_sparsekit_error(i_err,message)
      
      deallocate(iw1) 
      allocate(iw1(2*(ia_x22_elem(nrow+1)-ia_x22_elem(1))))
      iw1 = 0
      call csort(nrow,a_x22_elem,ja_x22_elem,ia_x22_elem,iw1,values)

      if(ndim == 3) then
        allocate(ia_x33_elem(nodesperelem+1))
        ia_x33_elem = 0

        ! DzDz
        n_degr = 0
        iw = 0
        call amubdg(nrow,ncol,ncolb,ja_x3_elem,ia_x3_elem,ja_x3_elem, &
          & ia_x3_elem,n_degr,nnz,iw)
        
        allocate(ja_x33_elem(nnz))
        ja_x33_elem = 0
        
        allocate(a_x33_elem(nnz))
        a_x33_elem = 0.0_wp
        
        iw = 0
        call amub(nrow,ncol,job,a_x3_elem,ja_x3_elem,ia_x3_elem,a_x3_elem, &
          & ja_x3_elem,ia_x3_elem,a_x33_elem,ja_x33_elem,ia_x33_elem,nnz,iw, &
          & i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)
        
        deallocate(iw1)
        allocate(iw1(2*(ia_x33_elem(nrow+1)-ia_x33_elem(1))))
        iw1 = 0
        call csort(nrow,a_x33_elem,ja_x33_elem,ia_x33_elem,iw1,values)
      endif

      ! If crossterms set the values of the CSR format
      if(crossterms) then

        allocate(ia_x12_elem(nodesperelem+1))
        ia_x12_elem = 0 

        if(ndim == 3) then
          allocate(ia_x13_elem(nodesperelem+1))
          ia_x13_elem = 0

          allocate(ia_x23_elem(nodesperelem+1))
          ia_x23_elem = 0 
        endif


        ! CSR stride and column pointers for Cross-terms
        ! DxDy
        n_degr = 0
        iw = 0
        call amubdg(nrow,ncol,ncolb,ja_x1_elem,ia_x1_elem,ja_x2_elem, &
          & ia_x2_elem,n_degr,nnz,iw)
        
        allocate(ja_x12_elem(nnz))
        ja_x12_elem = 0
        
        allocate(a_x12_elem(nnz))
        a_x12_elem = 0.0_wp

        iw = 0
        call amub(nrow,ncol,job,a_x1_elem,ja_x1_elem,ia_x1_elem,a_x2_elem, &
          & ja_x2_elem,ia_x2_elem,a_x12_elem,ja_x12_elem,ia_x12_elem,nnz,iw, &
          & i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        deallocate(iw)
        allocate(iw(max(nrow+1,2*nnz)))
        iw = 0
        call csort(nrow,a_x12_elem,ja_x12_elem,ia_x12_elem,iw,values)

        if(ndim == 3) then
          ! DxDz
          n_degr = 0
          iw = 0
          call amubdg(nrow,ncol,ncolb,ja_x1_elem,ia_x1_elem,ja_x3_elem, &
            & ia_x3_elem,n_degr,nnz,iw) 
          
          allocate(ja_x13_elem(nnz))
          ja_x13_elem = 0
          
          allocate(a_x13_elem(nnz))
          a_x13_elem = 0.0_wp
          
          deallocate(iw)
          allocate(iw(max(nrow+1,2*nnz))) 
          iw = 0
          call amub(nrow,ncol,job,a_x1_elem,ja_x1_elem,ia_x1_elem,a_x3_elem, &
            & ja_x3_elem,ia_x3_elem,a_x13_elem,ja_x13_elem,ia_x13_elem,nnz, &
            & iw,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
         
          deallocate(iw)
          allocate(iw(max(nrow+1,2*nnz)))
          iw = 0
          call csort(nrow,a_x13_elem,ja_x13_elem,ia_x13_elem,iw,values)
          
          ! DyDz
          n_degr = 0
          iw = 0
          call amubdg(nrow,ncol,ncolb,ja_x2_elem,ia_x2_elem,ja_x3_elem, &
            & ia_x3_elem,n_degr,nnz,iw)
          allocate(ja_x23_elem(nnz))
          ja_x23_elem = 0
          
          allocate(a_x23_elem(nnz))
          a_x23_elem = 0.0_wp

          deallocate(iw) ; allocate(iw(max(nrow+1,2*nnz))) 
          iw = 0
          call amub(nrow,ncol,job,a_x2_elem,ja_x2_elem,ia_x2_elem,a_x3_elem, &
            & ja_x3_elem,ia_x3_elem,a_x23_elem,ja_x23_elem,ia_x23_elem,nnz, &
            & iw,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
          
          deallocate(iw) ; allocate(iw(max(nrow+1,2*nnz)))
          iw = 0
          call csort(nrow,a_x23_elem,ja_x23_elem,ia_x23_elem,iw,values)
        endif

      endif
    endif

    ! Terms arising from I/dt
    ! I = identity matrix
    allocate(ia_0_elem(nodesperelem+1))
    allocate(ja_0_elem(nodesperelem))
    allocate(a_0_elem(nodesperelem))
    ia_0_elem = 0
    ja_0_elem = 0
    a_0_elem = 0.0_wp

    cnt = 1
    ia_0_elem(1) = 1

    do i=1,nodesperelem
      ia_0_elem(cnt+1) = ia_0_elem(cnt) + 1
      ja_0_elem(cnt) = cnt
      a_0_elem(cnt) = one
      cnt = cnt + 1 
    enddo

    return
  end subroutine csr_term_footprints

  !============================================================================

  !============================================================================
  ! csr_get_pointers - Constructs the CSR format of the differentiation 
  ! matrices for all the elements owned by a processor using the reference 
  ! standard "bricks" constructed in csr_term_footprints().

  subroutine csr_get_pointers()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use nsereferencevariables
    use unary_mod
    use controlvariables 

    ! Nothing is implicitly defined
    implicit none

    integer :: i, j, l
    integer :: cnt, nnz
    integer :: i_elem
    integer :: n_elems
    integer :: low_elem, high_elem
    integer :: global_shift

    continue

    ! Low volumetric element index
    low_elem = ihelems(1)

    ! High volumetric element index
    high_elem = ihelems(2)

    n_elems = 1 + high_elem - low_elem

    ! Shift for global indexing for parallel computation
    global_shift = (low_elem-1)*nodesperelem

    ! Diagonal terms stored in ia_0
    cnt = 0
    nnz = 0
    do l = 1, n_elems
      cnt = cnt + nodesperelem
      nnz = nnz + nodesperelem
    enddo

    allocate(ia_0(cnt+1))    ;   ia_0 = 0
    allocate(ja_0(nnz))      ;   ja_0 = 0
    allocate(ka_0(nnz))      ;   ka_0 = 0
    allocate( a_0(nnz))      ;    a_0 = 0.0_wp

    cnt = 1
    ia_0(cnt) = 1
    do i_elem = 1, n_elems
      do i = 1, nodesperelem
        ia_0(cnt+1) = ia_0(cnt) + 1
        ja_0(cnt) = cnt
        a_0(cnt) = one
        cnt = cnt + 1 
      enddo
    enddo

    ! Shift to global indexing for parallel computations 
    ja_0 = ja_0 + global_shift

    ! Terms in the xi direction
    if(ndim >= 1) then

      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x1(cnt+1))    ;   ia_x1 = 0
      allocate(ja_x1(nnz))      ;   ja_x1 = 0
      allocate(ka_x1(nnz))      ;   ka_x1 = 0
      allocate( a_x1(nnz))      ;    a_x1 = 0.0_wp

      nnz = 1
      cnt = 1
      ia_x1(1) = cnt
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x1(cnt) = ia_x1_elem(i) + (i_elem-1)*nodesperelem*nodesperedge 
          do j = ia_x1_elem(i), ia_x1_elem(i+1) - 1
            ja_x1(nnz) = ja_x1_elem(j) + (i_elem-1)*nodesperelem
            a_x1(nnz) = a_x1_elem(j)

            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x1 (cnt) = ia_x1_elem(nodesperelem+1) + (i_elem-1)*nodesperelem* &
          & nodesperedge
      end do

      ! Shift to global indexing for parallel computations 
      ja_x1 = ja_x1 + global_shift

    endif


    ! Terms in the eta direction
    if(ndim >= 2) then

      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x2(cnt+1))    ;   ia_x2 = 0
      allocate(ja_x2(nnz))      ;   ja_x2 = 0
      allocate(ka_x2(nnz))      ;   ka_x2 = 0
      allocate( a_x2(nnz))      ;    a_x2 = 0.0_wp

      nnz = 1 
      cnt = 1
      ia_x2(cnt) = cnt
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x2(cnt) = ia_x2_elem(i) + (i_elem-1)*nodesperelem*nodesperedge
          do j = ia_x2_elem(i), ia_x2_elem(i+1)-1
            ja_x2(nnz) = ja_x2_elem(j) + (i_elem-1)*nodesperelem
            a_x2(nnz) = a_x2_elem(j)

            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x2(cnt) = ia_x2_elem(nodesperelem+1) + (i_elem-1)*nodesperelem* &
          & nodesperedge
      end do

      ! Shift to global indexing for parallel computations 
      ja_x2 = ja_x2 + global_shift

    endif


    ! Terms in the zeta direction
    if(ndim >= 3) then

      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x3(cnt+1))    ;  ia_x3 = 0
      allocate(ja_x3(nnz))      ;  ja_x3 = 0
      allocate(ka_x3(nnz))      ;  ka_x3 = 0
      allocate( a_x3(nnz))      ;   a_x3 = 0.0_wp

      nnz = 1
      cnt = 1
      ia_x3(cnt) = cnt
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x3(cnt) = ia_x3_elem(i) + (i_elem-1)*nodesperelem*nodesperedge
          do j = ia_x3_elem(i), ia_x3_elem(i+1)-1
            ja_x3(nnz) = ja_x3_elem(j) + (i_elem-1)*nodesperelem
            a_x3(nnz) = a_x3_elem(j)

            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x3(cnt) = ia_x3_elem(nodesperelem+1) + (i_elem-1)*nodesperelem* &
          & nodesperedge
      end do

      ! Shift to global indexing for parallel computations 
      ja_x3 = ja_x3 + global_shift

    endif

    ! Second derivative path
    if(viscous) then
      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x11(cnt+1))   ;   ia_x11 = 0
      allocate(ja_x11(nnz))     ;   ja_x11 = 0
      allocate(ka_x11(nnz))     ;   ka_x11 = 0
      allocate( a_x11(nnz))     ;    a_x11 = 0.0_wp
      
      allocate(ia_x22(cnt+1))   ;   ia_x22 = 0
      allocate(ja_x22(nnz))     ;   ja_x22 = 0
      allocate(ka_x22(nnz))     ;   ka_x22 = 0
      allocate( a_x22(nnz))     ;    a_x22 = 0.0_wp

      if(ndim == 3) then
        allocate(ia_x33(cnt+1))  ;  ia_x33 = 0
        allocate(ja_x33(nnz))    ;  ja_x33 = 0
        allocate(ka_x33(nnz))    ;  ka_x33 = 0
        allocate( a_x33(nnz))    ;   a_x33 = 0.0_wp
      endif

      nnz = 1 
      cnt = 1 
      ia_x11(1) = cnt 
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x11(cnt) = ia_x11_elem(i) + (i_elem-1)*nodesperelem*nodesperedge
          do j = ia_x11_elem(i), ia_x11_elem(i+1)-1
            ja_x11(nnz) = ja_x11_elem(j) + (i_elem-1)*nodesperelem
            a_x11(nnz) = a_x11_elem(j)
            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x11(cnt) = ia_x11_elem(nodesperelem+1) + (i_elem-1)*nodesperelem* &
          & nodesperedge
      end do
      
      ! Shift to global indexing for parallel computations 
      ja_x11 = ja_x11 + global_shift

      
      nnz = 1 
      cnt = 1 
      ia_x22(1) = cnt 
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x22(cnt) = ia_x22_elem(i) + (i_elem-1)*nodesperelem*nodesperedge
          do j = ia_x22_elem(i), ia_x22_elem(i+1)-1
            ja_x22(nnz) = ja_x22_elem(j) + (i_elem-1)*nodesperelem
            a_x22(nnz) = a_x22_elem(j)
            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x22(cnt) = ia_x22_elem(nodesperelem+1) + (i_elem-1)*nodesperelem* &
          & nodesperedge
      end do

      ! Shift to global indexing for parallel computations 
      ja_x22 = ja_x22 + global_shift


      if(ndim == 3) then
        nnz = 1 
        cnt = 1 
        ia_x33(1) = cnt
        do i_elem = 1, n_elems
          do i = 1, nodesperelem
            ia_x33(cnt) = ia_x33_elem(i) + (i_elem-1)*nodesperelem*nodesperedge
            do j = ia_x33_elem(i),ia_x33_elem(i+1)-1
              ja_x33(nnz) = ja_x33_elem(j) + (i_elem-1)*nodesperelem
              a_x33(nnz) = a_x33_elem(j)
              nnz = nnz + 1
            enddo 
            cnt = cnt + 1
          enddo 
          ia_x33(cnt) = ia_x33_elem(nodesperelem+1) + (i_elem-1)*nodesperelem* &
            & nodesperedge
        end do

        ! Shift to global indexing for parallel computations 
        ja_x33 = ja_x33 + global_shift

      endif


      if(crossterms) then
        if(ndim >= 2) then
          cnt = 0
          nnz = 0
          do l = 1, n_elems
            cnt = cnt + nodesperelem
            nnz = nnz + nodesperelem*nodesperedge*nodesperedge
          enddo

          allocate(ia_x12(cnt+1))    ; ia_x12 = 0
          allocate(ja_x12(nnz))      ; ja_x12 = 0
          allocate(ka_x12(nnz))      ; ka_x12 = 0
          allocate( a_x12(nnz))      ;  a_x12 = 0.0_wp

          nnz = 1
          cnt = 1
          ia_x12(1) = cnt
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x12(cnt) = ia_x12_elem(i) + (i_elem-1)*nodesperelem* &
                & nodesperedge*nodesperedge
              do j = ia_x12_elem(i), ia_x12_elem(i+1)-1
                ja_x12(nnz) = ja_x12_elem(j) + (i_elem-1)*nodesperelem
                a_x12(nnz) = a_x12_elem(j)
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x12(cnt) = ia_x12_elem(nodesperelem+1) + (i_elem-1)* &
              & nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x12 = ja_x12 + global_shift

        endif


        if(ndim == 3) then

          cnt = 0
          nnz = 0
          do l = 1, n_elems
            cnt = cnt + nodesperelem
            nnz = nnz + nodesperelem*nodesperedge*nodesperedge
          enddo

          allocate(ia_x13(cnt+1))     ;   ia_x13 = 0
          allocate(ja_x13(nnz))       ;   ja_x13 = 0
          allocate(ka_x13(nnz))       ;   ka_x13 = 0
          allocate( a_x13(nnz))       ;    a_x13 = 0.0_wp
          
          allocate(ia_x23(cnt+1))     ;   ia_x23 = 0
          allocate(ja_x23(nnz))       ;   ja_x23 = 0
          allocate(ka_x23(nnz))       ;   ka_x23 = 0
          allocate( a_x23(nnz))       ;    a_x23 = 0.0_wp

          nnz = 1 
          cnt = 1
          ia_x13(1) = 1 
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x13(cnt) = ia_x13_elem(i) + (i_elem-1)*nodesperelem* &
                & nodesperedge*nodesperedge
              do j = ia_x13_elem(i), ia_x13_elem(i+1)-1
                ja_x13(nnz) = ja_x13_elem(j) + (i_elem-1)*nodesperelem
                 a_x13(nnz) =  a_x13_elem(j)
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x13(cnt) = ia_x13_elem(nodesperelem+1) + (i_elem-1)* &
              & nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x13 = ja_x13 + global_shift


          nnz = 1
          cnt = 1 
          ia_x23(1) = 1
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x23(cnt) = ia_x23_elem(i) + (i_elem-1)*nodesperelem* &
                & nodesperedge*nodesperedge
              do j = ia_x23_elem(i), ia_x23_elem(i+1)-1
                ja_x23(nnz) = ja_x23_elem(j) + (i_elem-1)*nodesperelem
                 a_x23(nnz) =  a_x23_elem(j)
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x23(cnt) = ia_x23_elem(nodesperelem+1) + (i_elem-1)* &
              & nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x23 = ja_x23 + global_shift

        endif

      endif

    endif

    ! Penalties:  Inviscid and Viscous
    if(IMEX_penalty == 'implicit') then

       ! penI = cI * (F_L^{I} - F_R^{I} )     :  Inviscid Penalty
       ! penV = cV * (F_L^{v} - F_R^{v} )     :  Viscous  Penalty

       call CSR_SAT_Penalty()

    endif

    return
  end subroutine csr_get_pointers

  !============================================================================
  
  !============================================================================
  subroutine CSR_SAT_Penalty()

    use CSRlocalvariables
    use nsereferencevariables
    use variables
    use referencevariables
    use controlvariables
    use unary_mod
!   use collocationvariables, only: iagrad,jagrad,dagrad,gradmat,pinv,pvol

    implicit none
    integer, parameter                    :: bigN = 200  ! Max number of NNZ per row for each operator

    !integer , allocatable, dimension(:)   :: iw, indu

    integer , allocatable, dimension(:)   :: ia_penI_T  , ia_penV_T  
    integer , allocatable, dimension(:)   :: ja_penI_T  , ja_penV_T  
    integer , allocatable, dimension(:)   :: la_penI_T  , la_penV_T  

    real(wp), allocatable, dimension(:)   ::  a_penI_T  ,  a_penV_T

    integer , allocatable, dimension(:)   :: ia_penI_tmp, ia_penV_tmp
    integer , allocatable, dimension(:,:) :: ja_penI_tmp, ja_penV_tmp
    integer , allocatable, dimension(:,:) :: la_penI_tmp, la_penV_tmp

    ! indices
    integer :: i_elem, n_elems
    integer :: inode, jnode, knode, kelem, iface
    integer :: low_elem, high_elem, mat_dim_proc, local_row
    integer :: i, j
    integer :: cntI, cntV
    integer :: cnt_penI_nnz, cnt_penV_nnz

    real(wp),              dimension(3)   :: n_v

    ! low : high volumetric element index
    low_elem = ihelems(1) ;  high_elem = ihelems(2)

    n_elems = 1 + high_elem - low_elem

    mat_dim_proc = n_elems * nodesperelem

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Compute boundary conditions pointers for CSR
    ! Inviscid and Viscous use same strategy
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    !  Temporary storage used to count terms

    allocate(ia_penI_tmp(mat_dim_proc))        ;  ia_penI_tmp(:)   = 0
    allocate(ja_penI_tmp(mat_dim_proc,bigN))   ;  ja_penI_tmp(:,:) = 0
    allocate(la_penI_tmp(mat_dim_proc,bigN))   ;  la_penI_tmp(:,:) = 0

    allocate(ia_penV_tmp(mat_dim_proc))        ;  ia_penV_tmp(:)   = 0
    allocate(ja_penV_tmp(mat_dim_proc,bigN))   ;  ja_penV_tmp(:,:) = 0
    allocate(la_penV_tmp(mat_dim_proc,bigN))   ;  la_penV_tmp(:,:) = 0

    ! First  pass over elements (elemloop) counts the number and position of all SAT entries 
    !    Stores data in two-dimensional arrays: ia_penI_tmp, ja_penI_tmp, la_penI_tmp

    ! Second pass, assembles data in CSR-format in  ia_penI_T, ja_penI_T

    ! Third  pass stores data in local arrays:  ia_penI_proc, ja_penI_proc, ka_penI_proc of appropriate size
    !    ia_penI_proc still contains redundent entries and is consistent with how SAT terms are visited in assembly

    ! Fourth pass: Clean up the CSR to remove duplicate column entries and reorder in assending order
    !  Data is squeezed down "in place" in the ia_penI_T . . . a_penI_T arrays.
    !  Compressed data is stored in ia_penI and ia_penV  arays

    cntI = 0 ; cntV = 0 ;

    ! Begin First pass : Loop over all elements
    elemloop: do i_elem = low_elem, high_elem

      ! loop over each face
      faceloop:do iface = 1,nfacesperelem
        ! check to see if face is a boundary face
        if (ef2e(1,iface,i_elem) < 0) then
          ! if we get here, we are on a boundary face
          ! specify the Boundary Condition procedure on the face
 
          ! loop over each node on the face
          do i = 1,nodesperface
            ! volumetric node index corresponding to face and node on face indices
            inode = kfacenodes(i,iface)
            
            ! facial index corresponding to face and node on face indices
            jnode = (iface-1)*nodesperface + i
            
            ! facial index corresponding to face and node on face indices
            knode = inode
            ! element index of partner node (which is itself because it's a boundary)
            kelem = i_elem

            ! Outward facing normal of facial node
            n_v = facenodenormal(:,jnode,i_elem)

            local_row = (i_elem-low_elem)*nodesperelem + inode

            
            cntI = cntI + 1
            ia_penI_tmp(local_row) = ia_penI_tmp(local_row) + 1
            ja_penI_tmp(local_row,ia_penI_tmp(local_row)) = (i_elem-1)*nodesperelem + inode
            la_penI_tmp(local_row,ia_penI_tmp(local_row)) = cntI

            cntI = cntI + 1
            ia_penI_tmp(local_row) = ia_penI_tmp(local_row) + 1
            ja_penI_tmp(local_row,ia_penI_tmp(local_row)) = (kelem-1)*nodesperelem + knode
            la_penI_tmp(local_row,ia_penI_tmp(local_row)) = cntI

            if(viscous) then
              !  Build ``On and Off  element'' contributions
              call Viscous_Pen_Initialize_CSR(iface,local_row,inode,i_elem,knode,kelem,cntV, &
                                              ia_penV_tmp, ja_penV_tmp, la_penV_tmp ) 
            endif

          end do

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       Off Processor Contributions to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        else if (ef2e(3,iface,i_elem) /= myprocid) then
          ! this is a parallel interface
          do i = 1, nodesperface
            jnode = nodesperface*(iface-1)+i
            ! volumetric node index corresponding to facial node index
            inode = ifacenodes(jnode)
            ! index in ghost
            knode = efn2efn(1,jnode,i_elem)
            ! element index of partner node
            kelem = efn2efn(2,jnode,i_elem)
            ! Outward facing normal of facial node
            n_v = facenodenormal(:,jnode,i_elem)

            local_row = (i_elem-low_elem)*nodesperelem + inode

!           Both On and Off process info is accumulated into process CSR

            cntI = cntI + 1
            ia_penI_tmp(local_row) = ia_penI_tmp(local_row) + 1
            ja_penI_tmp(local_row,ia_penI_tmp(local_row)) = (i_elem-1)*nodesperelem + inode
            la_penI_tmp(local_row,ia_penI_tmp(local_row)) = cntI

            cntI = cntI + 1
            ia_penI_tmp(local_row) = ia_penI_tmp(local_row) + 1
            ja_penI_tmp(local_row,ia_penI_tmp(local_row)) = ( kelem-1)*nodesperelem + knode
            la_penI_tmp(local_row,ia_penI_tmp(local_row)) = cntI

            if(viscous) then
              !  Build ``On and Off  element'' contributions
              call Viscous_Pen_Initialize_CSR(iface,local_row,inode,i_elem,knode,kelem,cntV, &
                                              ia_penV_tmp, ja_penV_tmp, la_penV_tmp ) 
            endif
          end do

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       ON-Processor Contributions to gsat
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        else
          do i = 1,nodesperface
            ! index in facial ordering
            jnode = nodesperface*(iface-1)+i
            ! volumetric node index corresponding to facial node index
            inode = ifacenodes(jnode)
            ! volumetric index of partner node
            knode = efn2efn(1,jnode,i_elem)
            ! element index of partner node
            kelem = efn2efn(2,jnode,i_elem)
            ! Outward facing normal of facial node
            n_v = facenodenormal(:,jnode,i_elem)

            local_row = (i_elem-low_elem)*nodesperelem + inode

            cntI = cntI + 1
            ia_penI_tmp(local_row) = ia_penI_tmp(local_row) + 1
            ja_penI_tmp(local_row,ia_penI_tmp(local_row)) = (i_elem-1)*nodesperelem + inode
            la_penI_tmp(local_row,ia_penI_tmp(local_row)) = cntI

            cntI = cntI + 1
            ia_penI_tmp(local_row) = ia_penI_tmp(local_row) + 1
            ja_penI_tmp(local_row,ia_penI_tmp(local_row)) = ( kelem-1)*nodesperelem + knode
            la_penI_tmp(local_row,ia_penI_tmp(local_row)) = cntI

            if(viscous) then
              !  Build ``On and Off  element'' contributions
              call Viscous_Pen_Initialize_CSR(iface,local_row,inode,i_elem,knode,kelem,cntV, &
                                              ia_penV_tmp, ja_penV_tmp, la_penV_tmp ) 
            endif
          end do
        end if

      enddo faceloop

    enddo elemloop

!   Flatten the 2-Dimensional Inviscid CSR into conventional CSR format

    ! Begin Second pass :    Count nnz and allocate ja_penI_T
    allocate(ia_penI_T(mat_dim_proc+1))      ;  ia_penI_T(:) = 0
    ia_penI_T(1) = 1
    do i = 1, mat_dim_proc
      ia_penI_T(i+1) = ia_penI_T(i) + ia_penI_tmp(i)
    enddo
    allocate(ja_penI_T(ia_penI_T(mat_dim_proc+1)-1))   ;  ja_penI_T(:) = 0
    allocate(la_penI_T(ia_penI_T(mat_dim_proc+1)-1))   ;  la_penI_T(:) = 0
    allocate( a_penI_T(ia_penI_T(mat_dim_proc+1)-1))   ;   a_penI_T(:) = zero

    !   form ja_penI_T in CSR format
    cnt_penI_nnz = 0
    do i = 1, mat_dim_proc
      if(ia_penI_tmp(i) >= 1) then
        do j = 1,ia_penI_tmp(i)
          cnt_penI_nnz = cnt_penI_nnz + 1
          ja_penI_T(cnt_penI_nnz) = ja_penI_tmp(i,j)
          la_penI_T(cnt_penI_nnz) = la_penI_tmp(i,j)
        enddo
      endif
    enddo
    
    !write(*,*) 'ja_penI_tmp, my procid', myprocid, ja_penI_tmp
    deallocate(ia_penI_tmp, ja_penI_tmp, la_penI_tmp )

    ! Begin Third pass
    allocate(ia_penI_proc(mat_dim_proc+1))                 ;  ia_penI_proc(:) = ia_penI_T(:)
    allocate(ja_penI_proc(ia_penI_T(mat_dim_proc+1)-1))    ;  ja_penI_proc(:) = ja_penI_T(:)
    allocate(ka_penI_proc(ia_penI_T(mat_dim_proc+1)-1))    ;  ka_penI_proc(:) = 0
    allocate(la_penI_proc(ia_penI_T(mat_dim_proc+1)-1))    ;  la_penI_proc(:) = la_penI_T(:)

   
    !  Begin Fourth pass:  clean CSR data and copy into final CSR arrays

 !  allocate(iw(mat_dim_proc+1),indu(mat_dim_proc))
 !  call clncsr(3,0,mat_dim_proc,a_penI_T,ja_penI_T,ia_penI_T,indu,iw)
 !  deallocate(iw,indu)

    call clean_CSR(mat_dim_proc, ja_penI_T, ia_penI_T)


    allocate(ia_penI(mat_dim_proc+1))
    allocate(ja_penI(ia_penI_T(mat_dim_proc+1)-1))
    allocate( a_penI(ia_penI_T(mat_dim_proc+1)-1))

    do i = 1,mat_dim_proc+1
      ia_penI(i) = ia_penI_T(i)
    enddo
    do i = 1,ia_penI(mat_dim_proc+1)-1
      ja_penI(i) = ja_penI_T(i)
    enddo
    deallocate(ia_penI_T, ja_penI_T, a_penI_T)

    !write(*,*) 'proc, ja_pen_proc', myprocid, ja_penI_proc
    !stop

    if(viscous) then

      ! Begin Second pass :    Count nnz and allocate ja_penV_T
      allocate(ia_penV_T(mat_dim_proc+1))      ;  ia_penV_T(:) = 0
      ia_penV_T(1) = 1
      do i = 1, mat_dim_proc
        ia_penV_T(i+1) = ia_penV_T(i) + ia_penV_tmp(i)
      enddo
      allocate(ja_penV_T(ia_penV_T(mat_dim_proc+1)-1))   ;  ja_penV_T(:) = 0
      allocate(la_penV_T(ia_penV_T(mat_dim_proc+1)-1))   ;  la_penV_T(:) = 0
      allocate( a_penV_T(ia_penV_T(mat_dim_proc+1)-1))   ;   a_penV_T(:) = zero

      cnt_penV_nnz = 0
      do i = 1, mat_dim_proc
        if(ia_penV_tmp(i) >= 1) then
          do j = 1,ia_penV_tmp(i)
            cnt_penV_nnz = cnt_penV_nnz + 1
            ja_penV_T(cnt_penV_nnz) = ja_penV_tmp(i,j)
            la_penV_T(cnt_penV_nnz) = la_penV_tmp(i,j)
          enddo
        endif
      enddo
      deallocate(ia_penV_tmp, ja_penV_tmp, la_penV_tmp)

      ! Begin Third pass
      allocate(ia_penV_proc(mat_dim_proc+1))                 ;  ia_penV_proc(:) = ia_penV_T(:)
      allocate(ja_penV_proc(ia_penV_T(mat_dim_proc+1)-1))    ;  ja_penV_proc(:) = ja_penV_T(:)
      allocate(ka_penV_proc(ia_penV_T(mat_dim_proc+1)-1))    ;  ka_penV_proc(:) = 0
      allocate(la_penV_proc(ia_penV_T(mat_dim_proc+1)-1))    ;  la_penV_proc(:) = la_penV_T(:)

      !  Clean up the CSR to remove duplicate column entries and reorder in assending order
!     allocate(iw(mat_dim_proc+1),indu(mat_dim_proc))
!     call clncsr(3,0,mat_dim_proc,a_penV_T,ja_penV_T,ia_penV_T,indu,iw)
!     deallocate(iw,indu)

      call clean_CSR(mat_dim_proc, ja_penV_T, ia_penV_T)

      allocate(ia_penV(mat_dim_proc+1))
      allocate(ja_penV(ia_penV_T(mat_dim_proc+1)-1))
      allocate(ka_penV(ia_penV_T(mat_dim_proc+1)-1))
      allocate( a_penV(ia_penV_T(mat_dim_proc+1)-1))

      do i = 1,mat_dim_proc+1
        ia_penV(i) = ia_penV_T(i)
      enddo
      do i = 1,ia_penV(mat_dim_proc+1)-1
        ja_penV(i) = ja_penV_T(i)
      enddo
      deallocate(ia_penV_T, ja_penV_T, la_penV_T, a_penV_T) 

    endif

    return
  end subroutine CSR_SAT_Penalty

  subroutine Viscous_Pen_Initialize_CSR(iface,local_row,inode,i_elem,knode,kelem,cntV, &
              ia_penV_tmp, ja_penV_tmp, la_penV_tmp )

    use CSRlocalvariables
    use nsereferencevariables
    use variables
    use referencevariables
    use controlvariables
    use collocationvariables, only: iagrad,jagrad,dagrad,dmat

    implicit none
    integer,                  intent(in   ) :: iface, local_row, inode, i_elem, knode, kelem
    integer,                  intent(inout) :: cntV
    integer, dimension(:),    intent(inout) :: ia_penV_tmp
    integer, dimension(:,:),  intent(inout) :: ja_penV_tmp
    integer, dimension(:,:),  intent(inout) :: la_penV_tmp

    integer              :: jdir,i, jnode

      ! On element  (use inode and i_elem)
      do jdir = 1,ndim

        do i = iagrad(inode), iagrad(inode+1)-1

          cntV = cntV + 1
          ia_penV_tmp(local_row) = ia_penV_tmp(local_row) + 1
          ja_penV_tmp(local_row,ia_penV_tmp(local_row)) = jagrad(jdir,i) + (i_elem-1) * nodesperelem
          la_penV_tmp(local_row,ia_penV_tmp(local_row)) = cntV
        end do
      end do

      ! Off element  (use knode and kelem)
      do jdir = 1,ndim

        if(jdir /= abs(facenormalcoordinate(iface))) then
          do i = iagrad(knode), iagrad(knode+1)-1
  
            cntV = cntV + 1
            ia_penV_tmp(local_row) = ia_penV_tmp(local_row) + 1
            ja_penV_tmp(local_row,ia_penV_tmp(local_row)) = jagrad(jdir,i) + (kelem-1) * nodesperelem
            la_penV_tmp(local_row,ia_penV_tmp(local_row)) = cntV
          end do
        else
          do i = iagrad(knode), iagrad(knode+1)-1
!           if(jnode == knode) then
!             jnode = inode + (i_elem-1) * nodesperelem
!           else
              jnode = jagrad(jdir,i) + (kelem-1) * nodesperelem
!           endif
  
            cntV = cntV + 1
            ia_penV_tmp(local_row) = ia_penV_tmp(local_row) + 1
            ja_penV_tmp(local_row,ia_penV_tmp(local_row)) = jnode
            la_penV_tmp(local_row,ia_penV_tmp(local_row)) = cntV
          end do
        endif
      end do

      return

  end subroutine Viscous_Pen_Initialize_CSR
  !============================================================================

  !============================================================================
  ! csr_combine_pointers - Combines the CSR formats of the differentiation 
  ! matrices for all the elements owned by a processor using the reference 
  ! standard "bricks" constructed in csr_get_pointers().
  
  subroutine csr_combine_pointers()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use nsereferencevariables
    use unary_mod
    use controlvariables
    use variables
    use collocationvariables

    ! Nothing is implicitly defined
    implicit none

    integer :: i,j,k
    integer :: i_elem, i_node
    integer :: n_elems

    integer , allocatable , dimension(:) :: n_degr
    integer , allocatable , dimension(:) :: iw
    integer :: nnz, n_rows_0, n_cols_0, iface, ii, jj, local_row, jdir
    integer :: inode, jnode, knode, kelem
    integer , dimension(10) :: status_allocation

    integer , allocatable , dimension(:) :: ia_tmp
    integer , allocatable , dimension(:) :: ja_tmp
    real(wp) , allocatable , dimension(:) :: a_tmp
    integer :: i_err
    real(wp) , allocatable , dimension(:) :: aS
    integer :: low_elem, high_elem
    character(120) :: message
    integer :: shift, ll, cntI, cntV
    integer :: ngrowsS

    continue

    ! Subroutine name for possible error message
    message = 'csr_combine_pointers'

    ! Low volumetric element index
    low_elem = ihelems(1)

    ! High volumetric element index
    high_elem = ihelems(2)

    ! Number of elements on the process
    n_elems = 1 + high_elem - low_elem

    ! Check iaS with another approach, i.e. using the sparskit
    n_rows_0 = size(ia_0)-1
    n_cols_0 = size(ia_0)-1

    ! Initialize i_err to normal value
    i_err = 0

    ! Number of global rows in the S matrix
    ngrowsS = ngrows / nequations
    allocate(iw(ngrowS)) ; iw = 0 

    ! Allocate memory for  n_degr and iw (sparkit variables)   
    status_allocation = 0
    allocate(n_degr(n_rows_0),stat=status_allocation(1)) ; n_degr = 0

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating n_degr and or iw in  &
        & csr_combine_pointers(): stopping' 
      stop
    endif


    ! Count number of nonzero element in (a_0 + a_x1)
    ! -----------------------------------------------
    call aplbdg(n_rows_0,ngrowsS,ja_0,ia_0,ja_x1,ia_x1,n_degr,nnz,iw)


    ! Allocate memory for temporary arrays
    allocate(ia_tmp(size(ia_0)),stat=status_allocation(1)) ; ia_tmp = 0
    allocate(ja_tmp(nnz),       stat=status_allocation(2)) ; ja_tmp = 0
    allocate( a_tmp(nnz),       stat=status_allocation(3)) ;  a_tmp = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating ia_tmp and or ja_tmp and or a_tmp in &
        & csr_combine_pointers(): stopping'
      stop
    endif

    ! Combine a_0 with a_x1
    call aplb1(n_rows_0,ngrowsS,0,a_0,ja_0,ia_0,a_x1,ja_x1,ia_x1,a_tmp, &
      & ja_tmp,ia_tmp,nnz,i_err)

    ! Check for error
    call check_sparsekit_error(i_err,message)

    ! Allocate memory for temporary arrays
    allocate(iaS(size(ia_0)),stat=status_allocation(1)) ; iaS = 0
    allocate(jaS(nnz),       stat=status_allocation(2)) ; jaS = 0
    allocate( aS(nnz),       stat=status_allocation(3)) ;  aS = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating iaS and or jaS and or aS in &
        & csr_combine_pointers(): stopping'
      stop
    endif

    ! Assign temporary values to iaS and jaS
    iaS(1) = ia_tmp(1)
    do i = 1,n_rows_0
      iaS(i+1) = ia_tmp(i+1)
      do j = ia_tmp(i),ia_tmp(i+1)-1
        jaS(j) = ja_tmp(j)
      enddo
    enddo


    ! Count number of nonzero element in (a_0 + a_x1) + a_x2
    ! ------------------------------------------------------
    n_degr = 0
    iw = 0
    call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x2,ia_x2,n_degr,nnz,iw) 

    ! Deallocate and allocate memory for new combination
    ia_tmp = 0
    deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
    deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp    
    
    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
        & csr_combine_pointers(): stopping'
      stop
    endif

    ! Combine (a_0 + a_x1) with a_x2
    call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x2,ja_x2,ia_x2,a_tmp,ja_tmp, &
      & ia_tmp,nnz,i_err)

    ! Check for error
    call check_sparsekit_error(i_err,message)

    ! Deallocate and allocate memory for new combination
    iaS = 0
    deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
    deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating jaS and or aS in &
        & csr_combine_pointers(): stopping'
      stop
    endif

    ! Assign temporary values to iaS and jaS
    iaS(1) = ia_tmp(1)
    do i = 1,n_rows_0
      iaS(i+1) = ia_tmp(i+1)
      do j = ia_tmp(i),ia_tmp(i+1)-1
        jaS(j) = ja_tmp(j)
      enddo
    enddo


    if(ndim == 3) then
      ! Count number of nonzero element in (a_0 + a_x1 + a_x2) + a_x3
      ! -------------------------------------------------------------
      n_degr = 0
      iw = 0
      call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x3,ia_x3,n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      ia_tmp = 0
      deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
      deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp    
    
      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Combine (a_0 + a_x1 + a_x2) with a_x3
      call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x3,ja_x3,ia_x3,a_tmp,ja_tmp, &
        & ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      iaS = 0
      deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
      deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp
      
      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating jaS and or aS in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Assign temporary values to iaS and jaS
      iaS(1) = ia_tmp(1)
      do i = 1,n_rows_0
        iaS(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i),ia_tmp(i+1)-1
          jaS(j) = ja_tmp(j)
        enddo
      enddo 

    endif


    if(viscous) then
      ! Count number of nonzero element
      ! -------------------------------
      n_degr = 0
      iw = 0
      call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x11,ia_x11,n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      ia_tmp = 0
      deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
      deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp   

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Combine terms
      call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x11,ja_x11,ia_x11,a_tmp, &
        & ja_tmp,ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      iaS = 0
      deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
      deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating jaS and or aS in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Assign temporary values to iaS and jaS
      iaS(1) = ia_tmp(1)
      do i = 1,n_rows_0
        iaS(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i),ia_tmp(i+1)-1
          jaS(j) = ja_tmp(j)
        enddo
      enddo 

      ! Count number of nonzero element 
      ! -------------------------------
      n_degr = 0
      iw = 0
      call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x22,ia_x22,n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      ia_tmp = 0
      deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
      deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp   

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Combine terms
      call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x22,ja_x22,ia_x22,a_tmp, &
        & ja_tmp,ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      iaS = 0
      deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
      deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp
      
      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating jaS and or aS in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Assign temporary values to iaS and jaS
      iaS(1) = ia_tmp(1)
      do i = 1,n_rows_0
        iaS(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i),ia_tmp(i+1)-1
          jaS(j) = ja_tmp(j)
        enddo
      enddo 

      if(ndim == 3) then
        ! Count number of nonzero element
        ! -------------------------------
        n_degr = 0
        iw = 0
        call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x33,ia_x33,n_degr,nnz,iw) 

        ! Deallocate and allocate memory for new combination
        ia_tmp = 0
        deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
        deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp   

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
            & csr_combine_pointers(): stopping'
          stop
        endif

        ! Combine terms
        call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x33,ja_x33,ia_x33,a_tmp, &
          & ja_tmp,ia_tmp,nnz,i_err) 

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Deallocate and allocate memory for new combination
        iaS = 0
        deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
        deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating jaS and or aS in &
            & csr_combine_pointers(): stopping'
          stop
        endif

        ! Assign temporary values to iaS and jaS
        iaS(1) = ia_tmp(1)
        do i = 1,n_rows_0
          iaS(i+1) = ia_tmp(i+1)
          do j = ia_tmp(i),ia_tmp(i+1)-1
            jaS(j) = ja_tmp(j)
          enddo
        enddo 

      endif

      if(crossterms) then
        ! Count number of nonzero element
        ! -------------------------------
        n_degr = 0
        iw = 0
        call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x12,ia_x12,n_degr,nnz,iw) 

        ! Deallocate and allocate memory for new combination
        ia_tmp = 0
        deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
        deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp   

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
            & csr_combine_pointers(): stopping'
          stop
        endif

        ! Combine terms
        call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x12,ja_x12,ia_x12,a_tmp, &
          & ja_tmp,ia_tmp,nnz,i_err) 
        
        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Deallocate and allocate memory for new combination
        iaS = 0
        deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
        deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating jaS and or aS in &
            & csr_combine_pointers(): stopping'
          stop
        endif

        ! Assign temporary values to iaS and jaS
        iaS(1) = ia_tmp(1)
        do i = 1,n_rows_0
          iaS(i+1) = ia_tmp(i+1)
          do j = ia_tmp(i),ia_tmp(i+1)-1
            jaS(j) = ja_tmp(j)
          enddo
        enddo 

        if(ndim == 3) then
          ! Count number of nonzero element
          ! -------------------------------
          n_degr = 0
          iw = 0
          call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x13,ia_x13,n_degr,nnz,iw) 

          ! Deallocate and allocate memory for new combination
          ia_tmp = 0
          deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
          deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp 
          
          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
              & csr_combine_pointers(): stopping'
            stop
          endif

          ! Combine terms
          call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x13,ja_x13,ia_x13,a_tmp, &
            & ja_tmp,ia_tmp,nnz,i_err) 

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Deallocate and allocate memory for new combination
          iaS = 0
          deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
          deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating jaS and or aS in &
              & csr_combine_pointers(): stopping'
            stop
          endif

          ! Assign temporary values to iaS and jaS
          iaS(1) = ia_tmp(1)
          do i = 1,n_rows_0
            iaS(i+1) = ia_tmp(i+1)
            do j = ia_tmp(i),ia_tmp(i+1)-1
              jaS(j) = ja_tmp(j)
            enddo
          enddo 

          ! Count number of nonzero element
          ! -------------------------------
          n_degr = 0
          iw = 0
          call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_x23,ia_x23,n_degr,nnz,iw) 

          ! Deallocate and allocate memory for new combination
          ia_tmp = 0
          deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
          deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp 

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
              & csr_combine_pointers(): stopping'
            stop
          endif

          ! Combine terms
          call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_x23,ja_x23,ia_x23,a_tmp, &
            & ja_tmp,ia_tmp,nnz,i_err) 

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Deallocate and allocate memory for new combination
          iaS = 0
          deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
          deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating jaS and or aS in &
              & csr_combine_pointers(): stopping'
            stop
          endif

          ! Assign temporary values to iaS and jaS
          iaS(1) = ia_tmp(1)
          do i = 1,n_rows_0
            iaS(i+1) = ia_tmp(i+1)
            do j = ia_tmp(i),ia_tmp(i+1)-1
              jaS(j) = ja_tmp(j)
            enddo
          enddo 

        endif ! End if 3D

      endif ! End cross terms

    endif ! End viscous


    if(IMEX_penalty == 'implicit') then

      ! Inviscid penalty contributions
      ! Count number of nonzero element
      ! -------------------------------

      n_degr = 0 ; iw     = 0

      call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_penI,ia_penI,n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      ia_tmp = 0
      deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
      deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp 
      
      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Combine terms
      call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_penI,ja_penI,ia_penI,a_tmp, &
        & ja_tmp,ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      iaS = 0
      deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
      deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating jaS and or aS in &
          & csr_combine_pointers(): stopping'
        stop
      endif

      ! Assign temporary values to iaS and jaS
      iaS(1) = ia_tmp(1)
      do i = 1,n_rows_0
        iaS(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i),ia_tmp(i+1)-1
          jaS(j) = ja_tmp(j)
        enddo
      enddo 

      ! Viscous penalty contributions
      if(viscous) then

        ! Count number of nonzero element
        ! -------------------------------
        n_degr = 0
        iw = 0
        call aplbdg(n_rows_0,ngrowsS,jaS,iaS,ja_penV,ia_penV,n_degr,nnz,iw) 
  
        ! Deallocate and allocate memory for new combination
        ia_tmp = 0
        deallocate(ja_tmp) ; allocate(ja_tmp(nnz),stat=status_allocation(2)) ; ja_tmp = 0
        deallocate( a_tmp) ; allocate( a_tmp(nnz),stat=status_allocation(3)) ;  a_tmp = 0.0_wp 
      
        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
            & csr_combine_pointers(): stopping'
          stop
        endif

        ! Combine terms
        call aplb1(n_rows_0,ngrowsS,0,aS,jaS,iaS,a_penV,ja_penV,ia_penV,a_tmp, &
          & ja_tmp,ia_tmp,nnz,i_err) 
  
        ! Check for error
        call check_sparsekit_error(i_err,message)
  
        ! Deallocate and allocate memory for new combination
        iaS = 0
        deallocate(jaS) ; allocate(jaS(nnz),stat=status_allocation(2)) ; jaS = 0
        deallocate( aS) ; allocate( aS(nnz),stat=status_allocation(3)) ;  aS = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating jaS and or aS in &
            & csr_combine_pointers(): stopping'
          stop
        endif

        ! Assign temporary values to iaS and jaS
        iaS(1) = ia_tmp(1)
        do i = 1,n_rows_0
          iaS(i+1) = ia_tmp(i+1)
          do j = ia_tmp(i),ia_tmp(i+1)-1
            jaS(j) = ja_tmp(j)
          enddo
        enddo 

      endif  !  End Viscous Contribution

    endif ! IMEX_penalty

    ! Deallocate temporary memory
    deallocate(ia_tmp,ja_tmp,a_tmp,n_degr,iw)

    ! Initialize i_node and shift
    i_node = 0
    shift = 0

    ! Assign term column pointers ka_j() that point to the jaS() matrix
    do i_elem = 1, n_elems
      do ll = 1, nodesperelem

        i_node = i_node + 1

        ! Shift for the ka_j pointers because the module
        ! jacobian_matrix_implicit_ts constructs the Jacobian matrix for a
        ! single element or cell. This shift makes the ka_j to start from 1.
        shift = iaS((i_elem-1)*nodesperelem + 1) - 1  

        totalja:do k = iaS(i_node),iaS(i_node+1) - 1
          do j = ia_0(i_node),ia_0(i_node+1) - 1
            if(ja_0(j) == jaS(k)) ka_0(j) = k - shift
          enddo

          if(ndim >= 1) then
            do j = ia_x1(i_node),ia_x1(i_node+1) - 1
              if(ja_x1(j) == jaS(k)) ka_x1(j) = k - shift
            enddo
            if(ndim >= 2) then
              do j = ia_x2(i_node),ia_x2(i_node+1) - 1
                if(ja_x2(j) == jaS(k)) ka_x2(j) = k - shift
              enddo
              if(ndim == 3) then
                do j = ia_x3(i_node),ia_x3(i_node+1) - 1
                  if(ja_x3(j) == jaS(k)) ka_x3(j) = k - shift
                enddo
              endif
            endif
          endif

          if(viscous) then
            if(ndim >= 1) then
              do j = ia_x11(i_node),ia_x11(i_node+1) - 1
                if(ja_x11(j) == jaS(k)) ka_x11(j) = k - shift
              enddo
              if(ndim >= 2) then
                do j = ia_x22(i_node),ia_x22(i_node+1) - 1
                  if(ja_x22(j) == jaS(k)) ka_x22(j) = k - shift
                enddo
                if(ndim == 3) then
                  do j = ia_x33(i_node),ia_x33(i_node+1) - 1
                    if(ja_x33(j) == jaS(k)) ka_x33(j) = k - shift
                  enddo
                endif
              endif
            endif
            if(crossterms) then
              if(ndim >= 2) then
                do j = ia_x12(i_node),ia_x12(i_node+1) - 1
                  if(ja_x12(j) == jaS(k)) ka_x12(j) = k - shift
                enddo
              endif
              if(ndim == 3) then
                do j = ia_x13(i_node),ia_x13(i_node+1) - 1
                  if(ja_x13(j) == jaS(k)) ka_x13(j) = k - shift
                enddo
                do j = ia_x23(i_node),ia_x23(i_node+1) - 1
                  if(ja_x23(j) == jaS(k)) ka_x23(j) = k - shift
                enddo
              endif
            endif
          endif

        enddo totalja
      enddo

    enddo

!   i_node = 0
!   do i_elem = 1, n_elems
!     do ll = 1, nodesperelem
!       i_node = i_node + 1
!       shift = iaS((i_elem-1)*nodesperelem + 1) - 1  
!        totaljT:do k = iaS(i_node),iaS(i_node+1) - 1
!        
!         if(IMEX_penalty == 'implicit') then
!           do j = ia_penI_proc(i_node),ia_penI_proc(i_node+1) - 1
!             cntI = la_penI_proc(j)
!             if(ja_penI_proc(j) == jaS(k)) ka_penI_proc(cntI) = k - shift
!           enddo
!           if(viscous) then
!             do j = ia_penV_proc(i_node),ia_penV_proc(i_node+1) - 1
!               cntV = la_penV_proc(j) 
!               if(ja_penV_proc(j) == jaS(k)) ka_penV_proc(cntV) = k - shift
!             enddo
!           endif
!         endif

!       enddo totaljT
!     enddo
!   enddo

    if(IMEX_penalty == 'implicit') then

      cntI = 0 ;
      elemloopI: do i_elem = low_elem, high_elem
        faceloopI:do iface = 1,nfacesperelem

          shift = iaS((i_elem-low_elem)*nodesperelem + 1) - 1  

          if (ef2e(1,iface,i_elem) < 0) then
            do jnode = 1,nodesperface
              inode = kfacenodes(jnode,iface)
              knode = inode
              kelem = i_elem
              local_row = (i_elem-low_elem)*nodesperelem + inode

              ! On process
              cntI = cntI + 1 ; jj   = (i_elem-1)*nodesperelem + inode

              do k = iaS(local_row),iaS(local_row+1)-1
                if(jaS(k) == jj) then
                  ka_penI_proc(cntI) = k - shift
                endif
              enddo

              ! Off process
              cntI = cntI + 1 ; jj   = (kelem-1)*nodesperelem + knode
              do k = iaS(local_row),iaS(local_row+1)-1
                if(jaS(k) == jj) then
                  ka_penI_proc(cntI) = k - shift
                endif
              enddo

            enddo

          else if (ef2e(3,iface,i_elem) /= myprocid) then
            do i = 1, nodesperface
              jnode = nodesperface*(iface-1)+i
              inode = ifacenodes(jnode)
              knode = efn2efn(1,jnode,i_elem)
              kelem = efn2efn(2,jnode,i_elem)
              local_row = (i_elem-low_elem)*nodesperelem + inode

              ! On process
              cntI = cntI + 1 ; jj   = (i_elem-1)*nodesperelem + inode
              do k = iaS(local_row),iaS(local_row+1)-1
                if(jaS(k) == jj) then
                  ka_penI_proc(cntI) = k - shift
                endif
              enddo

              ! Off process
              cntI = cntI + 1 ; jj   = (kelem-1)*nodesperelem + knode
              do k = iaS(local_row),iaS(local_row+1)-1
                if(jaS(k) == jj) then
                  ka_penI_proc(cntI) = k - shift
                endif
              enddo

            enddo

          else
            do i = 1,nodesperface
              jnode = nodesperface*(iface-1)+i
              inode = ifacenodes(jnode)
              knode = efn2efn(1,jnode,i_elem)
              kelem = efn2efn(2,jnode,i_elem)
              local_row = (i_elem-low_elem)*nodesperelem + inode

              ! On process
              cntI = cntI + 1 ; jj   = (i_elem-1)*nodesperelem + inode
              do k = iaS(local_row),iaS(local_row+1)-1
                if(jaS(k) == jj) then
                  ka_penI_proc(cntI) = k - shift
                endif
              enddo

              ! Off process
              cntI = cntI + 1 ; jj   = (kelem-1)*nodesperelem + knode
              do k = iaS(local_row),iaS(local_row+1)-1
                if(jaS(k) == jj) then
                  ka_penI_proc(cntI) = k - shift
                endif
              enddo

            enddo
          end if
        enddo faceloopI
      enddo elemloopI

      if(viscous) then
        cntV = 0 ;
        ka_penV_proc(:) = 0

        elemloopV: do i_elem = low_elem, high_elem
          shift = iaS((i_elem-low_elem)*nodesperelem + 1) - 1  
          faceloopV:do iface = 1,nfacesperelem


            if (ef2e(1,iface,i_elem) < 0) then
              do jnode = 1,nodesperface
                inode = kfacenodes(jnode,iface)
                knode = inode
                kelem = i_elem
                local_row = (i_elem-low_elem)*nodesperelem + inode
  
                ! On process
                do jdir = 1,ndim
  
                  do ii = iagrad(inode), iagrad(inode+1)-1
          
                    cntV = cntV + 1 ; jj   = jagrad(jdir,ii) + (i_elem-1) * nodesperelem
                    do k = iaS(local_row),iaS(local_row+1)-1
                      if(jaS(k) == jj) then
                        ka_penV_proc(cntV) = k - shift
                      endif
                    enddo
                  end do
  
                end do
  
                ! Off process
                do jdir = 1,ndim
          
                  do ii = iagrad(knode), iagrad(knode+1)-1
          
                    cntV = cntV + 1 ; jj   = jagrad(jdir,ii) + (kelem-1) * nodesperelem
                    do k = iaS(local_row),iaS(local_row+1)-1
                      if(jaS(k) == jj) then
                        ka_penV_proc(cntV) = k - shift
                      endif
                    enddo
                  end do
                end do
  
              enddo
  
            else if (ef2e(3,iface,i_elem) /= myprocid) then
              do i = 1, nodesperface
                jnode = nodesperface*(iface-1)+i
                inode = ifacenodes(jnode)
                knode = efn2efn(1,jnode,i_elem)
                kelem = efn2efn(2,jnode,i_elem)
                local_row = (i_elem-low_elem)*nodesperelem + inode

                ! On element  (use inode and i_elem)
                do jdir = 1,ndim
  
                  do ii = iagrad(inode), iagrad(inode+1)-1
          
                    cntV = cntV + 1 ; jj   = jagrad(jdir,ii) + (i_elem-1) * nodesperelem
                    do k = iaS(local_row),iaS(local_row+1)-1
                      if(jaS(k) == jj) then
                        ka_penV_proc(cntV) = k - shift
                      endif
                    enddo
                  end do
  
                end do

                ! Off element  (use knode and kelem)
                do jdir = 1,ndim
          
                  do ii = iagrad(knode), iagrad(knode+1)-1
          
                    cntV = cntV + 1 ; jj   = jagrad(jdir,ii) + (kelem-1) * nodesperelem
                    do k = iaS(local_row),iaS(local_row+1)-1
                      if(jaS(k) == jj) then
                        ka_penV_proc(cntV) = k - shift
                      endif
                    enddo
                  end do
  
                end do
  
              enddo
  
            else
              do i = 1,nodesperface
                jnode = nodesperface*(iface-1)+i
                inode = ifacenodes(jnode)
                knode = efn2efn(1,jnode,i_elem)
                kelem = efn2efn(2,jnode,i_elem)
                local_row = (i_elem-low_elem)*nodesperelem + inode
  
                ! On element  (use inode and i_elem)
                do jdir = 1,ndim
  
                  do ii = iagrad(inode), iagrad(inode+1)-1
          
                    cntV = cntV + 1 ; jj   = jagrad(jdir,ii) + (i_elem-1) * nodesperelem
                    do k = iaS(local_row),iaS(local_row+1)-1
                      if(jaS(k) == jj) then
                        ka_penV_proc(cntV) = k - shift
                      endif
                    enddo
                  end do
  
                end do
  
                ! Off element  (use knode and kelem)
                do jdir = 1,ndim
          
                  do ii = iagrad(knode), iagrad(knode+1)-1
          
                    cntV = cntV + 1 ; jj   = jagrad(jdir,ii) + (kelem-1) * nodesperelem
                    do k = iaS(local_row),iaS(local_row+1)-1
                      if(jaS(k) == jj) then
                        ka_penV_proc(cntV) = k - shift
                      endif
                    enddo
                  end do
  
                end do
              enddo
            end if
          enddo faceloopV
        enddo elemloopV
      endif
    endif

    return
  end subroutine csr_combine_pointers

  !============================================================================
   
  !============================================================================
  ! csr_combine_pointers_element - combines the CSR format of the 
  !                                differentiation matrices for one element.

  subroutine csr_combine_pointers_element()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use nsereferencevariables
    use unary_mod
    use jacobian_matrix_implicit_ts_variables
    use controlvariables, only : imex_element, imex_penalty

    ! Nothing is implicitly defined
    implicit none

    integer :: i,j
    integer , allocatable , dimension(:) :: n_degr
    integer , allocatable , dimension(:) :: iw
    integer :: nnz, n_rows, n_cols
    integer , dimension(10) :: status_allocation
    integer , allocatable , dimension(:) :: ia_tmp
    integer , allocatable , dimension(:) :: ja_tmp
    real(wp) , allocatable , dimension(:) :: a_tmp
    integer :: i_err
    real(wp) , allocatable , dimension(:) :: a_elem
    character(120) :: message

    continue

    ! Subroutine name for possible error message
    message = 'csr_combine_pointers_element'

    ! Number of rows and columns 
    n_rows = size(ia_0_elem)-1
    n_cols = size(ia_0_elem)-1

    ! Initialize i_err to normal value
    i_err = 0

    ! Allocate memory for  n_degr and iw (unary_mod.f90 variables)   
    status_allocation = 0
    allocate(n_degr(n_rows),stat=status_allocation(1))
    allocate(iw(n_cols),stat=status_allocation(2))
    n_degr = 0
    iw = 0

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating n_degr and or iw in &
        & csr_combine_pointers_element(): stopping' 
      stop
    endif


    ! Count number of nonzero element in (a_0_elem + a_x1_elem)
    ! ---------------------------------------------------------
    call aplbdg(n_rows,n_cols,ja_0_elem,ia_0_elem,ja_x1_elem, &
      & ia_x1_elem,n_degr,nnz,iw)

    ! Allocate memory for temporary arrays
    allocate(ia_tmp(size(ia_0_elem)),stat=status_allocation(1))
    allocate(ja_tmp(nnz),stat=status_allocation(2))
    allocate(a_tmp(nnz),stat=status_allocation(3))
    ia_tmp = 0
    ja_tmp = 0
    a_tmp = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating ia_tmp and or ja_tmp and or a_tmp in &
        & csr_combine_pointers_element(): stopping'
      stop
    endif

    ! Combine a_0_elem with a_x1_elem
    call aplb1(n_rows,n_cols,0,a_0_elem,ja_0_elem,ia_0_elem,a_x1_elem, &
      & ja_x1_elem,ia_x1_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err)

    ! Check for error
    call check_sparsekit_error(i_err,message)

    ! Allocate memory for temporary arrays
    allocate(ia_elem(size(ia_0_elem)),stat=status_allocation(1))
    allocate(ja_elem(nnz),stat=status_allocation(2))
    allocate(a_elem(nnz),stat=status_allocation(3))
    ia_elem = 0
    ja_elem = 0
    a_elem = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating ia_elem and or ja_elem and or a_elem &
        & in csr_combine_pointers_element(): stopping'
      stop
    endif

    ! Assign temporary values to ia_elem and ja_elem
    ia_elem(1) = ia_tmp(1)
    do i = 1, n_rows
      ia_elem(i+1) = ia_tmp(i+1)
      do j = ia_tmp(i), ia_tmp(i+1)-1
        ja_elem(j) = ja_tmp(j)
      enddo
    enddo

    ! Count number of nonzero element in (a_0_elem + a_x1_elem) + a_x2_elem
    ! ---------------------------------------------------------------------
    n_degr = 0
    iw = 0
    call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x2_elem,ia_x2_elem, &
      & n_degr,nnz,iw) 

    ! Deallocate and allocate memory for new combination
    deallocate(ja_tmp)
    allocate(ja_tmp(nnz),stat=status_allocation(2))
    deallocate(a_tmp)
    allocate(a_tmp(nnz),stat=status_allocation(3))
    ia_tmp = 0
    ja_tmp = 0
    a_tmp = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
        & csr_combine_pointers_element(): stopping'
      stop
    endif

    ! Combine (a_0_elem + a_x1_elem) with a_x2_elem
    call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x2_elem, &
      & ja_x2_elem,ia_x2_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

    ! Check for error
    call check_sparsekit_error(i_err,message)

    ! Deallocate and allocate memory for new combination
    deallocate(ja_elem)
    allocate(ja_elem(nnz),stat=status_allocation(2))
    deallocate(a_elem)
    allocate(a_elem(nnz),stat=status_allocation(3))
    ia_elem = 0
    ja_elem = 0
    a_elem = 0.0_wp

    ! Check allocation status
    if(sum(status_allocation) > 0) then
      write(*,*) 'Failure in allocating ja_element and or a_element in &
        & csr_combine_pointers_element(): stopping'
      stop
    endif

    ! Assign temporary values to ia_elem and ja_elem
    ia_elem(1) = ia_tmp(1)
    do i = 1, n_rows
      ia_elem(i+1) = ia_tmp(i+1)
      do j = ia_tmp(i), ia_tmp(i+1)-1
        ja_elem(j) = ja_tmp(j)
      enddo
    enddo


    if(ndim == 3) then
      ! Count number of nonzero element in (a_0_elem + a_x1_elem + a_x2_elem) 
      ! + a_x3_elem
      ! ---------------------------------------------------------------------
      n_degr = 0
      iw = 0
      call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x3_elem,ia_x3_elem, &
        & n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      deallocate(ja_tmp)
      allocate(ja_tmp(nnz),stat=status_allocation(2))
      deallocate(a_tmp)
      allocate(a_tmp(nnz),stat=status_allocation(3))
      ia_tmp = 0
      ja_tmp = 0
      a_tmp = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers_element(): stopping'
        stop
      endif

      ! Combine (a_0_elem + a_x1_elem + a_x2_elem) with a_x3_elem
      call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x3_elem, &
        & ja_x3_elem,ia_x3_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      deallocate(ja_elem)
      allocate(ja_elem(nnz),stat=status_allocation(2))
      deallocate(a_elem)
      allocate(a_elem(nnz),stat=status_allocation(3))
      ia_elem = 0
      ja_elem = 0
      a_elem = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_elem and or a_elem in &
          & csr_combine_pointers_element(): stopping'
        stop
      endif

      ! Assign temporary values to ia_elem and ja_elem
      ia_elem(1) = ia_tmp(1)
      do i = 1, n_rows
        ia_elem(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i), ia_tmp(i+1)-1
          ja_elem(j) = ja_tmp(j)
        enddo
      enddo 

    endif


    if(viscous) then
      ! Count number of nonzero element
      ! -------------------------------
      n_degr = 0
      iw = 0
      call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x11_elem,ia_x11_elem, &
        & n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      deallocate(ja_tmp)
      allocate(ja_tmp(nnz),stat=status_allocation(2))
      deallocate(a_tmp)
      allocate(a_tmp(nnz),stat=status_allocation(3))
      ia_tmp = 0
      ja_tmp = 0
      a_tmp = 0.0_wp
      
      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers_element(): stopping'
        stop
      endif

      ! Combine terms
      call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x11_elem, &
        & ja_x11_elem,ia_x11_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      deallocate(ja_elem)
      allocate(ja_elem(nnz),stat=status_allocation(2))
      deallocate(a_elem)
      allocate(a_elem(nnz),stat=status_allocation(3))
      ia_elem = 0
      ja_elem = 0
      a_elem = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_elem and or a_elem in &
          & csr_combine_pointers_element(): stopping'
        stop
      endif

      ! Assign temporary values to ia_elem and ja_elem
      ia_elem(1) = ia_tmp(1)
      do i = 1, n_rows
        ia_elem(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i), ia_tmp(i+1)-1
          ja_elem(j) = ja_tmp(j)
        enddo
      enddo 


      ! Count number of nonzero element 
      ! -------------------------------
      n_degr = 0
      iw = 0
      call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x22_elem,ia_x22_elem, &
        & n_degr,nnz,iw) 

      ! Deallocate and allocate memory for new combination
      deallocate(ja_tmp)
      allocate(ja_tmp(nnz),stat=status_allocation(2))
      deallocate(a_tmp)
      allocate(a_tmp(nnz),stat=status_allocation(3))
      ia_tmp = 0
      ja_tmp = 0
      a_tmp = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
          & csr_combine_pointers_element(): stopping'
        stop
      endif

      ! Combine terms
      call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x22_elem, &
        & ja_x22_elem,ia_x22_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Deallocate and allocate memory for new combination
      deallocate(ja_elem)
      allocate(ja_elem(nnz),stat=status_allocation(2))
      deallocate(a_elem)
      allocate(a_elem(nnz),stat=status_allocation(3))
      ia_elem = 0
      ja_elem = 0
      a_elem = 0.0_wp

      ! Check allocation status
      if(sum(status_allocation) > 0) then
        write(*,*) 'Failure in allocating ja_elem and or a_elem in &
          & csr_combine_pointers_element(): stopping'
        stop
      endif

      ! Assign temporary values to ia_elem and ja_elem
      ia_elem(1) = ia_tmp(1)
      do i = 1, n_rows
        ia_elem(i+1) = ia_tmp(i+1)
        do j = ia_tmp(i), ia_tmp(i+1)-1
          ja_elem(j) = ja_tmp(j)
        enddo
      enddo 

      if(ndim == 3) then
        ! Count number of nonzero element
        ! -------------------------------
        n_degr = 0
        iw = 0
        call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x33_elem,ia_x33_elem, &
          & n_degr,nnz,iw) 

        ! Deallocate and allocate memory for new combination
        deallocate(ja_tmp)
        allocate(ja_tmp(nnz),stat=status_allocation(2))
        deallocate(a_tmp)
        allocate(a_tmp(nnz),stat=status_allocation(3))
        ia_tmp = 0
        ja_tmp = 0
        a_tmp = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
            & csr_combine_pointers_element(): stopping'
          stop
        endif

        ! Combine terms
        call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x33_elem, &
          & ja_x33_elem,ia_x33_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Deallocate and allocate memory for new combination
        deallocate(ja_elem)
        allocate(ja_elem(nnz),stat=status_allocation(2))
        deallocate(a_elem)
        allocate(a_elem(nnz),stat=status_allocation(3))
        ia_elem = 0
        ja_elem = 0
        a_elem = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_elem and or a_elem in &
            & csr_combine_pointers_element(): stopping'
          stop
        endif

        ! Assign temporary values to ia_elem and ja_elem
        ia_elem(1) = ia_tmp(1)
        do i = 1, n_rows
          ia_elem(i+1) = ia_tmp(i+1)
          do j = ia_tmp(i), ia_tmp(i+1)-1
            ja_elem(j) = ja_tmp(j)
          enddo
        enddo 

      endif

      if(crossterms) then
        ! Count number of nonzero element
        ! -------------------------------
        n_degr = 0
        iw = 0
        call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x12_elem,ia_x12_elem, &
          & n_degr,nnz,iw) 

        ! Deallocate and allocate memory for new combination
        deallocate(ja_tmp)
        allocate(ja_tmp(nnz),stat=status_allocation(2))
        deallocate(a_tmp)
        allocate(a_tmp(nnz),stat=status_allocation(3))
        ia_tmp = 0
        ja_tmp = 0
        a_tmp = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
            & csr_combine_pointers_element(): stopping'
          stop
        endif

        ! Combine terms
        call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x12_elem, &
          & ja_x12_elem,ia_x12_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Deallocate and allocate memory for new combination
        deallocate(ja_elem)
        allocate(ja_elem(nnz),stat=status_allocation(2))
        deallocate(a_elem)
        allocate(a_elem(nnz),stat=status_allocation(3))
        ia_elem = 0
        ja_elem = 0
        a_elem = 0.0_wp

        ! Check allocation status
        if(sum(status_allocation) > 0) then
          write(*,*) 'Failure in allocating ja_elem and or a_elem in &
            & csr_combine_pointers_element(): stopping'
          stop
        endif

        ! Assign temporary values to ia_elem and ja_elem
        ia_elem(1) = ia_tmp(1)
        do i = 1, n_rows
          ia_elem(i+1) = ia_tmp(i+1)
          do j = ia_tmp(i), ia_tmp(i+1)-1
            ja_elem(j) = ja_tmp(j)
          enddo
        enddo 

        if(ndim == 3) then
          ! Count number of nonzero element
          ! -------------------------------
          n_degr = 0
          iw = 0
          call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x13_elem,ia_x13_elem, &
            & n_degr,nnz,iw) 

          ! Deallocate and allocate memory for new combination
          deallocate(ja_tmp)
          allocate(ja_tmp(nnz),stat=status_allocation(2))
          deallocate(a_tmp)
          allocate(a_tmp(nnz),stat=status_allocation(3))
          ia_tmp = 0
          ja_tmp = 0
          a_tmp = 0.0_wp

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
              & csr_combine_pointers_element(): stopping'
            stop
          endif

          ! Combine terms
          call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x13_elem, &
            & ja_x13_elem,ia_x13_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Deallocate and allocate memory for new combination
          deallocate(ja_elem)
          allocate(ja_elem(nnz),stat=status_allocation(2))
          deallocate(a_elem)
          allocate(a_elem(nnz),stat=status_allocation(3))
          ia_elem = 0
          ja_elem = 0
          a_elem = 0.0_wp

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating ja_elem and or a_elem in &
              & csr_combine_pointers_element(): stopping'
            stop
          endif

          ! Assign temporary values to ia_elem and ja_elem
          ia_elem(1) = ia_tmp(1)
          do i = 1, n_rows
            ia_elem(i+1) = ia_tmp(i+1)
            do j = ia_tmp(i), ia_tmp(i+1)-1
              ja_elem(j) = ja_tmp(j)
            enddo
          enddo 


          ! Count number of nonzero element
          ! -------------------------------
          n_degr = 0
          iw = 0
          call aplbdg(n_rows,n_cols,ja_elem,ia_elem,ja_x23_elem,ia_x23_elem, &
            & n_degr,nnz,iw) 

          ! Deallocate and allocate memory for new combination
          deallocate(ja_tmp)
          allocate(ja_tmp(nnz),stat=status_allocation(2))
          deallocate(a_tmp)
          allocate(a_tmp(nnz),stat=status_allocation(3))
          ia_tmp = 0
          ja_tmp = 0
          a_tmp = 0.0_wp

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating ja_tmp and or a_tmp in &
              & csr_combine_pointers_element(): stopping'
            stop
          endif

          ! Combine terms
          call aplb1(n_rows,n_cols,0,a_elem,ja_elem,ia_elem,a_x23_elem, &
            & ja_x23_elem,ia_x23_elem,a_tmp,ja_tmp,ia_tmp,nnz,i_err) 

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Deallocate and allocate memory for new combination
          deallocate(ja_elem)
          allocate(ja_elem(nnz),stat=status_allocation(2))
          deallocate(a_elem)
          allocate(a_elem(nnz),stat=status_allocation(3))
          ia_elem = 0
          ja_elem = 0
          a_elem = 0.0_wp

          ! Check allocation status
          if(sum(status_allocation) > 0) then
            write(*,*) 'Failure in allocating ja_elem and or a_elem in &
              & csr_combine_pointers_element(): stopping'
            stop
          endif

          ! Assign temporary values to ia_elem and ja_elem
          ia_elem(1) = ia_tmp(1)
          do i = 1, n_rows
            ia_elem(i+1) = ia_tmp(i+1)
            do j = ia_tmp(i), ia_tmp(i+1)-1
              ja_elem(j) = ja_tmp(j)
            enddo
          enddo 

        endif ! End if 3D

      endif ! End cross terms

    endif ! End viscous

    deallocate(ia_tmp,ja_tmp,a_tmp,n_degr,iw)

    return
  end subroutine csr_combine_pointers_element

  !============================================================================

  !============================================================================
  ! csr_on_element_operator - Constructs the CSR formats for computing the
  ! residual Jacobian on element using a matrix-matrix multiply approach. First 
  ! an element "brick" is constructed; then it is added to the array which
  ! holds the CSR differentiation matrix of all the elements owned by the
  ! process.

  subroutine csr_on_element_operator()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use nsereferencevariables
    use unary_mod
    use jacobian_matrix_implicit_ts_variables
    use controlvariables, only : imex_element, imex_penalty

    ! Nothing is implicitly defined
    implicit none
    
    integer :: i, j, k, l
    integer :: i_node
    integer :: i_err
    character(120) :: message
    real(wp), allocatable, dimension(:,:) :: eye_matrix
    real(wp), allocatable, dimension(:,:,:) :: dummy_matrix_elem
    integer :: low_elem, high_elem
    integer :: cnt, nnz, n_elems
    integer :: i_elem
    integer :: ll, shift
    integer :: global_shift

    continue

    ! Preliminary check on the dimensions in each direction
    if (ndim .eq. 2) then
      if (size(ia_x1_elem) .ne. size(ia_x2_elem)) then
        write(*,*) 'This routine is written assuming that the order of the &
            & polynomial reconstruction in the xi and eta directions is the &
            & same and all the cell use such a reconstruction! If this is not &
            & the case then a huge portion of the code needs to be re-written.'
        write(*,*) 'Stopping...'
        stop
      endif
    else
      if (size(ia_x1_elem) .ne. size(ia_x2_elem) .or. size(ia_x1_elem) .ne. &
        & size(ia_x3_elem)) then
        write(*,*) 'This routine is written assuming that the order of the &
            & polynomial reconstruction in the xi, eta and zeta directions is &
            & the same and all the cell use such a reconstruction! If this is &
            & not the case then a huge portion of the code needs to be &
            & re-written.'
        write(*,*) 'Stopping...'
        stop
      endif
    endif

    ! Low volumetric element index
    low_elem = ihelems(1)

    ! High volumetric element index
    high_elem = ihelems(2)

    ! Number of elements owned by the processor
    n_elems = 1 + high_elem - low_elem

    ! Subroutine name for possible error message
    message = 'csr_on_element_operator'

    ! Shift for global indexing for parallel computation
    global_shift = (low_elem-1)*nodesperelem

    ! Diagonal terms stored in ia_0
    cnt = 0
    nnz = 0
    do l = 1, n_elems
      cnt = cnt + nodesperelem
      nnz = nnz + nodesperelem
    enddo

    allocate(ia_0_matmul(cnt+1))    ;   ia_0_matmul = 0
    allocate(ja_0_matmul(nnz))      ;   ja_0_matmul = 0
    allocate(ka_0_matmul(nnz))      ;   ka_0_matmul = 0

    cnt = 1
    ia_0_matmul(cnt) = 1
    do i_elem = 1, n_elems
      do i = 1, nodesperelem
        ia_0_matmul(cnt+1) = ia_0_matmul(cnt) + 1
        ja_0_matmul(cnt) = cnt
        cnt = cnt + 1 
      enddo
    enddo

    ! Shift to global indexing for parallel computations 
    ja_0_matmul = ja_0_matmul + global_shift


    ! Construct a block identity matrix for the element
    ! -------------------------------------------------
    allocate(dummy_matrix_elem(nequations,nequations,nodesperelem))
    allocate(eye_matrix(nequations,nequations))
    ! Block identity matrix
    eye_matrix = 0.0_wp
    do i = 1, nequations
      eye_matrix(i,i) = 1.0_wp
    enddo
    ! Set the blocks identity matrices
    do i = 1, nodesperelem
      dummy_matrix_elem(:,:,i) = eye_matrix
    enddo
    
    ! Diagonal CSR format
    allocate(ja_diag_elem(nodesperelem))
    ja_diag_elem = 0
    
    allocate(ia_diag_elem(nodesperelem+1))
    ia_diag_elem = 0

    ia_diag_elem(1) = 1
    do i = 1, nodesperelem
      ia_diag_elem(i+1) = ia_diag_elem(i)+1
      ja_diag_elem(i) = i
    enddo
    

    ! x1
    ! ---
    ! Working array for the sparsekit.
    allocate(iw_inviscid(size(ia_x1_elem)-1))
    iw_inviscid = 0
 
    allocate(ia_1_matmul_elem(size(ia_x1_elem)))
    ia_1_matmul_elem = 0
    
    allocate(ja_1_matmul_elem(size(ja_x1_elem)))
    ja_1_matmul_elem = 0
    
    allocate(a_1_matmul_elem(nequations,nequations,size(a_x1_elem)))
    a_1_matmul_elem = 0.0_wp

    ! [D_x1] [dummy_matrix_elem]
    ! --------------------------
    ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of the
    ! element 
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
    a_1_matmul_elem = 0.0_wp
    iw_inviscid = 0
    call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
      & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
      & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
      & size(a_x1_elem),iw_inviscid,i_err)

    ! Check for error
    call check_sparsekit_error(i_err,message)

    ! Construct the array for all the elements owned by the processor
    cnt = 0
    nnz = 0
    do l = 1, n_elems
      cnt = cnt + nodesperelem
      nnz = nnz + nodesperelem*nodesperedge
    enddo

    allocate(ia_x1_matmul(cnt+1))
    ia_x1_matmul = 0
    
    allocate(ja_x1_matmul(nnz))
    ja_x1_matmul = 0

    allocate(ka_x1_matmul(nnz))
    ka_x1_matmul = 0

    nnz = 1
    cnt = 1
    ia_x1_matmul(1) = cnt
    do i_elem = 1, n_elems
      do i = 1, nodesperelem
        ia_x1_matmul(cnt) = ia_1_matmul_elem(i) + &
          & (i_elem-1)*nodesperelem*nodesperedge 
        do j = ia_1_matmul_elem(i), ia_1_matmul_elem(i+1) - 1
          ja_x1_matmul(nnz) = ja_1_matmul_elem(j) + (i_elem-1)*nodesperelem

          nnz = nnz + 1
        enddo 
        cnt = cnt + 1
      enddo 
      ia_x1_matmul(cnt) = ia_1_matmul_elem(nodesperelem+1) + &
        & (i_elem-1)*nodesperelem*nodesperedge
    end do
    
    ! Shift to global indexing for parallel computations 
    ja_x1_matmul = ja_x1_matmul + global_shift


    ! x2
    ! -- 
    ! [D_x2] [dummy_matrix_elem]
    ! --------------------------
    ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of the
    ! element 
    ia_1_matmul_elem = 0
    ja_1_matmul_elem = 0
    a_1_matmul_elem = 0.0_wp
    iw_inviscid = 0
    call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
      & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
      & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
      & size(a_x2_elem),iw_inviscid,i_err)

    ! Check for error
    call check_sparsekit_error(i_err,message)

    ! Construct the array for all the elements owned by the processor
    cnt = 0
    nnz = 0
    do l = 1, n_elems
      cnt = cnt + nodesperelem
      nnz = nnz + nodesperelem*nodesperedge
    enddo

    allocate(ia_x2_matmul(cnt+1))
    ia_x2_matmul = 0
    
    allocate(ja_x2_matmul(nnz))
    ja_x2_matmul = 0

    allocate(ka_x2_matmul(nnz))
    ka_x2_matmul = 0

    nnz = 1
    cnt = 1
    ia_x2_matmul(1) = cnt
    do i_elem = 1, n_elems
      do i = 1, nodesperelem
        ia_x2_matmul(cnt) = ia_1_matmul_elem(i) + &
          & (i_elem-1)*nodesperelem*nodesperedge 
        do j = ia_1_matmul_elem(i), ia_1_matmul_elem(i+1) - 1
          ja_x2_matmul(nnz) = ja_1_matmul_elem(j) + (i_elem-1)*nodesperelem

          nnz = nnz + 1
        enddo 
        cnt = cnt + 1
      enddo 
      ia_x2_matmul(cnt) = ia_1_matmul_elem(nodesperelem+1) + &
        & (i_elem-1)*nodesperelem*nodesperedge
    end do

    ! Shift to global indexing for parallel computations 
    ja_x2_matmul = ja_x2_matmul + global_shift


    if (viscous) then
      ! x11
      ! --
      ! Working array for the sparsekit.
      allocate(iw_viscous_1(size(ia_x11_elem)-1))
      iw_viscous_1 = 0

      allocate(ia_2_matmul_elem(size(ia_x1_elem)))
      ia_2_matmul_elem = 0
    
      allocate(ja_2_matmul_elem(size(ja_x1_elem)))
      ja_2_matmul_elem = 0
    
      allocate(a_2_matmul_elem(nequations,nequations,size(a_x1_elem)))
      a_2_matmul_elem = 0.0_wp
 
      allocate(ia_lap_matmul_elem(size(ia_x11_elem)))
      ia_lap_matmul_elem = 0
      
      allocate(ja_lap_matmul_elem(size(ja_x11_elem)))
      ja_lap_matmul_elem = 0
      
      allocate(a_lap_matmul_elem(nequations,nequations,size(a_x11_elem)))
      a_lap_matmul_elem = 0.0_wp

      ! [D_x1] [dummy_matrix_elem][D_x1][dummy_matrix_elem]
      ! ---------------------------------------------------
      ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of the
      ! element
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
        & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
        & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & size(a_x1_elem),iw_viscous_1,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of the
      ! element
      ia_2_matmul_elem = 0
      ja_2_matmul_elem = 0
      a_2_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
        & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
        & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
        & size(a_x1_elem),iw_viscous_1,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Compute [D_x1] [dummy_matrix_elem] times [D_x1] [dummy_matrix_elem]
      ia_lap_matmul_elem = 0
      ja_lap_matmul_elem = 0
      a_lap_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
        & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,a_lap_matmul_elem, &
        & ja_lap_matmul_elem,ia_lap_matmul_elem, &
        & size(a_x11_elem),iw_viscous_1,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Construct the array for all the elements owned by the processor
      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x11_matmul(cnt+1))
      ia_x11_matmul = 0
      
      allocate(ja_x11_matmul(nnz))  
      ja_x11_matmul = 0

      allocate(ka_x11_matmul(nnz))  
      ka_x11_matmul = 0

      nnz = 1 
      cnt = 1 
      ia_x11_matmul(1) = cnt 
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x11_matmul(cnt) = ia_lap_matmul_elem(i) + &
            & (i_elem-1)*nodesperelem*nodesperedge
          do j = ia_lap_matmul_elem(i), ia_lap_matmul_elem(i+1)-1
            ja_x11_matmul(nnz) = ja_lap_matmul_elem(j) + (i_elem-1)*nodesperelem

            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x11_matmul(cnt) = ia_lap_matmul_elem(nodesperelem+1) + &
          & (i_elem-1)*nodesperelem*nodesperedge
      end do

      ! Shift to global indexing for parallel computations 
      ja_x11_matmul = ja_x11_matmul + global_shift


      ! x22
      ! --
      ! [D_x2] [dummy_matrix_elem][D_x2][dummy_matrix_elem]
      ! ---------------------------------------------------
      ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of the
      ! element
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
        & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
        & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & size(a_x2_elem),iw_viscous_1,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of the
      ! element
      ia_2_matmul_elem = 0
      ja_2_matmul_elem = 0
      a_2_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
        & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
        & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
        size(a_x2_elem),iw_viscous_1, &
        & i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Compute [D_x2] [dummy_matrix_elem] times [D_x2] [dummy_matrix_elem]
      ia_lap_matmul_elem = 0
      ja_lap_matmul_elem = 0
      a_lap_matmul_elem = 0.0_wp
      iw_viscous_1 = 0
      call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
        & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem,a_lap_matmul_elem, &
        & ja_lap_matmul_elem,ia_lap_matmul_elem, &
        & size(a_x22_elem),iw_viscous_1,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Construct the array for all the elements owned by the processor
      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x22_matmul(cnt+1))
      ia_x22_matmul = 0
      
      allocate(ja_x22_matmul(nnz))  
      ja_x22_matmul = 0

      allocate(ka_x22_matmul(nnz))  
      ka_x22_matmul = 0
      
      nnz = 1 
      cnt = 1 
      ia_x22_matmul(1) = cnt 
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x22_matmul(cnt) = ia_lap_matmul_elem(i) + &
            & (i_elem-1)*nodesperelem*nodesperedge
          do j = ia_lap_matmul_elem(i), ia_lap_matmul_elem(i+1)-1
            ja_x22_matmul(nnz) = ja_lap_matmul_elem(j) + (i_elem-1)*nodesperelem

            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x22_matmul(cnt) = ia_lap_matmul_elem(nodesperelem+1) + &
          & (i_elem-1)*nodesperelem*nodesperedge
      end do
 
      ! Shift to global indexing for parallel computations 
      ja_x22_matmul = ja_x22_matmul + global_shift


      if (crossterms) then
        ! x12
        ! --
        allocate(iw_viscous_2(size(ia_x12_elem)-1))
        iw_viscous_2 = 0
        
        allocate(ia_ct_matmul_elem(size(ia_x12_elem)))
        ia_ct_matmul_elem = 0
        
        allocate(ja_ct_matmul_elem(size(ja_x12_elem)))
        ja_ct_matmul_elem = 0
        
        allocate(a_ct_matmul_elem(nequations,nequations,size(a_x12_elem)))
        a_ct_matmul_elem = 0.0_wp

        ! [D_x1] [dummy_matrix_elem] [D_x2] [dummy_matrix_elem]
        ! -----------------------------------------------------
        ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of the
        ! element
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
          & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x1_elem),iw_viscous_1,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of the
        ! element
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
          & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x2_elem),iw_viscous_1,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Compute [D_x1] [dummy_matrix_elem] times [D_x2] [dummy_matrix_elem]
        ia_ct_matmul_elem = 0
        ja_ct_matmul_elem = 0
        a_ct_matmul_elem = 0.0_wp
        iw_viscous_2 = 0
        call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
          & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
          & size(a_x12_elem),iw_viscous_2,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)
        
        ! Construct the array for all the elements owned by the processor
        cnt = 0
        nnz = 0
        do l = 1, n_elems
          cnt = cnt + nodesperelem
          nnz = nnz + nodesperelem*nodesperedge*nodesperedge
        enddo

        allocate(ia_x12_matmul(cnt+1))
        ia_x12_matmul = 0
        
        allocate(ja_x12_matmul(nnz))
        ja_x12_matmul = 0

        allocate(ka_x12_matmul(nnz))  
        ka_x12_matmul = 0
        
        nnz = 1
        cnt = 1
        ia_x12_matmul(1) = cnt
        do i_elem = 1, n_elems
          do i = 1, nodesperelem
            ia_x12_matmul(cnt) = ia_ct_matmul_elem(i) + &
              & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
            do j = ia_ct_matmul_elem(i), ia_ct_matmul_elem(i+1)-1
              ja_x12_matmul(nnz) = ja_ct_matmul_elem(j) + &
                & (i_elem-1)*nodesperelem
              
              nnz = nnz + 1
            enddo 
            cnt = cnt + 1
          enddo 
          ia_x12_matmul(cnt) = ia_ct_matmul_elem(nodesperelem+1) + &
            & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
        end do

        ! Shift to global indexing for parallel computations 
        ja_x12_matmul = ja_x12_matmul + global_shift


        ! x21
        ! --
        ! [D_x2] [dummy_matrix_elem] [D_x1] [dummy_matrix_elem]
        ! -----------------------------------------------------
        ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of the
        ! element
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_Sc_mu_b_Bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
          & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x2_elem),iw_viscous_1,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of the
        ! element
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_Sc_mu_b_Bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
          & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x1_elem),iw_viscous_1,i_err)
        
        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Compute [D_x2] [dummy_matrix_elem] [D_x1] [dummy_matrix_elem]
        ia_ct_matmul_elem = 0
        ja_ct_matmul_elem = 0
        a_ct_matmul_elem = 0.0_wp
        iw_viscous_2 = 0
        call a_Bl_mu_b_Bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
          & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
          & size(a_x12_elem),iw_viscous_2,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Construct the array for all the elements owned by the processor
        cnt = 0
        nnz = 0
        do l = 1, n_elems
          cnt = cnt + nodesperelem
          nnz = nnz + nodesperelem*nodesperedge*nodesperedge
        enddo

        allocate(ia_x21_matmul(cnt+1))
        ia_x21_matmul = 0
        
        allocate(ja_x21_matmul(nnz))
        ja_x21_matmul = 0

        allocate(ka_x21_matmul(nnz))  
        ka_x21_matmul = 0
        
        nnz = 1
        cnt = 1
        ia_x21_matmul(1) = cnt
        do i_elem = 1, n_elems
          do i = 1, nodesperelem
            ia_x21_matmul(cnt) = ia_ct_matmul_elem(i) + &
              & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
            do j = ia_ct_matmul_elem(i), ia_ct_matmul_elem(i+1)-1
              ja_x21_matmul(nnz) = ja_ct_matmul_elem(j) + &
                & (i_elem-1)*nodesperelem
              
              nnz = nnz + 1
            enddo 
            cnt = cnt + 1
          enddo 
          ia_x21_matmul(cnt) = ia_ct_matmul_elem(nodesperelem+1) + &
            & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
        end do

        ! Shift to global indexing for parallel computations 
        ja_x21_matmul = ja_x21_matmul + global_shift


      endif ! End if of crossterms
    endif ! End if od viscous


    if (ndim .eq. 3) then
      ! x3
      ! --
      ! [D_x3] [dummy_matrix_elem]
      ! --------------------------
      ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of the
      ! element
      ia_1_matmul_elem = 0
      ja_1_matmul_elem = 0
      a_1_matmul_elem = 0.0_wp
      iw_inviscid = 0
      call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
        & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
        & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
        & size(a_x3_elem),iw_inviscid,i_err)

      ! Check for error
      call check_sparsekit_error(i_err,message)

      ! Construct the array for all the elements owned by the processor
      cnt = 0
      nnz = 0
      do l = 1, n_elems
        cnt = cnt + nodesperelem
        nnz = nnz + nodesperelem*nodesperedge
      enddo

      allocate(ia_x3_matmul(cnt+1))
      ia_x3_matmul = 0
      
      allocate(ja_x3_matmul(nnz))
      ja_x3_matmul = 0

      allocate(ka_x3_matmul(nnz))  
      ka_x3_matmul = 0

      nnz = 1
      cnt = 1
      ia_x3_matmul(1) = cnt
      do i_elem = 1, n_elems
        do i = 1, nodesperelem
          ia_x3_matmul(cnt) = ia_1_matmul_elem(i) + &
            & (i_elem-1)*nodesperelem*nodesperedge 
          do j = ia_1_matmul_elem(i), ia_1_matmul_elem(i+1) - 1
            ja_x3_matmul(nnz) = ja_1_matmul_elem(j) + (i_elem-1)*nodesperelem

            nnz = nnz + 1
          enddo 
          cnt = cnt + 1
        enddo 
        ia_x3_matmul(cnt) = ia_1_matmul_elem(nodesperelem+1) + &
          & (i_elem-1)*nodesperelem*nodesperedge
      end do

      ! Shift to global indexing for parallel computations 
      ja_x3_matmul = ja_x3_matmul + global_shift


      if (viscous) then
        ! x33
        ! --
        ! [D_x3] [dummy_matrix_elem] [D_x3][dummy_matrix_elem]
        ! ----------------------------------------------------
        ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of 
        ! the element
        ia_1_matmul_elem = 0
        ja_1_matmul_elem = 0
        a_1_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
          & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)
 
        ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of
        ! the element
        ia_2_matmul_elem = 0
        ja_2_matmul_elem = 0
        a_2_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
          & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
          & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
          & size(a_x3_elem),iw_viscous_1,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Compute [D_x3] [dummy_matrix_elem] [D_x3][dummy_matrix_elem]
        ia_lap_matmul_elem = 0
        ja_lap_matmul_elem = 0
        a_lap_matmul_elem = 0.0_wp
        iw_viscous_1 = 0
        call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1, &
          & size(ia_2_matmul_elem)-1,nequations,a_1_matmul_elem, &
          & ja_1_matmul_elem,ia_1_matmul_elem,a_2_matmul_elem, &
          & ja_2_matmul_elem,ia_2_matmul_elem,a_lap_matmul_elem, &
          & ja_lap_matmul_elem,ia_lap_matmul_elem,size(a_x33_elem), &
          & iw_viscous_1,i_err)

        ! Check for error
        call check_sparsekit_error(i_err,message)

        ! Construct the array for all the elements owned by the processor
        cnt = 0
        nnz = 0
        do l = 1, n_elems
          cnt = cnt + nodesperelem
          nnz = nnz + nodesperelem*nodesperedge
        enddo

        allocate(ia_x33_matmul(cnt+1))
        ia_x33_matmul = 0
        
        allocate(ja_x33_matmul(nnz))  
        ja_x33_matmul = 0

        allocate(ka_x33_matmul(nnz))  
        ka_x33_matmul = 0
        
        nnz = 1 
        cnt = 1 
        ia_x33_matmul(1) = cnt 
        do i_elem = 1, n_elems
          do i = 1, nodesperelem
            ia_x33_matmul(cnt) = ia_lap_matmul_elem(i) + &
              & (i_elem-1)*nodesperelem*nodesperedge
            do j = ia_lap_matmul_elem(i), ia_lap_matmul_elem(i+1)-1
              ja_x33_matmul(nnz) = ja_lap_matmul_elem(j) + (i_elem-1)*nodesperelem

              nnz = nnz + 1
            enddo 
            cnt = cnt + 1
          enddo 
          ia_x33_matmul(cnt) = ia_lap_matmul_elem(nodesperelem+1) + &
            & (i_elem-1)*nodesperelem*nodesperedge
        end do

        ! Shift to global indexing for parallel computations 
        ja_x33_matmul = ja_x33_matmul + global_shift


        if (crossterms) then
          ! x13
          ! [D_x1] [dummy_matrix_elem] [D_x3][dummy_matrix_elem]
          ! ----------------------------------------------------
          ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_1_matmul_elem = 0
          ja_1_matmul_elem = 0
          a_1_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
            & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & size(a_x1_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
      
          ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_2_matmul_elem = 0
          ja_2_matmul_elem = 0
          a_2_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
            & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & size(a_x3_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Compute [D_x1] [dummy_matrix_elem] [D_x3][dummy_matrix_elem]
          ia_ct_matmul_elem = 0
          ja_ct_matmul_elem = 0
          a_ct_matmul_elem = 0.0_wp
          iw_viscous_2 = 0
          call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
            & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
            & size(a_x13_elem),iw_viscous_2,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
 
          ! Construct the array for all the elements owned by the processor
          cnt = 0
          nnz = 0
          do l = 1, n_elems
            cnt = cnt + nodesperelem
            nnz = nnz + nodesperelem*nodesperedge*nodesperedge
          enddo

          allocate(ia_x13_matmul(cnt+1))
          ia_x13_matmul = 0
          
          allocate(ja_x13_matmul(nnz))
          ja_x13_matmul = 0

          allocate(ka_x13_matmul(nnz))  
          ka_x13_matmul = 0
          
          nnz = 1
          cnt = 1
          ia_x13_matmul(1) = cnt
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x13_matmul(cnt) = ia_ct_matmul_elem(i) + &
                & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
              do j = ia_ct_matmul_elem(i), ia_ct_matmul_elem(i+1)-1
                ja_x13_matmul(nnz) = ja_ct_matmul_elem(j) + &
                  & (i_elem-1)*nodesperelem
                
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x13_matmul(cnt) = ia_ct_matmul_elem(nodesperelem+1) + &
              & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x13_matmul = ja_x13_matmul + global_shift


          ! x23
          ! --
          ! [D_x2] [dummy_matrix_elem] [D_x3][dummy_matrix_elem]
          ! ----------------------------------------------------
          ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_1_matmul_elem = 0
          ja_1_matmul_elem = 0
          a_1_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
            & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & size(a_x2_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
      
          ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_2_matmul_elem = 0
          ja_2_matmul_elem = 0
          a_2_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
            & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & size(a_x3_elem),iw_viscous_1,i_err)
          
          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Compute [D_x2] [dummy_matrix_elem] [D_x3][dummy_matrix_elem]
          ia_ct_matmul_elem = 0
          ja_ct_matmul_elem = 0
          a_ct_matmul_elem = 0.0_wp
          iw_viscous_2 = 0
          call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
            & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
            & size(a_x23_elem),iw_viscous_2,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Construct the array for all the elements owned by the processor
          cnt = 0
          nnz = 0
          do l = 1, n_elems
            cnt = cnt + nodesperelem
            nnz = nnz + nodesperelem*nodesperedge*nodesperedge
          enddo

          allocate(ia_x23_matmul(cnt+1))
          ia_x23_matmul = 0
          
          allocate(ja_x23_matmul(nnz))
          ja_x23_matmul = 0

          allocate(ka_x23_matmul(nnz))  
          ka_x23_matmul = 0
          
          nnz = 1
          cnt = 1
          ia_x23_matmul(1) = cnt
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x23_matmul(cnt) = ia_ct_matmul_elem(i) + &
                & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
              do j = ia_ct_matmul_elem(i), ia_ct_matmul_elem(i+1)-1
                ja_x23_matmul(nnz) = ja_ct_matmul_elem(j) + &
                  & (i_elem-1)*nodesperelem
                
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x23_matmul(cnt) = ia_ct_matmul_elem(nodesperelem+1) + &
              & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x23_matmul = ja_x23_matmul + global_shift


          ! x31
          ! --
          ! [D_x3] [dummy_matrix_elem] [D_x1][dummy_matrix_elem]
          ! ----------------------------------------------------
          ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_1_matmul_elem = 0
          ja_1_matmul_elem = 0
          a_1_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
            & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & size(a_x3_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
      
          ! Compute ([a_x1_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_2_matmul_elem = 0
          ja_2_matmul_elem = 0
          a_2_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x1_elem)-1,size(ia_x1_elem)-1,nequations, &
            & a_x1_elem,ja_x1_elem,ia_x1_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & size(a_x1_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Compute [D_x3] [\hat{C}] times [D_x1] [dW/dU]
          ia_ct_matmul_elem = 0
          ja_ct_matmul_elem = 0
          a_ct_matmul_elem = 0.0_wp
          iw_viscous_2 = 0
          call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
            & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
            & size(a_x13),iw_viscous_2,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Construct the array for all the elements owned by the processor
          cnt = 0
          nnz = 0
          do l = 1, n_elems
            cnt = cnt + nodesperelem
            nnz = nnz + nodesperelem*nodesperedge*nodesperedge
          enddo

          allocate(ia_x31_matmul(cnt+1))
          ia_x31_matmul = 0
          
          allocate(ja_x31_matmul(nnz))
          ja_x31_matmul = 0

          allocate(ka_x31_matmul(nnz))  
          ka_x31_matmul = 0
          
          nnz = 1
          cnt = 1
          ia_x31_matmul(1) = cnt
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x31_matmul(cnt) = ia_ct_matmul_elem(i) + &
                & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
              do j = ia_ct_matmul_elem(i), ia_ct_matmul_elem(i+1)-1
                ja_x31_matmul(nnz) = ja_ct_matmul_elem(j) + &
                  & (i_elem-1)*nodesperelem
                
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x31_matmul(cnt) = ia_ct_matmul_elem(nodesperelem+1) + &
              & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x31_matmul = ja_x31_matmul + global_shift


          ! x32
          ! --
          ! [D_x3] [dummy_matrix_elem] [D_x2][dummy_matrix_elem]
          ! ----------------------------------------------------
          ! Compute ([a_x3_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_1_matmul_elem = 0
          ja_1_matmul_elem = 0
          a_1_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x3_elem)-1,size(ia_x3_elem)-1,nequations, &
            & a_x3_elem,ja_x3_elem,ia_x3_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & size(a_x3_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)
      
          ! Compute ([a_x2_elem] tensor product [I_5]) times a dummy matrix of 
          ! the element
          ia_2_matmul_elem = 0
          ja_2_matmul_elem = 0
          a_2_matmul_elem = 0.0_wp
          iw_viscous_1 = 0
          call a_sc_mu_b_bl(size(ia_x2_elem)-1,size(ia_x2_elem)-1,nequations, &
            & a_x2_elem,ja_x2_elem,ia_x2_elem,dummy_matrix_elem,ja_diag_elem, &
            & ia_diag_elem,a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & size(a_x2_elem),iw_viscous_1,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Compute [D_x3] [dummy_matrix_elem] [D_x2][dummy_matrix_elem]
          ia_ct_matmul_elem = 0
          ja_ct_matmul_elem = 0
          a_ct_matmul_elem = 0.0_wp
          iw_viscous_2 = 0
          call a_bl_mu_b_bl(size(ia_1_matmul_elem)-1,size(ia_2_matmul_elem)-1, &
            & nequations,a_1_matmul_elem,ja_1_matmul_elem,ia_1_matmul_elem, &
            & a_2_matmul_elem,ja_2_matmul_elem,ia_2_matmul_elem, &
            & a_ct_matmul_elem,ja_ct_matmul_elem,ia_ct_matmul_elem, &
            & size(a_x23_elem),iw_viscous_2,i_err)

          ! Check for error
          call check_sparsekit_error(i_err,message)

          ! Construct the array for all the elements owned by the processor
          cnt = 0
          nnz = 0
          do l = 1, n_elems
            cnt = cnt + nodesperelem
            nnz = nnz + nodesperelem*nodesperedge*nodesperedge
          enddo

          allocate(ia_x32_matmul(cnt+1))
          ia_x32_matmul = 0
          
          allocate(ja_x32_matmul(nnz))
          ja_x32_matmul = 0

          allocate(ka_x32_matmul(nnz))  
          ka_x32_matmul = 0
          
          nnz = 1
          cnt = 1
          ia_x32_matmul(1) = cnt
          do i_elem = 1, n_elems
            do i = 1, nodesperelem
              ia_x32_matmul(cnt) = ia_ct_matmul_elem(i) + &
                & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
              do j = ia_ct_matmul_elem(i), ia_ct_matmul_elem(i+1)-1
                ja_x32_matmul(nnz) = ja_ct_matmul_elem(j) + &
                  & (i_elem-1)*nodesperelem
                
                nnz = nnz + 1
              enddo 
              cnt = cnt + 1
            enddo 
            ia_x32_matmul(cnt) = ia_ct_matmul_elem(nodesperelem+1) + &
              & (i_elem-1)*nodesperelem*nodesperedge*nodesperedge
          end do

          ! Shift to global indexing for parallel computations 
          ja_x32_matmul = ja_x32_matmul + global_shift


        endif ! End if crossterms

      endif ! End if viscous
    
    endif ! End if 3D

    ! Deallocate temporary memory
    if(allocated(dummy_matrix_elem)) deallocate(dummy_matrix_elem)
    if(allocated(eye_matrix)) deallocate(eye_matrix)

    ! Assign values to column pointers ka_j_matmul that point to the jaS vector
    ! -------------------------------------------------------------------------
    ! Initialize the shift
    shift = 0
    i_node = 0

    do i_elem = 1, n_elems
      do ll = 1, nodesperelem

        i_node = i_node + 1

        ! Shift for the ka_j pointers because the module
        ! jacobian_matrix_implicit_ts constructs the Jacobian matrix for a
        ! single element or cell. This shift makes the ka_j to start from 1.
        shift = iaS((i_elem-1)*nodesperelem + 1) - 1  

        totalja:do k = iaS(i_node),iaS(i_node+1) - 1
          do j = ia_0(i_node),ia_0(i_node+1) - 1
            if(ja_0_matmul(j) == jaS(k)) ka_0_matmul(j) = k - shift
          enddo

          if(ndim >= 1) then
            do j = ia_x1_matmul(i_node),ia_x1_matmul(i_node+1) - 1
              if(ja_x1_matmul(j) == jaS(k)) ka_x1_matmul(j) = k - shift
            enddo
            if(ndim >= 2) then
              do j = ia_x2_matmul(i_node),ia_x2_matmul(i_node+1) - 1
                if(ja_x2_matmul(j) == jaS(k)) ka_x2_matmul(j) = k - shift
              enddo
              if(ndim == 3) then
                do j = ia_x3_matmul(i_node),ia_x3_matmul(i_node+1) - 1
                  if(ja_x3_matmul(j) == jaS(k)) ka_x3_matmul(j) = k - shift
                enddo
              endif
            endif
          endif

          if(viscous) then
            if(ndim >= 1) then
              do j = ia_x11_matmul(i_node),ia_x11_matmul(i_node+1) - 1
                if(ja_x11_matmul(j) == jaS(k)) ka_x11_matmul(j) = k - shift
              enddo
              if(ndim >= 2) then
                do j = ia_x22_matmul(i_node),ia_x22_matmul(i_node+1) - 1
                  if(ja_x22_matmul(j) == jaS(k)) ka_x22_matmul(j) = k - shift
                enddo
                if(ndim == 3) then
                  do j = ia_x33_matmul(i_node),ia_x33_matmul(i_node+1) - 1
                    if(ja_x33_matmul(j) == jaS(k)) ka_x33_matmul(j) = k - shift
                  enddo
                endif
              endif
            endif
            if(crossterms) then
              if(ndim >= 2) then
                do j = ia_x12_matmul(i_node),ia_x12_matmul(i_node+1) - 1
                  if(ja_x12_matmul(j) == jaS(k)) ka_x12_matmul(j) = k - shift
                enddo
                do j = ia_x21_matmul(i_node),ia_x21_matmul(i_node+1) - 1
                  if(ja_x21_matmul(j) == jaS(k)) ka_x21_matmul(j) = k - shift
                enddo
              endif
              if(ndim == 3) then
                do j = ia_x13_matmul(i_node),ia_x13_matmul(i_node+1) - 1
                  if(ja_x13_matmul(j) == jaS(k)) ka_x13_matmul(j) = k - shift
                enddo
                do j = ia_x23_matmul(i_node),ia_x23_matmul(i_node+1) - 1
                  if(ja_x23_matmul(j) == jaS(k)) ka_x23_matmul(j) = k - shift
                enddo
                do j = ia_x31_matmul(i_node),ia_x31_matmul(i_node+1) - 1
                  if(ja_x31_matmul(j) == jaS(k)) ka_x31_matmul(j) = k - shift
                enddo
                do j = ia_x32_matmul(i_node),ia_x32_matmul(i_node+1) - 1
                  if(ja_x32_matmul(j) == jaS(k)) ka_x32_matmul(j) = k - shift
                enddo
              endif
            endif
          endif
        enddo totalja
      enddo

    enddo


    ! Allocate memory for the matrices stored in CSR format needed to compute 
    ! and assemble the residual Jacobian matrix of a single element (or cell).
    ! -------------------------------------------------------------------------
    ! Residual Jacobian of the element
    if (IMEX_penalty == 'implicit') then
      if (.not. viscous) then
        allocate(dfdu_a_elem(nequations,nequations,size(ja_elem) + &
          & nfacesperelem*nodesperface*2))
        dfdu_a_elem = 0.0_wp
      else
        allocate(dfdu_a_elem(nequations,nequations,size(ja_elem) + &
          & nfacesperelem*nodesperface*(1 + ndim*(nodesperedge-1))*2))
        dfdu_a_elem = 0.0_wp
      endif
    else
      allocate(dfdu_a_elem(nequations,nequations,size(ja_elem)))
      dfdu_a_elem = 0.0_wp
    endif

    ! Inviscid flux Jacobian of the element
    allocate(inviscid_flux_jacobian_elem(nequations,nequations,nodesperelem))
    inviscid_flux_jacobian_elem = 0.0_wp
    
    ! dwdu of the element
    allocate(dwdu_elem(nequations,nequations,nodesperelem))
    dwdu_elem = 0.0_wp

    ! hatc of the element
    allocate(hatc_elem(nequations,nequations,nodesperelem))
    hatc_elem = 0.0_wp
    
    ! dhatcdu*gradwj of the element
    allocate(dhatcdu_gradwj_elem(nequations,nequations,nodesperelem))
    dhatcdu_gradwj_elem = 0.0_wp

    ! Shift jaS for C indexing
    jaS = jaS - 1 
 
    return
  end subroutine csr_on_element_operator

  subroutine csr_on_element_Matrix_Multiply()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use nsereferencevariables
    use unary_mod
    use jacobian_matrix_implicit_ts_variables
    use controlvariables, only : imex_element, imex_penalty

    ! Nothing is implicitly defined
    implicit none
    
    integer :: npe, nq, n1, n12, n123, nFin
    integer, parameter :: bigN = 50000

    continue

    ! Preliminary check on the dimensions in each direction
    if (ndim .eq. 2) then
      if (size(ia_x1_elem) .ne. size(ia_x2_elem)) then
        write(*,*) 'This routine is written assuming that the order of the &
            & polynomial reconstruction in the xi and eta directions is the &
            & same and all the cell use such a reconstruction! If this is not &
            & the case then a huge portion of the code needs to be re-written.'
        write(*,*) 'Stopping...'
        stop
      endif
    else
      if (size(ia_x1_elem) .ne. size(ia_x2_elem) .or. size(ia_x1_elem) .ne. &
        & size(ia_x3_elem)) then
        write(*,*) 'This routine is written assuming that the order of the &
            & polynomial reconstruction in the xi, eta and zeta directions is &
            & the same and all the cell use such a reconstruction! If this is &
            & not the case then a huge portion of the code needs to be &
            & re-written.'
        write(*,*) 'Stopping...'
        stop
      endif
    endif

    npe  = nodesperelem
    nq   = nequations
    n1   = nodesperelem * nodesperedge
    n12  = n1 + n1 - nodesperelem
    n123 = n1 + n1 + n1 - nodesperelem - nodesperelem
    nFin = nodesperelem * nodesperedge * ((ndim-1)*nodesperedge - 1)

    allocate(ia_diag_tmp(npe+1))
    allocate(ja_diag_tmp(npe))
    allocate( a_diag_tmp(nq,nq,npe))

    allocate(ia_W1_elem(npe+1))
    allocate(ja_W1_elem(n1))
    allocate( a_W1_elem(nq,nq,n1))

    allocate(ia_W2_elem(npe+1))
    allocate(ja_W2_elem(n1))
    allocate( a_W2_elem(nq,nq,n1))

    allocate(ia_W3_elem(npe+1))
    allocate(ja_W3_elem(n1))
    allocate( a_W3_elem(nq,nq,n1))

    allocate(ia_W12_elem(npe+1))
    allocate(ja_W12_elem(n12))
    allocate( a_W12_elem(nq,nq,n12))

    allocate(ia_W123_elem(npe+1))
    allocate(ja_W123_elem(n123))
    allocate( a_W123_elem(nq,nq,n123))

    allocate(ia_dFvdU1_elem(npe+1))
    allocate(ja_dFvdU1_elem(npe))
    allocate( a_dFvdU1_elem(nq,nq,npe))

    allocate(ia_dFvdU2_elem(npe+1))
    allocate(ja_dFvdU2_elem(n123))
    allocate( a_dFvdU2_elem(nq,nq,n123))

    allocate(ia_dFvdU_elem(npe+1,3))
    allocate(ja_dFvdU_elem(n123,3))
    allocate( a_dFvdU_elem(nq,nq,n123,3))

    allocate(ia_dFvdUx_elem(npe+1))
    allocate(ja_dFvdUx_elem(nFin))
    allocate( a_dFvdUx_elem(nq,nq,nFin))

    allocate(ia_dFvdUy_elem(npe+1))
    allocate(ja_dFvdUy_elem(nFin))
    allocate( a_dFvdUy_elem(nq,nq,nFin))

    allocate(ia_dFvdUz_elem(npe+1))
    allocate(ja_dFvdUz_elem(nFin))
    allocate( a_dFvdUz_elem(nq,nq,nFin))

    allocate(ia_containerxy(npe+1))
    allocate(ja_containerxy(bigN))
    allocate( a_containerxy(nq,nq,bigN))

    allocate(ia_containerxyz(npe+1))
    allocate(ja_containerxyz(bigN))
    allocate( a_containerxyz(nq,nq,bigN))


   end subroutine csr_on_element_Matrix_Multiply


  !============================================================================
 
  !============================================================================
  ! eye - Assembles the identity matrix of dimension specified by the input
  !       parameter
  !
  ! Input parameter:
  ! n - size of the matrix

  pure function eye(n)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) ::n
    real, dimension(n,n) :: eye
    integer :: i
    
    !  Initialize the identity matrix to zero
    eye = 0.0_wp

    ! Set diagonal term
    do i = 1, n
      eye(i,i) = 1.0_wp
    enddo
    
    return
  end function

  !============================================================================

  subroutine quad_face_touch(Nvol,pp1,ndim,cnt,array)

  integer,               intent(in)  :: Nvol,pp1,ndim
  integer,               intent(out) :: cnt
  integer, dimension(:), intent(out) :: array

  integer      :: i,j,k
  integer      :: n1,n2,n3
  integer      :: stride1, stride2, remain


  if(ndim == 3) then
    stride2 = pp1 * pp1
    stride1 = pp1
    k = Nvol / stride2 + 1
    remain = mod(Nvol,stride2)
    j =         remain / stride1 + 1
    i =      mod(remain,stride1)
    if( (1 < i .and. i < pp1) .and.    &
        (1 < i .and. j < pp1) .and.    &
        (1 < i .and. k < pp1)) then
    write(*,*)'somethings up in quad_face_touch'
    write(*,*)'stopping'
    stop
  endif

  cnt = 0
  do n3 = 1,pp1
    do n2 = 1,pp1
      do n1 = 1,pp1
        if( (n3 == k) .or. (n2 == j) .or. (n1 == i)) then
          if( (n3 == k) .and. (n2 == j) .and. (n1 == i)) cycle
          cnt = cnt + 1
          array(cnt) = (k-1)*stride2 + (j-1)*stride1 + k
        endif
      enddo
    enddo
  enddo

  endif

  end subroutine quad_face_touch

  !============================================================================

  !============================================================================
  ! csr_testing - tests the CSR format of the differentiation matrices for all 
  !               elements.

  subroutine csr_testing()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use variables 
    use nsereferencevariables
    use matvec_module

    ! Nothing is implicitly defined
    integer :: cnt, n_tot_nodes, i_dir
    integer :: low_elem, high_elem, i_elem, j_node

    real(wp), allocatable, dimension(:) :: fn, exact
    real(wp), allocatable, dimension(:) :: dxi, deta, dzeta
    real(wp), allocatable, dimension(:) :: dxixi, dxieta ,dxizeta
    real(wp), allocatable, dimension(:) :: detaeta, detazeta
    real(wp), allocatable, dimension(:) :: dzetazeta
    real(wp) :: error

    continue

    ! Low volumetric element index
    low_elem = ihelems(1)

    ! High volumetric element index
    high_elem = ihelems(2)

    ! Total number of nodes owned by a processor
    n_tot_nodes = 0
    do i_elem = low_elem, high_elem
      n_tot_nodes = n_tot_nodes + nodesperelem
    enddo

    allocate(fn(n_tot_nodes))
    fn = 0.0_wp
    
    allocate(exact(n_tot_nodes))
    exact = 0.0_wp

    allocate(dxi(n_tot_nodes)) 
    dxi = 0.0_wp
    
    allocate(deta(n_tot_nodes))
    deta = 0.0_wp
    
    allocate(dzeta(n_tot_nodes))
    dzeta = 0.0_wp

    allocate(dxixi(n_tot_nodes))
    dxixi = 0.0_wp
    
    allocate(dxieta(n_tot_nodes))
    dxieta = 0.0_wp

    allocate(dxizeta(n_tot_nodes))
    dxizeta = 0.0_wp
    
    allocate(detaeta(n_tot_nodes))
    detaeta = 0.0_wp
    
    allocate(detazeta(n_tot_nodes))
    detazeta = 0.0_wp

    allocate(dzetazeta(n_tot_nodes))
    dzetazeta = 0.0_wp

    ! Initialize error to zero
    error = 0.0_wp

    !============   Testing First-derivative terms   =======
    !============ _x
    cnt = 0
    do i_elem = low_elem, high_elem
      do j_node = 1, nodesperelem 
        cnt = cnt + 1
        fn(cnt) = 1.0_wp*(xg(1,j_node,i_elem))**4
        exact(cnt) = 4.0_wp*(xg(1,j_node,i_elem))**3
      enddo
    enddo

    call amux(n_tot_nodes,fn,dxi,a_x1,ja_x1,ia_x1)
    call amux(n_tot_nodes,fn,deta,a_x2,ja_x2,ia_x2)

    if (ndim .eq. 3) then
      call amux(n_tot_nodes,fn,dzeta,a_x3,ja_x3,ia_x3)
    endif

    cnt = 0
    do i_elem = low_elem, high_elem
      do j_node = 1, nodesperelem 
        cnt = cnt + 1
        error = error + exact(cnt) - ( + dxi  (cnt)*r_x(1,1,j_node,i_elem)  &
                      + deta (cnt)*r_x(2,1,j_node,i_elem)  &
                      + dzeta(cnt)*r_x(3,1,j_node,i_elem)  )
      enddo
    enddo

    write(*,*)
    write(*,*) 'Processor ID:', myprocid
    write(*,*) 'Error D_x: ', error
    write(*,*) 
    
    !============ _y
    cnt = 0
    do i_elem = low_elem, high_elem
      do j_node = 1, nodesperelem 
        cnt = cnt + 1
        fn(cnt) = 1.0_wp*(xg(2,j_node,i_elem))**4
        exact(cnt) = 4.0_wp*(xg(2,j_node,i_elem))**3
      enddo
    enddo

    call amux(n_tot_nodes,fn,dxi,a_x1,ja_x1,ia_x1)
    call amux(n_tot_nodes,fn,deta,a_x2,ja_x2,ia_x2)

    if (ndim .eq. 3) then
      call amux(n_tot_nodes,fn,dzeta,a_x3,ja_x3,ia_x3)
    endif

    cnt = 0
    do i_elem = low_elem, high_elem
      do j_node = 1, nodesperelem 
        cnt = cnt + 1
        error = error + exact(cnt) - ( + dxi  (cnt)*r_x(1,2,j_node,i_elem)  &
                      + deta (cnt)*r_x(2,2,j_node,i_elem)  &
                      + dzeta(cnt)*r_x(3,2,j_node,i_elem)  )
      enddo
    enddo

    write(*,*)
    write(*,*) 'Processor ID:', myprocid
    write(*,*) 'Error D_y: ', error
    write(*,*) 

    !============ _z
    cnt = 0
    do i_elem = low_elem, high_elem
      do j_node = 1, nodesperelem 
        cnt = cnt + 1
        fn(cnt) = 1.0_wp*(xg(3,j_node,i_elem))**4
        exact(cnt) = 4.0_wp*(xg(3,j_node,i_elem))**3
      enddo
    enddo

    call amux(n_tot_nodes,fn,dxi,a_x1,ja_x1,ia_x1)
    call amux(n_tot_nodes,fn,deta,a_x2,ja_x2,ia_x2)

    if (ndim .eq. 3) then
      call amux(n_tot_nodes,fn,dzeta,a_x3,ja_x3,ia_x3)
    endif


    cnt = 0
    do i_elem = low_elem, high_elem
      do j_node = 1, nodesperelem 
        cnt = cnt + 1
        error = error + exact(cnt) - ( + dxi  (cnt)*r_x(1,3,j_node,i_elem)  &
                      + deta (cnt)*r_x(2,3,j_node,i_elem)  &
                      + dzeta(cnt)*r_x(3,3,j_node,i_elem)  )
      enddo
    enddo

    write(*,*)
    write(*,*) 'Processor ID:', myprocid
    write(*,*) 'Error D_z: ', error
    write(*,*) 
    
    !============

    if(viscous) then
    !============   Testing Second-derivative terms   =======
    !============ _xx, _yy, _zz
      do i_dir = 1,3
        cnt = 0
        do i_elem = low_elem, high_elem
          do j_node = 1, nodesperelem 
            cnt = cnt + 1
            fn(cnt) = 1.0_wp*(xg(i_dir,j_node,i_elem))**4
            exact(cnt) =12.0_wp*(xg(i_dir,j_node,i_elem))**2
          enddo
        enddo

        call amux(n_tot_nodes,fn,dxixi,a_x11,ja_x11,ia_x11)
        call amux(n_tot_nodes,fn,detaeta,a_x22,ja_x22,ia_x22)

        if (ndim .eq. 3) then
          call amux(n_tot_nodes,fn,dzetazeta,a_x33,ja_x33,ia_x33)
        endif


        cnt = 0
        do i_elem = low_elem, high_elem
          do j_node = 1, nodesperelem 
            cnt = cnt + 1
            error = error + exact(cnt) - ( + dxixi    (cnt)*r_x(1,i_dir,j_node,i_elem)*r_x(1,i_dir,j_node,i_elem)  &
                          + detaeta  (cnt)*r_x(2,i_dir,j_node,i_elem)*r_x(2,i_dir,j_node,i_elem)  &
                          + dzetazeta(cnt)*r_x(3,i_dir,j_node,i_elem)*r_x(3,i_dir,j_node,i_elem)  )
          enddo
        enddo
      enddo
    endif

    write(*,*)
    write(*,*) 'Processor ID:', myprocid
    write(*,*) 'Error Laplacian terms D_xx, D_yy, D_zz: ', error
    write(*,*) 

    !============   Testing Second-derivative cross-terms   =======
    !============ _xx, _xy, _yy, _xz, _yz, _zz

    if(viscous .and. crossterms) then

      ! U_{xy}
      cnt = 0
      do i_elem = low_elem, high_elem
        do j_node = 1, nodesperelem 
          cnt = cnt + 1
          fn(cnt) =  1.0_wp*xg(1,j_node,i_elem)**2 * xg(2,j_node,i_elem)**2 * xg(3,j_node,i_elem)**2
          exact(cnt) =  4.0_wp*xg(1,j_node,i_elem)**1 * xg(2,j_node,i_elem)**1 * xg(3,j_node,i_elem)**2
        enddo
      enddo

      call amux(n_tot_nodes,fn,dxixi,a_x11,ja_x11,ia_x11)
      call amux(n_tot_nodes,fn,dxieta,a_x12,ja_x12,ia_x12)
      call amux(n_tot_nodes,fn,detaeta,a_x22,ja_x22,ia_x22)

      if (ndim .eq. 3) then
        call amux(n_tot_nodes,fn,dxizeta,a_x13,ja_x13,ia_x13)
        call amux(n_tot_nodes,fn,detazeta,a_x23,ja_x23,ia_x23)
        call amux(n_tot_nodes,fn,dzetazeta,a_x33,ja_x33,ia_x33)
      endif


      cnt = 0
      do i_elem = low_elem, high_elem
        do j_node = 1, nodesperelem 
          cnt = cnt + 1
          error = error + exact(cnt) -                             &
                      ( + dxixi    (cnt)*r_x(1,1,j_node,i_elem)*r_x(1,2,j_node,i_elem)  &
                        + dxieta   (cnt)*r_x(2,1,j_node,i_elem)*r_x(1,2,j_node,i_elem)  &
                        + dxizeta  (cnt)*r_x(3,1,j_node,i_elem)*r_x(1,2,j_node,i_elem)  &

                        + dxieta   (cnt)*r_x(1,1,j_node,i_elem)*r_x(2,2,j_node,i_elem)  &
                        + detaeta  (cnt)*r_x(2,1,j_node,i_elem)*r_x(2,2,j_node,i_elem)  &
                        + detazeta (cnt)*r_x(3,1,j_node,i_elem)*r_x(2,2,j_node,i_elem)  &

                        + dxizeta  (cnt)*r_x(1,1,j_node,i_elem)*r_x(3,2,j_node,i_elem)  &
                        + detazeta (cnt)*r_x(2,1,j_node,i_elem)*r_x(3,2,j_node,i_elem)  &
                        + dzetazeta(cnt)*r_x(3,1,j_node,i_elem)*r_x(3,2,j_node,i_elem)  )
        enddo
      enddo


      write(*,*)
      write(*,*) 'Processor ID:', myprocid
      write(*,*) 'Error  D_xy: ', error
      write(*,*) 


      ! U_{xz}
      cnt = 0
      do i_elem = low_elem, high_elem
        do j_node = 1, nodesperelem 
          cnt = cnt + 1
          fn(cnt) =  1.0_wp*xg(1,j_node,i_elem)**2 * xg(2,j_node,i_elem)**2 * xg(3,j_node,i_elem)**2
          exact(cnt) =  4.0_wp*xg(1,j_node,i_elem)**1 * xg(2,j_node,i_elem)**2 * xg(3,j_node,i_elem)**1
        enddo
      enddo

      call amux(n_tot_nodes,fn,dxixi,a_x11,ja_x11,ia_x11)
      call amux(n_tot_nodes,fn,dxieta,a_x12,ja_x12,ia_x12)
      call amux(n_tot_nodes,fn,detaeta,a_x22,ja_x22,ia_x22)

      if (ndim .eq. 3) then
        call amux(n_tot_nodes,fn,dxizeta,a_x13,ja_x13,ia_x13)
        call amux(n_tot_nodes,fn,detazeta,a_x23,ja_x23,ia_x23)
        call amux(n_tot_nodes,fn,dzetazeta,a_x33,ja_x33,ia_x33)
      endif


      cnt = 0
      do i_elem = low_elem, high_elem
        do j_node = 1, nodesperelem 
          cnt = cnt + 1
          error = error + exact(cnt) -                             &
                      ( + dxixi    (cnt)*r_x(1,1,j_node,i_elem)*r_x(1,3,j_node,i_elem)  &
                        + dxieta   (cnt)*r_x(2,1,j_node,i_elem)*r_x(1,3,j_node,i_elem)  &
                        + dxizeta  (cnt)*r_x(3,1,j_node,i_elem)*r_x(1,3,j_node,i_elem)  &

                        + dxieta   (cnt)*r_x(1,1,j_node,i_elem)*r_x(2,3,j_node,i_elem)  &
                        + detaeta  (cnt)*r_x(2,1,j_node,i_elem)*r_x(2,3,j_node,i_elem)  &
                        + detazeta (cnt)*r_x(3,1,j_node,i_elem)*r_x(2,3,j_node,i_elem)  &

                        + dxizeta  (cnt)*r_x(1,1,j_node,i_elem)*r_x(3,3,j_node,i_elem)  &
                        + detazeta (cnt)*r_x(2,1,j_node,i_elem)*r_x(3,3,j_node,i_elem)  &
                        + dzetazeta(cnt)*r_x(3,1,j_node,i_elem)*r_x(3,3,j_node,i_elem)  )
        enddo
      enddo

      write(*,*)
      write(*,*) 'Processor ID:', myprocid
      write(*,*) 'Error  D_xz: ', error
      write(*,*) 

      ! U_{yz}
      cnt = 0
      do i_elem = low_elem, high_elem
        do j_node = 1, nodesperelem 
          cnt = cnt + 1
          fn(cnt) =  1.0_wp*xg(1,j_node,i_elem)**2 * xg(2,j_node,i_elem)**2 * xg(3,j_node,i_elem)**2
          exact(cnt) =  4.0_wp*xg(1,j_node,i_elem)**2 * xg(2,j_node,i_elem)**1 * xg(3,j_node,i_elem)**1
        enddo
      enddo

      call amux(n_tot_nodes,fn,dxixi,a_x11,ja_x11,ia_x11)
      call amux(n_tot_nodes,fn,dxieta,a_x12,ja_x12,ia_x12)
      call amux(n_tot_nodes,fn,detaeta,a_x22,ja_x22,ia_x22)

      if (ndim .eq. 3) then
        call amux(n_tot_nodes,fn,dxizeta,a_x13,ja_x13,ia_x13)
        call amux(n_tot_nodes,fn,detazeta,a_x23,ja_x23,ia_x23)
        call amux(n_tot_nodes,fn,dzetazeta,a_x33,ja_x33,ia_x33)
      endif



      cnt = 0
      do i_elem = low_elem, high_elem
        do j_node = 1, nodesperelem 
          cnt = cnt + 1
          error = error + exact(cnt) -                             &
                      ( + dxixi    (cnt)*r_x(1,2,j_node,i_elem)*r_x(1,3,j_node,i_elem)  &
                        + dxieta   (cnt)*r_x(2,2,j_node,i_elem)*r_x(1,3,j_node,i_elem)  &
                        + dxizeta  (cnt)*r_x(3,2,j_node,i_elem)*r_x(1,3,j_node,i_elem)  &

                        + dxieta   (cnt)*r_x(1,2,j_node,i_elem)*r_x(2,3,j_node,i_elem)  &
                        + detaeta  (cnt)*r_x(2,2,j_node,i_elem)*r_x(2,3,j_node,i_elem)  &
                        + detazeta (cnt)*r_x(3,2,j_node,i_elem)*r_x(2,3,j_node,i_elem)  &

                        + dxizeta  (cnt)*r_x(1,2,j_node,i_elem)*r_x(3,3,j_node,i_elem)  &
                        + detazeta (cnt)*r_x(2,2,j_node,i_elem)*r_x(3,3,j_node,i_elem)  &
                        + dzetazeta(cnt)*r_x(3,2,j_node,i_elem)*r_x(3,3,j_node,i_elem)  )
        enddo
      enddo

      write(*,*)
      write(*,*) 'Processor ID:', myprocid
      write(*,*) 'Error  D_yz: ', error
      write(*,*) 


    endif

    write(*,*) 
    write(*,*) 'Processor ID:', myprocid
    write(*,*) 'Total numerical error of the differentiation matrices:', error
    write(*,*)

    return
  end subroutine csr_testing

  !============================================================================

  !============================================================================
  ! csr_initialize_jacobian - Counts the number of rows of the Jacobian matrix 
  !                           owned by a process and its number of nonzero 
  !                           elements. 

  subroutine csr_initialize_jacobian()

    ! Load modules
    use CSRlocalvariables
    use referencevariables
    use variables 
    use mpimod 
    use nsereferencevariables

    integer :: nnz
    integer :: low_elem, high_elem, i_elem, i_face, i_node
    integer :: i_err

    continue

    ! Low volumetric element index
    low_elem = ihelems(1) 

    ! High volumetric element index
    high_elem = ihelems(2)

    ! Set to zero number of nonzero (nnz) elements
    nnz = 0

    ! Set to zero number of rows owned by a process and global number of rows    
    nprows = 0
    ngrows = 0

    ! Loop over elements
    do i_elem = low_elem, high_elem
      ! Interior contributions
      ! ----------------------

      ! Inviscid nonzero terms
      if(.not. viscous) then
        ! Each node will have 1 nnz element coming from itself plus 
        ! nodesperedge-1 contributions from the other nodes which lie in each 
        ! direction (tensor product cell), i.e nodesperedge-1. Hence, if the
        ! the problem is ndim-dimensional there are (nodesperedge-1)*ndim nnz 
        !contribution for one node
        nnz = nnz + (1+(nodesperedge-1)*ndim)*nodesperedge**ndim
      else 
        ! Viscous nonzero terms without cross derivatives, i.e. Laplacian 
        ! contributions. Note that the footprint of the Laplacian contribution 
        ! to the Jacobian matrix is exactly the same of that of the inviscid 
        ! contribution. Thus, the above discussion applies also here. 
        if(.not. crossterms) then
          nnz = nnz + (1 + (nodesperedge-1)*ndim) * nodesperedge**ndim
        else
          ! Cross contribution
          if (ndim == 2) then
            nnz = nnz + ( + (nodesperedge-0)*(nodesperedge-0))*nodesperedge**ndim
          elseif (ndim == 3) then
            nnz = nnz + ( + (nodesperedge-0)*(nodesperedge-0) &
                          + (nodesperedge-0)*(nodesperedge-1) &
                          + (nodesperedge-1)*(nodesperedge-1))*nodesperedge**ndim
          endif
        endif
      endif

      ! loop over faces
      do i_face = 1,nfacesperelem
        ! if on boundary, connect to self
        if (ef2e(1,i_face,i_elem) < 0) cycle
        ! loop over nodes on face
        do i_node = 1,nodesperface
          if(.not. viscous) then
            nnz = nnz + 1
          else
            ! assume no cross-derivative contributions in derivative penalty
            nnz = nnz + 1 + (nodesperedge-1)*ndim
          endif
        enddo
      enddo

      ! Count number of rows of the Jacobian matrix owned by the process
      nprows = nprows + nodesperelem*nequations

    enddo

    ! Reduce values on all processes to a single value
    call MPI_reduce(nprows,ngrows,1,mpi_integer,mpi_sum ,0,petsc_comm_world, &
      & i_err)

    ! Broadcast a message from the process with rank "root" to all other 
    ! processes of the communica
    call mpi_bcast(ngrows,1,mpi_integer,0,petsc_comm_world,i_err)

  end subroutine csr_initialize_jacobian

  !============================================================================

  !============================================================================
  ! check_sparsekit_error - Check error sent back by the sparsekit
  !
  ! Input parameters:
  ! i_err  - integer, if larger than 0 than there was an error

  subroutine check_sparsekit_error(i_err,message)

    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: i_err
    character(120), intent(in) :: message

    ! Check error status
    if(i_err > 0) then
      write(*,*) 'Failure in initialize_csr.f90; i_err > 0 in ' // trim(message)
      write(*,*) 'Stopping'
      stop
    endif

    return
  endsubroutine check_sparsekit_error

  !============================================================================

end module initialize_CSR

