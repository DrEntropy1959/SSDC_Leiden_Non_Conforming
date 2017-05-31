module matvec_module

  use precision_vars

  implicit none

  !     private
  !     public ::  amux , atmux , atmuxr , lsol , usol 

contains

  !=============================================================================80


  !----------------------------------------------------------------------c
  !                          S P A R S K I T                             c
  !----------------------------------------------------------------------c
  !          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
  !         Matrix-vector Mulitiplications and Triang. Solves            c
  !----------------------------------------------------------------------c
  ! contents: (as of Nov 18, 1991)                                       c
  !----------                                                            c
  ! 1) Matrix-vector products:                                           c
  !---------------------------                                           c
  ! amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
  ! atmux : Transp(A) times a vector. CSR format.                        c
  ! atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
  !                                                                      c
  ! 2) Triangular system solutions:                                      c
  !-------------------------------                                       c
  ! lsol  : Unit Lower Triang. solve. Compressed Sparse Row (CSR) format.c
  ! usol  : Unit Upper Triang. solve. Compressed Sparse Row (CSR) format.c
  ! udsolc: Upper Triang. solve.  Modified Sparse Column (MSC) format.   c
  !----------------------------------------------------------------------c
  ! 1)     M A T R I X    B Y    V E C T O R     P R O D U C T S         c
  !----------------------------------------------------------------------c
  subroutine amux (n, x, y, a,ja,ia) 
    real(wp)  x(*), y(*), a(*) 
    integer n, ja(*), ia(*)

    !-----------------------------------------------------------------------
    !         A times a vector
    !----------------------------------------------------------------------- 
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage.
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !
    !-----------------------------------------------------------------------
    ! local variables

    real(wp) t
    integer i, k

    do 100 i = 1,n

      !     compute the inner product of row i with vector x

      t = 0.0_wp
      do 99 k=ia(i), ia(i+1)-1 
        t = t + a(k)*x(ja(k))
        99      continue

        !     store result in y(i) 

        y(i) = t
        100  continue

      end subroutine amux

      !=============================================================================80

      subroutine atmux (n, x, y, a, ja, ia)
        real(wp) x(*), y(*), a(*) 
        integer n, ia(*), ja(*)
        !-----------------------------------------------------------------------
        !         transp( A ) times a vector
        !----------------------------------------------------------------------- 
        ! multiplies the transpose of a matrix by a vector when the original
        ! matrix is stored in compressed sparse row storage. Can also be
        ! viewed as the product of a matrix by a vector when the original
        ! matrix is stored in the compressed sparse column format.
        !-----------------------------------------------------------------------
        !
        ! on entry:
        !----------
        ! n     = row dimension of A
        ! x     = real array of length equal to the column dimension of
        !         the A matrix.
        ! a, ja,
        !    ia = input matrix in compressed sparse row format.
        !
        ! on return:
        !-----------
        ! y     = real array of length n, containing the product y=transp(A)*x
        !
        !-----------------------------------------------------------------------
        !     local variables 
        !
        integer i, k 
        !-----------------------------------------------------------------------
        !
        !     zero out output vector
        ! 
        do 1 i=1,n
          y(i) = 0.0
          1    continue

          ! loop over the rows

          do 100 i = 1,n
            do 99 k=ia(i), ia(i+1)-1 
              y(ja(k)) = y(ja(k)) + x(i)*a(k)
              99      continue
              100  continue

            end subroutine atmux

            !=============================================================================80

            subroutine atmuxr (m, n, x, y, a, ja, ia)
              real(wp) x(*), y(*), a(*) 
              integer m, n, ia(*), ja(*)

              !-----------------------------------------------------------------------
              !         transp( A ) times a vector, A can be rectangular
              !----------------------------------------------------------------------- 
              ! See also atmux.  The essential difference is how the solution vector
              ! is initially zeroed.  If using this to multiply rectangular CSC 
              ! matrices by a vector, m number of rows, n is number of columns.
              !-----------------------------------------------------------------------
              !
              ! on entry:
              !----------
              ! m     = column dimension of A
              ! n     = row dimension of A
              ! x     = real array of length equal to the column dimension of
              !         the A matrix.
              ! a, ja,
              !    ia = input matrix in compressed sparse row format.
              !
              ! on return:
              !-----------
              ! y     = real array of length n, containing the product y=transp(A)*x
              !
              !-----------------------------------------------------------------------
              !     local variables 
              !
              integer i, k 

              !     zero out output vector

              do 1 i=1,m
                y(i) = 0.0
                1    continue

                ! loop over the rows

                do 100 i = 1,n
                  do 99 k=ia(i), ia(i+1)-1 
                    y(ja(k)) = y(ja(k)) + x(i)*a(k)
                    99      continue
                    100  continue

                  end subroutine atmuxr

                  !=============================================================================80

                  !----------------------------------------------------------------------c
                  ! 2)     T R I A N G U L A R    S Y S T E M    S O L U T I O N S       c
                  !----------------------------------------------------------------------c
                  subroutine lsol (n,x,y,al,jal,ial)
                    integer n, jal(*),ial(n+1) 
                    real(wp)  x(n), y(n), al(*) 

                    !-----------------------------------------------------------------------
                    !   solves    L x = y ; L = lower unit triang. /  CSR format
                    !----------------------------------------------------------------------- 
                    ! solves a unit lower triangular system by standard (sequential )
                    ! forward elimination - matrix stored in CSR format. 
                    !-----------------------------------------------------------------------
                    !
                    ! On entry:
                    !---------- 
                    ! n      = integer. dimension of problem.
                    ! y      = real array containg the right side.
                    !
                    ! al,
                    ! jal,
                    ! ial,    = Lower triangular matrix stored in compressed sparse row
                    !          format. 
                    !
                    ! On return:
                    !----------- 
                    !	x  = The solution of  L x  = y.
                    !--------------------------------------------------------------------
                    ! local variables 

                    integer k, j 
                    real(wp)  t

                    x(1) = y(1) 
                    do 150 k = 2, n
                      t = y(k) 
                      do 100 j = ial(k), ial(k+1)-1
                        t = t-al(j)*x(jal(j))
                        100     continue
                        x(k) = t 
                        150  continue

                      end subroutine lsol

                      !=============================================================================80

                      subroutine ldsol (n,x,y,al,jal)
                        integer n, jal(*)
                        real(wp) x(n), y(n), al(*)
                        !-----------------------------------------------------------------------
                        !     Solves L x = y    L = triangular. MSR format
                        !-----------------------------------------------------------------------
                        ! solves a (non-unit) lower triangular system by standard (sequential)
                        ! forward elimination - matrix stored in MSR format
                        ! with diagonal elements already inverted (otherwise do inversion,
                        ! al(1:n) = 1.0/al(1:n),  before calling ldsol).
                        !-----------------------------------------------------------------------
                        !
                        ! On entry:
                        !----------
                        ! n      = integer. dimension of problem.
                        ! y      = real array containg the right hand side.
                        !
                        ! al,
                        ! jal,   = Lower triangular matrix stored in Modified Sparse Row
                        !          format.
                        !
                        ! On return:
                        !-----------
                        !       x = The solution of  L x = y .
                        !--------------------------------------------------------------------
                        ! local variables
                        !
                        integer k, j
                        real(wp) t
                        !-----------------------------------------------------------------------
                        x(1) = y(1)*al(1)
                        do 150 k = 2, n
                          t = y(k)
                          do 100 j = jal(k), jal(k+1)-1
                            t = t - al(j)*x(jal(j))
                            100     continue
                            x(k) = al(k)*t
                            150  continue
                            return

                          end subroutine ldsol

                          !=============================================================================80

                          subroutine usol (n,x,y,au,jau,iau)
                            integer n, jau(*),iau(n+1) 
                            real(wp)  x(n), y(n), au(*) 
                            !----------------------------------------------------------------------- 
                            !             Solves   U x = y    U = unit upper triangular. 
                            !-----------------------------------------------------------------------
                            ! solves a unit upper triangular system by standard (sequential )
                            ! backward elimination - matrix stored in CSR format. 
                            !-----------------------------------------------------------------------
                            !
                            ! On entry:
                            !---------- 
                            ! n      = integer. dimension of problem.
                            ! y      = real array containg the right side.
                            !
                            ! au,
                            ! jau,
                            ! iau,    = Lower triangular matrix stored in compressed sparse row
                            !          format. 
                            !
                            ! On return:
                            !----------- 
                            !	x = The solution of  U x = y . 
                            !-------------------------------------------------------------------- 
                            ! local variables 
                            !
                            integer k, j 
                            real(wp)  t

                            x(n) = y(n) 
                            do 150 k = n-1,1,-1 
                              t = y(k) 
                              do 100 j = iau(k), iau(k+1)-1
                                t = t - au(j)*x(jau(j))
                                100     continue
                                x(k) = t 
                                150  continue

                              end subroutine usol

                              !=============================================================================80
                              subroutine udsol (n,x,y,au,jau)

                                integer n, jau(*)
                                real(wp)  x(n), y(n),au(*)
                                !-----------------------------------------------------------------------
                                !             Solves   U x = y  ;   U = upper triangular in MSR format
                                !-----------------------------------------------------------------------
                                ! solves a non-unit upper triangular matrix by standard (sequential )
                                ! backward elimination - matrix stored in MSR format.
                                ! with diagonal elements already inverted (otherwise do inversion,
                                ! au(1:n) = 1.0/au(1:n),  before calling).
                                !-----------------------------------------------------------------------
                                !
                                ! On entry:
                                !----------
                                ! n      = integer. dimension of problem.
                                ! y      = real array containg the right side.
                                !
                                ! au,
                                ! jau,    = Lower triangular matrix stored in modified sparse row
                                !          format.
                                !
                                ! On return:
                                !-----------
                                !       x = The solution of  U x = y .
                                !--------------------------------------------------------------------
                                ! local variables
                                !
                                integer k, j
                                real(wp) t
                                !-----------------------------------------------------------------------
                                x(n) = y(n)*au(n)
                                do 150 k = n-1,1,-1
                                  t = y(k)
                                  do 100 j = jau(k), jau(k+1)-1
                                    t = t - au(j)*x(jau(j))
                                    100     continue
                                    x(k) = au(k)*t
                                    150  continue

                                  end subroutine udsol


                                end module matvec_module
