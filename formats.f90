!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                    FORMAT CONVERSION MODULE                          c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! csrdns  : converts a row-stored sparse matrix into the dense format. c
! dnscsr  : converts a dense matrix to a sparse storage format.        c
! coocsr  : converts coordinate to  to csr format                      c
! coicsr  : in-place conversion of coordinate to csr format            c
! csrcoo  : converts compressed sparse row to coordinate.              c
! csrssr  : converts compressed sparse row to symmetric sparse row     c
! ssrcsr  : converts symmetric sparse row to compressed sparse row     c
! csrell  : converts compressed sparse row to ellpack format           c
! ellcsr  : converts ellpack format to compressed sparse row format    c
! csrmsr  : converts compressed sparse row format to modified sparse   c
!           row format                                                 c
! msrcsr  : converts modified sparse row format to compressed sparse   c
!           row format.                                                c
! csrcsc  : converts compressed sparse row format to compressed sparse c
!           column format (transposition)                              c
! csrcsc2 : rectangular version of csrcsc                              c
! csrlnk  : converts compressed sparse row to linked list format       c
! lnkcsr  : converts linked list format to compressed sparse row fmt   c
! csrdia  : converts a compressed sparse row format into a diagonal    c
!           format.                                                    c
! diacsr  : converts a diagonal format into a compressed sparse row    c
!           format.                                                    c
! bsrcsr  : converts a block-row sparse format into a compressed       c
!           sparse row format.                                         c
! csrbsr  : converts a compressed sparse row format into a block-row   c
!           sparse format.                                             c
! csrbnd  : converts a compressed sparse row format into a banded      c
!           format (linpack style).                                    c
! bndcsr  : converts a banded format (linpack style) into a compressed c
!           sparse row storage.                                        c
! csrssk  : converts the compressed sparse row format to the symmetric c
!           skyline format                                             c
! sskssr  : converts symmetric skyline format to symmetric  sparse row c
!           format.                                                    c
! csrjad  : converts the csr format into the jagged diagonal format    c
! jadcsr  : converts the jagged-diagonal format into the csr format    c
! csruss  : Compressed Sparse Row to Unsymmetric Sparse Skyline        c
!           format                                                     c
! usscsr  : Unsymmetric Sparse Skyline format to Compressed Sparse Row c
! csrsss  : Compressed Sparse Row to Symmetric Sparse Skyline format   c
! ssscsr  : Symmetric Sparse Skyline format to Compressed Sparse Row   c
! csrvbr  : Converts compressed sparse row to var block row format     c
! vbrcsr  : Converts var block row to compressed sparse row format     c
! csorted : Checks if matrix in CSR format is sorted by columns        c
!--------- miscalleneous additions not involving the csr format--------c
! cooell  : converts coordinate to Ellpack/Itpack format               c
! dcsort  : sorting routine used by crsjad                             c
!----------------------------------------------------------------------c
      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr) 
      real*8 dns(ndns,*),a(*)
      integer ja(*),ia(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row    to    Dense 
!-----------------------------------------------------------------------
!
! converts a row-stored sparse matrix into a densely stored one
!
! On entry:
!---------- 
!
! nrow	= row-dimension of a
! ncol	= column dimension of a
! a, 
! ja, 
! ia    = input matrix in compressed sparse row format. 
!         (a=value array, ja=column array, ia=pointer array)
! dns   = array where to store dense matrix
! ndns	= first dimension of array dns 
!
! on return: 
!----------- 
! dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
! 
! ierr  = integer error indicator. 
!         ierr .eq. 0  means normal return
!         ierr .eq. i  means that the code has stopped when processing
!         row number i, because it found a column number .gt. ncol.
! 
!----------------------------------------------------------------------- 
      ierr = 0
      do 1 i=1, nrow
         do 2 j=1,ncol
	    dns(i,j) = 0.0d0
 2       continue
 1    continue
!     
      do 4 i=1,nrow
         do 3 k=ia(i),ia(i+1)-1
            j = ja(k) 
	    if (j .gt. ncol) then
               ierr = i
               return
	    endif
	    dns(i,j) = a(k)
 3       continue	   
 4    continue
      return
!---- end of csrdns ----------------------------------------------------
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
      real*8 dns(ndns,*),a(*)
      integer ia(*),ja(*)
!-----------------------------------------------------------------------
! Dense		to    Compressed Row Sparse 
!----------------------------------------------------------------------- 
!
! converts a densely stored matrix into a row orientied
! compactly sparse matrix. ( reverse of csrdns )
! Note: this routine does not check whether an element 
! is small. It considers that a(i,j) is zero if it is exactly
! equal to zero: see test below.
!-----------------------------------------------------------------------
! on entry:
!---------
!
! nrow	= row-dimension of a
! ncol	= column dimension of a
! nzmax = maximum number of nonzero elements allowed. This
!         should be set to be the lengths of the arrays a and ja.
! dns   = input nrow x ncol (dense) matrix.
! ndns	= first dimension of dns. 
!
! on return:
!---------- 
! 
! a, ja, ia = value, column, pointer  arrays for output matrix 
!
! ierr	= integer error indicator: 
!         ierr .eq. 0 means normal retur
!         ierr .eq. i means that the the code stopped while
!         processing row number i, because there was no space left in
!         a, and ja (as defined by parameter nzmax).
!----------------------------------------------------------------------- 
      ierr = 0
      next = 1
      ia(1) = 1
      do 4 i=1,nrow
         do 3 j=1, ncol 
            if (dns(i,j) .eq. 0.0d0) goto 3
            if (next .gt. nzmax) then
               ierr = i
               return
            endif
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
 3       continue	   
         ia(i+1) = next
 4    continue
      return
!---- end of dnscsr ---------------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!----------------------------------------------------------------------- 
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row 
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!--------- 
! nrow	= dimension of the matrix 
! nnz	= number of nonzero elements in matrix
! a,
! ir, 
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column 
!	  number. The order of the elements is arbitrary. 
!
! on return:
!----------- 
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao 
! 	continung the real values, jao containing the column indices, 
!	and iao being the pointer to the beginning of the row, 
!	in arrays ao, jao.
!
! Notes:
!------ This routine is NOT in place.  See coicsr
!
!------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
! determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
! starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
! go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
! shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!------------- end of coocsr ------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine coicsr (n,nnz,job,a,ja,ia,iwk)
      integer ia(nnz),ja(nnz),iwk(n+1) 
      real*8 a(*)
!------------------------------------------------------------------------
! IN-PLACE coo-csr conversion routine.
!------------------------------------------------------------------------
! this subroutine converts a matrix stored in coordinate format into 
! the csr format. The conversion is done in place in that the arrays 
! a,ja,ia of the result are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!--------- 
! n	= integer. row dimension of A.
! nnz	= integer. number of nonzero elements in A.
! job   = integer. Job indicator. when job=1, the real values in a are
!         filled. Otherwise a is not touched and the structure of the
!         array only (i.e. ja, ia)  is obtained.
! a	= real array of size nnz (number of nonzero elements in A)
!         containing the nonzero elements 
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer array of length nnz containing the row positions
! 	  of the corresponding elements in a.
! iwk	= integer work array of length n+1 
! on return:
!----------
! a
! ja 
! ia	= contains the compressed sparse row data structure for the 
!         resulting matrix.
! Note: 
!-------
!         the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use coocsr
!         if you want them sorted.
!----------------------------------------------------------------------c
!  Coded by Y. Saad, Sep. 26 1989                                      c
!----------------------------------------------------------------------c
      real*8 t,tnext
      logical values
!----------------------------------------------------------------------- 
      values = (job .eq. 1) 
! find pointer array for resulting matrix. 
      do 35 i=1,n+1
         iwk(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
 4    continue 
!------------------------------------------------------------------------
      iwk(1) = 1 
      do 44 i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
 44   continue 
!
!     loop for a cycle in chasing process. 
!
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
!------------------------------------------------------------------------
 6    k = k+1 		
!     current row number is i.  determine  where to go. 
      ipos = iwk(i)
!     save the chased element. 
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
!     then occupy its location.
      if (values) a(ipos)  = t
      ja(ipos) = j
!     update pointer information for next element to come in row i. 
      iwk(i) = ipos+1
!     determine  next element to be chased,
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext 
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
!     restart chasing --	
      goto 5
 70   do 80 i=1,n 
         ia(i+1) = iwk(i)
 80   continue
      ia(1) = 1
      return
!----------------- end of coicsr ----------------------------------------
!------------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine csrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
!-----------------------------------------------------------------------
      real*8 a(*),ao(*) 
      integer ir(*),jc(*),ja(*),ia(nrow+1) 
!----------------------------------------------------------------------- 
!  Compressed Sparse Row      to      Coordinate 
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry: 
!---------
! nrow	= dimension of the matrix.
! job   = integer serving as a job indicator. 
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao 
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return 
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) 
!         will write the output matrix in coordinate format on a, ja,ia.
!
! a,
! ja,
! ia    = matrix in compressed sparse row format.
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
! 
! on return:
!----------- 
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
! ierr       = integer error indicator.
!         ierr .eq. 0 means normal retur
!         ierr .eq. 1 means that the the code stopped 
!         because there was no space in ao, ir, jc 
!         (according to the value of  nzmax).
! 
! NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with 
!         ao being the same array as as a, and jc the same array as ja. 
!         but ir CANNOT be the same as ia. 
!         2) note the order in the output arrays, 
!------------------------------------------------------------------------
      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         return
      endif
!------------------------------------------------------------------------
      goto (3,2,1) job
 1    do 10 k=1,nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,nnz
         jc(k) = ja(k)
 11   continue
!
!     copy backward to allow for in-place processing. 
!
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
!------------- end-of-csrcoo ------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine csrssr (nrow,a,ja,ia,nzmax,ao,jao,iao,ierr)
      real*8 a(*), ao(*), t
      integer ia(*), ja(*), iao(*), jao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Symmetric Sparse Row
!----------------------------------------------------------------------- 
! this subroutine extracts the lower triangular part of a matrix.
! It can used as a means for converting a symmetric matrix for 
! which all the entries are stored in sparse format into one
! in which only the lower part is stored. The routine is in place in 
! that the output matrix ao, jao, iao can be overwritten on 
! the  input matrix  a, ja, ia if desired. Csrssr has been coded to
! put the diagonal elements of the matrix in the last position in
! each row (i.e. in position  ao(ia(i+1)-1   of ao and jao) 
!----------------------------------------------------------------------- 
! On entry
!-----------
! nrow  = dimension of the matrix a.
! a, ja, 
!    ia = matrix stored in compressed row sparse format
!
! nzmax = length of arrays ao,  and jao. 
!
! On return:
!----------- 
! ao, jao, 
!     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse 
!          row format format.
!  
! ierr   = integer error indicator. 
!          ierr .eq. 0  means normal return
!          ierr .eq. i  means that the code has stopped when processing
!          row number i, because there is not enough space in ao, jao
!          (according to the value of nzmax) 
!
!----------------------------------------------------------------------- 
      ierr = 0
      ko = 0
!-----------------------------------------------------------------------
      do  7 i=1, nrow
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            if (ko .gt. nzmax) then
               ierr = i
               return
            endif
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
!     
!     exchange
!     
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t
!     
         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
!     redefine iao(n+1)
      iao(nrow+1) = ko+1
      return
!--------- end of csrssr ----------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine ssrcsr(job, value2, nrow, a, ja, ia, nzmax, ao, jao, iao, indu, iwk, ierr)
!     .. Scalar Arguments ..
      integer            ierr, job, nrow, nzmax, value2
!     ..
!     .. Array Arguments ..
      integer            ia(nrow+1), iao(nrow+1), indu(nrow), iwk(nrow+1), ja(*), jao(nzmax)
      real*8             a(*), ao(nzmax)
!     ..
!-----------------------------------------------------------------------
!     Symmetric Sparse Row to Compressed Sparse Row format
!-----------------------------------------------------------------------
!     This subroutine converts a given matrix in SSR format to regular
!     CSR format by computing Ao = A + A' - diag(A), where A' is A
!     transpose.
!
!     Typically this routine is used to expand the SSR matrix of
!     Harwell Boeing matrices, or to obtain a symmetrized graph of
!     unsymmetric matrices.
!
!     This routine is inplace, i.e., (Ao,jao,iao) may be same as
!     (a,ja,ia).
!
!     It is possible to input an arbitrary CSR matrix to this routine,
!     since there is no syntactical difference between CSR and SSR
!     format. It also removes duplicate entries and perform a partial
!     ordering. The output matrix has an order of lower half, main
!     diagonal and upper half after the partial ordering.
!-----------------------------------------------------------------------
! on entry:
!---------
!
! job   = options
!         0 -- duplicate entries are not removed. If the input matrix is
!              SSR (not an arbitary CSR) matrix, no duplicate entry should
!              arise from this routine.
!         1 -- eliminate duplicate entries, zero entries.
!         2 -- eliminate duplicate entries and perform partial ordering.
!         3 -- eliminate duplicate entries, sort the entries in the
!              increasing order of clumn indices.
!              
! value2= will the values of A be copied?
!         0 -- only expand the graph (a, ao are not touched)
!         1 -- expand the matrix with the values.
!
! nrow  = column dimension of inout matrix
! a,
! ia,
! ja    = matrix in compressed sparse row format.
!
! nzmax = size of arrays ao and jao. SSRCSR will abort if the storage
!          provided in ao, jao is not sufficient to store A. See ierr.
!
! on return:
!----------
! ao, jao, iao
!       = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A'-D. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.
!
! indu  = integer array of length nrow. INDU will contain pointers
!         to the beginning of upper traigular part if job > 1.
!         Otherwise it is also used as a work array (size nrow).
!
! iwk   = integer work space (size nrow+1).
!
! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
!
!-----------------------------------------------------------------------
!     .. Local Scalars ..
      integer            i, ipos, j, k, kfirst, klast, ko, kosav, nnz
      real*8             tmp
!     ..
!     .. Executable Statements ..
      ierr = 0
      do 10 i = 1, nrow
         indu(i) = 0
         iwk(i) = 0
 10   continue
      iwk(nrow+1) = 0
!
!     .. compute number of elements in each row of (A'-D)
!     put result in iwk(i+1)  for row i.
!
      do 30 i = 1, nrow
         do 20 k = ia(i), ia(i+1) - 1
            j = ja(k)
            if (j.ne.i) then
              iwk(j+1) = iwk(j+1) + 1
            endif
 20      continue
 30   continue
!
!     .. find addresses of first elements of ouput matrix. result in iwk
!
      iwk(1) = 1
      do 40 i = 1, nrow
         indu(i) = iwk(i) + ia(i+1) - ia(i)
         iwk(i+1) = iwk(i+1) + indu(i)
         indu(i) = indu(i) - 1
 40   continue
!.....Have we been given enough storage in ao, jao ?
      nnz = iwk(nrow+1) - 1
      if (nnz.gt.nzmax) then
         ierr = nnz
         return
      endif
!
!     .. copy the existing matrix (backwards).
!
      kosav = iwk(nrow+1)
      do 60 i = nrow, 1, -1
         klast = ia(i+1) - 1
         kfirst = ia(i)
         iao(i+1) = kosav
         kosav = iwk(i)
         ko = iwk(i) - kfirst
         iwk(i) = ko + klast + 1
         do 50 k = klast, kfirst, -1
            if (value2.ne.0) then
              ao(k+ko) = a(k)
            endif
            jao(k+ko) = ja(k)
 50      continue
 60   continue
      iao(1) = 1
!
!     now copy (A'-D). Go through the structure of ao, jao, iao
!     that has already been copied. iwk(i) is the address
!     of the next free location in row i for ao, jao.
!
      do 80 i = 1, nrow
         do 70 k = iao(i), indu(i)
            j = jao(k)
            if (j.ne.i) then
               ipos = iwk(j)
               if (value2.ne.0) then
                 ao(ipos) = ao(k)
               endif
               jao(ipos) = i
               iwk(j) = ipos + 1
            endif
 70      continue
 80   continue
      if (job.le.0) return
!
!     .. eliminate duplicate entries --
!     array INDU is used as marker for existing indices, it is also the
!     location of the entry.
!     IWK is used to stored the old IAO array.
!     matrix is copied to squeeze out the space taken by the duplicated
!     entries.
!
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = iao(i)
 90   continue
      iwk(nrow+1) = iao(nrow+1)
      k = 1
      do 120 i = 1, nrow
         iao(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = jao(ipos)
            if (indu(j).eq.0) then
!     .. new entry ..
               if (value2.ne.0) then
                  if (ao(ipos) .ne. 0.0D0) then
                     indu(j) = k
                     jao(k) = jao(ipos)
                     ao(k) = ao(ipos)
                     k = k + 1
                  endif
               else
                  indu(j) = k
                  jao(k) = jao(ipos)
                  k = k + 1
               endif
            else if (value2.ne.0) then
!     .. duplicate entry ..
               ao(indu(j)) = ao(indu(j)) + ao(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
!     .. remove marks before working on the next row ..
         do 110 ipos = iao(i), k - 1
            indu(jao(ipos)) = 0
 110     continue
 120  continue
      iao(nrow+1) = k
      if (job.le.1) return
!
!     .. partial ordering ..
!     split the matrix into strict upper/lower triangular
!     parts, INDU points to the the beginning of the strict upper part.
!
      do 140 i = 1, nrow
         klast = iao(i+1) - 1
         kfirst = iao(i)
 130     if (klast.gt.kfirst) then
            if (jao(klast).lt.i .and. jao(kfirst).ge.i) then
!     .. swap klast with kfirst ..
               j = jao(klast)
               jao(klast) = jao(kfirst)
               jao(kfirst) = j
               if (value2.ne.0) then
                  tmp = ao(klast)
                  ao(klast) = ao(kfirst)
                  ao(kfirst) = tmp
               endif
            endif
            if (jao(klast).ge.i) then
              klast = klast - 1
            endif
            if (jao(kfirst).lt.i) then
              kfirst = kfirst + 1
            endif
            go to 130
         endif
!
         if (jao(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
      if (job.le.2) return
!
!     .. order the entries according to column indices
!     bubble-sort is used
!
      do 190 i = 1, nrow
         do 160 ipos = iao(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  if (value2.ne.0) then
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i), iao(i+1)-1
            do 170 j = iao(i+1)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  if (value2.ne.0) then
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
 170        continue
 180     continue
 190  continue
!
      return
!---- end of ssrcsr ----------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine xssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)
      integer ia(nrow+1),iao(nrow+1),ja(*),jao(nzmax),indu(nrow+1)
      real*8 a(*),ao(nzmax)
!-----------------------------------------------------------------------
! Symmetric Sparse Row   to    (regular) Compressed Sparse Row
!----------------------------------------------------------------------- 
! this subroutine converts  a symmetric  matrix in which only the lower 
! part is  stored in compressed sparse row format, i.e.,
! a matrix stored in symmetric sparse format, into a fully stored matrix
! i.e., a matrix where both the lower and upper parts are stored in 
! compressed sparse row format. the algorithm is in place (i.e. result 
! may be overwritten onto the input matrix a, ja, ia ----- ). 
! the output matrix delivered by ssrcsr is such that each row starts with
! the elements of the lower part followed by those of the upper part.
!----------------------------------------------------------------------- 
! on entry:
!--------- 
!	
! nrow  = row dimension of inout matrix
! a, 
! ia, 
! ja    = matrix in compressed sparse row format. This is assumed to be
!         a lower triangular matrix. 
!
! nzmax	= size of arrays ao and jao. ssrcsr will abort if the storage 
!	   provided in a, ja is not sufficient to store A. See ierr. 
!	
! on return:
!----------
! ao, iao, 
!   jao = output matrix in compressed sparse row format. The resulting 
!         matrix is symmetric and is equal to A+A**T - D, if
!         A is the original lower triangular matrix. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.
!      
! indu  = integer array of length nrow+1. If the input matrix is such 
!         that the last element in each row is its diagonal element then
!         on return, indu will contain the pointers to the diagonal 
!         element in each row of the output matrix. Otherwise used as
!         work array.
! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
! 
!----------------------------------------------------------------------- 
      ierr = 0
      do 1 i=1,nrow+1
         indu(i) = 0     
 1    continue
!     
!     compute  number of elements in each row of strict upper part. 
!     put result in indu(i+1)  for row i. 
!     
      do 3 i=1, nrow
         do 2 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .lt. i) indu(j+1) = indu(j+1)+1
 2       continue 
 3    continue
!-----------
!     find addresses of first elements of ouput matrix. result in indu
!-----------
      indu(1) = 1 
      do 4 i=1,nrow
         lenrow = ia(i+1)-ia(i)
         indu(i+1) = indu(i) + indu(i+1) + lenrow
 4    continue
!--------------------- enough storage in a, ja ? --------
      nnz = indu(nrow+1)-1 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
!
!     now copy lower part (backwards).
!     
      kosav = indu(nrow+1)
      do 6 i=nrow,1,-1
         klast = ia(i+1)-1
         kfirst = ia(i)
         iao(i+1) = kosav
         ko = indu(i) 
         kosav = ko
         do 5 k = kfirst, klast
            ao(ko) = a(k)
            jao(ko) = ja(k)
	    ko = ko+1
 5       continue
         indu(i) = ko 
 6    continue
      iao(1) = 1
!
!     now copy upper part. Go through the structure of ao, jao, iao
!     that has already been copied (lower part). indu(i) is the address
!     of the next free location in row i for ao, jao.
!     
      do 8 i=1,nrow
!     i-th row is now in ao, jao, iao structure -- lower half part
         do 9 k=iao(i), iao(i+1)-1 
            j = jao(k)
            if (j .ge. i)  goto 8
            ipos = indu(j)
            ao(ipos) = ao(k)
            jao(ipos) = i
            indu(j) = indu(j) + 1 
 9       continue
 8    continue
      return
!----- end of xssrcsr -------------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine csrell (nrow,a,ja,ia,maxcol,coef,jcoef,ncoef,ndiag,ierr)
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)  
      real*8 a(*), coef(ncoef,1)
!----------------------------------------------------------------------- 
! Compressed Sparse Row	    to    Ellpack - Itpack format 
!----------------------------------------------------------------------- 
! this subroutine converts  matrix stored in the general a, ja, ia 
! format into the coef, jcoef itpack format.
!
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! nrow 	  = row dimension of the matrix A.
!
! a, 
! ia, 
! ja      = input matrix in compressed sparse row format. 
!
! ncoef  = first dimension of arrays coef, and jcoef.
! 
! maxcol = integer equal to the number of columns available in coef.
!
! on return: 
!----------
! coef	= real array containing the values of the matrix A in 
!         itpack-ellpack format.
! jcoef = integer array containing the column indices of coef(i,j) 
!         in A.
! ndiag = number of active 'diagonals' found. 
!
! ierr 	= error message. 0 = correct return. If ierr .ne. 0 on
!	  return this means that the number of diagonals found
!         (ndiag) exceeds maxcol.
!
!----------------------------------------------------------------------- 
! first determine the length of each row of lower-part-of(A)
      ierr = 0
      ndiag = 0
      do 3 i=1, nrow
         k = ia(i+1)-ia(i)
         ndiag = max0(ndiag,k) 
 3    continue
!----- check whether sufficient columns are available. ----------------- 
      if (ndiag .gt. maxcol) then
         ierr = 1 
         return
      endif
!
! fill coef with zero elements and jcoef with row numbers.------------ 
!
      do 4 j=1,ndiag 
         do 41 i=1,nrow
            coef(i,j) = 0.0d0
            jcoef(i,j) = i
 41      continue
 4    continue
!     
!------- copy elements row by row.-------------------------------------- 
!     
      do 6 i=1, nrow
         k1 = ia(i)
         k2 = ia(i+1)-1
         do 5 k=k1,k2
            coef(i,k-k1+1) = a(k)
            jcoef(i,k-k1+1) = ja(k)
 5       continue
 6    continue
      return
!--- end of csrell------------------------------------------------------ 
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine ellcsr(nrow,coef,jcoef,ncoef,ndiag,a,ja,ia,nzmax,ierr)
      integer ia(nrow+1), ja(*), jcoef(ncoef,1) 
      real*8 a(*), coef(ncoef,1)
!----------------------------------------------------------------------- 
!  Ellpack - Itpack format  to  Compressed Sparse Row
!----------------------------------------------------------------------- 
! this subroutine converts a matrix stored in ellpack-itpack format 
! coef-jcoef into the compressed sparse row format. It actually checks
! whether an entry in the input matrix is a nonzero element before
! putting it in the output matrix. The test does not account for small
! values but only for exact zeros. 
!----------------------------------------------------------------------- 
! on entry:
!---------- 
!
! nrow 	= row dimension of the matrix A.
! coef	= array containing the values of the matrix A in ellpack format.
! jcoef = integer arraycontains the column indices of coef(i,j) in A.
! ncoef = first dimension of arrays coef, and jcoef.
! ndiag = number of active columns in coef, jcoef.
! 
! ndiag = on entry the number of columns made available in coef.
!
! on return: 
!----------
! a, ia, 
!    ja = matrix in a, ia, ja format where. 
! 
! nzmax	= size of arrays a and ja. ellcsr will abort if the storage 
!	   provided in a, ja is not sufficient to store A. See ierr. 
!
! ierr 	= integer. serves are output error message. 
!         ierr = 0 means normal return. 
!         ierr = 1 means that there is not enough space in
!         a and ja to store output matrix.
!----------------------------------------------------------------------- 
! first determine the length of each row of lower-part-of(A)
      ierr = 0
!-----check whether sufficient columns are available. ----------------- 
!
!------- copy elements row by row.-------------------------------------- 
      kpos = 1
      ia(1) = kpos
      do 6 i=1, nrow
         do 5 k=1,ndiag
            if (coef(i,k) .ne. 0.0d0) then
               if (kpos .gt. nzmax) then
                  ierr = kpos
                  return
               endif
               a(kpos) = coef(i,k)
               ja(kpos) = jcoef(i,k)
               kpos = kpos+1
	    endif
 5       continue
         ia(i+1) = kpos
 6    continue	
      return
!--- end of ellcsr ----------------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)
      real*8 a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)
!----------------------------------------------------------------------- 
! Compressed Sparse Row   to      Modified - Sparse Row 
!                                 Sparse row with separate main diagonal
!----------------------------------------------------------------------- 
! converts a general sparse matrix a, ja, ia into 
! a compressed matrix using a separated diagonal (referred to as
! the bell-labs format as it is used by bell labs semi conductor
! group. We refer to it here as the modified sparse row format.
! Note: this has been coded in such a way that one can overwrite
! the output matrix onto the input matrix if desired by a call of
! the form 
!
!     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!
! In case ao, jao, are different from a, ja, then one can
! use ao, jao as the work arrays in the calling sequence:
!
!     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!
!----------------------------------------------------------------------- 
!
! on entry :
!---------
! a, ja, ia = matrix in csr format. note that the 
!	     algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!	 
! on return :
!-----------
!
! ao, jao  = sparse matrix in modified sparse row storage format:
!	   +  ao(1:n) contains the diagonal of the matrix. 
!	   +  ao(n+2:nnz) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!	   +  jao(n+2:nnz) : their column indices
!	   +  jao(1:n+1) contains the pointer array for the nondiagonal
!             elements in ao(n+1:nnz) and jao(n+2:nnz).
!             i.e., for i .le. n+1 jao(i) points to beginning of row i 
!	      in arrays ao, jao.
!	       here nnz = number of nonzero elements+1 
! work arrays:
!------------
! wk	= real work array of length n
! iwk   = integer work array of length n+1
!
! notes: 
!------- 
!        Algorithm is in place.  i.e. both:
!
!          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!          (in which  ao, jao, are different from a, ja)
!           and
!          call csrmsr (n, a, ja, ia, a, ja, wk,iwk) 
!          (in which  wk, jwk, are different from a, ja)
!        are OK.
!--------
! coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
!-----------------------------------------------------------------------
      icount = 0
!
! store away diagonal elements and count nonzero diagonal elements.
!
      do 1 i=1,n
         wk(i) = 0.0d0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1 
               iwk(i+1) = iwk(i+1)-1
            endif
 2       continue
 1    continue
!     
! compute total length
!     
      iptr = n + ia(n+1) - icount
!     
!     copy backwards (to avoid collisions)
!     
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr-1
            endif
 100     continue
 500  continue
!
! compute pointer values and copy wk(*)
!
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i) 
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return	
!------------ end of subroutine csrmsr ---------------------------------
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine msrcsr (n,a,ja,ao,jao,iao,wk,iwk)
      real*8 a(*),ao(*),wk(n)
      integer ja(*),jao(*),iao(n+1),iwk(n+1)
!----------------------------------------------------------------------- 
!       Modified - Sparse Row  to   Compressed Sparse Row   
!
!----------------------------------------------------------------------- 
! converts a compressed matrix using a separated diagonal 
! (modified sparse row format) in the Compressed Sparse Row   
! format.
! does not check for zero elements in the diagonal.
!
!
! on entry :
!---------
! n          = row dimension of matrix
! a, ja      = sparse matrix in msr sparse storage format
!              see routine csrmsr for details on data structure 
!        
! on return :
!-----------
!
! ao,jao,iao = output matrix in csr format.  
!
! work arrays:
!------------
! wk       = real work array of length n
! iwk      = integer work array of length n+1
!
! notes:
!   The original version of this was NOT in place, but has
!   been modified by adding the vector iwk to be in place.
!   The original version had ja instead of iwk everywhere in
!   loop 500.  Modified  Sun 29 May 1994 by R. Bramley (Indiana).
!   
!----------------------------------------------------------------------- 
      logical added
      do 1 i=1,n
         wk(i) = a(i)
         iwk(i) = ja(i)
 1    continue
      iwk(n+1) = ja(n+1)
      iao(1) = 1
      iptr = 1
!---------
      do 500 ii=1,n 
         added = .false.
         idiag = iptr + (iwk(ii+1)-iwk(ii)) 
         do 100 k=iwk(ii),iwk(ii+1)-1
            j = ja(k)
            if (j .lt. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            elseif (added) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            else 
! add diag element - only reserve a position for it. 
               idiag = iptr
               iptr = iptr+1
               added = .true.
!     then other element
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            endif
 100     continue
         ao(idiag) = wk(ii)
         jao(idiag) = ii
         if (.not. added) iptr = iptr+1
         iao(ii+1) = iptr 
 500  continue
      return    
!------------ end of subroutine msrcsr --------------------------------- 
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place. 
!----------------------------------------------------------------------- 
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n	= dimension of A.
! job	= integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
!	  for any other normal usage, enter ipos=1.
! a	= real array of length nnz (nnz=number of nonzero elements in input 
!         matrix) containing the nonzero elements.
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of
!	  the beginning of the k-th row.
!
! on return:
! ---------- 
! output arguments:
! ao	= real array of size nzz containing the "a" part of the transpose
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the "ia" index array of
!	  the transpose. 
!
!----------------------------------------------------------------------- 
      call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
!-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place. 
!----------------------------------------------------------------------- 
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!----------------------------------------------------------------------- 
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n	= number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job	= integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
!	  for any other normal usage, enter ipos=1.
! a	= real array of length nnz (nnz=number of nonzero elements in input 
!         matrix) containing the nonzero elements.
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of
!	  the beginning of the k-th row.
!
! on return:
! ---------- 
! output arguments:
! ao	= real array of size nzz containing the "a" part of the transpose
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the "ia" index array of
!	  the transpose. 
!
!----------------------------------------------------------------------- 
!----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
!---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
!--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
!-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
!--------------- end of csrcsc2 ---------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine csrlnk (n,a,ja,ia,link) 
      real*8 a(*) 
      integer n, ja(*), ia(n+1), link(*)
!----------------------------------------------------------------------- 
!      Compressed Sparse Row         to    Linked storage format. 
!----------------------------------------------------------------------- 
! this subroutine translates a matrix stored in compressed sparse
! row into one with a linked list storage format. Only the link
! array needs to be obtained since the arrays a, ja, and ia may
! be unchanged and  carry the same meaning for the output matrix.
! in  other words a, ja, ia, link   is the output linked list data
! structure with a, ja, unchanged from input, and ia possibly 
! altered (in case therea re null rows in matrix). Details on
! the output array link are given below.
!----------------------------------------------------------------------- 
! Coded by Y. Saad, Feb 21, 1991.
!----------------------------------------------------------------------- 
!
! on entry:
!----------
! n	= integer equal to the dimension of A.	
!         
! a	= real array of size nna containing the nonzero elements
! ja	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1 containing the pointers to the beginning 
!         of each row. ia(k) contains the position in a, ja of the 
!         beginning of the k-th row.
!
! on return:
!---------- 
! a, ja, are not changed.
! ia    may be changed if there are null rows.
! 
! a     = nonzero elements.
! ja    = column positions. 
! ia    = ia(i) points to the first element of row i in linked structure.
! link	= integer array of size containing the linked list information.
!         link(k) points to the next element of the row after element 
!         a(k), ja(k). if link(k) = 0, then there is no next element,
!         i.e., a(k), jcol(k) is the last element of the current row.
!
!  Thus row number i can be accessed as follows:
!     next = ia(i) 
!     while(next .ne. 0) do 
!          value = a(next)      ! value a(i,j) 
!          jcol  = ja(next)     ! column index j
!          next  = link(next)   ! address of next element in row
!     endwhile
! notes:
! ------ ia may be altered on return.
!----------------------------------------------------------------------- 
! local variables
      integer i, k
!
! loop through all rows
!
      do 100 i =1, n
         istart = ia(i) 
         iend = ia(i+1)-1
         if (iend .gt. istart) then
            do 99  k=istart, iend-1 
               link(k) = k+1
 99         continue
            link(iend) = 0
         else
            ia(i) = 0
         endif
 100  continue
!     
      return
!-------------end-of-csrlnk --------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine lnkcsr (n, a, jcol, istart, link, ao, jao, iao) 
      real*8 a(*), ao(*) 
      integer n, jcol(*), istart(n), link(*), jao(*), iao(*) 
!----------------------------------------------------------------------- 
!     Linked list storage format   to      Compressed Sparse Row  format
!----------------------------------------------------------------------- 
! this subroutine translates a matrix stored in linked list storage 
! format into the compressed sparse row format. 
!----------------------------------------------------------------------- 
! Coded by Y. Saad, Feb 21, 1991.
!----------------------------------------------------------------------- 
!
! on entry:
!----------
! n	= integer equal to the dimension of A.	
!         
! a	= real array of size nna containing the nonzero elements
! jcol	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! istart= integer array of size n poiting to the beginning of the rows.
!         istart(i) contains the position of the first element of 
!         row i in data structure. (a, jcol, link).
!         if a row is empty istart(i) must be zero.
! link	= integer array of size nnz containing the links in the linked 
!         list data structure. link(k) points to the next element 
!         of the row after element ao(k), jcol(k). if link(k) = 0, 
!         then there is no next element, i.e., ao(k), jcol(k) is 
!         the last element of the current row.
!
! on return:
!-----------
! ao, jao, iao = matrix stored in csr format:
!
! ao    = real array containing the values of the nonzero elements of 
!         the matrix stored row-wise. 
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the pointers array to the 
!         beginning of each row. iao(i) is the address in ao,jao of
!         first element of row i.
!
!----------------------------------------------------------------------- 
! first determine individial bandwidths and pointers.
!----------------------------------------------------------------------- 
! local variables
      integer irow, ipos, next
!-----------------------------------------------------------------------
      ipos = 1
      iao(1) = ipos
!     
!     loop through all rows
!     
      do 100 irow =1, n
!     
!     unroll i-th row.
!     
         next = istart(irow)
 10      if (next .eq. 0) goto 99
         jao(ipos) = jcol(next)
         ao(ipos)  = a(next)
         ipos = ipos+1
         next = link(next) 
         goto 10
 99      iao(irow+1) = ipos 
 100  continue
!     
      return
!-------------end-of-lnkcsr ------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine csrdia (n,idiag,job,a,ja,ia,ndiag,diag,ioff,ao,jao,iao,ind)
      real*8 diag(ndiag,idiag), a(*), ao(*)
      integer ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)
!----------------------------------------------------------------------- 
! Compressed sparse row     to    diagonal format
!----------------------------------------------------------------------- 
! this subroutine extracts  idiag diagonals  from the  input matrix a,
! a, ia, and puts the rest of  the matrix  in the  output matrix ao,
! jao, iao.  The diagonals to be extracted depend  on the  value of job
! (see below for details.)  In  the first  case, the  diagonals to be
! extracted are simply identified by  their offsets  provided in ioff
! by the caller.  In the second case, the  code internally determines
! the idiag most significant diagonals, i.e., those  diagonals of the
! matrix which  have  the  largest  number  of  nonzero elements, and
! extracts them.
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! n	= dimension of the matrix a.
! idiag = integer equal to the number of diagonals to be extracted. 
!         Note: on return idiag may be modified.
! a, ja, 			
!    ia = matrix stored in a, ja, ia, format
! job	= integer. serves as a job indicator.  Job is better thought 
!         of as a two-digit number job=xy. If the first (x) digit
!         is one on entry then the diagonals to be extracted are 
!         internally determined. In this case csrdia exctracts the
!         idiag most important diagonals, i.e. those having the largest
!         number on nonzero elements. If the first digit is zero
!         then csrdia assumes that ioff(*) contains the offsets 
!         of the diagonals to be extracted. there is no verification 
!         that ioff(*) contains valid entries.
!         The second (y) digit of job determines whether or not
!         the remainder of the matrix is to be written on ao,jao,iao.
!         If it is zero  then ao, jao, iao is not filled, i.e., 
!         the diagonals are found  and put in array diag and the rest is
!         is discarded. if it is one, ao, jao, iao contains matrix
!         of the remaining elements.
!         Thus:
!         job= 0 means do not select diagonals internally (pick those
!                defined by ioff) and do not fill ao,jao,iao
!         job= 1 means do not select diagonals internally 
!                      and fill ao,jao,iao
!         job=10 means  select diagonals internally 
!                      and do not fill ao,jao,iao
!         job=11 means select diagonals internally 
!                      and fill ao,jao,iao
! 
! ndiag = integer equal to the first dimension of array diag.
!
! on return:
!----------- 
!
! idiag = number of diagonals found. This may be smaller than its value 
!         on entry. 
! diag  = real array of size (ndiag x idiag) containing the diagonals
!         of A on return
!          
! ioff  = integer array of length idiag, containing the offsets of the
!   	  diagonals to be extracted.
! ao, jao
!  iao  = remainder of the matrix in a, ja, ia format.
! work arrays:
!------------ 
! ind   = integer array of length 2*n-1 used as integer work space.
!         needed only when job.ge.10 i.e., in case the diagonals are to
!         be selected internally.
!
! Notes:
!-------
!    1) The algorithm is in place: ao, jao, iao can be overwritten on 
!       a, ja, ia if desired 
!    2) When the code is required to select the diagonals (job .ge. 10) 
!       the selection of the diagonals is done from left to right 
!       as a result if several diagonals have the same weight (number 
!       of nonzero elemnts) the leftmost one is selected first.
!-----------------------------------------------------------------------
      job1 = job/10
      job2 = job-job1*10
      if (job1 .eq. 0) goto 50
      n2 = n+n-1
      call infdia(n,ja,ia,ind,idum)
!----------- determine diagonals to  accept.---------------------------- 
!----------------------------------------------------------------------- 
      ii = 0
 4    ii=ii+1
      jmax = 0
      do 41 k=1, n2
         j = ind(k)
         if (j .le. jmax) goto 41
         i = k
         jmax = j
 41   continue
      if (jmax .le. 0) then
         ii = ii-1
         goto 42
      endif
      ioff(ii) = i-n
      ind(i) = - jmax
      if (ii .lt.  idiag) goto 4
 42   idiag = ii
!---------------- initialize diago to zero ----------------------------- 
 50   continue
      do 55 j=1,idiag
         do 54 i=1,n
            diag(i,j) = 0.0d0
 54      continue
 55   continue
!----------------------------------------------------------------------- 
      ko = 1
!----------------------------------------------------------------------- 
! extract diagonals and accumulate remaining matrix.
!----------------------------------------------------------------------- 
      do 6 i=1, n
         do 51 k=ia(i),ia(i+1)-1 
            j = ja(k)
            do 52 l=1,idiag
               if (j-i .ne. ioff(l)) goto 52
               diag(i,l) = a(k)
               goto 51
 52         continue
!--------------- append element not in any diagonal to ao,jao,iao ----- 
            if (job2 .eq. 0) goto 51
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko+1
 51      continue
         if (job2 .ne. 0 ) ind(i+1) = ko
 6    continue
      if (job2 .eq. 0) return
!     finish with iao
      iao(1) = 1
      do 7 i=2,n+1
         iao(i) = ind(i)
 7    continue
      return
!----------- end of csrdia ---------------------------------------------
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine diacsr (n,job,idiag,diag,ndiag,ioff,a,ja,ia)
      real*8 diag(ndiag,idiag), a(*), t
      integer ia(*), ja(*), ioff(*)
!----------------------------------------------------------------------- 
!    diagonal format     to     compressed sparse row     
!----------------------------------------------------------------------- 
! this subroutine extract the idiag most important diagonals from the 
! input matrix a, ja, ia, i.e, those diagonals of the matrix which have
! the largest number of nonzero elements. If requested (see job),
! the rest of the matrix is put in a the output matrix ao, jao, iao
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! n	= integer. dimension of the matrix a.
! job	= integer. job indicator with the following meaning.
!         if (job .eq. 0) then check for each entry in diag
!         whether this entry is zero. If it is then do not include
!         in the output matrix. Note that the test is a test for
!         an exact arithmetic zero. Be sure that the zeros are
!         actual zeros in double precision otherwise this would not
!         work.
!         
! idiag = integer equal to the number of diagonals to be extracted. 
!         Note: on return idiag may be modified.
!
! diag  = real array of size (ndiag x idiag) containing the diagonals
!         of A on return. 
! 
! ndiag = integer equal to the first dimension of array diag.
!
! ioff  = integer array of length idiag, containing the offsets of the
!   	  diagonals to be extracted.
!
! on return:
!----------- 
! a, 
! ja, 			
! ia    = matrix stored in a, ja, ia, format
!
! Note:
! ----- the arrays a and ja should be of length n*idiag.
!
!----------------------------------------------------------------------- 
      ia(1) = 1
      ko = 1
      do 80 i=1, n
         do 70 jj = 1, idiag
            j = i+ioff(jj) 
            if (j .lt. 1 .or. j .gt. n) goto 70
            t = diag(i,jj) 
            if (job .eq. 0 .and. t .eq. 0.0d0) goto 70
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 70      continue
         ia(i+1) = ko
 80   continue
      return
!----------- end of diacsr ---------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine bsrcsr (job, n, m, na, a, ja, ia, ao, jao, iao)
      implicit none
      integer job, n, m, na, ia(*), ja(*), jao(*), iao(n+1)
      real*8 a(na,*), ao(*)
!-----------------------------------------------------------------------
!             Block Sparse Row  to Compressed Sparse Row.
!----------------------------------------------------------------------- 
! NOTE: ** meanings of parameters may have changed wrt earlier versions
! FORMAT DEFINITION HAS CHANGED WRT TO EARLIER VERSIONS... 
!-----------------------------------------------------------------------
!
! converts a  matrix stored in block-reduced   a, ja, ia  format to the
! general  sparse row a,  ja, ia  format.  A matrix   that has  a block
! structure is a matrix whose entries are blocks  of the same size m
! (e.g.  3 x 3).   Then it is often preferred  to work with the reduced
! graph of the matrix. Instead of storing one element at a time one can
! store a whole block at a time.  In this storage scheme  an entry is a
! square array holding the m**2 elements of a block.
! 
!-----------------------------------------------------------------------
! on entry:
!----------
! job   = if job.eq.0 on entry, values are not copied (pattern only)
!
! n	= the block row dimension of the matrix.
!
! m     = the dimension of each block. Thus, the actual row dimension 
!         of A is n x m.  
!
! na	= first dimension of array a as declared in calling program.
!         This should be .ge. m**2.
!
! a	= real array containing the real entries of the matrix. Recall
!         that each entry is in fact an m x m block. These entries 
!         are stored column-wise in locations a(1:m*m,k) for each k-th
!         entry. See details below.
! 
! ja	= integer array of length n. ja(k) contains the column index 
!         of the leading element, i.e., the element (1,1) of the block
!         that is held in the column a(*,k) of the value array. 
!
! ia    = integer array of length n+1. ia(i) points to the beginning
!         of block row number i in the arrays a and ja. 
! 
! on return:
!-----------
! ao, jao, 
! iao   = matrix stored in compressed sparse row format. The number of
!         rows in the new matrix is n x m. 
! 
! Notes: THIS CODE IS NOT IN PLACE.
! 
!-----------------------------------------------------------------------
! BSR FORMAT.
!---------- 
! Each row of A contains the m x m block matrix unpacked column-
! wise (this allows the user to declare the array a as a(m,m,*) on entry 
! if desired). The block rows are stored in sequence just as for the
! compressed sparse row format. 
!
!-----------------------------------------------------------------------
!     example  with m = 2:
!                                                       1  2 3   
!    +-------|--------|--------+                       +-------+
!    | 1   2 |  0   0 |  3   4 |     Block             | x 0 x | 1
!    | 5   6 |  0   0 |  7   8 |     Representation:   | 0 x x | 2 
!    +-------+--------+--------+                       | x 0 0 | 3 
!    | 0   0 |  9  10 | 11  12 |                       +-------+ 
!    | 0   0 | 13  14 | 15  16 |  
!    +-------+--------+--------+   
!    | 17 18 |  0   0 |  0   0 |
!    | 22 23 |  0   0 |  0   0 |
!    +-------+--------+--------+
!
!    For this matrix:     n    = 3
!                         m    = 2
!                         nnz  = 5 
!-----------------------------------------------------------------------
! Data structure in Block Sparse Row format:      
!-------------------------------------------
! Array A:
!------------------------- 
!     1   3   9   11   17   <<--each m x m block is stored column-wise 
!     5   7   13  15   22       in a  column of the array A.
!     2   4   10  12   18      
!     6   8   14  16   23
!------------------------- 
! JA  1   3   2    3    1   <<-- column indices for each block. Note that
!-------------------------       these indices are wrt block matrix.
! IA  1   3   5    6        <<-- pointers to beginning of each block row 
!-------------------------       in arrays A and JA. 
!-----------------------------------------------------------------------
! locals 
! 
      integer i, i1, i2, ij, ii, irow, j, jstart, k, krow, no
      logical val
!     
      val = (job.ne.0)
      no = n * m 
      irow = 1	
      krow = 1
      iao(irow) = 1      
!-----------------------------------------------------------------------
      do 2 ii=1, n
!
!     recall: n is the block-row dimension
!
         i1 = ia(ii)
         i2 = ia(ii+1)-1
!
!     create m rows for each block row -- i.e., each k. 
!     
         do 23 i=1,m
            do 21 k=i1, i2
               jstart = m*(ja(k)-1)
               do 22  j=1,m
                  ij = (j-1)*m + i
                  if (val) ao(krow) = a(ij,k) 
                  jao(krow) = jstart+j
                  krow = krow+1
 22            continue	    
 21         continue
            irow = irow+1 
            iao(irow) = krow
 23      continue
 2    continue
      return
!-------------end-of-bsrcsr -------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine csrbsr (job,nrow,m,na,a,ja,ia,ao,jao,iao,iw,ierr)
      implicit none
      integer job,ierr,nrow,m,na,ia(nrow+1),ja(*),jao(na),iao(*),iw(*)
      real*8 a(*),ao(na,*)
!-----------------------------------------------------------------------
!     Compressed Sparse Row  to    Block Sparse Row
!-----------------------------------------------------------------------
!
! This  subroutine converts a matrix stored  in a general compressed a,
! ja, ia format into a a block  sparse row format a(m,m,*),ja(*),ia(*).
! See routine  bsrcsr  for more  details on  data   structure for block
! matrices. 
!
! NOTES: 1) the initial matrix does not have to have a block structure. 
! zero padding is done for general sparse matrices. 
!        2) For most practical purposes, na should be the same as m*m.
! 
!-----------------------------------------------------------------------
! 
! In what follows nr=1+(nrow-1)/m = block-row dimension of output matrix 
! 
! on entry:
!----------
! 
! job   =  job indicator.
!          job =  0 -> only the pattern of output matrix is generated 
!          job >  0 -> both pattern and values are generated. 
!          job = -1 -> iao(1) will return the number of nonzero blocks,
!            in the output matrix. In this case jao(1:nr) is used as 
!            workspace, ao is untouched, iao is untouched except iao(1)
!
! nrow	= integer, the actual row dimension of the matrix.
! 
! m     = integer equal to the dimension of each block. m should be > 0. 
! 
! na	= first dimension of array ao as declared in calling program.
!         na should be .ge. m*m. 
!
! a, ja, 
!    ia = input matrix stored in compressed sparse row format.
!
! on return:
!-----------
! 
! ao    = real  array containing the  values of the matrix. For details
!         on the format  see below. Each  row of  a contains the  m x m
!         block matrix  unpacked column-wise (this  allows the  user to
!         declare the  array a as ao(m,m,*) on  entry if desired).  The
!         block rows are stored in sequence  just as for the compressed
!         sparse row format. The block  dimension  of the output matrix
!         is  nr = 1 + (nrow-1) / m.
!         
! jao   = integer array. containing the block-column indices of the 
!         block-matrix. Each jao(k) is an integer between 1 and nr
!         containing the block column index of the block ao(*,k).  
!
! iao   = integer array of length nr+1. iao(i) points to the beginning
!         of block row number i in the arrays ao and jao. When job=-1
!         iao(1) contains the number of nonzero blocks of the output
!         matrix and the rest of iao is unused. This is useful for
!         determining the lengths of ao and jao. 
!
! ierr  = integer, error code. 
!              0 -- normal termination
!              1 -- m is equal to zero 
!              2 -- NA too small to hold the blocks (should be .ge. m**2)
!
! Work arrays:
!------------- 
! iw    = integer work array of dimension  nr = 1 + (nrow-1) / m
!
! NOTES: 
!-------
!     1) this code is not in place.
!     2) see routine bsrcsr for details on data sctructure for block 
!        sparse row format.
!     
!-----------------------------------------------------------------------
!     nr is the block-dimension of the output matrix.
!     
      integer nr, m2, io, ko, ii, len, k, jpos, j, i, ij, jr, irow   
      logical vals  
!----- 
      ierr = 0 
      if (m*m .gt. na) ierr = 2 
      if (m .eq. 0) ierr = 1 
      if (ierr .ne. 0) return
!----------------------------------------------------------------------- 
      vals = (job .gt. 0) 
      nr = 1 + (nrow-1) / m
      m2 = m*m 
      ko = 1 
      io = 1 
      iao(io) = 1 
      len = 0 
!     
!     iw determines structure of block-row (nonzero indicator) 
! 
         do j=1, nr
            iw(j) = 0
         enddo
!     
!     big loop -- leap by m rows each time.
!     
      do ii=1, nrow, m
         irow = 0
!
!     go through next m rows -- make sure not to go beyond nrow. 
!
         do while (ii+irow .le. nrow .and. irow .le. m-1) 
            do k=ia(ii+irow),ia(ii+irow+1)-1
!     
!     block column index = (scalar column index -1) / m + 1 
!
               j = ja(k)-1 
               jr = j/m + 1                
               j = j - (jr-1)*m 
               jpos = iw(jr) 
               if (jpos .eq. 0) then
!
!     create a new block
!     
                  iw(jr) = ko 
                  jao(ko) = jr 
                  if (vals) then
!     
!     initialize new block to zero -- then copy nonzero element
!     
                     do i=1, m2
                        ao(i,ko) = 0.0d0
                     enddo
                     ij = j*m + irow + 1 
                     ao(ij,ko) = a(k) 
                  endif
                  ko = ko+1
               else
!
!     copy column index and nonzero element 
!     
                  jao(jpos) = jr 
                  ij = j*m + irow + 1 
                  if (vals) ao(ij,jpos) = a(k) 
               endif 
            enddo  
            irow = irow+1
         enddo
!     
!     refresh iw
!                      
         do j = iao(io),ko-1 
            iw(jao(j)) = 0
         enddo
         if (job .eq. -1) then
            len = len + ko-1
            ko = 1
         else
            io = io+1 
            iao(io) = ko
         endif
      enddo
      if (job .eq. -1) iao(1) = len
!
      return
!--------------end-of-csrbsr-------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine csrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)
      real*8 a(*),abd(nabd,n)
      integer ia(n+1),ja(*)
!----------------------------------------------------------------------- 
!   Compressed Sparse Row  to  Banded (Linpack ) format.
!----------------------------------------------------------------------- 
! this subroutine converts a general sparse matrix stored in
! compressed sparse row format into the banded format. for the
! banded format,the Linpack conventions are assumed (see below).
!----------------------------------------------------------------------- 
! on entry:
!----------
! n	= integer,the actual row dimension of the matrix.
!
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.
!
! job	= integer. if job=1 then the values of the lower bandwith ml 
!         and the upper bandwidth mu are determined internally. 
!         otherwise it is assumed that the values of ml and mu 
!         are the correct bandwidths on input. See ml and mu below.
!
! nabd  = integer. first dimension of array abd.
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located. 
!         lowd should be  ( 1  .le.  lowd  .le. nabd).
!         if it is not known in advance what lowd should be
!         enter lowd = 0 and the default value lowd = ml+mu+1
!         will be chosen. Alternative: call routine getbwd from unary
!         first to detrermione ml and mu then define lowd accordingly.
!         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
!
! ml	= integer. equal to the bandwidth of the strict lower part of A
! mu	= integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than lowd then an error 
!         flag is raised (unless lowd = 0). see ierr.
!
! note:   ml and mu are assumed to have	 the correct bandwidth values
!         as defined above if job is set to zero on entry.
!
! on return:
!-----------
!
! abd   = real array of dimension abd(nabd,n).
!         on return contains the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
!         the bottom row (row lowd). See details below for this format.
!
! ml	= integer. equal to the bandwidth of the strict lower part of A
! mu	= integer. equal to the bandwidth of the strict upper part of A
!         if job=1 on entry then these two values are internally computed.
!
! lowd  = integer. row number in abd where the lowest diagonal 
!         (leftmost) of A is located on return. In case lowd = 0
!         on return, then it is defined to ml+mu+1 on return and the
!         lowd will contain this value on return. `
!
! ierr  = integer. used for error messages. On return:
!         ierr .eq. 0  :means normal return
!         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
!         or larger than nabd).
!         ierr .eq. -2 : means that lowd is not large enough and as 
!         result the matrix cannot be stored in array abd. 
!         lowd should be at least ml+mu+1, where ml and mu are as
!         provided on output.
!
!----------------------------------------------------------------------* 
! Additional details on banded format.  (this closely follows the      *
! format used in linpack. may be useful for converting a matrix into   *
! this storage format in order to use the linpack  banded solvers).    * 
!----------------------------------------------------------------------*
!             ---  band storage format  for matrix abd ---             * 
! uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           *
! a in rows of abd starting from the lowest (sub)-diagonal  which  is  *
! stored in row number lowd of abd. the minimum number of rows needed  *
! in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  *
! j-th  column  of  abd contains the elements of the j-th column of a, *
! from bottom to top: the element a(j+ml,j) is stored in  position     *
! abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   *
! Generally, the element a(j+k,j) of original matrix a is stored in    *
! position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               *
! The first dimension nabd of abd must be .ge. lowd                    *
!                                                                      *
!     example [from linpack ]:   if the original matrix is             *
!                                                                      *
!              11 12 13  0  0  0                                       *
!              21 22 23 24  0  0                                       *
!               0 32 33 34 35  0     original banded matrix            *
!               0  0 43 44 45 46                                       *
!               0  0  0 54 55 56                                       *
!               0  0  0  0 65 66                                       *
!                                                                      *
! then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   *
! if lowd = 5 for example, abd  should be:                             *
!                                                                      *
! untouched --> x  x  x  x  x  x                                       *
!               *  * 13 24 35 46                                       *
!               * 12 23 34 45 56    resulting abd matrix in banded     *
!              11 22 33 44 55 66    format                             *
!  row lowd--> 21 32 43 54 65  *                                       *
!                                                                      *
! * = not used                                                         *
!                                                                      
*
!----------------------------------------------------------------------*
! first determine ml and mu.
!----------------------------------------------------------------------- 
      ierr = 0
!-----------
      if (job .eq. 1) call getbwd(n,a,ja,ia,ml,mu)
      m = ml+mu+1
      if (lowd .eq. 0) lowd = m
      if (m .gt. lowd)  ierr = -2
      if (lowd .gt. nabd .or. lowd .lt. 0) ierr = -1
      if (ierr .lt. 0) return
!------------
      do 15  i=1,m
         ii = lowd -i+1
         do 10 j=1,n
	    abd(ii,j) = 0.0d0
 10      continue
 15   continue
!---------------------------------------------------------------------	   
      mdiag = lowd-ml
      do 30 i=1,n
         do 20 k=ia(i),ia(i+1)-1
            j = ja(k)
            abd(i-j+mdiag,j) = a(k) 
 20      continue
 30   continue
      return
!------------- end of csrbnd ------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)
      real*8 a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)
!----------------------------------------------------------------------- 
! Banded (Linpack ) format   to    Compressed Sparse Row  format.
!----------------------------------------------------------------------- 
! on entry:
!----------
! n	= integer,the actual row dimension of the matrix.
!
! nabd  = first dimension of array abd.
!
! abd   = real array containing the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix,comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
!         in row lowd (see below). 
!    
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located. 
!         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
!         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.
!
! ml	= integer. equal to the bandwidth of the strict lower part of A
! mu	= integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than nabd then an error 
!         message is set. see ierr.
!
! len   = integer. length of arrays a and ja. bndcsr will stop if the
!         length of the arrays a and ja is insufficient to store the 
!         matrix. see ierr.
!
! on return:
!-----------
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.
!
! lowd  = if on entry lowd was zero then lowd is reset to the default
!         value ml+mu+l. 
!
! ierr  = integer. used for error message output. 
!         ierr .eq. 0 :means normal return
!         ierr .eq. -1 : means invalid value for lowd. 
!	  ierr .gt. 0 : means that there was not enough storage in a and ja
!         for storing the ourput matrix. The process ran out of space 
!         (as indicated by len) while trying to fill row number ierr. 
!         This should give an idea of much more storage might be required. 
!         Moreover, the first irow-1 rows are correctly filled. 
!
! notes:  the values in abd found to be equal to zero
! -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
!         The resulting may not be identical to a csr matrix
!         originally transformed to a bnd format.
!          
!----------------------------------------------------------------------- 
      ierr = 0
!-----------
      if (lowd .gt. nabd .or. lowd .le. 0) then 
         ierr = -1
         return
      endif
!-----------
      ko = 1
      ia(1) = 1
      do 30 irow=1,n
!-----------------------------------------------------------------------
         i = lowd 
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j) 
             if (t .eq. 0.0d0) goto 19
             if (ko .gt. len) then 
               ierr = irow 
               return
            endif
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
!     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
!------------- end of bndcsr ------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine csrssk (n,imod,a,ja,ia,asky,isky,nzmax,ierr)
      real*8 a(*),asky(nzmax) 
      integer n, imod, nzmax, ierr, ia(n+1), isky(n+1), ja(*)
!----------------------------------------------------------------------- 
!      Compressed Sparse Row         to     Symmetric Skyline Format 
!  or  Symmetric Sparse Row        
!----------------------------------------------------------------------- 
! this subroutine translates a compressed sparse row or a symmetric
! sparse row format into a symmetric skyline format.
! the input matrix can be in either compressed sparse row or the 
! symmetric sparse row format. The output matrix is in a symmetric
! skyline format: a real array containing the (active portions) of the
! rows in  sequence and a pointer to the beginning of each row.
!
! This module is NOT  in place.
!----------------------------------------------------------------------- 
! Coded by Y. Saad, Oct 5, 1989. Revised Feb. 18, 1991.
!----------------------------------------------------------------------- 
!
! on entry:
!----------
! n	= integer equal to the dimension of A.	
! imod  = integer indicating the variant of skyline format wanted:
!         imod = 0 means the pointer isky points to the `zeroth' 
!         element of the row, i.e., to the position of the diagonal
!         element of previous row (for i=1, isky(1)= 0)
!         imod = 1 means that itpr points to the beginning of the row. 
!         imod = 2 means that isky points to the end of the row (diagonal
!                  element) 
!         
! a	= real array of size nna containing the nonzero elements
! ja	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of
!	  the beginning of the k-th row.
! nzmax = integer. must be set to the number of available locations
!         in the output array asky. 
!
! on return:
!---------- 
!
! asky    = real array containing the values of the matrix stored in skyline
!         format. asky contains the sequence of active rows from 
!         i=1, to n, an active row being the row of elemnts of 
!         the matrix contained between the leftmost nonzero element 
!         and the diagonal element. 
! isky	= integer array of size n+1 containing the pointer array to 
!         each row. The meaning of isky depends on the input value of
!         imod (see above). 
! ierr  =  integer.  Error message. If the length of the 
!         output array asky exceeds nzmax. ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
! 
! Notes:
!         1) This module is NOT  in place.
!         2) even when imod = 2, length of  isky is  n+1, not n.
!
!----------------------------------------------------------------------- 
! first determine individial bandwidths and pointers.
!----------------------------------------------------------------------- 
      ierr = 0
      isky(1) = 0
      do 3 i=1,n
         ml = 0
         do 31 k=ia(i),ia(i+1)-1 
            ml = max(ml,i-ja(k)+1) 
 31      continue
         isky(i+1) = isky(i)+ml
 3    continue
!
!     test if there is enough space  asky to do the copying.  
!
      nnz = isky(n+1) 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
!    
!   fill asky with zeros.
!     
      do 1 k=1, nnz 
         asky(k) = 0.0d0
 1    continue
!     
!     copy nonzero elements.
!     
      do 4 i=1,n
         kend = isky(i+1) 
         do 41 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .le. i) asky(kend+j-i) = a(k)
 41      continue
 4    continue
! 
! modify pointer according to imod if necessary.
!
      if (imod .eq. 0) return
      if (imod .eq. 1) then 
         do 50 k=1, n+1
            isky(k) = isky(k)+1
 50      continue
      endif
      if (imod .eq. 2) then
         do 60 k=1, n
            isky(k) = isky(k+1) 
 60      continue
      endif
!
      return
!------------- end of csrssk ------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!-----------------------------------------------------------------------
      subroutine sskssr (n,imod,asky,isky,ao,jao,iao,nzmax,ierr)
      real*8 asky(*),ao(nzmax) 
      integer n, imod,nzmax,ierr, isky(n+1),iao(n+1),jao(nzmax) 
!----------------------------------------------------------------------- 
!     Symmetric Skyline Format  to  Symmetric Sparse Row format.
!----------------------------------------------------------------------- 
!  tests for exact zeros in skyline matrix (and ignores them in
!  output matrix).  In place routine (a, isky :: ao, iao)
!----------------------------------------------------------------------- 
! this subroutine translates a  symmetric skyline format into a 
! symmetric sparse row format. Each element is tested to see if it is
! a zero element. Only the actual nonzero elements are retained. Note 
! that the test used is simple and does take into account the smallness 
! of a value. the subroutine filter (see unary module) can be used
! for this purpose. 
!----------------------------------------------------------------------- 
! Coded by Y. Saad, Oct 5, 1989. Revised Feb 18, 1991./
!----------------------------------------------------------------------- 
!
! on entry:
!----------
! n	= integer equal to the dimension of A.	
! imod  = integer indicating the variant of skyline format used:
!         imod = 0 means the pointer iao points to the `zeroth' 
!         element of the row, i.e., to the position of the diagonal
!         element of previous row (for i=1, iao(1)= 0)
!         imod = 1 means that itpr points to the beginning of the row. 
!         imod = 2 means that iao points to the end of the row 
!                  (diagonal element) 
! asky  = real array containing the values of the matrix. asky contains 
!         the sequence of active rows from i=1, to n, an active row 
!         being the row of elemnts of the matrix contained between the 
!         leftmost nonzero element and the diagonal element. 
! isky 	= integer array of size n+1 containing the pointer array to 
!         each row. isky (k) contains the address of the beginning of the
!         k-th active row in the array asky. 
! nzmax = integer. equal to the number of available locations in the 
!         output array ao.  
!
! on return:
! ---------- 
! ao	= real array of size nna containing the nonzero elements
! jao	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! iao	= integer of size n+1. iao(k) contains the position in a, ja of
!	  the beginning of the k-th row.
! ierr  = integer. Serving as error message. If the length of the 
!         output arrays ao, jao exceeds nzmax then ierr returns 
!         the row number where the algorithm stopped: rows
!         i, to ierr-1 have been processed succesfully.
!         ierr = 0 means normal return.
!         ierr = -1  : illegal value for imod
! Notes:
!------- 
! This module is in place: ao and iao can be the same as asky, and isky.
!-----------------------------------------------------------------------
! local variables
      integer next, kend, kstart, i, j 
      ierr = 0
!
! check for validity of imod
! 
      if (imod.ne.0 .and. imod.ne.1 .and. imod .ne. 2) then
         ierr =-1
         return
      endif 
!
! next  = pointer to next available position in output matrix
! kend  = pointer to end of current row in skyline matrix. 
!
      next = 1
!
! set kend = start position -1 in  skyline matrix.
! 
      kend = 0 
      if (imod .eq. 1) kend = isky(1)-1
      if (imod .eq. 0) kend = isky(1) 
!
! loop through all rows
!     
      do 50 i=1,n
!
! save value of pointer to ith row in output matrix
!
         iao(i) = next
!
! get beginnning and end of skyline  row 
!
         kstart = kend+1
         if (imod .eq. 0) kend = isky(i+1)
         if (imod .eq. 1) kend = isky(i+1)-1
         if (imod .eq. 2) kend = isky(i) 
! 
! copy element into output matrix unless it is a zero element.
! 
         do 40 k=kstart,kend
            if (asky(k) .eq. 0.0d0) goto 40
            j = i-(kend-k) 
            jao(next) = j
            ao(next)  = asky(k)
            next=next+1
            if (next .gt. nzmax+1) then
               ierr = i
               return
            endif 
 40      continue
 50    continue
      iao(n+1) = next
      return
!-------------end-of-sskssr -------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine csrjad (nrow, a, ja, ia, idiag, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(nrow+1), iperm(nrow), iao(nrow) 
      real*8 a(*), ao(*)
!-----------------------------------------------------------------------
!    Compressed Sparse Row  to   JAgged Diagonal storage. 
!----------------------------------------------------------------------- 
! this subroutine converts  matrix stored in the compressed sparse
! row format to the jagged diagonal format. The data structure
! for the JAD (Jagged Diagonal storage) is as follows. The rows of 
! the matrix are (implicitly) permuted so that their lengths are in
! decreasing order. The real entries ao(*) and their column indices 
! jao(*) are stored in succession. The number of such diagonals is idiag.
! the lengths of each of these diagonals is stored in iao(*).
! For more details see [E. Anderson and Y. Saad,
! ``Solving sparse triangular systems on parallel computers'' in
! Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
! or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
! SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! nrow 	  = row dimension of the matrix A.
!
! a, 
! ia, 
! ja      = input matrix in compressed sparse row format. 
!
! on return: 
!----------
! 
! idiag = integer. The number of jagged diagonals in the matrix.
!
! iperm = integer array of length nrow containing the permutation
!         of the rows that leads to a decreasing order of the
!         number of nonzero elements.
!
! ao    = real array containing the values of the matrix A in 
!         jagged diagonal storage. The j-diagonals are stored
!         in ao in sequence. 
!
! jao   = integer array containing the column indices of the 
!         entries in ao.
!
! iao   = integer array containing pointers to the beginning 
!         of each j-diagonal in ao, jao. iao is also used as 
!         a work array and it should be of length n at least.
!
!----------------------------------------------------------------------- 
!     ---- define initial iperm and get lengths of each row
!     ---- jao is used a work vector to store tehse lengths
!     
      idiag = 0
      ilo = nrow 
      do 10 j=1, nrow
         iperm(j) = j 
         len = ia(j+1) - ia(j)
         ilo = min(ilo,len) 
         idiag = max(idiag,len) 
         jao(j) = len
 10   continue 
!     
!     call sorter to get permutation. use iao as work array.
!    
      call dcsort (jao, nrow, iao, iperm, ilo, idiag) 
!     
!     define output data structure. first lengths of j-diagonals
!     
      do 20 j=1, nrow
         iao(j) = 0
 20   continue
      do 40 k=1, nrow
         len = jao(iperm(k)) 
         do 30 i=1,len
            iao(i) = iao(i)+1
 30      continue
 40   continue
!     
!     get the output matrix itself
!     
      k1 = 1
      k0 = k1
      do 60 jj=1, idiag
         len = iao(jj)
         do 50 k=1,len
            i = ia(iperm(k))+jj-1
            ao(k1) = a(i)
            jao(k1) = ja(i) 
            k1 = k1+1
 50      continue
         iao(jj) = k0
         k0 = k1
 60   continue
      iao(idiag+1) = k1
      return
!----------end-of-csrjad------------------------------------------------
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine jadcsr (nrow, idiag, a, ja, ia, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(idiag+1), iperm(nrow), iao(nrow+1) 
      real*8 a(*), ao(*)
!-----------------------------------------------------------------------
!     Jagged Diagonal Storage   to     Compressed Sparse Row  
!----------------------------------------------------------------------- 
! this subroutine converts a matrix stored in the jagged diagonal format
! to the compressed sparse row format.
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! nrow 	  = integer. the row dimension of the matrix A.
! 
! idiag   = integer. The  number of jagged diagonals in the data
!           structure a, ja, ia.
! 
! a, 
! ja,
! ia      = input matrix in jagged diagonal format. 
!
! iperm   = permutation of the rows used to obtain the JAD ordering. 
!        
! on return: 
!----------
! 
! ao, jao,
! iao     = matrix in CSR format.
!-----------------------------------------------------------------------
! determine first the pointers for output matrix. Go through the
! structure once:
!
      do 137 j=1,nrow
         jao(j) = 0
 137  continue
!     
!     compute the lengths of each row of output matrix - 
!     
      do 140 i=1, idiag
         len = ia(i+1)-ia(i) 
         do 138 k=1,len
            jao(iperm(k)) = jao(iperm(k))+1
 138     continue
 140  continue
!     
!     remember to permute
!     
      kpos = 1
      iao(1) = 1
      do 141 i=1, nrow 
         kpos = kpos+jao(i) 
         iao(i+1) = kpos
 141  continue
!     
!     copy elemnts one at a time.
!     
      do 200 jj = 1, idiag
         k1 = ia(jj)-1
         len = ia(jj+1)-k1-1 
         do 160 k=1,len
            kpos = iao(iperm(k))
            ao(kpos) = a(k1+k) 
            jao(kpos) = ja(k1+k) 
            iao(iperm(k)) = kpos+1
 160     continue
 200  continue
!     
!     rewind pointers
!     
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!----------end-of-jadcsr------------------------------------------------
!-----------------------------------------------------------------------
      end
      subroutine dcsort(ival, n, icnt, index, ilo, ihi)
!-----------------------------------------------------------------------
!     Specifications for arguments:
!     ----------------------------
      integer n, ilo, ihi, ival(n), icnt(ilo:ihi), index(n)
!-----------------------------------------------------------------------
!    This routine computes a permutation which, when applied to the
!    input vector ival, sorts the integers in ival in descending
!    order.  The permutation is represented by the vector index.  The
!    permuted ival can be interpreted as follows:
!      ival(index(i-1)) .ge. ival(index(i)) .ge. ival(index(i+1))
!
!    A specialized sort, the distribution counting sort, is used 
!    which takes advantage of the knowledge that
!        1)  The values are in the (small) range [ ilo, ihi ]
!        2)  Values are likely to be repeated often
!
!    contributed to SPARSKIT by Mike Heroux. (Cray Research) 
!    --------------------------------------- 
!----------------------------------------------------------------------- 
! Usage:
!------ 
!     call dcsort( ival, n, icnt, index, ilo, ihi )
!
! Arguments:
!----------- 
!    ival  integer array (input)
!          On entry, ia is an n dimensional array that contains
!          the values to be sorted.  ival is unchanged on exit.
!
!    n     integer (input)
!          On entry, n is the number of elements in ival and index.
!
!    icnt  integer (work)
!          On entry, is an integer work vector of length 
!          (ihi - ilo + 1).
!
!    index integer array (output)
!          On exit, index is an n-length integer vector containing
!          the permutation which sorts the vector ival.
!
!    ilo   integer (input)
!          On entry, ilo is .le. to the minimum value in ival.
!
!    ihi   integer (input)
!          On entry, ihi is .ge. to the maximum value in ival.
!
! Remarks:
!--------- 
!         The permutation is NOT applied to the vector ival.
!
!----------------------------------------------------------------
!
! Local variables:
!    Other integer values are temporary indices.
!
! Author: 
!-------- 
!    Michael Heroux
!    Sandra Carney
!       Mathematical Software Research Group
!       Cray Research, Inc.
!
! References:
!    Knuth, Donald E., "The Art of Computer Programming, Volume 3:
!    Sorting and Searching," Addison-Wesley, Reading, Massachusetts,
!    1973, pp. 78-79.
!
! Revision history:
!    05/09/90: Original implementation.  A variation of the 
!              Distribution Counting Sort recommended by
!              Sandra Carney. (Mike Heroux)
!
!-----------------------------------------------------------------
!     ----------------------------------
!     Specifications for local variables
!     ----------------------------------
      integer i, j, ivalj
!
!     --------------------------
!     First executable statement
!     --------------------------
      do 10 i = ilo, ihi
        icnt(i) = 0
 10   continue
!
      do 20 i = 1, n
        icnt(ival(i)) = icnt(ival(i)) + 1
 20   continue
!
      do 30 i = ihi-1,ilo,-1
        icnt(i) = icnt(i) + icnt(i+1)
 30   continue
!
      do 40 j = n, 1, -1
        ivalj = ival(j)
        index(icnt(ivalj)) = j
        icnt(ivalj) = icnt(ivalj) - 1
 40   continue
      return
      end
!-------end-of-dcsort---------------------------------------------------
!-----------------------------------------------------------------------
      subroutine cooell(job,n,nnz,a,ja,ia,ao,jao,lda,ncmax,nc,ierr)
      implicit none
      integer job,n,nnz,lda,ncmax,nc,ierr
      integer ja(nnz),ia(nnz),jao(lda,ncmax)
      real*8  a(nnz),ao(lda,ncmax)
!-----------------------------------------------------------------------
!     COOrdinate format to ELLpack format
!-----------------------------------------------------------------------
!     On entry:
!     job     -- 0 if only pattern is to be processed(AO is not touched)
!     n       -- number of rows in the matrix
!     a,ja,ia -- input matix in COO format
!     lda     -- leading dimension of array AO and JAO
!     ncmax   -- size of the second dimension of array AO and JAO
!
!     On exit:
!     ao,jao  -- the matrix in ELL format
!     nc      -- maximum number of nonzeros per row
!     ierr    -- 0 if convertion succeeded
!                -1 if LDA < N
!                nc if NC > ncmax
!
!     NOTE: the last column of JAO is used as work space!!
!-----------------------------------------------------------------------
      integer i,j,k,ip
      real*8  zero
      logical copyval
      parameter (zero=0.0D0)
!     .. first executable statement ..
      copyval = (job.ne.0)
      if (lda .lt. n) then
         ierr = -1
         return
      endif
!     .. use the last column of JAO as workspace
!     .. initialize the work space
      do i = 1, n
         jao(i,ncmax) = 0
      enddo
      nc = 0
!     .. go through ia and ja to find out number nonzero per row
      do k = 1, nnz
         i = ia(k)
         jao(i,ncmax) = jao(i,ncmax) + 1
      enddo
!     .. maximum number of nonzero per row
      nc = 0
      do i = 1, n
         if (nc.lt.jao(i,ncmax)) nc = jao(i,ncmax)
         jao(i,ncmax) = 0
      enddo
!     .. if nc > ncmax retrun now
      if (nc.gt.ncmax) then
         ierr = nc
         return
      endif
!     .. go through ia and ja to copy the matrix to AO and JAO
      do k = 1, nnz
         i = ia(k)
         j = ja(k)
         jao(i,ncmax) = jao(i,ncmax) + 1
         ip = jao(i,ncmax)
         if (ip.gt.nc) nc = ip
         if (copyval) ao(i,ip) = a(k)
         jao(i,ip) = j
      enddo
!     .. fill the unspecified elements of AO and JAO with zero diagonals
      do i = 1, n
         do j = ia(i+1)-ia(i)+1, nc
            jao(i,j)=i
            if(copyval) ao(i,j) = zero
         enddo
      enddo
      ierr = 0
!
      return
      end
!-----end-of-cooell-----------------------------------------------------
!-----------------------------------------------------------------------
      subroutine xcooell(n,nnz,a,ja,ia,ac,jac,nac,ner,ncmax,ierr)
!-----------------------------------------------------------------------
!   coordinate format to ellpack format.
!-----------------------------------------------------------------------
!
!   DATE WRITTEN: June 4, 1989. 
!
!   PURPOSE
!   -------
!  This subroutine takes a sparse matrix in coordinate format and
!  converts it into the Ellpack-Itpack storage.
!
!  Example:
!  -------
!       (   11   0   13    0     0     0  )
!       |   21  22    0   24     0     0  |
!       |    0  32   33    0    35     0  |
!   A = |    0   0   43   44     0    46  |
!       |   51   0    0   54    55     0  |
!       (   61  62    0    0    65    66  )
!
!   Coordinate storage scheme:
!
!    A  = (11,22,33,44,55,66,13,21,24,32,35,43,46,51,54,61,62,65)
!    IA = (1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6 )
!    JA = ( 1, 2, 3, 4, 5, 6, 3, 1, 4, 2, 5, 3, 6, 1, 4, 1, 2, 5)
!
!   Ellpack-Itpack storage scheme:
!
!       (   11  13    0    0   )          (   1   3   *    *  )
!       |   22  21   24    0   |          |   2   1   4    *  |
!  AC = |   33  32   35    0   |    JAC = |   3   2   5    *  |
!       |   44  43   46    0   |          |   4   3   6    *  |
!       |   55  51   54    0   |          |   5   1   4    *  |
!       (   66  61   62   65   )          (   6   1   2    5  )
!
!   Note: * means that you can store values from 1 to 6 (1 to n, where
!         n is the order of the matrix) in that position in the array.
!
!   Contributed by:
!   --------------- 
!   Ernest E. Rothman
!   Cornell Thoery Center/Cornell National Supercomputer Facility
!   e-mail address: BITNET:   EER@CORNELLF.BITNET
!                   INTERNET: eer@cornellf.tn.cornell.edu
!   
!   checked and modified  04/13/90 Y.Saad.
!
!   REFERENCES
!   ----------
!   Kincaid, D. R.; Oppe, T. C.; Respess, J. R.; Young, D. M. 1984.
!   ITPACKV 2C User's Guide, CNA-191. Center for Numerical Analysis,
!   University of Texas at Austin.
!
!   "Engineering and Scientific Subroutine Library; Guide and
!   Reference; Release 3 (SC23-0184-3). Pp. 79-86.
!
!-----------------------------------------------------------------------
!
!   INPUT PARAMETERS
!   ----------------
!  N       - Integer. The size of the square matrix.
!
!  NNZ     - Integer. Must be greater than or equal to the number of
!            nonzero elements in the sparse matrix. Dimension of A, IA 
!            and JA.
!
!  NCA     - Integer. First dimension of output arrays ca and jac.
!
!  A(NNZ)  - Real array. (Double precision)
!            Stored entries of the sparse matrix A.
!            NNZ is the number of nonzeros.
!
!  IA(NNZ) - Integer array.
!            Pointers to specify rows for the stored nonzero entries
!            in A.
!
!  JA(NNZ) - Integer array.
!            Pointers to specify columns for the stored nonzero
!            entries in A.
!
!  NER     - Integer. Must be set greater than or equal to the maximum
!            number of nonzeros in any row of the sparse matrix.
!
!  OUTPUT PARAMETERS
!  -----------------
!  AC(NAC,*)  - Real array. (Double precision)
!               Stored entries of the sparse matrix A in compressed
!               storage mode.
!
!  JAC(NAC,*) - Integer array.
!               Contains the column numbers of the sparse matrix
!               elements stored in the corresponding positions in
!               array AC.
!
!  NCMAX   -  Integer. Equals the maximum number of nonzeros in any
!             row of the sparse matrix.
!
!  IERR    - Error parameter is returned as zero on successful
!             execution of the subroutin<e.
!             Error diagnostics are given by means of positive values
!             of this parameter as follows:
!
!             IERR = -1   -  NER is too small and should be set equal
!                            to NCMAX. The array AC may not be large
!                            enough to accomodate all the non-zeros of
!                            of the sparse matrix.
!             IERR =  1   -  The array AC has a zero column. (Warning) 
!             IERR =  2   -  The array AC has a zero row.    (Warning)
!
!---------------------------------------------------------------------
      real*8 a(nnz), ac(nac,ner)
      integer ja(nnz), ia(nnz), jac(nac,ner), ierr, ncmax, icount
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   Initial error parameter to zero:
!
      ierr = 0
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   Initial output arrays to zero:
!
      do 4 in = 1,ner
         do 4 innz =1,n
            jac(innz,in) = n
            ac(innz,in) = 0.0d0
 4    continue
!     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   Assign nonzero elements of the sparse matrix (stored in the one
!   dimensional array A to the two dimensional array AC.
!   Also, assign the correct values with information about their
!   column indices to the two dimensional array KA. And at the same
!   time count the number of nonzeros in each row so that the
!   parameter NCMAX equals the maximum number of nonzeros in any row
!   of the sparse matrix.
!
      ncmax = 1
      do 10 is = 1,n
         k = 0
         do 30 ii = 1,nnz
            if(ia(ii).eq.is)then
               k = k + 1
               if (k .le. ner) then
                  ac(is,k) = a(ii)
                  jac(is,k) = ja(ii)
               endif 
            endif
 30      continue
         if (k.ge.ncmax) ncmax = k
 10   continue
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     
!     Perform some simple error checks:
!     
!heck maximum number of nonzeros in each row:
      if (ncmax.eq.ner) ierr = 0
      if (ncmax.gt.ner) then
         ierr = -1
         return
      endif
!     
!heck if there are any zero columns in AC:
!     
      do 45 in = 1,ncmax
         icount = 0
         do 44 inn =1,n
            if (ac(inn,in).ne.0.0d0) icount = 1
 44      continue
         if (icount.eq.0) then
            ierr = 1
            return
         endif
 45   continue
!     
!heck if there are any zero rows in AC:
!     
      do 55 inn = 1,n
         icount = 0
         do 54 in =1,ncmax
            if (ac(inn,in).ne.0.0d0) icount = 1
 54      continue
         if (icount.eq.0) then
            ierr = 2
            return
         endif
 55   continue
      return
!------------- end of xcooell ------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine csruss (nrow,a,ja,ia,diag,al,jal,ial,au,jau,iau) 
      real*8 a(*),al(*),diag(*),au(*) 
      integer nrow,ja(*),ia(nrow+1),jal(*),ial(nrow+1),jau(*),
     *     iau(nrow+1)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Unsymmetric Sparse Skyline format
!----------------------------------------------------------------------- 
! this subroutine converts a matrix stored in csr format into a nonsym. 
! sparse skyline format. This latter format does not assume
! that the matrix has a symmetric pattern and consists of the following 
! * the diagonal of A stored separately in diag(*);
! * The strict lower part of A is stored  in CSR format in al,jal,ial 
! * The strict upper part is stored in CSC format in au,jau,iau.
!----------------------------------------------------------------------- 
! On entry
!---------
! nrow  = dimension of the matrix a.
! a     = real array containing the nonzero values of the matrix 
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
! 
! On return
!----------
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in CSR format storing the strict lower 
!              trangular part of A.
! au,jau,iau = matrix in CSC format storing the strict upper
!              triangular part of A. 
!----------------------------------------------------------------------- 
      integer i, j, k, kl, ku 
!
! determine U's data structure first
! 
      do 1 i=1,nrow+1
         iau(i) = 0
 1    continue
      do 3 i=1, nrow
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)
            if (j .gt. i) iau(j+1) = iau(j+1)+1
 2       continue 
 3    continue
!
!     compute pointers from lengths
!
      iau(1) = 1
      do 4 i=1,nrow
         iau(i+1) = iau(i)+iau(i+1)
         ial(i+1) = ial(i)+ial(i+1)
 4    continue
!
!     now do the extractions. scan all rows.
!
      kl = 1
      ial(1) = kl
      do  7 i=1, nrow
!
!     scan all elements in a row
! 
         do 71 k = ia(i), ia(i+1)-1
            j = ja(k) 
!
!     if in upper part, store in row j (of transp(U) )
!     
            if (j  .gt. i) then
               ku = iau(j) 
               au(ku) = a(k)
               jau(ku) = i
               iau(j) = ku+1
            elseif (j  .eq. i) then
               diag(i) = a(k) 
            elseif (j .lt. i) then
               al(kl) = a(k)
               jal(kl) = j
               kl = kl+1
            endif
 71      continue
         ial(i+1) = kl 
 7    continue
!
! readjust iau
!
      do 8 i=nrow,1,-1
         iau(i+1) = iau(i)
 8    continue
      iau(1) = 1
!--------------- end-of-csruss ----------------------------------------- 
!-----------------------------------------------------------------------
      end 
!-----------------------------------------------------------------------
      subroutine usscsr (nrow,a,ja,ia,diag,al,jal,ial,au,jau,iau) 
      real*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1),jau(*),iau(nrow+1)
!-----------------------------------------------------------------------
! Unsymmetric Sparse Skyline   format   to Compressed Sparse Row 
!----------------------------------------------------------------------- 
! this subroutine converts a matrix stored in nonsymmetric sparse
! skyline format into csr format. The sparse skyline format is 
! described in routine csruss. 
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
! On entry
!-----------------------------------------------------------------------
! nrow  = dimension of the matrix a.
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in CSR format storing the strict lower 
!              trangular part of A.
! au,jau,iau = matrix in CSC format storing the strict upper
!              trangular part of A.
! On return
! --------- 
! a     = real array containing the nonzero values of the matrix 
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
! 
!-----------------------------------------------------------------------
!
! count elements in lower part + diagonal 
! 
      do 1 i=1, nrow
         ia(i+1) = ial(i+1)-ial(i)+1
 1    continue
!
! count elements in upper part
! 
      do 3 i=1, nrow
         do 2 k=iau(i), iau(i+1)-1 
            j = jau(k)
            ia(j+1) = ia(j+1)+1
 2       continue 
 3    continue
!---------- compute pointers from lengths ------------------------------
      ia(1) = 1
      do 4 i=1,nrow
         ia(i+1) = ia(i)+ia(i+1)
 4    continue
!
! copy lower part + diagonal 
! 
      do 6 i=1, nrow
         ka = ia(i) 
         do 5 k=ial(i), ial(i+1)-1
            a(ka) = al(k) 
            ja(ka) = jal(k) 
            ka = ka+1
 5       continue
         a(ka) = diag(i) 
         ja(ka) = i
         ia(i) = ka+1
 6    continue
!     
!     copy upper part
!     
      do 8 i=1, nrow
         do 7 k=iau(i), iau(i+1)-1
!
! row number
!
            jak = jau(k) 
!
! where element goes
!
            ka = ia(jak) 
            a(ka) = au(k) 
            ja(ka) = i
            ia(jak) = ka+1
 7       continue
 8    continue
!
! readjust ia
!
      do 9 i=nrow,1,-1
         ia(i+1) = ia(i)
 9    continue
      ia(1) = 1
!----------end-of-usscsr------------------------------------------------
      end 
!----------------------------------------------------------------------- 
      subroutine csrsss (nrow,a,ja,ia,sorted,diag,al,jal,ial,au)
      real*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1)
      logical sorted 
!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Symmetric Sparse Skyline   format 
!----------------------------------------------------------------------- 
! this subroutine converts a matrix stored in csr format into the 
! Symmetric sparse skyline   format. This latter format assumes that 
! that the matrix has a symmetric pattern. It consists of the following 
! * the diagonal of A stored separately in diag(*);
! * The strict lower part of A is stored  in csr format in al,jal,ial 
! * The values only of strict upper part as stored in csc format in au. 
!----------------------------------------------------------------------- 
! On entry
!-----------
! nrow  = dimension of the matrix a.
! a     = real array containing the nonzero values of the matrix 
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
! sorted= a logical indicating whether or not the elements in a,ja,ia
!         are sorted. 
! 
! On return
! --------- 
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in csr format storing the strict lower 
!              trangular part of A.
! au    = values of the strict upper trangular part of A, column wise.
!----------------------------------------------------------------------- 
! 
!     extract lower part and diagonal.
!
      kl = 1
      ial(1) = kl
      do  7 i=1, nrow
!
! scan all elements in a row
! 
         do 71 k = ia(i), ia(i+1)-1
            jak = ja(k) 
            if (jak  .eq. i) then
               diag(i) = a(k) 
            elseif (jak .lt. i) then
               al(kl) = a(k)
               jal(kl) = jak
               kl = kl+1
            endif
 71      continue
         ial(i+1) = kl 
 7    continue
!
! sort if not sorted
! 
      if (.not. sorted) then
!%%%%%---- incompatible arg list! 
         call csort (nrow, al, jal, ial, au, .true.) 
      endif
!
! copy u
! 
      do  8 i=1, nrow
!
! scan all elements in a row
! 
         do 81 k = ia(i), ia(i+1)-1
            jak = ja(k) 
            if (jak  .gt. i) then
               ku = ial(jak) 
               au(ku) = a(k)
               ial(jak) = ku+1
            endif
 81      continue
 8    continue
!   
! readjust ial
!
      do 9 i=nrow,1,-1
         ial(i+1) = ial(i)
 9    continue
      ial(1) = 1
!--------------- end-of-csrsss ----------------------------------------- 
!-----------------------------------------------------------------------
      end 
!
      subroutine ssscsr (nrow,a,ja,ia,diag,al,jal,ial,au) 
      real*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1) 
!-----------------------------------------------------------------------
! Unsymmetric Sparse Skyline   format   to Compressed Sparse Row 
!----------------------------------------------------------------------- 
! this subroutine converts a matrix stored in nonsymmetric sparse 
! skyline format into csr format. The sparse skyline format is 
! described in routine csruss. 
!----------------------------------------------------------------------- 
! On entry
!--------- 
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in csr format storing the strict lower 
!              trangular part of A.
! au    = values of strict upper part. 
!
! On return
! --------- 
! nrow  = dimension of the matrix a.
! a     = real array containing the nonzero values of the matrix 
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
! 
!-----------------------------------------------------------------------
!
! count elements in lower part + diagonal 
! 
      do 1 i=1, nrow
         ia(i+1) = ial(i+1)-ial(i)+1
 1    continue
!
! count elements in upper part
! 
      do 3 i=1, nrow
         do 2 k=ial(i), ial(i+1)-1 
            j = jal(k)
            ia(j+1) = ia(j+1)+1
 2       continue 
 3    continue
!---------- compute pointers from lengths ------------------------------
      ia(1) = 1
      do 4 i=1,nrow
         ia(i+1) = ia(i)+ia(i+1)
 4    continue
!
! copy lower part + diagonal 
! 
      do 6 i=1, nrow
         ka = ia(i) 
         do 5 k=ial(i), ial(i+1)-1
            a(ka) = al(k) 
            ja(ka) = jal(k) 
            ka = ka+1
 5       continue
         a(ka) = diag(i) 
         ia(i) = ka+1
 6    continue
!     
!     copy upper part
!     
      do 8 i=1, nrow
         do 7 k=ial(i), ial(i+1)-1
!
! row number
!
            jak = jal(k) 
!
! where element goes
!
            ka = ia(jak) 
            a(ka) = au(k) 
            ja(ka) = i
            ia(jak) = ka+1
 7       continue
 8    continue
!
! readjust ia
!
      do 9 i=nrow,1,-1
         ia(i+1) = ia(i)
 9    continue
      ia(1) = 1
!----------end-of-ssscsr------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine csrvbr(n,ia,ja,a,nr,nc,kvstr,kvstc,ib,jb,kb,
     &     b, job, iwk, nkmax, nzmax, ierr )
!-----------------------------------------------------------------------
      integer n, ia(n+1), ja(*), nr, nc, ib(*), jb(nkmax-1), kb(nkmax)
      integer kvstr(*), kvstc(*), job, iwk(*), nkmax, nzmax, ierr
      real*8  a(*), b(nzmax)
!-----------------------------------------------------------------------
!     Converts compressed sparse row to variable block row format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of matrix rows
!     ia,ja,a = input matrix in CSR format
!
!     job     = job indicator.
!               If job=0, kvstr and kvstc are used as supplied.
!               If job=1, kvstr and kvstc are determined by the code.
!               If job=2, a conformal row/col partitioning is found and
!               returned in both kvstr and kvstc.  In the latter two cases,
!               an optimized algorithm can be used to perform the
!               conversion because all blocks are full.
!
!     nkmax   = size of supplied jb and kb arrays
!     nzmax   = size of supplied b array
!
!     If job=0 then the following are input:
!     nr,nc   = matrix block row and block column dimension
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column.
!               (kvstr and kvstc may be the same array)
!
!     On return:
!---------------
!
!     ib,jb,kb,b = output matrix in VBR format
!
!     ierr    = error message
!               ierr = 0 means normal return
!               ierr = 1 out of space in jb and/or kb arrays
!               ierr = 2 out of space in b array
!               ierr = 3 nonsquare matrix used with job=2
!
!     If job=1,2 then the following are output:
!     nr,nc   = matrix block row and block column dimension
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column
!               If job=2, then kvstr and kvstc contain the same info.
!
!     Work space:
!----------------
!     iwk(1:ncol) = inverse kvstc array.  If job=1,2 then we also need:
!     iwk(ncol+1:ncol+nr) = used to help determine sparsity of each block row.
!     The workspace is not assumed to be initialized to zero, nor is it
!     left that way.
!
!     Algorithms:
!----------------
!     There are two conversion codes in this routine.  The first assumes
!     that all blocks are full (there is a nonzero in the CSR data
!     structure for each entry in the block), and is used if the routine
!     determines the block partitioning itself.  The second code makes
!     no assumptions about the block partitioning, and is used if the
!     caller provides the partitioning.  The second code is much less
!     efficient than the first code.
!
!     In the first code, the CSR data structure is traversed sequentially
!     and entries are placed into the VBR data structure with stride
!     equal to the row dimension of the block row.  The columns of the
!     CSR data structure are sorted first if necessary.
!
!     In the second code, the block sparsity pattern is first determined.
!     This is done by traversing the CSR data structure and using an
!     implied linked list to determine which blocks are nonzero.  Then
!     the VBR data structure is filled by mapping each individual entry
!     in the CSR data structure into the VBR data structure.  The columns
!     of the CSR data structure are sorted first if necessary.
!
!-----------------------------------------------------------------------
!     Local variables:
!---------------------
      integer ncol, nb, neqr, numc, a0, b0, b1, k0, i, ii, j, jj, jnew
      logical sorted
!
!     ncol = number of scalar columns in matrix
!     nb = number of blocks in conformal row/col partitioning
!     neqr = number of rows in block row
!     numc = number of nonzero columns in row
!     a0 = index for entries in CSR a array
!     b0 = index for entries in VBR b array
!     b1 = temp
!     k0 = index for entries in VBR kb array
!     i  = loop index for block rows
!     ii = loop index for scalar rows in block row
!     j  = loop index for block columns
!     jj = loop index for scalar columns in block column
!     jnew = block column number
!     sorted = used to indicate if matrix already sorted by columns
!
!-----------------------------------------------------------------------
      ierr = 0
!-----sort matrix by column indices
      call csorted(n, ia, ja, sorted)
      if (.not. sorted) then
         call csort (n, a, ja, ia, b, .true.)
      endif
      if (job .eq. 1 .or. job .eq. 2) then
!--------need to zero workspace; first find ncol
         ncol = 0
         do i = 2, n
            ncol = max0(ncol, ja(ia(i)-1))
         enddo
         do i = 1, ncol
            iwk(i) = 0
         enddo
         call csrkvstr(n, ia, ja, nr, kvstr)
         call csrkvstc(n, ia, ja, nc, kvstc, iwk)
      endif
!-----check if want conformal partitioning
      if (job .eq. 2) then
         if (kvstr(nr+1) .ne. kvstc(nc+1)) then
            ierr = 3
            return
         endif
!        use iwk temporarily
         call kvstmerge(nr, kvstr, nc, kvstc, nb, iwk)
         nr = nb
         nc = nb
         do i = 1, nb+1
            kvstr(i) = iwk(i)
            kvstc(i) = iwk(i)
         enddo
      endif
!-----------------------------------------------------------------------
!     inverse kvst (scalar col number) = block col number
!     stored in iwk(1:n)
!-----------------------------------------------------------------------
      do i = 1, nc
         do j = kvstc(i), kvstc(i+1)-1
            iwk(j) = i
         enddo
      enddo
      ncol = kvstc(nc+1)-1
!-----jump to conversion routine
      if (job .eq. 0) goto 400
!-----------------------------------------------------------------------
!     Fast conversion for computed block partitioning
!-----------------------------------------------------------------------
      a0 = 1
      b0 = 1
      k0 = 1
      kb(1) = 1
!-----loop on block rows
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ia(kvstr(i)+1) - ia(kvstr(i))
         ib(i) = k0
!--------loop on first row in block row to determine block sparsity
         j = 0
         do jj = ia(kvstr(i)), ia(kvstr(i)+1)-1
            jnew = iwk(ja(jj))
            if (jnew .ne. j) then
!--------------check there is enough space in kb and jb arrays
               if (k0+1 .gt. nkmax) then
                  ierr = 1
                  write (*,*) 'csrvbr: no space in kb for block row ', i
                  return
               endif
!--------------set entries for this block
               j = jnew
               b0 = b0 + neqr * (kvstc(j+1) - kvstc(j))
               kb(k0+1) = b0
               jb(k0) = j
               k0 = k0 + 1
            endif
         enddo
!--------loop on scalar rows in block row
         do ii = 0, neqr-1
            b1 = kb(ib(i))+ii
!-----------loop on elements in a scalar row
            do jj = 1, numc
!--------------check there is enough space in b array
               if (b1 .gt. nzmax) then
                  ierr = 2
                  write (*,*) 'csrvbr: no space in b for block row ', i
                  return
               endif
               b(b1) = a(a0)
               b1 = b1 + neqr
               a0 = a0 + 1
            enddo
         enddo
      enddo
      ib(nr+1) = k0
      return
!-----------------------------------------------------------------------
!     Conversion for user supplied block partitioning
!-----------------------------------------------------------------------
 400  continue
!-----initialize workspace for sparsity indicator
      do i = ncol+1, ncol+nc
         iwk(i) = 0
      enddo
      k0 = 1
      kb(1) = 1
!-----find sparsity of block rows
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ia(kvstr(i)+1) - ia(kvstr(i))
         ib(i) = k0
!--------loop on all the elements in the block row to determine block sparsity
         do jj = ia(kvstr(i)), ia(kvstr(i+1))-1
            iwk(iwk(ja(jj))+ncol) = 1
         enddo
!--------use sparsity to set jb and kb arrays
         do j = 1, nc
            if (iwk(j+ncol) .ne. 0) then
!--------------check there is enough space in kb and jb arrays
               if (k0+1 .gt. nkmax) then
                  ierr = 1
                  write (*,*) 'csrvbr: no space in kb for block row ', i
                  return
               endif
               kb(k0+1) = kb(k0) + neqr * (kvstc(j+1) - kvstc(j))
               jb(k0) = j
               k0 = k0 + 1
               iwk(j+ncol) = 0
            endif
         enddo
      enddo
      ib(nr+1) = k0
!-----Fill b with entries from a by traversing VBR data structure.
      a0 = 1
!-----loop on block rows
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
!--------loop on scalar rows in block row
         do ii = 0, neqr-1
            b0 = kb(ib(i)) + ii
!-----------loop on block columns
            do j = ib(i), ib(i+1)-1
!--------------loop on scalar columns within block column
               do jj = kvstc(jb(j)), kvstc(jb(j)+1)-1
!-----------------check there is enough space in b array
                  if (b0 .gt. nzmax) then
                     ierr = 2
                     write (*,*)'csrvbr: no space in b for blk row',i
                     return
                  endif
                  if (a0 .ge. ia(kvstr(i)+ii+1)) then
                     b(b0) = 0.d0
                  else
                     if (jj .eq. ja(a0)) then
                        b(b0) = a(a0)
                        a0 = a0 + 1
                     else
                        b(b0) = 0.d0
                     endif
                  endif
                  b0 = b0 + neqr
!--------------endloop on scalar columns
               enddo
!-----------endloop on block columns
            enddo
 2020       continue
         enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
!----------------------------end-of-csrvbr------------------------------
!----------------------------------------------------------------------c
      subroutine vbrcsr(ia, ja, a, nr, kvstr, kvstc, ib, jb, kb,
     &   b, nzmax, ierr)
!-----------------------------------------------------------------------
      integer ia(*), ja(*), nr, ib(nr+1), jb(*), kb(*)
      integer kvstr(nr+1), kvstc(*), nzmax, ierr
      real*8  a(*), b(nzmax)
!-----------------------------------------------------------------------
!     Converts variable block row to compressed sparse row format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     nr      = number of block rows
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column
!     ib,jb,kb,b = input matrix in VBR format
!     nzmax   = size of supplied ja and a arrays
!
!     On return:
!---------------
!     ia,ja,a = output matrix in CSR format
!
!     ierr    = error message
!               ierr = 0 means normal return
!               ierr = negative row number when out of space in
!                      ja and a arrays
!
!     Work space:
!----------------
!     None
!
!     Algorithm:
!---------------
!     The VBR data structure is traversed in the order that is required
!     to fill the CSR data structure.  In a given block row, consecutive
!     entries in the CSR data structure are entries in the VBR data
!     structure with stride equal to the row dimension of the block.
!     The VBR data structure is assumed to be sorted by block columns.
!
!-----------------------------------------------------------------------
!     Local variables:
!---------------------
      integer neqr, numc, a0, b0, i, ii, j, jj
!
!     neqr = number of rows in block row
!     numc = number of nonzero columns in row
!     a0 = index for entries in CSR a array
!     b0 = index for entries in VBR b array
!     i  = loop index for block rows
!     ii = loop index for scalar rows in block row
!     j  = loop index for block columns
!     jj = loop index for scalar columns in block column
!
!-----------------------------------------------------------------------
      ierr = 0
      a0 = 1
      b0 = 1
!-----loop on block rows
      do i = 1, nr
!--------set num of rows in block row, and num of nonzero cols in row
         neqr = kvstr(i+1) - kvstr(i)
         numc = ( kb(ib(i+1)) - kb(ib(i)) ) / neqr
!--------construct ja for a scalar row
         do j = ib(i), ib(i+1)-1
            do jj = kvstc(jb(j)), kvstc(jb(j)+1)-1
               ja(a0) = jj
               a0 = a0 + 1
            enddo
         enddo
!--------construct neqr-1 additional copies of ja for the block row
         do ii = 1, neqr-1
            do j = 1, numc
               ja(a0) = ja(a0-numc)
               a0 = a0 + 1
            enddo
         enddo
!--------reset a0 back to beginning of block row
         a0 = kb(ib(i))
!--------loop on scalar rows in block row
         do ii = 0, neqr-1
            ia(kvstr(i)+ii) = a0
            b0 = kb(ib(i)) + ii
!-----------loop on elements in a scalar row
            do jj = 1, numc
!--------------check there is enough space in a array
               if (a0 .gt. nzmax) then
                  ierr = -(kvstr(i)+ii)
                  write (*,*) 'vbrcsr: no space for row ', -ierr
                  return
               endif
               a(a0) = b(b0)
               a0 = a0 + 1
               b0 = b0 + neqr
            enddo
         enddo
!-----endloop on block rows
      enddo
      ia(kvstr(nr+1)) = a0
      return
      end
!-----------------------------------------------------------------------
!---------------------------end-of-vbrcsr-------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine csorted(n, ia, ja, sorted)
!-----------------------------------------------------------------------
      integer n, ia(n+1), ja(*)
      logical sorted
!-----------------------------------------------------------------------
!     Checks if matrix in CSR format is sorted by columns.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of rows in matrix
!     ia, ja  = sparsity structure of matrix in CSR format
!
!     On return:
!---------------
!     sorted  = indicates if matrix is sorted by columns
!
!-----------------------------------------------------------------------
!-----local variables
      integer i,j
!---------------------------------
      do i = 1, n
         do j = ia(i)+1, ia(i+1)-1
            if (ja(j-1) .ge. ja(j)) then
               sorted = .false.
               return
            endif
         enddo
      enddo
      sorted = .true.
      return
      end
!-----------------------------------------------------------------------
!------------------------end-of-csorted---------------------------------
