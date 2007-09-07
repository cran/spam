      subroutine backsolve(m,nsuper,nrhs,lindx,xlindx,lnz,
     &                   xlnz,xsuper,b)
c see below...
      implicit none

      integer m,nsuper,nrhs,lindx(*),xlindx(m+1),
     &        xlnz(m+1),xsuper(m+1)
      double precision lnz(*),b(m,nrhs)
      integer j
      do j = 1,nrhs
         call blkslb(nsuper,xsuper,xlindx,lindx,xlnz,lnz,b(1,j))
      enddo
      return
      end


      subroutine forwardsolve(m,nsuper,nrhs,lindx,xlindx,
     &  lnz,xlnz,xsuper,b)
c INPUT:
c     m -- the number of column in the matrix
c     lindx -- an nsub-vector of interger which contains, in 
c           column major oder, the row subscripts of the nonzero
c           entries in L in a compressed storage format
c     xlindx -- an nsuper-vector of integer of pointers for lindx
c     lnz -- First contains the non-zero entries of d; later
c            contains the entries of the Cholesky factor
c     xlnz -- column pointer for L stored in lnz
c     xsuper -- array of length m+1 containing the supernode
c               partitioning
c     b -- the rhs of the equality constraint
c OUTPUT:
c     b -- the solution
      
      implicit none

      integer m,nsuper,nrhs,lindx(*),xlindx(m+1),
     &        xlnz(m+1),xsuper(m+1)
      double precision lnz(*),b(m,nrhs)
      integer j
c
c      write(*,*)"forwardsolve"
      do j = 1,nrhs
         call blkslf(nsuper,xsuper,xlindx,lindx,xlnz,lnz,b(1,j))
      enddo
      return
      end
C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C                   Slight modification by Reinhard Furrer 
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C
C***********************************************************************
      subroutine bckslf(m,nsubmax,nsuper,nrhs,lindx,xlindx,nnzlmax,lnz,
     &                   xlnz,invp,perm,xsuper,newrhs,sol,b)
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c Sparse least squares solver via Ng-Peyton's sparse Cholesky 
c    factorization for sparse symmetric positive definite
c INPUT:
c     m -- the number of column in the design matrix X
c     nsubmax -- upper bound of the dimension of lindx
c     lindx -- an nsub-vector of interger which contains, in 
c           column major oder, the row subscripts of the nonzero
c           entries in L in a compressed storage format
c     xlindx -- an nsuper-vector of integer of pointers for lindx
c     nnzlmax -- the upper bound of the non-zero entries in
c                L stored in lnz, including the diagonal entries
c     lnz -- First contains the non-zero entries of d; later
c            contains the entries of the Cholesky factor
c     xlnz -- column pointer for L stored in lnz
c     invp -- an m-vector of integer of inverse permutation
c             vector
c     perm -- an m-vector of integer of permutation vector
c     xsuper -- array of length m+1 containing the supernode
c               partitioning
c     newrhs -- extra work vector for right-hand side and
c               solution
c     sol -- the least squares solution
c     b -- an m-vector, usualy the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c OUTPUT:
c     y -- an m-vector of least squares solution
c WORK ARRAYS:
c     b -- an m-vector, usually the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      integer nnzlmax,nsubmax,nsuper,nrhs,lindx(nsubmax),xlindx(m+1),
     &        invp(m),perm(m),xlnz(m+1),
     &        xsuper(m+1)
      double precision lnz(nnzlmax),b(m,nrhs),newrhs(m),sol(m,nrhs)
c Call blkslv: Numerical solution
c
c      write(*,*)"forwardsolve"
      do j = 1,nrhs
         do i = 1,m
            newrhs(i) = b(perm(i),j)
c            newrhs(i) = b(i,j)
         enddo
         call blkslf(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
         do i = 1,m
c            sol(i,j) = newrhs(i)
            sol(i,j) = newrhs(invp(i))
         enddo
      enddo
  100 continue
      return
      end
C***********************************************************************
      subroutine bckslb(m,nsubmax,nsuper,nrhs,lindx,xlindx,nnzlmax,lnz,
     &                   xlnz,invp,perm,xsuper,newrhs,sol,b)
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c Sparse least squares solver via Ng-Peyton's sparse Cholesky 
c    factorization for sparse symmetric positive definite
c INPUT:
c     m -- the number of column in the design matrix X
c     nsubmax -- upper bound of the dimension of lindx
c     lindx -- an nsub-vector of interger which contains, in 
c           column major oder, the row subscripts of the nonzero
c           entries in L in a compressed storage format
c     xlindx -- an nsuper-vector of integer of pointers for lindx
c     nnzlmax -- the upper bound of the non-zero entries in
c                L stored in lnz, including the diagonal entries
c     lnz -- First contains the non-zero entries of d; later
c            contains the entries of the Cholesky factor
c     xlnz -- column pointer for L stored in lnz
c     invp -- an m-vector of integer of inverse permutation
c             vector
c     perm -- an m-vector of integer of permutation vector
c     xsuper -- array of length m+1 containing the supernode
c               partitioning
c     newrhs -- extra work vector for right-hand side and
c               solution
c     sol -- the least squares solution
c     b -- an m-vector, usualy the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c OUTPUT:
c     y -- an m-vector of least squares solution
c WORK ARRAYS:
c     b -- an m-vector, usually the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      integer nnzlmax,nsubmax,nsuper,nrhs,lindx(nsubmax),xlindx(m+1),
     &        invp(m),perm(m),xlnz(m+1),
     &        xsuper(m+1)
      double precision lnz(nnzlmax),b(m,nrhs),newrhs(m),sol(m,nrhs)
c Call blkslv: Numerical solution
c
c      write(*,*)"backsolve"
      do j = 1,nrhs
         do i = 1,m
            newrhs(i) = b(perm(i),j)
c            newrhs(i) = b(i,j)
         enddo
         call blkslb(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
         do i = 1,m
            sol(i,j) = newrhs(invp(i))
c            sol(i,j) = newrhs(i)
         enddo
      enddo
  100 continue
      return
      end
C***********************************************************************
      subroutine bckslv(m,nsubmax,nsuper,nrhs,lindx,xlindx,nnzlmax,lnz,
     &                   xlnz,invp,perm,xsuper,newrhs,sol,b)
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c Sparse least squares solver via Ng-Peyton's sparse Cholesky 
c    factorization for sparse symmetric positive definite
c INPUT:
c     m -- the number of column in the design matrix X
c     nsubmax -- upper bound of the dimension of lindx
c     lindx -- an nsub-vector of interger which contains, in 
c           column major oder, the row subscripts of the nonzero
c           entries in L in a compressed storage format
c     xlindx -- an nsuper-vector of integer of pointers for lindx
c     nnzlmax -- the upper bound of the non-zero entries in
c                L stored in lnz, including the diagonal entries
c     lnz -- First contains the non-zero entries of d; later
c            contains the entries of the Cholesky factor
c     xlnz -- column pointer for L stored in lnz
c     invp -- an m-vector of integer of inverse permutation
c             vector
c     perm -- an m-vector of integer of permutation vector
c     xsuper -- array of length m+1 containing the supernode
c               partitioning
c     newrhs -- extra work vector for right-hand side and
c               solution
c     sol -- the least squares solution
c     b -- an m-vector, usualy the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c OUTPUT:
c     y -- an m-vector of least squares solution
c WORK ARRAYS:
c     b -- an m-vector, usually the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      integer nnzlmax,nsubmax,nsuper,nrhs,lindx(nsubmax),xlindx(m+1),
     &        invp(m),perm(m),xlnz(m+1),
     &        xsuper(m+1)
      double precision lnz(nnzlmax),b(m,nrhs),newrhs(m),sol(m,nrhs)

c Call blkslv: Numerical solution
c
c      write(*,*)"solve"
      do j = 1,nrhs
         do i = 1,m
            newrhs(i) = b(perm(i),j)
         enddo
         call blkslv(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
         do i = 1,m
            sol(i,j) = newrhs(invp(i))
         enddo
      enddo
  100 continue
      return
      end
