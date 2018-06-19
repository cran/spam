c *****************************************************************************
c
c     print function for errorflags from __aupd/__eupd; debugging resp verbose purpose
c
      subroutine errpr ( errf )
c
      implicit none
c
      integer errf
c
      if ( errf .eq. 0 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Normal exit', -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
      else if ( errf .gt. 0 ) then
c
         if ( errf .eq. 1 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Maximum number of iterations taken.', -1, 0, 0)
      call intpr(' Increase the argument nitr or ncv.', -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         end if
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Workload of the requested eigenvalues is too high.',
     &              -1, 0, 0)
      call intpr(' Increase nitr and decrease ncv.', -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
      else
c
      call intpr(' ', -1, 0, 0)
      call intpr(' ERROR in the Implicitly Restarted Arnoldi Process.',
     &           -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         if ( errf .eq. -1 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Dimension of input matrix is not positive.',
     &           -1, 0, 0)
      
      call intpr(' ', -1, 0, 0)
c
         else if ( errf .eq. -2 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Requested eigenvalues must be positive.',
     &            -1, 0, 0) 
      call intpr(' ', -1, 0, 0)
c
         else if ( errf .eq. -3 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' The number of requested eigenvalues is too high.',
     &            -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         else if ( errf .eq. -4 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' The maximum number of Arnoldi update',
     &           -1, 0, 0)
      call intpr(' iterations allowed must be greater than zero.',
     &           -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         else if ( errf .eq. -5 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' WHICH must be one of LA,SA,LM,SM or LR,SR,LI,SI.',
     &           -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         else if ( errf .eq. -14 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' The accuracy of the eigenvalues',
     &            -1, 0, 0)
      call intpr(' is not sufficent enough.',
     &            -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         else if ( errf .eq. -9999 ) then
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Could not build an Arnoldi factorization',
     &            -1, 0, 0)
      call intpr(' try to increase ncv and nitr.',
     &            -1, 0, 0)
      call intpr(' ', -1, 0, 0)
c
         else
c
      call intpr(' ', -1, 0, 0)
      call intpr(' Fatal error in fortran source, ',
     &            -1, 0, 0)
      call intpr(' The ERRORFLAG is = ', -1, errf, 4)
      call intpr(' ', -1, 0, 0)
c
         end if
c
      end if
c
      end
c
