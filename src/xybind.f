c     system("R CMD SHLIB ../src/xybind.f")

      subroutine cbindf(xncol,nrow, a,ia,ja, b,ib,jb,
     &     c,ic,jc)

      implicit none
      integer xncol, nrow
      integer ia(*), ja(*), ib(*), jb(*), ic(*), jc(*)

      double precision  a(*), b(*), c(*)

      integer j,j1,i,k

      k=1
      do j = 1,nrow
         jc(j)=ja(j)+jb(j)-1

         j1=j+1
         if (ja(j) .lt. ja(j1)) then
            do i=ja(j),ja(j1)-1
               c(k)=a(i)
               ic(k)=ia(i)
               k=k+1
c               if (k.gt.clen)  return
            enddo
         endif
         if (jb(j) .lt. jb(j1)) then
            do i=jb(j),jb(j1)-1
               c(k)=b(i)
               ic(k)=ib(i)+xncol
               k=k+1
c               if (k.gt.clen)  return
            enddo
         endif
      enddo
      j=nrow+1
      jc(j)=ja(j)+jb(j)-1


      return
      end

