C     closestdistXY,  distance between x and x or between x and y
C
C
C     We have four distances implemented:
C       c("euclidean", "maximum", "minkowski", "greatcircle")
C     
c     In case we need the distance matrix between x and x, then the
c     following parameters are used as well:
c     diag=0 include diagonal zero value, diag=1 no diagonal 
c     if part=-1, lower tri, part=0 the entire matrix
c        part= 1, upper tri only.
c     only values between eps and eta are considered.
c     p power for minkowski



      double precision function euclid(x,y,p)
      implicit none
      double precision x,y,p
      euclid=(x-y)**2
      return
      end



      double precision function minkowski(x,y,p)
      implicit none
      double precision x,y,p
      minkowski=abs(x-y)**p
      return
      end
      


      subroutine closestdist( ncol, x,nrowx, y, nrowy,
     &    diag, part, p, method, 
     &    eps, eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none
      
      
      double precision euclid, minkowski
      external euclid, minkowski

      integer ncol,nrowx, nrowy, nnz, method,  diag, part,  iflag
      integer colindices(nnz), rowpointers(nrowx+1)
      double precision p, x(nrowx,ncol),y(nrowy,ncol) 
      double precision eps, eta, entries(nnz)
      



      if (method.eq.1) then
         p=2.0
         call closestEdistXY( ncol, x,nrowx, y, nrowy,
     &        diag, part, p, euclid,
     &        eps, eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      if (method.eq.2) then
         p=1.0
         call closestMAXdistXY( ncol, x,nrowx, y, nrowy,
     &        diag, part, 
     &        eps, eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      if (method.eq.3) then
         call closestEdistXY( ncol, x,nrowx, y, nrowy,
     &        diag, part, p, minkowski,
     &        eps, eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      if (method.eq.4) then
         call closestGCdistXY( x,nrowx, y, nrowy,
     &        diag, part, p,
     &        eps, eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      return
      end


      
      subroutine closestEdistXY( ncol, x,xnrow, y, ynrow,
     &    diag, part, p, distfcn,
     &    eps, eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none
      
      double precision distfcn
      external distfcn

      integer ncol,xnrow, ynrow, nnz,  diag, part,iflag
      integer colindices(nnz),rowpointers(xnrow+1)
      double precision p,x(xnrow,ncol), y(ynrow,ncol)
      double precision eps, eta, entries(nnz)
      
c     local variables
      integer jja, i,j,k, ifrom,ito, jfrom, jto
      double precision epsp,etap, tmp,pinv
      

      epsp=eps**p
      etap=eta**p
      pinv=1/p
      
      jja=1
       
      rowpointers(1)=1
      jfrom = 1
      jto = ynrow
      
      
c cycle over all rows of x (independent of part and diag!)       
      
      do i= 1,xnrow
         
         if (part .lt. 0) then
            jto = i
         endif
         if (part .gt. 0) then
            jfrom = i
         endif
         
         do 10 j = jfrom,jto
            
c Start calculating the distance (until delta is exceeded)            
            tmp = 0.0
            do  k = 1, ncol
               tmp = tmp + distfcn(x(i,k),y(j,k),p)
               if( tmp.gt.etap) goto 10
            enddo
c Delta is not exceeded. 
c If we need to do not include the diagonal (diag==0)
c    eps < tmp < delta
c If we do include the diagonal,  
c    eps < tmp < delta for i=/=j only
 
            if (diag.eq.-1) then
               if (tmp.lt.epsp) goto 10
            else
               if (i.ne.j) then
                  if (tmp.lt.epsp) goto 10
               endif
            endif

c     (i,j) has a distance between eps and eta.
            
c     in case nnz was too small, recall line to get a better estimate
            if( jja .gt. nnz) then
               iflag = i 
               goto 20
            endif
               
            colindices(jja) = j
            if (p.eq.2) then
               entries(jja) = sqrt(tmp)
            else
               if (p.eq.1) then
                  entries(jja) = tmp
               else
                  entries(jja) = tmp**pinv
               endif
            endif
            jja = jja + 1
            
            
 10      continue
         rowpointers(i+1)=jja
      enddo
      
      
      if (part.gt.0) then
         rowpointers(xnrow+1)=jja
      endif
      nnz=jja-1
 20   continue
      
      return
      end
      
      subroutine closestMAXdistXY( ncol, x,xnrow, y, ynrow,
     &    diag, part, 
     &    eps, eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none
      
      integer ncol,xnrow, ynrow, nnz,  diag, part,iflag
      integer colindices(nnz),rowpointers(xnrow+1)
      double precision x(xnrow,ncol), y(ynrow,ncol)
      double precision eps, eta, entries(nnz)
      
c     local variables
      integer jja, i,j,k, ifrom,ito, jfrom, jto
      double precision  tmp
      

      
      jja=1
       
      rowpointers(1)=1
      jfrom = 1
      jto = ynrow
      
      
      do i= 1,xnrow
         
         if (part .lt. 0) then
            jto = i
         endif
         if (part .gt. 0) then
            jfrom = i
         endif
       
         do 10 j = jfrom,jto
            
c Start calculating the distance 
            tmp = 0.0
            do  k = 1, ncol
               tmp = max(tmp, abs(x(i,k)-y(j,k))) 
               if( tmp.gt.eta) goto 10
            enddo

c Delta is not exceeded. 
c If we need to do not include the diagonal (diag==0)
c    eps < tmp < delta
c If we do include the diagonal,  
c    eps < tmp < delta for i=/=j only
 
            if (diag.eq.-1) then
               if (tmp.lt.eps) goto 10
            else
               if (i.ne.j) then
                  if (tmp.lt.eps) goto 10
               endif
            endif

c     (i,j) has a distance between eps and eta.
      
                    
            
c     in case nnz was too small, recall line to get a better estimate
            if( jja .gt. nnz) then
               iflag = i 
               goto 20
            endif
               
            colindices(jja) = j
            entries(jja) = tmp
            jja = jja + 1
            
            
 10      continue
         rowpointers(i+1)=jja
      enddo
      
      
      if (part.gt.0) then
         rowpointers(xnrow+1)=jja
      endif
      nnz=jja-1
 20   continue
      
      return
      end
      


      subroutine closestGCdistXY( x,nx, y, ny,
     &    diag, part,p, 
     &    eps, eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none
      
      integer nx, ny, nnz, colindices(nnz),rowpointers(nx+1)
      integer diag, part, iflag
      double precision x(nx,2), y(ny,2), p, eps, eta, entries(nnz)
      
c     local variables
      logical equi
      integer jja, i,j,k, ifrom,ito, jfrom, jto
      double precision etap, epsp,tmp, rad, tmp1, tmp2
      double precision scy12(ny), ccy12(ny), sy2(ny)
      double precision scx12,     ccx12,     sx2
      
      parameter (rad = 0.01745329251994329)


c     Great savings if we know that x=y. This can be done by
c     multiplying diag by two.
      if (abs(diag).gt.1) then
         equi=.TRUE.
         diag=diag/2
      else 
         equi= .FALSE. 
      endif
      
      jja=1
       
      etap=cos(eta*rad)
      epsp=cos(eps*rad)
      rowpointers(1)=1
      jfrom = 1
      jto = ny
      

      DO j=1,ny
         tmp1=y(j,1)*rad
         tmp2=y(j,2)*rad
         ccy12(j)=dcos(tmp1)*dcos(tmp2)
         scy12(j)=dsin(tmp1)*dcos(tmp2)
         sy2(j)=dsin(tmp2)
      ENDDO

       
      
      do i= 1,nx
         
c     x2 is missing if equi=.TRUE. and we reuse the y stuff 
         if (equi .eqv. .TRUE.) then 
            ccx12=ccy12(i)
            scx12=scy12(i)
            sx2=sy2(i)
         else
            tmp1=x(i,1)*rad
            tmp2=x(i,2)*rad
            ccx12=dcos(tmp1)*dcos(tmp2)
            scx12=dsin(tmp1)*dcos(tmp2)
            sx2=dsin(tmp2)
         endif 

         if (part .lt. 0) then
            jto = i
         endif
         if (part .gt. 0) then
            jfrom = i
         endif
       
         do 10 j = jfrom,jto


c Start calculating the distance 
            tmp = ccx12 * ccy12(j) + scx12 * scy12(j)  + sx2*sy2(j)
            if (tmp .lt. etap) goto 10
            
            
c Delta is not exceeded. 
c If we need to do not include the diagonal (diag==0)
c    eps < tmp < delta
c If we do include the diagonal,  
c    eps < tmp < delta for i=/=j only
 
            if (diag.eq.-1) then
               if (tmp.ge.epsp) goto 10
            else
               if (i.ne.j) then
                  if (tmp.ge.epsp) goto 10
               endif
            endif

            tmp = dacos( tmp)


c     (i,j) has a distance between eps and eta.
            
c     In case nnz was too small, recall line to get a better estimate
            if( jja .gt. nnz) then
               iflag = i 
               goto 20
            endif
               
            colindices(jja) = j
            entries(jja) = tmp*p
            jja = jja + 1
            
            
 10      continue
         rowpointers(i+1)=jja
      enddo
      
      
      if (part.gt.0) then
         rowpointers(nx+1)=jja
      endif
      nnz=jja-1
 20   continue
      
      return
      end
      


 
