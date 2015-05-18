       subroutine REGREG(Tf,x0f,dxxf,dyxf,y0f,dxyf,dyyf,imaxf,jmaxf,Tt,x0t,dxxt,dyxt,y0t,dxyt,dyyt,imaxt,jmaxt)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Interpolates from a regular full field towards another
! regular full field
! Input geometry:
!   x= x0f + I dxxf + J dyxf
!   y= y0f + I dxyf + J dyyf
!
! Output geometry
!   x= x0t + I dxxt + J dyxt
!   y= y0t + I dxyt + J dyyt
!
!
! Interpolated field Tt computed from field Tf
!
!
! JMB 15/5/93
! fs 07/08/94 xt and yt computed by adding quantities at each loop step
!             and not any more by multiplying by i,j(new vars.: xbt, ybt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT NONE
      integer imaxf,jmaxf,imaxt,jmaxt,i,j,ii,jj
      real Tf(imaxf,jmaxf)
      real Tt(imaxt,jmaxt)
      real x0f,dxxf,dyxf,y0f,dxyf,dyyf
      real x0t,dxxt,dyxt,y0t,dxyt,dyyt
      real ri,rj,det,xt,yt,xi,yi,xj,yj
      real xbt,ybt
!     write(6,*) "Interpolating",x0f,dxxf,dyxf,y0f,dxyf,dyyf,imaxf,jmaxf
!     write(6,*) " to",x0t,dxxt,dyxt,y0t,dxyt,dyyt,imaxt,jmaxt

      det=dxxf*dyyf-dxyf*dyxf
      XI=dyyf/det
      YI=-dyxf/det
      XJ=-dxyf/det
      YJ=dxxf/det

      xbt=x0t
      ybt=y0t
      do j=1,jmaxt
       xbt=xbt+dyxt
       ybt=ybt+dyyt
       xt=xbt
       yt=ybt
       do i = 1,imaxt
       xt = xt+dxxt
       yt = yt+dxyt
       Ri = XI * ( xt - x0f ) + YI * ( yt - y0f )
       Rj = XJ * ( xt - x0f ) + YJ * ( yt - y0f )
       ii = Ri
       jj = Rj
       Ri = Ri - ii
       Rj = Rj - jj

       if( ii.lt.1 ) then
        if(jj.lt.1 ) jj=1
        if(jj.ge.jmaxf) jj=jmaxf
        Tt(i,j)=Tf(1,jj)
                     else
          if(ii.ge.imaxf) then
             if(jj.lt.1 ) jj=1
             if(jj.ge.jmaxf) jj=jmaxf
             Tt(i,j)=Tf(imaxf,jj)
                          else
             if(jj.lt.1 ) then
                 Tt(i,j)=Tf(ii,1)
                          else
                 if(jj.ge.jmaxf) then
                  Tt(i,j)=Tf(ii,jmaxf)
                                  else
! Interpolate...
                  Tt(i,j)= rj* ( ri * Tf(ii+1,jj+1) + ( 1 - ri )*Tf(ii,jj+1) ) + (1 -rj) * ( ri * Tf(ii+1,jj)  + ( 1 - ri )*Tf(ii,jj) )
                 endif
             endif
          endif
        endif
      end do
     end do
        end subroutine
