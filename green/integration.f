ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    This file includes two subroutines:
c       subroutine : DWIM(o,kmin,kmax,jj,dkrec,sums)
c                    SAIM(o,kmin,kmax,r0,cc,sums) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccc Discrete Wavenumber Integration Method cccccccccccccccccccc
       	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine DWIM_S(o,kmin,kmax,jj,dkrec,sums,nforce)

      include    "green_com.inc"

	integer    jj,flag0or1(12),nforce
      real*8     j01,j11,j21,j0d1,j1d1,j2d1
      real*8     j02,j12,j22,j0d2,j1d2,j2d2
	real*8     kmin,kmax,dkrec,restdk
	complex*16 fun1(7),fun2(7),integ1(3),integ2(3),sums(3)
      common     /flag0or1/ flag0or1

	jj=0
	kn=kmin
	call funval_S(kn,o,fun1,nforce)
	call bessj01d(r0*kn,j01,j11,j0d1,j1d1)
	do i=1,7
         if(flag0or1(i).eq.1) then
	      fun1(i)=0.
         end if
	end do
	integ1(1)=fun1(1)*j01+fun1(2)*j0d1+fun1(3)*j1d1
	integ1(2)=fun1(4)*j01+fun1(5)*j1d1
	integ1(3)=fun1(6)*j01+fun1(7)*j11
	do i=1,3
         sums(i)=integ1(i)*dkrec
	end do
!222   kn=kn+dkrec
!      if(kn.ge.kmax) goto 111
      do while(kn.lt.kmax)
	  kn=kn+dkrec
	  jj=jj+1
		call funval_S(kn,o,fun2,nforce)
        call bessj01d(r0*kn,j02,j12,j0d2,j1d2)
!	do i=1,7
!         if(flag0or1(i).eq.1) then
!	      fun2(i)=0.
!         end if
!	end do
	  integ2(1)=fun2(1)*j02+fun2(2)*j0d2+fun2(3)*j1d2
	  integ2(2)=fun2(4)*j02+fun2(5)*j1d2
	  integ2(3)=fun2(6)*j02+fun2(7)*j12
	  do i=1,3
	      sums(i)=sums(i)+(integ1(i)+integ2(i))/2.*dkrec
	  end do
	  do i=1,3
	      integ1(i)=integ2(i)
	  end do
!	goto 222
!111   restdk=dkrec-(kn-kmax)
      end do
      restdk=dkrec-(kn-kmax)
      kn=kmax
	jj=jj+1
	call funval_S(kn,o,fun2,nforce)
      call bessj01d(r0*kn,j02,j12,j0d2,j1d2)
!	do i=1,7
!         if(flag0or1(i).eq.1) then
!	      fun2(i)=0.
!         end if
!	end do
	integ2(1)=fun2(1)*j02+fun2(2)*j0d2+fun2(3)*j1d2
	integ2(2)=fun2(4)*j02+fun2(5)*j1d2
	integ2(3)=fun2(6)*j02+fun2(7)*j12
	do i=1,3
	   sums(i)=sums(i)+(integ1(i)+integ2(i))/2.*restdk
	end do

	return
	end 



