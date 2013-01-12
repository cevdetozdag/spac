cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This file includes the following subroutines or functions:
c      subroutine : calculate (o,k,f,r0,epsi,s2,flag)
c                   calculate_sub (k,fun,r0,s2)
c                   funval (k,o,fun)
c                   PQRScalc (n,x,P,Q,R,S)
c      functions  : poly2 (x1,x2,x3,f1,f2,f3,low,high)
c                   poly4bess2 (k,f,r0,n,m)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
      subroutine  calculate_S(o,k,f,r0,epsi,s2,flag,dkmin)
	
	implicit none
      integer    i,flag,judge(12),jj,flag0or1(12),nforce
	real*8     k(7),r0,epsi,rt(3),dkmin
	complex*16 o,f(12,7),poly2,fun3(12),fun5(12)
	complex*16 s2(3),s1d(12),s2d(12),s11(12),s12(12),s21(12),s22(12)
	common     /flag0or1/ flag0or1
	
	do i=1,7
	   judge(i)=0
	end do
	k(7)=k(5)
	k(6)=k(4)
	k(4)=k(3)
	k(5)=(k(4)+k(6))/2.
	k(3)=(k(2)+k(4))/2.
	call funval_S(k(3),o,fun3,nforce)
      call funval_S(k(5),o,fun5,nforce)
	do i=1,7
	   f(i,7)=f(i,5)
	   f(i,6)=f(i,4)
	   f(i,4)=f(i,3)
	   f(i,5)=fun5(i)
	   f(i,3)=fun3(i)
	end do

	if (k(4)-k(2).lt.dkmin) then
	
	   flag=1
	   call calculate_sub_S(k,f,r0,s2)
	
	else 

	   jj=1
	   do i=1,7
	      s1d(i)= poly2(k(2),k(4),k(6),f(i,2),f(i,4),f(i,6),k(2),k(6))
	      s2d(i)= poly2(k(2),k(3),k(4),f(i,2),f(i,3),f(i,4),k(2),k(4))
     &             +poly2(k(4),k(5),k(6),f(i,4),f(i,5),f(i,6),k(4),k(6))
	      s11(i)= poly2(k(2),k(3),k(4),f(i,2),f(i,3),f(i,4),k(2),k(4))
	      s12(i)= poly2(k(2),k(4),k(6),f(i,2),f(i,4),f(i,6),k(2),k(4))
	      s21(i)= poly2(k(4),k(5),k(6),f(i,4),f(i,5),f(i,6),k(4),k(6))
	      s22(i)= poly2(k(2),k(4),k(6),f(i,2),f(i,4),f(i,6),k(4),k(6))

	      if (flag0or1(i).eq.1) then
	         judge(i)=1
	      else 
	         rt(1)=cdabs(s12(i)-s11(i))/(cdabs(s11(i))+cdabs(s12(i)))
	         rt(2)=cdabs(s22(i)-s21(i))/(cdabs(s21(i))+cdabs(s22(i)))
	         rt(3)=cdabs(s2d(i)-s1d(i))/(cdabs(s1d(i))+cdabs(s2d(i)))
               if(rt(1).lt.epsi.and.rt(2).lt.epsi.and.rt(3).lt.epsi)then
	            judge(i)=1
	         end if
	      end if
	      jj=jj*judge(i)
	   end do

	   if (jj.eq.1) then
	      flag=1
	      call calculate_sub_S(k,f,r0,s2)
	   else 
	      flag=0
	   end if

      end if

	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine calculate_sub_S(k,fun,r0,s2) !?????????????!

	implicit none
      integer    i,j,flag0or1(12)
	real*8     k(7),r0,x1(5),x2(5)
	complex*16 sums(5),poly4bess2,s2(3),fun(12,7)
	complex*16 f1_1(5),f2_1(5),f1_2(5),f2_2(5),f1_3(5),f2_3(5)
	complex*16 f1_4(5),f2_4(5),f1_5(5),f2_5(5),f1_6(5),f2_6(5)	
      complex*16 f1_7(5),f2_7(5),f1_8(5),f2_8(5),f1_9(5),f2_9(5)
      complex*16 f1_10(5),f2_10(5),f1_11(5),f2_11(5),f1_12(5),f2_12(5)
      common     /flag0or1/ flag0or1

	do i=1,3
	   s2(i)=0
	end do
	do i=1,5
	   x1(i)=k(i)
 	   x2(i)=k(i+2)
	   f1_1(i)=fun(1,i)
	   f2_1(i)=fun(1,i+2)
	   f1_2(i)=fun(2,i)
	   f2_2(i)=fun(2,i+2)
	   f1_3(i)=fun(3,i)
	   f2_3(i)=fun(3,i+2)
	   f1_4(i)=fun(4,i)
	   f2_4(i)=fun(4,i+2)
	   f1_5(i)=fun(5,i)
	   f2_5(i)=fun(5,i+2)
	   f1_6(i)=fun(6,i)
	   f2_6(i)=fun(6,i+2)
	   f1_7(i)=fun(7,i)
	   f2_7(i)=fun(7,i+2)
	end do

	if (flag0or1(1).eq.1) then
	   sums(1)=0.0
	else 
	   sums(1)=poly4bess2(x1,f1_1,r0,0,0)+poly4bess2(x2,f2_1,r0,0,0)
	end if
      
	if (flag0or1(2).eq.1) then
	   sums(2)=0.0
	else 
	   sums(2)=poly4bess2(x1,f1_2,r0,0,1)+poly4bess2(x2,f2_2,r0,0,1)
	end if
	
	if (flag0or1(3).eq.1) then
	   sums(3)=0.0
	else 
	   sums(3)=poly4bess2(x1,f1_3,r0,1,0)+poly4bess2(x2,f2_3,r0,1,0)
	end if

	s2(1)=sums(1)+sums(2)+sums(3)
	
	if(flag0or1(4).eq.1) then
	   sums(1)=0.0
	else 
	   sums(1)=poly4bess2(x1,f1_4,r0,1,1)+poly4bess2(x2,f2_4,r0,1,1)
	end if
	
	if (flag0or1(5).eq.1) then
	   sums(2)=0.0
	else 
	   sums(2)=poly4bess2(x1,f1_5,r0,2,1)+poly4bess2(x2,f2_5,r0,2,1)
	end if
	
	s2(2)=sums(1)+sums(2)
	

	if (flag0or1(6).eq.1) then
	   sums(1)=0.0
	else 
	   sums(1)=poly4bess2(x1,f1_6,r0,0,0)+poly4bess2(x2,f2_6,r0,0,0)
      end if
	
	if (flag0or1(7).eq.1) then
	   sums(2)=0.0
	else 
	   sums(2)=poly4bess2(x1,f1_7,r0,1,0)+poly4bess2(x2,f2_7,r0,1,0)
	end if
	
	s2(3)=sums(1)+sums(2)

	return
	end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	complex*16 function poly2(x1,x2,x3,f1,f2,f3,low,high)

	real*8      x1,x2,x3,low,high,d
	complex*16	f1,f2,f3,a,b,c

	d=2./(x3-x1)**2
	a=d*(f1-2.*f2+f3)
	b=-d*((x2+x3)*f1-2.*(x1+x3)*f2+(x1+x2)*f3)
	c=d*(x2*x3*f1-2.*x1*x3*f2+x1*x2*f3)

	poly2=a/3.*(high**3-low**3)+b/2.*(high**2-low**2)+c*(high-low)
	
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine funval_S(k,o,fun,nforce)

	include "green_com.inc"

	real*8      k
	integer     nforce
	complex*16  e(4,4),a0(2),a1(2),a2(2),b1,b2,fun(7)
        common      /ab012/a0,a1,a2,b1,b2

	do lay = 1, nly
	   cpn(lay) = cdsqrt(k*k-(o/vp(lay))**2)
	   csn(lay) = cdsqrt(k*k-(o/vs(lay))**2)
	   if (real(cpn(lay)).lt.0.0) cpn(lay)=-cpn(lay)
	   if (real(csn(lay)).lt.0.0) csn(lay)=-csn(lay)
	end do
	call grt_coefs(k, o)	  
      call SourceVector_S(k,o) 
      call mtxe(lo, k, o, e)
	call Ydumtx(e)   	      
	call Uko_S

	fun(1)= b1*k*fpsv1(nforce)
	fun(2)= a0(1)*k*fpsv01(nforce)
	fun(3)=(a1(1)-b1)*k*fpsv1(nforce)  
	fun(4)= a1(1)*k*fsh1(nforce)
	fun(5)=(b1-a1(1))*k*fsh1(nforce)
	fun(6)=-a0(2)*k*fpsv01(nforce)
	fun(7)=-a1(2)*k*fpsv1(nforce)

	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      complex*16  function poly4bess2(k,f,r0,n,m)

c----------------------------------------------
c     switch=0 for integral(F(k)*Jn(kr0)*dk)
c           =1 for integral(F(k)*Jn'(kr0)*dk)
c----------------------------------------------  
      implicit none
	integer    i,j,n,m
	real*8     k(5),P(5),Q(5),R(5),S(5),a,b,c,d,ci(5)
      real*8     dk,r0,pp,sinp,cosp,sin2,cos2,sin4,cos4
      real*8     y1c,y1s,y2c,y2s,y3c,y3s,y4c,y4s,CC(5,5)
      real*8     i0c,i0s,i1c,i1s,i2c,i2s,i3c,i3s,i4c,i4s
     	complex*16 f(5),AA1(5),AA2(5)

	dk=k(4)-k(2)
	a=(k(1)-k(2))/dk
	b=(k(3)-k(2))/dk
	c=(k(4)-k(2))/dk
	d=(k(5)-k(2))/dk
	ci(1)=a*(a-b)*(a-c)*(a-d)
	ci(2)=a*b*c*d
	ci(3)=b*(a-b)*(b-c)*(b-d)
	ci(4)=c*(a-c)*(b-c)*(c-d)
	ci(5)=d*(a-d)*(b-d)*(c-d)

	CC(1,1)=0
	CC(2,1)=-(b*c*d)/ci(1)
	CC(3,1)=(b*(c+d)+c*d)/ci(1)
	CC(4,1)=-(b+c+d)/ci(1)
	CC(5,1)=1./ci(1)
	CC(1,2)=1
	CC(2,2)=-(a*(b*c+b*d+c*d)+b*c*d)/ci(2)
	CC(3,2)=(a*(b+c+d)+b*(c+d)+c*d)/ci(2)
	CC(4,2)=-(a+b+c+d)/ci(2)
	CC(5,2)=1./ci(2)
	CC(1,3)=0
	CC(2,3)=(a*c*d)/ci(3)
	CC(3,3)=-(a*(c+d)+c*d)/ci(3)
	CC(4,3)=(a+c+d)/ci(3)
	CC(5,3)=-1./ci(3)
	CC(1,4)=0
	CC(2,4)=-(a*b*d)/ci(4)
	CC(3,4)=(a*(b+d)+b*d)/ci(4)
	CC(4,4)=-(a+b+d)/ci(4)
	CC(5,4)=1./ci(4)
	CC(1,5)=0
	CC(2,5)=(a*b*c)/ci(5)
	CC(3,5)=-(a*(b+c)+b*c)/ci(5)
	CC(4,5)=(a+b+c)/ci(5)
	CC(5,5)=-1./ci(5)

  	call PQRScalc(n,k(1)*r0,P(1),Q(1),R(1),S(1))
	call PQRScalc(n,k(2)*r0,P(2),Q(2),R(2),S(2))
	call PQRScalc(n,k(3)*r0,P(3),Q(3),R(3),S(3))
	call PQRScalc(n,k(4)*r0,P(4),Q(4),R(4),S(4))
	call PQRScalc(n,k(5)*r0,P(5),Q(5),R(5),S(5))

	if (m.eq.0) then
         do i=1,5
	      AA1(i)=0.
	      AA2(i)=0.
            do j=1,5
	         AA1(i)=AA1(i)+CC(i,j)*f(j)*P(j)
	         AA2(i)=AA2(i)+CC(i,j)*f(j)*Q(j)
	      end do
	   end do
	else if(m.eq.1) then
         do i=1,5
	      AA1(i)=0.
	      AA2(i)=0.
            do j=1,5
	         AA1(i)=AA1(i)+CC(i,j)*f(j)*R(j)
	         AA2(i)=AA2(i)+CC(i,j)*f(j)*S(j)
	      end do
	   end do
	end if

	pp=dk*r0
	sinp=dsin(pp)
	cosp=dcos(pp)
	sin2=dsin(k(2)*r0)
	cos2=dcos(k(2)*r0)
	sin4=dsin(k(4)*r0)
	cos4=dcos(k(4)*r0)

	y1c= (sinp+(cosp-1)/pp)/pp
	y1s=-(cosp-sinp/pp)/pp
	y2c= (sinp-2*y1s)/pp
	y2s=-(cosp-2*y1c)/pp
	y3c= (sinp-3*y2s)/pp
	y3s=-(cosp-3*y2c)/pp
	y4c= (sinp-4*y3s)/pp
	y4s=-(cosp-4*y3c)/pp

	i0c= (sin4-sin2)/r0
	i0s=-(cos4-cos2)/r0
	i1c=dk*(cos2*y1c-sin2*y1s)
	i1s=dk*(cos2*y1s+sin2*y1c)
	i2c=dk*(cos2*y2c-sin2*y2s)
	i2s=dk*(cos2*y2s+sin2*y2c)
	i3c=dk*(cos2*y3c-sin2*y3s)
	i3s=dk*(cos2*y3s+sin2*y3c)
	i4c=dk*(cos2*y4c-sin2*y4s)
	i4s=dk*(cos2*y4s+sin2*y4c)


	poly4bess2=AA1(1)*i0c+AA1(2)*i1c+AA1(3)*i2c+AA1(4)*i3c+AA1(5)*i4c
     &         -(AA2(1)*i0s+AA2(2)*i1s+AA2(3)*i2s+AA2(4)*i3s+AA2(5)*i4s)

      return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine PQRScalc(n,x,P,Q,R,S)

	integer n,mu
	real*8  Pt,Qt,P,Q,Rt,St,R,S,pi,x
	real*8  cosv,sinv,coef,x2,x3,x4,x5,x6

	pi    = 4.0*atan(1.0)
	mu    = 4*n*n
	x2    = x*x
	x3    = x2*x
	x4    = x3*x
	x5    = x4*x
	x6    = x5*x
	cosv  = dcos((2*n+1)*pi/4)
	sinv  = dsin((2*n+1)*pi/4)
	coef  = dsqrt(2./(pi*x))

	Pt=1.-(mu-1)*(mu-9)/(128.*x2)
     %  +(mu-1)*(mu-9)*(mu-25)*(mu-49)/(98304.*x4)
     %  +(mu-1)*(mu-9)*(mu-25)*(mu-49)*(mu-81)*(mu-121)/(188743680.*x6)
	Qt=(mu-1)/(8.*x)-(mu-1)*(mu-9)*(mu-25)/(3072.*x3)
     %  +(mu-1)*(mu-9)*(mu-25)*(mu-49)*(mu-81)/(3932160.*x5)
	Rt=1.-(mu-1)*(mu+15)/(128.*x2)
     %  +(mu-1)*(mu-9)*(mu-25)*(mu+63)/(98304.*x4)  
     %  +(mu-1)*(mu-9)*(mu-25)*(mu-49)*(mu-81)*(mu+143)/(188743680.*x6)
	St=(mu+3)/(8.*x)-(mu-1)*(mu-9)*(mu+35)/(3072.*x3)
     %  +(mu-1)*(mu-9)*(mu-25)*(mu-49)*(mu+99)/(3932160.*x5)
	P=coef*(Pt*cosv+Qt*sinv)
	Q=coef*(Qt*cosv-Pt*sinv)
	R=coef*(Rt*sinv-St*cosv)
	S=coef*(Rt*cosv+St*sinv)
	
      return
	end

ccccccccccccccccccccccccccccccccc THE END ccccccccccccccccccccccccccccccc

