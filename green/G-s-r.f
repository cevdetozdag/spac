
      subroutine green_spectra_r(nf1, nf2,nforce)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Note: This is the core subroutine of the main program. The purpose c
c         of this subroutine is to calculate the Green spectrum for 	 c
c         the array-case.	    							 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	use         MSIMSL        
        include     "green_com.inc"

        integer     nf, k,jj(2),j31,j32,nn,flagflag,jz,flag(12),jflag
!	parameter   (MAX0=10)
	real*8      kmin,kmax,dk_DWIM,j0, j1, j2, j0d, j1d, j2d,r,dkk
  	real*8      ztop(MAX0,4),zbot(MAX0,4)
	real*8      ftop(MAX0,4),fbot(MAX0,4)
	real*8      rtop(MAX0,4),rbot(MAX0,4)
        complex*16  a0(2), a1(2), a2(2), b1, b2,dur,duf,duz
        complex*16  om,cs,cp,omax,otemp,fun(7)
	complex*16  sss1(3),sss2(3),sss(3),e(4,4),sums1(3)
        integer     counter,cnt(12)
	integer     judge, iii,nforce, xi
      
        common      /ab012/a0, a1, a2 , b1, b2
	common      /jz_nf/ nf,flag
	common      /rfztopbot/ztop,zbot,ftop,fbot,rtop,rbot
        common      /countercnt/ counter,cnt

!        open(2, file='CON', carriagecontrol='FORTRAN')

        dk_DWIM=dk(nforce)
	dkk=dk_DWIM
!      DO jz=1,n_z

        r0 = rst(nforce)
        print*,'nforce',nforce, 'r',r0
	print *, nf1, nf2, n
	do nf=nf1,nf2
           o = 2*pi*df(nforce)*(nf-1) - aj*oi(nforce)    ! complex*16 frequency
!	      om  = cmplx(2*pi*fc) 		      ! referenced frequency
           om=cmplx(2*pi)
	   do i=1,nly 			          ! complex*16 velocities
!	         cs = (1.0-cdlog(o/om)/(pi*Qs(i)))/(1.-aj/Qs(i)/2.)
!               cp = (1.0-cdlog(o/om)/(pi*Qp(i)))/(1.-aj/Qp(i)/2.)
              cs=(1.0+cdlog(o/(2*pi))/(pi*Qs(i))-aj/(2.0*Qs(i)))
              cp=(1.0+cdlog(o/(2*pi))/(pi*Qp(i))-aj/(2.0*Qp(i)))
              vs(i) = vs0(i)/cs
	      vp(i) = vp0(i)/cp
           end do	     

	   if (nf.le.3) then
	      otemp = 2*pi*df(nforce)*3 - aj*oi(nforce)
              call kmaxcalculate(otemp,kmax,3,nf2)
           else 
	      call kmaxcalculate(o,kmax,nf,nf2)
	   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            kmax=2.0*kmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           flagflag=1
	   counter=1
	   judge=1
	   do i=1,12
	      flag(i)=1
	      cnt(i)=0
           end do
	print *, '2a'            
           call DWIM_S(o,0d0,kmax,jj(1),dk_DWIM,sums1,i)
	    
           ur(nf)=sums1(1)
	   uf(nf)=sums1(2)
	   uz(nf)=sums1(3)

	   kn=kmax
	   
!	      jflag=1
            
           call grt_coefs(kn,o)
           call mtxe(lo,kn,o,e)
           call Ydumtx(e)
           call SourceVector_S(kn,o)
           call Uko_S
           if (cdabs(a0(1)).lt.1.0e-3.and.cdabs(a0(2)).lt.1.0e-3.and.
     &        cdabs(a1(1)).lt.1.0e-3.and.cdabs(a1(2)).lt.1.0e-3.and.
     &        cdabs(b1).lt.1.0e-3) then
	       jflag=0
	   else
	       jflag=1
	   end if
	print *, dkk
           do while(jflag.eq.1.and.flagflag.eq.1)
              kn=kn+dkk
              call funval_S(kn,o,fun,i)

              call bessj01d(r0*kn,j0,j1,j0d,j1d)

              ur(nf)=ur(nf)+(fun(1)*j0+fun(2)*j0d+fun(3)*j1d)*dkk
              uf(nf)=uf(nf)+(fun(4)*j0+fun(5)*j1d)*dkk
	      uz(nf)=uz(nf)+(fun(6)*j0+fun(7)*j1)*dkk
	      call JudgeFlag(kn,flagflag)
!	print *, xi
!	xi = xi+1
               
          end do      ! end of "k"-loop
	print *, '4a'
       end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!	END DO
       
      return
      end

