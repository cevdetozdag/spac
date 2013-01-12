
		include   "green_com.inc"

		integer   fff1, fff2
		real*8    f_1, f_2, amp0 ,win ,freal, window,tmp,ifftcoef
          real*8    urr(n), ufr(n), uzr(n),sur(n),suf(n),suz(n),t(n)
		complex*16   st0, sw, rrr(n), fff(n), zzz(n)
	    integer     minu,sec,nshift
	    real*8  time_begin,time_end

!	c  (0). Basic constants:
!	c  -----------------------------
		aj = cmplx(0d0, 1d0)
		pi = 4d0*atan(1d0)
		pi2= 2d0*pi

!	c  (1). Reading & Checking input:
!	c  -------------------------------   
		call green_input

!	c  (2). Basic parameters: 
!	c  ----------------------
		call green_basic
		amp0=s0

		open(11,file='Ur.dat')
		open(12,file='Uf.dat')
		open(13,file='Uz.dat')
!	c  (3). Calculate displacement spectra: Ur(r0,z0,fai0;o),

		do j=1,tmt
		   sur(j)=0.0d0
		   suf(j)=0.0d0
		   suz(j)=0.0d0
		end do

		call cpu_time(time_begin)
		do i=1,ns
		   m1=mt(i)/2
		   nf2=m1+1
		   nf1=1

!	c	fff2 is the total number of sample points in frequency domain         
		  if (WinSwitch.eq.'OFF') then
		     fff2 = m1+1 
		     fff1 = 1   
		  else 
		     fff2 = min(int(f4/df(i)), m1)+1
		     fff1 = int(f1/df(i))+1
		  end if 
		
		  call green_spectra_r(fff1, fff2,i)
	print *, '2'
!	c      pre-FFT process:  

		  ifftcoef=df(i)*mt(i)/pi2
		  do j = nf1,nf2
		     freal = df(i)*(j-1)         
		     o     = pi2*freal - aj*oi(i)
		     st0   = sw(o)
		     win   = window(freal)
		     rrr(j)=amp0*st0*win*ur(j)*ifftcoef
		     fff(j)=amp0*st0*win*uf(j)*ifftcoef
		     zzz(j)=amp0*st0*win*uz(j)*ifftcoef
		  end do
	print *, '3'
		  do j=nf2+1,mt(i)
		     rrr(j)=conjg(rrr(mt(i)+2-j))
		     fff(j)=conjg(fff(mt(i)+2-j))
		     zzz(j)=conjg(zzz(mt(i)+2-j))
		  end do

		  call fft(rrr,m,-1)
		  call fft(fff,m,-1)
		  call fft(zzz,m,-1)
	print *, '4'
!	c      post-FFT process:
		  do j=1,mt(i)
		     t(j)=(j-1)*dt
		     urr(j)=real(rrr(j))*exp(oi(i)*t(j))*1d-22
		     ufr(j)=real(fff(j))*exp(oi(i)*t(j))*1d-22
		     uzr(j)=real(zzz(j))*exp(oi(i)*t(j))*1d-22           
		  end do
		  
		  nshift=anint(shift(i)*(1.0/dt))
		  nshift=0
		  do j=1+nshift,mt(i)+nshift
		     sur(j)=sur(j)+urr(j-nshift)
		     suf(j)=suf(j)+ufr(j-nshift)
		     suz(j)=suz(j)+uzr(j-nshift)
		  end do		   
		end do
		call cpu_time(time_end)
	print *, '5'
          sec = time_end-time_begin
          minu=sec/60
	    sec=sec-minu*60

	    if (minu.eq.0)   then
	       print*,' *****************************************'
	       write(*,300) ' ** The Discrete Wavenumber Method takes '  
     &        ,sec,' seconds.   **'
	    else  
             print*,' *****************************************'
	       write(*,310) ' ** The Discrete Wavenumber Method takes '
     &         ,minu,' minutes ',sec, ' seconds.  **' 
	    end if
          pause
          
		do j=tmt,n
		   sur(j)=0.0
		   suf(j)=0.0
		   suz(j)=0.0
          end do
	
    	    do j=1,mt(i)
		   tmp=(j-1)*dt
		   write(11,111) tmp,urr(j)
		   write(12,111) tmp,ufr(j)
		   write(13,111) tmp,uzr(j)
		end do

111       format(F15.8,e20.12)
300       format(1x,a,i3,a)
310       format(1x,a,i5,a,i3,a)

		close (11)
		close (12)
		close (13)

		end   


