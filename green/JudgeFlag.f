	subroutine JudgeFlag(kn,flag_sum)

      include   "green_com.inc"

	integer   flag_sum,cnt(12),nf,counter
	real*8    aeps
	parameter (aeps=1.0E-8)
	real*8    ztop(MAX0,4),zbot(MAX0,4),zzr(3),zzi(3),azzr,azzi
	real*8    ftop(MAX0,4),fbot(MAX0,4),ffr(3),ffi(3),affr,affi
	real*8    rtop(MAX0,4),rbot(MAX0,4),rrr(3),rri(3),arrr,arri
	real*8    zrsum,rrsum,frsum,zisum,risum,fisum,kmax
	integer   flag(12),flagtemp, iii

	common    /jz_nf/ nf,flag
	common    /rfztopbot/ztop,zbot,ftop,fbot,rtop,rbot
      common    /countercnt/ counter,cnt

	flagtemp=0
	flag_sum=1

	if(counter.le.3) then
c     ---------------------
	  zzr(counter)=real(uz(nf))
	  zzi(counter)=imag(uz(nf))
	  rrr(counter)=real(ur(nf))
	  rri(counter)=imag(ur(nf))
	  ffr(counter)=real(uf(nf))
	  ffi(counter)=imag(uf(nf)) 
	else if (counter.gt.3) then
c     ----------------------
	  zzr(1)=zzr(2)
	  zzr(2)=zzr(3)
	  zzr(3)=real(uz(nf))
	  azzr  =abs(zzr(1))+abs(zzr(2))+abs(zzr(3))
	  zzi(1)=zzi(2)
	  zzi(2)=zzi(3)
	  zzi(3)=imag(uz(nf))
        azzi  =abs(zzi(1))+abs(zzi(2))+abs(zzi(3))
	  rrr(1)=rrr(2)
	  rrr(2)=rrr(3)
	  rrr(3)=real(ur(nf))
        arrr=abs(rrr(1))+abs(rrr(2))+abs(rrr(3))
	  rri(1)=rri(2)
	  rri(2)=rri(3)
	  rri(3)=imag(ur(nf))
        arri=abs(rri(1))+abs(rri(2))+abs(rri(3))
	  ffr(1)=ffr(2)
	  ffr(2)=ffr(3)
	  ffr(3)=real(uf(nf))
	  affr=abs(ffr(1))+abs(ffr(2))+abs(ffr(3))
	  ffi(1)=ffi(2)
	  ffi(2)=ffi(3)
	  ffi(3)=imag(uf(nf))
	  affi=abs(ffi(1))+abs(ffi(2))+abs(ffi(3))
	  if((cnt(1).lt.MAX0.and.zzr(2).gt.zzr(1).and.zzr(2).ge.zzr(3))
     & .OR.(cnt(1).lt.MAX0.and.azzr.lt.aeps)) then
	  	 cnt(1)=cnt(1)+1
	     if(cnt(1).eq.MAX0) flag(1)=0
		 ztop(cnt(1),1)=kn
		 ztop(cnt(1),2)=zzr(2)
	  end if
	  if((cnt(2).lt.MAX0.and.zzi(2).gt.zzi(1).and.zzi(2).ge.zzi(3))
     & .OR.(cnt(2).lt.MAX0.and.azzi.lt.aeps)) then
	  	 cnt(2)=cnt(2)+1
	     if(cnt(2).eq.MAX0) flag(2)=0
		 ztop(cnt(2),3)=kn
		 ztop(cnt(2),4)=zzi(2)
	  end if
	  if((cnt(3).lt.MAX0.and.zzr(2).lt.zzr(1).and.zzr(2).le.zzr(3))
     & .OR.(cnt(3).lt.MAX0.and.azzr.lt.aeps)) then
		 cnt(3)=cnt(3)+1
	     if(cnt(3).eq.MAX0-1) flag(3)=0
		 zbot(cnt(3),1)=kn
		 zbot(cnt(3),2)=zzr(2)
	  end if
  	  if((cnt(4).lt.MAX0.and.zzi(2).lt.zzi(1).and.zzi(2).le.zzi(3))
     & .OR.(cnt(4).lt.MAX0.and.azzi.lt.aeps)) then
	 	 cnt(4)=cnt(4)+1
	     if(cnt(4).eq.MAX0-1) flag(4)=0
		 zbot(cnt(4),3)=kn
	  	 zbot(cnt(4),4)=zzi(2)
	  end if

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if((cnt(5).lt.MAX0.and.rrr(2).gt.rrr(1).and.rrr(2).ge.rrr(3))
     & .OR.(cnt(5).lt.MAX0.and.arrr.lt.aeps)) then
		 cnt(5)=cnt(5)+1
	     if(cnt(5).eq.MAX0)   flag(5)=0
	  	 rtop(cnt(5),1)=kn
		 rtop(cnt(5),2)=rrr(2)
	  end if
	  if((cnt(6).lt.MAX0.and.rri(2).gt.rri(1).and.rri(2).ge.rri(3))
     & .OR.(cnt(6).lt.MAX0.and.arri.lt.aeps)) then
		 cnt(6)=cnt(6)+1
	     if(cnt(6).eq.MAX0)   flag(6)=0
	 	 rtop(cnt(6),3)=kn
		 rtop(cnt(6),4)=rri(2)
	  end if
	  if((cnt(7).lt.MAX0.and.rrr(2).lt.rrr(1).and.rrr(2).le.rrr(3))
     & .OR.(cnt(7).lt.MAX0.and.arrr.lt.aeps)) then
	     cnt(7)=cnt(7)+1
	     if(cnt(7).eq.MAX0-1)   flag(7)=0
		 rbot(cnt(7),1)=kn
		 rbot(cnt(7),2)=rrr(2)
	  end if
	  if((cnt(8).lt.MAX0.and.rri(2).lt.rri(1).and.rri(2).le.rri(3))
     & .OR.(cnt(8).lt.MAX0.and.arri.lt.aeps)) then
		 cnt(8)=cnt(8)+1
	     if(cnt(8).eq.MAX0-1)   flag(8)=0
		 rbot(cnt(8),3)=kn
		 rbot(cnt(8),4)=rri(2)
	  end if

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if((cnt(9).lt.MAX0.and.ffr(2).gt.ffr(1).and.ffr(2).ge.ffr(3))
     & .OR.(cnt(9).lt.MAX0.and.affr.lt.aeps)) then
		 cnt(9)=cnt(9)+1
	     if(cnt(9).eq.MAX0)   flag(9)=0
		 ftop(cnt(9),1)=kn
		 ftop(cnt(9),2)=ffr(2)
	  end if
	  if((cnt(10).lt.MAX0.and.ffi(2).gt.ffi(1).and.ffi(2).ge.ffi(3))
     & .OR.(cnt(10).lt.MAX0.and.affi.lt.aeps)) then
	  	 cnt(10)=cnt(10)+1
	     if(cnt(10).eq.MAX0)   flag(10)=0
		 ftop(cnt(10),3)=kn
		 ftop(cnt(10),4)=ffi(2)
	  end if
	  if((cnt(11).lt.MAX0.and.ffr(2).lt.ffr(1).and.ffr(2).le.ffr(3))
     & .OR.(cnt(11).lt.MAX0.and.affr.lt.aeps)) then
		 cnt(11)=cnt(11)+1
	     if(cnt(11).eq.MAX0-1)   flag(11)=0
		 fbot(cnt(11),1)=kn
		 fbot(cnt(11),2)=ffr(2)
	  end if
	  if((cnt(12).lt.MAX0.and.ffi(2).lt.ffi(1).and.ffi(2).le.ffi(3))
     & .OR.(cnt(12).lt.MAX0.and.affi.lt.aeps)) then
		 cnt(12)=cnt(12)+1
	     if(cnt(12).eq.MAX0-1)   flag(12)=0
		 fbot(cnt(12),3)=kn
		 fbot(cnt(12),4)=ffi(2)
	  end if
	end if

	do i=1,12
	   flagtemp=flagtemp+flag(i)
	end do	

c      write(*,33) 'counter=',counter,'   flagtemp=',flagtemp,
c	&             '   flag_sum=',flag_sum

	if(flagtemp.eq.0) then

	  zrsum=(ztop(1,2)+ztop(MAX0,2))/2.
	  if(ztop(1,1).gt.zbot(1,1)) then
	     do i=1,MAX0-2
	        zrsum=zrsum+zbot(i,2)+ztop(i+1,2)
           end do
	     zrsum=(zrsum+zbot(MAX0-1,2))/(2.*(MAX0-1))
	  else 
	     do i=2,MAX0-1
	        zrsum=zrsum+zbot(i,2)+ztop(i,2)
           end do
	     zrsum=(zrsum+zbot(MAX0,2))/(2.*(MAX0-1))
	  end if
c--------------------------------------------------
	  zisum=(ztop(1,4)+ztop(MAX0,4))/2.
        if(ztop(1,3).gt.zbot(1,3)) then
	     do i=1,MAX0-2
	        zisum=zisum+zbot(i,4)+ztop(i+1,4)
           end do
	     zisum=(zisum+zbot(MAX0-1,4))/(2.*(MAX0-1))
	  else 
	     do i=2,MAX0-1
	        zisum=zisum+zbot(i,4)+ztop(i,4)
           end do
	     zisum=(zisum+zbot(MAX0,2))/(2.*(MAX0-1))
	  end if
c---------------------------------------------------
	  frsum=(ftop(1,2)+ftop(MAX0,2))/2.
	  if(ftop(1,1).gt.fbot(1,1)) then
           do i=1,MAX0-2
	        frsum=frsum+fbot(i,2)+ftop(i+1,2)
           end do
	     frsum=(frsum+fbot(MAX0-1,2))/(2.*(MAX0-1))
	  else 
           do i=2,MAX0-1
	        frsum=frsum+fbot(i,2)+ftop(i,2)
           end do
	     frsum=(frsum+fbot(MAX0,2))/(2.*(MAX0-1))
	  end if
c----------------------------------------------------
	  fisum=(ftop(1,4)+ftop(MAX0,4))/2.
	  if(ftop(1,3).gt.fbot(1,3)) then
           do i=1,MAX0-2
	        fisum=fisum+fbot(i,4)+ftop(i+1,4)
           end do
	     fisum=(fisum+fbot(MAX0-1,4))/(2.*(MAX0-1))
        else 
           do i=2,MAX0-1
	        fisum=fisum+fbot(i,4)+ftop(i,4)
           end do
	     fisum=(fisum+fbot(MAX0,4))/(2.*(MAX0-1))
	  end if
c---------------------------------------------------
        rrsum=(rtop(1,2)+rtop(MAX0,2))/2.
	  if(rtop(1,1).gt.rbot(1,1)) then
	     do i=1,MAX0-2
	        rrsum=rrsum+rbot(i,2)+rtop(i+1,2)
           end do
	     rrsum=(rrsum+rbot(MAX0-1,2))/(2.*(MAX0-1))
	  else 
	     do i=2,MAX0-1
	        rrsum=rrsum+rbot(i,2)+rtop(i,2)
           end do
	     rrsum=(rrsum+rbot(MAX0,2))/(2.*(MAX0-1))
        end if
c---------------------------------------------------
	  risum=(rtop(1,4)+rtop(MAX0,4))/2.
	  if(rtop(1,3).gt.rbot(1,3)) then
           do i=1,MAX0-2
	        risum=risum+rbot(i,4)+rtop(i+1,4)
	     end do
	     risum=(risum+rbot(MAX0-1,4))/(2.*(MAX0-1))
        else 
           do i=2,MAX0-1
	        risum=risum+rbot(i,4)+rtop(i,4)
	     end do
	     risum=(risum+rbot(MAX0,4))/(2.*(MAX0-1))
        end if
c---------------------------------------------------

	  uz(nf)=zrsum+aj*zisum
	  uf(nf)=frsum+aj*fisum
	  ur(nf)=rrsum+aj*risum
	  flag_sum=0
	
	end if
				 
33    format('+',3(a,i5))
      counter=counter+1

	return 
	end
