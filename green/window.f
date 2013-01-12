	real*8 function window(f)

	character*10 WinSwitch,WinType
	real*8       f,a,pi,f1,f2,f3,f4
      common      /filter/ WinSwitch, WinType
	common      /f1234/ f1,f2,f3,f4

      pi=4.*atan(1.0)

c	print*,'pi=',pi
c	print*,'in window, f=',f
c	print*,'f1=',f1,'    f2=',f2
c	print*,'f3=',f3,'    f4=',f4
c	print*,'WinSwitch=',WinSwitch
c	print*,'WinType=',WinType

c   CASE 1 : All_pass filter :
      if (WinSwitch.eq.'OFF') then
	   window=1.
	   return
	end if

c   CASE 2 : Hamming window :
	if(WinSwitch.eq.'ON'.and.WinType.eq.'Hamming') then

		if (f.lt.f1.or.f.gt.f4) then
			window=0.
	        return
	    end if

		if (f.le.f3.and.f.ge.f2) then
			window=1.
	        return
	    end	if

	    if (f.ge.f1.and.f.lt.f2.and.f1.ne.f2) then
			a=(f-f1)/(f2-f1)*pi
			window=(0.46-0.46*cos(a))/0.92
	        return
	    end	 if

	    if (f.gt.f3.and.f.le.f4.and.f3.ne.f4) then
			a=(f-f3)/(f4-f3)*pi
			window=(0.46+0.46*cos(a))/0.92
	        return
		end if

	end if

	end