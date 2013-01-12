c
        subroutine green_basic

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Note: The purpose of this subroutine is to assign values to  c
c         some basic parameters, and calclulate values for some  c
c         other basic parameters.                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        include    "green_com.inc"	
	  real*8    xs(nsrc),ys(nsrc),L(nsrc),fai(nsrc)
	  real*8    xo,yo,vave,amin,amax

c  SETTING BASIC PARAMETERS:
c  -------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 call c_source(xs,ys)
	 call c_station(xo,yo)
	 call rfai(xs,ys,xo,yo,fai)
	 call fps(fai)

	 vmax = amax(nly,vp0)	  ! max velocity
	 vmin = amin(nly,vs0) 
	 vave = 0.5*(vmax+vmin)	
	 do i=1,ns
	    if (rst(i).lt.0.5) then
		m=10
		Twin(i)=6.5
	    else 
		m=11
		Twin(i)=13
	    end if
	    call mt_dk(m,i,vave)
	 end do

!	  print*,'!!!!!!!!!!!!!!!!!!! Basic Parameters !!!!!!!!!!!!!!!!!!'
!	  write(*,100)  '              Time window width : Twin =',Twin
!	  write(*,200)  '     Sample points(time domain) : mt   =',mt
!	  write(*,100)  '       Source-receiver distance : R    =',r0
!	  write(*,100)  '         Increment of time(sec) : dt   =',dt
!	  write(*,100)  '           Cutoff frequency(Hz) : fmax =',fmax
!	  write(*,100)  '     Increment of frequency(Hz) : df   =',df
!	  write(*,100)  'The period structure length(km) : L    =',L
!	  write(*,100)  ' Increment of wavenumber(km^-1) : dk   =',dk
!	  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!100     format(1x,a,f12.4)
!200     format(1x,a,i5)

c	  z(nly)=z(nly-1)+5.0

        do i = 1, nly 
           mu(i) = rho(i)*vs0(i)**2
        end do

c       Identity matrix:
        unit(1,1) = 1d0
        unit(1,2) = 0d0
        unit(2,1) = 0d0
        unit(2,2) = 1d0

        RETURN
        END

	subroutine c_source(xs,ys)

	 include "green_com.inc"
 	 real*8 xs(nsrc),ys(nsrc),lenr,a,rs(nsrc),fais(nsrc)

         pi = 4.0*atan(1.0)
	 lenr=rsmax-rsmin         
	call init_random_seed()
        do i=1,nsrc
            call random_number(a)
            rs(i)=rsmin+lenr*a
	end do
	do i=1, nsrc
	    call random_number(a)
            fais(i)=2*pi*a
            xs(i)=rs(i)*sin(fais(i))
            ys(i)=rs(i)*cos(fais(i))
	end do
	do i=1, nsrc
	    call random_number(a)
            shift(i)=tTwin*a
        end do
	
	return
	end subroutine

!-------------------------------------------------------------------------------

	subroutine c_station(xo,yo)

	 include "green_com.inc"
	 real*8 xo,yo

	 xo=0.0
	 yo=0.0
	
	 return
	end subroutine

!-------------------------------------------------------------------------------

	subroutine init_random_seed()
         integer :: n,clock
         integer,dimension(:),allocatable :: seed
         call random_seed(size=n)
         allocate(seed(n))
         call system_clock(count=clock)
         seed=clock+37*(/(i-1,i=1,n)/)
         call random_seed(put=seed)
         deallocate(seed)
         return
        end subroutine

!-------------------------------------------------------------------------------  

	subroutine rfai(xs,ys,xo,yo,fai)

	 include "green_com.inc"
	 real*8  xo,yo,xs(nsrc),ys(nsrc),fai(nsrc)		

	 pi=4.0*atan(1.0)
	 do i=1,ns
	    rst(i)=sqrt((xs(i)-xo)**2+(ys(i)-yo)**2)	 
	    if(xo.ge.xs(i)) then
	      fai(i)=acos((yo-ys(i))/rst(i))
            else if(yo.gt.ys(i)) then
       	      fai(i)=pi*2-acos((yo-ys(i))/rst(i))
            else
              fai(i)=pi+acos(-(yo-ys(i))/rst(i))
            end if
	 end do
	rst(1)=0.1
	 return
	end subroutine

!------------------------------------------------------------------------------

	subroutine fps(fai)

	 include "green_com.inc"
	 real*8 fai(nsrc)

	 do i=1,ns
            fpsv01(i) = nz
	    fpsv1(i)  = nx*dcos(fai(i)) + ny*dsin(fai(i))
	    fsh1(i)   = ny*dcos(fai(i)) - nx*dsin(fai(i))
	 end do

	 return
	end subroutine

!------------------------------------------------------------------------------

	subroutine mt_dk(m,k,vave)

	 include "green_com.inc"
	 integer  k
	 real*8 vave,L

	 pi=4.0*atan(1.0)
	 pi2=2*pi	
	 m1   = 2**(m-1)         ! Half of sampling points
       mt(k)   = 2*m1             ! Total sampling points
       dt   = Twin(k)/(mt(k)-1)
       fmax = 1d0/(2d0*dt)	      ! Nyquist frequency
       df(k)   = fmax/(m1-1)	  ! frequency interval
       oi(k)   = pi/Twin(k)   ! imaginary frequency
	  
       L=(vmax*Twin(k)+rst(k)+vmax/vave*dsqrt(rst(k)*rst(k)+(zs-z0)**2))
     &     +100      
       dk(k)   = pi2/(2*L)         ! k-integral step
	 tmt=tTwin/dt+1

	 if (f4.gt.fmax) then
             print*,'f4=',f4
	     print*,'fmax < f4 !   Please increase m to make fmax larger.'
	     pause
	     stop
	 end if

	 return
	end subroutine

