
        subroutine green_input

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Note: The purpose of this subroutine is to input parameters c
c         from the existed data files. So, before run the main  c
c         program, the input data files must be prepared in     c
c         appropriate forms.                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        include    "green_com.inc"      
        character  list*80
        integer    layernumb
        real*8     mus

c   Reading media parameters:

 101    format ( I10   )
 105    format ( a10   )
 106    format ( a90   )

        open(10,file='inmedia.dat',form='formatted',status='old')
        read(10, 106) list
        read(10, 106) list
        read(10, 106) list
        read(10, 106) list
        read(10, 101) nly
        read(10, 106) list
        do i=1,nly
           read(10, *) j,z(i-1),rho(i),vs0(i),vp0(i),qs(i),qp(i)
        end do
        close (10)
        
c   Reading BASIC, SOURCE & RECEIVER' parameters:

        open(10,file='inbasic.dat',form='formatted',status='old')
        read(10, 106) list
        read(10, 106) list
        read(10, 106) list    
        read(10, 106) list    
c    (1). Basic control parameters:
        read(10,  * ) tTwin
        read(10,  * ) coef
        read(10,  * ) criterion
c    (2). Source_type ('E':explosive_source, 'S':single_force,'D'double_couple)
        read(10, 106) list  
        read(10, 106) list
        read(10, 105) SourceType          
c    (3). Source parameters:
        read(10, 106) list
        read(10, 106) list
        read(10, 105) Type_STF
        read(10, 106) list
        read(10,  * ) fc
	read(10,  * ) tou
c    (4). Window swith and parameters:
        read(10, 106) list
        read(10, 106) list
        read(10, 105) WinSwitch
        read(10, 105) WinType
        read(10,  * ) f1
        read(10,  * ) f2
        read(10,  * ) f3
        read(10,  * ) f4
c    (5). Output Format:
        read(10, 106) list
        read(10, 106) list
        read(10, 105) OutDomain
        read(10, 105) Case
        close(10)

c       .....
c       single point source:          
c       .....
        open(10, file='in_s_s.dat')             
        read(10, 106) list
        read(10, 106) list
        read(10, 106) list
        read(10,  * ) s0
        read(10,  * ) nx
        read(10,  * ) ny
        read(10,  * ) nz
	  read(10,  * ) zs
        read(10,  * ) ns
	  read(10,  * ) rsmin
	  read(10,  * ) rsmax
        close(10)

c          ...... 

c     (5). Receiver's information:
        open(10, file='in_r_r.dat')
        read(10, 106) list   
        read(10, 106) list
        read(10, 106) list   
!	  read(10,  * ) x0
!	  read(10,  * ) y0
!	  read(10,  * ) z0
        read(10,  * ) nc
	  read(10,  * ) dc
        read(10,  * ) z0
        read(10,  * ) nr
        close(10)
c  CHECKING INPUTS:
c  ---------------------

        ls   = layernumb(zs,nly,z)
        lo   = layernumb(z0,nly,z)
       
        mus = (vs0(ls)**2)*rho(ls)

        if (nly.lt.ls.or.ls.lt.1) then
             print*,'Error: info@input-checking_1'
             stop
        end if
          
        if (nly.eq.ls) then                   
c# introduce a fictitious interface to prevent that source is within 
c  the last layer.
            if (zs.gt.z(nly-1)) then          
               print *, ' Add one fictious interface '
               z(nly)     = zs + 5.0  !!!!!!!!!                
               vs0(nly+1) = vs0(nly)   
               vp0(nly+1) = vp0(nly)
               qs(nly+1)  = qs(nly)
               qp(nly+1)  = qp(nly)
               rho(nly+1) = rho(nly)
               nly=nly+1
               print*,'nly=',nly,'     z(nly-1)=',z(nly-1)
             else
               print*,'Error: info@input-checking_2'
               stop
             end if
          end if

          if (nly .gt. NM) then
             print*,'nly>NM !!!'
             print*,'Current NM = ', NM, '  Current layer number=', nly
             stop
          endif

        END






