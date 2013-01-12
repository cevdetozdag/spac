

!c  "GREENMAIN_COM.INC"
!c       

        implicit  none
        integer   NM, N, MAX0, i, j, lay, nf1, nf2,nsrc,nst
        parameter ( NM=25, N=2**15+2, nsrc=50,nst=1,MAX0=5 )

!c  R/T coefficients:
        complex*16 rdu(0:NM,2,2), rud(0:NM,2,2), tdd(0:NM,2,2)
		complex*16 tuu(0:NM,2,2)
        complex*16 grdu(0:NM,2,2),grud(0:NM,2,2),gtu(0:NM,2,2)
		complex*16 gtd(0:NM,2,2)
        complex*16 rdu0(0:NM), rud0(0:NM), tdd0(0:NM), tuu0(0:NM)
        complex*16 exu(2), exd(2), ex(0:NM,2)
        complex*16 grdu0(0:NM), grud0(0:NM), gtu0(0:NM), gtd0(0:NM)
        common     /grt_r/grud, grdu, gtu, gtd, ex
        common     /grt_l/grud0, grdu0, gtd0, gtu0

!c  Basic quantities:
        complex*16 cpn(NM), csn(NM)
        complex*16 aj, o, unit(2,2)
        real*8     kn, pi, pi2
        common     /cn_ps/cpn, csn
        integer    m, m1, mt(nsrc), tmt
        common     /m_n1/tmt
		common     /m_n2/mt
        real*8     Twin(nsrc),oi(nsrc),df(nsrc),dt,dk(nsrc),tTwin
        common     /basic_1/Twin,oi,df,dk, tTwin, dt
        common     /basic_2/aj, pi, pi2, unit

!c  Quantities related to wave files:
        complex*16 yd(2,2), yu(2,2), yd0, yu0
        common     /ydu/yd, yu, yd0, yu0
        complex*16 ur2,ur1,ur0,ur01,ur02,uf2,uf1,uz2,uz1,uz01,uz02
        common     /urfz/ur2,ur1,ur0,ur01,ur02,uf2,uf1,uz2,uz1,uz01,uz02
        complex*16 su_0(2), sd_0(2), su_1(2), sd_1(2), su_2(2), sd_2(2)
        complex*16 su0_1, su0_2, sd0_1, sd0_2
        common     /SouVec_PSV/su_0, su_1, su_2, sd_0, sd_1, sd_2
        common     /SouVec_SH /su0_1, su0_2, sd0_1, sd0_2
        complex*16 cigma0(2),cigma1(2), cigma2(2), cigma1_sh, cigma2_sh
        common     /cigma/cigma0, cigma1, cigma2, cigma1_sh, cigma2_sh
!c        complex*16 ur_z(nvsp,n/2), uf_z(nvsp,n/2), uz_z(nvsp,n/2)
        complex*16 ur(0:n/2), uf(0:n/2), uz(0:n/2)
        integer    n_z
        real*8     dz0, z0_vsp
!c        common     /u_vsp/n_z, dz0, z0_vsp, ur_z, uf_z, uz_z
        common     /u_rfz/ur, uf, uz
        
!c  Source & receiver parameters:
        real*8     fc,zs,shift(nsrc),rsmin,rsmax,rst(nsrc),r0,tou
        real*8     z0,nc,dc,nr
        integer    lo, ls, nly,ns        
        common     /sou_rec/nly, lo, ls,ns
	common     /source_parameter/zs,r0,fc,shift,rst,rsmin,rsmax,tou
	common     /rec_parameter/ z0,nc,dc,nr
        real*8     nx,ny,nz,s0
        common     /single_force/ nx,ny,nz,s0
        complex*16 fpsv01(nsrc), fpsv1(nsrc),fsh1(nsrc)
        common    /radiation_pattern/fpsv01,fpsv1,fsh1
       
!c  Media parameters:
        real*8     mu(NM),z(0:NM),vs0(NM),vp0(NM),qs(NM),rho(NM),qp(NM)
        complex*16 vs(NM), vp(NM)
        common     /c_velocity/vs, vp
        common     /media/z, vp0, vs0, rho, mu, qs, qp

!c  Others:
        real*8      coef, criterion, vmin,vmax,fmax, f1,f2,f3,f4
        common      /others/coef,criterion,vmin,vmax,fmax
        common      /f1234/f1,f2,f3,f4
        character*10 SourceType, Case, Type_STF, OutDomain, OutFormat
        character*10 WinSwitch,WinType
        common     /strings/SourceType,Case,Type_STF,OutDomain,OutFormat
        common      /filter/ WinSwitch, WinType
