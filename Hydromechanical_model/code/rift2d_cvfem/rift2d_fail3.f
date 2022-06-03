!     last change:  map  22 jun 2004    5:19 pm
!     last change:  map  28 dec 1999    5:23 pm
!**************************************************************************
!**************************************************************************
!                                                                         *
!                             program rift2d                              *
!                                                                         *
!       "heat and brine transport in evolving rift basins"                *
!                                                                         *
!                                                                         *
!       this program is a finite element model that solves transient      *
!       groundwater flow with advective-dispersive heat transport         *
!       and petroleum generation within evolving rift and intra-          *
!       cratonic sedimentary basins. triangular elements that             *
!       employ linear interpolation functions are used to approximate     *
!       unknown hydraulic heads, temperatures, and fluid velocities.      *
!       the equation of fluid flow is formally coupled to heat and brine  *
!       transport through equations of state relating density and         *
!       viscosity to temperature, pressure, and salinity. transport       *
!       mechanisms that are represented in the program include:           *
!                                                                         *
!                                                                         *
!             1) gravity driven flow                                      *
!             2) density driven flow due to                               *
!                  a - temperature gradients                              *
!                  b - salinity gradients                                 *
!             3) compaction driven flow                                   *
!                                                                         *
!**************************************************************************
!**************************************************************************
!                                                                         *
!                               written by:                               *
!                                                                         *
!                       mark person, brian mailloux, elise bekele         *
!                       peter eadington, paul hseih, john swenson         *                                   *
!                                                                         *
!                                                                         *
!             we gratefully acknowledge contributions of:                 *
!                                                                         *
!             grant garven, denah toupin and jim wieck                    *
!             who worked on earlier versions of rift2d                    *
!                                                                         *
!                                                                         *					  *
!**************************************************************************
       program rift2d
!
!  declare and dimension all variables
!
      implicit real*8 (a-h,o-z)
! array sizes
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt
!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
! variables
!
      real*8 dt, theta, time, tja, tjb, tta,ttb
      real*8 tha,thb,sgrad, alphas, omega, delta
      real*8 domega, omegasum, total_time, msl
      real*8 dtp,dtpr,gamma,delzmn,tbrstr
      real*8 vztim1,vztim2,vztim3,vztim4,vztim5,vztim6
      real*8 vztim7,vztim8,vztim9,vztim10,vztim11,vztim12
!
      integer ibrine, iheat, ihchk, itchk, ntb, njb
      integer nfband, nhband, nnode, nelem, nmat
      integer icheck, icond, icoup, ncols
      integer nbnn, nr, iout, it, nflt, iflow, tncols,ipr
      integer tncf, nnr, rotflag, ijchk
      integer dummy,nhb,ioil,iperm,imtag
      integer ptec,pexp,psfc,pelem,iprfst
      integer ntimp,nump,numwel,maxit,ifltag,p_flag
!
! nodewise arrays
!
      real*8 xsp(maxnodes),zsp(maxnodes),x(maxnodes)
      real*8 phi(maxnodes),phiold(maxnodes),z(maxnodes)
      real*8 node_angle(maxnodes),temp(maxnodes),told(maxnodes)
      real*8 ttic(maxnodes),tti(maxnodes),tempc(maxnodes)
      real*8 vs(maxnodes),rvis(maxnodes),bb(maxnodes),head(maxnodes)
      real*8 conc(maxnodes),wrrn(maxnodes)
      real*8 rhof(maxnodes),vxp(maxnodes),vzp(maxnodes)
      real*8 xp(maxnodes),zp(maxnodes),concp(maxnodes)
      real*8 pload1(maxnodes),pload2(maxnodes)
      real*8 vsm(maxnodes),presl(maxnodes),prslold(maxnodes)
      real*8 pres(maxnodes),sigmax(maxnodes)
      real*8 phimin(maxnodes),sigmae(maxnodes)
      real*8 ss(maxnodes),perm(maxnodes),sfcel(maxnodes)
      real*8 delsfc(maxnodes),presh(maxnodes)
      real*8 heado(maxnodes),headg(maxnodes)
      real*8 rhoo(maxnodes),rhog(maxnodes)
      real*8 zmax(maxnodes),voxp(maxnodes),vozp(maxnodes)

      integer matp(maxnodes),ipflag(maxnodes),icbflg(maxnodes)
      integer ielneb(maxnodes,12),elhom(maxnodes),numpnb(maxnodes)
!
! fluid/rock property arrays
!
      real*8 perm1(0:maxmats),perm2(0:maxmats),anisop(0:maxmats)
      real*8 phi_o(0:maxmats),beta(0:maxmats)
      real*8 beta_ul(0:maxmats),phi_ir(0:maxmats)
      real*8 dif(0:maxmats),alpha(0:maxmats)
      real*8 tkf(0:maxmats),tks(0:maxmats),rhos(0:maxmats)
      real*8 cvf(0:maxmats),cvs(0:maxmats)
      real*8 ldis(0:maxmats),tdis(0:maxmats)
      real*8 perct1(0:maxmats),perct2(0:maxmats),perct3(0:maxmats)
      real*8 toc(0:maxmats),dgrn(0:maxmats)

      integer isourc(0:maxmats)
!
! boundary arrays
!
      integer njn(maxbnds),nj1(maxbnds),nj2(maxbnds),bec(maxbnds)
      integer ber(maxbnds),nh(maxbnds),nt(maxbnds),lbnd(maxbnds)
      integer rbnd(maxbnds),nrows(maxbnds),tnrows(maxbnds)
      integer nrold(maxbnds)
      integer nslp(maxbnds),nc(maxbnds),nevap(maxbnds)
      integer nclay(maxbnds),nbn(maxbnds),icase(maxbnds)

      real*8 j1(maxbnds),j2(maxbnds),j1base(maxbnds),j2base(maxbnds)
      real*8 j1inc(maxbnds),j2inc(maxbnds),tpbase(maxbnds)
      real*8 hdbase(maxbnds),hdinc(maxbnds),tpinc(maxbnds)
      real*8 cnbase(maxbnds),clay(maxbnds)
      real*8 vzt(maxbnds),zcheck(maxbnds),delz(maxbnds)
      real*8 cevap(maxbnds)
      real*8 zmin(maxbnds),vztinc(maxbnds)
      real*8 grdfac(maxbnds),ver(maxbnds),vztbase(maxbnds)

!
! elementwise arrays
!
      real*8 hflux_x(maxelems),hflux_z(maxelems),wrre(maxelems)
      real*8 twrre(maxelems),etemp(maxelems)
      real*8 area(maxelems),kx(maxelems),kz(maxelems)
      real*8 qx(maxelems),qz(maxelems),phitot(maxelems)
      real*8 tr(maxelems),rv(maxelems),masoil(maxelems)
      real*8 masgas(maxelems),mash2o(maxelems),masco2(maxelems)
      real*8 oilvol(maxelems),gasvol(maxelems),persat(maxelems)
      real*8 kgoil(maxelems), kggas(maxelems),fracoil(maxelems)
      real*8 trold(maxelems),qxo(maxelems),qzo(maxelems)

      integer mat(maxelems),ni(maxelems),nj(maxelems),nk(maxelems)
      integer eltp(maxelems),nl(maxelems)
!
! fault arrays
!
      real*8 bwidth(maxflts),centx(maxflts),fbwidth(maxflts)
      real*8 delcx(maxflts),centz(maxflts),xfault(maxflts)
      real*8 faultx(maxflts),faultz(maxflts)

      integer nftype(maxflts),nrel(maxflts),wellid(maxflts)

! character variables

      character*20 wellnm(maxflts)

!     ice loading arrays
      real*8 nu(maxbnds),nuold(maxbnds),nun(maxnodes),load(maxnodes)
      real*8 xls(maxbnds)

!
! multi-dimensional arrays
!
      real*8 aa(maxnodes,maxwdt),c(6,maxelems),b(6,maxelems)
      real*8 al(6,maxelems)

      character*40 tectlt
      character*20 file07,file08,file11,file12,file50,file65,file66
      character*20 file67,file68,file69,file70,file71
      character*20 file55,file56,file57,file58,file59,file60

!  deformation variable for lithosphere stress strain calc

       real*8 sig_kk(maxnodes),tpls(maxnodes),tplsh(maxnodes)
       real*8 xl(maxnodes),zl(maxnodes)
       real*8 dsig_dth(maxnodes),sig_1h(maxnodes),sig_3h(maxnodes)
       real*8 sigx(maxnodes),sigy(maxnodes),fc(maxnodes)
       real*8 sigxh(maxnodes),sigyh(maxnodes),fch(maxnodes)

       real*8 dsig_old(maxnodes)
       real*8 sig_1(maxnodes),sig_3(maxnodes),dsig_dt(maxnodes)
       real*8 sign_pp(maxnodes),tau_pp(maxnodes),fc_pp(maxnodes)

       integer ret(maxelems,9),elhom_l(maxelems)
       integer nnode_l,nelem_l,icount

       real*8 bl(3,maxelems),cl(3,maxelems),areal(maxelems)
       real*8 alll(3,maxelems)

! dynamic permeability variable MAP

      real*8   fc_time(maxelems),d_time,decay
      real*8   fc_crit,kx_out(maxnodes)
      integer  fc_flag(maxelems)
! ****************************************************************
!
!  assign input and output files:
!
!   #7 ---- input data file
!   #8 ---- output data file
!   #11 --- tecplot output file
!   #12 --- debugging data file
!   #50 --- surface fluid & heat flux data
!   #65 --- temperature, head, rv verse depth output data file
!   #66 --- temporal output for observation pt 1
!   #67 --- temporal output for observation pt 2
!   #68 --- temporal output for observation pt 3
!   #69 --- temporal output for observation pt 4
!   #70 --- temporal output for observation pt 5
!   #71 --- temporal output for observation pt 6


      open (unit=80,file='RIFTFILES.TXT')
      read (80,'(a)') file07
      read (80,'(a)') file08
      read (80,'(a)') file11
      read (80,'(a)') file50
      read (80,'(a)') file55
      read (80,'(a)') file56
      read (80,'(a)') file57
      read (80,'(a)') file58
      read (80,'(a)') file59
      read (80,'(a)') file60
      read (80,'(a)') file65
      read (80,'(a)') file66
      read (80,'(a)') file67
      read (80,'(a)') file68
      read (80,'(a)') file69
      read (80,'(a)') file70
      close(80)

      open (unit=7,file=file07)
      open (unit=8,file=file08)
!      open(unit=10,file='cvfem_export.txt')
      open(unit=23,file='ice_tec.dat')
      open(unit=24,file='nd1121_tec.dat')  
! node 1121 maximum under pressure in basal conf. unit node along
!  northern limb michigan basin
      open(unit=25,file='nd5931_tec.dat')
! node 1121 maximum overpressure node in basal conf. unit node along
!  southern limb illinois basin
      open(unit=26,file='nd1547_tec.dat')
! node 2366 overlying basal aquifer unit node along
      ! northern limb michigan basin
      open(unit=27,file='nd4993_tec.dat') 
! node 4993 overlying basal aquifer unit node along
!  southern limb illinois basin      
      open(unit=28,file='nd5971_tec.dat') 
      open(unit=29,file='head_debug_tec.dat')
      open(unit=30,file='nd3167_tec.dat')

! node 5971 crystallie basement node
! beyond limit of ice sheet
!      open(unit=10,file='cvfem_input_temp.txt')
      open(unit=10,file='cvfem_input_100ka.txt')
!      open(unit=10,file='cvfem_input_100ka_lund.txt')
      open(unit=12,file='femoc_input.txt')
      open (unit=11,file=file11)
      open (unit=50,file=file50)
      open (unit=55,file=file55)
      open (unit=56,file=file56)
      open (unit=57,file=file57)
      open (unit=58,file=file58)
      open (unit=59,file=file59)
      open (unit=60,file=file60)
      open (unit=65,file=file65)
      open (unit=66,file=file66)
      open (unit=67,file=file67)
      open (unit=68,file=file68)
      open (unit=69,file=file69)
      open (unit=70,file=file70)

      write (8,*) 'files opened by rift2d:'
      write (8,*) '  '
      write (8,*) '  '
      write (8,'(a)') file07
      write (8,'(a)') file08
      write (8,'(a)') file11
      write (8,'(a)') file50
      write (8,'(a)') file55
      write (8,'(a)') file56
      write (8,'(a)') file57
      write (8,'(a)') file58
      write (8,'(a)') file59
      write (8,'(a)') file60
      write (8,'(a)') file65
      write (8,'(a)') file66
      write (8,'(a)') file67
      write (8,'(a)') file68
      write (8,'(a)') file69
      write (8,'(a)') file70
      write (8,*) '  '
      write (8,*) '  '

! set up cvfem parameters

             nnode_l = 483
             nelem_l = 640
!
! set up simulation counters
!
        time =  0.0d+00
	omegasum=0.d0
	domega=0.d0
        dummy=0
        ipr=0
        it = 0
!
! initalize variables to zero
!

! dynamic permeability variable MAP

       fc_time(:) = 0
       d_time = 40
       fc_crit = 1e30 ! no failure 2.9e5
       fc_flag(:) = 0

       call zero (dt,theta,time,tja,tjb,tta,ttb,tha,thb,sgrad,
     $ alphas,omega,delta,domega,omegasum,total_time,msl,dtp,dtpr,
     $ gamma,delzmn,tbrstr,vztim1,vztim2,vztim3,vztim4,vztim5,vztim6,
     $ vztim7,vztim8,vztim9,vztim10,vztim11,vztim12,ibrine,iheat,ihchk,
     $ itchk,ntb,njb,nfband,nhband,nnode,nelem,nmat,icheck,icond,icoup,
     $ ncols,nbnn,nr,iout,it,nflt,iflow,tncols,ipr,tncf,nnr,rotflag,
     $ ijchk,dummy,nhb,ioil,iperm,imtag,ptec,pexp,psfc,pelem,iprfst,
     $ ntimp,nump,numwel,maxit,ifltag,p_flag,xsp,zsp,x,phi,phiold,z,
     $ node_angle,temp,told,ttic,tti,tempc,vs,rvis,bb,head,conc,wrrn,
     $ rhof,vxp,vzp,xp,zp,concp,pload1,pload2,vsm,presl,prslold,
     $ pres,sigmax,phimin,sigmae,ss,perm,sfcel,delsfc,presh,heado,
     $ headg,rhoo,rhog,zmax,voxp,vozp,matp,ipflag, icbflg,ielneb,elhom,
     $ numpnb,perm1,perm2,anisop,phi_o,beta,beta_ul,phi_ir,dif,alpha,
     $ tkf,tks,rhos,cvf,cvs,ldis,tdis,perct1,perct2,perct3,toc,dgrn,
     $ hflux_x,hflux_z,twrre,wrre,etemp,area,kx,kz,qx,qz,phitot,tr,rv,
     $ masoil,masgas,mash2o,masco2,oilvol,gasvol,persat,kgoil,kggas,
     $ fracoil,trold,qxo,qzo,mat,ni,nj,nk,eltp,nl,bwidth,centx,fbwidth,
     $ delcx,centz,xfault,faultx,faultz,nftype,nrel,wellid,aa,c,b,al,
     $ isourc,njn,nj1,nj2,bec,ber,nh,nt,lbnd,rbnd,nrows,tnrows,
     $ nrold,nslp,nc,nevap,nclay,nbn,icase,j1,j2,j1base,j2base,
     $ j1inc,j2inc,tpbase, hdbase,hdinc,tpinc,cnbase,clay,
     $ vzt,zcheck,delz,cevap,zmin,vztinc,grdfac,ver,vztbase)

!
!  read in input data and initialize arrays
!
      call  readin (xsp,zsp,rhos,ntime,iout,nh,dt,mat,head,tkf
     $   ,tks,cvf,cvs,ldis,tdis,phi,theta,dif,nt,iprint,temp,conc,iheat
     $   ,ibrine,icoup,ncols,nrows,nr,delz,zmin,j1base
     $   ,j2base,j1inc,j2inc,tja,tpbase,tpinc,tta
     $   ,itchk,vzt,vztbase,vztinc,vztim1,vztim2,vztim3
     $   ,vztim4,vztim5,vztim6,phiold,perm1,perm2
     $   ,anisop,toc,grdfac,isourc,sgrad,perct1,perct2
     $   ,perct3,icase,ver,iskip,tha,thb,hdbase,hdinc,ihchk
     $   ,nflt,bec,ber,iflow,nftype
     $   ,nrel,nnr,delta,domega,bwidth,centx,centz,alphas
     $   ,omega,faultx,faultz,zmax,rotflag
     $   ,node_angle,ptec,pexp,psfc,pelem,iprfst,msl
     $   ,phi_o,beta,phi_ir,beta_ul,tectlt,ioil,nc,ncb,cnbase
     $   ,numwel,wellnm,wellid,nevap,cevap,vztim7,vztim8,vztim9
     $   ,vztim10,vztim11,vztim12,gamma,tbrstr,delzmn,maxit
     $   ,nclay,clay,iperm,imtag,dgrn,ifltag,alpha,ijchk,tjb,ttb)


      write (8,*)
      write (8,*)
      write (8,*) '             **********************************'
      write (8,*) '             *  starting transient simulation *'
      write (8,*) '             **********************************'
      write (8,*)
!
!  march forward in time
!
      icount = 9

      total_time=dt*ntime
       do 100 it=1,ntime
          ipr=ipr+1
          icount = icount+1
          if(it.eq.9990.or.it.eq.19990.or.it.eq.29990) rewind 10
          if(it.eq.39990.or.it.eq.49990.or.it.eq.59990) rewind 10
          if(it.eq.69990.or.it.eq.79990.or.it.eq.89990) rewind 10
          if(it.eq.99990.or.it.eq.109990.or.it.eq.119990) rewind 10
          if(icount.eq.10) then
        call  read_cvfem(ret,xl,zl,nelem_l,nnode_l,sig_kk,
     $  sig_1,sig_3,dsig_dt,tpls,it,sigx,sigy,fc)
        if (it.eq.1) icount = 9
           if(it.gt.1) icount = 0
          end if

           write(8,*) 'after read_cvfem, xl(100) = ',xl(100)

       if(it.eq.1) then
               call areacv(xl,zl,ret,areal,alll,bl,cl,nelem_l)
       end if
                   write(8,*) 'after areacv, areal = ',areal(1)


       write (8,*) 'time step =',it
          if(ibrine.eq.0) write (*,1235)it,time
          if(ibrine.eq.1) write (*,1234)it,time,dt
          if(ibrine.eq.0) write (8,1235)it,time
          if(ibrine.eq.1) write (8,1234)it,time,dt
          if(mod(it,5).eq.0)then
          write(*,1009)nnode
          if(rotflag.eq.2) write(*,1010)omegasum
          endif
          time = time + dt
!
!  update basin geometry
!
       if(it.gt.1) call fbsep (xsp,zsp,x,z,head,vs,temp,told,tempc,tti
     $		  ,ttic,nnode,tncf,matp,ncols,nrows,phi,phiold
     $		  ,node_angle,sfcel,delsfc,conc,presl,prslold)
!
       call set_zmax (nrows,zsp,zmax,hdbase,nflt
     $                    ,bec,time,tha,hdinc,ihchk,dt)
!
! if no faults
!
      if(rotflag.eq.0)then
         call subsid (nrows,ncols,vs,zsp,dt,icase,ver,
     $	              head,rotflag,vzt)
         call angle(zsp,xsp,ncols,nrows,node_angle)
!
! if vertical faults
!
      else if(rotflag.eq.1) then
        call subsid (nrows,ncols,vs,zsp,dt,icase,ver,
     $	            head,rotflag,vzt)

        call angle(zsp,xsp,ncols,nrows,node_angle)
!
! if rotating faults
!
      else if(rotflag.eq.2)then

         call rotate (xsp,zsp,nrows,ncols,dt,icase,delta
     $	            ,omegasum,vs,bec,nflt,domega,faultx,faultz
     $              ,zmax,centx,centz,rotflag,bwidth
     $              ,fbwidth,alphas,delcx,node_angle)
       endif

!
! add or remove nodes from grid?
!
      call grdchk(ncols,nrows,zsp,zcheck,delz,icheck,grdfac,xsp,it
     $              ,vzt,dt,rotflag,icase,vs)

      if (icheck.eq.0.or.it.eq.1) then

       call gensis(nrows,ncols,nnode,xsp,zsp,vs,
     $      phi,phiold,head,temp,tempc,tti,ttic,told,delz,
     $      icase,zcheck,nrold,it,grdfac,nslp,
     $      -delta,-omegasum,vzt,dt,node_angle,sfcel,delsfc,
     $      conc,presl,prslold,rotflag,sigmax,phimin)


      end if

      call mesh(nrows,xsp,zsp,x,z,ni,nj,nk,nl,bec,ber,nelem,nnode
     $          ,nflt,eltp,mat,matp,nfband,nhband,ncols,tnrows,tncols
     $          ,lbnd,rbnd,xfault,vs,phi,phiold,head,temp,told,tempc,tti
     $ 		,ttic,tncf,dt,it,nftype,alphas,omegasum
     $     	,faultx,faultz,node_angle,sfcel,delsfc,rotflag,conc
     $          ,presl,prslold)
!
! update all boundary conditions
!
      call bound (tnrows,tncols,nbnn,nbn,nh,nt,nhb,ntb,vztbase,njn,nj1
     $            ,nj2,njb,nnode,nrmax,nhband,nflt,xfault,lbnd,rbnd
     $            ,x,z,zmax,nc,ncb)

      call bcread (time,vzt,vztbase,vztinc,vztim1,vztim2,vztim3
     & ,vztim4,vztim5,vztim6,icase,ver,ihchk,hdbase,hdinc,ncols,dt
     $ ,domega,omega,centz,ntime,nflt,rotflag,centx
     $ ,-delta,-omegasum,it,total_time,cnbase,nevap,cevap,tha,thb,vztim7
     $ ,vztim8,vztim9,vztim10,vztim11,vztim12,grdfac,mat,ber,bec
     $ ,iprint,imtag,msl)

       call bctime (time,head,nhb,nh,sgrad,
     $ j1,j2,njb,j1base,j2base,j1inc,j2inc,tja,temp,
     $ ntb,nt,tpbase,tpinc,tta,itchk,vzt,vztbase,vztinc,
     $ vztim1,vztim2,vztim3,vztim4,vztim5,vztim6
     & ,tncols,delz,dt,z,vs,grdfac,conc,hdbase,hdinc,nslp
     $ ,cnbase,nevap,cevap,nc,ncb,matp,tha,thb,vztim7,vztim8,vztim9
     $ ,vztim10,vztim11,vztim12,tnrows,nclay,clay,ihchk
     $ ,nflt,bec,msl,ijchk,tjb,ttb)

!
!  calculate elmental areas, derivatives of shape functions,
!  and zero aa matrix and b vector
!
       call areas (x,z,ni,nj,nk,area,nelem,al,b,c,nnode,time,
     $  tnrows,tncols,grdfac,vs,phi,phiold)
       call zero1 (aa,bb,nnode,nhband)
!
!  solve fluid flow equation
!
        if(iflow.eq.1) then
       if(it.eq.1) then
        dsig_old(:) = 0.0
       end if 
       if(icount.eq.0.and.it.gt.1) then
        dsig_old(:) = dsig_dth(:)
        end if

          write(*,*) 'before iceload'
          call iceld2(x,time,nu,load,nuold,it,dt,nun,z,head,
     $ nrows,ncols,xls)
         write(*,*) ' after iceload'
           call interp(ret,xl,zl,x,z,ni,nj,nk,nelem,elhom_l,
     $     dsig_dt,sig_1,sig_3,dsig_dth,sig_1h,sig_3h,
     $     alll,bl,cl,areal,it,nelem_l,nnode,
     $     head,sign_pp,tau_pp,fc_pp,tpls,tplsh,
     $     sigx,sigy,fc,sigxh,sigyh,fch,nun)
      

      call density(z,head,temp,rhof,nnode,conc)
          call viscos (z,nnode,temp,head,rvis,conc)
          do 555 iter=1,3
          call prscal(rhos,rhof,phi,z,presl,head,pres,sigmae,it,
     & matp,phi_o,beta,phi_ir,tncols,tnrows,phiold,sigmax,phimin,
     & prslold,iter)
           call porperm (perm,rhos,rhof,phi,icase,
     $ nmat,matp,phi_o,beta,phi_ir,beta_ul,perm1,perm2,sigmae,
     $ phiold,ss,sigmax,tnrows,tncols,phimin,presl,it,z,pres,vzt,
     $ p_flag,dgrn,iperm,presh,iter)
 555   continue
             call compact(z,nnode,phi,vzt,tncols,tnrows,dt,
     $               matp,phiold,vs,icase,it,rotflag,vsm,ver)
          call zero1 (aa,bb,nnode,nhband)
!          write(*,*) 'before iceload'
!          call iceld2(x,time,nu,load,nuold,it,dt,nun,z,head,
!     $ nrows,ncols,xls)
!          write(*,*) ' after iceload'
          call flow (aa,bb,x,z,ni,nj,nk,nl,rvis,rhof,kx,kz,
     $ nnode,nelem,iout,mat,area,al,b,c,phi,temp,told,dt,nfband,
     $ head,theta,phiold,perm1,perm2,anisop,vs,icoup,eltp,delta,
     $ omegasum,node_angle,ss,perm,vsm,rhos,itt,ntime,presl,prslold,
     $ beta,pload1,pload2,tr,trold,rhoo,beta_ul,icase,ifltag,alpha,
     $ rotflag,load,dsig_dth,fc_pp,dsig_old,icount,it,time,d_time,
     $ fc_time,fc_crit,kx_out,fc_flag)


          call update(temp,told,nnode)
          call febnd (aa,bb,head,nh,nhb,nnode,iout,nfband,it,ntime,
     $          hdbase,tplsh)
          call fgauss(aa,bb,nnode,nfband,head,nh,nhb,iout)

! write nodal output
       if(it.eq.1) then
       write(24,302)
       write(25,302)
       write(26,302)
       write(27,302)
       write(28,302)
       write(30,302)
       end if
          write (24,301) time,head(1131)-1.0e5,temp(1131),conc(1131),
     $  vxp(1131),vzp(1131),dsig_dth(1131),kx_out(1131)
          write (25,301) time,head(5931)-1.0e5,temp(5931),conc(5931),
     $  vxp(5931),vzp(5931),dsig_dth(5931),kx_out(5931)
          write (26,301) time,head(1547)-1.0e5,temp(1547),conc(1547),
     $  vxp(1547),vzp(1547),dsig_dth(1547),kx_out(1547)
          write (27,301) time,head(4993)-1.0e5,temp(4993),conc(4993),
     $  vxp(4993),vzp(4993),dsig_dth(4993),kx_out(4993)
          write (28,301) time,head(5971)-1.0e5,temp(5971),conc(5971),
     $  vxp(5971),vzp(5971),dsig_dth(5971),kx_out(5971)
          write (30,301) time,head(3167)-1.0e5,temp(3167),conc(3167),
     $  vxp(3167),vzp(3167),dsig_dth(3167),kx_out(3167)

 556   continue
          call velo (x,z,qx,qz,ni,nj,nk,nl,rvis,rhof,eltp,nelem
     $             ,mat,area,head,al,b,c,phi,perm,anisop
     $             ,icoup,delta,omegasum,node_angle
     $             ,qxo,qzo,rhoo,heado,ioil,pload1,ifltag,rotflag
     $             ,fc_pp,it,time,d_time,fc_time,fc_crit,fc_flag)
        endif
!
!  solve heat transport equation
!
       if (iheat.eq.1) then
      call zero1 (aa,bb,nnode,nhband)
          call heat(aa,bb,x,z,ni,nj,nk,nl,tkf,tks,qx,qz,cvs,cvf
     $             ,ldis,tdis,phi,rhof,nelem,iout,mat,area,icond,al
     $             ,b,c,rhos,dt,nhband,temp,theta,phiold,eltp
     $             ,it,ntime)
         call hnbnd(x,z,bb,nj1,nj2,njn,j1,j2,njb,it,iout,ntime)
         call hebnd(aa,bb,temp,nt,ntb,nnode,iout,nhband,it,ntime,tplsh)
         call hgauss(aa,bb,nnode,nhband,temp,nt,ntb)
        if(ioil.eq.1) then
         call kinetic (isourc,mat,it,nelem,temp,told,pres,
     $   perct1,perct2,perct3,dt,rv,tr,ni,nj,nk,toc,phi,area,head,
     $   oilvol,gasvol,persat,mash2o,masco2,masoil,masgas,fracoil,
     $   kgoil,kggas,trold)
         call petpot (head,heado,headg,rhof,nnode,z,rhoo,rhog,pres)
         call velo (x,z,qx,qz,ni,nj,nk,nl,rvis,rhof,eltp,nelem
     $             ,mat,area,head,al,b,c,phi,perm,anisop
     $             ,icoup,delta,omegasum,node_angle
     $             ,qxo,qzo,rhoo,heado,ioil,pload1,ifltag,rotflag
     $             ,fc_pp,it,time,d_time,fc_time,fc_crit,fc_flag)
      endif
      endif
!
!  solve solute transport equation using mmoc
!
      call pargen (xp,zp,nelem,nnode,nump,x,z,ni,nj,nk,nl,
     $ phi,qx,qz,vxp,vzp,numpnb,ielneb,area,nelem1,voxp,vozp,
     & qxo,qzo)


      if(ibrine.gt.0) then
      if(time.lt.tbrstr) then
      call concst (z,tncols,tnrows,conc,sgrad,cnbase)
      else
      call tstep (vxp,vzp,dtp,nelem,b,c,ni,nj,nk,dt,ntimp,dtpr
     $ ,gamma,delzmn,icbflg,nnode,maxit,total_time,ntime)

      call parmov (x,z,ni,nj,nk,nl,vxp,vzp,area,al,b,c,
     $ dtp,xp,zp,time,ielneb,numpnb,nump,concp,conc,njn,
     $ njb,ntp,nc,ncb,elhom,it,ipflag,ntime,nelem,tnrows,
     $ ncols,nnode)

      do 1000 ntp=1,ntimp+1
      timp=timp+dtp
      write(*,1111) ntp
      if(ntp.eq.ntimp+1) then
      dtp = dtpr
      call parmov (x,z,ni,nj,nk,nl,vxp,vzp,area,al,b,c,
     $ dtp,xp,zp,time,ielneb,numpnb,nump,concp,conc,njn,
     $ njb,ntp,nc,ncb,elhom,it,ipflag,ntime,nelem,tnrows,
     $ ncols,nnode)
      endif
      call parinterp (x,z,ni,nj,nk,area,al,b,c,elhom,conc
     $  ,concp,xp,zp,nump,nc,ncb,icbflg,tncols,tnrows,ntp
     $ ,numpnb,ielneb,njn,njb,nelem)

      call zero1 (aa,bb,nnode,nfband)
      call brine(aa,bb,x,z,ni,nj,nk,nl,dif,ldis,
     $ tdis,nelem,iout,mat,area,al,b,c,dtp,nfband,
     $ concp,theta,eltp,qx,qz,nnode,conc,it,rotflag,temp)
      call bebnd(aa,bb,conc,nc,ncb,nnode,iout,nfband)
      call bgauss(aa,bb,nnode,nfband,conc,nc,ncb)
 1000 continue
      endif
      endif

 302   format('variables = "time","head","temp",
     $    "conc","vx","vz","sig_kk","kx"',/,'zone')
 301   format (8(1pe16.8,2x))


!
!  calculate water rock ratio
!
      call wrratio(qx,qz,ni,nj,nk,nl,eltp,nelem,mat,area,phi
     $            ,wrre,twrre,dt,phitot,b,c,temp)
!
! calculate surface heatflow and recharge/discharge
!
         call heat_flow (hflux_x,hflux_z,ni,nj,nk,nl,eltp
     $                  ,nelem,mat,area,temp,b,c,phi
     $                  ,tkf,tks,faultx,matp)
!
!  check to see if we are at a slow printout time step
!
 456    continue
          if(pelem.eq.1)then
             call pelemout(nnr,nrel,it,mat,x,z,ni,nj,nk,qz,qx
     $         ,temp,head,ntime,iprint,msl,dt,nelem1)
          endif
       if(mod(it,iprint).eq.0)then
       write (*,*) 'calling tecout, it:',it,'  iprint: ',iprint
!
! generate graphical tecplot output
!
         if(ptec.eq.1.and.it.gt.iprfst)then
            call tecout (nnode,nelem,area,phi,ni,nj,nk,nl
     $  ,x,z,temp,head,matp,dummy,wrre,phitot,mat,tncols
     $  ,tnrows,time,ntime,iprint,it,nflt,wrrn,msl
     $  ,pres,presl,sigmae,vsm,ss,rhof,tectlt,tr
     $  ,masoil,masgas,rv,oilvol,gasvol,persat,mash2o,masco2
     $  ,heado,headg,rhoo,rhog,conc,perm,iflow,ibrine,iheat,ioil
     $  ,vxp,vzp,xp,zp,concp,numwel,wellnm,wellid,ipflag
     $  ,load,pload2,prslold,voxp,vozp,presh,node_angle
     $  ,anisop,ifltag,elhom,dsig_dth,sig_1h,sig_3h
     $  ,sign_pp,tau_pp,fc_pp,nu,xls,sigxh,sigyh,fch,nun
     $  ,d_time,fc_time,fc_crit,kx_out)
          write(*,*) ' after tecout'
             dummy=1
          endif
!
! general numerical output
!
        call output (x,z,head,temp,tti,qx,qz,
     $ nnode,it,iprint,dt,
     $ nelem,ipr,conc,iheat,ibrine,time,icoup,etemp,
     $ ni,nj,nk,nl,mat,nbn,nbnn,tempc,ttic,phi,kx,kz,rvis,
     $ tnrows,rhof,tr,masoil,masgas,rv,oilvol, gasvol,persat,mash2o,
     $ masco2,tncols,kgoil,kggas,iskip,grdfac,ber,bec,nflt,ntime,
     $ sigmae,pres,presl,ss,vsm,heado,headg,rhoo,rhog,ioil)

       endif
         if(mod(it,iprint).eq.0)then
!
! print surface flux output
!
	  if(psfc.eq.1.and.it.gt.iprfst)then
	     call sfcout (nnode,nelem,area,phi,ni,nj,nk,nl,qx,qz,x
     $                  ,z,temp,head,matp,wrre,phitot,iflow
     $                  ,mat,hflux_z,tncols,tnrows,time,ntime,iprint,it
     $                  ,bec,nflt,wrrn,iprfst,msl)
	  endif
        endif
 100   continue

c
      write (8,*)
      write (8,*)
      write (8,*) '             **********************************'
      write (8,*) '             *     simulation is complete     *'
      write (8,*) '             **********************************'
      write (8,*)

!
!  simulation complete
!
       close(7)
       close(8)
       close(11)
       close(12)
       close(18)
       close(50)
       close(65)
       close(66)
       close(67)
       close(68)
       close(69)
       close(70)
       close(71)
       write(*,*)'stop'
!***********************************************************************
!
!  format statements
!
 678  format(3i6,2(1pe12.4,1x))
 1111 format(3x,'  particle move = ',i10)
 1009 format(3x,'  number of nodes    = ',i6)
 1010 format(3x,'  amount of rotation = ',f6.2)
 1234 format(3x,'  time step # = ',i6,3x,'time(yr)=',1pe12.4,3x,
     $ 		'dt(yr)=',1pe12.4)
 1235 format(3x,'  it = ',i6,3x,'time(yr)=',1pe12.4)
      stop
      end
!**************************************
       subroutine areacv(xl,zl,ret,areal,alll,bl,cl,nelem_l)

       implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)

       integer  maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt
       integer maxcols,it,nnode,nelem

       include 'rift2d_dimensions.txt'

       real*8 xl(maxnodes),zl(maxnodes)
       real*8 bl(3,maxelems),cl(3,maxelems),areal(maxelems)
       real*8 alll(3,maxelems)

       integer ret(maxelems,9)
       integer i,j,k,m,nelem_l


       do 100 m=1,nelem_l
       i = ret(m,1)
       j = ret(m,2)
       k = ret(m,3)


       alll(1,m) = xl(j)*zl(k) - xl(k)*zl(j)
       alll(2,m) = xl(k)*zl(i) - xl(i)*zl(k)
       alll(3,m) = xl(i)*zl(j) - xl(j)*zl(i)


       bl(1,m)= zl(j)-zl(k)
       bl(2,m)= zl(k)-zl(i)
       bl(3,m)= zl(i)-zl(j)
       cl(1,m)= xl(k)-xl(j)
       cl(2,m)= xl(i)-xl(k)
       cl(3,m)= xl(j)-xl(i)

      areal(m) = 0.5d+00*((xl(i)*zl(j)-xl(j)*zl(i))
     $  + (xl(k)*zl(i)-xl(i)*zl(k))
     $  + (xl(j)*zl(k)-xl(k)*zl(j)))

  100 continue
      return
      end
!******************************************************
       subroutine read_cvfem(ret,xl,zl,nelem_l,nnode_l,sig_kk,
     $  sig_1,sig_3,dsig_dt,tpls,it,sigx,sigy,fc)


      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

       integer  maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt
       integer maxcols,it


      include 'rift2d_dimensions.txt'

      integer ret(maxelems,9)

      real*8 xl(maxnodes),zl(maxnodes),tpls(maxnodes)
      real*8 sig_kk(maxnodes),sigx(maxnodes),sigz(maxnodes)
      real*8 sigy(maxnodes),sigxy(maxnodes),sig_1(maxnodes)
      real*8 sig_3(maxnodes),dsig_dt(maxnodes)
      real*8 ux,uz,dx,dz
      real*8 signrm,tau,fc(maxnodes)
      character*80 title
      write(*,*) 'in read_cvfem'

!      if(it.eq.1) read(10,400) title
!      if(it.eq.1) write(8,400) title      
      read(10,400) title
      write(8,400) title
      write(*,*) 'starting read of nodal data'
      do n=1,nnode_l
      read(10,*) xl(n),zl(n),dx,dz,sig_kk(n),dsig_dt(n),
     $ ux,uz,sigx(n),sigy(n),sig_1(n),sig_3(n),tpls(n),
     $ fc(n),signrm,tau

      sigy(n) = -sigy(n)

      write(8,100) n,xl(n),zl(n),dx,dz,sig_kk(n),dsig_dt(n),
     $ ux,uz,sig_1(n),sig_3(n),tpls(n)

      end do
      write(*,*) 'starting read of element conectivity data'

      do m=1,nelem_l
      read(10,*) ret(m,1),ret(m,2),ret(m,3) ! ,ret(m,4),ret(m,5), 
!     $  ret(m,6),ret(m,7),ret(m,8),ret(m,9)         
      write(8,200) ret(m,1),ret(m,2),ret(m,3) ! ,ret(m,4),ret(m,5),
!     $  ret(m,6),ret(m,7),ret(m,8),ret(m,9)
      end do

  100  format('n=',i6,2x,'xl=',1pe12.4,2x,
     $    'yl=',1pe12.4,2x,'dx=',1pe12.4,2x,
     $    'dz=',1pe12.4,2x,'sig_kk=',1pe12.4,2x,'dsig_dt=',1pe12.4,2x
     $    'ux=',1pe12.4,'uz=',1pe12.4,2x,
     $    'sig1=',1pe12.4,2x,'sig_3=',1pe12.4)
  200  format('ret1=',i5,2x,'ret2=',i5,2x,'ret3=',i5)

!   ,2x,'ret4=',i5,2x,
!     $  'ret5=',i5,2x,'ret6=',i5,2x,'ret7=',i5,2x,'ret8=',i5,2x,
!     $ 'ret9=',i)
  400  format(a80)
      return
      end
! *****************************************************************

        subroutine interp(ret,xl,zl,x,z,ni,nj,nk,nelem,elhom_l,
     $  dsig_dt,sig_1,sig_3,dsig_dth,sig_1h,sig_3h,
     $  alll,bl,cl,areal,it,nelem_l,nnode,head,sign_pp,tau_pp,fc_pp,
     $  tpls,tplsh,sigx,sigy,fc,sigxh,sigyh,fch,nun)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

!
! array sizes
!
       integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt
       integer maxcols,it,maxrows

       include 'rift2d_dimensions.txt'

      integer ret(maxelems,9),nde(9),nde_h(3)
      integer elhom_l(maxelems),i,j,k
      integer ni(maxelems),nj(maxelems),nk(maxelems)
      integer q,m,n,ll,nelem_l,nnode_l,nnode,nelem
      integer node(3),ip,ii

      real*8 sigx(maxnodes),sigy(maxnodes),fc(maxnodes)
      real*8 sigxh(maxnodes),sigyh(maxnodes),fch(maxnodes)
      real*8 dist_h(maxelems),dist,head(maxnodes)
      real*8 xl(maxnodes),zl(maxnodes),x(maxnodes),z(maxnodes)
      real*8 sig_1(maxnodes),sig_3(maxnodes)
      real*8 areal(maxelems),bl(3,maxelems),cl(3,maxelems)
      real*8 alll(3,maxelems),tpls(maxnodes),tplsh(maxnodes)
      real*8 sig_1hp(maxnodes),sig_3hp(maxnodes)
      real*8 sign_pp(maxnodes),tau_pp(maxnodes),fc_pp(maxnodes)
      real*8 sig_1h(maxnodes),sig_3h(maxnodes)
      real*8 shapi,shapk,shapj
      real*8 dsig_dt(maxnodes),dsig_dth(maxnodes)

      real*8 min_dist,xxc(maxelems),zzc(maxelems)
      real*8 xc(maxelems),zc(maxelems),xxxc,zzzc
      real*8 dummy,nun(maxnodes)
!----------------------------------------------------------------------
!   map the source term (dsig_dt)(sig_1)(sig_3) from the i
!  deformation grid to the flow  gird
!  
!----------------------------------------------------------------------
         write(8,*) 'top interp: ret(1,1)=',ret(1,1)
   
         if(it.eq.1) then

          do ll=1,nelem_l   !  lithosphere elements
          xxc(ll) = (xl(ret(ll,1))+xl(ret(ll,2))+xl(ret(ll,3)))/3 ! xl(ret(ll,9))
          zzc(ll) = (zl(ret(ll,1))+zl(ret(ll,2))+zl(ret(ll,3)))/3
          end do

          do  m=1,nnode    ! hydrogeologic elements
           xc(m) = x(m) ! (x(ni(m))+x(nj(m))+x(nk(m)))/3
           zc(m) = z(m) ! (z(ni(m))+z(nj(m))+z(nk(m)))/3
           end do
           do m=1,nnode ! hydrogeologic elements
             min_dist = 1.0e+8
             do ll=1,nelem_l
               dist = sqrt( (xc(m)-xxc(ll))**2 + (zc(m)-zzc(ll))**2)
               if(dist.lt.min_dist) then
               min_dist = dist
               elhom_l(m) = ll
               dist_h(m) =  min_dist
               end if
             end do

             write(8,*) 'node=',m
             write(8,*) 'dist_h(m)=', dist_h(m)
             write(8,*) 'elhom_l(m)=',elhom_l(m)

           end do
!23456   
        end if
!        write(8,*) 'after mindist_loop'
        do q=1,nnode
        m = elhom_l(q)
!  i,j,k are indices of lithosphere element

        i = ret(m,1)
        j = ret(m,2)
        k = ret(m,3)
! node are indices of hydrologic element

!        node(1) = ni(q)
!        node(2) = nj(q)
!        node(3) = nk(q)

!        do ll=1,3      
!        ip = node(ll)

        xxxc = xc(q)
        zzzc = zc(q)

        shapi=0.5d+0*(alll(1,m)+bl(1,m)*xxxc+cl(1,m)*zzzc)/areal(m)
        shapj=0.5d+0*(alll(2,m)+bl(2,m)*xxxc+cl(2,m)*zzzc)/areal(m)
        shapk=0.5d+0*(alll(3,m)+bl(3,m)*xxxc+cl(3,m)*zzzc)/areal(m)

!  calculate positions of nodes of the hydrologic mesh within
!     lithosphere mesh using local coordinate system

        dsig_dth(q) = shapi*dsig_dt(i)+shapj*dsig_dt(j)
     $              + shapk*dsig_dt(k)
        sig_1h(q) = shapi*sig_1(i)+shapj*sig_1(j)
     $              + shapk*sig_1(k)
        sig_3h(q) = shapi*sig_3(i)+shapj*sig_3(j)
     $              + shapk*sig_3(k)
        tplsh(q) = shapi*tpls(i)+shapj*tpls(j)
     $              + shapk*tpls(k)
        sigxh(q) = shapi*sigx(i)+shapj*sigx(j)
     $              + shapk*sigx(k)
        sigyh(q) = shapi*sigy(i)+shapj*sigy(j)
     $              + shapk*sigy(k)
        fch(q) = shapi*fc(i)+shapj*fc(j)
     $              + shapk*fc(k)

!23456 


!        write(8,*) 'q=',q
!        write(8,*) 'm=',m
!        write(8,*) 'i,j,k',i,j,k
!        write(8,*) 'shapi,shapj,shapk',shapi,shapj,shapk
!        write(8,*) 'dsig_dt',dsig_dt(i),dsig_dt(j),dsig_dt(k)
!        write(8,*) 'sig_1',sig_1(i),sig_1(j),sig_1(k)
!        write(8,*) 'sig_3',sig_3(i),sig_3(j),sig_3(k)
!        write(8,*) 'disg_dth,sig_1h,sig_3h',dsig_dth(q),sig_1h(q),
!     $   sig_3h(q)
!            end do
        end do


        do n = 1,nnode
        dummy = head(n)-100000.0-0.9*nun(n)  ! convert to dev. head
        
        sig_1hp(n) = sig_1h(n)- dummy*10000
        sig_3hp(n) = sig_3h(n)- dummy*10000
        sign_pp(n) = 0.5*(sig_1hp(n)+sig_3hp(n))
     $   - 0.5*(sig_1hp(n)-sig_3hp(n))*cos(3.1415/3)
        tau_pp(n) = 0.5*(sig_1hp(n)-sig_3hp(n))*sin(3.1415/3)
        fc_pp(n) = tau_pp(n) - 0.6*sign_pp(n)
!        write(8,300) x(n),z(n),dsig_dth(n),sig_1h(n),sig_3h(n),
!     $   sign_pp(n),tau_pp(n),fc_pp(n)
!23456
         end do
 300   format('x=',1pe14.4,2x,'x=',1pe14.4,2x,'z=',
     $   1pe14.4,2x,'dsig_dth=',1pe14.4,2x,
     $  'sig_3h=',1pe14.4,2x,'sign_pp=',1pe14.4,2x,
     $  'tau_pp=',1pe14.4,2x,'fc_pp=',1pe14.4)
        return
        end


       subroutine iceld2(x,time,nu,load,nuold,it,dt,nun,y,head,
     $ nrows,ncols,xls)


      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

       !
       ! array sizes
       !
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt
      integer maxcols,it

      include 'rift2d_dimensions.txt'

       real*8 x(maxnodes),ho,lo,nu(maxbnds),nuold(maxbnds)
       real*8 load(maxnodes),lt,ht,delt0,delt1,delt2,delt3
       real*8 delx_h,time,hsl2,lsl2,hsl1,lsl1,xls(maxbnds)
       real*8 nun(maxnodes),head(maxnodes),y(maxnodes)
       real*8 newtime,maxtime

       integer nrows(maxbnds),ncols,m,n

       ho = 3000.0
       lo = 1.3e6
       delt0 = 0.0e4
       delt1 = 1.9e4
       delt2 = 2.0e4
       delt3 = 3.0e4
       delx_h = 16000.0
       maxtime = 3.0d+4
       newtime =mod(time,maxtime)
!       write(*,*) 'it = ',it
       write(*,*) ' ice load, time = ',newtime
       if(it.eq.1) then
               xx = 0.0
               do  m=1,ncols
               xls(m) = xx
               xx = xx + delx_h
               nuold(m) = 0.0
               nu(m) = 0.0
               end do
               icenodes = ncols
        end if

        if(newtime.le.delt0) then
                lt = 1.0
                ht = 0.0
        end if


       if (newtime.gt.delt0.and.newtime.le.delt1) then
        lsl1 = lo/(delt1-delt0)
        hsl1 = ho/(delt1-delt0)
        lt = lsl1*(newtime-delt0)
        ht = hsl1*(newtime-delt0)
        end if

        if(newtime.gt.delt1.and.newtime.le.delt2) then
                lt = lo
                ht = ho
        end if

        if(newtime.gt.delt2.and.newtime.le.delt3) then
        lt = lo
        ht = ho
        lsl2 = lo/(delt3-delt2)
        hsl2 = ho/(delt3-delt2)
        lt = lo - lsl2*(newtime-delt2)
        ht = ho - hsl2*(newtime-delt2)
        end if
   
        if(newtime.ge.delt3) then
                lt = 1
                ht = 0
        end if

       do n=1,ncols
       if(it.gt.1) nuold(n) = nu(n)
       arg = abs( 1.0 - ( xls(n)/lt )**2)
       nu(n)= ht*( arg**0.5 )
       if(xls(n).gt.lt) then
        nu(n) = 0.0
       end if
       write(12,404) n,xls(n),arg,ht,lt,nu(n)
 404   format('n=',i6,2x,'xls=',1pe14.3,3x,'arg=',1pe14.3,3x,'ht=',
     $  1pe14.3,2x,'lt=',1pe14.3,2x,'nu=',1pe14.2)
       if(nu(n).lt.0) nu(n) = 0.0
       if(nu(n).gt.ho) nu(n) = ho

!      write(*,*) 'xls=',xls(n)
!      write(*,*) 'nu=',nu(n)
!      write(*,*) 'n=',n
       end do

        nn = 0
        nnn = 0
        do 600 n=1,ncols
         nn = nn + nrows(n)
           do 400 m=1,nrows(n)
            nm = nn - m  + 1
              load(nm) = (nu(n)-nuold(n))/dt
              nun(nm) = nu(n)
 400   continue
               if(nu(n).le.0) head(nn) = y(nn)
               if(nu(n).gt.0) head(nn) = y(nn)+nu(n)*0.9

600    continue

       return
       end


!***********************************************************************
       subroutine angle(zsp,xsp,ncols,nrows,node_angle)
!
!     this subroutine calculates the elemental row slope based on
!     adjacent nodal column elevations.
!
! local variables used:
!
!  nn = bottom of current column
!  mm - bottom of column to the right
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 node_angle(maxnodes),zsp(maxnodes),xsp(maxnodes),delx,delz
      real*8 theta(maxbnds)
      real*8 atand

      integer nrows(maxbnds),ncols,nn,nnn,mm,n,m

      nn = 1
      nnn = 0
      do 100 n=1,ncols-1
      nnn=nnn+nrows(n)
      mm = nnn+1
      delz = zsp(mm) - zsp(nn)
      delx = xsp(mm) - xsp(nn)
      if(delx.gt.0.0d+0) then
      theta(n) = atan(delz/delx)
      else
      theta(n) = theta(n-1)
      endif
      do 200 m=1,nrows(n)
      node_angle(nn) = theta(n)
      nn=nn+1
  200 continue
      nn=nnn+1
  100 continue
      theta(ncols) = theta(ncols-1)
      do 500 m=1,nrows(ncols)
      node_angle(nn) = theta(ncols)
       nn=nn+1
  500 continue
!***********************************************************************
!  format statements
!
  201 format ('node1=',i5,3x,'node2=',i5,3x,'delx=',f9.1,3x,
     $   'delz=',f9.1,3x,'angle=',f12.4)

 250  return
      end
!***********************************************************************
      subroutine prscal(rhos,rhof,phi,z,presl,head,pres,sigmae,it,
     & matp,phi_o,beta,phi_ir,tncols,tnrows,phiold,sigmax,phimin,
     & prslold,iter)
!
!  this subroutine calculates the total stress (lithostatic pressure),
!  fluid pressure, and effective stress, for each node
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 rhos(0:maxmats),rhof(maxnodes),presl(maxnodes)
      real*8 phi(maxnodes),z(maxnodes),phimin(maxnodes),atemp
      real*8 g,depth,psum,sigmae(maxnodes),pres(maxnodes),head(maxnodes)
      real*8 phi_o(0:maxmats),beta(0:maxmats),phi_ir(0:maxmats)
      real*8 prslold(maxnodes),phiold(maxnodes),sigmax(maxnodes)
      real*8 dlz1,rhoff

      integer nn,nnn,m,n,tncols,tnrows(maxbnds),matp(maxnodes),it,q
      integer iter

      g = 9.81d+0
      rhoff = 1000.0d+0
!
!  if this is the first time step, set porosities assuming
!  hydrostatic initial conditions
!
      if(it.eq.1) then
      nn=0
       do 102 n=1,tncols
       do 202 m=1,tnrows(n)
        nn=nn+1
        nnn=nn+tnrows(n)-m
        q=matp(nn)
        depth = z(nnn)-z(nn)
        atemp =  -beta(q)*depth
        phi(nn) = phi_o(q)*dexp(atemp) + phi_ir(q)
        phiold(nn) = phi(nn)
        phimin(nn) = phi_o(q)

  202  continue
  102 continue

      do 1300 ii=1,10

      nnn = 0
      do 1100 n=1,tncols
      psum = 0.0d+0
      nnn = nnn + tnrows(n)
      nn = nnn
      do 1200 m=1,tnrows(n)-1
        nn=nn-1
        q=matp(nn)
        prslold(nn) = presl(nn)
        dlz1 = z(nn+1) - z(nn)
        psum  =  psum+ dlz1*((1.0-phi(nn))*rhos(q)*g
     $ +phi(nn)*rhof(nn)*g)
        presl(nn) = psum
        pres(nn) =  rhoff*g*(head(nn) - z(nn))
       	sigmae(nn) = presl(nn) - pres(nn)
       	if(pres(nn).gt.presl(nn)) sigmae(nn) = 0.0d+0
        sigmax(nn) = sigmae(nn)
        atemp = -beta(q)*sigmae(nn)/(g*(rhos(q)-rhof(nn)))
        phi(nn) = phi_o(q)*dexp(atemp) + phi_ir(q)
        phiold(nn) = phi(nn)
        if(phimin(nn).gt.phi(nn)) phimin(nn) = phi(nn)

 1200 continue

      q=matp(nnn)
      presl(nnn) = 0.0d+0
      prslold(nnn) = presl(nnn)
      pres(nnn) =  0.0d+0
      sigmae(nnn) = 0.0d+00
      sigmax(nnn) = sigmae(nnn)
      phi(nnn) = phi_o(q) + phi_ir(q)
      phiold(nnn) = phi(nnn)

 1100 continue
 1300 continue
 1301 format ('n=',i6,1x,'phi=',1pe12.4,1x,'pres=',1pe12.4,1x,
     $ 'presl=',1pe12.4,1x,'sigmae=',1pe12.4)
 1302 format ('rhos=',1pe12.4,1x,'rhof=',1pe12.4,1x,
     $ 'head=',1pe12.4,1x,'z=',1pe12.4)


      endif
      if(it.gt.1) then
      if(iter.eq.1) then
      nn = 0
      do 105 n=1,tncols
      do 205 m=1,tnrows(n)
        nn=nn+1
      	prslold(nn) = presl(nn)
  205 continue
  105 continue
      endif
      nnn = 0
      do 100 n=1,tncols
      psum = 0.0d+0
      nnn = nnn + tnrows(n)
      nn = nnn
      do 200 m=1,tnrows(n)-1
        nn=nn-1
        q=matp(nn)
        dlz1 = z(nn+1) - z(nn)
        psum  =  psum+ dlz1*g*((1.0-phi(nn))*rhos(q)
     $ +phi(nn)*rhof(nn))
      	presl(nn) = psum
        pres(nn) =  rhoff*g*(head(nn) - z(nn))
       	sigmae(nn) = presl(nn) - pres(nn)
       	if(pres(nn).gt.presl(nn)) sigmae(nn) = 0.0d+0
        if(sigmae(nn).gt.sigmax(nn)) sigmax(nn) = sigmae(nn)
  200 continue

      q=matp(nnn)
      presl(nnn) = 0.0d+0
      prslold(nnn) = presl(nnn)
      pres(nnn) =  0.0d+0
      sigmae(nnn) = 0.0d+00
      sigmax(nnn) = sigmae(nnn)
      phimin(nnn) = phi_o(q) + phi_ir(q)
      phiold(nnn) = phi(nnn)
  100 continue
      endif
      if(it.eq.1) then
      do 1105 n=1,tncols
      nn = 0
      do 1205 m=1,tnrows(n)
        nn=nn+1
      	prslold(nn) = presl(nn)
 1205 continue
 1105 continue
      endif
!***********************************************************************
!  format statements
!
  556 format(' node=',i7,' z=',f10.3,' presl=',1pe12.3,
     $ ' preslold=',1pe12.3)
  555 format(' node=',i7,' top=',i7,' phi=',f5.2,
     $ ' depth=',f6.0,' rhos=',f6.0,
     $ ' rhof=',f6.0,' matp',i7,' presl=',1pe12.3)

      return
      end
!***********************************************************************
      subroutine concst (z,tncols,tnrows,conc,sgrad,cnbase)
!
!  this subroutine imposes concentration gradients on variable conc
!  for simulation times less than tbrstr
!
!***********************************************************************

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 z(maxnodes),depth,conc(maxnodes),sgrad,cnbase(maxbnds)

      integer nn,nnn,m,n,tncols,tnrows(maxbnds)

      nn=0
      nnn=0
      do 100 n=1,tncols
        nnn=nnn+tnrows(n)
       do 200 m=1,tnrows(n)
        nn=nn+1
        depth = z(nnn)-z(nn)
        conc(nn) = cnbase(n) + depth*sgrad
  200  continue
  100 continue
      return
      end
!***********************************************************************
      subroutine parmov (x,z,ni,nj,nk,nl,vxp,vzp,area,al,b,c,
     $ dtp,xp,zp,time,ielneb,numpnb,nump,concp,conc,njn,
     $ njb,ntp,nc,ncb,elhom,it,ipflag,ntime,nelem,tnrows,ncols,
     $ nnode)
!
!
!   this subroutine locates particles in finite element mesh
!   and then moves particles to new locations
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 a1,l1,c1,e1,d1,l1test,diff1
      real*8 a2,l2,c2,e2,d2,l2test,diff2
      real*8 area,x,z
      real*8 dtp,vxp,vzp,xp,zp,conc,concp
      real*8 areai,areaj,areak,resid,totarea

      integer ni,nj,nk,nl,nump,ipflag(maxnodes),elhom(maxnodes),ntime
      integer ip,i,j,k,numpnb,ielneb(maxnodes,12),njn(maxbnds),njb,ntp
      integer tnrows(maxbnds),nelem

      dimension numpnb(maxnodes),conc(maxnodes),nl(maxelems)
      dimension area(maxelems),ni(maxelems),nj(maxelems),nk(maxelems)
      dimension x(maxnodes),z(maxnodes),vxp(maxnodes),vzp(maxnodes)
      dimension al(6,maxelems),b(6,maxelems),c(6,maxelems)
      dimension concp(maxnodes),xp(maxnodes),zp(maxnodes),nc(maxbnds)

!
!  move the particles backwards using dtp
!
       do 100 ip=1,nump
       xp(ip) = x(ip) - vxp(ip)*dtp
       zp(ip) = z(ip) - vzp(ip)*dtp
       ipflag(ip) = 0
       elhom(ip) = 0
  100 continue

!
! now find out what element the particle has moved into
! start by looking in elements surrounding the node
!
       do 200 ip=1,nump

  115  continue
       do 300 n=1,numpnb(ip)
       m = ielneb(ip,n)
       if(ni(m).eq.0) go to 300
       if(ipflag(ip).gt.0) go to 300

        i = ni(m)
        j = nj(m)
        k = nk(m)

        areai = 0.5d+00 * ((xp(ip)*z(j) - x(j)*zp(ip)) +
     $       (x(k)*zp(ip) - xp(ip)*z(k)) + (x(j)*z(k) - x(k)*z(j)))

        areaj = 0.5d+00 * ((x(i)*zp(ip) - xp(ip)*z(i)) + (x(k)*z(i) -
     $          x(i)*z(k)) + (xp(ip)*z(k) - x(k)*zp(ip)))

        areak = 0.5d+00 * ((x(i)*z(j) - x(j)*z(i)) + (xp(ip)*z(i) -
     $          x(i)*zp(ip)) + (x(j)*zp(ip) - xp(ip)*z(j)))

      if((areai.lt.0.0d+0).or.(areaj.lt.0.0d+0)
     $ .or.(areak.lt.0.0d+0)) go to 300
      totarea = areai + areaj + areak
      resid = abs(area(m) - totarea)
      if (resid .le. 0.001) then
!
! okay, we have found which element the particle resides in
! now determine its concentration using shape function weighting
! scheme
!
       ipflag(ip) = 1
       elhom(ip) = m
      endif
  300 continue

       if(ipflag(ip).gt.0) go to 200
       do 399 m=1,nelem
       if(ni(m).eq.0) go to 399
       if(ipflag(ip).gt.0) go to 399

        i = ni(m)
        j = nj(m)
        k = nk(m)

        areai = 0.5d+00 * ((xp(ip)*z(j) - x(j)*zp(ip)) +
     $       (x(k)*zp(ip) - xp(ip)*z(k)) + (x(j)*z(k) - x(k)*z(j)))

        areaj = 0.5d+00 * ((x(i)*zp(ip) - xp(ip)*z(i)) + (x(k)*z(i) -
     $          x(i)*z(k)) + (xp(ip)*z(k) - x(k)*zp(ip)))

        areak = 0.5d+00 * ((x(i)*z(j) - x(j)*z(i)) + (xp(ip)*z(i) -
     $          x(i)*zp(ip)) + (x(j)*zp(ip) - xp(ip)*z(j)))

      if((areai.lt.0.0d+0).or.(areaj.lt.0.0d+0)
     $ .or.(areak.lt.0.0d+0)) go to 399
      totarea = areai + areaj + areak
      resid = abs(area(m) - totarea)
      if (resid .le. 0.001) then
       ipflag(ip)=2
       elhom(ip) = m
      else
      endif
  399 continue
      if(ipflag(ip).gt.0) go to 200
!
! if were here, the particle hasn't been found in any
! element. it must be outside the solution domain.
! try to find it using njn (boundary node) array
!

       do 1700 nn=1,njb
       if(ip.eq.njn(nn)) then
       m2 = njn(nn)
       if(nn-1.eq.0) then
       m1 = njn(tnrows(1))
       else
       m1 = njn(nn-1)
       endif
       if(nn+1.gt.njb) then
       m3 = nnode
       else
       m3 = njn(nn+1)
       endif
       if(nn+1.eq.tnrows(1)) then
       m3 = tnrows(1)
       else
       endif
       endif
 1700 continue
      if(m1.le.0.or.m2.le.0) go to 232
      if(m1.gt.nnode.or.m2.gt.nnode) go to 232
      l1 = dsqrt( (x(m3)-x(m2))**2.0 +  (z(m3)-z(m2))**2.0 )
      c1 = dsqrt( (x(m3)-xp(ip))**2.0 +  (z(m3)-zp(ip))**2.0 )
      e1 = dsqrt( (x(m2)-xp(ip))**2.0 +  (z(m2)-zp(ip))**2.0 )
      d1 = (e1*e1+l1*l1-c1*c1)/(2.0*l1)
      a1 = l1-d1
      l1test = dabs(a1)+dabs(d1)
      diff1 = dabs(l1-l1test)
      if(diff1.lt.0.001) then
      xp(ip) = (x(m2) - x(m3))*(a1/l1) + x(m3)
      zp(ip) = (z(m2) - z(m3))*(a1/l1) + z(m3)
        do 701 n=1,numpnb(ip)
       m = ielneb(ip,n)
       if(m.le.0.or.m.gt.nelem) go to 701
       if(ni(m).eq.0) go to 701
        i = ni(m)
        j = nj(m)
        k = nk(m)
       if(m2.eq.i.and.m3.eq.j) elhom(ip) = m
       if(m2.eq.i.and.m3.eq.k) elhom(ip) = m
       if(m2.eq.j.and.m3.eq.i) elhom(ip) = m
       if(m2.eq.j.and.m3.eq.k) elhom(ip) = m
       if(m2.eq.k.and.m3.eq.i) elhom(ip) = m
       if(m2.eq.k.and.m3.eq.j) elhom(ip) = m
       if(m2.eq.i.and.m3.eq.j.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.i.and.m3.eq.k.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.j.and.m3.eq.i.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.j.and.m3.eq.k.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.k.and.m3.eq.i.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.k.and.m3.eq.j.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.i.and.m3.eq.j.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.i.and.m3.eq.k.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.j.and.m3.eq.i.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.j.and.m3.eq.k.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.k.and.m3.eq.i.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.k.and.m3.eq.j.and.ipflag(ip).eq.3)ipflag(ip) = 4

  701 continue
      endif
       if(ipflag(ip).gt.0) go to 200
      l2 = dsqrt( (x(m1)-x(m2))**2.0 +  (z(m1)-z(m2))**2.0 )
      e2 = dsqrt( (x(m1)-xp(ip))**2.0 +  (z(m1)-zp(ip))**2.0 )
      c2 = dsqrt( (x(m2)-xp(ip))**2.0 +  (z(m2)-zp(ip))**2.0 )
      d2 = (e2*e2+l2*l2-c2*c2)/(2.0*l2)
      a2 = l2-d2
      l2test = dabs(a2)+dabs(d2)
      diff2 = dabs(l2-l2test)
      if(diff2.lt.0.001) then
      xp(ip) = (x(m1) - x(m2))*(a2/l2) + x(m2)
      zp(ip) = (z(m1) - z(m2))*(a2/l2) + z(m2)

        do 702 n=1,numpnb(ip)
       m = ielneb(ip,n)
       if(m.le.0.or.m.gt.nelem) go to 702
       if(ni(m).eq.0) go to 702
       if(ipflag(ip).eq.3) go to 702
        i = ni(m)
        j = nj(m)
        k = nk(m)
       if(m2.eq.i.and.m1.eq.j) elhom(ip) = m
       if(m2.eq.i.and.m1.eq.k) elhom(ip) = m
       if(m2.eq.j.and.m1.eq.i) elhom(ip) = m
       if(m2.eq.j.and.m1.eq.k) elhom(ip) = m
       if(m2.eq.k.and.m1.eq.i) elhom(ip) = m
       if(m2.eq.k.and.m1.eq.j) elhom(ip) = m
       if(m2.eq.i.and.m1.eq.j.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.i.and.m1.eq.k.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.j.and.m1.eq.i.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.j.and.m1.eq.k.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.k.and.m1.eq.i.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.k.and.m1.eq.j.and.ipflag(ip).eq.0) ipflag(ip) = 3
       if(m2.eq.i.and.m1.eq.j.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.i.and.m1.eq.k.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.j.and.m1.eq.i.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.j.and.m1.eq.k.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.k.and.m1.eq.i.and.ipflag(ip).eq.3)ipflag(ip) = 4
       if(m2.eq.k.and.m1.eq.j.and.ipflag(ip).eq.3)ipflag(ip) = 4
  702 continue
      endif
      if(ipflag(ip).gt.0) go to 200
!
! since weve corrected the position of the particle, send it
! back to the top of the loop to find out if it resides in the
! surrounding elements
!
      if(ipflag(ip).eq.3) go to 115
      if(ipflag(ip).gt.0) go to 200
! if were here, the particle hasnt been found in any
! element. it must be outside the solution domain. we need
! to reduce the magnitude of the nodal velocity vector
! until the partile appears inside the solution domain
!
  232 continue
       velfrac = 1.0d+0
       do 700 iii=1,10
       velfrac = velfrac - 0.1
       xp(ip) = x(ip) - vxp(ip)*dtp*velfrac
       zp(ip) = z(ip) - vzp(ip)*dtp*velfrac
       do 1300 n=1,numpnb(ip)
       m = ielneb(ip,n)
       if(ni(m).eq.0) go to 1300
       if(ipflag(ip).gt.0) go to 1300

        i = ni(m)
        j = nj(m)
        k = nk(m)

        areai = 0.5d+00 * ((xp(ip)*z(j) - x(j)*zp(ip)) +
     $       (x(k)*zp(ip) - xp(ip)*z(k)) + (x(j)*z(k) - x(k)*z(j)))

        areaj = 0.5d+00 * ((x(i)*zp(ip) - xp(ip)*z(i)) + (x(k)*z(i) -
     $          x(i)*z(k)) + (xp(ip)*z(k) - x(k)*zp(ip)))

        areak = 0.5d+00 * ((x(i)*z(j) - x(j)*z(i)) + (xp(ip)*z(i) -
     $          x(i)*zp(ip)) + (x(j)*zp(ip) - xp(ip)*z(j)))

      if((areai.lt.0.0d+0).or.(areaj.lt.0.0d+0)
     $ .or.(areak.lt.0.0d+0)) go to 1300
      totarea = areai + areaj + areak
      resid = abs(area(m) - totarea)
      if (resid .le. 0.001) then
!
! okay, we have found which element the particle resides in
! now determine its concentration using shape function weighting
! scheme
!
       ipflag(ip) = 5
       elhom(ip) = m
      endif
 1300 continue
  700 continue

   15 continue
  200 continue
!
! check to see if particle is still lost. if it is
! then assign it the nodal concentration
!
      do 777 ip=1,nump
      if(ipflag(ip).eq.0) then
      write (8,101) ip,x(ip),xp(ip),z(ip),zp(ip)
      concp(ip) = conc(ip)
      else
      endif
  777 continue
!**********************************************************************
! format statements
  101 format ('unable to locate particle',i5,2x,'x=',
     $        1pe12.4,'     xp=',1pe12.4,'     z=',1pe12.4,
     $ '     zp=',1pe12.4)
      return
      end
!**********************************************************************
      subroutine parinterp (x,z,ni,nj,nk,area,al,b,c,elhom,conc
     $  ,concp,xp,zp,nump,nc,ncb,icbflg,tncols,tnrows,ntp
     $ ,numpnb,ielneb,njn,njb,nelem)
!
!
!   this subroutine calculates particle concentrations
!   using nodal concentrations
!   assuming three node triangular elements
!
!
!**********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 area,x,z,xc,zc,xp,zp,conc,concp
      real*8 shapi,shapj,shapk

      integer ni,nj,nk,nump,elhom(maxnodes),nc(maxbnds),ncb
      integer ip,i,j,k,icbflg,tncols,tnrows,ntp
      integer numpnb(maxnodes),ielneb(maxnodes,12)
      integer njb,njn(maxbnds),nelem

      dimension conc(maxnodes),area(maxelems),ni(maxelems),nj(maxelems)
      dimension x(maxnodes),z(maxnodes),al(6,maxelems),b(6,maxelems)
      dimension concp(maxnodes),xp(maxnodes),zp(maxnodes),c(6,maxelems)
      dimension nk(maxelems),tnrows(maxbnds),icbflg(maxnodes)
!
! now find out what element the particle has moved into
!
       do 100 ip=1,nump
       m = elhom(ip)
       if(m.gt.nelem.or.m.le.0) concp(ip) = conc(ip)
       if(m.gt.nelem.or.m.le.0) go to 100
       if(icbflg(ip).eq.1) go to 100
       if(ni(m).eq.0) go to 100

        i = ni(m)
        j = nj(m)
        k = nk(m)
!
! determine its concentration using shape function weighting
! scheme
!
       xc=xp(ip)
       zc=zp(ip)
       shapi=0.5d+0*(al(1,m)+b(1,m)*xc+c(1,m)*zc)/area(m)
       shapj=0.5d+0*(al(2,m)+b(2,m)*xc+c(2,m)*zc)/area(m)
       shapk=0.5d+0*(al(3,m)+b(3,m)*xc+c(3,m)*zc)/area(m)
       concp(ip) = shapi*conc(i)+shapj*conc(j)+shapk*conc(k)
       if(concp(ip).lt.0.0) concp(ip)=0.0d+0
  100  continue
!
! update boundary condition concetrations
!
      do 777 n=1,ncb
  777 concp(nc(n))=conc(nc(n))
!
! for nodes within thin (delz<delzmn) elements just below
! the land surface, set particle concentration equalt to
! boundary concentration
!
      nn = 0
      nnn = 0
      do 102 n=1,tncols
      nnn=nnn+tnrows(n)
      do 202 m=1,tnrows(n)
      nn=nn+1
      if(icbflg(nn).eq.1) concp(nn) = conc(nc(n))
 202  continue
 102  continue
       return
       end
!***********************************************************************
      subroutine pargen (xp,zp,nelem,nnode,nump,x,z,ni,nj,nk,nl,
     $ phi,qx,qz,vxp,vzp,numpnb,ielneb,area,nelem1,voxp,vozp,qxo,qzo)

!
!  this subroutine generates particles in triangular elements.
!  particles are distributed uniformly throughout each element
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 x,z,xp,zp,qx,qz,qxo,qzo,voxp,vozp

      integer ni,nj,nk,nl,nelem,nnode
      integer nump,i,j,k,l,ip,numpnb,ielneb(maxnodes,12)

      dimension x(maxnodes),z(maxnodes),numpnb(maxnodes)
      dimension xp(maxnodes),zp(maxnodes),area(maxelems)
      dimension ni(maxelems),nj(maxelems),nk(maxelems),nl(maxelems)
      dimension vxp(maxnodes),vzp(maxnodes),artot(maxnodes)
      dimension qx(maxelems),qz(maxelems),phi(maxnodes),voxp(maxnodes)
      dimension qxo(maxelems),qzo(maxelems),vozp(maxnodes)

      nump = nnode
!
! zero out arrays for current time step
!
      do 100, ip=1,nump
      vxp(ip) = 0.0d+0
      vzp(ip) = 0.0d+0
      voxp(ip) = 0.0d+0
      vozp(ip) = 0.0d+0
      artot(ip) = 0.0d+0
      numpnb(ip) = 0
  100 continue

      nelem1=0
      do 200 m = 1,nelem
!
! determine number of elements that share a common node
!
       i = ni(m)
       j = nj(m)
       k = nk(m)
       l = nl(m)
       if(i.eq.0) go to 200
       nelem1=nelem1+1
       numpnb(i) = numpnb(i) + 1
       numpnb(j) = numpnb(j) + 1
       numpnb(k) = numpnb(k) + 1
       if(l.gt.0) numpnb(l) = numpnb(l) + 1
       ielneb(i,numpnb(i)) = m
       ielneb(j,numpnb(j)) = m
       ielneb(k,numpnb(k)) = m
       if(l.gt.0) ielneb(l,numpnb(l)) = m
!
!  place particles at nodes
!
       xp(i) = x(i)
       zp(i) = z(i)
       xp(j) = x(j)
       zp(j) = z(j)
       xp(k) = x(k)
       zp(k) = z(k)
       if(l.gt.0) xp(l) = x(l)
       if(l.gt.0) zp(l) = z(l)
!
! calculate weighted nodal velocities
!
      if(phi(i).le.0.0) go to 200
      vxp(i) = qx(m)*area(m)/phi(i) + vxp(i)
      vzp(i) = qz(m)*area(m)/phi(i) + vzp(i)
      voxp(i) = qxo(m)*area(m)/phi(i) + voxp(i)
      vozp(i) = qzo(m)*area(m)/phi(i) + vozp(i)
      artot(i) = area(m) + artot(i)

      if(phi(j).le.0.0) go to 200
      vxp(j) = qx(m)*area(m)/phi(j) + vxp(j)
      vzp(j) = qz(m)*area(m)/phi(j) + vzp(j)
      voxp(j) = qxo(m)*area(m)/phi(j) + voxp(j)
      vozp(j) = qzo(m)*area(m)/phi(j) + vozp(j)
      artot(j) = area(m) + artot(j)

      if(phi(k).le.0.0) go to 200
      vxp(k) = qx(m)*area(m)/phi(k) + vxp(k)
      vzp(k) = qz(m)*area(m)/phi(k) + vzp(k)
      voxp(k) = qxo(m)*area(m)/phi(k) + voxp(k)
      vozp(k) = qzo(m)*area(m)/phi(k) + vozp(k)
      artot(k) = area(m) + artot(k)

      if(l.eq.0.or.phi(l).le.0.0) go to 200
      if(l.gt.0) vxp(l) = qx(m)*area(m)/phi(l) + vxp(l)
      if(l.gt.0) vzp(l) = qz(m)*area(m)/phi(l) + vzp(l)

      if(l.gt.0) voxp(l) = qxo(m)*area(m)/phi(l) + voxp(l)
      if(l.gt.0) vozp(l) = qzo(m)*area(m)/phi(l) + vozp(l)

      if(l.gt.0) artot(l) = area(m) + artot(l)
  200 continue
!
! account for differences in element area in weighted velocity
!
      do 300, ip=1,nump
      if(artot(ip).gt.0.0) then
      vxp(ip) = vxp(ip)/artot(ip)
      vzp(ip) = vzp(ip)/artot(ip)
      voxp(ip) = voxp(ip)/artot(ip)
      vozp(ip) = vozp(ip)/artot(ip)
      else
      vxp(ip) = 0.0
      vzp(ip) = 0.0
      voxp(ip) = 0.0
      vozp(ip) = 0.0
      endif

  300 continue
  301 format ('node=',i5,'num neb:',i4)
  302 format ('element neb:',12(i7,1x))
      return
      end
!***********************************************************************
      subroutine tstep (vxp,vzp,dtp,nelem,b,c,ni,nj,nk,dt,ntimp,dtpr,
     $ gamma,delzmn,icbflg,nnode,maxit,total_time,ntime)
!
!  this subroutine calculates the time step to be used
!  for the solute transport algorithm
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 total_time,delz,delx,deltz,deltx
      real*8 dtp,vxp,vzp,gamma,delzmn

      integer ni,nj,nk,nelem,nnode,icbflg,maxit

      dimension vxp(maxnodes),vzp(maxnodes),icbflg(maxnodes)
      dimension b(6,maxelems),c(6,maxelems)
      dimension ni(maxelems),nk(maxelems),nj(maxelems)
!
!  set up some initial parameter values
!
       dtp = 1.d+8
      do 455 n=1,nnode
      icbflg(n) = 0
  455 continue
!
!  determine maximum time step using element
!  velocities and characteristic distances of
!  each elements (ie. b ~ delz & c ~ delx)
!
        do 20 l = 1,nelem
            if(ni(l).eq.0) go to 20
!
!  first determine maximum dimensions of element in
!  x & z directions using derivatives of shape functions
!
       delx = dabs(c(1,l))
       if(dabs(c(2,l)).gt.delx) delx = dabs(c(2,l))
       if(dabs(c(3,l)).gt.delx) delx = dabs(c(3,l))
       delz = dabs(b(1,l))
       if(dabs(b(2,l)).gt.delz) delz = dabs(b(2,l))
       if(dabs(b(3,l)).gt.delz) delz = dabs(b(3,l))
       if(delz.le.delzmn)then
        icbflg(ni(l)) = 1
        icbflg(nj(l)) = 1
        icbflg(nk(l)) = 1
        go to 20
        else
        continue
        endif
!
!      calculate element velocity
!
       vxe = (vxp(ni(l))+ vxp(nj(l)) + vxp(nk(l)))/3.d+0
       vze = (vzp(ni(l))+ vzp(nj(l)) + vzp(nk(l)))/3.d+0
!
!  calculate delt for the element and
!  determine minimum elemental time step for
!  particle tracking
!
       if(vxe.eq.0.0d+0) go to 21
       deltx = gamma*delx/dabs(vxe)
       if(deltx.lt.dtp) dtp = deltx
   21 if(vze.eq.0.0d+0) go to 20
      deltz = gamma*delz/dabs(vze)
      if(deltx.le.0.0d+0.or.deltz.le.0.0d+0) go to 20
      if(deltz.lt.dtp) dtp = deltz
   20 continue
!
! calculate number of time steps required for the solute transport
!  loop
!
      if(dtp.gt.dt) dtp = dt/1.1
      ntimp = dint(dt/dtp)
      atemp = dble(ntimp)*dtp
      dtpr = dmod(dt,atemp)
      write(8,200) dtp,dtpr,ntimp,maxit
      write(*,201) dtp,dtpr
      if(ntimp.gt.maxit) ntimp=maxit
!***********************************************************************
!  format statement
!
  200 format (3x,'  dtp=',1pe12.4,3x,/,3x,
     $ '  dtpr=',1pe12.4,3x,/,'   ntimp=',i7,3x,/,'   maxit=',i7)
  201 format (3x,'  dtp=',1pe12.4,3x,/,3x,
     $ '  dtpr=',1pe12.4)
      return
      end
!***********************************************************************
      subroutine brine(aa,bb,x,z,ni,nj,nk,nl,dif,ldis,
     $ tdis,nelem,iout,mat,area,al,b,c,dtp,nfband,
     $ concp,theta,eltp,qx,qz,nnode,conc,it,rotflag,temp)
!
!     this subroutine computes the stiffness matrix and forcing
!     vectors for the brine transport  equation using three and four node
!     triangular elements using the modified method of characteristics.
!
!
!***********************************************************************

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 g,b,c,be,p,theta,x,z,area,dtp,kxe,kze,kxze
      real*8 c1,c2,c3,al,dd,ddfn,aa,bb,a,ae,qbar
      real*8 dif,dife,qx,qz,qx2,qz2,ldis,tdis,concp,conc
      real*8 temp(maxnodes),kint,tbar,cbar,gam,e1,ao

      integer nnode,nelem,eltp,i1,i2,i3,jj,ii,mat
      integer ni,nj,nk,nl,qq,iout,nn,ne,me,sw
      integer node,nfband,itf,i,j,k,l,m,rotflag

      dimension aa(maxnodes,maxwdt),bb(maxnodes),be(8),ddfn(8)
      dimension al(6,maxelems),b(6,maxelems),c(6,maxelems),mat(maxelems)
      dimension area(maxelems),ni(maxelems),nj(maxelems),nk(maxelems)
      dimension x(maxnodes),z(maxnodes),p(4,4),g(4,4),node(8)
      dimension dif(0:maxmats),qx(maxelems),qz(maxelems)
      dimension ldis(0:maxmats),tdis(0:maxmats),concp(maxnodes)
      dimension nl(maxelems),eltp(maxelems),conc(maxnodes)

      do 100 m = 1, nelem

      do 1402 lp=1,4
      do 1401 mp=1,4
 1401 g(lp,mp)=0.0
      be(lp)=0.0
 1402 continue

        if(eltp(m).eq.6)goto100
        ae = area(m)
        k = nk(m)
        i = ni(m)
        j = nj(m)
        k = nk(m)
        l = nl(m)
        if(i.eq.0) go to 100
        qq=mat(m)
        dife=dif(qq)
        qx2 = qx(m)*qx(m)
        qz2 = qz(m)*qz(m)
        qbar = dsqrt(qx2+qz2)
        if(qbar.eq.0.0d+00) go to 125
        if(rotflag.eq.0) then
        kxe = ((ldis(qq)*qx2+tdis(qq)*qz2))/qbar+dife
        kze = ((tdis(qq)*qx2+ldis(qq)*qz2))/qbar+dife
        kxze = (ldis(qq)-tdis(qq))*qx(m)*qz(m)/qbar
        else
        kxe = dife
        kze = dife
        kxze = 0.0d+00
        endif
        goto 135
  125 continue
        kxe = dife
        kze = dife
        kxze = 0.0d+00
  135 continue
        a = 1.d+0/ae
        d = ae/dtp
        if(eltp(m).eq.1) go to 300
!
!       determines which terms go with which nodes for the stiffness
!       and capacitence matrices, and which beta and gama terms are
!       used by the terms in the stiffnes matrix
!
          sw=4
          if(eltp(m).eq.2)then
            i1 = 1
            i2 = 2
            i3 = 3
          end if
          if((eltp(m).eq.4).or.(eltp(m).eq.5))then
            i1 = 2
            i2 = 3
            i3 = 1
            itf=i
            i=j
            j=k
            k=itf
          end if
          if(eltp(m).eq.3)then
            i1 = 3
            i2 = 2
            i3 = 1
            itf=i
            i=k
            k=itf
          end if
          node(1)=i
          node(2)=j
          node(3)=k
          node(4)=l
          b1=b(i1,m)
          b2=b(i2,m)
          b3=b(i3,m)
          g1=c(i1,m)
          g2=c(i2,m)
          g3=c(i3,m)
!
!     determination of values of the constants
!
      c2 = (al(i2,m)+b2*x(l)+g2*z(l))/(2.0d+00*ae)
      c3 = (al(i3,m)+b3*x(l)+g3*z(l))/(2.0d+00*ae)
      if(c2*c3.eq.0.00) go to 301
      c1 = 1.0d+00/(c3*c2)
!
!     calculate stiffness matrix terms
!
      g(1,1)=a*(((kxe*b1**2)+(kze*g1**2))/4.0d+00)
      g(1,2)=a*((kxe*((b1/4.0d+00)*(b2-(c1*c2/3.0d+00)*(b2+b3)))+(kze*
     $((g1/4.0d+00)*(g2-(c1*c2/3.0d+00)*(g2+g3))))))
      g(1,3)=a*((kxe*((b1/4.0d+00)*(b3-(c1*c3/3.0d+00)*(b2+b3)))+(kze*
     $((g1/4.0d+00)*(g3-(c1*c3/3.0d+00)*(g2+g3))))))
      g(1,4)=a*((kxe*((b1/12.0d+00)*c1*(b3+b2)))+(kze*((g1/12.0d+00)*
     $c1*(g3+g2))))
      g(2,1)=g(1,2)
      g(2,2)=a*(0.5d+00)*(((kxe*(b2**2)+kze*(g2**2))/2.0d+00)+
     $((c1*c2/3.0d+00)*(-kxe*(b2*b3+b2**2)-kze*(g2**2+g2*g3)+((c1*c2/
     $4.0d+00)*((kxe*(b3**2+b2*b3+b2**2)+kze*(g2**2+g3**2+g2*g3)))))))
      g(2,3)=a*(0.25d+00)*(kxe*((b2*b3))+(kze*(g2*g3))-(c1*c2/3.0d+00)*
     $(kxe*((b3**2+b2*b3))+(kze*(g2*g3+g3**2)))-(c1*c3/3.0d+00)*
     $(kxe*((b2**2+b2*b3))+(kze*(g2*g3+g2**2)))+(((c1**2)*c2*c3)/6.0d+
     $00)*(kxe*((b2**2+b2*b3+b3**2))+(kze*(g3**2+g2**2+g2*g3))))
      g(2,4)=a*((c1/12.0d+00)*(kxe*((b2*b3+b2**2))+(kze*(g2**2+g2*g3))-
     $(c1*c2/2.0d+00)*(kxe*((b2**2+b2*b3+b3**2))+(kze*(g3**2+g2**2+g2*
     $g3)))))
      g(3,1)=g(1,3)
      g(3,2)=g(2,3)
      g(3,3)=a*(0.5d+00)*(((kxe*(b3**2)+kze*(g3**2))/2.0d+00)+((c1*c3
     $/3.0d+00)*(-kxe*(b2*b3+b3**2)-kze*(g3**2+g2*g3)+((c1*c3/4.0d+00)
     $*((kxe*(b3**2+b2*b3+b2**2)+kze*(g2**2+g3**2+g2*g3)))))))
      g(3,4)=a*((c1/12.0d+00)*((kxe*(b2*b3+b3**2))+(kze*(g3**2+g2*g3))-
     $(c1*c3/2.0d+00)*(kxe*((b2**2+b2*b3+b3**2))+(kze*(g3**2+g2**2+g2*
     $g3)))))
      g(4,1)=g(1,4)
      g(4,2)=g(2,4)
      g(4,3)=g(3,4)
      g(4,4)=a*(((c1**2/24.0d+00)*((kxe*(b2**2+b2*b3+b3**2))+(kze*
     $(g2*g3+g3**2+g2**2)))))
!
! compute full capacitence matrix terms
!
      p(1,1) = d/6.0d+00
      p(1,2) = d/12.0d+00-c1*c2*d/60.0d+00
      p(1,3) = d/12.0d+00-c1*c3*d/60.0d+00
      p(1,4) = d*c1/60.0d+00
      p(2,1) = p(1,2)
      p(2,2) = d/6.0d+00-c2*c1*d/15.0d+00+(c2**2)*(c1**2)*d/90.0d+00
      p(2,3) = d/12.0d+00-c3*c1*d/30.0d+00-c2*c1*d/30.0d+00+c2*c3*
     $(c1**2)*d/90.0d+00
      p(2,4) = d*c1/30.0d+00-c2*(c1**2)*d/90.0d+00
      p(3,1) = p(1,3)
      p(3,2) = p(2,3)
      p(3,3) = d/6.0d+00-c3*c1*d/15.0d+00+(c3**2)*(c1**2)*d/90.0d+00
      p(3,4) = d*c1/30.0d+00-c3*(c1**2)*d/90.0d+00
      p(4,1) = p(1,4)
      p(4,2) = p(2,4)
      p(4,3) = p(3,4)
      p(4,4) = d*(c1**2)/90.0d+00

! create lumped capacitance matrix

      p(1,1) = p(1,1)+p(1,2)+p(1,3)+p(1,4)
      p(1,2) = 0.0d+00
      p(1,3) = 0.0d+00
      p(1,4) = 0.0d+00
      p(2,2) = p(2,1)+p(2,2)+p(2,3)+p(2,4)
      p(2,1) = 0.0d+00
      p(2,3) = 0.0d+00
      p(2,4) = 0.0d+00
      p(3,3) = p(3,1)+p(3,2)+p(3,3)+p(3,4)
      p(3,1) = 0.0d+00
      p(3,2) = 0.0d+00
      p(3,4) = 0.0d+00
      p(4,4) = p(4,1)+p(4,2)+p(4,3)+p(4,4)
      p(4,1) = 0.0d+00
      p(4,2) = 0.0d+00
      p(4,3) = 0.0d+00

!
!     detrmination of the weighted stiffness matrix
!
      do 200 n=1,4
        ddfn(n)=0.d+00
        do 250 nn=1,4
          ddfn(n) = ddfn(n) + ((1.d+00-theta)*g(n,nn) -
     $            p(n,nn))*concp(node(nn))

          g(n,nn) = theta*g(n,nn)+p(n,nn)

  250 continue

      be(n) = -ddfn(n)

  200 continue

  300 if(eltp(m).ne.1) go to 700
  301 continue
!
! calculate three node stiffnes matrix
!
      node(1)=ni(m)
      node(2)=nj(m)
      node(3)=nk(m)

      ao = 0.002
      eo = 1000.0
      cmax = 0.06
      tbar = (temp(ni(m))+temp(nj(m))+temp(nk(m)))/3.0
      kint = ao*exp(-eo/tbar)
      cbar = (concp(ni(m))+concp(nj(m))+concp(nk(m)))/3.0
      gam = kint*(cmax - cbar)
                                            

      do 400 n=1,3
        dd=0.0d+00
        do 500 nn=1,3
          if(n.eq.nn)then
             p(n,nn)=d/3.0d+00
          else if (n.ne.nn)then
             p(n,nn)=0.0d+00
          end if

          g(n,nn)=(a/4.0d+00)*(kxe*b(nn,m)*b(n,m)+kze*c(nn,m)*c(n,m)+
     $    kxze*(c(nn,m)*b(n,m)+b(nn,m)*c(n,m)))

          dd = dd + ((1.0d+00-theta)*g(n,nn)-p(n,nn))*concp(node(nn))

           g(n,nn) = theta*g(n,nn)+p(n,nn)

  500   continue

! add source term to crystalline basement

        if(qq.le.1) then
        be(n) =  -dd + dtp*gam*ae/3.0
        else
        be(n) =  -dd
        end if

  400 continue

!
!     form global a matrix and b vector
!
      sw=3

  700  continue

      if(iout.ne.2)  go to 401
      write(8,12) m,i,j,k,l
      do 402 lp=1,4
  402 write (8,21) (g(lp,mp),mp=1,4)
      write (8,13) (be(lp),lp=1,4)
  401 continue

        do 800 ii = 1,sw
          go to (35,40,45,50),ii
   35     ne = i
          go to 55
   40     ne = j
          go to 55
   45     ne = k
          go to 55
   50     ne = l
          go to 55

   55     do 900 jj=1,sw
            go to (60,65,70,75),jj
   60       me = i
            go to 80
   65       me = j
            go to 80
   70       me = k
            go to 80
   75       me = l
   80       kl=me-ne+1
            if(kl.le.0)go to 900
            aa(ne,kl)=aa(ne,kl)+g(ii,jj)
  900     continue
          bb(ne)=bb(ne)+be(ii)
  800   continue


  100 continue

!
!     print out stiffness matrix
!
      if (iout.ne.2) go to 902
      write(8,11)
      do 901 ll=1,nnode
      write(8,*) 'll=',ll
      write (8,21) (aa(ll,mm),mm=1,nfband)
  901 continue
      write (8,31) (bb(ll),ll=1,nnode)
  902 continue
!**********************************************************************
! format statements

   11 format (/1x,'brine stiffness matrix:')
   12 format (' local fstiff for elem:',i4,'i,j,k,l:',4(i5,1x),//)
   13 format (/,3x,'b vector:',4(1pe12.3))
   21 format (118(1pe12.3))
   31 format (//,1x,'b vector',10(1pe12.3))
      return
      end
!**********************************************************************
!  gaussian elimination for a symmetric,positive definite band matrix
!
!
      subroutine bgauss(aa,bb,nnode,nfband,conc,nc,ncb)
!
!   this subroutine solves solute transport aa matrix
!   and b vectore using gaussian elimination
!
! *********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8  aa(maxnodes,maxwdt),bb(maxnodes),conc(maxnodes)

      integer ncb,nnode,nfband,nc(maxbnds)
!
!     triangularization
!
      do 30 n=1,nnode
      do 20 l=2,nfband
      if(aa(n,l).eq.0.0)go to 20
      c=aa(n,l)/aa(n,1)
      i=n+l-1
      if(i.gt.nnode)go to 20
      j=0
      do 10 k=l,nfband
      j=j+1
   10 if(aa(n,k).ne.0.0)aa(i,j)=aa(i,j)-c*aa(n,k)
      aa(n,l)=c
      bb(i)=bb(i)-aa(n,l)*bb(n)
   20 continue
      bb(n)=bb(n)/aa(n,1)
   30 continue
!
! back substitution
!
      n=nnode
   40 do 50 k=2,nfband
      l=n+k-1
      if(l.gt.nnode)go to 60
      if(aa(n,k).ne.0.0)bb(n)=bb(n)-aa(n,k)*bb(l)
   50 continue
   60 n=n-1
      if(n.gt.0)go to 40
!
!     update concentrations
!
      do 700 i=1,nnode
      do 600 n=1,ncb
      if(i.eq.nc(n)) go to 700
  600 continue
      conc(i)=bb(i)
  700 continue

      return
      end
!***********************************************************************
      subroutine bebnd(aa,bb,conc,nc,ncb,nnode,iout,nfband)
!
!     this subroutine adjusts the aa matrix and bb vector to
!     account for constant concentration nodes
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 aa,bb,conc

      integer nc,ncb,nnode,iout,nfband

      dimension  aa(maxnodes,maxwdt),bb(maxnodes)
      dimension conc(maxnodes),nc(maxbnds)
!
!   alter aa and bb for constant concentration nodes, keeping symmetry
!
!
      do 100 n=1,ncb
      bb(nc(n))=conc(nc(n))
      aa(nc(n),1)=1.0d0

      do 140 l=2,nfband
      jl=nc(n)+l-1
      if(jl.gt.nnode)go to 130
      bb(jl)=bb(jl)-aa(nc(n),l)*bb(nc(n))
  130 aa(nc(n),l)=0.0d0
      jk=nc(n)-l+1
      if(jk.lt.1)go to 140
      bb(jk)=bb(jk)-aa(jk,l)*bb(nc(n))
      aa(jk,l)=0.0d0
  140 continue
  100 continue
!
!     print out aa matrix
!
      if (iout.ne.2) go to 300
      write(8,10)
      do 250 l=1,40
      write (8,20) (aa(l,m),m=1,20)
  250 continue
!
!     print out b vector
!

      write (8,30) (bb(l),l=1,40)
  300 continue
!
!     okay were done, back from wence we came ...
!
   10 format (/,1x,'subroutine bebnd',/,1x,'aa matrix:')
   20 format (11(1pe12.3))
   30 format (//,1x,'b vector',10(1pe12.3))
      return
      end
!***********************************************************************
      subroutine zero1 (aa,bb,nnode,nband)
!
!
!     this subroutine zero's the aa and bb array's and takes
!     care of some initial conditions
!
!
!************************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 aa(maxnodes,maxwdt),bb(maxnodes)

      integer nnode,nband

      do 100 l=1,nnode
      do 200 m=1,nband
      aa(l,m)=0.0d+00
  200 continue
      bb(l)=0.0d+00
  100 continue
  400 return
      end
!******************************************************************
!			subroutine porperm
!
       subroutine porperm (perm,rhos,rhof,phi,icase,
     $ nmat,matp,phi_o,beta,phi_ir,beta_ul,perm1,perm2,sigmae,
     $ phiold,ss,sigmax,tnrows,tncols,phimin,presl,it,z,pres,vzt,
     $ p_flag,dgrn,iperm,presh,iter)
!
!
!     this subroutine computes permeability and porostiy as a
!     function of effective stress. two relationships are used
!     depending on whether loading or unloading is occurring.
!
!******************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
       real*8 depth,z(maxnodes),phi(maxnodes),perm(maxnodes)
       real*8 phi_o(0:maxmats),rhos(0:maxmats),rhof(maxnodes)
       real*8 sigmae(maxnodes),atemp,perm1(0:maxmats),perm2(0:maxmats)
       real*8 phi_ir(0:maxmats),sigmax(maxnodes),phimin(maxnodes),a1,a2
       real*8 vzt(maxbnds),dgrn(0:maxmats),presh(maxnodes)
       real*8 presl(maxnodes),pres(maxnodes),beta(0:maxmats)
       real*8 beta_ul(0:maxmats),ss(maxnodes),phiold(maxnodes),g

       integer n,icase(maxbnds),matp(maxnodes),q,tnrows(maxbnds),tncols
       integer p_flag,it,iter

      units = 3.09d+14
      g=9.81d+0
      nn = 0
      nnn = 0
      do 100 n=1,tncols
      nnn=nn+tnrows(n)
      do 200 m=1,tnrows(n)
       nn=nn+1
       q = matp(nn)
       if(iter.eq.1) phiold(nn)=phi(nn)
       depth = z(nnn)-z(nn)
       presh(nn)=depth*g*rhof(nn)
	
            if(sigmae(nn).gt.sigmax(nn)) sigmax(nn)=sigmae(nn)
	    if(phimin(nn).gt.phi(nn))    phimin(nn)=phi(nn)
!
! use the basic porosity equation during deposition, uplift or a hiatus:
!
	if(icase(n).ne.3)then
         atemp = -beta(q)*sigmae(nn)/(g*(rhos(q)-rhof(nn)))
         phi(nn) = phi_o(q)*dexp(atemp) + phi_ir(q)
        else
         a1 = phimin(nn) - phi_ir(q)
         a2 =  sigmae(nn) - sigmax(nn)
         if(a2.lt.0) a2 = 0.0
         atemp = -beta_ul(q)*(a2)/(g*(rhos(q)-rhof(nn)))
         phi(nn) = a1*dexp(atemp) + phi_ir(q)
!         write(8,*) 'nn,phi,phimin,a1,a2',nn,phi(nn),phimin(nn),a1,a2
! project porosity at depth to land surface if erroding
         if(nn.eq.nnn) then
         dlzz =  z(nn) - z(nn-1)
         dlz =   z(nn-2) - z(nn-1)
         phi(nn) = phi(nn-1) - dlzz*(phi(nn-1)-phi(nn-2))/dlz
         endif
         endif

       if (iperm.eq.1) then
       perm(nn) = 10.d+0**(perm2(q)*phi(nn)+perm1(q))
       perm(nn) = perm(nn)*units
       endif
       if (iperm.eq.2) then
       perm(nn) = (dgrn(q)*dgrn(q)*phi(nn)*phi(nn)*phi(nn))
     $ /(180.d+0*(1.0d+0-phi(nn))*(1.0d+0-phi(nn)))
       perm(nn) = perm(nn)*units 
       endif    

	if(icase(n).ne.3)then   
           ss(nn) = beta(q)*phi(nn)
	else
           ss(nn) = beta_ul(q)*phi(nn)
	endif
    9  format(i5,2x,3(1pe12.3,2x))
  200  continue
  100  continue

 1112   format ('nn=',i5,2x,'phi=',1pe12.4,2x,
     $  'phimin=',1pe12.4,2x,'sigmae=',
     $  1pe12.4,2x,'sigmax=',1pe12.4,'beta_ul=',1pe12.4,2x,
     $ 'phi_ir=',1pe12.4,2x,'atemp=',1pe12.4,2x,'a1=',1pe12.4,
     $ 2x,'a2=',1pe12.4)
 1113   format ('nn=',i5,2x,'phi=',1pe12.4,2x,'phimin=',
     $   1pe12.4,2x,'sigmae=',
     $   1pe12.4,2x,'sigmax=',1pe12.4,'beta=',1pe12.4,2x,
     $   'phi_o=',1pe12.4,2x,'atemp=',1pe12.4)
       return
       end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine flow (aa,bb,x,z,ni,nj,nk,nl,rvis,rhof,kx,kz,
     $ nnode,nelem,iout,mat,area,al,b,c,phi,temp,told,dt,nfband,
     $ head,theta,phiold,perm1,perm2,anisop,vs,icoup,eltp,delta,
     $ omegasum,node_angle,ss,perm,vsm,rhos,itt,ntime,presl,prslold,
     $ beta,pload1,pload2,tr,trold,rhoo,beta_ul,icase,ifltag,alpha,
     $ rotflag,load,dsig_dth,fc_pp,dsig_old,icount,it,time,d_time,
     $ fc_time,fc_crit,kx_out,fc_flag)


!     this subroutine computes the stiffness matrix and forcing
!     vectors for the fluid flow equation using three and four node
!     triangular elements. three node elements use linear
!     interpolation functions summed over the element.  four node
!     elements have the fourth node positioned at any point allong
!     one side of the element and use both linear and quadratic
!     interpolation functions summed over the element.
!
!     latest revision:
!     transient working 5/28/92 james m. wieck
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 g(4,4),b(6,maxelems),c(6,maxelems),be(8),phi(maxnodes)
      real*8 p(4,4),alpha(0:maxmats),theta,x(maxnodes),z(maxnodes)
      real*8 area(maxelems),rhof(maxnodes),vs(maxnodes),dt
      real*8 anisop(0:maxmats),phiold(maxnodes),kxe,kze,kxze,rho0,c1,c2
      real*8 al(6,maxelems),dd,ddfn(8)
      real*8 aa(maxnodes,maxwdt),bb(maxnodes),a,ae
      real*8 ev(8),fv(8),b1,b2,b3,delta,omegasum,kx_temp,kz_temp
      real*8 told(maxnodes),head(maxnodes),kx(maxelems),kz(maxelems)
      real*8 perm1(0:maxmats),perm2(0:maxmats),c3,rvis(maxnodes)
      real*8 e_angle,node_angle(maxnodes),ss(maxnodes),vsm(maxnodes)
      real*8 sse,sedens,perm(maxnodes),presl(maxnodes),prslold(maxnodes)
      real*8 pload1(maxnodes),pload2(maxnodes),tr(maxelems)
      real*8 beta_ul(0:maxmats),rhos(0:maxmats),beta(0:maxmats)
      real*8 dcosd, dsind,trold(maxelems),rhoo(maxnodes),rhok,rhofe
      real*8 temp(maxnodes),load(maxnodes)
      real*8 dsig_dth(maxnodes),fc_pp(maxnodes),dsig_old(maxnodes),dummy
      integer icount

      integer nnode,nelem,eltp(maxelems),i1,i2,i3,jj,ii,mat(maxelems)
      integer ni(maxelems),nj(maxelems),nk(maxelems),nl(maxelems),qq
      integer j,k,l,node(8),nfband,ntime,it,itt,iout,icoup,nn
      integer icase(maxbnds),ifltag,rotflag,ne,me,sw

      real*8   fc_time(maxelems),time,kmin,kmax,d_time,decay
      real*8   fc_crit,kx_out(maxnodes) 
      integer  fc_flag(maxelems)

 
        rho0 = 998.2d+00
        rhok = 1100.0d+00
        grav = 9.82d+0



        do 100 m = 1, nelem

        if(eltp(m).eq.6)goto100

        ae = area(m)
        i = ni(m)
        j = nj(m)
        k = nk(m)
        l = nl(m)
        qq=mat(m)
        if(i.eq.0) go to 100
        if(qq.eq.0) then
!        write (*,*) 'ran out of material tags, program terminating'
        stop
        endif
        rvise = (rvis(i) + rvis(j) + rvis(k))/3.d+00
        rhofe = (rhof(i) + rhof(j) + rhof(k))/3.d+00
        sedens = (rhos(qq) - rhofe)/rhofe
        sse = 1.0e-6 ! (ss(i) + ss(j) + ss(k))/3.d+00
!
!     calculate uncoupled solution if icoup=0
!
        if(icoup.eq.0)then
          rvise=1.0d+00
          rhofe=rho0
        end if

        e_angle=(node_angle(i)+node_angle(j)+node_angle(k))/3.0
        if(rotflag.eq.2) e_angle=0.0
        rrho = (rhofe - rho0)/rho0
        if(qq.eq.ifltag) then
        kz_temp = (perm(i) + perm(j) + perm(k))/3.d+00
        kx_temp=kz_temp/anisop(qq)
        e_angle=delta-omegasum
        else
        kx_temp = (perm(i) + perm(j) + perm(k))/3.d+00
        kz_temp = kx_temp/anisop(qq)
        endif

         kxe  = kx_temp * (dcosd(e_angle)**2.0) +
     $         kz_temp * (dsind(e_angle)**2.0)
          kze  =kx_temp * (dsind(e_angle)**2.0)+
     $         kz_temp * (dcosd(e_angle)**2.0)
          kxze =(kx_temp-kz_temp)*dcosd(e_angle)*
     $                dsind(e_angle)


      kmax = kxe*100
      kmin = kxe

      if(fc_pp(i).ge.fc_crit.or.fc_pp(j).ge.fc_crit
     $      .or.fc_pp(k).ge.fc_crit) then
      fc_flag(m) = 1
      fc_time(m) = dble(it)
      decay = dble(it)-fc_time(m)
      kxe = (kmax-kmin)*exp(-decay/d_time)+kmin
      kze = kxe/anisop(qq)
      kxze = 0.0
      end if

      if(fc_flag(m).gt.0) then
      if(fc_pp(i).lt.fc_crit.or.fc_pp(j).lt.fc_crit
     $      .or.fc_pp(k).lt.fc_crit) then

      decay = dble(it)-fc_time(m)
      if(decay.lt.0) decay = d_time*100

      kxe = (kmax-kmin)*exp(-decay/d_time)+kmin
      kze = kxe/anisop(qq)
      kxze = 0.0
      end if
      else if (fc_flag(m).eq.0) then      
      kxe = kmin
      kze = kmin/anisop(qq)
      kxze = 0.0
      end if


      kx_out(i) = kxe
      kx_out(j) = kxe
      kx_out(k) = kxe

          a = (rhofe*rvise)/ae
          d = rho0*sse*ae/dt
          e = -rrho*rvise*rhofe
          f = -ae*rhofe/3.0d+0

!          if(m.lt.10) write (*,*) 'd=',d

          if(eltp(m).eq.1)goto300
!
!       determines which terms go with which nodes for the stiffness
!       and capacitence matrices, and which beta and gama terms are
!       used by the terms in the stiffnes matrix
!
          sw=4
          if(eltp(m).eq.2)then
            i1 = 1
            i2 = 2
            i3 = 3
          end if
          if((eltp(m).eq.4).or.(eltp(m).eq.5))then
            i1 = 2
            i2 = 3
            i3 = 1
            it=i
            i=j
            j=k
            k=it
          end if
          if(eltp(m).eq.3)then
            i1 = 3
            i2 = 2
            i3 = 1
            it=i
            i=k
            k=it
          end if

          node(1)=i
          node(2)=j
          node(3)=k
          node(4)=l
          b1=b(i1,m)
          b2=b(i2,m)
          b3=b(i3,m)
          g1=c(i1,m)
          g2=c(i2,m)
          g3=c(i3,m)
!
!     determination of values of the constants
!
         c2 = (al(i2,m)+b2*x(l)+g2*z(l))/(2.0d+00*ae)
         c3 = (al(i3,m)+b3*x(l)+g3*z(l))/(2.0d+00*ae)
         if(c2*c3.eq.0.00) go to 301
         c1 = 1.0d+00/(c3*c2)
!
!     stiffness matrix terms
!
      g(1,1)=a*(((kxe*b1**2)+(kze*g1**2))/4.0d+00+kxze*(b1*g1/2.0d+00))
      g(1,2)=a*((kxe*((b1/4.0d+00)*(b2-(c1*c2/3.0d+00)*(b2+b3)))+(kze*
     $((g1/4.0d+00)*(g2-(c1*c2/3.0d+00)*(g2+g3)))))+kxze*(((b1/4.0d+00)
     $*(g2-(c2*c1/3.0d+00)*(g3+g2)))+((g1/4.0d+00)*(b2-(c2*c1/3.0d+00)
     $*(b3+b2)))))
      g(1,3)=a*((kxe*((b1/4.0d+00)*(b3-(c1*c3/3.0d+00)*(b2+b3)))+(kze*
     $((g1/4.0d+00)*(g3-(c1*c3/3.0d+00)*(g2+g3)))))+kxze*(((b1/4.0d+00)
     $*(g3-(c3*c1/3.0d+00)*(g3+g2)))+((g1/4.0d+00)*(b3-(c3*c1/3.0d+00)
     $*(b3+b2)))))
      g(1,4)=a*((kxe*((b1/12.0d+00)*c1*(b3+b2)))+(kze*((g1/12.0d+00)*c1*
     $(g3+g2)))+kxze*(((c1*b1/12.0d+00)*(g2+g3))+(((c1*g1)/12.0d+00)*
     $(b2+b3))))
      g(2,1)=g(1,2)
      g(2,2)=a*(0.5d+00)*(((kxe*(b2**2)+kze*(g2**2))/2.0d+00)+((c1*c2/
     $3.0d+00)*(-kxe*(b2*b3+b2**2)-kze*(g2**2+g2*g3)+((c1*c2/4.0d+00)
     $ *((kxe*
     $(b3**2+b2*b3+b2**2)+kze*(g2**2+g3**2+g2*g3)))))))+kxze*(((1/(4.0d+
     $00*ae))*(b2*g2-(c1*c2/3.0d+00)*(b2*g3+2.0d+00*b2*g2+b3*g2-(c1*c2/
     $2.0d+00)*(b3*g3+b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2))))+((1.0d+00/
     $(4.0d+00*ae))*(g2*b2-(c1*c2/3.0d+00)*(g2*b3+2.0d+00*g2*b2+g3*b2-
     $(c1*c2/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2)))))
      g(2,3)=a*(0.25d+00)*(kxe*((b2*b3))+(kze*(g2*g3))-(c1*c2/3.0d+00)*
     $(kxe*((b3**2+b2*b3))+(kze*(g2*g3+g3**2)))-(c1*c3/3.0d+00)*
     $(kxe*((b2**2+b2*b3))+(kze*(g2*g3+g2**2)))+(((c1**2)*c2*c3)/6.0d+
     $00)*(kxe*((b2**2+b2*b3+b3**2))+(kze*(g3**2+g2**2+g2*g3))))+

     $kxze*(((1.0d+00/(4.0d+00*ae))*(b2*g3-(c3*c1/3.0d+00)*(b2*g3+b2*g2)
     $-(c2*c1/3.0d+00)*(b3*g3+b2*g3)+(c2*c3*c1**2/6.0d+00)*(b3*g3+b2*g2+
     $b3*g2/2.0d+00+b2*g3/2.0d+00)))+((1.0d+00/(4.0d+00*ae))*(g2*b3-(c3*
     $c1/3.0d+00)*(g2*b3+g2*b2)-(c2*c1/3.0d+00)*(g3*b3+g2*b3)+(c2*c3*c1*
     $*2/6.0d+00)*(g3*b3+g2*b2+g3*b2/2.0d+00+g2*b3/2.0d+00))))
      g(2,4)=a*((c1/12.0d+00)*(kxe*((b2*b3+b2**2))+(kze*(g2**2+g2*g3))-
     $(c1*c2/2.0d+00)*(kxe*((b2**2+b2*b3+b3**2))
     $+(kze*(g3**2+g2**2+g2*g3))
     $)))+kxze*(((c1/(12.0d+00*ae))*(b2*g3+b2*g2-(c2*c1/2.0d+00)*(b3*g3+
     $b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2)))+((c1/(12.0d+00*ae))*(g2*b3+g2
     $*b2-(c2*c1/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2))))
      g(3,1)=g(1,3)
      g(3,2)=g(2,3)
      g(3,3)=a*(0.5d+00)*(((kxe*(b3**2)+kze*(g3**2))/2.0d+00)+((c1*c3/3. 
     $d+00)*(-kxe*(b2*b3+b3**2)-kze*(g3**2+g2*g3)+((c1*c3/4.0d+00)
     $ *((kxe*
     $(b3**2+b2*b3+b2**2)+kze*(g2**2+g3**2+g2*g3)))))))+kxze*(((1.0d+00/
     $(4.0d+00*ae))*(b3*g3-(c1*c3/3.0d+00)*(b2*g3+2.0d+00*b3*g3+b3*g2-
     $(c1*c3/2.0d+00)*(b3*g3+b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2))))+((1.
     $0d+00/(4.0d+00*ae))*(g3*b3-(c1*c3/3.0d+00)*(g2*b3+2.0d+00*g3*b3+g
     $3*b2-(c1*c3/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2)))))
      g(3,4)=a*((c1/12.0d+00)*((kxe*(b2*b3+b3**2))+(kze*(g3**2+g2*g3))-
     $(c1*c3/2.0d+00)*(kxe*((b2**2+b2*b3+b3**2))
     $ +(kze*(g3**2+g2**2+g2*g3)
     $))))+kxze*(((c1/(12.0d+00*ae))*(b3*g3+b3*g2-(c3*c1/2.0d+00)*(b3*g3
     $+b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2)))+((c1/(12.0d+00*ae))*(g3*b3+
     $g3*b2-(c3*c1/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2))))
      g(4,1)=g(1,4)
      g(4,2)=g(2,4)
      g(4,3)=g(3,4)
      g(4,4)=a*(((c1**2/24.0d+00)*((kxe*(b2**2+b2*b3+b3**2))+(kze*(g2*g3
     $+g3**2+g2**2))))+kxze*(c1**2/12.0d+00)*(b3*g3+b2*g2+b3*g2/2.0d+00+
     $b2*g3/2.0d+00))
!
!     capacitence matrix terms
!
      p(1,1) = d/6.0d+00
      p(1,2) = d/12.0d+00-c1*c2*d/60.0d+00
      p(1,3) = d/12.0d+00-c1*c3*d/60.0d+00
      p(1,4) = d*c1/60.0d+00
      p(2,1) = p(1,2)
      p(2,2) = d/6.0d+00-c2*c1*d/15.0d+00+(c2**2)*(c1**2)*d/90.0d+00
      p(2,3) = d/12.0d+00-c3*c1*d/30.0d+00-c2*c1*d/30.0d+00+c2*c3*
     $(c1**2)*d/90.0d+00
      p(2,4) = d*c1/30.0d+00-c2*(c1**2)*d/90.0d+00
      p(3,1) = p(1,3)
      p(3,2) = p(2,3)
      p(3,3) = d/6.0d+00-c3*c1*d/15.0d+00+(c3**2)*(c1**2)*d/90.0d+00
      p(3,4) = d*c1/30.0d+00-c3*(c1**2)*d/90.0d+00
      p(4,1) = p(1,4)
      p(4,2) = p(2,4)
      p(4,3) = p(3,4)
      p(4,4) = d*(c1**2)/90.0d+00
!
!     detrmination of the weighted stiffness matrix
!
      do 200 n=1,4
        ddfn(n)=0.d+00
        do 250 nn=1,4
          ddfn(n) = ddfn(n) + ((1.d+00-theta)*g(n,nn) -
     $            p(n,nn))*head(node(nn))
          g(n,nn) = theta*g(n,nn)+p(n,nn)
  250 continue
!
!  for uncoupled solution; set e = 0.0
!
      if(icoup.eq.0) e = 0.0d+00

      fv(1)= (1.0d+00/3.0d+00)
      fv(2)= (1.0d+00/3.0d+00)*(1.0d+00-c2*c1/4.0d+00)
      fv(3)= (1.0d+00/3.0d+00)*(1.0d+00-c3*c1/4.0d+00)
      fv(4)= (c1/12.0d+00)

      ev(1)=c(i1,m)/2.0d+00
      ev(2)=c(i2,m)/2.0d+00-(c2*c1/6.0d+00*(c(i3,m)+c(i2,m)))
      ev(3)=c(i3,m)/2.0d+00-(c3*c1/6.0d+00*(c(i3,m)+c(i2,m)))
      ev(4)=c1/6.0d+00*(c(i3,m)+c(i2,m))
!
!  loading ala bethke and corbet
!

!
! calculate load teerm + oil generation source term
!
        if(rotflag.eq.2) f=0.0
        if(icase(1).ne.3) then
         be(n)=e*(kze*ev(n))
     $  + 3.0d+0*f*fv(n)*(-beta(qq)*phi(node(n))*
     $ (presl(node(n))-prslold(node(n)))/(rho0*grav*dt))
     $ + 3.0d+0*f*fv(n)*(phi(node(n))*alpha(qq)*(temp(node(n))
     $  -told(node(n)))/dt) - (tr(m)-trold(m))/(rhofe*dt)
     $  *3.0d+0*f*fv(n)*(rhok-rhoo(node(n))) -ddfn(n)


        pload1(node(n)) =(presl(node(n))-prslold(node(n)))

        pload2(node(n))=  f*fv(n)*(-beta(qq)*phi(node(n))
     $                    /(rho0*grav*dt))

       else
        be(n)=e*(kze*ev(n))
     $ + 3.0d+0*f*fv(n)*(-beta_ul(qq)*phi(node(n))*
     $(presl(node(n))-prslold(node(n)))/(rho0*grav*dt))
     $ + 3.0d+0*f*fv(n)*(phi(node(n))*alpha(qq)*(temp(node(n))
     $  -told(node(n)))/dt) - (tr(m)-trold(m))/(rhofe*dt)
     $  *3.0d+0*f*fv(n)*(rhok-rhoo(node(n))) -ddfn(n)


        pload1(node(n)) = (presl(node(n))-prslold(node(n)))

        pload2(node(n))= 3.0*f*fv(n)*(-beta_ul(qq)*phi(node(n))
     $   /(rho0*grav*dt))

       endif



  200 continue

  300 if(eltp(m).ne.1)goto 700
  301 continue
!
!     three node solution
!

!       write(12,494) m,ni(m),nj(m),nk(m),qq
  494  format(5(i7,2x))
      node(1)=ni(m)
      node(2)=nj(m)
      node(3)=nk(m)

      do 400 n=1,3

  556 format ('n=',i5,' rhofe=',1pe12.3,' tr=',1pe12.3,
     &   ' trold=',1pe12.3,' dt=',1pe12.3)
        dd=0.0d+00
        do 500 nn=1,3
          if(n.eq.nn)then
            p(n,nn)=d/3.0d+00
          else if (n.ne.nn)then
            p(n,nn)=0.0d+0     
          end if

         g(n,nn)=(a/4.0d+00)*(kxe*b(nn,m)*b(n,m)+kze*c(nn,m)*c(n,m)+
     $            kxze*(c(n,m)*b(nn,m)+b(n,m)*c(nn,m)))


          dd = dd + ((1.0d+00-theta)*g(n,nn)-p(n,nn))*head(node(nn))

          g(n,nn) = theta*g(n,nn)+p(n,nn)
  500   continue
!
!  for uncoupled solution; set e = 0.0
!
      if(icoup.eq.0) e = 0.0d+00


!
! calcluate ice sheet loading
!
!     f = 0.0

            dummy = dsig_dth(node(n))*dble(icount)/9.0d+0 +
     $ dsig_old(node(n))*(1-dble(icount)/9.0d+0) 
        if(m.eq.13819) write(*,*) 'sub flow, dummy = ',dummy

        be(n)=(e/2.0d+00)*(b(n,m)*kxze+c(n,m)*kze)
!     $ + f/3*(-0.91*load(node(n))*0.1*1.0e-6)-dd
     $ + f/3*(-0.91*load(node(n))*1.0e-6)-dd

!     $ + f/3*(-dummy)*1.0e-6/1.0e+4-dd
!     $ + f/3*(-dummy)*0.1*1.0e-6/1.0e+4-dd

!     $ + f/3*(-dsig_dth(node(n))*1.0e-6/1.0e+4)-dd
!     $ + f/3*(-dsig_dth(node(n))*1.0e-6*0.1/1.0e+4)-dd


  400 continue
!
!     form global a matrix and b vector
!
      sw=3

  700   do 800 ii = 1,sw
          go to (35,40,45,50),ii
   35     ne = i
          go to 55
   40     ne = j
          go to 55
   45     ne = k
          go to 55
   50     ne = l
          go to 55

   55     do 900 jj=1,sw
            go to (60,65,70,75),jj
   60       me = i
            go to 80
   65       me = j
            go to 80
   70       me = k
            go to 80
   75       me = l
   80       kl=me-ne+1
            if(kl.le.0)go to 900
            aa(ne,kl)=aa(ne,kl)+g(ii,jj)
  900     continue
          bb(ne)=bb(ne)+be(ii)
	  !write(*,*)bb(ne),be(ii)
  800   continue
c
      if(eltp(m).eq.1) then
      if(iout.ne.2) go to 401
      write(8,12) m
   12 format (' local fstiff for elem:',i4,//)
      do 402 lp=1,3
  402 write (8,21) (g(lp,mp),mp=1,3)
      write (8,13) (be(lp),lp=1,3)
   13 format (/,3x,'b vector:',3(1pe12.3))
  401 continue
      else if (eltp(m).gt.1) then
      if(iout.ne.2) go to 1401
      write(8,12) m
      do 1402 lp=1,4
 1402 write (8,21) (g(lp,mp),mp=1,4)
      write (8,13) (be(lp),lp=1,4)
 1401 continue
      endif
  100 continue

!
!     print out stiffness matrix
!
      if (iout.ne.2) go to 902
      write(8,11)
      do 901 ll=1,nnode
      write (8,21) (aa(ll,mm),mm=1,11)
  901 continue
      do ll=1,nnode
      write (8,31) ll,bb(ll)
      end do
  902 continue

   11 format (/1x,'sub. flow; stiffness matrix:')
   21 format (11(1pe12.3))
   31 format ('n=',i5,1x,'b:',1pe12.3)
      return
      end

!***********************************************************************
!
      subroutine depths (nrows,ncols,zsp,zmax,nflt,bec,sfcel,xsp,it
     $                  ,delsfc,head,vs,dt)
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
!
      real*8 zsp(maxnodes),zmax(maxnodes),sfcel(maxnodes),xsp(maxnodes)
      real*8 tempsfcel(maxnodes),head(maxnodes),vs(maxnodes),dt,slope
      real*8 delsfc(maxnodes)


      integer nrows(maxbnds),nflt,bec(maxbnds)
      integer sfcnode(maxbnds),sfccnt,ncols
      integer lnode,rnode,it,sfcflag




      nn=0
      sfccnt=0
      do 100 n=1,ncols
         do 50 m=1,nrows(n)
	 nn=nn+1
          if(m.eq.nrows(n))then
             sfccnt=sfccnt+1
	     sfcnode(sfccnt)=nn
	  endif
  50   continue
  100 continue


      nn=0
      do 600 n=1,ncols
         do 400 m=1,nrows(n)
	   nn=nn+1

	   sfcflag=0
           do 200 l=1,sfccnt

	    if(sfcflag.eq.0)then

	     lnode=sfcnode(l)
	     if(l.eq.sfccnt)then
	       rnode=lnode
	     else
	          rnode=sfcnode(l+1)
	     endif

	      if(xsp(nn).ge.xsp(lnode).and.xsp(nn).lt.xsp(rnode))then
	       slope=(zsp(lnode)-zsp(rnode))/(xsp(lnode)-xsp(rnode))
	       tempsfcel(nn)=zsp(lnode)+slope*(xsp(nn)-xsp(lnode))
	         if(it.gt.1)then
		  delsfc(nn)=sfcel(nn)-tempsfcel(nn)
		  sfcel(nn)=tempsfcel(nn)
	      else
	       delsfc(nn)=0.d0
	       sfcel(nn)=tempsfcel(nn)
	      endif

	       sfcflag=1
	      else if(l.eq.sfccnt.and.xsp(nn).ge.xsp(lnode)) then
	        tempsfcel(nn)=zsp(lnode)
		if(it.gt.1)then
		  delsfc(nn)=sfcel(nn)-tempsfcel(nn)
		  sfcel(nn)=tempsfcel(nn)
	      else
	       delsfc(nn)=0.d0
	       sfcel(nn)=tempsfcel(nn)
	      endif
		sfcflag=1
	      endif

	     endif

  200      continue
	  head(nn)=head(nn)-delsfc(nn)
          vs(nn)=vs(nn)+delsfc(nn)/dt
  400   continue
  600 continue

      iout = 0
      return
       end
!***********************************************************************
       subroutine compact(zsp,nnode,phi,vzt,ncols,nrows,dt,
     $                 matp,phiold,vs,icase,it,rotflag,vsm
     $                 ,ver)
!
!
!***********************************************************************

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 phi(maxnodes),phiold(maxnodes),phihlf,dt,vs(maxnodes)
      real*8 vzt(maxbnds),zsp(maxnodes),a1,a2,delzz,dz
      real*8 ver(maxbnds),delphi,vsm(maxnodes)

      integer matp(maxnodes),nnode,nrows(maxbnds),ncols
      integer nn,nnn,icase(maxbnds),rotflag

      nn = 1
      nnn = 0
      do 100 n=1,ncols

      nnn=nnn+nrows(n)

      if(rotflag.ne.2) vs(nn) = vzt(n)

      do 200 m=1,nrows(n)-1

!
! compute average delz
!
      dz = zsp(nn+1) - zsp(nn)
      a1 = dz*(1.d+00 - phiold(nn))/(1.d+00 - phi(nn))
      a2 = dz*(1.d+00 - phiold(nn+1))/(1.d+00 - phi(nn+1))
      delzz = (a1 + a2)/2.d+00

!
! compute average 1/(1-phi)
!

       a1 = 0.5d+00*( 1.0d+00/(1.0d+00-phiold(nn))
     $                     + 1.0d+00/(1.0d+00-phi(nn)) )
      a2 = 0.5d+00*( 1.0d+00/(1.0d+00-phiold(nn+1))
     $                     + 1.0d+00/(1.0d+00-phi(nn+1)) )
      phihlf = (a1 + a2)/2.d+00

!
! compute average phi-phi0ld
!
      a1 = phi(nn) - phiold(nn)
      a2 = phi(nn+1) - phiold(nn+1)
      delphi = (a1 + a2)/2.d+00

!
! compute velocity of solids and new z-coordinate locations
!
      if(rotflag.eq.0) then
      vs(nn+1) = vs(nn) + delzz*phihlf*(delphi)/dt
      else if (rotflag.eq.1) then
      vs(nn+1) = vzt(n)
      endif
      nn=nn+1
  200 continue
      nn=nnn+1
  100 continue
  500 continue
      nn = 0
      nnn = 0
      do 102 n=1,ncols
      nnn=nnn+nrows(n)
      do 202 m=1,nrows(n)
      nn=nn+1
      if(icase(n).eq.0.and.rotflag.ne.2) vsm(nn) = vzt(m)
      if(icase(n).eq.1.and.rotflag.ne.2) vsm(nn) = vs(nn)
      if(icase(n).eq.2.and.rotflag.ne.2) vsm(nn) = 0.0d+0
      if(icase(n).eq.3.and.rotflag.ne.2) vsm(nn) = ver(m)
      if(rotflag.eq.2) vsm(nn) = vs(nn)
 202  continue
 102  continue
 250  return
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      subroutine wrratio (qx,qz,ni,nj,nk,nl,eltp,nelem,mat,area,phi
     $ ,wrre,twrre,dt,phitot,b,c,temp)
!
!     this subroutine computes the water rock ration
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 qx(maxelems),qz(maxelems),area(maxelems),phi(maxnodes)
      real*8 phie,temp_b,temp_c,dt,wrre(maxelems),phitot(maxelems)
      real*8 temp(maxnodes),tempe,b(6,maxelems)
      real*8 c(6,maxelems),twrre(maxelems)


      integer nelem,eltp(maxelems), ni(maxelems),nj(maxelems)
      integer nk(maxelems),nl(maxelems), mat(maxelems)

      do 100 m=1,nelem

        if(eltp(m).eq.6)goto 100

        i=ni(m)
        j=nj(m)
        k=nk(m)

        phie=(phi(i) + phi(j) + phi(k))/3.d+00

        delx=0.0
        delz=0.0
        do 77 kk=1,3
          if(dabs(b(kk,m)).gt.dabs(delz)) delz=dabs(b(kk,m))
          if(dabs(c(kk,m)).gt.dabs(delx)) delx=dabs(c(kk,m))
  77    continue

        if(phie.gt.0.0) wrre(m)= wrre(m)+dabs(delx*qz(m)/phie
     $    + delz*qx(m)/phie ) * dt
  100 continue

      return
      end
!***********************************************************************
      subroutine sfcout (nnode,nelem,area,phi,ni,nj,nk,nl,qx,qz,
     $ x,z,temp,head,matp,wrre,phitot,iflow,mat,hflux_z
     $ ,tncols,tnrows,time,ntime,iprint,it,bec,nflt,wrrn
     $ ,iprfst,msl)


!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 area(maxelems),x(maxnodes),z(maxnodes),phi(maxnodes)
      real*8 vxp(maxnodes),vzp(maxnodes)
      real*8 artot(maxnodes),qx(maxelems),qz(maxelems),temp(maxnodes)
      real*8 wrre(maxelems),wrrn(maxnodes),phitot(maxelems)
      real*8 hfzn(maxnodes),time,msl,head(maxnodes), hflux_z(maxelems)

      integer ni(maxelems),nj(maxelems),nk(maxelems),nl(maxelems)
      integer nnode,nelem,ip,matp(maxnodes),iflow,mat(maxelems)
      integer tncols,tnrows(maxbnds),icnt,ntime,iprint,it,bec(maxbnds)
      integer nflt,iprfst,ihead

      data ihead /0/


      do 101, ip=1,nnode
      vxp(ip) = 0.0d+0
      vzp(ip) = 0.0d+0
      artot(ip) = 0.0d+0
      hfzn(ip)=0.0d+0
  101 continue

      do 201 m = 1,nelem
!
! determine number of elements that share a common node
!
       i = ni(m)
       j = nj(m)
       k = nk(m)
       l = nl(m)
      if(i.eq.0) go to 201

        vxp(i) = qx(m)*area(m)/phi(i) + vxp(i)
        vzp(i) = qz(m)*area(m)/phi(i) + vzp(i)
        hfzn(i)= hflux_z(m)*area(m)+hfzn(i)
        vxp(j) = qx(m)*area(m)/phi(j) + vxp(j)
        vzp(j) = qz(m)*area(m)/phi(j) + vzp(j)
        hfzn(j)= hflux_z(m)*area(m)+hfzn(j)
        vxp(k) = qx(m)*area(m)/phi(k) + vxp(k)
        vzp(k) = qz(m)*area(m)/phi(k) + vzp(k)
        hfzn(k)= hflux_z(m)*area(m)+hfzn(k)
        if(l.gt.0) vxp(l) = qx(m)*area(m)/phi(l) + vxp(l)
        if(l.gt.0) vzp(l) = qz(m)*area(m)/phi(l) + vzp(l)
        if(l.gt.0) hfzn(l)= hflux_z(m)*area(m)+hfzn(l)
        artot(i) = area(m) + artot(i)
        artot(j) = area(m) + artot(j)
        artot(k) = area(m) + artot(k)
        if(l.gt.0) artot(l) = area(m) + artot(l)

  201 continue
!
! account for differences in element area in weighted velocity
!
      do 303, ip=1,nnode
      vxp(ip) = vxp(ip)/artot(ip)
      vzp(ip) = vzp(ip)/artot(ip)
      hfzn(ip)=hfzn(ip)/artot(ip)
	  hfzn(ip)=hfzn(ip)/artot(ip)

  303 continue


c print out to scatter data
      if(ihead.eq.0)then
          write(50,*)'title = "surface data"'
       write(50,*)'variables = "x (m)","jz (mw/m*m)",
     $    "qx(m/yr)", "qz(m/yr)"'
      endif
          write(50,*)'zone f="point",i=',tncols
       icnt = 0
       do 62 jj=1,tncols
         do 42 mj=1,tnrows(jj)
	   icnt=icnt+1
 42      continue
         write(50,*)x(icnt),temp(icnt),hfzn(icnt),vxp(icnt),vzp(icnt)
 62    continue
       ihead=1

c
c  format statements
c
10    format ('title = ','"',a15,'"')
11    format ('variables = "x","z","vx","vz","temp"
     $ ,"head","strat","wrr","hfluxz"')
12    format ('zone n=',i5,', e=',i5,', f=fepoint, et=triangle ')
13    format (11(1pe12.4,1x))
15    format(4(i5,1x))
c
500   continue
      return
      end
!***********************************************************************

      subroutine tecout (nnode,nelem,area,phi,ni,nj,nk,nl
     $  ,x,z,temp,head,matp,dummy,wrre,phitot,mat,tncols
     $  ,tnrows,time,ntime,iprint,it,nflt,wrrn,msl
     $  ,pres,presl,sigmae,vsm,ss,rhof,tectlt,tr
     $  ,masoil,masgas,rv,oilvol,gasvol,persat,mash2o,masco2
     $  ,heado,headg,rhoo,rhog,conc,perm,iflow,ibrine,iheat,ioil
     $  ,vxp,vzp,xp,zp,concp,numwel,wellnm,wellid,ipflag
     $  ,pload1,pload2,prslold,voxp,vozp,presh,node_angle
     $  ,anisop,ifltag,elhom,dsig_dth,sig_1h,sig_3h
     $  ,sign_pp,tau_pp,fc_pp,nu,xls,sigxh,sigyh,fch,nun
     $  ,d_time,fc_time,fc_crit,kx_out)

!
!  this subroutine generates graphical output in tecplot format
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!


      real*8 area(maxelems),x(maxnodes),z(maxnodes),phi(maxnodes)
      real*8 vxp(maxnodes),vzp(maxnodes),xp(maxnodes),zp(maxnodes)
      real*8 artot(maxnodes),temp(maxnodes),head(maxnodes)
      real*8 wrrn(maxnodes),vsm(maxnodes),ss(maxnodes),perm(maxnodes)
      real*8 rhof(maxnodes),twrre(maxelems),wrre(maxelems)
      real*8 time,msl,phitot(maxelems),anisop(0:maxmats),presh(maxnodes)
      real*8 pres(maxnodes),presl(maxnodes),sigmae(maxnodes)
      real*8 masoil(maxelems),masgas(maxelems),tr(maxelems)
      real*8 gasvol(maxelems),persat(maxelems),mash2o(maxelems)
      real*8 masond(maxnodes),masgnd(maxnodes),masco2(maxelems)
      real*8 heado(maxnodes),headg(maxnodes),rhoo(maxnodes)
      real*8 conc(maxnodes),concp(maxnodes),rvnd(maxnodes),rv(maxelems)
      real*8 trnd(maxnodes),node_angle(maxnodes),rhog(maxnodes)
      real*8 pload1(maxnodes),pload2(maxnodes),prslold(maxnodes)
      real*8 voxp(maxnodes),vozp(maxnodes)
      real*8 kx(maxnodes),kz(maxnodes),kxz(maxnodes),kz_temp,kx_temp
      real*8 dcosd, dsind,kxxx(maxnodes),kzzz(maxnodes)
      real*8 nun(maxnodes)

! add cvfem variables

      real*8 dsig_dth(maxnodes),sig_1h(maxnodes),sig_3h(maxnodes)
      real*8 sign_pp(maxnodes),tau_pp(maxnodes),fc_pp(maxnodes)
      real*8 nu(maxbnds),xls(maxbnds),kxe
      real*8 sigxh(maxnodes),sigyh(maxnodes),fch(maxnodes)
      integer ni(maxelems),nj(maxelems),nk(maxelems),nl(maxelems)
      integer nnode,nelem,ip,matp(maxnodes),dummy,ipflag(maxnodes)
      integer nni,nnj,nnk,ntemp,kk,mat(maxelems)
      integer tncols,tnrows(maxbnds),ntime,iprint,it,ifltag
      integer nflt,iflow,iheat,ioil,ibrine,qq,wellid(maxflts),numwel
      integer elhom(maxnodes)

      character*15 title
      character*40 tectlt
      character*20 wellnm(maxflts)

      real*8   fc_time(maxelems),kmin,kmax,d_time,decay
      real*8   fc_crit,kx_out(maxnodes)


      title = 'rift2d basin  '

       ncount=0
       do 508 m=1,nelem
         twrre(m)=wrre(m)/(1-phitot(m))
        if(ni(m).eq.0) go to 508
          ncount=ncount+1
  508 continue
!
! calculate weighted nodal velocities
!
      do 101, ip=1,nnode
       artot(ip) = 0.0d+0
       wrrn(ip) = 0.0d+0
       masond(ip) = 0.0d+0
       masgnd(ip) = 0.0d+0
       rvnd(ip) = 0.0d+0
       trnd(ip) = 0.0d+0
  101 continue

      do 201 m = 1,nelem
!
! determine number of elements that share a common node
!
       i = ni(m)
       j = nj(m)
       k = nk(m)
       l = nl(m)
      if(i.eq.0) go to 201
        masond(i) = masoil(m)*area(m) + masond(i)
        masgnd(i) = masgas(m)*area(m) + masgnd(i)
	rvnd(i) = rv(m)*area(m) + rvnd(i)
	trnd(i) = tr(m)*area(m) + trnd(i)
        wrrn(i)= wrre(m)*area(m) + wrrn(i)

        masond(j) = masoil(m)*area(m) + masond(j)
        masgnd(j) = masgas(m)*area(m) + masgnd(j)
	rvnd(j) = rv(m)*area(m) + rvnd(j)
	trnd(j) = tr(m)*area(m) + trnd(j)
	wrrn(j)= wrre(m)*area(m) +wrrn(j)

        masond(k) = masoil(m)*area(m) + masond(k)
        masgnd(k) = masgas(m)*area(m) + masgnd(k)
	rvnd(k) = rv(m)*area(m) + rvnd(k)
	trnd(k) = tr(m)*area(m) + trnd(k)
	wrrn(k)= wrre(m)*area(m) +wrrn(k)

        if(l.gt.0) then
	wrrn(l)= twrre(m)*area(m) +wrrn(l)
        masond(l) = masoil(m)*area(m) + masond(l)
        masgnd(l) = masgas(m)*area(m) + masgnd(l)
	rvnd(l) = rv(m)*area(m) + rvnd(l)
	trnd(l) = tr(m)*area(m) + trnd(l)
        endif

        artot(i) = area(m) + artot(i)
        artot(j) = area(m) + artot(j)
        artot(k) = area(m) + artot(k)

        if(l.gt.0) artot(l) = area(m) + artot(l)

  201 continue
!
! account for differences in element area in weighted velocity
!
      do 303, ip=1,nnode
        masond(ip) = masond(ip)/artot(ip)
        masgnd(ip) = masgnd(ip)/artot(ip)
	rvnd(ip)   = rvnd(ip)/artot(ip)
	trnd(ip)   = trnd(ip)/artot(ip)
        wrrn(ip)   = wrrn(ip)/artot(ip)
        qq = matp(ip)


        if(qq.eq.ifltag) then
        kz_temp = perm(ip)
        kx_temp = perm(ip)/anisop(qq)
        else
        kx_temp = perm(ip)
        kz_temp = perm(ip)/anisop(qq)
        endif
        e_angle = node_angle(ip)

          kx(ip)  = kx_temp * (dcosd(e_angle)**2.0) +
     $         kz_temp * (dsind(e_angle)**2.0)
         kz(ip)  = kz_temp * (dcosd(e_angle)**2.0) +
     $         kx_temp * (dsind(e_angle)**2.0)
          kxz(ip) =(kx_temp-kz_temp)*dcosd(e_angle)*
     $                           dsind(e_angle)

!      kmax = kx(ip)*100
!      kmin = kx(ip)
!
!      if(fc_pp(i).ge.fc_crit.or.fc_pp(j).ge.fc_crit
!     $      .or.fc_pp(k).ge.fc_crit) then
!
!      decay = dble(it)-fc_time(m)
!      kxe = (kmax-kmin)*exp(-decay/d_time)+kmin
!      kze = kxe/anisop(qq)
!      kxze = 0.0
!      end if

!      if(fc_pp(i).lt.fc_crit.or.fc_pp(j).lt.fc_crit
!     $      .or.fc_pp(k).lt.fc_crit) then

!     decay = dble(it)-fc_time(m)
!     if(decay.lt.0) decay = d_time*100

!      kxe = (kmax-kmin)*exp(-decay/d_time)+kmin
!      kze = kxe/anisop(qq)
!      kxze = 0.0
!      end if

!       kxxx(ip) = kxe
       kzzz(ip) = kx_out(ip)/anisop(qq)

  303 continue
       if(dummy.eq.0) then
          write(23,*) 'variables = "x","nu"'
          write(11,10) title
          write(11,31)
	  if(iflow.eq.1) write(11,32)
	  if(iheat.eq.1) write(11,33)
	  if(ibrine.eq.1)write(11,34)
	  if(ioil.eq.1)  write(11,35)
       endif
       write(11,12) nnode,ncount
       write(23,*) 'zone'
       do ic = 1,tncols
       write(23,*) xls(ic),nu(ic)
       end do
!        1         2         3         4         5         6         7
!23456789012345678901234567890123456789012345678901234567890123456789012
       do 400 ip=1,nnode
          write (11,*)x(ip),z(ip),matp(ip),0.9*nun(ip),kx_out(ip),
     $ kzzz(ip),head(ip)
        if(iflow.eq.1) write(11,*)vxp(ip),vzp(ip),
     $  head(ip)-1.0e5-nun(ip)*0.9,
     $     dsig_dth(ip),sig_1h(ip),sig_3h(ip),
     $     sign_pp(ip),tau_pp(ip),fc_pp(ip),pload1(ip),
     $     sigxh(ip),sigyh(ip),fch(ip)
!     $       pres(ip),presh(ip),presl(ip),phi(ip),rhof(ip)
!     $       ,pload1(ip),ss(ip)
  
          if(iheat.eq.1) write(11,*)temp(ip)
          if(ibrine.eq.1)write(11,*)conc(ip),xp(ip),zp(ip),
     $            ipflag(ip),elhom(ip),concp(ip)
          if(ioil.eq.1)  write(11,*) masond(ip),masgnd(ip),
     $          trnd(ip),rvnd(ip),heado(ip),voxp(ip),vozp(ip)
 400   continue

!
! print out connectivity list after sortiung for matp
!
       do 501 m=1,nelem
          if(ni(m).eq.0) go to 501
           nni=ni(m)
           nnj=nj(m)
           nnk=nk(m)
          kk=1
         do while ((matp(nni).ne.mat(m)).and.(kk.lt.3))
           ntemp=nnk
           nnk=nnj
           nnj=nni
           nni=ntemp
           kk=kk+1
         enddo
         write(11,15) nni,nnj,nnk
  501 continue
      nn = 0
      nnn = 0
      if(dummy.eq.0)then
!
! print out profile data for each well
!
      do 1302 ll=1,numwel
              n1=ll+54
              write(n1,10) title
	      write(n1,51)
	  if(iflow.eq.1) write(n1,52)
	  if(iheat.eq.1) write(n1,53)
	  if(ibrine.eq.1)write(n1,54)
	  if(ioil.eq.1)  write(n1,55)
          write(n1,56) wellnm(ll)
1302  continue
      endif
      do 302 ll=1,numwel
       n1=ll+54
       ii = wellid(ll)
       nn=0

       write(n1,*)'zone f="point",i=',tnrows(ii)
       do 102 n=1,tncols
        nnn=nnn+tnrows(n)
        do 202 m=1,tnrows(n)
         nn=nn+1
         if(n.eq.ii) then
         depth = z(nnn) - z(nn)
         write(n1,509) depth
         if(iflow.eq.1) write(n1,*) head(nn)-1.0e5,pres(nn)
     $ ,presl(nn),presh(nn),vsm(nn),phi(nn),ss(nn),rhof(nn)
        if(iheat.eq.1) write(n1,*) temp(nn)
        if(ibrine.eq.1) write (n1,*) conc(nn)
        if(ioil.eq.1) write (n1,*) masond(nn),masgnd(nn)
     $ ,rvnd(nn),trnd(nn),rhoo(nn),rhog(nn)

       endif
 202    continue
 102   continue
 302  continue


10    format ('title = ','"',a15,'"')
12    format ('zone n=',i7,', e=',i7,', f=fepoint, et=triangle ')
13    format (11(1pe12.4,1x))
15    format(4(i7,1x))

31    format ('variables = "x","z","strat","node_angle",
     $ "kx","kz","kxz"')
32    format (',"vx","vz","head","dsig_dth","sig_1h",',
     $  '"sig_3h","sign_pp","tau_pp","fc_pp","pload1"',
     $  ',"sigx","sigy","fc"')
33    format (',"temp"')
34    format (',"conc","xp","zp","ipflag","elhom","concp"')
35    format (',"masoil","masgas","tr","rv","headoil",'
     $  ,'"vox","voz"')

41    format ('z	depth')
42    format ('z	depth	head	pres	presl	presh	',
     $ 'masoil	phi	ss	rhof	tr	temp	rv	',
     $ 'masgas')
43    format ('	temp')
44    format ('	conc')
45    format ('	masoil	masgas	rf')

51    format ('variables = "depth"')
52    format (',"head","pres","presl","presh","vsm","phi","ss","rhof"')
53    format (',"temp"')
54    format (',"conc"')
55    format (',"masoil","masgas","rf","tr","rhoo","rhog"')
56    format ('text  x=50.0,y=50.0,t=" ',a20,' "')
 509  format (20(1pe12.3,'	'))
 510    format ('z','	','depth','	','head','	',
     $  'pres','	','presl','	','presh','	',
     $  'vsm','	','phi','	','ss','	','temp',
     $  '	','masoil','	','rv')
 511    format ('variables="z","depth","head","pres"
     $  ,"presl","presh","vsm","phi","ss","temp","masoil","rv"')
 1120 format(6(1pe12.4,3x),i6)
 1511    format ('variables="x","z","conc",
     $ "xp","zp","concp","elhom"')
 1512 format ('zone i=',i7,',  f=point')
      return
      end
!***********************************************************************
!
      subroutine set_zmax (nrows,zsp,zmax,hdbase,nflt
     $                    ,bec,time,tha,hdinc,ihchk,dt)
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 zsp(maxnodes),zmax(maxnodes),hdbase(maxbnds),hdinc(maxbnds)
      real*8 tha,time

      integer nrows(maxbnds),colcount,nflt,bec(maxbnds),ihchk,add

       nn=0
       colcount=1
       add=0
       do 15 i=1,nflt+1
          add=i-1
            do 25 j=colcount,colcount+bec(i)
              do 35 k=1,nrows(j)
	        nn=nn+1
                if(ihchk.eq.2) then
                   zmax(nn)=hdbase(j-add)
     $                   + hdinc(j-add)*(time-tha)
                   if(k.eq.nrows(j))then
                   endif
                 else
                   zmax(nn)=hdbase(j-add)
		endif

 35           continue
 25         continue
          colcount=colcount+bec(i)+1
	  add=1
15    continue

      return
       end
c
c
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
       subroutine rotate(xsp,zsp,nrows,ncols,dt,icase,delta
     $	                ,omegasum,vs,bec,nflt,domega,faultx,faultz
     $                  ,zmax,centx,centz,rotflag,bwidth
     $                  ,fbwidth,alphas,delcx,node_angle)
!
! this subroutine rotates the mesh and accomadates extension
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
! 
       real*8   xsp(maxnodes),zsp(maxnodes),faultx(maxflts)
       real*8   fbwidth(maxflts),bwidth(maxflts),delcx(maxflts)
       real*8   delta,omegasum,dt,domega,faultz(maxflts),centx(maxflts)
       real*8   centz(maxflts)
       real*8   zmax(maxnodes),t_centx(20),vs(maxnodes)
       real*8   minx,alphas,xp,zp,node_angle(maxnodes)
       real*8   alpha1,alpha4,alpha6
       real*8   sidea,sideb,sidec,sided,dip
       real*8   dcosd, dsind, dtand
       real*4  theta1

       integer  colcount,rotflag,icase(maxbnds)
       integer  nn,nflt,nrows(maxbnds),bec(maxbnds),ncols,n

        minx=0.d0
!
! calculate fault block widths from original width, the
!  amount of rotation, and the original dip using trig.
!  relationships.
!
       omegasum=omegasum+domega
       dip=-1.d0*alphas
       alpha1=-1.d0*omegasum
       alpha6=dip-alpha1
       alpha4=90.d0-alpha6

! calculate block widths

       do 60 n = 1, nflt+1
          t_centx(n)=centx(n)
          sidea=bwidth(n)
          sideb=sidea*dsind(alpha1)
          sidec=sidea*dcosd(alpha1)
	  sided=sideb*dtand(alpha4)
          fbwidth(n)=sidec+sided
 60    continue

! calculate centroid locations

        centx(1)=fbwidth(1)+(10000.d0-centz(1)) * dtand(alpha4)
        do 70 n = 2, nflt+1
            centx(n) = centx(n-1)+fbwidth(n-1)+fbwidth(n)
     $              -(centz(n)-centz(n-1)) * dtand(alpha4)
 70      continue

       nn=1
       colcount=1
!
! loop over each fault block
!
       do 15 i=1,nflt+1
          do 25 j=colcount,colcount+bec(i)
             do 35 k=1,nrows(j)
!
! local cooridinates
!
               xsp(nn)=xsp(nn)-t_centx(i)
               zsp(nn)=zsp(nn)-centz(i)
               zmax(nn)=zmax(nn)-centz(i)
!
! rotation matrix
!
               xp =  xsp(nn)*dcosd(domega) + zsp(nn)*dsind(domega)
               zp = -xsp(nn)*dsind(domega) + zsp(nn)*dcosd(domega)
!
! uplift casesc  and increment node_angle
!
                if(k.eq.nrows(j))then
                 node_angle(nn)=0.0
                 if(icase(j).eq.1)then
                  if(alpha4.ne.0)then
                    xp=xp+((zp-zmax(nn))*dtand(alpha4))
                    vs(nn)=(zp-zmax(nn))/dt
                    zp=zmax(nn)
                   else
                     xp=xp
                     zp=zmax(nn)
                   endif
                 endif
		else
		   node_angle(nn)=node_angle(nn)+domega
		endif
!
! end uplift
!
               if(k.ne.nrows(j)) vs(nn)=(zp-zsp(nn))/dt
               xsp(nn)=xp
               zsp(nn)=zp

! global cooridinates

	     xsp(nn)=xsp(nn)+centx(i)
             zsp(nn)=zsp(nn)+centz(i)
             if(xsp(nn).le.minx) minx=xsp(nn)
             zmax(nn)=zmax(nn)+centz(i)
             nn=nn+1
35          continue
25        continue
          colcount=colcount+bec(i)+1
          faultx(i)=xsp(nn-1)
          faultz(i)=zsp(nn-1)
 15    continue

!
! keep in positive space
!
       do 195 i=1,nn-1
         xsp(i)=xsp(i)+dabs(minx)
 195    continue
        do 200 i=1,nflt+1
          faultx(i)=faultx(i)+dabs(minx)
 200    continue

       return
       end
!***********************************************************************
      subroutine bcread (time,vzt,vztbase,vztinc,vztim1,vztim2,vztim3
     & ,vztim4,vztim5,vztim6,icase,ver,ihchk,hdbase,hdinc,ncols,dt
     $ ,domega,omega,centz,ntime,nflt,rotflag,centx
     $ ,delta,omegasum,it,total_time,cnbase,nevap,cevap,tha,thb,vztim7
     $ ,vztim8,vztim9,vztim10,vztim11,vztim12,grdfac,mat,ber,bec
     $ ,iprint,imtag,msl)

!
!    this subroutine updates the boundary conditions for
!    each of the tectonic time periods
!***********************************************************************

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!


      real*8 time,vztbase(maxbnds),vztinc(maxbnds),ver(maxbnds)
      real*8 tha,thb,hdbase(maxbnds),hdinc(maxbnds),cnbase(maxbnds)
      real*8 dt,domega,omega,centz(maxflts),nexttime,tcentz(maxflts)
      real*8 delta,omegasum,cevap(maxbnds),grdfac(maxbnds)
      real*8 vztim1,vztim2,vztim3,vztim4,vztim5,vztim6,vzt(maxbnds)
      real*8 vztim7,vztim8,vztim9,vztim10,vztim11,vztim12
      real*8 time1,time2,time3,time4,time5,time6,centx(maxflts)
      real*8 time7,time8,time9,time10,time11,time12
      real*8 dtand,msl

      character*79 title

      integer dummy,ihchk,icase(maxbnds),ntime,nflt,rotflag
      integer it,nevap(maxbnds),mat(maxelems),ber(maxbnds),bec(maxbnds)
      integer icnt1,icnt2,iprtemp,imtag

      time1 = dabs(time-vztim1)
      time2 = dabs(time-vztim2)
      time3 = dabs(time-vztim3)
      time4 = dabs(time-vztim4)
      time5 = dabs(time-vztim5)
      time6 = dabs(time-vztim6)
      time7 = dabs(time-vztim7)
      time8 = dabs(time-vztim8)
      time9 = dabs(time-vztim9)
      time10 = dabs(time-vztim10)
      time11 = dabs(time-vztim11)
      time12 = dabs(time-vztim12)

      if ( time1.le.1.0e-03 .or. time2.le.1.0e-03 .or.
     $     time3.le.1.0e-03 .or.time4.le.1.0e-03 .or.  
     $     time5.le.1.0e-03 .or.time6.le.1.0e-03 .or.  
     $     time7.le.1.0e-03 .or.time8.le.1.0e-03 .or.  
     $     time9.le.1.0e-03 .or.time10.le.1.0e-03 .or.  
     $     time11.le.1.0e-03 .or.time12.le.1.0e-03 ) 
     $  then
!
!  read x-, z-spacing between element rows and tectonic
!  subsidence rate boundary conditions
!
      write(*,*) 'reading in x-, z-spacing between element rows'
      write(*,*) ' and tectonic subsidence rate boundary conditions'
      write(8,*) 'reading in x-, z-spacing between element rows'
      write(8,*) ' and tectonic subsidence rate boundary conditions'
      write(*,*)
      read (7,05) title
      read (7,05) title
      read (7,05) title
      write(8,800)
      do 936 n=1,ncols
      read (7,*) dummy,vztbase(n),vztinc(n),ver(n),icase(n),grdfac(n)
       write (8,803) n,vztbase(n),vztinc(n),ver(n),icase(n)
     $ ,grdfac(n)
       vzt(n) = vztbase(n)
  936 continue
!
!     read in boundary conditions for fluid flow and heat
!                           transport
!
      do 233 jj=1,480
           hdbase(jj)=0
	   hdinc(jj)=0
  233 continue

      write(*,*) 'reading in specified head boundary conditions'
      write(8,*) 'reading in specified head boundary conditions'
      write(*,*)
      read (7,05) title
      write(8,05) title
      read (7,05) title
      write(8,05) title
      read (7,*) ihchk,tha,thb,msl
       write (8,712) ihchk,tha,thb,msl,ncols
      read (7,05) title
      do 933 n=1,ncols-nflt
        read (7,*) dummy,hdbase(n),hdinc(n)
        write (8,*) dummy,n,hdbase(n),hdinc(n)
  933 continue

      write(*,*) 'reading in specified conc boundary conditions'
      write(8,*) 'reading in specified conc boundary conditions'
      write(*,*)
      write(8,912)
      read (7,05) title
      write(8,05) title
      read (7,05) title
      write(8,05) title
      do 1933 n=1,ncols-nflt
        read (7,*) dummy,cnbase(n),nevap(n),cevap(n)
        write (8,1803) n,cnbase(n),nevap(n),cevap(n)
 1933 continue

! rotation readin and setting

       read (7,05) title
       write(8,05) title
       read (7,05) title
       write(8,05) title
       read (7,*) omega,iprtemp
       if (iprtemp.gt.0) then
         iprint=iprtemp
       endif
       write(*,*)'reading in structure information'


! this section finds the next critical time by comparing the end of the simulations
! with the vztim's and then sets nextime

       if((rotflag.eq.2).or.(rotflag.eq.3))then
          if(time.eq.vztim1)then
             nexttime=vztim2
          else if(time2.le.1.0e-03)then
             nexttime=vztim3
          else if(time3.le.1.0e-03)then
             nexttime=vztim4
          else if(time4.le.1.0e-03)then
             nexttime=vztim5
          else if(time5.le.1.0e-03)then
             nexttime=vztim6
          else if(time6.le.1.0e-03)then
             nexttime=vztim7
          else if(time7.le.1.0e-03)then
             nexttime=vztim8
          else if(time8.le.1.0e-03)then
             nexttime=vztim9
          else if(time9.le.1.0e-03)then
             nexttime=vztim10
          else if(time10.le.1.0e-03)then
             nexttime=vztim11
          else if(time11.le.1.0e-03)then
             nexttime=total_time
          endif
          if(nexttime.gt.(total_time)) then
        	nexttime=total_time
          endif
!
! domega is set by taking the time to the next critical time, dividong by dt to get
! the incremental domega
!
       
          domega = -omega/((nexttime-time)/dt)
       else
          domega=0.d0
       endif
!
! this is where more info is read in for each block
!
         read (7,05) title
         do 222 i=1,nflt+1
            read(7,*)dummy,tcentz(i)
 222     continue
!
! the new centx is set by taking the present centx and centz and then
! using the angle to find the new centx depending on the spread of the centz's
!
          do 223 i=1,nflt+1
          icnt1 = dint(tcentz(i))
          icnt2 = dint(centz(i))
             if(icnt1.eq.icnt2)then
                continue
             else
             centx(i)=centx(i)+
     $                 ((centz(i)-tcentz(i))*dtand(delta+omegasum))
             endif
		centz(i)=tcentz(i)
 223     continue
!
! read material tags for new time period
!
      if(imtag.eq.3) then
      ll = 1
      oldlll = 0
      do 2251 j=1,nflt+1      
        read (7,05) title
        nr = ber(j)
        ncc = bec(j)
        mm = ll+1
        lll = (ncc*2) + oldlll
        write(8,2154) j
 2154 format (/,'fault block # ',i4,/)
        do 2215 n=1,nr
          read (7,2217) (mat(l),l=ll,lll-1,2)
          write (8,2217) (mat(l),l=ll,lll-1,2)
          do 2230 m = mm, lll, 2
              mat(m) = mat(m-1)
 2230     continue
          ll = ll + ncc*2
          lll = lll + ncc*2
          mm = mm + ncc*2
 2215   continue
        oldlll = oldlll+bec(j)*ber(j)*2
        ll=oldlll+1
 2251 continue

      else if(imtag.eq.2) then
      ll = 1
      oldlll = 0
      do 251 j=1,nflt+1      
        read (7,05) title
        nr = ber(j)
        ncc = bec(j)
        mm = ll+1
        lll = (ncc*2) + oldlll
        write(8,1154) j
 1154 format (/,'fault block # ',i4,/)
        do 15 n=1,nr
          read (7,17) (mat(l),l=ll,lll-1,2)
          write (8,17) (mat(l),l=ll,lll-1,2)
          do 30 m = mm, lll, 2
              mat(m) = mat(m-1)
   30     continue
          ll = ll + ncc*2
          lll = lll + ncc*2
          mm = mm + ncc*2
   15   continue
        oldlll = oldlll+bec(j)*ber(j)*2
        ll=oldlll+1
  251 continue
      else if (imtag.eq.1) then

      ll=1
      lll=0
      do j=1,nflt+1
        read (7,05) title
        nr=ber(j)
        ncc=bec(j)*2
        lll=ncc+lll
        write(8,1154) j
        do n=1,nr
          read (7,17) (mat(l),l=ll,lll)
          write (8,17) (mat(l),l=ll,lll)
          ll=ll+ncc
          lll=lll+ncc
        enddo
        lll=lll-ncc
      enddo     
      endif
            endif
 
!
!***********************************************************************
!                    format statements
!***********************************************************************
   05 format (a120)
   17 format (132i1)
 2217 format (132i2,/,132i2)
  712 format (//,5x,'specified head boundary condition data',
     $ //,10x,'ihchk=',i5,9x,'tha=',1pe10.3,9x,'thb=',1pe10.3,4x,
     $ 9x,'msl=',1pe10.3,4x'nhb=',i5,//)
  912 format (//,5x,'specified conc boundary condition data',
     $ //)
  800 format ('col#        vbase        vinc1         ver'
     $ ,'      icase		grdfac')
  803 format (8x,i4,3(4x,1pe10.3),4x,i4,10x,f12.2)
  815 format (8x,i4,4x,1pe10.3,3x,1pe10.3)
 1803 format (10x,'n=',i4,4x,'cnbase=',1pe10.3,4x,
     $ 'nevap=',i4,3x,'cevap=',1pe10.3)
      return
      end

!***********************************************************************
      subroutine bctime (time,head,nhb,nh,sgrad,
     $ j1,j2,njb,j1base,j2base,j1inc,j2inc,tja,temp,
     $ ntb,nt,tpbase,tpinc,tta,itchk,vzt,vztbase,vztinc,
     $ vztim1,vztim2,vztim3,vztim4,vztim5,vztim6
     & ,tncols,delz,dt,z,vs,grdfac,conc,hdbase,hdinc,nslp
     $ ,cnbase,nevap,cevap,nc,ncb,matp,tha,thb,vztim7,vztim8,vztim9
     $ ,vztim10,vztim11,vztim12,tnrows,nclay,clay,ihchk
     $ ,nflt,bec,msl,ijchk,tjb,ttb)
!
!
!  this subroutine computes the temporal changes in boundary conditions
!  for the model.
!
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 grdfac(maxbnds),conc(maxnodes),tha,thb,head(maxnodes),time
      real*8 vs(maxnodes),j1(maxbnds),j2(maxbnds),j1base(maxbnds)
      real*8 tja,j1inc(maxbnds),j2inc(maxbnds),temp(maxnodes)
      real*8 tta,vzt(maxbnds),vztbase(maxbnds),vztinc(maxbnds)
      real*8 vztim1,vztim2,vztim3,vztim4,vztim5,vztim6,dt,z(maxnodes)
      real*8 vztim7,vztim8,vztim9,vztim10,vztim11,vztim12
      real*8 delz(maxbnds),hdbase(maxbnds),hdinc(maxbnds)
      real*8 tpbase(maxbnds),tpinc(maxbnds),cevap(maxbnds),tjb,ttb
      real*8 j2base(maxbnds),msl,clay(maxbnds),cnbase(maxbnds)

      integer nh(maxbnds),nhb,ihchk,nslp(maxbnds),nc(maxbnds),ncb
      integer njb,nt(maxbnds),ntb,itchk,tncols,matp(maxnodes)
      integer nflt,bec(maxbnds),ijchk,nclay(maxbnds),tnrows(maxbnds)
      integer nevap(maxbnds)
!
!  compute temporal heat flow boundary conditions
!
      if (ijchk.eq.1) then
      do 51 n=1,njb
        j1(n) = (j1base(n) + j1inc(n)*(time-tja))*7538.57
        j2(n) = (j2base(n) + j2inc(n)*(time-tja))*7538.57
   51 continue
      endif
      if (ijchk.eq.2)then
      do 151 n=1,njb
       j1(n)=(j1base(n)+
     $  (1+dsin( ((time-tja)*2.*3.14)/tjb))*j1inc(n))*7538.57
       j2(n)=(j2base(n)+
     $  (1+dsin( ((time-tja)*2.*3.14)/tjb))*j2inc(n))*7538.57
  151 continue
      endif
      if (ijchk.eq.3)then
      do 152 n=1,njb
      j1(n)=(j1base(n)+
     $ j2inc(n)*(exp(j1inc(n)*(time-tja)/tjb)-1))*7538.57
      j2(n)=(j2base(n)+
     $ j2inc(n)*(exp(j1inc(n)*(time-tja)/tjb)-1))*7538.57
  152 continue
      endif

      if (ihchk.eq.1) then
        do 23 n=1,nhb
          head(nh(n))=z(nh(n))
   23   continue
      endif
      if (ihchk.eq.2) then
        if (time.lt.tha) go to 21
        do 24 n=1,nhb
           head(nh(n))=hdbase(n)+hdinc(n)*(time-tha)
   24   continue
   21   continue
      endif
! sea level fluctuation boundary condition
      if (ihchk.eq.3) then
        do 95 n=1,nhb
        head(nh(n))=msl+(1+dsin( ((time-tha)*2.*3.14)/thb))*hdinc(n)
        if(head(nh(n)).le.z(nh(n))) head(nh(n)) = z(nh(n))
   95   continue
      endif
   25 continue

!
!  compute temporal temperature boundary condition
!
      if (itchk.eq.0) go to 15
!
!  linear boundary condition
!
      if (itchk.eq.1) then
        if (time.lt.tta) go to 214
        do 13 n=1,ntb
         temp(nt(n)) = tpbase(n) + tpinc(n)*(time-tta)
   13   continue
  214   continue
      endif
   15 continue

! sine function boundary

      if (itchk.eq.2) then
        do 195 n=1,nhb
        temp(nt(n))=tpbase(n)
     $ +(1+dsin( ((time-tta)*2.*3.14)/ttb))*tpinc(n)
  195   continue
        endif
!
! compute specified concentration boundary conditions
!
       do 114 n=1,ncb
         if(nevap(n).eq.0)then
          conc(nc(n)) = cnbase(n)
	 else         
	 if(nevap(n).eq.matp(nc(n))) conc(nc(n)) = cevap(n)
	 if(nevap(n).ne.matp(nc(n))) conc(nc(n)) = cnbase(n)
	 endif
  114 continue
!
! specify conc. of sedimentary layer
!
      nn = 0
      do 1600 n=1,tncols
      do 1400 m=1,tnrows(n)
      nn = nn+1
      if(nclay(n).eq.matp(nn))then
      ncb=ncb+1
      conc(nn)=clay(n)
      nc(ncb) = nn
      else
      endif
 1400 continue
 1600 continue
!
!  specified seawater boundary condition if ihchk=3
!
      if (ihchk.eq.3) then
       do 1114 n=1,ncb
       if(head(nh(n)).le.z(nh(n))) then
            conc(nc(n)) =  0.00
            else
            conc(nc(n))  = 0.035
            endif

 1114   continue
       endif

!
!  compute tectonic subsidence boundary conditions
!
      if (time.gt.vztim1) then
      do 26 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim1)
   26 continue
      else if (time.gt.vztim2) then
      do 27 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim2)
   27 continue
      else if (time.gt.vztim3) then
      do 28 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim3)
   28 continue
      else if (time.gt.vztim4) then
      do 128 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim4)
  128 continue
      else if (time.gt.vztim5) then
      do 1328 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim5)
 1328 continue
      else if (time.gt.vztim6) then
      do 1228 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim6)
 1228 continue
      else if (time.gt.vztim7) then
      do 926 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim7)
  926 continue
      else if (time.gt.vztim8) then
      do 927 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim8)
  927 continue
      else if (time.gt.vztim9) then
      do 928 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim9)
  928 continue
      else if (time.gt.vztim10) then
      do 1128 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim10)
 1128 continue
      else if (time.gt.vztim11) then
      do 828 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim11)
  828 continue
      else if (time.gt.vztim12) then
      do 728 n=1,tncols
      vzt(n) = vztbase(n) + vztinc(n)*(time-vztim12)
  728 continue
      end if

   96 format ('ihchk=',i5,'   n=',i5,'   nh=',i5,'   headbc=',
     $   3x,1pe12.3,'     tha=',1pe12.3,3x,'    thb=',1pe12.3,
     $   3x,1pe12.3,'     hdbase=',1pe12.3,'    hdinc=',1pe12.3)
  196 format ('   n=',i5,'   nc=',i5,'   head=',
     $   3x,1pe12.3,'     z=',1pe12.3,3x,'    conc=',1pe12.3)

 3467 format ('col=',i5,2x,'node=',i5,2x,'conc=',1pe12.3)
c         write(8,3467) n,nc(n),conc(nc(n))
 3468 format ('solute b.c. data; time = ',1pe12.3)	 

      return
      end
!
!***********************************************************************
       subroutine readin (xsp,zsp,rhos,ntime,iout,nh,dt,mat,head,tkf
     $   ,tks,cvf,cvs,ldis,tdis,phi,theta,dif,nt,iprint,temp,conc,iheat
     $   ,ibrine,icoup,ncols,nrows,nr,delz,zmin,j1base
     $   ,j2base,j1inc,j2inc,tja,tpbase,tpinc,tta
     $   ,itchk,vzt,vztbase,vztinc,vztim1,vztim2,vztim3
     $   ,vztim4,vztim5,vztim6,phiold,perm1,perm2
     $   ,anisop,toc,grdfac,isourc,sgrad,perct1,perct2
     $   ,perct3,icase,ver,iskip,tha,thb,hdbase,hdinc,ihchk
     $   ,nflt,bec,ber,iflow,nftype
     $   ,nrel,nnr,delta,domega,bwidth,centx,centz,alphas
     $   ,omega,faultx,faultz,zmax,rotflag
     $   ,node_angle,ptec,pexp,psfc,pelem,iprfst,msl
     $   ,phi_o,beta,phi_ir,beta_ul,tectlt,ioil,nc,ncb,cnbase
     $   ,numwel,wellnm,wellid,nevap,cevap,vztim7,vztim8,vztim9
     $   ,vztim10,vztim11,vztim12,gamma,tbrstr,delzmn,maxit
     $   ,nclay,clay,iperm,imtag,dgrn,ifltag,alpha,ijchk,tjb,ttb)

!
! this subourtine reads in input data for simulation at start of simulation
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt
!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
c materials
      real*8 tkf(0:maxmats),tks(0:maxmats),toc(0:maxmats)
      real*8 cvf(0:maxmats),cvs(maxnodes),ldis(0:maxmats)
      real*8 perm1(0:maxmats),perm2(0:maxmats),anisop(0:maxmats)
      real*8 perct1(0:maxmats),perct2(0:maxmats),perct3(0:maxmats)
      real*8 bwidth(maxflts),centx(maxflts),centz(maxflts)
      real*8 phi_o(0:maxmats),beta(0:maxmats),beta_ul(0:maxmats)
      real*8 tdis(0:maxmats),dif(0:maxmats),rhos(0:maxmats)
      real*8 phi_ir(0:maxmats),faultx(maxflts),faultz(maxflts)
      integer nftype(maxflts),isourc(0:maxmats),nrel(maxflts)
c boundaries
      real*8 hdbase(maxbnds),hdinc(maxbnds),cnbase(maxbnds)
      real*8 j1base(maxbnds),j2base(maxbnds),j1inc(maxbnds)
      real*8 vzt(maxbnds),vztbase(maxbnds),tpinc(maxbnds)
      real*8 j2inc(maxbnds),tpbase(maxbnds),cevap(maxbnds)
      real*8 clay(maxbnds)

      integer nt(maxbnds),nrows(maxbnds),bec(maxbnds),ber(maxbnds)
      integer nevap(maxbnds),nh(maxbnds)
c variables
      real*8 tha,thb,sgrad,delz(maxbnds),tja,tta,tbrstr,delzmn
      real*8 dt,theta,alpha(0:maxmats),dummy
      real*8 dz,mpz,tz,alphas,omega,maxtime,delta,msl
      real*8 vztim1,vztim2,vztim3,vztim4,vztim5,vztim6,gamma
      real*8 vztim7,vztim8,vztim9,vztim10,vztim11,vztim12

      integer iprint,ntime
      integer icoup,ncols,nc(maxbnds),ncb,nclay(maxbnds)
      integer ihchk,itchk
      integer j,maxit,ifltag
      integer nflt,iflow,icount,tr1,rcount
      integer ic,nnr,iout
      integer bcounter,rotflag
      integer ptec,pexp,psfc,pelem,iprfst
      integer init_rows,ioil,numwel,iperm,imtag

! by node
      real*8 xsp(maxnodes),zsp(maxnodes),zmin(maxbnds),head(maxnodes)
      real*8 phi(maxnodes),vztinc(maxbnds),temp(maxnodes),conc(maxnodes)
      real*8 phiold(maxnodes),ver(maxbnds),xx(maxnodes),zmax(maxnodes)
      real*8 dgrn(0:maxmats),grdfac(maxbnds),node_angle(maxnodes)
      real*8 dcosd,  dtand

      integer icase(maxbnds)
! by element
      integer mat(maxelems)
! other
      character*120 title
      character*40 tectlt
      character*20 wellnm(20)

      integer wellid(maxflts)

      write (*,*)
      write (*,*)
      write (*,*) '             ***********************************'
      write (*,*) '             * reading in data and setting grid*'
      write (*,*) '             ***********************************'
      write (*,*)

      write (*,*) 'reading in simulation control parameters...'
      write (*,*)
      read (7,05) title
      tectlt=title(1:40)
      write (8,*) title
      read (7,05) title

      read (7,05)  title
      read  (7,*)  ncols,iout,iprint,iskip,iprfst,maxit
      write (8,20) ncols,iout,iprint,iskip,iprfst,maxit

      read (7,05)   title
      read (7,*)    iheat,ibrine,icoup,iflow,ioil,iperm,imtag,ifltag
      write (8,308) iheat,ibrine,icoup,iflow,ioil,iperm,imtag,ifltag

      read (7,05)  title
      read (7,*)   ntime,dt,theta,sgrad,nflt,gamma,tbrstr,delzmn
      write (8,08) ntime,dt,theta,sgrad,nflt,gamma,tbrstr,delzmn

      read (7,05) title
      write(8,12)
      do 100 l=1,nflt+1
         read(7,*)bec(l),ber(l)
	 write(8,11)l,bec(l),ber(l)
  100 continue
      read (7,05) title
      read (7,05) title
      write(8,113)
      read (7,*)ptec,pexp,psfc,pelem,nnr,numwel,msl
      write(8,114)ptec,pexp,psfc,pelem,nnr,numwel,msl

      if(nnr.gt.0) then
      read (7,05) title
      do 5550 l=1,nnr
 5550 read (7,*) dummy,nrel(l)
      write(8,115)
      write(8,116)(nrel(l),l=1,nnr)
      endif

      read (7,05) title
      read (7,122)(wellnm(l),wellid(l),l=1,numwel)
      write(8,120)
      write(8,119)(wellnm(l),wellid(l),l=1,numwel)

      write(*,*) 'reading in material property data'
      write(*,*)
      read (7,05)   title
      read (7,05)   title
      read (7,*)    nmat
      write (8,57)  nmat

      do 250 l=1,nmat

          read (7,05) title
          read (7,05) title
          read(7,*)   perm1(l),perm2(l),anisop(l),dgrn(l)
          write(8,53) l,perm1(l),perm2(l),anisop(l),dgrn(l)

          read(7,05)  title
          read(7,*)   tkf(l),tks(l),alpha(l)
          write(8,91) tkf(l),tks(l),alpha(l)
          units =  3.1558e+07/4.187
          tkf(l) = tkf(l)*units
          tks(l) = tks(l)*units

          read(7,05)  title
          read(7,*)   cvf(l),cvs(l),ldis(l),tdis(l)
          write(8,93) cvf(l),cvs(l),ldis(l),tdis(l)

          read(7,05)  title
          read(7,*)phi_o(l),beta(l),phi_ir(l),beta_ul(l),dif(l),rhos(l)
          write(8,94) phi_o(l),beta(l),phi_ir(l),beta_ul(l),
     $          dif(l),rhos(l)

          beta(l) = beta(l)*(rhos(l)-1000.0)*9.81
          beta_ul(l) = beta(l)*(rhos(l)-1000.0)*9.81

        read(7,05) title
        read(7,*) perct1(l),perct2(l),perct3(l),toc(l)
        if ((perct1(l)+perct2(l)+perct3(l)).gt.0)  then
        isourc(l) = 1
        write(8,41)
        else
        isourc(l) = 0
        endif
  250	continue
!
!  read in element material tag data (by fault block)
!

      ic=1
      do 200 k=1,nflt+1
        do 210 j=1,bec(k)
          do 220 l=1,ber(k)
            mat(ic)=1
            mat(ic+1)=1
            ic=ic+2
  220     continue
  210   continue
  200 continue
      ic=0
!
! this loop assigns material data for elements by assigning two elements
! for every property read in.
!
      write (*,*) 'reading in element material tag data..'
      write (8,*) 'reading in element material tag data..'
      write (8,*)

      if(imtag.eq.3) then
      ll = 1
      oldlll = 0
      do 2251 j=1,nflt+1      
        read (7,05) title
        nr = ber(j)
        ncc = bec(j)
        mm = ll+1
        lll = (ncc*2) + oldlll
        write(8,2154) j
 2154 format (/,'fault block # ',i4,/)
        do 2215 n=1,nr
          read (7,2217) (mat(l),l=ll,lll-1,2)
          write (8,2217) (mat(l),l=ll,lll-1,2)
          do 2230 m = mm, lll, 2
              mat(m) = mat(m-1)
 2230     continue
          ll = ll + ncc*2
          lll = lll + ncc*2
          mm = mm + ncc*2
 2215   continue
        oldlll = oldlll+bec(j)*ber(j)*2
        ll=oldlll+1
 2251 continue

      else if(imtag.eq.2) then
      ll = 1
      oldlll = 0
      do 251 j=1,nflt+1      
        read (7,05) title
        nr = ber(j)
        ncc = bec(j)
        mm = ll+1
        lll = (ncc*2) + oldlll
        write(8,1154) j
 1154 format (/,'fault block # ',i4,/)
        do 15 n=1,nr
          read (7,17) (mat(l),l=ll,lll-1,2)
          write (8,17) (mat(l),l=ll,lll-1,2)
          do 30 m = mm, lll, 2
              mat(m) = mat(m-1)
   30     continue
          ll = ll + ncc*2
          lll = lll + ncc*2
          mm = mm + ncc*2
   15   continue
        oldlll = oldlll+bec(j)*ber(j)*2
        ll=oldlll+1
  251 continue

       else if (imtag.eq.1) then

      ll=1
      lll=0
      do j=1,nflt+1
        read (7,05) title
        nr=ber(j)
        ncc=bec(j)*2
        lll=ncc+lll
        write(8,1154) j
        do n=1,nr
          read (7,17) (mat(l),l=ll,lll)
          write (8,17) (mat(l),l=ll,lll)
          ll=ll+ncc
          lll=lll+ncc
        enddo
        lll=lll-ncc
      enddo    
      endif
!
! read in rotation information
!
      read (7,05) title
      write (8,05) title
      read (7,05) title
      write (8,05) title
      read(7,*)rotflag,omega,alphas,delz(1),init_rows
      write(8,04)rotflag,omega,alphas,delz(1),init_rows
      delta = 90.d0-alphas
!
! read in centroid depth information
!
        read (7,05) title
        write (8,05) title
        do 226 i=1,nflt+1
            read(7,*)dummy,centz(i)
	    write (8,*)dummy,centz(i)
 226    continue
!
!  read x-, z-spacing between element rows and tectonic
!  subsidence rate boundary conditions
!
      write(*,*) 'reading in x-, z-spacing between element rows'
      write(*,*) ' and tectonic subsidence rate boundary conditions...'
      write(*,*)
      read (7,05) title
      write (8,05) title
      read (7,05) title
      read (7,*)  vztim1,vztim2,vztim3,vztim4,vztim5,vztim6
      write (8,123) vztim1,vztim2,vztim3,vztim4,vztim5,vztim6
      read (7,05) title
      read (7,*)  vztim7,vztim8,vztim9,vztim10,vztim11,vztim12
      write (8,124)  vztim7,vztim8,vztim9,vztim10,vztim11,vztim12
      write(8,800)
      ic=0
        read (7,05) title
      do 900 k=1,nflt+1
        read (7,05) title
        read (7,05) title
        do 936 n=1,bec(k)+1
          read (7,*) dummy,xx(n+ic),vztbase(n+ic),
     $    vztinc(n+ic),ver(n+ic),zmax(n+ic),zmin(n+ic),icase(n+ic),
     $    grdfac(n+ic)
          write (8,803) n+ic,xx(n+ic),vztbase(n+ic),
     $    vztinc(n+ic),ver(n+ic),zmax(n+ic),zmin(n+ic),icase(n+ic),
     $    grdfac(n+ic)
          vzt(n+ic) = vztbase(n+ic)
  936 continue
        ic=ic+n-1
  900 continue
!
! inititalize rotation conditions
! set incremental rotation rate
!
        if((dt*ntime).le.vztim1) then
        	maxtime=dt*ntime
        else
                maxtime=vztim1
        endif

       domega = -omega/(maxtime/dt)
!
! detemrine the block's half width
!
       bcounter = 1
       do 117 i=1, nflt+1
            bwidth(i)=xx(bec(i)+bcounter)-xx(bcounter)

              faultx(i)=xx(bec(i)+bcounter)
              faultz(i)=zmax(bec(i)+bcounter)

            bcounter = bcounter + bec(i)+1
            bwidth(i)=bwidth(i)/2.d0
117    continue

!
! for the first centroid find its position down from the surface of the block
! for each subsequent centroid find its position relative to the first
!
       do 112 i=1, nflt+1
             if(i.eq.1)then
                  centx(1) =bwidth(1)+
     $                  (10000.d0-centz(i))/dtand(-alphas)
              else
                  centx(i)=centx(i-1)+bwidth(i)+bwidth(i-1)
     $                   -(centz(i)-centz(i-1))/dtand(-alphas)
               endif

 112    continue

        if(rotflag.ne.0)then
          do 678, ke=1,ncols
	     delz(ke)=delz(1)
 678      continue
        endif

!
! generate initial grid coordinates and porosity estimate
!

       nn=0
       nnn=0
       do 901 n=1,ncols
          if(rotflag.eq.0)
     $         delz(n) = (zmax(n)-zmin(n))/init_rows
	  mpz=0.d0
	  tz=0.d0
	  dz=delz(n)
          tz=(zmax(n)-zmin(n))
          nn=nn+1
          zsp(nn)=zmin(n)
              npor=1
	  phi(nn) = phi_o(npor)*dexp(-beta(npor)*(tz))
	  phiold(nn)=phi(nn)

           if(rotflag.eq.2)then
              xsp(nn)=xx(n) + tz*dtand(-delta)
           else
              xsp(nn)=xx(n)
           endif
902       mpz=tz-delz(n)
          dz=delz(n)
          if(rotflag.eq.2) then
            dz=-1*delz(n)*dcosd(-delta)
          endif
         nn=nn+1
         tz=tz-dz
         zsp(nn)=zsp(nn-1)+dz
	 phi(nn) = phi_o(1)*dexp(-beta(1)*(zmax(n)-zsp(nn)))
	 phiold(nn)=phi(nn)
	 if(rotflag.eq.2)then
            xsp(nn) = xx(n) + tz*dtand(-delta)
         else
            xsp(nn) = xx(n)
         endif
         nnn=nnn+1
         idiff = dint(tz-delz(n))
         if(tz.gt.delz(n).and.idiff.gt.0)go to 902
         nn=nn+1
         zsp(nn)=zmax(n)
	 phi(nn) = phi_o(1)
	 phiold(nn)=phi(nn)
         xsp(nn)=xx(n)
         nrows(n)=nnn+2
         nnn=0
 901   continue
!
!     read in boundary conditions for fluid flow and heat
!                           transport

      write(*,*) 'reading in specified head boundary conditions'
      write(*,*)
      read (7,05) title
      read (7,05) title
      read (7,*) ihchk,tha,thb
      write (8,712) ihchk,tha,thb,ncols-nflt
      nn = 0
      read (7,05) title

      icount=0
      rcount=0
      tr1=0
      do 930 k=1,nflt+1
        do 933 n=1,bec(k)-tr1
          icount=icount+1
          rcount=rcount+1
          nn = nn + nrows(rcount)
          nh(icount) = nn
          nc(icount) = nn
          nt(icount) = nn
          read (7,*) dummy,hdbase(icount),hdinc(icount)
          write (8,815) icount,hdbase(icount),hdinc(icount)
          head(nh(icount)) = hdbase(icount)
  933   continue
        if(k.ne.nflt+1)then
          icount=icount+1
          rcount=rcount+1
          nn=nn+nrows(rcount)+nrows(rcount+1)
          nh(icount) = nn
          nc(icount) = nn
          nt(icount) = nn
          read (7,*) dummy,hdbase(icount),hdinc(icount)
          write (8,815) icount,hdbase(icount),hdinc(icount)
          head(nh(icount)) = hdbase(icount)
          tr1=1
          rcount=rcount+1
        else
          icount=icount+1
          rcount=rcount+1
          nn=nn+nrows(rcount)
          nh(icount) = nn
          nc(icount) = nn
          nt(icount) = nn
          read (7,*) dummy,hdbase(icount),hdinc(icount)
          write (8,815) icount,hdbase(icount),hdinc(icount)
          head(nh(icount)) = hdbase(icount)
        end if
  930 continue

      write (*,*) 'reading in boundary conditions for heat flux'
      write (*,*)
      read (7,05) title
      read (7,05) title
      read (7,*) ijchk,tja,tjb
      write (8,808) ijchk,tja,tjb,ncols
      read (7,05) title
      write(8,805)
      do 931 n=1,ncols
        read (7,*) dummy,j1base(n),j2base(n),j1inc(n),j2inc(n)
        write (8,809) n,j1base(n),j2base(n),j1inc(n),j2inc(n)
  931 continue

      write(*,*)'reading in specified temperature boundary conditions'
      write(*,*)
      read (7,05) title
      read (7,05) title
      read (7,*) itchk,tta,ttb
      write (8,714) itchk,tta,ttb,ncols-nflt
      read (7, 05) title
      write(8,819)
      do 934 n=1,ncols-nflt
        read (7,*) dummy,tpbase(n),tpinc(n)
        write (8,815) n,tpbase(n), tpinc(n)
        temp(nt(n)) = tpbase(n)
  934 continue

      write(*,*)'reading in specified concentration boundary conditions'
      write(*,*)
      read (7, 05) title
      read (7, 05) title
      write (8,1714) ncols-nflt
      write(8,1819)
      do 1934 n=1,ncols-nflt
        read (7,*) dummy,cnbase(n),nevap(n),cevap(n),nclay(n),clay(n)
        write (8,1815) n,cnbase(n),nevap(n),cevap(n),nclay(n),clay(n)
        conc(nc(n)) = cnbase(n)
 1934 continue
!
! assign initial heads, temperatures, and concentrations and porosity assuming
!  all consolidat like hydrostatic unit 1
!
      tr1=0
      icount=0
      rcount=0
      nn = 0
      nnn = 0
      do 13 k=1,nflt+1
        do 14 n=1,bec(k)-tr1
          icount=icount+1
          rcount=rcount+1
          nnn = nnn + nrows(rcount)
          do 16 m=1,nrows(rcount)
            nn = nn + 1
            head(nn) = head(nh(icount))
            temp(nn) = temp(nt(icount)) + (zsp(nnn)-zsp(nn))
     $ *(j1base(icount)*7538.57d+00)/
     $  (tks(1)*(1.d+0-phi(nn))+tkf(1)*phi(nn))
       conc(nn) = conc(nc(icount)) + (zsp(nnn)-zsp(nn))*sgrad
   16     continue
   14   continue
        if(k.ne.nflt+1)then
          icount=icount+1
          rcount=rcount+1
          nnn=nnn+nrows(rcount)+nrows(rcount+1)
          do 18 l=1,nrows(rcount)+nrows(rcount+1)
            nn=nn+1
           head(nn) = head(nh(icount))
           temp(nn) = temp(nt(icount)) + (zsp(nnn)-zsp(nn))
     $ *(j1base(icount)*7538.57d+00)/
     $  (tks(1)*(1.d+0-phi(nn))+tkf(1)*phi(nn))
         conc(nn) = conc(nc(icount)) + (zsp(nnn)-zsp(nn))*sgrad
   18     continue
          tr1=1
          rcount=rcount+1
        else
          icount=icount+1
          rcount=rcount+1
          nnn=nnn+nrows(rcount)
          do 19 l=1,nrows(rcount)
            nn=nn+1
            head(nn) = head(nh(icount))
           temp(nn) = temp(nt(icount)) + (zsp(nnn)-zsp(nn))
     $ *(j1base(icount)*7538.57d+00)/
     $  (tks(1)*(1.d+0-phi(nn))+tkf(1)*phi(nn))

      conc(nn) = conc(nc(icount)) + (zsp(nnn)-zsp(nn))*sgrad
   19     continue
        end if
   13 continue
      tr1=0
!
!***********************************************************************
!                    format statements
!***********************************************************************
   03 format (a120)
   04 format ('rotflag=',i5,2x,'omega=',f8.5,2x,'alphas='
     $ 1pe12.4,2x,'init. delz=',1pe12.4,2x,'init. rows=',i5)
   05 format (a120)
   06 format (//,30x,'***********************************',
     $         /,30x,'      p r o g r a m   r i f t      ',
     $         /,30x,'***********************************',
     $  /////,1x,a20)
   07 format (11x,i5,8x,e10.2,8x,f4.2,8x,f4.2)
   08 format  (10x,'ntime=',i5,8x,'delt=',1pe10.2,4x,'theta=',
     $ 1pe10.2,/,10x,'sgrad=',1pe10.2,5x,'nflt=',1i5,3x,'gamma=',
     $ 1pe10.2,/,10x,'tbrstr=',1pe10.2,10x,'delzmn=',1pe10.2)
   09 format (i4,5x,i10)
   10 format (3(10x,i4))
   11 format (10x,i5,6x,2i5)
   12 format (10x,'fault block   bec     ber')
   17 format (132i1)
 2217 format (132i2,/,132i2)
   20 format (///,20x,'****** rift_2d ******'///,
     $   20x,'a finite element model which simulates fluid ',
     $ /,20x,'flow, brine migration, heat transport and ',
     $ /,20x,'petroleum generation within evolving rift basins',
     $ ////,3x,' simulation control data:',//,10x,'ncols=',i4,9x,
     $ 'iout=',i4,9x,'iprint=',i4,/,10x,'iskip=',i4,'iprfst=',
     $ i4,8x,'maxit=',i4)
   41 format (/,'**this unit is a source bed**')
   53 format(/,4x,'hydrostratigraphic unit #',i2,//,10x,'perm1:'
     $ ,1pe10.2,3x,'perm2:',1pe10.2,3x,'anisop:',1pe10.2,
     $ 3x,'charac. grain size:',1pe10.2)
   57 format (/,1h1,//,1x,'number of material properties:',i2,
     $  3x,'fault tag1:',i2,3x,'fault tag2:',i2)
   91 format(10x,'tkf:',1pe10.2,5x,'tks:',1pe10.2,5x,'1:',
     $ 1pe10.2)
   93 format(10x,'cvf:',1pe10.2,5x,'cvs:',1pe10.2,5x,'ldis:',
     $ 1pe10.2,3x,'tdis:',1pe10.2)
   94 format(10x,'phi_o:',1pe10.2,4x,'beta:',1pe10.2,4x,
     $ 'phi_ir:',1pe10.2,/,'beta_ul:',1pe10.2,4x,'dif:',
     $ 1pe10.2,4x,'rhos:',1pe10.2)
  113 format (3x,'output control data:',/,
     $ 'ptec   pexp   psfc   pelem   nnr   numwel  msl')
  114 format (6(i5,5x),2x,1pe12.3)
  115 format (/, 'elements that will be outputed with ',
     $ 'time-dependent information:')
  116 format  (1('el=',i5,2x))
  120 format  ('columns which will be outputed')
  118 format  ('well name'12x,'column number')
  119 format  (1('wellnm=',a20,2x,'column #=',i5,2x))
  122 format  (a20,2x,i5)
  123 format ('vztim2=',1pe12.3,/,'vztim3=',1pe12.3,/,'vztim4='
     $ ,1pe12.3,/,'vztim5=',1pe12.3,/,'vztim6=',1pe12.3,/,
     $ 'vztim7=',1pe12.3)
  124 format ('vztim8=',1pe12.3,/,'vztim9=',1pe12.3,/,'vztim10='
     $ ,1pe12.3,/,'vztim11=',1pe12.3,/,'vztim12=',1pe12.3,/,
     $ 'vztim13=',1pe12.3)
  183 format (132i2,/,4i2)
  307 format (3(10x,i4))
  308 format (10x,'iheat = ',i3,8x,'ibrine=',i4,8x,
     $ 'icoup=',i4,8x,'iflow=',i4,8x,'ioil=',i4,8x,'iperm=',
     $ i4,8x,'imtag=',i4,8x,'ifltag=',i4)
  712 format (//,5x,'specified head boundary condition data',
     $ //,10x,'ihchk=',i5,3x,'tha=',1pe12.3,3x,'thb=',1pe12.3,3x,
     $ 'nhb=',i5,//)
  714 format (//,5x,'specified temperature boundary condition data',
     $ //,10x,'itchk=',i5,9x,'tta=',1pe10.3,5x,'ttb=',1pe10.3,
     $ 5x, 'ntb=',i5,//)
  800 format ('        col#       x            vztbase'
     $ ,'    vztinc             ver         zmax         zmin    icase'
     $ ,'      grdfac')
  801 format (1h1,//,5x,'tectonic subsidence data',//
     $ ,10x,'vztim1=',1pe10.3,4x,'vztim2=',1pe10.3,
     $ 3x,'vztim3=',1pe10.3)
  803 format (8x,i4,4x,1pe10.3,5(4x,1pe10.3),4x,i4,4x,1pe10.2)
  805 format(10x,'col',5x,'j1base',7x,'j2base',7x,'j1inc',7x,
     $ 'j2inc')
  808 format (1h1,//,5x,'specified heat flux boundary condition data',
     $ //,8x,'ijchk=',i2,8x,'tja=',1pe10.3,8x,
     $'tjb=',1pe10.3,5x,'njb=',i3,//)
  809 format (8x,i4,4x,1pe10.3,3x,1pe10.3,3x,1pe10.3,3x,1pe10.3)
  813 format (8x,i4,4x,1pe10.3,3x,1pe10.3,3x,1pe10.3,3x,1pe10.3)
  815 format (8x,i4,4x,1pe10.3,3x,1pe10.3)
  819 format (10x,'col',5x,'tpbase',7x,'tpinc')
  820 format (5x,i2,5x,i2)
 1121 format (2(1pe12.4,1x),i5)
 1714 format (//,5x,'specified concentration boundary condition data',
     $ //,5x,'ncb=',i5,//)
 1815 format (8x,i4,4x,1pe10.3,4x,i4,4x,1pe10.3,4x,i4,4x,1pe10.3)
 1819 format (10x,'col',5x,'cnbase')
      return
      end
!***********************************************************************
!
      subroutine bound(tnrows,tncols,nbnn,nbn,nh,nt,nhb,ntb
     $      ,vztbase,njn,nj1,nj2,njb,nnode,nrmax,nhband,nflt
     $      ,xfault,lbnd,rbnd,x,z,zmax,nc,ncb)
!
!  this subroutine determines boundary nodes for boundary condition
!  suboutines using the nbn array
!
!***********************************************************************

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      integer njb,ntb,nhb,nnode,ncb
      integer nbnn,nflt,count
      integer nflag(maxnodes),sflag

      integer nh(maxbnds),nt(maxbnds),nbn(maxbnds),nc(maxbnds)
      integer nj1(maxbnds),njn(maxbnds),nj2(maxbnds)
      integer tnrows(maxbnds),lbnd(maxbnds),rbnd(maxbnds),tncols

      real*8 vztbase(maxbnds)

      real*8 xfault(maxflts),x(maxnodes),z(maxnodes),zmax(maxnodes)
!
!  search for maximum number of nodes in a column
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     zeroing out boundary condition arrays
!
      do 100 n=1,tncols+nflt
        nt(n)=0
        nh(n)=0
        nc(n)=0
  100 continue

      nrmax = tnrows(1)
      do 110 n=2,tncols
      if(nrmax.lt.tnrows(n)) nrmax = tnrows(n)
  110 continue

      nbnn = tncols*2 + tnrows(1) + tnrows(tncols) - 3 + nflt
      nhb = tncols
      ntb = tncols
      ncb = tncols
!
!  make sure we are not past array dimensions
!
      if(nhband.gt.maxbnds) write (8,503) nhband
      if(nhband.gt.maxbnds) then
         write(*,*)'stop2'
         stop
       endif
!
!  impose new boundary conditions
!  these do the first boundary condition for each boundary
!
      count=2
      nt(1) = tnrows(1)
      nc(1) = tnrows(1)
      nh(1) = tnrows(1)
      nt(tncols) = nnode
      nh(tncols) = nnode
      nc(tncols) = nnode
      nn = tnrows(1)
!
! for each column this assigns the constant value to the top node
!
        do 50 n=2,tncols-1
        nn = nn + tnrows(n)
           nt(n) = nn
           nc(n) = nn
           nh(n) = nn
 50     continue

      count=0
         do 101,lf=1,nnode
	  nflag(lf)=0
 101    continue

        njb=0

        ll=0
        do 105 ncc=1,tncols
           do 115 nr=1,tnrows(ncc)
	   sflag=0
	   ll=ll+1
	   if(nflag(ll).eq.1)goto 115
!
! make sure the fault is not the next column
!
         do 108, kf=1,nflt
	    if(ll+tnrows(ncc).eq.rbnd(kf).and.
     $          rbnd(kf).lt.lbnd(kf))then
                sflag=1
            endif
!
! if two fault blocks equal flag it
!
	    if(lbnd(kf).eq.ll.and.ll.eq.rbnd(kf))then
	        sflag=2
            endif
 108     continue

!
! node 1
!
	  if(ll.eq.1)then
	     njb=njb+1
             njn(njb)=ll
             nj1(njb)=ll+1
             nj2(njb)=ll+tnrows(1)
	  endif
!
! lower right hand node
!
	  if(nr.eq.1.and.ncc.eq.tncols)then
	    njb=njb+1
            njn(njb)=ll
            nj1(njb)=ll-tnrows(ncc-1)
            nj2(njb)=ll+1
	  endif
!
! faults
!
	     do 125 nf=1,nflt
	       !write(*,*)rbnd(nf),lbnd(nf)
	       if(ll.eq.rbnd(nf).or.ll.eq.lbnd(nf))then
		    !nflag(ll)=1
!
! rbnd < lbnd
!
	            if(rbnd(nf).lt.lbnd(nf).and.
     $                 ll.eq.rbnd(nf))then
!
! very special case 1 left of f block when rbnd<lbnb
!
                       njb=njb+1
                       njn(njb)=rbnd(nf)-tnrows(ncc-1)
                       nj1(njb)=njn(njb)-tnrows(ncc-2)
                       nj2(njb)=lbnd(nf)
		       nflag(njn(njb))=1
!
! end special case
! rbnd ex.2
!
                       njb=njb+1
                       njn(njb)=ll
                       nj1(njb)=ll+1
                       nj2(njb)=ll+tnrows(ncc)
		       lltemp=rbnd(nf)+1
		       nflag(njn(njb))=1
!
! loop up the fault plane
!
		        do 1251 nd=rbnd(nf)+1,lbnd(nf)
		          if(nd.eq.lbnd(nf))then
			     njb=njb+1
                             njn(njb)=nd
			     nflag(njn(njb))=1
                             nj1(njb)=rbnd(nf)-tnrows(ncc-1)
                             nj2(njb)=nd-1
			  else
		             njb=njb+1
                             njn(njb)=nd
			     nflag(njn(njb))=1
                             nj1(njb)=nd+1
                             nj2(njb)=nd-1
	                  endif
		          lltemp=lltemp+1
 1251                    continue

		     else if(rbnd(nf).gt.lbnd(nf).and.
     $                 ll.eq.lbnd(nf))then

                       njb=njb+1
                       njn(njb)=lbnd(nf)+tnrows(ncc)
                       nj1(njb)=rbnd(nf)
                       nj2(njb)= njn(njb)+tnrows(ncc+1)
                       nflag(njn(njb))=1
		       njb=njb+1
                       njn(njb)=ll
                       nj1(njb)=ll-tnrows(ncc-1)
                       nj2(njb)=ll+1
		       lltemp=lbnd(nf)+1
		        do 1252 nd=lbnd(nf)+1,rbnd(nf)
		           if(nd.eq.rbnd(nf))then

			     njb=njb+1
                             njn(njb)=lltemp
                             nj1(njb)=lltemp-1
                             nj2(njb)=lbnd(nf)+tnrows(ncc)
			   else

			     njb=njb+1
                             njn(njb)=lltemp
                             nj1(njb)=lltemp-1
                             nj2(njb)=lltemp+1
			   endif
		         lltemp=lltemp+1
 1252                 continue
		     endif
		endif
 125         continue

	  if(ll.lt.tnrows(1).and.nr.ne.1)then
	     njb=njb+1
             njn(njb)=ll
	     nflag(njn(njb))=1
             nj1(njb)=ll+1
             nj2(njb)=ll-1
	   endif

           if(nr.eq.1.and.ncc.ne.1.and.ncc.ne.tncols
     $    .and.nflag(ll).eq.0.and.sflag.ne.1.or.sflag.eq.2)then
             njb=njb+1
             njn(njb)=ll
	     nflag(njn(njb))=1
             nj1(njb)=ll-tnrows(ncc-1)
             nj2(njb)=ll+tnrows(ncc)
	   endif

     	   if(ncc.eq.tncols.and.nr.ne.1.and.nr.ne.tnrows(ncc))then
	     njb=njb+1
             njn(njb)=ll
	     nflag(njn(njb))=1
             nj1(njb)=ll-1
             nj2(njb)=ll+1
	   endif

 115      continue
 105     continue


      do 85 l=1,nflt
       xdif1 =  dabs(x(nn)-xfault(l))
        if(xdif1.le.1.0e-3.and.((nn.ne.lbnd(l)).or.
     $  (nn.ne.rbnd(l))))then

            if(z(lbnd(l)).lt.z(rbnd(l)))then
               nbn(count)=rbnd(l)
             end if
          nbn(count+1)=lbnd(l)
          count=count+1
        endif

        if((nn.eq.lbnd(l)).and.(nn.eq.rbnd(l)))then
          nbnn=nbnn-1
        end if
   85 continue
      count=count+1
      m = m - 1
   90 continue

  503 format (1h1,//,3x,
     $ '**** warning nhband g.t. array dimensions ***',
     $ /,3x,'nhband = ',i5,5x,'aborting simulation')
  208 return
      end
!***********************************************************************
!
      subroutine subsid (nrows,ncols,vs,zsp,dt,icase,ver,
     $	                 head,rotflag,vzt)
!
!  this subroutines moves the mesh due to tectonic subsidence and
!  adjusts the heads to account for changes in potential energy due
!  to uplift (icase.eq.2)
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
       real*8 vs(maxnodes),zsp(maxnodes),ver(maxbnds),dt,head(maxnodes)
       real*8 vzt(maxbnds)

      integer nrows(maxbnds),ncols,icase(maxbnds),rotflag
!
!  icase = 0    pinned top node with subsidence/uplift, no grid deformation
!  icase = 1    pinned top node with subsidence/uplift with grid deformation
!  icase = 2    subsidence/uplift with unpinned top node
!  icase = 3    uplift with erosion applied to top node
!
      nn = 0
	
      do 600 n=1,ncols
       ii=1
       if(icase(n).gt.1) ii = 0
        nn = nn + nrows(n)
          do 400 m=1,nrows(n) - ii
          if(icase(n).eq.0) nm = nn - m
          if(icase(n).eq.1) nm = nn - m
          if(icase(n).eq.2) nm = nn - m + 1
          if(icase(n).eq.3) nm = nn - m + 1
          if(rotflag.eq.0.and.icase(n).eq.0)  zsp(nm) = zsp(nm)
     $    + vzt(n)*dt
          if(rotflag.eq.0.and.icase(n).eq.1)  zsp(nm) = zsp(nm)
     $    + vs(nm)*dt
          if(rotflag.eq.1)  zsp(nm) = zsp(nm) + vzt(n)*dt
         

         if(icase(n).eq.2) then
           head(nm) = head(nm) + vzt(n)*dt
           zsp(nm) = zsp(nm) + vzt(n)*dt
         endif

	 if(icase(n).eq.3) then
            write(8,*) 'nm,head,ver*dt',nm,head(nm),ver(n)*dt
	    head(nm) = head(nm)-ver(n)*dt
            write(8,*) 'head(nm)=',head(nm)
	 endif
400   continue
        if(icase(n).eq.3) then
        zsp(nn) = zsp(nn) - ver(n)*dt
	head(nn) = head(nn) - ver(n)*dt
       else
        continue
      end if
  600 continue
  111  format ('nn=',i5,'zsp=',1pe12.4,'ver=',1pe12.4,
     $  'icase=',i5)
   10	format(3(1pe15.6,5x))
      return
      end
!***********************************************************************
!
      subroutine gensis (nrows,ncols,nnode,xsp,zsp,vs,
     $ phi,phiold,head,temp,tempc,tti,ttic,told,delz,
     $ icase,zcheck,nrold,it,grdfac,nslp,
     $ delta,omega,vzt,dt,node_angle,sfcel,delsfc,conc,presl,
     $ prslold,rotflag,sigmax,phimin)
!
!  this subroutine adds/removes nodes due to subsidence/erosion and renumbers
!  nodal variables accordingly
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 xsp(maxnodes),zsp(maxnodes),xold(maxnodes),zold(maxnodes)
      real*8 vs(maxnodes),vsold(maxnodes),phi(maxnodes),phiold(maxnodes)
      real*8 head(maxnodes),temp(maxnodes),tempc(maxnodes),tti(maxnodes)
      real*8 told(maxnodes),zcheck(maxbnds),grdfac(maxbnds),vzt(maxbnds)
      real*8 phild1(maxnodes),phild2(maxnodes),hold(maxnodes)
      real*8 tld1(maxnodes),tcold(maxnodes),ttiold(maxnodes)
      real*8 tticld(maxnodes),tld2(maxnodes),ttic(maxnodes)
      real*8 delz(maxbnds), delta, omega
      real*8 node_angle(maxnodes),n_angleold(maxnodes)
      real*8 sfcel(maxnodes),lsfcelold(maxnodes)
      real*8 delsfc(maxnodes),delsfcold(maxnodes)
      real*8 conc(maxnodes),concld(maxnodes)
      real*8 presl(maxnodes),prslold(maxnodes)
      real*8 prsltmp1(maxnodes),prsltmp2(maxnodes)
      real*8 sigmax(maxnodes),sigmax_old(maxnodes),phimin(maxnodes)
      real*8 phimin_old(maxnodes)
      real*8 dcosd, dsind
      real*4 theta1


      integer nslp(maxbnds),m1,m2,m3,nrold(maxbnds),nrows(maxbnds)
      integer icase(maxbnds),rotflag

      if(it.eq.0) then
      nnold = 0
      nnode = 0
      do 901 n=1,ncols
      nrold(n) = nrows(n)
      nnode = nnode + nrows(n)
      nnold = nnold + nrold(n)
  901 continue
      endif
      if(it.eq.0) return
!
!  see where we need to add or eliminate nodes from a column
!
05    ll = 0
      do 903 n=1,ncols
      do 902 l=1,nrows(n)
      ll = ll + 1
      xold(ll) = xsp(ll)
      zold(ll) = zsp(ll)
      vsold(ll) = vs(ll)
      phild1(ll) = phi(ll)
      phild2(ll) = phiold(ll)
      hold(ll) = head(ll)
      concld(ll) =conc(ll)
      tld1(ll) = temp(ll)
      tld2(ll) = told(ll)
      tcold(ll) = tempc(ll)
      ttiold(ll) = tti(ll)
      tticld(ll) = ttic(ll)
      n_angleold(ll)= node_angle(ll)
      lsfcelold(ll)= sfcel(ll)
      delsfcold(ll)= delsfc(ll)
      prsltmp1(ll)= presl(ll)
      prsltmp2(ll)= prslold(ll)
      sigmax_old(ll)=sigmax(ll)
      phimin_old(ll)=phimin(ll)
  902 continue
  903 continue



      nnold = 0
      nnode = 0
      do 701 n=1,ncols
        nrold(n) = nrows(n)

          if(icase(n).le.1.and.zcheck(n).gt.delz(n)+1.0d-0)then
           nrows(n) = nrows(n) + 1
	   endif
	    if(zcheck(n).le.1.0d-0)  nrows(n)=nrows(n)-1


        nnode = nnode + nrows(n)
        nnold = nnold + nrold(n)
  701 continue

        if(nnode.gt.(nnold+ncols))then
       		write(*,*)'you are subsiding and are in trouble'
       		write(*,*)'stop3'
       		stop
        endif


      if(nnode.gt.maxnodes) write (8,501) nnode
      if(nnode.gt.maxnodes) then
        write(*,*)'stop4'
        stop
       endif

      mm = nnode
      mmm = nnold


      do 101 n=ncols,1,-1

      m = nrows(n) - nrold(n)
!
!  take action if there is a new node in this column
! if m is > 0 then the new column is > old column.
!
      if(m.gt.0) then

      m1 = mm
      m2 = mmm
!
! if there is a new node then the value for the new
! node is givenfrom the oldnode
!
      xsp(m1) = xold(m2)
      zsp(m1) = zold(m2)
      vs(m1) = vsold(m2)
      phi(m1) = phild1(m2)
      phiold(m1) = phild2(m2)
      head(m1) = hold(m2)
      conc(m1) = concld(m2)
      temp(m1) = tld1(m2)
      told(m1) = tld2(m2)
      tempc(m1) = tcold(m2)
      tti(m1) = ttiold(m2)
      ttic(m1) = tticld(m2)
      node_angle(m1)=0.0d0

      sfcel(m1)=lsfcelold(m2)
      delsfc(m1)=delsfcold(m2)
      presl(m1) = prsltmp1(m2)
      prslold(m1)= prsltmp2(m2)
      sigmax(m1) = sigmax_old(m2)
      phimin(m1) = phimin_old(m2)

      if(rotflag.eq.2) then

      zsp(m1-1) = zold(m2-1) -  delz(n) * dcosd(delta+omega)   
      xsp(m1-1) = xold(m2-1)  + (delz(n)*dsind(delta+omega))

       else
      zsp(m1-1)=zold(m2-1)+(zold(m2)-zold(m2-1))*0.99
      if(zsp(m1-1).gt.zold(m2)) zsp(m1-1) = zold(m2) - 1.0d+0
      xsp(m1-1)=xold(m2-1)
      endif
      m1 = m1 - 2
      m2 = m2 - 1

      do 201 l=1,nrows(n)-2

      xsp(m1) = xold(m2)
      zsp(m1) = zold(m2)
      vs(m1) = vsold(m2)
      phi(m1) = phild1(m2)
      phiold(m1) = phild2(m2)
      head(m1) = hold(m2)
      conc(m1) = concld(m2)
      temp(m1) = tld1(m2)
      told(m1) = tld2(m2)
      tempc(m1) = tcold(m2)
      tti(m1) = ttiold(m2)
      ttic(m1) = tticld(m2)
      node_angle(m1)=n_angleold(m2)
      sfcel(m1)=lsfcelold(m2)
      delsfc(m1)=delsfcold(m2)
      presl(m1) = prsltmp1(m2)
      prslold(m1)= prsltmp2(m2)
      sigmax(m1) = sigmax_old(m2)
      phimin(m1) = phimin_old(m2)

      m1 = m1 - 1
      m2 = m2 - 1
  201 continue
!
!  interpolate new heads, concenrations, temperatures, etc.
!  for the new  node that we just generated
!
      dlz = zsp(mm)-zsp(mm-2)
      dlzz = zsp(mm-1) - zsp(mm-2)
      vs(mm-1) = vs(mm-2)
      phi(mm-1) = phi(mm-2) + dlzz*(phi(mm)-phi(mm-2))/dlz
      phiold(mm-1) = phiold(mm-2) + dlzz*(phiold(mm)-phiold(mm-2))/dlz
      head(mm-1) = head(mm-2) + dlzz*(head(mm)-head(mm-2))/dlz
      conc(mm-1) = conc(mm-2) + dlzz*(conc(mm)-conc(mm-2))/dlz
      temp(mm-1) = temp(mm-2) + dlzz*(temp(mm)-temp(mm-2))/dlz
      told(mm-1) = told(mm-2) + dlzz*(told(mm)-told(mm-2))/dlz
      tempc(mm-1) = tempc(mm-2) + dlzz*(tempc(mm)-tempc(mm-2))/dlz
      tti(mm-1) = tti(mm-2) + dlzz*(tti(mm)-tti(mm-2))/dlz
      ttic(mm-1) = ttic(mm-2) + dlzz*(ttic(mm)-ttic(mm-2))/dlz
      node_angle(mm-1)=0.d0
      n_angleold(mm-1)=0.d0
      sfcel(mm-1)=sfcel(mm-2)
      lsfcelold(mm-1)=lsfcelold(mm-2)
      delsfc(mm-1)=delsfc(mm-2)
      delsfcold(mm-1)=delsfcold(mm-2)
      presl(mm-1) = presl(mm-2) + dlzz*(presl(mm)-presl(mm-2))/dlz
      prslold(mm-1) = prslold(mm-2)
     $   + dlzz*(prslold(mm)-prslold(mm-2))/dlz
     
      sigmax(mm-1) = sigmax(mm-2) + dlzz*(sigmax(mm)-sigmax(mm-2))/dlz
      phimin(mm-1) = phimin(mm-2) + dlzz*(phimin(mm)-phimin(mm-2))/dlz

      else
      continue
      end if
!
!  take action if there are no new nodes in this column
!
      if(m.eq.0.and.mm.ne.mmm) then

      m1 = mm
      m2 = mmm

      do 202 l=1,nrows(n)

      xsp(m1) = xold(m2)
      zsp(m1) = zold(m2)
      vs(m1) = vsold(m2)
      phi(m1) = phild1(m2)
      phiold(m1) = phild2(m2)
      head(m1) = hold(m2)
      conc(m1) =  concld(m2)
      temp(m1) = tld1(m2)
      told(m1) = tld2(m2)
      tempc(m1) = tcold(m2)
      tti(m1) = ttiold(m2)
      ttic(m1) = tticld(m2)
      node_angle(m1)=n_angleold(m2)
      sfcel(m1)=lsfcelold(m2)
      delsfc(m1)=delsfcold(m2)
      presl(m1) = prsltmp1(m2)
      prslold(m1)= prsltmp2(m2)
      sigmax(m1) = sigmax_old(m2)
      phimin(m1) = phimin_old(m2)

      m1 = m1 - 1
      m2 = m2 - 1
  202 continue

      else
      continue
      end if

      if(m.lt.0) then

      m1 = mm
      m2 = mmm
      m3 = m2-1

      xsp(m1) = xold(m2)
      zsp(m1) = zold(m2)
      vs(m1) = vsold(m2)
      head(m1) = hold(m2)
      conc(m1) = concld(m2)
      temp(m1) = tld1(m2)
      told(m1) = tld2(m2)
      tempc(m1) = tcold(m2)
      tti(m1) = ttiold(m3)
      ttic(m1) = tticld(m3)
      node_angle(m1)=n_angleold(m3)
      sfcel(m1)=lsfcelold(m3)
      delsfc(m1)=delsfcold(m3)
      presl(m1) = prsltmp1(m2)
      prslold(m1)= prsltmp2(m2)
      phimin(m1) = phimin_old(m2)
      sigmax(m1) = sigmax_old(m2)
      phi(m1) = phild1(m2)
      phiold(m1) = phild2(m2)

      m1 = mm - 1
      m2 = mmm - 2

      do 203 l=1,nrows(n) - 1

      xsp(m1) = xold(m2)
      zsp(m1) = zold(m2)
      vs(m1) = vsold(m2)
      phi(m1) = phild1(m2)
      phiold(m1) = phild2(m2)
      head(m1) = hold(m2)
      conc(m1) = concld(m2)
      temp(m1) = tld1(m2)
      told(m1) = tld2(m2)
      tempc(m1) = tcold(m2)
      tti(m1) = ttiold(m2)
      ttic(m1) = tticld(m2)
      node_angle(m1)=n_angleold(m2)
      sfcel(m1)=lsfcelold(m2)
      delsfc(m1)=delsfcold(m2)
      presl(m1) = prsltmp1(m2)
      prslold(m1)= prsltmp2(m2)
      sigmax(m1) = sigmax_old(m2)
      phimin(m1) = phimin_old(m2)

      m1 = m1 - 1
      m2 = m2 - 1
  203 continue

      else
      continue
      end if

      mm = mm - nrows(n)
      mmm = mmm - nrold(n)
101   continue

  501 format (1h1,//,3x,
     $ '**** warning nnode g.t. array dimensions ***',
     $ /,3x,'nnode = ',i7,5x,'aborting simulation')
 1122 format('node=',i5,2x,'nodeold=',i5,2x,'delta=',
     $ 1pe12.4,2x,'omega=',1pe12.4,/,'xold=',1pe12.4,2x,
     $ 'zold=',1pe12.4,2x,'xsp=',1pe12.4,2x,'zsp=',1pe12.4,2x,
     $ 'delz=',1pe12.4)


      return
      end
!***********************************************************************
!
      subroutine grdchk (ncols,nrows,zsp,zcheck,delz,icheck,grdfac,xsp
     $  ,it,vzt,dt,rotflag,icase,vs)
!
!  determine whether its time to add a new node and set flag to update
!  the mesh calling gensis and mesh
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 delz(maxbnds),vzt(maxbnds),zsp(maxnodes),xsp(maxnodes)
      real*8 grdfac(maxbnds),zcheck(maxbnds),vs(maxnodes)


      integer it,rotflag,nrows(maxbnds),icase(maxbnds)

      nn = 0
      do 601 n=1,ncols
      nn = nn + nrows(n)

       if(rotflag.eq.2) then
      zcheck(n) = ( (zsp(nn)-zsp(nn-1))**2.0d+00
     $         + (xsp(nn)-xsp(nn-1))**2.0d+00) **0.5
      if(zsp(nn).lt.zsp(nn-1))then
          zcheck(n)=-zcheck(n)
      endif
      else
        zcheck(n) = zsp(nn)-zsp(nn-1)
      endif

  601 continue
      icheck = 1
      nn = 0
      do 50 n=1,ncols
      nn = nn + nrows(n)
        if(it.gt.1.and.rotflag.eq.0)then
	       delz(n)=dabs(vzt(n)*grdfac(n)*dt)
!               write(8,*) 'delz=',delz(n),'grdfac=',grdfac(n)
	endif
! m.a.p. 12/26/99
!    if(icase(n).eq.1.and.zcheck(n).gt.(delz(n)+1.0d-0)) icheck = 0

      if(rotflag.ne.2) then
      if(icase(n).le.1.and.zcheck(n).gt.(delz(n)+1.0d-0)) icheck = 0
      if(zcheck(n).le.1.0d-0) icheck = 0
      else
      if(zcheck(n).gt.delz(n)) icheck = 0
      if(zcheck(n).le.0.0d+0) icheck = 0
      endif

   50 continue
      return
      end
!***********************************************************************
!
      subroutine areas (x,z,ni,nj,nk,area,nelem,al,b,c,nnode,time,
     $  tnrows,tncols,grdfac,vs,phi,phiold)
!
!     this subroutine calculates the triangular element areas
!    for the finite element mesh & and checks the integrity
!     of the input data set.
!
!***********************************************************************
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 x,z,area,b,c,time,al,vs,phi,phiold

      integer ni,nj,nk,nelem,nnode,flag,tncols,tnrows

      dimension x(maxnodes),z(maxnodes),area(maxelems)
      dimension b(6,maxelems),c(6,maxelems)
      dimension ni(maxelems),nj(maxelems),nk(maxelems)
      dimension tnrows(maxbnds)
      dimension al(6,maxelems),grdfac(maxbnds)
      dimension vs(maxnodes),phi(maxnodes),phiold(maxnodes)

      flag = 0
      if (nelem.gt.maxelems.or.nelem.lt.1) then
            write (8,115 ) nelem

      endif

      do 100 m=1,nelem

      i=ni(m)
      j=nj(m)
      k=nk(m)
      if(i.eq.0.and.j.eq.0.and.k.eq.0) go to 100
      if (i.gt.nnode.or.i.lt.1) write (*,116) nnode,nelem,i,j,k
      if (j.gt.nnode.or.j.lt.1) write (*,116) nnode,nelem,i,j,k
      if (k.gt.nnode.or.k.lt.1) write (*,116) nnode,nelem,i,j,k
!
!     compute derivatives
!

      al(1,m) = x(j)*z(k) - x(k)*z(j)
      al(2,m) = x(k)*z(i) - x(i)*z(k)
      al(3,m) = x(i)*z(j) - x(j)*z(i)

      b(1,m)= z(j)-z(k)
      b(2,m)= z(k)-z(i)
      b(3,m)= z(i)-z(j)
      c(1,m)= x(k)-x(j)
      c(2,m)= x(i)-x(k)
      c(3,m)= x(j)-x(i)
!
!     calculate element area & check for negative area
!
      area(m) = 0.5d+00*((x(i)*z(j)-x(j)*z(i))
     $  + (x(k)*z(i)-x(i)*z(k)) + (x(j)*z(k)-x(k)*z(j)))

      if (area(m).gt.0.0d+00) go to 100
      write (8,10) m,i,j,k,area(m)
      write(8,1101)x(i),z(i),x(j),z(j),x(k),z(k)
      flag = 1
  100 continue
!
!     if we have any negative element areas (ie. flag=1) exit program
!
      if (flag.eq.0) go to 101
      write (8,11) time
      write (8,63) nnode,nelem
      write (8,65)
      do 201 i = 1, tncols
      write (8,66) i, tnrows(i), grdfac(i)
201   continue
      write (8,62)
      write (8,20) (l,x(l),l=1,nnode)
      write (8,61)
      write (8,20) (l,z(l),l=1,nnode)
      write (8,162)
      write (8,20) (l,vs(l),l=1,nnode)
      write (8,161)
      write (8,20) (l,phi(l),l=1,nnode)
      write (8,163)
      write (8,20) (l,phiold(l),l=1,nnode)
  101 continue
      if (flag.eq.1) then
         write(*,*)'negative area cacluated, stop'
	 stop
       endif

   10 format (/,1x,'***warning***  negative area computed for',
     $ ' element # ',i4,/,1x,'i,j,k: ',3(i4,2x),5x,'area: ',
     $  1pe12.3)
   11 format (1h1,///,'  **** simulation aborting ****',
     $ /,4x,'time (yr) = ',1pe12.1)
   20 format (6(5x,i4,1x,1pe11.3))
   61 format (1h1,///,3x,'z coordinates',/,6x,'(m)',///)
   62 format (1h1,///,3x,'x coordinates',/,6x,'(m)',///)
   63 format (///,3x,'nnode = ',i4,5x,'nelem=',i4)
   65 format (1h1,///,3x,'no. of nodes per column',///, 1x, 'col #',
     $ 6x, 'no. of nodes', 6x, 'grdfac')
   66 format (3x, i4, 4x, i5, 2x,1pe12.3)
  161 format (1h1,///,3x,'phi ',/,6x,'(m^3/m^3)',///)
  162 format (1h1,///,3x,'vs',/,6x,'(m/yr)',///)
  163 format (1h1,///,3x,'phiold ',/,6x,'(m^3/m^3)',///)

 1101 format (3(/,'x=',1pe12.3,3x,'z=',1pe12.3))

  115 format (5x,'****** warning nelem out of range in sub. areas',
     $ //,' nelem: ',i5)
  116 format (5x,'****** warning i,j,k out of range in sub. areas',
     $ //,' nnode,nelem,i,j,k: ',5i5)
      return
      end
!***********************************************************************
      subroutine febnd(aa,bb,head,nh,nhb,nnode,iout,nfband,it,ntime,
     $   hdbase,tplsh)
!
!     this subroutine adjusts the aa matrix and bb vector to
!     account for constant head nodes
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 aa,bb,head,hdbase,tplsh(maxnodes)
      integer nh,nhb,nnode,iout,nfband,it,ntime
      dimension aa(maxnodes,maxwdt),bb(maxnodes)
      dimension head(maxnodes),nh(maxbnds),hdbase(maxbnds)
!
!   alter aa and bb for constant head nodes, keeping symmetry
!

      do 100 n=1,nhb
      if(tplsh(nh(n)).lt.0) go to 100
      bb(nh(n))=head(nh(n))
!            bb(nh(n))=hdbase(n)
      aa(nh(n),1)=1.0d0

      do 140 l=2,nfband
      jl=nh(n)+l-1
      if(jl.gt.nnode)go to 130
      bb(jl)=bb(jl)-aa(nh(n),l)*bb(nh(n))
  130 aa(nh(n),l)=0.0d0
      jk=nh(n)-l+1
      if(jk.lt.1)go to 140
      bb(jk)=bb(jk)-aa(jk,l)*bb(nh(n))
      aa(jk,l)=0.0d0
  140 continue
  100 continue
      if (iout.ne.2) go to 300
      write(8,10)
      do 250 l=1,40
      write (8,20) (aa(l,m),m=1,20)
  250 continue
      write (8,30) (bb(l),l=1,nnode)
  300 continue
   10 format (/,1x,'subroutine febnd',/,1x,'aa matrix:')
   20 format (11(1pe12.3))
   30 format (//,1x,'b vector',10(1pe12.3))
      return
      end
!***********************************************************************
      subroutine hebnd(aa,bb,temp,nt,ntb,nnode,iout,nhband,it,ntime,
     $  tplsh)
!
!     this subroutine adjusts the aa matrix and bb vector to
!     account for constant temperature nodes
!
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      integer nt(maxbnds),ntb,nnode,iout,nhband,it,ntime
      real*8 aa(maxnodes,maxwdt),bb(maxnodes),temp(maxnodes)
      real*8  tplsh(maxnodes)
!
!   alter aa and bb for constant temp nodes, keeping symmetry
!
!
      ij=(nhband-3)/2
      ij2=ij+2

      do 100 n=1,ntb
      if(tplsh(nt(n)).le.0.0) temp(nt(n))= 0.0 
      if(tplsh(nt(n)).gt.0.0) temp(nt(n)) = tplsh(nt(n))
      bb(nt(n))=temp(nt(n))
      aa(nt(n),ij2)=1.0d0
      do 140 l=1,nhband
      if (l.eq.ij2) go to 140
      aa(nt(n),l)=0.0d0
  140 continue
  100 continue
!
!     print out aa matrix if iout>1
!
      if (iout.ne.2.or.it.ne.ntime) go to 300
      write(8,10)
      do 250 l=1,nnode
      write (8,20) (aa(l,m),m=1,nhband)
  250 continue
      write (8,30) (bb(l),l=1,20)
  300 continue
   10 format (/,1x,'subroutine hebnd',/,1x,'aa matrix:')
   20 format (11(1pe12.3))
   30 format (//,1x,'b vector',10(1pe12.3))
      return
      end
c
c
!
!***********************************************************************
      subroutine hnbnd(x,z,bb,nj1,nj2,njn,j1,j2,njb,it,iout,ntime)
!
!     this subroutine evaluates the heat flux term, jflux on the lower
!     and side boundaries of the basin
!
!
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
!
      integer njb,iout,nj1(maxbnds),nj2(maxbnds),njn(maxbnds),it,ntime
      real*8 bb(maxnodes), j1(maxbnds),j2(maxbnds),x(maxnodes)
      real*8 z(maxnodes),l1,l2,jflux
!
!  evaluate flux vector q
!
      if(njb.eq.0) go to 200
      if(iout.eq.2.or.it.eq.ntime) write (8,10) it
      do 100 n=1,njb

        l1 = x(njn(n))-x(nj1(n))

        l2 = x(nj2(n))-x(njn(n))

        if((j1(n).le.0.0d+00).or.(j2(n).le.0.0d+00))then
            j1(n)=j1(1)
            j2(n)=j2(1)
          endif

      jflux = j1(n)*l1/2.d+00 + j2(n)*l2/2.d+00

       bb(njn(n)) = bb(njn(n)) + jflux

  100 continue
   10 format (//,2x,'natural boundary condition for time step:',i4,//)
   20 format (2x,'njn:',i4,2x,'nj1:',i4,2x,'nj2:',i4,3x,'j:',
     $1pe12.3,3x,'l1:',1pe12.3,3x,'l2:',1pe12.3,3x,1pe12.3,3x,i4,
     $3x,1pe12.3,3x,1pe12.3,3x,1pe12.3)
!
!     okay were done
!
  200 return
      end
!***********************************************************************
      subroutine output (x,z,head,temp,tti,qx,qz,nnode,it,iprint,dt,
     $ nelem,ipr,conc,iheat,ibrine,time,icoup,etemp,
     $ ni,nj,nk,nl,mat,nbn,nbnn,tempc,ttic,phi,kx,kz,rvis,
     $ tnrows,rhof,tr,masoil,masgas,rv,oilvol, gasvol,persat,mash2o,
     $ masco2,tncols,kgoil,kggas,iskip,grdfac,ber,bec,nflt,ntime,
     $ sigmae,pres,presl,ss,vsm,heado,headg,rhoo,rhog,ioil)
!
!
!     this subroutine prints output to unit 8
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 x(maxnodes),z(maxnodes),head(maxnodes),qx(maxelems)
      real*8 tti(maxnodes),conc(maxnodes),kx(maxelems),kz(maxelems)
      real*8 etemp(maxelems),tr(maxelems),rv(maxelems),masoil(maxelems)
      real*8 oilvol(maxelems),gasvol(maxelems),persat(maxelems)
      real*8 kgoil(maxelems),kggas(maxelems),tempc(maxnodes)
      real*8 sigmae(maxnodes),pres(maxnodes),presl(maxnodes)
      real*8 heado(maxnodes),headg(maxnodes),rhoo(maxnodes)
      real*8 qz(maxelems),rvis(maxnodes),rhof(maxnodes),temp(maxnodes)
      real*8 phi(maxnodes),masgas(maxelems),grdfac(maxbnds)
      real*8 mash2o(maxelems),masco2(maxelems),ttic(maxnodes)
      real*8 ss(maxnodes),vsm(maxnodes),rhog(maxnodes)

      integer nnode,it,iprint,ibrine,iheat,icoup,nl(maxelems),iskip,ipr
      integer nbnn,nbn(maxbnds),ni(maxelems),nj(maxelems),nk(maxelems)
      integer ioil,tncols,tnrows(maxbnds)
      integer ber(maxbnds),bec(maxbnds),nflt,ntime,mat(maxelems)


      if(ipr.ne.iprint) go to 300
      ipr = 0
      write (8,10) time
      write (8,63) nnode,nelem
      write (8,165)
      do 201 i = 1, tncols
         write (8,166) i, tnrows(i), grdfac(i)
201   continue

      write (8,51)
       write (8,20) (l,head(l),l=1,nnode,iskip)
      if(iheat.eq.1) write (8,50)
      if(iheat.eq.1) write (8,20) (l,temp(l),l=1,nnode,iskip)
      if(ibrine.eq.1) write (8,70)
      if(ibrine.eq.1) write (8,20) (l,conc(l),l=1,nnode,iskip)
      if(ibrine.eq.1.or.iheat.eq.1) write (8,79)
      if(ibrine.eq.1.or.iheat.eq.1) write (8,20)
     $   (l,rhof(l),l=1,nnode,iskip)
       write (8,30)
       write (8,40) (l,qx(l),qz(l),l=1,nelem,iskip*2)
       write (8,151)
       write (8,20) (l,phi(l),l=1,nnode,iskip)
       write (8,152)
       write (8,20) (l,pres(l)/1.d+6,l=1,nnode,iskip)
       write (8,153)
       write (8,20) (l,sigmae(l)/1.d+6,l=1,nnode,iskip)
       write (8,154)
       write (8,20) (l,vsm(l),l=1,nnode,iskip)
         if(ioil.eq.1) then
         write (8,161)
         write (8,20) (l,heado(l),l=1,nnode,iskip)
         write (8,162)
         write (8,20) (l,headg(l),l=1,nnode,iskip)
         write (8,163)
         write (8,20) (l,rhoo(l),l=1,nnode,iskip)
         write (8,164)
         write (8,20) (l,rhog(l),l=1,nnode,iskip)
         write (8,82)
         write (8,83) (l,tr(l),l=1,nelem,iskip*2)
         write (8,186)
         write (8,85) (l,persat(l),l=1,nelem,iskip*2)
         write (8,87)
         write (8,85) (l,rv(l),l=1,nelem,iskip*2)
         endif


c
   10 format (1h1,///,'  **** solution output ****',/,4x,'time (yr) = '
     $ ,1pe12.1)
   17 format (132i1,/,132i1)
   20 format (5(5x,i4,1x,1pe13.6))
   21 format (10(1pe13.6))
   30 format (1h1,///,'   elemental darcy velocities: ',
     $      /,5x,'(qx & qz in m/year)',///)
   31 format (1h1,///,'   elemental hydraulic conductivities: ',
     $      /,5x,'(qx & qz in m/year)',///)
   40 format (4(5x,i4,1x,1pe11.3,1x,1pe11.3))
   41 format (5(3x,1pe11.3,1x,1pe11.3))
   42 format (1x, 1pe12.2, 5x, 1pe12.3)
   50 format (1h1,///,3x,'nodal temperatures:',/,6x,'(centigrade)',///)
   51 format (1h1,///,10x,'nodal hydraulic heads:',/,15x,
     $ ' (meters)',///)
   53 format (1h1,///,10x,'conductive nodal temperatures:',/,15x,
     $ ' (meters)',///)
   55 format (1h1,///,3x,'nodal fluid density:',/,6x
     $ ,'(kg/m^3)',///)
   60 format (1h1,///,3x,'nodal maturation:',/,6x,'(tti)',///)
   61 format (1h1,///,3x,'updated z coordinates:',/,6x,'(m)',///)
   62 format (1h1,///,3x,'updated x coordinates:',/,6x,'(m)',///)
   63 format (///,3x,'nnode = ',i4,5x,'nelem=',i4)
   64 format (1h1,///,3x,'update connectivity matrix:',////)
   65 format (5(1x,'m=',i4,3(1x,i4),4x))
   73 format (3(3x,i8))
   75 format (20(2x,i4))
   66 format (////,3x,'number of boundary nodes:',i4)
   67 format (10(3x,i5))
   70 format (1h1,///,3x,'nodal brine concentrations',/,
     $ 6x,'(mass fraction):',///)
   71 format (1h1,///,3x,'nodal viscosity data:',/,6x,'(m)',///)
   79 format (1h1,///,3x,'nodal fluid density',/,
     $ 6x,'(kg/m**3):',///)
   80 format (1h1,///,3x,'elemental temperatures:',/,6x,'(centigrade)',
     $  ///)
   81 format (5(5x,i4,1x,1pe13.6))
   82 format (1h1,///,3x,'elemental maturation: transformation ratio',
     $  /,6x,'(first-order rate kinetics)',///)
   83 format (5(5x,i4,1x,1pe13.6))
   84 format (1h1,///,3x,'kilograms of oil per element',
     $  /,6x,'(first-order rate kinetics)',///)
   85 format (5(5x,i4,1x,1pe13.6))
   86 format (1h1,///,3x,'kilograms of gas per element',
     $  /,6x,'(first-order rate kinetics)',///)
   87 format (1h1,///,3x,'vitrinite reflectance (rv)',
     $  /,6x,'(for type 3 kerogen only)',///)
   90 format (1x, a)
  151 format (1h1,///,10x,'nodal porosity:',/,15x,
     $ ' (meters**3/meters**3)',///)
  152 format (1h1,///,10x,'nodal pressures:',/,15x,
     $ ' (mega-pascals)',///)
  153 format (1h1,///,10x,'nodal effective stress:',/,15x,
     $ ' (mega-pascals)',///)
  154 format (1h1,///,10x,'nodal subsidence/uplift:',/,15x,
     $ ' (m/yr)',///)
  161 format (1h1,///,10x,'nodal oil heads:',/,15x,
     $ ' (meters)',///)
  162 format (1h1,///,10x,'nodal gas heads:',/,15x,
     $ ' (meters)',///)
  163 format (1h1,///,10x,'nodal oil density:',/,15x,
     $ ' (kg/m**3)',///)
  164 format (1h1,///,10x,'nodal gas density:',/,15x,
     $ ' (kg/m**3)',///)
  165 format (1h1,///,3x,'no. of nodes per column',///, 1x, 'col #',
     $ 6x, 'no. of nodes', 6x, 'grdfac')
  166 format (3x, i4, 4x, i5, 2x,1pe12.3)
  186 format (1h1,///,3x,'percent saturation per element',///)

  301 format (1pe12.3,'	',1pe12.3,'	',1pe12.3)
  700 format (1x,'time (yrs):   ', 1pe12.3)
  701 format (1x, a)
  720 format (1x,2pe13.3,'	',2pe13.3,'	',1pe12.3)
  740 format (1x, 2pe13.3,8('	',1pe13.3))
  750 format (1x,/)
  300 return
      end
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      subroutine heat (aa,bb,x,z,ni,nj,nk,nl,tkf,tks,qx,qz,cvs,cvf,ldis,
     $tdis,phi,rhof,nelem,iout,mat,area,icond,al,b,c,rhos,
     $dt,nhband,temp,theta,phiold,eltp,itt,ntime)
!
!     this subroutine computes the stiffness matrix and forcing
!     vectors for the heat transfer equation using three and four node
!     triangular elements. three node elements use linear
!     interpolation functions summed over the element.  four node
!     elements have the fourth node positioned at any point allong
!     one side of the element and use both linear and quadratic
!     interpolation functions summed over the element.
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 s(4,4),a,b(6,maxelems)
      real*8 c(6,maxelems),qx(maxelems),qz(maxelems)
      real*8 tke,ex,ez,exz,tdis(0:maxmats)
      real*8 ldis(0:maxmats),phi(maxnodes)
      real*8 fix,cvs(0:maxmats),r(4,4)
      real*8 aa(maxnodes,maxwdt),bb(maxnodes)
      real*8 temp(maxnodes),area(maxelems),qbar,rhos(0:maxmats),phie
      real*8 c1,c2,c3,al(6,maxelems)
      real*8 dd,ddfn(8),ae,be(8),dt,rcf,theta,fv(8)
      real*8 tks(0:maxmats),qx2,qz2,tkf(0:maxmats)
      real*8 rhofe,rhof(maxnodes),cvf(0:maxmats)
      real*8 x(maxnodes),z(maxnodes),phiold(maxnodes),hw

      integer nelem,eltp(maxelems),i1,i2,i3,jj,ii,mat(maxelems)
      integer ni(maxelems),nj(maxelems),nk(maxelems),nl(maxelems),qq
      integer it,j,k,l,node(8),iout,nn,ne,me,sw
      integer nhband,icond,ij2,ntime,itt

       icond=1
      hw = 7.0d+04
      ij2 = ((nhband-3)/2)+ 2

      do 100 m = 1, nelem

        if(eltp(m).eq.6)goto100
        ae = area(m)
        i = ni(m)
        j = nj(m)
        k = nk(m)
        l = nl(m)
        qq=mat(m)
        if(i.eq.0) go to 100

        phie = (phi(i)+phi(j)+phi(k))/3.0d+00
        tke = (tkf(qq)**phie)*(tks(qq)**(1.0d+00-phie))
        rhofe = (rhof(i)+rhof(j)+rhof(k))/3.0d+00
        a = 1.0d+00/ae
        d = phie*rhofe*cvf(qq) + (1.0d+00-phie)*rhos(qq)*cvs(qq)
        e = ae*d/dt
        f = (-ae*rhofe*hw)/dt

        if(icond.eq.0)goto 125
        rhofe = (rhof(i)+rhof(j)+rhof(k))/3.0d+00
        qx2 = qx(m)*qx(m)
        qz2 = qz(m)*qz(m)
        qbar = dsqrt(qx2+qz2)
        if(qbar.eq.0.0d+00)goto 125
        ex = (rhofe*cvf(qq)*(ldis(qq)*qx2+tdis(qq)*qz2))/qbar+tke
        ez = (rhofe*cvf(qq)*(tdis(qq)*qx2+ldis(qq)*qz2))/qbar+tke
        exz = rhofe*cvf(qq)*(ldis(qq)-tdis(qq))*qx(m)*qz(m)/qbar
        fix = 1.0d+00
        goto 135
  125 continue
        ex = tke
        ez = tke
        exz = 0.0d+00
        fix = 0.0d+00
  135 continue

        rcf=rhofe*cvf(qq)*fix
!
!       assemble local stiffness matrix
!
        if(eltp(m).eq.1)goto300
          sw=4
          if(eltp(m).eq.2)then
            i1 = 1
            i2 = 2
            i3 = 3
          end if
          if((eltp(m).eq.4).or.(eltp(m).eq.5))then
            i1 = 2
            i2 = 3
            i3 = 1
            it=i
            i=j
            j=k
            k=it
          end if
          if(eltp(m).eq.3)then
            i1 = 3
            i2 = 2
            i3 = 1
            it=i
            i=k
            k=it
          end if

          node(1)=i
          node(2)=j
          node(3)=k
          node(4)=l
          b1=b(i1,m)
          b2=b(i2,m)
          b3=b(i3,m)
          g1=c(i1,m)
          g2=c(i2,m)
          g3=c(i3,m)
!
!     determination of values of shape function derivatives for
!     four node elements
!
      c2 = (al(i2,m)+b2*x(l)+g2*z(l))/(2.0d+00*ae)
      c3 = (al(i3,m)+b3*x(l)+g3*z(l))/(2.0d+00*ae)
      if(c2*c3.eq.0.00) go to 301
      c1 = 1.0d+00/(c3*c2)
!
!     calculate stiffness matrix for four node elements
!
      s(1,1)=a*(((ex*b1**2)+(ez*g1**2))/4.0d+00 + exz*(b1*g1/2.0d+00))
      s(1,2)=a*((ex*((b1/4.0d+00)*(b2-(c1*c2/3.0d+00)*(b2+b3)))+(ez*
     $((g1/4.0d+00)*(g2-(c1*c2/3.0d+00)*(g2+g3)))))+exz*(((b1/4.0d+00)
     $*(g2-(c2*c1/3.0d+00)*(g3+g2)))+((g1/4.0d+00)*(b2-(c2*c1/3.0d+00)
     $*(b3+b2)))))
      s(1,3)=a*((ex*((b1/4.0d+00)*(b3-(c1*c3/3.0d+00)*(b2+b3)))+(ez*
     $((g1/4.0d+00)*(g3-(c1*c3/3.0d+00)*(g2+g3)))))+exz*(((b1/4.0d+00)
     $*(g3-(c3*c1/3.0d+00)*(g3+g2)))+((g1/4.0d+00)*(b3-(c3*c1/3.0d+00)
     $*(b3+b2)))))
      s(1,4)=a*((ex*((b1/12.0d+00)*c1*(b3+b2)))+(ez*((g1/12.0d+00)*c1*
     $(g3+g2)))+exz*(((c1*b1/12.0d+00)*(g2+g3))+(((c1*g1)/12.0d+00)*
     $(b2+b3))))
      s(2,1)=s(1,2)
      s(2,2)=a*(0.5d+00)*(((ex*(b2**2)+ez*(g2**2))/2.0d+00)+((c1*c2/3.0
     $d+00)*(-ex*(b2*b3+b2**2)-ez*(g2**2+g2*g3)+((c1*c2/4.0d+00)*((ex*
     $(b3**2+b2*b3+b2**2)+ez*(g2**2+g3**2+g2*g3)))))))+exz*(((1/(4.0d+
     $00*ae))*(b2*g2-(c1*c2/3.0d+00)*(b2*g3+2.0d+00*b2*g2+b3*g2-(c1*c2/
     $2.0d+00)*(b3*g3+b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2))))+((1.0d+00/
     $(4.0d+00*ae))*(g2*b2-(c1*c2/3.0d+00)*(g2*b3+2.0d+00*g2*b2+g3*b2-
     $(c1*c2/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2)))))
      s(2,3)=a*(0.25d+00)*(ex*((b2*b3))+(ez*(g2*g3))-(c1*c2/3.0d+00)*
     $(ex*((b3**2+b2*b3))+(ez*(g2*g3+g3**2)))-(c1*c3/3.0d+00)*
     $(ex*((b2**2+b2*b3))+(ez*(g2*g3+g2**2)))+(((c1**2)*c2*c3)/6.0d+
     $00)*(ex*((b2**2+b2*b3+b3**2))+(ez*(g3**2+g2**2+g2*g3))))+
     $exz*(((1.0d+00/(4.0d+00*ae))*(b2*g3-(c3*c1/3.0d+00)*(b2*g3+b2*g2)-
     $(c2*c1/3.0d+00)*(b3*g3+b2*g3)+(c2*c3*c1**2/6.0d+00)*(b3*g3+b2*g2+
     $b3*g2/2.0d+00+b2*g3/2.0d+00)))+((1.0d+00/(4.0d+00*ae))*(g2*b3-(c3*
     $c1/3.0d+00)*(g2*b3+g2*b2)-(c2*c1/3.0d+00)*(g3*b3+g2*b3)+(c2*c3*c1*
     $*2/6.0d+00)*(g3*b3+g2*b2+g3*b2/2.0d+00+g2*b3/2.0d+00))))
      s(2,4)=a*((c1/12.0d+00)*(ex*((b2*b3+b2**2))+(ez*(g2**2+g2*g3))-
     $(c1*c2/2.0d+00)*(ex*((b2**2+b2*b3+b3**2))+(ez*(g3**2+g2**2+g2*g3))
     $)))+exz*(((c1/(12.0d+00*ae))*(b2*g3+b2*g2-(c2*c1/2.0d+00)*(b3*g3+
     $b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2)))+((c1/(12.0d+00*ae))*(g2*b3+g2
     $*b2-(c2*c1/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2))))
      s(3,1)=s(1,3)
      s(3,2)=s(2,3)
      s(3,3)=a*(0.5d+00)*(((ex*(b3**2)+ez*(g3**2))/2.0d+00)+((c1*c3/3.0
     $d+00)*(-ex*(b2*b3+b3**2)-ez*(g3**2+g2*g3)+((c1*c3/4.0d+00)*((ex*
     $(b3**2+b2*b3+b2**2)+ez*(g2**2+g3**2+g2*g3)))))))+exz*(((1.0d+00/
     $(4.0d+00*ae))*(b3*g3-(c1*c3/3.0d+00)*(b2*g3+2.0d+00*b3*g3+b3*g2-
     $(c1*c3/2.0d+00)*(b3*g3+b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2))))+((1.
     $0d+00/(4.0d+00*ae))*(g3*b3-(c1*c3/3.0d+00)*(g2*b3+2.0d+00*g3*b3+g
     $3*b2-(c1*c3/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2)))))
      s(3,4)=a*((c1/12.0d+00)*((ex*(b2*b3+b3**2))+(ez*(g3**2+g2*g3))-
     $(c1*c3/2.0d+00)*(ex*((b2**2+b2*b3+b3**2))+(ez*(g3**2+g2**2+g2*g3)
     $))))+exz*(((c1/(12.0d+00*ae))*(b3*g3+b3*g2-(c3*c1/2.0d+00)*(b3*g3
     $+b3*g2/2.0d+00+b2*g3/2.0d+00+b2*g2)))+((c1/(12.0d+00*ae))*(g3*b3+
     $g3*b2-(c3*c1/2.0d+00)*(g3*b3+g3*b2/2.0d+00+g2*b3/2.0d+00+g2*b2))))
      s(4,1)=s(1,4)
      s(4,2)=s(2,4)
      s(4,3)=s(3,4)
      s(4,4)=a*(((c1**2/24.0d+00)*((ex*(b2**2+b2*b3+b3**2))+(ez*(g2*g3+
     $g3**2+g2**2))))+exz*(c1**2/12.0d+00)*(b3*g3+b2*g2+b3*g2/2.0d+00+
     $b2*g3/2.0d+00))
c
c     advective terms
c
      s(1,1)=s(1,1)+((qx(m)*b1+qz(m)*g1)/6.0d+00)*rcf
      s(1,2)=s(1,2)+(qx(m)*(b2/6.0d+00-c2*c1/24.0d+00*(b3+b2))+
     $qz(m)*(g2/6.0d+00-c2*c1/24.0d+00*(g3+g2)))*rcf
      s(1,3)=s(1,3)+(qx(m)*(b3/6.0d+00-c3*c1/24.0d+00*(b3+b2))+
     $qz(m)*(g3/6.0d+00-c3*c1/24.0d+00*(g3+g2)))*rcf
      s(1,4)=s(1,4)+(qx(m)*((c1/24.0d+00)*(b3+b2))+qz(m)*((c1/24.0d+00)
     $*(g3+g2)))*rcf
      s(2,1)=s(2,1)+(qx(m)*((b1/6.0d+00)*(1-c2*c1*2))+qz(m)*((g1/6.0d+00
     $)*(1-c2*c1*2)))*rcf
      s(2,2)=s(2,2)+(qx(m)*(b2-c2*c1*(2.0d+00*b2+b3/2.0d+00+b2/4.0d+00-
     $c2*c1/2.0d+00*(b3+b2/2+b3/2+b2)))/6.0d+00+qz(m)*(g2-c2*c1*(2.0d+00
     $*g2+g3/2.0d+00+g2/4.0d+00-c2*c1/2.0d+00*(g3+g2/2.0d+00+g3/2.0d+00+
     $g2)))/6.0d+00)*rcf
      s(2,3)=s(2,3)+(qx(m)*(b3-2.0d+00*c2*c1*b3-c3*c1/2.0d+00*(b3+b2/
     $2.0d+00)+c2*c3*(c1**2)/2.0d+00*(b3+b2/2.0d+00+b3/2.0d+00+b2))/6.0d
     $+00+qz(m)*(g3-2.0d+00*c2*c1*g3-c3*c1/2.0d+00*(g3+g2/2.0d+00)+c2*c3
     $*(c1**2)/2.0d+00*(g3+g2/2.0d+00+g3/2.0d+00+g2))/6.0d+00)*rcf
      s(2,4)=s(2,4)+(qx(m)*((c1/12.0d+00)*(b3+b2/2.0d+00-c2*c1*(b3+b2/2.
     $0d+00+b3/2.0d+00+b2)))+qz(m)*((c1/12.0d+00)*(g3+g2/2.0d+00-c2*c1*
     $(g3+g2/2.0d+00+g3/2+g2))))*rcf
      s(3,1)=s(3,1)+(qx(m)*((b1/6.0d+00)*(1.0d+00-c3*c1*2.0d+00))+qz(m)*
     $((g1/6.0d+00)*(1.0d+00-c3*c1*2.0d+00)))*rcf
      s(3,2)=s(3,2)+(qx(m)*(b2-2.0d+00*c3*c1*b2-c2*c1/2.0d+00*(b2+b3/2.0
     $d+00)+c2*c3*(c1**2)/2.0d+00*(b3+b2/2.0d+00+b3/2.0d+00+b2))/6.0d+00
     $+qz(m)*(g2-2.0d+00*c3*c1*g2-c2*c1/2.0d+00*(g2+g3/2.0d+00)+c2*c3*
     $(c1**2)/2.0d+00*(g3+g2/2.0d+00+g3/2.0d+00+g2))/6.0d+00)*rcf
      s(3,3)=s(3,3)+(qx(m)*(b3-c3*c1*(2.0d+00*b3+b2/2.0d+00+b3/4.0d+00-
     $c3*c1/2.0d+00*(b3+b2/2.0d+00+b3/2.0d+00+b2)))/6.0d+00+qz(m)*(g3-c3
     $*c1*(2.0d+00*g3+g2/2+g3/4.0d+00-c3*c1/2.0d+00*(g3+g2/2.0d+00+
     $g3/2.0d+00+g2)))/6.0d+00)*rcf
      s(3,4)=s(3,4)+(qx(m)*((c1/12.0d+00)*(b2+b3/2.0d+00-c3*c1*
     $(b3+b2/2.0d+00+b3/2.0d+00+b2)))+qz(m)*((c1/12.0d+00)*(g2+g3/2.0d+
     $00-c3*c1*(g3+g2/2+g3/2.0d+00+g2))))*rcf
      s(4,1)=s(4,1)+(qx(m)*(c1*b1/3.0d+00)+qz(m)*(c1*g1/3.0d+00))*rcf
      s(4,2)=s(4,2)+(qx(m)*((c1/6.0d+00)*(2.0d+00*b2-c2*c1/2.0d+00*(b3+
     $b2/2.0d+00+b3/2.0d+00+b2)))+qz(m)*((c1/6.0d+00)*(2.0d+00*g2-c2*
     $c1/2.0d+00*(g3+g2/2.0d+00+g3/2.0d+00+g2))))*rcf
      s(4,3)=s(4,3)+(qx(m)*((c1/6.0d+00)*(2.0d+00*b3-c3*c1/2.0d+00*(b3+
     $b2/2.0d+00+b3/2.0d+00+b2)))+qz(m)*((c1/6.0d+00)*(2.0d+00*g3-c3*c1
     $/2.0d+00*(g3+g2/2.0d+00+g3/2.0d+00+g2))))*rcf
      s(4,4)=s(4,4)+(qx(m)*((c1**2*(b2+b2/2.0d+00+b3+b3/2.0d+00))/12.0d+
     $00)+qz(m)*((c1**2*(g2+g2/2.0d+00+g3+g3/2.0d+00))/12.0d+00))*rcf
!
!     calculate four node capacatence matrix terms
!
      r(1,1) = e/6.0d+00
      r(1,2) = e/12.0d+00-c1*c2*e/60.0d+00
      r(1,3) = e/12.0d+00-c1*c3*e/60.0d+00
      r(1,4) = e*c1/60.0d+00
      r(2,1) = r(1,2)
      r(2,2) = e/6.0d+00-c2*c1*e/15.0d+00+(c2**2)*(c1**2)*e/90.0d+00
      r(2,3) = e/12.0d+00-c3*c1*e/30.0d+00-c2*c1*e/30.0d+00+c2*c3*
     $(c1**2)*e/90.0d+00
      r(2,4) = e*c1/30.0d+00-c2*(c1**2)*e/90.0d+00
      r(3,1) = r(1,3)
      r(3,2) = r(2,3)
      r(3,3) = e/6.0d+00-c3*c1*e/15.0d+00+(c3**2)*(c1**2)*e/90.0d+00
      r(3,4) = e*c1/30.0d+00-c3*(c1**2)*e/90.0d+00
      r(4,1) = r(1,4)
      r(4,2) = r(2,4)
      r(4,3) = r(3,4)
      r(4,4) = e*(c1**2)/90.0d+00

      fv(1)=1.0d+00/3.0d+00
      fv(2)=1.0d+00/3.0d+00*(1.0d+00-c2*c1/4.0d+00)
      fv(3)=1.0d+00/3.0d+00*(1.0d+00-c3*c1/4.0d+00)
      fv(4)=c1/12.0d+00

      do 200 n=1,4
        ddfn(n)=0.0d+00
        do 250 nn=1,4
          ddfn(n) = ddfn(n) + ((1.0d+00-theta)*s(n,nn) -
     $            r(n,nn))*temp(node(nn))
          s(n,nn) = theta*s(n,nn)+r(n,nn)
  250 continue
          be(n) = f*fv(n)*(phi(node(n))-phiold(node(n)))/(1.0d+00-
     $          phi(node(n)))-ddfn(n)
  200 continue

  300 if(eltp(m).ne.1)goto 700
  301 continue
!
!    calculate three node a matrix
!
      node(1)=ni(m)
      node(2)=nj(m)
      node(3)=nk(m)

      do 400 n=1,3
        dd=0.0d+00
        do 500 nn=1,3
          if(n.eq.nn)then
            r(n,nn)=e/3.0d+00
          else if (n.ne.nn)then
            r(n,nn)= 0.0 ! e/12.0d+00
          end if

          s(n,nn)=(a/4.0d+00)*(ex*b(nn,m)*b(n,m)+ez*c(nn,m)*c(n,m)+
     $            exz*(c(n,m)*b(nn,m)+b(n,m)*c(nn,m)))+
     $            ((rhofe*cvf(qq))/6.0d+00)*(qx(m)*b(nn,m)+
     $            qz(m)*c(nn,m))*fix


          dd = dd + ((1.0d+00-theta)*s(n,nn)-r(n,nn))*temp(node(nn))

          s(n,nn) = theta*s(n,nn)+r(n,nn)
  500   continue

          be(n)=f/3.0d+00*(phi(node(n))-phiold(node(n)))/(1.0d+00-
     $          phi(node(n)))-dd

  400 continue

      if(eltp(m).eq.1) then
      if(iout.ne.2.or.itt.ne.ntime) go to 401
      write(8,12) m
   12 format (' local fstiff for elem:',i4,//)
      do 402 lp=1,3
  402 write (8,21) (s(lp,mp),mp=1,3)
      write (8,13) (be(lp),lp=1,3)
   13 format (/,3x,'b vector:',3(1pe12.3))
  401 continue
      else if (eltp(m).ne.1) then
      if(iout.ne.2.or.itt.ne.ntime) go to 1401
      write(8,12) m
      do 1402 lp=1,4
 1402 write (8,21) (s(lp,mp),mp=1,4)
      write (8,13) (be(lp),lp=1,4)
 1401 continue
      endif

      sw=3
!
!     form global a matrix and b vector
!

  700   do 800 ii = 1,sw
          go to (35,40,45,50),ii
   35     ne = i
          go to 55
   40     ne = j
          go to 55
   45     ne = k
          go to 55
   50     ne = l
          go to 55
c
   55     do 900 jj=1,sw
            go to (60,65,70,75),jj
   60       me = i
            go to 80
   65       me = j
            go to 80
   70       me = k
            go to 80
   75       me = l
   80       kl=me-ne+ij2
          aa(ne,kl)=aa(ne,kl)+s(ii,jj)
  900     continue
          bb(ne)=bb(ne)+be(ii)
  800   continue
  100 continue
!
!     print out stiffness matrix
!
      if (iout.ne.2.or.itt.ne.ntime) go to 903
      write(8,11)
      do 902 ll=1,35
      write (8,21) (aa(ll,mm),mm=1,nhband)
  902 continue
      write (8,31) (bb(ll),ll=1,35)
  903 continue

   11 format (/1x,'energy stiffness matrix:')
   21 format (11(1pe12.3))
   31 format (//,1x,'b vector',10(1pe12.3))
      return
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      subroutine velo (x,z,qx,qz,ni,nj,nk,nl,rvis,rhof,eltp,
     $ nelem,mat,area,head,al,b,c,phi,perm,anisop,icoup
     $ ,delta,omegasum,node_angle,qxo,qzo,rhoo,heado,ioil,pload1
     $ ,ifltag,rotflag,fc_pp,it,time,d_time,fc_time,fc_crit)

!
!     this subroutine computes the darcy velocites
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 fc_pp(maxnodes)
      real*8 a,al(6,maxelems),b(6,maxelems),c(6,maxelems)
      real*8 x(maxnodes),z(maxnodes),head(maxnodes),area(maxelems)
      real*8 kxe,kze,kxze,perm(maxnodes),anisop(0:maxmats),phi(maxnodes)
      real*8 c1,c2,c3,psi2,psi3,nb(8),ng(8),b1,b2,b3,g1,g2,g3
      real*8 centx,centz,delta,omegasum,kx_temp,kz_temp,rho0
      real*8 e_angle,node_angle(maxnodes),qxo(maxelems),qzo(maxelems)
      real*8 heado(maxnodes),kxoe,kzoe,kxzoe,visw,viso,rhooe
      real*8 dcosd, dsind,rhof(maxnodes),rvis(maxnodes),rhoo(maxnodes)
      real*8 pload1(maxnodes),qx(maxelems),qz(maxelems),rhofe


      integer nelem,node(8),eltp(maxelems),ifltag,icoup
      integer ni(maxelems),nj(maxelems),nk(maxelems),nl(maxelems),qq
      integer i1,i2,i3,mat(maxelems),ioil,rotflag

      real*8   fc_time(maxelems),time,kmin,kmax,d_time,decay
      real*8  fc_crit
      integer it
      integer  fc_flag(maxelems)


      rho0 = 998.2d+00
      visw = 1.d-03
      viso = 5.d-02

      do 100 m=1,nelem

        if(eltp(m).eq.6)goto 100

        i=ni(m)
        j=nj(m)
        k=nk(m)
        l=nl(m)
        qq=mat(m)
        ae=area(m)
        if(i.eq.0) go to 100
        node(1)=i
        node(2)=j
        node(3)=k

        rhofe = (rhof(i) + rhof(j) + rhof(k))/3.d+00
        rrho = (rhofe - rho0)/rho0
        rvise = (rvis(i) + rvis(j) + rvis(k))/3.d+00
        if(icoup.eq.0) rvise = 1.0d+00
        if(icoup.eq.0) rrho = 0.0d+00
        a = rvise/(2.0d+00*ae)
        e_angle=(node_angle(i)+node_angle(j)+node_angle(k))/3.0
        if(qq.eq.ifltag) then
        kz_temp = (perm(i) + perm(j) + perm(k))/3.d+00
        kx_temp=kz_temp/anisop(qq)
        e_angle=delta-omegasum
        else
        kx_temp = (perm(i) + perm(j) + perm(k))/3.d+00
        kz_temp = kx_temp/anisop(qq)
        endif

          kxe  =kx_temp * (dcosd(e_angle)**2.0) +
     $         kz_temp * (dsind(e_angle)**2.0)
          kze  =kx_temp * (dsind(e_angle)**2.0)+
     $         kz_temp * (dcosd(e_angle)**2.0)
          kxze =(kx_temp-kz_temp)*dcosd(e_angle)*
     $                           dsind(e_angle)

      kmax = kxe*100
      kmin = kxe

      if(fc_pp(i).ge.fc_crit.or.fc_pp(j).ge.fc_crit
     $      .or.fc_pp(k).ge.fc_crit) then

      decay = dble(it)-fc_time(m)
      kxe = (kmax-kmin)*exp(-decay/d_time)+kmin
      kze = kxe/anisop(qq)
      kxze = 0.0
      end if

      if(fc_flag(m).gt.0) then
      if(fc_pp(i).lt.fc_crit.or.fc_pp(j).lt.fc_crit
     $      .or.fc_pp(k).lt.fc_crit) then

      decay = dble(it)-fc_time(m)
      if(decay.lt.0) decay = d_time*100

      kxe = (kmax-kmin)*exp(-decay/d_time)+kmin
      kze = kxe/anisop(qq)
      kxze = 0.0
      end if
      else if (fc_flag(m).eq.0) then
      kxe = kmin
      kze = kmin/anisop(qq)
      kxze = 0.0
      end if

         if(ioil.eq.1) then
         rhooe = (rhoo(i) + rhoo(j) + rhoo(k))/3.d+00
         rrhoo = (rhooe - rhofe)/rhooe
         kxzoe = 0.0d+00
         kxoe = kxe*rhooe*visw/(rhofe*viso)
         kzoe = kxoe/anisop(qq)
         endif

      e_angle=(node_angle(i)+node_angle(j)+node_angle(k))/3.0

        qx(m)=0.0d+00
        qz(m)=0.0d+00
        qxo(m)=0.0d+00
        qzo(m)=0.0d+00

      do 400 n=1,3
          qx(m)=qx(m) -
     $          a*(kxe*b(n,m)*head(node(n))+kxze*c(n,m)*head(node(n)))
          qz(m)=qz(m) -
     $          a*(kxze*b(n,m)*head(node(n))+kze*c(n,m)*head(node(n)))
      if(ioil.eq.1) then
          qxo(m) =  qxo(m) -
     $      a*(kxoe*b(n,m)*heado(node(n))+kxzoe*c(n,m)*heado(node(n)))
          qzo(m) =  qzo(m) -
     $      a*(kxzoe*b(n,m)*heado(node(n))+kzoe*c(n,m)*heado(node(n)))
       endif

  400  continue


  500  qx(m)=qx(m) - kxze*rrho*rvise
       qz(m)=qz(m) - kze*rrho*rvise

  100 continue
      return
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       subroutine heat_flow (hflux_x,hflux_z,ni,nj,nk,nl,eltp,
     $ nelem,mat,area,tld1,b,c,phi,tkf,tks,faultx
     & ,matp)
!
! this subroutine calculates surface heat flow
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!


      real*8 a,b(6,maxelems),c(6,maxelems),area(maxelems),phi(maxnodes)
      real*8 tke,tkf(0:maxmats),tks(0:maxmats),hflux_x(maxelems)
      real*8 hflux_z(maxelems),tld1(maxnodes),units_tk,phie

      integer nelem,node(8),eltp(maxelems),ni(maxelems),nj(maxelems)
      integer mat(maxelems),matp(maxnodes),nk(maxelems),nl(maxelems),qq

          units_tk = 3.1558e+07/4.187

      do 100 m=1,nelem
        if(eltp(m).eq.6)goto 100
        i=ni(m)
        j=nj(m)
        k=nk(m)
        qq=mat(m)
        ae=area(m)
        if(i.eq.0) go to 100

        node(1)=i
        node(2)=j
        node(3)=k

        phie = (phi(i)+phi(j)+phi(k))/3.0d+00

        tke = (tkf(qq)/units_tk*phie)+(tks(qq)/units_tk*(1.0d+00-phie))

        a = tke/(2.0d+00*ae)

        hflux_x(m)=0.0d+00
        hflux_z(m)=0.0d+00

  300 if(eltp(m).eq.6)goto 500


      do 400 n=1,3
          hflux_x(m)=hflux_x(m) -
     $          a*(b(n,m)*tld1(node(n)))
          hflux_z(m)=hflux_z(m) -
     $          a*(c(n,m)*tld1(node(n)))

  400  continue

       hflux_x(m)=hflux_x(m)*1000.0
       hflux_z(m)=hflux_z(m)*1000.0
  500  continue

 100    continue

      return
      end



!***********************************************************************
      subroutine zero (dt,theta,time,tja,tjb,tta,ttb,tha,thb,sgrad,
     $ alphas,omega,delta,domega,omegasum,total_time,msl,dtp,dtpr,
     $ gamma,delzmn,tbrstr,vztim1,vztim2,vztim3,vztim4,vztim5,vztim6,
     $ vztim7,vztim8,vztim9,vztim10,vztim11,vztim12,ibrine,iheat,ihchk,
     $ itchk,ntb,njb,nfband,nhband,nnode,nelem,nmat,icheck,icond,icoup,
     $ ncols,nbnn,nr,iout,it,nflt,iflow,tncols,ipr,tncf,nnr,rotflag,
     $ ijchk,dummy,nhb,ioil,iperm,imtag,ptec,pexp,psfc,pelem,iprfst,
     $ ntimp,nump,numwel,maxit,ifltag,p_flag,xsp,zsp,x,phi,phiold,z,
     $ node_angle,temp,told,ttic,tti,tempc,vs,rvis,bb,head,conc,wrrn,
     $ rhof,vxp,vzp,xp,zp,concp,pload1,pload2,vsm,presl,prslold,
     $ pres,sigmax,phimin,sigmae,ss,perm,sfcel,delsfc,presh,heado,
     $ headg,rhoo,rhog,zmax,voxp,vozp,matp,ipflag,icbflg,ielneb,elhom,
     $ numpnb,perm1,perm2,anisop,phi_o,beta,beta_ul,phi_ir,dif,alpha,
     $ tkf,tks,rhos,cvf,cvs,ldis,tdis,perct1,perct2,perct3,toc,dgrn,
     $ hflux_x,hflux_z,twrre,wrre,etemp,area,kx,kz,qx,qz,phitot,tr,rv,
     $ masoil,masgas,mash2o,masco2,oilvol,gasvol,persat,kgoil,kggas,
     $ fracoil,trold,qxo,qzo,mat,ni,nj,nk,eltp,nl,bwidth,centx,fbwidth,
     $ delcx,centz,xfault,faultx,faultz,nftype,nrel,wellid,aa,c,b,al,
     $ isourc,njn,nj1,nj2,bec,ber,nh,nt,lbnd,rbnd,nrows,tnrows,
     $ nrold,nslp,nc,nevap,nclay,nbn,icase,j1,j2,j1base,j2base,
     $ j1inc,j2inc,tpbase, hdbase,hdinc,tpinc,cnbase,clay,
     $ vzt,zcheck,delz,cevap,zmin,vztinc,grdfac,ver,vztbase)

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

! variables
!
      real*8 dt, theta, time, tja, tjb, tta,ttb
      real*8 tha,thb,sgrad, alphas, omega, delta
      real*8 domega, omegasum, total_time, msl
      real*8 dtp,dtpr,gamma,delzmn,tbrstr
      real*8 vztim1,vztim2,vztim3,vztim4,vztim5,vztim6
      real*8 vztim7,vztim8,vztim9,vztim10,vztim11,vztim12
!
      integer ibrine, iheat, ihchk, itchk, ntb, njb
      integer nfband, nhband, nnode, nelem, nmat
      integer icheck, icond, icoup, ncols
      integer nbnn, nr, iout, it, nflt, iflow, tncols,ipr
      integer tncf, nnr, rotflag, ijchk
      integer dummy,nhb,ioil,iperm,imtag
      integer ptec,pexp,psfc,pelem,iprfst
      integer ntimp,nump,numwel,maxit,ifltag,p_flag
!
! nodewise arrays
!
      real*8 xsp(maxnodes),zsp(maxnodes),x(maxnodes)
      real*8 phi(maxnodes),phiold(maxnodes),z(maxnodes)
      real*8 node_angle(maxnodes),temp(maxnodes),told(maxnodes)
      real*8 ttic(maxnodes),tti(maxnodes),tempc(maxnodes)
      real*8 vs(maxnodes),rvis(maxnodes),bb(maxnodes),head(maxnodes)
      real*8 conc(maxnodes),wrrn(maxnodes)
      real*8 rhof(maxnodes),vxp(maxnodes),vzp(maxnodes)
      real*8 xp(maxnodes),zp(maxnodes),concp(maxnodes)
      real*8 pload1(maxnodes),pload2(maxnodes)
      real*8 vsm(maxnodes),presl(maxnodes),prslold(maxnodes)
      real*8 pres(maxnodes),sigmax(maxnodes)
      real*8 phimin(maxnodes),sigmae(maxnodes)
      real*8 ss(maxnodes),perm(maxnodes),sfcel(maxnodes)
      real*8 delsfc(maxnodes),presh(maxnodes)
      real*8 heado(maxnodes),headg(maxnodes)
      real*8 rhoo(maxnodes),rhog(maxnodes)
      real*8 zmax(maxnodes),voxp(maxnodes),vozp(maxnodes)

      integer matp(maxnodes),ipflag(maxnodes),icbflg(maxnodes)
      integer ielneb(maxnodes,12),elhom(maxnodes),numpnb(maxnodes)
!
! fluid/rock property arrays
!
      real*8 perm1(0:maxmats),perm2(0:maxmats),anisop(0:maxmats)
      real*8 phi_o(0:maxmats),beta(0:maxmats)
      real*8 beta_ul(0:maxmats),phi_ir(0:maxmats)
      real*8 dif(0:maxmats),alpha(0:maxmats)
      real*8 tkf(0:maxmats),tks(0:maxmats),rhos(0:maxmats)
      real*8 cvf(0:maxmats),cvs(0:maxmats)
      real*8 ldis(0:maxmats),tdis(0:maxmats)
      real*8 perct1(0:maxmats),perct2(0:maxmats),perct3(0:maxmats)
      real*8 toc(0:maxmats),dgrn(0:maxmats)

      integer isourc(0:maxmats)
!
! boundary arrays
!
      integer njn(maxbnds),nj1(maxbnds),nj2(maxbnds),bec(maxbnds)
      integer ber(maxbnds),nh(maxbnds),nt(maxbnds),lbnd(maxbnds)
      integer rbnd(maxbnds),nrows(maxbnds),tnrows(maxbnds)
      integer nrold(maxbnds)
      integer nslp(maxbnds),nc(maxbnds),nevap(maxbnds)
      integer nclay(maxbnds),nbn(maxbnds),icase(maxbnds)

      real*8 j1(maxbnds),j2(maxbnds),j1base(maxbnds),j2base(maxbnds)
      real*8 j1inc(maxbnds),j2inc(maxbnds),tpbase(maxbnds)
      real*8 hdbase(maxbnds),hdinc(maxbnds),tpinc(maxbnds)
      real*8 cnbase(maxbnds),clay(maxbnds)
      real*8 vzt(maxbnds),zcheck(maxbnds),delz(maxbnds)
      real*8 cevap(maxbnds)
      real*8 zmin(maxbnds),vztinc(maxbnds)
      real*8 grdfac(maxbnds),ver(maxbnds),vztbase(maxbnds)

!
! elementwise arrays
!
      real*8 hflux_x(maxelems),hflux_z(maxelems),wrre(maxelems)
      real*8 twrre(maxelems),etemp(maxelems)
      real*8 area(maxelems),kx(maxelems),kz(maxelems)
      real*8 qx(maxelems),qz(maxelems),phitot(maxelems)
      real*8 tr(maxelems),rv(maxelems),masoil(maxelems)
      real*8 masgas(maxelems),mash2o(maxelems),masco2(maxelems)
      real*8 oilvol(maxelems),gasvol(maxelems),persat(maxelems)
      real*8 kgoil(maxelems), kggas(maxelems),fracoil(maxelems)
      real*8 trold(maxelems),qxo(maxelems),qzo(maxelems)

      integer mat(maxelems),ni(maxelems),nj(maxelems),nk(maxelems)
      integer eltp(maxelems),nl(maxelems)
!
! fault arrays
!
      real*8 bwidth(maxflts),centx(maxflts),fbwidth(maxflts)
      real*8 delcx(maxflts),centz(maxflts),xfault(maxflts)
      real*8 faultx(maxflts),faultz(maxflts)

      integer nftype(maxflts),nrel(maxflts),wellid(maxflts)

!
! multi-dimensional arrays
!
      real*8 aa(maxnodes,maxwdt),c(6,maxelems),b(6,maxelems)
      real*8 al(6,maxelems)
!
! zero constants
!
! real*8
      dt=0.0d+0
      theta=0.0d+0
      time=0.0d+0
      tja=0.0d+0
      tjb=0.0d+0
      tta=0.0d+0
      ttb=0.0d+0
      tha=0.0d+0
      thb=0.0d+0
      sgrad=0.0d+0
      alphas=0.0d+0
      omega=0.0d+0
      delta=0.0d+0
      domega=0.0d+0
      omegasum=0.0d+0
      total_time=0.0d+0
      msl=0.0d+0
      dtp=0.0d+0
      dtpr=0.0d+0
      gamma=0.0d+0
      delzmn=0.0d+0
      tbrstr=0.0d+0
      vztim1=1.0d+20
      vztim2=1.0d+20
      vztim3=1.0d+20
      vztim4=1.0d+20
      vztim5=1.0d+20
      vztim6=1.0d+20
      vztim7=1.0d+20
      vztim8=1.0d+20
      vztim9=1.0d+20
      vztim10=1.0d+20
      vztim11=1.0d+20
      vztim12=1.0d+20
! integers
      ibrine=0
      iheat=0
      ihchk=0
      itchk=0
      ntb=0
      njb=0
      nfband=0
      nhband=0
      nnode=0
      nelem=0
      nmat=0
      icheck=0
      icond=0
      icoup=0
      ncols=0
      nbnn=0
      nr=0
      iout=0
      it=0
      nflt=0
      iflow=0
      tncols=0
      ipr=0
      tncf=0
      nnr=0
      rotflag=0
      ijchk=0
      dummy=0
      nhb=0
      ioil=0
      iperm=0
      imtag=0
      ptec=0
      pexp=0
      psfc=0
      pelem=0
      iprfst=0
      ntimp=0
      nump=0
      numwel=0
      maxit=0
      ifltag=0
      p_flag=0

!
! zero nodewise variables
!

      do 100 l=1,maxnodes
      do 200 m=1,maxwdt
      aa(l,m)=0.0d+00
  200 continue
      bb(l)=0.0d+00
      xsp(l)=0.0d+00
      zsp(l)=0.0d+00
      x(l)=0.0d+00
      phi(l)=0.0d+00
      phiold(l)=0.0d+00
      z(l)=0.0d+00
      node_angle(l)=0.0d+00
      temp(l)=0.0d+00
      told(l)=0.0d+00
      ttic(l)=0.0d+00
      tti(l)=0.0d+00
      tempc(l)=0.0d+00
      vs(l)=0.0d+00
      rvis(l)=0.0d+00
      bb(l)=0.0d+00
      head(l)=0.0d+00
      conc(l)=0.0d+00
      wrrn(l)=0.0d+00
      rhof(l)=0.0d+00
      vxp(l)=0.0d+00
      vzp(l)=0.0d+00
      xp(l)=0.0d+00
      zp(l)=0.0d+00
      concp(l)=0.0d+00
      pload1(l)=0.0d+00
      pload2(l)=0.0d+00
      vsm(l)=0.0d+00
      presl(l)=0.0d+00
      prslold(l)=0.0d+00
      pres(l)=0.0d+00
      sigmax(l)=0.0d+00
      phimin(l)=0.0d+00
      sigmae(l)=0.0d+00
      ss(l)=0.0d+00
      perm(l)=0.0d+00
      sfcel(l)=0.0d+00
      delsfc(l)=0.0d+00
      presh(l)=0.0d+00
      heado(l)=0.0d+00
      headg(l)=0.0d+00
      rhoo(l)=0.0d+00
      rhog(l)=0.0d+00
      zmax(l)=0.0d+00
      voxp(l)=0.0d+00
      vozp(l)=0.0d+00
      tti(l) = 0.0d+00
      ttic(l) = 0.0d+00
!
! integer variables
!
      do 201 m=1,12
      ielneb(l,m)=0
  201 continue
      matp(l)=0
      ipflag(l)=0
      icbflg(l)=0
      elhom(l)=0
      numpnb(l)=0
  100 continue
!
! zero elementwise variables
!
      do 300 m=1,nelem
      do 301 l=1,6
      c(l,m)=0.0d+00
      b(l,m)=0.0d+00
      al(l,m)=0.0d+00
  301 continue
! real*8
      qx(m)=0.0d+00
      qz(m)=0.0d+00
      hflux_x(m)=0.0d+00
      hflux_z(m)=0.0d+00
      wrre(m)=0.0d+00
      twrre(m)=0.0d+00
      etemp(m)=0.0d+00
      area(m)=0.0d+00
      kx(m)=0.0d+00
      kz(m)=0.0d+00
      qx(m)=0.0d+00
      qz(m)=0.0d+00
      phitot(m)=0.0d+00
      tr(m)=0.0d+00
      rv(m)=0.0d+00
      masoil(m)=0.0d+00
      masgas(m)=0.0d+00
      mash2o(m)=0.0d+00
      masco2(m)=0.0d+00
      oilvol(m)=0.0d+00
      gasvol(m)=0.0d+00
      persat(m)=0.0d+00
      kgoil(m)=0.0d+00
      kggas(m)=0.0d+00
      fracoil(m)=0.0d+00
      trold(m)=0.0d+00
      qxo(m)=0.0d+00
      qzo(m)=0.0d+00
! integer
      mat(m)=0
      ni(m)=0
      nj(m)=0
      nk(m)=0
      eltp(m)=0
      nl(m)=0
  300 continue

      do 400 n=0,maxmats
! real*8
      perm1(n)=0.0d+00
      perm2(n)=0.0d+00
      anisop(n)=0.0d+00
      phi_o(n)=0.0d+00
      beta(n)=0.0d+00
      beta_ul(n)=0.0d+00
      phi_ir(n)=0.0d+00
      dif(n)=0.0d+00
      alpha(n)=0.0d+00
      tkf(n)=0.0d+00
      tks(n)=0.0d+00
      rhos(n)=0.0d+00
      cvf(n)=0.0d+00
      cvs(n)=0.0d+00
      ldis(n)=0.0d+00
      tdis(n)=0.0d+00
      perct1(n)=0.0d+00
      perct2(n)=0.0d+00
      perct3(n)=0.0d+00
      toc(n)=0.0d+00
      dgrn(n)=0.0d+00
! integer
      isourc(n)=0
  400 continue
      do 500 n=1,maxflts
! real*8
      bwidth(n)=0.0d+00
      centx(n)=0.0d+00
      fbwidth(n)=0.0d+00
      delcx(n)=0.0d+00
      centz(n)=0.0d+00
      xfault(n)=0.0d+00
      faultx(n)=0.0d+00
      faultz(n)=0.0d+00
! integer
      nftype(n)=0
      nrel(n)=0
      wellid(n)=0
  500 continue
      do 600 n=1,maxbnds
      njn(n)=0
      nj1(n)=0
      nj2(n)=0
      bec(n)=0
      ber(n)=0
      nh(n)=0
      nt(n)=0
      lbnd(n)=0
      rbnd(n)=0
      nrows(n)=0
      tnrows(n)=0
      nrold(n)=0
      nslp(n)=0
      nc(n)=0
      nevap(n)=0
      nclay(n)=0
      nbn(n)=0
      icase(n)=0
      j1(n)=0.0d+0
      j2(n)=0.0d+0
      j1base(n)=0.0d+0
      j2base(n)=0.0d+0
      j1inc(n)=0.0d+0
      j2inc(n)=0.0d+0
      tpbase(n)=0.0d+0
      hdbase(n)=0.0d+0
      hdinc(n)=0.0d+0
      tpinc(n)=0.0d+0
      cnbase(n)=0.0d+0
      clay(n)=0.0d+0
      vzt(n)=0.0d+0
      zcheck(n)=0.0d+0
      delz(n)=0.0d+0
      cevap(n)=0.0d+0
      zmin(n)=0.0d+0
      vztinc(n)=0.0d+0
      grdfac(n)=0.0d+0
      ver(n)=0.0d+0
      vztbase(n)=0.0d+0
  600 continue
      return
      end
!**********************************************************************
!  gaussian elimination for a symmetric,positive definite band matrix
!
!
      subroutine fgauss(aa,bb,nnode,nfband,head,nh,nhb,iout)
!
!
! *********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      dimension aa(maxnodes,maxwdt),bb(maxnodes),head(maxnodes)
      integer nh(maxbnds),nhb,nnode,nfband,iout

      if (iout.ne.2) go to 902
      write(8,11)
      do 901 ll=1,40
      write (8,21) (aa(ll,mm),mm=1,20)
  901 continue
      write (8,31) (bb(ll),ll=1,40)
  902 continue


      do 30 n=1,nnode
      do 20 l=2,nfband
      if(aa(n,l).eq.0.0)go to 20
      c=aa(n,l)/aa(n,1)
      i=n+l-1
      if(i.gt.nnode)go to 20
      j=0
      do 10 k=l,nfband
      j=j+1
   10 if(aa(n,k).ne.0.0)aa(i,j)=aa(i,j)-c*aa(n,k)
      aa(n,l)=c
      bb(i)=bb(i)-aa(n,l)*bb(n)
   20 continue
      bb(n)=bb(n)/aa(n,1)
   30 continue

      n=nnode
   40 do 50 k=2,nfband
      l=n+k-1
      if(l.gt.nnode)go to 60
      if(aa(n,k).ne.0.0)bb(n)=bb(n)-aa(n,k)*bb(l)
   50 continue
   60 n=n-1
      if(n.gt.0)go to 40

      do 700 i=1,nnode
      do 600 n=1,nhb
      if(i.eq.nh(n)) go to 700
  600 continue
      head(i)=bb(i)
  700 continue
   11 format (/1x,'sub. fgaus; stiffness matrix:')
   21 format (11(1pe12.3))
   31 format (//,1x,'b vector',10(1pe12.3))
      return
      end
! **********************************************************************
! subroutine hgauss performs gaussian elimination for non-symmetric band
! matrices(diagonally dominant)
!

      subroutine hgauss(aa,bb,nnode,nhband,temp,nt,ntb)
!
!
! *********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      integer ntb,nt(maxbnds),nhband,nnode

      dimension aa(maxnodes,maxwdt),bb(maxnodes),temp(maxnodes)

      ij=(nhband-3)/2
      ij2=ij+2
      ij3=ij+3

      do 30 n=1,nnode
      if(aa(n,ij2).gt.1.0d-7)goto 5
      write(8,7) n
      write(*,7) n
      write(*,*)'stop6',temp(n)
      stop
   05 bb(n)=bb(n)/aa(n,ij2)
      if(n.ge.nnode)go to 30
      k=0
   10 do 25 j=ij3,nhband
      k=k+1
      aa(n,j)=aa(n,j)/aa(n,ij2)
      kk=0
      do 20 i=ij3,nhband
      kk=kk+1
      l=n+i-ij2
      if(l.gt.nnode)go to 20
      m=ij2+k-kk
      mm=ij2-kk
      aa(l,m)=aa(l,m)-aa(l,mm)*aa(n,j)
   20 continue
      jj=ij2-k
      ll=n+j-ij2
      if(ll.gt.nnode)go to 30
      bb(ll)=bb(ll)-aa(ll,jj)*bb(n)
   25 continue
   30 continue
      n=nnode
   40 n=n-1
      if(n.le.0) go to 400
      do 50 j=ij3,nhband
        l=n+j-ij2
        if(l.gt.nnode)go to 50
        bb(n)=bb(n)-aa(n,j)*bb(l)
   50 continue
      go to 40

  400 continue
      do 700 i=1,nnode
      do 600 n=1,ntb
      if(i.eq.nt(n)) go to 700
  600 continue
      temp(i)=bb(i)
  700 continue
   07 format(4x,'subroutine hgauss ',/,'divided by zero at node',i4)
      end
!*******************************************************************
!
      subroutine density(z,head,temp,rhof,nnode,conc)
!
!   this subroutine calculates fluid density based on temperature,
!   pressure, and concentration
!
!*******************************************************************
!
      implicit real*8 (a-h,o-z)

!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!
      real*8 z(maxnodes),rhof(maxnodes),temp(maxnodes),head(maxnodes)
      real*8 conc(maxnodes)
      real*8 am2,am1,a0,ap1,ap2,bm2,bm1,b0,bp1,bp2,c0,c1
      real*8 d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3,h1,h2,h3
      real*8 t,w,p,t2,a,b,c,d,e,f,g,h,sv

      am2= 1.006741d+02
      am1=-1.127522d+00
      a0 = 5.916365d-03
      ap1=-1.035794d-05
      ap2= 9.270048d-09

      bm2= 1.042948d+00
      bm1=-1.1933677d-2
      b0 = 5.307535d-05
      bp1=-1.0688768d-7
      bp2= 8.492739d-11

      c0= 1.23268d-09
      c1=-6.861928d-12
      d1=-2.5166d-03
      d2= 1.11766d-05
      d3=-1.70552d-08
      e1= 2.84851d-03
      e2=-1.54305d-05
      e3= 2.23982d-08
      f1=-1.5106d-05
      f2= 8.4605d-08
      f3=-1.2715d-10
      g1= 2.7676d-05
      g2=-1.5694d-07
      g3= 2.3102d-10
      h1= 6.4633d-08
      h2=-4.1671d-10
      h3= 6.8599d-13

      do 100 l=1,nnode
      t=temp(l)
      w=conc(l)
      if(t.lt.0) t=0.0
      if(w.lt.0.0d+00) w = 0.0d+00
      if(w.gt.0.3) w = 0.3d+00
      p=998.2d0*9.806d0*(head(l)-z(l))/1.0d+06
      if(p.lt.0.0d+0) p = 0.0d+0
      t=t + 273.15
      if(t.gt.473.15) t=473.15
      t2=t*t
      a=am2*(1.0d0/t**2)+am1*(1.0d0/t)+a0+ap1*t+ap2*t2
      b=bm2*(1.0d0/t**2)+bm1*(1.0d0/t)+b0+bp1*t+bp2*t2
      c=c0 + c1*t
      d=d1 + d2*t + d3*t2
      e=e1 + e2*t + e3*t2
      f=f1 + f2*t + f3*t2
      g=g1 + g2*t + g3*t2
      h=h1 + h2*t + h3*t2
      sv=a-b*p-c*p**2+w*d+w**2*e-w*f*p-w**2*g*p-0.5d0*h*p**2
      rhof(l)=1.0d0/sv
  100 continue
      return
      end
! **********************************************************************
!
      subroutine viscos (z,nnode,temp,head,rvis,conc)
!
! subroutine viscos calculates fluid viscosity (kg/m yr)
! as a function of pressure ,temperature,and salinity(moles/kg h20)
!
! **********************************************************************

      implicit real*8 (a-h,o-z)

!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      integer nnode
      real*8 muf,muop,muow20,ms,m0,m1,m2
      real*8 muo,lmuow,lzpv,mull

      real*8 z(maxnodes),head(maxnodes),rvis(maxnodes)
      real*8 temp(maxnodes),conc(maxnodes)

      muow20=1002.0d0

      muo=3.162d+04
      alpha1= 1.2378d0
      alpha2=-1.303d-03
      alpha3= 3.06d-06
      alpha4= 2.55d-08

      a1= 3.324d-02
      a2= 3.624d-03
      a3=-1.879d-04
      b1=-3.960d-02
      b2= 1.020d-02
      b3=-7.020d-04

      beta0=-1.297d+00
      beta1= 5.740d-02
      beta2=-6.970d-04
      beta3= 4.470d-06
      beta4=-1.050d-08

      gamma0=0.545d+00
      gamma1=2.800d-03

      m0=6.044d+00
      m1=2.80d-03
      m2=3.6d-05

      bata1=2.50d+00
      bata2=-2.00d+00
      bata3=0.5d+00

      do 100 l=1,nnode
      t=temp(l)
      s=conc(l)
      if(t.lt.0.0d+00) t = 0.0d+00
      if(t.gt.2.0d+02) t = 2.0d+02
      if(s.lt.0.0d+00) s = 0.0d+00
      if(s.gt.0.3d+00) s = 0.3d+00
      p=998.2d0*9.806d0*(head(l)-z(l))/1.0d+06
      if(p.lt.0.0d+0) p = 0.0d+0


      if(t.gt.20.0d0)go to 200
      yo=t + 133.15
      y1=248.37/yo
      y2=10.0d0**y1
      mull=2.40d-05 * y2 * 31.558d+06
      rvis(l) = muo/mull
      go to 105

  200 wnacl=s*1000.0d0/(1.0d0 - s)
      csalt=wnacl/58.443d0
      sum=alpha1*(20.0d0-t) + alpha2*((20.0d0-t)**2) +
     $alpha3*((20.0d0-t)**3) + alpha4*((20.0d0-t)**4)
      quot=sum/(96.0d0+t)
      lmuow=dlog10(muow20)+quot
      a=a1*csalt + a2*csalt**2 + a3*csalt**3
      b=b1*csalt + b2*csalt**2 + b3*csalt**3
      lzpv=lmuow + a +b*quot
      muop=10.0d0**lzpv
      ms=m0 + m1*t + m2*(t**2)
      bstar=bata1*(csalt/ms) + bata2*((csalt/ms)**2) +
     $      bata3*((csalt/ms)**3)
      bw=beta0+beta1*t+beta2*(t**2)+beta3*(t**3)+beta4*(t**4)
      bse=gamma0 + gamma1*t - bw
      beeta=bse*bstar + bw
      fp=p*1.00d-03
      muf=muop*(1.0d0 + beeta*fp)
      mull=muf*31.558
      rvis(l) = muo/mull

  105 continue
  100 continue
      return
      end
!***********************************************************************
      subroutine update(temp,told,nnode)
!
!   this subroutine updates the told variable for calculation
!   of thermal expansion in subroutine flow
!
!***********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      real*8 temp(maxnodes),told(maxnodes)

      integer nnode

      do 100 l=1,nnode
      told(l) = temp(l)
  100 continue
      return
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      subroutine mesh (nrows,xsp,zsp,x,z,ni,nj,nk,nl,bec,ber,
     $ nelem,nnode,nflt,eltp,mat,matp,nfband,nhband,ncols,tnrows,tncols,
     $ lbnd,rbnd,xfault,vs,phi,phiold,head,temp,told,tempc,tti,ttic,
     $ tncf,dt,it,nftype,delta,omega,faultx,faultz,
     $ node_angle,sfcel,delsfc,rotflag,conc,presl,prslold)
!
!     this subroutine develops the conectivity matrix for a
!     solution domaine containing a combination of three and four
!     node triangular elements.  the subroutine develops blocks of
!     three node elements then compaires the locations of nodes on
!     the adjoining sides of the faultblocks and creates four node
!     elements when appropriate.
!     the program also eliminates nodes occuping the same space
!     and redevelops the conectivity matrix as necessary.
!     the surfaces of the faultblocks can be erroded.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

      integer  count,ip,in(8,maxelems),nflt,nnode
      integer  sk,sk2,sw,coe,sw3,sof,a,b,c,roe,tval,ncf
      integer  ber(maxbnds),bec(maxbnds),nelem,ael,node,eltp(maxelems)
      integer  bae,bnode,p,k,h,elcount,tncf,baet,ncols
      integer  onnode,fth,nr(maxbnds,20),nfn,oor(maxbnds),nor(maxbnds)
      integer  lnode,last,fcel,anc,nrows(maxbnds)
      integer  ncel,cel,tcel,bwfw,s1,s2,ni(maxelems),nj(maxelems)
      integer  fcol,nk(maxelems),nl(maxelems)
      integer  tr1,tr2,tr3,tr6,tr7,tr8,tr9,tr10,tr11,tr12
      integer  mat(maxelems),matp(maxnodes),qq,maxnode,minnode,diff
      integer  nhband,tnrows(maxbnds),tncols,fnode,lbnd(maxbnds)
      integer  tdcount,bn,it,mc,lc,nftype(maxflts),tr50,fflag
      integer  rotflag,trip,maxdiff,nfband
      integer  tnc,fictel,rbnd(maxbnds)
c
      real*8  zsp(maxnodes),xsp(maxnodes),nlow,olow,zoor(maxbnds)
      real*8  x(maxnodes),z(maxnodes),xfault(maxflts),vs(maxnodes)
      real*8  phi(maxnodes),phild1(maxnodes),phiold(maxnodes),diflow
      real*8  head(maxnodes),hold(maxnodes),tld2(maxnodes)
      real*8  dt,temp(maxnodes),tld1(maxnodes),told(maxnodes)
      real*8  tti(maxnodes),ttiold(maxnodes),ttic(maxnodes)
      real*8  saf,tcold(maxnodes),phild2(maxnodes)
      real*8  delta,omega,faultx(maxflts),faultz(maxflts)
      real*8  node_angle(maxnodes),nd_agld(maxnodes),tempc(maxnodes)
      real*8  alpha4,sfcel(maxnodes),sfcelold(maxnodes),vsold(maxnodes)
      real*8  delsfc(maxnodes),delsfcold(maxnodes),xx,zz
      real*8  conc(maxnodes),concld(maxnodes),tdata(maxelems)
      real*8  presl(maxnodes),prslold(maxnodes),tticld(maxnodes)
      real*8  prsltmp1(maxnodes),prsltmp2(maxnodes)
      real*8  dtand

      real*4 axz,adtz,adif1


      nnode=nnode+tncf
!
!     this loop zeros out the previous conectivity matrix
!
      do 50 j=1,nelem
        do 50 k=1,4
          in(k,j)=0
   50 continue

      do 55 j=1,nnode
        z(j)=zsp(j)
        x(j)=xsp(j)
        vsold(j)=vs(j)
        phild1(j) = phi(j)
        phild2(j) = phiold(j)
        hold(j) = head(j)
        concld(j) = conc(j)
        tld1(j) = temp(j)
        tld2(j) = told(j)
        tcold(j) = tempc(j)
        ttiold(j) = tti(j)
        tticld(j) = ttic(j)
        nd_agld(j)=node_angle(j)
        sfcelold(j)=sfcel(j)
        delsfcold(j)=delsfc(j)
        prsltmp1(j)= presl(j)
        prsltmp2(j)= prslold(j)
   55 continue
c
      do 60 j=1,ncols
        tnrows(j)=nrows(j)
   60 continue
!
!     here we are developing the conectivity matricies for all
!     elements in a block wise pattern.
!
      count = 1
      l = 1
      sk = 1
      nelem = 0
      tnode=0
      nnode=0
!
!     converts from nrows(x) to nr(x,y)
!
      tncols=ncols
      do 75 k=1,nflt+1
        do 75 m=1,bec(k)+1
          nr(m,k)=nrows(count)
          count=count+1
   75   continue

      count=1
      do 100 k = 1, nflt + 1
        ip = tnode+1
        sk = 1
        sk2 =1
        fictel = 0
        tnr=ber(k)+1
        tnc=bec(k)+1
        bnode=0
        do 150 m=1,bec(k)+1
          bnode = bnode+nr(m,k)
  150   continue
        tnode=tnode+bnode
        do 200 m = 1,tnr-1
          sk = m
          sk2 = m
          do 250 p = 1,tnc-1
            if(mod(sk,2).gt.0.0) then
              s1 = 1
              s2 = 0
            else
              s1 = 0
              s2 = 1
             end if
            if((p.eq.1).or.(p.eq.tnc-1))then
              bn=0
            else
              bn=1
            end if
!
!     this section is developes left sided elements
!
            fictel = 0
            if(((nr(p,k)-2).gt.nr(p+bn,k)).and.(m.gt.nr(p+1,k)))then
              fictel=1
            end if
            if(nr(p,k).gt.m)then
              if((nr(p,k).eq.m+1).or.(m.eq.nr(p+1,k)))then
                if(nr(p,k).gt.nr(p+1,k))then
                  s1 = 0
                  s2 = 1
                end if
              end if
              if(fictel.eq.1)then
                eltp(l) = 6
                fictel = 0
              else
                if((nr(p,k)-nr(p+1,k).eq.2).and.(nr(p+1,k).eq.m-1))then
                  s1=-1
                end if
                in(1,l) = ip
                in(2,l) = ip + nr(p,k) + (1 * s1)
                in(3,l) = ip + 1
                eltp(l)=1
              end if
            else
              eltp(l)=6
            end if
!
!      this section is developes right sided elements
!
            sk2=sk2+1
            if(((nr(p+1,k).gt.m).and.(m.lt.(nr(p,k)+2))))then
              if((nr(p+1,k).eq.m+1).or.(nr(p,k).eq.m))then
                if((mod(sk2,2).gt.0.0).or.(nr(p,k).eq.m+1))then
                  if(p.eq.tnc)then
                    bn=0
                  else
                    bn=1
                  end if
                  if(nr(p+1,k).gt.nr(p+1-bn,k))then
                    s1 = 1
                    s2 = 0
                  end if
                end if
              end if
              in(1,l+1) = ip + (nr(p,k) * s2)
              in(2,l+1) = ip + nr(p,k) + (1 * s2)
              in(3,l+1) = ip + (nr(p,k) * s1) + 1
              eltp(l+1)=1
              if((nr(p+1,k)-nr(p,k).ge.2).and.(nr(p,k).eq.m-1))then
                in(1,l+1) = ip-1
              end if
            else
              eltp(l+1)=6
            end if
            ip = ip + nr(p,k)
            l = count + (p*2)
            sk = sk + 1
  250     continue
          count = count + (tnc-1) * 2
          ip = nnode+1 + m
  200   continue
        nelem = count-1
        nnode=tnode
  100 continue
!
!     now we determine which of the block boundary elements are
!     four node elements and redevelope the conectivity matrices
!     for these elements and flag them for the four node flow
!     and heat transfer subroutines.
!

      tncf=0

      fnode=1
      fcol=1
      tnode=0
      bae = 0
      onnode = nnode

      do 340 k=1,nflt
        xfault(k)=0
  340 continue

      do 350 k = 1,nflt
        tr1 = 0
        tr2 = 0
        tr6 = 0
        tncf=tncf+ncf
        ncf=0
        bnode = 0
        do 400 m=1,bec(k)+1
          bnode = bnode+nr(m,k)
  400   continue

        do 450 j = 1,2
          if(j.eq.1) then
            sof = 1
            roe = ber(k)
            coe = bec(k)
            tnode=tnode+bnode
            node = tnode + 1 - tncf
            ael = bae + bec(k) * 2
          else
            roe = ber(k + 1)
            coe = bec(k + 1)
            sof = 2
            node = tnode - nr(bec(k)+1,k) - tncf+1
            bae = bae + bec(k) * ber(k) * 2
            ael = bae + 1
          endif
          if(mod(coe,2).ne.0.0) then
            sw = 2
          else
            sw = 1
          endif
          tr7=2

          do 500 t = 1, roe
            if(sof.eq.2)then
              mc=0
              do 505 lc=1,2
                if(nftype(k).eq.1)then
                  mat(ael+mc)=2
                  mc=mc+1
                end if
  505         continue
            end if
            if(sof.eq.1)then
              mc=0
              do 506 lc=1,2
                if(nftype(k).eq.1)then
                  mat(ael-mc)=2
                  mc=mc+1
                end if
  506         continue
            end if
            tr7=tr7+1
            sw = sw + 1
            if(sof.eq.1) then
              if(mod(sw,2).eq.0.0) then
                a = 1
                b = 2
              else
                a = 2
                b = 3
              end if
            else
              if(mod(tr7,2).eq.0.0) then
                trip=3
              else
                trip=1
              end if
              a = 1
              b = 3
            end if
            if(sof.eq.1)then
              anc=1
            else
              anc=bec(k)+1
            end if
!
!     this section reassigns the conectivity matrix for the nodes of
!    elements which are located along the fault for elements in which
!     the lower of the two nodes located on the fault is = to the
!     node being looked at denoted by m
!   saf is a factor of saftey
!
            saf=0.09d+00
            do 550 m = node, node + nr(anc,k + 2 - sof)-1
              if(m.eq.node+nr(anc,k+2-sof)-1)then
                saf=0.00d+00
              end if
              if(sof.eq.2)then
                if(eltp(ael).ne.6)then
                  if((z(m).ge.z(in(a,ael))-saf).and.
     $            (z(m).le.z(in(b,ael))+saf))then
                  tr50=1
                    if(m.eq.node + ber(k + 2 - sof))then
                    end if
                    if(z(m).le.z(in(a,ael))+saf)then
                      if(tr1.eq.1)then
                        tr1=0
                      else
                        ncf=ncf+1
                        nnode=nnode-1
                      end if
                      do 560 h = 1,bae
                        if(in(4,h).eq.in(b,ael))then
                          in(4,h)=in(b,ael)-ncf
                        end if
  560                 continue
                      if(eltp(ael).ne.6)then
!
!     this section of the code preserves the node number and nodal
!     properties of eliminated nodes for rencorperation after flow
!     and heat subroutines
!
                        hold(m)=(hold(m)+hold(in(a,ael)))/2.0d+00
                        in(a,ael)=m
                        in(b,ael)=in(b,ael)-ncf
                        x(in(b,ael))=x(in(b,ael)+ncf)
                        z(in(b,ael))=z(in(b,ael)+ncf)
                        vsold(in(b,ael)) = vsold(in(b,ael)+ncf)
                        phild1(in(b,ael)) = phild1(in(b,ael)+ncf)
                        phild2(in(b,ael)) = phild2(in(b,ael)+ncf)
                        hold(in(b,ael)) = hold(in(b,ael)+ncf)
                        concld(in(b,ael)) = concld(in(b,ael)+ncf)
                        tld1(in(b,ael)) = tld1(in(b,ael)+ncf)
                        tld2(in(b,ael)) = tld2(in(b,ael)+ncf)
                        tcold(in(b,ael)) = tcold(in(b,ael)+ncf)
                        ttiold(in(b,ael)) = ttiold(in(b,ael)+ncf)
                        tticld(in(b,ael)) = tticld(in(b,ael)+ncf)
                        nd_agld(in(b,ael)) = nd_agld(in(b,ael)+ncf)
                        sfcelold(in(b,ael)) = sfcelold(in(b,ael)+ncf)
                      delsfcold(in(b,ael)) = delsfcold(in(b,ael)+ncf)
                        prsltmp1(in(b,ael))= prsltmp1(in(b,ael)+ncf)
                        prsltmp2(in(b,ael))= prsltmp2(in(b,ael)+ncf)

                      end if
                      if(eltp(ael+1).ne.6)then
                        if(tr2.ne.1)then
                          if(trip.eq.3)then
                            do 570 p=node,node+nr(anc,k+2-sof)-1
                         if(z(p).le.z(in(3,ael+1))+saf)then
                           tr8=1
                           tval=p
                          end if
  570                       continue
                            if(tr8.eq.1)then
                              in(3,ael+1) = tval
                            else
                              in(3,ael+1) = in(3,ael+1) - ncf
                            end if
                          else
                            in(1,ael+1) = m
                          end if
                        end if
                      end if
                      tr2=1
                    end if
!
!     this section reassigns the conectivity matrix for the nodes of
!     elements which are located along the fault for elements in which
!     the two nodes located on the fault are located on ither side of
!     of the node being looked at denoted by m
!
                    if((z(m).gt.z(in(a,ael))+saf).and.
     $              (z(m).lt.z(in(b,ael))-saf))then
                      if(tr2.ne.1)then
                        if(eltp(ael).ne.6)then
                          do 580 h = 1,bae
                            if(in(4,h).eq.in(a,ael))then
                              in(4,h)=in(4,h)-ncf
                            end if
  580                     continue

                          in(a,ael)=in(a,ael)-ncf
                          x(in(a,ael))=x(in(a,ael)+ncf)
                          z(in(a,ael))=z(in(a,ael)+ncf)
                          vsold(in(a,ael)) = vsold(in(a,ael)+ncf)
                          phild1(in(a,ael)) = phild1(in(a,ael)+ncf)
                          phild2(in(a,ael)) = phild2(in(a,ael)+ncf)
                          hold(in(a,ael)) = hold(in(a,ael)+ncf)
                          concld(in(a,ael)) = concld(in(a,ael)+ncf)
                          tld1(in(a,ael)) = tld1(in(a,ael)+ncf)
                          tld2(in(a,ael)) = tld2(in(a,ael)+ncf)
                          tcold(in(a,ael)) = tcold(in(a,ael)+ncf)
                          ttiold(in(a,ael)) = ttiold(in(a,ael)+ncf)
                          tticld(in(a,ael)) = tticld(in(a,ael)+ncf)
                          nd_agld(in(a,ael)) = nd_agld(in(a,ael)+ncf)
                          sfcelold(in(a,ael)) = sfcelold(in(a,ael)+ncf)
                        delsfcold(in(a,ael)) = delsfcold(in(a,ael)+ncf)
                       prsltmp1(in(a,ael))=   prsltmp1(in(a,ael)+ncf)
                       prsltmp2(in(a,ael))=   prsltmp2(in(a,ael)+ncf)

                          in(b,ael)=in(b,ael)-ncf
                          x(in(b,ael))=x(in(b,ael)+ncf)
                          z(in(b,ael))=z(in(b,ael)+ncf)
                          vsold(in(b,ael)) = vsold(in(b,ael)+ncf)
                          phild1(in(b,ael)) = phild1(in(b,ael)+ncf)
                          phild2(in(b,ael)) = phild2(in(b,ael)+ncf)
                          hold(in(b,ael)) = hold(in(b,ael)+ncf)
                          concld(in(b,ael)) = concld(in(b,ael)+ncf)
                          tld1(in(b,ael)) = tld1(in(b,ael)+ncf)
                          tld2(in(b,ael)) = tld2(in(b,ael)+ncf)
                          tcold(in(b,ael)) = tcold(in(b,ael)+ncf)
                          ttiold(in(b,ael)) = ttiold(in(b,ael)+ncf)
                          tticld(in(b,ael)) = tticld(in(b,ael)+ncf)
            nd_agld(in(b,ael)) = nd_agld(in(b,ael)+ncf)
            sfcelold(in(b,ael)) = sfcelold(in(b,ael)+ncf)
            delsfcold(in(b,ael)) = delsfcold(in(b,ael)+ncf)
                prsltmp1(in(b,ael))=   prsltmp1(in(b,ael)+ncf)
                prsltmp2(in(b,ael))=   prsltmp2(in(b,ael)+ncf)

                        end if
                        in(trip,ael+1) = in(trip,ael+1)-ncf
                        tr2=1
                      end if
                    end if
!
!     this section reassigns the conectivity matrix for the nodes of
!     elements which are located along the fault for elements in which
!     the upper of the two nodes located on the fault is = to the
!     node being looked at denoted by m
!
                    if(z(m).ge.z(in(b,ael))-saf)then
                      do 590 h = 1,bae
                        if(in(4,h).eq.in(a,ael))then
                          in(4,h)=in(a,ael)-ncf
                        end if
  590                   continue
                        if(tr2.ne.1)then
                          if(eltp(ael).ne.6)then
                            in(a,ael)=in(a,ael)-ncf
                            x(in(a,ael))=x(in(a,ael)+ncf)
                            z(in(a,ael))=z(in(a,ael)+ncf)
                            vsold(in(a,ael)) = vsold(in(a,ael)+ncf)
                            phild1(in(a,ael)) = phild1(in(a,ael)+ncf)
                            phild2(in(a,ael)) = phild2(in(a,ael)+ncf)
                            hold(in(a,ael)) = hold(in(a,ael)+ncf)
                            concld(in(a,ael)) = concld(in(a,ael)+ncf)
                            tld1(in(a,ael)) = tld1(in(a,ael)+ncf)
                            tld2(in(a,ael)) = tld2(in(a,ael)+ncf)
                            tcold(in(a,ael)) = tcold(in(a,ael)+ncf)
                            ttiold(in(a,ael)) = ttiold(in(a,ael)+ncf)
                            tticld(in(a,ael)) = tticld(in(a,ael)+ncf)
           nd_agld(in(a,ael)) = nd_agld(in(a,ael)+ncf)
	             sfcelold(in(a,ael)) = sfcelold(in(a,ael)+ncf)
		     delsfcold(in(a,ael)) = delsfcold(in(a,ael)+ncf)
                      prsltmp1(in(a,ael))=   prsltmp1(in(a,ael)+ncf)
                      prsltmp2(in(a,ael))=   prsltmp2(in(a,ael)+ncf)
                          end if
                        end if
                        ncf=ncf+1
                        nnode=nnode-1
                        if(eltp(ael).ne.6)then
                          hold(m)=(hold(m)+hold(in(b,ael)))/2.0d+00
                          in(b,ael)=m
                        end if
                        if(eltp(ael+1).ne.6)then
                          if(trip.eq.3)then
                            in(3,ael+1) = m
                          else
                          if(tr2.ne.1)then
                            in(1,ael+1) = in(1,ael+1) - ncf + 1
                          end if
                        end if
                      end if
                      tr1=1
                    end if
                  end if
                else
                  if((eltp(ael+1).ne.6).and.(tr6.ne.1))then
                    in(1,ael+1)=in(3,ael-bec(k+1)*2)
                    tr6=1
                  end if
                end if
              end if
!
!     this section determines the element type and adds in the
!     fourth node to some of the elements
!
              if(eltp(ael).eq.6)goto 550
              if(z(m).lt.z(in(a,ael))+saf)then
                goto 550
              else
                if(z(m).lt.z(in(b,ael))-saf)then
                  if(sof.eq.1)then
                    if(a.eq.2)then
                      eltp(ael) = 2
                    else
                      eltp(ael) = 3
                    end if
                  else
                    eltp(ael) = 4
                  end if
                  in(4,ael) = m
                endif
              end if
  550       continue
c
            if((sof.eq.2).and.(ncf.ne.0))then
            if((tr50.eq.0).and.(eltp(ael).ne.6))then
             in(a,ael)=in(a,ael)-ncf
             in(b,ael)=in(b,ael)-ncf
             in(trip,ael+1)=in(trip,ael+1)-ncf
             write(8,*)'warning 5 node element near',ael,
     $                  in(a,ael),ncf,
     $                  in(b,ael),in(2,ael),eltp(ael)
            end if
            end if
c
            ael = ael + coe * 2
            tr2=0
            tr50=0
  500     continue
  450   continue

!
!       this loop changes the node numbering of  elements in
!       the first elemental column in the faultblock located to
!       the right of the fault
!
        sw3=0
        tr3=1
        fth=4
        do 600 n = bae+1,ael+1-coe*2,coe*2
          do 620 p = 1,2
            do 640 h = 1,p
              if(mod(tr3,fth).eq.0.0)then
                sw3=1
              end if
              if(eltp(n+p-1).eq.6)goto640
              if(eltp(n+p-2).eq.6)sw3=0
              in(1+h-sw3,n+p-1)=in(1+h-sw3,n+p-1)-ncf
  640       continue
            tr3=tr3+1
            sw3=0
  620     continue
  600   continue
!
!      this loop changes the node numbering of the rest of the
!      elements in faultblock located to the right of the fault
!
        baet=bae+2
        coun=2
        do 680 a = 2,coe
          do 700 n = baet,ael+1-coe*2+coun,coe*2
            do 720 h = 1,2
              do 740 p = 1,4
                if(in(p,n+h).eq.0)goto 740
                if(eltp(n+h).eq.6)goto 740
                in(p,n+h)=in(p,n+h)-ncf
  740         continue
  720       continue
  700     continue
          baet=baet+2
          coun=coun+2
  680   continue
!
!     this section reassigns the conectivity matrix for elements
!     located in the faultblocks remaining to the right of the two
!     faultblocks currently being processed
!
        elcount=0
        do 760 n=1,k+1
          elcount=elcount+bec(n)*ber(n)*2
  760   continue
        do 780 n = elcount+1,nelem
          do 780 p = 1,4
            if(in(p,n).eq.0)goto 780
            if(in(p,n).lt.0)write(*,*)'lt 0'
            in(p,n)=in(p,n)-ncf
  780   continue
!
!     reassigns nodal properties to nodes located to the right of but
!     not located on the fault
!
        do 800 n = in(2,bae+2),nnode
          x(n)=x(n+ncf)
          z(n)=z(n+ncf)
          vsold(n) = vsold(n+ncf)
          phild1(n) = phild1(n+ncf)
          phild2(n) = phild2(n+ncf)
          hold(n) = hold(n+ncf)
          concld(n) = concld(n+ncf)
          tld1(n) = tld1(n+ncf)
          tld2(n) = tld2(n+ncf)
          tcold(n) = tcold(n+ncf)
          ttiold(n) = ttiold(n+ncf)
          nd_agld(n) = nd_agld(n+ncf)
          sfcelold(n) = sfcelold(n+ncf)
          delsfcold(n) = delsfcold(n+ncf)
          prsltmp1(n)=   prsltmp1(n+ncf)
          prsltmp2(n)=   prsltmp2(n+ncf)
  800   continue
c
c     zeros out the node numbers that are no longer used
c
        do 820 n = nnode+1,onnode
          x(n)=0.0
          z(n)=0.0
  820   continue
!
!     the remainder of the code sorts the node numbers out along the
!     fault columns and reassings the conectivity matrix and nodal
!     properties for nodes and elements located along the faults
!     this section creates a tmporary array of the node
!     numbers of nodes located allong the fault in question
!
!     the 830 loop determines the location of the fault
!
        do 830 l=1,bec(k)
          fnode=fnode+nr(l,k)
  830   continue
        xfault(k) = xsp(fnode)
        fnode=fnode+nr(l,k)
        nfn=0
        do 840 l=1,nnode
	        fflag=0
         	xx = x(l) - faultx(k)
          zz = faultz(k)-z(l)
          alpha4=90.d0-(-delta+omega)
           if(zz.eq.0.d0)then
	           if(xx.eq.0.d0)then
             fflag=1
             endif
           else if(sngl(xx).eq.0.d0.and.
     $          rotflag.eq.1)then
           fflag=1

           axz = sngl(xx/zz)
           adtz = sngl(dtand(alpha4))
           adif1 = abs(axz-adtz)
           else if(sngl(xx/zz).eq.sngl(dtand(alpha4)))then
	    fflag=1
            endif

              if(fflag.eq.1) then
                      nfn=nfn+1
                      zoor(nfn)=z(l)
                      oor(nfn)=l
                      nor(nfn)=l
                  end if
840       continue

        do 860 m=1,nfn
          lnode=1
          olow=1e+10
          do 880 l=1,nfn+1-m
            nlow=dmax1(zoor(l),zoor(lnode))
            diflow = dabs(nlow-olow)
            if(diflow.gt.1.e-03)then
              olow=nlow
              lnode=l
            end if
            last=l
  880     continue
          zoor(lnode)=zoor(last)
          zoor(last)=olow
          nortemp=nor(lnode)
          nor(lnode)=nor(last)
          nor(last)=nortemp
  860   continue
!
!     this section of is used to reassigne the nodal
!     properties for the nodes which have been rearanged
!	the nodes have to be rearaanfged for best ordering
!
        tdcount=1
        do 900 l=1,nfn
          tdata(tdcount)   = z(nor(l))
          tdata(tdcount+1) = vsold(nor(l))
          tdata(tdcount+2) = phild1(nor(l))
          tdata(tdcount+3) = phild2(nor(l))
          tdata(tdcount+4) = hold(nor(l))
          tdata(tdcount+5) = tld1(nor(l))
          tdata(tdcount+6) = tld2(nor(l))
          tdata(tdcount+7) = tcold(nor(l))
          tdata(tdcount+8) = ttiold(nor(l))
          tdata(tdcount+9) = tticld(nor(l))
          tdata(tdcount+10) = x(nor(l))
          tdata(tdcount+11) = nd_agld(nor(l))
          tdata(tdcount+13) = sfcelold(nor(l))
          tdata(tdcount+14) = delsfcold(nor(l))
          tdata(tdcount+15) = concld(nor(l))
          tdata(tdcount+16) = prsltmp1(nor(l))
          tdata(tdcount+17) = prsltmp2(nor(l))
          tdcount=tdcount+18
          if(nor(l).lt.0) write(*,*)'back to begin'
  900   continue
        tdcount=1
        do 910 l=1,nfn
          z(oor(l))   =  tdata(tdcount)
          vsold(oor(l))  =  tdata(tdcount+1)
          phild1(oor(l)) =  tdata(tdcount+2)
          phild2(oor(l)) =  tdata(tdcount+3)
          hold(oor(l))   =  tdata(tdcount+4)
          tld1(oor(l))   =  tdata(tdcount+5)
          tld2(oor(l))   =  tdata(tdcount+6)
          tcold(oor(l))  =  tdata(tdcount+7)
          ttiold(oor(l)) =  tdata(tdcount+8)
          tticld(oor(l)) =  tdata(tdcount+9)
          x(oor(l))   =  tdata(tdcount+10)
          nd_agld(oor(l))=tdata(tdcount+11)
          sfcelold(oor(l))=tdata(tdcount+13)
          delsfcold(oor(l))=tdata(tdcount+14)
          concld(oor(l))   =  tdata(tdcount+15)
          prsltmp1(oor(l))   =  tdata(tdcount+16)
          prsltmp2(oor(l))   =  tdata(tdcount+17)
          tdcount=tdcount+18

  910   continue
        tdcount=0

            xfault(k) = x(oor(1))
        do 920 m=1,2
          if(m.eq.2)then
            fcel=tcel+bec(k)*ber(k)*2+1
            anc=1
          else
            fcel=tcel+bec(k)*2
            anc=bec(k)+1
          end if
          ncel=ber(k-1+m)*2
          cel=fcel
          if(m.eq.2)then
            bwfw=+1
            sw=1
          else
            bwfw=-1
            if(mod(bec(k-1+m),2).eq.0)then
              sw=2
            else
              sw=1
            end if
          end if
          do 940 l=1,ncel,2
            if(m.eq.1)then
              if(mod(sw,2).eq.0)then
                a=1
                b=2
              else
                a=2
                b=3
              end if
              if(eltp(cel-1).eq.6)then
                if(nr(bec(k)+1,k).gt.nr(bec(k),k))then
                  a=2
                  b=3
                end if
              end if
              c=2
            else
              a=1
              b=3
              if(mod(sw,2).eq.0)then
                c=3
              else
                c=1
              end if
            end if
!
!     this section reassigns the old node numbers with the new node
!     numbers in the conectivity matrix for elements along the fault
!
            tr9=0
            tr10=0
            tr11=0
            tr12=0
            do 960 n=1,nfn
              if(tr9.ne.1)then
                if(in(a,cel).eq.nor(n))then
                  in(a,cel)=oor(n)
                  if((m.eq.1).and.(cel.eq.fcel))then
                    lbnd(k)=oor(n)
                  end if
                  if((m.eq.2).and.(cel.eq.fcel))then
                    rbnd(k)=oor(n)
                  end if
                  tr9=1
                end if
              end if
              if(tr10.ne.1)then
                if(in(b,cel).eq.nor(n))then
                  in(b,cel)=oor(n)
                  tr10=1
                end if
              end if
              if(tr12.ne.1)then
                if(in(4,cel).ne.0)then
                  if(in(4,cel).eq.nor(n))then
                    in(4,cel)=oor(n)
                    tr12=1
                  end if
                end if
              end if
              if(tr11.ne.1)then
                if(eltp(cel+bwfw).eq.6)goto 960
                if(in(c,cel+bwfw).eq.nor(n))then
                  in(c,cel+bwfw)=oor(n)
                  tr11=1
                end if
              end if
  960       continue
            cel = cel+(bec(k-1+m))*2
            sw=sw+1
  940     continue
  920   continue
        tcel=tcel+bec(k)*ber(k)*2
!
!       loop renumbers nrows and reduces ncols by 1
!    does it need to lower the col numbering to account for the fault col???
!
        tncols=tncols-1
        fcol=fcol+bec(k)
        tnrows(fcol)=tnrows(fcol)+tnrows(fcol+1)-ncf
        do 970 n=fcol+1,tncols
          tnrows(n)=tnrows(n+1)
  970   continue
  350 continue

      tcel=0
      tncf=tncf+ncf
      ncf=0

      do 980 k=1,nelem
      ni(k)=in(1,k)
      if(ni(k).lt.0)then
		write(*,*)'ni lt 0',k,x(ni(k)),z(ni(k))
  		ni(k)=0
      endif
      nj(k)=in(2,k)
       if(nj(k).lt.0) nj(k)=0
      nk(k)=in(3,k)
       if(nk(k).lt.0)  nk(k)=0
      nl(k)=in(4,k)
       if(nl(k).lt.0)  nl(k)=0
  980 continue
      maxdiff=0
      do 985 m=1,nelem
        i=in(1,m)
        j=in(2,m)
        k=in(3,m)
        l=in(4,m)
        if(l.eq.0)then
          l = k
        end if
        maxnode=max(i,j,k,l)
        minnode=min(i,j,k,l)
        diff=maxnode-minnode
        if(diff.gt.maxdiff)then
          maxdiff=diff
          if(2.0*maxdiff+1.gt.maxbnds)then
          write(8,*)i,j,k,l,m
          write(8,*)x(i),x(j),x(k)
          write(8,*)z(i),z(j),z(k)
          write(8,*)eltp(m)
          write(*,*)'stop7'
          end if
        end if
  985 continue
      nfband=maxdiff+1
      nhband=(maxdiff*2)+1
!
!     assign nodal material porosities
!
      do 990 l=1,nelem
        if(ni(l).eq.0) go to 990
        matp(ni(l))=mat(l)
        matp(nj(l))=mat(l)
        matp(nk(l))=mat(l)
        if(nl(l).gt.0) matp(nl(l))=mat(l)
  990 continue
!
!     check to see if the no. of elements exceeds the number of
!     assigned material properties.  stop the program if this occurs!
!
      qq = mat(nelem)
        if(qq.eq.0) then
         write(8,503) nelem
          write(*,*)'stop8'
         stop
      end if
!
! reseting variables back to their correct (non old) names
!
       do 575 j=1,nnode
         vs(j)=vsold(j)
         phi(j)=phild1(j)
         phiold(j)=phild2(j)
         head(j)=hold(j)
         conc(j)=concld(j)
         temp(j)=tld1(j)
        told(j)=tld2(j)
        tempc(j)=tcold(j)
        tti(j)=ttiold(j)
        ttic(j)=tticld(j)
        node_angle(j)=nd_agld(j)
        sfcel(j)=sfcelold(j)
        delsfc(j)=delsfcold(j)
        presl(j) = prsltmp1(j)
        prslold(j) = prsltmp2(j)
  575 continue

  503 format (1h1,//,3x,
     $ '**** warning:  number of elements exceeds no. of
     $  assigned material properties ***',/,3x,'nelem = ',i5,
     $  5x,'aborting simulation')
      return
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      subroutine fbsep (xsp,zsp,x,z,head,vs,temp,
     $told,tempc,tti,ttic,nnode,tncf,matp,ncols,nrows,phi,
     $phiold,node_angle,sfcel,delsfc,conc,presl,prslold)
!
!     this subroutine separates fault block and global spatial
!     coordinates
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!


      real*8 node_angle(maxnodes),nd_agtmp(maxnodes),xsp(maxnodes)
      real*8 head(maxnodes),htmp(maxnodes),vs(maxnodes),vstmp(maxnodes)
      real*8 ttmp1(maxnodes),told(maxnodes),ttmp2(maxnodes)
      real*8 tctmp(maxnodes),tti(maxnodes),ttitmp(maxnodes)
      real*8 ttic(maxnodes),ttictmp(maxnodes), phi(maxnodes)
      real*8 phitmp1(maxnodes),zsp(maxnodes),x(maxnodes),z(maxnodes)
      real*8 sfcel(maxnodes),sfceltmp(maxnodes),delsfc(maxnodes)
      real*8 conc(maxnodes),conctmp(maxnodes),temp(maxnodes)
      real*8 presl(maxnodes),prslold(maxnodes),prsltmp1(maxnodes)
      real*8 phitmp2(maxnodes),phiold(maxnodes),delsfctmp(maxnodes)
      real*8 prsltmp2(maxnodes),tempc(maxnodes)

      integer nnode,tncf,ncols
      integer matptmp(maxnodes),nrows(maxbnds),matp(maxnodes)

      do 546 m=1,nnode+tncf
            htmp(m)=head(m)
            conctmp(m)=conc(m)
            phitmp1(m)=phi(m)
            phitmp2(m)=phiold(m)
            vstmp(m)=vs(m)
            ttmp1(m)=temp(m)
            ttmp2(m)=told(m)
            tctmp(m)=tempc(m)
            ttitmp(m)=tti(m)
            ttictmp(m)=ttic(m)
            matptmp(m)=matp(m)
            nd_agtmp(m)= node_angle(m)
            sfceltmp(m)=sfcel(m)
            delsfctmp(m)=delsfc(m)
            prsltmp1(m)= presl(m)
            prsltmp2(m)= prslold(m)
  546 continue
      do 50 l=1,nnode
        matptmp(l)=matp(l)
   50 continue
c
c saf is a saftey factor to see if two nodes are in the same location
c
      l=0
      do 100 n=1,ncols
        saf=0.09d+00
        do 100 nn=1,nrows(n)
          l=l+1
          if(nn.eq.nrows(n))saf=0.00d+00
        do 100 m=1,nnode
         axdif = dabs(xsp(l)-x(m))
          if((axdif.le.1.0e-3).and.(zsp(l).ge.z(m)-saf).and.
     $      (zsp(l).le.z(m)+saf))then
            head(l)=htmp(m)
            conc(l)=conctmp(m)
            phi(l)=phitmp1(m)
            phiold(l)=phitmp2(m)
            vs(l)=vstmp(m)
            temp(l)=ttmp1(m)
            told(l)=ttmp2(m)
            tempc(l)=tctmp(m)
            tti(l)=ttitmp(m)
            ttic(l)=ttictmp(m)
            matp(l)=matptmp(m)
            node_angle(l)=nd_agtmp(m)
            sfcel(l)=sfceltmp(m)
            delsfc(l)=delsfctmp(m)
            presl(l)=prsltmp1(m)
            prslold(l)= prsltmp2(m)
          end if
  100 continue
      do 200 l=1,nnode+tncf
        htmp(l)=0.00d+00
        conctmp(l)=0.00d+00
        phitmp1(l)=0.00d+00
        phitmp2(l)=0.00d+00
        vstmp(l)=0.00d+00
        ttmp1(l)=0.00d+00
        ttmp2(l)=0.00d+00
        tctmp(l)=0.00d+00
        ttitmp(l)=0.00d+00
        ttictmp(l)=0.00d+00
        matptmp(l)=0
        nd_agtmp(l)=0.00d+00
        sfceltmp(l)=0.d+00
        delsfctmp(l)=0.d+00
        prsltmp1(m) = 0.d+00
        prsltmp2(m) = 0.d+00
  200 continue
      return
      end
!**************** element output ***********************************
      subroutine pelemout(nnr,nrel,it,mat,xold,zold,ni,nj,nk,qz,qx
     $   ,tld1,hold,ntime,iprint,msl,dt,nelem1)
!                                                                   *
!	this subroutine outputs values for specific elements        *
!                                                                   *
!********************************************************************
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

	integer mat(maxelems),nnr,nrel(maxflts),it
        integer ni(maxelems),nj(maxelems)
	integer ntime,iprint,nk(maxelems),ck

        real*8 xold(maxnodes),zold(maxnodes),qz(maxelems),qx(maxelems)
	real*8 tempc,heads,msl,dt,tld1(maxnodes),hold(maxnodes),x,z
        
          if(it.eq.1) then
          do 100 n=1,nnr
          lll = n+65
          write(lll,56) nrel(n),ntime
  100 continue
	  endif
       do 223 ck=1,nnr
          if(nrel(ck).gt.nelem1) go to 223
          if(ni(nrel(ck)).le.0) go to 223
          if(nj(nrel(ck)).le.0) go to 223
          if(nk(nrel(ck)).le.0) go to 223
          lll = ck+65
          x=xold(ni(nrel(ck)))/3+xold(nj(nrel(ck)))/3
     $      +xold(nk(nrel(ck)))/3
          z=zold(ni(nrel(ck)))/3+zold(nj(nrel(ck)))/3
     $      +zold(nk(nrel(ck)))/3
          tempc=tld1(ni(nrel(ck)))/3+tld1(nj(nrel(ck)))/3
     $      +tld1(nk(nrel(ck)))/3
          heads=hold(ni(nrel(ck)))/3+hold(nj(nrel(ck)))/3
     $      +hold(nk(nrel(ck)))/3


          time = dble(it)*dt
          write(lll,*)time,tempc,qx(nrel(ck)),qz(nrel(ck)),
     $          x,z,mat(nrel(ck)),heads
 223   continue
   56 format ('title = "elemental data"',/,'variables = "time (yr)", "te
     $mp (deg. c)","qx(m/yr)","qz(m/yr)","x(m)","z(m)","mat","head(m)"'
     $ ,/,'text  x=50.0,y=50.0,t="element #',i6,' "',/,
     $   'zone f="point",i=',i6)
        return
	end
!***********************************************************************
      subroutine kinetic (isourc,mat,it,nelem,temp,told,pres,
     $   perct1,perct2,perct3,dt,rv,tr,ni,nj,nk,toc,phi,area,head,
     $   oilvol,gasvol,persat,mash2o,masco2,masoil,masgas,fracoil,
     $   kgoil,kggas,trold)
!
!
!      this subroutine compute source rock maturation using first-order
!      rate kinetic theory.  this subroutine calculates the
!      transformation ratio (tr); this information can be used to
!      determine the "oil window" (10%<tr<40% is the oil window).  it also
!      calculates vitrinite reflectance (rv), the mass and volume of
!      oil and gas generated, and the percent saturation.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!

       integer isourc,mat,qq,i,j,it,l,nelem,ni,nk,nj
       integer rtn1,rtn2,rtnh2o,rtnco2
       integer rtnoil,rtnch4,rtntgs
       integer ioilck(maxelems)

       real*8 temp,newtmp,oldtmp,told,perct1,perct2,perct3
       real*8 fa,fab,fb,simpsn,dt,x1,x2,x3h2o,x3co2
       real*8 x3oil,x3ch4,xtgs,oldx1,oldx2,oldh2o,oldco2,oldoil
       real*8 oldtgs,oldch4,newx1,newx2,newh2o
       real*8 newco2,newoil,newch4,newtgs,e1,e2,e3h2o,e3co2
       real*8 e3ch4,etgs,e3oil,arr1,arr2,arrh2o,arrco2,arroil
       real*8 arrch4,arrtgs,r,trnum,trden,tr1,tr2
       real*8 tr3oil,tr
       real*8 masoil,mash2o,masco2,masch4,mastgs,masgas
       real*8 masol3,masgs3,newol1,newol2,newol3
       real*8 oldol1,oldol2,oldol3,newwat,newcd,newmet,newgas
       real*8 oldwat,oldcd,oldmet,oldgas
       real*8 genh2o,genco2,genoil3,gench4,genoil1,genoil2
       real*8 falpha,fbeta,fgamma,fdelta,hoc,ooc,rv
       real*8 moil,mgas,toc,phie,phi,rockdens,width,maselem
       real*8 area,mgoil,mggas,head,pres,pelem,oildens,kggas
       real*8 gasdens,oilvol,gasvol,porvol,persat,fracoil,kgoil
       real*8 trold(maxelems)

       dimension ni(maxelems),nk(maxelems),nj(maxelems),mat(maxelems)
       dimension phi(maxnodes),area(maxelems),pres(maxnodes)
       dimension perct1(0:maxmats),perct2(0:maxmats),perct3(0:maxmats)
       dimension masoil(maxelems),mash2o(maxelems),masco2(maxelems)
       dimension toc(0:maxmats),oilvol(maxelems),gasvol(maxelems)
       dimension head(maxnodes),kgoil(maxelems),kggas(maxelems)
       dimension temp(maxnodes),told(maxnodes),isourc(0:maxmats)
       dimension tr(maxelems),rv(maxelems),masgas(maxelems)
       dimension persat(maxelems),fracoil(maxelems)

       dimension x1(7),x2(10),x3h2o(8),x3co2(7),oldtgs(3)
       dimension x3oil(7),x3ch4(13),xtgs(3)
       dimension e1(7),e2(10),e3h2o(8),e3co2(7),e3oil(7)
       dimension e3ch4(13),etgs(3)
       dimension oldx1(maxelems,7),oldx2(maxelems,10)
       dimension oldco2(maxelems,7),oldoil(maxelems,7)
       dimension masch4(maxelems),mastgs(maxelems),oldh2o(maxelems,8)
       dimension masol3(maxelems),masgs3(maxelems),oldch4(maxelems,13)
       dimension oldol1(maxelems),oldol2(maxelems),oldol3(maxelems)
       dimension oldwat(maxelems),oldcd(maxelems),oldmet(maxelems)
       dimension newx1(7),newx2(10),newh2o(8),newco2(7)
       dimension newoil(7),newch4(13),newtgs(3)

!
!  assign values to variables first
!
       data rtn1,rtn2,rtnh2o,rtnco2,rtnoil / 7,10,8,7,7/
       data rtnch4,rtntgs /13, 3/
       data r / .001987/
       data arr1,arr2,arrh2o,arrco2/ 1e+14,1e+15,2*1e+13 /
       data arrch4 / 1e+13/
       data arroil,arrtgs/ 2e+13, 1e+15 /
!
!  genetic potentials, xi, for type i, ii, iii (h2o,co2,oil,and ch4)
!  and thermal gas
!
       data x1(1),x1(2),x1(3),x1(4) / 5.,10.,15.,20./
       data x1(5),x1(6),x1(7) / 800.,10.,15./
       data x2(1),x2(2),x2(3),x2(4) / 5.,3.,10.,30./
       data x2(5),x2(6),x2(7),x2(8) /150.,325.,60.,20./
       data x2(9),x2(10) / 10.,5./
       data x3h2o(1),x3h2o(2) /4.4,8.8/
       data x3h2o(3),x3h2o(4) /13.2,17.6/
       data x3h2o(5),x3h2o(6) /17.6,13.2/
       data x3h2o(7),x3h2o(8) /8.8,4.4 /
       data x3co2(1),x3co2(2) /1.45,4.35/
       data x3co2(3),x3co2(4) / 72.5,72.5/
       data x3co2(5),x3co2(6), x3co2(7) /43.5,29.0,14.5/
       data x3oil(1),x3oil(2) / .35,.7/
       data x3oil(3),x3oil(4) /1.4,2.1 /
       data x3oil(5),x3oil(6) /1.4,0.7/
       data x3oil(7),x3ch4(1) /0.35,7.56/
       data x3ch4(2), x3ch4(3),x3ch4(4)/12.96,15.12,14.04/
       data x3ch4(5),x3ch4(6) /12.96,11.88/
       data x3ch4(7),x3ch4(8) /9.72,7.56/
       data x3ch4(9),x3ch4(10) /5.4,4.32/
       data x3ch4(11),x3ch4(12) /3.24,2.16/
       data x3ch4(13) / 1.08 /
       data xtgs(1),xtgs(2),xtgs(3) / 0.5, 2.0, 7.5/
!
!  activation energies, ei, for type i, ii, iii (h2o,co2,oil,
!  and ch4) and thermal gas...
!
       data e1(1),e1(2),e1(3),e1(4) / 48.0,50.0,52.0,54.0 /
       data e1(5),e1(6),e1(7) / 56.0,58.0,60.0 /
       data e2(1),e2(2),e2(3),e2(4) / 40.0,42.0,48.0,50.0 /
       data e2(5),e2(6),e2(7),e2(8) / 52.0,54.0,56.0,58.0 /
       data e2(9),e2(10) / 60.0,62.0/
       data e3h2o(1),e3h2o(2), e3h2o(3),e3h2o(4)/ 38.0,40.0,42.0,44.0/
       data e3h2o(5),e3h2o(6), e3h2o(7),e3h2o(8)/46.0,48.0,50.0,52.0/
       data e3co2(1),e3co2(2), e3co2(3),e3co2(4)/ 42.0,44.0,46.0,48.0/
       data e3co2(5),e3co2(6), e3co2(7) /50.0,52.0,54.0/
       data e3oil(1),e3oil(2), e3oil(3),e3oil(4)/ 46.0,48.0,50.0,52.0/
       data e3oil(5),e3oil(6), e3oil(7) /54.0,56.0,58.0/
       data e3ch4(1),e3ch4(2), e3ch4(3),e3ch4(4)/50.0,52.0,54.0,56.0/
       data e3ch4(5),e3ch4(6), e3ch4(7),e3ch4(8)/58.0,60.0,62.0,64.0/
       data e3ch4(9),e3ch4(10)/66.0,68.0/
       data e3ch4(11),e3ch4(12)/70.0,72.0/
       data e3ch4(13)/ 74.0 /
       data etgs(1),etgs(2),etgs(3)/ 59.9,64.7,69.5 /

         if(it.eq.1) then
         do 455 kk=1,nelem
          ioilck(kk)=0
 455     continue
         endif



!
!  convert dt to units of seconds
       dt = dt*3.15576e+07
!
!  first, we need to determine the total genetic potentials
!  for type 3 (for vitrinite reflectance calculations...) and for
!  types 1 and 2 (only need to be computed on the first iteration.
!
       if (it .eq. 1) then
          do 91 j = 1,rtnh2o
             genh2o = genh2o + x3h2o(j)
91        continue
          do 92 j = 1,rtnco2
             genco2 = genco2 + x3co2(j)
92        continue
          do 93 j = 1,rtnoil
             genoil3 = genoil3 + x3oil(j)
93        continue
          do 94 j = 1,rtnch4
             gench4 = gench4 + x3ch4(j)
94        continue
          do 95 j = 1, rtn1
             genoil1 = genoil1 + x1(j)
95        continue
          do 96 j = 1, rtn2
             genoil2 = genoil2 + x2(j)
96        continue
       endif
!
!  set initial value of "old" number of elements
!  if we've added new elements at this time step, we need to assign
!  initial values to these elements before entering the loop.
!
          do 40 i = 1, nelem
             if(ni(i).eq.0) goto 40
             if(ioilck(i).eq.1) goto 40
             ioilck(i)=1
             qq = mat(i)
             if (perct1(qq) .gt. 0) then
                do 51 j = 1,rtn1
                   oldx1(i,j) = x1(j)
                   oldol1(i) = (oldol1(i) + x1(j))
51              continue
                oldol1(i) = oldol1(i) * perct1(qq)
             endif
             if (perct2(qq) .gt. 0) then
                do 52 j = 1,rtn2
                   oldx2(i,j) = x2(j)
                   oldol2(i) = (oldol2(i) + x2(j))
52              continue
                oldol2(i) = oldol2(i) * perct2(qq)
             endif
             if (perct3(qq) .gt. 0) then
                do 53 j = 1,rtnh2o
                   oldh2o(i,j) = x3h2o(j)
                   oldwat(i) = (oldwat(i) + x3h2o(j))
53              continue
                do 54 j = 1,rtnco2
                   oldco2(i,j) = x3co2(j)
                   oldcd(i) = (oldcd(i) + x3co2(j))
54              continue
                do 55 j = 1,rtnoil
                   oldoil(i,j) = x3oil(j)
                   oldol3(i) = (oldol3(i) + x3oil(j))
55              continue
                do 56 j = 1,rtnch4
                   oldch4(i,j) = x3ch4(j)
                   oldmet(i) = (oldmet(i) + x3ch4(j))
56              continue
                oldwat(i) = oldwat(i)*perct3(qq)
                oldcd(i) = oldcd(i)*perct3(qq)
                oldol3(i) = oldol3(i)*perct3(qq)
                oldmet(i) = oldmet(i)*perct3(qq)
             endif
40        continue
!
!             main body of subroutine (100-loop)
!  outer loop looks at each element ...
!  if the element is a source element, then the kinetics are computed.
!  { if the element doesn't yet exist, jump to end of loop}
!
      do 100 m=1,nelem
         qq=mat(m)
         if (isourc(qq).eq.0) go to 700
!
!  first, compute elemental temperature at current and past time step..
!
         i = ni(m)
         j = nj(m)
         k = nk(m)
         if(i.eq.0) go to 100
         newtmp = (temp(i)+temp(j)+temp(k))/3.d+0 + 273.15
         oldtmp = ((told(i)+told(j)+told(k))/3) + 273.15
!
!  type i:  next, compute mass fraction of oil generated
!
         if (perct1(qq) .gt. 0) then
!
!  next, perform kinetic calculation for oil generation
!
            do 120 l=1,rtn1
               fa=dexp(-e1(l) / (r*oldtmp))
               fb=dexp(-e1(l) / (r*newtmp))
               fab=dexp(-e1(l) / (r*((newtmp+oldtmp)/2)))
               simpsn= (dt/6)*(fa + 4*fab + fb)
               newx1(l)=oldx1(m,l)*dexp(-arr1*simpsn)
               trnum=trnum+(x1(l)-newx1(l))
               trden=trden+x1(l)
120     continue
!
!  now, compute transformation ratio and add to running total
!  for type 1 o.m. for this element...
!
            tr1=(trnum/trden)*perct1(qq)
            trnum=0
            trden=0
         endif
!
!  type ii:  next, compute the mass fraction of oil generated
!
           if (perct2(qq) .gt. 0) then
!
!  next, perform kinetic calculation for oil generation
!
            do 140 l=1,rtn2
               fa=dexp(-e2(l) / (r*oldtmp))
               fb=dexp(-e2(l) / (r*newtmp))
               fab=dexp(-e2(l) / (r*((newtmp+oldtmp)/2)))
               simpsn= (dt/6)*(fa + 4*fab + fb)
               newx2(l)=(oldx2(m,l)*dexp(-arr2*simpsn))
               trnum=trnum+(x2(l)-newx2(l))
               trden=trden+x2(l)
140     continue
!
!  now, compute transformation ratio and add to running total
!  for type 2 o.m. for this element...
!
            tr2=(trnum/trden)*perct2(qq)
            trnum=0
            trden=0
          endif
!
!  type iii: next, compute the mass fraction of h2o, co2, oil and ch4.
!
             if (perct3(qq) .gt. 0) then
!
!  next, perform kinetic calculation for h2o generation (type 3)
!
            do 190 l=1,rtnh2o
               fa=dexp(-e3h2o(l) / (r*oldtmp))
               fb=dexp(-e3h2o(l) / (r*newtmp))
               fab=dexp(-e3h2o(l) / (r*((newtmp+oldtmp)/2)))
               simpsn= (dt/6)*(fa + 4*fab + fb)
               newh2o(l)=oldh2o(m,l)*dexp(-arrh2o*simpsn)
               trnum=trnum+(x3h2o(l)-newh2o(l))
               trden=trden+x3h2o(l)
190    continue
!
!  now, compute transformation ratio and add to running total
!  for type 3 o.m. (h2o) for this element...
!
            trnum=0
            trden=0
!
!  next, perform kinetic calculation for co2 generation (type 3)
!
            do 200 l=1,rtnco2
               fa=dexp(-e3co2(l) / (r*oldtmp))
               fb=dexp(-e3co2(l) / (r*newtmp))
               fab=dexp(-e3co2(l) / (r*((newtmp+oldtmp)/2)))
               simpsn= (dt/6)*(fa + 4*fab + fb)
               newco2(l)=oldco2(m,l)*dexp(-arrco2*simpsn)
               trnum=trnum+(x3co2(l)-newco2(l))
               trden=trden+x3co2(l)
200         continue
!
!  now, compute transformation ratio and add to running total
!  for type 3 o.m. (co2) for this element...
!
            trnum=0
            trden=0
!
!  next, perform kinetic calculation for oil generation (type 3)
!
            do 210 l=1,rtnoil
               fa=dexp(-e3oil(l) / (r*oldtmp))
               fb=dexp(-e3oil(l) / (r*newtmp))
               fab=dexp(-e3oil(l) / (r*((newtmp+oldtmp)/2)))
               simpsn= (dt/6)*(fa + 4*fab + fb)
               newoil(l)=oldoil(m,l)*dexp(-arroil*simpsn)
               trnum=trnum+(x3oil(l)-newoil(l))
               trden=trden+x3oil(l)
210         continue
!
!  now, compute transformation ratio and add to running total
!  for type 3 o.m. for this element...
!
            tr3oil=(trnum/trden)*perct3(qq)
            trnum=0
            trden=0
!
!  finally, perform kinetic calculation for ch4 generation (type 3)
!
            do 220 l=1,rtnch4
               fa=dexp(-e3ch4(l) / (r*oldtmp))
               fb=dexp(-e3ch4(l) / (r*newtmp))
               fab=dexp(-e3ch4(l) / (r*((newtmp+oldtmp)/2)))
               simpsn= (dt/6)*(fa + 4*fab + fb)
               newch4(l)=oldch4(m,l)*dexp(-arrch4*simpsn)
               trnum=trnum+(x3ch4(l)-newch4(l))
               trden=trden+x3ch4(l)
220         continue
!
!  now, compute transformation ratio and add to running total
!  for type 3 o.m. (ch4) for this element...
!
            trnum=0
            trden=0
       endif
!
!  this marks the end of kinetic computations for primary cracking..
!  the next step is to add up all the oil formed during this time step and add it
!   to a running total of oil formed called masoil.
!
         do 221 l = 1, rtn1
             newol1 = newol1 + newx1(l)
221      continue
         do 222 l = 1, rtn2
             newol2 =  newol2 + newx2(l)
222      continue
         do 223 l = 1, rtnoil
             newol3 = newol3 + newoil(l)
223      continue
         newol1 = newol1*perct1(qq)
         newol2 = newol2*perct2(qq)
         newol3 = newol3*perct3(qq)
         masoil(m)=masoil(m) + ( (oldol1(m)-newol1) +
     $     (oldol2(m)-newol2)+ (oldol3(m) - newol3) )
!
!  now compute mash2o, masco2 and masch4 the same way...
!
         do 224 l = 1, rtnh2o
            newwat = newwat + newh2o(l)
224      continue
         newwat = newwat*perct3(qq)
         mash2o(m)=mash2o(m)+(oldwat(m)-newwat)
         do 225 l = 1, rtnco2
            newcd = newcd + newco2(l)
225      continue
         newcd = newcd*perct3(qq)
         masco2(m)=masco2(m)+(oldcd(m)-newcd)
         do 226 l = 1, rtnch4
            newmet = newmet + newch4(l)
226      continue
         newmet = newmet*perct3(qq)
         masch4(m)=masch4(m)+(oldmet(m)-newmet)
c
c  next, compute masol3 and masgs3 for vitrinite reflectance calcs.
c
         masol3(m) = masol3(m) + (oldol3(m) - newol3)
         masgs3(m) = masgs3(m) + (oldmet(m) - newmet)
!
! vitrinite reflectance calculations
!
      if (perct3(qq) .gt. 0.0) then
         falpha = mash2o(m) / (genh2o * perct3(qq))
         fbeta = masco2(m) / (genco2 * perct3(qq))
         fgamma = masol3(m) / (genoil3 * perct3(qq))
         fdelta = masgs3(m) / (gench4 * perct3(qq))
         hoc=(.9-(2*.35*.25*falpha)-(1.7*.01*fgamma)-(4*.125*fdelta))
     $       /(1-(.35*.7*fbeta)/2-.01*fgamma-.125*fdelta)
         ooc=.35*(1-(.25*falpha)-(.7*fbeta))/(1-(.35*.7*fbeta)/2
     $       -(.01*fgamma)-(.125*fdelta))
         rv(m) = 12 * dexp(-3.3 * hoc) - ooc
      endif
!
!  thermal gas calculations
!
         do 240 l=1,rtntgs
            fa=dexp(-etgs(l) / (r*oldtmp))
            fb=dexp(-etgs(l) / (r*newtmp))
            fab=dexp(-etgs(l) / (r*((newtmp+oldtmp)/2)))
            simpsn= (dt/6)*(fa + 4*fab + fb)
            oldtgs(l) =xtgs(l)*(masoil(m)*.9)
            newtgs(l)=oldtgs(l)*dexp(-arrtgs*simpsn)
240      continue
!
!  the next step is to add up all the thermal gas formed during
!  this time step (mastgs)
!
         do 241 l = 1, rtntgs
            newgas = newgas + newtgs(l)
            oldgas = oldgas + oldtgs(l)
241      continue
         mastgs(m)=mastgs(m) +(oldgas-newgas)
         masoil(m) = masoil(m) - (oldgas-newgas)
         trold(m) = tr(m)
         tr(m)=(tr1+tr2+tr3oil)
         masgas(m) = mastgs(m) + masch4(m)
!
!  compute mass of oil and gas per gram of source rock
!
          moil = masoil(m) * toc(qq)
          mgas = masgas(m) * toc(qq)
!
! compute mass of source rock in element (grams)
!
          i = ni(m)
          j = nj(m)
          k = nk(m)
          phie = (phi(i) + phi(j) + phi(k) ) / 3.0d+0
          rockdens = 2.65d+6
          width = 1.0d+0
          maselem = area(m) * width * (1.d+0 - phie ) * rockdens
!
!  compute mass of oil and gas in element (in milligrams)
!
          mgoil = moil * maselem
          mggas = mgas * maselem
!
!  compute a new array of masses of oil and gas in element for
!  use in plot file, kgoil, kggas {note:  all references to these
!  arrays have been deleted, but the computation is here...}
!
           kgoil(m) = mgoil / 1.0d+06
           kggas(m) = mggas / 1.0d+06
!
!  compute the relative fraction of kerogen that has been transformed
!  to oil and gas.
!
           fracoil(m) = 0.0
           if ( (perct1(qq).gt.0) .or. (perct2(qq).gt.0)
     $         .or. (perct3(qq).gt.0) )  then
              fracoil(m) = masoil(m) / ( (genoil1*perct1(qq)) +
     $           (genoil2*perct2(qq)) + (genoil3*perct3(qq)) )
           endif
!  compute density (kg/m^3) of oil and gas in element using equations
!  of state (relationships by england et. al, 1987)
!
          pelem = (pres(i) + pres(j) + pres(k)) / 3.0d+6
          oildens = 781.64 - 7.8944*pelem +
     $               4.1642*dexp(-2*(pelem**2))
          gasdens = -43.444 + (10.281*pelem)
!
!  change to mg/m3
!
          oildens = oildens * 1.d+06
          gasdens = gasdens * 1.d+06
!
! compute volume of oil and gas in element (in m^3)
!
          oilvol(m) = mgoil / oildens
          gasvol(m) = mggas / gasdens
!
!  compute percent saturation of oil+gas in element
!
          porvol = area(m) * width * phie
          persat(m) = ((oilvol(m) + gasvol(m))/porvol)*100.d+0
!
!  once each element is processed at current temperature, reassign
!  current-time variables to previous-time variables...
!
          oldol1(m) = newol1
          oldol2(m) = newol2
          oldol3(m) = newol3
          oldwat(m) = newwat
          oldcd(m) = newcd
          oldmet(m)= newmet
c
          newol1 = 0
          newol2 = 0
          newol3 = 0
          newwat = 0
          newcd = 0
          newmet = 0
          newgas = 0
          oldgas = 0
!
!  we need to switch arrays used in kinetic computations as well..
!
          do 400 l = 1, rtn1
               oldx1(m,l) = newx1(l)
400       continue
          do 410 l = 1, rtn2
               oldx2(m,l) = newx2(l)
410       continue
          do 420 l = 1, rtnh2o
               oldh2o(m,l) = newh2o(l)
420       continue
          do 430 l = 1, rtnco2
               oldco2(m,l) = newco2(l)
430       continue
          do 440 l = 1, rtnoil
               oldoil(m,l) = newoil(l)
440       continue
          do 450 l = 1, rtnch4
               oldch4(m,l) = newch4(l)
450       continue
!
!   now, move on to process next source element at current temperature
!
700    continue
100    continue
       dt = dt/3.15576e+07
       if (it .eq. 1) then

       write (8,704)
       write (8,705) 'summary of input values for kinetic calculations'
       write (8,705) 'type i kinetic parameters'
       write (8,730) 'ei', 'xi'
       do 760 j = 1,rtn1
            write (8,740) e1(j), x1(j)
760    continue
       write (8,720)
       write (8,705) 'type ii kinetic parameters'
       write (8,730) 'ei', 'xi'
       do 765 j = 1,rtn2
            write (8,740) e2(j), x2(j)
765    continue
       write (8,720)
       write (8,705) 'type iii kinetic parameters: h2o'
       write (8,730) 'ei', 'xi'
       do 770 j = 1,rtnh2o
            write (8,740) e3h2o(j), x3h2o(j)
770    continue
       write (8,720)
       write (8,705) 'type iii kinetic parameters: co2'
       write (8,730) 'ei', 'xi'
       do 775 j = 1,rtnco2
            write (8,740) e3co2(j), x3co2(j)
775    continue
       write (8,720)
       write (8,705) 'type iii kinetic parameters: oil'
       write (8,730) 'ei', 'xi'
       do 780 j = 1,rtnoil
            write (8,740) e3oil(j), x3oil(j)
780    continue
       write (8,720)
       write (8,705) 'type iii kinetic parameters: ch4'
       write (8,730) 'ei', 'xi'
       do 785 j = 1,rtnch4
            write (8,740) e3ch4(j), x3ch4(j)
785    continue
       write (8,720)
       write (8,705) 'kinetic parameters: thermal gas'
       write (8,730) 'ei', 'xi'
       do 790 j = 1,rtntgs
            write (8,740) etgs(j), xtgs(j)
790    continue
       write (8,722)
704    format (1x,//)
705    format (3x, a)
715    format (10x, a, 3(11x, a ))
720    format (1x, /)
730    format (9x ,a, 15x, a)
740    format (2 (3x, 1pe12.3 ))
722    format (1h1,
     $         //,20x,'*****************************************',
     $         /,20x,' t r a n s i e n t    s i m u l a t i o n',
     $         /,20x,'*****************************************',
     $  //)
       endif
       return
       end
!***********************************************************************
      subroutine petpot (head,heado,headg,rhof,nnode,z,rhoo,rhog
     $  ,pres)
!
!
!     this subroutine computes to oil heads
!
!***********************************************************************

      implicit real*8 (a-h,o-z)
!
! array sizes
!
      integer maxnodes,maxelems,maxbnds,maxmats,maxflts,maxwdt

!
! declaring array sizes
!
      include 'rift2d_dimensions.txt'
!


      integer nnode,n

      real*8 head(maxnodes),z(maxnodes),rhof(maxnodes),pres(maxnodes)
      real*8 heado(maxnodes),headg(maxnodes),rhoo(maxnodes)
      real*8 rhog(maxnodes)

      do 100 n=1,nnode
          pnode = pres(n) / 1.0d+06
          
        if(pnode.lt.0.0d+00)pnode=0.00d+00
          rhoo(n) = 781.64 - 7.8944*pnode +
     $               4.1642*dexp(-2*(pnode**2))
          rhog(n) = 43.444 + (10.281*pnode)
         if(rhoo(n).lt.450.d+00) rhoo(n) = 450.00d+00
         if(rhog(n).gt.400.d+00) rhog(n) = 400.00d+00

      heado(n) = rhof(n)*head(n)/rhoo(n)
     $ + z(n)*(rhoo(n)-rhof(n))/rhoo(n)
      headg(n) = rhoo(n)*heado(n)/rhog(n)
     $ + z(n)*(rhog(n)-rhoo(n))/rhog(n)
  100 continue
      return
      end
!*************************************************************************
      real*8 function atand(x)
      real*8 x
      atand = atan(x) * 57.2957795130823d0
      return
      end
!*************************************************************************
      real*8 function dcosd(x)
      real*8 x
      dcosd = dcos(x/57.2957795130823d0)
      return
      end
!*************************************************************************
      real*8 function dsind(x)
      real*8 x
      dsind = dsin(x/57.2957795130823d0)
      return
      end
!*************************************************************************
      real*8 function dtand(x)
      real*8 x
      dtand = dtan(x/57.2957795130823d0)
      return
      end
