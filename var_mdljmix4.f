************************************************************************
*    include file for the MDLJ.FOR program
*
************************************************************************
      implicit none
      real*8 tstar, pstar, plrc, upot, ulrc, rnp, xx, yy, zz, fx,
     &  fy, fz, vx, vy, vz, cutoff, cutsq, del_t, del_tby2, del_tsqby2,
     &  rho, treq, xdisp, ydisp, zdisp, vol, ukin, pshift, delr, 
     &  rdels, rmaxsq, accum, sigbox, sigsq, pi, virial, rlistsq, 
     &  rdiffsq, xldisp, yldisp, zldisp, free, sig,eps,mass,rmass,
     &  pshift6,utot,del_eps,del_sig,xmol,yij0,cpot,ryij,cpave,cpsq,
     &  yijave, yijsq,yijr,yrmax,upcs,rytest,eps11,sig11,mass11,
     &  ander

      integer npmax, np, nac, nco, nbinmax, nsteps, ngofr, nprint, 
     &  nconf, nsave, nbins, nsub, iseed, jseed, counts, hist, ncut,
     &  nlmax, list, point, nscale,ncomp,itype,ntype,ncmax,ngij,
     &  nwidom,ntp,ybins,yrcount,nangles,tconst

      logical contin, update, use_list, cut_and_shift, zero_momen

      PARAMETER (PI = 3.14159265358979323846)
      parameter (npmax = 20000, nac=40, nco=20, nbinmax=6000)
      parameter (nlmax = 1500*npmax)
      parameter (ncmax = 10)

      common /vects/ xx(npmax), yy(npmax), zz(npmax), fx(npmax), 
     &  fy(npmax), fz(npmax), vx(npmax), vy(npmax), vz(npmax),
     &  xdisp(npmax), ydisp(npmax), zdisp(npmax), accum(nac),
     &  xldisp(npmax), yldisp(npmax), zldisp(npmax),sig(ncmax,ncmax),
     &  eps(ncmax,ncmax),mass(ncmax),rmass(ncmax),xmol(ncmax),
     &  yij0(ncmax*(ncmax+1)/2),cpot(ncmax),ryij(0:nbinmax),
     &  cpave(ncmax),cpsq(ncmax),yijave(ncmax*(ncmax+1)/2),
     &  yijsq(ncmax*(ncmax+1)/2),yijr(0:nbinmax,ncmax*ncmax)


      common /param/ tstar, pstar, plrc, upot, ulrc, rnp, cutoff,
     &  cutsq, del_t, del_tby2, del_tsqby2, rho, treq, vol, ukin,
     &  delr, rdels, rmaxsq, sigbox, sigsq, virial, pshift, rlistsq,
     &  rdiffsq, free,pshift6,utot,del_eps,del_sig,yrmax,upcs,rytest,
     &  eps11,sig11,mass11,ander

      common /ivects/counts(nco),hist(nbinmax,ncmax*(ncmax+1)/2), 
     &  list(nlmax),point(npmax),itype(npmax),ntype(ncmax),
     &  yrcount(0:nbinmax,ncmax*ncmax)

      common /iparam/ np, nsteps, ngofr, nprint, nconf, nsave, nbins,
     &  nsub, iseed, jseed, ncut, nscale,ncomp,ngij,
     &  nwidom,ntp,ybins,nangles,tconst

      common /logi/contin,update,use_list,cut_and_shift,zero_momen

************************************************************************
