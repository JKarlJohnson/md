************************************************************************
*
      program mdljmix
*
*  This is a coding of the velocity Verlet algorithm for molecular
*  dynamics of Lennard-Jones atoms.  A cut and shifted LJ potential
*  used.  The potential shift is added back in and the long range
*  corrections are included at the end of the simulation.  The program
*  also does calculations using the leap-frog algorithm.  This code 
*  has been extended to mixtures.
*
************************************************************************
*                HISTORY
*
*  30 September 1991   Coding began
*  20 January 1992     Added the option of doing a leap-frog 
*                      algorithm simulation.
*  11 February 1992    Added neighbor list option
*  15 February 1992    Fixed error in neighbor list automatic update
*  17 February 1992    Added option to equilibrate and rezero counters
*   4 March 1992       Rezeros momentum if starting a new run from an
*                      old configuration.  Also scales initial velocities
*                      if temperature is constant.
*   5 March 1992       Changed the NVT code to perform simple velocity
*                      scaling at a specified interval.
*   6 Sept. 1993       Added a counter for the number of atoms within
*                      a distance of 1.0 sigma.
*  14 April 1995       Prints out mean-square displacement instead of
*                      RMS displacement.
*  26 June 1995   Made standard deviation calculations include the
*                 factor of num. subblocks for error of the mean.
*  11 Oct. 1995   Added in the Widom particle insertion technique
*                 to compute the chemical potential.
*  18 June 1996   Began modifying to do MD of LJ mixtures.
*  31 July 1996   Began coding yij(r) collection based on the 
*                 method of Lee & Shing, J.C.P. 91, 477-488 (1989).
*  2 Sept 1996    Coded alternate yij(r) collection method based
*                 on Lee & Shing Eqns (4.2) and (4.3)
*  30 July 1997   Fixed bug in init_vel, pointed out by Sivakumar.
*  17 Sept 1997   Modified code to do cut and shifted fluid calculaions.
*  18 March 1998  Completed code to collect g(r) over the entire box.
*                 Uses numerical integration rather than Walter's polynomial.
*  25 June 2002   Fixed bug in init_vel. Was it introduced on 30 July 1997?
*  27 June 2002   Zero momenta of each component separately in init_vel.
*  9 Jan. 2003    Added ability of adjusting the lattice to any integer 
*                 number of atoms by taking a larger lattice and randomly 
*                 deleting a number of atoms.
*
*  J Karl Johnson
************************************************************************
*
*                   Main variables
*
*  rho        Reduced number density.
*  np         Number of particles.
*  del_t      Time step in units of (m/epsilon)^(1/2)sigma (Changed
*             to box units in the program).
*  nsteps     Total number of time steps to make.
*  xx(i)      X-component of the position vector for the center of mass.
*  fx(i)      X-component of the force on particle i.
*  vx(i)      X-component of the velocity vector for particle i.
*  tstar      Current temperature.
*  treq       Required temperature for constant temperature runs.
*  cutoff     Potential cutoff in sigma units (changed to box units).
*  upot       Potential energy at current time (no l.r.c.).
*  upcs       Cut and shifted energy.
*  ukin       Kinetic energy per particle at current time.
*  pshift6    The 6th power shift in the potential energy for the cut and 
*             shifted potential.
*  pshift     Total potential energy shift at the current time step.
*  ncut	      Number of pairs inside the potential cutoff.
*  utot	      Total energy = ukin + upot - pshift
*  pstar      Current pressure from virial calculations.
*  ncomp      Number of components in the mixture.
*  sig(i,j)   LJ diameter between molecules of type i and j.
*  eps(i,j)   LJ well depth between molecules of type i and j.
*  itype(i)   The type of molecule i.
*  ntype(i)   Number of molecules of type i. np = sum_i ntype(i).
*  xmol(i)    Mole fraction of component i.
*  mass(i)    Mass of molecule of type i.
*  rmass(i)   Reciprocal mass of molecule of type i.
*  del_sig    Lorentz-Berthelot perturbation parameter in sigma.
*  del_eps    Lorentz-Berthelot perturbation parameter in epsilon.
*  nequil     Maximum number of equilibration steps.
*  dequil     Root mean square displacement for equilibration.
*  eps11      The epsilon for component 1 in real units.
*  sig11      The value of sigma for component 1 in real units.
*  mass11     The mass of component 1 in real units.
*
*  counts(i)  Integer counters
*  counts(1)  Total number of time steps taken thus far.
*  counts(2)  Number of calls to subroutine move in a subblock.
*  counts(3)  Number of subblocks used for error estimation.
*  counts(4)  Number of times the g(r) histograms were updated.
*  counts(5)  Number of updates to the neighbor list during the run.
*  counts(6)  Number of times the atoms within 1 sigma are counted.
*  counts(7)  Accumulator for the number of atoms within 1 sigma.
*  counts(8)  Number of times the chemical potential estimator was 
*                updated in a subblock.
*  counts(9)  Number of times the yij(0) accumulator was updated 
*                in a subblock.
*  counts(10) Number of times the yij(r) accumulator was updated
*                in a subblock.
*
*  accum(i)   Accumulators
*  accum(1)   upot in a subblock.
*  accum(2)   upot of subaverages.
*  accum(3)   upot**2 of subaverages.
*  accum(4)   ukin in a subblock.
*  accum(5)   ukin of subaverages.
*  accum(6)   ukin**2 of subaverages.
*  accum(7)   pstar in a subblock.
*  accum(8)   pstar of subaverages.
*  accum(9)   pstar**2 of subaverages.
*  accum(10)  tstar in a subblock.
*  accum(11)  tstar of subaverages.
*  accum(12)  tstar**2 of subaverages.
*  accum(13)  ncut in a subblock.
*  accum(14)  ncut of subaverages.
*  accum(15)  ncut**2 of subaverages.
*  accum(16)  utot.
*  accum(17)  utot**2.
*  accum(18)  upot**2 for fluctuation calculations.
*  accum(19)  pstar**2 for fluctuation calculations.
*  accum(20)  upcs in a subblock.
*  accum(21)  upcs of subaverages.
*  accum(22)  upcs**2 of subaverages.
*  cpot(i)    Chemial potential for component i in a subblock.
*  cpave(i)   Chemical potential for component i subaverage.
*  cpsq(i)    Chemical potential squared for component i subaverage.
*  yij0(l)    yij(0) estimator in a subblock.
*  yijave(l)  yij(0) estimator in a subaverage.
*  yijsq(l)   yij(0) squared estimator in a subaverage.
*  yijr(l,k)  yij(r) histograms.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 dequil, msd, rsig
      integer ip, nequil, ntemp, ntgofr, ntwidom
      character*1 label

      call startup(nequil, dequil)
      rsig = 1.0d0/sigbox
*  Build a neighbor list if we are using one and call the force routine
      update = use_list
      call forces
      if (nequil .gt. 0) then
         ntemp = nsteps
         nsteps = nequil
         ntgofr = ngofr
         ntwidom = nwidom
         ngofr = 0
         nwidom = 0
      endif
100   if (counts(1) .lt. nsteps) then
         counts(1) = counts(1) + 1
*  Use the velocity Verlet algorithm.
         call move_vv
*  Calculate the chemical potential if needed.
         if (nequil .eq. 0) then 
            if (nwidom .gt. 0) then 
               if (mod(counts(1),nwidom) .eq. 0) call widom(3)
            endif
         endif
         call collect_aves
*  Calculate subblock averages if end of subblock
         if (mod(counts(1), nsub) .eq. 0) call avergs
*  Print out intermediate results every nprint steps.
         if (nprint .gt. 0) then
            if (mod(counts(1), nprint) .eq. 0) call printit(msd)
            if ((nequil .gt. 0) .and. (sqrt(msd) .gt. dequil)) goto 500
         endif
         if (nconf .gt. 0) then
*  Write the current configuration to a file every nconf steps.
            if ((nequil .eq. 0) .and. 
     &          (mod(counts(1), nconf) .eq. 0)) then
               write(55) counts(1), np
               write(55) (xdisp(ip)*rsig, ydisp(ip)*rsig, 
     &                    zdisp(ip)*rsig, ip=1,np)
* Write out an xyz file in Angstrom units:
               write(44,*) np
               write(44,*) 
               do ip = 1, np
                  if (itype(ip) .eq. 1) then
                     label = 'C'
                  else
                     label = 'H'
                  endif
                  write(44,250) label, xx(ip)*sig11/sigbox,yy(ip)*
     &                 sig11/sigbox, zz(ip)*sig11/sigbox
 250              format(1x, 1a, 3(2x, e15.8))
               enddo
            endif
         endif
*  Write out the restart file every nsave steps.
         if (mod(counts(1), nsave) .eq. 0) call saveit
         goto 100
      endif
*  Test for equilibration.
500   if (nequil .gt. 0) then
         nequil = 0
         call rezero
         nsteps = ntemp
         ngofr = ntgofr
         nwidom = ntwidom
         goto 100
      endif

      call ending

      stop 
      end
************************* end of main program **************************
************************************************************************
*
      block data
*
************************************************************************
      implicit none
      logical test

      common /raset4/ test
      data test /.FALSE./

      end
********************* end of block data subroutine *********************
************************************************************************
*
      subroutine startup(nequil, dequil)
*
*  Reads in the initial configuration and the input parameters.  
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 xdum, ydum, zdum, rc3, rc9, rlist, 
     &    dequil,sigij3,trho
      integer ntest, istep, ip, jp, nequil, ntst(ncmax), ntoobig
      character title*79, latype*3
      logical newconf, newrho, adjust

      open(unit=1, file='in.mdljmix', status='old')
      open(unit=3, file='res.mdljmix', status='unknown')
* Open the file for the xyz coordinates:       
      open(unit=44, file='config.xyz', status='unknown')

      read(1,990) title
990   format(a79)
      read(1,*) np
      read(1,*) ncomp
      read(1,*) (ntype(ip),ip=1,ncomp)
      read(1,*) rho
      read(1,*) tconst
      read(1,*) zero_momen
      read(1,*) ander
      read(1,*) treq
      read(1,*) del_t
      read(1,*) cutoff
      read(1,*) nsteps
      read(1,*) nsub
      read(1,*) nprint
      read(1,*) nconf
      read(1,*) nsave
      read(1,*) ngofr
      read(1,*) delr
      read(1,*) newconf
      read(1,*) contin
      read(1,*) iseed
      read(1,*) jseed
      read(1,*) use_list
      read(1,*) rlist
      read(1,*) nequil
      read(1,*) dequil
      read(1,*) nscale
      read(1,*) ntp
      read(1,*) nwidom
      read(1,*) (eps(ip,ip),ip=1,ncomp)
      read(1,*) del_eps
      read(1,*) (sig(ip,ip),ip=1,ncomp)
      read(1,*) del_sig
      read(1,*) (mass(ip),ip=1,ncomp)
*  Make sure the parameters are reduced w.r.t. compontent 1:
      eps11 = eps(1,1)
      sig11 = sig(1,1)
      mass11 = mass(1)
      do ip = 1, ncomp
         eps(ip,ip) = eps(ip,ip)/eps11
         sig(ip,ip) = sig(ip,ip)/sig11
         mass(ip) = mass(ip)/mass11
      enddo
      if (abs(eps(1,1)-1.0d0) .gt. 1.0d-6 .or. 
     &    abs(sig(1,1)-1.0d0) .gt. 1.0d-6 .or.
     &    abs(mass(1)-1.0d0) .gt. 1.0d-6) then 
         write(*,*) '    Error:  eps(1,1) must be 1'
         write(3,*) '    Error:  eps(1,1) must be 1'
         write(*,*) '    Error:  sig(1,1) must be 1'
         write(3,*) '    Error:  sig(1,1) must be 1'
         write(*,*) '    Error:  mass(1) must be 1'
         write(3,*) '    Error:  mass(1) must be 1'
         stop
      endif
      read(1,*) yrmax
      read(1,*) nangles
      read(1,*) rytest
      read(1,*) cut_and_shift
      read(1,555) latype
 555  format(a3)
      read(1,*) newrho
      close(1)
      do ip = 1, ncomp
*  Store the reciprocal mass for calculating the velocities:
         rmass(ip) = 1.0d0/mass(ip)
      enddo

      rnp = real(np)
*  Compute the number of degrees of freedom:
      if (zero_momen) then
         free = 3.0d00*rnp - 3.0d00
      else
         free = 3.0d00*rnp
      endif
      call rmarin(iseed, jseed)
      if (newconf) then
         call lattice(np, xx, yy, zz, latype, ntoobig, adjust)
*  Adjust the number of atoms in the lattice if needed:
         if (adjust) call adjust_lattice(np, ntoobig, xx, yy, zz)
         call assign_type(np, ncomp, ntype, itype)
*  Check the type assignment.
         do ip = 1, ncomp
            ntst(ip) = 0
         enddo
         do ip = 1, np
            ntst(itype(ip)) = ntst(itype(ip)) + 1
         enddo
         do ip = 1, ncomp
            if (ntst(ip) .ne. ntype(ip)) then
               write(*,*) '   Error: ', ntst(ip), '.ne.', ntype(ip)
               write(3,*) '   Error: ', ntst(ip), '.ne.', ntype(ip)
            endif
         enddo
         call init_vel
      else
         open(unit=2, file='old.mdljmix', status='old')
*  Read in the variables for a restart or continuation
         read(2,*) trho
         if (newrho) then
            if (rho - trho .gt. 0.05) then
               write(3,*) '  Error: density in start file too high.'
               write(*,*) '  Error: density in start file too high.'
            endif
         else
            rho = trho
         endif
         read(2,*) (accum(ip), ip=1,nac)
         read(2,*) (counts(ip), ip=1,nco)
         read(2,*) ntest, ncomp
         if (ntest .ne. np) then
            write(3,*) '          >>>>>>>> ERROR <<<<<<<<'
            write(3,*) '   The number of particles specified is not ',
     &                  'the same as in the configuration file.'
            stop
         endif
*  Read in the LJ parameters
         read(2,*) (eps(ip,ip),ip=1,ncomp)
         read(2,*) del_eps
         read(2,*) (sig(ip,ip),ip=1,ncomp)
         read(2,*) del_sig
*  Read in the mass of each type:
         read(2,*) (mass(ip),ip=1,ncomp)
*  Read in the type of each molecule and compute the number of each type.
         do ip=1,ncomp
*  Store the reciprocal mass for calculating the velocities:
            rmass(ip) = 1.0d0/mass(ip)
            ntype(ip) = 0
         enddo
         read(2,*) (itype(ip), ip=1,np)
         do ip=1,np
            ntype(itype(ip)) = ntype(itype(ip))+1
         enddo
*  Read in the positions and velocities
         read(2,*) (xx(ip), yy(ip), zz(ip), vx(ip), vy(ip), vz(ip), 
     &             ip = 1, np)
         close(2)
*  Rezero the momentum if this is not a continuation run
         if (zero_momen .and. (.not. contin)) call p_zero
      endif
*  Compute the LJ cross terms using perturbed Lorentz-Berthelot mixing rules:
      do ip = 1, ncomp-1
         do jp = ip+1, ncomp
            eps(ip,jp) = sqrt(eps(ip,ip)*eps(jp,jp))*(1.0d0+del_eps)
            eps(jp,ip) = eps(ip,jp)
            sig(ip,jp) = 0.5d0*(sig(ip,ip)+sig(jp,jp))*(1.0d0+del_sig)
            sig(jp,ip) = sig(ip,jp)
         enddo
      enddo

*  Convert some variables
      vol = rnp/rho
      sigbox = (1.0d00/vol)**(1.0d00/3.0d00)
      sigsq = sigbox**2
      del_t = del_t*sigbox
      del_tsqby2 = 0.5d00*del_t**2
      del_tby2 = 0.5d00*del_t
      cutoff = cutoff*sigbox
      if (cutoff .gt. 0.5d00) cutoff = 0.5d00
      cutsq = cutoff**2
      rlist = rlist*sigbox
*  Make sure the list cutoff is less than half the box length
      if (rlist .gt. 0.5d00) rlist = 0.5d00
*  Use a list only if rlist > cutoff
      if (use_list) use_list = (rlist .gt. cutoff)
      rlistsq = rlist*rlist
      rdiffsq = (rlist - cutoff)**2
      rdels = 1.0d00/(delr*sigbox)
*  Calculate the 6th power shift in the potential.
      pshift6 = (sigbox/cutoff)**6
*  Calculate the long range corrections to the potential energy and
*  the pressure
      rc3 = (sigbox/cutoff)**3
      rc9 = rc3*rc3*rc3
*  ulrc is on a per particle basis.
      ulrc = 0.0
      plrc = 0.0
      if (.not. cut_and_shift) then
         do ip = 1, ncomp
            do jp = 1, ncomp
               sigij3 = sig(ip,jp)
               sigij3 = sigij3*sigij3*sigij3
               ulrc=ulrc+ntype(ip)*ntype(jp)*eps(ip,jp)*sigij3*
     &           sigij3*(rc9*sigij3*sigij3-3.0d0*rc3)
               plrc=plrc+ntype(ip)*ntype(jp)*eps(ip,jp)*sigij3*
     &           sigij3*(rc9*sigij3*sigij3-1.5d0*rc3)
            enddo
         enddo
         ulrc = 8.0d00*pi*rho*ulrc/(9.0d0*np*np)
         plrc = 32.0d00*pi*rho**2*plrc/(9.0d0*np*np)
      endif
*  Compute long range corrections for chemical potential and yij
      if (nwidom .gt. 0) call widom(0)
*  Compute the mole fractions:
      do ip = 1, ncomp
         xmol(ip) = float(ntype(ip))/rnp
      enddo
*  Variables for the g(r) subroutine:
*  Compute the number of distinct pair correlation functions:
      ngij = ncomp*(ncomp+1)/2
      if (ngofr .gt. 0) then
*  Collect g(r) upto sqrt(3)*half the box length, i.e., the semi-diagonal
*  of the simulation cell. Subtract off one sigma at the end.
         nbins = (sqrt(3.0)*0.5d00-sigbox)*rdels
         if (nbins .gt. nbinmax) nbins=nbinmax
         rmaxsq = (float(nbins)/rdels)**2
         do ip = 1, nbins
            do jp = 1, ngij 
               hist(ip,jp) = 0
            enddo
         enddo
      endif
*  Read in the other variables if this is a continuation run
      if (contin) then
         open(unit=2, file='msd.mdljmix', status='old', 
     &    form='unformatted')
         read(2) (xdisp(ip), ydisp(ip), zdisp(ip), ip = 1, np)
         close(2)
         if (ngofr .gt. 0) then
            open(unit=2, file='hist.mdljmix', status='unknown')
            read(2,*) nbins
            do jp=1,ngij
               read(2,*) (hist(ip,jp), ip = 1, nbins)
            enddo
            close(2)
         endif
         if (nwidom .gt. 0) then 
*  Read in the chemical potential and yij accumulators:
            open(unit=8, file='yij.cont', status='old')
            read(8,*) (cpot(ip), ip=1,ncomp)
            read(8,*) (cpave(ip), ip=1,ncomp)
            read(8,*) (cpsq(ip), ip=1,ncomp)
            read(8,*) (yij0(ip), ip=1,ngij)
            read(8,*) (yijave(ip), ip=1,ngij)
            read(8,*) (yijsq(ip), ip=1,ngij)
            read(8,*) ybins
            read(8,*) (ryij(ip), ip=1,ybins)
            read(8,*) ((yrcount(ip,jp), jp=1,ncomp*ncomp), ip=0,ybins)
            read(8,*) ((yijr(ip,jp), jp=1,ncomp*ncomp), ip=0,ybins)
            close(8)
         endif
         if (nprint .gt. 0) open(unit=66, file='prnt.mdljmix', 
     &         status='unknown')
*  I took the following off the end of the open statement
*    , access='append')
      else
         if (nprint .gt. 0) open(unit=66, file='prnt.mdljmix', 
     &         status='unknown')
      endif
      if ((newconf) .or. (.not. contin)) then
*  Zero the accumulators
         do 100 ip = 1, nac
            accum(ip) = 0.0
100      continue
         do 105 ip = 1, nco
            counts(ip) = 0
105      continue
*  Zero the displacement accumulators
         do 110 ip = 1, np
            xdisp(ip) = 0.0
            ydisp(ip) = 0.0
            zdisp(ip) = 0.0
110      continue
         if (nwidom .gt. 0) then
*  Fill the ryij(i) array
            call widom(1)
*  Zero the chemical potential and yij accumulators
            call widom(2)
         endif
      endif

*  Write out the input data to a file for checking
      write(3,900)
 900  format(15x, 'Velocity Verlet Molecular Dynamics Program',
     &     ' for Mixtures')
      write(3,1000) title, np, ncomp, rho
1000  format(5x, a79, 
     &       /, 10x, '------------ Echo of Input Data -------------',
     & /, 5x, 'Number of Atoms:',t54,i5,/,5x,'Number of components:', 
     & t54, i5, /, 5x, 'Density:',t56, f12.8)
      if (tconst .eq. 1) write(3,1010) treq, nscale
1010  format(5x, 'This is a constant temperature simulation', 
     &  /, 5x, 'Required temperature:', t54, f12.6,
     &  /, 5x, 'Velocity scaling used. ', 
     &  /, 5x, 'Number of steps between velocity scalings:', t54, i5)
      if (tconst .eq. 2) write(3,1011) treq, ander 
1011  format(5x, 'This is a constant temperature simulation', 
     &  /, 5x, 'Required temperature:', t54, f12.6,
     &  /, 5x, 'Andersen Coupling Parameter:', t54, f12.6)
      if (cut_and_shift) then 
         write(3,1015)
 1015    format(5x, 'This simulation is for a cut & shifted potential')
      endif
      write(3,1020) del_t/sigbox, cutoff/sigbox, nsteps, nsub, nprint, 
     &  nconf, nsave, ngofr, delr, iseed, jseed
1020  format(5x, 'Reduced time step:', t54, f12.6, 
     &  /, 5x, 'Cutoff value', t56, f10.6,
     &  /, 5x, 'Total number of time steps this run:', t51, i8,
     &  /, 5x, 'Number of time steps in a subaverage:', t54, i5,
     &  /, 5x, 'Number of time steps between print outs', t54, i5,
     &  /, 5x, 'Number of time steps between configuration writes', 
     &  t56, i3,
     &  /, 5x, 'Number of time steps between saves', t55, i4,
     &  /, 5x, 'Number of time steps between g(r) calls', t56, i3,
     &  /, 5x, 'Width of a g(r) shell in sigma units:', t56, f10.6,
     &  /, 5x, 'ISEED value:', t54, i5, /, 5x, 'JSEED value:', t54, i5)
      write(3,1021) ncomp
 1021 format(5x, 'Number of components:', t55, i4)
      write(3,1022) (ip,ip,eps(ip,ip),ip=1,ncomp)
 1022 format(5x, 'eps(',i2,i2,') = ', t56, f10.6)
      write(3,1023) del_eps
 1023 format(5x, 'LB modifier for epsilon:', t56, f10.6)
      write(3,1024) (ip,ip,sig(ip,ip),ip=1,ncomp)
 1024 format(5x, 'sig(',i2,i2,') = ', t56, f10.6)
      write(3,1025) del_sig
 1025 format(5x, 'LB modifier for sigma:', t56, f10.6)
      write(3,1026) (ip,ntype(ip),ip=1,ncomp)
 1026 format(5x, 'Number of molecules of type ', i2,':', t55, i4)
      write(3,1027) (ip,xmol(ip),ip=1,ncomp)
 1027 format(5x, 'Mole fraction of compontent ', i2,':', t56, f10.6)

      if (newconf) write(3,1030) tstar
*  Write out the initial temperature:
1030  format(5x, 'This simulation was started from a lattice',
     &  /, 10x, 'Initial Temperature =', 1x, f12.6)
      if (use_list) write(3,1040) rlist/sigbox
1040  format(5x, 'Neighbor List is used with r_list =', t56, f10.6)
      if (nwidom .gt. 0) write(3,1050) ntp, nwidom, yrmax, nangles, 
     &  rytest
 1050 format(5x, 'Number of test particles for Widom insertion:', 
     &   t54, i5, /, 5x, 
     &  'Number of time steps between calls to widom:', t54, i5, 
     &  /, 5x, 'Maximum distance to compute y(r):', t56, f10.6, 
     &  /, 5x, 'Number of random angles to use for y(r):', t54, i5, 
     &  /, 5x, 'Minimum distance for these angles:', t56, f10.6)
* Open the configuration file and position the file to the correct
* record
      if (nconf .gt. 0) then
         open(unit=55, file='conf.mdljmix', status='unknown', 
     &      form='unformatted')
         if ((contin) .and. (counts(1) .gt. 0)) then
200         read(55) istep, ntest
            do 210 ip = 1, np
               read(55) xdum, ydum, zdum
210         continue
            if ((istep + nconf) .le. counts(1)) goto 200
         endif
      endif

      return
      end
************************* end of subroutine startup ********************
************************************************************************
*
         subroutine rezero
*
*  Rezeros the accumulators after the system has equilibrated for
*  a specified number of steps or a specified root mean square 
*  displacement
************************************************************************
      include 'var_mdljmix4.f'
      real*8 msd
      integer ic

      call msdisp(msd)
      write(3,1000) counts(1), sqrt(msd)
1000  format(10x, 'The system has equilibrated --->'
     &  /, 15x, 'for ', i6,
     &  ' time steps and RMS displacement of ', f5.2, ' sigma')
      if (nconf .gt. 0) rewind(55)
      do 100 ic = 1, nac
         accum(ic) = 0.0
100   continue
      do 200 ic = 1, nco
         counts(ic) = 0
200   continue
*  Zero the displacement accumulators
      do 300 ic = 1, np
            xdisp(ic) = 0.0
            ydisp(ic) = 0.0
            zdisp(ic) = 0.0
300      continue

      return
      end
****************** end of subroutine rezero ***************************
************************************************************************
*
      subroutine move_vv
*
*  This procedure updates the positions, accelerations and velocities
*  and also calls the force subroutine in order to update the 
*  accelerations (forces).  The velocity Verlet algorithm is used
*  to move the particles.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 xxold, yyold,zzold,xxnew,yynew,zznew,scale,stdv,zeta,
     &       ranmar,dummy,gauss_dev
      integer ip, itp

      do 100 ip = 1, np
         itp = itype(ip)
*  Calculate the velocities at time t + 1/2 del_t:
         vx(ip) = vx(ip) + del_tby2*fx(ip)*rmass(itp)
         vy(ip) = vy(ip) + del_tby2*fy(ip)*rmass(itp)
         vz(ip) = vz(ip) + del_tby2*fz(ip)*rmass(itp)
*  Update the positions based on the current velocities and accelerations:
         xxold = xx(ip)
         yyold = yy(ip)
         zzold = zz(ip)
         xxnew = xxold + del_t*vx(ip)
         yynew = yyold + del_t*vy(ip)
         zznew = zzold + del_t*vz(ip)
*  Update the displacement vectors.  Note that this must be done before
*  periodic boundary conditions are applied.
         xdisp(ip) = xdisp(ip) + xxnew - xxold
         ydisp(ip) = ydisp(ip) + yynew - yyold
         zdisp(ip) = zdisp(ip) + zznew - zzold
*  If we are using a neighbor list then update the list displacement vector
         if (use_list) then
            xldisp(ip) = xldisp(ip) + xxnew - xxold
            yldisp(ip) = yldisp(ip) + yynew - yyold
            zldisp(ip) = zldisp(ip) + zznew - zzold
         endif
*  Apply periodic boundary conditions.  The box is of unit length.
         if (xxnew .gt. 0.5d00) then
             xxnew = xxnew - 1.0d00
         elseif (xxnew .lt. -0.5d00) then
             xxnew = xxnew + 1.0d00
         endif
         if (yynew .gt. 0.5d00) then
             yynew = yynew - 1.0d00
         elseif (yynew .lt. -0.5d00) then
             yynew = yynew + 1.0d00
         endif
         if (zznew .gt. 0.5d00) then
             zznew = zznew - 1.0d00
         elseif (zznew .lt. -0.5d00) then
             zznew = zznew + 1.0d00
         endif
*  End of periodic boundary conditions.  Now update position vectors.
         xx(ip) = xxnew
         yy(ip) = yynew
         zz(ip) = zznew
100   continue
*  Check to see if we need to build a list
      if (use_list) call check_list
*  Now update the focres (accelerations) from the velocity at the mid-point
*  Collect g(r) information into the histograms if taking g(r) data.
      if (ngofr .gt. 0) then
         if (mod(counts(1), ngofr) .eq. 0) then
            call fgr
         else
            call forces
         endif
      else
         call forces
      endif
*  Complete the velocity move to t + del_t and calculate the kinetic energy.
      ukin = 0.0
      do 200 ip = 1, np
         itp = itype(ip)
         vx(ip) = vx(ip) + del_tby2*fx(ip)*rmass(itp)
         vy(ip) = vy(ip) + del_tby2*fy(ip)*rmass(itp)
         vz(ip) = vz(ip) + del_tby2*fz(ip)*rmass(itp)
         ukin=ukin+mass(itp)*(vx(ip)*vx(ip)+vy(ip)*vy(ip)+
     &        vz(ip)*vz(ip))
200   continue
      ukin = 0.5d00*ukin/rnp
*  Calculate the temperature (based on the number of degrees of freedom).
      tstar = 2.0d00*ukin*rnp/free

*  Only for constant temperature runs: Scale velocities every nscale steps
      if ((tconst .eq. 1) .and. (mod(counts(1), nscale) .eq. 0)) then
*  Scale the velocities to get the requiried temperature.
         scale = sqrt(treq/tstar)
         ukin = 0.0
         do 210 ip = 1, np
            itp = itype(ip)
            vx(ip) = scale*vx(ip)
            vy(ip) = scale*vy(ip)
            vz(ip) = scale*vz(ip)
            ukin = ukin + mass(itp)*(vx(ip)**2 + vy(ip)**2 + 
     &           vz(ip)**2)
210      continue
         ukin = 0.5d00*ukin/rnp
*  Calculate the temperature (based on the degrees of freedom).
         tstar = 2.0d00*ukin*rnp/free
*  End of temperature scaling
      elseif (tconst .eq. 2) then 
*  This is the Andersen thermostat:
         ukin = 0.0
         stdv = sqrt(treq)
         do ip = 1, np
            itp = itype(ip)
            zeta = ranmar()
            if (zeta .lt. ander*del_t/sigbox) then
* Pick a new velocity for this atom from the Maxwell-Boltzmann distribution:
               vx(ip) = stdv*gauss_dev(dummy)/sqrt(mass(itp))
               vy(ip) = stdv*gauss_dev(dummy)/sqrt(mass(itp))
               vz(ip) = stdv*gauss_dev(dummy)/sqrt(mass(itp))
            endif
            ukin = ukin + mass(itp)*(vx(ip)**2 + vy(ip)**2 +
     &           vz(ip)**2)
         enddo
         ukin = 0.5d00*ukin/rnp
         tstar = 2.0d00*ukin*rnp/free
*   End of Andersen thermostat section.
      endif
*  Using the temperature calculate the pressure.
      pstar = rho*tstar + virial/vol

      return
      end
****************** end of subroutine move_vv  **************************
************************************************************************
*
      subroutine collect_aves
*
*  Collects the running averages
*
************************************************************************
      include 'var_mdljmix4.f'

*  Running totals.
      accum(1) = accum(1) + upot
      accum(18) = accum(18) + upot*upot
      accum(4) = accum(4) + ukin
      accum(7) = accum(7) + pstar
      accum(19) = accum(19) + pstar*pstar
      accum(10) = accum(10) + tstar
      accum(13) = accum(13) + ncut
      counts(2) = counts(2) + 1
*  Quantities for energy fluctuations.
      utot = upot - pshift/rnp + ukin
      accum(16) = accum(16) + utot
      accum(17) = accum(17) + utot*utot
      accum(20) = accum(20) + upcs

      return
      end
*********  end of subroutine collect_aves ******************************
************************************************************************
*
      subroutine check_list
*
*  This procedure checks whether the neighbor list needs to be
*  reconstructed or not.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 dispmax, tempx, tempy, tempz, tmax
      integer ip

*  Estimate the maximum displacement since the last update
      dispmax = 0.0
      do 100 ip = 1, np
         tempx = xldisp(ip)
         tempx = tempx*tempx
         tempy = yldisp(ip)
         tempy = tempy*tempy
         tempz = zldisp(ip)
         tempz = tempz*tempz
         tmax = tempx + tempy + tempz
         if (tmax .gt. dispmax) dispmax = tmax
100   continue
*  Use 2x the displacement of one molecule as a conservative estimate
*  of the relative displacement of a pair of molecules.
         dispmax = 4.0*dispmax
         update = (dispmax .gt. rdiffsq)

         return
         end

*************** end of subroutine check_list ***************************
************************************************************************
      subroutine forces
*
*  This procedure calculates the forces (accelerations) for each 
*  molecule and the energy and pressure for the system.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 xxi, yyi, zzi, fxi, fyi, fzi, fxij, fyij, fzij, xij, yij,
     &       zij, rij2, r_rij2, uij, virij, r6, r12, fij,sigij6,epsij
      integer ip, jp, jbeg, jend, jnabor, nlist,itp,jtp

*  Zero forces:
      do 100 ip = 1, np
         fx(ip) = 0.0
         fy(ip) = 0.0
         fz(ip) = 0.0
100   continue
*  Zero other variables:
      ncut = 0
      pshift = 0.0
      upot = 0.0
      virial = 0.0
*  Test to see if we are updating the neighbor list
      if (.not. use_list) goto 1000
      if (update) then
         update = .false.
*  Increment the number of updates:
         counts(5) = counts(5) + 1
*  Rezero displacement vector for the neighbor list.
         do 150 ip = 1, np
            xldisp(ip) = 0.0
            yldisp(ip) = 0.0
            zldisp(ip) = 0.0
150      continue
         nlist = 0
*  Begin the outer loop
         do 300 ip = 1, np - 1
            point(ip) = nlist + 1
            xxi = xx(ip)
            yyi = yy(ip)
            zzi = zz(ip)
            fxi = fx(ip)
            fyi = fy(ip)
            fzi = fz(ip)
            itp = itype(ip)
*  Begin the inner loop
            do 200 jp = ip + 1, np
               xij = xxi - xx(jp)
               yij = yyi - yy(jp)
               zij = zzi - zz(jp)
               jtp = itype(jp)
*  Apply the minimum image convention
               if (xij .lt. -0.5d00) then
                  xij = xij + 1.0d00
               elseif (xij .gt. 0.5d00) then
                  xij = xij - 1.0d00
               endif
               if (yij .lt. -0.5d00) then
                  yij = yij + 1.0d00
               elseif (yij .gt. 0.5d00) then
                  yij = yij - 1.0d00
               endif
               if (zij .lt. -0.5d00) then
                  zij = zij + 1.0d00
               elseif (zij .gt. 0.5d00) then
                  zij = zij - 1.0d00
               endif
               rij2 = xij*xij + yij*yij + zij*zij
*  End of minimum image convention
               if (rij2 .lt. rlistsq) then
                  nlist = nlist + 1
                  if (nlist .gt. nlmax) then
                     write(3,*) '  Neighbor list too small'
                     stop
                  endif
                  list(nlist) = jp
               endif
*  Check to see if the distance is within the cutoff
               if (rij2 .gt. cutsq) goto 200
                  sigij6=sig(itp,jtp)
                  sigij6=sigij6*sigij6*sigij6
                  sigij6=sigij6*sigij6
                  epsij = eps(itp,jtp)
                  r_rij2 = 1.0d00/rij2
                  r6 = sigsq*r_rij2
                  r6 = r6*r6*r6*sigij6
                  r12 = r6*r6
                  uij = epsij*(r12 - r6)
                  upot = upot + uij
                  virij = uij + epsij*r12
                  virial = virial + virij
*  Calculate the forces:
                  fij = virij*r_rij2
                  fxij = fij*xij
                  fyij = fij*yij
                  fzij = fij*zij
*  Force on particle i
                  fxi = fxi + fxij
                  fyi = fyi + fyij
                  fzi = fzi + fzij
*  Force on particle j
                  fx(jp) = fx(jp) - fxij
                  fy(jp) = fy(jp) - fyij
                  fz(jp) = fz(jp) - fzij
*  Update the number of pairs within the cutoff
                  ncut = ncut + 1
*  Compute the potential shift:
                  pshift=pshift+epsij*pshift6*sigij6*(pshift6*
     &                 sigij6-1.0d0)
*  End of inner loop
200         continue  
*  Update the force vector for particle ip
            fx(ip) = fxi
            fy(ip) = fyi
            fz(ip) = fzi
*  End of outter loop
300      continue
         point(np) = nlist + 1
      else
*  Use the neighbor lists
*  Begin the outer loop
         do 500 ip = 1, np - 1
            jbeg = point(ip)
            jend = point(ip+1) - 1
*  Make sure that molecule ip has neighbors
            if (jbeg .le. jend) then
               xxi = xx(ip)
               yyi = yy(ip)
               zzi = zz(ip)
               fxi = fx(ip)
               fyi = fy(ip)
               fzi = fz(ip)
               itp = itype(ip)
*  Begin the inner loop
               do 400 jnabor = jbeg, jend
                  jp = list(jnabor)
                  xij = xxi - xx(jp)
                  yij = yyi - yy(jp)
                  zij = zzi - zz(jp)
                  jtp = itype(jp)
*  Apply the minimum image convention
                  if (xij .lt. -0.5d00) then
                     xij = xij + 1.0d00
                  elseif (xij .gt. 0.5d00) then
                     xij = xij - 1.0d00
                  endif
                  if (yij .lt. -0.5d00) then
                     yij = yij + 1.0d00
                  elseif (yij .gt. 0.5d00) then
                     yij = yij - 1.0d00
                  endif
                  if (zij .lt. -0.5d00) then
                     zij = zij + 1.0d00
                  elseif (zij .gt. 0.5d00) then
                     zij = zij - 1.0d00
                  endif
                  rij2 = xij*xij + yij*yij + zij*zij
*  End of minimum image convention
*  Check to see if the distance is within the cutoff
                  if (rij2 .gt. cutsq) goto 400
                     sigij6=sig(itp,jtp)
                     sigij6=sigij6*sigij6*sigij6
                     sigij6=sigij6*sigij6
                     epsij = eps(itp,jtp)
                     r_rij2 = 1.0d00/rij2
                     r6 = sigsq*r_rij2
                     r6 = r6*r6*r6*sigij6
                     r12 = r6*r6
                     uij = epsij*(r12 - r6)
                     upot = upot + uij
                     virij = uij + epsij*r12
                     virial = virial + virij
*  Calculate the forces:
                     fij = virij*r_rij2
                     fxij = fij*xij
                     fyij = fij*yij
                     fzij = fij*zij
*  Force on particle i
                     fxi = fxi + fxij
                     fyi = fyi + fyij
                     fzi = fzi + fzij
*  Force on particle j
                     fx(jp) = fx(jp) - fxij
                     fy(jp) = fy(jp) - fyij
                     fz(jp) = fz(jp) - fzij
*  Update the number of pairs within the cutoff
                     ncut = ncut + 1
*  Compute the potential shift:
                     pshift=pshift+epsij*pshift6*sigij6*(pshift6*
     &                 sigij6-1.0d0)
*  End of inner loop
400            continue  
*  Update the force vector for particle i
               fx(ip) = fxi
               fy(ip) = fyi
               fz(ip) = fzi
*  End of if (jbeg .le. jend) block
            endif
*  End of outter loop
500      continue
*  End of if block for neighbor list
      endif
*  Skip the next block--it is for calculations without a neighbor list
      goto 1500
*  This section is for no neighbor list
1000  continue
*  Begin the outer loop
      do 1200 ip = 1, np - 1
         xxi = xx(ip)
         yyi = yy(ip)
         zzi = zz(ip)
         fxi = fx(ip)
         fyi = fy(ip)
         fzi = fz(ip)
         itp = itype(ip)
*  Begin the inner loop
         do 1100 jp = ip + 1, np
            xij = xxi - xx(jp)
            yij = yyi - yy(jp)
            zij = zzi - zz(jp)
            jtp = itype(jp)
*  Apply the minimum image convention
            if (xij .lt. -0.5d00) then
               xij = xij + 1.0d00
            elseif (xij .gt. 0.5d00) then
               xij = xij - 1.0d00
            endif
            if (yij .lt. -0.5d00) then
               yij = yij + 1.0d00
            elseif (yij .gt. 0.5d00) then
               yij = yij - 1.0d00
            endif
            if (zij .lt. -0.5d00) then
               zij = zij + 1.0d00
            elseif (zij .gt. 0.5d00) then
               zij = zij - 1.0d00
            endif
            rij2 = xij*xij + yij*yij + zij*zij
*  End of minimum image convention
*  Check to see if the distance is within the cutoff
            if (rij2 .lt. cutsq) then
               sigij6=sig(itp,jtp)
               sigij6=sigij6*sigij6*sigij6
               sigij6=sigij6*sigij6
               epsij = eps(itp,jtp)
               r_rij2 = 1.0d00/rij2
               r6 = sigsq*r_rij2
               r6 = r6*r6*r6*sigij6
               r12 = r6*r6
               uij = epsij*(r12 - r6)
               upot = upot + uij
               virij = uij + epsij*r12
               virial = virial + virij
*  Calculate the forces:
               fij = virij*r_rij2
               fxij = fij*xij
               fyij = fij*yij
               fzij = fij*zij
*  Force on particle i
               fxi = fxi + fxij
               fyi = fyi + fyij
               fzi = fzi + fzij
*  Force on particle j
               fx(jp) = fx(jp) - fxij
               fy(jp) = fy(jp) - fyij
               fz(jp) = fz(jp) - fzij
*  Update the number of pairs within the cutoff
               ncut = ncut + 1
*  Compute the potential shift:
               pshift=pshift+epsij*pshift6*sigij6*(pshift6*
     &              sigij6-1.0d0)
            endif
*  End of inner loop
1100     continue  
*  Update the force vector for particle i
         fx(ip) = fxi
         fy(ip) = fyi
         fz(ip) = fzi
*  End of outter loop
 1200 continue

*  Transfer to this point from the loops using a neighbor list
 1500 continue
*  Allow for a factor of 4 in the energy and 8 in the virial (24/3).
      upot = 4.0d00*upot/rnp
      virial = 8.0d00*virial
*  Allow for the factor of 24 in the forces.
      do 2000 ip = 1, np
         fx(ip) = 24.0d00*fx(ip)
         fy(ip) = 24.0d00*fy(ip)
         fz(ip) = 24.0d00*fz(ip)
 2000 continue
      pshift = 4.0*pshift
      upcs = upot - pshift/rnp

      return
      end
******************** end of subroutine forces **************************
************************************************************************
*
      subroutine fgr
*
*  This procedure calculates the forces (accelerations) for each 
*  molecule and the energy and pressure for the system, and also
*  updates the histograms for the g(r) calculations.  This subroutine
*  is slower than the forces subroutine and so it is only called
*  when the g(r) histograms need to be updated.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 xxi, yyi, zzi, fxi, fyi, fzi, fxij, fyij, fzij, xij, yij,
     &       zij,rij2,r_rij2,uij,virij,r6,r12,fij,rij,sigij6,epsij
      integer ip, jp, ibin, nlist,itp,jtp,kp

*  Increment the counter for the number of atoms within 1 sigma:
      counts(6) = counts(6) + 1
*  Zero forces:
      do 100 ip = 1, np
         fx(ip) = 0.0
         fy(ip) = 0.0
         fz(ip) = 0.0
100   continue
*  Zero other variables:
      ncut = 0
      pshift = 0.0
      upot = 0.0
      virial = 0.0
*  Update the neighbor list if we are using one
      if (use_list) then
         update = .false.
*  Increment the number of updates:
         counts(5) = counts(5) + 1
*  Rezero displacement vector for the neighbor list.
         do 150 ip = 1, np
            xldisp(ip) = 0.0
            yldisp(ip) = 0.0
            zldisp(ip) = 0.0
150      continue
         nlist = 0
      endif
*  Begin the outer loop
      do 300 ip = 1, np - 1
         point(ip) = nlist + 1
         xxi = xx(ip)
         yyi = yy(ip)
         zzi = zz(ip)
         fxi = fx(ip)
         fyi = fy(ip)
         fzi = fz(ip)
         itp = itype(ip)
*  Begin the inner loop
         do 200 jp = ip + 1, np
            xij = xxi - xx(jp)
            yij = yyi - yy(jp)
            zij = zzi - zz(jp)
            jtp = itype(jp)
*  Apply the minimum image convention
            if (xij .lt. -0.5d00) then
               xij = xij + 1.0d00
            elseif (xij .gt. 0.5d00) then
               xij = xij - 1.0d00
            endif
            if (yij .lt. -0.5d00) then
               yij = yij + 1.0d00
            elseif (yij .gt. 0.5d00) then
               yij = yij - 1.0d00
            endif
            if (zij .lt. -0.5d00) then
               zij = zij + 1.0d00
            elseif (zij .gt. 0.5d00) then
               zij = zij - 1.0d00
            endif
*  End of minimum image convention
            rij2 = xij*xij + yij*yij + zij*zij
            if ((use_list) .and. (rij2 .lt. rlistsq)) then
               nlist = nlist + 1
                  if (nlist .gt. nlmax) then
                     write(3,*) '  Neighbor list too small'
                     stop
                  endif
               list(nlist) = jp
            endif
*  Check to see if the distance is within the cutoff range 
*  for g(r) calculations
            if (rij2 .lt. rmaxsq) then
*  See if the pair are within 1.0 sigma of each other
               sigij6=sig(itp,jtp)
               sigij6=sigij6*sigij6
               if (rij2 .lt. sigij6*sigsq) 
     &              counts(7) = counts(7) + 2
               rij = sqrt(rij2)
*  Compute which bin to store the ij pair in.
               ibin = 1 + int(rij*rdels)
*  Compute which histogram to put this in:
               if (itp .le. jtp) then
                   kp = ncomp*(itp-1)-itp*(itp+1)/2+itp+jtp
               else
                   kp = ncomp*(jtp-1)-jtp*(jtp+1)/2+itp+jtp
               endif
               if (itp .eq. jtp) then
                  hist(ibin,kp) = hist(ibin,kp) + 2
               else
                  hist(ibin,kp) = hist(ibin,kp) + 1
               endif
*  Check to see if the ij pair is within the cutoff for force calculations.
               if (rij2 .lt. cutsq) then 
*  Note that sigij6 is computed above as sig(itp,jtp)**2
                  sigij6=sigij6*sigij6*sigij6
                  epsij = eps(itp,jtp)
                  r_rij2 = 1.0d00/rij2
                  r6 = sigsq*r_rij2
                  r6 = r6*r6*r6*sigij6
                  r12 = r6*r6
                  uij = epsij*(r12 - r6)
                  upot = upot + uij
                  virij = uij + epsij*r12
                  virial = virial + virij
*  Calculate the forces:
                  fij = virij*r_rij2
                  fxij = fij*xij
                  fyij = fij*yij
                  fzij = fij*zij
*  Force on particle i
                  fxi = fxi + fxij
                  fyi = fyi + fyij
                  fzi = fzi + fzij
*     Force on particle j
                  fx(jp) = fx(jp) - fxij
                  fy(jp) = fy(jp) - fyij
                  fz(jp) = fz(jp) - fzij
*  Update the number of pairs within the cutoff
                  ncut = ncut + 1
*  Compute the potential shift:
                  pshift=pshift+epsij*pshift6*sigij6*(pshift6*
     &                 sigij6-1.0d0)
               endif
*  End of inner loop
            endif
 200     continue  
*  Update the force vector for particle i
         fx(ip) = fxi
         fy(ip) = fyi
         fz(ip) = fzi
*  End of outter loop
300   continue
*  Finish the neighbor list update
      if (use_list) then
         point(np) = nlist + 1
      endif
*  Allow for a factor of 4 in the energy and 8 in the virial (24/3).
      upot = 4.0d00*upot/rnp
      virial = 8.0d00*virial
*  Allow for the factor of 24 in the forces.
      do 400 ip = 1, np
         fx(ip) = 24.0d00*fx(ip)
         fy(ip) = 24.0d00*fy(ip)
         fz(ip) = 24.0d00*fz(ip)
400   continue
      pshift = 4.0*pshift
      upcs = upot - pshift/rnp

*  Increment the counter.  counts(4) holds the number of times that
*  the g(r) histograms were updated.
      counts(4) = counts(4) + 1

      return
      end
******************** end of subroutine fgr *****************************
************************************************************************
*
      subroutine printit(msd)
*
*  Calculates and prints out intermediate results.
*
************************************************************************
      include 'var_mdljmix4.f'
*  Local variables
      real*8 msd

*  Calculate the RMS displacement.
      call msdisp(msd)
      write(66,1000) counts(1), upot, ukin, utot, pstar, tstar, msd
      write(*,1000) counts(1), upot, ukin, utot, pstar, tstar, msd
1000  format(i8, 6(1x, e14.7))

      return
      end
********************** end of printit **********************************
************************************************************************
*
      subroutine saveit
*
*  Writes out the current configuration to a file.
*  
************************************************************************
      include 'var_mdljmix4.f'
*  Local variables
      integer ip,jp

*  Write out the new configuration and related variables.
      open(unit=4, file='new.mdljmix', status='unknown')
      write(4,*) rho
      write(4,*) (accum(ip), ip=1,nac)
      write(4,*) (counts(ip), ip=1,nco)
      write(4,*) np, ncomp
      write(4,*) (eps(ip,ip),ip=1,ncomp)
      write(4,*) del_eps
      write(4,*) (sig(ip,ip),ip=1,ncomp)
      write(4,*) del_sig
      write(4,*) (mass(ip),ip=1,ncomp)
      write(4,*) (itype(ip), ip=1,np)
      write(4,*) (xx(ip), yy(ip), zz(ip), vx(ip), vy(ip), vz(ip), 
     &             ip = 1, np)
      close(4)
      if (ngofr .gt. 0) then
*  Write out the histogram data for a restart
         open(unit=2, file='hist.mdljmix', status='unknown')
         write(2,*) nbins
         do jp=1,ngij
            write(2,*) (hist(ip,jp), ip = 1, nbins)
         enddo
         close(2)
      endif
*  Write out the mean square displacement vectors.
      open(unit=2, file='msd.mdljmix', status='unknown', 
     &  form='unformatted')
      write(2) (xdisp(ip), ydisp(ip), zdisp(ip), ip=1,np)
      close(2)
*  Write out the chemical potential and yij accumulators:
      if (nwidom .gt. 0) then
         open(unit=8, file='yij.cont', status='unknown')
         write(8,*) (cpot(ip), ip=1,ncomp)
         write(8,*) (cpave(ip), ip=1,ncomp)
         write(8,*) (cpsq(ip), ip=1,ncomp)
         write(8,*) (yij0(ip), ip=1,ngij)
         write(8,*) (yijave(ip), ip=1,ngij)
         write(8,*) (yijsq(ip), ip=1,ngij)
         write(8,*) ybins
         write(8,*) (ryij(ip), ip=1,ybins)
         write(8,*) ((yrcount(ip,jp), jp=1,ncomp*ncomp), ip=0,ybins)
         write(8,*) ((yijr(ip,jp), jp=1,ncomp*ncomp), ip=0,ybins)
         close(8)
      endif

      return
      end
******************* end of subroutine saveit ***************************
************************************************************************
*
      subroutine msdisp(msd)
*
*  Calculates the current mean-square displacement.
*
************************************************************************
      include 'var_mdljmix4.f'
*  Local variables
      real*8 msd
      integer i
      msd = 0.0
      do 10 i = 1, np
         msd = msd + xdisp(i)**2 + ydisp(i)**2 + zdisp(i)**2
10    continue
      msd = msd/(rnp*sigsq)

      return
      end
************************ end of subroutine msdisp **********************
************************************************************************
*
      subroutine avergs
*
*  Procedure to collect averages used for error estimation.
*
************************************************************************
      include 'var_mdljmix4.f'

*  Local variables
      real*8 upave, ukave, pave, tave, ncave, rmoves, cpt,upcsave
      integer ix, jx, kx

*  Increment the number of subblocks used for error estimation.
      counts(3) = counts(3) + 1

      rmoves = 1.0d00/float(counts(2))
*  Calculate the average properties for this subblock.
      upave = accum(1)*rmoves
      ukave = accum(4)*rmoves
      pave = accum(7)*rmoves
      tave = accum(10)*rmoves
      ncave = accum(13)*rmoves
      upcsave = accum(20)*rmoves
      if (nwidom .gt. 0) then
*  Compute the chemical potential and yij terms:
         kx = 0
         do ix = 1, ncomp
            cpt = cpot(ix)/float(counts(8))
            cpave(ix) = cpave(ix) + cpt
            cpsq(ix) = cpsq(ix) + cpt*cpt
            do jx = ix, ncomp
               kx = kx + 1
               cpt = yij0(kx)/float(counts(9))
               yijave(kx)=yijave(kx)+cpt
               yijsq(kx)=yijsq(kx)+cpt*cpt
            enddo
         enddo
      endif

*  Add the averages to the accumulators
      accum(2) = accum(2) + upave
      accum(3) = accum(3) + upave*upave
      accum(5) = accum(5) + ukave
      accum(6) = accum(6) + ukave*ukave
      accum(8) = accum(8) + pave
      accum(9) = accum(9) + pave*pave
      accum(11) = accum(11) + tave
      accum(12) = accum(12) + tave*tave
      accum(14) = accum(14) + ncave
      accum(15) = accum(15) + ncave*ncave
      accum(21) = accum(21) + upcsave
      accum(22) = accum(22) + upcsave*upcsave
*  Rezero all the counters and accumulators for the next subblock.
      accum(1) = 0.0
      accum(4) = 0.0
      accum(7) = 0.0
      accum(10) = 0.0
      accum(13) = 0.0
      accum(20) = 0.0
      kx = 0
      do ix = 1, ncomp
         cpot(ix) = 0.0
         do jx = ix, ncomp
            kx = kx + 1
            yij0(kx) = 0.0
         enddo
      enddo
      counts(2) = 0
      counts(8) = 0
      counts(9) = 0

      return
      end
******************** end of subroutine avergs **************************
************************************************************************
*
      subroutine ending
*
*  Prints out final results and computes final averages.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 upave, upvar, ukave, ukvar, psave, psvar, tsave, tsvar,
     &  ncave, ncvar, msd, sub1, rsub1, rsub2, utave, utvar, upfluc,
     &  psfluc, upcsave, upcsvar
      integer ip, jp

*  Write out the new configuration and related variables.
      open(unit=4, file='new.mdljmix', status='unknown')
      write(4,*) rho
      write(4,*) (accum(ip), ip=1,nac)
      write(4,*) (counts(ip), ip=1,nco)
      write(4,*) np, ncomp
      write(4,*) (eps(ip,ip),ip=1,ncomp)
      write(4,*) del_eps
      write(4,*) (sig(ip,ip),ip=1,ncomp)
      write(4,*) del_sig
      write(4,*) (mass(ip),ip=1,ncomp)
      write(4,*) (itype(ip), ip=1,np)
      write(4,*) (xx(ip), yy(ip), zz(ip), vx(ip), vy(ip), vz(ip), 
     &             ip = 1, np)
      close(4)
      if (ngofr .gt. 0) then
*  Write out the histogram data for a restart
         open(unit=2, file='hist.mdljmix', status='unknown')
         write(2,*) nbins
         do jp=1,ngij
            write(2,*) (hist(ip,jp), ip = 1, nbins)
         enddo
         close(2)
      endif
*  Write out the mean square displacement vectors.
      open(unit=2, file='msd.mdljmix', status='unknown', 
     &  form='unformatted')
      write(2) (xdisp(ip), ydisp(ip), zdisp(ip), ip=1,np)
      close(2)
*  Write out the chemical potential and yij accumulators:
      if (nwidom .gt. 0) then
         open(unit=8, file='yij.cont', status='unknown')
         write(8,*) (cpot(ip), ip=1,ncomp)
         write(8,*) (cpave(ip), ip=1,ncomp)
         write(8,*) (cpsq(ip), ip=1,ncomp)
         write(8,*) (yij0(ip), ip=1,ngij)
         write(8,*) (yijave(ip), ip=1,ngij)
         write(8,*) (yijsq(ip), ip=1,ngij)
         write(8,*) ybins
         write(8,*) (ryij(ip), ip=1,ybins)
         write(8,*) ((yrcount(ip,jp), jp=1,ncomp*ncomp), ip=0,ybins)
         write(8,*) ((yijr(ip,jp), jp=1,ncomp*ncomp), ip=0,ybins)
         close(8)
      endif

      write(3,1000)
1000  format(/, 10x, '+++++++++++++++ Results +++++++++++++++')
*  Calculate the final RMS displacement
      call msdisp(msd)
*  Calculate averages and sample standard deviations.
      sub1 = float(counts(3))
      rsub1 = 1.0d00/sub1
      rsub2 = rsub1/float(counts(3) - 1) 
      upave = accum(2)*rsub1
      upvar = (accum(3) - sub1*upave**2)*rsub2
      if (upvar .gt. 0.0) upvar = sqrt(upvar)
      upfluc = accum(18)/real(counts(1)) - upave*upave
      ukave = accum(5)*rsub1
      ukvar = (accum(6) - sub1*ukave**2)*rsub2
      if (ukvar .gt. 0.0) ukvar = sqrt(ukvar)
      psave = accum(8)*rsub1
      psvar = (accum(9) - sub1*psave**2)*rsub2
      if (psvar .gt. 0.0) psvar = sqrt(psvar)
      psfluc = accum(19)/real(counts(1)) - psave*psave
      tsave = accum(11)*rsub1
      tsvar = (accum(12) - sub1*tsave**2)*rsub2
      if (tsvar .gt. 0.0) tsvar = sqrt(tsvar)
      ncave = accum(14)*rsub1
      ncvar = (accum(15) - sub1*ncave**2)*rsub2
      if (ncvar .gt. 0.0) ncvar = sqrt(ncvar)
      utave = accum(16)/real(counts(1))
      utvar = accum(17)/real(counts(1)) - utave*utave
      upcsave = accum(21)*rsub1
      upcsvar = (accum(22) - sub1*upcsave**2)*rsub2
      if (upcsvar .gt. 0) upcsvar = sqrt(upcsvar)
* Write out the radial distribution function:
      if (ngofr .gt. 0) call a_gofr(psave+plrc, upave+ulrc)
      if (nwidom .gt. 0) then
*  Compute the chemical potential and yij terms
         call cavity(sub1,rsub1,rsub2)
      endif

      if (utvar .gt. 0.0) utvar = sqrt(utvar)
      write(3,1010) counts(1), sqrt(msd), counts(4), counts(3)
1010  format(5x, 'Total number of time steps:', t49, i10,
     &  /, 5x, 'RMS displacement:', t57, e13.6,
     &  /, 5x, 'Number of times g(r) histograms updated:', t50, i9,
     &  /, 5x, 'Number of subblocks for error estimation:', t50, i9)

      write(3,1020) upave+ulrc, upvar, upfluc, ulrc, ukave, ukvar, 
     &     psave+plrc, psvar, psfluc, plrc, tsave, tsvar, ncave, 
     &     ncvar, utvar
1020  format(5x, 'Total potential energy (including l.r.c.):', t57, 
     &   e13.6, /, 10x, 'Sample standard deviation:', t57, e13.6,
     &   /, 10x, 'Potential energy fluctuations', t57, e13.6,
     &   /, 5x, 'Potential energy Long range corrections:', t57, e13.6,
     &   /, 5x, 'Average kinetic energy:', t57, e13.6,
     &   /, 10x, 'Sample standard deviation:', t57, e13.6,
     &   /, 5x, 'Average pressure (including l.r.c.):', t57, e13.6,
     &   /, 10x, 'Sample standard deviation:', t57, e13.6,
     &   /, 10x, 'Pressure fluctuations', t57, e13.6,
     &   /, 5x, 'Pressure Long range corrections:', t57, e13.6,
     &   /, 5x, 'Average temperature:', t57, e13.6, 
     &   /, 10x, 'Sample standard deviation:', t57, e13.6,
     &   /, 5x, 'Average number of pairs within the cutoff:', t57,
     &   e13.6, /, 10x, 'Sample standard deviation:', t57, e13.6,
     &   /, 5x, 'RMS deviations in the total energy:', t57, e13.6) 
      if (cut_and_shift) then
         write(3,1023) upcsave, upcsvar
 1023    format(5x, 'Cut and shifted potential:', t57, e13.6, 
     &   /, 10x, 'Sample standard deviation:', t57, e13.6)
      endif
      if (nwidom .gt. 0) then
         write(3,1027) (ip, -tstar*log(cpave(ip)), tstar*sqrt(cpsq(ip))/
     &        cpave(ip), ip=1,ncomp) 
      endif
 1027 format(5x, 'Chemical potential for component', i2, t57, e13.6, 
     &   /, 10x, 'Sample standard deviation:', t57, e13.6)
      if (use_list) write(3,1030) counts(5), float(counts(1) - 
     &   counts(5))/float(counts(5))
1030  format(5x, 'Number of list updates:', t50, i9,
     & /, 5x, 'Average No. of time steps between updates:', t57, e13.6)

      close(3)
      close(55)
      close(6)
* Close the file holding the xyz coordinates 
      close(44)

      return
      end
******************* end of subroutine ending ***************************
************************************************************************
*
      subroutine p_zero
*
*  Rezeros the momentum if we start a new run from an old configuration
*  Also scales the initial temperature if we are doing a constant 
*  temperature run.
*
*  Scale the momendum for each component separately, as this is required 
*  for mixtures that have very different masses.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 px(ncmax), py(ncmax), pz(ncmax), scale
      integer ip

      do ip = 1, ncomp
         px(ip) = 0.0
         py(ip) = 0.0
         pz(ip) = 0.0
      enddo
      do 100 ip = 1, np
         px(itype(ip)) = px(itype(ip)) + vx(ip)*mass(itype(ip))
         py(itype(ip)) = py(itype(ip)) + vy(ip)*mass(itype(ip))
         pz(itype(ip)) = pz(itype(ip)) + vz(ip)*mass(itype(ip))
100   continue
      ukin = 0.0
      do 200 ip = 1, np
         vx(ip) = vx(ip) - px(itype(ip))*rmass(itype(ip))/
     &        ntype(itype(ip))
         vy(ip) = vy(ip) - py(itype(ip))*rmass(itype(ip))/
     &        ntype(itype(ip))
         vz(ip) = vz(ip) - pz(itype(ip))*rmass(itype(ip))/
     &        ntype(itype(ip))
         ukin = ukin + mass(itype(ip))*(vx(ip)**2 + vy(ip)**2 + 
     &        vz(ip)**2)
200   continue
      ukin = 0.5d00*ukin/rnp
      tstar = 2.0d00*ukin*rnp/free
*  Scale the initial velocities if this is a constant temperature run
      if (tconst .eq. 1) then
         scale = sqrt(treq/tstar)
         do 300 ip = 1, np
            vx(ip) = scale*vx(ip)
            vy(ip) = scale*vy(ip)
            vz(ip) = scale*vz(ip)
300      continue
         tstar = treq
         ukin = treq*free/(2.0d0*rnp)
      endif

      return
      end
******************* end of subroutine p_zero ***************************
************************************************************************
*
      subroutine init_vel
*
*  This procedure picks initial velocities from a Gaussian distribution
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 stdv, gauss_dev, px(ncmax), py(ncmax), pz(ncmax), dummy, 
     &       tcomp(ncmax)
      integer ip
*  Calculate the standard deviation of the distribution you wish to
*  sample from.  This depends on the temperature: sigma = sqrt(kT/m).
      stdv = sqrt(treq)

      do ip = 1, ncomp
         tcomp(ip) = 0.0
         px(ip) = 0.0
         py(ip) = 0.0
         pz(ip) = 0.0
      enddo
      do 100 ip = 1, np
*  Note that it is the momenta that have a Maxwell-Boltzmann distribution, 
*  not the velocities.  Hence, v_i = p_i/m_i
         vx(ip) = stdv*gauss_dev(dummy)/sqrt(mass(itype(ip)))
         vy(ip) = stdv*gauss_dev(dummy)/sqrt(mass(itype(ip)))
         vz(ip) = stdv*gauss_dev(dummy)/sqrt(mass(itype(ip)))
100   continue
*  Make sure that the momentum is zero to begin with.
      if (zero_momen) call p_zero
      ukin = 0.0
      do ip = 1, np 
         tcomp(itype(ip)) = tcomp(itype(ip)) +
     &           (vx(ip)**2 + vy(ip)**2 + vz(ip)**2)*mass(itype(ip))
         ukin=ukin+mass(itype(ip))*(vx(ip)*vx(ip)+vy(ip)*vy(ip)+
     &        vz(ip)*vz(ip))
      enddo
      ukin = 0.5d00*ukin/rnp
*  Calculate the temperature (based on the number of degrees of freedom).
      tstar = 2.0d00*ukin*rnp/free
      do ip = 1, ncomp
         if (zero_momen) then
            tcomp(ip) = tcomp(ip)/(3.0*ntype(ip) -3.0)
         else
            tcomp(ip) = tcomp(ip)/(3.0*ntype(ip))
         endif
      enddo

      return
      end
******************* end of subroutine init_vel *************************
************************************************************************
*
      real*8function gauss_dev(dummy)
*
*  This is a routine to do a Box-Muller transformation on random deviates
*  to return a normally distributed deviate with zero mean and unit 
*  variance.  This was taken from Numerical Recipes p. 203.
************************************************************************
      implicit none
      real*8 zeta, v1, v2, radius, factor, ranmar, dummy
      integer iflag
      data iflag /0/
      save zeta
      save iflag
     

      if (iflag .eq. 0) then
*  Generate random deviates to use in the transform
100      v1 = 2.0d00*ranmar() - 1.0d00
         v2 = 2.0d00*ranmar() - 1.0d00
*  See if they are in the unit circle
         radius = v1*v1 + v2*v2
         if (radius .ge. 1.0d00) goto 100
*  Make the Box-Muller transformation to get two normal deviates
         factor = sqrt(-2.0d00*log(radius)/radius)
         gauss_dev = v1*factor
*  Save this deviate for use next time.
         zeta = v2*factor
         iflag = 1
      else
*  We have an extra deviate so use that one.
         gauss_dev = zeta
         iflag = 0
      endif
      return
      end
******************** end of function gauss_dev *************************
************************************************************************
*
      subroutine RMARIN(IJ,KL)
*
*
*  The limits of IJ and KL are:
* 	    0 <= IJ <= 31328
*           0 <= KJ <= 30081
*
      real*8 U(97), C, CD, CM
      real*8 s, t
      integer I97, J97, IJ, KL
      integer i, j, k, l, ii, jj, m
      logical test
      common /raset1/ U
      common /raset2/ C, CD, CM
      common /raset3/ I97, J97
      common /raset4/ test

      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.
     *    KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a',
     * ' value between 0 and 31328'
          print '(A)',' The second seed must have a value between 0',
     *  '  and 30081'
            stop
      endif

      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169) 

      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         U(ii) = s
2     continue

      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0

      I97 = 97
      J97 = 33

      test = .TRUE.
      return
      end
*
************** end of subroutine RMARIN ********************************
************************************************************************
*
      function ranmar()
*
C This is the random number generator proposed by George Marsaglia in 
C Florida State University Report: FSU-SCRI-87-50
C It was slightly modified by F. James to produce an array of pseudorandom
C numbers.

      real*8 U(97), C, CD, CM
      real*8 ranmar, uni
      integer I97, J97
      logical test
      common /raset1/ U
      common /raset2/ C, CD, CM
      common /raset3/ I97, J97
      common /raset4/ test
 
      if( .NOT. test ) then
         print '(A)',' Call the init routine (RMARIN) before calling', 
     *    ' RANMAR'  
         stop
      endif


         uni = U(I97) - U(J97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         U(I97) = uni
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         uni = uni - C
         if( uni .lt. 0.0 ) uni = uni + 1.0
         RANMAR = uni

      return
      end
*
************** end of function RANMAR **********************************
************************************************************************
*
      subroutine lattice(npart, rx, ry, rz, type, np, adjust)
*
*  This program generates a lattice structure for atoms
*
************************************************************************
*
*                 VARIABLE IDENTIFICATION
*
*  NAC		Number of atoms per unit cell
*  NAX		Number of unit cells in the X direction
*  NAY		Number of unit cells in the Y direction
*  NAZ		Number of unit cells in the Z direction
*               (For a cubic lattice NAX = NAY = NAZ)
*  NP		Total number of particles in the system
*		(NP = NAC*NAX*NAY*NAZ)
*  RX(I)        X position of the Ith atom
*  RY(I)        Y position of the Ith atom
*  RZ(I)        Z position of the Ith atom
*  XFRAC(I)	Fractional coordinate of the Ith atom in the unit cell
*  YFRAC(I)	Fractional coordinate of the Ith atom in the unit cell
*  ZFRAC(I)	Fractional coordinate of the Ith atom in the unit cell
*  TYPE 	Type of lattice to be generated: FCC, BCC or SC.
*
************************************************************************
*
      IMPLICIT NONE
      integer npart, np, ncx, ncy, ncz, nac, i, ix, iy, iz, j
      real*8 rx, ry, rz, xfrac, yfrac, zfrac
      CHARACTER TYPE*3
      logical adjust
      DIMENSION RX(npart), RY(npart), RZ(npart), 
     & XFRAC(4), YFRAC(4), ZFRAC(4)

      adjust = .false.
*
*  Assign varriables according to the type of matrix selected
*
 100  CONTINUE
      IF ((TYPE .EQ. 'FCC') .OR. (TYPE .EQ. 'fcc')) THEN
        NAC = 4
        XFRAC(1) = 0.0
        YFRAC(1) = 0.0
        ZFRAC(1) = 0.0
        XFRAC(2) = 0.5
        YFRAC(2) = 0.0
        ZFRAC(2) = 0.5
        XFRAC(3) = 0.5
        YFRAC(3) = 0.5
        ZFRAC(3) = 0.0
        XFRAC(4) = 0.0
        YFRAC(4) = 0.5
        ZFRAC(4) = 0.5
      ELSE IF ((TYPE .EQ. 'BCC') .OR. (TYPE .EQ. 'bcc')) THEN
        NAC = 2
        XFRAC(1) = 0.0
        YFRAC(1) = 0.0
        ZFRAC(1) = 0.0
        XFRAC(2) = 0.5
        YFRAC(2) = 0.5
        ZFRAC(2) = 0.5
      ELSE IF ((TYPE .EQ. 'SC') .OR. (TYPE .EQ. 'sc')) THEN
        NAC = 1           
        XFRAC(1) = 0.5
        YFRAC(1) = 0.5
        ZFRAC(1) = 0.5
      else
         write(*,*) '*** Unknown lattice type ***'
         write(*,*) 'Assuming FCC '
         TYPE = 'FCC'
         GOTO 100
      ENDIF
*
*  Check to see if this configuration is acceptable
*
      NCX = nint((npart/nac)**(1.0d00/3.0d00))
      NCY = NCX
      NCZ = NCX
      NP = NAC*NCX*NCY*NCZ
      if (npart .ne. np) then
*  We need to adjust the lattice:
         adjust = .true.
         if (npart .gt. np) then 
*     Make sure we have a lattice that is larger than the required number of atoms
  5         ncx = ncx + 1
            ncy = ncx
            ncz = ncx
            np = nac*ncx*ncy*ncz
            if (npart .gt. np) goto 5
         endif
      endif
*
*  Place the atoms in a lattice using the box (program) units of
*  -0.5 to 0.5 for the dimensions of the box
*
      I = 0
      DO 40 IZ = 1, NCZ
        DO 30 IY = 1, NCY
          DO 20 IX = 1, NCX
            DO 10 J = 1, NAC
              I = I + 1
              RX(I) = - 0.5 + (IX - 1 + XFRAC(J))/FLOAT(NCX)
              RY(I) = - 0.5 + (IY - 1 + YFRAC(J))/FLOAT(NCY)
              RZ(I) = - 0.5 + (IZ - 1 + ZFRAC(J))/FLOAT(NCZ)
10          CONTINUE
20        CONTINUE
30      CONTINUE
40    CONTINUE
      return
      END
**************** end of subroutine lattice *****************************
************************************************************************
*
      subroutine assign_type(np, ncomp, ntype, itype)
*
*  Randomly assigns a type to each of the molecules in the system.
*  This subroutine is called only for new configurations.
*
************************************************************************
      implicit none
      integer np,ncomp,ntype(ncomp),itype(np),ip,count,num_trys,i_try
      real*8 ranmar

*  Initialize the itype array:
      do ip = 1, np
         itype(ip) = 0
      enddo

*  Assign each type
      num_trys = 0
      do ip = 1, ncomp
         count = 0
 100     if (count .lt. ntype(ip)) then
            num_trys = num_trys + 1
            i_try = np*ranmar()+1
            if (itype(i_try) .eq. 0) then
               itype(i_try) = ip
               count = count + 1
            endif
            if (num_trys .gt. 100000*np)then
               write(*,*)  '     Error in assign_type'
               write(3,*)  '     Error in assign_type'
               stop
            endif
            goto 100
         endif
      enddo

      return
      end
************** end of subroutine assign_type ***************************
************************************************************************
*
      subroutine widom(flag)
*
*  This procedure calculates the chemical potential via Widom's 
*  test particle insertion technique.  Also computes the cavity 
*  correlation function, based on the method of Lee & Shing, 
*  J.C.P. 91, 477-488 (1989).
*
************************************************************************
*            History
*
*  31 July 1996  Created from the mdlj.f widom subroutine.
*  18 Sept 1997  Allows different number of random vectors zero and 
*                non-zero separations.
*
*  flag is a control variable:
*  flag == 0 --> calculate the long range corrections only
*  flag == 1 --> find the largest diameter sphere and compute ryij(i)
*  flag == 2 --> zero the accumulators
*  flag != (0,1,2) perform yij calculations
*  Computes the yij(0) values from Lee & Shing, JCP 91, 477 (1989), 
*         Eq. (4.2) or (4.3)
*  Computes the yij(r) values, possibley from the method of 
*         Llano-Restrepo & Chapman, JCP, 97, 2046 (1992), Eq. (12)
*      J Karl Johnson
************************************************************************
      include 'var_mdljmix4.f'
      integer flag

*  Local variables
      real*8 xt0,yt0,zt0,ranmar,big,lrcs(ncmax),rc3,rc9,sigij3,
     &     xtn,ytn,ztn,upot0(ncmax),upotn(ncmax),xix,xiy,xiz
      integer ic, ip, jp, kp, lp, itp, jtp, count
      save lrcs

      if (flag .eq. 0) then
         do ic = 1, ncomp
*  Zero the long-range corrections:
            lrcs(ic) = 0.0
         enddo
         if (cut_and_shift) then
*  Do not compute long-range corrections if using a cut and shifted potential
            return
         endif
*  Calculate the long range corrections for each component
         rc3 = (sigbox/cutoff)**3
         rc9 = rc3*rc3*rc3
         do ic = 1, ncomp
            ntype(ic) = ntype(ic)+1
            do ip = 1, ncomp
               do jp = 1, ncomp
                  sigij3 = sig(ip,jp)
                  sigij3 = sigij3*sigij3*sigij3
                  lrcs(ic)=lrcs(ic)+ntype(ip)*ntype(jp)*eps(ip,jp)*
     &                 sigij3*sigij3*(rc9*sigij3*sigij3-3.0d0*rc3)
               enddo
            enddo
            lrcs(ic) = lrcs(ic)*8.0*pi*rho/(9.0*np)-ulrc*np
            ntype(ic) = ntype(ic) - 1
         enddo
         return
      else if (flag .eq. 1) then
*  Find the largest value of sig(i,j)
         big = 0.0
         do ip = 1, ncomp
            do jp = 1, ncomp
               if (sig(ip,jp) .gt. big) big = sig(ip,jp)
            enddo
         enddo
*  Compute the number of r-points for yij(r)
         ybins = big*yrmax*sigbox*rdels
         if (ybins .gt. nbinmax) ybins = nbinmax
         do ic = 0, ybins
            ryij(ic) = delr*sigbox*ic
         enddo
         return
      else if (flag .eq. 2) then
*  Zero the yij0 accumulator
         do ip = 1, ngij
            yij0(ip) = 0.0
            yijave(ip) = 0.0
            yijsq(ip) = 0.0
         enddo
*  Zero the chemical potential accumulators
         do ip = 1, ncomp
            cpot(ip) = 0.0
            cpave(ip) = 0.0
            cpsq(ip) = 0.0
            do jp = 1, ncomp
               do lp = 0, ybins
                  yijr(lp,ncomp*(ip-1)+jp) = 0.0
                  yrcount(lp,ncomp*(ip-1)+jp) = 0
               enddo
            enddo
         enddo
         return
      endif


*  Enter the loop for the number of test particles.
      do ic = 1, ntp
*  Randomly pick positions
         xt0 = ranmar() - 0.5d00
         yt0 = ranmar() - 0.5d00
         zt0 = ranmar() - 0.5d00
*  Compute the test particle energy for all molecule types:
         call one_move(xt0, yt0, zt0, upot0, np+1)
*  Loop through each component
         do ip = 1, ncomp
*  Include the long range corrections
            upot0(ip) = upot0(ip) + lrcs(ip)
*  Collect the argument for the chemical potentials
            cpot(ip) = cpot(ip) + exp(-upot0(ip)/tstar)
         enddo
*  Update the counter for chemical potentials
         counts(8) = counts(8) + 1
*  Collect the yij(0) terms:
         kp = 0
         do ip = 1, ncomp
            do jp = ip, ncomp
               kp = kp + 1
               yij0(kp)=yij0(kp)+exp(-(upot0(ip)+upot0(jp))/tstar)
            enddo
         enddo
*  Update the counter for yij(0)
         counts(9) = counts(9) + 1
*  End of test particle loop
      enddo

*  Loop through all the particles in the system to collect yij(r).
      do jp = 1, np
*  Find the type of particle jp
         jtp = itype(jp)
*  Compute the energies for yij(0).
         call one_move(xx(jp), yy(jp), zz(jp), upot0, jp)
         do itp = 1, ncomp
            upot0(itp) = upot0(itp) + lrcs(itp)
*  Compute which histogram to put this in:
            kp = ncomp*(itp-1)+jtp
            yijr(0,kp)=yijr(0,kp)+exp(-upot0(itp)/tstar)
*  Update the yij(r) counter
            yrcount(0,kp) = yrcount(0,kp) + 1
         enddo
*  Loop through the number of r vectors:
         do lp = 1, ybins
*  If the bond length is greater than some value then 
*  loop through a number of randomly selected angles:
            if (ryij(lp) .gt. rytest) then
               count = nangles
            else
               count = 1
            endif   
            do ic = 1, count
               call vector3d(xix,xiy,xiz,ryij(lp))
               xtn = xx(jp) + xix
               ytn = yy(jp) + xiy
               ztn = zz(jp) + xiz
*  Compute the energies for all components:
               call one_move(xtn, ytn, ztn, upotn, jp)               
               do itp = 1, ncomp
*  Add the long-range corrections:
                  upotn(itp) = upotn(itp) + lrcs(itp)
*     Compute which histogram to put this in:
                  kp = ncomp*(itp-1)+jtp
                  yijr(lp,kp)=yijr(lp,kp)+exp(-upotn(itp)/tstar)
*  Update the yij(r) counter
                  yrcount(lp,kp) = yrcount(lp,kp) + 1
               enddo
            enddo
         enddo
      enddo

      return
      end
***************** end of subroutine widom ******************************
************************************************************************
*
      subroutine one_move(xi, yi, zi, ui, i_cur)
*
*  Calculate the potential energy change due to moving one molecule.
*
************************************************************************
      include 'var_mdljmix4.f'
*  Local variables
      real*8 xi,yi,zi,xij,yij,zij,rij2,ui(ncomp),r12,r6,sigij6,epsij
      integer i_cur, j,itp,jtp

      do itp = 1, ncomp
         ui(itp) = 0.0
      enddo

*  Loop over all other particles.
      do 120 j = 1, np
         if (j .ne. i_cur) then
*  Calculate the distance vector between the centers of mass for
*  molecules i_cur and j.
            xij = xi - xx(j)
            yij = yi - yy(j)
            zij = zi - zz(j)
            jtp = itype(j)
*  Invoke the Minimum Image Convention.
            if (xij .le. -0.5) then
                xij = xij + 1.0d00
            elseif (xij .ge. 0.5) then
                xij = xij - 1.0d00
            endif
            if (yij .le. -0.5) then
                yij = yij + 1.0d00
            elseif (yij .ge. 0.5) then
                yij = yij - 1.0d00
            endif
            if (zij .le. -0.5) then
                zij = zij + 1.0d00
            elseif (zij .ge. 0.5) then
                zij = zij - 1.0d00
            endif
*  End of minimum image convention
*  Calculate the square of the separation distance
            rij2 = xij*xij + yij*yij + zij*zij
*  Check to see if the distance is within the cutoff
            if (rij2 .lt. cutsq) then
               do itp = 1, ncomp
                  sigij6 = sig(itp,jtp)
                  sigij6=sigij6*sigij6*sigij6
                  sigij6=sigij6*sigij6
                  epsij = eps(itp,jtp)
                  r6 = sigsq/rij2
                  r6 = r6*r6*r6*sigij6
                  r12 = r6*r6
                  ui(itp) = ui(itp) + epsij*(r12 - r6)
                  if (cut_and_shift) then
*  Subtract off the shift
                     ui(itp) = ui(itp) - epsij*pshift6*sigij6*(
     &                    pshift6*sigij6 - 1.0d0)
                  endif
               enddo
            endif
*  End of the "if (j .ne. i_cur)" statement
         endif
120   continue
      do itp = 1, ncomp
         ui(itp) = 4.0d00*ui(itp)
      enddo

      return
      end
****************** end of subroutine one_move **************************
************************************************************************
*
      subroutine cavity(sub1,rsub1,rsub2)
*
*  Calculates the chemical potentials and the cavity correlation
*  functions with error bars.
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 sub1,rsub1,rsub2
*  Local variable
      integer ix, jx, kx, lx
      real*8 error_bars(ncmax*ncmax),lnyij(ncmax*(ncmax+1)/2)

      do ix = 1, ncomp
         cpave(ix) = cpave(ix)*rsub1
         cpsq(ix) = (cpsq(ix)-sub1*cpave(ix)**2)*rsub2
      enddo
      kx = 0
      do ix = 1, ncomp
         do jx = ix, ncomp
            kx = kx + 1
            yijave(kx)=yijave(kx)*rsub1
            yijsq(kx)=(yijsq(kx)-sub1*yijave(kx)**2)*rsub2
            yij0(kx) = log(yijave(kx)/(cpave(ix)*cpave(jx)))
         enddo
      enddo
*  Compute the error bars:
      kx = 0
      do ix = 1, ncomp
         do jx = ix, ncomp
            kx = kx + 1
            error_bars(kx)=sqrt(cpsq(ix))/cpave(ix)+
     &           sqrt(cpsq(jx))/cpave(jx)+sqrt(yijsq(kx))/
     &           yijave(kx)
         enddo
      enddo

*  Write out the results to a file

      open(unit=55, file='yij.dat', status='unknown')
      write(55,1000) 0.0, (exp(yij0(kx)),exp(yij0(kx)+error_bars(kx))-
     &  exp(yij0(kx)),exp(yij0(kx))-exp(yij0(kx)-error_bars(kx)),
     &  kx=1,ngij)
 1000 format(20(1x,e13.6))
      do lx = 0, ybins
         kx = 0
         do ix = 1, ncomp
            do jx = 1, ncomp
               kx = kx + 1
               if (yrcount(lx,kx) .gt. 0) then
                  lnyij(kx) = log(yijr(lx,kx)/(yrcount(lx,kx)*
     &                 cpave(ix)))
               else
*  Prevent divide by zero:
                  lnyij(kx) = 0
               endif
            enddo
         enddo
         write(55,1000) ryij(lx)/sigbox, (exp(lnyij(kx)),kx=1,
     &        ncomp*ncomp)
      enddo

      return
      end
********************* end of subroutine cavity *************************
************************************************************************
*
      subroutine vector3d(xx, yy, zz, rad)
*
*  Generate a random vector in 3-d, normalize to lie on a circle.
*
*  Corrected error: 8/20/97
*  Implements algorithm of Marsaglis taken from Allen & Tildesley.
*
************************************************************************
      implicit none
      real*8 xi1, xi2, zeta1, zeta2, xx, yy, zz, norm, rad, ranmar

* Find two random numbers whose length is less than one:
 100  xi1 = ranmar()
      xi2 = ranmar()
      zeta1 = 1-2.0*xi1
      zeta2 = 1-2.0*xi2
      norm = zeta1*zeta1 + zeta2*zeta2
      if (norm .gt. 1.0) goto 100
* Now compute the vector elements, with total length rad:
      norm = sqrt(1.0-norm)
      xx = 2.0*zeta1*norm
      yy = 2.0*zeta2*norm
      zz = norm*norm
      norm = xx*xx + yy*yy + zz*zz
      norm = rad/sqrt(norm)
      xx = xx*norm
      yy = yy*norm
      zz = zz*norm
      
      return
      end
****************** end of subroutine vector3d **************************
************************************************************************
*
*
      subroutine a_gofr(pr,up)
*
*  This procedure analyzes the data collected in the g(r) histograms
*  and writes the results to a file.
*
************************************************************************
      include 'var_mdljmix4.f'
*  Local variables
      real*8 rad_i, rad_ip1, vshell, radave, factor,term1, term2,
     &    pr, up, gofr(ncmax*(ncmax+1)/2) 
      integer ishell,ip,jp,kp
      real*8 volume_element
      real*8 coord(ncmax*(ncmax+1)/2)
**
*  Open the output file for the g(r) calculation results
*
      open(unit=2, file='gofr.mdljmix', status='unknown')
*
*  Open the output file for the coordination numbers
      open(unit=72, file='coord.mdljmix', status='unknown')
*
*  Zero the coordination counters
      do kp = 1, ngij
         coord(kp) = 0
      enddo
*  Calculate the conversion factor.
      factor = 1.0d0/(rho*rnp*float(counts(4)))
      write(2, *) rho, tstar, pr, up
      write(2, 1000) rho, tstar, cutoff/sigbox
1000  format(5x,'Results from the g(r) calculations',
     &       ' with rho* =', f10.3, 2x, 'T*=', f8.4, 1x, 'and cutoff=',
     &      f8.5, /, '!', 1x, 'Shell', 2x, 'R(ave)', 5x, 'g(r)')
      do 100 ishell = 1, nbins
         rad_i = float(ishell - 1)*delr
         rad_ip1 = rad_i + delr
         radave = 5.0d-01*(rad_i + rad_ip1)
* Compute the volume of the shell based on the probability density for 
* a cubic box.  The idea for this was taken from Walter Chapman's thesis.
         term1 = volume_element(rad_i)
         term2 = volume_element(rad_ip1)
         vshell = term2 - term1
         kp = 0
         do ip = 1, ncomp
            do jp = ip, ncomp
               kp = kp + 1
               if ((xmol(ip) .gt. 0.0) .and. (xmol(jp) .gt. 0.0)) then
                  term1 = factor/(vshell*xmol(ip)*xmol(jp))
               else
                  term1 = 0.0
               endif   
               gofr(kp) = hist(ishell,kp)*term1
*  Compute the coordination number:
               coord(kp) = coord(kp) + gofr(kp)*vshell*rho*xmol(ip)*
     &               xmol(jp)
            enddo
         enddo
         write(2,1010) ishell, radave, (gofr(ip),ip=1,ngij)
1010     format(1x, i4, 20(1x, f10.6))
         write(72,1020) ishell, radave, (coord(ip),ip=1,ngij)
1020     format(1x, i4, 20(1x, e13.6))
100   continue

      close(2)
      close(72)

      return
      end
*
****************** end of subroutine a_gofr ****************************
************************************************************************
*
      real*8 function volume_element(rr)
*
************************************************************************
      include 'var_mdljmix4.f'
      real*8 rr, alpha, sqrt2, sqrt3, lambda, T1, T2

      alpha = 0.5/sigbox
      sqrt2 = sqrt(2.0d0)
      sqrt3 = sqrt(3.0d0)

      if (rr .le. alpha) then
         volume_element = 4.0d0*pi*rr**3/3.0d0
      else if (rr .le. sqrt2*alpha) then
         volume_element = -2.0d0*pi*alpha**3 + 6.0d0*pi*alpha*rr**2 -
     &        8.0d0*pi*rr**3/3.0d0
      else if (rr .le. sqrt3*alpha) then
         lambda = rr/alpha
         T1 = sqrt(lambda*lambda-1)
         T2 = sqrt(lambda*lambda-2)
         volume_element = 8.0*alpha**3*(T2+(lambda*lambda-1.0/3.0)*(2.0*
     &        asin(1.0/T1)-asin(T2/T1)-0.25*pi)-lambda**3/3.0*(
     &        asin((lambda-T2)/(2.0*T1))-asin((lambda+T2)/(2.0*T1))
     &        - asin((lambda*lambda-lambda-1)/((lambda-1)*T1)) +
     &        asin((lambda*lambda+lambda-1)/((lambda+1)*T1))))
      else
         write(*,*) 'Error in a_gofr integration'
      endif
      return
      end
******************** End of function volume_element ********************
************************************************************************
*
      subroutine adjust_lattice(n1, n2, xx, yy, zz)
*
* Routine to delete atoms from the lattice at random and close up holes 
* in the resulting list.
*
************************************************************************
      implicit none
      integer n1, n2, ndel, ix, jx, ikill
      real*8 xx(n2), yy(n2), zz(n2), ranmar, xt, yt, zt

      ndel = n2 - n1
      do ix = 1, ndel
         ikill = n2*ranmar()+1
         n2 = n2 - 1
         if (ikill .lt. n2) then 
            do jx = ikill, n2
               xx(jx) = xx(jx+1)
               yy(jx) = yy(jx+1)
               zz(jx) = zz(jx+1)
            enddo
         endif
      enddo
      return
      end
***************** end of subroutine adjust_lattice *********************
************************************************************************
