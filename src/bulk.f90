! =============================================================================
! =============================================================================
!
! MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC IONIC SOLUTION
!
! =============================================================================
! =============================================================================


program bulk

   implicit double precision (a-h, o-z)
   real :: ran2
   include 'bulk_f90.inc'

   dimension :: uula(2), uuta(2), suu(25), tijd(10)
   dimension :: mtot(mxspec), macc(mxspec), menrj(mxspec), mhcrj(mxspec)

   data uuta  / 2*0.0 /
   data suu   / 25*0.0 /
   data mtot, macc, menrj, mhcrj / 40*0 /
   data kn1, kn2, key, kntot / 4*0 /
   data tijd  / 10*0.0 /

43 format( 'MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC    ', &
            'IONIC SYSTEM  ', /, /, 'Written by Peter Bolhuis, February', &
            '1992', / )

44 format( '************************************************************', &
            '************************', / )

   iclk = 0
   lclk = 0

   kkk = 12
   lll = 13
   mmm = 15
   jjj = 6

   write(jjj, 44)
   write(jjj, 43)
   write(jjj, 44)


   ! ==========================================================================
   ! Physical constants and conversion factors
   ! ==========================================================================

   pi   = acos(-1.0d0)
   epsx = 8.85418782d-12
   ech  = 1.60219d-19
   avno = 6.0223d23
   bk   = 8.314


   ! ==========================================================================
   ! Sample units are angstrom and elementary charge
   ! Kilojoules are used for a mole
   !
   ! MEANING OF THE INPUT VARIABLES
   !
   ! ink         = 0   from file
   ! ink         = 1   from slump
   ! nspec       = number of ionic components
   ! hion        = array contains the input info on the ions
   !               (1=nr, 2=number, 3=radius, 4=charge, 5=displacement)
   ! ny(1,2,3)   = number of MC loops
   ! dtemp       = temperature
   ! eps         = relative dielectric constant
   ! box         = box size
   ! nwint       = a widom/collision calc is performed every nwint cycle
   ! nwins       = number of widom particles inserted each calc.
   ! nfix        = where is the chem. pot. measured
   !               0 = averaged over the cell
   !
   ! ==========================================================================
   ! ==========================================================================

   open(lll, file='bulk.conf', form='formatted')
   open(mmm, file='bulk.inp',  form='formatted')
   open(kkk, file='bulk.macro', form='formatted')

   read(mmm, *)
   read(mmm, *) nspec
   read(mmm, *)
   read(mmm, *)
   read(mmm, *) (hion(l, 1), hion(l, 2), hion(l, 3), hion(l, 4), hion(l, 5), l = 1, nspec)
   read(mmm, *)
   read(mmm, *)
   read(mmm, *) ny1, ny2, ny3
   read(mmm, *)
   read(mmm, *)
   read(mmm, *) ink, islu
   islu = -islu
   read(mmm, *)
   read(mmm, *)
   read(mmm, *) box
   read(mmm, *)
   read(mmm, *)
   read(mmm, *) dtemp, eps
   read(mmm, *)
   read(mmm, *)
   read(mmm, *) nwins, nwint, nfix

   if (nfix > 2)  nfix = 2
   if (islu == 0) islu = 7
   if (nwint == 0) nwint = 2 * ny1


   ! ==========================================================================
   ! Boxsize variables
   ! ==========================================================================

   box2  = box / 2.0
   box2i = 1.0 / box2


   ! ==========================================================================
   ! Set up hard core and charge vectors
   !
   ! hc2v contains in the first column the particle type given
   ! by hion(i,1), whereas in the following columns the squared
   ! hard core distances between the particles are given.
   ! chv gives for every particle its charge
   ! npart is the total number of particles
   ! ==========================================================================

   npart = 0

   do i = 1, nspec
      num = int(hion(i, 2))

      do j = 1, num
         npart = npart + 1
         ispc(npart) = i

         do k = 1, nspec
            hc2v(npart, k) = (hion(i, 3) + hion(k, 3)) ** 2
         end do

         chv(npart) = hion(i, 4)
         dp(npart)  = hion(i, 5)
      end do
   end do

   do i = 1, nspec
      caver(i) = (hion(i, 2)) / (box ** 3)
   end do

   do i = 1, npart
      ei(i) = 0.0
   end do

   nkonf = ny1 * ny2 * ny3

   if (ny3 > 25) ny3 = 25

   t     = dtemp
   abeta = -1.0 / (bk * dtemp)


   ! ==========================================================================
   ! Read configuration from file
   ! ==========================================================================

   if (ink == 0) then
      read(lll, *) (x6(l), y6(l), z6(l), l = 1, npart)
      write(jjj, *) 'Configuration read from file'
   end if

   if (ink == 1) call slump
   call collision
   call widom


   ! ==========================================================================
   ! Write input data
   ! ==========================================================================

   nwtot = nwins * (int(ny1 / nwint) * ny2 * ny3)

   write(jjj, 800) npart
   write(jjj, 801) (l, int(hion(l, 2)), hion(l, 3), hion(l, 4), &
                    hion(l, 5), caver(l) * 1.0d27 / avno, l = 1, nspec)
   write(jjj, 802) eps, dtemp, box, ny1, ny2, ny3, nkonf, &
                   dble(nkonf / npart), nwtot

800 format (/, 'VARIABLES FROM INPUT FILE', /, &
            /, ' number of particles        =', i10, /, &
            /, ' species   number   radius    charge   dp   average conc', /)
801 format (i6, i10, f10.1, f10.1, f10.2, f10.6)
802 format (/, ' dielectric constant         =', f10.1, /, &
            ' temperature                 =', f10.1, /, &
            ' box size (AA)               =', f10.1, /, &
            ' number of configurations    =', i6, ' *', i6, ' *', i6, ' =', i8, /, &
            ' number of conf. per part    =', f10.1, /, &
            ' widom test per species      =', i10, /)

   ecf = ech ** 2 * 1.0d10 * avno / (eps * epsx * 4.0 * pi)


   ! ==========================================================================
   ! Energy evaluation initial configuration
   ! ==========================================================================

   call liv

   ytot1 = xww1 * ecf * 1.0d-3

   write(jjj, '(a)') 'Initial configuration'
   write(jjj, '(a, e12.5)') 'Coulomb energy (kJ/mol)  =', ytot1


   ! ==========================================================================
   ! Start of Monte Carlo loop
   !
   ! Loop over ny1*ny2*ny3 particle moves.
   ! il    = current particle
   ! ei    = new coulomb energy of particle il
   ! esa   = all coulomb interaction energies of the particle
   ! xww1  = total energy (also ulaa(1))
   ! dp(1..10)  displacement parameter for every ion
   !
   ! ==========================================================================

   do my3 = 1, ny3
      uula(1) = 0.0

      do my2 = 1, ny2
         do my1 = 1, ny1
            il    = 1 + mod(kntot, npart)
            ispec = ispc(il)
            kntot = kntot + 1
            mtot(ispec) = mtot(ispec) + 1

            tx6 = x6(il) + dp(il) * (ran2(islu) - 0.5)
            ty6 = y6(il) + dp(il) * (ran2(islu) - 0.5)
            tz6 = z6(il) + dp(il) * (ran2(islu) - 0.5)

            if (tx6 >  box2) tx6 = tx6 - box
            if (tx6 < -box2) tx6 = tx6 + box
            if (ty6 >  box2) ty6 = ty6 - box
            if (ty6 < -box2) ty6 = ty6 + box
            if (tz6 >  box2) tz6 = tz6 - box
            if (tz6 < -box2) tz6 = tz6 + box

            isos = 99
            call qin

            if (isos /= 0) go to 63

            deltae = 0.0
            do j = 1, npart
               deltae = deltae + ei(j) - esa(j, il)
            end do

            dzxp = abeta * (deltae * ecf)

            if (dzxp < -80.0) go to 64
            if (dzxp >   0.0) go to 62
            if (exp(dzxp) < ran2(islu)) go to 64

62          macc(ispec) = macc(ispec) + 1

            ! Trial config accepted
            x6(il) = tx6
            y6(il) = ty6
            z6(il) = tz6

            do j = 1, npart
               esa(il, j) = ei(j)
               esa(j, il) = ei(j)
            end do

            xww1 = xww1 + deltae
            go to 65

63          mhcrj(ispec) = mhcrj(ispec) + 1
            go to 65

64          menrj(ispec) = menrj(ispec) + 1

65          uula(1) = uula(1) + xww1

            if (mod(my1, nwint) == 0) then
               call collision1
               call widom1
            end if

         end do
      end do

      qww1 = xww1
      call liv

      if (xww1 /= 0) qww1 = (qww1 - xww1) / xww1

      write(kkk, *)
      write(kkk, 44)
      write(kkk, 73) my3, qww1, qfww1

      if (abs(qww1) > 0.001) stop

      call liv

      if (nwint <= ny1) then
         call collision2
         call widom2
      end if

      uula(1) = uula(1) / dble(ny1 * ny2)
      uuta(1) = uuta(1) + uula(1)
      uula(1) = uula(1) * ecf * 1.0d-3
      suu(my3) = uula(1)

      ! Calculate total acceptance and rejection
      key = 0
      kn1 = 0
      kn2 = 0

      do i = 1, nspec
         key = key + macc(i)
         kn1 = kn1 + mhcrj(i)
         kn2 = kn2 + menrj(i)
      end do

      write(kkk, 810) kntot, key, kn2, kn1
      write(kkk, '(a, e12.5)') 'Coulomb energy (kj/mol)  =', uula(1)
   end do

73  format(/ ' MACROSTEP  ', i10, /, /, 'Checkparameters : ', 3e10.3)
810 format(/, ' total number of configurations    =', i8, /, &
           ' accepted configurations           =', i8, /, &
           ' energy rejected configurations    =', i8, /, &
           ' hard core rejected configurations =', i8, /)
811 format(/, 'Divided per species', /, 'species      accepted   ', &
           ' energy rej.   hard core rej.', /, 10(i6, 3(6x, f8.4), /), /)


   ! ==========================================================================
   ! End of Monte Carlo loop
   !
   ! Calculation of the average energies and the pressures
   !
   ! ==========================================================================

   yn3     = 1.0 / float(ny3)
   uuta(1) = yn3 * uuta(1)
   uuta(1) = uuta(1) * ecf * 1.0d-3

   call earth(suu, uuvar, uuta(1), ny3)

   write(jjj, '(/)')
   write(jjj, 44)
   write(jjj, '(/, a, i3, a, /)') 'FINAL RESULTS AFTER ', ny3, &
                                  ' MACROSTEPS'
   write(jjj, 810) kntot, key, kn2, kn1
   write(jjj, 811) (l, dble(macc(l)) / mtot(l), dble(menrj(l)) &
                    / mtot(l), dble(mhcrj(l)) / mtot(l), l = 1, nspec)
   write(jjj, '(a, 2e12.5)') 'Coulomb energy (kj/mol)  =', uuta(1), uuvar
   write(jjj, *)

   if (nwint <= ny1) then
      call collision3
      call widom3
   end if

   rewind lll
   write(lll, 771) (x6(l), y6(l), z6(l), l = 1, npart)
771 format(5e16.8)

   close(lll)


   ! ==========================================================================
   !
   ! Calculation of pressures from the densities
   !
   ! ==========================================================================

   pid = 0.0
   do k = 1, nspec
      pid = pid + caver(k) * 1.0d27 / avno
   end do

   pexen  = uuta(1) * 1.0d30 / (3 * bk * dtemp * box ** 3 * avno)
   pexenv = uuvar   * 1.0d30 / (3 * bk * dtemp * box ** 3 * avno)
   ptot   = pid + pcollav + pexen
   ptotv  = sqrt(pcollv * pcollv + pexenv * pexenv)

   write(jjj, '(/, a, /)') 'TOTAL BULK PRESSURE '
   write(jjj, 729) 'Ideal pressure       ', pid, 0.0
   write(jjj, 729) 'Energy con. <E/3V>   ', pexen, pexenv
   write(jjj, 729) 'Collision press      ', pcollav, pcollv
   write(jjj, 731)
   write(jjj, 729) 'Bulk pressure        ' , ptot, ptotv

729 format(a, f12.4, f10.4)
731 format(70('_'))

   stop

end program bulk



! =============================================================================
! =============================================================================
!
! Subroutine: earth
!
! Calculates the average and deviation of nnn values stored in
! array a(1..25)
!
! =============================================================================
! =============================================================================

subroutine earth(a, bbbwww, xb, nnn)

   implicit double precision (a-h, o-z)

   dimension :: a(25)

   yak = 1.0 / (nnn * (nnn - 1))
   b   = 0.0
   xb  = 0.0

   do k = 1, nnn
      xb = xb + a(k)
   end do

   xb = xb / nnn

   do k = 1, nnn
      b = b + (a(k) - xb) * (a(k) - xb)
   end do

   b = sqrt(yak * b)
   bbbwww = b

   return

end subroutine earth


! =============================================================================
! =============================================================================
!
! Subroutine: earth2
!
! Calculates the average and deviation of nnn values stored in
! matrix a(1..25,1..num) and stores them in arrays xb and bbbwww
!
! =============================================================================
! =============================================================================

subroutine earth2(a, bbbwww, xb, nnn, num)

   implicit double precision (a-h, o-z)

   parameter (mxspec = 10)
   dimension :: a(25, mxspec), xb(mxspec), bbbwww(mxspec)

   yak = 1.0 / (nnn * (nnn - 1))

   do i = 1, num
      b = 0.0
      xb(i) = 0.0

      do k = 1, nnn
         xb(i) = xb(i) + a(k, i)
      end do

      xb(i) = xb(i) / nnn

      do k = 1, nnn
         b = b + (a(k, i) - xb(i)) * (a(k, i) - xb(i))
      end do

      b = sqrt(yak * b)
      bbbwww(i) = b
   end do

   return

end subroutine earth2


! =============================================================================
! =============================================================================
!
! Subroutine: slump
!
! Set up of the random initial configuration.
!
! =============================================================================
! =============================================================================

subroutine slump

   implicit double precision (a-h, o-z)
   real :: ran2
   include 'bulk_f90.inc'

   write(jjj, *) 'Random configuration generated'

   nsl    = 0
   nslmax = 100000
   k12    = 0

1  nsl = nsl + 1

   if (nsl > nslmax) then
      write(jjj, *)' too dense system'
      stop
   end if

   x6tt  = (ran2(islu) - 0.5) * box
   y6tt  = (ran2(islu) - 0.5) * box
   z6tt  = (ran2(islu) - 0.5) * box
   ispec = ispc(k12 + 1)

   do i = 1, k12
      ddx = x6tt - x6(i)
      ddy = y6tt - y6(i)
      ddz = z6tt - z6(i)
      ddx = ddx - aint(ddx * box2i) * box
      ddy = ddy - aint(ddy * box2i) * box
      ddz = ddz - aint(ddz * box2i) * box
      r2  = ddx ** 2 + ddy ** 2 + ddz ** 2

      if (r2 < hc2v(i, ispec)) go to 1
   end do

   k12 = k12 + 1
   x6(k12) = x6tt
   y6(k12) = y6tt
   z6(k12) = z6tt

   if (k12 == npart) return

   go to 1

end subroutine slump


! =============================================================================
! =============================================================================
!
! Subroutine: qin
!
! Essential part of the program. Checks for overlap and calculate
! the energy change after a MC step.
!
! =============================================================================
! =============================================================================

subroutine qin

   implicit double precision (a-h, o-z)
   include 'bulk_f90.inc'

   do k = 1, npart
      ddx = abs(tx6 - x6(k))
      ddy = abs(ty6 - y6(k))
      ddz = abs(tz6 - z6(k))

      if (ddx > box2) ddx = ddx - box
      if (ddy > box2) ddy = ddy - box
      if (ddz > box2) ddz = ddz - box

      rw2(k) = ddx * ddx + ddy * ddy + ddz * ddz

      ! Cancel out the hard core overlap with itself
      rw2(il) = 1000000

      if (rw2(k) < hc2v(k, ispec)) return
   end do

   isos = 0
   chil = chv(il)

   do k = 1, npart
      ei(k) = chil * chv(k) / sqrt(rw2(k))
   end do

   ei(il) = 0

   return

end subroutine qin


! =============================================================================
! =============================================================================
!
! Subroutine: liv
!
! Gives the total coulombic energy and the force by recalculating
! every particleinteraction. Used to check QIN
!
! =============================================================================
! =============================================================================

subroutine liv

   implicit double precision (a-h, o-z)
   include 'bulk_f90.inc'

   xww1 = 0.0

   do i = 1, npart - 1
      tx6 = x6(i)
      ty6 = y6(i)
      tz6 = z6(i)
      ip1 = i + 1

      do k = ip1, npart
         ddx = tx6 - x6(k)
         ddy = ty6 - y6(k)
         ddz = tz6 - z6(k)
         ddx = ddx - aint(ddx * box2i) * box
         ddy = ddy - aint(ddy * box2i) * box
         ddz = ddz - aint(ddz * box2i) * box
         rw2(k) = ddx * ddx + ddy * ddy + ddz * ddz
      end do

      do k = ip1, npart
         uj1 = chv(i) * chv(k) / sqrt(rw2(k))
         xww1 = xww1 + uj1
         esa(i, k) = uj1
      end do
   end do

   do i = 1, npart - 1
      do k = i + 1, npart
         esa(k, i) = esa(i, k)
      end do
   end do

   do i = 1, npart
      esa(i, i) = 0.0
   end do

   return

end subroutine liv


! =============================================================================
!
! Subroutine: collision
!
! This routine checkes if a particle of type i in one halfcel is
! able to collide with a patricle of type j in the other. The
! contribution to the pressure is calculation by intergration
! At a random position around the real particle of type i a ghost
! particle of type j is inserted and the potential with every
! other particle in the system is calculated. This gives the
! collision probability
!
! =============================================================================
! =============================================================================

subroutine collision

   implicit double precision (a-h, o-z)
   real :: ran2
   include 'bulk_f90.inc'

   dimension :: rel(mxspec)
   dimension :: scoll(25, mxspec, mxspec), coll(mxspec, mxspec)

   do i = 1, mxspec
      do k = 1, mxspec
         coll(i, k) = 0.0

         do j = 1, ny3
            scoll(j, i, k) = 0.0
         end do
      end do
   end do

   return


   ! -------------------------------------------------------------------------

   entry collision1

   do i = 1, nwins
      num = 1

      do isp = 1, nspec
         nisp = int(hion(isp, 2))

         do ksp = 1, nspec
            nnn = ran2(islu) * nisp + num
            dis = hion(isp, 3) + hion(ksp, 3)
            aa  = (2 * ran2(islu) - 1) * dis
            v   = 2 * pi * ran2(islu)
            wz6 = z6(nnn) + aa
            g2  = sqrt(dis * dis - aa * aa) + 0.00001
            wx6 = x6(nnn) + g2 * cos(v)
            wy6 = y6(nnn) + g2 * sin(v)
            wx6 = wx6 - aint(wx6 * box2i) * box
            wy6 = wy6 - aint(wy6 * box2i) * box
            wz6 = wz6 - aint(wz6 * box2i) * box
            urej = 0

            do k = 1, npart
               ddx = dabs(wx6 - x6(k))
               ddy = dabs(wy6 - y6(k))
               ddz = dabs(wz6 - z6(k))

               if (ddx > box2) ddx = ddx - box
               if (ddy > box2) ddy = ddy - box
               if (ddz > box2) ddz = ddz - box

               rw2(k) = ddx * ddx + ddy * ddy + ddz * ddz

               if (rw2(k) < hc2v(k, ksp)) go to 121
            end do

            do k = 1, npart
               rwi(k) = 1.0 / sqrt(rw2(k))
            end do

            wtot2 = 0

            do k = 1, npart
               wtot2 = wtot2 + rwi(k) * chv(k)
            end do

            coll(isp, ksp) = exp(abeta * wtot2 * hion(ksp, 4) * ecf) + coll(isp, ksp)

121         continue
         end do

         num = num + nisp
      end do
   end do

   return


   ! -------------------------------------------------------------------------

   entry collision2

   nwtot = nwins * int(ny1 / nwint) * ny2

   write (kkk, '(/, /, a)') 'COLLISION PRESSURE MATRIX'
   write (kkk, '(/, a, i6)') 'Total collision trials per species ', nwtot
   write (kkk, 2010) (i, i = 1, nspec)

   do k = 1, nspec
      do i = 1, nspec
         scoll(my3, i, k) = coll(i, k) / nwtot
         coll(i, k) = 0
      end do

      write(kkk, 2012) (scoll(my3, i, k), i = 1, nspec)
   end do

   return

2010 format('Species          ', 10(12x, i6), /)
2011 format(/, i4, '    Col. samples', 10('            ', i6))
2012 format('        pressure    ', 10e18.6, /)
2013 format(i4, '        ', 20e12.5)


   ! -------------------------------------------------------------------------

   entry collision3

   write(jjj, '(/, a, /)') 'COLLISION MATRIX  AND RELATIVE ERROR'
   write (jjj, 2010) (i, i = 1, nspec)

   do i = 1, nspec
      do k = 1, nspec
         collav = 0

         do j = 1, ny3
            cwi(j, i, k) = scoll(j, i, k)
            dum(j) = scoll(j, i, k)
         end do

         call earth(dum, collv, collav, ny3)

         rel(k) = 0
         if (collav /= 0) rel(k) = collv / collav
         coll(i, k) = caver(i) * collav
      end do

      write(jjj, 2013) i, (coll(i, k), rel(k), k = 1, nspec)
   end do

   write(jjj, '(/)')

   return

end subroutine collision


! =============================================================================
! =============================================================================
!
! Subroutine: widom
!
! Puts a ghost particle anywhere in the cell and calculates the
! energy with the rest of the system. Averaging of the exponent
! results in the chemical potential. A rescaling of the ion-
! charges is neccessary to assure neutrality. The integration
! is performed in entry 2 with a simpson's rule.
! The total chemical potential is divided in to the ideal (chid),
! hard core (chhc) and the electric contribution (chel). For
! each species the chem.pot. is given cellaverage, at the wall,
! and in the midplan.
!
! =============================================================================
! =============================================================================

subroutine widom

   implicit double precision (a-h, o-z)
   real :: ran2
   include 'bulk_f90.inc'

   dimension :: chel(25, mxspec), chhc(25, mxspec), chex(25, mxspec)
   dimension :: chto(25, mxspec), chexw(25, mxspec), dch1(25, mxspec)
   dimension :: dch2(25, mxspec), chint(mxspec, 11), ewnom(mxspec, 11)
   dimension :: ewden(mxspec, 11), chinta(mxspec, 11)
   dimension :: expuw(mxspec), chid(mxspec)
   dimension :: chelav(mxspec), chelv(mxspec), chhcav(mxspec), chhcv(mxspec)
   dimension :: chexav(mxspec), chexv(mxspec), chtoav(mxspec), chtov(mxspec)
   dimension :: chexwa(mxspec), chexwv(mxspec)
   dimension :: ihc(mxspec), irej(mxspec), mwcn(25), ihcall(0:5)

   character*80 str(mxspec)

   do i = 0, nfix
      do j = 1, nspec
         if (caver(j) /= 0) then
            chid(j + nspec * i) = dlog( caver(j) / avno * 1.0d27)
         else
            chid(j + nspec * i) = -77
         end if
      end do
   end do

   do j = 1, mxspec
      do i = 1, 25
         chel(i, j) = 0
         chhc(i, j) = 0
         chex(i, j) = 0
         chto(i, j) = 0
      end do

      chelav(j) = 0
      chhcav(j) = 0
      chexav(j) = 0
      chtoav(j) = 0
      expuw(j)  = 0
      ihc(j)    = 0

      do k = 1, 11
         ewden(j, k)  = 0
         ewnom(j, k)  = 0
         chinta(j, k) = 0
      end do
   end do

   do j = 1, 25
      mwcn(j) = 0
   end do

   do j = 0, 5
      ihcall(j) = 0
   end do

   ntel = 0

   str(1) = 'Chemical potentials averaged over the whole cell'
   str(2) = 'Chemical potentials at the wall'
   str(3) = 'Chemical potentials at the midplane'

   return


   entry widom1

   mwcn(my3) = mwcn(my3) + 1

   do mp = 0, nfix
      do i = 1, nwins
         x = box * (ran2(islu) - 0.5)
         y = box * (ran2(islu) - 0.5)
         z = 0

         if (mp <= 1) then
            z = box * (ran2(islu) - 0.5)
            if (mp == 1) z = sign(box2, z)
         end if

         do j = 1, npart
            ddx = abs(x - x6(j))
            if (ddx > box2) ddx = ddx - box
            ddy = abs(y - y6(j))
            if (ddy > box2) ddy = ddy - box
            ddz = abs(z - z6(j))
            if (ddz > box2) ddz = ddz - box
            rw2(j) = ddx * ddx + ddy * ddy + ddz * ddz
         end do

         irsum = 0

         do j = 1, nspec
            jm = mp * nspec + j
            irej(jm) = 0

            do k = 1, npart
               if (rw2(k) < hc2v(k, j)) irej(jm) = 1
            end do

            irsum = irsum + irej(jm)
         end do

         if (irsum == nspec) then
            ihcall(mp) = ihcall(mp) + 1
            go to 110
         end if

         do k = 1, npart
            rwi(k) = 1.0 / sqrt(rw2(k))
         end do

         wtot2 = 0
         wtot3 = 0
         nl    = 1

         do j = 1, nspec
            nh   = int(hion(j, 2))
            wsum = 0

            do k = nl, nl + nh - 1
               wsum = wsum + rwi(k)
            end do

            nl = nl + nh
            wtot2 = wtot2 + wsum * hion(j, 4)
            wtot3 = wtot3 + wsum
         end do

         wtot2 = wtot2 + uj1

         do j = 1, nspec
            jm = mp * nspec + j

            if (irej(jm) == 1) then
               ihc(jm) = ihc(jm) + 1
               go to 160
            end if

            expuw(jm) = expuw(jm) + exp(abeta * wtot2 * hion(j, 4) * ecf)

            do k1 = 0, 10
               k   = k1 + 1
               ew  = hion(j, 4) * (wtot2 - k1 * 0.1 * hion(j, 4) * wtot3 / npart)
               ewla = ew * k1 * 0.1
               ewd  = exp(abeta * ecf * ewla)
               ewden(jm, k) = ewden(jm, k) + ewd
               ewnom(jm, k) = ewnom(jm, k) - ew * abeta * ecf * ewd
            end do

160         continue
         end do

110      continue
      end do
   end do

   return


   entry widom2

   ntocp = nspec * (nfix + 1)

   do i = 1, ntocp
      do j = 1, 11
         if (ewden(i, j) == 0) then
            write(jjj, *) ' WIDOM DENOMINATOR EQUALS ZERO', i, j
         else
            chint(i, j)  = ewnom(i, j) / ewden(i, j)
            ewnom(i, j)  = 0
            ewden(i, j)  = 0
            chinta(i, j) = chinta(i, j) + 1.0 / ny3 * chint(i, j)
         end if
      end do

      aint4 = chint(i, 2) + chint(i, 4) + chint(i, 6) + chint(i, 8) + chint(i, 10)
      aint2 = chint(i, 3) + chint(i, 5) + chint(i, 7) + chint(i, 9)
      aint1 = chint(i, 1) + chint(i, 11)
      chel(my3, i) = 1.0 / 30.0 * (aint1 + 2 * aint2 + 4 * aint4)
   end do

   nwtot = mwcn(my3) * nwins
   nwtot = nwins * int(ny1 / nwint) * ny2

   do i = 1, ntocp
      ihc(i) = ihc(i) + ihcall(int((i - 1) / nspec))
      chhc(my3, i)  = -dlog(dble(nwtot - ihc(i)) / nwtot)
      chexw(my3, i) = -dlog(expuw(i) / nwtot)
      chex(my3, i)  = chel(my3, i) + chhc(my3, i)
      chto(my3, i)  = chex(my3, i) + chid(i)
      expuw(i) = 0
      ihc(i)   = 0
   end do

   do j = 1, nspec
      dch1(my3, j) = chto(my3, j) - chex(my3, nspec + j)
      dch2(my3, j) = chto(my3, j) - chex(my3, 2 * nspec + j)
   end do

   write(kkk, 2010) nwtot

   do j = 0, nfix
      ihcall(j) = 0
      write(kkk, '(/, a)') str(j + 1)
      write(kkk, 2015)
      write(kkk, 2020) (i, chid(i), chhc(my3, i), chel(my3, i), &
                        chex(my3, i), chto(my3, i), chexw(my3, i), &
                        i = j * nspec + 1, (j + 1) * nspec)
   end do

   return


   entry widom3

   nwtot = 0

   do i = 1, ny3
      nwtot = nwtot + mwcn(i) * nwins
   end do

   ntocp = nspec * (nfix + 1)
   i = 1

   call earth2(chel,  chelv,  chelav,  ny3, ntocp)
   call earth2(chexw, chexwv, chexwa,  ny3, ntocp)
   call earth2(chhc,  chhcv,  chhcav,  ny3, ntocp)
   call earth2(chex,  chexv,  chexav,  ny3, ntocp)
   call earth2(chto,  chtov,  chtoav,  ny3, ntocp)

   write (jjj, '(a, /)') 'CONTACT CORRELATION g(r)'
   write (jjj, 2001) (i, i = 1, nspec)

   pcollav = 0.0
   pcollv  = 0.0

   do i = 1, nspec
      do k = 1, nspec
         do j = 1, ny3
            dum(j) = cwi(j, i, k) * exp(chexw(j, k))
         end do

         call earth(dum, cwiv, cwiav, ny3)
         cwi(11, i, k) = cwiv
         cwi(12, i, k) = cwiav
      end do

      write(jjj, 2002) i, (cwi(12, i, k), cwi(11, i, k), k = 1, nspec)
   end do

   write (jjj, '(/, a, /)') 'CONTACT PRESSURE MATRIX'
   write (jjj, 2001) (i, i = 1, nspec)

   do i = 1, nspec
      do k = 1, nspec
         do j = 1, ny3
            dum(j) = 2.0 / 3.0 * pi * caver(i) * cwi(j, i, k) * &
                     exp(chexw(j, k) + chid(k)) * (hion(i, 3) + hion(k, 3)) ** 3
         end do

         call earth(dum, cwiv, cwiav, ny3)
         cwi(11, i, k) = cwiv
         cwi(12, i, k) = cwiav
         pcollav = pcollav + cwiav
         pcollv  = pcollv  + cwiv * cwiv
      end do

      write(jjj, 2002) i, (cwi(12, i, k), cwi(11, i, k), k = 1, nspec)
   end do

   pcollv = sqrt(pcollv)

   write(jjj, '(/, a, f12.5, f10.5, /)') 'Total Collision pressure   =        &
                                          ', pcollav, pcollv


   write(jjj, 44)
   write(jjj, 2034)
   write(jjj, 2035) ((i - 1.0) / 10.0, i = 1, 11)

   do i = 1, ntocp
      write(jjj, 2040) i, (chinta(i, j), j = 1, 11)
   end do

   write(jjj, 2010) nwtot

   do j = 0, nfix
      write(jjj, '(/, a)') str(j + 1)
      write(jjj, 2015)
      write(jjj, 2030) (mod(i - 1, nspec) + 1, chid(i), chhcav(i), chhcv(i) &
                        , chelav(i), chelv(i), chexav(i), chexv(i), chtoav(i), chtov(i) &
                        , chexwa(i), chexwv(i), i = j * nspec + 1, nspec * (j + 1))
      write(jjj, 2033) ( exp((2 * chexav(i - 2) + chexav(i - 1)) / 3) )
   end do


2001 format('Species      ', 10(i12, 12x))
2002 format(i4, '        ', 10(f12.5, f10.5))
2010 format(/, /, 'CHEMICAL POTENTIALS IN UNITS OF KT, MEASURED', i8, ' TIMES')
2015 format(/, 'Type  Ideal     Hard Core       El. Static       Excess', &
           '         Total         Uncorrected Excess')
2020 format(10(i3, f8.4, 5(f9.4, '       '), /))
2030 format(10(i3, f8.4, 5(f9.4, f7.4), /))
2033 format('2:1 ', 10(f8.4), /)
2034 format(/, /, 'INTEGRAND IN WIDOM CORRECTION FOR DIFFERENT CHARGE PARAMETERS')
2035 format('TYPE', 11('  ', f4.2, '  '), /)
2040 format(10(i3, 11f8.3), /)
44   format(/, '************************************************************', &
           '**********************************************')

   return

end subroutine widom
