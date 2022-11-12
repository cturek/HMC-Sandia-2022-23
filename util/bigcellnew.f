      PROGRAM BIGCELL
c
c###############################################################
c
c Program: BIGCELL
c
c Written: Peter A. Schultz, January 1999
c
c Purpose: generate NxMxL larger unit cells from basis cell
c
c Revision history:
c   3Aug04-PAS/1.3: bias distance sort, add greeting
c  11Jan04-PAS/1.2: offsets
c  19Nov03-PAS/1.1: general Quest geom reader
c  15Sep03-PAS/1.0a: fix bug in rprim expansion (thx A. Thompson)
c  18Nov02-PAS/1.0: add sorting options, input general atomline
c  19Jul00-PAS: add skewed cells, origins
c  21Apr00-PAS: add dimensionality
c###############################################################
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c natmd=maximum number of total atoms
      PARAMETER  (natmd=1000)
c
      CHARACTER  version*5,         lastup*11
      DATA       version /'1.3  '/, lastup/' 5Aug04-PAS'/
c
      DIMENSION  rprim(3,3), rprimb(3,3)
      DIMENSION  ratmuc(3,natmd),itypuc(natmd)
      DIMENSION  ratm(3,natmd),itypa(natmd)
      DIMENSION  rsqat(natmd), rsqbat(natmd), rbias(3)
      DIMENSION  r1(3), nlat(3), rlat(3,27)
      DIMENSION  rscale(3)
      DIMENSION  orig(3),offset(3),rvec(3),rvecb(3)
      DIMENSION  nlsc(3,3)
      CHARACTER  filnm*80, line*80
      CHARACTER  opt*1
      LOGICAL    incell
c
      PARAMETER  (natmnm=12)
      PARAMETER  (ntypd=10)
      CHARACTER*(24)  typnm(ntypd), atmnm(natmd), atmnmt
      CHARACTER*(24)  atmnmuc(natmd)
c
      DATA  IRD,IWR, IRV, IDT / 5,6, 7,8 /
      DATA  zero,one / 0.d0,1.d0 /
      DATA  rbias / 0.000035d0, 0.000027d0, 0.000021d0 /
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,9000)  version, lastup
 9000 format(
     $ 1x,'Program BIGCELL: generates supercells'
     $/1x,'Version:  ',a,'     Last revised: ',a
     $ /)
c Initialize origin, offset to zero's
      ioffset = 0
      do  i=1,3
        orig(i) = zero
        offset(i) = zero
        rscale(i) = one
      enddo
c
      write(IWR,*) 'Enter dimensionality:'
      read(IRD,*)  ndim
c
      write(IWR,*) 'Enter lattice vectors of basis cell:'
      read(IRD,*)  rprim
c
      write(IWR,*) 'Enter number of atoms in basis cell:'
      read(IRD,*)  natmuc
      if( natmuc.eq.0 ) goto 999
c
c How do we default order atoms? 1=block,2=interleave,3=R(atm) ...
      iatmorder = 3
c  ... and do we project them into the principal unit cell?
      incell = .false.
c
c Ask user for geometry file name:
      IRDR = IRD
      write(IWR,*) 'Enter geometry file name (i=direct input):'
      read(IRD,'(a)') filnm
      OPEN(unit=IRV,file=filnm,status='old',form='formatted',err=30)
      IRDR = IRV
c
c     Determine format of file from lead line:
      read(IRDR,'(a80)') line
      write(IWR,'(a80)') line
      call ATFMTIN( iatmfmt, line )
      write(IWR,*) 'format style=',iatmfmt
c     Tell it to pick up atom types without checking if have it:
      if( iatmfmt.eq.2 ) iatmfmt = 3
      goto 40
c
   30 continue
      write(IWR,*) 'Enter {atom type x y z} for atoms in basis cell:'
      iatmfmt = 3
c
   40 continue
c Use Quest configuration reader to extract data:
c  turn off write:
      IWRX = -1
      call CONFIGRD( IRDR,IWRX, iatmfmt, natmnm, ntyp,typnm,
     $ natmuc, atmnmuc, itypuc, ratmuc )
      if( iatmfmt .eq. 3 ) iatmfmt = 2
C     do  iatm=1,natmuc
C       read(IRDR,*)  jatm,itypuc(iatm),(ratmuc(i,iatm),i=1,3)
C     enddo
c
      if( IRDR.eq.IRV ) CLOSE(unit=IRDR)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Setup options
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Regurgitate defaults
      if( iatmorder.eq.1 )then
        write(IWR,*) 'Default atom order is by BLOCK'
      elseif( iatmorder.eq.2 )then
        write(IWR,*) 'Default atom order is INTERLEAVE'
      elseif( iatmorder.eq.3 )then
        write(IWR,*) 'Default atom order is distance to origin'
      else
        iatmorder = 3
      endif
      if( incell )then
        write(IWR,*) 'Projection into primary cell is ON/'
      else
        write(IWR,*) 'Projection into primary cell is OFF/'
      endif
c
      ilatx = 0
   50 continue
      write(IWR,*) 'Enter option: (g-o,s-ort,o-ffset,h-elp,c-enter):'
      read(IRD,'(a1)')  opt
      if(     opt .eq. 'g' .or. opt .eq. 'G' )then
        goto  90
      elseif( opt .eq. 'p' .or. opt .eq. 'P' )then
c       Toggle projection into cell
        incell = .not. incell
        if( incell )then
          write(IWR,*) 'Projection into cell is toggled ON/',incell
        else
          write(IWR,*) 'Projection into cell is toggled OFF/',incell
        endif
      elseif( opt .eq. 's' .or. opt .eq. 'S' )then
c       Specify atom ordering:
        write(IWR,*) 'Atom order (1=block,2=interleave,3=R(atm):'
        read(IRD,*) iopt
        if( iopt.gt.0 .and. iopt.le.3 ) iatmorder = iopt
      elseif( opt .eq. 'c' .or. opt .eq. 'C' )then
        write(IWR,*) 'Enter center for R-sort:'
        read(IWR,*)  orig
      elseif( opt .eq. 'o' .or. opt .eq. 'O' )then
        write(IWR,*) 'Enter position offset:'
        read(IWR,*)  offset
        ioffset = 1
      elseif( opt .eq. 'f' .or. opt .eq. 'F' )then
        write(IWR,*) 'Enter scaling factors:'
        read(IWR,*)  rscale
      elseif( opt .eq. 'l' .or. opt .eq. 'L' )then
        ilatx = 1
      else
        write(IWR,9050)
 9050   format(/1x,'Available options:'
     $   /5x,'h(elp)   = this message'
     $   /5x,'o(ffset) = offset coordinate'
     $   /5x,'c(enter) = set alternate center for R-sort'
     $   /5x,'l(v)     = don''t ask - lv units'
     $   /5x,'p(roj)   = project atom coords into cell?'
     $   /5x,'f(actors)= scaling *F*actors'
     $   /5x,'s(ort)   = specify sort of atom list -'
     $   /5x,'         (1) 1-N,1a-Na,...  [blocked by uc]'
     $   /5x,'         (2) 1,1a,...,2,2a,... [interleaved]'
     $   /5x,'         (3) R1<R2<R3... [distance from origin]'
     $   /)
      endif
      goto  50
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Completed setup, now go on to generating output
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   90 continue
c  First, install all scaling functions.
c     Install scaling into primitive cell vectors:
      do  j=1,ndim
        do  i=1,ndim
          rprim(i,j) = rscale(i)*rprim(i,j)
        enddo
      enddo
c     Install scaling into atomic positions:
      do  iatm=1,natmuc
        do  i=1,3
          ratmuc(i,iatm) = rscale(i)*ratmuc(i,iatm)
        enddo
      enddo
      do  i=1,3
        orig(i) = rscale(i)*orig(i)
        offset(i) = rscale(i)*offset(i)
      enddo
c
  100 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Select expanded cell to generate from primitive cell
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write(IWR,*) 'Enter nla,nlb,nlc (dims of big cell):'
      read(IRD,*)  nla,nlb,nlc
c
c Set non-periodic directions to single repeat
      if( ndim.lt.3 ) nlc = 1
      if( ndim.lt.2 ) nlb = 1
      if( ndim.lt.1 ) STOP 'ndim=0'
      write(IWR,'(a,3i4)') 'nla,nlb,nlc=',nla,nlb,nlc
      if( nla .le. 0 )then
        write(IWR,*) '>>>> nla.le.0 - all done'
        goto 999
      endif
c
cpas:start
      write(IWR,913) 'rprim=',(rprim(i,1),i=1,3)
      write(IWR,913) 'rprim=',(rprim(i,2),i=1,3)
      write(IWR,913) 'rprim=',(rprim(i,3),i=1,3)
      write(IWR,913) 'orig=',orig
      write(IWR,'(a,3i4)') 'nla,nlb,nlc=',nla,nlb,nlc
  913 format(1x,a,3f16.6)
C 914 format(1x,2i3,3f16.6)
C     do  iatm=1,natmuc
C       write(IWR,914) iatm,itypa(iatm),(ratmuc(i,iatm),i=1,3)
C     enddo
c
cpas:end
c
      nabig = nla*nlb*nlc*natmuc
      if( nabig.gt.natmd )then
        write(IWR,*) '***** ERROR: too many atoms'
        write(IWR,*) '>>>>> # atms,max=',nabig,natmd
        goto 100
      endif
c
      if( ilatx.eq.1 )then
        write(IWR,*) 'Enter supercell vectors, lv units:'
        write(IWR,*) '***** BE AFRAID, BE VERY AFRAID *****'
        read(IWR,*) nlsc
        write(IWR,*) '>>>> DISABLED - diagonal expansion ONLY'
      else
        do  j=1,3
          do  i=1,3
            nlsc(i,j) = 0
          enddo
        enddo
        nlsc(1,1) = nla
        nlsc(2,2) = nlb
        nlsc(3,3) = nlc
      endif
c
      write(IWR,*) 'Enter output file name:'
      read(IRD,'(a)')  filnm
c
      if( filnm(1:4).eq.'quit' ) goto 999
c
c Perform expansion based on input lattice expansion matrix
c  (error in original code discovered by A. Thompson
C      do  idim=1,3
C        do  j=1,3
C          rprimb(j,idim) = DBLE( nlsc(1,idim) )*rprim(j,1)
C     $                   + DBLE( nlsc(2,idim) )*rprim(j,2)
C     $                   + DBLE( nlsc(3,idim) )*rprim(j,3)
C        enddo
C      enddo
c The atom-coord code currently cannot use general form, so ...
c I will do a strictly diagonal expansion A'=lA,B'=mB,C'=nC
      do  idim=1,3
        do  j=1,3
          rprimb(j,idim) = DBLE( nlsc(idim,idim) ) * rprim(j,idim)
        enddo
      enddo
c
      nlat(1) = nla
      nlat(2) = nlb
      nlat(3) = nlc
c
c Copy primary cell over directly:
c
c If interleaved order specified, change stride:
      jdel = 1
      if( iatmorder.eq.2 ) jdel = nla*nlb*nlc
c
      jatm = 1
      do  iatm=1,natmuc
        atmnm(jatm) = atmnmuc(iatm)
        itypa(jatm) = itypuc(iatm)
        do  i=1,3
          ratm(i,jatm) = ratmuc(i,iatm)
        enddo
        jatm = jatm + jdel
      enddo
c
      natm = natmuc
      ncopys = 1
      jcell = 1
      do 300 idim=1,3
        jatm = natm
        r1(1) = zero
        r1(2) = zero
        r1(3) = zero
        ncopy = ncopys
        do 250 ilat=2,nlat(idim)
          r1(1) = r1(1) + rprim(1,idim)
          r1(2) = r1(2) + rprim(2,idim)
          r1(3) = r1(3) + rprim(3,idim)
c
          iatm = 1
          jatm = natm + 1
          do 210 icell=1,ncopy
            jcell = jcell + 1
            if( iatmorder.eq.2 )then
c             Interleaved, start at copy number
              iatm = icell
              jatm = jcell
            endif
            do 200 iatmuc=1,natmuc
              atmnm(jatm) = atmnm(iatm)
              itypa(jatm) = itypa(iatm)
              ratm(1,jatm) = ratm(1,iatm) + r1(1)
              ratm(2,jatm) = ratm(2,iatm) + r1(2)
              ratm(3,jatm) = ratm(3,iatm) + r1(3)
              jatm = jatm + jdel
              iatm = iatm + jdel
  200       continue
            natm = natm + natmuc
            ncopys = ncopys + 1
  210     continue 
c
  250   continue
c
  300 continue
c
c Have raw cell, now install desired offset
c
      if( ioffset.ne.0 )then
        write(IWR,*) 'Install offset=',offset
        do  iatm=1,natm
          ratm(1,iatm) = ratm(1,iatm) - offset(1)
          ratm(2,iatm) = ratm(2,iatm) - offset(2)
          ratm(3,iatm) = ratm(3,iatm) - offset(3)
        enddo
      endif
c
      if( iatmorder.eq.3 )then
c       We will be sorting atoms by distance from origin
c        so set up local lattice (search) vectors
c
c       Get nearby lattice vectors ...
c        ... first, restrict dimensionality
        jlc = 1
        jlb = 1
        jla = 1
        if( ndim.lt.3 )then
          jlc = 0
          if( ndim.lt.2 ) jlb = 0
        endif
c        ... and then create lattice vectors:
        ils = 1
        do 303 ilc=-jlc,jlc
          do 302 ilb=-jlb,jlb
            do 301 ila=-jla,jla
c             Put origin at beginning of list
              if( ila.eq.0 .and. ilb.eq.0 .and. ilc.eq.0 )then
                il = 1
              else
                ils = ils + 1
                il = ils
              endif
              do i=1,3
                rlat(i,il) = ila*rprimb(i,1) + ilb*rprimb(i,2)
     $                     + ilc*rprimb(i,3)
              enddo
  301       continue
  302     continue
  303   continue
      endif
c
c  Sort by distance from origin:
c
      do 400 iatm=1,natm
c
        if( incell .or. iatmorder.eq.3 )then
c
c         Project atom position onto cell vectors:
          call PROJLV( ratm(1,iatm),
     $     rprimb(1,1),rprimb(1,2),rprimb(1,3), r1 )
          edge = 0.0001d0
c
c         Put atom position inside primary cell:
          do 310 idim=1,ndim
            vchk = r1(idim) + edge
            ilat = ABS( vchk )
            if( r1(idim)+edge .lt. zero ) ilat = -ilat - 1
c
            if( ilat.ne.0 )then
              clat = DBLE( ilat )
              do  i=1,3
                ratm(i,iatm) = ratm(i,iatm) - clat*rprimb(i,idim)
              enddo
            endif
  310     continue
        endif
c
        if( iatmorder.eq.3 )then
c         Sort atoms in big cell by distance to origin
c
c         Find lattice image of atom iatm closest to origin:
c          (biased by small shift to get canonical order)
          do  i=1,3
            rvec(i) = ratm(i,iatm) - orig(i)
            rvecb(i) = ratm(i,iatm) - orig(i) - rbias(i)
          enddo
          rsqmn = zero
          rsqbmn = zero
          do  id=1,ndim
            rsqmn = rsqmn + rvec(id)**2
            rsqbmn = rsqbmn + rvecb(id)**2
          enddo
          ilmn = 1
          do  il=2,ils
            rsq = zero
            rsqb = zero
            do  id=1,ndim
              rsq = rsq + ( rvec(id) - rlat(id,il) )**2
              rsqb = rsqb + ( rvecb(id) - rlat(id,il) )**2
            enddo
            if( rsqb.lt.rsqbmn )then
              rsqmn = rsq
              rsqbmn = rsqb
              ilmn = il
            endif
          enddo
c
c         Shift atom coordinates to image closest to origin:
          ityp = itypa(iatm)
          do  i=1,3
            r1(i) = ratm(i,iatm) - rlat(i,ilmn)
          enddo
c
c         Insert atom iatm into list ordered by (biased) range:
          katm = iatm
          atmnmt = atmnm(iatm)
          do 350 jatm=iatm,2,-1
            if( rsqbmn.lt.rsqbat(jatm-1) )then
              atmnm(jatm) = atmnm(jatm-1)
              itypa(jatm) = itypa(jatm-1)
              ratm(1,jatm) = ratm(1,jatm-1)
              ratm(2,jatm) = ratm(2,jatm-1)
              ratm(3,jatm) = ratm(3,jatm-1)
              rsqat(jatm) = rsqat(jatm-1)
              rsqbat(jatm) = rsqbat(jatm-1)
              katm = jatm - 1
            else
              goto 351
            endif
  350     continue
  351     continue
          atmnm(katm) = atmnmt
          itypa(katm) = ityp
          ratm(1,katm) = r1(1)
          ratm(2,katm) = r1(2)
          ratm(3,katm) = r1(3)
          rsqat(katm) = rsqmn
          rsqbat(katm) = rsqbmn
c
c         End sorting this atom by distance to origin
        endif
c
c       end loop iatm=1,natm
  400 continue
c
      if( iatmorder .eq. 3 )then
        do  iatm=1,natm
          rsqat(iatm) = SQRT( rsqat(iatm) )
        enddo
        write(IWR,917)
        write(IWR,*) '>>>> R for range-sorted atoms:'
  915   format(a,5(i6,f8.3))
  916   format(5x,5(i7,f8.3))
        ncol = 5
        nline = ( natm + ncol - 1 ) / ncol
        do  iline=1,nline
          iatm1 = iline
C         write(IWR,915) 'range:', ( iatm,rsqat(iatm),
          write(IWR,916)  ( iatm,rsqat(iatm),
     $      iatm=iatm1,natm,nline )
        enddo
        write(IWR,917)
  917   format(1x, 72('-') )
      endif
c
  910 format(3f12.8)
  920 format(i3)
  930 format(i3,4x,i2,4x,4f12.8)
c
c Install scaling:
      do  j=1,ndim
        do  i=1,ndim
          rprimb(i,j) = rprimb(i,j) / rscale(i)
        enddo
      enddo
      do  iatm=1,natm
        do  i=1,3
          ratm(i,iatm) = ratm(i,iatm) / rscale(i)
        enddo
      enddo
c
      OPEN(unit=IDT,file=filnm,status='unknown',form='formatted')
      write(IDT,'(a,3i3)') 'primitive lattice vectors',nlat
      write(IDT,910)  rprimb
      write(IDT,'(a)') 'number of atoms in unit cell'
      write(IDT,920)  natm
      write(IDT,'(a)') 'atom, type, position vector'
C     do iatm=1,natm
C       write(IDT,930) iatm,itypa(iatm),(ratm(i,iatm),i=1,3)
C     enddo
c
c renumber the atoms (could make this an option instead ... )
        do  iatm=1,natm
          atmnm(iatm) = ' '
          write(atmnm(iatm),'(i6)') iatm
        enddo
c
      call CONFIGWR( IDT,      iatmfmt, natmnm, ntyp,typnm,
     $ natm  , atmnm  , itypa , ratm   )
      CLOSE(unit=IDT)
c
      goto 100
c
  999 STOP 'bye guy'
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PROJLV
c
c
      subroutine PROJLV( r, a,b,c, v )
c---------------------------------------------------------------
c Purpose: project vector r() onto linearly independent vectors
c          a(), b(), c(), putting result (Cramer's rule) in v().
c
c Written: Peter A. Schultz, 8-May-1991
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  r(3),a(3),b(3),c(3)
c output arrays:
      DIMENSION  v(3)
c
c >>>> EXECUTABLE CODE:
c
      denom =  a(1)*(b(2)*c(3)-b(3)*c(2)) +
     $         b(1)*(c(2)*a(3)-c(3)*a(2)) +
     $         c(1)*(a(2)*b(3)-a(3)*b(2))
c
      v(1) = ( r(1)*(b(2)*c(3)-b(3)*c(2)) +
     $         b(1)*(c(2)*r(3)-c(3)*r(2)) +
     $         c(1)*(r(2)*b(3)-r(3)*b(2)) ) / denom
      v(2) = ( a(1)*(r(2)*c(3)-r(3)*c(2)) +
     $         r(1)*(c(2)*a(3)-c(3)*a(2)) +
     $         c(1)*(a(2)*r(3)-a(3)*r(2)) ) / denom
      v(3) = ( a(1)*(b(2)*r(3)-b(3)*r(2)) +
     $         b(1)*(r(2)*a(3)-r(3)*a(2)) +
     $         r(1)*(a(2)*b(3)-a(3)*b(2)) ) / denom
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATFMT
c
c
      subroutine ATFMTIN( iatmfmt, line )
c---------------------------------------------------------------
c Purpose: connect type of atom-line format to atom-line label
c
c Written: Peter A. Schultz, 11-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c iatmfmt = format of atom description
c         = 0  "coordinates"
c         = 1  "atom, type, positions"
c
c Warning: routine assumes (does not check) that input "line"
c is at least six characters long, and that the output "line"
c has enough space to load the identifying label.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
      CHARACTER*(*)  readlbl*6
c
c >>>> EXECUTABLE CODE:
c
      iatmfmt = -1
c
      readlbl = line(1:6)
      if( readlbl.eq. 'coordi' )then
c       This will be position-only input (x,y,z)
        iatmfmt = 0
      elseif( readlbl.eq. 'atom, ' )then
c       This denotes atom,type,position (at,ty,x,y,z) input
        iatmfmt = 2
      endif
      RETURN
c
      entry ATFMTOUT( iatmfmt, line )
c
      if( iatmfmt.eq.0 )then
        line = 'coordinates (x y z) for atoms'
      elseif( iatmfmt.eq.1 .or. iatmfmt.eq.2 )then
        line = 'atom, type, position (at,ty,x,y,z) for atoms'
      else
C        call STOPXERR( 'ATFMTOUT: unknown atom output format' )
        write(*,*)     'ATFMTOUT: unknown atom output format'
        STOP 'atom-fmt'
      endif
      RETURN
c
c >>>> That's all, Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONFIGRD
c
c
      subroutine CONFIGRD( IDT,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c---------------------------------------------------------------
c Purpose: input atomic configuration
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  10Jan02-PAS/2.54: bare coord read enabled
c---------------------------------------------------------------
c
c iatmfmt = format of atom description
c         = 0  just x,y,z
c         = 1  atom, type(as number), coords
c         = 2  atom name, type(as name), coords
c         = 3  types previously unknown
c
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
      DIMENSION  itypa(*)
      DIMENSION  ratm(3,natm)
c
c >>>> EXECUTABLE CODE:
c
      do  iatm=1,natm
        call ATOMRD( IDT,IWR, iatmfmt, natmnm, ntyp, typnm,
     $   atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONFIGWR
c
c
      entry CONFIGWR( IDT, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c---------------------------------------------------------------
c Purpose: output atomic configuration
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  iatm=1,natm
        call ATOMWR( IDT, iatmfmt, natmnm, ntyp, typnm,
     $   atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COORDRD
c
c
      subroutine COORDRD( IDT, iatmfmt, natm,ratm )
c---------------------------------------------------------------
c Purpose: extract atomic coordinates from geometry file
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c atomic coordinates:
      DIMENSION  ratm(3,natm)
c Local
      DIMENSION  ratmi(3)
c
c >>>> EXECUTABLE CODE:
c
      jatmfmt = iatmfmt
      do  iatm=1,natm
        call RATMRD( IDT, jatmfmt, ratmi )
        ratm(1,iatm) = ratmi(1)
        ratm(2,iatm) = ratmi(2)
        ratm(3,iatm) = ratmi(3)
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COORDWR
c
c
      entry COORDWR( IDT,  natm,ratm )
c---------------------------------------------------------------
c Purpose: dump raw atomic coordinates to geometry file
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  iatm=1,natm
        ratmi(1) = ratm(1,iatm)
        ratmi(2) = ratm(2,iatm)
        ratmi(3) = ratm(3,iatm)
        call RATMWR( IDT, ratmi )
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMRD
c
c
      subroutine ATOMRD( IDT,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ atmnm, iatm,itypa, ratm )
c---------------------------------------------------------------
c Purpose: input atom coordinates
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Oh, the ineffable joys of character manipulation in fortran ...
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  itypa(*)
      DIMENSION  ratm(3)
c
c local declarations:
      PARAMETER  ( nline=80 )
      CHARACTER  line*(nline), strtmp*(nline)
      LOGICAL    printon
c
c >>>> EXECUTABLE CODE:
c
      printon = ( IWR .ge. 0 )
c Global format for output of atomic coordinates:
 9100 format(1x,3f16.10)
c
      if( iatmfmt .eq. 0 )then
c       Read the raw coordinates
        read(IDT,*)  ratm
        if( printon ) write(IWR,9100)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 .or.
     $        iatmfmt .eq. 3 )then
c
c       Read the line with the atomic information:
c
        read(IDT,'(a80)',end=1301,err=1301 )  line
        if( printon ) write(IWR,'(a)' )   line
        nl1 = 1
        nl2 = nline
c
c       Parse the atom number/name (first string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom name found' )
        if( n1 .eq. 0 ) STOP 'atm-name'
        nmlen = n2 - n1 + 1
C       if( nmlen .gt. natmnm ) call STOPXERR( 'atom name too long' )
        if( nmlen .gt. natmnm ) STOP 'longatom'
        strtmp = ' '
        strtmp(1:nmlen) = line(n1:n2)
        atmnm = strtmp(1:natmnm)
        nl1 = n2 + 1
c
c       First, check if atoms input as number sequence.
c       If they are, check that atoms are indeed in sequence.
        read( atmnm, * , err=120 ) jatm
C        if( jatm .ne. iatm ) call STOPXERR( 'atom out of sequence' )
        if( jatm .ne. iatm )  STOP 'atom-seq'
c       Shift atom name to right aligned 6-char deep
        atmnm = ' '
        write( atmnm, '(i6)' ) iatm
c
c       Parse the atom type:
c
  120   continue
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom type found on line' )
        if( n1 .eq. 0 )  STOP 'no type'
        nmlen = n2 - n1 + 1
C        if( nmlen .gt. natmnm ) call STOPXERR( 'type name too long' )
        if( nmlen .gt. natmnm ) STOP 'longtype'
        nl1 = n2 + 1
c
c       And now try and get the atom type:
c
        if( iatmfmt .eq. 2 .or. iatmfmt .eq. 3 )then
          itypa(iatm) = 0
          do  ityp=1,ntyp
            strtmp = typnm(ityp)
            call STRPARS( strtmp, 1,nline, nt1,nt2 )
            ntlen = nt2 - nt1 + 1
            if( ntlen .eq. nmlen .and.
     $          line(n1:n2) .eq. strtmp(nt1:nt2) )then
              itypa(iatm) = ityp
              goto 200
            endif
          enddo
        endif
        if( iatmfmt .eq. 3 )then
c         Declare a new type, and name it
c         **** WARNING **** ntyp dim not checked to be adequate
          ntyp = ntyp + 1
          write(*,9125) iatm,ntyp,line(n1:n2)
 9125     format(1x,'>>>> Atom#',i4,' - new type#',i2,' named >',a,'<')
          typnm(ntyp) = line(n1:n2)
          itypa(iatm) = ntyp
          goto 200
        endif
c       Did not find type in list of type names, perhaps as type #?
        read( line(n1:n2), * , err=130 ) ityp
C        if( ityp.lt.1 .or. ityp.gt.ntyp ) call STOPXERR( 'bad type' )
        if( ityp.lt.1 .or. ityp.gt.ntyp ) STOP 'bad type'
        itypa(iatm) = ityp
c
c       Atom type are input as numbers, not names:
        iatmfmt = 1
c
        goto 200
c
C  130   call STOPXERR( 'atom type unknown' )
  130   continue
        write(*,'(a,a)') 'type input=',line(n1:n2)
        write(*,*) 'number of types:',ntyp
        do  ityp=1,ntyp
           write(*,'(a,i2,1x,a,a)') 'type#',ityp,'=',line(n1:n2)
        enddo
        STOP 'unk-type'
c
  200   continue
c
c       Get the atom coordinates:
        read( line(nl1:nl2), * )  ratm
c
      else
c       Unknown geometry format
C        call STOPXERR( 'unknown geometry (input) format' )
        STOP 'geomform'
      endif
c
      RETURN
Cc
Cc >>>> ERROR HANDLING
Cc
C1301 call STOPXERR( 'error (or eof) reading atom coordinates' )
C     RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMWR
c
c
       entry ATOMWR( IDT, iatmfmt, natmnm, ntyp,typnm,
     $ atmnm, iatm,itypa, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( iatmfmt .eq. 0 )then
c       Write out raw coordinates without atomic id's
        write(IDT,9100)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c       Build line with the atomic information:
        line = ' '
        nl1 = 2
        nl2 = nline
c
c       Parse and incorporate the atom name:
        call STRPARS( atmnm, 1,natmnm, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'nameless atom in atom output' )
        if( n1 .eq. 0 ) STOP 'no_atmnm'
        nmlen = n2
        nl2 = nl1 + nmlen - 1
        line(nl1:nl2) = atmnm(1:n2)
        nl1 = nl1 + natmnm + 1
c
c       Parse and incorporate the atom type name:
        ityp = itypa(iatm)
        strtmp = ' '
        if( iatmfmt .eq. 1 )then
c         Output type as type #, not type name
          n2 = 6
          write( strtmp(1:6), '(i6)' ) ityp
        else
c         Output type as type name
          call STRPARS( typnm(ityp), 1,natmnm, n1,n2 )
C          if( n1 .eq. 0 ) call STOPXERR( 'nameless type in output' )
          if( n1 .eq. 0 ) STOP 'no_typnm'
          strtmp(1:natmnm) = typnm(ityp)
        endif
        nmlen = n2
        nl2 = nl1 + nmlen - 1
        line(nl1:nl2) = strtmp(1:n2)
        nl1 = nl1 + natmnm + 1
c
c       Incorporate the atom coordinates:
        nl2 = nl1 + 1 + 3*16
C        if( nl2 .gt. nline ) call STOPXERR( 'atom type/coord too long' )
        if( nl2 .gt. nline ) STOP 'longatll'
        write( line(nl1:nl2), 9100 )  ratm
c
        write(IDT,'(a)')  line(1:nl2)
c
      else
c       Unknown geometry format
C        call STOPXERR( 'unknown geometry (output) format' )
        STOP 'geom-fmt'
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RATMRD
c
c
      entry RATMRD( IDT, iatmfmt, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( iatmfmt .eq. 0 )then
c       Read the raw coordinates
        read(IDT,*)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c
c       Read the line with the atomic information:
c
        read(IDT,'(a80)',end=1301,err=1301 )  line
        nl1 = 1
        nl2 = nline
c
c       Skip the atom number/name (first string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom name found' )
        if( n1 .eq. 0 ) STOP 'no-atmnm'
        nl1 = n2 + 1
c
c       Skip the atom type (second string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom type found on line' )
        if( n1 .eq. 0 ) STOP 'no-typin'
        nl1 = n2 + 1
c
c       Get the atom coordinates:
        read( line(nl1:nl2), * )  ratm
c
      else
c       Unknown geometry format
C        call STOPXERR( 'unknown geometry (input) format' )
        STOP 'geom-fmt'
c
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RATMWR
c
c
      entry RATMWR( IDT, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      write(IDT,9100)  ratm
c
      RETURN
c
c >>>> ERROR HANDLING
c
C 1301 call STOPXERR( 'error (or eof) reading atom coordinates' )
 1301 STOP 'err-ratm'
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRPARS
c
c
      subroutine STRPARS( line, nl1,nl2, n1,n2 )
c---------------------------------------------------------------
c Purpose: parse bounds of first character string within line
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
c
c local declarations:
      CHARACTER*(1)  blank
      DATA           blank / ' ' /
c
c >>>> EXECUTABLE CODE:
c
      n1 = 0
      n2 = nl2
      do  il=nl1,nl2
        if( line(il:il) .ne. blank )then
c         Found the first non-blank character ...
          n1 = il
          n2 = n1
c          ... and find the last character in string:
          do  jl=n1+1,nl2
            if( line(jl:jl) .eq. blank ) goto 999
            n2 = jl
          enddo
          goto 999
        endif
      enddo
c
  999 continue
c
c    That's all Folks!
c
      RETURN
      END
