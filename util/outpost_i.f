      PROGRAM OUTPOST
c===========================================================
c Program OUTPOST
c Purpose: take data from Quest history file and
c          examine, output into another format.
c
c Written: Peter A. Schultz (PAS), 20-May-2003
c
c Contributors: Thomas Mattsson (TM)
c
c Revision history:
c  11Sep03-PAS: some cleanup
c  16Jul03-TM: JMOL-compatible output.
c===========================================================
c
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER  version*5,         lastup*11
      DATA       version /'2.56 '/, lastup/'16Jul03-TM '/
c
c *************** DIMENSIONING ***************
      INCLUDE    'lcao.params'
c
      COMMON     /SUMMSTUFF/ r_units
c Atom definition data arrays:
      PARAMETER  (natmnm=12)
c      CHARACTER*(natmnm)  typnm(ntypd)
      CHARACTER*12 typnm(ntypd)
      DIMENSION  itypa(natmd)
      DIMENSION  vatm(3,natmd), vec(3)
c
      CHARACTER  filenm*128, line*128,line12*12
      EQUIVALENCE  (line,line12)
      CHARACTER  opt*1
c
      LOGICAL  idtopen
      LOGICAL  interact
      DATA  IRD,IWR,IDT,IOUT / 5,6,7,8 /
c
c >>>> EXECUTABLE CODE:
c
c interactive or batch mode:
      interact = .true.
c
      if( interact )  write(IWR,9000)  version,lastup
 9000 format(1x,72('=')
     $ // 4x,'Program OUTPOST: Quest post-analysis code',
     $ // 4x,'Written by Peter A. Schultz, May-2003'
     $ // 4x,'       and Thomas Mattsson, July-2003'
     $ // 4x,'Version: ',a,'        Last revised: ',a
     $ // 1x,72('=') / )
c
c  Open the file, and do trivial processing.
c
   20 continue
      idtopen = .false.
      if( interact )then
        write(IWR,'(a)') 'Enter quest history file name'
        read(IRD,'(a)') filenm
        OPEN( IDT,file=filenm,status='old',form='formatted',err=1320)
        write(IWR,*) ' ... file successfully opened'
        idtopen = .true.
      else
        IDT = IRD
      endif
c
c  Get basic data from history file ...
c
      IWRX = -1
      call SUMM( IWRX, IDT, marker,
     $ natm,ntyp, itypa,natmnm,typnm,  ncell, ngeom, efinal, ierr )
      if( ierr.ne.0 )then
        write(IWR,*) 'Error processing Quest summary file'
        goto 999
      endif
c
c   ... and give a summary of what we found in the file
c
      if( interact )then
        write(IWR,'(a,i6)') 'Number of atoms found:',natm
        write(IWR,'(a,f20.8)') 'Last energy found (Ry)=',efinal
 9021   format(1x,'  out of',i4,' cells and',i5,' total geometries')
 9022   format(1x,'  out of',i5,' total geometries detected')
        if( ncell .gt. 0 )then
          write(IWR,9021)  ncell, ngeom
        else
          write(IWR,9022)  ngeom
        endif
      endif
c
c  Offer options to user:
c
  100 continue
      if( interact )
     $write(IWR,'(/a)')  'Enter option (h-elp,q-uit,m-olden,j-mol):'
      if( interact )then
        read(IRD,'(a1)')  opt
      else
        opt = 'm'
      endif
 9100 format(1x,'Available options:'
     $ /5x,'q = q(uit)'
     $ /5x,'h = h(elp)'
     $ /5x,'e = e(nergy) summary'
     $ /5x,'m = m(olden) output'
     $ /5x,'j = j(mol) output'
     $ /5x,'t = t(ype) names'
     $ / )
c
c  Execute selected option:
c
  101 continue
      if( opt .eq. 'q' .or. opt .eq. 'x' )then
        goto 999
      elseif( opt .eq. 'e' )then
c       Echo energy summary
        call SUMM( IWR, IDT, marker,
     $   natm,ntyp, itypa,natmnm,typnm,  ncell,ngeom, efinal, ierr )
        goto 100
      elseif( opt .eq. 'm' )then
c       Output Molden file geometry or cell.
        goto 200
      elseif( opt .eq. 'j' )then
c       Output jmol file geometry or cell.
        goto 250
      elseif( opt .eq. 'g' )then
c       Output Molden file based on geometries
        if( ngeom .eq. 0 )then
          write(IWR,*) 'No geometries available'
          goto 100
        endif
        goto 200
      elseif( opt .eq. 't' )then
c       Change atom type names
        goto 300
      else
c       Print the help message, and try again
        write(IWR,9100)
        goto 100
      endif
c
c ********** Output MOLDEN-specific files **********
c
  200 continue
c
      if( interact )then
        IWRM = IWR
        write(IWR,'(a)')  'Enter output file name:'
        read(IRD,'(a)')  filenm
        OPEN(IOUT,file=filenm,status='NEW',form='formatted',err=200)
      else
        IWRM = -1
        IOUT = IWR
      endif
c
      if( ncell.lt.1 .or. opt .eq. 'g' )then
        call MOLDG( IDT, IWRM, IOUT, marker,
     $   itypa, natmnm,typnm, natm,vatm, ierr )
      else
        call MOLDUC( IDT, IWRM, IOUT, marker,
     $   itypa, natmnm,typnm, natm,vatm, ierr )
      endif
c
      CLOSE( IOUT )
c
      if( interact ) goto 100
      goto 999
c
c
c ********** Output JMOL-specific files **********
c
  250 continue
c
      if( interact )then
        IWRM = IWR
        write(IWR,'(a)')  'Enter output file name:'
        read(IRD,'(a)')  filenm
        OPEN(IOUT,file=filenm,status='NEW',form='formatted',err=200)
      else
        IWRM = -1
        IOUT = IWR
      endif
c
      if( ncell.lt.1 .or. opt .eq. 'g' )then
        call JMOLOUT( IDT, IWRM, IOUT, marker,
     $   itypa, natmnm,typnm, natm,vatm, ierr )
      else
        call MOLDUC( IDT, IWRM, IOUT, marker,
     $   itypa, natmnm,typnm, natm,vatm, ierr )
      endif
c
      CLOSE( IOUT )
c
      if( interact ) goto 100
      goto 999
c
c ********** Offer user chance to change type names **********
c
  300 continue
      write(IWR,'(a)') '***** Change type names:'
      do  ityp=1,ntyp
        write(IWR,'(/a,i3,a)')  'Type #',ityp,' name=',typnm(ityp)
        write(IWR,'(a)')  '  Enter new type name (blank=keep):'
        read(IRD,'(a)')  line
        if( line12 .ne. ' ' )then
          typnm(ityp) = line(1:natmnm)
          write(IWR,'(a,a)')  '  New type name=',typnm(ityp)
        endif
      enddo
c
      goto 100
c
c ********** Error handling **********
c
 1320 continue
      write(IWR,*) '**** ERROR: could not open file'
      if( interact ) goto 20
      goto 999
c
c ********** Good-bye **********
c
  999 continue
      if( idtopen ) CLOSE( unit=IDT )
c
      STOP
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUMM
c
c
      subroutine SUMM( IWR, IDT, marker,
     $ natm,ntyp, itypa,natmnm,typnm,  ncell,ngeom, efinal, ierr )
c---------------------------------------------------------------
c Purpose: extract energy summary from Quest summary file.
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON     /SUMMSTUFF/ r_units
c
      DIMENSION  itypa(*)
      CHARACTER*12  typnm(7)
c      CHARACTER*(natmnm)  typnm(*)
      LOGICAL  printon
      CHARACTER  line*128, readlbl*4, dist*1
      CHARACTER  line1*1
      EQUIVALENCE  (line,line1)
      CHARACTER  marker*1
c
      DATA  one / 1.d0 /
      DATA  au_to_A / 0.5291771d0 /
c
c >>>> EXECUTABLE CODE:
c
      ierr = 0
      printon = .false.
      if( IWR.gt.0 )  printon = .true.
c
c Default is bohr units:
      r_units = au_to_A
c
      REWIND( unit=IDT )
c
c Identify key-character for tag
      read(IDT,9101)  line
      if( line(2:9) .eq. '=KEYCHAR' )then
        marker = line1
      else
        ierr = 1
        RETURN
      endif
c
 9120 format(1x,a,i6)
 9125 format(1x,'Geometry #',i5,' energy=',f20.10)
 9126 format(1x,'Unit cell #',i4, 'energy=',f20.10)
 9129 format(1x,'FINAL LISTED ENERGY=',f20.10)
c
      efinal = 0
      natm = 0
      ncell = 0
      ngeom = 0
  100 continue
      read(IDT,9101,end=900)  line
 9101 format(a128)
      if( line1 .ne. marker ) goto 100
      readlbl = line(2:5)
c
      if( readlbl .eq. 'NUMB' )then
c       Get number of atoms, and their flavors
        read(IDT,*)  natm
        if( printon )  write(IWR,9120)  'Number of atoms=',natm
        do  iatm0=0,natm-1,20
          iatm1 = iatm0 + 1
          iatms = iatm0 + 20
          if( natm.lt. (iatm0+20) ) iatms = natm
          read(IDT,*)  (itypa(iatm),iatm=iatm1,iatms)
        enddo
c
      elseif( readlbl .eq. 'DIST' )then
c       Get the units used for distances in the file
        read(IDT,9101)  line
        call STRPARS( line, 1,128, ic1,icn )
        dist = line(ic1:ic1)
        if( dist.eq.'a' .or. dist.eq.'A' )then
          r_units = one
          if( printon )  write(IWR,9120)  'Distance unit=A'
        elseif( dist.eq.'b' .or. dist.eq.'B' )then
          r_units = au_to_A
          if( printon )  write(IWR,9120)  'Distance unit=bohr'
        else
          write(*,*) 'Distance unit unknown:',line(ic1:icn)
          STOP
        endif
c
      elseif( readlbl .eq. 'TYPE' )then
c       Get the number of types, and the names of those types
        read(IDT,*)  ntyp
        if( printon )  write(IWR,9120)  'Number of types=',ntyp
        do  ityp=1,ntyp
          read(IDT,9101)  line
          call STRPARS( line, 1,128, ic1,icn )
          if( icn .gt. (ic1+natmnm-1) ) icn = (ic1+natmnm-1)
          typnm(ityp) = line(ic1:icn)
        enddo
c
      elseif( readlbl .eq. 'ESCF' )then
c       Get the scf energy for a given geometry
        read(IDT,*)  efinal
        ngeom = ngeom + 1
        if( printon )  write(IWR,9125)  ngeom,efinal
c
      elseif( readlbl .eq. 'ECEL' )then
c       Get the energy of a relaxed cell
        read(IDT,*)  efinal
        ncell = ncell + 1
        if( printon )  write(IWR,9126)  ncell,efinal
c
      elseif( readlbl .eq. 'END ' )then
c       Found the end of the file
        if( printon )  write(IWR,'(a)')  'Found end of data'
        goto 900
c
      endif
      goto 100
c
  900 continue
      if( natm.le.0 .or. (ngeom.eq.0 .and. ncell.eq.0) ) ierr = 1
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GREAD
c
c
      subroutine GREAD( IDT, natm,vatm )
c---------------------------------------------------------------
c Purpose: retrieve atomic coordinates from file
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON     /SUMMSTUFF/ r_units
c
      DIMENSION  vatm(3,*)
      DIMENSION  vec(3)
c
c >>>> EXECUTABLE CODE:
c
        do  iatm=1,natm
c         Retrieve coordinate vector from file, convert to Angstrom
          read(IDT,*)  vec
          vatm(1,iatm) = vec(1) * r_units
          vatm(2,iatm) = vec(2) * r_units
          vatm(3,iatm) = vec(3) * r_units
        enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MOLGWRIT
c
c
      subroutine MOLGWRIT( imoldfl, itypa, natmnm,typnm, natm,vatm )
c---------------------------------------------------------------
c Purpose: write out geometry in molden format
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c      CHARACTER*(natmnm)  typnm(*)
      CHARACTER*12 typnm(*)
      DIMENSION  itypa(*), vatm(3,*)
c
      CHARACTER  line*128
c
c >>>> EXECUTABLE CODE:
c
Cc Heading of section (not used)
C 9040 format('[GEOMETRIES]')
c Number of atoms
 9041 format( 1x, i4 )
c Blank line after number of atoms
 9042 format( 1x, ' ' )
c Format for atom-type, and atom-positions
 9043 format( 1x, a6, 3(5x,e20.10) )
c
c       Write out heading for geometry ...
        write(imoldfl,9041)  natm
        write(imoldfl,9042)
        do  iatm=1,natm
c         Construct atom type name
c          (Limit six characters to accommodate MOLDEN)
          ityp = itypa(iatm)
          line = typnm(ityp)
          call STRPARS( line, 1,natmnm, i1,i2 )
          if( i2 .gt. (i1+5) ) i2 = i1 + 5
c
c         Write atom line to Molden file
          write(imoldfl,9043)  line(i1:i2), (vatm(i,iatm),i=1,3)
        enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MOLGWRIT
c
c
      subroutine JMOLWRITE( imoldfl, itypa, natmnm,typnm, natm,vatm,
     &     imdstep, eneout )
c---------------------------------------------------------------
c Purpose: write out geometry in format for jmol
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c Thomas R Mattsson July 2003, adapting with energy for each step.
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c      CHARACTER*(natmnm)  typnm(*)
      CHARACTER*12 typnm(*) 
      DIMENSION  itypa(*), vatm(3,*)
c
      CHARACTER  line*128
c
c >>>> EXECUTABLE CODE:
c
Cc Heading of section (not used)
C 9040 format('[GEOMETRIES]')
c Number of atoms
 9041 format( 1x, i4 )
c Energy on the line after number of atoms.
 9044 format( 1x, i4, 1x, e20.10 )
c Format for atom-type, and atom-positions
 9043 format( 1x, a6, 3(5x,e20.10) )
c
c       Write out heading for geometry ...
        write(imoldfl,9041)  natm
        write(imoldfl,9044)  imdstep, eneout
        do  iatm=1,natm
c         Construct atom type name
c          (Limit six characters to accommodate MOLDEN)
          ityp = itypa(iatm)
          line = typnm(ityp)
          call STRPARS( line, 1,natmnm, i1,i2 )
          if( i2 .gt. (i1+5) ) i2 = i1 + 5
c
c         Write atom line to Molden file
          write(imoldfl,9043)  line(i1:i2), (vatm(i,iatm),i=1,3)
        enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>JMOLWRITE
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>JMOLOUT
c
c
      subroutine JMOLOUT( IDT,IWR, imoldfl, marker,
     $ itypa, natmnm,typnm, natm,vatm, ierr )
c---------------------------------------------------------------
c Purpose: assemble JMOL output file.
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c Thomas Mattsson July-11-2003.
c adaptation for JMOL format.
c---------------------------------------------------------------
c
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*12 typnm(12)
c      CHARACTER*(natmnm)  typnm(*)
      DIMENSION  itypa(*), vatm(3,*), energy(1000)
c
      CHARACTER  line*128, readlbl*4, wantlbl*4
      CHARACTER  line1*1
      EQUIVALENCE  (line,line1)
      CHARACTER  marker*1
      LOGICAL    printon
c
c >>>> EXECUTABLE CODE:
c
      printon = .false.
      if( IWR.gt.0 ) printon = .true.
      REWIND( unit=IDT )
c
      ierr = 0
c File heading
 9000 format(a)
Cc Headings of sections (not used)
C 9040 format('[GEOMETRIES]')
c Number of atoms
 9041 format( 1x, i4 )
c Blank line after number of atoms
 9042 format( 1x, ' ' )
c Format for atom-type, and atom-positions
 9043 format( 1x, a6, 3(5x,e20.10) )
c
      mode = -1
   90 continue
      mode = mode + 1
      npoint = 0
      if( mode.eq.0 )then
c       Apply global heading for file
c        write(imoldfl,9000) '[Molden Format]'
c        write(imoldfl,9000) '[GEOCONV]'
        goto 90
c
      elseif( mode.eq.1 )then
        wantlbl = 'ESCF'
c        write(imoldfl,9000)  'energy'
c
      elseif( mode.eq.2 )then
        wantlbl = 'FMAX'
c        write(imoldfl,9000)  'max-force'
c
      elseif( mode.eq.3 )then
        wantlbl = 'FRMS'
c        write(imoldfl,9000)  'rms-force'
c
      elseif( mode.eq.4 )then
        wantlbl = 'COOR'
c        write(imoldfl,9000)  '[GEOMETRIES]'
c
      else
c       End of the line for this file
        goto 999
      endif
      REWIND( unit=IDT )
c
c Process quest summary file:
c
  100 continue
      read(IDT,9101,end=900)  line
 9101 format(a128)
      if( line1 .ne. marker ) goto 100
      readlbl = line(2:5)
c
      if( readlbl .eq. 'ECEL' .or. readlbl .eq. 'END ' )then
c       This is signal to stop looking for this sequence
        goto 900
      elseif( readlbl .eq. wantlbl )then
c       Found the keyword we wanted
        goto 200
      endif
      goto 100
c
  200 continue
c
      if( wantlbl .eq. 'COOR' )then
c       We found a geometry, retrieve it ...
        call GREAD( IDT, natm,vatm )
c        ... and write it out in jmol format, including the frame and energy

        eprt = energy(npoint)
        if( npoint.lt.npointmax) then 
        call JMOLWRITE( imoldfl, itypa, natmnm,typnm, natm,vatm,
     $       (npoint+1), eprt)
        endif
        npoint = npoint + 1
c
      elseif( wantlbl .eq. 'ESCF'
     $   .or. wantlbl .eq. 'FRMS'
     $   .or. wantlbl .eq. 'FMAX' )then
c       Transfer a scalar quantity (energy, f-max, f-rms)
        read(IDT,*)  scalar
        if( wantlbl .eq. 'ESCF' ) then
           energy(npoint) = scalar
        endif
c        write(imoldfl,*)  scalar
c        write(imoldfl,*)  energy(npoint)
        npoint = npoint + 1
c
      else
c       Coding screwup
        write(IWR,*) '***** ERROR - unknown key *****', wantlbl
        STOP 'MOLDG'
      endif
c
      goto 100
c
  900 continue
      if( printon )
     $write(IWR,'(a,a,i5)') wantlbl,': number of points found=',npoint
      if(wantlbl.eq.'ESCF')then
         npointmax = npoint
      endif
      goto 90
c
  999 continue
      if( printon )
     $write(IWR,*) 'JMOL file successfully completed'

      RETURN
c
c    That's all Folks!
c
      END
c
c
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MOLDG
c
c
      subroutine MOLDG( IDT,IWR, imoldfl, marker,
     $ itypa, natmnm,typnm, natm,vatm, ierr )
c---------------------------------------------------------------
c Purpose: assemble Molden output file.
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*12 typnm(12)
c      CHARACTER*(natmnm)  typnm(*)
      DIMENSION  itypa(*), vatm(3,*)
c
      CHARACTER  line*128, readlbl*4, wantlbl*4
      CHARACTER  line1*1
      EQUIVALENCE  (line,line1)
      CHARACTER  marker*1
      LOGICAL    printon
c
c >>>> EXECUTABLE CODE:
c
      printon = .false.
      if( IWR.gt.0 ) printon = .true.
      REWIND( unit=IDT )
c
      ierr = 0
c File heading
 9000 format(a)
Cc Headings of sections (not used)
C 9040 format('[GEOMETRIES]')
c Number of atoms
 9041 format( 1x, i4 )
c Blank line after number of atoms
 9042 format( 1x, ' ' )
c Format for atom-type, and atom-positions
 9043 format( 1x, a6, 3(5x,e20.10) )
c
      mode = -1
   90 continue
      mode = mode + 1
      npoint = 0
      if( mode.eq.0 )then
c       Apply global heading for file
        write(imoldfl,9000) '[Molden Format]'
        write(imoldfl,9000) '[GEOCONV]'
        goto 90
c
      elseif( mode.eq.1 )then
        wantlbl = 'ESCF'
        write(imoldfl,9000)  'energy'
c
      elseif( mode.eq.2 )then
        wantlbl = 'FMAX'
        write(imoldfl,9000)  'max-force'
c
      elseif( mode.eq.3 )then
        wantlbl = 'FRMS'
        write(imoldfl,9000)  'rms-force'
c
      elseif( mode.eq.4 )then
        wantlbl = 'COOR'
        write(imoldfl,9000)  '[GEOMETRIES]'
c
      else
c       End of the line for this file
        goto 999
      endif
      REWIND( unit=IDT )
c
c Process quest summary file:
c
  100 continue
      read(IDT,9101,end=900)  line
 9101 format(a128)
      if( line1 .ne. marker ) goto 100
      readlbl = line(2:5)
c
      if( readlbl .eq. 'ECEL' .or. readlbl .eq. 'END ' )then
c       This is signal to stop looking for this sequence
        goto 900
      elseif( readlbl .eq. wantlbl )then
c       Found the keyword we wanted
        goto 200
      endif
      goto 100
c
  200 continue
c
      if( wantlbl .eq. 'COOR' )then
c       We found a geometry, retrieve it ...
        call GREAD( IDT, natm,vatm )
c        ... and write it out in molden format:
        call MOLGWRIT( imoldfl, itypa, natmnm,typnm, natm,vatm )
        npoint = npoint + 1
c
      elseif( wantlbl .eq. 'ESCF'
     $   .or. wantlbl .eq. 'FRMS'
     $   .or. wantlbl .eq. 'FMAX' )then
c       Transfer a scalar quantity (energy, f-max, f-rms)
        read(IDT,*)  scalar
        write(imoldfl,*)  scalar
        npoint = npoint + 1
c
      else
c       Coding screwup
        write(IWR,*) '***** ERROR - unknown key *****', wantlbl
        STOP 'MOLDG'
      endif
c
      goto 100
c
  900 continue
      if( printon )
     $write(IWR,'(a,a,i5)') wantlbl,': number of points found=',npoint
      goto 90
c
  999 continue
      if( printon )
     $write(IWR,*) 'Molden file successfully completed'
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MOLDUC
c
c
      subroutine MOLDUC( IDT,IWR, imoldfl, marker,
     $ itypa, natmnm,typnm, natm,vatm, ierr )
c---------------------------------------------------------------
c Purpose: assemble Molden output file from cell optimization
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c      CHARACTER*(natmnm)  typnm(*)
      CHARACTER*12 typnm(*)
      DIMENSION  itypa(*), vatm(3,*)
c
      CHARACTER  line*128, readlbl*4, wantlbl*4
      CHARACTER  line1*1
      EQUIVALENCE  (line,line1)
      CHARACTER  marker*1
      LOGICAL    gotg
      LOGICAL    printon
c
c >>>> EXECUTABLE CODE:
c
      printon = .false.
      if( IWR.gt.0 ) printon = .true.
      REWIND( unit=IDT )
c
      ierr = 0
c File heading
 9000 format(a)
Cc Headings of sections (not used)
C 9040 format('[GEOMETRIES]')
c
      mode = -1
   90 continue
      mode = mode + 1
      npoint = 0
      if( mode.eq.0 )then
c       Apply global heading for file
        write(imoldfl,9000) '[Molden Format]'
        write(imoldfl,9000) '[GEOCONV]'
        goto 90
c
      elseif( mode.eq.1 )then
        wantlbl = 'ECEL'
        write(imoldfl,9000)  'energy'
c
      elseif( mode.eq.2 )then
        wantlbl = 'MAXS'
        write(imoldfl,9000)  'max-stress'
c
      elseif( mode.eq.3 )then
        wantlbl = 'RMSS'
        write(imoldfl,9000)  'rms-stress'
c
      elseif( mode.eq.4 )then
        wantlbl = 'COOR'
        write(imoldfl,9000)  '[GEOMETRIES]'
c
      else
c       Write out final geometry in cell, if we have it:
        if( gotg ) call MOLGWRIT( imoldfl,itypa,natmnm,typnm,natm,vatm )
c       End of the line for this file
        goto 999
      endif
      REWIND( unit=IDT )
      gotg = .false.
c
c Process quest summary file:
c
  100 continue
      read(IDT,9101,end=900)  line
 9101 format(a128)
      if( line1 .ne. marker ) goto 100
      readlbl = line(2:5)
c
      if( readlbl .eq. 'END ' )then
c       This is signal to stop looking for this sequence
        goto 900
      elseif( readlbl .eq. wantlbl )then
c       Found the keyword we wanted
        goto 200
      endif
      goto 100
c
  200 continue
c
      if( wantlbl .eq. 'COOR' )then
c       We found a geometry, retrieve it ...
        call GREAD( IDT, natm,vatm )
        gotg = .true.
        npoint = 1
c
      elseif( wantlbl .eq. 'ECEL'
     $   .or. wantlbl .eq. 'RMSS'
     $   .or. wantlbl .eq. 'MAXS' )then
c       Transfer a scalar quantity (energy, stress-max/rms)
        read(IDT,*)  scalar
        write(imoldfl,*)  scalar
        npoint = npoint + 1
c
      else
c       Coding screwup
        write(IWR,*) '***** ERROR - unknown key *****', wantlbl
        STOP 'MOLDG'
      endif
c
      goto 100
c
  900 continue
      if( printon )
     $write(IWR,'(a,a,i5)') wantlbl,': number of points found=',npoint
      goto 90
c
  999 continue
      if( printon )
     $write(IWR,*) 'Molden file successfully completed'
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
