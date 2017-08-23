      SUBROUTINE EXTINCT(L,B,D,AVT,SAVT,AVC,JMAX,AV,SAV)
C----------------------------------------------------------------------
C   Subroutine EXTINCT (v. 2.0.2, June 2, 1997)
C     By Jon Hakkila, Jeannette Myers, Brett Stidham, & Dieter Hartmann.
C     Calculates visual interstellar extinction due to the Milky Way.
C     The origins of this code are described by the authors in a
C     manuscript scheduled to appear in the Nov. 1997 Astronomical Journal.
C
C   Changes implemented since version 1.0 (Dec. 2, 1996):
C     Updated error analysis for the FitzGerald subroutine.
C     fixed a coding error in the FitzGerald subroutine that caused
C       extinction out of the Galactic plane to be overestimated at
C       distances larger than 1 kpc.
C     fixed cell boundaries & extinction in Neckel & Klare subroutine.
C     More accurate angular cloud boundaries for the high-Galactic
C       clouds analyzed in the high-latitude study (formerly the Penprase
C       subroutine). This routine now also identifies additional high-
C       Galactic latitude clouds.
C     Systematic underrepresentation of extinction is corrected by using
C       an additional correction term. 
C----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,L-Z)
      DIMENSION AV(5), SAV(5)
C----------------------------------------------------------------------
C   L      = Galactic longitude (degrees)    0 <= L <  360
C   B      = Galactic latitude (degrees)   -90 <= B <=  90
C   D      = Source distance (in kpc)        0 <= D
C   AVT    = Total visual extinction (magnitudes), using available subroutines.
C   SAVT   = Total extinction error (magnitudes), using available subroutines.
C   AVC    = Visual extinction correction (mag), using available subroutines.
C            The corrected extinction AVT+AVC removes suspected systematics.
C   JMAX   = Number of subroutines used to calculate extinction.
C   AV(i)  = Extinction calculated by the iTH subroutine.     
C   SAV(i) = Error calculated by the iTH subroutine.      
C     Specific subroutines used to calculate the extinction are returned
C     via AV(i): If AV(i) has a value not equal to -99. then the ith 
C     study has been used in calculating extinction.
C     i=1 is from Fitzgerald, AJ 73, 983 (1968).
C     i=2 is from Neckel and Klare, A&A Supp 42, 251 (1980).
C     i=3 is from Berdnikov & Pavlovskaya, Sov. Astron. Lett. 17, 215 (1990).
C     i=4 is from Arenou et al., A&A 258, 104 (1992).
C     i=5 is from Penprase, ApJ. Supp. 83, 273 (1992),
C                 Magnani et al., ApJ. 295, 402 (1985),
C                 Keto & Myers, ApJ. 304, 466 (1986),
C                 Desert et al., Ap.J. 334, 815 (1988),
C                 Odenwald, Ap.J. 325, 320 (1988),
C                 Hughes et al., A.J. 105, 571 (1993),
C                 Kenyon et al., A.J. 108, 1872 (1994),
C                 Cernis, ApSS 166, 315 (1990),
C                 Cernis, Baltic Astron. 2, 214 (1993),
C                 Rossano, AJ, 83, 234 (1978),
C                 Rossano, AJ, 83, 241 (1978),
C                 Kutner et al., ApJ, 215, 521 (1977).
C
C   All studies have been modified to statistically account for
C     unsampled regions.
C
C   A0 = 1.5 magnitudes per kpc
C   R  = A0/E(B-V) = 3.0
C
C----------------------------------------------------------------------
      A0 = 1.5
      R  = 3.1
C----------------------------------------------------------------------
C     If source is beyond Milky Way edge (assume Sun is 8.5 kpc from
C     Galactic Center and that Disk has a 15 kpc radius), then limit
C     extinction to what would be seen through Milky Way.
C----------------------------------------------------------------------
      XD=D      !Distance
      XL=L      !Galactic Longitude
      XB=B      !Galactic Latitude
      DSUN=8.5
      RGAL=15.
      T1=COS(XB*3.141592/180.)*COS(XL*3.141592/180.)
      T2=1.-(RGAL/DSUN)*(RGAL/DSUN)
      DMAX=8.5*(T1+SQRT(T1*T1-T2))
      IF (D .GT. DMAX) XD=DMAX
C----------------------------------------------------------------------
C     Calculate extinction from desired subroutines.
C----------------------------------------------------------------------
      CALL FITZGERALD(XL,XB,XD,AV(1),SAV(1),A0,R)
      CALL NECKEL_KLARE(XL,XB,XD,AV(2),SAV(2),A0)
      CALL BERDNIKOV(XL,XB,XD,AV(3),SAV(3),A0)
      CALL ARENOU_ETAL(XL,XB,XD,AV(4),SAV(4),A0)
      CALL HIGH_LAT(XL,XB,XD,AV(5),SAV(5),IHFLG)
C----------------------------------------------------------------------
C     Combine outputs from subroutines.
C----------------------------------------------------------------------
      CALL COMBINE(AV,SAV,AVT,SAVT,JMAX,IHFLG)
      CALL CORRECT(XL,XB,XD,AVC)
      RETURN
      END

      SUBROUTINE FITZGERALD(L,B,DIST,AV,SAV,A0,R)
C----------------------------------------------------------------------
C     From Fitzgerald, AJ 73, 983 (1968).
C----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,L-Z)
      IF (ABS(B) .GT. 0.) THEN
         D = DIST*COS(B*3.141592/180.)
      ELSE
         D = DIST
      ENDIF
      A0R= A0/R
      IF (L .LT. 0.) L=L+360.
      IFLG=0
      IF (B .NE. 0.) IFLG=1
      BETA0=.114
      IF ((L .GE. 10.) .AND. (L .LT. 60.)) THEN
      BETA=.045
      GOTO 660
      ENDIF
      IF ((L .GE. 60.) .AND. (L .LT. 100.)) THEN
      BETA=.040
      GOTO 661
      ENDIF
      IF ((L .GE. 100.) .AND. (L .LT. 130.)) THEN
      BETA=.055
      GOTO 662
      ENDIF
      IF ((L .GE. 130.) .AND. (L .LT. 160.)) THEN
      BETA=.170
      GOTO 663
      ENDIF
      IF ((L .GE. 160.) .AND. (L .LT. 190.)) THEN
      BETA=.105
      GOTO 664
      ENDIF
      IF ((L .GE. 190.) .AND. (L .LT. 220.)) THEN
      BETA=.070
      GOTO 665
      ENDIF
      IF ((L .GE. 220.) .AND. (L .LT. 250.)) THEN
      BETA=.065
      GOTO 666
      ENDIF
      IF ((L .GE. 250.) .AND. (L .LT. 280.)) THEN
      BETA=.080
      GOTO 667
      ENDIF
      IF ((L .GE. 280.) .AND. (L .LT. 310.)) THEN
      BETA=.040
      GOTO 668
      ENDIF
      IF ((L .GE. 310.) .AND. (L .LT. 340.)) THEN
      BETA=.060
      GOTO 669
      ENDIF
      IF ((L .GE. 340.) .AND. (L .LT. 360.)) THEN
      BETA=.100
      GOTO 670
      ENDIF
      IF ((L .GE. 0.) .AND. (L .LT. 10.)) THEN
      BETA=.100
      GOTO 659
      ENDIF
 659   IF ((L .GE. 0.0) .AND.  (L .LT. 4.0)) THEN
        IF ((D .GE. 0.0) .AND. (D .LT. 0.5)) THEN
          EY = 0.0
          GOTO 1000
        ENDIF
        IF ((D .GE. 0.5) .AND. (D .LT. 0.8)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.0 + 2.7*(D-0.5)      
           ELSE
            EY=EYFCN(2.7,0.5,D,B,BETA)
           ENDIF    
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.8) .AND.(D .LT. 4.2)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.81
         ELSE
            EY=EYFCN(2.7,0.5,0.8,B,BETA0)
           ENDIF
          GOTO 1000
          ENDIF
          IF ((D .GE. 4.2)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.81 + A0R*(D-4.2)
           ELSE
            EY=EYFCN(2.7,0.5,0.8,B,BETA0)+EYFCN(A0R,4.2,D,B,BETA0)
           ENDIF
          GOTO 1000
        ENDIF
       ENDIF
       IF ((L .GE. 4.0) .AND. (L .LT. 7.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
          EY = 0.0
          GOTO 1000
        ENDIF
        IF ((D .GE. 0.1) .AND. (D .LT. 0.5)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.0 + 0.75*(D-0.1)
           ELSE
            EY=EYFCN(0.75,0.1,D,B,BETA)
           ENDIF  
        GOTO 1000
        ENDIF
        IF ((D .GE. 0.5) .AND. (D .LT. 1.2)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.3
         ELSE
            EY=EYFCN(0.75,0.1,0.5,B,BETA)
           ENDIF  
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.2) .AND. (D .LT. 1.8)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.3 + 1.38*(D-1.2)
           ELSE
            EY=EYFCN(0.75,0.1,0.5,B,BETA)+EYFCN(1.38,1.2,D,B,BETA0)
           ENDIF
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.8)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 1.128 + A0R*(D-1.8)
           ELSE
            EY=EYFCN(0.75,0.1,0.5,B,BETA)+EYFCN(1.38,1.2,1.8,B,BETA0)
     &         +EYFCN(A0R,1.8,D,B,BETA0)
           ENDIF
        GOTO 1000
        ENDIF
       ENDIF
       IF ((L .GE. 7.0) .AND. (L. LT. 10.0)) THEN
        IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
          EY = 0.0
          GOTO 1000
        ENDIF
        IF ((D .GE. 0.2) .AND. (D .LT. 0.5)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.0 + 1.17*(D-0.2)
          ELSE
            EY=EYFCN(1.17,0.2,D,B,BETA)
           ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 0.5) .AND. (D .LT. 4.0)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.351
         ELSE
          EY=EYFCN(1.17,0.2,0.5,B,BETA)
           ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.351 + A0R*(D-4.0)
             ELSE
          EY=EYFCN(1.17,0.2,0.5,B,BETA)+EYFCN(A0R,4.0,D,B,BETA0)
         ENDIF
        GOTO 1000
        ENDIF
       ENDIF
 660   IF ((L .GE. 10.0) .AND. (L .LT. 12.0)) THEN
        IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
          EY = 0.0
          GOTO 1000
        ENDIF
        IF ((D .GE. 0.1) .AND. (D .LT. 0.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.0 + 1.43*(D-0.1)
         ELSE
          EY=EYFCN(1.43,0.1,D,B,BETA)
         ENDIF       
        GOTO 1000
        ENDIF
        IF ((D .GE. 0.45) .AND. (D .LT. 1.9)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.5005
         ELSE
          EY=EYFCN(1.43,0.1,0.45,B,BETA)
         ENDIF       
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.9) .AND. (D .LT. 2.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.5005 + 1.2488*(D-1.9)
         ELSE
          EY=EYFCN(1.43,0.1,0.45,B,BETA)+EYFCN(1.2488,1.9,D,B,BETA0)
         ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 2.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.00002 + A0R*(D-2.3)
         ELSE
          EY=EYFCN(1.43,0.1,0.45,B,BETA)+EYFCN(1.2488,1.9,2.3,B,BETA0)
     &         +EYFCN(A0R,2.3,D,B,BETA0)
         ENDIF
        GOTO 1000
        ENDIF
       ENDIF
       IF ((L .GE. 12.0) .AND. (L .LT. 14.0)) THEN
        IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN      
          EY = 0.0
          GOTO 1000
        ENDIF
        IF ((D .GE. 0.25) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.0 + 2.*(D-0.25)
         ELSE
          EY=EYFCN(2.,0.25,D,B,BETA)
         ENDIF    
        GOTO 1000
        ENDIF
        IF ((D .GE. 0.5) .AND. (D .LT. 1.9)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.5
         ELSE
          EY=EYFCN(2.,0.25,0.5,B,BETA)
         ENDIF    
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.9) .AND. (D .LT. 2.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.5 + 1.11*(D-1.9)
         ELSE
          EY=EYFCN(2.,0.25,0.5,B,BETA)+EYFCN(1.11,1.9,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9995 + A0R*(D-2.35)
         ELSE
          EY=EYFCN(2.,0.25,0.5,B,BETA)+EYFCN(1.11,1.9,2.35,B,BETA0)
     &        +EYFCN(A0R,2.35,D,B,BETA0)
         ENDIF
        GOTO 1000
        ENDIF
       ENDIF
       IF ((L .GE. 14.0) .AND. (L .LT. 16.0)) THEN
        IF ((D .GE. 0.0) .AND. (D .LT. 0.6)) THEN
          EY = 0.0
          GOTO 1000
        ENDIF
        IF ((D .GE. 0.6) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.0 + 1.5*(D-0.6)
         ELSE
          EY=EYFCN(1.5,0.6,D,B,BETA)
         ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.0) .AND. (D .LT. 2.5)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.6
         ELSE
          EY=EYFCN(1.5,0.6,1.,B,BETA)
         ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 2.5) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.6 + 1.2*(D-2.5)
           ELSE
          EY=EYFCN(1.5,0.6,1.0,B,BETA)+EYFCN(1.2,2.5,D,B,BETA0)
         ENDIF
        GOTO 1000
        ENDIF
        IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.2 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.5,0.6,1.0,B,BETA)+EYFCN(1.2,2.5,3.,B,BETA0)
     &        +EYFCN(A0R,3.0,D,B,BETA0)
         ENDIF
        GOTO 1000
              ENDIF
       ENDIF
       IF ((L .GE. 16.0) .AND. (L .LT. 17.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.67*(D-0.1)
         ELSE
          EY=EYFCN(1.67,0.1,D,B,BETA)
         ENDIF        
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.5)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY=0.2505
         ELSE
          EY=EYFCN(1.67,0.1,0.25,B,BETA)
         ENDIF        
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 0.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.2505 + 1.874*(D-0.5)
         ELSE
          EY=EYFCN(1.67,0.1,0.25,B,BETA)+EYFCN(1.874,0.5,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                
        IF ((D .GE. 0.9) .AND. (D .LT. 2.5)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 1.0001
         ELSE
          EY=EYFCN(1.67,0.1,0.25,B,BETA)+EYFCN(1.874,0.5,0.9,B,BETA)
         ENDIF      
        GOTO 1000
        ENDIF
          IF ((D .GE. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.0001 + A0R*(D-2.5)
         ELSE
          EY=EYFCN(1.67,0.1,0.25,B,BETA)+EYFCN(1.874,0.5,0.9,B,BETA)
     &        +EYFCN(A0R,2.5,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 17.0) .AND. (L .LT. 18.5)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 2.*D 
         ELSE
          EY=EYFCN(2.,0.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.75)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.4
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.75) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4 + 2.8*(D-0.75)
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)+EYFCN(2.8,0.75,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF        
          IF ((D .GE. 1.0) .AND. (D .LT. 3.1)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 1.1
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)+EYFCN(2.8,0.75,1.,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.1 + A0R*(D-3.1)
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)+EYFCN(2.8,0.75,1.,B,BETA)
     &        +EYFCN(A0R,3.1,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 18.5) .AND. (L .LT. 20.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 2.*D
         ELSE
          EY=EYFCN(2.,0.,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF               
          IF ((D .GE. 0.2) .AND. (D .LT. 0.75)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.4
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.75) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4 + 2.8*(D-0.75)
           ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)+EYFCN(2.8,0.75,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF  
          IF ((D .GE. 1.0) .AND. (D .LT. 3.0)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 1.1
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)+EYFCN(2.8,0.75,1.,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.1 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(2.,0.,0.2,B,BETA)+EYFCN(2.8,0.75,1.,B,BETA)
     &        +EYFCN(A0R,3.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 20.0) .AND. (L .LT. 24.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 1.5*D
         ELSE
          EY=EYFCN(1.5,0.,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3 + 0.67*(D-0.2)
         ELSE
          EY=EYFCN(1.5,0.,0.2,B,BETA)+EYFCN(0.67,0.2,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.35) .AND. (D .LT. 0.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4005 + 1.599*(D-0.35)
         ELSE
          EY=EYFCN(1.5,0.,0.2,B,BETA)+EYFCN(0.67,0.2,0.35,B,BETA)
     &        +EYFCN(1.599,0.35,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 0.85) .AND. (D .LT. 2.75)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 1.2
         ELSE
          EY=EYFCN(1.5,0.,0.2,B,BETA)+EYFCN(0.67,0.2,0.35,B,BETA)
     &        +EYFCN(1.599,0.35,0.85,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.2 + A0R*(D-2.75)
           ELSE
          EY=EYFCN(1.5,0.,0.2,B,BETA)+EYFCN(0.67,0.2,0.35,B,BETA)
     &        +EYFCN(1.599,0.35,0.85,B,BETA)+EYFCN(A0R,2.75,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 24.0) .AND. (L .LT. 30.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.6*(D-0.2)
         ELSE
          EY=EYFCN(1.6,0.2,D,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF    
          IF ((D .GE. 0.45) .AND. (D .LT. 1.6)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.4
         ELSE
          EY=EYFCN(1.6,0.2,0.45,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.6) .AND. (D .LT. 2.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4 + 1.25*(D-1.6)
         ELSE
          EY=EYFCN(1.6,0.2,0.45,B,BETA)+EYFCN(1.25,1.6,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF              
          IF ((D .GE. 2.0) .AND. (D .LT. 4.0)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.9
         ELSE
          EY=EYFCN(1.6,0.2,0.45,B,BETA)+EYFCN(1.25,1.6,2.,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.0)) THEN        
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9 + A0R*(D-4.0)
         ELSE
          EY=EYFCN(1.6,0.2,0.45,B,BETA)+EYFCN(1.25,1.6,2.,B,BETA0)
     &        +EYFCN(A0R,4.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 30.0) .AND. (L .LT. 40.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.15) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 2.*(D-0.15)
         ELSE
          EY=EYFCN(2.,0.15,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF             
          IF ((D .GE. 0.6) .AND. (D .LT. 2.8)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.9
         ELSE
          EY=EYFCN(2.,0.15,0.6,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9 + A0R*(D-2.8)
         ELSE
          EY=EYFCN(2.,0.15,0.6,B,BETA)+EYFCN(A0R,2.8,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 40.0) .AND. (L .LT. 50.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + D
         ELSE
          EY=EYFCN(1.,0.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF              
          IF ((D .GE. 0.1) .AND. (D .LT. 0.25)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.1
         ELSE
          EY=EYFCN(1.,0.,0.1,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 1.4*(D-0.25)
         ELSE
          EY=EYFCN(1.,0.,0.1,B,BETA)+EYFCN(1.4,0.25,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.45 + 0.11*(D-0.5)
         ELSE
          EY=EYFCN(1.,0.,0.1,B,BETA)+EYFCN(1.4,0.25,0.5,B,BETA)
     &        +EYFCN(0.11,0.5,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.95) .AND. (D .LT. 1.3)) THEN  
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4995 + 0.86*(D-0.95)
         ELSE
          EY=EYFCN(1.,0.,0.1,B,BETA)+EYFCN(1.4,0.25,0.5,B,BETA)
     &        +EYFCN(0.11,0.5,0.95,B,BETA)+EYFCN(0.86,0.95,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF                                               
        IF ((D .GE. 1.3) .AND. (D .LT. 3.0)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.8005
         ELSE
          EY=EYFCN(1.,0.,0.1,B,BETA)+EYFCN(1.4,0.25,0.5,B,BETA)
     &        +EYFCN(0.11,0.5,0.95,B,BETA)+EYFCN(0.86,0.95,1.3,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.8005 + A0R*(D-3.0)
             ELSE
          EY=EYFCN(1.,0.,0.1,B,BETA)+EYFCN(1.4,0.25,0.5,B,BETA)
     &        +EYFCN(0.11,0.5,0.95,B,BETA)+EYFCN(0.86,0.95,1.3,B,BETA0)
     &        +EYFCN(A0R,3.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 50.0) .AND. (L .LT. 60.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 3.*(D-0.25)
         ELSE
          EY=EYFCN(3.,0.25,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.4) .AND. (D .LT. 0.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.45 + 0.11*(D-0.4)
         ELSE
          EY=EYFCN(3.,0.25,0.4,B,BETA)+EYFCN(0.11,0.4,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
        IF ((D .GE. 0.85) .AND. (D .LT. 1.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.4995 + 1.202*(D-0.85)
         ELSE
          EY=EYFCN(3.,0.25,0.4,B,BETA)+EYFCN(0.11,0.4,0.85,B,BETA)
     &         +EYFCN(1.202,0.85,D,B,BETA0)
         ENDIF  
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.1) .AND. (D .LT. 2.9)) THEN
           IF ((IFLG .EQ. 0)) THEN
          EY = 0.8
         ELSE
          EY=EYFCN(3.,0.25,0.4,B,BETA)+EYFCN(0.11,0.4,0.85,B,BETA)
     &         +EYFCN(1.202,0.85,1.1,B,BETA0)
         ENDIF  
        GOTO 1000
        ENDIF
          IF ((D .GE. 2.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.8 + A0R*(D-2.9)
         ELSE
          EY=EYFCN(3.,0.25,0.4,B,BETA)+EYFCN(0.11,0.4,0.85,B,BETA)
     &         +EYFCN(1.202,0.85,1.1,B,BETA0)+EYFCN(A0R,2.9,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
 661   IF ((L .GE. 60.0) .AND. (L .LT. 62.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.3)) THEN
             IF ((IFLG .EQ. 0)) THEN
                EY = 0.33*D
             ELSE
                EY=EYFCN(0.33,0.,D,B,BETA)
             ENDIF 
             GOTO 1000
          ENDIF
          IF ((D .GE. 0.3) .AND. (D .LT. 0.7)) THEN
             IF ((IFLG .EQ. 0)) THEN
                EY = 0.099 + 2.2525*(D-0.3)
             ELSE
                EY=EYFCN(0.33,0.,0.3,B,BETA)+EYFCN(2.2525,0.3,D,B,BETA)
             ENDIF   
             GOTO 1000
          ENDIF              
          IF ((D .GE. 0.7) .AND. (D .LT. 2.1)) THEN
              IF ((IFLG .EQ. 0)) THEN
                 EY = 1.0
              ELSE
               EY=EYFCN(0.33,0.,0.3,B,BETA)+EYFCN(2.2525,0.3,0.7,B,BETA)
              ENDIF   
              GOTO 1000
          ENDIF
          IF ((D .GE. 2.1)) THEN    
             IF ((IFLG .EQ. 0)) THEN
                EY = 1.0 + A0R*(D-2.1)
             ELSE
               EY=EYFCN(0.33,0.,0.3,B,BETA)+EYFCN(2.2525,0.3,0.7,B,BETA)
     &            +EYFCN(A0R,2.1,D,B,BETA0)
             ENDIF
             GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 62.0) .AND. (L .LT. 63.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.3)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.3) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.05*(D-0.3)
         ELSE
          EY=EYFCN(1.05,0.3,D,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 1.25) .AND. (D .LT. 3.0)) THEN
           IF ((IFLG .EQ. 0)) THEN
            EY = 0.9975
         ELSE
          EY=EYFCN(1.05,0.3,1.25,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9975 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.05,0.3,1.25,B,BETA)+EYFCN(A0R,3.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 63.0) .AND. (L .LT. 70.0)) THEN
        IF ((D .GE. 0.0) .AND. (D .LT. 0.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.38*D
         ELSE
          EY=EYFCN(0.38,0.,D,B,BETA)
         ENDIF
        GOTO 1000
        ENDIF
        IF ((D .GE. 0.4) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.152 + 1.49*(D-0.4)
         ELSE
          EY=EYFCN(0.38,0.,0.4,B,BETA)+EYFCN(1.49,0.4,D,B,BETA)
         ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 0.6) .AND. (D .LT. 1.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.45 + 0.22*(D-0.6)
         ELSE
          EY=EYFCN(0.38,0.,0.4,B,BETA)+EYFCN(1.49,0.4,0.6,B,BETA)
     &        +EYFCN(0.22,0.6,D,B,BETA)
         ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.5) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.648 + 0.102*(D-1.5)
         ELSE
          EY=EYFCN(0.38,0.,0.4,B,BETA)+EYFCN(1.49,0.4,0.6,B,BETA)
     &        +EYFCN(0.22,0.6,1.5,B,BETA)+EYFCN(0.102,1.5,D,B,BETA0)
         ENDIF 
        GOTO 1000
        ENDIF
        IF ((D .GE. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.75 + A0R*(D-2.5)
         ELSE
          EY=EYFCN(0.38,0.,0.4,B,BETA)+EYFCN(1.49,0.4,0.6,B,BETA)
     &        +EYFCN(0.22,0.6,1.5,B,BETA)+EYFCN(0.102,1.5,2.5,B,BETA0)
     &        +EYFCN(A0R,2.5,D,B,BETA0)
         ENDIF
        GOTO 1000
        ENDIF
       ENDIF
       IF ((L .GE. 70.0) .AND. (L .LT. 72.0) .AND. (B .GE. -0.3)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.17*D
         ELSE
          EY=EYFCN(0.17,0.,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 0.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.102 + 0.66*(D-0.6)
         ELSE
          EY=EYFCN(0.17,0.,0.6,B,BETA)+EYFCN(0.66,0.6,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.9) .AND. (D .LT. 2.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3 + 1.5*(D-0.9)
         ELSE
          EY=EYFCN(0.17,0.,0.6,B,BETA)+EYFCN(0.66,0.6,0.9,B,BETA)
     &        +EYFCN(1.5,0.9,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.95 + A0R*(D-2.0)
         ELSE
          EY=EYFCN(0.17,0.,0.6,B,BETA)+EYFCN(0.66,0.6,0.9,B,BETA)
     &        +EYFCN(1.5,0.9,2.,B,BETA)+EYFCN(A0R,2.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 70.0) .AND. (L .LT. 74.0) .AND. (B .LT. -0.3)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.6)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 1.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 3.2*(D-0.6)
         ELSE
          EY=EYFCN(3.2,0.6,D,B,BETA)
           ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.6 + A0R*(D-1.1)
         ELSE
          EY=EYFCN(3.2,0.6,1.1,B,BETA)+EYFCN(A0R,1.1,D,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 72.0) .AND. (L .LT. 74.0) .AND. (B .GE. -0.3)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.25*D
         ELSE
          EY=EYFCN(0.25,0.,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 0.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15 + 0.83*(D-0.6)
         ELSE
          EY=EYFCN(0.25,0.,0.6,B,BETA)+EYFCN(0.83,0.6,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF              
        IF ((D .GE. 0.9) .AND. (D .LT. 1.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.399
         ELSE
          EY=EYFCN(0.25,0.,0.6,B,BETA)+EYFCN(0.83,0.6,0.9,B,BETA)
         ENDIF      
        GOTO 1000
        ENDIF
        IF ((D .GE. 1.9) .AND. (D .LT. 2.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.399 + 2.005*(D-1.9)
         ELSE
          EY=EYFCN(0.25,0.,0.6,B,BETA)+EYFCN(0.83,0.6,0.9,B,BETA)
     &        +EYFCN(2.005,1.9,D,B,BETA0)
         ENDIF               
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.8 + A0R*(D-2.1)
         ELSE
          EY=EYFCN(0.25,0.,0.6,B,BETA)+EYFCN(0.83,0.6,0.9,B,BETA)
     &        +EYFCN(2.005,1.9,2.1,B,BETA0)+EYFCN(A0R,2.1,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 74.0) .AND. (L .LT. 75.5)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.33*D
         ELSE
          EY=EYFCN(0.33,0.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF   
          IF ((D .GE. 0.45) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1485
         ELSE
          EY=EYFCN(0.33,0.,0.45,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.0) .AND. (D .LT. 1.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1485 + 1.548*(D-1.0)
         ELSE
          EY=EYFCN(0.33,0.,0.45,B,BETA)+EYFCN(1.548,1.,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF                                             
          IF ((D .GE. 1.55) .AND. (D .LT. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.9999
         ELSE
          EY=EYFCN(0.33,0.,0.45,B,BETA)+EYFCN(1.548,1.,1.55,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9999 + A0R*(D-4.0)
         ELSE
          EY=EYFCN(0.33,0.,0.45,B,BETA)+EYFCN(1.548,1.,1.55,B,BETA0)
     &        +EYFCN(A0R,4.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 75.5) .AND. (L .LT. 77.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4*D
         ELSE
          EY=EYFCN(0.4,0.,D,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 0.25) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 1.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 1.08*(D-0.6)
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)+EYFCN(1.08,0.6,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.396 + A0R*(D-1.8)
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)+EYFCN(1.08,0.6,1.8,B,BETA)
     &        +EYFCN(A0R,1.8,D,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF          
       ENDIF
       IF ((L .GE. 77.0) .AND. (L .LT. 80.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 1.15)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.67*(D-0.25)
         ELSE
          EY=EYFCN(1.67,0.25,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.15)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.503 + A0R*(D-1.15)
         ELSE
          EY=EYFCN(1.67,0.25,1.15,B,BETA)+EYFCN(A0R,1.15,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 80.0) .AND. (L .LT. 84.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.15) .AND. (D .LT. 0.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 3.57*(D-0.15)
         ELSE
          EY=EYFCN(3.57,0.15,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 2.499 + A0R*(D-0.85)
         ELSE
          EY=EYFCN(3.57,0.15,0.85,B,BETA)+EYFCN(A0R,0.85,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 84.0) .AND. (L .LT. 87.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.7)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 1.4*(D-0.2)
         ELSE
          EY=EYFCN(1.4,0.2,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF 
          IF ((D .GE. 0.7) .AND. (D .LT. 2.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.7
         ELSE
          EY=EYFCN(1.4,0.2,0.7,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.7 + A0R*(D-2.65)
         ELSE
          EY=EYFCN(1.4,0.2,0.7,B,BETA)+EYFCN(A0R,2.65,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 87.0) .AND. (L .LT. 90.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.25)
         ELSE
          EY=EYFCN(1.,0.25,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 0.4) .AND. (D .LT. 0.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15
         ELSE
          EY=EYFCN(1.,0.25,0.4,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.8) .AND. (D .LT. 2.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15 + 0.42*(D-0.8)
         ELSE
          EY=EYFCN(1.,0.25,0.4,B,BETA)+EYFCN(0.42,0.8,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                                           
          IF ((D .GE. 2.0) .AND. (D .LT. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.654
         ELSE
          EY=EYFCN(1.,0.25,0.4,B,BETA)+EYFCN(0.42,0.8,2.,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.654 + A0R*(D-3.5)
         ELSE
          EY=EYFCN(1.,0.25,0.4,B,BETA)+EYFCN(0.42,0.8,2.,B,BETA0)
     &        +EYFCN(A0R,3.5,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 90.0) .AND. (L .LT. 97.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.45*(D-0.25)
         ELSE
          EY=EYFCN(1.45,0.25,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 0.8) .AND. (D .LT. 2.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.7975
         ELSE
          EY=EYFCN(1.45,0.25,0.8,B,BETA)
         ENDIF          
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.6)) THEN
            IF ((IFLG .EQ. 0)) THEN
          EY = 0.7975 + A0R*(D-2.6)
         ELSE
          EY=EYFCN(1.45,0.25,0.8,B,BETA)+EYFCN(A0R,2.6,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 97.0) .AND. (L .LT. 100.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.25)
         ELSE
          EY=EYFCN(1.,0.25,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                  
          IF ((D .GE. 0.95) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.7
         ELSE
          EY=EYFCN(1.,0.25,0.95,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.7 + A0R*(D-2.5)
         ELSE
          EY=EYFCN(1.,0.25,0.95,B,BETA)+EYFCN(A0R,2.5,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
 662   IF ((L .GE. 100.0) .AND. (L .LT. 103.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.33*(D-0.2)
         ELSE
          EY=EYFCN(1.33,0.2,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 0.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.399 + 1.02*(D-0.5)
         ELSE
          EY=EYFCN(1.33,0.2,0.5,B,BETA)+EYFCN(1.02,0.5,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.55) .AND. (D .LT. 1.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.45 + 0.6*(D-0.55)
         ELSE
          EY=EYFCN(1.33,0.2,0.5,B,BETA)+EYFCN(1.02,0.5,0.55,B,BETA)
     &        +EYFCN(0.6,0.55,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                  
        IF ((D .GE. 1.3) .AND. (D. LT. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9
         ELSE
          EY=EYFCN(1.33,0.2,0.5,B,BETA)+EYFCN(1.02,0.5,0.55,B,BETA)
     &        +EYFCN(0.6,0.55,1.3,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9 + A0R*(D-4.0)
         ELSE
          EY=EYFCN(1.33,0.2,0.5,B,BETA)+EYFCN(1.02,0.5,0.55,B,BETA)
     &        +EYFCN(0.6,0.55,1.3,B,BETA0)+EYFCN(A0R,4.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 103.0) .AND. (L .LT. 107.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.15) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 2.15*(D-0.15)
         ELSE
          EY=EYFCN(2.15,0.15,D,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF      
          IF ((D .GE. 0.5) .AND. (D .LT. 4.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.7525
         ELSE
          EY=EYFCN(2.15,0.15,0.5,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.7525 + A0R*(D-4.4)
         ELSE
          EY=EYFCN(2.15,0.15,0.5,B,BETA)+EYFCN(A0R,4.4,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 107.0) .AND. (L .LT. 110.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*D
         ELSE
          EY=EYFCN(0.5,0.,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25 + 1.33*(D-0.5)
         ELSE
          EY=EYFCN(0.5,0.,0.5,B,BETA)+EYFCN(1.33,0.5,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                                            
          IF ((D .GE. 0.95) .AND. (D .LT. 3.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.8485
         ELSE
          EY=EYFCN(0.5,0.,0.5,B,BETA)+EYFCN(1.33,0.5,0.95,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.8485 + A0R*(D-3.6)
         ELSE
          EY=EYFCN(0.5,0.,0.5,B,BETA)+EYFCN(1.33,0.5,0.95,B,BETA0)
     &        +EYFCN(A0R,3.6,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 110.0) .AND. (L .LT. 113.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*D
         ELSE
          EY=EYFCN(0.5,0.,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25 + 1.7*(D-0.5)
         ELSE
          EY=EYFCN(0.5,0.,0.5,B,BETA)+EYFCN(1.7,0.5,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF       
        IF ((D .GE. 1.0) .AND. (D .LT. 3.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.1 
         ELSE
          EY=EYFCN(0.5,0.,0.5,B,BETA)+EYFCN(1.7,0.5,1.,B,BETA0)
         ENDIF 
        GOTO 1000
        ENDIF
          IF ((D .GE. 3.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.1 + A0R*(D-3.3)
         ELSE
          EY=EYFCN(0.5,0.,0.5,B,BETA)+EYFCN(1.7,0.5,1.,B,BETA0)
     &        +EYFCN(A0R,3.3,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 113.0) .AND. (L .LT. 117.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.33*D
         ELSE
          EY=EYFCN(1.33,0.,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF              
          IF ((D .GE. 0.15) .AND. (D .LT. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1995
         ELSE
          EY=EYFCN(1.33,0.,0.15,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.65) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1995 + 1.287*(D-0.65)
         ELSE                
          EY=EYFCN(1.33,0.,0.15,B,BETA)+EYFCN(1.287,0.65,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF              
          IF ((D .GE. 1.0) .AND. (D .LT. 3.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.64995
         ELSE
          EY=EYFCN(1.33,0.,0.15,B,BETA)+EYFCN(1.287,0.65,1.,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.64995 + A0R*(D-3.9)
         ELSE
          EY=EYFCN(1.33,0.,0.15,B,BETA)+EYFCN(1.287,0.65,1.,B,BETA)
     &        +EYFCN(A0R,3.9,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 117.0) .AND. (L .LT. 120.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.63*D
         ELSE
          EY=EYFCN(0.63,0.,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                  
          IF ((D .GE. 0.95) .AND. (D .LT. 3.7)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5985
         ELSE
          EY=EYFCN(0.63,0.,0.95,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.7)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.5985 + A0R*(D-3.7)
             ELSE
          EY=EYFCN(0.63,0.,0.95,B,BETA)+EYFCN(A0R,3.7,D,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 120.0) .AND. (L .LT. 122.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.25*(D-0.1)
         ELSE
          EY=EYFCN(1.25,0.1,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF   
          IF ((D .GE. 0.3) .AND. (D .LT. 0.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25
         ELSE
          EY=EYFCN(1.25,0.1,0.3,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.75) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25 + 1.5*(D-0.75)
         ELSE
          EY=EYFCN(1.25,0.1,0.3,B,BETA)+EYFCN(1.5,0.75,D,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.25) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 1.0
         ELSE
          EY=EYFCN(1.25,0.1,0.3,B,BETA)+EYFCN(1.5,0.75,1.25,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.0 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.25,0.1,0.3,B,BETA)+EYFCN(1.5,0.75,1.25,B,BETA)
     &        +EYFCN(A0R,3.,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 122.0) .AND. (L .LT. 124.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.22*D
         ELSE
          EY=EYFCN(0.22,0.,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.45) .AND. (D .LT. 1.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.099 + 0.9424*(D-0.45)
         ELSE
          EY=EYFCN(0.22,0.,0.45,B,BETA)+EYFCN(0.9424,0.45,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF                                              
          IF ((D .GE. 1.3) .AND. (D .LT. 4.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.90004
         ELSE
          EY=EYFCN(0.22,0.,0.45,B,BETA)+EYFCN(0.9424,0.45,1.3,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.90004 + A0R*(D-4.1)
         ELSE
          EY=EYFCN(0.22,0.,0.45,B,BETA)+EYFCN(0.9424,0.45,1.3,B,BETA)
     &        +EYFCN(A0R,4.1,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 124.0) .AND. (L .LT. 127.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.36*D
         ELSE
          EY=EYFCN(0.36,0.,D,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.55) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.198 + 1.56*(D-0.55)
         ELSE
          EY=EYFCN(0.36,0.,0.55,B,BETA)+EYFCN(1.56,0.55,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF   
          IF ((D .GE. 1.0) .AND. (D .LT. 3.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.9
         ELSE
          EY=EYFCN(0.36,0.,0.55,B,BETA)+EYFCN(1.56,0.55,1.,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9 + A0R*(D-3.65)
         ELSE
          EY=EYFCN(0.36,0.,0.55,B,BETA)+EYFCN(1.56,0.55,1.,B,BETA)
     &        +EYFCN(A0R,3.65,D,B,BETA0)
         ENDIF    
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 127.0) .AND. (L .LT. 130.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.17*D
         ELSE
          EY=EYFCN(0.17,0.,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN      
            EY = 0.102 + 1.37*(D-0.6)
         ELSE
          EY=EYFCN(0.17,0.,0.6,B,BETA)+EYFCN(1.37,0.6,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF      
          IF ((D .GE. 1.0) .AND. (D .LT. 4.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.65
         ELSE
          EY=EYFCN(0.17,0.,0.6,B,BETA)+EYFCN(1.37,0.6,1.,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.65 + A0R*(D-4.45)
         ELSE
          EY=EYFCN(0.17,0.,0.6,B,BETA)+EYFCN(1.37,0.6,1.,B,BETA)
     &        +EYFCN(A0R,4.45,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
 663   IF ((L .GE. 130.0) .AND. (L .LT. 132.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.2*D
         ELSE
          EY=EYFCN(0.2,0.,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 1.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 1.07*(D-0.5)
         ELSE
          EY=EYFCN(0.2,0.,0.5,B,BETA)+EYFCN(1.07,0.5,D,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF              
          IF ((D .GE. 1.2) .AND. (D .LT. 4.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.849
         ELSE
          EY=EYFCN(0.2,0.,0.5,B,BETA)+EYFCN(1.07,0.5,1.2,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.849 + A0R*(D-4.25)
         ELSE
          EY=EYFCN(0.2,0.,0.5,B,BETA)+EYFCN(1.07,0.5,1.2,B,BETA)
     &        +EYFCN(A0R,4.25,D,B,BETA0)
         ENDIF
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 132.0) .AND. (L .LT. 134.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 2.5*(D-0.2)
         ELSE
          EY=EYFCN(2.5,0.2,D,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.4) .AND. (D .LT. 1.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5 + 0.81*(D-0.4)
         ELSE
          EY=EYFCN(2.5,0.2,0.4,B,BETA)+EYFCN(0.81,0.4,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF     
        IF ((D .GE. 1.2) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.148
         ELSE
          EY=EYFCN(2.5,0.2,0.4,B,BETA)+EYFCN(0.81,0.4,1.2,B,BETA)
         ENDIF   
        GOTO 1000
        ENDIF              
          IF ((D .GE. 2.5)) THEN      
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.148 + A0R*(D-2.5)
         ELSE
          EY=EYFCN(2.5,0.2,0.4,B,BETA)+EYFCN(0.81,0.4,1.2,B,BETA)
     &        +EYFCN(A0R,2.5,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 134.0) .AND. (L .LT. 137.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.5*(D-0.1)
         ELSE
          EY=EYFCN(1.5,0.1,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.3) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3 + 0.86*(D-0.3)
         ELSE
          EY=EYFCN(1.5,0.1,0.3,B,BETA)+EYFCN(0.86,0.3,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF           
          IF ((D .GE. 1.0) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.902
         ELSE
          EY=EYFCN(1.5,0.1,0.3,B,BETA)+EYFCN(0.86,0.3,1.,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.902 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.5,0.1,0.3,B,BETA)+EYFCN(0.86,0.3,1.,B,BETA)
     &        +EYFCN(A0R,3.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 137.0) .AND. (L .LT. 140.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.2)
         ELSE
          EY=EYFCN(1.,0.2,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.55) .AND. (D .LT. 0.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35 + 1.71*(D-0.55)
         ELSE
          EY=EYFCN(1.,0.2,0.55,B,BETA)+EYFCN(1.71,0.55,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF              
        IF ((D .GE. 0.9) .AND. (D .LT. 3.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9485
         ELSE
          EY=EYFCN(1.,0.2,0.55,B,BETA)+EYFCN(1.71,0.55,0.9,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9485 + A0R*(D-3.2)
         ELSE
          EY=EYFCN(1.,0.2,0.55,B,BETA)+EYFCN(1.71,0.55,0.9,B,BETA)
     &        +EYFCN(A0R,3.2,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF             
       IF ((L .GE. 140.0) .AND. (L .LT. 150.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*D 
         ELSE
          EY=EYFCN(0.5,0.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 0.7)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 3.5*(D-0.5)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(3.5,0.5,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF              
          IF ((D .GE. 0.7) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.8
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(3.5,0.5,0.7,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.5) .AND. (D .LT. 3.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.8 + 0.5*(D-2.5)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(3.5,0.5,0.7,B,BETA)
     &        +EYFCN(0.5,2.5,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.2 + A0R*(D-3.3)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(3.5,0.5,0.7,B,BETA)
     &        +EYFCN(0.5,2.5,3.3,B,BETA)+EYFCN(A0R,3.3,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 150.0) .AND. (L .LT. 160.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*D
         ELSE
          EY=EYFCN(0.5,0.,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 0.33*(D-0.2)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(0.33,0.2,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.35) .AND. (D .LT. 0.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1495 + 0.3763*(D-0.35)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(0.33,0.2,0.35,B,BETA)
     &        +EYFCN(0.3763,0.35,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.75) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3002 + 1.749*(D-0.75)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(0.33,0.2,0.35,B,BETA)
     &        +EYFCN(0.3763,0.35,0.75,B,BETA)
     &        +EYFCN(1.749,0.75,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.95) .AND. (D .LT. 1.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.65 + 0.4*(D-0.95)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(0.33,0.2,0.35,B,BETA)
     &        +EYFCN(0.3763,0.35,0.75,B,BETA)
     &        +EYFCN(1.749,0.75,0.95,B,BETA)
     &        +EYFCN(0.4,0.95,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 1.45) .AND. (D .LT. 2.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.85
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(0.33,0.2,0.35,B,BETA)
     &        +EYFCN(0.3763,0.35,0.75,B,BETA)
     &        +EYFCN(1.749,0.75,0.95,B,BETA)
     &        +EYFCN(0.4,0.95,1.45,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.85 + A0R*(D-2.85)
         ELSE
          EY=EYFCN(0.5,0.,0.2,B,BETA)+EYFCN(0.33,0.2,0.35,B,BETA)
     &        +EYFCN(0.3763,0.35,0.75,B,BETA)
     &        +EYFCN(1.749,0.75,0.95,B,BETA)
     &        +EYFCN(0.4,0.95,1.45,B,BETA)+EYFCN(A0R,2.85,D,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
 664   IF ((L .GE. 160.0) .AND. (L .LT. 170.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.6)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 0.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.33*(D-0.6)
         ELSE
          EY=EYFCN(1.33,0.6,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                  
          IF ((D .GE. 0.9) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.399
         ELSE
          EY=EYFCN(1.33,0.6,0.9,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.5) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.399 + 0.802*(D-2.5)
         ELSE
          EY=EYFCN(1.33,0.6,0.9,B,BETA)+EYFCN(0.802,2.5,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.8 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.33,0.6,0.9,B,BETA)+EYFCN(0.802,2.5,3.,B,BETA0)
     &        +EYFCN(A0R,3.,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 170.0) .AND. (L .LT. 174.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.1)
         ELSE
          EY=EYFCN(1.,0.1,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF                
          IF ((D .GE. 0.5) .AND. (D .LT. 1.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4
         ELSE
          EY=EYFCN(1.,0.1,0.5,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.5) .AND. (D .LT. 1.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4 + 0.43*(D-1.5)
         ELSE
          EY=EYFCN(1.,0.1,0.5,B,BETA)+EYFCN(0.43,1.5,D,B,BETA0)
         ENDIF        
        GOTO 1000
          ENDIF              
          IF ((D .GE. 1.85) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5505
         ELSE
          EY=EYFCN(1.,0.1,0.5,B,BETA)+EYFCN(0.43,1.5,1.85,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.5505 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.,0.1,0.5,B,BETA)+EYFCN(0.43,1.5,1.85,B,BETA0)
     &        +EYFCN(A0R,3.,D,B,BETA0)
         ENDIF        
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 174.0) .AND. (L .LT. 184.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.6*D
         ELSE
          EY=EYFCN(0.6,0.,D,B,BETA)
         ENDIF        
        GOTO 1000
          ENDIF              
          IF ((D .GE. 0.5) .AND. (D .LT. 1.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3 
         ELSE
          EY=EYFCN(0.6,0.,0.5,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.95) .AND. (D .LT. 2.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3 + (D-1.95)
         ELSE
          EY=EYFCN(0.6,0.,0.5,B,BETA)+EYFCN(1.,1.95,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.95 + A0R*(D-2.6)
         ELSE
          EY=EYFCN(0.6,0.,0.5,B,BETA)+EYFCN(1.,1.95,2.6,B,BETA0)
     &        +EYFCN(A0R,2.6,D,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 184.0) .AND. (L .LT. 190.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.05)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.05) .AND. (D .LT. 0.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*(D-0.05)
         ELSE
          EY=EYFCN(0.5,0.05,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 0.25) .AND. (D .LT. 0.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1
         ELSE
          EY=EYFCN(0.5,0.05,0.25,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.8) .AND. (D .LT. 1.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 1.7*(D-0.8)
         ELSE
          EY=EYFCN(0.5,0.05,0.25,B,BETA)+EYFCN(1.7,0.8,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF              
        IF ((D .GE. 1.3) .AND. (D .LT. 3.35)) THEN 
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.95
         ELSE
          EY=EYFCN(0.5,0.05,0.25,B,BETA)+EYFCN(1.7,0.8,1.3,B,BETA)
         ENDIF   
        GOTO 1000
        ENDIF
          IF ((D .GE. 3.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.95 + A0R*(D-3.35)
         ELSE
          EY=EYFCN(0.5,0.05,0.25,B,BETA)+EYFCN(1.7,0.8,1.3,B,BETA)
     &        +EYFCN(A0R,3.35,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
 665   IF ((L .GE. 190.0) .AND. (L .LT. 200.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*(D-0.1)
         ELSE
          EY=EYFCN(0.5,0.1,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.3) .AND. (D .LT. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1
         ELSE
          EY=EYFCN(0.5,0.1,0.3,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.65) .AND. (D .LT. 1.05)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 0.88*(D-0.65)
         ELSE
          EY=EYFCN(0.5,0.1,0.3,B,BETA)+EYFCN(0.88,0.65,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.05) .AND. (D .LT. 1.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.452
         ELSE
          EY=EYFCN(0.5,0.1,0.3,B,BETA)+EYFCN(0.88,0.65,1.05,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.8) .AND. (D .LT. 2.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.452 + 0.44*(D-1.8)
         ELSE
          EY=EYFCN(0.5,0.1,0.3,B,BETA)+EYFCN(0.88,0.65,1.05,B,BETA)
     &        +EYFCN(0.44,1.8,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF              
        IF ((D .GE. 2.25) .AND. (D .LT. 4.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.65
         ELSE
          EY=EYFCN(0.5,0.1,0.3,B,BETA)+EYFCN(0.88,0.65,1.05,B,BETA)
     &        +EYFCN(0.44,1.8,2.25,B,BETA0)
         ENDIF 
        GOTO 1000
        ENDIF
          IF ((D .GE. 4.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.65 + A0R*(D-4.4)
         ELSE
          EY=EYFCN(0.5,0.1,0.3,B,BETA)+EYFCN(0.88,0.65,1.05,B,BETA)
     &        +EYFCN(0.44,1.8,2.25,B,BETA0)+EYFCN(A0R,4.4,D,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 200.0) .AND. (L .LT. 205.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25*D
         ELSE
          EY=EYFCN(0.25,0.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF     
          IF ((D .GE. 0.2) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.05
         ELSE
          EY=EYFCN(0.25,0.,0.2,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.6) .AND. (D .LT. 0.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.05 + 0.4*(D-0.6)
         ELSE
          EY=EYFCN(0.25,0.,0.2,B,BETA)+EYFCN(0.4,0.6,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF              
          IF ((D .GE. 0.85) .AND. (D .LT. 1.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15
         ELSE
          EY=EYFCN(0.25,0.,0.2,B,BETA)+EYFCN(0.4,0.6,0.85,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.75) .AND. (D .LT. 2.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15 + (D-1.75)
         ELSE
          EY=EYFCN(0.25,0.,0.2,B,BETA)+EYFCN(0.4,0.6,0.85,B,BETA)
     &        +EYFCN(1.,1.75,D,B,BETA0)
         ENDIF    
        GOTO 1000
          ENDIF                
          IF ((D .GE. 2.0) .AND. (D .LT. 3.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4
         ELSE
          EY=EYFCN(0.25,0.,0.2,B,BETA)+EYFCN(0.4,0.6,0.85,B,BETA)
     &        +EYFCN(1.,1.75,2.,B,BETA0)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.4 + A0R*(D-3.65)
         ELSE
          EY=EYFCN(0.25,0.,0.2,B,BETA)+EYFCN(0.4,0.6,0.85,B,BETA)
     &        +EYFCN(1.,1.75,2.,B,BETA0)+EYFCN(A0R,3.65,D,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 205.0) .AND. (L .LT. 210.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.15*D
         ELSE
          EY=EYFCN(0.15,0.,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.65) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0975 + 1.1708*(D-0.65)
         ELSE
          EY=EYFCN(0.15,0.,0.65,B,BETA)+EYFCN(1.1708,0.65,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.79998 + A0R*(D-1.25)
         ELSE
          EY=EYFCN(0.15,0.,0.65,B,BETA)+EYFCN(1.1708,0.65,1.25,B,BETA)
     &        +EYFCN(A0R,1.25,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 210.0) .AND. (L .LT. 220.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.9)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.9) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.9)
         ELSE
          EY=EYFCN(1.,0.9,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 1.25) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35
         ELSE
          EY=EYFCN(1.,0.9,1.25,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0) .AND. (D .LT. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35 + 0.8*(D-3.)
         ELSE
          EY=EYFCN(1.,0.9,1.25,B,BETA)+EYFCN(0.8,3.,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.75 + A0R*(D-3.5)
         ELSE
          EY=EYFCN(1.,0.9,1.25,B,BETA)+EYFCN(0.8,3.,3.5,B,BETA0)
     &        +EYFCN(A0R,3.5,D,B,BETA0)
         ENDIF    
        GOTO 1000
          ENDIF
       ENDIF
 666   IF ((L .GE. 220.0) .AND. (L .LT. 230.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.3)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.3) .AND. (D .LT. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.3)
         ELSE
          EY=EYFCN(1.,0.3,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF              
          IF ((D .GE. 0.65) .AND. (D .LT. 1.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35
         ELSE
          EY=EYFCN(1.,0.3,0.65,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.75) .AND. (D .LT. 2.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35 + 0.43*(D-1.75)
          ELSE
          EY=EYFCN(1.,0.3,0.65,B,BETA)+EYFCN(0.43,1.75,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF         
          IF ((D .GE. 2.1) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5005      
         ELSE
          EY=EYFCN(1.,0.3,0.65,B,BETA)+EYFCN(0.43,1.75,2.1,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.5005 + A0R*(D-2.5)
         ELSE
          EY=EYFCN(1.,0.3,0.65,B,BETA)+EYFCN(0.43,1.75,2.1,B,BETA0)
     &        +EYFCN(A0R,2.5,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 230.0) .AND. (L .LT. 240.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 1.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1*D
         ELSE
          EY=EYFCN(0.1,0.,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.15 + A0R*(D-1.5)
         ELSE
          EY=EYFCN(0.1,0.,1.5,B,BETA)+EYFCN(A0R,1.5,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 240.0) .AND. (L .LT. 250.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.5*(D-0.1)
         ELSE
          EY=EYFCN(0.5,0.1,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.5) .AND. (D .LT. 3.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.2
         ELSE
          EY=EYFCN(0.5,0.1,0.5,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.4) .AND. (D .LT. 4.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.2 + 0.27*(D-3.4)
         ELSE
          EY=EYFCN(0.5,0.1,0.5,B,BETA)+EYFCN(0.27,3.4,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF  
          IF ((D .GE. 4.5) .AND. (D .LT. 9.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.497
         ELSE
          EY=EYFCN(0.5,0.1,0.5,B,BETA)+EYFCN(0.27,3.4,4.5,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 9.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.497 + A0R*(D-9.0)
         ELSE
          EY=EYFCN(0.5,0.1,0.5,B,BETA)+EYFCN(0.27,3.4,4.5,B,BETA0)
     &        +EYFCN(A0R,9.,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
 667   IF ((L .GE. 250.0) .AND. (L .LT. 257.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.4)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.4) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 1.25*(D-0.4)
         ELSE
          EY=EYFCN(1.25,0.4,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                  
          IF ((D .GE. 0.6) .AND. (D .LT. 2.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25
         ELSE
          EY=EYFCN(1.25,0.4,0.6,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.0) .AND. (D .LT. 2.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.25 + 0.47*(D-2.)
         ELSE
          EY=EYFCN(1.25,0.4,0.6,B,BETA)+EYFCN(0.47,2.,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF                                             
          IF ((D .GE. 2.75) .AND. (D .LT. 4.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.6025      
         ELSE
          EY=EYFCN(1.25,0.4,0.6,B,BETA)+EYFCN(0.47,2.,2.75,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.6025 + A0R*(D-4.25)
         ELSE
          EY=EYFCN(1.25,0.4,0.6,B,BETA)+EYFCN(0.47,2.,2.75,B,BETA0)
     &        +EYFCN(A0R,4.25,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 257.0) .AND. (L .LT. 264.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.65)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.65) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.86*(D-0.65)
         ELSE
          EY=EYFCN(1.86,0.65,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF                
          IF ((D .GE. 1.0) .AND. (D .LT. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.651
         ELSE
          EY=EYFCN(1.86,0.65,1.,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.651 + A0R*(D-3.0)
         ELSE
          EY=EYFCN(1.86,0.65,1.,B,BETA)+EYFCN(A0R,3.,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 264.0) .AND. (L .LT. 267.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.9)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.9) .AND. (D .LT. 1.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 2.*(D-0.9)
         ELSE
          EY=EYFCN(2.,0.9,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 1.35) .AND. (D .LT. 3.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.9
         ELSE
          EY=EYFCN(2.,0.9,1.35,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.9 + A0R*(D-3.1)
         ELSE
          EY=EYFCN(2.,0.9,1.35,B,BETA)+EYFCN(A0R,3.1,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 267.0) .AND. (L .LT. 270.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.105*D
         ELSE
          EY=EYFCN(0.105,0.,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.95) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.09975 + 2.67*(D-0.95)
         ELSE
          EY=EYFCN(0.105,0.,0.95,B,BETA)+EYFCN(2.67,0.95,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF          
          IF ((D .GE. 1.25) .AND. (D .LT. 3.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.90075
         ELSE
          EY=EYFCN(0.105,0.,0.95,B,BETA)+EYFCN(2.67,0.95,1.25,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.2)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.90075 + A0R*(D-3.2)
         ELSE
          EY=EYFCN(0.105,0.,0.95,B,BETA)+EYFCN(2.67,0.95,1.25,B,BETA)
     &        +EYFCN(A0R,3.2,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 270.0) .AND. (L .LT. 280.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 1.0)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 1.0) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 2.*(D-1.)
         ELSE
          EY=EYFCN(2.,1.,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                
          IF ((D .GE. 1.25) .AND. (D .LT. 3.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5
         ELSE
          EY=EYFCN(2.,1.,1.25,B,BETA)
         ENDIF
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.5 + A0R*(D-3.9)
         ELSE
          EY=EYFCN(2.,1.,1.25,B,BETA)+EYFCN(A0R,3.9,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
 668   IF ((L .GE. 280.0) .AND. (L .LT. 290.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.45)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.45) .AND. (D .LT. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5*(D-0.45)
         ELSE
          EY=EYFCN(0.5,0.45,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 0.65) .AND. (D .LT. 1.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 
         ELSE
          EY=EYFCN(0.5,0.45,0.65,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.3) .AND. (D .LT. 1.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + 1.5*(D-1.3)
         ELSE
          EY=EYFCN(0.5,0.45,0.65,B,BETA)+EYFCN(1.5,1.3,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.5) .AND. (D .LT. 2.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4
         ELSE
          EY=EYFCN(0.5,0.45,0.65,B,BETA)+EYFCN(1.5,1.3,1.5,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.4 + A0R*(D-2.65)
         ELSE
          EY=EYFCN(0.5,0.45,0.65,B,BETA)+EYFCN(1.5,1.3,1.5,B,BETA0)
     &        +EYFCN(A0R,2.65,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 290.0) .AND. (L .LT. 291.5)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.55)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.55) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.11*(D-0.55)
           ELSE
          EY=EYFCN(1.11,0.55,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.0) .AND. (D .LT. 2.9)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4995 
         ELSE
          EY=EYFCN(1.11,0.55,1.,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.9) .AND. (D .LT. 3.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4995 + 1.1546*(D-2.9)
         ELSE
          EY=EYFCN(1.11,0.55,1.,B,BETA)+EYFCN(1.1546,2.9,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.24999 + A0R*(D-3.55)
         ELSE
          EY=EYFCN(1.11,0.55,1.,B,BETA)+EYFCN(1.1546,2.9,3.55,B,BETA0)
     &        +EYFCN(A0R,3.55,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 291.5) .AND. (L .LT. 293.5)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4*D
         ELSE
          EY=EYFCN(0.4,0.,D,B,BETA)
         ENDIF         
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 0.25) .AND. (D .LT. 0.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.8) .AND. (D .LT. 1.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1 + (D-0.8)
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)+EYFCN(1.,0.8,D,B,BETA)
         ENDIF         
        GOTO 1000
          ENDIF                                          
          IF ((D .GE. 1.1) .AND. (D .LT. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.4      
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)+EYFCN(1.,0.8,1.1,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.4 + A0R*(D-4.0)
         ELSE
          EY=EYFCN(0.4,0.,0.25,B,BETA)+EYFCN(1.,0.8,1.1,B,BETA)
     &        +EYFCN(A0R,4.0,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF               
       ENDIF
       IF ((L .GE. 293.5) .AND. (L .LT. 300.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.25)
         ELSE
          EY=EYFCN(1.,0.25,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 0.6) .AND. (D .LT. 4.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35 
         ELSE
          EY=EYFCN(1.,0.25,0.6,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.35 + A0R*(D-4.5)
         ELSE
          EY=EYFCN(1.,0.25,0.6,B,BETA)+EYFCN(A0R,4.5,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 300.0) .AND. (L .LT. 310.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.15) .AND. (D .LT. 0.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 1.57*(D-0.15)
         ELSE
          EY=EYFCN(1.57,0.15,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF                   
          IF ((D .GE. 0.5) .AND. (D .LT. 1.6)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5495
         ELSE
          EY=EYFCN(1.57,0.15,0.5,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.6) .AND. (D .LT. 1.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5495 + 0.8586*(D-1.6)
         ELSE
          EY=EYFCN(1.57,0.15,0.5,B,BETA)+EYFCN(0.8586,1.6,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.85001 + A0R*(D-1.95)
         ELSE
          EY=EYFCN(1.57,0.15,0.5,B,BETA)+EYFCN(0.8586,1.6,1.95,B,BETA0)
     &        +EYFCN(A0R,1.95,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
 669   IF ((L .GE. 310.0) .AND. (L .LT. 320.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.15) .AND. (D .LT. 0.4)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.8*(D-0.15)
         ELSE
          EY=EYFCN(0.8,0.15,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF                  
          IF ((D .GE. 0.4) .AND. (D .LT. 0.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.2
         ELSE
          EY=EYFCN(0.8,0.15,0.4,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.75) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.2 + 0.8*(D-0.75)
         ELSE
          EY=EYFCN(0.8,0.15,0.4,B,BETA)+EYFCN(0.8,0.75,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.4 + A0R*(D-1.0)
         ELSE
          EY=EYFCN(0.8,0.15,0.4,B,BETA)+EYFCN(0.8,0.75,1.,B,BETA0)
     &        +EYFCN(A0R,1.,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 320.0) .AND. (L .LT. 330.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.2)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.2) .AND. (D .LT. 0.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.55*(D-0.2)
         ELSE
          EY=EYFCN(0.55,0.2,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF                  
          IF ((D .GE. 0.75) .AND. (D .LT. 3.1)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3025
         ELSE
          EY=EYFCN(0.55,0.2,0.75,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.1) .AND. (D .LT. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.3025 + 1.13*(D-3.1)
         ELSE
          EY=EYFCN(0.55,0.2,0.75,B,BETA)+EYFCN(1.13,3.1,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.7545 + A0R*(D-3.5)
         ELSE
          EY=EYFCN(0.55,0.2,0.75,B,BETA)+EYFCN(1.13,3.1,3.5,B,BETA0)
     &        +EYFCN(A0R,3.5,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 330.0) .AND. (L .LT. 333.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.7)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.7) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 3.33*(D-0.7)
         ELSE
          EY=EYFCN(3.33,0.7,D,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF  
          IF ((D .GE. 1.0) .AND. (D .LT. 4.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.999
         ELSE
          EY=EYFCN(3.33,0.7,1.0,B,BETA)
         ENDIF   
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.999 + A0R*(D-4.3)
              ELSE
          EY=EYFCN(3.33,0.7,1.0,B,BETA)+EYFCN(A0R,4.3,D,B,BETA0)
         ENDIF   
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 333.0) .AND. (L .LT. 337.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.15)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.15) .AND. (D .LT. 0.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.5*(D-0.15)
         ELSE
          EY=EYFCN(1.5,0.15,D,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF  
          IF ((D .GE. 0.25) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15
         ELSE
          EY=EYFCN(1.5,0.15,0.25,B,BETA)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.95) .AND. (D .LT. 1.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15 + 1.86*(D-0.95)
         ELSE
          EY=EYFCN(1.5,0.15,0.25,B,BETA)+EYFCN(1.86,0.95,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF          
          IF ((D .GE. 1.3) .AND. (D .LT. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.801 
         ELSE
          EY=EYFCN(1.5,0.15,0.25,B,BETA)+EYFCN(1.86,0.95,1.3,B,BETA0)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.801 + A0R*(D-3.5)
         ELSE
          EY=EYFCN(1.5,0.15,0.25,B,BETA)+EYFCN(1.86,0.95,1.3,B,BETA0)
     &        +EYFCN(A0R,3.5,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 337.0) .AND. (L .LT. 338.5)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.8)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.8) .AND. (D .LT. 1.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 2.75*(D-0.8)
         ELSE
          EY=EYFCN(2.75,0.8,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF        
          IF ((D .GE. 1.0) .AND. (D .LT. 2.05)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.55
         ELSE
          EY=EYFCN(2.75,0.8,1.,B,BETA)
         ENDIF    
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.05) .AND. (D .LT. 2.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.55 + 0.33*(D-2.05)
         ELSE
          EY=EYFCN(2.75,0.8,1.,B,BETA)+EYFCN(0.33,2.05,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF  
          IF ((D .GE. 2.5) .AND. (D .LT. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.6985
         ELSE
          EY=EYFCN(2.75,0.8,1.,B,BETA)+EYFCN(0.33,2.05,2.5,B,BETA0)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 4.0)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.6985 + A0R*(D-4.0)
         ELSE
          EY=EYFCN(2.75,0.8,1.,B,BETA)+EYFCN(0.33,2.05,2.5,B,BETA0)
     &        +EYFCN(A0R,4.0,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF     
       IF ((L .GE. 338.5) .AND. (L .LT. 340.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.1)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.1) .AND. (D .LT. 0.3)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 2.5*(D-0.1)
         ELSE
          EY=EYFCN(2.5,0.1,D,B,BETA)
         ENDIF        
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 0.3) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5
         ELSE
          EY=EYFCN(2.5,0.1,0.3,B,BETA)
         ENDIF  
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.25) .AND. (D .LT. 1.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.5 + 0.7*(D-1.25)
         ELSE
          EY=EYFCN(2.5,0.1,0.3,B,BETA)+EYFCN(0.7,1.25,D,B,BETA0)
         ENDIF        
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.75) .AND. (D .LT. 2.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.85
         ELSE
          EY=EYFCN(2.5,0.1,0.3,B,BETA)+EYFCN(0.7,1.25,1.75,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 2.85)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.85 + A0R*(D-2.85)
         ELSE
          EY=EYFCN(2.5,0.1,0.3,B,BETA)+EYFCN(0.7,1.25,1.75,B,BETA0)
     &        +EYFCN(A0R,2.85,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
       ENDIF
 670   IF ((L .GE. 340.0) .AND. (L .LT. 343.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.35)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.29*D
         ELSE
          EY=EYFCN(0.29,0.,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.35) .AND. (D .LT. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.1015 + 3.28*(D-0.35)
         ELSE
          EY=EYFCN(0.29,0.,0.35,B,BETA)+EYFCN(3.28,0.35,D,B,BETA)
         ENDIF      
        GOTO 1000
          ENDIF
          IF ((D .GE. 0.65)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.0855 + A0R*(D-0.65)
         ELSE
          EY=EYFCN(0.29,0.,0.35,B,BETA)+EYFCN(3.28,0.35,0.65,B,BETA)
     &        +EYFCN(A0R,0.65,D,B,BETA0)
         ENDIF      
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 343.0) .AND. (L .LT. 350.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.3)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.3) .AND. (D .LT. 0.95)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + (D-0.3)
         ELSE
          EY=EYFCN(1.,0.3,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF                 
          IF ((D .GE. 0.95) .AND. (D .LT. 1.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.65
         ELSE
          EY=EYFCN(1.,0.3,0.95,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.8)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 0.65 + A0R*(D-1.8)
         ELSE
          EY=EYFCN(1.,0.3,0.95,B,BETA)+EYFCN(A0R,1.8,D,B,BETA0)
         ENDIF       
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 350.0) .AND. (L .LT. 357.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 0.25)) THEN
            EY = 0.0
          GOTO 1000
          ENDIF
          IF ((D .GE. 0.25) .AND. (D .LT. 0.45)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 1.75*(D-0.25)
         ELSE
          EY=EYFCN(1.75,0.25,D,B,BETA)
         ENDIF       
        GOTO 1000
          ENDIF                    
          IF ((D .GE. 0.45) .AND. (D .LT. 1.15)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35
         ELSE
          EY=EYFCN(1.75,0.25,0.45,B,BETA)
         ENDIF     
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.15) .AND. (D .LT. 1.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.35 + 2.63*(D-1.15)
         ELSE
          EY=EYFCN(1.75,0.25,0.45,B,BETA)+EYFCN(2.63,1.15,D,B,BETA0)
         ENDIF        
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.55)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.402 + A0R*(D-1.55)
         ELSE
          EY=EYFCN(1.75,0.25,0.45,B,BETA)+EYFCN(2.63,1.15,1.55,B,BETA0)
     &        +EYFCN(A0R,1.55,D,B,BETA0)
         ENDIF        
        GOTO 1000
          ENDIF
       ENDIF
       IF ((L .GE. 357.0) .AND. (L .LT. 360.0)) THEN
          IF ((D .GE. 0.0) .AND. (D .LT. 1.25)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.0 + 0.12*D
         ELSE
          EY=EYFCN(0.12,0.,D,B,BETA)
         ENDIF        
        GOTO 1000
          ENDIF
          IF ((D .GE. 1.25) .AND. (D .LT. 1.75)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 0.15 + 2.*(D-1.25)
         ELSE
          EY=EYFCN(0.12,0.,1.25,B,BETA)+EYFCN(2.,1.25,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF      
          IF ((D .GE. 1.75) .AND. (D .LT. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
            EY = 1.15
         ELSE
          EY=EYFCN(0.12,0.,1.25,B,BETA)+EYFCN(2.,1.25,1.75,B,BETA0)
         ENDIF 
        GOTO 1000
          ENDIF
          IF ((D .GE. 3.5)) THEN
         IF ((IFLG .EQ. 0)) THEN
          EY = 1.15 + A0R*(D-3.5)
         ELSE
          EY=EYFCN(0.12,0.,1.25,B,BETA)+EYFCN(2.,1.25,1.75,B,BETA0)
     &        +EYFCN(A0R,3.5,D,B,BETA0)
         ENDIF     
        GOTO 1000
          ENDIF
       ENDIF
C---------------------------------------------------------------------------
C      Error calculation (correlation coefficient r = 0.998).
C---------------------------------------------------------------------------
c 1000  SEY= 0.059*EY*EY+0.218*EY+0.44
 1000  AV = EY*R
       SAV= SQRT(0.132*0.132+AV*AV/9.)
C       SAV= SEY*R
       IF (ABS(B) .EQ. 90.) SAV=0.0
       RETURN
       END

      FUNCTION EYFCN(S,D0,D,B,BETA)
      IMPLICIT REAL*4 (A-H,L-Z)
C---------------------------------------------------------------------------
C     This function calculates the color excess Ey for stars significantly
C     out of the Galactic plane (z > 55 pc), and is used in FITZGERALD. 
C     VARIABLE DEFINITIONS: S is the differential color excess per kpc,
C     D0 is the distance to the front edge of the cloud, D is the distance 
C     to the star (or distance to the back edge of the intervening cloud), 
C     b is galactic latitude, and beta is the half width of the material 
C     perpendicular to the galactic plane.
C---------------------------------------------------------------------------
      IF (ABS(B) .GT. 89.999) THEN
         EYFCN=0.
         RETURN
      ENDIF
      IF (D0 .GT. 0.) THEN
         E0=EXP(-D0*ABS(TAN(B*3.141592/180.))/BETA)
      ELSE
         E0=1.
      ENDIF
      E = EXP(-D * ABS(TAN(B * 3.141592 / 180.)) / BETA)
      EYFCN=BETA*S*(E0-E)/ABS(SIN(B*3.141592/180.))
      RETURN
      END

      SUBROUTINE NECKEL_KLARE(LL,BB,D,AV,SAV,A0)
C-----------------------------------------------------------------------
C     From Neckel and Klare, A&A Supp 42, 251 (1980).
C-----------------------------------------------------------------------     
C MAIN SUBROUTINE FOR THE NECKEL AND KLARE PAPER

C THERE ARE TWO FUNCTIONS AND 21 ADDITIONAL SUBROUTINES ASSOCIATED
C WITH THIS COMPUTER MODEL OF THE NECKEL AND KLARE PAPER.
C THE SUBROUTINES WITH THE NAMES MAP####(12 IN ALL) ARE THE SKYMAP
C SELECTION SUBROUTINES. ONE FOR EACH OF THE 12 MAPS(PAGES 258 THRU
C 261) GIVEN IN THE PAPER. THE SUBROUTINES WITH THE NAMES ENK#(8
C TOTAL)ARE THE EXTINCTION SUBROUTINES. THE VARIABLE CELL IN THE ENK#
C SUBROUTINES IS ASSOCIATED WITH THE EXTINCTIONS PLOTS GIVEN IN THE
C NK PAPER(PAGES 262 THRU 276). FINALLY THE SUBROUTINE ERR IS USED
C FOR FINDING THE SAV OF THE EXTINCTION FOR A GIVEN DISTANCE.
C-----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,L-Z)

C SETTING SENTINAL VALUES
      AV = -99.0
      SAV = -99.0

C FORCING L AND B COORDINATES INTO 0.05 STEPS
      L = FORCE(LL)
      B = FORCE(BB)

C SELECTING CORRECT SKYMAP ACCORDING TO L AND B COORDINATES
      IF ((B .GE. -10.0) .AND. (B .LT. 10.0)) THEN
        IF ((L .GE. 0.)   .AND. (L .LT. 30.))  CALL MAP258A(L,B,D,AV,A0)
        IF ((L .GE. 30.)  .AND. (L .LT. 60.))  CALL MAP258B(L,B,D,AV,A0)
        IF ((L .GE. 60.)  .AND. (L .LT. 90.))  CALL MAP258C(L,B,D,AV,A0)
        IF ((L .GE. 90.)  .AND. (L .LT. 120.)) CALL MAP259A(L,B,D,AV,A0)
        IF ((L .GE. 120.) .AND. (L .LT. 150.)) CALL MAP259B(L,B,D,AV,A0)
        IF ((L .GE. 150.) .AND. (L .LT. 180.)) CALL MAP259C(L,B,D,AV,A0)
        IF ((L .GE. 180.) .AND. (L .LT. 210.)) CALL MAP260A(L,B,D,AV,A0)
        IF ((L .GE. 210.) .AND. (L .LT. 240.)) CALL MAP260B(L,B,D,AV,A0)
        IF ((L .GE. 240.) .AND. (L .LT. 270.)) CALL MAP260C(L,B,D,AV,A0)
        IF ((L .GE. 270.) .AND. (L .LT. 300.)) CALL MAP261A(L,B,D,AV,A0)
        IF ((L .GE. 300.) .AND. (L .LT. 330.)) CALL MAP261B(L,B,D,AV,A0)
        IF ((L .GE. 330.) .AND. (L .LT. 360.)) CALL MAP261C(L,B,D,AV,A0)

c        AV = ROUND(AV,2)
        CALL ERR(D,SAV)
      ENDIF

      RETURN
      END
C ---------------------------------------------------------------

      SUBROUTINE ERR(D,SAV)
      IMPLICIT REAL*4 (A-H,L-Z)
C
C      Functional form for errors.
C
      SAV=0.214*(10**(0.1057*D))

C
C      Errors as obtained from figures, if the user prefers these.
C
C      IF ((D .GT. 0.0) .AND. (D .LE. 0.1)) SAV = 0.34
C      IF ((D .GT. 0.1) .AND. (D .LE. 0.2)) SAV = 0.17
C      IF ((D .GT. 0.2) .AND. (D .LE. 0.4)) SAV = 0.19
C      IF ((D .GT. 0.4) .AND. (D .LE. 0.6)) SAV = 0.26
C      IF ((D .GT. 0.6) .AND. (D .LE. 0.8)) SAV = 0.27
C      IF ((D .GT. 0.8) .AND. (D .LE. 1.0)) SAV = 0.19
C      IF ((D .GT. 1.0) .AND. (D .LE. 1.5)) SAV = 0.23
C      IF ((D .GT. 1.5) .AND. (D .LE. 2.0)) SAV = 0.28
C      IF ((D .GT. 2.0) .AND. (D .LE. 2.5)) SAV = 0.38
C      IF ((D .GT. 2.5) .AND. (D .LE. 3.0)) SAV = 0.36
C      IF ((D .GT. 3.0) .AND. (D .LE. 4.0)) SAV = 0.45
C      IF ((D .GT. 4.0) .AND. (D .LE. 5.0)) SAV = 0.71

      RETURN
      END

C ---------------------------------------------------------------
C THIS FUNCTION FORCES X INTO STEPS OF 0.05
      REAL FUNCTION FORCE (X)

      REAL DIGIT,X

C ROUNDING TO TWO PLACES
      X = ROUND(X,2)

C FINDING THE VALUE OF THE LAST DIGIT OF X
      DIGIT = 10*(10*X - INT(10*X))

C ROUNDING DIGIT OFF TO THE NEAREST STEP
      IF( (DIGIT .GE. 0) .AND. (DIGIT .LT. 3) ) THEN
        DIGIT = 0.0
      ELSE
        IF( (DIGIT .GE. 3) .AND. (DIGIT .LT. 8) ) THEN
          DIGIT = 0.5
        ELSE
          IF( (DIGIT .GE. 8) .AND. (DIGIT .LE. 9) ) THEN
            DIGIT = 1.0
          ENDIF
        ENDIF
      ENDIF

C COMPUTING NEW COORDINATE
      FORCE = (INT(10*X) + DIGIT)/10

      RETURN
      END

C ROUNDS THE VALUE OF TAU TO N PLACES

      REAL FUNCTION ROUND (TAU,N)

      INTEGER N

      REAL TAU,TEMP, POWER

      POWER = 10.0 ** N
      TEMP = TAU * POWER
      ROUND = NINT(TEMP)/POWER

      RETURN
      END

C ---------------------------------------------------------------

      SUBROUTINE MAP258A(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE.   0.00) .AND. (L .LT.   1.01)) .AND.
     &   ((B .GE.  -0.87) .AND. (B .LT.   0.67*L   -0.87))) THEN
        CALL ENK6(231,D,AV,A0)
        RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   1.03)) .AND.
     &   ((B .GE.  -1.89) .AND. (B .LT.  -0.87))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   1.03)) .AND.
     &   ((B .LT.  -1.89) .AND. (B .GE.   1.08*L   -2.99))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   2.20)) .AND.
     &   ((B .GE.   1.11) .AND. (B .LT.  -0.46*L +   2.18))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   2.22)) .AND.
     &   ((B .GE.  -0.21) .AND. (B .LT.   1.13))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.98) .AND. (L .LT.   2.23)) .AND.
     &   ((B .GE.  -0.87) .AND. (B .LT.  -0.18))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   1.01)) .AND.
     &   ((B .LT.  -0.23) .AND. (B .GE.   0.67*L   -0.87))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   1.00)) .AND.
     &   ((B .GE.   2.17) .AND. (B .LT.   7.74))) THEN
         CALL ENK6(235,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   1.01)) .AND.
     &   ((B .LT.   2.18) .AND. (B .GE.  -0.45*L +   2.15))) THEN
         CALL ENK6(235,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   5.66) .AND. (L .LT.   7.24)) .AND.
     &   ((B .GE.  -4.91) .AND. (B .LT.   0.63*L   -8.50))) THEN
         CALL ENK6(237,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   2.07) .AND. (L .LT.   4.99)) .AND.
     &   ((B .GE.  -4.93) .AND. (B .LT.  -0.56*L   -2.15))) THEN
         CALL ENK6(237,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   2.09)) .AND.
     &   ((B .GE.  -4.98) .AND. (B .LT.  -3.30))) THEN
         CALL ENK6(237,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   7.23)) .AND.
     &   ((B .GE.  -6.97) .AND. (B .LT.  -4.91))) THEN
         CALL ENK6(237,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.02) .AND. (L .LT.  11.87)) .AND.
     &   ((B .GE.  -6.97) .AND. (B .LT.  -1.94))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   7.19) .AND. (L .LT.  10.02)) .AND.
     &   ((B .GE.  -3.90) .AND. (B .LT.  -1.96))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   5.97)) .AND.
     &   ((B .GE.  -2.07) .AND. (B .LT.  -1.91))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   5.03) .AND. (L .LT.   5.95)) .AND.
     &   ((B .GE.  -2.12) .AND. (B .LT.   0.19*L   -3.06))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.87) .AND. (L .LT.   7.20)) .AND.
     &   ((B .GE.  -3.30) .AND. (B .LT.  -2.08))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   1.03)) .AND.
     &   ((B .GE.  -3.30) .AND. (B .LT.   1.08*L   -2.99))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   0.00) .AND. (L .LT.   0.88)) .AND.
     &   ((B .GE.  -3.32) .AND. (B .LT.  -3.01))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   2.07) .AND. (L .LT.   3.17)) .AND.
     &   ((B .LT.  -3.24) .AND. (B .GE.  -0.56*L   -2.15))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   3.16) .AND. (L .LT.   7.21)) .AND.
     &   ((B .GE.  -3.92) .AND. (B .LT.  -3.21))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   5.66) .AND. (L .LT.   7.24)) .AND.
     &   ((B .LT.  -3.89) .AND. (B .GE.   0.63*L   -8.50))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   4.98) .AND. (L .LT.   5.66)) .AND.
     &   ((B .GE.  -4.91) .AND. (B .LT.  -3.89))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   3.20) .AND. (L .LT.   4.99)) .AND.
     &   ((B .LT.  -3.92) .AND. (B .GE.  -0.56*L   -2.15))) THEN
         CALL ENK6(238,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   5.01) .AND. (L .LT.   6.52)) .AND.
     &   ((B .GE.  -1.91) .AND. (B .LT.  -0.91))) THEN
         CALL ENK6(239,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   5.03) .AND. (L .LT.   5.97)) .AND.
     &   ((B .LT.  -1.93) .AND. (B .GE.   0.16*L   -2.89))) THEN
         CALL ENK6(239,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   8.10) .AND. (L .LT.   9.24)) .AND.
     &   ((B .GE.   0.62) .AND. (B .LT.   0.35*L   -2.10))) THEN
         CALL ENK6(240,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   8.58) .AND. (L .LT.   9.25)) .AND.
     &   ((B .GE.   0.09) .AND. (B .LT.   0.70))) THEN
         CALL ENK6(240,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   6.49) .AND. (L .LT.   8.60)) .AND.
     &   ((B .GE.  -1.13) .AND. (B .LT.   0.69))) THEN
         CALL ENK6(240,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   8.60) .AND. (L .LT.   9.25)) .AND.
     &   ((B .LT.   0.12) .AND. (B .GE.   1.89*L  -17.42))) THEN
         CALL ENK6(240,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.   5.01) .AND. (L .LT.   6.51)) .AND.
     &   ((B .GE.  -0.92) .AND. (B .LT.   0.70))) THEN
         CALL ENK6(240,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.28) .AND. (L .LT.  11.68)) .AND.
     &   ((B .LT.  -0.12) .AND. (B .GE.   1.52*L  -17.89))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.92) .AND. (L .LT.  11.28)) .AND.
     &   ((B .GE.  -0.75) .AND. (B .LT.  -0.10))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.28) .AND. (L .LT.  10.92)) .AND.
     &   ((B .GE.  -0.75) .AND. (B .LT.   0.91*L  -10.07))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.26) .AND. (L .LT.  11.27)) .AND.
     &   ((B .GE.  -1.50) .AND. (B .LT.  -0.70))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.10) .AND. (L .LT.  10.35)) .AND.
     &   ((B .GE.  -1.52) .AND. (B .LT.  -0.87))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.06) .AND. (L .LT.  10.61)) .AND.
     &   ((B .GE.  -1.74) .AND. (B .LT.  -1.54))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.63) .AND. (L .LT.  10.84)) .AND.
     &   ((B .LT.  -1.54) .AND. (B .GE.   0.92*L  -11.48))) THEN
         CALL ENK7(241,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.05) .AND. (L .LT.  14.01)) .AND.
     &   ((B .GE.   1.87) .AND. (B .LT.   4.08))) THEN
         CALL ENK7(242,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.22) .AND. (L .LT.  13.08)) .AND.
     &   ((B .LT.   1.90) .AND. (B .GE.   2.02*L  -24.57))) THEN
         CALL ENK7(242,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.06) .AND. (L .LT.  12.22)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   1.90))) THEN
         CALL ENK7(242,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.28) .AND. (L .LT.  13.31)) .AND.
     &   ((B .GE.  -0.89) .AND. (B .LT.  -0.34*L +   3.62))) THEN
         CALL ENK7(243,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.64) .AND. (L .LT.  12.28)) .AND.
     &   ((B .GE.  -0.89) .AND. (B .LT.  -0.59))) THEN
         CALL ENK7(243,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.67) .AND. (L .LT.  13.29)) .AND.
     &   ((B .LT.  -0.89) .AND. (B .GE.  -0.31*L +   2.73))) THEN
         CALL ENK7(243,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.08) .AND. (L .LT.  14.75)) .AND.
     &   ((B .GE.   1.07) .AND. (B .LT.   1.89))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.69) .AND. (L .LT.  13.08)) .AND.
     &   ((B .GE.   1.07) .AND. (B .LT.   2.02*L  -24.57))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.85) .AND. (L .LT.  14.45)) .AND.
     &   ((B .GE.   0.10) .AND. (B .LT.  -1.61*L +  23.35))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.71) .AND. (L .LT.  13.89)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   1.10))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.22) .AND. (L .LT.  12.71)) .AND.
     &   ((B .GE.   0.06) .AND. (B .LT.   2.02*L  -24.57))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.29) .AND. (L .LT.  14.45)) .AND.
     &   ((B .GE.  -1.20) .AND. (B .LT.   0.08))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.65) .AND. (L .LT.  13.28)) .AND.
     &   ((B .GE.  -0.60) .AND. (B .LT.   0.09))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.28) .AND. (L .LT.  13.31)) .AND.
     &   ((B .LT.  -0.58) .AND. (B .GE.  -0.34*L +   3.62))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.20) .AND. (L .LT.  14.46)) .AND.
     &   ((B .LT.  -1.15) .AND. (B .GE.   0.28*L   -5.25))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.86) .AND. (L .LT.  12.88)) .AND.
     &   ((B .GE.  -1.33) .AND. (B .LT.   3.17*L  -42.03))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.86) .AND. (L .LT.  13.24)) .AND.
     &   ((B .LT.  -1.35) .AND. (B .GE.   0.28*L   -5.25))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.27) .AND. (L .LT.  12.86)) .AND.
     &   ((B .LT.  -1.28) .AND. (B .GE.  -0.39*L +   3.45))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.67) .AND. (L .LT.  13.29)) .AND.
     &   ((B .GE.  -1.34) .AND. (B .LT.  -0.31*L +   2.73))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.26) .AND. (L .LT.  11.67)) .AND.
     &   ((B .GE.  -1.32) .AND. (B .LT.  -0.88))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.27) .AND. (L .LT.  11.88)) .AND.
     &   ((B .GE.  -1.94) .AND. (B .LT.  -1.35))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.28) .AND. (L .LT.  11.68)) .AND.
     &   ((B .GE.  -0.88) .AND. (B .LT.   1.52*L  -17.89))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.04) .AND. (L .LT.  11.65)) .AND.
     &   ((B .GE.  -0.11) .AND. (B .LT.   0.07))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.28) .AND. (L .LT.  10.92)) .AND.
     &   ((B .LT.  -0.06) .AND. (B .GE.   0.91*L  -10.07))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  10.06) .AND. (L .LT.  10.26)) .AND.
     &   ((B .GE.  -0.88) .AND. (B .LT.  -0.11))) THEN
         CALL ENK7(244,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.65) .AND. (L .LT.  14.48)) .AND.
     &   ((B .GE.  -1.97) .AND. (B .LT.  -0.70*L +   8.14))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.86) .AND. (L .LT.  13.65)) .AND.
     &   ((B .GE.  -1.54) .AND. (B .LT.   0.28*L   -5.25))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.27) .AND. (L .LT.  12.86)) .AND.
     &   ((B .GE.  -1.57) .AND. (B .LT.  -0.39*L +   3.45))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.86) .AND. (L .LT.  12.29)) .AND.
     &   ((B .GE.  -2.95) .AND. (B .LT.  -1.34))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.29) .AND. (L .LT.  12.91)) .AND.
     &   ((B .GE.  -2.96) .AND. (B .LT.  -1.57))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  12.89) .AND. (L .LT.  13.67)) .AND.
     &   ((B .GE.  -2.94) .AND. (B .LT.  -1.54))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.65) .AND. (L .LT.  14.48)) .AND.
     &   ((B .GE.  -2.95) .AND. (B .LT.  -1.94))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.46) .AND. (L .LT.  15.07)) .AND.
     &   ((B .GE.  -2.92) .AND. (B .LT.  -1.61*L +  21.34))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.66) .AND. (L .LT.  15.07)) .AND.
     &   ((B .LT.  -2.92) .AND. (B .GE.   4.03*L  -63.58))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.45) .AND. (L .LT.  14.66)) .AND.
     &   ((B .GE.  -4.58) .AND. (B .LT.  -2.94))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.86) .AND. (L .LT.  13.53)) .AND.
     &   ((B .GE.  -3.48) .AND. (B .LT.  -2.93))) THEN
         CALL ENK7(245,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.45) .AND. (L .LT.  16.76)) .AND.
     &   ((B .GE.  -1.19) .AND. (B .LT.   0.10))) THEN
         CALL ENK7(246,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.45) .AND. (L .LT.  17.48)) .AND.
     &   ((B .LT.  -1.12) .AND. (B .GE.  -0.17*L +   1.25))) THEN
         CALL ENK7(246,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.72) .AND. (L .LT.  17.47)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.  -0.44))) THEN
         CALL ENK7(246,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  17.48) .AND. (L .LT.  18.08)) .AND.
     &   ((B .GE.  -0.91) .AND. (B .LT.  -0.44))) THEN
         CALL ENK7(246,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.05) .AND. (L .LT.  18.39)) .AND.
     &   ((B .GE.  -0.91) .AND. (B .LT.  -1.38*L +  24.52))) THEN
         CALL ENK7(246,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  15.05) .AND. (L .LT.  16.72)) .AND.
     &   ((B .GE.   3.06) .AND. (B .LT.  -0.61*L +  13.22))) THEN
         CALL ENK7(247,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.69) .AND. (L .LT.  15.06)) .AND.
     &   ((B .GE.   3.05) .AND. (B .LT.   4.09))) THEN
         CALL ENK7(247,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.02) .AND. (L .LT.  14.74)) .AND.
     &   ((B .GE.   1.87) .AND. (B .LT.   4.11))) THEN
         CALL ENK7(247,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.73) .AND. (L .LT.  16.73)) .AND.
     &   ((B .GE.   1.11) .AND. (B .LT.   3.05))) THEN
         CALL ENK7(247,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.46) .AND. (L .LT.  16.74)) .AND.
     &   ((B .GE.   0.10) .AND. (B .LT.   1.13))) THEN
         CALL ENK7(247,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.85) .AND. (L .LT.  14.45)) .AND.
     &   ((B .LT.   1.06) .AND. (B .GE.  -1.61*L +  23.35))) THEN
         CALL ENK7(247,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.07) .AND. (L .LT.  19.07)) .AND.
     &   ((B .GE.  -6.93) .AND. (B .LT.  -1.33*L +  18.49))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.12) .AND. (L .LT.  19.08)) .AND.
     &   ((B .GE.  -7.58) .AND. (B .LT.  -6.93))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  15.04) .AND. (L .LT.  16.12)) .AND.
     &   ((B .GE.  -6.94) .AND. (B .LT.  -2.95))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.10) .AND. (L .LT.  16.13)) .AND.
     &   ((B .GE.  -7.56) .AND. (B .LT.  -6.91))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.64) .AND. (L .LT.  15.07)) .AND.
     &   ((B .GE.  -4.57) .AND. (B .LT.   3.75*L  -59.38))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.86) .AND. (L .LT.  15.12)) .AND.
     &   ((B .GE.  -7.00) .AND. (B .LT.  -4.54))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  11.87) .AND. (L .LT.  13.47)) .AND.
     &   ((B .GE.  -4.55) .AND. (B .LT.  -3.46))) THEN
         CALL ENK7(248,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.75) .AND. (L .LT.  17.49)) .AND.
     &   ((B .GE.   0.45) .AND. (B .LT.   1.45))) THEN
         CALL ENK7(249,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.08) .AND. (L .LT.  21.07)) .AND.
     &   ((B .GE.   0.01) .AND. (B .LT.  -1.06*L +  22.35))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.07) .AND. (L .LT.  21.07)) .AND.
     &   ((B .GE.  -0.42) .AND. (B .LT.   0.05))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.10) .AND. (L .LT.  21.07)) .AND.
     &   ((B .GE.  -1.33) .AND. (B .LT.  -0.41))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  19.43) .AND. (L .LT.  20.09)) .AND.
     &   ((B .LT.  -0.42) .AND. (B .GE.  -1.42*L +  27.14))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  19.47) .AND. (L .LT.  20.10)) .AND.
     &   ((B .GE.  -0.43) .AND. (B .LT.   0.47))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  17.42) .AND. (L .LT.  20.09)) .AND.
     &   ((B .GE.   0.45) .AND. (B .LT.   1.14))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  19.06) .AND. (L .LT.  20.05)) .AND.
     &   ((B .GE.   1.13) .AND. (B .LT.  -2.92*L +  59.75))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.72) .AND. (L .LT.  19.06)) .AND.
     &   ((B .GE.   3.12) .AND. (B .LT.   0.44*L   -4.29))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  17.46) .AND. (L .LT.  19.10)) .AND.
     &   ((B .GE.   1.11) .AND. (B .LT.   3.10))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.75) .AND. (L .LT.  17.46)) .AND.
     &   ((B .GE.   1.45) .AND. (B .LT.   3.10))) THEN
         CALL ENK7(250,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  19.56) .AND. (L .LT.  20.09)) .AND.
     &   ((B .LT.  -1.57) .AND. (B .GE.   3.16*L  -65.09))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.95) .AND. (L .LT.  19.35)) .AND.
     &   ((B .GE.  -1.56) .AND. (B .LT.  -0.40*L +   6.25))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.40) .AND. (L .LT.  19.30)) .AND.
     &   ((B .GE.  -1.59) .AND. (B .LT.  -1.39))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.39) .AND. (L .LT.  18.71)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.  -1.38*L +  24.52))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.28) .AND. (L .LT.  19.63)) .AND.
     &   ((B .GE.  -2.91) .AND. (B .LT.  -1.60))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  17.44) .AND. (L .LT.  18.42)) .AND.
     &   ((B .GE.  -2.94) .AND. (B .LT.  -0.91))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  15.05) .AND. (L .LT.  17.49)) .AND.
     &   ((B .GE.  -2.91) .AND. (B .LT.  -1.82))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.45) .AND. (L .LT.  17.48)) .AND.
     &   ((B .GE.  -1.82) .AND. (B .LT.  -0.17*L +   1.25))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.65) .AND. (L .LT.  14.46)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.   0.48*L   -7.75))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  13.65) .AND. (L .LT.  14.48)) .AND.
     &   ((B .LT.  -1.36) .AND. (B .GE.  -0.70*L +   8.14))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.45) .AND. (L .LT.  15.20)) .AND.
     &   ((B .GE.  -1.94) .AND. (B .LT.  -1.71))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  14.48) .AND. (L .LT.  15.07)) .AND.
     &   ((B .LT.  -1.91) .AND. (B .GE.  -1.76*L +  23.51))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.11) .AND. (L .LT.  17.71)) .AND.
     &   ((B .LT.  -2.94) .AND. (B .GE.  -1.40*L +  19.61))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  17.68) .AND. (L .LT.  18.95)) .AND.
     &   ((B .GE.  -5.17) .AND. (B .LT.  -2.93))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.93) .AND. (L .LT.  19.56)) .AND.
     &   ((B .LT.  -2.92) .AND. (B .GE.   3.16*L  -65.09))) THEN
         CALL ENK7(251,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  16.74) .AND. (L .LT.  19.48)) .AND.
     &   ((B .GE.  -0.45) .AND. (B .LT.   0.49))) THEN
         CALL ENK7(252,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  19.36) .AND. (L .LT.  20.11)) .AND.
     &   ((B .GE.  -1.36) .AND. (B .LT.  -1.38*L +  26.48))) THEN
         CALL ENK7(253,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.70) .AND. (L .LT.  19.45)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.  -0.45))) THEN
         CALL ENK7(253,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.05) .AND. (L .LT.  18.71)) .AND.
     &   ((B .LT.  -0.46) .AND. (B .GE.  -1.38*L +  24.52))) THEN
         CALL ENK7(253,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.96) .AND. (L .LT.  22.06)) .AND.
     &   ((B .GE.   2.06) .AND. (B .LT.  -3.80*L +  85.89))) THEN
         CALL ENK7(254,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.01) .AND. (L .LT.  21.05)) .AND.
     &   ((B .GE.   1.14) .AND. (B .LT.   6.26))) THEN
         CALL ENK7(254,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.08) .AND. (L .LT.  21.07)) .AND.
     &   ((B .LT.   1.15) .AND. (B .GE.  -1.06*L +  22.35))) THEN
         CALL ENK7(254,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.04) .AND. (L .LT.  22.04)) .AND.
     &   ((B .LT.   2.06) .AND. (B .GE.   2.03*L  -42.75))) THEN
         CALL ENK7(254,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.09) .AND. (L .LT.  22.08)) .AND.
     &   ((B .LT.  -5.05) .AND. (B .GE.   0.07*L   -6.62))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.09) .AND. (L .LT.  22.08)) .AND.
     &   ((B .GE.  -5.05) .AND. (B .LT.  -3.57*L +  73.89))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.13) .AND. (L .LT.  22.06)) .AND.
     &   ((B .LT.  -4.97) .AND. (B .GE.   0.14*L   -8.03))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  20.09) .AND. (L .LT.  21.13)) .AND.
     &   ((B .GE.  -5.07) .AND. (B .LT.  -1.35))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  17.37) .AND. (L .LT.  18.37)) .AND.
     &   ((B .GE.  -1.60) .AND. (B .LT.  -1.38))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  19.33) .AND. (L .LT.  20.09)) .AND.
     &   ((B .LT.  -1.34) .AND. (B .GE.  -0.59*L +   9.77))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  18.93) .AND. (L .LT.  20.09)) .AND.
     &   ((B .GE.  -5.19) .AND. (B .LT.   3.16*L  -65.09))) THEN
         CALL ENK7(255,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.61) .AND. (L .LT.  25.06)) .AND.
     &   ((B .GE.  -3.40) .AND. (B .LT.   0.01))) THEN
         CALL ENK7(256,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.07) .AND. (L .LT.  21.61)) .AND.
     &   ((B .GE.  -1.35) .AND. (B .LT.   0.05))) THEN
         CALL ENK7(256,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.04) .AND. (L .LT.  21.63)) .AND.
     &   ((B .LT.  -1.34) .AND. (B .GE.  -3.40*L +  70.21))) THEN
         CALL ENK7(256,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  22.09) .AND. (L .LT.  24.04)) .AND.
     &   ((B .GE.   0.01) .AND. (B .LT.  -1.02*L +  24.62))) THEN
         CALL ENK7(257,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  21.04) .AND. (L .LT.  22.09)) .AND.
     &   ((B .GE.   0.01) .AND. (B .LT.   1.92*L  -40.32))) THEN
         CALL ENK7(257,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  25.02) .AND. (L .LT.  27.44)) .AND.
     &   ((B .GE.  -5.31) .AND. (B .LT.  -0.37))) THEN
         CALL ENK7(258,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  28.22) .AND. (L .LT.  29.06)) .AND.
     &   ((B .GE.  -3.96) .AND. (B .LT.  -4.77*L + 134.68))) THEN
         CALL ENK7(259,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  27.41) .AND. (L .LT.  28.30)) .AND.
     &   ((B .GE.  -3.97) .AND. (B .LT.   0.02))) THEN
         CALL ENK7(259,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  25.82) .AND. (L .LT.  27.40)) .AND.
     &   ((B .GE.  -0.39) .AND. (B .LT.   0.02))) THEN
         CALL ENK7(259,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  25.01) .AND. (L .LT.  30.91)) .AND.
     &   ((B .GE.   0.01) .AND. (B .LT.   7.12))) THEN
         CALL ENK7(260,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  28.22) .AND. (L .LT.  29.08)) .AND.
     &   ((B .LT.   0.03) .AND. (B .GE.  -4.59*L + 129.44))) THEN
         CALL ENK7(260,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  29.03) .AND. (L .LT.  30.96)) .AND.
     &   ((B .GE.  -3.98) .AND. (B .LT.   0.05))) THEN
         CALL ENK7(260,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  37.24) .AND. (L .LT.  38.45)) .AND.
     &   ((B .GE.  -1.93) .AND. (B .LT.  -0.50*L +  17.22))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  36.86) .AND. (L .LT.  37.24)) .AND.
     &   ((B .GE.  -1.91) .AND. (B .LT.   1.61*L  -61.19))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  36.86) .AND. (L .LT.  38.01)) .AND.
     &   ((B .LT.  -1.90) .AND. (B .GE.  -0.49*L +  16.04))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  38.01) .AND. (L .LT.  38.50)) .AND.
     &   ((B .LT.  -1.90) .AND. (B .GE.   1.18*L  -47.16))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  44.03) .AND. (L .LT.  45.08)) .AND.
     &   ((B .GE.  -2.07) .AND. (B .LT.  -1.46))) THEN
         CALL ENK7(262,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  44.06) .AND. (L .LT.  45.06)) .AND.
     &   ((B .GE.  -2.08) .AND. (B .LT.  -0.58*L +  24.20))) THEN
         CALL ENK7(262,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  40.05) .AND. (L .LT.  45.06)) .AND.
     &   ((B .GE.  -7.70) .AND. (B .LT.  -2.03))) THEN
         CALL ENK7(262,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  47.12) .AND. (L .LT.  50.10)) .AND.
     &   ((B .GE.   3.95) .AND. (B .LT.   4.38))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  45.08) .AND. (L .LT.  50.10)) .AND.
     &   ((B .GE.   0.36) .AND. (B .LT.   4.00))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  41.09) .AND. (L .LT.  45.09)) .AND.
     &   ((B .GE.   2.96) .AND. (B .LT.   3.95))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  41.08) .AND. (L .LT.  45.11)) .AND.
     &   ((B .LT.   3.01) .AND. (B .GE.  -0.49*L +  23.06))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.11) .AND. (L .LT.  57.29)) .AND.
     &   ((B .GE.   1.53) .AND. (B .LT.   5.01))) THEN
         CALL ENK7(266,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.09) .AND. (L .LT.  55.07)) .AND.
     &   ((B .GE.  -3.01) .AND. (B .LT.   1.60))) THEN
         CALL ENK7(266,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.11) .AND. (L .LT.  55.47)) .AND.
     &   ((B .GE.   0.02) .AND. (B .LT.   0.99))) THEN
         CALL ENK7(266,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.09) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.   4.95) .AND. (B .LT.   7.67))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  57.31) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.   2.04) .AND. (B .LT.   5.06))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.48) .AND. (L .LT.  58.25)) .AND.
     &   ((B .GE.   0.01) .AND. (B .LT.   1.62))) THEN
         CALL ENK7(268,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.10) .AND. (L .LT.  55.49)) .AND.
     &   ((B .GE.   1.01) .AND. (B .LT.   1.61))) THEN
         CALL ENK7(268,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  57.28) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.   1.63) .AND. (B .LT.   2.07))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  58.26) .AND. (L .LT.  60.98)) .AND.
     &   ((B .GE.   0.05) .AND. (B .LT.   1.67))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  59.44) .AND. (L .LT.  60.95)) .AND.
     &   ((B .GE.  -0.39) .AND. (B .LT.   0.09))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  58.27) .AND. (L .LT.  59.44)) .AND.
     &   ((B .LT.   0.06) .AND. (B .GE.  -0.37*L +  21.92))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.09) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.  -5.05) .AND. (B .LT.  -0.34))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  58.27) .AND. (L .LT.  59.44)) .AND.
     &   ((B .GE.  -0.36) .AND. (B .LT.  -0.33*L +  19.49))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.09) .AND. (L .LT.  58.25)) .AND.
     &   ((B .GE.  -0.43) .AND. (B .LT.   0.05))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.11) .AND. (L .LT.  60.02)) .AND.
     &   ((B .GE.  -7.58) .AND. (B .LT.  -5.00))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------
      SUBROUTINE MAP258B(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE.  29.05) .AND. (L .LT.  35.05)) .AND.
     &   ((B .GE.  -4.18) .AND. (B .LT.   6.91))) THEN
         CALL ENK7(260,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  37.24) .AND. (L .LT.  38.45)) .AND.
     &   ((B .GE.  -1.93) .AND. (B .LT.  -0.50*L +  17.22))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  36.86) .AND. (L .LT.  37.24)) .AND.
     &   ((B .GE.  -1.91) .AND. (B .LT.   1.61*L  -61.19))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  36.86) .AND. (L .LT.  38.01)) .AND.
     &   ((B .LT.  -1.90) .AND. (B .GE.  -0.49*L +  16.04))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  38.01) .AND. (L .LT.  38.50)) .AND.
     &   ((B .LT.  -1.90) .AND. (B .GE.   1.18*L  -47.16))) THEN
         CALL ENK7(261,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  44.06) .AND. (L .LT.  45.06)) .AND.
     &   ((B .GE.  -2.05) .AND. (B .LT.  -1.28))) THEN
         CALL ENK7(262,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  40.10) .AND. (L .LT.  44.06)) .AND.
     &   ((B .GE.  -2.05) .AND. (B .LT.   0.15*L   -7.89))) THEN
         CALL ENK7(262,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  40.10) .AND. (L .LT.  45.08)) .AND.
     &   ((B .GE.  -7.68) .AND. (B .LT.  -2.03))) THEN
         CALL ENK7(262,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  45.09) .AND. (L .LT.  50.10)) .AND.
     &   ((B .GE.  -3.04) .AND. (B .LT.   0.35))) THEN
         CALL ENK7(263,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  47.12) .AND. (L .LT.  50.10)) .AND.
     &   ((B .GE.   3.93) .AND. (B .LT.   4.38))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  45.08) .AND. (L .LT.  50.10)) .AND.
     &   ((B .GE.   0.36) .AND. (B .LT.   4.00))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  41.09) .AND. (L .LT.  45.09)) .AND.
     &   ((B .GE.   2.96) .AND. (B .LT.   3.95))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  41.08) .AND. (L .LT.  45.11)) .AND.
     &   ((B .LT.   3.01) .AND. (B .GE.  -0.49*L +  23.06))) THEN
         CALL ENK7(264,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  45.09) .AND. (L .LT.  54.11)) .AND.
     &   ((B .GE.  -7.66) .AND. (B .LT.  -3.00))) THEN
         CALL ENK7(265,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.11) .AND. (L .LT.  57.28)) .AND.
     &   ((B .GE.   1.53) .AND. (B .LT.   5.01))) THEN
         CALL ENK7(266,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.09) .AND. (L .LT.  55.07)) .AND.
     &   ((B .GE.  -3.01) .AND. (B .LT.   1.60))) THEN
         CALL ENK7(266,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.10) .AND. (L .LT.  55.47)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   1.01))) THEN
         CALL ENK7(266,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  50.09) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.   4.95) .AND. (B .LT.   7.67))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  57.28) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.   2.04) .AND. (B .LT.   5.06))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.48) .AND. (L .LT.  58.25)) .AND.
     &   ((B .GE.   0.01) .AND. (B .LT.   1.62))) THEN
         CALL ENK7(268,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.07) .AND. (L .LT.  55.49)) .AND.
     &   ((B .GE.   1.01) .AND. (B .LT.   1.61))) THEN
         CALL ENK7(268,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  57.28) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.   1.63) .AND. (B .LT.   2.07))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  58.20) .AND. (L .LT.  60.98)) .AND.
     &   ((B .GE.   0.05) .AND. (B .LT.   1.67))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  59.44) .AND. (L .LT.  60.98)) .AND.
     &   ((B .GE.  -0.39) .AND. (B .LT.   0.09))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  58.07) .AND. (L .LT.  59.44)) .AND.
     &   ((B .LT.   0.06) .AND. (B .GE.  -0.33*L +  19.49))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.09) .AND. (L .LT.  60.97)) .AND.
     &   ((B .GE.  -5.05) .AND. (B .LT.  -0.34))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  58.07) .AND. (L .LT.  59.44)) .AND.
     &   ((B .GE.  -0.36) .AND. (B .LT.  -0.33*L +  19.49))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.09) .AND. (L .LT.  58.25)) .AND.
     &   ((B .GE.  -0.43) .AND. (B .LT.   0.05))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  55.11) .AND. (L .LT.  60.02)) .AND.
     &   ((B .GE.  -7.58) .AND. (B .LT.  -5.00))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP258C(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE.  63.28) .AND. (L .LT.  64.58)) .AND.
     &   ((B .GE.   4.83) .AND. (B .LT.  -1.81*L + 122.00))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  59.15) .AND. (L .LT.  63.28)) .AND.
     &   ((B .GE.   4.80) .AND. (B .LT.   7.45))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  63.06) .AND. (L .LT.  64.56)) .AND.
     &   ((B .GE.   1.14) .AND. (B .LT.   4.85))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  59.16) .AND. (L .LT.  63.07)) .AND.
     &   ((B .GE.   1.77) .AND. (B .LT.   4.83))) THEN
         CALL ENK7(267,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  59.10) .AND. (L .LT.  63.05)) .AND.
     &   ((B .GE.  -0.60) .AND. (B .LT.   1.84))) THEN
         CALL ENK7(269,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  59.08) .AND. (L .LT.  63.25)) .AND.
     &   ((B .GE.  -5.23) .AND. (B .LT.  -0.56))) THEN
         CALL ENK7(270,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  66.11) .AND. (L .LT.  67.49)) .AND.
     &   ((B .GE.  -1.01) .AND. (B .LT.  -0.61*L +  40.00))) THEN
         CALL ENK7(271,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  65.09) .AND. (L .LT.  66.12)) .AND.
     &   ((B .GE.  -1.02) .AND. (B .LT.  -0.15))) THEN
         CALL ENK7(271,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  65.09) .AND. (L .LT.  67.49)) .AND.
     &   ((B .GE.  -2.15) .AND. (B .LT.  -1.00))) THEN
         CALL ENK7(271,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  63.07) .AND. (L .LT.  65.09)) .AND.
     &   ((B .GE.  -0.56) .AND. (B .LT.   1.14))) THEN
         CALL ENK7(271,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  63.25) .AND. (L .LT.  65.10)) .AND.
     &   ((B .GE.  -2.15) .AND. (B .LT.  -0.54))) THEN
         CALL ENK7(271,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  65.05) .AND. (L .LT.  71.67)) .AND.
     &   ((B .GE.  -6.20) .AND. (B .LT.  -2.11))) THEN
         CALL ENK7(272,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  63.25) .AND. (L .LT.  65.10)) .AND.
     &   ((B .GE.  -5.19) .AND. (B .LT.  -2.11))) THEN
         CALL ENK7(272,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  65.06) .AND. (L .LT.  66.25)) .AND.
     &   ((B .GE.   1.86) .AND. (B .LT.   2.85))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  66.25) .AND. (L .LT.  67.25)) .AND.
     &   ((B .GE.   0.55) .AND. (B .LT.  -1.29*L +  87.68))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  65.08) .AND. (L .LT.  66.25)) .AND.
     &   ((B .GE.   0.50) .AND. (B .LT.   1.87))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  64.59) .AND. (L .LT.  65.11)) .AND.
     &   ((B .GE.   1.13) .AND. (B .LT.   1.91))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  65.09) .AND. (L .LT.  67.25)) .AND.
     &   ((B .GE.  -0.16) .AND. (B .LT.   0.56))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  66.11) .AND. (L .LT.  67.41)) .AND.
     &   ((B .LT.  -0.14) .AND. (B .GE.  -0.61*L +  40.00))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  67.28) .AND. (L .LT.  67.47)) .AND.
     &   ((B .GE.  -0.83) .AND. (B .LT.  -4.80*L + 322.53))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  67.28) .AND. (L .LT.  67.47)) .AND.
     &   ((B .LT.  -0.83) .AND. (B .GE.  -4.97*L + 334.03))) THEN
         CALL ENK7(273,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  64.58) .AND. (L .LT.  69.70)) .AND.
     &   ((B .GE.   2.87) .AND. (B .LT.   6.53))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  63.70) .AND. (L .LT.  64.58)) .AND.
     &   ((B .LT.   6.52) .AND. (B .GE.  -1.81*L + 122.00))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  64.59) .AND. (L .LT.  65.07)) .AND.
     &   ((B .GE.   1.94) .AND. (B .LT.   2.88))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  66.25) .AND. (L .LT.  69.13)) .AND.
     &   ((B .GE.   1.93) .AND. (B .LT.   2.91))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  67.25) .AND. (L .LT.  69.13)) .AND.
     &   ((B .GE.   0.66) .AND. (B .LT.   1.95))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  66.28) .AND. (L .LT.  67.25)) .AND.
     &   ((B .LT.   1.93) .AND. (B .GE.  -1.29*L +  87.68))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  67.23) .AND. (L .LT.  68.09)) .AND.
     &   ((B .LT.   0.67) .AND. (B .GE.  -0.86*L +  58.50))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  68.09) .AND. (L .LT.  69.13)) .AND.
     &   ((B .GE.  -0.14) .AND. (B .LT.   0.68))) THEN
         CALL ENK7(274,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  70.11) .AND. (L .LT.  72.70)) .AND.
     &   ((B .GE.  -2.12) .AND. (B .LT.  -0.14))) THEN
         CALL ENK7(275,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.65) .AND. (L .LT.  75.10)) .AND.
     &   ((B .GE.  -6.85) .AND. (B .LT.  -2.12))) THEN
         CALL ENK7(276,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  80.08)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -4.75))) THEN
         CALL ENK7(276,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.16) .AND. (L .LT.  80.07)) .AND.
     &   ((B .GE.  -4.76) .AND. (B .LT.  -0.21*L +  11.99))) THEN
         CALL ENK7(276,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  76.16)) .AND.
     &   ((B .GE.  -4.75) .AND. (B .LT.  -3.67))) THEN
         CALL ENK7(276,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  76.16)) .AND.
     &   ((B .GE.  -3.66) .AND. (B .LT.  -1.45*L + 106.78))) THEN
         CALL ENK7(276,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.08) .AND. (L .LT.  72.40)) .AND.
     &   ((B .GE.   2.42) .AND. (B .LT.  -1.87*L + 137.91))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.83) .AND. (L .LT.  72.40)) .AND.
     &   ((B .LT.   2.42) .AND. (B .GE.   1.66*L  -117.66))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.83) .AND. (L .LT.  72.09)) .AND.
     &   ((B .GE.   1.15) .AND. (B .LT.  -1.27*L +  92.73))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.83) .AND. (L .LT.  72.12)) .AND.
     &   ((B .GE.   0.60) .AND. (B .LT.   1.16))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.09) .AND. (L .LT.  71.83)) .AND.
     &   ((B .GE.   0.55) .AND. (B .LT.   2.43))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  70.11) .AND. (L .LT.  71.09)) .AND.
     &   ((B .GE.   2.84) .AND. (B .LT.   4.89))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  69.10) .AND. (L .LT.  71.10)) .AND.
     &   ((B .GE.   0.58) .AND. (B .LT.   2.85))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  69.10) .AND. (L .LT.  70.44)) .AND.
     &   ((B .GE.  -0.14) .AND. (B .LT.   0.57))) THEN
         CALL ENK7(277,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.12) .AND. (L .LT.  72.57)) .AND.
     &   ((B .GE.   0.58) .AND. (B .LT.   0.71))) THEN
         CALL ENK7(278,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  70.43) .AND. (L .LT.  72.60)) .AND.
     &   ((B .GE.  -0.15) .AND. (B .LT.   0.59))) THEN
         CALL ENK7(278,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  78.09)) .AND.
     &   ((B .GE.  -2.10) .AND. (B .LT.  -0.87))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.26) .AND. (L .LT.  78.11)) .AND.
     &   ((B .GE.  -2.27) .AND. (B .LT.  -2.07))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.10) .AND. (L .LT.  75.26)) .AND.
     &   ((B .LT.  -2.10) .AND. (B .GE.  -1.00*L +  72.75))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  78.11)) .AND.
     &   ((B .GE.  -0.90) .AND. (B .LT.  -0.20*L +  15.13))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.71) .AND. (L .LT.  75.10)) .AND.
     &   ((B .GE.  -2.10) .AND. (B .LT.   0.75))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.57) .AND. (L .LT.  72.72)) .AND.
     &   ((B .GE.  -0.11) .AND. (B .LT.   0.71))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.59) .AND. (L .LT.  75.10)) .AND.
     &   ((B .GE.   0.74) .AND. (B .LT.   1.11))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  74.66) .AND. (L .LT.  75.10)) .AND.
     &   ((B .GE.   1.11) .AND. (B .LT.   0.52*L  -37.37))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.61) .AND. (L .LT.  74.09)) .AND.
     &   ((B .GE.   1.10) .AND. (B .LT.   1.51))) THEN
         CALL ENK7(279,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.35) .AND. (L .LT.  73.61)) .AND.
     &   ((B .GE.   1.36) .AND. (B .LT.  -0.81*L +  61.33))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.11) .AND. (L .LT.  73.38)) .AND.
     &   ((B .GE.   1.36) .AND. (B .LT.   1.68))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.12) .AND. (L .LT.  73.12)) .AND.
     &   ((B .GE.   1.35) .AND. (B .LT.   1.50))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.12) .AND. (L .LT.  72.57)) .AND.
     &   ((B .GE.   1.49) .AND. (B .LT.  -0.87*L +  64.76))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.83) .AND. (L .LT.  72.12)) .AND.
     &   ((B .GE.   1.49) .AND. (B .LT.   1.66*L  -117.66))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.83) .AND. (L .LT.  72.40)) .AND.
     &   ((B .LT.   1.49) .AND. (B .GE.  -1.25*L +  91.07))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.12) .AND. (L .LT.  73.62)) .AND.
     &   ((B .GE.   0.70) .AND. (B .LT.   1.36))) THEN
         CALL ENK7(280,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.57) .AND. (L .LT.  73.11)) .AND.
     &   ((B .GE.   1.50) .AND. (B .LT.   1.79))) THEN
         CALL ENK8(281,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.31) .AND. (L .LT.  72.59)) .AND.
     &   ((B .LT.   1.78) .AND. (B .GE.  -0.87*L +  64.76))) THEN
         CALL ENK8(281,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.31) .AND. (L .LT.  72.47)) .AND.
     &   ((B .GE.   1.87) .AND. (B .LT.   2.04))) THEN
         CALL ENK8(281,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.12) .AND. (L .LT.  73.12)) .AND.
     &   ((B .LT.   1.95) .AND. (B .GE.  -0.87*L +  64.76))) THEN
         CALL ENK8(281,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.12) .AND. (L .LT.  72.45)) .AND.
     &   ((B .GE.   1.95) .AND. (B .LT.   1.66*L  -117.66))) THEN
         CALL ENK8(281,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.51) .AND. (L .LT.  73.90)) .AND.
     &   ((B .GE.   2.32) .AND. (B .LT.   3.91))) THEN
         CALL ENK8(282,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.40) .AND. (L .LT.  73.70)) .AND.
     &   ((B .GE.   2.31) .AND. (B .LT.   3.90))) THEN
         CALL ENK8(282,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.63) .AND. (L .LT.  72.40)) .AND.
     &   ((B .LT.   3.92) .AND. (B .GE.  -1.87*L + 137.91))) THEN
         CALL ENK8(282,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.11) .AND. (L .LT.  73.52)) .AND.
     &   ((B .GE.   2.03) .AND. (B .LT.   2.32))) THEN
         CALL ENK8(282,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.11) .AND. (L .LT.  73.54)) .AND.
     &   ((B .LT.   2.03) .AND. (B .GE.   0.52*L  -35.89))) THEN
         CALL ENK8(282,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  72.45) .AND. (L .LT.  73.11)) .AND.
     &   ((B .GE.   1.78) .AND. (B .LT.   2.32))) THEN
         CALL ENK8(282,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.90) .AND. (L .LT.  75.09)) .AND.
     &   ((B .GE.   1.68) .AND. (B .LT.  -0.39*L +  30.94))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  74.08) .AND. (L .LT.  75.09)) .AND.
     &   ((B .GE.   1.32) .AND. (B .LT.   1.69))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.92) .AND. (L .LT.  74.11)) .AND.
     &   ((B .GE.   1.51) .AND. (B .LT.   1.68))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  74.66) .AND. (L .LT.  75.13)) .AND.
     &   ((B .LT.   1.33) .AND. (B .GE.   0.46*L  -33.58))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  74.08) .AND. (L .LT.  74.66)) .AND.
     &   ((B .GE.   1.10) .AND. (B .LT.   1.33))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.61) .AND. (L .LT.  73.92)) .AND.
     &   ((B .GE.   1.51) .AND. (B .LT.   2.14))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.52) .AND. (L .LT.  73.92)) .AND.
     &   ((B .GE.   2.12) .AND. (B .LT.   2.40))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.52) .AND. (L .LT.  73.61)) .AND.
     &   ((B .GE.   1.79) .AND. (B .LT.   2.04))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.11) .AND. (L .LT.  73.52)) .AND.
     &   ((B .GE.   1.79) .AND. (B .LT.   0.60*L  -42.41))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.11) .AND. (L .LT.  73.61)) .AND.
     &   ((B .GE.   1.62) .AND. (B .LT.   1.80))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.35) .AND. (L .LT.  73.61)) .AND.
     &   ((B .LT.   1.60) .AND. (B .GE.  -0.91*L +  68.04))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.11) .AND. (L .LT.  73.61)) .AND.
     &   ((B .GE.   1.62) .AND. (B .LT.   2.11))) THEN
         CALL ENK8(283,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.23) .AND. (L .LT.  76.60)) .AND.
     &   ((B .GE.   1.53) .AND. (B .LT.   2.30))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.87) .AND. (L .LT.  76.23)) .AND.
     &   ((B .GE.   2.06) .AND. (B .LT.   0.82*L  -60.15))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.87) .AND. (L .LT.  76.23)) .AND.
     &   ((B .GE.   1.53) .AND. (B .LT.   2.06))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  75.87)) .AND.
     &   ((B .GE.   1.52) .AND. (B .LT.   2.06))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  76.63)) .AND.
     &   ((B .LT.   1.52) .AND. (B .GE.   0.39*L  -28.25))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  75.87)) .AND.
     &   ((B .GE.   2.06) .AND. (B .LT.  -1.76*L + 135.60))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  74.29) .AND. (L .LT.  75.11)) .AND.
     &   ((B .GE.   2.19) .AND. (B .LT.   4.28))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.90) .AND. (L .LT.  75.13)) .AND.
     &   ((B .LT.   2.17) .AND. (B .GE.  -0.39*L +  30.94))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.90) .AND. (L .LT.  74.30)) .AND.
     &   ((B .GE.   2.17) .AND. (B .LT.   2.88))) THEN
         CALL ENK8(284,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.08) .AND. (L .LT.  80.16)) .AND.
     &   ((B .GE.   5.92) .AND. (B .LT.   7.54))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.14) .AND. (L .LT.  80.16)) .AND.
     &   ((B .LT.   5.95) .AND. (B .GE.   0.16*L   -6.90))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.07) .AND. (L .LT.  77.10)) .AND.
     &   ((B .GE.   3.61) .AND. (B .LT.   5.95))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.49) .AND. (L .LT.  76.72)) .AND.
     &   ((B .LT.   3.69) .AND. (B .GE.   0.89*L  -64.33))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.13) .AND. (L .LT.  75.48)) .AND.
     &   ((B .GE.   3.25) .AND. (B .LT.   3.61))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.10) .AND. (L .LT.  75.51)) .AND.
     &   ((B .LT.   3.25) .AND. (B .GE.  -1.76*L + 135.60))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  74.29) .AND. (L .LT.  75.08)) .AND.
     &   ((B .GE.   4.30) .AND. (B .LT.   5.93))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  73.87) .AND. (L .LT.  74.30)) .AND.
     &   ((B .GE.   2.90) .AND. (B .LT.   5.93))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.63) .AND. (L .LT.  73.91)) .AND.
     &   ((B .GE.   3.89) .AND. (B .LT.   5.92))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.11) .AND. (L .LT.  71.63)) .AND.
     &   ((B .GE.   4.86) .AND. (B .LT.   5.99))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  71.08) .AND. (L .LT.  71.63)) .AND.
     &   ((B .LT.   4.87) .AND. (B .GE.  -1.87*L + 137.92))) THEN
         CALL ENK8(285,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  76.60)) .AND.
     &   ((B .GE.   0.96) .AND. (B .LT.   0.39*L  -28.25))) THEN
         CALL ENK8(286,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  76.60)) .AND.
     &   ((B .GE.   0.59) .AND. (B .LT.   0.96))) THEN
         CALL ENK8(286,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.14) .AND. (L .LT.  76.60)) .AND.
     &   ((B .LT.   0.59) .AND. (B .GE.   0.22*L  -16.41))) THEN
         CALL ENK8(286,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.72) .AND. (L .LT.  77.53)) .AND.
     &   ((B .GE.   2.31) .AND. (B .LT.   3.63))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.49) .AND. (L .LT.  76.72)) .AND.
     &   ((B .GE.   2.54) .AND. (B .LT.   0.89*L  -64.33))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.63) .AND. (L .LT.  76.72)) .AND.
     &   ((B .GE.   2.32) .AND. (B .LT.   2.54))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.87) .AND. (L .LT.  76.23)) .AND.
     &   ((B .LT.   2.46) .AND. (B .GE.   0.82*L  -60.15))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.53) .AND. (L .LT.  75.87)) .AND.
     &   ((B .LT.   2.46) .AND. (B .GE.  -1.76*L + 135.60))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.60) .AND. (L .LT.  77.56)) .AND.
     &   ((B .GE.   1.61) .AND. (B .LT.   2.34))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.08) .AND. (L .LT.  77.53)) .AND.
     &   ((B .LT.   1.61) .AND. (B .GE.   0.86*L  -65.22))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.89) .AND. (L .LT.  77.08)) .AND.
     &   ((B .GE.   1.20) .AND. (B .LT.   1.61))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.60) .AND. (L .LT.  76.89)) .AND.
     &   ((B .GE.   1.21) .AND. (B .LT.   1.61))) THEN
         CALL ENK8(287,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  80.16) .AND. (L .LT.  80.73)) .AND.
     &   ((B .GE.   4.22) .AND. (B .LT.   5.96))) THEN
         CALL ENK8(288,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.10) .AND. (L .LT.  80.16)) .AND.
     &   ((B .GE.   4.22) .AND. (B .LT.   5.53))) THEN
         CALL ENK8(288,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.14) .AND. (L .LT.  80.16)) .AND.
     &   ((B .GE.   5.53) .AND. (B .LT.   0.16*L   -6.90))) THEN
         CALL ENK8(288,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.52) .AND. (L .LT.  78.50)) .AND.
     &   ((B .GE.   3.22) .AND. (B .LT.   4.22))) THEN
         CALL ENK8(288,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.10) .AND. (L .LT.  77.53)) .AND.
     &   ((B .GE.   3.63) .AND. (B .LT.   4.22))) THEN
         CALL ENK8(288,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  80.07) .AND. (L .LT.  85.09)) .AND.
     &   ((B .GE.  -6.46) .AND. (B .LT.  -4.04))) THEN
         CALL ENK8(289,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  80.09) .AND. (L .LT.  83.04)) .AND.
     &   ((B .GE.  -4.09) .AND. (B .LT.  -0.56*L +  42.37))) THEN
         CALL ENK8(289,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  78.08) .AND. (L .LT.  80.07)) .AND.
     &   ((B .GE.  -3.63) .AND. (B .LT.  -0.08))) THEN
         CALL ENK8(289,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.17) .AND. (L .LT.  78.12)) .AND.
     &   ((B .GE.  -3.67) .AND. (B .LT.  -2.29))) THEN
         CALL ENK8(289,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.17) .AND. (L .LT.  80.07)) .AND.
     &   ((B .LT.  -3.61) .AND. (B .GE.  -0.21*L +  11.99))) THEN
         CALL ENK8(289,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  76.16)) .AND.
     &   ((B .LT.  -2.29) .AND. (B .GE.  -1.45*L + 106.78))) THEN
         CALL ENK8(289,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.49) .AND. (L .LT.  78.11)) .AND.
     &   ((B .GE.   0.34) .AND. (B .LT.   0.64))) THEN
         CALL ENK8(290,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.14) .AND. (L .LT.  76.49)) .AND.
     &   ((B .GE.   0.31) .AND. (B .LT.   0.24*L  -17.74))) THEN
         CALL ENK8(290,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.14) .AND. (L .LT.  78.11)) .AND.
     &   ((B .GE.  -0.23) .AND. (B .LT.   0.33))) THEN
         CALL ENK8(290,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  75.09) .AND. (L .LT.  78.11)) .AND.
     &   ((B .LT.  -0.22) .AND. (B .GE.  -0.20*L +  15.13))) THEN
         CALL ENK8(290,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  78.93) .AND. (L .LT.  79.53)) .AND.
     &   ((B .GE.   1.34) .AND. (B .LT.  -1.28*L + 102.92))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.55) .AND. (L .LT.  78.93)) .AND.
     &   ((B .GE.   2.10) .AND. (B .LT.  -0.25*L +  21.94))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.55) .AND. (L .LT.  78.93)) .AND.
     &   ((B .GE.   1.33) .AND. (B .LT.   2.11))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  78.08) .AND. (L .LT.  79.53)) .AND.
     &   ((B .GE.  -0.10) .AND. (B .LT.   1.35))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.53) .AND. (L .LT.  78.08)) .AND.
     &   ((B .GE.   0.61) .AND. (B .LT.   1.34))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.08) .AND. (L .LT.  77.53)) .AND.
     &   ((B .GE.   1.21) .AND. (B .LT.   0.92*L  -69.37))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.49) .AND. (L .LT.  76.61)) .AND.
     &   ((B .GE.   1.42) .AND. (B .LT.   1.02*L  -76.78))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.49) .AND. (L .LT.  76.61)) .AND.
     &   ((B .GE.   1.20) .AND. (B .LT.   1.42))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  76.51) .AND. (L .LT.  77.54)) .AND.
     &   ((B .GE.   0.62) .AND. (B .LT.   1.21))) THEN
         CALL ENK8(291,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  79.55) .AND. (L .LT.  81.10)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   1.34))) THEN
         CALL ENK8(292,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  79.54) .AND. (L .LT.  82.88)) .AND.
     &   ((B .GE.   1.33) .AND. (B .LT.   4.24))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  78.93) .AND. (L .LT.  79.55)) .AND.
     &   ((B .GE.   2.11) .AND. (B .LT.   4.23))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  78.93) .AND. (L .LT.  79.59)) .AND.
     &   ((B .LT.   2.12) .AND. (B .GE.  -1.28*L + 102.92))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  78.46) .AND. (L .LT.  78.93)) .AND.
     &   ((B .GE.   2.47) .AND. (B .LT.   4.22))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.53) .AND. (L .LT.  78.48)) .AND.
     &   ((B .GE.   2.46) .AND. (B .LT.   3.23))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  77.55) .AND. (L .LT.  78.93)) .AND.
     &   ((B .LT.   2.47) .AND. (B .GE.  -0.25*L +  21.94))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  81.12) .AND. (L .LT.  82.88)) .AND.
     &   ((B .GE.  -0.05) .AND. (B .LT.   1.36))) THEN
         CALL ENK8(293,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  80.05) .AND. (L .LT.  85.64)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -6.46))) THEN
         CALL ENK8(294,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.11) .AND. (L .LT.  85.59)) .AND.
     &   ((B .GE.   0.85) .AND. (B .LT.  -1.50*L + 129.00))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  82.90) .AND. (L .LT.  85.11)) .AND.
     &   ((B .GE.   0.83) .AND. (B .LT.   1.55))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  82.90) .AND. (L .LT.  85.57)) .AND.
     &   ((B .GE.  -0.03) .AND. (B .LT.   0.86))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  83.02) .AND. (L .LT.  85.58)) .AND.
     &   ((B .GE.  -4.06) .AND. (B .LT.   0.02))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.57) .AND. (L .LT.  85.69)) .AND.
     &   ((B .GE.  -0.22) .AND. (B .LT.  -1.80*L + 153.82))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.57) .AND. (L .LT.  85.67)) .AND.
     &   ((B .GE.  -4.06) .AND. (B .LT.   0.03))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  80.06) .AND. (L .LT.  83.03)) .AND.
     &   ((B .GE.  -2.43) .AND. (B .LT.  -0.03))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  80.09) .AND. (L .LT.  83.04)) .AND.
     &   ((B .LT.  -2.41) .AND. (B .GE.  -0.56*L +  42.37))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.06) .AND. (L .LT.  85.66)) .AND.
     &   ((B .GE.  -6.46) .AND. (B .LT.  -4.03))) THEN
         CALL ENK8(295,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.37) .AND. (L .LT.  89.21)) .AND.
     &   ((B .GE.  -2.23) .AND. (B .LT.  -1.70))) THEN
         CALL ENK8(296,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.14) .AND. (L .LT.  88.38)) .AND.
     &   ((B .LT.  -1.68) .AND. (B .GE.  -2.30*L + 200.60))) THEN
         CALL ENK8(296,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.69) .AND. (L .LT.  88.38)) .AND.
     &   ((B .GE.  -2.26) .AND. (B .LT.  -0.30*L +  24.52))) THEN
         CALL ENK8(296,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.64) .AND. (L .LT.  89.25)) .AND.
     &   ((B .GE.  -7.01) .AND. (B .LT.  -2.23))) THEN
         CALL ENK8(296,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.24) .AND. (L .LT.  90.03)) .AND.
     &   ((B .GE.  -7.00) .AND. (B .LT.  -4.02))) THEN
         CALL ENK8(296,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.23) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.  -4.02) .AND. (B .LT.  -1.37*L + 119.50))) THEN
         CALL ENK8(296,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.01) .AND. (L .LT.  89.25)) .AND.
     &   ((B .GE.  -0.32) .AND. (B .LT.  -0.26*L +  22.85))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.68) .AND. (L .LT.  88.01)) .AND.
     &   ((B .GE.  -0.34) .AND. (B .LT.  -0.01))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.66) .AND. (L .LT.  85.68)) .AND.
     &   ((B .LT.  -0.03) .AND. (B .GE.   7.24*L  -620.78))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.16) .AND. (L .LT.  89.25)) .AND.
     &   ((B .GE.  -1.71) .AND. (B .LT.  -0.31))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.66) .AND. (L .LT.  88.16)) .AND.
     &   ((B .GE.  -1.42) .AND. (B .LT.  -0.32))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.69) .AND. (L .LT.  88.16)) .AND.
     &   ((B .LT.  -1.40) .AND. (B .GE.  -0.29*L +  23.51))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.16) .AND. (L .LT.  88.38)) .AND.
     &   ((B .GE.  -2.16) .AND. (B .LT.  -2.11*L + 184.09))) THEN
         CALL ENK8(297,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.60) .AND. (L .LT.  86.93)) .AND.
     &   ((B .LT.   1.03) .AND. (B .GE.   2.54*L  -219.75))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  84.38) .AND. (L .LT.  86.80)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   0.00*L +   0.02))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.65) .AND. (L .LT.  89.01)) .AND.
     &   ((B .GE.   1.17) .AND. (B .LT.  -1.93*L + 172.79))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.99) .AND. (L .LT.  88.65)) .AND.
     &   ((B .GE.   1.13) .AND. (B .LT.   1.92))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.99) .AND. (L .LT.  89.01)) .AND.
     &   ((B .GE.   0.99) .AND. (B .LT.   1.18))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.87) .AND. (L .LT.  88.87)) .AND.
     &   ((B .LT.   1.03) .AND. (B .GE.   0.02*L +   0.00))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.68) .AND. (L .LT.  88.86)) .AND.
     &   ((B .GE.   0.22) .AND. (B .LT.   1.01))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.01) .AND. (L .LT.  88.84)) .AND.
     &   ((B .LT.   0.21) .AND. (B .GE.   0.27*L  -23.73))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.68) .AND. (L .LT.  88.01)) .AND.
     &   ((B .GE.  -0.02) .AND. (B .LT.   0.28))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.56) .AND. (L .LT.  86.68)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   1.78))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.09) .AND. (L .LT.  85.58)) .AND.
     &   ((B .GE.   1.57) .AND. (B .LT.   1.77))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.11) .AND. (L .LT.  85.59)) .AND.
     &   ((B .LT.   1.58) .AND. (B .GE.  -1.50*L + 129.00))) THEN
         CALL ENK8(298,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.25) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.  -2.97) .AND. (B .LT.   0.02))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.01) .AND. (L .LT.  89.25)) .AND.
     &   ((B .LT.   0.00) .AND. (B .GE.  -0.26*L +  22.85))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.01) .AND. (L .LT.  88.84)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.   0.27*L  -23.73))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.87) .AND. (L .LT.  88.99)) .AND.
     &   ((B .GE.   0.22) .AND. (B .LT.   0.21*L  -18.62))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.87) .AND. (L .LT.  88.99)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.   0.22))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.01) .AND. (L .LT.  89.25)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.   1.18))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.01) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.   1.20) .AND. (B .LT.   1.91))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  88.65) .AND. (L .LT.  89.01)) .AND.
     &   ((B .LT.   1.92) .AND. (B .GE.  -1.93*L + 172.79))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.95) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.   1.90) .AND. (B .LT.   4.59))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  86.64) .AND. (L .LT.  86.96)) .AND.
     &   ((B .GE.   1.01) .AND. (B .LT.   4.54))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.12) .AND. (L .LT.  86.65)) .AND.
     &   ((B .GE.   1.76) .AND. (B .LT.   4.54))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  85.12) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.   4.57) .AND. (B .LT.   0.19*L  -11.99))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.25) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.  -0.53*L +  48.47))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.25) .AND. (L .LT.  90.00)) .AND.
     &   ((B .LT.   1.20) .AND. (B .GE.  -0.53*L +  48.47))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP259A(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 118.58) .AND. (L .LT. 119.60)) .AND.
     &   ((B .GE.   1.04) .AND. (B .LT.  -0.30*L +  36.85))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.59) .AND. (L .LT. 119.51)) .AND.
     &   ((B .LT.   1.04) .AND. (B .GE.   2.55*L  -303.30))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.96) .AND. (L .LT. 118.58)) .AND.
     &   ((B .GE.   0.66) .AND. (B .LT.   1.20*L  -140.95))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.65) .AND. (L .LT. 118.59)) .AND.
     &   ((B .GE.  -0.87) .AND. (B .LT.   0.68))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.98) .AND. (L .LT. 117.56)) .AND.
     &   ((B .GE.   0.68) .AND. (B .LT.   1.87))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.46) .AND. (L .LT. 116.98)) .AND.
     &   ((B .GE.   0.68) .AND. (B .LT.   0.81*L  -92.90))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.38) .AND. (L .LT. 116.68)) .AND.
     &   ((B .GE.  -1.30) .AND. (B .LT.   0.66))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.22) .AND. (L .LT. 116.37)) .AND.
     &   ((B .GE.   0.11) .AND. (B .LT.   0.68))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.46) .AND. (L .LT. 116.22)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   0.68))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.18) .AND. (L .LT. 116.35)) .AND.
     &   ((B .LT.   0.13) .AND. (B .GE.  -7.32*L + 850.36))) THEN
         CALL ENK1(  1,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.55) .AND. (L .LT. 118.58)) .AND.
     &   ((B .GE.   1.44) .AND. (B .LT.   1.87))) THEN
         CALL ENK1(  2,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.55) .AND. (L .LT. 117.96)) .AND.
     &   ((B .GE.   0.68) .AND. (B .LT.   1.44))) THEN
         CALL ENK1(  2,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.96) .AND. (L .LT. 118.58)) .AND.
     &   ((B .LT.   1.42) .AND. (B .GE.   1.20*L  -140.95))) THEN
         CALL ENK1(  2,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.23) .AND. (L .LT. 118.10)) .AND.
     &   ((B .GE.  -3.14) .AND. (B .LT.  -0.89))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.25) .AND. (L .LT. 117.82)) .AND.
     &   ((B .GE.  -3.81) .AND. (B .LT.  -3.14))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 117.80) .AND. (L .LT. 118.06)) .AND.
     &   ((B .LT.  -3.14) .AND. (B .GE.   2.59*L  -308.93))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.68) .AND. (L .LT. 117.24)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.  -0.89))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.35) .AND. (L .LT. 116.68)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.  -1.30))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.40) .AND. (L .LT. 117.25)) .AND.
     &   ((B .LT.  -2.68) .AND. (B .GE.  -1.57*L + 180.20))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.83) .AND. (L .LT. 116.50)) .AND.
     &   ((B .GE.  -1.34) .AND. (B .LT.   0.35*L  -41.90))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.83) .AND. (L .LT. 116.40)) .AND.
     &   ((B .LT.  -1.34) .AND. (B .GE.  -2.52*L + 290.65))) THEN
         CALL ENK1(  3,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.97) .AND. (L .LT. 119.98)) .AND.
     &   ((B .GE.   1.87) .AND. (B .LT.   6.65))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.45) .AND. (L .LT. 114.97)) .AND.
     &   ((B .GE.   1.85) .AND. (B .LT.   5.07))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.47) .AND. (L .LT. 115.98)) .AND.
     &   ((B .GE.   1.54) .AND. (B .LT.   1.87))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.47) .AND. (L .LT. 113.95)) .AND.
     &   ((B .GE.   1.10) .AND. (B .LT.   1.54))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.95) .AND. (L .LT. 115.00)) .AND.
     &   ((B .LT.   1.54) .AND. (B .GE.   0.36*L  -39.90))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.00) .AND. (L .LT. 115.98)) .AND.
     &   ((B .LT.   1.54) .AND. (B .GE.  -0.46*L +  54.40))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.98) .AND. (L .LT. 116.98)) .AND.
     &   ((B .LT.   1.87) .AND. (B .GE.   0.81*L  -92.90))) THEN
         CALL ENK1(  4,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.75) .AND. (L .LT. 117.25)) .AND.
     &   ((B .GE.  -3.81) .AND. (B .LT.  -1.57*L + 180.20))) THEN
         CALL ENK1(  5,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.39) .AND. (L .LT. 116.75)) .AND.
     &   ((B .GE.  -3.81) .AND. (B .LT.  -3.31))) THEN
         CALL ENK1(  5,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.36) .AND. (L .LT. 118.38)) .AND.
     &   ((B .GE.  -5.48) .AND. (B .LT.  -3.81))) THEN
         CALL ENK1(  5,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.35) .AND. (L .LT. 118.39)) .AND.
     &   ((B .LT.  -5.48) .AND. (B .GE.  -0.26*L +  24.27))) THEN
         CALL ENK1(  5,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 109.97) .AND. (L .LT. 114.39)) .AND.
     &   ((B .GE.  -7.46) .AND. (B .LT.  -3.31))) THEN
         CALL ENK1(  6,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.36) .AND. (L .LT. 116.99)) .AND.
     &   ((B .GE.  -7.46) .AND. (B .LT.  -6.17))) THEN
         CALL ENK1(  6,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.35) .AND. (L .LT. 117.00)) .AND.
     &   ((B .GE.  -6.17) .AND. (B .LT.  -0.26*L +  24.27))) THEN
         CALL ENK1(  6,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 119.57) .AND. (L .LT. 119.99)) .AND.
     &   ((B .GE.   1.09) .AND. (B .LT.   1.87))) THEN
         CALL ENK1(  7,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.57) .AND. (L .LT. 119.58)) .AND.
     &   ((B .GE.   1.40) .AND. (B .LT.   1.87))) THEN
         CALL ENK1(  7,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.58) .AND. (L .LT. 119.60)) .AND.
     &   ((B .LT.   1.40) .AND. (B .GE.  -0.30*L +  36.85))) THEN
         CALL ENK1(  7,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 119.56) .AND. (L .LT. 120.94)) .AND.
     &   ((B .GE.   0.51) .AND. (B .LT.   1.13))) THEN
         CALL ENK1(  7,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 119.99) .AND. (L .LT. 120.90)) .AND.
     &   ((B .GE.  -0.92) .AND. (B .LT.   0.51))) THEN
         CALL ENK1(  8,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.97) .AND. (L .LT. 119.99)) .AND.
     &   ((B .GE.   0.06) .AND. (B .LT.   0.44*L  -52.30))) THEN
         CALL ENK1(  8,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.97) .AND. (L .LT. 120.00)) .AND.
     &   ((B .GE.  -0.89) .AND. (B .LT.   0.06))) THEN
         CALL ENK1(  8,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.59) .AND. (L .LT. 118.97)) .AND.
     &   ((B .GE.  -0.89) .AND. (B .LT.   2.55*L  -303.30))) THEN
         CALL ENK1(  8,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 118.59) .AND. (L .LT. 119.98)) .AND.
     &   ((B .LT.  -0.89) .AND. (B .GE.  -0.45*L +  52.54))) THEN
         CALL ENK1(  8,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  90.00) .AND. (L .LT.  90.97)) .AND.
     &   ((B .GE.  -2.80) .AND. (B .LT.   0.20))) THEN
         CALL ENK8(299,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  93.85) .AND. (L .LT.  95.00)) .AND.
     &   ((B .GE.   0.20) .AND. (B .LT.   5.32))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  90.00) .AND. (L .LT.  93.85)) .AND.
     &   ((B .GE.   2.74) .AND. (B .LT.   0.67*L  -57.58))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  89.01) .AND. (L .LT.  90.00)) .AND.
     &   ((B .GE.   0.20) .AND. (B .LT.   1.42))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  90.00) .AND. (L .LT.  93.90)) .AND.
     &   ((B .GE.   0.20) .AND. (B .LT.   2.74))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  90.97) .AND. (L .LT.  94.63)) .AND.
     &   ((B .GE.  -2.47) .AND. (B .LT.   0.20))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  90.97) .AND. (L .LT.  95.00)) .AND.
     &   ((B .GE.  -2.80) .AND. (B .LT.  -2.46))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  90.00) .AND. (L .LT.  95.00)) .AND.
     &   ((B .GE.  -7.41) .AND. (B .LT.  -2.79))) THEN
         CALL ENK8(300,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.02) .AND. (L .LT. 100.64)) .AND.
     &   ((B .GE.  -4.88) .AND. (B .LT.  -2.47))) THEN
         CALL ENK8(301,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.05) .AND. (L .LT. 101.06)) .AND.
     &   ((B .GE.  -6.45) .AND. (B .LT.  -4.88))) THEN
         CALL ENK8(301,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  94.99) .AND. (L .LT. 100.05)) .AND.
     &   ((B .GE.  -6.45) .AND. (B .LT.  -2.47))) THEN
         CALL ENK8(301,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.05) .AND. (L .LT. 101.06)) .AND.
     &   ((B .LT.  -6.45) .AND. (B .GE.  -0.99*L +  92.55))) THEN
         CALL ENK8(301,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.04) .AND. (L .LT. 103.56)) .AND.
     &   ((B .GE.  -7.46) .AND. (B .LT.  -4.87))) THEN
         CALL ENK8(301,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.62) .AND. (L .LT. 103.54)) .AND.
     &   ((B .GE.  -4.88) .AND. (B .LT.  -0.83*L +  81.06))) THEN
         CALL ENK8(301,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.56) .AND. (L .LT. 101.31)) .AND.
     &   ((B .GE.  -1.87) .AND. (B .LT.  -0.84))) THEN
         CALL ENK8(302,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  97.00) .AND. (L .LT. 100.62)) .AND.
     &   ((B .GE.  -2.47) .AND. (B .LT.  -0.84))) THEN
         CALL ENK8(302,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.51) .AND. (L .LT. 105.04)) .AND.
     &   ((B .GE.  -3.32) .AND. (B .LT.  -2.65))) THEN
         CALL ENK8(303,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.96) .AND. (L .LT. 103.54)) .AND.
     &   ((B .LT.  -2.68) .AND. (B .GE.  -0.83*L +  81.06))) THEN
         CALL ENK8(303,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.51) .AND. (L .LT. 106.00)) .AND.
     &   ((B .GE.  -5.25) .AND. (B .LT.  -3.32))) THEN
         CALL ENK8(303,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 105.03) .AND. (L .LT. 105.97)) .AND.
     &   ((B .GE.  -3.32) .AND. (B .LT.  -0.71*L +  72.00))) THEN
         CALL ENK8(303,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.07) .AND. (L .LT. 103.43)) .AND.
     &   ((B .GE.  -1.48) .AND. (B .LT.  -2.04*L + 209.60))) THEN
         CALL ENK8(304,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.30) .AND. (L .LT. 103.07)) .AND.
     &   ((B .GE.  -1.49) .AND. (B .LT.  -0.74))) THEN
         CALL ENK8(304,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.31) .AND. (L .LT. 104.20)) .AND.
     &   ((B .GE.  -2.66) .AND. (B .LT.  -1.46))) THEN
         CALL ENK8(304,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.59) .AND. (L .LT. 101.31)) .AND.
     &   ((B .GE.  -2.44) .AND. (B .LT.  -1.87))) THEN
         CALL ENK8(304,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.95) .AND. (L .LT. 101.31)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.  -2.44))) THEN
         CALL ENK8(304,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.62) .AND. (L .LT. 100.96)) .AND.
     &   ((B .LT.  -2.44) .AND. (B .GE.  -0.83*L +  81.06))) THEN
         CALL ENK8(304,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.00) .AND. (L .LT. 101.53)) .AND.
     &   ((B .GE.   0.58) .AND. (B .LT.   1.95))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.53) .AND. (L .LT. 102.51)) .AND.
     &   ((B .GE.   0.59) .AND. (B .LT.   1.21))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.48) .AND. (L .LT. 103.63)) .AND.
     &   ((B .GE.   0.59) .AND. (B .LT.  -0.59*L +  61.72))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.63) .AND. (L .LT. 104.18)) .AND.
     &   ((B .GE.   0.59) .AND. (B .LT.   0.66*L  -67.73))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.00) .AND. (L .LT. 104.21)) .AND.
     &   ((B .GE.   0.11) .AND. (B .LT.   0.59))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 100.00) .AND. (L .LT. 101.28)) .AND.
     &   ((B .GE.  -0.89) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.30) .AND. (L .LT. 102.28)) .AND.
     &   ((B .GE.  -1.49) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.30) .AND. (L .LT. 103.04)) .AND.
     &   ((B .GE.  -0.75) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.04) .AND. (L .LT. 104.78)) .AND.
     &   ((B .GE.  -0.53) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.19) .AND. (L .LT. 104.78)) .AND.
     &   ((B .LT.  -0.53) .AND. (B .GE.   1.62*L  -170.32))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.40) .AND. (L .LT. 104.19)) .AND.
     &   ((B .GE.  -1.48) .AND. (B .LT.  -0.53))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.04) .AND. (L .LT. 103.40)) .AND.
     &   ((B .GE.  -0.74) .AND. (B .LT.  -0.53))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.07) .AND. (L .LT. 103.43)) .AND.
     &   ((B .LT.  -0.74) .AND. (B .GE.  -2.04*L + 209.60))) THEN
         CALL ENK8(305,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  97.00) .AND. (L .LT. 100.00)) .AND.
     &   ((B .GE.  -0.84) .AND. (B .LT.   0.20))) THEN
         CALL ENK8(306,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  94.64) .AND. (L .LT.  97.00)) .AND.
     &   ((B .GE.  -2.46) .AND. (B .LT.   0.20))) THEN
         CALL ENK8(306,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  95.00) .AND. (L .LT. 100.00)) .AND.
     &   ((B .GE.   0.20) .AND. (B .LT.   2.12))) THEN
         CALL ENK8(307,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  99.20) .AND. (L .LT.  99.49)) .AND.
     &   ((B .GE.   6.12) .AND. (B .LT.  -5.24*L + 527.30))) THEN
         CALL ENK8(308,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  95.00) .AND. (L .LT.  99.21)) .AND.
     &   ((B .GE.   6.11) .AND. (B .LT.   7.71))) THEN
         CALL ENK8(308,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  95.00) .AND. (L .LT. 100.00)) .AND.
     &   ((B .GE.   2.12) .AND. (B .LT.   6.11))) THEN
         CALL ENK8(308,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.99) .AND. (L .LT. 103.99)) .AND.
     &   ((B .GE.   7.30) .AND. (B .LT.   7.71))) THEN
         CALL ENK8(309,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.00) .AND. (L .LT. 103.99)) .AND.
     &   ((B .LT.   7.30) .AND. (B .GE.   0.60*L  -55.05))) THEN
         CALL ENK8(309,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  99.48) .AND. (L .LT. 102.00)) .AND.
     &   ((B .GE.   6.10) .AND. (B .LT.   7.71))) THEN
         CALL ENK8(309,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE.  99.20) .AND. (L .LT.  99.49)) .AND.
     &   ((B .LT.   7.71) .AND. (B .GE.  -5.24*L + 527.30))) THEN
         CALL ENK8(309,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.17) .AND. (L .LT. 104.48)) .AND.
     &   ((B .GE.   1.53) .AND. (B .LT.   2.43))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.20) .AND. (L .LT. 104.46)) .AND.
     &   ((B .LT.   1.53) .AND. (B .GE.   2.12*L  -220.01))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.59) .AND. (L .LT. 104.20)) .AND.
     &   ((B .GE.   0.95) .AND. (B .LT.   2.43))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.63) .AND. (L .LT. 104.18)) .AND.
     &   ((B .LT.   0.95) .AND. (B .GE.   0.66*L  -67.73))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.21) .AND. (L .LT. 103.60)) .AND.
     &   ((B .GE.   1.21) .AND. (B .LT.   2.43))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.50) .AND. (L .LT. 102.21)) .AND.
     &   ((B .GE.   1.90) .AND. (B .LT.   0.74*L  -73.19))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.52) .AND. (L .LT. 102.21)) .AND.
     &   ((B .GE.   1.23) .AND. (B .LT.   1.90))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.48) .AND. (L .LT. 103.63)) .AND.
     &   ((B .LT.   1.24) .AND. (B .GE.  -0.59*L +  61.72))) THEN
         CALL ENK8(310,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.00) .AND. (L .LT. 102.73)) .AND.
     &   ((B .GE.   6.10) .AND. (B .LT.   0.60*L  -55.05))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.00) .AND. (L .LT. 102.75)) .AND.
     &   ((B .GE.   3.52) .AND. (B .LT.   6.10))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.73) .AND. (L .LT. 103.35)) .AND.
     &   ((B .GE.   3.53) .AND. (B .LT.  -4.89*L + 508.63))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.00) .AND. (L .LT. 103.33)) .AND.
     &   ((B .GE.   2.40) .AND. (B .LT.   3.53))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.32) .AND. (L .LT. 104.48)) .AND.
     &   ((B .GE.   2.43) .AND. (B .LT.  -0.93*L +  99.69))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.50) .AND. (L .LT. 102.21)) .AND.
     &   ((B .LT.   2.43) .AND. (B .GE.   0.74*L  -73.19))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 101.00) .AND. (L .LT. 101.50)) .AND.
     &   ((B .GE.   1.90) .AND. (B .LT.   2.43))) THEN
         CALL ENK8(311,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.97) .AND. (L .LT. 104.97)) .AND.
     &   ((B .GE.   3.53) .AND. (B .LT.   7.30))) THEN
         CALL ENK8(312,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.73) .AND. (L .LT. 103.99)) .AND.
     &   ((B .GE.   6.54) .AND. (B .LT.   0.60*L  -55.05))) THEN
         CALL ENK8(312,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 102.73) .AND. (L .LT. 103.35)) .AND.
     &   ((B .LT.   6.54) .AND. (B .GE.  -4.89*L + 508.63))) THEN
         CALL ENK8(312,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 103.33) .AND. (L .LT. 104.02)) .AND.
     &   ((B .GE.   3.53) .AND. (B .LT.   6.54))) THEN
         CALL ENK8(312,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 109.00) .AND. (L .LT. 109.97)) .AND.
     &   ((B .GE.  -3.27) .AND. (B .LT.  -3.50*L + 381.61))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 105.95) .AND. (L .LT. 109.02)) .AND.
     &   ((B .GE.  -3.27) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.96) .AND. (L .LT. 106.97)) .AND.
     &   ((B .GE.   0.11) .AND. (B .LT.  -0.43*L +  46.07))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.18) .AND. (L .LT. 104.97)) .AND.
     &   ((B .GE.   0.12) .AND. (B .LT.   0.95))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.80) .AND. (L .LT. 105.96)) .AND.
     &   ((B .GE.  -1.48) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.19) .AND. (L .LT. 104.78)) .AND.
     &   ((B .GE.  -1.48) .AND. (B .LT.   1.62*L  -170.32))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 105.03) .AND. (L .LT. 105.97)) .AND.
     &   ((B .LT.  -2.65) .AND. (B .GE.  -0.71*L +  72.00))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 104.21) .AND. (L .LT. 105.96)) .AND.
     &   ((B .GE.  -2.67) .AND. (B .LT.  -1.48))) THEN
         CALL ENK8(313,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 109.59) .AND. (L .LT. 110.94)) .AND.
     &   ((B .GE.  -1.90) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(314,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 109.00) .AND. (L .LT. 109.61)) .AND.
     &   ((B .LT.   0.13) .AND. (B .GE.  -3.50*L + 381.61))) THEN
         CALL ENK8(314,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.73) .AND. (L .LT. 114.00)) .AND.
     &   ((B .GE.   0.49) .AND. (B .LT.   1.10))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 110.73) .AND. (L .LT. 111.74)) .AND.
     &   ((B .GE.   0.13) .AND. (B .LT.   1.10))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.74) .AND. (L .LT. 112.79)) .AND.
     &   ((B .LT.   0.49) .AND. (B .GE.   0.34*L  -37.90))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 110.94) .AND. (L .LT. 111.74)) .AND.
     &   ((B .GE.  -0.46) .AND. (B .LT.   0.11))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.74) .AND. (L .LT. 112.25)) .AND.
     &   ((B .GE.  -0.46) .AND. (B .LT.  -1.13*L + 126.30))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.37) .AND. (L .LT. 112.20)) .AND.
     &   ((B .LT.  -0.46) .AND. (B .GE.   0.61*L  -69.00))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 110.94) .AND. (L .LT. 111.39)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.  -0.46))) THEN
         CALL ENK8(315,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 108.96) .AND. (L .LT. 111.47)) .AND.
     &   ((B .GE.   2.78) .AND. (B .LT.   5.14))) THEN
         CALL ENK8(316,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 109.96) .AND. (L .LT. 111.47)) .AND.
     &   ((B .GE.   1.10) .AND. (B .LT.   2.78))) THEN
         CALL ENK8(316,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 108.96) .AND. (L .LT. 110.01)) .AND.
     &   ((B .LT.   2.80) .AND. (B .GE.  -0.75*L +  84.30))) THEN
         CALL ENK8(316,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.37) .AND. (L .LT. 112.20)) .AND.
     &   ((B .GE.  -0.94) .AND. (B .LT.   0.61*L  -69.00))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.24) .AND. (L .LT. 114.18)) .AND.
     &   ((B .GE.  -0.94) .AND. (B .LT.   0.13))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.74) .AND. (L .LT. 112.25)) .AND.
     &   ((B .LT.   0.13) .AND. (B .GE.  -1.13*L + 126.30))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.18) .AND. (L .LT. 115.01)) .AND.
     &   ((B .LT.   0.13) .AND. (B .GE.   1.21*L  -139.10))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.74) .AND. (L .LT. 112.79)) .AND.
     &   ((B .GE.   0.13) .AND. (B .LT.   0.34*L  -37.90))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.79) .AND. (L .LT. 113.98)) .AND.
     &   ((B .GE.   0.13) .AND. (B .LT.   0.49))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.95) .AND. (L .LT. 115.03)) .AND.
     &   ((B .GE.   0.13) .AND. (B .LT.   1.09))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.00) .AND. (L .LT. 115.46)) .AND.
     &   ((B .GE.   0.13) .AND. (B .LT.   1.09))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.46) .AND. (L .LT. 115.96)) .AND.
     &   ((B .LT.   1.09) .AND. (B .GE.   0.81*L  -92.90))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.00) .AND. (L .LT. 115.98)) .AND.
     &   ((B .GE.   1.09) .AND. (B .LT.  -0.46*L +  54.40))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.95) .AND. (L .LT. 115.00)) .AND.
     &   ((B .GE.   1.09) .AND. (B .LT.   0.36*L  -39.90))) THEN
         CALL ENK8(317,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.25) .AND. (L .LT. 112.99)) .AND.
     &   ((B .GE.  -1.25) .AND. (B .LT.  -0.89))) THEN
         CALL ENK8(318,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.37) .AND. (L .LT. 112.23)) .AND.
     &   ((B .GE.  -1.25) .AND. (B .LT.  -0.94))) THEN
         CALL ENK8(318,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.39) .AND. (L .LT. 111.78)) .AND.
     &   ((B .GE.  -1.59) .AND. (B .LT.  -1.25))) THEN
         CALL ENK8(318,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.78) .AND. (L .LT. 112.96)) .AND.
     &   ((B .LT.  -1.28) .AND. (B .GE.   1.19*L  -135.71))) THEN
         CALL ENK8(318,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 110.94) .AND. (L .LT. 111.76)) .AND.
     &   ((B .GE.  -3.31) .AND. (B .LT.  -1.56))) THEN
         CALL ENK8(318,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.96) .AND. (L .LT. 113.58)) .AND.
     &   ((B .GE.  -1.06) .AND. (B .LT.  -0.87))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.99) .AND. (L .LT. 113.16)) .AND.
     &   ((B .GE.  -1.06) .AND. (B .LT.  -0.87))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.14) .AND. (L .LT. 113.58)) .AND.
     &   ((B .LT.  -1.06) .AND. (B .GE.   1.70*L  -194.10))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.99) .AND. (L .LT. 113.16)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.  -1.08))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.56) .AND. (L .LT. 112.96)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.   1.19*L  -135.71))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.56) .AND. (L .LT. 113.14)) .AND.
     &   ((B .GE.  -2.30) .AND. (B .LT.  -1.83))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 112.78) .AND. (L .LT. 113.31)) .AND.
     &   ((B .LT.  -2.28) .AND. (B .GE.  -0.32*L +  33.74))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.14) .AND. (L .LT. 113.31)) .AND.
     &   ((B .GE.  -2.28) .AND. (B .LT.  -2.95*L + 332.00))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.78) .AND. (L .LT. 112.56)) .AND.
     &   ((B .GE.  -2.66) .AND. (B .LT.   1.19*L  -135.71))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 111.78) .AND. (L .LT. 112.59)) .AND.
     &   ((B .LT.  -2.66) .AND. (B .GE.  -0.24*L +  24.12))) THEN
         CALL ENK8(319,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.82) .AND. (L .LT. 116.22)) .AND.
     &   ((B .GE.  -0.08) .AND. (B .LT.  -0.47*L +  54.44))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.84) .AND. (L .LT. 116.23)) .AND.
     &   ((B .GE.  -1.16) .AND. (B .LT.  -0.11))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.22) .AND. (L .LT. 116.38)) .AND.
     &   ((B .GE.  -1.18) .AND. (B .LT.  -7.32*L + 850.36))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.83) .AND. (L .LT. 116.50)) .AND.
     &   ((B .LT.  -1.18) .AND. (B .GE.   0.35*L  -41.90))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.83) .AND. (L .LT. 116.43)) .AND.
     &   ((B .GE.  -2.83) .AND. (B .LT.  -2.52*L + 290.65))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.84) .AND. (L .LT. 116.46)) .AND.
     &   ((B .GE.  -3.31) .AND. (B .LT.  -2.83))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 116.43) .AND. (L .LT. 116.78)) .AND.
     &   ((B .GE.  -3.31) .AND. (B .LT.  -1.57*L + 180.20))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 115.01) .AND. (L .LT. 115.86)) .AND.
     &   ((B .GE.  -3.31) .AND. (B .LT.   0.13))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 114.18) .AND. (L .LT. 115.01)) .AND.
     &   ((B .GE.  -0.97) .AND. (B .LT.   1.21*L  -139.10))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.58) .AND. (L .LT. 113.58)) .AND.
     &   ((B .GE.  -0.99) .AND. (B .LT.  -0.87))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.58) .AND. (L .LT. 115.02)) .AND.
     &   ((B .GE.  -1.85) .AND. (B .LT.  -0.94))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.14) .AND. (L .LT. 113.58)) .AND.
     &   ((B .GE.  -1.85) .AND. (B .LT.   1.70*L  -194.10))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.14) .AND. (L .LT. 113.31)) .AND.
     &   ((B .LT.  -1.85) .AND. (B .GE.  -2.95*L + 332.00))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 113.31) .AND. (L .LT. 115.03)) .AND.
     &   ((B .GE.  -3.31) .AND. (B .LT.  -1.83))) THEN
         CALL ENK8(320,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP259B(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 120.00) .AND. (L .LT. 121.62)) .AND.
     &   ((B .GE.   0.68) .AND. (B .LT.   1.26))) THEN
         CALL ENK1(  7,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 120.00) .AND. (L .LT. 121.62)) .AND.
     &   ((B .GE.  -1.00) .AND. (B .LT.   0.68))) THEN
         CALL ENK1(  8,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 120.00) .AND. (L .LT. 121.62)) .AND.
     &   ((B .GE.   1.14) .AND. (B .LT.   2.48))) THEN
         CALL ENK1(  9,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.62) .AND. (L .LT. 121.81)) .AND.
     &   ((B .GE.   1.26) .AND. (B .LT.   2.51))) THEN
         CALL ENK1(  9,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.25) .AND. (L .LT. 125.07)) .AND.
     &   ((B .GE.  -1.52) .AND. (B .LT.   0.48*L  -60.79))) THEN
         CALL ENK1( 10,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 120.00) .AND. (L .LT. 125.07)) .AND.
     &   ((B .GE.  -3.30) .AND. (B .LT.  -1.52))) THEN
         CALL ENK1( 10,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.98) .AND. (L .LT. 125.62)) .AND.
     &   ((B .GE.   0.06) .AND. (B .LT.   0.52))) THEN
         CALL ENK1( 11,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 125.00)) .AND.
     &   ((B .GE.   0.04) .AND. (B .LT.   0.54*L  -66.94))) THEN
         CALL ENK1( 11,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.00) .AND. (L .LT. 125.62)) .AND.
     &   ((B .LT.   0.06) .AND. (B .GE.   1.36*L  -170.89))) THEN
         CALL ENK1( 11,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.25) .AND. (L .LT. 124.98)) .AND.
     &   ((B .GE.  -0.71) .AND. (B .LT.   0.08))) THEN
         CALL ENK1( 11,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.62) .AND. (L .LT. 123.24)) .AND.
     &   ((B .GE.  -1.52) .AND. (B .LT.   0.08))) THEN
         CALL ENK1( 11,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.25) .AND. (L .LT. 124.98)) .AND.
     &   ((B .LT.  -0.71) .AND. (B .GE.   0.48*L  -60.79))) THEN
         CALL ENK1( 11,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.83) .AND. (L .LT. 123.43)) .AND.
     &   ((B .GE.   0.88) .AND. (B .LT.  -0.28*L +  35.52))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.83) .AND. (L .LT. 123.41)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   0.88))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.21) .AND. (L .LT. 122.86)) .AND.
     &   ((B .GE.   0.40) .AND. (B .LT.   1.05*L  -127.83))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.19) .AND. (L .LT. 122.86)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   0.40))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.40) .AND. (L .LT. 124.00)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   0.92))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 124.48)) .AND.
     &   ((B .GE.   0.20) .AND. (B .LT.  -0.56*L +  69.86))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 124.48)) .AND.
     &   ((B .LT.   0.20) .AND. (B .GE.   0.54*L  -66.94))) THEN
         CALL ENK1( 12,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.21) .AND. (L .LT. 122.83)) .AND.
     &   ((B .GE.   1.09) .AND. (B .LT.  -0.28*L +  35.52))) THEN
         CALL ENK1( 13,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.21) .AND. (L .LT. 122.86)) .AND.
     &   ((B .LT.   1.09) .AND. (B .GE.   1.05*L  -127.83))) THEN
         CALL ENK1( 13,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.62) .AND. (L .LT. 122.21)) .AND.
     &   ((B .GE.   0.40) .AND. (B .LT.   1.31))) THEN
         CALL ENK1( 13,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.62) .AND. (L .LT. 122.21)) .AND.
     &   ((B .LT.   0.40) .AND. (B .GE.   0.53*L  -62.19))) THEN
         CALL ENK1( 13,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.22) .AND. (L .LT. 124.25)) .AND.
     &   ((B .GE.   3.90) .AND. (B .LT.  -0.22*L +  31.35))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.82) .AND. (L .LT. 123.22)) .AND.
     &   ((B .GE.   4.02) .AND. (B .LT.  -0.22*L +  31.35))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.82) .AND. (L .LT. 123.23)) .AND.
     &   ((B .GE.   2.12) .AND. (B .LT.   4.02))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.23) .AND. (L .LT. 124.25)) .AND.
     &   ((B .LT.   3.90) .AND. (B .GE.   2.35*L  -287.91))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 121.85) .AND. (L .LT. 123.22)) .AND.
     &   ((B .GE.   1.33) .AND. (B .LT.   2.10))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 122.16) .AND. (L .LT. 123.43)) .AND.
     &   ((B .LT.   1.31) .AND. (B .GE.  -0.28*L +  35.52))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.25) .AND. (L .LT. 124.40)) .AND.
     &   ((B .GE.   1.45) .AND. (B .LT.  -0.59*L +  74.75))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.43) .AND. (L .LT. 124.00)) .AND.
     &   ((B .GE.   0.90) .AND. (B .LT.   1.45))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 124.40)) .AND.
     &   ((B .LT.   1.45) .AND. (B .GE.   1.31*L  -161.51))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.14) .AND. (L .LT. 123.43)) .AND.
     &   ((B .GE.   1.31) .AND. (B .LT.   1.45))) THEN
         CALL ENK1( 14,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.39) .AND. (L .LT. 125.60)) .AND.
     &   ((B .GE.   0.92) .AND. (B .LT.   3.27))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.75) .AND. (L .LT. 124.40)) .AND.
     &   ((B .GE.   2.12) .AND. (B .LT.   3.30))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.23) .AND. (L .LT. 123.99)) .AND.
     &   ((B .GE.   2.12) .AND. (B .LT.   2.35*L  -287.91))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 123.25) .AND. (L .LT. 124.40)) .AND.
     &   ((B .LT.   2.12) .AND. (B .GE.  -0.59*L +  74.75))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 124.40)) .AND.
     &   ((B .GE.   0.92) .AND. (B .LT.   1.31*L  -161.51))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 125.62)) .AND.
     &   ((B .GE.   0.52) .AND. (B .LT.   0.92))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.00) .AND. (L .LT. 124.48)) .AND.
     &   ((B .LT.   0.52) .AND. (B .GE.  -0.56*L +  69.86))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.48) .AND. (L .LT. 125.00)) .AND.
     &   ((B .LT.   0.52) .AND. (B .GE.   0.54*L  -66.94))) THEN
         CALL ENK1( 15,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.58) .AND. (L .LT. 126.38)) .AND.
     &   ((B .GE.   1.57) .AND. (B .LT.   5.10))) THEN
         CALL ENK1( 16,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.00) .AND. (L .LT. 125.58)) .AND.
     &   ((B .GE.   3.28) .AND. (B .LT.   5.12))) THEN
         CALL ENK1( 16,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 128.23) .AND. (L .LT. 129.25)) .AND.
     &   ((B .GE.  -1.40) .AND. (B .LT.  -0.61*L +  77.49))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.06) .AND. (L .LT. 128.23)) .AND.
     &   ((B .GE.  -1.74) .AND. (B .LT.   0.40*L  -52.03))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.07) .AND. (L .LT. 126.04)) .AND.
     &   ((B .GE.  -2.91) .AND. (B .LT.   1.21*L  -154.11))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 129.25) .AND. (L .LT. 130.01)) .AND.
     &   ((B .GE.  -4.90) .AND. (B .LT.  -1.33))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 128.20) .AND. (L .LT. 129.26)) .AND.
     &   ((B .GE.  -2.91) .AND. (B .LT.  -1.40))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.04) .AND. (L .LT. 128.24)) .AND.
     &   ((B .GE.  -2.94) .AND. (B .LT.  -1.71))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.07) .AND. (L .LT. 129.29)) .AND.
     &   ((B .GE.  -4.90) .AND. (B .LT.  -2.91))) THEN
         CALL ENK1( 17,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 129.40) .AND. (L .LT. 130.00)) .AND.
     &   ((B .GE.   1.07) .AND. (B .LT.   2.31))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 127.67) .AND. (L .LT. 129.43)) .AND.
     &   ((B .GE.   1.07) .AND. (B .LT.   0.71*L  -89.56))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.81) .AND. (L .LT. 130.01)) .AND.
     &   ((B .GE.  -0.76) .AND. (B .LT.   1.07))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 128.23) .AND. (L .LT. 129.25)) .AND.
     &   ((B .LT.  -0.78) .AND. (B .GE.  -0.61*L +  77.49))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 129.25) .AND. (L .LT. 130.01)) .AND.
     &   ((B .GE.  -1.33) .AND. (B .LT.  -0.78))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.06) .AND. (L .LT. 128.26)) .AND.
     &   ((B .LT.  -0.73) .AND. (B .GE.   0.40*L  -52.03))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.62) .AND. (L .LT. 126.89)) .AND.
     &   ((B .GE.  -0.76) .AND. (B .LT.   0.11))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.01) .AND. (L .LT. 125.62)) .AND.
     &   ((B .GE.  -0.76) .AND. (B .LT.   1.36*L  -170.89))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.03) .AND. (L .LT. 126.06)) .AND.
     &   ((B .LT.  -0.76) .AND. (B .GE.  -0.91*L + 112.89))) THEN
         CALL ENK1( 18,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 127.69) .AND. (L .LT. 128.48)) .AND.
     &   ((B .LT.   1.62) .AND. (B .GE.   0.71*L  -89.56))) THEN
         CALL ENK1( 19,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.81) .AND. (L .LT. 127.67)) .AND.
     &   ((B .GE.   1.09) .AND. (B .LT.   1.62))) THEN
         CALL ENK1( 19,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 125.59) .AND. (L .LT. 126.91)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   1.62))) THEN
         CALL ENK1( 19,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.39) .AND. (L .LT. 135.00)) .AND.
     &   ((B .GE.   3.25) .AND. (B .LT.   7.66))) THEN
         CALL ENK1( 20,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 124.01) .AND. (L .LT. 126.39)) .AND.
     &   ((B .GE.   5.10) .AND. (B .LT.   7.66))) THEN
         CALL ENK1( 20,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 126.39) .AND. (L .LT. 128.50)) .AND.
     &   ((B .GE.   1.62) .AND. (B .LT.   3.25))) THEN
         CALL ENK1( 20,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 128.49) .AND. (L .LT. 130.00)) .AND.
     &   ((B .GE.   2.31) .AND. (B .LT.   3.25))) THEN
         CALL ENK1( 20,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 128.48) .AND. (L .LT. 129.43)) .AND.
     &   ((B .LT.   2.34) .AND. (B .GE.   0.71*L  -89.56))) THEN
         CALL ENK1( 20,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.16) .AND. (L .LT. 134.98)) .AND.
     &   ((B .LT.   3.27) .AND. (B .GE.  -0.31*L +  44.00))) THEN
         CALL ENK1( 20,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 129.01) .AND. (L .LT. 137.98)) .AND.
     &   ((B .GE.  -7.54) .AND. (B .LT.  -5.34))) THEN
         CALL ENK1( 21,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 134.10) .AND. (L .LT. 135.04)) .AND.
     &   ((B .LT.  -3.49) .AND. (B .GE.  -0.70*L +  90.41))) THEN
         CALL ENK1( 22,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 135.04) .AND. (L .LT. 135.46)) .AND.
     &   ((B .GE.  -4.06) .AND. (B .LT.  -3.63))) THEN
         CALL ENK1( 22,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 135.46) .AND. (L .LT. 137.97)) .AND.
     &   ((B .GE.  -5.36) .AND. (B .LT.  -3.15))) THEN
         CALL ENK1( 23,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.46) .AND. (L .LT. 135.46)) .AND.
     &   ((B .GE.  -5.36) .AND. (B .LT.  -4.16))) THEN
         CALL ENK1( 23,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 134.15) .AND. (L .LT. 135.50)) .AND.
     &   ((B .GE.  -4.18) .AND. (B .LT.  -4.04))) THEN
         CALL ENK1( 23,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 134.10) .AND. (L .LT. 135.04)) .AND.
     &   ((B .GE.  -4.10) .AND. (B .LT.  -0.70*L +  90.41))) THEN
         CALL ENK1( 23,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 133.03) .AND. (L .LT. 134.15)) .AND.
     &   ((B .GE.  -4.16) .AND. (B .LT.  -3.49))) THEN
         CALL ENK1( 23,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.46) .AND. (L .LT. 133.03)) .AND.
     &   ((B .GE.  -4.16) .AND. (B .LT.   0.46*L  -64.64))) THEN
         CALL ENK1( 23,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 133.73) .AND. (L .LT. 135.01)) .AND.
     &   ((B .GE.  -3.46) .AND. (B .LT.  -1.14))) THEN
         CALL ENK1( 24,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.00) .AND. (L .LT. 133.73)) .AND.
     &   ((B .GE.  -2.41) .AND. (B .LT.   0.21*L  -29.55))) THEN
         CALL ENK1( 24,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.00) .AND. (L .LT. 133.73)) .AND.
     &   ((B .GE.  -3.46) .AND. (B .LT.  -2.21))) THEN
         CALL ENK1( 24,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.00) .AND. (L .LT. 131.46)) .AND.
     &   ((B .GE.  -4.14) .AND. (B .LT.  -3.44))) THEN
         CALL ENK1( 24,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.46) .AND. (L .LT. 133.03)) .AND.
     &   ((B .LT.  -3.47) .AND. (B .GE.   0.46*L  -64.64))) THEN
         CALL ENK1( 24,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.04) .AND. (L .LT. 132.04)) .AND.
     &   ((B .GE.  -1.57) .AND. (B .LT.  -0.39*L +  50.20))) THEN
         CALL ENK1( 25,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.02) .AND. (L .LT. 132.04)) .AND.
     &   ((B .GE.  -1.98) .AND. (B .LT.  -1.57))) THEN
         CALL ENK1( 25,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.00) .AND. (L .LT. 132.04)) .AND.
     &   ((B .LT.  -1.98) .AND. (B .GE.   0.21*L  -29.55))) THEN
         CALL ENK1( 25,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 133.04) .AND. (L .LT. 133.73)) .AND.
     &   ((B .GE.  -1.16) .AND. (B .LT.  -0.96*L + 127.25))) THEN
         CALL ENK1( 26,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 132.04) .AND. (L .LT. 133.04)) .AND.
     &   ((B .GE.  -1.16) .AND. (B .LT.  -0.54))) THEN
         CALL ENK1( 26,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.92) .AND. (L .LT. 132.06)) .AND.
     &   ((B .GE.  -1.04) .AND. (B .LT.  -0.54))) THEN
         CALL ENK1( 26,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.90) .AND. (L .LT. 132.04)) .AND.
     &   ((B .LT.  -1.04) .AND. (B .GE.  -0.39*L +  50.20))) THEN
         CALL ENK1( 26,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 132.00) .AND. (L .LT. 133.73)) .AND.
     &   ((B .GE.  -1.64) .AND. (B .LT.  -1.16))) THEN
         CALL ENK1( 26,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 132.00) .AND. (L .LT. 133.73)) .AND.
     &   ((B .LT.  -1.64) .AND. (B .GE.   0.21*L  -29.55))) THEN
         CALL ENK1( 26,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.16) .AND. (L .LT. 134.98)) .AND.
     &   ((B .GE.   2.02) .AND. (B .LT.  -0.31*L +  44.00))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.00) .AND. (L .LT. 131.18)) .AND.
     &   ((B .GE.  -0.54) .AND. (B .LT.   3.27))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.04) .AND. (L .LT. 130.90)) .AND.
     &   ((B .GE.  -0.76) .AND. (B .LT.  -0.54))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.18) .AND. (L .LT. 135.00)) .AND.
     &   ((B .GE.  -0.54) .AND. (B .LT.   2.03))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 130.06) .AND. (L .LT. 130.90)) .AND.
     &   ((B .LT.  -0.80) .AND. (B .GE.  -0.39*L +  50.20))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 131.20) .AND. (L .LT. 135.00)) .AND.
     &   ((B .LT.   2.03) .AND. (B .GE.   0.68*L  -89.09))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 133.04) .AND. (L .LT. 133.73)) .AND.
     &   ((B .LT.  -0.52) .AND. (B .GE.  -0.96*L + 127.25))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 133.73) .AND. (L .LT. 134.97)) .AND.
     &   ((B .GE.  -1.14) .AND. (B .LT.  -0.53))) THEN
         CALL ENK1( 27,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 135.00) .AND. (L .LT. 138.49)) .AND.
     &   ((B .GE.   0.20) .AND. (B .LT.  -0.41*L +  56.97))) THEN
         CALL ENK1( 28,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 134.99) .AND. (L .LT. 138.46)) .AND.
     &   ((B .GE.  -3.13) .AND. (B .LT.   0.18))) THEN
         CALL ENK1( 28,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 135.03) .AND. (L .LT. 135.46)) .AND.
     &   ((B .GE.  -3.63) .AND. (B .LT.  -3.13))) THEN
         CALL ENK1( 28,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 134.98) .AND. (L .LT. 138.46)) .AND.
     &   ((B .GE.   1.64) .AND. (B .LT.   6.03))) THEN
         CALL ENK1( 29,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 135.01) .AND. (L .LT. 138.49)) .AND.
     &   ((B .LT.   1.64) .AND. (B .GE.  -0.41*L +  56.97))) THEN
         CALL ENK1( 29,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 138.44) .AND. (L .LT. 140.34)) .AND.
     &   ((B .GE.  -3.15) .AND. (B .LT.   0.06))) THEN
         CALL ENK1( 30,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 140.35) .AND. (L .LT. 144.97)) .AND.
     &   ((B .GE.  -5.00) .AND. (B .LT.   0.06))) THEN
         CALL ENK1( 31,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 140.37) .AND. (L .LT. 150.00)) .AND.
     &   ((B .GE.  -7.59) .AND. (B .LT.  -5.00))) THEN
         CALL ENK1( 31,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 139.95) .AND. (L .LT. 144.97)) .AND.
     &   ((B .GE.   0.06) .AND. (B .LT.   5.00))) THEN
         CALL ENK1( 32,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 148.01) .AND. (L .LT. 150.00)) .AND.
     &   ((B .GE.  -2.99) .AND. (B .LT.   0.01))) THEN
         CALL ENK1( 33,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 146.74) .AND. (L .LT. 148.01)) .AND.
     &   ((B .GE.  -1.62) .AND. (B .LT.   0.06))) THEN
         CALL ENK1( 33,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 144.92) .AND. (L .LT. 147.95)) .AND.
     &   ((B .GE.   0.06) .AND. (B .LT.   5.19))) THEN
         CALL ENK1( 34,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 147.97) .AND. (L .LT. 150.00)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   5.24))) THEN
         CALL ENK1( 35,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP259C(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 150.00) .AND. (L .LT. 158.00)) .AND.
     &   ((B .GE.  -6.10) .AND. (B .LT.   5.02))) THEN
         CALL ENK1( 35,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 150.00) .AND. (L .LT. 154.95)) .AND.
     &   ((B .GE.  -7.66) .AND. (B .LT.  -6.10))) THEN
         CALL ENK1( 35,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 150.00) .AND. (L .LT. 150.04)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   5.02))) THEN
         CALL ENK1( 35,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 158.00) .AND. (L .LT. 160.09)) .AND.
     &   ((B .GE.  -0.85) .AND. (B .LT.   0.04))) THEN
         CALL ENK1( 36,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 158.00) .AND. (L .LT. 164.45)) .AND.
     &   ((B .GE.  -2.65) .AND. (B .LT.  -0.85))) THEN
         CALL ENK1( 36,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.03) .AND. (L .LT. 165.67)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.   1.19))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.62) .AND. (L .LT. 166.91)) .AND.
     &   ((B .GE.  -0.08) .AND. (B .LT.  -0.99*L + 165.18))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.02) .AND. (L .LT. 166.92)) .AND.
     &   ((B .GE.  -5.53) .AND. (B .LT.  -0.08))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 164.45) .AND. (L .LT. 165.04)) .AND.
     &   ((B .GE.  -2.65) .AND. (B .LT.   0.06))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 160.04) .AND. (L .LT. 165.04)) .AND.
     &   ((B .GE.  -4.04) .AND. (B .LT.  -2.65))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 166.92) .AND. (L .LT. 170.41)) .AND.
     &   ((B .GE.  -1.07) .AND. (B .LT.  -0.42))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 166.90) .AND. (L .LT. 169.76)) .AND.
     &   ((B .GE.  -2.17) .AND. (B .LT.  -1.02))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 166.92) .AND. (L .LT. 169.78)) .AND.
     &   ((B .LT.  -2.19) .AND. (B .GE.   1.17*L  -200.85))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 168.69) .AND. (L .LT. 170.38)) .AND.
     &   ((B .GE.  -0.42) .AND. (B .LT.  -0.06))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 167.57) .AND. (L .LT. 168.67)) .AND.
     &   ((B .GE.  -0.44) .AND. (B .LT.   0.41*L  -69.05))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 166.93) .AND. (L .LT. 167.57)) .AND.
     &   ((B .GE.  -0.44) .AND. (B .LT.  -0.48*L +  80.05))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 168.67) .AND. (L .LT. 170.00)) .AND.
     &   ((B .GE.  -0.11) .AND. (B .LT.   2.00))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 169.98) .AND. (L .LT. 171.43)) .AND.
     &   ((B .GE.  -0.08) .AND. (B .LT.  -0.45*L +  77.16))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.38) .AND. (L .LT. 171.43)) .AND.
     &   ((B .GE.  -2.63) .AND. (B .LT.  -0.06))) THEN
         CALL ENK1( 37,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 161.09) .AND. (L .LT. 164.45)) .AND.
     &   ((B .GE.  -0.85) .AND. (B .LT.   0.08))) THEN
         CALL ENK1( 38,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 160.64) .AND. (L .LT. 161.83)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   1.40))) THEN
         CALL ENK1( 39,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 161.83) .AND. (L .LT. 163.16)) .AND.
     &   ((B .GE.   0.73) .AND. (B .LT.  -0.49*L +  80.78))) THEN
         CALL ENK1( 39,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 161.83) .AND. (L .LT. 163.16)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   0.73))) THEN
         CALL ENK1( 39,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 158.00) .AND. (L .LT. 165.03)) .AND.
     &   ((B .GE.   1.38) .AND. (B .LT.   7.61))) THEN
         CALL ENK1( 40,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 158.00) .AND. (L .LT. 160.62)) .AND.
     &   ((B .GE.   0.04) .AND. (B .LT.   1.38))) THEN
         CALL ENK1( 40,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 163.15) .AND. (L .LT. 165.03)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   1.38))) THEN
         CALL ENK1( 40,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 161.83) .AND. (L .LT. 163.16)) .AND.
     &   ((B .LT.   1.40) .AND. (B .GE.  -0.49*L +  80.78))) THEN
         CALL ENK1( 40,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 166.91) .AND. (L .LT. 168.69)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.   1.81))) THEN
         CALL ENK2( 41,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.62) .AND. (L .LT. 166.93)) .AND.
     &   ((B .GE.   1.19) .AND. (B .LT.   1.81))) THEN
         CALL ENK2( 41,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.62) .AND. (L .LT. 166.91)) .AND.
     &   ((B .LT.   1.19) .AND. (B .GE.  -0.99*L + 165.18))) THEN
         CALL ENK2( 41,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 166.93) .AND. (L .LT. 167.57)) .AND.
     &   ((B .LT.  -0.08) .AND. (B .GE.  -0.48*L +  80.05))) THEN
         CALL ENK2( 41,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 167.57) .AND. (L .LT. 168.67)) .AND.
     &   ((B .LT.  -0.08) .AND. (B .GE.   0.41*L  -69.05))) THEN
         CALL ENK2( 41,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.03) .AND. (L .LT. 168.70)) .AND.
     &   ((B .GE.   1.81) .AND. (B .LT.   4.50))) THEN
         CALL ENK2( 42,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 165.04) .AND. (L .LT. 168.73)) .AND.
     &   ((B .GE.   4.50) .AND. (B .LT.   7.61))) THEN
         CALL ENK2( 43,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 168.74) .AND. (L .LT. 170.04)) .AND.
     &   ((B .GE.   3.01) .AND. (B .LT.   7.61))) THEN
         CALL ENK2( 43,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.01) .AND. (L .LT. 172.01)) .AND.
     &   ((B .GE.   3.01) .AND. (B .LT.   3.78))) THEN
         CALL ENK2( 43,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.02) .AND. (L .LT. 172.04)) .AND.
     &   ((B .GE.   3.78) .AND. (B .LT.  -0.62*L + 110.68))) THEN
         CALL ENK2( 43,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 172.17) .AND. (L .LT. 174.59)) .AND.
     &   ((B .GE.  -4.04) .AND. (B .LT.  -1.83))) THEN
         CALL ENK2( 44,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 172.17) .AND. (L .LT. 173.24)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.  -1.40))) THEN
         CALL ENK2( 44,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 173.22) .AND. (L .LT. 173.50)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.  -1.60*L + 275.69))) THEN
         CALL ENK2( 44,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 172.02) .AND. (L .LT. 175.42)) .AND.
     &   ((B .GE.   2.00) .AND. (B .LT.   7.59))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.00) .AND. (L .LT. 172.02)) .AND.
     &   ((B .GE.   5.02) .AND. (B .LT.   7.57))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.00) .AND. (L .LT. 172.04)) .AND.
     &   ((B .LT.   5.02) .AND. (B .GE.  -0.62*L + 110.68))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.00) .AND. (L .LT. 172.01)) .AND.
     &   ((B .GE.   2.00) .AND. (B .LT.   3.01))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.00) .AND. (L .LT. 174.22)) .AND.
     &   ((B .GE.   0.56) .AND. (B .LT.   2.00))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 170.00) .AND. (L .LT. 171.43)) .AND.
     &   ((B .LT.   0.56) .AND. (B .GE.  -0.45*L +  77.16))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 171.43) .AND. (L .LT. 174.22)) .AND.
     &   ((B .GE.  -1.40) .AND. (B .LT.   0.56))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 174.23) .AND. (L .LT. 174.63)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.  -5.14*L + 897.75))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 174.22) .AND. (L .LT. 174.60)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.  -0.04))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 173.50) .AND. (L .LT. 174.22)) .AND.
     &   ((B .GE.  -1.86) .AND. (B .LT.  -1.40))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 173.24) .AND. (L .LT. 173.50)) .AND.
     &   ((B .GE.  -1.81) .AND. (B .LT.   0.18*L  -33.14))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 171.43) .AND. (L .LT. 172.21)) .AND.
     &   ((B .GE.  -2.60) .AND. (B .LT.  -1.40))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 173.22) .AND. (L .LT. 173.50)) .AND.
     &   ((B .LT.  -1.40) .AND. (B .GE.  -1.60*L + 275.69))) THEN
         CALL ENK2( 45,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 175.03) .AND. (L .LT. 178.02)) .AND.
     &   ((B .GE.  -4.23) .AND. (B .LT.  -2.00))) THEN
         CALL ENK2( 46,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 177.82) .AND. (L .LT. 179.47)) .AND.
     &   ((B .GE.  -0.61) .AND. (B .LT.  -0.32))) THEN
         CALL ENK2( 47,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 175.01) .AND. (L .LT. 179.50)) .AND.
     &   ((B .GE.  -2.03) .AND. (B .LT.  -0.61))) THEN
         CALL ENK2( 47,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 178.05) .AND. (L .LT. 179.50)) .AND.
     &   ((B .GE.  -3.01) .AND. (B .LT.  -2.03))) THEN
         CALL ENK2( 47,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 175.01) .AND. (L .LT. 176.17)) .AND.
     &   ((B .GE.  -0.64) .AND. (B .LT.  -0.41*L +  71.55))) THEN
         CALL ENK2( 47,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 179.50) .AND. (L .LT. 180.00)) .AND.
     &   ((B .GE.  -3.01) .AND. (B .LT.  -5.09*L + 913.00))) THEN
         CALL ENK2( 47,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 176.63) .AND. (L .LT. 178.66)) .AND.
     &   ((B .GE.   0.95) .AND. (B .LT.  -0.52*L +  93.85))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 176.18) .AND. (L .LT. 176.65)) .AND.
     &   ((B .GE.  -0.61) .AND. (B .LT.   2.00))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 174.63) .AND. (L .LT. 176.18)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.   2.00))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 174.23) .AND. (L .LT. 174.63)) .AND.
     &   ((B .LT.   1.98) .AND. (B .GE.  -5.14*L + 897.75))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 174.65) .AND. (L .LT. 176.17)) .AND.
     &   ((B .LT.  -0.06) .AND. (B .GE.  -0.41*L +  71.55))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 176.63) .AND. (L .LT. 177.82)) .AND.
     &   ((B .GE.  -0.64) .AND. (B .LT.   0.95))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 177.82) .AND. (L .LT. 180.00)) .AND.
     &   ((B .GE.  -0.04) .AND. (B .LT.   0.95))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 177.82) .AND. (L .LT. 179.50)) .AND.
     &   ((B .GE.  -0.32) .AND. (B .LT.  -0.04))) THEN
         CALL ENK2( 48,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 178.69) .AND. (L .LT. 179.04)) .AND.
     &   ((B .GE.   0.97) .AND. (B .LT.   5.00))) THEN
         CALL ENK2( 49,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 175.43) .AND. (L .LT. 178.66)) .AND.
     &   ((B .GE.   2.00) .AND. (B .LT.   5.00))) THEN
         CALL ENK2( 49,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 176.63) .AND. (L .LT. 178.66)) .AND.
     &   ((B .LT.   1.98) .AND. (B .GE.  -0.52*L +  93.85))) THEN
         CALL ENK2( 49,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 179.04) .AND. (L .LT. 180.00)) .AND.
     &   ((B .GE.   1.00) .AND. (B .LT.   5.00))) THEN
         CALL ENK2( 50,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 179.50) .AND. (L .LT. 180.00)) .AND.
     &   ((B .GE.  -0.64) .AND. (B .LT.  -0.04))) THEN
         CALL ENK2( 50,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 179.50) .AND. (L .LT. 180.00)) .AND.
     &   ((B .LT.  -0.64) .AND. (B .GE.  -5.09*L + 913.00))) THEN
         CALL ENK2( 50,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP260A(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 180.00) .AND. (L .LT. 183.07)) .AND.
     &   ((B .GE.  -2.65) .AND. (B .LT.   5.29))) THEN
         CALL ENK2( 50,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 180.00) .AND. (L .LT. 183.04)) .AND.
     &   ((B .GE.  -4.76) .AND. (B .LT.  -2.65))) THEN
         CALL ENK2( 50,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 183.04) .AND. (L .LT. 190.08)) .AND.
     &   ((B .GE.  -7.36) .AND. (B .LT.  -1.75))) THEN
         CALL ENK2( 51,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 182.99) .AND. (L .LT. 190.10)) .AND.
     &   ((B .GE.  -1.75) .AND. (B .LT.   1.25))) THEN
         CALL ENK2( 52,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 185.06) .AND. (L .LT. 187.08)) .AND.
     &   ((B .GE.   2.24) .AND. (B .LT.   3.88))) THEN
         CALL ENK2( 53,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 185.06) .AND. (L .LT. 189.08)) .AND.
     &   ((B .GE.   1.25) .AND. (B .LT.   2.24))) THEN
         CALL ENK2( 53,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 189.08) .AND. (L .LT. 190.03)) .AND.
     &   ((B .GE.   1.25) .AND. (B .LT.   1.42))) THEN
         CALL ENK2( 53,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 187.08) .AND. (L .LT. 188.25)) .AND.
     &   ((B .GE.   2.25) .AND. (B .LT.  -1.40*L + 265.83))) THEN
         CALL ENK2( 53,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 188.25) .AND. (L .LT. 189.70)) .AND.
     &   ((B .GE.   2.84) .AND. (B .LT.   3.67))) THEN
         CALL ENK2( 54,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 189.66) .AND. (L .LT. 192.84)) .AND.
     &   ((B .GE.   2.64) .AND. (B .LT.   5.22))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 188.25) .AND. (L .LT. 189.70)) .AND.
     &   ((B .GE.   3.70) .AND. (B .LT.   5.24))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 187.06) .AND. (L .LT. 188.27)) .AND.
     &   ((B .GE.   3.88) .AND. (B .LT.   5.25))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 187.08) .AND. (L .LT. 188.25)) .AND.
     &   ((B .LT.   3.87) .AND. (B .GE.  -1.40*L + 265.83))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 188.25) .AND. (L .LT. 189.78)) .AND.
     &   ((B .GE.   2.24) .AND. (B .LT.   2.85))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 189.08) .AND. (L .LT. 189.70)) .AND.
     &   ((B .GE.   1.25) .AND. (B .LT.   2.26))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 189.70) .AND. (L .LT. 190.05)) .AND.
     &   ((B .GE.   1.25) .AND. (B .LT.   2.66))) THEN
         CALL ENK2( 55,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.97) .AND. (L .LT. 191.53)) .AND.
     &   ((B .GE.  -0.91) .AND. (B .LT.   0.07))) THEN
         CALL ENK2( 56,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 191.54) .AND. (L .LT. 191.70)) .AND.
     &   ((B .GE.  -0.91) .AND. (B .LT.  -0.14))) THEN
         CALL ENK2( 56,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 191.74) .AND. (L .LT. 192.86)) .AND.
     &   ((B .GE.  -0.91) .AND. (B .LT.  -1.01*L + 193.71))) THEN
         CALL ENK2( 56,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.13) .AND. (L .LT. 192.55)) .AND.
     &   ((B .LT.   1.26) .AND. (B .GE.  -1.51*L + 293.80))) THEN
         CALL ENK2( 57,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.08) .AND. (L .LT. 192.86)) .AND.
     &   ((B .GE.   0.26) .AND. (B .LT.  -1.51*L + 291.60))) THEN
         CALL ENK2( 57,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.48) .AND. (L .LT. 192.76)) .AND.
     &   ((B .LT.   0.26) .AND. (B .GE.   1.23*L  -236.85))) THEN
         CALL ENK2( 57,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.05) .AND. (L .LT. 192.48)) .AND.
     &   ((B .LT.   0.26) .AND. (B .GE.  -0.72*L + 138.65))) THEN
         CALL ENK2( 57,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.03) .AND. (L .LT. 192.86)) .AND.
     &   ((B .LT.   2.62) .AND. (B .GE.   1.46*L  -279.00))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.04) .AND. (L .LT. 192.08)) .AND.
     &   ((B .GE.   1.44) .AND. (B .LT.   2.65))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.04) .AND. (L .LT. 192.08)) .AND.
     &   ((B .GE.   0.94) .AND. (B .LT.   1.43))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 191.70) .AND. (L .LT. 192.08)) .AND.
     &   ((B .GE.   0.93) .AND. (B .LT.  -1.51*L + 291.30))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.04) .AND. (L .LT. 192.08)) .AND.
     &   ((B .GE.   0.26) .AND. (B .LT.   0.93))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.99) .AND. (L .LT. 191.55)) .AND.
     &   ((B .GE.   0.05) .AND. (B .LT.   0.26))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 191.53) .AND. (L .LT. 192.05)) .AND.
     &   ((B .GE.  -0.12) .AND. (B .LT.   0.26))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.10) .AND. (L .LT. 190.98)) .AND.
     &   ((B .GE.  -2.80) .AND. (B .LT.   0.27))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.10) .AND. (L .LT. 195.81)) .AND.
     &   ((B .GE.  -7.41) .AND. (B .LT.  -2.80))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 190.97) .AND. (L .LT. 195.02)) .AND.
     &   ((B .GE.  -2.81) .AND. (B .LT.  -0.91))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.02) .AND. (L .LT. 195.83)) .AND.
     &   ((B .GE.  -2.80) .AND. (B .LT.  -2.70*L + 526.22))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.76) .AND. (L .LT. 195.02)) .AND.
     &   ((B .GE.  -0.91) .AND. (B .LT.  -0.46*L +  89.10))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.48) .AND. (L .LT. 192.76)) .AND.
     &   ((B .GE.  -0.15) .AND. (B .LT.   1.23*L  -236.85))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.05) .AND. (L .LT. 192.48)) .AND.
     &   ((B .GE.  -0.15) .AND. (B .LT.  -0.72*L + 138.65))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 191.74) .AND. (L .LT. 192.86)) .AND.
     &   ((B .LT.  -0.15) .AND. (B .GE.  -1.01*L + 193.71))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 191.53) .AND. (L .LT. 192.05)) .AND.
     &   ((B .GE.  -0.14) .AND. (B .LT.   0.26))) THEN
         CALL ENK2( 58,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.03) .AND. (L .LT. 192.46)) .AND.
     &   ((B .GE.   1.26) .AND. (B .LT.   1.46*L  -279.00))) THEN
         CALL ENK2( 59,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.03) .AND. (L .LT. 192.86)) .AND.
     &   ((B .LT.   1.26) .AND. (B .GE.  -1.51*L + 291.30))) THEN
         CALL ENK2( 59,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.46) .AND. (L .LT. 192.86)) .AND.
     &   ((B .GE.   1.26) .AND. (B .LT.  -1.51*L + 292.60))) THEN
         CALL ENK2( 59,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.76) .AND. (L .LT. 195.02)) .AND.
     &   ((B .LT.   0.21) .AND. (B .GE.  -0.46*L +  89.10))) THEN
         CALL ENK2( 59,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.86) .AND. (L .LT. 195.02)) .AND.
     &   ((B .GE.   0.21) .AND. (B .LT.  -0.46*L +  90.07))) THEN
         CALL ENK2( 59,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.02) .AND. (L .LT. 196.00)) .AND.
     &   ((B .GE.   0.19) .AND. (B .LT.   0.85*L  -164.90))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.02) .AND. (L .LT. 195.85)) .AND.
     &   ((B .GE.  -0.79) .AND. (B .LT.   0.19))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.02) .AND. (L .LT. 195.83)) .AND.
     &   ((B .LT.  -0.79) .AND. (B .GE.  -2.70*L + 526.22))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.83) .AND. (L .LT. 197.21)) .AND.
     &   ((B .GE.  -6.81) .AND. (B .LT.   0.19))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 197.24) .AND. (L .LT. 199.97)) .AND.
     &   ((B .GE.  -6.81) .AND. (B .LT.  -1.54*L + 301.20))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.83) .AND. (L .LT. 199.99)) .AND.
     &   ((B .GE.  -7.45) .AND. (B .LT.  -6.81))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.96) .AND. (L .LT. 200.94)) .AND.
     &   ((B .GE.  -7.48) .AND. (B .LT.  -0.71*L + 135.37))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 197.21) .AND. (L .LT. 199.99)) .AND.
     &   ((B .LT.  -1.83) .AND. (B .GE.   0.28*L  -58.02))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 197.22) .AND. (L .LT. 199.99)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.   0.16))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 196.00) .AND. (L .LT. 199.99)) .AND.
     &   ((B .GE.   0.19) .AND. (B .LT.  -0.23*L +  46.82))) THEN
         CALL ENK2( 60,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 196.03) .AND. (L .LT. 199.95)) .AND.
     &   ((B .GE.   1.04) .AND. (B .LT.   7.74))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.85) .AND. (L .LT. 196.02)) .AND.
     &   ((B .GE.   2.62) .AND. (B .LT.   7.75))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.86) .AND. (L .LT. 195.03)) .AND.
     &   ((B .LT.   1.27) .AND. (B .GE.  -1.51*L + 292.02))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.00) .AND. (L .LT. 196.03)) .AND.
     &   ((B .GE.   1.03) .AND. (B .LT.   2.62))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 195.02) .AND. (L .LT. 196.00)) .AND.
     &   ((B .LT.   1.05) .AND. (B .GE.   0.85*L  -164.90))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 196.00) .AND. (L .LT. 199.99)) .AND.
     &   ((B .LT.   1.01) .AND. (B .GE.  -0.23*L +  46.81))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.86) .AND. (L .LT. 195.00)) .AND.
     &   ((B .GE.   1.26) .AND. (B .LT.   2.60))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.46) .AND. (L .LT. 192.86)) .AND.
     &   ((B .GE.   2.02) .AND. (B .LT.   1.46*L  -279.00))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 192.46) .AND. (L .LT. 192.86)) .AND.
     &   ((B .LT.   2.02) .AND. (B .GE.  -1.51*L + 292.60))) THEN
         CALL ENK2( 61,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 201.98) .AND. (L .LT. 205.01)) .AND.
     &   ((B .GE.   0.90) .AND. (B .LT.   1.30))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.99) .AND. (L .LT. 202.01)) .AND.
     &   ((B .GE.   0.92) .AND. (B .LT.   0.22*L  -42.99))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.99) .AND. (L .LT. 205.01)) .AND.
     &   ((B .GE.  -6.79) .AND. (B .LT.   0.90))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.99) .AND. (L .LT. 200.94)) .AND.
     &   ((B .GE.  -6.49) .AND. (B .LT.  -0.71*L + 135.37))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 197.21) .AND. (L .LT. 199.99)) .AND.
     &   ((B .LT.  -2.63) .AND. (B .GE.  -1.54*L + 301.20))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 197.21) .AND. (L .LT. 199.99)) .AND.
     &   ((B .GE.  -2.63) .AND. (B .LT.   0.28*L  -58.02))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 200.93) .AND. (L .LT. 205.01)) .AND.
     &   ((B .GE.  -7.46) .AND. (B .LT.  -6.79))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 204.97) .AND. (L .LT. 207.01)) .AND.
     &   ((B .GE.  -7.49) .AND. (B .LT.  -2.42))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 205.01) .AND. (L .LT. 207.01)) .AND.
     &   ((B .GE.  -2.41) .AND. (B .LT.  -0.26*L +  51.40))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.92) .AND. (L .LT. 200.94)) .AND.
     &   ((B .LT.  -6.79) .AND. (B .GE.  -0.70*L + 133.16))) THEN
         CALL ENK2( 62,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.99) .AND. (L .LT. 202.01)) .AND.
     &   ((B .GE.   1.37) .AND. (B .LT.   2.34))) THEN
         CALL ENK2( 63,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.99) .AND. (L .LT. 202.01)) .AND.
     &   ((B .LT.   1.40) .AND. (B .GE.   0.22*L  -42.99))) THEN
         CALL ENK2( 63,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 207.10) .AND. (L .LT. 208.95)) .AND.
     &   ((B .GE.   0.10) .AND. (B .LT.   2.07))) THEN
         CALL ENK2( 64,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 205.01) .AND. (L .LT. 207.10)) .AND.
     &   ((B .GE.  -0.83) .AND. (B .LT.   2.08))) THEN
         CALL ENK2( 64,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 207.09) .AND. (L .LT. 207.81)) .AND.
     &   ((B .GE.  -0.83) .AND. (B .LT.  -0.48))) THEN
         CALL ENK2( 65,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 207.01) .AND. (L .LT. 207.79)) .AND.
     &   ((B .GE.  -2.40) .AND. (B .LT.  -0.83))) THEN
         CALL ENK2( 65,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 204.98) .AND. (L .LT. 207.01)) .AND.
     &   ((B .GE.  -1.84) .AND. (B .LT.  -0.85))) THEN
         CALL ENK2( 65,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 205.01) .AND. (L .LT. 207.01)) .AND.
     &   ((B .LT.  -1.84) .AND. (B .GE.  -0.26*L +  51.40))) THEN
         CALL ENK2( 65,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 207.79) .AND. (L .LT. 209.52)) .AND.
     &   ((B .GE.  -1.31) .AND. (B .LT.  -0.49*L + 100.97))) THEN
         CALL ENK2( 66,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 207.79) .AND. (L .LT. 209.57)) .AND.
     &   ((B .GE.  -2.91) .AND. (B .LT.  -1.32))) THEN
         CALL ENK2( 66,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 201.99) .AND. (L .LT. 208.91)) .AND.
     &   ((B .GE.   2.07) .AND. (B .LT.   7.68))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 199.97) .AND. (L .LT. 201.99)) .AND.
     &   ((B .GE.   2.34) .AND. (B .LT.   7.71))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 202.01) .AND. (L .LT. 205.01)) .AND.
     &   ((B .GE.   1.30) .AND. (B .LT.   2.10))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.91) .AND. (L .LT. 210.00)) .AND.
     &   ((B .GE.   3.08) .AND. (B .LT.  -0.87*L + 185.80))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.91) .AND. (L .LT. 210.90)) .AND.
     &   ((B .GE.   2.03) .AND. (B .LT.   3.08))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.93) .AND. (L .LT. 209.88)) .AND.
     &   ((B .GE.   0.09) .AND. (B .LT.   2.07))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 209.88) .AND. (L .LT. 210.90)) .AND.
     &   ((B .GE.  -0.25) .AND. (B .LT.   2.05))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 209.58) .AND. (L .LT. 210.93)) .AND.
     &   ((B .GE.  -7.53) .AND. (B .LT.  -1.33))) THEN
         CALL ENK2( 68,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.91) .AND. (L .LT. 210.86)) .AND.
     &   ((B .GE.   4.03) .AND. (B .LT.   7.69))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.91) .AND. (L .LT. 210.00)) .AND.
     &   ((B .LT.   4.02) .AND. (B .GE.  -0.87*L + 185.80))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 209.91) .AND. (L .LT. 210.93)) .AND.
     &   ((B .GE.   3.11) .AND. (B .LT.   4.01))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP260B(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 208.88) .AND. (L .LT. 209.84)) .AND.
     &  ((B .GE.  -0.43) .AND. (B .LT.   2.91))) THEN
        CALL ENK2( 67,D,AV,A0)
        RETURN
      ENDIF
      IF (((L .GE. 209.86) .AND. (L .LT. 210.73)) .AND.
     &   ((B .GE.   1.27) .AND. (B .LT.  -1.84*L + 388.91))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 209.85) .AND. (L .LT. 210.73)) .AND.
     &   ((B .GE.   0.29) .AND. (B .LT.   1.27))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 209.84) .AND. (L .LT. 212.01)) .AND.
     &   ((B .GE.  -0.41) .AND. (B .LT.   0.29))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 210.73) .AND. (L .LT. 212.02)) .AND.
     &   ((B .GE.   0.29) .AND. (B .LT.  -0.75*L + 159.31))) THEN
         CALL ENK2( 67,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.88) .AND. (L .LT. 219.94)) .AND.
     &   ((B .GE.  -7.76) .AND. (B .LT.  -0.41))) THEN
         CALL ENK2( 68,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 212.01) .AND. (L .LT. 219.94)) .AND.
     &   ((B .GE.  -0.48) .AND. (B .LT.  -0.07))) THEN
         CALL ENK2( 68,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 212.54) .AND. (L .LT. 213.83)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.  -0.32*L +  66.82))) THEN
         CALL ENK2( 68,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 212.02) .AND. (L .LT. 212.54)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.31))) THEN
         CALL ENK2( 68,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 213.83) .AND. (L .LT. 219.94)) .AND.
     &   ((B .GE.  -0.09) .AND. (B .LT.   7.60))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 208.89) .AND. (L .LT. 213.83)) .AND.
     &   ((B .GE.   2.93) .AND. (B .LT.   7.60))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 209.86) .AND. (L .LT. 210.73)) .AND.
     &   ((B .LT.   2.94) .AND. (B .GE.  -1.84*L + 388.91))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 210.73) .AND. (L .LT. 212.02)) .AND.
     &   ((B .LT.   1.25) .AND. (B .GE.  -0.75*L + 159.31))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 212.02) .AND. (L .LT. 212.54)) .AND.
     &   ((B .GE.   0.33) .AND. (B .LT.   2.94))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 210.73) .AND. (L .LT. 212.02)) .AND.
     &   ((B .GE.   1.25) .AND. (B .LT.   2.94))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 212.54) .AND. (L .LT. 213.83)) .AND.
     &   ((B .LT.   0.31) .AND. (B .GE.  -0.32*L +  66.82))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 212.54) .AND. (L .LT. 213.83)) .AND.
     &   ((B .GE.   0.31) .AND. (B .LT.   2.91))) THEN
         CALL ENK2( 69,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 221.95) .AND. (L .LT. 224.96)) .AND.
     &   ((B .GE.   0.12) .AND. (B .LT.   3.00))) THEN
         CALL ENK2( 70,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 220.97) .AND. (L .LT. 224.33)) .AND.
     &   ((B .GE.  -1.80) .AND. (B .LT.   0.12))) THEN
         CALL ENK2( 71,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 224.33) .AND. (L .LT. 224.96)) .AND.
     &   ((B .LT.   0.12) .AND. (B .GE.   2.86*L  -641.00))) THEN
         CALL ENK2( 71,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 222.46) .AND. (L .LT. 223.13)) .AND.
     &   ((B .LT.  -1.80) .AND. (B .GE.   0.86*L  -193.75))) THEN
         CALL ENK2( 71,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 220.94) .AND. (L .LT. 222.51)) .AND.
     &   ((B .GE.  -2.36) .AND. (B .LT.  -1.80))) THEN
         CALL ENK2( 71,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 223.13) .AND. (L .LT. 224.97)) .AND.
     &   ((B .GE.  -2.38) .AND. (B .LT.  -1.83))) THEN
         CALL ENK2( 72,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 222.46) .AND. (L .LT. 223.13)) .AND.
     &   ((B .GE.  -2.36) .AND. (B .LT.   0.86*L  -193.75))) THEN
         CALL ENK2( 72,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 220.93) .AND. (L .LT. 224.96)) .AND.
     &   ((B .GE.  -3.25) .AND. (B .LT.  -2.38))) THEN
         CALL ENK2( 72,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 219.94) .AND. (L .LT. 221.90)) .AND.
     &   ((B .GE.  -7.72) .AND. (B .LT.  -5.10))) THEN
         CALL ENK2( 73,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 224.96) .AND. (L .LT. 226.93)) .AND.
     &   ((B .GE.  -2.24) .AND. (B .LT.   0.38))) THEN
         CALL ENK2( 74,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 224.33) .AND. (L .LT. 224.96)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.   2.85*L  -641.00))) THEN
         CALL ENK2( 74,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 224.97) .AND. (L .LT. 228.32)) .AND.
     &   ((B .GE.  -5.48) .AND. (B .LT.  -2.24))) THEN
         CALL ENK2( 75,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 224.95) .AND. (L .LT. 230.52)) .AND.
     &   ((B .GE.  -7.74) .AND. (B .LT.  -5.43))) THEN
         CALL ENK2( 76,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 228.32) .AND. (L .LT. 233.96)) .AND.
     &   ((B .GE.  -5.43) .AND. (B .LT.  -1.20))) THEN
         CALL ENK2( 77,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 233.96) .AND. (L .LT. 234.73)) .AND.
     &   ((B .GE.  -5.43) .AND. (B .LT.  -4.09))) THEN
         CALL ENK2( 77,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 233.96) .AND. (L .LT. 234.73)) .AND.
     &   ((B .GE.  -4.09) .AND. (B .LT.  -3.72*L + 869.30))) THEN
         CALL ENK2( 77,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 230.50) .AND. (L .LT. 232.96)) .AND.
     &   ((B .GE.  -7.69) .AND. (B .LT.  -7.28))) THEN
         CALL ENK2( 78,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 230.50) .AND. (L .LT. 231.27)) .AND.
     &   ((B .GE.  -7.28) .AND. (B .LT.  -6.66))) THEN
         CALL ENK2( 78,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 231.27) .AND. (L .LT. 232.94)) .AND.
     &   ((B .GE.  -7.28) .AND. (B .LT.  -0.36*L +  76.67))) THEN
         CALL ENK2( 78,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 232.93) .AND. (L .LT. 233.94)) .AND.
     &   ((B .GE.  -7.24) .AND. (B .LT.  -5.41))) THEN
         CALL ENK2( 79,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 231.33) .AND. (L .LT. 232.92)) .AND.
     &   ((B .GE.  -6.66) .AND. (B .LT.  -5.41))) THEN
         CALL ENK2( 79,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 231.27) .AND. (L .LT. 232.94)) .AND.
     &   ((B .LT.  -6.66) .AND. (B .GE.  -0.36*L +  76.67))) THEN
         CALL ENK2( 79,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 228.69) .AND. (L .LT. 231.12)) .AND.
     &   ((B .GE.   2.57) .AND. (B .LT.   4.04))) THEN
         CALL ENK2( 80,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 228.73) .AND. (L .LT. 230.36)) .AND.
     &   ((B .GE.  -0.02) .AND. (B .LT.   2.57))) THEN
         CALL ENK2( 80,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 230.31) .AND. (L .LT. 231.19)) .AND.
     &   ((B .LT.   2.57) .AND. (B .GE.   2.96*L  -681.59))) THEN
         CALL ENK2( 80,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 231.19) .AND. (L .LT. 232.64)) .AND.
     &   ((B .GE.  -1.25) .AND. (B .LT.   2.60))) THEN
         CALL ENK3( 81,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 230.31) .AND. (L .LT. 231.19)) .AND.
     &   ((B .GE.  -0.02) .AND. (B .LT.   2.96*L  -681.59))) THEN
         CALL ENK3( 81,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 230.31) .AND. (L .LT. 231.16)) .AND.
     &   ((B .GE.  -1.23) .AND. (B .LT.  -0.02))) THEN
         CALL ENK3( 81,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.73) .AND. (L .LT. 235.75)) .AND.
     &   ((B .GE.  -1.23) .AND. (B .LT.   0.00))) THEN
         CALL ENK3( 82,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.32) .AND. (L .LT. 235.79)) .AND.
     &   ((B .GE.  -2.38) .AND. (B .LT.  -1.23))) THEN
         CALL ENK3( 83,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 233.96) .AND. (L .LT. 234.33)) .AND.
     &   ((B .LT.  -1.23) .AND. (B .GE.  -3.72*L + 869.30))) THEN
         CALL ENK3( 83,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.44) .AND. (L .LT. 235.76)) .AND.
     &   ((B .LT.  -2.40) .AND. (B .GE.   0.73*L  -174.42))) THEN
         CALL ENK3( 83,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.33) .AND. (L .LT. 234.73)) .AND.
     &   ((B .LT.  -2.36) .AND. (B .GE.  -3.72*L + 867.93))) THEN
         CALL ENK3( 83,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 236.95) .AND. (L .LT. 238.37)) .AND.
     &   ((B .GE.  -6.20) .AND. (B .LT.  -4.42))) THEN
         CALL ENK3( 84,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 238.36) .AND. (L .LT. 240.89)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -3.97))) THEN
         CALL ENK3( 85,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 236.93) .AND. (L .LT. 238.36)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -6.20))) THEN
         CALL ENK3( 85,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.96) .AND. (L .LT. 236.95)) .AND.
     &   ((B .GE.  -7.69) .AND. (B .LT.  -4.42))) THEN
         CALL ENK3( 85,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 236.93) .AND. (L .LT. 237.40)) .AND.
     &   ((B .LT.   3.00) .AND. (B .GE.   6.37*L  -1509.31))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.99) .AND. (L .LT. 236.93)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.   2.98))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 233.94) .AND. (L .LT. 234.97)) .AND.
     &   ((B .LT.   3.00) .AND. (B .GE.  -2.89*L + 678.97))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 236.94) .AND. (L .LT. 237.52)) .AND.
     &   ((B .GE.  -1.03) .AND. (B .LT.  -1.86*L + 440.62))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 235.72) .AND. (L .LT. 236.92)) .AND.
     &   ((B .GE.  -1.03) .AND. (B .LT.   0.00))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.99) .AND. (L .LT. 235.75)) .AND.
     &   ((B .GE.  -1.23) .AND. (B .LT.   0.00))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 235.70) .AND. (L .LT. 237.53)) .AND.
     &   ((B .GE.  -3.00) .AND. (B .LT.  -1.03))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.95) .AND. (L .LT. 235.76)) .AND.
     &   ((B .GE.  -3.03) .AND. (B .LT.   0.73*L  -174.42))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 234.92) .AND. (L .LT. 235.77)) .AND.
     &   ((B .GE.  -4.42) .AND. (B .LT.  -3.03))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 235.76) .AND. (L .LT. 237.54)) .AND.
     &   ((B .GE.  -4.42) .AND. (B .LT.  -3.03))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 237.55) .AND. (L .LT. 238.09)) .AND.
     &   ((B .GE.  -4.42) .AND. (B .LT.  -2.58*L + 610.00))) THEN
         CALL ENK3( 86,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 236.94) .AND. (L .LT. 237.40)) .AND.
     &   ((B .GE.  -0.02) .AND. (B .LT.   6.37*L  -1509.31))) THEN
         CALL ENK3( 87,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 237.40) .AND. (L .LT. 239.40)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.  -1.57*L + 375.80))) THEN
         CALL ENK3( 87,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 236.94) .AND. (L .LT. 237.52)) .AND.
     &   ((B .LT.  -0.02) .AND. (B .GE.  -1.86*L + 440.00))) THEN
         CALL ENK3( 87,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 237.52) .AND. (L .LT. 238.96)) .AND.
     &   ((B .GE.  -3.00) .AND. (B .LT.  -0.02))) THEN
         CALL ENK3( 87,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 238.96) .AND. (L .LT. 239.40)) .AND.
     &   ((B .LT.  -0.02) .AND. (B .GE.   7.31*L  -1750.20))) THEN
         CALL ENK3( 87,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 238.02) .AND. (L .LT. 240.96)) .AND.
     &   ((B .GE.   2.24) .AND. (B .LT.   0.41*L  -95.28))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 239.40) .AND. (L .LT. 240.90)) .AND.
     &   ((B .GE.  -2.96) .AND. (B .LT.   2.24))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 237.99) .AND. (L .LT. 239.40)) .AND.
     &   ((B .LT.   2.24) .AND. (B .GE.  -1.57*L + 375.80))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 238.96) .AND. (L .LT. 240.94)) .AND.
     &   ((B .GE.  -4.01) .AND. (B .LT.  -2.96))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 238.96) .AND. (L .LT. 239.40)) .AND.
     &   ((B .GE.  -2.96) .AND. (B .LT.   7.31*L  -1750.20))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP260C(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 240.00) .AND. (L .LT. 241.13)) .AND.
     &   ((B .GE.  -7.79) .AND. (B .LT.  -4.21))) THEN
         CALL ENK3( 85,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 241.13) .AND. (L .LT. 243.07)) .AND.
     &   ((B .GE.  -6.17) .AND. (B .LT.  -1.03*L + 244.16))) THEN
         CALL ENK3( 85,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 240.00) .AND. (L .LT. 245.18)) .AND.
     &   ((B .GE.  -7.79) .AND. (B .LT.  -6.17))) THEN
         CALL ENK3( 85,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.03) .AND. (L .LT. 243.05)) .AND.
     &   ((B .GE.   1.51) .AND. (B .LT.  -2.35*L + 572.49))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 241.99) .AND. (L .LT. 242.42)) .AND.
     &   ((B .GE.   0.57) .AND. (B .LT.  -0.96*L + 234.56))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.23) .AND. (L .LT. 243.05)) .AND.
     &   ((B .LT.   1.53) .AND. (B .GE.   0.79*L  -190.65))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.02) .AND. (L .LT. 242.23)) .AND.
     &   ((B .GE.   0.98) .AND. (B .LT.   1.51))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 240.42) .AND. (L .LT. 242.42)) .AND.
     &   ((B .LT.   0.58) .AND. (B .GE.   0.86*L  -208.15))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 240.42) .AND. (L .LT. 243.30)) .AND.
     &   ((B .GE.  -4.21) .AND. (B .LT.  -0.98*L + 234.25))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 239.03) .AND. (L .LT. 240.42)) .AND.
     &   ((B .GE.  -4.21) .AND. (B .LT.   2.86))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 240.42) .AND. (L .LT. 242.03)) .AND.
     &   ((B .GE.   0.58) .AND. (B .LT.   2.86))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 241.03) .AND. (L .LT. 242.03)) .AND.
     &   ((B .GE.   2.86) .AND. (B .LT.   1.16*L  -274.39))) THEN
         CALL ENK3( 88,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.42) .AND. (L .LT. 243.10)) .AND.
     &   ((B .GE.   0.21) .AND. (B .LT.  -0.59*L + 143.58))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.42) .AND. (L .LT. 243.10)) .AND.
     &   ((B .LT.   0.21) .AND. (B .GE.   0.54*L  -131.05))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 240.42) .AND. (L .LT. 242.42)) .AND.
     &   ((B .GE.  -1.13) .AND. (B .LT.   0.86*L  -208.15))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 240.45) .AND. (L .LT. 242.42)) .AND.
     &   ((B .LT.  -1.13) .AND. (B .GE.  -0.98*L + 234.25))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.42) .AND. (L .LT. 243.67)) .AND.
     &   ((B .GE.  -3.06) .AND. (B .LT.  -1.37*L + 330.81))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.42) .AND. (L .LT. 243.30)) .AND.
     &   ((B .LT.  -3.06) .AND. (B .GE.  -0.98*L + 234.25))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.29) .AND. (L .LT. 243.65)) .AND.
     &   ((B .LT.  -3.06) .AND. (B .GE.   3.06*L  -748.48))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.29) .AND. (L .LT. 243.65)) .AND.
     &   ((B .GE.  -3.06) .AND. (B .LT.  -1.37*L + 330.81))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.41) .AND. (L .LT. 243.45)) .AND.
     &   ((B .GE.  -1.33) .AND. (B .LT.  -1.21*L + 293.40))) THEN
         CALL ENK3( 89,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.65) .AND. (L .LT. 245.18)) .AND.
     &   ((B .GE.  -4.19) .AND. (B .LT.  -2.34))) THEN
         CALL ENK3( 90,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.15) .AND. (L .LT. 243.65)) .AND.
     &   ((B .LT.  -2.34) .AND. (B .GE.  -1.37*L + 330.81))) THEN
         CALL ENK3( 90,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.29) .AND. (L .LT. 243.65)) .AND.
     &   ((B .GE.  -4.19) .AND. (B .LT.   3.06*L  -748.48))) THEN
         CALL ENK3( 90,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 241.13) .AND. (L .LT. 243.07)) .AND.
     &   ((B .LT.  -4.19) .AND. (B .GE.  -1.03*L + 244.16))) THEN
         CALL ENK3( 90,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.09) .AND. (L .LT. 245.02)) .AND.
     &   ((B .GE.  -6.18) .AND. (B .LT.  -4.19))) THEN
         CALL ENK3( 90,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.02) .AND. (L .LT. 245.18)) .AND.
     &   ((B .GE.  -6.14) .AND. (B .LT.  -2.34))) THEN
         CALL ENK3( 90,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 244.29) .AND. (L .LT. 249.08)) .AND.
     &   ((B .GE.  -1.19) .AND. (B .LT.  -0.55*L + 135.90))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 244.01) .AND. (L .LT. 244.29)) .AND.
     &   ((B .GE.  -0.13) .AND. (B .LT.   5.77*L  -1407.90))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.49) .AND. (L .LT. 244.26)) .AND.
     &   ((B .GE.  -1.16) .AND. (B .LT.  -0.13))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.41) .AND. (L .LT. 243.45)) .AND.
     &   ((B .LT.  -0.10) .AND. (B .GE.  -1.21*L + 293.40))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.16) .AND. (L .LT. 243.44)) .AND.
     &   ((B .GE.  -2.34) .AND. (B .LT.  -1.33))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.42) .AND. (L .LT. 243.15)) .AND.
     &   ((B .LT.  -1.35) .AND. (B .GE.  -1.37*L + 330.81))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.44) .AND. (L .LT. 245.35)) .AND.
     &   ((B .GE.  -2.34) .AND. (B .LT.  -1.17))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.19) .AND. (L .LT. 249.06)) .AND.
     &   ((B .GE.  -3.02) .AND. (B .LT.  -1.20))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.19) .AND. (L .LT. 249.08)) .AND.
     &   ((B .LT.  -3.02) .AND. (B .GE.   0.28*L  -73.00))) THEN
         CALL ENK3( 91,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 244.01) .AND. (L .LT. 244.29)) .AND.
     &   ((B .LT.   1.51) .AND. (B .GE.   5.77*L  -1407.90))) THEN
         CALL ENK3( 92,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.05) .AND. (L .LT. 244.02)) .AND.
     &   ((B .GE.   0.81) .AND. (B .LT.   1.51))) THEN
         CALL ENK3( 92,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.23) .AND. (L .LT. 243.05)) .AND.
     &   ((B .GE.   0.81) .AND. (B .LT.   0.79*L  -190.65))) THEN
         CALL ENK3( 92,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.26) .AND. (L .LT. 244.01)) .AND.
     &   ((B .LT.   0.81) .AND. (B .GE.  -0.59*L + 143.58))) THEN
         CALL ENK3( 92,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.02) .AND. (L .LT. 250.09)) .AND.
     &   ((B .GE.  -7.75) .AND. (B .LT.  -6.11))) THEN
         CALL ENK3( 93,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.18) .AND. (L .LT. 249.08)) .AND.
     &   ((B .GE.  -4.18) .AND. (B .LT.   0.30*L  -76.59))) THEN
         CALL ENK3( 94,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.18) .AND. (L .LT. 249.08)) .AND.
     &   ((B .GE.  -6.17) .AND. (B .LT.  -4.18))) THEN
         CALL ENK3( 94,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 246.23) .AND. (L .LT. 249.10)) .AND.
     &   ((B .GE.   0.85) .AND. (B .LT.   1.95))) THEN
         CALL ENK3( 95,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.03) .AND. (L .LT. 246.23)) .AND.
     &   ((B .GE.   1.38) .AND. (B .LT.   0.48*L  -116.20))) THEN
         CALL ENK3( 95,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.05) .AND. (L .LT. 246.21)) .AND.
     &   ((B .GE.   0.88) .AND. (B .LT.   1.38))) THEN
         CALL ENK3( 95,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.07) .AND. (L .LT. 249.08)) .AND.
     &   ((B .LT.   0.87) .AND. (B .GE.  -0.55*L + 135.90))) THEN
         CALL ENK3( 95,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.05) .AND. (L .LT. 250.09)) .AND.
     &   ((B .GE.   1.94) .AND. (B .LT.   7.55))) THEN
         CALL ENK3( 96,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.09) .AND. (L .LT. 243.09)) .AND.
     &   ((B .GE.   3.94) .AND. (B .LT.   7.52))) THEN
         CALL ENK3( 96,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 242.03) .AND. (L .LT. 243.05)) .AND.
     &   ((B .LT.   3.94) .AND. (B .GE.  -2.35*L + 572.49))) THEN
         CALL ENK3( 96,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 243.03) .AND. (L .LT. 245.03)) .AND.
     &   ((B .GE.   1.48) .AND. (B .LT.   1.94))) THEN
         CALL ENK3( 96,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 245.03) .AND. (L .LT. 246.23)) .AND.
     &   ((B .LT.   1.94) .AND. (B .GE.   0.48*L  -116.20))) THEN
         CALL ENK3( 96,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 250.09) .AND. (L .LT. 257.48)) .AND.
     &   ((B .GE.  -7.76) .AND. (B .LT.  -5.09))) THEN
         CALL ENK3( 97,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 250.07) .AND. (L .LT. 255.10)) .AND.
     &   ((B .GE.  -5.10) .AND. (B .LT.   2.97))) THEN
         CALL ENK3( 98,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 255.10) .AND. (L .LT. 257.46)) .AND.
     &   ((B .GE.  -5.11) .AND. (B .LT.  -0.02))) THEN
         CALL ENK3( 99,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 255.10) .AND. (L .LT. 260.01)) .AND.
     &   ((B .GE.  -0.08) .AND. (B .LT.   4.98))) THEN
         CALL ENK3(100,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 257.48) .AND. (L .LT. 259.96)) .AND.
     &   ((B .GE.  -5.09) .AND. (B .LT.  -0.10))) THEN
         CALL ENK3(100,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 259.96) .AND. (L .LT. 260.60)) .AND.
     &   ((B .GE.  -2.81) .AND. (B .LT.  -4.88*L + 1268.65))) THEN
         CALL ENK3(100,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 259.96) .AND. (L .LT. 260.60)) .AND.
     &   ((B .LT.  -2.81) .AND. (B .GE.   3.39*L  -886.52))) THEN
         CALL ENK3(100,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.63) .AND. (L .LT. 264.96)) .AND.
     &   ((B .GE.  -3.47) .AND. (B .LT.  -1.54*L + 404.60))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.04) .AND. (L .LT. 264.63)) .AND.
     &   ((B .GE.  -3.47) .AND. (B .LT.  -2.92))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.43) .AND. (L .LT. 264.04)) .AND.
     &   ((B .GE.  -3.47) .AND. (B .LT.   0.88*L  -235.25))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.43) .AND. (L .LT. 264.21)) .AND.
     &   ((B .LT.  -3.47) .AND. (B .GE.  -0.74*L + 191.44))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.19) .AND. (L .LT. 265.04)) .AND.
     &   ((B .GE.  -5.56) .AND. (B .LT.  -3.47))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.01) .AND. (L .LT. 264.19)) .AND.
     &   ((B .GE.  -5.56) .AND. (B .LT.   1.26*L  -336.97))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.04) .AND. (L .LT. 268.02)) .AND.
     &   ((B .GE.  -5.04) .AND. (B .LT.  -0.55*L + 142.35))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.94) .AND. (L .LT. 267.99)) .AND.
     &   ((B .GE.  -5.56) .AND. (B .LT.  -5.04))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.00) .AND. (L .LT. 267.97)) .AND.
     &   ((B .GE.  -7.69) .AND. (B .LT.  -5.57))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 268.00) .AND. (L .LT. 269.23)) .AND.
     &   ((B .GE.  -7.69) .AND. (B .LT.  -2.10*L + 557.80))) THEN
         CALL ENK3(101,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.60) .AND. (L .LT. 262.63)) .AND.
     &   ((B .GE.  -3.25) .AND. (B .LT.  -2.79))) THEN
         CALL ENK3(102,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.80) .AND. (L .LT. 262.04)) .AND.
     &   ((B .GE.   0.43) .AND. (B .LT.   1.30*L  -338.46))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 259.96) .AND. (L .LT. 260.60)) .AND.
     &   ((B .GE.   0.43) .AND. (B .LT.  -1.50*L + 390.95))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 259.96) .AND. (L .LT. 260.60)) .AND.
     &   ((B .LT.   0.43) .AND. (B .GE.  -4.88*L + 1268.65))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.60) .AND. (L .LT. 262.01)) .AND.
     &   ((B .GE.  -2.77) .AND. (B .LT.   0.42))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.03) .AND. (L .LT. 262.32)) .AND.
     &   ((B .GE.  -2.75) .AND. (B .LT.   0.06))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.35) .AND. (L .LT. 262.79)) .AND.
     &   ((B .GE.  -2.78) .AND. (B .LT.  -7.04*L + 1847.10))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.63) .AND. (L .LT. 262.81)) .AND.
     &   ((B .GE.  -4.01) .AND. (B .LT.  -2.78))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.52) .AND. (L .LT. 262.65)) .AND.
     &   ((B .GE.  -4.04) .AND. (B .LT.  -3.19))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.82) .AND. (L .LT. 264.21)) .AND.
     &   ((B .GE.  -4.04) .AND. (B .LT.  -0.74*L + 191.44))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.01) .AND. (L .LT. 264.19)) .AND.
     &   ((B .LT.  -4.04) .AND. (B .GE.   1.26*L  -336.97))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.52) .AND. (L .LT. 263.01)) .AND.
     &   ((B .GE.  -5.57) .AND. (B .LT.  -4.04))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 259.96) .AND. (L .LT. 260.60)) .AND.
     &   ((B .GE.  -5.00) .AND. (B .LT.   3.39*L  -886.52))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.01) .AND. (L .LT. 260.51)) .AND.
     &   ((B .GE.  -5.57) .AND. (B .LT.  -5.00))) THEN
         CALL ENK3(103,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.01) .AND. (L .LT. 261.55)) .AND.
     &   ((B .GE.   2.00) .AND. (B .LT.   3.85))) THEN
         CALL ENK3(104,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 261.55) .AND. (L .LT. 262.04)) .AND.
     &   ((B .GE.   1.99) .AND. (B .LT.  -3.96*L + 1039.80))) THEN
         CALL ENK3(104,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.01) .AND. (L .LT. 260.43)) .AND.
     &   ((B .GE.   0.99) .AND. (B .LT.   1.99))) THEN
         CALL ENK3(104,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.01) .AND. (L .LT. 260.43)) .AND.
     &   ((B .LT.   0.98) .AND. (B .GE.  -1.50*L + 390.95))) THEN
         CALL ENK3(104,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.44) .AND. (L .LT. 260.80)) .AND.
     &   ((B .GE.   0.36) .AND. (B .LT.   1.99))) THEN
         CALL ENK3(104,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 260.80) .AND. (L .LT. 262.04)) .AND.
     &   ((B .LT.   2.02) .AND. (B .GE.   1.30*L  -338.46))) THEN
         CALL ENK3(104,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.47) .AND. (L .LT. 265.03)) .AND.
     &   ((B .LT.  -2.61) .AND. (B .GE.  -1.54*L + 404.60))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.04) .AND. (L .LT. 266.75)) .AND.
     &   ((B .GE.  -3.47) .AND. (B .LT.  -0.52*L + 135.21))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.04) .AND. (L .LT. 266.75)) .AND.
     &   ((B .LT.  -3.47) .AND. (B .GE.  -0.55*L + 142.35))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.75) .AND. (L .LT. 268.03)) .AND.
     &   ((B .GE.  -4.20) .AND. (B .LT.  -0.52*L + 135.21))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.74) .AND. (L .LT. 268.03)) .AND.
     &   ((B .GE.  -4.37) .AND. (B .LT.  -4.20))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.74) .AND. (L .LT. 268.02)) .AND.
     &   ((B .LT.  -4.37) .AND. (B .GE.  -0.55*L + 142.35))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 268.03) .AND. (L .LT. 268.97)) .AND.
     &   ((B .GE.  -5.04) .AND. (B .LT.  -0.90*L + 237.00))) THEN
         CALL ENK3(105,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.48) .AND. (L .LT. 262.79)) .AND.
     &   ((B .LT.  -0.82) .AND. (B .GE.  -7.04*L + 1847.10))) THEN
         CALL ENK3(106,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.79) .AND. (L .LT. 265.07)) .AND.
     &   ((B .GE.  -2.00) .AND. (B .LT.  -0.54*L + 141.17))) THEN
         CALL ENK3(106,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.44) .AND. (L .LT. 265.00)) .AND.
     &   ((B .LT.  -2.00) .AND. (B .GE.   0.88*L  -235.30))) THEN
         CALL ENK3(106,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.79) .AND. (L .LT. 263.44)) .AND.
     &   ((B .GE.  -3.01) .AND. (B .LT.  -2.00))) THEN
         CALL ENK3(106,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.82) .AND. (L .LT. 263.44)) .AND.
     &   ((B .LT.  -3.01) .AND. (B .GE.  -0.74*L + 191.44))) THEN
         CALL ENK3(106,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.25) .AND. (L .LT. 263.49)) .AND.
     &   ((B .GE.   1.97) .AND. (B .LT.   2.40))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.01) .AND. (L .LT. 262.28)) .AND.
     &   ((B .GE.   1.97) .AND. (B .LT.   1.44*L  -375.20))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.04) .AND. (L .LT. 263.47)) .AND.
     &   ((B .GE.   0.02) .AND. (B .LT.   1.96))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.35) .AND. (L .LT. 262.48)) .AND.
     &   ((B .LT.   0.02) .AND. (B .GE.  -7.04*L + 1847.10))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.47) .AND. (L .LT. 263.47)) .AND.
     &   ((B .GE.  -0.82) .AND. (B .LT.   0.02))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.47) .AND. (L .LT. 265.06)) .AND.
     &   ((B .GE.  -0.83) .AND. (B .LT.  -0.04))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.49) .AND. (L .LT. 264.08)) .AND.
     &   ((B .GE.   0.95) .AND. (B .LT.  -2.48*L + 655.95))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.48) .AND. (L .LT. 264.07)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.   0.95))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.08) .AND. (L .LT. 266.28)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.  -0.46*L + 122.40))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.04) .AND. (L .LT. 266.33)) .AND.
     &   ((B .GE.  -1.19) .AND. (B .LT.  -0.06))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.33) .AND. (L .LT. 266.53)) .AND.
     &   ((B .GE.  -1.19) .AND. (B .LT.  -3.97*L + 1056.98))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.07) .AND. (L .LT. 266.56)) .AND.
     &   ((B .LT.  -1.19) .AND. (B .GE.   0.60*L  -161.00))) THEN
         CALL ENK3(107,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.54) .AND. (L .LT. 265.07)) .AND.
     &   ((B .GE.   2.40) .AND. (B .LT.   3.85))) THEN
         CALL ENK3(108,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 262.48) .AND. (L .LT. 263.52)) .AND.
     &   ((B .GE.   2.37) .AND. (B .LT.   3.85))) THEN
         CALL ENK3(108,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 261.55) .AND. (L .LT. 262.47)) .AND.
     &   ((B .LT.   3.85) .AND. (B .GE.  -1.29*L + 341.33))) THEN
         CALL ENK3(108,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 263.49) .AND. (L .LT. 264.08)) .AND.
     &   ((B .LT.   2.39) .AND. (B .GE.  -2.48*L + 655.95))) THEN
         CALL ENK3(108,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.08) .AND. (L .LT. 265.05)) .AND.
     &   ((B .LT.   0.95) .AND. (B .GE.  -0.46*L + 122.40))) THEN
         CALL ENK3(108,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 264.08) .AND. (L .LT. 265.08)) .AND.
     &   ((B .GE.   0.95) .AND. (B .LT.   2.39))) THEN
         CALL ENK3(108,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.48) .AND. (L .LT. 266.56)) .AND.
     &   ((B .GE.  -1.84) .AND. (B .LT.   0.60*L  -161.00))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.05) .AND. (L .LT. 265.48)) .AND.
     &   ((B .GE.  -2.06) .AND. (B .LT.   0.60*L  -161.00))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.48) .AND. (L .LT. 266.26)) .AND.
     &   ((B .GE.  -2.06) .AND. (B .LT.  -1.84))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.05) .AND. (L .LT. 266.21)) .AND.
     &   ((B .GE.  -2.61) .AND. (B .LT.  -2.06))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.24) .AND. (L .LT. 266.71)) .AND.
     &   ((B .GE.  -2.61) .AND. (B .LT.  -1.23*L + 325.44))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.04) .AND. (L .LT. 266.71)) .AND.
     &   ((B .LT.  -2.61) .AND. (B .GE.  -0.52*L + 135.21))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.71) .AND. (L .LT. 267.32)) .AND.
     &   ((B .GE.  -3.48) .AND. (B .LT.  -1.23*L + 325.44))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.71) .AND. (L .LT. 267.32)) .AND.
     &   ((B .LT.  -3.48) .AND. (B .GE.  -0.52*L + 135.21))) THEN
         CALL ENK3(109,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.34) .AND. (L .LT. 270.00)) .AND.
     &   ((B .GE.  -0.38) .AND. (B .LT.   5.62))) THEN
         CALL ENK3(110,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.10) .AND. (L .LT. 266.34)) .AND.
     &   ((B .GE.   0.59) .AND. (B .LT.   5.59))) THEN
         CALL ENK3(110,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 265.08) .AND. (L .LT. 266.34)) .AND.
     &   ((B .LT.   0.59) .AND. (B .GE.  -0.46*L + 122.40))) THEN
         CALL ENK3(110,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 268.71) .AND. (L .LT. 270.96)) .AND.
     &   ((B .GE.  -1.61) .AND. (B .LT.  -0.36))) THEN
         CALL ENK3(111,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.53) .AND. (L .LT. 268.71)) .AND.
     &   ((B .GE.  -1.20) .AND. (B .LT.  -0.38))) THEN
         CALL ENK3(111,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.33) .AND. (L .LT. 266.53)) .AND.
     &   ((B .LT.  -0.38) .AND. (B .GE.  -3.97*L + 1056.98))) THEN
         CALL ENK3(111,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.53) .AND. (L .LT. 268.71)) .AND.
     &   ((B .LT.  -1.20) .AND. (B .GE.  -0.20*L +  52.25))) THEN
         CALL ENK3(111,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.27) .AND. (L .LT. 266.62)) .AND.
     &   ((B .GE.  -2.03) .AND. (B .LT.  -1.79))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.62) .AND. (L .LT. 267.12)) .AND.
     &   ((B .GE.  -2.04) .AND. (B .LT.  -0.45*L + 118.20))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.24) .AND. (L .LT. 267.11)) .AND.
     &   ((B .LT.  -2.04) .AND. (B .GE.  -1.23*L + 325.44))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 267.12) .AND. (L .LT. 269.28)) .AND.
     &   ((B .GE.  -3.00) .AND. (B .LT.  -0.45*L + 118.20))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 267.11) .AND. (L .LT. 267.32)) .AND.
     &   ((B .LT.  -3.00) .AND. (B .GE.  -1.23*L + 325.44))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 267.33) .AND. (L .LT. 268.03)) .AND.
     &   ((B .GE.  -3.77) .AND. (B .LT.  -3.02))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 267.32) .AND. (L .LT. 268.03)) .AND.
     &   ((B .LT.  -3.77) .AND. (B .GE.  -0.52*L + 135.21))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 268.04) .AND. (L .LT. 269.98)) .AND.
     &   ((B .GE.  -4.18) .AND. (B .LT.  -3.00))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 268.03) .AND. (L .LT. 269.97)) .AND.
     &   ((B .LT.  -4.18) .AND. (B .GE.  -0.90*L + 237.00))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 270.97)) .AND.
     &   ((B .GE.  -5.89) .AND. (B .LT.  -2.96))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.53) .AND. (L .LT. 268.71)) .AND.
     &   ((B .GE.  -1.63) .AND. (B .LT.  -0.20*L +  52.25))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.53) .AND. (L .LT. 269.28)) .AND.
     &   ((B .GE.  -1.79) .AND. (B .LT.  -1.59))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 266.62) .AND. (L .LT. 269.28)) .AND.
     &   ((B .LT.  -1.79) .AND. (B .GE.  -0.45*L + 118.20))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 269.28) .AND. (L .LT. 270.95)) .AND.
     &   ((B .GE.  -3.00) .AND. (B .LT.  -1.59))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP261A(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 270.00) .AND. (L .LT. 273.07)) .AND.
     &   ((B .GE.  -0.20) .AND. (B .LT.   0.77))) THEN
         CALL ENK3(111,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 273.07)) .AND.
     &   ((B .LT.  -0.20) .AND. (B .GE.   0.52*L  -142.21))) THEN
         CALL ENK3(111,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 272.05)) .AND.
     &   ((B .GE.  -3.74) .AND. (B .LT.  -0.28*L +  72.45))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 272.04)) .AND.
     &   ((B .GE.  -6.02) .AND. (B .LT.  -3.74))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 272.00)) .AND.
     &   ((B .LT.  -6.02) .AND. (B .GE.  -0.83*L + 218.00))) THEN
         CALL ENK3(112,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 275.07)) .AND.
     &   ((B .GE.   0.80) .AND. (B .LT.   3.79))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 273.07) .AND. (L .LT. 275.03)) .AND.
     &   ((B .GE.  -0.80) .AND. (B .LT.   0.77))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 273.07)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.   0.52*L  -142.21))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 273.04) .AND. (L .LT. 273.46)) .AND.
     &   ((B .GE.  -2.58) .AND. (B .LT.  -0.80))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 273.47) .AND. (L .LT. 274.44)) .AND.
     &   ((B .GE.  -2.55) .AND. (B .LT.  -1.80*L + 491.40))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 273.05)) .AND.
     &   ((B .GE.  -2.55) .AND. (B .LT.  -1.77))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 272.30)) .AND.
     &   ((B .GE.  -3.17) .AND. (B .LT.  -2.55))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 272.30) .AND. (L .LT. 273.04)) .AND.
     &   ((B .GE.  -3.17) .AND. (B .LT.  -0.83*L + 223.40))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 270.00) .AND. (L .LT. 273.04)) .AND.
     &   ((B .LT.  -3.17) .AND. (B .GE.  -0.28*L +  72.45))) THEN
         CALL ENK3(113,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 273.04) .AND. (L .LT. 274.44)) .AND.
     &   ((B .GE.  -3.34) .AND. (B .LT.  -2.58))) THEN
         CALL ENK3(114,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 272.30) .AND. (L .LT. 273.04)) .AND.
     &   ((B .LT.  -2.55) .AND. (B .GE.  -0.83*L + 223.40))) THEN
         CALL ENK3(114,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 276.06) .AND. (L .LT. 277.06)) .AND.
     &   ((B .GE.  -2.17) .AND. (B .LT.  -0.97*L + 266.80))) THEN
         CALL ENK3(115,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.00) .AND. (L .LT. 276.06)) .AND.
     &   ((B .GE.  -1.44) .AND. (B .LT.   0.26*L  -72.86))) THEN
         CALL ENK3(115,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.03) .AND. (L .LT. 276.06)) .AND.
     &   ((B .GE.  -2.20) .AND. (B .LT.  -1.44))) THEN
         CALL ENK3(115,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 274.44) .AND. (L .LT. 275.03)) .AND.
     &   ((B .GE.  -3.31) .AND. (B .LT.  -1.41))) THEN
         CALL ENK3(115,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 274.42) .AND. (L .LT. 275.01)) .AND.
     &   ((B .LT.  -3.31) .AND. (B .GE.  -1.39*L + 379.49))) THEN
         CALL ENK3(115,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.02) .AND. (L .LT. 277.06)) .AND.
     &   ((B .GE.  -4.14) .AND. (B .LT.  -2.17))) THEN
         CALL ENK3(115,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.94) .AND. (L .LT. 286.00)) .AND.
     &   ((B .GE.  -4.69) .AND. (B .LT.  -1.54*L + 435.65))) THEN
         CALL ENK3(116,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.94) .AND. (L .LT. 286.17)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -4.69))) THEN
         CALL ENK3(116,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 277.04) .AND. (L .LT. 284.94)) .AND.
     &   ((B .GE.  -7.71) .AND. (B .LT.  -3.27))) THEN
         CALL ENK3(116,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.01) .AND. (L .LT. 277.04)) .AND.
     &   ((B .GE.  -7.73) .AND. (B .LT.  -4.12))) THEN
         CALL ENK3(116,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.03) .AND. (L .LT. 277.06)) .AND.
     &   ((B .GE.  -0.77) .AND. (B .LT.   0.32*L  -88.72))) THEN
         CALL ENK3(117,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.03) .AND. (L .LT. 277.05)) .AND.
     &   ((B .GE.  -1.15) .AND. (B .LT.  -0.77))) THEN
         CALL ENK3(117,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.00) .AND. (L .LT. 276.06)) .AND.
     &   ((B .LT.  -1.15) .AND. (B .GE.   0.26*L  -72.86))) THEN
         CALL ENK3(117,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 276.06) .AND. (L .LT. 277.06)) .AND.
     &   ((B .LT.  -1.15) .AND. (B .GE.  -0.97*L + 266.80))) THEN
         CALL ENK3(117,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 277.06) .AND. (L .LT. 278.49)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.  -0.88*L + 243.78))) THEN
         CALL ENK3(117,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 277.06) .AND. (L .LT. 278.49)) .AND.
     &   ((B .GE.  -3.27) .AND. (B .LT.  -1.39))) THEN
         CALL ENK3(117,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.07) .AND. (L .LT. 276.61)) .AND.
     &   ((B .GE.   0.94) .AND. (B .LT.   4.98))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 276.63) .AND. (L .LT. 280.00)) .AND.
     &   ((B .GE.   0.91) .AND. (B .LT.  -1.20*L + 336.95))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 277.45) .AND. (L .LT. 279.97)) .AND.
     &   ((B .GE.  -0.13) .AND. (B .LT.   0.91))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 276.61) .AND. (L .LT. 277.42)) .AND.
     &   ((B .LT.   0.91) .AND. (B .GE.  -0.58*L + 161.30))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 277.06) .AND. (L .LT. 277.42)) .AND.
     &   ((B .GE.  -0.13) .AND. (B .LT.   1.50*L  -415.73))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 280.47)) .AND.
     &   ((B .GE.  -0.13) .AND. (B .LT.  -2.09*L + 586.22))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 278.50) .AND. (L .LT. 280.51)) .AND.
     &   ((B .GE.  -1.37) .AND. (B .LT.  -0.13))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 277.06) .AND. (L .LT. 278.49)) .AND.
     &   ((B .LT.  -0.13) .AND. (B .GE.  -0.88*L + 243.78))) THEN
         CALL ENK3(118,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 279.06) .AND. (L .LT. 280.00)) .AND.
     &   ((B .GE.   3.98) .AND. (B .LT.  -0.33*L +  96.45))) THEN
         CALL ENK3(119,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 278.07) .AND. (L .LT. 279.06)) .AND.
     &   ((B .GE.   3.27) .AND. (B .LT.   1.03*L  -283.05))) THEN
         CALL ENK3(119,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 279.06) .AND. (L .LT. 280.01)) .AND.
     &   ((B .GE.   3.27) .AND. (B .LT.   3.98))) THEN
         CALL ENK3(119,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 278.02) .AND. (L .LT. 280.00)) .AND.
     &   ((B .LT.   3.27) .AND. (B .GE.  -1.20*L + 336.95))) THEN
         CALL ENK3(119,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 275.03) .AND. (L .LT. 280.00)) .AND.
     &   ((B .GE.   4.93) .AND. (B .LT.   6.76))) THEN
         CALL ENK3(120,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 276.63) .AND. (L .LT. 278.07)) .AND.
     &   ((B .LT.   4.98) .AND. (B .GE.  -1.20*L + 336.95))) THEN
         CALL ENK3(120,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 278.08) .AND. (L .LT. 279.06)) .AND.
     &   ((B .GE.   4.29) .AND. (B .LT.   4.98))) THEN
         CALL ENK3(120,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 278.07) .AND. (L .LT. 279.06)) .AND.
     &   ((B .LT.   4.29) .AND. (B .GE.   1.03*L  -283.05))) THEN
         CALL ENK3(120,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 279.06) .AND. (L .LT. 280.00)) .AND.
     &   ((B .GE.   4.29) .AND. (B .LT.   4.98))) THEN
         CALL ENK3(120,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 279.06) .AND. (L .LT. 280.00)) .AND.
     &   ((B .LT.   4.29) .AND. (B .GE.  -0.33*L +  96.45))) THEN
         CALL ENK3(120,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 282.52)) .AND.
     &   ((B .GE.   3.03) .AND. (B .LT.  -0.87*L + 249.02))) THEN
         CALL ENK4(121,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 282.51)) .AND.
     &   ((B .GE.   1.94) .AND. (B .LT.   3.03))) THEN
         CALL ENK4(121,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 280.48)) .AND.
     &   ((B .GE.   0.91) .AND. (B .LT.   1.96))) THEN
         CALL ENK4(121,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 280.47)) .AND.
     &   ((B .LT.   0.91) .AND. (B .GE.  -2.09*L + 585.22))) THEN
         CALL ENK4(121,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.47) .AND. (L .LT. 282.51)) .AND.
     &   ((B .LT.   1.96) .AND. (B .GE.   1.03*L  -289.26))) THEN
         CALL ENK4(121,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 282.52) .AND. (L .LT. 285.01)) .AND.
     &   ((B .GE.   2.43) .AND. (B .LT.   6.83))) THEN
         CALL ENK4(122,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 282.52)) .AND.
     &   ((B .GE.   5.17) .AND. (B .LT.   6.78))) THEN
         CALL ENK4(122,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 280.00) .AND. (L .LT. 282.52)) .AND.
     &   ((B .LT.   5.17) .AND. (B .GE.  -0.87*L + 249.02))) THEN
         CALL ENK4(122,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 282.52) .AND. (L .LT. 285.01)) .AND.
     &   ((B .LT.   2.43) .AND. (B .GE.   0.19*L  -52.00))) THEN
         CALL ENK4(122,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 282.99) .AND. (L .LT. 283.90)) .AND.
     &   ((B .GE.  -3.27) .AND. (B .LT.  -1.67))) THEN
         CALL ENK4(123,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 282.60) .AND. (L .LT. 282.98)) .AND.
     &   ((B .LT.  -1.67) .AND. (B .GE.  -4.24*L + 1196.64))) THEN
         CALL ENK4(123,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 283.90) .AND. (L .LT. 284.97)) .AND.
     &   ((B .GE.  -2.53) .AND. (B .LT.  -0.79*L + 222.61))) THEN
         CALL ENK4(123,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 283.90) .AND. (L .LT. 284.97)) .AND.
     &   ((B .GE.  -3.27) .AND. (B .LT.  -2.53))) THEN
         CALL ENK4(123,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.77) .AND. (L .LT. 286.15)) .AND.
     &   ((B .GE.   0.08) .AND. (B .LT.   0.80*L  -228.51))) THEN
         CALL ENK4(124,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.33) .AND. (L .LT. 286.17)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.   0.11))) THEN
         CALL ENK4(124,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.33) .AND. (L .LT. 285.77)) .AND.
     &   ((B .LT.  -0.06) .AND. (B .GE.   0.82*L  -234.35))) THEN
         CALL ENK4(124,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.78) .AND. (L .LT. 285.33)) .AND.
     &   ((B .GE.  -0.39) .AND. (B .LT.   0.11))) THEN
         CALL ENK4(124,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 283.99) .AND. (L .LT. 284.78)) .AND.
     &   ((B .GE.  -0.39) .AND. (B .LT.  -0.08))) THEN
         CALL ENK4(124,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.96) .AND. (L .LT. 286.33)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.  -0.89))) THEN
         CALL ENK4(125,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.00) .AND. (L .LT. 284.97)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.   0.91*L  -260.15))) THEN
         CALL ENK4(125,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.96) .AND. (L .LT. 285.27)) .AND.
     &   ((B .GE.  -3.24) .AND. (B .LT.  -1.39))) THEN
         CALL ENK4(125,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.00) .AND. (L .LT. 284.96)) .AND.
     &   ((B .LT.  -1.77) .AND. (B .GE.  -0.79*L + 222.61))) THEN
         CALL ENK4(125,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.77) .AND. (L .LT. 286.17)) .AND.
     &   ((B .GE.  -0.23) .AND. (B .LT.  -0.08))) THEN
         CALL ENK4(126,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.76) .AND. (L .LT. 286.63)) .AND.
     &   ((B .LT.  -0.20) .AND. (B .GE.   0.68*L  -195.12))) THEN
         CALL ENK4(126,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.33) .AND. (L .LT. 285.77)) .AND.
     &   ((B .GE.  -0.44) .AND. (B .LT.   0.82*L  -234.35))) THEN
         CALL ENK4(126,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.99) .AND. (L .LT. 285.76)) .AND.
     &   ((B .GE.  -0.77) .AND. (B .LT.  -0.39))) THEN
         CALL ENK4(126,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.28) .AND. (L .LT. 286.38)) .AND.
     &   ((B .GE.  -2.32) .AND. (B .LT.  -0.84*L + 238.40))) THEN
         CALL ENK4(127,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.25) .AND. (L .LT. 286.35)) .AND.
     &   ((B .GE.  -3.55) .AND. (B .LT.  -2.32))) THEN
         CALL ENK4(127,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.26) .AND. (L .LT. 286.00)) .AND.
     &   ((B .LT.  -3.55) .AND. (B .GE.  -1.54*L + 435.65))) THEN
         CALL ENK4(127,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.98) .AND. (L .LT. 286.96)) .AND.
     &   ((B .GE.  -4.69) .AND. (B .LT.  -3.57))) THEN
         CALL ENK4(127,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.97) .AND. (L .LT. 287.93)) .AND.
     &   ((B .GE.  -4.26) .AND. (B .LT.  -3.55))) THEN
         CALL ENK4(127,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.96) .AND. (L .LT. 287.93)) .AND.
     &   ((B .LT.  -4.26) .AND. (B .GE.   0.42*L  -125.20))) THEN
         CALL ENK4(127,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.38) .AND. (L .LT. 288.62)) .AND.
     &   ((B .GE.  -2.32) .AND. (B .LT.   0.14*L  -42.38))) THEN
         CALL ENK4(128,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.96) .AND. (L .LT. 288.56)) .AND.
     &   ((B .GE.  -2.77) .AND. (B .LT.  -2.34))) THEN
         CALL ENK4(128,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.35) .AND. (L .LT. 287.96)) .AND.
     &   ((B .GE.  -3.55) .AND. (B .LT.  -2.32))) THEN
         CALL ENK4(128,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.59) .AND. (L .LT. 289.21)) .AND.
     &   ((B .GE.  -2.77) .AND. (B .LT.  -0.57*L + 162.18))) THEN
         CALL ENK4(128,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.93) .AND. (L .LT. 289.43)) .AND.
     &   ((B .LT.  -2.77) .AND. (B .GE.   0.53*L  -156.20))) THEN
         CALL ENK4(128,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.47) .AND. (L .LT. 289.80)) .AND.
     &   ((B .GE.  -0.68) .AND. (B .LT.  -0.64*L + 184.86))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.47) .AND. (L .LT. 289.80)) .AND.
     &   ((B .GE.  -1.03) .AND. (B .LT.  -0.68))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.60) .AND. (L .LT. 289.46)) .AND.
     &   ((B .GE.  -1.03) .AND. (B .LT.  -0.46))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.60) .AND. (L .LT. 289.01)) .AND.
     &   ((B .GE.  -1.44) .AND. (B .LT.  -1.03))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.09) .AND. (L .LT. 288.60)) .AND.
     &   ((B .GE.  -1.25) .AND. (B .LT.   0.36*L  -104.93))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.07) .AND. (L .LT. 288.59)) .AND.
     &   ((B .GE.  -2.01) .AND. (B .LT.  -1.25))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.38) .AND. (L .LT. 288.07)) .AND.
     &   ((B .GE.  -2.01) .AND. (B .LT.   0.14*L  -41.94))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.38) .AND. (L .LT. 288.62)) .AND.
     &   ((B .LT.  -2.01) .AND. (B .GE.   0.14*L  -42.38))) THEN
         CALL ENK4(129,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.36) .AND. (L .LT. 288.07)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.  -1.44))) THEN
         CALL ENK4(130,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.38) .AND. (L .LT. 288.07)) .AND.
     &   ((B .LT.  -1.77) .AND. (B .GE.   0.14*L  -41.94))) THEN
         CALL ENK4(130,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.17) .AND. (L .LT. 288.60)) .AND.
     &   ((B .GE.  -1.01) .AND. (B .LT.  -0.46))) THEN
         CALL ENK4(131,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.54) .AND. (L .LT. 288.32)) .AND.
     &   ((B .GE.  -0.77) .AND. (B .LT.   0.47*L  -135.85))) THEN
         CALL ENK4(131,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.35) .AND. (L .LT. 288.60)) .AND.
     &   ((B .LT.  -0.96) .AND. (B .GE.   0.36*L  -104.93))) THEN
         CALL ENK4(131,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.94) .AND. (L .LT. 287.35)) .AND.
     &   ((B .GE.  -1.46) .AND. (B .LT.  -0.96))) THEN
         CALL ENK4(131,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.94) .AND. (L .LT. 287.11)) .AND.
     &   ((B .GE.  -1.01) .AND. (B .LT.  -0.27))) THEN
         CALL ENK4(131,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.63) .AND. (L .LT. 286.94)) .AND.
     &   ((B .LT.  -0.30) .AND. (B .GE.  -3.89*L + 1114.65))) THEN
         CALL ENK4(131,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.13) .AND. (L .LT. 287.54)) .AND.
     &   ((B .GE.  -0.96) .AND. (B .LT.  -0.63))) THEN
         CALL ENK4(132,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.54) .AND. (L .LT. 288.17)) .AND.
     &   ((B .GE.  -0.96) .AND. (B .LT.  -0.77))) THEN
         CALL ENK4(132,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.54) .AND. (L .LT. 288.32)) .AND.
     &   ((B .LT.  -0.27) .AND. (B .GE.   0.47*L  -135.85))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.12) .AND. (L .LT. 288.27)) .AND.
     &   ((B .GE.  -0.27) .AND. (B .LT.  -0.45*L + 129.45))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.11) .AND. (L .LT. 287.57)) .AND.
     &   ((B .GE.  -0.63) .AND. (B .LT.   0.00))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.64) .AND. (L .LT. 287.12)) .AND.
     &   ((B .GE.   0.25) .AND. (B .LT.  -0.45*L + 129.45))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.63) .AND. (L .LT. 287.11)) .AND.
     &   ((B .GE.  -0.30) .AND. (B .LT.   0.25))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.16) .AND. (L .LT. 286.64)) .AND.
     &   ((B .GE.   0.42) .AND. (B .LT.  -0.45*L + 129.45))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.15) .AND. (L .LT. 286.63)) .AND.
     &   ((B .GE.  -0.20) .AND. (B .LT.   0.44))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.77) .AND. (L .LT. 286.15)) .AND.
     &   ((B .GE.   0.42) .AND. (B .LT.   0.68))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.01) .AND. (L .LT. 284.85)) .AND.
     &   ((B .GE.  -0.06) .AND. (B .LT.   0.68))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.75) .AND. (L .LT. 285.77)) .AND.
     &   ((B .GE.   0.11) .AND. (B .LT.   0.68))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 285.77) .AND. (L .LT. 286.15)) .AND.
     &   ((B .LT.   0.42) .AND. (B .GE.   0.80*L  -228.51))) THEN
         CALL ENK4(133,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.02) .AND. (L .LT. 290.00)) .AND.
     &   ((B .GE.   2.29) .AND. (B .LT.   4.69))) THEN
         CALL ENK4(134,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .GE.   1.70) .AND. (B .LT.  -1.05*L + 307.54))) THEN
         CALL ENK4(134,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .GE.   1.15) .AND. (B .LT.   1.70))) THEN
         CALL ENK4(134,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .LT.   1.15) .AND. (B .GE.   0.24*L  -68.68))) THEN
         CALL ENK4(134,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.97) .AND. (L .LT. 287.02)) .AND.
     &   ((B .GE.   3.84) .AND. (B .LT.   0.39*L  -107.20))) THEN
         CALL ENK4(134,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 284.97) .AND. (L .LT. 287.02)) .AND.
     &   ((B .GE.   0.68) .AND. (B .LT.   3.84))) THEN
         CALL ENK4(134,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.21) .AND. (L .LT. 290.00)) .AND.
     &   ((B .GE.   0.65) .AND. (B .LT.   2.32))) THEN
         CALL ENK4(135,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.00) .AND. (L .LT. 288.21)) .AND.
     &   ((B .GE.   0.65) .AND. (B .LT.   1.48))) THEN
         CALL ENK4(135,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.00) .AND. (L .LT. 288.21)) .AND.
     &   ((B .GE.   1.48) .AND. (B .LT.   0.67*L  -190.87))) THEN
         CALL ENK4(135,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.23) .AND. (L .LT. 290.00)) .AND.
     &   ((B .GE.  -0.20) .AND. (B .LT.   0.68))) THEN
         CALL ENK4(136,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.83) .AND. (L .LT. 290.60)) .AND.
     &   ((B .GE.  -0.61) .AND. (B .LT.  -0.20))) THEN
         CALL ENK4(136,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .GE.  -0.18) .AND. (B .LT.  -0.01))) THEN
         CALL ENK4(136,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.  -0.44*L + 128.31))) THEN
         CALL ENK4(136,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.16) .AND. (L .LT. 288.27)) .AND.
     &   ((B .LT.   0.65) .AND. (B .GE.  -0.45*L + 129.45))) THEN
         CALL ENK4(136,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.60) .AND. (L .LT. 291.10)) .AND.
     &   ((B .GE.  -0.87) .AND. (B .LT.  -0.47*L + 136.00))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.83) .AND. (L .LT. 290.60)) .AND.
     &   ((B .GE.  -1.03) .AND. (B .LT.  -0.61))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.81) .AND. (L .LT. 291.10)) .AND.
     &   ((B .LT.  -0.87) .AND. (B .GE.   1.14*L  -332.80))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.60) .AND. (L .LT. 290.81)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.  -0.89))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.46) .AND. (L .LT. 290.64)) .AND.
     &   ((B .GE.  -1.32) .AND. (B .LT.  -1.03))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.01) .AND. (L .LT. 289.46)) .AND.
     &   ((B .GE.  -1.46) .AND. (B .LT.  -1.03))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.97) .AND. (L .LT. 290.64)) .AND.
     &   ((B .LT.  -1.32) .AND. (B .GE.   0.42*L  -123.47))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.46) .AND. (L .LT. 289.97)) .AND.
     &   ((B .LT.  -1.32) .AND. (B .GE.  -0.61*L + 175.19))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.17) .AND. (L .LT. 289.46)) .AND.
     &   ((B .LT.  -1.46) .AND. (B .GE.   0.63*L  -183.79))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.60) .AND. (L .LT. 289.19)) .AND.
     &   ((B .GE.  -2.41) .AND. (B .LT.  -1.44))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.17) .AND. (L .LT. 289.46)) .AND.
     &   ((B .GE.  -2.27) .AND. (B .LT.  -2.41*L + 695.30))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.19) .AND. (L .LT. 289.40)) .AND.
     &   ((B .GE.  -2.77) .AND. (B .LT.  -2.27))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 288.59) .AND. (L .LT. 289.21)) .AND.
     &   ((B .LT.  -2.41) .AND. (B .GE.  -0.57*L + 162.18))) THEN
         CALL ENK4(137,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.96) .AND. (L .LT. 287.93)) .AND.
     &   ((B .GE.  -4.69) .AND. (B .LT.   0.42*L  -125.20))) THEN
         CALL ENK4(138,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.17) .AND. (L .LT. 287.94)) .AND.
     &   ((B .GE.  -6.02) .AND. (B .LT.  -4.69))) THEN
         CALL ENK4(138,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 286.17) .AND. (L .LT. 289.95)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -6.05))) THEN
         CALL ENK4(138,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.77) .AND. (L .LT. 290.41)) .AND.
     &   ((B .GE.  -3.57) .AND. (B .LT.  -2.55))) THEN
         CALL ENK4(139,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.93) .AND. (L .LT. 289.77)) .AND.
     &   ((B .GE.  -3.55) .AND. (B .LT.   0.53*L  -156.20))) THEN
         CALL ENK4(139,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 287.93) .AND. (L .LT. 290.37)) .AND.
     &   ((B .GE.  -6.02) .AND. (B .LT.  -3.55))) THEN
         CALL ENK4(139,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.97) .AND. (L .LT. 290.64)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.   0.42*L  -123.47))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.46) .AND. (L .LT. 289.97)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.  -0.61*L + 175.19))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.17) .AND. (L .LT. 289.46)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.   0.63*L  -183.79))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.17) .AND. (L .LT. 289.46)) .AND.
     &   ((B .LT.  -1.67) .AND. (B .GE.  -2.41*L + 695.30))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.92) .AND. (L .LT. 290.59)) .AND.
     &   ((B .GE.  -2.55) .AND. (B .LT.  -1.67))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.46) .AND. (L .LT. 289.91)) .AND.
     &   ((B .GE.  -2.24) .AND. (B .LT.  -1.67))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.43) .AND. (L .LT. 289.91)) .AND.
     &   ((B .GE.  -2.55) .AND. (B .LT.  -2.24))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.43) .AND. (L .LT. 289.77)) .AND.
     &   ((B .LT.  -2.55) .AND. (B .GE.   0.53*L  -156.20))) THEN
         CALL ENK4(140,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .GE.   0.87) .AND. (B .LT.   0.24*L  -68.68))) THEN
         CALL ENK4(141,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.16)) .AND.
     &   ((B .GE.   0.53) .AND. (B .LT.   0.87))) THEN
         CALL ENK4(141,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .LT.   0.53) .AND. (B .GE.  -0.44*L + 128.31))) THEN
         CALL ENK4(141,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.59) .AND. (L .LT. 291.45)) .AND.
     &   ((B .GE.  -2.60) .AND. (B .LT.  -1.22))) THEN
         CALL ENK4(142,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.39) .AND. (L .LT. 291.42)) .AND.
     &   ((B .GE.  -6.02) .AND. (B .LT.  -2.55))) THEN
         CALL ENK4(142,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.44) .AND. (L .LT. 291.57)) .AND.
     &   ((B .GE.  -2.24) .AND. (B .LT.  -1.22))) THEN
         CALL ENK4(142,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.53) .AND. (L .LT. 292.63)) .AND.
     &   ((B .GE.  -2.24) .AND. (B .LT.  -1.12*L + 325.48))) THEN
         CALL ENK4(142,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.43) .AND. (L .LT. 292.62)) .AND.
     &   ((B .GE.  -3.65) .AND. (B .LT.  -2.24))) THEN
         CALL ENK4(142,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.46) .AND. (L .LT. 291.18)) .AND.
     &   ((B .GE.  -0.58) .AND. (B .LT.  -0.20))) THEN
         CALL ENK4(143,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.60) .AND. (L .LT. 291.10)) .AND.
     &   ((B .LT.  -0.56) .AND. (B .GE.  -0.47*L + 136.00))) THEN
         CALL ENK4(143,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.10) .AND. (L .LT. 291.29)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.  -0.58))) THEN
         CALL ENK4(143,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.28) .AND. (L .LT. 291.56)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.  -2.47*L + 718.82))) THEN
         CALL ENK4(143,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.81) .AND. (L .LT. 291.10)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.   1.14*L  -332.80))) THEN
         CALL ENK4(143,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.54) .AND. (L .LT. 293.00)) .AND.
     &   ((B .GE.  -1.03) .AND. (B .LT.  -0.04))) THEN
         CALL ENK4(144,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.28) .AND. (L .LT. 291.51)) .AND.
     &   ((B .LT.  -0.61) .AND. (B .GE.  -2.47*L + 718.82))) THEN
         CALL ENK4(144,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 292.64) .AND. (L .LT. 293.33)) .AND.
     &   ((B .GE.   1.63) .AND. (B .LT.   2.53))) THEN
         CALL ENK4(145,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.20) .AND. (L .LT. 295.03)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   3.93))) THEN
         CALL ENK4(146,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.20)) .AND.
     &   ((B .GE.   2.93) .AND. (B .LT.   3.93))) THEN
         CALL ENK4(146,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 290.00) .AND. (L .LT. 291.18)) .AND.
     &   ((B .LT.   2.93) .AND. (B .GE.  -1.05*L + 307.54))) THEN
         CALL ENK4(146,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 292.63) .AND. (L .LT. 295.03)) .AND.
     &   ((B .GE.  -3.65) .AND. (B .LT.  -1.96))) THEN
         CALL ENK4(147,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 292.63) .AND. (L .LT. 293.02)) .AND.
     &   ((B .GE.  -1.96) .AND. (B .LT.  -1.03))) THEN
         CALL ENK4(147,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.53) .AND. (L .LT. 292.63)) .AND.
     &   ((B .LT.  -1.03) .AND. (B .GE.  -1.12*L + 325.48))) THEN
         CALL ENK4(147,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 291.42) .AND. (L .LT. 295.03)) .AND.
     &   ((B .GE.  -6.02) .AND. (B .LT.  -3.65))) THEN
         CALL ENK4(147,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 289.95) .AND. (L .LT. 295.03)) .AND.
     &   ((B .GE.  -7.64) .AND. (B .LT.  -6.02))) THEN
         CALL ENK4(147,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 293.00) .AND. (L .LT. 295.06)) .AND.
     &   ((B .GE.  -1.08) .AND. (B .LT.  -0.04))) THEN
         CALL ENK4(148,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 293.03) .AND. (L .LT. 293.84)) .AND.
     &   ((B .GE.  -1.96) .AND. (B .LT.  -1.08))) THEN
         CALL ENK4(148,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 293.84) .AND. (L .LT. 294.32)) .AND.
     &   ((B .LT.  -1.08) .AND. (B .GE.   1.71*L  -504.40))) THEN
         CALL ENK4(148,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 293.02) .AND. (L .LT. 293.84)) .AND.
     &   ((B .GE.  -1.96) .AND. (B .LT.  -1.08))) THEN
         CALL ENK4(149,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 294.32) .AND. (L .LT. 295.06)) .AND.
     &   ((B .GE.  -1.96) .AND. (B .LT.  -1.08))) THEN
         CALL ENK4(149,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 293.84) .AND. (L .LT. 294.32)) .AND.
     &   ((B .GE.  -1.96) .AND. (B .LT.   1.71*L  -504.40))) THEN
         CALL ENK4(149,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.05) .AND. (L .LT. 295.43)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.  -1.12*L + 329.38))) THEN
         CALL ENK4(149,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.02) .AND. (L .LT. 295.59)) .AND.
     &   ((B .GE.  -2.41) .AND. (B .LT.  -1.58))) THEN
         CALL ENK4(149,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.59) .AND. (L .LT. 297.00)) .AND.
     &   ((B .GE.  -2.41) .AND. (B .LT.  -2.03))) THEN
         CALL ENK4(150,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.04) .AND. (L .LT. 297.00)) .AND.
     &   ((B .GE.  -7.61) .AND. (B .LT.  -2.41))) THEN
         CALL ENK4(150,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 296.82) .AND. (L .LT. 297.11)) .AND.
     &   ((B .LT.  -0.96) .AND. (B .GE.   2.12*L  -630.95))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.03) .AND. (L .LT. 297.11)) .AND.
     &   ((B .GE.  -0.96) .AND. (B .LT.  -0.46*L + 135.70))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 296.39) .AND. (L .LT. 296.82)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.  -0.96))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.05) .AND. (L .LT. 296.39)) .AND.
     &   ((B .GE.  -1.15) .AND. (B .LT.  -0.96))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.84) .AND. (L .LT. 296.39)) .AND.
     &   ((B .LT.  -1.15) .AND. (B .GE.  -0.78*L + 229.59))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.43) .AND. (L .LT. 295.84)) .AND.
     &   ((B .LT.  -1.15) .AND. (B .GE.   1.04*L  -308.87))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.05) .AND. (L .LT. 295.43)) .AND.
     &   ((B .LT.  -1.15) .AND. (B .GE.  -1.12*L + 329.38))) THEN
         CALL ENK4(151,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.03) .AND. (L .LT. 298.03)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   1.34))) THEN
         CALL ENK4(152,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.03) .AND. (L .LT. 298.51)) .AND.
     &   ((B .LT.  -0.01) .AND. (B .GE.  -0.46*L + 135.70))) THEN
         CALL ENK4(152,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.52) .AND. (L .LT. 300.00)) .AND.
     &   ((B .GE.  -0.94) .AND. (B .LT.  -0.01))) THEN
         CALL ENK4(152,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.51) .AND. (L .LT. 300.00)) .AND.
     &   ((B .LT.  -0.94) .AND. (B .GE.   0.38*L  -114.97))) THEN
         CALL ENK4(152,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.03) .AND. (L .LT. 300.00)) .AND.
     &   ((B .GE.   2.96) .AND. (B .LT.   7.52))) THEN
         CALL ENK4(153,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.03) .AND. (L .LT. 298.04)) .AND.
     &   ((B .GE.   1.32) .AND. (B .LT.   2.96))) THEN
         CALL ENK4(153,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.04) .AND. (L .LT. 298.77)) .AND.
     &   ((B .GE.   1.94) .AND. (B .LT.  -1.40*L + 420.20))) THEN
         CALL ENK4(154,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.04) .AND. (L .LT. 298.77)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   1.94))) THEN
         CALL ENK4(154,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.51) .AND. (L .LT. 299.52)) .AND.
     &   ((B .GE.  -1.60) .AND. (B .LT.   0.38*L  -114.97))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 297.11) .AND. (L .LT. 298.51)) .AND.
     &   ((B .GE.  -1.60) .AND. (B .LT.  -0.46*L + 135.70))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 297.02) .AND. (L .LT. 299.54)) .AND.
     &   ((B .GE.  -2.46) .AND. (B .LT.  -1.60))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.59) .AND. (L .LT. 297.01)) .AND.
     &   ((B .GE.  -2.03) .AND. (B .LT.  -1.58))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 297.00) .AND. (L .LT. 299.05)) .AND.
     &   ((B .GE.  -3.48) .AND. (B .LT.  -1.58))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 296.82) .AND. (L .LT. 297.11)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.   2.12*L  -630.95))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.84) .AND. (L .LT. 296.39)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.  -0.78*L + 229.59))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 295.43) .AND. (L .LT. 295.84)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.   1.04*L  -308.87))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 299.05) .AND. (L .LT. 299.97)) .AND.
     &   ((B .GE.  -5.45) .AND. (B .LT.  -2.46))) THEN
         CALL ENK4(155,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 299.54) .AND. (L .LT. 300.00)) .AND.
     &   ((B .GE.  -1.20) .AND. (B .LT.   0.38*L  -114.97))) THEN
         CALL ENK4(156,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 299.54) .AND. (L .LT. 300.00)) .AND.
     &   ((B .GE.  -2.46) .AND. (B .LT.  -1.20))) THEN
         CALL ENK4(156,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.77) .AND. (L .LT. 300.00)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.   3.00))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 298.04) .AND. (L .LT. 298.77)) .AND.
     &   ((B .LT.   3.00) .AND. (B .GE.  -1.40*L + 420.20))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP261B(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 300.40) .AND. (L .LT. 301.00)) .AND.
     &  ((B .LT.  -0.11) .AND. (B .GE.   1.35*L  -406.57))) THEN
        CALL ENK4(156,D,AV,A0)
        RETURN
      ENDIF
      IF (((L .GE. 300.00) .AND. (L .LT. 300.40)) .AND.
     &   ((B .GE.  -1.11) .AND. (B .LT.  -0.11))) THEN
         CALL ENK4(156,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.40) .AND. (L .LT. 301.04)) .AND.
     &   ((B .GE.  -1.11) .AND. (B .LT.  -0.63*L + 189.00))) THEN
         CALL ENK4(156,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.00) .AND. (L .LT. 301.00)) .AND.
     &   ((B .GE.  -3.04) .AND. (B .LT.  -1.11))) THEN
         CALL ENK4(156,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.77) .AND. (L .LT. 304.05)) .AND.
     &   ((B .LT.   4.51) .AND. (B .GE.   3.67*L  -1111.50))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 302.65) .AND. (L .LT. 303.77)) .AND.
     &   ((B .GE.   2.54) .AND. (B .LT.   4.51))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.70) .AND. (L .LT. 302.65)) .AND.
     &   ((B .GE.   2.54) .AND. (B .LT.   4.51))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.00) .AND. (L .LT. 300.70)) .AND.
     &   ((B .GE.   0.70) .AND. (B .LT.   4.51))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.68) .AND. (L .LT. 301.14)) .AND.
     &   ((B .LT.   2.54) .AND. (B .GE.   2.46*L  -738.07))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.63) .AND. (L .LT. 301.54)) .AND.
     &   ((B .GE.   0.73) .AND. (B .LT.  -0.74*L + 223.90))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.98) .AND. (L .LT. 301.54)) .AND.
     &   ((B .LT.   0.73) .AND. (B .GE.   1.35*L  -406.57))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.00) .AND. (L .LT. 300.98)) .AND.
     &   ((B .GE.  -0.11) .AND. (B .LT.   0.73))) THEN
         CALL ENK4(157,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 301.54) .AND. (L .LT. 302.81)) .AND.
     &   ((B .LT.   2.56) .AND. (B .GE.   1.35*L  -406.57))) THEN
         CALL ENK4(158,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 301.14) .AND. (L .LT. 301.55)) .AND.
     &   ((B .GE.   1.42) .AND. (B .LT.   2.54))) THEN
         CALL ENK4(158,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.68) .AND. (L .LT. 301.14)) .AND.
     &   ((B .GE.   1.42) .AND. (B .LT.   2.46*L  -738.07))) THEN
         CALL ENK4(158,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.63) .AND. (L .LT. 301.54)) .AND.
     &   ((B .LT.   1.42) .AND. (B .GE.  -0.74*L + 223.90))) THEN
         CALL ENK4(158,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.25) .AND. (L .LT. 304.54)) .AND.
     &   ((B .GE.   2.53) .AND. (B .LT.  -1.09*L + 334.45))) THEN
         CALL ENK4(159,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.45) .AND. (L .LT. 304.25)) .AND.
     &   ((B .GE.   2.53) .AND. (B .LT.   0.39*L  -115.93))) THEN
         CALL ENK4(159,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 302.71) .AND. (L .LT. 304.54)) .AND.
     &   ((B .GE.   1.70) .AND. (B .LT.   2.53))) THEN
         CALL ENK4(159,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 302.20) .AND. (L .LT. 302.81)) .AND.
     &   ((B .GE.   1.70) .AND. (B .LT.   1.35*L  -406.57))) THEN
         CALL ENK4(159,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.53) .AND. (L .LT. 305.01)) .AND.
     &   ((B .GE.   1.70) .AND. (B .LT.   1.93))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.01) .AND. (L .LT. 305.34)) .AND.
     &   ((B .GE.   1.17) .AND. (B .LT.  -2.31*L + 706.53))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.01) .AND. (L .LT. 305.34)) .AND.
     &   ((B .GE.   0.86) .AND. (B .LT.   1.17))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.78) .AND. (L .LT. 305.01)) .AND.
     &   ((B .GE.   0.86) .AND. (B .LT.   1.69))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.78) .AND. (L .LT. 305.50)) .AND.
     &   ((B .GE.  -0.45) .AND. (B .LT.   0.86))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.76) .AND. (L .LT. 305.02)) .AND.
     &   ((B .LT.  -0.45) .AND. (B .GE.   0.98*L  -299.42))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 302.20) .AND. (L .LT. 303.77)) .AND.
     &   ((B .GE.  -0.71) .AND. (B .LT.   1.72))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.40) .AND. (L .LT. 302.20)) .AND.
     &   ((B .GE.  -0.71) .AND. (B .LT.   1.35*L  -406.57))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.04) .AND. (L .LT. 303.76)) .AND.
     &   ((B .GE.  -1.69) .AND. (B .LT.  -0.71))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.40) .AND. (L .LT. 303.04)) .AND.
     &   ((B .LT.  -0.71) .AND. (B .GE.  -0.63*L + 189.00))) THEN
         CALL ENK4(160,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.76) .AND. (L .LT. 305.02)) .AND.
     &   ((B .GE.  -1.65) .AND. (B .LT.   0.98*L  -299.42))) THEN
         CALL ENK5(161,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 301.00) .AND. (L .LT. 303.04)) .AND.
     &   ((B .GE.  -1.71) .AND. (B .LT.  -0.63*L + 189.00))) THEN
         CALL ENK5(162,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 301.00) .AND. (L .LT. 303.02)) .AND.
     &   ((B .GE.  -3.09) .AND. (B .LT.  -1.71))) THEN
         CALL ENK5(162,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.02) .AND. (L .LT. 305.99)) .AND.
     &   ((B .GE.  -2.12) .AND. (B .LT.  -1.43))) THEN
         CALL ENK5(163,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.82) .AND. (L .LT. 305.05)) .AND.
     &   ((B .GE.  -2.12) .AND. (B .LT.  -1.65))) THEN
         CALL ENK5(163,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.40) .AND. (L .LT. 304.84)) .AND.
     &   ((B .GE.  -2.12) .AND. (B .LT.   1.04*L  -318.73))) THEN
         CALL ENK5(163,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.43) .AND. (L .LT. 305.98)) .AND.
     &   ((B .GE.  -3.08) .AND. (B .LT.  -2.13))) THEN
         CALL ENK5(163,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 300.00) .AND. (L .LT. 306.36)) .AND.
     &   ((B .GE.  -7.67) .AND. (B .LT.  -3.04))) THEN
         CALL ENK5(164,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.02) .AND. (L .LT. 306.00)) .AND.
     &   ((B .GE.  -1.50) .AND. (B .LT.  -0.48))) THEN
         CALL ENK5(165,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 307.58) .AND. (L .LT. 308.05)) .AND.
     &   ((B .GE.   0.23) .AND. (B .LT.  -1.25*L + 385.31))) THEN
         CALL ENK5(166,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.50) .AND. (L .LT. 307.57)) .AND.
     &   ((B .GE.   0.24) .AND. (B .LT.   0.83))) THEN
         CALL ENK5(166,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 306.00) .AND. (L .LT. 307.97)) .AND.
     &   ((B .GE.  -0.98) .AND. (B .LT.   0.23))) THEN
         CALL ENK5(166,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.50) .AND. (L .LT. 306.00)) .AND.
     &   ((B .GE.  -0.46) .AND. (B .LT.   0.24))) THEN
         CALL ENK5(166,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.36) .AND. (L .LT. 306.21)) .AND.
     &   ((B .GE.   1.14) .AND. (B .LT.   3.53))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.54) .AND. (L .LT. 305.36)) .AND.
     &   ((B .GE.   1.93) .AND. (B .LT.   3.53))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.77) .AND. (L .LT. 304.56)) .AND.
     &   ((B .GE.   2.84) .AND. (B .LT.   3.53))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 307.58) .AND. (L .LT. 308.59)) .AND.
     &   ((B .GE.   0.25) .AND. (B .LT.   1.14))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 303.45) .AND. (L .LT. 304.25)) .AND.
     &   ((B .LT.   2.84) .AND. (B .GE.   0.39*L  -115.93))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 304.25) .AND. (L .LT. 304.56)) .AND.
     &   ((B .LT.   2.84) .AND. (B .GE.  -1.09*L + 334.45))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.01) .AND. (L .LT. 305.36)) .AND.
     &   ((B .LT.   1.93) .AND. (B .GE.  -2.31*L + 706.53))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 306.22) .AND. (L .LT. 308.59)) .AND.
     &   ((B .GE.   1.14) .AND. (B .LT.  -0.50*L + 155.51))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 305.34) .AND. (L .LT. 308.05)) .AND.
     &   ((B .GE.   0.83) .AND. (B .LT.   1.18))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 307.58) .AND. (L .LT. 308.05)) .AND.
     &   ((B .LT.   0.83) .AND. (B .GE.  -1.25*L + 385.31))) THEN
         CALL ENK5(167,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 306.23) .AND. (L .LT. 307.23)) .AND.
     &   ((B .GE.   2.31) .AND. (B .LT.   3.53))) THEN
         CALL ENK5(168,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 306.22) .AND. (L .LT. 307.23)) .AND.
     &   ((B .LT.   2.31) .AND. (B .GE.  -0.50*L + 155.51))) THEN
         CALL ENK5(168,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 308.00) .AND. (L .LT. 309.98)) .AND.
     &   ((B .GE.  -3.06) .AND. (B .LT.  -0.06))) THEN
         CALL ENK5(169,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 306.00) .AND. (L .LT. 307.96)) .AND.
     &   ((B .GE.  -3.08) .AND. (B .LT.  -0.96))) THEN
         CALL ENK5(169,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 306.39) .AND. (L .LT. 309.98)) .AND.
     &   ((B .GE.  -7.66) .AND. (B .LT.  -3.07))) THEN
         CALL ENK5(169,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 308.59) .AND. (L .LT. 310.00)) .AND.
     &   ((B .GE.   0.23) .AND. (B .LT.   6.91))) THEN
         CALL ENK5(170,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 307.23) .AND. (L .LT. 308.64)) .AND.
     &   ((B .GE.   1.88) .AND. (B .LT.   6.94))) THEN
         CALL ENK5(170,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 307.23) .AND. (L .LT. 308.59)) .AND.
     &   ((B .LT.   1.88) .AND. (B .GE.  -0.50*L + 155.51))) THEN
         CALL ENK5(170,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 310.00) .AND. (L .LT. 311.26)) .AND.
     &   ((B .GE.  -2.92) .AND. (B .LT.   3.78))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 311.16) .AND. (L .LT. 312.49)) .AND.
     &   ((B .GE.   0.96) .AND. (B .LT.  -2.21*L + 691.39))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 311.16) .AND. (L .LT. 312.85)) .AND.
     &   ((B .GE.  -1.85) .AND. (B .LT.   0.96))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 312.85) .AND. (L .LT. 313.39)) .AND.
     &   ((B .GE.  -0.02) .AND. (B .LT.   0.96))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 312.85) .AND. (L .LT. 313.63)) .AND.
     &   ((B .GE.  -1.00) .AND. (B .LT.  -0.02))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 311.16) .AND. (L .LT. 312.85)) .AND.
     &   ((B .LT.  -1.85) .AND. (B .GE.   0.63*L  -199.02))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 312.85) .AND. (L .LT. 313.63)) .AND.
     &   ((B .LT.  -1.00) .AND. (B .GE.   1.11*L  -349.15))) THEN
         CALL ENK5(171,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 311.16) .AND. (L .LT. 312.87)) .AND.
     &   ((B .GE.  -2.93) .AND. (B .LT.   0.63*L  -199.02))) THEN
         CALL ENK5(172,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 309.98) .AND. (L .LT. 312.87)) .AND.
     &   ((B .GE.  -7.65) .AND. (B .LT.  -2.93))) THEN
         CALL ENK5(172,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 313.63) .AND. (L .LT. 315.57)) .AND.
     &   ((B .GE.  -1.86) .AND. (B .LT.  -1.00))) THEN
         CALL ENK5(173,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 312.85) .AND. (L .LT. 313.63)) .AND.
     &   ((B .GE.  -1.86) .AND. (B .LT.   1.11*L  -349.15))) THEN
         CALL ENK5(173,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 312.85) .AND. (L .LT. 315.57)) .AND.
     &   ((B .GE.  -7.01) .AND. (B .LT.  -1.86))) THEN
         CALL ENK5(173,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 313.35) .AND. (L .LT. 314.59)) .AND.
     &   ((B .GE.   0.93) .AND. (B .LT.   0.64*L  -199.55))) THEN
         CALL ENK5(174,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 313.35) .AND. (L .LT. 314.60)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.93))) THEN
         CALL ENK5(174,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 314.60) .AND. (L .LT. 314.98)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.48))) THEN
         CALL ENK5(174,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 313.58) .AND. (L .LT. 314.97)) .AND.
     &   ((B .GE.  -1.00) .AND. (B .LT.  -0.07))) THEN
         CALL ENK5(174,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 314.98) .AND. (L .LT. 315.55)) .AND.
     &   ((B .GE.  -0.45) .AND. (B .LT.  -0.67*L + 211.05))) THEN
         CALL ENK5(174,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 314.98) .AND. (L .LT. 315.55)) .AND.
     &   ((B .GE.  -1.00) .AND. (B .LT.  -0.45))) THEN
         CALL ENK5(174,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 310.00) .AND. (L .LT. 315.40)) .AND.
     &   ((B .GE.   3.76) .AND. (B .LT.   5.94))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 313.35) .AND. (L .LT. 315.40)) .AND.
     &   ((B .GE.   1.72) .AND. (B .LT.   3.77))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 312.49) .AND. (L .LT. 313.35)) .AND.
     &   ((B .GE.   0.96) .AND. (B .LT.   3.77))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 311.16) .AND. (L .LT. 312.49)) .AND.
     &   ((B .LT.   3.78) .AND. (B .GE.  -2.21*L + 691.39))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 313.35) .AND. (L .LT. 314.59)) .AND.
     &   ((B .LT.   1.72) .AND. (B .GE.   0.64*L  -199.55))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 314.59) .AND. (L .LT. 315.40)) .AND.
     &   ((B .GE.   1.00) .AND. (B .LT.   1.72))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 314.58) .AND. (L .LT. 314.96)) .AND.
     &   ((B .GE.   0.48) .AND. (B .LT.   1.00))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 314.96) .AND. (L .LT. 316.59)) .AND.
     &   ((B .GE.  -0.05) .AND. (B .LT.   1.02))) THEN
         CALL ENK5(175,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 317.99) .AND. (L .LT. 318.82)) .AND.
     &   ((B .GE.  -1.25) .AND. (B .LT.  -1.04))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.58) .AND. (L .LT. 318.82)) .AND.
     &   ((B .GE.  -3.40) .AND. (B .LT.  -1.25))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.54) .AND. (L .LT. 318.12)) .AND.
     &   ((B .GE.  -5.64) .AND. (B .LT.  -3.40))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.12) .AND. (L .LT. 318.81)) .AND.
     &   ((B .GE.  -4.07) .AND. (B .LT.  -3.40))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.82) .AND. (L .LT. 319.97)) .AND.
     &   ((B .GE.  -4.07) .AND. (B .LT.  -0.56*L + 175.15))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.10) .AND. (L .LT. 319.97)) .AND.
     &   ((B .LT.  -4.07) .AND. (B .GE.   0.84*L  -272.85))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.52) .AND. (L .LT. 318.13)) .AND.
     &   ((B .LT.  -5.64) .AND. (B .GE.   0.15*L  -53.30))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.84) .AND. (L .LT. 319.00)) .AND.
     &   ((B .GE.  -2.14) .AND. (B .LT.  -1.04))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.00) .AND. (L .LT. 320.00)) .AND.
     &   ((B .GE.  -2.14) .AND. (B .LT.  -1.09*L + 346.72))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.82) .AND. (L .LT. 320.00)) .AND.
     &   ((B .LT.  -3.04) .AND. (B .GE.   0.32*L  -105.62))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.83) .AND. (L .LT. 320.00)) .AND.
     &   ((B .GE.  -3.04) .AND. (B .LT.  -2.14))) THEN
         CALL ENK5(176,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 328.42) .AND. (L .LT. 330.00)) .AND.
     &   ((B .GE.  -3.95) .AND. (B .LT.   0.49*L  -164.82))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 328.42) .AND. (L .LT. 330.00)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -3.95))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 326.06) .AND. (L .LT. 328.40)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -3.93))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.41) .AND. (L .LT. 326.06)) .AND.
     &   ((B .GE.  -4.43) .AND. (B .LT.   0.73*L  -241.98))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.33) .AND. (L .LT. 325.41)) .AND.
     &   ((B .GE.  -4.43) .AND. (B .LT.  -0.21*L +  63.89))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.99) .AND. (L .LT. 322.33)) .AND.
     &   ((B .GE.  -4.43) .AND. (B .LT.  -3.78))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.10) .AND. (L .LT. 319.97)) .AND.
     &   ((B .GE.  -5.64) .AND. (B .LT.   0.84*L  -272.85))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.13) .AND. (L .LT. 319.97)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -5.64))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.52) .AND. (L .LT. 318.13)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.   0.15*L  -53.30))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.97) .AND. (L .LT. 325.41)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -4.43))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.41) .AND. (L .LT. 326.06)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -4.40))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.55) .AND. (L .LT. 330.00)) .AND.
     &   ((B .LT.  -6.03) .AND. (B .GE.  -0.11*L +  28.68))) THEN
         CALL ENK5(177,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 320.00) .AND. (L .LT. 320.74)) .AND.
     &   ((B .GE.  -3.78) .AND. (B .LT.  -1.80))) THEN
         CALL ENK5(178,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.82) .AND. (L .LT. 320.00)) .AND.
     &   ((B .GE.  -3.42) .AND. (B .LT.   0.32*L  -105.62))) THEN
         CALL ENK5(178,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.61) .AND. (L .LT. 320.00)) .AND.
     &   ((B .GE.  -3.78) .AND. (B .LT.  -3.42))) THEN
         CALL ENK5(178,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.82) .AND. (L .LT. 319.63)) .AND.
     &   ((B .LT.  -3.42) .AND. (B .GE.  -0.56*L + 175.15))) THEN
         CALL ENK5(178,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.41) .AND. (L .LT. 320.01)) .AND.
     &   ((B .GE.   0.97) .AND. (B .LT.   3.95))) THEN
         CALL ENK5(179,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 318.01) .AND. (L .LT. 319.01)) .AND.
     &   ((B .GE.  -1.04) .AND. (B .LT.   0.99))) THEN
         CALL ENK5(179,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 316.60) .AND. (L .LT. 318.01)) .AND.
     &   ((B .GE.  -1.28) .AND. (B .LT.   0.99))) THEN
         CALL ENK5(179,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 315.60) .AND. (L .LT. 316.60)) .AND.
     &   ((B .GE.  -1.25) .AND. (B .LT.  -0.05))) THEN
         CALL ENK5(179,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.77) .AND. (L .LT. 321.67)) .AND.
     &   ((B .GE.  -1.83) .AND. (B .LT.   0.99))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.00) .AND. (L .LT. 319.77)) .AND.
     &   ((B .GE.  -1.04) .AND. (B .LT.   0.97))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 319.00) .AND. (L .LT. 319.77)) .AND.
     &   ((B .LT.  -1.04) .AND. (B .GE.  -1.09*L + 346.72))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 320.74) .AND. (L .LT. 321.67)) .AND.
     &   ((B .GE.  -2.11) .AND. (B .LT.  -1.80))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 321.69) .AND. (L .LT. 322.42)) .AND.
     &   ((B .GE.  -2.14) .AND. (B .LT.   0.01))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 321.69) .AND. (L .LT. 322.41)) .AND.
     &   ((B .LT.  -2.14) .AND. (B .GE.  -0.68*L + 216.64))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.42) .AND. (L .LT. 322.82)) .AND.
     &   ((B .GE.  -2.64) .AND. (B .LT.  -1.61))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.83) .AND. (L .LT. 323.35)) .AND.
     &   ((B .GE.  -2.64) .AND. (B .LT.  -1.97*L + 634.29))) THEN
         CALL ENK5(180,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.04) .AND. (L .LT. 325.41)) .AND.
     &   ((B .GE.  -2.92) .AND. (B .LT.  -2.47*L + 800.98))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.04) .AND. (L .LT. 325.41)) .AND.
     &   ((B .GE.  -3.78) .AND. (B .LT.  -2.92))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 324.63) .AND. (L .LT. 325.04)) .AND.
     &   ((B .GE.  -3.78) .AND. (B .LT.  -1.99))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 324.00) .AND. (L .LT. 324.63)) .AND.
     &   ((B .GE.  -2.64) .AND. (B .LT.   0.98*L  -320.15))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 320.74) .AND. (L .LT. 324.62)) .AND.
     &   ((B .GE.  -3.83) .AND. (B .LT.  -2.64))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.33) .AND. (L .LT. 325.41)) .AND.
     &   ((B .LT.  -3.78) .AND. (B .GE.  -0.21*L +  63.89))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 321.64) .AND. (L .LT. 322.41)) .AND.
     &   ((B .GE.  -2.64) .AND. (B .LT.  -0.68*L + 216.64))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 320.74) .AND. (L .LT. 321.64)) .AND.
     &   ((B .GE.  -2.64) .AND. (B .LT.  -2.14))) THEN
         CALL ENK5(181,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.45) .AND. (L .LT. 323.36)) .AND.
     &   ((B .GE.  -1.61) .AND. (B .LT.  -0.42*L + 134.18))) THEN
         CALL ENK5(182,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.83) .AND. (L .LT. 323.36)) .AND.
     &   ((B .LT.  -1.61) .AND. (B .GE.  -1.97*L + 634.29))) THEN
         CALL ENK5(182,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 323.36) .AND. (L .LT. 324.44)) .AND.
     &   ((B .GE.  -2.16) .AND. (B .LT.  -0.42*L + 134.18))) THEN
         CALL ENK5(182,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 323.35) .AND. (L .LT. 324.00)) .AND.
     &   ((B .GE.  -2.64) .AND. (B .LT.  -2.16))) THEN
         CALL ENK5(182,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 324.00) .AND. (L .LT. 324.44)) .AND.
     &   ((B .LT.  -2.16) .AND. (B .GE.   0.98*L  -320.15))) THEN
         CALL ENK5(182,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.43) .AND. (L .LT. 325.05)) .AND.
     &   ((B .GE.  -1.23) .AND. (B .LT.  -0.01))) THEN
         CALL ENK5(183,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 324.61) .AND. (L .LT. 325.05)) .AND.
     &   ((B .GE.  -1.99) .AND. (B .LT.  -1.23))) THEN
         CALL ENK5(183,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 324.44) .AND. (L .LT. 324.63)) .AND.
     &   ((B .LT.  -2.02) .AND. (B .GE.   0.98*L  -320.15))) THEN
         CALL ENK5(183,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 324.42) .AND. (L .LT. 324.63)) .AND.
     &   ((B .GE.  -2.01) .AND. (B .LT.  -1.23))) THEN
         CALL ENK5(183,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 322.45) .AND. (L .LT. 324.44)) .AND.
     &   ((B .LT.  -1.23) .AND. (B .GE.  -0.42*L + 134.18))) THEN
         CALL ENK5(183,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 321.69) .AND. (L .LT. 325.07)) .AND.
     &   ((B .GE.  -0.01) .AND. (B .LT.   5.98))) THEN
         CALL ENK5(184,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 320.00) .AND. (L .LT. 321.73)) .AND.
     &   ((B .GE.   0.94) .AND. (B .LT.   5.98))) THEN
         CALL ENK5(184,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 328.24) .AND. (L .LT. 329.09)) .AND.
     &   ((B .GE.  -2.35) .AND. (B .LT.   0.45*L  -150.10))) THEN
         CALL ENK5(185,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 328.43) .AND. (L .LT. 329.05)) .AND.
     &   ((B .GE.  -3.14) .AND. (B .LT.  -2.35))) THEN
         CALL ENK5(185,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 329.06) .AND. (L .LT. 330.00)) .AND.
     &   ((B .GE.  -3.14) .AND. (B .LT.  -0.84*L + 274.13))) THEN
         CALL ENK5(185,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 328.42) .AND. (L .LT. 330.00)) .AND.
     &   ((B .LT.  -3.14) .AND. (B .GE.   0.49*L  -164.82))) THEN
         CALL ENK5(185,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.41) .AND. (L .LT. 328.42)) .AND.
     &   ((B .GE.  -3.93) .AND. (B .LT.  -2.35))) THEN
         CALL ENK5(185,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 329.06) .AND. (L .LT. 330.00)) .AND.
     &   ((B .GE.  -2.37) .AND. (B .LT.  -1.97))) THEN
         CALL ENK5(186,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 329.06) .AND. (L .LT. 330.00)) .AND.
     &   ((B .LT.  -2.35) .AND. (B .GE.  -0.84*L + 274.12))) THEN
         CALL ENK5(186,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.07) .AND. (L .LT. 330.00)) .AND.
     &   ((B .GE.  -1.97) .AND. (B .LT.   1.04))) THEN
         CALL ENK5(187,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.02) .AND. (L .LT. 325.47)) .AND.
     &   ((B .LT.  -1.97) .AND. (B .GE.  -1.00*L + 323.15))) THEN
         CALL ENK5(187,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 325.47) .AND. (L .LT. 328.24)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -1.97))) THEN
         CALL ENK5(187,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 328.24) .AND. (L .LT. 329.09)) .AND.
     &   ((B .LT.  -1.97) .AND. (B .GE.   0.45*L  -150.10))) THEN
         CALL ENK5(187,D,AV,A0)
         RETURN
      ENDIF
      END
C ---------------------------------------------------------------

      SUBROUTINE MAP261C(L,B,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      IF (((L .GE. 330.00) .AND. (L .LT. 331.45)) .AND.
     &   ((B .GE.  -2.20) .AND. (B .LT.  -0.77))) THEN
         CALL ENK5(188,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.45) .AND. (L .LT. 332.03)) .AND.
     &   ((B .GE.  -1.24) .AND. (B .LT.  -0.96))) THEN
         CALL ENK5(189,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.45) .AND. (L .LT. 332.03)) .AND.
     &   ((B .LT.  -1.24) .AND. (B .GE.   0.55*L  -183.97))) THEN
         CALL ENK5(189,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 333.66) .AND. (L .LT. 335.03)) .AND.
     &   ((B .GE.   0.62) .AND. (B .LT.  -0.86*L + 288.77))) THEN
         CALL ENK5(190,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 332.18) .AND. (L .LT. 333.66)) .AND.
     &   ((B .GE.   0.65) .AND. (B .LT.   0.79*L  -261.70))) THEN
         CALL ENK5(190,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 332.18) .AND. (L .LT. 335.02)) .AND.
     &   ((B .GE.  -0.96) .AND. (B .LT.   0.65))) THEN
         CALL ENK5(190,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.45) .AND. (L .LT. 332.18)) .AND.
     &   ((B .GE.  -0.96) .AND. (B .LT.   0.29))) THEN
         CALL ENK5(190,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 330.45) .AND. (L .LT. 331.45)) .AND.
     &   ((B .GE.  -0.77) .AND. (B .LT.   0.26))) THEN
         CALL ENK5(190,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 330.45) .AND. (L .LT. 331.03)) .AND.
     &   ((B .GE.   0.29) .AND. (B .LT.   0.65))) THEN
         CALL ENK5(190,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 332.38) .AND. (L .LT. 333.01)) .AND.
     &   ((B .GE.  -2.44) .AND. (B .LT.  -1.91))) THEN
         CALL ENK5(191,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.23) .AND. (L .LT. 332.38)) .AND.
     &   ((B .GE.  -2.25) .AND. (B .LT.   0.29*L  -98.25))) THEN
         CALL ENK5(191,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.25) .AND. (L .LT. 332.38)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -2.20))) THEN
         CALL ENK5(191,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 333.02) .AND. (L .LT. 335.00)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -0.93))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 332.03) .AND. (L .LT. 333.01)) .AND.
     &   ((B .GE.  -1.91) .AND. (B .LT.  -0.93))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.45) .AND. (L .LT. 332.03)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.   0.55*L  -183.97))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.42) .AND. (L .LT. 332.05)) .AND.
     &   ((B .GE.  -1.91) .AND. (B .LT.  -1.58))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.42) .AND. (L .LT. 332.38)) .AND.
     &   ((B .LT.  -1.91) .AND. (B .GE.   0.29*L  -98.25))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 330.00) .AND. (L .LT. 335.00)) .AND.
     &   ((B .GE.  -7.73) .AND. (B .LT.  -2.42))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 330.00) .AND. (L .LT. 331.25)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -2.13))) THEN
         CALL ENK5(192,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 333.87) .AND. (L .LT. 335.01)) .AND.
     &   ((B .GE.   1.84) .AND. (B .LT.   3.47))) THEN
         CALL ENK5(193,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 333.66) .AND. (L .LT. 335.03)) .AND.
     &   ((B .LT.   1.84) .AND. (B .GE.  -0.86*L + 288.77))) THEN
         CALL ENK5(193,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 332.19) .AND. (L .LT. 333.87)) .AND.
     &   ((B .GE.   1.84) .AND. (B .LT.   0.79*L  -260.30))) THEN
         CALL ENK5(193,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 332.18) .AND. (L .LT. 333.66)) .AND.
     &   ((B .LT.   1.84) .AND. (B .GE.   0.79*L  -261.70))) THEN
         CALL ENK5(193,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.03) .AND. (L .LT. 332.19)) .AND.
     &   ((B .GE.   0.65) .AND. (B .LT.   0.79*L  -260.30))) THEN
         CALL ENK5(193,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.03) .AND. (L .LT. 332.18)) .AND.
     &   ((B .GE.   0.29) .AND. (B .LT.   0.65))) THEN
         CALL ENK5(193,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.04) .AND. (L .LT. 332.07)) .AND.
     &   ((B .GE.   1.75) .AND. (B .LT.  -3.00*L + 997.77))) THEN
         CALL ENK5(194,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 331.03) .AND. (L .LT. 332.07)) .AND.
     &   ((B .LT.   1.75) .AND. (B .GE.   0.79*L  -260.30))) THEN
         CALL ENK5(194,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 330.00) .AND. (L .LT. 331.03)) .AND.
     &   ((B .GE.   0.65) .AND. (B .LT.   4.83))) THEN
         CALL ENK5(194,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.93) .AND. (L .LT. 345.95)) .AND.
     &   ((B .GE.  -7.68) .AND. (B .LT.  -2.16*L + 739.60))) THEN
         CALL ENK5(195,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.00) .AND. (L .LT. 344.91)) .AND.
     &   ((B .GE.  -7.68) .AND. (B .LT.  -5.41))) THEN
         CALL ENK5(195,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.19) .AND. (L .LT. 340.00)) .AND.
     &   ((B .GE.  -7.68) .AND. (B .LT.  -6.03))) THEN
         CALL ENK5(195,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.90) .AND. (L .LT. 339.98)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -3.25))) THEN
         CALL ENK5(196,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.19) .AND. (L .LT. 336.90)) .AND.
     &   ((B .GE.  -6.03) .AND. (B .LT.  -3.09))) THEN
         CALL ENK5(196,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.61) .AND. (L .LT. 336.90)) .AND.
     &   ((B .GE.  -2.03) .AND. (B .LT.   0.13*L  -45.60))) THEN
         CALL ENK5(197,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.00) .AND. (L .LT. 336.90)) .AND.
     &   ((B .GE.  -3.09) .AND. (B .LT.  -2.03))) THEN
         CALL ENK5(197,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.77) .AND. (L .LT. 337.39)) .AND.
     &   ((B .GE.  -1.58) .AND. (B .LT.  -0.42*L + 140.13))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.77) .AND. (L .LT. 337.39)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.  -1.58))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.21) .AND. (L .LT. 336.77)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.   0.82*L  -277.40))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.61) .AND. (L .LT. 337.39)) .AND.
     &   ((B .LT.  -1.77) .AND. (B .GE.   0.13*L  -45.60))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.90) .AND. (L .LT. 336.24)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.  -1.36))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.59) .AND. (L .LT. 335.90)) .AND.
     &   ((B .GE.  -1.56) .AND. (B .LT.   0.61*L  -206.25))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.59) .AND. (L .LT. 335.90)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.  -1.56))) THEN
         CALL ENK5(198,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.68) .AND. (L .LT. 338.01)) .AND.
     &   ((B .GE.  -2.27) .AND. (B .LT.   0.46*L  -157.55))) THEN
         CALL ENK5(199,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.36) .AND. (L .LT. 338.01)) .AND.
     &   ((B .GE.  -2.49) .AND. (B .LT.  -2.27))) THEN
         CALL ENK5(199,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.40) .AND. (L .LT. 337.68)) .AND.
     &   ((B .GE.  -2.27) .AND. (B .LT.  -1.76*L + 591.98))) THEN
         CALL ENK5(199,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.89) .AND. (L .LT. 337.40)) .AND.
     &   ((B .GE.  -1.91) .AND. (B .LT.   0.13*L  -45.60))) THEN
         CALL ENK5(199,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.89) .AND. (L .LT. 337.40)) .AND.
     &   ((B .GE.  -2.46) .AND. (B .LT.  -1.91))) THEN
         CALL ENK5(199,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.21) .AND. (L .LT. 338.41)) .AND.
     &   ((B .LT.  -0.48) .AND. (B .GE.   5.29*L  -1790.95))) THEN
         CALL ENK5(200,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.08) .AND. (L .LT. 338.21)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.  -0.48))) THEN
         CALL ENK5(200,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.83) .AND. (L .LT. 338.07)) .AND.
     &   ((B .GE.  -1.29) .AND. (B .LT.  -0.48))) THEN
         CALL ENK5(200,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.68) .AND. (L .LT. 338.08)) .AND.
     &   ((B .LT.  -1.29) .AND. (B .GE.  -1.00*L + 336.48))) THEN
         CALL ENK5(200,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.68) .AND. (L .LT. 337.83)) .AND.
     &   ((B .GE.  -1.29) .AND. (B .LT.   5.42*L  -1831.42))) THEN
         CALL ENK5(200,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.97) .AND. (L .LT. 339.21)) .AND.
     &   ((B .LT.  -0.05) .AND. (B .GE.   1.76*L  -597.00))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.67) .AND. (L .LT. 338.97)) .AND.
     &   ((B .GE.  -0.48) .AND. (B .LT.  -0.05))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.68) .AND. (L .LT. 337.83)) .AND.
     &   ((B .LT.  -0.48) .AND. (B .GE.   5.42*L  -1831.42))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.75) .AND. (L .LT. 337.67)) .AND.
     &   ((B .GE.  -1.27) .AND. (B .LT.   1.30*L  -439.01))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.77) .AND. (L .LT. 337.39)) .AND.
     &   ((B .LT.  -1.27) .AND. (B .GE.  -0.42*L + 140.13))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.35) .AND. (L .LT. 337.68)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.  -1.27))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.68) .AND. (L .LT. 338.08)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.  -1.00*L + 336.48))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.40) .AND. (L .LT. 337.68)) .AND.
     &   ((B .LT.  -1.77) .AND. (B .GE.  -1.76*L + 591.98))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.68) .AND. (L .LT. 338.98)) .AND.
     &   ((B .LT.  -1.67) .AND. (B .GE.   0.46*L  -157.55))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.21) .AND. (L .LT. 338.41)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.   5.29*L  -1790.95))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.41) .AND. (L .LT. 338.98)) .AND.
     &   ((B .GE.  -1.67) .AND. (B .LT.  -0.48))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.97) .AND. (L .LT. 339.32)) .AND.
     &   ((B .GE.  -1.44) .AND. (B .LT.  -2.63*L + 891.05))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.98) .AND. (L .LT. 339.32)) .AND.
     &   ((B .LT.  -1.44) .AND. (B .GE.   0.46*L  -157.55))) THEN
         CALL ENK6(201,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.04) .AND. (L .LT. 337.67)) .AND.
     &   ((B .GE.  -0.05) .AND. (B .LT.  -0.76*L + 256.70))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.75) .AND. (L .LT. 337.67)) .AND.
     &   ((B .LT.  -0.05) .AND. (B .GE.   1.30*L  -439.01))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.02) .AND. (L .LT. 336.75)) .AND.
     &   ((B .GE.  -0.93) .AND. (B .LT.  -0.05))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.33) .AND. (L .LT. 336.74)) .AND.
     &   ((B .GE.  -1.42) .AND. (B .LT.  -0.93))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.02) .AND. (L .LT. 335.33)) .AND.
     &   ((B .LT.  -0.93) .AND. (B .GE.  -0.92*L + 307.35))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.30) .AND. (L .LT. 335.90)) .AND.
     &   ((B .LT.  -1.36) .AND. (B .GE.   0.61*L  -206.25))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.21) .AND. (L .LT. 336.77)) .AND.
     &   ((B .LT.  -1.32) .AND. (B .GE.   0.82*L  -277.40))) THEN
         CALL ENK6(202,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.37) .AND. (L .LT. 340.06)) .AND.
     &   ((B .GE.  -2.13) .AND. (B .LT.  -2.63*L + 891.05))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.99) .AND. (L .LT. 339.37)) .AND.
     &   ((B .GE.  -2.13) .AND. (B .LT.   0.46*L  -157.55))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.98) .AND. (L .LT. 341.40)) .AND.
     &   ((B .GE.  -3.06) .AND. (B .LT.  -2.13))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.77) .AND. (L .LT. 341.38)) .AND.
     &   ((B .LT.  -3.06) .AND. (B .GE.   1.69*L  -580.00))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.98) .AND. (L .LT. 340.77)) .AND.
     &   ((B .GE.  -4.09) .AND. (B .LT.  -3.06))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.02) .AND. (L .LT. 339.98)) .AND.
     &   ((B .GE.  -3.28) .AND. (B .LT.  -2.13))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 336.88) .AND. (L .LT. 338.02)) .AND.
     &   ((B .GE.  -3.28) .AND. (B .LT.  -2.49))) THEN
         CALL ENK6(203,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.02) .AND. (L .LT. 338.62)) .AND.
     &   ((B .GE.   1.29) .AND. (B .LT.   2.89))) THEN
         CALL ENK6(204,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.03) .AND. (L .LT. 338.63)) .AND.
     &   ((B .GE.   2.89) .AND. (B .LT.   4.93))) THEN
         CALL ENK6(205,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.04) .AND. (L .LT. 337.01)) .AND.
     &   ((B .GE.   1.87) .AND. (B .LT.   2.89))) THEN
         CALL ENK6(205,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 335.04) .AND. (L .LT. 337.02)) .AND.
     &   ((B .LT.   1.87) .AND. (B .GE.  -0.76*L + 256.70))) THEN
         CALL ENK6(205,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.03) .AND. (L .LT. 337.67)) .AND.
     &   ((B .GE.   0.36) .AND. (B .LT.   1.29))) THEN
         CALL ENK6(205,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.02) .AND. (L .LT. 337.67)) .AND.
     &   ((B .LT.   0.36) .AND. (B .GE.  -0.76*L + 256.70))) THEN
         CALL ENK6(205,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 337.61) .AND. (L .LT. 338.61)) .AND.
     &   ((B .GE.  -0.05) .AND. (B .LT.   1.29))) THEN
         CALL ENK6(205,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.50) .AND. (L .LT. 340.05)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.  -1.04*L + 353.60))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.21) .AND. (L .LT. 339.50)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.53))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.61) .AND. (L .LT. 339.21)) .AND.
     &   ((B .LT.   0.50) .AND. (B .GE.  -1.04*L + 352.70))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.97) .AND. (L .LT. 339.21)) .AND.
     &   ((B .GE.  -0.57) .AND. (B .LT.   1.76*L  -597.00))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.21) .AND. (L .LT. 339.83)) .AND.
     &   ((B .GE.  -0.60) .AND. (B .LT.  -0.07))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.84) .AND. (L .LT. 342.02)) .AND.
     &   ((B .GE.  -0.60) .AND. (B .LT.  -0.02))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.83) .AND. (L .LT. 340.05)) .AND.
     &   ((B .LT.  -0.60) .AND. (B .GE.  -1.34*L + 454.74))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.05) .AND. (L .LT. 340.34)) .AND.
     &   ((B .GE.  -0.93) .AND. (B .LT.  -0.60))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.34) .AND. (L .LT. 341.04)) .AND.
     &   ((B .LT.  -0.60) .AND. (B .GE.   0.48*L  -164.27))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.04) .AND. (L .LT. 342.02)) .AND.
     &   ((B .LT.  -0.60) .AND. (B .GE.  -0.71*L + 241.55))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.03) .AND. (L .LT. 342.98)) .AND.
     &   ((B .GE.  -0.81) .AND. (B .LT.  -0.85*L + 290.71))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.02) .AND. (L .LT. 344.98)) .AND.
     &   ((B .GE.  -1.29) .AND. (B .LT.  -0.81))) THEN
         CALL ENK6(206,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.35) .AND. (L .LT. 343.68)) .AND.
     &   ((B .GE.  -4.64) .AND. (B .LT.  -3.85))) THEN
         CALL ENK6(207,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.57) .AND. (L .LT. 343.35)) .AND.
     &   ((B .GE.  -4.64) .AND. (B .LT.   0.99*L  -343.80))) THEN
         CALL ENK6(207,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.00) .AND. (L .LT. 342.57)) .AND.
     &   ((B .GE.  -4.64) .AND. (B .LT.  -0.67*L + 224.85))) THEN
         CALL ENK6(207,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.39) .AND. (L .LT. 341.99)) .AND.
     &   ((B .GE.  -5.45) .AND. (B .LT.  -4.26))) THEN
         CALL ENK6(207,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.99) .AND. (L .LT. 343.67)) .AND.
     &   ((B .GE.  -5.45) .AND. (B .LT.  -4.64))) THEN
         CALL ENK6(207,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.68) .AND. (L .LT. 343.99)) .AND.
     &   ((B .GE.  -5.43) .AND. (B .LT.  -3.53*L + 1208.94))) THEN
         CALL ENK6(207,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.69) .AND. (L .LT. 344.96)) .AND.
     &   ((B .GE.  -4.38) .AND. (B .LT.  -3.04))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.99) .AND. (L .LT. 344.96)) .AND.
     &   ((B .GE.  -5.43) .AND. (B .LT.  -4.38))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.68) .AND. (L .LT. 343.99)) .AND.
     &   ((B .LT.  -4.38) .AND. (B .GE.  -3.53*L + 1208.94))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.58) .AND. (L .LT. 343.68)) .AND.
     &   ((B .GE.  -3.85) .AND. (B .LT.  -3.04))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.38) .AND. (L .LT. 342.58)) .AND.
     &   ((B .GE.  -4.26) .AND. (B .LT.  -3.04))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.57) .AND. (L .LT. 343.35)) .AND.
     &   ((B .LT.  -3.85) .AND. (B .GE.   0.99*L  -343.80))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.00) .AND. (L .LT. 342.57)) .AND.
     &   ((B .LT.  -4.26) .AND. (B .GE.  -0.67*L + 224.85))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.77) .AND. (L .LT. 341.38)) .AND.
     &   ((B .GE.  -4.09) .AND. (B .LT.   1.69*L  -580.00))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.00) .AND. (L .LT. 341.36)) .AND.
     &   ((B .GE.  -5.43) .AND. (B .LT.  -4.09))) THEN
         CALL ENK6(208,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.41) .AND. (L .LT. 344.97)) .AND.
     &   ((B .GE.  -3.06) .AND. (B .LT.  -1.65))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.62) .AND. (L .LT. 342.40)) .AND.
     &   ((B .GE.  -1.65) .AND. (B .LT.  -0.38*L + 128.42))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.07) .AND. (L .LT. 341.41)) .AND.
     &   ((B .GE.  -2.13) .AND. (B .LT.  -1.65))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.07) .AND. (L .LT. 340.62)) .AND.
     &   ((B .GE.  -1.65) .AND. (B .LT.  -0.98))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.06) .AND. (L .LT. 340.06)) .AND.
     &   ((B .LT.  -0.98) .AND. (B .GE.  -2.63*L + 890.85))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.99) .AND. (L .LT. 339.08)) .AND.
     &   ((B .LT.  -0.57) .AND. (B .GE.  -4.12*L + 1394.40))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.06) .AND. (L .LT. 339.83)) .AND.
     &   ((B .GE.  -0.98) .AND. (B .LT.  -0.53*L + 179.18))) THEN
         CALL ENK6(209,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.40) .AND. (L .LT. 344.98)) .AND.
     &   ((B .GE.  -1.65) .AND. (B .LT.  -1.24))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.51) .AND. (L .LT. 342.40)) .AND.
     &   ((B .LT.  -1.27) .AND. (B .GE.  -0.38*L + 128.42))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.62) .AND. (L .LT. 341.51)) .AND.
     &   ((B .LT.  -1.00) .AND. (B .GE.  -0.38*L + 128.42))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.52) .AND. (L .LT. 342.02)) .AND.
     &   ((B .GE.  -1.29) .AND. (B .LT.  -0.71*L + 241.55))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.04) .AND. (L .LT. 341.52)) .AND.
     &   ((B .GE.  -1.00) .AND. (B .LT.  -0.71*L + 241.55))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.34) .AND. (L .LT. 341.04)) .AND.
     &   ((B .GE.  -1.00) .AND. (B .LT.   0.48*L  -164.27))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.83) .AND. (L .LT. 340.05)) .AND.
     &   ((B .GE.  -0.86) .AND. (B .LT.  -1.34*L + 454.74))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.05) .AND. (L .LT. 340.33)) .AND.
     &   ((B .GE.  -0.98) .AND. (B .LT.  -0.86))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 339.06) .AND. (L .LT. 339.83)) .AND.
     &   ((B .LT.  -0.60) .AND. (B .GE.  -0.53*L + 179.18))) THEN
         CALL ENK6(210,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.81) .AND. (L .LT. 344.03)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   1.20))) THEN
         CALL ENK6(211,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.03) .AND. (L .LT. 344.24)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.  -5.98*L + 2058.22))) THEN
         CALL ENK6(211,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.24) .AND. (L .LT. 344.03)) .AND.
     &   ((B .GE.   1.20) .AND. (B .LT.   1.53))) THEN
         CALL ENK6(212,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.81) .AND. (L .LT. 343.24)) .AND.
     &   ((B .GE.   1.20) .AND. (B .LT.   0.77*L  -262.76))) THEN
         CALL ENK6(212,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.29) .AND. (L .LT. 345.00)) .AND.
     &   ((B .GE.   1.91) .AND. (B .LT.   3.73))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.79) .AND. (L .LT. 344.29)) .AND.
     &   ((B .GE.   1.75) .AND. (B .LT.   3.73))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.89) .AND. (L .LT. 342.79)) .AND.
     &   ((B .LT.   3.73) .AND. (B .GE.  -1.96*L + 673.74))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.24) .AND. (L .LT. 344.28)) .AND.
     &   ((B .GE.   1.51) .AND. (B .LT.   1.77))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 342.81) .AND. (L .LT. 343.24)) .AND.
     &   ((B .LT.   1.77) .AND. (B .GE.  -0.55*L + 190.25))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.03) .AND. (L .LT. 344.28)) .AND.
     &   ((B .GE.   1.17) .AND. (B .LT.   1.51))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 343.24) .AND. (L .LT. 344.28)) .AND.
     &   ((B .GE.   1.53) .AND. (B .LT.   1.77))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.24) .AND. (L .LT. 345.18)) .AND.
     &   ((B .GE.   0.00) .AND. (B .LT.   1.91))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.03) .AND. (L .LT. 344.24)) .AND.
     &   ((B .LT.   1.17) .AND. (B .GE.  -5.98*L + 2058.22))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.19) .AND. (L .LT. 345.17)) .AND.
     &   ((B .LT.   0.00) .AND. (B .GE.  -0.99*L + 340.75))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.17) .AND. (L .LT. 347.15)) .AND.
     &   ((B .LT.  -0.98) .AND. (B .GE.  -0.99*L + 340.75))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.14) .AND. (L .LT. 347.66)) .AND.
     &   ((B .LT.  -2.70) .AND. (B .GE.  -0.99*L + 340.75))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.19) .AND. (L .LT. 347.15)) .AND.
     &   ((B .GE.  -0.98) .AND. (B .LT.  -1.00*L + 346.30))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.15) .AND. (L .LT. 347.67)) .AND.
     &   ((B .GE.  -1.48) .AND. (B .LT.  -1.00*L + 346.30))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.67) .AND. (L .LT. 348.20)) .AND.
     &   ((B .GE.  -1.87) .AND. (B .LT.  -1.00*L + 346.30))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.14) .AND. (L .LT. 347.67)) .AND.
     &   ((B .GE.  -2.70) .AND. (B .LT.  -1.48))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.67) .AND. (L .LT. 348.26)) .AND.
     &   ((B .GE.  -3.41) .AND. (B .LT.  -1.87))) THEN
         CALL ENK6(213,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 341.19) .AND. (L .LT. 342.79)) .AND.
     &   ((B .GE.   1.72) .AND. (B .LT.  -1.96*L + 673.74))) THEN
         CALL ENK6(214,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 340.04) .AND. (L .LT. 341.17)) .AND.
     &   ((B .GE.   1.72) .AND. (B .LT.   4.97))) THEN
         CALL ENK6(214,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 338.60) .AND. (L .LT. 340.02)) .AND.
     &   ((B .GE.   1.60) .AND. (B .LT.   2.49))) THEN
         CALL ENK6(214,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 346.59) .AND. (L .LT. 347.39)) .AND.
     &   ((B .GE.   2.58) .AND. (B .LT.   4.98))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.39) .AND. (L .LT. 350.00)) .AND.
     &   ((B .GE.   1.91) .AND. (B .LT.   4.98))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 346.06) .AND. (L .LT. 346.59)) .AND.
     &   ((B .GE.   1.91) .AND. (B .LT.   4.98))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.22) .AND. (L .LT. 346.06)) .AND.
     &   ((B .GE.   1.91) .AND. (B .LT.   3.61*L  -1244.28))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.20) .AND. (L .LT. 346.56)) .AND.
     &   ((B .GE.   0.98) .AND. (B .LT.   1.91))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 346.03) .AND. (L .LT. 346.55)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.98))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.19) .AND. (L .LT. 346.02)) .AND.
     &   ((B .LT.   0.98) .AND. (B .GE.  -1.00*L + 346.30))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 346.02) .AND. (L .LT. 347.48)) .AND.
     &   ((B .LT.  -0.07) .AND. (B .GE.  -1.00*L + 346.30))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.49) .AND. (L .LT. 350.59)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.  -0.07))) THEN
         CALL ENK6(215,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.00) .AND. (L .LT. 347.66)) .AND.
     &   ((B .GE.  -3.21) .AND. (B .LT.  -0.99*L + 340.75))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.97) .AND. (L .LT. 348.26)) .AND.
     &   ((B .GE.  -3.73) .AND. (B .LT.  -3.21))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 346.40) .AND. (L .LT. 346.96)) .AND.
     &   ((B .LT.  -3.68) .AND. (B .GE.   0.82*L  -288.10))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.96) .AND. (L .LT. 346.40)) .AND.
     &   ((B .GE.  -4.14) .AND. (B .LT.  -3.68))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.25) .AND. (L .LT. 348.98)) .AND.
     &   ((B .GE.  -3.73) .AND. (B .LT.  -2.42))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.42) .AND. (L .LT. 348.99)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -1.51*L + 524.35))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.23) .AND. (L .LT. 348.41)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -1.91))) THEN
         CALL ENK6(216,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.99) .AND. (L .LT. 350.01)) .AND.
     &   ((B .LT.  -1.24) .AND. (B .GE.   1.15*L  -403.75))) THEN
         CALL ENK6(217,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.23) .AND. (L .LT. 348.99)) .AND.
     &   ((B .LT.  -1.24) .AND. (B .GE.  -1.51*L + 524.35))) THEN
         CALL ENK6(217,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 346.43) .AND. (L .LT. 347.96)) .AND.
     &   ((B .GE.  -5.02) .AND. (B .LT.  -0.58*L + 196.76))) THEN
         CALL ENK6(218,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.96) .AND. (L .LT. 346.40)) .AND.
     &   ((B .GE.  -5.02) .AND. (B .LT.  -4.14))) THEN
         CALL ENK6(218,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.95) .AND. (L .LT. 345.96)) .AND.
     &   ((B .GE.  -5.43) .AND. (B .LT.  -5.02))) THEN
         CALL ENK6(218,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 344.93) .AND. (L .LT. 345.96)) .AND.
     &   ((B .LT.  -5.41) .AND. (B .GE.  -2.16*L + 739.60))) THEN
         CALL ENK6(218,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 345.95) .AND. (L .LT. 347.96)) .AND.
     &   ((B .GE.  -7.68) .AND. (B .LT.  -5.02))) THEN
         CALL ENK6(218,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.96) .AND. (L .LT. 349.98)) .AND.
     &   ((B .GE.  -4.98) .AND. (B .LT.  -3.71))) THEN
         CALL ENK6(220,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 347.98) .AND. (L .LT. 349.26)) .AND.
     &   ((B .GE.  -7.68) .AND. (B .LT.  -4.98))) THEN
         CALL ENK6(220,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 349.26) .AND. (L .LT. 352.96)) .AND.
     &   ((B .GE.  -6.99) .AND. (B .LT.  -4.98))) THEN
         CALL ENK6(220,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.59) .AND. (L .LT. 351.01)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.  -3.45*L + 1208.29))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.01) .AND. (L .LT. 350.58)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.  -1.22))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.99) .AND. (L .LT. 350.01)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.   1.15*L  -403.75))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.99) .AND. (L .LT. 350.01)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.  -2.42))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 348.99) .AND. (L .LT. 349.99)) .AND.
     &   ((B .GE.  -3.73) .AND. (B .LT.  -2.68))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 353.45) .AND. (L .LT. 354.40)) .AND.
     &   ((B .GE.  -4.00) .AND. (B .LT.  -2.39))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.01) .AND. (L .LT. 353.45)) .AND.
     &   ((B .GE.  -2.68) .AND. (B .LT.   0.06*L  -23.70))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.00) .AND. (L .LT. 353.45)) .AND.
     &   ((B .GE.  -4.00) .AND. (B .LT.  -2.68))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 349.99) .AND. (L .LT. 352.99)) .AND.
     &   ((B .GE.  -4.98) .AND. (B .LT.  -3.97))) THEN
         CALL ENK6(221,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 353.40) .AND. (L .LT. 355.10)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.96))) THEN
         CALL ENK6(222,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.03) .AND. (L .LT. 355.81)) .AND.
     &   ((B .GE.  -1.77) .AND. (B .LT.  -0.07))) THEN
         CALL ENK6(222,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.02) .AND. (L .LT. 354.99)) .AND.
     &   ((B .LT.  -1.77) .AND. (B .GE.   0.23*L  -83.40))) THEN
         CALL ENK6(222,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.60) .AND. (L .LT. 351.02)) .AND.
     &   ((B .GE.  -1.22) .AND. (B .LT.  -0.07))) THEN
         CALL ENK6(222,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.59) .AND. (L .LT. 351.01)) .AND.
     &   ((B .LT.  -1.22) .AND. (B .GE.  -3.45*L + 1208.29))) THEN
         CALL ENK6(222,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.00) .AND. (L .LT. 351.73)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   1.27))) THEN
         CALL ENK6(223,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.72) .AND. (L .LT. 353.42)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   1.70))) THEN
         CALL ENK6(224,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.72) .AND. (L .LT. 353.42)) .AND.
     &   ((B .GE.   1.65) .AND. (B .LT.   2.18))) THEN
         CALL ENK6(225,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.00) .AND. (L .LT. 351.72)) .AND.
     &   ((B .GE.   1.29) .AND. (B .LT.   2.18))) THEN
         CALL ENK6(225,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 352.05) .AND. (L .LT. 353.42)) .AND.
     &   ((B .GE.   2.97) .AND. (B .LT.   3.76))) THEN
         CALL ENK6(226,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 353.42) .AND. (L .LT. 356.02)) .AND.
     &   ((B .GE.   1.00) .AND. (B .LT.   4.98))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 355.10) .AND. (L .LT. 356.04)) .AND.
     &   ((B .GE.  -0.07) .AND. (B .LT.   0.98))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 353.04) .AND. (L .LT. 353.42)) .AND.
     &   ((B .GE.   3.78) .AND. (B .LT.   5.00))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 352.07) .AND. (L .LT. 353.43)) .AND.
     &   ((B .GE.   2.18) .AND. (B .LT.   2.97))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.43) .AND. (L .LT. 352.06)) .AND.
     &   ((B .LT.   3.54) .AND. (B .GE.  -0.85*L + 301.50))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.46) .AND. (L .LT. 351.79)) .AND.
     &   ((B .GE.   3.54) .AND. (B .LT.   0.47*L  -161.20))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.78) .AND. (L .LT. 352.05)) .AND.
     &   ((B .GE.   3.78) .AND. (B .LT.  -1.40*L + 496.70))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 351.78) .AND. (L .LT. 352.05)) .AND.
     &   ((B .GE.   3.54) .AND. (B .LT.   3.78))) THEN
         CALL ENK6(227,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 350.00) .AND. (L .LT. 353.02)) .AND.
     &   ((B .GE.   4.95) .AND. (B .LT.   7.58))) THEN
         CALL ENK6(228,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.03) .AND. (L .LT. 358.38)) .AND.
     &   ((B .GE.   0.02) .AND. (B .LT.  -0.59*L + 211.52))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 355.82) .AND. (L .LT. 358.36)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.   0.02))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 358.38) .AND. (L .LT. 359.64)) .AND.
     &   ((B .GE.  -1.39) .AND. (B .LT.  -1.11*L + 397.90))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 358.89) .AND. (L .LT. 359.64)) .AND.
     &   ((B .GE.  -2.01) .AND. (B .LT.  -1.39))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.91) .AND. (L .LT. 358.89)) .AND.
     &   ((B .LT.  -1.39) .AND. (B .GE.  -0.63*L + 224.11))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 355.81) .AND. (L .LT. 357.42)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -1.39))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.02) .AND. (L .LT. 357.41)) .AND.
     &   ((B .GE.  -3.18) .AND. (B .LT.  -2.42))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 354.99) .AND. (L .LT. 355.81)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -1.79))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 353.40) .AND. (L .LT. 354.99)) .AND.
     &   ((B .GE.  -2.15) .AND. (B .LT.   0.23*L  -83.40))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 353.40) .AND. (L .LT. 354.99)) .AND.
     &   ((B .GE.  -2.42) .AND. (B .LT.  -2.15))) THEN
         CALL ENK6(229,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 358.07) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.  -5.48) .AND. (B .LT.   0.53*L  -195.46))) THEN
         CALL ENK6(230,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 352.97) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.  -7.56) .AND. (B .LT.  -5.48))) THEN
         CALL ENK6(230,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.00) .AND. (L .LT. 357.04)) .AND.
     &   ((B .GE.  -5.48) .AND. (B .LT.  -1.47*L + 519.40))) THEN
         CALL ENK6(230,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 352.98) .AND. (L .LT. 356.00)) .AND.
     &   ((B .GE.  -5.48) .AND. (B .LT.  -3.97))) THEN
         CALL ENK6(230,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 359.62) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.  -1.36) .AND. (B .LT.   1.05*L  -378.97))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 359.64) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.  -1.96) .AND. (B .LT.  -1.36))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.90) .AND. (L .LT. 359.98)) .AND.
     &   ((B .GE.  -4.47) .AND. (B .LT.  -1.96))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 358.07) .AND. (L .LT. 360.00)) .AND.
     &   ((B .LT.  -4.47) .AND. (B .GE.   0.53*L  -195.46))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.91) .AND. (L .LT. 358.89)) .AND.
     &   ((B .GE.  -2.01) .AND. (B .LT.  -0.63*L + 225.31))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.45) .AND. (L .LT. 357.91)) .AND.
     &   ((B .GE.  -3.18) .AND. (B .LT.  -1.39))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.04) .AND. (L .LT. 358.07)) .AND.
     &   ((B .GE.  -5.48) .AND. (B .LT.  -3.18))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.04) .AND. (L .LT. 357.05)) .AND.
     &   ((B .GE.  -3.97) .AND. (B .LT.  -3.16))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.00) .AND. (L .LT. 357.04)) .AND.
     &   ((B .LT.  -3.97) .AND. (B .GE.  -1.47*L + 519.40))) THEN
         CALL ENK6(231,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.13) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.   0.81) .AND. (B .LT.   0.43*L  -151.23))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 358.38) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.   0.02) .AND. (B .LT.   0.81))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.11) .AND. (L .LT. 358.38)) .AND.
     &   ((B .LT.   0.81) .AND. (B .GE.  -0.59*L + 211.52))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 359.72) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.  -0.93) .AND. (B .LT.   0.02))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 359.62) .AND. (L .LT. 360.00)) .AND.
     &   ((B .LT.  -0.93) .AND. (B .GE.   1.05*L  -378.97))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 359.26) .AND. (L .LT. 359.64)) .AND.
     &   ((B .LT.  -0.93) .AND. (B .GE.  -1.11*L + 397.90))) THEN
         CALL ENK6(233,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.03) .AND. (L .LT. 360.00)) .AND.
     &   ((B .GE.   2.08) .AND. (B .LT.   7.63))) THEN
         CALL ENK6(235,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 357.13) .AND. (L .LT. 360.00)) .AND.
     &   ((B .LT.   2.06) .AND. (B .GE.   0.43*L  -151.23))) THEN
         CALL ENK6(235,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.03) .AND. (L .LT. 357.13)) .AND.
     &   ((B .GE.   1.41) .AND. (B .LT.   2.08))) THEN
         CALL ENK6(235,D,AV,A0)
         RETURN
      ENDIF
      IF (((L .GE. 356.03) .AND. (L .LT. 357.11)) .AND.
     &   ((B .LT.   1.41) .AND. (B .GE.  -0.59*L + 211.52))) THEN
         CALL ENK6(235,D,AV,A0)
         RETURN
      ENDIF
      END

C     ===============================================================
      SUBROUTINE ENK1(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL
      
      IF(CELL .EQ.   1) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.02)) AV= 1.808*D
        IF((D .GE. 1.02).AND.(D .LT. 5.27)) AV= 1.844
        IF(D .GE. 5.27) AV= A0*(D- 5.27)+ 1.844
      RETURN
      ENDIF
      IF(CELL .EQ.   2) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.28)) AV= 0.000
        IF((D .GE. 0.28).AND.(D .LT. 0.99)) AV= 1.846*(D- 0.28)
        IF((D .GE. 0.99).AND.(D .LT. 4.67)) AV= 1.311
        IF(D .GE. 4.67) AV= A0*(D- 4.67)+ 1.311
      RETURN
      ENDIF
      IF(CELL .EQ.   3) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 0.256*D
        IF((D .GE. 0.14).AND.(D .LT. 0.47)) AV= 2.770*(D- 0.14)+0.036
        IF((D .GE. 0.47).AND.(D .LT. 0.77)) AV= 0.950
        IF((D .GE. 0.77).AND.(D .LT. 1.02)) AV= 4.515*(D- 0.77)+0.950
        IF((D .GE. 1.02).AND.(D .LT. 4.37)) AV= 2.079
        IF(D .GE. 4.37) AV= A0*(D- 4.37)+ 2.079
      RETURN
      ENDIF
      IF(CELL .EQ.   4) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 0.000
        IF((D .GE. 0.15).AND.(D .LT. 0.39)) AV= 2.649*(D- 0.15)
        IF((D .GE. 0.39).AND.(D .LT. 0.86)) AV= 0.641*(D- 0.39)+0.636
        IF((D .GE. 0.86).AND.(D .LT. 1.03)) AV=12.989*(D- 0.86)+0.937
        IF((D .GE. 1.03).AND.(D .LT. 3.54)) AV= 3.145
        IF(D .GE. 3.54) AV= A0*(D- 3.54)+ 3.145
      RETURN
      ENDIF
      IF(CELL .EQ.   5) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.31)) AV= 2.454*D
        IF((D .GE. 0.31).AND.(D .LT. 0.89)) AV= 0.761
        IF((D .GE. 0.89).AND.(D .LT. 0.95)) AV= 8.943*(D- 0.89)+0.761
        IF((D .GE. 0.95).AND.(D .LT. 3.36)) AV= 1.298
        IF(D .GE. 3.36) AV= A0*(D- 3.36)+ 1.298
      RETURN
      ENDIF
      IF(CELL .EQ.   6) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 0.97)) AV= 0.825*(D- 0.11)
        IF((D .GE. 0.97).AND.(D .LT. 3.62)) AV= 0.710
        IF(D .GE. 3.62) AV= A0*(D- 3.62)+ 0.710
      RETURN
      ENDIF
      IF(CELL .EQ.   7) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 1.00)) AV= 2.808*(D- 0.11)
        IF((D .GE. 1.00).AND.(D .LT. 4.26)) AV= 2.499
        IF(D .GE. 4.26) AV= A0*(D- 4.26)+ 2.499
      RETURN
      ENDIF
      IF(CELL .EQ.   8) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 1.07)) AV= 1.436*(D- 0.11)
        IF((D .GE. 1.07).AND.(D .LT. 2.72)) AV= 1.379
        IF((D .GE. 2.72).AND.(D .LT. 2.93)) AV= 2.564*(D- 2.72)+1.379
        IF((D .GE. 2.93).AND.(D .LT. 5.72)) AV= 1.917
        IF(D .GE. 5.72) AV= A0*(D- 5.72)+ 1.917
      RETURN
      ENDIF
      IF(CELL .EQ.   9) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.00)) AV= 3.524*D
        IF((D .GE. 1.00).AND.(D .LT. 3.85)) AV= 3.524
        IF(D .GE. 3.85) AV= A0*(D- 3.85)+ 3.524
      RETURN
      ENDIF
      IF(CELL .EQ.  10) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 0.067*D
        IF((D .GE. 0.16).AND.(D .LT. 1.08)) AV= 2.041*(D- 0.16)+0.011
        IF((D .GE. 1.08).AND.(D .LT. 4.70)) AV= 1.889
        IF(D .GE. 4.70) AV= A0*(D- 4.70)+ 1.889
      RETURN
      ENDIF
      IF(CELL .EQ.  11) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 0.50)) AV= 1.803*(D- 0.11)
        IF((D .GE. 0.50).AND.(D .LT. 1.00)) AV= 4.538*(D- 0.50)+0.703
        IF((D .GE. 1.00).AND.(D .LT. 3.79)) AV= 2.972
        IF(D .GE. 3.79) AV= A0*(D- 3.79)+ 2.972
      RETURN
      ENDIF
                IF(CELL .EQ.  12) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 0.000
        IF((D .GE. 0.10).AND.(D .LT. 1.00)) AV= 2.186*(D- 0.10)
        IF((D .GE. 1.00).AND.(D .LT. 4.57)) AV= 1.967
        IF(D .GE. 4.57) AV= A0*(D- 4.57)+ 1.967
      RETURN
      ENDIF
      IF(CELL .EQ.  13) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 0.000
        IF((D .GE. 0.12).AND.(D .LT. 1.09)) AV= 3.045*(D- 0.12)
        IF((D .GE. 1.09).AND.(D .LT. 3.94)) AV= 2.954
        IF(D .GE. 3.94) AV= A0*(D- 3.94)+ 2.954
      RETURN
      ENDIF
      IF(CELL .EQ.  14) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.23)) AV= 0.557*D
        IF((D .GE. 0.23).AND.(D .LT. 1.04)) AV= 2.923*(D- 0.23)+0.128
        IF((D .GE. 1.04).AND.(D .LT. 4.84)) AV= 2.496
        IF(D .GE. 4.84) AV= A0*(D- 4.84)+ 2.496
      RETURN
      ENDIF
      IF(CELL .EQ.  15) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 0.000
        IF((D .GE. 0.10).AND.(D .LT. 1.01)) AV= 3.416*(D- 0.10)
        IF((D .GE. 1.01).AND.(D .LT. 5.81)) AV= 3.109
        IF(D .GE. 5.81) AV= A0*(D- 5.81)+ 3.109
      RETURN
      ENDIF
      IF(CELL .EQ.  16) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 0.000
        IF((D .GE. 0.14).AND.(D .LT. 1.04)) AV= 3.007*(D- 0.14)
        IF((D .GE. 1.04).AND.(D .LT. 5.53)) AV= 2.706
        IF(D .GE. 5.53) AV= A0*(D- 5.53)+ 2.706
      RETURN
      ENDIF
      IF(CELL .EQ.  17) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 0.165*D
        IF((D .GE. 0.07).AND.(D .LT. 0.19)) AV= 3.715*(D- 0.07)+0.012
        IF((D .GE. 0.19).AND.(D .LT. 1.00)) AV= 1.359*(D- 0.19)+0.458
        IF((D .GE. 1.00).AND.(D .LT. 4.37)) AV= 1.559
        IF(D .GE. 4.37) AV= A0*(D- 4.37)+ 1.559
      RETURN
      ENDIF
      IF(CELL .EQ.  18) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 0.000
        IF((D .GE. 0.09).AND.(D .LT. 0.99)) AV= 2.773*(D- 0.09)
        IF((D .GE. 0.99).AND.(D .LT. 4.80)) AV= 2.496
        IF(D .GE. 4.80) AV= A0*(D- 4.80)+ 2.496
      RETURN
      ENDIF
      IF(CELL .EQ.  19) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.06)) AV= 0.199*D
        IF((D .GE. 0.06).AND.(D .LT. 0.98)) AV= 3.389*(D- 0.06)+0.012
        IF((D .GE. 0.98).AND.(D .LT. 4.59)) AV= 3.130
        IF(D .GE. 4.59) AV= A0*(D- 4.59)+ 3.130
      RETURN
      ENDIF
      IF(CELL .EQ.  20) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 0.249*D
        IF((D .GE. 0.09).AND.(D .LT. 0.99)) AV= 2.092*(D- 0.09)+0.022
        IF((D .GE. 0.99).AND.(D .LT. 4.98)) AV= 1.905
        IF(D .GE. 4.98) AV= A0*(D- 4.98)+ 1.905
      RETURN
      ENDIF
      IF(CELL .EQ.  21) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.02)) AV= 0.982*D
        IF((D .GE. 0.02).AND.(D .LT. 1.03)) AV= 1.075*(D- 0.02)+0.020
        IF((D .GE. 1.03).AND.(D .LT. 4.98)) AV= 1.106
        IF(D .GE. 4.98) AV= A0*(D- 4.98)+ 1.106
      RETURN
      ENDIF
      IF(CELL .EQ.  22) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.03)) AV= 1.736*D
        IF((D .GE. 1.03).AND.(D .LT. 3.90)) AV= 1.788
        IF(D .GE. 3.90) AV= A0*(D- 3.90)+ 1.788
      RETURN
      ENDIF
      IF(CELL .EQ.  23) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.02)) AV= 1.401*D
        IF((D .GE. 1.02).AND.(D .LT. 3.52)) AV= 1.429
        IF(D .GE. 3.52) AV= A0*(D- 3.52)+ 1.429
      RETURN
      ENDIF
      IF(CELL .EQ.  24) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 1.01)) AV= 2.313*(D- 0.11)
        IF((D .GE. 1.01).AND.(D .LT. 4.92)) AV= 2.082
        IF(D .GE. 4.92) AV= A0*(D- 4.92)+ 2.082
      RETURN
      ENDIF
      IF(CELL .EQ.  25) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.97)) AV= 2.001*D
        IF((D .GE. 0.97).AND.(D .LT. 1.67)) AV= 1.941
        IF((D .GE. 1.67).AND.(D .LT. 2.02)) AV= 1.706*(D- 1.67)+1.941
        IF((D .GE. 2.02).AND.(D .LT. 4.94)) AV= 2.538
        IF(D .GE. 4.94) AV= A0*(D- 4.94)+ 2.538
      RETURN
      ENDIF
      IF(CELL .EQ.  26) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.23)) AV= 4.769*D
        IF((D .GE. 0.23).AND.(D .LT. 1.02)) AV= 2.118*(D- 0.23)+1.097
        IF((D .GE. 1.02).AND.(D .LT. 1.76)) AV= 2.770
        IF((D .GE. 1.76).AND.(D .LT. 2.11)) AV= 3.840*(D- 1.76)+2.770
        IF(D .GE. 2.11) AV= A0*(D- 2.11)+ 4.114
      RETURN
      ENDIF
      IF(CELL .EQ.  27) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 5.145*D
        IF((D .GE. 0.14).AND.(D .LT. 0.88)) AV= 1.365*(D- 0.14)+0.720
        IF((D .GE. 0.88).AND.(D .LT. 1.02)) AV= 8.678*(D- 0.88)+1.730
        IF((D .GE. 1.02).AND.(D .LT. 2.01)) AV= 2.945
        IF(D .GE. 2.01) AV= A0*(D- 2.01)+ 2.945
      RETURN
      ENDIF
      IF(CELL .EQ.  28) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.21)) AV= 1.549*D
        IF((D .GE. 0.21).AND.(D .LT. 1.00)) AV= 2.844*(D- 0.21)+0.325
        IF((D .GE. 1.00).AND.(D .LT. 4.61)) AV= 2.572
        IF(D .GE. 4.61) AV= A0*(D- 4.61)+ 2.572
      RETURN
      ENDIF
      IF(CELL .EQ.  29) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.44)) AV= 3.414*D
        IF((D .GE. 0.44).AND.(D .LT. 1.03)) AV= 1.063*(D- 0.44)+1.502
        IF((D .GE. 1.03).AND.(D .LT. 4.35)) AV= 2.129
        IF(D .GE. 4.35) AV= A0*(D- 4.35)+ 2.129
      RETURN
      ENDIF
      IF(CELL .EQ.  30) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.02)) AV= 3.288*D
        IF((D .GE. 1.02).AND.(D .LT. 4.27)) AV= 3.354
        IF(D .GE. 4.27) AV= A0*(D- 4.27)+ 3.354
      RETURN
      ENDIF
      IF(CELL .EQ.  31) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 0.000
        IF((D .GE. 0.15).AND.(D .LT. 0.54)) AV= 5.354*(D- 0.15)
        IF((D .GE. 0.54).AND.(D .LT. 4.37)) AV= 2.088
        IF(D .GE. 4.37) AV= A0*(D- 4.37)+ 2.088
      RETURN
      ENDIF
      IF(CELL .EQ.  32) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.07)) AV= 2.130*D
        IF((D .GE. 1.07).AND.(D .LT. 2.43)) AV= 2.279
        IF((D .GE. 2.43).AND.(D .LT. 2.52)) AV=13.064*(D- 2.43)+2.279
        IF((D .GE. 2.52).AND.(D .LT. 4.71)) AV= 3.455
        IF(D .GE. 4.71) AV= A0*(D- 4.71)+ 3.455
      RETURN
      ENDIF
      IF(CELL .EQ.  33) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.70)) AV= 2.727*D
        IF((D .GE. 0.70).AND.(D .LT. 1.01)) AV= 5.163*(D- 0.70)+1.909
        IF((D .GE. 1.01).AND.(D .LT. 4.16)) AV= 3.510
        IF(D .GE. 4.16) AV= A0*(D- 4.16)+ 3.510
      RETURN
      ENDIF
      IF(CELL .EQ.  34) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.99)) AV= 2.364*D
        IF((D .GE. 0.99).AND.(D .LT. 3.92)) AV= 2.340
        IF(D .GE. 3.92) AV= A0*(D- 3.92)+ 2.340
      RETURN
      ENDIF
      IF(CELL .EQ.  35) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.64)) AV= 1.436*D
        IF(D .GE. 1.64) AV= A0*(D- 1.64)+ 2.355
      RETURN
      ENDIF
      IF(CELL .EQ.  36) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.72)) AV= 2.167*D
        IF(D .GE. 0.72) AV= A0*(D- 0.72)+ 1.560
      RETURN
      ENDIF
      IF(CELL .EQ.  37) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.67)) AV= 1.305*D
        IF((D .GE. 0.67).AND.(D .LT. 2.96)) AV= 0.874
        IF(D .GE. 2.96) AV= A0*(D- 2.96)+ 0.874
      RETURN
      ENDIF
      IF(CELL .EQ.  38) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.36)) AV= 0.845*D
        IF((D .GE. 0.36).AND.(D .LT. 1.67)) AV= 0.493*(D- 0.36)+0.304
        IF(D .GE. 1.67) AV= A0*(D- 1.67)+ 0.950
      RETURN
      ENDIF
      IF(CELL .EQ.  39) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 0.096*D
        IF((D .GE. 0.15).AND.(D .LT. 0.53)) AV= 2.651*(D- 0.15)+0.014
        IF(D .GE. 0.53) AV= A0*(D- 0.53)+ 1.021
      RETURN
      ENDIF
      IF(CELL .EQ.  40) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 0.000
        IF((D .GE. 0.16).AND.(D .LT. 0.97)) AV= 0.989*(D- 0.16)
        IF((D .GE. 0.97).AND.(D .LT. 2.75)) AV= 0.801
        IF((D .GE. 2.75).AND.(D .LT. 3.25)) AV= 3.608*(D- 2.75)+0.801
        IF((D .GE. 3.25).AND.(D .LT. 5.47)) AV= 2.605
        IF(D .GE. 5.47) AV= A0*(D- 5.47)+ 2.605
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK2(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ.  41) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.03)) AV= 1.219*D
        IF((D .GE. 1.03).AND.(D .LT. 2.53)) AV= 1.256
        IF(D .GE. 2.53) AV= A0*(D- 2.53)+ 1.256
      RETURN
      ENDIF
      IF(CELL .EQ.  42) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 5.473*D
        IF((D .GE. 0.05).AND.(D .LT. 0.82)) AV= 0.274
        IF((D .GE. 0.82).AND.(D .LT. 1.27)) AV= 0.969*(D- 0.82)+0.274
        IF((D .GE. 1.27).AND.(D .LT. 2.78)) AV= 0.710
        IF((D .GE. 2.78).AND.(D .LT. 3.33)) AV= 1.984*(D- 2.78)+0.710
        IF((D .GE. 3.33).AND.(D .LT. 5.57)) AV= 1.801
        IF(D .GE. 5.57) AV= A0*(D- 5.57)+ 1.801
      RETURN
      ENDIF
      IF(CELL .EQ.  43) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.39)) AV= 1.183*D
        IF((D .GE. 0.39).AND.(D .LT. 0.79)) AV= 0.799*(D- 0.39)+0.461
        IF(D .GE. 0.79) AV= A0*(D- 0.79)+ 0.781
      RETURN
      ENDIF
      IF(CELL .EQ.  44) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.82)) AV= 0.989*D
        IF((D .GE. 0.82).AND.(D .LT. 2.30)) AV= 0.811
        IF(D .GE. 2.30) AV= A0*(D- 2.30)+ 0.811
      RETURN
      ENDIF
      IF(CELL .EQ.  45) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.84)) AV= 1.652*D
        IF((D .GE. 0.84).AND.(D .LT. 2.20)) AV= 0.356*(D- 0.84)+1.388
        IF((D .GE. 2.20).AND.(D .LT. 4.57)) AV= 1.872
        IF(D .GE. 4.57) AV= A0*(D- 4.57)+ 1.872
      RETURN
      ENDIF
      IF(CELL .EQ.  46) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.98)) AV= 0.789*D
        IF((D .GE. 0.98).AND.(D .LT. 1.96)) AV= 0.452*(D- 0.98)+0.773
        IF(D .GE. 1.96) AV= A0*(D- 1.96)+ 1.216
      RETURN
      ENDIF
      IF(CELL .EQ.  47) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 0.19)) AV=18.698*(D- 0.11)
        IF((D .GE. 0.19).AND.(D .LT. 1.69)) AV= 0.457*(D- 0.19)+1.496
        IF(D .GE. 1.69) AV= A0*(D- 1.69)+ 2.182
      RETURN
      ENDIF
      IF(CELL .EQ.  48) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.06)) AV=18.527*D
        IF((D .GE. 0.06).AND.(D .LT. 3.87)) AV= 1.112
        IF(D .GE. 3.87) AV= A0*(D- 3.87)+ 1.112
      RETURN
      ENDIF
      IF(CELL .EQ.  49) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.03)) AV= 0.883*D
        IF((D .GE. 1.03).AND.(D .LT. 3.24)) AV= 0.909
        IF(D .GE. 3.24) AV= A0*(D- 3.24)+ 0.909
      RETURN
      ENDIF
      IF(CELL .EQ.  50) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.98)) AV= 1.314*(D- 0.08)
        IF((D .GE. 0.98).AND.(D .LT. 1.53)) AV= 0.464*(D- 0.98)+1.183
        IF((D .GE. 1.53).AND.(D .LT. 4.01)) AV= 1.438
        IF(D .GE. 4.01) AV= A0*(D- 4.01)+ 1.438
      RETURN
      ENDIF
      IF(CELL .EQ.  51) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.88)) AV= 0.000
        IF((D .GE. 0.88).AND.(D .LT. 1.23)) AV= 4.392*(D- 0.88)
        IF((D .GE. 1.23).AND.(D .LT. 4.95)) AV= 1.537
        IF(D .GE. 4.95) AV= A0*(D- 4.95)+ 1.537
      RETURN
      ENDIF
      IF(CELL .EQ.  52) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.00)) AV= 0.196*D
        IF((D .GE. 1.00).AND.(D .LT. 1.29)) AV= 8.103*(D- 1.00)+0.196
        IF((D .GE. 1.29).AND.(D .LT. 4.05)) AV= 2.546
        IF(D .GE. 4.05) AV= A0*(D- 4.05)+ 2.546
      RETURN
      ENDIF
      IF(CELL .EQ.  53) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.79)) AV= 0.304*D
        IF((D .GE. 0.79).AND.(D .LT. 1.02)) AV= 3.627*(D- 0.79)+0.240
        IF((D .GE. 1.02).AND.(D .LT. 3.84)) AV= 1.074
        IF(D .GE. 3.84) AV= A0*(D- 3.84)+ 1.074
      RETURN
      ENDIF
      IF(CELL .EQ.  54) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.99)) AV= 3.019*D
        IF((D .GE. 0.99).AND.(D .LT. 2.65)) AV= 2.989
        IF(D .GE. 2.65) AV= A0*(D- 2.65)+ 2.989
      RETURN
      ENDIF
      IF(CELL .EQ.  55) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 1.031*D
        IF((D .GE. 0.11).AND.(D .LT. 0.75)) AV= 0.113
        IF((D .GE. 0.75).AND.(D .LT. 1.00)) AV= 6.490*(D- 0.75)+0.113
        IF((D .GE. 1.00).AND.(D .LT. 3.35)) AV= 1.736
        IF(D .GE. 3.35) AV= A0*(D- 3.35)+ 1.736
      RETURN
      ENDIF
      IF(CELL .EQ.  56) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.83)) AV= 0.000
        IF((D .GE. 0.83).AND.(D .LT. 1.03)) AV= 9.452*(D- 0.83)
        IF((D .GE. 1.03).AND.(D .LT. 1.88)) AV= 1.890
        IF(D .GE. 1.88) AV= A0*(D- 1.88)+ 1.890
      RETURN
      ENDIF
      IF(CELL .EQ.  57) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.74)) AV= 0.000
        IF((D .GE. 0.74).AND.(D .LT. 1.01)) AV= 7.538*(D- 0.74)
        IF((D .GE. 1.01).AND.(D .LT. 3.55)) AV= 2.035
        IF(D .GE. 3.55) AV= A0*(D- 3.55)+ 2.035
      RETURN
      ENDIF
      IF(CELL .EQ.  58) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.36)) AV= 0.430*(D- 0.01)
        IF((D .GE. 0.36).AND.(D .LT. 0.92)) AV= 0.906*(D- 0.36)+0.155
        IF((D .GE. 0.92).AND.(D .LT. 1.00)) AV= 1.937*(D- 0.92)+0.662
        IF((D .GE. 1.00).AND.(D .LT. 2.03)) AV= 0.390*(D- 1.00)+0.817
        IF((D .GE. 2.03).AND.(D .LT. 2.48)) AV= 1.907*(D- 2.03)+1.219
        IF(D .GE. 2.48) AV= A0*(D- 2.48)+ 2.077
      RETURN
      ENDIF
      IF(CELL .EQ.  59) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.25)) AV= 0.197*D
        IF(D .GE. 1.25) AV= A0*(D- 1.25)+ 0.246
      RETURN
      ENDIF
      IF(CELL .EQ.  60) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.97)) AV= 0.340*D
        IF((D .GE. 0.97).AND.(D .LT. 1.18)) AV= 3.367*(D- 0.97)+0.330
        IF((D .GE. 1.18).AND.(D .LT. 2.94)) AV= 0.599*(D- 1.18)+1.037
        IF(D .GE. 2.94) AV= A0*(D- 2.94)+ 2.091
      RETURN
      ENDIF
      IF(CELL .EQ.  61) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.79)) AV= 0.720*D
        IF(D .GE. 0.79) AV= A0*(D- 0.79)+ 0.569
      RETURN
      ENDIF
      IF(CELL .EQ.  62) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.24)) AV= 0.000
        IF((D .GE. 0.24).AND.(D .LT. 0.59)) AV= 1.302*(D- 0.24)
        IF((D .GE. 0.59).AND.(D .LT. 2.11)) AV= 0.484*(D- 0.59)+0.456
        IF(D .GE. 2.11) AV= A0*(D- 2.11)+ 1.192
      RETURN
      ENDIF
      IF(CELL .EQ.  63) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 0.228*D
        IF((D .GE. 0.14).AND.(D .LT. 0.79)) AV= 0.363*(D- 0.14)+0.032
        IF((D .GE. 0.79).AND.(D .LT. 1.26)) AV= 3.987*(D- 0.79)+0.268
        IF(D .GE. 1.26) AV= A0*(D- 1.26)+ 2.142
      RETURN
      ENDIF
      IF(CELL .EQ.  64) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 3.768*D
        IF((D .GE. 0.05).AND.(D .LT. 0.75)) AV= 0.188
        IF((D .GE. 0.75).AND.(D .LT. 1.02)) AV= 2.793*(D- 0.75)+0.188
        IF((D .GE. 1.02).AND.(D .LT. 1.96)) AV= 0.942
        IF(D .GE. 1.96) AV= A0*(D- 1.96)+ 0.942
      RETURN
      ENDIF
      IF(CELL .EQ.  65) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.23)) AV= 1.084*D
        IF((D .GE. 0.23).AND.(D .LT. 0.63)) AV= 0.249
        IF((D .GE. 0.63).AND.(D .LT. 1.00)) AV= 2.935*(D- 0.63)+0.249
        IF((D .GE. 1.00).AND.(D .LT. 1.82)) AV= 0.540*(D- 1.00)+1.335
        IF(D .GE. 1.82) AV= A0*(D- 1.82)+ 1.778
      RETURN
      ENDIF
      IF(CELL .EQ.  66) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.92)) AV= 2.968*D
        IF((D .GE. 0.92).AND.(D .LT. 1.99)) AV= 2.731
        IF(D .GE. 1.99) AV= A0*(D- 1.99)+ 2.731
      RETURN
      ENDIF
      IF(CELL .EQ.  67) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.74)) AV= 0.341*D
        IF((D .GE. 0.74).AND.(D .LT. 3.50)) AV= 0.252
        IF(D .GE. 3.50) AV= A0*(D- 3.50)+ 0.252
      RETURN
      ENDIF
      IF(CELL .EQ.  68) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.24)) AV= 0.000
        IF((D .GE. 0.24).AND.(D .LT. 0.29)) AV= 4.092*(D- 0.24)
        IF((D .GE. 0.29).AND.(D .LT. 1.53)) AV= 0.205
        IF((D .GE. 1.53).AND.(D .LT. 1.74)) AV= 5.227*(D- 1.53)+0.205
        IF((D .GE. 1.74).AND.(D .LT. 3.55)) AV= 1.303
        IF((D .GE. 3.55).AND.(D .LT. 3.69)) AV= 7.056*(D- 3.55)+1.303
        IF(D .GE. 3.69) AV= A0*(D- 3.69)+ 2.291
      RETURN
      ENDIF
      IF(CELL .EQ.  69) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.58)) AV= 0.452*D
        IF((D .GE. 0.58).AND.(D .LT. 0.96)) AV= 1.341*(D- 0.58)+0.262
        IF((D .GE. 0.96).AND.(D .LT. 2.99)) AV= 0.772
        IF((D .GE. 2.99).AND.(D .LT. 3.22)) AV= 5.414*(D- 2.99)+0.772
        IF(D .GE. 3.22) AV= A0*(D- 3.22)+ 2.017
      RETURN
      ENDIF
      IF(CELL .EQ.  70) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.98)) AV= 0.610*(D- 0.08)
        IF((D .GE. 0.98).AND.(D .LT. 2.49)) AV= 0.549
        IF(D .GE. 2.49) AV= A0*(D- 2.49)+ 0.549
      RETURN
      ENDIF
      IF(CELL .EQ.  71) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.35)) AV= 0.757*D
        IF((D .GE. 0.35).AND.(D .LT. 0.91)) AV= 0.265
        IF((D .GE. 0.91).AND.(D .LT. 1.14)) AV= 2.481*(D- 0.91)+0.265
        IF((D .GE. 1.14).AND.(D .LT. 1.74)) AV= 1.080*(D- 1.14)+0.836
        IF((D .GE. 1.74).AND.(D .LT. 3.31)) AV= 1.484
        IF(D .GE. 3.31) AV= A0*(D- 3.31)+ 1.484
      RETURN
      ENDIF
      IF(CELL .EQ.  72) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.31)) AV= 1.166*D
        IF((D .GE. 0.31).AND.(D .LT. 1.19)) AV= 0.361
        IF((D .GE. 1.19).AND.(D .LT. 1.36)) AV= 3.841*(D- 1.19)+0.361
        IF(D .GE. 1.36) AV= A0*(D- 1.36)+ 1.014
      RETURN
      ENDIF
      IF(CELL .EQ.  73) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.19)) AV= 1.401*D
        IF((D .GE. 0.19).AND.(D .LT. 0.96)) AV= 0.266
        IF((D .GE. 0.96).AND.(D .LT. 1.16)) AV= 6.370*(D- 0.96)+0.266
        IF(D .GE. 1.16) AV= A0*(D- 1.16)+ 1.540
      RETURN
      ENDIF
      IF(CELL .EQ.  74) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.41)) AV= 0.638*D
        IF((D .GE. 0.41).AND.(D .LT. 0.85)) AV= 0.262
        IF((D .GE. 0.85).AND.(D .LT. 1.09)) AV= 1.529*(D- 0.85)+0.262
        IF((D .GE. 1.09).AND.(D .LT. 1.90)) AV= 0.418*(D- 1.09)+0.629
        IF(D .GE. 1.90) AV= A0*(D- 1.90)+ 0.968
      RETURN
      ENDIF
      IF(CELL .EQ.  75) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.25)) AV= 0.865*D
        IF((D .GE. 0.25).AND.(D .LT. 0.82)) AV= 0.216
        IF((D .GE. 0.82).AND.(D .LT. 0.94)) AV= 2.478*(D- 0.82)+0.216
        IF((D .GE. 0.94).AND.(D .LT. 1.31)) AV= 0.529*(D- 0.94)+0.513
        IF((D .GE. 1.31).AND.(D .LT. 2.61)) AV= 0.709
        IF(D .GE. 2.61) AV= A0*(D- 2.61)+ 0.709
      RETURN
      ENDIF
      IF(CELL .EQ.  76) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.39)) AV= 0.824*D
        IF((D .GE. 0.39).AND.(D .LT. 0.78)) AV= 0.321
        IF((D .GE. 0.78).AND.(D .LT. 0.87)) AV=13.240*(D- 0.78)+0.321
        IF((D .GE. 0.87).AND.(D .LT. 2.01)) AV= 1.513
        IF(D .GE. 2.01) AV= A0*(D- 2.01)+ 1.513
      RETURN
      ENDIF
      IF(CELL .EQ.  77) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 3.150*D
        IF((D .GE. 0.10).AND.(D .LT. 1.21)) AV= 0.315
        IF((D .GE. 1.21).AND.(D .LT. 1.79)) AV= 1.182*(D- 1.21)+0.315
        IF((D .GE. 1.79).AND.(D .LT. 5.76)) AV= 1.001
        IF(D .GE. 5.76) AV= A0*(D- 5.76)+ 1.001
      RETURN
      ENDIF
      IF(CELL .EQ.  78) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 3.619*D
        IF((D .GE. 0.14).AND.(D .LT. 1.42)) AV= 0.507
        IF(D .GE. 1.42) AV= A0*(D- 1.42)+ 0.507
      RETURN
      ENDIF
      IF(CELL .EQ.  79) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.22)) AV= 1.247*D
        IF((D .GE. 0.22).AND.(D .LT. 1.88)) AV= 0.274
        IF(D .GE. 1.88) AV= A0*(D- 1.88)+ 0.274
      RETURN
      ENDIF
      IF(CELL .EQ.  80) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.38)) AV= 0.000
        IF((D .GE. 0.38).AND.(D .LT. 0.70)) AV= 3.929*(D- 0.38)
        IF((D .GE. 0.70).AND.(D .LT. 4.85)) AV= 1.257
        IF(D .GE. 4.85) AV= A0*(D- 4.85)+ 1.257
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK3(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ.  81) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 4.525*D
        IF((D .GE. 0.12).AND.(D .LT. 1.01)) AV= 0.543
        IF((D .GE. 1.01).AND.(D .LT. 1.22)) AV= 4.566*(D- 1.01)+0.543
        IF((D .GE. 1.22).AND.(D .LT. 3.68)) AV= 1.502
        IF(D .GE. 3.68) AV= A0*(D- 3.68)+ 1.502
      RETURN
      ENDIF
      IF(CELL .EQ.  82) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 4.875*D
        IF((D .GE. 0.07).AND.(D .LT. 1.14)) AV= 0.341
        IF((D .GE. 1.14).AND.(D .LT. 2.60)) AV= 1.572*(D- 1.14)+0.341
        IF((D .GE. 2.60).AND.(D .LT. 5.45)) AV= 2.636
        IF(D .GE. 5.45) AV= A0*(D- 5.45)+ 2.636
      RETURN
      ENDIF
      IF(CELL .EQ.  83) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.00)) AV= 0.000
        IF((D .GE. 1.00).AND.(D .LT. 1.40)) AV= 2.708*(D- 1.00)
        IF((D .GE. 1.40).AND.(D .LT. 3.99)) AV= 1.083
        IF(D .GE. 3.99) AV= A0*(D- 3.99)+ 1.083
      RETURN
      ENDIF
      IF(CELL .EQ.  84) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.86)) AV= 0.000
        IF((D .GE. 0.86).AND.(D .LT. 1.18)) AV= 2.221*(D- 0.86)
        IF((D .GE. 1.18).AND.(D .LT. 3.43)) AV= 0.711
        IF(D .GE. 3.43) AV= A0*(D- 3.43)+ 0.711
      RETURN
      ENDIF
      IF(CELL .EQ.  85) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 1.287*D
        IF((D .GE. 0.10).AND.(D .LT. 1.28)) AV= 0.129
        IF((D .GE. 1.28).AND.(D .LT. 1.52)) AV= 0.754*(D- 1.28)+0.129
        IF((D .GE. 1.52).AND.(D .LT. 5.24)) AV= 0.310
        IF(D .GE. 5.24) AV= A0*(D- 5.24)+ 0.310
      RETURN
      ENDIF
      IF(CELL .EQ.  86) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 2.750*D
        IF((D .GE. 0.08).AND.(D .LT. 1.01)) AV= 0.220
        IF((D .GE. 1.01).AND.(D .LT. 1.81)) AV= 1.468*(D- 1.01)+0.220
        IF((D .GE. 1.81).AND.(D .LT. 4.90)) AV= 1.394
        IF(D .GE. 4.90) AV= A0*(D- 4.90)+ 1.394
      RETURN
      ENDIF
      IF(CELL .EQ.  87) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.38)) AV= 0.000
        IF((D .GE. 0.38).AND.(D .LT. 1.02)) AV= 1.542*(D- 0.38)
        IF((D .GE. 1.02).AND.(D .LT. 4.91)) AV= 0.987
        IF(D .GE. 4.91) AV= A0*(D- 4.91)+ 0.987
      RETURN
      ENDIF
      IF(CELL .EQ.  88) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.75)) AV= 0.372*D
        IF((D .GE. 0.75).AND.(D .LT. 3.50)) AV= 0.174*(D- 0.75)+0.279
        IF((D .GE. 3.50).AND.(D .LT. 3.78)) AV= 6.305*(D- 3.50)+0.758
        IF(D .GE. 3.78) AV= A0*(D- 3.78)+ 2.523
      RETURN
      ENDIF
      IF(CELL .EQ.  89) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.26)) AV= 1.884*(D- 0.08)
        IF((D .GE. 0.26).AND.(D .LT. 1.44)) AV= 0.339
        IF((D .GE. 1.44).AND.(D .LT. 1.96)) AV= 1.153*(D- 1.44)+0.339
        IF((D .GE. 1.96).AND.(D .LT. 4.99)) AV= 0.939
        IF((D .GE. 4.99).AND.(D .LT. 5.49)) AV= 2.266*(D- 4.99)+0.939
        IF(D .GE. 5.49) AV= A0*(D- 5.49)+ 2.072
      RETURN
      ENDIF
      IF(CELL .EQ.  90) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.25)) AV= 1.000*D
        IF((D .GE. 0.25).AND.(D .LT. 1.58)) AV= 0.250
        IF((D .GE. 1.58).AND.(D .LT. 2.08)) AV= 2.117*(D- 1.58)+0.250
        IF((D .GE. 2.08).AND.(D .LT. 5.52)) AV= 1.309
        IF(D .GE. 5.52) AV= A0*(D- 5.52)+ 1.309
      RETURN
      ENDIF
      IF(CELL .EQ.  91) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 1.004*D
        IF((D .GE. 0.15).AND.(D .LT. 1.78)) AV= 0.151
        IF((D .GE. 1.78).AND.(D .LT. 2.00)) AV= 6.681*(D- 1.78)+0.151
        IF((D .GE. 2.00).AND.(D .LT. 5.59)) AV= 1.621
        IF(D .GE. 5.59) AV= A0*(D- 5.59)+ 1.621
      RETURN
      ENDIF
      IF(CELL .EQ.  92) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.34)) AV= 0.505*D
        IF((D .GE. 0.34).AND.(D .LT. 1.74)) AV= 0.172
        IF((D .GE. 1.74).AND.(D .LT. 2.03)) AV= 5.743*(D- 1.74)+0.172
        IF((D .GE. 2.03).AND.(D .LT. 5.85)) AV= 1.837
        IF(D .GE. 5.85) AV= A0*(D- 5.85)+ 1.837
      RETURN
      ENDIF
      IF(CELL .EQ.  93) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 1.004*D
        IF((D .GE. 0.15).AND.(D .LT. 1.74)) AV= 0.151
        IF((D .GE. 1.74).AND.(D .LT. 2.03)) AV= 2.750*(D- 1.74)+0.151
        IF((D .GE. 2.03).AND.(D .LT. 5.71)) AV= 0.949
        IF(D .GE. 5.71) AV= A0*(D- 5.71)+ 0.949
      RETURN
      ENDIF
      IF(CELL .EQ.  94) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 0.449*D
        IF((D .GE. 0.10).AND.(D .LT. 0.82)) AV= 0.045
        IF((D .GE. 0.82).AND.(D .LT. 0.93)) AV= 1.956*(D- 0.82)+0.045
        IF((D .GE. 0.93).AND.(D .LT. 1.05)) AV= 0.260
        IF((D .GE. 1.05).AND.(D .LT. 1.22)) AV=11.801*(D- 1.05)+0.260
        IF(D .GE. 1.22) AV= A0*(D- 1.22)+ 2.266
      RETURN
      ENDIF
      IF(CELL .EQ.  95) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.09)) AV= 0.953*D
        IF((D .GE. 1.09).AND.(D .LT. 5.03)) AV= 1.039
        IF(D .GE. 5.03) AV= A0*(D- 5.03)+ 1.039
      RETURN
      ENDIF
      IF(CELL .EQ.  96) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.24)) AV= 1.159*D
        IF((D .GE. 0.24).AND.(D .LT. 1.85)) AV= 0.278
        IF((D .GE. 1.85).AND.(D .LT. 3.53)) AV= 0.245*(D- 1.85)+0.278
        IF((D .GE. 3.53).AND.(D .LT. 5.87)) AV= 0.690
        IF(D .GE. 5.87) AV= A0*(D- 5.87)+ 0.690
      RETURN
      ENDIF
      IF(CELL .EQ.  97) THEN
        IF((D .GE. 0.00).AND.(D .LT. 2.30)) AV= 0.573*D
        IF(D .GE. 2.30) AV= A0*(D- 2.30)+ 1.318
                RETURN
      ENDIF
      IF(CELL .EQ.  98) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 2.943*D
        IF((D .GE. 0.07).AND.(D .LT. 0.49)) AV= 0.206
        IF((D .GE. 0.49).AND.(D .LT. 1.04)) AV= 0.792*(D- 0.49)+0.206
        IF((D .GE. 1.04).AND.(D .LT. 1.98)) AV= 0.642
        IF((D .GE. 1.98).AND.(D .LT. 2.25)) AV= 5.016*(D- 1.98)+0.642
        IF((D .GE. 2.25).AND.(D .LT. 5.44)) AV= 1.996
        IF(D .GE. 5.44) AV= A0*(D- 5.44)+ 1.996
      RETURN
      ENDIF
      IF(CELL .EQ.  99) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.17)) AV= 1.233*D
        IF((D .GE. 0.17).AND.(D .LT. 1.19)) AV= 0.210
        IF((D .GE. 1.19).AND.(D .LT. 1.43)) AV= 7.109*(D- 1.19)+0.210
        IF((D .GE. 1.43).AND.(D .LT. 3.71)) AV= 1.916
        IF(D .GE. 3.71) AV= A0*(D- 3.71)+ 1.916
      RETURN
      ENDIF
      IF(CELL .EQ. 100) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.57)) AV= 0.475*D
        IF((D .GE. 0.57).AND.(D .LT. 0.99)) AV= 3.645*(D- 0.57)+0.271
        IF((D .GE. 0.99).AND.(D .LT. 2.46)) AV= 1.802
        IF(D .GE. 2.46) AV= A0*(D- 2.46)+ 1.802
      RETURN
      ENDIF
      IF(CELL .EQ. 101) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.46)) AV= 0.316*D
        IF((D .GE. 0.46).AND.(D .LT. 1.03)) AV= 1.282*(D- 0.46)+0.145
        IF((D .GE. 1.03).AND.(D .LT. 3.74)) AV= 0.876
        IF(D .GE. 3.74) AV= A0*(D- 3.74)+ 0.876
      RETURN
      ENDIF
      IF(CELL .EQ. 102) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 1.793*D
        IF((D .GE. 0.15).AND.(D .LT. 0.77)) AV= 0.269
        IF((D .GE. 0.77).AND.(D .LT. 1.17)) AV= 1.717*(D- 0.77)+0.269
        IF((D .GE. 1.17).AND.(D .LT. 4.74)) AV= 0.956
        IF(D .GE. 4.74) AV= A0*(D- 4.74)+ 0.956
      RETURN
      ENDIF
      IF(CELL .EQ. 103) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 0.000
        IF((D .GE. 0.07).AND.(D .LT. 0.33)) AV= 1.112*(D- 0.07)
        IF((D .GE. 0.33).AND.(D .LT. 0.77)) AV= 0.289
        IF((D .GE. 0.77).AND.(D .LT. 1.00)) AV= 8.498*(D- 0.77)+0.289
        IF(D .GE. 1.00) AV= A0*(D- 1.00)+ 2.244
      RETURN
      ENDIF
      IF(CELL .EQ. 104) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.04)) AV= 3.876*D
        IF((D .GE. 0.04).AND.(D .LT. 0.75)) AV= 0.155
        IF((D .GE. 0.75).AND.(D .LT. 1.02)) AV= 3.239*(D- 0.75)+0.155
        IF((D .GE. 1.02).AND.(D .LT. 2.11)) AV= 1.030
        IF(D .GE. 2.11) AV= A0*(D- 2.11)+ 1.030
      RETURN
      ENDIF
      IF(CELL .EQ. 105) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.000
        IF((D .GE. 0.11).AND.(D .LT. 1.01)) AV= 0.434*(D- 0.11)
        IF((D .GE. 1.01).AND.(D .LT. 1.26)) AV= 4.448*(D- 1.01)+0.391
        IF((D .GE. 1.26).AND.(D .LT. 4.85)) AV= 1.503
        IF(D .GE. 4.85) AV= A0*(D- 4.85)+ 1.503
      RETURN
      ENDIF
      IF(CELL .EQ. 106) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.98)) AV= 0.000
        IF((D .GE. 0.98).AND.(D .LT. 1.25)) AV= 5.805*(D- 0.98)
        IF((D .GE. 1.25).AND.(D .LT. 3.54)) AV= 1.567
        IF(D .GE. 3.54) AV= A0*(D- 3.54)+ 1.567
      RETURN
      ENDIF
      IF(CELL .EQ. 107) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.06)) AV= 0.000
        IF((D .GE. 0.06).AND.(D .LT. 0.42)) AV= 0.595*(D- 0.06)
        IF((D .GE. 0.42).AND.(D .LT. 0.77)) AV= 0.214
        IF((D .GE. 0.77).AND.(D .LT. 1.91)) AV= 2.627*(D- 0.77)+0.214
        IF(D .GE. 1.91) AV= A0*(D- 1.91)+ 3.209
      RETURN
      ENDIF
      IF(CELL .EQ. 108) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.04)) AV= 3.513*D
        IF((D .GE. 0.04).AND.(D .LT. 0.79)) AV= 0.141
        IF((D .GE. 0.79).AND.(D .LT. 1.28)) AV= 4.199*(D- 0.79)+0.141
        IF(D .GE. 1.28) AV= A0*(D- 1.28)+ 2.199
      RETURN
      ENDIF
      IF(CELL .EQ. 109) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.13)) AV= 0.962*D
        IF((D .GE. 0.13).AND.(D .LT. 0.83)) AV= 0.125
        IF((D .GE. 0.83).AND.(D .LT. 1.31)) AV= 1.777*(D- 0.83)+0.125
        IF((D .GE. 1.31).AND.(D .LT. 3.18)) AV= 0.978
        IF(D .GE. 3.18) AV= A0*(D- 3.18)+ 0.978
      RETURN
      ENDIF
      IF(CELL .EQ. 110) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 0.488*D
        IF((D .GE. 0.11).AND.(D .LT. 0.85)) AV= 0.124*(D- 0.11)+0.054
        IF((D .GE. 0.85).AND.(D .LT. 1.27)) AV= 5.639*(D- 0.85)+0.146
        IF((D .GE. 1.27).AND.(D .LT. 2.28)) AV= 2.514
        IF(D .GE. 2.28) AV= A0*(D- 2.28)+ 2.514
      RETURN
      ENDIF
      IF(CELL .EQ. 111) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 0.000
        IF((D .GE. 0.09).AND.(D .LT. 0.27)) AV= 1.019*(D- 0.09)
        IF((D .GE. 0.27).AND.(D .LT. 0.96)) AV= 0.183
        IF((D .GE. 0.96).AND.(D .LT. 1.29)) AV=10.321*(D- 0.96)+0.183
        IF((D .GE. 1.29).AND.(D .LT. 2.44)) AV= 3.589
        IF(D .GE. 2.44) AV= A0*(D- 2.44)+ 3.589
      RETURN
      ENDIF
      IF(CELL .EQ. 112) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.00)) AV= 0.362*D
        IF((D .GE. 1.00).AND.(D .LT. 2.23)) AV= 0.928*(D- 1.00)+0.362
        IF(D .GE. 2.23) AV= A0*(D- 2.23)+ 1.503
      RETURN
      ENDIF
      IF(CELL .EQ. 113) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 1.392*D
        IF((D .GE. 0.16).AND.(D .LT. 0.83)) AV= 0.223
        IF((D .GE. 0.83).AND.(D .LT. 1.32)) AV= 4.222*(D- 0.83)+0.223
        IF(D .GE. 1.32) AV= A0*(D- 1.32)+ 2.292
      RETURN
      ENDIF
      IF(CELL .EQ. 114) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.74)) AV= 0.324*D
        IF((D .GE. 0.74).AND.(D .LT. 2.45)) AV= 0.422*(D- 0.74)+0.240
        IF((D .GE. 2.45).AND.(D .LT. 3.22)) AV= 0.962
        IF(D .GE. 3.22) AV= A0*(D- 3.22)+ 0.962
      RETURN
      ENDIF
      IF(CELL .EQ. 115) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.01)) AV=10.769*D
        IF((D .GE. 0.01).AND.(D .LT. 0.90)) AV= 0.108
        IF((D .GE. 0.90).AND.(D .LT. 1.26)) AV= 5.363*(D- 0.90)+0.108
        IF((D .GE. 1.26).AND.(D .LT. 3.94)) AV= 2.039
        IF(D .GE. 3.94) AV= A0*(D- 3.94)+ 2.039
      RETURN
      ENDIF
      IF(CELL .EQ. 116) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.76)) AV= 0.865*D
        IF((D .GE. 0.76).AND.(D .LT. 7.76)) AV= 0.657
        IF(D .GE. 7.76) AV= A0*(D- 7.76)+ 0.657
      RETURN
      ENDIF
      IF(CELL .EQ. 117) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.32)) AV= 0.000
        IF((D .GE. 0.32).AND.(D .LT. 0.68)) AV= 1.543*(D- 0.32)
        IF((D .GE. 0.68).AND.(D .LT. 1.59)) AV= 0.555
        IF(D .GE. 1.59) AV= A0*(D- 1.59)+ 0.555
      RETURN
      ENDIF
      IF(CELL .EQ. 118) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.31)) AV= 0.532*D
        IF((D .GE. 0.31).AND.(D .LT. 0.73)) AV= 0.165
        IF((D .GE. 0.73).AND.(D .LT. 1.19)) AV= 2.133*(D- 0.73)+0.165
        IF((D .GE. 1.19).AND.(D .LT. 1.52)) AV= 1.434*(D- 1.19)+1.146
        IF((D .GE. 1.52).AND.(D .LT. 4.55)) AV= 1.619
        IF(D .GE. 4.55) AV= A0*(D- 4.55)+ 1.619
      RETURN
      ENDIF
      IF(CELL .EQ. 119) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.13)) AV= 2.710*D
        IF((D .GE. 0.13).AND.(D .LT. 0.72)) AV= 0.352
        IF((D .GE. 0.72).AND.(D .LT. 1.02)) AV= 1.410*(D- 0.72)+0.352
        IF(D .GE. 1.02) AV= A0*(D- 1.02)+ 0.775
      RETURN
      ENDIF
      IF(CELL .EQ. 120) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.13)) AV= 1.078*D
        IF((D .GE. 0.13).AND.(D .LT. 0.74)) AV= 0.140
        IF((D .GE. 0.74).AND.(D .LT. 1.26)) AV= 1.881*(D- 0.74)+0.140
        IF((D .GE. 1.26).AND.(D .LT. 3.50)) AV= 1.118
        IF(D .GE. 3.50) AV= A0*(D- 3.50)+ 1.118
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK4(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ. 121) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 6.086*D
        IF((D .GE. 0.07).AND.(D .LT. 0.82)) AV= 0.426
        IF((D .GE. 0.82).AND.(D .LT. 1.29)) AV= 1.423*(D- 0.82)+0.426
        IF((D .GE. 1.29).AND.(D .LT. 4.58)) AV= 1.095
        IF(D .GE. 4.58) AV= A0*(D- 4.58)+ 1.095
      RETURN
      ENDIF
      IF(CELL .EQ. 122) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 0.202*D
        IF((D .GE. 0.10).AND.(D .LT. 0.76)) AV= 0.790*(D- 0.10)+0.020
        IF((D .GE. 0.76).AND.(D .LT. 3.16)) AV= 0.541
        IF(D .GE. 3.16) AV= A0*(D- 3.16)+ 0.541
      RETURN
      ENDIF
      IF(CELL .EQ. 123) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 5.247*D
        IF((D .GE. 0.05).AND.(D .LT. 1.41)) AV= 0.262
        IF((D .GE. 1.41).AND.(D .LT. 2.09)) AV= 2.610*(D- 1.41)+0.262
        IF((D .GE. 2.09).AND.(D .LT. 5.61)) AV= 2.037
        IF(D .GE. 5.61) AV= A0*(D- 5.61)+ 2.037
      RETURN
      ENDIF
      IF(CELL .EQ. 124) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 1.782*D
        IF((D .GE. 0.16).AND.(D .LT. 1.54)) AV= 0.285
        IF((D .GE. 1.54).AND.(D .LT. 2.01)) AV= 1.904*(D- 1.54)+0.285
        IF((D .GE. 2.01).AND.(D .LT. 5.49)) AV= 1.180
        IF(D .GE. 5.49) AV= A0*(D- 5.49)+ 1.180
      RETURN
      ENDIF
      IF(CELL .EQ. 125) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.97)) AV= 0.000
        IF((D .GE. 0.97).AND.(D .LT. 3.17)) AV= 0.920*(D- 0.97)
        IF((D .GE. 3.17).AND.(D .LT. 5.69)) AV= 2.024
        IF(D .GE. 5.69) AV= A0*(D- 5.69)+ 2.024
      RETURN
      ENDIF
      IF(CELL .EQ. 126) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 2.492*D
        IF((D .GE. 0.10).AND.(D .LT. 1.22)) AV= 0.249
        IF((D .GE. 1.22).AND.(D .LT. 1.77)) AV= 2.289*(D- 1.22)+0.249
        IF((D .GE. 1.77).AND.(D .LT. 5.61)) AV= 1.508
        IF(D .GE. 5.61) AV= A0*(D- 5.61)+ 1.508
      RETURN
      ENDIF
      IF(CELL .EQ. 127) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 2.896*D
        IF((D .GE. 0.11).AND.(D .LT. 1.50)) AV= 0.319
        IF((D .GE. 1.50).AND.(D .LT. 2.03)) AV= 1.616*(D- 1.50)+0.319
        IF((D .GE. 2.03).AND.(D .LT. 5.78)) AV= 1.175
        IF(D .GE. 5.78) AV= A0*(D- 5.78)+ 1.175
      RETURN
      ENDIF
      IF(CELL .EQ. 128) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 3.719*D
        IF((D .GE. 0.08).AND.(D .LT. 1.50)) AV= 0.298
        IF((D .GE. 1.50).AND.(D .LT. 2.01)) AV= 0.999*(D- 1.50)+0.298
        IF((D .GE. 2.01).AND.(D .LT. 5.61)) AV= 0.807
        IF(D .GE. 5.61) AV= A0*(D- 5.61)+ 0.807
      RETURN
      ENDIF
      IF(CELL .EQ. 129) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.85)) AV= 1.608*D
        IF((D .GE. 0.85).AND.(D .LT. 5.65)) AV= 1.367
        IF(D .GE. 5.65) AV= A0*(D- 5.65)+ 1.367
      RETURN
      ENDIF
      IF(CELL .EQ. 130) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.96)) AV= 0.439*D
        IF((D .GE. 0.96).AND.(D .LT. 2.81)) AV= 0.421
        IF((D .GE. 2.81).AND.(D .LT. 3.28)) AV= 0.684*(D- 2.81)+0.421
        IF((D .GE. 3.28).AND.(D .LT. 5.05)) AV= 0.742
        IF(D .GE. 5.05) AV= A0*(D- 5.05)+ 0.742
      RETURN
      ENDIF
      IF(CELL .EQ. 131) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.50)) AV= 0.678*D
        IF((D .GE. 1.50).AND.(D .LT. 5.84)) AV= 1.017
        IF(D .GE. 5.84) AV= A0*(D- 5.84)+ 1.017
      RETURN
      ENDIF
      IF(CELL .EQ. 132) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.23)) AV= 1.409*D
        IF((D .GE. 1.23).AND.(D .LT. 3.22)) AV= 1.733
        IF(D .GE. 3.22) AV= A0*(D- 3.22)+ 1.733
      RETURN
      ENDIF
      IF(CELL .EQ. 133) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 1.998*D
        IF((D .GE. 0.11).AND.(D .LT. 0.74)) AV= 0.220
        IF((D .GE. 0.74).AND.(D .LT. 1.74)) AV= 1.570*(D- 0.74)+0.220
        IF((D .GE. 1.74).AND.(D .LT. 2.39)) AV= 1.790
        IF((D .GE. 2.39).AND.(D .LT. 3.02)) AV= 0.811*(D- 2.39)+1.790
        IF((D .GE. 3.02).AND.(D .LT. 4.70)) AV= 2.301
        IF(D .GE. 4.70) AV= A0*(D- 4.70)+ 2.301
      RETURN
      ENDIF
      IF(CELL .EQ. 134) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 1.107*D
        IF((D .GE. 0.14).AND.(D .LT. 1.18)) AV= 0.155
        IF((D .GE. 1.18).AND.(D .LT. 1.57)) AV= 2.033*(D- 1.18)+0.155
        IF((D .GE. 1.57).AND.(D .LT. 4.18)) AV= 0.948
        IF(D .GE. 4.18) AV= A0*(D- 4.18)+ 0.948
      RETURN
      ENDIF
      IF(CELL .EQ. 135) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.04)) AV= 1.567*D
        IF((D .GE. 0.04).AND.(D .LT. 1.23)) AV= 0.063
        IF((D .GE. 1.23).AND.(D .LT. 1.60)) AV= 4.639*(D- 1.23)+0.063
        IF((D .GE. 1.60).AND.(D .LT. 5.45)) AV= 1.779
        IF(D .GE. 5.45) AV= A0*(D- 5.45)+ 1.779
      RETURN
      ENDIF
      IF(CELL .EQ. 136) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.13)) AV= 2.054*D
        IF((D .GE. 0.13).AND.(D .LT. 0.77)) AV= 0.267
        IF((D .GE. 0.77).AND.(D .LT. 1.04)) AV= 3.505*(D- 0.77)+0.267
        IF((D .GE. 1.04).AND.(D .LT. 5.76)) AV= 1.213
        IF(D .GE. 5.76) AV= A0*(D- 5.76)+ 1.213
      RETURN
      ENDIF
      IF(CELL .EQ. 137) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 0.000
        IF((D .GE. 0.12).AND.(D .LT. 0.71)) AV= 1.310*(D- 0.12)
        IF((D .GE. 0.71).AND.(D .LT. 4.49)) AV= 0.773
        IF((D .GE. 4.49).AND.(D .LT. 4.66)) AV= 3.144*(D- 4.49)+0.773
        IF(D .GE. 4.66) AV= A0*(D- 4.66)+ 1.307
      RETURN
      ENDIF
      IF(CELL .EQ. 138) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.20)) AV= 1.644*D
        IF((D .GE. 0.20).AND.(D .LT. 0.66)) AV= 0.563*(D- 0.20)+0.329
        IF((D .GE. 0.66).AND.(D .LT. 0.90)) AV= 1.455*(D- 0.66)+0.588
        IF((D .GE. 0.90).AND.(D .LT. 5.32)) AV= 0.937
        IF(D .GE. 5.32) AV= A0*(D- 5.32)+ 0.937
      RETURN
      ENDIF
      IF(CELL .EQ. 139) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.21)) AV= 0.000
        IF((D .GE. 0.21).AND.(D .LT. 0.34)) AV= 6.427*(D- 0.21)
        IF((D .GE. 0.34).AND.(D .LT. 0.77)) AV= 0.836
        IF((D .GE. 0.77).AND.(D .LT. 1.05)) AV= 2.354*(D- 0.77)+0.836
        IF((D .GE. 1.05).AND.(D .LT. 3.73)) AV= 1.495
        IF(D .GE. 3.73) AV= A0*(D- 3.73)+ 1.495
      RETURN
      ENDIF
      IF(CELL .EQ. 140) THEN
        IF((D .GE. 0.00).AND.(D .LT. 2.66)) AV= 0.609*D
        IF((D .GE. 2.66).AND.(D .LT. 5.53)) AV= 1.620
        IF(D .GE. 5.53) AV= A0*(D- 5.53)+ 1.620
      RETURN
      ENDIF
      IF(CELL .EQ. 141) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 1.900*D
        IF((D .GE. 0.16).AND.(D .LT. 0.75)) AV= 0.304
        IF((D .GE. 0.75).AND.(D .LT. 1.03)) AV= 6.555*(D- 0.75)+0.304
        IF((D .GE. 1.03).AND.(D .LT. 2.56)) AV= 2.139
        IF(D .GE. 2.56) AV= A0*(D- 2.56)+ 2.139
      RETURN
      ENDIF
      IF(CELL .EQ. 142) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.70)) AV= 0.000
        IF((D .GE. 0.70).AND.(D .LT. 1.01)) AV= 2.951*(D- 0.70)
        IF((D .GE. 1.01).AND.(D .LT. 3.95)) AV= 0.915
        IF(D .GE. 3.95) AV= A0*(D- 3.95)+ 0.915
      RETURN
      ENDIF
      IF(CELL .EQ. 143) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.27)) AV= 4.867*D
        IF((D .GE. 0.27).AND.(D .LT. 3.86)) AV= 1.314
        IF((D .GE. 3.86).AND.(D .LT. 4.05)) AV= 4.304*(D- 3.86)+1.314
        IF(D .GE. 4.05) AV= A0*(D- 4.05)+ 2.137
      RETURN
      ENDIF
      IF(CELL .EQ. 144) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.75)) AV= 0.000
        IF((D .GE. 0.75).AND.(D .LT. 0.81)) AV= 2.033*(D- 0.75)
        IF((D .GE. 0.81).AND.(D .LT. 1.02)) AV= 2.755*(D- 0.81)+0.122
        IF((D .GE. 1.02).AND.(D .LT. 1.92)) AV= 0.701
        IF((D .GE. 1.92).AND.(D .LT. 2.12)) AV= 2.953*(D- 1.92)+0.701
        IF((D .GE. 2.12).AND.(D .LT. 4.73)) AV= 1.292
        IF(D .GE. 4.73) AV= A0*(D- 4.73)+ 1.292
      RETURN
      ENDIF
      IF(CELL .EQ. 145) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.29)) AV= 1.173*D
        IF((D .GE. 1.29).AND.(D .LT. 3.21)) AV= 1.513
        IF(D .GE. 3.21) AV= A0*(D- 3.21)+ 1.513
      RETURN
      ENDIF
      IF(CELL .EQ. 146) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.40)) AV= 2.081*D
        IF((D .GE. 0.40).AND.(D .LT. 1.89)) AV= 0.832
        IF((D .GE. 1.89).AND.(D .LT. 2.04)) AV= 1.433*(D- 1.89)+0.832
        IF((D .GE. 2.04).AND.(D .LT. 5.10)) AV= 1.047
        IF(D .GE. 5.10) AV= A0*(D- 5.10)+ 1.047
      RETURN
      ENDIF
      IF(CELL .EQ. 147) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.22)) AV= 0.000
        IF((D .GE. 0.22).AND.(D .LT. 0.39)) AV= 7.221*(D- 0.22)
        IF((D .GE. 0.39).AND.(D .LT. 1.76)) AV= 1.228
        IF((D .GE. 1.76).AND.(D .LT. 2.04)) AV= 0.570*(D- 1.76)+1.228
        IF((D .GE. 2.04).AND.(D .LT. 5.14)) AV= 1.388
        IF(D .GE. 5.14) AV= A0*(D- 5.14)+ 1.388
      RETURN
      ENDIF
      IF(CELL .EQ. 148) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.31)) AV= 2.466*D
        IF((D .GE. 0.31).AND.(D .LT. 5.62)) AV= 0.764
        IF(D .GE. 5.62) AV= A0*(D- 5.62)+ 0.764
      RETURN
      ENDIF
      IF(CELL .EQ. 149) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.39)) AV= 0.000
        IF((D .GE. 0.39).AND.(D .LT. 0.71)) AV= 3.095*(D- 0.39)
        IF((D .GE. 0.71).AND.(D .LT. 5.93)) AV= 0.990
        IF(D .GE. 5.93) AV= A0*(D- 5.93)+ 0.990
      RETURN
      ENDIF
      IF(CELL .EQ. 150) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 0.000
        IF((D .GE. 0.09).AND.(D .LT. 0.37)) AV= 4.176*(D- 0.09)
        IF((D .GE. 0.37).AND.(D .LT. 5.28)) AV= 1.169
        IF(D .GE. 5.28) AV= A0*(D- 5.28)+ 1.169
      RETURN
      ENDIF
      IF(CELL .EQ. 151) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.21)) AV= 0.435*D
        IF((D .GE. 0.21).AND.(D .LT. 0.52)) AV= 1.329*(D- 0.21)+0.091
        IF((D .GE. 0.52).AND.(D .LT. 0.63)) AV= 0.503
        IF((D .GE. 0.63).AND.(D .LT. 0.75)) AV= 2.701*(D- 0.63)+0.503
        IF((D .GE. 0.75).AND.(D .LT. 4.29)) AV= 0.827
        IF(D .GE. 4.29) AV= A0*(D- 4.29)+ 0.827
      RETURN
      ENDIF
      IF(CELL .EQ. 152) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 0.000
        IF((D .GE. 0.16).AND.(D .LT. 0.56)) AV= 1.301*(D- 0.16)
        IF((D .GE. 0.56).AND.(D .LT. 0.89)) AV= 1.937*(D- 0.56)+0.520
        IF((D .GE. 0.89).AND.(D .LT. 4.27)) AV= 1.159
        IF((D .GE. 4.27).AND.(D .LT. 4.91)) AV= 2.107*(D- 4.27)+1.159
        IF(D .GE. 4.91) AV= A0*(D- 4.91)+ 2.507
      RETURN
      ENDIF
      IF(CELL .EQ. 153) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.06)) AV= 0.000
        IF((D .GE. 0.06).AND.(D .LT. 0.73)) AV= 0.963*(D- 0.06)
        IF((D .GE. 0.73).AND.(D .LT. 4.32)) AV= 0.645
        IF((D .GE. 4.32).AND.(D .LT. 4.96)) AV= 2.127*(D- 4.32)+0.645
        IF(D .GE. 4.96) AV= A0*(D- 4.96)+ 2.006
      RETURN
      ENDIF
      IF(CELL .EQ. 154) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.88)) AV= 1.133*D
        IF((D .GE. 1.88).AND.(D .LT. 4.39)) AV= 2.130
        IF(D .GE. 4.39) AV= A0*(D- 4.39)+ 2.130
      RETURN
      ENDIF
      IF(CELL .EQ. 155) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 6.312*D
        IF((D .GE. 0.09).AND.(D .LT. 1.24)) AV= 0.568
        IF((D .GE. 1.24).AND.(D .LT. 1.55)) AV= 0.910*(D- 1.24)+0.568
        IF((D .GE. 1.55).AND.(D .LT. 4.01)) AV= 0.850
        IF((D .GE. 4.01).AND.(D .LT. 4.07)) AV=15.928*(D- 4.01)+0.850
        IF(D .GE. 4.07) AV= A0*(D- 4.07)+ 1.806
      RETURN
      ENDIF
      IF(CELL .EQ. 156) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.95)) AV= 1.241*D
        IF((D .GE. 0.95).AND.(D .LT. 3.47)) AV= 1.179
        IF(D .GE. 3.47) AV= A0*(D- 3.47)+ 1.179
      RETURN
      ENDIF
      IF(CELL .EQ. 157) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 0.000
        IF((D .GE. 0.12).AND.(D .LT. 0.29)) AV= 3.030*(D- 0.12)
        IF((D .GE. 0.29).AND.(D .LT. 1.23)) AV= 0.515
        IF((D .GE. 1.23).AND.(D .LT. 1.52)) AV= 2.083*(D- 1.23)+0.515
        IF((D .GE. 1.52).AND.(D .LT. 2.48)) AV= 1.119
        IF((D .GE. 2.48).AND.(D .LT. 2.79)) AV= 1.967*(D- 2.48)+1.119
        IF((D .GE. 2.79).AND.(D .LT. 6.90)) AV= 1.729
        IF(D .GE. 6.90) AV= A0*(D- 6.90)+ 1.729
      RETURN
      ENDIF
      IF(CELL .EQ. 158) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.05)) AV= 1.515*D
        IF((D .GE. 1.05).AND.(D .LT. 3.65)) AV= 1.591
        IF((D .GE. 3.65).AND.(D .LT. 3.75)) AV=10.402*(D- 3.65)+1.591
        IF((D .GE. 3.75).AND.(D .LT. 5.91)) AV= 2.631
        IF(D .GE. 5.91) AV= A0*(D- 5.91)+ 2.631
      RETURN
      ENDIF
      IF(CELL .EQ. 159) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.02)) AV= 4.343*D
        IF((D .GE. 0.02).AND.(D .LT. 0.24)) AV= 3.123*(D- 0.02)+0.087
        IF((D .GE. 0.24).AND.(D .LT. 0.75)) AV= 0.774
        IF((D .GE. 0.75).AND.(D .LT. 0.98)) AV= 1.682*(D- 0.75)+0.774
        IF((D .GE. 0.98).AND.(D .LT. 3.54)) AV= 1.161
        IF(D .GE. 3.54) AV= A0*(D- 3.54)+ 1.161
      RETURN
      ENDIF
      IF(CELL .EQ. 160) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.72)) AV= 0.000
        IF((D .GE. 0.72).AND.(D .LT. 1.03)) AV= 9.608*(D- 0.72)
        IF(D .GE. 1.03) AV= A0*(D- 1.03)+ 2.978
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK5(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ. 161) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.56)) AV= 1.621*D
        IF((D .GE. 1.56).AND.(D .LT. 4.91)) AV= 2.529
        IF(D .GE. 4.91) AV= A0*(D- 4.91)+ 2.529
      RETURN
      ENDIF
      IF(CELL .EQ. 162) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.53)) AV= 1.132*D
        IF((D .GE. 1.53).AND.(D .LT. 3.43)) AV= 1.732
        IF(D .GE. 3.43) AV= A0*(D- 3.43)+ 1.732
      RETURN
      ENDIF
      IF(CELL .EQ. 163) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.78)) AV= 0.814*D
        IF((D .GE. 0.78).AND.(D .LT. 2.50)) AV= 0.635
        IF((D .GE. 2.50).AND.(D .LT. 2.84)) AV= 0.635*(D- 2.50)+0.635
        IF((D .GE. 2.84).AND.(D .LT. 3.78)) AV= 0.851
        IF(D .GE. 3.78) AV= A0*(D- 3.78)+ 0.851
      RETURN
      ENDIF
      IF(CELL .EQ. 164) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 4.665*D
        IF((D .GE. 0.09).AND.(D .LT. 1.17)) AV= 0.420
        IF((D .GE. 1.17).AND.(D .LT. 2.04)) AV= 0.402*(D- 1.17)+0.420
        IF((D .GE. 2.04).AND.(D .LT. 5.91)) AV= 0.770
        IF(D .GE. 5.91) AV= A0*(D- 5.91)+ 0.770
      RETURN
      ENDIF
      IF(CELL .EQ. 165) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.48)) AV= 2.697*D
        IF((D .GE. 0.48).AND.(D .LT. 2.37)) AV= 1.295
        IF(D .GE. 2.37) AV= A0*(D- 2.37)+ 1.295
      RETURN
      ENDIF
      IF(CELL .EQ. 166) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.55)) AV= 2.888*D
        IF((D .GE. 0.55).AND.(D .LT. 2.72)) AV= 1.588
        IF((D .GE. 2.72).AND.(D .LT. 3.38)) AV= 2.981*(D- 2.72)+1.588
        IF(D .GE. 3.38) AV= A0*(D- 3.38)+ 3.555
      RETURN
      ENDIF
      IF(CELL .EQ. 167) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 0.160*D
        IF((D .GE. 0.07).AND.(D .LT. 0.28)) AV= 3.054*(D- 0.07)+0.011
        IF((D .GE. 0.28).AND.(D .LT. 0.86)) AV= 0.652
        IF((D .GE. 0.86).AND.(D .LT. 1.12)) AV= 3.317*(D- 0.86)+0.652
        IF((D .GE. 1.12).AND.(D .LT. 5.53)) AV= 1.514
        IF(D .GE. 5.53) AV= A0*(D- 5.53)+ 1.514
      RETURN
      ENDIF
      IF(CELL .EQ. 168) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.98)) AV= 1.523*D
        IF((D .GE. 0.98).AND.(D .LT. 3.59)) AV= 1.493
        IF(D .GE. 3.59) AV= A0*(D- 3.59)+ 1.493
      RETURN
      ENDIF
      IF(CELL .EQ. 169) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.04)) AV= 8.137*D
        IF((D .GE. 0.04).AND.(D .LT. 0.71)) AV= 0.325
        IF((D .GE. 0.71).AND.(D .LT. 0.98)) AV= 2.893*(D- 0.71)+0.325
        IF((D .GE. 0.98).AND.(D .LT. 3.87)) AV= 1.106
        IF((D .GE. 3.87).AND.(D .LT. 4.09)) AV= 5.526*(D- 3.87)+1.106
        IF(D .GE. 4.09) AV= A0*(D- 4.09)+ 2.322
      RETURN
      ENDIF
      IF(CELL .EQ. 170) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 4.263*D
        IF((D .GE. 0.11).AND.(D .LT. 0.48)) AV= 0.469
        IF((D .GE. 0.48).AND.(D .LT. 1.11)) AV= 0.946*(D- 0.48)+0.469
        IF((D .GE. 1.11).AND.(D .LT. 7.86)) AV= 1.065
        IF(D .GE. 7.86) AV= A0*(D- 7.86)+ 1.065
      RETURN
      ENDIF
      IF(CELL .EQ. 171) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.33)) AV= 4.980*D
        IF((D .GE. 0.33).AND.(D .LT. 5.70)) AV= 1.643
        IF(D .GE. 5.70) AV= A0*(D- 5.70)+ 1.643
      RETURN
      ENDIF
      IF(CELL .EQ. 172) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 4.665*D
        IF((D .GE. 0.09).AND.(D .LT. 0.60)) AV= 0.420
        IF((D .GE. 0.60).AND.(D .LT. 1.37)) AV= 1.150*(D- 0.60)+0.420
        IF((D .GE. 1.37).AND.(D .LT. 3.08)) AV= 1.306
        IF(D .GE. 3.08) AV= A0*(D- 3.08)+ 1.306
      RETURN
      ENDIF
      IF(CELL .EQ. 173) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.06)) AV= 8.503*D
        IF((D .GE. 0.06).AND.(D .LT. 0.73)) AV= 0.510
        IF((D .GE. 0.73).AND.(D .LT. 0.96)) AV= 3.769*(D- 0.73)+0.510
        IF((D .GE. 0.96).AND.(D .LT. 2.22)) AV= 1.377
        IF((D .GE. 2.22).AND.(D .LT. 2.50)) AV= 1.208*(D- 2.22)+1.377
        IF((D .GE. 2.50).AND.(D .LT. 5.12)) AV= 1.715
        IF(D .GE. 5.12) AV= A0*(D- 5.12)+ 1.715
      RETURN
      ENDIF
      IF(CELL .EQ. 174) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.03)) AV= 2.000*D
        IF((D .GE. 1.03).AND.(D .LT. 2.23)) AV= 2.060
        IF((D .GE. 2.23).AND.(D .LT. 2.53)) AV= 1.427*(D- 2.23)+2.060
        IF((D .GE. 2.53).AND.(D .LT. 3.60)) AV= 2.488
        IF(D .GE. 3.60) AV= A0*(D- 3.60)+ 2.488
      RETURN
      ENDIF
      IF(CELL .EQ. 175) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 6.000*D
        IF((D .GE. 0.10).AND.(D .LT. 2.24)) AV= 0.600
        IF((D .GE. 2.24).AND.(D .LT. 2.51)) AV= 1.662*(D- 2.24)+0.600
        IF((D .GE. 2.51).AND.(D .LT. 4.91)) AV= 1.049
        IF(D .GE. 4.91) AV= A0*(D- 4.91)+ 1.049
      RETURN
      ENDIF
      IF(CELL .EQ. 176) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 0.191*D
        IF((D .GE. 0.12).AND.(D .LT. 0.21)) AV= 7.426*(D- 0.12)+0.023
        IF((D .GE. 0.21).AND.(D .LT. 0.65)) AV= 0.691
        IF((D .GE. 0.65).AND.(D .LT. 1.09)) AV= 3.757*(D- 0.65)+0.691
        IF((D .GE. 1.09).AND.(D .LT. 4.53)) AV= 2.344
        IF((D .GE. 4.53).AND.(D .LT. 4.81)) AV= 0.956*(D- 4.53)+2.344
        IF(D .GE. 4.81) AV= A0*(D- 4.81)+ 2.612
      RETURN
      ENDIF
      IF(CELL .EQ. 177) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.18)) AV= 4.221*D
        IF((D .GE. 0.18).AND.(D .LT. 1.68)) AV= 0.760
        IF((D .GE. 1.68).AND.(D .LT. 2.03)) AV= 1.007*(D- 1.68)+0.760
        IF((D .GE. 2.03).AND.(D .LT. 5.37)) AV= 1.112
        IF(D .GE. 5.37) AV= A0*(D- 5.37)+ 1.112
      RETURN
      ENDIF
      IF(CELL .EQ. 178) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.41)) AV= 1.680*D
        IF((D .GE. 0.41).AND.(D .LT. 1.78)) AV= 0.689
        IF((D .GE. 1.78).AND.(D .LT. 2.06)) AV= 2.631*(D- 1.78)+0.689
        IF((D .GE. 2.06).AND.(D .LT. 4.31)) AV= 1.426
        IF(D .GE. 4.31) AV= A0*(D- 4.31)+ 1.426
      RETURN
      ENDIF
      IF(CELL .EQ. 179) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.28)) AV= 2.330*D
        IF((D .GE. 0.28).AND.(D .LT. 0.59)) AV= 0.652
        IF((D .GE. 0.59).AND.(D .LT. 0.94)) AV= 1.880*(D- 0.59)+0.652
        IF((D .GE. 0.94).AND.(D .LT. 1.04)) AV=10.994*(D- 0.94)+1.310
        IF((D .GE. 1.04).AND.(D .LT. 3.86)) AV= 2.409
        IF((D .GE. 3.86).AND.(D .LT. 4.24)) AV= 2.395*(D- 3.86)+2.409
        IF(D .GE. 4.24) AV= A0*(D- 4.24)+ 3.319
      RETURN
      ENDIF
      IF(CELL .EQ. 180) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.02)) AV= 8.541*D
        IF((D .GE. 0.02).AND.(D .LT. 0.12)) AV= 4.280*(D- 0.02)+0.171
        IF((D .GE. 0.12).AND.(D .LT. 0.69)) AV= 0.599
        IF((D .GE. 0.69).AND.(D .LT. 1.11)) AV= 2.914*(D- 0.69)+0.599
        IF((D .GE. 1.11).AND.(D .LT. 1.31)) AV= 8.397*(D- 1.11)+1.823
        IF((D .GE. 1.31).AND.(D .LT. 4.21)) AV= 3.502
        IF(D .GE. 4.21) AV= A0*(D- 4.21)+ 3.502
      RETURN
      ENDIF
      IF(CELL .EQ. 181) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.82)) AV= 0.000
        IF((D .GE. 0.82).AND.(D .LT. 1.27)) AV= 1.815*(D- 0.82)
        IF((D .GE. 1.27).AND.(D .LT. 3.80)) AV= 0.817
        IF((D .GE. 3.80).AND.(D .LT. 4.22)) AV= 2.046*(D- 3.80)+0.817
        IF(D .GE. 4.22) AV= A0*(D- 4.22)+ 1.676
      RETURN
      ENDIF
      IF(CELL .EQ. 182) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.38)) AV= 8.800*D
        IF((D .GE. 0.38).AND.(D .LT. 2.06)) AV= 3.344
        IF(D .GE. 2.06) AV= A0*(D- 2.06)+ 3.344
      RETURN
      ENDIF
      IF(CELL .EQ. 183) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.43)) AV= 1.385*D
        IF((D .GE. 0.43).AND.(D .LT. 0.89)) AV= 1.617*(D- 0.43)+0.596
        IF(D .GE. 0.89) AV= A0*(D- 0.89)+ 1.340
      RETURN
      ENDIF
      IF(CELL .EQ. 184) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV=10.000*D
        IF((D .GE. 0.05).AND.(D .LT. 1.21)) AV= 0.500
        IF((D .GE. 1.21).AND.(D .LT. 1.37)) AV= 4.964*(D- 1.21)+0.500
        IF((D .GE. 1.37).AND.(D .LT. 2.41)) AV= 1.294
        IF(D .GE. 2.41) AV= A0*(D- 2.41)+ 1.294
      RETURN
      ENDIF
      IF(CELL .EQ. 185) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 9.404*D
        IF((D .GE. 0.09).AND.(D .LT. 2.02)) AV= 0.837
        IF(D .GE. 2.02) AV= A0*(D- 2.02)+ 0.837
      RETURN
      ENDIF
      IF(CELL .EQ. 186) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.53)) AV= 2.126*D
        IF((D .GE. 0.53).AND.(D .LT. 3.49)) AV= 1.127
        IF(D .GE. 3.49) AV= A0*(D- 3.49)+ 1.127
      RETURN
      ENDIF
      IF(CELL .EQ. 187) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 5.225*D
        IF((D .GE. 0.10).AND.(D .LT. 0.64)) AV= 0.523
        IF((D .GE. 0.64).AND.(D .LT. 0.85)) AV= 1.052*(D- 0.64)+0.523
        IF((D .GE. 0.85).AND.(D .LT. 0.87)) AV=33.283*(D- 0.85)+0.744
        IF((D .GE. 0.87).AND.(D .LT. 1.20)) AV= 4.378*(D- 0.87)+1.410
        IF(D .GE. 1.20) AV= A0*(D- 1.20)+ 2.855
      RETURN
      ENDIF
      IF(CELL .EQ. 188) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 7.125*D
        IF((D .GE. 0.08).AND.(D .LT. 0.74)) AV= 0.570
        IF((D .GE. 0.74).AND.(D .LT. 1.49)) AV= 2.215*(D- 0.74)+0.570
        IF(D .GE. 1.49) AV= A0*(D- 1.49)+ 2.231
      RETURN
      ENDIF
      IF(CELL .EQ. 189) THEN
        IF((D .GE. 0.00).AND.(D .LT. 2.43)) AV= 1.394*D
        IF((D .GE. 2.43).AND.(D .LT. 3.49)) AV= 3.387
        IF(D .GE. 3.49) AV= A0*(D- 3.49)+ 3.387
      RETURN
      ENDIF
      IF(CELL .EQ. 190) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.20)) AV= 3.347*D
        IF((D .GE. 0.20).AND.(D .LT. 1.00)) AV= 0.669
        IF((D .GE. 1.00).AND.(D .LT. 1.20)) AV= 7.806*(D- 1.00)+0.669
        IF((D .GE. 1.20).AND.(D .LT. 5.57)) AV= 2.230
        IF(D .GE. 5.57) AV= A0*(D- 5.57)+ 2.230
      RETURN
      ENDIF
      IF(CELL .EQ. 191) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.57)) AV= 0.997*D
        IF((D .GE. 1.57).AND.(D .LT. 2.98)) AV= 1.565
        IF((D .GE. 2.98).AND.(D .LT. 3.46)) AV= 0.800*(D- 2.98)+1.565
        IF((D .GE. 3.46).AND.(D .LT. 5.09)) AV= 1.949
        IF(D .GE. 5.09) AV= A0*(D- 5.09)+ 1.949
      RETURN
      ENDIF
      IF(CELL .EQ. 192) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 2.837*D
        IF((D .GE. 0.16).AND.(D .LT. 0.92)) AV= 0.454
        IF((D .GE. 0.92).AND.(D .LT. 1.24)) AV= 2.645*(D- 0.92)+0.454
        IF((D .GE. 1.24).AND.(D .LT. 5.43)) AV= 1.300
        IF(D .GE. 5.43) AV= A0*(D- 5.43)+ 1.300
      RETURN
      ENDIF
      IF(CELL .EQ. 193) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 0.000
        IF((D .GE. 0.07).AND.(D .LT. 0.20)) AV= 3.644*(D- 0.07)
        IF((D .GE. 0.20).AND.(D .LT. 0.63)) AV= 0.474
        IF((D .GE. 0.63).AND.(D .LT. 1.00)) AV= 7.453*(D- 0.63)+0.474
        IF((D .GE. 1.00).AND.(D .LT. 1.93)) AV= 3.232
        IF(D .GE. 1.93) AV= A0*(D- 1.93)+ 3.232
      RETURN
      ENDIF
      IF(CELL .EQ. 194) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 7.870*D
        IF((D .GE. 0.10).AND.(D .LT. 0.75)) AV= 0.787
        IF((D .GE. 0.75).AND.(D .LT. 1.06)) AV= 6.391*(D- 0.75)+0.787
        IF((D .GE. 1.06).AND.(D .LT. 3.94)) AV= 2.768
        IF(D .GE. 3.94) AV= A0*(D- 3.94)+ 2.768
      RETURN
      ENDIF
      IF(CELL .EQ. 195) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.21)) AV= 1.748*D
        IF((D .GE. 0.21).AND.(D .LT. 1.15)) AV= 0.367
        IF(D .GE. 1.15) AV= A0*(D- 1.15)+ 0.367
      RETURN
      ENDIF
      IF(CELL .EQ. 196) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.78)) AV= 1.944*D
        IF((D .GE. 0.78).AND.(D .LT. 3.11)) AV= 1.516
        IF(D .GE. 3.11) AV= A0*(D- 3.11)+ 1.516
      RETURN
      ENDIF
      IF(CELL .EQ. 197) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.77)) AV= 2.797*D
        IF((D .GE. 0.77).AND.(D .LT. 1.91)) AV= 2.154
        IF(D .GE. 1.91) AV= A0*(D- 1.91)+ 2.154
      RETURN
      ENDIF
      IF(CELL .EQ. 198) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.83)) AV= 1.877*D
        IF((D .GE. 0.83).AND.(D .LT. 3.28)) AV= 1.558
        IF(D .GE. 3.28) AV= A0*(D- 3.28)+ 1.558
      RETURN
      ENDIF
      IF(CELL .EQ. 199) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.81)) AV= 0.000
        IF((D .GE. 0.81).AND.(D .LT. 1.02)) AV= 6.737*(D- 0.81)
        IF((D .GE. 1.02).AND.(D .LT. 3.00)) AV= 1.415
        IF(D .GE. 3.00) AV= A0*(D- 3.00)+ 1.415
      RETURN
      ENDIF
      IF(CELL .EQ. 200) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.02)) AV= 1.793*D
        IF((D .GE. 1.02).AND.(D .LT. 5.64)) AV= 1.829
        IF(D .GE. 5.64) AV= A0*(D- 5.64)+ 1.829
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK6(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ. 201) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.19)) AV= 3.127*D
        IF((D .GE. 0.19).AND.(D .LT. 0.62)) AV= 0.594
        IF((D .GE. 0.62).AND.(D .LT. 0.92)) AV= 5.287*(D- 0.62)+0.594
        IF((D .GE. 0.92).AND.(D .LT. 4.51)) AV= 2.180
        IF(D .GE. 4.51) AV= A0*(D- 4.51)+ 2.180
      RETURN
      ENDIF
      IF(CELL .EQ. 202) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 6.655*D
        IF((D .GE. 0.08).AND.(D .LT. 0.91)) AV= 0.532
        IF((D .GE. 0.91).AND.(D .LT. 1.28)) AV= 5.971*(D- 0.91)+0.532
        IF((D .GE. 1.28).AND.(D .LT. 3.37)) AV= 2.741
        IF(D .GE. 3.37) AV= A0*(D- 3.37)+ 2.741
      RETURN
      ENDIF
      IF(CELL .EQ. 203) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 3.517*D
        IF((D .GE. 0.12).AND.(D .LT. 0.83)) AV= 0.422
        IF((D .GE. 0.83).AND.(D .LT. 1.39)) AV= 4.980*(D- 0.83)+0.422
        IF((D .GE. 1.39).AND.(D .LT. 4.48)) AV= 3.211
        IF(D .GE. 4.48) AV= A0*(D- 4.48)+ 3.211
      RETURN
      ENDIF
      IF(CELL .EQ. 204) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.62)) AV= 0.000
        IF((D .GE. 0.62).AND.(D .LT. 1.04)) AV= 6.686*(D- 0.62)
        IF((D .GE. 1.04).AND.(D .LT. 2.57)) AV= 2.808
        IF(D .GE. 2.57) AV= A0*(D- 2.57)+ 2.808
      RETURN
      ENDIF
      IF(CELL .EQ. 205) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.82)) AV= 0.000
        IF((D .GE. 0.82).AND.(D .LT. 1.16)) AV= 5.196*(D- 0.82)
        IF((D .GE. 1.16).AND.(D .LT. 4.11)) AV= 1.767
        IF(D .GE. 4.11) AV= A0*(D- 4.11)+ 1.767
      RETURN
      ENDIF
      IF(CELL .EQ. 206) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.94)) AV= 3.076*D
        IF((D .GE. 0.94).AND.(D .LT. 3.74)) AV= 2.891
        IF(D .GE. 3.74) AV= A0*(D- 3.74)+ 2.891
      RETURN
      ENDIF
      IF(CELL .EQ. 207) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.64)) AV= 1.077*D
        IF((D .GE. 0.64).AND.(D .LT. 1.42)) AV= 1.142*(D- 0.64)+0.689
        IF(D .GE. 1.42) AV= A0*(D- 1.42)+ 1.580
      RETURN
      ENDIF
      IF(CELL .EQ. 208) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.26)) AV= 3.330*D
        IF((D .GE. 0.26).AND.(D .LT. 0.75)) AV= 0.866
        IF((D .GE. 0.75).AND.(D .LT. 1.16)) AV= 4.273*(D- 0.75)+0.866
        IF((D .GE. 1.16).AND.(D .LT. 2.96)) AV= 2.618
        IF(D .GE. 2.96) AV= A0*(D- 2.96)+ 2.618
      RETURN
      ENDIF
      IF(CELL .EQ. 209) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.22)) AV= 2.459*D
        IF((D .GE. 0.22).AND.(D .LT. 0.80)) AV= 0.541
        IF((D .GE. 0.80).AND.(D .LT. 1.11)) AV= 2.101*(D- 0.80)+0.541
        IF((D .GE. 1.11).AND.(D .LT. 1.97)) AV= 1.192
        IF((D .GE. 1.97).AND.(D .LT. 2.29)) AV= 8.367*(D- 1.97)+1.192
        IF(D .GE. 2.29) AV= A0*(D- 2.29)+ 3.869
      RETURN
      ENDIF
      IF(CELL .EQ. 210) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 3.830*D
        IF((D .GE. 0.12).AND.(D .LT. 0.82)) AV= 0.460
        IF((D .GE. 0.82).AND.(D .LT. 1.12)) AV= 5.140*(D- 0.82)+0.460
        IF((D .GE. 1.12).AND.(D .LT. 3.44)) AV= 2.002
        IF(D .GE. 3.44) AV= A0*(D- 3.44)+ 2.002
      RETURN
      ENDIF
      IF(CELL .EQ. 211) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 9.152*D
        IF((D .GE. 0.07).AND.(D .LT. 0.72)) AV= 0.641
        IF((D .GE. 0.72).AND.(D .LT. 0.96)) AV= 5.439*(D- 0.72)+0.641
        IF((D .GE. 0.96).AND.(D .LT. 3.05)) AV= 1.946
        IF(D .GE. 3.05) AV= A0*(D- 3.05)+ 1.946
      RETURN
      ENDIF
      IF(CELL .EQ. 212) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 5.467*D
        IF((D .GE. 0.09).AND.(D .LT. 0.75)) AV= 0.492
        IF((D .GE. 0.75).AND.(D .LT. 1.00)) AV= 4.268*(D- 0.75)+0.492
        IF((D .GE. 1.00).AND.(D .LT. 3.37)) AV= 1.559
        IF(D .GE. 3.37) AV= A0*(D- 3.37)+ 1.559
      RETURN
      ENDIF
      IF(CELL .EQ. 213) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 0.000
        IF((D .GE. 0.16).AND.(D .LT. 0.27)) AV= 6.521*(D- 0.16)
        IF((D .GE. 0.27).AND.(D .LT. 0.77)) AV= 0.717
        IF((D .GE. 0.77).AND.(D .LT. 0.95)) AV= 3.075*(D- 0.77)+0.717
        IF((D .GE. 0.95).AND.(D .LT. 3.94)) AV= 1.271
        IF(D .GE. 3.94) AV= A0*(D- 3.94)+ 1.271
      RETURN
      ENDIF
      IF(CELL .EQ. 214) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.19)) AV= 3.867*D
        IF((D .GE. 0.19).AND.(D .LT. 0.76)) AV= 0.735
        IF((D .GE. 0.76).AND.(D .LT. 0.88)) AV=10.612*(D- 0.76)+0.735
        IF((D .GE. 0.88).AND.(D .LT. 2.51)) AV= 2.008
        IF(D .GE. 2.51) AV= A0*(D- 2.51)+ 2.008
      RETURN
      ENDIF
      IF(CELL .EQ. 215) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.59)) AV= 0.930*D
        IF((D .GE. 0.59).AND.(D .LT. 0.68)) AV= 7.868*(D- 0.59)+0.549
        IF((D .GE. 0.68).AND.(D .LT. 0.99)) AV= 3.138*(D- 0.68)+1.257
        IF((D .GE. 0.99).AND.(D .LT. 1.61)) AV= 2.230
        IF((D .GE. 1.61).AND.(D .LT. 1.97)) AV= 3.354*(D- 1.61)+2.230
        IF(D .GE. 1.97) AV= A0*(D- 1.97)+ 3.437
      RETURN
      ENDIF
      IF(CELL .EQ. 216) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 6.912*D
        IF((D .GE. 0.14).AND.(D .LT. 0.79)) AV= 0.968
        IF((D .GE. 0.79).AND.(D .LT. 0.96)) AV= 7.904*(D- 0.79)+0.968
        IF((D .GE. 0.96).AND.(D .LT. 3.19)) AV= 2.312
        IF(D .GE. 3.19) AV= A0*(D- 3.19)+ 2.312
      RETURN
      ENDIF
      IF(CELL .EQ. 217) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.68)) AV= 0.000
        IF((D .GE. 0.68).AND.(D .LT. 1.06)) AV= 3.190*(D- 0.68)
        IF((D .GE. 1.06).AND.(D .LT. 1.92)) AV= 1.212
        IF(D .GE. 1.92) AV= A0*(D- 1.92)+ 1.212
      RETURN
      ENDIF
      IF(CELL .EQ. 218) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV= 5.000*D
        IF((D .GE. 0.10).AND.(D .LT. 0.99)) AV= 0.500
        IF(D .GE. 0.99) AV= A0*(D- 0.99)+ 0.500
      RETURN
      ENDIF
      IF(CELL .EQ. 220) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.17)) AV= 2.464*D
        IF((D .GE. 0.17).AND.(D .LT. 0.79)) AV= 0.419
        IF((D .GE. 0.79).AND.(D .LT. 1.02)) AV= 1.352*(D- 0.79)+0.419
        IF((D .GE. 1.02).AND.(D .LT. 2.80)) AV= 0.730
        IF(D .GE. 2.80) AV= A0*(D- 2.80)+ 0.730
      RETURN
      ENDIF
      IF(CELL .EQ. 221) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.84)) AV= 1.935*D
        IF((D .GE. 0.84).AND.(D .LT. 2.04)) AV= 1.625
        IF(D .GE. 2.04) AV= A0*(D- 2.04)+ 1.625
      RETURN
      ENDIF
      IF(CELL .EQ. 222) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.87)) AV= 0.926*D
        IF((D .GE. 0.87).AND.(D .LT. 1.13)) AV= 4.978*(D- 0.87)+0.806
        IF((D .GE. 1.13).AND.(D .LT. 1.26)) AV=12.680*(D- 1.13)+2.100
        IF((D .GE. 1.26).AND.(D .LT. 2.17)) AV= 3.748
        IF(D .GE. 2.17) AV= A0*(D- 2.17)+ 3.748
      RETURN
      ENDIF
      IF(CELL .EQ. 223) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.94)) AV= 0.000
        IF((D .GE. 0.94).AND.(D .LT. 1.11)) AV=11.628*(D- 0.94)
        IF((D .GE. 1.11).AND.(D .LT. 3.83)) AV= 1.977
        IF(D .GE. 3.83) AV= A0*(D- 3.83)+ 1.977
      RETURN
      ENDIF
      IF(CELL .EQ. 224) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.89)) AV= 0.000
        IF((D .GE. 0.89).AND.(D .LT. 1.15)) AV= 8.955*(D- 0.89)
        IF((D .GE. 1.15).AND.(D .LT. 2.27)) AV= 2.328
        IF(D .GE. 2.27) AV= A0*(D- 2.27)+ 2.328
      RETURN
      ENDIF
      IF(CELL .EQ. 225) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.87)) AV= 1.043*D
        IF((D .GE. 0.87).AND.(D .LT. 1.06)) AV=10.919*(D- 0.87)+0.907
        IF((D .GE. 1.06).AND.(D .LT. 2.96)) AV= 2.982
        IF(D .GE. 2.96) AV= A0*(D- 2.96)+ 2.982
      RETURN
      ENDIF
      IF(CELL .EQ. 226) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.68)) AV= 0.912*D
        IF((D .GE. 0.68).AND.(D .LT. 1.41)) AV= 0.620
        IF((D .GE. 1.41).AND.(D .LT. 1.57)) AV= 5.329*(D- 1.41)+0.620
        IF((D .GE. 1.57).AND.(D .LT. 2.54)) AV= 1.473
        IF(D .GE. 2.54) AV= A0*(D- 2.54)+ 1.473
      RETURN
      ENDIF
      IF(CELL .EQ. 227) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.10)) AV=10.000*D
        IF((D .GE. 0.10).AND.(D .LT. 1.49)) AV= 1.000
        IF((D .GE. 1.49).AND.(D .LT. 1.62)) AV= 0.935*(D- 1.49)+1.000
        IF((D .GE. 1.62).AND.(D .LT. 3.09)) AV= 1.122
        IF(D .GE. 3.09) AV= A0*(D- 3.09)+ 1.122
      RETURN
      ENDIF
      IF(CELL .EQ. 228) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.30)) AV= 2.430*D
        IF((D .GE. 0.30).AND.(D .LT. 1.95)) AV= 0.729
        IF(D .GE. 1.95) AV= A0*(D- 1.95)+ 0.729
      RETURN
      ENDIF
      IF(CELL .EQ. 229) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.11)) AV= 4.779*D
        IF((D .GE. 0.11).AND.(D .LT. 0.87)) AV= 0.526
        IF((D .GE. 0.87).AND.(D .LT. 1.24)) AV= 5.485*(D- 0.87)+0.526
        IF((D .GE. 1.24).AND.(D .LT. 3.33)) AV= 2.555
        IF(D .GE. 3.33) AV= A0*(D- 3.33)+ 2.555
      RETURN
      ENDIF
      IF(CELL .EQ. 230) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.07)) AV= 0.323*D
        IF((D .GE. 0.07).AND.(D .LT. 0.23)) AV= 1.882*(D- 0.07)+0.023
        IF((D .GE. 0.23).AND.(D .LT. 0.89)) AV= 0.324
        IF((D .GE. 0.89).AND.(D .LT. 1.16)) AV= 1.673*(D- 0.89)+0.324
        IF((D .GE. 1.16).AND.(D .LT. 3.98)) AV= 0.776
        IF(D .GE. 3.98) AV= A0*(D- 3.98)+ 0.776
      RETURN
      ENDIF
      IF(CELL .EQ. 231) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.91)) AV= 0.679*D
        IF((D .GE. 0.91).AND.(D .LT. 1.21)) AV= 3.872*(D- 0.91)+0.618
        IF((D .GE. 1.21).AND.(D .LT. 4.69)) AV= 1.780
        IF(D .GE. 4.69) AV= A0*(D- 4.69)+ 1.780
      RETURN
      ENDIF
      IF(CELL .EQ. 233) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.34)) AV= 1.319*D
        IF((D .GE. 0.34).AND.(D .LT. 0.63)) AV= 9.076*(D- 0.34)+0.448
        IF((D .GE. 0.63).AND.(D .LT. 4.70)) AV= 3.080
        IF(D .GE. 4.70) AV= A0*(D- 4.70)+ 3.080
      RETURN
      ENDIF
      IF(CELL .EQ. 235) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.82)) AV= 2.424*D
        IF((D .GE. 0.82).AND.(D .LT. 3.86)) AV= 1.998
        IF(D .GE. 3.86) AV= A0*(D- 3.86)+ 1.998
      RETURN
      ENDIF
      IF(CELL .EQ. 237) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.19)) AV= 2.724*D
        IF((D .GE. 0.19).AND.(D .LT. 1.43)) AV= 0.518
        IF((D .GE. 1.43).AND.(D .LT. 5.02)) AV= 0.105*(D- 1.43)+0.518
        IF((D .GE. 5.02).AND.(D .LT. 7.51)) AV= 0.895
        IF(D .GE. 7.51) AV= A0*(D- 7.51)+ 0.895
      RETURN
      ENDIF
      IF(CELL .EQ. 238) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.26)) AV= 0.970*D
        IF((D .GE. 1.26).AND.(D .LT. 1.96)) AV= 0.561*(D- 1.26)+1.222
        IF((D .GE. 1.96).AND.(D .LT. 5.62)) AV= 1.615
        IF(D .GE. 5.62) AV= A0*(D- 5.62)+ 1.615
      RETURN
      ENDIF
      IF(CELL .EQ. 239) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.00)) AV= 1.070*D
        IF((D .GE. 1.00).AND.(D .LT. 2.07)) AV= 1.070
        IF((D .GE. 2.07).AND.(D .LT. 2.28)) AV= 5.444*(D- 2.07)+1.070
        IF(D .GE. 2.28) AV= A0*(D- 2.28)+ 2.213
      RETURN
      ENDIF
      IF(CELL .EQ. 240) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.02)) AV= 0.857*D
        IF((D .GE. 1.02).AND.(D .LT. 2.48)) AV= 0.874
        IF(D .GE. 2.48) AV= A0*(D- 2.48)+ 0.874
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK7(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ. 241) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.98)) AV= 1.040*D
        IF((D .GE. 0.98).AND.(D .LT. 2.00)) AV= 1.019
        IF((D .GE. 2.00).AND.(D .LT. 2.24)) AV=11.576*(D- 2.00)+1.019
        IF(D .GE. 2.24) AV= A0*(D- 2.24)+ 3.797
      RETURN
      ENDIF
      IF(CELL .EQ. 242) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.71)) AV= 2.623*D
        IF((D .GE. 0.71).AND.(D .LT. 1.98)) AV= 1.862
        IF((D .GE. 1.98).AND.(D .LT. 2.21)) AV= 3.092*(D- 1.98)+1.862
        IF(D .GE. 2.21) AV= A0*(D- 2.21)+ 2.573
      RETURN
      ENDIF
      IF(CELL .EQ. 243) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.91)) AV= 1.184*D
        IF((D .GE. 0.91).AND.(D .LT. 4.74)) AV= 1.080
        IF(D .GE. 4.74) AV= A0*(D- 4.74)+ 1.080
      RETURN
      ENDIF
      IF(CELL .EQ. 244) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.06)) AV= 1.532*D
        IF((D .GE. 1.06).AND.(D .LT. 4.05)) AV= 1.624
        IF(D .GE. 4.05) AV= A0*(D- 4.05)+ 1.624
      RETURN
      ENDIF
      IF(CELL .EQ. 245) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.18)) AV= 2.584*D
        IF((D .GE. 0.18).AND.(D .LT. 0.36)) AV= 6.320*(D- 0.18)+0.465
        IF((D .GE. 0.36).AND.(D .LT. 1.75)) AV= 1.603
        IF((D .GE. 1.75).AND.(D .LT. 1.98)) AV= 3.358*(D- 1.75)+1.603
        IF(D .GE. 1.98) AV= A0*(D- 1.98)+ 2.375
      RETURN
      ENDIF
      IF(CELL .EQ. 246) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.58)) AV= 0.000
        IF((D .GE. 0.58).AND.(D .LT. 1.36)) AV= 4.712*(D- 0.58)
        IF((D .GE. 1.36).AND.(D .LT. 4.89)) AV= 3.675
        IF(D .GE. 4.89) AV= A0*(D- 4.89)+ 3.675
      RETURN
      ENDIF
      IF(CELL .EQ. 247) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.82)) AV= 0.000
        IF((D .GE. 0.82).AND.(D .LT. 2.49)) AV= 1.435*(D- 0.82)
        IF((D .GE. 2.49).AND.(D .LT. 4.28)) AV= 2.396
        IF(D .GE. 4.28) AV= A0*(D- 4.28)+ 2.396
      RETURN
      ENDIF
      IF(CELL .EQ. 248) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.09)) AV= 0.130*D
        IF((D .GE. 0.09).AND.(D .LT. 0.30)) AV= 5.872*(D- 0.09)+0.012
        IF((D .GE. 0.30).AND.(D .LT. 3.23)) AV= 1.245
        IF(D .GE. 3.23) AV= A0*(D- 3.23)+ 1.245
      RETURN
      ENDIF
      IF(CELL .EQ. 249) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.51)) AV= 1.622*D
        IF((D .GE. 1.51).AND.(D .LT. 2.95)) AV= 2.449
        IF(D .GE. 2.95) AV= A0*(D- 2.95)+ 2.449
      RETURN
      ENDIF
      IF(CELL .EQ. 250) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.48)) AV= 0.000
        IF((D .GE. 0.48).AND.(D .LT. 0.77)) AV=12.676*(D- 0.48)
        IF((D .GE. 0.77).AND.(D .LT. 3.08)) AV= 3.676
        IF(D .GE. 3.08) AV= A0*(D- 3.08)+ 3.676
      RETURN
      ENDIF
      IF(CELL .EQ. 251) THEN
        IF((D .GE. 0.00).AND.(D .LT. 3.49)) AV= 0.821*D
        IF((D .GE. 3.49).AND.(D .LT. 8.74)) AV= 2.865
        IF(D .GE. 8.74) AV= A0*(D- 8.74)+ 2.865
      RETURN
      ENDIF
      IF(CELL .EQ. 252) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.90)) AV= 0.000
        IF((D .GE. 0.90).AND.(D .LT. 1.87)) AV= 4.599*(D- 0.90)
        IF(D .GE. 1.87) AV= A0*(D- 1.87)+ 4.461
      RETURN
      ENDIF
      IF(CELL .EQ. 253) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.69)) AV= 3.591*D
        IF((D .GE. 0.69).AND.(D .LT. 5.20)) AV= 2.478
        IF(D .GE. 5.20) AV= A0*(D- 5.20)+ 2.478
      RETURN
      ENDIF
      IF(CELL .EQ. 254) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.13)) AV= 0.000
        IF((D .GE. 0.13).AND.(D .LT. 0.42)) AV= 9.201*(D- 0.13)
        IF((D .GE. 0.42).AND.(D .LT. 3.45)) AV= 2.668
        IF(D .GE. 3.45) AV= A0*(D- 3.45)+ 2.668
      RETURN
      ENDIF
      IF(CELL .EQ. 255) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.39)) AV= 7.726*D
        IF((D .GE. 0.39).AND.(D .LT. 3.15)) AV= 3.013
        IF(D .GE. 3.15) AV= A0*(D- 3.15)+ 3.013
      RETURN
      ENDIF
      IF(CELL .EQ. 256) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.24)) AV= 0.437*D
        IF((D .GE. 0.24).AND.(D .LT. 0.67)) AV= 5.991*(D- 0.24)+0.105
        IF((D .GE. 0.67).AND.(D .LT. 3.26)) AV= 2.681
        IF(D .GE. 3.26) AV= A0*(D- 3.26)+ 2.681
      RETURN
      ENDIF
      IF(CELL .EQ. 257) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.61)) AV= 5.490*D
        IF((D .GE. 0.61).AND.(D .LT. 3.64)) AV= 3.349
        IF(D .GE. 3.64) AV= A0*(D- 3.64)+ 3.349
      RETURN
      ENDIF
      IF(CELL .EQ. 258) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.34)) AV= 1.233*D
        IF((D .GE. 1.34).AND.(D .LT. 8.73)) AV= 1.652
        IF(D .GE. 8.73) AV= A0*(D- 8.73)+ 1.652
      RETURN
      ENDIF
      IF(CELL .EQ. 259) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.41)) AV= 0.000
        IF((D .GE. 0.41).AND.(D .LT. 0.73)) AV=10.210*(D- 0.41)
        IF((D .GE. 0.73).AND.(D .LT. 4.51)) AV= 3.267
        IF(D .GE. 4.51) AV= A0*(D- 4.51)+ 3.267
      RETURN
      ENDIF
      IF(CELL .EQ. 260) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.36)) AV= 0.313*D
        IF((D .GE. 0.36).AND.(D .LT. 0.58)) AV=23.712*(D- 0.36)+0.113
        IF(D .GE. 0.58) AV= A0*(D- 0.58)+ 5.330
      RETURN
      ENDIF
      IF(CELL .EQ. 261) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.68)) AV= 1.511*D
        IF((D .GE. 1.68).AND.(D .LT. 3.08)) AV= 2.538
        IF(D .GE. 3.08) AV= A0*(D- 3.08)+ 2.538
      RETURN
      ENDIF
      IF(CELL .EQ. 262) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.89)) AV= 0.980*D
        IF(D .GE. 1.89) AV= A0*(D- 1.89)+ 1.852
      RETURN
      ENDIF
      IF(CELL .EQ. 263) THEN
        IF((D .GE. 0.00).AND.(D .LT. 2.38)) AV= 1.345*D
        IF((D .GE. 2.38).AND.(D .LT. 3.95)) AV= 3.201
        IF(D .GE. 3.95) AV= A0*(D- 3.95)+ 3.201
      RETURN
      ENDIF
      IF(CELL .EQ. 264) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.50)) AV= 5.947*(D- 0.08)
        IF((D .GE. 0.50).AND.(D .LT. 5.94)) AV= 2.498
        IF(D .GE. 5.94) AV= A0*(D- 5.94)+ 2.498
      RETURN
      ENDIF
      IF(CELL .EQ. 265) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 0.000
        IF((D .GE. 0.16).AND.(D .LT. 0.47)) AV= 2.802*(D- 0.16)
        IF((D .GE. 0.47).AND.(D .LT. 4.84)) AV= 0.888
        IF(D .GE. 4.84) AV= A0*(D- 4.84)+ 0.888
      RETURN
      ENDIF
      IF(CELL .EQ. 266) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.78)) AV= 1.916*D
        IF((D .GE. 0.78).AND.(D .LT. 1.00)) AV= 7.052*(D- 0.78)+1.494
        IF(D .GE. 1.00) AV= A0*(D- 1.00)+ 3.045
      RETURN
      ENDIF
      IF(CELL .EQ. 267) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.29)) AV= 1.232*D
        IF((D .GE. 0.29).AND.(D .LT. 0.59)) AV= 3.591*(D- 0.29)+0.357
        IF((D .GE. 0.59).AND.(D .LT. 4.59)) AV= 1.434
        IF(D .GE. 4.59) AV= A0*(D- 4.59)+ 1.434
      RETURN
      ENDIF
      IF(CELL .EQ. 268) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.18)) AV= 0.000
        IF((D .GE. 0.18).AND.(D .LT. 0.49)) AV= 2.136*(D- 0.18)
        IF((D .GE. 0.49).AND.(D .LT. 3.09)) AV= 0.662
        IF(D .GE. 3.09) AV= A0*(D- 3.09)+ 0.662
      RETURN
      ENDIF
      IF(CELL .EQ. 269) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.43)) AV= 6.960*D
        IF((D .GE. 0.43).AND.(D .LT. 3.39)) AV= 2.993
        IF(D .GE. 3.39) AV= A0*(D- 3.39)+ 2.993
      RETURN
      ENDIF
      IF(CELL .EQ. 270) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.40)) AV= 1.669*(D- 0.08)
        IF((D .GE. 0.40).AND.(D .LT. 0.49)) AV=16.060*(D- 0.40)+0.534
        IF((D .GE. 0.49).AND.(D .LT. 2.36)) AV= 1.979
        IF(D .GE. 2.36) AV= A0*(D- 2.36)+ 1.979
      RETURN
      ENDIF
      IF(CELL .EQ. 271) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 0.000
        IF((D .GE. 0.15).AND.(D .LT. 0.28)) AV= 2.319*(D- 0.15)
        IF((D .GE. 0.28).AND.(D .LT. 0.69)) AV= 0.283*(D- 0.28)+0.311
        IF((D .GE. 0.69).AND.(D .LT. 1.01)) AV= 9.168*(D- 0.69)+0.425
        IF((D .GE. 1.01).AND.(D .LT. 4.39)) AV= 3.396
        IF(D .GE. 4.39) AV= A0*(D- 4.39)+ 3.396
      RETURN
      ENDIF
      IF(CELL .EQ. 272) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.97)) AV= 0.365*D
        IF((D .GE. 0.97).AND.(D .LT. 1.37)) AV= 5.972*(D- 0.97)+0.354
        IF(D .GE. 1.37) AV= A0*(D- 1.37)+ 2.743
      RETURN
      ENDIF
      IF(CELL .EQ. 273) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.49)) AV= 0.508*D
        IF((D .GE. 0.49).AND.(D .LT. 0.74)) AV= 3.288*(D- 0.49)+0.249
        IF((D .GE. 0.74).AND.(D .LT. 2.39)) AV= 1.071
        IF((D .GE. 2.39).AND.(D .LT. 2.65)) AV= 4.622*(D- 2.39)+1.071
        IF((D .GE. 2.65).AND.(D .LT. 5.55)) AV= 2.273
        IF(D .GE. 5.55) AV= A0*(D- 5.55)+ 2.273
      RETURN
      ENDIF
      IF(CELL .EQ. 274) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.22)) AV= 0.915*D
        IF((D .GE. 0.22).AND.(D .LT. 0.99)) AV= 0.201
        IF(D .GE. 0.99) AV= A0*(D- 0.99)+ 0.201
      RETURN
      ENDIF
      IF(CELL .EQ. 275) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.34)) AV= 0.000
        IF((D .GE. 0.34).AND.(D .LT. 0.38)) AV= 5.380*(D- 0.34)
        IF((D .GE. 0.38).AND.(D .LT. 0.84)) AV= 0.215
        IF((D .GE. 0.84).AND.(D .LT. 1.09)) AV= 7.339*(D- 0.84)+0.215
        IF(D .GE. 1.09) AV= A0*(D- 1.09)+ 2.050
      RETURN
      ENDIF
      IF(CELL .EQ. 276) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.25)) AV= 0.000
        IF((D .GE. 0.25).AND.(D .LT. 0.31)) AV= 2.793*(D- 0.25)
        IF((D .GE. 0.31).AND.(D .LT. 3.06)) AV= 0.302*(D- 0.31)+0.168
        IF(D .GE. 3.06) AV= A0*(D- 3.06)+ 0.999
      RETURN
      ENDIF
      IF(CELL .EQ. 277) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.24)) AV= 0.000
        IF((D .GE. 0.24).AND.(D .LT. 0.53)) AV= 0.517*(D- 0.24)
        IF((D .GE. 0.53).AND.(D .LT. 1.02)) AV= 0.653*(D- 0.53)+0.150
        IF((D .GE. 1.02).AND.(D .LT. 2.54)) AV= 1.209*(D- 1.02)+0.470
        IF((D .GE. 2.54).AND.(D .LT. 7.87)) AV= 0.271*(D- 2.54)+2.308
        IF(D .GE. 7.87) AV= A0*(D- 7.87)+ 3.752
      RETURN
      ENDIF
      IF(CELL .EQ. 278) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 4.935*D
        IF((D .GE. 0.05).AND.(D .LT. 2.02)) AV= 0.247
        IF(D .GE. 2.02) AV= A0*(D- 2.02)+ 0.247
      RETURN
      ENDIF
      IF(CELL .EQ. 279) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.04)) AV= 6.525*D
        IF((D .GE. 0.04).AND.(D .LT. 0.68)) AV= 0.261
        IF((D .GE. 0.68).AND.(D .LT. 0.94)) AV= 9.719*(D- 0.68)+0.261
        IF((D .GE. 0.94).AND.(D .LT. 2.76)) AV= 2.788
        IF(D .GE. 2.76) AV= A0*(D- 2.76)+ 2.788
      RETURN
      ENDIF
      IF(CELL .EQ. 280) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.01)) AV= 1.568*D
        IF((D .GE. 1.01).AND.(D .LT. 3.73)) AV= 1.584
        IF(D .GE. 3.73) AV= A0*(D- 3.73)+ 1.584
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE ENK8(CELL,D,AV,A0)
      IMPLICIT REAL*4 (A-H,L-Z)

      INTEGER CELL

      IF(CELL .EQ. 281) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.30)) AV= 1.839*D
        IF((D .GE. 0.30).AND.(D .LT. 3.65)) AV= 0.552
        IF(D .GE. 3.65) AV= A0*(D- 3.65)+ 0.552
      RETURN
      ENDIF
      IF(CELL .EQ. 282) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.27)) AV= 0.611*D
        IF((D .GE. 0.27).AND.(D .LT. 1.00)) AV= 1.865*(D- 0.27)+0.165
        IF((D .GE. 1.00).AND.(D .LT. 4.77)) AV= 1.526
        IF(D .GE. 4.77) AV= A0*(D- 4.77)+ 1.526
      RETURN
      ENDIF
      IF(CELL .EQ. 283) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.43)) AV= 1.168*(D- 0.08)
        IF((D .GE. 0.43).AND.(D .LT. 1.12)) AV= 0.409
        IF((D .GE. 1.12).AND.(D .LT. 1.32)) AV= 7.655*(D- 1.12)+0.409
        IF((D .GE. 1.32).AND.(D .LT. 4.76)) AV= 1.940
        IF(D .GE. 4.76) AV= A0*(D- 4.76)+ 1.940
      RETURN
      ENDIF
      IF(CELL .EQ. 284) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 2.088*D
        IF((D .GE. 0.05).AND.(D .LT. 0.46)) AV= 0.104
        IF((D .GE. 0.46).AND.(D .LT. 0.72)) AV= 5.550*(D- 0.46)+0.104
        IF((D .GE. 0.72).AND.(D .LT. 3.89)) AV= 1.547
        IF(D .GE. 3.89) AV= A0*(D- 3.89)+ 1.547
      RETURN
      ENDIF
      IF(CELL .EQ. 285) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.12)) AV= 1.955*D
        IF((D .GE. 0.12).AND.(D .LT. 0.69)) AV= 0.235
        IF((D .GE. 0.69).AND.(D .LT. 0.97)) AV= 2.475*(D- 0.69)+0.235
        IF((D .GE. 0.97).AND.(D .LT. 3.99)) AV= 0.928
        IF(D .GE. 3.99) AV= A0*(D- 3.99)+ 0.928
      RETURN
      ENDIF
      IF(CELL .EQ. 286) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.02)) AV= 4.314*D
        IF((D .GE. 0.02).AND.(D .LT. 0.72)) AV= 0.086
        IF((D .GE. 0.72).AND.(D .LT. 1.00)) AV= 8.890*(D- 0.72)+0.086
        IF((D .GE. 1.00).AND.(D .LT. 2.55)) AV= 2.575
        IF(D .GE. 2.55) AV= A0*(D- 2.55)+ 2.575
      RETURN
      ENDIF
      IF(CELL .EQ. 287) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.03)) AV= 6.645*D
        IF((D .GE. 0.03).AND.(D .LT. 0.49)) AV= 0.199
        IF((D .GE. 0.49).AND.(D .LT. 0.75)) AV= 8.076*(D- 0.49)+0.199
        IF((D .GE. 0.75).AND.(D .LT. 3.36)) AV= 2.299
        IF(D .GE. 3.36) AV= A0*(D- 3.36)+ 2.299
      RETURN
      ENDIF
      IF(CELL .EQ. 288) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.15)) AV= 1.952*D
        IF((D .GE. 0.15).AND.(D .LT. 1.01)) AV= 0.293
        IF((D .GE. 1.01).AND.(D .LT. 1.52)) AV= 2.407*(D- 1.01)+0.293
        IF((D .GE. 1.52).AND.(D .LT. 3.54)) AV= 1.521
        IF(D .GE. 3.54) AV= A0*(D- 3.54)+ 1.521
      RETURN
      ENDIF
      IF(CELL .EQ. 289) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.58)) AV= 0.386*D
        IF((D .GE. 0.58).AND.(D .LT. 0.87)) AV= 5.021*(D- 0.58)+0.224
        IF((D .GE. 0.87).AND.(D .LT. 3.72)) AV= 1.680
        IF(D .GE. 3.72) AV= A0*(D- 3.72)+ 1.680
      RETURN
      ENDIF
      IF(CELL .EQ. 290) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.78)) AV= 0.491*D
        IF((D .GE. 0.78).AND.(D .LT. 1.05)) AV=10.007*(D- 0.78)+0.383
        IF(D .GE. 1.05) AV= A0*(D- 1.05)+ 3.085
      RETURN
      ENDIF
      IF(CELL .EQ. 291) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 3.594*D
        IF((D .GE. 0.05).AND.(D .LT. 0.76)) AV= 0.180
        IF((D .GE. 0.76).AND.(D .LT. 1.02)) AV= 9.738*(D- 0.76)+0.180
        IF(D .GE. 1.02) AV= A0*(D- 1.02)+ 2.712
      RETURN
      ENDIF
      IF(CELL .EQ. 292) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.37)) AV= 0.068*D
        IF((D .GE. 0.37).AND.(D .LT. 1.24)) AV= 6.571*(D- 0.37)+0.025
        IF(D .GE. 1.24) AV= A0*(D- 1.24)+ 5.742
      RETURN
      ENDIF
      IF(CELL .EQ. 293) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.04)) AV= 6.820*D
        IF((D .GE. 0.04).AND.(D .LT. 0.61)) AV= 0.280
        IF((D .GE. 0.61).AND.(D .LT. 1.00)) AV= 8.180*(D- 0.61)+0.280
        IF(D .GE. 1.00) AV= A0*(D- 1.00)+ 3.486
      RETURN
      ENDIF
      IF(CELL .EQ. 294) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.05)) AV= 6.842*D
        IF((D .GE. 0.05).AND.(D .LT. 1.69)) AV= 0.342
        IF(D .GE. 1.69) AV= A0*(D- 1.69)+ 0.342
      RETURN
      ENDIF
      IF(CELL .EQ. 295) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.03)) AV= 4.690*D
        IF((D .GE. 0.03).AND.(D .LT. 0.57)) AV= 0.141
        IF((D .GE. 0.57).AND.(D .LT. 0.82)) AV=11.186*(D- 0.57)+0.141
        IF((D .GE. 0.82).AND.(D .LT. 2.49)) AV= 2.938
        IF(D .GE. 2.49) AV= A0*(D- 2.49)+ 2.938
      RETURN
      ENDIF
      IF(CELL .EQ. 296) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.22)) AV= 1.262*D
        IF((D .GE. 1.22).AND.(D .LT. 3.52)) AV= 1.540
        IF(D .GE. 3.52) AV= A0*(D- 3.52)+ 1.540
      RETURN
      ENDIF
      IF(CELL .EQ. 297) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.41)) AV= 0.398*D
        IF((D .GE. 0.41).AND.(D .LT. 0.63)) AV= 2.182*(D- 0.41)+0.163
        IF((D .GE. 0.63).AND.(D .LT. 1.73)) AV= 0.643
        IF(D .GE. 1.73) AV= A0*(D- 1.73)+ 0.643
      RETURN
      ENDIF
      IF(CELL .EQ. 298) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.14)) AV= 0.000
        IF((D .GE. 0.14).AND.(D .LT. 0.57)) AV= 1.592*(D- 0.14)
        IF((D .GE. 0.57).AND.(D .LT. 0.84)) AV= 0.685
        IF((D .GE. 0.84).AND.(D .LT. 1.09)) AV= 5.932*(D- 0.84)+0.685
        IF((D .GE. 1.09).AND.(D .LT. 3.14)) AV= 2.168
        IF(D .GE. 3.14) AV= A0*(D- 3.14)+ 2.168
      RETURN
      ENDIF
      IF(CELL .EQ. 299) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.18)) AV= 0.000
        IF((D .GE. 0.18).AND.(D .LT. 0.48)) AV= 7.556*(D- 0.18)
        IF((D .GE. 0.48).AND.(D .LT. 3.01)) AV= 2.267
        IF(D .GE. 3.01) AV= A0*(D- 3.01)+ 2.267
      RETURN
      ENDIF
      IF(CELL .EQ. 300) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.01)) AV= 0.734*D
        IF((D .GE. 1.01).AND.(D .LT. 1.26)) AV= 6.702*(D- 1.01)+0.741
        IF((D .GE. 1.26).AND.(D .LT. 3.67)) AV= 2.417
        IF(D .GE. 3.67) AV= A0*(D- 3.67)+ 2.417
      RETURN
      ENDIF
      IF(CELL .EQ. 301) THEN
        IF((D .GE. 0.00).AND.(D .LT. 1.00)) AV= 0.769*D
        IF((D .GE. 1.00).AND.(D .LT. 3.44)) AV= 0.769
        IF(D .GE. 3.44) AV= A0*(D- 3.44)+ 0.769
      RETURN
      ENDIF
      IF(CELL .EQ. 302) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.32)) AV= 0.076*D
        IF((D .GE. 0.32).AND.(D .LT. 0.89)) AV= 1.959*(D- 0.32)+0.024
        IF((D .GE. 0.89).AND.(D .LT. 1.01)) AV= 3.599*(D- 0.89)+1.141
        IF((D .GE. 1.01).AND.(D .LT. 5.89)) AV= 1.576
        IF(D .GE. 5.89) AV= A0*(D- 5.89)+ 1.576
      RETURN
      ENDIF
      IF(CELL .EQ. 303) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.83)) AV= 1.151*D
        IF((D .GE. 0.83).AND.(D .LT. 5.06)) AV= 0.955
        IF(D .GE. 5.06) AV= A0*(D- 5.06)+ 0.955
      RETURN
      ENDIF
      IF(CELL .EQ. 304) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.26)) AV= 0.000
        IF((D .GE. 0.26).AND.(D .LT. 0.88)) AV= 2.587*(D- 0.26)
        IF((D .GE. 0.88).AND.(D .LT. 5.48)) AV= 1.604
        IF(D .GE. 5.48) AV= A0*(D- 5.48)+ 1.604
      RETURN
      ENDIF
      IF(CELL .EQ. 305) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.27)) AV= 0.000
        IF((D .GE. 0.27).AND.(D .LT. 0.56)) AV= 7.742*(D- 0.27)
        IF((D .GE. 0.56).AND.(D .LT. 4.93)) AV= 2.245
        IF(D .GE. 4.93) AV= A0*(D- 4.93)+ 2.245
      RETURN
      ENDIF
      IF(CELL .EQ. 306) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.25)) AV= 1.186*D
        IF((D .GE. 0.25).AND.(D .LT. 0.59)) AV= 4.434*(D- 0.25)+0.297
        IF((D .GE. 0.59).AND.(D .LT. 3.06)) AV= 1.805
        IF(D .GE. 3.06) AV= A0*(D- 3.06)+ 1.805
      RETURN
      ENDIF
      IF(CELL .EQ. 307) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.38)) AV= 1.895*D
        IF((D .GE. 0.38).AND.(D .LT. 0.54)) AV=14.082*(D- 0.38)+0.722
        IF((D .GE. 0.54).AND.(D .LT. 3.41)) AV= 2.975
        IF(D .GE. 3.41) AV= A0*(D- 3.41)+ 2.975
      RETURN
      ENDIF
      IF(CELL .EQ. 308) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.26)) AV= 0.000
        IF((D .GE. 0.26).AND.(D .LT. 0.53)) AV= 7.162*(D- 0.26)
        IF((D .GE. 0.53).AND.(D .LT. 1.68)) AV= 1.934
        IF(D .GE. 1.68) AV= A0*(D- 1.68)+ 1.934
      RETURN
      ENDIF
      IF(CELL .EQ. 309) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.41)) AV= 3.801*D
        IF((D .GE. 0.41).AND.(D .LT. 2.58)) AV= 1.558
        IF(D .GE. 2.58) AV= A0*(D- 2.58)+ 1.558
      RETURN
      ENDIF
      IF(CELL .EQ. 310) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.48)) AV= 3.332*D
        IF((D .GE. 0.48).AND.(D .LT. 1.94)) AV= 1.599
        IF((D .GE. 1.94).AND.(D .LT. 2.26)) AV= 4.631*(D- 1.94)+1.599
        IF((D .GE. 2.26).AND.(D .LT. 4.00)) AV= 3.081
        IF(D .GE. 4.00) AV= A0*(D- 4.00)+ 3.081
      RETURN
      ENDIF
      IF(CELL .EQ. 311) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.53)) AV= 0.262*D
        IF((D .GE. 0.53).AND.(D .LT. 0.79)) AV= 6.686*(D- 0.53)+0.139
        IF(D .GE. 0.79) AV= A0*(D- 0.79)+ 1.877
      RETURN
      ENDIF
      IF(CELL .EQ. 312) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.32)) AV= 3.565*D
        IF((D .GE. 0.32).AND.(D .LT. 1.68)) AV= 1.141
        IF(D .GE. 1.68) AV= A0*(D- 1.68)+ 1.141
      RETURN
      ENDIF
      IF(CELL .EQ. 313) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.16)) AV= 0.962*D
        IF((D .GE. 0.16).AND.(D .LT. 0.33)) AV= 3.097*(D- 0.16)+0.154
        IF((D .GE. 0.33).AND.(D .LT. 1.66)) AV= 0.680
        IF((D .GE. 1.66).AND.(D .LT. 1.93)) AV= 5.390*(D- 1.66)+0.680
        IF((D .GE. 1.93).AND.(D .LT. 5.39)) AV= 2.135
        IF(D .GE. 5.39) AV= A0*(D- 5.39)+ 2.135
      RETURN
      ENDIF
      IF(CELL .EQ. 314) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.28)) AV= 0.832*D
        IF((D .GE. 0.28).AND.(D .LT. 0.52)) AV= 4.754*(D- 0.28)+0.233
        IF((D .GE. 0.52).AND.(D .LT. 1.75)) AV= 1.374
        IF((D .GE. 1.75).AND.(D .LT. 2.02)) AV= 7.029*(D- 1.75)+1.374
        IF((D .GE. 2.02).AND.(D .LT. 2.87)) AV= 3.272
        IF(D .GE. 2.87) AV= A0*(D- 2.87)+ 3.272
      RETURN
      ENDIF
      IF(CELL .EQ. 315) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.32)) AV= 0.758*D
        IF((D .GE. 0.32).AND.(D .LT. 0.55)) AV= 7.212*(D- 0.32)+0.243
        IF((D .GE. 0.55).AND.(D .LT. 5.18)) AV= 1.902
        IF(D .GE. 5.18) AV= A0*(D- 5.18)+ 1.902
      RETURN
      ENDIF
      IF(CELL .EQ. 316) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.51)) AV= 1.921*D
        IF((D .GE. 0.51).AND.(D .LT. 0.74)) AV= 4.306*(D- 0.51)+0.980
        IF((D .GE. 0.74).AND.(D .LT. 1.17)) AV= 2.712*(D- 0.74)+1.970
        IF(D .GE. 1.17) AV= A0*(D- 1.17)+ 3.136
      RETURN
      ENDIF
      IF(CELL .EQ. 317) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.28)) AV= 0.716*D
        IF((D .GE. 0.28).AND.(D .LT. 0.30)) AV= 8.567*(D- 0.28)+0.200
        IF((D .GE. 0.30).AND.(D .LT. 0.51)) AV= 0.124*(D- 0.30)+0.371
        IF((D .GE. 0.51).AND.(D .LT. 0.76)) AV= 8.271*(D- 0.51)+0.397
        IF((D .GE. 0.76).AND.(D .LT. 3.72)) AV= 2.465
        IF(D .GE. 3.72) AV= A0*(D- 3.72)+ 2.498
      RETURN
      ENDIF
      IF(CELL .EQ. 318) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.19)) AV= 1.593*D
        IF((D .GE. 0.19).AND.(D .LT. 0.64)) AV= 6.771*(D- 0.19)+0.301
        IF((D .GE. 0.64).AND.(D .LT. 2.00)) AV= 3.348
        IF(D .GE. 2.00) AV= A0*(D- 2.00)+ 3.348
      RETURN
      ENDIF
      IF(CELL .EQ. 319) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.60)) AV= 0.634*D
        IF((D .GE. 0.60).AND.(D .LT. 1.68)) AV= 0.380
        IF(D .GE. 1.68) AV= A0*(D- 1.68)+ 0.380
      RETURN
      ENDIF
      IF(CELL .EQ. 320) THEN
        IF((D .GE. 0.00).AND.(D .LT. 0.08)) AV= 0.000
        IF((D .GE. 0.08).AND.(D .LT. 0.18)) AV= 1.042*(D- 0.08)
        IF((D .GE. 0.18).AND.(D .LT. 0.38)) AV=11.668*(D- 0.18)+0.104
        IF((D .GE. 0.38).AND.(D .LT. 2.49)) AV= 2.438
        IF((D .GE. 2.49).AND.(D .LT. 2.73)) AV= 2.715*(D- 2.49)+2.438
        IF((D .GE. 2.73).AND.(D .LT. 5.04)) AV= 3.090
        IF(D .GE. 5.04) AV= A0*(D- 5.04)+ 3.090
      RETURN
      ENDIF

      RETURN
      END

      SUBROUTINE BERDNIKOV(L,B,D,AV,SAV,A)
C-----------------------------------------------------------------------
C     From Berdnikov and Pavlovskaya, Sov. Astron. Lett. 17, 215 (1990).
C-----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,L-Z)
C-----------------------------------------------------------------------
C     BETA= parameter characterizing the distribution over z       
C     A0 = parameter characterizing the extinction in the Galactic 
C         plane by Parenago, P. P. (1945). Astron. Zh. 22, 129.   
C     NOTE: abs(b) > 10 degrees            
C-----------------------------------------------------------------------
      AV=-99.
      SAV=-99.
      PHI=64.
      BETA=0.17
C-----------------------------------------------------------------------
C
C     A0 defined in this subroutine to be 1.09 mag/kpc
C
C-----------------------------------------------------------------------
      A0=1.09
C-----------------------------------------------------------------------
      B0=0.49
      IF (ABS(B) .EQ. 90.) L=115.
      IF ((L .GE. 65.) .AND. (L .LE. 165.) 
     &   .AND. (ABS(B) .GT. 10.)) THEN
       IF ((D .LE. 1.5)) THEN
            ABSZ=D*ABS(SIN(B*3.141592/180.))
            AVA=B0*SIN((L-PHI)*3.141592/180.)+A0
            AVB=1.-EXP(-(ABSZ/BETA))
            AVC=BETA/(ABS(SIN(B*3.141592/180.)))
            AV=AVA*AVB*AVC
            SAV=0.41
         ELSE
          ABSZ=1.5*ABS(SIN(B*3.141592/180.))
          AVA=B0*SIN((L-PHI)*3.141592/180.)+A0
          AVB=1.-EXP(-(ABSZ/BETA))
          AVC=BETA/(ABS(SIN(B*3.141592/180.)))
            AVFAR=BERDFCN(1.09,1.5,D,B,BETA)
            AV=AVA*AVB*AVC+AVFAR
          SAV=SQRT(0.1681+AVFAR*AVFAR*(D-1.5)*(D-1.5)/9.)
       ENDIF
      ENDIF
      RETURN
      END

      FUNCTION BERDFCN(S,D0,D,B,BETA)
      IMPLICIT REAL*4 (A-H,L-Z)
C---------------------------------------------------------------------------
C     This function calculates the extinction Av for stars beyond d=1.5 kpc,
C     and is used in BERDNIKOV & PAVLOVSKAYA. 
C     VARIABLE DEFINITIONS: S is the differential Galactic plane extinction,
C     D0 is 1.5 kpc, D is the distance for which extinction needs determined,
C     b is galactic latitude, and beta is the half width of the material 
C     perpendicular to the galactic plane.
C---------------------------------------------------------------------------
      E0=EXP(-D0*ABS(SIN(B*3.141592/180.))/BETA)
      E=EXP(-D*ABS(SIN(B*3.141592/180.))/BETA)
      IF (ABS(B) .GT. 89.999) THEN
         BERDFCN=0.
         RETURN
      ENDIF
      BERDFCN=BETA*S*(E0-E)/ABS(SIN(B*3.141592/180.))
      RETURN
      END

      SUBROUTINE ARENOU_ETAL(L,B,D,AV,SAV,A0)
C-----------------------------------------------------------------------
C     From Arenou et al., A&A 258, 104 (1992).
C-----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,L-Z)
C-----------------------------------------------------------------------
C     ALPHA, BETA, and GAMMA are the regression coefficients the authors 
C       obtained by using the least square method.      
C     R = Value where the plotted curve becomes linear.          
C     SAV4 = the relative error in percentages.                  
C-----------------------------------------------------------------------
      GAMMA=0.0    !  REASSIGN ONLY IF NON-ZERO
      IF ((B .GE. -90.) .AND. (B .LT. -60.)) THEN
         IF ((L .GE. 0.) .AND. (L .LT. 29.)) THEN
            ALPHA=2.22534
            BETA=-6.00212
            R=.052
          SAV4=13.
         ENDIF
         IF ((L .GE. 29.) .AND. (L .LT. 57.)) THEN
            ALPHA=3.35436
            BETA=-14.74567
            R=.035
          SAV4=40.
         ENDIF
         IF ((L .GE. 57.) .AND. (L .LT. 85.)) THEN
            ALPHA=2.77637
            BETA=-9.62706
            R=.042
          SAV4=15.
         ENDIF
         IF ((L .GE. 85.0) .AND. (L .LT. 110.0)) THEN
           ALPHA=4.44190
           BETA=-19.92097
           R=0.025
          SAV4=36.
         ENDIF
         IF ((L .GE. 110.0) .AND. (L .LT. 150.0)) THEN
           ALPHA=4.46685
           BETA=-26.07305
           R=0.026
         SAV4=28.
         ENDIF
         IF ((L .GE. 150.0) .AND. (L .LT. 180.0)) THEN
           ALPHA=7.63699
           BETA=-46.10856
           R=0.014
         SAV4=38.
         ENDIF
         IF ((L .GE. 180.0) .AND. (L .LT. 210.0)) THEN
           ALPHA=2.43412
           BETA=-8.69913
           R=0.050
         SAV4=36.
         ENDIF
         IF ((L .GE. 210.0) .AND. (L .LT. 240.0)) THEN
           ALPHA=3.34481
           BETA=-13.93228
           R=0.035
         SAV4=38.
         ENDIF
         IF ((L .GE. 240.0) .AND. (L .LT. 270.0)) THEN
           ALPHA=1.40733
           BETA=-3.43418
           R=0.091
         SAV4=30.
         ENDIF
         IF ((L .GE. 270.0) .AND. (L .LT. 300.0)) THEN
           ALPHA=1.64466
           BETA=-3.97380
           R=0.074
         SAV4=28.
         ENDIF
         IF ((L .GE. 300.0) .AND. (L .LT. 330.0)) THEN
           ALPHA=2.12696
           BETA=-6.05682
           R=0.056
         SAV4=14.
         ENDIF
         IF ((L .GE. 330.0) .AND. (L .LT. 360.0)) THEN
           ALPHA=2.34636
           BETA=-8.17784
           R=0.052
         SAV4=16.
         ENDIF
      ENDIF
      IF ((B .GE. -60.) .AND. (B .LT. -45.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 30.0)) THEN
          ALPHA=2.77060
          BETA=-9.52310
          R=0.145
        SAV4=16.
        ENDIF
        IF ((L .GE.30.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.96533
          BETA=-5.63251
          R=0.174
          SAV4=6.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 110)) THEN
          ALPHA=1.93622
          BETA=-13.31757
          R=0.073
        SAV4=26.
        ENDIF
        IF ((L .GE. 110.0) .AND. (L .LT. 180.0)) THEN
          ALPHA=1.05414
          BETA=-2.36540
          R=0.223
          SAV4=74.
        ENDIF
        IF ((L .GE.180.0) .AND. (L .LT. 210.0)) THEN
          ALPHA=1.39990
          BETA=-1.35325
          R=0.252
        SAV4=10.
        ENDIF
        IF ((L .GE. 210.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=2.73481
          BETA=-11.70266
          R=0.117
        SAV4=8.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 270.0)) THEN
          ALPHA=2.99784
          BETA=-11.64272
          R=0.129
          SAV4=3.
        ENDIF
        IF ((L .GE. 270.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=3.23802
          BETA=-11.63810
          R=0.139
        SAV4=7.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=1.72679
          BETA=-6.05085
          R=0.143
        SAV4=7.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=1.88890
          BETA=-5.51861
          R=0.171
        SAV4=14.
        ENDIF
      ENDIF
      IF ((B .GE. -45.) .AND. (B .LT. -30.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 30.0)) THEN
          ALPHA=1.98973
          BETA=-4.86206
          R=0.205
        SAV4=6.
        ENDIF
        IF ((L .GE. 30.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.49901
          BETA=-3.75837
          R=0.199
        SAV4=16.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 90.0)) THEN
          ALPHA=0.90091
          BETA=-1.30459
          R=0.329
        SAV4=73.
        ENDIF
        IF ((L .GE. 90.0) .AND. (L .LT. 120.0)) THEN
          ALPHA=1.94200
          BETA=-6.26833
          R=0.155
        SAV4=18.
        ENDIF
        IF ((L .GE. 120.0) .AND. (L .LT. 160.0)) THEN
          ALPHA=-0.37804
          BETA=10.77372
          R=0.210
        SAV4=100.
        ENDIF
        IF  ((L .GE. 160.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=-0.15710
          BETA=5.03190
          R=0.294
        SAV4=24.
        ENDIF
        IF  ((L .GE. 200.0) .AND. (L .LT. 235.0)) THEN
          ALPHA=3.20162
          BETA=-10.59297
          R=0.151
        SAV4=9.
        ENDIF
        IF ((L .GE. 235.0) .AND. (L .LT. 265.0)) THEN
          ALPHA=1.95079
          BETA=-4.73280
          R=0.206
        SAV4=21.
        ENDIF
        IF ((L .GE. 265.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=1.91200
          BETA=-4.97640
          R=0.192
        SAV4=17.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=2.50487
          BETA=-8.63106
          R=0.145
        SAV4=28.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=2.44394
          BETA=-9.17612
          R=0.133
        SAV4=7.
        ENDIF
      ENDIF
      IF ((B .GE. -30.) .AND. (B .LT. -15.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 20.0)) THEN
          ALPHA=2.82440
          BETA=-4.78217
          R=0.295
        SAV4=32.
        ENDIF
        IF ((L .GE. 20.0) .AND. (L .LT. 40.0)) THEN
          ALPHA=3.84362
          BETA=-8.04690
          R=0.239
        SAV4=46.
        ENDIF
        IF ((L .GE. 40.0) .AND. (L .LT. 80.0)) THEN
          ALPHA=0.60365
          BETA=0.07893
          R=0.523
        SAV4=22.
        ENDIF
        IF ((L .GE. 80.0) .AND. (L .LT. 100.0)) THEN
          ALPHA=0.58307
          BETA=-0.21053
          R=0.523
        SAV4=53.
        ENDIF
        IF ((L .GE. 100.0) .AND. (L .LT. 120.0)) THEN
          ALPHA=2.03861
          BETA=-4.40843
          R=0.231
        SAV4=60.
        ENDIF
        IF ((L .GE. 120.0) .AND. (L .LT. 140.0)) THEN
          ALPHA=1.14271
          BETA=-1.35635
          R=0.421
        SAV4=34.
        ENDIF
        IF ((L .GE. 140.0) .AND. (L .LT. 160.0)) THEN
          ALPHA=0.79908
          BETA=1.48074
          R=0.523
        SAV4=61.
        ENDIF
        IF ((L .GE. 160.0) .AND. (L .LT. 180.0)) THEN
          ALPHA=0.94260
          BETA=8.16346
          R=0.441
        SAV4=42.
        ENDIF
        IF ((L .GE. 180.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=1.66398
          BETA=0.26775
          R=0.523
        SAV4=42.
        ENDIF
        IF ((L .GE. 200.0) .AND. (L .LT. 220.0)) THEN
          ALPHA=1.08760
          BETA=-1.02443
          R=0.523
        SAV4=45.
        ENDIF
        IF ((L .GE. 220.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=1.20087
          BETA=-2.45407
          R=0.245
        SAV4=6.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 260.0)) THEN
          ALPHA=1.13147
          BETA=-1.87916
          R=0.301
        SAV4=16.
        ENDIF
        IF ((L .GE. 260.0) .AND. (L .LT. 280.0)) THEN
          ALPHA=1.97804
          BETA=-2.92838
          R=0.338
        SAV4=21.
        ENDIF
        IF ((L .GE. 280.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=1.40086
          BETA=-1.12403
          R=0.523
        SAV4=19.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 320.0)) THEN
          ALPHA=2.06355
          BETA=-3.68278
          R=0.280
        SAV4=42.
        ENDIF
        IF ((L .GE. 320.0) .AND. (L .LT. 340.0)) THEN
          ALPHA=1.59260
          BETA=-2.18754
          R=0.364
        SAV4=15.
        ENDIF
        IF ((L .GE. 340.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=1.45589
          BETA=-1.90598
          R=0.382
        SAV4=21.
        ENDIF
      ENDIF
      IF ((B .GE. -15.) .AND. (B .LT. -5.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 10.0)) THEN
          ALPHA=2.56438
          BETA=-2.31586
          R=0.554
        SAV4=37.
        ENDIF
        IF ((L .GE. 10.0) .AND. (L .LT. 20.0)) THEN
          ALPHA=3.24095
          BETA=-2.78217
          R=0.582
        SAV4=38.
        ENDIF
        IF ((L .GE. 20.0) .AND. (L .LT. 30.0)) THEN
          ALPHA=2.95627
          BETA=-2.57422
          R=0.574
          GAMMA=0.08336
        SAV4=41.
        ENDIF
        IF ((L .GE. 30.0) .AND. (L .LT. 40.0)) THEN
          ALPHA=1.85158
          BETA=-0.67716
          R=1.152
        SAV4=4.
        ENDIF
        IF ((L .GE. 40.0) .AND. (L .LT. 50.0)) THEN
          ALPHA=1.60770
          BETA=0.35279
          R=0.661
        SAV4=24.
        ENDIF
        IF ((L .GE. 50.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=0.69920
          BETA=-0.09146
          R=0.952
          GAMMA=0.12839
        SAV4=2.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 80.0)) THEN
          ALPHA=1.36189
          BETA=-1.05290
          R=0.647
          GAMMA=0.16258
        SAV4=45.
        ENDIF
        IF ((L .GE. 80.0) .AND. (L .LT. 90.0)) THEN
          ALPHA=0.33179
          BETA=0.37629
          R=1.152
        SAV4=62.
        ENDIF
        IF ((L .GE. 90.0) .AND. (L .LT. 100.0)) THEN
          ALPHA=1.70303
          BETA=-0.75246
          R=1.132
        SAV4=31.
        ENDIF
        IF ((L .GE. 100.0) .AND. (L .LT. 110.0)) THEN
          ALPHA=1.97414
          BETA=-1.59784
          R=0.618
          GAMMA=0.12847
        SAV4=35.
        ENDIF
        IF ((L .GE. 110.0) .AND. (L .LT. 120.0)) THEN
          ALPHA=1.07407
          BETA=-0.40066
          R=1.152
          GAMMA=0.17698
        SAV4=14.
        ENDIF
        IF ((L .GE. 120.0) .AND. (L .LT. 130.0)) THEN
          ALPHA=1.69495
          BETA=-1.00071
          R=0.847
          GAMMA=0.08567
        SAV4=28.
        ENDIF
        IF ((L .GE. 130.0) .AND. (L .LT. 140.0)) THEN
          ALPHA=1.51449
          BETA=-0.08441
          R=0.897
        SAV4=12.
        ENDIF
        IF ((L .GE. 140.0) .AND. (L .LT. 150.0)) THEN
          ALPHA=1.87894
          BETA=-0.73314
          R=1.152
        SAV4=23.
        ENDIF
        IF ((L .GE. 150.0) .AND. (L .LT. 160.0)) THEN
          ALPHA=1.43670
          BETA=0.67706
          R=0.778
          GAMMA=0.42624
        SAV4=3.
        ENDIF
        IF ((L .GE. 160.0) .AND. (L .LT. 180.0)) THEN
          ALPHA=6.84802
          BETA=-5.06864
          R=0.676
        SAV4=50.
        ENDIF
        IF ((L .GE. 180.0) .AND. (L .LT. 190.0)) THEN
          ALPHA=4.16321
          BETA=-5.80016
          R=0.359
          GAMMA=0.60387
        SAV4=51.
        ENDIF
        IF ((L .GE. 190.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=0.78135
          BETA=-0.27826
          R=1.152
        SAV4=4.
        ENDIF
        IF ((L .GE. 200.0) .AND. (L .LT. 210.0)) THEN
          ALPHA=0.85535
          BETA=0.20848
          R=1.152
        SAV4=17.
        ENDIF
        IF ((L .GE. 210.0) .AND. (L .LT. 220.0)) THEN
          ALPHA=0.52521
          BETA=0.65726
          R=1.152
        SAV4=7.
        ENDIF
        IF ((L .GE. 220.0) .AND. (L .LT. 230.0)) THEN
          ALPHA=0.88376
          BETA=-0.44519
          R=0.993
          GAMMA=0.06013
        SAV4=40.
        ENDIF
        IF ((L .GE. 230.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=0.42228
          BETA=-0.26304
          R=0.803
        SAV4=26.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 250.0)) THEN
          ALPHA=0.71318
          BETA=-0.67229
          R=0.530
          GAMMA=0.20994
        SAV4=55.
        ENDIF
        IF ((L .GE. 250.0) .AND. (L .LT. 260.0)) THEN
          ALPHA=0.99606
          BETA=-0.70103
          R=0.710
          GAMMA=0.01323
        SAV4=48.
        ENDIF
        IF ((L .GE. 260.0) .AND. (L .LT. 270.0)) THEN
          ALPHA=0.91519
          BETA=-0.39690
          R=1.152
          GAMMA=0.01961
        SAV4=48.
        ENDIF
        IF ((L .GE. 270.0) .AND. (L .LT. 280.0)) THEN
          ALPHA=0.85791
          BETA=-0.29115
          R=1.152
        SAV4=19.
        ENDIF
        IF ((L .GE. 280.0) .AND. (L .LT. 290.0)) THEN
          ALPHA=1.44226
          BETA=-1.09775
          R=0.657
        SAV4=39.
        ENDIF
        IF ((L .GE. 290.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=2.55486
          BETA=-1.68293
          R=0.759
        SAV4=31.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 310.0)) THEN
          ALPHA=3.18047
          BETA=-2.69796
          R=0.589
        SAV4=40.
        ENDIF
        IF ((L .GE. 310.0) .AND. (L .LT. 320.0)) THEN
          ALPHA=2.11235
          BETA=-1.77506
          R=0.595
        SAV4=29.
        ENDIF
        IF ((L .GE. 320.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=1.75884
          BETA=-1.38574
          R=0.635
        SAV4=25.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 340.0)) THEN
          ALPHA=1.97257
          BETA=-1.55545
          R=0.634
          GAMMA=0.00043
        SAV4=34.
        ENDIF
        IF ((L .GE. 340.0) .AND. (L .LT. 350.0)) THEN
          ALPHA=1.41497
          BETA=-1.05722
          R=0.669
          GAMMA=0.03264
        SAV4=46.
        ENDIF
        IF ((L .GE. 350.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=1.17795
          BETA=-0.95012
          R=0.620
          GAMMA=0.03339
        SAV4=46.
        ENDIF          
      ENDIF
      IF ((B .GE. -5.) .AND. (B .LT. 5.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 10.0)) THEN
          ALPHA=2.62556
          BETA=-1.11097
          R=1.182
          GAMMA=0.00340
        SAV4=57.
        ENDIF
        IF ((L .GE. 10.0) .AND. (L .LT. 20.0)) THEN
          ALPHA=3.14461
          BETA=-1.01140
          R=1.555
        SAV4=42.
        ENDIF
        IF ((L .GE. 20.0) .AND. (L .LT. 30.0)) THEN
          ALPHA=4.26624
          BETA=-1.61242
          R=1.323
        SAV4=34.
        ENDIF
        IF ((L .GE. 30.0) .AND. (L .LT. 40.0)) THEN
          ALPHA=2.54447
          BETA=-0.12771
          R=1.300
        SAV4=30.
        ENDIF
        IF ((L .GE. 40.0) .AND. (L .LT. 50.0)) THEN
          ALPHA=2.27030
          BETA=-0.68720
          R=1.652
          GAMMA=0.04928
        SAV4=52.
        ENDIF
        IF ((L .GE. 50.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.34359
          BETA=-0.05416
          R=2.000
        SAV4=32.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 70.0)) THEN
          ALPHA=1.76327
          BETA=-0.26407
          R=2.000
        SAV4=37.
        ENDIF
        IF ((L .GE. 70.0) .AND. (L .LT. 80.0)) THEN
          ALPHA=2.20666
          BETA=-0.41651
          R=2.000
        SAV4=36.
        ENDIF
        IF ((L .GE. 80.0) .AND. (L .LT. 90.0)) THEN
          ALPHA=1.50130
          BETA=-0.09855
          R=1.475
        SAV4=45.
        ENDIF
        IF ((L .GE. 90.0) .AND. (L .LT. 100.0)) THEN
          ALPHA=2.43965
          BETA=-0.82128
          R=1.485
          GAMMA=0.01959
        SAV4=36.
        ENDIF
        IF ((L .GE. 100.0) .AND. (L .LT. 110.0)) THEN
          ALPHA=3.35775
          BETA=-1.16400
          R=0.841
          GAMMA=0.00298
        SAV4=35.
        ENDIF
        IF ((L .GE. 110.0) .AND. (L .LT. 120.0)) THEN
          ALPHA=2.60621
          BETA=-0.68687
          R=1.897
        SAV4=36.
        ENDIF
        IF ((L .GE. 120.0) .AND. (L .LT. 130.0)) THEN
          ALPHA=2.90112
          BETA=-0.97988
          R=1.480
        SAV4=32.
        ENDIF
        IF ((L .GE. 130.0) .AND. (L .LT. 140.0)) THEN
          ALPHA=2.55377
          BETA=-0.71214
          R=1.793
        SAV4=38.
        ENDIF
        IF ((L .GE. 140.0) .AND. (L .LT. 150.0)) THEN
          ALPHA=3.12598
          BETA=-1.21437
          R=1.287
          GAMMA=0.15298
        SAV4=23.
        ENDIF
        IF ((L .GE. 150.0) .AND. (L .LT. 160.0)) THEN
          ALPHA=3.66930
          BETA=-2.29731
          R=0.799
          GAMMA=0.33473
        SAV4=40.
        ENDIF
        IF ((L .GE. 160.0) .AND. (L .LT. 170.0)) THEN
          ALPHA=2.15465
          BETA=-0.70690
          R=1.524
          GAMMA=0.14017
        SAV4=37.
        ENDIF
        IF ((L .GE. 170.0) .AND. (L .LT. 180.0)) THEN
          ALPHA= 1.82465
          BETA=-0.60223
          R=1.515
          GAMMA= 0.20730
        SAV4=29.
        ENDIF
        IF ((L .GE. 180.0) .AND. (L .LT. 190.0)) THEN
          ALPHA=1.76269
          BETA=-0.35945
          R=2.000
          GAMMA=0.08052
        SAV4=28.
        ENDIF
        IF ((L .GE. 190.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=1.06085
          BETA=-0.14211
          R=2.000
        SAV4=48.
        ENDIF
        IF ((L .GE. 200.0) .AND. (L .LT. 210.0)) THEN
          ALPHA=1.21333
          BETA=-0.23225
          R=2.000
        SAV4=57.
        ENDIF
        IF ((L .GE. 210.0) .AND. (L .LT. 220.0)) THEN
          ALPHA=0.58326
          BETA=-0.06097
          R=2.000
          GAMMA=0.36962
        SAV4=30.
        ENDIF
        IF ((L .GE. 220.0) .AND. (L .LT. 230.0)) THEN
          ALPHA=0.74200
          BETA=-0.19293
          R=1.923
          GAMMA=0.07459
        SAV4=48.
        ENDIF
        IF ((L .GE. 230.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=0.67520
          BETA=-0.21041
          R=1.604
          GAMMA=0.16602
        SAV4=49.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 250.0)) THEN
          ALPHA=0.62609
          BETA=-0.25312
          R=1.237
          GAMMA=0.14437
        SAV4=73.
        ENDIF
        IF ((L .GE. 250.0) .AND. (L .LT. 260.0)) THEN
          ALPHA=0.61415
          BETA=-0.13788
          R=2.000
          GAMMA=0.26859
        SAV4=42.
        ENDIF
        IF ((L .GE. 260.0) .AND. (L .LT. 270.0)) THEN
          ALPHA=0.58108
          BETA=0.01195
          R=2.000
          GAMMA=0.07661
        SAV4=40.
        ENDIF
        IF ((L .GE. 270.0) .AND. (L .LT. 280.0)) THEN
          ALPHA=0.68352
          BETA=-0.10743
          R=2.000
          GAMMA=0.00849
        SAV4=50.
        ENDIF
        IF ((L .GE. 280.0) .AND. (L .LT. 290.0)) THEN
          ALPHA=0.61747
          BETA=0.02675
          R=2.000
        SAV4=49.
        ENDIF
        IF ((L .GE. 290.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=1.06827
          BETA=-0.26290
          R=2.000
        SAV4=44.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 310.0)) THEN
          ALPHA=1.53631
          BETA=-0.36833
          R=2.000
          GAMMA=0.02960
        SAV4=37.
        ENDIF
        IF ((L .GE. 310.0) .AND. (L .LT. 320.0)) THEN
          ALPHA=1.94300
          BETA=-0.71445
          R=1.360
          GAMMA=0.15643
        SAV4=36.
        ENDIF
        IF ((L .GE. 320.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=1.22185
          BETA=-0.40185
          R=1.520
          GAMMA=0.07354
        SAV4=48.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 340.0)) THEN
          ALPHA=1.79429
          BETA=-0.48657
          R=1.844
        SAV4=43.
        ENDIF
        IF ((L .GE. 340.0) .AND. (L .LT. 350.0)) THEN
          ALPHA=2.29545
          BETA=-0.84096
          R=1.365
        SAV4=32.
        ENDIF
        IF ((L .GE. 350.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=2.07408
          BETA=-0.64745
          R=1.602
          GAMMA=0.12750
        SAV4=36.
        ENDIF
      ENDIF
      IF ((B .GE. 5.) .AND. (B .LT. 15.)) THEN
        IF ((L .GE . 0.0) .AND. (L .LT. 10.0)) THEN
          ALPHA=2.94213
          BETA=-2.09258
          R=0.703
          GAMMA=0.05490
        SAV4=41.
        ENDIF
        IF ((L .GE. 10.0) .AND. (L .LT. 30.0)) THEN
          ALPHA=3.04627
          BETA=7.71159
          R=0.355
        SAV4=45.
        ENDIF
        IF ((L .GE. 30.0) .AND. (L .LT. 40.0)) THEN
          ALPHA=3.78033
          BETA=-3.91956
          R=0.482
        SAV4=42.
        ENDIF
        IF  ((L .GE. 40.0) .AND. (L .LT. 50.0)) THEN
          ALPHA=2.18119
          BETA=-2.4050
          R=0.453
        SAV4=27.
        ENDIF
        IF  ((L .GE. 50.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.45372
          BETA=-0.49522
          R=1.152
        SAV4=31.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 70.0)) THEN
          ALPHA=1.05051
          BETA=-1.01704
          R=0.516
        SAV4=2.
        ENDIF
        IF ((L .GE. 70.0) .AND. (L .LT. 80.0)) THEN
          ALPHA=0.48416
          BETA=-0.27182
          R=0.891
          GAMMA=0.08639
        SAV4=94.
        ENDIF
        IF ((L .GE. 80.0) .AND. (L .LT. 90.0)) THEN
          ALPHA=0.61963
          BETA=0.41697
          R=1.152
          GAMMA=0.47171
        SAV4=35.
        ENDIF
        IF ((L .GE. 90.0) .AND. (L .LT. 100.0)) THEN
          ALPHA=4.40348
          BETA=-2.95611
          R=0.745
        SAV4=52.
        ENDIF
        IF ((L .GE. 100.0) .AND. (L .LT. 120.0)) THEN
          ALPHA=2.50938
          BETA=-0.56541
          R=1.152
        SAV4=27.
        ENDIF
        IF ((L .GE. 120.0) .AND. (L .LT. 130.0)) THEN
          ALPHA=0.44180
          BETA=1.58923
          R=0.949
        SAV4=4.
        ENDIF
        IF ((L .GE. 130.0) .AND. (L .LT. 140.0)) THEN
          ALPHA=3.96084
          BETA=-3.37374
          R=0.587
          GAMMA=0.34109
        SAV4=40.
        ENDIF
        IF ((L .GE. 140.0) .AND. (L .LT. 160.0)) THEN
          ALPHA=2.53335
          BETA=-0.40541
          R=1.152
        SAV4=38.
        ENDIF
        IF ((L .GE. 160.0) .AND. (L .LT. 170.0)) THEN
          ALPHA=2.03760
          BETA=-0.66317
          R=1.152
        SAV4=23.
        ENDIF
        IF ((L .GE. 170.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=1.06946
          BETA=-0.87395
          R=0.612
          GAMMA=0.29230
        SAV4=29.
        ENDIF
        IF ((L .GE. 200.0) .AND. (L .LT. 210.0)) THEN
          ALPHA=0.86348
          BETA=-0.65870
          R=0.655
          GAMMA=0.09089
        SAV4=79.
        ENDIF
        IF ((L .GE. 210.0) .AND. (L .LT. 230.0)) THEN
          ALPHA=0.30117
          BETA=-0.16136
          R=0.933
          GAMMA=0.07495
        SAV4=17.
        ENDIF
        IF ((L .GE. 230.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=0.75171
          BETA=-0.57143
          R=0.658
          GAMMA=0.00534
        SAV4=12.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 250.0)) THEN
          ALPHA=1.97427
          BETA=-2.02654
          R=0.487
        SAV4=67.
        ENDIF
        IF ((L .GE. 250.0) .AND. (L .LT. 260.0)) THEN
          ALPHA=1.25208
          BETA=-1.47763
          R=0.424
          GAMMA=0.31600
        SAV4=19.
        ENDIF
        IF ((L .GE. 260.0) .AND. (L .LT. 270.0)) THEN
          ALPHA=0.89448
          BETA=-0.43870
          R=1.019
        SAV4=5.
        ENDIF
        IF ((L .GE. 270.0) .AND. (L .LT. 280.0)) THEN
          ALPHA=0.81141
          BETA=-0.51001
          R=0.795
          GAMMA=0.03505
        SAV4=27.
        ENDIF
        IF ((L .GE. 280.0) .AND. (L .LT. 290.0)) THEN
          ALPHA=0.83781
          BETA=-0.44138
          R=0.949
          GAMMA=0.02820
        SAV4=50.
        ENDIF
        IF ((L .GE. 290.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=1.10600
          BETA=-0.86263
          R=0.641
          GAMMA=0.03402
        SAV4=28.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 310.0)) THEN
          ALPHA=1.37040
          BETA=-1.02779
          R=0.667
          GAMMA=0.05608
        SAV4=28.
        ENDIF
        IF ((L .GE. 310.0) .AND. (L .LT. 320.0)) THEN
          ALPHA=1.77590
          BETA=-1.26951
          R=0.699
          GAMMA=0.06972
        SAV4=37.
        ENDIF
        IF ((L .GE. 320.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=1.20865
          BETA=-0.70679
          R=0.855
          GAMMA=0.02902
        SAV4=35.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 340.0)) THEN
          ALPHA=2.28830
          BETA=-1.71890
          R=0.666
          GAMMA=0.22887
        SAV4=42.
        ENDIF
        IF ((L .GE. 340.0) .AND. (L .LT. 350.0)) THEN
          ALPHA=3.26278
          BETA=-0.94181
          R=1.152
        SAV4=38.
        ENDIF
        IF ((L .GE. 350.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=2.58100
          BETA=-1.69237
          R=0.763
        SAV4=53.
        ENDIF
      ENDIF
      IF ((B .GE. 15.) .AND. (B .LT. 30.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 20.0)) THEN
          ALPHA=6.23279
          BETA=-10.30384
          R=0.302
        SAV4=42.
        ENDIF
        IF ((L .GE. 20.0) .AND. (L .LT. 40.0)) THEN
          ALPHA=4.47693
          BETA=-7.28366
          R=0.307
        SAV4=29.
        ENDIF
        IF ((L .GE. 40.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.22938
          BETA=-1.19030
          R=0.516
        SAV4=5.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 80.0)) THEN
          ALPHA=0.84291
          BETA=-1.59338
          R=0.265
        SAV4=4.
        ENDIF
        IF ((L .GE. 80.0) .AND. (L .LT. 100.0)) THEN
          ALPHA=0.23996
          BETA=0.06304
          R=0.523
        SAV4=32.
        ENDIF
        IF ((L .GE. 100.0) .AND. (L .LT. 140.0)) THEN
          ALPHA=0.40062
          BETA=-1.75628
          R=0.114
        SAV4=16.
        ENDIF
        IF ((L .GE. 140.0) .AND. (L .LT. 180.0)) THEN
          ALPHA=0.56898
          BETA=-0.53331
          R=0.523
        SAV4=41.
        ENDIF
        IF ((L .GE. 180.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=-0.95721
          BETA=11.69217
          R=0.240
        SAV4=2.
        ENDIF
        IF ((L .GE. 200.0) .AND. (L .LT. 220.0)) THEN
          ALPHA=-0.19051
          BETA=1.45670
          R=0.376
        SAV4=1.
        ENDIF
        IF ((L .GE. 220.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=2.31305
          BETA=-7.82531
          R=0.148
        SAV4=95.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 260.0)) THEN
          ALPHA=1.39169
          BETA=-1.72984
          R=0.402
        SAV4=6.
        ENDIF
        IF ((L .GE. 260.0) .AND. (L .LT. 280.0)) THEN
          ALPHA=1.59418
          BETA=-1.28296
          R=0.523
        SAV4=36.
        ENDIF
        IF ((L .GE. 280.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=1.57082
          BETA=-1.97295
          R=0.398
        SAV4=10.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 320.0)) THEN
          ALPHA=1.95998
          BETA=-3.26159
          R=0.300
        SAV4=11.
        ENDIF
        IF ((L .GE. 320.0) .AND. (L .LT. 340.0)) THEN
          ALPHA=2.59567
          BETA=-4.84133
          R=0.268
        SAV4=37.
        ENDIF
        IF ((L .GE. 340.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=5.30273
          BETA=-7.43033
          R=0.357
        SAV4=37.
        ENDIF
      ENDIF
      IF ((B .GE. 30.) .AND. (B .LT. 45.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 20.0)) THEN
          ALPHA=2.93960
          BETA=-6.48049
          R=0.227
        SAV4=77.
        ENDIF
        IF ((L .GE. 20.0) .AND. (L .LT. 50.0)) THEN
          ALPHA=1.65864
          BETA=-9.99317
          R=0.083
        SAV4=100.
        ENDIF
        IF ((L .GE. 50.0) .AND. (L .LT. 80.0)) THEN
          ALPHA=1.71831
          BETA=-7.25286
          R=0.118
        SAV4=28.
        ENDIF
        IF ((L .GE. 80.0) .AND. (L .LT. 110.0)) THEN
          ALPHA=1.33617
          BETA=-10.39799
          R=0.064
        SAV4=100.
        ENDIF
        IF ((L .GE. 110.0) .AND. (L .LT. 160.0)) THEN
          ALPHA=-0.31330
          BETA=1.35622
          R=0.329
        SAV4=24.
        ENDIF
        IF ((L .GE. 160.0) .AND. (L .LT. 190.0)) THEN
          ALPHA=1.51984
          BETA=-8.69502
          R=0.087
        SAV4=100.
        ENDIF
        IF ((L .GE. 190.0) .AND. (L .LT. 220.0)) THEN
          ALPHA=-0.50758
          BETA=4.73320
          R=0.250
        SAV4=78.
        ENDIF
        IF ((L .GE. 220.0) .AND. (L .LT. 250.0)) THEN
          ALPHA=1.25864
          BETA=-12.59627
          R=0.050
        SAV4=70.
        ENDIF
        IF ((L .GE. 250.0) .AND. (L .LT. 280.0)) THEN
          ALPHA=1.54243
          BETA=-3.76065
          R=0.205
        SAV4=10.
        ENDIF
        IF ((L .GE. 280.0) .AND. (L .LT. 320.0)) THEN
          ALPHA=2.72258
          BETA=-7.47806
          R=0.182
        SAV4=5.
        ENDIF
        IF ((L .GE. 320.0) .AND. (L .LT. 340.0)) THEN
          ALPHA=2.81545
          BETA=-5.52139
          R=0.255
        SAV4=10.
        ENDIF
        IF ((L .GE. 340.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=2.23818
          BETA=0.81772
          R=0.329
        SAV4=19.
        ENDIF
      ENDIF
      IF ((B .GE. 45.) .AND. (B .LT. 60.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.38587
          BETA=-9.06536
          R=0.076
        SAV4=3.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 90.0)) THEN
          ALPHA=2.28570
          BETA=-9.88812
          R=0.116
        SAV4=3.
        ENDIF
        IF ((L .GE. 90.0) .AND. (L .LT. 110.0)) THEN
          ALPHA=1.36385
          BETA=-8.10127
          R=0.084
        SAV4=4.
        ENDIF
        IF ((L .GE. 110.0) .AND. (L .LT. 170.0)) THEN
          ALPHA=0.05943
          BETA=-1.08126
          R=0.027
        SAV4=50.
        ENDIF
        IF ((L .GE. 170.0) .AND. (L .LT. 200.0)) THEN
          ALPHA=1.40171
          BETA=-3.21783
          R=0.218
        SAV4=99.
        ENDIF
        IF ((L .GE. 200.0) .AND. (L .LT. 230.0)) THEN
          ALPHA=0.14718
          BETA=3.92670
          R=0.252
        SAV4=14.
        ENDIF
        IF ((L .GE. 230.0) .AND. (L .LT. 290.0)) THEN
          ALPHA=0.57124
          BETA=-4.30242
          R=0.066
        SAV4=10.
        ENDIF
        IF ((L .GE. 290.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=3.69891
          BETA=-19.62204
          R=0.094
        SAV4=5.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=1.19568
          BETA=-0.45043
          R=0.252
        SAV4=9.
        ENDIF
      ENDIF
      IF ((B .GE. 60.) .AND. (B .LE. 90.)) THEN
        IF ((L .GE. 0.0) .AND. (L .LT. 30.0)) THEN
          ALPHA=0.69443
          BETA=-0.27600
          R=0.153
        SAV4=100.
        ENDIF
        IF ((L .GE. 30.0) .AND. (L .LT. 60.0)) THEN
          ALPHA=1.11811
          BETA=0.71179
          R=0.085
        SAV4=73.
        ENDIF
        IF ((L .GE. 60.0) .AND. (L .LT. 90.0)) THEN
          ALPHA=1.10427
          BETA=-2.37654
          R=0.123
        SAV4=100.
        ENDIF
        IF ((L .GE. 90.0) .AND. (L .LT. 120.0)) THEN
          ALPHA=-0.42211
          BETA=5.24037
          R=0.184
        SAV4=12.
        ENDIF
        IF ((L .GE. 120.0) .AND. (L .LT. 150.0)) THEN
          ALPHA=0.87576
          BETA=-4.38033
          R=0.100
        SAV4=35.
        ENDIF
        IF ((L .GE. 150.0) .AND. (L .LT. 180.0)) THEN
          ALPHA=1.27477
          BETA=-4.98307
          R=0.128
        SAV4=72.
        ENDIF
        IF ((L .GE. 180.0) .AND. (L .LT. 210.0)) THEN
          ALPHA=1.19512
          BETA=-6.58464
          R=0.091
        SAV4=49.
        ENDIF
        IF ((L .GE. 210.0) .AND. (L .LT. 240.0)) THEN
          ALPHA=0.97581
          BETA=-4.89869
          R=0.100
        SAV4=95.
        ENDIF
        IF ((L .GE. 240.0) .AND. (L .LT. 270.0)) THEN
          ALPHA=0.54379
          BETA=-0.84403
          R=0.207
        SAV4=35.
        ENDIF
        IF ((L .GE. 270.0) .AND. (L .LT. 300.0)) THEN
          ALPHA=-0.85054
          BETA=13.01249
          R=0.126
        SAV4=39.
        ENDIF
        IF ((L .GE. 300.0) .AND. (L .LT. 330.0)) THEN
          ALPHA=0.74347
          BETA=-1.39825
          R=0.207
        SAV4=10.
        ENDIF
        IF ((L .GE. 330.0) .AND. (L .LT. 360.0)) THEN
          ALPHA=0.77310
          BETA=-4.45005
          R=0.087
        SAV4=16.
        ENDIF
      ENDIF
C----------------------------------------------------------------------
C     We are assuming the the extinction takes on the curve out to a 
C       value at R. Then it becomes a linear fit until it reaches the 
C       limiting distance of 2kpc.     
C     After 2kpc we are assuming the linear fit continues.   
C     NOTE: During the testing of this subroutine it was found that at 
C       distances of D less than R, the extinction takes on a negative
C       value. It was determined that the linear curve fit the authors
C       used was trying to fit points much farther out. Therefore, we
C       have replaced any negative AV value with 0.
C----------------------------------------------------------------------
      IF ((D .LE. R)) THEN
        AV=ALPHA*D+BETA*D*D
        IF ((AV .LT. 0.)) AV=0.
        SAV=SQRT(0.0225+(SAV4/100.*AV)*(SAV4/100.*AV))
      ENDIF
      IF ((D .GT. R) .AND. (D .LE. 2.)) THEN 
      AV=ALPHA*R+BETA*R*R+GAMMA*(D-R)
      SAV=SQRT(0.0225+(SAV4/100.*AV)*(SAV4/100.*AV))
      ENDIF
      IF ((D .GT. 2.)) THEN
      AV= ALPHA*R+BETA*R*R+GAMMA*(2.-R)
      SAV=SQRT(0.0225+(SAV4/100.*AV)*(SAV4/100.*AV))
      IF (D*ABS(SIN(B*3.141592/180.)) .LT. 0.1) THEN
          AV=AV+A0*(D-2.)
        SAV=SQRT(SAV*SAV+0.25*(D-2.)*(D-2.))
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
C     Assuming that the paper by Arenou et.al. has a limit of 2kpc.
C     Beyond 2kpc, used Ao + Parenago. Where ao = 1.5 mag/kpc and  
C           sigma ao = 50%.                                            
C-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE HIGH_LAT(L,B,D,AV,SAV,IHFLG)
      IMPLICIT REAL*4 (A-H,L-Z)
C-----------------------------------------------------------------------
C     From Penprase, ApJ. Supp. 83, 273 (1992),
C          Magnani et al., ApJ. 295, 402 (1985),
C          Keto & Myers, ApJ. 304, 466 (1986),
C          Odenwald, ApJ. 325, 320 (1988),
C          Desert et al., ApJ. 334, 815 (1988),
C          Hughes et al., AJ 105, 571 (1993),
C          Kenyon et al., AJ 108, 1872 (1994),
C          Cernis, ApSS 166, 315 (1990),
C          Cernis, Baltic Astron. 2, 214 (1993),
C          Rossano, AJ, 83, 234 (1978),
C          Rossano, AJ, 83, 241 (1978),
C          Kutner et al., ApJ, 215, 521 (1977).
C-----------------------------------------------------------------------
      AV=-99.
      SAV=-99.
      IHFLG=0
      IF (ABS(B) .LT. 7.) RETURN
C
C     Penprase values
C
      SAV=0.24
      IF ((L .GE. 4.6) .AND. (L .LT. 5.3) .AND. (B .GE. 30.3)
     &  .AND. (B .LT. 31.0)) THEN                 ! DBB03
         IHFLG=1
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. 0.15) .AND. (D .LE. 0.21)) AV=23.33*(D-0.15)
         IF (D .GT. 0.21) AV=1.3998
         RETURN
      ENDIF
      IF ((L .GE. 3.3) .AND. (L .LT. 6.3) .AND. (B .GE. 35.0)
     &  .AND. (B .LT. 37.)) THEN                 ! MBM57
         IHFLG=1
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. 0.15) .AND. (D .LE. 0.21)) AV=23.33*(D-0.15)
         IF (D .GT. 0.21) AV=1.3998
         RETURN
      ENDIF
      IF ((L .GE. 10.0) .AND. (L .LT. 13.5) .AND. (B .GE. 39.5)
     &  .AND. (B .LT. 41.0)) THEN                 ! DBB11,DBB13
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.18)) AV=6.25*(D-0.1)
         IF (D .GT. 0.18) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 12.0) .AND. (L .LT. 14.0) .AND. (B .GE. 38.5)
     &  .AND. (B .LT. 39.5)) THEN                 ! DBB11,DBB13
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.18)) AV=6.25*(D-0.1)
         IF (D .GT. 0.18) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 9.0) .AND. (L .LT. 12.0) .AND. (B .GE. 45.0)
     &  .AND. (B .LT. 46.5)) THEN                 ! DBB20
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.26)) AV=2.5*(D-0.1)
         IF (D .GT. 0.26) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 12.0) .AND. (L .LT. 17.0) .AND. (B .GE. 44.5)
     &  .AND. (B .LT. 46.0)) THEN                 ! DBB20
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.26)) AV=2.5*(D-0.1)
         IF (D .GT. 0.26) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 17.7) .AND. (L .LT. 19.8) .AND. (B .GE. 16.3)
     &  .AND. (B .LT. 19.9)) THEN                 ! DBB25
         IHFLG=1
         IF (D .LE. 0.06) AV=0.
         IF ((D .GT. 0.06) .AND. (D .LE. 0.16)) AV=9.*(D-0.06)
         IF (D .GT. 0.16) AV=0.9
         RETURN
      ENDIF
      IF ((L .GE. 23.5) .AND. (L .LT. 26.3) .AND. (B .GE. -22.5)
     &  .AND. (B .LT. -17.0)) THEN                ! DBB36
         IHFLG=1
         IF (D .LE. 0.12) AV=0.
         IF ((D .GT. 0.12) .AND. (D .LE. 0.22)) AV=10.*(D-0.12)
         IF (D .GT. 0.22) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 24.8) .AND. (L .LT. 30.0) .AND. (B .GE. 28.3)
     &  .AND. (B .LT. 31.2)) THEN                 ! DBB45
         IHFLG=1
         IF (D .LE. 0.07) AV=0.
         IF ((D .GT. 0.07) .AND. (D .LE. 0.17)) AV=20.*(D-0.07)
         IF (D .GT. 0.17) AV=2.0
         RETURN
      ENDIF
      IF ((L .GE. 37.0) .AND. (L .LT. 42.0) .AND. (B .GE. -35.5)
     &  .AND. (B .LT. -30.0)) THEN                ! MBM46,47,DBB55
         IHFLG=1
         IF (D .LE. 0.24) AV=0.
         IF ((D .GT. 0.24) .AND. (D .LE. 0.35)) AV=3.64*(D-0.24)
         IF (D .GT. 0.35) AV=0.4004
         RETURN
      ENDIF
      IF ((L .GE. 40.0) .AND. (L .LT. 46.0) .AND. (B .GE. -37.0)
     &  .AND. (B .LT. -32.0)) THEN                ! MBM46,47,DBB55
         IHFLG=1
         IF (D .LE. 0.24) AV=0.
         IF ((D .GT. 0.24) .AND. (D .LE. 0.35)) AV=3.64*(D-0.24)
         IF (D .GT. 0.35) AV=0.4004
         RETURN
      ENDIF
      IF ((L .GE. 89.5) .AND. (L .LT. 90.5) .AND. (B .GE. 38.0)
     &  .AND. (B .LT. 39.0)) THEN                 ! MBM41
         IHFLG=1
         IF (D .LE. 0.2) AV=0.
         IF ((D .GT. 0.2) .AND. (D .LE. 0.3)) AV=9.*(D-0.2)
         IF (D .GT. 0.3) AV=0.9
         RETURN
      ENDIF
      IF ((L .GE. 93.2) .AND. (L .LT. 94.7) .AND. (B .GE. -35.0)
     &  .AND. (B .LT. -30.5)) THEN                ! MBM53
         IHFLG=1
         IF (D .LE. 0.07) AV=0.
         IF ((D .GT. 0.07) .AND. (D .LE. 0.13)) AV=5.*(D-0.07)
         IF (D .GT. 0.13) AV=0.3
         RETURN
      ENDIF
      IF ((L .GE. 114.2) .AND. (L .LT. 116.0) .AND. (B .GE. -47.7)
     &  .AND. (B .LT. -45.0)) THEN                ! DBB153
         IHFLG=1
         IF (D .LE. 0.12) AV=0.
         IF ((D .GT. 0.12) .AND. (D .LE. 0.14)) AV=25.*(D-0.12)
         IF (D .GT. 0.14) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 114.2) .AND. (L .LT. 117.5) .AND. (B .GE. -45.0)
     &  .AND. (B .LT. -43.0)) THEN                ! DBB153
         IHFLG=1
         IF (D .LE. 0.12) AV=0.
         IF ((D .GT. 0.12) .AND. (D .LE. 0.14)) AV=25.*(D-0.12)
         IF (D .GT. 0.14) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 125.7) .AND. (L .LT. 132.0) .AND. (B .GE. -48.0)
     &  .AND. (B .LT. -45.8)) THEN                 ! MBM3
         IHFLG=1
         IF (D .LE. 0.09) AV=0.
         IF ((D .GT. 0.09) .AND. (D .LE. 0.18)) AV=3.33*(D-0.09)
         IF (D .GT. 0.18) AV=0.2997
         RETURN
      ENDIF
      IF ((L .GE. 132.0) .AND. (L .LT. 141.5) .AND. (B .GE. -45.8)
     &  .AND. (B .LT. -43.0)) THEN                 ! MBM4
         IHFLG=1
         IF (D .LE. 0.09) AV=0.
         IF ((D .GT. 0.09) .AND. (D .LE. 0.18)) AV=3.33*(D-0.09)
         IF (D .GT. 0.18) AV=0.2997
         RETURN
      ENDIF
      IF ((L .GE. 141.0) .AND. (L .LT. 143.0) .AND. (B .GE. 34.5)
     &  .AND. (B .LT. 38.5)) THEN                  ! MBM27-31
         IHFLG=1
         IF (D .LE. 0.13) AV=0.
         IF ((D .GT. 0.13) .AND. (D .LE. 0.22)) AV=12.22*(D-0.13)
         IF (D .GT. 0.22) AV=1.0998
         RETURN
      ENDIF
      IF ((L .GE. 143.0) .AND. (L .LT. 151.0) .AND. (B .GE. 38.5)
     &  .AND. (B .LT. 40.5)) THEN                  ! MBM27-31
         IHFLG=1
         IF (D .LE. 0.13) AV=0.
         IF ((D .GT. 0.13) .AND. (D .LE. 0.22)) AV=12.22*(D-0.13)
         IF (D .GT. 0.22) AV=1.0998
         RETURN
      ENDIF
      IF ((L .GE. 170.5) .AND. (L .LT. 174.5) .AND. (B .GE. -42.8)
     &  .AND. (B .LT. -40.0)) THEN                 ! DBB259
         IHFLG=1
         IF (D .LE. 0.08) AV=0.
         IF ((D .GT. 0.08) .AND. (D .LE. 0.2)) AV=4.17*(D-0.08)
         IF (D .GT. 0.2) AV=0.5004
         RETURN
      ENDIF
      IF ((L .GE. 175.0) .AND. (L .LT. 179.0) .AND. (B .GE. -40.0)
     &  .AND. (B .LT. -39.0)) THEN                 ! DBB267
         IHFLG=1
         IF (D .LE. 0.11) AV=0.
         IF ((D .GT. 0.11) .AND. (D .LE. 0.19)) AV=5.*(D-0.11)
         IF (D .GT. 0.19) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 186.0) .AND. (L .LT. 190.0) .AND. (B .GE. -39.5) .AND.
     &    (B .LT. -30.7)) THEN                      ! MBM18, DBB276
         IHFLG=1
         IF (D .LE. 0.14) AV=0.
         IF ((D .GT. 0.14) .AND. (D .LE. 0.2)) AV=25.*(D-0.14)
         IF (D .GT. 0.2) AV=1.5
         RETURN
      ENDIF
      IF ((L .GE. 187.0) .AND. (L .LT. 194.0) .AND. (B .GE. -53.5) .AND.
     &    (B .LT. -52.0)) THEN                      ! MBM15
         IHFLG=1
         IF (D .LE. 0.11) AV=0.
         IF ((D .GT. 0.11) .AND. (D .LE. 0.16)) AV=8.*(D-0.11)
         IF (D .GT. 0.16) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 208.0) .AND. (L .LT. 210.0) .AND. (B .GE. -30.0) .AND.
     &    (B .LT. -28.0)) THEN                      ! MBM21
         IHFLG=1
         IF (D .LE. 0.4) AV=0.
         IF ((D .GT. 0.4) .AND. (D .LE. 0.5)) AV=(D-0.4)
         IF (D .GT. 0.5) AV=1.5
         RETURN
      ENDIF
      IF ((L .GE. 206.2) .AND. (L .LT. 208.0) .AND. (B .GE. -28.0) .AND.
     &    (B .LT. -25.5)) THEN                      ! MBM22
         IHFLG=1
         IF (D .LE. 0.4) AV=0.
         IF ((D .GT. 0.4) .AND. (D .LE. 0.5)) AV=(D-0.4)
         IF (D .GT. 0.5) AV=1.5
         RETURN
      ENDIF
      IF ((L .GE. 210.7) .AND. (L .LT. 212.0) .AND. (B .GE. -37.0) .AND.
     &    (B .LT. -35.6)) THEN                      ! MBM20
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.13)) AV=30.*(D-0.1)
         IF (D .GT. 0.13) AV=0.9
         RETURN
      ENDIF
      IF ((L .GE. 212.33) .AND. (L .LT. 216.) .AND. (B .GE. 25.5) .AND.
     &    (B .LT. 26.5)) THEN                       ! G213+26
         IHFLG=1
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. 0.15) .AND. (D .LE. 0.22)) AV=7.14*(D-0.15)
         IF (D .GT. 0.22) AV=0.4998
         RETURN
      ENDIF
      IF ((L .GE. 259.4) .AND. (L .LT. 259.7) .AND. (B .GE. -16.7) .AND.
     &    (B .LT. -16.3)) THEN                      ! K259.5-16.5
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.45)) AV=3.14*(D-0.1)
         IF (D .GT. 0.45) AV=1.099
         RETURN
      ENDIF
      IF ((L .GE. 272.8) .AND. (L .LT. 273.1) .AND. (B .GE. 29.17) .AND.
     &    (B .LT. 29.45)) THEN                      ! K272.9+29.3
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.23)) AV=3.85*(D-0.1)
         IF (D .GT. 0.23) AV=0.5005
         RETURN
      ENDIF
      IF ((L .GE. 281.0) .AND. (L .LT. 282.5) .AND. (B .GE. 15.7) .AND.
     &    (B .LT. 16.6)) THEN                       ! DBB414
         IHFLG=1
         IF (D .LE. 0.2) AV=0.
         IF ((D .GT. 0.2) .AND. (D .LE. 0.3)) AV=1.*(D-0.2)
         IF (D .GT. 0.3) AV=0.1
         RETURN
      ENDIF
      IF ((L .GE. 292.5) .AND. (L .LT. 293.0) .AND. (B .GE. -20.0) .AND.
     &    (B .LT. -19.0)) THEN                      ! K292.7-19.8
         IHFLG=1
         IF (D .LE. 0.1) AV=0.
         IF ((D .GT. 0.1) .AND. (D .LE. 0.2)) AV=5.*(D-0.1)
         IF (D .GT. 0.2) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 293.8) .AND. (L .LT. 294.5) .AND. (B .GE. 14.6) .AND.
     &    (B .LT. -14.)) THEN                       ! K294.4-14.3
         IHFLG=1
         IF (D .LE. 0.11) AV=0.
         IF ((D .GT. 0.11) .AND. (D .LE. 0.21)) AV=12.*(D-0.11)
         IF (D .GT. 0.21) AV=1.2
         RETURN
      ENDIF
      SAV=0.3
C
C     Selected values from Magnani et al.
C
      IF ((L .GE. 116.75) .AND. (L .LT. 117.75) .AND. (B .GE. -52.9)
     &  .AND. (B .LT. -51.9)) THEN                 ! MBM2
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=300.*(D-0.099)
         IF (D .GT. 0.101) AV=0.6
         RETURN
      ENDIF
      IF ((L .GE. 69.0) .AND. (L .LT. 71.0) .AND. (B .GE. -22.75)
     &  .AND. (B .LT. -21.75)) THEN                ! MBM7,8
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=550.*(D-0.099)
         IF (D .GT. 0.101) AV=1.1
         RETURN
      ENDIF
      IF ((L .GE. 169.5) .AND. (L .LT. 172.5) .AND. (B .GE. -38.5)
     &  .AND. (B .LT. -38.5)) THEN                 ! MBM16
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=650.*(D-0.099)
         IF (D .GT. 0.101) AV=1.3
         RETURN
      ENDIF
      IF ((L .GE. 185.75) .AND. (L .LT. 186.0) .AND. (B .GE. -30.0)
     &  .AND. (B .LT. -29.5)) THEN                 ! MBM19
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=400.*(D-0.099)
         IF (D .GT. 0.101) AV=0.8
         RETURN
      ENDIF
      IF ((L .GE. 171.75) .AND. (L .LT. 172.) .AND. (B .GE. 26.5) .AND.
     &    (B .LT. 27.0)) THEN                      ! MBM23,24
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=350.*(D-0.099)
         IF (D .GT. 0.101) AV=0.7
         RETURN
      ENDIF
      IF ((L .GE. 173.25) .AND. (L .LT. 173.75) .AND.(B .GE. 31.) .AND.
     &    (B .LT. 32.0)) THEN                      ! MBM25
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=350.*(D-0.099)
         IF (D .GT. 0.101) AV=0.7
         RETURN
      ENDIF
      IF ((L .GE. 155.0) .AND. (L .LT. 157.0) .AND. (B .GE. 32.0) .AND.
     &    (B .LT. 34.0)) THEN                      ! MBM26
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=250.*(D-0.099)
         IF (D .GT. 0.101) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 146.75) .AND. (L .LT. 147.75) .AND.(B .GE. 40.5).AND.
     &    (B .LT. 41.0)) THEN                      ! MBM32
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=450.*(D-0.099)
         IF (D .GT. 0.101) AV=0.9
         RETURN
      ENDIF
      IF ((L .GE. 37.0) .AND. (L .LT. 38.0) .AND. (B .GE. 44.0) .AND.
     &    (B .LT. 45.0)) THEN                      ! MBM40
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=350.*(D-0.099)
         IF (D .GT. 0.101) AV=0.7
         RETURN
      ENDIF
      IF ((L .GE. 90.0) .AND. (L .LT. 92.0) .AND. (B .GE. 36.5) .AND.
     &    (B .LT. 38.0)) THEN                      ! MBM42-44
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=350.*(D-0.099)
         IF (D .GT. 0.101) AV=0.7
         RETURN
      ENDIF
      IF ((L .GE. 64.5) .AND. (L .LT. 65.0) .AND. (B .GE. -27.0).AND.
     &    (B .LT. -26.5)) THEN                     ! MBM49
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=250.*(D-0.099)
         IF (D .GT. 0.101) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 69.9) .AND. (L .LT. 70.1) .AND. (B .GE. -31.25) .AND.
     &    (B .LT. -31.0)) THEN                     ! MBM50
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=400.*(D-0.099)
         IF (D .GT. 0.101) AV=0.8
         RETURN
      ENDIF
      IF ((L .GE. 73.5) .AND. (L .LT. 74.0) .AND. (B .GE. -52.0) .AND.
     &    (B .LT. -51.0)) THEN                     ! MBM51,52
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=200.*(D-0.099)
         IF (D .GT. 0.101) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 84.0) .AND. (L .LT. 90.0) .AND. (B .GE. -41.75).AND.
     &    (B .LT. -40.0)) THEN                     ! MBM54,55
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=550.*(D-0.099)
         IF (D .GT. 0.101) AV=1.1
         RETURN
      ENDIF
      IF ((L .GE. 90.0) .AND. (L .LT. 95.0) .AND. (B .GE. -40.) .AND.
     &    (B .LT. -37.5)) THEN                     ! MBM54,55
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=550.*(D-0.099)
         IF (D .GT. 0.101) AV=1.1
         RETURN
      ENDIF
      IF ((L .GE.102.75) .AND. (L .LT.103.75) .AND.(B .GE. -26.5).AND.
     &    (B .LT. -25.5)) THEN                     ! MBM56
         IHFLG=1
         IF (D .LE. 0.099) AV=0.
         IF ((D .GT. 0.099) .AND. (D .LE. 0.101)) AV=800.*(D-0.099)
         IF (D .GT. 0.101) AV=1.6
         RETURN
      ENDIF
C
C     Selected values from Keto & Myers
C
      IF ((L .GE. 290.0) .AND. (L .LT. 293.5) .AND. (B .GE. -31.3)
     &  .AND. (B .LT. -29.3)) THEN                 ! K293.2-30.9
         IHFLG=1
         IF (D .LE. 0.0995) AV=0.
         IF ((D .GT. .0995) .AND. (D .LE. .1005)) AV=400.*(D-0.0995)
         IF (D .GT. 0.1005) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 141.50) .AND. (L .LT. 142.50) .AND. (B .GE. 34.00)
     &  .AND. (B .LT. 39.00)) THEN                  ! K225.3-66.3
         IHFLG=1
         IF (D .LE. 0.0995) AV=0.
         IF ((D .GT. .0995) .AND. (D .LE. .1005)) AV=500.*(D-0.0995)
         IF (D .GT. 0.1005) AV=0.5
         RETURN
      ENDIF
      IF ((L .GE. 298.1) .AND. (L .LT. 298.5) .AND. (B .GE. -29.1)
     &  .AND. (B .LT. -28.6)) THEN                  ! K298.3-28.9
         IHFLG=1
         IF (D .LE. 0.0997) AV=0.
         IF ((D .GT. .0997) .AND. (D .LE. .1003)) AV=666.*(D-0.0997)
         IF (D .GT. 0.1003) AV=0.3996
         RETURN
      ENDIF
      IF ((L .GE. 292.85) .AND. (L .LT. 293.25) .AND. (B .GE. -31.0)
     &  .AND. (B .LT. -30.6)) THEN                  ! K292.1-31.0
         IHFLG=1
         IF (D .LE. 0.0995) AV=0.
         IF ((D .GT. .0995) .AND. (D .LE. .1005)) AV=400.*(D-0.0995)
         IF (D .GT. 0.1005) AV=0.4
         RETURN
      ENDIF
      IF ((L .GE. 291.7) .AND. (L .LT. 292.35) .AND. (B .GE. -31.1)
     &  .AND. (B .LT. -30.8)) THEN                  ! K292.1-31.0
         IHFLG=1
         IF (D .LE. 0.0995) AV=0.
         IF ((D .GT. .0995) .AND. (D .LE. .1005)) AV=600.*(D-0.0995)
         IF (D .GT. 0.1005) AV=0.6
         RETURN
      ENDIF
      IF ((L .GE. 294.6) .AND. (L .LT. 295.4) .AND. (B .GE.-36.4).AND.
     &    (B .LT. -36.0)) THEN                      ! K295.3-36.2
         IHFLG=1
         IF (D .LE. 0.0996) AV=0.
         IF ((D .GT. .0996) .AND. (D .LE. .1004)) AV=750.*(D-0.0996)
         IF (D .GT. 0.1004) AV=0.6
         RETURN
      ENDIF
      IF ((L .GE. 308.6) .AND. (L .LT. 309.2) .AND. (B .GE.-24.6).AND.
     &    (B .LT. -24.1)) THEN                      ! K308.9-24.3
         IHFLG=1
         IF (D .LE. 0.0993) AV=0.
         IF ((D .GT. .0993) .AND. (D .LE. .1007)) AV=833.*(D-0.0993)
         IF (D .GT. 0.1007) AV=0.4998
         RETURN
      ENDIF
      IF ((L .GE. 300.4) .AND. (L .LT. 300.9) .AND. (B .GE.-32.2).AND.
     &    (B .LT. -31.3)) THEN                       ! K300.8-31.5
         IHFLG=1
         IF (D .LE. 0.0992) AV=0.
         IF ((D .GT. .0992) .AND. (D .LE. .1008)) AV=500.*(D-0.0992)
         IF (D .GT. 0.1008) AV=0.9
         RETURN
      ENDIF
      IF ((L .GE. 316.25) .AND. (L .LT. 316.8) .AND. (B .GE.20.5).AND.
     &    (B .LT. 21.3)) THEN                       ! K316.5+21.0
         IHFLG=1
         IF (D .LE. 0.0994) AV=0.
         IF ((D .GT. .0994) .AND. (D .LE. .1006)) AV=667.*(D-0.0994)
         IF (D .GT. 0.1006) AV=0.8004
         RETURN
      ENDIF
C
C     Selected values from Odenwald
C
      IF ((L .GE. 354.5) .AND. (L .LT. 356.7) .AND. (B .GE. 20.0)
     &  .AND. (B .LT. 22.5)) THEN                   ! G354+24
         IHFLG=1
         IF (D .LE. 0.13) AV=0.
         IF ((D .GT. .13) .AND. (D .LE. .132)) AV=1750.*(D-0.13)
         IF (D .GT. 0.132) AV=3.5
         RETURN
      ENDIF
      IF ((L .GE. 355.0) .AND. (L .LT. 356.0) .AND. (B .GE. 22.5)
     &  .AND. (B .LT. 24.2)) THEN                   ! G354+24
         IHFLG=1
         IF (D .LE. 0.13) AV=0.
         IF ((D .GT. .13) .AND. (D .LE. .132)) AV=6350.*(D-0.13)
         IF (D .GT. 0.132) AV=12.7
         RETURN
      ENDIF
C
C     Lupus Complex from Hughes et al.
C
      IF ((L .GE. 335.0) .AND. (L .LT. 341.0) .AND. (B .GE. 7.0)
     &  .AND. (B .LT. 17.0)) THEN                   ! Lupus complex
         IHFLG=1
         IF (D .LE. 0.14) AV=0.
         IF ((D .GT. .14) .AND. (D .LE. .24)) AV=15.*(D-0.14)
         IF (D .GT. 0.14) AV=1.5
         RETURN
      ENDIF
C
C     Taurus Cloud from Kenyon et al.
C
      IF ((L .GE. 168.0) .AND. (L .LT. 176.0) .AND. (B .GE. -17.5)
     &  .AND. (B .LT. -12.0)) THEN                  ! Taurus cloud
         IHFLG=1
         IF (D .LE. 0.14) AV=0.
         IF ((D .GT. .14) .AND. (D .LE. .18)) AV=38.*(D-0.14)
         IF (D .GT. 0.18) AV=1.9
         RETURN
      ENDIF
C
C     Perseus Clouds from Cernis (1990,1993)
C
      IF ((L .GE. 160.0) .AND. (L .LT. 161.0) .AND. (B .GE. -17.5)
     &  .AND. (B .LT. -16.5)) THEN                  ! Perseus clouds
         IHFLG=1
         SAV=0.1
         IF (D .LE. 0.34) AV=0.
         IF ((D .GT. .34) .AND. (D .LE. .39)) AV=20.*(D-0.34)
         IF (D .GT. 0.39) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 159.0) .AND. (L .LT. 160.5) .AND. (B .GE. -20.0)
     &  .AND. (B .LT. -17.5)) THEN                  ! Perseus clouds
         IHFLG=1
         SAV=0.3
         IF (D .LE. 0.16) AV=0.
         IF ((D .GT. .16) .AND. (D .LE. .21)) AV=10.*(D-0.16)
         IF ((D .GT. .21) .AND. (D .LE. .39)) AV=0.5
         IF ((D .GT. .39) .AND. (D .LE. .44)) AV=0.5+40.*(D-0.39)
         IF (D .GT. 0.39) SAV=0.5
         IF (D .GT. 0.44) AV=2.5
         RETURN
      ENDIF
      IF ((L .GE. 156.5) .AND. (L .LT. 160.5) .AND. (B .GE. -22.0)
     &  .AND. (B .LT. -20.0)) THEN                  ! Perseus clouds
         IHFLG=1
         SAV=0.3
         IF (D .LE. 0.16) AV=0.
         IF ((D .GT. .16) .AND. (D .LE. .21)) AV=10.*(D-0.16)
         IF ((D .GT. .21) .AND. (D .LE. .39)) AV=0.5
         IF ((D .GT. .39) .AND. (D .LE. .44)) AV=0.5+40.*(D-0.39)
         IF (D .GT. 0.39) SAV=0.5
         IF (D .GT. 0.44) AV=2.5
         RETURN
      ENDIF
C
C     Corona Australis Complex from Rossano, AJ, 83, 234 (1978).
C
      IF ((L .GE. 359.0) .AND. (L .LT. 360.0) .AND. (B .GE. -20.5)
     &  .AND. (B .LT. -17.0)) THEN                  ! Corona Aus. complex
         IHFLG=1
         SAV=1.0
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. .15) .AND. (D .LE. .16)) AV=200.*(D-0.15)
         IF (D .GT. 0.16) AV=2.0
         RETURN
      ENDIF
      IF ((L .GE. 0.0) .AND. (L .LT. 0.5) .AND. (B .GE. -20.5)
     &  .AND. (B .LT. -17.0)) THEN                  ! Corona Aus. complex
         IHFLG=1
         SAV=1.0
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. .15) .AND. (D .LE. .16)) AV=200.*(D-0.15)
         IF (D .GT. 0.16) AV=2.0
         RETURN
      ENDIF
      IF ((L .GE. 359.0) .AND. (L .LT. 360.0) .AND. (B .GE. -22.5)
     &  .AND. (B .LT. -20.5)) THEN                  ! Corona Aus. complex
         IHFLG=1
         SAV=0.4
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. .15) .AND. (D .LE. .16)) AV=80.*(D-0.15)
         IF (D .GT. 0.16) AV=0.8
         RETURN
      ENDIF
      IF ((L .GE. 0.0) .AND. (L .LT. 1.5) .AND. (B .GE. -21.5)
     &  .AND. (B .LT. -20.5)) THEN                  ! Corona Aus. complex
         IHFLG=1
         SAV=0.6
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. .15) .AND. (D .LE. .16)) AV=120.*(D-0.15)
         IF (D .GT. 0.16) AV=1.2
         RETURN
      ENDIF
      IF ((L .GE. 0.0) .AND. (L .LT. 3.0) .AND. (B .GE. -22.5)
     &  .AND. (B .LT. -21.5)) THEN                  ! Corona Aus. complex
         IHFLG=1
         SAV=0.5
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. .15) .AND. (D .LE. .16)) AV=100.*(D-0.15)
         IF (D .GT. 0.16) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 0.5) .AND. (L .LT. 3.0) .AND. (B .GE. -24.0)
     &  .AND. (B .LT. -22.5)) THEN                  ! Corona Aus. complex
         IHFLG=1
         SAV=0.4
         IF (D .LE. 0.15) AV=0.
         IF ((D .GT. .15) .AND. (D .LE. .16)) AV=80.*(D-0.15)
         IF (D .GT. 0.16) AV=0.8
         RETURN
      ENDIF
C
C     Ophiuchus/Scorpio Complex from Rossano, AJ, 83, 241 (1978).
C
      IF ((L .GE. 348.0) .AND. (L .LT. 354.0) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  18.0)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.7
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=210.*(D-0.20)
         IF (D .GT. 0.21) AV=2.1
         RETURN
      ENDIF
      IF ((L .GE. 354.0) .AND. (L .LT. 358.5) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  20.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.6
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=180.*(D-0.20)
         IF (D .GT. 0.21) AV=1.8
         RETURN
      ENDIF
      IF ((L .GE. 358.5) .AND. (L .LT. 360.0) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  23.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 0.0) .AND. (L .LT. 1.5) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  23.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 1.5) .AND. (L .LT. 2.5) .AND. (B .GE. 18.0)
     &  .AND. (B .LT.  21.0)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.4
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=120.*(D-0.20)
         IF (D .GT. 0.21) AV=1.2
         RETURN
      ENDIF
      IF ((L .GE. 2.5) .AND. (L .LT. 4.5) .AND. (B .GE. 18.0)
     &  .AND. (B .LT.  20.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 2.5) .AND. (L .LT. 5.0) .AND. (B .GE. 20.5)
     &  .AND. (B .LT.  23.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.27
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=80.*(D-0.20)
         IF (D .GT. 0.21) AV=0.8
         RETURN
      ENDIF
      IF ((L .GE. 4.0) .AND. (L .LT. 6.5) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  18.0)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.6
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=180.*(D-0.20)
         IF (D .GT. 0.21) AV=1.8
         RETURN
      ENDIF
      IF ((L .GE. 1.5) .AND. (L .LT. 4.0) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  18.0)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.27
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=80.*(D-0.20)
         IF (D .GT. 0.21) AV=0.8
         RETURN
      ENDIF
      IF ((L .GE. 4.5) .AND. (L .LT. 5.0) .AND. (B .GE. 18.0)
     &  .AND. (B .LT.  20.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.2
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=60.*(D-0.20)
         IF (D .GT. 0.21) AV=0.6
         RETURN
      ENDIF
      IF ((L .GE. 5.0) .AND. (L .LT. 6.5) .AND. (B .GE. 20.5)
     &  .AND. (B .LT.  22.0)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.2
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=60.*(D-0.20)
         IF (D .GT. 0.21) AV=0.6
         RETURN
      ENDIF
      IF ((L .GE. 5.0) .AND. (L .LT. 6.5) .AND. (B .GE. 22.0)
     &  .AND. (B .LT.  23.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 5.0) .AND. (L .LT. 11.0) .AND. (B .GE. 18.0)
     &  .AND. (B .LT.  20.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.5
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=150.*(D-0.20)
         IF (D .GT. 0.21) AV=1.5
         RETURN
      ENDIF
      IF ((L .GE. 6.0) .AND. (L .LT. 10.0) .AND. (B .GE. 20.5)
     &  .AND. (B .LT.  23.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 11.0) .AND. (L .LT. 18.0) .AND. (B .GE. 15.0)
     &  .AND. (B .LT.  20.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
      IF ((L .GE. 18.0) .AND. (L .LT. 20.0) .AND. (B .GE. 19.0)
     &  .AND. (B .LT.  20.5)) THEN                  ! Oph./Sco. complex
         IHFLG=1
         SAV=0.33
         IF (D .LE. 0.20) AV=0.
         IF ((D .GT. .20) .AND. (D .LE. .21)) AV=100.*(D-0.20)
         IF (D .GT. 0.21) AV=1.0
         RETURN
      ENDIF
C
C     Orion Complex estimated from Kutner et al., ApJ, 215, 521 (1977).
C
      IF ((L .GE. 196.5) .AND. (L .LT. 198.0) .AND. (B .GE. -27.5)
     &  .AND. (B .LT.  -26.5)) THEN                  ! Orion complex
         IHFLG=1
         SAV=0.9
         IF (D .LE. 0.50) AV=0.
         IF ((D .GT. .50) .AND. (D .LE. .51)) AV=180.*(D-0.50)
         IF (D .GT. 0.51) AV=1.8
         RETURN
      ENDIF
      IF ((L .GE. 198.0) .AND. (L .LT. 200.0) .AND. (B .GE. -30.0)
     &  .AND. (B .LT.  -27.5)) THEN                  ! Orion complex
         IHFLG=1
         SAV=1.5
         IF (D .LE. 0.50) AV=0.
         IF ((D .GT. .50) .AND. (D .LE. .51)) AV=300.*(D-0.50)
         IF (D .GT. 0.51) AV=3.0
         RETURN
      ENDIF
      IF ((L .GE. 199.5) .AND. (L .LT. 201.0) .AND. (B .GE. -34.0)
     &  .AND. (B .LT.  -32.0)) THEN                  ! Orion complex
         IHFLG=1
         SAV=1.5
         IF (D .LE. 0.50) AV=0.
         IF ((D .GT. .50) .AND. (D .LE. .51)) AV=300.*(D-0.50)
         IF (D .GT. 0.51) AV=3.0
         RETURN
      ENDIF
      IF ((L .GE. 201.0) .AND. (L .LT. 207.5) .AND. (B .GE. -34.0)
     &  .AND. (B .LT.  -32.0)) THEN                  ! Orion complex
         IHFLG=1
         SAV=0.4
         IF (D .LE. 0.50) AV=0.
         IF ((D .GT. .50) .AND. (D .LE. .51)) AV=80.*(D-0.50)
         IF (D .GT. 0.51) AV=0.8
         RETURN
      ENDIF
C
C
C     All other high-galactic latitude regions:
C     If abs(b) <= 34.75 degrees (less than 1 114-pc scale-height above
C       the Galactic plane so that the chances of observing multiple
C       clouds in the line-of-sight are negligible) then background
C       extinction is zero; otherwise do not calculate background
C       extinction in unsampled areas.
C
      IF (ABS(B) .GE. 34.75) AV=0.
      SAV=0.24
      RETURN
      END

      SUBROUTINE COMBINE(A,SA,AVT,SAVT,JMAX,IHFLG)
      IMPLICIT REAL*4 (A-H,L-Z)
      DIMENSION A(5),SA(5)
C------------------------------------------------------------------------
C     The COMBINE subroutine takes the values for AV(1to5) and SAV(1to5)
C       and combines them together to obtain the average extinction AVT and 
C       its associated error SAVT.              
C     NOTE: JMAX = 5 because there are 5 subroutines that will pass data 
C       points to COMBINE. This value will change in the future as more 
C       subroutines are created.             
C------------------------------------------------------------------------
      JMAX=5
      AVT=0.
      SAVT=0.
      AVP=0.
      SAVP=0.
      SAVP2=0.
      DO 1 J=1,JMAX
      IF ((A(J) .EQ. -99.) .OR. (SA(J) .EQ. -99.)) THEN
         JMAX=JMAX-1.
         AVP=0.
         SAVP=0.
      ELSE
         AVP=A(J)
         SAVP=SA(J)
      ENDIF
      AVT=AVT+AVP
      SAVP2=SAVP2+SAVP
C----------------------------------------------------------------------
C    AVP, SAVP, SAVP2 are only parts of the equations needed to
C       find the average extinction and the (statistically dependent) error.  
C    NOTE: AVP and SAVP are used in the loop to calculate the individual  
C       contributions to AVT and SAVP2.       
C    NOTE: JMAX is only a counter that tells you if all the data points
C       have values or not. It also tells us how many subroutines were 
C       used in calculating the final extinction and error.                
C----------------------------------------------------------------------
      IF ((JMAX .EQ. 0.)) THEN
         AVT=-99.
         SAVT=-99.
         RETURN
      ENDIF
  1   CONTINUE  
      AVT=AVT/JMAX
      SAVT=SAVP2/JMAX
C----------------------------------------------------------------------
C
C     Override average value for well-defined high-latitude clouds
C
C----------------------------------------------------------------------
      IF ((A(5) .GE. 0.) .AND. (IHFLG .EQ. 1))THEN
         AVT=A(5)
         SAVT=SA(5)
      ENDIF
      RETURN
      END

      SUBROUTINE CORRECT(L,B,D,AVC)
      IMPLICIT REAL*4 (A-H,L-Z)
      BETA=0.114
      IF (B .EQ. 0.) THEN
         IF ((D .LT. 1.) .OR. (D .GE. 5.)) THEN
            AVC=0.
            RETURN
         ENDIF
         IF ((D .GE. 1.0) .AND. (D .LT. 1.5)) THEN
            AVC=0.73*(D-1.0)
            RETURN
         ENDIF
         IF ((D .GE. 1.5) .AND. (D .LT. 2.0)) THEN
            AVC=0.365+1.09*(D-1.5)
            RETURN
         ENDIF
         IF ((D .GE. 2.0) .AND. (D .LT. 2.5)) THEN
            AVC=0.91+0.67*(D-2.0)
            RETURN
         ENDIF
         IF ((D .GE. 2.5) .AND. (D .LT. 3.0)) THEN
            AVC=1.245+0.55*(D-2.5)
            RETURN
         ENDIF
         IF ((D .GE. 3.0) .AND. (D .LT. 3.5)) THEN
            AVC=1.52+0.36*(D-3.0)
            RETURN
         ENDIF
         IF ((D .GE. 3.5) .AND. (D .LT. 4.0)) THEN
            AVC=1.70+0.30*(D-3.5)
            RETURN
         ENDIF
         IF ((D .GE. 4.0) .AND. (D .LT. 4.5)) THEN
            AVC=1.85+0.18*(D-4.0)
            RETURN
         ENDIF
         IF ((D .GE. 4.5) .AND. (D .LT. 5.0)) THEN
            AVC=1.94+0.10*(D-4.5)
            RETURN
         ENDIF
      ELSE
         DP=D*COS(B*3.141592/180.)
         IF ((DP .LT. 1.) .OR. (DP .GE. 5.)) THEN
            AVC=0.
            RETURN
         ENDIF
         IF ((DP .GE. 1.0) .AND. (DP .LT. 1.5)) THEN
            AVC=EYFCN(0.73,1.0,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 1.5) .AND. (DP .LT. 2.0)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 2.0) .AND. (DP .LT. 2.5)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,2.0,B,BETA)
     &         +EYFCN(0.67,2.0,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 2.5) .AND. (DP .LT. 3.0)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,2.0,B,BETA)
     &         +EYFCN(0.67,2.0,2.5,B,BETA)+EYFCN(0.55,2.5,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 3.0) .AND. (DP .LT. 3.5)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,2.0,B,BETA)
     &         +EYFCN(0.67,2.0,2.5,B,BETA)+EYFCN(0.55,2.5,3.0,B,BETA)
     &         +EYFCN(0.36,3.0,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 3.5) .AND. (DP .LT. 4.0)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,2.0,B,BETA)
     &         +EYFCN(0.67,2.0,2.5,B,BETA)+EYFCN(0.55,2.5,3.0,B,BETA)
     &         +EYFCN(0.36,3.0,3.5,B,BETA)+EYFCN(0.30,3.5,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 4.0) .AND. (DP .LT. 4.5)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,2.0,B,BETA)
     &         +EYFCN(0.67,2.0,2.5,B,BETA)+EYFCN(0.55,2.5,3.0,B,BETA)
     &         +EYFCN(0.36,3.0,3.5,B,BETA)+EYFCN(0.30,3.5,4.0,B,BETA)
     &         +EYFCN(0.18,4.0,D,B,BETA)
            RETURN
         ENDIF
         IF ((DP .GE. 4.5) .AND. (DP .LT. 5.0)) THEN
            AVC=EYFCN(0.73,1.0,1.5,B,BETA)+EYFCN(1.09,1.5,2.0,B,BETA)
     &         +EYFCN(0.67,2.0,2.5,B,BETA)+EYFCN(0.55,2.5,3.0,B,BETA)
     &         +EYFCN(0.36,3.0,3.5,B,BETA)+EYFCN(0.30,3.5,4.0,B,BETA)
     &         +EYFCN(0.18,4.0,4.5,B,BETA)+EYFCN(0.10,4.5,D,B,BETA)
            RETURN
         ENDIF
      ENDIF
      RETURN
      END
