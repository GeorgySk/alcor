C***********************************************************************
C           FUNCTION RAND
C***********************************************************************
        FUNCTION RAN(IDUMMY)
C
C     THIS IS AN ADAPTED VERSION OF SUBROUTINE RANECU WRITTEN BY
C     F. JAMES (COMPUT. PHYS. COMMUN. 60 (1990) 329-344, WHICH HAS
C     BEEN MODIFIED TO GIVE A SINGLE RANDOM NUMBER AT EACH CALL.
C     THE 'SEEDS' ISEED1 AND ISEED2 MUST BE INITIALIZED IN THE 
C     MAIN PROGRAM AND TRANSFERRED THROUGH THE NAMED COMMON BLOCK
C     /RSEED/.
        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real ran
      PARAMETER (USCALE=1.0D0/2.0D0**31)
      COMMON /RSEED/ ISEED1,ISEED2

      I1=ISEED1/53668
      ISEED1=40014*(ISEED1-I1*53668)-I1*12211
      
      IF (ISEED1.LT.0) ISEED2=ISEED1+2147483563
      
      I2=ISEED2/52774
      ISEED2=40692*(ISEED2-I2*52774)-I2*3791
      
      IF (ISEED2.LT.0) ISEED2=ISEED2+2147483399
      
      IZ=ISEED1-ISEED2
      
      IF(IZ.LT.1) IZ=IZ+2147483562
    
      RAN=IZ*USCALE

       
      RETURN
      END 
C***********************************************************************
        

            
C***********************************************************************
C     SUBROUTINE ODEINT
C     Subroutine for the numerical integration of a system of first  
C     order differential equations
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************      

      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,
     &           NBAD,DERIVS,RKQC,yscal,y,dydx)
      implicit double precision (a-h,o-z)
      PARAMETER (MAXSTP=10000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.E-30)
      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP(200),YP(10,200)
      DIMENSION YSTART(NVAR),YSCAL(Nvar),Y(Nvar),DYDX(Nvar)
      EXTERNAL DERIVS
      EXTERNAL RKQC
      X=X1
      H=dSIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE
      IF (KMAX.GT.0) XSAV=X-DXSAV*TWO
        DO 16 NSTP=1,MAXSTP
          CALL DERIVS(X,Y,DYDX,xpla,ypla)
          DO 12 I=1,NVAR
            YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
12        CONTINUE
          IF(KMAX.GT.0)THEN
            IF(DABS(X-XSAV).GT.DABS(DXSAV)) THEN
              IF(KOUNT.LT.KMAX-1)THEN
                KOUNT=KOUNT+1
                XP(KOUNT)=X
                DO 13 I=1,NVAR
                  YP(I,KOUNT)=Y(I)
13              CONTINUE
                XSAV=X
              ENDIF
            ENDIF
          ENDIF
          IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
          CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
          IF(HDID.EQ.H)THEN
            NOK=NOK+1
          ELSE
            NBAD=NBAD+1
          ENDIF
          IF((X-X2)*(X2-X1).GE.ZERO)THEN
            DO 14 I=1,NVAR
              YSTART(I)=Y(I)
14          CONTINUE
            IF(KMAX.NE.0)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO 15 I=1,NVAR
                YP(I,KOUNT)=Y(I)
15            CONTINUE
            ENDIF
            RETURN
          ENDIF
          IF(dABS(HNEXT).LT.HMIN) stop 'Stepsize smaller than minimum.'
            H=HNEXT
16      CONTINUE
        STOP 'Too many steps.'

      END
C***********************************************************************

  
C***********************************************************************
C     SUBROUTINE RKQC
C     Runge-Kutta Quality Control
C     Compares the results between one big step and two small steps and
C     decides the size of the next step according to how different the
C     results are
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      implicit double precision(a-h,o-z)
      PARAMETER (NMAX=10,FCOR=.0666666667,ONE=1.,SAFETY=0.9,
     &          ERRCON=6.E-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV) stop 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=dMAX1(ERRMAX,DABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      
      RETURN
      END
C***********************************************************************


C***********************************************************************
C     SUBROUTINE RK4
C     Runge-Kutta 4
C     Runge-Kutta step on a set of n differential equations. You input  
C     the values of the independent variables, and you get out new  
C     values which are stepped by a stepsize h (which can be positive or 
C     negative). The routine that has three calls to derivs
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************        
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      implicit double precision (a-h,o-z)
      
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      EXTERNAL DERIVS
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
C***********************************************************************        


C***********************************************************************
C     SUBROUTINE derivs
C     Sample derivatives routine for stiff 
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************              
      subroutine derivs(x,y,dydx)
C=======================================================================
C     This subroutine provides the values of the derivatives of the 
C     variables y_i
C-------------------------------------------------------------------
      implicit double precision (a-h,o-z)

C     ---   Dimensions   ---      
      dimension y(2),dydx(2)

C     ---   Calculating the derivatives    ---
      dydx(1)=y(2)
      call fuerza(y(1),f)
      dydx(2)=f
      
      return
      end     
C***********************************************************************   

      SUBROUTINE CHSTWO(BINS1,BINS2,NBINS,KNSTRN,DF,CHSQ,PROB)
      implicit double precision(a-h,o-z)
      integer knstrn
      DIMENSION BINS1(NBINS),BINS2(NBINS)
      DF=dfloat(NBINS-1-KNSTRN)
      CHSQ=0.
      DO 11 J=1,NBINS
        IF(BINS1(J).EQ.0..AND.BINS2(J).EQ.0.)THEN
          DF=DF-1.
        ELSE
          CHSQ=CHSQ+(BINS1(J)-BINS2(J))**2/(BINS1(J)+BINS2(J))
        ENDIF
11    CONTINUE
      PROB=GAMMQ(0.5*DF,0.5*CHSQ)
      RETURN
      END

      FUNCTION GAMMQ(A,X)
      implicit double precision (a-h,o-z)
      IF(X.LT.0..OR.A.LE.0.)read (*,*)
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
      ENDIF
      RETURN
      END
      
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      implicit double precision(a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.) read(*,*)
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(dABS(DEL).LT.dABS(SUM)*EPS)GO TO 1
11    CONTINUE
      write(6,*) 'A too large, ITMAX too small'
      read(*,*)
1     GAMSER=SUM*dEXP(-X+A*dLOG(X)-GLN)
      RETURN
      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      implicit double precision (a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=dFLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(dABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      write(6,*) 'A too large, ITMAX too small'
      read(*,*)
1     GAMMCF=dEXP(-X+A*dLOG(X)-GLN)*G
      RETURN
      END


      FUNCTION GAMMLN(XX)
      implicit double precision (a-h,o-z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*dLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+dLOG(STP*SER)
      RETURN
      END