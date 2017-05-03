C     TODO: give better names
      subroutine incoolea
C=======================================================================
C
C     This subroutine reads the cooling tables by Leandro Althaus
C
C-----------------------------------------------------------------------
C     Input parameters
C
C       none
C
C-----------------------------------------------------------------------
C     Output parameters
C
C       luminosity:   Luminosuty of the WD. [log(L/L0)]
C       mtrk:  Mass of the WD. [M0]
C       ttrk:  Cooling time. [Gyr]
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Declaration of variables   ---
      double precision mtrk,luminosity

C     ---   Parameters   ---
      parameter (ncol=3)
      parameter (nbank=3)
      parameter (nrow=900)

C     ---   Dimensions   ---
      dimension mtrk(ncol),ntrk(ncol,nbank)
      dimension ttrk(ncol,nrow,nbank),luminosity(ncol,nrow,nbank)
      dimension ginic(nbank)

C     ---   Commons   ---
C     NOTE: it is used only here
      common /tracks/ ginic,luminosity,mtrk,ttrk,ntrk
      
C     --- Values of the mass and G for A=-1.1d-11 ---
      mtrk(1)=0.52
      mtrk(2)=0.6
      mtrk(3)=1.0

      ginic(1)=1.020
      ginic(2)=1.050
      ginic(3)=1.100

C     ---   Initialization   ---
      irdr=299
 
C     ---   Reading the tables  Z=0.05 ---
      do 4 k=1,nbank 
        do 3 i=1,ncol
C         ---   Reading unit   ---
          irdr=irdr+1
C         ---   Reading the cooling curves   ---
          do 1 j=1,nrow
C           QUESTION: what does end=2 mean?            
            read(irdr,*,end=2) luminosity(i,j,k),a2,a3,a4,ttrk(i,j,k),
     &                         a6,a7,a8,a9,goverg
            ttrk(i,j,k)=10.0**(ttrk(i,j,k)-3.0)
    1     continue
    2     ntrk(i,k)=j-1
    3   continue 
    4 continue

      return
      end
C***********************************************************************
