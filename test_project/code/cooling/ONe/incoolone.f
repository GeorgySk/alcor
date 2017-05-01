C***********************************************************************
C     TODO: rewrite 
      SUBROUTINE incoolone
C=======================================================================
C
C       Reading the cooling tables of ONe WD's
C
C       Revised in 27.09.07 by S. Torres      
C
C-----------------------------------------------------------------------
C       mtabone:  mass of the WD of ONE 
C       ltabone:  logarithm of the luminosity
C       mvtabone: visual absolute magnitude visual
C       lgtabone: logarithm of the cooling time in years
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      integer nrow,ncol,nrow2,ncol2
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13 

C   ---  Parameter ---
      parameter(ncol=6)
      parameter(ncol2=5)
      parameter(nrow=300)
      parameter(nrow2=300)

C   ---  Dimension  ---
      integer i,j,ird
      integer ndatsone(ncol),ndatsone2(ncol2)
      double precision mtabone(ncol),mtabone2(ncol2)
      double precision lgtabone(ncol,nrow),ltabone(ncol,nrow)
      double precision mvtabone(ncol,nrow),lgtetabone(ncol,nrow)
      double precision bvtabone(ncol,nrow),vitabone(ncol,nrow)
      double precision vrtabone(ncol,nrow),uvtabone(ncol,nrow)
      double precision lgrtabone(ncol2,nrow2),lgt2tabone(ncol2,nrow2)

C  ---  Commons ---
      common /fredone/ lgtabone,ltabone,mvtabone,lgtetabone
      common /fredone2/ mtabone,ndatsone
      common /colorsone/ bvtabone,vitabone,vrtabone,uvtabone
      common /newone/ lgrtabone,lgt2tabone
      common /newone2/ mtabone2,ndatsone2

C   ---  ird is the unit of reading archives ---
      ird=121
 
C   ---  Masses of each archive ---
      mtabONE(1)=1.06d0
      mtabONE(2)=1.10d0
      mtabONE(3)=1.16d0
      mtabONE(4)=1.20d0
      mtabONE(5)=1.24d0
      mtabONE(6)=1.28d0
      mtabone2(1)=1.06d0
      mtabone2(2)=1.10d0
      mtabone2(3)=1.16d0
      mtabone2(4)=1.20d0
      mtabone2(5)=1.28d0

C   ---   Reading the files ---
      do  i=1,ncol
        do j=1,nrow
          read(ird,*,end=2) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13
          ltabone(i,j)=a1
          mvtabone(i,j)=a12
          lgtabone(i,j)=a13
          lgtetabone(i,j)=a2
          bvtabone(i,j)=a3
          vitabone(i,j)=a9
          vrtabone(i,j)=a4
          uvtabone(i,j)=a10
        end do      
2       ndatsONE(i)=j-1  
        ird=ird+1
      end do

C      ---  Reading the data of log Teff and log g ---
      ird=127      
      do  i=1,ncol2
        do j=1,nrow2
          read(ird,*,end=3) a1,a2,a3,a4,a5,a6
C         ---  Reconverting radii in cm to radii in solar radius ---
          a3=10.0**a3
          a3=a3/(6.96d10)
          a3=dlog10(a3)
          lgrtabone(i,j)=a3
          lgt2tabone(i,j)=a5 
        end do      
3       continue       
        ndatsONE2(i)=j-1  
        ird=ird+1
      end do

      return
      end
C***********************************************************************
