C***********************************************************************
C     TODO: rewrite      
      subroutine polar(iseed,minimumSectorRadius,maximumSectorRadius,
     &           angleCoveringSector,radiusOfSector,
     &           solarGalactocentricDistance,scaleLength)
C=======================================================================
C
C     This subrutina generates the positiones, in polar coordinates, of
C     the stars.
C   
C     Revised in 22.09.07 by S. Torres
C
C-----------------------------------------------------------------------
C     Input parameters:
C     minimumSectorRadius: minimum radius of the sector; in Kpc from the 
C                          Galactic Center
C     maximumSectorRadius: maximum radius
C     angleCoveringSector: angle covering the sector in degrees;
C                          from the GC
C     QUESTION: did i fail to name this variable correctly?
C     radiusOfSector: radial distance to the Sun
C     solarGalactocentricDistance: galactocentric distance of the Sun
C     scaleLength: scale length
C     numberOfWDs: total number of WDs
C=======================================================================
      implicit real (a-h,m,o-z)

C     --- Declaration of variables  ---
      external ran       
      real ran
      integer numberOfWDs,iseed,j
      real dospi,minimumSectorRadius,maximumSectorRadius,
     &                 angleCoveringSector,radiusOfSector,
     &                 scaleLength,asr
      double precision :: solarGalactocentricDistance
      real drsun2,dist,pi,xmin,xmax,zzz,zzr,zzy,zz,xx,xc,yc
      real xcte,xinc
       
C     ---  Parameters  ---
      integer numberOfStars
      parameter (numberOfStars=6000000)
           
C     ---  Dimensions  ---
      double precision :: coordinate_R(numberOfStars),
     &                    coordinate_Theta(numberOfStars),
     &                    coordinate_Zcylindr(numberOfStars)
      real x(numberOfStars),y(numberOfStars)
      real flagOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)

C     ---  Commons  ---
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /plano/ x,y
      common /index/ flagOfWD,numberOfWDs,disk_belonging

C     ---  Inicialization of pi and sigma ---
      pi=4.0*atan(1.0)
      dospi=2.0*pi
      drsun2=radiusOfSector*radiusOfSector
      
C     --- Calculating the angle in the sector
C         -angleCoveringSector/2 y +angleCoveringSector/2 degrees
C      and radius between minimumSectorRadius and maximumSectorRadius ---
      asr=(angleCoveringSector*pi)/180.0
      xmax=(maximumSectorRadius*maximumSectorRadius)
      xmin=(minimumSectorRadius*minimumSectorRadius)
      xcte=xmax-xmin
      xinc=maximumSectorRadius-minimumSectorRadius
                
      do 2 j=1,numberOfWDs
3       coordinate_Theta(j)=asr*ran(iseed)-(asr/2)
        if (coordinate_Theta(j).lt.0.0) then
          coordinate_Theta(j)=dospi+coordinate_Theta(j)
        endif
31      zzz=minimumSectorRadius+xinc*ran(iseed)
        zzy=0.16*ran(iseed)
        zzr=exp(-zzz/scaleLength)
        if (zzy.gt.zzr) then
          goto 31
        else
        endif
        zz=(zzz-minimumSectorRadius)/xinc
        xx=xmin+(xcte*zz)
        coordinate_R(j)=sqrt(xx)  
        xc=coordinate_R(j)*cos(coordinate_Theta(j))
        yc=coordinate_R(j)*sin(coordinate_Theta(j))
        dist=((xc-solarGalactocentricDistance)*
     &       (xc-solarGalactocentricDistance)+yc*yc)
C       QUESTION: what does it mean?       
C       --- Sol no hay más que uno ---
        if (dist.gt.drsun2.or.dist.lt.0.0000015) then 
          goto 3
        endif
        x(j)=xc
        y(j)=yc
 2    continue

      return
      end
C***********************************************************************
