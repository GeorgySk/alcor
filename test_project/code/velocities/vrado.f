C***********************************************************************
C     TODO:rewrite      
      subroutine vrado(u,v,w)
C***********************************************************************
C     This subroutine calculates the heliocentric velocities, starting
C     from the proper motions in galactic coordinates, making zero
C     the component of radial velocity
C***********************************************************************
      implicit real (a-h,m,o-z)
       
      integer numberOfStars,i,numberOfWDs
      real k,xcb,xsb,xcl,xsl,r
      real a1,a2,b1,b2,b3,c1,c2,c3
C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (k=4.740)
C     ---   Dimensiones   --- 
      real longitude_proper_motion(numberOfStars),
     &     latitude_proper_motion(numberOfStars),
     &     radial_velocity(numberOfStars)
      real rgac(numberOfStars)
      double precision :: lgac(numberOfStars),
     &                    bgac(numberOfStars)
      real flagOfWD(numberOfStars)
      real u(numberOfWDs),v(numberOfWDs),w(numberOfWDs)
      integer disk_belonging(numberOfStars)

C    ---   Commons  ---
      common /lb/ lgac,bgac
      common /paral/ rgac
      common /mopro/ latitude_proper_motion,
     &               longitude_proper_motion,
     &               radial_velocity
      common /index/ flagOfWD,numberOfWDs,disk_belonging

C     ---  Calculating the heliocentric velocity  
C          making zero the radial velocity ---
      do 1 i=1,numberOfWDs
        xcb=cos(bgac(i))
        xsb=sin(bgac(i))
        xcl=cos(lgac(i))
        xsl=sin(lgac(i)) 
        r=rgac(i)*1000.0
        a1=-k*xcb*xsl
        b1=-k*xsb*xcl
        c1=0.0
        u(i) = a1 * longitude_proper_motion(i) * r
     &         + b1 * latitude_proper_motion(i) * r
     &         + c1 * radial_velocity(i)
        a2=k*xcb*xcl
        b2=-k*xsb*xsl
        c2=0.0
        v(i) = a2 * longitude_proper_motion(i) * r
     &         + b2 * latitude_proper_motion(i) * r
     &         + c2 * radial_velocity(i)
        b3=k*xcb
        c3=0.0
        w(i) = b3 * latitude_proper_motion(i) * r
     &         + c3 * radial_velocity(i)
1     continue
         
      return
      end
C***********************************************************************
