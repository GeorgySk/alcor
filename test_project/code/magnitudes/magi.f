C     TODO: rewrite      
      subroutine magi(fractionOfDB,table)
      use external_types
C=======================================================================
C
C     This subroutine calculates ltc,cbv,cvi,cvr,cuv visual absolute 
C     and apparent magnitude of the WDs.
C     
C     Created by S.Torres
C     Introduced metallicity 08/2012 (ER Cojocaru)
C=======================================================================
      implicit real (a-h,m,o-z)

      integer numberOfStars,iseed,i,numberOfWDs,in
      real lum,teff,xlog,c1,c2,c3,c4,c5,n1,n2,n3,n4,n5
      real xg,xug,xgr,xri,xiz,xgi,fractionOfDB
      real mone

      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

      real luminosityOfWD(numberOfStars),
     &                 massOfWD(numberOfStars),
     &                 metallicityOfWD(numberOfStars),
     &                 effTempOfWD(numberOfStars)
      real flagOfWD(numberOfStars)
      real v(numberOfStars)
      real :: coolingTime(numberOfStars)
      integer typeOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)

C     TODO: no need to keep these
      real :: ugriz_ug(numberOfStars),
     &        ugriz_gr(numberOfStars),
     &        ugriz_ri(numberOfStars),
     &        ugriz_iz(numberOfStars)

      real :: UB(numberOfStars), BV(numberOfStars), 
     &        VRR(numberOfStars), RI(numberOfStars),
     &        ugriz_g_apparent(numberOfStars),
     &        rgac(numberOfStars)
      double precision :: lgac(numberOfStars),
     &                    bgac(numberOfStars)
      real :: AVT,
     &        SAVT,
     &        AVC,
     &        AV(5),
     &        SAV(5),
     &        extinction,
     &        ugriz_u,
     &        ugriz_g,
     &        ugriz_r,
     &        ugriz_i,
     &        ugriz_z,
     &        pi,
     &        ugriz_u_apparent,
     &        ugriz_r_apparent,
     &        ugriz_i_apparent,
     &        ugriz_z_apparent,
     &        ugriz_u_apparent_w_error,
     &        ugriz_g_apparent_w_error,
     &        ugriz_r_apparent_w_error,
     &        ugriz_i_apparent_w_error,
     &        ugriz_z_apparent_w_error
      integer :: JMAX

      TYPE(FileGroupInfo),DIMENSION(11) :: table

      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD
      common /index/ flagOfWD,numberOfWDs,disk_belonging
      common /paral/ rgac
C     TODO: no need to keep these
      common /photo/ ugriz_ug, ugriz_gr, ugriz_ri, ugriz_iz, 
     &               ugriz_g_apparent
      common /johnson/ v
      common /cool/ coolingTime
      common /indexdb/ typeOfWD
      common /ubvri/ UB, BV, VRR, RI
      common /lb/ lgac,bgac

      pi = 4.0 * atan(1.0)

      n1=0
      n2=0
      n3=0
      n4=0
      n5=0

C     ---  Interpolating Mv, luminosity, colors and other variables
C          from coolingTime and the mwd  ---
C     ---  Start DO on all the stars
      do i=1,numberOfWDs
C       ---  ATENTION! choosing only if .lt.1.1!!!  ---
C       ---  Start IF mass <1.4  ----
        if(massOfWD(i).le.1.4) then
C         ---  IF CO core ---
          if(massOfWD(i).lt.mone) then  
C           --- Atention We put the "old" ones to 0.6Msol ---
C           --- Distribucion DA/DB ---
            call dbd_fid(iseed,fractionOfDB,in)
C           --- End of distribution ---
C           ---  IF DA ---
            if(in.eq.0) then
              typeOfWD(i)=0
              n1=n1+1
              call interlumda(coolingTime(i),massOfWD(i),
     &             metallicityOfWD(i),lum,teff,xlog,c1,c2,c3,c4,c5,
     &             table)
C           ---  ELSE DB  ---
            else
              n3=n3+1
              typeOfWD(i)=1    
              call interlumdb(coolingTime(i),massOfWD(i),
     &             metallicityOfWD(i),lum,c1,c2,c3,c4,c5,teff,xlog,
     &             table)
              if(teff.lt.6000) n5=n5+1
            end if
C           ---  END IF DB/NON-DB
C         ---  ELSE ONe ---
          else
            n2=n2+1
            typeOfWD(i)=2
            call interlumone(coolingTime(i),massOfWD(i),lum,c1,c2,c3,c4,
     &           c5,teff,xlog)
          end if
C         ---  END IF CO/ONe ---
          if(teff.lt.6000) n4=n4+1
          luminosityOfWD(i)=-lum
          effTempOfWD(i)=teff            
          V(i)=c3
C         TODO: no need to keep indexes
C         TODO: rename VRR to VR
          UB(i)=c1-c2
          BV(i)=c2-c3
          VRR(i)=c3-c4
          RI(i)=c4-c5
          call chanco(V(i),UB(i),BV(i),VRR(i),RI(i),xg,xug,xgr,xri,xiz,
     &                xgi)
          ugriz_u = xug + xg
          ugriz_g = xg
          ugriz_r = xg - xgr
          ugriz_i = xg - xgi
          ugriz_z = xg - xgi - xiz

          call extinct(real(lgac(i) * 180.0 / pi),
     &                 real(bgac(i) * 180.0 / pi),
     &                 real(rgac(i)),
     &                 AVT,SAVT,AVC,JMAX,AV,SAV)
          extinction = AVT + AVC

          ugriz_u = ugriz_u + 1.579 * extinction
          ugriz_g = ugriz_g + 1.161 * extinction
          ugriz_r = ugriz_r + 0.843 * extinction
          ugriz_i = ugriz_i + 0.639 * extinction
          ugriz_z = ugriz_z + 0.453 * extinction

          ugriz_u_apparent = ugriz_u - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                              + 3.0)
          ugriz_g_apparent(i) = ugriz_g - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                              + 3.0)
          ugriz_r_apparent = ugriz_r - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                              + 3.0)
          ugriz_i_apparent = ugriz_i - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                              + 3.0)
          ugriz_z_apparent = ugriz_z - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                              + 3.0)

          call errfot(ugriz_u_apparent,
     &                ugriz_u_apparent_w_error,
     &                1)
          call errfot(ugriz_g_apparent(i),
     &                ugriz_g_apparent_w_error,
     &                2)
          call errfot(ugriz_r_apparent,
     &                ugriz_r_apparent_w_error,
     &                3)
          call errfot(ugriz_i_apparent,
     &                ugriz_i_apparent_w_error,
     &                4)
          call errfot(ugriz_z_apparent,
     &                ugriz_z_apparent_w_error,
     &                5) 

C         TODO: figure out what to do with commons
          ugriz_g_apparent(i) = ugriz_g_apparent_w_error

C         TODO: this is another way to calculate phot.error. do smth
C           ugriz_u_apparent_w_error = ugriz_u_apparent 
C      &                               + 0.1 * gasdev(iseed)
C           ugriz_g_apparent_w_error = ugriz_g_apparent(i) 
C      &                               + 0.025 * gasdev(iseed)
C           ugriz_r_apparent_w_error = ugriz_r_apparent 
C      &                               + 0.025 * gasdev(iseed)
C           ugriz_i_apparent_w_error = ugriz_i_apparent 
C      &                               + 0.05 * gasdev(iseed)
C           ugriz_z_apparent_w_error = ugriz_z_apparent 
C      &                               + 0.1 * gasdev(iseed)

          ugriz_ug(i) = ugriz_u_apparent_w_error 
     &                  - ugriz_g_apparent_w_error
          ugriz_gr(i) = ugriz_g_apparent_w_error 
     &                  - ugriz_r_apparent_w_error
          ugriz_ri(i) = ugriz_r_apparent_w_error 
     &                  - ugriz_i_apparent_w_error
          ugriz_iz(i) = ugriz_i_apparent_w_error 
     &                  - ugriz_z_apparent_w_error

C           ugriz_ug(i) = ugriz_u
C      &                  - ugriz_g
C           ugriz_gr(i) = ugriz_g
C      &                  - ugriz_r
C           ugriz_ri(i) = ugriz_r
C      &                  - ugriz_i
C           ugriz_iz(i) = ugriz_i
C      &                  - ugriz_z

C         ---  Making V apparent magnitude ---
          V(i) = V(i) - 5.0 + 5.0 * (log10(rgac(i)) + 3.0)
C       ---  ELSE mass >= 1.4  --- EXPLOTA, exceeding Chandrasekar limit
        else
          typeOfWD(i) = 5
        end if
      end do

      write(*,*) "DA CO ",n1
      write(*,*) "DA ONe ",n2
      write(*,*) "DB",n3
      write(*,*) "<6000 DA", n4
      write(*,*) "<6000 DB", n5

      return
      end
