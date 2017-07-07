C***********************************************************************
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
C
C-----------------------------------------------------------------------
C     Input parameters:
C       numberOfStarsInSample
C-----------------------------------------------------------------------
C     Output parameters
C       none
C=======================================================================
      implicit double precision (a-h,m,o-z)

      integer numberOfStars,iseed,i,numberOfWDs,in
      
C     ---   Variables  ---
      double precision lum,teff,xlog,c1,c2,c3,c4,c5,n1,n2,n3,n4,n5
      double precision UB,BV,VR,RI,xg,xug,xgr,xri,xiz,xgi,fractionOfDB
      double precision mone

C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

C     ---   Dimensions  ---
      double precision luminosityOfWD(numberOfStars),
     &                 massOfWD(numberOfStars),
     &                 metallicityOfWD(numberOfStars),
     &                 effTempOfWD(numberOfStars)
      double precision flagOfWD(numberOfStars)
      double precision rgac(numberOfStars)
      double precision g(numberOfWDs),v(numberOfStars)
      double precision coolingTime(numberOfStars)
      integer typeOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)

      double precision :: ugriz_ug(numberOfStars),
     &                    ugriz_gr(numberOfStars),
     &                    ugriz_ri(numberOfStars),
     &                    ugriz_iz(numberOfStars),
     &                    ugriz_g_apparent(numberOfStars)

      TYPE(FileGroupInfo),DIMENSION(11) :: table

C     ---   Commons   ---
      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD
      common /index/ flagOfWD,numberOfWDs,disk_belonging
      common /paral/ rgac
      common /photo/ ugriz_ug, ugriz_gr, ugriz_ri, ugriz_iz, 
     &               ugriz_g_apparent
      common /johnson/ v
      common /cool/ coolingTime
      common /indexdb/ typeOfWD

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
          UB=c1-c2
          BV=c2-c3
          VR=c3-c4
          RI=c4-c5
          call chanco(V(i),UB,BV,VR,RI,xg,xug,xgr,xri,xiz,xgi)
          g(i) = xg
          ugriz_ug(i) = xug
          ugriz_gr(i) = xgr
          ugriz_ri(i) = xri
          ugriz_iz(i) = xiz
C         ---  Making g and V apparent magnitude ---
          ugriz_g_apparent(i) = g(i)-5.0d0+5.0d0*(dlog10(rgac(i))+3.0d0)
          V(i)=V(i)-5.0d0+5.0d0*(dlog10(rgac(i))+3.0d0)
C       ---  ELSE mass >= 1.4  --- EXPLOTA, exceeding Chandrasekar limit
        else
          typeOfWD(i)=5
        end if
C       ---  END IF about WD mass ---
      end do
C     ---  END DO about all the stars ---

      write(*,*) "DA CO ",n1
      write(*,*) "DA ONe ",n2
      write(*,*) "DB",n3
      write(*,*) "<6000 DA", n4
      write(*,*) "<6000 DB", n5

      return
      end
C***********************************************************************
