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
      implicit real (a-h,m,o-z)

      integer numberOfStars,iseed,i,numberOfWDs,in
      
C     ---   Variables  ---
      real lum,teff,xlog,c1,c2,c3,c4,c5,n1,n2,n3,n4,n5
C       real UB,BV,VR,RI
      real xg,xug,xgr,xri,xiz,xgi,fractionOfDB
      real mone

C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

C     ---   Dimensions  ---
      real luminosityOfWD(numberOfStars),
     &                 massOfWD(numberOfStars),
     &                 metallicityOfWD(numberOfStars),
     &                 effTempOfWD(numberOfStars)
      real flagOfWD(numberOfStars)
      real rgac(numberOfStars)
      real g(numberOfWDs),v(numberOfStars)
      real :: coolingTime(numberOfStars)
      integer typeOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)

C       TODO: uncomment this after UBVRI diag-s plotted
      real :: ugriz_ug(numberOfStars),
     &                    ugriz_gr(numberOfStars),
     &                    ugriz_ri(numberOfStars),
     &                    ugriz_iz(numberOfStars)
C      &                    ugriz_g_apparent(numberOfStars)

      real :: UB(numberOfStars), BV(numberOfStars), 
     &                    VRR(numberOfStars), RI(numberOfStars),
     &                    ugriz_g_apparent(numberOfStars)

      TYPE(FileGroupInfo),DIMENSION(11) :: table

C     ---   Commons   ---
      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD
      common /index/ flagOfWD,numberOfWDs,disk_belonging
      common /paral/ rgac
C       TODO: uncomment this after UBVRI diag-s plotted
      common /photo/ ugriz_ug, ugriz_gr, ugriz_ri, ugriz_iz, 
     &               ugriz_g_apparent
      common /johnson/ v
      common /cool/ coolingTime
      common /indexdb/ typeOfWD
      common /ubvri/ UB, BV, VRR, RI
C      & , ugriz_g_apparent

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
C         TODO: delete indexes after UBVRI diag-s are plotted
C         TODO: rename VRR to VR also after diag-s
          UB(i)=c1-c2
          BV(i)=c2-c3
          VRR(i)=c3-c4
          RI(i)=c4-c5
          call chanco(V(i),UB(i),BV(i),VRR(i),RI(i),xg,xug,xgr,xri,xiz,
     &                xgi)
          g(i) = xg
C           TODO: uncomment this after UBVRI color diagrams are plotted
          ugriz_ug(i) = xug
          ugriz_gr(i) = xgr
          ugriz_ri(i) = xri
          ugriz_iz(i) = xiz
C         ---  Making g and V apparent magnitude ---
          ugriz_g_apparent(i) = g(i) - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                              + 3.0)
          V(i) = V(i) - 5.0 + 5.0 * (log10(rgac(i)) + 3.0)
C       ---  ELSE mass >= 1.4  --- EXPLOTA, exceeding Chandrasekar limit
        else
          typeOfWD(i) = 5
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
