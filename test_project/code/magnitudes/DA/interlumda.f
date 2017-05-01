C***********************************************************************
C     TODO: rewrite      
      subroutine interlumda(tcool,mass,Z,lum,teff,logg,c1,c2,c3,c4,c5,
     &           table)
      use external_types
C=======================================================================
C
C     This subroutine interpolates the luminosity of DA WD according
C     to its mass, cooling time and metallicity.
C     It uses the cooling tables by Althaus et al. (2009) and
C     Renedo et al. (2010)
C
C     Modifications in 07/2012 (ER Cojocaru)
C
C-------------------------------------------------------------------
C     Input parameters:
C       tcool: cooling time
C       mass: mass of WD
C-------------------------------------------------------------------
C     Output parameters
C       lum: luminosity
C       teff: effective temperature [K]
C       logg: logarithm of the superficial gravity
C       c1,c2,c3,c4,c5: UBVRI colors 
C=======================================================================
      use external_types
      implicit double precision (a-h,m,o-z)
      
C     ---   Parameters   ---
      integer nrow,model,modlog
      parameter (nrow=650)
      double precision zet1,zet2,zet3,zet4
      parameter (zet1=0.001)
      parameter (zet2=0.01)
      parameter (zet3=0.03)
      parameter (zet4=0.06)

C     ---   Declaration of variables   ---
      integer numberOfMassesWithColors,numberOfMassesWithCoolSeq_1,
     &        numberOfMassesWithCoolSeq_2,numberOfMassesWithCoolSeq_3,
     &        numberOfMassesWithCoolSeq_4
      integer ntrkda(10),numberOfRows_1(7),numberOfRows_2(10),
     &        numberOfRows_3(8),numberOfRows_4(8)
      double precision tcool,mass,Z,lum,teff,logg,c1,c2,c3,c4,c5
      double precision Z1,lum1,teff1,logg1
      double precision Z2,lum2,teff2,logg2
      double precision massOfWD(10),massOfWD_1(7),massOfWD_2(10),
     &                 massOfWD_3(8),massOfWD_4(8)
      double precision tprewdda1(7),tprewdda2(10),tprewdda3(8)
      double precision tprewdda4(8),coolingTimes_1(7,nrow),
     &                 coolingTimes_2(10,nrow)
      double precision coolingTimes_3(8,nrow),coolingTimes_4(8,nrow)
      double precision luminosity_1(7,nrow),
     &                 effectiveTemperature_1(7,nrow),
     &                 gravitationalAcceleration_1(7,nrow)
      double precision luminosity_2(10,nrow),
     &                 effectiveTemperature_2(10,nrow),
     &                 gravitationalAcceleration_2(10,nrow)
      double precision luminosity_3(8,nrow),
     &                 effectiveTemperature_3(8,nrow),
     &                 gravitationalAcceleration_3(8,nrow)
      double precision luminosity_4(8,nrow),
     &                 effectiveTemperature_4(8,nrow),
     &                 gravitationalAcceleration_4(8,nrow)
      double precision luminosity(10,nrow),color_U(10,nrow),
     &                 color_B(10,nrow),color_V(10,nrow)
      double precision color_R(10,nrow),color_I(10,nrow)

      TYPE(FileGroupInfo),DIMENSION(11) :: table

C     ---   Commons   ---
C     TODO: delete this
      common /nums/ numberOfMassesWithColors,
     &              numberOfMassesWithCoolSeq_1,
     &              numberOfMassesWithCoolSeq_2,
     &              numberOfMassesWithCoolSeq_3,
     &              numberOfMassesWithCoolSeq_4
      common /nums2/ ntrkda,numberOfRows_1,numberOfRows_2,
     &               numberOfRows_3,numberOfRows_4
      common /masses/ massOfWD,massOfWD_1,massOfWD_2,massOfWD_3,
     &                massOfWD_4
      common /datrks1/ coolingTimes_1,luminosity_1,
     &                 effectiveTemperature_1,
     &                 gravitationalAcceleration_1
      common /datrks2/ coolingTimes_2,luminosity_2,
     &                 effectiveTemperature_2,
     &                 gravitationalAcceleration_2
      common /datrks3/ coolingTimes_3,luminosity_3,
     &                 effectiveTemperature_3,
     &                 gravitationalAcceleration_3
      common /datrks4/ coolingTimes_4,luminosity_4,
     &                 effectiveTemperature_4,
     &                 gravitationalAcceleration_4
      common /datprewd/ tprewdda1,tprewdda2,tprewdda3,tprewdda4
      common /colors/ luminosity,color_U,color_B,color_V,color_R,color_I

      
      model=0
      modlog=0
C     QUESTION: what does it mean?
C     miro entre que valores de metalicidad se encuentra el Z de entrada
C     calculo lum,teff,logg para esas valores Z1 y Z2 y después hago la 
C     Interpolación para Z
      if(Z.ge.zet1.AND.Z.lt.zet2) then
        Z1=zet1
        call interp(model,modlog,tcool,mass,table(1)%ncol,
     &       table(1)%ntrk,table(1)%coolingTime,table(1)%prevTime,
     &       table(1)%mass,
     &       table(1)%luminosity,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,table(1)%ncol,
     &       table(1)%ntrk,table(1)%coolingTime,table(1)%prevTime,
     &       table(1)%mass,
     &       table(1)%effTemp,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,table(1)%ncol,
     &       table(1)%ntrk,table(1)%coolingTime,table(1)%prevTime,
     &       table(1)%mass,
     &       table(1)%gravAcc,logg1)
        Z2=zet2
        call interp(model,modlog,tcool,mass,
     &       table(2)%ncol,table(2)%ntrk,table(2)%coolingTime,
     &       table(2)%prevTime,table(2)%mass,table(2)%luminosity,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       table(2)%ncol,table(2)%ntrk,table(2)%coolingTime,
     &       table(2)%prevTime,table(2)%mass,table(2)%effTemp,
     &       teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       table(2)%ncol,table(2)%ntrk,table(2)%coolingTime,
     &       table(2)%prevTime,table(2)%mass,
     &       table(2)%gravAcc,logg2)
      end if
      
      if(Z.ge.zet2.AND.Z.lt.zet3) then
        Z1=zet2
        call interp(model,modlog,tcool,mass,
     &       table(2)%ncol,table(2)%ntrk,table(2)%coolingTime,
     &       table(2)%prevTime,table(2)%mass,table(2)%luminosity,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       table(2)%ncol,table(2)%ntrk,table(2)%coolingTime,
     &       table(2)%prevTime,table(2)%mass,table(2)%effTemp,
     &       teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       table(2)%ncol,table(2)%ntrk,table(2)%coolingTime,
     &       table(2)%prevTime,table(2)%mass,
     &       table(2)%gravAcc,logg1)
        Z2=zet3
        call interp(model,modlog,tcool,mass,
     &       table(3)%ncol,table(3)%ntrk,table(3)%coolingTime,
     &       table(3)%prevTime,table(3)%mass,table(3)%luminosity,
     &       lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       table(3)%ncol,table(3)%ntrk,table(3)%coolingTime,
     &       table(3)%prevTime,table(3)%mass,table(3)%effTemp,
     &       teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       table(3)%ncol,table(3)%ntrk,table(3)%coolingTime,
     &       table(3)%prevTime,table(3)%mass,
     &       table(3)%gravAcc,logg2)
      end if
      
      if(Z.ge.zet3.AND.Z.le.zet4) then
        Z1=zet3
        call interp(model,modlog,tcool,mass,
     &       table(3)%ncol,table(3)%ntrk,table(3)%coolingTime,
     &       table(3)%prevTime,table(3)%mass,table(3)%luminosity,
     &       lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       table(3)%ncol,table(3)%ntrk,table(3)%coolingTime,
     &       table(3)%prevTime,table(3)%mass,table(3)%effTemp,
     &       teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       table(3)%ncol,table(3)%ntrk,table(3)%coolingTime,
     &       table(3)%prevTime,table(3)%mass,
     &       table(3)%gravAcc,logg1)
        Z2=zet4
        call interp(model,modlog,tcool,mass,
     &       table(4)%ncol,table(4)%ntrk,table(4)%coolingTime,
     &       table(4)%prevTime,table(4)%mass,table(4)%luminosity,
     &       lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       table(4)%ncol,table(4)%ntrk,table(4)%coolingTime,
     &       table(4)%prevTime,table(4)%mass,table(4)%effTemp,
     &       teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       table(4)%ncol,table(4)%ntrk,table(4)%coolingTime,
     &       table(4)%prevTime,table(4)%mass,
     &       table(4)%gravAcc,logg2)
      end if
      
      lum=lum1+(lum2-lum1)*(Z-Z1)/(Z2-Z1)
      teff=teff1+(teff2-teff1)*(Z-Z1)/(Z2-Z1)
      logg=logg1+(logg2-logg1)*(Z-Z1)/(Z2-Z1)
      
C     interpolation of colors, here there is no metalicity
      call intermag(mass,lum,table(11)%ncol,table(11)%ntrk,
     &     table(11)%mass,
     &     table(11)%luminosity,table(11)%color_U,table(11)%color_B,
     &     table(11)%color_V,table(11)%color_R,table(11)%color_I,
     &     c1,c2,c3,
     &     c4,c5)

      return
      end
C***********************************************************************
