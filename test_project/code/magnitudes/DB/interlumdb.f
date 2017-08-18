      subroutine interlumdb(cooling_time,
     &                      wd_mass,
     &                      Z,
     &                      luminosity,
     &                      c1,
     &                      c2,
     &                      c3,
     &                      c4,
     &                      c5,
     &                      effective_temperature,
     &                      log_surface_gravity,
     &                      table)
C     Interpolating luminosity of DB WDs using cooling tables by 
C     Althaus et al. (2009) for metallicities 0.001, 0.01, 0.06
C         TODO: rename
C         c1,c2,c3,c4,c5: Johnson colors Johnson (UBVRI)
      use external_types
      implicit none

C     TODO: find out the meanings of variables and rename them
      integer, parameter :: NROWB = 400,
     &                      NROWB2 = 60  
      integer :: numberOfSequencesInGroup_1,
     &           numberOfSequencesInGroup_2,
     &           numberOfSequencesInGroup_3,
     &           numberOfSequences,
     &           model,
     &           modlog   
      real :: cooling_time,
     &        wd_mass,
     &        Z,
     &        luminosity,
     &        effective_temperature,
     &        log_surface_gravity,
     &        Z1,
     &        lum1,
     &        teff1,
     &        logg1,
     &        Z2,
     &        lum2,
     &        teff2,
     &        logg2,
     &        c1,
     &        c2,
     &        c3,
     &        c4,
     &        c5,
     &        zet1,
     &        zet2,
     &        zet3
      integer ::  vectorOfPointsNumberOfSeq_1(7),
     &            vectorOfPointsNumberOfSeq_2(9),
     &            vectorOfPointsNumberOfSeq_3(9),
     &            numberOfPointsInSequence(7)
      real :: vectorOfMasses_1(7),
     &        vectorOfPreviousTimes_1(7),
     &        matrixOfCoolingTimes_1(7,NROWB),
     &        matrixOfLuminosities_1(7,NROWB),
     &        matrixOfEffectiveTemperatures_1(7,NROWB),
     &        matrixOfLog_g_1(7,NROWB),
     &        vectorOfMasses_2(9),
     &        vectorOfPreviousTimes_2(9),
     &        matrixOfCoolingTimes_2(9,NROWB),
     &        matrixOfLuminosities_2(9,NROWB),
     &        matrixOfEffectiveTemperatures_2(9,NROWB),
     &        matrixOfLog_g_2(9,NROWB),
     &        vectorOfMasses_3(9),
     &        vectorOfPreviousTimes_3(9),
     &        matrixOfCoolingTimes_3(9,NROWB),
     &        matrixOfLuminosities_3(9,NROWB),
     &        matrixOfEffectiveTemperatures_3(9,NROWB),
     &        matrixOfLog_g_3(9,NROWB),
     &        massSequence(7),
     &        luminosityDB(7,NROWB2),
     &        colorDB_U(7,NROWB2),
     &        colorDB_B(7,NROWB2),
     &        colorDB_V(7,NROWB2),
     &        colorDB_R(7,NROWB2),
     &        colorDB_I(7,NROWB2)

      TYPE(FileGroupInfo),DIMENSION(11) :: table

C     TODO: delete this
      common /dbnums/ numberOfSequences,
     &                numberOfSequencesInGroup_1,
     &                numberOfSequencesInGroup_2,
     &                numberOfSequencesInGroup_3
      common /dbnums2/ vectorOfPointsNumberOfSeq_1,
     &                 vectorOfPointsNumberOfSeq_2,
     &                 vectorOfPointsNumberOfSeq_3,
     &                 numberOfPointsInSequence
      common /massesdb/ massSequence,
     &                  vectorOfMasses_1,
     &                  vectorOfMasses_2,
     &                  vectorOfMasses_3
      common /dbtrks1/ matrixOfLuminosities_1,
     &                 matrixOfEffectiveTemperatures_1,
     &                 matrixOfLog_g_1,
     &                 matrixOfCoolingTimes_1
      common /dbtrks2/ matrixOfLuminosities_2,
     &                 matrixOfEffectiveTemperatures_2,
     &                 matrixOfLog_g_2,
     &                 matrixOfCoolingTimes_2
      common /dbtrks3/ matrixOfLuminosities_3,
     &                 matrixOfEffectiveTemperatures_3,
     &                 matrixOfLog_g_3,
     &                 matrixOfCoolingTimes_3
      common /dbtprewd/ vectorOfPreviousTimes_1,
     &                  vectorOfPreviousTimes_2,
     &                  vectorOfPreviousTimes_3
      common /dbcolors/ luminosityDB,
     &                  colorDB_U,
     &                  colorDB_B,
     &                  colorDB_V,
     &                  colorDB_R,
     &                  colorDB_I
     
      zet1 = 0.001
      zet2 = 0.01
      zet3 = 0.06
     
C     Determinating luminosities, magnitudes for every star
      model = 0
      modlog = 0
      
C     Checking between what values of metallicity is the Z of input
C     calculating luminosity, effective_temperature, log_surface_gravity
C     for those values Z1 and Z2 and then interpolating for Z 
      if (Z >= zet1 .and. Z < zet2) then
          Z1 = zet1
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(5)%ncol,table(5)%ntrk,
     &         table(5)%coolingTime,table(5)%prevTime,
     &         table(5)%mass,table(5)%luminosity,lum1)
          modlog = 1
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(5)%ncol,table(5)%ntrk,
     &         table(5)%coolingTime,table(5)%prevTime,
     &         table(5)%mass,table(5)%effTemp,teff1)
          modlog = 0
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(5)%ncol,table(5)%ntrk,
     &         table(5)%coolingTime,table(5)%prevTime,
     &         table(5)%mass,table(5)%gravAcc,logg1)
  
          Z2 = zet2
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(6)%ncol,table(6)%ntrk,
     &         table(6)%coolingTime,table(6)%prevTime,
     &         table(6)%mass,table(6)%luminosity,lum2)
          modlog = 1
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(6)%ncol,table(6)%ntrk,
     &         table(6)%coolingTime,table(6)%prevTime,
     &         table(6)%mass,table(6)%effTemp,teff2)
          modlog = 0
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(6)%ncol,table(6)%ntrk,
     &         table(6)%coolingTime,table(6)%prevTime,
     &         table(6)%mass,table(6)%gravAcc,logg2)
      end if
      
      if (Z >= zet2 .and. Z < zet3) then
          Z1 = zet2
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(6)%ncol,table(6)%ntrk,
     &         table(6)%coolingTime,table(6)%prevTime,
     &         table(6)%mass,table(6)%luminosity,lum1)
          modlog = 1
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(6)%ncol,table(6)%ntrk,
     &         table(6)%coolingTime,table(6)%prevTime,
     &         table(6)%mass,table(6)%effTemp,teff1)
          modlog = 0
          call interp(model,modlog,cooling_time,wd_mass,
     &         table(6)%ncol,table(6)%ntrk,
     &         table(6)%coolingTime,table(6)%prevTime,
     &         table(6)%mass,table(6)%gravAcc,logg1)
          Z2 = zet3
          call interp(model,modlog,cooling_time,wd_mass,table(7)%ncol,
     &         table(7)%ntrk,table(7)%coolingTime,
     &         table(7)%prevTime,table(7)%mass,
     &         table(7)%luminosity,lum2)
          modlog = 1
          call interp(model,modlog,cooling_time,wd_mass,table(7)%ncol,
     &         table(7)%ntrk,table(7)%coolingTime,
     &         table(7)%prevTime,table(7)%mass,
     &         table(7)%effTemp,teff2)
          modlog = 0
          call interp(model,modlog,cooling_time,wd_mass,table(7)%ncol,
     &         table(7)%ntrk,table(7)%coolingTime,
     &         table(7)%prevTime,table(7)%mass,table(7)%gravAcc,
     &         logg2)
      end if
      
C     Interpolation 3D in function of Z
      luminosity = lum1 + (lum2 - lum1) * (Z - Z1) / (Z2 - Z1)
      effective_temperature = teff1 + (teff2 - teff1) * (Z - Z1) 
     &                                / (Z2 - Z1)
      log_surface_gravity = logg1 + (logg2 - logg1) * (Z - Z1) 
     &                              / (Z2 - Z1)
      
C     Colors interpolation, there is no metallicity here
      call intermag(wd_mass,luminosity,table(10)%ncol,
     &     table(10)%ntrk,table(10)%mass,table(10)%luminosity,
     &     table(10)%color_U,table(10)%color_B,table(10)%color_V,
     &     table(10)%color_R,table(10)%color_I,c1,c2,c3,c4,c5)

      end subroutine
