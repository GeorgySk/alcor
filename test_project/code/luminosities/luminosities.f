      subroutine lumx(numberOfStarsInSample,
     &                thick_disk_age,
     &                halo_age)
C         Determining what stars are WDs and calculating cooling time 
C         and luminosity for it
          implicit none

          integer, parameter :: MAX_STARS_COUNT = 6000000
C         Mass of Oxygen-Neon White Dwarf (ONe WD)
          real, parameter :: O_NE_WD_MASS = 1.14
          integer :: numberOfStarsInSample,
     &               numberOfWDs,
C                    Number of WDs of Oxygen-Neon type
     &               wd_ONe_count,
     &               i,
     &               k
          real :: thin_disk_age,
     &            parameterIFMR,
     &            fractionOfDB,
     &            parameterIMF,
     &            burst_age,
     &            thick_disk_age,
     &            halo_age,
     &            max_age

C       TODO: change to logical - is_WD
C       flagOfWD: 0 - it's not WD, 1 - it's a WD
C       TODO: rename
C       m: mass in the main sequence
          real :: starBirthTime(MAX_STARS_COUNT),
     &            main_sequence_lifetime(MAX_STARS_COUNT),
     &            coolingTime(MAX_STARS_COUNT),
     &            luminosityOfWD(MAX_STARS_COUNT),
     &            massOfWD(MAX_STARS_COUNT),
     &            metallicityOfWD(MAX_STARS_COUNT),
     &            effTempOfWD(MAX_STARS_COUNT),
     &            m(MAX_STARS_COUNT)
          integer disk_belonging(MAX_STARS_COUNT),
     &            flagOfWD(MAX_STARS_COUNT)
          double precision :: coordinate_R(MAX_STARS_COUNT),
     &                        coordinate_Theta(MAX_STARS_COUNT),
     &                        coordinate_Zcylindr(MAX_STARS_COUNT)

C         TODO: give names with one style
          common /tm/ starBirthTime, m
          common /enanas/ luminosityOfWD,
     &                    massOfWD,
     &                    metallicityOfWD,
     &                    effTempOfWD
          common /index/ flagOfWD,
     &                   numberOfWDs,
     &                   disk_belonging
          common /cool/ coolingTime
          common /param/ fractionOfDB,
     &                   thin_disk_age,
     &                   parameterIMF,
     &                   parameterIFMR,
     &                   burst_age
         common /coorcil/ coordinate_R,
     &                     coordinate_Theta,
     &                     coordinate_Zcylindr

          max_age = max(thin_disk_age, thick_disk_age, halo_age)

          numberOfWDs = 0
          wd_ONe_count = 0

C         Deciding if the star is a WD
          do i = 1, numberOfStarsInSample
              flagOfWD(i) = 0
C             Progenitor star that generates a ONe WD: 8.5 <M_MS < 10.5
C             WD of CO: m_WD <1.14; of ONe: m_wd>1.14
              if (m(i) <= 10.5) then
                  if (disk_belonging(i) == 3) then
                      metallicityOfWD(i) = 0.001
                  else
                      metallicityOfWD(i) = 0.01
                  end if
C                 Calculating the lifetime in the main sequence
C                 TODO: resolve the problem with halo stars
                  call get_main_sequence_lifetime(
     &                     m(i), 
     &                     metallicityOfWD(i), 
     &                     main_sequence_lifetime(i))
                  coolingTime(i) = max_age - starBirthTime(i) 
     &                             - main_sequence_lifetime(i)
                  if (coolingTime(i) > 0.0) then
C                     Initial-to-Final Mass Relation (IFMR)
                      call mmswd(m(i), massOfWD(i))
C                     Using Z solar z = 0.01
                      massOfWD(i) = parameterIFMR * massOfWD(i)
                      if (massOfWD(i) <= 1.4) then 
                          flagOfWD(i) = 1
                          numberOfWDs = numberOfWDs + 1
                          if (massOfWD(i) > O_NE_WD_MASS) then
                              wd_ONe_count = wd_ONe_count + 1
                          end if
                      else
                          flagOfWD(i) = 0
                      end if
                  else
                      flagOfWD(i) = 0
                  end if
              else
                  flagOfWD(i) = 0
              end if
          end do

          write (6,*) ' Number of WDs: ', numberOfWDs
          write (6,*) ' Number of ONe: ', wd_ONe_count
          if (numberOfWDs /= 0) then
              write (6,*) ' ONe percentage: ',
     &                    100.0 * float(wd_ONe_count) 
     &                    / float(numberOfWDs), '%'
          end if
    
          k = 0
C         Making the transfer
          do i = 1, numberOfStarsInSample
              if (flagOfWD(i) == 1) then
                  k = k + 1
                  coolingTime(k) = coolingTime(i)
                  massOfWD(k) = massOfWD(i)
                  metallicityOfWD(k) = metallicityOfWD(i)
                  starBirthTime(k) = starBirthTime(i)
                  disk_belonging(k) = disk_belonging(i)
                  coordinate_Zcylindr(k) = coordinate_Zcylindr(i)
              end if
          end do
    
          write(6,*) '      Total number of WDs=', numberOfWDs
          write(6,*) '      WDs of ONe=', wd_ONe_count
      end subroutine


      subroutine get_main_sequence_lifetime(stellar_mass, 
     &                                      metallicity,
     &                                      main_sequence_lifetime)
C         Calculates lifetime in the main sequence for a given 
C         metallicity Z € [0.01, 0.001] for standart helium content 
C         according to model by Leandro & Renedo et al.(2010)
C         Data in solar masses and Gyr
          implicit none

          real :: stellar_mass,
     &            metallicity,
     &            main_sequence_lifetime,
     &            mms(10),
     &            tms(10),
     &            mms2(7),
     &            tms2(7),
     &            pen,
     &            tsol,
     &            tsub
          integer :: k

C         Table of values Z solar --
          mms(1) = 1.00
          mms(2) = 1.50    
          mms(3) = 1.75   
          mms(4) = 2.00   
          mms(5) = 2.25 
          mms(6) = 2.50 
          mms(7) = 3.00 
          mms(8) = 3.50  
          mms(9) = 4.00 
          mms(10) = 5.00   

C         Althaus priv. comm X = 0.725, Y = 0.265
          tms(1) = 8.614
          tms(2) = 1.968
          tms(3) = 1.249  
          tms(4) = 0.865  
          tms(5) = 0.632
          tms(6) = 0.480
          tms(7) = 0.302
          tms(8) = 0.226 
          tms(9) = 0.149
          tms(10) = 0.088   

C         Interpolating
          if (stellar_mass < mms(1)) then
C             Mass less than the first, linear extrapolation, 2 last 
C             points 
              pen = (tms(2) - tms(1)) / (mms(2) - mms(1))
              tsol = pen * stellar_mass + (tms(1) - pen * mms(1))
          else
              if (stellar_mass > mms(10)) then
C                 Mass greater than the last, taking the fraction last 
C                 point
                  tsol = (mms(10) / stellar_mass) * tms(10)
              else
C                 TODO: eliminate goto
                  k = 1
  1               k = k + 1
                  if (stellar_mass < mms(k)) then
                      pen = (tms(k) - tms(k - 1)) 
     &                      / (mms(k) - mms(k - 1))
                      tsol = pen * stellar_mass 
     &                       + (tms(k) - pen * mms(k))
                  else
                      goto 1
                  end if 
              end if
          end if

C         Tabla of values Z Sub-Solar
          mms2(1) = 0.85
          mms2(2) = 1.00    
          mms2(3) = 1.25   
          mms2(4) = 1.50   
          mms2(5) = 1.75 
          mms2(6) = 2.00 
          mms2(7) = 3.00 

C         Althaus priv. comm X = 0.752, Y = 0.247
          tms2(1) = 10.34
          tms2(2) = 5.756
          tms2(3) = 2.623  
          tms2(4) = 1.412  
          tms2(5) = 0.905
          tms2(6) = 0.639
          tms2(7) = 0.245

C         Interpolating
          if (stellar_mass < mms2(1)) then
C             Mass less than the first, linear extrapolation, 2 last 
C             points 
              pen = (tms2(2) - tms2(1)) / (mms2(2) - mms2(1))
              tsub = pen * stellar_mass + (tms2(1) - pen * mms2(1))
          else
              if (stellar_mass > mms2(7)) then
C                 Masa greater than the last, extrapolating 2 last 
C                 points
                  tsub = (mms(7) / stellar_mass) * tms(7)
              else  
C                 TODO: eliminate goto 
                  k = 1
  2               k = k + 1
                  if (stellar_mass < mms2(k)) then
                      pen = (tms2(k) - tms2(k - 1)) 
     &                      / (mms2(k) - mms2(k - 1))
                      tsub = pen * stellar_mass 
     &                       + (tms2(k) - pen * mms2(k))
                  else
                      goto 2
                  end if 
              end if
          end if

C         Interpolating for the value of Z, z solar 10
          main_sequence_lifetime = tsub + ((tsol - tsub) 
     &                                     / (0.01 - 0.001)) 
     &                                    * (metallicity - 0.001)
      end subroutine
      

C     TODO: find out the meaning of the name of the function mmswd
      subroutine mmswd(mass,massOfWD)
C         IFMR according to model by Catalán et al.2008
C         combination with the model by Iben for Mi>6Mo
          implicit none

          real :: mass,
     &            massOfWD
      
          if (mass < 2.7) then
              massOfWD = 0.096 * mass + 0.429
          else if ((mass >= 2.7) .and. (mass <= 6.0)) then
C             Small correction for continuity
              massOfWD = 0.137 * mass + 0.3183
          else
C             Slope of Iben + continuity in Mi=6Mo
              massOfWD = 0.1057 * mass + 0.5061
          end if
      end subroutine
