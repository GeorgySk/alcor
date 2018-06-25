      subroutine velh(iseed,
     &                sample_stars_count,
     &                geometry)
C         Сalculating heliocentrical velocities. 
C         Velocity dispersion is based on scale height
C           TODO: get rid of this, implement another RNG
C           iseed: random number generator parameter (dummy)
          implicit none
          
C         TODO: this is numberOfStars, take it up as const
          integer, parameter :: MAX_STARS_COUNT = 6000000
C         A, B: Oort constants in km/sKpc (values by Kerr and 
C              Lynden-Bell 1986)
          real, parameter :: oort_const_A = 14.4,
     &                       oort_const_B = -12.8,
     &                       ro = 10.5,
     &                       l = 5.5,
     &                       sigo = 80.0,
     &                       sigm = 145.0,
     &                       pi = 4.0 * atan(1.0),
     &                       vc = 220.0,
     &                       gamma = 3.4
          double precision, parameter :: solar_galactocentric_distance 
     &                                   = 8.5d0
          integer :: iseed,
     &               sample_stars_count,
     &               i,
     &               k,
     &               numberOfWDs,
     &               disk_belonging(MAX_STARS_COUNT),
     &               flagOfWD(MAX_STARS_COUNT)
          real :: peculiar_solar_velocity_u,
     &            peculiar_solar_velocity_v,
     &            peculiar_solar_velocity_w,
     &            uop,
     &            vop,
     &            gasdev,
     &            uom,
     &            vom,
C         Heliocentric velocities, B3 system
     &            uu(MAX_STARS_COUNT),
     &            vv(MAX_STARS_COUNT),
     &            ww(MAX_STARS_COUNT),
     &            velocities_dispersions(3),
     &            stellar_galactocentric_distance,
     &            xx,
     &            sigo2,
     &            sigm2,
     &            sigr2,
     &            dsigr2dr,
     &            sigt2,
     &            sigr,
     &            sigt,
     &            vr,
     &            vth,
     &            vph,
     &            delta,
     &            sindel,
     &            cosdel,
     &            vx,
     &            vy,
     &            vz,
     &            sini,
     &            cosi
          character(len = 6) :: geometry
          double precision :: coordinate_R(MAX_STARS_COUNT),
     &                        coordinate_Theta(MAX_STARS_COUNT),
     &                        coordinate_Zcylindr(MAX_STARS_COUNT),
     &                        lgac(MAX_STARS_COUNT),
     &                        bgac(MAX_STARS_COUNT)
         real :: starBirthTime(MAX_STARS_COUNT),
     &           m(MAX_STARS_COUNT)
          
          common /vel/ uu, vv, ww
          common /coorcil/ coordinate_R,
     &                     coordinate_Theta,
     &                     coordinate_Zcylindr
          common /index/ flagOfWD,
     &                   numberOfWDs,
     &                   disk_belonging
          common /lb/ lgac, bgac
          common /tm/ starBirthTime, m
    
C         Peculiar solar velocities, values from Anguiano et al. 2017
C         TODO: find out why signs are different from values in Python
          peculiar_solar_velocity_u = -11.0
          peculiar_solar_velocity_v = -12.0
          peculiar_solar_velocity_w = -7.0
    
C         Making the transfer of z, r, theta
          k = 0      
          do i = 1, sample_stars_count
              if (flagOfWD(i) == 1) then
                  k = k + 1
                  if (geometry == 'cones') then
                      coordinate_R(k) = coordinate_R(i)
                      coordinate_Theta(k) = coordinate_Theta(i)
                  end if
              end if
          end do
    
          do i = 1, numberOfWDs
C             TODO: find these constants and make variables for them
              if (disk_belonging(i) == 1) then
                  velocities_dispersions(1) = 32.4
                  velocities_dispersions(2) = 23.0
                  velocities_dispersions(3) = 18.1 
              else if (disk_belonging(i) == 2) then
                  velocities_dispersions(1) = 50.0
                  velocities_dispersions(2) = 56.0
                  velocities_dispersions(3) = 34.0
              else if (disk_belonging(i) /= 3) then
                  write(6, *) "Error in velh:"
                  write(6, *) "    Number of WD's =", numberOfWDs
                  write(6, *) "    index =", i
                  write(6, *) "    disk_belonging =", disk_belonging(i)
                  stop
              end if
            
C             TODO: find out what is going on here (see article 
C             Simulating Gaia performances on white dwarfs by Torres)
              if (disk_belonging(i) == 3) then
C                 Radial and tangential dispersions
                  stellar_galactocentric_distance = real(sqrt(
     &                coordinate_R(i) ** 2 
     &                + coordinate_Zcylindr(i) ** 2))
                  xx = (stellar_galactocentric_distance - ro) / l
                  sigo2 = sigo * sigo
                  sigm2 = sigm * sigm
                  sigr2 = sigo2 + sigm2 * (0.5 - (atan(xx) / pi))
                  dsigr2dr = -(1.0 / l) * sigm2 / ((1.0 + xx * xx) * pi)
                  sigt2 = 0.5 * vc * vc + (1.0 - gamma / 2.0) * sigr2 
     &                    + (stellar_galactocentric_distance / 2.0) 
     &                      * dsigr2dr
                  sigr = sqrt(sigr2)
                  sigt = sqrt(sigt2)
                  sigr = vc / sqrt(2.0)
                  sigt = vc / sqrt(2.0)

C                 Spherical coordinates
                  vr = sigr * gasdev(iseed)
                  vth = sigt * gasdev(iseed)
                  vph = sigt * gasdev(iseed)

C                 Cartesian coordinates
                  delta = real(pi - lgac(i) - coordinate_Theta(i))
                  sindel = sin(delta)
                  cosdel = cos(delta)
                  vx = -cosdel * vr + sindel * vph
                  vy = sindel * vr + cosdel * vph
                  vz = vth

                  sini = sin(real(lgac(i)))
                  cosi = cos(real(lgac(i)))

                  uu(i) = vz * cosi + vx * sini
                  vv(i) = -vz * sini + vx * cosi
                  ww(i) = vy

                  vv(i) = vv(i) - vc

                  uu(i) = uu(i) + peculiar_solar_velocity_u
                  vv(i) = vv(i) + peculiar_solar_velocity_v
                  ww(i) = ww(i) + peculiar_solar_velocity_w
              else
C                 Calling function of gaussian distribution
C                 TODO: find out the meaning of uop
                  uop = uom(coordinate_R(i),
     &                      coordinate_Theta(i),
     &                      oort_const_A, 
     &                      oort_const_B,
     &                      solar_galactocentric_distance,
     &                      peculiar_solar_velocity_u)
                  uu(i) = velocities_dispersions(1) * gasdev(iseed) +uop
C                 TODO: find out the meaning of vop
                  vop = vom(coordinate_R(i),
     &                      coordinate_Theta(i),
     &                      oort_const_A, 
     &                      oort_const_B,
     &                      solar_galactocentric_distance,
     &                      peculiar_solar_velocity_v)
C                 TODO: find out what this 120.0 is
                  vv(i) = velocities_dispersions(2) * gasdev(iseed) +vop 
     &                    - velocities_dispersions(1) 
     &                      * velocities_dispersions(1)
     &                      / 120.0
                  ww(i) = velocities_dispersions(3) * gasdev(iseed) 
     &                    + peculiar_solar_velocity_w
              end if
          end do
      end subroutine


C     TODO: give a better name
      function uom(coordinate_R, 
     &             coordinate_Theta, 
     &             oort_const_A,
     &             oort_const_B, 
     &             solar_galactocentric_distance, 
     &             peculiar_solar_velocity_u)
C         Calculating uom taking into account the effect of gal. rotation
          implicit none
      
          real :: oort_const_A,
     &            oort_const_B,
     &            peculiar_solar_velocity_u,
     &            uom
          double precision :: coordinate_R, 
     &                        coordinate_Theta, 
     &                        solar_galactocentric_distance
C         TODO: what is going on here? 
          uom = real(peculiar_solar_velocity_u 
     &               + ((3.0 - (2.0 * coordinate_R) 
     &                         / solar_galactocentric_distance) 
     &                  * oort_const_A - oort_const_B) * coordinate_R 
     &                 * sin(coordinate_Theta))
      end function


C     TODO: give a better name
      function vom(coordinate_R,
     &             coordinate_Theta,
     &             oort_const_A,
     &             oort_const_B,
     &             solar_galactocentric_distance,
     &             peculiar_solar_velocity_v)
C         Calculating vom taking into account the effect of gal. rotation
          implicit none
          
          real :: oort_const_A, 
     &            oort_const_B,
     &            peculiar_solar_velocity_v,
     &            vom
          double precision :: coordinate_R, 
     &                        coordinate_Theta,
     &                        solar_galactocentric_distance
    
          vom = real(peculiar_solar_velocity_v 
     &               + ((3.0 - (2.0 * coordinate_R) 
     &                         / solar_galactocentric_distance) 
     &                  * oort_const_A - oort_const_B) * coordinate_R 
     &                 * cos(coordinate_Theta) 
     &               - (oort_const_A - oort_const_B) 
     &                 * solar_galactocentric_distance)
      end function
