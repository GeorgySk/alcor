      subroutine generate_cone_stars(cone_height_longitude,
     &                               cone_height_latitude,
     &                               numberOfStarsInSample,
     &                               iseed,
     &                               kinematicModel)
          implicit none
          double precision :: cone_height_longitude,
     &                        cone_height_latitude
          integer :: iseed,
     &               numberOfStarsInSample,
     &               kinematicModel

          integer, parameter :: numberOfStars = 6000000
          double precision :: starBirthTime(numberOfStars),
     &                        m(numberOfStars),
     &                        coordinate_Theta(numberOfStars),
     &                        coordinate_R(numberOfStars),
     &                        coordinate_Zcylindr(numberOfStars),
     &                        heightPattern(numberOfStars),
     &                        flagOfWD(numberOfStars)
          integer :: numberOfWDs,
     &               disk_belonging(numberOfStars)

          common /tm/ starBirthTime,
     &                m
          common /coorcil/ coordinate_R,
     &                     coordinate_Theta,
     &                     coordinate_Zcylindr
          common /patron/ heightPattern
          common /index/ flagOfWD,numberOfWDs,disk_belonging


          
      end subroutine
