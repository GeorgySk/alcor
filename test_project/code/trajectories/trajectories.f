      subroutine traject(thin_disk_age)
C     Calculating trajectories of WDs according to z-axis using the 4th 
C     order Runge-Kuttta.
      implicit none
      
C     TODO: take this up (ex numberOfStars)
C     TODO: find out the meaning of njumps, n
      integer, parameter :: MAX_STARS_COUNT = 6000000,
     &                      NJUMPS = 100,
     &                      N = 2
C     TODO: find out the meaning of hmin, wosun
      real, parameter :: WOSUN = -8.0,
     &                   HMIN = 0.0,
     &                   EPSILON = 1.0e-4,
     &                   SECONDS_IN_HOUR = 3600.0,
     &                   HOURS_IN_DAY = 24.0,
     &                   DAYS_IN_YEAR = 365.25,
     &                   YEARS_IN_GYR = (1.0e+9),
     &                   SECONDS_IN_GYR = SECONDS_IN_HOUR 
     &                                    * HOURS_IN_DAY
     &                                    * DAYS_IN_YEAR
     &                                    * YEARS_IN_GYR,
     &                   METERS_IN_PARSEC = 3.086e+16
      integer :: numberOfWDs,
     &           wd_index,
     &           NOK,
     &           NBAD,
     &           flagOfWD(MAX_STARS_COUNT)
      real :: thin_disk_age,
     &        final_time,
     &        xcar,
     &        ycar,
     &        wo,
     &        zo,
     &        ecini,
     &        time_increment,
     &        htry,
     &        initial_time,
     &        ecinf,
     &        f
      real :: uu(MAX_STARS_COUNT),
     &        vv(MAX_STARS_COUNT),
     &        ww(MAX_STARS_COUNT),
     &        starBirthTime(MAX_STARS_COUNT),
     &        m(MAX_STARS_COUNT),
     &        yscal(2),
     &        y(2),
     &        dydx(2),
     &        xpla(MAX_STARS_COUNT),
     &        ypla(MAX_STARS_COUNT)
      double precision :: coordinate_R(MAX_STARS_COUNT),
     &                    coordinate_Theta(MAX_STARS_COUNT),
     &                    coordinate_Zcylindr(MAX_STARS_COUNT)
      integer :: disk_belonging(MAX_STARS_COUNT)
        
      common /vel/ uu, vv, ww
      common /coorcil/ coordinate_R,
     &                 coordinate_Theta,
     &                 coordinate_Zcylindr
      common /tm/ starBirthTime, m
      common /index/ flagOfWD,
     &               numberOfWDs,
     &               disk_belonging         
C     TODO: rename these as in other subroutines and find out the mean..    
      common /plano/ xpla, ypla
      common /carte/ xcar, ycar

C     External - specifies procedures as external, and allows their 
C     symbolic names to be used as actual arguments.
      external DERIVS
      external RKQC
            
      final_time = thin_disk_age * SECONDS_IN_GYR 

C     Integrating trajectories
      do wd_index = 1, numberOfWDs
          xcar = xpla(wd_index)
          ycar = ypla(wd_index)
C         TODO: find out the meaning of wo and 8.0
          wo = ww(wd_index) + 8.0
          zo = real(coordinate_Zcylindr(wd_index) * METERS_IN_PARSEC)
C         TODO: find out the meaning of ecini
          ecini = 0.5 * wo * wo
          time_increment = (thin_disk_age - starBirthTime(wd_index)) 
     &                     / float(NJUMPS)
C         Time in seconds
C         TODO: find out the meaning of htry
          htry = time_increment * SECONDS_IN_GYR
C         Initial conditions
C         TODO: find out the meaning of y and dydx
          y(1) = zo
          y(2) = wo
          dydx(1) = wo
          call fuerza(zo, f)
          dydx(2) = f
          initial_time = starBirthTime(wd_index) * SECONDS_IN_GYR 
C         Calling to the Runge-Kutta integrator
          call ODEINT(y, N, initial_time, final_time, EPSILON, htry,
     &                HMIN, NOK, NBAD, DERIVS, RKQC, yscal, y, dydx)
C         TODO: find out the meaning of ecinf    
          ecinf = 0.5 * y(2) * y(2)
          coordinate_Zcylindr(wd_index) = y(1) / METERS_IN_PARSEC
          ww(wd_index) = y(2) + WOSUN
      end do
      end subroutine


      subroutine fuerza(z_coordinate, force)
C     Calculating the force along z-coordinate. z in km
      implicit none

C     TODO: find out the meaning of all the variables
      real, parameter :: ro = 8.5,
     &                   squared_ro = ro * ro,
     &                   vh = 220.0,
     &                   squared_vh = vh * vh,
     &                   rc1 = 2.7,
     &                   squared_rc1 = rc1 * rc1,
     &                   mc1 = 3.0d+09,
     &                   rc2 = 0.42,
     &                   squared_rc2 = rc2 * rc2,
     &                   mc2 = 1.6d+10,
     &                   b = 0.3,
     &                   squared_b = b * b,
     &                   md1 = 6.6e+10, 
     &                   a1 = 5.81, 
     &                   md2 = -2.9e+10,
     &                   a2 = 17.43,
     &                   md3 = 3.3d+09, 
     &                   a3 = 34.86,
     &                   g = 4.30026e-6,
C                        TODO: take this out
     &                   METERS_IN_PARSEC = 3.086e+16

      real :: xcar, 
     &        ycar, 
     &        xpla, 
     &        ypla,
     &        rpla,
     &        zsig, 
     &        z_coordinate,
     &        rpla2,
     &        r2,
     &        halo_force, 
     &        central_force_1,
     &        central_force_2,
     &        central_force,
     &        bzr,
     &        disk_force_1,
     &        disk_force_2,
     &        disk_force_3,
     &        disk_force,
     &        total_force,
     &        fcv,
     &        force

      common /carte/ xcar, ycar

      xpla = xcar
      ypla = ycar
      zsig = z_coordinate
      z_coordinate = abs(z_coordinate) / METERS_IN_PARSEC
      rpla2 = xpla * xpla + ypla * ypla
      rpla = sqrt(rpla2)
      r2 = rpla * rpla + z_coordinate * z_coordinate
                   
C     Calculating the forces
C     Dark halo
      halo_force = squared_vh * z_coordinate / (squared_ro + r2)

C     Central component      
      central_force_1 = g * mc1 * z_coordinate 
     &                  / ((squared_rc1 + r2) ** 1.5)   
      central_force_2 = g * mc2 * z_coordinate
     &                  / ((squared_rc2 + r2) ** 1.5)      
      central_force = central_force_1 + central_force_2

C     Disk
      bzr = sqrt(squared_b + z_coordinate * z_coordinate)
      disk_force_1 = g * md1 * z_coordinate * (a1 + bzr) 
     &       / (bzr * (rpla2 + (a1 + bzr) * (a1 + bzr)) ** 1.5)    
      disk_force_2 = g * md2 * z_coordinate * (a2 + bzr) 
     &       / (bzr * (rpla2 + (a2 + bzr) * (a2 + bzr)) ** 1.5)     
      disk_force_3 = g * md3 * z_coordinate * (a3 + bzr) 
     &       /(bzr * (rpla2 + (a3 + bzr) * (a3 + bzr)) ** 1.5)    
      disk_force = disk_force_1 + disk_force_2 + disk_force_3

C     Total force
      total_force = halo_force + central_force + disk_force 

C     If we want the result in km/s²
C     TODO: find out the meaning of the following const, conversion?
      fcv = 1.0 / METERS_IN_PARSEC
      total_force = fcv * abs(total_force)

C     The sign of z_coordinate will be
      force = -sign(total_force, zsig)
      z_coordinate = zsig
      
      end subroutine
