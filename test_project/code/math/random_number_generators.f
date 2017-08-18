C     NOTE: this all function is an old bunch of.. don't touch it
      function gasdev(iseed)
C         Returns a normally distributed deviate with zero mean and unit
C         variance
          implicit none

          external ran
          real ran

          save

          integer :: iseed,
     &               iset
          real :: v1,
     &            v2,
     &            r,
     &            fac,
     &            gset,
     &            gasdev
          
          data iset/0/
          if (iset == 0) then
  1           v1 = 2.0 * (ran(iseed)) - 1.0
              v2 = 2.0 * (ran(iseed)) - 1.0
              r = v1 * v1 + v2 * v2
              if (r >= 1.0 .or. r == 0.0) then 
                  goto 1
              end if
              fac = sqrt(-2.0e0 * log(r) / r)
              gset = v1 * fac
              gasdev = v2 * fac
              iset = 1
          else
              gasdev = gset
              iset = 0
          end if
      end function
