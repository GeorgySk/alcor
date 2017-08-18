C     TODO: give a better name
      subroutine dbd_fid(iseed,
     &                   DB_fraction,
     &                   in)
C     Determines if WD is DA or non-DA (DB) based on DB_fraction
C     TODO: change to logical type 
C     in = 0 DA
C     in = 1 non-DA
          implicit none

          external ran
          real ran
          integer :: iseed,
     &               in
          real :: DB_fraction
      
          if (ran(iseed) < DB_fraction) then 
              in = 1
          else
              in = 0
          end if
      end subroutine
