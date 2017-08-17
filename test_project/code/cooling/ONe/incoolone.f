      subroutine incoolone
C     Reading the cooling tables of ONe WDs
C     TODO: apply following names to commons
C           mtabone:  mass of the WD of ONE 
C           ltabone:  logarithm of the luminosity
C           mvtabone: visual absolute magnitude visual
C           lgtabone: logarithm of the cooling time in years
      implicit none
      
C     TODO: find out the meaning of following constants
      integer, parameter :: NROW = 300,
     &                      NCOL = 6,
     &                      NCOL2 = 5
      real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13 

      integer :: i,
     &           j,
     &           ird,
     &           ndatsone(NCOL),
     &           ndatsone2(NCOL2)
      real :: mtabone(NCOL),
     &        mtabone2(NCOL2),
     &        lgtabone(NCOL, NROW),
     &        ltabone(NCOL, NROW),
     &        mvtabone(NCOL, NROW),
     &        lgtetabone(NCOL, NROW),
     &        bvtabone(NCOL, NROW),
     &        vitabone(NCOL, NROW),
     &        vrtabone(NCOL, NROW),
     &        uvtabone(NCOL, NROW),
     &        lgrtabone(NCOL2, NROW),
     &        lgt2tabone(NCOL2, NROW)

      common /fredone/ lgtabone, ltabone, mvtabone, lgtetabone
      common /fredone2/ mtabone, ndatsone
      common /colorsone/ bvtabone, vitabone, vrtabone, uvtabone
      common /newone/ lgrtabone, lgt2tabone
      common /newone2/ mtabone2, ndatsone2

C     TODO: give a btter name to ird
C     ird is the unit of reading archives
      ird = 121
 
C     Masses of each archive ---
      mtabONE(1) = 1.06
      mtabONE(2) = 1.10
      mtabONE(3) = 1.16
      mtabONE(4) = 1.20
      mtabONE(5) = 1.24
      mtabONE(6) = 1.28
      mtabone2(1) = 1.06
      mtabone2(2) = 1.10
      mtabone2(3) = 1.16
      mtabone2(4) = 1.20
      mtabone2(5) = 1.28

C     Reading the files
      do i = 1, NCOL
        do j = 1, NROW
C         TODO: find out the way to avoid using end-goto
          read(ird, *, end=2) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13
          ltabone(i, j) = a1
          mvtabone(i, j) = a12
          lgtabone(i, j) = a13
          lgtetabone(i, j) = a2
          bvtabone(i, j) = a3
          vitabone(i, j) = a9
          vrtabone(i, j) = a4
          uvtabone(i, j) = a10
        end do      
2       ndatsONE(i) = j - 1  
        ird = ird + 1
      end do

C     Reading the data of log Teff and log g
      ird = 127      
      do i = 1, NCOL2
        do j = 1, NROW
C         TODO: find out the way to avoid using end-goto
          read(ird, *, end=3) a1, a2, a3, a4, a5, a6
C         Converting radii in cm to radii in solar radius
          a3 = 10.0 ** a3
          a3 = a3 / (6.96e10)
          a3 = log10(a3)
          lgrtabone(i, j) = a3
          lgt2tabone(i, j) = a5 
        end do      
3       continue       
        ndatsONE2(i) = j - 1  
        ird = ird + 1
      end do
      end subroutine
