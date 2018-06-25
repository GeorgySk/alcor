      subroutine magi(fractionOfDB,
     &                table,
     &                u_ugriz,
     &                g_ugriz,
     &                r_ugriz,
     &                i_ugriz,
     &                z_ugriz)
      use external_types
C     This subroutine calculates ltc,cbv,cvi,cvr,cuv visual absolute 
C     and apparent magnitude of the WDs.
      implicit real (a-h,m,o-z)

      integer numberOfStars,iseed,i,numberOfWDs,in
      real lum,teff,xlog,c1,c2,c3,c4,c5,c6,n1,n2,n3,n4,n5
      real fractionOfDB
      real mone

      real :: AVT,
     &        SAVT,
     &        AVC,
     &        AV(5),
     &        SAV(5),
     &        pi
      integer :: JMAX

      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

      real :: luminosityOfWD(numberOfStars),
     &        massOfWD(numberOfStars),
     &        metallicityOfWD(numberOfStars),
     &        effTempOfWD(numberOfStars),
     &        log_g(numberOfStars)
      integer :: flagOfWD(numberOfStars)
      real :: coolingTime(numberOfStars)
      integer typeOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)
      logical :: eliminated(numberOfStars)

      real :: u_ubvrij,
     &        b_ubvrij,
     &        v_ubvrij,
     &        r_ubvrij,
     &        i_ubvrij,
     &        j_ubvrij,
     &        extinction

      real :: ug, 
     &        gr,
     &        ri,
     &        iz

      real :: u_ugriz(numberOfStars), 
     &        g_ugriz(numberOfStars),
     &        r_ugriz(numberOfStars),
     &        i_ugriz(numberOfStars),
     &        z_ugriz(numberOfStars)

      real ::
     & u_ugriz_w_error, g_ugriz_w_error, r_ugriz_w_error, 
     & i_ugriz_w_error, z_ugriz_w_error,
     & ug_up, ug_low, gr_up, gr_low, top, left, right

      real :: rgac(numberOfStars)
      double precision :: lgac(numberOfStars),
     &                    bgac(numberOfStars)

      real :: d_ug, d_gr, d_ri, d_iz

      TYPE(FileGroupInfo),DIMENSION(11) :: table

      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD, log_g
      common /index/ flagOfWD,numberOfWDs,disk_belonging
      common /paral/ rgac
      common /cool/ coolingTime
      common /indexdb/ typeOfWD, eliminated
      common /lb/ lgac,bgac

      n1=0
      n2=0
      n3=0
      n4=0
      n5=0

      pi = 4.0 * atan(1.0)

C     ---  Interpolating Mv, luminosity, colors and other variables
C          from coolingTime and the mwd  ---
      do i = 1, numberOfWDs
C       ---  ATENTION! choosing only if .lt.1.1!!!  ---
        if (massOfWD(i) .le. 1.4) then
C         ---  IF CO core ---
          if (massOfWD(i) .lt. mone) then  
C           --- Atention We put the "old" ones to 0.6Msol ---
C           --- Distribucion DA/DB ---
            call dbd_fid(iseed,fractionOfDB,in)
C           --- End of distribution ---
C           ---  IF DA ---
            if(in.eq.0) then
              typeOfWD(i)=0
              n1=n1+1
              call interlumda(coolingTime(i),massOfWD(i),
     &             metallicityOfWD(i),lum,teff,xlog,c1,c2,c3,c4,c5,c6,
     &             table)
C           ---  ELSE DB  ---
            else
              n3=n3+1
              typeOfWD(i)=1    
              call interlumdb(coolingTime(i),massOfWD(i),
     &             metallicityOfWD(i),lum,c1,c2,c3,c4,c5,c6,teff,xlog,
     &             table)
              if(teff.lt.6000) n5=n5+1
            end if
C           ---  END IF DB/NON-DB
C         ---  ELSE ONe ---
          else
            n2=n2+1
            typeOfWD(i)=2
            call interlumone(coolingTime(i),massOfWD(i),lum,c1,c2,c3,c4,
     &           c5,c6,teff,xlog)
          end if
C         ---  END IF CO/ONe ---
          if (teff < 6000) then
            n4 = n4 + 1
          end if

          luminosityOfWD(i) = -lum
          effTempOfWD(i) = teff
          log_g(i) = xlog

C         Calculating extinction
          call EXTINCT(real(lgac(i) * 180.0 / pi),
     &                 real(bgac(i) * 180.0 / pi),
     &                 real(rgac(i)),
     &                 AVT,SAVT,AVC,JMAX,AV,SAV)
          extinction = AVT + AVC

          u_ubvrij = c1 + extinction * 1.664
          b_ubvrij = c2 + extinction * 1.321
          v_ubvrij = c3 + extinction * 1.015
          r_ubvrij = c4 + extinction * 0.819
          i_ubvrij = c5 + extinction * 0.594
          j_ubvrij = c6 + extinction * 0.276

          g_ugriz(i) = v_ubvrij + 0.63 * (b_ubvrij - v_ubvrij)
     &              - 0.124
          u_ugriz(i) = (g_ugriz(i)
     &               + 0.75 * (u_ubvrij - b_ubvrij)
     &               + 0.77 * (b_ubvrij - v_ubvrij)
     &               + 0.72)
          r_ugriz(i) = (g_ugriz(i) 
     &               - 1.646 * (v_ubvrij - r_ubvrij) + 0.139)
          i_ugriz(i) = (r_ugriz(i)
     &               - 1.007 * (r_ubvrij - i_ubvrij) + 0.236)
          z_ugriz(i) = (i_ugriz(i) 
     &               - (1.584 - 1.007) * (r_ubvrij - i_ubvrij)
     &               + 0.386 - 0.236)

          call errfot(u_ugriz(i),
     &                u_ugriz_w_error,
     &                1)
          call errfot(g_ugriz(i),
     &                g_ugriz_w_error,
     &                2)
          call errfot(r_ugriz(i),
     &                r_ugriz_w_error,
     &                3)
          call errfot(i_ugriz(i),
     &                i_ugriz_w_error,
     &                4)
          call errfot(z_ugriz(i),
     &                z_ugriz_w_error,
     &                5)

          u_ugriz(i) = u_ugriz_w_error
          g_ugriz(i) = g_ugriz_w_error
          r_ugriz(i) = r_ugriz_w_error
          i_ugriz(i) = i_ugriz_w_error
          z_ugriz(i) = z_ugriz_w_error

          ug = u_ugriz(i) - g_ugriz(i)
          gr = g_ugriz(i) - r_ugriz(i)
          ri = r_ugriz(i) - i_ugriz(i)
          iz = i_ugriz(i) - z_ugriz(i)

          ug_up = -24.384 * gr ** 5 - 19.0 * gr ** 4 
     &            + 3.497 * gr ** 3 + 1.193 * gr ** 2 
     &            + 0.083 * gr + 0.61
          ug_low = -20.653 * gr ** 5 + 10.816 * gr ** 4 
     &            + 15.718 * gr ** 3 - 1.294 * gr ** 2 
     &            - 0.084 * gr + 0.3
          gr_up = -0.6993 * ri ** 2 + 0.947 * ri + 0.192
          gr_low = -1.32 * ri ** 3 + 2.173 * ri ** 2 
     &             + 2.452 * ri - 0.07
          top = -0.56
          left = 0.176 * iz + 0.127
          right = -0.754 * iz + 0.11

          if ((ug > ug_low) .and. (ug < ug_up) 
     &        .and. (gr > gr_low) .and. (gr < gr_up) 
     &        .and. (ri > top) 
     &        .and. (ri < left) .and. (ri < right)
     &       ) then
            g_ugriz_apparent = g_ugriz(i) - 5.0 + 5.0 * (log10(rgac(i)) 
     &                                                   + 3.0)
C             i_ugriz_apparent = i_ugriz - 5.0 + 5.0 * (log10(rgac(i)) 
C      &                                                + 3.0)
C             if (i_ugriz_apparent < 15.0 
C      &          .or. i_ugriz_apparent > 19.1) then
            if (g_ugriz_apparent < 15.0 
     &          .or. g_ugriz_apparent > 22.0) then
              eliminated(i) = .true.
            else
              eliminated(i) = .false.
            end if
          else
              eliminated(i) = .true.
          end if

C       ---  ELSE mass >= 1.4  --- EXPLOTA, exceeding Chandrasekar limit
        else
          typeOfWD(i) = 5
        end if
      end do

      write(*,*) "DA CO ",n1
      write(*,*) "DA ONe ",n2
      write(*,*) "DB",n3
      write(*,*) "<6000 DA", n4
      write(*,*) "<6000 DB", n5

      return
      end
