      subroutine magi(fractionOfDB,
     &                table,
     &                u_ubvrij,
     &                b_ubvrij,
     &                v_ubvrij,
     &                r_ubvrij,
     &                i_ubvrij,
     &                j_ubvrij)
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

      real :: u_ubvrij(numberOfStars),
     &        b_ubvrij(numberOfStars),
     &        v_ubvrij(numberOfStars),
     &        r_ubvrij(numberOfStars),
     &        i_ubvrij(numberOfStars),
     &        j_ubvrij(numberOfStars),
     &        extinction

      real :: ug, 
     &        gr,
     &        ri,
     &        iz

      real :: u_ugriz, 
     &        g_ugriz,
     &        r_ugriz,
     &        i_ugriz,
     &        z_ugriz

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

          u_ubvrij(i) = c1 + extinction * 1.664
          b_ubvrij(i) = c2 + extinction * 1.321
          v_ubvrij(i) = c3 + extinction * 1.015
          r_ubvrij(i) = c4 + extinction * 0.819
          i_ubvrij(i) = c5 + extinction * 0.594
          j_ubvrij(i) = c6 + extinction * 0.276

C           g_ugriz = v_ubvrij(i) + 0.63 * (b_ubvrij(i) - v_ubvrij(i))
C      &              - 0.124
C           u_ugriz = (g_ugriz
C      &               + 0.75 * (u_ubvrij(i) - b_ubvrij(i))
C      &               + 0.77 * (b_ubvrij(i) - v_ubvrij(i))
C      &               + 0.72)
C           r_ugriz = (g_ugriz 
C      &               - 1.646 * (v_ubvrij(i) - r_ubvrij(i)) + 0.139)
C           i_ugriz = (r_ugriz
C      &               - 1.007 * (r_ubvrij(i) - i_ubvrij(i)) + 0.236)
C           z_ugriz = (i_ugriz 
C      &               - (1.584 - 1.007) * (r_ubvrij(i) - i_ubvrij(i))
C      &               + 0.386 - 0.236)

C           call errfot(u_ugriz,
C      &                u_ugriz_w_error,
C      &                1)
C           call errfot(g_ugriz,
C      &                g_ugriz_w_error,
C      &                2)
C           call errfot(r_ugriz,
C      &                r_ugriz_w_error,
C      &                3)
C           call errfot(i_ugriz,
C      &                i_ugriz_w_error,
C      &                4)
C           call errfot(z_ugriz,
C      &                z_ugriz_w_error,
C      &                5)

C           u_ugriz = u_ugriz_w_error
C           g_ugriz = g_ugriz_w_error
C           r_ugriz = r_ugriz_w_error
C           i_ugriz = i_ugriz_w_error
C           z_ugriz = z_ugriz_w_error

C           ug = u_ugriz - g_ugriz
C           gr = g_ugriz - r_ugriz
C           ri = r_ugriz - i_ugriz
C           iz = i_ugriz - z_ugriz

C           ug_up = -24.384 * gr ** 5 - 19.0 * gr ** 4 
C      &            + 3.497 * gr ** 3 + 1.193 * gr ** 2 
C      &            + 0.083 * gr + 0.61
C           ug_low = -20.653 * gr ** 5 + 10.816 * gr ** 4 
C      &            + 15.718 * gr ** 3 - 1.294 * gr ** 2 
C      &            - 0.084 * gr + 0.3
C           gr_up = -0.6993 * ri ** 2 + 0.947 * ri + 0.192
C           gr_low = -1.32 * ri ** 3 + 2.173 * ri ** 2 
C      &             + 2.452 * ri - 0.07
C           top = -0.56
C           left = 0.176 * iz + 0.127
C           right = -0.754 * iz + 0.11

C           if ((ug > ug_low) .and. (ug < ug_up) 
C      &        .and. (gr > gr_low) .and. (gr < gr_up) 
C      &        .and. (ri > top) 
C      &        .and. (ri < left) .and. (ri < right)
C      &       ) then
C             g_ugriz_apparent = g_ugriz(i) - 5.0 + 5.0 * (log10(rgac(i)) 
C      &                                                   + 3.0)
C C             i_ugriz_apparent = i_ugriz - 5.0 + 5.0 * (log10(rgac(i)) 
C C      &                                                + 3.0)
C C             if (i_ugriz_apparent < 15.0 
C C      &          .or. i_ugriz_apparent > 19.1) then
C             if (g_ugriz_apparent < 15.0 
C      &          .or. g_ugriz_apparent > 22.0) then
C               eliminated(i) = .true.
C             else
C               eliminated(i) = .false.
C             end if
C           else
C               eliminated(i) = .true.
C           end if

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
