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

      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

      real :: luminosityOfWD(numberOfStars),
     &        massOfWD(numberOfStars),
     &        metallicityOfWD(numberOfStars),
     &        effTempOfWD(numberOfStars)
      integer :: flagOfWD(numberOfStars)
      real :: coolingTime(numberOfStars)
      integer typeOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)

      real :: u_ubvrij(numberOfStars),
     &        b_ubvrij(numberOfStars),
     &        v_ubvrij(numberOfStars),
     &        r_ubvrij(numberOfStars),
     &        i_ubvrij(numberOfStars),
     &        j_ubvrij(numberOfStars)

      real :: rgac(numberOfStars)
      double precision :: lgac(numberOfStars),
     &                    bgac(numberOfStars)

      TYPE(FileGroupInfo),DIMENSION(11) :: table

      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD
      common /index/ flagOfWD,numberOfWDs,disk_belonging
      common /paral/ rgac
      common /cool/ coolingTime
      common /indexdb/ typeOfWD
      common /lb/ lgac,bgac

      n1=0
      n2=0
      n3=0
      n4=0
      n5=0

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

          u_ubvrij(i) = c1
          b_ubvrij(i) = c2
          v_ubvrij(i) = c3
          r_ubvrij(i) = c4
          i_ubvrij(i) = c5
          j_ubvrij(i) = c6

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
