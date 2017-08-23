      subroutine incoolda(table)
C         Reads cooling tables for sequences of WDs of DA type by 
C         Althaus et al. (2009) and Renedo et al. (2010) and
C         interpolates/extrapolates according to the vector of time 
C         reference. This subroutine can be used to read the tables for 
C         one metallicity in particular, according to the input index.
          use external_types
          implicit none
    
          integer :: j, k
          real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, 
     &            a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23,
     &            a24, a25, a26
          TYPE(FileGroupInfo) :: table
          
C         Reading masses, depending on the format of the input files,
C         because they are different
          if (table%flag == 1) then
C             Masses  (Z = 0.001)
              table%mass(1) = 0.505
              table%mass(2) = 0.553
              table%mass(3) = 0.593
              table%mass(4) = 0.627
              table%mass(5) = 0.660
              table%mass(6) = 0.692
              table%mass(7) = 0.863
C             Tprew_WD LPCODE
              table%prevTime(1) = 0.0
              table%prevTime(2) = 0.0
              table%prevTime(3) = 0.0
              table%prevTime(4) = 0.0
              table%prevTime(5) = 0.0
              table%prevTime(6) = 0.0
              table%prevTime(7) = 0.0
          end if
    
          if (table%flag == 2) then
C             Masses (Z = 0.01)
              table%mass(1) = 0.524
              table%mass(2) = 0.570
              table%mass(3) = 0.593
              table%mass(4) = 0.609
              table%mass(5) = 0.632
              table%mass(6) = 0.659
              table%mass(7) = 0.705
              table%mass(8) = 0.767
              table%mass(9) = 0.837
              table%mass(10) = 0.877
C             Tprew_WD LPCODE
              table%prevTime(1) = 0.0
              table%prevTime(2) = 0.0
              table%prevTime(3) = 0.0
              table%prevTime(4) = 0.0
              table%prevTime(5) = 0.0
              table%prevTime(6) = 0.0
              table%prevTime(7) = 0.0
              table%prevTime(8) = 0.0
              table%prevTime(9) = 0.0
              table%prevTime(10) = 0.0
          end if
    
          if (table%flag == 3) then
C             Masses (Z = 0.03/0.06)
              table%mass(1) = 0.524
              table%mass(2) = 0.570
              table%mass(3) = 0.593
              table%mass(4) = 0.609
              table%mass(5) = 0.632
              table%mass(6) = 0.659
              table%mass(7) = 0.705
              table%mass(8) = 1.000
C             Tprew_WD LPCODE
              table%prevTime(1) = 11.117
              table%prevTime(2) = 2.7004
              table%prevTime(3) = 1.699
              table%prevTime(4) = 1.2114
              table%prevTime(5) = 0.9892
              table%prevTime(6) = 0.7422
              table%prevTime(7) = 0.4431
              table%prevTime(8) = 0.0
          end if
    
C         Choosing the format of the input file
          if (table%flag == 1 .or. table%flag == 2) then
C             Reading the tables
              do k = 1, table%ncol
C                 Reading the cooling curves
                  do j = 1, table%nrow
                      read(table%initLink + k, *, end=15) a1, a2, a3, 
     &                                                    a4, a5, a6, 
     &                                                    a7, a8, a9,
     &                                                    a10, a11, a12,
     &                                                    a13
                      table%coolingTime(k, j) = a5 / 1000.0
                      table%effTemp(k, j) = 10 ** a2
                      table%gravAcc(k, j) = a12
                      table%luminosity(k, j) = a1
                  end do
15                table%ntrk(k) = j - 1
              end do
          else
C             Reading the tables
              do k = 1, table%ncol
C                 Reading the cooling curves
                  do j = 1, table%nrow
                      read(table%initLink + k, *, end=20) a1, a2, a3, 
     &                                                    a4, a5, a6, 
     &                                                    a7, a8, a9,
     &                                                    a10, a11, a12,
     &                                                    a13, a14, a15,
     &                                                    a16, a17, a18,
     &                                                    a19, a20, a21,
     &                                                    a22, a23, a24,
     &                                                    a25, a26
                      table%coolingTime(k, j) = 10.0 ** a9 / 1000.0
C                     Restricting the lifetime of pre-WD
                      table%coolingTime(k, j) = table%coolingTime(k, j)
     &                                          - table%prevTime(k)
                      table%effTemp(k, j) = 10 ** a2
                      table%gravAcc(k, j) = a23
                      table%luminosity(k, j) = a1   
                  end do
20                table%ntrk(k) = j - 1
              end do
          end if
      end subroutine
