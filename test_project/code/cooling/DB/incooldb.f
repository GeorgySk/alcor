      subroutine incooldb(table)
C         Reads cooling tables by:
C             Althaus for Z = 0.001
C             Althaus PG-DB, Z = 0.01
C             Althaus for Z = 0.06
          use external_types
          implicit none
    
          integer :: i, j
          real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, 
     &            a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23,
     &            a24, a25, a26
          TYPE(FileGroupInfo) :: table
    
          if (table%flag == 1) then
C             Mass (Z = 0.001) new by Leandro
              table%mass(1) = 0.5047
              table%mass(2) = 0.5527
              table%mass(3) = 0.59328
              table%mass(4) = 0.62738
              table%mass(5) = 0.6602
              table%mass(6) = 0.69289
              table%mass(7) = 0.8637
          else if (table%flag == 2) then
C             Mass (Z = 0.01) PG-DB
              table%mass(1) = 0.514
              table%mass(2) = 0.53
              table%mass(3) = 0.542
              table%mass(4) = 0.565
              table%mass(5) = 0.584
              table%mass(6) = 0.61
              table%mass(7) = 0.664
              table%mass(8) = 0.741
              table%mass(9) = 0.869
          else
C             Mass (Z = 0.06) LPCODE
              table%mass(1) = 0.524
              table%mass(2) = 0.570
              table%mass(3) = 0.593
              table%mass(4) = 0.61
              table%mass(5) = 0.632
              table%mass(6) = 0.659
              table%mass(7) = 0.70
              table%mass(8) = 0.76
              table%mass(9) = 0.87
          end if
    
          if (table%flag == 1) then
C             Tprew_WD (Z = 0.001) new
              table%prevTime(1) = 0.0
              table%prevTime(2) = 0.0
              table%prevTime(3) = 0.0
              table%prevTime(4) = 0.0
              table%prevTime(5) = 0.0
              table%prevTime(6) = 0.0
              table%prevTime(7) = 0.0     
          else
C             Tprew_WD (Z = 0.06) LPCODE
              table%prevTime(1) = 11.117
              table%prevTime(2) = 2.7004
              table%prevTime(3) = 1.699
              table%prevTime(4) = 1.2114
              table%prevTime(5) = 0.9892
              table%prevTime(6) = 0.7422
              table%prevTime(7) = 0.4431
              table%prevTime(8) = 0.287
              table%prevTime(9) = 0.114
          end if
     
          if (table%flag == 1 .or. table%flag == 3) then
C             Reading the tables
              do i = 1, table%ncol
C                 Reading unit
                  table%initLink = table%initLink + 1
C                 Reading the cooling curves
                  do j = 1, table%nrow
                      read(table%initLink, *, end=2) a1, a2, a3, a4, a5,
     &                                               a6, a7, a8, a9, 
     &                                               a10, a11, a12, a13,
     &                                               a14, a15, a16, a17,
     &                                               a18, a19, a20, a21,
     &                                               a22, a23, a24, a25,
     &                                               a26
                      table%coolingTime(i, j) = 10.0 ** a9 / 1000.0
C                     Restricting the lifetime of pre-WD's
                      table%coolingTime(i, j) = table%coolingTime(i, j)
     &                                          - table%prevTime(i)
                      table%effTemp(i, j) = 10.0 ** a2
                      table%gravAcc(i, j) = a23
                      table%luminosity(i, j) = a1      
                  end do
2                 table%ntrk(i) = j - 1       
              end do
          else
C             Reading the tables
              do i = 1, table%ncol
C                 Reading unit
                  table%initLink = table%initLink + 1 
C                 Reading the cooling curves
                  do j = 1, table%nrow
                      read (table%initLink, *, end=20) a1, a2, a3, a4      
                      table%coolingTime(i, j) = a3
                      table%effTemp(i, j) = 10.0 ** a2
                      table%gravAcc(i, j) = a4
                      table%luminosity(i, j) = a1
                  end do
20                table%ntrk(i) = j - 1
              end do
          end if
      end subroutine
