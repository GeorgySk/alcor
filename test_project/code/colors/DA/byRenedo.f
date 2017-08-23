      subroutine color(table)
C         Reading colors by Rene and interpolating according to the 
C         vector of reference time
          use external_types
          implicit none

          real, parameter :: VTANMIN = 30
          TYPE(FileGroupInfo) :: table
          integer :: i, k
          real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, 
     &            a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23,
     &            a24, a25, a26, a27, a28

C         Read masses
          table%mass(1) = 0.524
          table%mass(2) = 0.570
          table%mass(3) = 0.593
          table%mass(4) = 0.610
          table%mass(5) = 0.632
          table%mass(6) = 0.659
          table%mass(7) = 0.705
          table%mass(8) = 0.767
          table%mass(9) = 0.837
          table%mass(10) = 0.877
          
C         Read values from files
          do k = 1, table%ncol
              do i = 1, table%nrow
                  read(table%initLink + k, *, end=15) a1, a2, a3, a4, 
     &                                                a5, a6, a7, a8, 
     &                                                a9, a10, a11, a12,
     &                                                a13, a14, a15, 
     &                                                a16, a17, a18, 
     &                                                a19, a20, a21, 
     &                                                a22, a23, a24, 
     &                                                a25, a26, a27, a28
                  table%luminosity(k,i) = a3
                  table%color_U(k,i) = a24
                  table%color_B(k,i) = a25
                  table%color_V(k,i) = a26
                  table%color_R(k,i) = a27
                  table%color_I(k,i) = a28
              end do
15            table%ntrk(k) = i - 1
          end do
      end subroutine
