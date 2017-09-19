subroutine color(initlink, ncol, nrow, ntrk, mass, luminosity, color_u, &
                 color_b, color_v, color_r, color_i)
    implicit none

    integer, intent(in) :: initLink
    integer, intent(in) :: ncol
    integer, intent(in) :: nrow
    integer, intent(inout) :: ntrk(ncol)
    real, intent(inout) :: mass(ncol)
    real, intent(inout) :: luminosity(ncol, nrow)
    real, intent(inout) :: color_U(ncol, nrow)
    real, intent(inout) :: color_B(ncol, nrow)
    real, intent(inout) :: color_V(ncol, nrow)
    real, intent(inout) :: color_R(ncol, nrow)
    real, intent(inout) :: color_I(ncol, nrow)
    integer :: i, k
    real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, &
            a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, &
            a24, a25, a26, a27, a28

    mass(1) = 0.524
    mass(2) = 0.570
    mass(3) = 0.593
    mass(4) = 0.610
    mass(5) = 0.632
    mass(6) = 0.659
    mass(7) = 0.705
    mass(8) = 0.767
    mass(9) = 0.837
    mass(10) = 0.877

    do k = 1, ncol
        do i = 1, nrow
            read(initLink + k, *, end=15) a1, a2, a3, a4, &
                                          a5, a6, a7, a8,  &
                                          a9, a10, a11, a12, &
                                          a13, a14, a15, &
                                          a16, a17, a18, &
                                          a19, a20, a21, &
                                          a22, a23, a24,  &
                                          a25, a26, a27, a28
                  luminosity(k,i) = a3
                  color_U(k,i) = a24
                  color_B(k,i) = a25
                  color_V(k,i) = a26
                  color_R(k,i) = a27
                  color_I(k,i) = a28
        end do
15      ntrk(k) = i - 1
    end do
end subroutine
