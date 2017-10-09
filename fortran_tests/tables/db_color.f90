subroutine colordb(initlink, ncol, nrow, ntrk, mass, luminosity, color_u, &
                   color_b, color_v, color_r, color_i, color_j)
    implicit none
    
    integer :: i, j
    real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, &
            a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, &
            a24, a25, a26, a27
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
    real, intent(inout) :: color_J(ncol, nrow)
    
    mass(1) = 0.5
    mass(2) = 0.6
    mass(3) = 0.7
    mass(4) = 0.8
    mass(5) = 0.9
    mass(6) = 1.0
    mass(7) = 1.2
     
    do i = 1, ncol
        do j = 1, nrow
            read(initLink + i, *, end=2) a1, a2, a3, a4, a5, a6, &
                                         a7, a8, a9, a10, a11, &
                                         a12, a13, a14, a15, &
                                         a16, a17, a18, a19, &
                                         a20, a21, a22, a23, &
                                         a24, a25, a26, a27
            luminosity(i, j) = -(a3 - 4.72) / 2.5
            color_U(i, j) = a5
            color_B(i, j) = a6
            color_V(i, j) = a7
            color_R(i, j) = a8
            color_I(i, j) = a9
            color_J(i, j) = a10
        end do
2           ntrk(i) = j - 1    
    end do
end subroutine
