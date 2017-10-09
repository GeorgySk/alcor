subroutine incoolda(flag, initLink, ncol, nrow, &
                    ntrk, mass, coolingTime, prevTime, luminosity, effTemp, &
                    gravAcc)
    implicit none
    
    integer :: j, k
    real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, &
            a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, &
            a24, a25, a26
    integer, intent(in) :: flag
    integer, intent(in) :: initLink
    integer, intent(in) :: ncol
    integer, intent(in) :: nrow
    integer, intent(inout) :: ntrk(ncol)
    real, intent(inout) :: mass(ncol)
    real, intent(inout) :: coolingTime(ncol, nrow)
    real, intent(inout) :: prevTime(ncol)
    real, intent(inout) :: luminosity(ncol, nrow)
    real, intent(inout) :: effTemp(ncol, nrow)
    real, intent(inout) :: gravAcc(ncol, nrow)
          
    if (flag == 1) then
        mass(1) = 0.505
        mass(2) = 0.553
        mass(3) = 0.593
        mass(4) = 0.627
        mass(5) = 0.660
        mass(6) = 0.692
        mass(7) = 0.863
        prevTime(1) = 0.0
        prevTime(2) = 0.0
        prevTime(3) = 0.0
        prevTime(4) = 0.0
        prevTime(5) = 0.0
        prevTime(6) = 0.0
        prevTime(7) = 0.0
    end if
    
    if (flag == 2) then
        mass(1) = 0.524
        mass(2) = 0.570
        mass(3) = 0.593
        mass(4) = 0.609
        mass(5) = 0.632
        mass(6) = 0.659
        mass(7) = 0.705
        mass(8) = 0.767
        mass(9) = 0.837
        mass(10) = 0.877
        prevTime(1) = 0.0
        prevTime(2) = 0.0
        prevTime(3) = 0.0
        prevTime(4) = 0.0
        prevTime(5) = 0.0
        prevTime(6) = 0.0
        prevTime(7) = 0.0
        prevTime(8) = 0.0
        prevTime(9) = 0.0
        prevTime(10) = 0.0
    end if
    
    if (flag == 3) then
        mass(1) = 0.524
        mass(2) = 0.570
        mass(3) = 0.593
        mass(4) = 0.609
        mass(5) = 0.632
        mass(6) = 0.659
        mass(7) = 0.705
        mass(8) = 1.000
        prevTime(1) = 11.117
        prevTime(2) = 2.7004
        prevTime(3) = 1.699
        prevTime(4) = 1.2114
        prevTime(5) = 0.9892
        prevTime(6) = 0.7422
        prevTime(7) = 0.4431
        prevTime(8) = 0.0
    end if
    
    if (flag == 1 .or. flag == 2) then
        do k = 1, ncol
            do j = 1, nrow
                read(initLink + k, *, end=15) a1, a2, a3, &
                                              a4, a5, a6, &
                                              a7, a8, a9, &
                                              a10, a11, a12, &
                                              a13
                coolingTime(k, j) = a5 / 1000.0
                effTemp(k, j) = 10 ** a2
                gravAcc(k, j) = a12
                luminosity(k, j) = a1
            end do
15          ntrk(k) = j - 1
        end do
    else
        do k = 1, ncol
            do j = 1, nrow
                read(initLink + k, *, end=20) a1, a2, a3, &
                                              a4, a5, a6, &
                                              a7, a8, a9, &
                                              a10, a11, a12, &
                                              a13, a14, a15, &
                                              a16, a17, a18, &
                                              a19, a20, a21, &
                                              a22, a23, a24, &
                                              a25, a26
                coolingTime(k, j) = 10.0 ** a9 / 1000.0
                coolingTime(k, j) = coolingTime(k, j) - prevTime(k)
                effTemp(k, j) = 10 ** a2
                gravAcc(k, j) = a23
                luminosity(k, j) = a1   
            end do
20          ntrk(k) = j - 1
        end do
    end if
end subroutine
