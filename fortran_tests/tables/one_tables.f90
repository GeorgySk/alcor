subroutine incoolone(lgtabone, ltabone, mvtabone, lgtetabone, &
                     mtabone, ndatsone, &
                     bvtabone, vitabone, vrtabone, uvtabone, jone)
    implicit none

    integer, parameter :: NROW = 300, &
                          NCOL = 6, &
                          NCOL2 = 5
    real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13 

    integer :: i, &
               j, &
               file_unit, &
               ndatsone(NCOL)
    real :: mtabone(NCOL), &
            lgtabone(NCOL, NROW), &
            ltabone(NCOL, NROW), &
            mvtabone(NCOL, NROW), &
            lgtetabone(NCOL, NROW), &
            bvtabone(NCOL, NROW), &
            vitabone(NCOL, NROW), &
            vrtabone(NCOL, NROW), &
            uvtabone(NCOL, NROW), &
            jone(NCOL, NROW)
 
    mtabONE(1) = 1.06
    mtabONE(2) = 1.10
    mtabONE(3) = 1.16
    mtabONE(4) = 1.20
    mtabONE(5) = 1.24
    mtabONE(6) = 1.28

    file_unit = 121

    do i = 1, NCOL
        do j = 1, NROW
            read(file_unit, *, end=2) a1, a2, a3, a4, a5, a6, a7, a8, a9, &
                                    a10, a11, a12, a13
            ! luminosity
            ltabone(i, j) = a1
            ! apparent V from UBVRI
            mvtabone(i, j) = a12
            ! log cooling time
            lgtabone(i, j) = a13
            ! log effective temperature
            lgtetabone(i, j) = a2
            ! B - V from UBVRI
            bvtabone(i, j) = a3
            ! V - I from UBVRI
            vitabone(i, j) = a9
            ! V - R from UBVRI
            vrtabone(i, j) = a4
            ! U - V from UBVRI
            uvtabone(i, j) = a10
            jone(i, j) = a7 + a8 - a5 + a12
        end do      
        ! rows count
2       ndatsONE(i) = j - 1  
        file_unit = file_unit + 1
    end do
end subroutine
