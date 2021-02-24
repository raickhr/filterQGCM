module operators
    use kinds
    use gridMod
    use constants
    use fields
    implicit none

    contains

    subroutine calc_grad(nxpo, nypo, inField, gradX, gradY, &
                        boundTreatment)
        INTEGER(kind = i4) , INTENT(IN) :: nxpo, nypo
        REAL(kind=r8),DIMENSION(:,:), INTENT(IN) :: inField
        REAL(kind=r8),DIMENSION(:,:), INTENT(OUT) :: gradX, gradY
        CHARACTER(len=*), OPTIONAL :: boundTreatment
        
        REAL(kind=r8) :: dx,dy
        REAL(kind=r8), ALLOCATABLE :: left(:,:), right(:,:), & 
                         up(:,:), down(:,:)

        ALLOCATE(left(nxpo, nypo), right(nxpo, nypo), up(nxpo, nypo), down(nxpo, nypo))

        left = cshift(inField, SHIFT=-1, DIM=1)
        right = cshift(inField, SHIFT=1, DIM=1)

        up = cshift(inField, SHIFT=1, DIM=2)
        down = cshift(inField, SHIFT=-1, DIM=2)

        if (PRESENT(boundTreatment)) then
            if (boundTreatment == 'zero') then 
                left(1,:) = 0 
                right(nxpo,:) = 0
                up(:,nypo) = 0
                down(:,1) = 0
            else if (boundTreatment == 'mirror') then
                left(1,:) = inField(1,:) 
                right(nxpo,:) = inField(nxpo, :)
                up(:,nypo) = inField(:, nypo)
                down(:,1) = inField(:, 1)
                
            else
                print *, boundTreatment
                stop 'invalid bound condition for derivative'

            endif
        else !(use mirror)
            left(1,:) = inField(1,:) 
            right(nxpo,:) = inField(nxpo, :)
            up(:,nypo) = inField(:, nypo)
            down(:,1) = inField(:, 1)
        endif

        dx = getDxpo() * 1000    ! change from km to meters
        dy = getDypo() * 1000    ! change from km to meters

        gradX = (right - left)/(2*dx)
        gradY = (up-down)/(2*dy)

        DEALLOCATE(left, right, up, down)


    end subroutine
    
    subroutine calc_Ugos_Vgos_taux_tauy()

        !!!!   2020 Oct 27 !!!!
        !!!!   This subroutine is now modified to calulate 
        !!!!   the ocean mixed layer velocity 

        INTEGER(kind = i4) :: nxpo, nypo
        REAL(kind=r8), ALLOCATABLE, DIMENSION(:,:) :: &
                    P_field, & !!! Pressure field
                    ugos, & !!! goestrophic velocity field from pressure
                    vgos, & !!! goestrophic velocity field from pressure
                    umx, &  !!! mixed layer velocity (includes Ekman velocity)
                    vmx, &  !!! mixed layer velocity (includes Ekman velocity)
                    gradX, &
                    gradY, &
                    taux, &
                    tauy, &
                    PPA

        call setOcnPgridXYsizeto(nxpo, nypo)

        ALLOCATE(P_field(nxpo, nypo), & 
                    ugos(nxpo, nypo), &
                    vgos(nxpo, nypo), & 
                    umx(nxpo, nypo), &
                    vmx(nxpo, nypo), & 
                    gradX(nxpo, nypo), &
                    gradY(nxpo, nypo), &
                    taux(nxpo, nypo), &
                    tauy(nxpo, nypo), &
                    PPA(nxpo, nypo))

        call getInputFields(nxpo, nypo, ugos(:,:), vgos(:,:), &
                               taux(:,:), tauy(:,:), PPA(:,:))

        call getPressureField(nxpo,nypo, P_field)
        call calc_grad(nxpo, nypo, P_field, &
                       gradX, &
                       gradY, &
                       boundTreatment='mirror')

        ugos(:,:) = 0.0
        vgos(:,:) = 0.0

        ugos = -gradY/f0
        vgos= gradX/f0

        vmx = vgos - (taux)/(oHm * f0)
        umx = ugos + (tauy)/(oHm * f0)

        taux = rho_o * taux
        tauy = rho_o * tauy

        call saveInputFields(nxpo, nypo, umx(:,:), vmx(:,:), taux(:,:), tauy(:,:))

        print *, 'Saved mixed layer velocity and tau'

        DEALLOCATE(P_field, &
                    ugos, vgos, &
                    umx, vmx, &
                    gradX, gradY, &
                    taux, tauy, &
                    PPA)

    end subroutine

end module
