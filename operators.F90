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
    
    subroutine calc_Ugos_Vgos()

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
                    tauy

        call setOcnPgridXYsizeto(nxpo, nypo)

        ALLOCATE(P_field(nxpo, nypo), & 
                    ugos(nxpo, nypo), &
                    vgos(nxpo, nypo), & 
                    umx(nxpo, nypo), &
                    vmx(nxpo, nypo), & 
                    gradX(nxpo, nypo), &
                    gradY(nxpo, nypo), &
                    taux(nxpo, nypo), &
                    tauy(nxpo, nypo))

        call getTauxTauy(nxpo, nypo, taux(:,:), tauy(:,:))

        call getPressureField(nxpo,nypo, P_field)
        call calc_grad(nxpo, nypo, P_field, &
                       gradX, &
                       gradY, &
                       boundTreatment='mirror')

        ugos(:,:) = 0.0
        vgos(:,:) = 0.0

        ugos = -gradY/f0
        vgos= gradX/f0

        call saveInputOcnUgosVgos(nxpo, nypo, ugos(:,:), vgos(:,:))

        vmx = vgos - (taux)/(oHm * f0)
        umx = ugos + (tauy)/(oHm * f0)

        taux = rho_o * taux !! This tau will be updated later
        tauy = rho_o * tauy

        !! This sets the value of the ocean surface velocity 
        call saveInputOcnUmixVmix(nxpo, nypo, umx(:,:), vmx(:,:))

        print *, 'Saved ocean mixed layer velocity ' 

        DEALLOCATE(P_field, &
                    ugos, vgos, &
                    umx, vmx, &
                    gradX, gradY, &
                    taux, tauy)

    end subroutine

    subroutine calc_taux_tauy()

        INTEGER(kind = i4) :: nxpo, nypo
        REAL(kind=r8), ALLOCATABLE, DIMENSION(:,:) :: &
                    umx, vmx, &     !! mixed layer velocity (includes Ekman velocity)
                    umx_a, vmx_a, & !! atmospheric mixed layer velocity
                    relVelMag, &    !! relative velocity magnitude
                    taux,tauy       !! taux tauy

        call setOcnPgridXYsizeto(nxpo, nypo)

        ALLOCATE(   umx(nxpo, nypo), &
                    vmx(nxpo, nypo), & 
                    umx_a(nxpo, nypo), &
                    vmx_a(nxpo, nypo), &
                    taux(nxpo, nypo), &
                    tauy(nxpo, nypo), & 
                    relVelMag(nxpo, nypo) )

        call getAtmFilteredUmixVmix(nxpo, nypo, umx_a, vmx_a)

        call getOcnUmixVmix(nxpo, nypo, umx, vmx)

        relVelMag = SQRT((umx_a - umx)**2 + (vmx_a - vmx)**2)

        taux = Cd * relVelMag * (umx_a - umx)
        tauy = Cd * relVelMag * (vmx_a - vmx)

        call saveInputTauxTauy(nxpo, nypo, taux(:,:), tauy(:,:))

        print *, 'Calculated and saved taux and tauy .... ' 

        DEALLOCATE(umx, vmx, &
                   umx_a, vmx_a, &
                   taux, tauy, &
                   relVelMag)

    end subroutine

end module
