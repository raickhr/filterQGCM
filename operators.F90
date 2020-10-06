module operators
    use kinds
    use gridMod
    use constants
    implicit none

    contains

    subroutine calc_grad(inField, gradX, gradY, &
                        boundTreatment)

        REAL(kind=r8),INTENT(IN) :: inField(nxpo, nypo)
        REAL(kind=r8),INTENT(OUT) :: gradX(nxpo, nypo), gradY(nxpo, nypo)
        CHARACTER(len=*), OPTIONAL :: boundTreatment

        REAL(kind=r8) :: left(nxpo, nypo), right(nxpo, nypo), & 
                         up(nxpo, nypo), down(nxpo, nypo), &
                         dx, dy

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

        dx = (xpo(2) - xpo(1))*1000    ! change from km to meters
        dy = (ypo(2) - ypo(1))*1000    ! change from km to meters

        gradX = (right - left)/(2*dx)
        gradY = (up-down)/(2*dy)

    end subroutine
    
    subroutine calc_Ugos_Vgos(P_Field, ugos, vgos)

        REAL(kind=r8), INTENT(IN) :: P_field(nxpo, nypo) !!! Pressure field
        REAL(kind=r8), INTENT(OUT) :: ugos(nxpo, nypo), &
                                      vgos(nxpo, nypo) !!! goestrophic field from pressure
        
        REAL(kind=r8) :: gradX(nxpo, nypo), gradY(nxpo, nypo), m1(nxpo, nypo)
        call calc_grad(P_Field, &
                       gradX, &
                       gradY, &
                       boundTreatment='mirror')

        ugos(:,:) = 0.0
        vgos(:,:) = 0.0
        m1(:,:) = -1.0

        ugos = -gradY/f0
        vgos= gradX/f0
    end subroutine

end module
