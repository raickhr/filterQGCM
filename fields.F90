module fields
    use kinds
    use gridMod
    use configMod

    implicit none
    PRIVATE
    SAVE

    type field_info                       
        character (len=char_len_short):: fieldName
        character (len=char_len_short):: longName
        character (len=char_len_short):: units
    end type


    type(field_info), ALLOCATABLE, DIMENSION(:), PUBLIC:: input2DOcnFields
    type(field_info), ALLOCATABLE, DIMENSION(:), PUBLIC:: output2DOcnFields
    type(field_info), ALLOCATABLE, DIMENSION(:), PUBLIC:: output2DAtmOvrOcnFields


    REAL(kind=r8), ALLOCATABLE, DIMENSION(:,:), PRIVATE :: &
                    Pressure, UVEL, VVEL, &    !! ocean surface mixed layer velocity
                    UGOS, VGOS, &
                    TAUX, TAUY, &              !! wind stress
                    UATM, VATM, &              !! atmosphere surface mixed layer velocity 
                    PowerPerArea, &            !! power per area    
                    OL_UVEL, OL_VVEL, &        !! filtered ocean surface mixed layer velocity
                    OL_UATM, OL_VATM, &        !! filtered atmospheric surface mixed layer velocity
                    OL_TAUX, OL_TAUY, &        !! filtered surface wind stress
                    OL_PowerPerArea            !! filtered power per area

    PUBLIC ::   init_inputOcnFields,    &
                init_outputOcnFields,    &
                init_outputAtmOvrOcnFields,    &
                saveReadInputFields,    &
                getPressureField,    &
                getTauxTauy,    &
                getAtmUmixVmix,    &
                getAtmFilteredUmixVmix,    &
                reset_workFields,    &
                reset_FilteredFields,    &
                reset_FilteredAtmMixedLayerVel,    &
                saveInputOcnUmixVmix,    &
                saveInputOcnUgosVgos, &
                saveInputTauxTauy,    &
                saveInputAtmUmixVmix,    &
                saveUnfilteredFields,    &
                saveAtmFilteredUmixVmix,    &
                saveOutputFields,    &
                getInputFields,    &
                getOcnUmixVmix,    &
                getOcnUgosVgos,    &
                getOutputFields 

    contains

    subroutine init_inputOcnFields()
        INTEGER(kind=i4):: i, nVars2read
        CHARACTER (len = char_len_short), ALLOCATABLE, DIMENSION(:) :: varNameList

        allocate(varNameList(numVarsToRead()), &
                 input2DOcnFields(numVarsToRead()))

        call getVarNameListIn(varNameList)

        do i=1, numVarsToRead() 
            input2DOcnFields(i)%fieldName = trim(adjustl(varNameList(i)))
        enddo

        DEALLOCATE(varNameList)

    end subroutine

    subroutine init_outputOcnFields()
        INTEGER(kind=i4):: i, nOutputFields
        nOutputFields = 7
        allocate(output2DOcnFields(nOutputFields))

        output2DOcnFields(1)%fieldName = 'UVEL'
        output2DOcnFields(1)%longName = 'filtered zonal velocity'
        output2DOcnFields(1)%units = 'm/s'

        output2DOcnFields(2)%fieldName = 'VVEL'
        output2DOcnFields(2)%longName = 'filtered meridional velocity'
        output2DOcnFields(2)%units = 'm/s'

        output2DOcnFields(3)%fieldName = 'TAUX'
        output2DOcnFields(3)%longName = 'filtered zonal surface wind stress'
        output2DOcnFields(3)%units = 'Pascal'

        output2DOcnFields(4)%fieldName = 'TAUY'
        output2DOcnFields(4)%longName = 'filtered meriodional surface wind stress'
        output2DOcnFields(4)%units = 'Pascal'

        output2DOcnFields(5)%fieldName = 'TotalPowerPerArea'
        output2DOcnFields(5)%longName = '\overline{\tau . u}'
        output2DOcnFields(5)%units = 'Watt/m^2'

        output2DOcnFields(6)%fieldName = 'MeanPowerPerArea'
        output2DOcnFields(6)%longName = '\overline{\tau} . \overline{u}'
        output2DOcnFields(6)%units = 'Watt/m^2'

        output2DOcnFields(7)%fieldName = 'EddyPowerPerArea'
        output2DOcnFields(7)%longName = '\overline{\tau . u} - \overline{\tau} . \overline{u}'
        output2DOcnFields(7)%units = 'Watt/m^2'
        
        ! do i=1, 7
        !     allocate(output2DOcnFields(i)%fieldVal(nxpo,nypo))
        ! enddo
    end subroutine

    subroutine init_outputAtmOvrOcnFields()
        INTEGER(kind=i4):: i, nOutputFields
        nOutputFields = 2
        allocate(output2DAtmOvrOcnFields(nOutputFields))

        output2DAtmOvrOcnFields(1)%fieldName = 'UVEL'
        output2DAtmOvrOcnFields(1)%longName = 'filtered atm. zonal vel'
        output2DAtmOvrOcnFields(1)%units = 'm/s'

        output2DAtmOvrOcnFields(2)%fieldName = 'VVEL'
        output2DAtmOvrOcnFields(2)%longName = 'filtered atm. meridional vel'
        output2DAtmOvrOcnFields(2)%units = 'm/s'

    end subroutine

    subroutine saveReadInputFields(nxpo, nypo, inPressure, &
                               inTAUX, inTAUY, inUmix, inVmix)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inPressure(nxpo, nypo), &
                                        inTAUX(nxpo, nypo), inTAUY(nxpo, nypo), &
                                        inUmix(nxpo, nypo), inVmix(nxpo, nypo)

        Pressure(:,:) = inPressure(:,:)
        TAUX(:,:) = inTAUX(:,:)
        TAUY(:,:) = inTAUY(:,:)
        UATM(:,:) = inUmix(:,:)
        VATM(:,:) = inVmix(:,:)

    end subroutine

    subroutine getPressureField(nxpo, nypo, outPressure)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outPressure(nxpo, nypo)

        outPressure(:,:) = Pressure(:,:)

    end subroutine

    subroutine getTauxTauy(nxpo, nypo, outTaux, outTauy)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outTaux(nxpo, nypo), outTauy(nxpo, nypo)

        outTaux(:,:) = TAUX(:,:)
        outTauy(:,:) = TAUY(:,:)

    end subroutine

    subroutine getAtmUmixVmix(nxpo, nypo, outUATM, outVATM)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outUATM(nxpo, nypo), outVATM(nxpo, nypo)

        outUATM(:,:) = UATM(:,:)
        outVATM(:,:) = VATM(:,:)
    end subroutine

    subroutine getAtmFilteredUmixVmix(nxpo, nypo, outUATM, outVATM)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outUATM(nxpo, nypo), outVATM(nxpo, nypo)

        outUATM(:,:) = OL_UATM(:,:)
        outVATM(:,:) = OL_VATM(:,:)
    end subroutine

    subroutine reset_workFields()
        INTEGER(kind = i4) :: nxpo, nypo
        call setOcnPgridXYsizeto(nxpo, nypo)

        if (ALLOCATED(UVEL)) then
            DEALLOCATE(Pressure, UVEL, VVEL, &
                       UGOS, VGOS, &
                       TAUX, TAUY, &
                       UATM, VATM, &
                       PowerPerArea)
        endif

        allocate(Pressure(nxpo, nypo), &
                 UVEL(nxpo, nypo), VVEL(nxpo, nypo), &
                 UGOS(nxpo, nypo), VGOS(nxpo, nypo), &
                 TAUX(nxpo, nypo), TAUY(nxpo, nypo), &
                 UATM(nxpo, nypo), VATM(nxpo, nypo), &   
                 PowerPerArea(nxpo, nypo))

    end subroutine

    subroutine reset_FilteredFields()
        INTEGER(kind = i4) :: nxpo, nypo
        call setOcnPgridXYsizeto(nxpo, nypo)
        
        if (ALLOCATED(OL_UVEL)) then
            DEALLOCATE(OL_UVEL, OL_VVEL, &
                       OL_TAUX, OL_TAUY, &
                       OL_PowerPerArea)
        endif

        allocate(OL_UVEL(nxpo, nypo), OL_VVEL(nxpo, nypo), &
                 OL_TAUX(nxpo, nypo), OL_TAUY(nxpo, nypo), &
                 OL_PowerPerArea(nxpo, nypo))
    endsubroutine

    subroutine reset_FilteredAtmMixedLayerVel()
        INTEGER(kind = i4) :: nxpo, nypo
        call setOcnPgridXYsizeto(nxpo, nypo)
        
        if (ALLOCATED(OL_UATM)) then
            DEALLOCATE(OL_UATM, OL_VATM)
        endif

        allocate(OL_UATM(nxpo, nypo), OL_VATM(nxpo, nypo))
    end subroutine 
    
    subroutine saveInputOcnUmixVmix(nxpo, nypo, inUVEL, inVVEL)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUVEL(nxpo, nypo), inVVEL(nxpo, nypo)

        UVEL(:,:) = inUVEL(:,:)
        VVEL(:,:) = inVVEL(:,:)

    end subroutine

    subroutine saveInputOcnUgosVgos(nxpo, nypo, inUVEL, inVVEL)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUVEL(nxpo, nypo), inVVEL(nxpo, nypo)

        UGOS(:,:) = inUVEL(:,:)
        VGOS(:,:) = inVVEL(:,:)

    end subroutine

    subroutine saveInputTauxTauy(nxpo, nypo, &
                               inTAUX, inTAUY)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inTAUX(nxpo, nypo), inTAUY(nxpo, nypo)

        TAUX(:,:) = inTAUX(:,:)
        TAUY(:,:) = inTAUY(:,:)

    end subroutine

    subroutine saveInputAtmUmixVmix(nxpo, nypo, &
                               inUATM, inVATM)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUATM(nxpo, nypo), inVATM(nxpo, nypo)

        UATM(:,:) = inUATM(:,:)
        VATM(:,:) = inVATM(:,:)

    end subroutine

    subroutine saveUnfilteredFields(nxpo, nypo, inUVEL, inVVEL, &
                               inTAUX, inTAUY, inPPA)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUVEL(nxpo, nypo), inVVEL(nxpo, nypo), &
                                        inTAUX(nxpo, nypo), inTAUY(nxpo, nypo), &
                                        inPPA(nxpo, nypo)

        
        UVEL(:,:) = inUVEL(:,:)
        VVEL(:,:) = inVVEL(:,:)
        
        TAUX(:,:) = inTAUX(:,:)
        TAUY(:,:) = inTAUY(:,:)

        PowerPerArea(:,:) = inPPA(:,:)

    end subroutine
    
    subroutine saveAtmFilteredUmixVmix(nxpo, nypo, inUVEL, inVVEL)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUVEL(nxpo, nypo), inVVEL(nxpo, nypo)

        OL_UATM(:,:) = inUVEL(:,:)
        OL_VATM(:,:) = inVVEL(:,:)

    end subroutine

    subroutine saveOutputFields(nxpo, nypo, inUVEL, inVVEL, &
                               inTAUX, inTAUY, inPPA)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN), DIMENSION(:, :) ::  inUVEL, inVVEL, &
                                        inTAUX, inTAUY, &
                                        inPPA

        
        OL_UVEL = inUVEL
        OL_VVEL = inVVEL
        OL_TAUX = inTAUX
        OL_TAUY = inTAUY
        OL_PowerPerArea = inPPA
        
    end subroutine

    subroutine getInputFields(nxpo, nypo, outUVEL, outVVEL, &
                               outTAUX, outTAUY, outPPA)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outUVEL(nxpo, nypo), outVVEL(nxpo, nypo), &
                                        outTAUX(nxpo, nypo), outTAUY(nxpo, nypo), &
                                        outPPA(nxpo, nypo)

        
        outUVEL(:,:) = UVEL(:,:)
        outVVEL(:,:) = VVEL(:,:)
        outTAUX(:,:) = TAUX(:,:)
        outTAUY(:,:) = TAUY(:,:)
        outPPA(:,:) = PowerPerArea(:,:)

    end subroutine

    
    subroutine getOcnUmixVmix(nxpo, nypo, outUVEL, outVVEL)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outUVEL(nxpo, nypo), outVVEL(nxpo, nypo)

        outUVEL(:,:) = UVEL(:,:)
        outVVEL(:,:) = VVEL(:,:)
    end subroutine

    subroutine getOcnUgosVgos(nxpo, nypo, outUVEL, outVVEL)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outUVEL(nxpo, nypo), outVVEL(nxpo, nypo)

        outUVEL(:,:) = UGOS(:,:)
        outVVEL(:,:) = VGOS(:,:)
    end subroutine


    subroutine getOutputFields(nxpo, nypo, outUVEL, outVVEL, &
                               outTAUX, outTAUY, outPPA)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outUVEL(nxpo, nypo), outVVEL(nxpo, nypo), &
                                        outTAUX(nxpo, nypo), outTAUY(nxpo, nypo), &
                                        outPPA(nxpo, nypo)

        
        outUVEL(:,:) = OL_UVEL(:,:)
        outVVEL(:,:) = OL_VVEL(:,:)
        outTAUX(:,:) = OL_TAUX(:,:)
        outTAUY(:,:) = OL_TAUY(:,:)
        outPPA(:,:) = OL_PowerPerArea(:,:)

    end subroutine
    

end module

