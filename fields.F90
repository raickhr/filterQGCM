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


    REAL(kind=r8), ALLOCATABLE, DIMENSION(:,:), PRIVATE :: &
                    Pressure, UVEL, VVEL, &
                    TAUX, TAUY, &
                    PowerPerArea, &
                    OL_UVEL, OL_VVEL, &
                    OL_TAUX, OL_TAUY, &
                    OL_PowerPerArea

    PUBLIC :: init_inputOcnFields, &
                init_outputOcnFields, &
                reset_workFields, &
                reset_FilteredFields, &
                saveReadInputFields, &
                saveInputFields, &
                saveOutputFields, &
                getPressureField, &
                getInputFields, &
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


    subroutine reset_workFields()
        INTEGER(kind = i4) :: nxpo, nypo
        call setOcnPgridXYsizeto(nxpo, nypo)

        if (ALLOCATED(UVEL)) then
            DEALLOCATE(Pressure, UVEL, VVEL, &
                       TAUX, TAUY, &
                       PowerPerArea)
        endif

        allocate(Pressure(nxpo, nypo), &
                 UVEL(nxpo, nypo), VVEL(nxpo, nypo), &
                 TAUX(nxpo, nypo), TAUY(nxpo, nypo), &
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

    subroutine saveReadInputFields(nxpo, nypo, inPressure, &
                               inTAUX, inTAUY)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inPressure(nxpo, nypo), &
                                        inTAUX(nxpo, nypo), inTAUY(nxpo, nypo)

        Pressure(:,:) = inPressure(:,:)
        TAUX(:,:) = inTAUX(:,:)
        TAUY(:,:) = inTAUY(:,:)

    end subroutine

    subroutine getPressureField(nxpo, nypo, outPressure)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(OUT) ::  outPressure(nxpo, nypo)

        outPressure(:,:) = Pressure(:,:)

    end subroutine

    subroutine saveInputFields(nxpo, nypo, inUVEL, inVVEL, &
                               inTAUX, inTAUY)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUVEL(nxpo, nypo), inVVEL(nxpo, nypo), &
                                        inTAUX(nxpo, nypo), inTAUY(nxpo, nypo)

        
        UVEL(:,:) = inUVEL(:,:)
        VVEL(:,:) = inVVEL(:,:)
        TAUX(:,:) = inTAUX(:,:)
        TAUY(:,:) = inTAUY(:,:)
        PowerPerArea = UVEL * TAUX + VVEL * TAUY 

    end subroutine

    subroutine saveOutputFields(nxpo, nypo, inUVEL, inVVEL, &
                               inTAUX, inTAUY, inPPA)

        INTEGER(kind = r4), INTENT(IN) :: nxpo, nypo
    
        REAL(kind = r8), INTENT(IN) ::  inUVEL(nxpo, nypo), inVVEL(nxpo, nypo), &
                                        inTAUX(nxpo, nypo), inTAUY(nxpo, nypo), &
                                        inPPA(nxpo, nypo)

        
        OL_UVEL(:,:) = inUVEL(:,:)
        OL_VVEL(:,:) = inVVEL(:,:)
        OL_TAUX(:,:) = inTAUX(:,:)
        OL_TAUY(:,:) = inTAUY(:,:)
        OL_PowerPerArea(:,:) = inPPA(:,:)
        
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

