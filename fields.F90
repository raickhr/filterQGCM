module fields
    use kinds
    use gridMod
    use configMod

    implicit none
    type field_info                       
        character (len=char_len_short):: fieldName
        character (len=char_len_short):: longName
        character (len=char_len_short):: units
    end type

    type field2D
        type (field_info) info
        real(kind=r8), allocatable, dimension(:,:) :: fieldVal
    end type


    type(field2D), ALLOCATABLE, DIMENSION(:):: input2DOcnFields
    type(field2D), ALLOCATABLE, DIMENSION(:):: output2DOcnFields
    REAL(kind=r8) :: timeVal
    CHARACTER(len = char_len_short) :: timeUnits 

    contains

    subroutine init_inputOcnFields()
        INTEGER(kind=i4):: i
        allocate(input2DOcnFields(nVars2read))
        do i=1, nVars2read
            input2DOcnFields(i)%info%fieldName = trim(adjustl(varNameList(i)))
            allocate(input2DOcnFields(i)%fieldVal(nxpo,nypo))
        enddo
    end subroutine

    subroutine init_outputOcnFields()
        INTEGER(kind=i4):: i, nOutputFields
        nOutputFields = 7
        allocate(output2DOcnFields(nOutputFields))

        output2DOcnFields(1)%info%fieldName = 'UVEL'
        output2DOcnFields(1)%info%longName = 'filtered zonal velocity'
        output2DOcnFields(1)%info%units = 'm/s'

        output2DOcnFields(2)%info%fieldName = 'VVEL'
        output2DOcnFields(2)%info%longName = 'filtered meridional velocity'
        output2DOcnFields(2)%info%units = 'm/s'

        output2DOcnFields(3)%info%fieldName = 'TAUX'
        output2DOcnFields(3)%info%longName = 'filtered zonal surface wind stress'
        output2DOcnFields(3)%info%units = 'Pascal'

        output2DOcnFields(4)%info%fieldName = 'TAUY'
        output2DOcnFields(4)%info%longName = 'filtered meriodional surface wind stress'
        output2DOcnFields(4)%info%units = 'Pascal'

        output2DOcnFields(5)%info%fieldName = 'TotalPowerPerArea'
        output2DOcnFields(5)%info%longName = '\overline{\tau . u}'
        output2DOcnFields(5)%info%units = 'Watt/m^2'

        output2DOcnFields(6)%info%fieldName = 'MeanPowerPerArea'
        output2DOcnFields(6)%info%longName = '\overline{\tau} . \overline{u}'
        output2DOcnFields(6)%info%units = 'Watt/m^2'

        output2DOcnFields(7)%info%fieldName = 'EddyPowerPerArea'
        output2DOcnFields(7)%info%longName = '\overline{\tau . u} - \overline{\tau} . \overline{u}'
        output2DOcnFields(7)%info%units = 'Watt/m^2'
        
        do i=1, nVars2read
            allocate(output2DOcnFields(i)%fieldVal(nxpo,nypo))
        enddo
    end subroutine

end module

