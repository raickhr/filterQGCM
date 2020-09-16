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
        output2DOcnFields(2)%info%fieldName = 'VVEL'
        output2DOcnFields(3)%info%fieldName = 'TAUX'
        output2DOcnFields(4)%info%fieldName = 'TAUY'
        output2DOcnFields(5)%info%fieldName = 'TotalPowerPerArea'
        output2DOcnFields(6)%info%fieldName = 'MeanPowerPerArea'
        output2DOcnFields(6)%info%fieldName = 'EddyPowerPerArea'
        
        do i=1, nVars2read
            allocate(output2DOcnFields(i)%fieldVal(nxpo,nypo))
        enddo
    end subroutine

end module

