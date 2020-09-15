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

end module

