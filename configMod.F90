module configMod
    use kinds
    use constants
    use gatherScatter
    use mpiMod

    implicit none
    private 
    save

    !! file informations
    CHARACTER (len = char_len_long) :: inputLoc,   & !! input folder and path
                                       outputLoc     !! output folder and path

    INTEGER(kind = i4) :: nxto_c, nyto_c, nzo_c, nxta_c, nyta_c, nza_c
       
    INTEGER(kind = i4) :: nInputFiles,    & !! number of input files
                          nVars2read, &     !! number of variables to read from each file
                          nFilterLength    !! number of filter lengths 

    CHARACTER (len = char_len), ALLOCATABLE, DIMENSION(:) :: inputFileList    !! names of input Files

    CHARACTER (len = char_len) :: gridFileName

    CHARACTER (len = char_len_short), ALLOCATABLE, DIMENSION(:) :: &
                            varNameList                     !! names of input Variables


    REAL(kind =r8) , ALLOCATABLE, DIMENSION(:) :: filterLengthList   !! filter lengths in km

    INTEGER(kind=i4) :: startRecCount, endRecCount  !! record counts to be read from a file
    
    !! Public member function
    PUBLIC :: makeConfig, &
              setAtmXYto, &
              setOcnXYto, &
              getInputLoc, &
              getOutputLoc, &
              numFilesToRead, &
              numVarsToRead, &
              numFilterLength, &
              getInFileListIn, &
              getVarNameListIn, &
              getInFileNo, &
              getVarNameNo, &
              getFilterLenNo, &
              startRecIndx, &
              endRecIndx

    contains 

    subroutine makeConfig()

        INTEGER(kind = i4) :: i , &
                              err


        namelist /paths/      &
        inputLoc,       &
        outputLoc

        namelist /gridFile/ &
        gridFileName

        namelist /ocngridsize/ &
        nxto_c, nyto_c, nzo_c

        namelist /atmgridsize/ &
        nxta_c, nyta_c, nza_c

        namelist /for_allocation/ &
        nInputFiles, &
        nVars2read,    &
        nFilterLength, &
        startRecCount, &
        endRecCount

        namelist /fileName_VarName/ &
        inputFileList, &
        varNameList

        namelist /filterLenList/ &
        filterLengthList

        namelist /coriolis/ &
        f0

        if ( thisProc() == MASTER) then 
        

            print *, ''
            print *, 'CONFIG FILE'
            print *, ''
            print *, 'Input/Output Directories'
            read(*, paths)
            WRITE(* ,'(A25,A50)') 'Input Location : ', inputLoc
            WRITE(* ,'(A25,A50)') 'Output Location : ', outputLoc

            read(*, gridFile)
            print *, ''
            WRITE(* ,'(A25,A50)') 'GridFile : ', gridFileName

    100     format('     ocean grid (nx,ny,nz) : ','(',I4,', ',I4,', ',I4,')')
    101     format('atmosphere grid (nx,ny,nz) : ','(',I4,', ',I4,', ',I4,')')
            read(*, ocngridsize)
            WRITE(* ,100) nxto_c, nyto_c, nzo_c

            read(*, atmgridsize)
            WRITE(* ,101) nxta_c, nyta_c, nza_c

            read(*, for_allocation)
            print *, ''
            print *, 'Reading information ...'
            WRITE(* ,'(A30,I4)') 'number of input files =', nInputFiles
            WRITE(* ,'(A30,I4)') 'number of variables to read =', nVars2read
            WRITE(* ,'(A30,I4)') 'number of filterlengths =', nFilterLength
            WRITE(* ,'(A30,I4,A4,I4)') 'reading from record no ', startRecCount,' to ', endRecCount

            allocate(inputFileList(nInputFiles), &
                    varNameList(nVars2read), &
                    filterLengthList(nFilterLength))

            read(*, fileName_VarName)
            read(*, filterLenList)

            print *, ''
            print *, 'Files to be read ... '
            do i=1, nInputFiles
                inputFileList(i) = trim(adjustl(inputFileList(i)))
                WRITE(*,'(A50)') inputFileList(i)
            end do
            
            print *, ''
            print *, 'Variables to be read ... '
            do i=1, nVars2read
                varNameList(i) = trim(adjustl(varNameList(i)))
                print*,varNameList(i)
            end do

            print *, ''
            print *, 'Filterlengths to be filtered at ... '
            do i=1, nFilterLength
                WRITE(* ,'(f07.2,A3)') filterLengthList(i),'km'
            end do

            read(*, coriolis)
            print *, ''
            print *, 'Reading coriolis paramater ...'
            WRITE(* ,'(A30,E10.2)') 'coriolis parameter(f0) =', f0
        
        endif

        call broadCastInt(nInputFiles, MASTER, err)
        call broadCastInt(nFilterLength, MASTER, err)
        call broadCastInt(startRecCount, MASTER, err)
        call broadCastInt(endRecCount, MASTER, err)

        
    end subroutine makeConfig

    subroutine setAtmXYto(atmX, atmY, atmLyrs)
        integer(kind = i4), INTENT(OUT) :: atmX, atmY, atmLyrs
        atmX = nxta_c
        atmY = nyta_c
        atmLyrs = nza_c
    end subroutine 

    subroutine setOcnXYto(ocnX, ocnY, ocnLyrs)
        integer(kind = i4), INTENT(OUT) :: ocnX, ocnY, ocnLyrs
        ocnX = nxto_c
        ocnY = nyto_c
        ocnLyrs = nzo_c
    end subroutine 

    function getInputLoc() result(resultVal)
        CHARACTER (len = char_len_long) :: resultVal
        resultVal = inputLoc
    end function

    function getOutputLoc() result(resultVal)
        CHARACTER (len = char_len_long) :: resultVal
        resultVal = outputLoc
    end function

    function numFilesToRead() result(returnVal)
        INTEGER(kind = i4) :: returnVal
        returnVal = nInputFiles
    end function

    function numVarsToRead() result(returnVal)
        INTEGER(kind = i4) :: returnVal
        returnVal = nVars2read
    end function

    function numFilterLength() result(returnVal)
        INTEGER(kind = i4) :: returnVal
        returnVal = nFilterLength
    end function

    subroutine getInFileListIn(list)
        CHARACTER (len = char_len), INTENT(OUT) :: list(nInputFiles)
        list(:) = inputFileList(:)
    end subroutine

    subroutine getVarNameListIn(list)
        CHARACTER (len = char_len_short), INTENT(OUT) :: list(nVars2read)
        list(:) = varNameList(:)
    end subroutine

    function getInFileNo(number) result(returnVal)
        INTEGER(kind = i4) :: number
        CHARACTER (len = char_len):: returnVal
        returnVal = inputFileList(number)
    end function

    function getVarNameNo(number) result(returnVal)
        INTEGER(kind = i4) :: number
        CHARACTER (len = char_len_short):: returnVal
        returnVal = varNameList(number)
    end function

    function getFilterLenNo(number) result(returnVal)
        INTEGER(kind = i4) :: number
        REAL(kind=r4) :: returnVal
        returnVal = filterLengthList(number)
    end function

    function startRecIndx() result(returnVal)
        INTEGER(kind = i4) :: returnVal
        returnVal = startRecCount
    end function

    function endRecIndx() result(returnVal)
        INTEGER(kind = i4) :: returnVal
        returnVal = endRecCount
    end function

    

end module
