program filterQGCM
    use fields
    use configMod
    !use mpiMod
    use netCDFio
    use gridMod
    implicit none

    INTEGER(kind=i4) :: recDim, nRecs
    !! CHARACTER(len=char_len) :: fullFileName
    !! call startMPI()   ! initialize MPI

    call makeConfig() !!! Reads the i/o location and variables to read 
    call init_grid()
    call init_inputOcnFields()
    call init_outputOcnFields()


    !!! loop over record dimnesion for a file
    nRecs = 1
    do recDim = 1, nRecs
        call readInputOcnVars2D_P(inputLoc, &
                              inputFileList(1), & 
                              1, &
                              recDim)

        call writeOcnOutPut_Pfields(outputLoc, &
                                    'testOutput.nc')
    enddo

end program filterQGCM
