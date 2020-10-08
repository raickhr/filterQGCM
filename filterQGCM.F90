program filterQGCM
    use configMod
    use constants
    use fields
    use filter
    use gridMod
    use kinds
    use mpiMod
    use netCDFio
    use operators
    use workDiv    

    implicit none

    INTEGER(kind = i4) :: fileCounter, recDim, flenCounter
    REAL(kind = r8) :: prevTimeVal, currTimeVal
    ! initialize MPI
    call startMPI()

    !!! Reads the i/o location and variables to read
    call makeConfig()

    !! initialize grid
    call init_grid() 

    !!! divide the tasks to the workers
    call divide_task('P', ocnORatmGrid='O')

    !! initialize the input fields
    call init_inputOcnFields()

    !! initialize the output fields
    call init_outputOcnFields()

    !!! loop over file counter
    do fileCounter =1, numFilesToRead()

        !!! loop over record dimnesion for a file
        do recDim = startRecIndx(), endRecIndx()
            !! intialize the workFields
            call reset_workFields()
            
            if (thisProc() .EQ. MASTER ) then
                call readInputOcnVars2D_P(fileCounter, recDim)
                
                !!!! claculate the ugos, vgos from pressure and taux , tauy copy for filtering
                
                call calc_Ugos_Vgos_taux_tauy()

                ! print *, 'ugos, vgos, taux, tauy obtained ... '
                
            endif

            call bcastVarsAfterReadingInputFile()

            prevTimeVal = getPrevTimeVal()       !! this sets different time val to different processors
            call setPrevTimeVal(prevTimeVal)     !! setting same timeVel to all processors
            
            currTimeVal = getCurrTimeVal()       !! this sets different time val to different processors
            call setCurrTimeVal(currTimeVal)     !! setting same timeVal to all processors
            
            if (prevTimeVal == currTimeVal) then
                if (thisProc() .EQ. MASTER) then
                    print *, 'Found repeated time at rec count:',recDim
                endif
                continue
            else
                call setPrevTimeVal(currTimeVal)
            endif

            do flenCounter = 1, numFilterLength()

                if (taskid .EQ. MASTER ) then
                    write(*,'(A15,F5.0,A5)') 'Filtering at ', getFilterLenNo(flenCounter),' km'
                    print*, 'Making kernel'
                endif
                !!! calculate kernel
                call makeKernel(getFilterLenNo(flenCounter), 'P', ocnORatmGrid='O')
                !print * , flenCounter, 'Made Kernel at', taskid

                !!!! filter ugos, vgos, taux, tauy, PowerPerArea
                if (taskid == MASTER) then
                    print *, 'Filtering ...'
                endif

                call filterFields()

                if (taskid .EQ. MASTER ) then
                    !!!! update the output field
                    !call setOutputFields()
                    write(str_recDim,'(i4.4)') recDim
                    write(str_filterLen,'(i4.4)') int(workFilterLen)
                    outputFileName = 'recNo_'//trim(str_recDim)//'.filterLen_'//trim(str_filterLen)//'.nc'
                    call writeOcnOutPut_Pfields(outputLoc, &
                                                outputFileName)
                endif
                call MPI_Barrier(  MPI_COMM_WORLD, i_err)

            enddo !! filterLength
            
            
        enddo  !! recDim
    enddo

    call MPI_Finalize()

end program filterQGCM