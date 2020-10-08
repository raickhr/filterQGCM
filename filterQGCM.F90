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
    use gatherScatter    

    implicit none

    INTEGER(kind = i4) :: fileCounter, recDim, flenCounter, errorCode
    REAL(kind = r8) :: prevTimeVal, currTimeVal, filterLength
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

            call bcastVarsAfterReadingInputFile(errorCode)
            
            prevTimeVal = getPrevTimeVal()       !! this sets different time val to different processors
            call broadCastFloat(prevTimeVal,MASTER,errorCode)     
            call setPrevTimeVal(prevTimeVal)     !! setting same timeVel to all processors
            
            currTimeVal = getCurrTimeVal()       !! this sets different time val to different processors
            call broadCastFloat(currTimeVal,MASTER,errorCode)     
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

                
                
                if (thisProc() .EQ. MASTER ) then
                    filterLength = getFilterLenNo(flenCounter)
                    write(*,'(A15,F5.0,A5)') 'Filtering at ', filterLength,' km'
                    print*, 'Making kernel'
                endif
                call broadCastFloat(filterLength, MASTER, errorCode)

                call syncProcs(errorCode)

                if (thisProc() .EQ. MASTER) then 
                    print *, 'before making kernel'
                endif

                !!! calculate kernel
                call makeKernel(filterLength, 'P', ocnORatmGrid='O')
                !print * , flenCounter, 'Made Kernel at', taskid

                !!!! filter ugos, vgos, taux, tauy, PowerPerArea
                if (thisProc() == MASTER) then
                    print *, 'Filtering ...'
                endif

                call filterFields()

                if (thisProc() .EQ. MASTER ) then
                    call writeOcnOutPut_Pfields(recDim, filterLength )
                endif
                !call MPI_Barrier(  MPI_COMM_WORLD, i_err)

            enddo !! filterLength
            
            
        enddo  !! recDim
    enddo

    call MPI_Finalize()

end program filterQGCM