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

            call MPI_BCAST(xpo, nxpo, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(ypo, nypo, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(timeVal, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            
            ! print *,'at taskid',taskid,' received xpo, ypo, timeVal'

            call MPI_BARRIER(MPI_COMM_WORLD, i_err)

            msgSize = nxpo * nypo

            ! print *,'at taskid',taskid,'nxpo, nypo', nxpo, nypo

            !!!! broadcast ugos, vgos, taux, tauy, 
            call MPI_BCAST(UVEL, msgSize, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(VVEL, msgSize, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(TAUX, msgSize, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(TAUY, msgSize, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(PowerPerArea, msgSize, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

            ! print *,'at taskid',taskid,'received UVEL, VVEL, TAUX, TAUY, PowerPerArea'

            call MPI_BARRIER(MPI_COMM_WORLD, i_err)

            if (prevTimeVal == timeVal) then
                print *, 'Found repeated time at rec count:',recDim
                continue
            else
                prevTimeVal = timeVal
            endif

            do filterCounter = 1, nFilterLength

                if (taskid .EQ. MASTER ) then
                    workFilterLen = filterLengthList(filterCounter)
                    write(*,'(A15,F5.0,A5)') 'Filtering at ', workFilterLen,' km'
                    print*, 'Making kernel'
                endif
                call MPI_BCAST(workFilterLen, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BARRIER(MPI_COMM_WORLD, i_err)
                !!! calculate kernel
                call makeKernel(workFilterLen, 'P', ocnORatmGrid='O')
                

                !print * , filterCounter, 'Made Kernel at', taskid

                !!!! filter ugos, vgos, taux, tauy, PowerPerArea
                if (taskid == MASTER) then
                    print *, 'Filtering ...'
                endif

                call reset_FilteredFields()

                OL_UVEL(:,:) = 0.0
                !call getFilterFields(UVEL(:,:), OL_UVEL(:,:))

                OL_VVEL(:,:) = 0.0
                !call getFilterFields(VVEL(:,:), OL_VVEL(:,:))

                OL_TAUX(:,:) = 0.0
                !call getFilterFields(TAUX(:,:), OL_TAUX(:,:))

                OL_TAUY(:,:) = 0.0
                !call getFilterFields(TAUY(:,:), OL_TAUY(:,:))

                OL_PowerPerArea(:,:) = 0.0
                !call getFilterFields(PowerPerArea(:,:), OL_PowerPerArea(:,:))
       
                call filterFields()

                call MPI_BARRIER(MPI_COMM_WORLD, i_err)

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