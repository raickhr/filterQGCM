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

    INTEGER(kind=i4) :: recDim, &
                        fileCounter, filterCounter, &
                        glbGridSizeWithPadX, glbGridSizeWithPadY, &
                        msgSize, & ! mpi message size
                        knx, kny, &! kernel size
                        is, ie, js, je, &
                        dummmy
                        
    REAL(kind=r4) :: workFilterLen

    CHARACTER(len=256) :: outputFileName, str_recDim, str_filterLen
    ! initialize MPI
    call startMPI()

    !!! Reads the i/o location and variables to read
    if (taskid .EQ. MASTER) then 
        call makeConfig()
    endif

    call MPI_BCAST(nInputFiles, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(nFilterLength, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(startRecCount, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(endRecCount, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

    !! initialize grid
    call init_grid() 

    !!! divide the tasks to the workers
    call divide_task('P', ocnORatmGrid='O')

    !! initialize the input fields
    call init_inputOcnFields()

    !! initialize the output fields
    call init_outputOcnFields()


    !! intialize the workFields
    call init_workFields()

    do fileCounter =1, nInputFiles
        !!! loop over record dimnesion for a file
        do recDim = startRecCount, endRecCount
            if (taskid .EQ. MASTER ) then
                call readInputOcnVars2D_P(inputLoc, &
                                          inputFileList(1), & 
                                          1, &
                                          recDim)
                
                !!!! claculate the ugos, vgos from pressure and taux , tauy copy for filtering
                
                call calc_Ugos_Vgos(Pressure, UVEL, VVEL)

                TAUX(:,:) = rho_o * TAUX
                TAUY(:,:) = rho_o * TAUY

                PowerPerArea = UVEL * TAUX + VVEL * TAUY

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
                    write(*,'(A15,F4.0,A5)') 'Filtering at ', workFilterLen,' km'
                    print*, 'Making kernel'
                endif
                call MPI_BCAST(workFilterLen, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BARRIER(MPI_COMM_WORLD, i_err)
                !!! calculate kernel
                call makeKernel(workFilterLen, 'P', ocnORatmGrid='O')
                knx = size(kernel,1, kind=i4)
                kny = size(kernel,2, kind=i4)

                if (taskid .EQ. MASTER ) then !MASTER
                    write(*,'(A15i4i4)') 'kernel size', knx, kny
                    ! do dummmy =1, knx
                    !     write(*,'(13F11.5)') kernel(dummmy, :)
                    ! enddo
                    print *, 'sum Kernel=', sum(kernel)
                endif

                !print * , filterCounter, 'Made Kernel at', taskid

                !!!! filter ugos, vgos, taux, tauy, PowerPerArea
                if (taskid == MASTER) then
                    print *, 'Filtering ...'
                endif

                OL_UVEL(:,:) = 0
                OL_VVEL(:,:) = 0
                OL_TAUX(:,:) = 0
                OL_TAUY(:,:) = 0
                OL_PowerPerArea(:,:) = 0
                
                call getFilterFields(PowerPerArea, OL_PowerPerArea)

                call getFilterFields(TAUX, OL_TAUX)
                call getFilterFields(TAUY, OL_TAUY)

                call getFilterFields(UVEL, OL_UVEL)
                call getFilterFields(VVEL, OL_VVEL)


                !call getFilterFields()
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
