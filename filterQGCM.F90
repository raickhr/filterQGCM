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

    INTEGER(kind = i4) :: fileCounter, recDim, flenCounter, atmFlenCounter, errorCode
    REAL(kind = r8) :: prevTimeVal, currTimeVal, atmFilterLength, filterLength
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

    !! initialize the output field name and var for atm umix vmix
    call init_outputAtmOvrOcnFields()

    !!! loop over file counter
    do fileCounter =1, numFilesToRead()

        !!! loop over record dimnesion for a file
        do recDim = startRecIndx(), endRecIndx()
            !! intialize the workFields
            call reset_workFields()
            
            if (thisProc() .EQ. MASTER ) then
                call readInputOcnVars2D_P(fileCounter, recDim)
                
                !!!! claculate the ugos, vgos from pressure and taux , tauy copy for filtering    
                call calc_Ugos_Vgos()

                ! print *, 'ugos, vgos, taux, tauy obtained ... '

                ! call writeAfterReadingAtmUmixVmix(recDim)
                
                
            endif

            call syncProcs(errorCode)

            call bcastDimensionVars(errorCode)

            prevTimeVal = getPrevTimeVal()       
            currTimeVal = getCurrTimeVal()

            call syncProcs(errorCode)
            if (prevTimeVal == currTimeVal) then
                if (thisProc() .EQ. MASTER) then
                    print *, 'Found repeated time at rec count:',recDim
                endif
                continue
            endif
            call syncProcs(errorCode)

            call bcastUatmVatmBeforeFiltering(errorCode)

            !!! loop to coarsen the atmosphere
            do atmFlenCounter = 1, numAtmFilterLength()

                if (thisProc() .EQ. MASTER ) then
                    atmFilterLength = getAtmFilterLenNo(atmFlenCounter)
                    write(*,'(A25,F5.0,A5)') 'Filtering atmosphere at ', atmFilterLength,' km'
                    print*, 'Making kernel for atm filtering'
                endif
                call broadCastFloat(atmFilterLength, MASTER, errorCode)

                call syncProcs(errorCode)
        
                !!! calculate kernel for atmosphere 
                call makeAtmKernel(atmFilterLength, 'P', ocnORatmGrid='O')

                if (thisProc() == MASTER) then
                    print *, 'Filtering atmospheric mixed layered velocity ...'
                endif

                call filterAtmUmixVmix()

                call syncProcs(errorCode)

                if (thisProc() .EQ. MASTER ) then
                    call writeAppendFilteredAtmUmixVmix(recDim, atmFilterLength)
                    write(*,'(A20)') 'Calculating tau '
                    call calc_taux_tauy()
                endif

                !make atmKernel
                !filter atm mix
                !calculate tau in master
                !broadcast tau

                call bcastAllVarsBeforeFiltering(errorCode)

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
                            call writeAppendOcnOutPut_Pfields(recDim, atmFilterLength, filterLength )
                        endif
                        !call MPI_Barrier(  MPI_COMM_WORLD, i_err)

                    enddo !! filterLength
            enddo   !! atmfilterLength
        enddo  !! recDim
    enddo
    print *, 'Finalizing MPI'
    call endMPI()
    print *, 'end program'

end program filterQGCM
