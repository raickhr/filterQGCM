module workDiv
    use kinds
    use gridMod
    use mpiMod
    implicit none
    private 
    save

    !! variables for workDivison
    integer(kind = i4) :: start_col, &
                          end_col, &
                          ncol

    PUBLIC :: divide_task, &
              get_start_col, &
              get_end_col, &
              get_ncol
    contains

    subroutine divide_task(gridType, ocnORatmGrid)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !         THIS FUNCTION FIRST SHARES THE DIVIDES THE WORK BETWEEN THE WORKER NODES COLUMNWISE
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CHARACTER, INTENT(IN) :: gridType
        CHARACTER, OPTIONAL :: ocnORatmGrid

        INTEGER(kind=i4) :: nx, ny, avg_col, extra_col, &
                            dest, source, numworkers, &
                            err, msg_tag

        numworkers = totalWorkers()

        msg_tag = 0

        if (thisProc() .EQ. MASTER ) then

            print *, 'dividing task ...'
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !         MASTER DIVIDING THE WORK BETWEEN THE WORKERS
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (PRESENT(ocnORatmGrid)) then
                if ((ocnORatmGrid .EQ. 'A') .OR. &
                    (ocnORatmGrid .EQ. 'a')) then
                    if ( gridType .EQ. 'T') then 
                        call setAtmTgridXYsizeto(nx,ny)
                    else
                        call setAtmPgridXYsizeto(nx,ny)
                    endif
                else
                    if ( gridType .EQ. 'T') then 
                        call setOcnTgridXYsizeto(nx,ny)
                    else
                        call setOcnPgridXYsizeto(nx,ny)
                    endif
                endif
            else
                call setOcnPgridXYsizeto(nx,ny)
            endif

            extra_col = mod(ny,numworkers)
            avg_col = ny/numworkers
            start_col = 1

            do dest = 1, numworkers  !! calculating the data size for every workers and sending them to the respective workers
                if (dest .LE. extra_col) then
                    ncol = avg_col +1
                else
                    ncol = avg_col
                endif

                end_col = start_col + ncol -1
            
                call MPI_SEND( start_col, 1, MPI_INTEGER, dest, msg_tag + 1, &
                               MPI_COMM_WORLD, err )
                call MPI_SEND( end_col, 1, MPI_INTEGER, dest, msg_tag + 2, &
                               MPI_COMM_WORLD, err )
                call MPI_SEND( ncol, 1, MPI_INTEGER, dest, msg_tag + 3, &
                               MPI_COMM_WORLD, err )
                
                print *,'proc, start_col, end_col, ncol', dest, start_col, end_col, ncol
                
                start_col = start_col + ncol
            enddo

        else 

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !         RECEIVING THE WORK DIVISON FROM MASTER
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            call MPI_RECV( start_col, 1, MPI_INTEGER, MASTER, &
                            msg_tag + 1, MPI_COMM_WORLD, status, err )

            call MPI_RECV( end_col, 1, MPI_INTEGER, MASTER, &
                            msg_tag + 2, MPI_COMM_WORLD, status, err )

            call MPI_RECV( ncol, 1, MPI_INTEGER, MASTER, &
                            msg_tag + 3, MPI_COMM_WORLD, status, err )

            ! write(*,'(A14 I3 A11 I2 A9 I2 A7 I2)') &
            ! 'at processor ',taskid, ' start_col ', start_col, &
            ! ' end_col ', end_col, ' ncol ', ncol

        end if
    end subroutine divide_task

    function get_start_col() result(returnVal)
        INTEGER(kind=i4)  :: returnVal
        returnVal = start_col
    end function
    
    function get_end_col() result(returnVal)
        INTEGER(kind=i4)  :: returnVal
        returnVal = end_col
    end function

    function get_ncol() result(returnVal)
        INTEGER(kind=i4)  :: returnVal
        returnVal = ncol
    end function





end module
