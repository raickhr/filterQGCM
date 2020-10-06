        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !         THIS FUNCTION DIVIDES THE WORK BETWEEN THE WORKER NODES COLUMNWISE
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CHARACTER, INTENT(IN) :: gridType
        CHARACTER, OPTIONAL :: ocnORatmGrid

        if (taskid .EQ. MASTER ) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !         MASTER DIVIDING THE WORK BETWEEN THE WORKERS
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (PRESENT(ocnORatmGrid)) then
                if (ocnORatmGrid .EQ. 'A' .OR.
                    ocnORatmGrid .EQ. 'a' ) then
                        if ( gridType .EQ. 'T') then 
                            nx = nxta
                            ny = nyta
                        else
                            nx = nxpa
                            ny = nypa
                        endif
                else
                    if ( gridType .EQ. 'T') then 
                        nx = nxto
                        ny = nyto
                    else
                        nx = nxpo
                        ny = nypo
                    endif

                endif
            else
                nx = nxpo
                ny = nypo
            endif

            extra_col = mod(ny,numworkers)
            avg_col = ny/numworkers
            offset = 1

            do dest = 1, numworkers  !! calculating the data size for every 
                                     !! workers and sending them to the respective workers
                
                if (dest .LE. extra_col) then
                    ncol = avg_col +1
                else
                    ncol = avg_col
                endif

                call MPI_SEND( offset, 1, MPI_INTEGER, dest, msg_tag, &
                &              MPI_COMM_WORLD, i_err )
                call MPI_SEND( ncol, 1, MPI_INTEGER, dest, msg_tag, &
                &              MPI_COMM_WORLD, i_err )

                offset = offset + ncol

            !print *,'offset and ncol sent to taskid', dest


            enddo

        else 

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !         RECEIVING THE WORK DIVISON FROM MASTER
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
                call MPI_RECV( offset, 1, MPI_INTEGER, MASTER, &
                               msg_tag, MPI_COMM_WORLD, status, i_err )

                call MPI_RECV( ncol, 1, MPI_INTEGER, MASTER, &
                              msg_tag, MPI_COMM_WORLD, status, i_err )

                start_col = offset
                end_col = start_col + ncol - 1

                data_size = nx * ncol

        end if