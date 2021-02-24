module gatherScatter
    use kinds
    use mpiMod
    use fields
    use gridMod
    use workDiv

    implicit none
    PRIVATE
    save

    public :: broadCastInt, &
              broadCastFloat, &
              syncProcs, &
              bcastVarsAfterReadingInputFile, &
              gatherFilteredField

    contains

    subroutine broadCastInt(intVar, sourceProc, errorCode)
        INTEGER(kind = i4), INTENT(IN) :: intVar, sourceProc
        INTEGER(kind = i4), INTENT(OUT) :: errorCode
        call MPI_BCAST(intVar, 1, MPI_INTEGER , sourceProc, MPI_COMM_WORLD, errorCode)
    end subroutine

    subroutine broadCastFloat(fltVar, sourceProc, errorCode)
        REAL(kind = r4), INTENT(IN) :: fltVar
        INTEGER(kind = i4), INTENT(OUT) :: sourceProc, errorCode
        call MPI_BCAST(fltVar, 1, MPI_REAL , sourceProc, MPI_COMM_WORLD, errorCode)
    end subroutine

    subroutine syncProcs(errorCode)
        INTEGER(kind = i4), INTENT(OUT) :: errorCode
        call MPI_Barrier(MPI_COMM_WORLD, errorCode)
    end subroutine 


    subroutine bcastVarsAfterReadingInputFile(i_err)
        INTEGER(kind = i4), INTENT(OUT) :: i_err

        INTEGER(kind = i4) :: nxpo, nypo, msgSize

        REAL(kind = r8) :: timeVal

        REAL(kind = r8), ALLOCATABLE, DIMENSION(:) :: xpo, ypo

        REAL(kind = r8), ALLOCATABLE, DIMENSION(:,:) :: UVEL, VVEL, TAUX, TAUY, PPA

        call setOcnPgridXYsizeto(nxpo, nypo)

        ALLOCATE(xpo(nxpo), ypo(nypo))

        ALLOCATE(UVEL(nxpo, nypo), VVEL(nxpo, nypo), &
                 TAUX(nxpo, nypo), TAUY(nxpo, nypo), &
                 PPA(nxpo, nypo))

        call getInputFields(nxpo, nypo, UVEL, VVEL, &
                               TAUX, TAUY, PPA)

        call getXpoYpo(xpo, ypo)

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

        ! print *,'at taskid',taskid,'received UVEL, VVEL, TAUX, TAUY, PowerPerArea'

        call MPI_BARRIER(MPI_COMM_WORLD, i_err)

        call saveInputFields(nxpo, nypo, UVEL, VVEL, TAUX, TAUY)
        call saveReadXpoYpo(xpo, ypo)

        DEALLOCATE(xpo, ypo, &
                   UVEL, VVEL, TAUX, TAUY, PPA)
    end subroutine


    subroutine gatherFilteredField(startXoffset,  &
                                   startYoffset,  &
                                   field,         &
                                   refField,      &
                                   fieldName)

        INTEGER(kind = i4), INTENT(IN) :: startXoffset, startYoffset
        REAL(kind = r8), INTENT(INOUT), DIMENSION(:, :) :: field, refField
        CHARACTER(len = *), INTENT(IN), OPTIONAL :: fieldName

        INTEGER(kind = i4) :: source, dest, numworkers, msg_tag, &
                              nxpo, nypo, & 
                              startCol, endCol, nCol, &
                              xcolSize, datasize, i_err, &
                              i, j  !! dummy iterators

        REAL(KIND=r8), ALLOCATABLE, DIMENSION(:,:) :: buffer_msg
        msg_tag = 0

        call setOcnPgridXYsizeto(nxpo, nypo)

        numworkers = totalWorkers()

        if (thisProc() == MASTER) then
            field(:, :) = -9999
            do source = 1, numworkers
                call MPI_RECV( startCol, 1, MPI_INTEGER, source, msg_tag , &
                &                 MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()',thisProc(),' send start_col'

                call MPI_RECV( endCol, 1, MPI_INTEGER, source, msg_tag , &
                &                 MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()',thisProc(),' send end_col'

                call MPI_RECV( nCol, 1, MPI_INTEGER, source, msg_tag , &
                &                 MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()',thisProc(),' send end_col'

                call MPI_RECV( datasize, 1, MPI_INTEGER, source, msg_tag , &
                &                 MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()',thisProc(),' send data size'

                ALLOCATE(buffer_msg(nxpo + 2*startXoffset, ncol))
                
                call MPI_RECV( buffer_msg(:,:), datasize, MPI_REAL, source, &
                & msg_tag , MPI_COMM_WORLD, status, i_err )
                !print *, 'from thisProc()',thisProc(),' send outUVEL'
                field(1:nxpo+2*startXoffset, startCol:endCol) = buffer_msg(:,:)
                DEALLOCATE(buffer_msg)

                ! call MPI_RECV( field(1, startCol), datasize, MPI_REAL, source, &
                ! & msg_tag , MPI_COMM_WORLD, status, i_err )
                ! !print *, 'from thisProc()',thisProc(),' send outUVEL'
            
                !

                
            end do

            

        else
            dest = MASTER

            startCol = get_start_col() + startYoffset
            endCol = get_end_col() + startYoffset
            ncol = get_ncol()

            datasize = (2* startXoffset + nxpo) * ncol

            call MPI_SEND( startCol, 1, MPI_INTEGER, dest, msg_tag , &
            &                 MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send start_col'

            call MPI_SEND( endCol, 1, MPI_INTEGER, dest, msg_tag , &
            &                 MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send end_col'

            call MPI_SEND( nCol, 1, MPI_INTEGER, dest, msg_tag , &
            &                 MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send num_col'



            call MPI_SEND( datasize, 1, MPI_INTEGER, dest, msg_tag , &
            &                 MPI_COMM_WORLD, status, i_err )
            !print *, 'from thisProc()',thisProc(),' send data size'

            ALLOCATE(buffer_msg(nxpo + 2*startXoffset, ncol))

            buffer_msg(:,:) = field(1:nxpo + 2*startXoffset,startCol:endCol)

            call MPI_SEND( buffer_msg(:,:), datasize, MPI_REAL, dest, &
            & msg_tag , MPI_COMM_WORLD, status, i_err )
            ! print *, 'from thisProc()',thisProc(),' send outUVEL' 
            DEALLOCATE(buffer_msg)  

            ! call MPI_SEND( field(:,startCol:endCol), datasize, MPI_REAL, dest, &
            !  & msg_tag , MPI_COMM_WORLD, status, i_err )
            
            

        endif
        
        ! if (thisProc() == MASTER) then
        !     do j = startYoffset + 1 , startYoffset + nypo
        !         do i = startXoffset + 1, startXoffset + nxpo
        !             if (field(i,j) .NE. refField(i,j)) then 
        !                 if (PRESENT(fieldName)) then
        !                     print *, fieldName
        !                 endif
        !                 print *,i-startXoffset,j-startYoffset,' gather field error'
        !                 call mpiAbort()
        !             endif
        !         enddo
        !     enddo
        ! endif

    end subroutine

    

end module 
