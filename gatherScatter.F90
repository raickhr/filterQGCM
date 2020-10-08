module gatherScatter
    use kinds
    use mpiMod
    use fields
    use gridMod

    implicit none
    PRIVATE
    save

    public :: broadCastInt, &
              broadCastFloat, &
              syncProcs, &
              bcastVarsAfterReadingInputFile

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

    

end module 
