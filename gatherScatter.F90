module gatherScatter
    use kinds
    use mpiMod

    implicit none
    PRIVATE
    save

    public :: broadCastInt, &
              broadCastFloat, &
              syncProcs

    contains

    subroutine broadCastInt(intVar, sourceProc, errorCode)
        INTEGER(kind = i4), INTENT(IN) :: intVar, sourceProc
        INTEGER(kind = i4), INTENT(OUT) :: errorCode
        call MPI_BCAST(intVar, 1, MPI_INTEGER , sourceProc, MPI_COMM_WORLD, errorCode)
    end subroutine

    subroutine broadCastFloat(fltVar, sourceProc, errorCode)
        REAL(kind = r4), INTENT(IN) :: fltVar, sourceProc
        INTEGER(kind = i4), INTENT(OUT) :: errorCode
        call MPI_BCAST(fltVar, 1, MPI_REAL , sourceProc, MPI_COMM_WORLD, errorCode)
    end subroutine

    subroutine syncProcs(errorCode)
        INTEGER(kind = i4), INTENT(OUT) :: errorCode
        call MPI_Barrier(MPI_COMM_WORLD, errorCode)
    end subroutine 

    

end module 
