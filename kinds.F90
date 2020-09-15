module kinds
! !REFDOC:
!  This module is supplied to provide consistent data representation
!  across machine architectures.  It is meant to replace the old
!  Fortran double precision and real *X declarations that were
!  implementation-specific.
!  Users should not need to adjust anything in this module.  If various
!  character strings like long paths to files exceed the default
!  character length, the default value may be increased.

! !USES:
!  uses no other modules

  implicit none
  private
  save

! !DEFINED PARAMETERS:

  integer, parameter, public ::                &
    char_len_short   = 50                     ,&
    char_len         = 256                    ,&
    char_len_long    = 512                    ,&
    logic_kind       = kind(.true.)           ,&
    int_kind         = kind(1)                ,&
    i4               = selected_int_kind(6)   ,&
    i8               = selected_int_kind(13)  ,&
    r4               = selected_real_kind(6)  ,&
    r8               = selected_real_kind(13)


end module kinds


