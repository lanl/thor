module comm_types
  implicit none
  contains
    function get_mpi_op(oper) result (op)
      USE mpi_f08
        character(len=*), intent(in)  :: oper
        TYPE(MPI_OP)                  :: op
        SELECT CASE (oper)
            CASE ('null')
                op = MPI_OP_NULL
            CASE ('sum')
                op = MPI_SUM
            CASE ('min')
                op = MPI_MIN
            CASE ('max')
                op = MPI_MAX
            CASE ('prod')
                op = MPI_PROD
            CASE ('land')
                op = MPI_LAND
            CASE ('band')
                op = MPI_BAND
            CASE ('lor')
                op = MPI_LOR
            CASE ('bor')
                op = MPI_BOR
            CASE ('lxor')
                op = MPI_LXOR
            CASE ('bxor')
                op = MPI_BXOR
            CASE ('arg_min')
                op = MPI_MINLOC
            CASE ('arg_max')
                op = MPI_MAXLOC
            CASE ('min_loc')
                op = MPI_MINLOC
            CASE ('max_loc')
                op = MPI_MAXLOC
            CASE ('replace')
                op = MPI_REPLACE
            CASE DEFAULT
                op = MPI_SUM
        END SELECT
    end function get_mpi_op
    function get_mpi_type(T) result (tp)
        USE mpi_f08
        character(len=*), intent(in)  :: T
        TYPE(MPI_Datatype)            :: tp
        SELECT CASE (T)
            CASE ('MPI_i1')
                tp = MPI_INTEGER1
            CASE ('MPI_i2')
                tp = MPI_INTEGER2
            CASE ('MPI_i4')
                tp = MPI_INTEGER4
            CASE ('MPI_i8')
                tp = MPI_INTEGER8
            CASE ('MPI_i16')
                tp = MPI_INTEGER16
            CASE ('MPI_r2')
                tp = MPI_REAL2
            CASE ('MPI_r4')
                tp = MPI_REAL4
            CASE ('MPI_r8')
                tp = MPI_REAL8
            CASE ('MPI_r16')
                tp = MPI_REAL16
        END SELECT
    end function get_mpi_type

end module comm_types
