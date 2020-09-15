module netCDFio
    implicit none
    
    contains

    subroutine read_PscalarField2D(inputpath,filename,var_name,k_level)

        character(len=*), intent(in) :: inputpath,filename
        character(len=*), intent(in) :: var_name(:)
        integer, intent(in) :: k_level

        character(len=128):: path_n_file
        integer:: num_of_var, var_counter, var_id, file_id, field_id, ierr, endcount(4),startcount(4)
        real(kind =r8) :: field2D(imt,jmt)

        path_n_file = trim(adjustl(inputpath))//'/'//trim(adjustl(filename))

        print *,''
        print *, 'Opening file  ', path_n_file, ' ... '

        ierr = nf_open (trim(path_n_file), nf_nowrite, file_id)
        if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_open')

        print *,''
        print *, 'Opened  ', path_n_file

        num_of_var = size(var_name)

        do var_counter = 1,num_of_var

            print *,''
            print *, 'Trying to 2D field  ', var_name(var_counter), 'at k-level', k_level


            ierr = nf_inq_varid(file_id,trim(var_name(var_counter)),var_id)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_inq_varid')

            print *, 'VAR ID for ', var_name(var_counter), 'is  ', var_id

            startcount = (/1, 1, k_level, 1/)
            endcount = (/ imt, jmt, 1 , 1/)


            ierr = nf_get_vara_real(file_id,var_id,startcount,endcount,field2D)
            if ( ierr /= nf_noerr )  call handle_err(ierr, 'nf_get_field')

            print *, 'Field Obtained'

            inputFields2D(var_counter)%info%fieldName = trim(adjustl(var_name(var_counter)))
            inputFields2D(var_counter)%fieldVal(:,:) = field2D(:,:)

            print *, ''
            print *, var_name(var_counter), 'read SUCCESS'

        end do

        print *,''
        print *, 'Closing file  ', path_n_file, ' ... '

        ierr = nf_close(file_id)
        

    end subroutine
end module netCDFio