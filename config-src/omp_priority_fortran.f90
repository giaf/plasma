
program main
    use omp_lib
    implicit none
    integer, parameter :: n = 1000
    integer :: x(n), i, expect
    
    !$omp parallel
        !! todo: verify priority syntax
        !$omp task depend(inout:x(1:n)), priority(1)
        call task( n, x, 0 );
        !$omp end task

        !$omp task depend(inout:x(1:n)), priority(2)
        call task( n, x, 100 );
        !$omp end task
    !$omp end parallel
    
    do i = 1, n
        expect = 100 + i
        if (x(i) .ne. expect) then
            print '(a,i4,a,i4,a,i4)', 'openmp task priority failed, x(', &
                  i, ') = ', x(i), ', expected ', expect
            stop 1
        endif
    end do
    print '(a)', 'openmp task priority ok'
end program

subroutine task( n, x, id )
    integer :: n, x(n), id
    integer :: i
    do i = 1, n
        x(i) = id + i
    end do
end subroutine
