
program main
    use omp_lib
    implicit none
    integer, parameter :: n = 1000, last = 1000
    integer :: x(n), iter, i, expect

    do iter = 1, 100
        !! inserts last/10 tasks that update x
        !$omp parallel
        do i = 0, last, 10
            !$omp task depend(inout:x(1:n))
            call task( n, x, i )
            !$omp end task
        end do
        !$omp end parallel
        
        !! verify that updates worked
        do i = 1, n
            expect = last + i
            if (x(i) .ne. expect) then
                print '(a,i4,a,i4)', 'openmp task depend failed, x(', i, ') = ', x(i)
                stop 1
            endif
        end do
    end do
    print '(a)', 'openmp task depend seems ok'
end program
    
subroutine task( n, x, id )
    integer :: n, x(n), id
    integer :: i
    do i = 1, n
        x(i) = id + i
    end do
end subroutine task
