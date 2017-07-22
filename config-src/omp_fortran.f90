
program main
    use omp_lib
    implicit none
    integer :: x(10), nt, i

    !$omp parallel
    nt = omp_get_num_threads()
    !$omp end parallel

    !$omp parallel do
    do i = 1, 10
        x(i) = i
    end do
    print '(a,i3,a,i3)', 'openmp x(1)=', x(1), 'nt=', nt
end program
