
program main
    use iso_c_binding
    implicit none

interface
    subroutine foo( x ) &
    bind( C, name="foo" )
        use iso_c_binding
        integer(c_int), value :: x
    end subroutine foo
end interface

    integer(c_int), parameter :: x = 100
    print '(a,i3)', 'hello', x
    call foo( x )
end program main
