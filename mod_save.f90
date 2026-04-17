module mod_save

    use mod_param, only: rp

    implicit none

contains

    subroutine save(label, n, t, x, phi, phi_ex, nx, L)

        character(len=*), intent(in) :: label
        integer, intent(in) :: n, nx
        real(rp), intent(in) :: t, x(0:nx-1), phi(0:nx-1), phi_ex(0:nx-1), L

        character(len=30) :: fname
        integer :: i

        write(fname, '(A, "_t", I6.6, ".txt")') trim(label), n
        open(unit = 10, file = trim(fname), status = 'replace')

        write(10, '(A, F12.6)') 'Time: ', t
        write(10, '(A)') 'x  ; phi_num ; phi_exact'

        do i = 0, nx-1
            write(10, '(3ES16.7)') x(i), phi(i), phi_ex(i)
        end do

        write(10, '(3ES16.7)') L, phi(0), phi_ex(0)

        close(10)

    end subroutine save

    end module mod_save