module mod_param
    implicit none

    integer, parameter :: rp = selected_real_kind(15, 307) ! double precision

    real(rp), parameter :: pi = acos(-1.0_rp)

end module mod_param