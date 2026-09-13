module pcg_module
  use iso_c_binding  !imports the standard c-compatible data types c_long_long etc
  implicit none

  interface !translate the code so that fortran and C understands each other
    subroutine pcg_init(seed, seq) bind(C, name="pcg32_srandom")
      !maps fortran subroutine call to C's pcg32_srandom
      import :: c_long_long
      integer(c_long_long), value :: seed
      integer(c_long_long), value :: seq
    end subroutine

    function pcg_rand() bind(C, name="pcg32_random")
      !maps fortran subroutine call to C's pcg32_srandom
      import :: c_int32_t
      integer(c_int32_t) :: pcg_rand
    end function
  end interface

contains

  function pcg_random_real() result(r)
    !Because C's pcg32_random returns a raw 32-bit integer, pcg_random_real converts that integer into a standard double-precision float (real(8)) between $0.0$ and $1.0$.
    use iso_c_binding
    real(8) :: r
    integer(c_int32_t) :: x
    integer(c_int64_t) :: ux

    x = pcg_rand()

    ! Convert to unsigned safely
    ux = int(x, c_int64_t)
    if (ux < 0_c_int64_t) ux = ux + 4294967296_c_int64_t

    r = real(ux, kind=8) / 4294967296.0_8
  end function

end module