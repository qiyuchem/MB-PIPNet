module bemsa1b
  implicit none

contains
    subroutine evmono1b(x,m)
    implicit none
    real,dimension(1:3),intent(in)::x
    real,dimension(0:4),intent(out)::m
    !::::::::::::::::::::

    m(0) = 1.0D0
    m(1) = x(3)
    m(2) = x(2)
    m(3) = x(1)
    m(4) = m(1)*m(2)

    return
  end subroutine evmono1b

  subroutine evpoly1b(m,p)
    implicit none
    real,dimension(0:4),intent(in)::m
    real,dimension(0:49),intent(out)::p
    !::::::::::::::::::::

    p(0) = m(0)
    p(1) = m(1) + m(2)
    p(2) = m(3)
    p(3) = m(4)
    p(4) = p(2)*p(1)
    p(5) = p(1)*p(1) - p(3) - p(3)
    p(6) = p(2)*p(2)
    p(7) = p(2)*p(3)
    p(8) = p(3)*p(1)
    p(9) = p(2)*p(5)
    p(10) = p(2)*p(4)
    p(11) = p(1)*p(5) - p(8)
    p(12) = p(2)*p(6)
    p(13) = p(2)*p(8)
    p(14) = p(2)*p(7)
    p(15) = p(3)*p(3)
    p(16) = p(3)*p(5)
    p(17) = p(2)*p(11)
    p(18) = p(2)*p(9)
    p(19) = p(2)*p(10)
    p(20) = p(1)*p(11) - p(16)
    p(21) = p(2)*p(12)
    p(22) = p(2)*p(15)
    p(23) = p(2)*p(16)
    p(24) = p(2)*p(13)
    p(25) = p(2)*p(14)
    p(26) = p(3)*p(8)
    p(27) = p(3)*p(11)
    p(28) = p(2)*p(20)
    p(29) = p(2)*p(17)
    p(30) = p(2)*p(18)
    p(31) = p(2)*p(19)
    p(32) = p(1)*p(20) - p(27)
    p(33) = p(2)*p(21)
    p(34) = p(2)*p(26)
    p(35) = p(2)*p(27)
    p(36) = p(2)*p(22)
    p(37) = p(2)*p(23)
    p(38) = p(2)*p(24)
    p(39) = p(2)*p(25)
    p(40) = p(3)*p(15)
    p(41) = p(3)*p(16)
    p(42) = p(3)*p(20)
    p(43) = p(2)*p(32)
    p(44) = p(2)*p(28)
    p(45) = p(2)*p(29)
    p(46) = p(2)*p(30)
    p(47) = p(2)*p(31)
    p(48) = p(1)*p(32) - p(42)
    p(49) = p(2)*p(33)

    return
  end subroutine evpoly1b

  subroutine devmono1b(drdx,dm,m,flag,a1)
    implicit none
    real,dimension(9,3),intent(in)::drdx
    real,dimension(0:4),intent(out)::dm
    real,dimension(0:4),intent(in)::m
    integer::flag
    !::::::::::::::::::::
    real::a1

    dm(0) = 0.0D0
    dm(1) = -m(1)/a1*drdx(flag,3)
    dm(2) = -m(2)/a1*drdx(flag,2)
    dm(3) = -m(3)/a1*drdx(flag,1)
    dm(4) = dm(1)*m(2) + m(1)*dm(2)

    return
  end subroutine devmono1b

  subroutine devpoly1b(dm,p,dp)
    implicit none
    real,dimension(0:4),intent(in)::dm
    real,dimension(0:49),intent(in)::p
    real,dimension(0:49),intent(out)::dp
    !::::::::::::::::::::

    dp(0) = dm(0)
    dp(1) = dm(1) + dm(2)
    dp(2) = dm(3)
    dp(3) = dm(4)
    dp(4) = dp(2)*p(1) + p(2)*dp(1)
    dp(5) = dp(1)*p(1) + p(1)*dp(1) - dp(3) - dp(3)
    dp(6) = dp(2)*p(2) + p(2)*dp(2)
    dp(7) = dp(2)*p(3) + p(2)*dp(3)
    dp(8) = dp(3)*p(1) + p(3)*dp(1)
    dp(9) = dp(2)*p(5) + p(2)*dp(5)
    dp(10) = dp(2)*p(4) + p(2)*dp(4)
    dp(11) = dp(1)*p(5) + p(1)*dp(5) - dp(8)
    dp(12) = dp(2)*p(6) + p(2)*dp(6)
    dp(13) = dp(2)*p(8) + p(2)*dp(8)
    dp(14) = dp(2)*p(7) + p(2)*dp(7)
    dp(15) = dp(3)*p(3) + p(3)*dp(3)
    dp(16) = dp(3)*p(5) + p(3)*dp(5)
    dp(17) = dp(2)*p(11) + p(2)*dp(11)
    dp(18) = dp(2)*p(9) + p(2)*dp(9)
    dp(19) = dp(2)*p(10) + p(2)*dp(10)
    dp(20) = dp(1)*p(11) + p(1)*dp(11) - dp(16)
    dp(21) = dp(2)*p(12) + p(2)*dp(12)
    dp(22) = dp(2)*p(15) + p(2)*dp(15)
    dp(23) = dp(2)*p(16) + p(2)*dp(16)
    dp(24) = dp(2)*p(13) + p(2)*dp(13)
    dp(25) = dp(2)*p(14) + p(2)*dp(14)
    dp(26) = dp(3)*p(8) + p(3)*dp(8)
    dp(27) = dp(3)*p(11) + p(3)*dp(11)
    dp(28) = dp(2)*p(20) + p(2)*dp(20)
    dp(29) = dp(2)*p(17) + p(2)*dp(17)
    dp(30) = dp(2)*p(18) + p(2)*dp(18)
    dp(31) = dp(2)*p(19) + p(2)*dp(19)
    dp(32) = dp(1)*p(20) + p(1)*dp(20) - dp(27)
    dp(33) = dp(2)*p(21) + p(2)*dp(21)
    dp(34) = dp(2)*p(26) + p(2)*dp(26)
    dp(35) = dp(2)*p(27) + p(2)*dp(27)
    dp(36) = dp(2)*p(22) + p(2)*dp(22)
    dp(37) = dp(2)*p(23) + p(2)*dp(23)
    dp(38) = dp(2)*p(24) + p(2)*dp(24)
    dp(39) = dp(2)*p(25) + p(2)*dp(25)
    dp(40) = dp(3)*p(15) + p(3)*dp(15)
    dp(41) = dp(3)*p(16) + p(3)*dp(16)
    dp(42) = dp(3)*p(20) + p(3)*dp(20)
    dp(43) = dp(2)*p(32) + p(2)*dp(32)
    dp(44) = dp(2)*p(28) + p(2)*dp(28)
    dp(45) = dp(2)*p(29) + p(2)*dp(29)
    dp(46) = dp(2)*p(30) + p(2)*dp(30)
    dp(47) = dp(2)*p(31) + p(2)*dp(31)
    dp(48) = dp(1)*p(32) + p(1)*dp(32) - dp(42)
    dp(49) = dp(2)*p(33) + p(2)*dp(33)

    return
  end subroutine devpoly1b

end module bemsa1b
