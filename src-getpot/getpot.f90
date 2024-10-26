program main
use constants
use pes_shell
implicit none
real::x(3,9),eng
integer::i,j,k,ierr,ntm
character::symb(9)

real::gd(3,9),tmpx(3,9),pot(2)
call pes_init()

open(1,file="trimer.xyz")

do 
  read(1,*,iostat=ierr) ntm
  if (ierr < 0) exit
  read(1,*) 
  do i=1,ntm
   read(1,*) symb(i),x(:,i)
  end do
  x = x / auang

  call get_pg(x,eng,gd)
  write(*,*) "Potential energy (a.u.)"
  write(*,'(2F15.8)') eng
  write(*,*) "Gradient: Analytical V.S. Numerical, (a.u.)"
  do i=1,3
    do j=1,9
     tmpx=x
     tmpx(i,j)=tmpx(i,j)+0.0005
     call get_pot(tmpx,pot(1))
     tmpx=x
     tmpx(i,j)=tmpx(i,j)-0.0005
     call get_pot(tmpx,pot(2))
     write(*,'(3F15.8)') gd(i,j),(pot(1)-pot(2))/0.001d0
   end do
  end do

end do

end program 
