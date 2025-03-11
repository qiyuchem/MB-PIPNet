program main
use constants
use pes_shell
implicit none
real,dimension(:,:),allocatable::x,tmpx
character,dimension(:),allocatable::symb
real,dimension(:,:),allocatable::gd,gd1
real::eng,eng_ref,box(9),box0(9)
character(len=32)::filename

real::eng1
integer::i,j,k,ierr,ntm
real::rmse,mae,tend,tstart,rmse_gd
integer::cnt,cnt_gd
real::pot(2),vir0(3,3),delta

character(len=20)::date,time
real::t1,t2

call pes_init()

call getarg(1,filename)
open(1,file=trim(filename),status="old")

read(1,*) ntm
rewind(1)

allocate(x(3,ntm))
allocate(tmpx(3,ntm))
allocate(gd(3,ntm),gd1(3,ntm))
allocate(symb(ntm))
box=0.d0

rmse=0.d0
mae=0.d0
rmse_gd=0.d0
cnt=0
cnt_gd=0
do 
  read(1,*,iostat=ierr) ntm
  if (ierr < 0) exit
  read(1,*) eng_ref,box(1),box(5),box(9)
  do i=1,ntm
   read(1,*) symb(i),x(:,i)
  end do
  x = x / auang
  box = box / auang

  call get_pg(x,eng,gd,vir0,box)

  write(*,*) "Potential energy (a.u.)"
  write(*,'(2F15.8)') eng

  write(*,*) "Gradient: Analytical V.S. Numerical, (a.u.)"
  do i=1,3
    do j=1,768
     tmpx=x
     tmpx(i,j)=tmpx(i,j)+0.0005
     call get_pot(tmpx,pot(1),box)
     tmpx=x
     tmpx(i,j)=tmpx(i,j)-0.0005
     call get_pot(tmpx,pot(2),box)
     write(*,'(3F15.8)') gd(i,j),(pot(1)-pot(2))/0.001d0
   end do
  end do

end do
end program 
