module pes_shell
  use constants
  use NN_class
  use bemsa1b
  use bemsa2b
  implicit none
 
contains
  !=================================================!
  ! PES initialization                              !
  !=================================================!
  subroutine pes_init()
    implicit none

    call FFNN_init()
    
    return
  end subroutine pes_init

    !=================================================!
    ! Obtain basis function values                    !
    !=================================================!
    subroutine molbasis(geo,xbasis,nin1b,nin2b,n)
        implicit none
        real,dimension(:,:),intent(in)::geo
        real,dimension(:,:),intent(inout)::xbasis 
        integer::n,nin1b,nin2b,nin3b 

        !!!!!!!!
        real,dimension(3,size(geo)/3) :: xn0
        integer::cnt,i,j,k,l
        real::oodist(n,n)

        !!!!!!!!
        real,dimension(:,:),allocatable::x1b
        real,dimension(:,:,:),allocatable::x2b
        integer,dimension(:,:),allocatable::loc2b
        integer,dimension(:,:),allocatable::loc
        real,dimension(:,:),allocatable::xbasis0
        real,dimension(:),allocatable::xbasis1
        integer::n2b,iloc,jloc
        real::rmax2b(6),rmax,s
        real::p1b(nin1b),p2b(nin2b)
        real::factor
     
        factor=1.d0/3.d0
        xn0=geo

        xbasis=0.d0

        allocate(x1b(3,3))
        do i=1,n
          x1b(1,:)=xn0(:,2*i-1)
          x1b(2,:)=xn0(:,2*i)
          x1b(3,:)=xn0(:,2*n+i)
          call getbasis1b(x1b,p1b)
          xbasis(1:nin1b,i)=p1b(1:nin1b)
        end do

        do i=1,n
          do j=i+1,n
            oodist(i,j)=sqrt(sum((xn0(:,2*n+i)-xn0(:,2*n+j))**2))
          end do
        end do

        n2b=n*(n-1)/2
        allocate(x2b(n2b,7,3))
        allocate(loc2b(n,n))
        allocate(loc(n2b,2))

        cnt=0
        do i=1,n
          do j=i+1,n
            rmax2b(1)=oodist(i,j)
              cnt=cnt+1
              x2b(cnt,1,:)=xn0(:,2*i-1)
              x2b(cnt,2,:)=xn0(:,2*i)
              x2b(cnt,3,:)=xn0(:,2*j-1)
              x2b(cnt,4,:)=xn0(:,2*j)
              x2b(cnt,5,:)=xn0(:,2*n+i)
              x2b(cnt,6,:)=xn0(:,2*n+j)
              x2b(cnt,7,1)=rmax2b(1)
              loc2b(i,j)=cnt
              loc(cnt,1)=i
              loc(cnt,2)=j
          end do
        end do

        allocate(xbasis0(nin2b,cnt))

        !obtain PIP 2b bases
        do i=1,cnt
           call getbasis2b(x2b(i,:,:),p2b)
           xbasis0(:,i)=p2b(:)
        end do

        do i=1,cnt
           iloc=nin1b+1
           jloc=nin1b+nin2b
           xbasis(iloc:jloc,loc(i,1))=xbasis(iloc:jloc,loc(i,1))+xbasis0(:,i)
           xbasis(iloc:jloc,loc(i,2))=xbasis(iloc:jloc,loc(i,2))+xbasis0(:,i)
        end do

       return
    end subroutine molbasis


    !=================================================!
    ! Obtain basis function values and derivatives    !
    !=================================================!
    subroutine molbasis_gd(geo,xbasis,gdbasis,nin1b,nin2b,n)
        implicit none
        real,dimension(:,:),intent(in)::geo 
        real,dimension(:,:),intent(inout)::xbasis 
        real,dimension(:,:,:),intent(inout)::gdbasis 
        integer::n,nin1b,nin2b,nin3b
 
        !!!!!!!!!
        real,dimension(3,size(geo)/3) :: xn0
        integer::cnt,i,j,k,l,ij,ik,jk
        real::oodist(n,n)

        !!!!!!!!
        real,dimension(:,:),allocatable::x1b
        real,dimension(:,:,:),allocatable::x2b
        integer,dimension(:,:),allocatable::loc2b
        integer,dimension(:,:),allocatable::loc
        real,dimension(:,:),allocatable::xbasis0
        real,dimension(:),allocatable::xbasis1 
        real,dimension(:),allocatable::xbasis2
        real,dimension(:,:,:),allocatable::xbasis0_gd
        integer::n2b,iloc,jloc
        real::rmax2b(6),rmax,s
        real::p1b(nin1b),p2b(nin2b),gd_p1b(nin1b,9),gd_p2b(nin2b,18)
        real::factor,tmp1(6),tmp2(6),tmp3(6),tmp4(3),tmp5(3),tmp6(3)

        factor=1.d0/3.d0
        xn0=geo
 
        xbasis=0.d0
        gdbasis=0.d0

        allocate(x1b(3,3))
        !obtain PIP 1b bases, this part can be parallelized
        do i=1,n
          x1b(1,:)=xn0(:,2*i-1)
          x1b(2,:)=xn0(:,2*i)
          x1b(3,:)=xn0(:,2*n+i)
          call getbasis1b_gd(x1b,p1b,gd_p1b)
          xbasis(1:nin1b,i)=p1b(1:nin1b)
          gdbasis(1:nin1b,i,6*i-5:6*i-3)=gd_p1b(1:nin1b,1:3)
          gdbasis(1:nin1b,i,6*i-2:6*i)=gd_p1b(1:nin1b,4:6)
          gdbasis(1:nin1b,i,6*n+3*i-2:6*n+3*i)=gd_p1b(1:nin1b,7:9)
        end do

        !calculate all oo distances 
        do i=1,n
          do j=i+1,n
            oodist(i,j)=sqrt(sum((xn0(:,2*n+i)-xn0(:,2*n+j))**2))
          end do
        end do

        n2b=n*(n-1)/2
        allocate(x2b(n2b,7,3))
        allocate(loc2b(n,n))
        allocate(loc(n2b,2))

        !prepare 2b geometries     
        cnt=0
        do i=1,n
          do j=i+1,n
            rmax2b(1)=oodist(i,j)
              cnt=cnt+1
              x2b(cnt,1,:)=xn0(:,2*i-1)
              x2b(cnt,2,:)=xn0(:,2*i)
              x2b(cnt,3,:)=xn0(:,2*j-1)
              x2b(cnt,4,:)=xn0(:,2*j)
              x2b(cnt,5,:)=xn0(:,2*n+i)
              x2b(cnt,6,:)=xn0(:,2*n+j)
              x2b(cnt,7,1)=rmax2b(1)
              loc2b(i,j)=cnt
              loc(cnt,1)=i
              loc(cnt,2)=j
          end do
        end do

        allocate(xbasis0(nin2b,cnt))
        allocate(xbasis0_gd(nin2b,18,cnt))

        !obtain PIP 2b bases
        do i=1,cnt
           call getbasis2b_gd(x2b(i,:,:),p2b,gd_p2b)
           xbasis0(:,i)=p2b(:)
           xbasis0_gd(:,:,i)=gd_p2b
        end do

        do i=1,cnt
           iloc=nin1b+1
           jloc=nin1b+nin2b
           k=loc(i,1)
           l=loc(i,2)
           xbasis(iloc:jloc,k)=xbasis(iloc:jloc,k)+xbasis0(:,i)
           xbasis(iloc:jloc,l)=xbasis(iloc:jloc,l)+xbasis0(:,i)
           gdbasis(iloc:jloc,k,6*k-5:6*k)=gdbasis(iloc:jloc,k,6*k-5:6*k)+xbasis0_gd(1:nin2b,1:6,i)
           gdbasis(iloc:jloc,k,6*l-5:6*l)=gdbasis(iloc:jloc,k,6*l-5:6*l)+xbasis0_gd(1:nin2b,7:12,i)
           gdbasis(iloc:jloc,k,6*n+3*k-2:6*n+3*k)=gdbasis(iloc:jloc,k,6*n+3*k-2:6*n+3*k)+xbasis0_gd(1:nin2b,13:15,i)
           gdbasis(iloc:jloc,k,6*n+3*l-2:6*n+3*l)=gdbasis(iloc:jloc,k,6*n+3*l-2:6*n+3*l)+xbasis0_gd(1:nin2b,16:18,i)
           gdbasis(iloc:jloc,l,6*k-5:6*k)=gdbasis(iloc:jloc,l,6*k-5:6*k)+xbasis0_gd(1:nin2b,1:6,i)
           gdbasis(iloc:jloc,l,6*l-5:6*l)=gdbasis(iloc:jloc,l,6*l-5:6*l)+xbasis0_gd(1:nin2b,7:12,i)
           gdbasis(iloc:jloc,l,6*n+3*k-2:6*n+3*k)=gdbasis(iloc:jloc,l,6*n+3*k-2:6*n+3*k)+xbasis0_gd(1:nin2b,13:15,i)
           gdbasis(iloc:jloc,l,6*n+3*l-2:6*n+3*l)=gdbasis(iloc:jloc,l,6*n+3*l-2:6*n+3*l)+xbasis0_gd(1:nin2b,16:18,i)
        end do


        return
    end subroutine molbasis_gd

  subroutine get_pot(xx,pot)
    real,dimension(:,:),intent(in) :: xx
    real :: pot
    real :: x0(Ninput,Nmol),eng(Nout)
    integer :: i,j 

    pot=0.d0
   
    call molbasis(xx,x0,Ninput1b,Ninput2b,Nmol)

    do i=1,Ninput
      do j=1,Nmol
       x0(i,j)=2.d0*(x0(i,j)-xmin(i))/(xmax(i)-xmin(i))-1.d0
      end do
    end do

    do i=1,Nmol
      call FFNN_output(x0(:,i),eng,wt0,bias0,Neuron0,Activation0)
      pot=pot+eng(1)
    end do
    pot=(pot+1.d0)*(ymax-ymin)/2.d0+ymin
    return
  end subroutine get_pot


  subroutine get_pg(xx,pot,gd)
    real,dimension(:,:),intent(in) :: xx
    real :: pot
    real,dimension(:,:),intent(inout) :: gd
    real :: x0(Ninput,Nmol),gdx0(Ninput,Nmol,size(xx)),pg(Ninput+1,Nout)
    real :: eng(Nout)
    integer :: i,j,k,cnt
   
    pot=0.d0
    gd=0.d0
    call molbasis_gd(xx,x0,gdx0,Ninput1b,Ninput2b,Nmol)

    do i=1,Ninput
      do j=1,Nmol
        x0(i,j)=2.d0*(x0(i,j)-xmin(i))/(xmax(i)-xmin(i))-1.d0
        gdx0(i,j,:)=gdx0(i,j,:)*2.d0/(xmax(i)-xmin(i))
      end do
    end do

    do i=1,Nmol
      call FFNN_pot_gd(x0(:,i),pg,wt0,bias0,Neuron0,Activation0)
      pot=pot+pg(1,1)
      
      cnt=1
      do j=1,size(xx,2)
        do k=1,size(xx,1)
           gd(k,j)=gd(k,j)+dot_product(pg(2:Ninput+1,1),gdx0(1:Ninput,i,cnt))*(ymax-ymin)/2.d0
           cnt=cnt+1
        end do
      end do
    end do
    pot=(pot+1.d0)*(ymax-ymin)/2.d0+ymin

    return
  end subroutine get_pg


  subroutine geteng(x,eng)
    real,dimension(:),intent(in) ::x
    real::eng
    real::xx(3,size(x)/3)
  
    xx=reshape(x,(/3,size(x)/3/))
    call get_pot(xx,eng)
    return
  end subroutine geteng

  subroutine getenggd(x,eng,gd)
    real,dimension(:),intent(in) ::x
    real::eng
    real,dimension(:),intent(out) ::gd
    real::xx(3,size(x)/3),grad(3,size(x)/3)
    integer::i,j

    xx=reshape(x,(/3,size(x)/3/))
    call get_pg(xx,eng,grad)
    do i=1,size(x)/3
       gd(3*i-2:3*i)=grad(:,i)
    end do

    return
  end subroutine getenggd

  subroutine hessian(x,H)
    real,dimension(:),intent(in)::x
    real,dimension(:,:),intent(inout)::H

    real::eps,pot,f_ff,f_fb,f_bf,f_bb,fx
    real,dimension(1:size(x))::tx

    integer::i,j

    eps=0.001d0
    H=0.d0

    tx=x
    call geteng(tx,fx)
    do i=1,size(x)
      do j=i+1,size(x)
        tx=x;tx(i)=tx(i)+eps;tx(j)=tx(j)+eps;call geteng(tx,f_ff)
        tx=x;tx(i)=tx(i)+eps;tx(j)=tx(j)-eps;call geteng(tx,f_fb)
        tx=x;tx(i)=tx(i)-eps;tx(j)=tx(j)+eps;call geteng(tx,f_bf)
        tx=x;tx(i)=tx(i)-eps;tx(j)=tx(j)-eps;call geteng(tx,f_bb)
        H(i,j)=0.25*(f_ff-f_fb-f_bf+f_bb)/eps/eps
        H(j,i)=H(i,j)
      end do
    end do

    do i=1,size(x)
       tx=x;tx(i)=tx(i)+eps;call geteng(tx,f_ff)
       tx=x;tx(i)=tx(i)-eps;call geteng(tx,f_bb)
       H(i,i)=(f_ff-2*fx+f_bb)/eps/eps
    end do

    return
  end subroutine hessian

  subroutine getbasis1b(x1b,p1b)
    implicit none
    real,dimension(:,:),intent(in)::x1b
    real,dimension(:),intent(inout)::p1b
    real,dimension(0:4)::m
    real,dimension(0:49)::p
    real::xyz(3,3),x(3),a1b

    a1b=2.d0
    xyz=x1b
    x(1)=sqrt((xyz(1,1)-xyz(2,1))**2+ &
             (xyz(1,2)-xyz(2,2))**2+ &
             (xyz(1,3)-xyz(2,3))**2)
    x(2)=sqrt((xyz(1,1)-xyz(3,1))**2+ &
             (xyz(1,2)-xyz(3,2))**2+ &
             (xyz(1,3)-xyz(3,3))**2)
    x(3)=sqrt((xyz(2,1)-xyz(3,1))**2+ &
             (xyz(2,2)-xyz(3,2))**2+ &
             (xyz(2,3)-xyz(3,3))**2)
    x(1) = dexp(-x(1)/a1b)
    x(2) = dexp(-x(2)/a1b)
    x(3) = dexp(-x(3)/a1b)

    call evmono1b(x,m)
    call evpoly1b(m,p)
    p1b=p(1:49)
    return
  end subroutine

  subroutine getbasis1b_gd(x1b,p1b,gd_p1b)
    implicit none
    real,dimension(:,:),intent(in)::x1b
    real,dimension(:),intent(inout)::p1b
    real,dimension(:,:),intent(inout)::gd_p1b
    real,dimension(0:4)::m
    real,dimension(0:4)::dm
    real,dimension(0:49)::p
    real,dimension(0:49)::dp
    real::xyz(3,3),x(3),a1b
    real::drdx(9,3),dr(3)
    integer::cnt,i,j

    a1b=2.d0
    xyz=x1b

    cnt=1
    do i=1,3
     do j=i+1,3
       dr=xyz(i,:)-xyz(j,:)
       x(cnt)=dsqrt(dot_product(dr,dr))
       drdx(3*i-2:3*i,cnt)=dr/x(cnt)
       drdx(3*j-2:3*j,cnt)=-drdx(3*i-2:3*i,cnt)
       cnt=cnt+1
     end do
    end do

    x(1) = dexp(-x(1)/a1b)
    x(2) = dexp(-x(2)/a1b)
    x(3) = dexp(-x(3)/a1b)

    call evmono1b(x,m)
    call evpoly1b(m,p)
    p1b=p(1:49)

    !get derivatives of 1b bases
    do i=1,9
      call devmono1b(drdx,dm,m,i,a1b)
      call devpoly1b(dm,p,dp)
      gd_p1b(:,i)=dp(1:49)
    end do

    return
  end subroutine getbasis1b_gd

  subroutine getbasis2b(x2b,p2b)
    implicit none
    real,dimension(:,:),intent(in)::x2b
    real,dimension(:),intent(inout)::p2b
    real,dimension(0:1091)::m
    real,dimension(0:139)::p
    real,dimension(1:54)::q
    real::roo,xyz(6,3),x(15)
    real::a2b,dr(3),s
    integer::cnt,i,j

    xyz=x2b(1:6,1:3)
    roo=x2b(7,1)
    a2b=3.0d0

    cnt=1
    do i=1,6
     do j=i+1,6
       dr=xyz(i,:)-xyz(j,:)
       x(cnt)=dsqrt(dot_product(dr,dr))
       x(cnt)=dexp(-x(cnt)/a2b)
       cnt=cnt+1
     end do
    end do

    call evmono2b(x,m)
    call evpoly2b(m,p,q)
    p2b=p(1:139)
    return
  end subroutine

  subroutine getbasis2b_gd(x2b,p2b,gd_p2b)
    implicit none
    real,dimension(:,:),intent(in)::x2b
    real,dimension(:),intent(inout)::p2b
    real,dimension(:,:),intent(inout)::gd_p2b
    real,dimension(0:1091)::m
    real,dimension(0:1091)::dm
    real,dimension(0:139)::p
    real,dimension(0:139)::dp
    real,dimension(1:54)::q
    real,dimension(1:54)::dq
    real::roo,xyz(6,3),x(15)
    real::a2b,dr(3),drdx(18,15),s
    integer::cnt,i,j

    xyz=x2b(1:6,1:3)
    roo=x2b(7,1)

    a2b=3.0d0

    cnt=1
    do i=1,6
     do j=i+1,6
       dr=xyz(i,:)-xyz(j,:)
       x(cnt)=dsqrt(dot_product(dr,dr))
       drdx(3*i-2:3*i,cnt)=dr/x(cnt)
       drdx(3*j-2:3*j,cnt)=-drdx(3*i-2:3*i,cnt)
       x(cnt)=dexp(-x(cnt)/a2b)
       cnt=cnt+1
     end do
    end do

    call evmono2b(x,m)
    call evpoly2b(m,p,q)
    p2b=p(1:139)
    
    do i=1,18
      call devmono2b(drdx,dm,m,i,a2b)
      call devpoly2b(dm,p,q,dp)
      gd_p2b(:,i)=dp(1:139)
    end do
    return
  end subroutine getbasis2b_gd

  !==================================================!
  ! switching functions for 2b,3b,4b
  !==================================================!
  subroutine f_switch(s,r,ri,rf)
    real,intent(out)::s
    real,intent(in)::r,ri,rf
    !::::::::::::::::::::
    real::ra,ra2,ra3

    if(r.le.ri) then
        s = 1.0
    else if(r.le.rf) then
        ra=(r-ri)/(rf-ri)
        ra2=ra*ra
        ra3=ra2*ra
        s = 1.0-(10.0*ra3-15.0*ra*ra3+6.0*ra3*ra2)
    else
        s = 0.0
    end if
    return
  end subroutine f_switch

end module pes_shell
