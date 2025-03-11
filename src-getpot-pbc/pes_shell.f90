module pes_shell
  use constants
  use NN_class
  use bemsa1b
  use bemsa2b
  implicit none
 
    real,parameter::r2i = 7.0/0.5291772083
    real,parameter::r2f = 9.0/0.5291772083

    integer,parameter::max1b=2048
    integer,parameter::max2b=1000000
    real::m1b(max1b,0:4)
    real::p1b(max1b,0:49)
    real::drdx1b(max1b,9,3)
    real::x2b(max2b,7,3),tmp2bgd(max2b,3,12)
    integer::loc(max2b,3),cnt2b
 
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
    subroutine molbasis(geo,xbasis,nin1b,nin2b,n,box)
        implicit none
        real,dimension(:,:),intent(in)::geo
        real,dimension(:,:),intent(inout)::xbasis !include 1b and >=2b bases
        integer::n,nin1b,nin2b !no. of monomers,no.of 1b and >=2b bases
        real::box(9)

        !!!!!!!!
        real,dimension(3,size(geo)/3) :: xn0
        integer::ix,iy,iz,cnt
        integer::i,j,k,l
        integer::i1,j1,k1,l1

        !!!!!!!!
        real,dimension(:,:),allocatable::x1b
        real,dimension(:,:),allocatable::xbasis0
        real,dimension(:),allocatable::xbasis1
        real,dimension(:),allocatable::xbasis2
        integer::n2b,i1loc,j1loc,i2loc,j2loc,i3loc,j3loc
        real::rmax2b(6)
        real::pp1b(nin1b),pp2b(nin2b)
       
        real::ximage1(3,3),ximage2(3,3)
        integer::flag,flag1(3)

        !!!!!!!
        i1loc=nin1b+1
        j1loc=nin1b+nin2b

        xbasis=0.d0

        xn0=geo

        allocate(x1b(3,3))
!!$omp parallel do private(x1b, p1b) shared(xn0, xbasis, n, nin1b)
        do i=1,n
          x1b(1,:)=xn0(:,2*i-1)
          x1b(2,:)=xn0(:,2*i)
          x1b(3,:)=xn0(:,2*n+i)
          call getbasis1b(x1b,pp1b)
          xbasis(1:nin1b,i)=pp1b(1:nin1b)
        end do 
!!$omp end parallel do

        cnt=0
        do i=1,n
          do j=i+1,n
            ximage1(:,1:2)=xn0(:,i*2-1:i*2)
            ximage1(:,3)=xn0(:,2*n+i)
            ximage2(:,1:2) = xn0(:,j*2-1:j*2)
            ximage2(:,3) = xn0(:,2*n+j)
            call setPBC_image(ximage1,ximage2,flag,flag1(:),box)
            rmax2b(1)=sqrt(sum((ximage1(:,3)-ximage2(:,3))**2))
              if(rmax2b(1).le.r2f) then
                cnt=cnt+1
                x2b(cnt,1,:)=ximage1(:,1)
                x2b(cnt,2,:)=ximage1(:,2)
                x2b(cnt,3,:)=ximage2(:,1)
                x2b(cnt,4,:)=ximage2(:,2)
                x2b(cnt,5,:)=ximage1(:,3)
                x2b(cnt,6,:)=ximage2(:,3)
                x2b(cnt,7,1)=rmax2b(1)
                loc(cnt,1)=i
                loc(cnt,2)=j
              end if
          end do
        end do

        allocate(xbasis0(nin2b,cnt))

!$OMP PARALLEL DO 
        do i=1,cnt
           call getbasis2b(x2b(i,:,:),xbasis0(:,i))
        end do
!$OMP END PARALLEL DO

        do i=1,cnt
           xbasis(i1loc:j1loc,loc(i,1))=xbasis(i1loc:j1loc,loc(i,1))+xbasis0(:,i)
           xbasis(i1loc:j1loc,loc(i,2))=xbasis(i1loc:j1loc,loc(i,2))+xbasis0(:,i)
        end do

       return
    end subroutine molbasis


    !=================================================!
    ! Obtain basis function values and derivatives    !
    !=================================================!
    subroutine molbasis_gd(geo,xbasis,nin1b,nin2b,n,box)
        implicit none
        real,dimension(:,:),intent(in)::geo 
        real,dimension(:,:),intent(inout)::xbasis
        integer::n,nin1b,nin2b
        real::box(9)
 
        !!!!!!!!!
        real,dimension(3,size(geo)/3) :: xn0
        integer::ix,iy,iz,cnt
        integer::i,j,k,l
        integer::i1,j1,k1,l1


        !!!!!!!!
        real,dimension(:,:),allocatable::x1b
        real,dimension(:,:),allocatable::xbasis0
        real,dimension(:),allocatable::xbasis1 
        real,dimension(:),allocatable::xbasis2
        integer::n2b,i1loc,j1loc,i2loc,j2loc,i3loc,j3loc
        real::rmax2b(6)
        real::pp1b(nin1b),p2b(nin2b)

        character(len=20)::date,time

        real::ximage1(3,3),ximage2(3,3)
        integer::flag,flag1(3)

        !!!!!!!
        i1loc=nin1b+1
        j1loc=nin1b+nin2b
        xbasis=0.d0

        xn0=geo

        allocate(x1b(3,3))
!!$omp parallel do private(x1b, p1b) shared(xn0, xbasis, n, nin1b)
        do i=1,n
          x1b(1,:)=xn0(:,2*i-1)
          x1b(2,:)=xn0(:,2*i)
          x1b(3,:)=xn0(:,2*n+i)
          call getbasis1b_gd(x1b,pp1b,i)
          xbasis(1:nin1b,i)=pp1b(1:nin1b)
        end do
!!$omp end parallel do

        cnt=0
        do i=1,n
          do j=i+1,n
            ximage1(:,1:2)=xn0(:,i*2-1:i*2)
            ximage1(:,3)=xn0(:,2*n+i)
            ximage2(:,1:2) = xn0(:,j*2-1:j*2)
            ximage2(:,3) = xn0(:,2*n+j)
            call setPBC_image(ximage1,ximage2,flag,flag1(:),box)
            rmax2b(1)=sqrt(sum((ximage1(:,3)-ximage2(:,3))**2))

              if(rmax2b(1).le.r2f) then
                cnt=cnt+1
                x2b(cnt,1,:)=ximage1(:,1)
                x2b(cnt,2,:)=ximage1(:,2)
                x2b(cnt,3,:)=ximage2(:,1)
                x2b(cnt,4,:)=ximage2(:,2)
                x2b(cnt,5,:)=ximage1(:,3)
                x2b(cnt,6,:)=ximage2(:,3)
                x2b(cnt,7,1)=rmax2b(1)
                loc(cnt,1)=i
                loc(cnt,2)=j
              end if
          end do
        end do
        
        cnt2b=cnt
        allocate(xbasis0(nin2b,cnt))

!$OMP PARALLEL DO 
        do i=1,cnt
           call getbasis2b(x2b(i,:,:),xbasis0(:,i))
        end do
!$OMP END PARALLEL DO
        do i=1,cnt
           k=loc(i,1)
           l=loc(i,2)
           xbasis(i1loc:j1loc,k)=xbasis(i1loc:j1loc,k)+xbasis0(:,i)
           xbasis(i1loc:j1loc,l)=xbasis(i1loc:j1loc,l)+xbasis0(:,i)
        end do

 
        return
    end subroutine molbasis_gd

  subroutine get_pot(xx,pot,box)
    real,dimension(:,:),intent(in) :: xx
    real :: pot
    real :: box(9)
    real :: x0(Ninput,Nmol),eng(Nout)
    integer :: i,j
    real,dimension(3,size(xx)/3) :: xn0 

    pot=0.d0
    xn0=xx
    call setPBC(xn0,box)
    call molbasis(xn0,x0,Ninput1b,Ninput2b,Nmol,box)

    do i=1,Ninput
      do j=1,Nmol
       x0(i,j)=x0(i,j)*scale0(i)+scale1(i)
      end do
    end do 
  
!$OMP PARALLEL DO REDUCTION(+:pot) PRIVATE(eng)
    do i=1,Nmol
      call FFNN_output(x0(:,i),eng,wt0,bias0,Neuron0,Activation0)
      pot = pot + eng(1)
    end do
!$OMP END PARALLEL DO

    return
  end subroutine get_pot


  subroutine get_pg(xx,pot,gd,virial,box)
    real,dimension(:,:),intent(in) :: xx
    real :: pot
    real :: box(9)
    real,dimension(:,:),intent(inout) :: gd
    real,dimension(:,:),intent(inout)::virial
    real :: x0(Ninput,Nmol),tmppg(Ninput+1,Nout)!,pg(Ninput,Nout,Nmol)
    integer :: i,j
    real :: a1b,a2b
    real,dimension(:,:,:),allocatable::pg
    real,dimension(3,size(xx)/3) :: xn0

    real,dimension(0:4)::dm
    real,dimension(0:49)::dp
    real::gd1b(3,3)

    character(len=20)::date,time
    integer::iw1,iw2,i1loc,j1loc,i2loc,j2loc
   
    pot=0.d0
    gd=0.d0
    virial=0.d0
    i1loc=1
    j1loc=Ninput1b
    i2loc=1+Ninput1b
    j2loc=Ninput1b+Ninput2b  
    allocate(pg(Ninput,Nout,Nmol))

    xn0=xx
    call setPBC(xn0,box)


    call molbasis_gd(xn0,x0,Ninput1b,Ninput2b,Nmol,box)

    do i=1,Ninput
      do j=1,Nmol
        x0(i,j)=x0(i,j)*scale0(i)+scale1(i)
      end do
    end do

!$OMP PARALLEL DO REDUCTION(+:pot) PRIVATE(tmppg)
    do i=1,Nmol
      call FFNN_pot_gd(x0(:,i),tmppg,wt0,bias0,Neuron0,Activation0)
      pot=pot+tmppg(1,1)
      pg(1:Ninput,1,i)=tmppg(2:Ninput+1,1)*scale0(1:Ninput) 
    end do
!$OMP END PARALLEL DO


    a1b=2.d0
    do i=1,Nmol
      do j=1,3 
         call devmono1b(drdx1b(i,:,:),dm,m1b(i,:),j,a1b)
         call devpoly1b(dm,p1b(i,:),dp)
         gd1b(j,1)=sum(pg(i1loc:j1loc,1,i)*dp(1:49))
      end do
      do j=4,6 
         call devmono1b(drdx1b(i,:,:),dm,m1b(i,:),j,a1b)
         call devpoly1b(dm,p1b(i,:),dp)
         gd1b(j-3,2)=sum(pg(i1loc:j1loc,1,i)*dp(1:49))
      end do
      do j=7,9 
         call devmono1b(drdx1b(i,:,:),dm,m1b(i,:),j,a1b)
         call devpoly1b(dm,p1b(i,:),dp)
         gd1b(j-6,3)=sum(pg(i1loc:j1loc,1,i)*dp(1:49))
      end do
    
      gd(1:3,2*i-1)=gd(1:3,2*i-1)+gd1b(:,1)
      gd(1:3,2*i)=gd(1:3,2*i)+gd1b(:,2)
      gd(1:3,2*Nmol+i)=gd(1:3,2*Nmol+i)+gd1b(:,3)
      virial(1,1)=virial(1,1)+gd1b(1,1)*xn0(1,2*i-1)+gd1b(1,2)*xn0(1,2*i)+gd1b(1,3)*xn0(1,2*Nmol+i)
      virial(1,2)=virial(1,2)+gd1b(2,1)*xn0(1,2*i-1)+gd1b(2,2)*xn0(1,2*i)+gd1b(2,3)*xn0(1,2*Nmol+i)
      virial(1,3)=virial(1,3)+gd1b(3,1)*xn0(1,2*i-1)+gd1b(3,2)*xn0(1,2*i)+gd1b(3,3)*xn0(1,2*Nmol+i)
      virial(2,2)=virial(2,2)+gd1b(2,1)*xn0(2,2*i-1)+gd1b(2,2)*xn0(2,2*i)+gd1b(2,3)*xn0(2,2*Nmol+i)
      virial(2,3)=virial(2,3)+gd1b(3,1)*xn0(2,2*i-1)+gd1b(3,2)*xn0(2,2*i)+gd1b(3,3)*xn0(2,2*Nmol+i)
      virial(3,3)=virial(3,3)+gd1b(3,1)*xn0(3,2*i-1)+gd1b(3,2)*xn0(3,2*i)+gd1b(3,3)*xn0(3,2*Nmol+i)
    end do

!$OMP PARALLEL DO
    do i=1,cnt2b
      call gradient_mol2b(x2b(i,:,:),pg(i2loc:j2loc,1,loc(i,1)),pg(i2loc:j2loc,1,loc(i,2)),tmp2bgd(i,:,:))
    end do
!$OMP END PARALLEL DO

    do i=1,cnt2b
      iw1=loc(i,1)
      iw2=loc(i,2)
      gd(1:3,2*iw1-1)=gd(1:3,2*iw1-1)+tmp2bgd(i,:,1)
      gd(1:3,2*iw1)=gd(1:3,2*iw1)+tmp2bgd(i,:,2)
      gd(1:3,2*iw2-1)=gd(1:3,2*iw2-1)+tmp2bgd(i,:,3)
      gd(1:3,2*iw2)=gd(1:3,2*iw2)+tmp2bgd(i,:,4)
      gd(1:3,2*Nmol+iw1)=gd(1:3,2*Nmol+iw1)+tmp2bgd(i,:,5)
      gd(1:3,2*Nmol+iw2)=gd(1:3,2*Nmol+iw2)+tmp2bgd(i,:,6)

      gd(1:3,2*iw1-1)=gd(1:3,2*iw1-1)+tmp2bgd(i,:,7)
      gd(1:3,2*iw1)=gd(1:3,2*iw1)+tmp2bgd(i,:,8)
      gd(1:3,2*iw2-1)=gd(1:3,2*iw2-1)+tmp2bgd(i,:,9)
      gd(1:3,2*iw2)=gd(1:3,2*iw2)+tmp2bgd(i,:,10)
      gd(1:3,2*Nmol+iw1)=gd(1:3,2*Nmol+iw1)+tmp2bgd(i,:,11)
      gd(1:3,2*Nmol+iw2)=gd(1:3,2*Nmol+iw2)+tmp2bgd(i,:,12)

      virial(1,1)=virial(1,1)+tmp2bgd(i,1,1)*x2b(i,1,1)+tmp2bgd(i,1,2)*x2b(i,2,1)+tmp2bgd(i,1,3)*x2b(i,3,1)+&
                  tmp2bgd(i,1,4)*x2b(i,4,1)+tmp2bgd(i,1,5)*x2b(i,5,1)+tmp2bgd(i,1,6)*x2b(i,6,1)+&
                  tmp2bgd(i,1,7)*x2b(i,1,1)+tmp2bgd(i,1,8)*x2b(i,2,1)+tmp2bgd(i,1,9)*x2b(i,3,1)+&
                  tmp2bgd(i,1,10)*x2b(i,4,1)+tmp2bgd(i,1,11)*x2b(i,5,1)+tmp2bgd(i,1,12)*x2b(i,6,1)
      virial(1,2)=virial(1,2)+tmp2bgd(i,2,1)*x2b(i,1,1)+tmp2bgd(i,2,2)*x2b(i,2,1)+tmp2bgd(i,2,3)*x2b(i,3,1)+&
                  tmp2bgd(i,2,4)*x2b(i,4,1)+tmp2bgd(i,2,5)*x2b(i,5,1)+tmp2bgd(i,2,6)*x2b(i,6,1)+&
                  tmp2bgd(i,2,7)*x2b(i,1,1)+tmp2bgd(i,2,8)*x2b(i,2,1)+tmp2bgd(i,2,9)*x2b(i,3,1)+&
                  tmp2bgd(i,2,10)*x2b(i,4,1)+tmp2bgd(i,2,11)*x2b(i,5,1)+tmp2bgd(i,2,12)*x2b(i,6,1)
      virial(1,3)=virial(1,3)+tmp2bgd(i,3,1)*x2b(i,1,1)+tmp2bgd(i,3,2)*x2b(i,2,1)+tmp2bgd(i,3,3)*x2b(i,3,1)+&
                  tmp2bgd(i,3,4)*x2b(i,4,1)+tmp2bgd(i,3,5)*x2b(i,5,1)+tmp2bgd(i,3,6)*x2b(i,6,1)+&
                  tmp2bgd(i,3,7)*x2b(i,1,1)+tmp2bgd(i,3,8)*x2b(i,2,1)+tmp2bgd(i,3,9)*x2b(i,3,1)+&
                  tmp2bgd(i,3,10)*x2b(i,4,1)+tmp2bgd(i,3,11)*x2b(i,5,1)+tmp2bgd(i,3,12)*x2b(i,6,1)
      virial(2,2)=virial(2,2)+tmp2bgd(i,2,1)*x2b(i,1,2)+tmp2bgd(i,2,2)*x2b(i,2,2)+tmp2bgd(i,2,3)*x2b(i,3,2)+&
                  tmp2bgd(i,2,4)*x2b(i,4,2)+tmp2bgd(i,2,5)*x2b(i,5,2)+tmp2bgd(i,2,6)*x2b(i,6,2)+&
                  tmp2bgd(i,2,7)*x2b(i,1,2)+tmp2bgd(i,2,8)*x2b(i,2,2)+tmp2bgd(i,2,9)*x2b(i,3,2)+&
                  tmp2bgd(i,2,10)*x2b(i,4,2)+tmp2bgd(i,2,11)*x2b(i,5,2)+tmp2bgd(i,2,12)*x2b(i,6,2)
      virial(2,3)=virial(2,3)+tmp2bgd(i,3,1)*x2b(i,1,2)+tmp2bgd(i,3,2)*x2b(i,2,2)+tmp2bgd(i,3,3)*x2b(i,3,2)+&
                  tmp2bgd(i,3,4)*x2b(i,4,2)+tmp2bgd(i,3,5)*x2b(i,5,2)+tmp2bgd(i,3,6)*x2b(i,6,2)+&
                  tmp2bgd(i,3,7)*x2b(i,1,2)+tmp2bgd(i,3,8)*x2b(i,2,2)+tmp2bgd(i,3,9)*x2b(i,3,2)+&
                  tmp2bgd(i,3,10)*x2b(i,4,2)+tmp2bgd(i,3,11)*x2b(i,5,2)+tmp2bgd(i,3,12)*x2b(i,6,2)
      virial(3,3)=virial(3,3)+tmp2bgd(i,3,1)*x2b(i,1,3)+tmp2bgd(i,3,2)*x2b(i,2,3)+tmp2bgd(i,3,3)*x2b(i,3,3)+&
                  tmp2bgd(i,3,4)*x2b(i,4,3)+tmp2bgd(i,3,5)*x2b(i,5,3)+tmp2bgd(i,3,6)*x2b(i,6,3)+&
                  tmp2bgd(i,3,7)*x2b(i,1,3)+tmp2bgd(i,3,8)*x2b(i,2,3)+tmp2bgd(i,3,9)*x2b(i,3,3)+&
                  tmp2bgd(i,3,10)*x2b(i,4,3)+tmp2bgd(i,3,11)*x2b(i,5,3)+tmp2bgd(i,3,12)*x2b(i,6,3)
    end do

    virial(2,1)=virial(1,2)
    virial(3,1)=virial(1,3)
    virial(3,2)=virial(2,3)

    return
  end subroutine get_pg


  subroutine geteng(x,eng,box)
    real,dimension(:),intent(in) ::x
    real::eng,box(9)
    real::xx(3,size(x)/3)
  
    xx=reshape(x,(/3,size(x)/3/))
    call get_pot(xx,eng,box)
    return
  end subroutine geteng

  subroutine getenggd(x,eng,gd,virial,box)
    real,dimension(:),intent(in) ::x
    real::eng,box(9)
    real,dimension(:),intent(out) ::gd
    real,dimension(:,:),intent(inout)::virial
    real::xx(3,size(x)/3),grad(3,size(x)/3)
    integer::i,j

    xx=reshape(x,(/3,size(x)/3/))
    call get_pg(xx,eng,grad,virial,box)
    do i=1,size(x)/3
       gd(3*i-2:3*i)=grad(:,i)
    end do

  end subroutine getenggd

  subroutine hessian(x,H,box)
    real,dimension(:),intent(in)::x
    real,dimension(:,:),intent(inout)::H
    real::box

    real::eps,pot,f_ff,f_fb,f_bf,f_bb,fx
    real,dimension(1:size(x))::tx

    integer::i,j

    eps=0.001d0
    H=0.d0

    tx=x
    call geteng(tx,fx,box)
    do i=1,size(x)
      do j=i+1,size(x)
        tx=x;tx(i)=tx(i)+eps;tx(j)=tx(j)+eps;call geteng(tx,f_ff,box)
        tx=x;tx(i)=tx(i)+eps;tx(j)=tx(j)-eps;call geteng(tx,f_fb,box)
        tx=x;tx(i)=tx(i)-eps;tx(j)=tx(j)+eps;call geteng(tx,f_bf,box)
        tx=x;tx(i)=tx(i)-eps;tx(j)=tx(j)-eps;call geteng(tx,f_bb,box)
        H(i,j)=0.25*(f_ff-f_fb-f_bf+f_bb)/eps/eps
        H(j,i)=H(i,j)
      end do
    end do

    do i=1,size(x)
       tx=x;tx(i)=tx(i)+eps;call geteng(tx,f_ff,box)
       tx=x;tx(i)=tx(i)-eps;call geteng(tx,f_bb,box)
       H(i,i)=(f_ff-2*fx+f_bb)/eps/eps
    end do

    return
  end subroutine hessian

  subroutine getbasis1b(x1b,p1)
    implicit none
    real,dimension(:,:),intent(in)::x1b
    real,dimension(:),intent(inout)::p1
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
    p1=p(1:49)
    return
  end subroutine

  subroutine getbasis1b_gd(x1b,p1,i1b)
    implicit none
    real,dimension(:,:),intent(in)::x1b
    real,dimension(:),intent(inout)::p1
    integer::i1b
    real::xyz(3,3),x(3),a1b
    real::drdx(9,3),dr(3)
    integer::cnt,i,j

    a1b=2.d0
    xyz=x1b

    drdx=0.d0
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

    drdx1b(i1b,:,:)=drdx

    x(1) = dexp(-x(1)/a1b)
    x(2) = dexp(-x(2)/a1b)
    x(3) = dexp(-x(3)/a1b)

    call evmono1b(x,m1b(i1b,:))
    call evpoly1b(m1b(i1b,:),p1b(i1b,:))
    p1=p1b(i1b,1:49)

    return
  end subroutine getbasis1b_gd

  subroutine getbasis2b(x2b,p2)
    implicit none
    real,dimension(:,:),intent(in)::x2b
    real,dimension(:),intent(inout)::p2
    real,dimension(0:1091)::m
    real,dimension(0:139)::p
    real,dimension(1:54)::q
    real::roo,xyz(6,3),x(15)
    real::a2b,dr(3),s
    integer::cnt,i,j

    xyz=x2b(1:6,1:3)
    roo=x2b(7,1)
    a2b=2.5d0

    cnt=1
    do i=1,6
     do j=i+1,6
       dr=xyz(i,:)-xyz(j,:)
       x(cnt)=dsqrt(dot_product(dr,dr))
       x(cnt)=dexp(-x(cnt)/a2b)
       cnt=cnt+1
     end do
    end do

    call f_switch(s,roo,r2i,r2f)
    call evmono2b(x,m)
    call evpoly2b(m,p,q)
    p2=p(1:139)*s
    return
  end subroutine

  subroutine getbasis2b_gd(x2b,m,p,q,s,dsdx,drdx)
    implicit none
    real,dimension(:,:),intent(in)::x2b
    real,dimension(0:1091)::m
    real,dimension(0:139)::p
    real,dimension(1:54)::q
    real::s,dsdx(18),drdx(18,15)

    real::roo,xyz(6,3),x(15)
    real::a2b,dr(3)
    integer::cnt,i,j
    real::tmpx(3),tmps(2),eps,tmproo

    xyz=x2b(1:6,1:3)
    roo=x2b(7,1)
    a2b=2.5d0

    drdx=0.d0
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

    call f_switch(s,roo,r2i,r2f)
    call evmono2b(x,m)
    call evpoly2b(m,p,q)
    
    dsdx=0.d0
    eps=0.001d0
    if((roo.ge.r2i).and.(roo.le.r2f)) then
      do i=1,3
         tmpx(:)=xyz(5,:)
         tmpx(i)=tmpx(i)-eps
         dr=tmpx-xyz(6,:)
         tmproo=dsqrt(dot_product(dr,dr))
         call f_switch(tmps(1),tmproo,r2i,r2f)
         tmpx(:)=xyz(5,:)
         tmpx(i)=tmpx(i)+eps
         dr=tmpx-xyz(6,:)
         tmproo=dsqrt(dot_product(dr,dr))
         call f_switch(tmps(2),tmproo,r2i,r2f)  
         dsdx(12+i)=0.5d0*(tmps(2)-tmps(1))/eps

         tmpx(:)=xyz(6,:)
         tmpx(i)=tmpx(i)-eps
         dr=tmpx-xyz(5,:)
         tmproo=dsqrt(dot_product(dr,dr))
         call f_switch(tmps(1),tmproo,r2i,r2f)
         tmpx(:)=xyz(6,:)
         tmpx(i)=tmpx(i)+eps
         dr=tmpx-xyz(5,:)
         tmproo=dsqrt(dot_product(dr,dr))
         call f_switch(tmps(2),tmproo,r2i,r2f)
         dsdx(15+i)=0.5d0*(tmps(2)-tmps(1))/eps
      end do
    end if

    return
  end subroutine getbasis2b_gd

  subroutine gradient_mol2b(x,pg1,pg2,gd)
    implicit none
    real,dimension(:,:),intent(in)::x
    real,dimension(:),intent(in)::pg1
    real,dimension(:),intent(in)::pg2
    real,dimension(:,:),intent(inout)::gd

    real,dimension(0:139)::dc
    real::xxp(18)
    real,dimension(0:1091)::m
    real,dimension(0:139)::p
    real,dimension(1:54)::q
    real::s,dsdx(18),drdx(18,15)
    real::a2b

    integer::i,j

    a2b=2.5d0
    call getbasis2b_gd(x,m,p,q,s,dsdx,drdx)
    dc=0.d0
    dc(1:139)=pg1(:)
    call derivative_reverse(dc,m,p,q,drdx,xxp,a2b)
    gd(1:3,1)=xxp(1:3)*s
    gd(1:3,2)=xxp(4:6)*s
    gd(1:3,3)=xxp(7:9)*s
    gd(1:3,4)=xxp(10:12)*s
    gd(1:3,5)=xxp(13:15)*s+sum(dc(:)*p(:))*dsdx(13:15)
    gd(1:3,6)=xxp(16:18)*s+sum(dc(:)*p(:))*dsdx(16:18)

    dc=0.d0
    dc(1:139)=pg2(:)
    call derivative_reverse(dc,m,p,q,drdx,xxp,a2b)
    gd(1:3,7)=xxp(1:3)*s
    gd(1:3,8)=xxp(4:6)*s
    gd(1:3,9)=xxp(7:9)*s
    gd(1:3,10)=xxp(10:12)*s
    gd(1:3,11)=xxp(13:15)*s+sum(dc(:)*p(:))*dsdx(13:15)
    gd(1:3,12)=xxp(16:18)*s+sum(dc(:)*p(:))*dsdx(16:18)

    return
  end subroutine gradient_mol2b


      !==================================!!
      !Set PBC for water system
      !==================================!
      subroutine setPBC(x,box)
        real,dimension(:,:),intent(inout) :: x !x(3,natm)
        real :: box(9)
        real :: tempx(3,size(x)/3) !store new coordinates    
        real :: box2(size(box)) !half the box size
        real :: cen_z(3),cen_w(3),di
        integer :: i,j,k,nw

        if (size(box).ne.9) then
          write(*,*) "Current Box size is not acceptable"
          return
        end if

        box2 = 0.5*box

        nw = size(x)/9
        do i = 1,nw
           cen_w = x(:,2*nw+i)
           do j = 1,3
             if(cen_w(j).lt.-box2(4*j-3)) then
               do k = 1,2
               x(j,2*i-2+k) = x(j,2*i-2+k)+box(4*j-3)*nint(-cen_w(j)/box(4*j-3))
               end do
               x(j,2*nw+i) = x(j,2*nw+i) + box(4*j-3)*nint(-cen_w(j)/box(4*j-3))
             else if(cen_w(j).gt.box2(4*j-3)) then
               do k = 1,2
               x(j,2*i-2+k) = x(j,2*i-2+k) -box(4*j-3)*nint(cen_w(j)/box(4*j-3))
               end do
               x(j,2*nw+i) = x(j,2*nw+i) - box(4*j-3)*nint(cen_w(j)/box(4*j-3))
             end if
           end do
           do j = 1,3
             di = x(j,2*i-1) - x(j,2*nw+i)
             if(di.lt.-box2(4*j-3)) then
               x(j,2*i-1) = x(j,2*i-1) + box(4*j-3)
             else if(di.gt.box2(4*j-3)) then
               x(j,2*i-1) = x(j,2*i-1) - box(4*j-3)
             end if
           end do

           do j = 1,3
             di = x(j,2*i) - x(j,2*nw+i)
             if(di.lt.-box2(4*j-3)) then
               x(j,2*i) = x(j,2*i) + box(4*j-3)
             else if(di.gt.box2(4*j-3)) then
               x(j,2*i) = x(j,2*i) - box(4*j-3)
             end if
           end do
        end do

        return
      end subroutine setPBC

      subroutine setPBC_image(x1,x2,flag,flag0,box)
        real,dimension(:,:),intent(in) :: x1 !x(3,3)
        real,dimension(:,:),intent(inout) :: x2 !x(3,3)
        real,dimension(:),intent(inout) :: box
        real :: box2(size(box)) !half the box size
        real :: di
        integer :: i,j,k,flag0(3),flag

        if (size(box).ne.9) then
          write(*,*) "Current Box size is not acceptable"
          return
        end if

        box2 = 0.5*box
        flag = 0

        do i = 1,3
          di = x1(i,3)-x2(i,3)
          if(di.lt.-box2(4*i-3)) then
             x2(i,1) = x2(i,1) - box(4*i-3)
             x2(i,2) = x2(i,2) - box(4*i-3)
             x2(i,3) = x2(i,3) - box(4*i-3)
             flag0(i)=1
          else if(di.gt.box2(4*i-3)) then
             x2(i,1) = x2(i,1) + box(4*i-3)
             x2(i,2) = x2(i,2) + box(4*i-3)
             x2(i,3) = x2(i,3) + box(4*i-3)
             flag0(i)=3
          else
             flag0(i)=2
          end if
        end do
        flag=9*(flag0(1)-1)+3*(flag0(2)-1)+flag0(3)

        return
      end subroutine setPBC_image

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
