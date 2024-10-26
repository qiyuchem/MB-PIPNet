!======================================
!Module for Feedforward Neural Network
!Qi Yu 
!======================================

module NN_class
    implicit none
    !NN structure information
    character(len=100) :: NNstructure 
    integer :: Natm, Nmol, Ninput1b, Ninput2b, Ninput3b
    integer :: Ninput, Nout, Nlayer, maxNeuron, Nvar
    integer,dimension(:),allocatable :: Neuron0 
    integer,dimension(:),allocatable :: Activation0
    
    !NN weights and bias
    real,dimension(:,:,:),allocatable :: wt0 
    real,dimension(:,:),allocatable :: bias0 
    real,allocatable :: tmpwt(:,:,:),tmpbias(:,:)

    real,dimension(:),allocatable :: xmin,xmax
    real:: ymax,ymin

  contains

  !==================================================
  !Read and initialize the structure of FFNN
  !==================================================
  subroutine FFNN_init()
    implicit none
    integer :: file,ierr,i,j
    character(len=99)::filecoef

    !Read basic structure of NN
    call get_file(file)
    open(file,status='old',file='NN_struct.dat')

    read(file,*,iostat=ierr) Natm, Nmol, Ninput1b, Ninput2b

    read(file,*,iostat=ierr) Nlayer
    if(ierr.ne.0) stop "Error in reading NN structure"
    
    allocate(Neuron0(1:Nlayer+1))
    allocate(Activation0(1:Nlayer))

    read(file,*,iostat=ierr) Neuron0(1:Nlayer+1)
    if(ierr.ne.0) stop "Error in reading NN structure"
    
    Ninput = Neuron0(1)
    Nout = Neuron0(Nlayer+1)
    maxNeuron = maxval(Neuron0(1:Nlayer+1))

    if(Ninput.ne.(Ninput1b+Ninput2b)) stop "Error in no. of inputs"   

    !define number of weights 
    Nvar = 0
    do i = 1, Nlayer
      Nvar=Nvar+Neuron0(i+1)*Neuron0(i)+Neuron0(i+1)
      read(file,*,iostat=ierr) Activation0(i)
      if(ierr.ne.0) stop "Error in reading NN structure"
    end do
    close(file)

    !define FFNN parameters
    allocate(wt0(maxNeuron,maxNeuron,Nlayer))
    allocate(bias0(maxNeuron,Nlayer))
 
    allocate(tmpwt(maxNeuron,maxNeuron,Nlayer))
    allocate(tmpbias(maxNeuron,Nlayer))

    allocate(xmin(Ninput),xmax(Ninput))

    filecoef="WB01.txt"
    call init_wt_bias(-1.d0,1.d0,1,filecoef)    
    return
  end subroutine FFNN_init
  
  !==================================================
  !Deallocation the structure of FFNN
  !==================================================
  subroutine FFNN_clean()
    implicit none
    
    deallocate(Neuron0)
    deallocate(wt0)
    deallocate(bias0)
    deallocate(tmpwt)
    deallocate(tmpbias)
    deallocate(Activation0)
    return
  end subroutine FFNN_clean

  !==================================================
  !Initialize/copy/record/update the weight and bias
  !==================================================
  subroutine init_wt_bias(lower, upper, flag, filewt)
    implicit none
    real,intent(in) :: lower, upper
    character(len=99),intent(in),optional :: filewt
    integer::flag

    integer::cnt,iNeuron,i,j,k,iun

    if(flag==0) then
      call random_number(wt0)
      call random_number(bias0)
      wt0=wt0*(upper-lower)+lower
      bias0=bias0*(upper-lower)+lower
    end if

    if(flag==1) then
      call get_file(iun)
      open(iun,file=trim(filewt))
      read(iun,*)
      do i=1,Ninput
        read(iun,*) xmin(i),xmax(i) 
      end do
      read(iun,*)
      read(iun,*) ymin,ymax
      read(iun,*)
      cnt=0
      do i=1,Nlayer
        iNeuron=Neuron0(i+1)
        do j=1,Neuron0(i)
           do k=1,iNeuron
              read(iun,*) wt0(k,j,i)
              cnt=cnt+1
           end do
        end do
        do j=1,iNeuron
           read(iun,*) bias0(j,i)
           cnt=cnt+1
        end do
      end do
      if(cnt .ne. Nvar) then
        write(*,*) "Wrong number of parameters"
      end if
      close(iun)
     end if

    return
  end subroutine init_wt_bias

  subroutine copy_wt_bias()
    implicit none
    tmpwt = wt0
    tmpbias = bias0

    return
  end subroutine copy_wt_bias

  subroutine rec_wt_bias()
    implicit none
    wt0 = tmpwt
    bias0 = tmpbias

    return
  end subroutine rec_wt_bias

  subroutine update_wt_bias(dwb)
    implicit none
    real,dimension(:),intent(in)::dwb
    integer::i,j,cnt,Nbias,iNeuron

    cnt=0
    do i = 1,Nlayer
      iNeuron=Neuron0(i+1)
      do j = 1,Neuron0(i)
        wt0(1:iNeuron,j,i)=wt0(1:iNeuron,j,i)-dwb(cnt+1:cnt+iNeuron)
        cnt=cnt+iNeuron
      end do

      bias0(1:iNeuron,i)=bias0(1:iNeuron,i)-dwb(cnt+1:cnt+iNeuron)
      cnt=cnt+iNeuron
    end do

    return
  end subroutine update_wt_bias

  !==================================================
  !Calculate the output of FFNN
  !==================================================
  subroutine FFNN_output(x,pot,wt,bias,Neuron,Activation)
    implicit none
    real,intent(in) :: x(Ninput)
    real,dimension(:),intent(inout) :: pot
    real,dimension(:,:,:),intent(in) :: wt
    real,dimension(:,:),intent(in) :: bias
    integer,dimension(:),intent(in) :: Neuron
    integer,dimension(:),intent(in) :: Activation

    real :: input1(1:maxNeuron,0:Nlayer)
    real :: input0(1:maxNeuron,1:Nlayer)

    integer :: i,j,iNeuron,jNeuron

    input1(1:Ninput,0) = x(1:Ninput)
    do i=1,Nlayer
       iNeuron=Neuron(i+1)
       jNeuron=Neuron(i)
       input0(1:iNeuron,i)=bias(1:iNeuron,i)
       call dgemv('N', iNeuron, jNeuron, 1.d0, &
                wt(1:iNeuron,1:jNeuron,i), iNeuron, &
                input1(1:jNeuron,i-1), 1, 1.d0, input0(1:iNeuron,i),1)
       call factivate(input0(1:iNeuron,i),input1(1:iNeuron,i), Activation(i))
    end do

    pot(1:Nout)=input1(1:Nout,Nlayer)

    return
  end subroutine FFNN_output

  !====================================================
  !Calculate the output of FFNN and gradient over input
  !====================================================
  subroutine FFNN_pot_gd(x,pg,wt,bias,Neuron,Activation)
    implicit none
    real,intent(in) :: x(Ninput)
    real,dimension(:,:),intent(inout) :: pg 
    real,dimension(:,:,:),intent(in) :: wt
    real,dimension(:,:),intent(in) :: bias
    integer,dimension(:),intent(in) :: Neuron
    integer,dimension(:),intent(in) :: Activation

    real :: input1(1:maxNeuron,0:Nlayer)
    real :: input0(1:maxNeuron,1:Nlayer)
    real :: dinput1dx(1:Ninput,1:maxNeuron,0:Nlayer)
    real :: dinput0dx(1:Ninput,1:maxNeuron,1:Nlayer)
    real :: wk1(1:maxNeuron)
    integer :: i,j,iNeuron,jNeuron
    
    input1(1:Ninput,0)=x(1:Ninput)
    dinput1dx(1:Ninput,1:Ninput,0)=0.d0
    do i=1,Ninput
       dinput1dx(i,i,0)=1.d0
    end do

    do i=1,Nlayer
       iNeuron=Neuron(i+1)
       jNeuron=Neuron(i)

       input0(1:iNeuron,i)=bias(1:iNeuron,i)
       call dgemv('N', iNeuron, jNeuron, 1.d0, &
                wt(1:iNeuron,1:jNeuron,i), iNeuron, &
                input1(1:jNeuron,i-1), 1, 1.d0, input0(1:iNeuron,i),1)
       call factivate(input0(1:iNeuron,i),input1(1:iNeuron,i), Activation(i))

       call dgemm('N','T',Ninput,iNeuron,jNeuron,1.d0,dinput1dx(1:Ninput,1:jNeuron,i-1),&
            Ninput,wt(1:iNeuron,1:jNeuron,i),iNeuron,0.d0,dinput0dx(1:Ninput,1:iNeuron,i),Ninput)
       call dfactivate(input0(1:iNeuron,i),wk1(1:iNeuron),Activation(i))
       do j=1,Ninput
          dinput1dx(j,1:iNeuron,i)=dinput0dx(j,1:iNeuron,i)*wk1(1:iNeuron)
       end do
    end do

    pg(1,1:Nout)=input1(1:Nout,Nlayer)
    pg(2:Ninput+1,1:Nout)=dinput1dx(1:Ninput,1:Nout,Nlayer)

    return
  end subroutine FFNN_pot_gd

  !====================================================
  !Calculate the FFNN derivative over parameters
  !====================================================
  subroutine FFNN_wb(x,dydwb,wt,bias,Neuron,Activation)
    implicit none
    real,intent(in) :: x(Ninput)
    real,dimension(:,:),intent(inout) :: dydwb 
    real,dimension(:,:,:),intent(in) :: wt
    real,dimension(:,:),intent(in) :: bias
    integer,dimension(:),intent(in) :: Neuron
    integer,dimension(:),intent(in) :: Activation

    real :: input1(1:maxNeuron,0:Nlayer)
    real :: input0(1:maxNeuron,1:Nlayer)
    real :: dydinput0(maxNeuron,Nout,Nlayer)

    integer :: i,j,k,iNeuron,jNeuron,cnt
    real :: df(1:maxNeuron)


    input1(1:Ninput,0)=x(1:Ninput)  
    do i=1,Nlayer
       iNeuron=Neuron(i+1)
       jNeuron=Neuron(i)

       input0(1:iNeuron,i)=bias(1:iNeuron,i)
       
       call dgemv('N', iNeuron, jNeuron, 1.d0, &
                wt(1:iNeuron,1:jNeuron,i), iNeuron, &
                input1(1:jNeuron,i-1), 1, 1.d0, input0(1:iNeuron,i),1)
       call factivate(input0(1:iNeuron,i),input1(1:iNeuron,i), Activation(i))
    end do


    dydinput0(1:Nout,1:Nout,Nlayer)=0.d0
    do i=1,Nout
       call dfactivate(input0(i:i,Nlayer),dydinput0(i:i,i,Nlayer),Activation(Nlayer))
    end do
 
    do i=Nlayer,2,-1
       iNeuron=Neuron(i)
       jNeuron=Neuron(i+1)
       call dgemm('T','N',iNeuron,Nout,jNeuron,1.d0,wt(1:jNeuron,1:iNeuron,i),jNeuron,&
                         dydinput0(1:jNeuron,1:Nout,i),jNeuron,0.d0,dydinput0(1:iNeuron,1:Nout,i-1),iNeuron)

       call dfactivate(input0(1:iNeuron,i-1), df(1:iNeuron),Activation(i-1))

       do j=1,Nout
          dydinput0(1:iNeuron,j,i-1)=dydinput0(1:iNeuron,j,i-1)*df(1:iNeuron)
       end do
    end do

    do i=1,Nout
       cnt = 0
       do j=1,Nlayer
         iNeuron=Neuron(j+1)
         do k=1,Neuron(j)
           dydwb(cnt+1:cnt+iNeuron,i)=dydinput0(1:iNeuron,i,j)*input1(k,j-1)
           cnt=cnt+iNeuron
         end do
         dydwb(cnt+1:cnt+iNeuron,i)=dydinput0(1:iNeuron,i,j)
         cnt=cnt+iNeuron
       end do
       if(cnt.ne.Nvar) then
         write(*,*) "Wrong number of variables! Stop!"
         stop
       end if
    end do

    return
  end subroutine FFNN_wb


  !=================================================
  !Read or Write NN parameters
  !=================================================
  subroutine wt_bias(id,newwb)
    implicit none
    integer, intent(in) :: id
    real :: newwb(Nvar)
    integer :: i,j,cnt,iNeuron

    cnt=0
    do i = 1,Nlayer
      iNeuron=Neuron0(i+1)
      do j = 1,Neuron0(i)
        if(id.eq.0) then
          wt0(1:iNeuron,j,i)=newwb(cnt+1:cnt+iNeuron)
        else
          newwb(cnt+1:cnt+iNeuron)=wt0(1:iNeuron,j,i)
        end if
        cnt=cnt+iNeuron
      end do

      if(id.eq.0) then
          bias0(1:iNeuron,i)=newwb(cnt+1:cnt+iNeuron)
      else
          newwb(cnt+1:cnt+iNeuron)=bias0(1:iNeuron,i)
      end if
        cnt=cnt+iNeuron
    end do

    if(cnt.ne.Nvar) stop 'Wrong number of variables'

    return
  end subroutine wt_bias


  !=================================================
  !Save the NN parameters
  !=================================================
  subroutine nnsave(iun,io)
    implicit none
    integer, intent(in) :: io
    integer, intent(in) :: iun
    integer :: i,ios
    real,dimension(:),allocatable :: newwb

    allocate(newwb(Nvar))
    if(io.eq.0) then
      do i=1,Nvar
        read(iun,*,iostat=ios) newwb(i)
        if(ios.ne.0) stop 'Error reading NN parameter file!'
      end do
      call wt_bias(0,newwb)
    else
      call wt_bias(1,newwb)
      do i=1,Nvar
        write(iun,'(F15.8)') newwb(i)
      end do
    end if

    close(iun)
    deallocate(newwb)
    return
  end subroutine nnsave


  !==================================================
  ! Get a file handle
  !==================================================
  subroutine get_file(iun)
    integer,intent(out)::iun
    integer::k
    logical::b

    k=20
    inquire(unit=k,opened=b)

    do while(b .and. k<=100)
       k=k+1
       inquire(unit=k,opened=b)
    end do
    if(.not. b) then
       iun=k
    else
       stop "get_file: no free unit"
    end if

    return
  end subroutine get_file

 
  !==================================================
  ! Apply activation function
  !==================================================
  subroutine factivate(input0,input1,activation)
    implicit none
    real,dimension(:),intent(in)::input0
    real,dimension(:),intent(inout)::input1
    real,parameter :: fac=2.d0/acos(-1.d0)
    real,parameter :: fac1=0.5d0*acos(-1.d0)
    integer::activation
    integer::i,j

    !case 1: logsig
    !case 2: tansig
    !case 3: purelin
    !case 4: arctan
    if(activation .eq. 1) then
      input1 = exp(input0)
      input1 = input1/(1+input1)
    else if(activation .eq. 2) then
          input1 = tanh(input0)
    else if(activation .eq. 3) then
          input1 = input0
    else if(activation .eq. 4) then
          input1 = fac*atan(fac1*input0)
    else
        stop "Wrong Activation Function!"
    end if

    return
  end subroutine factivate

  !==================================================
  ! Derivative of activation function
  !==================================================
  subroutine dfactivate(input0,dydinput0,activation)
    implicit none
    real,dimension(:),intent(in)::input0
    real,dimension(:),intent(inout)::dydinput0
    real,parameter :: fac1=0.5d0*acos(-1.d0)
    integer::activation
    integer::i,j

    !case 1: logsig
    !case 2: tansig
    !case 3: purelin
    !case 4: arctan
    if(activation .eq. 1) then
      dydinput0 = exp(input0)
      dydinput0 = dydinput0/((1+dydinput0)**2)
    else if(activation .eq. 2) then
          dydinput0 = 1.d0/cosh(input0)**2
    else if(activation .eq. 3) then
          dydinput0 = 1.d0
    else if(activation .eq. 4) then
          dydinput0=1.d0/(1.d0+(fac1*input0)**2)
    else
        stop "Wrong Activation Function!"
    end if

    return
  end subroutine dfactivate

end module NN_class  



