

program inplementation
  implicit none

  real(8), parameter :: eps=1d-4    !! DSA RK45 Precision
  integer, parameter :: OMEGA=100   !! SSA Pecision (nb of part for 1uM) 
  integer, parameter :: OMEGA_max=4*OMEGA  !! SSA MAXIMUN PART NUMBER


  integer, parameter :: nbase=2        !! Number of state of one atom (eg c12 c13)
  integer, parameter :: n_cl_max=7        !! Maximum atom chain length  
  integer, parameter :: nstate=nbase**n_cl_max !! nb d'état possible au max (pour S7P))

  integer, parameter :: nv=13,nr=15    !! nb variable, nb of reactions 


  real(8) :: X(nv,nstate)      !! variable pour le deterministe
  real(8) :: Xinput(nv,nstate) !! variable d'initialisation (marquage)
  real(8) :: M(nv,n_cl_max+1)   !!   MID TABLE
  real(8) :: k(2*nr)            !! reaction flux (1:nr forward, nr+1:2nr backward)  


  type variab
     integer :: np
     integer :: s(n_cl_max)=0    
     integer :: nb=OMEGA                !! nb of particles
     integer :: pop(n_cl_max,OMEGA_max) !! isotopomer population for SSA 
     !integer :: mid(OMEGA_max)
     character(4) :: lab                !! label
  end type variab

  type reac !A+B<-> C+D
     integer :: nn
     integer :: na,nb                   !! reactant index
     integer :: nc,nd                   !! product index
     integer :: suf(2*n_cl_max,2)=0     !! permutation rule
     logical :: rev                     !! reaction is reversible
     logical :: input                   !! input reaction 
     real(8) :: kf,kb                   !! reaction rate
  end type reac

  type(variab) :: S(nv)  !! TABLE CHEMICAL SPECIES
  type(reac)   :: R(nr)  !! TABLE CHEMICAL REACTIONS

  real(8) :: initial_concentration(nv)  !! abondance relative des espèces

  integer :: iv,ir

  real(8) :: t,dt
  real(8), parameter :: Tmax=50d0  !! evolution duration (in sec)
  integer, parameter :: NT=500     !! number of output steps

  integer :: ncu=0

  character(300) :: filename




  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%% SYSTEM  & VARIABLES DEFINITIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  ! !!----------------------------
  !! §§ Glycolysis INTERNAL REACTIONS
  ! v1 G6P (abcdef)  -> F6P (abcdef)
  ! v2 F6P (abcdef) <-> FBP (abcdef)
  ! v3 FBP (acbdef) -> DHAP (cba) + GAP (def)
  ! v4 DHAP (abc) -> GAP (abc)

  ! S1 : G6P
  ! S2 : F6P
  ! S3 : FBP
  ! S4 : DHAP
  ! S5 : GAP

  S(1)%lab='G6P'
  S(2)%lab='F6P'
  S(3)%lab='FBP'
  S(4)%lab='DHAP'
  S(5)%lab='GAP'


  ! nb carbon variables
  S(1:3)%np=6 
  S(4:5)%np=3 



  ! v1 G6P (abcdef)  <-> F6P (abcdef)
  ir=1 
  R(ir)%nn=6
  R(ir)%na=1; R(ir)%nb=0;
  R(ir)%nc=2; R(ir)%nd=0
  R(ir)%suf(1:6,1)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%suf(1:6,2)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0



  ! v2 F6P (abcdef) <-> FBP (abcdef)
  ir=2
  R(ir)%nn=6
  R(ir)%na=2; R(ir)%nb=0;
  R(ir)%nc=3; R(ir)%nd=0
  R(ir)%suf(1:6,1)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%suf(1:6,2)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0

  ! v3 FBP (acbdef) <-> DHAP (cba) + GAP (def)
  ir=3
  R(ir)%nn=6
  R(ir)%na=3; R(ir)%nb=0;
  R(ir)%nc=4; R(ir)%nd=5
  R(ir)%suf(1:6,1)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%suf(1:6,2)=(/+3,+2,+1,-1,-2,-3/) 
  !  R(ir)%suf(1:6,2)=(/+1,+2,+3,-1,-2,-3/) 
  !  R(ir)%suf(1:6,2)=(/+1,+2,+3,-3,-2,-1/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0


  ! v4 DHAP (abc) <-> GAP (abc)
  ir=4
  R(ir)%nn=3
  R(ir)%na=4; R(ir)%nb=0;
  R(ir)%nc=5; R(ir)%nd=0
  R(ir)%suf(1:3,1)=(/+1,+2,+3/) 
  R(ir)%suf(1:3,2)=(/+1,+2,+3/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0




  ! !!----------------------------
  !! §§ PPP INTERNAL REACTIONS
  ! v5 G6P (abcdef) -> 6PG (abcdef)
  ! v6 6PG (abcdef) -> Ru5P (bcdef)
  ! v7 Ru5P (abcde) <-> R5P (abcde)
  ! v8 Ru5P (abcde) <-> X5P (abcde)
  ! v9 X5P (abcde) + R5P (ABCDE) <-> S7P (abABCDE) + GAP (cde)
  ! v10 S7P (abcdefg) + GAP (ABC) <-> F6P (abcABC) + E4P (defg)
  ! v11 X5P (abcde) + E4P (ABCD) <-> F6P (abABCD) + GAP (cde)


  ! S6 : 6PG
  ! S7 : Ru5P
  ! S8 : CO2
  ! S9 : X5P
  ! S10 : R5P
  ! S11 : E4P
  ! S12 : S7P

  S(6 )%lab='6PG'
  S(7 )%lab='Ru5P'
  S(8 )%lab='CO2'
  S(9 )%lab='X5P'
  S(10)%lab='R5P'
  S(11)%lab='E4P'
  S(12)%lab='S7P'



  ! nb carbon variables
  S(6)%np=6 
  S(7)%np=5 
  S(8)%np=1
  S(9:10)%np=5 
  S(11)%np=4 
  S(12)%np=7 



  ! v5  G6P (abcdef) -> 6PG (abcdef)
  ir=5
  R(ir)%nn=6
  R(ir)%na=1; R(ir)%nb=0;
  R(ir)%nc=6; R(ir)%nd=0
  R(ir)%suf(1:6,1)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%suf(1:6,2)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%rev=.FALSE.
  R(ir)%kf=1D0
  R(ir)%kb=0D0

  ! v6  6PG (abcdef) -> CO2 (a) + Ru5P (bcdef)
  ir=6
  R(ir)%nn=6
  R(ir)%na=6; R(ir)%nb=0;
  R(ir)%nc=8; R(ir)%nd=7
  R(ir)%suf(1:6,1)=(/+1,+2,+3,+4,+5,+6/) 
  R(ir)%suf(1:6,2)=(/+1,-1,-2,-3,-4,-5/) 
  R(ir)%rev=.FALSE.
  R(ir)%kf=1D0
  R(ir)%kb=0D0

  ! v8 Ru5P (abcde) <-> X5P (abcde)
  ir=8
  R(ir)%nn=5
  R(ir)%na=7; R(ir)%nb=0;
  R(ir)%nc=9; R(ir)%nd=0
  R(ir)%suf(1:5,1)=(/+1,+2,+3,+4,+5/) 
  R(ir)%suf(1:5,2)=(/+1,+2,+3,+4,+5/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0


  ! v7 Ru5P (abcde) <-> R5P (abcde)
  ir=7
  R(ir)%nn=5
  R(ir)%na=7; R(ir)%nb=0;
  R(ir)%nc=10; R(ir)%nd=0
  R(ir)%suf(1:5,1)=(/+1,+2,+3,+4,+5/) 
  R(ir)%suf(1:5,2)=(/+1,+2,+3,+4,+5/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0


  ! v9 X5P (abcde) + R5P(ABCDE) -> S7P (abABCDE) + GAP (cde)
  ir=9
  R(ir)%nn=10
  R(ir)%na=9; R(ir)%nb=10;
  R(ir)%nc=12; R(ir)%nd=5
  R(ir)%suf(1:10,1)=(/+1,+2,+3,+4,+5,-1,-2,-3,-4,-5/) 
  R(ir)%suf(1:10,2)=(/+1,+2,-1,-2,-3,+3,+4,+5,+6,+7/) 
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0


  ! v10 S7P (abcdefg) + GAP (ABC) <-> F6P (abcABC) + E4P (defg)
  ir=10
  R(ir)%nn=10
  R(ir)%na=12; R(ir)%nb=5;
  R(ir)%nc=2; R(ir)%nd=11
  R(ir)%suf(1:10,1)=(/+1,+2,+3,+4,+5,+6,+7,-1,-2,-3/) 
  R(ir)%suf(1:10,2)=(/+1,+2,+3,-1,-2,-3,-4,+4,+5,+6/)
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0

  ! v11 X5P (abcde) + E4P (ABCD) <-> F6P (abABCD) + GAP (cde)
  ir=11
  R(ir)%nn=9
  R(ir)%na=9; R(ir)%nb=11;
  R(ir)%nc=2; R(ir)%nd=5
  R(ir)%suf(1:9,1)=(/+1,+2,+3,+4,+5,-1,-2,-3,-4/)
  R(ir)%suf(1:9,2)=(/+1,+2,-1,-2,-3,+3,+4,+5,+6/)
  R(ir)%rev=.TRUE.
  R(ir)%kf=1D0
  R(ir)%kb=1D0


  ! !!----------------------------
  !! §§ OUTPUT REACTIONS
  ! v12 GAP (abc) -> exit1 (abc)
  ! v13 R5P (abcde) -> exit2 (abcde)
  ! v14 CO2 (a) -> exit3 (abcde)


  ! v12 GAP (abc) -> exit1 (abc)
  ir=12
  R(ir)%nn=3
  R(ir)%na=5; R(ir)%nb=0
  R(ir)%nc=0; R(ir)%nd=0
  R(ir)%suf(1:3,1)=(/+1,+2,+3/)
  R(ir)%suf(1:3,2)=(/+1,+2,+3/)
  R(ir)%rev=.FALSE.
  R(ir)%kf=1D0
  R(ir)%kb=0D0

  ! v13 R5P (abcde) -> exit2 (abcde)
  ir=13
  R(ir)%nn=5
  R(ir)%na=10; R(ir)%nb=0
  R(ir)%nc=0; R(ir)%nd=0
  R(ir)%suf(1:5,1)=(/+1,+2,+3,+4,+5/)
  R(ir)%suf(1:5,2)=(/+1,+2,+3,+4,+5/)
  R(ir)%rev=.FALSE.
  R(ir)%kf=1D0
  R(ir)%kb=0D0

  ! v14 CO2 (a) -> exit3 (abcde)
  ir=14
  R(ir)%nn=1
  R(ir)%na=8; R(ir)%nb=0
  R(ir)%nc=0; R(ir)%nd=0
  R(ir)%suf(1:1,1)=(/+1/)
  R(ir)%suf(1:1,2)=(/+1/)
  R(ir)%rev=.FALSE.
  R(ir)%kf=1D0
  R(ir)%kb=0D0

  ! !!----------------------------
  !! §§ INPUT REACTIONS
  ! v15 GLU (abcde) -> G6P (abcde)

  ! S13 : GLU
  S(13)%np=6 
  S(13)%lab='GLU'

  ir=15
  R(ir)%nn=6
  R(ir)%na=13; R(ir)%nb=0
  R(ir)%nc=1; R(ir)%nd=0
  R(ir)%suf(1:6,1)=(/+1,+2,+3,+4,+5,+6/)
  R(ir)%suf(1:6,2)=(/+1,+2,+3,+4,+5,+6/)
  R(ir)%rev=.FALSE.
  R(ir)%input=.TRUE.
  R(ir)%kf=2D0

  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%% END SYSTEM  & VARIABLES DEFINITIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




  !! def concentration relative des especes
  initial_concentration=1d0

  !! FIXED POINT if kf=kb=1
  initial_concentration(01)=0.983
  initial_concentration(02)=0.966
  initial_concentration(03)=0.727
  initial_concentration(04)=0.828
  initial_concentration(05)=0.589
  initial_concentration(06)=0.983
  initial_concentration(07)=1.411
  initial_concentration(08)=0.983
  initial_concentration(09)=1.189
  initial_concentration(10)=0.650
  initial_concentration(11)=0.571
  initial_concentration(12)=1.125
  initial_concentration(13)=1.000


  k(1:nr)=R(:)%kf
  k(nr+1:2*nr)=R(:)%kb

  !! ---------------------------------------
  !! LABELLING DEFINITION
  !! ---------------------------------------

  !! Xinput(iv,1)=(0,0)
  !! Xinput(iv,2)=(1,0)
  !! Xinput(iv,3)=(0,1)
  !! Xinput(iv,4)=(1,1)

  Xinput(:,:)=0
  Xinput(:,1)=1

  !! INPUT LABELLING ON GLU :: variab # 13
  Xinput(13,1)=0
  Xinput(13,2)=0.5
  Xinput(13,3)=0.5
  Xinput(13,4)=0

  !! ---------------------------------------
  !! END LABELLING DEFINITION
  !! ---------------------------------------




  !! ---------------------------------------
  !! OUTPUT FILES
  !! ---------------------------------------
  
  !! CONCENTRATION TIME EVOLUTION,  ONE LINE PER DATE
  !! time var_1 var_2 var_3 ... 
  open(20,file='DSA_CONC.dat')
  open(120,file='SSA_CONC.dat')

  !! ONE FILE PER VARIABLE, ONE LINE PER DATE
  !! time concentration m+0 m+1 m+2 ....
  do iv=1,nv
     write(filename,'(A,i2.2,A,A,A)') 'DSA_V',iv,'_',trim(S(iv)%lab),'.dat'
     open(20+iv,file=filename)
     write(filename,'(A,i2.2,A,A,A)') 'SSA_V',iv,'_',trim(S(iv)%lab),'.dat'
     open(120+iv,file=filename)
  end do
  !! ---------------------------------------
  !! END OUTPUT FILES
  !! ---------------------------------------


  
  dt=TMAX/NT
  call RUN_DSA   !! integration deterniste sur NT*dt

  dt=TMAX/NT
  call RUN_SSA  !! integration stochastique sur NT*dt





contains


  subroutine verif_flux
    integer :: iv,ir

    real(8) :: STOCH(nv,2*nr),var(nv)


    ! test des reactions
    !write(6,'(20i5)') 0,(iv,iv=1,nv)
    do ir=1,2*nr
       do iv=1,nv
          S(iv)%pop(:,:)=0
          S(iv)%nb=OMEGA
       end do
       call SSA_do_reac(ir)
       STOCH(:,ir)=S(:)%nb-OMEGA
    end do

    var=matmul(STOCH,k)
    do iv=1,nv
       if (abs(var(iv)) .GT. 1D-6) then
          print*,var(iv)
          print*,'IMBALANCE FLUX, PLEASE CORRECT'
          STOP
       end if
    end do
    print*,'FLUX BALANCE OK'

    !    read(5,*)

  end subroutine verif_flux




  !! ----------------------------------------------
  !! ----------------------------------------------
  !! RUN THE SSA
  !! ----------------------------------------------
  !! ----------------------------------------------
  subroutine RUN_SSA
    integer :: iv,ir,i,j,it,ip

    real(8) :: ksto(2*nr),ksto_old(2*nr),Rlim(2*nr)
    integer :: nb_reac(2*nr)
    real(8) :: next_time(2*nr)

    integer :: iit,iito

    integer :: ns,n1,n2



    !! initial condition
    do iv=1,nv
       S(iv)%pop(:,:)=0
       S(iv)%nb=int(OMEGA*initial_concentration(iv)) !ns(iv)

       ns=nbase**S(iv)%np
       n1=0
       do j=1,ns

          do i=0,n_cl_max-1
             S(iv)%s(i+1)=merge(1,0,btest(j-1,i))
          end do

          n2=int(Xinput(iv,j)/sum(Xinput(iv,:))*S(iv)%nb)

          if (n2 .GT.0) then
             do ip=n1+1,n1+n2
                S(iv)%pop(:,ip)=s(iv)%s
             end do
             n1=n1+n2
          end if

       end do
    end do
    !! end initial condition


    call cal_ksto(ksto)

    Rlim=0D0
    !! for random reaction time
    !call random_number(next_time)
    !! for fixed reaction time
    next_time=exp(-1d0)

    next_time=-(log(next_time))/max(ksto,1D-16)

    t=0
    iit=int(t/dt)
    iito=iit

    nb_reac=0
    do it=1,NT 
 
       do !! internal loop, evolution from t -> t+dt

          ncu=ncu+1


          ir=minloc(next_time,dim=1)
          nb_reac(ir)=nb_reac(ir)+1
          call SSA_do_reac(ir)

          t=next_time(ir)
          ksto_old=ksto
          call cal_ksto(ksto)

          !! for random reaction time
          !call random_number(next_time(ir))
          !! for fixed reaction time
          next_time(ir)=exp(-1d0)

          next_time(ir)=t-(log(next_time(ir)))/max(ksto(ir),1D-16)

          next_time=t+max(ksto_old,1D-16)/max(ksto,1D-16)*(next_time-t)



          iito=iit
          iit=int(t/dt)
          if (iit .NE. iito) exit

       end do !! end internal loop


       !! RECORD MID
       call getmid(M,'SSA')
       write(120,*) t,(sum(X(iv,:)),iv=1,nv)
       do iv=1,nv
          write(120+iv,*) t,sum(M(iv,:)),M(iv,1:S(iv)%np)
       end do
       !! END RECORD MID



       if (t .GT. Tmax) exit

    end do


  end subroutine RUN_SSA


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! COMPUTE THE RATE V_k(S) for SSA
  !! ----------------------------------------------
  !! ----------------------------------------------

  subroutine cal_ksto(ksto)
    real(8) :: ksto(2*nr)
    integer :: ir 
    integer :: va,vb,vc,vd


    do ir=1,nr
       va=R(ir)%na !! reactant #1 index 
       vb=R(ir)%nb !! reactant #2 index 
       vc=R(ir)%nc !! product #1 index 
       vd=R(ir)%nd !! product #2 index 


       ksto(ir)=k(ir)*OMEGA
       if (va .NE. 0) ksto(ir)=ksto(ir)* real(S(va)%nb)/real(OMEGA)
       if (vb .NE. 0) ksto(ir)=ksto(ir)* real(S(vb)%nb)/real(OMEGA)

       ksto(nr+ir)=k(nr+ir)*OMEGA
       if (vc .NE. 0) ksto(nr+ir)=ksto(nr+ir)* real(S(vc)%nb)/real(OMEGA)
       if (vd .NE. 0) ksto(nr+ir)=ksto(nr+ir)* real(S(vd)%nb)/real(OMEGA)

    end do


  end subroutine cal_ksto




  !! ----------------------------------------------
  !! ----------------------------------------------
  !! PERFORM A REACTION IN THE SSA
  !! ----------------------------------------------
  !! ----------------------------------------------

  subroutine SSA_do_reac(iir)
    integer, intent(in) :: iir
    integer :: ir, sens

    integer :: va,vb,vc,vd
    integer :: na,nb,nc,nd

    if (iir .LE. nr) then
       ir=iir
       sens=+1

       va=R(ir)%na !! reactant #1 index 
       vb=R(ir)%nb !! reactant #2 index 
       vc=R(ir)%nc !! product #1 index 
       vd=R(ir)%nd !! product #2 index 

    else
       ir=iir-nr
       sens=-1

       va=R(ir)%nc !! reactant #1 index 
       vb=R(ir)%nd !! reactant #2 index 
       vc=R(ir)%na !! product #1 index 
       vd=R(ir)%nb !! product #2 index 
    end if


    if (va.NE.0) then
       na=int(S(va)%nb*rand()+1)
       S(va)%s=S(va)%pop(:,na)
    end if

    if (vb.NE.0) then
       nb=int(S(vb)%nb*rand()+1)
       S(vb)%s=S(vb)%pop(:,nb)
    end if

    call doreac(ir,sens)

    if (.not. R(ir)%input) then 
       if (va.NE.0) then
          S(va)%pop(:,na)=S(va)%pop(:,S(va)%nb)
          S(va)%nb=S(va)%nb-1

          if (S(va)%nb .EQ. 0) then
             print*, 'Pb: stock à 0 sur variable ',va
             stop
          end if
       end if
       if (vb.NE.0) then
          S(vb)%pop(:,nb)=S(vb)%pop(:,S(vb)%nb)
          S(vb)%nb=S(vb)%nb-1

          if (S(vb)%nb .EQ. 0) then
             print*, 'Pb: stock à 0 sur variable ',vb
             stop
          end if

       end if
    end if

    if (vc.NE.0) then
       nc=S(vc)%nb+1
       S(vc)%pop(:,nc)=S(vc)%s
       S(vc)%nb=nc
    end if

    if (vd.NE.0) then
       nd=S(vd)%nb+1
       S(vd)%pop(:,nd)=S(vd)%s
       S(vd)%nb=nd
    end if

  end subroutine SSA_do_reac


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! COMPUTE MID FROM ISOTOPOMERS
  !! ----------------------------------------------
  !! ----------------------------------------------
  subroutine getMID(M,typ)
    real(8) ,intent(out) :: M(nv,n_cl_max+1)
    character(3), intent(in) :: typ
    integer :: iv , ns,nd,is,id,nmid,j
    M=0


    do iv=1,nv

       if (typ=='DSA') then 
          nd=S(iv)%np  !! nb de digit
          ns=nbase**n_cl_max !! nb état possible

          do is=1,ns !parcours les états
             nmid=1
             do id=0,nd-1
                j=merge(1,0,btest(is-1,id))
                nmid=nmid+j !! test le bit
             end do

             M(iv,nmid)=M(iv,nmid)+X(iv,is)

          end do
       else  if (typ=='SSA') then 
          ns=S(iv)%nb !! nb état possible
          do is=1,ns !parcour les états
             nmid=sum(S(iv)%pop(:,is))+1
             M(iv,nmid)=M(iv,nmid)+1d0/real(OMEGA)
          end do

       end if

    end do


  end subroutine GETMID


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! RUN THE DSA 
  !! ----------------------------------------------
  !! ----------------------------------------------

  subroutine RUN_DSA
    real(8) :: t

    !! INITIAL CONDITION
    t=0
    do iv=1,nv
       X(iv,:)=Xinput(iv,:)*initial_concentration(iv)
    end do

    !! TIME EVOLUTION
    do 

       !! RECORD MID
       call getmid(M,'DSA')
       write(20,*) t,(sum(X(iv,:)),iv=1,nv)
       do iv=1,nv
          write(20+iv,*) t,sum(M(iv,:)),M(iv,1:S(iv)%np)
       end do
       !! END RECORD MID

       call DSA_evol(t,dt)

       if (t .GT. Tmax) exit

    end do

  end subroutine RUN_DSA


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! RUNGE-KUTTA 4-5  FRO DSA
  !! ----------------------------------------------
  !! ----------------------------------------------
  subroutine DSA_evol(t,dt) 
    real(8), intent(inout) :: t,dt
    real(8), dimension(nv,nstate) :: k1,k2,k3,k4,k5,k6,errr
    real(8) :: te


    do
       ! RK4-5
       call DSA_RHS(t,X,k1)
       call DSA_RHS(t+0.5*dt , X+dt*0.5*k1 , k2)
       call DSA_RHS(t+0.5*dt , X+dt*(0.25*k1+0.25*k2) , k3)
       call DSA_RHS(t+1.0*dt , X+dt*(-1.0*k2+2.0*k3) , k4)
       call DSA_RHS(t+2.*dt/3. , X+dt*(7*k1+10*k2+1*k4)/27. , k5)
       call DSA_RHS(t+1.*dt/5. , X+dt*(28.*k1-125.*k2+546.*k3+54.*k4-378.*k5)/665. , k6)


       errr=dt*(1.*k1/360.-128.*k3/4275-2187.*k4/75240.+1.*k5/50.+2.*k6/55.)
       te=maxval(abs(errr))


       if (te .LT. eps) then 
          X=X+(1.*k1/24.+5.*k4/48.+27.*k5/56.+125.*k6/336.)*dt

          t=t+dt
          exit
       else

          dt=0.9*dt*(eps/te)**(1./5.)

       end if

    end do
    if (te .LT. eps/2.)  dt=min(2d0*dt,10d0)


!!$    ! RK4
!!$    call DSA_RHS(t,X,k1)
!!$    call DSA_RHS(t+0.5*dt , X+dt*0.5*k1 , k2)
!!$    call DSA_RHS(t+0.5*dt , X+dt*0.5*k2 , k3)
!!$    call DSA_RHS(t+1.0*dt , X+dt*1.0*k3 , k4)
!!$    X=X+(k1+2*k2+2*k3+k4)*dt/6.
!!$    t=t+dt


    !  EULER
!!$    call DSA_RHS(t,X,k1)
!!$    X=X+k1*dt
!!$    t=t+dt

  end subroutine DSA_evol


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! RHS FOR DSA
 !! ----------------------------------------------
  !! ----------------------------------------------

  subroutine DSA_RHS(t,X,DXDT)
    real(8), intent(in) :: t,X(nv,nstate)
    real(8), intent(out) :: DXDT(nv,nstate)

    integer :: ir

    ncu=ncu+1

    DXDT=0
    do ir=1,nr
       if (R(ir)%rev) then 
          if (R(ir)%kf .GT. 1d-6)  call DSA_DO_REAC(ir,+1,X,DXDT)
          if (R(ir)%kb .GT. 1d-6) call DSA_DO_REAC(ir,-1,X,DXDT)
       else
          if (R(ir)%kf .GT. 1d-6) call DSA_DO_REAC(ir,+1,X,DXDT)
       end if
    end do



  end subroutine DSA_RHS


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! RHS FOR ONE REACTION OF THE DSA 
  !! ----------------------------------------------
  !! ----------------------------------------------
  subroutine DSA_DO_REAC(ir,sens,X,DXDT)
    integer, intent(in) :: ir,sens
    real(8), intent(in) :: X(nv,nstate)
    real(8), intent(inout) :: DXDT(nv,nstate)

    integer :: i , v1,v2,v3,v4
    integer :: n1,n2,i1,i2

    integer :: va,vb,vc,vd

    integer :: ns1,ns2,ns3,ns4

    real(8) :: reac_rate,taux


    !print*,ir,real(R(ir)%k)


    va=R(ir)%na !! reactant #1 index 
    vb=R(ir)%nb !! reactant #2 index 
    vc=R(ir)%nc !! product #1 index 
    vd=R(ir)%nd !! product #2 index 



    if (sens==1) then
       v1=va
       v2=vb
       v3=vc
       v4=vd
       taux=R(ir)%kf
    else if (sens==-1) then
       v1=vc
       v2=vd
       v3=va
       v4=vb
       taux=R(ir)%kb
    end if

    n1=nbase**S(v1)%np
    if (v2==0) then
       n2=1
    else
       n2=nbase**S(v2)%np
    end if



    do i1=1,n1
       do i2=1,n2

          do i=0,n_cl_max-1
             S(v1)%s(i+1)=merge(1,0,btest(i1-1,i))
             if (v2 .NE. 0) S(v2)%s(i+1)=merge(1,0,btest(i2-1,i))
          end do


          call doreac(ir,sens)


          ns1=0
          ns2=0
          ns3=0
          ns4=0
          do i=1,n_cl_max
             if (S(v1)%s(i)==1) ns1=ibset(ns1,i-1)
             if (v2 .NE. 0) then 
                if (S(v2)%s(i)==1) ns2=ibset(ns2,i-1)
             end if
             if (v3 .NE. 0) then 
                if (S(v3)%s(i)==1) ns3=ibset(ns3,i-1)
             end if
             if (v4 .NE. 0) then 
                if (S(v4)%s(i)==1) ns4=ibset(ns4,i-1)
             end if
          end do
          ns1=ns1+1 !! numero de la configurtion
          ns2=ns2+1
          ns3=ns3+1
          ns4=ns4+1


          if (v2==0) then
             reac_rate=taux*X(v1,ns1)
          else
             reac_rate=taux*X(v1,ns1)*X(v2,ns2)
          end if



          if (.NOT. R(ir)%input)  DXDT(v1,ns1)=DXDT(v1,ns1)-reac_rate
          if (v2 .NE. 0) DXDT(v2,ns2)=DXDT(v2,ns2)-reac_rate
          if (v3 .NE. 0) DXDT(v3,ns3)=DXDT(v3,ns3)+reac_rate
          if (v4 .NE. 0) DXDT(v4,ns4)=DXDT(v4,ns4)+reac_rate

       end do
    end do



  end subroutine DSA_DO_REAC


  !! ----------------------------------------------
  !! ----------------------------------------------
  !! FERFORM ONE REACTION WITH SPECIFIC ISOTOPOMERS 
  !! ----------------------------------------------
  !! ----------------------------------------------
  subroutine doreac(ir,sens)
    integer :: ir,sens
    integer :: i , ni,v1,v2,l1,l2


    !! do the reaction
    ni=R(ir)%nn
    do i=1,ni
       if (R(ir)%suf(i,1) .GT. 0) then
          v1=R(ir)%na
       else 
          v1=R(ir)%nb
       end if
       l1=abs(R(ir)%suf(i,1))

       if (R(ir)%suf(i,2) .GT. 0) then
          v2=R(ir)%nc
       else 
          v2=R(ir)%nd
       end if
       l2=abs(R(ir)%suf(i,2))

       !! sens direct
       if (sens ==1 ) then 
          if (v2.NE.0) S(v2)%s(l2)=S(v1)%s(l1)
       else if (sens ==-1 ) then 
          !! sens reverse
          if (v2.NE.0) S(v1)%s(l1)=S(v2)%s(l2)
       end if
    end do
    !!


  end subroutine doreac


end program inplementation


