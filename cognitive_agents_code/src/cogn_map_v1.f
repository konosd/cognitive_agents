c Cognitive map model
c (c) Maksym Romenskyy
c  22 Dec 2018

c Variables                 Description

c ### INTEGER ###
c i, j                      general int for loops
c np                        number of particles
c st                        current step
c nsteps                    number of steps
c tr                        current vtraj
c nvtraj                    number of vtraj
c st_virt                   current step on vtraj
c tau                       number of steps in vtraj
c L2                        integer half size of the box
c start                     first timestep for analysis
c incr                      increment at which data is saved
c incr_coor                 increment at which coordinates are saved
c NBIN                      number of bins in analysis
c count                     used in VACF
c BL                        used in VACF
c seqn                      counter for simulation run
c tau_i                     increased tau in swarming
c max_tau                   artificial stop for tau_i
c last_step                 tau for a particle in single virtual trajectory
c extra_tau                 tau bonus in swarming

c ### REAL ###
c sigma                     size of the agents
c L                         size of the box
c noi                       used in Langevin sqrt(2.0*gamma*temp/dt)
c gamma                     drag
c temp                      temperature
c dt                        timestep
c Ar                        strength of repulsion
c pi                        pi constant as pi=4.0d0*atan(1.0d0)
c Fx, Fy                    total force
c Flx, Fly                  random force (from Langevin equation)
c Rtempx, Rtempy            coordinates after integration before box wrap adjustment
c Gaus_x, Gaus_y            noise for Langevin
c Rgyr_norm                 normalisation factore for Rgyr
c rgx, rgy                  radius of gyration (sum of mean square displacement from center)
c r_mean_x, r_mean_y        center of Rgyr
c R2virt_x, R2virt_y        x and y components of modulo current virt coordinates
c Rbx,  Rby                 x and y components of past coordinates
c Vvirt_x, Vvirt_y          virtual velocity components
c Frepx_virt, Frepy_virt    total repulsion after virtual step
c Frepx, Frepy              total repulsion after virtual step (in SUBROUTINE)
c skor                      temp skore
c Rijx, Rijy                distance between particles (per axis)
c Rij                       Euclidian distance between particles
c Rux, Ruy                  repulsion force modification after collision
c b                         force shape
c polar, apolar             polar and apolar order parameters
c tor                       angular order parameter
c dur                       time to reach steady state
c nl(7)                     pnf mean no of particles per cell
c dnl(7)                    pnf equation (SD of particle no)
c dnl1(7)                   average pnf
c Vmx, Vmy                  x and y components of the initial velocity in VACF
c Vm                        Euclidian of Vmx and Vmy
c Va                        velocity autocorrelation
c Va1(200)
c V21                       Euclidian of x and y from V @ st in VACF
c pull                      atraction distance modification for swarming

c ### REAL ARRAYS ###
c R                         (np, 2) current coordinates
c Rb                        (np, 2) previous coordinates
c R2                        (np, 2) modulo current coordinates
c V                         (np, 2) velocities
c Rvirt                     (max_tau+1, 2) virtual coordinates
c gnumb_x, gnumb_y          (max_tau) random numbers for virtual trajectories
c Rvirt1, Rvirt2            (max_tau) deprecated, now use r_mean_x, r_mean_y
c Gr                        (L2) current RDF
c Gr1                       (nsteps/incr, L2) stores RDFs
c CCF                       (L2) curretn CCF
c CCF1                      (nsteps/incr, L2) stores CCFs
c Gausv_x, Gausv_y          (nvtraj) x and y components of noise
c Rgyr                      (nvtraj) ratio of gyration

c ### CHARACTER ###
c repultype                 repulsion type [repul_hcore, repul_soft, repul_ahc]
c tau_char, np_char, nsteps_char, n_char, extra_tau_char, pull_char,
c file_name_char
c dir_name,file_name,file,path,file_path
c vtraj1, tau1, sigma1, L1, a1, seqn1, nsteps1, np1
c InfoFile, CoorFile, CoorlFile, PolFile, ApolFile, AngFile, RDFFile, VACFFile,
c CCFFile, PNFFile, TempFile, TauFile

c ### LOGICAL ###
c d_exists, f_exists        directory exists, file exists
c result                    making a new directory
c status                    deprecated change directory
c save_temp, save_coor, save_ord, save_vacf

       PROGRAM cognitive

         use IFPORT !for making new directory
         implicit none

        INTEGER i,j,np,st,nsteps,tr,nvtraj,st_virt,tau,L2
        INTEGER start,incr,incr_coor,count,NBIN,BL,seqn
        INTEGER tau_i,max_tau,last_step,extra_tau

        REAL sigma,L,noi,gamma,temp,dt,Ar,pi,Fx,Flx,Fy,Fly
        REAL Rtempx,Rtempy,Gaus_x,Gaus_y,Rgyr_norm,rgx,rgy,dnl(7)
        REAL r_mean_x,r_mean_y,R2virt_x,R2virt_y,Rbx,Rby,nl(7)
        REAL Vvirt_x,Vvirt_y,Frepx_virt,Frepy_virt,Frepx,Frepy,skor
        REAL Rijx,Rijy,Rij,Rux,Ruy,b,polar,apolar,tor,dur,dnl1(7)
        REAL Vmx,Vmy,Va,Va1(200),Vm,V21,pull

        REAL, DIMENSION(:,:), ALLOCATABLE :: R,Rb,R2,V,Rvirt
        REAL, DIMENSION(:), ALLOCATABLE :: gnumb_x,gnumb_y
        REAL, DIMENSION(:), ALLOCATABLE :: Rvirt1,Rvirt2,Gr,CCF
        REAL, DIMENSION(:,:), ALLOCATABLE :: Gr1,CCF1
        REAL, DIMENSION(:), ALLOCATABLE :: Gausv_x,Gausv_y,Rgyr

        CHARACTER(len=90)::nvtraj1,tau1,sigma1,L1,a1,seqn1,nsteps1
        CHARACTER(len=90)::np1,repultype,InfoFile,CoorFile,CoorlFile
        CHARACTER(len=90)::PolFile,ApolFile,AngFile,RDFFile,VACFFile
        CHARACTER(len=90)::CCFFile,PNFFile,TempFile,TauFile
        CHARACTER(len=90)::tau_char,np_char,nsteps_char,n_char
        CHARACTER(len=90)::dir_name,file_name,file,path,file_path
        CHARACTER(len=90)::extra_tau_char,pull_char,file_name_char

        LOGICAL save_temp,save_coor,save_ord,save_vacf
        LOGICAL d_exists,f_exists,result,status

c getting the command line arguments

        CALL GET_COMMAND_ARGUMENT(1,np_char)   !load the command line arguments
        CALL GET_COMMAND_ARGUMENT(2,tau_char)
        CALL GET_COMMAND_ARGUMENT(3,nsteps_char)
        CALL GET_COMMAND_ARGUMENT(4,file_name)
        !CALL GET_COMMAND_ARGUMENT(4,extra_tau_char)
        !CALL GET_COMMAND_ARGUMENT(5,pull_char)

        READ(np_char,*)np      !convert them from CHAR (default) to INTEGER
        READ(tau_char,*)tau
        READ(nsteps_char,*)nsteps
        !READ(extra_tau_char,*)extra_tau
        !READ(pull_char,*)pull

c Simulation parameters
        ! Real trajectories parameters
  !      np = 200            ! number of real particles
  !      nsteps = 5000       ! number of real steps

        ! Virtual trajectories parameters
        nvtraj = 20         ! number of virtual trajectories
  !      tau = 20            ! number of virtual steps
        vrepultype="h"      ! h - hard, ah - adjust hard, s - soft, sw - swarm

        ! General parameters
        sigma = 1.0         ! particle size
        Ar = 20.0           ! force strength
        b = 200.0           ! force shape
        L = 80.0*sigma      ! box size
        dt = 0.1            ! time step
        repultype="s"      ! h - hard, ah - adjust hard, s - soft, sw - swarm

        ! Swarming parameters
        max_tau = tau       ! tau limitation for swarming
        extra_tau = 10      ! in swarming
        pull = 2
        pull = sigma*pull     ! zone of atraction for swarming


        ! Thermostat
        gamma = 1.0
        temp = 1.0

c Analysis parameters

        start = 1          ! starting timestep for analysis - should be ss
        incr = 10          ! analyse every incr timestep
        incr_coor = 50     ! save coor and vel every incr timestep
        NBIN = 100         ! number of bins

c Types of Analysis

        save_temp=.true.        ! compute temperature
        save_coor=.true.        ! write coordinates every incr_coor step
        save_ord=.true.        ! compute order parameters
        save_vacf=.true.       ! compute VACF (velocity autocorr function)

c Setup
        ! Create unique filename with simulation parameters encoded
        write(np1,*)np
        write(nvtraj1,*)nvtraj
        write(tau1,*)tau
        write(sigma1,*)int(sigma)
        write(L1,*)int(L)
        write(a1,*)int(Ar)
        write(seqn1,*)seqn
        write(nsteps1,*)nsteps
        dir_name=trim(adjustl("n"))//trim(adjustl(np1))//
     *  trim(adjustl("nvt"))//trim(adjustl(nvtraj1))//
     *  trim(adjustl("tau"))//trim(adjustl(tau1))//
     *  trim(adjustl("s"))//trim(adjustl(sigma1))//
     *  trim(adjustl("L"))//trim(adjustl(L1))//
     *  trim(adjustl("a"))//trim(adjustl(a1))//
     *  trim(adjustl(repultype))//
     *  trim(adjustl("nsteps"))//trim(adjustl(nsteps1))

        INQUIRE (DIRECTORY=dir_name, EXIST=d_exists)

        if (.not.d_exists) then
            result = MakeDirQQ(dir_name)
        end if

        !do i=1,10
            !write(n_char,'(I5.5)')i ! formatting
            !file = trim(adjustl(n_char))//"_coor.dat"
            !status = CHANGEDIRQQ(dir_name)
            !INQUIRE (FILE=file, EXIST=f_exists)
            !if (.not.f_exists) then
            !file_name = trim(adjustl(n_char))
            !status = CHANGEDIRQQ('..')
            !exit
            !end if
        !end do

        path = trim(dir_name)//'/'//trim(file_name)

        ! Write info not encoded in filename into a *.txt file
        write(InfoFile,'(a,a)')trim(dir_name),"_info.txt"
        OPEN(51,FILE=InfoFile,form='formatted')
        write(51,*)"nsteps",nsteps
        write(51,*)"dt",dt
        write(51,*)"gamma",gamma
        write(51,*)"temp",temp
        write(51,*)"start",start
        write(51,*)"incr",incr
        write(51,*)"incr_coor",incr_coor
        CLOSE(51)

        ! Output setup
        IF(save_coor)write(CoorFile,'(a,a)')trim(path),"_coor.dat"
        write(CoorlFile,'(a,a)')trim(path),"_coorl.dat"
        IF(save_ord)write(PolFile,'(a,a)')trim(path),"_pol.dat"
        IF(save_ord)write(ApolFile,'(a,a)')trim(path),"_apol.dat"
        IF(save_ord)write(AngFile,'(a,a)')trim(path),"_ang.dat"
        write(RDFFile,'(a,a)')trim(path),"_rdf.dat"
        write(CCFFile,'(a,a)')trim(path),"_ccf.dat"
        write(PNFFile,'(a,a)')trim(path),"_pnf.dat"
        IF(save_temp)write(TempFile,'(a,a)')trim(path),"_temp.dat"
        IF(save_vacf)write(VACFFile,'(a,a)')trim(path),"_VACF.dat"
        write(TauFile,'(a,a)')trim(path),"_tau_i.dat"

        OPEN(49,FILE=CoorlFile,form='formatted')
        IF(save_coor)OPEN(50,FILE=CoorFile,form='formatted')
        IF(save_ord)OPEN(52,FILE=PolFile,form='formatted')
        IF(save_ord)OPEN(53,FILE=ApolFile,form='formatted')
        IF(save_ord)OPEN(54,FILE=AngFile,form='formatted')
        OPEN(55,FILE=RDFFile,form='formatted')
        OPEN(56,FILE=CCFFile,form='formatted')
        OPEN(57,FILE=PNFFile,form='formatted')
        IF(save_temp)OPEN(60,FILE=TempFile,form='formatted')
        IF(save_vacf)OPEN(61,FILE=VACFFile,form='formatted')
        OPEN(62,FILE=TauFile,form='formatted')

        ! Memory and randomization setup
        call init_random_seed()
        allocate(R(np,2))
        allocate(Rb(np,2))
        allocate(R2(np,2))
        allocate(V(np,2))
        allocate(Gausv_x(nvtraj))
        allocate(Gausv_y(nvtraj))
        allocate(Rgyr(nvtraj))

        allocate(Rvirt(max_tau+1,2))
        allocate(gnumb_x(max_tau))
        allocate(gnumb_y(max_tau))
        allocate(Rvirt1(max_tau))
        allocate(Rvirt2(max_tau))

        L2=int(L)*0.5
        allocate(Gr(L2))
        allocate(CCF(L2))
        allocate(Gr1(nsteps/incr,L2))
        allocate(CCF1(nsteps/incr,L2))

        pi=4.0d0*atan(1.0d0)
        count=0
        Gr1(1:(nsteps/incr),1:L2)=0.0
        CCF1(1:(nsteps/incr),1:L2)=0.0
        dnl1(1:7)=0.0
        skor=0.0
        BL=start
        Va1=0.0

c Initial conditions
        ! initialise particles
        call rand_nonoverlap_circle(R,np,sigma,sigma+sigma/10,L)
        noi=sqrt(2.0*gamma*temp/dt)
        Rb(:,1)=R(:,1)
        Rb(:,2)=R(:,2)
        R2(:,1)=R(:,1)
        R2(:,2)=R(:,2)
        DO i=1,np
            V(i,1)=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
            V(i,2)=sqrt(-2.0*log(rand(0)))*sin(2*pi*rand(0))
        END DO

c Main loop
        DO st=1,nsteps !for all timesteps
            DO i=1,np ! for all particles
                ! Virtual trajectory for particle i at timestep st
                DO tr=1,nvtraj ! for all virtual trajectories
                    ! Setup for the first virtual step
                    Rvirt(1,1)=R(i,1)
                    Rvirt(1,2)=R(i,2)
                    R2virt_x=R2(i,1)
                    R2virt_y=R2(i,2)
                    Rbx=Rb(i,1)
                    Rby=Rb(i,2)
                    Vvirt_x=V(i,1)
                    Vvirt_y=V(i,2)

        ! Generate random numbers for the whole virtual trajectory
              DO st_virt=1,max_tau
              gnumb_x(st_virt)=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
              gnumb_y(st_virt)=sqrt(-2.0*log(rand(0)))*sin(2*pi*rand(0))
              END DO
            ! Save the first-step-random-numbers
                Gausv_x(tr)=gnumb_x(1)
                Gausv_y(tr)=gnumb_y(1)

                tau_i=tau
                st_virt=1
                ! Virtual steps
                DO WHILE (st_virt.le.tau_i+1)
                    if (st_virt.eq.max_tau) exit
                    if (vrepultype.eq.'hs') then
                      call repul_hcore(np,L,sigma,Ar,R2virt_x,R2virt_y,R2,
                   *  Frepx_virt,Frepy_virt)
                    else if (vrepultype.eq.'sw') then
                       call swarm(np,L,R2virt_x,R2virt_y,R2,
                    *  st_virt,tau_i,extra_tau,pull)
                    else if (vrepultype.eq.'s') then
                      call repul_soft(np,L,sigma,Ar,R2virt_x,R2virt_y,R2,
                           *  Frepx_virt,Frepy_virt)
                    else if (vrepultype.eq.'ah')
                      call repul_ahc(np,L,sigma,Ar,b,R2virt_x,R2virt_y,R2,
                        *  Frepx_virt,Frepy_virt)
                    end if

                    ! Virtual Langevin equation
                    Flx=-gamma*Vvirt_x+noi*gnumb_x(st_virt)
                    Fly=-gamma*Vvirt_y+noi*gnumb_y(st_virt)

                    ! Virtual total force
                    Fx=Flx+Frepx_virt
                    Fy=Fly+Frepy_virt

                    ! Virtual integration
                    Rtempx=2*Rvirt(st_virt,1)-Rbx+Fx*dt**2
                    Rtempy=2*Rvirt(st_virt,2)-Rby+Fy*dt**2
                    Vvirt_x=(Rtempx-Rvirt(st_virt,1))/dt
                    Vvirt_y=(Rtempy-Rvirt(st_virt,2))/dt
                    Rbx=Rvirt(st_virt,1)
                    Rby=Rvirt(st_virt,2)
                    Rvirt(st_virt+1,1)=Rtempx
                    Rvirt(st_virt+1,2)=Rtempy
                    R2virt_x=MODULO(Rvirt(st_virt,1),L)
                    R2virt_y=MODULO(Rvirt(st_virt,2),L)

                    st_virt=st_virt+1
                END DO   ! End virtual steps for trajectory i

                last_step=st_virt-1
                write(62,*)last_step
                    ! Radius of gyration of trajectory i
        !            Rvirt1(:)=Rvirt(:,1)
        !            Rvirt2(:)=Rvirt(:,2)
        !            r_mean_x=sum(Rvirt1)/last_step
        !            r_mean_y=sum(Rvirt2)/last_step

                r_mean_x=0
                r_mean_y=0
                DO j=1,last_step+1
                    r_mean_x=r_mean_x+Rvirt(j,1)
                    r_mean_y=r_mean_y+Rvirt(j,2)
                END DO
                  r_mean_x=r_mean_x/(last_step+1)
                  r_mean_y=r_mean_y/(last_step+1)

                rgx=0
                rgy=0

                DO j=1,last_step+1
                    rgx=rgx+(Rvirt(j,1)-r_mean_x)**2
                    rgy=rgy+(Rvirt(j,2)-r_mean_y)**2
                END DO

                rgx=rgx/(last_step+1)
                rgy=rgy/(last_step+1)
                Rgyr(tr)=sqrt(rgx+rgy)
            END DO  ! End of all virtual trajectories

            Gaus_x=0
            Gaus_y=0

            DO tr=1,nvtraj
               Rgyr_norm=log(Rgyr(tr)/(sum(Rgyr)/size(Rgyr)))
               Gaus_x=Gaus_x+Gausv_x(tr)*Rgyr_norm
               Gaus_y=Gaus_y+Gausv_y(tr)*Rgyr_norm
            END DO

            Gaus_x=Gaus_x/nvtraj
            Gaus_y=Gaus_y/nvtraj

          if (vrepultype.eq.'h') then
            call repul_hcore(np,L,sigma,Ar,R2(i,1),R2(i,2),R2,Frepx,Frepy)
          else if (vrepultype.eq.'sw') then
             call swarm(np,L,R2virt_x,R2virt_y,R2,
          *  st_virt,tau_i,extra_tau,pull)
          else if (vrepultype.eq.'s') then
             call repul_soft(np,L,sigma,Ar,R2(i,1),R2(i,2),R2,Frepx,Frepy)
          else if (vrepultype.eq.'ah')
            call repul_ahc(np,L,sigma,Ar,b,R2(i,1),R2(i,2),R2,Frepx,Frepy)
          end if

                ! Langevin equation
            Flx=-gamma*V(i,1)+noi*Gaus_x
            Fly=-gamma*V(i,2)+noi*Gaus_y

            ! Total force
            Fx=Flx+Frepx
            Fy=Fly+Frepy

            ! Integration
            Rtempx=2*R(i,1)-Rb(i,1)+Fx*dt**2
            Rtempy=2*R(i,2)-Rb(i,2)+Fy*dt**2
            V(i,1)=(Rtempx-R(i,1))/dt
            V(i,2)=(Rtempy-R(i,2))/dt
            Rb(i,1)=R(i,1)
            Rb(i,2)=R(i,2)
            R(i,1)=Rtempx
            R(i,2)=Rtempy
            R2(i,1)=MODULO(R(i,1),L)
            R2(i,2)=MODULO(R(i,2),L)
        END DO

c Analysis
            IF ((st.ge.start).and.(mod(st,incr).eq.0)) THEN

                ! Order parameters
                IF(save_ord)THEN
                call ordparam(np,R2,V,polar,apolar,tor)
                write(52,*)polar
                write(53,*)apolar
                write(54,*)tor
                END IF
                ! RDF and CCF
                call rdfccf(np,L,L2,R2,V,pi,Gr,CCF)
                DO i=1,L2 ! r instead of d
                    Gr1((st/incr),i)=Gr(i)
                    CCF1((st/incr),i)=CCF(i)
                END DO
                ! Particle number fluctuations
                call pnf(pi,L,np,R2,nl,dnl)
                DO i=1,7
                    dnl1(i)=dnl1(i)+dnl(i)
                END DO

                DO i=1,7
                    write(57,*)nl(i),dnl1(i)
                END DO
            END IF
            !----Velocity autocorrelation function----!
            IF(save_vacf)THEN
            IF(st>BL)THEN
                IF(st.le.(BL+200))THEN
                count=count+1
                DO i=1,np !
                    IF(st==BL+1)THEN ! WHY USE THE SAME THHING TWICE
                        Vmx=V(i,1)
                        Vmy=V(i,2)
                        Vm=sqrt(Vmy**2+Vmx**2)
                    END IF
                    Va=Vmx*V(i,1)+Vmy*V(i,2) ! this is the same as VV
                    V21=sqrt(V(i,1)**2+V(i,2)**2)
                    ! IF((Va.eq.0.0.and.Vm.eq.0.0).or.
     ! *   (Va.eq.0.0.and.V21.eq.0.0))THEN
                    IF(Va.eq.0.0.and.(Vm.eq.0.0.or.V21.eq.0.0))THEN
                      Vm=0.000001 ! unsafe and unsavoury
                      V21=0.000001 ! unsafe and unsavoury
                    END IF
                    Va1(count)=Va1(count)+Va/(Vm*V21)
                END DO !
                write(61,*)Va1(count)
                    IF (st.eq.BL+200) THEN
                        BL=BL+200
                        count=0
                    END IF
                END IF
            END IF
            END IF

            IF (save_coor.and.(st.ge.start).and.
     *    (mod(st,incr_coor).eq.0)) THEN
                DO i=1,np
                    write(50,*)R2(i,1),R2(i,2),V(i,1),V(i,2)
                END DO
            END IF

            ! recording temperature at every time step
            IF(save_temp) THEN
            DO i=1,np !
                skor=skor+(V(i,1)**2+V(i,2)**2)
            END DO !
            skor=skor/(np)
            IF(mod(st,incr).eq.0)write(60,*)skor
            END IF

            END DO ! end of the main loop


            ! below should be measured at steady state only
            dur=real((st-start)/incr) ! how long it took to reach steady state
            ! RDF and CCF
            DO i=1,(nsteps/incr)
                    !write(55,*)i,Gr1(i)
                    !write(56,*)i,CCF1(i)
                    write(55,'(1000F14.7)') ( Gr1(i,j), j=1,L2)
                    write(56,'(1000F14.7)') ( CCF1(i,j), j=1,L2)
        END DO
        ! Coordinates and velocities for last step
        do i=1,np
            write(49,*)R2(i,1),R2(i,2),V(i,1),V(i,2)
        end do

        ! VACF
    !    IF(save_vacf)THEN
    !    do i=1,199
    !        write(61,*)Va1(i)/(np*((st-start)/200))
    !    end do
    !    END IF

        CLOSE(49)
        IF(save_coor)CLOSE(50)
        IF(save_ord)CLOSE(52)
        IF(save_ord)CLOSE(53)
        IF(save_ord)CLOSE(54)
        CLOSE(55)
        CLOSE(56)
        CLOSE(57)
        IF(save_temp)CLOSE(60)
        IF(save_vacf)CLOSE(61)
        CLOSE(62)

        deallocate(R,Rb,R2,V,Rvirt,Gr,CCF,Gr1,CCF1)
        deallocate(Gausv_x,Gausv_y,Rgyr)
        deallocate(gnumb_x,gnumb_y)
        deallocate(Rvirt1,Rvirt2)

       END PROGRAM cognitive

!- - - - - - Generate random nonoverlapping circles - - - - - -!

        SUBROUTINE rand_nonoverlap_circle(circles,nCircles,radii,r,L)

c ### INTEGER ###
c i, j                              general int for loops
c nCircles                np        number of agents
c theormax                          theor number of agents that can be placed
c slownumb                          2/3 theormax : int(theormax*2/3)

c ### REAL ###
c radii                   sigma     size of agent
c r                       1.1sigma  minimum dist between agents for init
c L                       L         size of box
c rand                              random number generator
c x, y                              coordinates
c dx, dy                            x and y components of dist from prev circle
c tstart, tcurrent                  timing variables
c distFromPrevCircles               Euclidian between current and others
c totdistFromPrevCircles            sum of distFromPrevCircles

c ### REAL ARRAYS ###
c circles                  R        (nCircles, 2) current coordinates

c ### LOGICAL ###
c newCircleFound                    new circle found

        INTEGER i,j,nCircles,theormax,slownumb

        REAL radii,r,L,rand,x,y,dx,dy,tstart,tcurrent
        REAL distFromPrevCircles,totdistFromPrevCircles
        REAL circles(nCircles,2)

        LOGICAL newCircleFound

        ! Setup
        call cpu_time(tstart)
        theormax = (L*0.5)**2
        IF (nCircles.ge.theormax) THEN
            write(*,*)"Number of circles exceeds theoretical maximum!"
            write(*,*)nCircles,"vs",theormax
            write(*,*)"TERMINATING THE PROGRAM"
            stop
        END IF

        slownumb = int(theormax*2/3)
        IF (nCircles.ge.slownumb) THEN
            write(*,*)"WARNING: Slow regime!"
            write(*,*)nCircles,"vs",slownumb
        END IF

        circles(1:nCircles,1:2) = 0.0

        ! Main calculations
        DO i=1,nCircles
            !write(*,*)"Inserting random circle", i
            ! Flag which holds true whenever a new circle was found
            newCircleFound= .false.

            ! Loop iteration which runs until finding a circle which doesn't intersect with previous ones
            DO WHILE (newCircleFound.eq. .false.)
            totdistFromPrevCircles = 0.0
            ! Draw random coordinates of a candidate circle
            x = radii+(L-radii*2)*rand(0)
            y = radii+(L-radii*2)*rand(0)

            ! Calculate distances from previous drawn circles
            DO j=1,i
                dx = circles(j,1)-x
                dy = circles(j,2)-y
                distFromPrevCircles = sqrt(dx**2+dy**2)
                IF (distFromPrevCircles.le.real(2)*r) THEN
                    totdistFromPrevCircles = totdistFromPrevCircles
     * + distFromPrevCircles
                END IF

            END DO

            ! If the distance is not to small - add the new circle to the list
            IF ((i.eq.1).or.(totdistFromPrevCircles.eq.0.0)) THEN
                newCircleFound = .true.
                circles(i,1) = x
                circles(i,2) = y
            END IF
            call cpu_time(tcurrent)
            IF(tcurrent-tstart.gt.600)THEN
                write(*,*)"This is taking too long..."
                write(*,*)"TERMINATING THE PROGRAM"
                stop
            END IF
            END DO
        END DO

        END SUBROUTINE

!- - - - - - - - - Random seed initialisation - - - - - - - - -!

        SUBROUTINE init_random_seed()
        INTEGER*8 :: i,n,clock,seedi
        INTEGER*8, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed=abs((clock+37*(/ (i-1,i=1,n)/))**2)
        CALL RANDOM_SEED(PUT=seed)
        seedi=rand(seed)
        DEALLOCATE(seed)

        END SUBROUTINE

!- - - - - - - - - Repulsion (hard core) - - - - - - - - -!

        SUBROUTINE repul_hcore(np,L,sigma,Ar,Ri_x,Ri_y,R2,Frepx,Frepy)
        implicit none

c ### INTEGER ###
c i, j                                  general int for loops
c np            np                      number of particles

c ### REAL ###
c L             L                       size of the box
c sigma         sigma                   size of the agents
c Ar            Ar                      strength of repulsion
c Ri_x, Ri_y    R2virt_x, R2virt_y      x and y comp of modulo curr virt coord
c Frepx, Frepy  Frepx_virt, Frepy_virt  total repulsion after virtual step
c Rrepul                                repulsion distance
c Rijx, Rijy                            distance between particles (per axis)
c Rij                                   Euclidian distance between particles
c Rux, Ruy                              repulsion force modif after collision

c ### REAL ARRAYS ###
c R2            R2                      (np, 2) modulo current coordinates

        INTEGER i,j,np
        REAL L,sigma,Ar,Ri_x,Ri_y,R2(np,2),Frepx,Frepy
        REAL Rrepul,Rijx,Rijy,Rij,Rux,Ruy
        Frepx=0.0
        Frepy=0.0
        Rrepul=2*sigma

        DO j=1,np
            if(i.eq.j)continue
            Rijx=Ri_x-R2(j,1)
            Rijy=Ri_y-R2(j,2)

            if(Rijx.gt.L*0.5)then
                Rijx=Rijx-L
            else if(Rijx.lt.-L*0.5)then
                Rijx=Rijx+L
            end if

            if(Rijy.gt.L*0.5)then
                Rijy=Rijy-L
            else if(Rijy.lt.-L*0.5)then
                Rijy=Rijy+L
            end if

            Rij=sqrt(Rijx**2+Rijy**2)

            if (Rij.lt.Rrepul)then
                if (Rijx.ne.0)then
                    Rux=Rijx/Rij
                else
                    Rux=0.0
                end if
                if (Rijy.ne.0)then
                    Ruy=Rijy/Rij
                else
                    Ruy=0.0
                end if
                Frepx=Frepx+Ar*Rux
                Frepy=Frepy+Ar*Ruy
            end if
        END DO

        END SUBROUTINE

!- - - - - - - - - Swarm - - - - - - - - -!

        SUBROUTINE swarm(np,L,Ri_x,Ri_y,R2,st_virt,tau_i,extra_tau,pull)
        implicit none

c ### INTEGER ###
c i, j                                  general int for loops
c np            np                      number of particles
c st_virt       st_virt                 current step on vtraj
c tau_i         tau_i                   increased tau in swarming
c extra_tau     extra_tau               tau bonus in swarming


c ### REAL ###
c L             L                       size of the box
c pull          pull                    atraction distance modif for swarming
c Ri_x, Ri_y    R2virt_x, R2virt_y      x and y comp of modulo curr virt coord
c Frepx, Frepy  Frepx_virt, Frepy_virt  total repulsion after virtual step
c Rrepul                                repulsion distance
c Rijx, Rijy                            distance between particles (per axis)
c Rij                                   Euclidian distance between particles
c Rux, Ruy                              repulsion force modif after collision

c ### REAL ARRAYS ###
c R2            R2                      (np, 2) modulo current coordinates

        INTEGER i,j,np,st_virt,tau_i,extra_tau
        REAL L,Ri_x,Ri_y,R2(np,2),pull,Rijx,Rijy,Rij

        DO j=1,np
            if (i.eq.j) continue   ! skip to avoid self attraction

            Rijx=Ri_x-R2(j,1)
            Rijy=Ri_y-R2(j,2)

            if(Rijx.gt.L*0.5)then
                Rijx=Rijx-L
            else if(Rijx.lt.-L*0.5)then
                Rijx=Rijx+L
            end if
            if(Rijy.gt.L*0.5)then
                Rijy=Rijy-L
            else if(Rijy.lt.-L*0.5)then
                Rijy=Rijy+L
            end if

            Rij=sqrt(Rijx**2+Rijy**2)

            if (Rij.lt.pull) then        ! if within atraction distance
                tau_i=st_virt+extra_tau  ! add extra tau steps
                exit
            end if

        END DO

        Frepx_virt=0
        Frepy_virt=0

        END SUBROUTINE

!- - - - - - - - - - - Repulsion (soft) - - - - - - - - - - -!

        SUBROUTINE repul_soft(np,L,sigma,Ar,Ri_x,Ri_y,R2,Frepx,Frepy)
        implicit none

c ### INTEGER ###
c i, j                                  general int for loops
c np            np                      number of particles

c ### REAL ###
c L             L                       size of the box
c sigma         sigma                   size of the agents
c Ar            Ar                      strength of repulsion
c Ri_x, Ri_y    R2virt_x, R2virt_y      x and y comp of modulo curr virt coord
c Frepx, Frepy  Frepx_virt, Frepy_virt  total repulsion after virtual step
c Rrepul                                repulsion distance
c Rijx, Rijy                            distance between particles (per axis)
c Rij                                   Euclidian distance between particles
c Rux, Ruy                              repulsion force modif after collision
c Lin                                   adjustment for repulsion strength

c ### REAL ARRAYS ###
c R2            R2                      (np, 2) modulo current coordinates


        INTEGER i,j,np
        REAL L,sigma,Ar,Ri_x,Ri_y,R2(np,2),Frepx,Frepy
        REAL Rrepul,Rijx,Rijy,Rij,Rux,Ruy,Lin

        Frepx=0.0
        Frepy=0.0
        Rrepul=2*sigma

        DO j=1,np
            if(i.eq.j)continue
            Rijx=Ri_x-R2(j,1)
            Rijy=Ri_y-R2(j,2)

            if(Rijx.gt.L*0.5)then
                Rijx=Rijx-L
            else if(Rijx.lt.-L*0.5)then
                Rijx=Rijx+L
            end if

            if(Rijy.gt.L*0.5)then
                Rijy=Rijy-L
            else if(Rijy.lt.-L*0.5)then
                Rijy=Rijy+L
            end if

            Rij=sqrt(Rijx**2+Rijy**2)

            if (Rij.lt.Rrepul)then
                Lin=1-Rij/Rrepul

                if (Rijx.ne.0.0)then
                    Rux=Rijx/abs(Rijx)
                else
                    Rux=0.0
                end if
                if (Rijy.ne.0.0)then
                    Ruy=Rijy/abs(Rijy)
                else
                    Ruy=0.0
                end if
                Frepx=Frepx+Ar*Lin*Rux
                Frepy=Frepy+Ar*Lin*Ruy
            end if
        END DO

        END SUBROUTINE

!- - - - - - - Repulsion (adjustable hard core) - - - - - - -!

        SUBROUTINE repul_ahc(np,L,sigma,Ar,b,Ri_x,Ri_y,R2,Frepx,Frepy)
        implicit none

c ### INTEGER ###
c i, j                                  general int for loops
c np            np                      number of particles

c ### REAL ###
c L             L                       size of the box
c sigma         sigma                   size of the agents
c Ar            Ar                      strength of repulsion
c b             b                       force shape
c Ri_x, Ri_y    R2virt_x, R2virt_y      x and y comp of modulo curr virt coord
c Frepx, Frepy  Frepx_virt, Frepy_virt  total repulsion after virtual step
c Rrepul                                repulsion distance
c Rijx, Rijy                            distance between particles (per axis)
c Rij                                   Euclidian distance between particles
c Rux, Ruy                              repulsion force modif after collision
c Lin                                   adjustment for repulsion strength

c ### REAL ARRAYS ###
c R2            R2                      (np, 2) modulo current coordinates

        INTEGER i,j,np
        REAL L,sigma,Ar,b,Ri_x,Ri_y,R2(np,2),Frepx,Frepy
        REAL Rrepul,Rijx,Rijy,Rij,Rux,Ruy,Lin

        Frepx=0.0
        Frepy=0.0
        Rrepul=2*sigma

        DO j=1,np
            if(i.eq.j)continue
            Rijx=Ri_x-R2(j,1)
            Rijy=Ri_y-R2(j,2)

            if(Rijx.gt.L*0.5)then
                Rijx=Rijx-L
            else if(Rijx.lt.-L*0.5)then
                Rijx=Rijx+L
            end if

            if(Rijy.gt.L*0.5)then
                Rijy=Rijy-L
            else if(Rijy.lt.-L*0.5)then
                Rijy=Rijy+L
            end if

            Rij=sqrt(Rijx**2+Rijy**2)

            if (Rij.lt.Rrepul)then
                Lin=(exp(b*Rij/Rrepul)-exp(b))/(1-exp(b))
                if (Rijx.ne.0.0)then
                    Rux=Rijx/abs(Rijx)
                else
                    Rux=0.0
                end if
                if (Rijy.ne.0.0)then
                    Ruy=Rijy/abs(Rijy)
                else
                    Ruy=0.0
                end if
                Frepx=Frepx+Ar*Lin*Rux
                Frepy=Frepy+Ar*Lin*Ruy
            end if
        END DO

        END SUBROUTINE

!- - - - - - - - - - Order parameters - - - - - - - - - -!

        SUBROUTINE ordparam(np,R2,V,polar,apolar,tor)
        implicit none

c ### INTEGER ###
c i                                 general int for loops
c np              np                number of particles

c ### REAL ###
c polar, apolar   polar, apolar     polar and apolar order parameters
c tor             tor               angular order parameter
c speed                             Euclidian of x and y from V
c dst                               Euclidian of x and y from R2
c rot                               rotation parameter

c ### REAL ARRAYS ###
c R2              R2                (np, 2) modulo current coordinates
c V               V                 (np, 2) velocities

c ### COMPLEX ###
c comp                              complex number built from V
c ord                               parameter for polar
c ord1                              parameter for apolar

        INTEGER i,np
        REAL polar,apolar,tor,speed,dst,rot,R2(np,2),V(np,2)
        COMPLEX ord,ord1,comp

        ord=0.0
        ord1=0.0
        speed=0.0
        dst=0.0
        rot=0.0
        tor=0.0

        DO i=1,np
            !Polar and apolar order parameters
            speed=sqrt(V(i,1)**2+V(i,2)**2)
            comp=cmplx(V(i,1),V(i,2))/speed
            ord=ord+comp
            ord1=ord1+comp**2
            ! Angular order parameter
            dst=sqrt(R2(i,1)**2+R2(i,2)**2)
            rot=R2(i,1)*V(i,2)-R2(i,2)*V(i,1)
            IF(rot.ne.0)tor=tor+rot/(speed*dst)
        END DO

        ! Polar order parameter
        ord=ord/np
        polar=sqrt(REAL(ord)**2+AIMAG(ord)**2)

        ! Apolar order parameter
        ord1=ord1/np
        apolar=sqrt(REAL(ord1)**2+AIMAG(ord1)**2)

        ! Angular order parameter
        tor=abs(tor/np)

        END SUBROUTINE

!- - - - - - - - - - - RDF and CCF - - - - - - - - - - -!

        SUBROUTINE rdfccf(np,L,L2,R2,V,pi,Gr,CCF)
        implicit none

c ### INTEGER ###
c i, j                      general int for loops
c np              np        number of particles
c nhis                      total number of bins
c delg                      bin size

c ### INTEGER ARRAYS ###
c NN1                       (L2) bins keeping tabs of agent number

c ### REAL ###
c L               L         size of the box
c pi                        pi constant as pi=4.0d0*atan(1.0d0)
c dx, dy                    distance in x and y
c dxy                       Euclidian of dx and dy
c VV                        velocity correlation
c Sfs                       total surface area
c Rhos                      general agent density
c SS                        area difference of bins
c Rho                       agent density in bin

c ### REAL ARRAYS ###
c R2              R2        (np, 2) modulo current coordinates
c V               V         (np, 2) velocities
c Gr              Gr        (L2) current RDF
c CCF             CCF       (L2) curretn CCF
c VCOR                      (L2) bins keeping tabs of VV

        INTEGER i,j,np,L2,nhis,delg,,NN1(L2)
        REAL L,R2(np,2),V(np,2),pi,Gr(L2),CCF(L2)
        REAL VCOR(L2),dx,dy,dxy,VV,Sfs,Rhos,SS,Rho

        nhis=L2 !Total number of bins
        delg=int(L)/(2*nhis) !Bin size

        Gr(1:nhis)=0.0
        CCF(1:nhis)=0.0

        NN1(1:nhis)=0
        VCOR(1:nhis)=0

        DO i=1,np
            DO j=1,np
                if(i.eq.j)continue ! if different agents
                dx=R2(i,1)-R2(j,1) ! distance in x
                dy=R2(i,2)-R2(j,2) ! distance in y
                dx=dx-L*nint(dx/int(L)) ! if dx > 0.5L, use other dx
                dy=dy-L*nint(dy/int(L)) ! if dy > 0.5L, use other dy
                dxy=(sqrt(dx**2+dy**2)) ! Euclidian distance
                VV=V(i,1)*V(j,1)+V(i,2)*V(j,2) ! Sum of product of V components
                IF(dxy<L2)THEN ! if distance less than half box (should be)
                    dxy=int(dxy/delg) ! pick destination box from distance
                    NN1(dxy)=NN1(dxy)+1 ! +1 particle in corresponding bin
                    VCOR(dxy)=VCOR(dxy)+VV ! add VV to corresponding VCOR value
                END IF
            END DO
        END DO

        Sfs=L**2
        Rhos=np/Sfs

        DO i=1,nhis
          SS=pi*real((i+1)**2-i**2)*real(delg**2) ! Marta says increase bin area
          Rho=real(NN1(i))/SS
          Gr(i)=Rho/real(Rhos*np)
          CCF(i)=VCOR(i)/real(NN1(i))
          IF((VCOR(i).eq.0).and.(NN1(i).eq.0))THEN
              CCF(i) = 0.0
          END IF
        END DO

        END SUBROUTINE

!- - - - - - - - Particle number fluctuations - - - - - - - -!

        SUBROUTINE pnf(pi,L,np,R2,nl,dnl)
        Implicit none
c ### INTEGER ###
c i, j                      general int for loops
c np              np        number of particles
c cellsize                  cells edge lenght
c rn                        cells edge lenght INT ?????????
c ncel                      no of cells per grid row
c icel                      cell index
c totalncel                 total no of cells
c icelx, icely              grid cell index per axis

c ### INTEGER ARRAYS ###
c sizes(7)                  grid cell edge lenght

c ### REAL ###
c L               L         size of the box
c pi                        pi constant as pi=4.0d0*atan(1.0d0)
c srkv                      distribution evenness
c kvsr                      should be np

c ### REAL ARRAYS ###
c R2                        (np, 2) modulo current coordinates
c nl(7)                     pnf mean no of particles per cell
c dnl(7)                    pnf equation (SD of particle no)

        INTEGER i,j,np,cellsize,rn,ncel,icel
        INTEGER totalncel,icelx,icely,sizes(7)
        REAL L,R2(np,2),nu1,pi,nl(7),nl1,dnl(7),srkv,kvsr
        INTEGER, DIMENSION(:), ALLOCATABLE :: N

        sizes(1)=20
        sizes(2)=10
        sizes(3)=8
        sizes(4)=5
        sizes(5)=4
        sizes(6)=2
        sizes(7)=1

        nl(1:size(sizes,1))=0.0
        dnl(1:size(sizes,1))=0.0

        DO j=1,size(sizes,1)

            rn=int(L)/int(L/sizes(j)) ! grid cell size from INTS
            ncel=int(L)/rn ! number of cells in row
            totalncel=ncel**2 ! number of cells

            srkv=0.0
            kvsr=0.0

            allocate(N(0:totalncel))
            N(0:totalncel)=0 ! initialise N for cells

            do i=1,np
                icelx=int(R2(i,1)/rn) ! find cell for x
                icely=int(R2(i,2)/rn) ! find cell for y
                icel=(icely)*ncel+icelx ! find cell index
                N(icel)=N(icel)+1 ! increment
            end do

            do icel=0,totalncel-1
                srkv=srkv+N(icel)**2 ! increase if higher density than expected
                kvsr=kvsr+N(icel) ! should increase with np
            end do

            deallocate(N)

            srkv=srkv/totalncel
            nl(j)=kvsr/totalncel
            nl1=nl(j)**2
            dnl(j)=sqrt(srkv-nl1)

        END DO

        END SUBROUTINE
