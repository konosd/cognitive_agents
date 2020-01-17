c Cognitive map model
c (c) Maksym Romenskyy
c  22 Dec 2018
c Hello World


       PROGRAM cognitive
         implicit none

        INTEGER i,j,np,st,nsteps,tr,nvtraj,st_virt,tau,L2
        INTEGER start,incr,incr_coor,count,NBIN,BL,seqn

        REAL sigma,rand,L,noi,gamma,temp,dt,Ar,pi,Fx,Flx,Fy,Fly
        REAL Rtempx,Rtempy,Gaus_x,Gaus_y,Rgyr_norm,rgx,rgy,dnl(7)
        REAL r_mean_x,r_mean_y,R2virt_x,R2virt_y,Rbx,Rby,nl(7)
        REAL Vvirt_x,Vvirt_y,Frepx_virt,Frepy_virt,Frepx,Frepy,skor
        REAL Rijx,Rijy,Rij,Rux,Ruy,b,polar,apolar,tor,dur,dnl1(7)
        REAL Vmx,Vmy,Va,Va1(200),Vm,V21

        REAL, DIMENSION(:,:), ALLOCATABLE :: R,Rb,R2,V,Rvirt
        REAL, DIMENSION(:), ALLOCATABLE :: gnumb_x,gnumb_y
        REAL, DIMENSION(:), ALLOCATABLE :: Rvirt1,Rvirt2,Gr,CCF
        REAL, DIMENSION(:), ALLOCATABLE :: Gr1,CCF1
        REAL, DIMENSION(:), ALLOCATABLE :: Gausv_x,Gausv_y,Rgyr

        CHARACTER(len=90)::name,nvtraj1,tau1,sigma1,L1,a1,seqn1
        CHARACTER(len=90)::np1,repultype,InfoFile,CoorFile,CoorlFile
        CHARACTER(len=90)::PolFile,ApolFile,AngFile,RDFFile,VACFFile
        CHARACTER(len=90)::CCFFile,PNFFile,TempFile

        LOGICAL save_temp,save_coor,save_ord,save_vacf

c Simulation parameters
        ! Real trajectrories
        np = 100            ! number of real particles
        nsteps = 2000       ! number of real steps

        ! Virtual trajectrories
        nvtraj = 20         ! number of virtual trajectories
        tau = 20            ! number of virtual steps

        ! General
        sigma = 1.0         ! particle size
        Ar=20.0             ! force strength
        b = 200.0            ! force shape
        L = 80.0*sigma      ! box size
        dt = 0.1            ! time step

        ! Thermostat
        gamma = 1.0
        temp = 1.0

c Analysis parameters

        seqn=1           ! file in sequence
        start = 1        ! starting timestep for analysis
        incr = 10           ! analyse every incr timestep
        incr_coor = 50     ! save coor and vel every incr timestep
        NBIN = 100

c Types of Analysis

        save_temp=.true.        ! compute temperature
        save_coor=.true.        ! write coordinates every incr_coor step
        save_ord=.true.        ! compute order parameters
        save_vacf=.true.       ! compute VACF

c Setup
        ! Create unique filename with simulation parameters encoded
        repultype="hs"      ! h - hard, ah - adjust hard, s - soft
        write(np1,*)np
        write(nvtraj1,*)nvtraj
        write(tau1,*)tau
        write(sigma1,*)int(sigma)
        write(L1,*)int(L)
        write(a1,*)int(Ar)
        write(seqn1,*)seqn
        name=trim(adjustl("n"))//trim(adjustl(np1))//
     *  trim(adjustl("nvt"))//trim(adjustl(nvtraj1))//
     *  trim(adjustl("tau"))//trim(adjustl(tau1))//
     *  trim(adjustl("s"))//trim(adjustl(sigma1))//
     *  trim(adjustl("L"))//trim(adjustl(L1))//
     *  trim(adjustl("a"))//trim(adjustl(a1))//
     *  trim(adjustl(repultype))//trim(adjustl("_"))//
     *  trim(adjustl(seqn1))

        ! Write info not encoded in filename into a *.txt file
        write(InfoFile,'(a,a)')trim(name),"_info.txt"
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
        IF(save_coor)write(CoorFile,'(a,a)')trim(name),"_coor.dat"
        write(CoorlFile,'(a,a)')trim(name),"_coorl.dat"
        IF(save_ord)write(PolFile,'(a,a)')trim(name),"_pol.dat"
        IF(save_ord)write(ApolFile,'(a,a)')trim(name),"_apol.dat"
        IF(save_ord)write(AngFile,'(a,a)')trim(name),"_ang.dat"
        write(RDFFile,'(a,a)')trim(name),"_rdf.dat"
        write(CCFFile,'(a,a)')trim(name),"_ccf.dat"
        write(PNFFile,'(a,a)')trim(name),"_pnf.dat"
        IF(save_temp)write(TempFile,'(a,a)')trim(name),"_temp.dat"
        IF(save_vacf)write(VACFFile,'(a,a)')trim(name),"_VACF.dat"

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


        ! Memory and randomization setup
        call init_random_seed()
        allocate(R(np,2))
        allocate(Rb(np,2))
        allocate(R2(np,2))
        allocate(V(np,2))
        allocate(Gausv_x(nvtraj))
        allocate(Gausv_y(nvtraj))
        allocate(Rgyr(nvtraj))
        allocate(Rvirt(tau+1,2))
        allocate(gnumb_x(tau))
        allocate(gnumb_y(tau))
        allocate(Rvirt1(tau))
        allocate(Rvirt2(tau))
        L2=int(L)*0.5
        allocate(Gr(L2))
        allocate(CCF(L2))
        allocate(Gr1(L2))
        allocate(CCF1(L2))
        pi=4.0d0*ATAN(1.0d0)
        count=0
        Gr1(1:L2)=0.0
        CCF1(1:L2)=0.0
        dnl1(1:7)=0.0
        skor=0.0
        BL=start
        Va1=0.0

c Initial conditions
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
        DO st=1,nsteps
            DO i=1,np
                ! Virtual trajectory for particle i
                DO tr=1,nvtraj
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
        DO st_virt=1,tau
            gnumb_x(st_virt)=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
            gnumb_y(st_virt)=sqrt(-2.0*log(rand(0)))*sin(2*pi*rand(0))
        END DO
            ! Save the first-step-random-numbers
            Gausv_x(tr)=gnumb_x(1)
            Gausv_y(tr)=gnumb_y(1)

                    ! Virtual steps
                    DO st_virt=1,tau

        ! Virtual repulsion (hard core)
        call repul_hcore(np,L,sigma,Ar,R2virt_x,R2virt_y,R2,
     *  Frepx_virt,Frepy_virt)
        ! Virtual repulsion (soft)
c        call repul_soft(np,L,sigma,Ar,R2virt_x,R2virt_y,R2,
c     *  Frepx_virt,Frepy_virt)
        ! Virtual repulsion (adjustable hard core)
c        call repul_ahc(np,L,sigma,Ar,b,R2virt_x,R2virt_y,R2,
c     *  Frepx_virt,Frepy_virt)

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
                    END DO

                    ! Radius of gyration of trajectory i
                    Rvirt1(:)=Rvirt(:,1)
                    Rvirt2(:)=Rvirt(:,2)
                    r_mean_x=sum(Rvirt1)/size(Rvirt1)
                    r_mean_y=sum(Rvirt2)/size(Rvirt2)
                    rgx=0
                    rgy=0

                    DO st_virt=1,tau+1
                        rgx=rgx+(Rvirt(st_virt,1)-r_mean_x)**2
                        rgy=rgy+(Rvirt(st_virt,2)-r_mean_y)**2
                    END DO

                    rgx=rgx/size(Rvirt,1)
                    rgy=rgy/size(Rvirt,1)
                    Rgyr(tr)=sqrt(rgx+rgy)
                END DO

                Gaus_x=0
                Gaus_y=0

                DO tr=1,nvtraj
                    Rgyr_norm=log(Rgyr(tr)/(sum(Rgyr)/size(Rgyr)))
                    Gaus_x=Gaus_x+Gausv_x(tr)*Rgyr_norm
                    Gaus_y=Gaus_y+Gausv_y(tr)*Rgyr_norm
                END DO

                Gaus_x=Gaus_x/nvtraj
                Gaus_y=Gaus_y/nvtraj

                ! Langevin equation
                Flx=-gamma*V(i,1)+noi*Gaus_x
                Fly=-gamma*V(i,2)+noi*Gaus_y

        ! Repulsion (hard core)
c        call repul_hcore(np,L,sigma,Ar,R2(i,1),R2(i,2),R2,Frepx,Frepy)
        ! Repulsion (soft)
        call repul_soft(np,L,sigma,Ar,R2(i,1),R2(i,2),R2,Frepx,Frepy)
        ! Repulsion (adjustable hard core)
c        call repul_ahc(np,L,sigma,Ar,b,R2(i,1),R2(i,2),R2,Frepx,Frepy)

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
                DO i=1,np
                    IF(save_temp)skor=skor+(V(i,1)**2+V(i,2)**2)
                END DO
                ! Order parameters
                IF(save_ord)THEN
                call ordparam(np,R2,V,polar,apolar,tor)
                write(52,*)polar
                write(53,*)apolar
                write(54,*)tor
                END IF
                ! RDF and CCF
                call rdfccf(np,L,L2,R2,V,pi,Gr,CCF)
                DO i=1,L2
                    Gr1(i)=Gr1(i)+Gr(i)
                    CCF1(i)=CCF1(i)+CCF(i)
                END DO
                ! Particle number fluctuations
                call pnf(pi,L,np,R2,nl,dnl)
                DO i=1,7
                    dnl1(i)=dnl1(i)+dnl(i)
                END DO
            END IF
            !----Velocity autocorrelation function----!
            IF(save_vacf)THEN
            IF(st>BL)THEN
                IF(st<BL+200)THEN
                count=count+1
                DO i=1,np
                    IF(st==BL+1)THEN
                        Vmx=V(i,1)
                        Vmy=V(i,2)
                        Vm=sqrt(Vmy**2+Vmx**2)
                    END IF
                    Va=Vmx*V(i,1)+Vmy*V(i,2)
                    V21=sqrt(V(i,1)**2+V(i,2)**2)
                    IF((Va.eq.0.0.and.Vm.eq.0.0).or.
     *    (Va.eq.0.0.and.V21.eq.0.0))THEN
                    Vm=0.000001
                    V21=0.000001
                    END IF
                    Va1(count)=Va1(count)+Va/(Vm*V21)
                END DO
                ELSE IF (st.eq.BL+200) THEN
                    BL=BL+200
                    count=0
                END IF
            END IF
            END IF

            IF (save_coor.and.(st.ge.start).and.
     *    (mod(st,incr_coor).eq.0)) THEN
                DO i=1,np
                    write(50,*)R2(i,1),R2(i,2),V(i,1),V(i,2)
                END DO
            END IF

        END DO

        dur=real((st-start)/incr)
        ! RDF and CCF
        DO i=1,L2
                write(55,*)i,Gr1(i)/dur
                write(56,*)i,CCF1(i)/dur
        END DO
        ! Particle number fluctuations
        DO i=1,7
            write(57,*)nl(i),dnl1(i)/dur
        END DO
        ! Coordinates and velocities for last step
        do i=1,np
            write(49,*)R2(i,1),R2(i,2),V(i,1),V(i,2)
        end do
        IF(save_temp)skor=skor/(np*dur)
        IF(save_temp)write(60,*)skor
        ! VACF
        IF(save_vacf)THEN
        do i=1,199
            write(61,*)Va1(i)/(np*((st-start)/200))
        end do
        END IF

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

        deallocate(R,Rb,R2,V,Rvirt,Gr,CCF,Gr1,CCF1)
        deallocate(Gausv_x,Gausv_y,Rgyr,gnumb_x,gnumb_y)
        deallocate(Rvirt1,Rvirt2)

       END PROGRAM cognitive

!- - - - - - Generate random nonoverlapping circles - - - - - -!

        SUBROUTINE rand_nonoverlap_circle(circles,nCircles,radii,r,L)

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
        IF (nCircles.ge. slownumb) THEN
            write(*,*)"WARNING: Slow regime!"
            write(*,*)nCircles,"vs",slownumb
        END IF

        circles(1:nCircles,1:2) = 0.0

        ! Main calculations
        DO i=1,nCircles
            !write(*,*)"Inserting random circle", i
            ! Flag which holds true whenever a new circle was found
            newCircleFound = .false.

            ! Loop iteration which runs until finding a circle which doesn't intersect with previous ones
            DO WHILE (newCircleFound .eq. .false.)
            totdistFromPrevCircles = 0.0
            ! Draw random coordinates of a candidate circle
            x = radii + (L-radii*2)*rand(0)
            y = radii + (L-radii*2)*rand(0)

            ! Calculate distances from previous drawn circles
            DO j=1,i
                dx = circles(j,1)-x
                dy = circles(j,2)-y
                distFromPrevCircles = sqrt(dx**2+dy**2)
                IF (distFromPrevCircles .le. real(2)*r) THEN
                    totdistFromPrevCircles = totdistFromPrevCircles
     * + distFromPrevCircles
                END IF

            END DO

            ! If the distance is not to small - add the new circle to the list
            IF (i .eq. 1 .or. totdistFromPrevCircles .eq. 0.0) THEN
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
        INTEGER*8 :: i, n, clock,seedi
        INTEGER*8, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed=abs((clock+37*(/ (i - 1, i = 1, n) /))**2)
        CALL RANDOM_SEED(PUT = seed)
        seedi=rand(seed)
        DEALLOCATE(seed)

        END SUBROUTINE

!- - - - - - - - - Repulsion (hard core) - - - - - - - - -!

        SUBROUTINE repul_hcore(np,L,sigma,Ar,Ri_x,Ri_y,R2,Frepx,Frepy)
        implicit none

        INTEGER i,j,np
        REAL Rijx,Rijy,Rij,Ri_x,Ri_y,R2(np,2),Rux,Ruy
        REAL sigma,Ar,Frepx,Frepy,Rrepul,L

        Frepx=0.0
        Frepy=0.0
        Rrepul=2*sigma

        DO j=1,np
            if (i .eq. j) cycle
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
            if (Rij .lt. Rrepul)then
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

!- - - - - - - - - - - Repulsion (soft) - - - - - - - - - - -!

        SUBROUTINE repul_soft(np,L,sigma,Ar,Ri_x,Ri_y,R2,Frepx,Frepy)
        implicit none

        INTEGER i,j,np
        REAL Rijx,Rijy,Rij,Ri_x,Ri_y,R2(np,2),Rux,Ruy
        REAL sigma,Ar,L,Frepx,Frepy,Lin,Rrepul

        Frepx=0.0
        Frepy=0.0
        Rrepul=2*sigma

        DO j=1,np
            if (i .eq. j) cycle
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

            if (Rij .lt. Rrepul)then
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

        INTEGER i,j,np
        REAL Rijx,Rijy,Rij,Ri_x,Ri_y,R2(np,2),Rux,Ruy
        REAL sigma,Ar,b,L,Frepx,Frepy,Lin,Rrepul

        Frepx=0.0
        Frepy=0.0
        Rrepul=2*sigma

        DO j=1,np
            if (i .eq. j) cycle
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

            if (Rij .lt. Rrepul)then
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

        INTEGER i,np
        REAL polar,apolar,TOR,speed,DST,rot,R2(np,2),V(np,2)
        COMPLEX ord,ord1,comp

        ord=0.0
        ord1=0.0
        speed=0.0
        DST=0.0
        rot=0.0
        TOR=0.0

        DO i=1,np
            !Polar and apolar order parameters
            speed=sqrt(V(i,1)**2+V(i,2)**2)
            comp=CMPLX(V(i,1),V(i,2))/speed
            ord=ord+comp
            ord1=ord1+comp**2
            ! Angular order parameter
            DST=SQRT(R2(i,1)**2+R2(i,2)**2)
            rot=R2(i,1)*V(i,2)-R2(i,2)*V(i,1)
            IF(rot.ne.0)TOR=TOR+rot/(speed*DST)
        END DO

        ! Polar order parameter
        ord=ord/np
        polar=SQRT(REAL(ord)**2+AIMAG(ord)**2)

        ! Apolar order parameter
        ord1=ord1/np
        apolar=SQRT(REAL(ord1)**2+AIMAG(ord1)**2)

        ! Angular order parameter
        TOR=abs(TOR/np)

        END SUBROUTINE

!- - - - - - - - - - - RDF and CCF - - - - - - - - - - -!

        SUBROUTINE rdfccf(np,L,L2,R2,V,pi,Gr,CCF)
        implicit none

        INTEGER i,j,np,L2,NN1(L2),nhis,delg
        REAL R2(np,2),V(np,2),L,Gr(L2),CCF(L2)
        REAL dx,dy,dxy,VV,VCOR(L2),Sfs,Rhos,SS,Rho,pi

        nhis=L2 !Total number of bins
        delg=int(L)/(2*nhis) !Bin size

        Gr(1:nhis)=0.0
        CCF(1:nhis)=0.0

        NN1(1:nhis)=0
        VCOR(1:nhis)=0

        DO i=1,np
            DO j=1,np
                if(i .eq. j)cycle
                dx=R2(i,1)-R2(j,1)
                dy=R2(i,2)-R2(j,2)
                dx=dx-L*NINT(dx/int(L))
                dy=dy-L*NINT(dy/int(L))
                dxy=(SQRT(dx**2+dy**2))
                VV=V(i,1)*V(j,1)+V(i,2)*V(j,2)
                IF(dxy<L2)THEN
                    dxy=int(dxy/delg)
                    NN1(dxy)=NN1(dxy)+1 !Number of particles in each circle
                    VCOR(dxy)=VCOR(dxy)+VV
                END IF
            END DO
        END DO

        Sfs=L**2
        Rhos=np/Sfs

        DO i=1,nhis
            SS=pi*real((i+1)**2-i**2)*real(delg**2)
            Rho=real(NN1(i))/SS
            Gr(i)=Rho/real(Rhos*np)
            CCF(i)=VCOR(i)/real(NN1(i))
            !write(74,*)i,Gr(i)
            !write(75,*)i,CCF(i)
        END DO

        END SUBROUTINE

!- - - - - - - - Particle number fluctuations - - - - - - - -!

        SUBROUTINE pnf(pi,L,np,R2,nl,dnl)
        Implicit none

        INTEGER rn,cellsize,icel,sizes(7)
        INTEGER ncel,np,i,j,totalncel,icelx,icely
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

            cellsize=sizes(j)

            rn=int(L)/int(L/cellsize)
            ncel=int(L)/rn
            totalncel=ncel**2

            srkv=0.0
            kvsr=0.0

            allocate(N(0:totalncel))
            N(0:totalncel)=0

            do i=1,np
                icelx=int(R2(i,1)/rn)
                icely=int(R2(i,2)/rn)
                icel=(icely)*ncel+icelx
                N(icel)=N(icel)+1
            end do

            do icel=0,totalncel-1
                srkv=srkv+N(icel)**2
                kvsr=kvsr+N(icel)
            end do

            deallocate(N)

            srkv=srkv/totalncel
            nl(j)=kvsr/totalncel
            nl1=nl(j)**2
            dnl(j)=SQRT(srkv-nl1)

        END DO

        END SUBROUTINE
