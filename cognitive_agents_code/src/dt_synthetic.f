c Cognitive map model
c (c) Maksym Romenskyy
c  22 Dec 2018

       PROGRAM cognitive
         use Ziggurat
         implicit none

        INTEGER i,j,np,st,nsteps,tr,nvtraj,st_virt,tau,L2,vrr
        INTEGER start,incr_coor,count,NBIN,seqn,t_precise
        !INTEGER incr,BL


        REAL sigma,rand,L,noi,gamma,temp,dt,Ar,Fx,Flx,Fy,Fly
        REAL Rtempx,Rtempy,Gaus_x,Gaus_y,Rgyr_norm,rgx,rgy!,dnl(7)
        REAL r_mean_x,r_mean_y,R2virt_x,R2virt_y,Rbx,Rby!,nl(7)
        REAL Vvirt_x,Vvirt_y,Frepx_virt,Frepy_virt,Frepx,Frepy,skor
        REAL Rijx,Rijy,Rij,Rux,Ruy,b





      
        REAL T, t_count, dtvirt, frac, Rv_x, Rv_y



        REAL inidist




        REAL, DIMENSION(:,:), ALLOCATABLE :: R,Rb,R2,V,Rvirt





        
        REAL, DIMENSION(:,:), ALLOCATABLE :: vrecord,rg
        
        
        
        
        
        
        REAL, DIMENSION(:), ALLOCATABLE :: gnumb_x,gnumb_y
        REAL, DIMENSION(:), ALLOCATABLE :: Rvirt1,Rvirt2,Gr,CCF
        REAL, DIMENSION(:), ALLOCATABLE :: Gr1,CCF1
        REAL, DIMENSION(:), ALLOCATABLE :: Gausv_x,Gausv_y,Rgyr

        REAL, DIMENSION(4):: r4
        REAL(4)::pi
        CHARACTER(len=90)::name,nvtraj1,tau1,sigma1,L1,a1,seqn1,nsteps1
        CHARACTER(len=90)::np1,repultype,InfoFile,CoorFile    
        !  CHARACTER(len=90)::TempFile,CoorlFile


        CHARACTER(len=90)::vrecordfile, rgfile, forcefile
        
        
        
        
        CHARACTER(len=90)::tau_char,np_char,nsteps_char,file_nr,
     *  inidist_char, frac_char
        CHARACTER(len=90)::dir_name,file_name,path

      !  LOGICAL save_temp,save_ord,save_vacf
        LOGICAL d_exists,save_coor,testing




        ! ###SYNTHETIC###
        ! FIX THE SEEDS
c       integer, allocatable :: seed(:)
c       integer seed_size
c       call random_seed(size=seed_size)
c       allocate(seed(seed_size))
c       seed(:) = 10
c       call random_seed(put=seed)
c       call zigset(10)
c       deallocate(seed)
c getting the command line arguments

                CALL GET_COMMAND_ARGUMENT(1,np_char)   !load the command line arguments
                CALL GET_COMMAND_ARGUMENT(2,tau_char)
                CALL GET_COMMAND_ARGUMENT(3,nsteps_char)
                CALL GET_COMMAND_ARGUMENT(4,file_nr)
                CALL GET_COMMAND_ARGUMENT(5,frac_char)
                CALL GET_COMMAND_ARGUMENT(6,inidist_char)
                READ(np_char,*)np      !then, convert them from CHAR (default) to INTEGER
                READ(tau_char,*)tau
                READ(nsteps_char,*)nsteps
                READ(frac_char,*)frac
                READ(inidist_char,*)inidist


c Simulation parameters




        ! For changing dt
        dt=0.1
        T=tau*dt
c       frac=1







        
        vrr=1  ! virtual traj record particle number (<np)
        ! Real trajectrories
        !np = 200            ! number of real particles
        !nsteps = 10000       ! number of real steps
        !file_nr = "1"

        ! Virtual trajectrories
        nvtraj = 360         ! number of virtual trajectories
        !tau = 20            ! number of virtual steps

        ! General
        sigma = 1.0         ! particle size
        Ar=20.0             ! force strength
        b = 200.0            ! force shape
        L = 80.0*sigma      ! box size

        ! Thermostat
        gamma = 1.0
        temp = 1.0

        testing = .false.

c Analysis parameters

        seqn=1           ! file in sequence
        start = 1        ! starting timestep for analysis
        !incr = 1           ! analyse every incr timestep
        incr_coor = 1     ! save coor and vel every incr timestep
        NBIN = 100

c Types of Analysis

    !    save_temp=.true.        ! compute temperature
        save_coor=.true.        ! write coordinates every incr_coor step
    !    save_ord=.true.        ! compute order parameters
    !    save_vacf=.true.       ! compute VACF

c Setup
        ! Create unique filename with simulation parameters encoded
        repultype="S"      ! h - hard, a - adjust hard, s - soft
        write(np1,*)np
        write(nvtraj1,*)nvtraj
        write(tau1,*)tau
        !write(sigma1,*)int(sigma)
        write(L1,*)int(L)
        !write(a1,*)int(Ar)
        !write(seqn1,*)seqn
        write(nsteps1,*)nsteps

        dir_name=trim(adjustl("na"))//trim(adjustl(np1))//
     *  trim(adjustl("_vs"))//trim(adjustl(tau1))//
     *  trim(adjustl("_st"))//trim(adjustl(nsteps1))//  !!!nsteps_char??
     *  trim(adjustl("_lb"))//trim(adjustl(L1))//
     *  trim(adjustl("_rt"))//trim(adjustl(repultype))

        INQUIRE(FILE=dir_name, EXIST=d_exists)

        if (.not.d_exists) then
          call system('mkdir '//dir_name)
        end if

        ! Write info not encoded in filename into a *.txt file
        write(InfoFile,'(a,a)')trim(dir_name),"_info.txt"
        OPEN(51,FILE=InfoFile,form='formatted')
        write(51,*)"nsteps",nsteps
        write(51,*)"dt",dt
        write(51,*)"gamma",gamma
        write(51,*)"temp",temp
        write(51,*)"start",start
        !write(51,*)"incr",incr
        write(51,*)"incr_coor",incr_coor
        CLOSE(51)

        ! Output setup
        file_name = trim(adjustl(file_nr))
        path = trim(dir_name)//'/'//trim(file_name)

        IF(save_coor)write(CoorFile,'(a,a)')trim(path),"_coor.dat"
        write(vrecordfile,'(a,a)')trim(path),"_vrecord.dat"
        write(rgfile,'(a,a)')trim(path),"_rg.dat"
        write(forcefile,'(a,a)')trim(path),"_force.dat"





    !    IF(save_temp)write(TempFile,'(a,a)')trim(name),"_temp.dat"

    !    OPEN(49,FILE=CoorlFile,form='formatted')
        IF(save_coor)OPEN(50,FILE=CoorFile,form='formatted')
        OPEN(49,FILE=vrecordfile,form='formatted')
        OPEN(48,FILE=rgfile,form='formatted')
        OPEN(47,FILE=forcefile,form='formatted')





    !    IF(save_temp)OPEN(60,FILE=TempFile,form='formatted')


        ! Memory and randomization setup
        call init_random_seed()
        allocate(R(np,2))
        allocate(Rb(np,2))
        allocate(R2(np,2))
        allocate(V(np,2))
        allocate(Gausv_x(nvtraj))
        allocate(Gausv_y(nvtraj))
        allocate(Rgyr(nvtraj))
        allocate(Rvirt(tau,2))
        allocate(gnumb_x(tau))
        allocate(gnumb_y(tau))
        allocate(Rvirt1(tau))
        allocate(Rvirt2(tau))

        allocate(vrecord(nvtraj/18,tau*2))
        allocate(rg(nvtraj,2))
        L2=int(L)*0.5
        pi=4.0d0*ATAN(1.0d0)

c Initial conditions
        ! call rand_nonoverlap_circle(R,np,sigma,sigma+sigma/10,L)
        R(:,1)=(/40.00-inidist,42.00,42.00,42.00,42.00,42.00,
     * 42.00,42.00/)
        R(:,2)=(/40.00,34.00,36.00,38.00,40.00,42.00,
     * 44.00,46.00/)
        noi=sqrt(2.0*gamma*temp/dt)
        Rb(:,1)=R(:,1)
        Rb(:,2)=R(:,2)
        R2(:,1)=R(:,1)
        R2(:,2)=R(:,2)
        DO i=1,np
        V(i,1)=rnor( )
        V(i,2)=rnor( )
        END DO

c Main loop
        DO st=1,nsteps
        rg(:,2)=0



            DO i=1,1
                ! Virtual trajectory for particle i
                IF (.not.testing)THEN
                DO tr=1,nvtraj
                    ! Setup for the first virtual step
                    ! Rvirt(1,1)=R(i,1)
                    ! Rvirt(1,2)=R(i,2)
                    R2virt_x=R2(i,1)
                    R2virt_y=R2(i,2)
                    Rbx=Rb(i,1)
                    Rby=Rb(i,2)
                    Vvirt_x=V(i,1)
                    Vvirt_y=V(i,2)






                    t_count=0
                    Rv_x=R(i,1)
                    Rv_y=R(i,2)
                    st_virt=floor(t_count/dt)





                    ! Generate random numbers for the whole virtual trajectory
                    ! Save the first-step-random-numbers
                    Gausv_x(tr)=rnor()
                    Gausv_y(tr)=rnor()

                    if ((Gausv_x(tr).gt.0).and.(i.eq.vrr)) then
                        rg(tr,2)=1
                    end if

                    ! Virtual Langevin equation
                    Flx=-gamma*Vvirt_x+noi*Gausv_x(tr)
                    Fly=-gamma*Vvirt_y+noi*Gausv_y(tr)
                    




                    ! Virtual steps
                    DO WHILE (t_count.lt.(T-0.5*dt*frac))

        ! Virtual repulsion (hard core)
        call repul_soft(i,np,L,sigma,Ar,R2virt_x,R2virt_y,R2,
     *  Frepx_virt,Frepy_virt)
        ! Virtual repulsion (soft)
!        call repul_soft(np,L,sigma,Ar,R2virt_x,R2virt_y,R2,
!     *  Frepx_virt,Frepy_virt)
        ! Virtual repulsion (adjustable hard core)
!        call repul_ahc(np,L,sigma,Ar,b,R2virt_x,R2virt_y,R2,
!     *  Frepx_virt,Frepy_virt)


                        ! Virtual total force
                        Fx=Flx+Frepx_virt
                        Fy=Fly+Frepy_virt





                        ! Distinguish encountering another agent or not
                        IF ((Frepx_virt.eq.0).and.
     *(Frepy_virt.eq.0)) THEN
                            dtvirt=dt
                        ELSE 
                            dtvirt=frac*dt
                        END IF 






                        ! Virtual integration
c                       Rtempx=2*Rv_x-Rbx+Fx*dtvirt**2
c                       Rtempy=2*Rv_y-Rby+Fy*dtvirt**2
c                       Vvirt_x=(Rtempx-Rv_x)/dtvirt
c                       Vvirt_y=(Rtempy-Rv_y)/dtvirt
                        
                        
                        
                        
                        Vvirt_x=Vvirt_x+Fx*dtvirt
                        Vvirt_y=Vvirt_y+Fy*dtvirt
                        Rtempx=Rv_x+Vvirt_x*dtvirt
                        Rtempy=Rv_y+Vvirt_y*dtvirt



                        Rbx=Rv_x
                        Rby=Rv_y
                        Rv_x=Rtempx
                        Rv_y=Rtempy
                        R2virt_x=MODULO(Rv_x,L)
                        R2virt_y=MODULO(Rv_y,L)




                        ! Update tau counting and record virtual
                        ! coordinates for radius of gyration
                        t_count=t_count+dtvirt
                        t_precise=nint(t_count/(dt*frac))
                        IF (floor(t_precise*frac).ne.st_virt) THEN
                            st_virt=floor(t_precise*frac)
                            Rvirt(st_virt,1)=Rv_x
                            Rvirt(st_virt,2)=Rv_y
                        END IF

                            IF ((st.eq.1).and.(i.eq.vrr).and.(tr.eq.2))
     *                      THEN
                            print *,Rv_x, Rv_y ,  dtvirt
     *                     ,floor(t_precise*frac),t_count
                            END IF

                        ! Langevin equation
                        Flx=-gamma*Vvirt_x+noi*rnor()
                        Fly=-gamma*Vvirt_y+noi*rnor()
                    END DO











                    ! Radius of gyration of trajectory i
                    Rvirt1(:)=Rvirt(:,1)
                    Rvirt2(:)=Rvirt(:,2)
                    r_mean_x=sum(Rvirt1)/(tau+1)
                    r_mean_y=sum(Rvirt2)/(tau+1)
                    rgx=0
                    rgy=0

                    DO st_virt=1,tau
                        rgx=rgx+(Rvirt(st_virt,1)-r_mean_x)**2
                        rgy=rgy+(Rvirt(st_virt,2)-r_mean_y)**2
                    END DO

                    rgx=rgx/(tau)
                    rgy=rgy/(tau)
                    Rgyr(tr)=sqrt(rgx+rgy)





                    if (i.eq.vrr) then
                        rg(tr,1)=Rgyr(tr)
                    end if

                    if ((i.eq.vrr).and.(MODULO(tr,18).eq.0)) then
                        do st_virt=1,tau
                        vrecord(tr/18,st_virt)=MODULO(Rvirt(st_virt
     *                  ,1),L)
                        vrecord(tr/18,tau+st_virt)=MODULO(Rvirt
     *                  (st_virt,2),L)
                        end do
                    end if






                END DO


                Gaus_x=0
                Gaus_y=0

                DO tr=1,nvtraj
                    Rgyr_norm=log(Rgyr(tr)/(sum(Rgyr)/size(Rgyr)))
                    Gaus_x=Gaus_x+Gausv_x(tr)*Rgyr_norm
                    Gaus_y=Gaus_y+Gausv_y(tr)*Rgyr_norm
                END DO
                END IF !testing
                Gaus_x=Gaus_x/nvtraj
                Gaus_y=Gaus_y/nvtraj




                write(47,*)Gaus_x, Gaus_y,R(i,1),R(i,2)




                ! Langevin equation
                Flx=-gamma*V(i,1)+noi*Gaus_x
                Fly=-gamma*V(i,2)+noi*Gaus_y


        ! Repulsion (hard core)
       ! call repul_hcore(np,L,sigma,Ar,R2(i,1),R2(i,2),R2,Frepx,Frepy)
        ! Repulsion (soft)
        call repul_soft(i,np,L,sigma,Ar,R2(i,1),R2(i,2),R2,Frepx,Frepy)
        ! Repulsion (adjustable hard core)
        !call repul_ahc(np,L,sigma,Ar,b,R2(i,1),R2(i,2),R2,Frepx,Frepy)

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

c Analysis (Saving)


            !IF ((st.ge.start).and.(mod(st,incr).eq.0)) THEN
            !    DO i=1,np
            !      IF(save_temp)skor=skor+(V(i,1)**2+V(i,2)**2)
            !    END DO
            !END IF



            IF (save_coor.and.(st.ge.start).and.
     *    (mod(st,incr_coor).eq.0)) THEN
                DO i=1,np
                    write(50,*)R2(i,1),R2(i,2),V(i,1),V(i,2)
                END DO
            END IF



            ! This tr is different from tr before, =tr/18
            do tr=1,nvtraj/18
                do st_virt=1,tau*2
                    write(49,*)vrecord(tr,st_virt)
                end do
            end do

            do tr=1,nvtraj
                write(48,*)rg(tr,1),rg(tr,2)
            end do




        END DO



        !dur=real((st-start)/incr)

        ! Coordinates and velocities for last step
      !  do i=1,np
      !      write(49,*)R2(i,1),R2(i,2),V(i,1),V(i,2)
      !  end do
        !IF(save_temp)skor=skor/(np*dur)
        !IF(save_temp)write(60,*)skor


      !  CLOSE(49)
        IF(save_coor)CLOSE(50)
        !IF(save_temp)CLOSE(60)

        close(49)
        close(49)
        close(48)

        deallocate(R,Rb,R2,V,Rvirt)
        deallocate(Gausv_x,Gausv_y,Rgyr,gnumb_x,gnumb_y)
        deallocate(Rvirt1,Rvirt2)

        deallocate(vrecord,rg)


       END PROGRAM cognitive
!- - - - - - - - - Random seed initialisation - - - - - - - - -!
        subroutine init_random_seed()
          use Ziggurat
          INTEGER :: i, n, clock, zigseed
          INTEGER, DIMENSION(:), ALLOCATABLE :: ranseed

          CALL RANDOM_SEED(size = n)
          ALLOCATE(ranseed(n))

          CALL SYSTEM_CLOCK(COUNT=clock)

          ranseed=abs((clock+37*[ (i - 1, i = 1, n) ])**2)
          CALL RANDOM_SEED(PUT = ranseed)
          zigseed=sum(ranseed)
          call zigset(zigseed)

          DEALLOCATE(ranseed)
        end subroutine
!- - - - - - Generate random nonoverlapping circles - - - - - -!

        SUBROUTINE rand_nonoverlap_circle(circles,nCircles,radii,r,L)

        INTEGER i,j,nCircles,theormax,slownumb

        REAL radii,r,L,rand,x,y,dx,dy,tstart,tcurrent
        REAL distFromPrevCircles,totdistFromPrevCircles
        REAL circles(nCircles,2)

        real, dimension(2)    :: rand2

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
            DO WHILE (newCircleFound .eqv. .false.)
            totdistFromPrevCircles = 0.0
            ! Draw random coordinates of a candidate circle
            call random_number(rand2)
            x = radii + (L-radii*2)*rand2(1)
            y = radii + (L-radii*2)*rand2(2)

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


!- - - - - - - - - Repulsion (hard core) - - - - - - - - -!

        SUBROUTINE repul_hcore(i,np,L,sigma,Ar,Ri_x,Ri_y,R2,Frepx,Frepy)
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

        SUBROUTINE repul_soft(i,np,L,sigma,Ar,Ri_x,Ri_y,R2,Frepx,Frepy)
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
                    Rux=Rijx/Rij
                else
                    Rux=0.0
                end if
                if (Rijy.ne.0.0)then
                    Ruy=Rijy/Rij
                else
                    Ruy=0.0
                end if
                Frepx=Frepx+Ar*Lin*Rux
                Frepy=Frepy+Ar*Lin*Ruy
            end if
        END DO

        END SUBROUTINE

!- - - - - - - Repulsion (adjustable hard core) - - - - - - -!

        SUBROUTINE repul_ahc(i,np,L,sigma,Ar,b,Ri_x,Ri_y,R2,Frepx,Frepy)
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
                    Rux=Rijx/Rij
                else
                    Rux=0.0
                end if
                if (Rijy.ne.0.0)then
                    Ruy=Rijy/Rij
                else
                    Ruy=0.0
                end if
                Frepx=Frepx+Ar*Lin*Rux
                Frepy=Frepy+Ar*Lin*Ruy
            end if
        END DO

        END SUBROUTINE
