program FD
    parameter (NX=81)
    common N,NM,FI(NX),AE(NX),AW(NX),AP(NX),Q(NX),X(NX)
    dimension FIEX(NX)
    character FILOUT*10
!  Open file
    print *,'Enter output file name:'
    read(*,1) FILOUT
    1 format(A10)
    open (UNIT=8,FILE=FILOUT)
! Read input data
    print *, 'Enter: DEN(SITY),VEL(OCITY),DIF(FUSION COEFF.)'
    read(*,*) DEN,VEL,DIF
    print *, 'Enter boundary values: FI0, FIN'
    read(*,*) FI0,FIN
    print *, 'Choose convection scheme: 1 - CDS, 2 - UDS'
    read(*,*) IC
!   print *, DEN,VEL,DIF,FI0,FIN

! Define the grid: EX - Expansion factor
! N - Number of nodes incl. boundary ones

    print *, 'Enter: XMIN, XMAX, EX, N'
    read(*,*) XMIN,XMAX,EX,N
    NM=N-1
    if(EX.EQ.1.) then
        DX=(XMAX-XMIN)/REAL(N-1)
    else
        DX=(XMAX-XMIN)*(1.-EX)/(1.-EX**(N-1))
    endif
    X(1)=XMIN
    do I=2,N
        X(I)=X(I-1)+DX
        DX=DX*EX
    end do
    
!  Initialize fields
    do I=1,N
        FI(I)=0.
    end do
    FI(1)=FI0
    FI(N)=FIN
    DENVEL=DEN*VEL
    ZERO=0.

    do I=2,NM
! Central difference convection approx. (CDS)
        if(IC.EQ.1) then
            AEC=DENVEL/(X(I+1)-X(I-1))
            AWC=-AEC
! Upwind convection approx. (UDS)
        elseif(IC.EQ.2) then
            AEC= MIN(DENVEL,ZERO)/(X(I+1)-X(I-1))
            AWC=-MAX(DENVEL,ZERO)/(X(I+1)-X(I-1))
        endif
! Central difference convection approx. (CDS)
        DXR=2./(X(I+1)-X(I-1))
        AED=-DIF*DXR/(X(I+1)-X(I-1))
        AWD=-DIF*DXR/(X(I)-X(I-1))
! Assemble coefficient matrix
        AE(I)=AEC+AED
        AW(I)=AWC+AWD
        AP(I)=-AW(I)-AE(I)
        Q(I)=0.
    end do
! Boundary conditions
    Q(2)=Q(2)-AW(2)*FI(1)
    AW(2)=0.
    Q(NM)=Q(NM)-AE(NM)*FI(N)
    AE(NM)=0.
! Solve equation system
    print *, 'Chosse solver: 1 - JACOBI, 2 - GS, 3 - GSOR, 4 - TDMA'
    read(*,) IS
    if(IS.EQ.1) then
        call JACOBI
    elseif(IS.EQ.2) then
        call GS
    elseif(IS.EQ.3) then
        call GSOR
    elseif(IS.EQ.4) then
        call TDMA
    endif


















end program FD
