program test
    real, dimension (100) :: X
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
        print *, X(I)
    end do
end program test
