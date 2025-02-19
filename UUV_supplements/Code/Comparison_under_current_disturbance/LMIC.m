function sys=LMIC(A,B1,B2,C1,D11,D12,gamma,N)
    setlmis([]);
    [X,n,Sx]=lmivar(1,[N,1]);
    [W,n,Sw]=lmivar(2,[6,N]);
   
    lmiterm([1,1,1,X],A,1,'s');  %AX+(AX)'
    lmiterm([1,1,1,W],B2,1,'s');  %B2*W+(B2*W)'
    lmiterm([1,1,2,0], B1);  %C1*X
    lmiterm([1,1,3,X],1,C1' ); %D12*W
    lmiterm([1,1,3,-W],1,D12' ); %W'*D12'
    lmiterm([1,2,2,0],-gamma^2);%-gamma
    lmiterm([1,2,3,0],D11'); %D12'
    lmiterm([1,3,3,0], -1 ); %-I
    %------------------------------------------
    lmiterm([-2,1,1,X],1,1);
    %------------------------------------------
    lmisys=getlmis ;
    
    [tmin, xfeas]=feasp(lmisys);
    if(tmin<0)
        disp('Feasible');
    else
        sys=[];
        return
        
    end
    X=dec2mat(lmisys, xfeas, X);
    W=dec2mat(lmisys, xfeas, W);
    
    K=W/X;
    
    sys=K;
    
end