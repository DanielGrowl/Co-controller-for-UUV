function [X,W,gamma]=LMID(A,B1,B2,C1,D11,D12,N)
    setlmis([]);
    X=lmivar(1,[N 1]);
    W=lmivar(2,[6 N]);
    g2=lmivar(1,[1 1]); 
    
    lmiterm([1,1,1,X],A,1,'s');  %AX+(AX)'
    lmiterm([1,1,1,W],B2,1,'s'); %B2*W+(B2*W)'
    lmiterm([1,1,2,0],B1); %B1
    lmiterm([1,1,3,-W],1,D12');%W*D2'
    lmiterm([1,1,3,X],1,C1');%X*C1'
    lmiterm([1,2,2,0],-1); %-I
    lmiterm([1,2,3,0],D11' );%D11'
    lmiterm([1,3,3,g2], -1,1);%-gamma^2*I
    
    lmiterm([-2,1,1,X],1,1);
    lmisys=getlmis ;
    n=decnbr(lmisys);
    c=zeros(n,1);
    for j=1:n
        bj=defcx(lmisys, j, g2);
        c(j)=bj;
    end
    [copt,xopt]=mincx(lmisys,c );
    X=dec2mat(lmisys, xopt, X);
    W=dec2mat(lmisys, xopt, W);
    gamma=sqrt(copt); 
end