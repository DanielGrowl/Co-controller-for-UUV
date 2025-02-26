function [Pack] = model_pack(A,B2,B1,C1,D12,D11)
   Pack{1,1} = double(A);
   Pack{2,1} = double(B2);
   Pack{3,1} = double(B1);
   Pack{4,1} = double(C1);
   Pack{5,1} = double(D12);
   Pack{6,1} = double(D11);
end