%%%%%%%%%%%
% Requires the following additional files:
% - PNLA_MATLAB_OCTAVE functions by Kim Batseiler: https://github.com/kbatseli/PNLA_MATLAB_OCTAVE
% - CVX: http://cvxr.com/cvx/download/


clear 
% x1^3 - x1 + x2 + x3 - 10;
polysys{1,1} = [1 -1 1 1 -10];
polysys{1,2} = [3 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
% x2^3 + x1 - x2 + x3  - 10;
polysys{2,1} = [1 1 -1 1 -10];
polysys{2,2} = [0 3 0; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
% x3^3 + x1 + x2 - x3  - 10;
polysys{3,1} = [1 1 1 -1 -10];
polysys{3,2} = [0 0 3; 1 0 0; 0 1 0; 0 0 1; 0 0 0];

%p = S(1,2) + S(2,3) 
p{1,1} = [-1,1,-2,-1,2,1,-1,1,10,-10];
p{1,2} = [4,0,0; 3,1,0; 1,3,0; 3,0,1; 0,3,1; 1,0,3; 0,1,3; 0,0,4; 3,0,0; 0,0,3];


n=size(polysys{1,2},2); % number of variables is length of some multiindex
dsmall = 4; % choose degree in which to represent polynomial p as a vector


testrangemin = dsmall;
testrangemax = dsmall+4;

data = zeros(testrangemax-testrangemin,3);

dataindex = 0;
for d = testrangemin:testrangemax
    
    dataindex = dataindex + 1;
    dtest = d;
   
    % represent p as vector in some degree dtest >= dsmall
    pvectest = polysys2vec(p,dtest) ;
    
    % make Macaulay matrix of degree dtest
    Mdtest = getM(polysys,dtest);

    %run 1-norm minimization in CVX to look for fewest # rows to represent p 
    cvx_begin quiet
    cvx_solver mosek
    variable x(size(Mdtest,1)); 
    pvectest' == Mdtest'*x;
    minimize( norm(x,1) )
    cvx_end
    
    x = x - rem(x,1e-5)

    data(dataindex,:) = [dtest , nnz(abs(x) > 1e-5), norm(x,1)];

    
end

array2table(data, 'VariableNames',{'Degree','numNonzeroX','oneNormX'})

 
