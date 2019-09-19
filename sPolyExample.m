%%%%%%%%%%%
% Requires additional files:
% - PNLA_MATLAB_OCTAVE functions by Kim Batseiler: https://github.com/kbatseli/PNLA_MATLAB_OCTAVE
% - CVX: http://cvxr.com/cvx/download/
%%%%%%%%%%%

clear 
format shortg

% x1^3 - x1 + x2 + x3 - 10;
polysys{1,1} = [1 -1 1 1 -10];
polysys{1,2} = [3 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
% x2^3 + x1 - x2 + x3  - 10;
polysys{2,1} = [1 1 -1 1 -10];
polysys{2,2} = [0 3 0; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
% x3^3 + x1 + x2 - x3  - 10;
polysys{3,1} = [1 1 1 -1 -10];
polysys{3,2} = [0 0 3; 1 0 0; 0 1 0; 0 0 1; 0 0 0];

n=size(polysys{1,2},2); % number of variables is length of some multiindex
dsmall = 5; % choose degree in which to represent polynomial p as a vector


testrangemin = dsmall;
testrangemax = dsmall+6;

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
    variable h(size(Mdtest,1)); 
    pvectest' == Mdtest'*h;
    minimize( norm(h,1) )
    cvx_end
    
%     h = h - rem(x,1e-5)

    data(dataindex,:) = [dtest , nnz(abs(h) > 1e-5), norm(h,1)];

    
end

array2table(data, 'VariableNames',{'P_d','numNonzeroH','oneNormH'})

 
