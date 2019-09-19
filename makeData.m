
numVars = 4;
syms x_1 x_2 x_3 x_4
allVars = [x_1, x_2 , x_3, x_4];

% read in generators from M2 
generators = polyArrToPolySys( readPolysFromFile('generators') , allVars); 

% read in reducibles from M2 (assume GRevLex for now to make these)
reducibles = polyArrToPolySys( readPolysFromFile('reducibles') , allVars);



% maximal degree among reducibles, used to determine inital degree,
% currently set manually
maxDegReducibles = 6;

numReducibles = size(reducibles,1);


testrangemin = maxDegReducibles;
testrangemax = testrangemin; %change to allow search in higher degree

data = zeros(testrangemax-testrangemin,3);

dataindex = 0;
for d = testrangemin:testrangemax
    
    dataindex = dataindex + 1;
    dtest = d;
   
    % represent p as vector in some degree dtest >= dsmall
    reduciblesVec = polysys2vec(reducibles,dtest) ;
    
    % make Macaulay matrix of degree dtest
    Mdtest = getM(generators,dtest,1);

    %run 1-norm minimization in CVX to look for fewest # rows to represent p 
    %currently solver not seeming to work - declares feasible problems
    %infeasible
    cvx_begin
    cvx_solver mosek
    variable h( size(Mdtest,1) , numReducibles); 
    reduciblesVec' ==  Mdtest' * h;
    minimize( norm( h,1) );
    cvx_end

    data(dataindex,:) = [dtest , nnz(abs( h ) > 1e-5), norm(h,1)];

    
end

array2table(data, 'VariableNames',{'P_d','numNonzeroH','oneNormH'})

 



