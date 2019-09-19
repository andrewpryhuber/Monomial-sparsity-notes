%%%%%%%%%%%%%%%%%
% set the number of variables manually depending on input system
% run commands in polysToFile.m2 to change input system 
numVars = 3;
syms x_1 x_2 x_3
allVars = [x_1, x_2 , x_3];
%%%%%%%%%%%%%%%%%

% read in generators from M2 
generators = polyArrToPolySys( readPolysFromFile('generators') , allVars); 

% read in reducibles from M2 (Note: these come from normal set to GRevLex GB)
reducibles = polyArrToPolySys( readPolysFromFile('reducibles') , allVars);


% maximal degree among reducibles, used to determine inital degree,
% currently set manually, can get max degree of reducibles in polysToFile.m2
maxDegReducibles = 7;
numReducibles = size(reducibles,1);


testrangemin = maxDegReducibles;
testrangemax = testrangemin + 3; %change to allow search in higher degree

data = zeros(testrangemax-testrangemin,3);

dataindex = 0;
for d = testrangemin:testrangemax
    
    dataindex = dataindex + 1;
    dtest = d;
   
    % represent reducibles as vectors in some degree dtest 
    reduciblesVec = polysys2vec(reducibles,dtest) ;
    
    % make Macaulay matrix of degree dtest
    Mdtest = getM(generators,dtest);

    % run 1-norm minimization in CVX to look for fewest # rows to represent reducibles
    % currently solver not seeming to work for larger problems 
    % declares feasible problems infeasible
    % maybe also need to tune objective function 
    cvx_begin quiet
    cvx_solver mosek
    variable h( size(Mdtest,1) , numReducibles); 
    reduciblesVec' ==  Mdtest' * h;
    obj = 0;
    for i = 1:numReducibles 
        obj = obj + norm( h(:,i) ,1);
    end
    minimize( obj );
    cvx_end
 
    data(dataindex,:) = [dtest , nnz( sum( abs(h) ,2)  > 1e-5), obj];

    
end

array2table(data, 'VariableNames',{'dtest','numNonzeroH','oneNormH'})

 



