clear
%%%%%%%%%%%%%%%%%
% set the number of variables manually depending on input system
% run commands in polysToFile.m2 to change input system
numVars = 2;
allVars = sym('x_%d', [1 numVars]);
%%%%%%%%%%%%%%%%%

% read in generators from M2
generators = polyArrToPolySys( readPolysFromFile('generators') , allVars);

% read in reducibles from M2 (Note: these come from normal set to GRevLex GB)
reducibles = polyArrToPolySys( readPolysFromFile('reducibles') , allVars);


% maximal degree among reducibles, used to determine inital degree,
% currently set manually, can get max degree of reducibles in polysToFile.m2
maxDegReducibles = 16;
numReducibles = size(reducibles,1);


testrangemin = maxDegReducibles;
testrangemax = testrangemin + 5; %change to allow search in higher degree

data = zeros(testrangemax-testrangemin,3);

dataindex = 0;
for d = testrangemin:testrangemax
    
    dataindex = dataindex + 1;
    dtest = d;
    
    % represent reducibles as vectors in some degree dtest
    %(last argument = 1 makes representation sparse)
    reduciblesVec = polysys2vec(reducibles,dtest,1) ;

    % make Macaulay matrix of degree dtest
    % (last argument = 1 makes representation sparse)
    Mdtest = getM(generators,dtest,1);
    [numRows, numCols] = size(Mdtest);
    
    % run 1-norm minimization in CVX to look for fewest # rows to represent reducibles
    cvx_begin quiet
%     cvx_solver sdpt3
    cvx_solver mosek
%     cvx_solver gurobi
    variable h( numRows , numReducibles);
    reduciblesVec' ==  Mdtest' * h;
    obj = 0;
    for i = 1:numReducibles
        obj = obj + norm( h(:,i) ,1);
    end
%         obj = norm(sum(h,2),1);
    minimize( obj );
    cvx_end
    
    
    data(dataindex,:) = [dtest ,  nnz( sum( abs(h) ,2)  > 1e-6), obj];
    
end


array2table(data, 'VariableNames',{'dtest','numNonzeroHcvx','obj'})



