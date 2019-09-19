% compute the multiindex of a monomial mon with respect to the
% set of indeterminates allVars, returns row vector

function mdeg = multiindex( mon, allVars)

numVars = size(allVars,1);
mdeg = zeros(1,numVars);

varIndex = 0;
for indet = allVars
    varIndex = varIndex +1;
    mdeg(1,varIndex ) = polynomialDegree(mon,indet);
end
    