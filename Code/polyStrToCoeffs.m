% convert a polynomial from representation as a sum of symbolic terms into
% a matrix representing coefficients and a matrix representing monomials 
%
% coeffMat is a single row vector
% monoMat has one row for every term. Rows of monoMat are multi-indices 


function [coeffMat , monoMat] = polyStrToCoeffs(poly, allVars)

[c, m] = coeffs(poly);

numTerms = size(c,2);
numVars = size(allVars,2); 

coeffMat = c;
monoMat = zeros( numTerms, numVars );

for i= 1: numTerms
    monoMat(i, :) = multiindex(m(1,i) , allVars);
end


end
