function polySystem = polyArrToPolySys(polys, allVars)

numEqs = size(polys{1} , 1 );

polySystem = cell(numEqs , 2);

for i = 1:numEqs
    currSymPoly = str2sym( polys{1}{i} );
    
    [c,m] = polyStrToCoeffs( currSymPoly , allVars ) ;  
    polySystem{i,1} = c ;
    polySystem{i,2} = m ;
end

end