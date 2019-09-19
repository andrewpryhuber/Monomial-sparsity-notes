-- an example system to play with
-- more systems available in polyscollection.m2
R = QQ[x_1..x_3];
f_0 = x_1^3 - x_1 + x_2 + x_3 - 10
f_1 = x_2^3 + x_1 - x_2 + x_3 - 10
f_2 = x_3^3 + x_1 + x_2 - x_3 - 10
I = ideal(f_0, f_1, f_2)


-- construct reducibles r_i = m * b_i where m is action monomial
-- returns matrix with rows (r_i - normalform(r_i)) , pruned to remove zero rows
makeReducibles = (m, I) -> transpose ( compress(( m*lift(basis(R/I), R) )  - ( ( m*lift(basis(R/I), R) )) % I ))

-- write generators to file for import to matlab
gensMatrix = transpose gens I ;
file = "generators" << ""
for i from 0 to numRows(gensMatrix) - 1 do (file << toString( gensMatrix_(i,0) ) << endl )
file << close;

-- write r_i - normalform(r_i) for all reducibles i
m = x_2; -- action monomial
redsMatrix =  makeReducibles(m,I) ;
file = "reducibles" << ""
for i from 0 to numRows(redsMatrix) - 1 do (file << toString( redsMatrix_(i,0) ) << endl )
file << close;


-- -- returns max degree polynomial in a matrix
-- getMaxDeg = M -> (L = degrees M; L = flatten flatten L; return ( max apply( L, abs) + 1))
-- getMaxDeg(redsMatrix)

-- -- columns of M generate syzygy module  Syz(f_i)
-- M = syz (gens I);

-- -- computes largest degree monomial encountered when multiplying basis syzygies with generators
-- numCol = numColumns M;
-- numRow = numRows M;
-- maxSyzDeg = 0;
-- for i from 0 to numRow - 1 do for j from 0 to numCol - 1 do maxSyzDeg = max(maxSyzDeg, first degree ( M_(i,j) * f_i ) );
-- maxSyzDeg
