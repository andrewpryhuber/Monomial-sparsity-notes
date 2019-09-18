R = QQ[x_1..x_3];

h_0 = x_1^3 - x_1 + x_2 + x_3 - 10
h_1 = x_2^3 + x_1 - x_2 + x_3 - 10
h_2 = x_3^3 + x_1 + x_2 - x_3 - 10

I = ideal(h_0, h_1, h_2)

-- p = S(h_0, h_1) + S(h_1, h_2) so can be represented using 4 monomial multiples of the h_i
p = (x_2^3*h_0 - x_1^3 * h_1) + (x_3^3*h_1 - x_2^3* h_2);
p // (gens I) -- returns representation which requires 10 monomial multiples of the h_i

-- columns of M generate syzygy module  Syz(h_0, h_1, h_2)
M = syz (gens I);


-- computes largest degree monomial encountered when multiplying basis syzygies with generators h_0, h_1, h_2
numCol = numColumns M;
numRow = numRows M;
maxSyzDeg = 0;
for i from 0 to numRow - 1 do for j from 0 to numCol - 1 do maxSyzDeg = max(maxSyzDeg, first degree ( M_(i,j) * h_i ) );
maxSyzDeg 

