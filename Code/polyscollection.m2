R = QQ[x_1,x_2];
f_0 = 2*x_1^2 + x_2^2 +3*x_2 -12;
f_1 = x_1^2 - x_2^2 +x_1 + 3*x_2 -4;
I = ideal(f_0, f_1);
 

R = QQ[x_1..x_3];
f_0 = x_1^3 - x_1 + x_2 + x_3 - 10
f_1 = x_2^3 + x_1 - x_2 + x_3 - 10
f_2 = x_3^3 + x_1 + x_2 - x_3 - 10
I = ideal(f_0, f_1, f_2)


R = QQ[x_1, x_2, x_3, x_4];
f_0 = x_2^2*x_3 + 2 * x_1*x_2*x_4 -2*x_1 -x_3
f_1 = -x_1^3*x_3+ 4*x_1*x_2^2*x_3 + 4*x_1^2*x_2*x_4 + 2*x_2^3*x_4+ 4*x_1^2 -10*x_2^2
f_2 = 2*x_2*x_3*x_4 + x_1*x_4^2 - x_1 -2*x_3
f_3 = -x_1*x_3^3 + 4*x_2*x_3^2*x_4 + 4 *x_1*x_3*x_4^2 + 2*x_2*x_4^3 + 4*x_1*x_3 + 4*x_3^2 -10*x_2*x_4 -10*x_4^2 +2
I = ideal(f_0,f_1,f_2,f_3)


R = QQ[x_1,x_2];
f_0 = 72*x_1^6*x_2^3 + 91*x_1^5*x_2^3 + 90*x_1^4*x_2^3 + 34*x_1^4*x_2^2 + 70*x_1^3*x_2^3 + 20*x_1^3*x_2^2 + 4*x_1^2*x_2^3 + 75*x_1^2*x_2^2 + 48*x_1*x_2^3 + 51*x_1^2*x_2 + 91*x_1*x_2^2 + 62*x_2^3 + 61*x_1*x_2 + 86*x_2^2 + 81*x_2 + 58;
f_1 = 19*x_1^6*x_2^3 + 24*x_1^5*x_2^3 + 89*x_1^4*x_2^3 + 3*x_1^4*x_2^2 + 49*x_1^3*x_2^3 + 17*x_1^3*x_2^2 + 98*x_1^2*x_2^3 + 72*x_1^2*x_2^2 + 48*x_1*x_2^3 + 51*x_1^2*x_2 + 6*x_1*x_2^2 + 5*x_2^3 + 69*x_1*x_2 + 8*x_2^2 + 53*x_2 + 10;
I = ideal (f_0,f_1);




-- 5pt rel pose - Demazure
R = QQ[x_1,x_2,x_3];
A0 = random(R^3,R^3), A1 = random(R^3,R^3), A2 = random(R^3,R^3), A3 = random(R^3,R^3);
E = A0 + x_1*A1 + x_2*A2 + x_3*A3;
I = ideal(det(E)) + ideal (2 * (E*transpose(E))*E - trace(E*transpose(E))*E);


-- 6pt relpose focal
R = QQ[x_1,x_2,x_3];
A0 = random(R^3,R^3), A1 = random(R^3,R^3), A2 = random(R^3,R^3);
F = A0 + x_1*A1 + x_2*A2;
Q = matrix{{1,0,0},{0,1,0},{0,0,x_3}};
I = ideal( det(F)) + ideal(2*(F*Q*transpose(F)*Q)*F - trace(F*Q*transpose(F)*Q)*F)


vec = mat -> ( vectorization = submatrix(mat, {0},); for i from 1 to numRows(mat)-1 do
        (vectorization = vectorization |   submatrix(mat, {i},) ) ;
	    return transpose vectorization)


-- abs pose quiver
R = QQ[x_1,x_2,x_3,x_4];
A = random(R^4, R^9);
Q = matrix{{1 +x_1^2 - x_2^2 -x_3^2 , 2*(x_1*x_2 - x_3), 2*(x_1*x_3+x_2)},
    {2*(x_1*x_2 + x_3) , 1-x_1^2+x_2^2-x_3^2, 2*(x_2*x_3-x_1)},
    {2*(x_1*x_3 - x_2), 2*(x_2*x_3+x_1) , 1-x_1^2-x_2^2 +x_3^2}};
K = matrix{{1,0,0},{0,1,0},{0,0,x_4}};
I = ideal(A*vec(Q*K))


