KK = ZZ / 30097;
R = KK[x_1,x_2,MonomialOrder=>GRevLex];
red = matrix({{x_1^5*x_2^3,x_1^3*x_2^4,x_1*x_2^5,x_1^5*x_2^2,x_1^3*x_2,x_1,x_1^4*x_2^4,x_1^3*x_2^4,x_1^2*x_2^5,x_1*x_2^5,x_2^6}});
b = matrix({{x_1^4*x_2^3_R,x_1^3*x_2^3_R,x_1^2*x_2^4_R,x_1^3*x_2^2_R,x_1^2*x_2^3_R,x_1*x_2^4_R,x_1^2*x_2^2_R,x_1*x_2^3_R,x_1*x_2^2_R,x_1*x_2_R,x_2^5_R,x_2^4_R,x_2^3_R,x_2^2_R,x_2_R,x_1^4*x_2^2_R,x_1^2*x_2_R,1_R}});
eq1 = 72*x_1^6*x_2^3 + 91*x_1^5*x_2^3 + 90*x_1^4*x_2^3 + 34*x_1^4*x_2^2 + 70*x_1^3*x_2^3 + 20*x_1^3*x_2^2 + 4*x_1^2*x_2^3 + 75*x_1^2*x_2^2 + 48*x_1*x_2^3 + 51*x_1^2*x_2 + 91*x_1*x_2^2 + 62*x_2^3 + 61*x_1*x_2 + 86*x_2^2 + 81*x_2 + 58;
eq2 = 19*x_1^6*x_2^3 + 24*x_1^5*x_2^3 + 89*x_1^4*x_2^3 + 3*x_1^4*x_2^2 + 49*x_1^3*x_2^3 + 17*x_1^3*x_2^2 + 98*x_1^2*x_2^3 + 72*x_1^2*x_2^2 + 48*x_1*x_2^3 + 51*x_1^2*x_2 + 6*x_1*x_2^2 + 5*x_2^3 + 69*x_1*x_2 + 8*x_2^2 + 53*x_2 + 10;
eqs = matrix({{eq1,eq2}});
I = ideal eqs;
gbTrace = 0;
Q = R/I;
b0 = lift(basis Q,R);
use R
S = (coefficients(b%I,Monomials => b0))_1;
if numcols b0 <= numcols b then (
  Sinv = transpose(S)*inverse(S*transpose(S));
) else (
  Sinv = inverse(transpose(S)*S)*transpose(S);
)
AM = Sinv*((coefficients(red%I,Monomials => b0))_1);
pp = red - b*AM;
A = pp // eqs;
gbRemove(I);
M = kernel eqs;
A = A % M;
-- "cache/matrix.test.1.txt" << toString A << close;
end
load "LAG-template-generation-example.m2"



redsMatrix = transpose (pp);
redsMatrix
