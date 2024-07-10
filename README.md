# Volterra_maps
Maple code for verifying tau function relations for integrable 4D Volterra map
Volterra map MAPLE CODE: 

This is the MAPLE code for verifying the computer algebra proof of Proposition 2.1 in the paper 
"A family of integrable maps associated with the Volterra lattice" by A.N.W. Hone, J.A.G. Roberts and P. Vanhaecke, 
which is freely available in preprint form here: https://arxiv.org/abs/2309.02336 (to appear in Nonlinearity). 

For the interested reader who does not have access to the MAPLE software package, the plain text of the MAPLE file is reproduced here: 

restart;
with(LinearAlgebra);
wrec := w[4]*w[3]*w[2] + w[2]*w[1]*w[0] + 2*w[2]^2*(w[3] + w[1]) + w[2]*(w[1]^2 + w[1]*w[3] + w[3]^2) + w[2]^3 + nu*w[2]*(w[3] + w[2] + w[1]) + b*w[2] + a;
;
;
tauseven := wrec;
for n from 0 to 4 do
    tauseven := simplify(subs(w[n] = tau[n]*tau[n + 3]/(tau[n + 1]*tau[n + 2]), tauseven));
end do;
tauseven := numer(tauseven);
NULL;
NULL;
NULL;
NULL;
NULL;
NULL;
NULL;
topsol := solve(tauseven, tau[7]);
for n from 7 to 10 do
    t[n] := topsol;
    for m from 0 to 6 do
        t[n] := factor(simplify(subs(tau[m] = t[m + n - 7], t[n])));
    end do;
    print(t[n]);
end do;
NULL;
NULL;
NULL;
botsol := solve(tauseven, tau[0]);
for n to 3 do
    t[-n] := botsol;
    for m from 0 to 6 do
        t[-n] := factor(simplify(subs(tau[m + 1] = t[m - n + 1], t[-n])));
    end do;
    print(t[-n]);
end do;
;
NULL;
f := (i, j) -> t[i + j - 5]*t[i - j + 6];
M := Matrix(5, f);
simplify(Determinant(M));
nullM := NullSpace(M);
nops(nullM);
kernelvec := op(1, nullM);
nullity := 5 - Rank(M);
for j to 5 do
    al[j] := kernelvec[j];
end do;
;
NULL;
NULL;
NULL;
for j to 4 do
    alw[j] := al[j];
    for m from 0 to 3 do
        alw[j] := subs(t[m] = w[m]*t[m + 1]*t[m + 2]/t[m + 3], alw[j]);
    end do;
end do;
for j to 4 do alw[j] := simplify(alw[j]); print(alw[j]); end do;
k[1] := denom(alw[4]);
NULL;
alw[5] := 1;
for j to 5 do
    alpha[j] := simplify(-k[1]*alw[j]*denom(alw[1]));
end do;
NULL;
k[2] := collect(collect(collect(simplify((k[1]^2 + alpha[2])/a), a), b), nu);
NULL;
simplify(alpha[1] - k[1]);
simplify(-a*k[2] + k[1]^2 + alpha[2]);
simplify(alpha[3] - a*(a*k[2] - 2*k[1]^2));
simplify(alpha[4] - a*(a^2*k[1] + b*k[1]^2 + nu*k[1]*k[2] + k[2]^2));
simplify(alpha[5] + k[1]*(a^2*k[1] + b*k[1]^2 + nu*k[1]*k[2] + k[2]^2));
;
