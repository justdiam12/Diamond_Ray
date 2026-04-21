

function [TL] = TL_geometrical_depth(r, r0, D)
w1=sqrt(1+(r./r0).^2);
w2=r.^2;

W=w2./w1;

TL=10 *log10(W)+40*log10(D/D(1));



