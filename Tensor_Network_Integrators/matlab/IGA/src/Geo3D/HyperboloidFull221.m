function [wn, ctpxn, ctpyn, ctpzn, knot1n, knot2n, knot3n] = HyperboloidFull221(rr, RR, r0, HH, nn1, nn2, nn3)
%% Define inital Geometry
p1 = 2; p2 = 2; p3 = 1;
n1 = 3; n2 = 9; n3 = 2;
%
pi = 3.141654;
ff = pi/2;
knot1 = NonUniformKnot(n1, p1);
knot2  = [0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1];
knot3 = NonUniformKnot(n3, p3);

ctpx2D = [-RR         -RR         0.0   RR       RR       RR       0.0  -RR       -RR;...
          -rr^2/RR    -rr^2/RR    0.0   rr^2/RR  rr^2/RR  rr^2/RR  0.0  -rr^2/RR  -rr^2/RR;... 
          -RR         -RR         0.0   RR       RR       RR       0.0  -RR       -RR];

ctpy2D = [0.0    RR       RR       RR       0.0   -RR       -RR       -RR       0.0;...
          0.0    rr^2/RR  rr^2/RR  rr^2/RR  0.0   -rr^2/RR  -rr^2/RR  -rr^2/RR  0.0;...
          0.0    RR       RR       RR       0.0   -RR       -RR       -RR       0.0];

ctpz2D = HH*[-0.5  -0.5  -0.5  -0.5  -0.5  -0.5  -0.5  -0.5  -0.5;...
              0.0   0.0   0.0   0.0   0.0  0.0   0.0   0.0   0.0;...
              0.5   0.5   0.5   0.5   0.5  0.5   0.5   0.5   0.5];
%%
ctpx = zeros(n1, n2, n3);

ctpx(:,:,1) = (r0/rr)*ctpx2D;
ctpx(:,:,2) = ctpx2D;
%%
ctpy = zeros(n1, n2, n3);

ctpy(:,:,1) = (r0/rr)*ctpy2D;
ctpy(:,:,2) = ctpy2D;
%%

ctpz = zeros(n1, n2, n3);

ctpz(:,:,1) = ctpz2D;
ctpz(:,:,2) = ctpz2D;
%PlotCtp3D(ctpx, ctpy, ctpz)


% control weight
wi = zeros(n1, n2, n3);

w2D = [1.0       cos(ff/2)         1.0      cos(ff/2)         1.0    cos(ff/2)        1.0     cos(ff/2)        1.0;...
       RR/rr     RR/rr*cos(ff/2)   RR/rr    RR/rr*cos(ff/2)   RR/rr  RR/rr*cos(ff/2)  RR/rr   RR/rr*cos(ff/2)  RR/rr;...
       1.0       cos(ff/2)         1.0      cos(ff/2)         1.0    cos(ff/2)        1.0     cos(ff/2)        1.0];

wi(:,:,1) = w2D;
wi(:,:,2) = w2D;

%% Refinement
[wn, ctpxn, ctpyn, ctpzn, knot1n]  = KnotInsertMW1_3D(wi, ctpx, ctpy, ctpz, knot1, sref(nn1, ReduceKnotUnbound(knot1)));
[wn, ctpxn, ctpyn, ctpzn, knot2n]  = KnotInsertMW2_3D(wn, ctpxn, ctpyn, ctpzn, knot2, sref([nn2 nn2 nn2 nn2], ReduceKnotUnbound(knot2)));
[wn, ctpxn, ctpyn, ctpzn, knot3n]  = KnotInsertMW3_3D(wn, ctpxn, ctpyn, ctpzn, knot3, sref(nn3, ReduceKnotUnbound(knot3)));

end