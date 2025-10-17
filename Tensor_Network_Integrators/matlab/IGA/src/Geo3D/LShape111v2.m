function [wn, ctpxn, ctpyn, ctpzn, knot1n, knot2n, knot3n] = LShape111v2(LL, nn1, nn2, nn3)
%
n1 = 3; n2 = 2; n3 = 2;
p1 = 1.0; p2 = 1.0; p3 = 1.0;
knot1 = NonUniformKnot(n1, p1);
knot2 = NonUniformKnot(n2, p2);
knot3 = NonUniformKnot(n3, p3);
%%
ctpx2D = LL*[-1.0  0.0;...
             -1.0  0.0;... 
              1.0  1.0];

ctpy2D = LL*[1.0   1.0;...
            -1.0   0.0;...
            -1.0   0.0];
%%
ctpx = zeros(n1, n2, n3);

ctpx(:,:,1) = ctpx2D;
ctpx(:,:,2) = ctpx2D;
%%
ctpy = zeros(n1, n2, n3);

ctpy(:,:,1) = ctpy2D;
ctpy(:,:,2) = ctpy2D;
%%

ctpz = zeros(n1, n2, n3);

ctpz(:,:,1) = zeros(n1, n2);
ctpz(:,:,2) = ones(n1, n2);

%PlotCtp3D(ctpx, ctpy, ctpz)

%weight
wi = ones(n1, n2, n3);

%% Refinement
[wn, ctpxn, ctpyn, ctpzn, knot1n]  = KnotInsertMW1_3D(wi, ctpx, ctpy, ctpz, knot1, sref([nn1 nn1], ReduceKnotUnbound(knot1)));
[wn, ctpxn, ctpyn, ctpzn, knot2n]  = KnotInsertMW2_3D(wn, ctpxn, ctpyn, ctpzn, knot2, sref(nn2, ReduceKnotUnbound(knot2)));
[wn, ctpxn, ctpyn, ctpzn, knot3n]  = KnotInsertMW3_3D(wn, ctpxn, ctpyn, ctpzn, knot3, sref(nn3, ReduceKnotUnbound(knot3)));

end