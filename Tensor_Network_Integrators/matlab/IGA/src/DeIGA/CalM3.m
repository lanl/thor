function Mt = CalM3(ctpxn, ctpyn, knot1n, knot2n)

[Gd1, Jd1, alphad1] = getIntegralData1D(knot1n);
[Gd2, Jd2, alphad2] = getIntegralData1D(knot2n);

MT = 0;
for ii = 1:length(Gd1)
    for jj = 1:length(Gd2)
        MT = MT + GetM3(knot1n, knot2n, Gd1(ii), Gd2(jj), Jd1(ii)*Jd2(jj), ctpxn, ctpyn, alphad1(ii), alphad2(jj));
    end
end

%% Return KT to 2D
[n1, n2, n3, n4] = size(MT);

Mt = reshape(permute(MT, [1, 2, 3, 4]), n1*n2, n3*n4);
% Mt = reshape(permute(MT, [2, 1, 4, 3]), n1*n2, n3*n4);
end