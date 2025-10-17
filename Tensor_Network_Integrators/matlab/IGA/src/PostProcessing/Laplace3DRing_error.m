function er = Laplace3DRing_error(ue, ctpxe, ctpye, ctpze, u_in, u_out, r_int, r_out)
[n1, n2] = size(ue);
er = zeros(n1, n2);
for i = 1:n1
    for j =1:n2
        r = sqrt(ctpxe(i,j)^2 + ctpye(i,j)^2);
        u_exact = (u_in*log(r_out/r)+ u_out*log(r/r_int))/log(r_out/r_int);
        er(i,j) = ue(i,j) - u_exact;
    end
end
end