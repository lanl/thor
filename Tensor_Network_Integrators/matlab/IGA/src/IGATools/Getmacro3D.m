function m3D = Getmacro3D(knot1, knot2, knot3)
k1 = ReduceKnotFull(knot1);
k2 = ReduceKnotFull(knot2);
k3 = ReduceKnotFull(knot3);
m1 = length(k1)-1;
m2 = length(k2)-1;
m3 = length(k3)-1;

m3D = zeros(m1*m2*m3, 6);
%% The order of extraction is 1st-2nd-3rd direction
for i1 = 1:m1
    for i2 = 1:m2
        for i3 = 1:m3
            idx = m1*(i2-1) + i1 + m1*m2*(i3-1);
            m3D(idx,:) = [k1(i1), k1(i1 + 1),...
                          k2(i2), k2(i2 + 1),...
                          k3(i3), k3(i3 + 1)];
        end
    end
end

end