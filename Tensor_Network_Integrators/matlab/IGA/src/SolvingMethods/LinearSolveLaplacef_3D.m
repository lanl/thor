function [ue, K11, Mt, f1, u1] = LinearSolveLaplacef_3D(ctpxe, ctpye, ctpze,...
                            knot1n, knot2n, knot3n, ...
                            pp1, pp2, pp3, ...
                            we, m3D, IndexE, ...
                            uindex1, ffunction)
[Ke, Me, fe] = GetKLaplace_3Df(ctpxe, ctpye, ctpze, knot1n, knot2n, knot3n, pp1, pp2, pp3, we, m3D, ffunction);

Kt = AssemblyK_sparse(Ke, IndexE, IndexE);
Mt = AssemblyK_sparse(Me, IndexE, IndexE);
ft = Assemblyf(fe, IndexE);

K11 = Kt(uindex1, uindex1);
f1 = ft(uindex1);


% Solve the system of equations
u1 = K11 \ f1; 

us = Getus1(u1, uindex1, IndexE);
ue = Getue(us, IndexE);
end