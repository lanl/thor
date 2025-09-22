function [ue, Kt, Mt, u1] = LinearSolveLaplace_3D(ctpxe, ctpye, ctpze,...
                            knot1n, knot2n, knot3n, ...
                            pp1, pp2, pp3, ...
                            we, m3D, IndexE, ...
                            uindex1, uindex2, ubc)

[Ke, Me] = GetKLaplace_3D(ctpxe, ctpye, ctpze, knot1n, knot2n, knot3n, pp1, pp2, pp3, we, m3D);

Kt = AssemblyK_sparse(Ke, IndexE, IndexE);
Mt = AssemblyK_sparse(Me, IndexE, IndexE);

K11 = Kt(uindex1, uindex1);
K12 = Kt(uindex1, uindex2);

u2 = ubc;
f1 = -K12 * u2';

% Solve the system of equations
u1 = K11 \ f1; 

us = Getus(u1, uindex1, u2, uindex2);
ue = Getue(us, IndexE);
end