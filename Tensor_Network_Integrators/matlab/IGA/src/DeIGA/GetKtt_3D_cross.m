function KT = GetKtt_3D_cross(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, tt_tol)
Kftt33 = GetKf_3D_tt_cross(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, tt_tol);
%[max_K, ~] = tt_max(Kftt33(:,:,:,:,:));
%fprintf('Max Jvm total: = %d\n', max_K);
%thresh = 1e-16;

%%
for i =1:3
  for j =1:3
    %[max_Kij, ~] = tt_max(Kftt33(:,:,:,i,j));
    %fprintf('Max Jvm_%d %d = %d \n', i, j, max_Kij);
    if 1 %max_K > (thresh * max_Kij)
      Kf_tt = round(Kftt33(:,:,:,i,j), tt_tol);
      R1 = Kf_tt.r(2);
      R2 = Kf_tt.r(3);
      fprintf('Jvm_%d%d: rank = %d %d\n', i, j, R1, R2);
      for r1 = 1:R1
        for r2 = 1:R2
          Kx = Kf_tt{1}(1,:,r1);
          Ky = Kf_tt{2}(r1,:,r2);
          Kz = Kf_tt{3}(r2,:,1);
          if i==1 && j==1
            KxT = CalDNDN_1D(knot1n, Kx);
            KyT = CalNN_1D(knot2n, Ky);
            KzT = CalNN_1D(knot3n, Kz);
          elseif i==1 && j==2
            KxT = CalNDN_1D(knot1n, Kx);
            KyT = CalDNN_1D(knot2n, Ky);
            KzT = CalNN_1D(knot3n, Kz);
          elseif i==1 && j==3
            KxT = CalNDN_1D(knot1n, Kx);
            KyT = CalNN_1D(knot2n, Ky);
            KzT = CalDNN_1D(knot3n, Kz);
          elseif i==2 && j==1
            KxT = CalDNN_1D(knot1n, Kx);
            KyT = CalNDN_1D(knot2n, Ky);
            KzT = CalNN_1D(knot3n, Kz);
          elseif i==2 && j==2
            KxT = CalNN_1D(knot1n, Kx);
            KyT = CalDNDN_1D(knot2n, Ky);
            KzT = CalNN_1D(knot3n, Kz);
          elseif i==2 && j==3
            KxT = CalNN_1D(knot1n, Kx);
            KyT = CalNDN_1D(knot2n, Ky);
            KzT = CalDNN_1D(knot3n, Kz);
          elseif i==3 && j==1
            KxT = CalDNN_1D(knot1n, Kx);
            KyT = CalNN_1D(knot2n, Ky);
            KzT = CalNDN_1D(knot3n, Kz);
          elseif i==3 && j==2
            KxT = CalNN_1D(knot1n, Kx);
            KyT = CalDNN_1D(knot2n, Ky);
            KzT = CalNDN_1D(knot3n, Kz);
          elseif i==3 && j==3
            KxT = CalNN_1D(knot1n, Kx);
            KyT = CalNN_1D(knot2n, Ky);
            KzT = CalDNDN_1D(knot3n, Kz);
          end
          if i*j*r1*r2 == 1
            KT = matrices_to_tt_matrix_fn({KxT, KyT, KzT});
          else
            KT = round(KT + matrices_to_tt_matrix_fn({KxT, KyT, KzT}), tt_tol);
          end
        end
      end
    end
  end
end

KT = round(KT, tt_tol);

end
