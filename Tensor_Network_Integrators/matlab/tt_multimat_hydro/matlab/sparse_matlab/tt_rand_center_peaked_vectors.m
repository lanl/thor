function ytt = tt_rand_center_peaked_vectors(N_list, shape)
% GENERATE_CENTER_PEAKED_VECTORS Generate multiple vectors with center maximum
%
% Inputs:
%   N_list - Array of vector lengths (e.g. [5 7 9])
%   shape  - Optional, one of: 'triangular' (default), 'gaussian', 'random'
%
% Output:
%   vecs - Cell array of 1D vectors with center element as maximum

  if nargin < 2
    shape = 'gaussian';
  end

  vecs = cell(1, numel(N_list));

  for idx = 1:numel(N_list)
    N = N_list(idx);
    mid = ceil(N / 2);

    switch lower(shape)
      case 'triangular'
        vecs{idx} = [1:mid, mid-1:-1:1];
        vecs{idx} = vecs{idx}(1:N);  % handle even N by trimming

      case 'gaussian'
        x = linspace(-1, 1, N);
        vecs{idx} = exp(-10 * x.^2);  % narrow Gaussian peak at center

      case 'random'
        v = rand(1, N);
        v(mid) = max(v) + 1;
        vecs{idx} = v;

      otherwise
        error('Unknown shape type: %s', shape);
    end
  end
  ytt = cell2core(tt_tensor,vecs);

end %function