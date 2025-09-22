function [L2, u_error] = CalH1L2E_Lshape( ...
    ue, ctpxe, ctpye, ctpze, ...
    k1, k2, k3, pp1, pp2, pp3, we, m2D)

  % Gauss points & weights
  [Gpx, wx] = GaussPoint(pp1+1);
  [Gpy, wy] = GaussPoint(pp2+1);
  [Gpz, wz] = GaussPoint(pp3+1);

  nGPx = numel(Gpx);
  nGPy = numel(Gpy);
  nGPz = numel(Gpz);
  nGP  = nGPx * nGPy * nGPz;      % total Gauss points per element

  [nE, ~] = size(ctpxe);          % number of elements
  u_error  = zeros(nE, nGP);      % preallocate error storage

  % Precompute weight products and parametric shifts
  Wq  = zeros(1, nGP);
  Pq  = zeros(nGP, 3);            % store [ξ, η, ζ] for each q
  q   = 0;
  for ix = 1:nGPx
    for iy = 1:nGPy
      for iz = 1:nGPz
        q         = q + 1;
        Wq(q)     = wx(ix) * wy(iy) * wz(iz);
        Pq(q, :)  = [Gpx(ix), Gpy(iy), Gpz(iz)];
      end
    end
  end

  L2n = 0;
  L2d = 0;

  % Loop over elements
  for e = 1:nE
    % extract element‐specific parametric ranges
    x1 = m2D(e,1); x2 = m2D(e,2);
    y1 = m2D(e,3); y2 = m2D(e,4);
    z1 = m2D(e,5); z2 = m2D(e,6);

    % midpoint & half‐range
    xm = (x1 + x2)/2; dx = (x2 - x1)/2;
    ym = (y1 + y2)/2; dy = (y2 - y1)/2;
    zm = (z1 + z2)/2; dz = (z2 - z1)/2;

    % constant Jacobian factor from parametric → reference
    J2 = dx * dy * dz;

    % now do all Gauss points in one flat loop
    for q = 1:nGP
      % map to physical ξ,η,ζ
      xi = xm + dx * Pq(q,1);
      yi = ym + dy * Pq(q,2);
      zi = zm + dz * Pq(q,3);

      % evaluate basis, jacobian and solution at this Gauss pt
      [NN, ~, J1] = GetB1D_3D( k1, k2, k3, ...
                      xi, yi, zi, ...
                      we(e,:), pp1, pp2, pp3, ...
                      ctpxe(e,:), ctpye(e,:), ctpze(e,:) );

      % physical coords & approx.‐value
      xe = ctpxe(e,:) * NN;
      ye = ctpye(e,:) * NN;
      ug = ue(e,:)    * NN;

      % analytical & error
      ua = (1/(2*pi^2)) * sin(pi*xe) * sin(pi*ye);
      du = abs(ug - ua);

      u_error(e, q) = du;

      % accumulate weighted norms
      wJ = Wq(q) * J1 * J2;
      L2n = L2n + wJ * du;
      L2d = L2d + wJ * sqrt(abs(ua));
    end
  end

  % final ratio
  L2 = L2n / L2d;
end
