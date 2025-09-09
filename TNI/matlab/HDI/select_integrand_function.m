function [ fun, Iexa ] = select_integrand_function( itest, d )

if itest==1
  % oscillatory
  fprintf("======================= genz function 1 (oscillatory) ===============================\n");
  fun  = @g1;
  Iexa = 0.5;
  %fun  = @g1_old;
  %Iexa = 1.;
elseif itest == 2
  % product peak
  fprintf("======================= genz function 2 (product peak-Lorentzian) ===================\n");
  fun  = @g2;
  Iexa = 1;
elseif itest == 3
  % corner peak
  fprintf("======================= genz function 3 (corner peak) ===============================\n");
  fun  = @g3;
  Iexa = 1/factorial(d+1);
elseif itest == 4
  % gaussian
  fprintf("======================= genz function 4 (Gaussian) ==================================\n");
  fun  = @g4;
  Iexa = ( sqrt(pi)/2 * erf(1) )^d;
elseif itest == 5
  % simple exponential (continuous)
  fprintf("======================= genz function 5 (continuous-exponential) ===================\n");
  fun  = @g5;
  Iexa = (1-1/exp(1))^d;
elseif itest == 6
  % ANOVA example
  fprintf("======================= ANOVA function (test for MC method) ========================\n");
  fun  = @anova;
  Iexa = 1;

elseif itest == 7
  % ANOVA example
  fprintf("======================= ANOVA 3/4 function (test for MC method) ========================\n");
  fun  = @anova34;
  Iexa = 1.0;
elseif itest == 8
  fprintf("======================= Chebishev polynomial degree 10 ========================\n");
  fun = @(x) (x(1)-1/2)^10*sign(x(1)-1/2) + ...
  (1/d)*sum(512.*x.^10 - 1280.*x.^8 + 1120.*x.^6 - 400.*x.^4 + 50.*x.^2 - 1);
  Iexa = -1/99;
end

end