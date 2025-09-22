function RHS=calculateRHS(t,x,eps_tt,u,A,F,B,isA,isF,n)
  %
  RHS = tt_zeros(n);
  %
  if(~isempty(B))
    error("u(t) is not yet implemented")
  end 
  %
  if(isA)
    %
    RHS = round(RHS + tt_contract(A,x,eps_tt), eps_tt); 
    %
  end
  %
  if(isF)
    %
    x2 = round(tkron(x,x),eps_tt);
    %
    RHS = round(RHS + tt_contract(F,x2,eps_tt), eps_tt); 
    %
  end
  %
end
%
