function c=wenokCoeff(rq,k)
  %
  c = zeros(k,k);
  for r=0:k-1
      for j=0:k-1
          c(r+1,j+1) = lagrangeInterp(rq,r,j,k);
      end
  end
end

function c=lagrangeInterp(rq,r,j,k)
  %
  c=0;
  %
  for m=j+1:k
      num   = 0;
      denom = 1;
      for l=0:k
          if(l~=m)
              prod = 1;
              for q=0:k
                  if(q~=m && q~=l)
                      prod = prod * (rq + r - q);
                  end
              end
              num   = num + prod;
              denom = denom * (m-l);  
          end
      end
      %
      c = c + num/denom;
      %
  end
end