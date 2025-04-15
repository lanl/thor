function Aqtt = convert_tt_matrix_to_qtt_matrix_fn(Att,tol,dims)

r = Att.r;
n = Att.n;
m = Att.m;
G = core2cell(Att);
numcore = numel(G);
Gnew = [];
for i = 1:numcore
  if ismember(i,dims) && (n(i)~=2) && (m(i)~=2)
    crc = reshape(G{i},r(i)*n(i),m(i)*r(i+1));
    %get the dimension for reshapes
    l = log2(n(i)); %number of qtt cores
    ln = 2*ones(1,l);ln(1) = ln(1)*r(i);
    lm = 2*ones(1,l);lm(end) = lm(end)*r(i+1);
    cqtt = tt_matrix(crc,tol,ln,lm);
    gqtt = core2cell(cqtt);
    % reshape the first and 
    crsize = [r(i),2,2,size(gqtt{1},4)];
    gqtt{1} = reshape(gqtt{1},crsize);
    %
    crsize = [size(gqtt{end},1),2,2,r(i+1)];
    gqtt{end} = reshape(gqtt{end},crsize);
    %
    Gnew = [Gnew; gqtt];
    
  else
    Gnew = [Gnew; G(i)];
  end
end
Aqtt = round(cell2core(tt_matrix,Gnew),tol);
end
