function y = funmap(n, mmap)
  %
  M    = size(n,1);
  dims = [mmap.rows, mmap.Nt];
  %
  j    = ones(M,1);
  mult = 1;
  for k = 1:size(n,2)
    j    = j + (n(:,k)-1)*mult;
    mult = mult * dims(k);
  end
  %
  rid = mod(j-1, mmap.rows) + 1;
  cid = floor((j-1)/mmap.rows) + 1;
  %
  edges    = [mmap.startCol(:); mmap.Nt+1];
  chunkID  = discretize(cid, edges);
  prevEnds = [0; cumsum(mmap.colCount(1:end-1))];
  localC   = cid - prevEnds(chunkID);
  %
  y = zeros(M,1);
  for c = 1:mmap.nChunks
    idx = (chunkID == c);
    if ~any(idx), continue, end
    % 
    S  = size(mmap.maps{c}.Data.Q);
    %
    lin = sub2ind(S, rid(idx), localC(idx));
    %
    y(idx) =  mmap.maps{c}.Data.Q(lin);
    %
  end
end