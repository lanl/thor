function Xtt = tt_batch_svd(X, num_batches, eps_tt)
  % Get number of snapshots
  Nt = size(X,2);
  % Split data into batches
  edges = round(linspace(1,Nt+1,num_batches+1));
  %
  X_batches = cell(1,num_batches);
  %
  for i=1:num_batches
    idx_start    = edges(i);
    idx_end      = edges(i+1)-1;
    X_batches{i} = X(:,idx_start:idx_end);
  end
  
  % Decompose each batch using TT-SVD and assemble
  for i=1:num_batches
    % TT-SVD on batch
    temp = tt_tensor(X_batches{i}, eps_tt);
    
    % Extract TT cores
    r = temp.r;
    G = core2cell(temp);

    % Expand middle core to full time length with zeros
    Gtemp                            = zeros(r(2), Nt, r(1));
    Gtemp(:,edges(i):edges(i+1)-1,:) = G{2};
    %
    G{2} = Gtemp;
    % Reassemble and add to final TT object
    if(i==1)
      Xtt = cell2core(tt_tensor,G);
    else
      Xtt = round(Xtt + cell2core(tt_tensor,G), eps_tt);
    end
  end
end