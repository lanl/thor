function [x,xtt,t,dt,K]=prepareTrainingData(foldername,tint_order,th,skip,verb)
  %
  if(~exist("verb","var"))
    verb = true;
  end
  %
  if(verb)
    fprintf("Starting reading snapshots\n")
  end
  %
  xtt = tt_read(foldername+"/snapshots.tt");
  %
  fid = fopen(foldername+"/snapshots.times","r");
  %
  t = fread(fid,"double"); 
  %
  fclose(fid);
  %
  if(skip>1)
    %
    G = core2cell(xtt);
    %
    G{5} = G{5}(:,1:skip:end,:);
    %
    xtt = cell2core(tt_tensor,G);
    %
    xtt = round(xtt,1e-14);
    %
    t = t(1:skip:end);
    %
  end
  %
  x = full(xtt,xtt.n');
  %
  % calculate time step for training data
  %
  dt = (t(end)-t(1))/(length(t)-1);
  %
  % calculate the timestep where training ends
  %
  [~,K]=min(abs(t-th));
  % 
  % Correct time step values according to the time integration scheme's order
  %
  K  =  K - tint_order;
  %
  if(verb)
    fprintf("Finished reading snapshots\n")
  end
  %
end