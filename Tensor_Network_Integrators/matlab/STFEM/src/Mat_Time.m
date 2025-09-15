function [MT] = Mat_Time(T_Loc)

%Compute Matrices for time interval
% T_Loc  contains coordinates for starting and ending
   
    t1=T_Loc(1); t2=T_Loc(2); % Starting and Ending Points
    ht=t2-t1;                 % Length of the time interval
    t_md=mean([t1,t2]);
    
    function fy= Phit1(t)
            fy= (t2-t)/ht;
    end 

    function fy= Phit2(t)
            fy= (t-t1)/ht;
    end 


  DPhi1=(-1/ht); DPhi2=(1/ht);
  
  MT=zeros(2,2);
  % Integration Matrices
  
  MT(1,1)=(ht/6)*(DPhi1*Phit1(t1) + ...
     4*DPhi1*Phit1(t_md)+...
     DPhi1*Phit1(t2));

 MT(1,2)=(ht/6)*(DPhi2*Phit1(t1) + ...
     4*DPhi2*Phit1(t_md)+...
     DPhi2*Phit1(t2));


 MT(2,1)=(ht/6)*(DPhi1*Phit2(t1) + ...
     4*DPhi1*Phit2(t_md)+...
     DPhi1*Phit2(t2));
 
 MT(2,2)= (ht/6)*(DPhi2*Phit2(t1) + ...
     4*DPhi2*Phit2(t_md)+...
     DPhi2*Phit2(t2));


end

