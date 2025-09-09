function [DMA] = dmchgl(N,ET)
% Compute the derivative 
if (N==1)
    DMA(1,1)=0;
end
DMA=zeros(N+1,N+1);
if (N>1)
    DN=N;
    CN=(2*DN*DN+1)/6;
    if mod(N,2)==0 % cioe' N e' pari
              N2=N/2; % Dibyendu; Original N/2-1; 
              SN=1+4*N2-2*N;
              
    else   % cioe' N e' dispari
              N2 = (N-1)/2;
              SN=1+4*N2-2*N;
    end 
    
   % SN=1+4*(N/2)-2*N;
    DMA(1,1)=-CN;
    DMA(N+1,N+1)=CN;
    DMA(1,N+1)=-0.5*SN;
    DMA(N+1,1)=0.5*SN;
    
    if (N==2)
        return;
    end
    
    SGN=-1;
    
    for J=2:N
        EJ=ET(J);
        DMA(1,J)=(-2*SGN)/(1+EJ);
        DMA(N+1,J)=(2*SGN*SN)/(1-EJ);
        DMA(J,1)=(.5*SGN)/(1+EJ);
        DMA(J,N+1)=-(.5*SGN*SN)/(1-EJ);
        SGN=-SGN;
    end
    
    SGNI=-1;
    
    for I=2:N
        EI=ET(I);
        SGNJ=-1;
        for J=2:N
            if (I ~=J)
                EJ=ET(J);
                DMA(I,J)=(SGNI*SGNJ)/(EI-EJ);
            else
                DMA(I,J)=-(.5*EI)/(1-EI*EI);
            end
            SGNJ=-SGNJ;
        end
        SGNI=-SGNI;
    end
end

end

