function Nout = BCIndexM(NI,BCI)
Nout=NI;
for k=1:length(BCI)
    Nout=BCIndex(Nout,BCI(k)-k+1);
end

    function Nout = BCIndex(NI,BCI)
        [n1,n2]=size(NI);
        Nout=zeros(n1,n2);
        for i=1:n1
            for j=1:n2
                if NI(i,j)<BCI
                    Nout(i,j)=NI(i,j);
                end
                if NI(i,j)==BCI
                    Nout(i,j)=-1;
                end
                if NI(i,j)>BCI
                    Nout(i,j)=NI(i,j)-1;
                end
            end
        end
    end
end