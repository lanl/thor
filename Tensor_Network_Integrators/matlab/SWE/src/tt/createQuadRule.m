function [x,w] = createQuadRule(Nq)
    %
    if(Nq==1)
        x = [0];
        w = [2];
    elseif(Nq==2)
        x = [-1 1]/sqrt(3);
        w = [ 1 1];
    elseif(Nq==3)
        x = [-sqrt(3/5) 0 sqrt(3/5)];
        w = [5/9 8/9 5/9];
    else
        error("Number of quaderature points must be less than 4")
    end
    % shift the quadrature points such that $x \in (0,1)$
    x = x'/2 + 0.5;
    % adjust the quadrature weights such that sum(w)=1
    w = w'/2;
    %
end