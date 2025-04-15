function X = cut_and_pad(X,dims,adding)
%%Inputs%%%%%%%
%X - tensor we want to cut and pad the cut dimentions with zeros
%dims - a cell list that contains which col we want to preserve or pad
%        ex: dims ={(1:nx-1),(2:ny),(2:nz-1)} 
%adding - how many col do we want to add per core
%        ex: adding = [0,2,2] adding zero to G{1} adding 2 to G{2} and
%        adding 2 to G{3}, where the new ones are at the begining and the
%        end
%%Outputs%%%%%%
%X - cut and padded with zeros

X{1} = X{1}(:,dims{1},:);
X{2} = X{2}(:,dims{2},:);
X{3} = X{3}(:,dims{3},:);
for j = 1: numel(dims)
    G = core2cell(X);
    [r1,a2,r2] = size(G{j});
    Gforward = zeros([r1,a2 + adding(j),r2]);
    Gforward(:,dims{j},:)= G{j}; 
    G{j} = Gforward;
    X = cell2core(X,G);
end
end