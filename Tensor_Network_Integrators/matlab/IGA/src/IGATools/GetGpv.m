function [Gpv1, Gpv2] = GetGpv(Gp1, Gp2)    
    len1 = length(Gp1);
    len2 = length(Gp2);
    Gpv1 = repelem(Gp1(:), len2);
    Gpv2 = [Gp2(:); Gp1(:)];
end
