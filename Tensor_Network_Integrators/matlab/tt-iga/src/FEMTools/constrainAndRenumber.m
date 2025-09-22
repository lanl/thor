function IndexE_new = constrainAndRenumber(IndexE, In1, In2)
%CONSTRAINANDRENAMBER  Apply master–slave constraints and then renumber DOFs
%
%   [IndexE_new, old2new] = constrainAndRenumber(IndexE, In1, In2)
%
%   Inputs:
%     IndexE   : [ne×npe] element-to-DOF connectivity array
%     In1,In2  : vectors of the same length, giving pairs
%                master DOF = In1(k), slave DOF = In2(k)
%
%   Outputs:
%     IndexE_new : same size as IndexE, with constraints applied and
%                  DOFs renumbered 1..Nnew
%     old2new    : a vector mapping original DOF i to new DOF old2new(i)
%
%   Behavior:
%     1) For each k, replaces every occurrence of slave In2(k) by In1(k).
%     2) Collects the set of all used DOFs in IndexE, sorts them,
%        and assigns them new labels 1..Nnew in ascending order.
%     3) Remaps IndexE to these new labels, returning IndexE_new.

  assert(numel(In1)==numel(In2), 'In1 and In2 must be same length');

  % 1) apply master<-slave replacement
  IndexE_constrained = IndexE;
  for k = 1:numel(In1)
    m = In1(k);  s = In2(k);
    IndexE_constrained(IndexE_constrained==s) = m;
  end

  % 2) find the sorted list of used DOFs
  used = unique(IndexE_constrained(:));
  % build old-to-new map (preallocate to max original DOF)
  maxOld = max(used);
  old2new = zeros(maxOld,1);
  old2new(used) = 1:numel(used);

  % 3) remap connectivity
  IndexE_new = old2new(IndexE_constrained);

end
