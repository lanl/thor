function att = convert_mesh_to_tt_fn(a,att,tol)

% att = TTMesh();
prop = properties(a);

for i = 1:length(prop)
  if ~isempty(a.(prop{i})) && ~isscalar(a.(prop{i})) && ~isvector(a.(prop{i}))
    % prop{i}
    att.(prop{i}) = tt_tensor(a.(prop{i}),tol);
  else
    att.(prop{i}) = a.(prop{i});
  end
end