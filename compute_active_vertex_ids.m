function vertex_ids = compute_active_vertex_ids(is_vertex_active)

% 1. initialize
vertex_ids = is_vertex_active*0;

% 2. find active vertices
loc = find(is_vertex_active == 1);

% 3. set IDs
vertex_ids(loc) = [1:length(loc)];
