function zisoi = compute_vertical_interface_depth_for_elm()

scalez = 0.025;
nlevgrnd = 15;

% node depth
j = [1:nlevgrnd];
zsoi = scalez*(exp(0.5*(j-0.5))-1.);

% z at interface
zisoi(1)          = 0.;
zisoi(2:nlevgrnd) = 0.5*(zsoi(1:end-1)+zsoi(2:end));
zisoi(nlevgrnd+1) = zsoi(nlevgrnd) + 0.5*(zsoi(nlevgrnd)-zsoi(nlevgrnd-1));

