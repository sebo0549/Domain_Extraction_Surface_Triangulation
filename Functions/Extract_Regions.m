function Connectivities_Regions=Extract_Regions(Initial_Connectivity)
%% function: EXTRACT_REGIONS3D
% determines the number and the connectivities of all
% regions of the geometry
%
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

% create planar graph from the initial connectivity
g = graph(Initial_Connectivity, Initial_Connectivity(:, [2 3 1]));

% get region index for each vertex
inds_vertices_regions=conncomp(g);
num_regions=length(unique(inds_vertices_regions));

Connectivities_Regions=cell(num_regions,1);
vec_vertices=1:length(inds_vertices_regions);
for n=1:num_regions
    inds_con=sum(ismember(Initial_Connectivity,vec_vertices(inds_vertices_regions==n)),2)==3;
    Connectivities_Regions{n}=Initial_Connectivity(inds_con,:);
end
end