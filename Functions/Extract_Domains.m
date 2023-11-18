function Domain_Connectivites_per_Region=Extract_Domains(Points,Connectivities_regions)
%% function: EXTRACT_DOMAINS3D 
% returns a cell array of cell arrays containing
% all connectivity matrices of all domains and the number of the
% domains in all regions
%
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

Domain_Connectivites_per_Region=cell(length(Connectivities_regions),1);
for k=1:length(Connectivities_regions)
    Domain_Connectivites_per_Region{k}=Extract_Domains_per_Region(Points,Connectivities_regions{k});
end
end

function  Domain_Connectivities=Extract_Domains_per_Region(Points,Connectivity)
%% function: EXTRACT_DOMAINS_PER_REGION3D 
% determines the connectivity of all domains for the
% given region
% 
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

% since we only handle the connectivity of a region suppress the warning 
% that not all vertices are referenced by the triangulation
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')

% build a matrix with all edge connectivities of the triangulation
all_edges=[Connectivity(:,1),Connectivity(:,2);
           Connectivity(:,2),Connectivity(:,3);
           Connectivity(:,3),Connectivity(:,1)];
all_edges=unique(sort(all_edges,2),'rows');

% save some useful information about the triangulation
num_Vertices=length(Points(:,1));
num_Edges=length(all_edges(:,1));
num_Triangles=length(Connectivity(:,1));
anz_wedges=6*num_Triangles;

%% 1. phase - setting up the wedge matrix
% extract all neighboring triangles for each edge
ID = edgeAttachments(triangulation(Connectivity,Points),all_edges(:,1),all_edges(:,2));

% get all neighbor vertices for each edge
neighbors=cell(num_Edges,1);
for n=1:num_Edges
    Mat_Neighbor_Triangles=Connectivity(ID{n},:);
    neighbors{n}=Mat_Neighbor_Triangles(Mat_Neighbor_Triangles~=all_edges(n,2) & Mat_Neighbor_Triangles~=all_edges(n,1));
end

% rotate and sort the neighbor vertices of each edge
edge_vecs=(Points(all_edges(:,2),:)-Points(all_edges(:,1),:))';
[theta,phi,~]=cart2sph(edge_vecs(1,:),edge_vecs(2,:),edge_vecs(3,:));
sorted_list=zeros(num_Triangles*3,4); j=1; 
for n=1:num_Edges
    % get rotation matrix
    R=rotYZ(phi(n),theta(n)); 
    
    act_neighbors=neighbors{n};
    act_edge=all_edges(n,:);
    
    points_shifted=(Points(act_neighbors,:)-Points(act_edge(1),:))';
    for k=1:length(neighbors{n}(:))
        % rotate all neighbors
        rotated_point=R*points_shifted(:,k);

        % get angle in the x,y-plane
        angle=atan2(rotated_point(3),rotated_point(2))+ 2*pi*(rotated_point(3)<0);
        
        % store edge connectivity, neighbors and angles in the sorted list
        % matrix
        sorted_list(j,1)=act_edge(1);       
        sorted_list(j,2)=act_edge(2);
        sorted_list(j,3)=act_neighbors(k);  
        sorted_list(j,4)=angle;
        j=j+1;
    end
end

% sort the list of edges,neighbors and angles to make the extration of the
% list per vertex easier
sorted_list=sortrows(sorted_list,[1,2,4]);

sorted_list_per_vertex=cell(num_Vertices,1);
[Group_Counts,Group_Nums] = groupcounts(sorted_list(:,1));
inds_sorted_per_vertex=[0,cumsum(Group_Counts')];
for n=1:length(Group_Nums)
    sorted_list_per_vertex{Group_Nums(n)}=sorted_list(inds_sorted_per_vertex(n)+1:inds_sorted_per_vertex(n+1),:);
end

% set up the wedge matrix
j=1;
all_wedges=zeros(num_Triangles*3,4);
for n=1:num_Edges
    actual_edge=all_edges(n,:);
    actual_list=sorted_list_per_vertex{actual_edge(1)};

    inds=actual_list(:,2)==actual_edge(2);
    num_Edges_act=sum(inds);
    temp_list=actual_list(inds,:);
    
    for m=1:num_Edges_act
        all_wedges(j,:)=[temp_list(mod(m,num_Edges_act)+1,3) temp_list(m,1:3)];
        j=j+1;
    end
end

% appending the wedge matrix with permuted columns
all_wedges=[all_wedges;all_wedges(:,[4,3,2,1])];

%% 2. phase - extraction of the domain connectivities
% sort wedge matrix to make it easier to extract the wedges per vertex 
all_wedges=sortrows(all_wedges,1); 

[Group_Counts,Group_Nums] = groupcounts(all_wedges(:,1));
inds_sorted_wedges_per_vertex=[0,cumsum(Group_Counts')];
sorted_list_wegdes_per_vertex=cell(num_Vertices,1);
ind_left_per_vertex=cell(num_Vertices,1);
for n=1:length(Group_Nums)
    sorted_list_wegdes_per_vertex{Group_Nums(n)}=all_wedges(inds_sorted_wedges_per_vertex(n)+1:inds_sorted_wedges_per_vertex(n+1),:);
    ind_left_per_vertex{Group_Nums(n)}=1:Group_Counts(n);
end

% initialize used arrays and variables
domain_wedges=[]; n=1; last_hit=1;
ind_non_empty=1:num_Vertices; 
non_empty_true_false=true(num_Vertices,1); 
ind_empty_changed=false;
while(true)
    % check if there are some unchecked wedges left --> if true: use this
    % wedge as the first wedge of the new domain connectivity
    new_domain_available=false; j=1; 
    for k=ind_non_empty
        if(~isempty(ind_left_per_vertex{k}))
            W=sorted_list_wegdes_per_vertex{k}(ind_left_per_vertex{k}(1),:);
            ind_left_per_vertex{k}(1)=[];
            new_domain_available=true;
            break;
        else
            non_empty_true_false(j)=false;
            ind_empty_changed=true;
        end
        j=j+1;
    end
    if(ind_empty_changed)
        ind_non_empty=ind_non_empty(non_empty_true_false);
        ind_empty_changed=false;
        non_empty_true_false=true(length(ind_non_empty),1);
    end
    
    if(~new_domain_available)
        % all wedges have been checked and there is no further domain left 
        break;
    end
    
    % find all continuous wedges
    subset_first_entry=W; new_domain=false; j=1; dom_wedges=zeros(anz_wedges,4); dom_wedges(1,:)=W;
    while(~new_domain)
        while(true)
            ind_list=ind_left_per_vertex{W(3)};
            wedge_list=sorted_list_wegdes_per_vertex{W(3)}(ind_list,:);
            
            ind_new= find(W(2)==wedge_list(:,2) & W(4)==wedge_list(:,3));
            W_new=wedge_list(ind_new,:);

            ind_left_per_vertex{W(3)}(ind_new)=[];
            W=W_new;
            
            j=j+1;
            dom_wedges(j,:)=W;
            
            if(subset_first_entry(1)==W_new(3) && subset_first_entry(2)==W_new(2) && subset_first_entry(3)==W_new(4))
                % closed loop in the wedges discovered
                break;
            end
        end
        
        cont=true;
        while(cont)
            % test whether there is another wedge that has not been checked yet 
            for k=last_hit:j
                A=dom_wedges(k,:);
                ind_list=ind_left_per_vertex{A(2)};
                wedge_list=sorted_list_wegdes_per_vertex{A(2)}(ind_list,:);
                ind_new=find(A(3)==wedge_list(:,2) & A(1)==wedge_list(:,3));
                
                if(~isempty(ind_new))
                    last_hit=k;
                    cont=false;
                    
                    W=wedge_list(ind_new,:);
                    
                    subset_first_entry=W;
                    j=j+1;
                    dom_wedges(j,:)=W;
                    
                    ind_left_per_vertex{A(2)}(ind_new)=[];
                    break;
                end
            end
            
            if(isempty(ind_new))
                % there is no further continuous wedge 
                % --> all continous wedges of the domain has been successfully found
                domain_wedges{n}=dom_wedges(1:j,:);
                dom_wedges=[];
                last_hit=1;
                n=n+1;
                cont=false;
                new_domain=true;
            end
        end
    end
end

% the domain connectivities are given by the first three columnes of all
% domain wedge lists
Domain_Connectivities=cell(length(domain_wedges),1);
for n=1:length(domain_wedges)
    Domain_Connectivities{n}=domain_wedges{n}(:,1:3);
    [~,i]=unique(sort(Domain_Connectivities{n},2),'rows');
    Domain_Connectivities{n}=Domain_Connectivities{n}(i,:);
end
end

function R=rotYZ(alpha,beta)
%% function: ROTYZ  
% returns thecombined rotation matrix for the rotation arround the y- and
% the z-axis
% rotate ccw arround z-axis -> theta is defined ccw regarding the
% z-axis
% rotate ccw arround y-axis -> phi is defined ccw regarding
% the y axis
% 
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

R=[cos(alpha)*cos(beta) cos(alpha)*sin(beta) sin(alpha);...
    -sin(beta) cos(beta) 0;...
    -sin(alpha)*cos(beta) -sin(alpha)*sin(beta) cos(alpha)];
end