function  Connectivities=Correct_Nested_Geometry_Parts(Points,Domain_Connectivites_per_Region)
%% function: CORRECT_NESTED_GEOMETRY_PARTS3D
% Algorithm:
% 1. step: Create a geometry information matrix that contains the information about which region
%  each domain belongs to and whether it is an internal or external domain (depending on the sign of the volume). 
% 2. step: if a region contains several external domains (as these are only connected at one edge), these are combined
% 3. step: extract outer hulls for each region i.e. the connectivity of the
% external domain
% 4. step: test if a point from another region lies inside of the hulls of the actual region 
% --> in this case the region is nested inside the other region
% 5. step: if multiple regions are nested inside one other region we have
% to test if they are also nested
% 6. step: determine in which subdomain the nested regions are nested,
% therefor we calculate the hulls of all possible nested subdomains
% 7. step: now we can correct the connectivities of the single subdomains,
% therefor we include the outer triangulations of an enclosed region in the
% triangulation of the specific subdomain
% 8. step: correct the geometry information matrix and delete the otherwise doubled
% entries in the domain connectivities
% 9. step: correct the connectivities of all regions and the geometry information matrix
% 10. step: combine all connectivites in a single list, correct the geometry information matrix and set the outer
% connectivity of all regions as the last entry
%
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

Number_of_regions=length(Domain_Connectivites_per_Region);

%% 1.step: create geometry information matrix
Reg_Dom_ExtInt=[];
for k=1:Number_of_regions
    for n=1:length(Domain_Connectivites_per_Region{k})
        Reg_Dom_ExtInt=[Reg_Dom_ExtInt;[k,n,...
            Determine_Orientation(Points,Domain_Connectivites_per_Region{k}{n})]];
    end
end

%% 2.step: combine multiple external domains
for k=1:Number_of_regions
    ind_act_region=find(Reg_Dom_ExtInt(:,1)==k);
    ind_external_domains=find(Reg_Dom_ExtInt(ind_act_region,3)==1);

    if(length(ind_external_domains)>1)
        for n=2:length(ind_external_domains)
            Domain_Connectivites_per_Region{k}{ind_external_domains(1)}=[Domain_Connectivites_per_Region{k}{ind_external_domains(1)};...
                Domain_Connectivites_per_Region{k}{ind_external_domains(n)}];
        end

        Reg_Dom_ExtInt(ind_act_region(ind_external_domains(2:end)),:)=[];
        Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,2)>ind_external_domains(2),2)=Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,2)>ind_external_domains(2),2)-1;
        Domain_Connectivites_per_Region{k}(ind_external_domains(2:end))=[];
    end
end

%% 3.step: extract outer hulls
tri_per_regions_outer_hull=cell(length(unique(Reg_Dom_ExtInt(:,1))),1);
for n=1:Number_of_regions
    entries_act_region=Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)==n,:);
    tri_per_regions_outer_hull{n}=Domain_Connectivites_per_Region{n}{entries_act_region(:,end)==1};
end

%% 4.step: test if regions are nested
% nested tells us which region is nested inside another region
% if (i,j)==1 --> region j is inside region i
nested_regions=zeros(Number_of_regions);
for n=1:Number_of_regions   
    inds_outer=tri_per_regions_outer_hull{n}(:);

    xmin_outer=min(Points(inds_outer,1)); xmax_outer=max(Points(inds_outer,1));
    ymin_outer=min(Points(inds_outer,2)); ymax_outer=max(Points(inds_outer,2));
    
    for k=1:Number_of_regions
        if(k==n)
            nested_regions(n,k)=0;
            continue
        else
            ind_P_test=tri_per_regions_outer_hull{k}(1,1);      
            
            if((Points(ind_P_test,1)>xmin_outer && Points(ind_P_test,1)<xmax_outer)...
                    && (Points(ind_P_test,2)>ymin_outer && Points(ind_P_test,2)<ymax_outer))

                zmin_outer=min(Points(inds_outer,3)); zmax_outer=max(Points(inds_outer,3));

                if((Points(ind_P_test,3)>zmin_outer && Points(ind_P_test,3)<zmax_outer))
                     % intersection possible since the bounding boxes
                     % overlap

                     inside_true_false=in_polyhedron(tri_per_regions_outer_hull{n},...
                         Points,Points(ind_P_test,:));

                     if(inside_true_false)
                         nested_regions(n,k)=1;
                     else
                         nested_regions(n,k)=0;
                     end
                else
                    nested_regions(n,k)=0;
                end
            else
                nested_regions(n,k)=0;
            end
        end
    end
end

%% 5. step: test if multiple regions are nested inside of another region
nested_sum=sum(nested_regions);
for n=1:Number_of_regions
    if(nested_sum(n)>1) % multiple regions are nested inside region n
        for k=1:Number_of_regions
            for j=1:Number_of_regions
                if(nested_regions(k,j)==1 && nested_regions(k,n)==1 && nested_regions(j,n)==1)
                    nested_regions(k,n)=0;
                end
            end
        end
    end
end

%% 6. step: test in which domains the regions are nested
num_nested_regions=sum(nested_regions,2)>0;
num_cons_per_region=zeros(Number_of_regions,1);
for n=1:Number_of_regions
    num_cons_per_region(n)=length(Domain_Connectivites_per_Region{n});
end

nested_domains=zeros(Number_of_regions,Number_of_regions,max(num_cons_per_region));
% nested_domains tells us which region is nested inside which subdomain
% if (i,k,j)==1 --> region j is inside subdomain k of region i

for n=1:Number_of_regions
    if(num_nested_regions(n)>0) 
        for k=1:Number_of_regions
            if(nested_regions(n,k)==1)                
                entries_act_region=Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)==n,:);
                vec_inds=1:length(entries_act_region(:,1)); % vec from 1 to #domains of the actual region
                inds_internal_domains=vec_inds(entries_act_region(:,end)==0); % internal domains of the actual region
                ind_P_test=tri_per_regions_outer_hull{k}(1,1);  
                
                for j=1:length(inds_internal_domains) 
                    % test for each internal domain                     
                    con_outer=Domain_Connectivites_per_Region{n}{entries_act_region(inds_internal_domains(j),2)};
                    inds_outer=con_outer(:);

                    xmin_outer=min(Points(inds_outer,1)); xmax_outer=max(Points(inds_outer,1));
                    ymin_outer=min(Points(inds_outer,2)); ymax_outer=max(Points(inds_outer,2));             

                    if((Points(ind_P_test,1)>xmin_outer && Points(ind_P_test,1)<xmax_outer)...
                            && (Points(ind_P_test,2)>ymin_outer && Points(ind_P_test,2)<ymax_outer))

                        zmin_outer=min(Points(inds_outer,3)); zmax_outer=max(Points(inds_outer,3));
                        
                        if((Points(ind_P_test,3)>zmin_outer && Points(ind_P_test,3)<zmax_outer))
                            % intersection possible since the bounding boxes
                            % overlap

                            inside_true_false=in_polyhedron(con_outer,...
                                Points,Points(ind_P_test,:));

                            if(inside_true_false)
                                nested_domains(n,k,inds_internal_domains(j))=1;
                            else
                                nested_domains(n,k,inds_internal_domains(j))=0;
                            end
                        else
                            nested_domains(n,k,inds_internal_domains(j))=0;
                        end
                    else
                        nested_domains(n,k,inds_internal_domains(j))=0;
                    end                           
                end
            else
                continue
            end
        end
    end
end

%% 6. step: correct the connectivities of all domains
for n=1:Number_of_regions
    for k=1:Number_of_regions
        for j=1:max(num_cons_per_region)
            if(nested_domains(n,k,j)==1)
                Domain_Connectivites_per_Region{n}{j}=[Domain_Connectivites_per_Region{n}{j};...
                    tri_per_regions_outer_hull{k}];
            end
        end
    end
end

%% 7.step: correct geometry information matrix and delete the double entries from the connectivity list
for n=1:Number_of_regions
    for k=1:Number_of_regions
        if(nested_regions(n,k)==1)
            % delete unused connectivity
            Domain_Connectivites_per_Region{k}(Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)==k,3)==1)=[];
            
            % index of the removed conenctivity
            ind_removed_external_domain=Reg_Dom_ExtInt(:,1)==k & Reg_Dom_ExtInt(:,3)==1;

            % domain number of the removed conenctivity
            dom_num_removed_external_domain=Reg_Dom_ExtInt(ind_removed_external_domain,2);
            
            % reduce all domain numbers that are larger than the removed
            % domain number by one
            Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)==k & Reg_Dom_ExtInt(:,2)>dom_num_removed_external_domain,2)=...
                Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)==k & Reg_Dom_ExtInt(:,2)>dom_num_removed_external_domain,2)-1;
            
            % remove the entry from the geometry information matrix
            Reg_Dom_ExtInt(ind_removed_external_domain,:)=[];
        end
    end
end

%% 8. step: correct the connectivities of all regions
for n=Number_of_regions:-1:1
    for k=Number_of_regions:-1:1
        if(nested_regions(n,k)==1)
            Domain_Connectivites_per_Region{n}=...
                [Domain_Connectivites_per_Region{n};...
                 Domain_Connectivites_per_Region{k}];

            entries_act_region=Reg_Dom_ExtInt(:,1)==k;
            entries_test_region=Reg_Dom_ExtInt(:,1)==n;
            
            number_subdomains=max(Reg_Dom_ExtInt(entries_test_region,2));
            vec_new_nums=(number_subdomains+1):(number_subdomains+sum(entries_act_region));
            
            Reg_Dom_ExtInt(entries_act_region,1)=n;
            Reg_Dom_ExtInt(entries_act_region,2)=vec_new_nums;
        end
    end
end

num_remaining_regions=length(unique(Reg_Dom_ExtInt(:,1)));
for n=1:num_remaining_regions-1
    m=min(Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)>n,1));
    inds=m==Reg_Dom_ExtInt(:,1);
    Reg_Dom_ExtInt(inds,1)=m-(m-n)+1;
end

% remove unnecessary region connectivities 
delete=sum(nested_regions)==1;
Domain_Connectivites_per_Region(delete)=[];

% final number of regions
Number_of_regions=length(Domain_Connectivites_per_Region);

% correct the region numbers in the geometry information matrix
reg_entries=unique(Reg_Dom_ExtInt(:,1));
for n=1:Number_of_regions
    vec_Reg_Dom_ExtInt(Reg_Dom_ExtInt(:,1)==reg_entries(n),1)=n;
end
Reg_Dom_ExtInt(:,1)=vec_Reg_Dom_ExtInt;


%% 9. step: combine connectivities
k=1;
num_internal_subdomains=sum(Reg_Dom_ExtInt(:,3)==0);
Connectivities=cell(num_internal_subdomains+1,1);
for n=1:length(Reg_Dom_ExtInt(:,1))
   if(Reg_Dom_ExtInt(n,3)==1)
       Connectivities{num_internal_subdomains+1}=[Connectivities{num_internal_subdomains+1};...
           Domain_Connectivites_per_Region{Reg_Dom_ExtInt(n,1)}{Reg_Dom_ExtInt(n,2)}];
   else
       Connectivities{k}=Domain_Connectivites_per_Region{Reg_Dom_ExtInt(n,1)}{Reg_Dom_ExtInt(n,2)};
       k=k+1;
   end   
end
end

function [External_true_false,V]=Determine_Orientation(Points,Connectivity)
%% function DETERMINE_ORIENTATION3D 
% determines the orientation and the volume of the triangulation
%
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

% all points of all vertices
x0=Points(Connectivity(:,1),1); x1=Points(Connectivity(:,2),1); x2=Points(Connectivity(:,3),1);
y0=Points(Connectivity(:,1),2); y1=Points(Connectivity(:,2),2); y2=Points(Connectivity(:,3),2);
z0=Points(Connectivity(:,1),3); z1=Points(Connectivity(:,2),3); z2=Points(Connectivity(:,3),3);

d12= [x0-x2,y0-y2,z0-z2]; d13= [x1-x2,y1-y2,z1-z2];
n = 0.5*CrossVectorized(d12,d13);

%calculate volume
V = 1/3*sum(dot([x0,y0,z0],n));

% test if external or not
External_true_false=false;
if(V<0)
    % normals point inward --> external domain
    External_true_false=true;
end
end

function N=CrossVectorized(a,b)
%% function: CROSS3D
% is a vectorized version of the cross product
%
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-11-2023
% License: MIT License

N = [a(:,2).*b(:,3) - a(:,3).*b(:,2)...
     a(:,3).*b(:,1) - a(:,1).*b(:,3)...
     a(:,1).*b(:,2) - a(:,2).*b(:,1)]; 
end
