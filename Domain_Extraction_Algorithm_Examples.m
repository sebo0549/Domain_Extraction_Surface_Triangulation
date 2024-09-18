function Domain_Extraction_Algorithm_Examples
%% function: DOMAIN_EXTRACTION_ALGORTIHM_EXAMPLES
% extracts all contained domains from a three-dimensional triangulation
%
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 18-05-2024
% License: MIT License

close all
addpath 'Functions' 'Meshfiles' % add folders to path

%% Parameter
file_name='6_domains_4_regions_nested_two_times'; % file name of the triangulation in the 'Meshfiles' folder
plot_geom_true_false=true; % if true: all domains are plotted in a single subplot.
% Should be used with care: gets slow if the number of domains is high!

fprintf('Extracting domains from the mesh file: ''%s.stl''\n\n',file_name)

% import mesh file
TR = stlread(sprintf('%s.stl',file_name));
Num_Elements=length(TR.ConnectivityList(:,1));

% store coordinates and initial connectivity
Points=TR.Points;
Initial_Connectivity=TR.ConnectivityList;

tic
% determine the regions
Region_Connectivities=Extract_Regions(Initial_Connectivity);

% extract the domains
Domain_Connectivites_per_Region=Extract_Domains(Points,Region_Connectivities);

% correct nested geometry parts
Connectivities=Correct_Nested_Geometry_Parts(TR.Points,Domain_Connectivites_per_Region);
t=toc;

% print runtime
fprintf('Geometry consists of %1.0f domains\n',length(Connectivities))
fprintf('Time needed: t = %1.3fs\n',t)
fprintf('Number of elements: M = %1.0f\n\n',Num_Elements)

if(plot_geom_true_false)
    % plot the initial geometry and all extracted domains
    Num_Doms=length(Connectivities);
    Colors=lines(Num_Doms);

    Xmin=min(Points(:,1)); Ymin=min(Points(:,2)); Zmin=min(Points(:,3));
    Xmax=max(Points(:,1)); Ymax=max(Points(:,2)); Zmax=max(Points(:,3));

    figure()
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    hold on
    p=numSubplots(Num_Doms);
      
    for n=1:Num_Doms
        subplot(p(1),p(2),n)
        hold on
        trisurf(Initial_Connectivity,...
            Points(:,1),Points(:,2),Points(:,3),...
            'FaceColor','none','EdgeAlpha',0.25)
        trisurf(Connectivities{n},...
            Points(:,1),Points(:,2),Points(:,3),...
            'FaceColor',Colors(n,:),...
            'FaceAlpha',0.4,'EdgeColor',Colors(n,:),'EdgeAlpha',1)
        axis equal
        set(gca,'XLim',1.1*[Xmin,Xmax],...
            'YLim',1.1*[Ymin,Ymax],'ZLim',1.1*[Zmin,Zmax])
        grid on
        view(30,10)
    end
end
end
