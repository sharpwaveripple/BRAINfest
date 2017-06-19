function programmatic_brain

close

%MAKE LINES TRANSLUCENT? add info about half brain; add front back rotator;
%legend in bottom
%make adjacency matrix
V = [0	1	0	1	0	0	0	0	1	0	0	0	0	1	0	0	0	0	1	1	0	0	1	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	0	0	0	0];

info

M = triu(ones(10)); %based on number of nodes; way to calculate this easier I'm sure
M(M==1) = V;
A = M+M'-diag(diag(M)); %can use this without vector method
B = any(A); %check to see if any elements = 0
D = distance_bin(A);
d1 = D(6:9,6:9); %motor network
d2 = D(3:6,3:6); %visual network
d3 = D(1:3,1:3); %language network
% d3 = D(1:3,1:3)
[lambda,E,~,~,d] = charpath(D,0,1);
[~,e1] = charpath(d1,0,1);
[~,e2] = charpath(d2,0,1);
[~,e3] = charpath(d3,0,1);

maxE=E;
maxe1=e1;
maxe2=e2;
maxe3=e3;

%load nodes - x, y, z, colour, size, label
fid=fopen('AAL10');
data = textscan(fid,'%f %f %f %f %f %s','CommentStyle','#');
fclose(fid);
surf.sphere = [cell2mat(data(1)) cell2mat(data(2)) cell2mat(data(3)) cell2mat(data(4)) cell2mat(data(5))];
surf.sphere(:,6) = 4;
surf.label = data{6};
surf.nsph = size(surf.sphere,1); %number nodes
labels = surf.label.';

f = figure('Name','Simulate a stroke','Visible','on',...
    'Position',[0 0 800 600]);
%position = [left bottom width height]

%  Construct the components.
p = uipanel(f,'Title','Select area(s)','Position',[.05 .33 .15 .4],...
    'BackgroundColor',[1 1 1],'FontSize',12);

lb = uicontrol(p,'Style','listbox','String',labels,...
    'Max',2,'Min',0,'Value',[],'Units','normalized',...
    'Position', [.02 .18 .95 .8],'BackgroundColor',[1 1 1],...
    'Callback',@lb_Callback,'FontSize',10);

p2 = uipanel(f,'Title','Network efficiency','Position',[.05 .1 .31 .2],...
    'BackgroundColor',[1 1 1],'FontSize',12);

p3 = uipanel(f,'Title','Area information','Position',[.21 .33 .15 .4],...
    'BackgroundColor',[1 1 1],'FontSize',12);

p4 = uipanel(f,'Title','How to play','Position',[.05 .75 .31 .2],...
    'BackgroundColor',[1 1 1],'FontSize',12);

visualheader = uipanel(p2,'Title','Visual Network','Position',[.03 .55 .3 .4],...
    'BackgroundColor',[1 1 1],'FontSize',11);

motorheader = uipanel(p2,'Title','Motor Network','Position',[.35 .55 .3 .4],...
    'BackgroundColor',[1 1 1],'FontSize',11);

languageheader = uipanel(p2,'Title','Language Network','Position',[.67 .55 .3 .4],...
    'BackgroundColor',[1 1 1],'FontSize',11);

brainheader = uipanel(p2,'Title','Whole brain','Position',[.03 .07 .94 .4],...
    'BackgroundColor',[1 1 1],'FontSize',11);

nodeinfo = uicontrol(p3,'Style','text','String',caudate,...
    'HorizontalAlignment','left','BackgroundColor',[1 1 1],...
        'Units','normalized','Position',[.02 .04 .95 .95]);
    
instructiontext = uicontrol(p4,'Style','text','String',instructions,...
        'HorizontalAlignment','left','FontSize',10,'FontName','Helvetica',...
        'BackgroundColor',[1 1 1],...
        'Units','normalized','Position',[.02 .04 .95 .95]);

motortext = uicontrol(motorheader,'Style','text','String',sprintf('%d%%',e1/maxe1*100),...
    'Units','normalized','Position',[.25 .3 .5 .5],...
    'BackgroundColor',[1 1 1],'FontSize',11);

visualtext = uicontrol(visualheader,'Style','text','String',sprintf('%d%%',e2/maxe2*100),...
    'Units','normalized','Position',[.25 .3 .5 .5],...
    'BackgroundColor',[1 1 1],'FontSize',11);

languagetext = uicontrol(languageheader,'Style','text','String',sprintf('%d%%',e3/maxe3*100),...
    'Units','normalized','Position',[.25 .3 .5 .5],...
    'BackgroundColor',[1 1 1],'FontSize',11);

braintext = uicontrol(brainheader,'Style','text','String',sprintf('%d%%',E/maxE*100),...
    'Units','normalized','Position',[.25 .4 .5 .5],...
    'BackgroundColor',[1 1 1],'FontSize',11);

resetbutton = uicontrol(p,'Style','pushbutton','String','Clear selection',...
    'Units','normalized','Position',[.02 .03 .95 .12],...
    'Callback',@pb_Callback);

ax = axes('Units','normalized','Position',[.3 .2 .8 .7]);

plot_mesh

[surf.ncyl,surf.cylinder] = prep_network(A);

hold on
for i=1:surf.nsph
    draw_sphere(surf,i,B);
end

for j=1:surf.ncyl
    draw_line(surf,j);
end

axis tight; axis vis3d off;daspect([1 1 1]);
eval(['material ','dull',';']);eval(['lighting ','gouraud',';']);
cam = camlight('headlight','infinite');
set(gcf,'color','w');

hold off

    function pb_Callback(hObject,eventdata)
        set(lb,'Value',[]);
        cla(ax);
        plot_mesh
        B = any(A); 
        [surf.ncyl,surf.cylinder] = prep_network(A);

        hold on
        for i = 1:surf.nsph
            draw_sphere(surf,i,B);
        end
        
        for j=1:surf.ncyl
            draw_line(surf,j);
        end
        
        axis tight; axis vis3d off;daspect([1 1 1]);
        eval(['material ','dull',';']);eval(['lighting ','flat',';']);
        cam = camlight('headlight');
        
        hold off
        
        set(braintext,'String',sprintf('%d%%',E/maxE*100));
        set(motortext,'String',sprintf('%d%%',e1/maxe1*100));
        set(visualtext,'String',sprintf('%d%%',e2/maxe2*100));
        set(languagetext,'String',sprintf('%d%%',e3/maxe3*100));


    end

    function lb_Callback(hObject,eventdata)
        Aprime = A;
        h = hObject.Value;
        Aprime(:,h) = 0;
        Aprime(h,:) = 0;
        Dprime = distance_bin(Aprime);
        d1prime = Dprime(6:9,6:9);
        d2prime = Dprime(3:6,3:6);
        d3prime = Dprime(1:3,1:3);
        [lambda,Eprime,~,~,d] = charpath(Dprime,0,1);
        [~,e1prime] = charpath(d1prime,0,1);
        [~,e2prime] = charpath(d2prime,0,1);    
        [~,e3prime] = charpath(d3prime,0,1);      

        B = any(Aprime);

        set(braintext,'String',sprintf('%d%%',round(Eprime/maxE*100)));
        set(motortext,'String',sprintf('%d%%',round(e1prime/maxe1*100)));
        set(visualtext,'String',sprintf('%d%%',round(e2prime/maxe2*100)));
        set(languagetext,'String',sprintf('%d%%',round(e3prime/maxe3*100)));

        cla(ax);
        
        plot_mesh
        [surf.ncyl,surf.cylinder] = prep_network(Aprime);
        
        hold on
        for i = 1:surf.nsph
            draw_sphere(surf,i,B);
        end
        
        for j=1:surf.ncyl
            draw_line(surf,j);
        end
        
        axis tight; axis vis3d off;daspect([1 1 1]);
        eval(['material ','dull',';']);eval(['lighting ','flat',';']);
        cam = camlight('headlight');
        
        hold off
        
    end

    function [ncyl,cylinder] = prep_network(A)
        index = find(triu(A) == 1);
        ncyl = length(index);
        cylinder=zeros(ncyl,7);
        [cylinder(:,1), cylinder(:,2)] = ind2sub(size(A),index);
        cylinder(:,3) = 1;
        cylinder(:,4) = 2; %size
        cylinder(:,5) = 1; %colour
        cylinder(:,6) = .7; %opacity
        cylinder(:,7) = 1;
    end

    function plot_mesh
        SurfFileName='BrainMesh_ICBM152Right.nv';
        fid=fopen(SurfFileName);
        data = textscan(fid,'%f','CommentStyle','#');
        surf.vertex_number = data{1}(1);
        surf.coord  = reshape(data{1}(2:1+3*surf.vertex_number),[3,surf.vertex_number]);
        surf.ntri = data{1}(3*surf.vertex_number+2);
        surf.tri = reshape(data{1}(3*surf.vertex_number+3:end),[3,surf.ntri])';
        Brain=trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),'EdgeColor','none');
%         whitebg(gcf,[1,1,1]);
%         set(gcf,'Color',[1,1,1],'InvertHardcopy','off');
        eval(['material ','dull',';']); eval(['shading ','interp',';']);axis off
        set(Brain,'FaceColor',[.95,.95,.95]);
        set(Brain,'FaceAlpha',.3);
        daspect([1 1 1])
        view(90,0);
    end

    function draw_sphere(surf,i,B)
        
        switch surf.sphere(i,4)
            case 1
                c = [1 1 0]; %yellow
            case 2
                c = [1 0 1]; %magenta
            case 3
                c = [0 1 1]; %cyan
            case 4
                c = [0 1 0]; %green
        end
        
        if B(1,i) == 0
            c = [1 0 0];
        end
        
        t = surf.sphere(i,6); %size of node
        [x,y,z]=sphere(100); %matrices of unit sphere
        x=x.*t+surf.sphere(i,1); %scales sphere by size, then adds coordinates
        y=y.*t+surf.sphere(i,2);
        z=z.*t+surf.sphere(i,3);
        Node=mesh(x,y,z,'EdgeColor','none');
        set(Node,'FaceColor',c);
        set(Node,'EdgeAlpha',0)
        PlotLabel(surf,i);
        
    end

    function PlotLabel(surf,i)
        text(surf.sphere(i,1),surf.sphere(i,2),surf.sphere(i,3)+surf.sphere(i,4)+5,surf.label{i},...
            'FontName','Arial','FontWeight','bold','FontAngle','normal','FontSize',11,...
            'FontUnits','points','HorizontalAlignment','center');
    end

    function draw_line(surf,j)
        sphere = surf.sphere;
        %the length of the cylinder is found by subtracting the distance between
        %the spheres at specified coordinate points
        length_cyl=norm(sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3))-sphere(surf.cylinder(j,2),5)-sphere(surf.cylinder(j,1),5); %change 5s to 7s
        det = 5; %detail; determines how many points on the line to make
        n = surf.cylinder(j,4);
        theta = (0:det) / det * 2 * pi;
        sintheta = sin(theta);
        sintheta(det + 1) = 0;
        n = ones(100,1) * 0.5 * n ;
        x = n * cos(theta);
        y = n * sintheta;
        w = length(n);
        z = (0:w-1)'/(w-1) * ones(1,det + 1);
        Line = mesh(x,y,z * length_cyl);
        unit_Vx=[0 0 1];
        angle_X1X2=acos( dot( unit_Vx,sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3) )/( norm(unit_Vx)*norm(sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3))) )*180/pi;
        axis_rot=cross([0 0 1],(sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3)) );
        if angle_X1X2~=0
            rotate(Line,axis_rot,angle_X1X2,[0 0 0])
        end
        set(Line,'XData',get(Line,'XData')+sphere(surf.cylinder(j,1),1) + (sphere(surf.cylinder(j,2),1) -sphere(surf.cylinder(j,1),1))/norm(sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3))*sphere(surf.cylinder(j,1),5)); %change 5's at end to 7s
        set(Line,'YData',get(Line,'YData')+sphere(surf.cylinder(j,1),2) + (sphere(surf.cylinder(j,2),2) -sphere(surf.cylinder(j,1),2))/norm(sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3))*sphere(surf.cylinder(j,1),5));
        set(Line,'ZData',get(Line,'ZData')+sphere(surf.cylinder(j,1),3) + (sphere(surf.cylinder(j,2),3) -sphere(surf.cylinder(j,1),3))/norm(sphere(surf.cylinder(j,2),1:3)-sphere(surf.cylinder(j,1),1:3))*sphere(surf.cylinder(j,1),5));
        set(Line,'EdgeColor','none');
        set(Line,'FaceColor',[0 0 0]);
        set(Line,'FaceAlpha',surf.cylinder(j,6));
        set(Line,'EdgeAlpha',0);
    end

    function D = distance_bin(A)
        %DISTANCE_BIN       Distance matrix
        %
        %   D = distance_bin(A);
        %
        %   The distance matrix contains lengths of shortest paths between all
        %   pairs of nodes. An entry (u,v) represents the length of shortest path
        %   from node u to node v. The average shortest path length is the
        %   characteristic path length of the network.
        %
        %   Input:      A,      binary directed/undirected connection matrix
        %
        %   Output:     D,      distance matrix
        %
        %   Notes:
        %       Lengths between disconnected nodes are set to Inf.
        %       Lengths on the main diagonal are set to 0.
        %
        %   Algorithm: Algebraic shortest paths.
        %
        %
        %   Mika Rubinov, U Cambridge
        %   Jonathan Clayden, UCL
        %   2007-2013
        
        % Modification history:
        % 2007: Original (MR)
        % 2013: Bug fix, enforce zero distance for self-connections (JC)
        
        A=double(A~=0);                 %binarize and convert to double format
        
        l=1;                            %path length
        Lpath=A;                        %matrix of paths l
        D=A;                            %distance matrix
        
        Idx=true;
        while any(Idx(:))
            l=l+1;
            Lpath=Lpath*A;
            Idx=(Lpath~=0)&(D==0);
            D(Idx)=l;
        end
        
        D(~D)=inf;                      %assign inf to disconnected nodes
        D(1:length(A)+1:end)=0;         %clear diagonal
    end

    function  [lambda,efficiency,ecc,radius,diameter] = charpath(D,diagonal_dist,infinite_dist)
        %CHARPATH       Characteristic path length, global efficiency and related statistics
        %
        %   lambda                                  = charpath(D);
        %   lambda                                  = charpath(D);
        %   [lambda,efficiency]                     = charpath(D);
        %   [lambda,efficiency,ecc,radius,diameter] = charpath(D,diagonal_dist,infinite_dist);
        %
        %   The network characteristic path length is the average shortest path
        %   length between all pairs of nodes in the network. The global efficiency
        %   is the average inverse shortest path length in the network. The nodal
        %   eccentricity is the maximal path length between a node and any other
        %   node in the network. The radius is the minimal eccentricity, and the
        %   diameter is the maximal eccentricity.
        %
        %   Input:      D,              distance matrix
        %               diagonal_dist   optional argument
        %                               include distances on the main diagonal
        %                                   (default: diagonal_dist=0)
        %               infinite_dist   optional argument
        %                               include infinite distances in calculation
        %                                   (default: infinite_dist=1)
        %
        %   Outputs:    lambda,         network characteristic path length
        %               efficiency,     network global efficiency
        %               ecc,            nodal eccentricity
        %               radius,         network radius
        %               diameter,       network diameter
        %
        %   Notes:
        %       The input distance matrix may be obtained with any of the distance
        %   functions, e.g. distance_bin, distance_wei.
        %       Characteristic path length is defined here as the mean shortest
        %   path length between all pairs of nodes, for consistency with common
        %   usage. Note that characteristic path length is also defined as the
        %   median of the mean shortest path length from each node to all other
        %   nodes.
        %       Infinitely long paths (i.e. paths between disconnected nodes) are
        %   included in computations by default. This behavior may be modified with
        %   via the infinite_dist argument.
        %
        %
        %   Olaf Sporns, Indiana University, 2002/2007/2008
        %   Mika Rubinov, U Cambridge, 2010/2015
        
        %   Modification history
        %   2002: original (OS)
        %   2010: incorporation of global efficiency (MR)
        %   2015: exclusion of diagonal weights by default (MR)
        %   2016: inclusion of infinite distances by default (MR)
        
        n = size(D,1);
        if any(any(isnan(D)))
            error('The distance matrix must not contain NaN values');
        end
        if ~exist('diagonal_dist','var') || ~diagonal_dist || isempty(diagonal_dist)
            D(1:n+1:end) = NaN;             % set diagonal distance to NaN
        end
        if  exist('infinite_dist','var') && ~infinite_dist
            D(isinf(D))  = NaN;             % ignore infinite path lengths
        end
        
        Dv = D(~isnan(D));                  % get non-NaN indices of D
        
        % Mean of entries of D(G)
        lambda     = mean(Dv);
        
        % Efficiency: mean of inverse entries of D(G)
        efficiency = mean(1./Dv);
        
        % Eccentricity for each vertex
        ecc        = nanmax(D,[],2);
        
        % Radius of graph
        radius     = min(ecc);
        
        % Diameter of graph
        diameter   = max(ecc);
    end

    function info
        caudate='The caudate nucleus is part of the basal ganglia, and plays an important role in motor memory and associative learning.';
        IFG='The inferior frontal gyrus is the home of Broca''s area, a famous structure. Damage to this area causes non-fluent aphasia, which impairs an individual''s ability to produce sentences.';
        instructions = sprintf('Our ability to see, move, and speak depends on brain connectivity. Connected brain regions form networks. Stroke can cause damage to many of these structures, causing network disconnection. \nBy clicking on a specific node in the box below, you can see how a stroke can damage network communication, causing a loss of function. Hold down the ctrl button (bottom left of keyboard) while clicking to select multiple areas at once.');

    end

end
