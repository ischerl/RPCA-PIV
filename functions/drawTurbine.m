function drawTurbine(geom)

% Draws a turbine!
% default is a 2d turbine, but geom.threeD=1 makes a 3D turbine!

% geom , a structure defining the geometry of the turbine

% Possible fields (all are optional):

%%%% Geometry Fields %%%%
% N - number of blades, default, 2
% center - center of turbine, default [0,0]
% radius - radius of turbine, default 0.5
% chord - chord of foil, default 4.06/17.2
% alpha_p - preset pitch angle, default 6 deg
% theta - current angular position of one blade, default, 0

%%%% 3D options %%%%
% threeD - default, 0. If 1, plots the turbine in 3d
% foil_span - Turbine height. Default, 1.4
% end_plates - If 1, plots the turbine end plates. Default 0

%%%% plot options %%%%
% fig - figure handle, default, current figure
% ax - axes handle, default, current axes
% rad_circ - if 1, plots turbine circumference, otherwise not. default, 1
% rad_circ_color - color of radius circle
% center_dot - if 1 (default), plots a dot at the center of the turbine
% center_dot_color - center dot color
% foil_color - color vector for foil, default, blue (empty for no fill)
% foil_line_color - default, black (empty ,[], for no line)
% foil_line_width - default 1.5
% norm - if 1 (default), plot absolute size, else normalize by value of
% 'norm'

%% Figure and axes handling

if ~exist('geom','var') || isempty(geom)
    geom=[];
end

if isfield(geom,'fig')
    figure(geom.fig)
end

if isfield(geom,'ax')
    ax = geom.ax;
else
    ax = gca;
end

if ~isfield(geom,'foil_color')
    geom.foil_color=[0    0.4470    0.7410];
elseif isempty(geom.foil_color)
    geom.foil_color='none';
end

if ~isfield(geom,'foil_line_color')
    geom.foil_line_color='k';
elseif isempty(geom.foil_line_color)
    geom.foil_line_color='none';
end

if ~isfield(geom,'foil_line_width')
    geom.foil_line_width=1.5;
end

if ~isfield(geom,'rad_circ')
    geom.rad_circ=1;
end

if ~isfield(geom,'rad_circ_color')
    geom.rad_circ_color = 'k';
end

if ~isfield(geom,'rad_circle_fill_color')
    geom.rad_circle_fill_color = 0;
end

if ~isfield(geom,'center_dot')
    geom.center_dot=1;
end

if ~isfield(geom,'center_dot_color')
    geom.center_dot_color='k';
end

if ~isfield(geom,'norm')
    geom.norm = 1;
end

if ~isfield(geom,'threeD')
    geom.threeD=0;
elseif geom.threeD==1
    if ~isfield(geom,'foil_span')
        geom.foil_span=1.4;
    end
    
    if ~isfield(geom,'end_plates')
        geom.end_plates=0;
    end
end
%% Apply default geometries as necessary

if ~isfield(geom,'center')
    geom.center=[0,0,0];
end

if ~isfield(geom,'chord')
    geom.chord=4.06/17.2;
end

if ~isfield(geom,'radius')
    geom.radius=0.5;
end
outerRad=geom.radius;
geom.radius=geom.radius-geom.chord*0.18/2;

if ~isfield(geom,'alpha_p')
    geom.alpha_p=6;
end

if ~isfield(geom,'theta')
    geom.theta=0;
end

if ~isfield(geom,'N')
    geom.N=2;
end

%% Plot stuff for 2d turbine
holdstate = ishold(ax);
hold(ax,'on');

if geom.threeD==0
    
        % plot curcumference line
    if geom.rad_circ==1
        thetas=0:0.001:2*pi;
        xc=outerRad/geom.norm.*sin(thetas)+geom.center(1)/geom.norm;
        yc=outerRad/geom.norm.*cos(thetas)+geom.center(2)/geom.norm;
        plot(ax,xc,yc,'-.','Color',geom.rad_circ_color,'LineWidth',1.5)
    end
    
    if geom.rad_circle_fill_color~=0
        thetas=0:0.001:2*pi;
        xc=outerRad/geom.norm.*sin(thetas)+geom.center(1)/geom.norm;
        yc=outerRad/geom.norm.*cos(thetas)+geom.center(2)/geom.norm;
        patch(ax,xc,yc,0,'FaceColor',geom.rad_circle_fill_color)
    end
    
    
    % Plot foils
    for i=1:geom.N
        thetaN=geom.theta+(i-1)*360/geom.N;
        [x,y]=draw_foil(geom.chord/geom.norm,geom.radius/geom.norm,thetaN,geom.alpha_p);
        x=x+geom.center(1)/geom.norm;
        y=y+geom.center(2)/geom.norm;
        %geom.foil_line_width;
        if ~isempty(geom.foil_color)
            patch(ax,x,y,0,'FaceColor',geom.foil_color,'EdgeColor',geom.foil_line_color,'LineWidth',geom.foil_line_width)
        end
    end
    

    % plot center dot
    if geom.center_dot==1
        plot(ax,geom.center(1)/geom.norm,geom.center(2)/geom.norm,'o','MarkerFaceColor',geom.center_dot_color,'MarkerEdgeColor',geom.center_dot_color)
    end
    
    % axis equal
    hold off
    
else
    
    % Plot foils
    for i=1:geom.N
        thetaN=geom.theta+(i-1)*360/geom.N;
        [x,y]=draw_foil(geom.chord/geom.norm,geom.radius/geom.norm,thetaN,geom.alpha_p);
        x=x+geom.center(1)/geom.norm;
        y=y+geom.center(2)/geom.norm;
        [X,Y,Z] = extrude([x;y],[0 0;0 0;-(geom.foil_span/2+geom.center(3)) (geom.foil_span/2+geom.center(3))]./geom.norm,0);
        h=surf(ax,X,Y,Z,'FaceColor',geom.foil_color,'LineStyle','none','AmbientStrength',0.2);
        h.AmbientStrength = 0.3;
        h.DiffuseStrength = 0.8;
        h.SpecularStrength = 0.9;
        h.SpecularExponent = 25;
        h.BackFaceLighting = 'unlit';
        % cap foil ends if no endplates
        if geom.end_plates==0
            zc=ones(size(x)).*(geom.foil_span/2+geom.center(3))./geom.norm+0.001;
            h=patch(ax,'XData',x,'YData',y,'ZData',zc,'FaceColor',geom.foil_color,'LineStyle','none');
            h.AmbientStrength = 0.3;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 25;
            h.BackFaceLighting = 'unlit';
            h=patch(ax,'XData',x,'YData',y,'ZData',-zc,'FaceColor',geom.foil_color,'LineStyle','none');
            h.AmbientStrength = 0.3;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 25;
            h.BackFaceLighting = 'unlit';
        end
        
    end
    
    % plot end caps, if necessary
    if geom.end_plates==1
        ang=0:0.002:2*pi;
        xc=geom.radius*1.05.*cos(ang)./geom.norm;
        yc=geom.radius*1.05.*sin(ang)./geom.norm;
        zc=ones(size(xc)).*(geom.foil_span/2+geom.center(3))./geom.norm+0.001;
        h=patch(ax,'XData',xc,'YData',yc,'ZData',zc,'FaceColor',geom.foil_color,'LineStyle','none');
        h.AmbientStrength = 0.3;
        h.DiffuseStrength = 0.8;
        h.SpecularStrength = 0.9;
        h.SpecularExponent = 25;
        h.BackFaceLighting = 'unlit';
        h=patch(ax,'XData',xc,'YData',yc,'ZData',-zc,'FaceColor',geom.foil_color,'LineStyle','none');
        h.AmbientStrength = 0.3;
        h.DiffuseStrength = 0.8;
        h.SpecularStrength = 0.9;
        h.SpecularExponent = 25;
        h.BackFaceLighting = 'unlit';
    end
    
    % touch up the lighting and camera veiw
    
    camproj(ax,'perspective')
    lighting gouraud
    light(ax,'Position',[-1 -5 2],'Style','local')
    light(ax,'Position',[0 0 10])
end

if holdstate
    hold(ax,'on');
else
    hold(ax,'off')
end

    function [x,y]=draw_foil(c,r,theta,alpha)
        t=c*0.18; % foil thickness for NACA0018
        theta=theta*pi/180;
        alpha=alpha*pi/180;
        
        % draw foil shape
        x=0:c/500:c;
        ytop=t/0.2*(0.2969*sqrt(x/c)-0.126*x/c-0.3516*(x/c).^2+0.2843*(x/c).^3-0.1036*(x/c).^4);
        x=[x,fliplr(x)];
        y=[ytop,fliplr(-ytop)];
        x=x-c/4; % shift so quarter chord is on origin
        
        % perform rotation
        xn=x.*cos(theta-alpha)-y.*sin(theta-alpha);
        yn=x.*sin(theta-alpha)+y.*cos(theta-alpha);
        
        % perform translation
        % calculate position in cartesian
        x=xn+r*sin(-theta);
        y=yn+r*cos(-theta);
    end


end
