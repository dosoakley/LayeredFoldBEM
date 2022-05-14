%This attempts to implement a 2D elastic BEM driven by remote strain, as in
%Huang and Johnson (2016).
%This version propagates according to a P/S ratio, where S is the maximum
%slip on the fault in each deformation increment.

%Define the fault geometry.
%Model domain is a half space with the free surface at y = 0.
dip = 40*pi/180;
xt = 200; %(xt, yt) is the position of the fault tip.
yt = -700;
fault_length = 2e3; %Length along dip.
seg_length = 10; %Length of each fault segment.
xb = xt-fault_length*cos(dip); %(xb, yb) is the position of the fault base.
yb = yt-fault_length*sin(dip);
PS_ratio = 2; %Fault propagation-to-slip ratio.
prop_at_dip = false; %If true, propagation is always at the dip angle above. If false, it is at the orientation of the current fault tip.

%Define the elastic constants
nu = 0.25; %Poisson's ratio
E = 3e10; %Young's modulus. Was 3e10.
lambda = (E*nu)/((1+nu)*(1-2*nu)); %Lame's first parameter.
mu = E/(2*(1+nu)); %Lame's second parameter / shear modulus
rho = 2500; %Density (kg/m^3)
mu_fric_fault = 0.2; %Coefficient of friction for fault.
mu_fric_beds = 0.6; %Coefficient of friction for beds. Was 0.6.

%Define the far-field stress field based on the far-field strains.
inc_strain = -0.002; %Strain per increment.
n_inc = 100; %Number of strain increments.
ffepsilon_xx = inc_strain;
ffepsilon_yy = -ffepsilon_xx*nu/(1-nu); %Huang and Johnson, Eqn. 2
ffsigma_xx = 2*mu*ffepsilon_xx+lambda*(ffepsilon_xx+ffepsilon_yy); %Far-field horizontal stress. Huang and Johnson, Eqn. 3.

%Define the fault geometry.
nseg = ceil(fault_length/seg_length);
b_dat = zeros(nseg,6);
xf = [xt-(0:nseg-1)*seg_length*cos(dip),xb]'; %Fault points.
yf = [yt-(0:nseg-1)*seg_length*sin(dip),yb]';

%Create the BEM model and add the fault to it.
model = BEMModel_frict(nu,mu,rho);
model.min_cut_dist = 5;
AddSurface(model,xf,yf,mu_fric_fault)

%Define the beds and add them to the model.
xmin = -2.5e3;
xmax = 2.5e3;
dx = 50;
ymin = -2000;
ymax = -50;
dy = 200; %Was 50. Then 200. Then 100.
xgrid = xmin:dx:xmax;
ygrid = ymin:dy:ymax;
[X,Y] = meshgrid(xgrid,ygrid);
% Add as active surfaces:
for i = 1:length(ygrid)
    AddSurface(model,X(i,:),Y(i,:),mu_fric_beds);
end
% AddObsPoints(model,X,Y)

%Cut the beds where they intersect the fault.
CutIntersections(model,1) %Fault is the first surface.

%Plot the initial state.
figure(1)
plot(model.surfaces(1).x,model.surfaces(1).y,'r-','LineWidth',2)
hold on
plot([xmin,xmax],[0,0],'k-','LineWidth',2)
for j = 2:model.nsurf
    mask = model.surfaces(j).active; %Active segments only.
    change_pts = [1;find(diff(mask)~=0)+1;model.surfaces(j).nvert]; %Vertex numbers separating active and inactive segments.
    for k = 1:length(change_pts)-1
        if mask(change_pts(k)) %Active
            plot(model.surfaces(j).x(change_pts(k):change_pts(k+1)),model.surfaces(j).y(change_pts(k):change_pts(k+1)),'k-')
        end
    end
end
if ~isempty(model.obsx)
   for j = 1:size(Y,1)
       plot(model.obsx(j,:),model.obsy(j,:),'k-')
   end
end
axis equal
hold off
title('Initial State')
xlabel('Distance (m)')
ylabel('Elevation (m)')

%Loop through the strain increments.
tic;
acc_prop = 0; %Accumulated fault propagation that hasn't yet been added on to it.
reached_surface = false; %Tells if the fault reached the surface, at which point it should stop propagating.
for i = 1:n_inc
    disp(i)
    %Impose the far-field horizontal strain and deform the model in response.
    slip = ImposeHorizStrain_ff(model,ffepsilon_xx,yb);
    if ~reached_surface
        %Propagate the fault:
%         disp([min(slip(1:model.surfaces(1).nel)),max(slip(1:model.surfaces(1).nel)),mean(slip(1:model.surfaces(1).nel))])
        prop_inc = PS_ratio*max(abs(slip(1:model.surfaces(1).nel)));
%         prop_inc=2;
        acc_prop = acc_prop+prop_inc;
        disp(['propagation = ',num2str(prop_inc)])
        while acc_prop>=seg_length
            [xfault,yfault]=deal(model.surfaces(1).x,model.surfaces(1).y);
            if prop_at_dip
                xnew = xfault(1)+seg_length*cos(dip);
                ynew = yfault(1)+seg_length*sin(dip);
            else
                tip_angle = atan(diff(yfault(1:2))/diff(xfault(1:2)));
                xnew = xfault(1)+seg_length*cos(tip_angle);
                ynew = yfault(1)+seg_length*sin(tip_angle);
            end
            AddPoint(model.surfaces(1),xnew,ynew,0) %Add the new point to the tip (becoming new element 1.
            model.nvert = model.nvert+1;
            model.nel = model.nel+1;
            if ynew>0
                CheckFreeSurface(model)
                reached_surface = true;
            end
            acc_prop = acc_prop-seg_length;
        end
    end
    %Check if any deformed beds are now crossing the fault, and cut them if so.
    CutIntersections(model,1) %Fault is the first surface.
    %Plot the deformed state.
    figure(2)
%     plot(model.surfaces(1).x,model.surfaces(1).y,'r-','LineWidth',2)
    plot([xmin,xmax],[0,0],'k-','LineWidth',2)
    hold on
%     plot([xmin,xmax],[0,0],'k-','LineWidth',2)
    for j = 1:model.nsurf
        mask = model.surfaces(j).active; %Active segments only.
        change_pts = [1;find(diff(mask)~=0)+1;model.surfaces(j).nvert]; %Vertex numbers separating active and inactive segments.
        for k = 1:length(change_pts)-1
            if mask(change_pts(k)) %Active
                if j == 1
                    plot(model.surfaces(j).x(change_pts(k):change_pts(k+1)),model.surfaces(j).y(change_pts(k):change_pts(k+1)),'r-','LineWidth',2)
                else
                    plot(model.surfaces(j).x(change_pts(k):change_pts(k+1)),model.surfaces(j).y(change_pts(k):change_pts(k+1)),'k-')
                end
            end
        end
    end
    if ~isempty(model.obsx)
       for j = 1:size(Y,1)
           plot(model.obsx(j,:),model.obsy(j,:),'k-')
       end
    end
    axis equal
    hold off
    xlabel('Distance (m)')
    ylabel('Elevation (m)')
end
t = toc;
disp(['Time = ',num2str(t),' seconds.'])

%Plot the total accumulated slip on each model element.
%I might want to plot element centers rather than vertices, since that's
%where the slip is measured.
figure(3)
plot([xmin,xmax],[0,0],'k-','LineWidth',2)
hold on
min_slip = min(cat(1,model.surfaces.slip));
max_slip = max(cat(1,model.surfaces.slip));
for j = 1:model.nsurf
    mask = model.surfaces(j).active; %Active segments only.
    change_pts = [1;find(diff(mask)~=0)+1;model.surfaces(j).nvert]; %Vertex numbers separating active and inactive segments.
    for k = 1:length(change_pts)-1
        if mask(change_pts(k)) %Active
            p = plot(model.surfaces(j).x(change_pts(k):change_pts(k+1)),model.surfaces(j).y(change_pts(k):change_pts(k+1)),'k-');
            drawnow
            cd = colormap('jet');
            element_slip = model.surfaces(j).slip(mask(change_pts(k):change_pts(k+1)-1));
            vertex_slip = [element_slip(1);(element_slip(1:end-1)+element_slip(2:end))/2;element_slip(end)]; %To plot, I need this by vertex rather than edge.
            cd = interp1(linspace(min_slip,max_slip,length(cd)),cd,vertex_slip);
            cd = uint8(cd'*255);
            cd(4,:) = 255;
            set(p.Edge,'ColorBinding','interpolated','ColorData',cd)
        end
    end
end
caxis([min_slip,max_slip])
c = colorbar();
c.Label.String = 'Slip (m)';
hold off
title('Total Accumulated Slip')
xlabel('Distance (m)')
ylabel('Elevation (m)')