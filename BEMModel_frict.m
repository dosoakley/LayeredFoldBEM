%This file contains functions necessary for Boundary Element Modeling
%(BEM) of Faults and Associated Deformation in a 2D Elastic Half-Space.
%Equations mostly come from Segall (2010) and Huang and Johnson (2016).

%This requires the BEMsurface_frict.m, EdgeDisloc_gpu.m,
%EdgeDisloc_parallel.m, and RemeshLine.m files to run.

%Sign Conventions:
%For slip: Positive is reverse and negative is normal, regardless of the
%direction in which the fault dips. If a slip is broken down into x and y
%components, then the vector (sx,sy) points in the direction of slip of the
%hanging wall. For perfectly horizontal elements, we define slip as
%positive to the right and negative to the left.

%For stress: Normal stresses are positive in tension and negative in
%compression. Shear stresses are positive if pointing updip in the hanging
%wall and negative if pointing downdip in the hanging wall.


%This version includes friction on the elements. It uses the
%BEMsurface_frict class for the surfaces.

classdef BEMModel_frict < handle
    properties
        nu %Poisson's ratio (scalar)
        mu %Shear modulus (scalar)
        rho %Density (scalar)
        g = 9.8 %Acceleration due to gravity.
        surfaces %Array of objects of the BEMsurface_frict class.
        nsurf %Number of surfaces.
        nel %Total number of fault elements across all surfaces.
        nvert %Total number of vertices across all surfaces.
        obsx %x coordinates of observation points.
        obsy %y coordinates of observation points.
        min_cut_dist = 5 %Minimum allowed distance of points on one surface from a surface that cuts it.
        min_el_len = 5 %Minimum allowed element length.
        max_el_len = 150 %Maximum allowed element length.
        max_len_ratio = 1.2 %1.1 %1.2 %Maximum allowed length ratio between adjacent elements.
        tol = 0.1 %Tolerance for convergence in the ImposeChange function.
        maxit = 1e3 %Maximum number of iterations for the ImposeChange function.
        water_depth = 0; %Thickness of water above the free surface. Used in calculating lithostatic stress.
        rho_water = 1000; %Density of water. Used if water_depth is not 0.
        use_GPU_EdgeDisloc = true; %Tells whether to use the GPU to speed up calculations of stress and displacement.
        use_GPU_slip = false; %Tells whether to use the GPU in calculating slip in ImposeChange. This can speed it up a bit but can also cause instability.
        gpu_source_block_size = 1000; %Number of sources to process at a time if use_GPU_EdgeDisloc is true. 
        xmax = 1e10; %Maximum allowed x coordinate of any point in the model. Used to catch problems and return an error.
    end
    methods
        function M = BEMModel_frict(nu,mu,rho)
            %Create a BEMModel object.
            %Adding the surfaces is done separately with AddSurface.
            [M.nu,M.mu,M.rho] = deal(nu,mu,rho);
            [M.nsurf,M.nel,M.nvert] = deal(0,0,0);
        end
        function AddSurface(M,x,y,mu,varargin)
            %Add a surface to the BEMModel.
            %x and y are coordinate vectors for the surface points.
            %mu and c can be scalars (to give the same value to each
            %element) or vectors.
            %varargin can optionally contain values for dsigma_s and
            %dsigma_n (see BEMsurface_frict.m).
            surface = BEMsurface_frict(x,y,mu,varargin{:});
            M.surfaces = [M.surfaces,surface];
            M.nsurf = M.nsurf+1;
            M.nel = M.nel+surface.nel;
            M.nvert = M.nvert+surface.nvert;
        end
        function RemoveSurface(M,ind)
            %Remove a surface from the BEMModel.
            %ind is the index of the surface in M.surfaces.
            M.nsurf = M.nsurf-1;
            M.nel = M.nel-M.surfaces(ind).nel;
            M.nvert = M.nvert-M.surfaces(ind).nvert;
            M.surfaces = M.surfaces([1:ind-1,ind+1:end]);
        end
        function CutIntersections(M,icut,varargin)
            %Find surfaces that are intersected by the surface cut_surf and
            %cut them at the intersection point. icut should be the
            %index of one surface in the model. Cutting is done by setting
            %the elements' active tag to false.
            %To avoid problems that arise when surfaces get too close to
            %each other, cutting is done with an envelope a distance
            %min_cut_dist to either side of cut_surf rather than with 
            %cut_surf itself.
            %vargin is an array of indices of surfaces to cut. If this is
            %absent, all surfaces (except icut) will be cut.
            nx = [M.surfaces(icut).nx(1),(M.surfaces(icut).nx(1:end-1)+M.surfaces(icut).nx(2:end))/2,M.surfaces(icut).nx(end)]; %Normals at the points.
            ny = [M.surfaces(icut).ny(1),(M.surfaces(icut).ny(1:end-1)+M.surfaces(icut).ny(2:end))/2,M.surfaces(icut).ny(end)];
            L = sqrt(nx.^2+ny.^2);
            [nx,ny] = deal(nx./L,ny./L);
            [env1x,env1y] = deal(M.surfaces(icut).x+M.min_cut_dist*nx,M.surfaces(icut).y+M.min_cut_dist*ny);
            [env2x,env2y] = deal(M.surfaces(icut).x-M.min_cut_dist*nx,M.surfaces(icut).y-M.min_cut_dist*ny);
            angle1 = atan2(diff(M.surfaces(icut).y(2:-1:1)),diff(M.surfaces(icut).x(2:-1:1))); %Direction first segment of icut points.
            angle2 = atan2(diff(M.surfaces(icut).y(end-1:end)),diff(M.surfaces(icut).x(end-1:end))); %Direction last segment of icut points.
            polyx = [env1x,M.surfaces(icut).x(end)+M.min_cut_dist*cos(angle2),...
                env2x((end:-1:1)),M.surfaces(icut).x(1)+M.min_cut_dist*cos(angle1),env1x(1)]; %Polygon coordinates.
            polyy = [env1y,M.surfaces(icut).y(end)+M.min_cut_dist*sin(angle2),...
                env2y((end:-1:1)),M.surfaces(icut).y(1)+M.min_cut_dist*sin(angle1),env1y(1)];
            if isempty(varargin)
                surfs_to_cut = 1:M.nsurf; %Cut all the surfaces. (icut will still be excluded by the if i~=icut line below)
            else
                surfs_to_cut = varargin{1};
            end
            for i = surfs_to_cut
                if i ~= icut %n give the original index position, which is what icut is based on.
                    [xcut,ycut,ii] = polyxpoly(M.surfaces(i).x,M.surfaces(i).y,polyx,polyy);
                    %Add new points at the intersections.
                    for k = 1:length(xcut)
                        ind = ii(k,1);
                        AddPoint(M.surfaces(i),xcut(k),ycut(k),ind)
                        [M.nel,M.nvert] = deal(M.nel+1,M.nvert+1);
                        if size(ii,1)>1
                            mask = ii(:,1)>ind | (ii(:,1)==ind & xcut>min(M.surfaces(i).x(ind+1:ind+2)) &...
                                xcut<max(M.surfaces(i).x(ind+1:ind+2))); %Mask of ii entries that need to be increased.
                            ii(mask,1) = ii(mask,1)+1;
                        end
                    end
                    %Deactivate elements inside the polygon. I could do
                    %this by checking both end points or by checking center
                    %points, but checking the end points sometimes misses
                    %ones right on the edge of the polygon, so I'm using
                    %center points. There's a very small chance of an
                    %element with its center outside a complicated polygon 
                    %and endpoints inside it, but that seems unlikely.
                    mask = inpolygon(M.surfaces(i).cx,M.surfaces(i).cy,polyx,polyy); %This gets points in the polygon or on its edges.
                    M.surfaces(i).active(mask) = false; %Make elements inside the polygon inactive.
                    n_removed = MergeInactiveElements(M.surfaces(i));
                    [M.nel,M.nvert] = deal(M.nel-n_removed,M.nvert-n_removed);
                end
            end
        end
        function AddObsPoints(M,x,y)
            %Add observation points to the BEMModel.
            %These aren't part of slipping surfaces, but their
            %displacements will be tracked.
            %If there are no existing observation points, x and y can be
            %any shape, but if you call this function more than once, then
            %x and y should be row vectors in all cases.
            if isempty(M.obsx)
                M.obsx = x;
                M.obsy = y;
            else
                M.obsx = [M.obsx,x];
                M.obsy = [M.obsy,y];
            end
        end
        function [Gs,Gn] = CalculateG(M)
            %Calculate the G matrices that relate stresses on each element
            %to slip on each element by the equation sigma=G_sigma*s as in Eqn. 1
            %of Huang and Johnson (2016) and similar equations for other streses.
            %Gs is for shear stress. Gn is for normal stress.
            %cx, cy = The x and y coordinates of each fault element centroid.
            %nx, ny = The upward-pointing unit normal vector to each fault element.
            %len = The length of each fault element.
            %Row indices of G tell which element sigma_s is acting on. Column indices
            %tell which element s (the slip) is occurring on.
            %Stresses are calculated for a unit slip, which can then
            %be multiplied by a constant to get the stress due to the unknown true
            %slip.
            x_source = zeros(2,M.nel);
            y_source = zeros(2,M.nel);
            ind_start = 1;
            for j =1:M.nsurf
                x_source(1,ind_start:ind_start+M.surfaces(j).nel-1) = M.surfaces(j).x(1:end-1); %I'm not sure if this is the right orientation. I need to check.
                x_source(2,ind_start:ind_start+M.surfaces(j).nel-1) = M.surfaces(j).x(2:end);
                y_source(1,ind_start:ind_start+M.surfaces(j).nel-1) = M.surfaces(j).y(1:end-1);
                y_source(2,ind_start:ind_start+M.surfaces(j).nel-1) = M.surfaces(j).y(2:end);
                ind_start = ind_start+M.surfaces(j).nel;
            end
            %Create arrays of observation element points.
            x_obs = [M.surfaces.cx];
            y_obs = [M.surfaces.cy];
            %Create arrays of slip values.
            sx_sign = -sign([M.surfaces.nx]);
            sx_sign(sx_sign==0) = 1; %For perfectly flat elements.
            sx_u = sx_sign.*abs([M.surfaces.ny]);
            sy_u = abs([M.surfaces.nx]);
            if M.use_GPU_EdgeDisloc
                [sigmaxx,sigmaxy,sigmayy] = EdgeDisloc_gpu.Stress(x_source,y_source,...
                    x_obs,y_obs,sx_u,sy_u,M.nu,M.mu,M.gpu_source_block_size);
            else
                [sigmaxx,sigmaxy,sigmayy] = EdgeDisloc_parallel.Stress(x_source,y_source,...
                    x_obs,y_obs,sx_u,sy_u,M.nu,M.mu);
            end
            [Gs,Gn] = ProjectStress(sigmaxx,sigmaxy,sigmayy,repmat([M.surfaces.nx]',[1,M.nel]),repmat([M.surfaces.ny]',[1,M.nel]));
            if M.use_GPU_EdgeDisloc && ~M.use_GPU_slip
                [Gs,Gn] = deal(gather(Gs),gather(Gn)); %Gs and Gn need to be converted to non-GPU arrays for ImposeChange.
            end
            if ~M.use_GPU_EdgeDisloc && M.use_GPU_slip
                [Gs,Gn] = deal(gpuArray(Gs),gpuArray(Gn)); %Gs and Gn need to be converted to GPU arrays for ImposeChange.
            end
        end
        
        function s = ImposeChange(M,sigma_xx_new,sigma_xy_new,sigma_yy_new,s_im,BCtype)
            %This function imposes a change on the model by imposing new (additional)
            %stresses onto the elements and / or imposing slip onto the elements.
            %First, any imposed slip is applied and new stresses are
            %calculated.
            %Next, the slip on each fault element necessary to reduce the
            %shear stress on all elements below the frictional resistance
            %is calculated, and new stresses are calculated. This second
            %stage is done iterativley until all shear stresses are below
            %the frictional resistance.
            %Finally, the final slip values on all elements are used to
            %move the model.
            %I'm solving for the slip until it converges before I move
            %anything, which means I only have to calculate the influence
            %coefficients once, but I might want to consider doing the
            %movement iteratively too.
            %sigma_xx_new, sigma_xy_new, and sigma_yy_new are the imposed shear
            %stress on each element. They should be vectors equal in length
            %to the number of fault elements.
            %s_im = imposed slip on each element. Should be of length equal to number of elements with BCtype = 2.
            %BCtype = boundary condition type on each element. 1 for
            %stress (frictional), 2 for prescribed slip.
            if isrow(BCtype)
                BCtype = BCtype'; %Needs to be a column vector for the matrix equation.
            end
            %Project the new stresses.
            [sigma_s_new,sigma_n_new] = ProjectStress(sigma_xx_new,...
                sigma_xy_new,sigma_yy_new,[M.surfaces.nx]',[M.surfaces.ny]');
            %Calculate the lithostatic stresses.
            %Keeping these separate makes it easy to recalculate them as points move.
            sigma_litho = -M.rho*M.g*abs([M.surfaces.cy]')-M.rho_water*M.g*M.water_depth; %Negative for compression in this sign convention.
            [sigma_s_litho,sigma_n_litho] = ProjectStress(sigma_litho,...
                zeros(M.nel,1),sigma_litho,[M.surfaces.nx]',[M.surfaces.ny]');
            %Get the old stresses.
            [dsigma_s,dsigma_n] = deal(cat(1,M.surfaces.dsigma_s),cat(1,M.surfaces.dsigma_n));
            %Combine the stresses: lithostatic stress + stress from previous increments + new stress.
            sigma_s = sigma_s_litho+dsigma_s+sigma_s_new;
            sigma_n = sigma_n_litho+dsigma_n+sigma_n_new;
            %Calculate the influence coefficients.
            [Gs,Gn] = CalculateG(M);
            %Combine all mu_fric and active together. This seems like
            %unnecessary memory overhead, but since I have to mask them
            %out, I'm not sure how else to do it.
            mu_fric_all = cat(1,M.surfaces.mu_fric);
            active_all = cat(1,M.surfaces.active);
            %Convert to gpuArray if needed:
            if M.use_GPU_slip
                [sigma_s,sigma_n,mu_fric_all,active_all] = deal(gpuArray(sigma_s),...
                    gpuArray(sigma_n),gpuArray(mu_fric_all),gpuArray(active_all));
            end
            %Initialize the slip on each element.
            s = zeros(M.nel,1,class(sigma_s)); %Accumulated slip on each element.
            %Deal with any slip boundary conditions first.
            if any(BCtype==2) %If slip BCs exist.
                s(BCtype==2) = s_im;
                sigma_s = sigma_s+Gs*s;
                sigma_n = sigma_n+Gn*s;
            end
            %Calculate slip on the elements that don't have slip boundary conditions.
            converged = false;
            i = 0; %Iteration number.
            while ~converged
                i = i+1;
                %Figure out which elements can slip.
                failure_criterion = - mu_fric_all.*sigma_n; %sigma_n is assumed to be negative in compression here.
                mask1 = (abs(sigma_s) > failure_criterion) & BCtype == 1 & sigma_n <= 0 & active_all; %Elements that fail in compression.
                mask2 = sigma_n > 0 & BCtype == 1 & abs(sigma_s) > 0 & active_all; %Elements that fail in tension.
                mask = mask1 | mask2;
                %Calculate the slip for only the elements that can slip.
                s_iter = zeros(M.nel,1,class(s)); %Slip in this iteration only.
                s_iter(mask)=(sign(sigma_s(mask)).*Gs(mask,mask)+[mu_fric_all(mask1).*Gn(mask1,mask);...
                    zeros(sum(mask2),sum(mask),class(s))])\([-mu_fric_all(mask1).*sigma_n(mask1);zeros(sum(mask2),1,class(s))]...
                    -abs(sigma_s(mask)));
                s = s+s_iter;
                %Calculate the change in stress.
                sigma_s_old = sigma_s; %Needed for the convergence test.
                sigma_s = sigma_s+Gs*s_iter;
                sigma_n = sigma_n+Gn*s_iter;
                if mean(abs(sigma_s-sigma_s_old)./(1+abs(sigma_s_old)))<M.tol || sum(mask)==0
%                 if max(abs(sigma_s-sigma_s_old)./(1+abs(sigma_s_old)))<M.tol || sum(mask)==0
                    %First part is the convergence criterion from Eqn. 2 of Cooke
                    %and Pollard (1997). They have abs(sigma_s) -
                    %abs(sigma_s_old), but abs(sigma_s-sigma_s_old) makes
                    %more sense to me.
                    %I'm not sure if I should require all to be <M.tol or
                    %just the mean. Mean leads to faster convergence, of
                    %course, and it doesn't seem to make a big difference
                    %to the final fold geometry.
                    converged = true;
                end
                if ~converged && i>=M.maxit
                    disp('Warning: Slip calculation did not converge.')
                    converged = true;
                end
            end
            clear mu_fric_all active_all mask mask1 mask2 mu_fric_all ...
                failure_criterion sigma_s_old
            %Calculate the displacements of all points in the model.
            %Slip:
            sx_sign = -sign([M.surfaces.nx]);
            sx_sign(sx_sign==0) = 1; %For perfectly flat elements.
            sx = s'.*sx_sign.*abs([M.surfaces.ny]);%Slip vector is perpendicular to normal vector.
            sy = s'.*abs([M.surfaces.nx]);
            %Source points (nel sets of 2 end points):
            [x_all,y_all] = deal([M.surfaces.x],[M.surfaces.y]);
            vert_inds = setdiff(1:M.nvert,cumsum([M.surfaces.nvert])); %Indices in x_all, y_all for first vertices for each element.
            [ex1,ey1] = deal(x_all(vert_inds),y_all(vert_inds)); %Starting points of each element.
            [ex2,ey2] = deal(x_all(vert_inds+1),y_all(vert_inds+1)); %Ending points of each element.
            [x_source,y_source] = deal(zeros(2,M.nel),zeros(2,M.nel));
            [x_source(1,:),x_source(2,:)] = deal(ex1,ex2);
            [y_source(1,:),y_source(2,:),] = deal(ey1,ey2);
            mask = s~=0; %Mask for the elements that slip
            if M.use_GPU_EdgeDisloc
                [ux_array,uy_array] = EdgeDisloc_gpu.Displacement(x_source(:,mask),...
                    y_source(:,mask),x_all,y_all,sx(mask),sy(mask),M.nu,M.gpu_source_block_size);
            else
                [ux_array,uy_array] = EdgeDisloc_parallel.Displacement(x_source(:,mask),...
                y_source(:,mask),x_all,y_all,sx(mask),sy(mask),M.nu);
            end
            ux_array(isnan(ux_array)) = 0; %NaNs occur at the tip points of the source segment.
            uy_array(isnan(uy_array)) = 0;
            [ux,uy] = deal(sum(ux_array,2)',sum(uy_array,2)'); %Sum over the sources (second dimension).
            if M.use_GPU_EdgeDisloc
                [ux,uy] = deal(gather(ux),gather(uy)); %Bring them back to the CPU.
            end
            clear ux_array uy_array
            if ~isempty(M.obsx) %If observation points exist.
                if M.use_GPU_EdgeDisloc
                    [uxobs_array,uyobs_array] = EdgeDisloc_gpu.Displacement(x_source(:,mask),...
                        y_source(:,mask),M.obsx(:)',M.obsy(:)',sx(mask),sy(mask),M.nu,M.gpu_source_block_size);
                else
                    [uxobs_array,uyobs_array] = EdgeDisloc_parallel.Displacement(x_source(:,mask),...
                        y_source(:,mask),M.obsx(:)',M.obsy(:)',sx(mask),sy(mask),M.nu);
                end
                uxobs_array(isnan(uxobs_array)) = 0; %This can happen when an observation point is at an element tip.
                uyobs_array(isnan(uyobs_array)) = 0;
                [uxobs,uyobs] = deal(sum(uxobs_array,2),sum(uyobs_array,2)); %Sum over the sources (second dimension).
                if M.use_GPU_EdgeDisloc
                    [uxobs,uyobs] = deal(gather(uxobs),gather(uyobs)); %Bring them back to the CPU.
                end
                clear uxobs_array uyobs_array
            end
            clear x_source y_source sx sy
            %Loop through and deform each surface.
            ind_start = 1; %Index of ux and uy for the first point in this surface.
            for j = 1:M.nsurf
                ind_end = ind_start+M.surfaces(j).nvert-1;
                DeformSurface(M.surfaces(j),ux(ind_start:ind_end),uy(ind_start:ind_end));
                ind_start = ind_end+1; %Starting index for the next surface.
            end
            %Move the observation points.
            if ~isempty(M.obsx)
                M.obsx(:) = M.obsx(:)+uxobs;
                M.obsy(:) = M.obsy(:)+uyobs;
            end
            %Calculate the modified xx, xy, and yy stresses, and assign
            %them to each element of each surface. Also update the
            %accumulated slip on each element.
            if M.use_GPU_slip
                s = gather(s); %Below expects s as non-GPU array.
            end
            ind_start = 1; %Index of ux and uy for the first point in this surface.
            [dsigma_s_new,dsigma_n_new] = deal(gather(sigma_s)-sigma_s_litho,gather(sigma_n)-sigma_n_litho);
            for j = 1:M.nsurf
                ind_end = ind_start+M.surfaces(j).nel-1;
                M.surfaces(j).dsigma_s = dsigma_s_new(ind_start:ind_end);
                M.surfaces(j).dsigma_n = dsigma_n_new(ind_start:ind_end);
                M.surfaces(j).slip = M.surfaces(j).slip+s(ind_start:ind_end);
                ind_start = ind_end+1; %Starting index for the next surface.
            end
            %Delete anything above the free surface.
            CheckFreeSurface(M)
            %Fix any self-intersections in the fault or beds:
            for j = 1:M.nsurf
                n_removed = FixSelfIntersections(M.surfaces(j));
                if n_removed~=0
                    [M.nel,M.nvert] = deal(M.nel-n_removed,M.nvert-n_removed);
                end
            end
            %Remesh, but first check that the model hasn't blown up huge.
            %If it has, return an error since something is wrong, and remeshing
            %can lead to a huge number of elements, causing the program to 
            %freeze up in the remeshing loop.
            if max(abs([M.surfaces.x]))>M.xmax
                error('Model exceeds xmax.')
            else
                %Remesh the lines.
                n_change = 0;
                for j = 1:M.nsurf
                    n_change_surf = Remesh(M.surfaces(j),M.min_el_len,M.max_el_len,M.max_len_ratio,M.maxit);
                    n_change = n_change+n_change_surf;
                end
                if n_change~=0
                    [M.nel,M.nvert] = deal(M.nel+n_change,M.nvert+n_change);
                end
            end
        end
        function s = ImposeSlip(M,s,inds)
            %Impose slip on fault elements and deform the model in response.
            %s should be a vector equal in length to the number of elements in inds.
            %inds = the indices of the elements on which to apply slip.
            BCtype = ones(M.nel,1);
            BCtype(inds) = 2;
            s = ImposeChange(M,zeros(M.nel,1),zeros(M.nel,1),zeros(M.nel,1),s,BCtype);
        end
        function s = ImposeStress_ff(M,ffsigma_xx,ffsigma_xy,ffsigma_yy)
            %Impose a far-field stress on the model (defined by its xx, xy,
            %and yy components), which applies to all fault elements. Then
            %deform the model in response to this stress.
            s = ImposeChange(M,ffsigma_xx*ones(M.nel,1),ffsigma_xy*ones(M.nel,1),ffsigma_yy*ones(M.nel,1),[],ones(M.nel,1));
        end
        function s = ImposeHorizStrain_ff(M,ffepsilon_xx,yref)
            %Impose a far-field horizontal strain on the model,
            %which is converted into stress. Then deform the model in
            %response to this stress. yref is the line of no vertical
            %strain.
            lambda = 2*M.mu*M.nu/(1-2*M.nu);
            ffepsilon_yy = -ffepsilon_xx*M.nu/(1-M.nu); %Huang and Johnson, Eqn. 2
            ffsigma_xx = 2*M.mu*ffepsilon_xx+lambda*(ffepsilon_xx+ffepsilon_yy); %Far-field horizontal stress. Huang and Johnson, Eqn. 3.
            s = ImposeStress_ff(M,ffsigma_xx,0,0);
            %Impose the far-field strain as in Huang and Johnson, Eqn. 6 and 7.
            for i=1:M.nsurf
                DeformSurface(M.surfaces(i),ffepsilon_xx*M.surfaces(i).x,ffepsilon_yy*(M.surfaces(i).y-yref))
            end
            if ~isempty(M.obsx)
                M.obsx = M.obsx+ffepsilon_xx*M.obsx;
                M.obsy = M.obsy+ffepsilon_yy*(M.obsy-yref);
            end
            %Strain may put points above the free surface that
            %ImposeStress_ff didn't.
            CheckFreeSurface(M)
        end
        function CheckFreeSurface(M)
            %This algorithm does not work if anything goes above the free
            %surface (y=0). Therefore, we assume that anything above the
            %free surface gets instantly eroded away.
%             %If observation points are above the free surface, set them to
%             %NaN.
%             if ~isempty(M.obsx)
%                 mask = M.obsy>0;
%                 M.obsx(mask) = NaN;
%                 M.obsy(mask) = NaN;
%             end
            for i = 1:M.nsurf
                if any(M.surfaces(i).y>0)
                    [xcut,ycut,ii] = polyxpoly(M.surfaces(i).x,M.surfaces(i).y,...
                        [min(M.surfaces(i).x)-1,max(M.surfaces(i).x)+1],[0.0,0.0]);
                    %Add new points at the surface where it cuts the beds.
                    for j = 1:length(xcut)
                        ind = ii(j,1);
                        if M.surfaces(i).active(ind) && (M.surfaces(i).x(ind)~=xcut(j) || M.surfaces(i).y(ind)~=ycut(j))
                            %If the segment was not previously inactive and the cut isn't right at a vertex,
                            %then add a new point.
                            AddPoint(M.surfaces(i),xcut(j),ycut(j),ind) %The two new segments will be numbers ind and ind+1.
                            [M.nel,M.nvert] = deal(M.nel+1,M.nvert+1);
                            ii(j+1:end,1) = ii(j+1:end,1)+1; %We've added an element so remaining indices to deal with increase.
                        end
                    end
                    %Inactivate all segments with a point above the free surface.
                    M.surfaces(i).active(M.surfaces(i).y(1:end-1)>0 | M.surfaces(i).y(2:end)>0) = false;
                    %Merge inactive elements.
                    if ~isempty(xcut)
                        n_removed = MergeInactiveElements(M.surfaces(i));
                        [M.nel,M.nvert] = deal(M.nel-n_removed,M.nvert-n_removed);
                    end
                end
            end
        end
    end
end

function [sigma_s,sigma_n] = ProjectStress(sigmaxx,sigmaxy,sigmayy,nx,ny)
%Project stresses in 2D Cartesian coordinates into shear and
%normal components (sigma_s and sigma_n) on fault elements.
%sigmaxx, sigmaxy, sigmayy = The xx, xy, and yy stress components on each fault element.
%nx, ny = The upward-pointing unit normal vector to each fault element.
sigma_s = sign(nx).*((sigmayy-sigmaxx).*nx.*ny+sigmaxy.*(nx.^2-ny.^2)); %The sign(nx) converts to the positive if pointing updip convention.
sigma_n = sigmaxx.*nx.^2+sigmayy.*ny.^2+2*sigmaxy.*nx.*ny;
end