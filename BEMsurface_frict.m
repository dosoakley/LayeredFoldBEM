
%This file defines a 2D surface class (really a line in 2D) for use in
%Boundary Element Modeling. A surface may be a fault or a bed along which
%slip can occur (which thus acts like a fault).
%A surface is defined by a set of points, which are assumed to be connected
%in the order they are specified. The line segments between the points are
%the surface elements, which are defined by their centroid and unit normal
%vector.

%This version includes friction, defined by the coefficient of friction. 
%This is defined for each element for maximum flexibility.
%If not set, it defaults to 0 as in the frictionless case.

classdef BEMsurface_frict < matlab.mixin.Copyable
    properties
        nel %number of elements
        nvert %number of vertices
        x %x coordinates of points, row vector
        y %y coordinates of points, row vector
        cx %x coordinates of element centroids, row vector
        cy %y coordinates of element centroids, row vector
        nx %x coordinates of unit normal vectors, row vector
        ny %y coordinates of unit normal vectors, row vector
        mu_fric; %Coefficient of friction, column vector
        dsigma_s %Shear stress on each element, difference from lithostatic, column vector
        dsigma_n %Normal stress on each element, difference from lithostatic, column vector
        active %Boolean. If false, the element will be ignored in calculations, column vector
        slip %Cumulative slip on the element.
        remesh %If true, the element can be modified during remeshing. If false, it cannot.
    end
    methods
        function S = BEMsurface_frict(x,y,mu,varargin)
            %Create a BEMsurface object.
            %mu (coefficient of friction) can be a scalar or a vector.
            %vargin should either be empty or be three vectors (one each
            %for dsigma_xx, dsigma_xy, and dsigma_yy).
            if ~isrow(x) %For using the BEMModel class, x and y need to be row vectors.
                x = x';
            end
            if ~isrow(y)
                y = y';
            end
            S.x = x;
            S.y = y;
            S.nvert = length(x);
            S.nel = S.nvert-1;
            S.active = true(S.nel,1);
            S.remesh = true(S.nel,1);
            S.slip = zeros(S.nel,1);
            if length(mu)>1
                S.mu_fric = mu;
            else
                S.mu_fric = mu*ones(S.nel,1);
            end
            if isempty(varargin)
                S.dsigma_s = zeros(S.nel,1);
                S.dsigma_n = zeros(S.nel,1);
            else
                S.dsigma_s = varargin{1};
                S.dsigma_n = varargin{2};
            end
            UpdateSurface(S);
        end
        function AddPoint(S,x,y,ind)
            %Add a point to the surface, splitting one element into two.
            %(x,y) are the coordinates of the point, and ind is the element
            %it should be added into. (That is, it will be added between 
            %points ind and ind+1, cutting element ind.) To add a new point
            %to the end use ind = 0 or ind = S.nel+1;
            [S.x,S.y] = deal([S.x(1:ind),x,S.x(ind+1:end)],[S.y(1:ind),y,S.y(ind+1:end)]); %This line works regardless.
            if (ind>=1 && ind<=S.nel) %Adding a point within the existing line.
                S.mu_fric = [S.mu_fric(1:ind);S.mu_fric(ind);S.mu_fric(ind+1:end)];
                S.dsigma_s = [S.dsigma_s(1:ind);S.dsigma_s(ind);S.dsigma_s(ind+1:end)];
                S.dsigma_n = [S.dsigma_n(1:ind);S.dsigma_n(ind);S.dsigma_n(ind+1:end)];
                S.active = [S.active(1:ind);S.active(ind);S.active(ind+1:end)];
                S.remesh = [S.remesh(1:ind);S.remesh(ind);S.remesh(ind+1:end)];
                S.slip = [S.slip(1:ind);S.slip(ind);S.slip(ind+1:end)];
            else %Adding a point to the end of the line.
                if ind==0
                    S.mu_fric = [S.mu_fric(1);S.mu_fric];
                    S.dsigma_s = [S.dsigma_s(1);S.dsigma_s];
                    S.dsigma_n = [S.dsigma_n(1);S.dsigma_n];
                    S.active = [true;S.active];
                    S.remesh = [true;S.remesh];
                    S.slip = [0;S.slip];
                elseif ind==S.nel+1
                    S.mu_fric = [S.mu_fric;S.mu_fric(end)];
                    S.dsigma_s = [S.dsigma_s;S.dsigma_s(end)];
                    S.dsigma_n = [S.dsigma_n;S.dsigma_n(end)];
                    S.active = [S.active;true];
                    S.remesh = [S.remesh;true];
                    S.slip = [S.slip;0];
                else
                    disp('Error: ind must be in the range 0 to S.nel+1')
                end
            end
            S.nvert = S.nvert+1;
            S.nel = S.nel+1;
            UpdateSurface(S);
        end
        function n = MergeInactiveElements(S)
            %Merge adjacent inactive elements together to reduce the total
            %number of elements.
            mask = (~S.active(1:end-1)) & (~S.active(2:end)); %Non-end points with both adjacent elements inactive.
            n = sum(mask);
            if n>0
                S.nvert = S.nvert-n;
                S.nel = S.nel-n;
                [S.x,S.y] = deal(S.x([true;~mask;true]'),S.y([true;~mask;true]'));
                %Mu and sigma values don't really matter since the elements
                %are inactive, so we're just taking the first one. Slip
                %could matter more, but taking the first is still probably
                %as good as anything.
                S.mu_fric = S.mu_fric([true;~mask]);
                [S.dsigma_s,S.dsigma_n] = deal(S.dsigma_s([true;~mask]),...
                    S.dsigma_n([true;~mask]));
                S.active = S.active([true;~mask]);
                S.remesh = S.remesh([true;~mask]);
                S.slip = S.slip([true;~mask]);
                UpdateSurface(S);
            end
        end
        function n = FixSelfIntersections(S)
            nel_old = S.nel;
            while 1
                [xi,yi,ii] = polyxpoly(S.x,S.y,S.x,S.y);
                inds = find(ii(:,1)~=ii(:,2));
                if ~isempty(inds)
                    [ii1,ii2] = deal(min(ii(inds(1),:)),max(ii(inds(1),:)));
                    %Add points to split elements at the intersection.
                    AddPoint(S,xi(inds(1)),yi(inds(1)),ii2) %Do the max index first so the other index number doesn't change.
                    AddPoint(S,xi(inds(1)),yi(inds(1)),ii1)
                    %Delete the elements and points inbetween, including one of
                    %the now duplicate points that were added.
                    vert_ind_keep = [1:ii1+1,ii2+3:S.nvert]; %+1 for vertex after element, +3 for that, point added, and don't keep duplicate.
                    el_ind_keep = [1:ii1,ii2+2:S.nel]';
                    [S.x,S.y] = deal(S.x(vert_ind_keep),S.y(vert_ind_keep));
                    S.mu_fric = S.mu_fric(el_ind_keep);
                    [S.dsigma_s,S.dsigma_n] = deal(S.dsigma_s(el_ind_keep),...
                        S.dsigma_n(el_ind_keep));
                    S.active = S.active(el_ind_keep);
                    S.remesh = S.remesh(el_ind_keep);
                    S.slip = S.slip(el_ind_keep);
                    %Update things.
                    [S.nvert,S.nel] = deal(length(vert_ind_keep),length(el_ind_keep));
                    UpdateSurface(S);
                end
                %If there are no more self-intersections break. If there
                %are, rerun polyxpoly since all the index numbers have
                %changed.
                if length(inds)<=1
                    break
                end
            end
            n = nel_old-S.nel; %Number of elements (and vertices) removed.
        end
        function n = Remesh(S,Lmin,Lmax,rmax,maxit)
            %This function remeshes a surface using the RemeshLine
            %function.
            n = 0; %Number of added (if positive) or removed (if negative) segments and vertices.
            change_pts = [1;find(diff(S.active)~=0)+1;S.nvert]; %Vertex numbers separating active and inactive segments.
            for k = 1:length(change_pts)-1
                if S.active(change_pts(k)) %Active
                    inds = change_pts(k):change_pts(k+1);
                    [cxold,cyold] = deal(S.cx(inds(1:end-1)),S.cy(inds(1:end-1)));
                    [xnew,ynew,remesh_out,change] = RemeshLine(S.x(inds),S.y(inds),Lmin,Lmax,rmax,maxit,S.remesh(inds(1:end-1)));
                    if change
                        [cxnew,cynew] = deal((xnew(1:end-1)+xnew(2:end))/2,(ynew(1:end-1)+ynew(2:end))/2);
                        S.x = [S.x(1:inds(1)-1),xnew,S.x(inds(end)+1:end)];
                        S.y = [S.y(1:inds(1)-1),ynew,S.y(inds(end)+1:end)];
                        field_names = {'mu_fric','dsigma_s','dsigma_n','slip'}; %The fields to reinterpolate.
                        %Note: Reinterpolating dsigma_s and disima_n
                        %assumes the change in orientation is small. If
                        %that's not true, it would be better to work in
                        %dsigma_xx, dsigma_xy, dsigma_yy.
                        for i = 1:length(field_names)
                            if length(cxold)>2 %scatteredInterpolant needs at least 3 values..
                                if min(cyold)==max(cyold) %This happens for perfectly flat beds, which don't work well with scatteredInterpolant.
                                    new_values = interp1(cxold',S.(field_names{i})(inds(1:end-1)),cxnew');
                                else
                                    F = scatteredInterpolant(cxold',cyold',S.(field_names{i})(inds(1:end-1)));
                                    new_values = F(cxnew',cynew');
                                end
                            else
                                new_values = mean(S.(field_names{i})(inds(1:end-1)))*ones(length(cxnew),1); %This isn't ideal but will do for now.
                            end
                            S.(field_names{i}) = [S.(field_names{i})(1:inds(1)-1);new_values;S.(field_names{i})(inds(end):end)];
                        end
                        nvert_new = length(S.x);
                        nchange = (nvert_new-S.nvert);
                        S.active = [S.active(1:inds(1)-1);true(length(cxnew),1);S.active(inds(end):end)];
                        S.remesh = [S.remesh(1:inds(1)-1);remesh_out;S.remesh(inds(end):end)];
                        n = n+nchange;
                        [S.nvert,S.nel] = deal(nvert_new,nvert_new-1);
                        change_pts(k+1:end) = change_pts(k+1:end)+nchange; %Update the vertex numbers.
                        UpdateSurface(S)
                    end
                end
            end
        end
        function UpdateSurface(S)
            %Update other things after changing x and y.
            %Calculate the centroid:
            S.cx = (S.x(1:end-1)+S.x(2:end))/2;
            S.cy = (S.y(1:end-1)+S.y(2:end))/2;
            %Calculate the normal vector:
            dx = S.x(1:end-1)-S.x(2:end);
            dy = S.y(1:end-1)-S.y(2:end);
            dxsign = sign(dx); %Using this ensures that the normal vector will always be the upward-pointing one.
            %Note: If an element is perfectly vertical, dxsign will be 0,
            %which will cause problems. That isn't likely, however.
            len = sqrt(dx.^2+dy.^2);
            S.nx = -dxsign.*dy./len;
            S.ny = dxsign.*dx./len;
        end
        function DeformSurface(S,ux,uy)
            %Apply displacements ux and uy to the elements of the surface.
            S.x = S.x+ux;
            S.y = S.y+uy;
            UpdateSurface(S);
        end
    end
end