function [x,y,active,change] = RemeshLine(x,y,Lmin,Lmax,rmax,maxit,active)
%Remesh a line in 2D to a specified desired element length.
%The goal is to maintain the shape as much as possible while keeping
%element sizes within a specified range and avoiding abrupt changes in
%size of adjacent elements.
%Arguments:
    %x, y = The point coordinates of the initial line.
    %Lmin = The minimum allowed element length.
    %Lmax = The maximum allowed element length. Should be >2*Lmin to prevent an infinite loop.
    %rmax = The maximum allowed size ratio of adjacent elements.
    %maxit = The mamximum number of iterations.
    %active = Tells whether a segment can be modified in the remeshing
    %   (true) or not (false).
%Returns:
    %x, y = The point coordinates of the final line.
    %active = The Boolean array for the final line telling which segments 
        %are allowed to be remeshed.
    %change = Boolean to tell whether or not anything has changed.
    
%I might want to add a de-spiking algorithm too.
%Also, I might want to put the check for self intersections into this
%function.

converged = [false,false,false]; %Tells whether any of the 3 tests have come back empty.
n = 0; %Number of iterations.
while any(~converged)
    %Split large elements.
    len = sqrt(diff(x).^2+diff(y).^2);
    inds_large = find(len>Lmax & active'); %These are element indices. They are bounded by vertices inds and inds+1.
    converged(1) = isempty(inds_large);
    for i = 1:length(inds_large)
        ind = inds_large(i);
        x = [x(1:ind),mean(x(ind:ind+1)),x(ind+1:end)];
        y = [y(1:ind),mean(y(ind:ind+1)),y(ind+1:end)];
        active = [active(1:ind);true;active(ind+1:end)];
        inds_large = inds_large+1; %Because we've added an element, the remaining indices increase.
    end
    %Remove small elements.
    len = sqrt(diff(x).^2+diff(y).^2);
    if length(len)==1 %If there's only one element, don't remove it.
        inds_small = [];
    else
        inds_small = find(len<Lmin & active'); %These are element indices. They are bounded by vertices inds and inds+1.
    end
    converged(2) = isempty(inds_small);
    for i = 1:length(inds_small)
        ind = inds_small(i);
        if len(ind)<Lmin && length(len)>1 %We have to check again, because removing a previous one could have lengthened this one enough to make it > min_len.
            x = [x(1:ind-1),mean(x(ind:ind+1)),x(ind+2:end)];
            y = [y(1:ind-1),mean(y(ind:ind+1)),y(ind+2:end)];
            len = sqrt(diff(x).^2+diff(y).^2); %Update the length vector.
            active = [active(1:ind-1);active(ind+1:end)];
            inds_small = inds_small-1; %Because we've removed an element, the remaining indices decrease.
        end
    end
    %Shift points connecting adjacent elements with large length ratios so the element lengths are equal.
    %See 25/2/27 notes for derivation of these equations.
    ratio = len(2:end)./len(1:end-1); %Ratio of element length to previous.
    both_active = active(1:end-1)' & active(2:end)';
    inds_ratio = find((ratio>rmax | 1./ratio>rmax) & both_active); %These give the index of the 1st of the 2 elements involved.
    converged(3) = isempty(inds_ratio);
    for i = 1:length(inds_ratio)
        ind = inds_ratio(i);
        if (ratio(ind)>rmax || 1/ratio(ind)>rmax) %We have to check again because a previous move could have changed this.
            [L1,L2] = deal(len(ind),len(ind+1));
            X21 = [x(ind)-x(ind+1),y(ind)-y(ind+1)]; %Vector from middle point to previous point.
            X23 = [x(ind+2)-x(ind+1),(y(ind+2)-y(ind+1))];%Vector from middle point to next point.
            uhat = (X21+X23)/sqrt(sum((X21+X23).^2)); %Unit vector for displacement of middle points.
            u_mag = abs((L1^2-L2^2)/(2*dot(X21-X23,uhat))); %Magnitude of the displacement.
            [x(ind+1),y(ind+1)] = deal(x(ind+1)+u_mag*uhat(1),y(ind+1)+u_mag*uhat(2)); %ind+1 is the middle point.
            len = sqrt(diff(x).^2+diff(y).^2); %Update the length vector.
            ratio = len(2:end)./len(1:end-1); %Update the ratio vector.
        end
    end
    n = n+1;
%     disp([n,length(inds_large),length(inds_small),length(inds_ratio)])
    if n>=maxit && any(~converged)
        disp('Warning: Remeshing failed to converge.')
        break
    end
end
if n==1 %First iteration was already converged, so nothing changed.
    change = false;
else
    change = true;
end