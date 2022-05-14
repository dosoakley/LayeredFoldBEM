%This file contains functions for modeling of edge dislocations in a 2D
%Elastic Half-Space.
%Equations come from Segall (2010), with a correction from point 8 of the
%errata to that book.

%The Stress and Displacement functions are for an edge dislocation that is
%infinitely long in one direction.
%The Stress2 and Displacement2 functions are for slip over a finite
%interval, which is represented by edge dislocations of opposite sign at
%either end of the invertval.

%A note on sign convention:
%For s1, and s2, the sign convention should be such that the vector (s1,s2)
%points in the direction of slip on the hanging wall of the fault.
%The equations that I have from Segall (2010) use a sign convention for
%s1 and s2 that is the same as this for faults dipping to the right but
%opposite to it for faults dipping to the left. Therefore, in the Stress
%and Displacement functions, I have added code to flip the signs in this
%case before evaluating the equations from Segall (2010).
%For stresses calculated with the Stress function, the sign convention is
%that tensile stresses are positive and compressional stresses are negative.

classdef EdgeDisloc_gpu
    methods (Static = true)
        function [sigma11,sigma12,sigma22] = Stress(xi1_vec,xi2_vec,x1_vec,x2_vec,s1_vec,s2_vec,nu,mu,block)
            %This function models stress due to slip over a confined
            %interval, represented by edge dislocations at each end with
            %opposite signs. %It uses Equantions 3.77, 3.78, and 3.79 of
            %Segall (2010) as corrected in the Errata to that book.
            %Incides 1 and 2 correspond to the x and y directions respectively.
            %xi1_vec and xi2_vec are 2xN arrays of source segment endpoints.
            %x1_vec and x2_vec are 1xM arrays of observation points. s1_vec
            %and s2_vec are 1xM arrays of slip on each source segment.
            %nu is Poisson's ratio. mu is shear modulus.
            %block is the number of sources to work with in a block at one
            %time.
            %The output variables are M by N arrays of the stresses at each
            %of the M observation points produced by slip on each of the N
            %sources.
            
            %Figure out which point is the top (or right if horizontal) one, and which is the bottom and switch if needed.
            top_first = xi2_vec(1,:)>xi2_vec(2,:) | (xi2_vec(1,:)==xi2_vec(2,:) & xi1_vec(1,:)>xi1_vec(2,:));
            if any(~top_first) %Switch so the top (or right) entries are always first.
                [temp1,temp2] = deal(xi1_vec(1,~top_first),xi2_vec(1,~top_first));
                [xi1_vec(1,~top_first),xi2_vec(1,~top_first)] = deal(xi1_vec(2,~top_first),xi2_vec(2,~top_first));
                [xi1_vec(2,~top_first),xi2_vec(2,~top_first)] = deal(temp1,temp2);
                clear temp1 temp2
            end
            %Prepare the sigma arrays:
            N = size(xi1_vec,2); %Number of sources.
            M = length(x1_vec); %Number of observation points.
            sigma11 = zeros(M,N,'gpuArray');
            sigma12 = zeros(M,N,'gpuArray');
            sigma22 = zeros(M,N,'gpuArray');
            %Loop through the blocks.
            nblocks = ceil(N/block);
            for n = 1:nblocks
                ind = (n-1)*block+1; %Starting index.
                nsources = min(block,N-ind+1); %Number of sources in this block.
                %Set up the arrays of each source with each observation point.
                x1_unshifted = repmat(gpuArray(x1_vec'),[1,nsources]);
                x2 = repmat(gpuArray(x2_vec'),[1,nsources]);
                s1 = repmat(gpuArray(s1_vec(ind:ind+nsources-1)),[M,1]);
                s2 = repmat(gpuArray(s2_vec(ind:ind+nsources-1)),[M,1]);
                %Calculate the stresses as the sum of two edge dislocations, which cancel except along the element:
                for i = 1:2 %First edge dislocation at the tops of each element, then at the bottoms with opposite slip.
                    if i == 1
                        xi1 = repmat(gpuArray(xi1_vec(1,ind:ind+nsources-1)),[M,1]);
                        xi2 = repmat(gpuArray(xi2_vec(1,ind:ind+nsources-1)),[M,1]);
                        %Flipping the signs on s if the fault dips left only has to be done once.
                        left_mask = sign(s1)==sign(s2) | s2==0; %s2==0 is horizontal. Since the right point is considered "top", these all dip left.
                        [s1(left_mask),s2(left_mask)] = deal(-s1(left_mask),-s2(left_mask)); %Flip to the sign convention that Segall's equations use.
                        clear left_mask
                    else
                        xi1 = repmat(gpuArray(xi1_vec(2,ind:ind+nsources-1)),[M,1]);
                        xi2 = repmat(gpuArray(xi2_vec(2,ind:ind+nsources-1)),[M,1]);
                        [s1,s2] = deal(-s1,-s2);
                    end
                    %Solve the Edge Dislocation stress equations.
                    %For gpu it's better to solve the equations all in one step
                    %to use less memory.
                    x1 = x1_unshifted-xi1;%These equations are written for xi1 = 0. If that's not the case, this makes it equivalent to that case.
                    clear xi1
                    %Solve using Segall (2010) Eqns. 3.77, 3.78, and 3.79.
                    %The equations have been simplified (using Matlab's
                    %Symbolic Math Toolbox) to remove the intermediate r1
                    %and r2 variables and reduce GPU memory usage.
                    r1_sq = x1.^2+(x2-xi2).^2;
                    r2_sq = x1.^2+(x2+xi2).^2; %Eqn. 3.75 of Segall (2010).
                    sigma11(:,ind:ind+nsources-1) = sigma11(:,ind:ind+nsources-1)+...
                        ((mu*s2/(2*pi*(1-nu))).*(x1.*((x2-xi2).^2-x1.^2)./(r1_sq.^2)-...
                        x1.*((x2+xi2).^2-x1.^2)./r2_sq.^2+(4*xi2.*x1./(r2_sq.^3)).*((2*xi2-x2)...
                        .*(x2+xi2).^2+(3*x2+2*xi2).*x1.^2))+(mu*s1/(2*pi*(1-nu))).*((x2-xi2)...
                        .*((x2-xi2).^2+3*x1.^2)./r1_sq.^2-(x2+xi2).*((x2+xi2).^2+3*x1.^2)./r2_sq.^2+...
                        (2*xi2./r2_sq.^3).*(6*x2.*(x2+xi2).*x1.^2-(x2-xi2).*(x2+xi2).^3-x1.^4))); %Eqn. 3.77
                    sigma22(:,ind:ind+nsources-1) = sigma22(:,ind:ind+nsources-1)+...
                        ((-mu*s2/(2*pi*(1-nu))).*(x1.*(3*(x2-xi2).^2+x1.^2)./r1_sq.^2-...
                        x1.*(3*(x2+xi2).^2+x1.^2)./r2_sq.^2-(4*xi2.*x2.*x1./r2_sq.^3).*(3*(x2+xi2).^2-x1.^2))...
                        +(mu*s1/(2*pi*(1-nu))).*((x2-xi2).*((x2-xi2).^2-x1.^2)./r1_sq.^2-...
                        (x2+xi2).*((x2+xi2).^2-x1.^2)./r2_sq.^2-(2*xi2./r2_sq.^3).*(6*x2.*(x2+xi2).*x1.^2-...
                        (3*x2+xi2).*(x2+xi2).^3+x1.^4))); %Eqn. 3.78
                    sigma12(:,ind:ind+nsources-1) = sigma12(:,ind:ind+nsources-1)+...
                        ((mu*s2./(2*pi*(1-nu))).*((x2-xi2).*((x2-xi2).^2-x1.^2)./r1_sq.^2-...
                        (x2+xi2).*((x2+xi2).^2-x1.^2)./r2_sq.^2+(2*xi2./r2_sq.^3).*(6*x2.*(x2+xi2).*x1.^2-...
                        x1.^4+(xi2-x2).*(x2+xi2).^3))+(mu*s1./(2*pi*(1-nu))).*(x1.*((x2-xi2).^2-x1.^2)./r1_sq.^2-...
                        x1.*((x2+xi2).^2-x1.^2)./r2_sq.^2+(4*xi2.*x2.*x1./r2_sq.^3).*(3*(x2+xi2).^2-x1.^2))); %Eqn. 3.79
                end
            end
        end
        
        function [u1,u2] = Displacement(xi1_vec,xi2_vec,x1_vec,x2_vec,s1_vec,s2_vec,nu,block)
            %This function models displacement due to slip over a confined
            %interval, represented by edge dislocations at each end with
            %opposite signs. %It uses Equantion 3.74 of
            %Segall (2010) as corrected in the Errata to that book.
            %Incides 1 and 2 correspond to the x and y directions respectively.
            %xi1_vec and xi2_vec are 2xN arrays of source segment endpoints.
            %x1_vec and x2_vec are 1xM arrays of observation points. s1_vec
            %and s2_vec are 1xM arrays of slip on each source segment.
            %nu is Poisson's ratio.
            %block is the number of sources to work with in a block at one
            %time.
            %The output variables are M by N arrays of the displacements at 
            %each of the M observation points produced by slip on each of 
            %the N sources.
            
            %Figure out which point is the top (or right if horizontal) one, and which is the bottom and switch if needed.
            top_first = xi2_vec(1,:)>xi2_vec(2,:) | (xi2_vec(1,:)==xi2_vec(2,:) & xi1_vec(1,:)>xi1_vec(2,:));
            if any(~top_first) %Switch so the top (or right) entries are always first.
                [temp1,temp2] = deal(xi1_vec(1,~top_first),xi2_vec(1,~top_first));
                [xi1_vec(1,~top_first),xi2_vec(1,~top_first)] = deal(xi1_vec(2,~top_first),xi2_vec(2,~top_first));
                [xi1_vec(2,~top_first),xi2_vec(2,~top_first)] = deal(temp1,temp2);
                clear temp1 temp2
            end
            %Prepare the u arrays:
            N = size(xi1_vec,2); %Number of sources.
            M = length(x1_vec); %Number of observation points.
            u1 = zeros(M,N,'gpuArray');
            u2 = zeros(M,N,'gpuArray');
            slip0_mask = false(M,N,'gpuArray');
            %Loop through the blocks.
            nblocks = ceil(N/block);
            for n = 1:nblocks %I should probably use parfor if not gpu and for otherwise, but not sure how to do that.
                ind = (n-1)*block+1; %Starting index.
                nsources = min(block,N-ind+1); %Number of sources in this block.
                %Set up the arrays of each source with each observation point.
                x1_unshifted = repmat(gpuArray(x1_vec'),[1,nsources]);
                x2 = repmat(gpuArray(x2_vec'),[1,nsources]);
                s1 = repmat(gpuArray(s1_vec(ind:ind+nsources-1)),[M,1]);
                s2 = repmat(gpuArray(s2_vec(ind:ind+nsources-1)),[M,1]);
                for i = 1:2 %First edge dislocation at the tops of each element, then at the bottoms with opposite slip.
                    if i == 1
                        xi1 = repmat(gpuArray(xi1_vec(1,ind:ind+nsources-1)),[M,1]);
                        xi2 = repmat(gpuArray(xi2_vec(1,ind:ind+nsources-1)),[M,1]);
                        %Flipping the signs on s if the fault dips left only has to be done once.
                        left_mask = sign(s1)==sign(s2) | s2==0; %s2==0 is horizontal. Since the right point is considered "top", these all dip left.
                        [s1(left_mask),s2(left_mask)] = deal(-s1(left_mask),-s2(left_mask)); %Flip to the sign convention that Segall's equations use.
                        clear left_mask
                    else
                        xi1 = repmat(gpuArray(xi1_vec(2,ind:ind+nsources-1)),[M,1]);
                        xi2 = repmat(gpuArray(xi2_vec(2,ind:ind+nsources-1)),[M,1]);
                        [s1,s2] = deal(-s1,-s2);
                    end
                    %Solve the Edge Dislocation stress equations.
                    x1 = x1_unshifted-xi1;%These equations are written for xi1 = 0. If that's not the case, this makes it equivalent to that case.
                    clear xi1
                    theta1 = atan2(x1,(x2-xi2)); %Calculate the angles relative to the vertical axis.
                    theta2 = atan2(x1,(x2+xi2));
                    dip = atan(s2./s1);
                    dip(dip<0) = dip(dip<0)+pi; %This puts dip in the range [0,pi], which is how I've written things below to expect it. Actually, this may not be necessary.
                    %The following puts the atan branch cuts along the fault by rotating theta1
                    %to point down the fault and theta2 to point up (above the half space)
                    %along it. This means that theta2 is shifted by -pi from theta1. This
                    %matches Fig. 3.16 in which theta1 and theta2 are measured in opposite
                    %direction (clockwise vs. counterclockwise) from the x2 axis, but it
                    %doesn't seem to match Eqn.  3.76, in which theta2 appears to be measured
                    %counterclockwise. In any case, this way works, when tested for surface
                    %points against Eqn. 3.70 and 3.73, and trying to correct it afterwords for
                    %that -pi factor doesn't work.
                    theta1 = theta1+pi/2+dip; %Shift theta1 to be measured up from the fault.
                    theta2 = theta2+dip-pi/2; %Shift theta2 to be measured from a line pointing up opposite the fault from the image dislocation.
                    clear dip
                    theta1 = mod(theta1,2*pi);
                    theta2 = mod(theta2,2*pi);
                    %Make a correction for rounding errors that can occur when
                    %theta1 or theta2 is very close to 0 or 2*pi.
                    if any(2*pi-theta1<1e-10)
                        theta1(2*pi-theta1<1e-10) = 0;
                    end
                    if any(2*pi-theta2<1e-10)
                        theta2(2*pi-theta2<1e-10) = 0;
                    end
                    theta_diff = theta1-theta2;
                    clear theta1 theta2
                    %For gpu it's better to solve the equations all in one step
                    %to use less memory.
                    %Solve using Segall (2010) Eqn. 3.74.
                    %The equations have been simplified (using Matlab's
                    %Symbolic Math Toolbox) to remove the intermediate r1
                    %and r2 variables and reduce GPU memory usage.
                    r1_sq = x1.^2+(x2-xi2).^2;
                    r2_sq = x1.^2+(x2+xi2).^2; %Eqn. 3.75 of Segall (2010).
                    log_r2_r1 = log(sqrt(r2_sq)./sqrt(r1_sq));
                    %Solve Eqn. 3.74 (Segall, 2010) for the displacements:
                    u1(:,ind:ind+nsources-1) = u1(:,ind:ind+nsources-1)+...
                        ((-s1/(pi*(1-nu))).*(((1-nu)/2).*(-theta_diff)+x1.*(x2-xi2)./(4*r1_sq)-...
                        x1.*(x2+(3-4*nu).*xi2)./(4*r2_sq)+xi2.*x2.*x1.*(x2+xi2)./r2_sq.^2)...
                        +(s2/(pi*(1-nu))).*(((1-2*nu)/4).*log_r2_r1-(x2-xi2).^2./(4*r1_sq)+...
                        (x2.^2+xi2.^2-4*(1-nu)*xi2.*(x2+xi2))./(4*r2_sq)+x2.*xi2.*(x2+xi2).^2./r2_sq.^2));
                    u2(:,ind:ind+nsources-1) = u2(:,ind:ind+nsources-1)+...
                        ((-s1/(pi*(1-nu))).*(((1-2*nu)/4).*log_r2_r1+(x2-xi2).^2./(4*r1_sq)-...
                        ((x2+xi2).^2-2*xi2.^2-2*(1-2*nu)*xi2.*(x2+xi2))./(4*r2_sq)+...
                        x2.*xi2.*(x2+xi2).^2./r2_sq.^2)+(s2/(pi*(1-nu))).*(((1-nu)/2).*(theta_diff)+...
                        x1.*(x2-xi2)./(4*r1_sq)-x1.*(x2+(3-4*nu).*xi2)./(4*r2_sq)-xi2.*x2.*x1.*(x2+xi2)./r2_sq.^2));
                end
                slip0_mask(:,ind:ind+nsources-1) = (s1==0 & s2==0); %If slip is 0, displacement is 0, but calling the equations gives NaNs.
            end
            u1(slip0_mask) = 0;
            u2(slip0_mask) = 0;
        end
    end
end