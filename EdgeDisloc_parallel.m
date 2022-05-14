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

classdef EdgeDisloc_parallel
    methods (Static = true)
        function [sigma11,sigma12,sigma22] = Stress(xi1_vec,xi2_vec,x1_vec,x2_vec,s1_vec,s2_vec,nu,mu)
            %This function models stress due to slip over a confined
            %interval, represented by edge dislocations at each end with
            %opposite signs. %It uses Equantions 3.77, 3.78, and 3.79 of
            %Segall (2010) as corrected in the Errata to that book.
            %Incides 1 and 2 correspond to the x and y directions respectively.
            %xi1_vec and xi2_vec are 2xN arrays of source segment endpoints.
            %x1_vec and x2_vec are 1xM arrays of observation points. s1_vec
            %and s2_vec are 1xM arrays of slip on each source segment.
            %nu is Poisson's ratio. mu is shear modulus.
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
            sigma11 = zeros(M,N);
            sigma12 = zeros(M,N);
            sigma22 = zeros(M,N);
            [xi1_t,xi2_t] = deal(xi1_vec(1,:),xi2_vec(1,:)); %Top point coordinates.
            [xi1_b,xi2_b] = deal(xi1_vec(2,:),xi2_vec(2,:)); %Bottom point coordinates.
            parfor n = 1:N
                %Set up the arrays of each source with each observation point.
                x1_unshifted = x1_vec';
                x2 = x2_vec';
                s1 = s1_vec(n);
                s2 = s2_vec(n);
                %Calculate the stresses as the sum of two edge dislocations, which cancel except along the element:
                for i = 1:2 %First edge dislocation at the tops of each element, then at the bottoms with opposite slip.
                    if i == 1
                        xi1 = xi1_t(n);
                        xi2 = xi2_t(n);
                        %Flipping the signs on s if the fault dips left only has to be done once.
                        if sign(s1)==sign(s2) || s2==0 
                            %This occurs if the fault dips to the left.
                            %s2==0 is horizontal, which is considered
                            %left-dipping because the point on the right is
                            %considered the top point (see above).
                            [s1,s2] = deal(-s1,-s2); %Flip to the sign convention that Segall's equations use.
                        end
                    else
                        xi1 = xi1_b(n);
                        xi2 = xi2_b(n);
                        [s1,s2] = deal(-s1,-s2);
                    end
                    %Solve the Edge Dislocation stress equations.
                    %For cpu, it's better to use the
                    %optimized versions with intermediate variables.
                    x1 = x1_unshifted-xi1;%These equations are written for xi1 = 0. If that's not the case, this makes it equivalent to that case.
                    [sigma11_part,sigma22_part,sigma12_part] = EdgeDisloc_parallel.OptimizedStressEqns(xi2,x1,x2,s1,s2,nu,mu);
                    [sigma11(:,n),sigma22(:,n),sigma12(:,n)] = ...
                        deal(sigma11(:,n)+sigma11_part,sigma22(:,n)+sigma22_part,sigma12(:,n)+sigma12_part);
                end
            end
        end
        
        function [u1,u2] = Displacement(xi1_vec,xi2_vec,x1_vec,x2_vec,s1_vec,s2_vec,nu)
            %This function models displacement due to slip over a confined
            %interval, represented by edge dislocations at each end with
            %opposite signs. %It uses Equantion 3.74 of
            %Segall (2010) as corrected in the Errata to that book.
            %Incides 1 and 2 correspond to the x and y directions respectively.
            %xi1_vec and xi2_vec are 2xN arrays of source segment endpoints.
            %x1_vec and x2_vec are 1xM arrays of observation points. s1_vec
            %and s2_vec are 1xM arrays of slip on each source segment.
            %nu is Poisson's ratio.
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
            u1 = zeros(M,N);
            u2 = zeros(M,N);
            [xi1_t,xi2_t] = deal(xi1_vec(1,:),xi2_vec(1,:)); %Top point coordinates.
            [xi1_b,xi2_b] = deal(xi1_vec(2,:),xi2_vec(2,:)); %Bottom point coordinates.
            parfor n = 1:N
                %Set up the arrays of each source with each observation point.
                x1_unshifted = x1_vec';
                x2 = x2_vec';
                s1 = s1_vec(n);
                s2 = s2_vec(n);
                if s1==0 && s2==0 %If slip is 0, displacement is 0, but calling the equations gives NaNs.
                    [u1(:,n),u2(:,n)] = deal(0,0);
                else
                    for i = 1:2 %First edge dislocation at the tops of each element, then at the bottoms with opposite slip.
                        if i == 1
                            xi1 = xi1_t(n);
                            xi2 = xi2_t(n);
                            %Flipping the signs on s if the fault dips left only has to be done once.
                            if sign(s1)==sign(s2) || s2==0
                                %This occurs if the fault dips to the left.
                                %s2==0 is horizontal, which is considered
                                %left-dipping because the point on the right is
                                %considered the top point (see above).
                                [s1,s2] = deal(-s1,-s2); %Flip to the sign convention that Segall's equations use.
                            end
                        else
                            xi1 = xi1_b(n);
                            xi2 = xi2_b(n);
                            [s1,s2] = deal(-s1,-s2);
                        end
                        %Solve the Edge Dislocation stress equations.
                        x1 = x1_unshifted-xi1;%These equations are written for xi1 = 0. If that's not the case, this makes it equivalent to that case.
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
                        %For cpu, it's better to use the
                        %optimized versions with intermediate variables.
                        [u1_part,u2_part] = EdgeDisloc_parallel.OptimizedDisplacementEqns(xi2,x1,x2,s1,s2,nu,theta_diff);
                        [u1(:,n),u2(:,n)] = deal(u1(:,n)+u1_part,u2(:,n)+u2_part);
                    end
                end
            end
        end
        
        function [sigma11,sigma22,sigma12] = OptimizedStressEqns(xi2,x1,x2,s1,s2,nu,mu)
            %STRESSEQNS
            %    [SIGMA11,SIGMA22,SIGMA12] = STRESSEQNS(XI2,X1,X2,S1,S2,NU,MU)
            
            %    This function was generated by the Symbolic Math Toolbox version 8.6.
            %    03-May-2021 11:56:50
            
            t2 = x2+xi2;
            t4 = x1.^2;
            t6 = xi2.*2.0;
            t7 = 1.0./pi;
            t8 = nu-1.0;
            t9 = -xi2;
            t5 = t4.^2;
            t10 = t4.*3.0;
            t11 = t2.^2;
            t12 = t2.^3;
            t14 = t9+x2;
            t16 = 1.0./t8;
            t19 = t2.*t4.*x2.*6.0;
            t15 = t11.*3.0;
            t17 = t14.^2;
            t18 = t4+t11;
            t20 = -t19;
            t25 = t12.*t14;
            t21 = t4+t17;
            t23 = 1.0./t18.^2;
            t24 = 1.0./t18.^3;
            t30 = t5+t20+t25;
            t28 = 1.0./t21.^2;
            t29 = -t23.*x1.*(t4-t11);
            t31 = -t2.*t23.*(t4-t11);
            t32 = t23.*x1.*(t4-t11);
            t37 = t6.*t24.*t30;
            t34 = -t28.*x1.*(t4-t17);
            sigma11 = (mu.*s1.*t7.*t16.*(t37+t2.*t23.*(t10+t11)-t14.*t28.*(t10+t17)))./2.0-(mu.*s2.*t7.*t16.*(t32+t34+t24.*x1.*xi2.*(t4.*(t2.*2.0+x2)+t11.*(t6-x2)).*4.0))./2.0;
            if nargout > 1
                t36 = t14.*t28.*(t4-t17);
                sigma22 = (mu.*s1.*t7.*t16.*(t31+t36+t6.*t24.*(t5+t19-t12.*(t2+x2.*2.0))))./2.0+(mu.*s2.*t7.*t16.*(-t23.*x1.*(t4+t15)+t28.*x1.*(t17.*2.0+t21)+t24.*x1.*x2.*xi2.*(t4-t15).*4.0))./2.0;
            end
            if nargout > 2
                sigma12 = (mu.*s1.*t7.*t16.*(t29+t28.*x1.*(t4-t17)+t24.*x1.*x2.*xi2.*(t4-t15).*4.0))./2.0+(mu.*s2.*t7.*t16.*(t31+t36+t37))./2.0;
            end
        end
        
        function [u1,u2] = OptimizedDisplacementEqns(xi2,x1,x2,s1,s2,nu,theta_diff)
            %DISPLACEMENTEQNS
            %    [U1,U2] = DISPLACEMENTEQNS(XI2,X1,X2,S1,S2,NU,THETA_DIFF)
            
            %    This function was generated by the Symbolic Math Toolbox version 8.6.
            %    03-May-2021 11:57:39
            
            t2 = x2+xi2;
            t3 = nu.*4.0;
            t4 = x1.^2;
            t5 = xi2.^2;
            t6 = 1.0./pi;
            t7 = nu-1.0;
            t8 = -xi2;
            t11 = nu./2.0;
            t9 = t4.*4.0;
            t10 = t2.^2;
            t12 = t3-3.0;
            t13 = t8+x2;
            t15 = 1.0./t7;
            t20 = t11-1.0./2.0;
            t21 = t11-1.0./4.0;
            t14 = t10.*4.0;
            t17 = t13.^2;
            t18 = t4+t10;
            t19 = t8.*t12;
            t23 = t20.*theta_diff;
            t22 = t17.*4.0;
            t24 = t4+t17;
            t25 = 1.0./t18.^2;
            t26 = t19+x2;
            t27 = sqrt(t18);
            t28 = t9+t14;
            t29 = t9+t22;
            t30 = 1.0./t28;
            t31 = 1.0./sqrt(t24);
            t32 = t2.*t25.*x1.*x2.*xi2;
            t35 = t8.*t10.*t25.*x2;
            t33 = 1.0./t29;
            t38 = t26.*t30.*x1;
            t39 = t27.*t31;
            t36 = t13.*t33.*x1;
            t37 = t17.*t33;
            t40 = log(t39);
            t41 = t21.*t40;
            u1 = s2.*t6.*t15.*(t35+t37+t41-t30.*(t5+x2.^2+t2.*xi2.*(t12-1.0)))+s1.*t6.*t15.*(t23+t32+t36-t38);
            if nargout > 1
                u2 = s2.*t6.*t15.*(t23+t32-t36+t38)-s1.*t6.*t15.*(t35-t37+t41+t30.*(t5.*-2.0+t10+t2.*xi2.*(t3-2.0)));
            end
        end
    end
end