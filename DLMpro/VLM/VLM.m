function D_vlm = VLM(Panel,Mach,dsp)
% VLM calulates the steady influence of the AIC
%% Description:
%   The VLM calculates the steady part of the influence coefficient.
%   First, the Prandtl-Glauert-Transformation is applied by dividing the
%   x-coordinates of points with beta. This is used as correction of the compressibility
%   of air. In further progress several geometries are calculated to determine
%   the horseshoe vortex and its influence. D1 is the calulation for the bound
%   vortex, D2 and D3 are the calculations for the border vortices. Each of this
%   influences is weighted by a normalized normalvector for use with dihedral.
%   The influences are summed up. In order to the relation of the influence matrix
%   to pressure coefficients the influences have to be multiplied with panel
%   half-chord. 
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input variables:
%   Panel: Object of panel information
%   Mach: Mach number
%   
%   Optional input variables:
%   dsp: Flag for display status in command window
%
%   Name-Value-Pairs:
%   Variable Explanation:
%   beta: Prandtl-Glauert-Transformation
%   CP_b = Control point coordinates
%   D1_b and D5_b: border points of bound vortex line
%   N_y and N_z: normalized normal vectors in y- and z-direction
%   epsilon: value for recognition of other values as 0 to prevent singularities
%   r0: distance of border points
%   r1 and r2: distances from control point to border points
%   r1Xr2: crossproduct of r1 and r2
%   abs_r1Xr2: amount of vector r1Xr2
%   abs_r1: amount of vector r1
%   abs_r2: amount of vector r2
%   r0r1 and r0r2: scalar products of r0 to r1 and r2
%   D1: influences of bound vortex
%   D2 and D3: influences of border vortices
%   D: steady influence matrix related to circulation
%   F: Factor to relate influences to pressure coefficients
%
%% Output:
%   D_vlm: Resulting steady influence matrix
%
%% References:
%   [Ref1]
%   [Ref2]
%
%% Disclaimer:
%
% Last editor:  Marc Bangel
% Last edit on: 22.12.2021
% Code version: X.Y.Z
%
%% Disclaimer
%   Copyright (c) 2021 Nils Böhnisch, Marc Bangel.
%
%   This file is part of DLMPro.
%
%   DLMPro is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version. Also see the file "License".

%% ToDo / Changelog
% -

%% Prandtl-Glauert-Transformation

beta = sqrt(1-Mach^2);

%% Multiply y-coordinates by beta for compressibility consideration

%Corrected control point
CP_b = Panel.CP;
CP_b(:,1) = CP_b(:,1)./beta;
%corrected border point
D1_b = Panel.D1;
D1_b(:,1) = Panel.D1(:,1)./beta;
%corrected border point
D5_b = Panel.D5;
D5_b(:,1) = Panel.D5(:,1)./beta;

%% Normal vector initialization

N_y = repmat(Panel.N(:,2),1,length(Panel.N(:,2)));
N_z = repmat(Panel.N(:,3),1,length(Panel.N(:,3)));

%% Calculate r0, r1 and r2

epsilon = 10e-6;

%vector from inner to outer point on bound vortex line
r0x = repmat(D5_b(:,1)-D1_b(:,1),1,length(D5_b))';
r0y = repmat(D5_b(:,2)-D1_b(:,2),1,length(D5_b))';
r0z = repmat(D5_b(:,3)-D1_b(:,3),1,length(D5_b))';

%vector from inner bound vortex point to control point
r1x = CP_b(:,1)-D1_b(:,1)';
r1y = CP_b(:,2)-D1_b(:,2)';
r1z = CP_b(:,3)-D1_b(:,3)';

%vector from outer bound vortex point to control point
r2x = CP_b(:,1)-D5_b(:,1)';
r2y = CP_b(:,2)-D5_b(:,2)';
r2z = CP_b(:,3)-D5_b(:,3)';

%% Calculate D1

%scalar products of vectors from border voints to control point
r1Xr2_x = r1y.*r2z-r1z.*r2y;
r1Xr2_y = r1z.*r2x-r1x.*r2z;
r1Xr2_z = r1x.*r2y-r1y.*r2x;

%amount of vectors
abs_r1Xr2 = sqrt(r1Xr2_x.^2+r1Xr2_y.^2+r1Xr2_z.^2);
abs_r1 = sqrt(r1x.^2+r1y.^2+r1z.^2);
abs_r2 = sqrt(r2x.^2+r2y.^2+r2z.^2);

%scalar products of vectors border points to control point and bound vortex
%length vector
r0r1 = r0x.*r1x+r0y.*r1y+r0z.*r1z;
r0r2 = r0x.*r2x+r0y.*r2y+r0z.*r2z;

%calculate first part of influence matrix
one = ones(size(r1Xr2_x));
D1_base = (1/(4*pi))*(one./(abs_r1Xr2.^2)).*((r0r1./abs_r1) - (r0r2./abs_r2));
D1_u = r1Xr2_x.*D1_base;
D1_v = r1Xr2_y.*D1_base;
D1_w = r1Xr2_z.*D1_base;

ind =  find(abs_r1<epsilon);
D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

ind =  find(abs_r2<epsilon);
D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

ind =  find(abs_r1Xr2<epsilon);
D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

D1 = D1_v.*N_y + D1_w.*N_z;

%% Calculate D2

d2 = sqrt(r1y.^2+r1z.^2); %distance from control point to edge vortex
cos_b1 = 1;
cos_b2 = -r1x./abs_r1;
cosGamma = r1y./d2;
sinGamma = -r1z./d2;

D2_base = (-1)/(4*pi)*(cos_b1-cos_b2)./d2;
D2_v = D2_base.*sinGamma;
D2_w = D2_base.*cosGamma;

ind =  find(abs_r1<epsilon);
D2_v(ind) = 0;
D2_w(ind) = 0;

ind =  find(d2<epsilon);
D2_v(ind) = 0;
D2_w(ind) = 0;

D2 = D2_v.*N_y + D2_w.*N_z;

%% Calculate D3

d3 = sqrt(r2y.^2+r2z.^2); %distance from control point to edge vortex
cos_b1 = r2x./abs_r2;
cos_b2 = -1;
cosGammaD3 = -r2y./d3;
sinGammaD3 = r2z./d3;

D3_base = (-1)/(4*pi)*(cos_b1-cos_b2)./d3;
D3_v = D3_base.*sinGammaD3;
D3_w = D3_base.*cosGammaD3;

ind =  find(abs_r2<epsilon);
D3_v(ind) = 0;
D3_w(ind) = 0;

ind =  find(d3<epsilon);
D3_v(ind) = 0;
D3_w(ind) = 0;

D3 = D3_v.*N_y + D3_w.*N_z;

%% Summarize D1, D2, D3

D = D1 + D2 + D3;

F = 0.5*Panel.chord;
F = F';
F = repmat(F,length(F),1);

D_vlm = D.*F;

if (dsp)
    dt = datestr(now,'HH:MM:SS.FFF');
    disp(dt+" VLM Matrix: Done")
end

end