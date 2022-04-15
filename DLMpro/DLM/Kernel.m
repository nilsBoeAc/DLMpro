function [P1,P2] = Kernel(Mach,k,x,y,z,e,tanLambda,relGamma,method,dsp)
% This function solves the Kernel function for a specific
% doublet-line-point and control point pair
%% Description:
%   First the coefficients for the specific point pair are determined. The
%   integrals are deteminded are dependent on the desired approximateion
%   method. At least the solutions P1 and P2 are calculated.
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input variables:
%   Mach: Mach number
%   k: Reduced frequency
%   x, y, z: Matrix of distances between sending and receiving point on
%   each panel
%   e: Matrix of panel half spans
%   with zeros used for kernel evaluation
%   tanLambda: Vector of sweep angles for each panel
%   relGamma: Matrix of relative dihedrals of panels
%   
%   Optional input variables
%   method: Method of integral approximation inside kernel function.
%   Watkins, Laschka or Desmarais
%   dsp: Flag if status shall be shown in command window
%
%% Output:
%   P1 & P2: Solutions of Kernel function for specific point
%
%% References:
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

%% Coefficients

r1 = sqrt((y-e).^2+z.^2);
beta2 = 1-Mach^2;
R = sqrt((x-e.*tanLambda').^2+beta2*r1.^2);
u1 = (Mach*R-x+e.*tanLambda')./(beta2*r1);
k1 = k*r1;
i=sqrt(-1);
eiku = exp(-i.*k1.*u1);

%% Calculation
T1 = cos(relGamma);
T2 = z.*(z.*cos(relGamma)+(y-e).*sin(relGamma));

%selecting approximation methods
switch method
    case "Watkins"
        [I1, I2] = WatkinsApprox(u1,k1,dsp);
    case "Laschka"
        [I1, I2] = LaschkaApprox(u1,k1,dsp);
    case "Desmarais"
        [I1, I2] = DesmaraisApprox(u1,k1,dsp);
end

K1 = -I1-(eiku.*Mach.*r1)./(R.*sqrt(1+u1.^2));
K2 = +3*I2 + ((i.*k1.*eiku.*(Mach.^2).*(r1.^2))./(R.^2.*sqrt(1+u1.^2))) + (eiku.*Mach.*r1)./(R.*(1+u1.^2).^(3/2))...
     .*((1+u1.^2).*beta2.*r1.^2./R.^2+2+Mach.*r1.*u1./R);

K10 = -1-(x-e.*tanLambda')./R;
K20 = 2+((x-e.*tanLambda')./R).*(2+beta2.*r1.^2./(R.^2));

ind1 = r1==0 & x>=0;
ind2 = r1==0 & x<0;

K1(ind1==1) = -2;
K2(ind1==1) = 4;

K1(ind2==1) = 0;
K2(ind2==1) = 0;

P1 = -(K1.*exp(-i.*k.*(x-e.*tanLambda'))-K10).*T1;
P2 = -(K2.*exp(-i.*k.*(x-e.*tanLambda'))-K20).*T2;

if dsp
    dt = datestr(now,'HH:MM:SS.FFF');
    disp(dt+" Calculation of P1 and P2: Done")
end
end

