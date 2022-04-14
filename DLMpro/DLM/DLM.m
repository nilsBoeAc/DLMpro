function [D_dlm] = DLM(Panel,Mach,k,Integration,Approximation,pkern,dsp)
% First level function of DLM influence matrix calculation
%% Description:
%   First important geometric values are calculated and and transformation
%   of panel coordinates from global to panel coordinate system is
%   performed. Depending on desired Integration method the corresponding
%   function is called.
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input values:
%   Panel: Object with panel information
%   Mach: Mach number
%   k: Specific reduced frequency for DLM
% 
%   Optional input values
%   Integration: Integration method (Parabolic, Quartic)
%   Approximation: Approximation method (Laschka, Watkins, Desmarais)
%   pkern: flag if kernel function shall be plotted.
%   dsp: flag if status shall be shown in command window
% 
%% Output:
%   D_dlm: influence matrix of DLM part
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

%% Get distances

%distances from midspan point of doublet line to control point
x0= Panel.CP(:,1)-Panel.SP(:,1)';
y0 = Panel.CP(:,2)-Panel.SP(:,2)';
z0 = Panel.CP(:,3)-Panel.SP(:,3)';

sinGamma = (Panel.D5(:,3)-Panel.D1(:,3))./Panel.s; %dihedral sinus of each panel
cosGamma = (Panel.D5(:,2)-Panel.D1(:,2))./Panel.s; %dihedral cosinus of each panel
tanLambda = (Panel.D5(:,1)-Panel.D1(:,1))./Panel.s; %tanget of sweep of each panel
relGamma = asin(sinGamma)-asin(sinGamma)'; %matrix of relative dihedrals

%transformation of distances to element coordinate system(sending element)
x = x0;
y = y0.*cosGamma'+z0.*sinGamma';
z = z0.*cosGamma'-y0.*sinGamma';

e = Panel.s./2; %halfspan of panel
e_m = repmat(e,1,length(e)); %matrix of panel halfspan
e_2 = e_m./2; %halfspan of panel halfspan for quartic integration
e_0 = zeros(size(e_m));
%% Integration

switch Integration
    case "Parabolic"
        D_dlm = ParabolicApprox(Panel,Mach,k,x,y,z,e_m,e_0,tanLambda,relGamma,Approximation,pkern,dsp);
    case "Quartic"
        D_dlm = QuarticApprox(Panel,Mach,k,x,y,z,e_m,e_0,e_2,tanLambda,relGamma,Approximation,pkern,dsp);
end

%% Status display in command window
if (dsp)
    dt = datestr(now,'HH:MM:SS.FFF');
    disp(dt+" Calculation of D_dlm: Done")
end

end