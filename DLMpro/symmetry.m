function Dnew = symmetry(obj,D)  
% symmetry calculates a new influence coefficients matrix due to symmetry effects
%% Description:
%   Symmetry is used to calculate a new influence coefficients matrix
%   that is used to determine the aerodynamic influence coefficients.
%   This is needed to consider symmetry effects and the influences from one
%   wing to another if only one wing shall be analyzed.
%   
%% Syntax:
%   Dnew = symmetry(obj,D);
%
%% Input:
%   Required input variables:
%     D: influence coefficients matrx
%   
%% Output:
%   Dnew: new influence coefficients matrix considered of symmetry effects
%
%% Disclaimer:
%
% Last editor:  Marc Bangel
% Last edit on: 17.01.2022
%   Copyright (c) 2021 Nils Böhnisch, Marc Bangel.
%
%   This file is part of DLMPro.
%
%   DLMPro is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version. Also see the file "License".


%% ToDo / Changelog
% -Improve kernel function plotting behaviour
NS_hw = obj.wingProp.NS;
NC = obj.wingProp.NC;
D_r = D(NS_hw*NC+1:end,NS_hw*NC+1:end);   % influence from right wing on itself
D_lr = D(NS_hw*NC+1:end,1:NS_hw*NC);      % influence from left wing on right wing
sub_D_lr = zeros(NS_hw*NC,NC);

%sorting the influence matrix from left wing on
%right wing so that it can easily summed up with
%influnece from right wing on itself
m=1;
l = NS_hw;
for k=1:NS_hw
    sub_D_lr(:,:,m) = D_lr(:,(k-1)*NC+1:k*NC);
    m = m+1;
end

for k=1:NS_hw
    new_D_lr(:,(k-1)*NC+1:k*NC) = sub_D_lr(:,:,l);
    l = l - 1;
end

Dnew = D_r+new_D_lr;