
% Example 1 
% 
%% Description
% This script, contains an example for DLMpro.
% It contains:
%   - define wing and create object
%   - calc AIC with various Methods
%   - access results
%
%% Note
%   Please make sure, that DLMpro is installed.
%
%% Disclaimer
%   Copyright (c) 2021 Nils Böhnisch, Marc Bangel.
%
%   This file is part of DLMPro.
%
%   DLMPro is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version. Also see the file "License"

%% Initialize
clearvars;

%% Wing definition
span  = 5.7*2;  % [m]
chord = 1.25;   % [m]
NS    = 30;
NC    = 5;
sweep = 0;      % [degree]
dihedral = 0;   % [degree]
taperRatio =1;  % [-]

%% Object creation
offset = [-0.625,0,0];
dlm = DLMpro(span,chord,NS,NC,sweep,dihedral,taperRatio,offset);

%% Calc AIC with given red. frequencies and mach numbers
k = 0.1:0.1:1;      % reduced Frequency
ma = 0.1:0.1:0.6;  % mach Number
%ma = 0.2;

res1 = dlm.calcAIC(k,ma,'verbose',1);

%% Adjust Method (Integration and Approximation)
int = "Parabolic";
app = "Watkins";
res2 = dlm.calcAIC(k,ma,'Integration',int,'Approximation',app,'verbose',0);

int = "Parabolic";
app = "Desmarais";
res3 = dlm.calcAIC(k,ma,'Integration',int,'Approximation',app,'verbose',0);

int = "Quartic";
app = "Laschka";
res4 = dlm.calcAIC(k,ma,'Integration',int,'Approximation',app,'verbose',0);

int = "Quartic";
app = "Watkins";
res5 = dlm.calcAIC(k,ma,'Integration',int,'Approximation',app,'verbose',0);

int = "Quartic";
app = "Desmarais";
res6 = dlm.calcAIC(k,ma,'Integration',int,'Approximation',app,'verbose',1);


%% Access results
dlm.resultOverview
% Result No: Number of Result
% k-M-List : Pair for which the AIC are calculated
% AIC      : Aerodynamic Influence Coefficients
% Qss      : Aerodyanmic Influence with respect to 2 DoF per panel (pitch and plunge) - in structure coordinate System
% Qkk      : Aerodyanmic Influence with respect to 2 DoF per panel (pitch and plunge) - in aero coordinate System
% other    : information on used methods

% Note: Since GSA (Transformation matrix between structural and aeor grid)
%       is not given, the size and values of Qss and Qkk are equal.




