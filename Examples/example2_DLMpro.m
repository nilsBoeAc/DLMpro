
% Example 2
% 
%% Description
% This script, contains an example for DLMpro.
% It contains:
%   - define wing and create object
%   - calc AIC with one pair of method
%   - define Transformatino matrix in case of different orientation of
%     structural grid
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
NS    = 20;
NC    = 8;
sweep = 0;      % [degree]
dihedral = 0;   % [degree]
taperRatio =1;  % [-]

%% Spline / Transformation Matrix
% Assume that a structural beam model is available wiht 50 nodes and each
% has 3DoF. By using an appropriate spline approach, the following spline
% matrix is generated.
% Note: This Matrix assumes the coordinate systems of aero and structure
%       have the same orientations.
GSA = load('GSA.mat').GAK;

% The different orientation of the coordinate systems can be taken into
% account by defining the "CSChange" property using Coordinate rotation
% In this example the structural CS defined as follows (total 180 degree around y-axis):
%   x_struct = -x_aero
%   y_struct =  y_aero
%   z_struct = -z_aero
% Thus, a rotation about the z-axis of 180 degree and then a 180 degree
% rotation around the x'-axis is required -> see coordinate rotation
CSChange = [180,180,0]; 

%% Object creation
offset = [-0.625,0,0];
dlm = DLMpro(span,chord,NS,NC,sweep,dihedral,taperRatio,offset,'CSChange',CSChange);

%% Calc AIC with given red. frequencies and mach numbers
k = 0.1:0.1:1;     % reduced Frequency
ma = 0.1:0.1:0.4;  % mach Number

res1 = dlm.calcAIC(k,ma,"GSA",GSA,"verbose",0);

%% Adjust Method (Integration and Approximation)
int = "Parabolic";
app = "Watkins";
res2 = dlm.calcAIC(k,ma,'GSA',GSA,'Integration',int,'Approximation',app);

int = "Parabolic";
app = "Desmarais";
res3 = dlm.calcAIC(k,ma,'GSA',GSA,'Integration',int,'Approximation',app);

int = "Quartic";
app = "Laschka";
res4 = dlm.calcAIC(k,ma,'GSA',GSA,'Integration',int,'Approximation',app);

int = "Quartic";
app = "Watkins";
res5 = dlm.calcAIC(k,ma,'GSA',GSA,'Integration',int,'Approximation',app);

int = "Quartic";
app = "Desmarais";
res6 = dlm.calcAIC(k,ma,'GSA',GSA,'Integration',int,'Approximation',app);


%% Access results
dlm.resultOverview

% Result No: Number of Result
% k-M-List : Pair for which the AIC are calculated
% AIC      : Aerodynamic Influence Coefficients
% Qss      : Aerodyanmic Influence with respect to 2 DoF per panel (pitch and plunge) - in structure coordinate System
% Qkk      : Aerodyanmic Influence with respect to 2 DoF per panel (pitch and plunge) - in aero coordinate System
% other    : information on used methods
%
% Note: Since GSA (Transformation matrix between structural and aeor grid)
%       is given, the size and values of Qss and Qkk are not equal.
%       Also "CSChange" Parameter is given, so that some of the coefficients changes in sign, due to the different orientation
%       Equation : Qss = [GSA]*[TR]*[Qkk]*[TR]*[GSA]'; TR: Matrix for orientation change
%
% Comment: The differentiation of GSA and CSChange is done, due to the
%          capability of visualization of the results to a later stage and due to requirements
%          by the auhtors internal codes.
%          The influence of "CSChange" is currently kind of "hard-code" - use this with caution!