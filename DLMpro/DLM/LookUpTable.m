function [a] = LookUpTable(n, Approx)
% Coefficients of Laschka and Desmarais Approximation are selected
%% Description:
%   Depending on the iteration counter 'n' an approximation coefficient
%   a is selected. 
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input values:
%   n: Iteration number
%   Approx: Selected approximation
%
%% Output:
%   a: constant for Laschka approximation depending on iteration number
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
% Copyright (c) 2021

%% ToDo / Changelog
% -

%% Switch of approximation method

switch Approx
    case 'Laschka'
        switch n
            case 1
                a = 0.24186198;
            case 2
                a = -2.7918027;
            case 3
                a = 24.991079;
            case 4
                a = -111.59196;
            case 5
                a = 271.43549;
            case 6
                a = -305.75288;
            case 7
                a = -41.183630;
            case 8
                a = 545.98537;
            case 9
                a = -644.78155;
            case 10
                a = 328.72755;
            case 11
                a = -64.279511;
        end
    case 'Desmarais'
        switch n
            case 1
                a = 0.000319759140;
            case 2
                a = -0.000055461471;
            case 3
                a = 0.002726074362;
            case 4
                a = 0.005749551566;
            case 5
                a = 0.031455895072;
            case 6
                a = 0.106031126212;
            case 7
                a = 0.406838011567;
            case 8
                a = 0.798112357155;
            case 9
                a = -0.417749229098;
            case 10
                a = 0.077480713894;
            case 11
                a = -0.012677284771;
            case 12
                a = 0.001787032960;
        end
end

