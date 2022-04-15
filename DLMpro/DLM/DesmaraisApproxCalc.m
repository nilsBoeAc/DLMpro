function [I1, I2] = DesmaraisApproxCalc(u1,k1)
% DesmaraisApproxCalc calculates integrals I1 and I2 with Desmarais
% Approximation
%% Description:
%   Desmarais approximation has 12 iteration steps. Each step has a different
%   constant a. This constant is taken from function LookUpTable() where the
%   constants are saved. With iteration the sub integrals I0 and J0 are
%   calculated. These sub integrals are used to calculate integrals I1 and I2
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input values:
%   u1, k1: Factors of kernel function
%
%   Name-Value-Pairs:
%   i: imgainary unit
%   b & m: constant for Desmarais approximation
%   I0: sub integral for calculation of whole integral I1 (and I2)
%   J0: sub integral for calculation of whole integral I2
%   n: iteration of Laschka approximation
%   a: constants for Laschka approximation depending on iteration number
%   I1 and I2: approximated integrals
%
%% Output:
%   I1 and I2: approximated integrals
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

%% ToDo / Changelog
% -

%% Preallocation

i = sqrt(-1);
b = 0.009054814793;
m = 1;
I0 = zeros(size(u1));
J0 = zeros(size(u1));

%% Calculation

for n=1:12
    a = LookUpTable(n,'Desmarais');
    I0_step = (a.*exp(-2.^(n./m).*b.*u1))./((2.^(n./m)).^2.*b.^2+k1.^2).*(2.^(n./m).*b-i.*k1);
    I0 = I0+I0_step;
    J0_step = (a.*exp(-2.^(n./m).*b.*u1))./(((2.^(n./m)).^2.*b.^2+k1.^2).^2).*((2.^(n./m)).^2.*b.^2-k1.^2+2.^(n/m).*b.*u1.*((2.^(n/m)).^2.*b.^2+k1.^2)-i.*k1.*(2*2.^(n/m).*b+u1.*((2.^(n/m)).^2.*b.^2+k1.^2)));
    J0 = J0+J0_step;
end

I1 = (1-(u1./(1+u1.^2).^0.5)-i.*k1.*I0).*exp(-i.*u1.*k1);
I2 = ((2+i.*k1.*u1).*(1-(u1./((1+u1.^2).^0.5)))-(u1./((1+u1.^2).^1.5))-i.*k1.*I0+k1.^2.*J0).*(1/3).*exp(-i.*u1.*k1);

end

