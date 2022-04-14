function [I1,I2] = WatkinsApproxCalc(u1,k1)
% WatkinsApproxCalc calculates integrals I1 and I2 with Watkins
% Approximation
%% Description:
%   Watkins approximation is done by substituting the term in the integral by
%   a function with 6 constants. The integration of this function leads to
%   the sub integrals I0 and J0 that are used to calculate integrals I1 and
%   I2.
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
%   a's & b's: constants for Watkins approximation
%   I0: sub integral for calculation of whole integral I1 (and I2)
%   J0: sub integral for calculation of whole integral I2
%
%% Output:
%   I1 and I2: approximated integrals
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

%% Constants

a1 = 0.101;
a2 = 0.899;
a3 = 0.09480933;
b1 = 0.329;
b2 = 1.4067;
b3 = 2.90;
i = sqrt(-1);

%% Calculation

I0 = ((a1*exp((-b1-(i*k1)).*u1)./(b1+(i*k1))) + (a2*exp((-b2-(i*k1)).*u1)./(b2+(i*k1)))+...
     (a3./(((b3 + (i*k1)).^2) + (pi^2))).*(((b3+(i*k1)).*sin(pi.*u1)) + (pi*cos(pi.*u1))).*exp((-b3-(i*k1)).*u1)).*exp(i.*k1.*u1);

J0 = (((a1.*((i.*k1+b1).*u1+1).*exp((-b1-i.*k1).*u1))./(i.*k1+b1).^2)+((a2.*((i.*k1+b2).*u1+1).*exp((-b2-i.*k1).*u1))./(i.*k1+b2).^2)...
     +(a3.*exp(-(k1.*i+b3).*u1).*(((k1.*i+b3).*(((k1.*i+b3).^2+pi^2).*u1+k1.*i+b3)-pi^2).*sin(pi*u1)+pi*(((k1.*i+b3).^2+pi^2).*u1+2*(k1.*i+b3)).*cos(pi.*u1)))./((k1.*i+b3).^2+pi^2).^2).*exp(i.*u1.*k1);
 
I1 = (1-(u1./(1+u1.^2).^0.5)-i.*k1.*I0).*exp(-i.*u1.*k1);
I2 = ((2+i.*k1.*u1).*(1-(u1./((1+u1.^2).^0.5)))-(u1./((1+u1.^2).^1.5))-i.*k1.*I0+k1.^2.*J0).*(1/3).*exp(-i.*u1.*k1);

end