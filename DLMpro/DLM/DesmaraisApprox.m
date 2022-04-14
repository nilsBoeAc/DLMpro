function [I1, I2] = DesmaraisApprox(u1,k1,dsp)
% DesmaraisApprox calculates integrals I1 and I2 with Desmarais
% Approximation
%% Description:
%   Desmarais approximation is calculated in DesmaraisApproxCalc(). This
%   approximation is only valid for u1>=0. To approximate the integrals if
%   u1<0 the approximation has to be done for different values and has to be
%   calculated in the end.
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input values:
%   u1, k1: Factors of kernel function
% 
%   Optional input variables:
%   dsp: Flag for display status in command window
% 
%   Name-Value-Pairs:
%   i: imgainary unit
%   ind1: index of positiv u1 values in matrix
%   ind2: index of negativ u1 values in matrix
%   u1_t1: matrix only with positiv u1 values
%   u1_t2: matrix only with negativ u1 values
%   k1_t1: matrix of k1 if u1 is positiv
%   k1_t2: matrix of k1 if u1 is negativ
%   u1_0: matrix of zeros for u1
%   I1_pos & I2_pos: approximated integral I1 for positiv u1 values
%   I1_0 & I2_0: approximated integral I1 for u1 = 0
%   I1_neg & I2_neg: approximated integrals for negative u1 values
% 
%% Output:
%   I1 and I2: approximated integrals for u1>=0 and u1<0
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

%% Timestamp in command window

if dsp
    dt = datestr(now,'HH:MM:SS.FFF');
    disp(dt+" Start Desmarais approximation...")
end

%% Preallocation

I1 = zeros(size(u1));
I2 = zeros(size(u1));

u1_t1 = zeros(size(u1));
u1_t2 = zeros(size(u1));
k1_t1 = zeros(size(u1));
k1_t2 = zeros(size(u1));
u1_0 = zeros(size(u1));
i = sqrt(-1);

%% Calculate I1 & I2 for u >= 0

ind1 = find(u1>=0);
u1_t1(ind1) = u1(ind1);
k1_t1(ind1) = k1(ind1);

[I1_pos,I2_pos] = DesmaraisApproxCalc(u1_t1,k1_t1);

I1(ind1) = I1_pos(ind1);
I2(ind1) = I2_pos(ind1);


%% Calculate I1 & I2 for u < 0

ind2 = find(u1<0);
u1_t2(ind2) = -u1(ind2);
k1_t2(ind2) = k1(ind2);

[I1_0,I2_0] = DesmaraisApproxCalc(u1_0,k1_t2);
[I1_neg,I2_neg] = DesmaraisApproxCalc(u1_t2,k1_t2);

I1(ind2) = (2*real(I1_0(ind2))) - real(I1_neg(ind2)) + (i*imag(I1_neg(ind2)));
I2(ind2) = (2*real(I2_0(ind2))) - real(I2_neg(ind2)) + (i*imag(I2_neg(ind2))); 

if dsp
    dt = datestr(now,'HH:MM:SS.FFF');
    disp(dt+" Desmarais Approximation: Done")
end

end