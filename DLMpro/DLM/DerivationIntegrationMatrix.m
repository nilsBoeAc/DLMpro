function [W,B] = DerivationIntegrationMatrix(k,RefChord,Panel)
% DerivationIntegraionMatrix calculates Derivation matrix W and
% Integration matrix B
%% Description:
%   The Derivation matrix describes the local deflections at center point corresponding to the
%   local aerodynamic degrees of freedom at AIC control point(3/4 chord).
%   The degrees of freedom for a wing without control surfaces are plunging and
%   pitching. The value for plunging is specified as: plg = i*(kr/br) with br as
%   half-chord of wing. The value for pitching is specified as:
%   alpha = 1+i*(kr/br)*(xcp-xc) where (xcp-xc) is the distance of control point
%   and center point in x direction.
% 
%   The Integration matrix relates the pressures to the center point forces and
%   moments. Deflections on 1/4 chord line leads to forces and moments. Plunging
%   causes forces and pitching causes moments. The value for plunging is specified
%   as the panel area A. The value for pitching is specified as: alpha = A*(xc-xq)
%   where (xc-xq) is the distance of control point and 1/4 cord line in x direction.
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input values:
%   k: Reduced frequency
%   RefChord: Reference chord of the wing
%   Panel: Object with panel information
%
%   Name-Value-Pairs:
%   plg: value for plunging
%   alpha: value for pitching
%   c: combined matrix of plunging and pitching for Derivation matrix
%   W2: reshaped c matrix
%   i: loop iteration counter
%   l: column counter
%   A: panel Area
%   B_v: combined matrix of plunging and pitching for Integration matrix
%
%% Output:
%   W: derivation matrix
%   B: Integration Matrix
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

%% Derivation Matrix

W = zeros(length(Panel.A)); %preallocation for W matrix
W = [W;W];  %increase W matrix to degrees of freedom
[~,n]=size(W);
plg = 1j*(k*2/RefChord); %plunging for W matrix
plg(1:length(Panel.A)) = plg; %build vector
alpha = 1+plg.*0.25.*Panel.chord';  %pitching for W matrix
c=[plg;alpha];  %concatenate pitching and plunging
W2 = reshape(c,[],1);   %reshape c to correct order of values

l=1;

for i=1:n
    %set values from W2 vector to the correct place in W matrix
    W(l,i)= -W2(l);
    W(l+1,i)=W2(l+1);
    l=l+2;
end

W = W.';

%% Integration Matrix

B = zeros(length(Panel.A));  %preallocation of B matrix
B = [B;B];  %increase B matrix to degrees of freedom
[~,n]=size(B);

A = Panel.A'; %Panel area is equal to influence of plunging
alpha = (A.*0.25.*Panel.chord');  %pitching influence
B_v = reshape([A;alpha],[],1);  %concatenated plunging and pitching influence

l = 1;

for i=1:n
    %set values from B_v vector to the correct place in B matrix
    B(l,i)= B_v(l);
    B(l+1,i)=B_v(l+1);
    l=l+2;
end
end

