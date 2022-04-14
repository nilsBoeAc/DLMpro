function [D_dlm] = ParabolicApprox(Panel,Mach,k,x,y,z,e_m,e_0,tanLambda,relGamma,method,pkern,dsp)
% ParabolicApprox calculates the influence matrix if parabolic integration
% is selected
%% Description:
%   In the parabolic integration the kernel function is calculated three
%   times for three different points on the doublet line. The coefficients
%   for the approximation of the integral of the kernel function are
%   calculated. The approximation describes a parabolic equation. At least
%   the planar and non-planar parts of the influence matrix are calculated
%   ant put together.
%   
%% Syntax:
%   Here are the syntax
%
%% Input:
%   Required input variables:
%   Panel: Object with all panel information
%   Mach: Mach number
%   k: reduced frequency
%   x, y, z: matrix of distances between sending and receiving point on
%   each panel
%   e_m, e_0: matrix of panel half spans and a panel with halfspan matrix
%   with zeros used for kernel evaluation
%   tanLambda: vector of sweep angles for each panel
%   relGamma: matrix of relative dihedrals of panels
%   
%   Optional input variables
%   method: method of integral approximation inside kernel function.
%   Watkins, Laschka or Desmarais
%   pkern: flag if kernel function shall be plotted.
%   dsp: flag if status shall be shown in command window
%
%% Output:
%   D_dlm: influence matrix of DLM part
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
% -Improve kernel function plotting behaviour

%% Kernel function

[P1_e,P2_e] = Kernel(Mach,k,x,y,z,e_m,tanLambda,relGamma,method,dsp);
[P1_0,P2_0] = Kernel(Mach,k,x,y,z,e_0,tanLambda,relGamma,method,dsp);
[P1_ne,P2_ne] = Kernel(Mach,k,x,y,z,-e_m,tanLambda,relGamma,method,dsp);

%% Preallocation

a = y.^2+z.^2-e_m.^2;
ratio = 2*e_m.*abs(z)./a;
d1 = zeros(size(a));
d2 = zeros(size(a));
Fp = zeros(size(a));
D2_cs = zeros(size(a));
chord=repmat(Panel.chord',length(Panel.chord'),1);

ind1=a>=10e-6;
ind2=and(a<10e-6,a>-10e-6);
ind3=a<=-10e-6;

d1(ind1) = 1;
d1(ind3) = 1;
d1(ind2) = 0;

d2(ind1) = 0;
d2(ind3) = 1;
d2(ind2) = 0.5;

%% Approximation

%Planar case
i0 = abs(z)./e_m <= 0.001;
Fp_i0 = 2*e_m./(y.^2-e_m.^2);

%Near-planar case
ia = ratio <= 0.3 & abs(z)./e_m > 0.001;
s=0;
for n=2:7
    s_step=(-1).^n./(2.*n-1).*ratio.^(2.*n-4);
    s=s+s_step;
end
eps_ia = (4.*e_m.^4./a.^2).*s;
Fp_ia = (d1.*2.*e_m./a).*(1-eps_ia.*(z.^2./e_m.^2))+d2.*(pi./z);

%rest cases
ir = ratio > 0.3 & abs(z)./e_m > 0.001;
eps_ir = (e_m.^2./z.^2).*(1-(a./(2.*e_m.*abs(z)))).*atan(2.*e_m.*abs(z)./a);
Fp_ir = (d1.*2.*e_m./a).*(1-eps_ir.*(z.^2./e_m.^2))+d2.*(pi./z);

Fp(i0)=Fp_i0(i0);
Fp(ia)=Fp_ia(ia);
Fp(ir)=Fp_ir(ir);

%% D1

A1 = (P1_ne-2*P1_0+P1_e)./(2.*e_m.^2);
B1 = (P1_e-P1_ne)./(2.*e_m);
C1 = P1_0;

D1_cs = chord.*(((y.^2-z.^2).*A1+y.*B1+C1).*Fp+(0.5*B1+y.*A1).*log(((y-e_m).^2+z.^2)./((y+e_m).^2+z.^2))+2*e_m.*A1)/(8*pi);

%% Plot Kernel function
if pkern
    eta = -0.5:0.01:0.5;
    eta2 = eta.^2;
    
    y_ydel = 0:0.01:1;
    
    K_total = A1(1,1)*eta2 + B1(1,1)*eta + C1(1,1);
    
    figure(1);
    clf;
%     ax1 = subplot(1,2,1); hold on;
    plot(y_ydel,real(K_total),'b-')
%     title('Real part of Kernel function')
%     xlabel('Dimensionless Panel span')
%     ylabel('Planar Kernel Increment')
    title('Realteil der Kernelfunktion')
    xlabel('Dimensionslose Elementspannweite')
    ylabel('Planarer Inkrement der Kernelfunktion')

    figure(2);
    clf;
%     ax2 = subplot(1,2,2); hold on;
    plot(y_ydel,imag(K_total),'b-')
%     title('Imaginary part of Kernel function')
%     xlabel('Dimensionless Panel span')
%     ylabel('Planar Kernel Increment')
    title('Imaginärteil der Kernelfunktion')
    xlabel('Dimensionslose Elementspannweite')
    ylabel('Planarer Inkrement der Kernelfunktion')
end

%% D2

A2 = (P2_ne-2*P2_0+P1_e)./(2.*e_m.^2);
B2 = (P2_e-P2_ne)./(2.*e_m.^2);
C2 = P2_0;

ib = abs(1./ratio) <= 0.1 & abs(z)./e_m > 0.001;
ic = abs(1./ratio) > 0.1 & abs(z)./e_m > 0.001;

D2_cs_ib = (chord./(16.*pi.*z.^2)).*(((y.^2+z.^2).*A2+y.*B2+C2).*Fp+(1./((y+e_m).^2+z.^2)).*(((y.^2+z.^2).*y+(y.^2-z.^2).*e_m).*A2+(y.^2+z.^2+y.*e_m).*B2+(y+e_m).*C2)...
            -(1./((y-e_m).^2.*z.^2)).*(((y.^2+z.^2).*y-(y.^2-z.^2).*e_m).*A2+(y.^2+z.^2-y.*e_m).*B2+(y-e_m).*C2));

D2_cs_ic = (chord.*e_m./(8.*pi.*a)).*(((2.*(y.^2+z.^2+e_m.^2).*(e_m.^2.*A2+C2)+4.*y.*e_m.^2.*B2)./(((y+e_m).^2+z.^2).*((y-e_m).^2+z.^2)))-(eps_ia./e_m.^2).*((y.^2+z.^2).*A2+y.*B2+C2));

D2_cs(ib) = D2_cs_ib(ib);
D2_cs(ic) = D2_cs_ic(ic);

%% Addition of D's

D_dlm = D1_cs;%+D2_cs;

end

