function [D_dlm] = QuarticApprox(Panel,Mach,k,x,y,z,e_m,e_0,e_2,tanLambda,relGamma,method,pkern,dsp)
% QuarticApprox calculates the influence matrix if quartic integration
% is selected
%% Description:
%   In the quartic integration the kernel function is calculated five
%   times for five different points on the doublet line. The coefficients
%   for the approximation of the integral of the kernel function are
%   calculated. The approximation describes a quartic equation. At least
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
%   e_2: matrix of semispans of panel half spans
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

[Q1_e,Q2_e] = Kernel(Mach,k,x,y,z,e_m,tanLambda,relGamma,method,dsp);
[Q1_e2,Q2_e2] = Kernel(Mach,k,x,y,z,e_2,tanLambda,relGamma,method,dsp);
[Q1_0,Q2_0] = Kernel(Mach,k,x,y,z,e_0,tanLambda,relGamma,method,dsp);
[Q1_ne2,Q2_ne2] = Kernel(Mach,k,x,y,z,-e_2,tanLambda,relGamma,method,dsp);
[Q1_ne,Q2_ne] = Kernel(Mach,k,x,y,z,-e_m,tanLambda,relGamma,method,dsp);

%% Preallocation

a = y.^2+z.^2-e_m.^2;
ratio = 2*e_m.*abs(z)./a;
d1 = zeros(size(a));
d2 = zeros(size(a));
Fq = zeros(size(a));
D2_cs = zeros(size(a));
chord=repmat(Panel.chord,1,length(Panel.chord)); %matrix of panel chords

ind1=a>=10e-6;
ind2=and(a<10e-6,a>-10e-6);
ind3=a<=-10e-6;

d1(ind1) = 1;
d1(ind3) = 1;
d1(ind2) = 0;

d2(ind1) = 0;
d2(ind3) = 1;
d2(ind2) = 0.5;

%% Calculations of cases

%Planar case
i0 = abs(z)./e_m <= 0.001;
Fq_i0 = d1.*2.*e_m./(y.^2-e_m.^2);

%Near-planar case
ia = ratio <= 0.3 & abs(z)./e_m > 0.001;
s=0;
for n=2:7
    s_step=(-1).^n./(2.*n-1).*ratio.^(2.*n-4);
    s=s+s_step;
end
eps_ia = (4.*e_m.^4./a.^2).*s;
Fq_ia = (d1.*2.*e_m./a).*(1-eps_ia.*(z.^2./e_m.^2))+d2.*(pi./z);

%rest cases
ir = ratio > 0.3 & abs(z)./e_m > 0.001;
eps_ir = (e_m.^2./z.^2).*(1-(a./(2.*e_m.*abs(z)))).*atan(2.*e_m.*abs(z)./a);
Fq_ir = (d1.*2.*e_m./a).*(1-eps_ir.*(z.^2./e_m.^2))+d2.*(pi./z);

Fq(i0)=Fq_i0(i0);
Fq(ia)=Fq_ia(ia);
Fq(ir)=Fq_ir(ir);

%% D1

A1 = (-1./(6.*e_m.^2)).*(Q1_ne-16.*Q1_ne2+30.*Q1_0-16.*Q1_e2+Q1_e);
B1 = (1./(6.*e_m)).*(Q1_ne-8.*Q1_ne2+8.*Q1_e2-Q1_e);
C1 = Q1_0;
D1 = (-2./(3.*e_m.^3)).*(Q1_ne-2.*Q1_ne2+2.*Q1_e2-Q1_e);
E1 = (2./(3.*e_m.^4)).*(Q1_ne-4.*Q1_ne2+6.*Q1_0-4.*Q1_e2+Q1_e);

D1_cs = (chord/(8*pi)).*(((y.^2-z.^2).*A1+y.*B1+C1+y.*(y.^2-3.*z.^2).*D1+(y.^4-6.*y.^2.*z.^2+z.^4).*E1).*Fq...
        +(y.*A1+0.5.*B1+0.5.*(3.*y.^2-z.^2).*D1+2.*y.*(y.^2-z.^2).*E1).*log(((y-e_m).^2+z.^2)./((y+e_m).^2+z.^2))...
        +2.*e_m.*(A1+2.*y.*D1+(3.*y.^2-z.^2+(1/3).*e_m.^2).*E1));

%% Plot Kernel function

if pkern
    eta = -0.5:0.01:0.5;
    eta2 = eta.^2;
    eta3 = eta.^3;
    eta4 = eta.^4;
    
    y_ydel = 0:0.01:1;
    
    K_total = A1(1,1)*eta2 + B1(1,1)*eta + C1(1,1) + D1(1,1)*eta3 + E1(1,1)*eta4;
    
    fig = figure(1);
    hold on
    ax1 = fig.Children(1);
%     ax1.Title.String = "Real part of Kernel function";
%     ax1.XLabel.String = "Dimensionless Panel span";
%     ax1.YLabel.String = "Planar Kernel Increment";
    ax1.Title.String = "Realteil der Kernelfunktion";
    ax1.XLabel.String = "Dimensionslose Elementspannweite";
    ax1.YLabel.String = "Planarer Inkrement der Kernelfunktion";
    plot(ax1,y_ydel,real(K_total),'r-')
    legend('Parabolisch','Quartisch')

    fig2 = figure(2);
    hold on
    ax2 = fig2.Children(1);
%     ax2.Title.String = "Imaginary part of Kernel function";
%     ax2.XLabel.String = "Dimensionless Panel span";
%     ax2.YLabel.String = "Planar Kernel Increment";
    ax2.Title.String = "Imaginärteil der Kernelfunktion";
    ax2.XLabel.String = "Dimensionslose Elementspannweite";
    ax2.YLabel.String = "Planarer Inkrement der Kernelfunktion";
    plot(ax2,y_ydel,imag(K_total),'r-')
    legend('Parabolisch','Quartisch')
end

%% D2

A2 = (-1./(6*e_m.^2)).*(Q2_ne-16*Q2_ne2+30*Q2_0-16*Q2_e2+Q2_e);
B2 = (1./(6*e_m)).*(Q2_ne-8*Q2_ne2+8*Q2_e2-Q2_e);
C2 = Q2_0;
D2 = (-2./(3*e_m.^3)).*(Q2_ne-2*Q2_ne2+2*Q2_e2-Q2_e);
E2 = (2./(3*e_m.^4)).*(Q2_ne-4*Q2_ne2+6*Q2_0-4*Q2_e2+Q2_e);

ib = abs(1./ratio) <= 0.1 & abs(z)./e_m > 0.001;
ic = abs(1./ratio) > 0.1 & abs(z)./e_m > 0.001;
del = (e_m./abs(z)).^2.*(1-d1-d2.*(pi./ratio));

D2_cs_ib = (chord./(16.*pi.*z.^2)).*(((y.^2+z.^2).*A2+y.*B2+C2+y.*(y.^2+3.*z.^2).*D2+(y.^4+6.*y.^2.*z.^2-3.*z.^4).*E2).*Fq...
            +(1./((y+e_m).^2+z.^2)).*(((y.^2+z.^2).*y+(y.^2-z.^2).*e_m).*A2+(y.^2+z.^2+y.*e_m).*B2+(y+e_m).*C2+(y.^4-z.^4+(y.^2-3.*z.^2).*y.*e_m).*D2...
            +((y.^4-2.*y.^2.*z.^2-3.*z.^4).*y+(y.^4-6.*y.^2.*z.^2+z.^4).*e_m).*E2)-(1./((y-e_m).^2+z.^2)).*(((y.^2+z.^2).*y-(y.^2-z.^2).*e_m).*A2...
            +(y.^2+z.^2-y.*e_m).*B2+(y-e_m).*C2+(y.^4-z.^4-(y.^2-3.*z.^2).*y.*e_m).*D2+((y.^4-2.*y.^2.*z.^2-3.*z.^4).*y-(y.^4-6.*y.^2.*z.^2+z.^4).*e_m).*E2)...
            +(z.^2.*log(((y-e_m).^2+z.^2)./((y+e_m).^2+z.^2))).*D2+4.*z.^2.*(e_m+y.*log(((y-e_m).^2+z.^2)./((y+e_m).^2+z.^2))).*E2);

D2_cs_ic = (e_m.*chord./8.*pi.*a).*((1./(((y+e_m).^2+z.^2).*((y-e_m).^2+z.^2))).*(2.*(y.^2+z.^2+e_m.^2).*(e_m.^2.*A2+C2)+4.*y.*e_m.^2.*B2...
            +2.*y.*(y.^4-2.*e_m.^2.*y.^2+2.*y.^2.*z.^2+3.*e_m.^4+2.*e_m.^2.*z.^2+z.^4).*D2+2.*(3.*y.^6-7.*e_m.^2.*y.^4+5.*y.^4.*z.^2 ...
            +6.*e_m.^4.*y.^2+6.*e_m.^2.*y.^2.*z.^2-3.*e_m.^2.*z.^4-z.^6+y.^2.*z.^4-2.*e_m.^4.*z.^2).*E2)...
            -((d1.*eps_ir+del)./e_m.^2).*((y.^2+z.^2).*A2+y.*B2+C2+y.*(y.^2+3.*z.^2).*D2+(y.^4+6.*y.^2.*z.^2-3.*z.^4).*E2))...
            +(chord./(8.*pi)).*((D2./2).*log(((y-e_m).^2+z.^2)./((y+e_m).^2+z.^2))+2.*(e_m+y.*log(((y-e_m).^2+z.^2)./((y+e_m).^2+z.^2))).*E2);

D2_cs(ib)=D2_cs_ib(ib);
D2_cs(ic)=D2_cs_ic(ic);

%% Addition of D's

D_dlm = D1_cs+D2_cs;

end

