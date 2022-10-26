classdef Panel < handle
    %% Disclaimer
    %   Copyright (c) 2021 Nils Böhnisch, Marc Bangel.
    %
    %   This file is part of DLMPro.
    %
    %   DLMPro is free software: you can redistribute it and/or modify
    %   it under the terms of the GNU General Public License as published by
    %   the Free Software Foundation, either version 3 of the License, or
    %   any later version. Also see the file "License".
    
    properties
        Pts(:,3)    %All points
        NC double   %Number of chordwise boxes
        NS double   %Number of spanwise boxes
        PtsIdx      %Index of corner points 1 in Pts matrix
        C1(:,3)     %Corner point 1
        C2(:,3)     %Corner point 2
        C3(:,3)     %Corner point 3
        C4(:,3)     %Corner point 4
        CP(:,3)     %Control point on 3/4 line
        SP(:,3)     %Sending point on 1/4 line midspan
        D1(:,3)     %inner border point on doublet line (1/4 line)
        D2(:,3)     %inner midhalfspan point on doublet line
        D4(:,3)     %outer midhalfspan point on doublet line
        D5(:,3)     %outer border point on doublet line
        RP(:,3)     %reference point of panel(center point)
        chord double%average chord
        s double    %span of box
        N(:,3)      %normal vector
        A double    %Area of box
        AR double   %Aspect ratio of box
    end
    
    methods
        function obj = Panel(Points,NC,NS)
            obj.Pts = Points;
            obj.NC = NC;
            obj.NS = NS;

            if any(obj.AR>=3)
                disp("--- Warning in DLMpro: Some Panels have an Aspect Ratio higher 3!");
                disp("   -> it is recommended to remesh the Geometry!");
                disp("   -> and/or use 'Quartic' Integration!");
            end
        end
        function PtsIdx = get.PtsIdx(obj)
            a = 1:length(obj.Pts);
            a(end-(obj.NC+1):end) = [];
            delPts = (1:obj.NS-1).*(obj.NC+1);
            a(delPts) = [];
            PtsIdx = a;
        end
        function C1 = get.C1(obj)
            C1 = obj.Pts(obj.PtsIdx,:);
        end
        function C2=get.C2(obj)
            a = obj.PtsIdx+(obj.NC+1);
            C2 = obj.Pts(a,:);
        end
        function C3=get.C3(obj)
            a = obj.PtsIdx+1;
            C3 = obj.Pts(a,:);
        end
        function C4=get.C4(obj)
            a = obj.PtsIdx+1+(obj.NC+1);
            C4 = obj.Pts(a,:);
        end
        function D5=get.D5(obj)
            D5=obj.C2+(obj.C4-obj.C2)*0.25;
        end
        function D1=get.D1(obj)
            D1=obj.C1+(obj.C3-obj.C1)*0.25;
        end
        function SP=get.SP(obj)
            SP=obj.D1+(obj.D5-obj.D1)*0.5;
        end
        function D2=get.D2(obj)
            D2=obj.D1+(obj.D5-obj.D1)*0.25;
        end
        function D4=get.D4(obj)
            D4=obj.D1+(obj.D5-obj.D1)*0.75;
        end
        function CP=get.CP(obj)
            MP=obj.C3+(obj.C4-obj.C3)*0.5;
            CP=obj.SP+(MP-obj.SP)*(2/3);
        end
        function RP=get.RP(obj)
            MP=obj.C3+(obj.C4-obj.C3)*0.5;
            RP=obj.SP+(MP-obj.SP)*(1/3);
        end
        function chord=get.chord(obj)
            r=obj.C3(:,1)-obj.C1(:,1);
            t=obj.C4(:,1)-obj.C2(:,1);
            chord=0.5*(r+t);
        end
        function s=get.s(obj)
            dist=obj.C4(:,2:3)-obj.C3(:,2:3);
            s=sqrt(dist(:,1).^2+dist(:,2).^2);
        end
        function A=get.A(obj)
            PChord=obj.C3(:,1)-obj.C1(:,1);
            PTip=obj.C4(:,1)-obj.C2(:,1);
            A=0.5*(PTip+PChord).*obj.s;
        end
        function AR=get.AR(obj)
            AR=obj.s.^2./obj.A;
        end
        function N=get.N(obj)
            Np=cross((obj.C3-obj.C1),(obj.C2-obj.C1));
            norm=sqrt(Np(:,1).^2+Np(:,2).^2+Np(:,3).^2);
            N_x=Np(:,1)./norm;
            N_y=Np(:,2)./norm;
            N_z=Np(:,3)./norm;
            N=[N_x,N_y,N_z];
        end
    end
end