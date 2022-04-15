classdef DLMpro < handle
    %% Description
    %   This class set up a Doublet-Lattice-Model for a given Wing
    %   Several Methods are provided for calculating the Aerodyanmic
    %   Influence Coefficients
    %
    %   Constructor:
    %       self = DLMPro(Span,Chord,NS,NC,Sweep,Dihedral,TR)
    %
    %   Methods:
    %       res = self.calcAIC(kVect,machVect)
    %
    %   This Tool will be updated continualy with further modules. (e.g. to
    %   visualize the results, etc.)
    %
    %
    %   For more detailed Information on methods see Documentation 
    %       -> "supplemental Software"
    %
    %% Info for Output
    %   The created object contains a property ("resultOverview"). This
    %   gives an short overview of the conducted analyses. Since several
    %   method are available, by reanalysing the AIC, a new result
    %   structure will be created. This allows to compare different
    %   calculated AIC for the same configuration.
    %
    %   Output Coordinate Systems:
    %       Qkk - in aero       Coordinate System 
    %       Qss - in structural Coordinate System (after Transformation with "CSChange")
    %
    %% References
    %   [1]     - Voß, A., “An Implementation ft he Vortex Lattice and the Doublet Lat-
    %             tice Method", Institut für Aeroelastik, Deutsches Zentrum für Luft- und 
    %             Raumfahrt, Göttingen, Oktober 2020. 
    %   [2]     - Albano, E. und Rodden, W. P., “A Doublet Lattice Method For Calculat-
    %             ing Lift Distributions on Oscillation Surfaces in Subsonic Flows," in AIAA 
    %             6th Aerospace Sciences Meeting, New York, 1968. 
    %   [3]     - Kotikalpudi, A., "Body Freedom Flutter (BFF) Doublet Lattice Method 
    %             (DLM)," University of Minnesota Digital Conservancy, 09-Sep-2014. 
    %             [Online]. Verfügbar: http://hdl.handle.net/11299/165566.
    %   [4]     - Kotikalpudi, A., Pfifer, H., und Balas, G. J., "Unsteady Aerodynamics
    %             Modeling for a Flexible Unmanned Air Vehicle," in AIAA Atmospheric
    %             Flight Mechanics Conference, Dallas, Texas, 2015.
    %
    %% Examples
    %   Examples can be found in the "Example_Folder"
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
        
    %% Properties
    properties (Constant)
        version = 1.0.0   % Version of Code
    end
    properties
        wingProp        % wing properties
        panelProp       % panel properies
        SYM             % symmetry flag (1: symmetry, 2: no symmetry)
        geo             % geometry for visualization
        cT              % parameter for coordinate System transformation
        resultArray={}  % result array
    end
    properties (Dependent)
        resultOverview  % result overview
        T_CS            % Transform_CS
    end

    %% Normal Methods
    methods
        % Constructor
        function self = DLMpro(Span,Chord,NS,NC,Sweep,Dihedral,TR,varargin)
            % creates DLMpro Object
            %% Syntax
            %   obj = DLMPro(Span,Chord,NS,NC,Sweep,Dihedral,TR)
            %   obj = DLMPro(__,varargin)
            %
            %% Input
            %   Mandatory:
            %       Span        [double]    - Span of Wing
            %       Chord       [double]    - Chord of Wing @ Root
            %       NS          [double]    - Number of Panels in Span Directions
            %       NC          [double]    - Number of Panels in Chord Directions
            %       Sweep       [double]    - Sweep of Wing (Measured at LE) [in degree]
            %       Dihedral    [double]    - Dihedral of Wing [in degree]
            %       TR          [double]    - Taper Ratio of Wing
            %
            %   Optional (Name-Value Pair)
            %       shift       [double]    - absolut shift of cs [x,y,z]
            %       kVect       [double]    - Vector containing reduced Freq.
            %       MachVect    [double]    - Vector containing mach number
            %       GSA         [double]    - Transformation matrix between Structure and AeroGrid
            %       Integration [string]    - Integration method: choose between
            %                                   ["Parabolic" (default), "Quartic"]
            %       Approximation [string]  - Approximation method: choose between
            %                                   ["Watkins","Laschka"(default),"Desmarais"]
            %
            %       SYM         [boolean]   - is wing symmetric - boolean
            %       CSChange    [double]    - coordinate system rotation describing change of
            %                                 coordinate system between structural
            %                                 and aerodynamic [0,0,0]
            %                                 (default) - [in degrees]
            %                               - This is purely experimently - use with caution!!
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

            %% Add Path
            if(~exist('DLM.m','file')||~exist('VLM.m','file'))
                addpath('.\DLM')
                addpath('.\VLM')
            end
            
            %% Input Parser
            p = inputParser;
            addOptional(p,'shift',[0,0,0]);
            addOptional(p,'kVect',[]);
            addOptional(p,'MachVect',0);
            addOptional(p,'GSA',[]);
            addOptional(p,'Integration',"Parabolic",@(s) any(strcmpi(s,["Parabolic","Quartic"])));
            addOptional(p,'Approximation',"Laschka",@(s) any(strcmpi(s,["Watkins","Laschka","Desmarais"])));
            addOptional(p,'SYM',1);
            addOptional(p,'CSChange',[])
            p.parse(varargin{:});

            shift    = p.Results.shift;
            kVect    = p.Results.kVect;
            machVect = p.Results.MachVect;
            gsa     = p.Results.GSA;
            int      = p.Results.Integration;
            app      = p.Results.Approximation;
            sym      = p.Results.SYM;
            ct       = p.Results.CSChange;

            if isempty(ct)
                ct = [0,0,0];
            end
            self.cT = ct;

            %% Parameter definition
            % Create Wing Object
            wg = Wing(shift(1),shift(2),shift(3));

            wg.Span = Span;         % Wing Span
            wg.Chord = Chord;       % wing chord
            wg.NS = NS;             % Number of spanwise boxes
            wg.NC = NC;             % Number of chordwise boxes
            wg.Sweep  = Sweep;      % sweep in degrees
            wg.Dihedral  = Dihedral;% dihedral in degrees
            wg.TR = TR;             % Taper ratio
            self.wingProp  = wg;    % add object to DLMpro Object

            % check for symmetry
            if sym == 1
                symPoints = wg.Symmetry; % points of wing with symmetry are calculated
                symNS = NS*2;            % double the number of spanwise panel in case of symmetry
            else
                symPoints = wg.Points;
                symNS = NS;
            end
            
            self.wingProp  = wg;
            self.panelProp = Panel(symPoints, wg.NC, symNS);
            self.SYM = sym;
            self.createGeo;
            
            % CalcAIC if direct requested
            if(~isempty(kVect))
                self.calcAIC(kVect,machVect,int,app,gsa)
            end         
        end
       
        function res = calcAIC(self,kVect,machVect,varargin)
            %  self.calcAIC function
            %  calculates AIC for given reduced Frequencies and Mach
            %  Numbers
            %% Syntax
            %   res = self.calcAIC(kVect,machVect)
            %   res = self.calcAIC(__,varargin)
            %% Input
            %   Mandatory:
            %       kVect       [double]    - Vector containing reduced Freq.
            %       MachVect    [double]    - Vector containing mach number
            %
            %   Optional (Name-Value-Pair)
            %       
            %       Integration [string]    - Integration method: choose between
            %                                   ["Parabolic" (default), "Quartic"]
            %       Approximation [string]  - Approximation method: choose between
            %                                   ["Watkins","Laschka"(default),"Desmarais"]
            %       GSA         [double]    - Transformation matrix between Structure and AeroGrid 
            %                               - only for describing the splines
            %                               - orientation of coordinate system will not be taken into account in this matrix 
            %                               - use CSChange property
            %       plotKernel  [boolean]   - plots the kernel function 
            %       verbose     [boolean]   - give detailed outout in console
            %
            %% Output
            %   The output will be stored in a result structure and also be
            %   added to the object's property obj.resultArray
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
            
            %% Input Parser
            p = inputParser;
            addOptional(p,'Integration',"Parabolic",@(s) any(strcmpi(s,["Parabolic","Quartic"])));
            addOptional(p,'Approximation',"Laschka",@(s) any(strcmpi(s,["Watkins","Laschka","Desmarais"])));
            addOptional(p,'GSA',[]);
            addOptional(p,'plotKernel',0)
            addOptional(p,'verbose',0)
            p.parse(varargin{:});

            gsa      = p.Results.GSA;
            int      = p.Results.Integration;
            app      = p.Results.Approximation;
            pkern    = p.Results.plotKernel;
            dsp      = p.Results.verbose;
            
            %% Initialize Variables
            Ref_chord = self.wingProp.Ref_chord;  % reference chord of wing
            totalNumberCases = length(kVect)*length(machVect);
            if(isempty(gsa))
                ndofAero = (self.wingProp.NC*self.wingProp.NS)*2;
                gsa = eye(ndofAero);
                Qss = zeros(ndofAero,ndofAero,totalNumberCases);
                Qkk = zeros(ndofAero,ndofAero,totalNumberCases);
            else
                ndofAero = (self.wingProp.NC*self.wingProp.NS)*2;
                [sn,an] = size(gsa);
                Qss = zeros(sn,sn,totalNumberCases);
                Qkk = zeros(an,an,totalNumberCases);
            end
                        
            %% Calc Values
            disp("--- START Calculation ---");
            disp("   --- Integration: "+int);
            disp("   --- Approximation: "+app);
            pa = self.panelProp;
            combi = zeros(totalNumberCases,2);
            AIC = zeros(self.wingProp.NC*self.wingProp.NS,self.wingProp.NC*self.wingProp.NS,totalNumberCases);
            n = 0;
            for ma = 1:length(machVect)
                for ki=1:length(kVect)
                    n = n+1;
                    mach = machVect(ma);
                    kr = kVect(ki);
                    combi(n,1) = kr;
                    combi(n,2) = mach;
                    
                    kred = kr*2/Ref_chord; %form of the reduced frequency used for the DLM
                    D_vlm = VLM(pa,mach,dsp); %Get VLM downwash effect (steady state effect)
                    D_dlm = DLM(pa,mach,kred,int,app,pkern,dsp); %Get DLM downwash effect (oscillatory effect)

                    % combine VLM and DLM matrices
                    D = -(D_vlm+D_dlm);

                    % concatenate results to one wing
                    if self.SYM == 1
                        NS_hw = self.wingProp.NS;
                        NC = self.wingProp.NC;
                        D_r = D(NS_hw*NC+1:end,NS_hw*NC+1:end);   % influence from right wing on itself
                        D_lr = D(NS_hw*NC+1:end,1:NS_hw*NC);      % influence from left wing on right wing
                        sub_D_lr = zeros(NS_hw*NC,NC,n);

                        %sorting the influence matrix from left wing on
                        %right wing so that it can easily summed up with
                        %influnece from right wing on itself
                        m=1;
                        l = NS_hw;
                        for k=1:NS_hw
                            sub_D_lr(:,:,m) = D_lr(:,(k-1)*NC+1:k*NC);
                            m = m+1;
                        end
                        
                        for k=1:NS_hw
                            new_D_lr(:,(k-1)*NC+1:k*NC) = sub_D_lr(:,:,l);
                            l = l - 1;
                        end
                        
                        D = D_r+new_D_lr;
                    end
                    AIC(:,:,n) = inv(D);

                    %Calculate Derivation matrix W and Integration matrix B
                    [W,B] = DerivationIntegrationMatrix(kr,Ref_chord,pa);
                    if self.SYM == 1
                        W = W(ndofAero/2+1:end,ndofAero+1:end);
                        B = B(ndofAero+1:end,ndofAero/2+1:end);
                    end

                    W_total(:,:,n) = W;
                    % Calc QKK and Qhh
                    Qkk(:,:,n) = B*(D\W); % B*inv(D)*W;

                    % create coordinate Transform matrix to account 
                    % for different CS-Orientation in DLMpro(aero) and Struture

                    T = [self.T_CS self.T_CS; self.T_CS self.T_CS];
                    T([1,2,4,6],:) = [];
                    T(:,[1,2,4,6]) = [];
                    
                    TR = zeros(size(gsa,2)/2);
                    for i = 1:size(TR,1)
                        s = 2*(i-1)+1;
                        TR(s:s+1,s:s+1) = T;
                        i=i+1;
                    end

                    QkkT = TR*Qkk(:,:,n)*TR;
                    Qss(:,:,n) = gsa*QkkT*gsa.';
                end
            end
            
            res.Nr              = length(self.resultArray)+1;
            res.kMList          = combi;
            res.Qss             = Qss;
            res.Qkk             = Qkk;
            res.integration     = int;
            res.approximation   = app;
            res.VLM             = D_vlm;
            res.AIC             = AIC;
            res.GSA             = gsa;
            res.B               = B;
            res.W               = W_total;

            self.resultArray{end+1,1} = res;
            disp("--- Calculation Finished ---"+newline);
        end
        
        function createGeo(obj)
            % creates geometry for visualization
            %% Note:
            % The Visualization module is not published, yet!
            %

            %% check for Sym
            if(obj.SYM)
                wg = obj.wingProp;
                Pts = obj.panelProp.Pts;
                RP  = obj.panelProp.RP;
                CP  = obj.panelProp.CP;
                SP  = obj.panelProp.SP;
                relNr = (wg.NS+1)*(wg.NC+1)-wg.NC;
                pa.Pts = Pts(relNr:end,:);
                relNr = wg.NS*wg.NC+1;
                pa.RP = RP(relNr:end,:);
                pa.SP = SP(relNr:end,:);
                pa.CP = CP(relNr:end,:);
            else
                pa = obj.panelProp;
                wg = obj.wingProp;
            end
            %% Code
%             pa = obj.panelProp;
%             wg = obj.wingProp;
            
            numberGridpoints = size(pa.Pts,1)+size(pa.RP,1)+size(pa.CP,1)+size(pa.SP,1);
            id = (1:numberGridpoints)';
            startRP = size(pa.Pts,1)+1;
            endRP   = size(pa.Pts,1)+size(pa.RP,1);
            startCP = endRP+1;
            endCP   = size(pa.Pts,1)+size(pa.RP,1)+size(pa.CP,1);
            startSP = endCP+1;
            endSP   = size(pa.Pts,1)+size(pa.RP,1)+size(pa.CP,1)+size(pa.SP,1);
            
            coord = [pa.Pts;pa.RP;pa.CP;pa.SP];
            % grid Points
            gridTable = table(id,coord);
            gridTable.Properties.VariableNames = {'ID','coord'};
            
            % Elements Panel
            numberElem = wg.NC*wg.NS+size(pa.RP,1)+size(pa.CP,1)+size(pa.SP,1);
            
            gridID_AERO = zeros(wg.NC*wg.NS,4);
            ni = 1; k = 1; c = 0;
            pX = wg.NC;  pY = wg.NS;
            for i=1:(pX)*(pY)
                if(c==(pX))
                    c=0;
                    k = k + 1;
                end
                gridID_AERO(ni,1) = k ;
                gridID_AERO(ni,2) = k+1 ;
                gridID_AERO(ni,3) = k+1+(pX+1);
                gridID_AERO(ni,4) = k+pX+1;
                
                k = k + 1;
                ni = ni+1;
                c = c+1;
            end
            type_AERO = repmat("AERO",wg.NC*wg.NS,1);
            
            % Elements RP
            gridID_RP = zeros(size(pa.RP,1),4);
            gridID_RP(:,1) = (startRP:endRP)';
            type_RP   = repmat("RP",size(pa.RP,1),1);
            
            % Elemetns CP
            gridID_CP = zeros(size(pa.CP,1),4);
            gridID_CP(:,1) = (startCP:endCP)';
            type_CP   = repmat("CP",size(pa.CP,1),1);
            
            % Elemetns SP
            gridID_SP = zeros(size(pa.SP,1),4);
            gridID_SP(:,1) = (startSP:endSP)';
            type_SP   = repmat("SP",size(pa.SP,1),1);
            
            % Elements total
            idELEM = (1:numberElem)';
            type = [type_AERO;type_RP;type_CP;type_SP];
            gridID = [gridID_AERO;gridID_RP;gridID_CP;gridID_SP];
            
            elemTable = table(type,idELEM,gridID);
            elemTable.Properties.VariableNames = {'type','idELEM','gridID'};
            
            % create object and store data
            try
                obj.geo = sdb_geometry3D;
                obj.geo.name = "AERO";
                obj.geo.viewSettings.viewGeo = [3.999198187806840e+04,24.40483056985801];
                obj.geo.addDesign('AERO','o',0.5,'k',2,	'#4DBEEE',1,'jet');
                obj.geo.addDesign('RP','o',0.5,'b',2.5,	'#4DBEEE');
                obj.geo.addDesign('CP','o',0.5,'g',2.5,	'#4DBEEE');
                obj.geo.addDesign('SP','o',0.5,'r',2.5,	'#4DBEEE');
            catch
                obj.setWarning(1);
            end

            % change coordinates
            for i = 1:size(gridTable,1)
                currentCoord = gridTable{i,'coord'};
                newCoord = obj.T_CS*currentCoord';
                gridTable{i,'coord'} = newCoord';
            end

            obj.geo.gridTable = gridTable;
            obj.geo.elementTable = elemTable;
            obj.geo.constraintMatrix = 1;
        end

        function plotDeformedMesh(obj,ax,disp)
            % plot a deformed mesh
            %% Note:
            % The Visualization module is not published, yet!
            %

            if(~isa(obj.geo,'sdb_geometry3D'))
                obj.setWarning(1);
                return;
            end

            % plotDeformed Mesh
            %% 
            %
            x = obj.panelProp.RP(:,1);
            y = obj.panelProp.RP(:,2);
            z = disp(:,3);
            
            coord = obj.geo.gridTable{:,3};
            pts_x = coord(:,1);
            pts_y = coord(:,2);
            dispNew = griddata(x,y,z,pts_x,pts_y,'v4');
            obj.geo.plotGeo('ax',ax,'disp',dispNew);   
        end
        
        function plotPressCoeff(obj,alpha,res,n,s,varargin)
            % plots the press Coefficient
            %% Note:
            % The Visualization module is not published, yet!
            %

            if(~isa(obj.geo,'sdb_geometry3D'))
                obj.setWarning(1);
                return;
            end
            % plotPressCoeff

            p = inputParser;
            addOptional(p,'axC',[]);
            addOptional(p,'axSpan',[]);
            addOptional(p,'axStrip',[]);
            p.parse(varargin{:});
            axC = p.Results.axC;
            axS = p.Results.axSpan;
            axSt = p.Results.axStrip;
                        
            
            NC = obj.wingProp.NC;
            NS = obj.wingProp.NS;
            aVect(1:NC*NS) = deg2rad(alpha);
            cp = real(res.AIC(:,:,n))*aVect'; %pressure coefficient as vector
            cp_m = reshape(cp,[NC,NS]); %shape to matrix like panel on wing
            cp_s = cp_m(:,s);
            cp_s2 = cp_m(:,NS);
%             if obj.SYM ==1
%                 PArea = obj.panelProp.A;
%                 PArea(1:NC*NS)=[];
%             else
%                 PArea = obj.panelProp.A;
%             end

            ca = mean(cp_m);    %calculated lift coefficient of each strip
            ca_w = mean(ca);
            
            %% plot
            if(isempty(axC))
                figure(200); clf;
                axC = nexttile;
                colorbar
                axC.Title.String = "Panel pressure coefficients";
            end
            obj.geo.plotGeo('ax',axC,'colorData',cp);

            if(isempty(axSt))
                figure(100); clf;
                axSt = nexttile;
            end
            
            x=(1/NC)*0.25:1/NC:1;
            x(end+1)=1;
            cp_s(end+1)=0;
            cp_s2(end+1)=0;
            ix =(1/NC)*0.25:0.01:1;
            iy = interp1(x,cp_s,ix,'spline');
            iy2 = interp1(x,cp_s2,ix,'spline');
            hold on
            plot(axSt,ix,iy,'b',ix,iy2,'r')
            legend('Wurzel','Spitze','AutoUpdate','off')
            plot(axSt,x(1:end-1),cp_s(1:end-1),'bo',x(1:end-1),cp_s2(1:end-1),'ro')
            axSt.XLabel.String = "Dimensionslose Tiefe";
            axSt.YLabel.String = "Cp,res";
            axSt.Title.String = "Druckverteilung über die Flügeltiefe";
            
            if(isempty(axS))
                figure(101); clf;
                axS = nexttile;
            end
            plot(axS,0:1/(NS-1):1,ca)
            axS.XLabel.String = "Dimensionless Span";
            axS.YLabel.String = "Ca";
            axS.Title.String = "Lift distribution over span";
            
            x = (1/NS)*0.5:(1/NS):1;
            x(end+1)=1;
            ca(end+1)=0;
            ix = 0:0.001:1;
            iy = interp1(x,ca,ix);
            plot(axS,ix,iy)
            axS.XLabel.String = "Dimensionslose Spannweite";
            axS.YLabel.String = "Cp,res";
            axS.Title.String = "Druckverteilung über die Spannweite";
            
            if(isempty(axC))
                figure(200); clf;
                axC = nexttile;
                colorbar
                axC.Title.String = "Druckbeiwerte der Elemente";
            end

            obj.geo.plotGeo('ax',axC,'xRev',0,'colorData',cp);
        end

        function createGeoT(obj,a,b,c)
            % transforms coordinates into new coordinate System defines by rotations parameter a,b,c
            % only affects the geo-object and thus the visualization
            %% Note:
            % The Visualization module is not published, yet!
            %

            T = asi_csRot(a,b,c);
            obj.csT = [a,b,c];
            for i = 1:obj.geo.noNodes
                currentCoord = obj.geo.gridTable{i,3};
                newCoord = T*currentCoord';
                newCoord(abs(newCoord)<10E-8) = 0;
                obj.geo.gridTable{i,3} = newCoord';
            end
        end

        function plotPhaseLag(obj,alpha,res,n,strip)
            % plots the phase lag
            %% Note:
            % The Visualization module is not published, yet!
            %
            if(~isa(obj.geo,'sdb_geometry3D'))
                obj.setWarning(1);
                return;
            end
            % plotPhaseLag

            t = 0:0.01:10;
            ar = deg2rad(alpha);
            a = ar*exp(1j*2.*t);
            NC = obj.wingProp.NC;
            NS = obj.wingProp.NS;
            cp = sum(res.AIC(:,:,n),2); %pressure coefficient as vector
            cp_m = reshape(cp,[NC,NS]); %shape to matrix like panel on wing
            ca = mean(cp_m);    %calculated lift coefficient of each strip
            ca_st = ca(strip);
            q = ca_st.*a;

            D_vlm = symmetry(obj,res.VLM);
            AIC_vlm = inv(D_vlm);
            cpr = sum(AIC_vlm,2); %pressure coefficient as vector
            cp_mr = reshape(cpr,[NC,NS]); %shape to matrix like panel on wing
            car = mean(cp_mr);
            car_st = car(strip);
            qs = abs(car_st)*ar;
            qs1 = qs*exp(1j*2.*t);

            figure(500); hold on;
            plot(t,real(q))
            plot(t,real(qs1))
            xlabel("time")
            ylabel("cp")
            legend("unsteady","quasi-steady")
        end
        
    end
    
    %% Get Methods
    methods
        function ovTable = get.resultOverview(obj)
            % creates Table to have an overview over all computations made 
            varNames = {'Result No.', 'k-M-List','AIC','Qss','Qkk','Integration','Approximation','GSA'};
            
            nResults = length(obj.resultArray);
            % Initialize Vectors for table
            ResultNo = 1:nResults;
            kMList   = {1:nResults};
            Qss      = {1:nResults};
            Qkk      = {1:nResults};
            Int      = repmat("",nResults,1);
            App      = repmat("",nResults,1);
            
            for rr = 1:nResults
                ResultNo(rr) = obj.resultArray{rr}.Nr;
                kMList{rr}   = obj.resultArray{rr}.kMList;
                AIC{rr}      = obj.resultArray{rr}.AIC;
                Qss{rr}      = obj.resultArray{rr}.Qss;
                Qkk{rr}      = obj.resultArray{rr}.Qkk;
                Int(rr)      = obj.resultArray{rr}.integration;
                App(rr)      = obj.resultArray{rr}.approximation;
                GSA{rr}      = obj.resultArray{rr}.GSA;
            end
            
            ovTable = table(ResultNo', kMList',AIC', Qss', Qkk', Int, App, GSA', 'VariableNames', varNames);
        end

        function T = get.T_CS(obj)
            % if available, use extra function
            try
                T = asi_csRot(obj.cT(1),obj.cT(2),obj.cT(3));
            catch
                a = deg2rad(obj.cT(1));
                b = deg2rad(obj.cT(2));
                c = deg2rad(obj.cT(3));
                
                %% Transformation
                T(1,1) = cos(c)*cos(a)-sin(c)*cos(b)*sin(a);
                T(1,2) = cos(c)*sin(a)+sin(c)*cos(b)*cos(a);
                T(1,2) = sin(c)*sin(b);
                T(2,1) = -sin(c)*cos(a)-cos(c)*cos(b)*sin(a);
                T(2,2) = -sin(c)*sin(a)+cos(c)*cos(b)*cos(a);
                T(2,3) = cos(c)*sin(b);
                T(3,1) = sin(b)*sin(a);
                T(3,2) = -sin(b)*cos(a);
                T(3,3) = cos(b);
            end
                T(abs(T)<10E-8) = 0;
        end

    end

    %% Private Methods
    methods (Access=private)
        function setWarning(obj,c)
            switch c
                case 1
                    warning("The visualization module is not published, yet! Coming soon...");
                otherwise
                    warning("An unexpected warning occurs. :(")
            end
        end
    end
end