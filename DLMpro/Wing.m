classdef Wing < handle
    properties
        Span double     %Span of Wing
        Chord double    %Chord of Wing
        Sweep double    %Sweep of Wing in deg
        Dihedral double %Dihedral of Wing in deg
        TR double       %Taper Ratio
        NC double       %Number of chordwise panels
        NS double       %Number of spanwise panels
        x_off double    %offset of x-coordinate
        y_off double    %offset of y-coordinate
        z_off double    %offset of z-coordinate
        CT    double    %CS-Transform with respect to Reference System
    end

    properties(Dependent)
        Semispan double     %Semispan of Wing
        Sweepr double       %Sweep in rad
        Dihedralr double    %Dihedral in rad
        Tiplength double    %Tiplength of Wing
        Area double         %Area of one wing
        AR double           %Aspect Ratio
        Points(:,3)         %Points of the aerodynamic grid
        Ref_chord double    %reference chord
        T_CT double         %Transformmationmatrix CS Transform
    end

    methods
        function obj = Wing(x_off,y_off,z_off)
            obj.x_off=x_off;
            obj.y_off=y_off;
            obj.z_off=z_off;
        end

        function symPoints = Symmetry(obj)
            sym = obj.Points;
            sym(1:obj.NC+1,:) = [];
            sym(:,2) = -sym(:,2);
            sym_f = [];
            for m=1:obj.NS
                sub_sym = sym(1:obj.NC+1,:);
                sym(1:obj.NC+1,:) = [];
                sym_f = [sub_sym;sym_f];
            end
            symPoints = [sym_f;obj.Points];
        end
    end

    % get-Methods
    methods
        function T = get.T_CT(obj)
            T = asi_csRot(obj.CT(1),obj.CT(2),obj.CT(3));
            T(abs(T)<10E-8) = 0;
        end

        function Sweepr=get.Sweepr(obj)
            Sweepr=deg2rad(obj.Sweep);    %converts sweep from deg in rad
        end

        function Dihedralr=get.Dihedralr(obj)
            Dihedralr=deg2rad(obj.Dihedral);    %converts dihedral from deg in rad
        end
        
        function Semispan=get.Semispan(obj)
            Semispan=obj.Span/2;   %calculates semispan
        end

        function Tiplength=get.Tiplength(obj)
            Tiplength=obj.Chord*obj.TR;     %calculates tiplength dependent on taper ratio and chord
        end

        function Ref_chord=get.Ref_chord(obj)
            Ref_chord = (obj.Chord+obj.Tiplength)/2;
        end

        function Area=get.Area(obj)
            Area=0.5*(obj.Tiplength+obj.Chord)*(obj.Semispan/cos(obj.Dihedralr)); %calculates area of one wing
        end

        function AR=get.AR(obj)
            AR=obj.Span^2/(2*obj.Area);  %calculates aspect ratio with overall wing area (2 wings)
        end

        function Points=get.Points(obj)
            TE_ang = atan((obj.Semispan*tan(obj.Sweepr)+obj.Tiplength-obj.Chord)/obj.Semispan); %angle of trailing edge
            i = 1;  %initialize point number
            y = 0;  %starting y coordinate
            z = 0;  %starting z coordinate
            dy = obj.Semispan/obj.NS;
            GridPoints=zeros((obj.NC+1)*(obj.NS+1),3);  %preallocation of point matrix for performance
            for k=1:obj.NS+1
                x = y*tan(obj.Sweepr);  %starting x coordinate
                segLength = obj.Chord + y*tan(TE_ang) - x;  %length of chord depending on y- coordinate
                dx = segLength/obj.NC;  %difference of x points
                for m=1:obj.NC+1
                    GridPoints(i,1:3) = [x, y,  z+tan(obj.Dihedralr)*y];   %points on wing for grid generation
                    x = x+dx;    %increment y coordinate for next point
                    i = i+1;    %increment point number
                end
                y = y+dy;   %increment x coordinate for next point row
            end
            Points=GridPoints;
            
            %add offset to points
            offset = [obj.x_off;obj.y_off;obj.z_off];
            Points(:,1) = Points(:,1)+offset(1);
            Points(:,2) = Points(:,2)+offset(2);
            Points(:,3) = Points(:,3)+offset(3);
        end
    end

end