function [Force, Torque] = recoverForces(Force, Torque, Markers)
    %RECOVERFORCES: Returns forces and torques in Vicon Global Coordinates from Force/Torque Sensor Readings
    %
    %   [Force, Torque] = recoverForces(Force, Torque, Markers) takes the
    %   raw Force and Torque measurements (in V) from the Nano25 F/T Sensor
    %   and converts them into Force (N) and Torque (Nm) measurements in
    %   the Vicon Global Coordinate System. recoverForces converts V to N
    %   and Nm using predefined constants. To convert from the force sensor
    %   coordinates to Vicon Coordinates, Markers must be specified. 
    %   
    %   The input syntax for recoverForces is:
    %       Force: a 3x1 vector of force measurements in (V)
    %       Torque: a 3x1 vector of torque (or moment) measurements in V
    %       Markers: a structure of force handle markers containing the
    %       following fields:
    %           frontleft: 3x1 vector
    %           frontright: 3x1 vector
    %           frontmiddle: 3x1 vector
    %           backmiddle: 3x1 vector
    %
    %   The outputs from recoverForces are:
    %       Force: a 3x1 vector of forces (N) in the Vicon Coordinate
    %       System
    %       Torque: a 3x1 vector of torques (Nm) in the Vicon Coordinate
    %       System
    %
    %   Note: if Markers is not specified, recoverForces will return the
    %   Forces and Torques in the proper units, but the vectors will be
    %   expressed in the sensor coordinate system.
    
    %   Luke Drnach
    %   October 27, 2017
    %
    %   Updated December 4, 2018. Update allows RECOVERFORCES to operate on
    %   arrays of forces and marker data, instead of solely on point data.
    
    % Scale Factors (Convert V to N or Nm)
    F_Scale = [12.5, 12.5, 50]';     %Coordinates are scaling factors in [X, Y, Z]
    T_Scale = [0.3, 0.3, 0.3]';
    % Sensor Biases (Sensor offsets in V)
%     F_Bias = [0 0 0]'; %Force bias in [X, Y, Z] (V) - set to zeros for test
    F_Bias = [-0.2567, -0.1931, -0.1590]'; %Force bias in [X, Y, Z] (V)
%     coordinates % Original code uses one set of values for all subjects
%     and not sure how they got these values
%     F_Bias = [-0.2002    0.0455   -0.2017]'; %Force bias in [X, Y, Z] (V) from mean of Assist Solo trials for HHI08
    T_Bias = [-0.6627, 0.6189, 0.7339]'; %Torque bias
    % Convert sensor output (V) to forces (N) and torques (Nm) assuming a
    % linear relationship: (Force/Torque) = (Scale)[(Sensor Reading) -
    % (Bias)]. This conversion is performed in the Force/Torque Sensor's
    % Reference frame
    Force = F_Scale.*(Force - F_Bias);
    Torque = T_Scale.*(Torque - T_Bias);
       
    if nargin == 3
        N = size(Force,2);
        for n = 1:N
            % Convert the [Force, Torque] from the sensor's reference frame to the
            % global reference frame using the markers
            
            % First, check that all the markers are reasonably where they should be.
            FR = Markers.frontright(:,n);
            FL = Markers.frontleft(:,n);
            FM = Markers.frontmiddle(:,n);
            BM = Markers.backmiddle(:,n);
            % Calculate distances between the markers
            tolerance = 150; % (mm) this is the only part of this code that cares about units of Markers
            dRL = norm(FR-FL);  %Distance between Front Right and Front Left
            dRM = norm(FR-FM);  %Distance between Front Right and Front Middle
            dRB = norm(FR-BM);  %Distance between Front Right and Back Middle
            dLM = norm(FL-FM);  %Distance between Front Left and Front Middle
            dLB = norm(FL-BM);  %Distance between Front Left and Back Middle
            dMB = norm(FM-BM);  %Distance between Front Middle and Back Middle
            % Check that the distances are reasonable. Large distances or nan
            % indicate a marker has been dropped
            % Then, we build a representation of the Sensor Reference Frame
            if all([dRL, dRM, dLM]<tolerance)
                xAxis = FR - FL;             %Unnormalized x-axis
                xAxis = xAxis./norm(xAxis);  %Normalize to get the x-axis unit vector
                vecML = FM - FL;             %Vector with x and y components
                xProj = dot(xAxis, vecML)*xAxis;
                yAxis = vecML - xProj;
                yAxis = yAxis./norm(yAxis);  %Normalize to get y-axis unit vector
                zAxis = cross(xAxis, yAxis); %Cross product to get the z-axis unit vector
                zAxis = zAxis./norm(zAxis);
                
%                 hold on;
%                 quiver(0,0,xAxis(1),xAxis(2));
%                 quiver(0,0,yAxis(1),yAxis(2));
%                 axis square;
            elseif all([dRL, dRB, dLB] < tolerance)
                xAxis = FR - FL;             %Unnormalized x-axis
                xAxis = xAxis./norm(xAxis);  %Normalize to get the x-axis unit vector
                vecML = BM - FL;             %Vector with x and y components
                xProj = dot(xAxis, vecML)*xAxis;
                yAxis = -1*(vecML - xProj);  %The vector we get here points in -Y. The -1* up front flips the orientation to +Y
                yAxis = yAxis./norm(yAxis);  %Normalize to get y-axis unit vector
                zAxis = cross(xAxis, yAxis); %Cross product to get the z-axis unit vector
                zAxis = zAxis./norm(zAxis);
            elseif all([dMB, dLM, dLB] < tolerance)
                yAxis = FM - BM;            %Compute the y-axis first, as part of the x-axis is no good
                yAxis = yAxis./norm(yAxis);
                vecML = FM - FL;
                yProj = dot(yAxis, vecML)*yAxis;    %Now we take the projection along the y-axis, instead of along the x-axis
                xAxis = vecML - yProj;
                xAxis = xAxis./norm(xAxis);
                zAxis = cross(xAxis, yAxis);
                zAxis = zAxis./norm(zAxis);
            elseif all([dMB, dRM, dRB] < tolerance)
                yAxis = FM - BM;            %Compute the y-axis first, as part of the x-axis is no good
                yAxis = yAxis./norm(yAxis);
                vecML = FM - FR;
                yProj = dot(yAxis, vecML)*yAxis;    %Now we take the projection along the y-axis, instead of along the x-axis
                xAxis = yProj - vecML;
                xAxis = xAxis./norm(xAxis);
                zAxis = cross(xAxis, yAxis);
                zAxis = zAxis./norm(zAxis);
            else
                xAxis = [1 0 0]';
                yAxis = [0 1 0]';
                zAxis = [0 0 1]';
            end
            %Now the axes of the force sensor are defined in terms of Vicon
            %Coordinates. We can define a rotation from the Sensor frame to the
            %Vicon Frame using the axes. 
            %
            %The axes we've calculated define the axes of the force sensor
            %relative to Vicon coordinates. Combined as columns, they form
            %a rotation from Vicon coordinates to the sensor coordinates.
            %As rows, they form the rotation from sensor coordinates to
            %Vicon coordinates.
%             Rot = [xAxis, yAxis, zAxis]'; % old code incorrect!
            Rot = [xAxis, yAxis, zAxis];
            % Now we can return the forces and torques expressed in the Vicon
            % frame:
            Force(:,n) = Rot*Force(:,n);
            Torque(:,n) = Rot*Torque(:,n);
        end
    end
end

