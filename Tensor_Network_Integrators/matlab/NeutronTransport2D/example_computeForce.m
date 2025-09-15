function output = example_computeForce(mass, acceleration)
    % computeForce - Calculates force using Newton's Second Law.
    %
    % Syntax:
    %   output = computeForce(mass, acceleration)
    %
    % Description:
    %   This function computes the force exerted on an object based on its
    %   mass and acceleration, using Newton's Second Law:
    %
    %       **F = m × a**
    %
    %   where:
    %       *F* is force in newtons (N),
    %       *m* is mass in kilograms (kg),
    %       *a* is acceleration in meters per second squared (m/s²).
    %
    % Inputs:
    %   mass - Numeric scalar. The mass of the object in kilograms.
    %   acceleration - Numeric scalar or vector. The acceleration applied
    %                  to the object in m/s².
    %
    % Outputs:
    %   output - Numeric scalar or vector. The resulting force in newtons (N).
    %
    % Examples:
    %   % Example 1: Compute force with scalar input
    %   f = computeForce(10, 9.81)
    %
    %   % Example 2: Vectorized acceleration
    %   f = computeForce(2, [0 1 2 3])
    %
    % See also: computeAcceleration, computeMass
    %
    % .. note::
    %    Make sure that units are consistent when passing inputs.
    %
    % .. tip::
    %    Use ``acceleration = [0:0.5:10]`` to analyze force response over a range.
    %
    
        output = mass * acceleration;
    end