classdef c_eos < handle
    % c_eos   Mie–Grüneisen equation of state 
    %
    %    P = { rho0 * C0^2 (eta - 1) [ eta - (Gamma0/2) (eta - 1)] } /
    %        { eta - s (eta -1) }^2 
    %      + Gamma0 *(e - e0)  for rho >= rho0        // from Wikipedia
    %
    %    P = rho0 * c0^2 (eta - 1) + Gamma0 (e - e0)  for rho < rho0  // from xRAGE
    %    E = rho0 * cv (T - T0) = (e - e0) = (e - e0)
    %    e0: reference internal enery density = rho0 * cv * T0   
    %
    %    For copper 
    %    rho0 = 8.96 g/cc 
    %    cv   = 3.9e+06 erg /(g-K)
    %    C0   = 3.933e+05 cm/sec 
    %    s    = 1.5
    %    Gamma0 = 1.99 (for T < T1)
    %      = 2.12 (for T >= T1)
    %
    %    T1 = 700 K  
    %    eta = rho /rho0   
    %
    %    solid sound speed 
    %        cp = sqrt( E (1-nu)/[ rho (1 + nu)(1 - 2 nu)] ) 
    %        cs = sqrt (G/rho) 
    %        E: Young's modulus
    %        nu: Poisson's ratio  
    %        G: shear modules of elastic material 
    %
    %   for copper: 
    %       E = 110 Gpa = 1.1e12 dyne/cm^2  
    %       nu = 0.33  
    
    properties (Constant)
        % Constants for Copper
        rho0_copper = 8.96;                     % Density of copper (g/cm^3)
        c0_copper   = 3.933e+05;                % Speed of sound in copper (cm/s)
        cv_copper   = 3.9e+06;                  % Specific heat at constant volume (erg/gK)
        s_copper    = 1.5;                      % Dimensionless parameter
        T1_copper   = 700.0;                    % Temperature T1 for copper (K)
        T0_copper   = 290.0;                    % Room temperature, reference (K)
        
        % Gamma constants
        gamma0_less = 1.99;                     % Gamma value for lower pressure/temperature
        gamma0_more = 2.12;                     % Gamma value for higher pressure/temperature
        
        % Mechanical properties of copper
        E_copper  = 1.1e+12;                    % Young's modulus for copper (dyne/cm^2)
        nu_copper = 0.33;                       % Poisson's ratio for copper

        e_T1_copper = c_eos.rho0_copper * c_eos.cv_copper * c_eos.T1_copper; % Internal energy at T1 (erg/cm^3)
        e_T0_copper = c_eos.rho0_copper * c_eos.cv_copper * c_eos.T0_copper; % Internal energy at T0 (erg/cm^3)

    end
    
    methods
%%%%%%%%%%%%%%% implemented and tested
        function cs = sound_speed_solid(obj, rho)
            % Function: sound_speed_solid
            % Inputs:
            %   rho - Density of the material (g/cm^3)
            %
            % Outputs:
            %   cs  - Calculated sound speed in the solid material (cm/s)
            cs = obj.E_copper * (1.0 - obj.nu_copper) / ...
                 (rho * (1 + obj.nu_copper) * (1.0 - obj.nu_copper - obj.nu_copper));
            cs = sqrt(cs);
        end

%%%%%%%%%%%%%%% implemented but not tested        
        function p = p_mie_gruneisen(obj, rho, ei)
            % Function: p_mie_gruneisen
            % Inputs:
            %   rho - Density of the material (g/cm^3)
            %   ei  - Internal energy density (erg/cm^3)
            %
            % Outputs:
            %   p   - Pressure (dyne/cm^2) calculated using the Mie-Grüneisen
            %         equation of state

            % Calculate the compression ratio
            eta = rho / obj.rho0_copper;
            eta1 = eta - 1.0;

            % Determine the appropriate gamma0 value based on internal energy
            gamma0 = obj.gamma0_less;

            % Calculate the denominator for the pressure equation
            denom = eta - obj.s_copper * eta1;
            denom = denom * denom;

            % Compute the pressure using the Mie-Grüneisen equation of state
            p = obj.rho0_copper * obj.c0_copper^2 * eta1 * ...
                (eta - 0.5 * gamma0 * eta1) / denom + ...
                gamma0 * (ei - obj.e_T0_copper);
        end

        function ei = e_mie_gruneisen(obj, rho, p)
            % Function: e_mie_gruneisen
            % Inputs:
            %   rho - Density of the material (g/cm^3)
            %   p   - Pressure (dyne/cm^2)
            %
            % Outputs:
            %   ei  - Internal energy density (erg/cm^3) calculated using the
            %         Mie-Grüneisen equation of state

            % Calculate the compression ratio
            eta = rho / obj.rho0_copper;
            eta1 = eta - 1.0;

            % Calculate the denominator and intermediate variables
            denom = eta - obj.s_copper * eta1;
            denom = denom * denom;
            rhocc02eta1 = obj.rho0_copper * obj.c0_copper^2 * eta1;

            % Calculate internal energy densities for both gamma values
            eiless = (p - rhocc02eta1 * (eta - 0.5 * obj.gamma0_less * eta1) / denom) / ...
                obj.gamma0_less + obj.e_T0_copper;
            eimore = (p - rhocc02eta1 * (eta - 0.5 * obj.gamma0_more * eta1) / denom) / ...
                obj.gamma0_more + obj.e_T0_copper;

            ei = eiless;
        end

    end % methods
end
