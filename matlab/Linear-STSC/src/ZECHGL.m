function[ET]=ZECHGL(N)
    % MATLAB equivalent of the ZECHGL subroutine
    

    % Check if N is zero, if so, return
    if N == 0
        return;
    end

    % Initialize ET(0) to -1
    ET(1) = -1.0;

    % Check if N is 1, if so, return
    if N == 1
        return;
    end

    % Initialize ET(N) to 1
    ET(N+1) = 1.0;

    % Check if N is 2, if so, return
    if N == 2
        return;
    end

    % Calculate N2
    N2 = round((N - 1)/(2));

    % Initialize ET(N2+1) to 0
    ET(N2 + 1) = 0.0;

    % Check if N is 2, if so, return
    if N == 2
        return;
    end

    % Define PI
    PI = 3.14159265358979323846;

    % Calculate constant C
    C = PI / double(N);

    % Loop to calculate ET values
    for I = 1:N2
        ETX = cos(C * double(I));
        ET(I+1) = -ETX;
        ET(N - I+1) = ETX;
    end

    % Note: MATLAB indexing starts from 1, so the loop index I is used directly.

end
