function printTextWithOK(text,width)

    % Define total width for each line (including the text, padding, and "OK")
    totalWidth = width;  % You can change this value to control the total length
    
    % Define a fixed number of '-' characters to add in front of the text
    frontPadding = '----- ';
    
    % Calculate the length of the padding based on the text length
    textLength = strlength(text) + length(frontPadding);  % Account for the front padding
    paddingLength = totalWidth - textLength - 3;  % Subtract 3 for space and 'OK'
    
    % Ensure padding is non-negative
    if paddingLength < 0
        paddingLength = 0;
    end
    
    % Create the padding string with '-' characters
    padding = repmat('-', 1, paddingLength);
    
    % Print the formatted line with front padding, text, and "OK"
    fprintf('%s%s %s OK\n', frontPadding, text, padding);
end
