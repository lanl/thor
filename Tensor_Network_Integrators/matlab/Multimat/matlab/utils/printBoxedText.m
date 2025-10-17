function printBoxedText(text,padding)

    % Define the padding size
    % padding = 20;

    % Length of the text
    textLength = length(text);
    
    % Create the top and bottom border with '-' characters
    border = repmat('-', 1, textLength + padding * 2 + 2); % Add 2 for the '| ' and ' |' borders
    
    % Print the top border
    fprintf('%s\n', border);
    
    % Print the text inside the box with padding
    fprintf('|%s%s%s|\n', repmat(' ', 1, padding), text, repmat(' ', 1, padding));
    
    % Print the bottom border
    fprintf('%s\n', border);
end