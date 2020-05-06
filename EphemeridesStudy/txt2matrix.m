function [matrix, time] = txt2matrix(filename)

%TXT2MATRIX convert a file in the text format '*.txt' to a matrix
% containing the ephemerides of the target body.
%
% INPUTS: 
%   filename:      '*.txt'        
%   Character vector containing the name of the file, followed by the
%   specificator of format ".txt".
%
% OUTPUT: 
%   matrix:     [Nx6] [Km ~ deg deg deg]        
%   Matrix containing the ephemerides [a e i Om om theta] of the target object.
%
%   time:       [Nx1] [s]
%   Time vector containing the time step of the variation of keplerian
%   elements.
%
% FUNCTION CALLS:
%   None
%

% Read the file
text = fileread(filename);

% Split the text in rows
text = splitlines(text);

% Find the index of the matching results
row1 = find(~cellfun('isempty', regexp(text,'\$\$SOE','match')));
row2 = find(~cellfun('isempty', regexp(text,'\$\$EOE','match')));

% Set to empty cell the unwanted text rows
N = length(text);
for k = 1:row1
    text{k} = [];
end

for k = row2:N
    text{k} = [];
end

% Remove empty cells from cell array
text = text(~cellfun('isempty',text));

% Convert to string array
text = string(text);

% Remove unnecessary characters
text = replace(text,[",","A.D."]," ");

% Neglect spaces at the begin/end of the row
text = strip(text);

% Divide the string array in columns at white spaces
text = split(text);

% Keplerian elements  columns from the full matrix, plus the time vector
text = text(:,[1,13,4,6,7,8,12]);    % [s km ~ deg deg deg deg]

% Conversion from string array to table
text = array2table(text);

% Conversion of each string in numeric value
text.text1 = str2num(char(text.text1));
text.text2 = str2num(char(text.text2));
text.text3 = str2num(char(text.text3));
text.text4 = str2num(char(text.text4));
text.text5 = str2num(char(text.text5));
text.text6 = str2num(char(text.text6));
text.text7 = str2num(char(text.text7));

% Vector containing the time span 
time = (text.text1(:,1)'-text.text1(1,1))*86400;

% Conversion from table to matrix
matrix = table2array(text);

% Neglect the first column containing the days 
matrix(:,1) = [];

end







