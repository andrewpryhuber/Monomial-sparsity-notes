% reads the file polyfile line by line and constructs a cell array of 
% strings representing each polynomial

function polyarray = readPolysFromFile(filename)

file = fopen( filename );

polyarray = textscan(file, '%s');

fclose(file);

end



