function [matrix, rowNames, colNames] = readMtxFile(filename)
    % readMtxFile Reads a .mtx file and returns the matrix, row, and column names
    %
    % Input:
    %   filename - Name of the .mtx file to read
    %
    % Outputs:
    %   matrix - Sparse matrix stored in the .mtx file
    %   rowNames - Optional row identifiers (if available in the file)
    %   colNames - Optional column identifiers (if available in the file)

    % Open the file
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('File cannot be opened: %s', filename);
    end
    
    % Initialize names (if metadata is available in the file)
    rowNames = {};
    colNames = {};
    
    % Read the first line to check if it contains metadata or dimensions
    firstLine = fgetl(fileID);
    if startsWith(firstLine, '%%MatrixMarket')
        % Skip header and comments
        while startsWith(firstLine, '%')
            firstLine = fgetl(fileID);
        end
    end
    
    % Read matrix dimensions and non-zero count
    dims = sscanf(firstLine, '%d %d %d', 3);
    numRows = dims(1);
    numCols = dims(2);
    numEntries = dims(3);
    
    % Initialize sparse matrix
    i = zeros(numEntries, 1);
    j = zeros(numEntries, 1);
    s = zeros(numEntries, 1);
    
    % Read non-zero entries
    for k = 1:numEntries
        data = fscanf(fileID, '%d %d %f', [1, 3]);
        i(k) = data(1);
        j(k) = data(2);
        s(k) = data(3);
    end
    
    % Close the file
    fclose(fileID);
    
    % Create sparse matrix
    matrix = sparse(i, j, s, numRows, numCols);
    
    % If the .mtx file contains row and column names or other metadata, additional parsing logic would be needed here.
    
    % Output
    fprintf('Matrix of size %dx%d with %d non-zero entries loaded.\n', numRows, numCols, numEntries);
end
