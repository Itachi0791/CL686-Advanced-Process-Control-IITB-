function M = createSymbolicMatrix(name, m, n)
    % Create a symbolic matrix of size (m, n) with the given name
    M = sym(name, [m, n], 'real');
end

