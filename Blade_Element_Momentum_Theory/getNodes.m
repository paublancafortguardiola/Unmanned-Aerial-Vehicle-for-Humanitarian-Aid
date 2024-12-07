function [nodes, dy] = getNodes(y0,y,n)
    % returns a vector with the position of the centered nodes of each element given the start position, the end position,
    % and the number of elements
    dy = (y-y0)/n;
    nodes = y0:dy:y;
    nodes = nodes(2:end) - dy/2;
    dy = ones(1,n)*dy;
end