%IMAN_XYSHIFT
%   Correct for shift in XY position during imaging

function [xyin] = iman_xyshift(xyin, p, shd)

%Assess direction of shift desired (no scaling allowed)
if ~exist('shd', 'var'); shd = 1; else shd = sign(shd); end

%Perform correction
if isstruct(xyin) && isfield(xyin, 'xCoord')            %IF a structure
    for ts = 1:numel(p.frame) %FOR each shift event
        %Adjust x Coordinates
        temp = arrayfun(@(x)x.xCoord + shd*p.dx(ts), ...
            xyin(p.frame(ts):end), 'UniformOutput', false);
        [xyin(p.frame(ts):end).xCoord] = deal(temp{:});
        %Adjust y Coordinates
        temp = arrayfun(@(x)x.yCoord + shd*p.dy(ts), ...
            xyin(p.frame(ts):end), 'UniformOutput', false);
        [xyin(p.frame(ts):end).yCoord] = deal(temp{:});
    end
    %Eliminate negative valued coordinates (u-Track will error)
    ispos = arrayfun(@(x)x.xCoord >= 1 & x.yCoord >= 1, xyin,...
        'UniformOutput', false);
    xyin = arrayfun(@(x,xispos)structfun(@(xx)xx(xispos{1}), ...
                                 x, 'UniformOutput', false), xyin, ispos);
    
elseif isnumeric(xyin)                                  %IF an array
    for ts = 1:numel(p.frame) %FOR each shift event
        %Adjust x Coordinates
        xyin(:, p.frame(ts):end, 1) = ...
            xyin(:, p.frame(ts):end, 1) + shd*p.dx(ts);
        %Adjust y Coordinates
        xyin(:, p.frame(ts):end, 2) = ...
            xyin(:, p.frame(ts):end, 2) + shd*p.dy(ts);
    end
end




end
