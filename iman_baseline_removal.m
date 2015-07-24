%IMAN_BASELINE_REMOVAL
%   Remove scalar baseline from images and trucate negative values (noise)

function im = iman_baseline_removal(im, bval)
%Remove baseline value
im = im - bval;
%Truncate negative values to zero (assume noise in transmission)
im(im < 0) = 0;
end