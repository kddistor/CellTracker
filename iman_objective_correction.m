%IMAN_OBJECTIVE_CORRECTION
%	Correct for objective view area (geometric collection efficiency)
%	across an image.
%
%	[imout, e_inv] = iman_objective_correction(im, p, e_inv)
%		produces the adjusted image 'imout' and the computed transform
%		'e_inv' that is the nominal collection efficiency divided by the
%		spatially biased estimate.  The procedure must be provided an image
%		'im' and associated metadata parameters 'p'.  Additional inputs may
%		be passed to provide a pre-computed transform, 'e_inv'.
%

%Usage example
%	cim = zeros(size(im_stack));
%	[cim(:,:,1), e_inv] = iman_objective_correction(im_stack(:,:,1), mdata); 
%	for s = 2:size(im_stack,3)
%		cim(:,:,s) = iman_objective_correction(im_stack(:,:,s), mdata, e_inv);
%	end

function [imout, e_inv] = iman_objective_correction(im, p, e_inv)
%Uses base geometric model for correction

%Interpolate the estimated collection efficiency to match image (if not provided)
if ~exist('e_inv','var') || isempty(e_inv)	
	%Estimate collection efficiency over the image (if not provided)
	[dx, e_col] = estimate_collection_efficiency(p);

	%Compile image information
	%Get pixel and image sizes (at the sample, 'object space')
	pix_dx = p.cam.PixSizeX;     pix_dy = p.cam.PixSizeY;
	imsz = size(im);  %Image size
	if numel(imsz) < 3; imsz(3) = 1; end
	
	%Define displacement of center of each pixel from image center
	dxmat = abs((1:imsz(1))' - (imsz(1)+1)/2)*pix_dx;   %x displacements
	dxmat = dxmat(:,ones(imsz(2),1),ones(imsz(3),1));   %Replicate image size
	dymat = abs((1:imsz(2)) - (imsz(2)+1)/2)*pix_dy;    %y displacements
	dymat = dymat(ones(imsz(1),1),:,ones(imsz(3),1));   %Replicate image size
	imdx =  sqrt(dxmat.^2 + dymat.^2);  clear dxmat dymat;
	
	%Perform interpolation
	e_col_im = interp1(dx, e_col, imdx, 'linear', 'extrap');
	
	%Modify original image by spatially modified collecion efficiency estimate
	%	Get nominal collection efficiency
	ceff_nom = ( 1 - sqrt(1 - (p.obj.NA./p.obj.RefIndex).^2) )/2;
	%	Get image correction factors
	e_inv = ceff_nom./e_col_im;
end

%Perform image correction; ceff_nom is then the collection efficiency throughout
imout = im.*e_inv;

end



%% Subfunction to estimate collection efficiency over the image
function [dx, e_col] = estimate_collection_efficiency(p)

%Prepare geometric collection efficiency estimate
%To determine the minor axis of the apparent ellipse of the objective from
%the object point of view, evaluate the vectors to each edge as a function
%of displacement from the center of view.

%Check Working Distance scale (often in mm, not um)
if p.obj.WkDist < 100; p.obj.WkDist = p.obj.WkDist*1000; end

thet = asin(p.obj.NA./p.obj.RefIndex);  %View half-angle
lens_half = tan(thet).*p.obj.WkDist;  %In um (as must be the Working Distance)
%   Maximum displacement expected (from camera, objective)
dx_max = sqrt( (p.cam.PixNumX.*p.cam.PixSizeX)^2 ...
              + (p.cam.PixNumY.*p.cam.PixSizeY)^2) / 2;  %In um          

%Set number of points to evaluate by pixel dimension
npts = ceil(sqrt(p.cam.PixNumX^2 + p.cam.PixNumY^2)/2);

%Prepare points to evaluate
dx = (0:dx_max/(npts-1):dx_max)';
%	Positions of viewing edges
A = [dx + lens_half, ones(npts,1)*p.obj.WkDist];	%Outer
B = [dx - lens_half, ones(npts,1)*p.obj.WkDist];	%Inner

%Major axis angle (half) (Corrects for increased distance from objective)
ta = atan(lens_half./sqrt(p.obj.WkDist.^2 + dx.^2));
%Minor axis angle (half) (Corrects for distance and angle)
tb = acos( dot(A,B,2)./(sqrt(sum(A.^2,2)).*sqrt(sum(B.^2,2))) )./2;

%Ellipse radius, function of phi
%--------------------------------------------------------------------------
%r = a*b/sqrt((b*cos(phi))^2 + (a*sin(phi))^2);
%OR r = b^2 / (a - sqrt(a^2 - b^2)*cos(phi));
%   Substituting the arc length approximating of phi*(focal length)
%   gives the same, with r, a, b being angles, rather than radii
%   Take integral and divide by max to get actual collection efficiency.
%   Area = Int_(-pi)^(pi) Int_0^r(phi) R^2 sin(thet) d_thet d_phi
%    = r^2 * (2*pi - Int_(-pi)^(pi) cos( r(phi) ) d_phi
%   Divide by nominal collection efficiency to give normalized transform.
%--------------------------------------------------------------------------

%I was not clever enough the solve the integral, so evaluating numerically.
%   r(phi) = tb.^2 / (ta - sqrt(ta^2 - tb^2)*cos(phi));
%   Here, already divided by R^2
e_col = ( 2*pi*ones(npts,1) - ...
    integral(@(x)cos(el_rang(ta,tb,x)), -pi, pi, 'ArrayValued', true)...
                    )./(4*pi);
                
%Adjustment for illumination uniformity
%   Formulated for point source illumination focused beyond sample
%   Applicable for focus before sample, as well as 
%Working Distance for illumination focal point
WD_il = p.obj.WkDist.*lens_half./(lens_half - dx_max);

%Illumination spherical radius at center
R_0 = WD_il - p.obj.WkDist;
%Illumination spherical radius at dx from center
R = sqrt(R_0.^2 + dx.^2);

%Surface area of illumination sphere patch
SA = 2*pi*R.^2 .* (1 - cos(thet));
%Relative illumination (to peak at center)
RI = 2*pi*R_0.^2 .* (1 - cos(thet)) ./ SA;

%Total effective collection/illumination efficiency
e_col = e_col.*RI;
end

%% Function for ellipse radial angle in spherical coords
function out = el_rang(a,b,phi)
%with Major Axis (a) and Minor Axis (b), centered on a Focus
% out = b.^2 ./ (a - sqrt(a.^2 - b.^2).*cos(phi)); %Major/Minor uncertain
out = a.*b./sqrt((b.*cos(phi)).^2 + (a.*sin(phi)).^2);
end