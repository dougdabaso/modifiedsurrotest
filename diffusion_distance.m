function [val, varargout] = diffusion_distance(im1, im2, sig, dim, threshold, pad_type)
% DIFFUSION_DISTANCE calculate the diffusion distance between matrices
% dist = diffusion_distance(im1, im2, [sig, dim], threshold, pad_type]) 
% [dist, iter] = diffusion_distance(im1, im2) compute the diffusion distance
% between im1 and im2 (histograms). 
% 
% Optionally returns the number of iterations. This number will be limited by
% the size of the histograms, and the threshold value, which can be set
% manually, or will take the value 1e-4.
% 
% The padding method used by the filtering stages can be set using the 4th
% argument, and defaults to circular. For available options, see help imfilter. 
% 
% Reference:
% Diffusion Distance for Histogram Comparison
% Ling and Okada, Proc. CVPR, 2006
%
% Matt Foster <ee1mpf@bath.ac.uk>

error(nargchk(2, 6, nargin, 'struct'))

if nargin < 4
  sig = 2;
  dim = 2;
end

if nargin == 3
  error('Strange number of arguments.');
end

if nargin < 5
  threshold = 1e-6;
end

if nargin < 6
  pad_type = 'symmetric'; % see help imfilter
  %pad_type = 'circular'; % see help imfilter 
end

kernel = fspecial('gaussian', dim, sig);


dist = abs(im1 - im2);

%eu acho que aqui seria a integracao
 %val = sum(abs(dist(:)));
 %val = norm(dist(:));
 val = norm(dist(:),'fro');
 %val = norm(dist(:),5);
 
level_val = 100000; % This need initialising to larger than the threshold.
iter = 0;

while level_val > threshold && all(size(dist) >= size(kernel))
  iter = iter + 1;
  % filter, and downsample
  dist = imresize(imfilter(dist, kernel, pad_type), 0.5);
  level_val = sum(abs(dist(:)));
%  level_val = norm(dist(:));
  val = val + level_val;
end

if nargout >= 2
  varargout{1} = iter;
end

