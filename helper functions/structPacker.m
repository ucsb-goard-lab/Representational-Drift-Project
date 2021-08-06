

function datastruct = structPacker(datastruct, varargin)
% data structure packer
	%elements of varargin are ordered as: variable, name, variable name, etc
	if isempty(datastruct)
	    datastruct = struct();
	end
	for ii = 1:2:nargin-1
	    datastruct.(varargin{ii+1}) = varargin{ii};
	end
end