function ret = create_currentspread(ret)

%#codegen
function y = fcn(S)

% Specify the class of the input S as struct.
assert(isstruct(S));

% Specify the size of the fields r and i
% based on the first element of the array.
assert(all(size(S) == [1 2]));
assert(isa(S(1).r,'double'));
assert(isa(S(1).i,'int8'));


% takes in parameter values describing electrodes and
% calculates the resulting current field

ret.I = zeros(size(ret.e(1).R));
for i = 1:length(ret.x) % for each electrode
    ret.e(i).I = ret.eI./(1+((ret.k.*ret.e(i).R).^ret.a));
    ret.I = ret.I + ret.e(i).I;
end
ret.I(ret.Z<max([ret.z ret.z])) = 0; % remove all current above the electrode
eIxy = ret.I(:, :, round(size(ret.I, 3)/2));
if isfield(ret, 'loc') % if choosing a particular location on the retina at which to know amplitude
    disp('calculating threshold at a specific location')
    ret.eIxy = interpn(unique(ret.Y), unique(ret.X),eIxy, ret.loc(2), ret.loc(1));
else % find the max amplitude
    ret.eIxy = max(eIxy(:));
end
end