% ===================================================
% *** FUNCTION luv2xyz
% ***
% *** function XYZ = luv2xyz(luv, up, vp, obs, xyzw)
% *** computes XYZ from Luv

% *** *** luv, up, vp can be obtained from function 'xyz2luv'

% *** Luv is an n by 3 matrix for L*u*v* (Ls, us, vs)
% *** e.g. set obs to 'd65_64' for D65 and 1964
% *** set obs to 'user' to use optional argument   
% *** xyzw as the XYZ of the white point
% *** This function is created by Qianqian Pan from University of Leeds,
% United Kingdom

    % E.g.,
        % XYZ =[29.54 25.28 19.59; 28.24 24.17 17.43];
        % [luv, up, vp] = xyz2luv(XYZ, 'd65_31');
        % XYZ = luv2xyz(luv, up, vp, 'd65_31')

        % XYZ =

        % 29.5400   25.2800   19.5900
        % 28.2400   24.1700   17.4300


% ===================================================
function XYZ = luv2xyz(luv,up, vp, obs, xyzw)

if (size(luv,2)~=3)
   disp('luv must be n by 3'); return;   
end
XYZ = zeros(size(luv,1),size(luv,2));

if strcmp('a_64',obs)
    white=[111.144 100.00 35.200];
elseif strcmp('a_31', obs)
    white=[109.850 100.00 35.585];
elseif strcmp('c_64', obs)
    white=[97.285 100.00 116.145];
elseif strcmp('c_31', obs)
    white=[98.074 100.00 118.232];
elseif strcmp('d50_64', obs)
    white=[96.720 100.00 81.427];
elseif strcmp('d50_31', obs)
    white=[96.422 100.00 82.521];
elseif strcmp('d55_64', obs)
    white=[95.799 100.00 90.926];
elseif strcmp('d55_31', obs)
    white=[95.682 100.00 92.149];
elseif strcmp('d65_64', obs)
    white=[94.811 100.00 107.304];
elseif strcmp('d65_31', obs)
    white=[95.047 100.00 108.883];
elseif strcmp('d75_64', obs)
    white=[94.416 100.00 120.641];
elseif strcmp('d75_31', obs)
    white=[94.072 100.00 122.638];
elseif strcmp('f2_64', obs)
    white=[103.279 100.00 69.027];
elseif strcmp('f2_31', obs)
    white=[99.186 100.00 67.393];
elseif strcmp('f7_64', obs)
    white=[95.792 100.00 107.686];
elseif strcmp('f7_31', obs)
    white=[95.041 100.00 108.747];
elseif strcmp('f11_64', obs)
    white=[103.863 100.00 65.607]; 
elseif strcmp('f11_31', obs)
    white=[100.962 100.00 64.350];
elseif strcmp('user', obs)
    white=xyzw;
else
   disp('unknown option obs'); 
   disp('use d65_64 for D65 and 1964 observer'); return;
end


Ls = luv(:,1);
us = luv(:,2);
vs = luv(:,3);

var_Y = (Ls + 16)/116;

for i = 1:size(var_Y,1)

    if (var_Y(i)^3 > (6/29)^3)
        var_Y(i) = var_Y(i)^3;
    else
        var_Y(i) = (var_Y(i) - 16/116)/7.787;
    end
end

ReferenceX = white(1);
ReferenceY = white(2);
ReferenceZ = white(3);

ref_U = (4 * ReferenceX)/(ReferenceX + (15 * ReferenceY) + (3 * ReferenceZ));
ref_V = (9 * ReferenceY)/(ReferenceX + (15 * ReferenceY) + (3 * ReferenceZ));

var_U = us./(13 * Ls) + ref_U;
var_V = vs./(13 * Ls) + ref_V;

Y = var_Y * 100;
X = -(9 * Y.* var_U)./((var_U - 4).* var_V - var_U.* var_V);
Z = (9 * Y - (15 * var_V.* Y) - (var_V .* X))./(3 * var_V);

XYZ = [X Y Z];

end


















