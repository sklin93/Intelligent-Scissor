function L = getCost(Img)
%%%%%%%%% Calculate local cost at node %%%%%%%%

Wz = 0.5;
Wg = 0.5;
Wd = 0.43;

% The zero-crossing cost Fz
Fz = ~edge(Img, 'zerocross');

% The gradient strength cost Fg
Img = double(Img);
[dX,dY] = imgradientxy(Img);
Fg = sqrt(dX.^2 + dY.^2);
Fg = 1 - Fg./max(Fg(:));

% Fd

%Sum:
L = Wz.*double(Fz)+ Wg.*Fg;

end

