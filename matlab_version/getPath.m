function [xPath, yPath] = getPath(xPathMap, yPathMap, Xpos, Ypos)

%get int position
Xpos = int16(Xpos);
Ypos = int16(Ypos);

xPathMap  = int16(xPathMap);
yPathMap  = int16(yPathMap);

MAXSTEP = 2000;
%initialize paths
xPath = zeros(MAXSTEP, 1, 'int16');
yPath = zeros(MAXSTEP, 1, 'int16');

step = 1;
xPath(step) = Xpos;
yPath(step) = Ypos;

% Not seed point: go back through path map
while xPathMap(Ypos, Xpos) ~= 0 || yPathMap(Ypos, Xpos) ~= 0
    step = step + 1;
    XposCopy = Xpos;
    Xpos = Xpos + xPathMap(Ypos, Xpos);
    Ypos = Ypos + yPathMap(Ypos, Xpos);
    
    xPath(step) = Xpos;
    yPath(step) = Ypos;
end

% revert to forward path
xPath = xPath(step - 1:-1:1);
yPath = yPath(step - 1:-1:1);

end