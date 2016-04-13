function [ r_x, r_y ] = setNeighbor( q_x, q_y, step, i )
%SETNEIGHBOR Summary of this function goes here
%   Detailed explanation goes here
switch i
    case 1
        r_x = q_x + step;
        r_y = q_y;
    case 2
        r_x = q_x + step;
        r_y = q_y - step;
    case 3
        r_x = q_x;
        r_y = q_y - step;
    case 4
        r_x = q_x - step;
        r_y = q_y - step;
    case 5
        r_x = q_x - step;
        r_y = q_y;
    case 6
        r_x = q_x - step;
        r_y = q_y + step;
    case 7
        r_x = q_x;
        r_y = q_y + step;
    case 8
        r_x = q_x + step;
        r_y = q_y + step;
end

end

