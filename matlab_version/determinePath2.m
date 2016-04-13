function [ xPathMap, yPathMap ] = determinePath2( localCost, Xpos, Ypos )

% imageSize
sizeX = size(localCost,1);
sizeY = size(localCost,2);
% expanded/processed state
state = zeros(sizeX, sizeY); % INITIAL=0, ACTIVE=1, EXPANDED=2
% empty active list (priority queue)
global pq
pq = zeros(3,0); 
% total cost matrix
tCost = Inf(sizeX, sizeY);
tCost(Xpos, Ypos) = 0.0; %seed with total cost 0
%insert seed to active list
insert(Xpos, Ypos, 0.0);
ExpandedCount = 0;

%%% --- finish initializing, start algorithm --- %%%

while ~isempty(pq)
    [q_x,q_y,q_cost]=extractmin; % remove minimumcost pixel q from active list
    state(q_x, q_y) = 2; % set as expanded
    ExpandedCount = ExpandedCount +1;
    step = 1;
    for i = 1:8 % loop around neighbors
        [r_x, r_y] = setNeighbor(q_x, q_y, step, i);
        if state(r_x,r_y) ~= 2 % r not expanded
            temptCost = tCost + localCost(q_x,q_y,i);
            if (state(r_x, r_y) == 1) && (temptCost < tCost(r_x, r_y))
                state(r_x, r_y) = 0;
            end
            if state(r_x, r_y) == 0
                tCost(r_x, r_y) = temptCost;
                insert(r_x, r_y,temptCost);
                state(r_x, r_y) = 1;
            end
        end
    end
end

end

