% Arguments:   
%           xcoord - x coordinate of the point
%           ycoord - y coordinate of the point
%           value  - value of the point
%
% Returns:
%            changed    - is the value of the point changed (used to update
%            the predecessor pointer). 0->not changed, 1->changed

function changed=insert(xcoord,ycoord,value)

    global pq
    % Check if the point already exists in the priority queue
    ans=find( (pq(1,:)==xcoord) & (pq(2,:)==ycoord));
    if ans %point is already in the priority queue
        % now we should check if the new value is greater that the existing
        % one inside the priority queue
        if pq(3,ans(1))> value %Change it 
            pq(3,ans(1))=value;
            changed=1;
            return;
        else % Don't change it
            changed=0;
            return;
        end
    else %point not found on pq, we should insert it
        pq_size=size(pq,2);
        pq(:,pq_size+1)=[xcoord; ycoord; value];
        changed=1;
        return;
    end
