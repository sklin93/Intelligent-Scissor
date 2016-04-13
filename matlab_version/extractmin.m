function [xcoord,ycoord,cost] = extractmin

global pq
    
    [cost,index]=min(pq(3,:));
    
    if index % the priority queue is NOT empty
        xcoord=pq(1,index);
        ycoord=pq(2,index);
        pq(:,index)=[]; %delete the minimum element
    else
        xcoord=0;
        ycoord=0;
    end
    return
