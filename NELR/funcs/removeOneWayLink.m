function [newW] = removeOneWayLink(W)
    newW = zeros(size(W));
    newW(W'>0) = W(W'>0);
end


