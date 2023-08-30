function [ Y ] = soft_thresh( X, thresh )
%SOFT_THRESH Summary of this function goes here
%   Detailed explanation goes here

idx = find(abs(X) > thresh);
Y = zeros(size(X));
Y(idx) = sign(X(idx)).*(abs(X(idx))-thresh);

end

