function [ status ] = check_N(t,y,flag,~)
%check_N Summary of this function goes here
%   Detailed explanation goes here
global K
if y(2)==0
    K=K+1;
end
if K==5
    status=1;
else
    status=0;
end
% l=length(ye);
% if l==6
%     status=1;
% else
%     status=0;
% end

end

