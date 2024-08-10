function [ branchOut ] = branch_2H( channelH )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
numChan = length(channelH);
b(1) = -0.7071-0.7071i;
b(2) = -0.7071+0.7071i;
b(3) = 0.7071+0.7071i;
b(4) = 0.7071-0.7071i;

v(1,:) = [b(1),b(1)];
v(2,:) = [b(1),b(2)];
v(3,:) = [b(1),b(4)];

v(4,:) = [b(2),b(1)];
v(5,:) = [b(2),b(2)];
v(6,:) = [b(2),b(3)];

v(7,:) = [b(3),b(2)];
v(8,:) = [b(3),b(3)];
v(9,:) = [b(3),b(4)];

v(10,:) = [b(4),b(1)];
v(11,:) = [b(4),b(3)];
v(12,:) = [b(4),b(4)];

branch = zeros(12,2);

for kk = 1:12
    branch(kk,:) = filter(channelH,1,v(kk,:));
end
branchOut = branch(:,numChan:end);

end

