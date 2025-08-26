data_in = data;
data_bodyFrame = labToBodyFrame(data_in, true);
data_out = bodyToLabFrame(data_bodyFrame);

% compare data_in chord values and data_out chord vals
res_ind = 3000000;
resFrameNum = data_in.res(:,1);
frameNum = 49;

res_in = data_in.res(resFrameNum==frameNum,2:4);
res_out = data_out.res(resFrameNum==frameNum,2:4);
res_bodyFrame = data_bodyFrame.res(resFrameNum==frameNum,2:4);

scatter3(res_in(:,1),res_in(:,2),res_in(:,3))
hold on
scatter3(res_out(:,1),res_out(:,2),res_out(:,3))
% scatter3(res_bodyFrame(:,1),res_bodyFrame(:,2),res_bodyFrame(:,3))
hold off

axis equal