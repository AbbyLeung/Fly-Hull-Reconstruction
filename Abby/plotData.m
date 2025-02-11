% convert cines to mp4 test
% this takes ~70s, much faster than our current script?!

%% plot some data from data_cleaned
flyAngles = data_cleaned.anglesBodyFrameSmooth;
frames = data_cleaned.params.startTrackingTime:data_cleaned.params.endTrackingTime;
t_ms = frames*(1/data_cleaned.params.fps)*1000;

t_cond = t_ms> -50 & t_ms < 100;

angleInd = 3;
currAngle1 = flyAngles(t_cond,angleInd);
currAngle2 = flyAngles(t_cond,angleInd+3);

% fig params
x0=100;
y0=100;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])

plot(t_ms(t_cond),currAngle1)
hold on
plot(t_ms(t_cond),currAngle2)
hold off

legend('Right','Left')
xlabel('Time (ms)')
ylabel('Stroke (deg)')
% ylim([-10,230])

xlim([-50,100])

%% plot left minus right
flyAngles = data_cleaned.anglesBodyFrameSmooth;
frames = data_cleaned.params.startTrackingTime:data_cleaned.params.endTrackingTime;
t_ms = frames*(1/data_cleaned.params.fps)*1000;

t_cond = t_ms> -50 & t_ms < 100;

% angleInd = 5;
strokeR = flyAngles(t_cond,3)*(-1);
strokeL = flyAngles(t_cond,6);

% LRminus = strokeL-strokeR;
% 
[upperEnvL,lowEnvL] = envelope(strokeL,30,'peak');
[upperEnvR,lowEnvR] = envelope(strokeR,30,'peak');

plot(upperEnvL)
hold on
plot(strokeL)
plot(lowEnvL)
hold off
plot(t_ms(t_cond),(lowEnvL-lowEnvR))
xline([0,50])

% plot(t_ms(t_cond),(strokeL-strokeR))
% xlabel('Time (ms)')
% ylabel('Left - right stroke')

% plot(t_ms(t_cond),strokeL)


%% plot the manually corrected data
flyAngles = data.anglesBodyFrameSmooth;
frames = data.params.startTrackingTime:data.params.endTrackingTime;
t_ms = frames*(1/data.params.fps)*1000;

t_cond = t_ms> -50 & t_ms < 100;

angleInd = 3;
currAngle1 = flyAngles(t_cond,angleInd)*(-1);
currAngle2 = flyAngles(t_cond,angleInd+3);

% fig params
x0=100;
y0=100;
width=900;
height=400;
set(gcf,'position',[x0,y0,width,height])

plot(t_ms(t_cond),currAngle1)
hold on
plot(t_ms(t_cond),currAngle2)
hold off

legend('Right','Left')
xlabel('Time (ms)')
ylabel('Stroke (deg)')
% ylim([-10,230])

xlim([-20,100])

%% plot from all_fly_bw

% mat = all_fly_bw.mat;
dims = all_fly_bw.dim;

camInd = 3;
frameNum = 100;

bw = getImage4D(all_fly_bw, camInd, 105);
imshow(bw)


