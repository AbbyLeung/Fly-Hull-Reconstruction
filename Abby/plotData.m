% convert cines to mp4 test
% this takes ~70s, much faster than our current script?!

tic

exprPath = 'Y:\Abby2\02_20012025';
camNames = {'xz','xy','yz'};
cineFilenames = strcat(camNames,'_005.cine');
cinePaths = fullfile(exprPath,cineFilenames);
saveToPath = fullfile(exprPath,'mp4');

if ~isfolder(saveToPath)
    mkdir(saveToPath)
end

% cine2mp4(cinePaths,1,5,saveToPath)

% toc

%% plot some data from data_cleaned
flyAngles = data_cleaned.anglesBodyFrameSmooth;
frames = data_cleaned.params.startTrackingTime:data_cleaned.params.endTrackingTime;
t_ms = frames*(1/data_cleaned.params.fps)*1000;

t_cond = t_ms> -50 & t_ms < 100;

angleInd = 4;
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


