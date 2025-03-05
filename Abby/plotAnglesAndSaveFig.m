% script that loops through files and generates figs
expFolder = 'Y:\Abby2\Feb7\Control haltereless';
flyDir = dir(fullfile(expFolder,'fly*'));

% analysisFolder = 'Y:\Abby2\Feb7\S01\fly 18\Analysis';
for i = 1:length(flyDir)
    analysisFolder = fullfile(expFolder,flyDir(i).name,'Analysis');
    analysisDir = dir(fullfile(analysisFolder,'Expr*'));
    saveToFolder = fullfile(analysisFolder,'figs');
    
    if ~isfolder(saveToFolder)
        mkdir(saveToFolder)
    end
    
    for folderInd = 1:length(analysisDir)
        currFolderName = analysisDir(folderInd).name;
        dataFilename = [currFolderName,'_cleaned.mat'];
    
        if ~isfile(fullfile(analysisFolder,currFolderName,dataFilename))
            disp(['Cannot find ',dataFilename])
            continue
        end
    
        currData = load(fullfile(analysisFolder,currFolderName,dataFilename));
        data = currData.data_cleaned;
    
        flyAngles = data.anglesBodyFrameSmooth;
        flyAngles(:,3) = flyAngles(:,3)*(-1); % flip because defined as negative
    
        frames = data.params.startTrackingTime:data.params.endTrackingTime;
        t_ms = frames*(1/data.params.fps)*1000;
        t_cond = t_ms> -50 & t_ms < 100;
    
        anglesFig = tiledlayout(3,1);
        
        nexttile
        angleInd = 3;
        currAngle1 = flyAngles(t_cond,angleInd);
        currAngle2 = flyAngles(t_cond,angleInd+3);
        
        plot(t_ms(t_cond),currAngle1)
        hold on
        plot(t_ms(t_cond),currAngle2)
        hold off
        
        legend('Right','Left')
        xlabel('Time (ms)')
        ylabel('Stroke (deg)')
        
        xlim([-50,100])
        ylim([0,220])
    
        nexttile
        angleInd = 4;
        currAngle1 = flyAngles(t_cond,angleInd);
        currAngle2 = flyAngles(t_cond,angleInd+3);
        
        plot(t_ms(t_cond),currAngle1)
        hold on
        plot(t_ms(t_cond),currAngle2)
        hold off
        
    %     legend('Right','Left')
        xlabel('Time (ms)')
        ylabel('Deviation (deg)')
        xlim([-50,100])
        ylim([-30,90])
    
        nexttile % pitch angle
        angleInd = 5;
        currAngle1 = flyAngles(t_cond,angleInd);
        currAngle2 = flyAngles(t_cond,angleInd+3);
        
        plot(t_ms(t_cond),currAngle1)
        hold on
        plot(t_ms(t_cond),currAngle2)
        hold off
        
    %     legend('Right','Left')
        xlabel('Time (ms)')
        ylabel('Pitch (deg)')
        xlim([-50,100])
        ylim([0,200])
        
        title(anglesFig,strrep(currFolderName,'_',' '))
        saveas(gcf,fullfile(saveToFolder,[currFolderName,'.png']))
        close
    end
end
