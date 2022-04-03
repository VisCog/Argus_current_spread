function [data,electrode] = csv2struct(filename)

% load in the csv file

% electrode names
count =0;
ename = {};
xval={'A' ;'B'; 'C'; 'D'; 'E'; 'F'};
for i=1:length(xval)
    for j=0:9
        count = count+1;
        eName{count} = sprintf('%s0%d',xval{i},j);
    end
end


[a,b,c] = xlsread(filename);

electrodeCol = find(strcmp('PTS_Electrodes',b(1,:)));
ampCol = find(strcmp('PTS_Amplitude',b(1,:)));
respCol = find(strcmp('KEY',b(1,:)));
trialCol = find(strcmp('Trial No.',b(1,:)));
subjectCol = find(strcmp('subject_id',b(1,:)));


electrodesUsed = unique(c(3:end,electrodeCol));

data.subject = c{2,subjectCol};

for i=1:size(c,1)
    trialNum = c{i,trialCol};
    if ~isempty(trialNum) && ~ischar(trialNum) && ~isnan(trialNum)
        curElectrode = c{i,electrodeCol};
        if strcmp(curElectrode,'none')
            data.electrodeNum(trialNum)=0;
            data.amp(trialNum) = 0;  % amplitude is zero for a catch trial
        else
            data.electrodeNum(trialNum) = find(strcmp(eName,curElectrode));
            data.amp(trialNum) = c{i,ampCol};
        end
        resp = c{i,respCol};
        switch resp
            case 'Y'
                data.resp(trialNum) = 1;
            case 'N'
                data.resp(trialNum) = 0;
        end
    end
end
        
for i=1:length(eName)
    electrode.x(i) = find(strcmp(eName{i}(1),xval));
    electrode.y(i) = str2num(eName{i}(2:end));
end


% Calculate P(yes | catch trial)

electrode.pCatch = mean(data.resp(data.electrodeNum==0));
electrode.pThresh = .5;
electrode.pLapse = 0;

[electrode.xPlot,electrode.yPlot] = meshgrid(0:.5:(max(electrode.x)+1),0:.5:(max(electrode.y)+1));

data.x = electrode.x;
data.y = electrode.y;



