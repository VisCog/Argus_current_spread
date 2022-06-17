close all
trl.dur = 100/1000;
trl.trialdur = 100/1000;
trl.dt = .0001;
trl.t = 0:trl.dt:trl.dur-trl.dt;
trl.freq = 20;
trl.lag = 10/1000;
trl.ip = 1/1000;
trl.pw = .46/1000;
trl.amp = 1;

[ trl ] = MakePulseTrain(trl)
subplot(2,1,1)
plot(trl.t, trl.pt(1:length(trl.t)));

trl.lag = 35/1000;
[ trl ] = MakePulseTrain(trl)
subplot(2,1,2)
plot(trl.t, trl.pt(1:length(trl.t)));

figure
trl.dur = 5/1000;
trl.trialdur = trl.dur; 
trl.dt = .00001;
trl.t = 0:trl.dt:trl.dur-trl.dt;
trl.lag = 1.5/1000;

[ trl ] = MakePulseTrain(trl)
subplot(2,1,1)
plot(trl.t, trl.pt(1:length(trl.t)));