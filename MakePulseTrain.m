function [ trl ] = MakePulseTrain(trl)
% generates a pulse train based on the parameters in input structure s

on =  mod(trl.t,1/trl.freq) < trl.pw;
delay =  trl.pw+trl.ip;
lag = round(trl.lag/trl.dt);
off = mod(trl.t-delay,1/trl.freq) < trl.pw;
tmp  = trl.amp.*(on-off);
trl.pt= zeros(1, lag+length(tmp));
trl.pt(lag+1:lag+length(tmp))=tmp;

if trl.dur<trl.trialdur
    trl.pt((end+1):round((trl.trialdur/trl.dt))) = 0;
    trl.t = 0:trl.dt:trl.trialdur-trl.dt;
end
