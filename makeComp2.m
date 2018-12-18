clear all;

[y,fs] = audioread('Audio 05_124.wav');
%[y,fs] = audioread('15K1103_hakumei_7_1.wav');
%[y, fs] = audioread('voicechanger.wav');
comp_y = y;
threshold = -15;
TS = 10.^(threshold/20);
ratio = 1/2;
attack = 0.002 * fs;
N = 256;
L = 80;
FR = frameindex(N,floor(N-L),length(comp_y));
comp_FR = comp_y(FR).*hann(N);
fnum = size(FR);
E = zeros(1,fnum(2));

%�G�����M�[�̌v�Z
for j = 1:fnum(2)
    k = 0;
    for i = 1:fnum(1)
       k = k + comp_FR(i,j).^2;       
    end
    E(j) = k;
end

RMS = zeros(1,fnum(2));
%RMS�ϊ�
for a = 1:length(E)
    RMS(a) = sqrt(E(a)/N);
end

%Treshold����
RP = 0;
count = 0;
for k = 1:length(RMS)
    if RMS(k) > TS
       RP = L*k - 1;
       %diff = comp_y(RP) - TS;
       %comp_y(RP) = TS + diff*ratio;
       %comp_y(RP) = comp_y(RP)*ratio;
       count = count + 1;
    end
end
CP = zeros(1,count);
p = 1;
for k = 1:length(RMS)
    if RMS(k) > TS
       RP = L*k - 1;
       CP(p) = RP;
       %diff = comp_y(RP) - TS;
       %comp_y(RP) = TS + diff*ratio;
       %comp_y(RP) = comp_y(RP)*ratio;
       p = p + 1;
    end
end

%L�{
S = zeros(count,L+1);
for k = 1:length(CP)
    S(k,:) = CP(k):CP(k)+L;
end

u = 0;
%�g�`�ɃR���v
for w = 1:length(CP)
    comp_y(S(w,:)) = comp_y(S(w,:))*ratio;
    u = u + 1;
end

%�R���v��̃G�����M�[
N = 256;
L = 80;
FRC = frameindex(N,floor(N-L),length(comp_y));
comp_FRC = comp_y(FRC).*hann(N);
fnumC = size(FRC);
EC = zeros(1,fnumC(2));
for j = 1:fnumC(2)
    k = 0;
    for i = 1:fnumC(1)
       k = k + comp_FRC(i,j).^2;       
    end
    EC(j) = k;
end

RMSC = zeros(1,fnumC(2));
%�R���vRMS�ϊ�
for a = 1:length(EC)
    RMSC(a) = sqrt(EC(a)/N);
end
    
%�v���b�g
%RMS�v���b�g
subplot(4,1,1);
t1 = 0:L/fs:(fnum(2)-1)*L/fs; %RMS����
plot(t1,RMS); %RMS
hold on;
RTL = ones(1,length(RMS))*TS;%RMStreshold
line(t1,RTL);
hold on;
RPL = ones(1,length(RMS))*RP/fs;
line(RPL,RMS);
hold off;

%Y�g�`�v���b�g
subplot(4,1,2);
t2 =  0:1/fs:(length(y)-1)/fs; %�g����
plot(t2,y); %�g�`
hold on;
th = ones(1,length(y))*TS; 
%line(t2,th); %threshold
hold on;
%line(t2,-th); %threshold
hold on;
RPL = ones(1,length(y))*RP/fs;
line(RPL,y);
hold off;

%RMS�v���b�g
subplot(4,1,3);
t1 = 0:L/fs:(fnumC(2)-1)*L/fs; %RMS����
plot(t1,RMSC); %RMS
hold on;
RTL = ones(1,length(RMSC))*TS;%RMStreshold
line(t1,RTL);
hold on;
RPL = ones(1,length(RMSC))*RP/fs;
line(RPL,RMSC);
hold off;

%comp�g�`�v���b�g
subplot(4,1,4);
t2 =  0:1/fs:(length(comp_y)-1)/fs; %�g����
plot(t2,comp_y); %�g�`
hold on;
th = ones(1,length(comp_y))*TS; 
%line(t2,th); %threshold
hold on;
%line(t2,-th); %threshold
hold on;
RPL = ones(1,length(comp_y))*RP/fs;
line(RPL,comp_y);
hold off;
%soundsc(comp_y,fs);
