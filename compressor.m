function comp_y = compressor(y,fs,threshold,ratio)
comp_y = y;
threshold = 10.^(threshold/20);%閾値;
N = 256;
L = 80;
FR = frameindex(N,floor(N-L),length(comp_y));
comp_FR = comp_y(FR).*hann(N);
fnum = size(FR);
E = zeros(1,fnum(2));

%エレルギーの計算
for j = 1:fnum(2)
    k = 0;
    for i = 1:fnum(1)
       k = k + comp_FR(i,j).^2;       
    end
    E(j) = k;
end

RMS = zeros(1,fnum(2));
%RMS変換
for i = 1:length(E)
    RMS(i) = sqrt(E(i)/N);
end

%閾値を超えたインデクスを保存
RP = 0;
count = 0;
for k = 1:length(RMS)
    if RMS(k) > threshold
        RP = L*k - 1;
       %diff = comp_y(RP) - TS;
       %comp_y(RP) = TS + diff*ratio;
       %comp_y(RP) = comp_y(RP)*ratio;
        count = count + 1;
    end
end
CP = zeros(2,count);
p = 1;
for k = 1:length(RMS)
    if RMS(k) > threshold
       RP = L*k - 1;
       CP(1,p) = RP;%インデクス
       CP(2,p) = RMS(k);%RMS値
       p = p + 1;
    end
end

%L倍し、エネルギーと波形のインデクスを合わせる。
cidx = zeros(count,L+1);
for i = 1:count
    cidx(i,:) = CP(1,i):L+CP(1,i);
end

%波形にコンプをかける。
for i = 1:count
    k = CP(2,i) - threshold;
    r = k*ratio + threshold;
    CompRatio = r/CP(2,i);
    comp_y(cidx(i,:)) = comp_y(cidx(i,:)) * CompRatio;
end

%コンプ後のエレルギー
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
%コンプRMS変換
for i = 1:length(EC)
    RMSC(i) = sqrt(EC(i)/N);
end
    
%プロット
%Y波形プロット
subplot(4,1,1);
t2 =  0:1/fs:(length(y)-1)/fs; %波時間
plot(t2,y); %波形

%RMSプロット
subplot(4,1,2);
t1 = 0:L/fs:(fnum(2)-1)*L/fs; %RMS時間
plot(t1,RMS); %RMS
hold on;
RTL = ones(1,length(RMS))*threshold;%RMStreshold
line(t1,RTL);
hold off;

%comp波形プロット
subplot(4,1,3);
t2 =  0:1/fs:(length(comp_y)-1)/fs; %波時間
plot(t2,comp_y); %波形

%RMSプロット
subplot(4,1,4);
t1 = 0:L/fs:(fnumC(2)-1)*L/fs; %RMS時間
plot(t1,RMSC); %RMS
hold on;
RTL = ones(1,length(RMSC))*threshold;%RMStreshold
line(t1,RTL);
hold off;
end