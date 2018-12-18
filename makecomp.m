clear all;
[y, fs] = audioread('15K1103_0_1.wav');
y = resample(y,16000,fs); fs = 16000;
comp_y = y;
threshold = 0.03;
ratio = 1/8;
attack = 0.00002 * fs;
inv_ratio = 1/ratio; 
count = 0;
GS = 0;
a = 0;
for i = 1:length(comp_y)
    if comp_y(i) > threshold
        if count >= attack
            diff = comp_y(i) - threshold;
            comp_y(i) = threshold + diff*ratio;
            GS = GS + 1;
        else
            diff = comp_y(i) - threshold;
            AttackRatio = ratio * (inv_ratio - (inv_ratio - 1) * count/attack);
            comp_y(i) = threshold + diff*AttackRatio;
            count = count + 1;
            a = a + 1;
        end
    elseif comp_y(i) < -threshold
        if count >= attack
            diff = comp_y(i) + threshold;
            comp_y(i) = -threshold + diff*ratio;
            GS = GS + 1;
        else
            diff = comp_y(i) + threshold;
            AttackRatio = ratio * (inv_ratio - (inv_ratio - 1) * count/attack);
            comp_y(i) = -threshold + diff* AttackRatio;
            count = count + 1;
        end
    else
        count = 0;
    end
end

subplot(2,1,1);
t = 0:1/fs:(length(comp_y)-1)/fs;
th = ones(1,length(comp_y))*threshold;
plot(t,y);
hold on;
plot(t,th);
hold on;
plot(t,-th);
hold off;
subplot(2,1,2);
plot(t,comp_y);
hold on;
plot(t,th);
hold on;
plot(t,-th);
hold off;

%soundsc(y,fs);
soundsc(comp_y,fs);