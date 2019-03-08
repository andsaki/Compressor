clear all;

threshold = -20;
ratio = 1/4;
[y, fs] = audioread('Tears in Heaven.wav');
comp_y = compressor(y,fs,threshold,ratio);
%soundsc(comp_y,fs); 

