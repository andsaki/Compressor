function findex = frameindex(framelength, noverlap, signallength)
sft=framelength-noverlap;
n=fix((signallength-framelength)/sft+1);
findex=repmat((1:framelength)',1,n)+repmat((0:(n-1))*sft,framelength,1);
end