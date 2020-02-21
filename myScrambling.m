
function seq = myScrambling(initState,lenSeq)

LFSR = initState;%[1 0 1 1 1 0 1].';%ones(7,1); 
seq = NaN(lenSeq,1);

for ii = 1:length(lenSeq)
    nxt = mod(LFSR(end)+LFSR(4),2);
    seq(ii) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
end