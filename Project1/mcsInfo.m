
function rateStruct = mcsInfo(RATE)
% mcs stands for Modulation and Coding Scheme
% RATE is the value of the supported data rates in 802.11-2012
% it can have the value of 6,9,12,18,24,36,48, or 54
% scheme and convolutional code of a given rate as defined in 802.11-2012
% eccRate is the rate of the convolutional code
% hConvEnc is the convolutional encoder obj
% R1R4 are the bits that get placed in the SIGNAL frame
% NBPSC is the number of bits per subcarrier
% NCBPS is the number of coded bits per OFDM symbol
% NDBPS is the number of data bits per OFDM symbol
% NMSPOS is the number of modulated symbols per OFDM symbol, constant 48

k = 7; %constraint length of both convolutional codes
trel = poly2trellis(k ,[133 171]);
punc34 = [1 1 1 0 0 1].';
punc23 = [1 1 1 0].';

% hConvEnc = comm.ConvolutionalEncoder(trel,'TerminationMethod','Truncated',...
%     'PuncturePatternSource','Property','PuncturePattern',punc34);

hCE12 = comm.ConvolutionalEncoder(trel,'TerminationMethod','Truncated');
hCE23 = comm.ConvolutionalEncoder(trel,'TerminationMethod','Truncated',...
    'PuncturePatternSource','Property','PuncturePattern',punc23);
hCE34 = comm.ConvolutionalEncoder(trel,'TerminationMethod','Truncated',...
    'PuncturePatternSource','Property','PuncturePattern',punc34);

hVD12 = comm.ViterbiDecoder(trel,...
        'InputFormat','Hard','TerminationMethod','Truncated');
hVD23 = comm.ViterbiDecoder(trel,...
        'InputFormat','Hard','TerminationMethod','Truncated',...
        'PuncturePatternSource','Property','PuncturePattern',punc23);
hVD34 = comm.ViterbiDecoder(trel,...
        'InputFormat','Hard','TerminationMethod','Truncated',...
        'PuncturePatternSource','Property','PuncturePattern',punc34);    

if isequal(RATE,6)
    rateStruct.eccRate = 1/2;
    rateStruct.hConvEnc = hCE12;
    rateStruct.hVitDec = hVD12;
    rateStruct.hVitDec.TracebackDepth = 5*(k-1);
    rateStruct.hMod = comm.BPSKModulator('PhaseOffset',pi);
    rateStruct.hDemod = comm.BPSKDemodulator('PhaseOffset',pi);
    rateStruct.R1R4 = [1 1 0 1];
    rateStruct.NBPSC = 1;
    rateStruct.NCBPS = 48;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '6 Mbps';
elseif isequal(RATE,9)
    rateStruct.eccRate = 3/4;
    rateStruct.hConvEnc = hCE34;
    rateStruct.hVitDec = hVD34;
    rateStruct.hVitDec.TracebackDepth = 10*(k-1);
    rateStruct.hMod = comm.BPSKModulator('PhaseOffset',pi);
    rateStruct.hDemod = comm.BPSKDemodulator('PhaseOffset',pi);
    rateStruct.R1R4 = [1 1 1 1];
    rateStruct.NBPSC = 1;
    rateStruct.NCBPS = 48;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '9 Mbps';
elseif isequal(RATE,12)
    rateStruct.eccRate = 1/2;
    rateStruct.hConvEnc = hCE12;
    rateStruct.hVitDec = hVD12;
    rateStruct.hVitDec.TracebackDepth = 5*(k-1);
    rateStruct.hMod = comm.QPSKModulator;
    rateStruct.hDemod = comm.QPSKDemodulator;
    rateStruct.R1R4 = [0 1 0 1];
    rateStruct.NBPSC = 2;
    rateStruct.NCBPS = 96;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '12 Mbps';
elseif isequal(RATE,18)
    rateStruct.eccRate = 3/4;
    rateStruct.hConvEnc = hCE34;
    rateStruct.hVitDec = hVD34;
    rateStruct.hVitDec.TracebackDepth = 10*(k-1);
    rateStruct.hMod = comm.QPSKModulator;
    rateStruct.hDemod = comm.QPSKDemodulator;
    rateStruct.R1R4 = [0 1 1 1];
    rateStruct.NBPSC = 2;
    rateStruct.NCBPS = 96;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '18 Mbps';
elseif isequal(RATE,24)
    rateStruct.eccRate = 1/2;
    rateStruct.hConvEnc = hCE12;
    rateStruct.hVitDec = hVD12;
    rateStruct.hVitDec.TracebackDepth = 5*(k-1);
    const = (-3:2:3)+1j*(3:-2:-3).';
    const = const(:)/sqrt(10);
    rateStruct.hMod = comm.GeneralQAMModulator(const);
    rateStruct.hDemod = comm.GeneralQAMDemodulator(const);
    rateStruct.R1R4 = [1 0 0 1];
    rateStruct.NBPSC = 4;
    rateStruct.NCBPS = 192;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '24 Mbps';
elseif isequal(RATE,36)
    rateStruct.eccRate = 3/4;
    rateStruct.hConvEnc = hCE34;
    rateStruct.hVitDec = hVD34;
    rateStruct.hVitDec.TracebackDepth = 10*(k-1);
    const = (-3:2:3)+1j*(3:-2:-3).';
    const = const(:)/sqrt(10);
    rateStruct.hMod = comm.GeneralQAMModulator(const);
    rateStruct.hDemod = comm.GeneralQAMDemodulator(const);
    rateStruct.R1R4 = [1 0 1 1];
    rateStruct.NBPSC = 4;
    rateStruct.NCBPS = 192;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '36 Mbps';
elseif isequal(RATE,48)
    rateStruct.eccRate = 2/3;
    rateStruct.hConvEnc = hCE23;
    rateStruct.hVitDec = hVD23;
    rateStruct.hVitDec.TracebackDepth = 8*(k-1);
    const = (-7:2:7)+1j*(7:-2:-7).';
    const = const(:)/sqrt(42); 
    rateStruct.hMod = comm.GeneralQAMModulator(const);
    rateStruct.hDemod = comm.GeneralQAMDemodulator(const);
    rateStruct.R1R4 = [0 0 0 1];
    rateStruct.NBPSC = 6;
    rateStruct.NCBPS = 288;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '48 Mbps';
elseif isequal(RATE, 54)
    rateStruct.eccRate = 3/4;
    rateStruct.hConvEnc = hCE34;
    rateStruct.hVitDec = hVD34;
    rateStruct.hVitDec.TracebackDepth = 10*(k-1);
    const = (-7:2:7)+1j*(7:-2:-7).';
    const = const(:)/sqrt(42);
    rateStruct.hMod = comm.GeneralQAMModulator(const);
    rateStruct.hDemod = comm.GeneralQAMDemodulator(const);
    rateStruct.R1R4 = [0 0 1 1];
    rateStruct.NBPSC = 6;
    rateStruct.NCBPS = 288;
    rateStruct.NDBPS = rateStruct.NCBPS*rateStruct.eccRate;
    rateStruct.ttl = '54 Mbps';
end
rateStruct.NMSPOS = rateStruct.NCBPS/rateStruct.NBPSC;
rateStruct.hOFDMmod = comm.OFDMModulator;
rateStruct.hOFDMdemod = comm.OFDMDemodulator;
end