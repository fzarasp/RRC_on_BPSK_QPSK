close all
clear all
clc
%% initial parameters
N = 1e5;
SPS = 5 ;
span = 10;
rolloff = 0.2 ;
N0 = 1;
Pt = 1 : 0.1 :16;
EbN = 10 * log10(Pt);
Rs = 1e3;
Fs = Rs  * SPS ;
BW  = (1 + rolloff) * Rs;
Period = 1;

%% Part 1 :  base band signals 
rbpsk = BB_signal(2 , N);
rqpsk = BB_signal(4 , N);
bpsks = 6 * pskmod( rbpsk, 2);
qpsks = 6 * pskmod(rqpsk,4,pi/4)/log2(4);

%% Part 2 : pulse shaping
h = rcosdesign(rolloff,span,SPS);
sbpsks = upfirdn(bpsks,h,SPS);
sqpsks = upfirdn(qpsks,h,SPS);
%% Part 2 : ploting  real signals
w = 1 ;
figure(w);
w = w + 1 ;
plot(real(sbpsks));
xlim([0 500]);
title('BPSK after pule shaping');

figure(w);
w = w + 1 ;
plot(real(sqpsks));
xlim([0 500]);
title('QPSK after pule shaping');

%% Part 3 : at reciever
% add noise 
nbpsk =  sbpsks + noise_generator(numel(sbpsks) , N0);
nqpsk =  sqpsks + noise_generator(numel(sqpsks) , N0);

% for Eb/N0 = inf
n0bpsk = sbpsks ;
n0qpsk = sqpsks ;

% down sample 
fbpsk = upfirdn(nbpsk,h,1 ,SPS);
fqpsk = upfirdn(nqpsk,h,1 ,SPS);
dsbpsk = fbpsk(span+1:end-span);
dsqpsk = fqpsk(span+1:end-span);

%demodulation
dbpsk = pskdemod(dsbpsk, 2);
dqpsk = pskdemod(dsqpsk, 4 , pi/4);

figure(w);
w = w +1 ;
% P.S.D
[pbpsk,f] = pwelch(nbpsk(1:2048),[],[],[],Fs,'centered');
[pqpsk,f] = pwelch(nqpsk(1:2048),[],[],[],Fs,'centered');
plot(f , pbpsk , 'r');
hold on;
plot(f , pqpsk , 'b');
hold off;
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend({'PSD of BPSK', 'PSD of QPSK'} ,'Location','northwest');
maxp = max(pbpsk);
minp = min(pqpsk);
rectangle('Position',[-BW/2 minp BW maxp])

%% Part 4 : eye diagram
Offset = 0;
eyediagram(nbpsk(1:1000),SPS,Period,Offset)
title('noisy BPSK')

eyediagram(n0bpsk(1:1000),SPS,Period,Offset)
title('noise free BPSK')

eyediagram(nqpsk(1:1000),SPS,Period,Offset)
title('noisy QPSK')

eyediagram(n0qpsk(1:1000),SPS,Period,Offset)
title('noise free QPSK')
w = w + 4;

%% Part 5 : BER analysis
scatterplot(dsbpsk)
title('BPSK');
scatterplot(dsqpsk)
title('QPSK');
w = w + 2;
inx = 1 ;
berror = zeros(1 , numel(Pt));
qerror = zeros(1 , numel(Pt));

for p = Pt
    % random signal generator
    rbpskt = BB_signal(2 , N);
    rqpskt = BB_signal(4 , N);
    
    %modulation
    bpskst =  p * pskmod( rbpskt, 2);
    qpskst =  p *pskmod(rqpskt,4,pi/4)/log2(4);
    
    % pulse shaping at transmiter
    sbpskst = upfirdn(bpskst,h,SPS);
    sqpskst = upfirdn(qpskst,h,SPS);
    
    % adding noise
    nbpskt = sbpskst + noise_generator(numel(sbpskst) , N0);
    nqpskt = sqpskst + noise_generator(numel(sqpskst) , N0);
    
    % down sampling
    fbpskt = upfirdn(nbpskt,h,1 ,SPS);
    fqpskt = upfirdn(nqpskt,h,1 ,SPS);
    dsbpskt = fbpskt(span+1:end-span);
    dsqpskt = fqpskt(span+1:end-span);
    
    % demodulation
    dbpskt = pskdemod(dsbpskt, 2);
    dqpskt = pskdemod(dsqpskt, 4 , pi/4);
    
    % SER analysis
    berror(inx) = symerr(rbpskt , dbpskt)/N;
    qerror(inx) = 0.5 * symerr(rqpskt , dqpskt)/N;
    
    inx = inx + 1 ;
end
tberror = qfunc(sqrt(2 * Pt));
% ps(QPSK) = 2 pb(QPSK)
tqerror = (1 - (1- qfunc(sqrt(Pt))).^2)/2;
figure(w);
w = w + 1 ;
semilogy(EbN , berror , 'r*-');
hold on;
semilogy(EbN , tberror , 'y-');
hold on;
semilogy(EbN , qerror , 'b*-');
hold on;
semilogy(EbN , tqerror , 'c-');
hold off;
xlabel('Eb / N0 (dB)')
ylabel('BER')
title('BER analysis for BPSK & QPSK ');
legend({'BPSK : simulation', 'BPSK : theorical ','QPSK : simulation','QPSK : theorical'} ,'Location','southwest');
grid on;