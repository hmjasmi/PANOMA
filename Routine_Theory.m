%%
clc
clear
close all
a = [0.05,0.95];
PL = 10.^([0 6]./10);
SNR_dB = -5:1:35;
for i_snr = 1:1:length(SNR_dB)
    [BER_avg(i_snr),BERth_NU(i_snr),BERth_FU(i_snr),sumRate(i_snr)] = CNOMA2UEavgBER(a,PL,SNR_dB(i_snr));
    [BER_avg_PANOMA(i_snr),BERth_NU_PANOMA(i_snr),BERth_FU_PANOMA(i_snr),sumRate_PANOMA(i_snr)] = PANOMA2UEavgBER(a,PL,SNR_dB(i_snr));
end
%%
clc
clear
close all
a1 = 0:0.001:0.5;
a2 = 1-a1;
a = [a1.' a2.'];
PL = 10.^([0 6]./10);
SNR_dB = [15,25,35];
for i_snr = 1:1:length(SNR_dB)
    for i_alpha = 1:1:size(a,1)
        [BERth_avg(i_snr,i_alpha),BER1th(i_snr,i_alpha),BER2th(i_snr,i_alpha),sumRate(i_snr,i_alpha)] = CNOMA2UEavgBER(a(i_alpha,:),PL,SNR_dB(i_snr));
        [BERth_avg_PANOMA(i_snr,i_alpha),BER1th_PANOMA(i_snr,i_alpha),BER2th_PANOMA(i_snr,i_alpha),sumRate_PANOMA(i_snr,i_alpha)] = PANOMA2UEavgBER(a(i_alpha,:),PL,SNR_dB(i_snr));
    end
end

%%
clc
clear
close all
a = [0.05 0.15 0.80;...
     0.10 0.20 0.70];
PL = 10.^([0 6 12]./10);
SNR_dB = -5:1:35;
for i_snr = 1:1:length(SNR_dB)
    for i_alpha = 1:1:size(a,1)
        [BERth_avg(i_alpha,i_snr),BER1th(i_alpha,i_snr),BER2th(i_alpha,i_snr),BER3th(i_alpha,i_snr),sumRate(i_alpha,i_snr)] = CNOMA3UEavgBER(a(i_alpha,:),PL,SNR_dB(i_snr));
        [BERth_avg_PANOMA(i_alpha,i_snr),BER1th_PANOMA(i_alpha,i_snr),BER2th_PANOMA(i_alpha,i_snr),BER3th_PANOMA(i_alpha,i_snr),sumRate_PANOMA(i_alpha,i_snr)] = PANOMA3UEavgBER(a(i_alpha,:),PL,SNR_dB(i_snr));
    end
end

