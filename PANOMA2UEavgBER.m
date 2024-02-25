function [BER_avg,BERth_NU_PANOMA,BERth_FU_PANOMA,sumRate_P] = PANOMA2UEavgBER(a,PL,SNR_dB)
N=2;
Va1 = [a(1) 0.5];               %NU power coefficients vector
Va2 = [a(2) 0.5];               %FU power coefficients vector
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
SNR = 10.^(SNR_dB./10);
%Lookup table for PANOMA
MLD_LUT_P = [(1*SVa1(2)+1*SVa2(2)) ...
    (1*SVa1(1)-1*SVa2(1)) ...
    (-1*SVa1(1)+1*SVa2(1)) ...
    (-1*SVa1(2)-1*SVa2(2))];
normFactor = (MLD_LUT_P*MLD_LUT_P')/4;
MLD_LUT_P = MLD_LUT_P/sqrt(normFactor);
% Theoretical BER calculations
xA = MLD_LUT_P(1);
xB = MLD_LUT_P(3);
Sd1 = sqrt(PL(1)./(2*N*SNR));
Sd2 = sqrt(PL(2)./(2*N*SNR));
gamma1PANOMA_u1 = (((xA - xB)/2)./(Sd1)).^2;
gamma2PANOMA_u1 = ((xB)./(Sd1)).^2;
gamma1cPANOMA_u1 = (((3*xA + xB)/2)./(Sd1)).^2;
gamma2cPANOMA_u1 = ((xA)./(Sd1)).^2;
gamma3cPANOMA_u1 = (((xA + 3*xB)/2)./(Sd1)).^2;

gamma1PANOMA_u2 = ((xA)./(Sd2)).^2;
gamma2PANOMA_u2 = ((xB)./(Sd2)).^2;

BERth_NU_PANOMA = 0.5*(1 - sqrt(0.5*gamma1PANOMA_u1./(1 + 0.5*gamma1PANOMA_u1))) + 0.25*(1 - sqrt(0.5*gamma2PANOMA_u1./(1 + 0.5*gamma2PANOMA_u1)))...
    - 0.25*(1 - sqrt(0.5*gamma2cPANOMA_u1./(1 + 0.5*gamma2cPANOMA_u1))) - 0.25*(1 - sqrt(0.5*gamma3cPANOMA_u1./(1 + 0.5*gamma3cPANOMA_u1)))...
    + 0.25*(1 - sqrt(0.5*gamma1cPANOMA_u1./(1 + 0.5*gamma1cPANOMA_u1)));
BERth_FU_PANOMA = 0.25*(1 - sqrt(0.5*gamma1PANOMA_u2./(1 + 0.5*gamma1PANOMA_u2))) + 0.25*(1 - sqrt(0.5*gamma2PANOMA_u2./(1 + 0.5*gamma2PANOMA_u2)));
BER_avg = 0.5*(BERth_NU_PANOMA+BERth_FU_PANOMA);
c1=-((1-BERth_NU_PANOMA).*log2(1-BERth_NU_PANOMA)) ...
        - ((BERth_NU_PANOMA).*log2(BERth_NU_PANOMA));
c1 = 1 - c1;
c2=-((1-BERth_FU_PANOMA).*log2(1-BERth_FU_PANOMA)) ...
        - ((BERth_FU_PANOMA).*log2(BERth_FU_PANOMA));
c2 = 1 - c2;
sumRate_P = c1+c2;