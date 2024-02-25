function [BER_avg,BERth_NU_PANOMA,BERth_MU_PANOMA,BERth_FU_PANOMA,sumRate] = PANOMA3UEavgBER(a,PL,SNR_dB)
N = 3;
Va1 = [0.15/2 a(1) 0.050 1/3];    %NU power coefficients vector
Va2 = [0.15/2 a(2) 0.475 1/3];    %MU power coefficients vector
Va3 = [0.8500 a(3) 0.475 1/3];    %FU power coefficients vector
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
SVa3 = sqrt(Va3);
SNR = 10.^(SNR_dB./10);

%Lookup table for PANOMA
MLD_LUT_P = [(1*SVa1(4)+1*SVa2(4)+1*SVa3(4)) ...
             (1*SVa1(1)+1*SVa2(1)-1*SVa3(1)) ...
             (1*SVa1(2)-1*SVa2(2)+1*SVa3(2)) ...
             (1*SVa1(3)-1*SVa2(3)-1*SVa3(3)) ...
             (-1*SVa1(3)+1*SVa2(3)+1*SVa3(3)) ...
             (-1*SVa1(2)+1*SVa2(2)-1*SVa3(2)) ...
             (-1*SVa1(1)-1*SVa2(1)+1*SVa3(1)) ...
             (-1*SVa1(4)-1*SVa2(4)-1*SVa3(4))];
normFactor = (MLD_LUT_P*MLD_LUT_P')/8;
MLD_LUT_P = MLD_LUT_P/sqrt(normFactor);

% Theoretical BER calculations
Sd1 = sqrt(PL(1)./(2*N*SNR));
Sd2 = sqrt(PL(2)./(2*N*SNR));
Sd3 = sqrt(PL(3)./(2*N*SNR));
xA = MLD_LUT_P(1);
xB = MLD_LUT_P(5);
xC = MLD_LUT_P(3);
xD = MLD_LUT_P(7);

gamma1PANOMA_u1 = (((xA - xB)/2)./(Sd1)).^2;
gamma2PANOMA_u1 = (((xB - xC)/2)./(Sd1)).^2;
gamma3PANOMA_u1 = (((xC - xD)/2)./(Sd1)).^2;
gamma4PANOMA_u1 = ((xD)./(Sd1)).^2;
gamma1cPANOMA_u1 = (((2*xA - xC - xD)/2)./(Sd1)).^2;
gamma2cPANOMA_u1 = (((2*xA + xC + xD)/2)./(Sd1)).^2;
gamma3cPANOMA_u1 = (((3*xA + xB)/2)./(Sd1)).^2;
gamma4cPANOMA_u1 = ((xB)./(Sd1)).^2;
gamma5cPANOMA_u1 = (((3*xB + xC)/2)./(Sd1)).^2;
gamma6cPANOMA_u1 = (((3*xC + xD)/2)./(Sd1)).^2;
gamma7cPANOMA_u1 = (((2*xC + xA + xB)/2)./(Sd1)).^2;
gamma8cPANOMA_u1 = (((xA + xB - 2*xD)/2)./(Sd1)).^2;
gamma9cPANOMA_u1 = (((2*xD + xB + xC)/2)./(Sd1)).^2;
gamma10cPANOMA_u1 = (((2*xA - xB - xC)/2)./(Sd1)).^2;
gamma11cPANOMA_u1 = ((xA)./(Sd1)).^2;
gamma12cPANOMA_u1 = (((2*xA + xB + xC)/2)./(Sd1)).^2;
gamma13cPANOMA_u1 = (((2*xB - xC - xD)/2)./(Sd1)).^2;
gamma14cPANOMA_u1 = (((2*xB + xC + xD)/2)./(Sd1)).^2;
gamma15cPANOMA_u1 = (((3*xB + xA)/2)./(Sd1)).^2;
gamma16cPANOMA_u1 = (((xA + xB - 2*xC)/2)./(Sd1)).^2;
gamma17cPANOMA_u1 = ((xC)./(Sd1)).^2;
gamma18cPANOMA_u1 = (((3*xC + xB)/2)./(Sd1)).^2;
gamma19cPANOMA_u1 = (((xB + xC - 2*xD)/2)./(Sd1)).^2;
gamma20cPANOMA_u1 = (((3*xD + xC)/2)./(Sd1)).^2;
gamma21cPANOMA_u1 = (((2*xD + xA + xB)/2)./(Sd1)).^2;

gamma1PANOMA_u2 = (((2*xA - xB - xC)/2)./(Sd2)).^2;
gamma2PANOMA_u2 = (((xB - xC)/2)./(Sd2)).^2;
gamma3PANOMA_u2 = ((xC)./(Sd2)).^2;
gamma4PANOMA_u2 = (((xB + xC - 2*xD)/2)./(Sd2)).^2;
gamma5PANOMA_u2 = ((xD)./(Sd2)).^2;
gamma1cPANOMA_u2 = (((3*xB + xC)/2)./(Sd2)).^2;
gamma2cPANOMA_u2 = (((2*xA + xB + xC)/2)./(Sd2)).^2;
gamma3cPANOMA_u2 = ((xA)./(Sd2)).^2;
gamma4cPANOMA_u2 = ((xB)./(Sd2)).^2;
gamma5cPANOMA_u2 = (((3*xC + xB)/2)./(Sd2)).^2;
gamma6cPANOMA_u2 = (((2*xD + xB + xC)/2)./(Sd2)).^2;

gamma1PANOMA_u3 = ((xA)./(Sd3)).^2;
gamma2PANOMA_u3 = ((xB)./(Sd3)).^2;
gamma3PANOMA_u3 = ((xC)./(Sd3)).^2;
gamma4PANOMA_u3 = ((xD)./(Sd3)).^2;

BERth_NU_PANOMA = 0.25*(1 - sqrt(0.5*gamma1PANOMA_u1./(1 + 0.5*gamma1PANOMA_u1))) + 0.25*(1 - sqrt(0.5*gamma2PANOMA_u1./(1 + 0.5*gamma2PANOMA_u1))) ...
    + 0.25*(1 - sqrt(0.5*gamma3PANOMA_u1./(1 + 0.5*gamma3PANOMA_u1))) + 0.125*(1 - sqrt(0.5*gamma4PANOMA_u1./(1 + 0.5*gamma4PANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma1cPANOMA_u1./(1 + 0.5*gamma1cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma10cPANOMA_u1./(1 + 0.5*gamma10cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma2cPANOMA_u1./(1 + 0.5*gamma2cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma11cPANOMA_u1./(1 + 0.5*gamma11cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma3cPANOMA_u1./(1 + 0.5*gamma3cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma12cPANOMA_u1./(1 + 0.5*gamma12cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma4cPANOMA_u1./(1 + 0.5*gamma4cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma13cPANOMA_u1./(1 + 0.5*gamma13cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma5cPANOMA_u1./(1 + 0.5*gamma5cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma14cPANOMA_u1./(1 + 0.5*gamma14cPANOMA_u1))) ...
    - 0.125*(1 - sqrt(0.5*gamma15cPANOMA_u1./(1 + 0.5*gamma15cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma16cPANOMA_u1./(1 + 0.5*gamma16cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma6cPANOMA_u1./(1 + 0.5*gamma6cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma17cPANOMA_u1./(1 + 0.5*gamma17cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma7cPANOMA_u1./(1 + 0.5*gamma7cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma18cPANOMA_u1./(1 + 0.5*gamma18cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma8cPANOMA_u1./(1 + 0.5*gamma8cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma19cPANOMA_u1./(1 + 0.5*gamma19cPANOMA_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma9cPANOMA_u1./(1 + 0.5*gamma9cPANOMA_u1))) - 0.125*(1 - sqrt(0.5*gamma20cPANOMA_u1./(1 + 0.5*gamma20cPANOMA_u1))) ...
    - 0.125*(1 - sqrt(0.5*gamma21cPANOMA_u1./(1 + 0.5*gamma21cPANOMA_u1)));

BERth_MU_PANOMA = 0.125*(1 - sqrt(0.5*gamma1PANOMA_u2./(1 + 0.5*gamma1PANOMA_u2))) ...
    + 0.25*(1 - sqrt(0.5*gamma2PANOMA_u2./(1 + 0.5*gamma2PANOMA_u2))) + 0.125*(1 - sqrt(0.5*gamma3PANOMA_u2./(1 + 0.5*gamma3PANOMA_u2))) ...
    + 0.125*(1 - sqrt(0.5*gamma4PANOMA_u2./(1 + 0.5*gamma4PANOMA_u2))) + + 0.125*(1 - sqrt(0.5*gamma5PANOMA_u2./(1 + 0.5*gamma5PANOMA_u2))) ...
    + 0.125*(1 - sqrt(0.5*gamma1cPANOMA_u2./(1 + 0.5*gamma1cPANOMA_u2))) - 0.125*(1 - sqrt(0.5*gamma3cPANOMA_u2./(1 + 0.5*gamma3cPANOMA_u2))) ...
     + 0.125*(1 - sqrt(0.5*gamma2cPANOMA_u2./(1 + 0.5*gamma2cPANOMA_u2))) - 0.125*(1 - sqrt(0.5*gamma4cPANOMA_u2./(1 + 0.5*gamma4cPANOMA_u2)))...
     - 0.125*(1 - sqrt(0.5*gamma5cPANOMA_u2./(1 + 0.5*gamma5cPANOMA_u2)))...
     - 0.125*(1 - sqrt(0.5*gamma6cPANOMA_u2./(1 + 0.5*gamma6cPANOMA_u2)));

BERth_FU_PANOMA = 0.125*(1 - sqrt(0.5*gamma1PANOMA_u3./(1 + 0.5*gamma1PANOMA_u3))) ...
    + 0.125*(1 - sqrt(0.5*gamma2PANOMA_u3./(1 + 0.5*gamma2PANOMA_u3))) + 0.125*(1 - sqrt(0.5*gamma3PANOMA_u3./(1 + 0.5*gamma3PANOMA_u3))) ...
    + 0.125*(1 - sqrt(0.5*gamma4PANOMA_u3./(1 + 0.5*gamma4PANOMA_u3))); 


BER_avg = 1/3*(BERth_NU_PANOMA+BERth_MU_PANOMA+BERth_FU_PANOMA);
c1=-((1-BERth_NU_PANOMA).*log2(1-BERth_NU_PANOMA)) ...
        - ((BERth_NU_PANOMA).*log2(BERth_NU_PANOMA));
c1 = 1 - c1;
c2=-((1-BERth_MU_PANOMA).*log2(1-BERth_MU_PANOMA)) ...
        - ((BERth_MU_PANOMA).*log2(BERth_MU_PANOMA));
c2 = 1 - c2;
c3=-((1-BERth_FU_PANOMA).*log2(1-BERth_FU_PANOMA)) ...
        - ((BERth_FU_PANOMA).*log2(BERth_FU_PANOMA));
c3 = 1 - c3;

sumRate = c1+c2+c3;