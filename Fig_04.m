clc
clear 
close all
BERth_avg = [0.38147963956440417 0.36768281756329624 0.35239160817530479 ...
           0.33565769631191567 0.31760228500827947 0.29840850234578109 ...
           0.27830811019882346 0.25756561324625793 0.23646283771287319 ...
           0.21528612058677055 0.19431684851849393 0.17382476042856518 ...
           0.15406261667673679 0.13526073400621366 0.11762044088979005 ...
           0.1013064916568788 0.086439524514177177 0.073090329180620014 ...
           0.061277673390942144 0.050970689184895 0.042095661503406934 ...
           0.034545998673074452 0.02819359658306303 0.022899842247475863 ...
           0.018524973272628031 0.014935141198752402 0.012007087468279568 ...
           0.0096307153666203849 0.0077100190672552693 0.00616285940134198 ...
           0.0049200171399426534];
BERth_avg_PANOMA = [0.375552664215781 0.36101301651916839 0.34480528529253263 ...
                  0.32697860137618573 0.30766769649462938 0.28708705355732828 ...
                  0.26551794947500756 0.24329132385574415 0.22076934176763882 ...
                  0.19832750684318806 0.17633779879821537 0.15515235011546685 ...
                  0.1350872054462301 0.11640667794636644 0.099310058968039142 ...
                  0.083923076756304546 0.070296028736010777 0.058409099162809633 ...
                  0.048183702955456764 0.039497500007442843 0.03220037862886363 ...
                  0.026129134051361871 0.02111941277258314 0.017014375656784514 ...
                  0.013670208112059656 0.010958996234534751 0.0087696299112877962 ...
                  0.007007370368463911 0.0055926120075789482 0.0044592355386600807 ...
                  0.0035528252865042603];
BER2th = [0.32055838656025515 0.30231256046267629 0.28300950091344684 ...
            0.26281458234602573 0.24194884976318284 0.22068420687791734 ...
            0.1993325288603234 0.17822921917381779 0.15771282293491293 ...
            0.13810310002064932 0.11968026745971866 0.1026679532557902 ...
            0.087221867257717856 0.073425369941550439 0.061292074173872696 ...
            0.050774504247887015 0.041776932793517713 0.034170083616831604 ...
            0.027805519599503059 0.022528101058188144 0.01818564522305377 ...
            0.014635598997124355 0.011749015117108602 0.0094123688524462568 ...
            0.00752781348357337 0.0060124178556764662 0.0047968211649842 ...
            0.0038236219312602071 0.0030457132505500861 0.0024246939271041212 ...
            0.001929425509166055];

BER2th_PANOMA = [0.32468460030695734 0.30698713475292649 0.28828473726530024 ...
                   0.26873154671408717 0.24852860193961451 0.22791808671761724 ...
                   0.20717321140827349 0.18658484789847871 0.16644656401059627 ...
                   0.14703966846015237 0.1286194042992295 0.11140292557967346 ...
                   0.095559529991334824 0.081203818632216984 0.068392664521673974 ...
                   0.057126689841123152 0.047356271546199313 0.038991190271346621 ...
                   0.031912338938348106 0.02598370021201768 0.021063087870904051 ...
                   0.01701071898674808 0.013695290047211844 0.010997698308477161 ...
                   0.008812819519692372 0.0070498469651252837 0.0056316723188548357 ...
                   0.00449370445505945 0.0035824208349132802 0.0028538520736524353 ...
                   0.0022721244120070816];
BER1th = [0.44240089256855319 0.43305307466391618 0.42177371543716269 ...
            0.4085008102778056 0.39325572025337613 0.37613279781364489 ...
            0.35728369153732353 0.336902007318698 0.31521285249083342 0.2924691411528918 ...
            0.2689534295772692 0.24498156760134016 0.22090336609575573 ...
            0.19709609807087689 0.17394880760570741 0.15183847906587059 ...
            0.13110211623483664 0.11201057474440843 0.09474982718238123 ...
            0.07941327731160186 0.0660056777837601 0.054456398349024548 ...
            0.044638178049017457 0.036387315642505469 0.029522133061682693 ...
            0.023857864541828339 0.019217353771574935 0.015437808801980563 ...
            0.012374324883960452 0.0099010248755798391 0.0079106087707192518 ...
            ];
BER1th_PANOMA = [0.42642072812460463 0.4150388982854103 0.401325833319765 ...
                   0.38522565603828435 0.36680679104964425 0.34625602039703929 ...
                   0.32386268754174163 0.29999779981300956 0.27509211952468138 ...
                   0.24961534522622375 0.22405619329720125 0.19890177465126022 ...
                   0.1746148809011254 0.1516095372605159 0.13022745341440431 ...
                   0.11071946367148594 0.09323578592582224 0.077827008054272645 ...
                   0.064455066972565422 0.053011299802868006 0.043337669386823208 ...
                   0.035247549115975663 0.028543535497954436 0.023031053005091867 ...
                   0.018527596704426941 0.014868145503944219 0.011907587503720757 ...
                   0.0095210362818683725 0.0076028031802446161 0.0060646190036677261 ...
                   0.0048335261610014391];
PL = [1 3.9810717055349722];
SNR_dB = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 ...
          19 20 21 22 23 24 25];
a = [0.05 0.95];
sumRate = [0.10460637194047528 0.12887151060637114 0.15821033319484767 ...
           0.19332946312398391 0.23486537297154553 0.28332023376646231 ...
           0.3390035341645119 0.40198888502071695 0.47209084414927982 0.54886116243227334 ...
           0.63160003324652991 0.71937703021110311 0.81105825010756516 ...
           0.90533950263355556 1.0007884925511175 1.0959000300197632 1.1891661165218035 ...
           1.2791575010108107 1.3646071755143732 1.4444825323550234 1.5180338117870802 ...
           1.5848117955027112 1.6446549686694483 1.6976525624100991 1.7440931178120556 ...
           1.7844082964426158 1.8191196022947225 1.8487927951919407 1.8740020815152771 ...
           1.8953041757291715 1.913221130451662];
sumRate_PANOMA = [0.10627489263355305 0.1312643162446715 0.1617839156359473 ...
                  0.19871213097574913 0.24287365632151947 0.29495578615355589 ...
                  0.35542279235866536 0.42443948412750843 0.50181191129323266 ...
                  0.58694988669354542 0.67885516416059155 0.77614064418958106 ...
                  0.87708696103727091 0.97973966050697625 1.0820417948564258 ...
                  1.1819859864477562 1.2777621640037131 1.3678765875684769 ...
                  1.4512250067996602 1.5271145199937397 1.5952400504691528 ...
                  1.6556286726282856 1.7085671990223021 1.7545265483411803 ...
                  1.7940924359910304 1.8279076529395288 1.8566277050389224 ...
                  1.880889250795138 1.9012895194110788 1.9183744337314079 ...
                  1.9326332021096175];
BER1 = ...
  [0.44245790014033287 0.37630603731320894 0.26872288029873231 0.15224662792790691 ...
   0.06600274874083753 0.023844009061636463 0.0079615099199669328;
   0.4214891158696138 0.33734035990713362 0.21217898023673254 0.10108984532551558 ...
   0.038327450366832112 0.013168445688514371 0.0041830850147166173;
   0.395542665899447 0.30280596981343394 0.17285168945270182 0.073473400921996931 ...
   0.026074288085706381 0.00865148157839474 0.0027972823423921922;
   0.37859149844500517 0.29537362479625068 0.17846435303548988 0.082611521503261662 ...
   0.031106510894963683 0.010635849963833453 0.0034253271239095872;
   0.36421764010786628 0.30063522184092717 0.21662888727870908 0.13549428793570689 ...
   0.069296644011186626 0.028337249292502358 0.0100494717100943];

BER1_PANOMA = ...
  [0.42650123457921807 0.34644402747824171 0.22375914996950011 0.11110681714394285 ...
   0.043316392070359769 0.014870262932456892 0.0049078742654191151;
   0.41300765976613413 0.3239590243032523 0.19494830850563832 0.088134628092906356 ...
   0.032257991431695228 0.010929651067829773 0.0034506395395348683;
   0.39308150222832589 0.30075589852200491 0.1715703135156216 0.072919001728327576 ...
   0.025860278382405393 0.0085850234665884439 0.0027782459475135081;
   0.3779925160666131 0.29539271327428906 0.17936783960720132 0.08365961696794344 ...
   0.031691404778650741 0.010862047126509579 0.0035002747907506974;
   0.36414303098156342 0.30068483105056315 0.21688826141412862 0.13592475525081582 ...
   0.069694559351468829 0.028553811070629764 0.010137205792647357];

BER2 = ...
  [0.32069747267509108 0.22074317043943187 0.11977814240619197 0.051005533106556311 ...
   0.018234470468431772 0.0059953445988513371 0.0019480403815320616;
   0.3261365691281029 0.23064509055803148 0.13258825595581347 0.059988055248149173 ...
   0.022464899075336414 0.0076092715107616308 0.002409965925113583;
   0.33837509604134652 0.25261663190289368 0.16194823621421262 0.085880130399565338 ...
   0.036092041151529496 0.012937743332522225 0.0042590222615924616;
   0.35066617486275048 0.27471416241112528 0.19790267365775446 0.12618996999176668 ...
   0.065076788285705719 0.026588322830590566 0.009338614704617651;
   0.36339472618424606 0.29846272908256971 0.23791251424995249 0.1858109170886097 ...
   0.13056628352905492 0.073965664906116974 0.03208434617717941];

BER2_PANOMA = ...
  [0.32482865682114392 0.22799619313768954 0.12874592501358328 0.057357439016869946 ...
   0.021108705679314402 0.0070276328245572513 0.0022937423541921527;
   0.329324865792114 0.23632142580358065 0.14010849130502898 0.065938373955420154 ...
   0.025390566406445311 0.00868609083803054 0.0027698605587981374;
   0.33989425764414116 0.25536847168842769 0.166070488098373 0.09002139263702455 ...
   0.038619350435498551 0.013999328335572214 0.0046310522714924285;
   0.35115349094669684 0.27561514899117 0.19938699162669457 0.12813498955003483 ...
   0.066745220224265919 0.0274919396102013 0.0096934572718090934;
   0.36346048117339608 0.29858489013369954 0.23812181042729857 0.18615255657481142 ...
   0.13103576113079624 0.074408658221139262 0.032345569264769114];

BER_su = ...
  [0.25528841987193374 0.14666401632827891 0.0640175209832634 0.023475468623437922 ...
   0.007672031718227606 0.0024715282199059335 0.00084166386112046289;
   0.25524930541898194 0.14627917386108713 0.064056557311475623 0.023284922383592054 ...
   0.0076193235605881314 0.0024701219745934182 0.00078624737917540278;
   0.25484698905170317 0.14611050254832483 0.064367441691861024 0.023194740392532026 ...
   0.0076732296309012307 0.002427465866780444 0.00080033587388042039;
   0.25480464544284853 0.14650805330648897 0.064203666196112674 0.023266953693487687 ...
   0.0076714588034706554 0.0024809031886560379 0.000788278622404592;
   0.2545645420765264 0.14617349713000957 0.064081479103402991 0.023267917231942561 ...
   0.0076489067953106823 0.0024616324195586015 0.00080840876363745453];

BERth_su = [0.254921913794786,0.233191431476236,0.211096653092735,0.189008368341342,...
    0.167325847773791,0.146446609406726,0.126733461932933,0.108484732049584,0.0919131757263162,...
    0.0771369160563911,0.0641826854495230,0.0529988839256388,0.0434744067460626,0.0354590676278381,...
    0.0287823671004334,0.0232687053772038,0.0187483912020588,0.0150646803703528,0.0120775474099916,...
    0.00966503894843851,0.00772300227202255,0.00616383474636451,0.00491473054446961,...
    0.00391574864793115,0.00311790514611077,0.00248140489500542,0.00197406790668275,...
    0.00156996788451691,0.00124827788184229,0.000992306076136817,0.000788699342467891];
BER_avg_PANOMA = 0.5*(BER1_PANOMA + BER2_PANOMA);
BER_avg = 0.5*(BER1 + BER2);
%a = ...
%  [0.05 0.95;
%   0.1 0.9;
%   0.2 0.8;
%   0.3 0.7;
%   0.4 0.6];

%% Plotting
figure
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.5, 4.5], 'PaperUnits', 'Inches', 'PaperSize', [7.5, 4.5]);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gcf,'color','w');
ax = gca;
h1 = subplot(1,3,1);
set(h1, 'Position', [0.05 0.1425 0.3 0.8375]);
semilogy(SNR_dB,BER1th,'-b','LineWidth',0.75); hold on
semilogy(SNR_dB,BER1th_PANOMA,'--r','LineWidth',0.75); hold on
semilogy(SNR_dB,BERth_su,'-.k','LineWidth',0.75); hold on
semilogy(SNR_dB(1:5:end),BER1(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER1_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok'); hold on
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.45,'$\mathrm{First \; User}\;(U_{1})$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$'},'Interpreter','latex');
label_h = ylabel('BER','interpreter','latex');
label_h.Position(1) = -7; % change horizontal position of ylabel
label_h.Position(2) = 3e-2; % change vertical position of ylabel
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;

h2 = subplot(1,3,2);
set(h2, 'Position', [0.37 0.1425 0.3 0.8375]);
semilogy(SNR_dB,BER2th,'-b','LineWidth',0.75); hold on
semilogy(SNR_dB,BER2th_PANOMA,'--r','LineWidth',0.75); hold on
semilogy(SNR_dB,BERth_su,'-.k','LineWidth',0.75); hold on
semilogy(SNR_dB(1:5:end),BER2(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER2_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok'); hold on
set(gca,'YTickLabel',[]);
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.45,'$\mathrm{Second \; User}\;(U_{2})$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$'},'Interpreter','latex');
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;

h3 = subplot(1,3,3);
set(h3, 'Position', [0.69 0.1425 0.3 0.8375]);

semilogy(-5,-5,'-b','LineWidth',0.75); hold on
semilogy(-5,-5,'--r','LineWidth',0.75); hold on
semilogy(-5,-5,'-.k','LineWidth',0.75); hold on
semilogy(-5,-5,'ok'); hold on

semilogy(SNR_dB,BERth_avg,'-b','LineWidth',0.75); hold on
semilogy(SNR_dB,BERth_avg_PANOMA,'--r','LineWidth',0.75); hold on
semilogy(SNR_dB,BERth_su(1,:),'-.k','LineWidth',0.75); hold on
semilogy(SNR_dB(1:5:end),BER_avg(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER_avg_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok'); hold on
set(gca,'YTickLabel',[]);
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.45,'$\mathrm{Average}$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$'},'Interpreter','latex');
legend('C-NOMA','PANOMA','SU','Simulations','Location','SouthWest','interpreter','latex','FontSize',9)
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;
% print(gcf,'2UE_ber_new.eps','-depsc','-r600');
