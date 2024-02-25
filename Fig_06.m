clc
clear
close all
SNR_dB = -5:1:35;
a = ...
  [0.05 0.15 0.8;
   0.1 0.2 0.7];

PL = [1 3.9810717055349722 15.848931924611133];

BERth_su = [0.254921913794786,0.233191431476236,0.211096653092735,...
    0.189008368341342,0.167325847773791,0.146446609406726,0.126733461932933,...
    0.108484732049584,0.0919131757263162,0.0771369160563911,0.0641826854495230,...
    0.0529988839256388,0.0434744067460626,0.0354590676278381,0.0287823671004334,...
    0.0232687053772038,0.0187483912020588,0.0150646803703528,0.0120775474099916,...
    0.00966503894843851,0.00772300227202255,0.00616383474636451,0.00491473054446961,...
    0.00391574864793115,0.00311790514611077,0.00248140489500542,0.00197406790668275,...
    0.00156996788451691,0.00124827788184229,0.000992306076136817,0.000788699342467891,...
    0.000626791033821517,0.000498070366772352,0.000395752940590410,...
    0.000314434496891669,0.000249812656113402,0.000198463832435103,...
    0.000157664730252138,0.000125249730198340,9.94970857463029e-05,...
    7.90381964438924e-05];

BER1th = ...
  [0.45809304328427486 0.45496654098286832 0.45076131227762167 0.44504003627451788 ...
   0.43732595968179655 0.42715804768822008 0.41415527485120723 0.39807791975856577 ...
   0.3788744083120919 0.35670559163373122 0.33194289182281406 0.30514117672690833 ...
   0.27699081372026718 0.24825592867700025 0.21970744865102176 0.19205989330282991 ...
   0.16591988289924831 0.14175195298924026 0.11986394807634837 0.10041080609165572 ...
   0.08341278298025416 0.068782663079351075 0.056356375489834951 0.045922404208606271 ...
   0.037246928096762028 0.030093234672837046 0.024235252514839933 0.019465885506906491 ...
   0.015601227598153328 0.012481803106723102 0.0099718496583948713 0.0079574499622172556 ...
   0.0063440979262489455 0.0050540923730609533 0.0040240016028987857 0.0032023339792784239 ...
   0.0025474769026565536 0.0020259205928484514 0.0016107562346922288 0.0012804239722890298 ...
   0.0010176803051444444;
   0.444793438986228 0.44176540780903212 0.43777985831052524 0.43248488667694451 ...
   0.42552060843196771 0.41656375023773473 0.40537460806176545 0.39183582595713395 ...
   0.375973985188835 0.35795931961589789 0.33808453976086417 0.31672891235277034 ...
   0.294316748759022 0.27127941111074366 0.24802711437630032 0.22493253977622812 ...
   0.20232438889282414 0.18048688388080208 0.15966122218753048 0.140046480039376 ...
   0.12179930566298863 0.10503302389312079 0.08981717088165446 0.07617824481943973 ...
   0.064102050384260811 0.053537714645815845 0.044403268604397722 0.036592508488624026 ...
   0.029982632684547081 0.0244419704612706 0.019837082389893984 0.016038652506488851 ...
   0.012925850510066297 0.010389119441893413 0.00833156002864488 0.0066692034035039094 ...
   0.0053304963366275632 0.0042552951435635006 0.0033936063114042292 0.0027042465885470657 ...
   0.0021535363466256008];

BER1th_PANOMA = ...
  [0.45813980030184853 0.45428539833950687 0.4493296996249041 0.44288458674121622 ...
   0.43452779916360851 0.42385298991832332 0.41052789441398729 0.3943491285974825 ...
   0.37528302360626642 0.353485427312271 0.32929805862117434 0.30322331016620319 ...
   0.27588242681217284 0.24796350092561142 0.22016600192822433 0.19314807094870112 ...
   0.16748191644230598 0.1436214330181963 0.12188454471224792 0.10245073062015991 ...
   0.085371993683488176 0.070593703525333173 0.057980785464477963 0.047344834991546556 ...
   0.0384687067300912 0.031126507364193395 0.02509825601394218 0.020179472231366583 ...
   0.016186524178054593 0.01295877633101071 0.01035853755969722 0.00826964536662153 ...
   0.0065953178088347664 0.0052557123879369561 0.0041854739144349118 0.0033314360746666827 ...
   0.0026505601008254137 0.0021081411385279469 0.0016762809541563672 0.0013326079979273098 ...
   0.0010592175308903745;
   0.45878377627762945 0.45486984367684469 0.44984158809574687 0.44331408369447234 ...
   0.43487073540438503 0.42411336315008474 0.41071987894096168 0.39449805767160151 ...
   0.37542490365520548 0.3536647149876076 0.32956365276193889 0.30362299113829955 ...
   0.27645633968650013 0.24873769664415041 0.22114737935499718 0.19432207349252356 ...
   0.16881385235034166 0.14506135788528737 0.12337460374095813 0.10393317157534054 ...
   0.086796035750170961 0.07192003806621601 0.059183353860853175 0.048410282692801176 ...
   0.039394318684336488 0.031917475395169467 0.0257649415249942 0.02073506213472151 ...
   0.016645242022812909 0.013334649774779078 0.010664638286735265 0.008517688476010879 ...
   0.0067955116152395995 0.0054167685086432271 0.0043147101929102427 0.0034349255557990926 ...
   0.0027332957520782158 0.0021741980774551972 0.0017289662269340117 0.0013745933305535085 ...
   0.0010926537704008255];

BER2th = ...
  [0.43957569264164276 0.43500879913974594 0.42961866591221876 0.42311776190912476 ...
   0.41522394611366664 0.40569681518607137 0.39437120207819731 0.38118090708667424 ...
   0.3661678267811953 0.34947510199673881 0.3313267409103694 0.31199923582719236 ...
   0.29179204368991518 0.27100307820555797 0.24991295691350307 0.22877864569098505 ...
   0.207834473148695 0.18729701357370224 0.16737020724289581 0.148247941742158 ...
   0.13011263135606382 0.1131297059318245 0.097439145756868684 0.083146118088700982 ...
   0.070313171348362513 0.058956136832285358 0.049044910671384029 0.040508984861615582 ...
   0.033246463873261009 0.027134718571354202 0.022040871204024648 0.017830789292060598 ...
   0.014375907000410318 0.011557757818962733 0.0092704790242864876 0.0074217280909345124 ...
   0.0059324824442135715 0.0047361385787471066 0.00377723554774477 0.0030100337777092673 ...
   0.00239709948666747;
   0.42952568822328024 0.42503480190124926 0.42015103845111695 0.41471037506508 ...
   0.408549679313135 0.40152932000686148 0.39355291134351267 0.38457958309825813 ...
   0.37462581528666589 0.36375646864471806 0.35206747257322168 0.33966478886702 ...
   0.32664498569387163 0.31308182331457096 0.29902116756185915 0.2844842171171057 ...
   0.26947725500378406 0.25400524631018545 0.23808647307973665 0.22176568612109304 ...
   0.20512368979624618 0.18828178359894279 0.17140014588069485 0.15467015132028084 ...
   0.13830170127287611 0.12250763625540649 0.10748781514104717 0.093415256400538882 ...
   0.080425950405783625 0.06861296356542787 0.058024676581779053 0.048666597369114717 ...
   0.040506054792991283 0.033479012249910761 0.027498145833267873 0.022461261358467668 ...
   0.018259173870375886 0.014782375493353844 0.011926118116633022 0.0095938420118401591 ...
   0.00769911404497893];

BER2th_PANOMA = ...
  [0.43635050871644465 0.43140524494067312 0.42560499835230081 0.41865293705869416 ...
   0.41025216734251968 0.40014273276030504 0.38813734099969188 0.37414885300191575 ...
   0.35820429513808938 0.34044327053191742 0.32110220943092704 0.30048884849875857 ...
   0.27895283824066525 0.25685814520523836 0.23456125119611193 0.21239679616637741 ...
   0.19067008503977093 0.16965436504204531 0.14959020541748158 0.13068461787239533 ...
   0.11310852989574391 0.096992547121920514 0.082422190470232817 0.069434508107130047 ...
   0.058017847784459189 0.048115703620688496 0.03963435083378361 0.032452988800748428 ...
   0.026434674916072168 0.02143647603411486 0.017317777320166636 0.013946293654502617 ...
   0.011201825009579156 0.0089781020683061125 0.0071831918323966987 0.0057389288937515076 ...
   0.0045797654576021224 0.0036513372652730858 0.00290895039374138 0.00231611797152792 ...
   0.0018432192430313105;
   0.43501730511018566 0.42997999898846057 0.424085378180557 0.41703306970017168 ...
   0.40852078726121877 0.39828160833831538 0.38612037914554326 0.37194224325212744 ...
   0.35576795355410418 0.33773369474041337 0.31807666645275823 0.2971105965968967 ...
   0.27519686236605861 0.25271671065943047 0.23004849455981091 0.20755159086124647 ...
   0.18555653216784279 0.16435946602852691 0.14421859523284519 0.12535069977066138 ...
   0.10792691238985144 0.092068187147606242 0.07784186825347196 0.065261034288934608 ...
   0.054287758225115779 0.04484036469310515 0.036803690155446822 0.030040702826509044 ...
   0.024403781045554351 0.019744356047897774 0.01592023188541633 0.012800460098439118 ...
   0.010268037781872852 0.0082208918691353133 0.0065716515716029433 0.005246656308811562 ...
   0.004184551376528442 0.0033347233349186045 0.0026557400414786853 0.0021138929045810828 ...
   0.00168189105373788];

BER3th = ...
  [0.39482428603084108 0.38311023014706846 0.37036496517716633 0.35658651745234848 ...
   0.34180103040902882 0.32606753936009292 0.30948096010769455 0.29217238014577429 ...
   0.27430601016015227 0.25607265399097234 0.23768017255645663 0.219341975784616 ...
   0.20126492494599429 0.18363810926057686 0.16662381460461517 0.15035170170416287 ...
   0.13491679754060126 0.12038137481779239 0.10678015832677791 0.0941276437409043 ...
   0.08242583463499592 0.07167061761554902 0.061855424796760358 0.052971698525851513 ...
   0.045006691740056615 0.037939929953849594 0.031739931408272776 0.026362467597291245 ...
   0.021750944582207318 0.0178387320584395 0.014552750138560902 0.0118174586880756 ...
   0.0095585224363545218 0.0077057012136448749 0.0061947997707648972 0.0049687241583756536 ...
   0.0039778095809230724 0.0031796233281337904 0.0025384343667896048 0.0020245048488762157 ...
   0.0016133164335538075;
   0.40213022853265878 0.39135461646030933 0.37967056901838192 0.36708887427384179 ...
   0.3536471495075556 0.33941334553573832 0.3244874617432541 0.30900091706154548 ...
   0.29311334120939292 0.27700692771239188 0.26087877957702 0.24493178900792098 ...
   0.22936456728336002 0.21436093130863018 0.20007958876936655 0.1866449243715777 ...
   0.1741399992804204 0.16260282348933655 0.1520265513499312 0.14236357676970091 ...
   0.13353280131215423 0.12542886563723171 0.11793201249903072 0.11091746676487561 ...
   0.1042636324423494 0.097858843342454177 0.0916067294708312 0.085430413395948007 ...
   0.079275729731870456 0.0731135068261513 0.06694072830122913 0.060780189219672778 ...
   0.054678178399064206 0.04869984877441258 0.042922321840193609 0.03742614221468879 ...
   0.03228626145126319 0.027564032468948008 0.023301557132267378 0.019519164768397143 ...
   0.016216028223972348];

BER3th_PANOMA = ...
  [0.39520566442068511 0.38355066714122288 0.37087230220548961 0.35716754510352489 ...
   0.3424597966092604 0.32680290226045666 0.31028339105611 0.29302023619207629 ...
   0.27516178000983482 0.2568799691556547 0.23836256427048091 0.21980434576115349 ...
   0.20139842457415344 0.18332859488806877 0.16576333961506762 0.14885174747162466 ...
   0.13272130774131879 0.11747733249251348 0.10320358320953016 0.089963527577341079 ...
   0.077801572003274713 0.0667437162761903 0.056797423438617281 0.047950998421750116 ...
   0.040173182894315027 0.033413764500970361 0.0276057060485191 0.022668783409145549 ...
   0.018514241833988795 0.015049735801171207 0.012183849224851925 0.0098297132060062437 ...
   0.0079075164716420066 0.0063459378181555115 0.00508267542817932 0.0040643065318751309 ...
   0.003245707347828114 0.00258922611560182 0.0020637535601238516 0.0016437887904737769 ...
   0.0013085606418312551;
   0.39551844106817668 0.38390555912031166 0.37127518619007222 0.35762483820786195 ...
   0.34297825562937895 0.32738926298219245 0.31094381156195211 0.29375953041368086 ...
   0.27598246798857845 0.25778122657893088 0.23933921757988241 0.22084611547952143 ...
   0.2024896472912851 0.18444862923161509 0.16688776472752681 0.14995429604566179 ...
   0.13377627636914466 0.11846204915146571 0.10410045854330371 0.090761308685467626 ...
   0.078495625790892254 0.067335413140888242 0.0572928822943769 0.048359537555532184 ...
   0.040505796485287626 0.033681848730044331 0.027820140624402037 0.022839373351922987 ...
   0.018649460218246081 0.015156682125487486 0.01226834013144093 0.00989643922278885 ...
   0.0079602194271718024 0.0063875821297857821 0.0051156000850586664 0.0040903532788243713 ...
   0.0032663254096871319 0.0026055560624400342 0.0020766936570538325 0.0016540471305860616 ...
   0.0013166959713661192];

BERth_avg = ...
  [0.43083100731891955 0.42436185675656085 0.4169149811223356 0.40824810521199695 ...
   0.39811697873483065 0.38630746741146149 0.372669145679033 0.35714373566367141 ...
   0.33978274841781314 0.32075111587381411 0.30031660176321334 0.27882746277957221 ...
   0.25668259411872552 0.23429903871437835 0.21208140672304668 0.19039674689932595 ...
   0.1695570511961815 0.14981011379357828 0.131338104548674 0.114262130524906 ...
   0.0986504163237713 0.0845276622089082 0.071883648681154669 0.060680073607719584 ...
   0.050855597061727049 0.042329767152990666 0.035006698198165577 0.028779112655271104 ...
   0.023532878684540549 0.0191517512455056 0.015521823666993473 0.01253523264745115 ...
   0.010092842454337928 0.0081058504685561864 0.0064964267993167235 0.0051975954095295294 ...
   0.0041525896425977319 0.0033138941665764492 0.0026421420497422009 0.0021049875329581707 ...
   0.0016760320751219072;
   0.4254831185807223 0.41938494205686361 0.41253382192667465 0.40476137867195539 ...
   0.39590581241755274 0.38583547192677814 0.374471660382844 0.36180544203897919 ...
   0.34790438056163131 0.33290757199100263 0.31701026397036858 0.30044183007590375 ...
   0.28344210057875119 0.26624072191131493 0.24904262356917531 0.23202056042163721 ...
   0.21531388105900953 0.19903165122677469 0.18325808220573281 0.16805858097672333 ...
   0.153485265590463 0.13958122437643175 0.12638310975379335 0.11392195430153207 ...
   0.10222246136649545 0.0913013980812255 0.0811659377387587 0.071812726095036972 ...
   0.063228104274067054 0.055389480284283255 0.04826749575763406 0.041828479698425444 ...
   0.036036694567373927 0.030855993488738916 0.026250675900702118 0.022185535658886788 ...
   0.018625310552755545 0.01553390103528845 0.012873760520101543 0.010605751122928123 ...
   0.0086895595385256263];

BERth_avg_PANOMA = ...
  [0.42989865781299275 0.42308043680713425 0.41526900006089812 0.40623502296781178 ...
   0.39574658770512949 0.38359954164636162 0.3696495421565964 0.35383940593049151 ...
   0.33621636625139689 0.316936222333281 0.29625427744086075 0.27450550147537173 ...
   0.25207789654233048 0.22938341367297282 0.20683019757980128 0.18479887152890107 ...
   0.16362443640779856 0.14358437685091835 0.12489277777975322 0.10769962535663209 ...
   0.09209403186083559 0.07810998897448132 0.065733466457776016 0.05491011384014223 ...
   0.045553245802955139 0.037551991828617418 0.03077943763208163 0.025100414813753517 ...
   0.020378480309371852 0.016481662722098923 0.013286721368238594 0.01068188407571013 ...
   0.0085682197633519756 0.0068599174247995264 0.00548378039167031 0.0043782238334311068 ...
   0.0034920109687518832 0.0027829015064676174 0.0022163283026738661 0.0017641715866430021 ...
   0.00140366580525098;
   0.42977317415199723 0.42291846726187227 0.41506738415545869 0.405990663867502 ...
   0.39545659276499423 0.38326141149019755 0.36926135654948566 0.35339994377913664 ...
   0.33572510839929603 0.3163932121023173 0.29565984559819314 0.27385990107157254 ...
   0.25138094978128123 0.228634345511732 0.20602787954744498 0.18394265346647726 ...
   0.16271555362910972 0.14262762435509335 0.12389788583903567 0.1066817266771565 ...
   0.091072857976971539 0.07710787945157016 0.064772701469567345 0.054010284845755985 ...
   0.04472929113157996 0.036813229606106314 0.03012959076828102 0.024538379437717847 ...
   0.019899494428871112 0.01607856264938811 0.012951070101197507 0.010404862599079615 ...
   0.0083412562747614168 0.0066750808358547742 0.0053339872831906169 0.0042573117144783417 ...
   0.0033947241794312633 0.0027048258249379453 0.00215379997515551 0.001714177788573551 ...
   0.0013637469318349416];


% Simulation
% a = ...
%  [0.05 0.15 0.8;
%   0.1 0.2 0.7;
%   0.05 0.3 0.65;
%   0.05 0.1 0.85];

BER1 = ...
  [0.45811072440355122 0.42725992673014657 0.33195863295773409 0.19223864677270647 ...
   0.08362553587392825 0.030140111594776812 0.010105136039727921 0.003135993728012544 ...
   0.000987748024503951;
   0.44470315746868505 0.41643208901082196 0.33796082407835182 0.22492444077611845 ...
   0.1218516312967374 0.053551939771120459 0.019837054075891848 0.006768283338433323 ...
   0.0021795737658524684;
   0.45878530117939764 0.41741410267179463 0.32570947358105284 0.21486439839620322 ...
   0.13052220770558459 0.076497175130649742 0.038690391369217259 0.015859452656094689 ...
   0.0056212856324287354;
   0.4523186891126218 0.42061500251999495 0.3449681069387861 0.22944647860704279 ...
   0.12171906906186188 0.051399022201955594 0.018610634653730692 0.0061391595966808063 ...
   0.0019841679066641867];

BER1_PANOMA = ...
  [0.45812623999752 0.42394215211569575 0.32936151315197371 0.1933161602426795 ...
   0.085565203869592266 0.031189234496531007 0.010500494624010752 0.0032587434825130349 ...
   0.0010299198151603698;
   0.45876773871452259 0.42393771462457075 0.32952318470363057 0.19432811134377731 ...
   0.086820795108409787 0.031967482940034118 0.010634009981980035 0.0035178054643890712 ...
   0.0011156540186919627;
   0.46158782682434635 0.42376305872388254 0.34629616678266645 0.25081979523540954 ...
   0.1662359956530087 0.095573090103819786 0.0433859288531423 0.016139514595970808 ...
   0.0054506609736780526;
   0.45482077785844427 0.422118749512501 0.33133411858176282 0.19826941596116807 ...
   0.090205397714204577 0.033187589874820247 0.011190180744638512 0.0035945396859206283 ...
   0.0011423727152545696];

BER2 = ...
  [0.43964637070725859 0.40579867277765447 0.33146224332551333 0.22867943326613346 ...
   0.12992255265489469 0.058892054090891821 0.022038565297869406 0.0074103758042483913 ...
   0.0024128076743846513;
   0.42960114079771838 0.401430619013762 0.35194070236859526 0.28453529030441937 ...
   0.20528546442907114 0.12247441130117739 0.057959384081231836 0.022435002004995989 ...
   0.0076408753432493136;
   0.4120330665588669 0.37663118423763153 0.31092514377471248 0.22540379919240161 ...
   0.15405744188511622 0.10723847302305395 0.071758715857568281 0.040542122040755919 ...
   0.017798558152883694;
   0.45095589496321009 0.42035086242327513 0.35216974878550245 0.26151807071385857 ...
   0.17679386516226966 0.10363351148297703 0.048533777932444136 0.018467900564198871 ...
   0.0062709093331813339];

BER2_PANOMA = ...
  [0.4364284083931832 0.40021965268569465 0.32126252934994132 0.21230101289797421 ...
   0.11291966478567043 0.048083716332567333 0.017323934102131795 0.0057243166763666468 ...
   0.0018584806580386839;
   0.43512864536770929 0.3982100004549991 0.31789050484399028 0.20752138183223634 ...
   0.10810220567058866 0.044768972962054072 0.015816905866188266 0.0052560207379585243 ...
   0.0016568873112253776;
   0.42275151387197224 0.38330796775906451 0.29866054330391339 0.17943509425481149 ...
   0.081353087293825413 0.030110721028557944 0.010018183088633822 0.0032566028617942765 ...
   0.0010563885122229755;
   0.44155694501111 0.40772521579956839 0.3342703314593371 0.23427712519574961 ...
   0.13876476934546131 0.066532070060859871 0.025983557407885185 0.00896748206503587 ...
   0.002911244177511645];

BER3 = ...
  [0.39490985080529839 0.32609555093389814 0.23776494634510731 0.1503650586448827 ...
   0.082325600973798055 0.038078783217433562 0.014512533474933051 0.004963021323957352 ...
   0.0015838874572250855;
   0.40226861733776531 0.33929558703382595 0.260885821978356 0.18659654868190265 ...
   0.13370621696256607 0.09791463229573541 0.066955397339205316 0.037388159598680804 ...
   0.016230764413471172;
   0.405938547497905 0.3464090103069794 0.27385192104615791 0.2063645403959192 ...
   0.15534378306243388 0.12236756776486447 0.097983522782954441 0.072612198525602953 ...
   0.0444745048009904;
   0.39150813885872227 0.31980539163921673 0.22636114102771795 0.13440741868516262 ...
   0.066278633067733864 0.026938305498389003 0.0095264496971006057 0.0031436187127625746 ...
   0.0010092948564102871];

BER3_PANOMA = ...
  [0.39529577190845616 0.32683204946090105 0.23843446063107873 0.14885046792406414 ...
   0.077694063361873275 0.033533276683446633 0.012152803819392361 0.0040592887564224875 ...
   0.0012838411823176354;
   0.39566864616270769 0.32728278293443414 0.23934070881858235 0.14988763772472455 ...
   0.078661764551470892 0.033793323038353924 0.012279662940674119 0.0040796480907038185 ...
   0.0013001380247239505;
   0.39893417088165822 0.33370744196011609 0.25036842113815772 0.16429931202637596 ...
   0.089408399308201381 0.039181281012437973 0.014313908872182255 0.0047324592850814294 ...
   0.0015542000165999668;
   0.39431588324323352 0.3250526311447377 0.23514529533440934 0.14540794355911288 ...
   0.075852973294053411 0.032622887879224238 0.011945616733766533 0.0039945857608284785 ...
   0.0012930130389739222];
BER_su = ...
  [0.25494044324411352 0.14650948823102353 0.064418121163757669 0.02335457829084342 ...
   0.0077369376511246977 0.0024603857042285915 0.00079660778178443641 0.000236249527500945 ...
   7.9406091187817624E-5;
   0.25485775590948817 0.14611827026345947 0.064229871540256916 0.023190844243311512 ...
   0.0076567346865306267 0.002483713782572435 0.00081896711206577591 0.00025215574568850861 ...
   8.4890455219089564E-5;
   0.25490992768014464 0.14657726934546131 0.064154215441569115 0.023273078453843091 ...
   0.0076522190705618585 0.0025232762034475931 0.00079084216831566332 0.0002517963714072572 ...
   8.0359214281571431E-5;
   0.25486714651570697 0.14651936321127357 0.064238262148475708 0.023302906519186962 ...
   0.0077732032035935927 0.0024421044907910186 0.00077999844000311994 0.0002473745052509895 ...
   7.9531090937818129E-5];

BER_avg_PANOMA = 1/3*(BER1_PANOMA + BER2_PANOMA + BER3_PANOMA);
BER_avg = 1/3*(BER1 + BER2 + BER3);

%% Plotting
figure
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.5, 4.5], 'PaperUnits', 'Inches', 'PaperSize', [7.5, 4.5]);
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(gcf,'color','w');
ax = gca;
h1 = subplot(1,4,1);
set(h1, 'Position', [0.05 0.1425 0.2175 0.8375]);
semilogy(SNR_dB(1:5:end),BER1(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER1_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER1(2,:),'om'); hold on
semilogy(SNR_dB(1:5:end),BER1_PANOMA(2,:),'oc'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok')
semilogy(SNR_dB,BER1th(1,:),'-b'); hold on
semilogy(SNR_dB,BER1th_PANOMA(1,:),':r'); hold on
semilogy(SNR_dB,BER1th(2,:),'-.m'); hold on
semilogy(SNR_dB,BER1th_PANOMA(2,:),'--c'); hold on
semilogy(SNR_dB,BERth_su(1,:),'-.k','LineWidth',1.25)
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.5,'$\mathrm{First \; User}\;(U_{1})$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$','(a)'},'Interpreter','latex');
ylabel('$\mathrm{BER}$','interpreter','latex');
label_h = ylabel('$\mathrm{BER}$','interpreter','latex');
label_h.Position(1) = -7; % change horizontal position of ylabel
label_h.Position(2) = 3e-2; % change vertical position of ylabel
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;

h2 = subplot(1,4,2);
set(h2, 'Position', [0.29 0.1425 0.2175 0.8375]);
semilogy(SNR_dB(1:5:end),BER2(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER2_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER2(2,:),'om'); hold on
semilogy(SNR_dB(1:5:end),BER2_PANOMA(2,:),'oc'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok')
semilogy(SNR_dB,BER2th(1,:),'-b'); hold on
semilogy(SNR_dB,BER2th_PANOMA(1,:),':r'); hold on
semilogy(SNR_dB,BER2th(2,:),'-.m'); hold on
semilogy(SNR_dB,BER2th_PANOMA(2,:),'--c'); hold on
semilogy(SNR_dB,BERth_su(1,:),'-.k','LineWidth',1.25)
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.5,'$\mathrm{Second \; User}\;(U_{2})$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$','(b)'},'Interpreter','latex');
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
set(gca,'YTickLabel',[]);
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;

h3 = subplot(1,4,3);
set(h3, 'Position', [0.53 0.1425 0.2175 0.8375]);
semilogy(SNR_dB(1:5:end),BER3(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER3_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER3(2,:),'om'); hold on
semilogy(SNR_dB(1:5:end),BER3_PANOMA(2,:),'oc'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok')
semilogy(SNR_dB,BER3th(1,:),'-b'); hold on
semilogy(SNR_dB,BER3th_PANOMA(1,:),':r'); hold on
semilogy(SNR_dB,BER3th(2,:),'-.m'); hold on
semilogy(SNR_dB,BER3th_PANOMA(2,:),'--c'); hold on
semilogy(SNR_dB,BERth_su(1,:),'-.k','LineWidth',1.25)
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.5,'$\mathrm{Third \; User}\;(U_{3})$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$','(c)'},'Interpreter','latex');
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
set(gca,'YTickLabel',[]);
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;

h4 = subplot(1,4,4);
set(h4, 'Position', [0.77 0.1425 0.2175 0.8375]);
semilogy(-5,-5,'-b');hold on
semilogy(-5,-5,':r');hold on
semilogy(-5,-5,'-.m');hold on
semilogy(-5,-5,'--c');hold on
semilogy(-5,-5,'-.k','LineWidth',1.25);hold on
semilogy(-5,-5,'ok');hold on

semilogy(SNR_dB(1:5:end),BER_avg(1,:),'ob'); hold on
semilogy(SNR_dB(1:5:end),BER_avg_PANOMA(1,:),'or'); hold on
semilogy(SNR_dB(1:5:end),BER_avg(2,:),'om'); hold on
semilogy(SNR_dB(1:5:end),BER_avg_PANOMA(2,:),'oc'); hold on
semilogy(SNR_dB(1:5:end),BER_su(1,:),'ok')
semilogy(SNR_dB,BERth_avg(1,:),'-b'); hold on
semilogy(SNR_dB,BERth_avg_PANOMA(1,:),':r'); hold on
semilogy(SNR_dB,BERth_avg(2,:),'-.m'); hold on
semilogy(SNR_dB,BERth_avg_PANOMA(2,:),'--c'); hold on
semilogy(SNR_dB,BERth_su(1,:),'-.k','LineWidth',1.25)
axp = get(gca,'Position');
text(axp(1)-3,axp(3)+0.5,'$\mathrm{Average}$','FontSize',11,'FontName', 'Times New Roman','FontWeight','bold', ...
    'Color', [0, 0, 0],'interpreter','latex')
hold off
xlabel({'$E_{b}/N_{0}\, (\mathrm{dB})$','(d)'},'Interpreter','latex');
legend('$\mathrm{C-NOMA}:\;\mathcal{P}_1$','$\mathrm{PANOMA}:\;\mathcal{P}_1$',...
    '$\mathrm{C-NOMA}:\;\mathcal{P}_2$','$\mathrm{PANOMA}:\;\mathcal{P}_2$','SU',...
    'Simulations','interpreter','latex','FontSize',9,'location','southwest')
axis([min(SNR_dB) max(SNR_dB) 10^-3 1])
set(gca,'YTickLabel',[]);
grid on
ax = gca;
ax.Color = 'white';
ax.FontSize = 11;
ax.LineWidth = 0.75;
% print(gcf,'3UE_ber_new.eps','-depsc','-r600');

