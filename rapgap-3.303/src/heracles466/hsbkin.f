C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C...INPUT FOR BADELEK/KWIECINSKI PARAMETRIZATION
C   NEW VERSION VALID FOR EXTENDED KINEMATIC RANGE (DEC. 93)

      SUBROUTINE HSBKIN
      COMMON //H1(2,20,20)
      DIMENSION H(2,20,20)
      DATA (H(1,1,N),N=1,20)/
     .316.3296644740772,   485.8975656910705,  259.6136606319242,
     .117.3307786807664,    32.16818321591835,
     . 10.64768876117214,    2.848065566825837, -0.7215608659036665,
     .  1.187429709239199,  -0.7932390539536663,
     .  0.4714339692858772, -0.1831410742125502, 1.0170687772727414E-03,
     .  8.4179584992745884E-02,-0.1117787295555190,
     .  9.1866157500536292E-02,-6.6935139980532398E-02,
     .  4.6752933908006442E-02,-2.0748727933986192E-02,
     .  1.1678824404468792E-02/
      DATA (H(1,2,N),N=1,20)/
     .254.3089293501125,   418.7477285501654,  230.8904850564271,
     . 93.64447846743268,   32.48257976526893,
     .  6.879896644486573,   2.936391479892251, -4.9063167715756220E-02,
     . -2.4767278893035123E-02, 0.3955014763248987,
     . -0.4533908636720629,  0.3829556835451522,-0.2505685158029276,
     .  0.1182806533021176, -2.6630772372001804E-02,
     . -2.9048052301551008E-02,4.5916476409516690E-02,
     . -4.2200780634649624E-02,3.1875262968723642E-02,
     . -1.6984452043858281E-02/
      DATA (H(1,3,N),N=1,20)/
     . -5.205517203742071,  -7.613923888380566, -2.139391763920144,
     . -0.8778185462976319, -0.6064977299189174,
     .  0.1964926533236207, -0.3638238753020571, 0.2493255700900238,
     . -0.1134190355714483, -1.6085746861335640E-02,
     .  9.3563287149740028E-02,-0.1181647883494726,
     .  0.1050190308667975, -7.1955618163729553E-02,
     .  3.9320613351923985E-02,
     . -1.2009815438625413E-02,-2.7604838214640864E-03,
     .  9.1599608346141909E-03,-9.3558304809775159E-03,
     .  6.7374042415057438E-03/
      DATA (H(1,4,N),N=1,20)/
     .  1.790899045505892,   3.053839222466831,   1.716760795147843,
     .  0.9389845862757266,  0.3552022319675952,
     .  3.0907987227430154E-02, 0.1089572387552265,
     . -8.7460054261984849E-02, 6.0878745078748803E-02,
     . -1.9558925664357729E-02,
     . -1.3279064772305173E-02, 3.0702202030882964E-02,
     . -3.4410989446204764E-02, 2.8108807989698516E-02,
     . -1.8858381016765550E-02,
     .  9.2738410819645021E-03,-2.9301436541826585E-03,
     . -1.0451584670259943E-03, 2.0371195301416922E-03,
     . -2.1740196816710108E-03/
      DATA (H(1,5,N),N=1,20)/
     . -1.125432508063247,  -1.868590239295427, -1.068872737618474,
     . -0.5180069312658861, -0.1529311946681504,
     . -3.1686389286933733E-02,-3.8870789237635146E-02,
     .  3.3116940661247520E-02,-2.6234759141101501E-02,
     .  1.1450631650217661E-02,
     .  1.7056679752367680E-03,-8.8379478159300044E-03,
     .  1.1203033071420430E-02,-9.9215808824230287E-03,
     .  7.0031753409176388E-03,
     . -4.1486685438180383E-03, 1.8415654248478509E-03,
     . -6.4946774876729796E-06,-2.4179162736877685E-04,
     .  6.9752055132026684E-04/
      DATA (H(1,6,N),N=1,20)/
     .  0.7859860842825152,  1.268601328031006,  0.7198376963883954,
     .  0.3337274972420689,  8.3403187833354106E-02,
     .  2.6098164746585593E-02, 2.0937160406365934E-02,
     . -1.6449241435357990E-02, 1.2500976388980837E-02,
     . -4.9348555980043545E-03,
     . -1.2847024081706998E-03, 3.7857560154279305E-03,
     . -4.2505327971354598E-03, 3.5652230238075717E-03,
     . -2.2208002904650896E-03,
     .  1.4900325694985228E-03,-6.9601213564879949E-04,
     . -8.4776803337891749E-05,-7.4556062576086525E-05,
     . -2.5171560143955646E-04/
      DATA (H(1,7,N),N=1,20)/
     . -0.3458275355070987, -0.5560796396633200, -0.3178436561855143,
     . -0.1432259256486573, -3.1191476091947454E-02,
     . -1.3781369721746755E-02,-8.7187118647218875E-03,
     .  7.6580822379931938E-03,-5.8729154797341954E-03,
     .  2.2832713852619971E-03,
     .  5.9895272062973058E-04,-1.6405192637641120E-03,
     .  1.7203622377909946E-03,-1.3468184125611596E-03,
     .  7.4291499573687678E-04,
     . -4.9990725304443526E-04, 2.2481399476398239E-04,
     .  6.2691047846582758E-05, 4.9961656278537477E-05,
     .  9.0458495053711477E-05/
      DATA (H(1,8,N),N=1,20)/
     . -0.1085455474187870, -0.1655087767157225,-8.9218106237787184E-02,
     . -4.3775058549888573E-02,-1.2630632154353540E-02,
     . -1.2351333903456096E-03,-2.5217085547781231E-03,
     . -4.0266867779742427E-04, 1.2853944564771603E-03,
     . -1.2973691451185036E-03,
     .  8.2431922683277483E-04,-1.3689526356241208E-04,
     . -3.2085154615627093E-04, 4.3498770490028007E-04,
     . -4.6939545426076912E-04,
     .  2.4900431886548259E-04,-1.5150459753039519E-04,
     .  1.0209815802393671E-04, 1.4812796521087268E-05,
     . -7.7146939394380492E-06/
      DATA (H(1,9,N),N=1,20)/
     .  0.3438552509916081,  0.5389308757927751,  0.3013583141188885,
     .  0.1389264048545785,  3.2556417603769212E-02,
     .  1.0326741528094614E-02, 8.2785611717296069E-03,
     . -3.8201384005378077E-03, 1.4216865195395041E-03,
     .  7.3202674486942658E-04,
     . -1.6676915387988519E-03, 1.1644228111795703E-03,
     . -4.2808616876797423E-04,-1.5871321459509406E-05,
     .  4.2473000724665851E-04,
     . -2.0456601380010781E-04, 1.7306727773199913E-04,
     . -2.1646329348452131E-04,-4.2131609511711239E-05,
     . -3.1138105795789064E-05/
      DATA (H(1,10,N),N=1,20)/
     . -0.2666284983331815, -0.4193036547186229, -0.2361279192119010,
     . -0.1072196731354234, -2.3509715113688884E-02,
     . -9.1084669244087511E-03,-6.4858604017878814E-03,
     .  3.6723705418607977E-03,-1.7590346279415437E-03,
     . -2.4900690073546540E-04,
     .  1.2769022045740275E-03,-1.0431994895220903E-03,
     .  4.9766983254701150E-04,-1.0532742411831749E-04,
     . -2.6526758757535447E-04,
     .  1.4025948676539830E-04,-1.3392153823798330E-04,
     .  1.7337473250470669E-04, 2.5950665371892346E-05,
     .  2.9269620241400226E-05/
      DATA (H(1,11,N),N=1,20)/
     . -2.1641655953025585E-03,-2.7280006177934300E-03,
     . -7.6672594403688454E-05,-1.2905786933555744E-03,
     . -1.8107522578205785E-03,
     .  7.6144454771083392E-04, 7.9082559269037676E-05,
     . -4.9312363482609612E-04, 3.6864351208980488E-04,
     . -1.0675557950985583E-04,
     . -7.1443297748211056E-05, 1.2751039991787558E-04,
     . -9.5767267168048481E-05, 3.7682934529567924E-05,
     .  6.3722345818114810E-06,
     . -2.7053378158035031E-05, 2.6923218234557537E-05,
     . -1.6266385457451837E-05, 9.6887873322577205E-06,
     . -2.0352265252877557E-06/
      DATA (H(1,12,N),N=1,20)/
     .  0.2318584783547141,     0.3640874804517640,
     .  0.2033123115140910,     9.3531266922188376E-02,
     .  2.2266520869560401E-02,
     .  7.0983081770853293E-03, 5.4250929901558780E-03,
     . -2.6499440057001076E-03, 1.1986084342017536E-03,
     .  2.4374299263204615E-04,
     . -9.6913977487532241E-04, 7.4712172998354589E-04,
     . -3.4173275196392854E-04, 7.7988600742127338E-05,
     .  1.9221963234209065E-04,
     . -7.1059036124745883E-05, 7.3105232343063069E-05,
     . -1.2575779554997444E-04,-3.4568904139992196E-05,
     . -2.3956753201861310E-05/
      DATA (H(1,13,N),N=1,20)/
     . -0.2551126027054218,    -0.4011298971323523,
     . -0.2249831008736070,    -0.1025452462911845,
     . -2.3367491092718078E-02,
     . -8.4353929707860059E-03,-6.0355086331645017E-03,
     .  3.2918676298121958E-03,-1.6157986916136473E-03,
     . -1.6732176165457908E-04,
     .  1.1080362119742959E-03,-9.1425250635234992E-04,
     .  4.5147914913725408E-04,-1.2111266682289053E-04,
     . -2.0931542116845973E-04,
     .  9.2481020975093035E-05,-9.6705156647598560E-05,
     .  1.4932068209479170E-04, 3.1787023323479374E-05,
     .  2.8251344672428605E-05/
      DATA (H(1,14,N),N=1,20)/
     .  7.8075093211252500E-02, 0.1230915801936174,
     .  6.9635656708936525E-02, 3.1196152781560376E-02,
     .  6.4933981650981401E-03,
     .  2.9483403022493872E-03, 1.8859674419684562E-03,
     . -1.2311906604126394E-03, 6.7465998701432196E-04,
     . -1.3952585110695949E-05,
     . -3.5970934183974964E-04, 3.3258202427922218E-04,
     . -1.8329315501507838E-04, 5.9591646797594392E-05,
     .  6.0591765883761565E-05,
     . -3.4972881170466804E-05, 3.7860821889947154E-05,
     . -5.1447163347695427E-05,-6.3823604161686913E-06,
     . -9.7702386617232954E-06/
      DATA (H(1,15,N),N=1,20)/
     .  0.1438319789368712,     0.2258070543013799,
     .  0.1259919888371545,     5.7927903372604848E-02,
     .  1.3837074640170796E-02,
     .  4.4006172533111952E-03, 3.3508418267598121E-03,
     . -1.6264681715208084E-03, 7.2990997529672634E-04,
     .  1.5856831464849014E-04,
     . -6.0311392583147860E-04, 4.6228050942256978E-04,
     . -2.0937636824327253E-04, 4.5772587830595487E-05,
     .  1.2175076251223997E-04,
     . -4.5887893430150516E-05, 4.6800637373385023E-05,
     . -7.8846348562105379E-05,-2.0974026719137824E-05,
     . -1.4884327106003285E-05/
      DATA (H(1,16,N),N=1,20)/
     . -0.2352880318790727,    -0.3698028202050644,
     . -0.2070853502842321,    -9.4516184750886790E-02,
     . -2.1804421877123153E-02,
     . -7.6657171316025243E-03,-5.5279236751719818E-03,
     .  2.9434083313969099E-03,-1.4231371166309967E-03,
     . -1.7553506755579768E-04,
     .  1.0117564418676352E-03,-8.2273032895695508E-04,
     .  3.9978820619553932E-04,-1.0373709599029986E-04,
     . -1.9443139412460696E-04,
     .  8.3282390876259020E-05,-8.6901792439193437E-05,
     .  1.3622111793675436E-04, 3.0072538256514326E-05,
     .  2.5791484490536205E-05/
      DATA (H(1,17,N),N=1,20)/
     .  0.1340924543603120,     0.2109023755965740,
     .  0.1183732053308312,     5.3772887774393171E-02,
     .  1.2123308000737075E-02,
     .  4.5401943808229436E-03, 3.1665582152355960E-03,
     . -1.7804878814693143E-03, 8.9478244297346624E-04,
     .  6.9306597750372190E-05,
     . -5.8572627300951668E-04, 4.9315828092124907E-04,
     . -2.4876897171992880E-04, 6.9666256248087523E-05,
     .  1.0909715545727975E-04,
     . -5.0493693698681596E-05, 5.3342317261339052E-05,
     . -8.0306058696077234E-05,-1.5569325013133603E-05,
     . -1.5231393378166980E-05/
      DATA (H(1,18,N),N=1,20)/
     .  6.9332981920138960E-02, 0.1088471578287316,
     .  6.0712918472132376E-02, 2.7907682015571262E-02,
     .  6.6789757841758525E-03,
     .  2.1203903268691982E-03, 1.6112807904989852E-03,
     . -7.8054412276863090E-04, 3.5061304367886554E-04,
     .  7.5895844327584269E-05,
     . -2.8959938417618394E-04, 2.2195241775818285E-04,
     . -1.0058448101794640E-04, 2.2138700937261226E-05,
     .  5.8425362973301124E-05,
     . -2.1867229809486621E-05, 2.2394644758000066E-05,
     . -3.7938907355849618E-05,-1.0128219344834088E-05,
     . -7.1864603050829368E-06/
      DATA (H(1,19,N),N=1,20)/
     . -0.2094335438293607,    -0.3291457167597506,
     . -0.1842463995640689,    -8.4099819707091857E-02,
     . -1.9455355646261942E-02,
     . -6.8077143644552556E-03,-4.9091727524013329E-03,
     .  2.6028771767728045E-03,-1.2566865356079515E-03,
     . -1.5769782415207982E-04,
     .  8.9727152319202583E-04,-7.2825324762207682E-04,
     .  3.5325983584691317E-04,-9.1551284431029190E-05,
     . -1.7270491484759411E-04,
     .  7.3384593677747381E-05,-7.6736948361805104E-05,
     .  1.2092337784544062E-04, 2.6905347839721798E-05,
     .  2.2944945547767758E-05/
      DATA (H(1,20,N),N=1,20)/
     .  0.1764241083215092,     0.2773161943401925,
     .  0.1553194347090352,     7.0814195640269395E-02,
     .  1.6292501141539833E-02,
     .  5.7896162067337869E-03, 4.1404075170619767E-03,
     . -2.2254584723765772E-03, 1.0853581943524741E-03,
     .  1.2300838389880883E-04,
     . -7.5875689530633292E-04, 6.2123569222668283E-04,
     . -3.0428300054637019E-04, 8.0505847821641867E-05,
     .  1.4493110464378429E-04,
     . -6.2785749290831476E-05, 6.5866331486169732E-05,
     . -1.0272471275959145E-04,-2.2158585089963858E-05,
     . -1.9501578701913792E-05/
      DATA (H(2,1,N),N=1,20)/
     . 315.3555650808810,      487.1523194474932,   259.6465107900499,
     . 116.4940461036096,      33.21540403363760,
     .  9.831397315937661,      3.218098014278122,
     . -0.6817154010258749,     0.9200946379397350,
     . -0.4866104630931506,
     .  0.2383907776459113,    -5.5988745449716497E-02,
     . -3.9178743903936125E-02, 7.2407546263734249E-02,
     . -7.9561173025909792E-02,
     .  5.8214987013573908E-02,-4.0877887493026279E-02,
     .  2.9363377988324933E-02,-1.1066386701897477E-02,
     .  6.9549005586796583E-03/
      DATA (H(2,2,N),N=1,20)/
     .254.3818289740736,      418.5275186398955,   231.1498494922002,
     . 93.49937336431225,      32.41140777734330,
     .  7.152067200268769,      2.587646605884788,
     .  0.2299337074089342,    -0.1499072195943851,
     .  0.3716051170379300,
     . -0.3395322792371598,     0.2481939288576829,
     . -0.1407742546302435,     5.1804775191776264E-02,
     .  9.2148510769340189E-04,
     . -2.9919145279846772E-02, 3.4634000756206363E-02,
     . -2.7546874461335447E-02, 2.0486183161881020E-02,
     . -9.8246558595267586E-03/
      DATA (H(2,3,N),N=1,20)/
     . -5.222753288349486,     -7.569075206092052,
     . -2.195918358368663,     -0.8299124712324074,
     . -0.6210124599979283,
     .  0.1616486397900743,    -0.2883521656087849,
     .  0.1622022683048591,    -4.7157807034083544E-02,
     . -4.4675383118636211E-02,
     .  8.5796938584630368E-02,-8.8589980315254132E-02,
     .  6.9274423787098518E-02,-4.2103020391994220E-02,
     .  1.9617402469487285E-02,
     . -2.7872258929548117E-03,-5.1104801020305181E-03,
     .  7.2200152134024414E-03,-6.7037539882323981E-03,
     .  4.1539590637957897E-03/
      DATA (H(2,4,N),N=1,20)/
     .  1.795655239403132,      3.042305197786875,
     .  1.731510650345955,      0.9245925080343382,
     .  0.3635697359701548,
     .  3.4770334331769240E-02, 9.1904811816968064E-02,
     . -6.2496362250700219E-02, 3.7414104058495009E-02,
     . -4.6536547847940985E-03,
     . -1.6938004389858978E-02, 2.5626445644071382E-02,
     . -2.4666614104332649E-02, 1.8140589733304925E-02,
     . -1.0864519578660768E-02,
     .  4.4550266994742695E-03,-5.5627123898424732E-04,
     . -1.4020650074635405E-03, 1.7417533968980771E-03,
     . -1.4283687571557594E-03/
      DATA (H(2,5,N),N=1,20)/
     . -1.126853083480120,     -1.865303194677924,
     . -1.073102970239374,     -0.5135741672826667,
     . -0.1562121942761721,
     . -3.1632271398583714E-02,-3.4828615756586772E-02,
     .  2.5849725963632295E-02,-1.8438367321697530E-02,
     .  5.5489037460329903E-03,
     .  4.2475786663977816E-03,-8.3539999305615269E-03,
     .  8.6737311523375181E-03,-6.8302270862958138E-03,
     .  4.1674577161149381E-03,
     . -2.2165421159097054E-03, 6.7951339675775276E-04,
     .  3.7790330032298208E-04,-3.2587438691072173E-04,
     .  5.0242900463093502E-04/
      DATA (H(2,6,N),N=1,20)/
     .  0.7864312243380061,      1.267607980471678,
     .  0.7211224905345020,      0.3323237488854016,
     .  8.4593623825662339E-02,
     .  2.5804678659169554E-02, 1.9951615107583772E-02,
     . -1.4277573204395454E-02, 9.9344292404631553E-03,
     . -2.7679879559438641E-03,
     . -2.4419729832566975E-03, 3.9342797751438191E-03,
     . -3.6125020449757010E-03, 2.6267220260744386E-03,
     . -1.2588652671254420E-03,
     .  7.8114114261018274E-04,-2.1923412548235279E-04,
     . -2.8171486406146958E-04, 5.5576500129959715E-06,
     . -2.0435510342071948E-04/
      DATA (H(2,7,N),N=1,20)/
     . -0.3459714042534403,    -0.5557681331266985,
     . -0.3182485416066988,    -0.1427718698835976,
     . -3.1613247839279500E-02,
     . -1.3618499107698309E-02,-8.4773043724060808E-03,
     .  6.9945056072709231E-03,-5.0253027746090357E-03,
     .  1.5091899675752218E-03,
     .  1.0657727843081655E-03,-1.7743351883985245E-03,
     .  1.5658858287012425E-03,-1.0636890094487570E-03,
     .  4.2176733617930925E-04,
     . -2.4920864200303785E-04, 4.2659968471616453E-05,
     .  1.4690131523422416E-04, 8.2876816011666943E-06,
     .  8.0107476865693879E-05/
      DATA (H(2,8,N),N=1,20)/
     . -0.1084980810772992,    -0.1656088807985460,
     . -8.9087376386210287E-02,-4.3924074265384052E-02,
     . -1.2482577849853671E-02,
     . -1.3065344936064699E-03,-2.5793676278996006E-03,
     . -1.9668365825810749E-04, 1.0041118427874141E-03,
     . -1.0240666464075246E-03,
     .  6.4550834159022521E-04,-6.7816159392562099E-05,
     . -2.8648780311074406E-04, 3.4992809259848387E-04,
     . -3.6282081763720117E-04,
     .  1.6185832607705037E-04,-8.4279742396267165E-05,
     .  6.8736558185330570E-05, 3.3157017208563402E-05,
     . -5.9345829973611333E-06/
      DATA (H(2,9,N),N=1,20)/
     .  0.3438393579778490,     0.5389635985445374,
     .  0.3013154171255550,     0.1389757547543309,
     .  3.2504722326641548E-02,
     .  1.0355376472440389E-02, 8.2913668159813562E-03,
     . -3.8847017032601163E-03, 1.5153449687520053E-03,
     .  6.3614404340656736E-04,
     . -1.6010333670451369E-03, 1.1339006548746402E-03,
     . -4.3426387082593887E-04, 9.4917241579444531E-06,
     .  3.8947712947007999E-04,
     . -1.7455910817432596E-04, 1.4873062177360617E-04,
     . -2.0375990653644737E-04,-4.9620779049301657E-05,
     . -3.1206664618377173E-05/
      DATA (H(2,10,N),N=1,20)/
     . -0.2666231040709893,    -0.4193145044508198,
     . -0.2361136868986530,    -0.1072360910235140,
     . -2.3491759328072356E-02,
     . -9.1194205360529169E-03,-6.4882354386895188E-03,
     .  3.6927084796933584E-03,-1.7902600447898725E-03,
     . -2.1552806369962549E-04,
     .  1.2524701436853841E-03,-1.0306463554775114E-03,
     .  4.9803908648807727E-04,-1.1278777536861013E-04,
     . -2.5364485080527975E-04,
     .  1.2999604384653926E-04,-1.2521969763617665E-04,
     .  1.6865130501965398E-04, 2.8883403189738144E-05,
     .  2.9123344756608840E-05/
      DATA (H(2,11,N),N=1,20)/
     . -2.1660351807644072E-03,-2.7243364001492458E-03,
     . -8.1435272540315466E-05,-1.2851179429450414E-03,
     . -1.8169381076812112E-03,
     .  7.6549364276813217E-04, 7.9308966343408952E-05,
     . -4.9954504267920127E-04, 3.7904763654213230E-04,
     . -1.1838963598946141E-04,
     . -6.2601422280735362E-05, 1.2255928196267087E-04,
     . -9.5336761775792050E-05, 3.9830327528115821E-05,
     .  2.5565310038369642E-06,
     . -2.3563899255824005E-05, 2.3840774617620322E-05,
     . -1.4538768861509441E-05, 8.5705688071448483E-06,
     . -1.9310428292372785E-06/
      DATA (H(2,12,N),N=1,20)/
     .  0.2318591581225479,     0.3640861941488209,
     .  0.2033139209296905,     9.3529466274897100E-02,
     .  2.2268614225486147E-02,
     .  7.0968688618638989E-03, 5.4251737911641374E-03,
     . -2.6479100630039140E-03, 1.1951494982517930E-03,
     .  2.4776160118547664E-04,
     . -9.7230098779385421E-04, 7.4901618532146762E-04,
     . -3.4205697579785480E-04, 7.7389287135517535E-05,
     .  1.9346583836977588E-04,
     . -7.2237973141837467E-05, 7.4188032169166135E-05,
     . -1.2638124806418650E-04,-3.4150899972068985E-05,
     . -2.4009607069959709E-05/
      DATA (H(2,13,N),N=1,20)/
     . -0.2551128802029598,    -0.4011294003881104,
     . -0.2249836551117676,    -0.1025446709226060,
     . -2.3368165547998303E-02,
     . -8.4349197469735719E-03,-6.0355616009768629E-03,
     .  3.2912148917943074E-03,-1.6146531788159599E-03,
     . -1.6869615025303481E-04,
     .  1.1091480240063034E-03,-9.1495592957547284E-04,
     .  4.5164275442913885E-04,-1.2095138684732984E-04,
     . -2.0972035015917538E-04,
     .  9.2876388617450081E-05,-9.7082065154004003E-05,
     .  1.4954255380035639E-04, 3.1633576812386239E-05,
     .  2.8274623990525821E-05/
      DATA (H(2,14,N),N=1,20)/
     .  7.8075235376370969E-02, 0.1230913449954655,
     .  6.9635856975744840E-02, 3.1195988326370767E-02,
     .  6.4935814396565525E-03,
     .  2.9482200390526297E-03, 1.8859683953848872E-03,
     . -1.2309708030527523E-03, 6.7428314752218884E-04,
     . -1.3492726150342946E-05,
     . -3.6008765732794689E-04, 3.3283120260798764E-04,
     . -1.8336057445542679E-04, 5.9548638500715326E-05,
     .  6.0723096244574487E-05,
     . -3.5104120561733622E-05, 3.7990281705310375E-05,
     . -5.1524502850817283E-05,-6.3273339575380029E-06,
     . -9.7794431641205986E-06/
      DATA (H(2,15,N),N=1,20)/
     .  0.1438318812767813,     0.2258072042788636,
     .  0.1259919072418009,     5.7927929457090424E-02,
     .  1.3837061884805410E-02,
     .  4.4006090730575940E-03, 3.3508729547767256E-03,
     . -1.6265535015216533E-03, 7.3003221412595660E-04,
     .  1.5842369175691944E-04,
     . -6.0299655557464298E-04, 4.6220200068652175E-04,
     . -2.0935573839129444E-04, 4.5786152006797103E-05,
     .  1.2170767562643602E-04,
     . -4.5845125254479122E-05, 4.6757378864284945E-05,
     . -7.8820587763527001E-05,-2.0992870871451951E-05,
     . -1.4881256790120993E-05/
      DATA (H(2,16,N),N=1,20)/
     . -0.2352879476493046,    -0.3698029441522917,
     . -0.2070853081732636,    -9.4516163586214274E-02,
     . -2.1804469324534519E-02,
     . -7.6656614736123753E-03,-5.5279712518447083E-03,
     .  2.9434525948627788E-03,-1.4231753036239052E-03,
     . -1.7549895555224489E-04,
     .  1.0117316711568970E-03,-8.2271541058852500E-04,
     .  3.9978841328645314E-04,-1.0374450655089577E-04,
     . -1.9441647123489869E-04,
     .  8.3269069670810831E-05,-8.6888408161558546E-05,
     .  1.3621365092123308E-04, 3.0078193672955222E-05,
     .  2.5790964380796467E-05/
      DATA (H(2,17,N),N=1,20)/
     .  0.1340923727789033,     0.2109024936138492,
     .  0.1183731759744742,     5.3772849717184081E-02,
     .  1.2123377974822865E-02,
     .  4.5401199766745579E-03, 3.1666142102052796E-03,
     . -1.7805202206368379E-03, 8.9479294892310841E-04,
     .  6.9308027891171892E-05,
     . -5.8573458849584875E-04, 4.9316714692979459E-04,
     . -2.4877815808025492E-04, 6.9673037398030878E-05,
     .  1.0909108718445778E-04,
     . -5.0490121639002391E-05, 5.3339261680800218E-05,
     . -8.0305072546186575E-05,-1.5570187753166633E-05,
     . -1.5231913605119485E-05/
      DATA (H(2,18,N),N=1,20)/
     .  6.9333064756997220E-02, 0.1088470386760217,
     .  6.0712944150708274E-02, 2.7907726960858526E-02,
     .  6.6788959412898632E-03,
     .  2.1204735436808044E-03, 1.6112198446343334E-03,
     . -7.8051463146365807E-04, 3.5061166757438840E-04,
     .  7.5881089234121166E-05,
     . -2.8957883517329916E-04, 2.2193434253996885E-04,
     . -1.0057132561036559E-04, 2.2131486682815538E-05,
     .  5.8428725478745573E-05,
     . -2.1867576938304706E-05, 2.2394106721910335E-05,
     . -3.7937567302082411E-05,-1.0129117078619442E-05,
     . -7.1855069091663770E-06/
      DATA (H(2,19,N),N=1,20)/
     . -0.2094336295158042,    -0.3291455937230821,
     . -0.1842464246997682,    -8.4099868361146289E-02,
     . -1.9455269978078640E-02,
     . -6.8078031852446442E-03,-4.9091081178749376E-03,
     .  2.6028476632793045E-03,-1.2566882207296329E-03,
     . -1.5767791759710169E-04,
     .  8.9724593667888015E-04,-7.2823118838405930E-04,
     .  3.5324469912063817E-04,-9.1543537856009860E-05,
     . -1.7270751195724047E-04,
     .  7.3383856366975696E-05,-7.6735115913733591E-05,
     .  1.2092115497396109E-04, 2.6906925838931395E-05,
     .  2.2943793016116745E-05/
      DATA (H(2,20,N),N=1,20)/
     .  0.1764241977012391,     0.2773160660648814,
     .  0.1553194604372024,     7.0814247132839963E-02,
     .  1.6292410696607291E-02,
     .  5.7897098459946868E-03, 4.1403394690611027E-03,
     . -2.2254279622666902E-03, 1.0853609872971168E-03,
     .  1.2298599076304767E-04,
     . -7.5872864510608314E-04, 6.2121142480548227E-04,
     . -3.0426659804017821E-04, 8.0497599697570645E-05,
     .  1.4493354353069614E-04,
     . -6.2784615612362718E-05, 6.5863971103739356E-05,
     . -1.0272209329061569E-04,-2.2160470633835152E-05,
     . -1.9500311080302588E-05/
      DO 200 I=1,2
      DO 200 IQ=1,20
      DO 200 N=1,20
      H1(I,IQ,N)=H(I,IQ,N)
200   CONTINUE
      END