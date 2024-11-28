* Last file till now
* My Biogas Supply Chain- Hamed Semsarian
* In most lines of code below, the letter "p" at the end of words and symbols, represents character prime (') e.g. b'. In the same way, "pp" represents double prime character (") e.g. b".

Sets
i_org Biomass source points of main suppliers /abadeh, ardakan, arsenjan, bakhtegan, bavanat, beyza, darab, eqlid, estahban, evaz, farashband, fasa, firuzabad, gerash, ghir, jahrom, kavar, kazerun, kerehee, khafr, kharameh, khonj, lamerd, lar, marvdasht, mohr, neyriz, nurabad, pasargad, qaemyeh, rostam, sadra, safashahr, sarvestan, shiraz, zarindasht, zarqan, eqlid_WWT, jahrom_WWT, kazerun_WWT, shiraz_industrial_WWT, shiraz_WWT_No_1, shiraz_WWT_No_2/
i(i_org) Biomass source points of main suppliers /abadeh, ardakan, arsenjan, bakhtegan, bavanat, beyza, darab, eqlid, estahban, evaz, farashband, fasa, firuzabad, gerash, ghir, jahrom, kavar, kazerun, kerehee, khafr, kharameh, khonj, lamerd, lar, marvdasht, mohr, neyriz, nurabad, pasargad, qaemyeh, rostam, sadra, safashahr, sarvestan, shiraz, zarindasht, zarqan, eqlid_WWT, jahrom_WWT, kazerun_WWT, shiraz_industrial_WWT, shiraz_WWT_No_1, shiraz_WWT_No_2/
ip(i) wastewater treatment centers /eqlid_WWT, jahrom_WWT, kazerun_WWT, shiraz_industrial_WWT, shiraz_WWT_No_1, shiraz_WWT_No_2/
ipp(i) cities /abadeh, ardakan, arsenjan, bakhtegan, bavanat, beyza, darab, eqlid, estahban, evaz, farashband, fasa, firuzabad, gerash, ghir, jahrom, kavar, kazerun, kerehee, khafr, kharameh, khonj, lamerd, lar, marvdasht, mohr, neyriz, nurabad, pasargad, qaemyeh, rostam, sadra, safashahr, sarvestan, shiraz, zarindasht, zarqan/
j_org Types of biomasses: Animal_manure. Agricultural_waste. Waste_water /Chicken, Cow, Horse, Agri_waste, Urban, Urban_ready/
j(j_org) Types of biomasses: Animal_manure. Agricultural_waste. Waste_water /Chicken, Cow, Horse, Agri_waste, Urban, Urban_ready/
jnp(j) Types of livestock biomass /Chicken, Cow, Horse, Agri_waste/
jp(j) Types of wastewaters available /Urban, Urban_ready/
j1(jp) Types of wastewaters available /Urban/
j2(jp) /Urban_ready/
b Locations related to biorefineries /point1, point4, point5, point6, point8, point10/
bp(b) Candidate points /point1, point4, point5, point6, point8, point10/
*bpp(b) Existing points //
p Locations related to power plants /hafez_power_plant, jahrom_power_plant, kazerun_power_plant, shiraz_power_plant/
d Candidate points of distribution centers /point1, point4, point5, point6, point8, point10/
u Capacity scales for warehouses /small, medium, large/
n Types of products /Raw_Bg, Biofert, SCP/
np(n) Products except Raw Biogas /Biofert, SCP/
o Market zones /abadeh, ardakan, arsenjan, bakhtegan, bavanat, beyza, darab, eqlid, estahban, evaz, farashband, fasa, firuzabad, gerash, ghir, jahrom, kavar, kazerun, kerehee, khafr, kharameh, khonj, lamerd, lar, marvdasht, mohr, neyriz, nurabad, pasargad, qaemyeh, rostam, sadra, safashahr, sarvestan, shiraz, zarindasht, zarqan/
w Power and Energy demand zones /abadeh, darab, jahrom, kazerun, lamerd, lar, shiraz, hafez_power_plant, jahrom_power_plant, kazerun_power_plant, shiraz_power_plant/
m Market zones of products /abadeh, ardakan, arsenjan, bakhtegan, bavanat, beyza, darab, eqlid, estahban, evaz, farashband, fasa, firuzabad, gerash, ghir, jahrom, kavar, kazerun, kerehee, khafr, kharameh, khonj, lamerd, lar, marvdasht, mohr, neyriz, nurabad, pasargad, qaemyeh, rostam, sadra, safashahr, sarvestan, shiraz, zarindasht, zarqan/
k Production technologies of biogas: Chamber_Floating_Drum_Type_Biogas_Plant. Digester_Floating_Drum_Type_Biogas_Holder. Fixed_Dome_Type_Biogas_Plant /CFDTBP, DFDTBH, FDTBP/
f Capacity levels for production technology k /low, medium, high/
h Purification technologies of biogas: Pressure_Swing_Adsorption. High_Pressure_Water_Scrubbing. Organic_Physical_Solvent. Cryogenic_Separation_Process. Membrane_Separation /PSA, HPWS, OPS, CSP, MS/
g Capacity levels for purification technology h /low, medium, high/
l Modes of transportation /pipe, rail, truck/
lpg(l) /pipe, rail, truck/
lp(lpg) pipeline /pipe/
lg(l) /rail, truck/
lg1(lg) /rail/
lg2(lg) /truck/
v Job classes /Class1, Class2, Class3/
t_org /t0*t10/
t(t_org) Time periods /t1*t10/
*q Kinds of markets /elecrticity_market, biological_products_market/
r Railway stations /station1, station2, station3, station4, station5, station6, station7/
s Types of fossil fuels used in power plants /shiraz_oil_refinery, farashband_gas_refinery/
c Capacity scales for warehouses /small, medium, large/
e Types of energy that can be produced /electricity, thermal/
e1(e) /electricity/
e2(e) /thermal/
counter2 /1*15/
counter3 /1*15/
alias (r, rp)
;

********************************************************************************

parameter
loopCount2
loopCount3
;

Scalar
*epsilon2 /10E-15/
*epsilon3 /10E-15/
EPSILON /10E-6/
epsilon1
epsilon2
epsilon3
step2 /224994/
step3 /505673/
LowerBound2 /60943190/
LowerBound3 /241171500/
q2 /14/
q3 /14/
;

*epsilon1 = LowerBound1 + q * step1;
*epsilon2 = LowerBound2 + q2 * step2;
*epsilon3 = LowerBound3 + q3 * step3;

Scalar
horizon /10/
pi /3.141592653589/
phi Interest rate /0.24/
alpha Annual discount rate /0.24/
maxDist Maximum distance allowed for pipeline /70/
bigM Big M /10E+15/
Diameter pipe diameter /0.5/
Velociy Fluid velocity in pipe /1.5/
lamda /31536000/
rho /1/
Cross_section
Discharge
cap_pipe
;

Cross_section = (pi * POWER(Diameter,2))/4;
Discharge = Velociy * Cross_section;
cap_pipe = lamda * rho * Discharge;

********************************************************************************
* The factor of social benefit
Parameters
chi_cons(bp,k,h,c,v)
chi_pipe(v)
chi_src(i,j,v)
chi_prod(b,k,h,c,v)
chi_PowP(p,v)
chi_trans_i_b(i, b, lg, v)
chi_trans_b_p(b, p, lg, v)
chi_trans_b_m(b, m, lg, v)
chi_trans_s_p(s, p, lg, v)
;

chi_cons(bp,k,h,c,v) = 1;
chi_pipe(v) = 1;
chi_src(i,j,v) = 1;
chi_prod(b,k,h,c,v) = 1;
chi_PowP(p,v) = 1;

chi_trans_i_b(i, b, lg1, "Class3") = 5;
chi_trans_i_b(i, b, lg2, "Class3") = 2;

chi_trans_b_p(b, p, lg1, "Class3") = 5;
chi_trans_b_p(b, p, lg2, "Class3") = 2;

chi_trans_b_m(b, m, lg1, "Class3") = 5;
chi_trans_b_m(b, m, lg2, "Class3") = 2;

chi_trans_s_p(s, p, lg1, "Class3") = 5;
chi_trans_s_p(s, p, lg2, "Class3") = 2;

********************************************************************************

Parameters
*time(lg) Travel time for each km (minutes) rail 50km.h^-1 truck 80km.h^-1 /rail 1.2, truck 0.75/
time(lg) Travel time for each km (hours) rail 50km.h^-1 truck 80km.h^-1 /rail 0.02, truck 0.0125/
sigma(b) Cost factor according to the climate
sigma_pipe(i,b)
;

sigma(b) = 0;
sigma(bp) = 1;
sigma_pipe(i,b) = 1;

********************************************************************************
* Road efficiency factors

Parameter
omega_i_b(i, b)
omega_i_r(i, r)

omega_b_p(b, p)
omega_b_d(b, d)
omega_b_r(b, r)

omega_d_m(d, m)
omega_d_r(d, r)

omega_s_p(s, p)
omega_s_r(s, r)

omega_r_b(r, b)
omega_r_d(r, d)
omega_r_p(r, p)
omega_r_m(r, m)

omega_r_rp(r, rp)
;

omega_i_b(i, b) = 1;
omega_i_r(i, r) = 1;

omega_b_p(b, p) = 1;
omega_b_d(b, d) = 1;
omega_b_r(b, r) = 1;

omega_d_m(d, m) = 1;
omega_d_r(d, r) = 1;

omega_s_p(s, p) = 1;
omega_s_r(s, r) = 1;

omega_r_b(r, b) = 1;
omega_r_d(r, d) = 1;
omega_r_p(r, p) = 1;
omega_r_m(r, m) = 1;

omega_r_rp(r, rp) = 1;


********************************************************************************
* Loading capacity parameters

Parameters
lc_Bm(lg, j)
lc_prdct(lg, n)
lc_PBg(lg) Loading capacity  /rail 800, truck 20/
lc_Fuel(lg, s) Loading capacity
mpc_Bm(lg, j) Minimum percentage of capacity allowed to be loaded
mpc_prdct(lg, n) Minimum percentage of capacity allowed to be loaded
mpc_PBg(lg) Minimum load percentage
mpc_Fuel(lg,s) Minimum load percentage
;

lc_Bm("rail", j) = 1000;
lc_Bm("truck", j) = 30;

$CALL GDXXRW BSC.xlsx par=lc_prdct rng=Sheet1!B7:E9 rdim=1 cdim=1
$gdxin BSC.gdx
$load lc_prdct
$gdxin

$CALL GDXXRW BSC.xlsx par=lc_Fuel rng=Sheet1!K2:M4 rdim=1 cdim=1
$gdxin BSC.gdx
$load lc_Fuel
$gdxin

mpc_Bm("rail", j) = 0.5;
mpc_Bm("truck", j) = 0.8;
mpc_prdct("rail", n) = 0.5;
mpc_prdct("truck", n) = 0.8;
mpc_PBg("rail") = 0.5;
mpc_PBg("truck") = 0.8;
mpc_Fuel("rail",s) = 0.5;
mpc_Fuel("truck",s) = 0.8;


********************************************************************************
* Distances

Parameter
d_i_b(i, b)
*d_i_r(i, r)

d_b_p(b, p)
*d_b_r(b, r)

d_b_m(b, m)

d_s_p(s, p)
*d_s_r(s, r)

*d_r_b(r, b)
*d_r_p(r, p)
*d_r_m(r, m)
*d_r_rp(r, rp)
;


$CALL GDXXRW BSC.xlsx par=d_i_b rng=distance_i_b!A1:G44 rdim=1 cdim=1
$gdxin BSC.gdx
$load d_i_b
$gdxin

*$CALL GDXXRW BSC.xlsx par=d_i_r rng=distance_i_r!A1:H44 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_i_r
*$gdxin

$CALL GDXXRW BSC.xlsx par=d_b_p rng=distance_b_p!A1:E7 rdim=1 cdim=1
$gdxin BSC.gdx
$load d_b_p
$gdxin

*$CALL GDXXRW BSC.xlsx par=d_b_r rng=distance_b_r!A1:H7 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_b_r
*$gdxin

$CALL GDXXRW BSC.xlsx par=d_b_m rng=distance_d_m!A1:AL7 rdim=1 cdim=1
$gdxin BSC.gdx
$load d_b_m
$gdxin


$CALL GDXXRW BSC.xlsx par=d_s_p rng=distance_s_p!A1:E3 rdim=1 cdim=1
$gdxin BSC.gdx
$load d_s_p
$gdxin

*$CALL GDXXRW BSC.xlsx par=d_s_r rng=distance_s_r!A1:H3 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_s_r
*$gdxin

*$CALL GDXXRW BSC.xlsx par=d_r_b rng=distance_r_b!A1:G8 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_r_b
*$gdxin

*$CALL GDXXRW BSC.xlsx par=d_r_p rng=distance_r_p!A1:E8 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_r_p
*$gdxin

*$CALL GDXXRW BSC.xlsx par=d_r_m rng=distance_r_m!A1:AL8 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_r_m
*$gdxin

*$CALL GDXXRW BSC.xlsx par=d_r_rp rng=distance_r_rp!A1:H8 rdim=1 cdim=1
*$gdxin BSC.gdx
*$load d_r_rp
*$gdxin

********************************************************************************
* Binary Parameter

Parameter
wbd(i, b) maxDist radius: source point to biorefinery;

loop((i,b),
if(
d_i_b(i, b) <= maxDist,
wbd(i, b) = 1;
else
wbd(i, b) = 0;
);
);

********************************************************************************
* Binary parameter
Parameter
z_src(i, j);

$CALL GDXXRW BSC.xlsx par=z_src rng=Sheet4!P2:V45 rdim=1 cdim=1
$gdxin BSC.gdx
$load z_src
$gdxin

********************************************************************************
* Binary parameter

Table tcb(k, j) Suitable technology for biomass
         Chicken Cow Horse Agri_waste Urban Urban_ready
CFDTBP   1       1   1     1          0     0
DFDTBH   1       1   1     1          1     1
FDTBP    1       1   0     1          1     1
;

********************************************************************************
* Water requirements

Parameter
c_purch_Bm(j, i, t)
wr(j, k) Water requirement for biomass processing
c_water(t) water cost per unit;

c_purch_Bm(j, i, t) = 2;
c_purch_Bm("Urban", i, t) = 0;
c_purch_Bm("Urban_ready", i, t) = 0;

wr(j, k) = 1000;
c_water(t) = 1;

********************************************************************************
* Transporting costs
Parameters
c_trans_Bm(j, l, t)
c_trans_prdct(n, lg, t)
c_trans_PBg(lg, t)
c_trans_Fuel(s, lg, t)
;

c_trans_Bm(j, "pipe", t) = 0;
c_trans_Bm(j, "rail", t) = 50;
c_trans_Bm(j, "truck", t) = 100;

c_trans_prdct(n, "rail", t) = 100;
c_trans_prdct(n, "truck", t) = 200;

c_trans_PBg("rail", t) = 100;
c_trans_PBg("truck", t) = 200;

c_trans_Fuel(s, "rail", t) = 100;
c_trans_Fuel(s, "truck", t) = 200;

********************************************************************************
* Process and production costs
Parameters
c_pre_Bm(j, b, t) Pre-processing cost of biomass
c_prod(n, k, t)
c_PBg_prod(h,t)
c_elec(p, t)
c_Bg_CHP(c, t)
c_PBg_CHP(p, t)
c_Fuel_CHP(s, p, t)
;

*c_prod(n, k, t) = 1;

*c_prod("Raw_Bg", "CFDTBP", t) = 0.47;
*c_prod("Raw_Bg", "DFDTBH", t) = 0.48;
*c_prod("Raw_Bg", "FDTBP", t) = 0.46;
c_prod("Raw_Bg", "CFDTBP", t) = 405.61;
c_prod("Raw_Bg", "DFDTBH", t) = 414.24;
c_prod("Raw_Bg", "FDTBP", t) = 396.98;

c_prod("Biofert", k, t) = 30;
c_prod("SCP", k, t) = 50;


c_pre_Bm(j, b, t) = 1;
c_pre_Bm("Urban_ready", b, t) = 0;


*c_PBg_prod("PSA",t) = 0.4;
*c_PBg_prod("HPWS",t) = 0.13;
*c_PBg_prod("OPS",t) = 0.11;
*c_PBg_prod("CSP",t) = 0.17;
*c_PBg_prod("MS",t) = 0.12;
c_PBg_prod("PSA",t) = 345.2;
c_PBg_prod("HPWS",t) = 112.19;
c_PBg_prod("OPS",t) = 94.93;
c_PBg_prod("CSP",t) = 146.71;
c_PBg_prod("MS",t) = 103.56;

*c_elec(p, t) = 0.13;
c_elec(p, t) = 0.05;

c_Bg_CHP(c, t) = 0.05;
c_PBg_CHP(p, t) = 0.05;
c_Fuel_CHP(s , p, t) = 0.05;

********************************************************************************
* Maintenance costs

Parameters
c_maint_Bref(b, k, h, c, t)
c_maint_PowP(p, t)
;

*c_maint_Bref(b, k, h, c, t) = 1;
c_maint_Bref(b, k, "PSA", c, t) = 106000;
c_maint_Bref(b, k, "HPWS", c, t) = 65000;
c_maint_Bref(b, k, "OPS", c, t) = 89000;
c_maint_Bref(b, k, "CSP", c, t) = 109000;
c_maint_Bref(b, k, "MS", c, t) = 75000;

c_maint_PowP(p, t) = 150000;

********************************************************************************
* Establishment costs

Scalar
c_pipe Establishment cost /30000/
;

Parameters
c_invst(k, h, c)
;
c_invst(k, h, c) = 1;
c_invst(k, "PSA", c) = 9092075;
c_invst(k, "HPWS", c) = 9499350;
c_invst(k, "OPS", c) = 10414000;
c_invst(k, "CSP", c) = 10659067;
c_invst(k, "MS", c) = 10382000;

*c_invst(k, "PSA", c) = 4092075;
*c_invst(k, "HPWS", c) = 4499350;
*c_invst(k, "OPS", c) = 5414000;
*c_invst(k, "CSP", c) = 5659067;
*c_invst(k, "MS", c) = 5382000;

********************************************************************************
* Demands
Parameters
q_min_prdct(n, m, t)
q_max_prdct(n, m, t)

q_min_PBg(m, t)
q_max_PBg(m, t)

q_min_nrg(e, b, t)
q_max_nrg(e, b, t)

qp_min_nrg(e, w, t)
qp_max_nrg(e, w, t)
;

$CALL GDXXRW Demands.xlsx par=q_max_prdct rng=Demands!B1:M112 rdim=2 cdim=1
$gdxin Demands.gdx
$load q_max_prdct
$gdxin
*q_min_prdct("Raw_Bg", m, t) = 150000;
*q_min_prdct("Biofert", m, t) = 50;
*q_min_prdct("SCP", m, t) = 20;
q_min_prdct("Raw_Bg", m, t) = 0;
q_min_prdct("Biofert", m, t) = 0;
q_min_prdct("SCP", m, t) = 0;

$CALL GDXXRW Demands.xlsx par=q_max_PBg rng=Demands!P31:Z68 rdim=1 cdim=1
$gdxin Demands.gdx
$load q_max_PBg
$gdxin
*q_min_PBg(m, t) = 30000;
q_min_PBg(m, t) = 0;

$CALL GDXXRW Demands.xlsx par=q_max_nrg rng=Demands!P15:AA27 rdim=2 cdim=1
$gdxin Demands.gdx
$load q_max_nrg
$gdxin
q_min_nrg(e, b, t) = 0;

$CALL GDXXRW Demands.xlsx par=qp_max_nrg rng=Demands!P1:AA12 rdim=2 cdim=1
$gdxin Demands.gdx
$load qp_max_nrg
$gdxin
qp_min_nrg(e, w, t) = 0;


********************************************************************************
* Environmental parameters: GHG emissions

Scalars
e_PBg_conv /0.0000125/
*e_Bg_CHP /4.4013/
e_Bg_CHP /0.0051/
*e_usg_PBg /2.7616/
e_usg_PBg /0.0032/
;

Parameters
e_Bm_src(j, i)
e_Bm_disp(j, i)
e_Bm_conv(j, k)
e_Bg_conv(h)
e_Fuel(s)

e_trans_Bm(j, l)
e_trans_prdct(n, lg)
e_trans_PBg(lg)
e_trans_Fuel(lg, s)

*e_usg_prdct(n) /Raw_Bg 4.4013, Biofert 0.0015, SCP 0.0015/
e_usg_prdct(n) /Raw_Bg 0.0051, Biofert 0.0015, SCP 0.0015/

e_LC_prdct(n, m)
e_LC_PBg(m)
ep_LC_PBg(p)
e_LC_BE(b)
e_LC_nrg(p)

beta_src(i)
beta_Bref(b)
beta_PowP(p)
beta_mrkt(m)
;

e_Bm_src(j, i) = 0.565;
e_Bm_disp(j, i) = 0.565;
e_Bm_conv(j, k) = 0.565;
e_Bg_conv(h) = 0.62;
e_Fuel(s) = 0.0011;

e_trans_Bm(j, l) = 0.594;
e_trans_prdct("Raw_Bg", lg) = 0.4;
e_trans_prdct("Biofert", lg) = 0.36;
e_trans_prdct("SCP", lg) = 0.36;
e_trans_PBg(lg) = 0.36;
e_trans_Fuel(lg, s) = 0.4;

e_LC_prdct(n, m) = 4;
e_LC_PBg(m) = 4;
e_LC_BE(b) = 0.01;
*e_LC_nrg(p) = 0.01;
*ep_LC_PBg(p) = 4;

beta_src(i) = 1;
beta_Bref(b) = 1;
beta_PowP(p) = 1;
beta_mrkt(m) = 1;

********************************************************************************
* Technical parameters: maximum production capacities and minimum capacity utilization rates
Parameter
ucap_Bm(j, i, t)

cap_tech_PBg(h) /PSA 912500, HPWS 912500, OPS 912500, CSP 912500, MS 912500/
capp_CHP(p, e)

delta_tech(k, n)
deltap_tech(h)
delta_PowP(p)
delta_CHP(c, e)
deltap_CHP(p, e)
;

Table cap_tech_prod(k, n)
         Raw_Bg       Biofert       SCP
CFDTBP   1268010      360           100
DFDTBH   1268010      360           100
FDTBP    1268010      360           100
;

Table cap_PowP(p, e)
                         electricity       thermal
hafez_power_plant        8514720000        8514720000
jahrom_power_plant       12649440000       12649440000
kazerun_power_plant      12027480000       12027480000
shiraz_power_plant       9198000000        9198000000
;

Table cap_CHP(c, e)
           electricity       thermal
small      87600000          87600000
medium     262800000         262800000
large      438000000         438000000
;

capp_CHP(p, "electricity") = 350400000;
capp_CHP(p, "thermal") = 350400000;

delta_tech(k, n) = 0.3;
deltap_tech(h) = 0.3;
delta_PowP(p) = 0.3;
delta_CHP(c, e) = 0.3;
deltap_CHP(p, e) = 0.3;

$CALL GDXXRW BSC.xlsx par=ucap_Bm rng=Sheet4!B3:M261 rdim=2 cdim=1
$gdxin BSC.gdx
$load ucap_Bm
$gdxin

********************************************************************************
* Other technical parameters: disposal rate of biomass and conversion factors
Parameters
varsigma_Bm(j, i)
kappa_Bm(e, j, k)
kappa_PBg(e, h)
kappa_Bg(c, e)
lambda_nrg(e)
lambdap_nrg(e)
mu_Fuel(s, e)
mup_Fuel(s, e)

theta_prdct(j, k, n)
phi_PBg(h) /PSA 0.75, HPWS 0.65, OPS 0.6, CSP 0.8, MS 0.8/
;

varsigma_Bm(j, i) = 0.1;

kappa_Bm(e, j, k) = 1000;
kappa_PBg(e, h) = 12945;
kappa_Bg(c, e) = 12945;

lambda_nrg(e) = 12945;
lambdap_nrg(e) = 12945;

mu_Fuel(s, e) = 10000;
mup_Fuel(s, e) = 10000;

theta_prdct(j, k, "Raw_Bg") = 0.5;
*theta_prdct(j, k, "Raw_Bg") = 0.579;
theta_prdct(j, k, "Biofert") = 0.05;
theta_prdct(j, k, "SCP") = 0.01;

********************************************************************************
* Energy requirements for productions and processes
Parameters
er_Bm(e, j, k)
er_Bg(e, h)
er_Bg_CHP(e, c)
er_PBg(e, p)
er_PBg_CHP(e, p)
er_Fuel(e, s, p)
er_Fuel_CHP(e, s, p)
;

er_Bm("electricity", j, k) = 26.78;
er_Bm("thermal", j, k) = 85.69;

*er_Bg(e, "PSA") = 0.335;
*er_Bg(e, "HPWS") = 0.43;
*er_Bg(e, "OPS") = 0.49;
*er_Bg(e, "CSP") = 0.646;
**er_Bg(e, "MS") = 0.769;

*er_PBg(e, p) = 0.5;
*er_Bg_CHP(e, c) = 0.5;
*er_PBg_CHP(e, p) = 0.5;

er_Bg(e, "PSA") = 289.105;
er_Bg(e, "HPWS") = 371.09;
er_Bg(e, "OPS") = 422.87;
er_Bg(e, "CSP") = 557.498;
er_Bg(e, "MS") = 663.647;

er_PBg(e, p) = 431.5;
er_Bg_CHP(e, c) = 431.5;
er_PBg_CHP(e, p) = 431.5;

*er_Fuel(e, s, p) = 20;
*er_Fuel_CHP(e, s, p) = 20;

er_Fuel(e, s, p) = 1000;
er_Fuel_CHP(e, s, p) = 1000;

********************************************************************************
* Working hours requirements

Parameters
w_Bref(bp, k, h, c, v)
w_pipe(v) /Class1 0, Class2 1035, Class3 1725/
w_src(j, i, v)
w_prod(j, n, k, v)
w_prod_PBg(h, v)
w_CHP(b, c, v)
w_PowP(p, v)
w_Fuel(s, p, v)
wp_CHP(p, v)
wpp_CHP(p, s, v)
;

$CALL GDXXRW BSC2.xlsx par=w_Bref rng=Sheet4!A1:G271 rdim=4 cdim=1
$gdxin BSC2.gdx
$load w_Bref
$gdxin

w_src("Chicken", i, "Class1") = 0.5;
w_src("Cow", i, "Class1") = 0.5;
w_src("Horse", i, "Class1") = 0.5;
w_src("Agri_waste", i, "Class1") = 0.5;
w_src("Urban", i, "Class1") = 0.5;
w_src("Urban_ready", i, "Class1") = 0.5;
w_src("Chicken", i, "Class2") = 1;
w_src("Cow", i, "Class2") = 1;
w_src("Horse", i, "Class2") = 1;
w_src("Agri_waste", i, "Class2") = 0;
w_src("Urban", i, "Class2") = 1;
w_src("Urban_ready", i, "Class2") = 1;
w_src("Chicken", i, "Class3") = 2.5;
w_src("Cow", i, "Class3") = 10;
w_src("Horse", i, "Class3") = 5;
w_src("Agri_waste", i, "Class3") = 7.5;
w_src("Urban", i, "Class3") = 7.5;
w_src("Urban_ready", i, "Class3") = 7.5;

w_prod(j, "Raw_Bg", k, "Class1") = 200;
w_prod(j, "Raw_Bg", k, "Class1") = 240;
w_prod(j, "Raw_Bg", k, "Class1") = 120;
w_prod(j, "Biofert", k, "Class2") = 200;
w_prod(j, "Biofert", k, "Class2") = 240;
w_prod(j, "Biofert", k, "Class2") = 120;
w_prod(j, "SCP", k, "Class3") = 200;
w_prod(j, "SCP", k, "Class3") = 240;
w_prod(j, "SCP", k, "Class3") = 120;

w_prod_PBg(h, "Class1") = 200;
w_prod_PBg(h, "Class2") = 240;
w_prod_PBg(h, "Class3") = 120;

w_CHP(b, c, "Class2") = 10;
w_CHP(b, c, "Class3") = 10;

w_PowP(p, "Class1") = 100;
w_PowP(p, "Class2") = 120;
w_PowP(p, "Class3") = 60;

w_Fuel(s, p, "Class1") = 100;
w_Fuel(s, p, "Class2") = 120;
w_Fuel(s, p, "Class3") = 60;

wp_CHP(p, "Class2") = 10;
wp_CHP(p, "Class3") = 10;

wpp_CHP(p, s, "Class2") = 10;
wpp_CHP(p, s, "Class3") = 10;

********************************************************************************
* Prices: "p" in the beginning of parameters below= Producer's perspective, "c" in the beginning of parameters below= Consumer's perspective

Parameters
c_cont(i, j, t)
c_Fuel(s, p, t)
c_PBg(b, p, t)
c_prdct(n, b, t)
cp_PBg(b, t)
c_nrg(e,b,t)
cp_nrg(e, p, t)
p_cont(i, j, t)
p_Fuel(s, p, t)
p_PBg(b, p, t)
p_prdct(n, b, t)
pp_PBg(b, t)
p_nrg(e,b,t)
pp_nrg(e, p, t)
;

c_cont(i, j, t) = 0;
p_cont(i, j, t) = 0;
c_cont(ipp, "Urban", t) = 1;
p_cont(ipp, "Urban", t) = 1;
c_cont(ip, "Urban_ready", t) = 1.5;
p_cont(ip, "Urban_ready", t) = 1.5;

*c_prdct("Raw_Bg", b, t) = 474.65;
*p_prdct("Raw_Bg", b, t) = 474.65;
*c_prdct("Raw_Bg", b, t) = 0.7;
*p_prdct("Raw_Bg", b, t) = 0.7;
p_prdct("Raw_Bg", b, t) = 8630;
p_prdct("Biofert", b, t) = 1500;
p_prdct("SCP", b, t) = 5000;

c_Fuel("shiraz_oil_refinery", p, t) = 560;
c_Fuel("farashband_gas_refinery", p, t) = 560;
p_Fuel("shiraz_oil_refinery", p, t) = 560;
p_Fuel("farashband_gas_refinery", p, t) = 560;

*cp_PBg(b, t) = 1.3;
*pp_PBg(b, t) = 1.3;

pp_PBg(b, t) = 9500;

*c_PBg(b, p, t) = 1.05;
*p_PBg(b, p, t) = 1.05;

c_PBg(b, p, t) = 9000;
p_PBg(b, p, t) = 9000;

c_nrg(e, b, t) = 1;
p_nrg(e, b, t) = 1;

cp_nrg(e, p, t) = 1;
pp_nrg(e, p, t) = 1;

********************************************************************************
Parameter
c_purch_nrg(e,t)
cp_purch_nrg(e,t);
c_purch_nrg(e,t) = 1;
cp_purch_nrg(e,t) = 1;
********************************************************************************

Parameter
discount_fraction(t);

loop(t$(ord(t) <= horizon),
discount_fraction(t) = (1/(POWER(1+alpha,(ord(t)-1))));
);

*////////////////////////////////////////////////////////////////////////////////

Variables
NPV
GrandTotalRevenues
GrandTotalRevenues_Bref
GrandTotalRevenues_PowP
GrandTotalCosts_Bref
GrandTotalCosts_PowP
GrandTotalCosts

TotalRevenues(t)
Brefs_Revenues(t)
PowPs_Revenues(t)

TotalCosts(t)
Brefs_costs(t)
dstrbs_costs(t)
PowPs_costs(t)

****************************************

Net_GHG_emission_savings
Avoided_GHG_emissions
Total_GHG_emissions_of_the_supply_chain

****************************************

Social_Benefit
Constructions
Transfers
Productions
;

Binary Variables
Z_Bref(b, k, h, c)
Z_dstrb(d, u)
Z_src_Bref(i, j, b, l)
Z_Bref_PowP(b, p, lg)
Z_Bref_dstrb(b, d, lg)
Z_dstrb_mrkt(d, m, lg)
;

Positive Variables
* Number of trucks /trains fleet required
FL_Bm(l, j, i, b, t)
FL_prdct(l, n, b, m, t)
FL_PBg(l, b, p, t)
FLp_PBg(l, b, m, t)
FL_Fuel(l, s, p, t)

Positive Variables
* Energy requirements
R_Bref(e, b, t)
R_PowP(e, p, t)

* Energy that must be supplied
E_Bref(e, b, t)
E_PowP(e, p, t)

* Productions
*X_Bm(j, i, t)
XT_Bm(j, i, b, l, k, t)
X_prdct(n, j, b, k, t)
X_PBg(b, h, t)
X_nrg(e, b, t)
Xp_nrg(e, p, t)

* Available amounts
*Y_Bm_conv(j, b, k, f, t)
Y_Bg_conv(b, h, t)
Y_Bg(b, c, t)
Y_PBg_conv(p, t)
Y_Fuel_conv(s, p, t)
Y_PBg_CHP(p, t)
Y_Fuel_CHP(s, p, t)

* Transfers
*T_trans_Bm(j, i, b, l, t)
*T_trans_prdct(n, b, d, lg, t)
T_trans_PBg(b, p, lg, t)
*Tp_trans_PBg(b, d, lg, t)
TS_prdct(n, b, m, lg, t)
TS_PBg(b, m, lg, t)

* Sales
S_prdct(n, d, m, t)
S_PBg(d, m, t)
S_nrg(e, b, t)
Sp_nrg(e, p, w, t)

* Fossil fuel purchases
V_Fuel(s, p, l, t)
;

Variable
OOF

Positive Variable
Surplus1
Surplus2
Surplus3
;

********************************************************************************

Equations
OOF_obj

NPV_obj, TotalRevenues_obj, Brefs_Revenues_obj, PowPs_Revenues_obj
TotalCosts_obj, Brefs_costs_obj, PowPs_costs_obj
GrandTotalRevenues_obj, GrandTotalCosts_obj
GrandTotalRevenues_Bref_obj, GrandTotalRevenues_PowP_obj
GrandTotalCosts_Bref_obj, GrandTotalCosts_PowP_obj
Net_GHG_emission_savings_obj, Avoided_GHG_emissions_obj, Total_GHG_emissions_of_the_supply_chain_obj
Social_Benefit_obj, Constructions_obj, Productions_obj, Transfers_obj

eps_co1, eps_co2,
co1,co2,co5,co6,co7,co8,co9,co10,co11,co12,co13,co15,co16,co18,co20,co21,co22_1,co22_2,co23,co24,co25,co26,co28,
co29,co30,co31_1,co31_2,co32_1,co32_2,co33_1,co33_2,co34,co35,co36,co37,co38,co39,co40,co41,co42,co43,co44
co_pipe
;

********************************************************************************

OOF_obj                     .. OOF =e= NPV + EPSILON * (Surplus2 + 0.1 * Surplus3);
*OOF_obj                     .. OOF =e= Avoided_GHG_emissions + EPSILON * (Surplus1 + 0.1 * Surplus3);

*$ontext

NPV_obj                     .. NPV =e= Sum((t),discount_fraction(t)*(TotalRevenues(t)-TotalCosts(t)));
*GrandTotalRevenues_obj     .. GrandTotalRevenues =e= Sum(t,Brefs_Revenues(t));
*GrandTotalCosts_obj        .. GrandTotalCosts =e= Sum(t,Brefs_costs(t));
*GrandTotalRevenues_Bref_obj .. GrandTotalRevenues_Bref =e= Sum(t,Brefs_Revenues(t));
*GrandTotalRevenues_PowP_obj .. GrandTotalRevenues_PowP =e= Sum(t,PowPs_Revenues(t));
*GrandTotalCosts_Bref_obj    .. GrandTotalCosts_Bref =e= Sum(t,Brefs_costs(t));
*GrandTotalCosts_PowP_obj    .. GrandTotalCosts_PowP =e= Sum(t,PowPs_costs(t));

*TotalRevenues_obj(t)        .. TotalRevenues(t) =e= Brefs_Revenues(t)+PowPs_Revenues(t);
TotalRevenues_obj(t)        .. TotalRevenues(t) =e= Brefs_Revenues(t);

Brefs_Revenues_obj(t)       .. Brefs_Revenues(t) =e= Sum((n,b,m,lg),(p_prdct(n,b,t)*TS_prdct(n,b,m,lg,t)+pp_PBg(b,t)*TS_PBg(b,m,lg,t)))
*                               + Sum((b,p,lg),p_PBg(b,p,t)*T_trans_PBg(b,p,lg,t))
                               + Sum((b,e),p_nrg(e,b,t)*S_nrg(e,b,t));

*PowPs_Revenues_obj(t)       .. PowPs_Revenues(t) =e= Sum((e,p,w),pp_nrg(e,p,t)*Sp_nrg(e,p,w,t));

*TotalCosts_obj(t)           .. TotalCosts(t) =e= Brefs_costs(t) + PowPs_costs(t);
TotalCosts_obj(t)           .. TotalCosts(t) =e= Brefs_costs(t);

Brefs_costs_obj(t)          .. Brefs_costs(t) =e= Sum((b,h,c,k),(sigma(b)*c_invst(k,h,c)+c_maint_Bref(b,k,h,c,t))*Z_Bref(b,k,h,c))
                               + Sum((i,jp,b,lp),c_pipe*d_i_b(i,b)*sigma_pipe(i,b)*Z_src_Bref(i,jp,b,lp))
                               + Sum((j,i,b,l,k),(c_purch_Bm(j,i,t)+c_cont(i,j,t)+c_water(t)*wr(j,k)+c_pre_Bm(j,b,t))*XT_Bm(j,i,b,l,k,t))
                               + Sum((jp,i,b,lp,k),c_trans_Bm(jp,lp,t)*XT_Bm(jp,i,b,lp,k,t))
                               + Sum((i,j,b,lg),c_trans_Bm(j,lg,t)*FL_Bm(lg,j,i,b,t)*d_i_b(i,b))
                               + Sum((n,j,b,h,c,k),(c_prod(n,k,t)*X_prdct(n,j,b,k,t)+c_PBg_prod(h,t)*X_PBg(b, h, t)+c_Bg_CHP(c,t)*Y_Bg(b,c,t)))
                               + Sum((n,b,m,lg),c_trans_prdct(n,lg,t)*FL_prdct(lg,n,b,m,t)*d_b_m(b,m))
                               + Sum((b,m,lg),c_trans_PBg(lg,t)*FLp_PBg(lg,b,m,t)*d_b_m(b,m))
                               + Sum((e,b),c_purch_nrg(e,t)*E_Bref(e,b,t));

*PowPs_costs_obj(t)          .. PowPs_costs(t) =e= Sum((b,p,s,lg,e1),c_maint_PowP(p,t)+c_PBg(b,p,t)*T_trans_PBg(b,p,lg,t)+c_Fuel(s,p,t)*V_Fuel(s,p,lg,t)+c_elec(p,t)*Xp_nrg(e1,p,t)+c_PBg_CHP(p,t)*Y_PBg_CHP(p,t)+c_Fuel_CHP(s,p,t)*Y_Fuel_CHP(s,p,t))
*                               + Sum((b,p,lg),c_trans_PBg(lg,t)*FL_PBg(lg,b,p,t)*d_b_p(b,p))
*                               + Sum((lg,p,s),c_trans_Fuel(s,lg,t)*FL_Fuel(lg,s,p,t)*d_s_p(s,p))
*                               + Sum((e,p),cp_purch_nrg(e,t)*E_PowP(e,p,t));

*$offtext



*$ontext


Net_GHG_emission_savings_obj   .. Net_GHG_emission_savings =e= Avoided_GHG_emissions - Total_GHG_emissions_of_the_supply_chain;
Avoided_GHG_emissions_obj .. Avoided_GHG_emissions =e= Sum((b,i,j,k,l,t),beta_src(i)*e_Bm_disp(j,i)*varsigma_Bm(j,i)*XT_Bm(j,i,b,l,k,t))
                                  + Sum((b,lg,m,n,t),beta_mrkt(m)*(e_LC_prdct(n,m)*TS_prdct(n,b,m,lg,t)+e_LC_PBg(m)*TS_PBg(b,m,lg,t)))
                                  + Sum((b,e,t),beta_Bref(b)*e_LC_BE(b)*S_nrg(e,b,t));
*                                 + Sum((m,e,p,t),beta_PowP(p)*e_LC_nrg(p)*(lambda_nrg(e)*Y_PBg_conv(p,t)+lambdap_nrg(e)*Y_PBg_CHP(p,t)));

Total_GHG_emissions_of_the_supply_chain_obj .. Total_GHG_emissions_of_the_supply_chain =e= Sum((b,i,j,k,l,t),(beta_src(i)*e_Bm_src(j,i)+e_trans_Bm(j,l)*d_i_b(i,b))*XT_Bm(j,i,b,l,k,t))
                                                    + Sum((b,c,i,j,h,k,lg,t),beta_Bref(b)*(e_Bm_conv(j,k)*XT_Bm(j,i,b,lg,k,t)+e_Bg_conv(h)*Y_Bg_conv(b,h,t)+e_Bg_CHP*Y_Bg(b,c,t)))
                                                    + Sum((b,lg,m,n,t),(e_trans_prdct(n,lg)*TS_prdct(n,b,m,lg,t)+e_trans_PBg(lg)*TS_PBg(b,m,lg,t))*d_b_m(b,m))
*                                                    + Sum((b,lg,p,s,t),e_trans_PBg(lg)*d_b_p(b,p)*T_trans_PBg(b,p,lg,t)+e_trans_Fuel(lg,s)*d_s_p(s,p)*V_Fuel(s,p,lg,t))
*                                                    + Sum((e,p,s,t),beta_PowP(p)*(e_PBg_conv*lambda_nrg(e)*Y_PBg_conv(p,t)+e_Fuel(s)*mu_Fuel(s,e)*Y_Fuel_conv(s,p,t)))
                                                    + Sum((b,lg,m,n,t),beta_mrkt(m)*(e_usg_prdct(n)*TS_prdct(n,b,m,lg,t)+e_usg_PBg*TS_PBg(b,m,lg,t)))
                                                    + Sum((b,e,t),e_LC_BE(b)*E_Bref(e,b,t));
*                                                    + Sum((e,p,t),e_LC_nrg(p)*E_PowP(e,p,t));

*$offtext



*$ontext

Social_Benefit_obj   .. Social_Benefit =e= Constructions + Productions + Transfers;

Constructions_obj .. Constructions =e= Sum((bp,c,h,k,v),w_Bref(bp,k,h,c,v)*Z_Bref(bp,k,h,c))
                          + Sum((b,i,jp,lp,v),w_pipe(v)*d_i_b(i,b)*Z_src_Bref(i,jp,b,lp));

Productions_obj .. Productions =e= Sum((b,i,j,l,k,t,v),w_src(j,i,v)*XT_Bm(j,i,b,l,k,t))
                        + Sum((b,c,t,v),(Sum((j,k,n),w_prod(j,n,k,v)*X_prdct(n,j,b,k,t))+Sum((h),w_prod_PBg(h,v)*X_PBg(b,h,t))+w_CHP(b,c,v)*Y_Bg(b,c,t)));
*                        + Sum((p,lg,s,t,v),(w_PowP(p,v)*Y_PBg_conv(p,t)+w_Fuel(s,p,v)*V_Fuel(s,p,lg,t))+wp_CHP(p,v)*Y_PBg_CHP(p,t)+wpp_CHP(p,s,v)*Y_Fuel_CHP(s,p,t));

*Transfers_obj .. Transfers =e= Sum((lg,t,v),time(lg)*(Sum((b,i,j),d_i_b(i,b)*FL_Bm(lg,j,i,b,t)*chi_trans_i_b(i,b,lg,v))+Sum((b,p),d_b_p(b,p)*FL_PBg(lg,b,p,t)*chi_trans_b_p(b,p,lg,v))
*                      + Sum((b,m,n),d_b_m(b,m)*(FL_prdct(lg,n,b,m,t)+FLp_PBg(lg,b,m,t))*chi_trans_b_m(b,m,lg,v))
*                      + Sum((p,s),d_s_p(s,p)*FL_Fuel(lg,s,p,t)*chi_trans_s_p(s,p,lg,v))));

*Transfers_obj .. Transfers =e= Sum((lg,t,v),time(lg)*(Sum((b,i,j,k),d_i_b(i,b)*XT_Bm(j,i,b,lg,k,t)*chi_trans_i_b(i,b,lg,v))+Sum((b,p),d_b_p(b,p)*T_trans_PBg(b,p,lg,t)*chi_trans_b_p(b,p,lg,v))
*                      + Sum((b,m,n),d_b_m(b,m)*(TS_prdct(n,b,m,lg,t)+TS_PBg(b,m,lg,t))*chi_trans_b_m(b,m,lg,v))
*                      + Sum((p,s),d_s_p(s,p)*V_Fuel(s,p,lg,t)*chi_trans_s_p(s,p,lg,v))));

Transfers_obj .. Transfers =e= Sum((lg,t,v),time(lg)*(Sum((b,i,j,k),d_i_b(i,b)*(XT_Bm(j,i,b,lg,k,t)/lc_Bm(lg,j))*chi_trans_i_b(i,b,lg,v))
*                              + Sum((b,p),d_b_p(b,p)*(T_trans_PBg(b,p,lg,t)/lc_PBg(lg))*chi_trans_b_p(b,p,lg,v))
                              + Sum((b,m,n),d_b_m(b,m)*((TS_prdct(n,b,m,lg,t)/lc_prdct(lg,n))+(TS_PBg(b,m,lg,t)/lc_PBg(lg)))*chi_trans_b_m(b,m,lg,v))));
*                             + Sum((p,s),d_s_p(s,p)*(V_Fuel(s,p,lg,t)/lc_Fuel(lg,s))*chi_trans_s_p(s,p,lg,v))));

*$offtext

****************************************************************************************

*eps_co1        ..  NPV - Surplus1 =e= epsilon1;
eps_co1        ..  Avoided_GHG_emissions - Surplus2 =e= epsilon2;
eps_co2        ..  Social_Benefit - Surplus3 =e= epsilon3;

****************************************************************************************

co1(i,j,t)                .. Sum((b,l,k),XT_Bm(j,i,b,l,k,t)*(1-varsigma_Bm(j,i))) =l= ucap_Bm(j,i,t)*z_src(i,j);
co2(b,i,j,lg,t,k)           .. XT_Bm(j,i,b,lg,k,t) =l= lc_Bm(lg,j)*FL_Bm(lg,j,i,b,t);
co5(jp,i,b)             .. Sum((lp,k,t),XT_Bm(jp,i,b,lp,k,t)) =l= wbd(i,b)*z_src(i,jp)*bigM;
co6(jnp,i,b,lp,k,t)        .. XT_Bm(jnp,i,b,lp,k,t) =e= 0;
co7(b,k,j)                    .. Sum((i,l,t),XT_Bm(j,i,b,l,k,t)) =l= Sum((h,c),Z_Bref(b,k,h,c)*tcb(k,j)*bigM);
co8(b,j,n,t,k)              .. X_prdct(n,j,b,k,t) =e= Sum((i,l),theta_prdct(j,k,n)*XT_Bm(j,i,b,l,k,t));
co9(b,t,h)                  .. X_PBg(b,h,t) =e= phi_PBg(h)*Y_Bg_conv(b,h,t);
co10(b,n,t)                .. Sum((j,k),X_prdct("Raw_Bg",j,b,k,t)) =e= Sum(h,Y_Bg_conv(b,h,t))+Sum(c,Y_Bg(b,c,t))+Sum((m,lg),TS_prdct("Raw_Bg",b,m,lg,t));
co11(b,n,t)               .. Sum(h,Y_Bg_conv(b,h,t))-Sum(c,Y_Bg(b,c,t)) =e= Sum((m,lg),TS_prdct("Raw_Bg",b,m,lg,t));
co12(b,np,t)              .. Sum((m,lg),TS_prdct(np,b,m,lg,t)) =e= Sum((j,k),X_prdct(np,j,b,k,t));
co13(b,m,lg,n,t)          .. TS_prdct(n,b,m,lg,t) =l= lc_prdct(lg,n)*FL_prdct(lg,n,b,m,t);
*co15(b,t)                 .. Sum(h,X_PBg(b,h,t)) =e= Sum((lg,p),T_trans_PBg(b,p,lg,t)) + Sum((lg,m),TS_PBg(b,m,lg,t));
co15(b,t)                 .. Sum(h,X_PBg(b,h,t)) =e= Sum((lg,m),TS_PBg(b,m,lg,t));
*co16(b,lg,p,t)            .. T_trans_PBg(b,p,lg,t) =l= lc_PBg(lg)*FL_PBg(lg,b,p,t);
co18(b,m,lg,t)            .. TS_PBg(b,m,lg,t) =l= lc_PBg(lg)*FLp_PBg(lg,b,m,t);
co20(b,e,t)               .. X_nrg(e,b,t) =e= Sum((j,i,l,k),kappa_Bm(e,j,k)*XT_Bm(j,i,b,l,k,t))+Sum(h,kappa_PBg(e,h)*Y_Bg_conv(b,h,t))+Sum(c,kappa_Bg(c,e)*Y_Bg(b,c,t));
co21(b,e,t)               .. X_nrg(e,b,t) =e= R_Bref(e,b,t)+S_nrg(e,b,t);
co22_1(b,e,t)             .. S_nrg(e,b,t) =l= q_max_nrg(e,b,t);
co22_2(b,e,t)             .. S_nrg(e,b,t) =g= q_min_nrg(e,b,t);
co23(b,e,t)               .. Sum((j,i,l,k),er_Bm(e,j,k)*XT_Bm(j,i,b,l,k,t))+Sum(h,er_Bg(e,h)*Y_Bg_conv(b,h,t))+Sum(c,er_Bg_CHP(e,c)*Y_Bg(b,c,t)) =e= R_Bref(e,b,t)+E_Bref(e,b,t);
*co24(p,t)                 .. Y_PBg_conv(p,t)+Y_PBg_CHP(p,t) =e= Sum((b,lg),T_trans_PBg(b,p,lg,t));
*co25(p,s,t)               .. Y_Fuel_conv(s,p,t)+Y_Fuel_CHP(s,p,t) =e= Sum(lg,V_Fuel(s,p,lg,t));
*co26(lg,p,s,t)            .. V_Fuel(s,p,lg,t) =l= lc_Fuel(lg,s)*FL_Fuel(lg,s,p,t);
*co28(e,p,t)               .. Xp_nrg(e,p,t) =e= lambda_nrg(e)*Y_PBg_conv(p,t)+Sum((s),mu_Fuel(s,e)*Y_Fuel_conv(s,p,t))+lambdap_nrg(e)*Y_PBg_CHP(p,t)+Sum((s),mup_Fuel(s,e)*Y_Fuel_CHP(s,p,t));
*co29(e,p,t)               .. Xp_nrg(e,p,t) =e= R_PowP(e,p,t)+Sum((w),Sp_nrg(e,p,w,t));
*co30(e,p,t)               .. er_PBg(e,p)*Y_PBg_conv(p,t)+Sum((s),er_Fuel(e,s,p)*Y_Fuel_conv(s,p,t))+er_PBg_CHP(e,p)*Y_PBg_CHP(p,t)+Sum((s),er_Fuel_CHP(e,s,p)*Y_Fuel_CHP(s,p,t)) =e= R_PowP(e,p,t)+E_PowP(e,p,t);
co31_1(m,n,t)             .. Sum((b,lg),TS_prdct(n,b,m,lg,t)) =l= q_max_prdct(n,m,t);
co31_2(m,n,t)             .. Sum((b,lg),TS_prdct(n,b,m,lg,t)) =g= q_min_prdct(n,m,t);
co32_1(m,t)               .. Sum((b,lg),TS_PBg(b,m,lg,t)) =l= q_max_PBg(m,t);
co32_2(m,t)               .. Sum((b,lg),TS_PBg(b,m,lg,t)) =g= q_min_PBg(m,t);
*co33_1(e,t,w)             .. Sum((p),Sp_nrg(e,p,w,t)) =l= qp_max_nrg(e,w,t);
*co33_2(e,t,w)             .. Sum((p),Sp_nrg(e,p,w,t)) =g= qp_min_nrg(e,w,t);

co34(b,n,k)             .. Sum((j,t),X_prdct(n,j,b,k,t)) =l= horizon*Sum((h,c),cap_tech_prod(k,n)*Z_Bref(b,k,h,c));
*co35(b,n,k)                 .. Sum((j,t),X_prdct(n,j,b,k,t)) =g= horizon*Sum((h,c),cap_tech_prod(k,n)*delta_tech(k,n)*Z_Bref(b,k,h,c));

co36(b,h)                   .. Sum((t),X_PBg(b,h,t)) =l= horizon*Sum((c,k),1.158*cap_tech_PBg(h)*Z_Bref(b,k,h,c));
*co37(b,h)                   .. Sum((t),X_PBg(b,h,t)) =g= horizon*Sum((c,k),cap_tech_PBg(h)*deltap_tech(h)*Z_Bref(b,k,h,c));

co38(b,e)                 .. Sum((t),X_nrg(e,b,t)) =l= horizon*Sum((k,h,c),cap_CHP(c,e)*Z_Bref(b,k,h,c));
*co39(b,e)                 .. Sum((t),X_nrg(e,b,t)) =g= horizon*Sum((k,h,c),cap_CHP(c,e)*delta_CHP(c,e)*Z_Bref(b,k,h,c));

*co40(e,p)                 .. Sum((t),Xp_nrg(e,p,t)) =l= horizon*(cap_PowP(p,e)+capp_CHP(p,e));
*co41(e,p)                 .. Sum((t),Xp_nrg(e,p,t)) =g= horizon*(cap_PowP(p,e)*delta_PowP(p)+capp_CHP(p,e)*deltap_CHP(p,e));

co_pipe(jp,i,b,lp,k,t)    .. XT_Bm(jp,i,b,lp,k,t) =l= cap_pipe;

co42(b)                  .. Sum((k,h,c),Z_Bref(b,k,h,c)) =l= 1;
co43(b)                .. Sum((i,jp,lp),Z_src_Bref(i,jp,b,lp)) =l= Sum((k,h,c),Z_Bref(b,k,h,c));
co44(lp,b,jp,i,k)         .. Sum(t,XT_Bm(jp,i,b,lp,k,t)) =l= Z_src_Bref(i,jp,b,lp)*bigM;


********************************************************************************

*Option reslim=2000, optca=0.1, optcr=0.1, mip=BARON;
Option optca=0, optcr=0, MIP=CPLEX;
Model BSC /
OOF_obj,

*GrandTotalRevenues_obj, GrandTotalCosts_obj,
*GrandTotalRevenues_Bref_obj,
*GrandTotalRevenues_PowP_obj,
*GrandTotalCosts_Bref_obj,
*GrandTotalCosts_PowP_obj,
NPV_obj, TotalRevenues_obj, Brefs_Revenues_obj,
*PowPs_Revenues_obj
TotalCosts_obj, Brefs_costs_obj,
*PowPs_costs_obj,

Net_GHG_emission_savings_obj,
Avoided_GHG_emissions_obj,
Total_GHG_emissions_of_the_supply_chain_obj

Social_Benefit_obj,
Constructions_obj,
Productions_obj,
Transfers_obj

eps_co1
eps_co2

co1,
co5,
co6,
co7,
co8,
co9,
co10,
co11,
co12,
co15,
co20,
co21,
co22_1,
co22_2,
co23,
*co24,
*co25,
*co28,
*co29,
*co30,
co31_1,
co31_2,
co32_1,
co32_2,
*co33_1,
*co33_2,
co34,
co36,
co38,
*co40,
co42,
co43,
co44
co_pipe

*co35,
*co37,
*co39,
*co41,

co2,
co13,
*co16,
co18,
*co26
/
;



loop(counter2$(ord(counter2) = q2),
    epsilon2 = LowerBound2 + ord(counter2)*step2;
    loopCount2 = ord(counter2);

*    loop(counter3$(ord(counter3) <= q3),
    loop(counter3$(ord(counter3) = q3),
        epsilon3 = LowerBound3 + ord(counter3)*step3;
        loopCount3 = ord(counter3);

Solve BSC using MIP max OOF;
*Solve BSC using MIP max NPV;
*Solve BSC using MIP max Avoided_GHG_emissions;
*Solve BSC using MIP max Social_Benefit;

Display
loopCount2, loopCount3,
*NPV.l, Avoided_GHG_emissions.l, Social_Benefit.l





*$ontext

*Display
discount_fraction, wbd,
NPV.l,
*GrandTotalRevenues_Bref.l,
*GrandTotalRevenues_PowP.l,
*GrandTotalCosts_Bref.l,
*GrandTotalCosts_PowP.l,
*GrandTotalCosts.l,
TotalRevenues.l,
TotalCosts.l,
Brefs_Revenues.l,
*PowPs_Revenues.l,
Brefs_costs.l,
*PowPs_costs.l,

Net_GHG_emission_savings.l,
Avoided_GHG_emissions.l,
Total_GHG_emissions_of_the_supply_chain.l,

Social_Benefit.l,
Constructions.l,
Productions.l,
Transfers.l,

Z_Bref.l,
Z_src_Bref.l,
R_Bref.l, E_Bref.l,
*R_PowP.l, E_PowP.l,
XT_Bm.l, X_prdct.l, X_PBg.l,
X_nrg.l,
*Xp_nrg.l,
Y_Bg_conv.l, Y_Bg.l,
*Y_PBg_conv.l,
*Y_Fuel_conv.l,
*Y_PBg_CHP.l, Y_Fuel_CHP.l,
*T_trans_PBg.l,
*V_Fuel.l,
TS_prdct.l, TS_PBg.l,
S_nrg.l,
*Sp_nrg.l
FL_Bm.l,
FL_prdct.l,
*FL_PBg.l,
FLp_PBg.l
*FL_Fuel.l
;

*$offtext



);
);









