/*  q type constants used by high precision check routines */
/* Verified using mathematica and pari GP. jacob */

#include "qhead.h"

/*                                     -1.0 */
Qfloat qminusone[1] = {
{-1,EXPONE,{0x8000000000000000ULL,0,0,0,0,0,0}}};
/*                                    0.0E0 */
Qfloat qzero[1] = {{0}};

/* 5.0E-1 */
Qfloat qhalf[1] = {
{0,EXPONE-1,{0x8000000000000000ULL,0,0,0,0,0,0}}};

/* 1.0E0 */
Qfloat qone[1] = { {0,EXPONE,{0x8000000000000000ULL,0,0,0,0,0,0}}};

/* 2.0E0 */
Qfloat qtwo[1] = {
{0,EXPONE+1,{0x8000000000000000ULL,0,0,0,0,0,0}}};

/* 3.0E0 */
Qfloat qthree[1] = {
{0,EXPONE+1,{0xc000000000000000ULL,0,0,0,0,0,0}}};

/* -3.0E0 */
Qfloat qminus_three[1] = {
{-1,EXPONE+1,{0xc000000000000000ULL,0,0,0,0,0,0}}};

/* 5.0E0 */
Qfloat qfive[1] = {{0,EXPONE+2,{0xa000000000000000ULL,0,0,0,0,0,0,}}};
/* 9.0e0 */
Qfloat qnine[1] = {
{0,EXPONE+3,{0x9000000000000000ULL,0,0,0,0,0,0}}};

/* 3.2E1 */
Qfloat q32[1] = {
{0,EXPONE+5,{0x8000000000000000ULL,0,0,0,0,0,0}}};

/* 1/3 */
Qfloat oneThird[1] = {
{0,0x7ffff,{0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,
0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,
0xaaaaaaaaaaaaaaabULL}}};

/*6.9314718055994530941723212145817656807550013436025525412068000949339362E-1*/
Qfloat qlog2[1] = {
{0,EXPONE-1,{0xb17217f7d1cf79abULL,0xc9e3b39803f2f6afULL,
0x40f343267298b62dULL,0x8a0d175b8baafa2bULL,0xe7b876206debac98ULL,
0x559552fb4afa1b10ULL,0xed2eae35c1382144ULL}}};

Qfloat qinv_log2[1]={ // Rounded
{0,0x00080001,{0xb8aa3b295c17f0bbULL,0xbe87fed0691d3e88ULL,0xeb577aa8dd695a58ULL,
0x8b25166cd1a13247ULL, 0xde1c43f755176cd6ULL,0x24d92f75c16be0b3ULL,0xea90b9e60c4a909fULL}}};

/*  1.41421356237309504880168872420969807856967187537694
 *    80731766797379907324784621070388503875343276415727e0
 */
Qfloat qsqrt2[1] = {
{0,EXPONE,{0xb504f333f9de6484ULL,0x597d89b3754abe9fULL,
0x1d6f60ba893ba84cULL,0xed17ac8583339915ULL,0x4afc83043ab8a2c3ULL,
0xa8b1fe6fdc83db39ULL,0x0f74a85e439c7b4aULL}}};

Qfloat qinv_sqrt2[1] = {
{0,EXPONE-1,{0xb504f333f9de6484ULL, 0x597d89b3754abe9fULL, 
0x1d6f60ba893ba84cULL, 0xed17ac8583339915ULL, 0x4afc83043ab8a2c3ULL, 
0xa8b1fe6fdc83db39ULL, 0x0f74a85e439c7b4a}}};

/* 2/sqrt(PI) =
 *1.12837916709551257389615890312154517168810125865799771368817144342128493688298682897348732040421472688605669581272341470337986298965E0*/
Qfloat oneopi[1] = {
{0,EXPONE,{0x906eba8214db688dULL,0x71d48a7f6bfec344ULL,
0x1409a0ebac3e7517ULL,0x39a15830cce620b0ULL,0xc0759cf859270f11ULL,
0x40c036096cc79aebULL,0xbd1f4eee48e1ca79ULL}}};

/* 3.14159265358979323846264338327950288419716939937510
 *  582097494459230781640628620899862803482534211706798e0
 */
Qfloat qpi[1] = {
{0,EXPONE+1,{0xc90fdaa22168c234ULL,0xc4c6628b80dc1cd1ULL,
0x29024e088a67cc74ULL,0x020bbea63b139b22ULL,0x514a08798e3404ddULL,
0xef9519b3cd3a431bULL,0x302b0a6df25f1437ULL}}};

QfloatAccum qpiAccum[1]={
{0,EXPONE+1,{0xc90fdaa22168c234ULL,0xc4c6628b80dc1cd1ULL,
 0x29024e088a67cc74ULL,0x020bbea63b139b22ULL,0x514a08798e3404ddULL,
 0xef9519b3cd3a431bULL,0x302b0a6df25f1437ULL,0x4fe1356d6d51c245ULL, 
0xe485b576625e7ec7ULL,0xf44c42e9a637ed6bULL}}};
/*
1.57079632679489661923132169163975144209858469968755291048747229615390820314310449
*/
Qfloat qPi_Div_2[1] = {{
// exponent is 1 less, i.e. PI/2
0,EXPONE,{0xc90fdaa22168c234ULL,0xc4c6628b80dc1cd1ULL,
0x29024e088a67cc74ULL,0x020bbea63b139b22ULL,0x514a08798e3404ddULL,
0xef9519b3cd3a431bULL,0x302b0a6df25f1437ULL}}};

Qfloat qinv_pi[1]={
{0,0x0007ffff,{0xa2f9836e4e441529ULL,0xfc2757d1f534ddc0ULL,0xdb6295993c439041ULL,
0xfe5163abdebbc561ULL, 0xb7246e3a424dd2e0ULL,0x06492eea09d1921cULL,0xfe1deb1cb129a73fULL}}};
/*  Euler's constant
5.772156649015328606065120900824024310421593359399235988
     057672348848677267776646709369470632917467495146314472  
*/
Qfloat qeul[1] = {{
0,EXPONE-1,{0x93c467e37db0c7a4ULL,0xd1be3f810152cb56ULL,
0xa1cecc3af65cc019ULL,0x0c03df34709affbdULL,0x8e4b59fa03a9f0eeULL,
0xd0649ccb621057d1ULL,0x1056ae9132135a08ULL}}}; // next 2: e43b4673d74bafea 58deb878cc86d734

/* log(10)=2.3025850929940456840179914546843642076011014886287729760333279009\
67572609677352480235997205089598298341967784042286248633409525465083 */
Qfloat qlog10c[1] = {{
0,EXPONE+1,{0x935d8dddaaa8ac16ULL,0xea56d62b82d30a28ULL,
0xe28fecf9da5df90eULL,0x83c61e8201f02d72ULL,0x962f02d7b1a8105cULL,
0xcc70cbc02c5f0d68ULL,0x2c622418410be2dbULL }}}; //, 0xfb8f788402e516d6ULL,0x782cf8a28a8c911eULL }}}; 

Qfloat qinv_log10[1]={
{0,EXPONE-2,{0xde5bd8a937287195ULL,0x355baaafad33dc32ULL,
 0x3ee3460245c9a202ULL,0x3a3f2d44f78ea53cULL,0x75424efa1402f3f2ULL,
 0x92235592c6464a15ULL,0x18ce3bd9fd38dcbcULL,/* 0x6fa2b8d2c8cda7b3ULL,0x4356bd1948d06ffaULL*/}}};


// -exp(-1.0) Used in Lambert's W function.
// -0.367879441171442321595523770161460867445811131031767834507836801697461495744899803357147274345919643746627
Qfloat qmem1[1] = {{
-1,0x0007ffff,{0xbc5ab1b16779be35ULL,0x75bd8f0520a9f21bULL,
0xb5300b556ad8ee66ULL,0x604973a14a0fb5dbULL,0x62c8017e8e56842bULL,
0x7048b6cd3b21a4f4ULL,0xb5d4aaa12d1b8724ULL}}};

//exp(1)
// 2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742747
Qfloat qexp1[1] = {{
0,0x00080002,{0xadf85458a2bb4a9aULL,0xafdc5620273d3cf1ULL,
0xd8b9c583ce2d3695ULL,0xa9e13641146433fbULL,0xcc939dce249b3ef9ULL,
0x7d2fe363630c75d8ULL,0xf681b202aec4617bULL}}};

Qfloat invSqrt2pi[1] = {{
0,0x0007ffff,{0xcc42299ea1b28468ULL,0x7e59e2805d5c717fULL,
0xa7053e60a24e0db1ULL,0x808e369004c73a5dULL,0x3bf158aa88fa7b6eULL,
0x47eba75766e9d11bULL,0xd6284ef2e17c2a26ULL}}};

// pow(2,-NBITS)
// 1.3758210268297397763667897526170065161360224664295482687804518709187165186228035782046256351
Qfloat qepsilon[1] = {{ 0,0x0007fe41,{0x8000000000000000ULL,0,0,0,0,0,0}}};

Qfloat qten[1]={{0,EXPONE+3,{0xa000000000000000ULL,0,0,0,0,0,0}}};
Qfloat qfifteen[1]={{0,EXPONE+3,{0xf000000000000000ULL,0,0,0,0,0,0}}};
Qfloat qeight[1]={{0,EXPONE+3,{0x8000000000000000ULL,0,0,0,0,0,0}}};

/* 0.6931471805599453093744803655607000791860627941787242889404296875 */
Qfloat C1[1] = {{0,EXPONE-1,{0xb17217f7d1cf79abULL}}};
/*
4.27517558974764888894373401815309651802503219933936219696947156058633269964186875420014810205706857336855202357581305570326707516350E-20
*/
Qfloat C2[1]= {
{0x00000000,EXPONE-65,{0xc9e3b39803f2f6afULL,0x40f343267298b62dULL,
0x8a0d175b8baafa2bULL,0xe7b876206debac98ULL,0x559552fb4afa1b10ULL,
0xed2eae35c1382144ULL,0x27573b291169b825ULL}}}; 
/*

static Qfloat C1[1] = {{0,EXPONE-1,{0xb17217f700000000ULL}}};
static Qfloat C2[1] = {
{0,EXPONE-33,{0xd1cf79abc9e3b398ULL,0x03f2f6af40f34326ULL,
0x7298b62d8a0d175bULL,0x8baafa2be7b87620ULL,0x6debac98559552fbULL,
0x4afa1b10ed2eae35ULL,0xc138214427573b29ULL}}};
*/
