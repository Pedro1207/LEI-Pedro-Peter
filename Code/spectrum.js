var A = 3.48e-3

var windDir = np.array([1, 0])
var windSpeed = 20

var FW = 256
var L = 256

var g = 9.81
var fetchValue = 120000
var depth = 20
var swell = 1

var Hs = 10

function resetGlobalVariables(){

    A = 3.48e-3
    
    windDir = np.array([1, 0])
    windDir = windDir / norm(windDir)
    
    windSpeed = 20
    
    FW = 256
    
    L = 256
    
    fetchValue = 120000
    
    depth = 20
    
    swell = 1

    Hs = 10
}
    
    
    // from https://github.com/libretro/glsl-shaders/tree/master/crt/shaders/crt-royale/port-helpers
function gamma_1(s){
    var s_inv = 1/s
    var g = 1.12906830989
    var c0 = 0.8109119309638332633713423362694399653724431 
    var c1 = 0.4808354605142681877121661197951496120000040
    var e = 2.71828182845904523536028747135266249775724709
    var sph = s + 0.5
    var lanczos_sum = c0 + c1/(s + 1.0)
    var base = (sph + g)/e
    return (pow(base, sph) * lanczos_sum) * s_inv
}

// Stirling approximation
// https://rosettacode.org/wiki/Gamma_function#C
function gamma_2(x){
    return sqrt(2.0 * pi/x)*pow(x/e, x)
}


// Spouge's approximation
// https://rosettacode.org/wiki/Gamma_function#C
function gamma_sp(z){

    
    var a = 12; // number of steps - dictates precision
    var c = np.zeros(a);
    
    // this should be done only once
    var k1_factrl = 1.0; // (k - 1)!*(-1)^k with 0!==1
    c[0] = sqrt(2.0*pi);
    for(var i = 1; i < a; i++){
        c[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
        k1_factrl *= -k;
    }
    
    // computation based on table c
    accm = c[0];
    for(var i = 1; i < a; i++){
        accm += c[k] / ( z + k );
    }
    
    accm *= exp(-(z+a)) * pow(z+a, z+0.5); // Gamma(z+1) 
    return accm/z;
}



// Dispersion
DISP_DEEP = 0
DISP_SHALLOW = 1
DISP_CAPILLARY = 2

function dispersion(dispersionMode, k){

    
    var tension = 0.074 // N/m
    var density = 1000.0 // Kg/m3
    
    if (dispersionMode == DISP_DEEP){
        return sqrt(k * g);
    }
    else if (dispersionMode ==  DISP_SHALLOW){
        return sqrt(k * g * tanh(k*depth));
    }
    else if (dispersionMode ==  DISP_CAPILLARY){
        return sqrt((g * k + pow(k, 3) * tension / density) * tanh(k*depth));
    }
    else{
        return 0;
    }
    
}
    
function dispersionDeriv(dispersionMode, k){

    var tension = 0.074; // N/m
    var density = 1000.0; // Kg/m3
    var w = dispersion(dispersionMode, k);
        
    if (dispersionMode == DISP_DEEP){
        return 0.5 * g / w;
    }
    else if (dispersionMode ==  DISP_SHALLOW){
        var dk = depth * k;
        var th = tanh(dk);
        return 0.5 * sqrt(g/(k*th)) * (th + dk*(1 - th*th)) 
    }
    else if (dispersionMode ==  DISP_CAPILLARY){
        var dk = depth * k;
        var th = tanh(dk);
        var b = tension / density;
        return 0.5 * ((9.81 + 3 * b*k*k)*th + dk * (b*k*k + 9.81) * pow(1.0 / cosh(dk), 2)) / w;
    } 
    else{
        return 0;
    } 
}



// Bretschneider Spectrum defined as a function of significant wave height
// https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/lec6_wavespectra.pdf


function Bretschneider(w){
    var wm = 0.4 * sqrt(g/Hs)
    return 1.25 / 4 * pow (wm,4) / pow (w,5) * Hs * Hs * exp(-1.25 * pow(wm/w, 4))
}


// Pierson-Moskowitz defined as a function of significant wave height
// https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/lec6_wavespectra.pdf

function PiersonMoskowitz(w){    
    var wp = 0.879 * g / windSpeed
    var expTerm = -1.25 * pow(wp / w, 4)
    var result = 8.1e-3 * pow(g, 2) * pow(w, -5) * exp(expTerm)
    return result
}


function PiersonMoskowitzPeakW(){
    return 0.879 * g / windSpeed;
}


function PiersonMoskowitzHs(w){    
    var expTerm = -0.032 * pow(g / Hs , 2) * pow(w, -4)
    var result = 8.1e-3 * pow(g, 2) * pow(w, -5) * exp(expTerm)
    return result
}


function PiersonMoskowitzHsPeakW(){
    return 0.4 * sqrt(g/Hs)
}


//JONSWAP

var gamma_js = 3.3
var sigma_a = 0.07
var sigma_b = 0.09
var alpha = 0.076 * pow(g * fetchValue / pow(windSpeed, 2), -0.22)



function JONSWAP(w){

    var dimensionlessFetch = g * fetchValue / pow(windSpeed, 2)    
    var wp = 22 * (g / windSpeed) * pow(dimensionlessFetch, -0.33)
    var expTerm = -1.25 * pow(wp / w, 4)
    
    if (w <= wp){
        sigma = sigma_a
    }
    else{
        sigma = sigma_b
    }
    
    var gammaExp = -0.5 * pow((w - wp) / (sigma*wp), 2)    
    
    var result = alpha * g * g * pow(w, -5) * exp(expTerm) * pow(gamma_js, exp(gammaExp))
    return result
}      


function JONSWAPPeakW(){
    return 22 * (g / windSpeed) * pow(windSpeed*windSpeed / (g * fetchValue), 0.33)
}


// JONSWAP with Hasselmann et.al. suggested average parameter values
function JONSWAP_Hasselmann(w){

    gamma_js = 3.3
    sigma_a = 0.07
    sigma_b = 0.09
    var dimensionlessFetch = g * fetchValue / pow(windSpeed, 2)   
    var alpha = 0.076 * pow(dimensionlessFetch, -0.22)
    return JONSWAP(w)
}
    
    
// alpha and gamma settings suggested by Mitsuyasu (1985)- see Ochi, 1998
function JONSWAP_Mitsuyasu(w){
    var dimensionlessFetch = g * fetchValue / pow(windSpeed, 2)
    var gamma_js = 7.0 * pow(dimensionlessFetch, -0.143)
    var alpha = 8.17e-2 * pow(dimensionlessFetch, -0.283)
    return JONSWAP(w)
}


//TMA

function TMA_JONSWAP(w){

    
    var wp = getPeakW(SPECTRUM_JONSWAP)
    // this is a simplification since it is assumed that we are using the deep water dispersion relationship
    var km = wp*wp/g
    var k = windSpeed * windSpeed * km / g
    alpha = 0.0078 * pow(k, 0.49)
    gamma_js = 2.47 * pow(k, 0.39)
    
    return JONSWAP(w);
}
    

function TMA(w){

    var result = TMA_JONSWAP(w)
    
    var wh = np.clip(w * sqrt(depth / 9.81), 0.0, 2.0);
    if (wh <= 1){
        var theta = 0.5 * wh*wh
    }
    else{
        var theta = 1 - 0.5 * pow(2 - wh, 2);
        np.clip(theta, 0, 1)
    }
    
    return result * theta;
}


function TMA_Hasselmann(w){

    
    var result = JONSWAP_Hasselmann(w);
    
    var wh = np.clip(w * sqrt(depth / 9.81), 0.0, 2.0);
    if (wh <= 1){
        var theta = 0.5 * wh*wh
    }
    else{
        var theta = 1 - 0.5 * pow(2 - wh, 2);
        np.clip(theta, 0, 1)

    }
    
    return result * theta;
}


//Donelan

function DonelanJONSWAP_IWA(w,iwa){
    var Omega = iwa

    var wp = Omega * g / windSpeed
    var expTerm = -pow(wp/w, 4)
    
    var sigma = 0.08 * (1 + 4 / pow(Omega,3))
        
    var gammaExp = -0.5 * pow(w - wp, 2) / (pow(sigma * wp, 2))
    
    // beta replaces alpha
    var beta = 0.006 * pow(Omega, 0.55);

    if (Omega >= 1){
        var gammat = 1.7 + 6 * log(Omega);
    }
    else{
        var gammat = 1.7
    }

    var result = beta * g * g * pow(wp, -1) * pow(w, -4) * exp(expTerm) * pow(gammat, exp(gammaExp));
    return result;    

}


function DonelanJONSWAP(w){

    
    var dimensionlessFetch = g * fetchValue / pow(windSpeed, 2)    
    
    var Omega = 11.6 * pow(dimensionlessFetch,-0.23)
    
    var wp = Omega * g / windSpeed
    var expTerm = -pow(wp/w, 4)
    
    var sigma = 0.08 * (1 + 4 / pow(Omega,3))
    
    var gammaExp = -0.5 * pow(w - wp, 2) / (pow(sigma * wp, 2))
    
    // beta replaces alpha
    var beta = 0.006 * pow(Omega, 0.55);
    
    if (Omega >= 1){
        var gammat = 1.7 + 6 * log(Omega);
    }
    else{
        var gammat = 1.7
    }
    
    var result = beta * g * g * pow(wp, -1) * pow(w, -4) * exp(expTerm) * pow(gammat, exp(gammaExp));
    return result;    
    
}

function DonelanJONSWAPPeakW(){
    
    var dimensionlessFetch = g * fetchValue / pow(windSpeed, 2)    
    var Omega = 11.6 * pow(dimensionlessFetch,-0.23)
    return Omega * g / windSpeed
}



//tessendorf
function Tessendorf(kvec){

    var k = np.sqrt(kvec.dot(kvec))
    if (k == 0.0){
        return 0.0
    }
    
    var l = windSpeed * windSpeed / g;
    var khat = kvec/k;
    var wd = windDir / np.sqrt(windDir.dot(windDir));
    var dot_k_w = khat.dot(wd);
    var result = A / pow(k, 4) * exp(-1 / (k*l*k*l)) * pow(dot_k_w, 2);
    result *= exp(-k*k*0.1*0.1);
    return result;
}

// From A unified directional spectrumfor long and short wind-driven waves
// T. Elfouhaily,B. Chapron,and K. Katsaros
// https://archimer.ifremer.fr/doc/00091/20226/

function omega(k){
    var km = 370
    return (sqrt(g * k * (1 + pow(k/km,2)))) // eq 24.
}


function Unified(directional, kvec){

    var k = np.sqrt(kvec.dot(kvec))
    if (k == 0.0){
        return 0.0
    }
    
    
    var X = g * fetch / pow(windSpeed, 2) // after eq. 4
    
    var X0 = 2.2e4
    var Omega = 0.84 * pow(tanh(pow(X/X0, 0.4)), -0.75) //eq 37
    
    var kp = g * pow(Omega/windSpeed,2) // below eq. 3
    
    var c = omega(k) / k // wave phase speed
    var cp = omega(kp) / kp // phase speed for peak frequency
    
    var alphap = 6e-3 * sqrt(Omega) //34
    
    var Lpm = exp( - 1.25 * pow(kp/k, 2)) // eq. 2
    
    var sigma = 0.08 * (1 + 4 * pow(Omega, -3))  // below eq. 3
    var Gamma = exp(- 0.5 * pow(sqrt(k/kp) - 1, 2) / pow(sigma,2)) 
    if (Omega < 1){
        var gammat = 1.7
    }
    else{
        var gammat = 1.7 + 6 * log(Omega)
    }
    
    var Jp = pow(gammat, Gamma) // eq. 3
    
    var Fp = Lpm * Jp * exp( - Omega / sqrt(10) * (sqrt(k/kp) - 1)) // eq. 32
    
    var Bl = 0.5 * alphap * (cp/c) * Fp // eq. 31
    
    var cm = 0.23 // eq. 59
    
    var z0 = 3.7e-5 * pow(windSpeed,2) / g * pow(windSpeed/cp, 0.9) // eq. 66
    // https://en.wikipedia.org/wiki/Von_K%C3%A1rm%C3%A1n_constant
    var K = 0.41
    varustar = K * windSpeed / log(10.0/z0) // eq. 60
    
    if (ustar < cm){
        alpham = 1 + log(ustar/cm) // eq 44 
    } 
    else{
        alpham = 1 + 3 * log(ustar/cm)     
        alpham = 0.01 * alpham 
    }
    
    var km = 370
    var Fm = exp(-0.25 * pow(k/km - 1, 2))
    
    var Bh = 0.5 * alpham * (cm/c) * Fm  * Lpm // eq. 40 (fixed?)
    
    if (!directional){
        return pow(k,-3)  * (Bl + Bh) 
    }
    
    var a0 = log(2)/4
    
    var ap = 4
    var am = 0.13 * ustar / cm
    var delta = tanh(a0 + ap * pow(c/cp, 2.5) + am * (pow(cm/c, 2.5))) // eq. 57
    
    var cos_twice_angle = 2 * pow(np.dot(windDir,kvec/norm(kvec)),2) -1
    return 1 / (2 * pi) * pow(k,-4) * (Bl + Bh) * (1 + delta * cos_twice_angle) //eq. 67
}


//ochi

var l1 = 3.0
var wm1 = 0.7 * exp(-0.046 * Hs);
var hs1 = 0.84 * Hs;
var l2 = 1.54 * exp(-0.062 *  Hs);
var wm2 = 1.15 * exp(-0.039 * Hs);
var hs2 = 0.54 * Hs;


function OchiSingle(w, l, wm, hs){
    var aux = (4 * l + 1)/4
    
    var s = 0.25 * pow(aux * pow(wm, 4), l)/gamma(l) * pow(hs,2)/pow(w, 4*l + 1) * exp(-aux * pow(wm/w, 4))
    return s
}



function Ochi(w, lambda1, wm1, hs1, lambda2, wm2, hs2){
    var s = OchiSingle(w, lambda1, wm1,hs1) + OchiSingle(w, lambda2, wm2, hs2)
    return s
}
    


function OchiHs(w, hs){
    return Ochi(w, 3, 0.7 * exp(-0.046 * hs), 0.84 * hs, 1.54 * exp(-0.062 * hs), 1.15 * exp(-0.039 * hs), 0.54 * hs)
}



// Pierson et. al, "Practical methods for observing and forecasting ocean waves by means of wave spectra and statistics", 1955

function cosPower(waveDir){
    var cos_angle = np.dot(windDir,waveDir/norm(waveDir))
    if (cos_angle > 0){
        return (2/pi  * pow(cos_angle,2))
    }
    else{
        return 0;
    }
}


// Mitsuyasu 1975
function mitsuyasu(w, waveDir, wPeak){

    var cos_angle = np.dot(windDir/norm(windDir),waveDir/norm(waveDir));
    var half_angle_cos = sqrt((1+ cos_angle)/2);
    var sp = 11.5 * pow((wPeak * windSpeed /g),-2.5);
    if (w <= wPeak){
        var s = sp * pow(w/wPeak, 5);
    }
    else{
        var s = sp * pow(w/wPeak, -2.5)  ;
    }
    
    var qs = pow(2, 2*s-1) / pi * pow(gamma(s+1), 2)/ gamma(2*s +1);
    
    return qs * pow(half_angle_cos, 2*s);
}


// Hasselmann 1980
function hasselmann(w, waveDir, wPeak){
   
    var cos_angle = np.dot(windDir,waveDir/norm(waveDir))
    var half_angle_cos = sqrt((1+ cos_angle)/2)
    var sp = 11.5 * pow((wPeak * windSpeed /g),-2.5)
    if (w <= wPeak){
        var s = 6.97 * pow(w/wPeak, 4.06)
    }
    else{
        var s = 9.77 * pow(w/wPeak, -2.33-1.45*((windSpeed * wPeak / g) - 1.17))    
        var qs = pow(2, 2*s-1) / pi * pow(gamma(s+1), 2)/ gamma(2*s +1)
        return qs * pow(abs(half_angle_cos), 2*s)
    }
}



// Donelan Hamilton and Hui 1985
function donelan(w, waveDir, wPeak){

    var cos_angle = np.dot(windDir,waveDir/norm(waveDir))
    var wave_angle = acos(cos_angle)
    var wcoeff = w/wPeak
    var epsilon = -0.4 + 0.8393 * exp( -0.567 * log(pow(wcoeff,2)));
    if (wcoeff < 0.95){
        var beta = 2.61 * pow(wcoeff, 1.3);
        
    }
    else if (wcoeff < 1.6){
        var beta = 2.28 * pow(wcoeff, -1.3);
    }
    else{
        var beta = pow(10,epsilon);
    }
    var res = beta / (2 * tanh(beta * pi)) * pow(1/cosh(beta * wave_angle),2);
    return res;
    
}



//numerical integration to compute normalization factor

function integrateCosine(power){

    var s = 0 
    var theta = 0;
    for(var i = 0; i < 30; i++){
        theta = -pi + (i+1)*2*pi/32
        s += pow(abs(cos(theta/2)), power);
    }  
    
    theta = -pi/2  
    var s0 = pow(abs(cos(theta/2)), power)
    
    theta = pi/2  
    var sn = pow(abs(cos(theta/2)), power)
    
    var res = (s + 0.5 * (s0 + sn))  * (2 * pi/32)   
    return res 
}
    

var dirSpreadUseNumIntegration = true

// Mitsuyasu + Horvath, 2015 
function mitsuyasuHorvath(w, waveDir, wPeak){

    var cos_angle = np.dot(windDir/norm(windDir),waveDir/norm(waveDir))
    var half_angle_cos = sqrt((1+ cos_angle)/2)
    var sp = 11.5 * pow((wPeak * windSpeed /g),-2.5)
    if (w <= wPeak){
        var s = sp * pow(w/wPeak, 5)
    }
    else{
        var s = sp * pow(w/wPeak, -2.5)  
    }
    
    s += 16 * tanh(wPeak/w) * swell   
    
    if (dirSpreadUseNumIntegration){
        return 1 / integrateCosine(2*s) * pow(abs(half_angle_cos), 2*s)
    }
    
    else{
        var qs = pow(2, 2*s-1) / pi * pow(gamma(s+1), 2)/ gamma(2*s +1)
        return qs * pow(abs(half_angle_cos), 2*s)
    } 
}
    
    
//Hasselmann + Horvath, 2015

function hasselmannHorvath(w, waveDir, wPeak){
    var cos_angle = np.dot(windDir/norm(windDir),waveDir/norm(waveDir))
    var half_angle_cos = sqrt((1+ cos_angle)/2)
    var sp = 11.5 * pow((wPeak * windSpeed /g),-2.5)
    if (w <= wPeak){
        var s = 6.97 * pow(w/wPeak, 4.06)
    }
    else{
        var s = 9.77 * pow(w/wPeak, -2.33-1.45*((windSpeed * wPeak / g) - 1.17))  
    }

    s += 16 * tanh(wPeak/w) * swell
    if (dirSpreadUseNumIntegration){
        return 1 / integrateCosine(2*s) * pow(abs(half_angle_cos), 2*s)
    }
    else{
        qs = pow(2, 2*s-1) / pi * pow(gamma(s+1), 2)/ gamma(2*s +1)
        return qs * pow(abs(half_angle_cos), 2*s)
    }

}



// Donelan-Banner + Horvath, 2015

// aux fuction to compute Horvath
function horvath(w, waveDir, wPeak){
    var cos_angle = np.dot(windDir/norm(windDir),waveDir/norm(waveDir))
    var half_angle_cos = sqrt((1+ cos_angle)/2)
    var s = 16 * tanh(wPeak/w) * swell
    var qs = pow(2, 2*s-1) / pi * pow(gamma(s+1), 2)/ gamma(2*s +1)
    return qs * pow(abs(half_angle_cos), 2*s)
}
    
    
//numerical integration to compute normalization factor
function integrateDonelanHorvath(w, wPeak){
    var s = 0
    var theta, wDir;

    for(var i = 0; i < 30; i++){
        theta = -pi/2 + (i+1)*pi/32
        wDir = np.array((cos(theta), sin(theta)))
        s += donelan(w, wDir, wPeak) * horvath(w, wDir, wPeak);
    }
    
    var wave0 =  np.array((cos(-pi/2), sin(-pi/2)))  
    var f0 = donelan(w, wave0, wPeak) * horvath(w, wave0, wPeak)
    
    var waven =  np.array((cos(pi/2), sin(pi/2)))  
    var fn = donelan(w, waven, wPeak) * horvath(w, waven, wPeak)
    
    var res = (s + 0.5 * (f0 + fn))  * (2 * pi/32)   
    return res 
}
        
    
    
function donelanHorvath(w, waveDir, wPeak){
    var res = donelan(w, waveDir, wPeak) * horvath(w, waveDir, wPeak)
    var normFactor = integrateDonelanHorvath(w, wPeak)
    
    var res = res / normFactor;
    return res
}    
        
                   