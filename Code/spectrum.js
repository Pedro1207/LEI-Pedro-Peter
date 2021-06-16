export var A = 3.48e-3

export var windDir = [1, 0]
export var windSpeed = 20

export var FW = 256
export var L = 256

export var g = 9.81
export var fetchValue = 120000
export var depth = 20
export var swell = 1

export var Hs = 10


export var dot = (a, b) => a.map((x, i) => a[i] * b[i]).reduce((m, n) => m + n);

export function clamp(x, a, b){
    return Math.max( a, Math.min(x, b))
} 


export function updateVariables(windSpeed_, windDir_, L_, fetchValue_, depth_, swell_, Hs_, gamma_js_, sigma_a_, sigma_b_, l1_, l2_, wm1_, wm2_, hs1_, hs2_){
    windSpeed = windSpeed_;
    windDir = windDir_;
    L = L_;
    fetchValue = fetchValue_;
    depth = depth_;
    swell = swell_;
    Hs = Hs_;
    gamma_js = gamma_js_;
    sigma_a = sigma_a_;
    sigma_b = sigma_b_;
    l1 = l1_;
    l2 = l2_;
    wm1 = wm1_;
    wm2 = wm2_;
    hs1 = hs1_;
    hs2 = hs2_;
}


export function resetGlobalVariables(){

    A = 3.48e-3
    
    windDir = [1, 0]
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
export function gamma_1(s){
    var s_inv = 1/s
    var g = 1.12906830989
    var c0 = 0.8109119309638332633713423362694399653724431 
    var c1 = 0.4808354605142681877121661197951496120000040
    var e = 2.71828182845904523536028747135266249775724709
    var sph = s + 0.5
    var lanczos_sum = c0 + c1/(s + 1.0)
    var base = (sph + g)/e
    return (Math.pow(base, sph) * lanczos_sum) * s_inv
}

// Stirling approximation
// https://rosettacode.org/wiki/Gamma_function#C
export function gamma_2(x){
    return Math.sqrt(2.0 * pi/x)*Math.pow(x/e, x)
}


// Spouge's approximation
// https://rosettacode.org/wiki/Gamma_function#C
export function gamma_sp(z){

    
    var a = 12; // number of steps - dictates precision
    var c = new Array(a).fill(0);
    
    // this should be done only once
    var k1_factrl = 1.0; // (k - 1)!*(-1)^k with 0!==1
    c[0] = Math.sqrt(2.0*pi);
    for(var i = 1; i < a; i++){
        c[k] = Math.exp(a-k) * Math.pow(a-k, k-0.5) / k1_factrl;
        k1_factrl *= -k;
    }
    
    // computation based on table c
    accm = c[0];
    for(var i = 1; i < a; i++){
        accm += c[k] / ( z + k );
    }
    
    accm *= Math.exp(-(z+a)) * Math.pow(z+a, z+0.5); // Gamma(z+1) 
    return accm/z;
}



// Dispersion
export const DISP_DEEP = 0
export const DISP_SHALLOW = 1
export const DISP_CAPILLARY = 2

export function dispersion(dispersionMode, k){

    
    var tension = 0.074 // N/m
    var density = 1000.0 // Kg/m3
    
    if (dispersionMode == DISP_DEEP){
        return Math.sqrt(k * g);
    }
    else if (dispersionMode ==  DISP_SHALLOW){
        return Math.sqrt(k * g * tanh(k*depth));
    }
    else if (dispersionMode ==  DISP_CAPILLARY){
        return Math.sqrt((g * k + Math.pow(k, 3) * tension / density) * tanh(k*depth));
    }
    else{
        return 0;
    }
    
}
    
export function dispersionDeriv(dispersionMode, k){

    var tension = 0.074; // N/m
    var density = 1000.0; // Kg/m3
    var w = dispersion(dispersionMode, k);
        
    if (dispersionMode == DISP_DEEP){
        return 0.5 * g / w;
    }
    else if (dispersionMode ==  DISP_SHALLOW){
        var dk = depth * k;
        var th = tanh(dk);
        return 0.5 * Math.sqrt(g/(k*th)) * (th + dk*(1 - th*th)) 
    }
    else if (dispersionMode ==  DISP_CAPILLARY){
        var dk = depth * k;
        var th = tanh(dk);
        var b = tension / density;
        return 0.5 * ((9.81 + 3 * b*k*k)*th + dk * (b*k*k + 9.81) * Math.pow(1.0 / cosh(dk), 2)) / w;
    } 
    else{
        return 0;
    } 
}



// Bretschneider Spectrum defined as a function of significant wave height
// https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/lec6_wavespectra.pdf


export function Bretschneider(w){
    var wm = 0.4 * Math.sqrt(g/Hs)
    return 1.25 / 4 * Math.pow (wm,4) / Math.pow (w,5) * Hs * Hs * Math.exp(-1.25 * Math.pow(wm/w, 4))
}


// Pierson-Moskowitz defined as a function of significant wave height
// https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/lec6_wavespectra.pdf

export function PiersonMoskowitz(w){    
    var wp = 0.879 * g / windSpeed
    var expTerm = -1.25 * Math.pow(wp / w, 4)
    var result = 8.1e-3 * Math.pow(g, 2) * Math.pow(w, -5) * Math.exp(expTerm)
    return result
}


export function PiersonMoskowitzPeakW(){
    return 0.879 * g / windSpeed;
}


export function PiersonMoskowitzHs(w){    
    var expTerm = -0.032 * Math.pow(g / Hs , 2) * Math.pow(w, -4)
    var result = 8.1e-3 * Math.pow(g, 2) * Math.pow(w, -5) * Math.exp(expTerm)
    return result
}


export function PiersonMoskowitzHsPeakW(){
    return 0.4 * Math.sqrt(g/Hs)
}


//JONSWAP

export var gamma_js = 3.3
export var sigma_a = 0.07
export var sigma_b = 0.09
export var alpha = 0.076 * Math.pow(g * fetchValue / Math.pow(windSpeed, 2), -0.22)



export function JONSWAP(w){

    var dimensionlessFetch = g * fetchValue / Math.pow(windSpeed, 2)    
    var wp = 22 * (g / windSpeed) * Math.pow(dimensionlessFetch, -0.33)
    var expTerm = -1.25 * Math.pow(wp / w, 4)
    
    if (w <= wp){
        var sigma = sigma_a
    }
    else{
        var sigma = sigma_b
    }
    
    var gammaExp = -0.5 * Math.pow((w - wp) / (sigma*wp), 2)    
    
    var result = alpha * g * g * Math.pow(w, -5) * Math.exp(expTerm) * Math.pow(gamma_js, Math.exp(gammaExp))
    return result
}      


export function JONSWAPPeakW(){
    return 22 * (g / windSpeed) * Math.pow(windSpeed*windSpeed / (g * fetchValue), 0.33)
}


// JONSWAP with Hasselmann et.al. suggested average parameter values
export function JONSWAP_Hasselmann(w){

    gamma_js = 3.3
    sigma_a = 0.07
    sigma_b = 0.09
    var dimensionlessFetch = g * fetchValue / Math.pow(windSpeed, 2)   
    var alpha = 0.076 * Math.pow(dimensionlessFetch, -0.22)
    return JONSWAP(w)
}
    
    
// alpha and gamma settings suggested by Mitsuyasu (1985)- see Ochi, 1998
export function JONSWAP_Mitsuyasu(w){
    var dimensionlessFetch = g * fetchValue / Math.pow(windSpeed, 2)
    var gamma_js = 7.0 * Math.pow(dimensionlessFetch, -0.143)
    var alpha = 8.17e-2 * Math.pow(dimensionlessFetch, -0.283)
    return JONSWAP(w)
}


//TMA

export function TMA_JONSWAP(w){

    
    var wp = getPeakW(SPECTRUM_JONSWAP)
    // this is a simplification since it is assumed that we are using the deep water dispersion relationship
    var km = wp*wp/g
    var k = windSpeed * windSpeed * km / g
    alpha = 0.0078 * Math.pow(k, 0.49)
    gamma_js = 2.47 * Math.pow(k, 0.39)
    
    return JONSWAP(w);
}
    

export function TMA(w){

    var result = TMA_JONSWAP(w)
    
    var wh = clamp(w * Math.sqrt(depth / 9.81), 0.0, 2.0);
    if (wh <= 1){
        var theta = 0.5 * wh*wh
    }
    else{
        var theta = 1 - 0.5 * Math.pow(2 - wh, 2);
        clamp(theta, 0, 1)
    }
    
    return result * theta;
}


export function TMA_Hasselmann(w){

    
    var result = JONSWAP_Hasselmann(w);
    
    var wh = clamp(w * Math.sqrt(depth / 9.81), 0.0, 2.0);
    if (wh <= 1){
        var theta = 0.5 * wh*wh
    }
    else{
        var theta = 1 - 0.5 * Math.pow(2 - wh, 2);
        clamp(theta, 0, 1)

    }
    
    return result * theta;
}


//Donelan

export function DonelanJONSWAP_IWA(w,iwa){
    var Omega = iwa

    var wp = Omega * g / windSpeed
    var expTerm = -Math.pow(wp/w, 4)
    
    var sigma = 0.08 * (1 + 4 / Math.pow(Omega,3))
        
    var gammaExp = -0.5 * Math.pow(w - wp, 2) / (Math.pow(sigma * wp, 2))
    
    // beta replaces alpha
    var beta = 0.006 * Math.pow(Omega, 0.55);

    if (Omega >= 1){
        var gammat = 1.7 + 6 * log(Omega);
    }
    else{
        var gammat = 1.7
    }

    var result = beta * g * g * Math.pow(wp, -1) * Math.pow(w, -4) * Math.exp(expTerm) * Math.pow(gammat, Math.exp(gammaExp));
    return result;    

}


export function DonelanJONSWAP(w){

    
    var dimensionlessFetch = g * fetchValue / Math.pow(windSpeed, 2)    
    
    var Omega = 11.6 * Math.pow(dimensionlessFetch,-0.23)
    
    var wp = Omega * g / windSpeed
    var expTerm = -Math.pow(wp/w, 4)
    
    var sigma = 0.08 * (1 + 4 / Math.pow(Omega,3))
    
    var gammaExp = -0.5 * Math.pow(w - wp, 2) / (Math.pow(sigma * wp, 2))
    
    // beta replaces alpha
    var beta = 0.006 * Math.pow(Omega, 0.55);
    
    if (Omega >= 1){
        var gammat = 1.7 + 6 * log(Omega);
    }
    else{
        var gammat = 1.7
    }
    
    var result = beta * g * g * Math.pow(wp, -1) * Math.pow(w, -4) * Math.exp(expTerm) * Math.pow(gammat, Math.exp(gammaExp));
    return result;    
    
}

export function DonelanJONSWAPPeakW(){
    
    var dimensionlessFetch = g * fetchValue / Math.pow(windSpeed, 2)    
    var Omega = 11.6 * Math.pow(dimensionlessFetch,-0.23)
    return Omega * g / windSpeed
}



//tessendorf
export function Tessendorf(kvec){

    var k = (kvec.dot(kvec)).map(Math.sqrt)
    if (k == 0.0){
        return 0.0
    }
    
    var l = windSpeed * windSpeed / g;
    var khat = kvec/k;
    var wd = windDir / (windDir.dot(windDir)).map(Math.sqrt);
    var dot_k_w = khat.dot(wd);
    var result = A / Math.pow(k, 4) * Math.exp(-1 / (k*l*k*l)) * Math.pow(dot_k_w, 2);
    result *= Math.exp(-k*k*0.1*0.1);
    return result;
}

// From A unified directional spectrumfor long and short wind-driven waves
// T. Elfouhaily,B. Chapron,and K. Katsaros
// https://archimer.ifremer.fr/doc/00091/20226/

export function omega(k){
    var km = 370
    return (Math.sqrt(g * k * (1 + Math.pow(k/km,2)))) // eq 24.
}


export function Unified(directional, kvec){

    var k = (kvec.dot(kvec)).map(Math.sqrt);
    if (k == 0.0){
        return 0.0
    }
    
    
    var X = g * fetch / Math.pow(windSpeed, 2) // after eq. 4
    
    var X0 = 2.2e4
    var Omega = 0.84 * Math.pow(tanh(Math.pow(X/X0, 0.4)), -0.75) //eq 37
    
    var kp = g * Math.pow(Omega/windSpeed,2) // below eq. 3
    
    var c = omega(k) / k // wave phase speed
    var cp = omega(kp) / kp // phase speed for peak frequency
    
    var alphap = 6e-3 * Math.sqrt(Omega) //34
    
    var Lpm = Math.exp( - 1.25 * Math.pow(kp/k, 2)) // eq. 2
    
    var sigma = 0.08 * (1 + 4 * Math.pow(Omega, -3))  // below eq. 3
    var Gamma = Math.exp(- 0.5 * Math.pow(Math.sqrt(k/kp) - 1, 2) / Math.pow(sigma,2)) 
    if (Omega < 1){
        var gammat = 1.7
    }
    else{
        var gammat = 1.7 + 6 * log(Omega)
    }
    
    var Jp = Math.pow(gammat, Gamma) // eq. 3
    
    var Fp = Lpm * Jp * Math.exp( - Omega / Math.sqrt(10) * (Math.sqrt(k/kp) - 1)) // eq. 32
    
    var Bl = 0.5 * alphap * (cp/c) * Fp // eq. 31
    
    var cm = 0.23 // eq. 59
    
    var z0 = 3.7e-5 * Math.pow(windSpeed,2) / g * Math.pow(windSpeed/cp, 0.9) // eq. 66
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
    var Fm = Math.exp(-0.25 * Math.pow(k/km - 1, 2))
    
    var Bh = 0.5 * alpham * (cm/c) * Fm  * Lpm // eq. 40 (fixed?)
    
    if (!directional){
        return Math.pow(k,-3)  * (Bl + Bh) 
    }
    
    var a0 = log(2)/4
    
    var ap = 4
    var am = 0.13 * ustar / cm
    var delta = tanh(a0 + ap * Math.pow(c/cp, 2.5) + am * (Math.pow(cm/c, 2.5))) // eq. 57
    
    var cos_twice_angle = 2 * Math.pow(dot(windDir,kvec/norm(kvec)),2) -1
    return 1 / (2 * pi) * Math.pow(k,-4) * (Bl + Bh) * (1 + delta * cos_twice_angle) //eq. 67
}


//ochi

export var l1 = 3.0
export var wm1 = 0.7 * Math.exp(-0.046 * Hs);
export var hs1 = 0.84 * Hs;
export var l2 = 1.54 * Math.exp(-0.062 *  Hs);
export var wm2 = 1.15 * Math.exp(-0.039 * Hs);
export var hs2 = 0.54 * Hs;


export function OchiSingle(w, l, wm, hs){
    var aux = (4 * l + 1)/4
    
    var s = 0.25 * Math.pow(aux * Math.pow(wm, 4), l)/gamma(l) * Math.pow(hs,2)/Math.pow(w, 4*l + 1) * Math.exp(-aux * Math.pow(wm/w, 4))
    return s
}



export function Ochi(w, lambda1, wm1, hs1, lambda2, wm2, hs2){
    var s = OchiSingle(w, lambda1, wm1,hs1) + OchiSingle(w, lambda2, wm2, hs2)
    return s
}
    


export function OchiHs(w, hs){
    return Ochi(w, 3, 0.7 * Math.exp(-0.046 * hs), 0.84 * hs, 1.54 * Math.exp(-0.062 * hs), 1.15 * Math.exp(-0.039 * hs), 0.54 * hs)
}



// Pierson et. al, "Practical methods for observing and forecasting ocean waves by means of wave spectra and statistics", 1955

export function cosPower(waveDir){
    var cos_angle = dot(windDir,waveDir/norm(waveDir))
    if (cos_angle > 0){
        return (2/pi  * Math.pow(cos_angle,2))
    }
    else{
        return 0;
    }
}


// Mitsuyasu 1975
export function mitsuyasu(w, waveDir, wPeak){

    var cos_angle = dot(windDir/norm(windDir),waveDir/norm(waveDir));
    var half_angle_cos = Math.sqrt((1+ cos_angle)/2);
    var sp = 11.5 * Math.pow((wPeak * windSpeed /g),-2.5);
    if (w <= wPeak){
        var s = sp * Math.pow(w/wPeak, 5);
    }
    else{
        var s = sp * Math.pow(w/wPeak, -2.5)  ;
    }
    
    var qs = Math.pow(2, 2*s-1) / pi * Math.pow(gamma(s+1), 2)/ gamma(2*s +1);
    
    return qs * Math.pow(half_angle_cos, 2*s);
}


// Hasselmann 1980
export function hasselmann(w, waveDir, wPeak){
   
    var cos_angle = dot(windDir,waveDir/norm(waveDir))
    var half_angle_cos = Math.sqrt((1+ cos_angle)/2)
    var sp = 11.5 * Math.pow((wPeak * windSpeed /g),-2.5)
    if (w <= wPeak){
        var s = 6.97 * Math.pow(w/wPeak, 4.06)
    }
    else{
        var s = 9.77 * Math.pow(w/wPeak, -2.33-1.45*((windSpeed * wPeak / g) - 1.17))    
        var qs = Math.pow(2, 2*s-1) / pi * Math.pow(gamma(s+1), 2)/ gamma(2*s +1)
        return qs * Math.pow(abs(half_angle_cos), 2*s)
    }
}



// Donelan Hamilton and Hui 1985
export function donelan(w, waveDir, wPeak){

    var cos_angle = dot(windDir,waveDir/norm(waveDir))
    var wave_angle = acos(cos_angle)
    var wcoeff = w/wPeak
    var epsilon = -0.4 + 0.8393 * Math.exp( -0.567 * log(Math.pow(wcoeff,2)));
    if (wcoeff < 0.95){
        var beta = 2.61 * Math.pow(wcoeff, 1.3);
        
    }
    else if (wcoeff < 1.6){
        var beta = 2.28 * Math.pow(wcoeff, -1.3);
    }
    else{
        var beta = Math.pow(10,epsilon);
    }
    var res = beta / (2 * tanh(beta * pi)) * Math.pow(1/cosh(beta * wave_angle),2);
    return res;
    
}



//numerical integration to compute normalization factor

export function integrateCosine(power){

    var s = 0 
    var theta = 0;
    for(var i = 0; i < 30; i++){
        theta = -pi + (i+1)*2*pi/32
        s += Math.pow(abs(cos(theta/2)), power);
    }  
    
    theta = -pi/2  
    var s0 = Math.pow(abs(cos(theta/2)), power)
    
    theta = pi/2  
    var sn = Math.pow(abs(cos(theta/2)), power)
    
    var res = (s + 0.5 * (s0 + sn))  * (2 * pi/32)   
    return res 
}
    

export var dirSpreadUseNumIntegration = true

// Mitsuyasu + Horvath, 2015 
export function mitsuyasuHorvath(w, waveDir, wPeak){

    var cos_angle = dot(windDir/norm(windDir),waveDir/norm(waveDir))
    var half_angle_cos = Math.sqrt((1+ cos_angle)/2)
    var sp = 11.5 * Math.pow((wPeak * windSpeed /g),-2.5)
    if (w <= wPeak){
        var s = sp * Math.pow(w/wPeak, 5)
    }
    else{
        var s = sp * Math.pow(w/wPeak, -2.5)  
    }
    
    s += 16 * tanh(wPeak/w) * swell   
    
    if (dirSpreadUseNumIntegration){
        return 1 / integrateCosine(2*s) * Math.pow(abs(half_angle_cos), 2*s)
    }
    
    else{
        var qs = Math.pow(2, 2*s-1) / pi * Math.pow(gamma(s+1), 2)/ gamma(2*s +1)
        return qs * Math.pow(abs(half_angle_cos), 2*s)
    } 
}
    
    
//Hasselmann + Horvath, 2015

export function hasselmannHorvath(w, waveDir, wPeak){
    var cos_angle = dot(windDir/norm(windDir),waveDir/norm(waveDir))
    var half_angle_cos = Math.sqrt((1+ cos_angle)/2)
    var sp = 11.5 * Math.pow((wPeak * windSpeed /g),-2.5)
    if (w <= wPeak){
        var s = 6.97 * Math.pow(w/wPeak, 4.06)
    }
    else{
        var s = 9.77 * Math.pow(w/wPeak, -2.33-1.45*((windSpeed * wPeak / g) - 1.17))  
    }

    s += 16 * tanh(wPeak/w) * swell
    if (dirSpreadUseNumIntegration){
        return 1 / integrateCosine(2*s) * Math.pow(abs(half_angle_cos), 2*s)
    }
    else{
        qs = Math.pow(2, 2*s-1) / pi * Math.pow(gamma(s+1), 2)/ gamma(2*s +1)
        return qs * Math.pow(abs(half_angle_cos), 2*s)
    }

}



// Donelan-Banner + Horvath, 2015

// aux fuction to compute Horvath
export function horvath(w, waveDir, wPeak){
    var cos_angle = dot(windDir/norm(windDir),waveDir/norm(waveDir))
    var half_angle_cos = Math.sqrt((1+ cos_angle)/2)
    var s = 16 * tanh(wPeak/w) * swell
    var qs = Math.pow(2, 2*s-1) / pi * Math.pow(gamma(s+1), 2)/ gamma(2*s +1)
    return qs * Math.pow(abs(half_angle_cos), 2*s)
}
    
    
//numerical integration to compute normalization factor
export function integrateDonelanHorvath(w, wPeak){
    var s = 0
    var theta, wDir;

    for(var i = 0; i < 30; i++){
        theta = -pi/2 + (i+1)*pi/32
        wDir = [cos(theta), sin(theta)]
        s += donelan(w, wDir, wPeak) * horvath(w, wDir, wPeak);
    }
    
    var wave0 =  [cos(-pi/2), sin(-pi/2)] 
    var f0 = donelan(w, wave0, wPeak) * horvath(w, wave0, wPeak)
    
    var waven =  [cos(pi/2), sin(pi/2)]
    var fn = donelan(w, waven, wPeak) * horvath(w, waven, wPeak)
    
    var res = (s + 0.5 * (f0 + fn))  * (2 * pi/32)   
    return res 
}
        
    
    
export function donelanHorvath(w, waveDir, wPeak){
    var res = donelan(w, waveDir, wPeak) * horvath(w, waveDir, wPeak)
    var normFactor = integrateDonelanHorvath(w, wPeak)
    
    var res = res / normFactor;
    return res
}    
        
                   