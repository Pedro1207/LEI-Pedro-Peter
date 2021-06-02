
// ------------------------------
// Dispersion
// ------------------------------

float getDispersionW(float k, float depth, int dispersionMode) {

	float tension = 0.074; // N/m
	float density = 1000; // Kg/m3
	
	// deep waters
	if (dispersionMode == 0) 
		return sqrt(k * G);
	// shallow waters	
	else if (dispersionMode == 1)	
		return sqrt(k * G * tanh(k*depth));
	else 
		return sqrt((G * k + pow(k,3) * tension/density)* tanh(k*depth));
}


float getDispersionDerivative(float k, float depth, int dispersionMode) {

	float tension = 0.074; // N/m
	float density = 1000; // Kg/m3
	float w = getDispersionW(k, depth, dispersionMode);
	
	// deep waters
	if (dispersionMode == 1) 
		return 0.5 * w  / k;
	// shallow waters	
	else if (dispersionMode == 2)	
		return 0.5 * (w/k) * (1 + k*depth * (1-pow(tanh(k*depth),2.0))/tanh(k*depth));
	else {
		float th = tanh(k*depth);
		float b = tension/density;
		float aTerm = G * k + pow(k,3) * tension/density;
		return ( (G + 3*b*k*k)*th + aTerm * k * (1-th*th))/sqrt(aTerm*th);
	}
}

// ------------------------------
// Directional Spreaading
// ------------------------------

float gamma_impl(const float s) {
	float s_inv = 1.0/s;
    //  Float version:
	const float g = 1.12906830989;
	const float c0 = 0.8109119309638332633713423362694399653724431;
	const float c1 = 0.4808354605142681877121661197951496120000040;
	const float e = 2.71828182845904523536028747135266249775724709;
	const float sph = s + 0.5;
	const float lanczos_sum = c0 + c1/(s + 1.0);
	const float base = (sph + g)/e;
	return (pow(base, sph) * lanczos_sum) * s_inv;
}

float getDirectionalSpreading(int directionalMode, float w, float wPeak, vec2 waveDir, vec2 windDir, float windSpeed, float swell) {

	switch(directionalMode) {
		// Tessendorf
		case 0:	return 2/M_PI * pow(dot(waveDir, windDir),2);
		// Occhi 1998 Ocean Waves, Cambridge Press
		// see Emprirical Directional Wave Spectra
		case 1: 
			float aux = dot(waveDir, windDir);
			if (aux > 0)
				return max(0, 2/M_PI * pow(aux,2));
			else
				return 0;
		// Mitsuyasu: see Emprirical Directional Wave Spectra
		case 2:
			float swellShape = 0;
			if (swell > 0)
				swellShape = w==0?0:16.1 * tanh(wPeak/w) * pow(swell,2);
			float s, sp = 11.5 * pow((wPeak * windSpeed /G),-2.5);
			if (w <= wPeak)
				s = sp * pow(w/wPeak, 5);
			else
				s = sp * pow(w/wPeak, -2.5);
			s += swellShape;
			float qs = pow(2, 2*s - 1)/ M_PI; 
			qs *= pow(gamma_impl(s + 1),2.0) / gamma_impl(2*s + 1);
			vec2 aux2 = waveDir + windDir;
			if (length(aux2) != 0)
				normalize(aux2);
			float c = 2 * pow(dot(waveDir, windDir),2) - 1;
			float res = qs * pow(abs(c),2*s);
			return res;
		// Hasselmann : see Emprirical Directional Wave Spectra
		case 3:
			swellShape = 0;
			if (swell > 0)
				swellShape = w==0?0:16.1 * tanh(wPeak/w) * pow(swell,2);
			sp = 11.5 * pow((wPeak * windSpeed /G),-2.5);
			if (w <= wPeak)
				s = 6.97 * pow(w/wPeak, 4.06);
			else
				s = 9.77 * pow(w/wPeak, -2.33-1.45*((windSpeed * wPeak / G) - 1.17));
			s += swellShape;
			qs = pow(2, 2*s - 1)/ M_PI; 
			qs *= pow(gamma_impl(s + 1),2.0) / gamma_impl(2*s + 1);
			aux2 = waveDir + windDir;
			if (length(aux2) != 0)
				normalize(aux2);
			c = 2 * pow(dot(waveDir, windDir),2) - 1;
			res = qs * pow(abs(c),2*s);
			return res;
		// Donelan-Banner : see Emprirical Directional Wave Spectra
		case 4:
			float beta;
			float wcoeff = w / wPeak;
			float epsilon = -0.4 + 0.8393 * exp( -0.567 * log(pow(wcoeff,2)));
			if (wcoeff < 0.95)
				beta = 2.61 * pow(wcoeff, 1.3);
			else if (wcoeff < 1.6)
				beta = 2.28 * pow(wcoeff, -1.3);
			else {
				float epsilon = -0.4 + 0.8393 * exp( -0.567 * log(pow(wcoeff,2)));
				beta = pow(10,epsilon);
			}
			res = beta / (2 * tanh(beta * M_PI)) * pow(1/cosh(beta * dot(waveDir, windDir)),2);
			return res;
			
	}	
}

// -----------------------------
// Spectra in wave number format
// -----------------------------

float Phillips(vec2 vec_k, float windSpeed, vec2 windDir) 
{
	if (vec_k == vec2(0.0f, 0.0f))
		return 0.0f;
        
    // Largest possible waves arising from a continuous wind of speed V
	float L = windSpeed*windSpeed / G; 

	float k = length(vec_k);
	vec2 k_hat = normalize(vec_k);
	float dot_k_w = dot(k_hat, normalize(windDir));
//	float result = 3e-6 * exp(-1 / (k*L*k*L)) / pow(k, 4) * pow(dot_k_w, 2);
	float result = 3.48e-3 * exp(-1 / (k*L*k*L)) / pow(k, 4) * pow(dot_k_w, 2);

    // suppressing very small waves ( l << L
	float l = 0.1;//0.00001 * L * L;
	result *= exp(-k*k*l*l);  // Eq24

	return result;
}

// ------------------------------
// Frequency Spectra
// ------------------------------


// http://dutch360hdr.com/downloads/other/EmpiricalDirectionalWaveSpectra_DigiPro2015_Jun20_2015.pdf
float Pierson_Moskowitz(vec2 vec_k, float windSpeed, float w) 
{
	// s(w)
	float w0 = G / (1.026 * windSpeed);
	float expTerm = -0.74 * pow(w0/w,4) ;
	float result = 8.1e-3 * G * G * pow(w, -5) * exp(expTerm);
	return result;
}


// https://wikiwaves.org/Ocean-Wave_Spectra
float JONSWAP(vec2 vec_k, float windSpeed, float w, float fetch) 
{
	if (vec_k == vec2(0.0f, 0.0f))
		return 0.0f;
        
	float aux1 = G * fetch / pow(windSpeed,2);

	float alpha = 0.076 * pow(aux1, -0.22);
	
	float wp = 22 * (G/windSpeed) * pow(aux1, -0.33);
	float expTerm = -1.25 * pow(wp/w, 4);
		
	float delta = 0.09;
	if (w <= wp)
		delta = 0.07;
	
	float gamaExp = pow(w - wp, 2) / (2 * pow(delta * wp, 2)); 
	
	float result = alpha * G * G * pow(w, -5) * exp(expTerm) * pow(3.3, exp(-gamaExp));
	
	return result;
}


float ModifiedJONSWAP(vec2 vec_k, float windSpeed, vec2 windDir, float depth, int dispersionMode, float fetch, float swell) 
{
	if (vec_k == vec2(0.0f, 0.0f))
		return 0.0f;
        
	float k = length(vec_k);
	vec2 k_hat = normalize(vec_k);

	float w = getDispersionW(k, depth, dispersionMode);
	
	float aux1 = G * fetch / pow(windSpeed,2);
	float wp = 22 * (G/windSpeed) * pow(aux1, -0.33);
	float v = wp * windSpeed  / (2* M_PI * G);
	float beta = 0.006 * pow(v, 0.55);
	
	float delta = 0.09;
	if (w <= wp)
		delta = 0.07;
		
	float gamma;	
	if (v >= 1)
		gamma = 6.489 + 6 * log(v);
	else
		gamma = 1.7;
	float gamaExp = pow(w - wp, 2) / (2 * pow(delta * wp, 2)); 
	float expTerm = pow(w/wp, 4);
	
	float result = beta * G * G * pow(wp, -1) * pow(w, -4) * exp(-expTerm) * pow(gamma, exp(-gamaExp));
	
	return result;
}

// http://dutch360hdr.com/downloads/other/EmpiricalDirectionalWaveSpectra_DigiPro2015_Jun20_2015.pdf
float TMA(vec2 vec_k, float windSpeed, float w, float depth, float fetch) {


	float result = JONSWAP(vec_k, windSpeed, w, fetch);

	float wh = clamp(w * sqrt(depth/G), 0.0, 2.0);
	float theta;
	if (wh <= 1)
		theta = 0.5 * wh*wh;
	else
		theta = 1-0.5*pow(2-wh, 2);
		
	clamp( theta,0,1);	
	
	return result * theta;
}

// Real-time Animation and Rendering of Ocean Whitecaps
// https://github.com/jdupuy/whitecaps
#define OMEGA 2 // sea state
#define km 370.0
#define cm 0.23
#define omnispectrum false

float sqr(float x)
{
    return x * x;
}

float omega(float k)
{
    return sqrt(9.81 * k * (1.0 + sqr(k / km))); // Eq 24
}

float WC(vec2  k_vec, float windSpeed, vec2 windDir, float depth, int dispersionMode, float swell) {

	if (length(k_vec) == 0)
		return 0;
	float A = 2;
    float U10 = windSpeed;
    float Omega = OMEGA;

    // phase speed
    float k = length(k_vec);
    float c = omega(k) / k;

    // spectral peak
    float kp = 9.81 * sqr(Omega / U10); // after Eq 3
    float cp = omega(kp) / kp;

    // friction velocity
    float z0 = 3.7e-5 * sqr(U10) / 9.81 * pow(U10 / cp, 0.9f); // Eq 66
    float u_star = 0.41 * U10 / log(10.0 / z0); // Eq 60

    float Lpm = exp(- 5.0 / 4.0 * sqr(kp / k)); // after Eq 3
    float gamma = Omega < 1.0 ? 1.7 : 1.7 + 6.0 * log(Omega); // after Eq 3 // log10 or log??
    float sigma = 0.08 * (1.0 + 4.0 / pow(Omega, 3.0f)); // after Eq 3
    float Gamma = exp(-1.0 / (2.0 * sqr(sigma)) * sqr(sqrt(k / kp) - 1.0));
    float Jp = pow(gamma, Gamma); // Eq 3
    float Fp = Lpm * Jp * exp(- Omega / sqrt(10.0) * (sqrt(k / kp) - 1.0)); // Eq 32
    float alphap = 0.006 * sqrt(Omega); // Eq 34
    float Bl = 0.5 * alphap * cp / c * Fp; // Eq 31

    float alpham = 0.01 * (u_star < cm ? 1.0 + log(u_star / cm) : 1.0 + 3.0 * log(u_star / cm)); // Eq 44
    float Fm = exp(-0.25 * sqr(k / km - 1.0)); // Eq 41
    float Bh = 0.5 * alpham * cm / c * Fm; // Eq 40

    Bh *= Lpm; 

    if (omnispectrum)
    {
        return A * (Bl + Bh) / (k * sqr(k)); // Eq 30
    }

    float a0 = log(2.0) / 4.0;
    float ap = 4.0;
    float am = 0.13 * u_star / cm; // Eq 59
    float Delta = tanh(a0 + ap * pow(c / cp, 2.5f) + am * pow(cm / c, 2.5f)); // Eq 57

    float phi = atan(k_vec.y, k_vec.x);

//    if (propagate == 1)
//    {
//        if (k_vec.x < 0.0)
//        {
//            return 0.0;
//        }
//        else
//        {
//            Bl *= 2.0;
//            Bh *= 2.0;
//        }
//    }

	// remove waves perpendicular to wind dir
	float tweak = sqrt(max(k_vec.x/k,0.0f));
	tweak = 1.0f;
    return A * (Bl + Bh) * (1.0 + Delta * cos(2.0 * phi)) / (2.0 * M_PI * sqr(sqr(k))) * tweak; // Eq 67

}


float getSpectrum(int spectrum, vec2 k_vec, float L, vec2 windDir, float windSpeed, int directionalMode, int dispersionMode, int propagate, float depth, float fetch, float swell) {


	if (k_vec == vec2(0,0))
		return 0;
		
	float res;
	float amplitudeCorrection = pow(2 * M_PI/L, 2);//2 * M_PI/L;//
	
	if (propagate == 1) {
		if (dot(-k_vec, windDir) <= 0.0)
		{
			return 0.0;
		}
		else
		{
			amplitudeCorrection *= 2;
		}
	}
	
	float k = length(k_vec);
	vec2 k_hat = normalize(k_vec);
	float w = getDispersionW(k, depth, dispersionMode);
	if (spectrum == 0)
		res = Phillips(k_vec, windSpeed, windDir);
	else if (spectrum == 1) {
		res = Pierson_Moskowitz(k_vec, windSpeed, w);
		res *= getDirectionalSpreading(directionalMode, w, 0.855*G/windSpeed, k_hat, windDir, windSpeed, swell);
		res *= getDispersionDerivative(k, depth, dispersionMode) / k;
		}
	else if (spectrum == 2) {
		res = JONSWAP(k_vec, windSpeed, w, fetch);
		float aux1 = G * fetch / pow(windSpeed,2);
		res *= getDirectionalSpreading(directionalMode, w, 22*(G/windSpeed) * pow(aux1, -0.33), k_hat, windDir, windSpeed, swell);
		res *= getDispersionDerivative(k, depth, dispersionMode) / k;
	}
	else if (spectrum == 3) {
		res = ModifiedJONSWAP(k_vec, windSpeed, windDir, depth, dispersionMode, fetch, swell);
		float aux1 = G * fetch / pow(windSpeed,2);
		res *= getDirectionalSpreading(directionalMode, w, 22*(G/windSpeed) * pow(aux1, -0.33), k_hat, windDir, windSpeed, swell); 
		res *= getDispersionDerivative(k, depth, dispersionMode) / k;
	}
	
		
	else if (spectrum == 4) {
		res = TMA(k_vec, windSpeed, w, depth, fetch);
		res *= getDirectionalSpreading(directionalMode, w, 22 * pow(G*G/(windSpeed*fetch), 0.33), k_hat, windDir, windSpeed, swell);
		res *= getDispersionDerivative(k, depth, dispersionMode) / k;
	}
	else if (spectrum == 5)
		res = WC(k_vec, windSpeed, windDir, depth, dispersionMode, swell);
		
	//amplitudeCorrection = 1;	
	return amplitudeCorrection * res ;
}


