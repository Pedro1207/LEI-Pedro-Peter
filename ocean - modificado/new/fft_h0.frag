#version 430

precision highp float;

in vec2 texCoordV;

#define M_PI 3.1415926535897932384626433832795
#define G 9.81


uniform int width;
uniform int L;

uniform vec2 windDir;
uniform float windSpeed;
uniform float depth;

uniform int randomDistribution; // Normal distribution
uniform int spectrum;
uniform int dispersionMode;
uniform int directionalMode;
uniform float spectrumScale;

uniform float fetch;
uniform float swell;
uniform float Hs;

// JONSWAP parameters
uniform float JONSWAP_gamma;
uniform float JONSWAP_sigmaA;
uniform float JONSWAP_sigmaB; 

// Bretschneider Param
uniform float Bretschneider_wm;

// Ochi parameters
uniform float Ochi_lambda1, Ochi_wm1, Ochi_Hs1, Ochi_lambda2, Ochi_wm2, Ochi_Hs2;


uniform int propagate;

uniform sampler2D noise_i0, texRnd;
layout (binding = 0, rgba32f) uniform image2D h0k; 
layout (binding = 1, rgba32f) uniform image2D texF; 

layout(std430, binding = 0) buffer oceanInfo{
	vec4 info[];
};


float noise(vec2 k) {
	return texture(noise_i0, k).r;
}

//http://amindforeverprogramming.blogspot.com/2013/07/random-floats-in-glsl-330.html
uint hash( uint x ) {
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}

// Compound versions of the hashing algorithm I whipped together.
uint hash( uvec2 v ) { return hash( v.x ^ hash(v.y)                         ); }

// Construct a float with half-open range [0:1] using low 23 bits.
// All zeroes yields 0.0, all ones yields the next smallest representable value below 1.0.
float floatConstruct( uint m ) {
    const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
    const uint ieeeOne      = 0x3F800000u; // 1.0 in IEEE binary32

    m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
    m |= ieeeOne;                          // Add fractional part to 1.0

    float  f = uintBitsToFloat( m );       // Range [1:2]
    return f - 1.0;                        // Range [0:1]
}

// Pseudo-random value in half-open range [0:1].
float random( uvec2  v ) { return floatConstruct(hash(v)); }


//#define TEXTURE_GRADS

//#define USE_NOISE

#define COMPUTE_SKY_FOR_REFLECTION

// foam
#define NO_FOAM 0
#define USE_JACOBIAN 1
#define USE_VERTICAL_ACCELERATION 2
#define FOAM USE_JACOBIAN

//#define TESTING_COLOR_DEPTHS
//#define TESTING_COLORS	

#define M_PI 3.1415926535897932384626433832795
#define G 9.81

#define RANDOM_UNIFORM 0
// normal dist can be computed with Box-Muller transform
#define RANDOM_NORMAL 1
// if X is from a normal distribution, then exp(X) follows a log distribution
#define RANDOM_LOG 2
// if X is uniform on [0,1] then −log(X) follows an exponential distribution
#define RANDOM_EXP 3


#define DISPERSION_DEEP 0
#define DISPERSION_SHALLOW 1
#define DISPERSION_ CAPILLARY 2 


#define DIRECTIONAL_COS_POWER 0
#define DIRECTIONAL_MITSUYASU 1
#define DIRECTIONAL_HASSELMANN 2
#define DIRECTIONAL_DONNELAN_BANNER 3
#define DIRECTIONAL_HORVATH_MITSUYASU 4
#define DIRECTIONAL_HORVATH_HASSELMANN 5
#define DIRECTIONAL_HORVATH_DONNELAN_BANNER 6

#define SPECTRUM_PHILLIPS 0
#define SPECTRUM_PIERSON_MOSKOWITZ 1
#define SPECTRUM_JONSWAP 2
#define SPECTRUM_DONNELAN_JONSWAP 3
#define SPECTRUM_TMA 4
#define SPECTRUM_UNIFIED 5
#define SPECTRUM_PIERSON_MOSKOWITZ_HS 6
#define SPECTRUM_BRETSCHNEIDER 7
#define SPECTRUM_OCHI 8
#define SPECTRUM_OCHI_HS 9



// ------------------------------
// Dispersion
// ------------------------------

float getDispersionW(float k, float depth, int dispersionMode) {

	float tension = 0.074; // N/m
	float density = 1000; // Kg/m3
	
	// deep waters
	if (dispersionMode == DISPERSION_DEEP) 
		return sqrt(k * G);
	// shallow waters	
	else if (dispersionMode == DISPERSION_SHALLOW)	
		return sqrt(k * G * tanh(k*depth));
	else  // DISPERSION_ CAPILLARY
		return sqrt((G * k + pow(k,3) * tension/density) * tanh(k*depth));
}


float getDispersionDerivative(float k, float depth, int dispersionMode) {

	float tension = 0.074; // N/m
	float density = 1000; // Kg/m3
	float w = getDispersionW(k, depth, dispersionMode);
	
	// deep waters
	if (dispersionMode == DISPERSION_DEEP) 
		return G * 0.5f / w; 
	// shallow waters	
	else if (dispersionMode == DISPERSION_SHALLOW) {
		float dk = depth * k;
		float th = tanh(dk);
		return 0.5f * sqrt(9.81f / (k*th)) * (th + dk*(1 - th*th));
	}
	else { // DISPERSION_ CAPILLARY
		float dk = depth * k;
		float th = tanh(dk);
		float b = tension/density;
		return 0.5f * ((9.81f + 3 * b*k*k)*th + dk * (k*k*b + 9.81f) * pow(1.0f / cosh(dk), 2)) / w;
	}
}






// ------------------------------
// Directional Spreaading
// ------------------------------

// Taken from:
// https://buildbot.libretro.com/assets/frontend/bundle/shaders/shaders_glsl/crt/shaders/crt-royale/port-helpers/special-functions.h

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


float numericalIntegrationCos(float two_s) {

	float theta;
	float steps = 64;
	float sum = 0.0f;
	for (int i = 1; i < steps - 1; ++i) {
		theta = -M_PI + i * (2 * M_PI) / steps;
		sum += pow(abs(cos(theta* 0.5)), two_s);
	}

	float s0 = 0;//pow(abs(cos(-M_PI * 0.5)), two_s);
	float sn = 0;//pow(abs(cos(M_PI * 0.5)), two_s);
	sum = (sum + 0.5 * (s0 + sn))  * (2 * M_PI) / steps;

	return sum;
}


float auxDonelan(float w, vec2 waveDir, vec2 windDir, float wPeak, float swell) {

	float cos_angle = dot(waveDir,windDir);
	if (cos_angle < -1)
		cos_angle = -1;
	else if (cos_angle > 1)
		cos_angle = 1;	
	float wave_angle = acos(cos_angle);
	float wcoeff = w / wPeak;
	float epsilon = -0.4f + 0.8393f * exp(-0.567f * log(pow(wcoeff, 2)));
	float beta;
	if (wcoeff < 0.95f)
		beta = 2.61f * pow(wcoeff, 1.3f);
	else if (wcoeff < 1.6f)
		beta = 2.28f * pow(wcoeff, -1.3f);
	else
		beta = pow(10.0f, epsilon);
	return beta / (2 * tanh(beta * M_PI)) * pow(1 / cosh(beta * wave_angle), 2);
}


float auxHorvath(float w, vec2 waveDir, vec2 windDir, float wPeak, float swell) {			// compute Horvath factor
	float cos_angle = dot(waveDir,windDir);
	if (cos_angle < -1)
		cos_angle = -1;
	float half_angle_cos = sqrt((1 + cos_angle) / 2);
	float s = 16 * tanh(wPeak / w) * swell;
	float qs = pow(2, 2 * s - 1) / M_PI * pow(gamma_impl(s + 1), 2) / gamma_impl(2 * s + 1);
	return qs * pow(abs(half_angle_cos), 2 * s);
}


float getDirectionalSpreading(int directionalMode, float w, float wPeak, vec2 waveDir, vec2 windDir, float windSpeed, float swell) {

	float cos_angle, half_angle_cos, s, sp, qs, wave_angle, wcoeff, epsilon, beta, norm, theta;

	switch(directionalMode) {
		// Massel (2018) - Ocean surface waves: their physics and precdiction
		//	"Historically, the first attempt to model directional energy spreading was
		//	suggested by Pierson et al. (1955) in the form of the cosine type function"
		case DIRECTIONAL_COS_POWER: 
			float cos_angle = dot(waveDir, windDir);
			if (cos_angle > 0)
				return 2/M_PI * pow(cos_angle,2);
			else
				return 0;

		// Mitsuyasu (1975)
		case DIRECTIONAL_MITSUYASU:
			cos_angle = dot(waveDir, windDir);
			if (cos_angle < -1)
				cos_angle = -1;
			half_angle_cos = sqrt(( 1 + cos_angle)/2);
			sp = 11.5 * pow((wPeak * windSpeed /G),-2.5);
			if (w <= wPeak)
				s = sp * pow(w / wPeak, 5);
			else
				s = sp * pow(w / wPeak, -2.5f);
			norm = numericalIntegrationCos(2*s);
			return 	pow(abs(half_angle_cos), 2*s) / (norm + 0.00001);
			//qs = pow(2, 2*s-1) / M_PI * pow(gamma_impl(s+1), 2)/ gamma_impl(2*s +1);
			//return qs * pow(abs(half_angle_cos), 2*s);

		// Mitsuyasu: see Emprirical Directional Wave Spectra, Horvath, 2015
		case DIRECTIONAL_HORVATH_MITSUYASU:
			cos_angle = dot(windDir,waveDir);
			if (cos_angle < -1)
				cos_angle = -1;
			half_angle_cos = sqrt((1+ cos_angle)/2);
			sp = 11.5 * pow((wPeak * windSpeed /G),-2.5);
			if (w <= wPeak)
				s = sp * pow(w / wPeak, 5);
			else
				s = sp * pow(w / wPeak, -2.5f);
			// swell term	
			s += 16 * tanh(wPeak/w) * swell;
			norm = numericalIntegrationCos(2*s);
			return 	pow(abs(half_angle_cos), 2*s) / norm;
			//qs = pow(2, 2*s-1) / M_PI * pow(gamma_impl(s+1), 2)/ gamma_impl(2*s +1);
			//return qs * pow(abs(half_angle_cos), 2*s);

		// Hasselmann 1980
		case DIRECTIONAL_HASSELMANN:

			cos_angle = dot(windDir,waveDir);
			if (cos_angle < -1)
				cos_angle = -1;
			half_angle_cos = sqrt((1+ cos_angle) * 0.5);
			sp = 11.5 * pow((wPeak * windSpeed /G),-2.5);
			if (w <= wPeak)
				s = 6.97f * pow(w / wPeak, 4.06f);
			else
				s = 9.77f * pow(w / wPeak, -2.33f - 1.45f*((windSpeed * wPeak / 9.81f) - 1.17f));
			norm = numericalIntegrationCos(2*s);
			return 	pow(abs(half_angle_cos), 2*s) / norm;
			//qs = pow(2, 2*s-1) / M_PI * pow(gamma_impl(s+1), 2)/ gamma_impl(2*s +1);
			//return qs * pow(abs(half_angle_cos), 2*s);
			
		// Hasselmann: see Emprirical Directional Wave Spectra, Horvath, 2015
		case DIRECTIONAL_HORVATH_HASSELMANN:

			cos_angle = dot(windDir,waveDir);
			if (cos_angle < -1)
				cos_angle = -1;
			half_angle_cos = sqrt((1+ cos_angle)/2);
			sp = 11.5 * pow((wPeak * windSpeed /G),-2.5);
			if (w <= wPeak)
				s = 6.97f * pow(w / wPeak, 4.06f);
			else
				s = 9.77f * pow(w / wPeak, -2.33f - 1.45f*((windSpeed * wPeak / 9.81f) - 1.17f));
			// swell term	
			s += 16 * tanh(wPeak/w) * swell;
			norm = numericalIntegrationCos(2*s);
			return 	pow(abs(half_angle_cos), 2*s) / norm;
			//qs = pow(2, 2*s-1) / M_PI * pow(gamma_impl(s+1), 2)/ gamma_impl(2*s +1);
			//return qs * pow(abs(half_angle_cos), 2*s);
			
		// Donelan-Banner : see Emprirical Directional Wave Spectra
		case DIRECTIONAL_DONNELAN_BANNER:
			cos_angle = dot(windDir,waveDir);
			if (cos_angle < -1)
				cos_angle = -1;
			else if (cos_angle > 1)
				cos_angle = 1;	
			wave_angle = acos(cos_angle);
			wcoeff = w/wPeak;
			epsilon = -0.4 + 0.8393 * exp( -0.567 * log(pow(wcoeff,2)));


			if (wcoeff < 0.95)
				beta = 2.61 * pow(wcoeff, 1.3);
			else if (wcoeff < 1.6)
				beta = 2.28 * pow(wcoeff, -1.3);
			else
				beta = pow(10,epsilon);
			return beta / (2 * tanh(beta * M_PI)) * pow(1 / cosh(beta * wave_angle), 2);;

		// Donelan-Banner: see Emprirical Directional Wave Spectra, Horvath, 2015	
		case DIRECTIONAL_HORVATH_DONNELAN_BANNER:
			vec2 wDir; float sum = 0;
			for (int i = 1; i < 30; ++i) {
				theta = -M_PI + (i + 1)* 2 * M_PI / 32;
				wDir = vec2(cos(theta), sin(theta));
				sum += auxDonelan(w, wDir, windDir, wPeak, swell) * auxHorvath(w, wDir, windDir, wPeak, swell);
			}
			wDir = vec2(cos(-M_PI), sin(-M_PI));
			float f0 = auxDonelan(w, wDir, windDir, wPeak, swell) * auxHorvath(w, wDir, windDir, wPeak, swell);

			wDir = vec2(cos(M_PI), sin(M_PI));
			float fn = auxDonelan(w, wDir, windDir, wPeak, swell) * auxHorvath(w, wDir, windDir, wPeak, swell);
			
			norm = (sum + 0.5f * (f0 + fn))  * (2 * M_PI / 32);

			return auxDonelan(w, waveDir, windDir, wPeak, swell) * auxHorvath(w, waveDir, windDir, wPeak, swell) / norm;
	
		
	}
	
}

// -----------------------------
// Spectra in wave number format
// -----------------------------
// Chen, Li-ning and Jin, Yi - cheng and Yin, Yong and Ren, Hong-xiang
// On the Wave Spectrum Selection in Ocean Wave Scene Simulation of the Maritime Simulator
// AsiaSim 2013     A = 3.48e-3f

float Phillips(vec2 vec_k, float windSpeed, vec2 windDir) 
{
	if (vec_k == vec2(0.0f, 0.0f))
		return 0.0f;
        
    // Largest possible waves arising from a continuous wind of speed V
	float L = windSpeed*windSpeed / G; 

	float k = length(vec_k);
	vec2 k_hat = normalize(vec_k);
	float dot_k_w = dot(k_hat, normalize(windDir));
	float result = 3.48e-3 * exp(-1 / (k*L*k*L)) / pow(k, 4) * pow(dot_k_w, 2);

    // suppressing very small waves ( l << L
	float l = 0.1;//0.00001 * L * L;
	result *= exp(-k*k*l*l);  // Eq24

	return result;
}

// ------------------------------
// Frequency Spectra
// ------------------------------


// From original article by Pierson and Moskowitz (1964) Journal of Geophysical Research
// The scaling coefficient 1.026 is to compensate the fact that in PM the wind is 
// measured at 19.5 mt, whereas in all other works a height of 10 mts is considered
// U19.5 ≈ 1.026 U10 (Stewart, 2006)
// the function below takes winds at 10 mt above sea level
float Pierson_Moskowitz(float windSpeed, float w) 
{
	// s(w)
	float wp = 0.879 * G / (1.026 * windSpeed);
	float expTerm = -1.25 * pow(wp/w,4) ;
	float result = 8.1e-3 * G * G * pow(w, -5) * exp(expTerm);
	return result;
}


// Pierson - Moskowitz defined as a function of significant wave height
// Ochi, pg 35
// https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/lec6_wavespectra.pdf
float PiersonMoskowitzHs(float w) {

	float expTerm = -0.032 * pow(G / (Hs), 2) * pow(w, -4);
	float result = 8.1e-3 * pow(G, 2) * pow(w, -5) * exp(expTerm);
	return result;
}


// Bretschneider two parameter spectrum: Hs and wm
// setting wm = 0.4 sqrt(g/Hs) reverts to Pierson-Moskowitz
// Ochi, pg 36
float Bretschneider(float w) {

	float expTerm = -1.25f * pow(Bretschneider_wm / w, 4);
	float result = 1.25f/4.0f * pow(Bretschneider_wm,4) * Hs*Hs * pow(w, -5) * exp(expTerm);
	return result;
}


// Ochi, Six-Parameter wave spectra, 1976
float OchiSingle(float w, float l, float wm, float hs) {

	float aux = (4 * l + 1) / 4.0f;
	return 0.25f * pow(aux * pow(wm, 4), l) / gamma_impl(l) * pow(hs, 2) / pow(w, 4 * l + 1) * exp(-aux* pow(wm / w, 4));
}


float Ochi(float w, float l1, float wm1, float hs1, float l2, float wm2, float hs2) {

	return OchiSingle(w, l1, wm1, hs1) + OchiSingle(w, l2, wm2, hs2);
}


float OchiHs(float w) {
	return Ochi(w, 3, 0.7f * exp(-0.046f * Hs), 0.84f * Hs,
		1.54f * exp(-0.062f *  Hs), 1.15f * exp(-0.039f * Hs), 0.54f* Hs);
}


// From Hasselman et. al 1973
float JONSWAP(float windSpeed, float w, float fetch) 
{    
	float dimensionlessFetch = G * fetch / pow(windSpeed,2);
	float wp = 22 * (G/windSpeed) * pow(dimensionlessFetch, -0.33);
	float expTerm = -1.25 * pow(wp/w, 4);
	
	float sigma = JONSWAP_sigmaB;
	if (w <= wp)
		sigma = JONSWAP_sigmaA;	
	float gamaExp = 0.5 * pow((w - wp) /(sigma * wp), 2); 
	
	float alpha = 0.076 * pow(dimensionlessFetch, -0.22);
	
	float result = alpha * G * G * pow(w, -5) * exp(expTerm) * pow(JONSWAP_gamma, exp(-gamaExp));
	
	return result;
}


// Donelan et al (1985) proposed a modification to the JONSWAP spectral form
float DonnelanJONSWAP(float windSpeed, float w, float fetch) 
{
	float dimensionlessFetch = G * fetch / pow(windSpeed, 2);
	float Omega = 11.6f * pow(dimensionlessFetch, -0.23f);

	// wp from JONSWAP
	float wp = 22 * (G/windSpeed) * pow(dimensionlessFetch, -0.33);
	Omega = wp * windSpeed / G;
	//float wp = Omega * (G / windSpeed);
	float expTerm = -pow(wp / w, 4);

	float sigma = 0.08f * (1 + 4 / pow(Omega,3));

	float gamaExp = 0.5 * pow((w - wp) /(sigma * wp), 2); 

	float beta = 0.006f * pow(Omega, 0.55f);

	float gamma;
	if (Omega > 1)
		gamma = 1.7f + 6 * log(Omega);
	else
		gamma = 1.7f;

	float result = beta * G * G * pow(wp, -1) * pow(w, -4) * exp(expTerm) * pow(gamma, exp(-gamaExp));

	return result;
}


// From TMA, Hughes et al 1984
float TMA(float windSpeed, float w, float depth, float fetch) {


	float result = JONSWAP(windSpeed, w, fetch);

	float wh = clamp(w * sqrt(depth/G), 0.0, 2.0);
	float theta;
	if (wh <= 1)
		theta = 0.5 * wh*wh;
	else
		theta = 1-0.5*pow(2-wh, 2);
		
	clamp( theta,0,1);	
	
	return result * theta;
}


// From A unified directional spectrumfor long and short wind - driven waves
// T.Elfouhaily, B.Chapron, and K.Katsaros
// https://archimer.ifremer.fr/doc/00091/20226/
float omega(float k) {
	float km = 370;
	return (sqrt(G * k * (1 + pow(k / km, 2)))); // eq 24.
}


float Unified(vec2 kvec, vec2 khat) {

	float k = length(kvec);
	if (k == 0.0f)
		return 0.0f;

	float X = G * fetch / pow(windSpeed, 2); // after eq. 4

	float X0 = 2.2e4f;
	float Omega = 0.84f * pow(tanh(pow(X / X0, 0.4f)), -0.75f); // eq. 37

	float kp = G * pow(Omega / windSpeed, 2); // below eq. 3

	float c = omega(k) / k; // wave phase speed
	float cp = omega(kp) / kp; // phase speed for peak frequency

	float alphap = 6e-3f * sqrt(Omega); // eq. 34

	float Lpm = exp(-1.25f * pow(kp / k, 2));  // eq. 2

	float sigma = 0.08f * (1 + 4 * pow(Omega, -3)); // below eq. 3
	float Gamma = exp(-0.5f * pow(sqrt(k / kp) - 1, 2) / pow(sigma, 2));
	float gamma = 1.7f;
	if (Omega > 1)
		gamma +=  6 * log(Omega);

	float Jp = pow(gamma, Gamma); // eq. 3

	float Fp = Lpm * Jp * exp(-Omega / sqrt(10) * (sqrt(k / kp) - 1)); // eq. 32

	float Bl = 0.5f * alphap * (cp / c) * Fp; // eq. 31

	float cm = 0.23f; // between eq. 40 and 41

	float z0 = 3.7e-5f * pow(windSpeed, 2) / G * pow(windSpeed / cp, 0.9f); // eq. 66

	// https://en.wikipedia.org/wiki/Von_K%C3%A1rm%C3%A1n_constant
	float K = 0.41f;
	float ustar = K * windSpeed / log(10.0f / z0); // eq. 60

	float alpham = 0.01f;
	if (ustar < cm) // eq 44
		alpham *= (1 + log(ustar / cm));
	else
		alpham *= (1 + 3 * log(ustar / cm));

	float km = 370;
	float Fm = exp(-0.25f * pow(k / km - 1, 2));

	float Bh = 0.5f * alpham * (cm / c) * Fm * Lpm; // eq. 40 (fixed ? )

	float a0 = log(2) / 4;
	float ap = 4;
	float am = 0.13f * ustar / cm;
	float delta = tanh(a0 + ap * pow(c / cp, 2.5f) + am * (pow(cm / c, 2.5f))); // eq. 57

	float cos_twice_angle = 2 * pow(dot(windDir,khat), 2) - 1;
	return 1 / (2 * M_PI) * pow(k, -4) * (Bl + Bh) * (1 + delta * cos_twice_angle); // eq. 67
}



float getSpectrum(int spectrum, vec2 k_vec, float L, vec2 windDir, float windSpeed, int directionalMode, int dispersionMode, int propagate, float depth, float fetch, float swell) {


	if (k_vec == vec2(0,0))
		return 0;
		
	float res;
	float amplitudeCorrection = pow(2 * M_PI/L, 2);//2 * M_PI/L;//
	
	if (propagate == 1) {
		if (dot(k_vec, windDir) <= 0.0) {
			return 0.0;
		}
		else {
			amplitudeCorrection *= 2;
		}
	}
	
	float k = length(k_vec);
	vec2 k_hat = normalize(k_vec);
	float aux1;
	float w = getDispersionW(k, depth, dispersionMode);
		
	//imageStore(texF, ivec2(gl_GlobalInvocationID.xy), vec4(directionalMode, spectrum, 3, 4));

	switch(spectrum) {
		case SPECTRUM_PHILLIPS:
			res = Phillips(k_vec, windSpeed, windDir);
			break;
		case SPECTRUM_PIERSON_MOSKOWITZ:
			res = Pierson_Moskowitz(windSpeed, w);
			res *= getDirectionalSpreading(directionalMode, w, 0.877*G/windSpeed, k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_JONSWAP:
			res = JONSWAP(windSpeed, w, fetch);
			aux1 = G * fetch / pow(windSpeed,2);
			res *= getDirectionalSpreading(directionalMode, w, 22*(G/windSpeed) * pow(aux1, -0.33), k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_DONNELAN_JONSWAP:
			res = DonnelanJONSWAP(windSpeed, w, fetch);
			aux1 = G * fetch / pow(windSpeed,2);
			res *= getDirectionalSpreading(directionalMode, w, 22*(G/windSpeed) * pow(aux1, -0.33), k_hat, windDir, windSpeed, swell); 
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_TMA:
			res = TMA(windSpeed, w, depth, fetch);
			res *= getDirectionalSpreading(directionalMode, w, 22 * pow(G*G/(windSpeed*fetch), 0.33), k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_UNIFIED:
			res = Unified(k_vec, k_hat);
			break;
		case SPECTRUM_PIERSON_MOSKOWITZ_HS:
			res = PiersonMoskowitzHs(w);
			res *= getDirectionalSpreading(directionalMode, w, 0.4f * sqrt(G / Hs), k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_BRETSCHNEIDER:
			res = Bretschneider(w);
			res *= getDirectionalSpreading(directionalMode, w, 0.4f * sqrt(G / Hs), k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_OCHI:
			res = Ochi(w,Ochi_lambda1, Ochi_wm1, Ochi_Hs1, Ochi_lambda2, Ochi_wm2, Ochi_Hs2);
			res *= getDirectionalSpreading(directionalMode, w, 0.4f * sqrt(G / Hs), k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;
		case SPECTRUM_OCHI_HS:
			res = OchiHs(w);
			res *= getDirectionalSpreading(directionalMode, w, 0.4f * sqrt(G / Hs), k_hat, windDir, windSpeed, swell);
			res *= getDispersionDerivative(k, depth, dispersionMode) / k;
			break;

	}	
	//amplitudeCorrection = 1;	
	return amplitudeCorrection * res ;
}



void main(void) {

	vec2 hK, hminusK;
	// reset max and min depth;
	info[0] = vec4(0);
	
	ivec2 pos = ivec2(texCoordV * width);
	
	int kx = pos.x >= width/2 ? pos.x - width: pos.x;
	int kz = pos.y >= width/2 ? pos.y - width: pos.y;
	float whalf = width / 2.0;

	vec2 k = vec2(kx, kz) * 2 * M_PI/L;
	
	float aux = getSpectrum(spectrum, k, L, windDir, windSpeed, directionalMode, 			
					dispersionMode, propagate, depth, fetch, swell) * spectrumScale;
	
	float p = sqrt(aux * 0.5f);

	float z1, z2;

	// clamp to prevent values equal to zero
	vec4 z = clamp(texelFetch(texRnd, ivec2(texCoordV * width), 0), 0.000001, 1.0);
	z1 = z.x-0.5;
	z2 = z.y-0.5;

	if (randomDistribution == RANDOM_NORMAL || randomDistribution == RANDOM_LOG) {
		// Box-Muller Transform
		// http://mathworld.wolfram.com/Box-MullerTransformation.html
		z1 = sqrt(-2*log(z.x))*cos(2 * M_PI * z.y);
		z2 = sqrt(-2*log(z.x))*sin(2 * M_PI * z.y);

		// https://www.quora.com/How-do-I-transform-between-log-normal-distribution-and-normal-distribution
		if (randomDistribution == RANDOM_LOG) {
			z1 = exp(z1);
			z2 = exp(z2);
		}
	}
	// https://stats.stackexchange.com/questions/234544/from-uniform-distribution-to-exponential-distribution-and-vice-versa
	else if (randomDistribution == RANDOM_EXP) {
		z1 = -log(z.x);
		z2 = -log(z.y);
	}
	hK = vec2(z1, z2) * p;
	//hK = sqrt(0.5) * vec2(z.x, z.y) * p;
	
	// Box-Muller Transform
	float aux2 = getSpectrum(spectrum, -k, L, windDir, windSpeed, directionalMode, dispersionMode, 
					propagate, depth, fetch, swell) * spectrumScale;
	p = sqrt(aux2 * 0.5);
	
	hminusK = vec2(z1, z2) * p;
	//hminusK = sqrt(0.5) * vec2(z.z, z.w) * p;
	hminusK.y = -hminusK.y;
	
	imageStore(h0k, pos ,vec4(hK, hminusK));
    discard;

	//imageStore(h0k, pos ,teste);
	//imageStore(h0k, pos ,vec4(dispersionMode, directionalMode, spectrum, spectrumScale));
}


