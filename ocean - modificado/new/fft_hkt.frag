#version 430

layout (binding = 1, rgba32f)   writeonly uniform image2D texF; 

// layers 
// VA is JXY
#define LAYER_Y_JXY_JXX_JYY 0
#define LAYER_DX_DZ_SX_SZ 1

// cascade layers
#define LAYER_Y	        0
#define LAYER_DX	    1
#define LAYER_DZ        2
#define LAYER_SX	    3
#define LAYER_SZ	    4
#define LAYER_VA_JXY	5
#define LAYER_JXX   	6
#define LAYER_JYY   	7

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


in vec2 texCoordV;

layout (binding = 0, rgba32f) writeonly uniform image2DArray tilde_hkt; 

uniform sampler2D tilde_h0k;

uniform int width;
uniform int L;
uniform float timer;
uniform int dispersionMode;
uniform float depth;


vec2 mult(vec2 v0, vec2 v1) {
	return vec2(v0.x * v1.x - v0.y * v1.y,
				v0.x * v1.y + v0.y * v1.x);
}

vec2 conj(vec2 v){
	return vec2(v.x, -v.y);
}

void main(void) {

//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_DY), 
//		texelFetch(tilde_h0k, ivec2(gl_GlobalInvocationID.xy),0));
//	return;	
	float t = (timer+1e7) / 1000;
	//float t = timer/1000;
	ivec2 pos = ivec2(texCoordV * 512);

	int kx = pos.x >= width/2 ? pos.x - width: pos.x;
	int kz = pos.y >= width/2 ? pos.y - width: pos.y;
	float whalf = width / 2.0;

	vec2 k = vec2(kx, kz) * 2 * M_PI/L;
	
	float magnitude = length(k);
//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_DY), vec4(magnitude,depth,k));
//	return;	
	if (magnitude < 0.00000000001) magnitude = 0.00000000001;
	
	float w = getDispersionW(magnitude, depth, dispersionMode);
	
	vec4 spectrum = texelFetch(tilde_h0k, ivec2(texCoordV * 512),0);
	int x,y;
	x = (width - int(texCoordV.x * 512)) % width;
	y = (width - int(texCoordV.y * 512)) % width;
	vec4 spectrumC = texelFetch(tilde_h0k, ivec2(x,y),0);

	vec2 fourier_amp = spectrum.xy;
	vec2 fourier_amp_conj = spectrumC.xy;
	fourier_amp_conj.y = -fourier_amp_conj.y;
	
	float cosinus = cos(w*t);
	float sinus   = sin(w*t);
	vec2 exp_iwt = vec2(cosinus, sinus);
	vec2 exp_iwt_inv = vec2(cosinus, -sinus);

	// dy
	vec2 h_k_t_dy = mult(fourier_amp, exp_iwt) + mult(fourier_amp_conj, exp_iwt_inv);
//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_DY), vec4(h_k_t_dy, cosinus, sinus));
//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_DY), vec4(fourier_amp, fourier_amp_conj));
//	return;	

	// dx
	vec2 dx = vec2(0.0, -k.x/magnitude);
	vec2 h_k_t_dx = mult(dx, h_k_t_dy);
	
	// dz
	vec2 dz = vec2(0.0, -k.y/magnitude);
	vec2 h_k_t_dz = mult(dz, h_k_t_dy);
	
	// sx
	dx = vec2(0.0, k.x);
	vec2 sx = mult(dx, h_k_t_dy);

	// sz
	dz = vec2(0.0, k.y);
	vec2 sz = mult(dz, h_k_t_dy);
	
#if (FOAM == USE_VERTICAL_ACCELERATION)

	vec2 a0k = fourier_amp * w * w;
	vec2 a0minusk = fourier_amp_conj * w * w;
	vec2 a = mult(a0k, exp_iwt) + mult(a0minusk, exp_iwt_inv);
	vec2 b = vec2(0);
	
#elif (FOAM == USE_JACOBIAN)

	vec2 daux = vec2(0.0, -k.x);
	vec2 jxx = mult(daux, h_k_t_dx);
	vec2 a = mult(daux, h_k_t_dz);
	daux = vec2(0.0, -k.y);
	vec2 jyy = mult(daux, h_k_t_dz);
	vec2 b = jxx +vec2(-jyy.y, jyy.x);
//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_JXXYY), vec4(jxx, jyy));
	
#else
	vec2 a = vec2(0);
#endif
	vec2 dy = h_k_t_dy + vec2(-a.y, a.x);
	vec2 dxz = h_k_t_dx + vec2(-h_k_t_dz.y, h_k_t_dz.x);
	vec2 sxz = sx + vec2(-sz.y, sz.x);
	imageStore(tilde_hkt, ivec3(texCoordV * 512, LAYER_Y_JXY_JXX_JYY), vec4(dy, b));
	imageStore(tilde_hkt, ivec3(texCoordV * 512, LAYER_DX_DZ_SX_SZ), vec4(dxz,sxz));
    discard;
//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_DXZ), vec4(h_k_t_dx, h_k_t_dz));
//	imageStore(tilde_hkt, ivec3(gl_GlobalInvocationID.xy, LAYER_SXZ), vec4(sx, sz));
	
}