
uniform sampler2DArray htk;
uniform sampler2D voronoi, foam;
uniform vec3 camPos;

uniform float choppyFactor;

const float indAir = 1.000293; //air refraction index
const float indWater = 1.333; //water index of refraction
const float Eta = indAir/indWater;
uniform float power = 5.0;


in Data {
	vec3 l_dir;
	vec3 pos;
	vec2 texCoord;
	vec3 normal;
} DataIn;

out vec4 outputF;


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


// Sky stuff

uniform sampler2D sky;

#define LINEAR 0
#define EXPONENTIAL 1
uniform int sampling = LINEAR;

uniform int divisions = 8;
uniform int divisionsLightRay = 4;
uniform float exposure = 1.5;

uniform vec3 betaR = vec3(3.67044e-07, 1.11688e-06, 1.80601e-06);
uniform float betaMf = 5.76e-07;
uniform float Hr = 7994;
uniform float Hm = 1200;
uniform float g = 0.99;
uniform vec2 sunAngles;

const float PI = 3.14159265358979323846;
const float earthRadius = 6360000;
const float atmosRadius = 6420000;
const float fourPI = 4.0 * PI;


///////////////////////////////////////////////////////////////////////////
// 					Sky stuff 
///////////////////////////////////////////////////////////////////////////

vec3 skyColor(vec3 dir, vec3 sunDir, vec3 origin);


vec4 computeSkyReflection(vec3 refl) {

	vec2 sunAnglesRad = vec2(sunAngles.x, sunAngles.y) * vec2(PI/180);
	vec3 sunDir = vec3(cos(sunAnglesRad.y) * sin(sunAnglesRad.x),
							 sin(sunAnglesRad.y),
							-cos(sunAnglesRad.y) * cos(sunAnglesRad.x));
							
#ifdef COMPUTE_SKY_FOR_REFLECTION		
	return vec4(skyColor(refl, sunDir, vec3(0.0, earthRadius+100, 0.0)),1);
#else	
	float phi = atan(refl.z, refl.x);
	float theta = acos(refl.y);
	float aux = tan(phi);
	float x = sqrt((1-cos(theta))/(1+aux*aux));
	float y = aux*x;
	vec2 tcSky = vec2(x, y);
	float ka = length(tcSky);
	if (ka >= 0.99) 
		tcSky *= 0.99/ka;
//	tcSky.x = 1 - tcSky.x;
	tcSky = tcSky * 0.5 + 0.5;
	return texture(sky, tcSky );
#endif
}



float distToTopAtmosphere(vec3 origin, vec3 dir) {

	// project the center of the earth on to the ray
	vec3 u = vec3(-origin);
	// k is the signed distance from the origin to the projection
	float k = dot(dir,u);
	vec3 proj = origin + k * dir;
	
	// compute the distance from the projection to the atmosphere
	float aux = length(proj); 
	float dist = sqrt(atmosRadius * atmosRadius - aux*aux);
	
	dist += k;	
	return dist;
}


void initSampling(in float dist, in int div, out float quotient, out float segLength) {

	if (sampling == EXPONENTIAL) {
		quotient =  pow(dist, 1.0/(div));
		//segLength = quotient - 1;
	}
	else { // linear sampling
		segLength = dist/div;
	}
}


void computeSegLength(float quotient, float current, inout float segLength) {

	if (sampling == EXPONENTIAL) {
		segLength = current * quotient - current;
	}
	else { // linear sampling
	}
}


vec3 skyColor(vec3 dir, vec3 sunDir, vec3 origin) {

	float dist = distToTopAtmosphere(origin, dir);

	float quotient, quotientLight, segLengthLight, segLength;
	
	float cosViewSun = dot(dir, sunDir);
	
	vec3 betaM = vec3(betaMf);
	
	vec3 rayleigh = vec3(0);
	vec3 mie = vec3(0);
	
	float opticalDepthRayleigh = 0;
	float opticalDepthMie = 0;

	// phase functions
	float phaseR = 0.75 * (1.0 + cosViewSun * cosViewSun);

	float aux = 1.0 + g*g - 2.0*g*cosViewSun;
	float phaseM = 3.0 * (1 - g*g) * (1 + cosViewSun * cosViewSun) / 
					(2.0 * (2 + g*g) * pow(aux, 1.5)); 

	float current = 1;
	initSampling(dist, divisions, quotient, segLength);
	float height;
	for(int i = 0; i < divisions; ++i) {
		computeSegLength(quotient, current, segLength);
		vec3 samplePos = origin + (current + segLength * 0.5) * dir;
		height = length(samplePos) - earthRadius;
		if (height < 0) {
			break;
		}
		float hr = exp(-height / Hr) * segLength;
		float hm = exp(-height / Hm) * segLength;
		opticalDepthRayleigh += hr;
		opticalDepthMie += hm;
		
		float distLightRay = distToTopAtmosphere(samplePos, sunDir);
		initSampling(distLightRay, divisionsLightRay, quotientLight, segLengthLight);
		float currentLight = 1;
		float opticalDepthLightR = 0;
		float opticalDepthLightM = 0;
		int j = 0;
		
		for (; j < divisionsLightRay; ++j) {
			computeSegLength(quotientLight, currentLight, segLengthLight);
			vec3 sampleLightPos = samplePos + (currentLight + segLengthLight * 0.5) * sunDir;
			float heightLight = length(sampleLightPos) - earthRadius;
			if (heightLight < 0){
				break;
			}

			opticalDepthLightR += exp(-heightLight / Hr) * segLengthLight;
			opticalDepthLightM += exp(-heightLight / Hm) * segLengthLight;
			currentLight += segLengthLight;

		}

		if (j == divisionsLightRay) {
			vec3 tau = fourPI * betaR * (opticalDepthRayleigh + opticalDepthLightR) + 
					   fourPI * 1.1 * betaM *  (opticalDepthMie + opticalDepthLightM);
			vec3 att = exp(-tau);
			rayleigh += att * hr;
			mie += att * hm;
		}

		current += segLength;
	}
	vec3 result = (rayleigh *betaR * phaseR + mie * betaM * phaseM) * 20;
	vec3 white_point = vec3(1.0);
	result = pow(vec3(1.0) - exp(-result / white_point * exposure), vec3(1.0 / 2.2));

	return result;
}



float schlickRatio (vec3 rayDirection, vec3 normal) {

	float f =  pow((1.0 - indWater) / (1.0 + indWater) , 2);
	float schlick = f + (1 - f) * pow(1 - dot(-rayDirection,normal), power);
	
	return clamp(schlick, 0 ,1);
}


// From white-caps master Dupuy and Bruneton
vec3 hdr(vec3 L) {
    L = L * 1.05;//hdrExposure;
    L.r = L.r < 1.413 ? pow(L.r * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.r);
    L.g = L.g < 1.413 ? pow(L.g * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.g);
    L.b = L.b < 1.413 ? pow(L.b * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.b);
    return L;
}





float computeFoamFactor() {

	float f = 0;
	
#if (FOAM == USE_VERTICAL_ACCELERATION)

#define minFoam 0
#define maxFoam 7
	float whiteCap = texture(htk, vec3(DataIn.texCoord, LAYER_Y_JXY_JXX_JYY)).y * choppyFactor;
	vec4 foamV = texture(foam, DataIn.texCoord*2);
	f = pow(smoothstep(1,7, whiteCap), 2.0);
	f = 2*f;
	outputF = outputF * (1-f) + foamV * f;

#elif (FOAM == USE_JACOBIAN)
	float jxx= 1, jyy = 1, jxy = 0;
	jxx += texture(htk, vec3(DataIn.texCoord, LAYER_Y_JXY_JXX_JYY)).z * choppyFactor;
	jyy += texture(htk, vec3(DataIn.texCoord, LAYER_Y_JXY_JXX_JYY)).w * choppyFactor;
	jxy += texture(htk, vec3(DataIn.texCoord, LAYER_Y_JXY_JXX_JYY)).y * choppyFactor;

	float det = jxx * jyy - jxy*jxy;
	float whiteCap = det;
	
	vec4 foamV = texture(foam, DataIn.texCoord*2);
	f = 1-smoothstep(0.0, 0.7, whiteCap);
	if (whiteCap < 0.0)
		f = 1;
#endif	
	return f;
}


vec4 computeOceanColor(vec3 wn) {

	vec2 sunAnglesRad = vec2(sunAngles.x, sunAngles.y) * vec2(M_PI/180);
	vec3 sunDir = vec3(cos(sunAnglesRad.y) * sin(sunAnglesRad.x),
							 sin(sunAnglesRad.y),
							-cos(sunAnglesRad.y) * cos(sunAnglesRad.x));

	vec3 viewDir = normalize(DataIn.pos - camPos);
	vec3 reflDir = normalize(reflect(viewDir, wn));	
	if (reflDir.y < 0)
		reflDir.y = -reflDir.y;

	float spec = pow (max(0, dot(reflDir, sunDir)), 128);
	
	vec4 skyC = computeSkyReflection(reflDir);	
	
	vec4 shallowColor = vec4(0.0, 0.64, 0.68, 1);
	vec4 deepColor = vec4(0.02, 0.05, 0.10, 1);

	float relativeHeight = clamp((DataIn.pos.y - (-40)) / (80), 0.0, 1.0);
	vec4 heightColor = (relativeHeight * shallowColor + (1 - relativeHeight) * deepColor) * 0.8;
	float ratio = schlickRatio(viewDir, wn);
	float refCoeff = pow(max(dot(wn, -viewDir), 0.0), 0.3);	// Smaller power will have more concentrated reflect.
	vec4 reflectColor = (ratio) * skyC;
	float specCoef = min(0.1, pow(max(dot(viewDir, reflDir), 0.0), 64) * 3);

	//return skyC;
	vec4 c =  (1 - specCoef) * (heightColor+reflectColor) + vec4(spec) ;
	return c;// * skyC;

}



vec3 intersect(vec3 origin, vec3 direction, vec3 planeNormal, float D) {

	float calpha = dot(normalize(direction),normalize(-planeNormal));
	if (calpha > 0) {
	
		float k = (origin.y - D) * calpha;
		float x = origin.x + k * direction.x;
		float z = origin.z + k * direction.z;
		return vec3(x,D,z);
	}
	//caso o vetor refratado não intercete o fundo do mar
	else return vec3(0,-D*100,0);
}



void main() {
	
	vec2 slope = texture(htk, vec3(DataIn.texCoord, LAYER_DX_DZ_SX_SZ)).zw;
	vec3 wn = normalize(vec3( -slope.x,  1,  -slope.y));//
	
	vec4 color = computeOceanColor(wn);

#if (FOAM != NO_FOAM)
	vec4 foamV = texture(foam, DataIn.texCoord);
	float f = computeFoamFactor();
	outputF = color * (1-f) + foamV * f;
#else
	outputF = color;
#endif	
}
	
