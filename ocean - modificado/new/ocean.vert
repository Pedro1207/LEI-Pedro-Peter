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



uniform sampler2DArray htk;
uniform int L;
uniform uint heightMapSize;
uniform float choppyFactor;

uniform mat4 m_pvm;
uniform mat3 m_normal;
uniform mat4 m_view;
uniform vec4 l_dir;

in vec4 position;
in vec2 texCoord0;

out Data {
	vec3 l_dir;
	vec3 pos;
	vec2 texCoord;
	vec3 normal;
} DataOut;

void main() {

	//vec2 tc = texCoord0;//position.xz * 0.5 + 0.5 ;
	vec2 tc = (vec2(position.xz) + vec2(float(heightMapSize) * 0.5)) / heightMapSize;
	
	DataOut.l_dir = vec3(normalize( (m_view * l_dir)));
	DataOut.l_dir = vec3(normalize( (l_dir)));

	float scaleFactor = float(heightMapSize)/L;
	vec4 dy = texture(htk, vec3(tc, LAYER_Y_JXY_JXX_JYY)) * scaleFactor;
	float h = dy.x  ;
	vec4 ds = texture(htk, vec3(tc, LAYER_DX_DZ_SX_SZ))  * scaleFactor;
	vec2 dxz = ds.xy * choppyFactor;// * scaleFactor;
	
	DataOut.texCoord = tc;

	vec2 slope = ds.zw;
	slope = vec2(slope.x/(1 + choppyFactor * dy.z), slope.y/(1 + choppyFactor * dy.w));
	DataOut.normal = normalize(vec3(-slope.x, 1, -slope.y));

	//DataOut.pos = vec3(position.x * L * 0.5f - dxz.x, h, position.z * L*0.5f - dxz.y);
	DataOut.pos = vec3(position.x  - dxz.x, h, position.z - dxz.y);
	gl_Position = m_pvm * vec4(DataOut.pos, 1);
}
