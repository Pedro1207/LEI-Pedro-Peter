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


in vec2 texCoordV;

// ping pong textures
layout (binding = 0, rgba32f) uniform image2DArray pingpong0;
layout (binding = 1, rgba32f) uniform image2DArray pingpong1;
layout (binding = 2, rgba32f) uniform image2DArray test;

uniform int pingpong;
uniform int log_width;
uniform int stage;


vec2 complexMult(vec2 v0, vec2 v1) {
	return vec2(v0.x * v1.x - v0.y * v1.y,
				v0.x * v1.y + v0.y * v1.x);
}

 vec4 complexMultTwice(vec2 v0, vec4 v1) {
	return vec4(v0.x * v1.x - v0.y * v1.y,
				v0.x * v1.y + v0.y * v1.x,
				v0.x * v1.z - v0.y * v1.w,
				v0.x * v1.w + v0.y * v1.z);
}


vec2 w(int k, int nn) {
	float it =  2 * k * M_PI / nn;
	return vec2(cos(it), sin(it));
}


void main() {

	vec2 aux, aux1, aux2, raux, rs1, rs2;

	uint line = uint(texCoordV.x * 512);
	int column = int(texCoordV.y * 512);
	
	int iter = int(pow(2,log_width-1));

	int halfGroupSize = int(pow(2, stage));
	int groupSize = 2 * halfGroupSize;
	int k = column % halfGroupSize;
	int group = column / halfGroupSize;
	int shift = int(pow(2, stage));
	int groupShift = shift * 2;

	int index = k + group * groupShift;

	vec2 ww;
	vec4 elemk, elemks;
	vec4 elems4, elemss4, elemxz, elemxzs, elemj, elemjs;
	
	// alternate between textures
	if (pingpong == 0) {
		// when stage = 0 use bit reverse indices
		if (stage == 0) {
			uint br = bitfieldReverse(uint(index));
			br = bitfieldExtract(br, 32 - log_width, log_width);
			uint brs = bitfieldReverse(uint(index + shift));
			brs = bitfieldExtract(brs, 32 - log_width, log_width);
			elemk = imageLoad(pingpong0, ivec3(line, br, LAYER_Y_JXY_JXX_JYY));
			elemks = imageLoad(pingpong0, ivec3(line, brs, LAYER_Y_JXY_JXX_JYY));
			elemxz = imageLoad(pingpong0, ivec3(line, br, LAYER_DX_DZ_SX_SZ));
			elemxzs = imageLoad(pingpong0, ivec3(line, brs, LAYER_DX_DZ_SX_SZ));
		}
		else {
			elemk = imageLoad(pingpong0, ivec3(line, index, LAYER_Y_JXY_JXX_JYY));
			elemks = imageLoad(pingpong0, ivec3(line, index + shift, LAYER_Y_JXY_JXX_JYY));
			elemxz = imageLoad(pingpong0, ivec3(line, index, LAYER_DX_DZ_SX_SZ));
			elemxzs = imageLoad(pingpong0, ivec3(line, index + shift, LAYER_DX_DZ_SX_SZ));
		}
		
		// compute the twiddle factor
		ww = w(k, groupShift);
			
		// write the outputs
		vec4 cm = complexMultTwice(ww, elemks);
		imageStore(pingpong1, ivec3(line, index, LAYER_Y_JXY_JXX_JYY), elemk + cm);
		imageStore(pingpong1, ivec3(line, index + shift, LAYER_Y_JXY_JXX_JYY), elemk - cm);
		imageStore(test, ivec3(line, index, LAYER_Y_JXY_JXX_JYY), elemk + cm);
		imageStore(test, ivec3(line, index + shift, LAYER_Y_JXY_JXX_JYY), elemk - cm);

		cm = complexMultTwice(ww,elemxzs);
		imageStore(pingpong1, ivec3(line, index, LAYER_DX_DZ_SX_SZ), elemxz + cm);
		imageStore(pingpong1, ivec3(line, index + shift, LAYER_DX_DZ_SX_SZ), elemxz - cm);
			
	}
	else {
		if (stage == 0) {
			uint br = bitfieldReverse(uint(index));
			br = bitfieldExtract(br, 32 - log_width, log_width);
			uint brs = bitfieldReverse(uint(index + shift));
			brs = bitfieldExtract(brs, 32 - log_width, log_width);
			elemk = imageLoad(pingpong1, ivec3(line, br, LAYER_Y_JXY_JXX_JYY));
			elemks = imageLoad(pingpong1, ivec3(line, brs, LAYER_Y_JXY_JXX_JYY));
			elemxz = imageLoad(pingpong1, ivec3(line, br, LAYER_DX_DZ_SX_SZ));
			elemxzs = imageLoad(pingpong1, ivec3(line, brs, LAYER_DX_DZ_SX_SZ));
		}
		else {	
			elemk = imageLoad(pingpong1, ivec3(line, index, LAYER_Y_JXY_JXX_JYY));
			elemks = imageLoad(pingpong1, ivec3(line, index + shift, LAYER_Y_JXY_JXX_JYY));
			elemxz = imageLoad(pingpong1, ivec3(line, index, LAYER_DX_DZ_SX_SZ));
			elemxzs = imageLoad(pingpong1, ivec3(line, index + shift, LAYER_DX_DZ_SX_SZ));
		}
		ww = w(k, groupShift);
			
		vec4 cm = complexMultTwice(ww, elemks);
		imageStore(pingpong0, ivec3(line, index, LAYER_Y_JXY_JXX_JYY), elemk + cm);
		imageStore(pingpong0, ivec3(line, index + shift, LAYER_Y_JXY_JXX_JYY), elemk - cm);
		imageStore(test, ivec3(line, index, LAYER_Y_JXY_JXX_JYY), elemk + cm);
		imageStore(test, ivec3(line, index + shift, LAYER_Y_JXY_JXX_JYY), elemk - cm);

		cm = complexMultTwice(ww, elemxzs);
		imageStore(pingpong0, ivec3(line, index, LAYER_DX_DZ_SX_SZ), elemxz + cm);
		imageStore(pingpong0, ivec3(line, index + shift, LAYER_DX_DZ_SX_SZ), elemxz - cm);
		
	}
    discard;
}