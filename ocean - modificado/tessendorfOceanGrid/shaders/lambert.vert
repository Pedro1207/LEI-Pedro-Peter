
uniform sampler2DArray htk;
uniform int L;
uniform float choppyFactor, windSpeed;
uniform uint heightMapSize;

uniform mat4 m_pvm;
uniform mat3 m_normal;

in vec4 position;
in vec3 normal;

out vec3 normalV;

void main() {

	vec2 disp = vec2(50,50);
	vec2 tc = (vec2(disp) + vec2(float(heightMapSize) * 0.5)) / heightMapSize;
	float scaleFactor = float(heightMapSize)/L;
	vec4 dy = texture(htk, vec3(tc, LAYER_Y_JXY_JXX_JYY)) * scaleFactor;
	float h = dy.x  ;
	vec4 ds = texture(htk, vec3(tc, LAYER_DX_DZ_SX_SZ))  * scaleFactor;
	vec2 dxz = ds.xy * choppyFactor;// * scaleFactor;


	vec4 pos = position*scaleFactor;
	pos.xz += disp - dxz.xy  ;//( 1+choppyFactor/(1+(exp(-windSpeed+20)/5)));
	
	pos.y += h ;
	pos.xyz *= scaleFactor;
	normalV = normalize(m_normal * normal);
	
	gl_Position = m_pvm * pos;
}

