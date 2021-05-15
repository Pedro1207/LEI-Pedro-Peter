
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