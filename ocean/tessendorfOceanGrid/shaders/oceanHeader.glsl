
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


