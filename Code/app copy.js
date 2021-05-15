var vertexShaderText = 
`#version 300 es
precision highp float;

in vec3 vertPosition;
in vec3 vertColor;
out vec3 fragColor;
uniform mat4 mWorld;
uniform mat4 mView;
uniform mat4 mProj;

void main()
{
  fragColor = vertColor;
  gl_Position = mProj * mView * mWorld * vec4(vertPosition, 1.0);
}`

var fragmentShaderText =
`#version 300 es
precision highp float;

in vec3 fragColor;
out vec4 fragColorOut;

void main()
{
	fragColorOut = vec4(fragColor, 1.0);
}
`


var InitDemo = function () {
	console.log('This is working');

	var canvas = document.getElementById('game-surface');
	var gl = canvas.getContext('webgl2');

	if (!gl) {
		alert('Your browser does not support WebGL');
	}

	gl.clearColor(0.75, 0.85, 0.8, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.frontFace(gl.CCW);
	gl.cullFace(gl.BACK);

	//
	// Create shaders
	// 
	var vertexShader = gl.createShader(gl.VERTEX_SHADER);
	var fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);

	gl.shaderSource(vertexShader, vertexShaderText);
	gl.shaderSource(fragmentShader, fragmentShaderText);

	gl.compileShader(vertexShader);
	if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS)) {
		console.error('ERROR compiling vertex shader!', gl.getShaderInfoLog(vertexShader));
		return;
	}

	gl.compileShader(fragmentShader);
	if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS)) {
		console.error('ERROR compiling fragment shader!', gl.getShaderInfoLog(fragmentShader));
		return;
	}

	var program = gl.createProgram();
	gl.attachShader(program, vertexShader);
	gl.attachShader(program, fragmentShader);
	gl.linkProgram(program);
	if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
		console.error('ERROR linking program!', gl.getProgramInfoLog(program));
		return;
	}
	gl.validateProgram(program);
	if (!gl.getProgramParameter(program, gl.VALIDATE_STATUS)) {
		console.error('ERROR validating program!', gl.getProgramInfoLog(program));
		return;
	}

	//
	// Create buffer
	//
	var boxVertices = 
	[ // X, Y, Z           R, G, B
		// Top
		-1.0, 1.0, -1.0,   0.5, 0.5, 0.5,
		-1.0, 1.0, 1.0,    0.5, 0.5, 0.5,
		1.0, 1.0, 1.0,     0.5, 0.5, 0.5,
		1.0, 1.0, -1.0,    0.5, 0.5, 0.5,

		// Left
		-1.0, 1.0, 1.0,    0.75, 0.25, 0.5,
		-1.0, -1.0, 1.0,   0.75, 0.25, 0.5,
		-1.0, -1.0, -1.0,  0.75, 0.25, 0.5,
		-1.0, 1.0, -1.0,   0.75, 0.25, 0.5,

		// Right
		1.0, 1.0, 1.0,    0.25, 0.25, 0.75,
		1.0, -1.0, 1.0,   0.25, 0.25, 0.75,
		1.0, -1.0, -1.0,  0.25, 0.25, 0.75,
		1.0, 1.0, -1.0,   0.25, 0.25, 0.75,

		// Front
		1.0, 1.0, 1.0,    1.0, 0.0, 0.15,
		1.0, -1.0, 1.0,    1.0, 0.0, 0.15,
		-1.0, -1.0, 1.0,    1.0, 0.0, 0.15,
		-1.0, 1.0, 1.0,    1.0, 0.0, 0.15,

		// Back
		1.0, 1.0, -1.0,    0.0, 1.0, 0.15,
		1.0, -1.0, -1.0,    0.0, 1.0, 0.15,
		-1.0, -1.0, -1.0,    0.0, 1.0, 0.15,
		-1.0, 1.0, -1.0,    0.0, 1.0, 0.15,

		// Bottom
		-1.0, -1.0, -1.0,   0.5, 0.5, 1.0,
		-1.0, -1.0, 1.0,    0.5, 0.5, 1.0,
		1.0, -1.0, 1.0,     0.5, 0.5, 1.0,
		1.0, -1.0, -1.0,    0.5, 0.5, 1.0,
	];

	var boxIndices =
	[
		// Top
		0, 1, 2,
		0, 2, 3,

		// Left
		5, 4, 6,
		6, 4, 7,

		// Right
		8, 9, 10,
		8, 10, 11,

		// Front
		13, 12, 14,
		15, 14, 12,

		// Back
		16, 17, 18,
		16, 18, 19,

		// Bottom
		21, 20, 22,
		22, 20, 23
	];

	var boxVertexBufferObject = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER, boxVertexBufferObject);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(boxVertices), gl.STATIC_DRAW);

	var boxIndexBufferObject = gl.createBuffer();
	gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, boxIndexBufferObject);
	gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(boxIndices), gl.STATIC_DRAW);

	var positionAttribLocation = gl.getAttribLocation(program, 'vertPosition');
	var colorAttribLocation = gl.getAttribLocation(program, 'vertColor');
	gl.vertexAttribPointer(
		positionAttribLocation, // Attribute location
		3, // Number of elements per attribute
		gl.FLOAT, // Type of elements
		gl.FALSE,
		6 * Float32Array.BYTES_PER_ELEMENT, // Size of an individual vertex
		0 // Offset from the beginning of a single vertex to this attribute
	);
	gl.vertexAttribPointer(
		colorAttribLocation, // Attribute location
		3, // Number of elements per attribute
		gl.FLOAT, // Type of elements
		gl.FALSE,
		6 * Float32Array.BYTES_PER_ELEMENT, // Size of an individual vertex
		3 * Float32Array.BYTES_PER_ELEMENT // Offset from the beginning of a single vertex to this attribute
	);

	gl.enableVertexAttribArray(positionAttribLocation);
	gl.enableVertexAttribArray(colorAttribLocation);

	// Tell OpenGL state machine which program should be active.
	gl.useProgram(program);

	var matWorldUniformLocation = gl.getUniformLocation(program, 'mWorld');
	var matViewUniformLocation = gl.getUniformLocation(program, 'mView');
	var matProjUniformLocation = gl.getUniformLocation(program, 'mProj');

	var worldMatrix = new Float32Array(16);
	var viewMatrix = new Float32Array(16);
	var projMatrix = new Float32Array(16);
	glMatrix.mat4.identity(worldMatrix);
	glMatrix.mat4.lookAt(viewMatrix, [0, 0, -8], [0, 0, 0], [0, 1, 0]);
	glMatrix.mat4.perspective(projMatrix, glMatrix.glMatrix.toRadian(45), canvas.clientWidth / canvas.clientHeight, 0.1, 1000.0);

	gl.uniformMatrix4fv(matWorldUniformLocation, gl.FALSE, worldMatrix);
	gl.uniformMatrix4fv(matViewUniformLocation, gl.FALSE, viewMatrix);
	gl.uniformMatrix4fv(matProjUniformLocation, gl.FALSE, projMatrix);
	

	/*
	// FFT texture dim -->
	var fftTextureDim = gl.getUniformLocation(program, 'fftTexDim');
	gl.uniform1i(fftTextureDim,512);

	// log of width - not set by lua script -->	
	var butterflyStage = gl.getUniformLocation(program, 'butterStage');
	gl.uniform1i(butterflyStage,9);

	//butterfly stage -->
	var logOfWidth = gl.getUniformLocation(program, 'logOfWidth');
	gl.uniform1i(logOfWidth,0);

	// pingpong status - set not by lua script 
	var pingpongStatus = gl.getUniformLocation(program, 'pingpongStatus');
	gl.uniform1i(pingpongStatus,0);

	// attribute for mipampping: current level - not set by lua script -->
	var currentLevelMipmap = gl.getUniformLocation(program, 'currentLevel');
	gl.uniform1i(currentLevelMipmap,0);

	// random distribution -->
	var randomDistribution = gl.getUniformLocation(program, 'randomDistribution');
	gl.uniform1i(randomDistribution,1);
	
	// Ocean attributes -->
	var oceanDepth = gl.getUniformLocation(program , 'oceanDepth');
	gl.uniform1f(oceanDepth,20);
	var L = gl.getUniformLocation(program, 'L');
	gl.uniform1i(L,512);
	var spectrum = gl.getUniformLocation(program, 'spectrum');
	gl.uniform1i(spectrum,1);
	var dispersion = gl.getUniformLocation(program, 'dispersion');
	gl.uniform1i(dispersion,0);
	var dirSpreading = gl.getUniformLocation(program, 'dirSpreading');
	gl.uniform1i(dirSpreading,0);
	var spectrumScale = gl.getUniformLocation(program, 'spectrumScale');
	gl.uniform1f(spectrumScale,1);
	var choppyFactor = gl.getUniformLocation(program, 'choppyFactor');
	gl.uniform1f(choppyFactor,1);

	//-- Spectrum parameters -->
	var fetch = gl.getUniformLocation(program, 'fetch');
	gl.uniform1f(fetch,1200000.0);
	var windDir = gl.getUniformLocation(program, 'windDir');
	gl.uniform2f(windDir,1,0);
	var windSpeed = gl.getUniformLocation(program, 'windSpeed');
	gl.uniform1f(windSpeed, 10.0);
	
	//Directional spreading parameters -->
	var propagate = gl.getUniformLocation(program, 'propagate');
	gl.uniform1i(propagate, 1);
	var swell = gl.getUniformLocation(program, 'swell');
	gl.uniform1f(swell, 0.0);

	//Occhi and Bretschneider -->
	var Hs = gl.getUniformLocation(program, 'Hs');
	gl.uniform1f(Hs, 10.0);
	
	//JONSWAP parameters
	var JONSWAP_gamma = gl.getUniformLocation(program, 'JONSWAP_gamma');
	gl.uniform1f(JONSWAP_gamma, 3.3);
	var JONSWAP_sigmaA = gl.getUniformLocation(program, 'JONSWAP_sigmaA');
	gl.uniform1f(JONSWAP_sigmaA, 0.07);
	var JONSWAP_sigmaB = gl.getUniformLocation(program, 'JONSWAP_sigmaB');
	gl.uniform1f(JONSWAP_sigmaB, 0.09);

	//Bretschneider parameters
	var Bretschneider_wm = gl.getUniformLocation(program, 'Bretschneider_wm');
	gl.uniform1f(Bretschneider_wm, 0.0);
			
	//Ochi parameters
	var Ochi_lambda1 = gl.getUniformLocation(program, 'Ochi_lambda1');
	gl.uniform1f(Ochi_lambda1, 3.0);
	var Ochi_lambda2 = gl.getUniformLocation(program, 'Ochi_lambda2');
	gl.uniform1f(Ochi_lambda2, 0.0);
	var Ochi_wm1 = gl.getUniformLocation(program, 'Ochi_wm1');
	gl.uniform1f(Ochi_wm1, 0.0);
	var Ochi_wm2 = gl.getUniformLocation(program, 'Ochi_wm2');
	gl.uniform1f(Ochi_wm2, 0.0);
	var Ochi_Hs1 = gl.getUniformLocation(program, 'Ochi_Hs1');
	gl.uniform1f(Ochi_Hs1, 0.0);
	var Ochi_Hs2 = gl.getUniformLocation(program, 'Ochi_Hs2');
	gl.uniform1f(Ochi_Hs2, 0.0);

	//Color information
	var oceanType = gl.getUniformLocation(program, 'oceanType');
	gl.uniform1i(oceanType, 1);
	var oceanTrans = gl.getUniformLocation(program, 'oceanTrans');
	gl.uniform3f(oceanTrans, 98.2, 95.8, 57.0);
	var oceanFloorColor = gl.getUniformLocation(program, 'oceanFloorColor');
	gl.uniform3f(oceanFloorColor, 0.956, 0.925, 0.925)
	var debug = gl.getUniformLocation(program, 'debug');
	gl.uniform3f(debug, 0, 0, 0)
	*/

	var xRotationMatrix = new Float32Array(16);
	var yRotationMatrix = new Float32Array(16);

	var angleDelta = 0.05;
	//set up sliders
    var slider = document.getElementById("slide");
	slider.oninput = function(){
		angleDelta = slider.value / 100;
	}

	//
	// Main render loop
	//
	var identityMatrix = new Float32Array(16);
	glMatrix.mat4.identity(identityMatrix);
	var angle = 0;
	var loop = function () {
		angle = angle + angleDelta;
		glMatrix.mat4.rotate(yRotationMatrix, identityMatrix, angle, [0, 1, 0]);
		glMatrix.mat4.rotate(xRotationMatrix, identityMatrix, angle / 4, [1, 0, 0]);
		glMatrix.mat4.mul(worldMatrix, yRotationMatrix, xRotationMatrix);
		gl.uniformMatrix4fv(matWorldUniformLocation, gl.FALSE, worldMatrix);

		gl.clearColor(0.75, 0.85, 0.8, 1.0);
		gl.clear(gl.DEPTH_BUFFER_BIT | gl.COLOR_BUFFER_BIT);
		gl.drawElements(gl.TRIANGLES, boxIndices.length, gl.UNSIGNED_SHORT, 0);

		requestAnimationFrame(loop);
	};
	requestAnimationFrame(loop);
};