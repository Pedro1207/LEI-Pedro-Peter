
-- sets common variables prior to starting both forward and inverse FFT
function resetCommon()

	local stage = {-1}
	
	-- we set the stage to -1 so it gets to zero on the first iteration
	setAttr("RENDERER", "CURRENT", "current_stage", 0, stage)

	
end




-- sets the variables prior to starting the inverse FFT
function resetInverseFFT()

	local pingpong = {1}
	setAttr("RENDERER", "CURRENT", "pingpong", 0, pingpong)
	
	resetCommon()
end


-- sets the variables prior to starting the inverse FFT
function resetInverseVertical()

	resetCommon()
end



-- controls the stage loop cycle
function fft()

	local pingpong = {}
	local stage = {}
	local log_2 = {}

	
	-- control the stage loop and pingpong variable
	getAttr("RENDERER", "CURRENT", "log_width", 0, log_2)
	getAttr("RENDERER", "CURRENT", "current_stage", 0, stage)
	stage[1] = stage[1] + 1
	
	if (stage[1] < log_2[1]) then
		-- flip the ping pong variable
		getAttr("RENDERER", "CURRENT", "pingpong", 0, pingpong)
		pingpong[1] = (pingpong[1] + 1) % 2
		setAttr("RENDERER", "CURRENT", "pingpong", 0, pingpong)
		setAttr("RENDERER", "CURRENT", "current_stage", 0, stage)
		return true
	else
		return false
	end
end

-- controls the mipmap loop cycle
testMipmap = function() 

	local f = {0}
	local g = {0}
	local h = {0}
	local k = {0}
	getAttr("RENDERER", "CURRENT", "log_width", 0, h);
	getAttr("RENDERER", "CURRENT", "currentLevel",0,f)
	
	if f[1] == h[1] then
		f[1] = 0
		setAttr("RENDERER", "CURRENT", "currentLevel",0,f)
		return false
	else
		f[1] = f[1] + 1
		k[1] = h[1] - f[1];
		if k[1] == 0 then
			g[1] = 1;
		else g[1] = 2;	
			for i=2,k[1] do
				g[1] = g[1] * 2
			end	
		end
		
		setAttr("RENDERER", "CURRENT", "currentLevel",0,f)
		setAttr("PASS", "CURRENT", "DIM_X", 0, g);
		setAttr("PASS", "CURRENT", "DIM_Y", 0, g);
		
		setAttr("IMAGE_TEXTURE", "Ocean::mipmap", "LEVEL", 1, f);
		return true
	end
end
