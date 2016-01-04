local v = require 'sliceth.vectors'
local c = require 'sliceth.config'

local module = {}

local F216 = c.DetectorConfig()

F216.mtfA = 1
F216.mtfB = 0
F216.mtfC = 5
F216.mtfD = 0
F216.DwellTimeMsec = 1 
F216.type = "Scintillator"
F216.n = v.int2(2048,2048)
F216.MaxElectronCounts = 5e4

module.TvipsF216 = F216
return module