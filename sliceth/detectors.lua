
local c = require 'sliceth.config'

local module = {}

local F216 = c.DetectorConfig()

dc.mtfA = 1
dc.mtfB = 0
dc.mtfC = 5
dc.mtfD = 0
dc.DwellTimeMsec = 1 
dc.type = "Scintillator"
dc.n = s.int2(2048,2048)
dc.MaxElectronCounts = 5e4

module.TvipsF216 = F216
return module