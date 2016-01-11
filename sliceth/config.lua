local C = require 'sliceth.sliceFFI'
local ffi = require 'ffi'

local c = {}

local StructureConfig_mt = {
  __call = function ()
              local sc = ffi.new("StructureConfig")
--              ffi.gc(sc,C.StructureConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("StructureConfig")
              ffi.copy(sc, self, ffi.sizeof("StructureConfig"))
              return sc
           end
}
local ModelConfig_mt = {
  __call = function ()
              local sc = ffi.new("ModelConfig")
--              ffi.gc(sc,C.ModelConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("ModelConfig")
              ffi.copy(sc, self, ffi.sizeof("ModelConfig"))
              return sc
           end
}
local WaveConfig_mt = {
  __call = function ()
              local sc = ffi.new("WaveConfig")
--              ffi.gc(sc,C.WaveConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("WaveConfig")
              ffi.copy(sc, self, ffi.sizeof("WaveConfig"))
              return sc
           end
}
local OutputConfig_mt = {
  __call = function ()
              local sc = ffi.new("OutputConfig")
--              ffi.gc(sc,C.OutputConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("OutputConfig")
              ffi.copy(sc, self, ffi.sizeof("OutputConfig"))
              return sc
           end
}
local DetectorConfig_mt = {
  __call = function ()
              local sc = ffi.new("DetectorConfig")
--              ffi.gc(sc,C.DetectorConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("DetectorConfig")
              ffi.copy(sc, self, ffi.sizeof("DetectorConfig"))
              return sc
           end
}
local ScanConfig_mt = {
  __call = function ()
              local sc = ffi.new("sScanConfig")
--              ffi.gc(sc,C.sScanConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("sScanConfig")
              ffi.copy(sc, self, ffi.sizeof("sScanConfig"))
              return sc
           end
}

c.StructureConfig = ffi.metatype("StructureConfig",StructureConfig_mt)
c.ModelConfig = ffi.metatype("ModelConfig",ModelConfig_mt)
c.WaveConfig = ffi.metatype("WaveConfig",WaveConfig_mt)
c.OutputConfig = ffi.metatype("OutputConfig",OutputConfig_mt)
c.DetectorConfig = ffi.metatype("DetectorConfig",DetectorConfig_mt)
c.ScanConfig = ffi.metatype("sScanConfig",ScanConfig_mt)

local Config_mt = {
  __call = function()
    local conf = ffi.new("sConfig")
    local sc = c.StructureConfig()
    local mc = c.ModelConfig()
    local wc = c.WaveConfig() 
    local oc = c.OutputConfig()
    local dc = c.DetectorConfig()
    local scanc = c.ScanConfig()
    
    conf.Structure = sc
    conf.Model = mc
    conf.Output = oc
    conf.Wave = wc
    conf.Scan = scanc
    conf.Detector = dc
    
--    ffi.gc(conf,C.sConfig_delete)
    return conf
  end,
  __index = {
    copy =  function (self)
              local sc = ffi.new("sConfig")
              ffi.copy(sc, self, ffi.sizeof("sConfig"))
              return sc
            end  
  }
}
c.Config = ffi.metatype("sConfig",Config_mt)

return c

