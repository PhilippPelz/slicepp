local C = require 'sliceth.sliceFFI'
local ffi = require 'ffi'

local c = {}

local StructureConfig_mt = {
  __call = function ()
              local sc = ffi.new("StructureConfig")
              ffi.gc(sc,C.StructureConfig_delete)
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
              ffi.gc(sc,C.StructureConfig_delete)
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
              ffi.gc(sc,C.StructureConfig_delete)
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
              ffi.gc(sc,C.StructureConfig_delete)
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
              ffi.gc(sc,C.StructureConfig_delete)
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
              local sc = ffi.new("ScanConfig")
              ffi.gc(sc,C.StructureConfig_delete)
              return sc
           end,                
  copy = function (self)
              local sc = ffi.new("ScanConfig")
              ffi.copy(sc, self, ffi.sizeof("ScanConfig"))
              return sc
           end
}

c.StructureConfig = ffi.metatype("StructureConfig",StructureConfig_mt)
c.ModelConfig = ffi.metatype("ModelConfig",ModelConfig_mt)
c.WaveConfig = ffi.metatype("WaveConfig",WaveConfig_mt)
c.OutputConfig = ffi.metatype("OutputConfig",OutputConfig_mt)
c.DetectorConfig = ffi.metatype("DetectorConfig",DetectorConfig_mt)
c.ScanConfig = ffi.metatype("ScanConfig",ScanConfig_mt)

local Config_mt = {
  __call = function()
    local c = ffi.new("c_Config")
    local sc = ffi.new("StructureConfig")
    local mc = ffi.new("ModelConfig")
    local wc = ffi.new("WaveConfig")
    local oc = ffi.new("OutputConfig")
    local dc = ffi.new("DetectorConfig")
    local scanc = ffi.new("ScanConfig")
    
    c.Structure = sc
    c.Model = mc
    c.Output = oc
    c.Wave = wc
    c.Scan = scanc
    c.Detector = dc
    
    ffi.gc(c,C.c_Config_delete)
    ffi.gc(sc,C.StructureConfig_delete)
    ffi.gc(mc,C.ModelConfig_delete)
    ffi.gc(wc,C.WaveConfig_delete)
    ffi.gc(oc,C.OutputConfig_delete)
    ffi.gc(dc,C.DetectorConfig_delete)
    ffi.gc(scanc,C.ScanConfig_delete)
    
    return c  
  end,
  __index = {
    copy =  function (self)
              local sc = ffi.new("c_Config")
              ffi.copy(sc, self, ffi.sizeof("c_Config"))
              return sc
            end  
  }
}
c.Config = ffi.metatype("c_Config",Config_mt)

return c

