local C = require 'sliceth.sliceFFI'
local ffi = require 'ffi'

local c = {}

local StructureConfig_mt = {}
local ModelConfig_mt = {}
local WaveConfig_mt = {}
local OutputConfig_mt = {}
local DetectorConfig_mt = {}
local ScanConfig_mt = {}

c.StructureConfig = ffi.metatype("StructureConfig",StructureConfig_mt)
c.ModelConfig = ffi.metatype("ModelConfig",ModelConfig_mt)
c.WaveConfig = ffi.metatype("WaveConfig",WaveConfig_mt)
c.OutputConfig = ffi.metatype("OutputConfig",OutputConfig_mt)
c.DetectorConfig = ffi.metatype("DetectorConfig",DetectorConfig_mt)
c.ScanConfig = ffi.metatype("ScanConfig",ScanConfig_mt)

local Config_mt = {}

function Config_mt.__call()
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
end

--function Config_mt.__newindex(self,index,value)
--  print(a)
--  print(index)
--  print(value)  
--  return self
--end

c.Config = ffi.metatype("c_Config",Config_mt)

return c

