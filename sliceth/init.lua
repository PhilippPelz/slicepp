require 'torch'
require 'ffi'
local sliceth = {}

sliceth.runner = require 'sliceth.runner'
sliceth.config = require 'sliceth.config'
sliceth.detectors = require 'sliceth.detectors'
sliceth.Config = sliceth.config.Config

function sliceth.float3(one,two,three)
  return ffi.new("float[3]",{one,two,three})
end 

function sliceth.int3(one,two,three)
  return ffi.new("int[3]",{one,two,three})
end

function sliceth.int2(one,two)
  return ffi.new("int[2]",{one,two})
end

return sliceth
