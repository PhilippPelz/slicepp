vectors = {}
local ffi = require 'ffi'
function vectors.float3(one,two,three)
  return ffi.new("float[3]",{one,two,three})
end 

function vectors.int3(one,two,three)
  return ffi.new("int[3]",{one,two,three})
end

function vectors.int2(one,two)
  return ffi.new("int[2]",{one,two})
end
return vectors