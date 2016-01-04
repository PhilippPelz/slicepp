local s = require 'sliceth'

local b = s.configbuilder.Builder():File('/home/philipp/projects/slicepp/Examples/superstructures/goldball')

local c = b:getConfig()

s.runner.run(c)


