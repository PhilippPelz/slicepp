require 'torch'
require 'ffi'
local sliceth = {}

sliceth.runner = require 'sliceth.runner'
sliceth.config = require 'sliceth.config'
sliceth.detectors = require 'sliceth.detectors'
sliceth.Config = sliceth.config.Config
sliceth.configbuilder = require 'sliceth.configbuilder'
sliceth.vectors = require 'sliceth.vectors'
sliceth.v = require 'sliceth.vectors'

return sliceth
