require "paths"
local c = require 'sliceth.config'
local v = require 'sliceth.vectors'
local module = {}

local default = c.Config()
default.ExpType = "CBED"
default.nThreads = 1
default.Device = 0

local mc = c.ModelConfig()
mc.TiltBack = false
mc.CenterSample = false

mc.UseTDS = false
mc.TDSRuns = 1

mc.n = v.int3(1024,1024,10)
mc.d = v.float3(0.5,0.5,0.5)
mc.offset = v.float3(0,0,0)

mc.beamTiltX = 0.0
mc.beamTiltY = 0.0
mc.SourceDiameterAngstrom = 0.0
mc.BeamCurrentpA = 0.0

mc.PlotVrr = false
mc.periodicXY = false
mc.periodicZ = false
mc.DoZInterpolation = false
mc.UseQPotentialOffsets = false

mc.StructFactorType = "Rez"
mc.SliceCalcType = "SliceThickness"
mc.ResCalcType = "FILLRES"
mc.DisplaceType = "None"

mc.ratom = 5.0
mc.PotType = "CUDA2D"
mc.EnergykeV = 200
mc.ImagPot = 0.0

local oc = c.OutputConfig()
oc.LogLevel = 1
oc.SaveWaveIterations = 1
oc.SavePotential = false
oc.SaveProjectedPotential = false
oc.WriteLogFile = false
oc.SaveProbe = true
oc.SaveWaveAfterTransmit = false
oc.SaveWaveAfterTransform = false
oc.SaveWaveAfterPropagation = false
oc.SaveWaveAfterSlice = true
oc.SaveAtomicPotential = false
oc.ComputeFromProjectedPotential = false
oc.SaveAtomDeltas = false
oc.SaveAtomConv = false

oc.LogFileName = "./slicepp.log"
oc.SavePath = "save.h5"
oc.ConfigPath = ""

oc.PendelloesungPlot = false
oc.readPotential = false

local wc = c.WaveConfig()
wc.Cs = 0
wc.C5 = 0
wc.Cc = 0

wc.alpha = 5.0
wc.Defocus = 500
wc.Astigmatism = 0
wc.AstigmatismAngle = 0

wc.a_33 = 0
wc.a_31 = 0
wc.a_44 = 0
wc.a_42 = 0
wc.a_55 = 0
wc.a_53 = 0
wc.a_51 = 0
wc.a_66 = 0
wc.a_64 = 0
wc.a_62 = 0
wc.phi_33 = 0
wc.phi_31 = 0
wc.phi_44 = 0
wc.phi_42 = 0
wc.phi_55 = 0
wc.phi_53 = 0
wc.phi_51 = 0
wc.phi_66 = 0
wc.phi_64 = 0
wc.phi_62 = 0
wc.gaussScale = 1.0

wc.dI_I = 0
wc.dE_E = 0
wc.dV_V = 0

wc.AISaperture = 0
wc.tiltX = 0
wc.tiltY = 0

wc.Smooth = false
wc.Gaussian = false
wc.type = "Plane"
wc.n = v.int2(512,512)

local dc = c.DetectorConfig()

dc.mtfA = 1
dc.mtfB = 0
dc.mtfC = 5
dc.mtfD = 0
dc.DwellTimeMsec = 1 
dc.type = "Scintillator"
dc.n = v.int2(1024,1024)
dc.MaxElectronCounts = 5e4

local sc = c.StructureConfig()
sc.StructureFilename = "test"
sc.T_Kelvin = 300
sc.crystalTilt = v.float3(0,0,0)
sc.box = v.float3(0,0,0)
sc.zoneAxis = v.int3(1,0,0)
sc.nCells = v.int3(1,1,1)

sc.isBoxed = false
sc.rotateToZoneAxis = false

local scanc = c.ScanConfig()
scanc.xPos = 0;
scanc.yPos = 0;
scanc.xStep = 1;
scanc.yStep = 1;
scanc.nSteps = v.int2(5,5);
scanc.type = "Raster";
--scanc.yaml = ;

default.Model = mc
default.Output = oc
default.Structure = sc 
default.Wave = wc
default.Detector = dc
default.Scan = scanc

local Builder = {}
local Builder_mt = {}
function Builder_mt:__call()
--  print(module.DefaultConfig)
  self.c = module.DefaultConfig:copy()
  return self
end

function Builder:File(filename)
  local base = paths.basename(filename) 
  local dir = paths.dirname(filename)
  if(paths.filep(filename .. '.cif')) then
    self.c.Structure.StructureFilename = filename .. '.cif'
  else 
    self.c.Structure.StructureFilename = filename .. '.gbm'
  end
  self.c.Output.SavePath = paths.concat(dir,base) .. '.h5'
  self.c.Output.LogFileName = paths.concat(dir,base) .. '.log'
  self.c.Output.ConfigPath = dir
  return self
end

function Builder:Wave(w)
  self.c.Wave = w
  return self
end
function Builder:Model(w)
  self.c.Model = w
  return self
end
function Builder:Output(w)
  self.c.Output = w
  return self
end
function Builder:Structure(w)
  self.c.Structure = w
  return self
end
function Builder:Detector(w)
  self.c.Detector = w
  return self
end

function Builder:getConfig()
  return self.c
end

setmetatable(Builder,Builder_mt)

module.DefaultConfig = default
module.Builder = Builder
return module