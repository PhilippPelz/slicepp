local s = require 'sliceth'

local c = s.Config()
c.ExperimentType = "CBED"

local sc = c.Structure 
sc.structureFilename = "test"
sc.T_Kelvin = 300
sc.crystalTilt = s.float3(0,0,0)
sc.box = s.float3(0,0,0)
sc.zoneAxis = s.int3(1,0,0)
sc.nCells = s.int3(1,1,1)

sc.isBoxed = false
sc.rotateToZoneAxis = false
  
local mc = c.Model 


  
local oc = c.Output

  
local wc = c.Wave 


local scanc = c.Scan 



local dc = c.Detector



slice.runner.run(c)

print(torch)
print(c)


--{
--   nthreads =  1 ,
--   scan = {
--     scanType =  1 ,
--     x Start Position =  0 ,
--     yStep =  128 ,
--     y Start Position =  0 ,
--     xStep =  128 ,
--     scanx =  256 ,
--     scany =  256 
--  },
--   wave = {

--  },
--   beam = {
--     sourceDiameterAngstrom =  0 ,
--     dwellTimeMsec =  1.6021773e-4 ,
--     beamCurrentpA =  1 ,
--     energy_keV =  200.000000 
--  },
--   mode =  1 ,
--   output = {
--     showProbe =  false ,
--     writeLogFile =  false ,
--     loglevel =  0 ,
--     SaveAtomDeltas =  true ,
--     SaveWaveAfterTransform =  false ,
--     SaveWaveAfterNSlices =  1 ,
--     logFileName =  slicelog.log ,
--     pendelloesungPlot =  false ,
--     ComputeFromProjectedPotential =  false ,
--     SaveAtomicPotential =  true ,
--     savePath =  gold.h5 ,
--     readPotential =  false ,
--     SaveWaveAfterPropagation =  false ,
--     SaveAtomConv =  true ,
--     SaveWaveAfterTransmit =  false ,
--     folder =  CBED ,
--     saveProjectedPotential =  true ,
--     saveProbe =  true ,
--     SaveWaveAfterSlice =  true ,
--     savePotential =  true 
--  },
--   model = {
--     slices =  80 ,
--     tiltBack =  false ,
--     centerSample =  false ,
--     3D =  false ,
--     periodicXY =  true ,
--     atomRadiusAngstrom =  5.0 ,
--     nx =  320 ,
--     ny =  320 ,
--     UseQPotentialOffsets =  false ,
--     structureFactors =  1 ,
--     CUDAOnTheFly =  false ,
--     type =  cuda ,
--     resolutionYAngstrom =  0.25 ,
--     plotVr_r =  false ,
--     beamTiltY =  0.000000 ,
--     beamTiltX =  0.000000 ,
--     resolutionCalculation =  2 ,
--     DoZ erpolation =  false ,
--     sliceThicknessCalculation =  2 ,
--     tds =  false ,
--     periodicZ =  false ,
--     FFT =  false ,
--     sliceThicknessAngstrom =  2.039125 ,
--     centerSlices =  false ,
--     resolutionXAngstrom =  0.25 ,
--     displacementType =  3 ,
--     tdsRuns =  1 
--  },
--   detector = {
--     mtfA =  10 ,
--     mtfC =  500 ,
--     mtfB =  200 ,
--     type =  1 ,
--     mtfD =  5 
--  },
--   structure = {
--     zOffset =  0 ,
--     ncellz =  7 ,
--     rotateToZoneAxis =  false ,
--     isBoxed =  false ,
--     structure_filename =  ../../Examples/superstructures/goldball.gbm ,
--     crystalTiltY =  0.000000 ,
--     crystalTiltX =  0.000000 ,
--     yOffset =  0.000000 ,
--     zoneAxis =  0,1,0 ,
--     ncellx =  20 ,
--     boxX =  40 ,
--     temperatureK =  300.000000 ,
--     boxZ =  40 ,
--     boxY =  40 ,
--     crystalTiltZ =  0.000000 ,
--     xOffset =  0.000000 ,
--     ncelly =  20 
--  }
--}
