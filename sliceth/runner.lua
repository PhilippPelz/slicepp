
local C = require("sliceth.sliceFFI")

local runner = {}

function runner.run(c)
  print(c.ExperimentType)
  C.run_simulation(c)
end

return runner
