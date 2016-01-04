
local C = require("sliceth.sliceFFI")

local runner = {}

function runner.run(c)
  C.run_simulation(c)
end

return runner
