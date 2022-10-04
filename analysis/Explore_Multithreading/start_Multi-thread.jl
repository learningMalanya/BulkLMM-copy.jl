using IJulia
IJulia.installkernel("Julia 64 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "64",
))