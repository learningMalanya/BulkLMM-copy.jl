using IJulia
IJulia.installkernel("Julia 4 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "4",
))