using IJulia
IJulia.installkernel("Julia 8 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "8",
))