using IJulia
IJulia.installkernel("Julia 16 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "16",
))