using IJulia
IJulia.installkernel("Julia 12 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "12",
))