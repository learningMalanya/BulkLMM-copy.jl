using IJulia
IJulia.installkernel("Julia 24 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "24",
))