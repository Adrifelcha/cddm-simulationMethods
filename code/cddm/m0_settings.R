settings <- list(
                iterations = 200,
                sampleSize.list  =  c(500),
                driftAngle.list  =  c(0.0,  2.0, 4.0),
                driftLength.list =  c(0.01, 1.0, 2.0),
                bound.list       =  c(1.5, 2.0, 2.5),
                nondecision.list =  c(0.1, 0.2, 0.3))
settings <- c(settings,
              s.topIdx  = length(settings$sampleSize.list),
              a.topIdx  = length(settings$driftAngle.list),
              m.topIdx  = length(settings$driftLength.list),   # 'm' for magnitude
              b.topIdx  = length(settings$bound.list),
              n.topIdx  = length(settings$nondecision.list))
settings <- c(settings,
              possible.combinations = settings$s.topIdx * settings$a.topIdx * settings$m.topIdx * settings$b.topIdx * settings$n.topIdx,
              n.iter    = 2500,
              n.burnin  = 500,
              n.chains  = 4,
              n.thin    = 1,
              output.folder = "./output_polar_n500/",
              simstudy.Name = "sim_")