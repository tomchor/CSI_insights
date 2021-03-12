using Parameters

@with_kw struct InteriorJetSimulations
    f0 = 1e-4
    Ny = 2^13
    Nz = 2^9
    Ly = 15_000 # m
    Lz = 500 # m

    JD15exp = (name = "JD15exp",
             f0 = f0,
             u₀ = 0.35e0,
             N2_inf = 4.9e-5,
             N2_pyc = 4.9e-5,
             Ny = Ny,
             Nz = Nz,
             Ly = 5000,
             Lz = Lz,
             σy = 1600,
             σz = 80,
             y₀ = Ly/2,
             z₀ = -Lz/2,
          )

    CIjet01 = (name = "CIintjet01",
              f0 = f0,
              u₀ = -0.4, # m/s
              N2_inf = 4e-5, # 1/s²
              N2_pyc = 4e-5, # 1/s²
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 1600, # m
              σz = 80, # m
              y₀ = Ly/3, # m
              z₀ = -Lz/2, # m
              )

end


@with_kw struct SurfaceJetSimulations
    f0 = 1e-4
    Ny = 2^14
    Nz = 2^8
    Ly = 15_000 # m
    Lz = 80 # m
    N2_pyc = 1e-6 # 1/s²

    CIjet1 = (name = "CIsurfjet1",
              f0 = f0,
              u₀ = -0.2,
              N2_inf = 1e-5,
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 800,
              σz = 80,
              y₀ = +Ly/2,
              z₀ = 0,
              N2_pyc = 1e-5,
             )

    SIjet1 = (name = "SIsurfjet1",
              f0 = f0,
              u₀ = -0.2, # m/s
              N2_inf = 1e-5, # 1/s²
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 1600, # m
              σz = 50, # m
              y₀ = +Ly/2, # m
              z₀ = 0, # m
              N2_pyc = 1e-5,
              )

    SIjet2 = (name = "SIsurfjet2",
              f0 = f0,
              u₀ = -0.2,
              N2_inf = 1e-6,
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 800,
              σz = 80,
              y₀ = +Ly/2,
              z₀ = 0,
              N2_pyc = 1e-6,
             )


    SIjet3 = (name = "SIsurfjet3",
              f0 = f0,
              u₀ = -0.2,
              N2_inf = 1e-6,
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 1400,
              σz = 80,
              y₀ = +Ly/2,
              z₀ = 0,
              N2_pyc = 1e-6,
             )



    SIjet4 = (name = "SIsurfjet4",
              f0 = f0,
              u₀ = -0.2,
              N2_inf = 1e-6,
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 1600,
              σz = 80,
              y₀ = +Ly/2,
              z₀ = 0,
              N2_pyc = 1e-6,
             )



    SIjet5 = (name = "SIsurfjet5",
              f0 = f0,
              u₀ = -0.1,
              N2_inf = 2.5e-7,
              Ny = Ny,
              Nz = Nz,
              Ly = Ly,
              Lz = Lz,
              σy = 800,
              σz = 80,
              y₀ = +Ly/2,
              z₀ = 0,
              N2_pyc = 2.5e-7,
             )



    Stabjet1 = (name = "stabsurfjet1",
               f0 = f0,
               u₀ = -0.08,
               N2_inf = 1e-5,
               Ny = Ny,
               Nz = Nz,
               Ly = Ly,
               Lz = Lz,
               σy = 800,
               σz = 80,
               y₀ = +Ly/2,
               z₀ = 0,
               N2_pyc = 1e-5,
             )


end

