using Parameters

@with_kw struct SurfaceJetSimulations
    f_0 = 1e-4
    Ny = 10*2^9
    Nz = 2^9
    Ny_2d = 100*2^9
    Nz_2d = 2^9
    Ly = 8_000 # m
    Lz = 80 # m
    νz = 5e-4
    sponge_frac = 1/16
    ThreeD = false

    CItest1 = (name = "CItest1",
               f_0 = f_0,
               u_0 = -0.2,
               N2_inf = 1e-5,
               Ny = ThreeD ? Ny : Ny_2d,
               Nz = ThreeD ? Nz : Nz_2d,
               Ly = Ly,
               Lz = Lz,
               σ_y = 800,
               σ_z = 80,
               y_0 = +Ly/2,
               z_0 = 0,
               νz = 5e-4,
               sponge_frac = sponge_frac,
               )

    CIfront1 = (name = "CIfront1",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 1e-5,
                  Ny = ThreeD ? Ny : Ny_2d,
                  Nz = ThreeD ? Nz : Nz_2d,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 800,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = 5e-4,
                  sponge_frac = sponge_frac,
                 )

    CIfront2 = (name = "CIfront2",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 5e-5,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 800,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = νz,
                  sponge_frac = sponge_frac,
                 )

    CIfront3 = (name = "CIfront3",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 5e-6,
                  Ny = ThreeD ? Ny : Ny_2d,
                  Nz = ThreeD ? Nz : Nz_2d,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 800,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = νz,
                  sponge_frac = sponge_frac,
                 )

    CIfront4 = (name = "CIfront4",
                  f_0 = 5e-5,
                  u_0 = -0.2,
                  N2_inf = 5e-6,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 800,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = νz,
                  sponge_frac = sponge_frac,
                 )


    CIfront5 = (name = "CIfront5",
                  f_0 = 7e-5,
                  u_0 = -0.2,
                  N2_inf = 1.4e-6,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 600,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = νz,
                  sponge_frac = sponge_frac,
                 )


    SIfront1 = (name = "SIfront1",
                  f_0 = f_0,
                  u_0 = -0.23, # m/s
                  N2_inf = 5e-6, # 1/s²
                  Ny = ThreeD ? Ny : 100*2^9,
                  Nz = ThreeD ? Nz : 2^9,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 1600, # m
                  σ_z = 80, # m
                  y_0 = +Ly/2, # m
                  z_0 = 0, # m
                  νz = νz,
                  sponge_frac = sponge_frac,
                  )
    
    SIfront2 = (name = "SIfront2",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 1e-6,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 800,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = νz,
                  sponge_frac = sponge_frac,
                 )
    
    
    SIfront3 = (name = "SIfront3",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 1.4e-6,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 1400,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = 8e-4,
                  sponge_frac = sponge_frac,
                 )
    
    
    
    SIfront4 = (name = "SIfront4",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 1e-6,
                  Ny = ThreeD ? Ny : Ny_2d,
                  Nz = ThreeD ? Nz : Nz_2d,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 1600,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = 1e-3,
                  sponge_frac = sponge_frac,
                 )
    
    
    
    SIfront5 = (name = "SIfront5",
                  f_0 = f_0,
                  u_0 = -0.1,
                  N2_inf = 2.5e-7,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 800,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = 1e-3,
                  sponge_frac = sponge_frac,
                  )
    
    
    
    SIfront6 = (name = "SIfront6",
                  f_0 = f_0,
                  u_0 = -0.2,
                  N2_inf = 2.5e-6,
                  Ny = Ny,
                  Nz = Nz,
                  Ly = Ly,
                  Lz = Lz,
                  σ_y = 1200,
                  σ_z = 80,
                  y_0 = +Ly/2,
                  z_0 = 0,
                  νz = νz,
                  sponge_frac = sponge_frac,
                  )
    
end

