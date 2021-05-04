
#Compute reduced quantities for error checking without resorting to the "ncdiff" tool
function reductions(model, grid)
    mass = 0.0
    te = 0.0
    nx, nz = grid.nx, grid.nz
    dx, dz = grid.dx, grid.dz

    for k = 1:nz, i = 1:nx
        r = model.state[hs+i, hs+k, ID_DENS] + model.hy_dens_cell[hs+k]       # Density
        u = model.state[hs+i, hs+k, ID_UMOM] / r                           # U-wind
        w = model.state[hs+i, hs+k, ID_WMOM] / r                           # W-wind
        th = (model.state[hs+i, hs+k, ID_RHOT] + model.hy_dens_theta_cell[hs+k]) / r # Potential Temperature (theta)
        p = C0 * (r * th)^gamma      # Pressure
        t = th / (p0 / p)^(rd / cp)  # Temperature
        ke = r * (u * u + w * w)           # Kinetic Energy
        ie = r * cv * t                # Internal Energy
        mass = mass + r * dx * dz # Accumulate domain mass
        te = te + (ke + ie) * dx * dz # Accumulate domain total energy
        #println(r, u, w, th)
      end

    #    double glob[2], loc[2];
    #    loc[0] = mass;
    #    loc[1] = te;
    #    int ierr = MPI_Allreduce(loc,glob,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    #    mass = glob[0];
    #    te   = glob[1];
    return mass, te
end

function createSnapshot(model, grid)

    snapshot = Array{Float64}(undef, 4, grid.nx, grid.nz)
    #Store perturbed values in the temp arrays for output
    for k = 1:grid.nz
        for i = 1:grid.nx
            snapshot[1, i, k] =  model.state[i, hs+k, ID_DENS]
            snapshot[2, i, k] =
                model.state[i, hs+k, ID_UMOM] /
                (model.hy_dens_cell[hs+k] + model.state[i, hs+k, ID_DENS])
            snapshot[3, i, k] =
                model.state[i, hs+k, ID_WMOM] /
                (model.hy_dens_cell[hs+k] + model.state[i, hs+k, ID_DENS])
            snapshot[4, i, k] =
                (model.state[i, hs+k, ID_RHOT] + model.hy_dens_theta_cell[hs+k]) /
                (model.hy_dens_cell[hs+k] + model.state[i, hs+k, ID_DENS]) -
                model.hy_dens_theta_cell[hs+k] / model.hy_dens_cell[hs+k]
        end
    end

    return snapshot
end


#
#Output the fluid state (state) to a NetCDF file at a given elapsed model time (etime)
#The file I/O uses parallel-netcdf, the only external library required for this mini-app.
#If it's too cumbersome, you can comment the I/O out, but you'll miss out on some potentially cool graphics
function output(snapshots, grid, etimes, ncfile="output.nc")
  #Create new file
  isfile(ncfile) && rm(ncfile)
  _, nx, nz, _ = size(snapshots)
  nt = length(etimes)
  nccreate(ncfile, "dens", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))
  nccreate(ncfile, "uwnd", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))
  nccreate(ncfile, "wwnd", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))
  nccreate(ncfile, "theta", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))

                           #nccreate(ncfile, "uwnd", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
  #nccreate(ncfile, "wwnd", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
  #nccreate(ncfile, "theta", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))    
  #exit()
  #Write the grid data to file with all the processes writing collectively
#  ncwrite(etimes, ncfile, "time")
  ncwrite(snapshots[1, :, :, 1:nt], ncfile, "dens")
  ncwrite(snapshots[2, :, :, 1:nt], ncfile, "uwnd")
  ncwrite(snapshots[3, :, :, 1:nt], ncfile, "wwnd")
  ncwrite(snapshots[4, :, :, 1:nt], ncfile, "theta")

end


