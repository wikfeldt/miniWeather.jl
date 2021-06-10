function run!(model, grid, output_freq, ncfilename)

    etime = 0.0
    sim_time = grid.nt * grid.dt

    if output_freq > 0.0
        nsnaps = floor(Int, sim_time/output_freq) + 1
        snapshots = Array{Float64}(undef, 4, grid.nx, grid.nz, nsnaps)
        etimes = Vector{Float64}(undef, nsnaps)    
        # output initial state
        snapshots[:,:,:,1] = createSnapshot(model, grid)
    end   

    direction_switch = true
    output_counter = 0.0
    counter = 1
    @showprogress for i in 1:grid.nt
        if output_freq > 0.0 && output_counter >= output_freq
            counter += 1
            output_counter = output_counter - output_freq
            snapshots[:,:,:,counter] = createSnapshot(model, grid)
            etimes[counter] = etime
        end

        perform_timestep!(model, grid, direction_switch)
        direction_switch = !direction_switch
        etime += grid.dt
        output_counter += grid.dt
        #output_gif!(model, etime, grid, anim);
    end


    if output_freq > 0
        output(snapshots, grid, etimes, ncfilename)
    end
    

end

#Performs a single dimensionally split time step using a simple low-storate three-stage Runge-Kutta time integrator
#The dimensional splitting is a second-order-accurate alternating Strang splitting in which the
#order of directions is alternated each time step.
#The Runge-Kutta method used here is defined as follows:
# q*     = q[n] + dt/3 * rhs(q[n])
# q^    = q[n] + dt/2 * rhs(q*  )
# q[n+1] = q[n] + dt/1 * rhs(q^ )
function perform_timestep!(model, grid, direction_switch)
    if direction_switch
        #x-direction first
        semi_discrete_step!(model, grid, DIR_X, 1)
        semi_discrete_step!(model, grid, DIR_X, 2)
        semi_discrete_step!(model, grid, DIR_X, 3)
        #z-direction second
        semi_discrete_step!(model, grid, DIR_Z, 1)
        semi_discrete_step!(model, grid, DIR_Z, 2)
        semi_discrete_step!(model, grid, DIR_Z, 3)
    else
        #z-direction second
        semi_discrete_step!(model, grid, DIR_Z, 1)
        semi_discrete_step!(model, grid, DIR_Z, 2)
        semi_discrete_step!(model, grid, DIR_Z, 3)
        #x-direction first
        semi_discrete_step!(model, grid, DIR_X, 1)
        semi_discrete_step!(model, grid, DIR_X, 2)
        semi_discrete_step!(model, grid, DIR_X, 3)
    end
end


#Perform a single semi-discretized step in time with the form:
#state_out = state_init + dt * rhs(state_forcing)
#Meaning the step starts from state_init, computes the rhs using state_forcing, and stores the result in state_out
function semi_discrete_step!(model, grid, dir, mode)
    # mode=1 sets model.state=state_init & state_forcing and model.state_tmp=state_out
    # mode=2 sets model.state=initial-state and model.state_tmp=state_forcing & state_out
    # mode=3 sets model.state=initial state & state_out and model.state_tmp=state_forcing      
    # time step depends on mode
    dt = grid.dt / (4.0 - mode)
    if dir == DIR_X
        #Set the halo values for this MPI task's fluid state in the x-direction
        if mode == 1
            set_halo_values_x!(
                model.state,
                model.hy_dens_cell,
                model.hy_dens_theta_cell,
                grid,
                model.weather_type,
            )
            compute_tendencies_x!(
                model.state,
                model.flux,
                model.tend,
                model.hy_dens_cell,
                model.hy_dens_theta_cell,
                grid,
            )
        elseif mode == 2 || mode == 3
            set_halo_values_x!(
                model.state_tmp,
                model.hy_dens_cell,
                model.hy_dens_theta_cell,
                grid,
                model.weather_type,
            )
            compute_tendencies_x!(
                model.state_tmp,
                model.flux,
                model.tend,
                model.hy_dens_cell,
                model.hy_dens_theta_cell,
                grid,
            )
        else
            throw(ArgumentError("mode must be either 1, 2 or 3"))
        end
        #Compute the time tendencies for the fluid state in the x-direction
    elseif dir == DIR_Z
        #Set the halo values for this MPI task's fluid state in the z-direction
        if mode == 1
            set_halo_values_z!(model.state, grid, model.weather_type)
            compute_tendencies_z!(
                model.state,
                model.flux,
                model.tend,
                model.hy_dens_int,
                model.hy_dens_theta_int,
                model.hy_pressure_int,
                grid,
            )
        elseif mode == 2 || mode == 3
            set_halo_values_z!(model.state_tmp, grid, model.weather_type)
            compute_tendencies_z!(
                model.state_tmp,
                model.flux,
                model.tend,
                model.hy_dens_int,
                model.hy_dens_theta_int,
                model.hy_pressure_int,
                grid,
            )
        else
            throw(ArgumentError("mode must be either 1, 2 or 3"))
        end
        #Compute the time tendencies for the fluid state in the z-direction
    else
        throw(ArgumentError("dir must be either $DIR_X or $DIR_Z, got $dir"))
    end

    #################################################
    ## TODO: THREAD ME
    #################################################
    #Apply the tendencies to the fluid state
    #state_out = zeros(size(model.state))
    state_out = similar(model.state)
    #    state_out = zeros(grid.nx + 2 * hs, grid.nz + 2 * hs, NUM_VARS)    
    for ll = 1:NUM_VARS
        for k = 1:grid.nz, i = 1:grid.nx         
            state_out[i+hs, k+hs, ll] =
                model.state[i+hs, k+hs, ll] + dt * model.tend[i, k, ll]            
        end
    end
    if mode == 1 || mode == 2
        model.state_tmp .= state_out
    elseif mode == 3
        model.state .= state_out
    else
        throw(ArgumentError("mode must be either 1, 2 or 3"))
    end
end


#Compute the time tendencies of the fluid state using forcing in the x-direction
#Since the halos are set in a separate routine, this will not require MPI
#First, compute the flux vector at each cell interface in the x-direction (including hyperviscosity)
#Then, compute the tendencies using those fluxes
function compute_tendencies_x!(state, flux, tend, hy_dens_cell, hy_dens_theta_cell, grid)
    #    flux = zeros(grid.nx + 1, grid.nz + 1, NUM_VARS)
    #    tend = zeros(grid.nx, grid.nz, NUM_VARS)
    d3_vals = zeros(NUM_VARS)
    vals = zeros(NUM_VARS)
    stencil = zeros(4)
    #Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * grid.dx / (16 * grid.dt)

    #################################################
    ## TODO: THREAD ME
    #################################################
    #Compute fluxes in the x-direction for each cell
    for k = 1:grid.nz, i = 1:grid.nx+1
        #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        for ll = 1:NUM_VARS
            for s = 1:sten_size
                stencil[s] = state[i+s-1, k+hs, ll]
            end
            #Fourth-order-accurate interpolation of the state
            vals[ll] =
                -stencil[1] / 12 + 7 * stencil[2] / 12 + 7 * stencil[3] / 12 -
                stencil[4] / 12
            #First-order-accurate interpolation of the third spatial derivative of the state (for artificial viscosity)
            d3_vals[ll] = -stencil[1] + 3 * stencil[2] - 3 * stencil[3] + stencil[4]
        end

        #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals[ID_DENS] + hy_dens_cell[k+hs]
        u = vals[ID_UMOM] / r
        w = vals[ID_WMOM] / r
        t = (vals[ID_RHOT] + hy_dens_theta_cell[k+hs]) / r
        p = C0 * (r * t)^gamma

        #Compute the flux vector
        flux[i, k, ID_DENS] = r * u - hv_coef * d3_vals[ID_DENS]
        flux[i, k, ID_UMOM] = r * u * u + p - hv_coef * d3_vals[ID_UMOM]
        flux[i, k, ID_WMOM] = r * u * w - hv_coef * d3_vals[ID_WMOM]
        flux[i, k, ID_RHOT] = r * u * t - hv_coef * d3_vals[ID_RHOT]
    end

    ######
    # TODO: THREAD ME
    #####
    #Use the fluxes to compute tendencies for each cell
    for ll = 1:NUM_VARS 
        for k = 1:grid.nz, i = 1:grid.nx
            tend[i, k, ll] = -(flux[i+1, k, ll] - flux[i, k, ll]) / grid.dx
        end
    end
end


#Compute the time tendencies of the fluid state using forcing in the z-direction
#Since the halos are set in a separate routine, this will not require MPI
#First, compute the flux vector at each cell interface in the z-direction (including hyperviscosity)
#Then, compute the tendencies using those fluxes
function compute_tendencies_z!(
    state,
    flux,
    tend,
    hy_dens_int,
    hy_dens_theta_int,
    hy_pressure_int,
    grid,
)
    #    flux = zeros(grid.nx + 1, grid.nz + 1, NUM_VARS)
    #    tend = zeros(grid.nx, grid.nz, NUM_VARS)
    d3_vals = zeros(NUM_VARS)
    vals = zeros(NUM_VARS)
    stencil = zeros(4)
    #Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * grid.dz / (16 * grid.dt)
    ####
    ## TODO: THREAD ME
    ####
    #Compute fluxes in the x-direction for each cell
    for k = 1:grid.nz+1, i = 1:grid.nx
        #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        for ll = 1:NUM_VARS
            for s = 1:sten_size
                stencil[s] = state[i+hs, k-1+s, ll]
            end
            #Fourth-order-accurate interpolation of the state
            vals[ll] =
                -stencil[1] / 12 + 7 * stencil[2] / 12 + 7 * stencil[3] / 12 -
                stencil[4] / 12
            #First-order-accurate interpolation of the third spatial derivative of the state
            d3_vals[ll] = -stencil[1] + 3 * stencil[2] - 3 * stencil[3] + stencil[4]
        end

        #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals[ID_DENS] + hy_dens_int[k]
        u = vals[ID_UMOM] / r
        w = vals[ID_WMOM] / r
        t = (vals[ID_RHOT] + hy_dens_theta_int[k]) / r
        p = C0 * (r * t)^gamma - hy_pressure_int[k]
        #Enforce vertical boundary condition and exact mass conservation
        if k == 1 || k == grid.nz + 1
            w = 0
            d3_vals[ID_DENS] = 0
        end

        #Compute the flux vector with hyperviscosity
        flux[i, k, ID_DENS] = r * w - hv_coef * d3_vals[ID_DENS]
        flux[i, k, ID_UMOM] = r * w * u - hv_coef * d3_vals[ID_UMOM]
        flux[i, k, ID_WMOM] = r * w * w + p - hv_coef * d3_vals[ID_WMOM]
        flux[i, k, ID_RHOT] = r * w * t - hv_coef * d3_vals[ID_RHOT]
    end

    ####
    ## TODO: THREAD ME
    ####
    #Use the fluxes to compute tendencies for each cell
    for ll = 1:NUM_VARS 
        for k = 1:grid.nz, i = 1:grid.nx
            tend[i, k, ll] = -(flux[i, k+1, ll] - flux[i, k, ll]) / grid.dz
            if ll == ID_WMOM
                tend[i, k, ll] -= state[i+hs, k+hs, ID_DENS] * grav
            end
        end
    end
end



#Set this MPI task's halo values in the x-direction. This routine will require MPI
function set_halo_values_x!(state, hy_dens_cell, hy_dens_theta_cell, grid, weather_type)
    #real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    ###
    # TODO: EXCHANGE HALO VALUES WITH NEIGHBORING MPI TASKS
    # (1) give    state(1:hs,1:nz,1:NUM_VARS)       to   my left  neighbor
    # (2) receive state(1-hs:0,1:nz,1:NUM_VARS)     from my left  neighbor
    # (3) give    state(nx-hs+1:nx,1:nz,1:NUM_VARS) to   my right neighbor
    # (4) receive state(nx+1:nx+hs,1:nz,1:NUM_VARS) from my right neighbor
    ###

    ###
    # DELETE THE SERIAL CODE BELOW AND REPLACE WITH MPI
    ###
    for ll = 1:NUM_VARS
        for k = 1:grid.nz
            state[1, hs+k, ll] = state[grid.nx+hs-1, hs+k, ll]
            state[2, hs+k, ll] = state[grid.nx+hs, hs+k, ll]
            state[grid.nx+hs+1, hs+k, ll] = state[hs+1, hs+k, ll]
            state[grid.nx+hs+2, hs+k, ll] = state[hs+2, hs+k, ll]
        end
    end
    ###

    if (weather_type == "injection")
        #if (myrank == 0)
        for k = 1:grid.nz
            z = (k_beg - 1 + k - 0.5) * grid.dz
            if (abs(z - 3 * zlen / 4) <= zlen / 16)
                state[1:2, hs+k, ID_UMOM] =
                    (state[1:2, hs+k, ID_DENS] + hy_dens_cell[k+hs]) * 50.0
                state[1:2, hs+k, ID_RHOT] =
                    (state[1:2, hs+k, ID_DENS] + hy_dens_cell[k+hs]) * 298.0 -
                    hy_dens_theta_cell[k+hs]
            end
        end
        #end
    end
end


#Set this MPI task's halo values in the z-direction. This does not require MPI because there is no MPI
#decomposition in the vertical direction
function set_halo_values_z!(state, grid, weather_type)
    #real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    mnt_width = xlen / 8
    ###
    # TODO: THREAD ME
    ###
    for ll = 1:NUM_VARS
        for i = 1:grid.nx+2*hs
            if (ll == ID_WMOM)
                state[i, 1, ll] = 0
                state[i, 2, ll] = 0
                state[i, grid.nz+hs+1, ll] = 0
                state[i, grid.nz+hs+2, ll] = 0
                #Impose the vertical momentum effects of an artificial cos^2 mountain at the lower boundary
                if (weather_type == "mountain_waves")
                    x = (i_beg - 1 + i - 0.5) * grid.dx
                    if (abs(x - xlen / 4) < mnt_width)
                        xloc = (x - (xlen / 4)) / mnt_width
                        #Compute the derivative of the fake mountain
                        mnt_deriv =
                            -pi * cos(pi * xloc / 2) * sin(pi * xloc / 2) * 10 / grid.dx
                        #w = (dz/dx)*u
                        state[i, 1, ID_WMOM] = mnt_deriv * state[i, 3, ID_UMOM]
                        state[i, 2, ID_WMOM] = mnt_deriv * state[i, 3, ID_UMOM]
                    end
                end
            else
                state[i, 1, ll] = state[i, 3, ll]
                state[i, 2, ll] = state[i, 3, ll]
                state[i, grid.nz+3, ll] = state[i, grid.nz+2, ll]
                state[i, grid.nz+4, ll] = state[i, grid.nz+2, ll]
            end
        end
    end
end

