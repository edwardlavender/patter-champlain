# Batch particle algorithms to handle memory requirements

# Memory available: see top
# Memory required:
# * Filter starts with a vector of .xinit state instances.
# * At each time step, we hold past and current state instances (vector).
# * The filter accumulates a .nr * .nt matrix of state instances in memory.
# * We need two filter runs and one smoother run.
# * I.e., we need to hold three .nr * .nt matrix in memory.

# Compute approximate memory used to store outputs of filtering & smoothing
p_mem <- function(.nr, .nt, .ns, .nc) {
  3 * (.nr * .nt * .ns * 8) * .nc / 1e6
}

# Compute the number of batches required
p_batch <- function(.mem_req, .mem_avail) {
  max(c(1L, ceiling(.mem_req / .mem_avail)))
}

# Compute approximate memory used to store outputs of filtering & smoothing, if:
# * 1000,1500, 2000 particles recorded
# * 22320 time steps
# * 4 state dimensions (map_value, x, y, heading)
# * 100 cores
p_mem(1000, 22320, 4, 100) # 214272 MB = 214 GB
p_mem(1500, 22320, 4, 100) # 321408 MB = 321 GB 
p_mem(2000, 22320, 4, 100) # 428544 MB = 428 GB 

# Compute required number of batches on siam-linux20: 
# * Double required memory (safety buffer)
# * Double required memory for iterative applications 
#   (we may start on while nearly finishing another)
# * Assume 75 % of ~400 000 MB server memory available
#   (~400 000 MB is max available assuming server not in use)
# > This suggests 3 batches is sufficient. 
p_batch(214272 * 4, 4e5 * 0.75) # 3 batches if 1000 particles recorded
p_batch(321408 * 4, 4e5 * 0.75) # 5 batches if 1500 particles recorded
p_batch(428544 * 4, 4e5 * 0.75) # 6 batches if 2000 particles recorded

# Compute number of batches if N GB memory available:
# > Using 1 core
p_batch(p_mem(1000, 22320, 4, 1), 4e3)
# > Using 2 cores
p_batch(p_mem(1000, 22320, 4, 2), 4e3)
# > Using more cores
p_batch(p_mem(1000, 22320, 4, 3), 4e3)
p_batch(p_mem(1000, 22320, 4, 4), 4e3)
p_batch(p_mem(1000, 22320, 4, 5), 4e3)
