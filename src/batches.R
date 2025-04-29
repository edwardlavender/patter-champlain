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
# * 1000 particles
# * 20,000 time steps
# * 4 state dimensions (map_value, x, y, heading)
# * 100 cores
p_mem(1000, 20000, 4, 100) # 192000 MB = 192 GB 

# Compute required number of batches on siam-linux20: 
# * Double required memory (safety buffer)
# * Double required memory for iterative applications (we may start on while nearly finishing another)
# * Assume 500,000 MB available memory (if server not in use)
# > This suggests 3 batches is sufficient. 
p_batch(192000 * 4, 5e5)
