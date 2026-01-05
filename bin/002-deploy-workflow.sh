# Set analysis = "sim" or analysis = "real" in scripts: 
# * 015-prepare-analysis.R
# * 001-run-algorithms.jl
# * 001-workflow.sh 

# Unlink old data/output/{analysis}/main/ (if needed):
# file.path("data", "output", "sim", "main", "runs")
# file.path("data", "output", "real", "main", "runs")

# Rebuild data/output/{analysis}/main/ directories:
# * Run 015-prepare-analysis.R

# Set NROW in 001-workflow.sh
# NROW=210
# NROW=2723

# Open tmux and run workflow:
# cd ~/documents/projects/patter-champlain
# tmux ls
# tmux new -s ch
# chmod +x ./bin/001-workflow.sh
# ./bin/001-workflow.sh