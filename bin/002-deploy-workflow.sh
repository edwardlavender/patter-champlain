# Set analysis = "sim" or analysis = "real" in scripts: 
# * 015-prepare-analysis.R
# * 001-run-algorithms.jl
# * 001-workflow.sh 

# Unlink old data/output/{analysis}/main/ (if needed):
# file.path("data", "output", "sim", "main", "blocks")
# file.path("data", "output", "real", "main", "blocks")

# Rebuild data/output/{analysis}/main/ directories:
# * Run 015-prepare-analysis.R

# Set NROW in 001-workflow.sh
# NROW=210
# NROW=2723

# Open tmux and run workflow:
# cd ~/Documents/work/projects/move-smc/patter/projects/patter-trout/patter-champlain
# cd ~/documents/projects/patter-champlain
# tmux ls
# tmux new -s ch
# chmod +x ./bin/001-workflow.sh
# ./bin/001-workflow.sh