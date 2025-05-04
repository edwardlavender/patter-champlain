# Deploy analysis-patter.R on siam-linux20

# Instructions
# * Ensure test <- FALSE in constructor_ac_core() (!)
# * (optional) Clean up output folders, including logs/ and runs/
# * Customise arguments below as required
# * Run tmux code below to deploy script

# tmux cheatsheet:
# tmux ls
# tmux new -s $session_name
# tmux a                  (attach to single session)
# tmux a -t champlain-1   (attach to named session)
# tmux ctr-b d            (detach)

# Start tmux session
# cd ~/documents/projects/patter-champlain
# tmux new -s ch-1

# Define arguments
analysis="real"
mobility="216"
dev="FALSE"
DIRECTORY_LOG="data/output/$analysis/main/logs/R-CMD-BATCH"
mkdir -p "$DIRECTORY_LOG"

# Run R code
Rscript --verbose ./R/011-analysis-patter.R \
  "$analysis" "$mobility" "$dev" \
  > "$DIRECTORY_LOG/log-$mobility.Rout" 2>&1

