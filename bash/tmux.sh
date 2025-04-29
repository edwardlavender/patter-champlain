# Deploy analysis-patter.R on siam-linux20

# tmux cheatsheet:
# tmux ls
# tmux new -s $session_name
# tmux a                  (attach to single session)
# tmux a -t champlain-1   (attach to named session)
# tmux ctr-b d            (detach)

cd ~/documents/projects/patter-champlain
tmux new -s champlain-1

analysis="sim"
mobility="116"
dev="FALSE"
DIRECTORY_LOG="data/output/$analysis/main/logs/R-CMD-BATCH"

mkdir -p "$DIRECTORY_LOG"

R CMD BATCH \
  --no-save \
  --no-restore \
  ./R/012-analysis-patter.R \
  "$DIRECTORY_LOG/log-$mobility.Rout" \
  --args "$analysis" "$mobility" "$dev"

