# Deploy analysis-patter.R on siam-linux20

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
  --args "$analysis" "$mobility" "$dev" \
  ./R/012-analyse-patter.R \
  "$DIRECTORY_LOG/log-$mobility.Rout"
