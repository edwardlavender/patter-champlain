# Deploy Julia workflows on siam-linux20

# Instructions
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
# tmux ls
# tmux new -s ch-1

# TO DO Run Julia scripts

# for i in {1..50}; do
#   Rscript run_model.R "$i" &
#   (( $(jobs -r | wc -l) >= 10 )) && wait -n
# done
# wait

# parallel -j 10 Rscript run.R ::: args

