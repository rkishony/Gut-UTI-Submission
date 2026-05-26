#!/usr/bin/env bash
set -euo pipefail

SESSION_NAME="matlab"
MATLAB_BIN="/usr/local/MATLAB/R2021b/bin/matlab"
MATLAB_CMD="${MATLAB_BIN} -nodisplay -nosplash -nodesktop"

if tmux has-session -t "${SESSION_NAME}" 2>/dev/null; then
  # Session exists. Check if MATLAB is still running inside it.
  matlab_pid=$(tmux list-panes -t "${SESSION_NAME}" -F '#{pane_pid}' 2>/dev/null | head -1)
  matlab_running=false
  if [ -n "${matlab_pid}" ]; then
    # Check if any child of the tmux pane is MATLAB
    if pgrep -P "${matlab_pid}" -f MATLAB >/dev/null 2>&1; then
      matlab_running=true
    fi
  fi

  if ${matlab_running}; then
    echo "MATLAB is running in tmux session '${SESSION_NAME}'."
    echo "If fonts are wrong (stale DISPLAY), type 'exit' in MATLAB, then re-run this script."
    exec tmux attach -t "${SESSION_NAME}"
  else
    # Session exists but MATLAB exited. Start fresh MATLAB with current DISPLAY.
    tmux send-keys -t "${SESSION_NAME}" "${MATLAB_CMD}" Enter
    exec tmux attach -t "${SESSION_NAME}"
  fi
else
  # No session. Create one and start MATLAB with current DISPLAY.
  exec tmux new -s "${SESSION_NAME}" "${MATLAB_CMD}"
fi
