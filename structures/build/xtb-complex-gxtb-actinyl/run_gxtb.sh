#!/bin/bash
# Step 3: g-xTB gas-phase optimisation of every ETKDG seed in tasks.txt -> runs/<idx>_<El>_s<seed>/.
# 9 workers, 1 thread each (2 cores free). Resumable: a run whose log says the geometry optimisation
# converged is skipped. usage: conda activate env-gxtb && bash xtb-complex-gxtb-ox4/run_gxtb.sh [nproc]
set -u
cd "$(dirname "$0")"
case "$(command -v xtb)" in
  */env-gxtb/opt/gxtb/bin/xtb) ;;
  *) echo "xtb is not the g-xTB binary ($(command -v xtb)); activate env-gxtb" >&2; exit 1 ;;
esac
NPROC=${1:-9}
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_STACKSIZE=4G

one() {
  name=$1; uhf=$2; chrg=$3; d=runs/$name
  if [ -f "$d/xtbopt.xyz" ] && grep -q 'GEOMETRY OPTIMIZATION CONVERGED' "$d/xtb.log" 2>/dev/null; then
    echo "$name uhf=$uhf chrg=$chrg rc=0 wall_s=0 skipped=1" >> runs/status.txt; return
  fi
  mkdir -p "$d"; cp "starts/$name.xyz" "$d/in.xyz"
  t0=$(date +%s)
  (cd "$d" && xtb in.xyz --gxtb --opt --chrg "$chrg" --uhf "$uhf" > xtb.log 2>&1)
  echo "$name uhf=$uhf chrg=$chrg rc=$? wall_s=$(( $(date +%s) - t0 )) skipped=0" >> runs/status.txt
}
export -f one

mkdir -p runs
touch runs/status.txt
# tasks.txt lines: name uhf chrg (written by build_list.py)
xargs -P "$NPROC" -L 1 bash -c 'one "$0" "$1" "$2"' < tasks.txt
echo done
