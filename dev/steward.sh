#!/bin/sh
# The merge steward: one lane per day, main only moves on a passing gate.
#
#   dev/steward.sh sync    every 4 hours: merge origin/main into the work branches
#   dev/steward.sh gate N  after the day-N gate passes: write its marker
#   dev/steward.sh merge   at 00:05: merge the work branch into main when the
#                          previous day's marker is on it, then drop merged lanes
#
# The marker is `data/metrics/gate-green-YYYY-MM-DD` on the work branch.  A merge
# without yesterday's marker is refused: the gate is a report, the marker is the
# permission.
set -eu

cd "$(dirname "$0")/.."
work=program/3day-20261127
chore=chore/hygiene-20261127
mode=${1:-sync}

git fetch origin --prune

exists() {
  git rev-parse --verify --quiet "refs/heads/$1" >/dev/null
}

push_branch() {
  git push origin "HEAD:refs/heads/$1"
}

sync_one() {
  branch=$1
  exists "$branch" || return 0
  git checkout "$branch"
  if ! git merge --no-edit origin/main; then
    echo "steward: merge of origin/main into $branch conflicted; left in place" >&2
    exit 1
  fi
  push_branch "$branch"
}

case "$mode" in
sync)
  sync_one "$work"
  sync_one "$chore"
  ;;
gate)
  day=${2:?usage: dev/steward.sh gate N}
  Rscript dev/check.R --quick --no-render
  Rscript dev/check-user-api.R
  Rscript dev/check-mcp-parity.R
  git checkout "$work"
  mkdir -p data/metrics
  marker="data/metrics/gate-green-$(date +%F)"
  touch "$marker"
  git add "$marker"
  git commit -m "gate(day $day): marker for the merge steward"
  push_branch "$work"
  ;;
merge)
  day=${2:-$(date -d yesterday +%F)}
  marker="data/metrics/gate-green-$day"
  if ! git cat-file -e "origin/$work:$marker" 2>/dev/null; then
    echo "steward: no $marker on origin/$work; main stays where it is"
    exit 0
  fi
  git checkout main
  git pull --ff-only origin main
  git merge --no-ff -m "merge: gate $day passed" "origin/$work"
  if exists "$chore" && git rev-parse --verify --quiet "origin/$chore" >/dev/null; then
    git merge --no-ff -m "merge: hygiene lane for $day" "origin/$chore"
  fi
  push_branch main
  # Keep the two named lanes: they carry the next day.  Remove every other
  # program/chore branch that is now fully merged.
  for branch in $(git branch --merged main --format='%(refname:short)' |
    grep -E '^(program|chore)/' | grep -vx "$work" | grep -vx "$chore"); do
    git branch -d "$branch"
    git push origin --delete "$branch" || true
  done
  git checkout "$work"
  git merge --no-edit origin/main
  push_branch "$work"
  ;;
*)
  echo "usage: dev/steward.sh sync|gate N|merge" >&2
  exit 2
  ;;
esac
