#!/usr/bin/env bash
set -euo pipefail

wiki_dir=${1:-wiki}
repository=${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is required}
server_url=${GITHUB_SERVER_URL:-https://github.com}
source_sha=${GITHUB_SHA:?GITHUB_SHA is required}
wiki_token=${WIKI_TOKEN:?WIKI_TOKEN is required}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
"$script_dir/wiki_check.sh" "$wiki_dir"

temp_root=$(mktemp -d "${RUNNER_TEMP:-/tmp}/z-fasta-wiki.XXXXXX")
trap 'rm -rf -- "$temp_root"' EXIT

wiki_repo="$temp_root/wiki"
wiki_url="${server_url%/}/${repository}.wiki.git"
auth_header=$(printf 'x-access-token:%s' "$wiki_token" | base64 | tr -d '\n')

git -c "http.extraheader=AUTHORIZATION: basic ${auth_header}" \
  clone --quiet "$wiki_url" "$wiki_repo"

wiki_branch=$(git -C "$wiki_repo" symbolic-ref --quiet --short HEAD) || {
  echo "wiki repository has no checked-out default branch; create its first page in GitHub" >&2
  exit 1
}

rsync --archive --delete --exclude='.git/' "$wiki_dir/" "$wiki_repo/"

git -C "$wiki_repo" add --all
if git -C "$wiki_repo" diff --cached --quiet --exit-code; then
  echo "wiki is already synchronized"
  exit 0
fi

echo "wiki changes to publish:"
git -C "$wiki_repo" diff --cached --name-status

git -C "$wiki_repo" config user.name 'github-actions[bot]'
git -C "$wiki_repo" config user.email '41898282+github-actions[bot]@users.noreply.github.com'
git -C "$wiki_repo" commit --quiet -m "Sync wiki from ${source_sha}"

git -c "http.extraheader=AUTHORIZATION: basic ${auth_header}" \
  -C "$wiki_repo" push --quiet origin "HEAD:${wiki_branch}"

local_sha=$(git -C "$wiki_repo" rev-parse HEAD)
remote_sha=$(git -c "http.extraheader=AUTHORIZATION: basic ${auth_header}" \
  -C "$wiki_repo" ls-remote origin "refs/heads/${wiki_branch}" | awk 'NR == 1 { print $1 }')

if [[ -z "$remote_sha" || "$remote_sha" != "$local_sha" ]]; then
  echo "wiki push verification failed: local=$local_sha remote=${remote_sha:-missing}" >&2
  exit 1
fi

echo "published wiki commit $local_sha to $repository ($wiki_branch)"
