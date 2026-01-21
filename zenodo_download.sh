#!/bin/bash

download_file() {
  local url=$1
  local output_file=$2
  wget --progress=dot "$url" -O "$output_file" 2>&1 | \
  grep --line-buffered "%" | \
  sed -u -e "s,\.,,g" | awk '{printf("\rDownloading %s: %s", "'$output_file'", $2)}'
  echo " Done."
}

extract_tar_gz() {
  local tar_file=$1
  local output_dir=$2
  echo -n "Extracting $tar_file..."
  tar -xvf "$tar_file" -C "$output_dir" > /dev/null
  echo "Done."
}

remove_file() {
  local file=$1
  echo -n "Removing $file..."
  rm -f "$file"
  echo "Done."
}

remove_directory() {
  local dir=$1
  if [ -d "$dir" ]; then
    echo -n "Removing existing directory $dir..."
    rm -rf "$dir"
    echo "Done."
  fi
}

ZENODO_BASE_URL="https://zenodo.org/record/18312238/files"
FILES=("data.tar.gz" "models.tar.gz")
DIRS=("data" "models")
REPO_DIR=$(pwd)

for i in "${!FILES[@]}"; do
  file="${FILES[$i]}"
  dir="${DIRS[$i]}"
  
  download_file "$ZENODO_BASE_URL/$file" "$REPO_DIR/$file"
  remove_directory "$REPO_DIR/$dir"
  extract_tar_gz "$REPO_DIR/$file" "$REPO_DIR"
  remove_file "$REPO_DIR/$file"
done

echo "Data and models were downloaded from Zenodo successfully."

