#! /bin/bash

set -e

cd $(dirname $0)/..

if [ -f build_gz/.checksum.txt ]; then
  rm -rf build_gz_new 2> /dev/null || true
  mkdir build_gz_new
  pushd build > /dev/null
  md5sum -c ../build_gz/.checksum.txt 2>/dev/null | \grep ': OK$' | while read file; do
    file=${file%: OK}
    if [ -f ../build_gz/$file ]; then
      mkdir -p ../build_gz_new/$(dirname $file)
      mv ../build_gz/$file ../build_gz_new/$file.gz
    fi
  done
  popd > /dev/null
  rm -rf build_gz
  mv build_gz_new build_gz
elif [ -d build_gz ]; then
  # corrupted
  rm -rf build_gz
fi

pushd build > /dev/null
find . -type f | while read file; do
  if [ -f ../build_gz/$file.gz ]; then
    mv ../build_gz/$file.gz ../build_gz/$file
    continue
  fi
  mkdir -p ../build_gz/$(dirname $file)
  ln ../build/$file ../build_gz/$file
  pigz -9 -f ../build_gz/$file
  echo "compress: build_gz/$file"
  mv ../build_gz/$file.gz ../build_gz/$file
done
popd > /dev/null

pushd build > /dev/null
find . -type f -print0 | xargs -0 md5sum -b > ../build_gz/.checksum.txt
popd > /dev/null
