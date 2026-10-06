#!/bin/bash
# Compile-only artifact extraction; does not launch a GPU application.
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof
llvm_root=/opt/rocm-6.4.2/llvm/bin
for variant in old new; do
  "$llvm_root/llvm-objcopy" --dump-section=".hip_fatbin=hip/$variant.fatbin" "hip/$variant.o"
  "$llvm_root/clang-offload-bundler" --unbundle --type=bc \
    --input="hip/$variant.fatbin" --targets=hipv4-amdgcn-amd-amdhsa--gfx90a \
    --output="hip/$variant.hsaco"
  for section in text rodata; do
    "$llvm_root/llvm-objcopy" --dump-section=".$section=hip/$variant-device.$section" "hip/$variant.hsaco"
  done
  "$llvm_root/llvm-objdump" --triple=amdgcn-amd-amdhsa --mcpu=gfx90a -d \
    "hip/$variant.hsaco" > "hip/$variant-device.disassembly"
done
python3 audit_proof.py
