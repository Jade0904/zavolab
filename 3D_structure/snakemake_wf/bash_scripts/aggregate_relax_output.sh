#!/bin/bash

out_paths=(${snakemake_input[out_path]})
sc_paths=(${snakemake_input[sc_path]})
aggregated_out_path=${snakemake_output[aggregated_out_path]}
aggregated_sc_path=${snakemake_output[aggregated_sc_path]}

echo "out_paths: ${out_paths[@]}"
echo "sc_paths: ${sc_paths[@]}"
echo "aggregated_out_path: $aggregated_out_path"
echo "aggregated_sc_path: $aggregated_sc_path"

out_header=""
out_content=""
sc_header=""
sc_content=""

for out_path in "${out_paths[@]}"; do
    current_out_header=$(awk '/^REMARK BINARY SILENTFILE/{print; exit} {print}' "$out_path")
    current_out_content=$(awk '/^REMARK BINARY SILENTFILE/{found=1; next} found' "$out_path")
    if [ -z "$out_header" ]; then
        out_header="$current_out_header"
    fi
    out_content+="${current_out_content}"$'\n'
done

printf "%s\n%s" "$out_header" "$out_content" > "$aggregated_out_path"

for sc_path in "${sc_paths[@]}"; do
    current_sc_header=$(awk '/^SCORE: total_score/{print; exit} {print}' "$sc_path")
    current_sc_content=$(awk '/^SCORE: total_score/{found=1; next} found' "$sc_path")
    if [ -z "$sc_header" ]; then
        sc_header="$current_sc_header"
    fi
    sc_content+="${current_sc_content}"$'\n'
    echo "sc_content: $sc_content"
done

printf "%s\n%s" "$sc_header" "$sc_content" > "$aggregated_sc_path"

