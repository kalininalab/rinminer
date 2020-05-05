#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset
export LANG=C
export LC_ALL=C

# number of entries in a family required to be considered
min_entries='3'

# increment used for SMLF support percentage labels
increment='5'

# maximum number of selected results for each run
max_subgraph_results='15'

# parameteres to limit resource usage
# 80% of system memory
max_memory_usage="$(echo $(awk -F" " '$1~/MemTotal:/{print $2}' /proc/meminfo) \* 9 \/ 10 | bc)"
# all available cores
threads="$(nproc)"
min_output_size='3'
max_output_file_size='300M'
max_run_time='1h'

default_support_percentage_start='50'
default_support_percentage_end='100'

# rinminer parameters
default_edge_inclusion='200'
default_min_score_difference='900'

# significance trheeshold before multiple testing correction
default_significance_base_threshold='0.0001'

run_script_dir="$(dirname "$(realpath -s "$0")")"
rinminer_path="$(realpath -m -s "${run_script_dir}/../rinminer/rinminer")"
scripts_path="$(realpath -s "${run_script_dir}/..")"

rerun=0
rerun_crashed=0
rerun_select=0
rerun_eval=0
use_patterns=0
use_significance=0
use_smlf=0
use_sm_classes=0

edge_inclusion=${default_edge_inclusion}
min_score_difference=${default_min_score_difference}
significance_base_threshold=${default_significance_base_threshold}
support_percentage_start=${default_support_percentage_start}
support_percentage_end=${default_support_percentage_end}


positional=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --rerun)
    rerun=1
    shift
    ;;
    --rerun-crashed)
    rerun_crashed=1
    shift
    ;;
    --rerun-select)
    rerun_select=1
    shift
    ;;
    --rerun-eval)
    rerun_eval=1
    shift
    ;;
    --use-patterns)
    use_patterns=1
    shift
    ;;
    --smlf)
    use_smlf=1
    shift
    ;;
    --structman)
    use_sm_classes=1
    shift
    sm_classes_file="$1"
    shift
    ;;
    --memory)
    shift
    max_memory_usage="$1"
    shift
    ;;
    --threads)
    shift
    threads="$1"
    shift
    ;;
    --min-size)
    shift
    min_output_size="$1"
    shift
    ;;
    --file-size)
    shift
    max_output_file_size="$1"
    shift
    ;;
    --run-time)
    shift
    max_run_time="$1"
    shift
    ;;
    --edge-match-dist)
    shift
    edge_inclusion="$1"
    shift
    ;;
    --min-score-diff)
    shift
    min_score_difference="$1"
    shift
    ;;
    --sig-threshold)
    shift
    significance_base_threshold="$1"
    shift
    ;;
    --sup-start)
    shift
    support_percentage_start="$1"
    shift
    ;;
    --sup-end)
    shift
    support_percentage_end="$1"
    shift
    ;;
    *)
    positional+=("$1")
    shift
    ;;
esac
done
set -- "${positional[@]}"

if [ $# -lt 3 -o $# -gt 4 ]
then
    echo "$0 [OPTIONS] /path/to/dataset /path/to/results /path/to/report [/path/to/classification-database]"
    echo "Running on a data set with existing results will not regenerate them unless specified by an option."
    echo "Options:"
    echo -e "--rerun\t\t\t\trerun for previously successfully run databases"
    echo -e "--rerun-crashed\t\t\trerun for previously failed databases"
    echo -e "--rerun-select\t\t\trerun subgraph selection"
    echo -e "--rerun-eval\t\t\trerun subgraph evaluation and analysis"
    echo -e "--use-patterns\t\t\tcompare found subgraphs to PROSITE style patterns"
    echo -e "--smlf\t\t\t\tgenerate SMLF files for analysis with StructMAn"
    echo -e "--structman simple_classes.tsv\tuse the provided simple class distribution generate by structman to annotate results"
    echo -e "--memory x[k,M,G]\t\tmaximum memory usage of RINminer"
    echo -e "--threads n\t\t\tnumber of threads used by RINminer"
    echo -e "--min-size n\t\t\tminimum number of edges for a subgraph to be included in the output"
    echo -e "--file-size n\t\t\tmaximum file size for the output of a single RINminer run"
    echo -e "--run-time x[m,h]\t\tmaximum run time for a single RINminer run"
    echo -e "--edge-match-dist n\t\tmaximum edge label distance difference considered for matches. As per mille of average distance. Default: ${default_edge_inclusion}"
    echo -e "--min-score-diff n\t\tminimum score increase compared to parent graphs. Default ${default_min_score_difference}"
    echo -e "--sig-threshold n\t\threshold for classification significance (before multiple testing correction). Requires classification database. Default: ${default_significance_base_threshold}"
    echo -e "--sup-start n\t\support range start. As integer of percentages. Default: ${default_support_percentage_start}"
    echo -e "--sup-end n\t\support range end. As integer of percentages. Default: ${default_support_percentage_end}"
    exit 1
fi

data_path="$1"
shift
out_path="$1"
shift
report_path="$1"
shift

if [ $# -eq 1 ]
then
    significance_path="$1"
    if [[ -f "${significance_path}" ]]
    then
        use_significance=1
    fi
else
    significance_path=""
fi

tmp_path="$(mktemp -d -t rinminer.XXXXX)"
function cleanup {
    rm -rf "${tmp_path}"
}
trap cleanup EXIT

max_memory_usage_kbytes=$(echo "${max_memory_usage}" | sed -e 's/[iIbB]*//g' -e 's/[kK]//g' -e 's/[mM]/\*1024/g' -e 's/[gG]/\*1024\*1024/g' -e 's/[tT]/\*1024\*1024\*1024/g' | bc)
ulimit -Sv "${max_memory_usage_kbytes}"

# disable core dumps
ulimit -Sc 0

echo -e "Memory limit:\t\t$(echo "${max_memory_usage_kbytes} / 1024" | bc)MiB"
echo -e "Threads:\t\t${threads}"
echo -e "Max file size:\t\t${max_output_file_size}"
echo -e "Max run time:\t\t${max_run_time}"
echo -e "Min graph size:\t\t${min_output_size}"
echo -e "Num selected graphs:\t${max_subgraph_results}"
echo -e "Max edge distance diff:\t${edge_inclusion}"
echo -e "Min score increase:\t${min_score_difference}"
if [[ "${use_significance}" -eq "1" ]]
then
    echo -e "Significance threshold:\t${significance_base_threshold} (before multiple testing correction)"
fi
echo -e "Support range:\t${support_percentage_start}% - ${support_percentage_end}%"

if [[ ! -f "${rinminer_path}" ]]
then
    echo "Error: no rinminer installation found under ${rinminer_path}"
    exit 1
fi

pretty_names=()
results_files=()
tmp_results_files=()
output_files=()
tmp_output_files=()
error_files=()
select_files=()
index_files=()
graph_dirs=()
commands=()
tags=()

family_names=()
family_ids=()
family_patterns=()
range_start=()
range_end=()

### report output stuff ########################################################

function join_by {
    local IFS="$1"
    shift
    echo "$*"
}

function report_header {
    local report_file="$1"
    local significance_parameter=""

    if [[ "${use_significance}" -eq "1" ]]
    then
        significance_parameter="<b>Significance threshold:</b> ${significance_threshold_string} = ${significance_threshold}<br>"
    fi

    cat > "${report_file}" <<EOF
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"> 
<title>Results Overview</title>
<style>
    body { font-family: sans-serif; color: black; background-color: white; }
    table { border-collapse: collapse; table-layout: fixed; }
    table, th, td { border: 1px solid #bbb; }
    th, td { min-width: 6.5em; padding: 4px; white-space: nowrap; }
    th { background-color: #f3f3f3; width: auto; }
    td { text-align: center; }
    td.interesting { background-color: rgb(136, 235, 140); }
    td.interesting.bad-score { background: repeating-linear-gradient(135deg,
                                                                     rgb(136, 235, 140), rgb(136, 235, 140) 14.2857%,
                                                                     transparent 14.2857%, transparent 28.5714%); }
    td.significant { box-shadow: 0px 0px 0 3px #fe3 inset; }
    td.failed { background-color: rgb(235, 136, 140); }
    td.empty { background-color: rgb(220, 220, 220); }
    td a { color: #449; }
    td a:visited { color: #b2a; }
    td.no-graphs { background-color: #f8f8f8; color: #888; }
    td.no-graphs a { color: #ccc; }
    td.no-graphs a:visited { color: #eae; }
</style>
</head>

<body>
<h1>Found Subgraphs</h1>
<span style="font-size: larger;"><b>Parameters:</b></span><br>
<b>Edge inclusion threshold (approximate isomorphism):</b> $(printf "%.1lf" "$(echo "${edge_inclusion}/10" | bc -q -l)")% of average distance<br>
<b>Minimum score difference:</b> ${min_score_difference}<br>
${significance_parameter}
Minimum database size: ${min_entries}<br>
Minimum subgraph size in output: ${min_output_size}<br>
Maximum output file size: ${max_output_file_size}<br>
Maximum run time: ${max_run_time}<br>
Maximum memory usage: ${max_memory_usage}<br>
Threads: ${threads}<br>
<p>
EOF
}

function report_footer {
    local report_file="$1"
    cat >> "${report_file}" <<EOF
</body>

</html> 
EOF
}

function scores_table_header {
    local report_file="$1"
    local min_range="$2"
    local max_range="$3"

cat >> "${report_file}" <<EOF
<span style="font-size: larger;"><b>Scores:</b></span><br>
<table>
<tr>
<th align="right">Support threshold</th>
EOF

    for s in $(seq $min_range $max_range)
    do
        echo "<th>${s}</th>" >> "${report_file}"
    done

    if [[ "${use_patterns}" -eq "1" ]]
    then
        echo "<th>Best Score</th>" >> "${report_file}"
        echo "<th>Strictly conserved overlap</th>" >> "${report_file}"
    fi

    echo "</tr>" >> "${report_file}"
}

### result overview ############################################################

function parse_result_csv {
    local overview_csv="$1"

    return_best_pattern_score_percent='0.0'
    return_best_pattern_score_binary='0.0'
    return_strictly_conserved=""
    return_strictly_conserved_matched=""
    return_style=""

    local first_row=1
    local result_significant=0
    local result_interesting=0
    local result_first_interesting=0

    if [[ -f "${overview_csv}" ]]
    then
        while IFS='|' read 'id' 'name' 'size' 'vertices' 'support' 'score' 'interesting' 'significant' 'tp' 'fn' 'fp' 'tn' 'mcc' 'odds_ratio' 'p_value' 'pattern_score_percent' 'pattern_score_binary' 'pattern_exceed' 'pattern_conserved' 'pattern_conserved_matched' 'rmsd' 'core' 'lig' 'prot' 'surf' 'dna' 'rna' 'metal' 'ion'
        do
            if [[ $first_row == 1 ]]
            then
                if [[ "${interesting}" == "True" ]]
                then
                    result_first_interesting=1
                fi
                # Only consider values from the best scoring result
                return_strictly_conserved="${pattern_conserved}"
                return_strictly_conserved_matched="${pattern_conserved_matched}"
                first_row=0
            fi
            if [[ "${interesting}" == "True" ]]
            then
                result_interesting=1
            fi
            if [[ "${significant}" == "True" ]]
            then
                result_significant=1
            fi
            if [[ "$(awk -va="${pattern_score_percent}" -vb="${return_best_pattern_score_percent}" "BEGIN{print(a>b)?1:0}")" == "1" ]]
            then
                return_best_pattern_score_percent="${pattern_score_percent}"
            fi
            if [[ "$(awk -va="${pattern_score_binary}" -vb="${return_best_pattern_score_binary}" "BEGIN{print(a>b)?1:0}")" == "1" ]]
            then
                return_best_pattern_score_binary="${pattern_score_binary}"
            fi
        done < "${overview_csv}"
    fi

    local style_classes=()

    if [[ $result_significant == 1 ]]
    then
        style_classes+=("significant")
    fi
    if [[ $result_interesting == 1 ]]
    then
        style_classes+=("interesting")
    fi
    if [[ $result_interesting == 1 && $result_first_interesting == 0 ]]
    then
        style_classes+=("bad-score")
    fi

    if [[ "${#style_classes[@]}" -gt 0 ]]
    then
        return_style="$(join_by ' ' "${style_classes[@]}")"
    fi
}

### output parser ##############################################################

function get_num_found_subgs {
    # Can fail if rinminer crashed
    grep -e 'Found .* supported subgraphs' "$1" | sed -e 's/Found \([0-9]*\) .*/\1/g' || echo "-"
}

function get_max_found_size {
    # Can fail if no subgraphs were found
    found_size="$(grep -e '.*|.*|.*|.*|.*' "$1" | tail -n +2 | wc -l)"
    if [[ "${found_size}" != '0' ]]
    then
        echo "${found_size}"
    else
        echo "-"
    fi
}

function get_best_score_and_size_and_mappings {
    # Can fail if no subgraphs were found
    grep -e '.*|.*|.*|.*|.*' "$1" | tail -n +2 | sed -e 's/\s*\([0-9]*\)\s*|.*|.*|.*|\s*\([0-9]*\) ([0-9]*[^0-9]*\([0-9]*\).*).*/\2 \1 \3/g' | sort -n | tail -n 1 || echo "- - -"
}

### main script ################################################################

mkdir -p "${out_path}"
mkdir -p "${report_path}"

report_file="${report_path}/report.html"
sig_report_file="${report_path}/significance.html"
mcc_report_file="${report_path}/mcc.html"
rmsd_report_file="${report_path}/rmsd.html"
skipped_file="${out_path}/skipped_commands.txt"
pattern_score_distribution_file_percent="${report_path}/pattern_score_distribution_percent.svg"
pattern_score_distribution_file_binary="${report_path}/pattern_score_distribution_binary.svg"

rm -f "${skipped_file}"

skipped=0

# Prepare command list
while read -d '' -r db_file;
do
    num_entries=$(grep -e "^t .*" "${db_file}" | wc -l)

    if [[ "$num_entries" -lt "$min_entries" ]]
    then
        echo "Skipping \"${db_file}\": Not enough entries (less than ${min_entries})" | tee -a "${skipped_file}"
        skipped=1
        continue
    fi

    fam_dir="$(dirname "${db_file}")"
    fam_id="$(basename "${fam_dir}")"
    fam_name="$(cat "$(dirname "${db_file}")/family_name.txt" 2>/dev/null || echo "${fam_id}")"
    fam_pattern=""
    pattern_file="${fam_dir}/pattern.txt"
    if [[ "${use_patterns}" -eq "1" && -f "${pattern_file}" ]]
    then
        fam_pattern="$(cat "${pattern_file}")"
    fi

    support_min=""
    support_max=""

    support_range_file="${fam_dir}/support_range.txt"
    if [[ -f "${support_range_file}" ]]
    then
        IFS=',;-' read support_min support_max < "${support_range_file}"
    fi

    if [[ -z "$support_min" ]]
    then
        let support_min=$num_entries*$support_percentage_start/100
    fi

    if [[ "$support_min" == "_" || "$support_min" -lt "$min_entries" ]]
    then
        support_min="$min_entries"
    fi

    if [[ -z "$support_max" || "$support_max" == "_" || "$support_max" -gt "$num_entries" ]]
    then
        let support_max=$num_entries*$support_percentage_end/100
    fi

    family_dirs+=("$fam_dir")
    family_ids+=("$fam_id")
    family_names+=("$fam_name")
    family_patterns+=("$fam_pattern")

    range_start+=("${support_min}")
    range_end+=("${support_max}")

    for support in $(seq $support_min $support_max)
    do
        out_name="${fam_id}_sup_${support}"

        pretty_name="${fam_name} (${fam_id}), support threshold: ${support}"
        results_file="${out_path}/${out_name}_results"
        tmp_results_file="${tmp_path}/${out_name}_results"
        output_file="${out_path}/${out_name}_output.txt"
        tmp_output_file="${tmp_path}/${out_name}_output.txt"
        error_file="${out_path}/${out_name}_error.txt"
        select_file="${out_path}/${out_name}_selected.txt"
        index_file="${out_path}/${out_name}_index.txt"
        graph_dir="${report_path}/graphs/${out_name}_graphs"
        if [ -n "${significance_path}" ]
        then
            command="${rinminer_path} -o ${min_output_size} -f ${max_output_file_size} -t ${threads} -i ${edge_inclusion} -s ${min_score_difference} -c ${fam_id} -l ${significance_path} ${support} ${db_file} ${tmp_results_file}"
        else
            command="${rinminer_path} -o ${min_output_size} -f ${max_output_file_size} -t ${threads} -i ${edge_inclusion} -s ${min_score_difference} ${support} ${db_file} ${tmp_results_file}"
        fi

        percentage_min=$(echo "y=100.0/$increment; perc=$support/$support_max; rem=$support%$support_max; scale=0; x=perc*y; x=x/1; if (rem) { x=x+1; }; x=x*100/y; x=x/1; if (100.0*perc-x<0.0000001) { x; } else { x+$increment; }" | bc -l)
        percentage_max=$(echo "y=100.0/$increment; perc=($support+1)/$support_max; rem=($support+1)%$support_max; scale=0; x=perc*y; x=x/1; if (rem) { x=x+1; }; x=x*100/y; x=x/1; if (100.0*perc-x<0.0000001) { x=x-$increment; }; if (x>100) { x=100; }; x"| bc -l)
        tag="$(echo $(seq $percentage_min $increment $percentage_max | sed -e 's/$/-sup/g') | tr ' ' ',')"

        pretty_names+=("$pretty_name")
        results_files+=("$results_file")
        tmp_results_files+=("$tmp_results_file")
        output_files+=("$output_file")
        tmp_output_files+=("$tmp_output_file")
        error_files+=("$error_file")
        select_files+=("$select_file")
        index_files+=("$index_file")
        graph_dirs+=("$graph_dir")
        commands+=("$command")
        tags+=("$tag")
    done
done < <(find -L "${data_path}" -mindepth 2 -maxdepth 2 -type f -name database -print0 | sort -z -V)

# base-threshold / number of runs / max number of results per run
significance_threshold_string="${significance_base_threshold} / ${#commands[@]} / ${max_subgraph_results}"
significance_threshold="$(awk "BEGIN{print ${significance_threshold_string}}")"

# Find min range start and max range end
min_range=99999999
max_range=0
num_ranges=${#range_start[@]}
for (( i=0; i<${num_ranges}; i++ ));
do
    start="${range_start[i]}"
    end="${range_end[i]}"
    if [[ "${start}" -lt "${min_range}" ]]
    then
        min_range="${start}"
    fi
    if [[ "${end}" -gt "${max_range}" ]]
    then
        max_range="${end}"
    fi
done

# Prepare the report
report_header "${report_file}"
scores_table_header "${report_file}" "${min_range}" "${max_range}"

# The command counter
i=0
num_commands="${#commands[@]}"
num_fams=${#family_names[@]}

pattern_scores_percent=()
pattern_scores_binary=()

for (( f=0; f<${num_fams}; f++ ));
do
    family_dir="${family_dirs[f]}"
    family_id="${family_ids[f]}"
    family_name="${family_names[f]}"
    family_pattern="${family_patterns[f]}"

    start="${range_start[f]}"
    end="${range_end[f]}"

    db_size="${end}"

    echo "<tr>" >> "${report_file}"
    if [[ "${family_name}" == "${family_id}" ]]
    then
        row_title="${family_name}"
    else
        row_title="${family_name} (${family_id})"
    fi
    echo "<th align=\"left\">${row_title}</th>" >> "${report_file}"

    # Report pre padding
    for (( x=min_range; x<${start}; x++ ));
    do
        echo "<td class=\"empty\"></td>" >> "${report_file}"
    done

    best_row_pattern_score_percent=""
    best_row_pattern_score_binary=""
    best_row_strictly_conserved=""
    best_row_strictly_conserved_matched=""
    best_row_score_found=0

    # Run the commands for this family and add them to the report
    for (( x=start; x<=${end}; x++ ));
    do
        command="${commands[i]}"
        tag="${tags[i]}"
        pretty_name="${pretty_names[i]}"
        results_file="${results_files[i]}"
        tmp_results_file="${tmp_results_files[i]}"
        output_file="${output_files[i]}"
        tmp_output_file="${tmp_output_files[i]}"
        error_file="${error_files[i]}"
        select_file="${select_files[i]}"
        index_file="${index_files[i]}"
        graph_dir="${graph_dirs[i]}"
        unique_prefix="${family_id}_sup_${x}_"
        overview_csv="${graph_dir}/overview.csv"
        overview_html="${graph_dir}/overview.html"

        num_found_subgs=0
        max_found_size=0
        best_score=""
        best_score_size=""
        best_score_mappings=""

        skip_reason=""
        exit_reason=""

        rerun_eval_needed="0"

        let current=${i}+1
        echo "($current/$num_commands) $pretty_name"

        if [[ -f "${error_file}" ]]
        then
            exit_reason="$(cat "${error_file}")"
            if [[ "${rerun_crashed}" -ne "1" ]]
            then
                skip_reason="previously run and crashed: ${exit_reason}"
            fi
        else
            if [[ -f "${results_file}" &&  "${rerun}" -ne "1" ]]
            then
                skip_reason="previously run"
            fi
        fi


        if [[ -n "${skip_reason}" ]]
        then
            echo "Skipped: ${skip_reason}"
        else
            # Run the command, it is allowed to fail, but fails are recorded
            echo "$command" | tee "${tmp_output_file}"

            exit_status=("0" "0")
            timeout $max_run_time $command 2>&1 | tee -a "${tmp_output_file}" || exit_status=("${PIPESTATUS[@]}") && true
            exit_reason=""

            # manually killed, stop the script        
            if [[ "${exit_status[0]}" -eq "130" || "${exit_status[1]}" -eq "130" ]]
            then
                exit 130
            # segfault from dereferencing NULL or kill signal from OOM killer
            elif [[ "${exit_status[0]}" -eq "139" || "${exit_status[0]}" -eq "137" ]]
            then
                exit_reason="out of memory"
            # timeouts only affect the timeout exit status
            elif [[ "${exit_status[0]}" -eq "124" ]]
            then
                exit_reason="timeout"
            # unknown error exist status
            elif [[ "${exit_status[0]}" -ne "0" ]]
            then
                exit_reason="unknown (${exit_status[0]})"
            fi

            if [[ -n "${exit_reason}" ]]
            then
                echo "FAILED: ${exit_reason}" | tee -a "${tmp_output_file}"
                echo "${exit_reason}" > "${error_file}"
                rm -f "${tmp_results_file}"
            else
                rm -f "${error_file}"
                mv -f "${tmp_results_file}" "${results_file}"
            fi
            mv -f "${tmp_output_file}" "${output_file}"

            rerun_eval_needed="1"
        fi

        # Select best subgraphs if required
        if [[ -f "${results_file}" && ( "${rerun_select}" == "1" || ! -f "${select_file}" || "${results_file}" -nt "${select_file}" || ! -f "${index_file}" || "${results_file}" -nt "${index_file}") ]]
        then
            options=('-f' '-b' '-l' "${max_subgraph_results}")
            "${scripts_path}/analyze-results/select-graphs.py" "${options[@]}" -- "${results_file}" "${index_file}" "${select_file}"
            rerun_eval_needed="1"
        fi

        num_found_subgs="$(get_num_found_subgs "${output_file}")"
        max_found_size="$(get_max_found_size "${output_file}")"
        best_score_and_size_and_mappings="$(get_best_score_and_size_and_mappings "${output_file}")"
        best_score="$(echo "${best_score_and_size_and_mappings}" | cut -d\  -f 1)"
        best_score_size="$(echo "${best_score_and_size_and_mappings}" | cut -d\  -f 2)"
        best_score_mappings="$(echo "${best_score_and_size_and_mappings}" | cut -d\  -f 3)"

        # Generate the svg-, pml-, csv- and the smlf-files and the overview for this result
        if [[ ! -f "${overview_csv}" || "${rerun_eval_needed}" == "1" || "${rerun_eval}" == "1" ]]
        then
            rm -rf "${graph_dir}"
            mkdir -p "${graph_dir}"

            selected_graphs=()
            if [[ -f "${select_file}" ]]
            then
                readarray -t selected_graphs < "${select_file}"
            fi

            options=('-f' '-b' '-i' '-p' '-c' '-y' "${family_dir}"  '-d' "${min_score_difference}" '-t' "${pretty_name}" '-o' "${output_file}" '-x' "${unique_prefix}" '-l' "${family_id}")

            if [[ "${use_significance}" -eq "1" ]]
            then
                options+=('-n' "${significance_threshold}")
            fi

            if [[ "${use_patterns}" -eq "1" ]]
            then
                options+=('-a' "${family_pattern}")
            fi

            if [[ "${use_smlf}" -eq "1" ]]
            then
                options+=('-m')
            fi

            if [[ "${use_sm_classes}" -eq "1" ]]
            then
                options+=('-w' "${sm_classes_file}")
            fi

            if [[ "${best_row_score_found}" -ne "1" ]]
            then
                options+=('-z')
            fi

            if [[ -n "$tag" ]]
            then
                options+=('-u' "$tag")
            fi

            "${scripts_path}/analyze-results/analyze-graphs.py" "${options[@]}" -- "${results_file}" "${index_file}" "${graph_dir}" "${selected_graphs[@]}"
        fi

        parse_result_csv "${overview_csv}"
        style="${return_style}"

        # Ignore pattern scores from timeouts/crashes
        if [[ -z "${exit_reason}" ]]
        then
            best_pattern_score_percent="${return_best_pattern_score_percent}"
            best_pattern_score_binary="${return_best_pattern_score_binary}"
        else
            best_pattern_score_percent=''
            best_pattern_score_binary=''
        fi


        if [[ -z "${best_row_pattern_score_percent}" || "$(awk -va="${best_pattern_score_percent}" -vb="${best_row_pattern_score_percent}" "BEGIN{print(a>b)?1:0}")" == "1" ]]
        then
            best_row_pattern_score_percent="${best_pattern_score_percent}"
        fi

        if [[ -z "${best_row_pattern_score_binary}" || "$(awk -va="${best_pattern_score_binary}" -vb="${best_row_pattern_score_binary}" "BEGIN{print(a>b)?1:0}")" == "1" ]]
        then
            best_row_pattern_score_binary="${best_pattern_score_binary}"
        fi

        # Only consider the first column with an entry, those are the largest subgraphs and therefore the highest scoring ones
        if [[ -z "${best_row_strictly_conserved}" && -n ${return_strictly_conserved} ]]
        then
            best_row_strictly_conserved=${return_strictly_conserved}
            best_row_strictly_conserved_matched=${return_strictly_conserved_matched}
        fi

        # Write the report data
        pretty_name_quoted="$(echo "${pretty_name}" | sed -e "s/\"/\&quot;/g")"
        title="${pretty_name_quoted}
Best score size: ${best_score_size}, best score support: ${best_score_mappings}, total subgraps: ${num_found_subgs}, maximum size: ${max_found_size}"

        text="${best_score}"
        if [[ "${use_patterns}" -eq "1" ]]
        then
            if [[ -n "${best_pattern_score_percent}" ]]
            then
                text="$(printf "%.3f" "${best_pattern_score_percent}")"
            else
                text="-"
            fi

            if [[ -n "${best_pattern_score_binary}" ]]
            then
                text="${text} / $(printf "%.3f" "${best_pattern_score_binary}")"
            else
                text="${text} / -"
            fi
        fi

        if [[ ${max_found_size} == '-' || ${max_found_size} -lt ${min_output_size} ]]
        then
            style="no-graphs"
        else
            best_row_score_found=1
        fi

        if [[ -n "${exit_reason}" ]]
        then
            style="failed"
            text="${exit_reason}"
        fi

        style_html=""
        if [[ -n "${style}" ]]
        then
            style_html=" class=\"${style}\""
        fi

        overview_html="$(realpath --relative-to="${report_path}" "${overview_html}")"
        cell_content="<a href=\"${overview_html}\"><span title=\"${title}\">${text}</span></a>"

        echo "<td${style_html}>${cell_content}</td>" >> "${report_file}"

        let i=i+1

        echo "-------------------------------------------------------------------------------------------"
    done

    # Report post padding
    for (( x=end+1; x<=${max_range}; x++ ));
    do
        echo "<td class=\"empty\"></td>" >> "${report_file}"
    done

    # Best row pattern score in pattern mode
    if [[ "${use_patterns}" -eq "1" ]]
    then
        if [[ -n "${best_row_pattern_score_percent}" ]]
        then
            text="$(printf "%.3f" "${best_row_pattern_score_percent}")"
            pattern_scores_percent+=("${best_row_pattern_score_percent}")
        else
            text="-"
        fi

        if [[ -n "${best_row_pattern_score_binary}" ]]
        then
            text="${text} / $(printf "%.3f" "${best_row_pattern_score_binary}")"
            pattern_scores_binary+=("${best_row_pattern_score_binary}")
        else
            text="${text} / -"
        fi
        echo "<td>${text}</td>" >> "${report_file}"


        text="- / -"
        if [[ -n ${best_row_strictly_conserved} ]]
        then
            text="${best_row_strictly_conserved_matched} / ${best_row_strictly_conserved}"
        fi
        echo "<td>${text}</td>" >> "${report_file}"
    fi

    echo "</tr>" >> "${report_file}"
done

echo "</table>" >> "${report_file}"

cat >> "${report_file}" <<EOF
<p>
<span style="font-size: larger;"><b>Legend:</b></span><br>
<table>
EOF

if [[ "${use_patterns}" -eq "1" ]]
then
    echo "<tr><td>Percent score / Binary score</td></tr>" >> "${report_file}"
fi

cat >> "${report_file}" <<EOF
<tr><td class="interesting">interesting</td></tr>
<tr><td class="interesting bad-score">interesting, but not highest score</td></tr>
<tr><td class="significant">significant</td></tr>
<tr><td class="no-graphs">no graphs within threshold</td></tr>
<tr><td class="failed">failed</td></tr>
<tr><td class="empty">not tested</td></tr>
</table>
EOF

# Add score distribution plots for pattern mode
function plot_score_distribution {
    local svg_file="$1"
    shift
    local title="$1"
    shift
    local scores=($*)
    Rscript - <<EOF
v <- c($(join_by ', ' "${scores[@]}"))
svg('${svg_file}')
h <- hist(v, probability=TRUE, col='blue', xlab='Score', main='${title}')
dev.off()
EOF
}

if [[ "${use_patterns}" -eq "1" ]]
then
    plot_score_distribution "${pattern_score_distribution_file_percent}" "Percentage score distribution" ${pattern_scores_percent[@]}
    pattern_score_distribution_file_percent="$(realpath --relative-to="${report_path}" "${pattern_score_distribution_file_percent}")"
    echo "<img src=\"${pattern_score_distribution_file_percent}\">" >> "${report_file}"

    plot_score_distribution "${pattern_score_distribution_file_binary}" "Binary score distribution" ${pattern_scores_binary[@]}
    pattern_score_distribution_file_binary="$(realpath --relative-to="${report_path}" "${pattern_score_distribution_file_binary}")"
    echo "<img src=\"${pattern_score_distribution_file_binary}\">" >> "${report_file}"
fi

# Finalize the report
report_footer "${report_file}"

# Combine all SMLF files
if [[ "${use_smlf}" -eq "1" ]]
then
    combined_smlf="${report_path}/all_mapped_residues.smlf"
    rm -f "${combined_smlf}"
    touch "${combined_smlf}"
    combined_random_smlf="${report_path}/all_mapped_residues_random.smlf"
    rm -f "${combined_random_smlf}"
    touch "${combined_random_smlf}"
    for (( i=0; i<${num_commands}; i++ ));
    do
        graph_dir="${graph_dirs[i]}"
        overview_smlf="${graph_dir}/overview.smlf"
        if [[ -f "${overview_smlf}" ]]
        then
            cat "${overview_smlf}" >> "${combined_smlf}"
        fi
        overview_random_smlf="${graph_dir}/overview_random.smlf"
        if [[ -f "${overview_random_smlf}" ]]
        then
            cat "${overview_random_smlf}" >> "${combined_random_smlf}"
        fi
    done
fi

### significance report ########################################################

function stat_report {
    local report_file="$1"
    local stat_csv_file="$2"
    local title="$3"

    report_header "${report_file}"

    cat >> "${report_file}" <<EOF
<span style="font-size: larger;"><b>${title}:</b></span><br>
<table>
<tr>
<th>Dataset</th>
<th>Motif</th>
<th>Size</th>
<th>Score</th>
<th>RMSD</th>
<th>TP</th>
<th>FN</th>
<th>FP</th>
<th>TN</th>
<th>MCC</th>
<th>Odds Ratio</th>
<th>p-Value</th>
EOF

    if [[ "${use_sm_classes}" -eq "1" ]]
    then
        echo "<th>Surface</th>" >> "${report_file}"
        echo "<th>Core</th>" >> "${report_file}"
        echo "<th>DNA</th>" >> "${report_file}"
        echo "<th>RNA</th>" >> "${report_file}"
        echo "<th>Ligand</th>" >> "${report_file}"
        echo "<th>Protein</th>" >> "${report_file}"
        echo "<th>Metal</th>" >> "${report_file}"
        echo "<th>Ion</th>" >> "${report_file}"
    fi

    echo "</tr>" >> "${report_file}"

    prev_name=""
    #                   1    2      3      4          5         6       7             8             9    10   11   12   13    14           15        16                      17                     18               19                  20                          21     22     23    24     25     26    27    28      29    30       31            32
    while IFS='|' read 'id' 'name' 'size' 'vertices' 'support' 'score' 'interesting' 'significant' 'tp' 'fn' 'fp' 'tn' 'mcc' 'odds_ratio' 'p_value' 'pattern_score_percent' 'pattern_score_binary' 'pattern_exceed' 'pattern_conserved' 'pattern_conserved_matched' 'rmsd' 'core' 'lig' 'prot' 'surf' 'dna' 'rna' 'metal' 'ion' 'number' 'pretty_name' 'graph_dir'
    do
        rel_graph_dir="$(realpath --relative-to="${report_path}" "${graph_dir}")"
        overview_html="${rel_graph_dir}/overview.html"
        pml_file="${rel_graph_dir}/${name}.pml"
        svg_file="${rel_graph_dir}/${name}.svg"

        # Hide repeating pretty names
        if [[ ${pretty_name} == ${prev_name} ]]
        then
            dataset="..."
        else
            prev_name="${pretty_name}"
            dataset="<a href=\"${overview_html}\">${pretty_name}</a>"
        fi

        style=()
        if [[ ${significant} == "True" ]]
        then
            style+=("significant")
        fi
        if [[ ${interesting} == "True" ]]
        then
            style+=("interesting")
            if [[ ${number} -ne 1 ]]
            then
                style+=("bad-score")
            fi
        fi

        style_string=""
        if [[ ${#style[@]} -gt 0 ]]
        then
            style_string=" class=\"$(join_by ' ' "${style[@]}")\""
        fi

        cat >> "${report_file}" <<EOF
<tr>
<td style="text-align: left;">${dataset}</td>
<td><a href="${svg_file}">SVG</a> <a href="${pml_file}">PyMol</a></td>
<td>${size}</td>
<td>${score}</td>
<td>$(printf "%.5f" "${rmsd}")</td>
<td>${tp}</td>
<td>${fn}</td>
<td>${fp}</td>
<td>${tn}</td>
<td>$(printf "%.5f" "${mcc}")</td>
<td>$(printf "%.3f" "${odds_ratio}")</td>
<td${style_string}>$(printf "%.3e" "${p_value}")</td>
EOF

        if [[ "${use_sm_classes}" -eq "1" ]]
        then
            echo "<td>$(printf "%.3lf" "${surf}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${core}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${dna}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${rna}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${lig}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${prot}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${metal}")</td>" >> "${report_file}"
            echo "<td>$(printf "%.3lf" "${ion}")</td>" >> "${report_file}"
        fi

        echo "</tr>" >> "${report_file}"

    done < "${stat_csv_file}"

    echo "</table>" >> "${report_file}"

    cat >> "${report_file}" <<EOF
<p>
<span style="font-size: larger;"><b>Legend:</b></span><br>
<table>
<tr><td class="interesting">interesting</td></tr>
<tr><td class="interesting bad-score">interesting, but not best score</td></tr>
<tr><td class="significant">significant</td></tr>
</table>
EOF

    report_footer "${report_file}"
}

tmp_csv="${tmp_path}/tmp.csv"
rm -f "${tmp_csv}"
touch "${tmp_csv}"

for (( i=0; i<${num_commands}; i++ ));
do
    pretty_name="${pretty_names[i]}"
    graph_dir="${graph_dirs[i]}"
    overview_csv="${graph_dir}/overview.csv"
    if [[ -f "${overview_csv}" ]]
    then
        number=0
        while read line
        do
            let number=$number+1
            echo -n "${line}|" >>"${tmp_csv}"
            echo "${number}|${pretty_name}|${graph_dir}" >>"${tmp_csv}"
        done < "${overview_csv}"
    fi
done

# sort by rmsd (asc), score (desc) and pretty name (asc)
rmsd_csv="${tmp_path}/rmsd.csv"
sort -t '|' -k 21,21g -k 6,6r -k31,31 --parallel=4 -o "${rmsd_csv}" "${tmp_csv}"
stat_report "${rmsd_report_file}" "${rmsd_csv}" "Subgraphs sorted by RMSD"

# Generate reports sorted by structman classification
if [[ "${use_sm_classes}" -eq "1" ]]
then
    # sort by core (desc), score (desc) and pretty name (asc)
    sm_report_file="${report_path}/sm-core.html"
    sm_csv="${tmp_path}/sm.csv"
    sort -t '|' -k 22,22rg -k 6,6r -k 31,31 --parallel=4 -o "${sm_csv}" "${tmp_csv}"
    stat_report "${sm_report_file}" "${sm_csv}" "Subgraphs sorted by 'core' percentage"

    # sort by ligand (desc), score (desc) and pretty name (asc)
    sm_report_file="${report_path}/sm-lig.html"
    sm_csv="${tmp_path}/sm.csv"
    sort -t '|' -k 23,23rg -k 6,6r -k 31,31 --parallel=4 -o "${sm_csv}" "${tmp_csv}"
    stat_report "${sm_report_file}" "${sm_csv}" "Subgraphs sorted by 'ligand' percentage"

    # sort by protein (desc), score (desc) and pretty name (asc)
    sm_report_file="${report_path}/sm-prot.html"
    sm_csv="${tmp_path}/sm.csv"
    sort -t '|' -k 24,24rg -k 6,6r -k 31,31 --parallel=4 -o "${sm_csv}" "${tmp_csv}"
    stat_report "${sm_report_file}" "${sm_csv}" "Subgraphs sorted by 'protein' percentage"

    # sort by surface (desc), score (desc) and pretty name (asc)
    sm_report_file="${report_path}/sm-surf.html"
    sm_csv="${tmp_path}/sm.csv"
    sort -t '|' -k 25,25rg -k 6,6r -k 31,31 --parallel=4 -o "${sm_csv}" "${tmp_csv}"
    stat_report "${sm_report_file}" "${sm_csv}" "Subgraphs sorted by 'surface' percentage"

    # sort by dna (desc), score (desc) and pretty name (asc)
    sm_report_file="${report_path}/sm-dna.html"
    sm_csv="${tmp_path}/sm.csv"
    sort -t '|' -k 26,26rg -k 6,6r -k 31,31 --parallel=4 -o "${sm_csv}" "${tmp_csv}"
    stat_report "${sm_report_file}" "${sm_csv}" "Subgraphs sorted by 'dna' percentage"

    # sort by rna (desc), score (desc) and pretty name (asc)
    sm_report_file="${report_path}/sm-rna.html"
    sm_csv="${tmp_path}/sm.csv"
    sort -t '|' -k 27,27rg -k 6,6r -k 31,31 --parallel=4 -o "${sm_csv}" "${tmp_csv}"
    stat_report "${sm_report_file}" "${sm_csv}" "Subgraphs sorted by 'rna' percentage"
fi

# Generate significance report sorted by p-value and MCC
if [[ "${use_significance}" -eq "1" ]]
then
    # sort by p-value (asc), score (desc) and pretty name (asc)
    sig_csv="${tmp_path}/sig.csv"
    sort -t '|' -k 15,15g -k 6,6r -k 31,31 --parallel=4 -o "${sig_csv}" "${tmp_csv}"
    stat_report "${sig_report_file}" "${sig_csv}" "Subgraphs sorted by p-value"

    # sort by mcc (desc), score (desc) and pretty name (asc)
    mcc_csv="${tmp_path}/mcc.csv"
    sort -t '|' -k 13,13gr -k 6,6r -k 31,31 --parallel=4 -o "${mcc_csv}" "${tmp_csv}"
    stat_report "${mcc_report_file}" "${mcc_csv}" "Subgraphs sorted by MCC"
fi
