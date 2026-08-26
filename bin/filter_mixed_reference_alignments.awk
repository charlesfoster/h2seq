BEGIN {
    FS = OFS = "\t"
    current_name = ""
    reference_count = 0
}

function flag_is_set(flag, bit) {
    return int(flag / bit) % 2
}

function is_primary_mapped(flag) {
    return !flag_is_set(flag, 4) && !flag_is_set(flag, 256) && !flag_is_set(flag, 2048)
}

function flush_fragment(    i, flag, primary_count, reference, assigned_reference, ambiguous) {
    if (record_count == 0) return

    primary_count = 0
    assigned_reference = ""
    ambiguous = 0
    for (i = 1; i <= record_count; i++) {
        flag = record_flag[i]
        if (!is_primary_mapped(flag)) continue
        primary_count++
        reference = record_reference[i]
        if (assigned_reference == "") assigned_reference = reference
        else if (reference != assigned_reference) ambiguous = 1
        if (reference_count > 1 && record_mapq[i] < min_mapq) ambiguous = 1
    }

    if (primary_count == 0) {
        unassigned_fragments++
        if (reference_count == 1) {
            for (i = 1; i <= record_count; i++) print record_line[i]
        }
    } else if (reference_count > 1 && ambiguous) {
        ambiguous_fragments++
        ambiguous_primary_records += primary_count
    } else {
        assigned_fragments[assigned_reference]++
        assigned_primary_records[assigned_reference] += primary_count
        for (i = 1; i <= record_count; i++) {
            flag = record_flag[i]
            if (reference_count == 1 || (!flag_is_set(flag, 256) && !flag_is_set(flag, 2048))) {
                print record_line[i]
            }
        }
    }

    delete record_line
    delete record_flag
    delete record_reference
    delete record_mapq
    record_count = 0
}

/^@/ {
    print
    if ($1 == "@SQ") {
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) {
                reference_count++
                reference_order[reference_count] = substr($i, 4)
                break
            }
        }
    }
    next
}

{
    query_name = $1
    if (current_name != "" && query_name != current_name) flush_fragment()
    current_name = query_name
    record_count++
    record_line[record_count] = $0
    record_flag[record_count] = $2 + 0
    record_reference[record_count] = $3
    record_mapq[record_count] = $5 + 0
}

END {
    flush_fragment()
    print "sample_id", "read_type", "assignment", "reference_name", "fragments", "primary_alignment_records", "min_mapq", "reference_count" > summary_output
    for (i = 1; i <= reference_count; i++) {
        reference = reference_order[i]
        if (assigned_fragments[reference] > 0) {
            print sample_id, read_type, "assigned", reference, assigned_fragments[reference], assigned_primary_records[reference], min_mapq, reference_count > summary_output
        }
    }
    if (ambiguous_fragments > 0) {
        print sample_id, read_type, "ambiguous", "", ambiguous_fragments, ambiguous_primary_records, min_mapq, reference_count > summary_output
    }
    if (unassigned_fragments > 0) {
        print sample_id, read_type, "unassigned", "", unassigned_fragments, 0, min_mapq, reference_count > summary_output
    }
}
