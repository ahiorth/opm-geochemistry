function sample_key(species, temperature)
{
    return species SUBSEP sprintf("%.2f", temperature)
}

FNR == NR {
    if (FNR > 1)
    {
        reference[sample_key($1, $2)] = $3
    }
    next
}

FNR > 1 {
    key = sample_key($1, $2)
    if (key in reference)
    {
        difference = $4 - reference[key]
        absolute_difference = difference < 0 ? -difference : difference
        compared++

        if (absolute_difference > maximum)
        {
            maximum = absolute_difference
            maximum_species = $1
            maximum_temperature = $2
            maximum_signed_difference = difference
        }
    }
}

END {
    if (compared == 0)
    {
        print "No common samples found." > "/dev/stderr"
        exit 1
    }

    printf "compared=%d max_abs_difference=%.7f species=%s T_C=%.2f signed_difference=%+.7f\n", \
        compared, maximum, maximum_species, maximum_temperature, maximum_signed_difference
}
