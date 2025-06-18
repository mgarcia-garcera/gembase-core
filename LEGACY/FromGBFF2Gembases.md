    #### GEMBASES from GBFF



- step 1.
- core

```{bash}
awk '

BEGIN {
    FS="\""
    OFS=""
    position=0
}
#extract coordinates
/\s+gene\s+/&&!/\/function=/ {
    position +=1
    c[position] = 1
    if (/complement\(/ ) {
        c[position] = 2
    }
    if(match($0, /[0-9]+..[0-9]+/)) {
        split(substr($0, RSTART, RLENGTH), coords, /\.\./)
        start[position] = coords[1]
        end[position] = coords[2]
    }
}
#extract translations
/\/translation="/ {
    translation[position] = $2
    if ($0 !~ /"$/) {
        while (getline line) {
            if (line ~ /"$/) {
                sub(/".*/, "", line)
                translation[position] = translation[position] line
                break
            }
            translation[position] = translation[position] line
        }
    }
}
#extract gene descriptions
/\/product="/ {
    product[position] = $2
    if ($0 !~ /"$/) {
        while (getline line) {
            if (line ~ /"$/) {
                sub(/".*/, "", line)
                prodcut[position] = product[position] line
                break
            }
            translation[position] = translation[position] line
        }
    }
}
#extract gene names
/\/gene="/ {
    gene[position] = $2
}
END {
    for (i=1; i<=position; i++) {
        gsub(/[\t\n\r]+/,"", translation[i])
        gsub(/[\t\n\r]+/,"", product[i])
        print i" "start[i]" "end[i]" "c[i]" "gene[i]" "product[i]
    }
}
' $file
```