#!/usr/bin/awk -f
BEGIN { FS = "|" }
{
    if ($1 ~ /^>ENST/) {
        print $1
    } else {
        print
    }
}
