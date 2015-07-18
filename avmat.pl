#!/usr/bin/perl

$total = 0.0;
$nvalues = 0;
while(<>)
{
    chomp;
    @values = split;
    shift @values;
    foreach $value (@values)
    {
        $total += $value;
        $nvalues++;
    }
}

print "Total = $total  NValues = $nvalues\n";
printf "Mean = %f\n", $total / $nvalues;


