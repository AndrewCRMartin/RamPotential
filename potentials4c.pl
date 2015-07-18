#!/usr/bin/perl -d

#Perl program to find the potentials between all pairs of amino acid atoms in a PDB file

$pdbfile = $ARGV[0];
$counta = 0;
$TotalEnergy = 0;

#Open the PDB file
open(PFILE, "$pdbfile") || die "Can't read $pdbfile";


#Opens the potentials datafile and reads it into an array
ReadDatafile();


#Getting coordinates from PDB file
GetCoordinates();


#for each pair of amino acid atom groups
for($i = 0; $i < $counta; $i++)
{
    #compare with all amino acids atom groups that come after it
    for ($j = $i + 1; $j < $counta; $j++)
    {
	#Using pythag to calculate the distance between the two atom groups
	$distSquared = ($x[$i] - $x[$j]) * ($x[$i] - $x[$j]) + 
	               ($y[$i] - $y[$j]) * ($y[$i] - $y[$j]) + 
		       ($z[$i] - $z[$j]) * ($z[$i] - $z[$j]);
	{
	    $distance = sqrt $distSquared;
		
	    
	    #Uses the $position and $matrix values to find the correct potential for the pair from the datafile and print it out
	    $potential = MatrixLine();
	    
            $TotalEnergy = $TotalEnergy + $potential;
	}	
    }
}

print "Total energy potential = ";
print "$TotalEnergy";



sub GetCoordinates
#for the length of the pdb file
{
    while(<PFILE>)
    {
	($type, $atom, $atomname, $resname) = split;
	if ($type eq "ATOM")
	{			
	    #finds the amino acid type, the group type and the coordiantes of the group
	    $resid[$counta] = substr ($_, 17, 3);
	    $atomtype[$counta] = substr ($_, 13, 3);
	    $x[$counta] = substr ($_, 30, 8);
	    $y[$counta] = substr ($_, 38, 8);
	    $z[$counta] = substr ($_, 46, 8);	    
	 
            #Converting residue id and atom type into a number that can be used to pull the potential out of the datafile
	    $position[$counta] = Convert();

	    $counta++;
	}
    }
}


sub MatrixNumber
#Take distance worked out from Pythag
#Use distance to work out which matrix to use
{
    my($distance) = @_;

    if ($distance < 3.0) 
    {
	$matrix = 0;
    }	
    elsif ($distance < 4.0) 
    {
	$matrix = 1;
    }
    elsif ($distance < 5.0) 
    {
	$matrix = 2;
    }
    elsif ($distance < 6.0) 
    {
	$matrix = 3;
    }
    elsif ($distance < 7.0) 
    {
	$matrix = 4;
    }
    elsif ($distance < 8.0) 
    {
	$matrix = 5;
    }
    elsif ($distance < 9.0) 
    {
	$matrix = 6;
    }
    elsif ($distance < 10.0) 
    {
	$matrix = 7;
    }
    elsif ($distance < 11.0) 
    {	
	$matrix = 8;
    }
    elsif ($distance < 12.0) 
    {
	$matrix = 9;
    }
    elsif ($distance < 13.0) 
    {
	$matrix = 10;
    }
    elsif ($distance < 14.0) 
    {
	$matrix = 11;
    }
    elsif ($distance < 15.0) 
    {
	$matrix = 12;
    }
    elsif ($distance < 16.0) 
    {
	$matrix = 13;
    }
    elsif ($distance < 17.0) 
    {
	$matrix = 14;
    }
    elsif ($distance < 18.0) 
    {
	$matrix = 15;
    }
    elsif ($distance < 19.0) 
    {
	$matrix = 16;
    }
    elsif ($distance < 20.0) 
    {
	$matrix = 17;
    }
    elsif ($distance > 20.0)
    {
	$matrix = "distance too great";
    }
    return($matrix);
}

sub Convert
#changing the residue type and atome type to find which matrix to use and position within the matrix
{
     if ($resid[$counta] eq "ALA")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 1;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 2;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 3;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 4;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 5;
	 }
     }
     elsif ($resid[$counta] eq "CYS")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {	
	     $position[$counta] = 6;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 7;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 8;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 9;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 10;
	 }
	 elsif ($atomtype[$counta] eq "SG ")
	 {
	     $position[$counta] = 11;
	 }
     }
     elsif ($resid[$counta] eq "ASP")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 12;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 13;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 14;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 15;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 16;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 17;
	 }
	 elsif ($atomtype[$counta] eq "OD1")
	 {
	     $position[$counta] = 18;
	 }
	 elsif ($atomtype[$counta] eq "OD2")
	 {
	     $position[$counta] = 19;
	 }
     }
     elsif ($resid[$counta] eq "GLU")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 20;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 21;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 22;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 23;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 24;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 25;
	 }
	 elsif ($atomtype[$counta] eq "CD ")
	 {
	     $position[$counta] = 26;
	 }
	 elsif ($atomtype[$counta] eq "OE1")
	 {
	     $position[$counta] = 27;
	 }
	 elsif ($atomtype[$counta] eq "OE2")
	 {
	     $position[$counta] = 28;
	 }
     }
     elsif ($resid[$counta] eq "PHE")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 29;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 30;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 31;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 32;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 33;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 34;
	 }
	 elsif ($atomtype[$counta] eq "CD1")
	 {
	     $position[$counta] = 35;
	 }
	 elsif ($atomtype[$counta] eq "CD2")
	 {
	     $position[$counta] = 36;
	 }
	 elsif ($atomtype[$counta] eq "CE1")
	 {
	     $position[$counta] = 37;
	 }
	 elsif ($atomtype[$counta] eq "CE2")
	 {
	     $position[$counta] = 38;
	 }
	 elsif ($atomtype[$counta] eq "CZ ")
	 {
	     $position[$counta] = 39;
	 }
     }
     elsif ($resid[$counta] eq "GLY")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 40;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 41;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 42;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 43;
	 }
     }
     elsif ($resid[$counta] eq "HIS")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 44;
	 } 
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 45;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 46;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 47;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 48;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 49;
	 }
	 elsif ($atomtype[$counta] eq "ND1")
	 {
	     $position[$counta] = 50;
	 }
	 elsif ($atomtype[$counta] eq "CD2")
	 {
	     $position[$counta] = 51;
	 }
	 elsif ($atomtype[$counta] eq "CE1")
	 {
	     $position[$counta] = 52;
	 }
	 elsif ($atomtype[$counta] eq "NE2")
	 {
	     $position[$counta] = 53;
	 }
     }
     elsif ($resid[$counta] eq "ILE")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 54;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 55;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 56;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 57;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 58;
	 }
	 elsif ($atomtype[$counta] eq "CG1")
	 {
	     $position[$counta] = 59;
	 }
	 elsif ($atomtype[$counta] eq "CG2")
	 {
	     $position[$counta] = 60;
	 }
	 elsif ($atomtype[$counta] eq "CD1")
	 {
	     $position[$counta] = 61;
	 }
     }
     elsif ($resid[$counta] eq "LYS")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 62;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 63;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 64;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 65;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 66;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 67;
	 }
	 elsif ($atomtype[$counta] eq "CD ")
	 {
	     $position[$counta] = 68;
	 }
	 elsif ($atomtype[$counta] eq "CE ")
	 {
	     $position[$counta] = 69;
	 }	
	 elsif ($atomtype[$counta] eq "NZ ")
	 {
	     $position[$counta] = 70;
	 }
	 elsif ($atomtype[$counta] eq "CZ ")
	 {
	     $position[$counta] = 70;
	 }
     }
     elsif ($resid[$counta] eq "LEU")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 71;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 72;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 73;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 74;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 75;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 76;
	 }
	 elsif ($atomtype[$counta] eq "CD1")
	 {
	     $position[$counta] = 77;
	 }
	 elsif ($atomtype[$counta] eq "CD2")
	 {
	     $position[$counta] = 78;
	 }
     }
     elsif ($resid[$counta] eq "MET")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 79;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 80;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 81;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 82;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 83;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 84;
	 }
	 elsif ($atomtype[$counta] eq "SD ")
	 {
	     $position[$counta] = 85;
	 }
	 elsif ($atomtype[$counta] eq "CE ")
	 {
	     $position[$counta] = 86;
	 }
     }
     elsif ($resid[$counta] eq "ASN")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 87;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 88;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 89;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 90;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 91;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 92;
	 }
	 elsif ($atomtype[$counta] eq "OD1")
	 {
	     $position[$counta] = 93;
	 }
	 elsif ($atomtype[$counta] eq "ND2")
	 {
	     $position[$counta] = 94;
	 }
     }
     elsif ($resid[$counta] eq "PRO")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 95;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 96;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 97;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 98;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 99;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 100;
	 }
	 elsif ($atomtype[$counta] eq "CD ")
	 {
	     $position[$counta] = 101;
	 }
     }
     elsif ($resid[$counta] eq "GLN")
     {
	 if ($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 102;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 103;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 104;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 105;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 106;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 107;
	 }
	 elsif ($atomtype[$counta] eq "CD ")
	 {
	     $position[$counta] = 108;
	 }
	 elsif ($atomtype[$counta] eq "OE1")
	 {
	     $position[$counta] = 109;
	 }
	 elsif ($atomtype[$counta] eq "NE2")
	 {
	     $position[$counta] = 110;
	 }
     }
     elsif ($resid[$counta] eq "ARG")
     {
	 if($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 111;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 112;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 113;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 114;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 115;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 116;
	 }
	 elsif ($atomtype[$counta] eq "CD ")
	 {
	     $position[$counta] = 117;
	 }
	 elsif ($atomtype[$counta] eq "NE ")
	 {
	     $position[$counta] = 118;
	 }
	 elsif ($atomtype[$counta] eq "CZ ")
	 {
	     $position[$counta] = 119;
	 }
	 elsif ($atomtype[$counta] eq "NH1")
	 {
	     $position[$counta] = 120;
	 }
	 elsif ($atomtype[$counta] eq "NH2")
	 {
	     $position[$counta] = 121;
	 }
     }
     elsif ($resid[$counta] eq "SER")
     {
	 if($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 122;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 123;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 124;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 125;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 126;
	 }
	 elsif ($atomtype[$counta] eq "OG ")
	 {
	     $position[$counta] = 127;
	 }
     }
     elsif ($resid[$counta] eq "THR")
     {
	 if($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 128;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 129;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 130;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 131;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 132;
	 }
	 elsif ($atomtype[$counta] eq "OG1")
	 {
	     $position[$counta] = 133;
	 }
	 elsif ($atomtype[$counta] eq "CG2")
	 {
	     $position[$counta] = 134;
	 }
     }
     elsif ($resid[$counta] eq "VAL")
     {
	 if($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 135;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 136;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 137;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 138;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 139;
	 }
	 elsif ($atomtype[$counta] eq "CG1")
	 {
	     $position[$counta] = 140;
	 }
	 elsif ($atomtype[$counta] eq "CG2")
	 {
	     $position[$counta] = 141;
	 }
     }
     elsif ($resid[$counta] eq "TRP")
     {
	 if($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 142;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 143;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 144;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 145;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 146;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 147;
	 }
	 elsif ($atomtype[$counta] eq "CD1")
	 {
	     $position[$counta] = 148;
	 }
	 elsif ($atomtype[$counta] eq "CD2")
	 {
	     $position[$counta] = 149;
	 }
	 elsif ($atomtype[$counta] eq "NE1")
	 {
	     $position[$counta] = 150;
	 }
	 elsif ($atomtype[$counta] eq "CE2")
	 {
	     $position[$counta] = 151;
	 }
	 elsif ($atomtype[$counta] eq "CE3")
	 {
	     $position[$counta] = 152;
	 }
	 elsif ($atomtype[$counta] eq "CZ2")
	 {
	     $position[$counta] = 153;
	 }
	 elsif ($atomtype[$counta] eq "CZ3")
	 {
	     $position[$counta] = 154;
	 }
	 elsif ($atomtype[$counta] eq "CH2")
	 {
	     $position[$counta] = 155;
	 }
     }
     elsif ($resid[$counta] eq "TYR")
     {
	 if($atomtype[$counta] eq "N  ")
	 {
	     $position[$counta] = 156;
	 }
	 elsif ($atomtype[$counta] eq "CA ")
	 {
	     $position[$counta] = 157;
	 }
	 elsif ($atomtype[$counta] eq "C  ")
	 {
	     $position[$counta] = 158;
	 }
	 elsif ($atomtype[$counta] eq "O  ")
	 {
	     $position[$counta] = 159;
	 }
	 elsif ($atomtype[$counta] eq "CB ")
	 {
	     $position[$counta] = 160;
	 }
	 elsif ($atomtype[$counta] eq "CG ")
	 {
	     $position[$counta] = 161;
	 }
	 elsif ($atomtype[$counta] eq "CD1")
	 {
	     $position[$counta] = 162;
	 }
	 elsif ($atomtype[$counta] eq "CD2")
	 {
	     $position[$counta] = 163;
	 }
	 elsif ($atomtype[$counta] eq "CE1")
	 {
	     $position[$counta] = 164;
	 }
	 elsif ($atomtype[$counta] eq "CE2")
	 {
	     $position[$counta] = 165;
	 }
	 elsif ($atomtype[$counta] eq "CZ ")
	 {
	     $position[$counta] = 166;
	 }
	 elsif ($atomtype[$counta] eq "OH ")
	 {
	     $position[$counta] = 167;
	 }
     }
     elsif ($atomtype[$counta] eq "OXT")
     {
         $position[$counta] = 0;
     }

}





sub ReadDatafile
#Reads the potentials datafile into an array
{
    #open the datafile or kill the program with a message
    open(DFILE, "/nfs/home/danielle/RAM/RAM.gen.par") || die "Can't read datafile";
    {
	while(<DFILE>)
	{
	    @data = <DFILE>;
	    close (DFILE);
	}
    }
}


sub MatrixLine
#Calculates the line of the array needed from the $position of the first amino acid and the matrix number
#Splits the correct line by white-space
#Counts using the $position of the second amino acid and finds the correct potential for that pair
{
    if ($matrix eq "distance too great")
    {
	$potential = 0;
	return ($potential);
    }	
    $linenumber = ($matrix * 168) + ($position[$i] - 1);

    {
	@datal = split ('\s+', $data[$linenumber]);
	{
	    $potential = $datal[$position[$j] + 1];
	}
    }
    return ($potential);
}

