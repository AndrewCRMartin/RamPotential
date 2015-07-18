#!/usr/bin/perl -s

#Perl program to find the potentials between two residues specified in the command line

if(defined($aa))
{
    $pdbfile = $ARGV[0];
    $residueAname = $ARGV[1];
    $residueAnumber = $ARGV[2];
    $residueBname = $ARGV[3];
    $residueBnumber = $ARGV[4];
    
    $TotalEnergy = 0;
    
    #Open the PDB file
    open(PFILE, "$pdbfile") || die "Can't read $pdbfile";
    
    
    #Opens the potentials datafile and reads it into an array
    ReadDatafile();
    
    
    #Getting coordinates from PDB file
    $natoms = GetCoordinates();
    
    
    #for each pair of amino acid atom groups
    for($i = 0; $i < $natoms; $i++)
    {
	$j = $i;
	while (($ChainNameNumber[$j] eq $ChainNameNumber[$i]) && ($j < $natoms))
	{
	    $j++;
	} 
	$k = $j;
	while (($ChainNameNumber[$k] eq $ChainNameNumber[$j]) && ($k < $natoms))
	{
	    $k++;
        }
	
	#compare with all amino acids atom groups that come after it
	for ($j = $k; $j < $natoms; $j++)
	{
	    #Using pythag to calculate the distance between the two atom groups
	    $distSquared = ($gX[$i] - $gX[$j]) * ($gX[$i] - $gX[$j]) + 
		           ($gY[$i] - $gY[$j]) * ($gY[$i] - $gY[$j]) + 
		           ($gZ[$i] - $gZ[$j]) * ($gZ[$i] - $gZ[$j]);
	    {
		$distance = sqrt $distSquared;
		
		#Using the distance to find out the first part of the table number in the data file
		$matrix = MatrixNumber($distance);
		
		if ($resid[$i] eq $residueAname)
		{
		    if ($residueNumber[$i] == $residueAnumber)
		    {
			if ($resid[$j] eq $residueBname)
			{
			    if ($residueNumber[$j] == $residueBnumber)
			    {
				
				#Uses the $position and $matrix values to find the correct potential for the pair from the datafile and print it out
				$potential = MatrixLine($matrix, $position[$i], $position[$j]);
				
				$TotalEnergy = $TotalEnergy + $potential;
				
			    }
			}
		    } 
		}
	    }
	}
    }
    
    print "Total energy for amino acid pair  = ";
    print $TotalEnergy;		    
}
elsif(defined($t))
{
    $pdbfile = $ARGV[0];
    $TotalEnergy = 0;
	      
#Open the PDB file
    open(PFILE, "$pdbfile") || die "Can't read $pdbfile";
    
    
    #Opens the potentials datafile and reads it into an array
    ReadDatafile();
    
    
    #Getting coordinates from PDB file
    $natoms = GetCoordinates();
    
    
    #for each pair of amino acid atom groups
    for($i = 0; $i < $natoms; $i++)
    {
	$j = $i;
	while (($ChainNameNumber[$j] eq $ChainNameNumber[$i]) && ($j < $natoms))
	{
	    $j++;
	} 
	$k = $j;
	while (($ChainNameNumber[$k] eq $ChainNameNumber[$j]) && ($k < $natoms))
	{
	    $k++
	    }
	
	#compare with all amino acids atom groups that come after it
	for ($j = $k; $j < $natoms; $j++)
	{
	    #Using pythag to calculate the distance between the two atom groups
	    $distSquared = ($gX[$i] - $gX[$j]) * ($gX[$i] - $gX[$j]) + 
		           ($gY[$i] - $gY[$j]) * ($gY[$i] - $gY[$j]) + 
		           ($gZ[$i] - $gZ[$j]) * ($gZ[$i] - $gZ[$j]);
	    {
		$distance = sqrt $distSquared;
		
		#Using the distance to find out the first part of the table number in the data file
		$matrix = MatrixNumber($distance);
		
                #Uses the $position and $matrix values to find the correct potential for the pair from the datafile and print it out
		$potential = MatrixLine($matrix, $position[$i], $position[$j]);
				
		$TotalEnergy = $TotalEnergy + $potential;
	    }
	}
    }
    print "Total energy = $TotalEnergy\n";		    
}
elsif(defined($v))
{
    $pdbfile = $ARGV[0];
    $TotalEnergy = 0;
	      
    #Open the PDB file
    open(PFILE, "$pdbfile") || die "Can't read $pdbfile";
    
    
    #Opens the potentials datafile and reads it into an array
    ReadDatafile();
    
    
    #Getting coordinates from PDB file
    $natoms = GetCoordinates();
    
    
    #for each pair of amino acid atom groups
    for($i = 0; $i < $natoms; $i++)
    {
	$j = $i;
	while (($ChainNameNumber[$j] eq $ChainNameNumber[$i]) && ($j < $natoms))
	{
	    $j++;
	} 
	$k = $j;
	while (($ChainNameNumber[$k] eq $ChainNameNumber[$j]) && ($k < $natoms))
	{
	    $k++
	    }
	
	#compare with all amino acids atom groups that come after it
	for ($j = $k; $j < $natoms; $j++)
	{
	    #Using pythag to calculate the distance between the two atom groups
	    $distSquared = ($gX[$i] - $gX[$j]) * ($gX[$i] - $gX[$j]) + 
		           ($gY[$i] - $gY[$j]) * ($gY[$i] - $gY[$j]) + 
		           ($gZ[$i] - $gZ[$j]) * ($gZ[$i] - $gZ[$j]);
	    {
		$distance = sqrt $distSquared;
		
		#Using the distance to find out the first part of the table number in the data file
		$matrix = MatrixNumber($distance);
	
                #prints the first and second amino acid atom groups	
		print "$residueNumber[$i] $resid[$i] $atomtype[$i]"; 
		printf ("%4d", "$position[$i]");
		print " ";
		print "$residueNumber[$j] $resid[$j] $atomtype[$j]"; 
		printf ("%4d", "$position[$j]");
		
	        #prints the matrix number
		print " ";
		print "$matrix";
		print "\n";
		
                #Uses the $position and $matrix values to find the correct potential for the pair from the datafile and print it out
		$potential = MatrixLine($matrix, $position[$i], $position[$j]);
		
		print "$potential\n";
	    }
	}
    }
}
		
#################################################################################
sub GetCoordinates
#for the length of the pdb file
{
    my $counta = 0;

    while(<PFILE>)
    {
	($type, $atom, $atomname, $resname) = split;
	if ($type eq "ATOM")
	{			
	    #finds the amino acid type, the group type and the coordiantes of the group
	    $resid[$counta] = substr ($_, 17, 3);
	    $atomtype[$counta] = substr ($_, 13, 3);
	    $gX[$counta] = substr ($_, 30, 8);
	    $gY[$counta] = substr ($_, 38, 8);
	    $gZ[$counta] = substr ($_, 46, 8);	
            $ChainNameNumber[$counta] = substr ($_, 21, 6);
	    
	    $ChainNameNumber[$counta] =~ s/\s+//;

	    $residueNumber[$counta] = substr ($_, 23, 3);    

	   
	 
            #Converting residue id and atom type into a number that can be used to pull the potential out of the datafile
	    $position[$counta] = Convert($resid[$counta], $atomtype[$counta]);

	    $counta++;
	}
    }
    return ($counta);
}


##############################################################################
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
	$matrix = -1;
    }
    return($matrix);
}


##########################################################################
sub Convert
#changing the residue type and atome type to find which matrix to use and position within the matrix
{
     my ($resid, $atomtype) = @_;

     $atomtype = "O  " if($atomtype eq "OXT");
 
     if ($resid eq "ALA")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 1;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 2;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 3;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 4;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 5;
	 }
     }
     elsif ($resid eq "CYS")
     {
	 if ($atomtype eq "N  ")
	 {	
	     return 6;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 7;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 8;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 9;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 10;
	 }
	 elsif ($atomtype eq "SG ")
	 {
	     return 11;
	 }
     }
     elsif ($resid eq "ASP")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 12;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 13;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 14;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 15;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 16;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 17;
	 }
	 elsif ($atomtype eq "OD1")
	 {
	     return 18;
	 }
	 elsif ($atomtype eq "OD2")
	 {
	     return 19;
	 }
     }
     elsif ($resid eq "GLU")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 20;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 21;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 22;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 23;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 24;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 25;
	 }
	 elsif ($atomtype eq "CD ")
	 {
	     return 26;
	 }
	 elsif ($atomtype eq "OE1")
	 {
	     return 27;
	 }
	 elsif ($atomtype eq "OE2")
	 {
	     return 28;
	 }
     }
     elsif ($resid eq "PHE")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 29;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 30;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 31;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 32;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 33;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 34;
	 }
	 elsif ($atomtype eq "CD1")
	 {
	     return 35;
	 }
	 elsif ($atomtype eq "CD2")
	 {
	     return 36;
	 }
	 elsif ($atomtype eq "CE1")
	 {
	     return 37;
	 }
	 elsif ($atomtype eq "CE2")
	 {
	     return 38;
	 }
	 elsif ($atomtype eq "CZ ")
	 {
	     return 39;
	 }
     }
     elsif ($resid eq "GLY")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 40;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 41;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 42;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 43;
	 }
     }
     elsif ($resid eq "HIS")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 44;
	 } 
	 elsif ($atomtype eq "CA ")
	 {
	     return 45;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 46;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 47;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 48;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 49;
	 }
	 elsif ($atomtype eq "ND1")
	 {
	     return 50;
	 }
	 elsif ($atomtype eq "CD2")
	 {
	     return 51;
	 }
	 elsif ($atomtype eq "CE1")
	 {
	     return 52;
	 }
	 elsif ($atomtype eq "NE2")
	 {
	     return 53;
	 }
     }
     elsif ($resid eq "ILE")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 54;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 55;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 56;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 57;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 58;
	 }
	 elsif ($atomtype eq "CG1")
	 {
	     return 59;
	 }
	 elsif ($atomtype eq "CG2")
	 {
	     return 60;
	 }
	 elsif ($atomtype eq "CD1")
	 {
	     return 61;
	 }
     }
     elsif ($resid eq "LYS")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 62;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 63;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 64;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 65;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 66;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 67;
	 }
	 elsif ($atomtype eq "CD ")
	 {
	     return 68;
	 }
	 elsif ($atomtype eq "CE ")
	 {
	     return 69;
	 }	
	 elsif ($atomtype eq "NZ ")
	 {
	     return  70;
	 }
	 elsif ($atomtype eq "CZ ")
	 {
	     return 70;
	 }
     }
     elsif ($resid[$counta] eq "LEU")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 71;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 72;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 73;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 74;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 75;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 76;
	 }
	 elsif ($atomtype eq "CD1")
	 {
	     return 77;
	 }
	 elsif ($atomtype eq "CD2")
	 {
	     return 78;
	 }
     }
     elsif ($resid eq "MET")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 79;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 80;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 81;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 82;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 83;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 84;
	 }
	 elsif ($atomtyp eq "SD ")
	 {
	     return 85;
	 }
	 elsif ($atomtype eq "CE ")
	 {
	     return 86;
	 }
     }
     elsif ($resid eq "ASN")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 87;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 88;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 89;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 90;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 91;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 92;
	 }
	 elsif ($atomtype eq "OD1")
	 {
	     return 93;
	 }
	 elsif ($atomtype eq "ND2")
	 {
	     return 94;
	 }
     }
     elsif ($resid eq "PRO")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 95;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 96;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 97;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 98;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 99;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 100;
	 }
	 elsif ($atomtype eq "CD ")
	 {
	     return 101;
	 }
     }
     elsif ($resid eq "GLN")
     {
	 if ($atomtype eq "N  ")
	 {
	     return 102;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 103;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 104;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 105;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 106;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 107;
	 }
	 elsif ($atomtype eq "CD ")
	 {
	     return 108;
	 }
	 elsif ($atomtype eq "OE1")
	 {
	     return 109;
	 }
	 elsif ($atomtype eq "NE2")
	 {
	     return 110;
	 }
     }
     elsif ($resid eq "ARG")
     {
	 if($atomtype eq "N  ")
	 {
	     return 111;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 112;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 113;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 114;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 115;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 116;
	 }
	 elsif ($atomtype eq "CD ")
	 {
	     return 117;
	 }
	 elsif ($atomtype eq "NE ")
	 {
	     return 118;
	 }
	 elsif ($atomtype eq "CZ ")
	 {
	     return 119;
	 }
	 elsif ($atomtype[$counta] eq "NH1")
	 {
	     return 120;
	 }
	 elsif ($atomtype eq "NH2")
	 {
	     return 121;
	 }
     }
     elsif ($resid eq "SER")
     {
	 if($atomtype eq "N  ")
	 {
	     return 122;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 123;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 124;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 125;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 126;
	 }
	 elsif ($atomtype eq "OG ")
	 {
	     return 127;
	 }
     }
     elsif ($resid eq "THR")
     {
	 if($atomtype eq "N  ")
	 {
	     return 128;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 129;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 130;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 131;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 132;
	 }
	 elsif ($atomtype eq "OG1")
	 {
	     return 133;
	 }
	 elsif ($atomtype eq "CG2")
	 {
	     return 134;
	 }
     }
     elsif ($resid eq "VAL")
     {
	 if($atomtype eq "N  ")
	 {
	     return 135;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 136;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 137;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 138;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 139;
	 }
	 elsif ($atomtype eq "CG1")
	 {
	     return 140;
	 }
	 elsif ($atomtype eq "CG2")
	 {
	     return 141;
	 }
     }
     elsif ($resid eq "TRP")
     {
	 if($atomtype eq "N  ")
	 {
	     return 142;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 143;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 144;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 145;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 146;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 147;
	 }
	 elsif ($atomtype eq "CD1")
	 {
	     return 148;
	 }
	 elsif ($atomtype eq "CD2")
	 {
	     return 149;
	 }
	 elsif ($atomtype eq "NE1")
	 {
	     return 150;
	 }
	 elsif ($atomtype eq "CE2")
	 {
	     return 151;
	 }
	 elsif ($atomtype eq "CE3")
	 {
	     return 152;
	 }
	 elsif ($atomtype eq "CZ2")
	 {
	     return 153;
	 }
	 elsif ($atomtype eq "CZ3")
	 {
	     return 154;
	 }
	 elsif ($atomtype eq "CH2")
	 {
	     return 155;
	 }
     }
     elsif ($resid eq "TYR")
     {
	 if($atomtype eq "N  ")
	 {
	     return 156;
	 }
	 elsif ($atomtype eq "CA ")
	 {
	     return 157;
	 }
	 elsif ($atomtype eq "C  ")
	 {
	     return 158;
	 }
	 elsif ($atomtype eq "O  ")
	 {
	     return 159;
	 }
	 elsif ($atomtype eq "CB ")
	 {
	     return 160;
	 }
	 elsif ($atomtype eq "CG ")
	 {
	     return 161;
	 }
	 elsif ($atomtype eq "CD1")
	 {
	     return 162;
	 }
	 elsif ($atomtype eq "CD2")
	 {
	     return 163;
	 }
	 elsif ($atomtype eq "CE1")
	 {
	     return 164;
	 }
	 elsif ($atomtype eq "CE2")
	 {
	     return 165;
	 }
	 elsif ($atomtype eq "CZ ")
	 {
	     return 166;
	 }
	 elsif ($atomtype eq "OH ")
	 {
	     return 167;
	 }
     }
     elsif ($atomtype eq "OXT ")
     {
         return 0;
     }
}




###########################################################################
sub ReadDatafile
#Reads the potentials datafile into an array
{
    #open the datafile or kill the program with a message
    open(DFILE, "/nfs/home/danielle/RAM/RAM.gen.par") || die "Can't read datafile";
    
    @gData = <DFILE>;
    close (DFILE);
}


###########################################################################
sub MatrixLine
#Calculates the line of the array needed from the $position of the first amino acid and the matrix number
#Splits the correct line by white-space
#Counts using the $position of the second amino acid and finds the correct potential for that pair
{
    my($matrix, $posi, $posj) = @_;
    my($linenumber, $dataline, @datal);

    if (($matrix == -1) || ($posi == 0) || ($posj == 0))
    {
	return 0;
    }	
    $linenumber = ($matrix * 168) + ($posi-1);
    $dataline = $gData[$linenumber];
    $dataline =~ s/^\s+//;
    @datal = split(/\s+/, $dataline);

    return ($datal[$posj]);
}

