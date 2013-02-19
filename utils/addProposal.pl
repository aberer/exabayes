
#WARNING needs to access /src and /examples from working directory. i.e. must be called as 'perl utils/addProposal.pl'
#This scipt modifies all relevant files to allow for the new proposal. Parameters/Names need to be changed to the required values in this file.
#
#The following still need to be taken care of manualy: actually implement apply (and reset) function in proposals.c. printf the new values in output.c, add values to config files

#################################################################################
$proposalName="UPDATE_SINGLE_BL_BIUNIF"; 
$weightName="initSingleBrancBiunifhWeight";
$initWeight="0.0";

$configName="initSingleBranchBiunifWeight";
$configWeight="1";

$nexConfigName="initSingleBranchBiunifWeight";#may be same as $configName
$nexConfigWeight="1";

$applyName="biunif_branch_length_proposal_apply";
$resetName="random_branch_length_proposal_reset";#may already exist
$priorName="get_branch_length_prior";#may already exist


###################################################################################

#bayes.c###########################
$file="src/bayes.c";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD addInitParameters/){
print WRITEFILE "curstate->proposalWeights[".$proposalName."] = initParams->".$weightName.";\n" ;

}elsif($line =~ /PROPOSALADD initDefaultValues/){
print WRITEFILE "theState->proposalWeights[".$proposalName."] = ".$initWeight.";\n" ;

}
print WRITEFILE $line;


}

close READFILE;
close WRITEFILE;

system ("mv $temp tempBayes.c");

#configParser.c#####################

$file="src/configParser.c";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD parseLine/){
print WRITEFILE "else if( ! strcmp(key, \"$configName\"))\n";
print WRITEFILE "theState->proposalWeights[$proposalName] = atof(value);\n" ;
}

print WRITEFILE $line;

}

close READFILE;
close WRITEFILE;

system ("mv $temp tempConfigParser.c");

#nclConfigReader.cpp###########################
$file="src/nclConfigReader.cpp";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD paramBadInit/){
print WRITEFILE "initParam->$nexConfigName  = -1 ;\n" ;

}elsif($line =~ /PROPOSALADD read/){
print WRITEFILE "else if (key.EqualsCaseInsensitive(\"$nexConfigName\"))\n" ;
print WRITEFILE "initParam->$nexConfigName = value.ConvertToDouble();\n" ;

}elsif($line =~ /PROPOSALADD assertInitialized/){
print WRITEFILE "assert(initParam->$nexConfigName != -1);\n" ;

}
print WRITEFILE $line;
}

close READFILE;
close WRITEFILE;

system ("mv $temp tempNclConfigReader.cpp");

#nclConfigReader.h###########################
$file="src/nclConfigReader.h";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD initParamStruct/){
print WRITEFILE "double $nexConfigName;\n" ;
}
print WRITEFILE $line;
}

close READFILE;
close WRITEFILE;

system ("mv $temp tempNclConfigReader.h");


#proposals.c###########################
$file="src/proposals.c";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD prop_funcs/){
print WRITEFILE "{ $proposalName, $applyName, $resetName, $priorName},\n" ;
}
print WRITEFILE $line;
}

close READFILE;
close WRITEFILE;

system ("mv $temp tempProposals.c");

#proposalStructs.h###########################
$file="src/proposalStructs.h";
$temp="temp.txt";

$numProposals=-1;

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD proposal_type/){
$string="//PROPOSALADD";
@splitline=split(/$string/,$line);
print WRITEFILE "$splitline[0],\n" ;
print WRITEFILE "$proposalName = ".($numProposals-1) ;
print WRITEFILE "//PROPOSALADD$splitline[1]\n" ;
}elsif($line =~ /PROPOSALADD NUM_PROPOSALS/){
@splitline=split(/\(/,$line);
print WRITEFILE "$splitline[0](";
@splitline=split(/\)/,$splitline[1]);
$numProposals=$splitline[0]+1;
print WRITEFILE "".$numProposals.")".$splitline[1];


}else{
print WRITEFILE $line;
}
}

close READFILE;
close WRITEFILE;

system ("mv $temp tempProposalStructs.h");

#test.nex###########################
$file="examples/test.nex";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD nexusConfig/){
print WRITEFILE "      $nexConfigName $nexConfigWeight\n" ;
}
print WRITEFILE $line;
}

close READFILE;
close WRITEFILE;

system ("mv $temp temptTest.nex");

#test.conf###########################
$file="examples/test.conf";
$temp="temp.txt";

open READFILE, "<$file" or die $!;
open WRITEFILE, ">$temp" or die $!;
while (my $line = <READFILE>) {
if ($line =~ /PROPOSALADD confConfig/){
print WRITEFILE "$configName=$configWeight\n" ;
}
print WRITEFILE $line;
}

close READFILE;
close WRITEFILE;

system ("mv $temp temptTest.conf");

