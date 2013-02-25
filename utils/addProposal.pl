
#WARNING needs to access /src and /examples from working directory. i.e. must be called as 'perl utils/addProposal.pl'
#This scipt modifies all relevant files to allow for the new proposal. Parameters/Names need to be changed to the required values in this file.
#
#If run with the argument 'test', the original files are kept and temp* files are kept in the working directory. 
#
#The following still need to be taken care of manualy: actually implement apply (and reset) function in proposals.c. printf the new values in output.c, add values to config files

#################################################################################
$proposalName="UPDATE_SINGLE_BL_BIUNIF"; 
$weightName="initSingleBranchBiunifWeight";
$initWeight="0.0";

#$configName="initSingleBranchBiunifWeight";
$configName=$weightName;
$configWeight="4.0";

#$nexConfigName="initSingleBranchBiunifWeight";#may be same as $configName
$nexConfigName=$weightName;
$nexConfigWeight="4.0";

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

#########################Move#And#Copy#Files################
$backupdir="utils/backup";
if($ARGV[0]!~ /test/){
system("mkdir $backupdir");

system ("mv src/bayes.c $backupdir/bayes.c");
system ("mv src/configParser.c $backupdir/configParser.c");
system ("mv src/nclConfigReader.cpp $backupdir/nclConfigReader.cpp");
system ("mv src/nclConfigReader.h $backupdir/nclConfigReader.h");
system ("mv src/proposals.c $backupdir/proposals.c");
system ("mv src/proposalStructs.h $backupdir/proposalStructs.h");
system ("mv examples/test.nex $backupdir/test.nex");
system ("mv examples/test.conf $backupdir/test.conf");

system ("mv tempBayes.c src/bayes.c");
system ("mv tempConfigParser.c src/configParser.c");
system ("mv tempNclConfigReader.cpp src/nclConfigReader.cpp");
system ("mv tempNclConfigReader.h src/nclConfigReader.h");
system ("mv tempProposals.c src/proposals.c");
system ("mv tempProposalStructs.h src/proposalStructs.h");
system ("mv temptTest.nex examples/test.nex");
system ("mv temptTest.conf examples/test.conf");
}