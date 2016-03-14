#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_makefile_EMMCMC

=head1 SYNOPSIS

generate_makefile_EMMCMC [options] 

Options:

  -h      brief help message
  --man   full documentation
  --mf    output make file

Example usage: ./generate_makefile_EMMCMC.pl -h

=head1 OPTIONS

=over 8

=item B<-help>

  -w         work directory : location of all output files
  -t         gemma tool directory
  --vd       vcf file directory
  --ad       annotation file directory
  --ac       annotation code file
  --pheno    phenotype file
  --hyp      initial hyper parameter value file
  -f         file with a list of vcf file heads
  -G         genotype format: GT(genotype data 0/0, 0/1, 1/1) or EC (dosage data)
  --maf      maf threshold: default 0.5% 
  --pp       specify prior for the causal probability: default 1e-6
  --abgamma  specify inverse gamma prior for the effect size variance: default 0.1
  --win      window size of the neighborhood: default 100
  --em       number of EM iterations: default 5
  -b         number of burn ins: default 50,000
  -N         number of MCMC iterations: default 50,000
  --NL       number of MCMC iterations for the last EM iteration: default 50,000
  -c         compress genotype data or not: default 0 (do not compress)
  --initype  specify initial model: default 3 (stepwise selection)
  --rv       fixed residual variance value: default 1
  --smin     minimum number of variates per block in the model: default 0
  --smax     maximum number of variates per block in the model: default 5
  --mem      specify maximum memory usage: default 2Gb
  --time     specify time of runing MCMC per block per MCMC iteration: default 24hr
  -l         specify how job will be submited: default srun
  --mf       output make file

=item B<-man>

prints the manual page and exits.

=back

=head1 DESCRIPTION

B<generate_makefile_EMMCMC.pl> will generate a makefile to be run in the Unix/Linux/OS/Window system with the Bayesian EM_MCMC algorithm for GWASs. 

=cut

# define default option variables
my $help;
my $verbose;
my $debug;
my $man;
my $launchMethod = "slurm";
my $wkDir=getcwd();
my $makeFile = "Makefile_EM_MCMC.mk";

my $toolDir="/net/wonderland/home/yjingj/AMD";
my $annoDir="/net/wonderland/home/yjingj/Data/AMD/ImputeLoci_BlockVCFs/GVS_Anno";
my $vcfDir = "/net/wonderland/home/yjingj/Data/AMD/ImputeLoci_BlockVCFs";
my $pheno="/net/wonderland/home/yjingj/AMD/amd_pheno_scaled.txt.gz";
my $annoCode="/net/wonderland/home/yjingj/AMD/AnnoPartition/Func7.txt";
my $hyppar="/net/wonderland/home/yjingj/AMD/EM_Initial/init_hyper_group7";
my $filelist = "/net/wonderland/home/yjingj/Data/AMD/ImputeLoci_BlockVCFs/vcf_filehead.txt";
my $rs="/net/fantasia/home/yjingj/My_Rcode/Mstep_multgroup.r";

my $EM=5;
my $GTfield="GT";
my $maf="0.005";
my $rho="1";
my $smin="0";
my $smax="5";
my $win="100";
my $burnin="50000";
my $Nmcmc="50000";
my $NmcmcLast="50000";
my $compress="0";
my $initype="3";
my $rv="1";
my $pp="1e-7";
my $abgamma="0.1";
my $gemma="gemma_emblock_multgroup";
my $saveSNP="0";

my $maxmem = "3000";
my $time = "24:00:00";
my $nice = "100";
my $jobid;
my $xnode="";
my $wnode="";
my $part="nomosix";

# "twins-mc[01-04],1000g-mc[01-04],amd-mc[02-04],c43,got2d-mc[01-04],psoriasis-mc[01-04],r6313"; # excluded host list

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug, 'm'=>\$man,
                'w:s'=>\$wkDir, 't:s' =>\$toolDir, 'ad:s'=>\$annoDir, 
                'ac:s'=>\$annoCode, 'vd:s'=>\$vcfDir, 'pheno:s'=>\$pheno,
                'hyp:s'=>\$hyppar, 'G:s'=>\$GTfield, 'maf:s'=>\$maf,
                'smin:s'=>\$smin, 'smax:s'=>\$smax, 'win:s'=>\$win,
                'b:s'=>\$burnin, 'N:s'=>\$Nmcmc, 'NL:s'=>\$NmcmcLast,
                'c:s'=>\$compress, 'initype:s'=>\$initype, 'rv:s'=>\$rv,
                'pp:s'=>\$pp, 'abgamma:s'=>\$abgamma, 
                'mem:s'=>\$maxmem, 'saveSNP:s'=>\$saveSNP,
                'time:s'=>\$time, 'f:s'=>\$filelist, 'em:i'=>\$EM, 'rs:s'=>\$rs,
                'l:s'=>\$launchMethod, 'mf:s'=>\$makeFile, 'nice:s'=>\$nice, 'gemma:s'=>\$gemma, 'j:s' =>\$jobid, 'xnode:s'=>\$xnode, 
                'wnode:s'=>\$wnode, 'part:s'=>\$part)
  || !defined($wkDir) || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(1);
        exit(0);
    }
    elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }
    else
    {
        pod2usage(1);
        exit(0);
    }
}

if ($help)
    {
        pod2usage(1);
        exit(0);
    }
elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }


if ($launchMethod ne "local" && $launchMethod ne "slurm" && $launchMethod ne "mosix" && $launchMethod ne "qsub")
{
    print STDERR "Launch method has to be local, slurm, mosix, or qsub\! \n";
    exit(1);
}

##############
#print options
##############
printf("Options\n");
printf("\n");
printf("launch method : %s\n", $launchMethod);
printf("work directory : %s\n", $wkDir);
print "vcfDir: ", $vcfDir, "\n", "annoDir: ", $annoDir, "\n",
        "pheno: ", $pheno, "\nannoCode: ", $annoCode, "\n", 
        "hyppar: ", $hyppar, "\nfileheads: ", $filelist, "\n", 
        "rscript: ", $rs, "\n",
        "GTfield ", $GTfield, "; maf ", $maf, "; smin ", $smin, "\n", 
        "smax ", $smax, "; win ", $win, "; burnin ", $burnin, "; Nmcmc ", $Nmcmc, "\n",
        "NmcmcLast ", $NmcmcLast, "; compress ", $compress, "; initype ", $initype, "\n",
        "rv ", $rv, "; pp ", $pp, "; abgamma ", $abgamma, "\n"; 
printf("\n");


#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();

#temporary variables
my $tgt;
my $dep;
my @cmd;

mkpath($wkDir);

########################################
# Initial Set up before EM iterations
########################################

### prepare files before EM_MCMC
my $hypcurrent="$wkDir/hypval.current";
$tgt = "$wkDir/pre_em.OK";
$dep = "";
@cmd = "rm -f -r $wkDir/output $wkDir/Eoutput $wkDir/OUT";
push(@cmd, "mkdir -p $wkDir/output $wkDir/Eoutput $wkDir/OUT");
push(@cmd, "cp -f $hyppar $hypcurrent");
push(@cmd, "> $wkDir/Eoutput/EM\_result.txt");
push(@cmd, "> $wkDir/Rout.txt");
makeJob("local", $tgt, $dep, $wkDir, @cmd);  


###### EM step 0 without dependencies ###########
my $i;
my $premcmcOK="";
my @filehead;
my $line;

open(my $FILELIST, $filelist)
    or die "Can not open $filelist \!";
 while ($line = <$FILELIST>) {
    chop $line;
    push(@filehead, $line);
}
close $FILELIST;

if(@filehead == 0) {
    print STDERR "file list is empty\! \n";
    exit(1);
} 
else{ print "Total \# of fileheads: ", scalar(@filehead), "\n \n"; }


for(my $j=0; $j< @filehead; ++$j)
    {
        $i=0;
        $line = $filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.$i.OK ";
        $tgt = "$wkDir/OUT/$line.$i.OK";
        $dep = "$wkDir/pre_em.OK";
        @cmd = "$toolDir/$gemma -vcf $vcfDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -vcfp $pheno -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -bslmm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -n 1 -o $line -w $burnin -s $Nmcmc -comp $compress -saveSNP $saveSNP -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

my $paramfile="$wkDir/Eoutput/paramtemp$i.txt";
my $hypfile="$wkDir/Eoutput/hyptemp$i.txt";
my $logfile="$wkDir/Eoutput/log$i.txt";

$tgt = "$wkDir/Eoutput/cp_param$i.OK";
$dep = "$premcmcOK";
@cmd = "cat \`ls -d -1 $wkDir/output/** | grep paramtemp | sort\` > $paramfile";
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep hyptemp | sort\` > $hypfile");
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
makeJob("local", $tgt, $dep, $wkDir, @cmd);

$tgt = "$wkDir/R$i.OK";
$dep = "$wkDir/Eoutput/cp_param$i.OK $wkDir/pre_em.OK";
@cmd = "Rscript --no-save --no-restore --verbose $rs $hypfile $i $pp $abgamma $wkDir/Eoutput/EM_result.txt $hypcurrent >> $wkDir/Rout.txt";
makeJob("local", $tgt, $dep, $wkDir, @cmd);


####### With dependencies of previous output 
my $ipre="";

for $i (1..$EM){

    $ipre=$i-1; $premcmcOK="";

    for(my $j=0; $j< @filehead; ++$j){
        $line=$filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.$i.OK ";
        $tgt = "$wkDir/OUT/$line.$i.OK";
        $dep = "$wkDir/R$ipre.OK";
        if($i < $EM){
            @cmd = "$toolDir/gemma\_emblock\_multgroup -vcf $vcfDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -vcfp $pheno -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -bslmm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -n 1 -o $line -w $burnin -s $Nmcmc -comp $compress -saveSNP $saveSNP -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt ";
            } elsif ($i == $EM){
            @cmd = "$toolDir/gemma\_emblock\_multgroup -vcf $vcfDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -vcfp $pheno -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -bslmm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -n 1 -o $line -w $burnin -s $NmcmcLast -comp $compress -saveSNP $saveSNP -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt ";
            }
        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

    $paramfile="$wkDir/Eoutput/paramtemp$i.txt";
    $hypfile="$wkDir/Eoutput/hyptemp$i.txt";
    $logfile="$wkDir/Eoutput/log$i.txt";

  $tgt = "$wkDir/Eoutput/cp_param$i.OK";
  $dep = "$premcmcOK";
  @cmd = "cat \`ls -d -1 $wkDir/output/** | grep paramtemp | sort \` > $paramfile";
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep hyptemp | sort\` > $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
  makeJob("local", $tgt, $dep, $wkDir, @cmd);

  $tgt = "$wkDir/R$i.OK";
  $dep = "$wkDir/Eoutput/cp_param$i.OK";
  @cmd = "Rscript --no-save --no-restore --verbose $rs $hypfile $i $pp $abgamma $wkDir/Eoutput/EM_result.txt $hypcurrent >> $wkDir/Rout.txt";
  makeJob("local", $tgt, $dep, $wkDir, @cmd);

}


#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
#this tells GNU Make to delete files that were generated during a command that did not end successfully.
print MAK ".DELETE_ON_ERROR:\n\n";
#this ensures that all the targets are executed; exclude clean
print MAK "all: @tgts\n\n";

######## Create clean jobs command #######
mkpath("$wkDir/slurm_err");
push(@tgts, "clean_err");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/slurm_err/*.err");


push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/*.OK $wkDir/Eoutput/*.OK $wkDir/OUT/*.OK $wkDir/slurm_err/*.err");


for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, $wkd, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    elsif ($method eq "slurm")
    {
        makeSlurm($tgt, $dep, $wkd, @cmd);
    }
    elsif($method eq "mosix")
    {
      makeMosix($tgt, $dep, $wkd, @cmd);
    }
}

#run mosix jobs
sub makeMosix
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -b -q1 " . $c . "\n";
    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run slurm jobs
sub makeSlurm
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
        $cmdtemp .= "\tsrun --exclude=$xnode --partition=$part --mem-per-cpu\=$maxmem --time\=$time --nice\=$nice --error\=$wkDir/slurm_err/\%N.\%j.err -J $jobid -D $wkd $c \n";

        #$cmdtemp .= "\tsrun --exclude=$xnode --nodelist=$wnode --mem-per-cpu\=$maxmem --time\=$time --nice\=$nice --error\=$wkDir/slurm_err/\%N.\%j.err -J $jobid -D $wkd $c \n";
    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

# Check empty directories
sub dir_is_empty
{
    my ($path) = @_;
    opendir DIR, $path;
    while(my $entry = readdir DIR) {
        next if($entry =~ /^\.\.?$/);
        closedir DIR;
        return 0;
    }
    closedir DIR;
    return 1;
}
