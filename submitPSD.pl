#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use YAML::Tiny 'LoadFile';
use YAML::Tiny;
use Data::Dumper;
use POSIX qw(ceil floor);

use constant false => 0;
use constant true  => 1;

my $localdir = $ENV{PWD};

#Open the configuration file
#my $configfile = YAML::Tiny::LoadFile( $localdir . "/config.yml");

my $firstchn = $ARGV[0];#$configfile->{FirstChannel};
my $nchn = $ARGV[1];#$configfile->{NChannels};
my $lastchn = $firstchn + $nchn;
my $tier2file = $ARGV[2];

my $myUser = $ENV{'USER'};

# save current umask
my $old_umask = umask;
umask 0000;

for(my $chn = $firstchn; $chn < $lastchn; $chn++){
    my $thisdir = $localdir . "/chn" . $chn;
    mkdir "$thisdir", 0770 unless -d "$thisdir";
    
    my @cmd = ( "!psdcommand!");
    
    my $psdprogram = $localdir . "/processPSD";
    my $psdout = $thisdir . "/psd_chn" . $chn . ".out";
    my $psdcommand = $psdprogram . " " . $thisdir . " " . $tier2file . " " . $chn . " > " . $psdout; 
    
    ### submit noise script
    
    my $scripttobecopied = $localdir . "/script.template.sh";
    open IN, $scripttobecopied or die "Can't read source file $scripttobecopied: $!\n";
    
    my $currscript = $thisdir . "/script_chn" . $chn . ".sh";
    open OUT, ">$currscript" or die "Can't write on file $currscript: $!\n";
    
    while (<IN>) {
	s/$cmd[0]/$psdcommand/g;
	print OUT $_;
    }
    close IN;
    close OUT;
    
    ### submit PSD job to queue
    my $scriptlog = $thisdir . "/psd_" . $chn . ".out";
    my $scripterr = $thisdir . "/psd_" . $chn . ".err";
    my $jobName = "psd_" . $chn;
    
    my $QUEUEcmd = "qsub -N " . $jobName . " -q gerda -V -d " . $thisdir . " -m abe -e localhost:". $scripterr . " -o localhost:" . $scriptlog . " -l mem=4000mb " . $currscript;
    system($QUEUEcmd);
    print $QUEUEcmd . "\n";
    
}

my $countJob = "qstat -u " . $myUser . " |  grep psd | wc -l ";
print $countJob . "\n";
my $inQueue = 1;

while ($inQueue) {
    sleep 30;
    my $actualJob = `$countJob`;
    print $actualJob;
    if($actualJob==0) {
	$inQueue=0;
	
        # Reset umask to old value
	umask $old_umask;
	
	my $from = 'psdtest';
	my $to = 'valerio.dandrea\@lngs.infn.it';
	my $subject = "psd";
	#my $fwhmFile = $scandir . "/FWHM_chn" . $chn . ".txt";
	my $body = "Hi, \n this is an automatic message from the test of the PSD";
	my $sendMail = "echo \"". $body . "\" | mailx -s \"". $subject . "\" -r " . $from . " " . $to;
	system($sendMail);
    }
}


