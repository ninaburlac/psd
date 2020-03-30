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

my $cycle = $ARGV[0];
my $run = $ARGV[1];

my $myUser = $ENV{'USER'};

# save current umask
my $old_umask = umask;
umask 0000;

my $resdir = $localdir . "/analysis";
mkdir "$resdir", 0770 unless -d "$resdir";

my @cmd = ( "!readcommand!");


my $readprogram = $localdir . "/readTier";
my $readout = $resdir . "/read.out";
my $readcommand = $readprogram . " " . $resdir . " " . $cycle . " " . $run . " > " . $readout; 

### submit noise script
my $scripttobecopied = $localdir . "/scriptread.template.sh";
open IN, $scripttobecopied or die "Can't read source file $scripttobecopied: $!\n";

my $currscriptr = $resdir . "/scriptread.sh";
open OUT, ">$currscriptr" or die "Can't write on file $currscriptr: $!\n";

while (<IN>) {
    s/$cmd[0]/$readcommand/g;
    print OUT $_;
}
close IN;
close OUT;

### submit READ job to queue
my $scriptlogr = $resdir . "/read_1.out";
my $scripterrr = $resdir . "/read_1.err";
my $jobNamer = "read";

my $QUEUEcmdr = "qsub -N " . $jobNamer . " -q gerda -V -d " . $resdir . " -m abe -e localhost:". $scripterrr . " -o localhost:" . $scriptlogr . " -l mem=4000mb " . $currscriptr;
system($QUEUEcmdr);
print $QUEUEcmdr . "\n";

my $countJobr = "qstat -u " . $myUser . " |  grep read | wc -l ";
print $countJobr . "\n";
my $inQueuer = 1;

while ($inQueuer) {
    sleep 30;
    my $actualJobr = `$countJobr`;
    print $actualJobr;
    if($actualJobr==0) {
	$inQueuer=0;
	
        # Reset umask to old value
	umask $old_umask;
	
	my $from = 'psdtest';
	my $to = 'ninaburlac.nb\@gmail.com';
	my $subject = "read root files";
	#my $fwhmFile = $scandir . "/FWHM_chn" . $chn . ".txt";
	my $body = "Hi, \n this is an automatic message from the test of the PSD";
	my $sendMail = "echo \"". $body . "\" | mailx -s \"". $subject . "\" -r " . $from . " " . $to;
	system($sendMail);
    }
}


